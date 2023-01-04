module initInterpMultilevelContainer

  ! ------------------------------------------------------------------------
  ! !DESCRIPTION:
  ! This module defines a class that contains one instance of each interp_multilevel
  ! type. This class is responsible for:
  !
  ! (1) Constructing each interp_multilevel object
  !
  ! (2) Determining which interp_multilevel object should be used for each multi-level
  !     field (based on the field's level dimension)
  !
  ! !USES:
#include "shr_assert.h" 
  use shr_kind_mod               , only : r8 => shr_kind_r8
  use initInterpBounds           , only : interp_bounds_type
  use initInterpMultilevelBase   , only : interp_multilevel_type
  use initInterpMultilevelCopy   , only : interp_multilevel_copy_type
  use initInterpMultilevelInterp , only : interp_multilevel_interp_type
  use initInterpMultilevelSnow   , only : interp_multilevel_snow_type
  use initInterpMultilevelSplit  , only : interp_multilevel_split_type, create_interp_multilevel_split_type
  use initInterp2dvar            , only : interp_2dvar_type
  use initInterp1dData           , only : interp_1d_data
  use ncdio_pio                  , only : file_desc_t, var_desc_t, check_var, ncd_io
  use clm_varctl                 , only : iulog
  use abortutils                 , only : endrun
  use shr_log_mod                , only : errMsg => shr_log_errMsg
  use spmdMod                    , only : masterproc
  use array_utils                , only : transpose_wrapper

  implicit none
  private
  save

  ! Public types

  public :: interp_multilevel_container_type

  type :: interp_multilevel_container_type
     private
     ! Components need to be pointers so that we can return pointers to them.
     !
     ! (Components of a derived type cannot have the target attribute, but rather take on
     ! the target attribute from their parent object. So the alternative to making these
     ! pointers would be to require all instances of this derived type to have the target
     ! attribute.)
     type(interp_multilevel_copy_type), pointer   :: interp_multilevel_copy
   contains
     procedure :: find_interpolator
  end type interp_multilevel_container_type

  interface interp_multilevel_container_type
     module procedure constructor
  end interface interp_multilevel_container_type

  character(len=*), parameter, private :: sourcefile = &
       __FILE__

contains

  ! ========================================================================
  ! Constructors
  ! ========================================================================

  !-----------------------------------------------------------------------
  function constructor(ncid_source, ncid_dest, bounds_source, bounds_dest, &
       pftindex, colindex) result(this)
    !
    ! !DESCRIPTION:
    ! Create an interp_multilevel_container_type instance.
    !
    ! !USES:
    use ncdio_pio, only : file_desc_t
    ! 
    ! !ARGUMENTS:
    type(interp_multilevel_container_type) :: this  ! function result
    type(file_desc_t), target, intent(inout) :: ncid_source ! netcdf ID for source file
    type(file_desc_t), target, intent(inout) :: ncid_dest   ! netcdf ID for dest file
    type(interp_bounds_type), intent(in) :: bounds_source
    type(interp_bounds_type), intent(in) :: bounds_dest

    ! The following give mappings from source to dest for pft and col-level variables.
    ! e.g., colindex(i) gives source col corresponding to dest col i.
    integer, intent(in) :: pftindex(:)
    integer, intent(in) :: colindex(:)

    ! !LOCAL VARIABLES:

    character(len=*), parameter :: subname = 'constructor'
    !-----------------------------------------------------------------------

    allocate(this%interp_multilevel_copy)
    this%interp_multilevel_copy = interp_multilevel_copy_type()

  end function constructor

  ! ========================================================================
  ! Public methods
  ! ========================================================================

  !-----------------------------------------------------------------------
  function find_interpolator(this, lev_dimname, vec_dimname) result(interpolator)
    !
    ! !DESCRIPTION:
    ! Given the name of the level dimension and the vector dimension, return a pointer to
    ! an interpolator that is appropriate for this multi-level variable.
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    class(interp_multilevel_type), pointer :: interpolator  ! function result

    class(interp_multilevel_container_type), intent(in) :: this
    character(len=*), intent(in) :: lev_dimname  ! name of level dimension
    character(len=*), intent(in) :: vec_dimname  ! name of vector dimension
    !
    ! !LOCAL VARIABLES:

    character(len=*), parameter :: subname = 'find_interpolator'
    !-----------------------------------------------------------------------

    select case (lev_dimname)
    case ('levgrnd')
       select case (vec_dimname)
       case default
          call error_not_found(subname, lev_dimname, vec_dimname)
       end select
    case default
       interpolator => this%interp_multilevel_copy
    end select

  contains
    subroutine error_not_found(subname, lev_dimname, vec_dimname)
      ! Write an error message and abort
      character(len=*), intent(in) :: subname
      character(len=*), intent(in) :: lev_dimname
      character(len=*), intent(in) :: vec_dimname

      write(iulog,*) subname//' ERROR: no multi-level interpolator found for:'
      write(iulog,*) 'lev_dimname = ', trim(lev_dimname)
      write(iulog,*) 'vec_dimname = ', trim(vec_dimname)
      call endrun(msg='ERROR: no multi-level interpolator found '//errMsg(sourcefile, __LINE__))
    end subroutine error_not_found

  end function find_interpolator

end module initInterpMultilevelContainer
