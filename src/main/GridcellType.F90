module GridcellType

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Gridcell data type allocation 
  ! -------------------------------------------------------- 
  ! gridcell types can have values of 
  ! -------------------------------------------------------- 
  !   1 => default
  !
  use shr_kind_mod   , only : r8 => shr_kind_r8
  use shr_infnan_mod , only : nan => shr_infnan_nan, assignment(=)
  use clm_varcon     , only : ispval
  !
  ! !PUBLIC TYPES:
  implicit none
  save
  private
  !
  type, public :: gridcell_type

     ! topological mapping functionality, local 1d gdc arrays
     integer , pointer :: gindex           (:) ! global index
     real(r8), pointer :: area             (:) ! total land area, gridcell (km^2)
     real(r8), pointer :: lat              (:) ! latitude (radians)
     real(r8), pointer :: lon              (:) ! longitude (radians)
     real(r8), pointer :: latdeg           (:) ! latitude (degrees)
     real(r8), pointer :: londeg           (:) ! longitude (degrees)
     logical , pointer :: active           (:) ! just needed for symmetry with other subgrid types

     ! Daylength
     real(r8) , pointer :: max_dayl        (:) ! maximum daylength for this grid cell (s)
     real(r8) , pointer :: dayl            (:) ! daylength (seconds)
     real(r8) , pointer :: prev_dayl       (:) ! daylength from previous timestep (seconds)

   contains

     procedure, public :: Init
     procedure, public :: Clean
     
  end type gridcell_type
  type(gridcell_type), public, target :: grc    !gridcell data structure
  !------------------------------------------------------------------------

contains
  
  !------------------------------------------------------------------------
  subroutine Init(this, begg, endg)
    !
    ! !ARGUMENTS:
    class(gridcell_type) :: this
    integer, intent(in)  :: begg, endg
    !------------------------------------------------------------------------

    ! The following is set in InitGridCells
    allocate(this%gindex    (begg:endg)) ; this%gindex    (:) = ispval
    allocate(this%area      (begg:endg)) ; this%area      (:) = nan
    allocate(this%lat       (begg:endg)) ; this%lat       (:) = nan
    allocate(this%lon       (begg:endg)) ; this%lon       (:) = nan
    allocate(this%latdeg    (begg:endg)) ; this%latdeg    (:) = nan
    allocate(this%londeg    (begg:endg)) ; this%londeg    (:) = nan
    allocate(this%active    (begg:endg)) ; this%active    (:) = .true.

    ! This is initiailized in module DayLength
    allocate(this%max_dayl  (begg:endg)) ; this%max_dayl  (:) = nan
    allocate(this%dayl      (begg:endg)) ; this%dayl      (:) = nan
    allocate(this%prev_dayl (begg:endg)) ; this%prev_dayl (:) = nan

  end subroutine Init

  !------------------------------------------------------------------------
  subroutine Clean(this)
    !
    ! !ARGUMENTS:
    class(gridcell_type) :: this
    !------------------------------------------------------------------------

    deallocate(this%gindex           )
    deallocate(this%area             )
    deallocate(this%lat              )
    deallocate(this%lon              )
    deallocate(this%latdeg           )
    deallocate(this%londeg           )
    deallocate(this%active           )
    deallocate(this%max_dayl         )
    deallocate(this%dayl             )
    deallocate(this%prev_dayl        )

  end subroutine Clean

end module GridcellType
