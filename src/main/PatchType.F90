module PatchType

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Patch data type allocation 
  ! -------------------------------------------------------- 
  ! patch types can have values of
  ! -------------------------------------------------------- 
  !   0  => not_vegetated
  !   1  => needleleaf_evergreen_temperate_tree
  !   2  => needleleaf_evergreen_boreal_tree
  !   3  => needleleaf_deciduous_boreal_tree
  !   4  => broadleaf_evergreen_tropical_tree
  !   5  => broadleaf_evergreen_temperate_tree
  !   6  => broadleaf_deciduous_tropical_tree
  !   7  => broadleaf_deciduous_temperate_tree
  !   8  => broadleaf_deciduous_boreal_tree
  !   9  => broadleaf_evergreen_shrub
  !   10 => broadleaf_deciduous_temperate_shrub
  !   11 => broadleaf_deciduous_boreal_shrub
  !   12 => c3_arctic_grass
  !   13 => c3_non-arctic_grass
  !   14 => c4_grass
  !   15 => c3_crop
  !   16 => c3_irrigated
  ! -------------------------------------------------------- 
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
  type, public :: patch_type

     ! g/l/c/p hierarchy, local g/l/c/p cells only
     integer , pointer :: column   (:) ! index into column level quantities
     real(r8), pointer :: wtcol    (:) ! weight (relative to column) 
     integer , pointer :: landunit (:) ! index into landunit level quantities
     real(r8), pointer :: wtlunit  (:) ! weight (relative to landunit) 
     integer , pointer :: gridcell (:) ! index into gridcell level quantities
     real(r8), pointer :: wtgcell  (:) ! weight (relative to gridcell) 

     ! Non-ED only 
     integer , pointer :: itype    (:) ! patch vegetation 
     integer , pointer :: mxy      (:) ! m index for laixy(i,j,m),etc. (undefined for special landunits)
     logical , pointer :: active   (:) ! true=>do computations on this patch

   contains

     procedure, public :: Init
     procedure, public :: Clean
     
  end type patch_type
  type(patch_type), public, target :: patch  ! patch type data structure
  !------------------------------------------------------------------------

contains
  
  !------------------------------------------------------------------------
  subroutine Init(this, begp, endp)
    !
    ! !ARGUMENTS:
    class(patch_type)   :: this
    integer, intent(in) :: begp,endp
    !
    ! LOCAL VARAIBLES:
    !------------------------------------------------------------------------

    ! The following is set in InitGridCells

    allocate(this%gridcell      (begp:endp)); this%gridcell   (:) = ispval
    allocate(this%wtgcell       (begp:endp)); this%wtgcell    (:) = nan

    allocate(this%landunit      (begp:endp)); this%landunit   (:) = ispval
    allocate(this%wtlunit       (begp:endp)); this%wtlunit    (:) = nan

    allocate(this%column        (begp:endp)); this%column     (:) = ispval
    allocate(this%wtcol         (begp:endp)); this%wtcol      (:) = nan

    allocate(this%mxy           (begp:endp)); this%mxy        (:) = ispval
    allocate(this%active        (begp:endp)); this%active     (:) = .false.

    ! TODO (MV, 10-17-14): The following must be commented out because
    ! currently the logic checking if patch%itype(p) is not equal to noveg
    ! is used in RootBiogeophysMod in zeng2001_rootfr- a filter is not used
    ! in that routine - which would elimate this problem

    allocate(this%itype      (begp:endp)); this%itype      (:) = ispval

  end subroutine Init

  !------------------------------------------------------------------------
  subroutine Clean(this)
    !
    ! !ARGUMENTS:
    class(patch_type) :: this
    !------------------------------------------------------------------------

    deallocate(this%gridcell)
    deallocate(this%wtgcell )
    deallocate(this%landunit)
    deallocate(this%wtlunit )
    deallocate(this%column  )
    deallocate(this%wtcol   )
    deallocate(this%itype   )
    deallocate(this%mxy     )
    deallocate(this%active  )

  end subroutine Clean

end module PatchType
