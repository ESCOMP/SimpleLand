module WaterfluxType

  !------------------------------------------------------------------------------
  ! !DESCRIPTION:
  !
  ! !USES:
  use shr_kind_mod   , only: r8 => shr_kind_r8
  use shr_infnan_mod , only : nan => shr_infnan_nan, assignment(=)
  use clm_varpar     , only : nlevsno, nlevsoi
  use clm_varcon     , only : spval
  use decompMod      , only : bounds_type
  use LandunitType   , only : lun                
  use ColumnType     , only : col                
  use PatchType      , only : patch   
  !
  implicit none
  private
  !
  ! !PUBLIC TYPES:
  type, public :: waterflux_type

     ! water fluxes are in units or mm/s

     real(r8), pointer :: qflx_evap_tot_patch      (:)   ! patch pft_qflx_evap_soi + pft_qflx_evap_veg + qflx_tran_veg

     real(r8), pointer :: qflx_surf_col            (:)   ! col surface runoff (mm H2O /s)
     real(r8), pointer :: qflx_drain_col           (:)   ! col sub-surface runoff (mm H2O /s)
     real(r8), pointer :: qflx_h2osfc_surf_col     (:)   ! col surface water runoff
     real(r8), pointer :: qflx_drain_perched_col   (:)   ! col sub-surface runoff from perched wt (mm H2O /s)
     real(r8), pointer :: qflx_qrgwl_col           (:)   ! col qflx_surf at glaciers, wetlands, lakes
     real(r8), pointer :: qflx_runoff_col          (:)   ! col total runoff (qflx_drain+qflx_surf+qflx_qrgwl) (mm H2O /s)

     real(r8), pointer :: qflx_liq_dynbal_grc      (:)   ! grc liq dynamic land cover change conversion runoff flux

   contains
 
     procedure, public  :: Init
     procedure, private :: InitAllocate 
     procedure, private :: InitCold     

  end type waterflux_type
  !------------------------------------------------------------------------

contains

  !------------------------------------------------------------------------
  subroutine Init(this, bounds)

    class(waterflux_type) :: this
    type(bounds_type), intent(in)    :: bounds  

    call this%InitAllocate(bounds) ! same as "call initAllocate_type(hydro, bounds)"
    call this%InitCold(bounds)

  end subroutine Init

  !------------------------------------------------------------------------
  subroutine InitAllocate(this, bounds)
    !
    ! !DESCRIPTION:
    ! Initialize module data structure
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    class(waterflux_type) :: this
    type(bounds_type), intent(in) :: bounds  
    !
    ! !LOCAL VARIABLES:
    integer :: begp, endp
    integer :: begc, endc
    integer :: begg, endg
    !------------------------------------------------------------------------

    begp = bounds%begp; endp= bounds%endp
    begc = bounds%begc; endc= bounds%endc
    begg = bounds%begg; endg= bounds%endg

    allocate(this%qflx_evap_tot_patch      (begp:endp))              ; this%qflx_evap_tot_patch      (:)   = nan

    allocate(this%qflx_surf_col            (begc:endc))              ; this%qflx_surf_col            (:)   = nan
    allocate(this%qflx_drain_col           (begc:endc))              ; this%qflx_drain_col           (:)   = nan
    allocate(this%qflx_h2osfc_surf_col     (begc:endc))              ; this%qflx_h2osfc_surf_col     (:)   = nan
    allocate(this%qflx_qrgwl_col           (begc:endc))              ; this%qflx_qrgwl_col           (:)   = nan
    allocate(this%qflx_drain_perched_col   (begc:endc))              ; this%qflx_drain_perched_col   (:)   = nan
    allocate(this%qflx_runoff_col          (begc:endc))              ; this%qflx_runoff_col          (:)   = nan

    allocate(this%qflx_liq_dynbal_grc      (begg:endg))              ; this%qflx_liq_dynbal_grc      (:)   = nan

  end subroutine InitAllocate

  !-----------------------------------------------------------------------
  subroutine InitCold(this, bounds)
    !
    ! !ARGUMENTS:
    class(waterflux_type) :: this
    type(bounds_type) , intent(in) :: bounds
    !-----------------------------------------------------------------------

    this%qflx_h2osfc_surf_col(bounds%begc:bounds%endc) = 0._r8
    this%qflx_drain_col(bounds%begc:bounds%endc) = 0._r8
    this%qflx_surf_col(bounds%begc:bounds%endc)  = 0._r8

  end subroutine InitCold

end module WaterfluxType
