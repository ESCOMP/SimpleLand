module SurfaceRadiationMod

  !------------------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Calculate solar fluxes absorbed by vegetation and ground surface
  !
  ! !USES:
  use shr_kind_mod      , only : r8 => shr_kind_r8
  use shr_log_mod       , only : errMsg => shr_log_errMsg
  use clm_varctl        , only : use_snicar_frc, use_fates
  use decompMod         , only : bounds_type
  use clm_varcon        , only : namec
  use atm2lndType       , only : atm2lnd_type
  use WaterstateType    , only : waterstate_type
  use CanopyStateType   , only : canopystate_type
  use SurfaceAlbedoType , only : surfalb_type
  use SolarAbsorbedType , only : solarabs_type
  use GridcellType      , only : grc                
  use LandunitType      , only : lun                
  use ColumnType        , only : col                
  use PatchType         , only : patch

  !
  ! !PRIVATE TYPES:
  implicit none
  private

  ! !PUBLIC MEMBER FUNCTIONS:

  !
  ! !PRIVATE DATA:
  type, public :: surfrad_type
     real(r8), pointer, private  :: sfc_frc_aer_patch     (:) ! patch surface forcing of snow with all aerosols (patch) [W/m2]
     real(r8), pointer, private  :: sfc_frc_bc_patch      (:) ! patch surface forcing of snow with BC (patch) [W/m2]
     real(r8), pointer, private  :: sfc_frc_oc_patch      (:) ! patch surface forcing of snow with OC (patch) [W/m2]
     real(r8), pointer, private  :: sfc_frc_dst_patch     (:) ! patch surface forcing of snow with dust (patch) [W/m2]
     real(r8), pointer, private  :: sfc_frc_aer_sno_patch (:) ! patch surface forcing of snow with all aerosols, averaged only when snow is present (patch) [W/m2]
     real(r8), pointer, private  :: sfc_frc_bc_sno_patch  (:) ! patch surface forcing of snow with BC, averaged only when snow is present (patch) [W/m2]
     real(r8), pointer, private  :: sfc_frc_oc_sno_patch  (:) ! patch surface forcing of snow with OC, averaged only when snow is present (patch) [W/m2]
     real(r8), pointer, private  :: sfc_frc_dst_sno_patch (:) ! patch surface forcing of snow with dust, averaged only when snow is present (patch) [W/m2]

     real(r8), pointer, private  :: parveg_ln_patch       (:) ! patch  absorbed par by vegetation at local noon (W/m**2)

     real(r8), pointer, private  :: fsr_sno_vd_patch      (:) ! patch reflected direct beam vis solar radiation from snow (W/m**2)
     real(r8), pointer, private  :: fsr_sno_nd_patch      (:) ! patch reflected direct beam NIR solar radiation from snow (W/m**2)
     real(r8), pointer, private  :: fsr_sno_vi_patch      (:) ! patch reflected diffuse vis solar radiation from snow (W/m**2)
     real(r8), pointer, private  :: fsr_sno_ni_patch      (:) ! patch reflected diffuse NIR solar radiation from snow (W/m**2)

     real(r8), pointer, private  :: fsr_vis_d_patch       (:) ! patch reflected direct beam vis solar radiation (W/m**2)
     real(r8), pointer, private  :: fsr_vis_i_patch       (:) ! patch reflected diffuse vis solar radiation (W/m**2)
     real(r8), pointer, private  :: fsr_vis_d_ln_patch    (:) ! patch reflected direct beam vis solar radiation at local noon (W/m**2)

     real(r8), pointer, private  :: fsds_sno_vd_patch     (:) ! patch incident visible, direct radiation on snow  (for history files)  [W/m2]
     real(r8), pointer, private  :: fsds_sno_nd_patch     (:) ! patch incident near-IR, direct radiation on snow  (for history files)  [W/m2]
     real(r8), pointer, private  :: fsds_sno_vi_patch     (:) ! patch incident visible, diffuse radiation on snow (for history files) [W/m2]
     real(r8), pointer, private  :: fsds_sno_ni_patch     (:) ! patch incident near-IR, diffuse radiation on snow (for history files) [W/m2]

     real(r8), pointer, private  :: fsds_vis_d_patch      (:) ! patch incident direct beam vis solar radiation (W/m**2)
     real(r8), pointer, private  :: fsds_vis_i_patch      (:) ! patch incident diffuse vis solar radiation (W/m**2)
     real(r8), pointer, private  :: fsds_vis_d_ln_patch   (:) ! patch incident direct beam vis solar radiation at local noon (W/m**2)
     real(r8), pointer, private  :: fsds_vis_i_ln_patch   (:) ! patch incident diffuse beam vis solar radiation at local noon (W/m**2)

   contains

     procedure, public  :: Init
     procedure, private :: InitAllocate 
     procedure, private :: InitHistory  
     procedure, private :: InitCold     

  end type surfrad_type

  character(len=*), parameter, private :: sourcefile = &
       __FILE__
  !-----------------------------------------------------------------------

contains

  !------------------------------------------------------------------------
  subroutine Init(this, bounds)

    class(surfrad_type) :: this
    type(bounds_type), intent(in) :: bounds  

    call this%InitAllocate(bounds)
    call this%InitHistory(bounds)
    call this%InitCold(bounds)

  end subroutine Init

  !-----------------------------------------------------------------------
  subroutine InitAllocate(this, bounds)
    !
    ! !USES:
    use shr_infnan_mod, only : nan => shr_infnan_nan, assignment(=)
    !
    ! !ARGUMENTS:
    class(surfrad_type) :: this
    type(bounds_type), intent(in) :: bounds  
    !
    ! !LOCAL VARIABLES:
    integer :: begp, endp
    !---------------------------------------------------------------------

    begp = bounds%begp; endp = bounds%endp

    allocate(this%sfc_frc_aer_patch     (begp:endp))              ; this%sfc_frc_aer_patch     (:)   = nan
    allocate(this%sfc_frc_bc_patch      (begp:endp))              ; this%sfc_frc_bc_patch      (:)   = nan
    allocate(this%sfc_frc_oc_patch      (begp:endp))              ; this%sfc_frc_oc_patch      (:)   = nan
    allocate(this%sfc_frc_dst_patch     (begp:endp))              ; this%sfc_frc_dst_patch     (:)   = nan
    allocate(this%sfc_frc_aer_sno_patch (begp:endp))              ; this%sfc_frc_aer_sno_patch (:)   = nan
    allocate(this%sfc_frc_bc_sno_patch  (begp:endp))              ; this%sfc_frc_bc_sno_patch  (:)   = nan
    allocate(this%sfc_frc_oc_sno_patch  (begp:endp))              ; this%sfc_frc_oc_sno_patch  (:)   = nan
    allocate(this%sfc_frc_dst_sno_patch (begp:endp))              ; this%sfc_frc_dst_sno_patch (:)   = nan

    allocate(this%parveg_ln_patch       (begp:endp))              ; this%parveg_ln_patch       (:)   = nan

    allocate(this%fsr_vis_d_patch       (begp:endp))              ; this%fsr_vis_d_patch       (:)   = nan
    allocate(this%fsr_vis_d_ln_patch    (begp:endp))              ; this%fsr_vis_d_ln_patch    (:)   = nan
    allocate(this%fsr_vis_i_patch       (begp:endp))              ; this%fsr_vis_i_patch       (:)   = nan
    allocate(this%fsr_sno_vd_patch      (begp:endp))              ; this%fsr_sno_vd_patch      (:)   = nan
    allocate(this%fsr_sno_nd_patch      (begp:endp))              ; this%fsr_sno_nd_patch      (:)   = nan
    allocate(this%fsr_sno_vi_patch      (begp:endp))              ; this%fsr_sno_vi_patch      (:)   = nan
    allocate(this%fsr_sno_ni_patch      (begp:endp))              ; this%fsr_sno_ni_patch      (:)   = nan

    allocate(this%fsds_vis_d_patch      (begp:endp))              ; this%fsds_vis_d_patch      (:)   = nan
    allocate(this%fsds_vis_i_patch      (begp:endp))              ; this%fsds_vis_i_patch      (:)   = nan
    allocate(this%fsds_vis_d_ln_patch   (begp:endp))              ; this%fsds_vis_d_ln_patch   (:)   = nan
    allocate(this%fsds_vis_i_ln_patch   (begp:endp))              ; this%fsds_vis_i_ln_patch   (:)   = nan
    allocate(this%fsds_sno_vd_patch     (begp:endp))              ; this%fsds_sno_vd_patch     (:)   = nan
    allocate(this%fsds_sno_nd_patch     (begp:endp))              ; this%fsds_sno_nd_patch     (:)   = nan
    allocate(this%fsds_sno_vi_patch     (begp:endp))              ; this%fsds_sno_vi_patch     (:)   = nan
    allocate(this%fsds_sno_ni_patch     (begp:endp))              ; this%fsds_sno_ni_patch     (:)   = nan

  end subroutine InitAllocate

  !-----------------------------------------------------------------------
  subroutine InitHistory(this, bounds)
    !
    ! History fields initialization
    !
    ! !USES:
    use shr_infnan_mod, only : nan => shr_infnan_nan, assignment(=)
    use clm_varcon    , only : spval
    use histFileMod   , only : hist_addfld1d, hist_addfld2d
    !
    ! !ARGUMENTS:
    class(surfrad_type) :: this
    type(bounds_type), intent(in) :: bounds  
    !
    ! !LOCAL VARIABLES:
    integer :: begp, endp
    integer :: begc, endc
    real(r8), pointer :: data2dptr(:,:) ! temp. pointers for slicing larger arrays
    !---------------------------------------------------------------------

    begp = bounds%begp; endp = bounds%endp
    begc = bounds%begc; endc = bounds%endc

    if (use_snicar_frc) then
       this%sfc_frc_aer_patch(begp:endp) = spval
       call hist_addfld1d (fname='SNOAERFRCL', units='W/m^2', &
            avgflag='A', long_name='surface forcing of all aerosols in snow (land) ', &
            ptr_patch=this%sfc_frc_aer_patch, set_urb=spval, default='inactive')

       this%sfc_frc_aer_sno_patch(begp:endp) = spval
       call hist_addfld1d (fname='SNOAERFRC2L', units='W/m^2', &
            avgflag='A', long_name='surface forcing of all aerosols in snow, averaged only when snow is present (land)', &
            ptr_patch=this%sfc_frc_aer_sno_patch, set_urb=spval, default='inactive')

       this%sfc_frc_bc_patch(begp:endp) = spval
       call hist_addfld1d (fname='SNOBCFRCL', units='W/m^2', &
            avgflag='A', long_name='surface forcing of BC in snow (land) ', &
            ptr_patch=this%sfc_frc_bc_patch, set_urb=spval, default='inactive')

       this%sfc_frc_bc_sno_patch(begp:endp) = spval
       call hist_addfld1d (fname='SNOBCFRC2L', units='W/m^2', &
            avgflag='A', long_name='surface forcing of BC in snow, averaged only when snow is present (land)', &
            ptr_patch=this%sfc_frc_bc_sno_patch, set_urb=spval, default='inactive')

       this%sfc_frc_oc_patch(begp:endp) = spval
       call hist_addfld1d (fname='SNOOCFRCL', units='W/m^2', &
            avgflag='A', long_name='surface forcing of OC in snow (land) ', &
            ptr_patch=this%sfc_frc_oc_patch, set_urb=spval, default='inactive')

       this%sfc_frc_oc_sno_patch(begp:endp) = spval
       call hist_addfld1d (fname='SNOOCFRC2L', units='W/m^2', &
            avgflag='A', long_name='surface forcing of OC in snow, averaged only when snow is present (land)', &
            ptr_patch=this%sfc_frc_oc_sno_patch, set_urb=spval, default='inactive')

       this%sfc_frc_dst_patch(begp:endp) = spval
       call hist_addfld1d (fname='SNODSTFRCL', units='W/m^2', &
            avgflag='A', long_name='surface forcing of dust in snow (land) ', &
            ptr_patch=this%sfc_frc_dst_patch, set_urb=spval, default='inactive')

       this%sfc_frc_dst_sno_patch(begp:endp) = spval
       call hist_addfld1d (fname='SNODSTFRC2L', units='W/m^2', &
            avgflag='A', long_name='surface forcing of dust in snow, averaged only when snow is present (land)', &
            ptr_patch=this%sfc_frc_dst_sno_patch, set_urb=spval, default='inactive')
    end if

    this%fsds_vis_d_patch(begp:endp) = spval
    call hist_addfld1d (fname='FSDSVD', units='W/m^2',  &
         avgflag='A', long_name='direct vis incident solar radiation', &
         ptr_patch=this%fsds_vis_d_patch, default='inactive')

    this%fsds_vis_i_patch(begp:endp) = spval
    call hist_addfld1d (fname='FSDSVI', units='W/m^2',  &
         avgflag='A', long_name='diffuse vis incident solar radiation', &
         ptr_patch=this%fsds_vis_i_patch, default='inactive')

    this%fsr_vis_d_patch(begp:endp) = spval
    call hist_addfld1d (fname='FSRVD', units='W/m^2',  &
         avgflag='A', long_name='direct vis reflected solar radiation', &
         ptr_patch=this%fsr_vis_d_patch, c2l_scale_type='urbanf', default='inactive')

    this%fsr_vis_i_patch(begp:endp) = spval
    call hist_addfld1d (fname='FSRVI', units='W/m^2',  &
         avgflag='A', long_name='diffuse vis reflected solar radiation', &
         ptr_patch=this%fsr_vis_i_patch, c2l_scale_type='urbanf', default='inactive')

    this%fsds_vis_d_ln_patch(begp:endp) = spval
    call hist_addfld1d (fname='FSDSVDLN', units='W/m^2',  &
         avgflag='A', long_name='direct vis incident solar radiation at local noon', &
         ptr_patch=this%fsds_vis_d_ln_patch, default='inactive')

    this%fsds_vis_i_ln_patch(begp:endp) = spval
    call hist_addfld1d (fname='FSDSVILN', units='W/m^2',  &
         avgflag='A', long_name='diffuse vis incident solar radiation at local noon', &
         ptr_patch=this%fsds_vis_i_ln_patch, default='inactive')

    this%parveg_ln_patch(begp:endp) = spval
    call hist_addfld1d (fname='PARVEGLN', units='W/m^2',  &
         avgflag='A', long_name='absorbed par by vegetation at local noon', &
         ptr_patch=this%parveg_ln_patch, default='inactive')

    this%fsr_vis_d_ln_patch(begp:endp) = spval
    call hist_addfld1d (fname='FSRVDLN', units='W/m^2',  &
         avgflag='A', long_name='direct vis reflected solar radiation at local noon', &
         ptr_patch=this%fsr_vis_d_ln_patch, c2l_scale_type='urbanf', default='inactive')

    this%fsds_sno_vd_patch(begp:endp) = spval
    call hist_addfld1d (fname='SNOFSDSVD', units='W/m^2',  &
         avgflag='A', long_name='direct vis incident solar radiation on snow', &
         ptr_patch=this%fsds_sno_vd_patch, default='inactive')

    this%fsds_sno_nd_patch(begp:endp) = spval
    call hist_addfld1d (fname='SNOFSDSND', units='W/m^2',  &
         avgflag='A', long_name='direct nir incident solar radiation on snow', &
         ptr_patch=this%fsds_sno_nd_patch, default='inactive')

    this%fsds_sno_vi_patch(begp:endp) = spval
    call hist_addfld1d (fname='SNOFSDSVI', units='W/m^2',  &
         avgflag='A', long_name='diffuse vis incident solar radiation on snow', &
         ptr_patch=this%fsds_sno_vi_patch, default='inactive')

    this%fsds_sno_ni_patch(begp:endp) = spval
    call hist_addfld1d (fname='SNOFSDSNI', units='W/m^2',  &
         avgflag='A', long_name='diffuse nir incident solar radiation on snow', &
         ptr_patch=this%fsds_sno_ni_patch, default='inactive')

    this%fsr_sno_vd_patch(begp:endp) = spval
    call hist_addfld1d (fname='SNOFSRVD', units='W/m^2',  &
         avgflag='A', long_name='direct vis reflected solar radiation from snow', &
         ptr_patch=this%fsr_sno_vd_patch, default='inactive')

    this%fsr_sno_nd_patch(begp:endp) = spval
    call hist_addfld1d (fname='SNOFSRND', units='W/m^2',  &
         avgflag='A', long_name='direct nir reflected solar radiation from snow', &
         ptr_patch=this%fsr_sno_nd_patch, default='inactive')

    this%fsr_sno_vi_patch(begp:endp) = spval
    call hist_addfld1d (fname='SNOFSRVI', units='W/m^2',  &
         avgflag='A', long_name='diffuse vis reflected solar radiation from snow', &
         ptr_patch=this%fsr_sno_vi_patch, default='inactive')

    this%fsr_sno_ni_patch(begp:endp) = spval
    call hist_addfld1d (fname='SNOFSRNI', units='W/m^2',  &
         avgflag='A', long_name='diffuse nir reflected solar radiation from snow', &
         ptr_patch=this%fsr_sno_ni_patch, default='inactive')


  end subroutine InitHistory

  !------------------------------------------------------------------------
  subroutine InitCold(this, bounds)
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    class(surfrad_type) :: this
    type(bounds_type), intent(in) :: bounds  
    !
    ! !LOCAL VARIABLES:
    integer :: p,l
    !-----------------------------------------------------------------------

    ! nothing for now
    
  end subroutine InitCold
  
end module SurfaceRadiationMod
