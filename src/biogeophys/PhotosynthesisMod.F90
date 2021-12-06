module  PhotosynthesisMod

#include "shr_assert.h"

  !------------------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Leaf photosynthesis and stomatal conductance calculation as described by
  ! Bonan et al (2011) JGR, 116, doi:10.1029/2010JG001593 and extended to
  ! a multi-layer canopy
  !
  ! !USES:
  use shr_sys_mod         , only : shr_sys_flush
  use shr_kind_mod        , only : r8 => shr_kind_r8
  use shr_log_mod         , only : errMsg => shr_log_errMsg
  use shr_infnan_mod      , only : nan => shr_infnan_nan, assignment(=)
  use abortutils          , only : endrun
  use clm_varctl          , only : use_cn, use_cndv, use_fates, use_luna, use_hydrstress
  use clm_varctl          , only : iulog
  use clm_varpar          , only : nlevcan, nvegwcs, mxpft
  use clm_varcon          , only : namep, spval
  use decompMod           , only : bounds_type
  use QuadraticMod        , only : quadratic
  use pftconMod           , only : pftcon
  use atm2lndType         , only : atm2lnd_type
  use CanopyStateType     , only : canopystate_type
  use WaterStateType      , only : waterstate_type
  use WaterfluxType       , only : waterflux_type
  use SoilStateType       , only : soilstate_type
  use TemperatureType     , only : temperature_type
  use SolarAbsorbedType   , only : solarabs_type
  use SurfaceAlbedoType   , only : surfalb_type
  use OzoneBaseMod        , only : ozone_base_type
  use LandunitType        , only : lun
  use PatchType           , only : patch
  use GridcellType        , only : grc
  !
  implicit none
  private
  !
  ! !PRIVATE DATA:
  integer, parameter, private :: leafresp_mtd_ryan1991  = 1  ! Ryan 1991 method for lmr25top
  integer, parameter, private :: leafresp_mtd_atkin2015 = 2  ! Atkin 2015 method for lmr25top
  integer, parameter, private :: sun=1     ! index for sunlit
  integer, parameter, private :: sha=2     ! index for shaded
  integer, parameter, private :: xyl=3     ! index for xylem
  integer, parameter, private :: root=4    ! index for root
  integer, parameter, private :: veg=0     ! index for vegetation
  integer, parameter, private :: soil=1    ! index for soil
  integer, parameter, private :: stomatalcond_mtd_bb1987     = 1   ! Ball-Berry 1987 method for photosynthesis
  integer, parameter, private :: stomatalcond_mtd_medlyn2011 = 2   ! Medlyn 2011 method for photosynthesis
  ! !PUBLIC VARIABLES:

  type :: photo_params_type
     real(r8), allocatable, public  :: krmax              (:)
     real(r8), allocatable, private :: kmax               (:,:)
     real(r8), allocatable, private :: psi50              (:,:)
     real(r8), allocatable, private :: ck                 (:,:)
     real(r8), allocatable, public  :: psi_soil_ref       (:)
     real(r8), allocatable, private :: lmr_intercept_atkin(:)
  contains
     procedure, private :: allocParams
  end type photo_params_type
  !
  type(photo_params_type), public, protected :: params_inst  ! params_inst is populated in readParamsMod 

  type, public :: photosyns_type

     logical , pointer, private :: c3flag_patch      (:)   ! patch true if C3 and false if C4
     ! Plant hydraulic stress specific variables
     real(r8), pointer, private :: ac_phs_patch      (:,:,:) ! patch Rubisco-limited gross photosynthesis (umol CO2/m**2/s)
     real(r8), pointer, private :: aj_phs_patch      (:,:,:) ! patch RuBP-limited gross photosynthesis (umol CO2/m**2/s)
     real(r8), pointer, private :: ap_phs_patch      (:,:,:) ! patch product-limited (C3) or CO2-limited (C4) gross photosynthesis (umol CO2/m**2/s)
     real(r8), pointer, private :: ag_phs_patch      (:,:,:) ! patch co-limited gross leaf photosynthesis (umol CO2/m**2/s)
     real(r8), pointer, private :: an_sun_patch      (:,:)   ! patch sunlit net leaf photosynthesis (umol CO2/m**2/s)
     real(r8), pointer, private :: an_sha_patch      (:,:)   ! patch shaded net leaf photosynthesis (umol CO2/m**2/s)
     real(r8), pointer, private :: vcmax_z_phs_patch (:,:,:) ! patch maximum rate of carboxylation (umol co2/m**2/s)
     real(r8), pointer, private :: kp_z_phs_patch    (:,:,:) ! patch initial slope of CO2 response curve (C4 plants)
     real(r8), pointer, private :: tpu_z_phs_patch   (:,:,:) ! patch triose phosphate utilization rate (umol CO2/m**2/s)
     real(r8), pointer, private :: gs_mol_sun_patch  (:,:) ! patch sunlit leaf stomatal conductance (umol H2O/m**2/s)
     real(r8), pointer, private :: gs_mol_sha_patch  (:,:) ! patch shaded leaf stomatal conductance (umol H2O/m**2/s)

     real(r8), pointer, private :: ac_patch          (:,:) ! patch Rubisco-limited gross photosynthesis (umol CO2/m**2/s)
     real(r8), pointer, private :: aj_patch          (:,:) ! patch RuBP-limited gross photosynthesis (umol CO2/m**2/s)
     real(r8), pointer, private :: ap_patch          (:,:) ! patch product-limited (C3) or CO2-limited (C4) gross photosynthesis (umol CO2/m**2/s)
     real(r8), pointer, private :: ag_patch          (:,:) ! patch co-limited gross leaf photosynthesis (umol CO2/m**2/s)
     real(r8), pointer, private :: an_patch          (:,:) ! patch net leaf photosynthesis (umol CO2/m**2/s)
     real(r8), pointer, private :: vcmax_z_patch     (:,:) ! patch maximum rate of carboxylation (umol co2/m**2/s)
     real(r8), pointer, private :: cp_patch          (:)   ! patch CO2 compensation point (Pa)
     real(r8), pointer, private :: kc_patch          (:)   ! patch Michaelis-Menten constant for CO2 (Pa)
     real(r8), pointer, private :: ko_patch          (:)   ! patch Michaelis-Menten constant for O2 (Pa)
     real(r8), pointer, private :: qe_patch          (:)   ! patch quantum efficiency, used only for C4 (mol CO2 / mol photons)
     real(r8), pointer, private :: tpu_z_patch       (:,:) ! patch triose phosphate utilization rate (umol CO2/m**2/s)
     real(r8), pointer, private :: kp_z_patch        (:,:) ! patch initial slope of CO2 response curve (C4 plants)
     real(r8), pointer, private :: theta_cj_patch    (:)   ! patch empirical curvature parameter for ac, aj photosynthesis co-limitation
     real(r8), pointer, private :: bbb_patch         (:)   ! patch Ball-Berry minimum leaf conductance (umol H2O/m**2/s)
     real(r8), pointer, private :: mbb_patch         (:)   ! patch Ball-Berry slope of conductance-photosynthesis relationship
     real(r8), pointer, private :: gs_mol_patch      (:,:) ! patch leaf stomatal conductance       (umol H2O/m**2/s)
     real(r8), pointer, private :: gb_mol_patch      (:)   ! patch leaf boundary layer conductance (umol H2O/m**2/s)
     real(r8), pointer, private :: rh_leaf_patch     (:)   ! patch fractional humidity at leaf surface (dimensionless)

     real(r8), pointer, private :: alphapsnsun_patch (:)   ! patch sunlit 13c fractionation ([])
     real(r8), pointer, private :: alphapsnsha_patch (:)   ! patch shaded 13c fractionation ([])

     real(r8), pointer, public  :: psnsun_patch      (:)   ! patch sunlit leaf photosynthesis     (umol CO2/m**2/s)
     real(r8), pointer, public  :: psnsha_patch      (:)   ! patch shaded leaf photosynthesis     (umol CO2/m**2/s)

     real(r8), pointer, private :: psnsun_z_patch    (:,:) ! patch canopy layer: sunlit leaf photosynthesis   (umol CO2/m**2/s)
     real(r8), pointer, private :: psnsha_z_patch    (:,:) ! patch canopy layer: shaded leaf photosynthesis   (umol CO2/m**2/s)
     real(r8), pointer, private :: psnsun_wc_patch   (:)   ! patch Rubsico-limited sunlit leaf photosynthesis (umol CO2/m**2/s)
     real(r8), pointer, private :: psnsha_wc_patch   (:)   ! patch Rubsico-limited shaded leaf photosynthesis (umol CO2/m**2/s)
     real(r8), pointer, private :: psnsun_wj_patch   (:)   ! patch RuBP-limited sunlit leaf photosynthesis    (umol CO2/m**2/s)
     real(r8), pointer, private :: psnsha_wj_patch   (:)   ! patch RuBP-limited shaded leaf photosynthesis    (umol CO2/m**2/s)
     real(r8), pointer, private :: psnsun_wp_patch   (:)   ! patch product-limited sunlit leaf photosynthesis (umol CO2/m**2/s)
     real(r8), pointer, private :: psnsha_wp_patch   (:)   ! patch product-limited shaded leaf photosynthesis (umol CO2/m**2/s)

     real(r8), pointer, public  :: fpsn_patch        (:)   ! patch photosynthesis                 (umol CO2/m**2 ground/s)
     real(r8), pointer, private :: fpsn_wc_patch     (:)   ! patch Rubisco-limited photosynthesis (umol CO2/m**2 ground/s)
     real(r8), pointer, private :: fpsn_wj_patch     (:)   ! patch RuBP-limited photosynthesis    (umol CO2/m**2 ground/s)
     real(r8), pointer, private :: fpsn_wp_patch     (:)   ! patch product-limited photosynthesis (umol CO2/m**2 ground/s)

     real(r8), pointer, public  :: lnca_patch        (:)   ! top leaf layer leaf N concentration (gN leaf/m^2)

     real(r8), pointer, public  :: lmrsun_patch      (:)   ! patch sunlit leaf maintenance respiration rate               (umol CO2/m**2/s)
     real(r8), pointer, public  :: lmrsha_patch      (:)   ! patch shaded leaf maintenance respiration rate               (umol CO2/m**2/s)
     real(r8), pointer, private :: lmrsun_z_patch    (:,:) ! patch canopy layer: sunlit leaf maintenance respiration rate (umol CO2/m**2/s)
     real(r8), pointer, private :: lmrsha_z_patch    (:,:) ! patch canopy layer: shaded leaf maintenance respiration rate (umol CO2/m**2/s)

     real(r8), pointer, public  :: cisun_z_patch     (:,:) ! patch intracellular sunlit leaf CO2 (Pa)
     real(r8), pointer, public  :: cisha_z_patch     (:,:) ! patch intracellular shaded leaf CO2 (Pa)

     real(r8), pointer, private :: rssun_z_patch     (:,:) ! patch canopy layer: sunlit leaf stomatal resistance (s/m)
     real(r8), pointer, private :: rssha_z_patch     (:,:) ! patch canopy layer: shaded leaf stomatal resistance (s/m)
     real(r8), pointer, public  :: rssun_patch       (:)   ! patch sunlit stomatal resistance (s/m)
     real(r8), pointer, public  :: rssha_patch       (:)   ! patch shaded stomatal resistance (s/m)
     real(r8), pointer, public  :: luvcmax25top_patch (:)   ! vcmax25 !     (umol/m2/s)
     real(r8), pointer, public  :: lujmax25top_patch  (:)   ! vcmax25 (umol/m2/s)
     real(r8), pointer, public  :: lutpu25top_patch   (:)   ! vcmax25 (umol/m2/s)
!!


     ! LUNA specific variables
     real(r8), pointer, public  :: vcmx25_z_patch    (:,:) ! patch  leaf Vc,max25 (umol CO2/m**2/s) for canopy layer 
     real(r8), pointer, public  :: jmx25_z_patch     (:,:) ! patch  leaf Jmax25 (umol electron/m**2/s) for canopy layer 
     real(r8), pointer, public  :: pnlc_z_patch      (:,:) ! patch proportion of leaf nitrogen allocated for light capture for canopy layer
     real(r8), pointer, public  :: enzs_z_patch      (:,:) ! enzyme decay status 1.0-fully active; 0-all decayed during stress
     real(r8), pointer, public  :: fpsn24_patch      (:)   ! 24 hour mean patch photosynthesis (umol CO2/m**2 ground/day)

     ! Logical switches for different options
     logical, public  :: rootstem_acc                      ! Respiratory acclimation for roots and stems
     logical, private :: light_inhibit                     ! If light should inhibit respiration
     integer, private :: leafresp_method                   ! leaf maintencence respiration at 25C for canopy top method to use
     integer, private :: stomatalcond_mtd                  ! Stomatal conduction method type
     logical, private :: modifyphoto_and_lmr_forcrop       ! Modify photosynthesis and LMR for crop
   contains

     ! Public procedures
     procedure, public  :: Init
     procedure, public  :: Restart
     procedure, public  :: ReadNML
     procedure, public  :: ReadParams

     ! Private procedures
     procedure, private :: InitAllocate
     procedure, private :: InitHistory
     procedure, private :: InitCold

  end type photosyns_type

  character(len=*), parameter, private :: sourcefile = &
       __FILE__
  !------------------------------------------------------------------------

contains

  !------------------------------------------------------------------------
  subroutine Init(this, bounds)

    class(photosyns_type) :: this
    type(bounds_type), intent(in) :: bounds

    call this%InitAllocate (bounds)
    call this%InitHistory  (bounds)
    call this%InitCold     (bounds)

  end subroutine Init

  !------------------------------------------------------------------------
  subroutine InitAllocate(this, bounds)
    !
    ! !ARGUMENTS:
    class(photosyns_type) :: this
    type(bounds_type), intent(in) :: bounds
    !
    ! !LOCAL VARIABLES:
    integer :: begp, endp
    integer :: begc, endc
    !------------------------------------------------------------------------

    begp = bounds%begp; endp= bounds%endp
    begc = bounds%begc; endc= bounds%endc

    allocate(this%c3flag_patch      (begp:endp))             ; this%c3flag_patch      (:)     =.false.
    allocate(this%ac_phs_patch      (begp:endp,2,1:nlevcan)) ; this%ac_phs_patch      (:,:,:) = nan
    allocate(this%aj_phs_patch      (begp:endp,2,1:nlevcan)) ; this%aj_phs_patch      (:,:,:) = nan
    allocate(this%ap_phs_patch      (begp:endp,2,1:nlevcan)) ; this%ap_phs_patch      (:,:,:) = nan
    allocate(this%ag_phs_patch      (begp:endp,2,1:nlevcan)) ; this%ag_phs_patch      (:,:,:) = nan
    allocate(this%an_sun_patch      (begp:endp,1:nlevcan))   ; this%an_sun_patch      (:,:)   = nan
    allocate(this%an_sha_patch      (begp:endp,1:nlevcan))   ; this%an_sha_patch      (:,:)   = nan
    allocate(this%vcmax_z_phs_patch (begp:endp,2,1:nlevcan)) ; this%vcmax_z_phs_patch (:,:,:) = nan
    allocate(this%tpu_z_phs_patch   (begp:endp,2,1:nlevcan)) ; this%tpu_z_phs_patch   (:,:,:) = nan
    allocate(this%kp_z_phs_patch    (begp:endp,2,1:nlevcan)) ; this%kp_z_phs_patch    (:,:,:) = nan
    allocate(this%gs_mol_sun_patch  (begp:endp,1:nlevcan))   ; this%gs_mol_sun_patch  (:,:)   = nan
    allocate(this%gs_mol_sha_patch  (begp:endp,1:nlevcan))   ; this%gs_mol_sha_patch  (:,:)   = nan
    allocate(this%ac_patch          (begp:endp,1:nlevcan)) ; this%ac_patch          (:,:) = nan
    allocate(this%aj_patch          (begp:endp,1:nlevcan)) ; this%aj_patch          (:,:) = nan
    allocate(this%ap_patch          (begp:endp,1:nlevcan)) ; this%ap_patch          (:,:) = nan
    allocate(this%ag_patch          (begp:endp,1:nlevcan)) ; this%ag_patch          (:,:) = nan
    allocate(this%an_patch          (begp:endp,1:nlevcan)) ; this%an_patch          (:,:) = nan
    allocate(this%vcmax_z_patch     (begp:endp,1:nlevcan)) ; this%vcmax_z_patch     (:,:) = nan
    allocate(this%tpu_z_patch       (begp:endp,1:nlevcan)) ; this%tpu_z_patch       (:,:) = nan
    allocate(this%kp_z_patch        (begp:endp,1:nlevcan)) ; this%kp_z_patch        (:,:) = nan
    allocate(this%gs_mol_patch      (begp:endp,1:nlevcan)) ; this%gs_mol_patch      (:,:) = nan
    allocate(this%cp_patch          (begp:endp))           ; this%cp_patch          (:)   = nan
    allocate(this%kc_patch          (begp:endp))           ; this%kc_patch          (:)   = nan
    allocate(this%ko_patch          (begp:endp))           ; this%ko_patch          (:)   = nan
    allocate(this%qe_patch          (begp:endp))           ; this%qe_patch          (:)   = nan
    allocate(this%theta_cj_patch    (begp:endp))           ; this%theta_cj_patch    (:)   = nan
    allocate(this%bbb_patch         (begp:endp))           ; this%bbb_patch         (:)   = nan
    allocate(this%mbb_patch         (begp:endp))           ; this%mbb_patch         (:)   = nan
    allocate(this%gb_mol_patch      (begp:endp))           ; this%gb_mol_patch      (:)   = nan
    allocate(this%rh_leaf_patch     (begp:endp))           ; this%rh_leaf_patch     (:)   = nan

    allocate(this%psnsun_patch      (begp:endp))           ; this%psnsun_patch      (:)   = nan
    allocate(this%psnsha_patch      (begp:endp))           ; this%psnsha_patch      (:)   = nan

    allocate(this%psnsun_z_patch    (begp:endp,1:nlevcan)) ; this%psnsun_z_patch    (:,:) = nan
    allocate(this%psnsha_z_patch    (begp:endp,1:nlevcan)) ; this%psnsha_z_patch    (:,:) = nan
    allocate(this%psnsun_wc_patch   (begp:endp))           ; this%psnsun_wc_patch   (:)   = nan
    allocate(this%psnsha_wc_patch   (begp:endp))           ; this%psnsha_wc_patch   (:)   = nan
    allocate(this%psnsun_wj_patch   (begp:endp))           ; this%psnsun_wj_patch   (:)   = nan
    allocate(this%psnsha_wj_patch   (begp:endp))           ; this%psnsha_wj_patch   (:)   = nan
    allocate(this%psnsun_wp_patch   (begp:endp))           ; this%psnsun_wp_patch   (:)   = nan
    allocate(this%psnsha_wp_patch   (begp:endp))           ; this%psnsha_wp_patch   (:)   = nan
    allocate(this%fpsn_patch        (begp:endp))           ; this%fpsn_patch        (:)   = nan
    allocate(this%fpsn_wc_patch     (begp:endp))           ; this%fpsn_wc_patch     (:)   = nan
    allocate(this%fpsn_wj_patch     (begp:endp))           ; this%fpsn_wj_patch     (:)   = nan
    allocate(this%fpsn_wp_patch     (begp:endp))           ; this%fpsn_wp_patch     (:)   = nan
    
    allocate(this%lnca_patch        (begp:endp))           ; this%lnca_patch        (:)   = nan

    allocate(this%lmrsun_z_patch    (begp:endp,1:nlevcan)) ; this%lmrsun_z_patch    (:,:) = nan
    allocate(this%lmrsha_z_patch    (begp:endp,1:nlevcan)) ; this%lmrsha_z_patch    (:,:) = nan
    allocate(this%lmrsun_patch      (begp:endp))           ; this%lmrsun_patch      (:)   = nan
    allocate(this%lmrsha_patch      (begp:endp))           ; this%lmrsha_patch      (:)   = nan

    allocate(this%alphapsnsun_patch (begp:endp))           ; this%alphapsnsun_patch (:)   = nan
    allocate(this%alphapsnsha_patch (begp:endp))           ; this%alphapsnsha_patch (:)   = nan

    allocate(this%cisun_z_patch     (begp:endp,1:nlevcan)) ; this%cisun_z_patch     (:,:) = nan
    allocate(this%cisha_z_patch     (begp:endp,1:nlevcan)) ; this%cisha_z_patch     (:,:) = nan

    allocate(this%rssun_z_patch     (begp:endp,1:nlevcan)) ; this%rssun_z_patch     (:,:) = nan
    allocate(this%rssha_z_patch     (begp:endp,1:nlevcan)) ; this%rssha_z_patch     (:,:) = nan
    allocate(this%rssun_patch       (begp:endp))           ; this%rssun_patch       (:)   = nan
    allocate(this%rssha_patch       (begp:endp))           ; this%rssha_patch       (:)   = nan
    allocate(this%luvcmax25top_patch(begp:endp))           ; this%luvcmax25top_patch(:) = nan
    allocate(this%lujmax25top_patch (begp:endp))           ; this%lujmax25top_patch(:)  = nan
    allocate(this%lutpu25top_patch  (begp:endp))           ; this%lutpu25top_patch(:)   = nan
!!
!    allocate(this%psncanopy_patch   (begp:endp))           ; this%psncanopy_patch   (:)   = nan
!    allocate(this%lmrcanopy_patch   (begp:endp))           ; this%lmrcanopy_patch   (:)   = nan
    if(use_luna)then
      ! NOTE(bja, 2015-09) because these variables are only allocated
      ! when luna is turned on, they can not be placed into associate
      ! statements.
      allocate(this%vcmx25_z_patch  (begp:endp,1:nlevcan)) ; this%vcmx25_z_patch    (:,:) = 30._r8
      allocate(this%jmx25_z_patch   (begp:endp,1:nlevcan)) ; this%jmx25_z_patch     (:,:) = 60._r8 
      allocate(this%pnlc_z_patch    (begp:endp,1:nlevcan)) ; this%pnlc_z_patch      (:,:) = 0.01_r8
      allocate(this%fpsn24_patch    (begp:endp))           ; this%fpsn24_patch      (:)   = nan
      allocate(this%enzs_z_patch    (begp:endp,1:nlevcan)) ; this%enzs_z_patch      (:,:) = 1._r8
    endif

  end subroutine InitAllocate

  !-----------------------------------------------------------------------
  subroutine InitHistory(this, bounds)
    !
    ! !USES:
    use histFileMod   , only: hist_addfld1d, hist_addfld2d
    !
    ! !ARGUMENTS:
    class(photosyns_type) :: this
    type(bounds_type), intent(in) :: bounds
    real(r8), pointer  :: ptr_1d(:)  ! pointer to 1d patch array
    !
    ! !LOCAL VARIABLES:
    integer :: begp, endp
    !---------------------------------------------------------------------

    begp = bounds%begp; endp= bounds%endp

    this%rh_leaf_patch(begp:endp) = spval
    call hist_addfld1d (fname='RH_LEAF', units='fraction', &
         avgflag='A', long_name='fractional humidity at leaf surface', &
         ptr_patch=this%rh_leaf_patch, set_spec=spval, default='inactive')
    this%lnca_patch(begp:endp) = spval
    call hist_addfld1d (fname='LNC', units='gN leaf/m^2', &
         avgflag='A', long_name='leaf N concentration', &
         ptr_patch=this%lnca_patch, set_spec=spval)

    ! Don't output photosynthesis variables when FATES is on as they aren't calculated
    if (.not. use_fates) then
       this%fpsn_patch(begp:endp) = spval
       call hist_addfld1d (fname='FPSN', units='umol/m2s',  &
            avgflag='A', long_name='photosynthesis', &
            ptr_patch=this%fpsn_patch, set_lake=0._r8, set_urb=0._r8)

       ! Don't by default output this rate limiting step as only makes sense if you are outputing
       ! the others each time-step
       this%fpsn_wc_patch(begp:endp) = spval
       call hist_addfld1d (fname='FPSN_WC', units='umol/m2s',  &
            avgflag='I', long_name='Rubisco-limited photosynthesis', &
            ptr_patch=this%fpsn_wc_patch, set_lake=0._r8, set_urb=0._r8, &
            default='inactive')

       ! Don't by default output this rate limiting step as only makes sense if you are outputing
       ! the others each time-step
       this%fpsn_wj_patch(begp:endp) = spval
       call hist_addfld1d (fname='FPSN_WJ', units='umol/m2s',  &
            avgflag='I', long_name='RuBP-limited photosynthesis', &
            ptr_patch=this%fpsn_wj_patch, set_lake=0._r8, set_urb=0._r8, &
            default='inactive')

       ! Don't by default output this rate limiting step as only makes sense if you are outputing
       ! the others each time-step
       this%fpsn_wp_patch(begp:endp) = spval
       call hist_addfld1d (fname='FPSN_WP', units='umol/m2s',  &
            avgflag='I', long_name='Product-limited photosynthesis', &
            ptr_patch=this%fpsn_wp_patch, set_lake=0._r8, set_urb=0._r8, &
            default='inactive')
    end if

    if (use_cn) then
       this%psnsun_patch(begp:endp) = spval
       call hist_addfld1d (fname='PSNSUN', units='umolCO2/m^2/s', &
            avgflag='A', long_name='sunlit leaf photosynthesis', &
            ptr_patch=this%psnsun_patch)

       this%psnsha_patch(begp:endp) = spval
       call hist_addfld1d (fname='PSNSHA', units='umolCO2/m^2/s', &
            avgflag='A', long_name='shaded leaf photosynthesis', &
            ptr_patch=this%psnsha_patch)
    end if

    this%rssun_patch(begp:endp) = spval
    call hist_addfld1d (fname='RSSUN', units='s/m',  &
         avgflag='M', long_name='sunlit leaf stomatal resistance', &
         ptr_patch=this%rssun_patch, set_lake=spval, set_urb=spval)

    this%rssha_patch(begp:endp) = spval
    call hist_addfld1d (fname='RSSHA', units='s/m',  &
         avgflag='M', long_name='shaded leaf stomatal resistance', &
         ptr_patch=this%rssha_patch, set_lake=spval, set_urb=spval)

    this%gs_mol_sun_patch(begp:endp,:) = spval
    this%gs_mol_sha_patch(begp:endp,:) = spval
    if (nlevcan>1) then 
       call hist_addfld2d (fname='GSSUN', units='umol H20/m2/s', type2d='nlevcan', &
          avgflag='A', long_name='sunlit leaf stomatal conductance', &
          ptr_patch=this%gs_mol_sun_patch, set_lake=spval, set_urb=spval)

       call hist_addfld2d (fname='GSSHA', units='umol H20/m2/s', type2d='nlevcan', &
          avgflag='A', long_name='shaded leaf stomatal conductance', &
          ptr_patch=this%gs_mol_sha_patch, set_lake=spval, set_urb=spval)
    else
       ptr_1d => this%gs_mol_sun_patch(begp:endp,1)
       call hist_addfld1d (fname='GSSUN', units='umol H20/m2/s', &
          avgflag='A', long_name='sunlit leaf stomatal conductance', &
          ptr_patch=ptr_1d)

       ptr_1d => this%gs_mol_sha_patch(begp:endp,1)
       call hist_addfld1d (fname='GSSHA', units='umol H20/m2/s', &
          avgflag='A', long_name='shaded leaf stomatal conductance', &
          ptr_patch=ptr_1d)

    endif

    if(use_luna)then  
       if(nlevcan>1)then
         call hist_addfld2d (fname='Vcmx25Z', units='umol/m2/s', type2d='nlevcan', &
            avgflag='A', long_name='canopy profile of vcmax25 predicted by LUNA model', &
            ptr_patch=this%vcmx25_z_patch)
 
         call hist_addfld2d (fname='Jmx25Z', units='umol/m2/s', type2d='nlevcan', &
            avgflag='A', long_name='canopy profile of  vcmax25 predicted by LUNA model', &
            ptr_patch=this%jmx25_z_patch)

         call hist_addfld2d (fname='PNLCZ', units='unitless', type2d='nlevcan', &
            avgflag='A', long_name='Proportion of nitrogen allocated for light capture', &
            ptr_patch=this%pnlc_z_patch,default='inactive')
       else
         ptr_1d => this%vcmx25_z_patch(:,1)
         call hist_addfld1d (fname='Vcmx25Z', units='umol/m2/s',&
            avgflag='A', long_name='canopy profile of vcmax25 predicted by LUNA model', &
            ptr_patch=ptr_1d)
         ptr_1d => this%jmx25_z_patch(:,1)
         call hist_addfld1d (fname='Jmx25Z', units='umol/m2/s',&
            avgflag='A', long_name='canopy profile of  vcmax25 predicted by LUNA model', &
            ptr_patch=ptr_1d)
         ptr_1d => this%pnlc_z_patch(:,1)
         call hist_addfld1d (fname='PNLCZ', units='unitless', &
            avgflag='A', long_name='Proportion of nitrogen allocated for light capture', &
            ptr_patch=ptr_1d,default='inactive')

         this%luvcmax25top_patch(begp:endp) = spval
         call hist_addfld1d (fname='VCMX25T', units='umol/m2/s',  &
            avgflag='M', long_name='canopy profile of vcmax25', &
            ptr_patch=this%luvcmax25top_patch, set_lake=spval, set_urb=spval)

         this%lujmax25top_patch(begp:endp) = spval
         call hist_addfld1d (fname='JMX25T', units='umol/m2/s',  &
            avgflag='M', long_name='canopy profile of jmax', &
            ptr_patch=this%lujmax25top_patch, set_lake=spval, set_urb=spval)

            this%lutpu25top_patch(begp:endp) = spval
            call hist_addfld1d (fname='TPU25T', units='umol/m2/s',  &
            avgflag='M', long_name='canopy profile of tpu', &
            ptr_patch=this%lutpu25top_patch, set_lake=spval, set_urb=spval)

       endif
       this%fpsn24_patch = spval 
       call hist_addfld1d (fname='FPSN24', units='umol CO2/m**2 ground/day',&
           avgflag='A', long_name='24 hour accumulative patch photosynthesis starting from mid-night', &
           ptr_patch=this%fpsn24_patch, default='inactive')
   
    endif

  end subroutine InitHistory

  !-----------------------------------------------------------------------
  subroutine InitCold(this, bounds)
    !
    ! !ARGUMENTS:
    class(photosyns_type) :: this
    type(bounds_type), intent(in) :: bounds
    !
    ! !LOCAL VARIABLES:
    integer :: p,l                        ! indices
    !-----------------------------------------------------------------------

    do p = bounds%begp,bounds%endp
       l = patch%landunit(p)

       this%alphapsnsun_patch(p) = spval
       this%alphapsnsha_patch(p) = spval

       if (lun%ifspecial(l)) then
          this%psnsun_patch(p) = 0._r8
          this%psnsha_patch(p) = 0._r8
       end if
    end do

  end subroutine InitCold

  !-----------------------------------------------------------------------
  subroutine allocParams ( this )
    !
    implicit none

    ! !ARGUMENTS:
    class(photo_params_type) :: this
    !
    ! !LOCAL VARIABLES:
    character(len=32)  :: subname = 'allocParams'
    !-----------------------------------------------------------------------

    ! allocate parameters

    allocate( this%krmax       (0:mxpft) )          ; this%krmax(:)        = nan
    allocate( this%kmax        (0:mxpft,nvegwcs) )  ; this%kmax(:,:)       = nan
    allocate( this%psi50       (0:mxpft,nvegwcs) )  ; this%psi50(:,:)      = nan
    allocate( this%ck          (0:mxpft,nvegwcs) )  ; this%ck(:,:)         = nan
    allocate( this%psi_soil_ref(0:mxpft) )          ; this%psi_soil_ref(:) = nan

    if ( use_hydrstress .and. nvegwcs /= 4 )then
       call endrun(msg='Error:: the Plant Hydraulics Stress methodology is for the spacA function is hardcoded for nvegwcs==4' &
                   //errMsg(__FILE__, __LINE__))
    end if

  end subroutine allocParams

  !-----------------------------------------------------------------------
  subroutine readParams ( this, ncid )
    !
    ! !USES:
    use ncdio_pio , only : file_desc_t,ncd_io
    implicit none

    ! !ARGUMENTS:
    class(photosyns_type) :: this
    type(file_desc_t),intent(inout) :: ncid   ! pio netCDF file id
    !
    ! !LOCAL VARIABLES:
    character(len=32)  :: subname = 'readParams'
    character(len=100) :: errCode = '-Error reading in parameters file:'
    logical            :: readv ! has variable been read in or not
    real(r8)           :: temp1d(0:mxpft) ! temporary to read in parameter
    real(r8)           :: temp2d(0:mxpft,nvegwcs) ! temporary to read in parameter
    character(len=100) :: tString ! temp. var for reading
    !-----------------------------------------------------------------------

    ! read in parameters


    call params_inst%allocParams()

    tString = "krmax"
    call ncd_io(varname=trim(tString),data=temp1d, flag='read', ncid=ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))
    params_inst%krmax=temp1d
    tString = "psi_soil_ref"
    call ncd_io(varname=trim(tString),data=temp1d, flag='read', ncid=ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))
    params_inst%psi_soil_ref=temp1d
    tString = "lmr_intercept_atkin"
    call ncd_io(varname=trim(tString),data=temp1d, flag='read', ncid=ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))
    params_inst%lmr_intercept_atkin=temp1d
    tString = "kmax"
    call ncd_io(varname=trim(tString),data=temp2d, flag='read', ncid=ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))
    params_inst%kmax=temp2d
    tString = "psi50"
    call ncd_io(varname=trim(tString),data=temp2d, flag='read', ncid=ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))
    params_inst%psi50=temp2d
    tString = "ck"
    call ncd_io(varname=trim(tString),data=temp2d, flag='read', ncid=ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))
    params_inst%ck=temp2d

  end subroutine readParams


  !------------------------------------------------------------------------
  subroutine ReadNML(this, NLFilename)
    !
    ! !DESCRIPTION:
    ! Read the namelist for Photosynthesis
    !
    ! !USES:
    use fileutils      , only : getavu, relavu, opnfil
    use shr_nl_mod     , only : shr_nl_find_group_name
    use spmdMod        , only : masterproc, mpicom
    use shr_mpi_mod    , only : shr_mpi_bcast
    use clm_varctl     , only : iulog
    !
    ! !ARGUMENTS:
    class(photosyns_type) :: this
    character(len=*), intent(IN) :: NLFilename ! Namelist filename
    !
    ! !LOCAL VARIABLES:
    integer :: ierr                 ! error code
    integer :: unitn                ! unit for namelist file

    character(len=*), parameter :: subname = 'Photosyn::ReadNML'
    character(len=*), parameter :: nmlname = 'photosyns_inparm'
    logical :: rootstem_acc    = .false.                     ! Respiratory acclimation for roots and stems
    logical :: light_inhibit   = .false.                     ! If light should inhibit respiration
    integer :: leafresp_method = leafresp_mtd_ryan1991       ! leaf maintencence respiration at 25C for canopy top method to use
    logical :: modifyphoto_and_lmr_forcrop = .false.            ! Modify photosynthesis and LMR for crop
    character(len=50) :: stomatalcond_method = 'Ball-Berry1987' ! Photosynthesis method string
    !-----------------------------------------------------------------------

    namelist /photosyns_inparm/ leafresp_method, light_inhibit, &
              rootstem_acc, stomatalcond_method, modifyphoto_and_lmr_forcrop

    ! Initialize options to default values, in case they are not specified in
    ! the namelist

    if (masterproc) then
       unitn = getavu()
       write(iulog,*) 'Read in '//nmlname//'  namelist'
       call opnfil (NLFilename, unitn, 'F')
       call shr_nl_find_group_name(unitn, nmlname, status=ierr)
       if (ierr == 0) then
          read(unitn, nml=photosyns_inparm, iostat=ierr)
          if (ierr /= 0) then
             call endrun(msg="ERROR reading "//nmlname//"namelist"//errmsg(sourcefile, __LINE__))
          end if
       else
          call endrun(msg="ERROR could NOT find "//nmlname//"namelist"//errmsg(sourcefile, __LINE__))
       end if
       call relavu( unitn )
       this%rootstem_acc    = rootstem_acc
       this%leafresp_method = leafresp_method
       this%light_inhibit   = light_inhibit
       this%modifyphoto_and_lmr_forcrop = modifyphoto_and_lmr_forcrop
       if (      trim(stomatalcond_method) == 'Ball-Berry1987' ) then
          this%stomatalcond_mtd = stomatalcond_mtd_bb1987
       else if ( trim(stomatalcond_method) == 'Medlyn2011'     ) then
          this%stomatalcond_mtd = stomatalcond_mtd_medlyn2011
       else
          call endrun(msg="ERROR bad value for stomtalcond_method in "//nmlname//"namelist"//errmsg(sourcefile, __LINE__))
       end if
    end if

    call shr_mpi_bcast (this%rootstem_acc   , mpicom)
    call shr_mpi_bcast (this%leafresp_method, mpicom)
    call shr_mpi_bcast (this%light_inhibit  , mpicom)
    call shr_mpi_bcast (this%stomatalcond_mtd, mpicom)
    call shr_mpi_bcast (this%modifyphoto_and_lmr_forcrop, mpicom)

    if (masterproc) then
       write(iulog,*) ' '
       write(iulog,*) nmlname//' settings:'
       write(iulog,nml=photosyns_inparm)
       write(iulog,*) ' '
    end if

  end subroutine ReadNML

  !------------------------------------------------------------------------
  subroutine Restart(this, bounds, ncid, flag)
    !
    ! !USES:
    use ncdio_pio  , only : file_desc_t, ncd_defvar, ncd_io, ncd_double, ncd_int, ncd_inqvdlen
    use restUtilMod
    !
    ! !ARGUMENTS:
    class(photosyns_type) :: this
    type(bounds_type), intent(in)    :: bounds
    type(file_desc_t), intent(inout) :: ncid   ! netcdf id
    character(len=*) , intent(in)    :: flag   ! 'read' or 'write'
    !
    ! !LOCAL VARIABLES:
    integer :: j,c ! indices
    logical :: readvar      ! determine if variable is on initial file
    !-----------------------------------------------------------------------

    call restartvar(ncid=ncid, flag=flag, varname='GSSUN', xtype=ncd_double,  &
         dim1name='pft', dim2name='levcan', switchdim=.true., &
         long_name='sunlit leaf stomatal conductance', units='umol H20/m2/s', &
         interpinic_flag='interp', readvar=readvar, data=this%gs_mol_sun_patch)
    
    call restartvar(ncid=ncid, flag=flag, varname='GSSHA', xtype=ncd_double,  &
         dim1name='pft', dim2name='levcan', switchdim=.true., &
         long_name='shaded leaf stomatal conductance', units='umol H20/m2/s', &
         interpinic_flag='interp', readvar=readvar, data=this%gs_mol_sha_patch)
    
    call restartvar(ncid=ncid, flag=flag, varname='lnca', xtype=ncd_double,  &
       dim1name='pft', long_name='leaf N concentration', units='gN leaf/m^2', &
       interpinic_flag='interp', readvar=readvar, data=this%lnca_patch)

    if(use_luna) then
      call restartvar(ncid=ncid, flag=flag, varname='vcmx25_z', xtype=ncd_double,  &
         dim1name='pft', dim2name='levcan', switchdim=.true., &
         long_name='Maximum carboxylation rate at 25 celcius for canopy layers', units='umol CO2/m**2/s', &
         interpinic_flag='interp', readvar=readvar, data=this%vcmx25_z_patch)
      call restartvar(ncid=ncid, flag=flag, varname='jmx25_z', xtype=ncd_double,  &
         dim1name='pft', dim2name='levcan', switchdim=.true., &
         long_name='Maximum carboxylation rate at 25 celcius for canopy layers', units='umol CO2/m**2/s', &
         interpinic_flag='interp', readvar=readvar, data=this%jmx25_z_patch)
      call restartvar(ncid=ncid, flag=flag, varname='pnlc_z', xtype=ncd_double,  &
         dim1name='pft', dim2name='levcan', switchdim=.true., &
         long_name='proportion of leaf nitrogen allocated for light capture', units='unitless', &
         interpinic_flag='interp', readvar=readvar, data=this%pnlc_z_patch )
      call restartvar(ncid=ncid, flag=flag, varname='enzs_z', xtype=ncd_double,  &
         dim1name='pft', dim2name='levcan', switchdim=.true., &
         long_name='enzyme decay status during stress: 1.0-fully active; 0.0-all decayed', units='unitless', &
         interpinic_flag='interp', readvar=readvar, data=this%enzs_z_patch )
      call restartvar(ncid=ncid, flag=flag, varname='gpp24', xtype=ncd_double,  &
            dim1name='pft', long_name='accumulative gross primary production', units='umol CO2/m**2 ground/day', &
            interpinic_flag='interp', readvar=readvar, data=this%fpsn24_patch)    
   endif
   call restartvar(ncid=ncid, flag=flag, varname='vcmx25t', xtype=ncd_double,  &
         dim1name='pft', long_name='canopy profile of vcmax25', &
         units='umol/m2/s', &
         interpinic_flag='interp', readvar=readvar, data=this%luvcmax25top_patch)    

   call restartvar(ncid=ncid, flag=flag, varname='jmx25t', xtype=ncd_double,  &
         dim1name='pft', long_name='canopy profile of jmax', &
         units='umol/m2/s', &
         interpinic_flag='interp', readvar=readvar, data=this%lujmax25top_patch)    

   call restartvar(ncid=ncid, flag=flag, varname='tpu25t', xtype=ncd_double,  &
         dim1name='pft', long_name='canopy profile of tpu', &
         units='umol/m2/s', &
         interpinic_flag='interp', readvar=readvar, data=this%lutpu25top_patch)    

  end subroutine Restart

end module PhotosynthesisMod
