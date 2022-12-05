module EnergyFluxType

#include "shr_assert.h"

  !------------------------------------------------------------------------------
  ! Energy flux data structure
  !
  use shr_kind_mod   , only : r8 => shr_kind_r8
  use shr_log_mod    , only : errMsg => shr_log_errMsg
  use clm_varcon     , only : spval
  use decompMod      , only : bounds_type
  use LandunitType   , only : lun                
  use ColumnType     , only : col                
  use PatchType      , only : patch                
  !
  implicit none
  save
  private
  !
  type, public :: energyflux_type

     ! Fluxes
     real(r8), pointer :: eflx_sh_grnd_patch      (:)   ! patch sensible heat flux from ground (W/m**2) [+ to atm]
     real(r8), pointer :: eflx_sh_veg_patch       (:)   ! patch sensible heat flux from leaves (W/m**2) [+ to atm]
     real(r8), pointer :: eflx_sh_snow_patch      (:)   ! patch sensible heat flux from snow (W/m**2) [+ to atm]
     real(r8), pointer :: eflx_sh_soil_patch      (:)   ! patch sensible heat flux from soil  (W/m**2) [+ to atm]
     real(r8), pointer :: eflx_sh_h2osfc_patch    (:)   ! patch sensible heat flux from surface water (W/m**2) [+ to atm]
     real(r8), pointer :: eflx_sh_tot_patch       (:)   ! patch total sensible heat flux (W/m**2) [+ to atm]
     real(r8), pointer :: eflx_sh_precip_conversion_col(:) ! col sensible heat flux from precipitation conversion (W/m**2) [+ to atm]
     real(r8), pointer :: eflx_lh_tot_patch       (:)   ! patch total latent heat flux (W/m**2) [+ to atm]
     real(r8), pointer :: eflx_lh_vegt_patch      (:)   ! patch transpiration heat flux from veg (W/m**2) [+ to atm]
     real(r8), pointer :: eflx_lh_vege_patch      (:)   ! patch evaporation heat flux from veg (W/m**2) [+ to atm]
     real(r8), pointer :: eflx_lh_grnd_patch      (:)   ! patch evaporation heat flux from ground (W/m**2) [+ to atm]
     real(r8), pointer :: eflx_soil_grnd_patch    (:)   ! patch soil heat flux (W/m**2) [+ = into soil] 
     real(r8), pointer :: eflx_lwrad_net_patch    (:)   ! patch net infrared (longwave) rad (W/m**2) [+ = to atm]
     real(r8), pointer :: eflx_lwrad_out_patch    (:)   ! patch emitted infrared (longwave) radiation (W/m**2)
     real(r8), pointer :: eflx_snomelt_col        (:)   ! col snow melt heat flux (W/m**2)
     real(r8), pointer :: eflx_gnet_patch         (:)   ! patch net heat flux into ground  (W/m**2)
     real(r8), pointer :: eflx_grnd_lake_patch    (:)   ! patch net heat flux into lake / snow surface, excluding light transmission (W/m**2)
     real(r8), pointer :: eflx_dynbal_grc         (:)   ! grc dynamic land cover change conversion energy flux (W/m**2)
     real(r8), pointer :: eflx_bot_col            (:)   ! col heat flux from beneath the soil or ice column (W/m**2)
     real(r8), pointer :: eflx_fgr12_col          (:)   ! col ground heat flux between soil layers 1 and 2 (W/m**2)
     real(r8), pointer :: eflx_fgr_col            (:,:) ! col (rural) soil downward heat flux (W/m2) (1:nlevgrnd)  (pos upward; usually eflx_bot >= 0)

     ! Derivatives of energy fluxes
     real(r8), pointer :: dgnetdT_patch           (:)   ! patch derivative of net ground heat flux wrt soil temp  (W/m**2 K)
     real(r8), pointer :: netrad_patch            (:)   ! col net radiation (W/m**2) [+ = to sfc]
     real(r8), pointer :: cgrnd_patch             (:)   ! col deriv. of soil energy flux wrt to soil temp [W/m2/k]
     real(r8), pointer :: cgrndl_patch            (:)   ! col deriv. of soil latent heat flux wrt soil temp  [W/m**2/k]
     real(r8), pointer :: cgrnds_patch            (:)   ! col deriv. of soil sensible heat flux wrt soil temp [W/m2/k]

     ! Canopy radiation
     real(r8), pointer :: dlrad_patch             (:)   ! col downward longwave radiation below the canopy [W/m2]
     real(r8), pointer :: ulrad_patch             (:)   ! col upward longwave radiation above the canopy [W/m2]

     ! Wind Stress
     real(r8), pointer :: taux_patch              (:)   ! patch wind (shear) stress: e-w (kg/m/s**2)
     real(r8), pointer :: tauy_patch              (:)   ! patch wind (shear) stress: n-s (kg/m/s**2)

     ! Conductance
     real(r8), pointer :: canopy_cond_patch       (:)   ! patch tracer conductance for canopy [m/s] 

     ! Transpiration
     real(r8), pointer :: btran_patch             (:)   ! patch transpiration wetness factor (0 to 1)
     real(r8), pointer :: btran_min_patch         (:)   ! patch daily minimum transpiration wetness factor (0 to 1)
     real(r8), pointer :: btran_min_inst_patch    (:)   ! patch instantaneous daily minimum transpiration wetness factor (0 to 1)
     real(r8), pointer :: bsun_patch              (:)   ! patch sunlit canopy transpiration wetness factor (0 to 1)
     real(r8), pointer :: bsha_patch              (:)   ! patch shaded canopy transpiration wetness factor (0 to 1)

     ! Roots
     real(r8), pointer :: btran2_patch            (:)   ! patch root zone soil wetness factor (0 to 1) 
     real(r8), pointer :: rresis_patch            (:,:) ! patch root resistance by layer (0-1)  (nlevgrnd)

     ! Latent heat
     real(r8), pointer :: htvp_col                (:)   ! latent heat of vapor of water (or sublimation) [j/kg]

     ! Balance Checks
     real(r8), pointer :: errsoi_patch            (:)   ! soil/lake energy conservation error   (W/m**2)
     real(r8), pointer :: errsoi_col              (:)   ! soil/lake energy conservation error   (W/m**2)
     real(r8), pointer :: errseb_patch            (:)   ! surface energy conservation error     (W/m**2)
     real(r8), pointer :: errseb_col              (:)   ! surface energy conservation error     (W/m**2)
     real(r8), pointer :: errsol_patch            (:)   ! solar radiation conservation error    (W/m**2)
     real(r8), pointer :: errsol_col              (:)   ! solar radiation conservation error    (W/m**2)
     real(r8), pointer :: errlon_patch            (:)   ! longwave radiation conservation error (W/m**2)
     real(r8), pointer :: errlon_col              (:)   ! longwave radiation conservation error (W/m**2)

   contains

     procedure, public  :: Init            ! Public initialization method
     procedure, private :: InitAllocate    ! initialize/allocate
     procedure, private :: InitHistory     ! setup history fields
     procedure, private :: InitCold        ! initialize for cold start
     procedure, public  :: Restart         ! setup restart fields

  end type energyflux_type

  character(len=*), parameter, private :: sourcefile = &
       __FILE__
  !------------------------------------------------------------------------

contains

  !------------------------------------------------------------------------
  subroutine Init(this, bounds, t_grnd_col)
    !
    ! !DESCRIPTION:
    !    Allocate and initialize the data type and setup history, and initialize for cold-start.
    ! !USES:
    implicit none
    ! !ARGUMENTS:
    class(energyflux_type)         :: this
    type(bounds_type) , intent(in) :: bounds  
    real(r8)          , intent(in) :: t_grnd_col( bounds%begc: )

    SHR_ASSERT_ALL((ubound(t_grnd_col) == (/bounds%endc/)), errMsg(sourcefile, __LINE__))

    call this%InitAllocate ( bounds )
    call this%InitHistory ( bounds)
    call this%InitCold ( bounds, t_grnd_col) 

  end subroutine Init

  !------------------------------------------------------------------------
  subroutine InitAllocate(this, bounds)
    !
    ! !DESCRIPTION:
    ! Initialize and allocate data structure
    !
    ! !USES:
    use shr_infnan_mod , only : nan => shr_infnan_nan, assignment(=)
    use clm_varpar     , only : nlevsno, nlevgrnd, nlevlak
    implicit none
    !
    ! !ARGUMENTS:
    class(energyflux_type) :: this
    type(bounds_type), intent(in) :: bounds  
    !
    ! !LOCAL VARIABLES:
    integer :: begp, endp
    integer :: begc, endc
    integer :: begl, endl
    integer :: begg, endg
    !------------------------------------------------------------------------

    begp = bounds%begp; endp= bounds%endp
    begc = bounds%begc; endc= bounds%endc
    begl = bounds%begl; endl= bounds%endl
    begg = bounds%begg; endg= bounds%endg

    allocate( this%eflx_sh_snow_patch      (begp:endp))             ; this%eflx_sh_snow_patch      (:)   = nan
    allocate( this%eflx_sh_soil_patch      (begp:endp))             ; this%eflx_sh_soil_patch      (:)   = nan
    allocate( this%eflx_sh_h2osfc_patch    (begp:endp))             ; this%eflx_sh_h2osfc_patch    (:)   = nan
    allocate( this%eflx_sh_tot_patch       (begp:endp))             ; this%eflx_sh_tot_patch       (:)   = nan
    allocate( this%eflx_sh_grnd_patch      (begp:endp))             ; this%eflx_sh_grnd_patch      (:)   = nan
    allocate( this%eflx_sh_veg_patch       (begp:endp))             ; this%eflx_sh_veg_patch       (:)   = nan
    allocate( this%eflx_sh_precip_conversion_col(begc:endc))        ; this%eflx_sh_precip_conversion_col(:) = nan
    allocate( this%eflx_lh_tot_patch       (begp:endp))             ; this%eflx_lh_tot_patch       (:)   = nan
    allocate( this%eflx_lh_grnd_patch      (begp:endp))             ; this%eflx_lh_grnd_patch      (:)   = nan
    allocate( this%eflx_lh_vege_patch      (begp:endp))             ; this%eflx_lh_vege_patch      (:)   = nan
    allocate( this%eflx_lh_vegt_patch      (begp:endp))             ; this%eflx_lh_vegt_patch      (:)   = nan
    allocate( this%eflx_soil_grnd_patch    (begp:endp))             ; this%eflx_soil_grnd_patch    (:)   = nan
    allocate( this%eflx_lwrad_net_patch    (begp:endp))             ; this%eflx_lwrad_net_patch    (:)   = nan
    allocate( this%eflx_lwrad_out_patch    (begp:endp))             ; this%eflx_lwrad_out_patch    (:)   = nan
    allocate( this%eflx_gnet_patch         (begp:endp))             ; this%eflx_gnet_patch         (:)   = nan
    allocate( this%eflx_grnd_lake_patch    (begp:endp))             ; this%eflx_grnd_lake_patch    (:)   = nan
    allocate( this%eflx_dynbal_grc         (begg:endg))             ; this%eflx_dynbal_grc         (:)   = nan
    allocate( this%eflx_bot_col            (begc:endc))             ; this%eflx_bot_col            (:)   = nan
    allocate( this%eflx_snomelt_col        (begc:endc))             ; this%eflx_snomelt_col        (:)   = nan
    allocate( this%eflx_fgr12_col          (begc:endc))             ; this%eflx_fgr12_col          (:)   = nan
    allocate( this%eflx_fgr_col            (begc:endc, 1:nlevgrnd)) ; this%eflx_fgr_col            (:,:) = nan

    allocate( this%dgnetdT_patch           (begp:endp))             ; this%dgnetdT_patch           (:)   = nan
    allocate( this%cgrnd_patch             (begp:endp))             ; this%cgrnd_patch             (:)   = nan
    allocate( this%cgrndl_patch            (begp:endp))             ; this%cgrndl_patch            (:)   = nan
    allocate( this%cgrnds_patch            (begp:endp))             ; this%cgrnds_patch            (:)   = nan
    allocate( this%dlrad_patch             (begp:endp))             ; this%dlrad_patch             (:)   = nan
    allocate( this%ulrad_patch             (begp:endp))             ; this%ulrad_patch             (:)   = nan
    allocate( this%netrad_patch            (begp:endp))             ; this%netrad_patch            (:)   = nan  

    allocate( this%taux_patch              (begp:endp))             ; this%taux_patch              (:)   = nan
    allocate( this%tauy_patch              (begp:endp))             ; this%tauy_patch              (:)   = nan

    allocate( this%canopy_cond_patch       (begp:endp))             ; this%canopy_cond_patch       (:)   = nan

    allocate( this%htvp_col                (begc:endc))             ; this%htvp_col                (:)   = nan

    allocate(this%rresis_patch             (begp:endp,1:nlevgrnd))  ; this%rresis_patch            (:,:) = nan
    allocate(this%btran_patch              (begp:endp))             ; this%btran_patch             (:)   = nan
    allocate(this%btran_min_patch          (begp:endp))             ; this%btran_min_patch         (:)   = nan
    allocate(this%btran_min_inst_patch     (begp:endp))             ; this%btran_min_inst_patch    (:)   = nan
    allocate(this%btran2_patch             (begp:endp))             ; this%btran2_patch            (:)   = nan
    allocate( this%bsun_patch              (begp:endp))             ; this%bsun_patch              (:)   = nan
    allocate( this%bsha_patch              (begp:endp))             ; this%bsha_patch              (:)   = nan
    allocate( this%errsoi_patch            (begp:endp))             ; this%errsoi_patch            (:)   = nan
    allocate( this%errsoi_col              (begc:endc))             ; this%errsoi_col              (:)   = nan
    allocate( this%errseb_patch            (begp:endp))             ; this%errseb_patch            (:)   = nan
    allocate( this%errseb_col              (begc:endc))             ; this%errseb_col              (:)   = nan
    allocate( this%errsol_patch            (begp:endp))             ; this%errsol_patch            (:)   = nan
    allocate( this%errsol_col              (begc:endc))             ; this%errsol_col              (:)   = nan
    allocate( this%errlon_patch            (begp:endp))             ; this%errlon_patch            (:)   = nan
    allocate( this%errlon_col              (begc:endc))             ; this%errlon_col              (:)   = nan

  end subroutine InitAllocate
    
  !------------------------------------------------------------------------
  subroutine InitHistory(this, bounds)
    !
    ! !DESCRIPTION:
    ! Setup fields that can be output to history files
    !
    ! !USES:
    use shr_infnan_mod , only : nan => shr_infnan_nan, assignment(=)
    use clm_varpar     , only : nlevsno, nlevgrnd
    use histFileMod    , only : hist_addfld1d, hist_addfld2d, no_snow_normal
    use ncdio_pio      , only : ncd_inqvdlen
    implicit none
    !
    ! !ARGUMENTS:
    class(energyflux_type) :: this
    type(bounds_type), intent(in) :: bounds  
    !
    ! !LOCAL VARIABLES:
    integer           :: begp, endp
    integer           :: begc, endc
    integer           :: begl, endl
    integer           :: begg, endg
    integer           :: dimlen
    integer           :: err_code
    logical           :: do_io
    character(10)     :: active
    real(r8), pointer :: data2dptr(:,:), data1dptr(:) ! temp. pointers for slicing larger arrays
    !------------------------------------------------------------------------

    begp = bounds%begp; endp= bounds%endp
    begc = bounds%begc; endc= bounds%endc
    begl = bounds%begl; endl= bounds%endl
    begg = bounds%begg; endg= bounds%endg


    this%eflx_dynbal_grc(begg:endg) = spval 
    call hist_addfld1d (fname='EFLX_DYNBAL',  units='W/m^2',  &
         avgflag='A', long_name='dynamic land cover change conversion energy flux', &
         ptr_lnd=this%eflx_dynbal_grc, default='inactive')

    this%eflx_snomelt_col(begc:endc) = spval
    call hist_addfld1d (fname='FSM',  units='W/m^2',  &
         avgflag='A', long_name='snow melt heat flux', &
         ptr_col=this%eflx_snomelt_col, c2l_scale_type='urbanf', default='inactive')

    call hist_addfld1d (fname='FSM_ICE', units='W/m^2',  &
         avgflag='A', long_name='snow melt heat flux (ice landunits only)', &
         ptr_col=this%eflx_snomelt_col, c2l_scale_type='urbanf', l2g_scale_type='ice', &
         default='inactive')

    this%eflx_lwrad_net_patch(begp:endp) = spval
    call hist_addfld1d (fname='FIRA', units='W/m^2',  &
         avgflag='A', long_name='net infrared (longwave) radiation', &
         ptr_patch=this%eflx_lwrad_net_patch, c2l_scale_type='urbanf', default='inactive')

    call hist_addfld1d (fname='FIRA_ICE', units='W/m^2',  &
         avgflag='A', long_name='net infrared (longwave) radiation (ice landunits only)', &
         ptr_patch=this%eflx_lwrad_net_patch, c2l_scale_type='urbanf', l2g_scale_type='ice',&
         default='inactive')

    this%eflx_lwrad_out_patch(begp:endp) = spval 
    call hist_addfld1d (fname='FIRE', units='W/m^2',  &
         avgflag='A', long_name='emitted infrared (longwave) radiation', &
         ptr_patch=this%eflx_lwrad_out_patch, c2l_scale_type='urbanf', default='inactive')

    call hist_addfld1d (fname='FIRE_ICE', units='W/m^2',  &
         avgflag='A', long_name='emitted infrared (longwave) radiation (ice landunits only)', &
         ptr_patch=this%eflx_lwrad_out_patch, c2l_scale_type='urbanf', l2g_scale_type='ice', &
         default='inactive')

    this%eflx_lh_grnd_patch(begp:endp) = spval
    call hist_addfld1d (fname='FGEV', units='W/m^2',  &
         avgflag='A', long_name='ground evaporation', &
         ptr_patch=this%eflx_lh_grnd_patch, c2l_scale_type='urbanf', default='inactive') 

    this%eflx_sh_tot_patch(begp:endp) = spval
    call hist_addfld1d (fname='FSH', units='W/m^2',  &
         avgflag='A', long_name='sensible heat not including correction for land use change and rain/snow conversion', &
         ptr_patch=this%eflx_sh_tot_patch, c2l_scale_type='urbanf', default='inactive')

    call hist_addfld1d (fname='FSH_ICE', units='W/m^2',  &
         avgflag='A', &
         long_name='sensible heat not including correction for land use change and rain/snow conversion (ice landunits only)', &
         ptr_patch=this%eflx_sh_tot_patch, c2l_scale_type='urbanf', l2g_scale_type='ice', &
         default='inactive')

    this%eflx_sh_tot_patch(begp:endp) = spval
    call hist_addfld1d (fname='Qh', units='W/m^2',  &
         avgflag='A', long_name='sensible heat', &
         ptr_patch=this%eflx_sh_tot_patch, c2l_scale_type='urbanf', &
         default = 'inactive')

    this%eflx_lh_tot_patch(begp:endp) = spval
    call hist_addfld1d (fname='Qle', units='W/m^2',  &
         avgflag='A', long_name='total evaporation', &
         ptr_patch=this%eflx_lh_tot_patch, c2l_scale_type='urbanf', &
         default = 'inactive')

    this%eflx_lh_tot_patch(begp:endp) = spval
    call hist_addfld1d (fname='EFLX_LH_TOT', units='W/m^2', &
         avgflag='A', long_name='total latent heat flux [+ to atm]', &
         ptr_patch=this%eflx_lh_tot_patch, c2l_scale_type='urbanf', default='inactive')

    call hist_addfld1d (fname='EFLX_LH_TOT_ICE', units='W/m^2',  &
         avgflag='A', long_name='total latent heat flux [+ to atm] (ice landunits only)', &
         ptr_patch=this%eflx_lh_tot_patch, c2l_scale_type='urbanf', l2g_scale_type='ice', &
         default='inactive')

    this%eflx_soil_grnd_patch(begp:endp) = spval
    call hist_addfld1d (fname='Qstor', units='W/m^2',  &
         avgflag='A', long_name='storage heat flux (includes snowmelt)', &
         ptr_patch=this%eflx_soil_grnd_patch, c2l_scale_type='urbanf', &
         default = 'inactive')

    this%eflx_sh_grnd_patch(begp:endp) = spval
    call hist_addfld1d (fname='FSH_G', units='W/m^2',  &
         avgflag='A', long_name='sensible heat from ground', &
         ptr_patch=this%eflx_sh_grnd_patch, c2l_scale_type='urbanf', default='inactive')

    this%eflx_soil_grnd_patch(begp:endp) = spval
    call hist_addfld1d (fname='FGR', units='W/m^2',  &
         avgflag='A', long_name='heat flux into soil/snow including snow melt and lake / snow light transmission', &
         ptr_patch=this%eflx_soil_grnd_patch, c2l_scale_type='urbanf', default='inactive')

    call hist_addfld1d (fname='FGR_ICE', units='W/m^2',  &
         avgflag='A', &
         long_name='heat flux into soil/snow including snow melt and lake / snow light transmission (ice landunits only)', &
         ptr_patch=this%eflx_soil_grnd_patch, c2l_scale_type='urbanf', l2g_scale_type='ice', &
         default='inactive')

    this%eflx_soil_grnd_patch(begp:endp) = spval
    call hist_addfld1d (fname='EFLX_SOIL_GRND', units='W/m^2', &
         avgflag='A', long_name='soil heat flux [+ into soil]', &
         ptr_patch=this%eflx_soil_grnd_patch, default='inactive', c2l_scale_type='urbanf')

    this%eflx_sh_precip_conversion_col(begc:endc) = spval
    call hist_addfld1d (fname = 'FSH_PRECIP_CONVERSION', units='W/m^2', &
         avgflag='A', long_name='Sensible heat flux from conversion of rain/snow atm forcing', &
         ptr_col=this%eflx_sh_precip_conversion_col, c2l_scale_type='urbanf', default='inactive')

    this%netrad_patch(begp:endp) = spval
    call hist_addfld1d (fname='Rnet', units='W/m^2',  &
         avgflag='A', long_name='net radiation', &
         ptr_patch=this%netrad_patch, c2l_scale_type='urbanf', &
         default='inactive')

    this%dgnetdT_patch(begp:endp) = spval
    call hist_addfld1d (fname='DGNETDT', units='W/m^2/K', &
         avgflag='A', long_name='derivative of net ground heat flux wrt soil temp', &
         ptr_patch=this%dgnetdT_patch, default='inactive', c2l_scale_type='urbanf')

    this%eflx_fgr12_col(begc:endc) = spval
    call hist_addfld1d (fname='FGR12',  units='W/m^2',  &
         avgflag='A', long_name='heat flux between soil layers 1 and 2', &
         ptr_col=this%eflx_fgr12_col, set_lake=spval, default='inactive')

    this%taux_patch(begp:endp) = spval
    call hist_addfld1d (fname='TAUX', units='kg/m/s^2',  &
         avgflag='A', long_name='zonal surface stress', &
         ptr_patch=this%taux_patch, default='inactive')

    this%tauy_patch(begp:endp) = spval
    call hist_addfld1d (fname='TAUY', units='kg/m/s^2',  &
         avgflag='A', long_name='meridional surface stress', &
         ptr_patch=this%tauy_patch, default='inactive')

    this%errsoi_col(begc:endc) = spval
    call hist_addfld1d (fname='ERRSOI',  units='W/m^2',  &
         avgflag='A', long_name='soil/lake energy conservation error', &
         ptr_col=this%errsoi_col, default='inactive')

    this%errseb_patch(begp:endp) = spval
    call hist_addfld1d (fname='ERRSEB',  units='W/m^2',  &
         avgflag='A', long_name='surface energy conservation error', &
         ptr_patch=this%errseb_patch, default='inactive')

    this%errsol_patch(begp:endp) = spval
    call hist_addfld1d (fname='ERRSOL',  units='W/m^2',  &
         avgflag='A', long_name='solar radiation conservation error', &
         ptr_patch=this%errsol_patch, set_urb=spval, default='inactive')

  end subroutine InitHistory

  !-----------------------------------------------------------------------
  subroutine InitCold(this, bounds, t_grnd_col)
    !
    ! !DESCRIPTION:
    ! Initialize cold start conditions for module variables
    !
    ! !USES:
    use shr_kind_mod    , only : r8 => shr_kind_r8
    use shr_const_mod   , only : SHR_CONST_TKFRZ
    use clm_varpar      , only : nlevsoi, nlevgrnd, nlevsno, nlevlak
    use clm_varcon      , only : denice, denh2o, sb
    use landunit_varcon , only : istwet, istsoil, istdlak
    use clm_varctl      , only : iulog
    implicit none
    !
    ! !ARGUMENTS:
    class(energyflux_type)         :: this
    type(bounds_type) , intent(in) :: bounds  
    real(r8)          , intent(in) :: t_grnd_col( bounds%begc: )
    !
    ! !LOCAL VARIABLES:
    integer  :: j,l,c,p,levs,lev
    !-----------------------------------------------------------------------

    SHR_ASSERT_ALL((ubound(t_grnd_col) == (/bounds%endc/)), errMsg(sourcefile, __LINE__))

    ! Patches
    do p = bounds%begp, bounds%endp 
       c = patch%column(p)
       l = patch%landunit(p)

       this%eflx_lwrad_out_patch(p) = sb * (t_grnd_col(c))**4
    end do

    ! initialize rresis, for use in ecosystemdyn
    do p = bounds%begp,bounds%endp
       do lev = 1,nlevgrnd
          this%rresis_patch(p,lev) = 0._r8
       end do
    end do 

  end subroutine InitCold

  !------------------------------------------------------------------------
  subroutine Restart(this, bounds, ncid, flag)
    ! 
    ! !DESCRIPTION:
    ! Read/Write module information to/from restart file.
    !
    ! !USES:
    use shr_log_mod, only : errMsg => shr_log_errMsg
    use spmdMod    , only : masterproc
    use abortutils , only : endrun
    use ncdio_pio  , only : file_desc_t, ncd_defvar, ncd_io, ncd_double, ncd_int, &
                            ncd_inqvdlen
    use restUtilMod
    use decompMod      , only : get_proc_global
    implicit none
    !
    ! !ARGUMENTS:
    class(energyflux_type) :: this
    type(bounds_type), intent(in)    :: bounds 
    type(file_desc_t), intent(inout) :: ncid   
    character(len=*) , intent(in)    :: flag   
    !
    ! !LOCAL VARIABLES:
    integer :: j,c ! indices
    integer :: dimlen
    integer :: err_code
    integer :: numl_global
    logical :: readvar      ! determine if variable is on initial file
    logical :: do_io
    !-----------------------------------------------------------------------

    call get_proc_global(nl=numl_global)
    call restartvar(ncid=ncid, flag=flag, varname='EFLX_LWRAD_OUT', xtype=ncd_double,  & 
         dim1name='pft', &
         long_name='emitted infrared (longwave) radiation', units='watt/m^2', &
         interpinic_flag='interp', readvar=readvar, data=this%eflx_lwrad_out_patch)

    call restartvar(ncid=ncid, flag=flag, varname='btran2', xtype=ncd_double,  &
         dim1name='pft', &
         long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%btran2_patch) 

    call restartvar(ncid=ncid, flag=flag, varname='BTRAN_MIN', xtype=ncd_double,  &
         dim1name='pft', &
         long_name='daily minimum of transpiration wetness factor', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%btran_min_patch) 

    call restartvar(ncid=ncid, flag=flag, varname='BTRAN_MIN_INST', xtype=ncd_double,  &
         dim1name='pft', &
         long_name='instantaneous daily minimum of transpiration wetness factor', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%btran_min_inst_patch) 

    call restartvar(ncid=ncid, flag=flag, varname='eflx_grnd_lake', xtype=ncd_double,  &
         dim1name='pft', &
         long_name='net heat flux into lake/snow surface, excluding light transmission', units='W/m^2', &
         interpinic_flag='interp', readvar=readvar, data=this%eflx_grnd_lake_patch)

  end subroutine Restart

end module EnergyFluxType
