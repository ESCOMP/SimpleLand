module lnd2atmType

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Handle atm2lnd, lnd2atm mapping
  !
  ! !USES:
  use shr_kind_mod  , only : r8 => shr_kind_r8
  use shr_infnan_mod, only : nan => shr_infnan_nan, assignment(=)
  use shr_log_mod   , only : errMsg => shr_log_errMsg
  use abortutils    , only : endrun
  use decompMod     , only : bounds_type
  use clm_varpar    , only : numrad, ndst, nlevgrnd !ndst = number of dust bins. 	! MML: ndst = 4 from clm varpar
  use clm_varcon    , only : spval
  use clm_varctl    , only : iulog
  !
  ! !PUBLIC TYPES:
  implicit none
  private

  ! ----------------------------------------------------
  ! land -> atmosphere variables structure
  !----------------------------------------------------
  type, public :: lnd2atm_type
     ! lnd->atm
     real(r8), pointer :: t_rad_grc          (:)   => null() ! radiative temperature (Kelvin)
     ! MML check tech note for examples on how to calculate this; use MO theory
     real(r8), pointer :: t_ref2m_grc        (:)   => null() ! 2m surface air temperature (Kelvin)
     real(r8), pointer :: q_ref2m_grc        (:)   => null() ! 2m surface specific humidity (kg/kg)
     real(r8), pointer :: u_ref10m_grc       (:)   => null() ! 10m surface wind speed (m/sec)
     real(r8), pointer :: h2osno_grc         (:)   => null() ! snow water (mm H2O)
     ! MML: change this so when its allocated it is size (:,mml_nsoi) ... in which case the dust would have to be the same size... hmm...
     real(r8), pointer :: h2osoi_vol_grc     (:,:) => null() ! volumetric soil water (0~watsat, m3/m3, nlevgrnd) (for dust model)
     ! MML: albedo (:,:) -> albd is direct, albd(:,1) direct vis, albd(:,2) direct nir
     ! 					 -> albi is diffuse, albi(:,1) diffuse vis, albi(:,2) diffuse nir (I THINK) 
     real(r8), pointer :: albd_grc           (:,:) => null() ! (numrad) surface albedo (direct)
     real(r8), pointer :: albi_grc           (:,:) => null() ! (numrad) surface albedo (diffuse)
     real(r8), pointer :: taux_grc           (:)   => null() ! wind stress: e-w (kg/m/s**2)
     real(r8), pointer :: tauy_grc           (:)   => null() ! wind stress: n-s (kg/m/s**2)
     real(r8), pointer :: eflx_lh_tot_grc    (:)   => null() ! total latent HF (W/m**2)  [+ to atm]
     real(r8), pointer :: eflx_sh_tot_grc    (:)   => null() ! total sensible HF (W/m**2) [+ to atm]
     real(r8), pointer :: eflx_sh_ice_to_liq_col(:) => null() ! sensible HF generated from conversion of ice runoff to liquid (W/m**2) [+ to atm]
     real(r8), pointer :: eflx_lwrad_out_grc (:)   => null() ! IR (longwave) radiation (W/m**2)
     real(r8), pointer :: qflx_evap_tot_grc  (:)   => null() ! qflx_evap_soi + qflx_evap_can + qflx_tran_veg
     real(r8), pointer :: fsa_grc            (:)   => null() ! solar rad absorbed (total) (W/m**2)
     real(r8), pointer :: net_carbon_exchange_grc(:) => null() ! net CO2 flux (kg CO2/m**2/s) [+ to atm]
     real(r8), pointer :: nem_grc            (:)   => null() ! gridcell average net methane correction to CO2 flux (g C/m^2/s)
     real(r8), pointer :: ram1_grc           (:)   => null() ! aerodynamical resistance (s/m)
     real(r8), pointer :: fv_grc             (:)   => null() ! friction velocity (m/s) (for dust model)
     real(r8), pointer :: flxdst_grc         (:,:) => null() ! dust flux (size bins)
     real(r8), pointer :: ddvel_grc          (:,:) => null() ! dry deposition velocities
     ! lnd->rof
     real(r8), pointer :: qflx_rofliq_grc         (:)   => null() ! rof liq forcing
     real(r8), pointer :: qflx_rofliq_qsur_grc    (:)   => null() ! rof liq -- surface runoff component
     real(r8), pointer :: qflx_rofliq_qsub_grc    (:)   => null() ! rof liq -- subsurface runoff component
     real(r8), pointer :: qflx_rofliq_qgwl_grc    (:)   => null() ! rof liq -- glacier, wetland and lakes water balance residual component
     real(r8), pointer :: qflx_rofliq_h2osfc_grc  (:)   => null() ! rof liq -- surface water runoff component
     real(r8), pointer :: qflx_rofliq_drain_perched_grc    (:)   => null() ! rof liq -- perched water table runoff component
     real(r8), pointer :: qflx_rofice_grc    (:)   => null() ! rof ice forcing
     real(r8), pointer :: qflx_liq_from_ice_col(:) => null() ! liquid runoff from converted ice runoff

   contains

     procedure, public  :: Init
     procedure, private :: InitAllocate 
     procedure, private :: InitHistory  

  end type lnd2atm_type
  !------------------------------------------------------------------------

  character(len=*), parameter, private :: sourcefile = &
       __FILE__
  !------------------------------------------------------------------------

contains

  !------------------------------------------------------------------------
  subroutine Init(this, bounds)

    class(lnd2atm_type) :: this
    type(bounds_type), intent(in) :: bounds  

    call this%InitAllocate(bounds)
    call this%InitHistory(bounds)
    
  end subroutine Init

  !------------------------------------------------------------------------
  subroutine InitAllocate(this, bounds)
    !
    ! !DESCRIPTION:
    ! Initialize lnd2atm derived type
    !
    ! !USES
    use clm_varcon, only: sb, tfrz
    !
    ! !ARGUMENTS:
    class (lnd2atm_type) :: this
    type(bounds_type), intent(in) :: bounds  
    !
    ! !LOCAL VARIABLES:
    real(r8) :: ival  = 0.0_r8  ! initial value
    integer  :: begc, endc
    integer  :: begg, endg
    !------------------------------------------------------------------------

    begc = bounds%begc; endc = bounds%endc
    begg = bounds%begg; endg = bounds%endg

    allocate(this%t_rad_grc          (begg:endg))            ; this%t_rad_grc          (:)   = tfrz + 2._r8
    allocate(this%t_ref2m_grc        (begg:endg))            ; this%t_ref2m_grc        (:)   =ival
    allocate(this%q_ref2m_grc        (begg:endg))            ; this%q_ref2m_grc        (:)   =ival
    allocate(this%u_ref10m_grc       (begg:endg))            ; this%u_ref10m_grc       (:)   =ival
    allocate(this%h2osno_grc         (begg:endg))            ; this%h2osno_grc         (:)   = 0._r8
    allocate(this%h2osoi_vol_grc     (begg:endg,1:nlevgrnd)) ; this%h2osoi_vol_grc     (:,:) =ival
    allocate(this%albd_grc           (begg:endg,1:numrad))   ; this%albd_grc           (:,:) = 0.2_r8
    allocate(this%albi_grc           (begg:endg,1:numrad))   ; this%albi_grc           (:,:) = 0.2_r8
    allocate(this%taux_grc           (begg:endg))            ; this%taux_grc           (:)   =ival
    allocate(this%tauy_grc           (begg:endg))            ; this%tauy_grc           (:)   =ival
    allocate(this%eflx_lwrad_out_grc (begg:endg))            ; this%eflx_lwrad_out_grc (:)   = sb * tfrz**4
    allocate(this%eflx_sh_tot_grc    (begg:endg))            ; this%eflx_sh_tot_grc    (:)   =ival
    allocate(this%eflx_sh_ice_to_liq_col(begc:endc))         ; this%eflx_sh_ice_to_liq_col(:) = ival
    allocate(this%eflx_lh_tot_grc    (begg:endg))            ; this%eflx_lh_tot_grc    (:)   =ival
    allocate(this%qflx_evap_tot_grc  (begg:endg))            ; this%qflx_evap_tot_grc  (:)   =ival
    allocate(this%fsa_grc            (begg:endg))            ; this%fsa_grc            (:)   =ival
    allocate(this%net_carbon_exchange_grc(begg:endg))        ; this%net_carbon_exchange_grc(:) =ival
    allocate(this%nem_grc            (begg:endg))            ; this%nem_grc            (:)   =ival
    allocate(this%ram1_grc           (begg:endg))            ; this%ram1_grc           (:)   =ival
    allocate(this%fv_grc             (begg:endg))            ; this%fv_grc             (:)   =ival
    allocate(this%flxdst_grc         (begg:endg,1:ndst))     ; this%flxdst_grc         (:,:) =ival
    allocate(this%qflx_rofliq_grc    (begg:endg))            ; this%qflx_rofliq_grc    (:)   =ival
    allocate(this%qflx_rofliq_qsur_grc    (begg:endg))       ; this%qflx_rofliq_qsur_grc    (:)   =ival
    allocate(this%qflx_rofliq_qsub_grc    (begg:endg))       ; this%qflx_rofliq_qsub_grc    (:)   =ival
    allocate(this%qflx_rofliq_qgwl_grc    (begg:endg))       ; this%qflx_rofliq_qgwl_grc    (:)   =ival
    allocate(this%qflx_rofliq_h2osfc_grc  (begg:endg))       ; this%qflx_rofliq_h2osfc_grc    (:)   =ival
    allocate(this%qflx_rofliq_drain_perched_grc    (begg:endg))       ; this%qflx_rofliq_drain_perched_grc    (:)   =ival
    allocate(this%qflx_rofice_grc    (begg:endg))            ; this%qflx_rofice_grc    (:)   =ival
    allocate(this%qflx_liq_from_ice_col(begc:endc))          ; this%qflx_liq_from_ice_col(:) = ival

  end subroutine InitAllocate

  !-----------------------------------------------------------------------
  subroutine InitHistory(this, bounds)
    !
    ! !USES:
    use histFileMod, only : hist_addfld1d
    !
    ! !ARGUMENTS:
    class(lnd2atm_type) :: this
    type(bounds_type), intent(in) :: bounds  
    !
    ! !LOCAL VARIABLES:
    integer  :: begg, endg
    !---------------------------------------------------------------------

    begg = bounds%begg; endg = bounds%endg

    this%eflx_sh_tot_grc(begg:endg) = 0._r8
    call hist_addfld1d (fname='FSH_TO_COUPLER', units='W/m^2',  &
         avgflag='A', &
         long_name='sensible heat sent to coupler &
              &(includes corrections for land use change, rain/snow conversion and conversion of ice runoff to liquid)', &
         ptr_lnd=this%eflx_sh_tot_grc)

  end subroutine InitHistory

end module lnd2atmType
