module mml_mainMod
  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  !    Module which will contain a simple bucket land model to overwrite
  !	   CLM output. For now, just try to get it to compile with the CLM 
  !	   code. 
  !
  !  Interface with the CLM code in main/clm_driver.F90, before the call
  !  to lnd2atm (and before the history files are written), such that we 
  !  can overwrite the CLM data that gets passed to the coupler and 
  !  written into the CLM history files. 
  !
  !  Operates on gridscale level atmospheric data to do a simplified 
  !  land energy and hydrology budget, then return values to the atmosphere.
  !
  ! !Author: Marysa Lague 
  !	!Date created: 2016.01.13
  !
  !-----------------------------------------------------------------------

  
  ! !USES:
  ! MML: bounds & data type
  use decompMod , 	only : bounds_type
  use atm2lndType,	only : atm2lnd_type
  use lnd2atmType, 	only : lnd2atm_type  ! MML: probably going to need a lnd2atm type
  ! to hand to the coupler as data coming from the land going to the atmosphere (l2x) 
  ! Though actually, I suppose one could store them in atm2lnd_inst, then when the coupler
  ! code lnd_import_exportMod asks for lnd2atm data, just grab it from the atm2lnd_inst 
  ! structure

  !MML: do need these:
  use shr_log_mod     , only : errMsg => shr_log_errMsg
  use abortutils      , only : endrun
  use clm_varctl      , only : single_column, iulog
  use shr_kind_mod    , only : r8 => shr_kind_r8
  use shr_infnan_mod  , only : isnan => shr_infnan_isnan
  
  use QSatMod		  , only : QSat
  use perf_mod			! for t_startf and t_stopf
   
  ! For using month-dependent values from forcing files
  use clm_time_manager, only : get_curr_date, get_nstep, get_step_size
  
  ! For namelist var
  use clm_varctl       , only: mml_surdat
  
  implicit none
  
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: mml_main
  
  ! !PRIVATE MEMBER FUNCTIONS:
  private :: nc_import
  
  ! mml: can I store the subroutines in other files? (just to keep this one from getting outrageously long?
  ! try it... (move nc_import, for starters... )
  ! Well, that didn't work (moving nc_import to its own .F90 file), so I'll stick with piling all
  ! the subroutines at the end of this file, but ask someone who knows Fortran better if I can
  ! get around this... without making a whole bunch of module
  ! (or, do I need to say "use nc_import" up at the top? 
  private	:: mo_hybrid 	! follow Gordon's matlab example functino hybrid.m to est the obu root
  private	:: most			! follow Gordon's matlab function for monin-obukhov similarity theory; solve 
  								! this function for a given set of values
  private 	:: psi_m
  private 	:: psi_h
  private 	:: zbrent
  
  private	:: soil_thermal_properties 	! Again, following Gordon's matlab script; account for
  								! soil type and how much of the freezable water in each layer is frozen
  							!	 to get a tk(i) and cv(i) for each layer i at each time step
  private	:: phase_change	! after we get the soil temperature profile, do phase change in that layer
  								! and update soil temperatures accordingly
 ! private	:: satvap		! calculate saturation vapour pressure and deriv at given Ts
 			! oops, this was an old polynomial (just have looked at Gordon's old LSM code to figure out fortran polynomials,
 			! then used those coeffs instead of the more recent one!) Instead, I'm using the equivalent, but newer, clm 
 			! function QSat, and doing the lhflx calculations with specific humidity rather than saturation vapour pressure
  
 

contains
  
  !-----------------------------------------------------------------------
  subroutine mml_main (bounds, atm2lnd_inst, lnd2atm_inst) !lnd2atm_inst
  
    implicit none
  
    type(bounds_type), intent(in) :: bounds
    type(atm2lnd_type), intent(inout) :: atm2lnd_inst 
    type(lnd2atm_type), intent(inout) :: lnd2atm_inst 
     
    
	!-----------------------------------------------------------------------
    ! !LOCAL VARIABLES:
    
    integer begg, endg
    integer g
    
    real(r8)alpha, stefan		! temporary albedo placeholder
     
    character(len=256) :: lfsurdat	! surface prescribed data file name !MML make flexible length!
    ! len=64, mml changed to see if it gets upset about tailing white space
 
    character(len=32) :: subname = 'testMod_mml'
   
    real(r8) 	:: dummy
    
    ! Define a whole bunch of physical constants... here, or in main... maybe in main...
   ! some are site-specific though, and they need to be defined here (eg bucket capacity)
    real(r8) vkc, grav, tfrz, sigma, mmdry, mmh2o, cpd, cpw, rgas, cwat, cice, &
    	rhowat, rhoice, cvwat, cvice, tkwat, tkice, hfus, hvap, hsub
    
	integer mml_nsoi	! numer of soil layers... make sure its consistent with atm2lndType and .nc file. Right now forcing 5
	
	!real(r8), allocatable, dimension (:) 	:: lambda, gamma, beta, snow_melt, t_to_snow, wat2snow ! hvap or hsub depending if it is below freezing or not
    !real(r8) gamma	! psychometric constant [Pa/K]
    !real(r8) beta	! soil wetness factor
    !real(r8) snow_melt, t_to_snow, wat2snow	! snow melting temporary vars
    real(r8) mms2kgm	! conversion factor to get precip from [mm/s] to [kg/m2] 
    
    ! soil temperature tri-diagonal solving placeholder variables
   ! real(r8) m		! placeholder for soil heat capacity
   ! real(r8) aa, bb, cc, dd, num, den ! tridiagonal local elements 
   ! real(r8), allocatable, dimension (:) :: cp, dp	! tridiagonal vector elements  (cprime, dprime)
   
 !   real(r8), allocatable, dimension (:) 	:: m, aa, bb, cc, dd, num, den
   ! real(r8), allocatable, dimension (:,:)	:: cp, dp
   
    !real(r8) epc, edif, err	! placeholder variables for soil energy tracking
   !	real(r8), allocatable, dimension (:) 	:: epc, edif, err, snow0, water0, obu_root, taux, tauy, Va, zeta
   
    ! Try NOT allocating, but specifically saying the variables are of size begg:endg
    real(r8)	::	esrf(bounds%begg:bounds%endg), desrf(bounds%begg:bounds%endg), &
    				dqsrf(bounds%begg:bounds%endg) !, qsrf(bounds%begg:bounds:endg) 	! qsrf already points to an atm2lnd var	
    
   
   	integer :: year    ! year (0, ...) for nstep+1
    integer :: mon     ! month (1, ..., 12) for nstep+1
    integer :: day     ! day of month (1, ..., 31) for nstep+1
    integer :: sec     ! seconds into current date for nstep+1
    integer :: mcdate  ! Current model date (yyyymmdd)
    
    real(r8)	:: dt	   ! length of time step, in seconds
   
   	integer :: i, j, ind , dew, A2, B2, P0	! index integers for looping
   	
   	real(r8)	tol, obu0, obu1
   	
!    	!Force-set a maximum snow value
	real(r8)	:: snowcap
!    	
   	
   	! Formerly "allocate" "deallocate" variables:
   	
   	!1d
   	real(r8)	::	mmair(bounds%begg:bounds%endg)	, 	&	! molar mass of air
   					res(bounds%begg:bounds%endg)	, 	&	! surface resistance [s/m] from .nc
   					dshflx(bounds%begg:bounds%endg)	, 	&	! change in sensible heat flux
   					dlhflx(bounds%begg:bounds%endg)	, 	&	! change in latent heat flux
   					dlwrad(bounds%begg:bounds%endg)	, 	&	! change in long wave radiation
   					f0(bounds%begg:bounds%endg)	, 	&		! energy flux into the soil (check?)
   					df0(bounds%begg:bounds%endg)	, 	&	! change in energy flux into the soil (check thats what f0 actually is)
   					esat(bounds%begg:bounds%endg)	, 	&	! saturation vapour pressure
   					desat(bounds%begg:bounds%endg)	, 	&	! change in sat vp pressure
   					temp(bounds%begg:bounds%endg)	, 	&	! temporary placeholder variable
   					alb_nir_dir(bounds%begg:bounds%endg)	, 	&	! albedo, near IR, direct, from .nc file
   					alb_vis_dir(bounds%begg:bounds%endg)	, 	&	! albedo, visible, direct, from .nc file
   					alb_nir_dif(bounds%begg:bounds%endg)	, 	&	! albedo, near IR, diffuse, from .nc file
   					alb_vis_dif(bounds%begg:bounds%endg)	, 	&	! albedo, visible, diffuse, from .nc file
   					! this needs to come from the netcdf file, and I need to remember it over multiple time steps, so 
   					! make it a part of atm2lnd_type !!!
   					!emiss(bounds%begg:bounds%endg)	, 	&			! emissivity, from .nc file
   					!glc_mask(bounds%begg:bounds%endg)	, 	&		! mask on glaciated points (1 = glacier, 0 = not), from .nc file (aiming for greenalnd + antarctica)
   					!dust(bounds%begg:bounds%endg,4)	, 	&			! dust flux, land to atm, from .nc file  lat x lon x 3 dust bins
   					zref_t(bounds%begg:bounds%endg)	, 	&			! reference height temperature for lnd2atm
   					zref_u(bounds%begg:bounds%endg)	, 	&			! reference height wind speed for lnd2atm
   					zref_q(bounds%begg:bounds%endg)	, 	&			! reference height humidity for lnd2atm
   					lwrad(bounds%begg:bounds%endg)	, 	&			! incoming longwave radiation from atm
   					lw_abs(bounds%begg:bounds%endg)	, 	&			! absorbed longwave radiation (emissivity*incoming)
   					lambda(bounds%begg:bounds%endg)	, 	&			! latent heat of vaporization, or fusion, depending on phase
   					gamma(bounds%begg:bounds%endg)	, 	&			! 
   					beta(bounds%begg:bounds%endg)	, 	&			!
   					!snow_melt(bounds%begg:bounds%endg)	, 	&		! how much snow melted
   					t_to_snow(bounds%begg:bounds%endg)	, 	&		! temperature change to be applied to phase change instead
   					wat2snow(bounds%begg:bounds%endg)	, 	&		! how much water to steal from snow bucket (or snow from water...)
   					! soil temperature tri-diagonal solving placeholder variables
   					m(bounds%begg:bounds%endg,10)	, 	&	
   					aa(bounds%begg:bounds%endg,10)	, 	&
   					bb(bounds%begg:bounds%endg,10)	, 	&
   					cc(bounds%begg:bounds%endg,10)	, 	&
   					dd(bounds%begg:bounds%endg,10)	, 	&
   					ee(bounds%begg:bounds%endg,10)	, 	&
   					ff(bounds%begg:bounds%endg,10)	, 	&
   					num(bounds%begg:bounds%endg)	, 	&
   					den(bounds%begg:bounds%endg)	, 	&
   					! Energy conservation
   					epc(bounds%begg:bounds%endg)	, 	&			! 
   					edif(bounds%begg:bounds%endg)	, 	&			!
   					err(bounds%begg:bounds%endg)	, 	&			!
   					! snow melt rate tracking vars
   					ptl_snow_melt(bounds%begg:bounds%endg)	, 	&	
   					max_snow_melt(bounds%begg:bounds%endg)	, 	&	
   					! Other vars
   					theta_srf(bounds%begg:bounds%endg)	, 	&
   					snow0(bounds%begg:bounds%endg)	, 	&
   					water0(bounds%begg:bounds%endg)	, 	&
   					obu_old(bounds%begg:bounds%endg)	, 	&
   					obu_new(bounds%begg:bounds%endg)	, 	&
   					fx(bounds%begg:bounds%endg)	, 	&
   					obu_root(bounds%begg:bounds%endg)	, 	&
   					!taux(bounds%begg:bounds%endg)	, 	&
   					!tauy(bounds%begg:bounds%endg)	, 	&
   					Va(bounds%begg:bounds%endg)	, 	&
   					zeta(bounds%begg:bounds%endg)	
   	
   	! For the lhflx limitation				
   	real(r8)	::	pbot(bounds%begg:bounds%endg)	,	&		! [Pa] midpoint of bottom layer (from atm)
   					p2(bounds%begg:bounds%endg)	,	&		! [Pa] top boundary of bottom layer (calculate using hybrid coords)
   					qbot(bounds%begg:bounds%endg)	,	&		! [kg/kg] specific humidity in lowest level of atm (check units?)
   					dpbot(bounds%begg:bounds%endg)	,	&		! thickness in pressure of bottom layer, approximating as dpbot = 2*(psrf-pbot)
   					q_avail(bounds%begg:bounds%endg)	,	&	! water available in lowest level
   					lh_avail(bounds%begg:bounds%endg)	,	&	! latent heat available in lowest level
   					!lh_excess(bounds%begg:bounds%endg)	,	&	! excess latent heat 
   					!q_excess(bounds%begg:bounds%endg)	,	&	! excess water (qflux units) (kg/m2/s vs w/m2)
   					sh_init(bounds%begg:bounds%endg)	,	&
   					q_min(bounds%begg:bounds%endg)	,	&		! minimum water allowed to be in bottom layer of atm. Setting = 0 b/c I don't know what else to use... [kg/kg]
   					ztodt(bounds%begg:bounds%endg)	,	&		! 2 delta t
   					srfrpdel(bounds%begg:bounds%endg)			! 1./(pint(K+1)-pint(K)) = 1/(dpbot)
   	
   	logical		:: 	is_qneg										! if q goes negative anywhere (ie dew exceeds atm water), toggle this
   	
   	!real(r8)	::	grav	! 9.81 m/s2
   	!2d
   	! relies on knowing that mml_nsoi = 10
   	real(r8)	::	tsoi0(bounds%begg:bounds%endg,10)	, 	&
   					dtsoi(bounds%begg:bounds%endg,10)	, 	&
   					cp(bounds%begg:bounds%endg,10)	, 	&
   					dp(bounds%begg:bounds%endg,10)	, 	&
   					fsds_dir(bounds%begg:bounds%endg,2)	, 	&
   					fsds_dif(bounds%begg:bounds%endg,2)	, 	&
   					sw_abs_dir(bounds%begg:bounds%endg,2)	, 	&
   					sw_abs_dif(bounds%begg:bounds%endg,2)	
   real(r8) :: fsds_tot      ! Total solar
   
   !-----------------------------------------------------------------------
   ! MML: associate the simple land model variables with their counterparts in atm2lnd
     ! This gives them easier to work with names for use in this subroutine.
     ! right now I'm giving values that are coming straight from the atmosphere the tag
     ! *_atm, and variables that are calculated with the simple land model the tag _mml
     ! until I come up with an acronym for the simple model. 
   associate (  &
     ! atm vars; pass these to their land model equivalents
     !tbot_atm      		=> atm2lnd_inst%forc_t_not_downscaled_grc  		, 	&
     !rain_in_atm 		=> atm2lnd_inst%forc_rain_not_downscaled_grc	,	&
     !snow_in_atm		=> atm2lnd_inst%forc_snow_not_downscaled_grc	, 	&
     !lwdown_atm			=> atm2lnd_inst%forc_lwrad_not_downscaled_grc	, 	&
     !swdown_atm			=> atm2lnd_inst%forc_solar_grc					,	&
     ! atm vars (need to grab them from atm portion, though... once this is written, can simplify by grabbing them right away)
     fsds	   	=> atm2lnd_inst%mml_atm_fsds_grc		,	&
      fsdsnd	   	=> atm2lnd_inst%mml_atm_fsdsnd_grc		,	&	! incoming shortwave nir direct
      fsdsvd	   	=> atm2lnd_inst%mml_atm_fsdsvd_grc		,	&
      fsdsni	   	=> atm2lnd_inst%mml_atm_fsdsni_grc		,	&
      fsdsvi	   	=> atm2lnd_inst%mml_atm_fsdsvi_grc		,	&
     lwdn		=> atm2lnd_inst%mml_atm_lwdn_grc    	,	&
     zref		=> atm2lnd_inst%mml_atm_zref_grc    	,	&
     tref		=> atm2lnd_inst%mml_atm_tbot_grc    	,	&
     thref		=> atm2lnd_inst%mml_atm_thref_grc    	,	&
     qref		=> atm2lnd_inst%mml_atm_qbot_grc    	,	&
     uref		=> atm2lnd_inst%mml_atm_uref_grc    	,	&
     eref		=> atm2lnd_inst%mml_atm_eref_grc    	,	&
     pref		=> atm2lnd_inst%mml_atm_pbot_grc    	,	&
     psrf		=> atm2lnd_inst%mml_atm_psrf_grc    	,	&
     rhomol		=> atm2lnd_inst%mml_atm_rhomol_grc    	,	&
     rhoair		=> atm2lnd_inst%mml_atm_rhoair_grc    	,	&
     cpair		=> atm2lnd_inst%mml_atm_cp_grc    		,	& ! MML: this is in 
     pco2		=> atm2lnd_inst%mml_atm_pco2			,	&
     prec_liq	=> atm2lnd_inst%mml_atm_prec_liq_grc    ,   &	! MML: in mm/s
     prec_frz	=> atm2lnd_inst%mml_atm_prec_frz_grc    ,   &
     ! lnd variables
     tsrf 		=> atm2lnd_inst%mml_lnd_ts_grc			,   &
     qsrf 		=> atm2lnd_inst%mml_lnd_qs_grc			,   &
     radforc	=> atm2lnd_inst%mml_lnd_qa_grc			,   &
     sw_abs		=> atm2lnd_inst%mml_lnd_swabs_grc		,   &
     fsr		=> atm2lnd_inst%mml_lnd_fsr_grc		,   &
      fsrnd		=> atm2lnd_inst%mml_lnd_fsrnd_grc		,   &
      fsrni		=> atm2lnd_inst%mml_lnd_fsrni_grc		,   &
      fsrvd		=> atm2lnd_inst%mml_lnd_fsrvd_grc		,   &
      fsrvi		=> atm2lnd_inst%mml_lnd_fsrvi_grc		,   &
     lwup		=> atm2lnd_inst%mml_lnd_lwup_grc		,   &
     fsns		=> atm2lnd_inst%mml_lnd_fsns_grc		,   &
     flns		=> atm2lnd_inst%mml_lnd_flns_grc		,   &
     shflx		=> atm2lnd_inst%mml_lnd_shflx_grc		,   &
     lhflx		=> atm2lnd_inst%mml_lnd_lhflx_grc		,   &
     gsoi		=> atm2lnd_inst%mml_lnd_gsoi_grc		,   &
     gsnow		=> atm2lnd_inst%mml_lnd_gsnow_grc		,   &
     evap		=> atm2lnd_inst%mml_lnd_evap_grc		,   &
     ustar		=> atm2lnd_inst%mml_lnd_ustar_grc		,   &
     tstar		=> atm2lnd_inst%mml_lnd_tstar_grc		,   &
     qstar		=> atm2lnd_inst%mml_lnd_qstar_grc		,   &
     tvstar		=> atm2lnd_inst%mml_lnd_tvstar_grc		,   &
     obu		=> atm2lnd_inst%mml_lnd_obu_grc			,   &
     ram		=> atm2lnd_inst%mml_lnd_ram_grc			,   &
     rah		=> atm2lnd_inst%mml_lnd_rah_grc			,   &
     h_disp		=> atm2lnd_inst%mml_lnd_disp_grc		,   &
     z0m		=> atm2lnd_inst%mml_lnd_z0m_grc			,   &
     z0h		=> atm2lnd_inst%mml_lnd_z0h_grc			,   &
     albedo_fin	=> atm2lnd_inst%mml_lnd_alb_grc		,   & 
     snow_melt	=> atm2lnd_inst%mml_lnd_snowmelt		,	&
     taux		=> atm2lnd_inst%mml_out_taux			,	&
     tauy		=> atm2lnd_inst%mml_out_tauy			,	&
     ! over-large dew:
     lh_excess	=> atm2lnd_inst%mml_lh_excess		,   & 
     q_excess	=> atm2lnd_inst%mml_q_excess		,   & 
     lh_demand	=> atm2lnd_inst%mml_lh_demand		,   &
     q_demand	=> atm2lnd_inst%mml_q_demand		,   &
     ! soil variables
     tsoi		=> atm2lnd_inst%mml_soil_t_grc		,   &
     soil_liq	=> atm2lnd_inst%mml_soil_liq_grc		,   &
     soil_ice	=> atm2lnd_inst%mml_soil_ice_grc		,   &
     soil_dz	=> atm2lnd_inst%mml_soil_dz_grc		,   &
     soil_zh	=> atm2lnd_inst%mml_soil_zh_grc		,   &
     soil_tk	=> atm2lnd_inst%mml_soil_tk_grc		,   &
     soil_tk_1d	=> atm2lnd_inst%mml_soil_tk_1d_grc		,   &
     soil_tkh	=> atm2lnd_inst%mml_soil_tkh_grc		,   &
     soil_dtsoi	=> atm2lnd_inst%mml_soil_dtsoi_grc		,   &
     soil_cv	=> atm2lnd_inst%mml_soil_cv_grc		,   &
     soil_cv_1d	=> atm2lnd_inst%mml_soil_cv_1d_grc		,   &
	 glc_tk_1d	=> atm2lnd_inst%mml_glc_tk_1d_grc		,   &
	 glc_cv_1d	=> atm2lnd_inst%mml_glc_cv_1d_grc		,   &
     water		=> atm2lnd_inst%mml_soil_water_grc		,   &
     snow		=> atm2lnd_inst%mml_soil_snow_grc		,   &
     runoff		=> atm2lnd_inst%mml_soil_runoff_grc		,   &
     ! values from .nc file
     albedo_gvd	=> atm2lnd_inst%mml_nc_alb_gvd_grc			,	&
     albedo_svd	=> atm2lnd_inst%mml_nc_alb_svd_grc			,	&
     albedo_gnd	=> atm2lnd_inst%mml_nc_alb_gnd_grc			,	&
     albedo_snd	=> atm2lnd_inst%mml_nc_alb_snd_grc			,	&
     albedo_gvf	=> atm2lnd_inst%mml_nc_alb_gvf_grc			,	&
     albedo_svf	=> atm2lnd_inst%mml_nc_alb_svf_grc			,	&
     albedo_gnf	=> atm2lnd_inst%mml_nc_alb_gnf_grc			,	&
     albedo_snf	=> atm2lnd_inst%mml_nc_alb_snf_grc			,	&
     snowmask		=> atm2lnd_inst%mml_nc_snowmask_grc				,	&
     evaprs			=> atm2lnd_inst%mml_nc_evaprs_grc				,	&
     bucket_cap		=> atm2lnd_inst%mml_nc_bucket_cap_grc			,	&
     soil_maxice	=> atm2lnd_inst%mml_nc_soil_maxice_grc			,	&
     soil_z			=> atm2lnd_inst%mml_nc_soil_levels_grc			,	&
     soil_type		=> atm2lnd_inst%mml_nc_soil_type_grc			,	&
     roughness		=> atm2lnd_inst%mml_nc_roughness_grc			,	&
     emiss			=> atm2lnd_inst%mml_nc_emiss_grc			,	&
     glc_mask		=> atm2lnd_inst%mml_nc_glcmask_grc			,	&
     dust			=> atm2lnd_inst%mml_nc_dust_grc			,	&
     ! temporary diagnostics
     diag1_1d		=> atm2lnd_inst%mml_diag1_1d_grc			,	&
     diag2_1d		=> atm2lnd_inst%mml_diag2_1d_grc			,	&
     diag3_1d		=> atm2lnd_inst%mml_diag3_1d_grc			,	&
     diag1_2d		=> atm2lnd_inst%mml_diag1_2d_grc			,	&
     diag2_2d		=> atm2lnd_inst%mml_diag2_2d_grc			,	&
     diag3_2d		=> atm2lnd_inst%mml_diag3_2d_grc				&
     !ddvel_grc		=> lnd2atm_inst%ddvel_grc						&	! lat x lon x 3 dust bins
   )
     !-----------------------------------------------------------------------
     
     
     !-----------------------------------------------------------------------
     ! Assign local values

     begg = bounds%begg
     endg = bounds%endg
     mml_nsoi = 10
     
     ! Maximum allowed snow:
     snowcap = 5000.0_r8    ! somewhat arbitrary... thats 10m of snow at 500 kg/m3 (mid-value for a firn)
     ! Paterson, W.S.B. 1994. The Physics of Glaciers
     ! kg/m3
     ! 
     ! New snow (immediately after falling in calm)	50-70
     ! Damp new snow	100-200
     ! Settled snow	200-300
     ! Depth hoar	100-300
     ! Wind packed snow	350-400
     ! Firn	400-830
     ! Very wet snow and firn	700-800
     ! Glacier ice	830-917
     !
     
     
     
!    ! GBB: You probably do not have to allocate memory if these variabls are local
!	! to this routine. You should be able to use: 
!	! real(r8) :: mmair(bounds%begg:bounds%endg)
!	! or you can pass bounds%begg and bounds%endg into this routine through the
!	! arguement list as begg and endg. Then use:
!	! real(r8) :: mmair(begg:endg)
!	!
  	 
     !lfsurdat = '/glade/u/home/mlague/cesmruns/simple_land/inputs/mml_prescribed_vals_current.nc'
     !lfsurdat = '/glade/p/work/mlague/CESM/simple_land/inputs/mml_prescribed_vals_current.nc'
	 lfsurdat = mml_surdat
	 !lfsurdat = '/home/mlague/cesm_runs/inputs/surfdat/mml_ceres_monthly_mean_albedo.nc'
	 !write(iulog,*)subname, 'MML lfsurdat = mml_surdat = ', mml_surdat
  	 
  	 !-----------------------------------------------------------------------
	 ! Physical constants
	 
	 ! GBB: CESM physical constants are in: cime/share/csm_share/shr/shr_const_mod.F90
	 ! and CLM physical constants are in: components/clm/src/main/clm_varcon.F90
	 !
	 ! MML: during code clean-up consider taking values from shr_const_mod.F90 (clm_varcon.f90 points to that)
	 ! 			NOTE though that Rgas in those modules is per kmol, not per mol, ie instead of 8.3 its 8314... 
	 !			adjust accordingly! (also check the others)
	 vkc = 0.4_r8								! von Karman constant
	 grav = 9.80665_r8						! gravity [m/s2]
	 tfrz = 273.15_r8							! freezing temperature [K]
	 sigma = 5.67e-08			!where do we put _r8 here?			! Stefan-Boltzmann constant [ W/m2/K4]
	 mmdry = 28.97_r8/1000.0_r8					! Molecular mass of dry air (kg/mol)
	 mmh2o = 18.02_r8/1000.0_r8					! Molecular mass of water vapour (kg/mol)
	 cpd = 1005.0_r8              			! Specific heat of dry air at constant pressure (J/kg/K)
	 cpw = 1846.0_r8             				! Specific heat of water vapor at constant pressure (J/kg/K)
	 rgas = 8.31446_r8            				! Universal gas constant (J/K/mol)
	 cwat = 4188.0_r8             				! Specific heat of water (J/kg/K)
	 cice = 2117.27_r8            				! Specific heat ice (J/kg/K)
	 rhowat = 1000.0_r8                       	! Density of water (kg/m3)
	 rhoice = 917.0_r8                        	! Density of ice (kg/m3)
	 cvwat = cwat * rhowat 					! Heat capacity of water (J/m3/K)
	 cvice = cice * rhoice 					! Heat capacity of ice (J/m3/K)
	 tkwat = 0.57_r8                          	! Thermal conductivity of water (W/m/K)
	 tkice = 2.29_r8                         	! Thermal conductivity of ice (W/m/K)
 	 hfus = 0.3337e6_r8                 	    ! Heat of fusion for water at 0 C (J/kg)
	 hvap = 2.501e6_r8                       	! Latent heat of evaporation (J/kg)
	 hsub = hfus + hvap	    				! Latent heat of sublimation (J/kg)
	 !-----------------------------------------------------------------------
	 
	 !emiss = 1._r8 ! 0.98 	! constant surface emissivity for longwave radiation
	 ! Need to add diurnal control on albedo
	 ! Need to add nir and vis fields, I think, for passing stuff back to atm
	 
	 
     !k = 1 ! start time = 1 when reading in .nc file? (I'm actually not sure what the nt=k is for)
  	 ! MML: loop over k=1:12 and read in each month's value, then store in a larger matrix
  	 
  	 !-----------------------------------------------------------------------
!  	 temp = 0
  	 
  	 !-----------------------------------------------------------------------
  	 ! Re-assign atmospheric forcing data to simple land model equivalent
  	 ! (this is all the "forcing" data
  	 fsds 		= atm2lnd_inst%forc_solar_grc
  	 fsds_dir	= atm2lnd_inst%forc_solad_grc
     fsds_dif	= atm2lnd_inst%forc_solai_grc
  	 lwdn 		= atm2lnd_inst%forc_lwrad_not_downscaled_grc
  	 zref		= atm2lnd_inst%forc_hgt_grc   ! Note, there is a u, t , and q height in atm2lnd... compare? 
 ! GBB: No need to use the separate values for t, u, q; only need zref
 ! MML: Keith said there are 3 separate ones for historical reasons, but all three should be the same as zref    
     zref_t		= atm2lnd_inst%forc_hgt_t_grc
     zref_u		= atm2lnd_inst%forc_hgt_u_grc
     zref_q		= atm2lnd_inst%forc_hgt_q_grc
     tref		= atm2lnd_inst%forc_t_not_downscaled_grc ! is this right? or does atm have a ref height value?
     uref		= atm2lnd_inst%forc_wind_grc
     eref		= atm2lnd_inst%forc_vp_grc
     qref		= atm2lnd_inst%forc_q_not_downscaled_grc
     pref		= atm2lnd_inst%forc_pbot_not_downscaled_grc
     rhoair		= atm2lnd_inst%forc_rho_not_downscaled_grc 
     prec_liq	= atm2lnd_inst%forc_rain_not_downscaled_grc
     prec_frz	= atm2lnd_inst%forc_snow_not_downscaled_grc
     pco2 		= atm2lnd_inst%forc_pco2_grc
     ! For checking the big neg lhflx:
     psrf		= atm2lnd_inst%forc_psrf_grc  ! surface pressure (Pa)
     pbot		= atm2lnd_inst%forc_pbot_not_downscaled_grc ! not downscaled atm pressure (Pa)
     qbot		= atm2lnd_inst%forc_q_not_downscaled_grc   ! not downscaled atm specific humidity (kg/kg) 
     ! NOTE: this is NOT going to be consistent with CAM, still, if I use pbot and psrf as the "edges" of 
     ! my lowest atm layer; cam uses the actual pressure levels at the edges of the lowermost 
     ! atmospheric layer, but all I've got is pbot (which is likely in the middle of the lowest layer)
     ! ... can I extrapolate the upper and lower pressures from pbot, somehow? eg based off 
     ! temperature? Is mass conserved within a layer (ie knowing tbot, qbot, and pbot, could I get
     ! the thickness of the layer?)
     ! or use pbot - psrf as = 1/2 layer thickness, in pressure, which would give me the 
     ! total layer thickness, in pressure, which I should then be able to plug in to their equation.
     ! Yes, lets do it that way!  
     
     ! Put direct/diffuse fsds vis/nir into right variable to be output:
     fsdsnd = fsds_dir(:,2)
     fsdsvd = fsds_dir(:,1)
     fsdsni = fsds_dif(:,2)
     fsdsvi = fsds_dif(:,1)	! I think? check...

	 ! Theta = T + 0.0098 * z  (Gamma = 0.0098)
	 thref = tref + 0.0098_r8 * zref
	 
	 ! Have to calculate rhomol from the vapor pressure, actual pressure, and actual temperature
	 rhomol	= pref / (rgas*tref);
	  ! rho_mol = (pd + forcvar.eref)/(physcon.rgas * forcvar.tref)
	 ! rho_kg = ((pref - eref)*mmdry + eref*mmh2o)/(rgas*tref)
	 
	 ! MML: might need to move into g loop if I can't figure out how to allocate a matrix of size 
	 ! begg:endg before I know begg and endg ... 
	 ! calculate heat capacity based off specific humidity:
	 mmair = rhomol / rhoair 					! mol/kg
	 cpair = cpd * (1._r8 + (cpw/cpd - 1._r8)*qref)	! J/kg/K
	! cpair = mmair * cpair_kg					! J/mol/K
	! physcon.cpd * (1.0 + (physcon.cpw/physcon.cpd - 1.0) * forcvar.qref) (* mmair); 
	 
  	 ! compare atm reference height to t,u,q reference height. 
  	! write(iulog,*)subname, 'MML hgt = ', atm2lnd_inst%forc_hgt_grc , ' u_hgt = ', &
  	! 		atm2lnd_inst%forc_hgt_u_grc, ' t_hgt = ', atm2lnd_inst%forc_hgt_t_grc , &
  	! 		'q_hgt = ', atm2lnd_inst%forc_hgt_q_grc
  	!write(iulog,*)subname, 'MML in the land model...'
  	 
  	 !-----------------------------------------------------------------------
  	 
  	 
  	 	 
  	 !-----------------------------------------------------------------------
  	 ! Get outside data 
  	 
     !MML: Grab the current model time so we know what month we're in
     call get_curr_date(year, mon, day, sec)   ! Actually all I need for now is mon
     mcdate = year*10000 + mon*100 + day
  
  	!write(iulog,*)subname, 'MML month = ', mon
  	!write(iulog,*)subname, 'MML day = ', day
  	
  	 ! Get time step length (in seconds)
  	 dt = real( get_step_size(), r8 )
  	 
  	 ! conversion factor for precipitation (don't put it in the loop because it doesn't change)
  	 ! rain [mm/s] * [1 m / 1000 mm ] * rhowat [kg/m3] * dt [s] = rain [kg/m2]
  	 mms2kgm = rhowat * dt / 1000._r8
  	 ! note: I'm using the same conversion for snow and water, as snow is reported also in mm/s water equivalence units
  	 ! consider moving this up and multiplying prec_liq and prec_frz by this factor as we read them in
  	 
  	 
  	 !-----------------------------------------------------------------------
  	 ! read in .nc data

 ! 	 do i = 1, mml_nsoi
			
!		soil_z(:,i) =  -0.025 * (exp(0.5*(i-0.5)) - 1.0);
		
!		if (i == 1) then
!			soil_maxice(begg:endg,i) = 0._r8
!		else
!			soil_maxice(begg:endg,i) = 300._r8
!		end if
		
!	 enddo
	

  	! write(iulog,*) 'MML: Yikes! Pre-nc reading, albedo_gvd at some point begg = ', albedo_gvd(begg)
  	call t_startf('mml_nc_import')
  	
  		! ONLY actually run nc_import if we're on the first timestep of the first day of the month...
  		!if (sec <= 1800) then !( day == 1 .and. sec <= 1800) then
  		if ( day == 1 .and. sec .le. 1800 ) then
  			! <= 1800 will read it in both first 2 time steps... but after a restart it 
  			! seems to start on 1800, not 0, so it needs to be able to read them then, too...
  			! Is there a better way to say "if you haven't still got the last values, read these in?"
  			!
  			! Added the nc vars to the restart file, so maybe now I can revert to just saying if sec = 0? 
  			! (sec <1800) -> as long as that instance HAPPENS that would work... I think...
  			write(iulog,*)'reading netcdf data for mon=',mon,', day=',day,', sec=',sec,')'
  			
  			call nc_import(begg, endg, mml_nsoi, lfsurdat, mon, &
 					albedo_gvd(begg:endg), albedo_svd(begg:endg), &
 					albedo_gnd(begg:endg), albedo_snd(begg:endg), &
 					albedo_gvf(begg:endg), albedo_svf(begg:endg), &
 					albedo_gnf(begg:endg), albedo_snf(begg:endg), &
 					snowmask(begg:endg), evaprs(begg:endg), &
 					bucket_cap(begg:endg), & 
 					soil_type(begg:endg), roughness(begg:endg), &
 					emiss(begg:endg), glc_mask(begg:endg), dust(begg:endg,:), &
 					soil_tk_1d(begg:endg), soil_cv_1d(begg:endg), &
 					glc_tk_1d(begg:endg), glc_cv_1d(begg:endg)   ) !, &

                       write(iulog,*)'read netcdf'
 					
		end if
	call t_stopf('mml_nc_import')
	

		! Hard code snowmask and see if it'll run with the new files using that
		!snowmask(begg:endg) = 100.0_r8
			
     ! *************************************************************
     ! ***       Start the simple model (science part)  		 ***
     ! *************************************************************
     
     ! -------------------------------------------------------------
     ! -------- Some setup:
     !				- soil interface depths and dz soil thicknesses (based off soil_z)
     !				- albedo
     !				- update water bucket with precip here? (or at end? hmm.)
     !						Check with Gordon and GFDL
     !				- get radiative forcing (sw in and lw in, with albedo accounted for)
     ! -------------------------------------------------------------
 
     !write(iulog,*)'MML: Commence actually running the model!'    
 
     ! displacement height
     		! for now, set equal to 0.7 * canopy height
     h_disp = 0.7_r8 * roughness
     
     ! Roughness length for momentum
     		! for now, set equal to 0.1 * canopy height
     z0m	= 0.1_r8 * roughness
     
     ! Roughness length for heat
     		! for now, set equal to 0.1 * momentum roughness length
     z0h	= 0.1_r8 * z0m
     
     
     
     ! snow masking factor
     ! SHOULD ALWAYS BE BETWEEN 0 AND 1!!!!!
     
     ! If snow is negative (it shouldn't be, but if it went a bit neg), set temp = 0
     
     !temp(begg:endg) = snow(begg:endg)/(snow(begg:endg) + snowmask(begg:endg)) ! snow masking factor
     !diag3_1d = temp
    
     
 
     do g = begg, endg
          ! MML 2021.09.29: initialize temp as all zeros, otherwise it might just not have a value in some places!
          temp(g) = 0.0_r8	

          if (snowmask(g) < 0.0_r8) then
               ! this should not happen. Never feed in negative snowmask! But a person technically could do so, so catch it here:
               write(iulog,*)'warning: user provided snowmask(g)<0 (snowmask(g) = ',snowmask(g),'), setting snowmask(g)=100.0'
               snowmask(g) = 100.0_r8
          end if
 
  	  if ( snow(g) <= 0.0_r8 ) then
  	       temp(g) = 0.0_r8
  	       write(iulog,*)'warning: snow<0, setting snowmasking factor to zero. (snow(g) = ',snow(g),', overwriting so snow(g)=0.0)'
               snow(g) = 0.0_r8
  	  else
               if ( snow(g) + snowmask(g) == 0.0_r8 ) then
  	          temp(g) = 0.0_r8
               else
  	          temp(g) = snow(g) / ( snow(g) + snowmask(g) )
               end if
  	  end if
  	 
  	  diag3_1d(g) = temp(g)
  	 
  	  if ( temp(g) < 0 ) then
                 write(iulog,*)'Error: snow masking factor < 0 (should be between 0 and 1) \n ', &
                       	'Instead, snowmasking factor = ',temp(g)
                 call endrun(msg=errmsg(__FILE__, __LINE__))
          elseif ( temp(g) > 1 ) then
                 write(iulog,*)'Error: snow masking factor > 1 (should be between 0 and 1)  \n ', &
                  'Instead, snowmasking factor = ',temp(g)
                 call endrun(msg=errmsg(__FILE__, __LINE__))
        end if
  	 
  	 end do
     
     
     ! -------------------------------------------------------------
	 ! Albedo stuff
	 
     ! Direct/Diffuse Visible/NIR
     ! for consistent coding, shove vis and nir into a (:,2) sized matrix
     alb_vis_dir(begg:endg) = (1._r8 - temp(begg:endg)) * albedo_gvd(begg:endg) + &
     							temp(begg:endg) * albedo_svd(begg:endg)
     alb_nir_dir(begg:endg) = (1._r8 - temp(begg:endg)) * albedo_gnd(begg:endg) + &
	 							temp(begg:endg) * albedo_snd(begg:endg)
     alb_vis_dif(begg:endg) = (1._r8 - temp(begg:endg)) * albedo_gvf(begg:endg) + &
	 							temp(begg:endg) * albedo_svf(begg:endg)
     alb_nir_dif(begg:endg) = (1._r8 - temp(begg:endg)) * albedo_gnf(begg:endg) + &
	 							temp(begg:endg) * albedo_snf(begg:endg)
	 
	 ! for now, output one of these as albedo_fin just so there is a value:
	 !albedo_fin = alb_vis_dir
	 diag2_1d = alb_vis_dir
	 
	 !diag3_1d = alb_vis_dif	! why is the albedo going to 1e22 in h0? try this one...
	 
	 ! Do something special for albedo where there is a glacier? 
	 ! at present, I'm just feeding in albedos that already "make sense" for a glacier
	 
	
	 !albedo_fin = 0.3_r8 ! see if that overwrites... 
	 
	 ! -------------------------------------------------------------
	 ! Net radiation
     ! variables from atm: 			lwdn, fsds, fsds_dir, fsds_dif ! sw is handed as total, direct, and diffuse
     ! variables to end up with: 	sw_abs, fsr, radforc (into ground) 
     !
     ! for lw, emissivity = absorptivity
     ! alpha (albed0) = reflected, so (1-alpha) = absorbed
     
     ! longwave
     lw_abs(begg:endg) = emiss(begg:endg)*lwdn(begg:endg)
     lwup(begg:endg) = (1._r8 - emiss(begg:endg)) * lwdn(begg:endg) ! reflected longwave. Later, add surface emission 
     ! Shortwave direct visible
     sw_abs_dir(begg:endg,1) = (1._r8 - alb_vis_dir(begg:endg)) * fsds_dir(begg:endg,1)
     ! Shortwave direct NIR
     sw_abs_dir(begg:endg,2) = (1._r8 - alb_nir_dir(begg:endg)) * fsds_dir(begg:endg,2)
     ! Shortwave diffuse visible
     sw_abs_dif(begg:endg,1) = (1._r8 - alb_vis_dif(begg:endg)) * fsds_dif(begg:endg,1)
     ! Shortwave diffuse NIR
     sw_abs_dif(begg:endg,2) = (1._r8 - alb_nir_dif(begg:endg)) * fsds_dif(begg:endg,2)
     
     !fsr(begg:endg)  = alb_vis_dir(begg:endg) * fsds_dir(begg:endg,1) + &
    ! 					alb_nir_dir(begg:endg) * fsds_dir(begg:endg,2) + &
    ! 		 			alb_vis_dif(begg:endg) * fsds_dif(begg:endg,1) + &
    ! 		 			alb_nir_dif(begg:endg) * fsds_dif(begg:endg,2)
     
     ! fsr by  vis/nir/dir/dif
     fsrnd	=	alb_nir_dir(begg:endg) * fsds_dir(begg:endg,2)
     fsrni	=	alb_nir_dif(begg:endg) * fsds_dif(begg:endg,2)
     fsrvd	=	alb_vis_dir(begg:endg) * fsds_dir(begg:endg,1)
     fsrvi	=	alb_vis_dif(begg:endg) * fsds_dif(begg:endg,1)
     
     ! put sum of these in diag2, should equal fsr... well, it will. thats math. don't bother. 
     
     		 
     sw_abs(begg:endg) = sw_abs_dir(begg:endg,1) + sw_abs_dir(begg:endg,2) + &
     						sw_abs_dif(begg:endg,1) + sw_abs_dif(begg:endg,2)
     						
     						
     ! should be able to write like:
	 fsr(:)  = alb_vis_dir * fsds_dir(:,1) + &
     					alb_nir_dir * fsds_dir(:,2) + &
     		 			alb_vis_dif * fsds_dif(:,1) + &
     		 			alb_nir_dif * fsds_dif(:,2)
     		 
     sw_abs(:) = sw_abs_dir(:,1) + sw_abs_dir(:,2) + &
     						sw_abs_dif(:,1) + sw_abs_dif(:,2)
     
     
     ! Make output albedo to be a combination of all 4 albedo streams:
     albedo_fin(:) = 1.0e36_r8
     do g = begg, endg
        fsds_tot = fsds_dir(g,1) + fsds_dir(g,2) + fsds_dif(g,1) +  fsds_dif(g,2)
        if ( fsds_tot > 0.0_r8 )then
           albedo_fin(g) = fsr(g) / fsds_tot
        end if
     end do
     ! temporary fix:
     !lw_abs(begg:endg) = lwdn(begg:endg)
     !sw_abs(begg:endg) = 0.7*fsds(begg:endg)
     
     
     radforc(begg:endg) = lw_abs(begg:endg) + sw_abs(begg:endg)
     
     
     !-----------------------------------------------------------------------
  	 ! Initial Checks -> crash run if these fail
  	 
  	 do g = begg, endg
  	 
  	       if ( zref(g) < h_disp(g) ) then
            write(iulog,*)'Error: Forcing height is below canopy displacement height (zref < h_disp) '
            call endrun(msg=errmsg(__FILE__, __LINE__))
         end if
  	 
  	 end do
     
     
     ! -------------------------------------------------------------
     ! -------- Monin-Obukhov Stuff
     ! -------------------------------------------------------------
     
     ! call version of matlab model's hybrid.m and feval.m on the function most.m to get 
     ! an obu_root and a delta obu, as well as ustar, tsrat
     
     ! For testing purposes, make some dummy params ustar, tstar to get the obukhov bit,
     ! just so I know if the fluxes part is compiling.
     
     !ustar(begg:endg) = 0.002		! should I be doing begg:endg here, or can I just say (:)?
     !tstar(:) = 0.004				! Ahmed says (:) should work. In general, begg:endg is 
     								! good clm practice, but this is kind of a rogue piece 
     								! of code anyhow... 
     								!
     								! Check: is tstar theta_star ? I hope so. If not, convert.
     !qstar(:) = 1.0e-08				! kind of arbitraily assigned, lets see if it'll run
     
     
     obu0 = 100.0_r8
     obu1 = -100.0_r8
     tol = 0.01_r8
     
     ! loop over g to get a new obu length for each point g 
     do g = begg, endg
     
     	! solve for the obukhov length
    	call mo_hybrid ( zref(g), uref(g), thref(g), qref(g), h_disp(g), &	! in
    	 				z0m(g), z0h(g), tsrf(g), qsrf(g), vkc, grav, &	! in
  						ustar(g), tstar(g), tvstar(g), qstar(g), obu(g), & ! out
  						obu_root(g), obu0, obu1, tol)	! out
  		
  		! take that obukhov length (obu_root(g)) and evaluate the most function once more
  		! at that value of the obukhov length (to update ustar, tstar, etc?)			
  		
  		 !GBB: This is not essential. If the solver works fine, L, u*, and T* are all
		! consistent. Sometimes a root is not found, so this call evaluates u* and T*
		! with the final root estimate. 
  		call most ( obu_root(g), zref(g), uref(g), thref(g), qref(g), h_disp(g), &	! in
  					z0m(g), z0h(g), tsrf(g), qsrf(g), vkc, grav, &					! in
  					ustar(g), tstar(g), tvstar(g), qstar(g), obu(g), fx(g))		! out
     
     
     end do
     
     
     ! reset evaprs to something bigger for now
     ! evaprs(:) = 100
     
	! calculate aerodynamic resistances for momentum (ram) and heat (rah) in [s/m], and 
	! the effective resistance combining ram with the canopy resistance (res)	
     ram(:)		=	uref / (ustar * ustar)				! [s/m] = [m/s] / ([m/s] * [m/s])
     rah(:)		=	(thref - tsrf) / (ustar * tstar)	! [s/m] = [K] / ([m/s] * [K])
     res(:)		=  	(evaprs + rah)						! [s/m]
     
     ! cap res at 100,000 ()
     where ( res > 100000. )
		res(:) = 100000.0_r8
	 end where


     ! GBB: See what GFDL does for its evaporative resistance; should be a function
	 ! of stomatal conductance and LAI
	 
 	 ! Save initial temperature profile for energy conservation check:
 	 tsoi0(:,:) = tsoi
 	 
 	 ! Call soil thermal properties for this time step: (right now, it doesn't matter b/c 
 	 ! it doesn't have water dependence, or soil type dependence, for that matter, 
 	 ! but eventually I should add those
 	 call soil_thermal_properties (begg, endg, glc_mask, &
 	 							   soil_type, soil_z, soil_zh, soil_dz, &
 	 							   soil_liq, soil_ice, &
 	 							   soil_tk_1d, soil_cv_1d, &
  								   glc_tk_1d, glc_cv_1d, &
 	 							   soil_tk, soil_cv, soil_tkh, &
  								   mml_nsoi)
    
	! need to pass tsrf-tfrz to the satvap function, but I don't think I can do the calculation
	! in the actual subroutine call
	!temp = tsrf - tfrz;
	
    do g = begg, endg
     	
     	
   		! Call satvap. It would be nice to do this vectorized (whole matrix at once, not 
       	! looped, but my first stab at that didn't work (had trouble saying "if t>0" where t 
       	! is a matrix), so for now I'll put it in the loop just to get this working, then 
       	! try and code it better later.
     	
     	! MML: see qsat mod, and actually use that instead. This implementation is using the wrong polynomial!
     	! (out of date polynomial, 6 degree rather than 9). Also, change it to a specific humidity
     	! so that I'm consistently using a specific humidity everywhere 
     	
     	!call satvap( temp(g) , esat(g), desat(g)) 	! [K], [Pa], [Pa/K] 
     	
     	! I'm after surface values here, so use tsrf, and psrf (assume psrf ~ pbot ? that
     	! looks like that CanopyTemperatureMod is doing - check with Gordon)
     	call QSat (tsrf(g), pref(g), esrf(g), desrf(g), qsrf(g), dqsrf(g))
     	
     	! Okay, this is giving an updated qsrf - so should I call it before MO section? 
     	! otherwise MO uses the qsrf from the last time step, which is fine I guess, since this is
     	! using the tsrf from the last time step to get qsrf anyhow... think about this.
     	!
     	! ... also, I probably shouldn't have a variable called qsat if I'm also calling a function
     	! called QSat (is fortran case-sensitive?)... for now, go rename qsat to qsrf 
     	
     	! call QSat (T, p, es, esdT, qs, qsdT) ! in: T,p ; out: es, esdT, qs, qsdT
     	! T = temperature (K)
    	! p = surface atmospheric pressure (pa)
    	! In CanopyFluxesMod:
    	! call QSat(t_ref2m(p), forc_pbot(c), e_ref2m, de2mdT, qsat_ref2m, dqsat2mdT)
    	! 
    	! In CanopyTemperatureMod:
    	! call QSat(t_grnd(c), forc_pbot(c), eg, degdT, qsatg, qsatgdT)
    	! where eg = ! water vapor pressure at temperature T [pa]
     	
    	!evap(g) = esat(g) ! lets see if satvap is our culprit...
     	
    end do
     
   ! lwrad(:) = emiss * sigma * tsrf**4 
     
   ! lwup = lwup + lwrad ! (reflected longwave + sigma T^4)
     
    ! Now lets introduce some realistic physics...
     
    ! -----------------------------------------------------
    ! setup some parameters
     
    !Latent heat of Vapourization [J/kg] (sublimation if below freezing)
    ! GBB: hsub is used if snow is on the ground (check GFDL code). Or CLM uses hvap
    ! (gfdl says sublimation if snow, CLM says sublimation if frozen... check and make sure,
    ! then choose one and run with it)
    lambda(:) = hvap
	where ( tsrf < tfrz) lambda(:) = hsub	
	
	
	! Psychometric Constant [Pa/K]
	gamma(:) = cpair(:) * pref(:) / lambda(:)		! [J/kg/K] * [Pa] / [J/kg]
	
	!lhflx(:) = lambda
	
    ! --------------------------------------------------
	! ---- Surface Fluxes
	! -------------------------------------------------- 
     
    ! Emitted longwave radiation from surface [W/m2] and temperature derivative [W/m2/K]
	lwrad(:) 	=	emiss * sigma * tsrf**4
	dlwrad(:) 	=	4.0_r8 * emiss * sigma * tsrf**3
	! GBB: dlwrad(:) = 4.0_r8 * emiss * sigma * tsrf**3
	! The exponents do not need to be real; but the factor 4 should be real
	
	! Sensible heat flux [W/m2] and temperature derivative [W/m2/K]
	! GBB: Need to multiply by rhoair: J/s/m2 = kg/m3 * J/kg/K * K * m/s
	shflx(:) 	=	cpair * (tsrf - thref) / rah * rhoair	! [W/m2] = [J/kg/K] * [K] / [s/m] * [kg/m3]
	dshflx(:)	=	cpair / rah	* rhoair					! [W/m2/K] = [J/kg/K] / [s/m] * [ kg/m3] 
	
	! Latent heat flux [W/m2] and temperature derivative [W/m2/K]
	! (check if lhflx > water available in snow and soil, in which case limit lhflx
	! to available water; also, if there is snow, don't use soil moisture as a factor)
		! GBB: Need to multiply by rhoair: J/s/m2 = kg/m3 * J/kg/K * K/Pa * Pa * m/s
		! Or you could simplify to: lhflx = rhoair * (qsat - qref) / res 
		! but then you have to convert esat to qsat and desat to dqsat. CLM does this in
		! QSatMod.F90 using:
		!
		!   vp    = 1.0_r8   / (p - 0.378_r8*es)
		!   vp1   = 0.622_r8 * vp
		!   vp2   = vp1   * vp
		!
		!   qs    = es    * vp1             ! kg/kg
		!   qsdT  = esdT  * vp2 * p         ! 1 / K
		!
		! where es is the saturation vapor pressure (Pa), esdT is the temperature
		! derivative (Pa/K) and P is the reference height pressure (Pa). You could copy
		! this entire routine (see my comment below about satvap). See also the CLM4.5
		! technote, equations 5.157, 5.158
		!
		! MML: plan - use qsat instead of esat, by calling CLM function QSat. Modify these
		! equations accordingly (and check units!!!!) 
	
	! Initialize beta = 1.0 (no extra bucket resistance) everywhere. Overwrite with smaller values where appropriate.
	beta(:) = 1.0_r8

	! similarly initialize mml_lnd_effective_res_grc and mml_lnd_res_grc to avoid nans
	atm2lnd_inst%mml_lnd_effective_res_grc = 1.0_r8 !9999.99_r8
	atm2lnd_inst%mml_lnd_res_grc = 1.0_r8 ! 9999.99_r8 
	
	where ( snow <= 0 )
		beta(:) = min ( water/(.75 * bucket_cap) , 1.0_r8 )		! scaling factor [unitless]
		! OH I bet the problem is that I only end up defining beta in places where snow<0 -- hence the nan problem!!! So I should initialize
		! a starting beta matrix where everywhere is 1.0 or something! 
		! add minimum beta value in case water is negative?
		!lhflx(:) 	= cpair / gamma * (esat - eref) / res * beta * rhoair 	! [W/m2] = [J/kg/K] / [Pa/K] * [Pa] / [s/m] * [unitless] * [kg/m3] 
		!dlhflx(:) 	= cpair / gamma * desat / res * beta * rhoair			! [W/m2/K]
		lhflx(:)	= rhoair * lambda * (qsrf - qref) * beta / res 	! [W/m2] = [kg/m3] * [J/kg] * [kg/kg] * [unitless] / [s/m] -> kg/m3 * J/kg * m/s = kg/kg J/s 1/m2 = W/m2
		dlhflx(:) 	= rhoair * lambda * dqsrf * beta / res			! [W/m2/K] = [kg/m3] * [J/kg] * [kg/kg/K] * [unitless] / [s/m] -> kg/m3 * J/kg * 1/K * m/s -> J/s /K /m2 = W/m2/K
		! got here doing unit analysis - make sure this is actually the right equation!!!  
	end where
	
	! make sure beta isn't negative (if neg, set equal to 0)
	where ( beta <= 0.0 )
		beta(:) = 0.0_r8
	end where
	
	where ( snow > 0 ) ! go where there is snow and overwrite the value of lhflx and dlhflx
		!lhflx(:) = cpair / gamma * ( esat - eref ) / res * rhoair			! [W/m2]
		!dlhflx(:) = cpair / gamma * desat / res * rhoair					! [W/m2]
		lhflx(:)	= rhoair * lambda * (qsrf - qref) / res 	! [W/m2] = [kg/m3] * [J/kg] * [kg/kg] * [unitless] / [s/m] -> kg/m3 * J/kg * m/s = kg/kg J/s 1/m2 = W/m2
		dlhflx(:) 	= rhoair * lambda * dqsrf / res			! [W/m2/K] = [kg/m3] * [J/kg] * [kg/kg/K] * [unitless] / [s/m] -> kg/m3 * J/kg * 1/K * m/s -> J/s /K /m2 = W/m2/K
	end where

	! Check if we tried to evaporate more water than is available
	! ... probably isn't the sneakiest way to do this... what if dlhflx is <0? then we might
	! be okay - would have to check at end of time step...
	where ( lhflx * dt / lambda > ( water + snow ) ) 	! [W/m2] * [s] / [J/kg] -> W * [s/J] * kg/m2 = kg/m2
		!write(iulog,*)subname, 'MML tried to evaporate more water than there is in snow + water, adjusting accordingly'
		!lhflx(:) = lambda / dt * ( water + snow ) * rhoair					! [W/m2]
		!dlhflx(:) = 0._r8														! [W/m2]
		lhflx(:)	= lambda / dt * ( water + snow )	! [W/m2] = [J/kg] / [s] * [kg/m2] ->  J/s * kg/kg/m2 = W/m2
		dlhflx(:) 	= 0._r8								! [W/m2/K]
	end where
	



	
	! Net flux of energy into soil [W/m2] and temperature derivative [W/m2/K] from the 
	! surface energy imbalance given other fluxes:
	f0(:) 	=	radforc - ( lwrad + lhflx + shflx )							! [W/m2]
	df0(:) 	= 	- ( dlwrad + dlhflx + dshflx )								! [W/m2]
    
    ! lets temporarily save this value out as gsoi (not the real gsoi, but the right "family"
    gsoi(:) = f0							! [W/m2]
     
     
	! -------------------------------------------------------------
	! Initial pass at soil temperatures
    ! -------------------------------------------------------------
    
    ! Initial change in soil temperatures = 0
    dtsoi(:,:) = 0.0_r8 ! see if this helps?
    
    ! -------------------------------------------------------------
    ! Set up tri-diagonal matrix
    
    ! surface
    i = 1
    
    aa(:,i) = 0.0_r8
    cc(:,i) = -soil_tkh(:,i) / ( soil_z(:,i) - soil_z(:,i+1) )
    bb(:,i) = soil_cv(:,i) * soil_dz(:,i) / dt - cc(:,i) - df0
    dd(:,i) = -soil_tkh(:,i) * ( tsoi(:,i) - tsoi(:,i+1) ) / ( soil_z(:,i) - soil_z(:,i+1) ) + f0
    
    ! layers 2 to nsoi-1
    dummy = mml_nsoi - 1
    do i = 2, dummy
    	aa(:,i) = -soil_tkh(:,i-1) / ( soil_z(:,i-1) - soil_z(:,i) )
    	cc(:,i) = -soil_tkh(:,i) / ( soil_z(:,i) - soil_z(:,i+1) )
    	bb(:,i) = soil_cv(:,i) * soil_dz(:,i) / dt - aa(:,i) - cc(:,i)
    	dd(:,i) = soil_tkh(:,i-1) * ( tsoi(:,i-1) - tsoi(:,i) ) / (soil_z(:,i-1) - soil_z(:,i)) &
    			- soil_tkh(:,i) * (tsoi(:,i) - tsoi(:,i+1)) / (soil_z(:,i) - soil_z(:,i+1))
    end do 
    
    ! Bottom soil layer
    i = mml_nsoi
    aa(:,i) = -soil_tkh(:,i-1) / (soil_z(:,i-1) - soil_z(:,i))
    cc(:,i) = 0.0_r8
    bb(:,i) = soil_cv(:,i) * soil_dz(:,i) / dt - aa(:,i)
    dd(:,i) = soil_tkh(:,i-1) * (tsoi(:,i-1) - tsoi(:,i)) / (soil_z(:,i-1) - soil_z(:,i))
    
    ! ----------------------------------------------------------
    ! Begin forward (upward) sweep of tridiagonal matrix from layer N to 1
    
    ! Bottom soil layer
    i = mml_nsoi
    ee(:,i) = aa(:,i) / bb(:,i)
    ff(:,i) = dd(:,i) / bb(:,i)
    
    ! Layers nsoi-1 to 2
    dummy = mml_nsoi-1
    do i = dummy, 2, -1
    	den = bb(:,i) - cc(:,i)*ee(:,i+1)
    	ee(:,i) = aa(:,i) / den
    	ff(:,i) = (dd(:,i) - cc(:,i)*ff(:,i+1)) / den
    end do
    
    ! Complete tridiagonal sol'n to get initial temperature guess for top soil layer
    i = 1
    num = dd(:,i) - cc(:,i) * ff(:,i+1)
    den = bb(:,i) - cc(:,i) * ee(:,i+1)
    tsrf = tsoi0(:,i) + num/den
    
	
	!write(iulog,*)subname, 'MML new tridiagonal solver IS being used'
	
	! -------------------------------------------------------------
    ! Snow accounting: 
    ! if tsrf>freezing and there is snow on the ground, melt some snow!
    ! -------------------------------------------------------------
        	
    !t_to_snow(:) = soil_cv(:,1) * soil_dz(:,1) / hfus	! factor to convert a change in temperature to snow melt
    
    ! how much snow can we melt given the temperature? 
    snow_melt = 0.0_r8
    !where ( snow > 0.0_r8 .and. tsrf > tfrz) snow_melt(:) = (tsrf(:) - tfrz) * den(:) * t_to_snow(:)
    
    ! Maximum snow melt RATE based on temperature above freezing:
    ptl_snow_melt(:) = max(0.0 , (tsrf(:) - tfrz) * den(:) / hfus)
    !where ( snow > 0.0_r8 .and. tsrf > tfrz) snow_melt(:) = (tsrf(:) - tfrz) * den(:) / hfus 
    
    ! Maximum melt RATE is the rate it would take to melt all the snow that is currently present:
    max_snow_melt(:) = snow / dt
    
    ! Set actual snow melt RATE to either the total the potential (if enough snow is present) or the total (if enoguh energy is present)
    snow_melt(:) = min( max_snow_melt(:) , ptl_snow_melt(:) )
    
    ! Energy flux associated with realized snow melt
    gsnow(:) = snow_melt(:) * hfus		! [kg/m2/s]*[J/kg] = [J/s/m2] = [W/m2]
    
    ! Recalculate melt based off how much snow is actually present (can't melt more
    ! than what is actually present)
    ! If we have more energy than snow to melt, update surface temperature accordingly
    !where ( snow > 0.0_r8 .and. snow_melt > 0.0_r8 .and. snow_melt <= snow ) tsrf(:) = tfrz	! where snow_melt < snow, temperature stays at freezing
    !where ( snow > 0.0_r8 .and. snow_melt > 0.0_r8 .and. snow_melt > snow )
    !	snow_melt(:) = snow	! melt all available snow
    !	tsrf(:) = tsoi(:,1) + (num(:) - snow_melt(:)/t_to_snow(:))/den(:)
    !end where
        
    ! Update snow and water buckets accordingly -> convert to water units, not rates
    snow(:)  = snow  - snow_melt*dt			! [kg/m2] = [kg/m2] - [kg/m2/s]*[s]
    water(:) = water + snow_melt*dt
       		
    ! Update surface temperature to reflect snow melt:
    !	If there is no snow melt, tsoi(1) = tsrf as above, unmodified
    !	While snow is actively melting, tsrf should be tfrz
    ! 	If snow melt was less than the total energy, tsrf should be > trfz but less tahn tsrf above
    tsoi(:,1) = tsoi(:,1) + (num - gsnow) / den;
    dtsoi(:,1) = tsoi(:,1) - tsoi0(:,1)
    
    
  	! -------------------------------------------------------------
  	! Complete the tri-diagonal solver for soil temperature given we now know the 
  	! surface temperature after snow melting
    ! -------------------------------------------------------------	
        	
    !dtsoi(:,1) 	= tsrf(:) - tsoi(:,1)	! save change in top soil layer
    !tsoi(:,1) 	= tsrf(:)					! update top soil layer to be surface temperature
    
    !------ Complete tri-diagonal solver (downwards sweep)
	do i = 2,mml_nsoi
		dtsoi(:,i) = ff(:,i) - ee(:,i)*dtsoi(:,i-1)
		tsoi(:,i) = tsoi(:,i) + dtsoi(:,i)
	end do
    
    !dummy = mml_nsoi - 1
    !do i = 1, dummy
    !    dtsoi(:,i+1) = dp(:,i) + cp(:,i)*dtsoi(:,i)	! ah, this hsould have been i+1
    !    tsoi(:,i+1) = tsoi(:,i+1) + dtsoi(:,i+1) 	! old tsoi + dtsoi
    !end do
	
	
	! -------------------------------------------------------------
    ! Update surface energy fluxes based on the change in surface temperature
    ! -------------------------------------------------------------	
		
	lwrad(:) = lwrad + dlwrad * dtsoi(:,1)
	lhflx(:) = lhflx + dlhflx * dtsoi(:,1)	! if lhflx = snow+water, dlhflx = 0
	shflx(:) = shflx + dshflx * dtsoi(:,1)
	! and the ground energy flux:
	gsoi(:)	 = f0 + df0 * dtsoi(:,1)
	
	! split energy flux into ground into flux into soil (gsoi) and snow (gsnow)
	gsoi(:) = gsoi(:) - gsnow(:)
	!gsoi(:)  = gsoi - snow_melt / dt * hfus
	!gsnow(:) = snow_melt / dt * hfus
	
	
	! Energy conservation check:
	! Sum change in energy (W/m2)
	edif(:) = 0._r8
	do i = 1,mml_nsoi
		edif(:) = edif(:) + soil_cv(:,i) * soil_dz(:,i) * ( tsoi(:,i) - tsoi0(:,i) ) / dt
	end do
	! Energy conservation check:
	err(:) = 0._r8
	err(:) = edif(:) - gsoi(:)
	do g = begg,endg
		if ( abs( err(g) ) > 1.0e-06 ) then
			write(iulog,*)subname, 'MML ERROR: Soil temperature energy conservation error: pre-phase change'
			call endrun(msg=errmsg(__FILE__, __LINE__))
		end if
	end do
	
	! Maybe should be checking lhflx HERE for if it is larger than water+snow
	
	
	lwup(:) = lwup + lwrad	! reflected longwave (0 at the moment) plus sigma*T^4
	
	
	! -------------------------------------------------------------
	! TO DO:
	! If lhflx < 0 and the total amount of water the land tries to suck out of the atmosphere is
	! larger than the total water available in the lowest level of the atmosphere, cap the negative LHFLX
	! at the amount of water in the atm_bot and put the excess energy into SHFLX (cam has a check
	! that does this (qneg4.f90)
	
	! check 1: if evap*dt > water + snow at this point, take excess and put into sensible heat flux?
	do g = begg, endg
		if ( lhflx(g) * dt / lambda(g) > (water(g) + snow(g)) ) then
	!where ( lhflx * dt / lambda > (water + snow) )
			temp(g) = lhflx(g) - (water(g) + snow(g)) * lambda(g) / dt	!excess energy that we don't have water for
			lhflx(g) = lhflx(g) - temp(g)	! remove the excess from lh
			shflx(g) = shflx(g) + temp(g)	! give it to shflx 	...  ask Gordon about a better way to do this...
			write(iulog,*)subname, 'MML Warning: lhflx > available water; put excess in shflx'
	!end where 
		end if	! put in an if loop just so I could get it to write the warning
	end do
	
	
        ! MML 2021.09.13: move update of evap (in water units) to AFTER the lh/sh check - otherwise lh and evap won't match (once put into proper units)
        
        ! LHFLX in water units [kg/m2/s = mm/s]
        ! update evap(g) 
        !evap(:) = lhflx * dt / lambda
        evap(:) = lhflx / lambda        ! kg/m2/s or mm/s, NOT times dt!!!!


	! -------------------------------------------------------------
	!	Check that dew doesn't exceed water available in lowest atm level
	! -------------------------------------------------------------
	! check 2: if evap*dt < 0 and requires more water than is available in the bottom of the atmosphere,
	! that is bad... the atmosphere corrects for it, but I want the atm and land to be self-consistent...
	! TODO STILL!
	! GBB: CLM does not do this
	!
	! MML: implement a check for this (go back to CAM QNEG3 OR QNEG4 to check how CAM does it)
	!	Then limit the CLM LHFLX to whatever CAM is going to adjust it to. Also, print out how
	!	big that energy difference is and save it somewhere - it'll be big in the first couple
	! 	of time steps, but I'm not sure how big/negligible it is after the model is sort of spun
	! 	up. Gordon said there was O(1) W/m2 of energy that sort of gets lost in the coupled 
	!	model - I'm curious if this contributes to that, or if this is totally negligible once
	!	the models spins up. 
	! 	(What CAM does is takes the excess energy that was in LHFLX (but there isn't enough water available
	! 	in the lower level of the atmosphere for) and adds it to the SHFLX, so its still conserving ENERGY
	! 	(ie shouldn't be a source of an energy leak), but its changing the PATHWAY the energy takes.
	
!	! Method:
!	! Following that of the CAM routine qneg4.F90 in cam/src/physics
!	! 
!	! Compute maximum available water in lowest level of the atmosphere (in lhflx units), 
!	! and compare to the lhflx calculated above. If dew formation calculated in CLM is larger
!	! than the total available water in the lowest atm level, set dew formation (neg lhflx)
!	! equal to the amount of water in the lowest atm level, and put excess energy into the shflx.
!	! To do this, we need to calculate how much water is available in the lowest atm level.
!	! Since I don't actually have the pressure bounds for the lowest atm level, I'll approximate
!	! them as follows: pbot = pressure in midpoint of lowest layer; psrf = pressure at bottom boundary.
!	! dpbot = 2*(psrf-pbot) (ie the width in pressure units is 2x the difference between the surface
!	! and the midpoint). 
!	! I'll also assume that the minimum water content of the lowest atm level is 0 (the 
!	! qneg4 fn has a qmin term, but I don't know what else to set it to without passing that
!	! value out of cam...
!	!
!	! Basis for using dp = psrf - pbot is from equation _______ in the cam (4?) tech note...
!	! read more into this, not sure I'm doing it right...
!	
!	! First, save initial lhflx and evap demand to variables to be written out (before modification)
!	q_demand = evap		! pos upward, neg downward, so need q_demand - q_avail to be negative for there to be a problem (ie demand is too large)
!	lh_demand = lhflx
!	sh_init   = shflx
!	
!	! Set constants used by these equations:
!	!grav = 9.81_r8	! m/s2
!	ztodt = 2*dt	! 2 delta t
!	
!	A2 = 0._r8
!	B2 =  0.985112190246582_r8
!	P0 = 100000
!	
!	p2 = A2*P0 + B2*psrf
!	!dpbot = 2 * (psrf - pbot)	! [Pa] , thickness in pressure of lowest level
!	dpbot = psrf - p2
!	write(iulog,*)subname, 'MML compare delta p: psrf-p2 = ', dpbot(begg), & 
!						' while 2*(psrf-pbot)=',2*(psrf(begg)-pbot(begg))
!	
!	srfrpdel = 1 / dpbot		! [1/Pa]
!	
!	! start off assuming there is no excess dew formation:
!	is_qneg = 0
!	
!	! set minimum value of water in lowest atm q_min = 0 [kg/kg] 
!	q_min = 0._r8	! [kg/kg], minimum amount of water allowed in the lowest atm level is 0 kg/kg (stop if going negative)
!	
!	! q_excess []    =    evap [kg/m2/s] - (qmin - qbot [kg / kg] ) / ( ztodt [ s ] * grav [ m/s2 ] * srfrpdel [ 1/ pa =  m s2 / kg]  )
!	! 				=	[kg/m2/s]	- 	([])/[m/s * m s2 / kg ]
!	!				=	[kg/m2/s]	-	1/[m2 s / kg]
!	!				=	[kg/m2/s]	- [kg/m2/s]
!	!				= 	[kg/m2/s]
!	q_avail = (q_min - qbot) / (ztodt * grav * srfrpdel)
!	q_excess = evap + q_avail 
!	! pos upward, neg downward, so need q_demand - q_avail to be negative for there to be a problem 
!	
!	is_qneg = any(q_excess < 0)	! set is_qneg = 1 if any of the points in q_excess are > 0, ie, if dew demand (evap) was larger than atm water (q_avail*const) at any of the tested points
!	
!	! into lh units:
!	! lh_excess = evap [kg/m2/s] * lambda [ J / kg ] = J/m2/s = W*s / m2 / s  = W / m2
!	!			= [kg/m2/s]	* [ N * m / kg]
!	!			=	[kg/m2/s] * [kg m2 / s2 / kg]
!	!			= 	[kg/s] * [ / s2]
!	!			= 	[ kg / s3 ]	* [m2/m2]		! W = kg m2 / s3
!	!			=	[kg m2 /s3] / [m2]
!	!			= 	[W/m2]
!	lh_excess 	=	q_excess * lambda	! [W/m2]
!	
!	! Take q_excess away from evap only if q_excess is negative
!	do g = begg, endg
!	
!		if ( q_excess(g) .lt. 0 ) then ! If we have negative q_excess, we demanded too much h20 from atm
!			write(iulog,*)subname, 'MML WARNING: dew demand by land exceeds water avaialable in lowest amt level at point g=',g, &
!									'evap_demand (pos upward) = ',q_demand(g),' q_available = ', q_avail(g)
!			! in this case, lh_excess would be negative (check:)
!			write(iulog,*)subname, 'MML WARNING: dew demand lh_excess is negative? lh_excess = ',lh_excess(g)
!			! so add lh_excess to lhflx, this should reduce lhflx
!			lhflx(g) = lhflx(g) + lh_excess(g)
!			
!			! do the same with evap (where q_excess is negative, so we should add)
!			evap(g) = evap(g) + q_excess(g)
!			
!			! put energy lost from lhflx into shflx:
!			shflx(g) = shflx(g) - lh_excess(g)
!		
!		end if
!	
!	end do
!	
!	
!	! Check that q_demand (initial evap) = evap (modified) + q_excess 
!	! and that initial lhflx+shflx = new lhflx+shflx
!	epc(:) = 0.0 ! for evap check
!	epc = q_demand - (evap + q_excess)
!	dew = 0
!	do g = begg, endg
!		if ( abs(epc(g)) .gt. 1.0e-06 ) then
!			
!			dew = 1
!		end if
!	end do
!	
!	if (dew .eq. 1) then
!		write(iulog,*)subname, 'MML dew water conservation error '
!		! send a warning that dew formation was larger than atm water availability
!		write(iulog,*)subname, 'MML Warning: lhflx > available water in lowest level of atm!'
!		write(iulog,*)subname, 'MML Warning: initial lhflx = ', lhflx(endg)
!		write(iulog,*)subname, 'MML Warning: available atm lhflx = ', lh_avail(endg)
!		write(iulog,*)subname, 'MML Warning: excess energy in non-existant dew:', lh_excess(endg)
!		write(iulog,*)subname, 'MML Warning: initial shflx = ', shflx(endg)
!		!call endrun(msg=errmsg(__FILE__, __LINE__))
!	end if
!	
!	epc(:) = 0.0 ! for evap check
!	epc = (lhflx + shflx) - (lh_demand + sh_init)
!	dew = 0 
!	do g = begg, endg
!		if ( abs(epc(g)) .gt. 1.0e-06 ) then
!			
!			dew = 1
!		end if
!	end do
!	
!	if (dew .eq. 1) then
!		write(iulog,*)subname, 'MML turbulent energy (lhflx + shflx) dew conservation error'
!		! send a warning that dew formation was larger than atm water availability
!		write(iulog,*)subname, 'MML Warning: lhflx > available water in lowest level of atm!'
!		write(iulog,*)subname, 'MML Warning: initial lhflx = ', lhflx(endg)
!		write(iulog,*)subname, 'MML Warning: available atm lhflx = ', lh_avail(endg)
!		write(iulog,*)subname, 'MML Warning: excess energy in non-existant dew:', lh_excess(endg)
!		write(iulog,*)subname, 'MML Warning: initial shflx = ', shflx(endg)
!		!call endrun(msg=errmsg(__FILE__, __LINE__))
!	end if
	
	
	
	! -------------------------------------------------------------
	! Update fsns and flns 
	fsns = fsds - fsr
	! compare to sw_abs, should be the same. Put in diag3_1d
	!diag3_1d = sw_abs
	
	flns = lwdn - lwup
	


	! -------------------------------------------------------------
	! Adjust soil temperatures for phase change (freezing/thawing in soil)
    ! -------------------------------------------------------------	
	
	! have to translate that function first :p
	! returns new tsoi and epc, where epc is the energy used in phase change [W/m2]
	epc(:) = 0.0 ! for now
	
	call phase_change (begg, endg, tsoi, soil_cv, soil_dz, &
  							soil_maxice, soil_liq, soil_ice, &
  							mml_nsoi, dt, hfus, tfrz, epc &
  							!diag1_1d, diag1_2d, diag2_2d, diag3_2d			& ! temporary diagnostics
  							)
	
	! -------------------------------------------------------------
  	! Check soil temperature energy conservation
 	! -------------------------------------------------------------	
	edif(:) = 0.0		! change in energy in each layer
	do i = 1, mml_nsoi
		edif(:) = edif(:) + soil_cv(:,i) * soil_dz(:,i) * (tsoi(:,i) - tsoi0(:,i)) / dt
	end do
	
	err(:) = edif(:) - gsoi(:) - epc(:)	! not counting gsnow here, because it didn't heat/cool soil
	
	do g = begg, endg
		if ( abs(err(g)) .gt. 1.0e-06 ) then
			write(iulog,*)subname, 'MML Soil Temperature Conservation Error :( at g = ', g, &
								'err(g) = ', err(g), ', edif(g) = ', edif(g),', gsoi(g) = ', gsoi(g)
			call endrun(msg=errmsg(__FILE__, __LINE__))
		end if
	end do
		
	! -------------------------------------------------------------
    ! Bucket hydrology!
    ! Remove water that evaporated via LHFLX from ye-old water and snow buckets
    ! Also add rain/snow falling in from the great-big-sometimes-blue sky
    ! Then calculate runoff if the bucket overflowed 
    !
    ! Ask Gordon - should I be raining into the bucket at the start of the time step?
    ! then let the bucket exceed capacity, do evaporation, and only if there is excess water
    ! at the end of the time step send it to runoff? 
    ! (right now, I'm raining after LHFLX is calculated, so if it was dry then rains,
    ! we have small lhflx, but it could catch up next time step...
    ! ... probably doesn't matter much on the monthly mean scale, but if doing it one
    ! way vs the other results in wibbly-wobbly surface fluxes from time step to time 
    ! step which can be avoided, should do it right... 
    !
    ! GBB: This is how I would do it (calculate latent heat flux on current soil
	! water) and then update the soil water. See what GFDL did.
    ! -------------------------------------------------------------	
        	
    !write(iulog,*)subname, 'MML welcome to bucket hydrology land!'
        	
    ! If there is snow on the ground, sublimate that to get lhflx
    ! If there isn't enough snow to accomodate evap(g) when there is snow, steal it from 
    ! the water bucket (without accounting for hvap or soil wetness or anything like that - 
    ! treating the snow like it has a magic straw into the soil pool)
    ! If there isn't snow, take the water in evap(g) right from the soil water bucket
    
    
    !------------------------------------
    ! Rain into buckets 
    
    ! (should I do this at the start of the time step? would up the amount of lh possible...)
    water = water + mms2kgm * prec_liq			! water in bucket [kg/m2]
    snow  = snow  + mms2kgm * prec_frz			! snow in bucket  [kg/m2]
    
            	
    ! -------------------------------------------------------------	
    ! Evaporation
    
    ! shouldn't ever be in a case where evap > snow + water, it checks that when calculating lhflx
    ! though its possible if lhflx was close to snow + water, that when we update with dTsrf, it goes negative... hmm...
    ! (allow it for now?) 
    
    ! Snow Evaporation:
    snow0  = snow
    water0 = water
    
    where (snow0 > 0 .and. evap*dt <= snow0)
    	! where snow is enough to cover all evaporation, take lhflx out of snow bucket
    	snow(:) = snow0(:) - evap(:)*dt	! here I need to say evap*dt to get kg/m2 not kg/m2/s
    	! NOTE: IF lhflx < 0, then evap < 0, so this will ADD snow to snow bucket (sucking water out of atm)
    end where
   
    ! MML 2021.09.21: changed from using snow to using snow0 in the where statments, otherwise I'm going to evaproate twice, aren't I?  
    where (snow0 > 0 .and. evap*dt > snow0)
    	! where snow isn't enough to cover all evaporation
    	
    	! steal excess water we need from soil bucket
    	wat2snow(:) = evap*dt - snow0
    	! remove wat2snow from water bucket
    	water(:) = water0 - wat2snow			! POSSIBLE that this could go negative at one time step, but shouldn't blow up
    	! give snow wat2snow and remove evap (should equal zero)
    	snow(:) = snow0 + wat2snow - evap*dt

		! NOTE: IF lhflx < 0, then evap < 0, so this will ADD water to the bucket (sucking it out of the atmosphere)
   		! 		... shouldn't actually happen in this case b/c evap*dt < 0 shouldn't also be > snow
    end where
    
    ! Snow-free Evaporation:
    
    where (snow0 <= 0 )
    	water(:) = water0 - evap*dt
    	! NOTE: IF lhflx < 0, then evap < 0, so this will ADD water to the bucket (sucking it out of the atmosphere)
    end where
    
    ! Check water and snow buckets
    do g = begg, endg
    	! Check that snow = 0 (ish) if we were supposed to evaporate it all
    	! if   evap more than init snow    .and.    still have snow    -> error
    	if ( evap(g)*dt > snow0(g) .and. abs(snow(g)) > 1.0e-06) then
    		write(iulog,*)subname, 'MML WARNING evaporation - snow should be 0, and its not! snow = ', snow(g)
    	end if
    	! Check if water went negative, if so, say something!
    	if ( water(g) < 0.0_r8 ) then
    		! changed from < 0 since it was tripping with values of -6e-18 ... 
    		write(iulog,*)subname, 'MML WARNING evaporation - water(g) < 0; water(g) = ', water(g), 'setting water(g) = 0.0'
                water(g) = 0.0_r8
    	end if
    	! Check if snow went negative, if so, say something!
    	if ( snow(g) < 0.0_r8 ) then
    		write(iulog,*)subname, 'MML WARNING evaporation - snow(g) < 0; snow(g) = ', snow(g),', setting snow(g)=0.0'
                snow(g) = 0.0_r8
    	end if
    end do   

    
    
    
    !------------------------------------
    ! Runoff: check if bucket overflowed
    
    where (water > bucket_cap)
    	runoff = water - bucket_cap	! excess h20
    	water = bucket_cap
    end where
    
	! Check we didn't let snow or water go negative
	do g = begg, endg
       	if ( snow(g) < 0.0_r8  .or. water(g) < 0.0_r8  ) then
			write(iulog,*)subname, 'MML WARNING snow or water bucket went negative, uhoh (after runoff)'
		end if
		
! 		if( isnan(snow(g)) ) then
!     		write(iulog,*)subname, 'MML ERROR:snow is a nan \n', &
!     		!call endrun(msg=errmsg(__FILE__, __LINE__))
!     	end if
!     
	end do
	
	!---------------------------------------
	! Snow build-up: if there is too much snow, send extra to runoff 
	! 	(note, not using any energy to melt it or anything - its getting sent as ice to runoff)
do g = begg, endg 
	if (snow(g) > snowcap ) then  
		runoff(g) = runoff(g) + (snow(g) - snowcap) 
		snow(g) = snowcap 
	end if 
	 
	if (snow(g) <  0.0) then 
		write(iulog,*)subname, 'MML WARNING snow went negative after implementing snow cap... ' 
		snow(g) = 0.0 
	end if 
	 
	if (snow(g) > snowcap) then 
		write(iulog,*)subname, 'MML WARNING snow exceeds snow cap after implementing snow cap... ' 
	end if 
	 
       if (water(g) <  0.0) then
                write(iulog,*)subname, 'MML WARNING water went negative, set to zero...  '
                water(g) = 0.0
        end if
end do 
	
	! -------------------------------------------------------------
    ! Now what now what now what? This is so exciting :)
    ! -------------------------------------------------------------
!           
    ! -------------------------------------------------------------
    ! Check energy balance
    ! -------------------------------------------------------------
          ! radforc = swabs + lwabs
          ! lwup = lw_reflected + lwrad
    err = radforc - (lwup + lhflx + shflx + gsoi + gsnow)
    do g = begg, endg
    	if( abs(err(g)) > 1.0e-06) then
    		write(iulog,*)subname, 'MML ERROR: Not conserving energy (surface fluxes) \n', &
    					'err = ', err(g), &
    					'\n radforc = ', radforc(g), &
    					'\n lwup = ', lwup(g), &
    					'\n lhflx = ', lhflx(g), &
    					'\n shflx = ', shflx(g), &
    					'\n gsoi = ', gsoi(g), &
    					'\n gsnow = ', gsnow(g)
    		call endrun(msg=errmsg(__FILE__, __LINE__))
    	end if
    end do
    	
    
    
 
     
     
     
     ! -------------------------------------------------------------
     ! Update qs (surface specific humidity - need it for next round's MO calculations)
     ! -------------------------------------------------------------
     
     ! instead of direct calculation, re-evaluate QSat on the new surface temperature to get qsrf
     qsrf(:) = qref + evap*dt * res / (dt * rhoair)
     ! Gordon says leave it with the above equation (the below is the inversion to calculate it...)
     !do g = begg, endg
     !	call QSat (tsrf(g), pref(g), esrf(g), desrf(g), qsrf(g), dqsrf(g))
     !end do
     ! updates qsrf for next time step using current tsrf
     
     
     ! do I need to do this if I'm using qsrf for my lhflx now? probably... or move the call
     ! to QSat to before the MO calculation... that sounds better... except I still want to
     ! print out qsrf to the h0 file. Hmm. Well, they SHOULD be consistent, right? Actually no,
     ! because I calculate the first pass at qsrf using tsrf from last time step, and I want to
     ! get the updated version to pass up to the atm. 
     ! Go through after implementing the call to QSat instead of satvap and make sure I'm
     ! being self-consistent within the module re: qsrf that I'm using / passing out / using on next time step.
     
     ! Trying to follow CLM 4.5 tech note, but there they're using uatm - us (u surface? not u star?) AND
     ! vatm and vs, which I don't think I've got...
     ! Definitely check these...
     !
     ! GBB: 
	 ! atm2lnd_inst%forc_wind_grc(g)  = sqrt(atm2lnd_inst%forc_u_grc(g)**2 + atm2lnd_inst%forc_v_grc(g)**2)
	 ! but for taux and tauy you want to preserve the zonal and meridonal components
	 ! taux = -rhoair * atm2lnd_inst%forc_u_grc(g) / ram
	 ! tauy = -rhoair * atm2lnd_inst%forc_v_grc(g) / ram
     taux = -rhoair * (atm2lnd_inst%forc_u_grc - 0._r8) / ram		! [kg/m/s2] = [kg/m3] * [m/s] / [s/m]	
     tauy = -rhoair * (atm2lnd_inst%forc_v_grc - 0._r8) / ram		! [kg/m/s2] = [kg/m3] * [m/s] / [s/m]	
     ! the - 0._r8 should be removed later, this is to remind myself I'm saying u_ref - u_srf, where u_srf = 0 by def'n
     
     
     ! -------------------------------------------------------------
     ! Set T2m, q2m, and u10m based off ustar and tstar 
     ! -------------------------------------------------------------
     
     ! Loop over g, post-clm meeting attempt at u10, q2, t2
     
     do g = begg,endg
     
     	! U10m
     	! u2 - u1 = ustar/vkc * ( ln((z2-d)/(z1-d)) - psi_m((z2-d)/L) + psi_m((z1-d)/L))
     	! u1 -> u_10m, u2 -> uref
     	! => u10 = uref - ustar/vkc * ( ln((zref - d)/(z0m+10)) - psi_m((zref-d)/L) + psi_m((z0m+10)/L)  )
     	atm2lnd_inst%mml_out_uref10m_grc(g) = uref(g) - ustar(g) / vkc * &
     			( log((zref(g) - h_disp(g))/(z0m(g) + 10)) - &
     			  psi_m((zref(g) - h_disp(g))/obu(g)) + psi_m((z0m(g) + 10)/obu(g)) )
     	
     	! T2m
     	! Note: this calculation is for theta_2m, NOT T_2m... convert! 
     	! conversion: theta = T * (p/p0) ^ (R/cp)
     	! need p2m, theta2m, to get T2m... for now, assume p0 ~ p2m so these are roughly equal
     	!
     	! GBB: Use thref, because this is used for a flux calculations
     	atm2lnd_inst%mml_out_tref2m_grc(g) = thref(g) - tstar(g) / vkc * &
     			( log((zref(g) - h_disp(g))/(z0h(g) + 2)) - &
     			  psi_h((zref(g) - h_disp(g))/obu(g)) + psi_h((z0h(g) + 2)/obu(g))  )
     	! MML: check this again (T's vs Theta's in all the places), make sure I've got 
     	! them set right (is tref2m a temperature or a potential temperature? Does this eq'n 
     	! GIVE a temperature or a potential temperature? etc... check!)
     			  
     	
     	! q2m
     	! GBB: This is the same as for temperature: use psi_h and z0h
        atm2lnd_inst%mml_out_qref2m_grc(g) = qref(g) - qstar(g) / vkc * &
     			( log((zref(g) - h_disp(g))/(z0h(g) + 2)) - &
     			  psi_h((zref(g) - h_disp(g))/obu(g)) + psi_h((z0h(g) + 2)/obu(g))  )
     	! MML problems with passing a nan here to the coupler - check if nan, if so print more info
     	
     	if( isnan(atm2lnd_inst%mml_out_qref2m_grc(g)) ) then
    		write(iulog,*)subname, 'MML ERROR: mml_out_qref2m_grc is a nan \n', &
    					'err = ', err(g), &
    					'\n qref = ', qref(g), &
    					'\n qstar = ', qstar(g), &
    					'\n zref = ', zref(g), &
    					'\n h_disp = ', h_disp(g), &
    					'\n z0h = ', z0h(g), &
    					'\n obu = ', obu(g)
    		call endrun(msg=errmsg(__FILE__, __LINE__))
    	end if
     	
     	
     
     	! GBB: Did you check this?
     	! MML: at one point, bad things happened. But now that I've got the equations right,
     	!		check again and compare! SHOULD give the same answers. CLM interpolates from surface
     	!		values instead of ref height values. Check that they're the same.
     	! compare to using surface values us = 0, ts, qs, and u*, t*, q*:
     	!
     	! still crashing the model when I include these (where it writes h0 files, interestingly enough).
     	! ... so, is it the equations that are wrong, or the fact that I'm overwriting a predefined
     	! variable diag?_1d that makes it unhappy? Check both.
     	
     	!u10m:
!		diag1_1d(g) = ustar(g) / vkc * &
!				( log( (z0m(g) + 10) / z0m(g) ) &
!				- psi_m( (z0m(g) + 10) / obu(g) ) + psi_m( z0m(g) / obu(g) ) )
!     	
!     	! t2m
!     	! also compare to using tsrf...
		! this differes from thref version -> maybe I need to calculate a theta_srf using a reference pressure p0? 
		! ( theta = T * (p0 / p ) ^ (R/cp) where R/cp ~ 0.286 for air... ... use a p0 = 1000 hPa = 100000 Pa 
		theta_srf(g) = tsrf(g) * (100000._r8 / pref(g)	)**(0.286_r8)
		! (using pref which is pbot, assuming psrf ~ pref)
		
!		diag2_1d(g) = tsrf(g) + tstar(g) / vkc * &		! aha! I bet that was the problem - I never defined theta_srf
		!diag2_1d(g) = theta_srf(g) + tstar(g) / vkc * &
!				( log( (z0h(g) + 2) / z0h(g) ) &
!				 - psi_h( (z0h(g) + 2) / obu(g) ) + psi_h( (z0h(g) ) / obu(g) ) )
				 ! Nope, that didn't fix it, still disagrees with calculation using thref
		
!     
!    	! q2m
!		diag3_1d(g) = qsrf(g) + qstar(g) / vkc * &		
!				( log( (z0h(g) + 2) / z0h(g) ) &
!				 - psi_h( (z0h(g) + 2) / obu(g) ) + psi_h( (z0h(g) ) / obu(g) ) )
				 
		! Still not matching using zref...
    
     	! Answer:
     	! u10 is the same, but t2m and q2m aren't. In particular, t2m is much cooler using the 
     	! surface interpolation vs using the zref interpolation. In q2m, the negative values are 
     	! larger (-.02 in surface interpolation vs -.002 in ref height interpolation).
     	! Look in to why this is. Bad values of qsrf / tsrf? (by bad, I mean inconsistent)
     end do
     
     !*********************************************
     ! Variables needed for use of diagnostics package
     ! fsds
     ! flds
     ! FSDSND	! NIR direct incident shortwave
     ! FSDSNI	! NIR diffuse incident shortwave
     ! FSRND	! NIR direct reflected shortwave
     ! FSRNI	! NIR diffuse reflected
     ! Q2M
     ! RH2M
     ! P-E
     ! 
     
     ! How about I start mapping my variables onto the CLM fields, now that CLM is turned off?
     ! then I won't have to attack the file with nco to get it into a CLM-format for the diagnostics!
     ! (this will only work on fields without any depth dimension).
     
     
     
     !	! save beta out for netcdf
!	atm2lnd_inst%mml_lnd_beta_grc(:) = beta(:)
!        	
!	! and 1/beta * (rs + rah)  1/beta * res , the effective resistnace
!	atm2lnd_inst%mml_lnd_effective_res_grc(:) = res(:) / beta(:)

! save beta out for netcdf
    do g = begg,endg

         atm2lnd_inst%mml_lnd_effective_res_grc(g) = 1.0_r8 ! this is not the actual value, but a nan is showing up in here so for now, just comment out the calculation - should be res/beta
         atm2lnd_inst%mml_lnd_beta_grc(g) = beta(g) 
         atm2lnd_inst%mml_lnd_res_grc(g) = res(g) 

        if (beta(g) < 1.0e-3) then
           atm2lnd_inst%mml_lnd_effective_res_grc(g) = 999999.0_r8
        else
           atm2lnd_inst%mml_lnd_effective_res_grc(g) = res(g) / beta(g)
        end if

        !        atm2lnd_inst%mml_lnd_beta_grc(g) = beta(g) !beta(:)
        !        if(isnan(atm2lnd_inst%mml_lnd_beta_grc(g))) then
        !                atm2lnd_inst%mml_lnd_beta_grc(g) = 1.0e-8_r8 !1.0e36_r8 ! something very small
        !        end if
        !        ! if beta smaller than 0.01 set it larger 
        !        !if(atm2lnd_inst%mml_lnd_beta_grc(g)<0.01) then
        !        !        atm2lnd_inst%mml_lnd_beta_grc(g) = 0.01 ! something very small
        !        !end if
        !        
        !        atm2lnd_inst%mml_lnd_effective_res_grc(g) = res(g) / beta(g) 
        !        if(isnan(atm2lnd_inst%mml_lnd_effective_res_grc(g))) then
        !                atm2lnd_inst%mml_lnd_effective_res_grc(g) = 1.0e-8_r8 !1.!0e36_r8
        !        end if
        !        !if(atm2lnd_inst%mml_lnd_effective_res_grc(g)>10000.) then
        !        !        atm2lnd_inst%mml_lnd_effective_res_grc(g) = 10001.0
        !        !end if
        !        !if(atm2lnd_inst%mml_lnd_effective_res_grc(g)>10000.) then
        !        !        atm2lnd_inst%mml_lnd_effective_res_grc(g) = 10000.0
        !        !end if
        !        
        !        atm2lnd_inst%mml_lnd_res_grc(g) = res(g)
    	!      	if( isnan(atm2lnd_inst%mml_lnd_res_grc(g)) ) then
    	!        	atm2lnd_inst%mml_lnd_res_grc(g) = 1.0e-8_r8 ! 1.e36_r8
     !		    end if
    ! 		    if( atm2lnd_inst%mml_lnd_res_grc(g)>10000. ) then
   ! 	        	atm2lnd_inst%mml_lnd_res_grc(g) = 1.0e-8_r8 !1.e36_r8
   !  		    end if
                
     end do
     
          ! save out res for the netcdf 
     !atm2lnd_inst%mml_lnd_res_grc(:) = res(:)
     
!     do g = begg,endg
!           
!     end do
     !*********************************************
     ! Export my land data to the atmosphere
     ! (send these to lnd2atm)
     
     ! To be sure of what is actually being sent to the atmosphere, see subroutine lnd_export
     ! in /glade/p/work/mlague/cesm_source/cesm1_5_beta05_mml_land/components/clm/src/cpl/lnd_import_export.F90 
    
 
     ! lnd -> atm
     lnd2atm_inst%t_rad_grc = tsrf									! radiative temperature (Kelvin)
     lnd2atm_inst%t_ref2m_grc = atm2lnd_inst%mml_out_tref2m_grc 	! 2m surface air temperature (Kelvin)
     !atm2lnd_inst%mml_lnd_ts_grc = tsrf ! dunno what its saving out now... 
     !lnd2atm_inst%q_ref2m_grc = atm2lnd_inst%mml_out_qref2m_grc 	! 2m surface specific humidity (kg/kg)
     !lnd2atm_inst%u_ref10m_grc = atm2lnd_inst%mml_out_uref10m_grc 	! 10m surface wind speed (m/sec)
     
     lnd2atm_inst%q_ref2m_grc = atm2lnd_inst%mml_out_qref2m_grc 	! 2m surface specific humidity (kg/kg)
     lnd2atm_inst%u_ref10m_grc = atm2lnd_inst%mml_out_uref10m_grc 	! 10m surface wind speed (m/sec)
     
     ! note: mm h20 snow if using rhowat to convert should be the same as kg/m2
     lnd2atm_inst%h2osno_grc = snow / rhowat * 1000	! [kg/m2] / [kg/m3]	* 1000[mm/m]! snow water (mm H2O)
     !lnd2atm_inst%h2osoi_vol_grc									! volumetric soil water (0~watsat, m3/m3, nlevgrnd) (for dust model)
     
     ! MML: albedo (:,:) -> albd is direct, albd(:,1) direct vis, albd(:,2) direct nir
     ! 					 -> albi is diffuse, albi(:,1) diffuse vis, albi(:,2) diffuse nir (I THINK)
     ! GBB: yes
     lnd2atm_inst%albd_grc(:,1) 		= 	alb_vis_dir				! (numrad=1, vis) surface albedo (direct)
     lnd2atm_inst%albd_grc(:,2) 		= 	alb_nir_dir				! (numrad=2, nir) surface albedo (direct)
     
     lnd2atm_inst%albi_grc(:,1) 		= 	alb_vis_dif				! (numrad=1, vis) surface albedo (diffuse)
     lnd2atm_inst%albi_grc(:,2) 		= 	alb_nir_dif				! (numrad=2, nir) surface albedo (diffuse)
     								
     lnd2atm_inst%taux_grc				=	taux					! wind stress: e-w (kg/m/s**2)
     lnd2atm_inst%tauy_grc				= 	tauy					! wind stress: n-s (kg/m/s**2)
     lnd2atm_inst%eflx_lh_tot_grc		= 	lhflx					! total latent HF (W/m**2)  [+ to atm]
     lnd2atm_inst%eflx_sh_tot_grc		= 	shflx					! total sensible HF (W/m**2) [+ to atm]
!     lnd2atm_inst%eflx_sh_precip_conversion_grc			! sensible HF from precipitation conversion (W/m**2) [+ to atm]				
     			! Land group says (a) this is new (sh_precip_converstion) and I can set it to 0 since I don't have multiple levels on my (currently nonexistent) ice sheets
     lnd2atm_inst%eflx_lwrad_out_grc	= 	lwup					! IR (longwave) radiation (W/m**2)
     lnd2atm_inst%qflx_evap_tot_grc		= 	evap					! (mm H2O/s) ! qflx_evap_soi + qflx_evap_can + qflx_tran_veg		
     lnd2atm_inst%fsa_grc				= 	sw_abs									! solar rad absorbed (total) (W/m**2)
     ! MML: not running interactive BGC, set CO2/Methane fluxes to 0
     !lnd2atm_inst%net_carbon_exchange_grc				= 	0._r8					! net CO2 flux (kg CO2/m**2/s) [+ to atm]
     lnd2atm_inst%net_carbon_exchange_grc				= 	0._r8	
     lnd2atm_inst%nem_grc				= 	0._r8					! gridcell average net methane correction to CO2 flux (g C/m^2/s)
     lnd2atm_inst%ram1_grc				= 	ram						! aerodynamical resistance (s/m)
     			! MML: check if it is ram (vs res) that I should be exporting here
     !lnd2atm_inst%fv_grc				= 							! friction velocity (m/s) (for dust model)
     			! MML: should be able to calculate this from MO theory... is this ustar? 
     
     ! Need to put the dust fluxes I read from the .nc file into the right size
     lnd2atm_inst%flxdst_grc			=   dust				! dust flux (size bins)
     !lnd2atm_inst%flxdst_grc			=	0._r8			! (:,ndust) where ndust=4, so I need a 4th dust flux field! and I think the ones I had were wrong...
     			! MML: need some sort of forcing file - see what the aquaplanet people are using
     			! 	currently borrowing the value from CLM by running the whole CLM model first... 
     ! MML: set all these to 0, not running interactive BGC
     !lnd2atm_inst%ddvel_grc				= 0._r8						! dry deposition velocities
     !lnd2atm_inst%flxvoc_grc			= 0._r8						! VOC flux (size bins)
     		! hmm, cesm crashes saying "Attempt to use pointer FLXVOC_GRC when it is not associated with a target"
     		! so for now, skip these ones... (maybe they aren't defined at all when I build as CLM%SP? 
     !lnd2atm_inst%fireflx_grc			= 0._r8						! Wild Fire Emissions
     !lnd2atm_inst%fireztop_grc			= 0._r8						! Wild Fire Emissions vertical distribution top
     !lnd2atm_inst%flux_ch4_grc			= 0._r8						! net CH4 flux (kg C/m**2/s) [+ to atm]
     
     
     ! lnd -> rof (not for now? or all zeros?) 
     !lnd2atm_inst%
     !lnd2atm_inst%
     
     
     ! Save the top-layer soil_cv and soil_tk to the temporary diag fields, to make sure the glacier mask is working 
     diag1_1d(begg:endg) = soil_tk(begg:endg,1)
     diag2_1d(begg:endg) = soil_cv(begg:endg,1)
     ! Excellent, yes, it appears to be working :) 
     
     
     !write(iulog,*)subname, 'MML Map some of my surface variables onto CLM fields for easy running of diagnostics'
     ! (currently doing this in post-processing, just making a new .nc with the CLM-named vars made out of my vars)
     ! Sensible Heat, Latent Heat, Surface T, 2m T, albedo, evaporation rate (kg/m2/s)
     
     
     
     !write(iulog,*)subname, 'MML done simple land. Cleanup time!'
     ! deallocate everything
     ! tidy code: shouldn't have to allocate/deallocate all of these. Instead, 
     ! define them up top as real(r8)	:: myvar(begg:endg)	(assuming I pass in begg and endg, else use bounds%begg)
     ! rather than real, allocatable, (:) then allocate with (begg:endg)
  

   end associate
   

  end subroutine mml_main
  
  
  
 ! TODO!!!
 ! MML: move prescribed value reading from .nc files to a new subroutine
 ! subroutine read_nc_data ( bounds )
 ! 
 ! end subroutine read_nc_data
 
 subroutine nc_import (begg, endg, mml_nsoi, lfsurdat, mon, &
 					albedo_gvd, albedo_svd, albedo_gnd, albedo_snd,&
 					albedo_gvf, albedo_svf, albedo_gnf, albedo_snf,&
 					snowmask, evaprs, &
 					bucket_cap, &
 					!soil_maxice, soil_z, 
 					soil_type, roughness, emiss, glc_mask, dust, &
 					soil_tk_1d, soil_cv_1d, &
 					glc_tk_1d, glc_cv_1d     ) !, &
 					!emiss, glc_mask ) !, soil_tk, soil_cv ) ! combine soil tk and glc tk into "soil_tk" that accounts for ground type

  
   	! For importing .nc data
   	use fileutils   , only : getfil  
   	!use domainMod   , only : domain_type, domain_init, domain_clean
   	use ncdio_pio
   	use clm_varcon      , only : grlnd
  
    implicit none
    
    ! !ARGUMENTS:
    integer, intent(in) 	:: begg
    integer, intent(in) 	:: endg
    integer, intent(in)		:: mml_nsoi
	character(len=256), intent(in) :: lfsurdat	! surface prescribed data file name !MML make flexible length!
	integer, intent(in)		:: mon		! current month
	! Direct Albedos
	real(r8), intent(out)	:: albedo_gvd(begg:endg)
	real(r8), intent(out)	:: albedo_svd(begg:endg)
	real(r8), intent(out)	:: albedo_gnd(begg:endg)
	real(r8), intent(out)	:: albedo_snd(begg:endg)
	! Diffuse Albedos
	real(r8), intent(out)	:: albedo_gvf(begg:endg)
	real(r8), intent(out)	:: albedo_svf(begg:endg)
	real(r8), intent(out)	:: albedo_gnf(begg:endg)
	real(r8), intent(out)	:: albedo_snf(begg:endg)
	! Other land params
	real(r8), intent(out)	:: snowmask(begg:endg)
	real(r8), intent(out)	:: evaprs(begg:endg)
	real(r8), intent(out)	:: bucket_cap(begg:endg)
	!real(r8), intent(out)	:: soil_maxice(:,:)
	!real(r8), intent(out)	:: soil_z(:,:)
	real(r8), intent(out)	:: soil_type(begg:endg)
	real(r8), intent(out)	:: roughness(begg:endg)
	real(r8), intent(out)	:: emiss(begg:endg)
	real(r8), intent(out)	:: glc_mask(begg:endg)
	real(r8), intent(out)	:: dust(begg:endg,4)
	real(r8), intent(out)	:: soil_tk_1d(begg:endg)
	real(r8), intent(out)	:: soil_cv_1d(begg:endg)
	real(r8), intent(out)	:: glc_tk_1d(begg:endg)
	real(r8), intent(out)	:: glc_cv_1d(begg:endg)

    ! !LOCAL VARIABLES:
    character(len=32) :: subname = 'nc_import_sub_mml'
    ! MML allocation variables to read from .nc file
    real(r8), pointer  :: nc_alb_gvd(:) 	=> null()       	! ground albedo read from .nc file
    real(r8), pointer  :: nc_alb_svd(:)     => null()   	! snow albedo read from .nc file
    real(r8), pointer  :: nc_alb_gnd(:)     => null()   	! 
    real(r8), pointer  :: nc_alb_snd(:)		=> null()
    real(r8), pointer  :: nc_alb_gvf(:)     => null()   	! ground albedo read from .nc file
    real(r8), pointer  :: nc_alb_svf(:)     => null()   	! snow albedo read from .nc file
    real(r8), pointer  :: nc_alb_gnf(:)     => null()   	! 
    real(r8), pointer  :: nc_alb_snf(:)		=> null()
    real(r8), pointer  :: nc_snowmask(:)    => null()    	! snow masking depth read from .nc file
    real(r8), pointer  :: nc_evaprs(:)      => null()  	! evap resistance from .nc file
    real(r8), pointer  :: nc_bucket(:)      => null()  	! soil bucket depth from .nc file
    real(r8), pointer  :: nc_ice(:,:)       => null() 	! freezeable water in each soil layer from .nc file
    real(r8), pointer  :: nc_z(:,:)        	=> null()	! depth from surf to each soil layer from .nc file
    real(r8), pointer  :: nc_type(:)        => null()		! soil type from .nc file
    real(r8), pointer  :: nc_rough(:)       => null() 	! roughness length from .nc file
    real(r8), pointer  :: nc_soil_tk(:)     => null() 	! soil thermal conductivity from .nc file
    real(r8), pointer  :: nc_glc_tk(:)      => null() 	! glacier thermal conductivity from .nc file
    real(r8), pointer  :: nc_soil_cv(:)     => null() 	! soil heat capacity from .nc file
    real(r8), pointer  :: nc_glc_cv(:)      => null() 	! glacier heat capacity from .nc file
    real(r8), pointer  :: nc_glc_mask(:)    => null() 	! glacier mask from .nc file
    real(r8), pointer  :: nc_emiss(:)       => null() 	! emissivity (for LW) from .nc file
    real(r8), pointer  :: nc_dust(:)       => null() 	! dust flux (clm5 climatology for now) from .nc file
    ! note: doing allocatable, pointer won't compile, says variable has already been 
    ! assigned the allocatbale tribute ... so does being a pointer encompass being allocatable? 
    ! same error if I do pointer, allocatable instead:
    !
    !/glade/u/home/mlague/cesmruns/simple_land/SLIMv1_11/SourceMods/src.clm/mml_main.F90(361): error #6448: This name has already bee
	! n assigned the ALLOCATABLE attribute.   [NC_ALB_BARE]
    ! real(r8), pointer, allocatable  :: nc_alb_bare(:)           ! ground albedo read from .nc file
	!---------------------------------------^
	! 
	! Just allocatable = says it doesn't have a function ncd_io that works 
	! ... so going back to pointer, but this time lets give it a dummy value.
    
    
   	! real(r8), pointer :: temp_nc(:)	! temporary placeholder for a single month's data
    integer :: k
 	! for getting current month
    integer :: year    ! year (0, ...) for nstep+1
   ! integer :: mon     ! month (1, ..., 12) for nstep+1
    integer :: day     ! day of month (1, ..., 31) for nstep+1
    integer :: sec     ! seconds into current date for nstep+1
    integer :: mcdate  ! Current model date (yyyymmdd)
        
    character(len=256)	:: locfn                ! local file name
    logical           	:: readvar              ! true => variable is on dataset
    type(file_desc_t) 	:: ncid                 ! netcdf id
    !type(file_desc_t)	:: ncid0
    integer			  	:: ier
    
   integer :: i 	! index for hard-coded looping over soil layers
   ! integer :: mml_nsoi	
   
    real(r8)	:: dummy, ival
    
   	
   	ival = 0.0_r8
   	
   	allocate( nc_alb_gvd	(begg:endg)		)		; nc_alb_gvd(:) = ival
    allocate( nc_alb_svd	(begg:endg)		)		; nc_alb_svd(:) = ival
    allocate( nc_alb_gnd	(begg:endg)		)		; nc_alb_gnd(:) = ival
    allocate( nc_alb_snd	(begg:endg)		)		; nc_alb_snd(:) = ival
    allocate( nc_alb_gvf	(begg:endg)		)		; nc_alb_gvf(:) = ival
    allocate( nc_alb_svf	(begg:endg)		)		; nc_alb_svf(:) = ival
	allocate( nc_alb_gnf	(begg:endg)		)		; nc_alb_gnf(:) = ival
    allocate( nc_alb_snf	(begg:endg)		)		; nc_alb_snf(:) = ival
    allocate( nc_snowmask	(begg:endg)		)		; nc_snowmask(:) = ival
    allocate( nc_evaprs		(begg:endg)		)		; nc_evaprs(:) 	= ival
    allocate( nc_bucket		(begg:endg)		)		; nc_bucket(:) 	= ival
    allocate( nc_ice		(begg:endg,mml_nsoi) )	; nc_ice(:,:)	= ival
    allocate( nc_z			(begg:endg,mml_nsoi) )  ; nc_z(:,:)		= ival
    allocate( nc_type		(begg:endg)		)		; nc_type(:) 	= ival
    allocate( nc_rough		(begg:endg)		)		; nc_rough(:) 	= ival
    allocate( nc_soil_tk	(begg:endg)		)		; nc_soil_tk(:) 	= ival
    allocate( nc_glc_tk		(begg:endg)		)		; nc_glc_tk(:) 	= ival
    allocate( nc_soil_cv	(begg:endg)		)		; nc_soil_cv(:) 	= ival
    allocate( nc_glc_cv		(begg:endg)		)		; nc_glc_cv(:) 	= ival
    allocate( nc_glc_mask	(begg:endg)		)		; nc_glc_mask(:) 	= ival
    allocate( nc_emiss		(begg:endg)		)		; nc_emiss(:) 	= ival
	allocate( nc_dust		(begg:endg)		)		; nc_dust(:) 	= ival	! keep overwriting this for each dust bin


!   	if (ier /= 0) then
!       write(iulog,*)subname, 'allocation error MML one of these died...'
!       call endrun(msg=errMsg(__FILE__, __LINE__))
!	end if
   	
    
	!write(iulog,*)subname, 'MML in nc_import subroutine start'
	
	
	!write(iulog,*)subname, 'MML assigned dummy value to nc_data'
	
	!MML: Grab the current model time so we know what month we're in
  !  call get_curr_date(year, mon, day, sec)   ! Actually all I need for now is mon
    !mcdate = year*10000 + mon*100 + day
  	k=mon  ! use k=day to check if its reading in short, 11 day trial runs
    
  	!write(iulog,*)subname, 'MML grabbed the date; mon = ',k
  	
    ! Open .nc file specified by lfsurdat
    call getfil( lfsurdat, locfn, 0 ) 
    call ncd_pio_openfile (ncid, trim(locfn), 0)
    
    !write(iulog,*)subname, 'MML opened the .nc file (I think...)'
    
    
    ! Can't read directly to a pointer variable, instead read to a temporary array that
    ! I've allocated memory to, then set the pointer equal to the temporary array
    
!    write(iulog,*)subname, 'MML failing at .nc reading now, move on with dummy values'
    !write(iulog,*)subname, 'MML begin reading .nc vars'
    
    ! Climatological dust flux
    ! Just in case I can't read in a dust flux, set dust = 0 everywhere to start
    dust(begg:endg,:) = 0.0_r8
    
    ! first dust bin:
    call ncd_io(ncid=ncid, varname='l2xavg_Fall_flxdst1', flag='read', data=nc_dust, &
              dim1name=grlnd, nt=k, readvar=readvar)
    if ( .NOT. readvar) then
		write(iulog,*)subname, 'MML tried to read dust-1, failed ', readvar
    else
	dust(begg:endg,1) = nc_dust
	end if 
	
	
	! second dust bin:
    call ncd_io(ncid=ncid, varname='l2xavg_Fall_flxdst2', flag='read', data=nc_dust, &
              dim1name=grlnd, nt=k, readvar=readvar)
    if ( .NOT. readvar) then
		write(iulog,*)subname, 'MML tried to read dust-2, failed ', readvar
    else
	dust(begg:endg,2) = nc_dust
	end if 
	
	
	! third dust bin:
    call ncd_io(ncid=ncid, varname='l2xavg_Fall_flxdst3', flag='read', data=nc_dust, &
              dim1name=grlnd, nt=k, readvar=readvar)
    if ( .NOT. readvar) then
		write(iulog,*)subname, 'MML tried to read dust-3, failed ', readvar
    else
	dust(begg:endg,3) = nc_dust
	end if 
	

	! fourth dust bin:
    call ncd_io(ncid=ncid, varname='l2xavg_Fall_flxdst4', flag='read', data=nc_dust, &
              dim1name=grlnd, nt=k, readvar=readvar)
    if ( .NOT. readvar) then
		write(iulog,*)subname, 'MML tried to read dust-4, failed ', readvar
    else
	dust(begg:endg,4) = nc_dust
	end if 
	
	
	
    ! Albedo Direct
    call ncd_io(ncid=ncid, varname='alb_gvd', flag='read', data=nc_alb_gvd, &
              dim1name=grlnd, nt=k, readvar=readvar)
    if ( .NOT. readvar) then
		write(iulog,*)subname, 'MML tried to read alb_gvd, failed ', readvar
	end if 
	
	call ncd_io(ncid=ncid, varname='alb_svd', flag='read', data=nc_alb_svd, &
              dim1name=grlnd, nt=k, readvar=readvar)
    if ( .NOT. readvar) then
		write(iulog,*)subname, 'MML tried to read alb_svd, failed ', readvar
	end if 
	
	call ncd_io(ncid=ncid, varname='alb_gnd', flag='read', data=nc_alb_gnd, &
              dim1name=grlnd, nt=k, readvar=readvar)
    if ( .NOT. readvar) then
		write(iulog,*)subname, 'MML tried to read alb_gnd, failed ', readvar
	end if 
	
	call ncd_io(ncid=ncid, varname='alb_snd', flag='read', data=nc_alb_snd, &
              dim1name=grlnd, nt=k, readvar=readvar)
    if ( .NOT. readvar) then
		write(iulog,*)subname, 'MML tried to read alb_snd, failed ', readvar
	end if 
	
	! Albedo Diffuse
	call ncd_io(ncid=ncid, varname='alb_gvf', flag='read', data=nc_alb_gvf, &
              dim1name=grlnd, nt=k, readvar=readvar)
    if ( .NOT. readvar) then
		write(iulog,*)subname, 'MML tried to read alb_gvf, failed ', readvar
	end if 
	
	call ncd_io(ncid=ncid, varname='alb_svf', flag='read', data=nc_alb_svf, &
              dim1name=grlnd, nt=k, readvar=readvar)
    if ( .NOT. readvar) then
		write(iulog,*)subname, 'MML tried to read alb_svf, failed ', readvar
	end if 
	
	call ncd_io(ncid=ncid, varname='alb_gnf', flag='read', data=nc_alb_gnf, &
              dim1name=grlnd, nt=k, readvar=readvar)
    if ( .NOT. readvar) then
		write(iulog,*)subname, 'MML tried to read alb_gnf, failed ', readvar
	end if 
	
	call ncd_io(ncid=ncid, varname='alb_snf', flag='read', data=nc_alb_snf, &
              dim1name=grlnd, nt=k, readvar=readvar)
    if ( .NOT. readvar) then
		write(iulog,*)subname, 'MML tried to read alb_snf, failed ', readvar
	end if 
	
	! The rest:
	
    call ncd_io(ncid=ncid, varname='snowmask', flag='read', data=nc_snowmask, &
             dim1name=grlnd, nt=k)
    if ( .NOT. readvar) then
		write(iulog,*)subname, 'MML tried to read snowmask, failed ', readvar
	end if 
	
	call ncd_io(ncid=ncid, varname='bucketdepth', flag='read', data=nc_bucket, &
             dim1name=grlnd, nt=k)
    if ( .NOT. readvar) then
		write(iulog,*)subname, 'MML tried to read bucketdepth, failed ', readvar
	end if 
	
!	! I don't need these if I prescribe them directly below (if loading them causes problems)
!	call ncd_io(ncid=ncid, varname='max_soil_ice', flag='read', data=nc_ice, &
!              nt=k)
!    if ( .NOT. readvar) then
!		write(iulog,*)subname, 'MML tried to read max_soil_ice, failed ', readvar
!	end if 
!	
!	call ncd_io(ncid=ncid, varname='soilz', flag='read', data=nc_z, &
!              nt=k)
!    if ( .NOT. readvar) then
!		write(iulog,*)subname, 'MML tried to read soilz, failed ', readvar
!	end if 

    call ncd_io(ncid=ncid, varname='roughness', flag='read', data=nc_rough, &
             dim1name=grlnd, nt=k)
    if ( .NOT. readvar) then
		write(iulog,*)subname, 'MML tried to read roughness, failed ', readvar
	end if 
	
	call ncd_io(ncid=ncid, varname='evap_res', flag='read', data=nc_evaprs, &
             dim1name=grlnd, nt=k)
    if ( .NOT. readvar) then
		write(iulog,*)subname, 'MML tried to read evap_res, failed ', readvar
	end if 
	
	call ncd_io(ncid=ncid, varname='soil_type', flag='read', data=nc_type, &
             dim1name=grlnd, nt=k)
    if ( .NOT. readvar) then
		write(iulog,*)subname, 'MML tried to read soil_type, failed ', readvar
	end if 
	
	call ncd_io(ncid=ncid, varname='glc_mask', flag='read', data=nc_glc_mask, &
             dim1name=grlnd, nt=k)
    if ( .NOT. readvar) then
		write(iulog,*)subname, 'MML tried to read glc_mask, failed ', readvar
	end if 
	
	call ncd_io(ncid=ncid, varname='emissivity', flag='read', data=nc_emiss, &
             dim1name=grlnd, nt=k)
    if ( .NOT. readvar) then
		write(iulog,*)subname, 'MML tried to read emissivity, failed ', readvar
	end if 
	
	call ncd_io(ncid=ncid, varname='soil_tk_1d', flag='read', data=nc_soil_tk, &
             dim1name=grlnd, nt=k)
    if ( .NOT. readvar) then
		write(iulog,*)subname, 'MML tried to read soil_tk_1d, failed ', readvar
	end if 
	
	call ncd_io(ncid=ncid, varname='soil_cv_1d', flag='read', data=nc_soil_cv, &
             dim1name=grlnd, nt=k)
    if ( .NOT. readvar) then
		write(iulog,*)subname, 'MML tried to read soil_cv_1d, failed ', readvar
	end if 
	
	call ncd_io(ncid=ncid, varname='glc_tk_1d', flag='read', data=nc_glc_tk, &
             dim1name=grlnd, nt=k)
    if ( .NOT. readvar) then
		write(iulog,*)subname, 'MML tried to read glc_tk_1d, failed ', readvar
	end if 
	
	call ncd_io(ncid=ncid, varname='glc_cv_1d', flag='read', data=nc_glc_cv, &
             dim1name=grlnd, nt=k)
    if ( .NOT. readvar) then
		write(iulog,*)subname, 'MML tried to read glc_cv_1d, failed ', readvar
	end if 
	
!	call ncd_io(ncid=ncid, varname='glc_tk_2d', flag='read', data=nc_glc_tk, &
!             dim1name=grlnd, nt=k)
!    if ( .NOT. readvar) then
!		write(iulog,*)subname, 'MML tried to read glc_tk_2d, failed ', readvar
!	end if 
!	
!	call ncd_io(ncid=ncid, varname='glc_cv_2d', flag='read', data=nc_glc_cv, &
!             dim1name=grlnd, nt=k)
!    if ( .NOT. readvar) then
!		write(iulog,*)subname, 'MML tried to read glc_cv_2d, failed ', readvar
!	end if 
	
!	! MML:  having trouble reading in 3d data from netcdf file. Hard code them here until 
!	! 		I figure out how to fix it... 
!	write(iulog,*)subname, 'MML nc_import: manually implement column vars '
!	do i = 1, mml_nsoi
!		
!		nc_z(:,i) =  -0.025 * (exp(0.5*(i-0.5)) - 1.0);
!		
!		if (i == 1) then
!			nc_ice(begg:endg,i) = 0
!		else
!			nc_ice(begg:endg,i) = 300
!		end if
!		
!	enddo
	
	

	
	
    !write(iulog,*)subname, 'MML start pointing nc data at real variables '
    
    
    ! initially didn't have begg:endg; adding to see if it stops the crazy land-moving
    ! stuff happening in the h0 files (and prob the other .nc files, I haven't checked)
!    glc_mask(begg:endg) = nc_glc_mask(begg:endg)
 !   emiss(begg:endg) = nc_emiss(begg:endg)

    ! Albedos (direct)
	albedo_gvd(begg:endg) 	= nc_alb_gvd(begg:endg) 
	albedo_svd(begg:endg) 	= nc_alb_svd(begg:endg)
	albedo_gnd(begg:endg) 	= nc_alb_gnd(begg:endg) 
	albedo_snd(begg:endg) 	= nc_alb_snd(begg:endg)
	! Albedos (diffuse)
	albedo_gvf(begg:endg) 	= nc_alb_gvf(begg:endg) 
	albedo_svf(begg:endg) 	= nc_alb_svf(begg:endg)
	albedo_gnf(begg:endg) 	= nc_alb_gnf(begg:endg) 
	albedo_snf(begg:endg) 	= nc_alb_snf(begg:endg)
    ! Rest
	snowmask(begg:endg)		=	nc_snowmask(begg:endg)
	evaprs(begg:endg)		=	nc_evaprs(begg:endg)
	bucket_cap(begg:endg)	=	nc_bucket(begg:endg)
	!soil_maxice(begg:endg,:)	=	nc_ice(begg:endg,:)
	!soil_z(begg:endg,:)		=	nc_z(begg:endg,:)  ! 1:mml_nsoi ? 
	soil_type(begg:endg)	=	nc_type(begg:endg)
	roughness(begg:endg)	=	nc_rough(begg:endg)
	emiss(begg:endg)		= 	nc_emiss(begg:endg)
	glc_mask(begg:endg)		= 	nc_glc_mask(begg:endg)
	soil_tk_1d(begg:endg)	=	nc_soil_tk
	soil_cv_1d(begg:endg)	=	nc_soil_cv
	glc_tk_1d(begg:endg)	=	nc_glc_tk
	glc_cv_1d(begg:endg)	=	nc_glc_cv
	!dust(begg:endg)			= 	nc_dust(begg:endg)

!	
    !write(iulog,*)subname, 'MML done pointing nc data at real variables '
    
   ! write(iulog,*)subname, 'MML nc_alb_bare shape is', SHAPE(nc_alb_bare)
   !      write(iulog,*)subname, 'MML nc_alb_bare(begg,day) = ', nc_alb_bare(begg)
   !      write(iulog,*)subname, 'MML day = ', day
         !write(iulog,*)subname, 'albedo_bare(1) = ', albedo_bare(1)
    
    call ncd_pio_closefile(ncid) ! use ncid0 in case ncid got modified by one of the imports...
   	!write(iulog,*)subname, 'MML nc_import actually closed netcdf, yay '
   	


	deallocate(&
            nc_alb_gvd, 	&
            nc_alb_svd, 	&
            nc_alb_gnd, 	&
            nc_alb_snd, 	&
            nc_alb_gvf, 	&
            nc_alb_svf, 	&
            nc_alb_gnf, 	&
            nc_alb_snf, 	&
            nc_snowmask, 	&
            nc_evaprs,   	&
            nc_bucket,   	&
            nc_ice, 		&
            nc_z,   		&
            nc_type,     	&
            nc_rough, 		&
            nc_soil_tk, 	&
    		nc_glc_tk, 		&
    		nc_soil_cv, 	&
   			nc_glc_cv, 		&
    		nc_glc_mask, 	&
    		nc_emiss 	,	&
    		nc_dust			&
		)


    !write(iulog,*)subname, 'MML deallocated, hurray! nc_import subroutine end '
    
 
  end subroutine nc_import
  
  
  
  !***********************************************
  !	phase change
  !***********************************************
  subroutine phase_change (begg, endg, tsoi, soil_cv, soil_dz, &
  							soil_maxice, soil_liq, soil_ice, &
  							mml_nsoi, dt, hfus, tfrz, epc &
  							!diag1_1d, diag1_2d, diag2_2d, diag3_2d			& ! temporary diagnostics
  							)
	!% -------------------------------------------------------------------------
	! Given the initial soil temperature calculation, go check if we should be 
	! freezing/thawing any of the available freezeable water in that layer. 
	! Adjust temperature accordingly, if so.
	!% -------------------------------------------------------------------------
	
	implicit none
	
	! ----- Input Variables --------
  	integer	, intent(in)	:: begg, endg	! spatial bounds
  	integer	, intent(in)	:: mml_nsoi
  	real(r8), intent(in)	:: dt
  	real(r8), intent(in)	:: hfus
  	real(r8), intent(in)	:: tfrz
  	real(r8), intent(in)	:: soil_cv(:,:)
  	real(r8), intent(in)	:: soil_dz(:,:)
  	real(r8), intent(in)	:: soil_maxice(:,:)	! Not using this right now, instead using presc. soil_liq and soil_ice vals
  	
  	! ----- Output Variables --------
  	real(r8), intent(inout)		:: tsoi(:,:)	
  		! tsoi(begg:endg,:)	! try defining them this way instead, to avoid the dummy vars and keep the correct g indices
  	real(r8), intent(inout) 	:: soil_liq(:,:)	!
  	real(r8), intent(inout) 	:: soil_ice(:,:)	! 
  	
  	real(r8), intent(out) 	:: epc(:)	! (:,:) derivative of sat vapour pressure at ta [Pa/K]
  	
  !	real(r8), intent(out)	:: diag1_1d(:)		! put alhf here
  !	real(r8), intent(out)	:: diag1_2d(:,:)	! put rfm here
  !	real(r8), intent(out)	:: diag2_2d(:,:)	! put hfm here
  !	real(r8), intent(out)	:: diag3_2d(:,:)	! put hfmx here
  	
  	
  	! ----- Local Variables --------
	integer	:: i, g	! indexing variable
	! alternatley, could loop over numg = endg-begg + 1 
	! but you should be able to grab "numg" from another routine
	real(r8) :: ival
	character(len=32) :: subname = 'phase_change'
	
	real(r8)	::	alhf(begg:endg)
	
	real(r8)	::	wliq0(begg:endg,mml_nsoi)	,	&
					wice0(begg:endg,mml_nsoi)	,	&
					wmass0(begg:endg,mml_nsoi)	,	&
					tsoi0(begg:endg,mml_nsoi)	,	&
					imelt(begg:endg,mml_nsoi)	,	&
					hfm(begg:endg,mml_nsoi)	,	&
					hfmx(begg:endg,mml_nsoi)	,	&
					rfm(begg:endg,mml_nsoi)	,	&
					phase_cv(begg:endg,mml_nsoi)	,	&
					phase_dz(begg:endg,mml_nsoi)	,	&
					phase_maxice(begg:endg,mml_nsoi)	,	&
					phase_tsoi(begg:endg,mml_nsoi)	,	&
					phase_liq(begg:endg,mml_nsoi)	,	&
					phase_ice(begg:endg,mml_nsoi)	
	

	
	phase_cv = soil_cv
	phase_dz = soil_dz
	phase_maxice = soil_maxice
	phase_tsoi = tsoi
	phase_liq = soil_liq
	phase_ice = soil_ice
	
	!------------------------------------------------------
	
	!-----------------------------
	! Initialization
	wliq0 = phase_liq		! [kg/m2] per layer
	wice0 = phase_ice
	wmass0 = wliq0 + wice0		! should equal 300/dz in all but top layer, where it should be 0
	tsoi0 = phase_tsoi
	
	!-----------------------------
	! Identify if layers should be melting or freezing
		!	imelt = 0 -> no phase change
		!	imelt = 1 -> melt
		! 	imelt = 2 -> freeze
	imelt(:,:) = 0._r8
	
	!do i = 1, mml_nsoi
	do i = 1, mml_nsoi	! should be no freezeable water in top layer... ie phase_ice and phase_liq should both ==0
		
		! Melting: if there is ice and phase_tsoi > 0
		where (phase_ice(:,i) > 0._r8 .and. phase_tsoi(:,i) > tfrz)
			imelt(:,i) = 1
			phase_tsoi(:,i) = tfrz
		end where
		
		! Freezing: if there is water and phase_tsoi < 0
		where (phase_liq(:,i) > 0._r8 .and. phase_tsoi(:,i) < tfrz)
			imelt(:,i) = 2
			phase_tsoi(:,i) = tfrz
		end where
		
		! otherwise, leave phase_tsoi as is and don't put energy into phase change
		
	end do
	
	!-----------------------------
	! Energy available for freezing or melting comes from difference between phase_tsoi(:,i) and
	! tfreeze
	! 
	
	do i = 1, mml_nsoi
		
		hfm(:,i) = 0._r8	! all the palces imelt=0, no phase change
		
		! Energy for freezing or melting [W/m2]; hfm > 0 freezing, hfm < 0 melting
		where (imelt(:,i) > 0)
			hfm(:,i) = ( phase_tsoi(:,i) - tsoi0(:,i) ) * phase_cv(:,i) * phase_dz(:,i) / dt
			! if I accounted for cv water/ice here, too, would that fix part of the problem?
			
			! how much energy for freezing or melting based only off Delta T (if you've got excess, use for T change)
			! maybe I need to include water in cv to conserve energy? Hmm. Don't think gfdl does, though...
		end where
		
		! Melting: maximum energy available for freezing or melting [W/m2]
		where (imelt(:,i) .eq. 1)	! Melting case
			hfmx(:,i) = - phase_ice(:,i) * hfus / dt	! total meltable = depends how much ice you've got
		end where
		
		! Freezing: maximum energy available for freezing or melting [W/m2]
		where (imelt(:,i) .eq. 2)	! freezing case
			hfmx(:,i) = phase_liq(:,i) * hfus / dt	! total freezable = depends how much water you've got
		end where
		
	end do
  	
  	
  	!-----------------------------
	! Calculate phase change
	
	epc(:) = 0._r8
	
	do i = 1, mml_nsoi
		
		where( imelt(:,i) > 0 )
		
			! Freeze or melt ice
			rfm(:,i) = hfm(:,i) / hfus								! change in ice (>0 freeze, <0 melt) [kg/m2/s]
			phase_ice(:,i) = wice0(:,i) + rfm(:,i) * dt				! update ice [kg/m2]
			phase_ice(:,i) = max( 0.0 , phase_ice(:,i) )				! can't melt more ice than is present
			phase_ice(:,i) = min( wmass0(:,i) , phase_ice(:,i) )		! can't exceed total water than is present (300*dz, should be)
			phase_liq(:,i) = max( 0.0 , ( wmass0(:,i) - phase_ice(:,i) ) )	! update liquid water (kg/m2)
			alhf(:) = hfus * (phase_ice(:,i) - wice0(:,i)) / dt		! actual heat flux from phase change [w/m2]
			epc(:) = epc + alhf										! sum of heat flux from phase change over soil column [w/m2]
			
			! If there is energy left over, use it to change soil layer temperature
			phase_tsoi(:,i) = phase_tsoi(:,i) - (hfm(:,i) - alhf(:)) * dt / (phase_cv(:,i) * phase_dz(:,i))
		
		end where
	
	
	!---------------------
	! Error Checking:
		
		! Do error checking outside of where loop (can't seem to put if inside where)
		! MML: I have a feeling I'm doing this looping wrong, and maybe thats why the errors
		! are getting triggered... ie, I'm not sure this subroutine knows the indices
		! of the matrices that got passed in are begg:endg vs 1:(size(begg:endg)), so maybe
		! we're updating the wrong g indices? 
		do g = begg, endg
		
			! only check at places were we tried to change something
			if (imelt(g,i) > 0 ) then
			
				! Error check: make sure actual phase change isn't bigger than maximum energy available for change
				if ( abs(alhf(g)) > abs(hfmx(g,i)) ) then
					write(iulog,*)subname, 'MML ERROR: Soil temperature energy conservation error: phase change'
				end if
					
				! Freezing: make sure actual phase change does not exceed permissible phase change
      			! and that the change in ice does not exceed permissible change
				if (imelt(g,i) .eq. 2) then
					
					! Case 1: actual heat used for freezing differs from that possible by delta T and the max possible given amount of liquid water
					if ( abs( alhf(g) - min( hfm(g,i), hfmx(g,i) ) ) > 1.0e-03 ) then
						write(iulog,*)subname, 'MML ERROR: Soil temperature energy conservation error: freezing 01 ... ', &
						'alhf = ', alhf(g),', hfm(g,i) = ',hfm(g,i),', hfmx(g,i) = ',hfmx(g,i),&
						'abs(alhf - min(hfm,hfmx)) = ', abs( alhf(g) - min( hfm(g,i), hfmx(g,i) ) ), &
						' layer i = ', i
					end if
					
					if ( abs( ( phase_ice(g,i) - wice0(g,i) ) - ( min( hfm(g,i), hfmx(g,i) ) / hfus * dt) ) > 1.0e-03 ) then
						write(iulog,*)subname, 'MML ERROR: Soil temperature energy conservation error: freezing 02 ... ', &
						'check if < 1.0e-03: actually = ', &
						abs( ( phase_ice(g,i) - wice0(g,i) ) - ( min( hfm(g,i), hfmx(g,i) ) / hfus * dt) ), &
						' layer i = ', i, ' , phase_ice(g,i) - wice0(g,i) = ', phase_ice(g,i) - wice0(g,i) ,&
						' , min( hfm(g,i), hfmx(g,i) ) / hfus * dt) = ', min( hfm(g,i), hfmx(g,i) ) / hfus * dt
					end if
					
				end if
				
				
				! Thawing: make sure actual phase change does not exceed permissible phase change
      			! and that the change in ice does not exceed permissible change
				if (imelt(g,i) .eq. 1) then
					
					if ( abs(alhf(g) - max(hfm(g,i),hfmx(g,i))) > 1.0e-03 ) then
						write(iulog,*)subname, 'MML ERROR: Soil temperature energy conservation error: thawing 01 ... ', &
						'alhf = ', alhf(g),', hfm(g,i) = ',hfm(g,i),', hfmx(g,i) = ',hfmx(g,i), &
						'abs(alhf - max(hfm,hfmx)) = ', abs( alhf(g) - max( hfm(g,i), hfmx(g,i) ) ), &
						' layer i = ', i
					end if
					
					if ( abs( (phase_ice(g,i)-wice0(g,i)) - (max(hfm(g,i),hfmx(g,i)) / hfus * dt) ) > 1.0e-03 ) then
						write(iulog,*)subname, 'MML ERROR: Soil temperature energy conservation error: thawing 02 ... ', &
						'check if < 1.0e-03: actually = ', &
						abs( ( phase_ice(g,i) - wice0(g,i) ) - ( max( hfm(g,i), hfmx(g,i) ) / hfus * dt)) , &
						' layer i = ', i, ' , phase_ice(g,i) - wice0(g,i) = ', phase_ice(g,i) - wice0(g,i) ,&
						' , max( hfm(g,i), hfmx(g,i) ) / hfus * dt) = ', max( hfm(g,i), hfmx(g,i) ) / hfus * dt
					end if
					
				end if
			
			end if
				
		end do
	!------------------------
	
	end do
  	
  	! update out vars
  	soil_liq = phase_liq
  	soil_ice = phase_ice
  	tsoi	 = phase_tsoi
  	



  end subroutine phase_change
  
  
  !***********************************************
  !	soil_thermal_properties
  !***********************************************
  subroutine soil_thermal_properties (	begg, endg, glc_mask, &
  										soil_type, soil_z, soil_zh, soil_dz, &
  										soil_liq, soil_ice, &
  										soil_tk_1d, soil_cv_1d, &
  										glc_tk_1d, glc_cv_1d, &
 	 							   		soil_tk, soil_cv, soil_tkh, &
  								   		mml_nsoi)
	!% -------------------------------------------------------------------------
	!% From Gordon's hybrid.m matlab script: 
	!
	!% Calculate soil thermal conductivity and heat capacity
	! % (given the soil type and how much of the prescribed freezable water is 
	! frozen in that layer)
	!% -------------------------------------------------------------------------

  	implicit none
  	
	! ----- Input Variables --------
  	real(r8), intent(in)	:: soil_type(:)	! silt/sand/clay identified from a table (in theory... not yet :p )
  	real(r8), intent(in)	:: soil_z(:,:)	! soil depth (mid point of soil layer)
  	real(r8), intent(in)	:: soil_zh(:,:)	! soil depth (bottom interface of soil layer)
  	real(r8), intent(in)	:: soil_dz(:,:)	! soil layer thickness
  	real(r8), intent(in)	:: soil_liq(:,:)	! soil layer water content (kg/m2)
  	real(r8), intent(in)	:: soil_ice(:,:)	! soil layer ice content (kg/m2)
  	
  	real(r8), intent(inout)	:: soil_tk_1d(:)	! nc prescribed soil tk (for all layers)
  	real(r8), intent(inout)	:: soil_cv_1d(:)	! nc prescribed soil cv (for all layers)
  	real(r8), intent(inout)	:: glc_tk_1d(:)	! nc prescribed soil tk (for all layers)
  	real(r8), intent(inout)	:: glc_cv_1d(:)	! nc prescribed soil cv (for all layers)
  	
  	
  	integer	, intent(in)	:: mml_nsoi
  	integer	, intent(in)	:: begg, endg	! spatial bounds
  	real(r8), intent(in)	:: glc_mask(:)	! mask of glaciated cells, use ice properties here.
  	  	
  	! ----- Output Variables --------
  	real(r8), intent(out)	:: soil_tk(:,:)	! soil thermal resistance at each layer
  	real(r8), intent(out) 	:: soil_cv(:,:)	! soil heat capacity at each layer
  	real(r8), intent(out) 	:: soil_tkh(:,:)! soil thermal resistance at the boundary (bottom) of each layer
  	
  	! ----- Local Variables --------
	integer	:: i	! indexing variable

	real(r8)	::	watliq(begg:endg,mml_nsoi), watice(begg:endg,mml_nsoi)
	
	real(r8)	:: ival, cv_wat, cv_ice, rho_wat, rho_ice, cwat, cice, tk_ice
	
	
	
	! local constants:
	ival = 1e22
	rho_wat = 1000.0_r8	! density of water [kg/m3]
	rho_ice = 917.0_r8	! density of ice [kg/m3]
	cwat = 4188.0_r8		! specific heat of water [J/kg/K]
	cice = 2117.27_r8		! specific heat of ice [J/kg/K]
	!cv_wat = cwat * rho_wat	! [J/m3/K]	! = specific heat of water [J/kg/K] * density of water [kg/m3]
	!		= 4.188e6 J/kg/K
	!cv_ice = cice * rho_ice	! [J/m3/K]  ! = specific heat of ice [J/kg/K] * density of ice [kg/m3]
	!		= 1.942e6 J/kg/K
	cv_wat = 4.188e6_r8		! [J/m3/K]	! from LSM
	cv_ice = 2.094e6_r8		! [J/m3/K]	! from LSM
	tk_ice = 2.2			! [W/m/K]	! Near 0 C ; taking from LSM ... down to 2.22 at 0, up to 3.48 at -100 C
  	
  	! this gives a cv_ice of 1941536.59 = 1.942e6 J/kg/K
  	! or cv_ice = 2.094 x 10^6 J/kg/K from LSM
  	! and the calc
  	! Later: add dependence on soil type. Now I'm setting all soils to have the same tk and cv
	! (consider using the table-implementation in Gordon's code and in the GFDL code)
  	
  	! calculate the volumetric liquid / ice water content in each soil layer:
  	watliq = soil_liq / (rho_wat * soil_dz)		! [kg/m2] / ([kg/m3] * [m]) -> unitless? hmm... or m3/m3 I guess
  	watice = soil_ice / (rho_ice * soil_dz)		! m3/m3 ?
  	
  	! I'm assuming matrix addition works as I expect in Fortran?
	! ie I don't have to loop over g = begg,endg, do I? (I might if it goes spatially 
	! varying and the equations aren't the same have to check. But for now, implement like this)
	
	do i = 1, mml_nsoi
		
		! For soil points (non-glaciated), use these values:
		
		!soil_tk(:,i) = 1.5_r8 		! [W/m/K]	! in the ballpark of that for various soils in LaD
		soil_tk(:,i) = soil_tk_1d(:)
		
		!soil_cv(:,i) = 2.0e06_r8	! [J/m3/K]	! that used for "medium" soil in LaD
		soil_cv(:,i) = soil_cv_1d(:)
		
		! If the point is a glacier (glc_mask=1), use these values instead:
		where(glc_mask .eq. 1)	! really, I should make glc_mask a logical...
		
			! Using heat capacity and thermal resistance of ice
			!soil_tk(:,i) = tk_ice	! [W/m/K]	
			soil_tk(:,i) = glc_tk_1d(:)
			
			!2.3_r8 		! [W/m/K]
			! value somewhat arbitrarily taken from:
			! http://www.engineeringtoolbox.com/ice-thermal-properties-d_576.html
			! ... find a more supportable value to use in the end
			
			!soil_cv(:,i) = cv_ice	! [J/m3/K]
			soil_cv(:,i) =	glc_cv_1d(:)
			
			!1.8e06_r8	! [J/m3/K]
			! value somewhat arbitrarily taken from:
			! http://www.engineeringtoolbox.com/ice-thermal-properties-d_576.html
			! ... find a more supportable value to use in the end
			
		end where
		
		
		! later: add water to thermal resistance? or no? 
		! is this right? soil_cv = actual_soil_cv + water_cv + ice_cv ?
		!soil_cv(:,i) = 1.926e06 + cv_wat*watliq(:,i) + cv_ice*watice(:,i)	! [J/m3/K]
		
	enddo ! loop over all soil layers and assign them the 1d value
	
	
	soil_tkh(:,:) = 0.0_r8 	! for now, just so each entry has a value (it should really be size (:, mml_soi-1), not (:,mml_nsoi)
	! now find tkh
	do i = 1, mml_nsoi-1 ! no heat diffusion through bottom layer
		soil_tkh(:,i) = soil_tk(:,i) * soil_tk(:,i+1) * ( soil_z(:,i) - soil_z(:,i+1) ) / &
						( soil_tk(:,i) * (soil_zh(:,i) - soil_z(:,i+1)) + &
						  soil_tk(:,i+1) * (soil_z(:,i) - soil_zh(:,i)) )
						  
						  ! This LOOKS the same as the matlab eq'n... add zh to 
						  ! output and see if that looks okay...
						  
		! NOTE: tk and tkh not currently dependent on water/ice content of layer!
		! ... but I'll keep it like that, for now anyhow. More straightforward. 
	enddo
	

  	
   
  end subroutine soil_thermal_properties
  
  
  
  !***********************************************
  !	mo_hybrid: solve for obukhov root
  !***********************************************
  
  ! OR call most and mo_hybrid FROM within a g loop, so both of these also just act on single points... 
  ! lets do that for now. Sigh. 
  subroutine mo_hybrid ( &
  						! variables for most fn
  						!begg, endg, 
  						 zref, uref, thref, qref, disp, z0m, z0h, tsrf, qsrf, vkc, grav, &
  						ustar, tstar, tvstar, qstar, obu_new, & ! fx, &
  						! specific hybrid vars
  						obu_root, obu0, obu1, tol &
  					)

  	implicit none
  	
  	! Input variables
  	!integer, intent(in)		:: begg, endg
  	
  	!real(r8), intent(in)	:: obu_old(:), zref(:), uref(:), thref(:), disp(:), &
  	!						   z0m(:), z0h(:), tsrf(:)
  	real(r8), intent(in)	:: zref, uref, thref, qref, disp, &
  							   z0m, z0h, tsrf, qsrf, vkc, grav, &
  							   obu0, obu1, tol
  							   						   
  	
  	! Output variables
  	real(r8), intent(out)	:: obu_root, ustar, tstar, tvstar, qstar, obu_new

  	
  	! Local variables
  	!real(r8), allocatable, dimension(:)	:: root, x0, x1, minx, minf, f0, f1, dx, x
  	real(r8)	:: root, x0, x1, minx, minf, f0, f1, dx, x, ival 	
  	! (f0 and f1 will be past back from most as fx)
 
  	integer		:: i, itmax, g, tag
  	
  	character(len=32) :: subname = 'MML_mo_hybrid'
  	
  	ival = 0._r8

  	
  	!---------------------------------------------------------------------------------
  	
  	! temporary out vals
  	obu_root = ival;
  	ustar = ival;
  	tstar = ival;
  	tvstar = ival;
  	qstar = ival;
  	obu_new = ival;
  	
  	!---------------------------------------------------------------------------------
  	! Already assuming this is being called from a loop in g, so everything is just point value
  	
  		i = 0 	! value for i if we don't hit do loop
  		!----------------------------------------------------------
	 	! (1) Evaluate func at xa and see if this is the root
  		
  		x0 = obu0
  		call most(  x0, zref, uref, thref, qref, disp, z0m, z0h, tsrf, qsrf, vkc, grav, &
  					ustar, tstar, tvstar, qstar, obu_new, f0)
  		if (f0 == 0.0) then
  			obu_root = x0
  			tag = 0
  			!write(iulog,*)subname, 'MML hybrid exit (hopefully) 0'
  			!exit	! I want to exit the whole subroutine, not just the if statment
  		end if
  		
  		! for now, instead of exit statments, just wrap it in if/else statments
  		
  		if (f0 .ne. 0.0_r8) then
  		
  			!----------------------------------------------------------
	  		! (2) Evaluate func at xb and see if this is the root
			!
			! (only want to do this if above failed...)
			
			x1 = obu1
			call most(  x1, zref, uref, thref, qref, disp, z0m, z0h, tsrf, qsrf, vkc, grav, &
  						ustar, tstar, tvstar, qstar, obu_new, f1)
			if (f1 == 0.0_r8) then
				obu_root = x1
				tag = 1
				!write(iulog,*)subname, 'MML hybrid exit (hopefully) 1'
				!exit
			end if
			
			if (f1 .ne. 0.0_r8) then
			
				!----------------------------------------------------------
	    		! (3) Order initial root estimates correctly
				!    
    			! if neither of the above holds, do something more fancy
				
				if (f1 < f0) then
					minx = x1
					minf = f1
				else
					minx = x0
					minf = f0
				end if
		
				!----------------------------------------------------------
	    		! (4) Iterative root calculation. Use the secant method, and Brent's method as a backup
				
				itmax = 40
				do i = 1, itmax
					
					dx = -f1 * (x1 - x0) / (f1 - f0)
					x = x1 + dx
					
					! Check if x is the root. If so, exit the iteration
					if( abs(dx) < tol) then
						x0 = x
						tag = 2
						!write(iulog,*)subname, 'MML hybrid exit 2'
						exit		! I think this will actually exit the do loop... not sure
					end if	
					
					! Otherwise keep going, and evaluate the fn at x
					x0 = x1
					f0 = f1
					x1 = x
					call most(  x1, zref, uref, thref, qref, disp, z0m, z0h, tsrf, qsrf, vkc, grav, &
  						ustar, tstar, tvstar, qstar, obu_new, f1)
					
					if( f1 < minf) then
						minx = x1
						minf = f1
					end if
					
					! if a root zone is found, use Brent's method for a robust backup strategy:
					
					if( f1*f0 < 0.0_r8) then
						! call zbrent(x, x0, x1)	!x0, x1 in, x out, plus a new obu length, non?
						!... needs all most vars; xa, xb, and tol
						call zbrent( &
								 ! most vars in: (no obu_old in, we'll use x0 and x1 to guess them)
								 zref, uref, thref, qref, disp, z0m, z0h, tsrf, qsrf, vkc, grav, &
								 ! most vars out:
  								ustar, tstar, tvstar, qstar, obu_new, &
  								 ! zbrent specific vars 
  								 x0, x1, tol, x)	! x = root from zbrent
						x0 = x
						tag = 3
						!write(iulog,*)subname, 'MML hybrid exit 3'
						exit
					end if
					
					! If we fail to converge within itmax iterations, stop with the minimum value fn
					
					if( i == itmax ) then
						call most(  minx, zref, uref, thref, qref, disp, z0m, z0h, tsrf, qsrf, vkc, grav, &
  							ustar, tstar, tvstar, qstar, obu_new, f1)
  						x0 = minx	! what is the point of calling most on this? I guess it updates obu_new, which we'll pass out
						tag = 4
						write(iulog,*)subname, 'MML hybrid exit 4 (caution, itmax reached)'
						exit
					end if
					
					
					
					
				end do	! loop i to itmax
				
				!write(iulog,*)subname, 'MML obu hybrid: exited loop at i = ',i,', tag = ', tag
				
				
				obu_root = x0
				
				
				
			
			end if 	! f1 ~= 0.0, root = x1
		
		end if	! f0 ~= 0.0 , root = x1
		
	
	!obu_root = root
	
  	!---------------------------------------------------------------------------------

    
  end subroutine mo_hybrid
 
 
  !***********************************************
  !	most: Monin-Obukhov Stability Function
  !***********************************************
  
  ! as 1-d point function
  subroutine most ( & !begg, endg, 
  					obu_old, zref, uref, thref, qref, disp, z0m, z0h, tsrf, qsrf, vkc, grav, &
  					ustar, tstar, tvstar, qstar, obu_new, fx)
	!% -------------------------------------------------------------------------
	!% From Gordon's hybrid.m matlab script: 
	!
	!% Use Monin-Obukhov similarity theory to obtain the Obukhov length (obu).
	!% This is the function to solve for the Obukhov length. For the current estimate
	!% of the Obukhov length (x), calculate u* and T* and then the new length (obu).
	!% The function value is the change in Obukhov length: fx = x - obu.
	!% -------------------------------------------------------------------------

  	implicit none
  	
  	! Game plan (2016.03.23):
  	! Loop it all over g (ugly, but might understand how to do it better than putting
  	! function calls inside where statements)
  	! Rewrite a zbrent subroutine following the matlab version such that we can update
  	! the obu length while minimizing x (the exiting fortran one can only work on functions
  	! with one output, which our current most subroutine isn't)
  	! Write the psi functions as actual functions, just for kicks (rather than subroutines).
  	! This will force them to only have one output, and they'll only be able to work 1d,
  	! but assumably there is a reason clm did it that way... (right, we can call them directly,
  	! rather than saying call psi_m(lalala) on its own line, we can put a*psi_m(lalala) right
  	! into a statment). 
  	! 
  	
  	! CALL THIS FROM WITHIN A g LOOP!!!
  	
  	
  	! ----- Input Variables --------
  	!integer	, intent(in)	:: begg, endg	! spatial bounds
  	 
  	real(r8), intent(in)	:: obu_old	! current estimate of obukhov length [m]
  	real(r8), intent(in)	:: zref			! reference height [m]
  	real(r8), intent(in)	:: uref			! wind speed at reference height [m/s]
	real(r8), intent(in)	:: thref		! potential temperature at reference height [m/s]
  	real(r8), intent(in)	:: qref			! specific humidity at reference height [kg/kg]
  	real(r8), intent(in)	:: disp			! displacement height [m]
  	real(r8), intent(in)	:: z0m			! roughness length for momentum [m]
  	real(r8), intent(in)	:: z0h			! roughness length for heat [m]
  	real(r8), intent(in)	:: tsrf			! surface temperature [K]
  	real(r8), intent(in)	:: qsrf			! surface humidity [kg/kg]
  	real(r8), intent(in)	:: vkc				! von Karman constant
  	real(r8), intent(in)	:: grav				! gravitational constant	! if I define vkc and grav at top of module, can I avoid passing them around?
  	 	
  	! ----- Output Variables --------
  	real(r8), intent(out)	:: ustar	! soil thermal resistance at each layer
  	real(r8), intent(out) 	:: tstar	! soil heat capacity at each layer
  	real(r8), intent(out) 	:: tvstar! soil thermal resistance at the boundary (bottom) of each layer
	real(r8), intent(out) 	:: qstar! soil thermal resistance at the boundary (bottom) of each layer
  	real(r8), intent(out) 	:: obu_new! soil thermal resistance at the boundary (bottom) of each layer
  	real(r8), intent(out) 	:: fx! soil thermal resistance at the boundary (bottom) of each layer
  	
  	
  	! ----- Local Variables --------
	integer	:: i, g	! indexing variable
	real(r8)	ival
	
	!real(r8), allocatable, dimension(:) 	::	psi_m_zref, psi_m_z0m, psi_h_zref, psi_h_z0h, &
	!										z_minus_d, zlog_psim, zlog_psih, obu_temp
											! psi_m, psi_h, &
	real(r8) 	:: psi_m_zref, psi_m_z0m, psi_h_zref, psi_h_z0h, &
				   z_minus_d, zlog_psim, zlog_psih, obu_temp, thv, &
				   zeta, zldis	
				   								
	character(len=32) :: subname = 'MML_most'
	
	ival = 0._r8

  	
  	obu_temp = obu_old
  	
  	! for now, just to compile
  	ustar = ival
  	tstar = ival
  	tvstar = ival
  	qstar = ival
  	obu_new = ival
  	fx = ival
  	
  	
  	!-----------------------------------------------------------------------
  	! ASSUME WE'RE ALREADY IN A g LOOP
  	
  	!-----------------------------
  	! Prevent near-zero values of obukhov length by imposing a minimum of 0.1 m
  	if ( abs(obu_temp) <= 0.1_r8) then
  		obu_temp = 0.1_r8
  	end if
  	
  	! check if zref < disp
  	if( zref - disp < 0.0_r8 ) then
  		write(iulog,*)subname, 'MML ERROR; zref < disp! zref = ',zref,' disp = ',disp
  	end if
  	
  	
  	!-----------------------------
  	! get psi_m and psi_h for zref and z0
  	
  	! Evaluate psi for momentum at ref height zref and surface (disp+z0m)
  	psi_m_zref 	= psi_m( (zref - disp) / obu_temp)	! what is zref < disp? ie zref = 2m, disp = 10m
  	! force zref to be above disp
  	!psi_m_zref 	= psi_m( (zref ) / obu_temp)	! what is zref < disp? ie zref = 2m, disp = 10m
  	psi_m_z0m	= psi_m( z0m / obu_temp )
  	
  	! Evaluate psi for heat at ref height zref and surface (disp+z0h)
  	psi_h_zref 	= psi_h( (zref - disp) / obu_temp)	! what is zref < disp? ie zref = 2m, disp = 10m
  	! force zref to be above disp
  	!psi_h_zref 	= psi_h( (zref ) / obu_temp)	! what is zref < disp? ie zref = 2m, disp = 10m
  	psi_h_z0h	= psi_h( z0h / obu_temp )
  	
  	!-----------------------------
  	! Calculate ustar [m/s] and tstar [K]
  	
  	z_minus_d 	= 	zref - disp	! what if disp > zref???? 
  	! temporarily just set to zref:
  	!z_minus_d = zref	
  	zlog_psim	=	log( z_minus_d / z0m ) - psi_m_zref + psi_m_z0m
  	zlog_psih	= 	log( z_minus_d / z0h ) - psi_h_zref + psi_h_z0h
  	
  	ustar	= 	uref * vkc / zlog_psim
  	tstar	=	(thref - tsrf) * vkc / zlog_psih
  	
  	!-----------------------------
  	! Calculate qstar using sfc and atm specific humidity
  	qstar 	=	(qref - qsrf) * vkc / zlog_psih
  	
  	!-----------------------------
  	! Calculate theta_v_star using new value of qstar 
  	tvstar 	=	tstar + 0.61_r8 * tsrf * qstar
  	
  	!-----------------------------
  	! Need theta_v for calculation (CHECK THIS EQUATION!!!)
  	thv	=	thref * (1 + 0.61_r8*qref)
  	
  	!-----------------------------
  	! Calculate obukhov length [m] 
  	obu_new = max( 0.1_r8, ustar**2 * thv / ( vkc * grav * tvstar ) )
  	
  	
  	! -------------------------------------------------------------
    ! limit the obukhov length by putting bounds on zeta (inside obu subroutine, probably)
    ! zeta = (z-d)/obu, zldis = z-d
    zldis = zref - disp
    zeta = (zref - disp) / obu_new
    
    if( zeta >= 0._r8) then 	! stable
    	zeta = min(2._r8,max(zeta,0.01_r8))
    else
    	zeta = max(-2._r8,min(zeta,-0.01_r8))	! or -100 instead of -2
    end if
    
    obu_new = zldis / zeta
  	
  	
  	!-----------------------------
  	! Calculate change in obukhov length [m] 
  	
  	fx = obu_temp - obu_new		! not obu_old - obu_new? potentially modify obu_temp if obu_old is too small, but this is what hte matlab code did

  	!-----------------------------------------------------------------------

  end subroutine most
  
 
  
   !=============================================================
  !*** Psi functions for most

  
  ! Monin-Obukhov psi function for momentum 
   real(r8) function psi_m (zeta)
	!% -------------------------------------------------------------------------
	!% From Gordon's psi_m_monin_obukhov.m
	!% -------------------------------------------------------------------------
	
	! written like this, it can probably only act on a single value, not on an
	! array, ie I need to loop it in g... pity...
  	implicit none
  	
  	
  	! ----- Input Variables --------
  	real(r8), intent(in)	:: zeta	! current estimate of obukhov length [m]
  	!integer, intent(in)		:: begg, endg
  	! ----- Output Variables --------
  	!real(r8), intent(out)	:: psi_m(:)	! 
  	
  	! ----- Local Variables --------
	!real(r8), allocatable, dimension(:)	:: y!, obu_temp
	!real(r8), allocatable, dimension(:)	:: psi_m!, obu_temp
	real(r8)	::	ival, pi
	real(r8)	y
	
	ival = 0._r8
	pi = 4._r8*atan(1._r8)		! this is defined with actual numbers in shr_const_mod
	

	!----------------------------
	
	if ( zeta < 0.0_r8 ) then
   		y = (1.0_r8 - 16.0_r8 * zeta)**0.25_r8
   		psi_m = 2.0_r8 * log((1.0_r8 + y)/2.0_r8) + log((1.0_r8 + y**2)/2.0_r8) - 2.0_r8*atan(y) + pi/2.0_r8
   	else
   		!if ( zeta >= 0.0 ) then !else
    	psi_m = -5.0_r8 * zeta
    end if

  end function psi_m
  
  
  ! Monin-Obukhov psi function for heat+water 
  real(r8) function psi_h (zeta)
	!% -------------------------------------------------------------------------
	!% From Gordon's psi_m_monin_obukhov.m
	!% -------------------------------------------------------------------------

  	implicit none
  	
  	
  	! ----- Input Variables --------
  	real(r8), intent(in)	:: zeta	! current estimate of obukhov length [m]
  	!integer, intent(in)		:: begg, endg
  	! ----- Output Variables --------
  	!real(r8), intent(out)	:: psi_h(:)	! 
  	
  	! ----- Local Variables --------
	!real(r8), allocatable, dimension(:)	:: y
	real(r8)	:: 	y
	

	!----------------------------
	if ( zeta < 0.0 ) then
   		y = (1.0_r8 - 16.0_r8 * zeta)**0.25_r8
   		psi_h = 2.0_r8 * log((1.0_r8 + y**2)/2.0_r8) 
   	else 
   		! if ( zeta >= 0.0 ) !else
    	psi_h = -5.0_r8 * zeta
    end if
   

   
  end function psi_h

  !*** End psi for most functions
  !=============================================================
  
  
!=============================================================
! zbrent solver
  subroutine zbrent( &
					! most vars in: (no obu_old in, we'll use x0 and x1 to guess them)
					zref, uref, thref, qref, disp, z0m, z0h, tsrf, qsrf, vkc, grav, &
					! most vars out:
  					ustar, tstar, tvstar, qstar, obu_new, &! fx, &
  					! zbrent specific vars 
  					xa, xb, tol, root)
  
  	implicit none 
  	
  	! ----- Input Variables --------
  	real(r8), intent(in)	:: zref	
  	real(r8), intent(in)	:: uref	
  	real(r8), intent(in)	:: thref	
  	real(r8), intent(in)	:: qref	
  	real(r8), intent(in)	:: disp	
  	real(r8), intent(in)	:: z0m	
  	real(r8), intent(in)	:: z0h	
  	real(r8), intent(in)	:: tsrf	
  	real(r8), intent(in)	:: qsrf	
  	real(r8), intent(in)	:: vkc	
  	real(r8), intent(in)	:: grav
  	
  	real(r8), intent(in)	:: xa, xb, tol
  		
  	! ----- Output Variables --------
  	real(r8), intent(out)	:: ustar
  	real(r8), intent(out)	:: tstar
  	real(r8), intent(out)	:: tvstar
	real(r8), intent(out)	:: qstar
	real(r8), intent(out)	:: obu_new
	!real(r8), intent(out)	:: fx
	real(r8), intent(out)	:: root
  	
  	! ----- Local Variables --------
	!real(r8), allocatable, dimension(:)	:: y
	real(r8)	:: 	a, b, c, d, e, &
					fa, fb, fc, &
					p, q, r, s, &
					xm, &
					eps1, tol1
	integer		:: 	itmax, iter
	real(r8)	::	ival, pi, i
	character(len=32) :: subname = 'MML_zbrent'
	
	ival = 0._r8

	
	!--------------------------------------------------------------
	
	a = xa
	b = xb
	call most ( a, zref, uref, thref, qref, disp, z0m, z0h, tsrf, qsrf, vkc, grav, &
  					ustar, tstar, tvstar, qstar, obu_new, fa)
  	call most ( b, zref, uref, thref, qref, disp, z0m, z0h, tsrf, qsrf, vkc, grav, &
  					ustar, tstar, tvstar, qstar, obu_new, fb)
	
	if( ( fa > 0.0_r8 .and. fb > 0.0_r8 ) .or. ( fa < 0.0_r8 .and. fb < 0.0_r8 ) ) then
		! error: 
		write(iulog,*)subname, 'MML brent error: root must be bracketed'
	end if
	
	itmax = 50		! max number of iterations
	eps1 = 1.0e-08_r8	! relative error tolerance
	
	c = b
	fc = fb
	
	do iter = 1, itmax
		
		if( (fb > 0.0_r8 .and. fc > 0.0_r8) .or. (fb < 0.0_r8 .and. fc < 0.0_r8) ) then
			c = a
			fc = fa
			d = b - a
			e = d
		end if
		
		if( abs(fc) < abs(fb) ) then
			a = b
			b = c
			c = a
			fa = fb
			fb = fc
			fc = fa
		end if
		
		tol1 = 2.0_r8 * eps1 * abs(b) + 0.5_r8 * tol
		xm = 0.5_r8 * (c - b)
		
		if( abs(xm) <= tol1 .or. fb .eq. 0.0_r8) then
			!write(iulog,*)subname, 'MML zbrent exit 1'
			exit	! stop iterating?
		end if
		
		if( abs(e) >= tol1 .and. abs(fa) > abs(fb) ) then
		
			s = fb/ fa
			
			if (a == c) then
				p = 2.0_r8 * xm * s
				q = 1.0_r8 - s
			else
				q = fa / fc
				r = fb / fc
				p = s * (2.0_r8 * xm * q * (q-r) - (b-a) * (r - 1.0_r8) )
				q = (q - 1.0_r8) * (r - 1.0_r8) * (s - 1.0_r8)
			end if
			
			if( p > 0.0_r8) then
				q = -q
			end if
			
			p = abs(p)
			
			if( 2.0_r8*p < min( 3.0_r8*xm*q-abs(tol1*q) , abs(e*q) ) ) then
				e = d
				d = p/q
			else
				d = xm
				e = d
			end if
			
		else
			
			d = xm
			e = d
			
		end if
		
		a = b
		fa = fb
		
		if( abs(d) > tol1 ) then
			b = b + d
		else
			if(	xm >= 0.0_r8 ) then
				b = b + abs(tol1)
			else
				b = b - abs(tol1)
			end if
		end if
		! call most
		call most ( b, zref, uref, thref, qref, disp, z0m, z0h, tsrf, qsrf, vkc, grav, &
  					ustar, tstar, tvstar, qstar, obu_new, fb)
		if( fb == 0.0_r8 ) then
			!write(iulog,*)subname, 'MML zbrent exit 2'
			exit;
		end if
		
		if( iter == itmax ) then
			!write(iulog,*)subname, 'MML brent ERROR: maximum number of iterations exceeded '
		end if 
	
	end do
	!write(iulog,*)subname, 'MML brent iteration exit with iter = ',iter
	
	root = b
	
	
	!most ( & !begg, endg, 
  !					obu_old, zref, uref, thref, qref, disp, z0m, z0h, tsrf, qsrf, vkc, grav, &
  !					ustar, tstar, tvstar, qstar, obu_new, fx)
	
	!--------------------------------------------------------------
  
  end subroutine zbrent
  

end module mml_mainMod
