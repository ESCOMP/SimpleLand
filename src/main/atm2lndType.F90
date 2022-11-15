module atm2lndType

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Handle atm2lnd, lnd2atm mapping
  !
  ! !USES:
  use shr_kind_mod  , only : r8 => shr_kind_r8
  use shr_infnan_mod, only : nan => shr_infnan_nan, assignment(=)
  use shr_log_mod   , only : errMsg => shr_log_errMsg
  use clm_varpar    , only : numrad, ndst, nlevgrnd !ndst = number of dust bins.  ! MML: numrad = 2, 1=vis, 2=nir
  use clm_varcon    , only : rair, grav, cpair, hfus, tfrz, spval
  use clm_varctl    , only : iulog
  use decompMod     , only : bounds_type
  use abortutils    , only : endrun
  use PatchType     , only : patch
  !
  ! !PUBLIC TYPES:
  implicit none
  private
  save
  !
  ! !PUBLIC DATA TYPES:

  type, public :: atm2lnd_params_type
     ! true => repartition rain/snow from atm based on temperature
     logical :: repartition_rain_snow

     ! Surface temperature lapse rate (K m-1)
     real(r8) :: lapse_rate

     ! Rain-snow ramp for glacier landunits
     ! frac_rain = (temp - all_snow_t) * frac_rain_slope
     ! (all_snow_t is in K)
     real(r8) :: precip_repartition_glc_all_snow_t
     real(r8) :: precip_repartition_glc_frac_rain_slope

     ! Rain-snow ramp for non-glacier landunits
     ! frac_rain = (temp - all_snow_t) * frac_rain_slope
     ! (all_snow_t is in K)
     real(r8) :: precip_repartition_nonglc_all_snow_t
     real(r8) :: precip_repartition_nonglc_frac_rain_slope
  end type atm2lnd_params_type

  !----------------------------------------------------
  ! atmosphere -> land variables structure
  !
  ! NOTE:
  ! IF there are forcing variables that are downscaled - then the
  ! non-downscaled versions SHOULD NOT be used in the code. Currently
  ! the non-downscaled versions are only used n a handful of places in
  ! the code (and needs to be used in lnd_import_export and the
  ! downscaling routines), but in general should NOT be used in new
  ! code. Instead use the datatype variables that have a _col suffix
  ! which gives the downscaled versions of these fields.
  !
  ! MML: I don't think this applies to me... I'm working at the grc level, not the col level...
  !----------------------------------------------------
  type, public :: atm2lnd_type
     type(atm2lnd_params_type) :: params

     ! atm->lnd not downscaled
     real(r8), pointer :: forc_u_grc                    (:)   => null() ! atm wind speed, east direction (m/s)
     real(r8), pointer :: forc_v_grc                    (:)   => null() ! atm wind speed, north direction (m/s)
     real(r8), pointer :: forc_wind_grc                 (:)   => null() ! atmospheric wind speed
     real(r8), pointer :: forc_hgt_grc                  (:)   => null() ! atmospheric reference height (m)
     real(r8), pointer :: forc_topo_grc                 (:)   => null() ! atmospheric surface height (m)
     real(r8), pointer :: forc_hgt_u_grc                (:)   => null() ! obs height of wind [m] (new)
     real(r8), pointer :: forc_hgt_t_grc                (:)   => null() ! obs height of temperature [m] (new)
     real(r8), pointer :: forc_hgt_q_grc                (:)   => null() ! obs height of humidity [m] (new)
     ! mml maybe use these
     real(r8), pointer :: forc_vp_grc                   (:)   => null() ! atmospheric vapor pressure (Pa)
     real(r8), pointer :: forc_rh_grc                   (:)   => null() ! atmospheric relative humidity (%)
     real(r8), pointer :: forc_psrf_grc                 (:)   => null() ! surface pressure (Pa)
     real(r8), pointer :: forc_pco2_grc                 (:)   => null() ! CO2 partial pressure (Pa)
     real(r8), pointer :: forc_pco2_240_patch           (:)   => null() ! 10-day mean CO2 partial pressure (Pa)
     real(r8), pointer :: forc_solad_grc                (:,:) => null() ! direct beam radiation (numrad) (vis=forc_sols , nir=forc_soll )
     real(r8), pointer :: forc_solai_grc                (:,:) => null() ! diffuse radiation (numrad) (vis=forc_solsd, nir=forc_solld)
     real(r8), pointer :: forc_solar_grc                (:)   => null() ! incident solar radiation
     real(r8), pointer :: forc_pc13o2_grc               (:)   => null() ! C13O2 partial pressure (Pa)
     real(r8), pointer :: forc_po2_grc                  (:)   => null() ! O2 partial pressure (Pa)
     real(r8), pointer :: forc_po2_240_patch            (:)   => null() ! 10-day mean O2 partial pressure (Pa)
     real(r8), pointer :: forc_aer_grc                  (:,:) => null() ! aerosol deposition array

     real(r8), pointer :: forc_t_not_downscaled_grc     (:)   => null() ! not downscaled atm temperature (Kelvin)       
     real(r8), pointer :: forc_th_not_downscaled_grc    (:)   => null() ! not downscaled atm potential temperature (Kelvin)    
     real(r8), pointer :: forc_q_not_downscaled_grc     (:)   => null() ! not downscaled atm specific humidity (kg/kg)  
     			! MML: I think this is the q I need to check if the negative LH is too big. 
     real(r8), pointer :: forc_pbot_not_downscaled_grc  (:)   => null() ! not downscaled atm pressure (Pa)   
     real(r8), pointer :: forc_pbot240_downscaled_patch (:)   => null() ! 10-day mean downscaled atm pressure (Pa)           
     real(r8), pointer :: forc_rho_not_downscaled_grc   (:)   => null() ! not downscaled atm density (kg/m**3)                      
     real(r8), pointer :: forc_rain_not_downscaled_grc  (:)   => null() ! not downscaled atm rain rate [mm/s]                       
     real(r8), pointer :: forc_snow_not_downscaled_grc  (:)   => null() ! not downscaled atm snow rate [mm/s]                       
     real(r8), pointer :: forc_lwrad_not_downscaled_grc (:)   => null() ! not downscaled atm downwrd IR longwave radiation (W/m**2) 
     
     ! MML: 2016.01.14 Adding a variable 2TBOT to see if I can print it to the h0 file
	 !real(r8), pointer :: forc_2t_not_downscaled_grc     (:)   => null() ! not downscaled atm temperature (Kelvin)       
	 
	 ! MML: 2016.01.14 Simple Land Energy and Hydrology variables (gridscale)
	 ! Instead of being passed x2l from coupler, these are going to be defined by the
	 ! simple land model. But they're being attached to the atm2lndType (for now at least)
!	 real(r8), pointer :: mml_swdn_grc     				(:)   => null() ! gridcell level incoming SW (from atm)
!	 real(r8), pointer :: mml_swup_grc     				(:)   => null() ! gridcell level reflected SW (calculated in simple model)
!	 real(r8), pointer :: mml_swab_grc     				(:)   => null() ! gridcell level absorbed SW (calculated in simple model)
!	 real(r8), pointer :: mml_lwup_grc     			(:)   => null() ! gridcell level emitted LW (calculated in simple model)
!	 real(r8), pointer :: mml_lwdn_grc    				(:)   => null() ! gridcell level incoming LW (from atm)
!	 real(r8), pointer :: mml_ts_grc     				(:)   => null() ! gridcell level surface temperature (calculated in simple model)
!	 real(r8), pointer :: mml_h2o_grc     				(:)   => null() ! gridcell level soil bucket water content (calculated in simple model)
!	 real(r8), pointer :: mml_sno_grc     				(:)   => null() ! gridcell level snow bucket water content (calculated in simple model)
!	 real(r8), pointer :: mml_alb_grc     				(:)   => null() ! gridcell level albedo (combination bare ground and snow; calculated in simple model)
!	 
!	 ! MML: 2016.01.25 these variables are read from a .nc file, prescribed by user 
!	 ! (there should be the ability to have a different value for each month, but here 
!	 ! I'll just try and save out the value corresponding to the current month... have to figure out how
!	 ! in the actual model subroutine. 
!	 real(r8), pointer :: mml_nc_alb_grc     			(:)   => null() ! gridcell level bare ground albedo (prescribed)
!	 real(r8), pointer :: mml_nc_soild_grc     			(:)   => null() ! gridcell level soil bucket depth... or volume... (prescribed)
!	 real(r8), pointer :: mml_nc_rough_grc     			(:)   => null() ! gridcell level surface roughness length (prescribed)
!	 real(r8), pointer :: mml_nc_evaprs_grc    			(:)   => null() ! gridcell level evaporative resistance (prescribed)

	! ------------------------------------------------------------------------------------
	! MML: 2016.02.29 Variables that I'm actually going to use (delete above once these are implemented
	! ------------------------------------------------------------------------------------
	
	! --- data to be read from .nc file ---
	! direct albedo:
	real(r8), pointer :: mml_nc_alb_gvd_grc     			(:)   => null() ! albedo ground visible direct (prescribed)
	real(r8), pointer :: mml_nc_alb_svd_grc     			(:)   => null() ! albedo snow visible direct (prescribed)
	real(r8), pointer :: mml_nc_alb_gnd_grc     			(:)   => null() ! albedo ground nir direct (prescribed)
	real(r8), pointer :: mml_nc_alb_snd_grc     			(:)   => null() ! albedo snow nir direct (prescribed)
	! diffuse albedo:
	real(r8), pointer :: mml_nc_alb_gvf_grc     			(:)   => null() ! albedo ground visible diffuse (prescribed)
	real(r8), pointer :: mml_nc_alb_svf_grc     			(:)   => null() ! albedo snow visible diffuse (prescribed)
	real(r8), pointer :: mml_nc_alb_gnf_grc     			(:)   => null() ! albedo ground nir diffuse (prescribed)
	real(r8), pointer :: mml_nc_alb_snf_grc     			(:)   => null() ! albedo snow nir diffuse (prescribed)
	! the rest:
	real(r8), pointer :: mml_nc_snowmask_grc     		(:)   => null() ! snow masking depth [kg/m2] (to weight albedo when snow is present)
	real(r8), pointer :: mml_nc_evaprs_grc     			(:)   => null() ! effective stomatal resistance [s/m] <- using Gordon's units for now, possibly convert later [units?]
	real(r8), pointer :: mml_nc_bucket_cap_grc     		(:)   => null() ! soil bucket water capacity [kg/m2]
	real(r8), pointer :: mml_nc_soil_maxice_grc     	(:,:)   => null() ! space x soil depth; maximum freezable water in each soil layer, [kg/m3] (x dz for actual value in layer)
	real(r8), pointer :: mml_nc_soil_levels_grc     	(:,:)   => null() ! depth from surface to each soil level [m]
	real(r8), pointer :: mml_nc_soil_type_grc			(:)	=> null() ! soil type at gridcell. Needed for thermal calculations.
																		  ! for simplicity, set all layers and all grid cells equal. For flexibility, program this way. 
	real(r8), pointer :: mml_nc_roughness_grc			(:)   => null() ! surface roughness length [m], like canopy height.
	real(r8), pointer :: mml_nc_emiss_grc				(:)   => null() ! surface roughness length [m], like canopy height.
	real(r8), pointer :: mml_nc_glcmask_grc				(:)   => null() ! surface roughness length [m], like canopy height.
	real(r8), pointer :: mml_nc_dust_grc				(:,:)   => null() ! surface roughness length [m], like canopy height.
	
	
	! Note: this .nc data really only needs to be read once at the start of the run... for now,
	! I'm going to read it in every time step (slow?) ... theoretically, albedo could be time varying... but on 
	! a monthly scale, not every 30 min! Probably a clever way to say "if month changed, do this ___" else don't bother. 
	
	! --- data from atmosphere to be carried over ---
	! Energy:
	real(r8), pointer :: mml_atm_fsds_grc     	(:)   => null() ! downwelling shortwave radiation [W/m2]
	 real(r8), pointer :: mml_atm_fsdsnd_grc     	(:)   => null() ! direct nir incident solar radiation
	 real(r8), pointer :: mml_atm_fsdsvd_grc     	(:)   => null() ! direct vis incident solar radiation
	 real(r8), pointer :: mml_atm_fsdsni_grc     	(:)   => null() ! diffuse nir incident solar radiation
	 real(r8), pointer :: mml_atm_fsdsvi_grc     	(:)   => null() ! diffuse visible incident solar radiation
	real(r8), pointer :: mml_atm_lwdn_grc     	(:)   => null() ! downwelling longwave radiation [W/m2]
	real(r8), pointer :: mml_atm_zref_grc		(:)	  => null() ! height of reference level of atm [m]
	real(r8), pointer :: mml_atm_tbot_grc		(:)	  => null() ! temperature at reference height [K]
	real(r8), pointer :: mml_atm_thref_grc		(:)	  => null() ! reference temperature theta at reference height [K]
	real(r8), pointer :: mml_atm_qbot_grc		(:)	  => null() ! specific humidity at reference height [kg/kg]
	real(r8), pointer :: mml_atm_uref_grc		(:)	  => null() ! wind speed at reference height [m/s]
	real(r8), pointer :: mml_atm_eref_grc		(:)	  => null() ! vapor pressure at ref height [Pa]
	real(r8), pointer :: mml_atm_pbot_grc		(:)	  => null() ! atmospheric pressure at ref height [Pa]
	real(r8), pointer :: mml_atm_psrf_grc		(:)	  => null() ! atmospheric pressure at surface [Pa]
	real(r8), pointer :: mml_atm_rhomol_grc		(:)	  => null() ! molar density of air at ref height [mol/m3]
	real(r8), pointer :: mml_atm_rhoair_grc		(:)	  => null() ! density of air at ref height [kg/m3]
	real(r8), pointer :: mml_atm_cp_grc			(:)	  => null() ! specific heat of air at const pressure + ref height [J/kg/K]	
	real(r8), pointer :: mml_atm_pco2			(:)	  => null() ! partial pressure of co2
	! Hydrology:
	real(r8), pointer :: mml_atm_prec_liq_grc    (:)   => null() ! liquid precipitation (rain) [mm/s] ! MML 20180615 - bug: used to say m/s, changing to mm/s
	real(r8), pointer :: mml_atm_prec_frz_grc	 (:)   => null() ! frozen precipitation (snow) [mm/s]
	
	! --- land model data that needs to be saved from timestep to timestep
	! land fluxes
	real(r8), pointer :: mml_lnd_ts_grc	 		(:)   => null() ! surface skin temperature [K]
	real(r8), pointer :: mml_lnd_qs_grc	 		(:)   => null() ! surface specific humidity [kg/kg] or [mol/mol]
	real(r8), pointer :: mml_lnd_qa_grc	 		(:)   => null() ! radiative forcing [W/m2] (calculated from swin, lwin, and albedo)
	real(r8), pointer :: mml_lnd_swabs_grc	 	(:)   => null() ! absorbed shortwave radiation [W/m2]
	real(r8), pointer :: mml_lnd_fsr_grc	 	(:)   => null() ! reflected shortwave radiation [W/m2]
	 real(r8), pointer :: mml_lnd_fsrnd_grc	 	(:)   => null() ! reflected shortwave nir direct radiation [W/m2]
	 real(r8), pointer :: mml_lnd_fsrni_grc	 	(:)   => null() ! reflected shortwave nir diffuse radiation [W/m2]
	 real(r8), pointer :: mml_lnd_fsrvd_grc	 	(:)   => null() ! reflected shortwave visible direct radiation [W/m2]
	 real(r8), pointer :: mml_lnd_fsrvi_grc	 	(:)   => null() ! reflected shortwave visible diffuse radiation [W/m2]
	real(r8), pointer :: mml_lnd_lwup_grc	 	(:)   => null() ! emitted longwave radiation [W/m2]
	real(r8), pointer :: mml_lnd_shflx_grc	 	(:)   => null() ! sensible heat flux [W/m2]
	real(r8), pointer :: mml_lnd_lhflx_grc	 	(:)   => null() ! latent heat flux [W/m2]
	real(r8), pointer :: mml_lnd_gsoi_grc	 	(:)   => null() ! heat flux into soil [W/m2]
	real(r8), pointer :: mml_lnd_gsnow_grc	 	(:)   => null() ! heat flux into snow [W/m2] (do I use this/save this?)
	real(r8), pointer :: mml_lnd_evap_grc	 	(:)   => null() ! [kg H2O/ m2 / s = mm/s] -> check!!! 
	real(r8), pointer :: mml_lnd_ustar_grc	 	(:)   => null() ! friction velocity [m/s]
	real(r8), pointer :: mml_lnd_tstar_grc	 	(:)   => null() ! Temperature scale [K]
	real(r8), pointer :: mml_lnd_qstar_grc	 	(:)   => null() ! friction velocity [kg/kg]
	real(r8), pointer :: mml_lnd_tvstar_grc	 	(:)   => null() ! Theta_v_star, virtual potential temperature flux [K]
	real(r8), pointer :: mml_lnd_obu_grc	 	(:)   => null() ! obukhov length [m]
	real(r8), pointer :: mml_lnd_ram_grc	 	(:)   => null() ! aerodynamic resistance for momentum (and moisture) [s/m]
	real(r8), pointer :: mml_lnd_rah_grc	 	(:)   => null() ! aerodynamic resistance for heat [s/m]
	real(r8), pointer :: mml_lnd_res_grc	 	(:)   => null() ! evap_rs (lid resistance) + rah (aerodynamic resistance) [s/m]!
	real(r8), pointer :: mml_lnd_effective_res_grc	 	(:)   => null() ! 1/ beta * ( evap_rs + rah)  [s/m]
	real(r8), pointer :: mml_lnd_beta_grc	 	(:)   => null() ! beta [unitless] from 0 to 1, resistance due to bucket fullness
	real(r8), pointer :: mml_lnd_disp_grc	 	(:)   => null() ! displacement height [m] (MML - what is this? do I use this?)
	real(r8), pointer :: mml_lnd_z0m_grc	 	(:)   => null() ! roughness length for momentum [m]
	real(r8), pointer :: mml_lnd_z0h_grc	 	(:)   => null() ! roughness length for heat [m]
	real(r8), pointer :: mml_lnd_alb_grc	 	(:)   => null() ! actual albedo (visible, direct?) (accounting for snow) [unitless]
	real(r8), pointer :: mml_lnd_fsns_grc	 	(:)   => null() ! flux of shortwave at the surface (in - out), pos into surface
	real(r8), pointer :: mml_lnd_flns_grc	 	(:)   => null() ! flux of longwave at the surface (out - in), pos out of surface
	real(r8), pointer :: mml_lnd_snowmelt	 	(:)   => null() ! flux of longwave at the surface (out - in), pos out of surface
	! soil energy
	real(r8), pointer :: mml_soil_t_grc	 		(:,:)   => null() ! soil temperature in each layer [K]
	real(r8), pointer :: mml_soil_liq_grc	 	(:,:)   => null() ! amount of liquid thermodynamic water [kg/m2]
	real(r8), pointer :: mml_soil_ice_grc	 	(:,:)   => null() ! amount of frozen thermodynamic water [kg/m2]
	real(r8), pointer :: mml_soil_dz_grc     	(:,:)   => null() ! thickness of each soil layer (calculated from depth field) [m]
	real(r8), pointer :: mml_soil_zh_grc     	(:,:)   => null() ! soil depth at interface between each layer [m]
	real(r8), pointer :: mml_soil_tk_grc     	(:,:)   => null() ! thermal conductivity of each soil layer [W/m/K] (this depends on soil properties...)
	real(r8), pointer :: mml_soil_tk_1d_grc     	(:)   => null() ! thermal conductivity of every soil layer [W/m/K] (this depends on soil properties...)
	real(r8), pointer :: mml_soil_tkh_grc     	(:,:)   => null() ! thermal conductivity at base of each soil layer [W/m/K] (this depends on soil properties...)
	real(r8), pointer :: mml_soil_dtsoi_grc     	(:,:)   => null() ! thermal conductivity at base of each soil layer [W/m/K] (this depends on soil properties...)
	real(r8), pointer :: mml_soil_cv_grc     	(:,:)   => null() ! heat capacity of each soil layer [J/m3/K] (this depends on soil properties...)
	real(r8), pointer :: mml_soil_cv_1d_grc     	(:)   => null() ! heat capacity of every soil layer [J/m3/K] (this depends on soil properties...)
	real(r8), pointer :: mml_glc_tk_1d_grc     	(:)   => null() ! thermal conductivity of every ice layer [W/m/K] (this depends on soil properties...)
	real(r8), pointer :: mml_glc_cv_1d_grc     	(:)   => null() ! thermal conductivity of every ice layer [W/m/K] (this depends on soil properties...)
	
	! soil hydrology
	real(r8), pointer :: mml_soil_water_grc	 	(:)   => null() ! soil bucket water content [kg/m2]
	real(r8), pointer :: mml_soil_snow_grc	 	(:)   => null() ! snow bucket water content [kg/m2]
	real(r8), pointer :: mml_soil_runoff_grc	(:)   => null() ! runoff water if soil bucket water > bucket capacity [kg/m2]
	
	! check over-large dew fluxes
	real(r8), pointer :: mml_lh_excess	(:)   => null() ! runoff water if soil bucket water > bucket capacity [kg/m2]
	real(r8), pointer :: mml_q_excess	(:)   => null() ! runoff water if soil bucket water > bucket capacity [kg/m2]
	real(r8), pointer :: mml_lh_demand	(:)   => null() ! runoff water if soil bucket water > bucket capacity [kg/m2]
	real(r8), pointer :: mml_q_demand	(:)   => null() ! runoff water if soil bucket water > bucket capacity [kg/m2]
	
	! --- values passed to atmosphere
	real(r8), pointer :: mml_out_tref2m_grc	(:)   => null() ! 2m air temperature calculated from tsrf and tstar (mo theory)
	real(r8), pointer :: mml_out_qref2m_grc	(:)   => null() ! 2m humidity calculated from qsrf and qstar (mo theory)
	real(r8), pointer :: mml_out_uref10m_grc	(:)   => null() ! 10m wind calculated from ustar (mo theory)
	real(r8), pointer :: mml_out_taux	(:)   => null() ! 10m wind calculated from ustar (mo theory)
	real(r8), pointer :: mml_out_tauy	(:)   => null() ! 10m wind calculated from ustar (mo theory)
	
	! --- temproary diagnostic vars
	real(r8), pointer :: mml_diag1_1d_grc	(:)   => null() ! 1d diagnostic
	real(r8), pointer :: mml_diag2_1d_grc	(:)   => null() ! 1d diagnostic
	real(r8), pointer :: mml_diag3_1d_grc	(:)   => null() ! 1d diagnostic
	real(r8), pointer :: mml_diag1_2d_grc	(:,:)   => null() ! 2d diagnostic
	real(r8), pointer :: mml_diag2_2d_grc	(:,:)   => null() ! 2d diagnostic
	real(r8), pointer :: mml_diag3_2d_grc	(:,:)   => null() ! 2d diagnostic
	
	! --- error values for fluxes/balances
	real(r8), pointer :: mml_err_h2o	(:)   => null() ! total water conservation error
	real(r8), pointer :: mml_err_h2osno	(:)   => null() ! imbalance in snow depth (liquid water)
	real(r8), pointer :: mml_err_seb	(:)   => null() ! surface energy conservation error
	real(r8), pointer :: mml_err_soi	(:)   => null() ! soil/lake energy conservation error
	real(r8), pointer :: mml_err_sol	(:)   => null() ! solar radiation conservation error
	
	
	! Note: variables I'm going to have to add to couple to cam:
	! 2m T, 10m T, 10m U, dust flux, others?
	
	 
	! ------------------------------------------------------------------------------------
	
	 	 
     ! atm->lnd downscaled
     real(r8), pointer :: forc_t_downscaled_col         (:)   => null() ! downscaled atm temperature (Kelvin)
     real(r8), pointer :: forc_th_downscaled_col        (:)   => null() ! downscaled atm potential temperature (Kelvin)
     real(r8), pointer :: forc_q_downscaled_col         (:)   => null() ! downscaled atm specific humidity (kg/kg)
     real(r8), pointer :: forc_pbot_downscaled_col      (:)   => null() ! downscaled atm pressure (Pa)
     real(r8), pointer :: forc_rho_downscaled_col       (:)   => null() ! downscaled atm density (kg/m**3)
     real(r8), pointer :: forc_rain_downscaled_col      (:)   => null() ! downscaled atm rain rate [mm/s]
     real(r8), pointer :: forc_snow_downscaled_col      (:)   => null() ! downscaled atm snow rate [mm/s]
     real(r8), pointer :: forc_lwrad_downscaled_col     (:)   => null() ! downscaled atm downwrd IR longwave radiation (W/m**2)

     !  rof->lnd
     real(r8), pointer :: forc_flood_grc                (:)   => null() ! rof flood (mm/s)
     real(r8), pointer :: volr_grc                      (:)   => null() ! rof volr total volume (m3)
     real(r8), pointer :: volrmch_grc                   (:)   => null() ! rof volr main channel (m3)

     ! anomaly forcing
     real(r8), pointer :: af_precip_grc                 (:)   => null() ! anomaly forcing 
     real(r8), pointer :: af_uwind_grc                  (:)   => null() ! anomaly forcing 
     real(r8), pointer :: af_vwind_grc                  (:)   => null() ! anomaly forcing 
     real(r8), pointer :: af_tbot_grc                   (:)   => null() ! anomaly forcing 
     real(r8), pointer :: af_pbot_grc                   (:)   => null() ! anomaly forcing 
     real(r8), pointer :: af_shum_grc                   (:)   => null() ! anomaly forcing 
     real(r8), pointer :: af_swdn_grc                   (:)   => null() ! anomaly forcing 
     real(r8), pointer :: af_lwdn_grc                   (:)   => null() ! anomaly forcing 
     real(r8), pointer :: bc_precip_grc                 (:)   => null() ! anomaly forcing - add bias correction

     ! time averaged quantities
     real(r8) , pointer :: fsd24_patch                  (:)   => null() ! patch 24hr average of direct beam radiation 
     real(r8) , pointer :: fsd240_patch                 (:)   => null() ! patch 240hr average of direct beam radiation 
     real(r8) , pointer :: fsi24_patch                  (:)   => null() ! patch 24hr average of diffuse beam radiation 
     real(r8) , pointer :: fsi240_patch                 (:)   => null() ! patch 240hr average of diffuse beam radiation 
     real(r8) , pointer :: prec365_col                  (:)   => null() ! col 365-day running mean of tot. precipitation (see comment in UpdateAccVars regarding why this is col-level despite other prec accumulators being patch-level)
     real(r8) , pointer :: prec60_patch                 (:)   => null() ! patch 60-day running mean of tot. precipitation (mm/s) 
     real(r8) , pointer :: prec10_patch                 (:)   => null() ! patch 10-day running mean of tot. precipitation (mm/s) 
     real(r8) , pointer :: rh30_patch                   (:)   => null() ! patch 30-day running mean of relative humidity 
     real(r8) , pointer :: prec24_patch                 (:)   => null() ! patch 24-hour running mean of tot. precipitation (mm/s) 
     real(r8) , pointer :: rh24_patch                   (:)   => null() ! patch 24-hour running mean of relative humidity
     real(r8) , pointer :: wind24_patch                 (:)   => null() ! patch 24-hour running mean of wind
     real(r8) , pointer :: t_mo_patch                   (:)   => null() ! patch 30-day average temperature (Kelvin)
     real(r8) , pointer :: t_mo_min_patch               (:)   => null() ! patch annual min of t_mo (Kelvin)

   contains

     procedure, public  :: Init
     procedure, public  :: InitForTesting  ! version of Init meant for unit testing
     procedure, private :: ReadNamelist
     procedure, private :: InitAllocate
     procedure, private :: InitHistory  
     procedure, private :: InitCold    		! MML 2016.01.15 adding InitCold to give accumulating variables a starting point 
     procedure, public  :: InitAccBuffer
     procedure, public  :: InitAccVars
     procedure, public  :: UpdateAccVars
     procedure, public  :: Restart
     procedure, public  :: Clean

  end type atm2lnd_type

  interface atm2lnd_params_type
     module procedure atm2lnd_params_constructor
  end interface atm2lnd_params_type

  character(len=*), parameter, private :: sourcefile = &
       __FILE__
  !----------------------------------------------------

contains

  !-----------------------------------------------------------------------
  function atm2lnd_params_constructor(repartition_rain_snow, &
       lapse_rate, &
       precip_repartition_glc_all_snow_t, precip_repartition_glc_all_rain_t, &
       precip_repartition_nonglc_all_snow_t, precip_repartition_nonglc_all_rain_t) &
       result(params)
    !
    ! !DESCRIPTION:
    ! Creates a new instance of atm2lnd_params_type
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    type(atm2lnd_params_type) :: params  ! function result
    logical, intent(in) :: repartition_rain_snow

    ! Surface temperature lapse rate (K m-1)
    real(r8), intent(in) :: lapse_rate

    ! End-points of the rain-snow ramp for glacier landunits (degrees C)
    ! Must be present if repartition_rain_snow is true; ignored otherwise
    real(r8), intent(in), optional :: precip_repartition_glc_all_snow_t
    real(r8), intent(in), optional :: precip_repartition_glc_all_rain_t

    ! End-points of the rain-snow ramp for non-glacier landunits (degrees C)
    ! Must be present if repartition_rain_snow is true; ignored otherwise
    real(r8), intent(in), optional :: precip_repartition_nonglc_all_snow_t
    real(r8), intent(in), optional :: precip_repartition_nonglc_all_rain_t
    !
    ! !LOCAL VARIABLES:

    character(len=*), parameter :: subname = 'atm2lnd_params_constructor'
    !-----------------------------------------------------------------------

    params%repartition_rain_snow = repartition_rain_snow

    params%lapse_rate = lapse_rate

    if (repartition_rain_snow) then

       ! Make sure all of the repartitioning-related parameters are present

       if (.not. present(precip_repartition_glc_all_snow_t)) then
          call endrun(subname // &
               ' ERROR: For repartition_rain_snow true, precip_repartition_glc_all_snow_t must be provided')
       end if
       if (.not. present(precip_repartition_glc_all_rain_t)) then
          call endrun(subname // &
               ' ERROR: For repartition_rain_snow true, precip_repartition_glc_all_rain_t must be provided')
       end if
       if (.not. present(precip_repartition_nonglc_all_snow_t)) then
          call endrun(subname // &
               ' ERROR: For repartition_rain_snow true, precip_repartition_nonglc_all_snow_t must be provided')
       end if
       if (.not. present(precip_repartition_nonglc_all_rain_t)) then
          call endrun(subname // &
               ' ERROR: For repartition_rain_snow true, precip_repartition_nonglc_all_rain_t must be provided')
       end if

       ! Do some other error checking

       if (precip_repartition_glc_all_rain_t <= precip_repartition_glc_all_snow_t) then
          call endrun(subname // &
               ' ERROR: Must have precip_repartition_glc_all_snow_t < precip_repartition_glc_all_rain_t')
       end if

       if (precip_repartition_nonglc_all_rain_t <= precip_repartition_nonglc_all_snow_t) then
          call endrun(subname // &
               ' ERROR: Must have precip_repartition_nonglc_all_snow_t < precip_repartition_nonglc_all_rain_t')
       end if

       ! Convert to the form of the parameters we want for the main code

       call compute_ramp_params( &
            all_snow_t_c = precip_repartition_glc_all_snow_t, &
            all_rain_t_c = precip_repartition_glc_all_rain_t, &
            all_snow_t_k = params%precip_repartition_glc_all_snow_t, &
            frac_rain_slope = params%precip_repartition_glc_frac_rain_slope)

       call compute_ramp_params( &
            all_snow_t_c = precip_repartition_nonglc_all_snow_t, &
            all_rain_t_c = precip_repartition_nonglc_all_rain_t, &
            all_snow_t_k = params%precip_repartition_nonglc_all_snow_t, &
            frac_rain_slope = params%precip_repartition_nonglc_frac_rain_slope)

    else  ! .not. repartition_rain_snow
       params%precip_repartition_glc_all_snow_t = nan
       params%precip_repartition_glc_frac_rain_slope = nan
       params%precip_repartition_nonglc_all_snow_t = nan
       params%precip_repartition_nonglc_frac_rain_slope = nan
    end if

  contains
    subroutine compute_ramp_params(all_snow_t_c, all_rain_t_c, &
         all_snow_t_k, frac_rain_slope)
      real(r8), intent(in)  :: all_snow_t_c  ! Temperature at which precip falls entirely as rain (deg C)
      real(r8), intent(in)  :: all_rain_t_c  ! Temperature at which precip falls entirely as snow (deg C)
      real(r8), intent(out) :: all_snow_t_k  ! Temperature at which precip falls entirely as snow (K)
      real(r8), intent(out) :: frac_rain_slope ! Slope of the frac_rain vs. T relationship

      frac_rain_slope = 1._r8 / (all_rain_t_c - all_snow_t_c)
      all_snow_t_k = all_snow_t_c + tfrz
    end subroutine compute_ramp_params

  end function atm2lnd_params_constructor


  !------------------------------------------------------------------------
  subroutine Init(this, bounds, NLFilename)

    class(atm2lnd_type) :: this
    type(bounds_type), intent(in) :: bounds
    character(len=*), intent(in) :: NLFilename ! namelist filename

    call this%InitAllocate(bounds)
    call this%ReadNamelist(NLFilename)
    call this%InitHistory(bounds)
    
    ! MML 2016.01.15 adding call to InitCold (make sure it doesn't keep using the 
    !   coldstart value after the first time step!)
    call this%InitCold(bounds)	
    
  end subroutine Init

  !-----------------------------------------------------------------------
  subroutine InitForTesting(this, bounds, params)
    !
    ! !DESCRIPTION:
    ! Does initialization needed for unit testing. Allows caller to prescribe parameter
    ! values (bypassing the namelist read)
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    class(atm2lnd_type) :: this
    type(bounds_type), intent(in) :: bounds

    ! If params isn't provided, we use default values
    type(atm2lnd_params_type), intent(in), optional :: params
    !
    ! !LOCAL VARIABLES:
    type(atm2lnd_params_type) :: l_params

    character(len=*), parameter :: subname = 'InitForTesting'
    !-----------------------------------------------------------------------

    if (present(params)) then
       l_params = params
    else
       ! Use arbitrary values
       l_params = atm2lnd_params_type( &
            repartition_rain_snow = .false., &
            lapse_rate = 0.01_r8)
    end if

    call this%InitAllocate(bounds)
    this%params = l_params

  end subroutine InitForTesting


  !-----------------------------------------------------------------------
  subroutine ReadNamelist(this, NLFilename)
    !
    ! !DESCRIPTION:
    ! Read the atm2lnd namelist
    !
    ! !USES:
    use fileutils      , only : getavu, relavu, opnfil
    use shr_nl_mod     , only : shr_nl_find_group_name
    use spmdMod        , only : masterproc, mpicom
    use shr_mpi_mod    , only : shr_mpi_bcast
    !
    ! !ARGUMENTS:
    character(len=*), intent(in) :: NLFilename ! Namelist filename
    class(atm2lnd_type), intent(inout) :: this
    !
    ! !LOCAL VARIABLES:

    ! temporary variables corresponding to the components of atm2lnd_params_type
    logical :: repartition_rain_snow
    real(r8) :: lapse_rate
    real(r8) :: precip_repartition_glc_all_snow_t
    real(r8) :: precip_repartition_glc_all_rain_t
    real(r8) :: precip_repartition_nonglc_all_snow_t
    real(r8) :: precip_repartition_nonglc_all_rain_t

    integer :: ierr                 ! error code
    integer :: unitn                ! unit for namelist file
    character(len=*), parameter :: nmlname = 'atm2lnd_inparm'

    character(len=*), parameter :: subname = 'ReadNamelist'
    !-----------------------------------------------------------------------

    namelist /atm2lnd_inparm/ repartition_rain_snow, &
         lapse_rate, &
         precip_repartition_glc_all_snow_t, precip_repartition_glc_all_rain_t, &
         precip_repartition_nonglc_all_snow_t, precip_repartition_nonglc_all_rain_t

    ! Initialize namelist variables to defaults
    repartition_rain_snow = .false.
    lapse_rate = nan
    precip_repartition_glc_all_snow_t = nan
    precip_repartition_glc_all_rain_t = nan
    precip_repartition_nonglc_all_snow_t = nan
    precip_repartition_nonglc_all_rain_t = nan

    if (masterproc) then
       unitn = getavu()
       call opnfil (NLFilename, unitn, 'F')
       call shr_nl_find_group_name(unitn, nmlname, status=ierr)
       if (ierr == 0) then
          read(unitn, nml=atm2lnd_inparm, iostat=ierr)
          if (ierr /= 0) then
             call endrun(msg="ERROR reading "//nmlname//"namelist"//errmsg(sourcefile, __LINE__))
          end if
       else
          call endrun(msg="ERROR could NOT find "//nmlname//"namelist"//errmsg(sourcefile, __LINE__))
       end if
       call relavu( unitn )
    end if

    call shr_mpi_bcast(repartition_rain_snow, mpicom)
    call shr_mpi_bcast(lapse_rate, mpicom)
    call shr_mpi_bcast(precip_repartition_glc_all_snow_t, mpicom)
    call shr_mpi_bcast(precip_repartition_glc_all_rain_t, mpicom)
    call shr_mpi_bcast(precip_repartition_nonglc_all_snow_t, mpicom)
    call shr_mpi_bcast(precip_repartition_nonglc_all_rain_t, mpicom)

    if (masterproc) then
       write(iulog,*) ' '
       write(iulog,*) nmlname//' settings:'
       ! Write settings one-by-one rather than with a nml write because some settings may
       ! be NaN if certain options are turned off.
       write(iulog,*) 'repartition_rain_snow = ', repartition_rain_snow
       write(iulog,*) 'lapse_rate = ', lapse_rate
       if (repartition_rain_snow) then
          write(iulog,*) 'precip_repartition_glc_all_snow_t = ', precip_repartition_glc_all_snow_t
          write(iulog,*) 'precip_repartition_glc_all_rain_t = ', precip_repartition_glc_all_rain_t
          write(iulog,*) 'precip_repartition_nonglc_all_snow_t = ', precip_repartition_nonglc_all_snow_t
          write(iulog,*) 'precip_repartition_nonglc_all_rain_t = ', precip_repartition_nonglc_all_rain_t
       end if
       write(iulog,*) ' '
    end if

    this%params = atm2lnd_params_type( &
         repartition_rain_snow = repartition_rain_snow, &
         lapse_rate = lapse_rate, &
         precip_repartition_glc_all_snow_t = precip_repartition_glc_all_snow_t, &
         precip_repartition_glc_all_rain_t = precip_repartition_glc_all_rain_t, &
         precip_repartition_nonglc_all_snow_t = precip_repartition_nonglc_all_snow_t, &
         precip_repartition_nonglc_all_rain_t = precip_repartition_nonglc_all_rain_t)

  end subroutine ReadNamelist


  !------------------------------------------------------------------------
  subroutine InitAllocate(this, bounds)
    !
    ! !DESCRIPTION:
    ! Initialize atm2lnd derived type
    !
    ! !ARGUMENTS:
    class(atm2lnd_type) :: this
    type(bounds_type), intent(in) :: bounds  
    !
    ! !LOCAL VARIABLES:
    real(r8) :: ival  = 0.0_r8  ! initial value
    integer  :: begg, endg
    integer  :: begc, endc
    integer  :: begp, endp
    !------------------------------------------------------------------------

    begg = bounds%begg; endg= bounds%endg
    begc = bounds%begc; endc= bounds%endc
    begp = bounds%begp; endp= bounds%endp

    ! atm->lnd
    allocate(this%forc_u_grc                    (begg:endg))        ; this%forc_u_grc                    (:)   = ival
    allocate(this%forc_v_grc                    (begg:endg))        ; this%forc_v_grc                    (:)   = ival
    allocate(this%forc_wind_grc                 (begg:endg))        ; this%forc_wind_grc                 (:)   = ival
    allocate(this%forc_rh_grc                   (begg:endg))        ; this%forc_rh_grc                   (:)   = ival
    allocate(this%forc_hgt_grc                  (begg:endg))        ; this%forc_hgt_grc                  (:)   = ival
    allocate(this%forc_topo_grc                 (begg:endg))        ; this%forc_topo_grc                 (:)   = ival
    allocate(this%forc_hgt_u_grc                (begg:endg))        ; this%forc_hgt_u_grc                (:)   = ival
    allocate(this%forc_hgt_t_grc                (begg:endg))        ; this%forc_hgt_t_grc                (:)   = ival
    allocate(this%forc_hgt_q_grc                (begg:endg))        ; this%forc_hgt_q_grc                (:)   = ival
    allocate(this%forc_vp_grc                   (begg:endg))        ; this%forc_vp_grc                   (:)   = ival
    allocate(this%forc_psrf_grc                 (begg:endg))        ; this%forc_psrf_grc                 (:)   = ival
    allocate(this%forc_pco2_grc                 (begg:endg))        ; this%forc_pco2_grc                 (:)   = ival
    allocate(this%forc_solad_grc                (begg:endg,numrad)) ; this%forc_solad_grc                (:,:) = ival
    allocate(this%forc_solai_grc                (begg:endg,numrad)) ; this%forc_solai_grc                (:,:) = ival
    allocate(this%forc_solar_grc                (begg:endg))        ; this%forc_solar_grc                (:)   = ival
    allocate(this%forc_pc13o2_grc               (begg:endg))        ; this%forc_pc13o2_grc               (:)   = ival
    allocate(this%forc_po2_grc                  (begg:endg))        ; this%forc_po2_grc                  (:)   = ival
    allocate(this%forc_aer_grc                  (begg:endg,14))     ; this%forc_aer_grc                  (:,:) = ival

    ! atm->lnd not downscaled
    allocate(this%forc_t_not_downscaled_grc     (begg:endg))        ; this%forc_t_not_downscaled_grc     (:)   = ival
    allocate(this%forc_q_not_downscaled_grc     (begg:endg))        ; this%forc_q_not_downscaled_grc     (:)   = ival
    allocate(this%forc_pbot_not_downscaled_grc  (begg:endg))        ; this%forc_pbot_not_downscaled_grc  (:)   = ival
    allocate(this%forc_th_not_downscaled_grc    (begg:endg))        ; this%forc_th_not_downscaled_grc    (:)   = ival
    allocate(this%forc_rho_not_downscaled_grc   (begg:endg))        ; this%forc_rho_not_downscaled_grc   (:)   = ival
    allocate(this%forc_lwrad_not_downscaled_grc (begg:endg))        ; this%forc_lwrad_not_downscaled_grc (:)   = ival
    allocate(this%forc_rain_not_downscaled_grc  (begg:endg))        ; this%forc_rain_not_downscaled_grc  (:)   = ival
    allocate(this%forc_snow_not_downscaled_grc  (begg:endg))        ; this%forc_snow_not_downscaled_grc  (:)   = ival
    ! MML: 2016.01.14 Adding a variable 2TBOT to see if I can print it to the h0 file
	!allocate(this%forc_2t_not_downscaled_grc    (begg:endg))        ; this%forc_2t_not_downscaled_grc     (:)   = ival       
	 

	! ---------------------------------------
	! MML: 2016.02.29 Allocate space for my land model variables
	
	! .nc data 
	! alb direct
	allocate(this%mml_nc_alb_gvd_grc    	 (begg:endg))     ; this%mml_nc_alb_gvd_grc     	(:)   = ival
	allocate(this%mml_nc_alb_svd_grc    	 (begg:endg))     ; this%mml_nc_alb_svd_grc     	(:)   = ival  
	allocate(this%mml_nc_alb_gnd_grc    	 (begg:endg))     ; this%mml_nc_alb_gnd_grc     	(:)   = ival
	allocate(this%mml_nc_alb_snd_grc    	 (begg:endg))     ; this%mml_nc_alb_snd_grc     	(:)   = ival  
	! alb diffuse
	allocate(this%mml_nc_alb_gvf_grc    	 (begg:endg))     ; this%mml_nc_alb_gvf_grc     	(:)   = ival
	allocate(this%mml_nc_alb_svf_grc    	 (begg:endg))     ; this%mml_nc_alb_svf_grc     	(:)   = ival  
	allocate(this%mml_nc_alb_gnf_grc    	 (begg:endg))     ; this%mml_nc_alb_gnf_grc     	(:)   = ival
	allocate(this%mml_nc_alb_snf_grc    	 (begg:endg))     ; this%mml_nc_alb_snf_grc     	(:)   = ival  
	! the rest
	allocate(this%mml_nc_snowmask_grc    	(begg:endg))     	; this%mml_nc_snowmask_grc    	(:)   = ival 
	allocate(this%mml_nc_evaprs_grc      	(begg:endg))     	; this%mml_nc_evaprs_grc   		(:)   = ival 
	allocate(this%mml_nc_bucket_cap_grc  	(begg:endg))     	; this%mml_nc_bucket_cap_grc   	(:)   = ival 
	allocate(this%mml_nc_soil_maxice_grc    (begg:endg,10))    ; this%mml_nc_soil_maxice_grc   (:,:)   = ival  ! hard coding for 5 soil layers!!!
	allocate(this%mml_nc_soil_levels_grc    (begg:endg,10))    ; this%mml_nc_soil_levels_grc   (:,:)   = ival 
	allocate(this%mml_nc_soil_type_grc    	(begg:endg))        ; this%mml_nc_soil_type_grc   	(:)   = ival  ! one kind of soil per grid cell (all layers equal)
	allocate(this%mml_nc_roughness_grc    	(begg:endg))        ; this%mml_nc_roughness_grc   	(:)   = ival  ! one kind of soil per grid cell (all layers equal)
	allocate(this%mml_nc_emiss_grc    		(begg:endg))        ; this%mml_nc_emiss_grc   		(:)   = ival  ! 
	allocate(this%mml_nc_glcmask_grc    	(begg:endg))        ; this%mml_nc_glcmask_grc   	(:)   = ival  ! 
	allocate(this%mml_nc_dust_grc    		(begg:endg,4))        ; this%mml_nc_dust_grc   		(:,:)   = ival  ! 
	
	! from atm data
	allocate(this%mml_atm_fsds_grc    	 	(begg:endg))     	; this%mml_atm_fsds_grc     	(:)   = ival
	 allocate(this%mml_atm_fsdsnd_grc    	 	(begg:endg))     ; this%mml_atm_fsdsnd_grc     	(:)   = ival
	 allocate(this%mml_atm_fsdsvd_grc    	 	(begg:endg))     ; this%mml_atm_fsdsvd_grc     	(:)   = ival
	 allocate(this%mml_atm_fsdsni_grc    	 	(begg:endg))     ; this%mml_atm_fsdsni_grc     	(:)   = ival
	 allocate(this%mml_atm_fsdsvi_grc    	 	(begg:endg))     ; this%mml_atm_fsdsvi_grc     	(:)   = ival
	allocate(this%mml_atm_lwdn_grc    	 	(begg:endg))     	; this%mml_atm_lwdn_grc     	(:)   = ival
	allocate(this%mml_atm_zref_grc    	 	(begg:endg))     	; this%mml_atm_zref_grc     	(:)   = ival
	allocate(this%mml_atm_tbot_grc    	 	(begg:endg))     	; this%mml_atm_tbot_grc     	(:)   = ival
	allocate(this%mml_atm_thref_grc    	 	(begg:endg))     	; this%mml_atm_thref_grc     	(:)   = ival
	allocate(this%mml_atm_qbot_grc    	 	(begg:endg))     	; this%mml_atm_qbot_grc     	(:)   = ival
	allocate(this%mml_atm_uref_grc    	 	(begg:endg))     	; this%mml_atm_uref_grc     	(:)   = ival
	allocate(this%mml_atm_eref_grc    	 	(begg:endg))     	; this%mml_atm_eref_grc     	(:)   = ival
	allocate(this%mml_atm_pbot_grc    	 	(begg:endg))     	; this%mml_atm_pbot_grc     	(:)   = ival
	allocate(this%mml_atm_psrf_grc    	 	(begg:endg))     	; this%mml_atm_psrf_grc     	(:)   = ival
	allocate(this%mml_atm_rhomol_grc    	(begg:endg))     	; this%mml_atm_rhomol_grc     	(:)   = ival
	allocate(this%mml_atm_rhoair_grc    	(begg:endg))     	; this%mml_atm_rhoair_grc     	(:)   = ival
	allocate(this%mml_atm_cp_grc    	 	(begg:endg))     	; this%mml_atm_cp_grc	     	(:)   = ival
	allocate(this%mml_atm_pco2    	 		(begg:endg))     	; this%mml_atm_pco2	     		(:)   = ival
	allocate(this%mml_atm_prec_liq_grc    	(begg:endg))     	; this%mml_atm_prec_liq_grc     (:)   = ival
	allocate(this%mml_atm_prec_frz_grc    	(begg:endg))     	; this%mml_atm_prec_frz_grc     (:)   = ival
	
	! land flux data
	allocate(this%mml_lnd_ts_grc    	(begg:endg))     	; this%mml_lnd_ts_grc     (:)   = ival
	allocate(this%mml_lnd_qs_grc    	(begg:endg))     	; this%mml_lnd_qs_grc     (:)   = ival
	allocate(this%mml_lnd_qa_grc    	(begg:endg))     	; this%mml_lnd_qa_grc     (:)   = ival
	allocate(this%mml_lnd_swabs_grc    	(begg:endg))     	; this%mml_lnd_swabs_grc  (:)   = ival
	allocate(this%mml_lnd_fsr_grc    	(begg:endg))     	; this%mml_lnd_fsr_grc   (:)   = ival
	 allocate(this%mml_lnd_fsrnd_grc    	(begg:endg))     	; this%mml_lnd_fsrnd_grc   (:)   = ival
	 allocate(this%mml_lnd_fsrni_grc    	(begg:endg))     	; this%mml_lnd_fsrni_grc   (:)   = ival
	 allocate(this%mml_lnd_fsrvd_grc    	(begg:endg))     	; this%mml_lnd_fsrvd_grc   (:)   = ival
	 allocate(this%mml_lnd_fsrvi_grc    	(begg:endg))     	; this%mml_lnd_fsrvi_grc   (:)   = ival
	allocate(this%mml_lnd_lwup_grc    	(begg:endg))     	; this%mml_lnd_lwup_grc   (:)   = ival
	allocate(this%mml_lnd_shflx_grc    	(begg:endg))     	; this%mml_lnd_shflx_grc  (:)   = ival
	allocate(this%mml_lnd_lhflx_grc    	(begg:endg))     	; this%mml_lnd_lhflx_grc  (:)   = ival
	allocate(this%mml_lnd_gsoi_grc    	(begg:endg))     	; this%mml_lnd_gsoi_grc   (:)   = ival
	allocate(this%mml_lnd_gsnow_grc    	(begg:endg))     	; this%mml_lnd_gsnow_grc  (:)   = ival
	allocate(this%mml_lnd_evap_grc    	(begg:endg))     	; this%mml_lnd_evap_grc   (:)   = ival
	allocate(this%mml_lnd_ustar_grc    	(begg:endg))     	; this%mml_lnd_ustar_grc  (:)   = ival
	allocate(this%mml_lnd_tstar_grc    	(begg:endg))     	; this%mml_lnd_tstar_grc  (:)   = ival
	allocate(this%mml_lnd_qstar_grc    	(begg:endg))     	; this%mml_lnd_qstar_grc  (:)   = ival
	allocate(this%mml_lnd_tvstar_grc    (begg:endg))     	; this%mml_lnd_tvstar_grc  (:)   = ival
	allocate(this%mml_lnd_obu_grc    	(begg:endg))     	; this%mml_lnd_obu_grc    (:)   = ival
	allocate(this%mml_lnd_ram_grc    	(begg:endg))     	; this%mml_lnd_ram_grc    (:)   = ival
	allocate(this%mml_lnd_rah_grc    	(begg:endg))     	; this%mml_lnd_rah_grc    (:)   = ival
	allocate(this%mml_lnd_res_grc    	(begg:endg))     	; this%mml_lnd_res_grc    (:)   = ival
	allocate(this%mml_lnd_effective_res_grc    	(begg:endg))     	; this%mml_lnd_effective_res_grc    (:)   = ival
	allocate(this%mml_lnd_beta_grc    	(begg:endg))     	; this%mml_lnd_beta_grc    (:)   = ival
	allocate(this%mml_lnd_disp_grc    	(begg:endg))     	; this%mml_lnd_disp_grc   (:)   = ival
	allocate(this%mml_lnd_z0m_grc    	(begg:endg))     	; this%mml_lnd_z0m_grc    (:)   = ival
	allocate(this%mml_lnd_z0h_grc    	(begg:endg))     	; this%mml_lnd_z0h_grc    (:)   = ival
	allocate(this%mml_lnd_alb_grc    	(begg:endg))     	; this%mml_lnd_alb_grc    (:)   = ival
	allocate(this%mml_lnd_fsns_grc    	(begg:endg))     	; this%mml_lnd_fsns_grc    (:)   = ival
	allocate(this%mml_lnd_flns_grc    	(begg:endg))     	; this%mml_lnd_flns_grc    (:)   = ival
	allocate(this%mml_lnd_snowmelt    	(begg:endg))     	; this%mml_lnd_snowmelt    (:)   = ival
	
	! soil data
	allocate(this%mml_soil_t_grc    	(begg:endg,10))    ; this%mml_soil_t_grc      (:,:)   = ival 
	allocate(this%mml_soil_liq_grc    	(begg:endg,10))    ; this%mml_soil_liq_grc    (:,:)   = ival 
	allocate(this%mml_soil_ice_grc    	(begg:endg,10))    ; this%mml_soil_ice_grc    (:,:)   = ival 
	allocate(this%mml_soil_dz_grc    	(begg:endg,10))    ; this%mml_soil_dz_grc     (:,:)   = ival 
	allocate(this%mml_soil_zh_grc    	(begg:endg,10))    ; this%mml_soil_zh_grc    (:,:)   = ival 
	allocate(this%mml_soil_tk_grc    	(begg:endg,10))    ; this%mml_soil_tk_grc     (:,:)   = ival 
	allocate(this%mml_soil_tk_1d_grc    	(begg:endg))    ; this%mml_soil_tk_1d_grc     (:)   = ival 
	allocate(this%mml_soil_tkh_grc    	(begg:endg,10))    ; this%mml_soil_tkh_grc     (:,:)   = ival 
	allocate(this%mml_soil_dtsoi_grc    	(begg:endg,10))    ; this%mml_soil_dtsoi_grc     (:,:)   = ival 
	allocate(this%mml_soil_cv_grc    	(begg:endg,10))    ; this%mml_soil_cv_grc     (:,:)   = ival 
	allocate(this%mml_soil_cv_1d_grc    	(begg:endg))    ; this%mml_soil_cv_1d_grc     (:)   = ival 
	allocate(this%mml_soil_water_grc    (begg:endg))     	; this%mml_soil_water_grc  (:)   = ival
	allocate(this%mml_soil_snow_grc    	(begg:endg))     	; this%mml_soil_snow_grc   (:)   = ival
	allocate(this%mml_soil_runoff_grc   (begg:endg))     	; this%mml_soil_runoff_grc (:)   = ival
	allocate(this%mml_glc_tk_1d_grc    	(begg:endg))   		; this%mml_glc_tk_1d_grc     (:)   = ival 
	allocate(this%mml_glc_cv_1d_grc    	(begg:endg))    	; this%mml_glc_cv_1d_grc     (:)   = ival 
	
	! outbound lnd2atm variables
	allocate(this%mml_out_tref2m_grc	(begg:endg))     	; this%mml_out_tref2m_grc (:)   = ival
	allocate(this%mml_out_qref2m_grc	(begg:endg))     	; this%mml_out_qref2m_grc (:)   = ival
	allocate(this%mml_out_uref10m_grc	(begg:endg))     	; this%mml_out_uref10m_grc (:)   = ival
	allocate(this%mml_out_taux			(begg:endg))     	; this%mml_out_taux (:)   = ival
	allocate(this%mml_out_tauy			(begg:endg))     	; this%mml_out_tauy (:)   = ival
	
	! check over-large dew fluxes
	allocate(this%mml_lh_excess			(begg:endg))     	; this%mml_lh_excess (:)   = ival
	allocate(this%mml_q_excess			(begg:endg))     	; this%mml_q_excess (:)   = ival
	allocate(this%mml_lh_demand			(begg:endg))     	; this%mml_lh_demand (:)   = ival
	allocate(this%mml_q_demand			(begg:endg))     	; this%mml_q_demand (:)   = ival
	
	! Add a few temporary diagnostic variables (no units, get rid of them later)
	allocate(this%mml_diag1_1d_grc	(begg:endg))     	; this%mml_diag1_1d_grc (:)   = ival
	allocate(this%mml_diag2_1d_grc	(begg:endg))     	; this%mml_diag2_1d_grc (:)   = ival
	allocate(this%mml_diag3_1d_grc	(begg:endg))     	; this%mml_diag3_1d_grc (:)   = ival
	
	allocate(this%mml_diag1_2d_grc	(begg:endg,10))     ; this%mml_diag1_2d_grc (:,:)   = ival
	allocate(this%mml_diag2_2d_grc	(begg:endg,10))     ; this%mml_diag2_2d_grc (:,:)   = ival
	allocate(this%mml_diag3_2d_grc	(begg:endg,10))     ; this%mml_diag3_2d_grc (:,:)   = ival
	
	! Error/balance check vars
	allocate(this%mml_err_h2o	(begg:endg))     	; this%mml_err_h2o (:)   = ival
	allocate(this%mml_err_h2osno(begg:endg))    	; this%mml_err_h2osno (:)   = ival
	allocate(this%mml_err_seb	(begg:endg))     	; this%mml_err_seb (:)   = ival
	allocate(this%mml_err_soi	(begg:endg))   	 	; this%mml_err_soi (:)   = ival
	allocate(this%mml_err_sol	(begg:endg))     	; this%mml_err_sol (:)   = ival
	
	
	
	
	! ---------------------------------------

    
    ! atm->lnd downscaled
    allocate(this%forc_t_downscaled_col         (begc:endc))        ; this%forc_t_downscaled_col         (:)   = ival
    allocate(this%forc_q_downscaled_col         (begc:endc))        ; this%forc_q_downscaled_col         (:)   = ival
    allocate(this%forc_pbot_downscaled_col      (begc:endc))        ; this%forc_pbot_downscaled_col      (:)   = ival
    allocate(this%forc_th_downscaled_col        (begc:endc))        ; this%forc_th_downscaled_col        (:)   = ival
    allocate(this%forc_rho_downscaled_col       (begc:endc))        ; this%forc_rho_downscaled_col       (:)   = ival
    allocate(this%forc_lwrad_downscaled_col     (begc:endc))        ; this%forc_lwrad_downscaled_col     (:)   = ival
    allocate(this%forc_rain_downscaled_col      (begc:endc))        ; this%forc_rain_downscaled_col      (:)   = ival
    allocate(this%forc_snow_downscaled_col      (begc:endc))        ; this%forc_snow_downscaled_col      (:)   = ival

    ! rof->lnd
    allocate(this%forc_flood_grc                (begg:endg))        ; this%forc_flood_grc                (:)   = ival
    allocate(this%volr_grc                      (begg:endg))        ; this%volr_grc                      (:)   = ival
    allocate(this%volrmch_grc                   (begg:endg))        ; this%volrmch_grc                   (:)   = ival

    ! anomaly forcing
    allocate(this%bc_precip_grc                 (begg:endg))        ; this%bc_precip_grc                 (:)   = ival
    allocate(this%af_precip_grc                 (begg:endg))        ; this%af_precip_grc                 (:)   = ival
    allocate(this%af_uwind_grc                  (begg:endg))        ; this%af_uwind_grc                  (:)   = ival
    allocate(this%af_vwind_grc                  (begg:endg))        ; this%af_vwind_grc                  (:)   = ival
    allocate(this%af_tbot_grc                   (begg:endg))        ; this%af_tbot_grc                   (:)   = ival
    allocate(this%af_pbot_grc                   (begg:endg))        ; this%af_pbot_grc                   (:)   = ival
    allocate(this%af_shum_grc                   (begg:endg))        ; this%af_shum_grc                   (:)   = ival
    allocate(this%af_swdn_grc                   (begg:endg))        ; this%af_swdn_grc                   (:)   = ival
    allocate(this%af_lwdn_grc                   (begg:endg))        ; this%af_lwdn_grc                   (:)   = ival

    allocate(this%fsd24_patch                   (begp:endp))        ; this%fsd24_patch                   (:)   = nan
    allocate(this%fsd240_patch                  (begp:endp))        ; this%fsd240_patch                  (:)   = nan
    allocate(this%fsi24_patch                   (begp:endp))        ; this%fsi24_patch                   (:)   = nan
    allocate(this%fsi240_patch                  (begp:endp))        ; this%fsi240_patch                  (:)   = nan
    allocate(this%prec10_patch                  (begp:endp))        ; this%prec10_patch                  (:)   = nan
    allocate(this%prec60_patch                  (begp:endp))        ; this%prec60_patch                  (:)   = nan
    allocate(this%rh30_patch                    (begp:endp))        ; this%rh30_patch                    (:)   = nan 
    allocate(this%prec365_col                   (begc:endc))        ; this%prec365_col                   (:)   = nan
    allocate(this%t_mo_patch                    (begp:endp))        ; this%t_mo_patch               (:)   = nan
    allocate(this%t_mo_min_patch                (begp:endp))        ; this%t_mo_min_patch           (:)   = spval ! TODO - initialize this elsewhere

  end subroutine InitAllocate

  !------------------------------------------------------------------------
  subroutine InitHistory(this, bounds)
    !
    ! !USES:
    ! use histFileMod, only : hist_addfld1d 
    ! MML:
    use histFileMod, only : hist_addfld1d, hist_addfld2d
    !
    ! !ARGUMENTS:
    class(atm2lnd_type) :: this
    type(bounds_type), intent(in) :: bounds  
    !
    ! !LOCAL VARIABLES:
    integer  :: begg, endg
    integer  :: begc, endc
    integer  :: begp, endp
    
    integer  :: mml_nsoi ! number of soil levels
    !---------------------------------------------------------------------

	mml_nsoi = 10
	
    begg = bounds%begg; endg= bounds%endg
    begc = bounds%begc; endc= bounds%endc
    begp = bounds%begp; endp= bounds%endp
	
	!write(iulog,*)  'MML trying write h0 - start'
	
    this%forc_flood_grc(begg:endg) = spval
    call hist_addfld1d (fname='QFLOOD',  units='mm/s',  &
         avgflag='A', long_name='runoff from river flooding', &
         ptr_lnd=this%forc_flood_grc)

    this%volr_grc(begg:endg) = spval
    call hist_addfld1d (fname='VOLR',  units='m3',  &
         avgflag='A', long_name='river channel total water storage', &
         ptr_lnd=this%volr_grc)

    this%volrmch_grc(begg:endg) = spval
    call hist_addfld1d (fname='VOLRMCH',  units='m3',  &
         avgflag='A', long_name='river channel main channel water storage', &
         ptr_lnd=this%volrmch_grc)

    this%forc_wind_grc(begg:endg) = spval
    call hist_addfld1d (fname='WIND', units='m/s',  &
         avgflag='A', long_name='atmospheric wind velocity magnitude', &
         ptr_lnd=this%forc_wind_grc)
    ! Rename of WIND for Urban intercomparision project
    call hist_addfld1d (fname='Wind', units='m/s',  &
         avgflag='A', long_name='atmospheric wind velocity magnitude', &
         ptr_gcell=this%forc_wind_grc, default = 'inactive')

    this%forc_hgt_grc(begg:endg) = spval
    call hist_addfld1d (fname='ZBOT', units='m',  &
         avgflag='A', long_name='atmospheric reference height', &
         ptr_lnd=this%forc_hgt_grc)

    this%forc_topo_grc(begg:endg) = spval
    call hist_addfld1d (fname='ATM_TOPO', units='m', &
         avgflag='A', long_name='atmospheric surface height', &
         ptr_lnd=this%forc_topo_grc)

    this%forc_solar_grc(begg:endg) = spval
    call hist_addfld1d (fname='FSDS', units='W/m^2',  &
         avgflag='A', long_name='atmospheric incident solar radiation', &
         ptr_lnd=this%forc_solar_grc)

    this%forc_pco2_grc(begg:endg) = spval
    call hist_addfld1d (fname='PCO2', units='Pa',  &
         avgflag='A', long_name='atmospheric partial pressure of CO2', &
         ptr_lnd=this%forc_pco2_grc)

    this%forc_solar_grc(begg:endg) = spval
    call hist_addfld1d (fname='SWdown', units='W/m^2',  &
         avgflag='A', long_name='atmospheric incident solar radiation', &
         ptr_gcell=this%forc_solar_grc, default='inactive')

    this%forc_rh_grc(begg:endg) = spval
    call hist_addfld1d (fname='RH', units='%',  &
         avgflag='A', long_name='atmospheric relative humidity', &
         ptr_gcell=this%forc_rh_grc, default='inactive')

    this%forc_t_not_downscaled_grc(begg:endg) = spval
    call hist_addfld1d (fname='Tair_from_atm', units='K',  &
         avgflag='A', long_name='atmospheric air temperature received from atmosphere (pre-downscaling)', &
         ptr_gcell=this%forc_t_not_downscaled_grc, default='inactive')

    this%forc_t_downscaled_col(begc:endc) = spval
    call hist_addfld1d (fname='TBOT', units='K',  &
         avgflag='A', long_name='atmospheric air temperature (downscaled to columns in glacier regions)', &
         ptr_col=this%forc_t_downscaled_col)
    call hist_addfld1d (fname='Tair', units='K', &
         avgflag='A', long_name='atmospheric air temperature (downscaled to columns in glacier regions)', &
         ptr_col=this%forc_t_downscaled_col, default='inactive')

    this%forc_pbot_downscaled_col(begc:endc) = spval
    call hist_addfld1d (fname='PBOT', units='Pa',  &
         avgflag='A', long_name='atmospheric pressure at surface (downscaled to columns in glacier regions)', &
         ptr_col=this%forc_pbot_downscaled_col)
    call hist_addfld1d (fname='PSurf', units='Pa',  &
         avgflag='A', long_name='atmospheric pressure at surface (downscaled to columns in glacier regions)', &
         ptr_col=this%forc_pbot_downscaled_col, default='inactive')

    this%forc_lwrad_downscaled_col(begc:endc) = spval
    call hist_addfld1d (fname='FLDS', units='W/m^2',  &
         avgflag='A', long_name='atmospheric longwave radiation (downscaled to columns in glacier regions)', &
         ptr_col=this%forc_lwrad_downscaled_col)
    call hist_addfld1d (fname='LWdown', units='W/m^2',  &
         avgflag='A', long_name='atmospheric longwave radiation (downscaled to columns in glacier regions)', &
         ptr_col=this%forc_lwrad_downscaled_col, default='inactive')

    this%forc_rain_not_downscaled_grc(begg:endg) = spval
    call hist_addfld1d (fname='RAIN_FROM_ATM', units='mm/s',  &
         avgflag='A', long_name='atmospheric rain received from atmosphere (pre-repartitioning)', &
         ptr_lnd=this%forc_rain_not_downscaled_grc)

    this%forc_snow_not_downscaled_grc(begg:endg) = spval
    call hist_addfld1d (fname='SNOW_FROM_ATM', units='mm/s',  &
         avgflag='A', long_name='atmospheric snow received from atmosphere (pre-repartitioning)', &
         ptr_lnd=this%forc_snow_not_downscaled_grc)

    this%forc_rain_downscaled_col(begc:endc) = spval
    call hist_addfld1d (fname='RAIN', units='mm/s',  &
         avgflag='A', long_name='atmospheric rain, after rain/snow repartitioning based on temperature', &
         ptr_col=this%forc_rain_downscaled_col)
    call hist_addfld1d (fname='Rainf', units='mm/s',  &
         avgflag='A', long_name='atmospheric rain, after rain/snow repartitioning based on temperature', &
         ptr_col=this%forc_rain_downscaled_col, default='inactive')

    this%forc_snow_downscaled_col(begc:endc) = spval
    call hist_addfld1d (fname='SNOW', units='mm/s',  &
         avgflag='A', long_name='atmospheric snow, after rain/snow repartitioning based on temperature', &
         ptr_col=this%forc_snow_downscaled_col)

    this%forc_th_downscaled_col(begc:endc) = spval
    call hist_addfld1d (fname='THBOT', units='K',  &
         avgflag='A', long_name='atmospheric air potential temperature (downscaled to columns in glacier regions)', &
         ptr_col=this%forc_th_downscaled_col)
         
!	! MML: 2016.01.14 Try and add a new history field variable equal to 2TBOT
!	! (just to see if it will print)
!	this%forc_2t_not_downscaled_grc(begg:endg) = spval
!    call hist_addfld1d (fname='T2BOT', units='K',  &
!         avgflag='A', long_name='2x atmospheric air temperature MML Test', &
!         ptr_lnd=this%forc_2t_not_downscaled_grc)
	
    this%forc_q_downscaled_col(begc:endc) = spval
    call hist_addfld1d (fname='QBOT', units='kg/kg',  &
         avgflag='A', long_name='atmospheric specific humidity (downscaled to columns in glacier regions)', &
         ptr_col=this%forc_q_downscaled_col)
    ! Rename of QBOT for Urban intercomparison project
    call hist_addfld1d (fname='Qair', units='kg/kg',  &
         avgflag='A', long_name='atmospheric specific humidity (downscaled to columns in glacier regions)', &
         ptr_col=this%forc_q_downscaled_col, default='inactive')

        
    !-----------------------------------------------------------------------
    ! MML: 2016.01.14 Simple Land Energy and Hydrology variables (gridscale)
    
    ! write(iulog,*)  'MML write to h0: netcdf vars '
    
    ! From .nc file:
    ! (don't typically print these - waste of space. But for now, could be useful... 
    ! Skipping because I'm lazy; add later
!    this%mml_nc_alb_grc(begg:endg) = spval
!    call hist_addfld1d (fname='MML_albdeo', units='unitless',  &
!         avgflag='A', long_name='MML prescribed snow-free surface albedo', &
!         ptr_lnd=this%mml_nc_alb_grc)
!         
!	this%mml_nc_snoalb_grc(begg:endg) = spval
!    call hist_addfld1d (fname='MML_snow_albdeo', units='unitless',  &
!         avgflag='A', long_name='MML prescribed deep-snow albedo', &
!         ptr_lnd=this%mml_nc_snoalb_grc)
    
         
	this%mml_nc_snowmask_grc(begg:endg) = spval
    call hist_addfld1d (fname='MML_snowmaskdepth', units='kg/m2',  &
         avgflag='A', long_name='MML snow required to toggle deep snow albedo', &
         ptr_lnd=this%mml_nc_snowmask_grc)
    
    this%mml_nc_evaprs_grc(begg:endg) = spval
    call hist_addfld1d (fname='MML_evap_rs', units='s/m',  &
         avgflag='A', long_name='MML like stomatal resistance of soil', &
         ptr_lnd=this%mml_nc_evaprs_grc)
    
    this%mml_nc_bucket_cap_grc(begg:endg) = spval
    call hist_addfld1d (fname='MML_bucket_cap', units='kg/m2',  &
         avgflag='A', long_name='MML soil water bucket capacity (maximum water soil can hold)', &
         ptr_lnd=this%mml_nc_bucket_cap_grc)
    
    !this%mml_nc_soil_maxice_grc(begg:endg,:) = spval
    !call hist_addfld1d (fname='MML_maxice', units='kg/m3',  &
    !     avgflag='A', long_name='MML maximum freezable water in each soil layer; for thermal calculations', &
    !     ptr_lnd=this%mml_nc_soil_maxice_grc)
    !data2dptr => this%mml_nc_soil_maxice_grc(begg:endg,:)
    !fieldname = 'MML_maxice'
    !longname =  'MML maximum freezable water in each soil layer; for thermal calculations'
    
    this%mml_nc_soil_type_grc(begg:endg) = spval
    call hist_addfld1d (fname='MML_soiltype', units='unitless',  &
         avgflag='A', long_name='MML Soil type (sand/clay), of 11 possible types; for thermal calculations', &
         ptr_lnd=this%mml_nc_soil_type_grc)
    
    this%mml_nc_roughness_grc(begg:endg) = spval
    call hist_addfld1d (fname='MML_roughness', units='m',  &
         avgflag='A', long_name='MML surface roughness length (e.g. canopy height) ', &
         ptr_lnd=this%mml_nc_roughness_grc)
    
    this%mml_nc_emiss_grc(begg:endg) = spval
    call hist_addfld1d (fname='MML_emiss', units='unitless',  &
         avgflag='A', long_name='MML surface emissivity ', &
         ptr_lnd=this%mml_nc_emiss_grc)
         
    this%mml_nc_glcmask_grc(begg:endg) = spval
    call hist_addfld1d (fname='MML_glcmask', units='unitless',  &
         avgflag='A', long_name='MML logical mask saying which cells should be treated as glaciers', &
         ptr_lnd=this%mml_nc_glcmask_grc)     
    
   
         
    ! Carried from atmosphere:
   ! write(iulog,*)  'MML write to h0: atm vars '
    
    this%mml_atm_fsds_grc(begg:endg) = spval
    call hist_addfld1d (fname='MML_fsds', units='W/m2',  &
         avgflag='A', long_name='MML incoming shortwave radiation', &
         ptr_lnd=this%mml_atm_fsds_grc)
    
    this%mml_atm_fsdsnd_grc(begg:endg) = spval
    call hist_addfld1d (fname='MML_fsdsnd', units='W/m2',  &
         avgflag='A', long_name='MML incoming shortwave nir direct radiation', &
         ptr_lnd=this%mml_atm_fsdsnd_grc)
	
	this%mml_atm_fsdsvd_grc(begg:endg) = spval
    call hist_addfld1d (fname='MML_fsdsvd', units='W/m2',  &
         avgflag='A', long_name='MML incoming shortwave visible direct radiation', &
         ptr_lnd=this%mml_atm_fsdsvd_grc)
    
    this%mml_atm_fsdsni_grc(begg:endg) = spval
    call hist_addfld1d (fname='MML_fsdsni', units='W/m2',  &
         avgflag='A', long_name='MML incoming shortwave nir diffuse radiation', &
         ptr_lnd=this%mml_atm_fsdsni_grc)
    
    this%mml_atm_fsdsvi_grc(begg:endg) = spval
    call hist_addfld1d (fname='MML_fsdsvi', units='W/m2',  &
         avgflag='A', long_name='MML incoming shortwave visible diffuse radiation', &
         ptr_lnd=this%mml_atm_fsdsvi_grc)
         
    this%mml_atm_lwdn_grc(begg:endg) = spval
    call hist_addfld1d (fname='MML_lwdn', units='W/m2',  &
         avgflag='A', long_name='MML incoming longwave radiation', &
         ptr_lnd=this%mml_atm_lwdn_grc)
    
    this%mml_atm_zref_grc(begg:endg) = spval
    call hist_addfld1d (fname='MML_zref', units='m',  &
         avgflag='A', long_name='MML height of atm reference level', &
         ptr_lnd=this%mml_atm_zref_grc)
    
    this%mml_atm_tbot_grc(begg:endg) = spval
    call hist_addfld1d (fname='MML_tbot', units='K',  &
         avgflag='A', long_name='MML temperature midpoint of lowest atm layer', &
         ptr_lnd=this%mml_atm_tbot_grc)
    
    this%mml_atm_thref_grc(begg:endg) = spval
    call hist_addfld1d (fname='MML_thref', units='K',  &
         avgflag='A', long_name='MML potential temperature theta at reference height', &
         ptr_lnd=this%mml_atm_thref_grc)
    
    this%mml_atm_qbot_grc(begg:endg) = spval
    call hist_addfld1d (fname='MML_qbot', units='kg/kg',  &
         avgflag='A', long_name='MML specific humidity midpoint of lowest atm layer', &
         ptr_lnd=this%mml_atm_qbot_grc)
    
    this%mml_atm_uref_grc(begg:endg) = spval
    call hist_addfld1d (fname='MML_uref', units='m/s',  &
         avgflag='A', long_name='MML wind speed at reference height', &
         ptr_lnd=this%mml_atm_uref_grc)
    
    this%mml_atm_eref_grc(begg:endg) = spval
    call hist_addfld1d (fname='MML_eref', units='Pa',  &
         avgflag='A', long_name='MML vapor pressure at reference height', &
         ptr_lnd=this%mml_atm_eref_grc)
    
    this%mml_atm_pbot_grc(begg:endg) = spval
    call hist_addfld1d (fname='MML_pbot', units='Pa',  &
         avgflag='A', long_name='MML atmospheric pressure midpoint of lowest atm layer', &
         ptr_lnd=this%mml_atm_pbot_grc)
    
    this%mml_atm_psrf_grc(begg:endg) = spval
    call hist_addfld1d (fname='MML_psrf', units='Pa',  &
         avgflag='A', long_name='MML atmospheric pressure surface', &
         ptr_lnd=this%mml_atm_psrf_grc)
    
    this%mml_atm_rhomol_grc(begg:endg) = spval
    call hist_addfld1d (fname='MML_rhomol', units='mol/m3',  &
         avgflag='A', long_name='MML molar density of air at reference height', &
         ptr_lnd=this%mml_atm_rhomol_grc)
    
    this%mml_atm_rhoair_grc(begg:endg) = spval
    call hist_addfld1d (fname='MML_rhoair', units='kg/m3',  &
         avgflag='A', long_name='MML mass density of air at reference height', &
         ptr_lnd=this%mml_atm_rhoair_grc)
         
    this%mml_atm_cp_grc(begg:endg) = spval
    call hist_addfld1d (fname='MML_cpair', units='J/kg/K',  &
         avgflag='A', long_name='MML specific heat of air at constant pressure at ref height', &
         ptr_lnd=this%mml_atm_cp_grc)
    
    this%mml_atm_pco2(begg:endg) = spval
    call hist_addfld1d (fname='MML_pco2', units='Pa',  &
         avgflag='A', long_name='MML partial pressure of co2', &
         ptr_lnd=this%mml_atm_pco2)
         
    this%mml_atm_prec_liq_grc(begg:endg) = spval
    call hist_addfld1d (fname='MML_prec_liq', units='mm/s',  &	! or mm/s? 
         avgflag='A', long_name='MML rate of liquid precipitation (rain)', &
         ptr_lnd=this%mml_atm_prec_liq_grc)
    
    this%mml_atm_prec_frz_grc(begg:endg) = spval
    call hist_addfld1d (fname='MML_prec_frz', units='mm/s',  &
         avgflag='A', long_name='MML rate of frozen precipitation (snow)', &
         ptr_lnd=this%mml_atm_prec_frz_grc)
    
    ! Land calculated surface variables
    !write(iulog,*)  'MML write to h0: 1d land vars '
    
    this%mml_lnd_ts_grc(begg:endg) = spval
    call hist_addfld1d (fname='MML_ts', units='K',  &
         avgflag='A', long_name='MML surface skin temperature', &
         ptr_lnd=this%mml_lnd_ts_grc)
    
    this%mml_lnd_qs_grc(begg:endg) = spval
    call hist_addfld1d (fname='MML_qs', units='kg/kg',  &
         avgflag='A', long_name='surface specific humidity [kg/kg] or [mol/mol]', &
         ptr_lnd=this%mml_lnd_qs_grc)
    
    this%mml_lnd_qa_grc(begg:endg) = spval
    call hist_addfld1d (fname='MML_qa', units='W/m2',  &
         avgflag='A', long_name='MML radiative forcing (SWin*albedo + LWin)', &
         ptr_lnd=this%mml_lnd_qa_grc)
    
    this%mml_lnd_swabs_grc(begg:endg) = spval
    call hist_addfld1d (fname='MML_swabs', units='W/m2',  &
         avgflag='A', long_name='MML absorbed shortwave radiation', &
         ptr_lnd=this%mml_lnd_swabs_grc)
    
    this%mml_lnd_fsr_grc(begg:endg) = spval
    call hist_addfld1d (fname='MML_fsr', units='W/m2',  &
         avgflag='A', long_name='MML reflected shortwave radation', &
         ptr_lnd=this%mml_lnd_fsr_grc)
    
    this%mml_lnd_fsrnd_grc(begg:endg) = spval
    call hist_addfld1d (fname='MML_fsrnd', units='W/m2',  &
         avgflag='A', long_name='MML reflected shortwave nir direct radation', &
         ptr_lnd=this%mml_lnd_fsrnd_grc)
    
    this%mml_lnd_fsrni_grc(begg:endg) = spval
    call hist_addfld1d (fname='MML_fsrni', units='W/m2',  &
         avgflag='A', long_name='MML reflected shortwave nir diffuse radation', &
         ptr_lnd=this%mml_lnd_fsrni_grc)
    
    this%mml_lnd_fsrvd_grc(begg:endg) = spval
    call hist_addfld1d (fname='MML_fsrvd', units='W/m2',  &
         avgflag='A', long_name='MML reflected shortwave visible direct radation', &
         ptr_lnd=this%mml_lnd_fsrvd_grc)
    
    this%mml_lnd_fsrvi_grc(begg:endg) = spval
    call hist_addfld1d (fname='MML_fsrvi', units='W/m2',  &
         avgflag='A', long_name='MML reflected shortwave visible diffuse radation', &
         ptr_lnd=this%mml_lnd_fsrvi_grc)
    
    this%mml_lnd_lwup_grc(begg:endg) = spval
    call hist_addfld1d (fname='MML_lwup', units='W/m2',  &
         avgflag='A', long_name='MML emitted longwave radiation', &
         ptr_lnd=this%mml_lnd_lwup_grc)
    
    this%mml_lnd_shflx_grc(begg:endg) = spval
    call hist_addfld1d (fname='MML_shflx', units='W/m2',  &
         avgflag='A', long_name='MML sensible heat flux', &
         ptr_lnd=this%mml_lnd_shflx_grc)
    
    this%mml_lnd_lhflx_grc(begg:endg) = spval
    call hist_addfld1d (fname='MML_lhflx', units='W/m2',  &
         avgflag='A', long_name='MML latent heat flux', &
         ptr_lnd=this%mml_lnd_lhflx_grc)
    
    this%mml_lnd_gsoi_grc(begg:endg) = spval
    call hist_addfld1d (fname='MML_gsoi', units='W/m2',  &
         avgflag='A', long_name='MML flux of energy into the soil', &
         ptr_lnd=this%mml_lnd_gsoi_grc)
    
    this%mml_lnd_gsnow_grc(begg:endg) = spval
    call hist_addfld1d (fname='MML_gsnow', units='W/m2',  &
         avgflag='A', long_name='MML flux of energy into snowmelt', &
         ptr_lnd=this%mml_lnd_gsnow_grc)
    
    this%mml_lnd_evap_grc(begg:endg) = spval
    call hist_addfld1d (fname='MML_evap', units='kg H20 / m2 / s = mm/s',  &
         avgflag='A', long_name='MML evapotranspiration (in kg water over whole time step)', &
         ptr_lnd=this%mml_lnd_evap_grc)
    
    this%mml_lnd_ustar_grc(begg:endg) = spval
    call hist_addfld1d (fname='MML_ustar', units='m/s',  &
         avgflag='A', long_name='MML friction velocity from MO theory', &
         ptr_lnd=this%mml_lnd_ustar_grc)
    
    this%mml_lnd_tstar_grc(begg:endg) = spval
    call hist_addfld1d (fname='MML_tstar', units='K',  &
         avgflag='A', long_name='MML temperature scale from MO theory', &
         ptr_lnd=this%mml_lnd_tstar_grc)
         
    this%mml_lnd_qstar_grc(begg:endg) = spval
    call hist_addfld1d (fname='MML_qstar', units='kg/kg',  &
         avgflag='A', long_name='MML humidity scale (?) from MO theory', &
         ptr_lnd=this%mml_lnd_qstar_grc)
         
    this%mml_lnd_tvstar_grc(begg:endg) = spval
    call hist_addfld1d (fname='MML_tvstar', units='K',  &
         avgflag='A', long_name='MML virtual potential temperature scale from MO theory', &
         ptr_lnd=this%mml_lnd_tvstar_grc)
    
    this%mml_lnd_obu_grc(begg:endg) = spval
    call hist_addfld1d (fname='MML_obu', units='m',  &
         avgflag='A', long_name='MML Obukhov length from MO theory', &
         ptr_lnd=this%mml_lnd_obu_grc)
    
    this%mml_lnd_ram_grc(begg:endg) = spval
    call hist_addfld1d (fname='MML_ram', units='s/m',  &
         avgflag='A', long_name='MML aerodynamic resistance for momentum (and moisture; from MO theory)', &
         ptr_lnd=this%mml_lnd_ram_grc)
    
    this%mml_lnd_rah_grc(begg:endg) = spval
    call hist_addfld1d (fname='MML_rah', units='s/m',  &
         avgflag='A', long_name='MML aerodynamic resistance for heat', &
         ptr_lnd=this%mml_lnd_rah_grc)
         
    this%mml_lnd_res_grc(begg:endg) = spval
    call hist_addfld1d (fname='MML_res_tot', units='s/m',  &
         avgflag='A', long_name='MML lid resistance + aerodynamic resistance for heat (MML_evap_rs + MML_rah)', &
         ptr_lnd=this%mml_lnd_res_grc)

    this%mml_lnd_effective_res_grc(begg:endg) = spval
    call hist_addfld1d (fname='MML_res_effective', units='s/m',  &
         avgflag='A', long_name='MML effective surface resistance = 1/beta * (MML_evap_rs + MML_rah)', &
         ptr_lnd=this%mml_lnd_effective_res_grc)

    this%mml_lnd_beta_grc(begg:endg) = spval
    call hist_addfld1d (fname='MML_beta', units='unitless',  &
         avgflag='A', long_name='MML beta factor for resistance due to bucket emptiness (between 0 and 1)', &
         ptr_lnd=this%mml_lnd_beta_grc)
                     
    this%mml_lnd_z0m_grc(begg:endg) = spval
    call hist_addfld1d (fname='MML_z0m', units='m',  &
         avgflag='A', long_name='MML roughness length for momentum', &
         ptr_lnd=this%mml_lnd_z0m_grc)
    
    this%mml_lnd_z0h_grc(begg:endg) = spval
    call hist_addfld1d (fname='MML_z0h', units='m',  &
         avgflag='A', long_name='MML roughness length for heat', &
         ptr_lnd=this%mml_lnd_z0h_grc)
    
    this%mml_lnd_alb_grc(begg:endg) = spval
    call hist_addfld1d (fname='MML_alb', units='unitless',  &
         avgflag='A', long_name='MML actual albedo (accounting for snow) used', &
         ptr_lnd=this%mml_lnd_alb_grc)
    
    this%mml_lnd_fsns_grc(begg:endg) = spval
    call hist_addfld1d (fname='MML_fsns', units='W/m2',  &
         avgflag='A', long_name='MML net flux of shortwave at surface (in - out), pos into land', &
         ptr_lnd=this%mml_lnd_fsns_grc)
    
    this%mml_lnd_flns_grc(begg:endg) = spval
    call hist_addfld1d (fname='MML_flns', units='W/m2',  &
         avgflag='A', long_name='MML net flux of longwave at surface (out-in), pos out of land', &
         ptr_lnd=this%mml_lnd_flns_grc)
         
    this%mml_lnd_snowmelt(begg:endg) = spval
    call hist_addfld1d (fname='MML_snowmelt', units='kg/m2',  &
         avgflag='A', long_name='MML snow that melted into water bucket', &
         ptr_lnd=this%mml_lnd_snowmelt)
         
    ! Soil variables
! start 2d
  
    ! I wanted to add an mml case to the type2d, but for now change it back, since its crashing
 
 	!write(iulog,*)  'MML write to h0: 2d soil vars '
 	
 	this%mml_nc_dust_grc(begg:endg,:) = spval
    call hist_addfld2d (fname='MML_dust_2atm', units='unknown',  type2d='mml_dust', &
         avgflag='A', long_name='MML surface dust flux to atmosphere ', &
         ptr_lnd=this%mml_nc_dust_grc)
 	
    this%mml_nc_soil_maxice_grc(begg:endg,:) = spval
    call hist_addfld2d (fname='MML_maxice', units='kg/m3',  type2d='mml_lev', &
         avgflag='A', long_name='MML maximum freezable water in each soil layer; for thermal calculations', &
         ptr_lnd=this%mml_nc_soil_maxice_grc)     


    this%mml_nc_soil_levels_grc(begg:endg,:) = spval
    call hist_addfld2d (fname='MML_soilz', units='m',  type2d='mml_lev', &
         avgflag='A', long_name='MML depth (negative) from surface of midpoint of each soil layer', &
         ptr_lnd=this%mml_nc_soil_levels_grc, mml_dim=mml_nsoi)
    
    this%mml_soil_t_grc(begg:endg,:) = spval
    call hist_addfld2d (fname='MML_soil_t', units='K', type2d='mml_lev', &
         avgflag='A', long_name='MML soil temperature at each layer', &
         ptr_lnd=this%mml_soil_t_grc, mml_dim=mml_nsoi)
         
    
    this%mml_soil_liq_grc(begg:endg,:) = spval
    call hist_addfld2d (fname='MML_soil_liq', units='kg/m2',  type2d='mml_lev', &
         avgflag='A', long_name='MML kg of liquid water in each soil layer (Thermodynamic ONLY)', &
         ptr_lnd=this%mml_soil_liq_grc, mml_dim=mml_nsoi)
    
    this%mml_soil_ice_grc(begg:endg,:) = spval
    call hist_addfld2d (fname='MML_soil_ice', units='kg/m2', type2d='mml_lev',  &
         avgflag='A', long_name='MML kg of frozen water in each soil layer (Thermodynamic ONLY)', &
         ptr_lnd=this%mml_soil_ice_grc, mml_dim=mml_nsoi)
    
    this%mml_soil_dz_grc(begg:endg,:) = spval
    call hist_addfld2d (fname='MML_dz', units='m', type2d='mml_lev',  &
         avgflag='A', long_name='MML thickness of each soil layer', &
         ptr_lnd=this%mml_soil_dz_grc, mml_dim=mml_nsoi)
    
    this%mml_soil_zh_grc(begg:endg,:) = spval
    call hist_addfld2d (fname='MML_zh', units='m', type2d='mml_lev',  &
         avgflag='A', long_name='MML soil depth at interface between each soil layer', &
         ptr_lnd=this%mml_soil_zh_grc, mml_dim=mml_nsoi)
    
    this%mml_soil_tk_grc(begg:endg,:) = spval
    call hist_addfld2d (fname='MML_tk', units='W/m/K', type2d='mml_lev',  &
         avgflag='A', long_name='MML thermal conductivity of each soil layer', &
         ptr_lnd=this%mml_soil_tk_grc, mml_dim=mml_nsoi)
    
    this%mml_soil_tk_1d_grc(begg:endg) = spval
    call hist_addfld1d (fname='MML_tk_1d', units='W/m/K',  &
         avgflag='A', long_name='MML thermal resistance of every soil layer', &
         ptr_lnd=this%mml_soil_tk_1d_grc)
    
    this%mml_soil_tkh_grc(begg:endg,:) = spval
    call hist_addfld2d (fname='MML_tkh', units='W/m/K', type2d='mml_lev',  &
         avgflag='A', long_name='MML thermal conductivity at bottom boundary of each soil layer', &
         ptr_lnd=this%mml_soil_tkh_grc, mml_dim=mml_nsoi)
    
    this%mml_soil_dtsoi_grc(begg:endg,:) = spval
    call hist_addfld2d (fname='MML_dtsoi', units='K', type2d='mml_lev',  &
         avgflag='A', long_name='MML temperature tendency in each soil layer', &
         ptr_lnd=this%mml_soil_dtsoi_grc, mml_dim=mml_nsoi)
    
    this%mml_soil_cv_grc(begg:endg,:) = spval
    call hist_addfld2d (fname='MML_cv', units='J/m3/K', type2d='mml_lev',  &
         avgflag='A', long_name='MML heat capacity of each soil layer (depends on soil type)', &
         ptr_lnd=this%mml_soil_cv_grc, mml_dim=mml_nsoi)
	
	this%mml_soil_cv_1d_grc(begg:endg) = spval
    call hist_addfld1d (fname='MML_cv_1d', units='J/m3/K',  &
         avgflag='A', long_name='MML heat capacity of every soil layer', &
         ptr_lnd=this%mml_soil_cv_1d_grc)
	
	this%mml_glc_tk_1d_grc(begg:endg) = spval
    call hist_addfld1d (fname='MML_glc_tk_1d', units='W/m/K',  &
         avgflag='A', long_name='MML thermal resistance of every ice layer where glaciated', &
         ptr_lnd=this%mml_glc_tk_1d_grc)
    
    this%mml_glc_cv_1d_grc(begg:endg) = spval
    call hist_addfld1d (fname='MML_glc_cv_1d', units='J/m3/K',  &
         avgflag='A', long_name='MML heat capacity of every ice layer where glaciated', &
         ptr_lnd=this%mml_glc_cv_1d_grc)
         
! end 2d
    
    !write(iulog,*)  'MML write to h0: 1d soil vars '
    
    this%mml_soil_water_grc(begg:endg) = spval
    call hist_addfld1d (fname='MML_water', units='kg/m2',   &
         avgflag='A', long_name='MML total amount of liquid water in soil bucket (hydrology)', &
         ptr_lnd=this%mml_soil_water_grc)
    
    this%mml_soil_snow_grc(begg:endg) = spval
    call hist_addfld1d (fname='MML_snow', units='kg/m2',  &
         avgflag='A', long_name='MML total amount of snow in snow bucket (hydrology)', &
         ptr_lnd=this%mml_soil_snow_grc)
    
    this%mml_soil_runoff_grc(begg:endg) = spval
    call hist_addfld1d (fname='MML_runoff', units='kg/m2',  &
         avgflag='A', long_name='MML water in excess of bucket capacity (runoff, but it disappears)', &
         ptr_lnd=this%mml_soil_runoff_grc)
    
    ! lnd2atm MML vars
    this%mml_out_tref2m_grc(begg:endg) = spval
    call hist_addfld1d (fname='MML_l2a_tref2m', units='K',  &
         avgflag='A', long_name='MML 2m ref height temperature calculated from tsrf and tstar', &
         ptr_lnd=this%mml_out_tref2m_grc)
    
    this%mml_out_qref2m_grc(begg:endg) = spval
    call hist_addfld1d (fname='MML_l2a_qref2m', units='kg/kg',  &
         avgflag='A', long_name='MML 2m ref height humidity calculated from qsrf and qstar', &
         ptr_lnd=this%mml_out_qref2m_grc)
    
    this%mml_out_uref10m_grc(begg:endg) = spval
    call hist_addfld1d (fname='MML_l2a_uref10m', units='m/s',  &
         avgflag='A', long_name='MML 10m ref wind calculated from ustar', &
         ptr_lnd=this%mml_out_uref10m_grc)
    
    this%mml_out_taux(begg:endg) = spval
    call hist_addfld1d (fname='MML_l2a_taux', units='m/s',  &
         avgflag='A', long_name='MML zonal surface stress	', &
         ptr_lnd=this%mml_out_taux)
    
    this%mml_out_tauy(begg:endg) = spval
    call hist_addfld1d (fname='MML_l2a_tauy', units='m/s',  &
         avgflag='A', long_name='MML meridional surface stress	', &
         ptr_lnd=this%mml_out_tauy)
    
    ! MML check if latent heat flux is larger than atm can support (giant dew)
    this%mml_q_excess(begg:endg) = spval
    call hist_addfld1d (fname='MML_q_excess', units='kg/m2/s',  &
         avgflag='A', long_name='MML over-demand of dew (positive downwards) by land from atmosphere', &
         ptr_lnd=this%mml_q_excess)
    
    this%mml_lh_excess(begg:endg) = spval
    call hist_addfld1d (fname='MML_lh_excess', units='W/m2',  &
         avgflag='A', long_name='MML over-demand of latent heat flux (dew; positive downwards) by land from atmosphere', &
         ptr_lnd=this%mml_lh_excess)
    
    this%mml_q_demand(begg:endg) = spval
    call hist_addfld1d (fname='MML_q_demand', units='kg/m2/s',  &
         avgflag='A', long_name='MML initial demand of water flux by land from atmosphere (before correction for excess dew)', &
         ptr_lnd=this%mml_q_demand)
    
    this%mml_lh_demand(begg:endg) = spval
    call hist_addfld1d (fname='MML_lh_demand', units='W/m2',  &
         avgflag='A', long_name='MML initial demand of latent heat flux by land from atmosphere (before correction for excess dew)', &
         ptr_lnd=this%mml_lh_demand)
         
         
    ! mml diagnostic vars (temproary)
    
    this%mml_diag1_1d_grc(begg:endg) = spval
    call hist_addfld1d (fname='MML_diag1_1d', units='n/a',  &
         avgflag='A', long_name='MML temporary 1d diagnostic var 1', &
         ptr_lnd=this%mml_diag1_1d_grc)
    
    this%mml_diag2_1d_grc(begg:endg) = spval
    call hist_addfld1d (fname='MML_diag2_1d', units='n/a',  &
         avgflag='A', long_name='MML temporary 1d diagnostic var 2', &
         ptr_lnd=this%mml_diag2_1d_grc)
    
    this%mml_diag3_1d_grc(begg:endg) = spval
    call hist_addfld1d (fname='MML_diag3_1d', units='n/a',  &
         avgflag='A', long_name='MML temporary 1d diagnostic var 3', &
         ptr_lnd=this%mml_diag3_1d_grc)
    
    this%mml_diag1_2d_grc(begg:endg,:) = spval
    call hist_addfld2d (fname='MML_diag1_2d', units='n/a', type2d='mml_lev',  &
         avgflag='A', long_name='MML temporary 2d diagnostic var 1', &
         ptr_lnd=this%mml_diag1_2d_grc, mml_dim=mml_nsoi)
    
    this%mml_diag2_2d_grc(begg:endg,:) = spval
    call hist_addfld2d (fname='MML_diag2_2d', units='n/a', type2d='mml_lev',  &
         avgflag='A', long_name='MML temporary 2d diagnostic var 2', &
         ptr_lnd=this%mml_diag2_2d_grc, mml_dim=mml_nsoi)
    
    this%mml_diag3_2d_grc(begg:endg,:) = spval
    call hist_addfld2d (fname='MML_diag3_2d', units='n/a', type2d='mml_lev',  &
         avgflag='A', long_name='MML temporary 2d diagnostic var 3', &
         ptr_lnd=this%mml_diag3_2d_grc, mml_dim=mml_nsoi)
    
    
    ! mml error flux/balance vars
    
    this%mml_err_h2o(begg:endg) = spval
    call hist_addfld1d (fname='mml_err_h2o', units='n/a',  &
         avgflag='A', long_name='MML total water conservation error', &
         ptr_lnd=this%mml_err_h2o)
         
    this%mml_err_h2osno(begg:endg) = spval
    call hist_addfld1d (fname='mml_err_h2osno', units='n/a',  &
         avgflag='A', long_name='MML imbalance in snow depth (liquid water)', &
         ptr_lnd=this%mml_err_h2osno) 
    
    this%mml_err_seb(begg:endg) = spval
    call hist_addfld1d (fname='mml_err_seb', units='n/a',  &
         avgflag='A', long_name='MML surface energy conservation error', &
         ptr_lnd=this%mml_err_seb)
         
    this%mml_err_soi(begg:endg) = spval
    call hist_addfld1d (fname='mml_err_soi', units='n/a',  &
         avgflag='A', long_name='MML soil/lake energy conservation error', &
         ptr_lnd=this%mml_err_soi)
         
    this%mml_err_sol(begg:endg) = spval
    call hist_addfld1d (fname='mml_err_sol', units='n/a',  &
         avgflag='A', long_name='MML solar radiation conservation error', &
         ptr_lnd=this%mml_err_sol)
         
                           
    !write(iulog,*)  'MML write to h0: end of my vars '
    
                     
	! End MML simple land model added variables
	!-----------------------------------------------------------------------


    ! Time averaged quantities
    this%fsi24_patch(begp:endp) = spval
    call hist_addfld1d (fname='FSI24', units='K',  &
         avgflag='A', long_name='indirect radiation (last 24hrs)', &
         ptr_patch=this%fsi24_patch, default='inactive')

    this%fsi240_patch(begp:endp) = spval
    call hist_addfld1d (fname='FSI240', units='K',  &
         avgflag='A', long_name='indirect radiation (last 240hrs)', &
         ptr_patch=this%fsi240_patch, default='inactive')

    this%fsd24_patch(begp:endp) = spval
    call hist_addfld1d (fname='FSD24', units='K',  &
         avgflag='A', long_name='direct radiation (last 24hrs)', &
         ptr_patch=this%fsd24_patch, default='inactive')

    this%fsd240_patch(begp:endp) = spval
    call hist_addfld1d (fname='FSD240', units='K',  &
         avgflag='A', long_name='direct radiation (last 240hrs)', &
         ptr_patch=this%fsd240_patch, default='inactive')

  end subroutine InitHistory

!-----------------------------------------------------------------------
!	MML: 2016.01.15 Adding InitCold subroutine
!-----------------------------------------------------------------------
  subroutine InitCold(this, bounds)
  !MML: how do I keep this from running every time step? 
  
  ! !DESCRIPTION:
  ! Initialize time constant variables (like bucket size?) and cold start conditions 
  ! (like bucket water amount)
  !
  ! !USES:
  
  ! !ARGUMENTS:
   class(atm2lnd_type) :: this
   type(bounds_type), intent(in) :: bounds  
  !
  !
  ! !LOCAL VARIABLES:
  integer 		:: g   	! loop over grid indices
  integer  		:: begg, endg
  real(r8)		:: max_freezable, tfreeze		! bucket capacity (constant), 
  integer		:: i 	! index
  integer		:: mml_nsoi
  ! 												! water in bucket, snow in bucket

  real(r8) 		:: ival

  begg = bounds%begg
  endg = bounds%endg
  
  mml_nsoi = 10
  
  ival = 1e22	! some stupidly large value to initialize fields that are going to get math done to them.
  				! if a number this big shows up... it means somewhere in the calculation of hte variable we're missing smoething
  tfreeze = 273.15_r8
  max_freezable = 300._r8 ! [kg/m3] !if I make this vary, I'm going to have to change how I initialize this!!!
  
  ! Just initialize variables that it would need from the previous time step, not 
  ! things that would be read in from the .nc file. 
  
  this%mml_soil_water_grc(:) 	= 150._r8  ! 150.0		! start bucket (typical capacity 200) with 150 kg/m2 = 150 mm of water in it.
  this%mml_soil_snow_grc(:)		= 0._r8
  this%mml_soil_runoff_grc(:)	= 0._r8
  
  this%mml_soil_t_grc(:,:)		= tfreeze+2._r8
  
  ! need this to be in kg/m2, not kg/m3 => must multiply by dz ... but theoretically don't know dz yet... 
  ! So I'll hard code it in to here, but thats not a very robust way of doing things (just want to 
  ! try and fix the soil temperature energy conservation error for now)
  do i = 1, mml_nsoi	
		this%mml_nc_soil_levels_grc(:,i) =  -0.025 * (exp(0.5*(i-0.5)) - 1.0);
  enddo
  
  this%mml_soil_dz_grc(:,1) = 0.5_r8 * abs( this%mml_nc_soil_levels_grc(:,1) + this%mml_nc_soil_levels_grc(:,2) )
  do i = 2, mml_nsoi-1
     this%mml_soil_dz_grc(:,i) =  0.5_r8 * ( this%mml_nc_soil_levels_grc(:,i-1) - this%mml_nc_soil_levels_grc(:,i+1) )
  enddo
  this%mml_soil_dz_grc(:,mml_nsoi) = this%mml_nc_soil_levels_grc(:,mml_nsoi-1) - this%mml_nc_soil_levels_grc(:,mml_nsoi)
  
  this%mml_nc_soil_maxice_grc(:,1) = 0._r8 * this%mml_soil_dz_grc(:,1)
  this%mml_soil_zh_grc(:,1) = -this%mml_soil_dz_grc(:,1)
  do i = 2, mml_nsoi
		this%mml_nc_soil_maxice_grc(:,i) = 300._r8 * this%mml_soil_dz_grc(:,i)
		this%mml_soil_zh_grc(:,i) = this%mml_soil_zh_grc(:,i-1) - this%mml_soil_dz_grc(:,i) ! aha, zh, not z
  end do
  
  
  this%mml_soil_liq_grc(:,:)	= this%mml_nc_soil_maxice_grc	
 ! this%mml_soil_liq_grc(:,1) 	= 0._r8
  this%mml_soil_ice_grc(:,:)	= 0._r8
  
  
  this%mml_lnd_ts_grc(:)		= tfreeze+2._r8
  this%mml_lnd_qs_grc(:)		= 0.02_r8
  
  ! or should I just cold-start ALL the MML variables to zero? 
     
	! init-cold nc vars with something *reasonable* in case they're not read from a .nc file right away?
	this%mml_nc_alb_gvd_grc(:)   = 0.23_r8
	this%mml_nc_alb_svd_grc(:)   = 0.23_r8
	this%mml_nc_alb_gnd_grc(:)  = 0.23_r8
	this%mml_nc_alb_snd_grc(:)  = 0.23_r8
	this%mml_nc_alb_gvf_grc(:)  = 0.23_r8
	this%mml_nc_alb_svf_grc(:)  = 0.23_r8
	this%mml_nc_alb_gnf_grc(:)  = 0.23_r8
	this%mml_nc_alb_snf_grc(:)   = 0.23_r8
	this%mml_nc_snowmask_grc(:)   = 100._r8
	this%mml_nc_evaprs_grc(:)   = 20._r8
	this%mml_nc_bucket_cap_grc(:)   = 200._r8
	!this%mml_nc_soil_maxice_grc(:,:)   
	!this%mml_nc_soil_levels_grc(:,:)  
	!this%mml_nc_soil_type_grc(:)	
																		  ! for simplicity, set all layers and all grid cells equal. For flexibility, program this way. 
	this%mml_nc_roughness_grc(:)   = 10._r8
	this%mml_nc_emiss_grc(:)  = 1._r8
	this%mml_nc_glcmask_grc(:)   = 0.0_r8	! no glc mask
	this%mml_nc_dust_grc(:,:)		= 0.0_r8	! no dust
	
	! **********************************
	! TEMPORARY FOR DEVELOPEMENT ONLY
	! **********************************
	! 
	! Init Cold-ing all my vars that aren't either grabbed from the .nc file or 
	! copied from the atmosphere, because starting them at ival in initAllocate doesn't
	! seem to be cutting it (getting weird spatially memory problems that I don't really 
	! understand, but doing math with those vars without giving them a new value first 
	! seems to cause problems...
	!
	! (except vars I initialized on purpose above :p )
	
	! land flux data
	
	! these 2 are commented because we already initialize them with somethign sensible above
	!this%mml_lnd_ts_grc     (:)   = ival
	!this%mml_lnd_qs_grc     (:)   = ival
	
	this%mml_lnd_qa_grc     (:)   = ival
	this%mml_lnd_swabs_grc  (:)   = ival
	this%mml_lnd_fsr_grc   (:)   = ival
	this%mml_lnd_lwup_grc   (:)   = ival
	this%mml_lnd_shflx_grc  (:)   = ival
	this%mml_lnd_lhflx_grc  (:)   = ival
	this%mml_lnd_gsoi_grc   (:)   = ival
	this%mml_lnd_gsnow_grc  (:)   = ival
	this%mml_lnd_evap_grc   (:)   = ival
	this%mml_lnd_ustar_grc  (:)   = ival
	this%mml_lnd_tstar_grc  (:)   = ival
	this%mml_lnd_qstar_grc  (:)   = ival
	this%mml_lnd_tvstar_grc  (:)   = ival
	this%mml_lnd_obu_grc    (:)   = ival
	this%mml_lnd_ram_grc    (:)   = ival
	this%mml_lnd_rah_grc    (:)   = ival
	this%mml_lnd_res_grc    (:)   = ival
	this%mml_lnd_effective_res_grc    (:)   = ival
	this%mml_lnd_beta_grc    (:)   = ival
	this%mml_lnd_disp_grc   (:)   = ival
	this%mml_lnd_z0m_grc    (:)   = ival
	this%mml_lnd_z0h_grc    (:)   = ival
	this%mml_lnd_alb_grc    (:)   = ival
	this%mml_lnd_fsns_grc    (:)   = ival
	this%mml_lnd_flns_grc    (:)   = ival
	this%mml_lnd_snowmelt    (:)   = 0.0_r8
	
	! soil data 
	this%mml_soil_tk_grc     (:,:)   = ival 
	this%mml_soil_tk_1d_grc     (:)   = ival 
	this%mml_soil_tkh_grc     (:,:)   = ival 
	this%mml_soil_dtsoi_grc     (:,:)   = ival 
	this%mml_soil_cv_grc     (:,:)   = ival 
	this%mml_soil_cv_1d_grc     (:)   = ival 
	this%mml_glc_tk_1d_grc     (:)   = ival
	this%mml_glc_cv_1d_grc     (:)   = ival
	! these were commented:
!	this%mml_soil_t_grc      (:,:)   = ival 
! 	this%mml_soil_dz_grc     (:,:)   = ival 
! 	this%mml_soil_zh_grc     (:,:)   = ival
! 	this%mml_soil_water_grc  (:)   	 = ival
! 	this%mml_soil_snow_grc   (:)   	 = ival
! 	this%mml_soil_runoff_grc (:)   	 = ival
  	
  	this%mml_out_tref2m_grc	 (:)		= ival
	this%mml_out_qref2m_grc	 (:)		= ival
	this%mml_out_uref10m_grc (:)		= ival
	this%mml_out_taux (:)		= ival
	this%mml_out_tauy (:)		= ival
	
	! 300 kg/m3, (0 in top layer) but need in kg/m2 units, ie x dz(:,i)
	! so for now I'll manually prescribe the initial values, but this is 
	! NOT a robust way of doing things!!!
	! we'll put it all in liquid, for now
	! (was commented... )
	this%mml_soil_ice_grc    (:,:)   = 0.0 
	
	this%mml_soil_liq_grc    (:,1)   = 0.0 
	this%mml_soil_liq_grc    (:,2)   = 0.082736907779029 
	this%mml_soil_liq_grc    (:,3)   = 0.136410099727240
	this%mml_soil_liq_grc    (:,4)   = 0.224902232958626
	this%mml_soil_liq_grc    (:,5)   = 0.370801095306842
	this%mml_soil_liq_grc    (:,6)   = 0.611347653031295
	this%mml_soil_liq_grc    (:,7)   = 1.007941879345298 
	this%mml_soil_liq_grc    (:,8)   = 1.661815216106055
	this%mml_soil_liq_grc    (:,9)   = 2.739870094767183
	this%mml_soil_liq_grc    (:,10)  = 3.410915413537486

  
  end subroutine InitCold

! MML: InitAccBuffer sounds like what I actually want... unless it only initializes a 
! value with the necessity of having a r0 file overwrite it. 
  !-----------------------------------------------------------------------
  subroutine InitAccBuffer (this, bounds)
    !
    ! !DESCRIPTION:
    ! Initialize accumulation buffer for all required module accumulated fields
    ! This routine set defaults values that are then overwritten by the
    ! restart file for restart or branch runs
    !
    ! !USES 
    use clm_varcon  , only : spval
    use accumulMod  , only : init_accum_field
    !
    ! !ARGUMENTS:
    class(atm2lnd_type) :: this
    type(bounds_type), intent(in) :: bounds  
    !---------------------------------------------------------------------

    this%fsd24_patch(bounds%begp:bounds%endp) = spval
    call init_accum_field (name='FSD24', units='W/m2',                                             &
         desc='24hr average of direct solar radiation',  accum_type='runmean', accum_period=-1,    &
         subgrid_type='pft', numlev=1, init_value=0._r8)

    this%fsd240_patch(bounds%begp:bounds%endp) = spval
    call init_accum_field (name='FSD240', units='W/m2',                                            &
         desc='240hr average of direct solar radiation',  accum_type='runmean', accum_period=-10,  &
         subgrid_type='pft', numlev=1, init_value=0._r8)

    this%fsi24_patch(bounds%begp:bounds%endp) = spval
    call init_accum_field (name='FSI24', units='W/m2',                                             &
         desc='24hr average of diffuse solar radiation',  accum_type='runmean', accum_period=-1,   &
         subgrid_type='pft', numlev=1, init_value=0._r8)

    this%fsi240_patch(bounds%begp:bounds%endp) = spval
    call init_accum_field (name='FSI240', units='W/m2',                                            &
         desc='240hr average of diffuse solar radiation',  accum_type='runmean', accum_period=-10, &
         subgrid_type='pft', numlev=1, init_value=0._r8)

  end subroutine InitAccBuffer

  !-----------------------------------------------------------------------
  subroutine InitAccVars(this, bounds)
    !
    ! !DESCRIPTION:
    ! Initialize module variables that are associated with
    ! time accumulated fields. This routine is called for both an initial run
    ! and a restart run (and must therefore must be called after the restart file 
    ! is read in and the accumulation buffer is obtained)
    !
    ! !USES 
    use accumulMod       , only : extract_accum_field
    use clm_time_manager , only : get_nstep
    !
    ! !ARGUMENTS:
    class(atm2lnd_type) :: this
    type(bounds_type), intent(in) :: bounds  
    !
    ! !LOCAL VARIABLES:
    integer  :: begp, endp
    integer  :: begc, endc
    integer  :: nstep
    integer  :: ier
    real(r8), pointer :: rbufslp(:)  ! temporary
    real(r8), pointer :: rbufslc(:)  ! temporary
    !---------------------------------------------------------------------

    begp = bounds%begp; endp = bounds%endp
    begc = bounds%begc; endc = bounds%endc

    ! Allocate needed dynamic memory for single level patch field
    allocate(rbufslp(begp:endp), stat=ier)
    if (ier/=0) then
       write(iulog,*)' in '
       call endrun(msg="InitAccVars allocation error for rbufslp"//&
            errMsg(sourcefile, __LINE__))
    endif
    ! Allocate needed dynamic memory for single level col field
    allocate(rbufslc(begc:endc), stat=ier)
    if (ier/=0) then
       write(iulog,*)' in '
       call endrun(msg="InitAccVars allocation error for rbufslc"//&
            errMsg(sourcefile, __LINE__))
    endif

    ! Determine time step
    nstep = get_nstep()

    call extract_accum_field ('FSD24', rbufslp, nstep)
    this%fsd24_patch(begp:endp) = rbufslp(begp:endp)

    call extract_accum_field ('FSD240', rbufslp, nstep)
    this%fsd240_patch(begp:endp) = rbufslp(begp:endp)

    call extract_accum_field ('FSI24', rbufslp, nstep)
    this%fsi24_patch(begp:endp) = rbufslp(begp:endp)

    call extract_accum_field ('FSI240', rbufslp, nstep)
    this%fsi240_patch(begp:endp) = rbufslp(begp:endp)

    deallocate(rbufslp)
    deallocate(rbufslc)

  end subroutine InitAccVars

  !-----------------------------------------------------------------------
  subroutine UpdateAccVars (this, bounds)
    !
    ! USES
    use clm_time_manager, only : get_nstep
    use accumulMod      , only : update_accum_field, extract_accum_field
    !
    ! !ARGUMENTS:
    class(atm2lnd_type)                 :: this
    type(bounds_type)      , intent(in) :: bounds  
    !
    ! !LOCAL VARIABLES:
    integer :: g,c,p                     ! indices
    integer :: dtime                     ! timestep size [seconds]
    integer :: nstep                     ! timestep number
    integer :: ier                       ! error status
    integer :: begp, endp
    integer :: begc, endc
    real(r8), pointer :: rbufslp(:)      ! temporary single level - patch level
    real(r8), pointer :: rbufslc(:)      ! temporary single level - column level
    !---------------------------------------------------------------------

    begp = bounds%begp; endp = bounds%endp
    begc = bounds%begc; endc = bounds%endc

    nstep = get_nstep()

    ! Allocate needed dynamic memory for single level patch field
    allocate(rbufslp(begp:endp), stat=ier)
    if (ier/=0) then
       write(iulog,*)'UpdateAccVars allocation error for rbufslp'
       call endrun(msg=errMsg(sourcefile, __LINE__))
    endif
    ! Allocate needed dynamic memory for single level col field
    allocate(rbufslc(begc:endc), stat=ier)
    if (ier/=0) then
       write(iulog,*)'UpdateAccVars allocation error for rbufslc'
       call endrun(msg=errMsg(sourcefile, __LINE__))
    endif

    ! Accumulate and extract forc_solad24 & forc_solad240 
    do p = begp,endp
       g = patch%gridcell(p)
       rbufslp(p) = this%forc_solad_grc(g,1)
    end do
    call update_accum_field  ('FSD240', rbufslp               , nstep)
    call extract_accum_field ('FSD240', this%fsd240_patch     , nstep)
    call update_accum_field  ('FSD24' , rbufslp               , nstep)
    call extract_accum_field ('FSD24' , this%fsd24_patch      , nstep)

    ! Accumulate and extract forc_solai24 & forc_solai240 
    do p = begp,endp
       g = patch%gridcell(p)
       rbufslp(p) = this%forc_solai_grc(g,1)
    end do
    call update_accum_field  ('FSI24' , rbufslp               , nstep)
    call extract_accum_field ('FSI24' , this%fsi24_patch      , nstep)
    call update_accum_field  ('FSI240', rbufslp               , nstep)
    call extract_accum_field ('FSI240', this%fsi240_patch     , nstep)

    ! Precipitation accumulators

    do p = begp,endp
       c = patch%column(p)
       rbufslp(p) = this%forc_rain_downscaled_col(c) + this%forc_snow_downscaled_col(c)
       rbufslc(c) = this%forc_rain_downscaled_col(c) + this%forc_snow_downscaled_col(c)
    end do

    deallocate(rbufslp)
    deallocate(rbufslc)

  end subroutine UpdateAccVars

  !------------------------------------------------------------------------
  subroutine Restart(this, bounds, ncid, flag)
    ! 
    ! !USES:
    use restUtilMod
    use ncdio_pio
    !
    ! !ARGUMENTS:
    class(atm2lnd_type) :: this
    type(bounds_type), intent(in) :: bounds  
    type(file_desc_t), intent(inout) :: ncid   
    character(len=*) , intent(in)    :: flag   
    !
    ! !LOCAL VARIABLES:
    logical            :: readvar 
    !------------------------------------------------------------------------

    call restartvar(ncid=ncid, flag=flag, varname='qflx_floodg', xtype=ncd_double, &
         dim1name='gridcell', &
         long_name='flood water flux', units='mm/s', &
         interpinic_flag='skip', readvar=readvar, data=this%forc_flood_grc)
    if (flag == 'read' .and. .not. readvar) then
       ! initial run, readvar=readvar, not restart: initialize flood to zero
       this%forc_flood_grc = 0._r8
    endif

    ! -----------------------------------------------------------------------
    ! Start MML simple land model restart variables section		! MML 2016.01.15
    
     write(iulog,*)  ' MML trying to write r 1d restart vars '
    
    ! MML: surface
    call restartvar(ncid=ncid, flag=flag, varname='mml_lnd_ts_grc', xtype=ncd_double, &
         dim1name='gridcell',  &
         long_name='Surface Temperature for MO', units='K', &
         interpinic_flag='skip', readvar=readvar, data=this%mml_lnd_ts_grc) 
         
    call restartvar(ncid=ncid, flag=flag, varname='mml_lnd_qs_grc', xtype=ncd_double, &
         dim1name='gridcell',  &
         long_name='surface specific humidity for MO', units='kg/kg', &
         interpinic_flag='skip', readvar=readvar, data=this%mml_lnd_qs_grc) 
    
  	! MML soil:
  	! MML Hydrology variables:
    call restartvar(ncid=ncid, flag=flag, varname='mml_soil_water_grc', xtype=ncd_double, &
         dim1name='gridcell', &
         long_name='soil bucket water content', units='kg', &
         interpinic_flag='skip', readvar=readvar, data=this%mml_soil_water_grc)
         
    call restartvar(ncid=ncid, flag=flag, varname='mml_soil_snow_grc', xtype=ncd_double, &
         dim1name='gridcell', &
         long_name='snow bucket snow content', units='kg', &
         interpinic_flag='skip', readvar=readvar, data=this%mml_soil_snow_grc)
         
    call restartvar(ncid=ncid, flag=flag, varname='mml_soil_runoff_grc', xtype=ncd_double, &
         dim1name='gridcell', &
         long_name='water runoff', units='kg', &
         interpinic_flag='skip', readvar=readvar, data=this%mml_soil_runoff_grc)
    
   !  write(iulog,*)  'MML trying to write r 2d restart vars '
    ! MML Thermodynamic vars for each soil level (3d)
    call restartvar(ncid=ncid, flag=flag, varname='mml_soil_liq_grc', xtype=ncd_double, &
         dim1name='gridcell',  dim2name='mml_lev', switchdim=.true., &  ! dim2 mml_lev?
         long_name='amount of liquid water in each soil layer', units='kg', &
         interpinic_flag='skip', readvar=readvar, data=this%mml_soil_liq_grc)
         
    call restartvar(ncid=ncid, flag=flag, varname='mml_soil_ice_grc', xtype=ncd_double, &
         dim1name='gridcell', dim2name='mml_lev', switchdim=.true., &
         long_name='amount of frozen water in each soil layer', units='kg', &
         interpinic_flag='skip', readvar=readvar, data=this%mml_soil_ice_grc)          
   
    call restartvar(ncid=ncid, flag=flag, varname='mml_soil_t_grc', xtype=ncd_double, &
         dim1name='gridcell', dim2name='mml_lev', switchdim=.true., &
         long_name='MML soil temperature at each layer', units='K', &
         interpinic_flag='skip', readvar=readvar, data=this%mml_soil_t_grc)   

    call restartvar(ncid=ncid, flag=flag, varname='mml_soil_dtsoi_grc', xtype=ncd_double, &
         dim1name='gridcell', dim2name='mml_lev', switchdim=.true., &
         long_name='MML temperature tendency in each soil layer', units='K', &
         interpinic_flag='skip', readvar=readvar, data=this%mml_soil_dtsoi_grc) 

 
    ! MML nc vars, so if I stop mid-month / mid-day I can still know what that month's nc params are
    call restartvar(ncid=ncid, flag=flag, varname='mml_nc_alb_gvd_grc', xtype=ncd_double, &
         dim1name='gridcell',  &
         long_name='Ground visible direct albedo (from netcdf file)', units='none', &
         interpinic_flag='skip', readvar=readvar, data=this%mml_nc_alb_gvd_grc) 
         
     call restartvar(ncid=ncid, flag=flag, varname='mml_nc_alb_svd_grc', xtype=ncd_double, &
         dim1name='gridcell',  &
         long_name='Snow visible direct albedo (from netcdf file)', units='none', &
         interpinic_flag='skip', readvar=readvar, data=this%mml_nc_alb_svd_grc) 
       
    call restartvar(ncid=ncid, flag=flag, varname='mml_nc_alb_gnd_grc', xtype=ncd_double, &
         dim1name='gridcell',  &
         long_name='Ground NIR direct albedo (from netcdf file)', units='none', &
         interpinic_flag='skip', readvar=readvar, data=this%mml_nc_alb_gnd_grc) 
       
    call restartvar(ncid=ncid, flag=flag, varname='mml_nc_alb_snd_grc', xtype=ncd_double, &
         dim1name='gridcell',  &
         long_name='Snow NIR direct albedo (from netcdf file)', units='none', &
         interpinic_flag='skip', readvar=readvar, data=this%mml_nc_alb_snd_grc) 
     
    call restartvar(ncid=ncid, flag=flag, varname='mml_nc_alb_gvf_grc', xtype=ncd_double, &
         dim1name='gridcell',  &
         long_name='Ground visible diffuse albedo (from netcdf file)', units='none', &
         interpinic_flag='skip', readvar=readvar, data=this%mml_nc_alb_gvf_grc) 
    
    call restartvar(ncid=ncid, flag=flag, varname='mml_nc_alb_svf_grc', xtype=ncd_double, &
         dim1name='gridcell',  &
         long_name='snow visible diffuse albedo (from netcdf file)', units='none', &
         interpinic_flag='skip', readvar=readvar, data=this%mml_nc_alb_svf_grc) 
         
    call restartvar(ncid=ncid, flag=flag, varname='mml_nc_alb_gnf_grc', xtype=ncd_double, &
         dim1name='gridcell',  &
         long_name='Ground NIR diffuse albedo (from netcdf file)', units='K', &
         interpinic_flag='skip', readvar=readvar, data=this%mml_nc_alb_gnf_grc) 
         
    call restartvar(ncid=ncid, flag=flag, varname='mml_nc_alb_snf_grc', xtype=ncd_double, &
         dim1name='gridcell',  &
         long_name='Snow NIR diffuse albedo (from netcdf file)', units='none', &
         interpinic_flag='skip', readvar=readvar, data=this%mml_nc_alb_snf_grc) 
    
    call restartvar(ncid=ncid, flag=flag, varname='mml_nc_snowmask_grc', xtype=ncd_double, &
         dim1name='gridcell',  &
         long_name='Amount of snow required to fully mask ground albedo (from netcdf file)', units='kg/m2', &
         interpinic_flag='skip', readvar=readvar, data=this%mml_nc_snowmask_grc) 
    
    call restartvar(ncid=ncid, flag=flag, varname='mml_nc_evaprs_grc', xtype=ncd_double, &
         dim1name='gridcell',  &
         long_name='Evaporative resistance (from netcdf file)', units='s/m', &
         interpinic_flag='skip', readvar=readvar, data=this%mml_nc_evaprs_grc) 
    
    call restartvar(ncid=ncid, flag=flag, varname='mml_nc_bucket_cap_grc', xtype=ncd_double, &
         dim1name='gridcell',  &
         long_name='Bucket Capacity (from netcdf file)', units='kg/m2', &
         interpinic_flag='skip', readvar=readvar, data=this%mml_nc_bucket_cap_grc) 
    
    call restartvar(ncid=ncid, flag=flag, varname='mml_nc_roughness_grc', xtype=ncd_double, &
         dim1name='gridcell',  &
         long_name='Surface roughness (vegetation height) (from netcdf file)', units='m', &
         interpinic_flag='skip', readvar=readvar, data=this%mml_nc_roughness_grc) 
    
    call restartvar(ncid=ncid, flag=flag, varname='mml_nc_emiss_grc', xtype=ncd_double, &
         dim1name='gridcell',  &
         long_name='Surface emissivity (from netcdf file)', units='none', &
         interpinic_flag='skip', readvar=readvar, data=this%mml_nc_emiss_grc) 
    
    call restartvar(ncid=ncid, flag=flag, varname='mml_nc_glcmask_grc', xtype=ncd_double, &
         dim1name='gridcell',  &
         long_name='Mask of glaciated points (from netcdf file)', units='none', &
         interpinic_flag='skip', readvar=readvar, data=this%mml_nc_glcmask_grc) 
    
    call restartvar(ncid=ncid, flag=flag, varname='mml_nc_dust_grc', xtype=ncd_double, &
         dim1name='gridcell',  &
         long_name='Dust flux to atm (from netcdf file)', units='unknown', &
         interpinic_flag='skip', readvar=readvar, data=this%mml_nc_dust_grc) 
    
    write(iulog,*)  ' MML end of 1d restart vars '
    
    ! 3d restart var example (soilbiogeochem carbon mod):
    !call restartvar(ncid=ncid, flag=flag, varname='col_ctrunc_vr', xtype=ncd_double,  &
    !           dim1name='column', dim2name='levgrnd', switchdim=.true., &
    !           long_name='',  units='', fill_value=spval, &
    !           interpinic_flag='interp', readvar=readvar, data=ptr2d)
                         
	! End MML simple land model added variables
	!-----------------------------------------------------------------------

  end subroutine Restart

  !-----------------------------------------------------------------------
  subroutine Clean(this)
    !
    ! !DESCRIPTION:
    ! Finalize this instance
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    class(atm2lnd_type), intent(inout) :: this
    !
    ! !LOCAL VARIABLES:

    character(len=*), parameter :: subname = 'Clean'
    !-----------------------------------------------------------------------

    ! atm->lnd
    deallocate(this%forc_u_grc)
    deallocate(this%forc_v_grc)
    deallocate(this%forc_wind_grc)
    deallocate(this%forc_rh_grc)
    deallocate(this%forc_hgt_grc)
    deallocate(this%forc_topo_grc)
    deallocate(this%forc_hgt_u_grc)
    deallocate(this%forc_hgt_t_grc)
    deallocate(this%forc_hgt_q_grc)
    deallocate(this%forc_vp_grc)
    deallocate(this%forc_psrf_grc)
    deallocate(this%forc_pco2_grc)
    deallocate(this%forc_solad_grc)
    deallocate(this%forc_solai_grc)
    deallocate(this%forc_solar_grc)
    deallocate(this%forc_pc13o2_grc)
    deallocate(this%forc_po2_grc)
    deallocate(this%forc_aer_grc)

    ! atm->lnd not downscaled
    deallocate(this%forc_t_not_downscaled_grc)
    deallocate(this%forc_q_not_downscaled_grc)
    deallocate(this%forc_pbot_not_downscaled_grc)
    deallocate(this%forc_th_not_downscaled_grc)
    deallocate(this%forc_rho_not_downscaled_grc)
    deallocate(this%forc_lwrad_not_downscaled_grc)
    deallocate(this%forc_rain_not_downscaled_grc)
    deallocate(this%forc_snow_not_downscaled_grc)
    
    ! atm->lnd downscaled
    deallocate(this%forc_t_downscaled_col)
    deallocate(this%forc_q_downscaled_col)
    deallocate(this%forc_pbot_downscaled_col)
    deallocate(this%forc_th_downscaled_col)
    deallocate(this%forc_rho_downscaled_col)
    deallocate(this%forc_lwrad_downscaled_col)
    deallocate(this%forc_rain_downscaled_col)
    deallocate(this%forc_snow_downscaled_col)

    ! rof->lnd
    deallocate(this%forc_flood_grc)
    deallocate(this%volr_grc)
    deallocate(this%volrmch_grc)

    ! anomaly forcing
    deallocate(this%bc_precip_grc)
    deallocate(this%af_precip_grc)
    deallocate(this%af_uwind_grc)
    deallocate(this%af_vwind_grc)
    deallocate(this%af_tbot_grc)
    deallocate(this%af_pbot_grc)
    deallocate(this%af_shum_grc)
    deallocate(this%af_swdn_grc)
    deallocate(this%af_lwdn_grc)

    deallocate(this%fsd24_patch)
    deallocate(this%fsd240_patch)
    deallocate(this%fsi24_patch)
    deallocate(this%fsi240_patch)
    deallocate(this%prec10_patch)
    deallocate(this%prec60_patch)
    deallocate(this%prec365_col)
    deallocate(this%t_mo_patch)
    deallocate(this%t_mo_min_patch)
    
    ! MML: deallocate mml vars:
    
    deallocate(this%mml_nc_alb_gvd_grc    	 )
	deallocate(this%mml_nc_alb_svd_grc    	 )
	deallocate(this%mml_nc_alb_gnd_grc    	 )
	deallocate(this%mml_nc_alb_snd_grc    	 )
	! alb diffuse
	deallocate(this%mml_nc_alb_gvf_grc    	)
	deallocate(this%mml_nc_alb_svf_grc    	) 
	deallocate(this%mml_nc_alb_gnf_grc    	)
	deallocate(this%mml_nc_alb_snf_grc    	)
	! the rest
	deallocate(this%mml_nc_snowmask_grc    	)
	deallocate(this%mml_nc_evaprs_grc      	)
	deallocate(this%mml_nc_bucket_cap_grc  	)
	deallocate(this%mml_nc_soil_maxice_grc  )
	deallocate(this%mml_nc_soil_levels_grc  )
	deallocate(this%mml_nc_soil_type_grc    )
	deallocate(this%mml_nc_roughness_grc    )
	deallocate(this%mml_nc_emiss_grc    )
	deallocate(this%mml_nc_glcmask_grc    )
	deallocate(this%mml_nc_dust_grc    )
	
	! from atm data
	deallocate(this%mml_atm_fsds_grc    )
	 deallocate(this%mml_atm_fsdsnd_grc    )
	 deallocate(this%mml_atm_fsdsvd_grc    )
	 deallocate(this%mml_atm_fsdsni_grc    )
	 deallocate(this%mml_atm_fsdsvi_grc    )
	deallocate(this%mml_atm_lwdn_grc    )
	deallocate(this%mml_atm_zref_grc    )
	deallocate(this%mml_atm_tbot_grc    )
	deallocate(this%mml_atm_thref_grc   )
	deallocate(this%mml_atm_qbot_grc    )
	deallocate(this%mml_atm_uref_grc    )
	deallocate(this%mml_atm_eref_grc    )
	deallocate(this%mml_atm_pbot_grc     )
	deallocate(this%mml_atm_psrf_grc     )
	deallocate(this%mml_atm_rhomol_grc   )
	deallocate(this%mml_atm_rhoair_grc   )
	deallocate(this%mml_atm_cp_grc    	 )
	deallocate(this%mml_atm_pco2    	 )
	deallocate(this%mml_atm_prec_liq_grc )
	deallocate(this%mml_atm_prec_frz_grc )
	
	! land flux data
	deallocate(this%mml_lnd_ts_grc     )
	deallocate(this%mml_lnd_qs_grc     )
	deallocate(this%mml_lnd_qa_grc     )
	deallocate(this%mml_lnd_swabs_grc  )
	deallocate(this%mml_lnd_fsr_grc   )
	 deallocate(this%mml_lnd_fsrnd_grc   )
	 deallocate(this%mml_lnd_fsrni_grc   )
	 deallocate(this%mml_lnd_fsrvd_grc   )
	 deallocate(this%mml_lnd_fsrvi_grc   )
	deallocate(this%mml_lnd_lwup_grc   )
	deallocate(this%mml_lnd_shflx_grc  )
	deallocate(this%mml_lnd_lhflx_grc  )
	deallocate(this%mml_lnd_gsoi_grc   )
	deallocate(this%mml_lnd_gsnow_grc  )
	deallocate(this%mml_lnd_evap_grc   )
	deallocate(this%mml_lnd_ustar_grc  )
	deallocate(this%mml_lnd_tstar_grc  )
	deallocate(this%mml_lnd_qstar_grc  )
	deallocate(this%mml_lnd_tvstar_grc  )
	deallocate(this%mml_lnd_obu_grc    )
	deallocate(this%mml_lnd_ram_grc    )
	deallocate(this%mml_lnd_rah_grc    )
	deallocate(this%mml_lnd_res_grc    )
	deallocate(this%mml_lnd_effective_res_grc    )
	deallocate(this%mml_lnd_beta_grc    )
	deallocate(this%mml_lnd_disp_grc   )
	deallocate(this%mml_lnd_z0m_grc    )
	deallocate(this%mml_lnd_z0h_grc    )
	deallocate(this%mml_lnd_alb_grc    )
	deallocate(this%mml_lnd_fsns_grc    )
	deallocate(this%mml_lnd_flns_grc    )
	deallocate(this%mml_lnd_snowmelt    )
	
	! soil data
	deallocate(this%mml_soil_t_grc     )
	deallocate(this%mml_soil_liq_grc   )
	deallocate(this%mml_soil_ice_grc   ) 
	deallocate(this%mml_soil_dz_grc    )
	deallocate(this%mml_soil_zh_grc    )
	deallocate(this%mml_soil_tk_grc    )
	deallocate(this%mml_soil_tk_1d_grc    )
	deallocate(this%mml_soil_tkh_grc    )
	deallocate(this%mml_soil_dtsoi_grc    )
	deallocate(this%mml_soil_cv_grc    )
	deallocate(this%mml_soil_cv_1d_grc    )
	deallocate(this%mml_soil_water_grc )
	deallocate(this%mml_soil_snow_grc  )
	deallocate(this%mml_soil_runoff_grc)
	deallocate(this%mml_glc_tk_1d_grc    )
	deallocate(this%mml_glc_cv_1d_grc    )
	
	! over-large dew
	deallocate(this%mml_lh_excess )
	deallocate(this%mml_q_excess )
	deallocate(this%mml_lh_demand )
	deallocate(this%mml_q_demand )
	
	! lnd2atm data
	deallocate(this%mml_out_tref2m_grc )
	deallocate(this%mml_out_qref2m_grc )
	deallocate(this%mml_out_uref10m_grc )
	deallocate(this%mml_out_taux )
	deallocate(this%mml_out_tauy )
	
	! diagnostic temporary vars
	deallocate(this%mml_diag1_1d_grc )
	deallocate(this%mml_diag2_1d_grc )
	deallocate(this%mml_diag3_1d_grc )
	
	deallocate(this%mml_diag1_2d_grc )
	deallocate(this%mml_diag2_2d_grc )
	deallocate(this%mml_diag3_2d_grc )
	
	! error tracking vars
	deallocate(this%mml_err_h2o )
	deallocate(this%mml_err_h2osno )
	deallocate(this%mml_err_seb )
	deallocate(this%mml_err_soi )
	deallocate(this%mml_err_sol )
	
  end subroutine Clean


end module atm2lndType
