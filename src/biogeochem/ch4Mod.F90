module ch4Mod

#include "shr_assert.h"

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Module holding routines to calculate methane fluxes
  ! The driver averages up to gridcell, weighting by finundated, and checks for balance errors.
  ! Sources, sinks, "competition" for CH4 & O2, & transport are resolved in ch4_tran.
  !
  ! !USES:
  use shr_kind_mod                   , only : r8 => shr_kind_r8
  use shr_infnan_mod                 , only : nan => shr_infnan_nan, assignment(=), shr_infnan_isnan
  use shr_log_mod                    , only : errMsg => shr_log_errMsg
  use abortutils                     , only : endrun
  !
  implicit none
  private

  ! Non-tunable constants
  real(r8) :: rgasm  ! J/mol.K; rgas / 1000; will be set below
  real(r8), parameter :: rgasLatm = 0.0821_r8 ! L.atm/mol.K

  ! !PUBLIC MEMBER FUNCTIONS:

  type, private :: params_type
     ! ch4 production constants
     real(r8) :: q10ch4               ! additional Q10 for methane production ABOVE the soil decomposition temperature relationship
     real(r8) :: q10ch4base           ! temperature at which the effective f_ch4 actually equals the constant f_ch4
     real(r8) :: f_ch4                ! ratio of CH4 production to total C mineralization
     real(r8) :: rootlitfrac          ! Fraction of soil organic matter associated with roots
     real(r8) :: cnscalefactor        ! scale factor on CN decomposition for assigning methane flux
     real(r8) :: redoxlag             ! Number of days to lag in the calculation of finundated_lag
     real(r8) :: lake_decomp_fact     ! Base decomposition rate (1/s) at 25C
     real(r8) :: redoxlag_vertical    ! time lag (days) to inhibit production for newly unsaturated layers
     real(r8) :: pHmax                ! maximum pH for methane production(= 9._r8)
     real(r8) :: pHmin                ! minimum pH for methane production(= 2.2_r8)
     real(r8) :: oxinhib              ! inhibition of methane production by oxygen (m^3/mol)

     ! ch4 oxidation constants
     real(r8) :: vmax_ch4_oxid        ! oxidation rate constant (= 45.e-6_r8 * 1000._r8 / 3600._r8) [mol/m3-w/s];
     real(r8) :: k_m                  ! Michaelis-Menten oxidation rate constant for CH4 concentration 
     real(r8) :: q10_ch4oxid          ! Q10 oxidation constant
     real(r8) :: smp_crit             ! Critical soil moisture potential
     real(r8) :: k_m_o2               ! Michaelis-Menten oxidation rate constant for O2 concentration
     real(r8) :: k_m_unsat            ! Michaelis-Menten oxidation rate constant for CH4 concentration
     real(r8) :: vmax_oxid_unsat      ! (= 45.e-6_r8 * 1000._r8 / 3600._r8 / 10._r8) [mol/m3-w/s]

     ! ch4 aerenchyma constants
     real(r8) :: aereoxid             ! fraction of methane flux entering aerenchyma rhizosphere that will be

     ! oxidized rather than emitted
     real(r8) :: scale_factor_aere    ! scale factor on the aerenchyma area for sensitivity tests
     real(r8) :: nongrassporosratio   ! Ratio of root porosity in non-grass to grass, used for aerenchyma transport
     real(r8) :: unsat_aere_ratio     ! Ratio to multiply upland vegetation aerenchyma porosity by compared to inundated systems (= 0.05_r8 / 0.3_r8)
     real(r8) :: porosmin             ! minimum aerenchyma porosity (unitless)(= 0.05_r8) 

     ! ch4 ebbulition constants
     real(r8) :: vgc_max              ! ratio of saturation pressure triggering ebullition

     ! ch4 transport constants
     real(r8) :: satpow               ! exponent on watsat for saturated soil solute diffusion
     real(r8) :: scale_factor_gasdiff ! For sensitivity tests; convection would allow this to be > 1
     real(r8) :: scale_factor_liqdiff ! For sensitivity tests; convection would allow this to be > 1
     real(r8) :: capthick             ! min thickness before assuming h2osfc is impermeable (mm) (= 100._r8)

     ! additional constants
     real(r8) :: f_sat                ! volumetric soil water defining top of water table or where production is allowed (=0.95)
     real(r8) :: qflxlagd             ! days to lag qflx_surf_lag in the tropics (days) ( = 30._r8)
     real(r8) :: highlatfact          ! multiple of qflxlagd for high latitudes	(= 2._r8)	
     real(r8) :: q10lakebase          ! (K) base temperature for lake CH4 production (= 298._r8)
     real(r8) :: atmch4               ! Atmospheric CH4 mixing ratio to prescribe if not provided by the atmospheric model (= 1.7e-6_r8) (mol/mol)
     real(r8) :: rob                  ! ratio of root length to vertical depth ("root obliquity") (= 3._r8)
  end type params_type
  type(params_type), private ::  params_inst

  type, public :: ch4_type
     real(r8), pointer, private :: ch4_prod_depth_sat_col     (:,:) ! col CH4 production rate from methanotrophs (mol/m3/s) (nlevsoi)
     real(r8), pointer, private :: ch4_prod_depth_unsat_col   (:,:) ! col CH4 production rate from methanotrophs (mol/m3/s) (nlevsoi)
     real(r8), pointer, private :: ch4_prod_depth_lake_col    (:,:) ! col CH4 production rate from methanotrophs (mol/m3/s) (nlevsoi)
     real(r8), pointer, private :: ch4_oxid_depth_sat_col     (:,:) ! col CH4 consumption rate via oxidation in each soil layer (mol/m3/s) (nlevsoi)
     real(r8), pointer, private :: ch4_oxid_depth_unsat_col   (:,:) ! col CH4 consumption rate via oxidation in each soil layer (mol/m3/s) (nlevsoi)
     real(r8), pointer, private :: ch4_oxid_depth_lake_col    (:,:) ! col CH4 consumption rate via oxidation in each soil layer (mol/m3/s) (nlevsoi)
     real(r8), pointer, private :: ch4_aere_depth_sat_col     (:,:) ! col CH4 loss rate via aerenchyma in each soil layer (mol/m3/s) (nlevsoi)
     real(r8), pointer, private :: ch4_aere_depth_unsat_col   (:,:) ! col CH4 loss rate via aerenchyma in each soil layer (mol/m3/s) (nlevsoi)
     real(r8), pointer, private :: ch4_tran_depth_sat_col     (:,:) ! col CH4 loss rate via transpiration in each soil layer (mol/m3/s) (nlevsoi)
     real(r8), pointer, private :: ch4_tran_depth_unsat_col   (:,:) ! col CH4 loss rate via transpiration in each soil layer (mol/m3/s) (nlevsoi)
     real(r8), pointer, private :: ch4_ebul_depth_sat_col     (:,:) ! col CH4 loss rate via ebullition in each soil layer (mol/m3/s) (nlevsoi)
     real(r8), pointer, private :: ch4_ebul_depth_unsat_col   (:,:) ! col CH4 loss rate via ebullition in each soil layer (mol/m3/s) (nlevsoi)
     real(r8), pointer, private :: ch4_ebul_total_sat_col     (:)   ! col Total col CH4 ebullition (mol/m2/s)
     real(r8), pointer, private :: ch4_ebul_total_unsat_col   (:)   ! col Total col CH4 ebullition (mol/m2/s)
     real(r8), pointer, private :: ch4_surf_aere_sat_col      (:)   ! col CH4 aerenchyma flux to atmosphere (after oxidation) (mol/m2/s)
     real(r8), pointer, private :: ch4_surf_aere_unsat_col    (:)   ! col CH4 aerenchyma flux to atmosphere (after oxidation) (mol/m2/s)
     real(r8), pointer, private :: ch4_surf_ebul_sat_col      (:)   ! col CH4 ebullition flux to atmosphere (after oxidation) (mol/m2/s)
     real(r8), pointer, private :: ch4_surf_ebul_unsat_col    (:)   ! col CH4 ebullition flux to atmosphere (after oxidation) (mol/m2/s)
     real(r8), pointer, private :: ch4_surf_ebul_lake_col     (:)   ! col CH4 ebullition flux to atmosphere (after oxidation) (mol/m2/s)
     real(r8), pointer, private :: co2_aere_depth_sat_col     (:,:) ! col CO2 loss rate via aerenchyma in each soil layer (mol/m3/s) (nlevsoi)
     real(r8), pointer, private :: co2_aere_depth_unsat_col   (:,:) ! col CO2 loss rate via aerenchyma in each soil layer (mol/m3/s) (nlevsoi)
     real(r8), pointer, private :: o2_oxid_depth_sat_col      (:,:) ! col O2 consumption rate via oxidation in each soil layer (mol/m3/s) (nlevsoi)
     real(r8), pointer, private :: o2_oxid_depth_unsat_col    (:,:) ! col O2 consumption rate via oxidation in each soil layer (mol/m3/s) (nlevsoi)
     real(r8), pointer, private :: o2_aere_depth_sat_col      (:,:) ! col O2 gain rate via aerenchyma in each soil layer (mol/m3/s) (nlevsoi)
     real(r8), pointer, private :: o2_aere_depth_unsat_col    (:,:) ! col O2 gain rate via aerenchyma in each soil layer (mol/m3/s) (nlevsoi)
     real(r8), pointer, private :: co2_decomp_depth_sat_col   (:,:) ! col CO2 production during decomposition in each soil layer (nlevsoi) (mol/m3/s)
     real(r8), pointer, private :: co2_decomp_depth_unsat_col (:,:) ! col CO2 production during decomposition in each soil layer (nlevsoi) (mol/m3/s)
     real(r8), pointer, private :: co2_oxid_depth_sat_col     (:,:) ! col CO2 production rate via oxidation in each soil layer (mol/m3/s) (nlevsoi)
     real(r8), pointer, private :: co2_oxid_depth_unsat_col   (:,:) ! col CO2 production rate via oxidation in each soil layer (mol/m3/s) (nlevsoi)
     real(r8), pointer, private :: conc_o2_lake_col           (:,:) ! col O2 conc in each soil layer (mol/m3) (nlevsoi)
     real(r8), pointer, private :: conc_ch4_sat_col           (:,:) ! col CH4 conc in each soil layer (mol/m3) (nlevsoi)
     real(r8), pointer, private :: conc_ch4_unsat_col         (:,:) ! col CH4 conc in each soil layer (mol/m3) (nlevsoi)
     real(r8), pointer, private :: conc_ch4_lake_col          (:,:) ! col CH4 conc in each soil layer (mol/m3) (nlevsoi)
     real(r8), pointer, private :: ch4_surf_diff_sat_col      (:)   ! col CH4 surface flux (mol/m2/s)
     real(r8), pointer, private :: ch4_surf_diff_unsat_col    (:)   ! col CH4 surface flux (mol/m2/s)
     real(r8), pointer, private :: ch4_surf_diff_lake_col     (:)   ! col CH4 surface flux (mol/m2/s)
     real(r8), pointer, private :: ch4_dfsat_flux_col         (:)   ! col CH4 flux to atm due to decreasing fsat (kg C/m^2/s) [+]

     real(r8), pointer, private :: zwt_ch4_unsat_col          (:)   ! col depth of water table for unsaturated fraction (m)
     real(r8), pointer, private :: lake_soilc_col             (:,:) ! col total soil organic matter found in level (g C / m^3) (nlevsoi)
     real(r8), pointer, private :: totcolch4_col              (:)   ! col total methane found in soil col (g C / m^2)
     real(r8), pointer, private :: totcolch4_bef_col          (:)   ! col total methane found in soil col, start of timestep (g C / m^2)
     real(r8), pointer, private :: annsum_counter_col         (:)   ! col seconds since last annual accumulator turnover
     real(r8), pointer, private :: tempavg_somhr_col          (:)   ! col temporary average SOM heterotrophic resp. (gC/m2/s)
     real(r8), pointer, private :: annavg_somhr_col           (:)   ! col annual average SOM heterotrophic resp. (gC/m2/s)
     real(r8), pointer, private :: tempavg_finrw_col          (:)   ! col respiration-weighted annual average of finundated
     real(r8), pointer, private :: annavg_finrw_col           (:)   ! col respiration-weighted annual average of finundated
     real(r8), pointer, private :: sif_col                    (:)   ! col (unitless) ratio applied to sat. prod. to account for seasonal inundation
     real(r8), pointer, private :: ch4stress_unsat_col        (:,:) ! col Ratio of methane available to the total per-timestep methane sinks (nlevsoi)
     real(r8), pointer, private :: ch4stress_sat_col          (:,:) ! col Ratio of methane available to the total per-timestep methane sinks (nlevsoi)
     real(r8), pointer, private :: qflx_surf_lag_col          (:)   ! col time-lagged surface runoff (mm H2O /s)
     real(r8), pointer, private :: finundated_lag_col         (:)   ! col time-lagged fractional inundated area
     real(r8), pointer, private :: layer_sat_lag_col          (:,:) ! col Lagged saturation status of soil layer in the unsaturated zone (1 = sat)
     real(r8), pointer, private :: zwt0_col                   (:)   ! col coefficient for determining finundated (m)
     real(r8), pointer, private :: f0_col                     (:)   ! col maximum inundated fraction for a gridcell (for methane code)
     real(r8), pointer, private :: p3_col                     (:)   ! col coefficient for determining finundated (m)
     real(r8), pointer, private :: pH_col                     (:)   ! col pH values for methane production
     !
     real(r8), pointer, private :: dyn_ch4bal_adjustments_col (:)   ! adjustments to each column made in this timestep via dynamic column area adjustments (only makes sense at the column-level: meaningless if averaged to the gridcell-level) (g C / m^2)
     !
     real(r8), pointer, private :: c_atm_grc                  (:,:) ! grc atmospheric conc of CH4, O2, CO2 (mol/m3)
     real(r8), pointer, private :: ch4co2f_grc                (:)   ! grc CO2 production from CH4 oxidation (g C/m**2/s)
     real(r8), pointer, private :: ch4prodg_grc               (:)   ! grc average CH4 production (g C/m^2/s)
     !
     ! for aerenchyma calculations 
     real(r8), pointer, private :: annavg_agnpp_patch         (:)   ! patch (gC/m2/s) annual average aboveground NPP
     real(r8), pointer, private :: annavg_bgnpp_patch         (:)   ! patch (gC/m2/s) annual average belowground NPP
     real(r8), pointer, private :: tempavg_agnpp_patch        (:)   ! patch (gC/m2/s) temp. average aboveground NPP
     real(r8), pointer, private :: tempavg_bgnpp_patch        (:)   ! patch (gC/m2/s) temp. average belowground NPP
     !
     ! The following variable reports whether this is the first timestep that includes
     ! ch4. It is true in the first timestep of the run, and remains true until the
     ! methane code is first run - at which point it becomes false, and remains
     ! false. This could be a scalar, but scalars cause problems with threading, so we use
     ! a column-level array (column-level for convenience, because it is referenced in
     ! column-level loops).
     logical , pointer, private :: ch4_first_time_col         (:)   ! col whether this is the first time step that includes ch4
     !
     real(r8), pointer, public :: finundated_col             (:)   ! col fractional inundated area (excluding dedicated wetland cols)
     real(r8), pointer, public :: finundated_pre_snow_col    (:)   ! col fractional inundated area (excluding dedicated wetland cols) before snow
     real(r8), pointer, public :: o2stress_unsat_col         (:,:) ! col Ratio of oxygen available to that demanded by roots, aerobes, & methanotrophs (nlevsoi)
     real(r8), pointer, public :: o2stress_sat_col           (:,:) ! col Ratio of oxygen available to that demanded by roots, aerobes, & methanotrophs (nlevsoi)
     real(r8), pointer, public :: conc_o2_sat_col            (:,:) ! col O2 conc in each soil layer (mol/m3) (nlevsoi)
     real(r8), pointer, public :: conc_o2_unsat_col          (:,:) ! col O2 conc in each soil layer (mol/m3) (nlevsoi)
     real(r8), pointer, public :: o2_decomp_depth_sat_col    (:,:) ! col O2 consumption during decomposition in each soil layer (nlevsoi) (mol/m3/s)
     real(r8), pointer, public :: o2_decomp_depth_unsat_col  (:,:) ! col O2 consumption during decomposition in each soil layer (nlevsoi) (mol/m3/s)
     real(r8), pointer, public :: ch4_surf_flux_tot_col      (:)   ! col CH4 surface flux (to atm) (kg C/m**2/s)

     real(r8), pointer, public :: grnd_ch4_cond_patch        (:)   ! patch tracer conductance for boundary layer [m/s]
     real(r8), pointer, public :: grnd_ch4_cond_col          (:)   ! col tracer conductance for boundary layer [m/s]

  end type ch4_type

  character(len=*), parameter, private :: sourcefile = &
       __FILE__
  !------------------------------------------------------------------------

end module ch4Mod

