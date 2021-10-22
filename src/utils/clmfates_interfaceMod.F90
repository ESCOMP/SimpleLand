module CLMFatesInterfaceMod
   
   ! -------------------------------------------------------------------------------------
   ! This module contains various functions and definitions to aid in the
   ! coupling of the FATES library/API with the CLM/ALM/ATS/etc model driver.  
   ! All connections between the two models should occur in this file alone.  
   ! 
   ! This is also the only location where CLM code is allowed to see FATES memory 
   ! structures.
   ! The routines here, that call FATES library routines, will not pass any types defined
   ! by the driving land model (HLM).
   ! 
   ! either native type arrays (int,real,log, etc) or packed into fates boundary condition
   ! structures.
   !
   ! Note that CLM/ALM does use Shared Memory Parallelism (SMP), where processes such as 
   ! the update of state variables are forked.  However, IO is not assumed to be 
   ! threadsafe and therefore memory spaces reserved for IO must be continuous vectors,
   ! and moreover they must be pushed/pulled from history IO for each individual 
   ! bounds_proc memory space as a unit.
   !
   ! Therefore, the state variables in the clm_fates communicator is vectorized by
   ! threadcount, and the IO communication arrays are not.
   !
   !
   ! Conventions:
   ! keep line widths within 90 spaces
   ! HLM acronym = Host Land Model
   !
   ! -------------------------------------------------------------------------------------

   !  use ed_driver_interface, only: 
   
   ! Used CLM Modules
   use PatchType         , only : patch
   use shr_kind_mod      , only : r8 => shr_kind_r8
   use decompMod         , only : bounds_type
   use WaterStateType    , only : waterstate_type
   use WaterFluxType     , only : waterflux_type
   use CanopyStateType   , only : canopystate_type
   use TemperatureType   , only : temperature_type
   use EnergyFluxType    , only : energyflux_type

   use SoilStateType     , only : soilstate_type 
   use clm_varctl        , only : iulog
   use clm_varctl        , only : use_vertsoilc 
   use clm_varctl        , only : use_fates_spitfire
   use clm_varctl        , only : use_fates_planthydro
   use clm_varctl        , only : use_fates_ed_st3
   use clm_varctl        , only : use_fates_ed_prescribed_phys
   use clm_varctl        , only : use_fates_logging
   use clm_varctl        , only : use_fates_inventory_init
   use clm_varctl        , only : fates_inventory_ctrl_filename
 
   use clm_varcon        , only : tfrz
   use clm_varcon        , only : spval 
   use clm_varcon        , only : denice
   use clm_varcon        , only : ispval

   use clm_varpar        , only : natpft_size
   use clm_varpar        , only : numrad
   use clm_varpar        , only : ivis
   use clm_varpar        , only : inir
   use clm_varpar        , only : nlevgrnd
   use clm_varpar        , only : nlevsoi
   use clm_varpar        , only : nlevdecomp
   use clm_varpar        , only : nlevdecomp_full
   use PhotosynthesisMod , only : photosyns_type
   use atm2lndType       , only : atm2lnd_type
   use SurfaceAlbedoType , only : surfalb_type
   use SolarAbsorbedType , only : solarabs_type
   use SoilBiogeochemCarbonFluxType, only :  soilbiogeochem_carbonflux_type
   use SoilBiogeochemCarbonStateType, only : soilbiogeochem_carbonstate_type
   use FrictionVelocityMod  , only : frictionvel_type
   use clm_time_manager  , only : is_restart
   use ncdio_pio         , only : file_desc_t, ncd_int, ncd_double
   use restUtilMod,        only : restartvar
   use clm_time_manager  , only : get_days_per_year, &
                                  get_curr_date,     &
                                  get_ref_date,      &
                                  timemgr_datediff,  &
                                  is_beg_curr_day,   &
                                  get_step_size,     &
                                  get_nstep
   use spmdMod           , only : masterproc
   use decompMod         , only : get_proc_bounds,   &
                                  get_proc_clumps,   &
                                  get_clump_bounds
   use GridCellType      , only : grc
   use ColumnType        , only : col
   use LandunitType      , only : lun
   use landunit_varcon   , only : istsoil
   use abortutils        , only : endrun
   use shr_log_mod       , only : errMsg => shr_log_errMsg    
   use clm_varcon        , only : dzsoi_decomp
!   use SoilWaterPlantSinkMod, only : Compute_EffecRootFrac_And_VertTranSink_Default

   ! Used FATES Modules
   use FatesInterfaceMod     , only : fates_interface_type
   use FatesHistoryInterfaceMod, only : fates_history_interface_type
   use FatesRestartInterfaceMod, only : fates_restart_interface_type

   implicit none
   
   type, public :: f2hmap_type

      ! This is the associated column index of each FATES site
      integer, allocatable :: fcolumn (:) 

      ! This is the associated site index of any HLM columns
      ! This vector may be sparse, and non-sites have index 0
      integer, allocatable :: hsites  (:)

   end type f2hmap_type
   

   type, public :: hlm_fates_interface_type
      
      !      private
      

      ! See above for descriptions of the sub-types populated
      ! by thread.  This type is somewhat self-explanatory, in that it simply
      ! breaks up memory and process by thread.  Each thread will have its
      ! own list of sites, and boundary conditions for those sites

      type(fates_interface_type), allocatable :: fates (:)
      

      ! This memory structure is used to map fates sites
      ! into the host model.  Currently, the FATES site
      ! and its column number matching are its only members

      type(f2hmap_type), allocatable  :: f2hmap(:)

      ! fates_hist is the interface class for the history output
      type(fates_history_interface_type) :: fates_hist

      ! fates_restart is the inteface calss for restarting the model
      type(fates_restart_interface_type) :: fates_restart

   contains
      
      !procedure, public :: init
      !procedure, public :: check_hlm_active
      !procedure, public :: restart
      !procedure, public :: init_coldstart
      !procedure, public :: dynamics_driv
      !procedure, public :: wrap_sunfrac
      !procedure, public :: wrap_btran
      !procedure, public :: wrap_photosynthesis
      !procedure, public :: wrap_accumulatefluxes
      !procedure, public :: prep_canopyfluxes
      !procedure, public :: wrap_canopy_radiation
      !procedure, public :: wrap_bgc_summary
      !procedure, public :: TransferZ0mDisp
      !procedure, private :: init_history_io
      !procedure, private :: wrap_update_hlmfates_dyn
      !procedure, private :: init_soil_depths
      !procedure, public  :: ComputeRootSoilFlux
      !procedure, public  :: wrap_hydraulics_drive

   end type hlm_fates_interface_type

   logical :: DEBUG  = .false.

   character(len=*), parameter, private :: sourcefile = &
        __FILE__
   
end module CLMFatesInterfaceMod
