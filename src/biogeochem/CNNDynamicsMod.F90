module CNNDynamicsMod

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Module for mineral nitrogen dynamics (deposition, fixation, leaching)
  ! for coupled carbon-nitrogen code.
  !
  ! !USES:
  use shr_kind_mod                    , only : r8 => shr_kind_r8
  use decompMod                       , only : bounds_type
  use clm_varcon                      , only : dzsoi_decomp, zisoi
  use clm_varctl                      , only : nfix_timeconst
  use subgridAveMod                   , only : p2c
  use atm2lndType                     , only : atm2lnd_type
  use CNVegStateType                  , only : cnveg_state_type
  use CNVegCarbonFluxType             , only : cnveg_carbonflux_type
  use CNVegNitrogenStateType	      , only : cnveg_nitrogenstate_type
  use CNVegNitrogenFluxType	      , only : cnveg_nitrogenflux_type
  use SoilBiogeochemStateType         , only : soilbiogeochem_state_type
  use SoilBiogeochemNitrogenStateType , only : soilbiogeochem_nitrogenstate_type
  use SoilBiogeochemNitrogenFluxType  , only : soilbiogeochem_nitrogenflux_type
  use WaterStateType                  , only : waterstate_type
  use WaterFluxType                   , only : waterflux_type
  use CropType                        , only : crop_type
  use ColumnType                      , only : col                
  use PatchType                       , only : patch                
  use perf_mod                        , only : t_startf, t_stopf
  !
  implicit none
  private
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: CNNDeposition               ! Update N deposition rate from atm forcing
  public :: CNNFixation                 ! Update N Fixation rate
  public :: CNNFert                     ! Update N fertilizer for crops
  public :: CNSoyfix                    ! N Fixation for soybeans
  public :: CNFreeLivingFixation        ! N free living fixation

  !
  ! !PRIVATE DATA MEMBERS:
  type, private :: params_type
     real(r8) :: freelivfix_intercept   ! intercept of line of free living fixation with annual ET
     real(r8) :: freelivfix_slope_wET   ! slope of line of free living fixation with annual ET
  end type params_type
  type(params_type) :: params_inst
  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine CNNDeposition( bounds, &
       atm2lnd_inst, soilbiogeochem_nitrogenflux_inst )
    !
    ! !DESCRIPTION:
    ! On the radiation time step, update the nitrogen deposition rate
    ! from atmospheric forcing. For now it is assumed that all the atmospheric
    ! N deposition goes to the soil mineral N pool.
    ! This could be updated later to divide the inputs between mineral N absorbed
    ! directly into the canopy and mineral N entering the soil pool.
    !
    ! !USES:
    ! !ARGUMENTS:
    type(bounds_type)        , intent(in)    :: bounds  
    type(atm2lnd_type)       , intent(in)    :: atm2lnd_inst
    type(soilbiogeochem_nitrogenflux_type) , intent(inout) :: soilbiogeochem_nitrogenflux_inst
    !
    ! !LOCAL VARIABLES:
    integer :: g,c                    ! indices
    !-----------------------------------------------------------------------
    
    associate(                                                                & 
         forc_ndep     =>  atm2lnd_inst%forc_ndep_grc ,                       & ! Input:  [real(r8) (:)]  nitrogen deposition rate (gN/m2/s)
         ndep_to_sminn =>  soilbiogeochem_nitrogenflux_inst%ndep_to_sminn_col & ! Output: [real(r8) (:)]  atmospheric N deposition to soil mineral N (gN/m2/s)
         )
      
      ! Loop through columns
      do c = bounds%begc, bounds%endc
         g = col%gridcell(c)
         ndep_to_sminn(c) = forc_ndep(g)

      end do

    end associate

  end subroutine CNNDeposition

  !-----------------------------------------------------------------------
  subroutine CNFreeLivingFixation(num_soilc, filter_soilc, &
       waterflux_inst, soilbiogeochem_nitrogenflux_inst)


    use clm_time_manager , only : get_days_per_year, get_step_size
    use shr_sys_mod      , only : shr_sys_flush
    use clm_varcon       , only : secspday, spval
 
    integer                                , intent(in)    :: num_soilc       ! number of soil columns in filter                                                                                                                     
    integer                                , intent(in)    :: filter_soilc(:) ! filter for soil columns                                                                                                                                  
   
    type(soilbiogeochem_nitrogenflux_type) , intent(inout) :: soilbiogeochem_nitrogenflux_inst 
    type(waterflux_type)                   , intent(inout) :: waterflux_inst 
    !
    ! !LOCAL VARIABLES:                                                                                                                                                                                                           
    integer  :: c,fc            !indices     
    real(r8) :: dayspyr         !days per year 
    real(r8) :: secs_per_year   !seconds per year   

       associate(                                                                        &
                  AnnET            => waterflux_inst%AnnET,                              & ! Input:  [real(:)  ] : Annual average ET flux mmH20/s
                  freelivfix_slope => params_inst%freelivfix_slope_wET,                  & ! Input:  [real     ] : slope of fixation with ET
                  freelivfix_inter => params_inst%freelivfix_intercept,                  & ! Input:  [real     ] : intercept of fixation with ET
                  ffix_to_sminn    => soilbiogeochem_nitrogenflux_inst%ffix_to_sminn_col & ! Output: [real(:)  ] : free living N fixation to soil mineral N (gN/m2/s)
                ) 
       
       dayspyr = get_days_per_year()
       secs_per_year = dayspyr*24_r8*3600_r8

       do fc = 1,num_soilc
           c = filter_soilc(fc)
          ffix_to_sminn(c) = (freelivfix_slope*(max(0._r8,AnnET(c))*secs_per_year) + freelivfix_inter )/secs_per_year !(units g N m-2 s-1)  

       end do

  end associate
  end subroutine CNFreeLivingFixation

  !-----------------------------------------------------------------------
  subroutine CNNFixation(num_soilc, filter_soilc, &
       cnveg_carbonflux_inst, soilbiogeochem_nitrogenflux_inst)
    !
    ! !DESCRIPTION:
    ! On the radiation time step, update the nitrogen fixation rate
    ! as a function of annual total NPP. This rate gets updated once per year.
    ! All N fixation goes to the soil mineral N pool.
    !
    ! !USES:
    use clm_time_manager , only : get_days_per_year, get_step_size
    use shr_sys_mod      , only : shr_sys_flush
    use clm_varcon       , only : secspday, spval
    !
    ! !ARGUMENTS:
    integer                                , intent(in)    :: num_soilc       ! number of soil columns in filter
    integer                                , intent(in)    :: filter_soilc(:) ! filter for soil columns
    type(cnveg_carbonflux_type)            , intent(inout) :: cnveg_carbonflux_inst
    type(soilbiogeochem_nitrogenflux_type) , intent(inout) :: soilbiogeochem_nitrogenflux_inst 
    !
    ! !LOCAL VARIABLES:
    integer  :: c,fc                  ! indices
    real(r8) :: t                     ! temporary
    real(r8) :: dayspyr               ! days per year
    !-----------------------------------------------------------------------

    associate(                                                                & 
         cannsum_npp    => cnveg_carbonflux_inst%annsum_npp_col ,             & ! Input:  [real(r8) (:)]  nitrogen deposition rate (gN/m2/s)                
         col_lag_npp    => cnveg_carbonflux_inst%lag_npp_col    ,             & ! Input: [real(r8) (:)]  (gC/m2/s) lagged net primary production           

         nfix_to_sminn  => soilbiogeochem_nitrogenflux_inst%nfix_to_sminn_col & ! Output: [real(r8) (:)]  symbiotic/asymbiotic N fixation to soil mineral N (gN/m2/s)
         )

      dayspyr = get_days_per_year()

      if ( nfix_timeconst > 0._r8 .and. nfix_timeconst < 500._r8 ) then
         ! use exponential relaxation with time constant nfix_timeconst for NPP - NFIX relation
         ! Loop through columns
         do fc = 1,num_soilc
            c = filter_soilc(fc)         

            if (col_lag_npp(c) /= spval) then
               ! need to put npp in units of gC/m^2/year here first
               t = (1.8_r8 * (1._r8 - exp(-0.003_r8 * col_lag_npp(c)*(secspday * dayspyr))))/(secspday * dayspyr)  
               nfix_to_sminn(c) = max(0._r8,t)
            else
               nfix_to_sminn(c) = 0._r8
            endif
         end do
      else
         ! use annual-mean values for NPP-NFIX relation
         do fc = 1,num_soilc
            c = filter_soilc(fc)

            t = (1.8_r8 * (1._r8 - exp(-0.003_r8 * cannsum_npp(c))))/(secspday * dayspyr)
            nfix_to_sminn(c) = max(0._r8,t)
         end do
      endif

    end associate

  end subroutine CNNFixation
 
  !-----------------------------------------------------------------------
  subroutine CNNFert(bounds, num_soilc, filter_soilc, &
       cnveg_nitrogenflux_inst, soilbiogeochem_nitrogenflux_inst)
    !
    ! !DESCRIPTION:
    ! On the radiation time step, update the nitrogen fertilizer for crops
    ! All fertilizer goes into the soil mineral N pool.
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    type(bounds_type)                      , intent(in)    :: bounds  
    integer                                , intent(in)    :: num_soilc       ! number of soil columns in filter
    integer                                , intent(in)    :: filter_soilc(:) ! filter for soil columns
    type(cnveg_nitrogenflux_type)          , intent(in)    :: cnveg_nitrogenflux_inst
    type(soilbiogeochem_nitrogenflux_type) , intent(inout) :: soilbiogeochem_nitrogenflux_inst 
    !
    ! !LOCAL VARIABLES:
    integer :: c,fc                 ! indices
    !-----------------------------------------------------------------------

    associate(                                                                  &   
         fert          =>    cnveg_nitrogenflux_inst%fert_patch ,               & ! Input:  [real(r8) (:)]  nitrogen fertilizer rate (gN/m2/s)                
         fert_to_sminn =>    soilbiogeochem_nitrogenflux_inst%fert_to_sminn_col & ! Output: [real(r8) (:)]                                                    
         )
      
      call p2c(bounds, num_soilc, filter_soilc, &
           fert(bounds%begp:bounds%endp), &
           fert_to_sminn(bounds%begc:bounds%endc))

    end associate

  end subroutine CNNFert

  !-----------------------------------------------------------------------
  subroutine CNSoyfix (bounds, num_soilc, filter_soilc, num_soilp, filter_soilp, &
       waterstate_inst, crop_inst, cnveg_state_inst, cnveg_nitrogenflux_inst , &
       soilbiogeochem_state_inst, soilbiogeochem_nitrogenstate_inst, soilbiogeochem_nitrogenflux_inst)
    !
    ! !DESCRIPTION:
    ! This routine handles the fixation of nitrogen for soybeans based on
    ! the EPICPHASE model M. Cabelguenne et al., Agricultural systems 60: 175-196, 1999
    ! N-fixation is based on soil moisture, plant growth phase, and availibility of
    ! nitrogen in the soil root zone.
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    type(bounds_type)                       , intent(in)    :: bounds  
    integer                                 , intent(in)    :: num_soilc       ! number of soil columns in filter
    integer                                 , intent(in)    :: filter_soilc(:) ! filter for soil columns
    integer                                 , intent(in)    :: num_soilp       ! number of soil patches in filter
    integer                                 , intent(in)    :: filter_soilp(:) ! filter for soil patches
    type(waterstate_type)                   , intent(in)    :: waterstate_inst
    type(crop_type)                         , intent(in)    :: crop_inst
    type(cnveg_state_type)                  , intent(in)    :: cnveg_state_inst
    type(cnveg_nitrogenflux_type)           , intent(inout) :: cnveg_nitrogenflux_inst
    type(soilbiogeochem_state_type)         , intent(in)    :: soilbiogeochem_state_inst
    type(soilbiogeochem_nitrogenstate_type) , intent(in)    :: soilbiogeochem_nitrogenstate_inst
    type(soilbiogeochem_nitrogenflux_type)  , intent(inout) :: soilbiogeochem_nitrogenflux_inst 
    !
    ! !LOCAL VARIABLES:
    integer :: fp,p,c
    real(r8):: fxw,fxn,fxg,fxr             ! soil water factor, nitrogen factor, growth stage factor
    real(r8):: soy_ndemand                 ! difference between nitrogen supply and demand
    real(r8):: GDDfrac
    real(r8):: sminnthreshold1, sminnthreshold2
    real(r8):: GDDfracthreshold1, GDDfracthreshold2
    real(r8):: GDDfracthreshold3, GDDfracthreshold4
    !-----------------------------------------------------------------------

    associate(                                                                      & 
         wf               =>  waterstate_inst%wf_col                      ,         & ! Input:  [real(r8) (:) ]  soil water as frac. of whc for top 0.5 m          

         hui              =>  crop_inst%gddplant_patch                    ,         & ! Input:  [real(r8) (:) ]  gdd since planting (gddplant)                    
         croplive         =>  crop_inst%croplive_patch                    ,         & ! Input:  [logical  (:) ]  true if planted and not harvested                  

         gddmaturity      =>  cnveg_state_inst%gddmaturity_patch          ,         & ! Input:  [real(r8) (:) ]  gdd needed to harvest                             

         plant_ndemand    =>  cnveg_nitrogenflux_inst%plant_ndemand_patch ,         & ! Input:  [real(r8) (:) ]  N flux required to support initial GPP (gN/m2/s)  
         soyfixn          =>  cnveg_nitrogenflux_inst%soyfixn_patch       ,         & ! Output: [real(r8) (:) ]  nitrogen fixed to each soybean crop               

         fpg              =>  soilbiogeochem_state_inst%fpg_col           ,         & ! Input:  [real(r8) (:) ]  fraction of potential gpp (no units)              

         sminn            =>  soilbiogeochem_nitrogenstate_inst%sminn_col ,         & ! Input:  [real(r8) (:) ]  (kgN/m2) soil mineral N                           
         soyfixn_to_sminn =>  soilbiogeochem_nitrogenflux_inst%soyfixn_to_sminn_col & ! Output: [real(r8) (:) ]                                                    
         )

      sminnthreshold1 = 30._r8
      sminnthreshold2 = 10._r8
      GDDfracthreshold1 = 0.15_r8
      GDDfracthreshold2 = 0.30_r8
      GDDfracthreshold3 = 0.55_r8
      GDDfracthreshold4 = 0.75_r8

      do fp = 1,num_soilp
         p = filter_soilp(fp)
         soyfixn(p) = 0._r8
      end do

      call p2c(bounds, num_soilc, filter_soilc, &
           soyfixn(bounds%begp:bounds%endp), &
           soyfixn_to_sminn(bounds%begc:bounds%endc))

    end associate

  end subroutine CNSoyfix

end module CNNDynamicsMod
