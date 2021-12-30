module SoilBiogeochemDecompCascadeBGCMod

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Sets the coeffiecients used in the decomposition cascade submodel.  
  ! This uses the CENTURY/BGC parameters
  !
  ! !USES:
  use shr_kind_mod                       , only : r8 => shr_kind_r8
  use shr_const_mod                      , only : SHR_CONST_TKFRZ
  use shr_log_mod                        , only : errMsg => shr_log_errMsg
  use clm_varpar                         , only : nlevsoi, nlevgrnd, nlevdecomp, ndecomp_cascade_transitions, ndecomp_pools
  use clm_varpar                         , only : i_met_lit, i_cel_lit, i_lig_lit, i_cwd
  use clm_varctl                         , only : iulog, spinup_state, anoxia, use_vertsoilc, use_fates
  use clm_varcon                         , only : zsoi
  use decompMod                          , only : bounds_type
  use spmdMod                            , only : masterproc
  use abortutils                         , only : endrun
  use CNSharedParamsMod                  , only : CNParamsShareInst, anoxia_wtsat, nlev_soildecomp_standard 
  use SoilBiogeochemDecompCascadeConType , only : decomp_cascade_con
  use SoilBiogeochemStateType            , only : soilbiogeochem_state_type
  use SoilBiogeochemCarbonFluxType       , only : soilbiogeochem_carbonflux_type
  use SoilStateType                      , only : soilstate_type
  use CanopyStateType                    , only : canopystate_type
  use TemperatureType                    , only : temperature_type 
  use ch4Mod                             , only : ch4_type
  use ColumnType                         , only : col                
  use GridcellType                       , only : grc
  use SoilBiogeochemStateType            , only : get_spinup_latitude_term

  !
  implicit none
  private
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: readParams                      ! Read in parameters from params file
  public :: init_decompcascade_bgc          ! Initialization
  !
  ! !PUBLIC DATA MEMBERS 
  logical , public :: normalize_q10_to_century_tfunc = .true.! do we normalize the century decomp. rates so that they match the CLM Q10 at a given tep?
  logical , public :: use_century_tfunc = .false.
  real(r8), public :: normalization_tref = 15._r8            ! reference temperature for normalizaion (degrees C)
  !
  ! !PRIVATE DATA MEMBERS 

  integer, private            :: i_soil1   = -9         ! Soil Organic Matter (SOM) first pool
  integer, private            :: i_soil2   = -9         ! SOM second pool
  integer, private            :: i_soil3   = -9         ! SOM third pool
  integer, private, parameter :: nsompools = 3          ! Number of SOM pools
  integer, private, parameter :: i_litr1   = i_met_lit  ! First litter pool, metobolic
  integer, private, parameter :: i_litr2   = i_cel_lit  ! Second litter pool, cellulose
  integer, private, parameter :: i_litr3   = i_lig_lit  ! Third litter pool, lignin

  type, private :: params_type
     real(r8):: cn_s1_bgc     !C:N for SOM 1
     real(r8):: cn_s2_bgc     !C:N for SOM 2
     real(r8):: cn_s3_bgc     !C:N for SOM 3

     real(r8):: rf_l1s1_bgc   !respiration fraction litter 1 -> SOM 1
     real(r8):: rf_l2s1_bgc
     real(r8):: rf_l3s2_bgc

     real(r8):: rf_s2s1_bgc    
     real(r8):: rf_s2s3_bgc    
     real(r8):: rf_s3s1_bgc    

     real(r8):: rf_cwdl2_bgc 
     real(r8):: rf_cwdl3_bgc

     real(r8):: tau_l1_bgc    ! turnover time of  litter 1 (yr)
     real(r8):: tau_l2_l3_bgc ! turnover time of  litter 2 and litter 3 (yr)
     real(r8):: tau_s1_bgc    ! turnover time of  SOM 1 (yr)
     real(r8):: tau_s2_bgc    ! turnover time of  SOM 2 (yr)
     real(r8):: tau_s3_bgc    ! turnover time of  SOM 3 (yr)
     real(r8):: tau_cwd_bgc   ! corrected fragmentation rate constant CWD

     real(r8) :: cwd_fcel_bgc !cellulose fraction for CWD
     real(r8) :: cwd_flig_bgc !

     real(r8) :: k_frag_bgc   !fragmentation rate for CWD
     real(r8) :: minpsi_bgc   !minimum soil water potential for heterotrophic resp
     real(r8) :: maxpsi_bgc   !maximum soil water potential for heterotrophic resp

     real(r8) :: initial_Cstocks(nsompools) ! Initial Carbon stocks for a cold-start
     real(r8) :: initial_Cstocks_depth      ! Soil depth for initial Carbon stocks for a cold-start
     
  end type params_type
  !
  type(params_type), private :: params_inst

  character(len=*), parameter, private :: sourcefile = &
       __FILE__

  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine readParams ( ncid )
    !
    ! !DESCRIPTION:
    !
    ! !USES:
    use ncdio_pio    , only: file_desc_t,ncd_io
    !
    ! !ARGUMENTS:
    type(file_desc_t),intent(inout) :: ncid   ! pio netCDF file id
    !
    ! !LOCAL VARIABLES:
    character(len=32)  :: subname = 'CNDecompBgcParamsType'
    character(len=100) :: errCode = 'Error reading in CN const file '
    logical            :: readv   ! has variable been read in or not
    real(r8)           :: tempr   ! temporary to read in constant
    character(len=100) :: tString ! temp. var for reading
    !-----------------------------------------------------------------------

    ! Read off of netcdf file
    tString='tau_l1'
    call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))
    params_inst%tau_l1_bgc=tempr

    tString='tau_l2_l3'
    call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))
    params_inst%tau_l2_l3_bgc=tempr

    tString='tau_s1'
    call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))
    params_inst%tau_s1_bgc=tempr

    tString='tau_s2'
    call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))
    params_inst%tau_s2_bgc=tempr

    tString='tau_s3'
    call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))
    params_inst%tau_s3_bgc=tempr

    tString='tau_cwd'
    call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))
    params_inst%tau_cwd_bgc=tempr

    tString='cn_s1_bgc'
    call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))
    params_inst%cn_s1_bgc=tempr

    tString='cn_s2_bgc'
    call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))
    params_inst%cn_s2_bgc=tempr

    tString='cn_s3_bgc'
    call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))
    params_inst%cn_s3_bgc=tempr

    tString='rf_l1s1_bgc'
    call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))
    params_inst%rf_l1s1_bgc=tempr

    tString='rf_l2s1_bgc'
    call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))
    params_inst%rf_l2s1_bgc=tempr

    tString='rf_l3s2_bgc'
    call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))
    params_inst%rf_l3s2_bgc=tempr   

    tString='rf_s2s1_bgc'
    call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))
    params_inst%rf_s2s1_bgc=tempr

    tString='rf_s2s3_bgc'
    call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))
    params_inst%rf_s2s3_bgc=tempr

    tString='rf_s3s1_bgc'
    call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))
    params_inst%rf_s3s1_bgc=tempr

    tString='rf_cwdl2_bgc'
    call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))
    params_inst%rf_cwdl2_bgc=tempr

    tString='rf_cwdl3_bgc'
    call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))
    params_inst%rf_cwdl3_bgc=tempr

    tString='cwd_fcel'
    call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))
    params_inst%cwd_fcel_bgc=tempr

    tString='k_frag'
    call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))
    params_inst%k_frag_bgc=tempr

    tString='minpsi_hr'
    call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))
    params_inst%minpsi_bgc=tempr 

    tString='maxpsi_hr'
    call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))
    params_inst%maxpsi_bgc=tempr 

    tString='cwd_flig'
    call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))
    params_inst%cwd_flig_bgc=tempr 
    
  end subroutine readParams

  !-----------------------------------------------------------------------
  subroutine init_decompcascade_bgc(bounds, soilbiogeochem_state_inst, soilstate_inst )
    !
    ! !DESCRIPTION:
    !  initialize rate constants and decomposition pathways following the decomposition cascade of the BGC model.
    !  written by C. Koven 
    !
    ! !USES:
    use clm_time_manager , only : get_step_size
    !
    ! !ARGUMENTS:
    type(bounds_type)               , intent(in)    :: bounds  
    type(soilbiogeochem_state_type) , intent(inout) :: soilbiogeochem_state_inst
    type(soilstate_type)            , intent(in)    :: soilstate_inst
    !
    ! !LOCAL VARIABLES
    !-- properties of each decomposing pool
    real(r8) :: rf_l1s1
    real(r8) :: rf_l2s1
    real(r8) :: rf_l3s2
    !real(r8) :: rf_s1s2(bounds%begc:bounds%endc,1:nlevdecomp)
    !real(r8) :: rf_s1s3(bounds%begc:bounds%endc,1:nlevdecomp)
    real(r8), allocatable :: rf_s1s2(:,:)
    real(r8), allocatable :: rf_s1s3(:,:)
    real(r8) :: rf_s2s1
    real(r8) :: rf_s2s3
    real(r8) :: rf_s3s1
    real(r8) :: rf_cwdl2
    real(r8) :: rf_cwdl3
    real(r8) :: cwd_fcel
    real(r8) :: cwd_flig
    real(r8) :: cn_s1
    real(r8) :: cn_s2
    real(r8) :: cn_s3
    !real(r8) :: f_s1s2(bounds%begc:bounds%endc,1:nlevdecomp)
    !real(r8) :: f_s1s3(bounds%begc:bounds%endc,1:nlevdecomp)
    real(r8), allocatable :: f_s1s2(:,:)
    real(r8), allocatable :: f_s1s3(:,:)
    real(r8) :: f_s2s1
    real(r8) :: f_s2s3

    integer :: i_l1s1
    integer :: i_l2s1
    integer :: i_l3s2
    integer :: i_s1s2
    integer :: i_s1s3
    integer :: i_s2s1
    integer :: i_s2s3
    integer :: i_s3s1
    integer :: i_cwdl2
    integer :: i_cwdl3
    real(r8):: speedup_fac                  ! acceleration factor, higher when vertsoilc = .true.

    integer  :: c, j    ! indices
    real(r8) :: t       ! temporary variable
    !-----------------------------------------------------------------------

    associate(                                                                                     &
         rf_decomp_cascade              => soilbiogeochem_state_inst%rf_decomp_cascade_col       , & ! Input:  [real(r8)          (:,:,:) ]  respired fraction in decomposition step (frac)       
         pathfrac_decomp_cascade        => soilbiogeochem_state_inst%pathfrac_decomp_cascade_col , & ! Input:  [real(r8)          (:,:,:) ]  what fraction of C leaving a given pool passes through a given transition (frac)

         cellsand                       => soilstate_inst%cellsand_col                           , & ! Input:  [real(r8)          (:,:)   ]  column 3D sand                                         
         
         cascade_step_name              => decomp_cascade_con%cascade_step_name                  , & ! Output: [character(len=8)  (:)     ]  name of transition                               
         cascade_donor_pool             => decomp_cascade_con%cascade_donor_pool                 , & ! Output: [integer           (:)     ]  which pool is C taken from for a given decomposition step 
         cascade_receiver_pool          => decomp_cascade_con%cascade_receiver_pool              , & ! Output: [integer           (:)     ]  which pool is C added to for a given decomposition step   
         floating_cn_ratio_decomp_pools => decomp_cascade_con%floating_cn_ratio_decomp_pools     , & ! Output: [logical           (:)     ]  TRUE => pool has fixed C:N ratio                          
         decomp_pool_name_restart       => decomp_cascade_con%decomp_pool_name_restart           , & ! Output: [character(len=8)  (:)     ]  name of pool for restart files                   
         decomp_pool_name_history       => decomp_cascade_con%decomp_pool_name_history           , & ! Output: [character(len=8)  (:)     ]  name of pool for history files                   
         decomp_pool_name_long          => decomp_cascade_con%decomp_pool_name_long              , & ! Output: [character(len=20) (:)     ]  name of pool for netcdf long names              
         decomp_pool_name_short         => decomp_cascade_con%decomp_pool_name_short             , & ! Output: [character(len=8)  (:)     ]  name of pool for netcdf short names              
         is_litter                      => decomp_cascade_con%is_litter                          , & ! Output: [logical           (:)     ]  TRUE => pool is a litter pool                             
         is_soil                        => decomp_cascade_con%is_soil                            , & ! Output: [logical           (:)     ]  TRUE => pool is a soil pool                               
         is_cwd                         => decomp_cascade_con%is_cwd                             , & ! Output: [logical           (:)     ]  TRUE => pool is a cwd pool                                
         initial_cn_ratio               => decomp_cascade_con%initial_cn_ratio                   , & ! Output: [real(r8)          (:)     ]  c:n ratio for initialization of pools                    
         initial_stock                  => decomp_cascade_con%initial_stock                      , & ! Output: [real(r8)          (:)     ]  initial concentration for seeding at spinup              
         initial_stock_soildepth        => decomp_cascade_con%initial_stock_soildepth            , & ! Output: [real(r8)          (:)     ]  soil depth for initial concentration for seeding at spinup              
         is_metabolic                   => decomp_cascade_con%is_metabolic                       , & ! Output: [logical           (:)     ]  TRUE => pool is metabolic material                        
         is_cellulose                   => decomp_cascade_con%is_cellulose                       , & ! Output: [logical           (:)     ]  TRUE => pool is cellulose                                 
         is_lignin                      => decomp_cascade_con%is_lignin                          , & ! Output: [logical           (:)     ]  TRUE => pool is lignin                                    
         spinup_factor                  => decomp_cascade_con%spinup_factor                        & ! Output: [real(r8)          (:)     ]  factor for AD spinup associated with each pool           

         )

      allocate(rf_s1s2(bounds%begc:bounds%endc,1:nlevdecomp))
      allocate(rf_s1s3(bounds%begc:bounds%endc,1:nlevdecomp))
      allocate(f_s1s2(bounds%begc:bounds%endc,1:nlevdecomp))
      allocate(f_s1s3(bounds%begc:bounds%endc,1:nlevdecomp))

      !------- time-constant coefficients ---------- !
      ! set soil organic matter compartment C:N ratios
      cn_s1 = params_inst%cn_s1_bgc
      cn_s2 = params_inst%cn_s2_bgc
      cn_s3 = params_inst%cn_s3_bgc

      ! set respiration fractions for fluxes between compartments
      rf_l1s1 = params_inst%rf_l1s1_bgc
      rf_l2s1 = params_inst%rf_l2s1_bgc
      rf_l3s2 = params_inst%rf_l3s2_bgc
      rf_s2s1 = params_inst%rf_s2s1_bgc
      rf_s2s3 = params_inst%rf_s2s3_bgc
      rf_s3s1 = params_inst%rf_s3s1_bgc

      rf_cwdl2 = params_inst%rf_cwdl2_bgc
      rf_cwdl3 = params_inst%rf_cwdl3_bgc

      ! set the cellulose and lignin fractions for coarse woody debris
      cwd_fcel = params_inst%cwd_fcel_bgc
      cwd_flig = params_inst%cwd_flig_bgc

      ! set path fractions
      f_s2s1 = 0.42_r8/(0.45_r8)
      f_s2s3 = 0.03_r8/(0.45_r8)

      ! some of these are dependent on the soil texture properties
      do c = bounds%begc, bounds%endc
         do j = 1, nlevdecomp
            t = 0.85_r8 - 0.68_r8 * 0.01_r8 * (100._r8 - cellsand(c,j))
            f_s1s2(c,j) = 1._r8 - .004_r8 / (1._r8 - t)
            f_s1s3(c,j) = .004_r8 / (1._r8 - t)
            rf_s1s2(c,j) = t
            rf_s1s3(c,j) = t
         end do
      end do
      initial_stock_soildepth = params_inst%initial_Cstocks_depth

      !-------------------  list of pools and their attributes  ------------
      floating_cn_ratio_decomp_pools(i_litr1) = .true.
      decomp_pool_name_restart(i_litr1) = 'litr1'
      decomp_pool_name_history(i_litr1) = 'LITR1'
      decomp_pool_name_long(i_litr1) = 'litter 1'
      decomp_pool_name_short(i_litr1) = 'L1'
      is_litter(i_litr1) = .true.
      is_soil(i_litr1) = .false.
      is_cwd(i_litr1) = .false.
      initial_cn_ratio(i_litr1) = 90._r8
      initial_stock(i_litr1) = 0._r8
      is_metabolic(i_litr1) = .true.
      is_cellulose(i_litr1) = .false.
      is_lignin(i_litr1) = .false.

      floating_cn_ratio_decomp_pools(i_litr2) = .true.
      decomp_pool_name_restart(i_litr2) = 'litr2'
      decomp_pool_name_history(i_litr2) = 'LITR2'
      decomp_pool_name_long(i_litr2) = 'litter 2'
      decomp_pool_name_short(i_litr2) = 'L2'
      is_litter(i_litr2) = .true.
      is_soil(i_litr2) = .false.
      is_cwd(i_litr2) = .false.
      initial_cn_ratio(i_litr2) = 90._r8
      initial_stock(i_litr2) = 0._r8
      is_metabolic(i_litr2) = .false.
      is_cellulose(i_litr2) = .true.
      is_lignin(i_litr2) = .false.

      floating_cn_ratio_decomp_pools(i_litr3) = .true.
      decomp_pool_name_restart(i_litr3) = 'litr3'
      decomp_pool_name_history(i_litr3) = 'LITR3'
      decomp_pool_name_long(i_litr3) = 'litter 3'
      decomp_pool_name_short(i_litr3) = 'L3'
      is_litter(i_litr3) = .true.
      is_soil(i_litr3) = .false.
      is_cwd(i_litr3) = .false.
      initial_cn_ratio(i_litr3) = 90._r8
      initial_stock(i_litr3) = 0._r8
      is_metabolic(i_litr3) = .false.
      is_cellulose(i_litr3) = .false.
      is_lignin(i_litr3) = .true.

      if (.not. use_fates) then
         ! CWD
         floating_cn_ratio_decomp_pools(i_cwd) = .true.
         decomp_pool_name_restart(i_cwd) = 'cwd'
         decomp_pool_name_history(i_cwd) = 'CWD'
         decomp_pool_name_long(i_cwd) = 'coarse woody debris'
         decomp_pool_name_short(i_cwd) = 'CWD'
         is_litter(i_cwd) = .false.
         is_soil(i_cwd) = .false.
         is_cwd(i_cwd) = .true.
         initial_cn_ratio(i_cwd) = 90._r8
         initial_stock(i_cwd) = 0._r8
         is_metabolic(i_cwd) = .false.
         is_cellulose(i_cwd) = .false.
         is_lignin(i_cwd) = .false.
      endif

      if (.not. use_fates) then
         i_soil1 = 5
      else
         i_soil1 = 4
      endif
      floating_cn_ratio_decomp_pools(i_soil1) = .false.
      decomp_pool_name_restart(i_soil1) = 'soil1'
      decomp_pool_name_history(i_soil1) = 'SOIL1'
      decomp_pool_name_long(i_soil1) = 'soil 1'
      decomp_pool_name_short(i_soil1) = 'S1'
      is_litter(i_soil1) = .false.
      is_soil(i_soil1) = .true.
      is_cwd(i_soil1) = .false.
      initial_cn_ratio(i_soil1) = cn_s1
      initial_stock(i_soil1) = params_inst%initial_Cstocks(1)
      is_metabolic(i_soil1) = .false.
      is_cellulose(i_soil1) = .false.
      is_lignin(i_soil1) = .false.

      if (.not. use_fates) then
         i_soil2 = 6
      else
         i_soil2 = 5
      endif
      floating_cn_ratio_decomp_pools(i_soil2) = .false.
      decomp_pool_name_restart(i_soil2) = 'soil2'
      decomp_pool_name_history(i_soil2) = 'SOIL2'
      decomp_pool_name_long(i_soil2) = 'soil 2'
      decomp_pool_name_short(i_soil2) = 'S2'
      is_litter(i_soil2) = .false.
      is_soil(i_soil2) = .true.
      is_cwd(i_soil2) = .false.
      initial_cn_ratio(i_soil2) = cn_s2
      initial_stock(i_soil2) = params_inst%initial_Cstocks(2)
      is_metabolic(i_soil2) = .false.
      is_cellulose(i_soil2) = .false.
      is_lignin(i_soil2) = .false.

      if (.not. use_fates) then
         i_soil3 = 7
      else
         i_soil3 = 6
      endif
      floating_cn_ratio_decomp_pools(i_soil3) = .false.
      decomp_pool_name_restart(i_soil3) = 'soil3'
      decomp_pool_name_history(i_soil3) = 'SOIL3'
      decomp_pool_name_long(i_soil3) = 'soil 3'
      decomp_pool_name_short(i_soil3) = 'S3'
      is_litter(i_soil3) = .false.
      is_soil(i_soil3) = .true.
      is_cwd(i_soil3) = .false.
      initial_cn_ratio(i_soil3) = cn_s3
      initial_stock(i_soil3) = params_inst%initial_Cstocks(3)
      is_metabolic(i_soil3) = .false.
      is_cellulose(i_soil3) = .false.
      is_lignin(i_soil3) = .false.


      speedup_fac = 1._r8

      !lit1
      spinup_factor(i_litr1) = 1._r8
      !lit2,3
      spinup_factor(i_litr2) = 1._r8
      spinup_factor(i_litr3) = 1._r8
      !CWD
      if (.not. use_fates) then
         spinup_factor(i_cwd) = max(1._r8, (speedup_fac * params_inst%tau_cwd_bgc / 2._r8 ))
      end if
      !som1
      spinup_factor(i_soil1) = 1._r8
      !som2,3
      spinup_factor(i_soil2) = max(1._r8, (speedup_fac * params_inst%tau_s2_bgc))
      spinup_factor(i_soil3) = max(1._r8, (speedup_fac * params_inst%tau_s3_bgc))

      if ( masterproc ) then
         write(iulog,*) 'Spinup_state ',spinup_state
         write(iulog,*) 'Spinup factors ',spinup_factor
      end if

      !----------------  list of transitions and their time-independent coefficients  ---------------!
      i_l1s1 = 1
      cascade_step_name(i_l1s1) = 'L1S1'
      rf_decomp_cascade(bounds%begc:bounds%endc,1:nlevdecomp,i_l1s1) = rf_l1s1
      cascade_donor_pool(i_l1s1) = i_litr1
      cascade_receiver_pool(i_l1s1) = i_soil1
      pathfrac_decomp_cascade(bounds%begc:bounds%endc,1:nlevdecomp,i_l1s1) = 1.0_r8

      i_l2s1 = 2
      cascade_step_name(i_l2s1) = 'L2S1'
      rf_decomp_cascade(bounds%begc:bounds%endc,1:nlevdecomp,i_l2s1) = rf_l2s1
      cascade_donor_pool(i_l2s1) = i_litr2
      cascade_receiver_pool(i_l2s1) = i_soil1
      pathfrac_decomp_cascade(bounds%begc:bounds%endc,1:nlevdecomp,i_l2s1)= 1.0_r8

      i_l3s2 = 3
      cascade_step_name(i_l3s2) = 'L3S2'
      rf_decomp_cascade(bounds%begc:bounds%endc,1:nlevdecomp,i_l3s2) = rf_l3s2
      cascade_donor_pool(i_l3s2) = i_litr3
      cascade_receiver_pool(i_l3s2) = i_soil2
      pathfrac_decomp_cascade(bounds%begc:bounds%endc,1:nlevdecomp,i_l3s2) = 1.0_r8

      i_s1s2 = 4
      cascade_step_name(i_s1s2) = 'S1S2'
      rf_decomp_cascade(bounds%begc:bounds%endc,1:nlevdecomp,i_s1s2) = rf_s1s2(bounds%begc:bounds%endc,1:nlevdecomp)
      cascade_donor_pool(i_s1s2) = i_soil1
      cascade_receiver_pool(i_s1s2) = i_soil2
      pathfrac_decomp_cascade(bounds%begc:bounds%endc,1:nlevdecomp,i_s1s2) = f_s1s2(bounds%begc:bounds%endc,1:nlevdecomp)

      i_s1s3 = 5
      cascade_step_name(i_s1s3) = 'S1S3'
      rf_decomp_cascade(bounds%begc:bounds%endc,1:nlevdecomp,i_s1s3) = rf_s1s3(bounds%begc:bounds%endc,1:nlevdecomp)
      cascade_donor_pool(i_s1s3) = i_soil1
      cascade_receiver_pool(i_s1s3) = i_soil3
      pathfrac_decomp_cascade(bounds%begc:bounds%endc,1:nlevdecomp,i_s1s3) = f_s1s3(bounds%begc:bounds%endc,1:nlevdecomp)

      i_s2s1 = 6
      cascade_step_name(i_s2s1) = 'S2S1'
      rf_decomp_cascade(bounds%begc:bounds%endc,1:nlevdecomp,i_s2s1) = rf_s2s1
      cascade_donor_pool(i_s2s1) = i_soil2
      cascade_receiver_pool(i_s2s1) = i_soil1
      pathfrac_decomp_cascade(bounds%begc:bounds%endc,1:nlevdecomp,i_s2s1) = f_s2s1

      i_s2s3 = 7 
      cascade_step_name(i_s2s3) = 'S2S3'
      rf_decomp_cascade(bounds%begc:bounds%endc,1:nlevdecomp,i_s2s3) = rf_s2s3
      cascade_donor_pool(i_s2s3) = i_soil2
      cascade_receiver_pool(i_s2s3) = i_soil3
      pathfrac_decomp_cascade(bounds%begc:bounds%endc,1:nlevdecomp,i_s2s3) = f_s2s3

      i_s3s1 = 8
      cascade_step_name(i_s3s1) = 'S3S1'
      rf_decomp_cascade(bounds%begc:bounds%endc,1:nlevdecomp,i_s3s1) = rf_s3s1
      cascade_donor_pool(i_s3s1) = i_soil3
      cascade_receiver_pool(i_s3s1) = i_soil1
      pathfrac_decomp_cascade(bounds%begc:bounds%endc,1:nlevdecomp,i_s3s1) = 1.0_r8

      if (.not. use_fates) then
         i_cwdl2 = 9
         cascade_step_name(i_cwdl2) = 'CWDL2'
         rf_decomp_cascade(bounds%begc:bounds%endc,1:nlevdecomp,i_cwdl2) = rf_cwdl2
         cascade_donor_pool(i_cwdl2) = i_cwd
         cascade_receiver_pool(i_cwdl2) = i_litr2
         pathfrac_decomp_cascade(bounds%begc:bounds%endc,1:nlevdecomp,i_cwdl2) = cwd_fcel
         
         i_cwdl3 = 10
         cascade_step_name(i_cwdl3) = 'CWDL3'
         rf_decomp_cascade(bounds%begc:bounds%endc,1:nlevdecomp,i_cwdl3) = rf_cwdl3
         cascade_donor_pool(i_cwdl3) = i_cwd
         cascade_receiver_pool(i_cwdl3) = i_litr3
         pathfrac_decomp_cascade(bounds%begc:bounds%endc,1:nlevdecomp,i_cwdl3) = cwd_flig
      end if

      deallocate(rf_s1s2)
      deallocate(rf_s1s3)
      deallocate(f_s1s2)
      deallocate(f_s1s3)

    end associate

  end subroutine init_decompcascade_bgc

end module SoilBiogeochemDecompCascadeBGCMod
