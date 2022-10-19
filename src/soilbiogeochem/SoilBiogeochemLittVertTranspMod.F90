module SoilBiogeochemLittVertTranspMod

  !-----------------------------------------------------------------------
  ! calculate vertical mixing of all decomposing C and N pools
  !
  use shr_kind_mod                       , only : r8 => shr_kind_r8
  use shr_log_mod                        , only : errMsg => shr_log_errMsg
  use clm_varctl                         , only : iulog, spinup_state, use_cn
  use clm_varcon                         , only : secspday
  use decompMod                          , only : bounds_type
  use abortutils                         , only : endrun
  use CanopyStateType                    , only : canopystate_type
  use SoilBiogeochemStateType            , only : soilbiogeochem_state_type
  use SoilBiogeochemCarbonFluxType       , only : soilbiogeochem_carbonflux_type
  use SoilBiogeochemCarbonStateType      , only : soilbiogeochem_carbonstate_type
  use SoilBiogeochemNitrogenFluxType     , only : soilbiogeochem_nitrogenflux_type
  use SoilBiogeochemNitrogenStateType    , only : soilbiogeochem_nitrogenstate_type
  use SoilBiogeochemDecompCascadeConType , only : decomp_cascade_con
  use ColumnType                         , only : col                
  use GridcellType                       , only : grc
  use SoilBiogeochemStateType            , only : get_spinup_latitude_term
  !
  implicit none
  private
  !
  public :: readParams
  public :: SoilBiogeochemLittVertTransp

  type, private :: params_type
     real(r8) :: som_diffus                 ! Soil organic matter diffusion
     real(r8) :: cryoturb_diffusion_k       ! The cryoturbation diffusive constant cryoturbation to the active layer thickness
     real(r8) :: max_altdepth_cryoturbation ! (m) maximum active layer thickness for cryoturbation to occur
  end type params_type

  type(params_type), private :: params_inst
  !
  real(r8), public :: som_adv_flux =  0._r8
  real(r8), public :: max_depth_cryoturb = 3._r8   ! (m) this is the maximum depth of cryoturbation
  real(r8) :: som_diffus                   ! [m^2/sec] = 1 cm^2 / yr
  real(r8) :: cryoturb_diffusion_k         ! [m^2/sec] = 5 cm^2 / yr = 1m^2 / 200 yr
  real(r8) :: max_altdepth_cryoturbation   ! (m) maximum active layer thickness for cryoturbation to occur

  character(len=*), parameter, private :: sourcefile = &
       __FILE__
  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------  
  subroutine readParams ( ncid )
    !
    use ncdio_pio   , only : file_desc_t,ncd_io
    !
    type(file_desc_t),intent(inout) :: ncid   ! pio netCDF file id
    !
    character(len=32)  :: subname = 'SoilBiogeochemLittVertTranspType'
    character(len=100) :: errCode = '-Error reading in parameters file:'
    logical            :: readv ! has variable been read in or not
    real(r8)           :: tempr ! temporary to read in constant
    character(len=100) :: tString ! temp. var for reading
    !-----------------------------------------------------------------------
    !
    ! read in parameters
    !

     tString='som_diffus'
     call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
     if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))
     !soilbiogeochem_litt_verttransp_params_inst%som_diffus=tempr
     ! FIX(SPM,032414) - can't be pulled out since division makes things not bfb
     params_inst%som_diffus = 1e-4_r8 / (secspday * 365._r8)  

     tString='cryoturb_diffusion_k'
     call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
     if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))
     !soilbiogeochem_litt_verttransp_params_inst%cryoturb_diffusion_k=tempr
     !FIX(SPM,032414) Todo.  This constant cannot be on file since the divide makes things
     !SPM Todo.  This constant cannot be on file since the divide makes things
     !not bfb
     params_inst%cryoturb_diffusion_k = 5e-4_r8 / (secspday * 365._r8)  ! [m^2/sec] = 5 cm^2 / yr = 1m^2 / 200 yr

     tString='max_altdepth_cryoturbation'
     call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
     if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))
     params_inst%max_altdepth_cryoturbation=tempr
    
   end subroutine readParams

  !-----------------------------------------------------------------------
  subroutine SoilBiogeochemLittVertTransp(bounds, num_soilc, filter_soilc,      &
       canopystate_inst, soilbiogeochem_state_inst,                     &
       soilbiogeochem_carbonstate_inst, soilbiogeochem_carbonflux_inst, &
       c13_soilbiogeochem_carbonstate_inst, c13_soilbiogeochem_carbonflux_inst, &
       c14_soilbiogeochem_carbonstate_inst, c14_soilbiogeochem_carbonflux_inst, &
       soilbiogeochem_nitrogenstate_inst, soilbiogeochem_nitrogenflux_inst)
    !
    ! !DESCRIPTION:
    ! Calculate vertical mixing of soil and litter pools.  Also reconcile sources and sinks of these pools 
    ! calculated in the CStateUpdate1 and NStateUpdate1 subroutines.
    ! Advection-diffusion code based on algorithm in Patankar (1980)
    ! Initial code by C. Koven and W. Riley
    !
    ! !USES:
    use clm_time_manager , only : get_step_size
    use clm_varpar       , only : nlevdecomp, ndecomp_pools, nlevdecomp_full
    use clm_varcon       , only : zsoi, dzsoi_decomp, zisoi
    use TridiagonalMod   , only : Tridiagonal
    use ColumnType       , only : col
    use clm_varctl       , only : use_bedrock

    !
    ! !ARGUMENTS:
    type(bounds_type)                       , intent(in)    :: bounds 
    integer                                 , intent(in)    :: num_soilc        ! number of soil columns in filter
    integer                                 , intent(in)    :: filter_soilc(:)  ! filter for soil columns
    type(canopystate_type)                  , intent(in)    :: canopystate_inst
    type(soilbiogeochem_state_type)         , intent(inout) :: soilbiogeochem_state_inst
    type(soilbiogeochem_carbonstate_type)   , intent(inout) :: soilbiogeochem_carbonstate_inst
    type(soilbiogeochem_carbonflux_type)    , intent(inout) :: soilbiogeochem_carbonflux_inst
    type(soilbiogeochem_carbonstate_type)   , intent(inout) :: c13_soilbiogeochem_carbonstate_inst
    type(soilbiogeochem_carbonflux_type)    , intent(inout) :: c13_soilbiogeochem_carbonflux_inst
    type(soilbiogeochem_carbonstate_type)   , intent(inout) :: c14_soilbiogeochem_carbonstate_inst
    type(soilbiogeochem_carbonflux_type)    , intent(inout) :: c14_soilbiogeochem_carbonflux_inst
    type(soilbiogeochem_nitrogenstate_type) , intent(inout) :: soilbiogeochem_nitrogenstate_inst
    type(soilbiogeochem_nitrogenflux_type)  , intent(inout) :: soilbiogeochem_nitrogenflux_inst
    !
    ! !LOCAL VARIABLES:
    real(r8) :: diffus (bounds%begc:bounds%endc,1:nlevdecomp+1)                    ! diffusivity (m2/s)  (includes spinup correction, if any)
    real(r8) :: adv_flux(bounds%begc:bounds%endc,1:nlevdecomp+1)                   ! advective flux (m/s)  (includes spinup correction, if any)
    real(r8) :: aaa                                                                ! "A" function in Patankar
    real(r8) :: pe                                                                 ! Pe for "A" function in Patankar
    real(r8) :: w_m1, w_p1                                                         ! Weights for calculating harmonic mean of diffusivity
    real(r8) :: d_m1, d_p1                                                         ! Harmonic mean of diffusivity
    real(r8) :: a_tri(bounds%begc:bounds%endc,0:nlevdecomp+1)                      ! "a" vector for tridiagonal matrix
    real(r8) :: b_tri(bounds%begc:bounds%endc,0:nlevdecomp+1)                      ! "b" vector for tridiagonal matrix
    real(r8) :: c_tri(bounds%begc:bounds%endc,0:nlevdecomp+1)                      ! "c" vector for tridiagonal matrix
    real(r8) :: r_tri(bounds%begc:bounds%endc,0:nlevdecomp+1)                      ! "r" vector for tridiagonal solution
    real(r8) :: d_p1_zp1(bounds%begc:bounds%endc,1:nlevdecomp+1)                   ! diffusivity/delta_z for next j  (set to zero for no diffusion)
    real(r8) :: d_m1_zm1(bounds%begc:bounds%endc,1:nlevdecomp+1)                   ! diffusivity/delta_z for previous j (set to zero for no diffusion)
    real(r8) :: f_p1(bounds%begc:bounds%endc,1:nlevdecomp+1)                       ! water flux for next j
    real(r8) :: f_m1(bounds%begc:bounds%endc,1:nlevdecomp+1)                       ! water flux for previous j
    real(r8) :: pe_p1(bounds%begc:bounds%endc,1:nlevdecomp+1)                      ! Peclet # for next j
    real(r8) :: pe_m1(bounds%begc:bounds%endc,1:nlevdecomp+1)                      ! Peclet # for previous j
    real(r8) :: dz_node(1:nlevdecomp+1)                                            ! difference between nodes
    real(r8) :: epsilon_t (bounds%begc:bounds%endc,1:nlevdecomp+1,1:ndecomp_pools) !
    real(r8) :: conc_trcr(bounds%begc:bounds%endc,0:nlevdecomp+1)                  !
    real(r8) :: a_p_0
    real(r8) :: deficit
    integer  :: ntype
    integer  :: i_type,s,fc,c,j,l                                                  ! indices
    integer  :: jtop(bounds%begc:bounds%endc)                                      ! top level at each column
    real(r8) :: dtime                                                              ! land model time step (sec)
    integer  :: zerolev_diffus
    real(r8) :: spinup_term                                                        ! spinup accelerated decomposition factor, used to accelerate transport as well
    real(r8) :: epsilon                                                            ! small number
    real(r8), pointer :: conc_ptr(:,:,:)                                           ! pointer, concentration state variable being transported
    real(r8), pointer :: source(:,:,:)                                             ! pointer, source term
    real(r8), pointer :: trcr_tendency_ptr(:,:,:)                                  ! poiner, store the vertical tendency (gain/loss due to vertical transport)
    !-----------------------------------------------------------------------

    ! Set statement functions
    aaa (pe) = max (0._r8, (1._r8 - 0.1_r8 * abs(pe))**5)  ! A function from Patankar, Table 5.2, pg 95
  
    associate(                                                             &
         is_cwd           => decomp_cascade_con%is_cwd                  ,  & ! Input:  [logical (:)    ]  TRUE => pool is a cwd pool                                
         spinup_factor    => decomp_cascade_con%spinup_factor           ,  & ! Input:  [real(r8) (:)   ]  spinup accelerated decomposition factor, used to accelerate transport as well

         altmax           => canopystate_inst%altmax_col                ,  & ! Input:  [real(r8) (:)   ]  maximum annual depth of thaw                             
         altmax_lastyear  => canopystate_inst%altmax_lastyear_col       ,  & ! Input:  [real(r8) (:)   ]  prior year maximum annual depth of thaw                  

         som_adv_coef     => soilbiogeochem_state_inst%som_adv_coef_col ,  & ! Output: [real(r8) (:,:) ]  SOM advective flux (m/s)                               
         som_diffus_coef  => soilbiogeochem_state_inst%som_diffus_coef_col & ! Output: [real(r8) (:,:) ]  SOM diffusivity due to bio/cryo-turbation (m2/s)       
         )

      !Set parameters of vertical mixing of SOM
      som_diffus                 = params_inst%som_diffus 
      cryoturb_diffusion_k       = params_inst%cryoturb_diffusion_k 
      max_altdepth_cryoturbation = params_inst%max_altdepth_cryoturbation 

      dtime = get_step_size()

      ntype = 2
      spinup_term = 1._r8
      epsilon = 1.e-30

      !------ loop over litter/som types
      do i_type = 1, ntype

         select case (i_type)
         case (1)  ! C
            conc_ptr          => soilbiogeochem_carbonstate_inst%decomp_cpools_vr_col
            source            => soilbiogeochem_carbonflux_inst%decomp_cpools_sourcesink_col
            trcr_tendency_ptr => soilbiogeochem_carbonflux_inst%decomp_cpools_transport_tendency_col
         case (2)  ! N
            if (use_cn ) then
               conc_ptr          => soilbiogeochem_nitrogenstate_inst%decomp_npools_vr_col
               source            => soilbiogeochem_nitrogenflux_inst%decomp_npools_sourcesink_col
               trcr_tendency_ptr => soilbiogeochem_nitrogenflux_inst%decomp_npools_transport_tendency_col
            endif
         case (3)
            write(iulog,*) 'error.  ncase = 4, but c13 and c14 not both enabled.'
            call endrun(msg=errMsg(sourcefile, __LINE__))
         case (4)
            write(iulog,*) 'error.  ncase = 4, but c13 and c14 not both enabled.'
            call endrun(msg=errMsg(sourcefile, __LINE__))
         end select

         !! for single level case, no transport; just update the fluxes calculated in the StateUpdate1 subroutines
         do l = 1, ndecomp_pools
            do j = 1,nlevdecomp
               do fc = 1, num_soilc
                  c = filter_soilc (fc)
                  conc_ptr(c,j,l) = conc_ptr(c,j,l) + source(c,j,l)
                  trcr_tendency_ptr(c,j,l) = 0._r8
               end do
            end do
         end do


      end do  ! i_type
   
    end associate

  end subroutine SoilBiogeochemLittVertTransp
 
end module SoilBiogeochemLittVertTranspMod
