module SoilStateType

  !------------------------------------------------------------------------------
  ! !USES:
  use shr_kind_mod    , only : r8 => shr_kind_r8
  use shr_log_mod     , only : errMsg => shr_log_errMsg
  use shr_infnan_mod  , only : nan => shr_infnan_nan, assignment(=)
  use decompMod       , only : bounds_type
  use abortutils      , only : endrun
  use clm_varpar      , only : nlevsoi, nlevgrnd, nlevlak, nlayer, nlevsno
  use clm_varcon      , only : spval
  use clm_varctl      , only : iulog
  use LandunitType    , only : lun                
  use ColumnType      , only : col                
  use PatchType       , only : patch                
  !
  implicit none
  save
  private
  !
  ! !PUBLIC TYPES:
  type, public :: soilstate_type

     ! hydraulic properties
     real(r8), pointer :: hk_l_col             (:,:) ! col hydraulic conductivity (mm/s)
     real(r8), pointer :: smp_l_col            (:,:) ! col soil matric potential (mm)
     real(r8), pointer :: smpmin_col           (:)   ! col restriction for min of soil potential (mm) 
     real(r8), pointer :: watsat_col           (:,:) ! col volumetric soil water at saturation (porosity) 
     real(r8), pointer :: dsl_col              (:)   ! col dry surface layer thickness (mm)
     real(r8), pointer :: soilresis_col        (:)   ! col soil evaporative resistance S&L14 (s/m)
     real(r8), pointer :: soilbeta_col         (:)   ! col factor that reduces ground evaporation L&P1992(-)
     real(r8), pointer :: soilalpha_col        (:)   ! col factor that reduces ground saturated specific humidity (-)
     real(r8), pointer :: soilpsi_col          (:,:) ! col soil water potential in each soil layer (MPa) (CN)
     real(r8), pointer :: wtfact_col           (:)   ! col maximum saturated fraction for a gridcell
     real(r8), pointer :: msw_col              (:,:) ! col vanGenuchtenClapp "m"
     real(r8), pointer :: nsw_col              (:,:) ! col vanGenuchtenClapp "n"
     real(r8), pointer :: alphasw_col          (:,:) ! col vanGenuchtenClapp "nalpha"
     real(r8), pointer :: watres_col           (:,:) ! residual soil water content
     ! thermal conductivity / heat capacity
     real(r8), pointer :: thk_col              (:,:) ! col thermal conductivity of each layer [W/m-K] 

     ! roots
     real(r8), pointer :: rootr_patch          (:,:) ! patch effective fraction of roots in each soil layer (nlevgrnd)
     real(r8), pointer :: rootr_col            (:,:) ! col effective fraction of roots in each soil layer (nlevgrnd)  
     real(r8), pointer :: rootfr_col           (:,:) ! col fraction of roots in each soil layer (nlevgrnd) 
     real(r8), pointer :: rootfr_patch         (:,:) ! patch fraction of roots for water in each soil layer (nlevgrnd)
     real(r8), pointer :: crootfr_patch        (:,:) ! patch fraction of roots for carbon in each soil layer (nlevgrnd)
     real(r8), pointer :: root_depth_patch     (:)   ! root depth
     real(r8), pointer :: k_soil_root_patch    (:,:) ! patch soil-root interface conductance [mm/s]
     real(r8), pointer :: root_conductance_patch(:,:) ! patch root conductance [mm/s]
     real(r8), pointer :: soil_conductance_patch(:,:) ! patch soil conductance [mm/s]

   contains

     procedure, public  :: Init         
     procedure, public  :: Restart
     procedure, private :: InitAllocate 
     procedure, private :: InitHistory  
     procedure, private :: InitCold     

  end type soilstate_type
  !------------------------------------------------------------------------

contains

  !------------------------------------------------------------------------
  subroutine Init(this, bounds)

    class(soilstate_type) :: this
    type(bounds_type), intent(in) :: bounds  

    call this%InitAllocate(bounds)
    call this%InitHistory(bounds)
    call this%InitCold(bounds)

  end subroutine Init

  !------------------------------------------------------------------------
  subroutine InitAllocate(this, bounds)
    !
    ! !DESCRIPTION:
    ! Initialize module data structure
    !
    ! !ARGUMENTS:
    class(soilstate_type) :: this
    type(bounds_type), intent(in) :: bounds  
    !
    ! !LOCAL VARIABLES:
    integer :: begp, endp
    integer :: begc, endc
    integer :: begg, endg
    !------------------------------------------------------------------------

    begp = bounds%begp; endp= bounds%endp
    begc = bounds%begc; endc= bounds%endc
    begg = bounds%begg; endg= bounds%endg

    allocate(this%hk_l_col             (begc:endc,nlevgrnd))            ; this%hk_l_col             (:,:) = nan   
    allocate(this%smp_l_col            (begc:endc,nlevgrnd))            ; this%smp_l_col            (:,:) = nan   
    allocate(this%smpmin_col           (begc:endc))                     ; this%smpmin_col           (:)   = nan

    allocate(this%watsat_col           (begc:endc,nlevgrnd))            ; this%watsat_col           (:,:) = nan
    allocate(this%dsl_col              (begc:endc))                     ; this%dsl_col         (:)   = spval!nan   
    allocate(this%soilresis_col        (begc:endc))                     ; this%soilresis_col         (:)   = spval!nan   
    allocate(this%soilbeta_col         (begc:endc))                     ; this%soilbeta_col         (:)   = nan   
    allocate(this%soilalpha_col        (begc:endc))                     ; this%soilalpha_col        (:)   = nan
    allocate(this%soilpsi_col          (begc:endc,nlevgrnd))            ; this%soilpsi_col          (:,:) = nan
    allocate(this%wtfact_col           (begc:endc))                     ; this%wtfact_col           (:)   = nan

    allocate(this%thk_col              (begc:endc,-nlevsno+1:nlevgrnd)) ; this%thk_col              (:,:) = nan

    allocate(this%rootr_patch          (begp:endp,1:nlevgrnd))          ; this%rootr_patch          (:,:) = nan
    allocate(this%root_depth_patch     (begp:endp))                     ; this%root_depth_patch     (:)   = nan
    allocate(this%rootr_col            (begc:endc,nlevgrnd))            ; this%rootr_col            (:,:) = nan
    allocate(this%rootfr_patch         (begp:endp,1:nlevgrnd))          ; this%rootfr_patch         (:,:) = nan
    allocate(this%crootfr_patch        (begp:endp,1:nlevgrnd))          ; this%crootfr_patch        (:,:) = nan
    allocate(this%rootfr_col           (begc:endc,1:nlevgrnd))          ; this%rootfr_col           (:,:) = nan 
    allocate(this%k_soil_root_patch    (begp:endp,1:nlevsoi))           ; this%k_soil_root_patch (:,:) = nan
    allocate(this%root_conductance_patch(begp:endp,1:nlevsoi))          ; this%root_conductance_patch (:,:) = nan
    allocate(this%soil_conductance_patch(begp:endp,1:nlevsoi))          ; this%soil_conductance_patch (:,:) = nan
    allocate(this%msw_col              (begc:endc,1:nlevgrnd))          ; this%msw_col              (:,:) = nan
    allocate(this%nsw_col              (begc:endc,1:nlevgrnd))          ; this%nsw_col              (:,:) = nan
    allocate(this%alphasw_col          (begc:endc,1:nlevgrnd))          ; this%alphasw_col          (:,:) = nan
    allocate(this%watres_col           (begc:endc,1:nlevgrnd))          ; this%watres_col           (:,:) = nan
  end subroutine InitAllocate

  !-----------------------------------------------------------------------
  subroutine InitHistory(this, bounds)
    !
    ! History fields initialization
    !
    ! !USES:
    use histFileMod   , only: hist_addfld1d, hist_addfld2d, no_snow_normal
    !
    ! !ARGUMENTS:
    class(soilstate_type) :: this
    type(bounds_type), intent(in) :: bounds  
    !
    ! !LOCAL VARIABLES:
    integer :: begc, endc
    integer :: begp, endp
    character(10)     :: active
    real(r8), pointer :: data2dptr(:,:), data1dptr(:) ! temp. pointers for slicing larger arrays
    !---------------------------------------------------------------------

    begp = bounds%begp; endp= bounds%endp
    begc = bounds%begc; endc= bounds%endc


    active = "inactive"

    call hist_addfld2d (fname='SMP',  units='mm', type2d='levgrnd',  &
         avgflag='A', long_name='soil matric potential (vegetated landunits only)', &
         ptr_col=this%smp_l_col, set_spec=spval, l2g_scale_type='veg', default='inactive')

       this%soil_conductance_patch(begp:endp,:) = spval
       call hist_addfld2d (fname='KSOIL', units='1/s', type2d='levsoi', &
          avgflag='A', long_name='soil conductance in each soil layer', &
          ptr_patch=this%soil_conductance_patch, default='inactive')

    this%thk_col(begc:endc,-nlevsno+1:0) = spval
    data2dptr => this%thk_col(:,-nlevsno+1:0)
    call hist_addfld2d (fname='SNO_TK', units='W/m-K', type2d='levsno', &
         avgflag='A', long_name='Thermal conductivity', &
         ptr_col=data2dptr, no_snow_behavior=no_snow_normal, default='inactive')

    call hist_addfld2d (fname='SNO_TK_ICE', units='W/m-K', type2d='levsno', &
         avgflag='A', long_name='Thermal conductivity (ice landunits only)', &
         ptr_col=data2dptr, no_snow_behavior=no_snow_normal, &
         l2g_scale_type='ice', default='inactive')

    this%hk_l_col(begc:endc,:) = spval
    call hist_addfld2d (fname='HK',  units='mm/s', type2d='levgrnd',  &
         avgflag='A', long_name='hydraulic conductivity (vegetated landunits only)', &
         ptr_col=this%hk_l_col, set_spec=spval, l2g_scale_type='veg', default='inactive')

    this%soilalpha_col(begc:endc) = spval
    call hist_addfld1d (fname='SoilAlpha',  units='unitless',  &
         avgflag='A', long_name='factor limiting ground evap', &
         ptr_col=this%soilalpha_col, set_urb=spval, default='inactive' )

    this%soilresis_col(begc:endc) = spval
    call hist_addfld1d (fname='SOILRESIS',  units='s/m',  &
         avgflag='A', long_name='soil resistance to evaporation', &
         ptr_col=this%soilresis_col, default='inactive')

  end subroutine InitHistory

  !-----------------------------------------------------------------------
  subroutine InitCold(this, bounds)
    !
    ! Initialize module soil state variables to reasonable values
    !
    ! !USES:
    use clm_varpar      , only : nlevgrnd
    !
    ! !ARGUMENTS:
    class(soilstate_type) :: this
    type(bounds_type), intent(in) :: bounds  
    !
    ! !LOCAL VARIABLES:
    !-----------------------------------------------------------------------

    this%smp_l_col(bounds%begc:bounds%endc,1:nlevgrnd) = -1000._r8
    this%hk_l_col(bounds%begc:bounds%endc,1:nlevgrnd) = 0._r8

  end subroutine InitCold

  !------------------------------------------------------------------------
  subroutine Restart(this, bounds, ncid, flag)
    ! 
    ! !DESCRIPTION:
    ! Read/Write module information to/from restart file.
    !
    ! !USES:
    use ncdio_pio        , only : file_desc_t, ncd_io, ncd_double
    use restUtilMod
    use spmdMod          , only : masterproc
    !
    ! !ARGUMENTS:
    class(soilstate_type) :: this
    type(bounds_type), intent(in)    :: bounds 
    type(file_desc_t), intent(inout) :: ncid   ! netcdf id
    character(len=*) , intent(in)    :: flag   ! 'read' or 'write'
    !
    ! !LOCAL VARIABLES:
    integer  :: c
    logical  :: readvar
    logical  :: readrootfr = .false.
    !------------------------------------------------------------------------

    call restartvar(ncid=ncid, flag=flag, varname='DSL', xtype=ncd_double,  &
         dim1name='column', long_name='dsl thickness', units='mm', &
         interpinic_flag='interp', readvar=readvar, data=this%dsl_col)
    
    call restartvar(ncid=ncid, flag=flag, varname='SOILRESIS', xtype=ncd_double,  &
         dim1name='column', long_name='soil resistance', units='s/m', &
         interpinic_flag='interp', readvar=readvar, data=this%soilresis_col)

    call restartvar(ncid=ncid, flag=flag, varname='SMP', xtype=ncd_double,  &
         dim1name='column', dim2name='levgrnd', switchdim=.true., &
         long_name='soil matric potential', units='mm', &
         interpinic_flag='interp', readvar=readvar, data=this%smp_l_col)

    call restartvar(ncid=ncid, flag=flag, varname='HK', xtype=ncd_double,  &
         dim1name='column', dim2name='levgrnd', switchdim=.true., &
         long_name='hydraulic conductivity', units='mm/s', &
         interpinic_flag='interp', readvar=readvar, data=this%hk_l_col)

  end subroutine Restart

end module SoilStateType
