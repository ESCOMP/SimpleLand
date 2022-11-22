module SolarAbsorbedType

  !------------------------------------------------------------------------------
  ! !USES:
  use shr_kind_mod , only: r8 => shr_kind_r8
  use shr_log_mod  , only: errMsg => shr_log_errMsg
  use decompMod    , only : bounds_type
  use clm_varcon   , only : spval
  !
  ! !PUBLIC TYPES:
  implicit none
  save
  private
  !
  ! !PUBLIC DATA MEMBERS:
  type, public :: solarabs_type

     ! Solar reflected
     real(r8), pointer :: fsr_patch              (:)   ! patch solar radiation reflected (W/m**2)         
     
     ! Solar Absorbed
     real(r8), pointer :: fsa_patch              (:)   ! patch solar radiation absorbed (total) (W/m**2)  
     real(r8), pointer :: parsun_z_patch         (:,:) ! patch absorbed PAR for sunlit leaves in canopy layer (W/m**2) 
     real(r8), pointer :: parsha_z_patch         (:,:) ! patch absorbed PAR for shaded leaves in canopy layer (W/m**2) 
     real(r8), pointer :: par240d_z_patch        (:,:) ! 10-day running mean of daytime patch absorbed PAR for leaves in canopy layer (W/m**2) 
     real(r8), pointer :: par240x_z_patch        (:,:) ! 10-day running mean of maximum patch absorbed PAR for leaves in canopy layer (W/m**2)
     real(r8), pointer :: par24d_z_patch         (:,:) ! daily accumulated  absorbed PAR for leaves in canopy layer from midnight to current step(J/m**2) 
     real(r8), pointer :: par24x_z_patch         (:,:) ! daily max of patch absorbed PAR for  leaves in canopy layer from midnight to current step(W/m**2)
     real(r8), pointer :: sabg_soil_patch        (:)   ! patch solar radiation absorbed by soil (W/m**2)       
     real(r8), pointer :: sabg_snow_patch        (:)   ! patch solar radiation absorbed by snow (W/m**2)       
     real(r8), pointer :: sabg_patch             (:)   ! patch solar radiation absorbed by ground (W/m**2)     
     real(r8), pointer :: sabg_chk_patch         (:)   ! patch fsno weighted sum (W/m**2)                                   
     real(r8), pointer :: sabg_lyr_patch         (:,:) ! patch absorbed radiation in each snow layer and top soil layer (pft,lyr) [W/m2]
     real(r8), pointer :: sabg_pen_patch         (:)   ! patch (rural) shortwave radiation penetrating top soisno layer [W/m2]

     real(r8), pointer :: sub_surf_abs_SW_patch  (:)   ! patch fraction of solar radiation absorbed below first snow layer
     real(r8), pointer :: sabv_patch             (:)   ! patch solar radiation absorbed by vegetation (W/m**2) 

     ! Currently needed by lake code 
     ! TODO (MV 8/20/2014) should be moved in the future
     real(r8), pointer :: fsds_nir_d_patch       (:)   ! patch incident direct beam nir solar radiation (W/m**2)
     real(r8), pointer :: fsds_nir_i_patch       (:)   ! patch incident diffuse nir solar radiation (W/m**2)    
     real(r8), pointer :: fsds_nir_d_ln_patch    (:)   ! patch incident direct beam nir solar radiation at local noon (W/m**2)
     real(r8), pointer :: fsr_nir_d_patch        (:)   ! patch reflected direct beam nir solar radiation (W/m**2) 
     real(r8), pointer :: fsr_nir_i_patch        (:)   ! patch reflected diffuse nir solar radiation (W/m**2)     
     real(r8), pointer :: fsr_nir_d_ln_patch     (:)   ! patch reflected direct beam nir solar radiation at local noon (W/m**2)

   contains

     procedure, public  :: Init         
     procedure, private :: InitAllocate 
     procedure, private :: InitHistory  

  end type solarabs_type
  !-----------------------------------------------------------------------

contains

  !------------------------------------------------------------------------
  subroutine Init(this, bounds)

    class(solarabs_type) :: this
    type(bounds_type), intent(in) :: bounds  

    call this%InitAllocate(bounds)
    call this%InitHistory(bounds)

  end subroutine Init

  !-----------------------------------------------------------------------
  subroutine InitAllocate(this, bounds)
    !
    ! Allocate module variables and data structures
    !
    ! !USES:
    use shr_infnan_mod, only : nan => shr_infnan_nan, assignment(=)
    use clm_varpar    , only : nlevcan, nlevcan, numrad, nlevsno
    !
    ! !ARGUMENTS:
    class(solarabs_type) :: this
    type(bounds_type), intent(in) :: bounds  
    !
    ! !LOCAL VARIABLES:
    integer :: begp, endp
    integer :: begc, endc
    integer :: begl, endl
    !---------------------------------------------------------------------

    begp = bounds%begp; endp = bounds%endp
    begc = bounds%begc; endc = bounds%endc
    begl = bounds%begl; endl = bounds%endl

    allocate(this%fsa_patch              (begp:endp))              ; this%fsa_patch              (:)   = nan
    allocate(this%parsun_z_patch         (begp:endp,1:nlevcan))    ; this%parsun_z_patch         (:,:) = nan
    allocate(this%parsha_z_patch         (begp:endp,1:nlevcan))    ; this%parsha_z_patch         (:,:) = nan 
    allocate(this%sabv_patch             (begp:endp))              ; this%sabv_patch             (:)   = nan
    allocate(this%sabg_patch             (begp:endp))              ; this%sabg_patch             (:)   = nan
    allocate(this%sabg_lyr_patch         (begp:endp,-nlevsno+1:1)) ; this%sabg_lyr_patch         (:,:) = nan
    allocate(this%sabg_pen_patch         (begp:endp))              ; this%sabg_pen_patch         (:)   = nan
    allocate(this%sabg_soil_patch        (begp:endp))              ; this%sabg_soil_patch        (:)   = nan
    allocate(this%sabg_snow_patch        (begp:endp))              ; this%sabg_snow_patch        (:)   = nan
    allocate(this%sabg_chk_patch         (begp:endp))              ; this%sabg_chk_patch         (:)   = nan
    allocate(this%sub_surf_abs_SW_patch  (begp:endp))              ; this%sub_surf_abs_SW_patch  (:)   = nan
    allocate(this%fsr_patch              (begp:endp))              ; this%fsr_patch              (:)   = nan
    allocate(this%fsr_nir_d_patch        (begp:endp))              ; this%fsr_nir_d_patch        (:)   = nan
    allocate(this%fsr_nir_i_patch        (begp:endp))              ; this%fsr_nir_i_patch        (:)   = nan
    allocate(this%fsr_nir_d_ln_patch     (begp:endp))              ; this%fsr_nir_d_ln_patch     (:)   = nan
    allocate(this%fsds_nir_d_patch       (begp:endp))              ; this%fsds_nir_d_patch       (:)   = nan
    allocate(this%fsds_nir_i_patch       (begp:endp))              ; this%fsds_nir_i_patch       (:)   = nan
    allocate(this%fsds_nir_d_ln_patch    (begp:endp))              ; this%fsds_nir_d_ln_patch    (:)   = nan

  end subroutine InitAllocate

  !-----------------------------------------------------------------------
  subroutine InitHistory(this, bounds)
    !
    ! History fields initialization
    !
    ! !USES:
    use shr_infnan_mod, only : nan => shr_infnan_nan, assignment(=)
    use clm_varpar    , only : nlevsno
    use histFileMod   , only : hist_addfld1d, hist_addfld2d
    use histFileMod   , only : no_snow_normal
    !
    ! !ARGUMENTS:
    class(solarabs_type) :: this
    type(bounds_type), intent(in) :: bounds  

    !
    ! !LOCAL VARIABLES:
    integer :: begp, endp
    integer :: begc, endc
    real(r8), pointer :: data2dptr(:,:) ! temp. pointers for slicing larger arrays
    real(r8), pointer :: ptr_1d(:)      ! pointer to 1d patch array
    !---------------------------------------------------------------------

    begp = bounds%begp; endp = bounds%endp
    begc = bounds%begc; endc = bounds%endc

    this%fsa_patch(begp:endp) = spval
    call hist_addfld1d (fname='FSA', units='W/m^2',  &
         avgflag='A', long_name='absorbed solar radiation', &
         ptr_patch=this%fsa_patch, c2l_scale_type='urbanf', default='inactive')

    this%fsr_patch(begp:endp) = spval
    call hist_addfld1d (fname='FSR', units='W/m^2',  &
         avgflag='A', long_name='reflected solar radiation', &
         ptr_patch=this%fsr_patch, c2l_scale_type='urbanf', default='inactive')
    ! Rename of FSR for Urban intercomparision project
    call hist_addfld1d (fname='SWup', units='W/m^2',  &
         avgflag='A', long_name='upwelling shortwave radiation', &
         ptr_patch=this%fsr_patch, c2l_scale_type='urbanf', default='inactive')

    this%sabg_patch(begp:endp) = spval
    call hist_addfld1d (fname='SABG', units='W/m^2',  &
         avgflag='A', long_name='solar rad absorbed by ground', &
         ptr_patch=this%sabg_patch, c2l_scale_type='urbanf', default='inactive')

     ! Currently needed by lake code - TODO should not be here
    this%fsds_nir_d_patch(begp:endp) = spval
    call hist_addfld1d (fname='FSDSND', units='W/m^2',  &
         avgflag='A', long_name='direct nir incident solar radiation', &
         ptr_patch=this%fsds_nir_d_patch, default='inactive')

    this%fsds_nir_i_patch(begp:endp) = spval
    call hist_addfld1d (fname='FSDSNI', units='W/m^2',  &
         avgflag='A', long_name='diffuse nir incident solar radiation', &
         ptr_patch=this%fsds_nir_i_patch, default='inactive')

    this%fsds_nir_d_ln_patch(begp:endp) = spval
    call hist_addfld1d (fname='FSDSNDLN', units='W/m^2',  &
         avgflag='A', long_name='direct nir incident solar radiation at local noon', &
         ptr_patch=this%fsds_nir_d_ln_patch, default='inactive')

    this%fsr_nir_d_patch(begp:endp) = spval
    call hist_addfld1d (fname='FSRND', units='W/m^2',  &
         avgflag='A', long_name='direct nir reflected solar radiation', &
         ptr_patch=this%fsr_nir_d_patch, c2l_scale_type='urbanf', default='inactive')

    this%fsr_nir_i_patch(begp:endp) = spval
    call hist_addfld1d (fname='FSRNI', units='W/m^2',  &
         avgflag='A', long_name='diffuse nir reflected solar radiation', &
         ptr_patch=this%fsr_nir_i_patch, c2l_scale_type='urbanf', default='inactive')

    this%fsr_nir_d_ln_patch(begp:endp) = spval
    call hist_addfld1d (fname='FSRNDLN', units='W/m^2',  &
         avgflag='A', long_name='direct nir reflected solar radiation at local noon', &
         ptr_patch=this%fsr_nir_d_ln_patch, c2l_scale_type='urbanf', default='inactive')

  end subroutine InitHistory

end module SolarAbsorbedType
