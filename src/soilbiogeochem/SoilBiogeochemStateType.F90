module SoilBiogeochemStateType

  use shr_kind_mod   , only : r8 => shr_kind_r8
  use shr_log_mod    , only : errMsg => shr_log_errMsg
  use shr_infnan_mod , only : nan => shr_infnan_nan, assignment(=)
  use decompMod      , only : bounds_type
  use abortutils     , only : endrun
  use spmdMod        , only : masterproc
  use clm_varpar     , only : ndecomp_cascade_transitions, nlevdecomp, nlevdecomp_full
  use clm_varcon     , only : spval
  use landunit_varcon, only : istsoil, istcrop
  use clm_varctl     , only : iulog
  use LandunitType   , only : lun                
  use ColumnType     , only : col                
  ! 
  ! !PUBLIC TYPES:
  implicit none
  private
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: get_spinup_latitude_term  
  !
  ! !PUBLIC TYPES:
  type, public :: soilbiogeochem_state_type

     real(r8) , pointer :: leaf_prof_patch             (:,:)   ! (1/m) profile of leaves (vertical profiles for calculating fluxes)
     real(r8) , pointer :: froot_prof_patch            (:,:)   ! (1/m) profile of fine roots (vertical profiles for calculating fluxes)
     real(r8) , pointer :: croot_prof_patch            (:,:)   ! (1/m) profile of coarse roots (vertical profiles for calculating fluxes)
     real(r8) , pointer :: stem_prof_patch             (:,:)   ! (1/m) profile of stems (vertical profiles for calculating fluxes)
     real(r8) , pointer :: fpi_vr_col                  (:,:)   ! (no units) fraction of potential immobilization 
     real(r8) , pointer :: fpi_col                     (:)     ! (no units) fraction of potential immobilization 
     real(r8),  pointer :: fpg_col                     (:)     ! (no units) fraction of potential gpp 
     real(r8) , pointer :: rf_decomp_cascade_col       (:,:,:) ! (frac) respired fraction in decomposition step 
     real(r8) , pointer :: pathfrac_decomp_cascade_col (:,:,:) ! (frac) what fraction of C leaving a given pool passes through a given transition 
     real(r8) , pointer :: nfixation_prof_col          (:,:)   ! (1/m) profile for N fixation additions 
     real(r8) , pointer :: ndep_prof_col               (:,:)   ! (1/m) profile for N fixation additions 
     real(r8) , pointer :: som_adv_coef_col            (:,:)   ! (m2/s) SOM advective flux 
     real(r8) , pointer :: som_diffus_coef_col         (:,:)   ! (m2/s) SOM diffusivity due to bio/cryo-turbation 
     real(r8) , pointer :: plant_ndemand_col           (:)     ! column-level plant N demand

   contains

     procedure, public  :: Init         
     procedure, public  :: Restart      
     procedure, private :: InitAllocate 
     procedure, private :: InitCold     

  end type soilbiogeochem_state_type
  !------------------------------------------------------------------------

contains

  !------------------------------------------------------------------------
  subroutine Init(this, bounds)

    class(soilbiogeochem_state_type) :: this
    type(bounds_type), intent(in) :: bounds  

    call this%InitAllocate ( bounds )
    call this%InitCold ( bounds ) 

  end subroutine Init

  !------------------------------------------------------------------------
  subroutine InitAllocate(this, bounds)
    !
    ! !DESCRIPTION:
    ! Initialize module data structure
    !
    ! !USES:
    use shr_infnan_mod , only : nan => shr_infnan_nan, assignment(=)
    !
    ! !ARGUMENTS:
    class(soilbiogeochem_state_type) :: this
    type(bounds_type), intent(in) :: bounds  
    !
    ! !LOCAL VARIABLES:
    integer :: begp, endp
    integer :: begc, endc
    !------------------------------------------------------------------------

    begp = bounds%begp; endp= bounds%endp
    begc = bounds%begc; endc= bounds%endc

    allocate(this%leaf_prof_patch     (begp:endp,1:nlevdecomp_full)) ; this%leaf_prof_patch     (:,:) = spval
    allocate(this%froot_prof_patch    (begp:endp,1:nlevdecomp_full)) ; this%froot_prof_patch    (:,:) = spval
    allocate(this%croot_prof_patch    (begp:endp,1:nlevdecomp_full)) ; this%croot_prof_patch    (:,:) = spval
    allocate(this%stem_prof_patch     (begp:endp,1:nlevdecomp_full)) ; this%stem_prof_patch     (:,:) = spval
    allocate(this%fpi_vr_col          (begc:endc,1:nlevdecomp_full)) ; this%fpi_vr_col          (:,:) = nan
    allocate(this%fpi_col             (begc:endc))                   ; this%fpi_col             (:)   = nan
    allocate(this%fpg_col             (begc:endc))                   ; this%fpg_col             (:)   = nan
    allocate(this%nfixation_prof_col  (begc:endc,1:nlevdecomp_full)) ; this%nfixation_prof_col  (:,:) = spval
    allocate(this%ndep_prof_col       (begc:endc,1:nlevdecomp_full)) ; this%ndep_prof_col       (:,:) = spval
    allocate(this%som_adv_coef_col    (begc:endc,1:nlevdecomp_full)) ; this%som_adv_coef_col    (:,:) = spval
    allocate(this%som_diffus_coef_col (begc:endc,1:nlevdecomp_full)) ; this%som_diffus_coef_col (:,:) = spval
    allocate(this%plant_ndemand_col   (begc:endc))                   ; this%plant_ndemand_col   (:)   = nan

    allocate(this%rf_decomp_cascade_col(begc:endc,1:nlevdecomp_full,1:ndecomp_cascade_transitions)); 
    this%rf_decomp_cascade_col(:,:,:) = nan

    allocate(this%pathfrac_decomp_cascade_col(begc:endc,1:nlevdecomp_full,1:ndecomp_cascade_transitions));     
    this%pathfrac_decomp_cascade_col(:,:,:) = nan

  end subroutine InitAllocate

  !------------------------------------------------------------------------
  subroutine initCold(this, bounds)
    !
    ! !USES:
    use spmdMod    , only : masterproc
    use fileutils  , only : getfil
    use ncdio_pio
    !
    ! !ARGUMENTS:
    class(soilbiogeochem_state_type) :: this
    type(bounds_type), intent(in) :: bounds   
    !
    ! !LOCAL VARIABLES:
    integer               :: g,l,c,p,n,j,m            ! indices
    integer               :: dimid                    ! dimension id
    integer               :: ier                      ! error status
    logical               :: readvar 
    integer               :: begc, endc
    !-----------------------------------------------------------------------

    begc = bounds%begc; endc= bounds%endc

    ! --------------------------------------------------------------------
    ! Initialize terms needed for dust model
    ! --------------------------------------------------------------------
       
    do c = bounds%begc, bounds%endc
       l = col%landunit(c)
       if (lun%ifspecial(l)) then
          this%fpi_col (c) = spval
          this%fpg_col (c) = spval
          do j = 1,nlevdecomp_full
             this%fpi_vr_col(c,j) = spval
          end do
       end if

       if (lun%itype(l) == istsoil .or. lun%itype(l) == istcrop) then
          ! initialize fpi_vr so that levels below nlevsoi are not nans
          this%fpi_vr_col(c,1:nlevdecomp_full)          = 0._r8 
          this%som_adv_coef_col(c,1:nlevdecomp_full)    = 0._r8 
          this%som_diffus_coef_col(c,1:nlevdecomp_full) = 0._r8 

          ! initialize the profiles for converting to vertically resolved carbon pools
          this%nfixation_prof_col(c,1:nlevdecomp_full)  = 0._r8 
          this%ndep_prof_col(c,1:nlevdecomp_full)       = 0._r8 
       end if
    end do

  end subroutine initCold

  !------------------------------------------------------------------------
  subroutine Restart(this, bounds, ncid, flag)
    !
    ! !USES:
    use shr_log_mod, only : errMsg => shr_log_errMsg
    use spmdMod    , only : masterproc
    use abortutils , only : endrun
    use restUtilMod
    use ncdio_pio
    !
    ! !ARGUMENTS:
    class(soilbiogeochem_state_type) :: this
    type(bounds_type), intent(in)    :: bounds 
    type(file_desc_t), intent(inout) :: ncid   
    character(len=*) , intent(in)    :: flag   
    !
    ! !LOCAL VARIABLES:
    integer, pointer  :: temp1d(:)  ! temporary
    integer           :: p,j,c,i    ! indices
    logical           :: readvar    ! determine if variable is on initial file
    real(r8), pointer :: ptr2d(:,:) ! temp. pointers for slicing larger arrays
    real(r8), pointer :: ptr1d(:)   ! temp. pointers for slicing larger arrays
    !-----------------------------------------------------------------------
  
    ptr1d => this%fpi_vr_col(:,1) ! nlevdecomp = 1; so treat as 1D variable
    call restartvar(ncid=ncid, flag=flag, varname='fpi', xtype=ncd_double,  &
         dim1name='column', &
         long_name='fraction of potential immobilization',  units='unitless', &
         interpinic_flag='interp' , readvar=readvar, data=ptr1d)

    call restartvar(ncid=ncid, flag=flag, varname='fpg', xtype=ncd_double,  &
         dim1name='column', &
         long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%fpg_col) 

  end subroutine Restart


  function get_spinup_latitude_term(latitude) result(ans)

    !!DESCRIPTION:
    ! calculate a logistic function to scale spinup factors so that spinup is more accelerated in high latitude regions
    !
    ! !REVISION HISTORY
    ! charlie koven, nov. 2015
    !
    ! !ARGUMENTS:
    real(r8), intent(in) :: latitude
    !
    ! !LOCAL VARIABLES:
    real(r8) :: ans

    ans = 1._r8 + 50._r8 / ( 1._r8 + exp(-0.15_r8 * (abs(latitude) - 60._r8) ) )
    
    return
  end function get_spinup_latitude_term

end module SoilBiogeochemStateType
