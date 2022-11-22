module CanopyStateType

  !------------------------------------------------------------------------------
  ! !USES:
  use shr_kind_mod    , only : r8 => shr_kind_r8
  use shr_infnan_mod  , only : nan => shr_infnan_nan, shr_infnan_isnan, assignment(=)
  use shr_log_mod     , only : errMsg => shr_log_errMsg
  use abortutils      , only : endrun
  use decompMod       , only : bounds_type
  use landunit_varcon , only : istsoil, istcrop
  use clm_varpar      , only : nlevcan
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
  type, public :: CanopyState_type

     integer  , pointer :: frac_veg_nosno_patch     (:)   ! patch fraction of vegetation not covered by snow (0 OR 1) [-] 
     integer  , pointer :: frac_veg_nosno_alb_patch (:)   ! patch fraction of vegetation not covered by snow (0 OR 1) [-] 

     real(r8) , pointer :: tlai_patch               (:)   ! patch canopy one-sided leaf area index, no burying by snow
     real(r8) , pointer :: tsai_patch               (:)   ! patch canopy one-sided stem area index, no burying by snow
     real(r8) , pointer :: elai_patch               (:)   ! patch canopy one-sided leaf area index with burying by snow
     real(r8) , pointer :: esai_patch               (:)   ! patch canopy one-sided stem area index with burying by snow
     real(r8) , pointer :: laisun_patch             (:)   ! patch patch sunlit projected leaf area index  
     real(r8) , pointer :: laisha_patch             (:)   ! patch patch shaded projected leaf area index  
     real(r8) , pointer :: laisun_z_patch           (:,:) ! patch patch sunlit leaf area for canopy layer 
     real(r8) , pointer :: laisha_z_patch           (:,:) ! patch patch shaded leaf area for canopy layer 
     real(r8) , pointer :: mlaidiff_patch           (:)   ! patch difference between lai month one and month two (for dry deposition of chemical tracers)
     real(r8) , pointer :: annlai_patch             (:,:) ! patch 12 months of monthly lai from input data set (for dry deposition of chemical tracers) 
     real(r8) , pointer :: htop_patch               (:)   ! patch canopy top (m)
     real(r8) , pointer :: hbot_patch               (:)   ! patch canopy bottom (m)
     real(r8) , pointer :: displa_patch             (:)   ! patch displacement height (m)
     real(r8) , pointer :: fsun_patch               (:)   ! patch sunlit fraction of canopy         

     real(r8) , pointer :: alt_col                  (:)   ! col current depth of thaw 
     integer  , pointer :: alt_indx_col             (:)   ! col current depth of thaw 
     real(r8) , pointer :: altmax_col               (:)   ! col maximum annual depth of thaw 
     real(r8) , pointer :: altmax_lastyear_col      (:)   ! col prior year maximum annual depth of thaw 
     integer  , pointer :: altmax_indx_col          (:)   ! col maximum annual depth of thaw 
     integer  , pointer :: altmax_lastyear_indx_col (:)   ! col prior year maximum annual depth of thaw 

     real(r8) , pointer :: dewmx_patch              (:)   ! patch maximum allowed dew [mm] 
     real(r8) , pointer :: dleaf_patch              (:)   ! patch characteristic leaf width (diameter) [m]
                                                          ! same as pftcon%dleaf()

   contains

     procedure, public  :: Init         
     procedure, private :: InitAllocate 
     procedure, private :: InitCold     

  end type CanopyState_type

  character(len=*), parameter, private :: sourcefile = &
       __FILE__
  !------------------------------------------------------------------------

contains   

  !------------------------------------------------------------------------
  subroutine Init(this, bounds)

    class(canopystate_type) :: this
    type(bounds_type), intent(in) :: bounds  

    call this%InitAllocate(bounds)
    call this%InitCold(bounds)

  end subroutine Init

  !------------------------------------------------------------------------
  subroutine InitAllocate(this, bounds)
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    class(canopystate_type) :: this
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

    allocate(this%frac_veg_nosno_patch     (begp:endp))           ; this%frac_veg_nosno_patch     (:)   = huge(1)
    allocate(this%frac_veg_nosno_alb_patch (begp:endp))           ; this%frac_veg_nosno_alb_patch (:)   = 0
    allocate(this%tlai_patch               (begp:endp))           ; this%tlai_patch               (:)   = nan
    allocate(this%tsai_patch               (begp:endp))           ; this%tsai_patch               (:)   = nan
    allocate(this%elai_patch               (begp:endp))           ; this%elai_patch               (:)   = nan
    allocate(this%esai_patch               (begp:endp))           ; this%esai_patch               (:)   = nan
    allocate(this%laisun_patch             (begp:endp))           ; this%laisun_patch             (:)   = nan
    allocate(this%laisha_patch             (begp:endp))           ; this%laisha_patch             (:)   = nan
    allocate(this%laisun_z_patch           (begp:endp,1:nlevcan)) ; this%laisun_z_patch           (:,:) = nan
    allocate(this%laisha_z_patch           (begp:endp,1:nlevcan)) ; this%laisha_z_patch           (:,:) = nan
    allocate(this%mlaidiff_patch           (begp:endp))           ; this%mlaidiff_patch           (:)   = nan
    allocate(this%annlai_patch          (12,begp:endp))           ; this%annlai_patch             (:,:) = nan
    allocate(this%htop_patch               (begp:endp))           ; this%htop_patch               (:)   = nan
    allocate(this%hbot_patch               (begp:endp))           ; this%hbot_patch               (:)   = nan
    allocate(this%displa_patch             (begp:endp))           ; this%displa_patch             (:)   = nan
    allocate(this%fsun_patch               (begp:endp))           ; this%fsun_patch               (:)   = nan

    allocate(this%alt_col                  (begc:endc))           ; this%alt_col                  (:)   = spval     
    allocate(this%altmax_col               (begc:endc))           ; this%altmax_col               (:)   = spval
    allocate(this%altmax_lastyear_col      (begc:endc))           ; this%altmax_lastyear_col      (:)   = spval
    allocate(this%alt_indx_col             (begc:endc))           ; this%alt_indx_col             (:)   = huge(1)
    allocate(this%altmax_indx_col          (begc:endc))           ; this%altmax_indx_col          (:)   = huge(1)
    allocate(this%altmax_lastyear_indx_col (begc:endc))           ; this%altmax_lastyear_indx_col (:)   = huge(1)

    allocate(this%dewmx_patch              (begp:endp))           ; this%dewmx_patch              (:)   = nan
    allocate(this%dleaf_patch              (begp:endp))           ; this%dleaf_patch              (:)   = nan

  end subroutine InitAllocate

  !-----------------------------------------------------------------------
  subroutine InitCold(this, bounds)
    !
    ! !ARGUMENTS:
    class(canopystate_type) :: this
    type(bounds_type), intent(in) :: bounds  
    !
    ! !LOCAL VARIABLES:
    integer  :: p,l,c,g 
    !-----------------------------------------------------------------------

    do p = bounds%begp, bounds%endp
       l = patch%landunit(p)

       this%frac_veg_nosno_patch(p) = 0._r8
       this%tlai_patch(p)       = 0._r8
       this%tsai_patch(p)       = 0._r8
       this%elai_patch(p)       = 0._r8
       this%esai_patch(p)       = 0._r8
       this%htop_patch(p)       = 0._r8
       this%hbot_patch(p)       = 0._r8
       this%dewmx_patch(p)      = 0.1_r8

       if (lun%itype(l) == istsoil .or. lun%itype(l) == istcrop) then
          this%laisun_patch(p) = 0._r8
          this%laisha_patch(p) = 0._r8
       end if

       ! needs to be initialized to spval to avoid problems when averaging for the accum
       ! field
       this%fsun_patch(p) = spval
    end do

    do c = bounds%begc, bounds%endc
       l = col%landunit(c)

       if (lun%itype(l) == istsoil .or. lun%itype(l) == istcrop) then
          this%alt_col(c)               = 0._r8 !iniitialized to spval for all columns
          this%altmax_col(c)            = 0._r8 !iniitialized to spval for all columns
          this%altmax_lastyear_col(c)   = 0._r8 !iniitialized to spval for all columns
          this%alt_indx_col(c)          = 0     !initiialized to huge  for all columns
          this%altmax_indx_col(c)       = 0     !initiialized to huge  for all columns
          this%altmax_lastyear_indx_col = 0     !initiialized to huge  for all columns
       end if
    end do

  end subroutine InitCold

end module CanopyStateType
