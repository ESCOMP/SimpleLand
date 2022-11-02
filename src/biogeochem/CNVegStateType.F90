module CNVegStateType

  use shr_kind_mod   , only : r8 => shr_kind_r8
  use shr_log_mod    , only : errMsg => shr_log_errMsg
  use shr_infnan_mod , only : nan => shr_infnan_nan, assignment(=)
  use decompMod      , only : bounds_type
  use abortutils     , only : endrun
  use spmdMod        , only : masterproc
  use clm_varpar     , only : nlevsno, nlevgrnd, nlevlak, nlevsoi
  use clm_varctl     , only : iulog, fsurdat
  use clm_varcon     , only : spval, ispval, grlnd
  use landunit_varcon, only : istsoil, istcrop
  use LandunitType   , only : lun                
  use ColumnType     , only : col                
  use PatchType      , only : patch                
  !
  ! !PUBLIC TYPES:
  implicit none
  private
  !
  ! !PUBLIC TYPES:
  type, public :: cnveg_state_type

     integer  , pointer :: burndate_patch              (:)     ! patch crop burn date
     real(r8) , pointer :: dwt_smoothed_patch          (:)     ! change in patch weight (-1 to 1) on the gridcell in this time step; changes in first time step of year are smoothed (dribbled) over the whole year

     ! Prognostic crop model
     !
     ! TODO(wjs, 2016-02-22) Most / all of these crop-specific state variables should be
     ! moved to CropType
     real(r8) , pointer :: hdidx_patch                 (:)     ! patch cold hardening index?
     real(r8) , pointer :: cumvd_patch                 (:)     ! patch cumulative vernalization d?ependence?
     real(r8) , pointer :: gddmaturity_patch           (:)     ! patch growing degree days (gdd) needed to harvest (ddays)
     real(r8) , pointer :: huileaf_patch               (:)     ! patch heat unit index needed from planting to leaf emergence
     real(r8) , pointer :: huigrain_patch              (:)     ! patch heat unit index needed to reach vegetative maturity
     real(r8) , pointer :: aleafi_patch                (:)     ! patch saved leaf allocation coefficient from phase 2
     real(r8) , pointer :: astemi_patch                (:)     ! patch saved stem allocation coefficient from phase 2
     real(r8) , pointer :: aleaf_patch                 (:)     ! patch leaf allocation coefficient
     real(r8) , pointer :: astem_patch                 (:)     ! patch stem allocation coefficient
     real(r8) , pointer :: htmx_patch                  (:)     ! patch max hgt attained by a crop during yr (m)
     integer  , pointer :: peaklai_patch               (:)     ! patch 1: max allowed lai; 0: not at max

     integer  , pointer :: idop_patch                  (:)     ! patch date of planting

     real(r8) , pointer :: gdp_lf_col                  (:)     ! col global real gdp data (k US$/capita)
     real(r8) , pointer :: peatf_lf_col                (:)     ! col global peatland fraction data (0-1)
     integer  , pointer :: abm_lf_col                  (:)     ! col global peak month of crop fire emissions 

     real(r8) , pointer :: lgdp_col                    (:)     ! col gdp limitation factor for fire occurrence (0-1)
     real(r8) , pointer :: lgdp1_col                   (:)     ! col gdp limitation factor for fire spreading (0-1)
     real(r8) , pointer :: lpop_col                    (:)     ! col pop limitation factor for fire spreading (0-1)

     real(r8) , pointer :: tempavg_t2m_patch           (:)     ! patch temporary average 2m air temperature (K)
     real(r8) , pointer :: annavg_t2m_patch            (:)     ! patch annual average 2m air temperature (K)
     real(r8) , pointer :: annavg_t2m_col              (:)     ! col annual average of 2m air temperature, averaged from patch-level (K)
     real(r8) , pointer :: annsum_counter_col          (:)     ! col seconds since last annual accumulator turnover

     ! Fire
     real(r8) , pointer :: nfire_col                   (:)     ! col fire counts (count/km2/sec), valid only in Reg. C
     real(r8) , pointer :: fsr_col                     (:)     ! col fire spread rate at column level (m/s)
     real(r8) , pointer :: fd_col                      (:)     ! col fire duration at column level (hr)
     real(r8) , pointer :: lfc_col                     (:)     ! col conversion area fraction of BET and BDT that haven't burned before (/timestep)
     real(r8) , pointer :: lfc2_col                    (:)     ! col conversion area fraction of BET and BDT that burned (/sec)
     real(r8) , pointer :: dtrotr_col                  (:)     ! col annual decreased fraction coverage of BET on the gridcell (0-1)
     real(r8) , pointer :: trotr1_col                  (:)     ! col patch weight of BET on the column (0-1)
     real(r8) , pointer :: trotr2_col                  (:)     ! col patch weight of BDT on the column (0-1)
     real(r8) , pointer :: cropf_col                   (:)     ! col crop fraction in veg column (0-1)
     real(r8) , pointer :: baf_crop_col                (:)     ! col baf for cropland(/sec)
     real(r8) , pointer :: baf_peatf_col               (:)     ! col baf for peatland (/sec)
     real(r8) , pointer :: fbac_col                    (:)     ! col total burned area out of conversion (/sec)
     real(r8) , pointer :: fbac1_col                   (:)     ! col burned area out of conversion region due to land use fire (/sec)
     real(r8) , pointer :: wtlf_col                    (:)     ! col fractional coverage of non-crop Patches (0-1)
     real(r8) , pointer :: lfwt_col                    (:)     ! col fractional coverage of non-crop and non-bare-soil Patches (0-1)
     real(r8) , pointer :: farea_burned_col            (:)     ! col fractional area burned (/sec) 

     real(r8), pointer :: dormant_flag_patch           (:)     ! patch dormancy flag
     real(r8), pointer :: days_active_patch            (:)     ! patch number of days since last dormancy
     real(r8), pointer :: onset_flag_patch             (:)     ! patch onset flag
     real(r8), pointer :: onset_counter_patch          (:)     ! patch onset days counter
     real(r8), pointer :: onset_gddflag_patch          (:)     ! patch onset flag for growing degree day sum
     real(r8), pointer :: onset_fdd_patch              (:)     ! patch onset freezing degree days counter
     real(r8), pointer :: onset_gdd_patch              (:)     ! patch onset growing degree days
     real(r8), pointer :: onset_swi_patch              (:)     ! patch onset soil water index
     real(r8), pointer :: offset_flag_patch            (:)     ! patch offset flag
     real(r8), pointer :: offset_counter_patch         (:)     ! patch offset days counter
     real(r8), pointer :: offset_fdd_patch             (:)     ! patch offset freezing degree days counter
     real(r8), pointer :: offset_swi_patch             (:)     ! patch offset soil water index
     real(r8), pointer :: grain_flag_patch             (:)     ! patch 1: grain fill stage; 0: not
     real(r8), pointer :: lgsf_patch                   (:)     ! patch long growing season factor [0-1]
     real(r8), pointer :: bglfr_patch                  (:)     ! patch background litterfall rate (1/s)
     real(r8), pointer :: bgtr_patch                   (:)     ! patch background transfer growth rate (1/s)
     real(r8), pointer :: c_allometry_patch            (:)     ! patch C allocation index (DIM)
     real(r8), pointer :: n_allometry_patch            (:)     ! patch N allocation index (DIM)

     real(r8), pointer :: tempsum_potential_gpp_patch  (:)     ! patch temporary annual sum of potential GPP
     real(r8), pointer :: annsum_potential_gpp_patch   (:)     ! patch annual sum of potential GPP
     real(r8), pointer :: tempmax_retransn_patch       (:)     ! patch temporary annual max of retranslocated N pool (gN/m2)
     real(r8), pointer :: annmax_retransn_patch        (:)     ! patch annual max of retranslocated N pool (gN/m2)
     real(r8), pointer :: downreg_patch                (:)     ! patch fractional reduction in GPP due to N limitation (DIM)
     real(r8), pointer :: leafcn_offset_patch          (:)     ! patch leaf C:N used by FUN
     real(r8), pointer :: plantCN_patch                (:)     ! patch plant C:N used by FUN

   contains

     procedure, public  :: Init         
     procedure, public  :: Restart      
     procedure, private :: InitAllocate 
     procedure, private :: InitCold     

  end type cnveg_state_type
  !------------------------------------------------------------------------

  character(len=*), parameter, private :: sourcefile = &
       __FILE__

contains

  !------------------------------------------------------------------------
  subroutine Init(this, bounds)

    class(cnveg_state_type) :: this
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
    class(cnveg_state_type) :: this
    type(bounds_type), intent(in) :: bounds  
    !
    ! !LOCAL VARIABLES:
    integer :: begp, endp
    integer :: begc, endc
    logical :: allows_non_annual_delta
    !------------------------------------------------------------------------

    begp = bounds%begp; endp= bounds%endp
    begc = bounds%begc; endc= bounds%endc

    allocate(this%burndate_patch      (begp:endp))                   ; this%burndate_patch      (:)   = ispval
    allocate(this%dwt_smoothed_patch  (begp:endp))                   ; this%dwt_smoothed_patch  (:)   = nan

    allocate(this%hdidx_patch         (begp:endp))                   ; this%hdidx_patch         (:)   = nan
    allocate(this%cumvd_patch         (begp:endp))                   ; this%cumvd_patch         (:)   = nan
    allocate(this%gddmaturity_patch   (begp:endp))                   ; this%gddmaturity_patch   (:)   = spval
    allocate(this%huileaf_patch       (begp:endp))                   ; this%huileaf_patch       (:)   = nan
    allocate(this%huigrain_patch      (begp:endp))                   ; this%huigrain_patch      (:)   = 0.0_r8
    allocate(this%aleafi_patch        (begp:endp))                   ; this%aleafi_patch        (:)   = nan
    allocate(this%astemi_patch        (begp:endp))                   ; this%astemi_patch        (:)   = nan
    allocate(this%aleaf_patch         (begp:endp))                   ; this%aleaf_patch         (:)   = nan
    allocate(this%astem_patch         (begp:endp))                   ; this%astem_patch         (:)   = nan
    allocate(this%htmx_patch          (begp:endp))                   ; this%htmx_patch          (:)   = 0.0_r8
    allocate(this%peaklai_patch       (begp:endp))                   ; this%peaklai_patch       (:)   = 0

    allocate(this%idop_patch          (begp:endp))                   ; this%idop_patch          (:)   = huge(1)

    allocate(this%gdp_lf_col          (begc:endc))                   ;
    allocate(this%peatf_lf_col        (begc:endc))                   ; 
    allocate(this%abm_lf_col          (begc:endc))                   ; 

    allocate(this%lgdp_col            (begc:endc))                   ; 
    allocate(this%lgdp1_col           (begc:endc))                   ; 
    allocate(this%lpop_col            (begc:endc))                   ;  

    allocate(this%tempavg_t2m_patch   (begp:endp))                   ; this%tempavg_t2m_patch   (:)   = nan
    allocate(this%annsum_counter_col  (begc:endc))                   ; this%annsum_counter_col  (:)   = nan
    allocate(this%annavg_t2m_col      (begc:endc))                   ; this%annavg_t2m_col      (:)   = nan
    allocate(this%annavg_t2m_patch    (begp:endp))                   ; this%annavg_t2m_patch    (:)   = nan

    allocate(this%nfire_col           (begc:endc))                   ; this%nfire_col           (:)   = spval
    allocate(this%fsr_col             (begc:endc))                   ; this%fsr_col             (:)   = nan
    allocate(this%fd_col              (begc:endc))                   ; this%fd_col              (:)   = nan
    allocate(this%lfc_col             (begc:endc))                   ; this%lfc_col             (:)   = spval
    allocate(this%lfc2_col            (begc:endc))                   ; this%lfc2_col            (:)   = 0._r8
    allocate(this%dtrotr_col          (begc:endc))                   ; this%dtrotr_col          (:)   = 0._r8
    allocate(this%trotr1_col          (begc:endc))                   ; this%trotr1_col          (:)   = 0._r8
    allocate(this%trotr2_col          (begc:endc))                   ; this%trotr2_col          (:)   = 0._r8
    allocate(this%cropf_col           (begc:endc))                   ; this%cropf_col           (:)   = nan
    allocate(this%baf_crop_col        (begc:endc))                   ; this%baf_crop_col        (:)   = nan
    allocate(this%baf_peatf_col       (begc:endc))                   ; this%baf_peatf_col       (:)   = nan
    allocate(this%fbac_col            (begc:endc))                   ; this%fbac_col            (:)   = nan
    allocate(this%fbac1_col           (begc:endc))                   ; this%fbac1_col           (:)   = nan
    allocate(this%wtlf_col            (begc:endc))                   ; this%wtlf_col            (:)   = nan
    allocate(this%lfwt_col            (begc:endc))                   ; this%lfwt_col            (:)   = nan
    allocate(this%farea_burned_col    (begc:endc))                   ; this%farea_burned_col    (:)   = nan

    allocate(this%dormant_flag_patch          (begp:endp)) ;    this%dormant_flag_patch          (:) = nan
    allocate(this%days_active_patch           (begp:endp)) ;    this%days_active_patch           (:) = nan
    allocate(this%onset_flag_patch            (begp:endp)) ;    this%onset_flag_patch            (:) = nan
    allocate(this%onset_counter_patch         (begp:endp)) ;    this%onset_counter_patch         (:) = nan
    allocate(this%onset_gddflag_patch         (begp:endp)) ;    this%onset_gddflag_patch         (:) = nan
    allocate(this%onset_fdd_patch             (begp:endp)) ;    this%onset_fdd_patch             (:) = nan
    allocate(this%onset_gdd_patch             (begp:endp)) ;    this%onset_gdd_patch             (:) = nan
    allocate(this%onset_swi_patch             (begp:endp)) ;    this%onset_swi_patch             (:) = nan
    allocate(this%offset_flag_patch           (begp:endp)) ;    this%offset_flag_patch           (:) = nan
    allocate(this%offset_counter_patch        (begp:endp)) ;    this%offset_counter_patch        (:) = nan
    allocate(this%offset_fdd_patch            (begp:endp)) ;    this%offset_fdd_patch            (:) = nan
    allocate(this%offset_swi_patch            (begp:endp)) ;    this%offset_swi_patch            (:) = nan
    allocate(this%grain_flag_patch            (begp:endp)) ;    this%grain_flag_patch            (:) = nan
    allocate(this%lgsf_patch                  (begp:endp)) ;    this%lgsf_patch                  (:) = nan
    allocate(this%bglfr_patch                 (begp:endp)) ;    this%bglfr_patch                 (:) = nan
    allocate(this%bgtr_patch                  (begp:endp)) ;    this%bgtr_patch                  (:) = nan
    allocate(this%c_allometry_patch           (begp:endp)) ;    this%c_allometry_patch           (:) = nan
    allocate(this%n_allometry_patch           (begp:endp)) ;    this%n_allometry_patch           (:) = nan
    allocate(this%tempsum_potential_gpp_patch (begp:endp)) ;    this%tempsum_potential_gpp_patch (:) = nan
    allocate(this%annsum_potential_gpp_patch  (begp:endp)) ;    this%annsum_potential_gpp_patch  (:) = nan
    allocate(this%tempmax_retransn_patch      (begp:endp)) ;    this%tempmax_retransn_patch      (:) = nan
    allocate(this%annmax_retransn_patch       (begp:endp)) ;    this%annmax_retransn_patch       (:) = nan
    allocate(this%downreg_patch               (begp:endp)) ;    this%downreg_patch               (:) = nan
    allocate(this%leafcn_offset_patch         (begp:endp)) ;    this%leafcn_offset_patch         (:) = nan
    allocate(this%plantCN_patch               (begp:endp)) ;    this%plantCN_patch               (:) = nan

  end subroutine InitAllocate

  !------------------------------------------------------------------------
  subroutine initCold(this, bounds)
    !
    ! !USES:
    use spmdMod    , only : masterproc
    use fileutils  , only : getfil
    use clm_varctl , only : nsrest, nsrStartup
    use ncdio_pio
    !
    ! !ARGUMENTS:
    class(cnveg_state_type) :: this
    type(bounds_type), intent(in) :: bounds   
    !
    ! !LOCAL VARIABLES:
    integer               :: g,l,c,p,n,j,m            ! indices
    real(r8) ,pointer     :: gdp (:)                  ! global gdp data (needs to be a pointer for use in ncdio)
    real(r8) ,pointer     :: peatf (:)                ! global peatf data (needs to be a pointer for use in ncdio)
    integer  ,pointer     :: abm (:)                  ! global abm data (needs to be a pointer for use in ncdio)
    real(r8) ,pointer     :: gti (:)                  ! read in - fmax (needs to be a pointer for use in ncdio)
    integer               :: dimid                    ! dimension id
    integer               :: ier                      ! error status
    type(file_desc_t)     :: ncid                     ! netcdf id
    logical               :: readvar 
    character(len=256)    :: locfn                    ! local filename
    integer               :: begc, endc
    integer               :: begg, endg
    !-----------------------------------------------------------------------

    begc = bounds%begc; endc= bounds%endc
    begg = bounds%begg; endg= bounds%endg

    ! --------------------------------------------------------------------
    ! Open surface dataset
    ! --------------------------------------------------------------------

    call getfil (fsurdat, locfn, 0)
    call ncd_pio_openfile (ncid, locfn, 0)

    ! --------------------------------------------------------------------
    ! Read in GDP data 
    ! --------------------------------------------------------------------

    allocate(gdp(bounds%begg:bounds%endg))
    call ncd_io(ncid=ncid, varname='gdp', flag='read', data=gdp, dim1name=grlnd, readvar=readvar)
    if (.not. readvar) then
       call endrun(msg=' ERROR: gdp NOT on surfdata file'//errMsg(sourcefile, __LINE__)) 
    end if
    do c = bounds%begc, bounds%endc
       g = col%gridcell(c)
       this%gdp_lf_col(c) = gdp(g)
    end do
    deallocate(gdp)

    ! --------------------------------------------------------------------
    ! Read in peatf data 
    ! --------------------------------------------------------------------

    allocate(peatf(bounds%begg:bounds%endg))
    call ncd_io(ncid=ncid, varname='peatf', flag='read', data=peatf, dim1name=grlnd, readvar=readvar)
    if (.not. readvar) then
       call endrun(msg=' ERROR: peatf NOT on surfdata file'//errMsg(sourcefile, __LINE__)) 
    end if
    do c = bounds%begc, bounds%endc
       g = col%gridcell(c)
       this%peatf_lf_col(c) = peatf(g)
    end do
    deallocate(peatf)

    ! --------------------------------------------------------------------
    ! Read in ABM data 
    ! --------------------------------------------------------------------

    allocate(abm(bounds%begg:bounds%endg))
    call ncd_io(ncid=ncid, varname='abm', flag='read', data=abm, dim1name=grlnd, readvar=readvar)
    if (.not. readvar) then
       call endrun(msg=' ERROR: abm NOT on surfdata file'//errMsg(sourcefile, __LINE__)) 
    end if
    do c = bounds%begc, bounds%endc
       g = col%gridcell(c)
       this%abm_lf_col(c) = abm(g)
    end do
    deallocate(abm)

    ! Close file

    call ncd_pio_closefile(ncid)

    if (masterproc) then
       write(iulog,*) 'Successfully read fmax, soil color, sand and clay boundary data'
       write(iulog,*)
    endif
    
    ! --------------------------------------------------------------------
    ! Initialize terms needed for dust model
    ! TODO - move these terms to DUSTMod module variables 
    ! --------------------------------------------------------------------
       
    do c = bounds%begc, bounds%endc
       l = col%landunit(c)
       if (lun%ifspecial(l)) then
          this%annsum_counter_col (c) = spval
          this%annavg_t2m_col     (c) = spval
          this%nfire_col          (c) = spval
          this%baf_crop_col       (c) = spval
          this%baf_peatf_col      (c) = spval
          this%fbac_col           (c) = spval
          this%fbac1_col          (c) = spval
          this%farea_burned_col   (c) = spval
       end if

       if (lun%itype(l) == istsoil .or. lun%itype(l) == istcrop) then
          this%annsum_counter_col(c) = 0._r8   
          this%annavg_t2m_col(c)     = 280._r8 

          ! fire related variables 
          this%baf_crop_col(c)       = 0._r8 
          this%baf_peatf_col(c)      = 0._r8 
          this%fbac_col(c)           = 0._r8 
          this%fbac1_col(c)          = 0._r8 
          this%farea_burned_col(c)   = 0._r8 
          this%nfire_col(c)          = 0._r8
       end if
    end do

    ! ecophysiological and phenology variables

    do p = bounds%begp,bounds%endp
       l = patch%landunit(p)

       if (lun%ifspecial(l)) then
          this%annavg_t2m_patch  (p)          = spval
          this%tempavg_t2m_patch (p)          = spval
          this%dormant_flag_patch(p)          = spval
          this%days_active_patch(p)           = spval
          this%onset_flag_patch(p)            = spval
          this%onset_counter_patch(p)         = spval
          this%onset_gddflag_patch(p)         = spval
          this%onset_fdd_patch(p)             = spval
          this%onset_gdd_patch(p)             = spval
          this%onset_swi_patch(p)             = spval
          this%offset_flag_patch(p)           = spval
          this%offset_counter_patch(p)        = spval
          this%offset_fdd_patch(p)            = spval
          this%offset_swi_patch(p)            = spval
          this%grain_flag_patch(p)            = spval
          this%lgsf_patch(p)                  = spval
          this%bglfr_patch(p)                 = spval
          this%bgtr_patch(p)                  = spval
          this%c_allometry_patch(p)           = spval
          this%n_allometry_patch(p)           = spval
          this%tempsum_potential_gpp_patch(p) = spval
          this%annsum_potential_gpp_patch(p)  = spval
          this%tempmax_retransn_patch(p)      = spval
          this%annmax_retransn_patch(p)       = spval
          this%downreg_patch(p)               = spval
          this%leafcn_offset_patch(p)         = spval
          this%plantCN_patch(p)               = spval
       end if

       if (lun%itype(l) == istsoil .or. lun%itype(l) == istcrop) then
          ! phenology variables
          this%dormant_flag_patch(p)   = 1._r8
          this%days_active_patch(p)    = 0._r8
          this%onset_flag_patch(p)     = 0._r8
          this%onset_counter_patch(p)  = 0._r8
          this%onset_gddflag_patch(p)  = 0._r8
          this%onset_fdd_patch(p)      = 0._r8
          this%onset_gdd_patch(p)      = 0._r8
          this%onset_swi_patch(p)      = 0._r8
          this%offset_flag_patch(p)    = 0._r8
          this%offset_counter_patch(p) = 0._r8
          this%offset_fdd_patch(p)     = 0._r8
          this%offset_swi_patch(p)     = 0._r8
          this%lgsf_patch(p)           = 0._r8
          this%bglfr_patch(p)          = 0._r8
          this%bgtr_patch(p)           = 0._r8
          this%annavg_t2m_patch(p)     = 280._r8
          this%tempavg_t2m_patch(p)    = 0._r8
          this%grain_flag_patch(p)     = 0._r8

          ! non-phenology variables
          this%c_allometry_patch(p)           = 0._r8
          this%n_allometry_patch(p)           = 0._r8
          this%tempsum_potential_gpp_patch(p) = 0._r8
          this%annsum_potential_gpp_patch(p)  = 0._r8
          this%tempmax_retransn_patch(p)      = 0._r8
          this%annmax_retransn_patch(p)       = 0._r8
          this%downreg_patch(p)               = 0._r8
          this%leafcn_offset_patch(p)         = spval 
          this%plantCN_patch(p)               = spval 
       end if

    end do

    ! fire variables

    do c = bounds%begc,bounds%endc
       this%lfc2_col(c) = 0._r8
    end do

  end subroutine initCold

  !------------------------------------------------------------------------
  subroutine Restart(this, bounds, ncid, flag, cnveg_carbonstate, &
             cnveg_nitrogenstate, filter_reseed_patch, num_reseed_patch)
    !
    ! !USES:
    use shr_log_mod, only : errMsg => shr_log_errMsg
    use spmdMod    , only : masterproc
    use abortutils , only : endrun
    use CNVegNitrogenStateType, only: cnveg_nitrogenstate_type
    use CNVegCarbonStateType  , only: cnveg_carbonstate_type
    use restUtilMod
    use ncdio_pio
    use pftconMod , only : pftcon
    !
    ! !ARGUMENTS:
    class(cnveg_state_type) :: this
    type(bounds_type), intent(in)    :: bounds 
    type(file_desc_t), intent(inout) :: ncid   
    character(len=*) , intent(in)    :: flag   
    type(cnveg_nitrogenstate_type), intent(in) :: cnveg_nitrogenstate
    type(cnveg_carbonstate_type)  , intent(in) :: cnveg_carbonstate
    integer                       , intent(out), optional :: filter_reseed_patch(:)
    integer                       , intent(out), optional :: num_reseed_patch
    !
    ! !LOCAL VARIABLES:
    integer          :: j,c,i,p ! indices
    logical          :: readvar   ! determine if variable is on initial file
    real(r8), pointer :: ptr2d(:,:) ! temp. pointers for slicing larger arrays
    real(r8), pointer :: ptr1d(:)   ! temp. pointers for slicing larger arrays
    !-----------------------------------------------------------------------

    call restartvar(ncid=ncid, flag=flag, varname='dormant_flag', xtype=ncd_double,  &
         dim1name='pft', &
         long_name='dormancy flag', units='unitless', &
         interpinic_flag='interp', readvar=readvar, data=this%dormant_flag_patch) 

    call restartvar(ncid=ncid, flag=flag, varname='days_active', xtype=ncd_double,  &
         dim1name='pft', &
         long_name='number of days since last dormancy', units='days' , &
         interpinic_flag='interp', readvar=readvar, data=this%days_active_patch) 

    call restartvar(ncid=ncid, flag=flag, varname='onset_flag', xtype=ncd_double,  &
         dim1name='pft', &
         long_name='flag if critical growing degree-day sum is exceeded', units='unitless' , &
         interpinic_flag='interp', readvar=readvar, data=this%onset_flag_patch) 

    call restartvar(ncid=ncid, flag=flag, varname='onset_counter', xtype=ncd_double,  &
         dim1name='pft', &
         long_name='onset days counter', units='sec' , &
         interpinic_flag='interp', readvar=readvar, data=this%onset_counter_patch) 

    call restartvar(ncid=ncid, flag=flag, varname='onset_gddflag', xtype=ncd_double,  &
         dim1name='pft', &
         long_name='onset flag for growing degree day sum', units='' , &
         interpinic_flag='interp', readvar=readvar, data=this%onset_gddflag_patch) 

    call restartvar(ncid=ncid, flag=flag, varname='onset_fdd', xtype=ncd_double,  &
         dim1name='pft', &
         long_name='onset freezing degree days counter', units='days' , &
         interpinic_flag='interp', readvar=readvar, data=this%onset_fdd_patch) 

    call restartvar(ncid=ncid, flag=flag, varname='onset_gdd', xtype=ncd_double,  &
         dim1name='pft', &
         long_name='onset growing degree days', units='days' , &
         interpinic_flag='interp', readvar=readvar, data=this%onset_gdd_patch) 

    call restartvar(ncid=ncid, flag=flag, varname='onset_swi', xtype=ncd_double,  &
         dim1name='pft', &
         long_name='onset soil water index', units='days' , &
         interpinic_flag='interp', readvar=readvar, data=this%onset_swi_patch) 

    call restartvar(ncid=ncid, flag=flag, varname='offset_flag', xtype=ncd_double,  &
         dim1name='pft', &
         long_name='offset flag', units='unitless' , &
         interpinic_flag='interp', readvar=readvar, data=this%offset_flag_patch) 

    call restartvar(ncid=ncid, flag=flag, varname='offset_counter', xtype=ncd_double,  &
         dim1name='pft', &
         long_name='offset days counter', units='sec' , &
         interpinic_flag='interp', readvar=readvar, data=this%offset_counter_patch) 

    call restartvar(ncid=ncid, flag=flag, varname='offset_fdd', xtype=ncd_double,  &
         dim1name='pft', &
         long_name='offset freezing degree days counter', units='days' , &
         interpinic_flag='interp', readvar=readvar, data=this%offset_fdd_patch) 

    call restartvar(ncid=ncid, flag=flag, varname='offset_swi', xtype=ncd_double,  &
         dim1name='pft', &
         long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%offset_swi_patch) 

    call restartvar(ncid=ncid, flag=flag, varname='lgsf', xtype=ncd_double,  &
         dim1name='pft', &
         long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%lgsf_patch) 

    call restartvar(ncid=ncid, flag=flag, varname='bglfr', xtype=ncd_double,  &
         dim1name='pft', &
         long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%bglfr_patch) 

    call restartvar(ncid=ncid, flag=flag, varname='bgtr', xtype=ncd_double,  &
         dim1name='pft', &
         long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%bgtr_patch) 

    call restartvar(ncid=ncid, flag=flag, varname='annavg_t2m', xtype=ncd_double,  &
         dim1name='pft', &
         long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%annavg_t2m_patch) 

    call restartvar(ncid=ncid, flag=flag, varname='tempavg_t2m', xtype=ncd_double,  &
         dim1name='pft', &
         long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%tempavg_t2m_patch) 

    call restartvar(ncid=ncid, flag=flag, varname='c_allometry', xtype=ncd_double,  &
         dim1name='pft', &
         long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%c_allometry_patch) 

    call restartvar(ncid=ncid, flag=flag, varname='n_allometry', xtype=ncd_double,  &
         dim1name='pft', &
         long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%n_allometry_patch) 

    call restartvar(ncid=ncid, flag=flag, varname='tempsum_potential_gpp', xtype=ncd_double,  &
         dim1name='pft', &
         long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%tempsum_potential_gpp_patch) 

    call restartvar(ncid=ncid, flag=flag, varname='annsum_potential_gpp', xtype=ncd_double,  &
         dim1name='pft', &
         long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%annsum_potential_gpp_patch) 

    call restartvar(ncid=ncid, flag=flag, varname='tempmax_retransn', xtype=ncd_double,  &
         dim1name='pft', &
         long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%tempmax_retransn_patch) 

    call restartvar(ncid=ncid, flag=flag, varname='annmax_retransn', xtype=ncd_double,  &
         dim1name='pft', &
         long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%annmax_retransn_patch) 

    call restartvar(ncid=ncid, flag=flag, varname='downreg', xtype=ncd_double,  &
         dim1name='pft', &
         long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%downreg_patch) 

    call restartvar(ncid=ncid, flag=flag, varname='leafcn_offset', xtype=ncd_double,  &
         dim1name='pft', &
         long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%leafcn_offset_patch)
     
    call restartvar(ncid=ncid, flag=flag, varname='plantCN', xtype=ncd_double,  &
         dim1name='pft', &
         long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%plantCN_patch)

    call restartvar(ncid=ncid, flag=flag, varname='annsum_counter', xtype=ncd_double,  &
         dim1name='column', &
         long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%annsum_counter_col) 

    call restartvar(ncid=ncid, flag=flag, varname='burndate', xtype=ncd_int,  &
         dim1name='pft', &
         long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%burndate_patch) 

    call restartvar(ncid=ncid, flag=flag, varname='lfc', xtype=ncd_double,  &
         dim1name='column', &
         long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%lfc_col) 

    call restartvar(ncid=ncid, flag=flag, varname='cannavg_t2m', xtype=ncd_double,  &
         dim1name='column', &
         long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%annavg_t2m_col) 

    if ( flag == 'read' .and. num_reseed_patch > 0 )then
       if ( masterproc ) write(iulog, *) 'Reseed dead plants for CNVegState'
       do i = 1, num_reseed_patch
          p = filter_reseed_patch(i)
          ! phenology variables
          this%dormant_flag_patch(p)   = 1._r8
          this%days_active_patch(p)    = 0._r8
          this%onset_flag_patch(p)     = 0._r8
          this%onset_counter_patch(p)  = 0._r8
          this%onset_gddflag_patch(p)  = 0._r8
          this%onset_fdd_patch(p)      = 0._r8
          this%onset_gdd_patch(p)      = 0._r8
          this%onset_swi_patch(p)      = 0._r8
          this%offset_flag_patch(p)    = 0._r8
          this%offset_counter_patch(p) = 0._r8
          this%offset_fdd_patch(p)     = 0._r8
          this%offset_swi_patch(p)     = 0._r8
          this%lgsf_patch(p)           = 0._r8
          this%bglfr_patch(p)          = 0._r8
          this%bgtr_patch(p)           = 0._r8
          this%annavg_t2m_patch(p)     = 280._r8
          this%tempavg_t2m_patch(p)    = 0._r8
          this%grain_flag_patch(p)     = 0._r8

          this%c_allometry_patch(p)           = 0._r8
          this%n_allometry_patch(p)           = 0._r8
          this%tempsum_potential_gpp_patch(p) = 0._r8
          this%annsum_potential_gpp_patch(p)  = 0._r8
          this%tempmax_retransn_patch(p)      = 0._r8
          this%annmax_retransn_patch(p)       = 0._r8
          this%downreg_patch(p)               = 0._r8
          this%leafcn_offset_patch(p) = spval
          this%plantCN_patch(p)       = spval
       end do
    end if

  end subroutine Restart

end module CNVegStateType
