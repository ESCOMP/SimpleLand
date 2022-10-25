module CNVegCarbonStateType

#include "shr_assert.h"

  !-----------------------------------------------------------------------
  !
  ! !USES:
  use shr_kind_mod   , only : r8 => shr_kind_r8
  use shr_infnan_mod , only : nan => shr_infnan_nan, assignment(=)
  use shr_const_mod  , only : SHR_CONST_PDB
  use shr_log_mod    , only : errMsg => shr_log_errMsg
  use pftconMod	     , only : noveg, pftcon
  use clm_varcon     , only : spval, c3_r2, c4_r2
  use clm_varctl     , only : iulog
  use decompMod      , only : bounds_type
  use abortutils     , only : endrun
  use spmdMod        , only : masterproc 
  use LandunitType   , only : lun                
  use ColumnType     , only : col                
  use PatchType      , only : patch
  use CNSpeciesMod   , only : species_from_string, CN_SPECIES_C12
  use CNVegComputeSeedMod, only : ComputeSeedAmounts
  ! 
  ! !PUBLIC TYPES:
  implicit none
  private
  !

  type, public :: cnveg_carbonstate_type

     integer :: species  ! c12

     real(r8), pointer :: grainc_patch             (:) ! (gC/m2) grain C (crop model)
     real(r8), pointer :: grainc_storage_patch     (:) ! (gC/m2) grain C storage (crop model)
     real(r8), pointer :: grainc_xfer_patch        (:) ! (gC/m2) grain C transfer (crop model)
     real(r8), pointer :: leafc_patch              (:) ! (gC/m2) leaf C
     real(r8), pointer :: leafc_storage_patch      (:) ! (gC/m2) leaf C storage
     real(r8), pointer :: leafc_xfer_patch         (:) ! (gC/m2) leaf C transfer
     real(r8), pointer :: leafc_storage_xfer_acc_patch   (:) ! (gC/m2) Accmulated leaf C transfer
     real(r8), pointer :: storage_cdemand_patch          (:) ! (gC/m2)       C use from the C storage pool 
     real(r8), pointer :: frootc_patch             (:) ! (gC/m2) fine root C
     real(r8), pointer :: frootc_storage_patch     (:) ! (gC/m2) fine root C storage
     real(r8), pointer :: frootc_xfer_patch        (:) ! (gC/m2) fine root C transfer
     real(r8), pointer :: livestemc_patch          (:) ! (gC/m2) live stem C
     real(r8), pointer :: livestemc_storage_patch  (:) ! (gC/m2) live stem C storage
     real(r8), pointer :: livestemc_xfer_patch     (:) ! (gC/m2) live stem C transfer
     real(r8), pointer :: deadstemc_patch          (:) ! (gC/m2) dead stem C
     real(r8), pointer :: deadstemc_storage_patch  (:) ! (gC/m2) dead stem C storage
     real(r8), pointer :: deadstemc_xfer_patch     (:) ! (gC/m2) dead stem C transfer
     real(r8), pointer :: livecrootc_patch         (:) ! (gC/m2) live coarse root C
     real(r8), pointer :: livecrootc_storage_patch (:) ! (gC/m2) live coarse root C storage
     real(r8), pointer :: livecrootc_xfer_patch    (:) ! (gC/m2) live coarse root C transfer
     real(r8), pointer :: deadcrootc_patch         (:) ! (gC/m2) dead coarse root C
     real(r8), pointer :: deadcrootc_storage_patch (:) ! (gC/m2) dead coarse root C storage
     real(r8), pointer :: deadcrootc_xfer_patch    (:) ! (gC/m2) dead coarse root C transfer
     real(r8), pointer :: gresp_storage_patch      (:) ! (gC/m2) growth respiration storage
     real(r8), pointer :: gresp_xfer_patch         (:) ! (gC/m2) growth respiration transfer
     real(r8), pointer :: cpool_patch              (:) ! (gC/m2) temporary photosynthate C pool
     real(r8), pointer :: xsmrpool_patch           (:) ! (gC/m2) abstract C pool to meet excess MR demand
     real(r8), pointer :: ctrunc_patch             (:) ! (gC/m2) patch-level sink for C truncation
     real(r8), pointer :: woodc_patch              (:) ! (gC/m2) wood C
     real(r8), pointer :: leafcmax_patch           (:) ! (gC/m2) ann max leaf C
     real(r8), pointer :: totc_patch               (:) ! (gC/m2) total patch-level carbon, including cpool
     real(r8), pointer :: rootc_col                (:) ! (gC/m2) root carbon at column level (fire)
     real(r8), pointer :: leafc_col                (:) ! (gC/m2) column-level leafc (fire)
     real(r8), pointer :: deadstemc_col            (:) ! (gC/m2) column-level deadstemc (fire)
     real(r8), pointer :: fuelc_col                (:) ! fuel load outside cropland
     real(r8), pointer :: fuelc_crop_col           (:) ! fuel load for cropland
     real(r8), pointer :: cropseedc_deficit_patch  (:) ! (gC/m2) pool for seeding new crop growth; this is a NEGATIVE term, indicating the amount of seed usage that needs to be repaid

     ! pools for dynamic landcover
     real(r8), pointer :: seedc_grc                (:) ! (gC/m2) gridcell-level pool for seeding new PFTs via dynamic landcover

     ! summary (diagnostic) state variables, not involved in mass balance
     real(r8), pointer :: dispvegc_patch           (:) ! (gC/m2) displayed veg carbon, excluding storage and cpool
     real(r8), pointer :: storvegc_patch           (:) ! (gC/m2) stored vegetation carbon, excluding cpool
     real(r8), pointer :: totvegc_patch            (:) ! (gC/m2) total vegetation carbon, excluding cpool
     real(r8), pointer :: totvegc_col              (:) ! (gC/m2) total vegetation carbon, excluding cpool averaged to column (p2c)

     ! Total C pools       
     real(r8), pointer :: totc_p2c_col             (:) ! (gC/m2) totc_patch averaged to col
     real(r8), pointer :: totc_col                 (:) ! (gC/m2) total column carbon, incl veg and cpool
     real(r8), pointer :: totecosysc_col           (:) ! (gC/m2) total ecosystem carbon, incl veg but excl cpool 

   contains

     procedure , public  :: Init   
     procedure , public  :: SetValues
     procedure , public  :: Restart
     
     procedure , private :: InitAllocate    ! Allocate arrays
     procedure , private :: InitReadNML     ! Read in namelist
     procedure , private :: InitHistory     ! Initialize history
     procedure , private :: InitCold        ! Initialize arrays for a cold-start

  end type cnveg_carbonstate_type

  ! !PRIVATE DATA:

  type, private :: cnvegcarbonstate_const_type
      ! !PRIVATE MEMBER DATA:
      real(r8) :: initial_vegC = 20._r8    ! Initial vegetation carbon for leafc/frootc and storage
  end type
  type(cnvegcarbonstate_const_type), private :: cnvegcstate_const    ! Constants used here
  character(len=*), parameter :: sourcefile = &
       __FILE__

  !------------------------------------------------------------------------

contains

  !------------------------------------------------------------------------
  subroutine Init(this, bounds, carbon_type, ratio, NLFilename, &
                  c12_cnveg_carbonstate_inst)

    class(cnveg_carbonstate_type)                       :: this
    type(bounds_type)            , intent(in)           :: bounds  
    real(r8)                     , intent(in)           :: ratio
    character(len=*)             , intent(in)           :: carbon_type                ! Carbon isotope type C12, C13 or C1
    character(len=*)             , intent(in)           :: NLFilename                 ! Namelist filename
    type(cnveg_carbonstate_type) , intent(in), optional :: c12_cnveg_carbonstate_inst ! cnveg_carbonstate for C12 (if C13 or C14)
    !-----------------------------------------------------------------------

    this%species = species_from_string(carbon_type)

    call this%InitAllocate ( bounds)
    call this%InitReadNML  ( NLFilename )
    call this%InitHistory ( bounds, carbon_type)
    if (present(c12_cnveg_carbonstate_inst)) then
       call this%InitCold  ( bounds, ratio, carbon_type, c12_cnveg_carbonstate_inst )
    else
       call this%InitCold  ( bounds, ratio, carbon_type )
    end if

  end subroutine Init

  !------------------------------------------------------------------------
  subroutine InitReadNML(this, NLFilename)
    !
    ! !DESCRIPTION:
    ! Read the namelist for CNVegCarbonState
    !
    !USES:
    use fileutils      , only : getavu, relavu, opnfil
    use shr_nl_mod     , only : shr_nl_find_group_name
    use spmdMod        , only : masterproc, mpicom
    use shr_mpi_mod    , only : shr_mpi_bcast
    use clm_varctl     , only : iulog
    !
    ! !ARGUMENTS:
    class(cnveg_carbonstate_type)                       :: this
    character(len=*)             , intent(in)           :: NLFilename                 ! Namelist filename
    !
    ! !LOCAL VARIABLES:
    integer :: ierr                 ! error code
    integer :: unitn                ! unit for namelist file

    character(len=*), parameter :: subname = 'InitReadNML'
    character(len=*), parameter :: nmlname = 'cnvegcarbonstate'   ! MUST match what is in namelist below
    !-----------------------------------------------------------------------
    real(r8) :: initial_vegC
    namelist /cnvegcarbonstate/ initial_vegC

    initial_vegC = cnvegcstate_const%initial_vegC

    if (masterproc) then
       unitn = getavu()
       write(iulog,*) 'Read in '//nmlname//'  namelist'
       call opnfil (NLFilename, unitn, 'F')
       call shr_nl_find_group_name(unitn, nmlname, status=ierr)
       if (ierr == 0) then
          read(unitn, nml=cnvegcarbonstate, iostat=ierr)
          if (ierr /= 0) then
             call endrun(msg="ERROR reading "//nmlname//"namelist"//errmsg(sourcefile, __LINE__))
          end if
       else
          call endrun(msg="ERROR could NOT find "//nmlname//"namelist"//errmsg(sourcefile, __LINE__))
       end if
       call relavu( unitn )
    end if

    call shr_mpi_bcast (initial_vegC            , mpicom)

    cnvegcstate_const%initial_vegC = initial_vegC

    if (masterproc) then
       write(iulog,*) ' '
       write(iulog,*) nmlname//' settings:'
       write(iulog,nml=cnvegcarbonstate)    ! Name here MUST be the same as in nmlname above!
       write(iulog,*) ' '
    end if

    !-----------------------------------------------------------------------

  end subroutine InitReadNML

  !------------------------------------------------------------------------
  subroutine InitAllocate(this, bounds)
    !
    ! !ARGUMENTS:
    class (cnveg_carbonstate_type) :: this
    type(bounds_type), intent(in) :: bounds  
    !
    ! !LOCAL VARIABLES:
    integer           :: begp,endp
    integer           :: begc,endc
    integer           :: begg,endg
    !------------------------------------------------------------------------

    begp = bounds%begp; endp = bounds%endp
    begc = bounds%begc; endc = bounds%endc
    begg = bounds%begg; endg = bounds%endg

    allocate(this%leafc_patch              (begp:endp)) ; this%leafc_patch              (:) = nan
    allocate(this%leafc_storage_patch      (begp:endp)) ; this%leafc_storage_patch      (:) = nan
    allocate(this%leafc_xfer_patch         (begp:endp)) ; this%leafc_xfer_patch         (:) = nan
    allocate(this%leafc_storage_xfer_acc_patch (begp:endp)) ; this%leafc_storage_xfer_acc_patch (:) = nan
    allocate(this%storage_cdemand_patch        (begp:endp)) ; this%storage_cdemand_patch        (:) = nan
    allocate(this%frootc_patch             (begp:endp)) ; this%frootc_patch             (:) = nan
    allocate(this%frootc_storage_patch     (begp:endp)) ; this%frootc_storage_patch     (:) = nan
    allocate(this%frootc_xfer_patch        (begp:endp)) ; this%frootc_xfer_patch        (:) = nan
    allocate(this%livestemc_patch          (begp:endp)) ; this%livestemc_patch          (:) = nan
    allocate(this%livestemc_storage_patch  (begp:endp)) ; this%livestemc_storage_patch  (:) = nan
    allocate(this%livestemc_xfer_patch     (begp:endp)) ; this%livestemc_xfer_patch     (:) = nan
    allocate(this%deadstemc_patch          (begp:endp)) ; this%deadstemc_patch          (:) = nan
    allocate(this%deadstemc_storage_patch  (begp:endp)) ; this%deadstemc_storage_patch  (:) = nan
    allocate(this%deadstemc_xfer_patch     (begp:endp)) ; this%deadstemc_xfer_patch     (:) = nan
    allocate(this%livecrootc_patch         (begp:endp)) ; this%livecrootc_patch         (:) = nan
    allocate(this%livecrootc_storage_patch (begp:endp)) ; this%livecrootc_storage_patch (:) = nan
    allocate(this%livecrootc_xfer_patch    (begp:endp)) ; this%livecrootc_xfer_patch    (:) = nan
    allocate(this%deadcrootc_patch         (begp:endp)) ; this%deadcrootc_patch         (:) = nan
    allocate(this%deadcrootc_storage_patch (begp:endp)) ; this%deadcrootc_storage_patch (:) = nan
    allocate(this%deadcrootc_xfer_patch    (begp:endp)) ; this%deadcrootc_xfer_patch    (:) = nan
    allocate(this%gresp_storage_patch      (begp:endp)) ; this%gresp_storage_patch      (:) = nan
    allocate(this%gresp_xfer_patch         (begp:endp)) ; this%gresp_xfer_patch         (:) = nan
    allocate(this%cpool_patch              (begp:endp)) ; this%cpool_patch              (:) = nan
    allocate(this%xsmrpool_patch           (begp:endp)) ; this%xsmrpool_patch           (:) = nan
    allocate(this%ctrunc_patch             (begp:endp)) ; this%ctrunc_patch             (:) = nan
    allocate(this%dispvegc_patch           (begp:endp)) ; this%dispvegc_patch           (:) = nan
    allocate(this%storvegc_patch           (begp:endp)) ; this%storvegc_patch           (:) = nan
    allocate(this%leafcmax_patch           (begp:endp)) ; this%leafcmax_patch           (:) = nan
    allocate(this%totc_patch               (begp:endp))  ; this%totc_patch               (:) = nan
    allocate(this%grainc_patch             (begp:endp)) ; this%grainc_patch             (:) = nan
    allocate(this%grainc_storage_patch     (begp:endp)) ; this%grainc_storage_patch     (:) = nan
    allocate(this%grainc_xfer_patch        (begp:endp)) ; this%grainc_xfer_patch        (:) = nan
    allocate(this%woodc_patch              (begp:endp)) ; this%woodc_patch              (:) = nan     

    allocate(this%cropseedc_deficit_patch  (begp:endp)) ; this%cropseedc_deficit_patch  (:) = nan
    allocate(this%seedc_grc                (begg:endg)) ; this%seedc_grc                (:) = nan
    allocate(this%rootc_col                (begc:endc)) ; this%rootc_col                (:) = nan
    allocate(this%leafc_col                (begc:endc)) ; this%leafc_col                (:) = nan
    allocate(this%deadstemc_col            (begc:endc)) ; this%deadstemc_col            (:) = nan
    allocate(this%fuelc_col                (begc:endc)) ; this%fuelc_col                (:) = nan
    allocate(this%fuelc_crop_col           (begc:endc)) ; this%fuelc_crop_col           (:) = nan

    allocate(this%totvegc_patch            (begp:endp)) ; this%totvegc_patch            (:) = nan
    allocate(this%totvegc_col              (begc:endc)) ; this%totvegc_col              (:) = nan

    allocate(this%totc_p2c_col             (begc:endc)) ; this%totc_p2c_col             (:) = nan
    allocate(this%totc_col                 (begc:endc)) ; this%totc_col                 (:) = nan
    allocate(this%totecosysc_col           (begc:endc)) ; this%totecosysc_col           (:) = nan

  end subroutine InitAllocate

  !------------------------------------------------------------------------
  subroutine InitHistory(this, bounds, carbon_type)
    !
    ! !DESCRIPTION:
    ! add history fields for all CN variables, always set as default='inactive'
    !
    ! !USES:
    use histFileMod, only : hist_addfld1d, hist_addfld2d, hist_addfld_decomp 
    !
    ! !ARGUMENTS:
    class (cnveg_carbonstate_type) :: this
    type(bounds_type)         , intent(in) :: bounds 
    character(len=*)          , intent(in) :: carbon_type ! one of ['c12', c13','c14']
    !
    ! !LOCAL VARIABLES:
    integer           :: k,l,ii,jj 
    character(10)     :: active
    integer           :: begp,endp
    integer           :: begc,endc
    integer           :: begg,endg 
    character(24)     :: fieldname
    character(100)    :: longname
    real(r8), pointer :: data1dptr(:)   ! temp. pointer for slicing larger arrays
    real(r8), pointer :: data2dptr(:,:) ! temp. pointer for slicing larger arrays
    !---------------------------------------------------------------------

    begp = bounds%begp; endp = bounds%endp
    begc = bounds%begc; endc = bounds%endc
    begg = bounds%begg; endg = bounds%endg

    !-------------------------------
    ! C12 state variables
    !-------------------------------

       this%woodc_patch(begp:endp) = spval
       call hist_addfld1d (fname='WOODC', units='gC/m^2', &
            avgflag='A', long_name='wood C', &
            ptr_patch=this%woodc_patch, default='inactive')

       this%leafc_patch(begp:endp) = spval
       call hist_addfld1d (fname='LEAFC', units='gC/m^2', &
            avgflag='A', long_name='leaf C', &
            ptr_patch=this%leafc_patch, default='inactive')

       this%leafc_storage_patch(begp:endp) = spval
       call hist_addfld1d (fname='LEAFC_STORAGE', units='gC/m^2', &
            avgflag='A', long_name='leaf C storage', &
            ptr_patch=this%leafc_storage_patch, default='inactive')    

       this%leafc_xfer_patch(begp:endp) = spval
       call hist_addfld1d (fname='LEAFC_XFER', units='gC/m^2', &
            avgflag='A', long_name='leaf C transfer', &
            ptr_patch=this%leafc_xfer_patch, default='inactive')    

       this%leafc_storage_xfer_acc_patch(begp:endp) = spval
       call hist_addfld1d (fname='LEAFC_STORAGE_XFER_ACC', units='gC/m^2', &
            avgflag='A', long_name='Accumulated leaf C transfer', &
            ptr_patch=this%leafc_storage_xfer_acc_patch, default='inactive')

       this%storage_cdemand_patch(begp:endp) = spval
       call hist_addfld1d (fname='STORAGE_CDEMAND', units='gC/m^2', &
            avgflag='A', long_name='C use from the C storage pool', &
            ptr_patch=this%storage_cdemand_patch, default='inactive')

       this%frootc_patch(begp:endp) = spval
       call hist_addfld1d (fname='FROOTC', units='gC/m^2', &
            avgflag='A', long_name='fine root C', &
            ptr_patch=this%frootc_patch, default='inactive')

       this%frootc_storage_patch(begp:endp) = spval
       call hist_addfld1d (fname='FROOTC_STORAGE', units='gC/m^2', &
            avgflag='A', long_name='fine root C storage', &
            ptr_patch=this%frootc_storage_patch, default='inactive')   

       this%frootc_xfer_patch(begp:endp) = spval
       call hist_addfld1d (fname='FROOTC_XFER', units='gC/m^2', &
            avgflag='A', long_name='fine root C transfer', &
            ptr_patch=this%frootc_xfer_patch, default='inactive')    

       this%livestemc_patch(begp:endp) = spval
       call hist_addfld1d (fname='LIVESTEMC', units='gC/m^2', &
            avgflag='A', long_name='live stem C', &
            ptr_patch=this%livestemc_patch, default='inactive')

       this%livestemc_storage_patch(begp:endp) = spval
       call hist_addfld1d (fname='LIVESTEMC_STORAGE', units='gC/m^2', &
            avgflag='A', long_name='live stem C storage', &
            ptr_patch=this%livestemc_storage_patch, default='inactive')    

       this%livestemc_xfer_patch(begp:endp) = spval
       call hist_addfld1d (fname='LIVESTEMC_XFER', units='gC/m^2', &
            avgflag='A', long_name='live stem C transfer', &
            ptr_patch=this%livestemc_xfer_patch, default='inactive')     

       this%deadstemc_patch(begp:endp) = spval
       call hist_addfld1d (fname='DEADSTEMC', units='gC/m^2', &
            avgflag='A', long_name='dead stem C', &
            ptr_patch=this%deadstemc_patch, default='inactive')

       this%deadstemc_storage_patch(begp:endp) = spval
       call hist_addfld1d (fname='DEADSTEMC_STORAGE', units='gC/m^2', &
            avgflag='A', long_name='dead stem C storage', &
            ptr_patch=this%deadstemc_storage_patch, default='inactive')    

       this%deadstemc_xfer_patch(begp:endp) = spval
       call hist_addfld1d (fname='DEADSTEMC_XFER', units='gC/m^2', &
            avgflag='A', long_name='dead stem C transfer', &
            ptr_patch=this%deadstemc_xfer_patch, default='inactive')    

       this%livecrootc_patch(begp:endp) = spval
       call hist_addfld1d (fname='LIVECROOTC', units='gC/m^2', &
            avgflag='A', long_name='live coarse root C', &
            ptr_patch=this%livecrootc_patch, default='inactive')

       this%livecrootc_storage_patch(begp:endp) = spval
       call hist_addfld1d (fname='LIVECROOTC_STORAGE', units='gC/m^2', &
            avgflag='A', long_name='live coarse root C storage', &
            ptr_patch=this%livecrootc_storage_patch, default='inactive')     

       this%livecrootc_xfer_patch(begp:endp) = spval
       call hist_addfld1d (fname='LIVECROOTC_XFER', units='gC/m^2', &
            avgflag='A', long_name='live coarse root C transfer', &
            ptr_patch=this%livecrootc_xfer_patch, default='inactive')    

       this%deadcrootc_patch(begp:endp) = spval
       call hist_addfld1d (fname='DEADCROOTC', units='gC/m^2', &
            avgflag='A', long_name='dead coarse root C', &
            ptr_patch=this%deadcrootc_patch, default='inactive')

       this%deadcrootc_storage_patch(begp:endp) = spval
       call hist_addfld1d (fname='DEADCROOTC_STORAGE', units='gC/m^2', &
            avgflag='A', long_name='dead coarse root C storage', &
            ptr_patch=this%deadcrootc_storage_patch, default='inactive')   

       this%deadcrootc_xfer_patch(begp:endp) = spval
       call hist_addfld1d (fname='DEADCROOTC_XFER', units='gC/m^2', &
            avgflag='A', long_name='dead coarse root C transfer', &
            ptr_patch=this%deadcrootc_xfer_patch, default='inactive')   

       this%gresp_storage_patch(begp:endp) = spval
       call hist_addfld1d (fname='GRESP_STORAGE', units='gC/m^2', &
            avgflag='A', long_name='growth respiration storage', &
            ptr_patch=this%gresp_storage_patch, default='inactive')    

       this%gresp_xfer_patch(begp:endp) = spval
       call hist_addfld1d (fname='GRESP_XFER', units='gC/m^2', &
            avgflag='A', long_name='growth respiration transfer', &
            ptr_patch=this%gresp_xfer_patch, default='inactive')     

       this%cpool_patch(begp:endp) = spval
       call hist_addfld1d (fname='CPOOL', units='gC/m^2', &
            avgflag='A', long_name='temporary photosynthate C pool', &
            ptr_patch=this%cpool_patch, default='inactive')

       this%xsmrpool_patch(begp:endp) = spval
       call hist_addfld1d (fname='XSMRPOOL', units='gC/m^2', &
            avgflag='A', long_name='temporary photosynthate C pool', &
            ptr_patch=this%xsmrpool_patch, default='inactive')

       this%ctrunc_patch(begp:endp) = spval
       call hist_addfld1d (fname='PFT_CTRUNC', units='gC/m^2', &
            avgflag='A', long_name='patch-level sink for C truncation', &
            ptr_patch=this%ctrunc_patch, default='inactive')

       this%dispvegc_patch(begp:endp) = spval
       call hist_addfld1d (fname='DISPVEGC', units='gC/m^2', &
            avgflag='A', long_name='displayed veg carbon, excluding storage and cpool', &
            ptr_patch=this%dispvegc_patch, default='inactive')

       this%storvegc_patch(begp:endp) = spval
       call hist_addfld1d (fname='STORVEGC', units='gC/m^2', &
            avgflag='A', long_name='stored vegetation carbon, excluding cpool', &
            ptr_patch=this%storvegc_patch, default='inactive')

       this%totvegc_patch(begp:endp) = spval
       call hist_addfld1d (fname='TOTVEGC', units='gC/m^2', &
            avgflag='A', long_name='total vegetation carbon, excluding cpool', &
            ptr_patch=this%totvegc_patch, default='inactive')

       this%totc_patch(begp:endp) = spval
       call hist_addfld1d (fname='TOTPFTC', units='gC/m^2', &
            avgflag='A', long_name='total patch-level carbon, including cpool', &
            ptr_patch=this%totc_patch, default='inactive')

       this%seedc_grc(begg:endg) = spval
       call hist_addfld1d (fname='SEEDC', units='gC/m^2', &
            avgflag='A', long_name='pool for seeding new PFTs via dynamic landcover', &
            ptr_gcell=this%seedc_grc, default='inactive')

       this%fuelc_col(begc:endc) = spval
       call hist_addfld1d (fname='FUELC', units='gC/m^2', &
            avgflag='A', long_name='fuel load', &
            ptr_col=this%fuelc_col, default='inactive')

       this%totc_col(begc:endc) = spval
       call hist_addfld1d (fname='TOTCOLC', units='gC/m^2', &
            avgflag='A', long_name='total column carbon, incl veg and cpool but excl product pools', &
            ptr_col=this%totc_col, default='inactive')

       this%totecosysc_col(begc:endc) = spval
       call hist_addfld1d (fname='TOTECOSYSC', units='gC/m^2', &
            avgflag='A', long_name='total ecosystem carbon, incl veg but excl cpool and product pools', &
            ptr_col=this%totecosysc_col, default='inactive')

  end subroutine InitHistory

  !-----------------------------------------------------------------------
  subroutine InitCold(this, bounds, ratio, carbon_type, c12_cnveg_carbonstate_inst)
    !
    ! !DESCRIPTION:
    ! Initializes time varying variables used only in coupled carbon-nitrogen mode (CN):
    !
    ! !USES, default='inactive':
    use landunit_varcon	 , only : istsoil, istcrop 
    use clm_time_manager , only : is_restart, get_nstep
    !
    ! !ARGUMENTS:
    class(cnveg_carbonstate_type)                       :: this 
    type(bounds_type)            , intent(in)           :: bounds  
    real(r8)                     , intent(in)           :: ratio              ! Standard isotope ratio
    character(len=*)             , intent(in)           :: carbon_type        ! 'c12' or 'c13' or 'c14'
    type(cnveg_carbonstate_type) , optional, intent(in) :: c12_cnveg_carbonstate_inst
    !
    ! !LOCAL VARIABLES:
    integer  :: p,c,l,g,j,k,i
    integer  :: fc                                       ! filter index
    integer  :: num_special_col                          ! number of good values in special_col filter
    integer  :: num_special_patch                        ! number of good values in special_patch filter
    integer  :: special_col(bounds%endc-bounds%begc+1)   ! special landunit filter - columns
    integer  :: special_patch(bounds%endp-bounds%begp+1) ! special landunit filter - patches
    !-----------------------------------------------------------------------

    ! Set column filters

    num_special_col = 0
    do c = bounds%begc, bounds%endc
       l = col%landunit(c)
       if (lun%ifspecial(l)) then
          num_special_col = num_special_col + 1
          special_col(num_special_col) = c
       end if
    end do

    ! Set patch filters

    num_special_patch = 0
    do p = bounds%begp,bounds%endp
       l = patch%landunit(p)
       if (lun%ifspecial(l)) then
          num_special_patch = num_special_patch + 1
          special_patch(num_special_patch) = p
       end if
    end do

    !-----------------------------------------------
    ! initialize patch-level carbon state variables
    !-----------------------------------------------

    do p = bounds%begp,bounds%endp

       this%leafcmax_patch(p) = 0._r8

       l = patch%landunit(p)
       if (lun%itype(l) == istsoil .or. lun%itype(l) == istcrop) then

          if (patch%itype(p) == noveg) then
             this%leafc_patch(p)          = 0._r8
             this%leafc_storage_patch(p)  = 0._r8
             this%frootc_patch(p)         = 0._r8            
             this%frootc_storage_patch(p) = 0._r8    
          else
             if (pftcon%evergreen(patch%itype(p)) == 1._r8) then
                this%leafc_patch(p)          = cnvegcstate_const%initial_vegC * ratio     
                this%leafc_storage_patch(p)  = 0._r8
                this%frootc_patch(p)         = cnvegcstate_const%initial_vegC * ratio           
                this%frootc_storage_patch(p) = 0._r8    
             else
                this%leafc_patch(p)          = 0._r8
                this%leafc_storage_patch(p)  = cnvegcstate_const%initial_vegC * ratio   
                this%frootc_patch(p)         = 0._r8            
                this%frootc_storage_patch(p) = cnvegcstate_const%initial_vegC * ratio   
             end if
          end if
          this%leafc_xfer_patch(p) = 0._r8
          this%leafc_storage_xfer_acc_patch(p)  = 0._r8
          this%storage_cdemand_patch(p)         = 0._r8

          this%frootc_patch(p)            = 0._r8 
          this%frootc_storage_patch(p)    = 0._r8 
          this%frootc_xfer_patch(p)       = 0._r8 

          this%livestemc_patch(p)         = 0._r8 
          this%livestemc_storage_patch(p) = 0._r8 
          this%livestemc_xfer_patch(p)    = 0._r8 

          if (pftcon%woody(patch%itype(p)) == 1._r8) then
             this%deadstemc_patch(p) = 0.1_r8 * ratio
          else
             this%deadstemc_patch(p) = 0._r8 
          end if
          this%deadstemc_storage_patch(p)  = 0._r8 
          this%deadstemc_xfer_patch(p)     = 0._r8 

          this%livecrootc_patch(p)         = 0._r8 
          this%livecrootc_storage_patch(p) = 0._r8 
          this%livecrootc_xfer_patch(p)    = 0._r8 

          this%deadcrootc_patch(p)         = 0._r8 
          this%deadcrootc_storage_patch(p) = 0._r8 
          this%deadcrootc_xfer_patch(p)    = 0._r8 

          this%gresp_storage_patch(p)      = 0._r8 
          this%gresp_xfer_patch(p)         = 0._r8 

          this%cpool_patch(p)              = 0._r8 
          this%xsmrpool_patch(p)           = 0._r8 
          this%ctrunc_patch(p)             = 0._r8 
          this%dispvegc_patch(p)           = 0._r8 
          this%storvegc_patch(p)           = 0._r8 
          this%woodc_patch(p)              = 0._r8
          this%totc_patch(p)               = 0._r8 

       endif

    end do

    ! -----------------------------------------------
    ! initialize column-level variables
    ! -----------------------------------------------

    do c = bounds%begc, bounds%endc
       l = col%landunit(c)
       if (lun%itype(l) == istsoil .or. lun%itype(l) == istcrop) then
!          this%totgrainc_col(c)  = 0._r8

          ! total carbon pools
          this%totecosysc_col(c) = 0._r8
          this%totc_p2c_col(c)   = 0._r8
          this%totc_col(c)       = 0._r8
       end if
    end do


    do g = bounds%begg, bounds%endg
       this%seedc_grc(g) = 0._r8
    end do

    if ( .not. is_restart() .and. get_nstep() == 1 ) then

       do p = bounds%begp,bounds%endp
          if (pftcon%c3psn(patch%itype(p)) == 1._r8) then
             this%grainc_patch(p)            = c12_cnveg_carbonstate_inst%grainc_patch(p)         * c3_r2
             this%grainc_storage_patch(p)    = c12_cnveg_carbonstate_inst%grainc_storage_patch(p) * c3_r2
             this%grainc_xfer_patch(p)       = c12_cnveg_carbonstate_inst%grainc_xfer_patch(p)    * c3_r2
             this%dispvegc_patch(p)          = c12_cnveg_carbonstate_inst%dispvegc_patch(p)       * c3_r2
             this%storvegc_patch(p)          = c12_cnveg_carbonstate_inst%storvegc_patch(p)       * c3_r2
             this%totvegc_patch(p)           = c12_cnveg_carbonstate_inst%totvegc_patch(p)        * c3_r2
             this%totc_patch(p)              = c12_cnveg_carbonstate_inst%totc_patch(p)           * c3_r2
             this%woodc_patch(p)             = c12_cnveg_carbonstate_inst%woodc_patch(p)          * c3_r2
          else
             this%grainc_patch(p)            = c12_cnveg_carbonstate_inst%grainc_patch(p)         * c4_r2
             this%grainc_storage_patch(p)    = c12_cnveg_carbonstate_inst%grainc_storage_patch(p) * c4_r2
             this%grainc_xfer_patch(p)       = c12_cnveg_carbonstate_inst%grainc_xfer_patch(p)    * c4_r2
             this%dispvegc_patch(p)          = c12_cnveg_carbonstate_inst%dispvegc_patch(p)       * c4_r2
             this%storvegc_patch(p)          = c12_cnveg_carbonstate_inst%storvegc_patch(p)       * c4_r2
             this%totvegc_patch(p)           = c12_cnveg_carbonstate_inst%totvegc_patch(p)        * c4_r2
             this%totc_patch(p)              = c12_cnveg_carbonstate_inst%totc_patch(p)           * c4_r2
             this%woodc_patch(p)             = c12_cnveg_carbonstate_inst%woodc_patch(p)          * c4_r2
          end if
       end do
    end if

    ! initialize fields for special filters

    call this%SetValues (&
         num_patch=num_special_patch, filter_patch=special_patch, value_patch=0._r8, &
         num_column=num_special_col, filter_column=special_col, value_column=0._r8)

  end subroutine InitCold

  !-----------------------------------------------------------------------
  subroutine Restart ( this,  bounds, ncid, flag, carbon_type, reseed_dead_plants, &
                       c12_cnveg_carbonstate_inst, filter_reseed_patch, &
                       num_reseed_patch)
    !
    ! !DESCRIPTION: 
    ! Read/write CN restart data for carbon state
    !
    ! !USES:
    use shr_infnan_mod   , only : isnan => shr_infnan_isnan, nan => shr_infnan_nan, assignment(=)
    use clm_varctl       , only : spinup_state
    use clm_time_manager , only : get_nstep, is_restart, get_nstep
    use landunit_varcon	 , only : istsoil, istcrop 
    use spmdMod          , only : mpicom
    use shr_mpi_mod      , only : shr_mpi_sum
    use restUtilMod
    use ncdio_pio
    !
    ! !ARGUMENTS:
    class (cnveg_carbonstate_type)                               :: this
    type(bounds_type)                     , intent(in)           :: bounds 
    type(file_desc_t)                     , intent(inout)        :: ncid   ! netcdf id
    character(len=*)                      , intent(in)           :: flag   !'read' or 'write'
    character(len=*)                      , intent(in)           :: carbon_type ! 'c12' or 'c13' or 'c14'
    logical                               , intent(in)           :: reseed_dead_plants
    type (cnveg_carbonstate_type)         , intent(in), optional :: c12_cnveg_carbonstate_inst 
    integer                               , intent(out), optional :: filter_reseed_patch(:)
    integer                               , intent(out), optional :: num_reseed_patch
    !
    ! !LOCAL VARIABLES:
    integer            :: i,j,k,l,c,p
    real(r8)           :: ratio
    character(len=128) :: varname   ! temporary
    logical            :: readvar
    integer            :: idata
    logical            :: exit_spinup  = .false.
    logical            :: enter_spinup = .false.
    ! flags for comparing the model and restart decomposition cascades
    integer            :: decomp_cascade_state, restart_file_decomp_cascade_state 
    ! spinup state as read from restart file, for determining whether to enter or exit spinup mode.
    integer            :: restart_file_spinup_state
    integer            :: total_num_reseed_patch      ! Total number of patches to reseed across all processors

    !------------------------------------------------------------------------

    ratio = 1._r8

    if ( (      present(num_reseed_patch) .and. .not. present(filter_reseed_patch)) &
    .or. (.not. present(num_reseed_patch) .and.       present(filter_reseed_patch) ) )then
       call endrun(msg=' ERROR: filter_reseed_patch and num_reseed_patch both need to be entered ' //&
       errMsg(sourcefile, __LINE__))
    end if
    if ( present(num_reseed_patch) )then
       num_reseed_patch = 0
       filter_reseed_patch(:) = -1
    end if

    !--------------------------------
    ! patch carbon state variables (c12)
    !--------------------------------

       call restartvar(ncid=ncid, flag=flag, varname='leafc', xtype=ncd_double,  &
            dim1name='pft', long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%leafc_patch) 

       call restartvar(ncid=ncid, flag=flag, varname='leafc_storage', xtype=ncd_double,  &
            dim1name='pft', long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%leafc_storage_patch) 

       call restartvar(ncid=ncid, flag=flag, varname='leafc_xfer', xtype=ncd_double,  &
            dim1name='pft', long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%leafc_xfer_patch) 

       call restartvar(ncid=ncid, flag=flag, varname='leafc_storage_xfer_acc', xtype=ncd_double,  &
            dim1name='pft', long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%leafc_storage_xfer_acc_patch)
 
       call restartvar(ncid=ncid, flag=flag, varname='storage_cdemand', xtype=ncd_double,  &
            dim1name='pft', long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%storage_cdemand_patch)

       call restartvar(ncid=ncid, flag=flag, varname='frootc', xtype=ncd_double,  &
            dim1name='pft', long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%frootc_patch) 

       call restartvar(ncid=ncid, flag=flag, varname='frootc_storage', xtype=ncd_double,  &
            dim1name='pft', long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%frootc_storage_patch) 

       call restartvar(ncid=ncid, flag=flag, varname='frootc_xfer', xtype=ncd_double,  &
            dim1name='pft', long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%frootc_xfer_patch) 

       call restartvar(ncid=ncid, flag=flag, varname='livestemc', xtype=ncd_double,  &
            dim1name='pft', long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%livestemc_patch) 

       call restartvar(ncid=ncid, flag=flag, varname='livestemc_storage', xtype=ncd_double,  &
            dim1name='pft', long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%livestemc_storage_patch) 

       call restartvar(ncid=ncid, flag=flag, varname='livestemc_xfer', xtype=ncd_double,  &
            dim1name='pft', long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%livestemc_xfer_patch) 

       call restartvar(ncid=ncid, flag=flag, varname='deadstemc', xtype=ncd_double,  &
            dim1name='pft', long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%deadstemc_patch) 

       call restartvar(ncid=ncid, flag=flag, varname='deadstemc_storage', xtype=ncd_double,  &
            dim1name='pft', long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%deadstemc_storage_patch) 

       call restartvar(ncid=ncid, flag=flag, varname='deadstemc_xfer', xtype=ncd_double,  &
            dim1name='pft', long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%deadstemc_xfer_patch) 

       call restartvar(ncid=ncid, flag=flag, varname='livecrootc', xtype=ncd_double,  &
            dim1name='pft', long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%livecrootc_patch) 

       call restartvar(ncid=ncid, flag=flag, varname='livecrootc_storage', xtype=ncd_double,  &
            dim1name='pft', long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%livecrootc_storage_patch) 

       call restartvar(ncid=ncid, flag=flag, varname='livecrootc_xfer', xtype=ncd_double,  &
            dim1name='pft', long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%livecrootc_xfer_patch) 

       call restartvar(ncid=ncid, flag=flag, varname='deadcrootc', xtype=ncd_double,  &
            dim1name='pft', long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%deadcrootc_patch) 

       call restartvar(ncid=ncid, flag=flag, varname='deadcrootc_storage', xtype=ncd_double,  &
            dim1name='pft', long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%deadcrootc_storage_patch) 

       call restartvar(ncid=ncid, flag=flag, varname='deadcrootc_xfer', xtype=ncd_double,  &
            dim1name='pft', long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%deadcrootc_xfer_patch) 

       call restartvar(ncid=ncid, flag=flag, varname='gresp_storage', xtype=ncd_double,  &
            dim1name='pft', long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%gresp_storage_patch) 

       call restartvar(ncid=ncid, flag=flag, varname='gresp_xfer', xtype=ncd_double,  &
            dim1name='pft', long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%gresp_xfer_patch) 

       call restartvar(ncid=ncid, flag=flag, varname='cpool', xtype=ncd_double,  &
            dim1name='pft', long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%cpool_patch) 

       call restartvar(ncid=ncid, flag=flag, varname='xsmrpool', xtype=ncd_double,  &
            dim1name='pft', long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%xsmrpool_patch) 

       call restartvar(ncid=ncid, flag=flag, varname='pft_ctrunc', xtype=ncd_double,  &
            dim1name='pft', long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%ctrunc_patch) 

       call restartvar(ncid=ncid, flag=flag, varname='leafcmax', xtype=ncd_double,  &
            dim1name='pft', &
            long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%leafcmax_patch)

       if (flag == 'read') then
          call restartvar(ncid=ncid, flag=flag, varname='spinup_state', xtype=ncd_int, &
            long_name='Spinup state of the model that wrote this restart file: ' &
            // ' 0 = normal model mode, 1 = AD spinup, 2 = AAD spinup', units='', &
            interpinic_flag='copy', readvar=readvar,  data=idata)

          if (readvar) then
             restart_file_spinup_state = idata
          else
             restart_file_spinup_state = spinup_state
             if ( masterproc ) then
                write(iulog,*) ' CNRest: WARNING!  Restart file does not contain info ' &
                      // ' on spinup state used to generate the restart file. '
                 write(iulog,*) '   Assuming the same as current setting: ', spinup_state
             end if
          end if
       end if

       if (flag == 'read' .and. spinup_state /= restart_file_spinup_state) then
          if ( masterproc ) write(iulog, *) 'exit_spinup ',exit_spinup,' restart_file_spinup_state ',restart_file_spinup_state
          if (spinup_state <= 1 .and. restart_file_spinup_state == 2 ) then
             if ( masterproc ) write(iulog,*) ' CNRest: taking Dead wood C pools out of AD spinup mode'
             exit_spinup = .true.
             if ( masterproc ) write(iulog, *) 'Multiplying stemc and crootc by 10 for exit spinup'
             do i = bounds%begp,bounds%endp
                this%deadstemc_patch(i) = this%deadstemc_patch(i) * 10._r8
                this%deadcrootc_patch(i) = this%deadcrootc_patch(i) * 10._r8
             end do
          else if (spinup_state == 2 .and. restart_file_spinup_state <= 1 )then
             if (spinup_state == 2 .and. restart_file_spinup_state <= 1 )then
                if ( masterproc ) write(iulog,*) ' CNRest: taking Dead wood C pools into AD spinup mode'
                enter_spinup = .true.
                if ( masterproc ) write(iulog, *) 'Dividing stemc and crootc by 10 for enter spinup '
                do i = bounds%begp,bounds%endp
                   this%deadstemc_patch(i) = this%deadstemc_patch(i) / 10._r8
                   this%deadcrootc_patch(i) = this%deadcrootc_patch(i) / 10._r8
                end do
             end if
          end if

       !--------------------------------
       ! C12 carbon state variables
       !--------------------------------

          call restartvar(ncid=ncid, flag=flag, varname='totvegc', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%totvegc_patch) 
          ! totvegc_col needed for resetting soil carbon stocks during AD spinup exit
          call restartvar(ncid=ncid, flag=flag, varname='totvegc_col', xtype=ncd_double,  &
               dim1name='column', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%totvegc_col)


       if (  flag == 'read' .and. (enter_spinup .or. (reseed_dead_plants .and. .not. is_restart()))) then
             if ( masterproc ) write(iulog, *) 'Reseeding dead plants for CNVegCarbonState'
             ! If a pft is dead (indicated by totvegc = 0) then we reseed that
             ! pft according to the cold start protocol in the InitCold subroutine.
             ! Thus, the variable totvegc is required to be read before here
             ! so that if it is zero for a given pft, the pft can be reseeded.
             do i = bounds%begp,bounds%endp
                if (this%totvegc_patch(i) .le. 0.0_r8) then
                   !-----------------------------------------------
                   ! initialize patch-level carbon state variables
                   !-----------------------------------------------

                   this%leafcmax_patch(i) = 0._r8

                   l = patch%landunit(i)
                   if (lun%itype(l) == istsoil )then
                      if ( present(num_reseed_patch) ) then
                         num_reseed_patch = num_reseed_patch + 1
                         filter_reseed_patch(num_reseed_patch) = i
                      end if

                      if (patch%itype(i) == noveg) then
                         this%leafc_patch(i)          = 0._r8
                         this%leafc_storage_patch(i)  = 0._r8
                         this%frootc_patch(i)         = 0._r8            
                         this%frootc_storage_patch(i) = 0._r8    
                      else
                         if (pftcon%evergreen(patch%itype(i)) == 1._r8) then
                            this%leafc_patch(i)          = cnvegcstate_const%initial_vegC * ratio     
                            this%leafc_storage_patch(i)  = 0._r8
                            this%frootc_patch(i)         = cnvegcstate_const%initial_vegC * ratio           
                            this%frootc_storage_patch(i) = 0._r8    
                         else
                            this%leafc_patch(i)          = 0._r8
                            this%leafc_storage_patch(i)  = cnvegcstate_const%initial_vegC * ratio   
                            this%frootc_patch(i)         = 0._r8            
                            this%frootc_storage_patch(i) = cnvegcstate_const%initial_vegC * ratio   
                         end if
                      end if
                      this%leafc_xfer_patch(i) = 0._r8
                      this%leafc_storage_xfer_acc_patch(i)  = 0._r8
                      this%storage_cdemand_patch(i)         = 0._r8

                      this%frootc_patch(i)            = 0._r8 
                      this%frootc_storage_patch(i)    = 0._r8 
                      this%frootc_xfer_patch(i)       = 0._r8 

                      this%livestemc_patch(i)         = 0._r8 
                      this%livestemc_storage_patch(i) = 0._r8 
                      this%livestemc_xfer_patch(i)    = 0._r8 

                      if (pftcon%woody(patch%itype(i)) == 1._r8) then
                         this%deadstemc_patch(i) = 0.1_r8 * ratio
                      else
                         this%deadstemc_patch(i) = 0._r8 
                      end if
                      this%deadstemc_storage_patch(i)  = 0._r8 
                      this%deadstemc_xfer_patch(i)     = 0._r8 

                      this%livecrootc_patch(i)         = 0._r8 
                      this%livecrootc_storage_patch(i) = 0._r8 
                      this%livecrootc_xfer_patch(i)    = 0._r8 

                      this%deadcrootc_patch(i)         = 0._r8 
                      this%deadcrootc_storage_patch(i) = 0._r8 
                      this%deadcrootc_xfer_patch(i)    = 0._r8 

                      this%gresp_storage_patch(i)      = 0._r8 
                      this%gresp_xfer_patch(i)         = 0._r8 

                      this%cpool_patch(i)              = 0._r8 
                      this%xsmrpool_patch(i)           = 0._r8 
                      this%ctrunc_patch(i)             = 0._r8 
                      this%dispvegc_patch(i)           = 0._r8 
                      this%storvegc_patch(i)           = 0._r8 
                      this%woodc_patch(i)              = 0._r8
                      this%totc_patch(i)               = 0._r8 

                      ! calculate totvegc explicitly so that it is available for the isotope 
                      ! code on the first time step.

                      this%totvegc_patch(i) = &
                           this%leafc_patch(i)              + &
                           this%leafc_storage_patch(i)      + &
                           this%leafc_xfer_patch(i)         + &
                           this%frootc_patch(i)             + &
                           this%frootc_storage_patch(i)     + &
                           this%frootc_xfer_patch(i)        + &
                           this%livestemc_patch(i)          + &
                           this%livestemc_storage_patch(i)  + &
                           this%livestemc_xfer_patch(i)     + &
                           this%deadstemc_patch(i)          + &
                           this%deadstemc_storage_patch(i)  + &
                           this%deadstemc_xfer_patch(i)     + &
                           this%livecrootc_patch(i)         + &
                           this%livecrootc_storage_patch(i) + &
                           this%livecrootc_xfer_patch(i)    + &
                           this%deadcrootc_patch(i)         + &
                           this%deadcrootc_storage_patch(i) + &
                           this%deadcrootc_xfer_patch(i)    + &
                           this%gresp_storage_patch(i)      + &
                           this%gresp_xfer_patch(i)         + &
                           this%cpool_patch(i)
                   endif
                end if
             end do
             if ( .not. is_restart() .and. get_nstep() == 1 ) then

                do p = bounds%begp,bounds%endp
                  if (this%leafc_patch(p) .lt. 0.01_r8) then
                   if (pftcon%c3psn(patch%itype(p)) == 1._r8) then
                      this%grainc_patch(p)         = c12_cnveg_carbonstate_inst%grainc_patch(p)         * c3_r2
                      this%grainc_storage_patch(p) = c12_cnveg_carbonstate_inst%grainc_storage_patch(p) * c3_r2
                      this%grainc_xfer_patch(p)    = c12_cnveg_carbonstate_inst%grainc_xfer_patch(p)    * c3_r2
                      this%dispvegc_patch(p)       = c12_cnveg_carbonstate_inst%dispvegc_patch(p)       * c3_r2
                      this%storvegc_patch(p)       = c12_cnveg_carbonstate_inst%storvegc_patch(p)       * c3_r2
                      this%totvegc_patch(p)        = c12_cnveg_carbonstate_inst%totvegc_patch(p)        * c3_r2
                      this%totc_patch(p)           = c12_cnveg_carbonstate_inst%totc_patch(p)           * c3_r2
                      this%woodc_patch(p)          = c12_cnveg_carbonstate_inst%woodc_patch(p)          * c3_r2
                   else
                      this%grainc_patch(p)         = c12_cnveg_carbonstate_inst%grainc_patch(p)         * c4_r2
                      this%grainc_storage_patch(p) = c12_cnveg_carbonstate_inst%grainc_storage_patch(p) * c4_r2
                      this%grainc_xfer_patch(p)    = c12_cnveg_carbonstate_inst%grainc_xfer_patch(p)    * c4_r2
                      this%dispvegc_patch(p)       = c12_cnveg_carbonstate_inst%dispvegc_patch(p)       * c4_r2
                      this%storvegc_patch(p)       = c12_cnveg_carbonstate_inst%storvegc_patch(p)       * c4_r2
                      this%totvegc_patch(p)        = c12_cnveg_carbonstate_inst%totvegc_patch(p)        * c4_r2
                      this%totc_patch(p)           = c12_cnveg_carbonstate_inst%totc_patch(p)           * c4_r2
                      this%woodc_patch(p)          = c12_cnveg_carbonstate_inst%woodc_patch(p)          * c4_r2
                   end if
                  end if
                end do
             end if
             if ( present(num_reseed_patch) ) then
                call shr_mpi_sum( num_reseed_patch, total_num_reseed_patch, mpicom )
                if ( masterproc ) write(iulog,*) 'Total num_reseed, over all tasks = ', total_num_reseed_patch
             end if
       end if

    end if

  end subroutine Restart

  !-----------------------------------------------------------------------
  subroutine SetValues ( this, &
       num_patch, filter_patch, value_patch, &
       num_column, filter_column, value_column)
    !
    ! !DESCRIPTION:
    ! Set carbon state variables
    !
    ! !ARGUMENTS:
    class (cnveg_carbonstate_type) :: this
    integer , intent(in) :: num_patch
    integer , intent(in) :: filter_patch(:)
    real(r8), intent(in) :: value_patch
    integer , intent(in) :: num_column
    integer , intent(in) :: filter_column(:)
    real(r8), intent(in) :: value_column
    !
    ! !LOCAL VARIABLES:
    integer :: fi,i,j,k,l     ! loop index
    !------------------------------------------------------------------------

    do fi = 1,num_patch
       i  = filter_patch(fi)
       this%leafc_patch(i)              = value_patch
       this%leafc_storage_patch(i)      = value_patch
       this%leafc_xfer_patch(i)         = value_patch
       this%leafc_storage_xfer_acc_patch(i) = value_patch
       this%storage_cdemand_patch(i)        = value_patch        
       this%frootc_patch(i)             = value_patch
       this%frootc_storage_patch(i)     = value_patch
       this%frootc_xfer_patch(i)        = value_patch
       this%livestemc_patch(i)          = value_patch
       this%livestemc_storage_patch(i)  = value_patch
       this%livestemc_xfer_patch(i)     = value_patch
       this%deadstemc_patch(i)          = value_patch
       this%deadstemc_storage_patch(i)  = value_patch
       this%deadstemc_xfer_patch(i)     = value_patch
       this%livecrootc_patch(i)         = value_patch
       this%livecrootc_storage_patch(i) = value_patch
       this%livecrootc_xfer_patch(i)    = value_patch
       this%deadcrootc_patch(i)         = value_patch
       this%deadcrootc_storage_patch(i) = value_patch
       this%deadcrootc_xfer_patch(i)    = value_patch
       this%gresp_storage_patch(i)      = value_patch
       this%gresp_xfer_patch(i)         = value_patch
       this%cpool_patch(i)              = value_patch
       this%xsmrpool_patch(i)           = value_patch
       this%ctrunc_patch(i)             = value_patch
       this%dispvegc_patch(i)           = value_patch
       this%storvegc_patch(i)           = value_patch
       this%woodc_patch(i)              = value_patch
       this%totvegc_patch(i)            = value_patch
       this%totc_patch(i)               = value_patch
    end do

    do fi = 1,num_column
       i  = filter_column(fi)
       this%rootc_col(i)                = value_column
       this%leafc_col(i)                = value_column
       this%deadstemc_col(i)            = value_column
       this%fuelc_col(i)                = value_column
       this%fuelc_crop_col(i)           = value_column
       this%totvegc_col(i)              = value_column
       this%totc_p2c_col(i)             = value_column
       this%totc_col(i)                 = value_column
       this%totecosysc_col(i)           = value_column
    end do

  end subroutine SetValues

end module CNVegCarbonStateType
