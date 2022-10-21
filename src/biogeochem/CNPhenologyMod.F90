module CNPhenologyMod

#include "shr_assert.h"

  !-----------------------------------------------------------------------
  ! !MODULE: CNPhenologyMod
  !
  ! !DESCRIPTION:
  ! Module holding routines used in phenology model for coupled carbon
  ! nitrogen code.
  !
  ! !USES:
  use shr_kind_mod                    , only : r8 => shr_kind_r8
  use shr_log_mod                     , only : errMsg => shr_log_errMsg
  use shr_sys_mod                     , only : shr_sys_flush
  use decompMod                       , only : bounds_type
  use clm_varpar                      , only : numpft, nlevdecomp_full
  use clm_varctl                      , only : iulog
  use clm_varcon                      , only : tfrz
  use abortutils                      , only : endrun
  use CanopyStateType                 , only : canopystate_type
  use CNVegstateType                  , only : cnveg_state_type
  use CNVegCarbonStateType            , only : cnveg_carbonstate_type
  use CNVegCarbonFluxType             , only : cnveg_carbonflux_type
  use CNVegnitrogenstateType          , only : cnveg_nitrogenstate_type
  use CNVegnitrogenfluxType           , only : cnveg_nitrogenflux_type
  use CropType                        , only : crop_type
  use pftconMod                       , only : pftcon
  use SoilStateType                   , only : soilstate_type
  use TemperatureType                 , only : temperature_type
  use WaterstateType                  , only : waterstate_type
  use ColumnType                      , only : col                
  use GridcellType                    , only : grc                
  use PatchType                       , only : patch   
  use atm2lndType                     , only : atm2lnd_type             
  use atm2lndType                     , only : atm2lnd_type
  !
  implicit none
  private
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: readParams           ! Read parameters
  public :: CNPhenologyInit      ! Initialization
  !
  ! !PRIVATE DATA MEMBERS:
  type, private :: params_type
     real(r8) :: crit_dayl       ! critical day length for senescence
     real(r8) :: ndays_on     	 ! number of days to complete leaf onset
     real(r8) :: ndays_off	 ! number of days to complete leaf offset
     real(r8) :: fstor2tran      ! fraction of storage to move to transfer for each onset
     real(r8) :: crit_onset_fdd  ! critical number of freezing days to set gdd counter
     real(r8) :: crit_onset_swi  ! critical number of days > soilpsi_on for onset
     real(r8) :: soilpsi_on      ! critical soil water potential for leaf onset
     real(r8) :: crit_offset_fdd ! critical number of freezing days to initiate offset
     real(r8) :: crit_offset_swi ! critical number of water stress days to initiate offset
     real(r8) :: soilpsi_off     ! critical soil water potential for leaf offset
     real(r8) :: lwtop   	 ! live wood turnover proportion (annual fraction)
  end type params_type

  type(params_type) :: params_inst

  real(r8) :: dt                            ! radiation time step delta t (seconds)
  real(r8) :: fracday                       ! dtime as a fraction of day
  real(r8) :: crit_dayl                     ! critical daylength for offset (seconds)
  real(r8) :: ndays_on                      ! number of days to complete onset
  real(r8) :: ndays_off                     ! number of days to complete offset
  real(r8) :: fstor2tran                    ! fraction of storage to move to transfer on each onset
  real(r8) :: crit_onset_fdd                ! critical number of freezing days
  real(r8) :: crit_onset_swi                ! water stress days for offset trigger
  real(r8) :: soilpsi_on                    ! water potential for onset trigger (MPa)
  real(r8) :: crit_offset_fdd               ! critical number of freezing degree days to trigger offset
  real(r8) :: crit_offset_swi               ! water stress days for offset trigger
  real(r8) :: soilpsi_off                   ! water potential for offset trigger (MPa)
  real(r8) :: lwtop                         ! live wood turnover proportion (annual fraction)

  ! CropPhenology variables and constants
  real(r8) :: p1d, p1v                      ! photoperiod factor constants for crop vernalization
  real(r8) :: hti                           ! cold hardening index threshold for vernalization
  real(r8) :: tbase                         ! base temperature for vernalization

  integer, parameter :: NOT_Planted   = 999 ! If not planted   yet in year
  integer, parameter :: NOT_Harvested = 999 ! If not harvested yet in year
  integer, parameter :: inNH       = 1      ! Northern Hemisphere
  integer, parameter :: inSH       = 2      ! Southern Hemisphere
  integer, pointer   :: inhemi(:)           ! Hemisphere that patch is in 

  integer, allocatable :: minplantjday(:,:) ! minimum planting julian day
  integer, allocatable :: maxplantjday(:,:) ! maximum planting julian day
  integer              :: jdayyrstart(inSH) ! julian day of start of year

  real(r8), private :: initial_seed_at_planting = 3._r8 ! Initial seed at planting

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

    ! !ARGUMENTS:
    implicit none
    type(file_desc_t),intent(inout) :: ncid   ! pio netCDF file id
    !
    ! !LOCAL VARIABLES:
    character(len=32)  :: subname = 'CNPhenolParamsType'
    character(len=100) :: errCode = '-Error reading in parameters file:'
    logical            :: readv ! has variable been read in or not
    real(r8)           :: tempr ! temporary to read in parameter
    character(len=100) :: tString ! temp. var for reading
    !-----------------------------------------------------------------------

    !
    ! read in parameters
    !   
    tString='crit_dayl'
    call ncd_io(varname=trim(tString),data=tempr, flag='read', ncid=ncid, readvar=readv)
    if ( .not. readv ) call endrun( msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))
    params_inst%crit_dayl=tempr

    tString='ndays_on'
    call ncd_io(varname=trim(tString),data=tempr, flag='read', ncid=ncid, readvar=readv)
    if ( .not. readv ) call endrun( msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))
    params_inst%ndays_on=tempr

    tString='ndays_off'
    call ncd_io(varname=trim(tString),data=tempr, flag='read', ncid=ncid, readvar=readv)
    if ( .not. readv ) call endrun( msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))
    params_inst%ndays_off=tempr

    tString='fstor2tran'
    call ncd_io(varname=trim(tString),data=tempr, flag='read', ncid=ncid, readvar=readv)
    if ( .not. readv ) call endrun( msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))
    params_inst%fstor2tran=tempr

    tString='crit_onset_fdd'
    call ncd_io(varname=trim(tString),data=tempr, flag='read', ncid=ncid, readvar=readv)
    if ( .not. readv ) call endrun( msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))
    params_inst%crit_onset_fdd=tempr

    tString='crit_onset_swi'
    call ncd_io(varname=trim(tString),data=tempr, flag='read', ncid=ncid, readvar=readv)
    if ( .not. readv ) call endrun( msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))
    params_inst%crit_onset_swi=tempr

    tString='soilpsi_on'
    call ncd_io(varname=trim(tString),data=tempr, flag='read', ncid=ncid, readvar=readv)
    if ( .not. readv ) call endrun( msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))
    params_inst%soilpsi_on=tempr

    tString='crit_offset_fdd'
    call ncd_io(varname=trim(tString),data=tempr, flag='read', ncid=ncid, readvar=readv)
    if ( .not. readv ) call endrun( msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))
    params_inst%crit_offset_fdd=tempr

    tString='crit_offset_swi'
    call ncd_io(varname=trim(tString),data=tempr, flag='read', ncid=ncid, readvar=readv)
    if ( .not. readv ) call endrun( msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))
    params_inst%crit_offset_swi=tempr

    tString='soilpsi_off'
    call ncd_io(varname=trim(tString),data=tempr, flag='read', ncid=ncid, readvar=readv)
    if ( .not. readv ) call endrun( msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))
    params_inst%soilpsi_off=tempr

    tString='lwtop_ann'
    call ncd_io(varname=trim(tString),data=tempr, flag='read', ncid=ncid, readvar=readv)
    if ( .not. readv ) call endrun( msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))
    params_inst%lwtop=tempr   

  end subroutine readParams

  !-----------------------------------------------------------------------
  subroutine CNPhenologyInit(bounds)
    !
    ! !DESCRIPTION:
    ! Initialization of CNPhenology. Must be called after time-manager is
    ! initialized, and after pftcon file is read in.
    !
    ! !USES:
    use clm_time_manager, only: get_step_size
    use clm_varcon      , only: secspday
    !
    ! !ARGUMENTS:
    type(bounds_type), intent(in) :: bounds  
    !------------------------------------------------------------------------

    !
    ! Get time-step and what fraction of a day it is
    !
    dt      = real( get_step_size(), r8 )
    fracday = dt/secspday

    ! set constants for CNSeasonDecidPhenology 
    ! (critical daylength from Biome-BGC, v4.1.2)
    crit_dayl=params_inst%crit_dayl

    ! Set constants for CNSeasonDecidPhenology and CNStressDecidPhenology
    ndays_on=params_inst%ndays_on
    ndays_off=params_inst%ndays_off

    ! set transfer parameters
    fstor2tran=params_inst%fstor2tran

    ! -----------------------------------------
    ! Constants for CNStressDecidPhenology
    ! -----------------------------------------

    ! onset parameters
    crit_onset_fdd=params_inst%crit_onset_fdd
    ! critical onset gdd now being calculated as a function of annual
    ! average 2m temp.
    ! crit_onset_gdd = 150.0 ! c3 grass value
    ! crit_onset_gdd = 1000.0   ! c4 grass value
    crit_onset_swi=params_inst%crit_onset_swi
    soilpsi_on=params_inst%soilpsi_on

    ! offset parameters
    crit_offset_fdd=params_inst%crit_offset_fdd
    crit_offset_swi=params_inst%crit_offset_swi
    soilpsi_off=params_inst%soilpsi_off    

    ! -----------------------------------------
    ! Constants for CNLivewoodTurnover
    ! -----------------------------------------

    ! set the global parameter for livewood turnover rate
    ! define as an annual fraction (0.7), and convert to fraction per second
    lwtop=params_inst%lwtop/31536000.0_r8 !annual fraction converted to per second

  end subroutine CNPhenologyInit

end module CNPhenologyMod
