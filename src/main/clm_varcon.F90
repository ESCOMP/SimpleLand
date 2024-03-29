module clm_varcon

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Module containing various model constants.
  !
  ! !USES:
  use shr_kind_mod , only: r8 => shr_kind_r8
  use shr_const_mod, only: SHR_CONST_G,SHR_CONST_STEBOL,SHR_CONST_KARMAN,     &
                           SHR_CONST_RWV,SHR_CONST_RDAIR,SHR_CONST_CPFW,      &
                           SHR_CONST_CPICE,SHR_CONST_CPDAIR,SHR_CONST_LATVAP, &
                           SHR_CONST_LATSUB,SHR_CONST_LATICE,SHR_CONST_RHOFW, &
                           SHR_CONST_RHOICE,SHR_CONST_TKFRZ,SHR_CONST_REARTH, &
                           SHR_CONST_PDB, SHR_CONST_PI, SHR_CONST_CDAY,       &
                           SHR_CONST_RGAS, SHR_CONST_PSTD,                    &
                           SHR_CONST_MWDAIR, SHR_CONST_MWWV
  !
  ! !PUBLIC TYPES:
  implicit none
  save
  !-----------------------------------------------------------------------
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: clm_varcon_init  ! initialize constants in clm_varcon
  !
  ! !REVISION HISTORY:
  ! Created by Mariana Vertenstein
  ! 27 February 2008: Keith Oleson; Add forcing height and aerodynamic parameters
  !-----------------------------------------------------------------------

  !------------------------------------------------------------------
  ! Initialize mathmatical constants
  !------------------------------------------------------------------

  real(r8) :: rpi    = SHR_CONST_PI

  !------------------------------------------------------------------
  ! Initialize physical constants
  !------------------------------------------------------------------

  real(r8), parameter :: pc = 0.4                           ! threshold probability
  real(r8), parameter :: secsphr = 3600._r8                 ! Seconds in an hour
  integer,  parameter :: isecsphr = int(secsphr)            ! Integer seconds in an hour
  integer,  parameter :: isecspmin= 60                      ! Integer seconds in a minute
  real(r8) :: grav   = SHR_CONST_G                          ! gravity constant [m/s2]
  real(r8) :: sb     = SHR_CONST_STEBOL                     ! stefan-boltzmann constant  [W/m2/K4]
  real(r8) :: vkc    = SHR_CONST_KARMAN                     ! von Karman constant [-]
  real(r8) :: rwat   = SHR_CONST_RWV                        ! gas constant for water vapor [J/(kg K)]
  real(r8) :: rair   = SHR_CONST_RDAIR                      ! gas constant for dry air [J/kg/K]
  real(r8) :: roverg = SHR_CONST_RWV/SHR_CONST_G*1000._r8   ! Rw/g constant = (8.3144/0.018)/(9.80616)*1000. mm/K
  real(r8) :: cpliq  = SHR_CONST_CPFW                       ! Specific heat of water [J/kg-K]
  real(r8) :: cpice  = SHR_CONST_CPICE                      ! Specific heat of ice [J/kg-K]
  real(r8) :: cpair  = SHR_CONST_CPDAIR                     ! specific heat of dry air [J/kg/K]
  real(r8) :: hvap   = SHR_CONST_LATVAP                     ! Latent heat of evap for water [J/kg]
  real(r8) :: hsub   = SHR_CONST_LATSUB                     ! Latent heat of sublimation    [J/kg]
  real(r8) :: hfus   = SHR_CONST_LATICE                     ! Latent heat of fusion for ice [J/kg]
  real(r8) :: denh2o = SHR_CONST_RHOFW                      ! density of liquid water [kg/m3]
  real(r8) :: denice = SHR_CONST_RHOICE                     ! density of ice [kg/m3]
  real(r8) :: rgas   = SHR_CONST_RGAS                       ! universal gas constant [J/K/kmole]
  real(r8) :: pstd   = SHR_CONST_PSTD                       ! standard pressure [Pa]

  ! TODO(wjs, 2016-04-08) The following should be used in place of hard-coded constants
  ! of 0.622 and 0.378 (which is 1 - 0.622) in various places in the code:
  real(r8), parameter :: wv_to_dair_weight_ratio = SHR_CONST_MWWV/SHR_CONST_MWDAIR ! ratio of molecular weight of water vapor to that of dry air [-]

  real(r8) :: tkair  = 0.023_r8                             ! thermal conductivity of air   [W/m/K]
  real(r8) :: tkice  = 2.290_r8                             ! thermal conductivity of ice   [W/m/K]
  real(r8) :: tkwat  = 0.57_r8                              ! thermal conductivity of water [W/m/K]
  real(r8), parameter :: tfrz   = SHR_CONST_TKFRZ           ! freezing temperature [K]
  real(r8), parameter :: tcrit  = 2.5_r8                    ! critical temperature to determine rain or snow
  real(r8) :: o2_molar_const = 0.209_r8                     ! constant atmospheric O2 molar ratio (mol/mol)
  real(r8) :: oneatm = 1.01325e5_r8                         ! one standard atmospheric pressure [Pa]
  real(r8) :: bdsno = 250._r8                               ! bulk density snow (kg/m**3)
  real(r8) :: alpha_aero = 1.0_r8                           ! constant for aerodynamic parameter weighting
  real(r8) :: tlsai_crit = 2.0_r8                           ! critical value of elai+esai for which aerodynamic parameters are maximum
  real(r8) :: watmin = 0.01_r8                              ! minimum soil moisture (mm)

  real(r8) :: re = SHR_CONST_REARTH*0.001_r8                ! radius of earth (km)

  real(r8), public, parameter :: degpsec = 15._r8/3600.0_r8 ! Degree's earth rotates per second
  real(r8), public, parameter :: secspday= SHR_CONST_CDAY   ! Seconds per day
  integer,  public, parameter :: isecspday= secspday        ! Integer seconds per day

  ! ------------------------------------------------------------------------
  ! Special value flags
  ! ------------------------------------------------------------------------

  ! NOTE(wjs, 2015-11-23) The presence / absence of spval should be static in time for
  ! multi-level fields.  i.e., if a given level & column has spval at initialization, it
  ! should remain spval throughout the run (e.g., indicating that this level is not valid
  ! for this column type); similarly, if it starts as a valid value, it should never
  ! become spval. This is needed for init_interp to work correctly on multi-level fields.
  ! For more details, see the note near the top of initInterpMultilevelInterp.
  real(r8), public, parameter ::  spval = 1.e36_r8          ! special value for real data

  ! Keep this negative to avoid conflicts with possible valid values
  integer , public, parameter :: ispval = -9999             ! special value for int data

  integer, private :: i  ! loop index

  !------------------------------------------------------------------
  ! Set subgrid names
  !------------------------------------------------------------------

  character(len=16), parameter :: grlnd  = 'lndgrid'      ! name of lndgrid
  character(len=16), parameter :: nameg  = 'gridcell'     ! name of gridcells

  !------------------------------------------------------------------
  ! Soil depths are constants for now; lake depths can vary by gridcell
  ! zlak and dzlak correspond to the default 50 m lake depth.
  ! The values for the following arrays are set in routine iniTimeConst
  !------------------------------------------------------------------

  real(r8), allocatable :: zsoi(:)         !soil z  (layers)

contains

  !------------------------------------------------------------------------------
  subroutine clm_varcon_init()
    !
    ! !DESCRIPTION:
    ! This subroutine initializes constant arrays in clm_varcon. 
    ! MUST be called  after clm_varpar_init.
    !
    ! !USES:
    use clm_varpar, only: nlevgrnd
    !
    ! !ARGUMENTS:
    implicit none
    !
    ! !REVISION HISTORY:
    !   Created by E. Kluzek
!------------------------------------------------------------------------------

    allocate( zsoi(1:nlevgrnd                ))

  end subroutine clm_varcon_init

end module clm_varcon
