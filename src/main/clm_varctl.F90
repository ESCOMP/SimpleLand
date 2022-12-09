module clm_varctl

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Module containing run control variables
  !
  ! !USES:
  use shr_kind_mod, only: r8 => shr_kind_r8, SHR_KIND_CL
  use shr_sys_mod , only: shr_sys_abort ! cannot use endrun here due to circular dependency
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  implicit none
  public :: clm_varctl_set    ! Set variables
  !
  private
  save
  !
  ! !PUBLIC TYPES:
  !
  integer , parameter, public ::  iundef = -9999999
  real(r8), parameter, public ::  rundef = -9999999._r8
  integer , parameter, public ::  fname_len = SHR_KIND_CL   ! max length of file names in this module
  !----------------------------------------------------------
  !
  ! Run control variables
  !
  ! case id
  character(len=256), public :: caseid  = ' '                            

  ! case title
  character(len=256), public :: ctitle  = ' '                            

  ! Type of run
  integer, public :: nsrest                = iundef                         
  logical, public :: is_cold_start         = .false.
  logical, public :: is_interpolated_start = .false. ! True if we're starting from initial conditions that have been run through init_interp

  ! Startup from initial conditions
  integer, public, parameter :: nsrStartup  = 0                          

  ! Continue from restart files
  integer, public, parameter :: nsrContinue = 1                          

  ! Branch from restart files
  integer, public, parameter :: nsrBranch   = 2                          

  ! true => allow case name to remain the same for branch run
  ! by default this is not allowed
  logical, public :: brnch_retain_casename = .false.                     

  !true => no valid land points -- do NOT run
  logical, public :: noland = .false.                                    

  ! Hostname of machine running on
  character(len=256), public :: hostname = ' '                           

  ! username of user running program
  character(len=256), public :: username = ' '                           

  ! description of this source
  character(len=256), public :: source   = "Community Land Model CLM4.0" 

  ! version of program
  character(len=256), public :: version  = " "                           

  ! dataset conventions
  character(len=256), public :: conventions = "CF-1.0"                   

  !----------------------------------------------------------
  ! Unit Numbers
  !----------------------------------------------------------
  !
  integer, public :: iulog = 6        ! "stdout" log file unit number, default is 6

  !----------------------------------------------------------
  ! Output NetCDF files
  !----------------------------------------------------------

  logical, public :: outnc_large_files = .true.         ! large file support for output NetCDF files

  !----------------------------------------------------------
  ! Run input files
  !----------------------------------------------------------

  character(len=fname_len), public :: finidat    = ' '        ! initial conditions file name
  character(len=fname_len), public :: fsurdat    = ' '        ! surface data file name
  character(len=fname_len), public :: fatmgrid   = ' '        ! atm grid file name
  character(len=fname_len), public :: fatmlndfrc = ' '        ! lnd frac file on atm grid
  character(len=fname_len), public :: nrevsn     = ' '        ! restart data file name for branch run

  !----------------------------------------------------------
  ! MML input files
  !----------------------------------------------------------
  character(len=fname_len), public :: mml_surdat   = ' '      ! MML surface data file for simple model
  	

  !----------------------------------------------------------
  ! Interpolation of finidat if requested
  !----------------------------------------------------------

  ! If finidat_interp_source is non-blank and finidat is blank then interpolation will be
  ! done from finidat_interp_source to finidat_interp_dest. Note that
  ! finidat_interp_source is not read in directly from the namelist - rather, it is set
  ! from finidat if use_init_interp is .true.

  character(len=fname_len), public :: finidat_interp_source = ' '
  character(len=fname_len), public :: finidat_interp_dest   = 'finidat_interp_dest.nc'     

  !----------------------------------------------------------
  ! BGC logic and datasets
  !----------------------------------------------------------

  ! values of 'prognostic','diagnostic','constant'
  character(len=16), public :: co2_type = 'constant'    

  !----------------------------------------------------------
  ! Physics
  !----------------------------------------------------------

  ! true => write global average diagnostics to std out
  logical,  public :: wrtdia       = .false.            

  ! atmospheric CO2 molar ratio (by volume) (umol/mol)
  real(r8), public :: co2_ppmv     = 355._r8            !

  !----------------------------------------------------------
  ! single column control variables
  !----------------------------------------------------------

  logical,  public :: single_column = .false. ! true => single column mode
  real(r8), public :: scmlat        = rundef  ! single column lat
  real(r8), public :: scmlon        = rundef  ! single column lon

  !----------------------------------------------------------
  ! instance control
  !----------------------------------------------------------

  integer, public :: inst_index
  character(len=16), public :: inst_name
  character(len=16), public :: inst_suffix

  !----------------------------------------------------------
  ! Decomp control variables
  !----------------------------------------------------------

  ! number of segments per clump for decomp
  integer, public :: nsegspc = 20                       

  !----------------------------------------------------------
  ! Derived variables (run, history and restart file)
  !----------------------------------------------------------

  ! directory name for local restart pointer file
  character(len=256), public :: rpntdir = '.'            

  ! file name for local restart pointer file
  character(len=256), public :: rpntfil = 'rpointer.lnd' 

  !----------------------------------------------------------
  ! Migration of CPP variables
  !----------------------------------------------------------
  logical, public :: use_noio            = .false.

  !----------------------------------------------------------
  ! To retrieve namelist
  !----------------------------------------------------------
  character(len=SHR_KIND_CL), public :: NLFilename_in ! Namelist filename
  !
  logical, private :: clmvarctl_isset = .false.
 !-----------------------------------------------------------------------

contains

  !---------------------------------------------------------------------------
  subroutine clm_varctl_set( caseid_in, ctitle_in, brnch_retain_casename_in,    &
       single_column_in, scmlat_in, scmlon_in, nsrest_in, &
       version_in, hostname_in, username_in)
    !
    ! !DESCRIPTION:
    ! Set input control variables.
    !
    ! !ARGUMENTS:
    character(len=256), optional, intent(IN) :: caseid_in                ! case id
    character(len=256), optional, intent(IN) :: ctitle_in                ! case title
    logical,            optional, intent(IN) :: brnch_retain_casename_in ! true => allow case name to remain the 
                                                                         ! same for branch run
    logical,            optional, intent(IN) :: single_column_in         ! true => single column mode
    real(r8),           optional, intent(IN) :: scmlat_in                ! single column lat
    real(r8),           optional, intent(IN) :: scmlon_in                ! single column lon
    integer,            optional, intent(IN) :: nsrest_in                ! 0: initial run. 1: restart: 3: branch
    character(len=256), optional, intent(IN) :: version_in               ! model version
    character(len=256), optional, intent(IN) :: hostname_in              ! hostname running on
    character(len=256), optional, intent(IN) :: username_in              ! username running job
    !-----------------------------------------------------------------------

    if ( clmvarctl_isset )then
       call shr_sys_abort(' ERROR:: control variables already set, cannot call this routine')
    end if

    if ( present(caseid_in       ) ) caseid        = caseid_in
    if ( present(ctitle_in       ) ) ctitle        = ctitle_in
    if ( present(single_column_in) ) single_column = single_column_in
    if ( present(scmlat_in       ) ) scmlat        = scmlat_in
    if ( present(scmlon_in       ) ) scmlon        = scmlon_in
    if ( present(nsrest_in       ) ) nsrest        = nsrest_in
    if ( present(brnch_retain_casename_in) ) brnch_retain_casename = brnch_retain_casename_in
    if ( present(version_in      ) ) version       = version_in
    if ( present(username_in     ) ) username      = username_in
    if ( present(hostname_in     ) ) hostname      = hostname_in

  end subroutine clm_varctl_set

end module clm_varctl
