module controlMod

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Module which initializes run control variables. The following possible
  ! namelist variables are set default values and possibly read in on startup
  !
  ! Note: For definitions of namelist variables see
  !       ../../bld/namelist_files/namelist_definition.xml
  !       Display the file in a browser to see it neatly formatted in html.
  !
  ! !USES:
  use shr_kind_mod                     , only: r8 => shr_kind_r8, SHR_KIND_CL
  use shr_nl_mod                       , only: shr_nl_find_group_name
  use shr_const_mod                    , only: SHR_CONST_CDAY
  use shr_log_mod                      , only: errMsg => shr_log_errMsg
  use abortutils                       , only: endrun
  use spmdMod                          , only: masterproc
  use decompMod                        , only: clump_pproc
  use clm_varpar                       , only: numrad, nlevsno
  use histFileMod                      , only: max_tapes, max_namlen 
  use histFileMod                      , only: hist_empty_htapes, hist_dov2xy, hist_avgflag_pertape, hist_type1d_pertape 
  use histFileMod                      , only: hist_nhtfrq, hist_ndens, hist_mfilt, hist_fincl1, hist_fincl2, hist_fincl3
  use histFileMod                      , only: hist_fincl4, hist_fincl5, hist_fincl6, hist_fexcl1, hist_fexcl2, hist_fexcl3
  use histFileMod                      , only: hist_fexcl4, hist_fexcl5, hist_fexcl6
  use initInterpMod                    , only: initInterp_readnl
  use clm_varctl                       , only: iundef, rundef, nsrest, caseid, ctitle, nsrStartup, nsrContinue
  use clm_varctl                       , only: nsrBranch, brnch_retain_casename, hostname, username, source, version, conventions
  use clm_varctl                       , only: iulog, outnc_large_files, finidat, fsurdat, fatmgrid, fatmlndfrc, nrevsn
  use clm_varctl                       , only: mml_surdat, finidat_interp_source, finidat_interp_dest, co2_type
  use clm_varctl                       , only: wrtdia, co2_ppmv, nsegspc, rpntdir, rpntfil
  use clm_varctl                       , only: use_noio, NLFilename_in
  use clm_varctl                       , only: clm_varctl_set
  use clm_varctl                       , only: single_column
  !
  ! !PUBLIC TYPES:
  implicit none
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: control_setNL ! Set namelist filename
  public :: control_init  ! initial run control information
  public :: control_print ! print run control information
  !
  !
  ! !PRIVATE MEMBER FUNCTIONS:
  private :: apply_use_init_interp  ! apply the use_init_interp namelist option, if set
  !
  ! !PRIVATE TYPES:
  character(len=  7) :: runtyp(4)                        ! run type
  character(len=SHR_KIND_CL) :: NLFilename = 'lnd.stdin' ! Namelist filename

#if (defined _OPENMP)
  integer, external :: omp_get_max_threads  ! max number of threads that can execute concurrently in a single parallel region
#endif

  character(len=*), parameter, private :: sourcefile = &
       __FILE__
  !-----------------------------------------------------------------------

contains

  !------------------------------------------------------------------------
  subroutine control_setNL( NLfile )
    !
    ! !DESCRIPTION:
    ! Set the namelist filename to use
    !
    ! !ARGUMENTS:
    character(len=*), intent(IN) :: NLFile ! Namelist filename
    !
    ! !LOCAL VARIABLES:
    character(len=32) :: subname = 'control_setNL'  ! subroutine name
    logical :: lexist                               ! File exists
    !------------------------------------------------------------------------

    ! Error checking...
    if ( len_trim(NLFile) == 0 )then
       call endrun(msg=' error: nlfilename entered is not set'//errMsg(sourcefile, __LINE__))
    end if
    inquire (file = trim(NLFile), exist = lexist)
    if ( .not. lexist )then
       call endrun(msg=' error: NLfilename entered does NOT exist:'//&
            trim(NLFile)//errMsg(sourcefile, __LINE__))
    end if
    if ( len_trim(NLFile) > len(NLFilename) )then
       call endrun(msg=' error: entered NLFile is too long'//errMsg(sourcefile, __LINE__))
    end if
    ! Set the filename
    NLFilename = NLFile
    NLFilename_in = NLFilename   ! For use in external namelists and to avoid creating dependencies on controlMod
  end subroutine control_setNL

  !------------------------------------------------------------------------
  subroutine control_init( )
    !
    ! !DESCRIPTION:
    ! Initialize CLM run control information
    !
    ! !USES:
    use clm_time_manager                 , only : set_timemgr_init
    use fileutils                        , only : getavu, relavu
    !
    ! !LOCAL VARIABLES:
    integer :: i                    ! loop indices
    integer :: ierr                 ! error code
    integer :: unitn                ! unit for namelist file
    integer :: dtime                ! Integer time-step
    integer :: override_nsrest      ! If want to override the startup type sent from driver
    logical :: use_init_interp      ! Apply initInterp to the file given by finidat
    !------------------------------------------------------------------------

    ! ----------------------------------------------------------------------
    ! Namelist Variables
    ! ----------------------------------------------------------------------

    ! Time step
    namelist / clm_inparm/ &
    dtime

    ! CLM namelist settings

    namelist /clm_inparm / &
         fatmlndfrc, finidat, nrevsn, &
         finidat_interp_dest, &
         use_init_interp

    ! Input datasets

    namelist /clm_inparm/ fsurdat

	! MML Input datasets for simple model
    namelist /clm_inparm/ &
    	 mml_surdat			
    	 ! MML forcing file w/ albedo, roughness, etc
    	 ! /glade/p/work/mlague/cesm_source/cesm1_5_beta05_mml_land/components/clm/bld/namelist_files/namelist_defaults.xml
    	 ! I think I need to modify one of the namelis_defaults xml files in the above folder in order for 
    	 ! the model to know to accept my new namelist var... 

    ! History, restart options

    namelist /clm_inparm/  &
         hist_empty_htapes, hist_dov2xy, &
         hist_avgflag_pertape, hist_type1d_pertape, &
         hist_nhtfrq,  hist_ndens, hist_mfilt, &
         hist_fincl1,  hist_fincl2, hist_fincl3, &
         hist_fincl4,  hist_fincl5, hist_fincl6, &
         hist_fexcl1,  hist_fexcl2, hist_fexcl3, &
         hist_fexcl4,  hist_fexcl5, hist_fexcl6

    ! BGC info

    namelist /clm_inparm / &
         co2_type

    ! Glacier_mec info
    namelist /clm_inparm/ nlevsno

    ! Other options

    namelist /clm_inparm/  &
         clump_pproc, wrtdia, &
         nsegspc, co2_ppmv, override_nsrest

    ! All old cpp-ifdefs are below and have been converted to namelist variables 
    namelist /clm_inparm/ use_noio

    ! Items not really needed, but do need to be properly set as they are used
    namelist / clm_inparm/ single_column

    ! ----------------------------------------------------------------------
    ! Default values
    ! ----------------------------------------------------------------------

    if (masterproc) then
       write(iulog,*) 'Attempting to initialize run control settings .....'
    endif

    runtyp(:)               = 'missing'
    runtyp(nsrStartup  + 1) = 'initial'
    runtyp(nsrContinue + 1) = 'restart'
    runtyp(nsrBranch   + 1) = 'branch '

    ! Set clumps per procoessor

#if (defined _OPENMP)
    clump_pproc = omp_get_max_threads()
    if ( clump_pproc > 1 )then
       call endrun(msg=' error: number of OMP threads is NOT 1 -- SLIM can NOT run with more than one thread'//errMsg(sourcefile, __LINE__))
    end if
#else
    clump_pproc = 1
#endif

    override_nsrest = nsrest

    use_init_interp = .false.

    if (masterproc) then

       ! ----------------------------------------------------------------------
       ! Read namelist from standard input. 
       ! ----------------------------------------------------------------------

       if ( len_trim(NLFilename) == 0  )then
          call endrun(msg=' error: nlfilename not set'//errMsg(sourcefile, __LINE__))
       end if
       unitn = getavu()
       write(iulog,*) 'Read in clm_inparm namelist from: ', trim(NLFilename)
       open( unitn, file=trim(NLFilename), status='old' )
       call shr_nl_find_group_name(unitn, 'clm_inparm', status=ierr)
       if (ierr == 0) then
          read(unitn, clm_inparm, iostat=ierr)
          if (ierr /= 0) then
             call endrun(msg='ERROR reading clm_inparm namelist'//errMsg(sourcefile, __LINE__))
          end if
       else
          call endrun(msg='ERROR finding clm_inparm namelist'//errMsg(sourcefile, __LINE__))
       end if

       call relavu( unitn )

       ! ----------------------------------------------------------------------
       ! Process some namelist variables, and perform consistency checks
       ! ----------------------------------------------------------------------

       call set_timemgr_init( dtime_in=dtime )

       ! Check for namelist variables that SLIM can NOT use
       if ( single_column )then
          call endrun(msg='ERROR SLIM can NOT run with single_column on'//errMsg(sourcefile, __LINE__))
       end if

       if (use_init_interp) then
          call apply_use_init_interp(finidat, finidat_interp_source)
       end if

       ! History and restart files

       do i = 1, max_tapes
          if (hist_nhtfrq(i) == 0) then
             hist_mfilt(i) = 1
          else if (hist_nhtfrq(i) < 0) then
             hist_nhtfrq(i) = nint(-hist_nhtfrq(i)*SHR_CONST_CDAY/(24._r8*dtime))
          endif
       end do

       ! Override start-type (can only override to branch (3)  and only 
       ! if the driver is a startup type
       if ( override_nsrest /= nsrest )then
           if ( override_nsrest /= nsrBranch .and. nsrest /= nsrStartup )then
              call endrun(msg= ' ERROR: can ONLY override clm start-type ' // &
                   'to branch type and ONLY if driver is a startup type'// &
                   errMsg(sourcefile, __LINE__))
           end if
           call clm_varctl_set( nsrest_in=override_nsrest )
       end if

       ! If nlevsno are equal to their junk
       ! default value, then they were not specified by the user namelist and we generate
       ! an error message. Also check nlevsno for bounds.
       if (nlevsno < 3 .or. nlevsno > 12)  then
          write(iulog,*)'ERROR: nlevsno = ',nlevsno,' is not supported, must be in range 3-12.'
          call endrun(msg=' ERROR: invalid value for nlevsno in CLM namelist. '//&
               errMsg(sourcefile, __LINE__))
       endif
    endif   ! end of if-masterproc if-block

    ! ----------------------------------------------------------------------
    ! Read in other namelists for other modules
    ! ----------------------------------------------------------------------

    call initInterp_readnl( NLFilename )

    ! ----------------------------------------------------------------------
    ! Broadcast all control information if appropriate
    ! ----------------------------------------------------------------------

    call control_spmd()
    
    ! ----------------------------------------------------------------------
    ! consistency checks
    ! ----------------------------------------------------------------------

    ! Consistency settings for co2 type
    if (co2_type /= 'constant' .and. co2_type /= 'prognostic' .and. co2_type /= 'diagnostic') then
       write(iulog,*)'co2_type = ',co2_type,' is not supported'
       call endrun(msg=' ERROR:: choices are constant, prognostic or diagnostic'//&
            errMsg(sourcefile, __LINE__))
    end if

    ! Check on run type
    if (nsrest == iundef) then
       call endrun(msg=' ERROR:: must set nsrest'//& 
            errMsg(sourcefile, __LINE__))
    end if
    if (nsrest == nsrBranch .and. nrevsn == ' ') then
       call endrun(msg=' ERROR: need to set restart data file name'//&
            errMsg(sourcefile, __LINE__))
    end if

    ! Consistency settings for co2_ppvm
    if ( (co2_ppmv <= 0.0_r8) .or. (co2_ppmv > 3000.0_r8) ) then
       call endrun(msg=' ERROR: co2_ppmv is out of a reasonable range'//& 
            errMsg(sourcefile, __LINE__))
    end if

    ! Consistency settings for nrevsn

    if (nsrest == nsrStartup ) nrevsn = ' '
    if (nsrest == nsrContinue) nrevsn = 'set by restart pointer file file'
    if (nsrest /= nsrStartup .and. nsrest /= nsrContinue .and. nsrest /= nsrBranch ) then
       call endrun(msg=' ERROR: nsrest NOT set to a valid value'//&
            errMsg(sourcefile, __LINE__))
    end if

    if (masterproc) then
       write(iulog,*) 'Successfully initialized run control settings'
       write(iulog,*)
    endif

  end subroutine control_init

  !------------------------------------------------------------------------
  subroutine control_spmd()
    !
    ! !DESCRIPTION:
    ! Distribute namelist data all processors. All program i/o is 
    ! funnelled through the master processor. Processor 0 either 
    ! reads restart/history data from the disk and distributes 
    ! it to all processors, or collects data from
    ! all processors and writes it to disk.
    !
    ! !USES:
    use spmdMod,    only : mpicom, MPI_CHARACTER, MPI_INTEGER, MPI_LOGICAL, MPI_REAL8
    !
    ! !ARGUMENTS:
    !
    ! !LOCAL VARIABLES:
    integer ier       !error code
    !-----------------------------------------------------------------------

    ! run control variables
    call mpi_bcast (caseid, len(caseid), MPI_CHARACTER, 0, mpicom, ier)
    call mpi_bcast (ctitle, len(ctitle), MPI_CHARACTER, 0, mpicom, ier)
    call mpi_bcast (version, len(version), MPI_CHARACTER, 0, mpicom, ier)
    call mpi_bcast (hostname, len(hostname), MPI_CHARACTER, 0, mpicom, ier)
    call mpi_bcast (username, len(username), MPI_CHARACTER, 0, mpicom, ier)
    call mpi_bcast (nsrest, 1, MPI_INTEGER, 0, mpicom, ier)

    call mpi_bcast (use_noio, 1, MPI_LOGICAL, 0, mpicom, ier)

    ! initial file variables
    call mpi_bcast (nrevsn, len(nrevsn), MPI_CHARACTER, 0, mpicom, ier)
    call mpi_bcast (finidat, len(finidat), MPI_CHARACTER, 0, mpicom, ier)
    call mpi_bcast (finidat_interp_source, len(finidat_interp_source), MPI_CHARACTER, 0, mpicom, ier)
    call mpi_bcast (finidat_interp_dest, len(finidat_interp_dest), MPI_CHARACTER, 0, mpicom, ier)
    call mpi_bcast (fsurdat, len(fsurdat), MPI_CHARACTER, 0, mpicom, ier)
    call mpi_bcast (fatmlndfrc,len(fatmlndfrc),MPI_CHARACTER, 0, mpicom, ier)

	! mml input file vars for simple model
	call mpi_bcast (mml_surdat,  len(mml_surdat),   MPI_CHARACTER, 0, mpicom, ier)
	
    call mpi_bcast (co2_type, len(co2_type), MPI_CHARACTER, 0, mpicom, ier)

    ! physics variables
    call mpi_bcast (nsegspc, 1, MPI_INTEGER, 0, mpicom, ier)
    call mpi_bcast (wrtdia, 1, MPI_LOGICAL, 0, mpicom, ier)
    call mpi_bcast (single_column,1, MPI_LOGICAL, 0, mpicom, ier)
    call mpi_bcast (co2_ppmv, 1, MPI_REAL8,0, mpicom, ier)

    ! snow pack variables
    call mpi_bcast (nlevsno, 1, MPI_INTEGER, 0, mpicom, ier)

    ! history file variables
    call mpi_bcast (hist_empty_htapes, 1, MPI_LOGICAL, 0, mpicom, ier)
    call mpi_bcast (hist_dov2xy, size(hist_dov2xy), MPI_LOGICAL, 0, mpicom, ier)
    call mpi_bcast (hist_nhtfrq, size(hist_nhtfrq), MPI_INTEGER, 0, mpicom, ier)
    call mpi_bcast (hist_mfilt, size(hist_mfilt), MPI_INTEGER, 0, mpicom, ier)
    call mpi_bcast (hist_ndens, size(hist_ndens), MPI_INTEGER, 0, mpicom, ier)
    call mpi_bcast (hist_avgflag_pertape, size(hist_avgflag_pertape), MPI_CHARACTER, 0, mpicom, ier)
    call mpi_bcast (hist_type1d_pertape, max_namlen*size(hist_type1d_pertape), MPI_CHARACTER, 0, mpicom, ier)
    call mpi_bcast (hist_fexcl1, max_namlen*size(hist_fexcl1), MPI_CHARACTER, 0, mpicom, ier)
    call mpi_bcast (hist_fexcl2, max_namlen*size(hist_fexcl2), MPI_CHARACTER, 0, mpicom, ier)
    call mpi_bcast (hist_fexcl3, max_namlen*size(hist_fexcl3), MPI_CHARACTER, 0, mpicom, ier)
    call mpi_bcast (hist_fexcl4, max_namlen*size(hist_fexcl4), MPI_CHARACTER, 0, mpicom, ier)
    call mpi_bcast (hist_fexcl5, max_namlen*size(hist_fexcl5), MPI_CHARACTER, 0, mpicom, ier)
    call mpi_bcast (hist_fexcl6, max_namlen*size(hist_fexcl6), MPI_CHARACTER, 0, mpicom, ier)
    call mpi_bcast (hist_fincl1, (max_namlen+2)*size(hist_fincl1), MPI_CHARACTER, 0, mpicom, ier)
    call mpi_bcast (hist_fincl2, (max_namlen+2)*size(hist_fincl2), MPI_CHARACTER, 0, mpicom, ier)
    call mpi_bcast (hist_fincl3, (max_namlen+2)*size(hist_fincl3), MPI_CHARACTER, 0, mpicom, ier)
    call mpi_bcast (hist_fincl4, (max_namlen+2)*size(hist_fincl4), MPI_CHARACTER, 0, mpicom, ier)
    call mpi_bcast (hist_fincl5, (max_namlen+2)*size(hist_fincl5), MPI_CHARACTER, 0, mpicom, ier)
    call mpi_bcast (hist_fincl6, (max_namlen+2)*size(hist_fincl6), MPI_CHARACTER, 0, mpicom, ier)

    ! restart file variables

    call mpi_bcast (rpntfil, len(rpntfil), MPI_CHARACTER, 0, mpicom, ier)

    ! clump decomposition variables

    call mpi_bcast (clump_pproc, 1, MPI_INTEGER, 0, mpicom, ier)

  end subroutine control_spmd

  !------------------------------------------------------------------------
  subroutine control_print ()
    !
    ! !DESCRIPTION:
    ! Write out the clm namelist run control variables
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    !
    ! !LOCAL VARIABLES:
    integer i  !loop index
    !------------------------------------------------------------------------

    write(iulog,*) 'define run:'
    write(iulog,*) '   source                = ',trim(source)
    write(iulog,*) '   model_version         = ',trim(version)
    write(iulog,*) '   run type              = ',runtyp(nsrest+1)
    write(iulog,*) '   case title            = ',trim(ctitle)
    write(iulog,*) '   username              = ',trim(username)
    write(iulog,*) '   hostname              = ',trim(hostname)
    write(iulog,*) 'process control parameters:'
    write(iulog,*) '    use_noio = ', use_noio

    write(iulog,*) 'input data files:'
    if (fsurdat == ' ') then
       write(iulog,*) '   fsurdat, surface dataset not set'
    else
       write(iulog,*) '   surface data   = ',trim(fsurdat)
    end if
    if (fatmlndfrc == ' ') then
       write(iulog,*) '   fatmlndfrc not set, setting frac/mask to 1'
    else
       write(iulog,*) '   land frac data = ',trim(fatmlndfrc)
    end if
    if (mml_surdat == ' ') then
       write(iulog,*) '   mml_surdat NOT set, check that we are using the default'
    else
       write(iulog,*) '   mml_surdat IS set, and = ',trim(mml_surdat)
    end if
    write(iulog,*) '   Number of snow layers =', nlevsno

    if (nsrest == nsrStartup) then
       if (finidat /= ' ') then
          write(iulog,*) '   initial data: ', trim(finidat)
       else if (finidat_interp_source /= ' ') then
          write(iulog,*) '   initial data interpolated from: ', trim(finidat_interp_source)
       else
          write(iulog,*) '   initial data created by model (cold start)'
       end if
    else
       write(iulog,*) '   restart data   = ',trim(nrevsn)
    end if

    write(iulog,*) '   atmospheric forcing data is from cesm atm model'
    write(iulog,*) 'Restart parameters:'
    write(iulog,*)'   restart pointer file directory     = ',trim(rpntdir)
    write(iulog,*)'   restart pointer file name          = ',trim(rpntfil)
    write(iulog,*) 'model physics parameters:'

    if ( trim(co2_type) == 'constant' )then
       write(iulog,*) '   CO2 volume mixing ratio   (umol/mol)   = ', co2_ppmv
    else
       write(iulog,*) '   CO2 volume mixing ratio                = ', co2_type
    end if

    if (nsrest == nsrContinue) then
       write(iulog,*) 'restart warning:'
       write(iulog,*) '   Namelist not checked for agreement with initial run.'
       write(iulog,*) '   Namelist should not differ except for ending time step and run type'
    end if
    if (nsrest == nsrBranch) then
       write(iulog,*) 'branch warning:'
       write(iulog,*) '   Namelist not checked for agreement with initial run.'
       write(iulog,*) '   Surface data set and reference date should not differ from initial run'
    end if
    write(iulog,*) '   nsegspc              = ',nsegspc

  end subroutine control_print


  !-----------------------------------------------------------------------
  subroutine apply_use_init_interp(finidat, finidat_interp_source)
    !
    ! !DESCRIPTION:
    ! Applies the use_init_interp option, setting finidat_interp_source to finidat
    !
    ! Should be called if use_init_interp is true.
    !
    ! Does error checking to ensure that it is valid to set use_init_interp to true,
    ! given the values of finidat and finidat_interp_source.
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    character(len=*), intent(inout) :: finidat
    character(len=*), intent(inout) :: finidat_interp_source
    !
    ! !LOCAL VARIABLES:

    character(len=*), parameter :: subname = 'apply_use_init_interp'
    !-----------------------------------------------------------------------

    if (finidat == ' ') then
       call endrun(msg=' ERROR: Can only set use_init_interp if finidat is set')
    end if

    if (finidat_interp_source /= ' ') then
       call endrun(msg=' ERROR: Cannot set use_init_interp if finidat_interp_source is &
            &already set')
    end if

    finidat_interp_source = finidat
    finidat = ' '

  end subroutine apply_use_init_interp



end module controlMod
