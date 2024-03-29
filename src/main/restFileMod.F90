module restFileMod

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Reads from or writes to/ the CLM restart file.
  !
  ! !USES:
#include "shr_assert.h"
  use shr_kind_mod     , only : r8 => shr_kind_r8
  use decompMod        , only : bounds_type, get_proc_clumps, get_clump_bounds
  use decompMod        , only : BOUNDS_LEVEL_PROC
  use spmdMod          , only : masterproc, mpicom
  use abortutils       , only : endrun
  use shr_log_mod      , only : errMsg => shr_log_errMsg
  use clm_time_manager , only : timemgr_restart_io, get_nstep
  use subgridRestMod   , only : subgridRestWrite
  use accumulMod       , only : accumulRest
  use clm_instMod      , only : clm_instRest
  use histFileMod      , only : hist_restart_ncd
  use clm_varctl       , only : iulog
  use clm_varcon       , only : nameg
  use ncdio_pio        , only : file_desc_t, ncd_pio_createfile, ncd_pio_openfile, ncd_global
  use ncdio_pio        , only : ncd_pio_closefile, ncd_defdim, ncd_putatt, ncd_enddef, check_dim
  use ncdio_pio        , only : check_att, ncd_getatt
  !
  ! !PUBLIC TYPES:
  implicit none
  private
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: restFile_read
  public :: restFile_write
  public :: restFile_open
  public :: restFile_close
  public :: restFile_getfile
  public :: restFile_filename        ! Sets restart filename
  !
  ! !PRIVATE MEMBER FUNCTIONS:
  private :: restFile_read_pfile
  private :: restFile_write_pfile       ! Writes restart pointer file
  private :: restFile_closeRestart      ! Close restart file and write restart pointer file
  private :: restFile_dimset
  private :: restFile_add_flag_metadata ! Add global metadata for some logical flag
  private :: restFile_dimcheck
  private :: restFile_enddef
  private :: restFile_check_consistency   ! Perform consistency checks on the restart file
  private :: restFile_read_consistency_nl ! Read namelist associated with consistency checks
  private :: restFile_check_year          ! Check consistency of year on the restart file
  !
  ! !PRIVATE TYPES:

  character(len=*), parameter, private :: sourcefile = &
       __FILE__
  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine restFile_write( bounds, file, rdate, noptr)
    !
    ! !DESCRIPTION:
    ! Define/write CLM restart file.
    !
    ! !ARGUMENTS:
    type(bounds_type) , intent(in)           :: bounds          
    character(len=*)  , intent(in)           :: file  ! output netcdf restart file
    character(len=*)  , intent(in), optional :: rdate ! restart file time stamp for name
    logical           , intent(in), optional :: noptr ! if should NOT write to the restart pointer file
    !
    ! !LOCAL VARIABLES:
    type(file_desc_t) :: ncid ! netcdf id
    integer :: i       ! index
    logical :: ptrfile ! write out the restart pointer file
    !-----------------------------------------------------------------------

    if ( present(noptr) )then
       ptrfile = .not. noptr
    else
       ptrfile = .true.
    end if

    ! Open file

    call restFile_open( flag='write', file=file, ncid=ncid )

    ! Define dimensions and variables

    call restFile_dimset ( ncid )

    call timemgr_restart_io(ncid, flag='define')

    call subgridRestWrite(bounds, ncid, flag='define' )

    call accumulRest( ncid, flag='define' )

    call clm_instRest(bounds, ncid, flag='define')

    if (present(rdate)) then 
       call hist_restart_ncd (bounds, ncid, flag='define', rdate=rdate )
    end if

    call restFile_enddef( ncid )

    ! Write variables
    
    call timemgr_restart_io( ncid, flag='write' )

    call subgridRestWrite(bounds, ncid, flag='write' )

    call accumulRest( ncid, flag='write' )

    call clm_instRest(bounds, ncid, flag='write')

    call hist_restart_ncd (bounds, ncid, flag='write' )

    ! Close file 
    
    call restFile_close( ncid )
    call restFile_closeRestart( file )
    
    ! Write restart pointer file
    
    if ( ptrfile ) call restFile_write_pfile( file )
    
    ! Write out diagnostic info

    if (masterproc) then
       write(iulog,*) 'Successfully wrote out restart data at nstep = ',get_nstep()
       write(iulog,'(72a1)') ("-",i=1,60)
    end if
    
  end subroutine restFile_write

  !-----------------------------------------------------------------------
  subroutine restFile_read( bounds_proc, file)
    !
    ! !DESCRIPTION:
    ! Read a CLM restart file.
    !
    ! !ARGUMENTS:
    type(bounds_type) , intent(in) :: bounds_proc      ! processor-level bounds
    character(len=*)  , intent(in) :: file             ! output netcdf restart file
    !
    ! !LOCAL VARIABLES:
    type(file_desc_t) :: ncid    ! netcdf id
    integer           :: i       ! index
    integer           :: nclumps ! number of clumps on this processor
    integer           :: nc      ! clump index
    type(bounds_type) :: bounds_clump

    character(len=*), parameter :: subname = 'restFile_read'
    !-----------------------------------------------------------------------

    ! The provided bounds need to be proc-level bounds. This is in part because of logic
    ! below that divides this into clump-level bounds for the sake of reweight_wrapup.
    ! But it *MAY* also be necessary to have proc-level bounds for these i/o routines.
    SHR_ASSERT(bounds_proc%level == BOUNDS_LEVEL_PROC, subname // ': argument must be PROC-level bounds')

    ! Open file

    call restFile_open( flag='read', file=file, ncid=ncid )

    ! Read file

    call restFile_dimcheck( ncid )

    call accumulRest( ncid, flag='read' )

    call clm_instRest( bounds_proc, ncid, flag='read' )

    call hist_restart_ncd (bounds_proc, ncid, flag='read' )

    ! Do error checking on file
    
    call restFile_check_consistency(bounds_proc, ncid)

    ! Close file 

    call restFile_close( ncid )

    ! Write out diagnostic info

    if (masterproc) then
       write(iulog,'(72a1)') ("-",i=1,60)
       write(iulog,*) 'Successfully read restart data for restart run'
       write(iulog,*)
    end if

  end subroutine restFile_read

  !-----------------------------------------------------------------------
  subroutine restFile_getfile( file, path )
    !
    ! !DESCRIPTION:
    ! Determine and obtain netcdf restart file
    !
    ! !USES:
    use clm_varctl, only : caseid, nrevsn, nsrest, brnch_retain_casename
    use clm_varctl, only : nsrContinue, nsrBranch
    use fileutils , only : getfil
    !
    ! !ARGUMENTS:
    character(len=*), intent(out) :: file  ! name of netcdf restart file
    character(len=*), intent(out) :: path  ! full pathname of netcdf restart file
    !
    ! !LOCAL VARIABLES:
    integer :: status                      ! return status
    integer :: length                      ! temporary          
    character(len=256) :: ftest,ctest      ! temporaries
    !-----------------------------------------------------------------------

    ! Continue run:
    ! Restart file pathname is read restart pointer file 

    if (nsrest==nsrContinue) then
       call restFile_read_pfile( path )
       call getfil( path, file, 0 )
    end if

    ! Branch run: 
    ! Restart file pathname is obtained from namelist "nrevsn"
    ! Check case name consistency (case name must be different for branch run, 
    ! unless namelist specification states otherwise)

    if (nsrest==nsrBranch) then
       length = len_trim(nrevsn)
       if (nrevsn(length-2:length) == '.nc') then
          path = trim(nrevsn) 
       else
          path = trim(nrevsn) // '.nc'
       end if
       call getfil( path, file, 0 )

       ! tcraig, adding xx. and .slim makes this more robust
       ctest = 'xx.'//trim(caseid)//'.slim'
       ftest = 'xx.'//trim(file)
       status = index(trim(ftest),trim(ctest))
       if (status /= 0 .and. .not.(brnch_retain_casename)) then
          if (masterproc) then
             write(iulog,*) 'Must change case name on branch run if ',&
                  'brnch_retain_casename namelist is not set'
             write(iulog,*) 'previous case filename= ',trim(file),&
                  ' current case = ',trim(caseid), &
                  ' ctest = ',trim(ctest), &
                  ' ftest = ',trim(ftest)
          end if
          call endrun(msg=errMsg(sourcefile, __LINE__)) 
       end if
    end if

  end subroutine restFile_getfile

  !-----------------------------------------------------------------------
  subroutine restFile_read_pfile( pnamer )
    !
    ! !DESCRIPTION:
    ! Setup restart file and perform necessary consistency checks
    !
    ! !USES:
    use fileutils , only : opnfil, getavu, relavu
    use clm_varctl, only : rpntfil, rpntdir, inst_suffix
    !
    ! !ARGUMENTS:
    character(len=*), intent(out) :: pnamer ! full path of restart file
    !
    ! !LOCAL VARIABLES:
    !EOP
    integer :: i                  ! indices
    integer :: nio                ! restart unit
    integer :: status             ! substring check status
    character(len=256) :: locfn   ! Restart pointer file name
    !-----------------------------------------------------------------------

    ! Obtain the restart file from the restart pointer file. 
    ! For restart runs, the restart pointer file contains the full pathname 
    ! of the restart file. For branch runs, the namelist variable 
    ! [nrevsn] contains the full pathname of the restart file. 
    ! New history files are always created for branch runs.

    if (masterproc) then
       write(iulog,*) 'Reading restart pointer file....'
    endif

    nio = getavu()
    locfn = trim(rpntdir) //'/'// trim(rpntfil)//trim(inst_suffix)
    call opnfil (locfn, nio, 'f')
    read (nio,'(a256)') pnamer
    call relavu (nio)

    if (masterproc) then
       write(iulog,*) 'Reading restart data.....'
       write(iulog,'(72a1)') ("-",i=1,60)
    end if

  end subroutine restFile_read_pfile

  !-----------------------------------------------------------------------
  subroutine restFile_closeRestart( file )
    !
    ! !DESCRIPTION:
    ! Close restart file and write restart pointer file if
    ! in write mode, otherwise just close restart file if in read mode
    !
    ! !USES:
    use clm_time_manager, only : is_last_step
    !
    ! !ARGUMENTS:
    character(len=*) , intent(in) :: file  ! local output filename
    !
    ! !CALLED FROM:
    ! subroutine restart in this module
    !
    ! !REVISION HISTORY:
    ! Author: Mariana Vertenstein
    !
    !
    ! !LOCAL VARIABLES:
    !EOP
    integer :: i                   !index
    !-----------------------------------------------------------------------

    if (masterproc) then
       write(iulog,*) 'Successfully wrote local restart file ',trim(file)
       write(iulog,'(72a1)') ("-",i=1,60)
       write(iulog,*)
    end if

  end subroutine restFile_closeRestart

  !-----------------------------------------------------------------------
  subroutine restFile_write_pfile( fnamer )
    !
    ! !DESCRIPTION:
    ! Open restart pointer file. Write names of current netcdf restart file.
    !
    ! !USES:
    use clm_varctl, only : rpntdir, rpntfil, inst_suffix
    use fileutils , only : relavu
    use fileutils , only : getavu, opnfil
    !
    ! !ARGUMENTS:
    character(len=*), intent(in) :: fnamer
    !
    ! !LOCAL VARIABLES:
    integer :: m                    ! index
    integer :: nio                  ! restart pointer file
    character(len=256) :: filename  ! local file name
    !-----------------------------------------------------------------------

    if (masterproc) then
       nio = getavu()
       filename= trim(rpntdir) //'/'// trim(rpntfil)//trim(inst_suffix)
       call opnfil( filename, nio, 'f' )

       write(nio,'(a)') fnamer
       call relavu( nio )
       write(iulog,*)'Successfully wrote local restart pointer file'
    end if

  end subroutine restFile_write_pfile

  !-----------------------------------------------------------------------
  subroutine restFile_open( flag, file, ncid )

    use clm_time_manager, only : get_nstep

    character(len=*),  intent(in) :: flag ! flag to specify read or write
    character(len=*),  intent(in) :: file ! filename
    type(file_desc_t), intent(out):: ncid ! netcdf id

    integer :: omode                              ! netCDF dummy variable
    character(len= 32) :: subname='restFile_open' ! subroutine name

    if (flag == 'write') then

       ! Create new netCDF file (in define mode) and set fill mode
       ! to "no fill" to optimize performance

       if (masterproc) then	
          write(iulog,*)
          write(iulog,*)'restFile_open: writing restart dataset at ',&
               trim(file), ' at nstep = ',get_nstep()
          write(iulog,*)
       end if
       call ncd_pio_createfile(ncid, trim(file))

    else if (flag == 'read') then

       ! Open netcdf restart file

       if (masterproc) then
          write(iulog,*) 'Reading restart dataset'
       end if
       call ncd_pio_openfile (ncid, trim(file), 0)

    end if

  end subroutine restFile_open

  !-----------------------------------------------------------------------
  character(len=256) function restFile_filename( rdate )
    !
    ! !DESCRIPTION:
    !
    ! !USES:
    use clm_varctl, only : caseid, inst_suffix
    !
    ! !ARGUMENTS:
    character(len=*), intent(in) :: rdate   ! input date for restart file name 
    !-----------------------------------------------------------------------

    restFile_filename = "./"//trim(caseid)//".slim"//trim(inst_suffix)//&
         ".r."//trim(rdate)//".nc"
    if (masterproc) then
       write(iulog,*)'writing restart file ',trim(restFile_filename),' for model date = ',rdate
    end if

  end function restFile_filename

  !------------------------------------------------------------------------
  subroutine restFile_dimset( ncid )
    !
    ! !DESCRIPTION:
    ! Read/Write initial data from/to netCDF instantaneous initial data file
    !
    ! !USES:
    use clm_time_manager     , only : get_nstep
    use clm_varctl           , only : caseid, ctitle, version, username, hostname, mml_surdat
    use clm_varctl           , only : conventions, source
    use clm_varpar           , only : numrad, nlevgrnd
    use decompMod            , only : get_proc_global
    !
    ! !ARGUMENTS:
    type(file_desc_t), intent(inout) :: ncid
    !
    ! !LOCAL VARIABLES:
    integer :: dimid               ! netCDF dimension id
    integer :: numg                ! total number of gridcells across all processors
    integer :: ier                 ! error status
    integer :: strlen_dimid        ! string dimension id
    character(len=  8) :: curdate  ! current date
    character(len=  8) :: curtime  ! current time
    character(len=256) :: str
    character(len= 32) :: subname='restFile_dimset' ! subroutine name
    !------------------------------------------------------------------------

    call get_proc_global(ng=numg)

    ! Define dimensions

    call ncd_defdim(ncid , nameg      , numg           ,  dimid)

    call ncd_defdim(ncid , 'levgrnd' , nlevgrnd       ,  dimid)
    call ncd_defdim(ncid , 'numrad'  , numrad         ,  dimid)
    call ncd_defdim(ncid , 'string_length', 64        ,  dimid)
	
	! mml add my soil dimension
    call ncd_defdim(ncid , 'mml_lev'  , 10        ,  dimid) ! mml: hard coded for six soil layers
    ! mml add my dust dimension
    call ncd_defdim(ncid , 'mml_dust'  , 4        ,  dimid) ! mml: hard coded for six soil layers
    
    ! Define global attributes

    call ncd_putatt(ncid, NCD_GLOBAL, 'Conventions', trim(conventions))
    call getdatetime(curdate, curtime)
    str = 'created on ' // curdate // ' ' // curtime
    call ncd_putatt(ncid, NCD_GLOBAL, 'history' , trim(str))
    call ncd_putatt(ncid, NCD_GLOBAL, 'username', trim(username))
    call ncd_putatt(ncid, NCD_GLOBAL, 'host'    , trim(hostname))
    call ncd_putatt(ncid, NCD_GLOBAL, 'version' , trim(version))
    call ncd_putatt(ncid, NCD_GLOBAL, 'source'  , trim(source))
    str = '$Id: restFileMod.F90 41292 2012-10-26 13:51:45Z erik $'
    call ncd_putatt(ncid, NCD_GLOBAL, 'revision_id'    , trim(str))
    call ncd_putatt(ncid, NCD_GLOBAL, 'case_title'     , trim(ctitle))
    call ncd_putatt(ncid, NCD_GLOBAL, 'case_id'        , trim(caseid))
    call ncd_putatt(ncid, NCD_GLOBAL, 'surface_dataset', trim(mml_surdat))
    call ncd_putatt(ncid, NCD_GLOBAL, 'title', 'CLM Restart information')

  end subroutine restFile_dimset

  !-----------------------------------------------------------------------
  subroutine restFile_add_flag_metadata(ncid, flag, flag_name)
    !
    ! !DESCRIPTION:
    ! Add global metadata for some logical flag
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    type(file_desc_t), intent(inout) :: ncid ! local file id
    logical          , intent(in)    :: flag ! logical flag
    character(len=*) , intent(in)    :: flag_name ! name of netcdf attribute
    !
    ! !LOCAL VARIABLES:

    character(len=*), parameter :: subname = 'restFile_add_flag_metadata'
    !-----------------------------------------------------------------------

    if (flag) then
       call ncd_putatt(ncid, ncd_global, flag_name, 'true')
    else
       call ncd_putatt(ncid, ncd_global, flag_name, 'false')
    end if

  end subroutine restFile_add_flag_metadata


  !-----------------------------------------------------------------------
  subroutine restFile_dimcheck( ncid )
    !
    ! !DESCRIPTION:
    ! Check dimensions of restart file
    !
    ! !USES:
    use decompMod,  only : get_proc_global
    use clm_varpar, only : nlevgrnd
    use clm_varctl, only : single_column, nsrest, nsrStartup
    !
    ! !ARGUMENTS:
    type(file_desc_t), intent(inout) :: ncid
    !
    ! !LOCAL VARIABLES:
    integer :: numg      ! total number of gridcells across all processors
    character(len=:), allocatable :: msg  ! diagnostic message
    character(len=32) :: subname='restFile_dimcheck' ! subroutine name
    !-----------------------------------------------------------------------

    ! Get relevant sizes

    if ( .not. single_column .or. nsrest /= nsrStartup )then
       call get_proc_global(ng=numg)
       msg = 'Did you mean to set use_init_interp = .true. in user_nl_clm?' // &
            new_line('x') // &
            '(Setting use_init_interp = .true. is needed when doing a' // &
            new_line('x') // &
            'transient run with crops using an initial conditions file from a non-transient run,' // &
            new_line('x') // &
            'or a non-transient run with crops using an initial conditions file from a transient run,' // &
            new_line('x') // &
            'or when running a resolution or configuration that differs from the initial conditions.)'
       call check_dim(ncid, nameg, numg, msg=msg)
    end if
    call check_dim(ncid, 'levgrnd' , nlevgrnd)
    ! mml add check for my dim?
    call check_dim(ncid, 'mml_lev'  , 10) 
    ! mml add check for my dust dim?
    call check_dim(ncid, 'mml_dust'  , 4)

  end subroutine restFile_dimcheck

  !-----------------------------------------------------------------------
  subroutine restFile_enddef( ncid )
    !
    ! !DESCRIPTION:
    ! Read a CLM restart file.
    ! !ARGUMENTS:
    type(file_desc_t), intent(inout) :: ncid
    !-----------------------------------------------------------------------

    call ncd_enddef(ncid)

  end subroutine restFile_enddef

  !-----------------------------------------------------------------------
  subroutine restFile_close( ncid )
    !
    ! !DESCRIPTION:
    ! Read a CLM restart file.
    !
    ! !ARGUMENTS:
    type(file_desc_t), intent(inout) :: ncid
    !
    ! !LOCAL VARIABLES:
    character(len=32) :: subname='restFile_close' ! subroutine name
    !-----------------------------------------------------------------------

    call ncd_pio_closefile(ncid)

  end subroutine restFile_close

  !-----------------------------------------------------------------------
  subroutine restFile_check_consistency(bounds, ncid)
    !
    ! !DESCRIPTION:
    ! Perform some consistency checks on the restart file
    !
    ! !ARGUMENTS:
    type(bounds_type), intent(in)    :: bounds  ! bounds
    type(file_desc_t), intent(inout) :: ncid    ! netcdf id
    !
    ! !LOCAL VARIABLES:
    logical :: check_finidat_year_consistency    ! whether to check consistency between year on finidat file and current year
    
    character(len=*), parameter :: subname = 'restFile_check_consistency'
    !-----------------------------------------------------------------------
    
    call restFile_read_consistency_nl( &
         check_finidat_year_consistency)

    if (check_finidat_year_consistency) then
       call restFile_check_year(ncid)
    end if

  end subroutine restFile_check_consistency

  !-----------------------------------------------------------------------
  subroutine restFile_read_consistency_nl( &
       check_finidat_year_consistency)

    !
    ! !DESCRIPTION:
    ! Read namelist settings related to finidat consistency checks
    !
    ! !USES:
    use fileutils      , only : getavu, relavu
    use clm_nlUtilsMod , only : find_nlgroup_name
    use controlMod     , only : NLFilename
    use shr_mpi_mod    , only : shr_mpi_bcast
    !
    ! !ARGUMENTS:
    logical, intent(out) :: check_finidat_year_consistency
    !
    ! !LOCAL VARIABLES:
    integer :: nu_nml    ! unit for namelist file
    integer :: nml_error ! namelist i/o error flag

    character(len=*), parameter :: subname = 'restFile_read_consistency_nl'
    !-----------------------------------------------------------------------

    namelist /finidat_consistency_checks/ &
         check_finidat_year_consistency

    ! Set default namelist values
    check_finidat_year_consistency = .true.

    ! Read namelist
    if (masterproc) then
       nu_nml = getavu()
       open( nu_nml, file=trim(NLFilename), status='old', iostat=nml_error )
       call find_nlgroup_name(nu_nml, 'finidat_consistency_checks', status=nml_error)
       if (nml_error == 0) then
          read(nu_nml, nml=finidat_consistency_checks,iostat=nml_error)
          if (nml_error /= 0) then
             call endrun(msg='ERROR reading finidat_consistency_checks namelist'//errMsg(sourcefile, __LINE__))
          end if
       else
          call endrun(msg='ERROR Could not find finidat_consistency_checks namelist'//errMsg(sourcefile, __LINE__))
       end if
       close(nu_nml)
       call relavu( nu_nml )
    endif

    call shr_mpi_bcast (check_finidat_year_consistency, mpicom)

    if (masterproc) then
       write(iulog,*) ' '
       write(iulog,*) 'finidat_consistency_checks settings:'
       write(iulog,nml=finidat_consistency_checks)
       write(iulog,*) ' '
    end if

  end subroutine restFile_read_consistency_nl

  !-----------------------------------------------------------------------
  subroutine restFile_check_year(ncid)
    !
    ! !DESCRIPTION:
    ! Make sure year on the restart file is consistent with the current model year
    !
    ! !USES:
    use clm_time_manager     , only : get_curr_date, get_rest_date
    use clm_varctl           , only : fname_len
    !
    ! !ARGUMENTS:
    type(file_desc_t), intent(inout) :: ncid    ! netcdf id
    !
    ! !LOCAL VARIABLES:
    logical                  :: att_found                ! whether the attribute was found on the netcdf file
    character(len=fname_len) :: flanduse_timeseries_rest ! flanduse_timeseries from the restart file
    integer                  :: year                     ! current model year
    integer                  :: mon                      ! current model month
    integer                  :: day                      ! current model day of month
    integer                  :: tod                      ! current model time of day
    integer                  :: rest_year                ! year from restart file

    character(len=*), parameter :: subname = 'restFile_check_year'
    !-----------------------------------------------------------------------
    
  end subroutine restFile_check_year

end module restFileMod



