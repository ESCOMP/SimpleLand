module clm_initializeMod

  !-----------------------------------------------------------------------
  ! Performs land model initialization
  !
  use shr_kind_mod    , only : r8 => shr_kind_r8
  use shr_sys_mod     , only : shr_sys_flush
  use shr_log_mod     , only : errMsg => shr_log_errMsg
  use spmdMod         , only : masterproc
  use decompMod       , only : bounds_type, get_proc_bounds, get_proc_clumps, get_clump_bounds
  use abortutils      , only : endrun
  use clm_varctl      , only : nsrest, nsrStartup, nsrContinue, nsrBranch
  use clm_varctl      , only : is_cold_start, is_interpolated_start
  use clm_varctl      , only : iulog
  use clm_instur      , only : wt_lunit, urban_valid, wt_nat_patch, wt_cft, wt_glc_mec, topo_glc_mec
  use perf_mod        , only : t_startf, t_stopf
  use ncdio_pio       , only : file_desc_t
  use GridcellType    , only : grc           ! instance     
  use LandunitType    , only : lun           ! instance          
  use ColumnType      , only : col           ! instance          
  use PatchType       , only : patch         ! instance            
  use reweightMod     , only : reweight_wrapup
  use filterMod       , only : allocFilters, filter

  use clm_instMod       
  ! 
  implicit none
  public   ! By default everything is public 
  !
  public :: initialize1  ! Phase one initialization
  public :: initialize2  ! Phase two initialization
  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine initialize1( )
    !
    ! !DESCRIPTION:
    ! CLM initialization first phase 
    !
    ! !USES:
    use clm_varpar       , only: clm_varpar_init, natpft_lb, natpft_ub, cft_lb, cft_ub
    use clm_varcon       , only: clm_varcon_init
    use landunit_varcon  , only: landunit_varcon_init, max_lunit
    use clm_varctl       , only: fsurdat, fatmlndfrc, noland, version, mml_surdat  
    use decompInitMod    , only: decompInit_lnd, decompInit_clumps, decompInit_glcp
    use domainMod        , only: domain_check, ldomain, domain_init
    use surfrdMod        , only: surfrd_get_globmask, surfrd_get_grid, surfrd_get_data 
    use controlMod       , only: control_init, control_print, NLFilename
    use ncdio_pio        , only: ncd_pio_init
    use initGridCellsMod , only: initGridCells
    use UrbanParamsType  , only: UrbanInput
    !
    ! !LOCAL VARIABLES:
    integer           :: ier                     ! error status
    integer           :: i,j,n,k,c,l,g           ! indices
    integer           :: nl                      ! gdc and glo lnd indices
    integer           :: ns, ni, nj              ! global grid sizes
    integer           :: begg, endg              ! processor bounds
    type(bounds_type) :: bounds_proc             
    type(bounds_type) :: bounds_clump
    integer           :: nclumps                 ! number of clumps on this processor
    integer           :: nc                      ! clump index
    integer ,pointer  :: amask(:)                ! global land mask
    character(len=32) :: subname = 'initialize1' ! subroutine name
    !-----------------------------------------------------------------------

    call t_startf('clm_init1')

    ! ------------------------------------------------------------------------
    ! Initialize run control variables, timestep
    ! ------------------------------------------------------------------------

    if ( masterproc )then
       write(iulog,*) trim(version)
       write(iulog,*)
       write(iulog,*) 'Attempting to initialize the land model .....'
       write(iulog,*)
       call shr_sys_flush(iulog)
    endif

    call control_init()
    call clm_varpar_init()
    call clm_varcon_init()
    call landunit_varcon_init()
    call ncd_pio_init()

    if (masterproc) call control_print()

    ! ------------------------------------------------------------------------
    ! Read in global land grid and land mask (amask)- needed to set decomposition
    ! ------------------------------------------------------------------------

    ! global memory for amask is allocate in surfrd_get_glomask - must be
    ! deallocated below
    if (masterproc) then
       write(iulog,*) 'Attempting to read global land mask from ',trim(fatmlndfrc)
       call shr_sys_flush(iulog)
    endif
    call surfrd_get_globmask(filename=fatmlndfrc, mask=amask, ni=ni, nj=nj)

    ! Exit early if no valid land points
    if ( all(amask == 0) )then
       if (masterproc) write(iulog,*) trim(subname)//': no valid land points do NOT run clm'
       noland = .true.
       return
    end if

    ! ------------------------------------------------------------------------
    ! Determine clm gridcell decomposition and processor bounds for gridcells
    ! ------------------------------------------------------------------------

    call decompInit_lnd(ni, nj, amask)
    deallocate(amask)

    ! *** Get JUST gridcell processor bounds ***
    ! Remaining bounds (landunits, columns, patches) will be determined 
    ! after the call to decompInit_glcp - so get_proc_bounds is called
    ! twice and the gridcell information is just filled in twice

    call get_proc_bounds(begg, endg)

    ! ------------------------------------------------------------------------
    ! Get grid and land fraction (set ldomain)
    ! ------------------------------------------------------------------------

    if (masterproc) then
       write(iulog,*) 'Attempting to read ldomain from ',trim(fatmlndfrc)
       call shr_sys_flush(iulog)
    endif
    call surfrd_get_grid(begg, endg, ldomain, fatmlndfrc)
    if (masterproc) then
       call domain_check(ldomain)
    endif
    ldomain%mask = 1  !!! TODO - is this needed?

    ! Initialize glc behavior
    call glc_behavior%Init(begg, endg, NLFilename)

    ! Initialize urban model input (initialize urbinp data structure)
    ! This needs to be called BEFORE the call to surfrd_get_data since
    ! that will call surfrd_get_special which in turn calls check_urban

    call UrbanInput(begg, endg, mode='initialize')

    ! Allocate surface grid dynamic memory (just gridcell bounds dependent)

    allocate (wt_lunit     (begg:endg, max_lunit           ))
    allocate (urban_valid  (begg:endg                      ))
    allocate (wt_nat_patch (begg:endg, natpft_lb:natpft_ub ))
    allocate (wt_cft       (begg:endg, cft_lb:cft_ub       ))
    allocate (wt_glc_mec  (begg:endg, 10))
    allocate (topo_glc_mec(begg:endg, 10))

    ! Read surface dataset and set up subgrid weight arrays
    
    call surfrd_get_data(begg, endg, ldomain, fsurdat)

    ! ------------------------------------------------------------------------
    ! Determine decomposition of subgrid scale landunits, columns, patches
    ! ------------------------------------------------------------------------

    call decompInit_clumps(ns, ni, nj, glc_behavior)

    ! *** Get ALL processor bounds - for gridcells, landunit, columns and patches ***

    call get_proc_bounds(bounds_proc)
    
    ! Allocate memory for subgrid data structures
    ! This is needed here BEFORE the following call to initGridcells
    ! Note that the assumption is made that none of the subgrid initialization
    ! can depend on other elements of the subgrid in the calls below

    call grc%Init  (bounds_proc%begg, bounds_proc%endg)
    call lun%Init  (bounds_proc%begl, bounds_proc%endl)
    call col%Init  (bounds_proc%begc, bounds_proc%endc)
    call patch%Init(bounds_proc%begp, bounds_proc%endp)

    ! Build hierarchy and topological info for derived types
    ! This is needed here for the following call to decompInit_glcp

    call initGridCells(glc_behavior)

    ! Set global seg maps for gridcells, landlunits, columns and patches

    call decompInit_glcp(ns, ni, nj, glc_behavior)

    ! Set filters

    call allocFilters()

    nclumps = get_proc_clumps()
    !$OMP PARALLEL DO PRIVATE (nc, bounds_clump)
    do nc = 1, nclumps
       call get_clump_bounds(nc, bounds_clump)
       call reweight_wrapup(bounds_clump)
    end do
    !$OMP END PARALLEL DO

    ! Deallocate surface grid dynamic memory for variables that aren't needed elsewhere.
    ! Some things are kept until the end of initialize2; urban_valid is kept through the
    ! end of the run for error checking.

    deallocate (wt_lunit, wt_cft, wt_glc_mec)

    call t_stopf('clm_init1')

  end subroutine initialize1


  !-----------------------------------------------------------------------
  subroutine initialize2( )
    !
    ! !DESCRIPTION:
    ! CLM initialization - second phase
    !
    ! !USES:
    use shr_scam_mod          , only : shr_scam_getCloseLatLon
    use clm_varpar            , only : nlevsno
    use clm_varcon            , only : spval
    use clm_varctl            , only : finidat, finidat_interp_source, finidat_interp_dest, fsurdat, mml_surdat
    use clm_varctl            , only : single_column, scmlat, scmlon
    use clm_time_manager      , only : get_curr_date, get_nstep, advance_timestep 
    use clm_time_manager      , only : timemgr_init, timemgr_restart_io, timemgr_restart
    use fileutils             , only : getfil
    use initInterpMod         , only : initInterp
    use subgridWeightsMod     , only : init_subgrid_weights_mod
    use histFileMod           , only : hist_htapes_build, htapes_fieldlist, hist_printflds
    use histFileMod           , only : hist_addfld1d, hist_addfld2d, no_snow_normal
    use restFileMod           , only : restFile_getfile, restFile_open, restFile_close
    use restFileMod           , only : restFile_read, restFile_write 
    use lnd2atmMod            , only : lnd2atm_minimal
    use controlMod            , only : NLFilename
    !
    ! !ARGUMENTS    
    !
    ! !LOCAL VARIABLES:
    integer               :: c,i,g,j,k,l,p! indices
    integer               :: yr           ! current year (0, ...)
    integer               :: mon          ! current month (1 -> 12)
    integer               :: day          ! current day (1 -> 31)
    integer               :: ncsec        ! current time of day [seconds]
    integer               :: nc           ! clump index
    integer               :: nclumps      ! number of clumps on this processor
    character(len=256)    :: fnamer       ! name of netcdf restart file 
    character(len=256)    :: pnamer       ! full pathname of netcdf restart file
    character(len=256)    :: locfn        ! local file name
    type(file_desc_t)     :: ncid         ! netcdf id
    integer               :: nstep        ! model time step
    type(bounds_type)     :: bounds_proc  ! processor bounds
    type(bounds_type)     :: bounds_clump ! clump bounds
    logical               :: lexist
    integer               :: closelatidx,closelonidx
    real(r8)              :: closelat,closelon
    integer               :: begp, endp
    integer               :: begc, endc
    integer               :: begl, endl
    real(r8), pointer     :: data2dptr(:,:) ! temp. pointers for slicing larger arrays
    character(len=32)     :: subname = 'initialize2' 
    !----------------------------------------------------------------------

    call t_startf('clm_init2')

    ! ------------------------------------------------------------------------
    ! Determine processor bounds and clumps for this processor
    ! ------------------------------------------------------------------------

    call get_proc_bounds(bounds_proc)
    nclumps = get_proc_clumps()

    ! ------------------------------------------------------------------------
    ! Initialize time manager
    ! ------------------------------------------------------------------------

    if (nsrest == nsrStartup) then  
       call timemgr_init()
    else
       call restFile_getfile(file=fnamer, path=pnamer)
       call restFile_open( flag='read', file=fnamer, ncid=ncid )
       call timemgr_restart_io( ncid=ncid, flag='read' )
       call restFile_close( ncid=ncid )
       call timemgr_restart()
    end if

    ! ------------------------------------------------------------------------
    ! Initialize component data structures 
    ! ------------------------------------------------------------------------

    ! Note: new logic is in place that sets all the history fields to spval so
    ! there is no guesswork in the initialization to nans of the allocated variables

    ! First put in history calls for subgrid data structures - these cannot appear in the
    ! module for the subgrid data definition due to circular dependencies that are introduced

    data2dptr => col%dz(:,-nlevsno+1:0)
    col%dz(bounds_proc%begc:bounds_proc%endc,:) = spval
    call hist_addfld2d (fname='SNO_Z', units='m', type2d='levsno',  &
         avgflag='A', long_name='Snow layer thicknesses', &
         ptr_col=data2dptr, no_snow_behavior=no_snow_normal, default='inactive')

    call hist_addfld2d (fname='SNO_Z_ICE', units='m', type2d='levsno',  &
         avgflag='A', long_name='Snow layer thicknesses (ice landunits only)', &
         ptr_col=data2dptr, no_snow_behavior=no_snow_normal, &
         l2g_scale_type='ice', default='inactive')

    col%zii(bounds_proc%begc:bounds_proc%endc) = spval
    call hist_addfld1d (fname='ZII', units='m', &
         avgflag='A', long_name='convective boundary height', &
         ptr_col=col%zii, default='inactive')

    ! If single-column determine closest latitude and longitude

    if (single_column) then
       call getfil (fsurdat, locfn, 0)
       call shr_scam_getCloseLatLon(locfn, scmlat, scmlon, &
            closelat, closelon, closelatidx, closelonidx)
    end if

    ! Initialize instances of all derived types as well as time constant variables

    call clm_instInit(bounds_proc)

    call hist_printflds()

    ! ------------------------------------------------------------------------
    ! Initializate dynamic subgrid weights (for prescribed transient Patches
    ! and/or dynamic landunits); note that these will be overwritten in a
    ! restart run
    ! ------------------------------------------------------------------------

    call t_startf('init_subgrid_weights')
    call init_subgrid_weights_mod(bounds_proc)
    call t_stopf('init_subgrid_weights')

    ! ------------------------------------------------------------------------
    ! Initialize modules (after time-manager initialization in most cases)
    ! ------------------------------------------------------------------------

    ! ------------------------------------------------------------------------
    ! On restart only - process the history namelist. 
    ! ------------------------------------------------------------------------

    ! Later the namelist from the restart file will be used.  This allows basic
    ! checking to make sure you didn't try to change the history namelist on restart.

    if (nsrest == nsrContinue ) then
       call htapes_fieldlist()
    end if

    ! ------------------------------------------------------------------------
    ! Read restart/initial info 
    ! ------------------------------------------------------------------------

    is_cold_start = .false.
    is_interpolated_start = .false.

    if (nsrest == nsrStartup) then

       if (finidat == ' ') then
          if (finidat_interp_source == ' ') then
             is_cold_start = .true.
             if (masterproc) then
                write(iulog,*)'Using cold start initial conditions '
             end if
          else 
             if (masterproc) then
                write(iulog,*)'Interpolating initial conditions from ',trim(finidat_interp_source),&
                     ' and creating new initial conditions ', trim(finidat_interp_dest)
             end if
          end if
       else 
          if (masterproc) then
             write(iulog,*)'Reading initial conditions from ',trim(finidat)
          end if
          call getfil( finidat, fnamer, 0 )
          call restFile_read(bounds_proc, fnamer, glc_behavior)
       end if

    else if ((nsrest == nsrContinue) .or. (nsrest == nsrBranch)) then

       if (masterproc) then
          write(iulog,*)'Reading restart file ',trim(fnamer)
       end if
       call restFile_read(bounds_proc, fnamer, glc_behavior)

    end if

    ! ------------------------------------------------------------------------
    ! If appropriate, create interpolated initial conditions
    ! ------------------------------------------------------------------------

    if (nsrest == nsrStartup .and. finidat_interp_source /= ' ') then

       is_interpolated_start = .true.

       ! Check that finidat is not cold start - abort if it is
       if (finidat /= ' ') then
          call endrun(msg='ERROR clm_initializeMod: '//&
               'finidat and finidat_interp_source cannot both be non-blank')
       end if

       ! Create new template file using cold start
       call restFile_write(bounds_proc, finidat_interp_dest)

       ! Interpolate finidat onto new template file
       call getfil( finidat_interp_source, fnamer,  0 )
       call initInterp(filei=fnamer, fileo=finidat_interp_dest, bounds=bounds_proc)

       ! Read new interpolated conditions file back in
       call restFile_read(bounds_proc, finidat_interp_dest, glc_behavior)

       ! Reset finidat to now be finidat_interp_dest 
       ! (to be compatible with routines still using finidat)
       finidat = trim(finidat_interp_dest)

    end if

    ! ------------------------------------------------------------------------
    ! Initialize active history fields. 
    ! ------------------------------------------------------------------------

    ! This is only done if not a restart run. If a restart run, then this 
    ! information has already been obtained from the restart data read above. 
    ! Note that routine hist_htapes_build needs time manager information,
    ! so this call must be made after the restart information has been read.

    if (nsrest /= nsrContinue) then
       call hist_htapes_build()
    end if

    ! ------------------------------------------------------------------------
    ! Initialize variables that are associated with accumulated fields.
    ! ------------------------------------------------------------------------

    ! The following is called for both initial and restart runs and must
    ! must be called after the restart file is read 
    call atm2lnd_inst%initAccVars(bounds_proc)

    !------------------------------------------------------------       
    ! Determine gridcell averaged properties to send to atm
    !------------------------------------------------------------       

    if (nsrest == nsrStartup) then
       call t_startf('init_map2gc')
       call lnd2atm_minimal(bounds_proc, waterstate_inst, lnd2atm_inst)
       call t_stopf('init_map2gc')
    end if

    !------------------------------------------------------------       
    ! Deallocate wt_nat_patch
    !------------------------------------------------------------       

    ! wt_nat_patch was allocated in initialize1, but needed to be kept around through
    ! initialize2 for some consistency checking; now it can be deallocated

    deallocate(wt_nat_patch)

    ! topo_glc_mec was allocated in initialize1, but needed to be kept around through
    ! initialize2 because it is used to initialize other variables; now it can be
    ! deallocated

    deallocate(topo_glc_mec)

    !------------------------------------------------------------       
    ! Write log output for end of initialization
    !------------------------------------------------------------       

    call t_startf('init_wlog')
    if (masterproc) then
       write(iulog,*) 'Successfully initialized the land model'
       if (nsrest == nsrStartup) then
          write(iulog,*) 'begin initial run at: '
       else
          write(iulog,*) 'begin continuation run at:'
       end if
       call get_curr_date(yr, mon, day, ncsec)
       write(iulog,*) '   nstep= ',get_nstep(), ' year= ',yr,' month= ',mon,&
            ' day= ',day,' seconds= ',ncsec
       write(iulog,*)
       write(iulog,'(72a1)') ("*",i=1,60)
       write(iulog,*)
    endif
    call t_stopf('init_wlog')

    call t_stopf('clm_init2')

  end subroutine initialize2

end module clm_initializeMod
