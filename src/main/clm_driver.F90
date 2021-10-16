module clm_driver

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! This module provides the main SLIM driver physics calling sequence.  Most
  ! computations occurs over ``clumps'' of gridcells (and associated subgrid
  ! scale entities) assigned to each MPI process. Computation is further
  ! parallelized by looping over clumps on each process using shared memory OpenMP.
  !
  ! !USES:
  use shr_kind_mod           , only : r8 => shr_kind_r8
  use clm_varctl             , only : wrtdia, iulog
  use clm_varctl             , only : use_cn, use_noio
  use clm_time_manager       , only : get_nstep
  use spmdMod                , only : masterproc, mpicom
  use decompMod              , only : get_proc_clumps, get_clump_bounds, get_proc_bounds, bounds_type
  use filterMod              , only : filter_inactive_and_active
  use histFileMod            , only : hist_update_hbuf, hist_htapes_wrapup
  use restFileMod            , only : restFile_write, restFile_filename
  use abortutils             , only : endrun
  !
  use SoilBiogeochemVerticalProfileMod   , only : SoilBiogeochemVerticalProfile
  use ActiveLayerMod         , only : alt_calc
  !
  use perf_mod				! MML: this is where t_startf and t_stopf are 
  !
  use clm_instMod            , only : temperature_inst, canopystate_inst
  use clm_instMod            , only : soilstate_inst, soilbiogeochem_state_inst
  use clm_instMod            , only : bgc_vegetation_inst
  use clm_instMod            , only : atm2lnd_inst, lnd2atm_inst
  use clm_instMod            , only : soilstate_inst

  ! MML: add use simple land model module
  use mml_mainMod			 , only : mml_main	! MML if I don't say "only", it'll be fine, yes?
										!  GBB: The 'only' is not required. It means that only 'mml_main' in the
										! mml_mainMod is available to clm_driver. It is a good coding practice,
										! so that it is clear what is being used here.


  !
  ! !PUBLIC TYPES:
  implicit none
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: clm_drv            ! Main clm driver 
  !
  ! !PRIVATE MEMBER FUNCTIONS:
  private :: write_diagnostic  ! Write diagnostic information to log file

  character(len=*), parameter, private :: sourcefile = &
       __FILE__
  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine clm_drv(doalb, nextsw_cday, declinp1, declin, rstwr, nlend, rdate, rof_prognostic)
    !
    ! !DESCRIPTION:
    !
    ! First phase of the SLIM driver calling the SLIM physics.
    !
    ! !USES:
    use clm_time_manager, only : get_curr_date    
    !
    ! !ARGUMENTS:
    implicit none
    logical ,        intent(in) :: doalb       ! true if time for surface albedo calc
    real(r8),        intent(in) :: nextsw_cday ! calendar day for nstep+1
    real(r8),        intent(in) :: declinp1    ! declination angle for next time step
    real(r8),        intent(in) :: declin      ! declination angle for current time step
    logical,         intent(in) :: rstwr       ! true => write restart file this step
    logical,         intent(in) :: nlend       ! true => end of run on this step
    character(len=*),intent(in) :: rdate       ! restart file time stamp for name

    ! Whether we're running with a prognostic ROF component. This shouldn't change from
    ! timestep to timestep, but we pass it into the driver loop because it isn't available
    ! in initialization.
    logical,         intent(in) :: rof_prognostic  ! whether we're running with a prognostic ROF component
    !
    ! !LOCAL VARIABLES:
    integer              :: nstep                   ! time step number
    integer              :: nc, c, p, l, g          ! indices
    integer              :: nclumps                 ! number of clumps on this processor
    character(len=256)   :: filer                   ! restart file name
    integer              :: ier                     ! error code
    type(bounds_type)    :: bounds_clump    
    type(bounds_type)    :: bounds_proc     

    ! COMPILER_BUG(wjs, 2014-11-29, pgi 14.7) Workaround for internal compiler error with
    ! pgi 14.7 ('normalize_forall_array: non-conformable'), which appears in the call to
    ! CalcIrrigationNeeded. Simply declaring this variable makes the ICE go away.
    real(r8), allocatable :: dummy1_to_make_pgi_happy(:)
    !-----------------------------------------------------------------------

    ! Determine processor bounds and clumps for this processor

    call get_proc_bounds(bounds_proc)
    nclumps = get_proc_clumps()

! MML: I think I need this bit
    !$OMP PARALLEL DO PRIVATE (nc,bounds_clump)
    do nc = 1,nclumps
       call get_clump_bounds(nc, bounds_clump)

       ! BUG(wjs, 2014-12-15, bugz 2107) Because of the placement of the following
       ! routines (alt_calc and SoilBiogeochemVerticalProfile) in the driver sequence -
       ! they are called very early in each timestep, before weights are adjusted and
       ! filters are updated - it may be necessary for these routines to compute values
       ! over inactive as well as active points (since some inactive points may soon
       ! become active) - so that's what is done now. Currently, it seems to be okay to do
       ! this, because the variables computed here seem to only depend on quantities that
       ! are valid over inactive as well as active points.

       call t_startf("decomp_vert")
       call alt_calc(filter_inactive_and_active(nc)%num_soilc, filter_inactive_and_active(nc)%soilc, &
            temperature_inst, canopystate_inst) 

       if (use_cn) then
          call SoilBiogeochemVerticalProfile(bounds_clump                                       , &
               filter_inactive_and_active(nc)%num_soilc, filter_inactive_and_active(nc)%soilc   , &
               filter_inactive_and_active(nc)%num_soilp, filter_inactive_and_active(nc)%soilp   , &
               canopystate_inst, soilstate_inst, soilbiogeochem_state_inst)
       end if

       call t_stopf("decomp_vert")
    end do
    !$OMP END PARALLEL DO

	! ============================================================================
    ! MML: Simple Land Model Override
    ! ============================================================================
    
    ! Moved down after lnd2atm and lnd2glc because I'm already computing gridcell-averaged
    ! values, and I don't want them to get overwritten with agregated patch data, which
    ! I THINK is what lnd2atm is doing
    
    ! MML: Put call to my module here, before creating lnd2atm? Or let it create 
    ! it then change it? Or I can just make my own, but in that case I have to make sure 
    ! I give it everything it needs. I think lnd2atm (but check!) actually hands the data
    ! off to the coupler, so if thats the case I need to make my changes before hand.

	call t_startf('mml_main')
		call mml_main(bounds_clump, atm2lnd_inst, lnd2atm_inst)
	call t_stopf('mml_main')
	
	!write(iulog,*) 'MML: done with simple model, back at clm_driver'
	
	
    ! ============================================================================
    ! Write global average diagnostics to standard output
    ! ============================================================================

    nstep = get_nstep()
    if (wrtdia) call mpi_barrier(mpicom,ier)
    call t_startf('wrtdiag')
    call write_diagnostic(bounds_proc, wrtdia, nstep, lnd2atm_inst)
    call t_stopf('wrtdiag')

    ! ============================================================================
    ! Update history buffer
    ! ============================================================================


    call t_startf('hbuf')
    call hist_update_hbuf(bounds_proc)
    call t_stopf('hbuf')

    ! ============================================================================
    ! History/Restart output
    ! ============================================================================

    if (.not. use_noio) then

       call t_startf('clm_drv_io')

       ! Create history and write history tapes if appropriate
       call t_startf('clm_drv_io_htapes')

       !write(iulog,*)'MML: about to call htapes_wrapup, prepare to die, my name is inigio montoya... also wtf does it want the soilstate for? '

       ! MML workaround to try and avoid the soilstate leading to crashing - this is CLM's soil state, not SLIM's, so the values shouldn't be meaningful anyhow
       !soilstate_inst%watsat_col(bounds_proc%begc:bounds_proc%endc, :) = 0.0_r8
       !soilstate_inst%sucsat_col(bounds_proc%begc:bounds_proc%endc, :) = 0.0_r8
       !soilstate_inst%bsw_col(bounds_proc%begc:bounds_proc%endc, :) = 0.0_r8
       !soilstate_inst%hksat_col(bounds_proc%begc:bounds_proc%endc, :) = 0.0_r8

       !write(iulog,*)'MML: clobbered the soilstate_inst values, call hist_htapes_wrapup now'

       !write(iulog,*)'MML: rstwr = ',rstwr,', nlend = ',nlend

       call hist_htapes_wrapup( rstwr, nlend, bounds_proc,                    &
            soilstate_inst%watsat_col(bounds_proc%begc:bounds_proc%endc, 1:), &
            soilstate_inst%sucsat_col(bounds_proc%begc:bounds_proc%endc, 1:), &
            soilstate_inst%bsw_col(bounds_proc%begc:bounds_proc%endc, 1:),    &
            soilstate_inst%hksat_col(bounds_proc%begc:bounds_proc%endc, 1:))

       !write(iulog,*)'MML: back from wrapup, yet we are still running'

       call t_stopf('clm_drv_io_htapes')

       if (use_cn) then
          call bgc_vegetation_inst%WriteHistory(bounds_proc)
       end if

       ! Write restart/initial files if appropriate
       if (rstwr) then
          !write(iulog,*)'MML: write restart file'
          call t_startf('clm_drv_io_wrest')
          filer = restFile_filename(rdate=rdate)

          call restFile_write( bounds_proc, filer, rdate=rdate )

          call t_stopf('clm_drv_io_wrest')
          
          ! MML:
         ! write(iulog,*)  'MML: end of restart if statment '
          
       end if
       call t_stopf('clm_drv_io')

	   ! MML:
       !write(iulog,*)  'MML: after restart call '
	    
    end if
    
    ! MML:
    !write(iulog,*)  'MML: end clm_drv routine '

  end subroutine clm_drv

  !------------------------------------------------------------------------
  subroutine write_diagnostic (bounds, wrtdia, nstep, lnd2atm_inst)
    !
    ! !DESCRIPTION:
    ! Write diagnostic surface temperature output each timestep.  Written to
    ! be fast but not bit-for-bit because order of summations can change each
    ! timestep.
    !
    ! !USES:
    use decompMod  , only : get_proc_global
    use spmdMod    , only : masterproc, npes, MPI_REAL8
    use spmdMod    , only : MPI_STATUS_SIZE, mpicom, MPI_SUM
    use shr_sys_mod, only : shr_sys_flush
    use abortutils , only : endrun
    use shr_log_mod, only : errMsg => shr_log_errMsg
    use lnd2atmType, only : lnd2atm_type
    !
    ! !ARGUMENTS:
    type(bounds_type)  , intent(in) :: bounds        
    logical            , intent(in) :: wrtdia     !true => write diagnostic
    integer            , intent(in) :: nstep      !model time step
    type(lnd2atm_type) , intent(in) :: lnd2atm_inst
    !
    ! !REVISION HISTORY:
    ! Created by Mariana Vertenstein
    !
    ! !LOCAL VARIABLES:
    integer :: p                       ! loop index
    integer :: numg                    ! total number of gridcells across all processors
    integer :: ier                     ! error status
    real(r8):: psum                    ! partial sum of ts
    real(r8):: tsum                    ! sum of ts
    real(r8):: tsxyav                  ! average ts for diagnostic output
    integer :: status(MPI_STATUS_SIZE) ! mpi status
    !------------------------------------------------------------------------

    call get_proc_global(ng=numg)

    if (wrtdia) then

       call t_barrierf('sync_write_diag', mpicom)
       psum = sum(lnd2atm_inst%t_rad_grc(bounds%begg:bounds%endg))
       call mpi_reduce(psum, tsum, 1, MPI_REAL8, MPI_SUM, 0, mpicom, ier)
       if (ier/=0) then
          write(iulog,*) 'write_diagnostic: Error in mpi_reduce()'
          call endrun(msg=errMsg(sourcefile, __LINE__))
       end if
       if (masterproc) then
          tsxyav = tsum / numg
          write(iulog,1000) nstep, tsxyav
          call shr_sys_flush(iulog)
       end if

    else

       if (masterproc) then
          write(iulog,*)'slim: completed timestep ',nstep
          call shr_sys_flush(iulog)
       end if

    endif

1000 format (1x,'nstep = ',i10,'   TS = ',f21.15)

  end subroutine write_diagnostic

end module clm_driver
