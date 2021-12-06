module IrrigationMod

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Calculates irrigation flux.
  !
  ! Usage:
  !
  !   - Call CalcIrrigationNeeded in order to compute whether and how much irrigation is
  !     needed for the next call to ApplyIrrigation. This should be called once per
  !     timestep.
  ! 
  !   - Call ApplyIrrigation in order to calculate qflx_irrig. This should be called
  !     exactly once per time step, before the first time qflx_irrig is needed by other
  !     parts of the code. It is acceptable for this to be called earlier in the timestep
  !     than CalcIrrigationNeeded.
  !
  !   - Access the timestep's irrigation flux via qflx_irrig_patch or
  !     qflx_irrig_col. These should be treated as read-only.
  !
  ! Design notes:
  !
  !   In principle, ApplyIrrigation and CalcIrrigationNeeded could be combined into a
  !   single routine. Their separation is largely for historical reasons: In the past,
  !   irrigation depended on btran, and qflx_irrig is needed earlier in the driver loop
  !   than when btran becomes available. (And qflx_irrig is also used late in the driver
  !   loop - so it wouldn't work, for example, to calculate qflx_irrig after btran is
  !   computed, and then save it on the restart file for the next iteration of the driver
  !   loop: then the uses of qflx_irrig early and late in the driver loop would be
  !   inconsistent.)
  !
  !   Now that we no longer have a dependency on btran, we could call CalcIrrigationNeeded
  !   before the first time qflx_irrig is needed. Thus, there might be some advantage to
  !   combining ApplyIrrigation and CalcIrrigationNeeded - or at least calling these two
  !   routines from the same place.  In particular: this separation of the irrigation
  !   calculation into two routines that are done at different times in the driver loop
  !   makes it harder and less desirable to nest the irrigation object within some other
  !   object: Doing so might make it harder to do the two separate steps at the right
  !   time, and would lead to less clarity about how these two steps are ordered with
  !   respect to the rest of the driver loop. So if we start trying to create a hierarchy
  !   of objects in CLM, we may want to rework this design. And we may want to rework it
  !   to keep things simpler anyway, even if we aren't nesting objects. Note that this
  !   rework would change behavior slightly, because irrigation would be applied in the
  !   same time step that CalcIrrigationNeeded first determines it's needed, rather than
  !   waiting until the following time step.
  !
  ! !USES:
#include "shr_assert.h"
  use shr_kind_mod     , only : r8 => shr_kind_r8
  use decompMod        , only : bounds_type, get_proc_global
  use shr_log_mod      , only : errMsg => shr_log_errMsg
  use abortutils       , only : endrun
  use clm_varctl       , only : iulog
  use clm_varcon       , only : isecspday, degpsec, denh2o, spval, namec
  use clm_varpar       , only : nlevsoi, nlevgrnd
  use clm_time_manager , only : get_step_size
  use SoilWaterRetentionCurveMod, only : soil_water_retention_curve_type
  use GridcellType     , only : grc                
  use ColumnType       , only : col                
  use PatchType        , only : patch                
  use subgridAveMod    , only : p2c, c2g
  use filterColMod     , only : filter_col_type, col_filter_from_logical_array
  !
  implicit none
  private

  ! !PUBLIC TYPES:
  
  ! This type is public (and its components are public, too) to aid unit testing
  type, public :: irrigation_params_type
     ! Minimum LAI for irrigation
     real(r8) :: irrig_min_lai

     ! Time of day to check whether we need irrigation, seconds (0 = midnight). 
     ! We start applying the irrigation in the time step FOLLOWING this time, 
     ! since we won't begin irrigating until the next call to ApplyIrrigation
     integer  :: irrig_start_time

     ! Desired amount of time to irrigate per day (sec). Actual time may 
     ! differ if this is not a multiple of dtime. Irrigation won't work properly 
     ! if dtime > secsperday
     integer  :: irrig_length

     ! Target soil matric potential for irrigation (mm)
     !
     ! When we irrigate, we aim to bring the total soil moisture in the top (irrig_depth)
     ! m of soil up to this level.
     !
     ! (Note: When we convert this to a relative saturation, we ensure that the relative
     ! saturation target is bounded by [0,1].)
     real(r8) :: irrig_target_smp

     ! Soil depth to which we measure for irrigation (m)
     real(r8) :: irrig_depth

     ! Determines soil moisture threshold at which we irrigate. If
     ! h2osoi_liq_wilting_point is the soil moisture level at wilting point and
     ! h2osoi_liq_target is the soil moisture level at the target irrigation level (given
     ! by irrig_target_smp), then the threshold at which we irrigate is
     !     h2osoi_liq_wilting_point +
     !          irrig_threshold_fraction*(h2osoi_liq_target - h2osoi_liq_wilting_point)
     ! A value of 1 means that we irrigate whenever soil moisture falls below the target
     ! A value of 0 means that we only irrigate when soil moisture falls below the
     ! wilting point
     real(r8) :: irrig_threshold_fraction

     ! Threshold for river water volume below which irrigation is shut off, if
     ! limit_irrigation is .true. (fraction of available river water). A threshold of 0
     ! means allow all river water to be used; a threshold of 0.1 means allow 90% of the
     ! river volume to be used; etc.
     real(r8) :: irrig_river_volume_threshold

     ! Whether irrigation is limited based on river storage. This only applies if ROF is
     ! enabled (i.e., rof_prognostic is .true.) - otherwise we don't limit irrigation,
     ! regardless of the value of this flag.
     logical :: limit_irrigation_if_rof_enabled

  end type irrigation_params_type


  type, public :: irrigation_type
     private
     ! Public data members
     ! Note: these should be treated as read-only by other modules
     real(r8), pointer, public :: qflx_irrig_patch(:) ! patch irrigation flux (mm H2O/s)
     real(r8), pointer, public :: qflx_irrig_col  (:) ! col irrigation flux (mm H2O/s)

     ! Private data members; set in initialization:
     type(irrigation_params_type) :: params
     integer :: dtime                ! land model time step (sec)
     integer :: irrig_nsteps_per_day ! number of time steps per day in which we irrigate
     real(r8), pointer :: relsat_wilting_point_col(:,:) ! relative saturation at which smp = wilting point [col, nlevsoi]
     real(r8), pointer :: relsat_target_col(:,:)        ! relative saturation at which smp is at the irrigation target [col, nlevsoi]

     ! Private data members; time-varying:
     real(r8), pointer :: irrig_rate_patch            (:) ! current irrigation rate [mm/s]
     real(r8), pointer :: irrig_rate_demand_patch     (:) ! current irrigation rate, neglecting surface water source limitation [mm/s]
     integer , pointer :: n_irrig_steps_left_patch    (:) ! number of time steps for which we still need to irrigate today (if 0, ignore)
     real(r8), pointer :: qflx_irrig_demand_patch     (:) ! irrigation flux neglecting surface water source limitation [mm/s]

   contains
     ! Public routines
     ! COMPILER_BUG(wjs, 2014-10-15, pgi 14.7) Add an "Irrigation" prefix to some  generic routines like "Init"
     ! (without this workaround, pgi compilation fails in restFileMod)
     procedure, public :: Init => IrrigationInit
     procedure, public :: Restart
     procedure, public :: Clean => IrrigationClean ! deallocate memory

     ! Private routines
     procedure, private :: ReadNamelist
     procedure, private :: CheckNamelistValidity   ! Check for validity of input parameters
     procedure, private :: InitAllocate => IrrigationInitAllocate
     procedure, private :: InitHistory => IrrigationInitHistory
     procedure, private :: InitCold => IrrigationInitCold
  end type irrigation_type

  interface irrigation_params_type
     module procedure irrigation_params_constructor
  end interface irrigation_params_type

  ! Soil matric potential at wilting point (mm)
  !
  ! There is no reason to make this a tunable parameter, because the behavior it governs
  ! (the trigger for irrigation) can be tuned via other parameters.
  !
  ! TODO(wjs, 2016-09-08) It looks like there is other code in CLM that also uses an
  ! assumed wilting point (CNRootDynMod, maybe others). We should probably make this a
  ! shared parameter, e.g., in clm_varcon.
  real(r8), parameter, private :: wilting_point_smp = -150000._r8

  ! Conversion factors
  real(r8), parameter :: m3_over_km2_to_mm = 1.e-3_r8

  character(len=*), parameter, private :: sourcefile = &
       __FILE__
  
contains

  ! ========================================================================
  ! Constructors
  ! ========================================================================

  !-----------------------------------------------------------------------
  function irrigation_params_constructor(irrig_min_lai, &
       irrig_start_time, irrig_length, &
       irrig_target_smp, &
       irrig_depth, irrig_threshold_fraction, irrig_river_volume_threshold, &
       limit_irrigation_if_rof_enabled) &
       result(this)
    !
    ! !DESCRIPTION:
    ! Create an irrigation_params instance
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    type(irrigation_params_type) :: this  ! function result
    real(r8), intent(in) :: irrig_min_lai
    integer , intent(in) :: irrig_start_time
    integer , intent(in) :: irrig_length
    real(r8), intent(in) :: irrig_target_smp
    real(r8), intent(in) :: irrig_depth
    real(r8), intent(in) :: irrig_threshold_fraction
    real(r8), intent(in) :: irrig_river_volume_threshold
    logical , intent(in) :: limit_irrigation_if_rof_enabled
    !
    ! !LOCAL VARIABLES:
    
    character(len=*), parameter :: subname = 'irrigation_params_constructor'
    !-----------------------------------------------------------------------
    
    this%irrig_min_lai = irrig_min_lai
    this%irrig_start_time = irrig_start_time
    this%irrig_length = irrig_length
    this%irrig_depth = irrig_depth
    this%irrig_target_smp = irrig_target_smp
    this%irrig_threshold_fraction = irrig_threshold_fraction
    this%irrig_river_volume_threshold = irrig_river_volume_threshold
    this%limit_irrigation_if_rof_enabled = limit_irrigation_if_rof_enabled

  end function irrigation_params_constructor


  ! ========================================================================
  ! Infrastructure routines (initialization, restart, etc.)
  ! ========================================================================
  
  !------------------------------------------------------------------------
  subroutine IrrigationInit(this, bounds, NLFilename, &
       soilstate_inst, soil_water_retention_curve)
    use SoilStateType , only : soilstate_type

    class(irrigation_type) , intent(inout) :: this
    type(bounds_type)      , intent(in)    :: bounds
    character(len=*)       , intent(in)    :: NLFilename ! Namelist filename
    type(soilstate_type)   , intent(in)    :: soilstate_inst
    class(soil_water_retention_curve_type), intent(in) :: soil_water_retention_curve

    call this%ReadNamelist(NLFilename)
    call this%InitAllocate(bounds)
    call this%InitHistory(bounds)
    call this%InitCold(bounds, soilstate_inst, soil_water_retention_curve)
  end subroutine IrrigationInit

  !-----------------------------------------------------------------------
  subroutine ReadNamelist(this, NLFilename)
    !
    ! !DESCRIPTION:
    ! Read the irrigation namelist
    !
    ! !USES:
    use fileutils      , only : getavu, relavu, opnfil
    use shr_nl_mod     , only : shr_nl_find_group_name
    use spmdMod        , only : masterproc, mpicom
    use shr_mpi_mod    , only : shr_mpi_bcast
    use shr_infnan_mod , only : nan => shr_infnan_nan, assignment(=)
    !
    ! !ARGUMENTS:
    character(len=*), intent(in) :: NLFilename ! Namelist filename
    class(irrigation_type) , intent(inout) :: this
    !
    ! !LOCAL VARIABLES:

    ! temporary variables corresponding to the components of irrigation_params_type
    real(r8) :: irrig_min_lai
    integer  :: irrig_start_time
    integer  :: irrig_length
    real(r8) :: irrig_target_smp
    real(r8) :: irrig_depth
    real(r8) :: irrig_threshold_fraction
    real(r8) :: irrig_river_volume_threshold
    logical  :: limit_irrigation_if_rof_enabled

    integer :: ierr                 ! error code
    integer :: unitn                ! unit for namelist file
    character(len=*), parameter :: nmlname = 'irrigation_inparm'

    character(len=*), parameter :: subname = 'ReadNamelist'
    !-----------------------------------------------------------------------

    namelist /irrigation_inparm/ irrig_min_lai, irrig_start_time, irrig_length, &
         irrig_target_smp, irrig_depth, irrig_threshold_fraction, &
         irrig_river_volume_threshold, limit_irrigation_if_rof_enabled

    ! Initialize options to garbage defaults, forcing all to be specified explicitly in
    ! order to get reasonable results
    irrig_min_lai = nan
    irrig_start_time = 0
    irrig_length = 0
    irrig_target_smp = nan
    irrig_depth = nan
    irrig_threshold_fraction = nan
    irrig_river_volume_threshold = nan
    limit_irrigation_if_rof_enabled = .false.

    if (masterproc) then
       unitn = getavu()
       write(iulog,*) 'Read in '//nmlname//'  namelist'
       call opnfil (NLFilename, unitn, 'F')
       call shr_nl_find_group_name(unitn, nmlname, status=ierr)
       if (ierr == 0) then
          read(unitn, nml=irrigation_inparm, iostat=ierr)
          if (ierr /= 0) then
             call endrun(msg="ERROR reading "//nmlname//"namelist"//errmsg(sourcefile, __LINE__))
          end if
       else
          call endrun(msg="ERROR could NOT find "//nmlname//"namelist"//errmsg(sourcefile, __LINE__))
       end if
       call relavu( unitn )
    end if

    call shr_mpi_bcast(irrig_min_lai, mpicom)
    call shr_mpi_bcast(irrig_start_time, mpicom)
    call shr_mpi_bcast(irrig_length, mpicom)
    call shr_mpi_bcast(irrig_target_smp, mpicom)
    call shr_mpi_bcast(irrig_depth, mpicom)
    call shr_mpi_bcast(irrig_threshold_fraction, mpicom)
    call shr_mpi_bcast(irrig_river_volume_threshold, mpicom)
    call shr_mpi_bcast(limit_irrigation_if_rof_enabled, mpicom)

    this%params = irrigation_params_type( &
         irrig_min_lai = irrig_min_lai, &
         irrig_start_time = irrig_start_time, &
         irrig_length = irrig_length, &
         irrig_target_smp = irrig_target_smp, &
         irrig_depth = irrig_depth, &
         irrig_threshold_fraction = irrig_threshold_fraction, &
         irrig_river_volume_threshold = irrig_river_volume_threshold, &
         limit_irrigation_if_rof_enabled = limit_irrigation_if_rof_enabled)

    if (masterproc) then
       write(iulog,*) ' '
       write(iulog,*) nmlname//' settings:'
       ! Write settings one-by-one rather than with a nml write because
       ! irrig_river_volume_threshold may be NaN
       write(iulog,*) 'irrig_min_lai = ', irrig_min_lai
       write(iulog,*) 'irrig_start_time = ', irrig_start_time
       write(iulog,*) 'irrig_length = ', irrig_length
       write(iulog,*) 'irrig_target_smp = ', irrig_target_smp
       write(iulog,*) 'irrig_depth = ', irrig_depth
       write(iulog,*) 'irrig_threshold_fraction = ', irrig_threshold_fraction
       write(iulog,*) 'limit_irrigation_if_rof_enabled = ', limit_irrigation_if_rof_enabled
       if (limit_irrigation_if_rof_enabled) then
          write(iulog,*) 'irrig_river_volume_threshold = ', irrig_river_volume_threshold
       end if
       write(iulog,*) ' '

       call this%CheckNamelistValidity()
    end if

  end subroutine ReadNamelist

  !-----------------------------------------------------------------------
  subroutine CheckNamelistValidity(this)
    !
    ! !DESCRIPTION:
    ! Check for validity of input parameters.
    !
    ! Assumes that the inputs have already been packed into 'this%params'.
    !
    ! Only needs to be called by the master task, since parameters are the same for all
    ! tasks.
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    class(irrigation_type), intent(in) :: this
    !
    ! !LOCAL VARIABLES:

    character(len=*), parameter :: subname = 'CheckNamelistValidity'
    !-----------------------------------------------------------------------

    associate( &
         irrig_min_lai => this%params%irrig_min_lai, &
         irrig_start_time => this%params%irrig_start_time, &
         irrig_length => this%params%irrig_length, &
         irrig_target_smp => this%params%irrig_target_smp, &
         irrig_depth => this%params%irrig_depth, &
         irrig_threshold_fraction => this%params%irrig_threshold_fraction, &
         irrig_river_volume_threshold => this%params%irrig_river_volume_threshold, &
         limit_irrigation_if_rof_enabled => this%params%limit_irrigation_if_rof_enabled)

    if (irrig_min_lai < 0._r8) then
       write(iulog,*) ' ERROR: irrig_min_lai must be >= 0'
       write(iulog,*) ' irrig_min_lai = ', irrig_min_lai
       call endrun(msg=' ERROR: irrig_min_lai must be >= 0 ' // errMsg(sourcefile, __LINE__))
    end if

    if (irrig_start_time < 0 .or. irrig_start_time >= isecspday) then
       write(iulog,*) ' ERROR: irrig_start_time must be >= 0 and < ', isecspday
       write(iulog,*) ' irrig_start_time = ', irrig_start_time
       call endrun(msg=' ERROR: irrig_start_time out of bounds ' // errMsg(sourcefile, __LINE__))
    end if

    if (irrig_length <= 0 .or. irrig_length > isecspday) then
       write(iulog,*) ' ERROR: irrig_length must be > 0 and <= ', isecspday
       write(iulog,*) ' irrig_length = ', irrig_length
       call endrun(msg=' ERROR: irrig_length out of bounds ' // errMsg(sourcefile, __LINE__))
    end if

    if (irrig_target_smp >= 0._r8) then
       write(iulog,*) ' ERROR: irrig_target_smp must be negative'
       write(iulog,*) ' irrig_target_smp = ', irrig_target_smp
       call endrun(msg=' ERROR: irrig_target_smp must be negative ' // errMsg(sourcefile, __LINE__))
    end if

    if (irrig_target_smp < wilting_point_smp) then
       write(iulog,*) ' ERROR: irrig_target_smp must be >= wilting_point_smp'
       write(iulog,*) ' irrig_target_smp (from namelist) = ', irrig_target_smp
       write(iulog,*) ' wilting_point_smp (hard-coded) = ', wilting_point_smp
       call endrun(msg=' ERROR: irrig_target_smp must be >= wilting_point_smp ' // errMsg(sourcefile, __LINE__))
    end if

    if (irrig_depth < 0._r8) then
       write(iulog,*) ' ERROR: irrig_depth must be > 0'
       write(iulog,*) ' irrig_depth = ', irrig_depth
       call endrun(msg=' ERROR: irrig_depth must be > 0 ' // errMsg(sourcefile, __LINE__))
    end if

    if (irrig_threshold_fraction < 0._r8 .or. irrig_threshold_fraction > 1._r8) then
       write(iulog,*) ' ERROR: irrig_threshold_fraction must be between 0 and 1'
       write(iulog,*) ' irrig_threshold_fraction = ', irrig_threshold_fraction
       call endrun(msg=' ERROR: irrig_threshold_fraction must be between 0 and 1 ' // &
            errMsg(sourcefile, __LINE__))
    end if

    if (limit_irrigation_if_rof_enabled) then
       if (irrig_river_volume_threshold < 0._r8 .or. irrig_river_volume_threshold > 1._r8) then
          write(iulog,*) ' ERROR: irrig_river_volume_threshold must be between 0 and 1'
          write(iulog,*) ' irrig_river_volume_threshold = ', irrig_river_volume_threshold
          call endrun(msg=' ERROR: irrig_river_volume_threshold must be between 0 and 1 ' // &
               errMsg(sourcefile, __LINE__))
       end if
    end if

    end associate

  end subroutine CheckNamelistValidity


  !-----------------------------------------------------------------------
  subroutine IrrigationInitAllocate(this, bounds)
    !
    ! !DESCRIPTION:
    ! Initialize irrigation data structure
    !
    ! !USES:
    use shr_infnan_mod , only : nan => shr_infnan_nan, assignment(=)
    !
    ! !ARGUMENTS:
    class(irrigation_type) , intent(inout) :: this
    type(bounds_type)      , intent(in)    :: bounds
    !
    ! !LOCAL VARIABLES:
    integer :: begp, endp
    integer :: begc, endc

    character(len=*), parameter :: subname = 'InitAllocate'
    !-----------------------------------------------------------------------

    begp = bounds%begp; endp= bounds%endp
    begc = bounds%begc; endc= bounds%endc

    allocate(this%qflx_irrig_patch         (begp:endp))          ; this%qflx_irrig_patch         (:)   = nan
    allocate(this%qflx_irrig_demand_patch  (begp:endp))          ; this%qflx_irrig_demand_patch  (:)   = nan
    allocate(this%qflx_irrig_col           (begc:endc))          ; this%qflx_irrig_col           (:)   = nan
    allocate(this%relsat_wilting_point_col (begc:endc,nlevsoi)) ; this%relsat_wilting_point_col (:,:) = nan
    allocate(this%relsat_target_col        (begc:endc,nlevsoi)) ; this%relsat_target_col        (:,:) = nan
    allocate(this%irrig_rate_patch         (begp:endp))          ; this%irrig_rate_patch         (:)   = nan
    allocate(this%irrig_rate_demand_patch  (begp:endp))          ; this%irrig_rate_demand_patch  (:)   = nan
    allocate(this%n_irrig_steps_left_patch (begp:endp))          ; this%n_irrig_steps_left_patch (:)   = 0

  end subroutine IrrigationInitAllocate

  !-----------------------------------------------------------------------
  subroutine IrrigationInitHistory(this, bounds)
    !
    ! !DESCRIPTION:
    ! Initialize irrigation history fields
    !
    ! !USES:
    use histFileMod  , only : hist_addfld1d
    !
    ! !ARGUMENTS:
    class(irrigation_type) , intent(inout) :: this
    type(bounds_type)      , intent(in)    :: bounds
    !
    ! !LOCAL VARIABLES:
    integer :: begp, endp
    
    character(len=*), parameter :: subname = 'InitHistory'
    !-----------------------------------------------------------------------

    begp = bounds%begp; endp= bounds%endp

    this%qflx_irrig_patch(begp:endp) = spval
    call hist_addfld1d (fname='QIRRIG', units='mm/s', &
         avgflag='A', long_name='water added through irrigation', &
         ptr_patch=this%qflx_irrig_patch)

    this%qflx_irrig_demand_patch(begp:endp) = spval
    call hist_addfld1d (fname='QIRRIG_DEMAND', units='mm/s', &
         avgflag='A', long_name='irrigation demand', &
         ptr_patch=this%qflx_irrig_demand_patch, default='inactive')

  end subroutine IrrigationInitHistory

  !-----------------------------------------------------------------------
  subroutine IrrigationInitCold(this, bounds, soilstate_inst, soil_water_retention_curve)
    !
    ! !DESCRIPTION:
    ! Do cold-start initialization for irrigation data structure
    !
    ! !USES:
    use SoilStateType , only : soilstate_type
    !
    ! !ARGUMENTS:
    class(irrigation_type) , intent(inout) :: this
    type(bounds_type)      , intent(in)    :: bounds
    type(soilstate_type)   , intent(in)    :: soilstate_inst
    class(soil_water_retention_curve_type), intent(in) :: soil_water_retention_curve
    !
    ! !LOCAL VARIABLES:
    integer :: c ! col index
    integer :: j ! level index

    character(len=*), parameter :: subname = 'InitCold'
    !-----------------------------------------------------------------------

    do j = 1, nlevsoi
       do c = bounds%begc, bounds%endc
          call soil_water_retention_curve%soil_suction_inverse( &
               c = c, &
               j = j, &
               smp_target = wilting_point_smp, &
               soilstate_inst = soilstate_inst, &
               s_target = this%relsat_wilting_point_col(c,j))

          call soil_water_retention_curve%soil_suction_inverse( &
               c = c, &
               j = j, &
               smp_target = this%params%irrig_target_smp, &
               soilstate_inst = soilstate_inst, &
               s_target = this%relsat_target_col(c,j))

          ! Make sure relative saturation targets are bounded by [0,1]
          !
          ! NOTE(wjs, 2016-11-17) These targets can be > 1 if smp_target is too small of
          ! a negative value; we want to force the target to a relative saturation of 1
          ! in that case. I don't see how these targets could end up negative, though I
          ! had a note that I saw negative values at one point. In practice, for
          ! reasonable parameter values, it seems that these min and max functions have
          ! no effect.
          this%relsat_wilting_point_col(c,j) = min(this%relsat_wilting_point_col(c,j), 1._r8)
          this%relsat_wilting_point_col(c,j) = max(this%relsat_wilting_point_col(c,j), 0._r8)
          this%relsat_target_col(c,j) = min(this%relsat_target_col(c,j), 1._r8)
          this%relsat_target_col(c,j) = max(this%relsat_target_col(c,j), 0._r8)
       end do
    end do

    this%dtime = get_step_size()
    this%irrig_nsteps_per_day = 0

  end subroutine IrrigationInitCold

  !-----------------------------------------------------------------------
  subroutine Restart(this, bounds, ncid, flag)
    !
    ! !DESCRIPTION:
    ! Handle restart of irrigation variables
    !
    ! !USES:
    use ncdio_pio        , only : file_desc_t, ncd_inqvdlen, ncd_double, ncd_int
    use restUtilMod
    !
    ! !ARGUMENTS:
    class(irrigation_type) :: this
    type(bounds_type), intent(in)    :: bounds 
    type(file_desc_t), intent(inout) :: ncid   ! netcdf id
    character(len=*) , intent(in)    :: flag   ! 'read', 'write' or 'define'
    !
    ! !LOCAL VARIABLES:
    logical :: do_io
    integer :: dimlen       ! dimension length
    integer :: nump_global  ! total number of patchs, globally
    integer :: err_code     ! error code
    logical :: readvar      ! determine if variable is on initial file

    character(len=*), parameter :: subname = 'Restart'
    !-----------------------------------------------------------------------
    
    ! Get expected total number of points, for later error checks
    call get_proc_global(np=nump_global)

    do_io = .true.
    readvar = .false.
    if (flag == 'read') then
       ! BACKWARDS_COMPATIBILITY
       ! On a read, confirm that this variable has the expected size; if not, don't read
       ! it (instead give it a default value). This is needed to support older initial
       ! conditions for which this variable had a different size.
       call ncd_inqvdlen(ncid, 'n_irrig_steps_left', 1, dimlen, err_code)
       if (dimlen /= nump_global) then
          do_io = .false.
       end if
    end if
    if (do_io) then
       call restartvar(ncid=ncid, flag=flag, varname='n_irrig_steps_left', xtype=ncd_int,  &
            dim1name='pft', &
            long_name='number of irrigation time steps left', units='#', &
            interpinic_flag='interp', readvar=readvar, data=this%n_irrig_steps_left_patch)
    end if
    if (flag=='read' .and. .not. readvar) then
       this%n_irrig_steps_left_patch = 0
    end if

    do_io = .true.
    readvar = .false.
    if (flag == 'read') then
       ! BACKWARDS_COMPATIBILITY
       ! On a read, confirm that this variable has the expected size; if not, don't read
       ! it (instead give it a default value). This is needed to support older initial
       ! conditions for which this variable had a different size.
       call ncd_inqvdlen(ncid, 'irrig_rate', 1, dimlen, err_code)
       if (dimlen /= nump_global) then
          do_io = .false.
       end if
    end if
    if (do_io) then
       call restartvar(ncid=ncid, flag=flag, varname='irrig_rate', xtype=ncd_double,  &
            dim1name='pft', &
            long_name='irrigation rate', units='mm/s', &
            interpinic_flag='interp', readvar=readvar, data=this%irrig_rate_patch)
    end if
    if (flag=='read' .and. .not. readvar) then
       this%irrig_rate_patch = 0.0_r8
    end if

    ! BACKWARDS_COMPATIBILITY(wjs, 2016-11-23) To support older restart files without an
    ! irrig_rate_demand field, get irrig_rate_demand from irrig_rate. I'm abusing the
    ! capability to specify multiple variable names here (even though irrig_rate isn't
    ! really an old version of irrig_rate_demand), rather than having code like 'if
    ! (.not. readvar) then ...', because the latter doesn't work when using init_interp.
    call restartvar(ncid=ncid, flag=flag, varname='irrig_rate_demand:irrig_rate', &
         xtype=ncd_double,  &
         dim1name='pft', &
         long_name='irrigation rate demand, neglecting surface water source limitation', &
         units='mm/s', &
         interpinic_flag='interp', readvar=readvar, data=this%irrig_rate_demand_patch)

  end subroutine Restart

  !-----------------------------------------------------------------------
  subroutine IrrigationClean(this)
    !
    ! !DESCRIPTION:
    ! Deallocate memory
    !
    ! !ARGUMENTS:
    class(irrigation_type), intent(inout) :: this
    !
    ! !LOCAL VARIABLES:
    
    character(len=*), parameter :: subname = 'Clean'
    !-----------------------------------------------------------------------
    
    deallocate(this%qflx_irrig_patch)
    deallocate(this%qflx_irrig_demand_patch)
    deallocate(this%qflx_irrig_col)
    deallocate(this%relsat_wilting_point_col)
    deallocate(this%relsat_target_col)
    deallocate(this%irrig_rate_patch)
    deallocate(this%irrig_rate_demand_patch)
    deallocate(this%n_irrig_steps_left_patch)

  end subroutine IrrigationClean

end module IrrigationMod
