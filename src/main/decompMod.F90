module decompMod

  !------------------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Module provides a descomposition into a clumped data structure which can
  ! be mapped back to atmosphere physics chunks.
  !
  ! !USES:
  use shr_kind_mod, only : r8 => shr_kind_r8
  ! Must use shr_sys_abort rather than endrun here to avoid circular dependency
  use shr_sys_mod , only : shr_sys_abort 
  use clm_varctl  , only : iulog
  use clm_varcon  , only : grlnd, nameg
  use mct_mod     , only : mct_gsMap
  !
  ! !PUBLIC TYPES:
  implicit none
  integer, public :: clump_pproc ! number of clumps per MPI process

  ! Define possible bounds subgrid levels
  integer, parameter, public :: BOUNDS_SUBGRID_GRIDCELL = 1

  ! Define possible bounds levels
  integer, parameter, public :: BOUNDS_LEVEL_PROC  = 1
  integer, parameter, public :: BOUNDS_LEVEL_CLUMP = 2
  !
  ! !PUBLIC MEMBER FUNCTIONS:

  public get_beg            ! get beg bound for a given subgrid level
  public get_end            ! get end bound for a given subgrid level
  public get_proc_clumps    ! number of clumps for this processor
  public get_proc_total     ! total no. of gridcells for any processor
  public get_proc_global    ! total gridcells across all processors
  public get_clmlevel_gsize ! get global size associated with clmlevel
  public get_clmlevel_gsmap ! get gsmap associated with clmlevel

  interface get_clump_bounds
     module procedure get_clump_bounds_old
     module procedure get_clump_bounds_new
  end interface
  public get_clump_bounds   ! clump beg and end gridcell

  interface get_proc_bounds
     module procedure get_proc_bounds_old
     module procedure get_proc_bounds_new
  end interface
  public get_proc_bounds    ! this processor beg and end gridcell

  ! !PRIVATE MEMBER FUNCTIONS:
  !
  ! !PRIVATE TYPES:
  private  ! (now mostly public for decompinitmod)

  integer,public :: nclumps     ! total number of clumps across all processors
  integer,public :: numg        ! total number of gridcells on all procs

  type bounds_type
     integer :: begg, endg       ! beginning and ending gridcell index
     integer :: level            ! whether defined on the proc or clump level
     integer :: clump_index      ! if defined on the clump level, this gives the clump index
  end type bounds_type
  public bounds_type

  !---global information on each pe
  type processor_type
     integer :: nclumps          ! number of clumps for processor_type iam
     integer,pointer :: cid(:)   ! clump indices
     integer :: ncells           ! number of gridcells in proc
     integer :: begg, endg       ! beginning and ending gridcell index
  end type processor_type
  public processor_type
  type(processor_type),public :: procinfo

  !---global information on each pe
  type clump_type
     integer :: owner            ! process id owning clump
     integer :: ncells           ! number of gridcells in clump
     integer :: begg, endg       ! beginning and ending gridcell index
  end type clump_type
  public clump_type
  type(clump_type),public, allocatable :: clumps(:)

  !---global information on each pe
  !--- glo = 1d global sn ordered
  !--- gdc = 1d global dc ordered compressed
  type decomp_type
     integer,pointer :: gdc2glo(:)    ! 1d gdc to 1d glo
  end type decomp_type
  public decomp_type
  type(decomp_type),public,target :: ldecomp

  type(mct_gsMap)  ,public,target :: gsMap_lnd_gdc2glo
  type(mct_gsMap)  ,public,target :: gsMap_gce_gdc2glo
  !------------------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  pure function get_beg(bounds, subgrid_level) result(beg_index)
    !
    ! !DESCRIPTION:
    ! Get beginning bounds for a given subgrid level
    !
    ! subgrid_level should be one of the constants defined in this module:
    ! BOUNDS_SUBGRID_GRIDCELL, etc.
    !
    ! Returns -1 for invalid subgrid_level (does not abort in this case, in order to keep
    ! this function pure).
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    integer :: beg_index  ! function result
    type(bounds_type), intent(in) :: bounds
    integer, intent(in) :: subgrid_level
    !
    ! !LOCAL VARIABLES:

    character(len=*), parameter :: subname = 'get_beg'
    !-----------------------------------------------------------------------

    select case (subgrid_level)
    case (BOUNDS_SUBGRID_GRIDCELL)
       beg_index = bounds%begg
    case default
       beg_index = -1
    end select

  end function get_beg

  !-----------------------------------------------------------------------
  pure function get_end(bounds, subgrid_level) result(end_index)
    !
    ! !DESCRIPTION:
    ! Get end bounds for a given subgrid level
    !
    ! subgrid_level should be one of the constants defined in this module:
    ! BOUNDS_SUBGRID_GRIDCELL, etc.
    !
    ! Returns -1 for invalid subgrid_level (does not abort in this case, in order to keep
    ! this function pure).
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    integer :: end_index  ! function result
    type(bounds_type), intent(in) :: bounds
    integer, intent(in) :: subgrid_level
    !
    ! !LOCAL VARIABLES:

    character(len=*), parameter :: subname = 'get_end'
    !-----------------------------------------------------------------------

    select case (subgrid_level)
    case (BOUNDS_SUBGRID_GRIDCELL)
       end_index = bounds%endg
    case default
       end_index = -1
    end select

  end function get_end

  !------------------------------------------------------------------------------
   subroutine get_clump_bounds_new (n, bounds)
     !
     ! !DESCRIPTION:
     ! Determine clump bounds
     !
     ! !ARGUMENTS:
     integer, intent(in)  :: n                ! processor clump index
     type(bounds_type), intent(out) :: bounds ! clump bounds
     !
     ! !LOCAL VARIABLES:
     character(len=32), parameter :: subname = 'get_clump_bounds'  ! Subroutine name
     integer :: cid                                                ! clump id
#ifdef _OPENMP
     integer, external :: OMP_GET_MAX_THREADS
     integer, external :: OMP_GET_NUM_THREADS
     integer, external :: OMP_GET_THREAD_NUM
#endif
     !------------------------------------------------------------------------------
     !    Make sure this IS being called from a threaded region
#ifdef _OPENMP
     if ( OMP_GET_NUM_THREADS() == 1 .and. OMP_GET_MAX_THREADS() > 1 )then
        call shr_sys_abort( trim(subname)//' ERROR: Calling from inside a non-threaded region)')
     end if
#endif

     cid  = procinfo%cid(n)
     bounds%begg = clumps(cid)%begg
     bounds%endg = clumps(cid)%endg
     
     bounds%level = BOUNDS_LEVEL_CLUMP
     bounds%clump_index = n

   end subroutine get_clump_bounds_new

   !------------------------------------------------------------------------------
   subroutine get_clump_bounds_old (n, begg, endg)
     integer, intent(in)  :: n           ! proc clump index
     integer, intent(out) :: begg, endg  ! clump beg and end gridcell indices
     integer :: cid                                                ! clump id
     !------------------------------------------------------------------------------

     cid  = procinfo%cid(n)
     begg = clumps(cid)%begg
     endg = clumps(cid)%endg
   end subroutine get_clump_bounds_old

   !------------------------------------------------------------------------------
   subroutine get_proc_bounds_new (bounds)
     !
     ! !DESCRIPTION:
     ! Retrieve processor bounds
     !
     ! !ARGUMENTS:
     type(bounds_type), intent(out) :: bounds ! processor bounds bounds
     !
     ! !LOCAL VARIABLES:
#ifdef _OPENMP
     integer, external :: OMP_GET_NUM_THREADS
     integer, external :: OMP_GET_MAX_THREADS
     integer, external :: OMP_GET_THREAD_NUM
#endif
     character(len=32), parameter :: subname = 'get_proc_bounds'  ! Subroutine name
     !------------------------------------------------------------------------------
     !    Make sure this is NOT being called from a threaded region
#ifdef _OPENMP
     if ( OMP_GET_NUM_THREADS() > 1 )then
        call shr_sys_abort( trim(subname)//' ERROR: Calling from inside  a threaded region')
     end if
#endif

     bounds%begg = procinfo%begg
     bounds%endg = procinfo%endg

     bounds%level = BOUNDS_LEVEL_PROC
     bounds%clump_index = -1           ! irrelevant for proc, so assigned a bogus value

   end subroutine get_proc_bounds_new

   !------------------------------------------------------------------------------
   subroutine get_proc_bounds_old (begg, endg)

     integer, optional, intent(out) :: begg, endg  ! proc beg and end gridcell indices
     !------------------------------------------------------------------------------

     if (present(begg)) begg = procinfo%begg
     if (present(endg)) endg = procinfo%endg
   end subroutine get_proc_bounds_old

   !------------------------------------------------------------------------------
   subroutine get_proc_total(pid, ncells)
     !
     ! !DESCRIPTION:
     ! Count up gridcells on process.
     !
     ! !ARGUMENTS:
     integer, intent(in)  :: pid     ! proc id
     integer, intent(out) :: ncells  ! total number of gridcells on the processor
     !
     ! !LOCAL VARIABLES:
     integer :: cid       ! clump index
     !------------------------------------------------------------------------------

     ncells   = 0
     do cid = 1,nclumps
        if (clumps(cid)%owner == pid) then
           ncells  = ncells    + clumps(cid)%ncells
        end if
     end do
   end subroutine get_proc_total

   !------------------------------------------------------------------------------
   subroutine get_proc_global(ng)
     !
     ! !DESCRIPTION:
     ! Return number of gridcells across all processes.
     !
     ! !ARGUMENTS:
     integer, optional, intent(out) :: ng        ! total number of gridcells across all processors
     !------------------------------------------------------------------------------

     if (present(ng)) ng             = numg

   end subroutine get_proc_global

   !------------------------------------------------------------------------------
   integer function get_proc_clumps()
     !
     ! !DESCRIPTION:
     ! Return the number of clumps.
     !------------------------------------------------------------------------------

     get_proc_clumps = procinfo%nclumps

   end function get_proc_clumps

   !-----------------------------------------------------------------------
   integer function get_clmlevel_gsize (clmlevel)
     !
     ! !DESCRIPTION:
     ! Determine 1d size from clmlevel
     !
     ! !USES:
     use domainMod , only : ldomain
     !
     ! !ARGUMENTS:
     character(len=*), intent(in) :: clmlevel    !type of clm 1d array
     !-----------------------------------------------------------------------

     select case (clmlevel)
     case(grlnd)
        get_clmlevel_gsize = ldomain%ns
     case(nameg)
        get_clmlevel_gsize = numg
     case default
        write(iulog,*) 'get_clmlevel_gsize does not match clmlevel type: ', trim(clmlevel)
        call shr_sys_abort()
     end select

   end function get_clmlevel_gsize

   !-----------------------------------------------------------------------
   subroutine get_clmlevel_gsmap (clmlevel, gsmap)
     !
     ! !DESCRIPTION:
     ! Compute arguments for gatherv, scatterv for vectors
     !
     ! !ARGUMENTS:
     character(len=*), intent(in) :: clmlevel     ! type of input data
     type(mct_gsmap) , pointer    :: gsmap
     !----------------------------------------------------------------------

    select case (clmlevel)
    case(grlnd)
       gsmap => gsmap_lnd_gdc2glo
    case(nameg)
       gsmap => gsmap_gce_gdc2glo
    case default
       write(iulog,*) 'get_clmlevel_gsmap: Invalid expansion character: ',trim(clmlevel)
       call shr_sys_abort()
    end select

  end subroutine get_clmlevel_gsmap

end module decompMod
