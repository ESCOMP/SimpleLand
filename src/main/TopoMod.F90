module TopoMod

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Handles topographic height of each column
  !
  ! !USES:
  use shr_kind_mod   , only : r8 => shr_kind_r8
  use shr_infnan_mod , only : nan => shr_infnan_nan, assignment(=)
  use decompMod      , only : bounds_type
  use PatchType      , only : patch
  use ColumnType     , only : col
  use LandunitType   , only : lun
  use glcBehaviorMod , only : glc_behavior_type
  use landunit_varcon, only : istice_mec
  use filterColMod   , only : filter_col_type, col_filter_from_logical_array_active_only
  !
  ! !PUBLIC TYPES:
  implicit none
  save
  private

  type, public :: topo_type
     private

     ! Public member data

     real(r8), pointer, public :: topo_col(:)  ! surface elevation (m)

     ! Private member data

     logical, pointer :: needs_downscaling_col(:)  ! whether a column needs to be downscaled

   contains
     procedure, public :: Init
     procedure, public :: Restart
     procedure, public :: Clean

     procedure, private :: InitAllocate
     procedure, private :: InitCold
  end type topo_type

contains

  !-----------------------------------------------------------------------
  subroutine Init(this, bounds)
    ! !ARGUMENTS:
    class(topo_type), intent(inout) :: this
    type(bounds_type), intent(in) :: bounds
    !
    ! !LOCAL VARIABLES:

    character(len=*), parameter :: subname = 'Init'
    !-----------------------------------------------------------------------

    call this%InitAllocate(bounds)
    call this%InitCold(bounds)

  end subroutine Init

  !-----------------------------------------------------------------------
  subroutine InitAllocate(this, bounds)
    ! !ARGUMENTS:
    class(topo_type), intent(inout) :: this
    type(bounds_type), intent(in) :: bounds
    !
    ! !LOCAL VARIABLES:
    integer :: begc, endc

    character(len=*), parameter :: subname = 'InitAllocate'
    !-----------------------------------------------------------------------

    begc = bounds%begc
    endc = bounds%endc

    allocate(this%topo_col(begc:endc))
    this%topo_col(:) = nan

    allocate(this%needs_downscaling_col(begc:endc))
    this%needs_downscaling_col(:) = .false.

  end subroutine InitAllocate

  !-----------------------------------------------------------------------
  subroutine InitCold(this, bounds)
    ! !USES:
    use column_varcon    , only: col_itype_to_icemec_class
    use clm_instur, only : topo_glc_mec
    ! !ARGUMENTS:
    class(topo_type), intent(inout) :: this
    type(bounds_type), intent(in) :: bounds
    !
    ! !LOCAL VARIABLES:
    integer :: c, l, g
    integer :: icemec_class            ! current icemec class (1..maxpatch_glcmec)

    character(len=*), parameter :: subname = 'InitCold'
    !-----------------------------------------------------------------------

    do c = bounds%begc, bounds%endc
       l = col%landunit(c)
       g = col%gridcell(c)

       if (lun%itype(l) == istice_mec) then
          ! For ice_mec landunits, initialize topo_col based on surface dataset; this
          ! will get overwritten in the run loop by values sent from CISM
          icemec_class = col_itype_to_icemec_class(col%itype(c))
          this%topo_col(c) = topo_glc_mec(g, icemec_class)
          this%needs_downscaling_col(c) = .true.
       else
          ! For other landunits, arbitrarily initialize topo_col to 0 m; for landunits
          ! where this matters, this will get overwritten in the run loop by values sent
          ! from CISM
          this%topo_col(c) = 0._r8
          this%needs_downscaling_col(c) = .false.
       end if
    end do

  end subroutine InitCold

  !-----------------------------------------------------------------------
  subroutine Restart(this, bounds, ncid, flag)
    ! !USES:
    use ncdio_pio, only : file_desc_t, ncd_double
    use restUtilMod
    !
    ! !ARGUMENTS:
    class(topo_type), intent(inout) :: this
    type(bounds_type), intent(in)    :: bounds 
    type(file_desc_t), intent(inout) :: ncid   ! netcdf id
    character(len=*) , intent(in)    :: flag   ! 'read', 'write' or 'define'
    !
    ! !LOCAL VARIABLES:
    integer :: p, c
    real(r8), pointer :: rparr(:)
    logical :: readvar

    character(len=*), parameter :: subname = 'Restart'
    !-----------------------------------------------------------------------

    allocate(rparr(bounds%begp:bounds%endp))

    call restartvar(ncid=ncid, flag=flag, varname='cols1d_topoglc', xtype=ncd_double,   &
         dim1name='column',                                                             &
         long_name='mean elevation on glacier elevation classes', units='m',            &
         interpinic_flag='skip', readvar=readvar, data=this%topo_col)

    if (flag /= 'read') then
       do p=bounds%begp,bounds%endp
          c = patch%column(p)
          rparr(p) = this%topo_col(c)
       enddo
       call restartvar(ncid=ncid, flag=flag, varname='pfts1d_topoglc', xtype=ncd_double,   &
            dim1name='pft',                                                             &
            long_name='mean elevation on glacier elevation classes', units='m',            &
            interpinic_flag='skip', readvar=readvar, data=rparr)
    end if

    deallocate(rparr)

  end subroutine Restart

  !-----------------------------------------------------------------------
  subroutine Clean(this)
    ! !ARGUMENTS:
    class(topo_type), intent(inout) :: this
    !
    ! !LOCAL VARIABLES:

    character(len=*), parameter :: subname = 'Clean'
    !-----------------------------------------------------------------------

    deallocate(this%topo_col)
    deallocate(this%needs_downscaling_col)

  end subroutine Clean

end module TopoMod
