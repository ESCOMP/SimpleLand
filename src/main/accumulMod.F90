module accumulMod

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! This module contains generic subroutines that can be used to
  ! define, accumulate and extract  user-specified fields over
  ! user-defined intervals. Each interval  and accumulation type is
  ! unique to each field processed.
  ! Subroutine [init_accumulator] defines the values of the accumulated
  ! field data structure. Subroutine [update_accum_field] does
  ! the actual accumulation for a given field.
  ! Three types of accumulations are possible:
  ! - Average over time interval. Time average fields are only
  !   valid at the end of the averaging interval.
  ! - Running mean over time interval. Running means are valid once the
  !   length of the simulation exceeds the
  ! - Running accumulation over time interval. Accumulated fields are
  !   continuously accumulated. The trigger value "-99999." resets
  !   the accumulation to zero.
  !
  ! !USES:
#include "shr_assert.h"
  use shr_kind_mod, only: r8 => shr_kind_r8
  use shr_sys_mod , only: shr_sys_abort
  use shr_log_mod , only: errMsg => shr_log_errMsg
  use abortutils  , only: endrun
  use clm_varctl  , only: iulog, nsrest, nsrStartup
  use clm_varcon  , only: spval, ispval
  use GridcellType, only : grc
  !
  ! !PUBLIC TYPES:
  implicit none
  save
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: accumulRest          ! Write/read restart of accumulation fields
  public :: init_accum_field     ! Initialize an accumulator field
  public :: print_accum_fields   ! Print info about accumulator fields
  public :: extract_accum_field  ! Extracts the current value of an accumulator field
  public :: update_accum_field   ! Update the current value of an accumulator field
  public :: clean_accum_fields   ! Deallocate space and reset accum fields list

  interface extract_accum_field
     module procedure extract_accum_field_sl ! Extract current val of single-level accumulator field
     module procedure extract_accum_field_ml ! Extract current val of multi-level accumulator field
  end interface
  interface update_accum_field               ! Updates the current value of an accumulator field
     module procedure update_accum_field_sl  ! Update single-level accumulator field
     module procedure update_accum_field_ml  ! Update multi-level accumulator field
  end interface
  private
  !
  ! !PRIVATE MEMBER FUNCTIONS:
  private :: extract_accum_field_basic   ! Extract values for one level of the given field
  private :: extract_accum_field_timeavg ! Extract values for one level of the given timeavg field
  private :: update_accum_field_timeavg  ! Update values for one level of the given timeavg field
  private :: update_accum_field_runmean  ! Update values for one level of the given runmean field
  private :: update_accum_field_runaccum ! Update values for one level of the given runaccum field
  private :: find_field                  ! Find field index given field name
  private :: acctype_to_string  ! Return a string representation of an ACCTYPE parameter
  !
  ! !PRIVATE TYPES:

  type accum_field
     character(len=  8) :: name     !field name
     character(len=128) :: desc     !field description
     character(len=  8) :: units    !field units
     integer            :: acctype  !accumulation type (see ACCTYPE parameters below)
     character(len=  8) :: type1d   !subgrid type: ["gridcell","landunit","column" or "pft"]
     character(len=  8) :: type2d   !type2d ('','levsoi','numrad',..etc. )
     integer            :: beg1d    !subgrid type beginning index
     integer            :: end1d    !subgrid type ending index
     integer            :: numlev   !number of vertical levels in field
     logical, pointer   :: active(:)!whether each point (patch, column, etc.) is active
     real(r8)           :: initval  !initial value of accumulated field
     real(r8), pointer  :: val(:,:) !accumulated field
     integer            :: period   !field accumulation period (in model time steps)

     ! In most cases, we could use a 1-d nsteps variable. However, that's awkward within
     ! nested loops (with level as the outer loop); also, runaccum can theoretically have
     ! different reset points for different levels.
     integer, pointer   :: nsteps(:,:)!number of steps each point has accumulated, since last reset time

     ! NOTE(wjs, 2017-12-03) We should convert this to fully object-oriented (with
     ! inheritance / polymorphism). For now, in the interest of time, I'm going with a
     ! semi-object-oriented solution of using procedure pointers.
     procedure(extract_accum_field_interface), pointer :: extract_accum_field_func
     procedure(update_accum_field_interface) , pointer :: update_accum_field_func
  end type accum_field

  abstract interface
     subroutine extract_accum_field_interface(this, level, nstep, field)
       use shr_kind_mod, only: r8 => shr_kind_r8
       import :: accum_field
       class(accum_field), intent(in) :: this
       integer, intent(in) :: level      ! level index to extract (1 for a 1-d field)
       integer, intent(in) :: nstep      ! timestep index
       real(r8), intent(inout) :: field(:) ! field values for current time step
     end subroutine extract_accum_field_interface

     subroutine update_accum_field_interface(this, level, nstep, field)
       use shr_kind_mod, only: r8 => shr_kind_r8
       import :: accum_field
       class(accum_field), intent(in) :: this
       integer, intent(in) :: level      ! level index to update (1 for a 1-d field)
       integer, intent(in) :: nstep      ! timestep index
       real(r8), intent(in) :: field(:)  ! field values for current time step
     end subroutine update_accum_field_interface
  end interface

  real(r8), parameter, public :: accumResetVal = -99999._r8 ! used to do an annual reset ( put in for bug 1858)

  integer, parameter :: ACCTYPE_TIMEAVG  = 1
  integer, parameter :: ACCTYPE_RUNMEAN  = 2
  integer, parameter :: ACCTYPE_RUNACCUM = 3

  integer, parameter :: max_accum = 100    !maximum number of accumulated fields
  type (accum_field) :: accum(max_accum)   !array accumulated fields
  integer :: naccflds = 0                  !accumulator field counter

  character(len=*), parameter, private :: sourcefile = &
       __FILE__

  !------------------------------------------------------------------------

contains

  !------------------------------------------------------------------------
  subroutine init_accum_field (name, units, desc, &
       accum_type, accum_period, numlev, subgrid_type, init_value, type2d)
    !
    ! !DESCRIPTION:
    ! Initialize accumulation fields. This subroutine sets:
    ! o name  of accumulated field
    ! o units of accumulated field
    ! o accumulation type of accumulated field
    ! o description of accumulated fields: accdes
    ! o accumulation period for accumulated field (in iterations)
    ! o initial value of accumulated field
    !
    ! Note about initial value: This must be 0 for a timeavg or runaccum field. For a
    ! runmean field, initializing to a non-zero value often won't accomplish anything, but
    ! can be used to start with a reasonable value in situations such as: (1) the field is
    ! extracted before the first update call; (2) edge case situations such as doing a
    ! branch run from a restart file that did not contain this field (though it's
    ! possible that init_value doesn't matter even in this case).
    !
    ! !USES:
    use shr_const_mod, only: SHR_CONST_CDAY
    use clm_time_manager, only : get_step_size
    use decompMod, only : get_proc_bounds
    !
    ! !ARGUMENTS:
    implicit none
    character(len=*), intent(in)           :: name         !field name
    character(len=*), intent(in)           :: units        !field units
    character(len=*), intent(in)           :: desc         !field description
    character(len=*), intent(in)           :: accum_type   !field type: timeavg, runmean, runaccum
    integer , intent(in)                   :: accum_period !field accumulation period
    character(len=*), intent(in)           :: subgrid_type !["gridcell"]
    integer , intent(in)                   :: numlev       !number of vertical levels
    real(r8), intent(in)                   :: init_value   !field initial or reset value
    character(len=*), intent(in), optional :: type2d       !level type (optional) - needed if numlev > 1
    !
    ! !LOCAL VARIABLES:
    integer :: nf           ! field index
    integer :: beg1d,end1d  ! beggining and end subgrid indices
    integer :: begg, endg   ! per-proc gridcell ending gridcell indices
    !------------------------------------------------------------------------

    ! Determine necessary indices

    call get_proc_bounds(begg, endg)

    ! update field index
    ! Consistency check that number of accumulated does not exceed maximum.

    naccflds = naccflds + 1
    if (naccflds > max_accum) then
       write(iulog,*) 'ACCUMULINIT error: user-defined accumulation fields ', &
            'equal to ',naccflds,' exceeds max_accum'
       call shr_sys_abort()
    end if
    nf = naccflds

    select case (trim(accum_type))
    case ('timeavg')
       accum(nf)%acctype = ACCTYPE_TIMEAVG
       accum(nf)%extract_accum_field_func => extract_accum_field_timeavg
       accum(nf)%update_accum_field_func  => update_accum_field_timeavg
    case ('runmean')
       accum(nf)%acctype = ACCTYPE_RUNMEAN
       accum(nf)%extract_accum_field_func => extract_accum_field_basic
       accum(nf)%update_accum_field_func  => update_accum_field_runmean
    case ('runaccum')
       accum(nf)%acctype = ACCTYPE_RUNACCUM
       accum(nf)%extract_accum_field_func => extract_accum_field_basic
       accum(nf)%update_accum_field_func  => update_accum_field_runaccum
    case default
       write(iulog,*) 'init_accum_field ERROR: unknown accum_type ', accum_type
       call shr_sys_abort('init_accum_field: unknown accum_type')
    end select

    accum(nf)%name    = trim(name)
    accum(nf)%units   = trim(units)
    accum(nf)%desc    = trim(desc)
    accum(nf)%initval = init_value

    ! Note accumulation period must be converted from days
    ! to number of iterations
    accum(nf)%period  = accum_period
    if (accum(nf)%period < 0) then
       accum(nf)%period = -accum(nf)%period * nint(SHR_CONST_CDAY) / get_step_size()
    end if

    select case (trim(subgrid_type))
    case ('gridcell')
       beg1d = begg
       end1d = endg
       accum(nf)%active => grc%active
    case default
       write(iulog,*)'init_accum_field: unknown subgrid type ',subgrid_type
       call shr_sys_abort ()
    end select

    accum(nf)%type1d = trim(subgrid_type)
    accum(nf)%beg1d = beg1d
    accum(nf)%end1d = end1d
    accum(nf)%numlev = numlev

    if (present(type2d)) then
       accum(nf)%type2d = type2d
    else
       accum(nf)%type2d = ' '
    end if
    
    ! Allocate and initialize accumulation field

    allocate(accum(nf)%val(beg1d:end1d,numlev))
    if (accum(nf)%acctype == ACCTYPE_TIMEAVG .or. &
         accum(nf)%acctype == ACCTYPE_RUNACCUM) then
       if (init_value /= 0._r8) then
          write(iulog,*) 'init_accum_field ERROR: for field ', trim(name)
          write(iulog,*) 'init_value must be 0 for timeavg and runaccum fields'
          call shr_sys_abort('init_accum_field: init_value must be 0 for timeavg and runaccum fields')
       end if
    end if
    accum(nf)%val(beg1d:end1d,1:numlev) = init_value

    allocate(accum(nf)%nsteps(beg1d:end1d,numlev))
    accum(nf)%nsteps(beg1d:end1d,1:numlev) = 0

  end subroutine init_accum_field

  !------------------------------------------------------------------------
  subroutine print_accum_fields()
    !
    ! !DESCRIPTION:
    ! Diagnostic printout of accumulated fields
    !
    ! !USES:
    use spmdMod, only : masterproc
    !
    ! !ARGUMENTS:
    implicit none
    !
    integer :: i,nf   !indices
    !------------------------------------------------------------------------

    if (masterproc) then
       write(iulog,*)
       write(iulog,*) 'Initializing variables for time accumulation .....'
       write(iulog,'(72a1)') ("-",i=1,60)
       write(iulog,*) 'Accumulated fields'
       write(iulog,1002)
       write(iulog,'(72a1)') ("_",i=1,71)
       do nf = 1, naccflds
          if (accum(nf)%period /= huge(1)) then
             write(iulog,1003) nf,accum(nf)%name,accum(nf)%units,&
                  acctype_to_string(accum(nf)%acctype), &
                  accum(nf)%period, accum(nf)%initval, &
                  accum(nf)%desc
          else
             write(iulog,1004) nf,accum(nf)%name,accum(nf)%units,&
                  acctype_to_string(accum(nf)%acctype), &
                  accum(nf)%initval, accum(nf)%desc
          endif
       end do
       write(iulog,'(72a1)') ("_",i=1,71)
       write(iulog,*)
       write(iulog,'(72a1)') ("-",i=1,60)
       write(iulog,*) 'Successfully initialized variables for accumulation'
       write(iulog,*)
    endif

1002 format(' No',' Name    ',' Units   ',' Type    ','Period',' Inival',' Description')
1003 format((1x,i2),(1x,a8),(1x,a8),(1x,a8), (1x,i5),(1x,f4.0),(1x,a40))
1004 format((1x,i2),(1x,a8),(1x,a8),(1x,a8),'  N.A.',(1x,f4.0),(1x,a40))

  end subroutine print_accum_fields

  !------------------------------------------------------------------------
  subroutine extract_accum_field_sl (name, field, nstep)
    !
    ! !DESCRIPTION:
    ! Extract single-level accumulated field.
    ! This routine extracts the field values from the multi-level
    ! accumulation field. It extracts the current value except if
    ! the field type is a time average. In this case, an absurd value
    ! is assigned to  indicate the time average is not yet valid.
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    implicit none
    character(len=*), intent(in) :: name     !field name
    real(r8), pointer, dimension(:) :: field !field values for current time step
    integer , intent(in) :: nstep            !timestep index
    !
    ! !LOCAL VARIABLES:
    integer :: nf     ! field index
    integer :: numlev ! number of vertical levels

    character(len=*), parameter :: subname = 'extract_accum_field_sl'
    !------------------------------------------------------------------------

    call find_field(field_name=name, caller_name=subname, field_index=nf)

    numlev = accum(nf)%numlev
    SHR_ASSERT(numlev == 1, errMsg(sourcefile, __LINE__))

    call accum(nf)%extract_accum_field_func( &
         level = 1, &
         nstep = nstep, &
         field = field)

  end subroutine extract_accum_field_sl

  !------------------------------------------------------------------------
  subroutine extract_accum_field_ml (name, field, nstep)
    !
    ! !DESCRIPTION:
    ! Extract mutli-level accumulated field.
    ! This routine extracts the field values from the multi-level
    ! accumulation field. It extracts the current value except if
    ! the field type is a time average. In this case, an absurd value
    ! is assigned to  indicate the time average is not yet valid.
    !
    ! !ARGUMENTS:
    implicit none
    character(len=*), intent(in) :: name       !field name
    real(r8), pointer, dimension(:,:) :: field !field values for current time step
    integer, intent(in) :: nstep               !timestep index
    !
    ! !LOCAL VARIABLES:
    integer :: j,nf            !indices
    integer :: numlev          !number of vertical levels

    character(len=*), parameter :: subname = 'extract_accum_field_ml'
    !------------------------------------------------------------------------

    call find_field(field_name=name, caller_name=subname, field_index=nf)

    numlev = accum(nf)%numlev
    SHR_ASSERT((size(field, 2) == numlev), errMsg(sourcefile, __LINE__))

    do j = 1,numlev
       call accum(nf)%extract_accum_field_func( &
            level = j, &
            nstep = nstep, &
            field = field(:,j))
    end do

  end subroutine extract_accum_field_ml

  !-----------------------------------------------------------------------
  subroutine extract_accum_field_basic(this, level, nstep, field)
    ! !DESCRIPTION:
    ! Extract values for one level of the given field
    !
    ! !ARGUMENTS:
    class(accum_field), intent(in) :: this
    integer, intent(in) :: level      ! level index to extract (1 for a 1-d field)
    integer, intent(in) :: nstep      ! timestep index
    real(r8), intent(inout) :: field(:) ! field values for current time step
    !
    ! !LOCAL VARIABLES:
    integer :: begi,endi         !subgrid beginning,ending indices
    integer :: k, kf

    character(len=*), parameter :: subname = 'extract_accum_field_basic'
    !-----------------------------------------------------------------------

    begi = this%beg1d
    endi = this%end1d
    SHR_ASSERT((size(field) == endi-begi+1), errMsg(sourcefile, __LINE__))

    do k = begi, endi
       kf = k - begi + 1
       field(kf) = this%val(k,level)
    end do

  end subroutine extract_accum_field_basic

  !-----------------------------------------------------------------------
  subroutine extract_accum_field_timeavg(this, level, nstep, field)
    ! !DESCRIPTION:
    ! Extract values for one level of the given timeavg field
    !
    ! !ARGUMENTS:
    class(accum_field), intent(in) :: this
    integer, intent(in) :: level      ! level index to extract (1 for a 1-d field)
    integer, intent(in) :: nstep      ! timestep index
    real(r8), intent(inout) :: field(:) ! field values for current time step
    !
    ! !LOCAL VARIABLES:
    integer :: begi,endi         !subgrid beginning,ending indices
    integer :: k, kf

    character(len=*), parameter :: subname = 'extract_accum_field_basic'
    !-----------------------------------------------------------------------

    begi = this%beg1d
    endi = this%end1d
    SHR_ASSERT((size(field) == endi-begi+1), errMsg(sourcefile, __LINE__))

    if (mod(nstep,this%period) == 0) then
       do k = begi, endi
          kf = k - begi + 1
          field(kf) = this%val(k,level)
       end do
    else
       do k = begi, endi
          kf = k - begi + 1
          field(kf) = spval
       end do
    end if

  end subroutine extract_accum_field_timeavg

  !------------------------------------------------------------------------
  subroutine update_accum_field_sl (name, field, nstep)
    !
    ! !DESCRIPTION:
    ! Accumulate single level field over specified time interval.
    ! The appropriate field is accumulated in the array [accval].
    !
    ! Values of 'field' are ignored in inactive points, so it's safe for 'field' to
    ! contain garbage in inactive points.
    !
    ! !ARGUMENTS:
    implicit none
    character(len=*), intent(in) :: name     !field name
    real(r8), pointer, dimension(:) :: field !field values for current time step
    integer , intent(in) :: nstep            !time step index
    !
    ! !LOCAL VARIABLES:
    integer :: nf     ! field index
    integer :: numlev ! number of vertical levels

    character(len=*), parameter :: subname = 'update_accum_field_sl'
    !------------------------------------------------------------------------

    call find_field(field_name=name, caller_name=subname, field_index=nf)

    numlev = accum(nf)%numlev
    SHR_ASSERT(numlev == 1, errMsg(sourcefile, __LINE__))

    call accum(nf)%update_accum_field_func( &
         level = 1, &
         nstep = nstep, &
         field = field)

  end subroutine update_accum_field_sl

  !------------------------------------------------------------------------
  subroutine update_accum_field_ml (name, field, nstep)
    !
    ! !DESCRIPTION:
    ! Accumulate multi level field over specified time interval.
    !
    ! Values of 'field' are ignored in inactive points, so it's safe for 'field' to
    ! contain garbage in inactive points.
    !
    ! !ARGUMENTS:
    implicit none
    character(len=*), intent(in) :: name       !field name
    real(r8), pointer, dimension(:,:) :: field !field values for current time step
    integer , intent(in) :: nstep              !time step index
    !
    ! !LOCAL VARIABLES:
    integer :: j,nf            !indices
    integer :: numlev          !number of vertical levels

    character(len=*), parameter :: subname = 'update_accum_field_ml'
    !------------------------------------------------------------------------

    call find_field(field_name=name, caller_name=subname, field_index=nf)

    numlev = accum(nf)%numlev
    SHR_ASSERT((size(field, 2) == numlev), errMsg(sourcefile, __LINE__))

    do j = 1,numlev
       call accum(nf)%update_accum_field_func( &
            level = j, &
            nstep = nstep, &
            field = field(:,j))
    end do

  end subroutine update_accum_field_ml

  !-----------------------------------------------------------------------
  subroutine update_accum_field_timeavg(this, level, nstep, field)
    !
    ! !DESCRIPTION:
    ! Update values for one level of the given timeavg field
    !
    ! !ARGUMENTS:
    class(accum_field), intent(in) :: this
    integer, intent(in) :: level      ! level index to update (1 for a 1-d field)
    integer, intent(in) :: nstep      ! timestep index
    real(r8), intent(in) :: field(:)  ! field values for current time step
    !
    ! !LOCAL VARIABLES:
    integer :: begi,endi         !subgrid beginning,ending indices
    integer :: k, kf

    character(len=*), parameter :: subname = 'update_accum_field_timeavg'
    !-----------------------------------------------------------------------

    begi = this%beg1d
    endi = this%end1d
    SHR_ASSERT((size(field) == endi-begi+1), errMsg(sourcefile, __LINE__))

    ! time average field: reset every accumulation period; normalize at end of
    ! accumulation period

    if ((mod(nstep,this%period) == 1 .or. this%period == 1) .and. (nstep /= 0))then
       do k = begi,endi
          if (this%active(k)) then
             this%val(k,level) = 0._r8
             this%nsteps(k,level) = 0
          end if
       end do
    end if

    do k = begi,endi
       if (this%active(k)) then
          kf = k - begi + 1
          this%val(k,level) =  this%val(k,level) + field(kf)
          this%nsteps(k,level) = this%nsteps(k,level) + 1
       end if
    end do

    if (mod(nstep,this%period) == 0) then
       do k = begi,endi
          if (this%active(k)) then
             this%val(k,level) = this%val(k,level) / this%nsteps(k,level)
          end if
       end do
    end if

  end subroutine update_accum_field_timeavg

  !-----------------------------------------------------------------------
  subroutine update_accum_field_runmean(this, level, nstep, field)
    !
    ! !DESCRIPTION:
    ! Update values for one level of the given runmean field
    !
    ! !ARGUMENTS:
    class(accum_field), intent(in) :: this
    integer, intent(in) :: level      ! level index to update (1 for a 1-d field)
    integer, intent(in) :: nstep      ! timestep index
    real(r8), intent(in) :: field(:)  ! field values for current time step
    !
    ! !LOCAL VARIABLES:
    integer :: begi,endi         !subgrid beginning,ending indices
    integer :: k, kf
    integer :: accper  ! accumulation period

    character(len=*), parameter :: subname = 'update_accum_field_runmean'
    !-----------------------------------------------------------------------

    begi = this%beg1d
    endi = this%end1d
    SHR_ASSERT((size(field) == endi-begi+1), errMsg(sourcefile, __LINE__))

    do k = begi,endi
       if (this%active(k)) then
          kf = k - begi + 1
          this%nsteps(k,level) = this%nsteps(k,level) + 1
          ! Cap nsteps at 'period' - partly to avoid overflow, but also because it doesn't
          ! help us to accumulate nsteps beyond a value of 'period' (because nsteps is
          ! just used to determine accper, and accper needs to be capped at 'period'). A
          ! side-benefit of this capping of nsteps is that accper (the accumulation
          ! period) is always equal to nsteps.
          this%nsteps(k,level) = min(this%nsteps(k,level), this%period)
          accper = this%nsteps(k,level)
          this%val(k,level) = &
               ((accper-1)*this%val(k,level) + field(kf)) / accper
       end if
    end do
       
  end subroutine update_accum_field_runmean

  !-----------------------------------------------------------------------
  subroutine update_accum_field_runaccum(this, level, nstep, field)
    !
    ! !DESCRIPTION:
    ! Update values for one level of the given runaccum field
    !
    ! !ARGUMENTS:
    class(accum_field), intent(in) :: this
    integer, intent(in) :: level      ! level index to update (1 for a 1-d field)
    integer, intent(in) :: nstep      ! timestep index
    real(r8), intent(in) :: field(:)  ! field values for current time step
    !
    ! !LOCAL VARIABLES:
    integer :: begi,endi         !subgrid beginning,ending indices
    integer :: k, kf

    character(len=*), parameter :: subname = 'update_accum_field_runaccum'
    !-----------------------------------------------------------------------

    begi = this%beg1d
    endi = this%end1d
    SHR_ASSERT((size(field) == endi-begi+1), errMsg(sourcefile, __LINE__))

    !running accumulation field; reset at trigger -99999

    do k = begi,endi
       if (this%active(k)) then
          kf = k - begi + 1
          if (nint(field(kf)) == -99999) then
             this%val(k,level) = 0._r8
             this%nsteps(k,level) = 0
          else
             this%val(k,level) = &
                  min(max(this%val(k,level) + field(kf), 0._r8), 99999._r8)
             this%nsteps(k,level) = this%nsteps(k,level) + 1
          end if
       end if
    end do

  end subroutine update_accum_field_runaccum


  !------------------------------------------------------------------------
  subroutine accumulRest( ncid, flag )
    !
    ! !DESCRIPTION:
    ! Read/write accumulation restart data
    !
    ! !USES:
    use restUtilMod     , only : restartvar
    use ncdio_pio       , only : file_desc_t, ncd_double, ncd_int
    !
    ! !ARGUMENTS:
    implicit none
    type(file_desc_t), intent(inout) :: ncid   !netcdf unit
    character(len=*) , intent(in) :: flag   !'define','read', or 'write'
    !
    ! !LOCAL VARIABLES:
    integer :: nf                             ! indices
    logical :: readvar                        ! determine if variable is on initial file
    character(len=128) :: varname             ! temporary
    character(len= 32) :: subname='AccumRest' ! subroutine name
    !------------------------------------------------------------------------

    do nf = 1,naccflds

       varname = trim(accum(nf)%name) // '_VALUE'
       if (accum(nf)%numlev == 1) then
          call restartvar(ncid=ncid, flag=flag, varname=varname, xtype=ncd_double, &
               dim1name=accum(nf)%type1d, &
               long_name=accum(nf)%desc, units=accum(nf)%units, &
               interpinic_flag='interp', &
               data=accum(nf)%val, readvar=readvar)
       else
          call restartvar(ncid=ncid, flag=flag, varname=varname, xtype=ncd_double, &
               dim1name=accum(nf)%type1d, dim2name=accum(nf)%type2d, &
               long_name=accum(nf)%desc, units=accum(nf)%units, &
               interpinic_flag='interp', &
               data=accum(nf)%val, readvar=readvar)
       end if

       varname = trim(accum(nf)%name) // '_NSTEPS'
       if (accum(nf)%numlev == 1) then
          call restartvar(ncid=ncid, flag=flag, varname=varname, xtype=ncd_int, &
               dim1name=accum(nf)%type1d, &
               long_name='number of accumulated steps for '//trim(accum(nf)%name), &
               units='-', &
               interpinic_flag='interp', &
               data=accum(nf)%nsteps, readvar=readvar)
       else
          call restartvar(ncid=ncid, flag=flag, varname=varname, xtype=ncd_int, &
               dim1name=accum(nf)%type1d, dim2name=accum(nf)%type2d, &
               long_name='number of accumulated steps for '//trim(accum(nf)%name), &
               units='-', &
               interpinic_flag='interp', &
               data=accum(nf)%nsteps, readvar=readvar)
       end if

       varname = trim(accum(nf)%name) // '_PERIOD'
       call restartvar(ncid=ncid, flag=flag, varname=varname, xtype=ncd_int, &
            long_name='', units='time steps', &
            imissing_value=ispval, ifill_value=huge(1), &
            interpinic_flag='copy', &
            data=accum(nf)%period, readvar=readvar)

    end do

  end subroutine accumulRest

  !-----------------------------------------------------------------------
  subroutine clean_accum_fields
    !
    ! !DESCRIPTION:
    ! Deallocate space and reset accum fields list
    !
    ! !ARGUMENTS:
    !
    ! !LOCAL VARIABLES:
    integer :: i

    character(len=*), parameter :: subname = 'clean_accum_fields'
    !-----------------------------------------------------------------------

    do i = 1, naccflds
       if (associated(accum(i)%val)) then
          deallocate(accum(i)%val)
       end if
       if (associated(accum(i)%nsteps)) then
          deallocate(accum(i)%nsteps)
       end if
    end do

    naccflds = 0

  end subroutine clean_accum_fields


  !-----------------------------------------------------------------------
  subroutine find_field(field_name, caller_name, field_index)
    !
    ! !DESCRIPTION:
    ! Find field index given field name
    !
    ! Aborts if the given field name isn't found
    !
    ! !ARGUMENTS:
    character(len=*), intent(in)  :: field_name  ! field name to find
    character(len=*), intent(in)  :: caller_name ! name of calling routine (for more informative error messages)
    integer         , intent(out) :: field_index ! index of given field in accum array
    !
    ! !LOCAL VARIABLES:
    integer :: i

    character(len=*), parameter :: subname = 'find_field'
    !-----------------------------------------------------------------------

    field_index = 0
    do i = 1, naccflds
       if (field_name == accum(i)%name) then
          field_index = i
          exit
       end if
    end do
    if (field_index == 0) then
       write(iulog,*) trim(caller_name), 'ERROR: field name ',trim(field_name),' not found'
       call endrun('Accumulation field not found')
    end if

  end subroutine find_field


  !-----------------------------------------------------------------------
  pure function acctype_to_string(acctype) result(acctype_str)
    !
    ! !DESCRIPTION:
    ! Return a string representation of an ACCTYPE parameter
    !
    ! !ARGUMENTS:
    character(len=32) :: acctype_str  ! function result
    integer, intent(in) :: acctype
    !
    ! !LOCAL VARIABLES:

    character(len=*), parameter :: subname = 'acctype_to_string'
    !-----------------------------------------------------------------------

    select case (acctype)
    case (ACCTYPE_TIMEAVG)
       acctype_str = 'timeavg'
    case (ACCTYPE_RUNMEAN)
       acctype_str = 'runmean'
    case (ACCTYPE_RUNACCUM)
       acctype_str = 'runaccum'
    case default
       acctype_str = 'unknown'
    end select

  end function acctype_to_string


end module accumulMod
