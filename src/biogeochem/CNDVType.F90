module CNDVType

#include "shr_assert.h"

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Module containing routines to drive the annual dynamic vegetation
  ! that works with CN, reset related variables,
  ! and initialize/reset time invariant variables
  !
  ! !USES:
  use shr_kind_mod , only : r8 => shr_kind_r8
  use shr_log_mod  , only : errMsg => shr_log_errMsg
  use abortutils   , only : endrun
  use decompMod    , only : bounds_type
  use clm_varctl   , only : iulog
  !
  ! !PUBLIC TYPES:
  implicit none
  private
  !
  ! !PUBLIC DATA TYPES:
  !
  ! DGVM-specific ecophysiological constants structure (patch-level)
  type, public :: dgv_ecophyscon_type
     real(r8), pointer :: crownarea_max(:)   ! patch tree maximum crown area [m2]
     real(r8), pointer :: tcmin(:)           ! patch minimum coldest monthly mean temperature [units?]
     real(r8), pointer :: tcmax(:)           ! patch maximum coldest monthly mean temperature [units?]
     real(r8), pointer :: gddmin(:)          ! patch minimum growing degree days (at or above 5 C)
     real(r8), pointer :: twmax(:)           ! patch upper limit of temperature of the warmest month [units?]
     real(r8), pointer :: reinickerp(:)      ! patch parameter in allometric equation
     real(r8), pointer :: allom1(:)          ! patch parameter in allometric
     real(r8), pointer :: allom2(:)          ! patch parameter in allometric
     real(r8), pointer :: allom3(:)          ! patch parameter in allometric
  end type dgv_ecophyscon_type
  type(dgv_ecophyscon_type), public :: dgv_ecophyscon
  !
  ! DGVM state variables structure
  type, public :: dgvs_type
     real(r8), pointer, public :: agdd_patch        (:) ! patch accumulated growing degree days above 5
     real(r8), pointer, public :: agddtw_patch      (:) ! patch accumulated growing degree days above twmax
     real(r8), pointer, public :: agdd20_patch      (:) ! patch 20-yr running mean of agdd
     real(r8), pointer, public :: tmomin20_patch    (:) ! patch 20-yr running mean of tmomin
     logical , pointer, public :: present_patch     (:) ! patch whether PATCH present in patch
     logical , pointer, public :: pftmayexist_patch (:) ! patch if .false. then exclude seasonal decid patches from tropics
     real(r8), pointer, public :: nind_patch        (:) ! patch number of individuals (#/m**2)
     real(r8), pointer, public :: lm_ind_patch      (:) ! patch individual leaf mass
     real(r8), pointer, public :: lai_ind_patch     (:) ! patch LAI per individual
     real(r8), pointer, public :: fpcinc_patch      (:) ! patch foliar projective cover increment (fraction) 
     real(r8), pointer, public :: fpcgrid_patch     (:) ! patch foliar projective cover on gridcell (fraction)
     real(r8), pointer, public :: fpcgridold_patch  (:) ! patch last yr's fpcgrid
     real(r8), pointer, public :: crownarea_patch   (:) ! patch area that each individual tree takes up (m^2)
     real(r8), pointer, public :: greffic_patch     (:)
     real(r8), pointer, public :: heatstress_patch  (:)

   contains

     procedure , public  :: Init   
     procedure , public  :: Restart
     procedure , private :: InitAllocate 
  end type dgvs_type

  character(len=*), parameter, private :: sourcefile = &
       __FILE__
  !-----------------------------------------------------------------------

contains

  !------------------------------------------------------------------------
  subroutine Init(this, bounds)

    class(dgvs_type) :: this
    type(bounds_type), intent(in) :: bounds  

    ! Note - need allocation so that associate statements can be used
    ! at run time for NAG (allocation of variables is needed)

    call this%InitAllocate (bounds)

  end subroutine Init

  !-----------------------------------------------------------------------
  subroutine InitAllocate(this, bounds)
    !
    ! !USES:
    use shr_infnan_mod , only : nan => shr_infnan_nan, assignment(=)
    use clm_varpar     , only : numpft
    use pftconMod      , only : allom1s, allom2s, allom1, allom2, allom3, reinickerp
    use pftconMod      , only : ntree, nbrdlf_dcd_brl_shrub
    use pftconMod      , only : pftcon
    !
    ! !ARGUMENTS:
    class(dgvs_type) :: this
    type(bounds_type), intent(in) :: bounds  
    !
    ! !LOCAL VARIABLES:
    integer  :: begp, endp
    integer  :: m
    !---------------------------------------------------------------------

    begp = bounds%begp; endp = bounds%endp
  
    allocate(this%agdd_patch        (begp:endp)) ;     this%agdd_patch        (:) = nan
    allocate(this%agddtw_patch      (begp:endp)) ;     this%agddtw_patch      (:) = nan
    allocate(this%agdd20_patch      (begp:endp)) ;     this%agdd20_patch      (:) = nan
    allocate(this%tmomin20_patch    (begp:endp)) ;     this%tmomin20_patch    (:) = nan
    allocate(this%present_patch     (begp:endp)) ;     this%present_patch     (:) = .false.
    allocate(this%pftmayexist_patch (begp:endp)) ;     this%pftmayexist_patch (:) = .true.
    allocate(this%nind_patch        (begp:endp)) ;     this%nind_patch        (:) = nan
    allocate(this%lm_ind_patch      (begp:endp)) ;     this%lm_ind_patch      (:) = nan
    allocate(this%lai_ind_patch     (begp:endp)) ;     this%lai_ind_patch     (:) = nan
    allocate(this%fpcinc_patch      (begp:endp)) ;     this%fpcinc_patch      (:) = nan
    allocate(this%fpcgrid_patch     (begp:endp)) ;     this%fpcgrid_patch     (:) = nan
    allocate(this%fpcgridold_patch  (begp:endp)) ;     this%fpcgridold_patch  (:) = nan
    allocate(this%crownarea_patch   (begp:endp)) ;     this%crownarea_patch   (:) = nan
    allocate(this%greffic_patch     (begp:endp)) ;     this%greffic_patch     (:) = nan
    allocate(this%heatstress_patch  (begp:endp)) ;     this%heatstress_patch  (:) = nan

    allocate(dgv_ecophyscon%crownarea_max (0:numpft)) 
    allocate(dgv_ecophyscon%tcmin         (0:numpft))         
    allocate(dgv_ecophyscon%tcmax         (0:numpft))         
    allocate(dgv_ecophyscon%gddmin        (0:numpft))        
    allocate(dgv_ecophyscon%twmax         (0:numpft))         
    allocate(dgv_ecophyscon%reinickerp    (0:numpft))    
    allocate(dgv_ecophyscon%allom1        (0:numpft))        
    allocate(dgv_ecophyscon%allom2        (0:numpft))        
    allocate(dgv_ecophyscon%allom3        (0:numpft))        

    do m = 0,numpft
       dgv_ecophyscon%crownarea_max(m) = pftcon%pftpar20(m)
       dgv_ecophyscon%tcmin(m)         = pftcon%pftpar28(m)
       dgv_ecophyscon%tcmax(m)         = pftcon%pftpar29(m)
       dgv_ecophyscon%gddmin(m)        = pftcon%pftpar30(m)
       dgv_ecophyscon%twmax(m)         = pftcon%pftpar31(m)
       dgv_ecophyscon%reinickerp(m)    = reinickerp
       dgv_ecophyscon%allom1(m)        = allom1
       dgv_ecophyscon%allom2(m)        = allom2
       dgv_ecophyscon%allom3(m)        = allom3
       ! modification for shrubs by X.D.Z
       if (m > ntree .and. m <= nbrdlf_dcd_brl_shrub ) then 
          dgv_ecophyscon%allom1(m) = allom1s
          dgv_ecophyscon%allom2(m) = allom2s
       end if
    end do

  end subroutine InitAllocate

  !-----------------------------------------------------------------------
  subroutine Restart(this, bounds, ncid, flag)
    ! 
    ! !DESCRIPTION:
    ! Read/Write module information to/from restart file.
    !
    ! !USES:
    use clm_varcon , only : spval  
    use spmdMod    , only : masterproc
    use decompMod  , only : get_proc_global
    use restUtilMod
    use ncdio_pio
    use pio
    !
    ! !ARGUMENTS:
    class(dgvs_type) :: this
    type(bounds_type), intent(in)    :: bounds 
    type(file_desc_t), intent(inout) :: ncid   ! netcdf id
    character(len=*) , intent(in)    :: flag   ! 'read' or 'write'
    !
    ! !LOCAL VARIABLES:
    integer           :: j,c,p		! indices
    logical           :: readvar	! determine if variable is on initial file
    logical           :: do_io		! whether to do i/o for the given variable
    integer           :: nump_global	! total number of patches, globally
    integer           :: dimlen         ! dimension length
    integer           :: ier            ! error status
    integer           :: itemp		! temporary 
    integer , pointer :: iptemp(:)	! pointer to memory to be allocated
    integer           :: err_code       ! error code
    !-----------------------------------------------------------------------

    ! Get expected total number of points, for later error checks
    call get_proc_global(np=nump_global)

    call restartvar(ncid=ncid, flag=flag, varname='CROWNAREA', xtype=ncd_double,  &
         dim1name='pft', &
         long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%crownarea_patch)

    call restartvar(ncid=ncid, flag=flag, varname='nind', xtype=ncd_double,  &
         dim1name='pft', &
         long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%nind_patch)

    call restartvar(ncid=ncid, flag=flag, varname='fpcgrid', xtype=ncd_double,  &
         dim1name='pft', &
         long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%fpcgrid_patch)

    call restartvar(ncid=ncid, flag=flag, varname='fpcgridold', xtype=ncd_double,  &
         dim1name='pft', &
         long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%fpcgridold_patch)

    ! tmomin20
    do_io = .true.
    if (flag == 'read') then
       ! On a read, confirm that this variable has the expected size; if not, don't
       ! read it (instead leave it at its arbitrary initial value). This is needed to
       ! support older initial conditions for which this variable had a different size.
       call ncd_inqvdlen(ncid, 'TMOMIN20', 1, dimlen, err_code)
       if (dimlen /= nump_global) then
          do_io = .false.
       end if
    end if
    if (do_io) then
       call restartvar(ncid=ncid, flag=flag, varname='TMOMIN20', xtype=ncd_double,  &
            dim1name='pft', &
            long_name='',units='', &
            interpinic_flag='interp', readvar=readvar, data=this%tmomin20_patch)
    end if

    ! agdd20
    do_io = .true.
    if (flag == 'read') then
       ! On a read, confirm that this variable has the expected size; if not, don't
       ! read it (instead leave it at its arbitrary initial value). This is needed to
       ! support older initial conditions for which this variable had a different size.
       call ncd_inqvdlen(ncid, 'AGDD20', 1, dimlen, err_code)
       if (dimlen /= nump_global) then
          do_io = .false.
       end if
    end if
    if (do_io) then
       call restartvar(ncid=ncid, flag=flag, varname='AGDD20', xtype=ncd_double,  &
            dim1name='pft',&
            long_name='',units='', &
            interpinic_flag='interp', readvar=readvar, data=this%agdd20_patch)
    end if

    ! present  
    if (flag == 'read' .or. flag == 'write') then
       allocate (iptemp(bounds%begp:bounds%endp), stat=ier)
    end if
    if (flag == 'write') then
       do p = bounds%begp,bounds%endp
          iptemp(p) = 0
          if (this%present_patch(p)) iptemp(p) = 1
       end do
    end if
    call restartvar(ncid=ncid, flag=flag, varname='present', xtype=ncd_int,  &
         dim1name='pft',&
         long_name='',units='', &
         interpinic_flag='interp', readvar=readvar, data=iptemp)
    if (flag=='read' .and. readvar) then
       do p = bounds%begp,bounds%endp
          this%present_patch(p) = .false.
          if (iptemp(p) == 1) this%present_patch(p) = .true.
       end do
    end if
    if (flag == 'read' .or. flag == 'write') then
       deallocate (iptemp)
    end if

    call restartvar(ncid=ncid, flag=flag, varname='heatstress', xtype=ncd_double,  &
         dim1name='pft', &
         long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%heatstress_patch)

    call restartvar(ncid=ncid, flag=flag, varname='greffic', xtype=ncd_double,  &
         dim1name='pft', &
         long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%greffic_patch)

  end subroutine Restart

end module CNDVType
