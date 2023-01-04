module clm_instMod

  !------------------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Instances and definitions of all data types
  !
  ! !USES:
  use shr_kind_mod    , only : r8 => shr_kind_r8
  use decompMod       , only : bounds_type
  use clm_varcon      , only : bdsno
  use clm_varctl      , only : iulog
  use perf_mod        , only : t_startf, t_stopf

  !-----------------------------------------
  ! Constants
  !-----------------------------------------

  !-----------------------------------------
  ! Definition of component types 
  !-----------------------------------------

  use atm2lndType                     , only : atm2lnd_type
  use lnd2atmType                     , only : lnd2atm_type
  use GridcellType                    , only : grc
  !
  implicit none
  public   ! By default everything is public 
  !
  !-----------------------------------------
  ! Instances of component types
  !-----------------------------------------

  ! Physics types 
  type(atm2lnd_type)                      :: atm2lnd_inst
  type(lnd2atm_type)                      :: lnd2atm_inst

  public :: clm_instInit       ! Initialize
  public :: clm_instRest       ! Setup restart
  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine clm_instInit(bounds)
    !
    ! !USES: 
    use clm_varpar                         , only : nlevsno
    !
    ! !ARGUMENTS    
    type(bounds_type), intent(in) :: bounds  ! processor bounds
    !
    ! !LOCAL VARIABLES:
    !----------------------------------------------------------------------

    ! Initialize clm->drv and drv->clm data structures

    call atm2lnd_inst%Init(bounds)
    call lnd2atm_inst%Init(bounds)

    ! Initialization of public data types

! TODO SLIM: slevis keeping an example of an accumulated field as template
!   ! ------------------------------------------------------------------------
!   ! Initialize accumulated fields
!   ! ------------------------------------------------------------------------

!   ! The time manager needs to be initialized before this called is made, since
!   ! the step size is needed. 

!   call t_startf('init_accflds')
!   call atm2lnd_inst%InitAccBuffer(bounds)
!   call t_stopf('init_accflds')

  end subroutine clm_instInit

  !-----------------------------------------------------------------------
  subroutine clm_instRest(bounds, ncid, flag)
    !
    ! !USES:
    use ncdio_pio       , only : file_desc_t
    !
    ! !DESCRIPTION:
    ! Define/write/read CLM restart file.
    !
    ! !ARGUMENTS:
    type(bounds_type) , intent(in)    :: bounds          
    type(file_desc_t) , intent(inout) :: ncid ! netcdf id
    character(len=*)  , intent(in)    :: flag ! 'define', 'write', 'read' 
    !-----------------------------------------------------------------------

    call atm2lnd_inst%restart (bounds, ncid, flag=flag)

 end subroutine clm_instRest

end module clm_instMod

