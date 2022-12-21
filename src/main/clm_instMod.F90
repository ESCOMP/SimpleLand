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
  use landunit_varcon , only : istice_mec, istsoil
  use perf_mod        , only : t_startf, t_stopf

  !-----------------------------------------
  ! Constants
  !-----------------------------------------

  use UrbanParamsType                    , only : urbanparams_type   ! Constants 
  !-----------------------------------------
  ! Definition of component types 
  !-----------------------------------------

  use UrbanParamsType                 , only : urbanparams_type
  use atm2lndType                     , only : atm2lnd_type
  use lnd2atmType                     , only : lnd2atm_type
  use glcBehaviorMod                  , only : glc_behavior_type
  use TopoMod                         , only : topo_type
  use GridcellType                    , only : grc
  use LandunitType                    , only : lun                
  use ColumnType                      , only : col                
  use PatchType                       , only : patch                
  !
  implicit none
  public   ! By default everything is public 
  !
  !-----------------------------------------
  ! Instances of component types
  !-----------------------------------------

  ! Physics types 
  type(urbanparams_type)                  :: urbanparams_inst
  type(atm2lnd_type)                      :: atm2lnd_inst
  type(lnd2atm_type)                      :: lnd2atm_inst
  type(glc_behavior_type), target         :: glc_behavior
  type(topo_type)                         :: topo_inst

  public :: clm_instInit       ! Initialize
  public :: clm_instRest       ! Setup restart
  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine clm_instInit(bounds)
    !
    ! !USES: 
    use clm_varpar                         , only : nlevsno
    use initVerticalMod                    , only : initVertical
    !
    ! !ARGUMENTS    
    type(bounds_type), intent(in) :: bounds  ! processor bounds
    !
    ! !LOCAL VARIABLES:
    integer               :: c,l,g
    integer               :: begp, endp
    integer               :: begc, endc
    integer               :: begl, endl
    !----------------------------------------------------------------------

    begp = bounds%begp; endp = bounds%endp 
    begc = bounds%begc; endc = bounds%endc 
    begl = bounds%begl; endl = bounds%endl

    ! Initialize urban constants

    call urbanparams_inst%Init(bounds)

    ! Initialize vertical data components 

    call initVertical(bounds)

    ! Initialize clm->drv and drv->clm data structures

    call atm2lnd_inst%Init(bounds)
    call lnd2atm_inst%Init(bounds)

    ! Initialization of public data types

    call topo_inst%Init(bounds)

    ! ------------------------------------------------------------------------
    ! Initialize accumulated fields
    ! ------------------------------------------------------------------------

    ! The time manager needs to be initialized before this called is made, since
    ! the step size is needed. 

    call t_startf('init_accflds')
    call atm2lnd_inst%InitAccBuffer(bounds)
    call t_stopf('init_accflds')

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

    call topo_inst%restart (bounds, ncid, flag=flag)

 end subroutine clm_instRest

end module clm_instMod

