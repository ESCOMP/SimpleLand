module clm_instMod

  !------------------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Instances and definitions of all data types
  !
  ! !USES:
  use shr_kind_mod    , only : r8 => shr_kind_r8
  use decompMod       , only : bounds_type
  use clm_varcon      , only : bdsno
  use landunit_varcon , only : istice_mec, istsoil
  use perf_mod        , only : t_startf, t_stopf
  use controlMod      , only : NLFilename

  !-----------------------------------------
  ! Constants
  !-----------------------------------------

  use UrbanParamsType                    , only : urbanparams_type   ! Constants 
  !-----------------------------------------
  ! Definition of component types 
  !-----------------------------------------

  use EnergyFluxType                  , only : energyflux_type
  use FrictionVelocityMod             , only : frictionvel_type
  use SoilHydrologyType               , only : soilhydrology_type  
  use SoilStateType                   , only : soilstate_type
  use SolarAbsorbedType               , only : solarabs_type
  use SurfaceAlbedoType               , only : surfalb_type
  use TemperatureType                 , only : temperature_type
  use WaterFluxType                   , only : waterflux_type
  use WaterStateType                  , only : waterstate_type
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
  use SoilStateInitTimeConstMod       , only : SoilStateInitTimeConst
  use SoilHydrologyInitTimeConstMod   , only : SoilHydrologyInitTimeConst
  !
  implicit none
  public   ! By default everything is public 
  !
  !-----------------------------------------
  ! Instances of component types
  !-----------------------------------------

  ! Physics types 
  type(energyflux_type)                   :: energyflux_inst
  type(frictionvel_type)                  :: frictionvel_inst
  type(soilstate_type)                    :: soilstate_inst
  type(soilhydrology_type)                :: soilhydrology_inst
  type(solarabs_type)                     :: solarabs_inst
  type(surfalb_type)                      :: surfalb_inst
  type(temperature_type)                  :: temperature_inst
  type(urbanparams_type)                  :: urbanparams_inst
  type(waterflux_type)                    :: waterflux_inst
  type(waterstate_type)                   :: waterstate_inst
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
    use clm_varpar                         , only : nlevsno, numpft
    use controlMod                         , only : nlfilename, fsurdat
    use domainMod                          , only : ldomain
    use initVerticalMod                    , only : initVertical
    use accumulMod                         , only : print_accum_fields 
    use decompMod                          , only : get_proc_bounds
    !
    ! !ARGUMENTS    
    type(bounds_type), intent(in) :: bounds  ! processor bounds
    !
    ! !LOCAL VARIABLES:
    integer               :: c,l,g
    integer               :: nclumps,nc
    integer               :: begp, endp
    integer               :: begc, endc
    integer               :: begl, endl
    type(bounds_type)     :: bounds_clump
    real(r8), allocatable :: h2osno_col(:)
    real(r8), allocatable :: snow_depth_col(:)

    integer :: dummy_to_make_pgi_happy
    !----------------------------------------------------------------------

    ! Note: h2osno_col and snow_depth_col are initialized as local variable 
    ! since they are needed to initialize vertical data structures  

    begp = bounds%begp; endp = bounds%endp 
    begc = bounds%begc; endc = bounds%endc 
    begl = bounds%begl; endl = bounds%endl

    allocate (h2osno_col(begc:endc))
    allocate (snow_depth_col(begc:endc))

    ! snow water
    do c = begc,endc
       l = col%landunit(c)
       g = col%gridcell(c)

       ! In areas that should be snow-covered, it can be problematic to start with 0 snow
       ! cover, because this can affect the long-term state through soil heating, albedo
       ! feedback, etc. On the other hand, we would introduce hysteresis by putting too
       ! much snow in places that are in a net melt regime, because the melt-albedo
       ! feedback may not activate on time (or at all). So, as a compromise, we start with
       ! a small amount of snow in places that are likely to be snow-covered for much or
       ! all of the year.
       if (lun%itype(l)==istice_mec) then
          h2osno_col(c) = 100._r8
       else if (lun%itype(l)==istsoil .and. abs(grc%latdeg(g)) >= 60._r8) then 
          h2osno_col(c) = 100._r8
       else
          h2osno_col(c) = 0._r8
       endif
       snow_depth_col(c)  = h2osno_col(c) / bdsno
    end do

    ! Initialize urban constants

    call urbanparams_inst%Init(bounds)

    ! Initialize vertical data components 

    call initVertical(bounds,               &
         glc_behavior, &
         snow_depth_col(begc:endc),              &
         urbanparams_inst%thick_wall(begl:endl), &
         urbanparams_inst%thick_roof(begl:endl))

    ! Initialize clm->drv and drv->clm data structures

    call atm2lnd_inst%Init( bounds, NLFilename )
    call lnd2atm_inst%Init( bounds, NLFilename )

    ! Initialization of public data types

    call temperature_inst%Init(bounds)

    call soilstate_inst%Init(bounds)
    call SoilStateInitTimeConst(bounds, soilstate_inst, nlfilename) ! sets hydraulic and thermal soil properties

    call waterstate_inst%Init(bounds,         &
         h2osno_col(begc:endc),                    &
         snow_depth_col(begc:endc),                &
         soilstate_inst%watsat_col(begc:endc, 1:), &
         temperature_inst%t_soisno_col(begc:endc, -nlevsno+1:) )

    call waterflux_inst%Init(bounds)

    ! COMPILER_BUG(wjs, 2014-11-29, pgi 14.7) Without the following assignment, the
    ! assertion in energyflux_inst%Init fails with pgi 14.7 on yellowstone, presumably due
    ! to a compiler bug.
    dummy_to_make_pgi_happy = ubound(temperature_inst%t_grnd_col, 1)
    call energyflux_inst%Init(bounds, temperature_inst%t_grnd_col(begc:endc))

    call frictionvel_inst%Init(bounds)

    call soilhydrology_inst%Init(bounds, nlfilename)
    call SoilHydrologyInitTimeConst(bounds, soilhydrology_inst) ! sets time constant properties

    call solarabs_inst%Init(bounds)

    call surfalb_inst%Init(bounds)

    call topo_inst%Init(bounds)

    deallocate (h2osno_col)
    deallocate (snow_depth_col)

    ! ------------------------------------------------------------------------
    ! Initialize accumulated fields
    ! ------------------------------------------------------------------------

    ! The time manager needs to be initialized before this called is made, since
    ! the step size is needed. 

    call t_startf('init_accflds')

    call atm2lnd_inst%InitAccBuffer(bounds)

    call temperature_inst%InitAccBuffer(bounds)
    
    call energyflux_inst%InitAccBuffer(bounds)

    call t_stopf('init_accflds')

  end subroutine clm_instInit

  !-----------------------------------------------------------------------
  subroutine clm_instRest(bounds, ncid, flag)
    !
    ! !USES:
    use ncdio_pio       , only : file_desc_t
    use decompMod       , only : get_proc_bounds, get_proc_clumps, get_clump_bounds

    !
    ! !DESCRIPTION:
    ! Define/write/read CLM restart file.
    !
    ! !ARGUMENTS:
    type(bounds_type) , intent(in)    :: bounds          
    
    type(file_desc_t) , intent(inout) :: ncid ! netcdf id
    character(len=*)  , intent(in)    :: flag ! 'define', 'write', 'read' 

    ! Local variables
    integer                           :: nc, nclumps
    type(bounds_type)                 :: bounds_clump

    !-----------------------------------------------------------------------

    call atm2lnd_inst%restart (bounds, ncid, flag=flag)

    call energyflux_inst%restart (bounds, ncid, flag=flag)

    call frictionvel_inst% restart (bounds, ncid, flag=flag)

    call soilhydrology_inst%restart (bounds, ncid, flag=flag)

    call temperature_inst%restart (bounds, ncid, flag=flag)

    call soilstate_inst%restart (bounds, ncid, flag=flag)

    call waterflux_inst%restart (bounds, ncid, flag=flag)

    call waterstate_inst%restart (bounds, ncid, flag=flag, &
         watsat_col=soilstate_inst%watsat_col(bounds%begc:bounds%endc,:)) 

    call surfalb_inst%restart (bounds, ncid, flag=flag)

    call topo_inst%restart (bounds, ncid, flag=flag)

 end subroutine clm_instRest

end module clm_instMod

