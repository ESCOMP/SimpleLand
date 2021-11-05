module CNVegetationFacade

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Facade for the CN Vegetation subsystem.
  !
  ! (A "facade", in software engineering terms, is a unified interface to a set of
  ! interfaces in a subsystem. The facade defines a higher-level interface that makes the
  ! subsystem easier to use.)
  !
  ! NOTE(wjs, 2016-02-19) I envision that we will introduce an abstract base class
  ! (VegBase). Then both CNVeg and EDVeg will extend VegBase. The rest of the CLM code can
  ! then have an instance of VegBase, which depending on the run, can be either a CNVeg or
  ! EDVeg instance.
  !
  ! In addition, we probably want an implementation when running without CN or fates - i.e.,
  ! an SPVeg inst. This would provide implementations for get_leafn_patch,
  ! get_downreg_patch, etc., so that we don't need to handle the non-cn case here (note
  ! that, currently, we return NaN for most of these getters, because these arrays are
  ! invalid and shouldn't be used when running in SP mode). Also, in its EcosystemDynamics
  ! routine, it would call SatellitePhenology (but note that the desired interface for
  ! EcosystemDynamics would be quite different... could just pass everything needed by any
  ! model, and ignore unneeded arguments). Then we can get rid of comments in this module
  ! like, "only call if use_cn is true", as well as use_cn conditionals in this module.
  !
  ! NOTE(wjs, 2016-02-23) Currently, SatellitePhenology is called even when running with
  ! CN, for the sake of dry deposition. This seems weird to me, and my gut feeling -
  ! without understanding it well - is that this should be rewritten to depend on LAI from
  ! CN rather than from satellite phenology. Until that is done, the separation between SP
  ! and other Veg modes will be messier.
  !
  ! NOTE(wjs, 2016-02-23) Currently, this class coordinates calls to soil BGC routines as
  ! well as veg BGC routines (even though it doesn't contain any soil BGC types). This is
  ! because CNDriver coordinates both the veg & soil BGC. We should probably split up
  ! CNDriver so that there is a cleaner separation between veg BGC and soil BGC, to allow
  ! easier swapping of (for example) CN and ED. At that point, this class could
  ! coordinate just the calls to veg BGC routines, with a similar facade class
  ! coordinating the calls to soil BGC routines.
  !
  ! !USES:
#include "shr_assert.h"
  use shr_kind_mod                    , only : r8 => shr_kind_r8
  use shr_infnan_mod                  , only : nan => shr_infnan_nan, assignment(=)
  use shr_log_mod                     , only : errMsg => shr_log_errMsg
  use perf_mod                        , only : t_startf, t_stopf
  use decompMod                       , only : bounds_type
  use clm_varctl                      , only : iulog, use_cn
  use abortutils                      , only : endrun
  use spmdMod                         , only : masterproc
  use CNBalanceCheckMod               , only : cn_balance_type
  use CNVegStateType                  , only : cnveg_state_type
  use CNVegCarbonFluxType             , only : cnveg_carbonflux_type
  use CNVegCarbonStateType            , only : cnveg_carbonstate_type
  use CNVegNitrogenFluxType           , only : cnveg_nitrogenflux_type
  use CNVegNitrogenStateType          , only : cnveg_nitrogenstate_type
  use CNFireMethodMod                 , only : cnfire_method_type
  use CNProductsMod                   , only : cn_products_type
  use SpeciesIsotopeType              , only : species_isotope_type
  use SpeciesNonIsotopeType           , only : species_non_isotope_type
  use CNDriverMod                     , only : CNDriverInit
  !
  implicit none
  private

  ! !PUBLIC TYPES:

  type, public :: cn_vegetation_type
     ! FIXME(bja, 2016-06) These need to be public for use when fates is
     ! turned on. Should either be moved out of here or create some ED
     ! version of the facade....
     type(cnveg_state_type)         :: cnveg_state_inst
     type(cnveg_carbonstate_type)   :: cnveg_carbonstate_inst
     type(cnveg_carbonflux_type)    :: cnveg_carbonflux_inst

     !X!private

     type(cnveg_carbonstate_type)   :: c13_cnveg_carbonstate_inst
     type(cnveg_carbonstate_type)   :: c14_cnveg_carbonstate_inst
     type(cnveg_carbonflux_type)    :: c13_cnveg_carbonflux_inst
     type(cnveg_carbonflux_type)    :: c14_cnveg_carbonflux_inst
     type(cnveg_nitrogenstate_type) :: cnveg_nitrogenstate_inst
     type(cnveg_nitrogenflux_type)  :: cnveg_nitrogenflux_inst

     type(cn_products_type)         :: c_products_inst
     type(cn_products_type)         :: c13_products_inst
     type(cn_products_type)         :: c14_products_inst
     type(cn_products_type)         :: n_products_inst

     type(cn_balance_type)          :: cn_balance_inst
     class(cnfire_method_type), allocatable :: cnfire_method

     ! Control variables
     logical, private :: reseed_dead_plants    ! Flag to indicate if should reseed dead plants when starting up the model

     ! TODO(wjs, 2016-02-19) Evaluate whether some other variables should be moved in
     ! here. Whether they should be moved in depends on how tightly they are tied in with
     ! the other CN Vegetation stuff. A question to ask is: Is this module used when
     ! running with SP or ED? If so, then it should probably remain outside of CNVeg.
     !
     ! From the clm_instMod section on "CN vegetation types":
     ! - nutrient_competition_method
     !   - I'm pretty sure this should be moved into here; it's just a little messy to do
     !     so, because of how it's initialized (specifically, the call to readParameters
     !     in clm_initializeMod).
     !
     ! From the clm_instMod section on "general biogeochem types":
     ! - ch4_inst
     !   - probably not: really seems to belong in soilbiogeochem
     ! - crop_inst
     ! - dust_inst
     ! - vocemis_inst
     ! - fireemis_inst
     ! - drydepvel_inst
     
   contains
     procedure, public :: Init
     procedure, public :: InitAccBuffer
     procedure, public :: InitAccVars
     procedure, public :: Restart

     procedure, public :: Init2                         ! Do initialization in initialize phase, after subgrid weights are determined
     procedure, public :: WriteHistory                  ! Do any history writes that are specific to veg dynamics

     procedure, public :: get_totvegc_col               ! Get column-level total vegetation carbon array

     procedure, private :: CNReadNML                    ! Read in the CN general namelist
  end type cn_vegetation_type

  character(len=*), parameter, private :: sourcefile = &
       __FILE__

contains

  !-----------------------------------------------------------------------
  subroutine Init(this, bounds, NLFilename)
    !
    ! !DESCRIPTION:
    ! Initialize a CNVeg object.
    !
    ! Should be called regardless of whether use_cn is true
    !
    ! !USES:
    use CNFireFactoryMod , only : create_cnfire_method
    use clm_varcon       , only : c13ratio, c14ratio
    !
    ! !ARGUMENTS:
    class(cn_vegetation_type), intent(inout) :: this
    type(bounds_type), intent(in)    :: bounds
    character(len=*) , intent(in)    :: NLFilename ! namelist filename
    !
    ! !LOCAL VARIABLES:
    integer :: begp, endp

    character(len=*), parameter :: subname = 'Init'
    !-----------------------------------------------------------------------

    begp = bounds%begp
    endp = bounds%endp

    ! Note - always initialize the memory for cnveg_state_inst (used in biogeophys/)
    call this%cnveg_state_inst%Init(bounds)
    
    if (use_cn) then

       ! Read in the general CN namelist
       call this%CNReadNML( NLFilename )    ! MUST be called first as passes down control information to others

       call this%cnveg_carbonstate_inst%Init(bounds, carbon_type='c12', ratio=1._r8, NLFilename=NLFilename)
       call this%cnveg_carbonflux_inst%Init(bounds, carbon_type='c12')
       call this%cnveg_nitrogenstate_inst%Init(bounds,                   &
            this%cnveg_carbonstate_inst%leafc_patch(begp:endp),          &
            this%cnveg_carbonstate_inst%leafc_storage_patch(begp:endp),  &
            this%cnveg_carbonstate_inst%frootc_patch(begp:endp),         &
            this%cnveg_carbonstate_inst%frootc_storage_patch(begp:endp), &
            this%cnveg_carbonstate_inst%deadstemc_patch(begp:endp) )
       call this%cnveg_nitrogenflux_inst%Init(bounds) 

       call this%c_products_inst%Init(bounds, species_non_isotope_type('C'))
       call this%n_products_inst%Init(bounds, species_non_isotope_type('N'))

       call this%cn_balance_inst%Init(bounds)

    end if

    allocate(this%cnfire_method, &
         source=create_cnfire_method(NLFilename))

  end subroutine Init

  !-----------------------------------------------------------------------
  subroutine CNReadNML( this, NLFilename )
    !
    ! !DESCRIPTION:
    ! Read in the general CN control namelist
    !
    ! !USES:
    use fileutils      , only : getavu, relavu, opnfil
    use shr_nl_mod     , only : shr_nl_find_group_name
    use spmdMod        , only : masterproc, mpicom
    use shr_mpi_mod    , only : shr_mpi_bcast
    use clm_varctl     , only : iulog
    !
    ! !ARGUMENTS:
    class(cn_vegetation_type), intent(inout) :: this
    character(len=*)         , intent(in)    :: NLFilename                 ! Namelist filename
    !
    ! !LOCAL VARIABLES:
    integer :: ierr                 ! error code
    integer :: unitn                ! unit for namelist file

    character(len=*), parameter :: subname = 'CNReadNML'
    character(len=*), parameter :: nmlname = 'cn_general'   ! MUST match what is in namelist below
    !-----------------------------------------------------------------------
    logical :: reseed_dead_plants
    namelist /cn_general/ reseed_dead_plants

    reseed_dead_plants = this%reseed_dead_plants

    if (masterproc) then
       unitn = getavu()
       write(iulog,*) 'Read in '//nmlname//'  namelist'
       call opnfil (NLFilename, unitn, 'F')
       call shr_nl_find_group_name(unitn, nmlname, status=ierr)
       if (ierr == 0) then
          read(unitn, nml=cn_general, iostat=ierr)   ! Namelist name here MUST be the same as in nmlname above!
          if (ierr /= 0) then
             call endrun(msg="ERROR reading "//nmlname//"namelist"//errmsg(sourcefile, __LINE__))
          end if
       else
          call endrun(msg="ERROR could NOT find "//nmlname//"namelist"//errmsg(sourcefile, __LINE__))
       end if
       call relavu( unitn )
    end if

    call shr_mpi_bcast (reseed_dead_plants      , mpicom)

    this%reseed_dead_plants = reseed_dead_plants

    if (masterproc) then
       write(iulog,*) ' '
       write(iulog,*) nmlname//' settings:'
       write(iulog,nml=cn_general)    ! Name here MUST be the same as in nmlname above!
       write(iulog,*) ' '
    end if

    !-----------------------------------------------------------------------

  end subroutine CNReadNML


  !-----------------------------------------------------------------------
  subroutine InitAccBuffer(this, bounds)
    !
    ! !DESCRIPTION:
    ! Initialize accumulation buffer for types contained here 
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    class(cn_vegetation_type), intent(inout) :: this
    type(bounds_type), intent(in)    :: bounds
    !
    ! !LOCAL VARIABLES:

    character(len=*), parameter :: subname = 'InitAccBuffer'
    !-----------------------------------------------------------------------

  end subroutine InitAccBuffer

  !-----------------------------------------------------------------------
  subroutine InitAccVars(this, bounds)
    !
    ! !DESCRIPTION:
    ! Initialize variables that are associated with accumulated fields
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    class(cn_vegetation_type), intent(inout) :: this
    type(bounds_type), intent(in)    :: bounds
    !
    ! !LOCAL VARIABLES:

    character(len=*), parameter :: subname = 'InitAccVars'
    !-----------------------------------------------------------------------

  end subroutine InitAccVars

  !-----------------------------------------------------------------------
  subroutine Restart(this, bounds, ncid, flag)
    !
    ! !DESCRIPTION:
    ! Handle restart (read / write) for CNVeg
    !
    ! Should be called regardless of whether use_cn is true
    !
    ! !USES:
    use ncdio_pio, only : file_desc_t
    use clm_varcon, only : c3_r2, c14ratio
    !
    ! !ARGUMENTS:
    class(cn_vegetation_type), intent(inout) :: this
    type(bounds_type), intent(in)    :: bounds 
    type(file_desc_t), intent(inout) :: ncid   
    character(len=*) , intent(in)    :: flag   
    integer  :: reseed_patch(bounds%endp-bounds%begp+1)
    integer  :: num_reseed_patch
    !
    ! !LOCAL VARIABLES:

    integer :: begp, endp

    character(len=*), parameter :: subname = 'Restart'
    !-----------------------------------------------------------------------

    if (use_cn) then
       begp = bounds%begp
       endp = bounds%endp
       call this%cnveg_carbonstate_inst%restart(bounds, ncid, flag=flag, carbon_type='c12', &
               reseed_dead_plants=this%reseed_dead_plants, filter_reseed_patch=reseed_patch, &
               num_reseed_patch=num_reseed_patch )
       if ( flag /= 'read' .and. num_reseed_patch /= 0 )then
          call endrun(msg="ERROR num_reseed should be zero and is not"//errmsg(sourcefile, __LINE__))
       end if
       call this%cnveg_carbonflux_inst%restart(bounds, ncid, flag=flag, carbon_type='c12')

       call this%cnveg_nitrogenstate_inst%restart(bounds, ncid, flag=flag,  &
            leafc_patch=this%cnveg_carbonstate_inst%leafc_patch(begp:endp),         &
            leafc_storage_patch=this%cnveg_carbonstate_inst%leafc_storage_patch(begp:endp), &
            frootc_patch=this%cnveg_carbonstate_inst%frootc_patch(begp:endp), &
            frootc_storage_patch=this%cnveg_carbonstate_inst%frootc_storage_patch(begp:endp), &
            deadstemc_patch=this%cnveg_carbonstate_inst%deadstemc_patch(begp:endp), &
            filter_reseed_patch=reseed_patch, num_reseed_patch=num_reseed_patch)
       call this%cnveg_nitrogenflux_inst%restart(bounds, ncid, flag=flag)
       call this%cnveg_state_inst%restart(bounds, ncid, flag=flag, &
            cnveg_carbonstate=this%cnveg_carbonstate_inst, &
            cnveg_nitrogenstate=this%cnveg_nitrogenstate_inst, &
            filter_reseed_patch=reseed_patch, num_reseed_patch=num_reseed_patch)

       call this%c_products_inst%restart(bounds, ncid, flag)
       call this%n_products_inst%restart(bounds, ncid, flag)

    end if

  end subroutine Restart

  !-----------------------------------------------------------------------
  subroutine Init2(this, bounds, NLFilename)
    !
    ! !DESCRIPTION:
    ! Do initialization that is needed in the initialize phase, after subgrid weights are
    ! determined
    !
    ! Should only be called if use_cn is true
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    class(cn_vegetation_type) , intent(inout) :: this
    type(bounds_type) , intent(in)    :: bounds
    character(len=*)  , intent(in)    :: NLFilename ! namelist filename
    !
    ! !LOCAL VARIABLES:

    character(len=*), parameter :: subname = 'Init2'
    !-----------------------------------------------------------------------

    call CNDriverInit(bounds, NLFilename, this%cnfire_method)

  end subroutine Init2


  !-----------------------------------------------------------------------
  subroutine WriteHistory(this, bounds)
    !
    ! !DESCRIPTION:
    ! Do any history writes that are specific to vegetation dynamics
    !
    ! NOTE(wjs, 2016-02-23) This could probably be combined with
    ! EndOfTimeStepVegDynamics, except for the fact that (currently) history writes are
    ! done with proc bounds rather than clump bounds. If that were changed, then the body
    ! of this could be moved into EndOfTimeStepVegDynamics, inside a "if (.not.
    ! use_noio)" conditional.
    !
    ! Should only be called if use_cn is true
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    class(cn_vegetation_type), intent(in) :: this
    type(bounds_type)  , intent(in) :: bounds                  
    !
    ! !LOCAL VARIABLES:

    character(len=*), parameter :: subname = 'WriteHistory'
    !-----------------------------------------------------------------------

  end subroutine WriteHistory

  !-----------------------------------------------------------------------
  function get_totvegc_col(this, bounds) result(totvegc_col)
    !
    ! !DESCRIPTION:
    ! Get column-level total vegetation carbon array
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    class(cn_vegetation_type), intent(in) :: this
    type(bounds_type), intent(in) :: bounds
    real(r8) :: totvegc_col(bounds%begc:bounds%endc)  ! function result: (gC/m2)
    !
    ! !LOCAL VARIABLES:

    character(len=*), parameter :: subname = 'get_totvegc_col'
    !-----------------------------------------------------------------------

    if (use_cn) then
       totvegc_col(bounds%begc:bounds%endc) = &
            this%cnveg_carbonstate_inst%totvegc_col(bounds%begc:bounds%endc)
    else
       totvegc_col(bounds%begc:bounds%endc) = nan
    end if

  end function get_totvegc_col


end module CNVegetationFacade
