module CanopyFluxesMod

#include "shr_assert.h"

  !------------------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Performs calculation of leaf temperature and surface fluxes.
  ! SoilFluxes then determines soil/snow and ground temperatures and updates the surface 
  ! fluxes for the new ground temperature.
  !
  ! !USES:
  use shr_sys_mod           , only : shr_sys_flush
  use shr_kind_mod          , only : r8 => shr_kind_r8
  use shr_log_mod           , only : errMsg => shr_log_errMsg
  use abortutils            , only : endrun
  use clm_varctl            , only : iulog, use_cn, use_cndv, &
                                     use_luna, use_hydrstress
  use clm_varpar            , only : nlevgrnd, nlevsno
  use clm_varcon            , only : namep 
  use pftconMod             , only : pftcon
  use decompMod             , only : bounds_type
  !
  ! !PUBLIC TYPES:
  implicit none
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: CanopyFluxesReadNML     ! Read in namelist settings
  !
  ! !PUBLIC DATA MEMBERS:
  ! true => btran is based only on unfrozen soil levels
  logical,  public :: perchroot     = .false.  

  ! true  => btran is based on active layer (defined over two years); 
  ! false => btran is based on currently unfrozen levels
  logical,  public :: perchroot_alt = .false.  
  !
  ! !PRIVATE DATA MEMBERS:
  ! Snow in vegetation canopy namelist options.
  logical, private :: snowveg_on     = .false.                ! snowveg_flag = 'ON'
  logical, private :: snowveg_onrad  = .true.                 ! snowveg_flag = 'ON_RAD'
  logical, private :: use_undercanopy_stability = .true.      ! use undercanopy stability term or not

  character(len=*), parameter, private :: sourcefile = &
       __FILE__
  !------------------------------------------------------------------------------

contains

  !------------------------------------------------------------------------
  subroutine CanopyFluxesReadNML(NLFilename)
    !
    ! !DESCRIPTION:
    ! Read the namelist for Canopy Fluxes
    !
    ! !USES:
    use fileutils      , only : getavu, relavu, opnfil
    use shr_nl_mod     , only : shr_nl_find_group_name
    use spmdMod        , only : masterproc, mpicom
    use shr_mpi_mod    , only : shr_mpi_bcast
    use clm_varctl     , only : iulog
    !
    ! !ARGUMENTS:
    character(len=*), intent(IN) :: NLFilename ! Namelist filename
    !
    ! !LOCAL VARIABLES:
    integer :: ierr                 ! error code
    integer :: unitn                ! unit for namelist file

    character(len=*), parameter :: subname = 'CanopyFluxeseadNML'
    character(len=*), parameter :: nmlname = 'canopyfluxes_inparm'
    !-----------------------------------------------------------------------

    namelist /canopyfluxes_inparm/ use_undercanopy_stability

    ! Initialize options to default values, in case they are not specified in
    ! the namelist

    if (masterproc) then
       unitn = getavu()
       write(iulog,*) 'Read in '//nmlname//'  namelist'
       call opnfil (NLFilename, unitn, 'F')
       call shr_nl_find_group_name(unitn, nmlname, status=ierr)
       if (ierr == 0) then
          read(unitn, nml=canopyfluxes_inparm, iostat=ierr)
          if (ierr /= 0) then
             call endrun(msg="ERROR reading "//nmlname//"namelist"//errmsg(sourcefile, __LINE__))
          end if
       else
          call endrun(msg="ERROR could NOT find "//nmlname//"namelist"//errmsg(sourcefile, __LINE__))
       end if
       call relavu( unitn )
    end if

    call shr_mpi_bcast (use_undercanopy_stability, mpicom)

    if (masterproc) then
       write(iulog,*) ' '
       write(iulog,*) nmlname//' settings:'
       write(iulog,nml=canopyfluxes_inparm)
       write(iulog,*) ' '
    end if

  end subroutine CanopyFluxesReadNML

end module CanopyFluxesMod

