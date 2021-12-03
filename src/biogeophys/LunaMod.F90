module LunaMod

#include "shr_assert.h"
  
  !********************************************************************************************************************************************************************** 
  ! !DESCRIPTION:
  ! Calculates the photosynthetic capacities based on a prescribed leaf nitrogen content, using the LUNA model, developed by Chonggang Xu, Ashehad Ali and Rosie Fisher
  ! Currently only works for C3 plants. See Xu et al 2012; Ali et al 2015a. Ecological Applications. http://dx.doi.org/10.1890/14-2111.1.  and Ali et al 2015b.In Review.
  ! !USES:
  use shr_kind_mod          , only : r8  => shr_kind_r8
  use shr_log_mod           , only : errMsg  => shr_log_errMsg
  use clm_varcon            , only : rgas, tfrz,spval
  use abortutils            , only : endrun
  use clm_varctl            , only : iulog
  use clm_varcon            , only : namep 
  use clm_varpar            , only : nlevcan
  use decompMod             , only : bounds_type
  use pftconMod             , only : pftcon
  use FrictionvelocityMod   , only : frictionvel_type 
  use atm2lndType           , only : atm2lnd_type
  use CanopyStateType       , only : canopystate_type
  use PhotosynthesisMod     , only : photosyns_type
  use TemperatureType       , only : temperature_type
  use PatchType             , only : patch
  use GridcellType          , only : grc     
  use SolarAbsorbedType     , only : solarabs_type
  use SurfaceAlbedoType     , only : surfalb_type
  use WaterstateType        , only : waterstate_type
  !use EDPhotosynthesisMod  , only : vcmaxc, jmaxc
  
  
  implicit none
  save
  
  !------------------------------------------------------------------------------
  ! PRIVATE MEMBER FUNCTIONS:
  public  :: LunaReadNML                                   !subroutine to read in the Luna namelist

  !------------------------------------------------------------------------------ 
  !Constants  
  real(r8), parameter :: Cv = 1.2e-5_r8 * 3600.0           ! conversion factor from umol CO2 to g carbon
  real(r8), parameter :: Kc25 = 40.49_r8                   ! Mechanis constant of CO2 for rubisco(Pa), Bernacchi et al (2001) Plant, Cell and Environment 24:253-259
  real(r8), parameter :: Ko25 = 27840_r8                   ! Mechanis constant of O2 for rubisco(Pa), Bernacchi et al (2001) Plant, Cell and Environment 24:253-259
  real(r8), parameter :: Cp25 = 4.275_r8                   ! CO2 compensation point at 25C (Pa), Bernacchi et al (2001) Plant, Cell and Environment 24:253-259
  real(r8), parameter :: Fc25 = 294.2_r8                   ! Fc25 = 6.22*47.3 #see Rogers (2014) Photosynthesis Research 
  real(r8), parameter :: Fj25 = 1257.0_r8                  ! Fj25 = 8.06*156 # #see COSTE 2005 and Xu et al 2012
  real(r8), parameter :: NUEr25 = 33.69_r8                 ! nitrogen use efficiency for respiration, see Xu et al 2012
  real(r8), parameter :: Cb = 1.78_r8                      ! nitrogen use effiency for choloraphyll for light capture, see Evans 1989  
  real(r8), parameter :: O2ref = 209460.0_r8                 ! ppm of O2 in the air
  real(r8), parameter :: CO2ref = 380.0_r8                   ! reference CO2 concentration for calculation of reference NUE. 
  real(r8), parameter :: forc_pbot_ref = 101325.0_r8       ! reference air pressure for calculation of reference NUE
  real(r8), parameter :: Q10Enz = 2.0_r8                   ! Q10 value for enzyme decay rate
  real(r8), parameter :: Jmaxb0 = 0.0311_r8                ! the baseline proportion of nitrogen allocated for electron transport (J)     
  real(r8)            :: Jmaxb1 = 0.1_r8                   ! the baseline proportion of nitrogen allocated for electron transport (J)    
  real(r8), parameter :: Wc2Wjb0 = 0.8054_r8               ! the baseline ratio of rubisco limited rate vs light limited photosynthetic rate (Wc:Wj) 
  real(r8), parameter :: relhExp = 6.0999_r8               ! electron transport parameters related to relative humidity
  real(r8), parameter :: Enzyme_turnover_daily = 0.1_r8    ! the daily turnover rate for photosynthetic enzyme at 25oC in view of ~7 days of half-life time for Rubisco (Suzuki et al. 2001)
  real(r8), parameter :: NMCp25 = 0.715_r8                 ! estimated by assuming 80% maintenance respiration is used for photosynthesis enzyme maintenance
  real(r8), parameter :: Trange1 = 5.0_r8                  ! lower temperature limit (oC) for nitrogen optimization  
  real(r8), parameter :: Trange2 = 42.0_r8                 ! upper temperature limit (oC) for nitrogen optimization
  real(r8), parameter :: SNC = 0.004_r8                    ! structural nitrogen concentration (g N g-1 dry mass carbon)
  real(r8), parameter :: mp = 9.0_r8                       ! slope of stomatal conductance; this is used to estimate model parameter, but may need to be updated from the physiology file, 
  real(r8), parameter :: PARLowLim = 200.0_r8              ! minimum photosynthetically active radiation for nitrogen optimization
  real(r8), parameter :: minrelh = 0.25_r8                 ! minimum relative humdity for nitrogen optimization

  character(len=*), parameter, private :: sourcefile = &
       __FILE__
  !------------------------------------------------------------------------------
  
  contains

  !********************************************************************************************************************************************************************** 
  ! Read in LUNA namelist
  subroutine LunaReadNML( NLFilename )
    !
    ! !DESCRIPTION:
    ! Read the namelist for LUNA
    !
    ! !USES:
    use fileutils      , only : getavu, relavu, opnfil
    use shr_nl_mod     , only : shr_nl_find_group_name
    use spmdMod        , only : masterproc, mpicom
    use shr_mpi_mod    , only : shr_mpi_bcast
    use clm_varctl     , only : iulog
    use shr_log_mod    , only : errMsg => shr_log_errMsg
    use abortutils     , only : endrun
    !
    ! !ARGUMENTS:
    character(len=*), intent(in) :: NLFilename ! Namelist filename
    !
    ! !LOCAL VARIABLES:
    integer :: ierr                 ! error code
    integer :: unitn                ! unit for namelist file

    character(len=*), parameter :: subname = 'lunaReadNML'
    character(len=*), parameter :: nmlname = 'luna'
    !-----------------------------------------------------------------------
    namelist /luna/ Jmaxb1

    ! Initialize options to default values, in case they are not specified in
    ! the namelist


    if (masterproc) then
       unitn = getavu()
       write(iulog,*) 'Read in '//nmlname//'  namelist'
       call opnfil (NLFilename, unitn, 'F')
       call shr_nl_find_group_name(unitn, nmlname, status=ierr)
       if (ierr == 0) then
          read(unitn, nml=luna, iostat=ierr)
          if (ierr /= 0) then
             call endrun(msg="ERROR reading "//nmlname//"namelist"//errmsg(__FILE__, __LINE__))
          end if
       else
          call endrun(msg="ERROR could NOT find "//nmlname//"namelist"//errmsg(__FILE__, __LINE__))
       end if
       call relavu( unitn )
    end if

    call shr_mpi_bcast (Jmaxb1, mpicom)

    if (masterproc) then
       write(iulog,*) ' '
       write(iulog,*) nmlname//' settings:'
       write(iulog,nml=luna)
       write(iulog,*) ' '
    end if

  end subroutine lunaReadNML

end module LunaMod

