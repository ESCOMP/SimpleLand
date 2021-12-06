module SurfaceAlbedoMod

#include "shr_assert.h"

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Performs surface albedo calculations
  !
  ! !PUBLIC TYPES:
  use shr_kind_mod      , only : r8 => shr_kind_r8
  use shr_log_mod       , only : errMsg => shr_log_errMsg
  use decompMod         , only : bounds_type
  use landunit_varcon   , only : istsoil, istcrop
  use clm_varcon        , only : grlnd, namep
  use clm_varpar        , only : numrad, nlevcan, nlevsno, nlevcan
  use clm_varctl        , only : fsurdat, iulog, use_snicar_frc
  use pftconMod         , only : pftcon
  use ColumnType        , only : col                
  !
  implicit none
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: SurfaceAlbedoInitTimeConst
  !
  ! !PUBLIC DATA MEMBERS:
  ! The CLM default albice values are too high.
  ! Full-spectral albedo for land ice is ~0.5 (Paterson, Physics of Glaciers, 1994, p. 59)
  ! This is the value used in CAM3 by Pritchard et al., GRL, 35, 2008.

  ! albedo land ice by waveband (1=vis, 2=nir)
  real(r8), public  :: albice(numrad) = (/ 0.80_r8, 0.55_r8 /)

  ! namelist default setting for inputting alblakwi
  real(r8), public  :: lake_melt_icealb(numrad) = (/ 0.10_r8, 0.10_r8/)

  ! Coefficient for calculating ice "fraction" for lake surface albedo
  ! From D. Mironov (2010) Boreal Env. Research
  real(r8), parameter :: calb = 95.6_r8   

  !
  ! !PRIVATE DATA MEMBERS:

  ! !PRIVATE DATA FUNCTIONS:
  real(r8), allocatable, private :: albsat(:,:) ! wet soil albedo by color class and waveband (1=vis,2=nir)
  real(r8), allocatable, private :: albdry(:,:) ! dry soil albedo by color class and waveband (1=vis,2=nir)
  integer , allocatable, private :: isoicol(:)  ! column soil color class

  character(len=*), parameter, private :: sourcefile = &
       __FILE__
  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine SurfaceAlbedoInitTimeConst(bounds)
    !
    ! !DESCRIPTION:
    ! Initialize module time constant variables
    !
    ! !USES:
    use shr_log_mod, only : errMsg => shr_log_errMsg
    use fileutils  , only : getfil
    use abortutils , only : endrun
    use ncdio_pio  , only : file_desc_t, ncd_defvar, ncd_io, ncd_pio_openfile, ncd_pio_closefile
    use spmdMod    , only : masterproc
    !
    ! !ARGUMENTS:
    type(bounds_type), intent(in) :: bounds  
    !
    ! !LOCAL VARIABLES:
    integer            :: c,g          ! indices
    integer            :: mxsoil_color ! maximum number of soil color classes
    type(file_desc_t)  :: ncid         ! netcdf id
    character(len=256) :: locfn        ! local filename
    integer            :: ier          ! error status
    logical            :: readvar 
    integer  ,pointer  :: soic2d (:)   ! read in - soil color 
    !---------------------------------------------------------------------

    ! Allocate module variable for soil color

    allocate(isoicol(bounds%begc:bounds%endc)) 

    ! Determine soil color and number of soil color classes 
    ! if number of soil color classes is not on input dataset set it to 8

    call getfil (fsurdat, locfn, 0)
    call ncd_pio_openfile (ncid, locfn, 0)

    call ncd_io(ncid=ncid, varname='mxsoil_color', flag='read', data=mxsoil_color, readvar=readvar)
    if ( .not. readvar ) mxsoil_color = 8  

    allocate(soic2d(bounds%begg:bounds%endg)) 
    call ncd_io(ncid=ncid, varname='SOIL_COLOR', flag='read', data=soic2d, dim1name=grlnd, readvar=readvar)
    if (.not. readvar) then
       call endrun(msg=' ERROR: SOIL_COLOR NOT on surfdata file'//errMsg(sourcefile, __LINE__)) 
    end if
    do c = bounds%begc, bounds%endc
       g = col%gridcell(c)
       isoicol(c) = soic2d(g)
    end do
    deallocate(soic2d)

    call ncd_pio_closefile(ncid)

    ! Determine saturated and dry soil albedos for n color classes and 
    ! numrad wavebands (1=vis, 2=nir)

    allocate(albsat(mxsoil_color,numrad), albdry(mxsoil_color,numrad), stat=ier)
    if (ier /= 0) then
       write(iulog,*)'allocation error for albsat, albdry'
       call endrun(msg=errMsg(sourcefile, __LINE__)) 
    end if

    if (masterproc) then
       write(iulog,*) 'Attempting to read soil colo data .....'
    end if
    
    if (mxsoil_color == 8) then
       albsat(1:8,1) = (/0.12_r8,0.11_r8,0.10_r8,0.09_r8,0.08_r8,0.07_r8,0.06_r8,0.05_r8/)
       albsat(1:8,2) = (/0.24_r8,0.22_r8,0.20_r8,0.18_r8,0.16_r8,0.14_r8,0.12_r8,0.10_r8/)
       albdry(1:8,1) = (/0.24_r8,0.22_r8,0.20_r8,0.18_r8,0.16_r8,0.14_r8,0.12_r8,0.10_r8/)
       albdry(1:8,2) = (/0.48_r8,0.44_r8,0.40_r8,0.36_r8,0.32_r8,0.28_r8,0.24_r8,0.20_r8/)
    else if (mxsoil_color == 20) then
       albsat(1:20,1) = (/0.25_r8,0.23_r8,0.21_r8,0.20_r8,0.19_r8,0.18_r8,0.17_r8,0.16_r8,&
            0.15_r8,0.14_r8,0.13_r8,0.12_r8,0.11_r8,0.10_r8,0.09_r8,0.08_r8,0.07_r8,0.06_r8,0.05_r8,0.04_r8/)
       albsat(1:20,2) = (/0.50_r8,0.46_r8,0.42_r8,0.40_r8,0.38_r8,0.36_r8,0.34_r8,0.32_r8,&
            0.30_r8,0.28_r8,0.26_r8,0.24_r8,0.22_r8,0.20_r8,0.18_r8,0.16_r8,0.14_r8,0.12_r8,0.10_r8,0.08_r8/)
       albdry(1:20,1) = (/0.36_r8,0.34_r8,0.32_r8,0.31_r8,0.30_r8,0.29_r8,0.28_r8,0.27_r8,&
            0.26_r8,0.25_r8,0.24_r8,0.23_r8,0.22_r8,0.20_r8,0.18_r8,0.16_r8,0.14_r8,0.12_r8,0.10_r8,0.08_r8/)
       albdry(1:20,2) = (/0.61_r8,0.57_r8,0.53_r8,0.51_r8,0.49_r8,0.48_r8,0.45_r8,0.43_r8,&
            0.41_r8,0.39_r8,0.37_r8,0.35_r8,0.33_r8,0.31_r8,0.29_r8,0.27_r8,0.25_r8,0.23_r8,0.21_r8,0.16_r8/)
    else
       write(iulog,*)'maximum color class = ',mxsoil_color,' is not supported'
       call endrun(msg=errMsg(sourcefile, __LINE__)) 
    end if

    ! Set alblakwi
    !alblakwi(:) = lake_melt_icealb(:)

  end subroutine SurfaceAlbedoInitTimeConst

end module SurfaceAlbedoMod
