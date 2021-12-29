module SnowSnicarMod

#include "shr_assert.h"

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Calculate albedo of snow containing impurities 
  ! and the evolution of snow effective radius
  !
  ! !USES:
  use shr_kind_mod    , only : r8 => shr_kind_r8
  use shr_sys_mod     , only : shr_sys_flush
  use shr_log_mod     , only : errMsg => shr_log_errMsg
  use clm_varctl      , only : iulog
  use clm_varcon      , only : namec , tfrz
  use shr_const_mod   , only : SHR_CONST_RHOICE
  use abortutils      , only : endrun
  use decompMod       , only : bounds_type
  use AerosolMod      , only : snw_rds_min
  !
  implicit none
  save
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: SnowAge_init     ! Initial read in of snow-aging file
  public :: SnowOptics_init  ! Initial read in of snow-optics file
  !
  ! !PUBLIC DATA MEMBERS:
  integer,  public, parameter :: sno_nbr_aer =   8        ! number of aerosol species in snowpack
                                                          ! (indices described above) [nbr]
  logical,  public, parameter :: DO_SNO_OC =    .false.   ! parameter to include organic carbon (OC)
                                                          ! in snowpack radiative calculations
  logical,  public, parameter :: DO_SNO_AER =   .true.    ! parameter to include aerosols in snowpack radiative calculations

  ! !PRIVATE DATA MEMBERS:
  integer,  parameter :: numrad_snw  =   5               ! number of spectral bands used in snow model [nbr]
  integer,  parameter :: nir_bnd_bgn =   2               ! first band index in near-IR spectrum [idx]
  integer,  parameter :: nir_bnd_end =   5               ! ending near-IR band index [idx]

  integer,  parameter :: idx_Mie_snw_mx = 1471           ! number of effective radius indices used in Mie lookup table [idx]
  integer,  parameter :: idx_T_max      = 11             ! maxiumum temperature index used in aging lookup table [idx]
  integer,  parameter :: idx_T_min      = 1              ! minimum temperature index used in aging lookup table [idx]
  integer,  parameter :: idx_Tgrd_max   = 31             ! maxiumum temperature gradient index used in aging lookup table [idx]
  integer,  parameter :: idx_Tgrd_min   = 1              ! minimum temperature gradient index used in aging lookup table [idx]
  integer,  parameter :: idx_rhos_max   = 8              ! maxiumum snow density index used in aging lookup table [idx]
  integer,  parameter :: idx_rhos_min   = 1              ! minimum snow density index used in aging lookup table [idx]

  integer,  parameter :: snw_rds_max_tbl = 1500          ! maximum effective radius defined in Mie lookup table [microns]
  integer,  parameter :: snw_rds_min_tbl = 30            ! minimium effective radius defined in Mie lookup table [microns]
  integer,  parameter :: snw_rds_min_int = nint(snw_rds_min) ! minimum allowed snow effective radius as integer [microns]
  real(r8), parameter :: snw_rds_max     = 1500._r8      ! maximum allowed snow effective radius [microns]
  real(r8), parameter :: snw_rds_refrz   = 1000._r8      ! effective radius of re-frozen snow [microns]

  real(r8), parameter :: min_snw = 1.0E-30_r8            ! minimum snow mass required for SNICAR RT calculation [kg m-2]

  !real(r8), parameter :: C1_liq_Brun89 = 1.28E-17_r8    ! constant for liquid water grain growth [m3 s-1],
                                                         ! from Brun89
  real(r8), parameter :: C1_liq_Brun89 = 0._r8           ! constant for liquid water grain growth [m3 s-1],
                                                         ! from Brun89: zeroed to accomodate dry snow aging
  real(r8), parameter :: C2_liq_Brun89 = 4.22E-13_r8     ! constant for liquid water grain growth [m3 s-1],
                                                         ! from Brun89: corrected for LWC in units of percent

  real(r8), parameter :: tim_cns_bc_rmv  = 2.2E-8_r8     ! time constant for removal of BC in snow on sea-ice
                                                         ! [s-1] (50% mass removal/year)
  real(r8), parameter :: tim_cns_oc_rmv  = 2.2E-8_r8     ! time constant for removal of OC in snow on sea-ice
                                                         ! [s-1] (50% mass removal/year)
  real(r8), parameter :: tim_cns_dst_rmv = 2.2E-8_r8     ! time constant for removal of dust in snow on sea-ice
                                                         ! [s-1] (50% mass removal/year)

  ! scaling of the snow aging rate (tuning option):
  logical :: flg_snoage_scl    = .false.                 ! flag for scaling the snow aging rate by some arbitrary factor
  real(r8), parameter :: xdrdt = 1.0_r8                  ! arbitrary factor applied to snow aging rate

  ! snow and aerosol Mie parameters:
  ! (arrays declared here, but are set in iniTimeConst)
  ! (idx_Mie_snw_mx is number of snow radii with defined parameters (i.e. from 30um to 1500um))
  
  ! direct-beam weighted ice optical properties
  real(r8) :: ss_alb_snw_drc(idx_Mie_snw_mx,numrad_snw)
  real(r8) :: asm_prm_snw_drc(idx_Mie_snw_mx,numrad_snw)
  real(r8) :: ext_cff_mss_snw_drc(idx_Mie_snw_mx,numrad_snw)

  ! diffuse radiation weighted ice optical properties
  real(r8) :: ss_alb_snw_dfs(idx_Mie_snw_mx,numrad_snw)
  real(r8) :: asm_prm_snw_dfs(idx_Mie_snw_mx,numrad_snw)
  real(r8) :: ext_cff_mss_snw_dfs(idx_Mie_snw_mx,numrad_snw)

  ! hydrophiliic BC
  real(r8) :: ss_alb_bc1(numrad_snw)
  real(r8) :: asm_prm_bc1(numrad_snw)
  real(r8) :: ext_cff_mss_bc1(numrad_snw)

  ! hydrophobic BC
  real(r8) :: ss_alb_bc2(numrad_snw)
  real(r8) :: asm_prm_bc2(numrad_snw)
  real(r8) :: ext_cff_mss_bc2(numrad_snw)

  ! hydrophobic OC
  real(r8) :: ss_alb_oc1(numrad_snw)
  real(r8) :: asm_prm_oc1(numrad_snw)
  real(r8) :: ext_cff_mss_oc1(numrad_snw)

  ! hydrophilic OC
  real(r8) :: ss_alb_oc2(numrad_snw)
  real(r8) :: asm_prm_oc2(numrad_snw)
  real(r8) :: ext_cff_mss_oc2(numrad_snw)

  ! dust species 1:
  real(r8) :: ss_alb_dst1(numrad_snw)
  real(r8) :: asm_prm_dst1(numrad_snw)
  real(r8) :: ext_cff_mss_dst1(numrad_snw)

  ! dust species 2:
  real(r8) :: ss_alb_dst2(numrad_snw)
  real(r8) :: asm_prm_dst2(numrad_snw)
  real(r8) :: ext_cff_mss_dst2(numrad_snw)

  ! dust species 3:
  real(r8) :: ss_alb_dst3(numrad_snw)
  real(r8) :: asm_prm_dst3(numrad_snw)
  real(r8) :: ext_cff_mss_dst3(numrad_snw)

  ! dust species 4:
  real(r8) :: ss_alb_dst4(numrad_snw)
  real(r8) :: asm_prm_dst4(numrad_snw)
  real(r8) :: ext_cff_mss_dst4(numrad_snw)

  ! best-fit parameters for snow aging defined over:
  !  11 temperatures from 225 to 273 K
  !  31 temperature gradients from 0 to 300 K/m
  !   8 snow densities from 0 to 350 kg/m3
  ! (arrays declared here, but are set in iniTimeConst)
  real(r8), pointer :: snowage_tau(:,:,:) ! (idx_rhos_max,idx_Tgrd_max,idx_T_max)
  real(r8), pointer :: snowage_kappa(:,:,:) ! (idx_rhos_max,idx_Tgrd_max,idx_T_max)
  real(r8), pointer :: snowage_drdt0(:,:,:) ! idx_rhos_max,idx_Tgrd_max,idx_T_max)
  !
  ! !REVISION HISTORY:
  ! Created by Mark Flanner

  character(len=*), parameter, private :: sourcefile = &
       __FILE__
  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
   subroutine SnowOptics_init( )
     
     use fileutils  , only : getfil
     use CLM_varctl , only : fsnowoptics
     use spmdMod    , only : masterproc
     use ncdio_pio  , only : file_desc_t, ncd_io, ncd_pio_openfile, ncd_pio_closefile

     type(file_desc_t)  :: ncid                        ! netCDF file id
     character(len=256) :: locfn                       ! local filename
     character(len= 32) :: subname = 'SnowOptics_init' ! subroutine name
     integer            :: ier                         ! error status

     return    ! return early
     !
     ! Open optics file:
     if(masterproc) write(iulog,*) 'Attempting to read snow optical properties .....'
     call getfil (fsnowoptics, locfn, 0)
     call ncd_pio_openfile(ncid, locfn, 0)
     if(masterproc) write(iulog,*) subname,trim(fsnowoptics)

     ! direct-beam snow Mie parameters:
     call ncd_io('ss_alb_ice_drc', ss_alb_snw_drc,            'read', ncid, posNOTonfile=.true.)
     call ncd_io( 'asm_prm_ice_drc',asm_prm_snw_drc,          'read', ncid, posNOTonfile=.true.)
     call ncd_io( 'ext_cff_mss_ice_drc', ext_cff_mss_snw_drc, 'read', ncid, posNOTonfile=.true.)

     ! diffuse snow Mie parameters
     call ncd_io( 'ss_alb_ice_dfs', ss_alb_snw_dfs,           'read', ncid, posNOTonfile=.true.)
     call ncd_io( 'asm_prm_ice_dfs', asm_prm_snw_dfs,         'read', ncid, posNOTonfile=.true.)
     call ncd_io( 'ext_cff_mss_ice_dfs', ext_cff_mss_snw_dfs, 'read', ncid, posNOTonfile=.true.)

     ! BC species 1 Mie parameters
     call ncd_io( 'ss_alb_bcphil', ss_alb_bc1,           'read', ncid, posNOTonfile=.true.)
     call ncd_io( 'asm_prm_bcphil', asm_prm_bc1,         'read', ncid, posNOTonfile=.true.)
     call ncd_io( 'ext_cff_mss_bcphil', ext_cff_mss_bc1, 'read', ncid, posNOTonfile=.true.)

     ! BC species 2 Mie parameters
     call ncd_io( 'ss_alb_bcphob', ss_alb_bc2,           'read', ncid, posNOTonfile=.true.)
     call ncd_io( 'asm_prm_bcphob', asm_prm_bc2,         'read', ncid, posNOTonfile=.true.)
     call ncd_io( 'ext_cff_mss_bcphob', ext_cff_mss_bc2, 'read', ncid, posNOTonfile=.true.)

     ! OC species 1 Mie parameters
     call ncd_io( 'ss_alb_ocphil', ss_alb_oc1,           'read', ncid, posNOTonfile=.true.)
     call ncd_io( 'asm_prm_ocphil', asm_prm_oc1,         'read', ncid, posNOTonfile=.true.)
     call ncd_io( 'ext_cff_mss_ocphil', ext_cff_mss_oc1, 'read', ncid, posNOTonfile=.true.)

     ! OC species 2 Mie parameters
     call ncd_io( 'ss_alb_ocphob', ss_alb_oc2,           'read', ncid, posNOTonfile=.true.)
     call ncd_io( 'asm_prm_ocphob', asm_prm_oc2,         'read', ncid, posNOTonfile=.true.)
     call ncd_io( 'ext_cff_mss_ocphob', ext_cff_mss_oc2, 'read', ncid, posNOTonfile=.true.)

     ! dust species 1 Mie parameters
     call ncd_io( 'ss_alb_dust01', ss_alb_dst1,           'read', ncid, posNOTonfile=.true.)
     call ncd_io( 'asm_prm_dust01', asm_prm_dst1,         'read', ncid, posNOTonfile=.true.)
     call ncd_io( 'ext_cff_mss_dust01', ext_cff_mss_dst1, 'read', ncid, posNOTonfile=.true.)

     ! dust species 2 Mie parameters
     call ncd_io( 'ss_alb_dust02', ss_alb_dst2,           'read', ncid, posNOTonfile=.true.)
     call ncd_io( 'asm_prm_dust02', asm_prm_dst2,         'read', ncid, posNOTonfile=.true.)
     call ncd_io( 'ext_cff_mss_dust02', ext_cff_mss_dst2, 'read', ncid, posNOTonfile=.true.)

     ! dust species 3 Mie parameters
     call ncd_io( 'ss_alb_dust03', ss_alb_dst3,           'read', ncid, posNOTonfile=.true.)
     call ncd_io( 'asm_prm_dust03', asm_prm_dst3,         'read', ncid, posNOTonfile=.true.)
     call ncd_io( 'ext_cff_mss_dust03', ext_cff_mss_dst3, 'read', ncid, posNOTonfile=.true.)

     ! dust species 4 Mie parameters
     call ncd_io( 'ss_alb_dust04', ss_alb_dst4,           'read', ncid, posNOTonfile=.true.)
     call ncd_io( 'asm_prm_dust04', asm_prm_dst4,         'read', ncid, posNOTonfile=.true.)
     call ncd_io( 'ext_cff_mss_dust04', ext_cff_mss_dst4, 'read', ncid, posNOTonfile=.true.)


     call ncd_pio_closefile(ncid)
     if (masterproc) then

        write(iulog,*) 'Successfully read snow optical properties'
        ! print some diagnostics:
        write (iulog,*) 'SNICAR: Mie single scatter albedos for direct-beam ice, rds=100um: ', &
             ss_alb_snw_drc(71,1), ss_alb_snw_drc(71,2), ss_alb_snw_drc(71,3),     &
             ss_alb_snw_drc(71,4), ss_alb_snw_drc(71,5)
        write (iulog,*) 'SNICAR: Mie single scatter albedos for diffuse ice, rds=100um: ',     &
             ss_alb_snw_dfs(71,1), ss_alb_snw_dfs(71,2), ss_alb_snw_dfs(71,3),     &
             ss_alb_snw_dfs(71,4), ss_alb_snw_dfs(71,5)
        if (DO_SNO_OC) then
           write (iulog,*) 'SNICAR: Including OC aerosols from snow radiative transfer calculations'
        else
           write (iulog,*) 'SNICAR: Excluding OC aerosols from snow radiative transfer calculations'
        endif
        write (iulog,*) 'SNICAR: Mie single scatter albedos for hydrophillic BC: ', &
             ss_alb_bc1(1), ss_alb_bc1(2), ss_alb_bc1(3), ss_alb_bc1(4), ss_alb_bc1(5)
        write (iulog,*) 'SNICAR: Mie single scatter albedos for hydrophobic BC: ', &
             ss_alb_bc2(1), ss_alb_bc2(2), ss_alb_bc2(3), ss_alb_bc2(4), ss_alb_bc2(5)
        if (DO_SNO_OC) then
           write (iulog,*) 'SNICAR: Mie single scatter albedos for hydrophillic OC: ', &
                ss_alb_oc1(1), ss_alb_oc1(2), ss_alb_oc1(3), ss_alb_oc1(4), ss_alb_oc1(5)
           write (iulog,*) 'SNICAR: Mie single scatter albedos for hydrophobic OC: ', &
                ss_alb_oc2(1), ss_alb_oc2(2), ss_alb_oc2(3), ss_alb_oc2(4), ss_alb_oc2(5)
        endif
        write (iulog,*) 'SNICAR: Mie single scatter albedos for dust species 1: ', &
             ss_alb_dst1(1), ss_alb_dst1(2), ss_alb_dst1(3), ss_alb_dst1(4), ss_alb_dst1(5)
        write (iulog,*) 'SNICAR: Mie single scatter albedos for dust species 2: ', &
             ss_alb_dst2(1), ss_alb_dst2(2), ss_alb_dst2(3), ss_alb_dst2(4), ss_alb_dst2(5)
        write (iulog,*) 'SNICAR: Mie single scatter albedos for dust species 3: ', &
             ss_alb_dst3(1), ss_alb_dst3(2), ss_alb_dst3(3), ss_alb_dst3(4), ss_alb_dst3(5)
        write (iulog,*) 'SNICAR: Mie single scatter albedos for dust species 4: ', &
             ss_alb_dst4(1), ss_alb_dst4(2), ss_alb_dst4(3), ss_alb_dst4(4), ss_alb_dst4(5)
        write(iulog,*)
     end if

   end subroutine SnowOptics_init

   !-----------------------------------------------------------------------
   subroutine SnowAge_init( )
     use CLM_varctl      , only : fsnowaging
     use fileutils       , only : getfil
     use spmdMod         , only : masterproc
     use ncdio_pio       , only : file_desc_t, ncd_io, ncd_pio_openfile, ncd_pio_closefile

     type(file_desc_t)  :: ncid                        ! netCDF file id
     character(len=256) :: locfn                       ! local filename
     character(len= 32) :: subname = 'SnowOptics_init' ! subroutine name
     integer            :: varid                       ! netCDF id's
     integer            :: ier                         ! error status

     ! Open snow aging (effective radius evolution) file:
     allocate(snowage_tau(idx_rhos_max,idx_Tgrd_max,idx_T_max))
     allocate(snowage_kappa(idx_rhos_max,idx_Tgrd_max,idx_T_max))
     allocate(snowage_drdt0(idx_rhos_max,idx_Tgrd_max,idx_T_max))

     return    ! return early
     if(masterproc)  write(iulog,*) 'Attempting to read snow aging parameters .....'
     call getfil (fsnowaging, locfn, 0)
     call ncd_pio_openfile(ncid, locfn, 0)
     if(masterproc) write(iulog,*) subname,trim(fsnowaging)

     ! snow aging parameters

     call ncd_io('tau', snowage_tau,       'read', ncid, posNOTonfile=.true.)
     call ncd_io('kappa', snowage_kappa,   'read', ncid, posNOTonfile=.true.)
     call ncd_io('drdsdt0', snowage_drdt0, 'read', ncid, posNOTonfile=.true.)

     call ncd_pio_closefile(ncid)
     if (masterproc) then

        write(iulog,*) 'Successfully read snow aging properties'

        ! print some diagnostics:
        write (iulog,*) 'SNICAR: snowage tau for T=263K, dTdz = 100 K/m, rhos = 150 kg/m3: ', snowage_tau(3,11,9)
        write (iulog,*) 'SNICAR: snowage kappa for T=263K, dTdz = 100 K/m, rhos = 150 kg/m3: ', snowage_kappa(3,11,9)
        write (iulog,*) 'SNICAR: snowage dr/dt_0 for T=263K, dTdz = 100 K/m, rhos = 150 kg/m3: ', snowage_drdt0(3,11,9)
     endif

   end subroutine SnowAge_init
   
 end module SnowSnicarMod
