module readParamsMod

  !-----------------------------------------------------------------------
  !
  ! Read parameters
  ! module used to read parameters for individual modules and/or for some 
  ! well defined functionality (eg. ED).
  !
  ! ! USES:
  use clm_varctl , only : paramfile, iulog, use_fates, use_cn
  use spmdMod    , only : masterproc
  use fileutils  , only : getfil
  use ncdio_pio  , only : ncd_pio_closefile, ncd_pio_openfile
  use ncdio_pio  , only : file_desc_t , ncd_inqdid, ncd_inqdlen

  implicit none
  private
  !
  public :: readParameters

  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine readParameters (photosyns_inst)
    !
    ! ! USES:
    use CNSharedParamsMod                 , only : CNParamsReadShared
    use CNGapMortalityMod                 , only : readCNGapMortParams                    => readParams
    use CNMRespMod                        , only : readCNMRespParams                      => readParams
    use CNPhenologyMod                    , only : readCNPhenolParams                     => readParams
    use SoilBiogeochemCompetitionMod      , only : readSoilBiogeochemCompetitionParams    => readParams
    use SoilBiogeochemNLeachingMod        , only : readSoilBiogeochemNLeachingParams      => readParams
    use SoilBiogeochemNitrifDenitrifMod   , only : readSoilBiogeochemNitrifDenitrifParams => readParams
    use SoilBiogeochemLittVertTranspMod   , only : readSoilBiogeochemLittVertTranspParams => readParams
    use SoilBiogeochemPotentialMod        , only : readSoilBiogeochemPotentialParams      => readParams
    use SoilBiogeochemDecompMod           , only : readSoilBiogeochemDecompParams         => readParams
    use SoilBiogeochemDecompCascadeBGCMod , only : readSoilBiogeochemDecompBgcParams      => readParams
    use SoilBiogeochemDecompCascadeCNMod  , only : readSoilBiogeochemDecompCnParams       => readParams
    !use ch4Mod                            , only : readCH4Params                          => readParams
    use clm_varctl,                         only : NLFilename_in
    use PhotosynthesisMod                 , only : photosyns_type
    !
    ! !ARGUMENTS:
    type(photosyns_type)                   , intent(in) :: photosyns_inst
    !
    ! !LOCAL VARIABLES:
    character(len=256) :: locfn ! local file name
    type(file_desc_t)  :: ncid  ! pio netCDF file id
    integer            :: dimid ! netCDF dimension id
    integer            :: npft  ! number of pfts on pft-physiology file
    character(len=32)  :: subname = 'readParameters'
    !-----------------------------------------------------------------------

    if (masterproc) then
       write(iulog,*) 'paramMod.F90::'//trim(subname)//' :: reading CLM '//' parameters '
    end if

    call getfil (paramfile, locfn, 0)
    call ncd_pio_openfile (ncid, trim(locfn), 0)
    call ncd_inqdid(ncid,'pft',dimid) 
    call ncd_inqdlen(ncid,dimid,npft) 

    !
    ! Above ground biogeochemistry...
    !
    if (use_cn) then
       call readCNGapMortParams(ncid)
       call readCNMRespParams(ncid)
       call readCNPhenolParams(ncid)
    end if

    !
    ! Soil biogeochemistry...
    !
    if (use_cn .or. use_fates) then
       call readSoilBiogeochemCompetitionParams(ncid)
       call readSoilBiogeochemDecompBgcParams(ncid)
       call readSoilBiogeochemDecompCnParams(ncid)
       call readSoilBiogeochemDecompParams(ncid)
       call readSoilBiogeochemLittVertTranspParams(ncid)
       call readSoilBiogeochemNitrifDenitrifParams(ncid)
       call readSoilBiogeochemNLeachingParams(ncid)
       call readSoilBiogeochemPotentialParams(ncid)
       call CNParamsReadShared(ncid, NLFilename_in)  ! this is called CN params but really is for the soil biogeochem parameters

       !call readCH4Params (ncid)
    end if

    !
    ! Biogeophysics
    !
    call photosyns_inst%ReadParams( ncid )


    !
    call ncd_pio_closefile(ncid)

  end subroutine readParameters

end module readParamsMod
