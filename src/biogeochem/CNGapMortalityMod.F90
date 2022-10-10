module CNGapMortalityMod

#include "shr_assert.h"

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Module holding routines used in gap mortality for coupled carbon
  ! nitrogen code.
  !
  ! !USES:
  use shr_kind_mod            , only : r8 => shr_kind_r8
  use decompMod               , only : bounds_type
  use abortutils              , only : endrun
  use shr_log_mod             , only : errMsg => shr_log_errMsg
  use pftconMod               , only : pftcon
  use CNVegCarbonStateType    , only : cnveg_carbonstate_type
  use CNVegCarbonFluxType     , only : cnveg_carbonflux_type
  use CNVegNitrogenStateType  , only : cnveg_nitrogenstate_type
  use CNVegNitrogenFluxType   , only : cnveg_nitrogenflux_type
  use CanopyStateType         , only : canopystate_type            
  use ColumnType              , only : col                
  use PatchType               , only : patch                
  !
  implicit none
  private
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: readParams

  type, private :: params_type
     real(r8):: am     ! mortality rate based on annual rate, fractional mortality (1/yr)
     real(r8):: k_mort ! coeff. of growth efficiency in mortality equation
  end type params_type
  !
  type(params_type), private :: params_inst
  !
  ! !PRIVATE MEMBER FUNCTIONS:
  private :: CNGap_PatchToColumn

  character(len=*), parameter, private :: sourcefile = &
       __FILE__
  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine readParams ( ncid )
    !
    ! !DESCRIPTION:
    ! Read in parameters
    !
    ! !USES:
    use ncdio_pio  , only : file_desc_t,ncd_io
    !
    ! !ARGUMENTS:
    implicit none
    type(file_desc_t),intent(inout) :: ncid   ! pio netCDF file id
    !
    ! !LOCAL VARIABLES:
    character(len=32)  :: subname = 'CNGapMortParamsType'
    character(len=100) :: errCode = '-Error reading in parameters file:'
    logical            :: readv ! has variable been read in or not
    real(r8)           :: tempr ! temporary to read in constant
    character(len=100) :: tString ! temp. var for reading
    !-----------------------------------------------------------------------

    tString='r_mort'
    call ncd_io(varname=trim(tString),data=tempr, flag='read', ncid=ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))
    params_inst%am=tempr

    tString='k_mort'
    call ncd_io(varname=trim(tString),data=tempr, flag='read', ncid=ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))
    params_inst%k_mort=tempr   
    
  end subroutine readParams

  !-----------------------------------------------------------------------
  subroutine CNGap_PatchToColumn (bounds, num_soilc, filter_soilc, &
       cnveg_carbonflux_inst, cnveg_nitrogenflux_inst, &
       leaf_prof_patch, froot_prof_patch, croot_prof_patch, stem_prof_patch)
    !
    ! !DESCRIPTION:
    ! gathers all patch-level gap mortality fluxes to the column level and
    ! assigns them to the three litter pools
    !
    ! !USES:
    use clm_varpar , only : maxpatch_pft, nlevdecomp, nlevdecomp_full
    !
    ! !ARGUMENTS:
    type(bounds_type)               , intent(in)    :: bounds
    integer                         , intent(in)    :: num_soilc       ! number of soil columns in filter
    integer                         , intent(in)    :: filter_soilc(:) ! soil column filter
    type(cnveg_carbonflux_type)     , intent(inout) :: cnveg_carbonflux_inst
    type(cnveg_nitrogenflux_type)   , intent(inout) :: cnveg_nitrogenflux_inst
    real(r8)                        , intent(in)    :: leaf_prof_patch(bounds%begp:,1:)
    real(r8)                        , intent(in)    :: froot_prof_patch(bounds%begp:,1:)
    real(r8)                        , intent(in)    :: croot_prof_patch(bounds%begp:,1:)
    real(r8)                        , intent(in)    :: stem_prof_patch(bounds%begp:,1:)
    !
    ! !LOCAL VARIABLES:
    integer :: fc,c,pi,p,j               ! indices
    !-----------------------------------------------------------------------

    SHR_ASSERT_ALL((ubound(leaf_prof_patch)   == (/bounds%endp,nlevdecomp_full/)), errMsg(sourcefile, __LINE__))
    SHR_ASSERT_ALL((ubound(froot_prof_patch)  == (/bounds%endp,nlevdecomp_full/)), errMsg(sourcefile, __LINE__))
    SHR_ASSERT_ALL((ubound(croot_prof_patch)  == (/bounds%endp,nlevdecomp_full/)), errMsg(sourcefile, __LINE__))
    SHR_ASSERT_ALL((ubound(stem_prof_patch)   == (/bounds%endp,nlevdecomp_full/)), errMsg(sourcefile, __LINE__))

    associate(                                                                                                 & 
         leaf_prof                           => leaf_prof_patch                                              , & ! Input:  [real(r8) (:,:) ]  (1/m) profile of leaves                         
         froot_prof                          => froot_prof_patch                                             , & ! Input:  [real(r8) (:,:) ]  (1/m) profile of fine roots                     
         croot_prof                          => croot_prof_patch                                             , & ! Input:  [real(r8) (:,:) ]  (1/m) profile of coarse roots                   
         stem_prof                           => stem_prof_patch                                              , & ! Input:  [real(r8) (:,:) ]  (1/m) profile of stems                          

         ivt                                 => patch%itype                                                    , & ! Input:  [integer  (:)   ]  patch vegetation type                                
         wtcol                               => patch%wtcol                                                    , & ! Input:  [real(r8) (:)   ]  patch weight relative to column (0-1)               
         
         lf_flab                             => pftcon%lf_flab                                               , & ! Input:  [real(r8) (:)   ]  leaf litter labile fraction                       
         lf_fcel                             => pftcon%lf_fcel                                               , & ! Input:  [real(r8) (:)   ]  leaf litter cellulose fraction                    
         lf_flig                             => pftcon%lf_flig                                               , & ! Input:  [real(r8) (:)   ]  leaf litter lignin fraction                       
         fr_flab                             => pftcon%fr_flab                                               , & ! Input:  [real(r8) (:)   ]  fine root litter labile fraction                  
         fr_fcel                             => pftcon%fr_fcel                                               , & ! Input:  [real(r8) (:)   ]  fine root litter cellulose fraction               
         fr_flig                             => pftcon%fr_flig                                               , & ! Input:  [real(r8) (:)   ]  fine root litter lignin fraction                  
         
         m_leafc_to_litter                   => cnveg_carbonflux_inst%m_leafc_to_litter_patch                , & ! Input:  [real(r8) (:)   ]                                                    
         m_frootc_to_litter                  => cnveg_carbonflux_inst%m_frootc_to_litter_patch               , & ! Input:  [real(r8) (:)   ]                                                    
         m_livestemc_to_litter               => cnveg_carbonflux_inst%m_livestemc_to_litter_patch            , & ! Input:  [real(r8) (:)   ]                                                    
         m_deadstemc_to_litter               => cnveg_carbonflux_inst%m_deadstemc_to_litter_patch            , & ! Input:  [real(r8) (:)   ]                                                    
         m_livecrootc_to_litter              => cnveg_carbonflux_inst%m_livecrootc_to_litter_patch           , & ! Input:  [real(r8) (:)   ]                                                    
         m_deadcrootc_to_litter              => cnveg_carbonflux_inst%m_deadcrootc_to_litter_patch           , & ! Input:  [real(r8) (:)   ]                                                    
         m_leafc_storage_to_litter           => cnveg_carbonflux_inst%m_leafc_storage_to_litter_patch        , & ! Input:  [real(r8) (:)   ]                                                    
         m_frootc_storage_to_litter          => cnveg_carbonflux_inst%m_frootc_storage_to_litter_patch       , & ! Input:  [real(r8) (:)   ]                                                    
         m_livestemc_storage_to_litter       => cnveg_carbonflux_inst%m_livestemc_storage_to_litter_patch    , & ! Input:  [real(r8) (:)   ]                                                    
         m_deadstemc_storage_to_litter       => cnveg_carbonflux_inst%m_deadstemc_storage_to_litter_patch    , & ! Input:  [real(r8) (:)   ]                                                    
         m_livecrootc_storage_to_litter      => cnveg_carbonflux_inst%m_livecrootc_storage_to_litter_patch   , & ! Input:  [real(r8) (:)   ]                                                    
         m_deadcrootc_storage_to_litter      => cnveg_carbonflux_inst%m_deadcrootc_storage_to_litter_patch   , & ! Input:  [real(r8) (:)   ]                                                    
         m_gresp_storage_to_litter           => cnveg_carbonflux_inst%m_gresp_storage_to_litter_patch        , & ! Input:  [real(r8) (:)   ]                                                    
         m_leafc_xfer_to_litter              => cnveg_carbonflux_inst%m_leafc_xfer_to_litter_patch           , & ! Input:  [real(r8) (:)   ]                                                    
         m_frootc_xfer_to_litter             => cnveg_carbonflux_inst%m_frootc_xfer_to_litter_patch          , & ! Input:  [real(r8) (:)   ]                                                    
         m_livestemc_xfer_to_litter          => cnveg_carbonflux_inst%m_livestemc_xfer_to_litter_patch       , & ! Input:  [real(r8) (:)   ]                                                    
         m_deadstemc_xfer_to_litter          => cnveg_carbonflux_inst%m_deadstemc_xfer_to_litter_patch       , & ! Input:  [real(r8) (:)   ]                                                    
         m_livecrootc_xfer_to_litter         => cnveg_carbonflux_inst%m_livecrootc_xfer_to_litter_patch      , & ! Input:  [real(r8) (:)   ]                                                    
         m_deadcrootc_xfer_to_litter         => cnveg_carbonflux_inst%m_deadcrootc_xfer_to_litter_patch      , & ! Input:  [real(r8) (:)   ]                                                    
         m_gresp_xfer_to_litter              => cnveg_carbonflux_inst%m_gresp_xfer_to_litter_patch           , & ! Input:  [real(r8) (:)   ]                                                    
         gap_mortality_c_to_litr_met_c       => cnveg_carbonflux_inst%gap_mortality_c_to_litr_met_c_col      , & ! Output: [real(r8) (:,:) ]  C fluxes associated with gap mortality to litter metabolic pool (gC/m3/s)
         gap_mortality_c_to_litr_cel_c       => cnveg_carbonflux_inst%gap_mortality_c_to_litr_cel_c_col      , & ! Output: [real(r8) (:,:) ]  C fluxes associated with gap mortality to litter cellulose pool (gC/m3/s)
         gap_mortality_c_to_litr_lig_c       => cnveg_carbonflux_inst%gap_mortality_c_to_litr_lig_c_col      , & ! Output: [real(r8) (:,:) ]  C fluxes associated with gap mortality to litter lignin pool (gC/m3/s)
         gap_mortality_c_to_cwdc             => cnveg_carbonflux_inst%gap_mortality_c_to_cwdc_col            , & ! Output: [real(r8) (:,:) ]  C fluxes associated with gap mortality to CWD pool (gC/m3/s)
         
         m_leafn_to_litter                   => cnveg_nitrogenflux_inst%m_leafn_to_litter_patch              , & ! Input:  [real(r8) (:)   ]                                                    
         m_frootn_to_litter                  => cnveg_nitrogenflux_inst%m_frootn_to_litter_patch             , & ! Input:  [real(r8) (:)   ]                                                    
         m_livestemn_to_litter               => cnveg_nitrogenflux_inst%m_livestemn_to_litter_patch          , & ! Input:  [real(r8) (:)   ]                                                    
         m_deadstemn_to_litter               => cnveg_nitrogenflux_inst%m_deadstemn_to_litter_patch          , & ! Input:  [real(r8) (:)   ]                                                    
         m_livecrootn_to_litter              => cnveg_nitrogenflux_inst%m_livecrootn_to_litter_patch         , & ! Input:  [real(r8) (:)   ]                                                    
         m_deadcrootn_to_litter              => cnveg_nitrogenflux_inst%m_deadcrootn_to_litter_patch         , & ! Input:  [real(r8) (:)   ]                                                    
         m_retransn_to_litter                => cnveg_nitrogenflux_inst%m_retransn_to_litter_patch           , & ! Input:  [real(r8) (:)   ]                                                    
         m_leafn_storage_to_litter           => cnveg_nitrogenflux_inst%m_leafn_storage_to_litter_patch      , & ! Input:  [real(r8) (:)   ]                                                    
         m_frootn_storage_to_litter          => cnveg_nitrogenflux_inst%m_frootn_storage_to_litter_patch     , & ! Input:  [real(r8) (:)   ]                                                    
         m_livestemn_storage_to_litter       => cnveg_nitrogenflux_inst%m_livestemn_storage_to_litter_patch  , & ! Input:  [real(r8) (:)   ]                                                    
         m_deadstemn_storage_to_litter       => cnveg_nitrogenflux_inst%m_deadstemn_storage_to_litter_patch  , & ! Input:  [real(r8) (:)   ]                                                    
         m_livecrootn_storage_to_litter      => cnveg_nitrogenflux_inst%m_livecrootn_storage_to_litter_patch , & ! Input:  [real(r8) (:)   ]                                                    
         m_deadcrootn_storage_to_litter      => cnveg_nitrogenflux_inst%m_deadcrootn_storage_to_litter_patch , & ! Input:  [real(r8) (:)   ]                                                    
         m_leafn_xfer_to_litter              => cnveg_nitrogenflux_inst%m_leafn_xfer_to_litter_patch         , & ! Input:  [real(r8) (:)   ]                                                    
         m_frootn_xfer_to_litter             => cnveg_nitrogenflux_inst%m_frootn_xfer_to_litter_patch        , & ! Input:  [real(r8) (:)   ]                                                    
         m_livestemn_xfer_to_litter          => cnveg_nitrogenflux_inst%m_livestemn_xfer_to_litter_patch     , & ! Input:  [real(r8) (:)   ]                                                    
         m_deadstemn_xfer_to_litter          => cnveg_nitrogenflux_inst%m_deadstemn_xfer_to_litter_patch     , & ! Input:  [real(r8) (:)   ]                                                    
         m_livecrootn_xfer_to_litter         => cnveg_nitrogenflux_inst%m_livecrootn_xfer_to_litter_patch    , & ! Input:  [real(r8) (:)   ]                                                    
         m_deadcrootn_xfer_to_litter         => cnveg_nitrogenflux_inst%m_deadcrootn_xfer_to_litter_patch    , & ! Input:  [real(r8) (:)   ]                                                    
         gap_mortality_n_to_litr_met_n       => cnveg_nitrogenflux_inst%gap_mortality_n_to_litr_met_n_col    , & ! Output: [real(r8) (:,:) ]  N fluxes associated with gap mortality to litter metabolic pool (gN/m3/s)
         gap_mortality_n_to_litr_cel_n       => cnveg_nitrogenflux_inst%gap_mortality_n_to_litr_cel_n_col    , & ! Output: [real(r8) (:,:) ]  N fluxes associated with gap mortality to litter cellulose pool (gN/m3/s)
         gap_mortality_n_to_litr_lig_n       => cnveg_nitrogenflux_inst%gap_mortality_n_to_litr_lig_n_col    , & ! Output: [real(r8) (:,:) ]  N fluxes associated with gap mortality to litter lignin pool (gN/m3/s)
         gap_mortality_n_to_cwdn             => cnveg_nitrogenflux_inst%gap_mortality_n_to_cwdn_col            & ! Output: [real(r8) (:,:) ]  N fluxes associated with gap mortality to CWD pool (gN/m3/s)
         )

      do j = 1,nlevdecomp
         do pi = 1,maxpatch_pft
            do fc = 1,num_soilc
               c = filter_soilc(fc)

               if (pi <=  col%npatches(c)) then
                  p = col%patchi(c) + pi - 1

                  if (patch%active(p)) then

                     ! leaf gap mortality carbon fluxes
                     gap_mortality_c_to_litr_met_c(c,j) = gap_mortality_c_to_litr_met_c(c,j) + &
                          m_leafc_to_litter(p) * lf_flab(ivt(p)) * wtcol(p) * leaf_prof(p,j)
                     gap_mortality_c_to_litr_cel_c(c,j) = gap_mortality_c_to_litr_cel_c(c,j) + &
                          m_leafc_to_litter(p) * lf_fcel(ivt(p)) * wtcol(p) * leaf_prof(p,j)
                     gap_mortality_c_to_litr_lig_c(c,j) = gap_mortality_c_to_litr_lig_c(c,j) + &
                          m_leafc_to_litter(p) * lf_flig(ivt(p)) * wtcol(p) * leaf_prof(p,j)

                     ! fine root gap mortality carbon fluxes
                     gap_mortality_c_to_litr_met_c(c,j) = gap_mortality_c_to_litr_met_c(c,j) + &
                          m_frootc_to_litter(p) * fr_flab(ivt(p)) * wtcol(p) * froot_prof(p,j)
                     gap_mortality_c_to_litr_cel_c(c,j) = gap_mortality_c_to_litr_cel_c(c,j) + &
                          m_frootc_to_litter(p) * fr_fcel(ivt(p)) * wtcol(p) * froot_prof(p,j)
                     gap_mortality_c_to_litr_lig_c(c,j) = gap_mortality_c_to_litr_lig_c(c,j) + &
                          m_frootc_to_litter(p) * fr_flig(ivt(p)) * wtcol(p) * froot_prof(p,j)

                     ! wood gap mortality carbon fluxes
                     gap_mortality_c_to_cwdc(c,j)  = gap_mortality_c_to_cwdc(c,j)  + &
                          (m_livestemc_to_litter(p) + m_deadstemc_to_litter(p))  * wtcol(p) * stem_prof(p,j)
                     gap_mortality_c_to_cwdc(c,j) = gap_mortality_c_to_cwdc(c,j) + &
                          (m_livecrootc_to_litter(p) + m_deadcrootc_to_litter(p)) * wtcol(p) * croot_prof(p,j)

                     ! storage gap mortality carbon fluxes
                     gap_mortality_c_to_litr_met_c(c,j) = gap_mortality_c_to_litr_met_c(c,j) + &
                          (m_leafc_storage_to_litter(p) + m_gresp_storage_to_litter(p)) * wtcol(p) * leaf_prof(p,j)
                     gap_mortality_c_to_litr_met_c(c,j) = gap_mortality_c_to_litr_met_c(c,j) + &
                          m_frootc_storage_to_litter(p) * wtcol(p) * froot_prof(p,j)
                     gap_mortality_c_to_litr_met_c(c,j) = gap_mortality_c_to_litr_met_c(c,j)  + &
                          (m_livestemc_storage_to_litter(p) + m_deadstemc_storage_to_litter(p)) * wtcol(p) * stem_prof(p,j)
                     gap_mortality_c_to_litr_met_c(c,j) = gap_mortality_c_to_litr_met_c(c,j) + &
                          (m_livecrootc_storage_to_litter(p) + m_deadcrootc_storage_to_litter(p)) * wtcol(p) * croot_prof(p,j)

                     ! transfer gap mortality carbon fluxes
                     gap_mortality_c_to_litr_met_c(c,j) = gap_mortality_c_to_litr_met_c(c,j) + &
                          (m_leafc_xfer_to_litter(p) + m_gresp_xfer_to_litter(p)) * wtcol(p) * leaf_prof(p,j)
                     gap_mortality_c_to_litr_met_c(c,j) = gap_mortality_c_to_litr_met_c(c,j) + &
                          m_frootc_xfer_to_litter(p) * wtcol(p) * froot_prof(p,j)
                     gap_mortality_c_to_litr_met_c(c,j)  = gap_mortality_c_to_litr_met_c(c,j)  + &
                          (m_livestemc_xfer_to_litter(p) + m_deadstemc_xfer_to_litter(p))  * wtcol(p) * stem_prof(p,j)
                     gap_mortality_c_to_litr_met_c(c,j) = gap_mortality_c_to_litr_met_c(c,j) + &
                          (m_livecrootc_xfer_to_litter(p) + m_deadcrootc_xfer_to_litter(p)) * wtcol(p) * croot_prof(p,j)

                     ! leaf gap mortality nitrogen fluxes
                     gap_mortality_n_to_litr_met_n(c,j) = gap_mortality_n_to_litr_met_n(c,j) + &
                          m_leafn_to_litter(p) * lf_flab(ivt(p)) * wtcol(p) * leaf_prof(p,j)
                     gap_mortality_n_to_litr_cel_n(c,j) = gap_mortality_n_to_litr_cel_n(c,j) + &
                          m_leafn_to_litter(p) * lf_fcel(ivt(p)) * wtcol(p) * leaf_prof(p,j)
                     gap_mortality_n_to_litr_lig_n(c,j) = gap_mortality_n_to_litr_lig_n(c,j) + &
                          m_leafn_to_litter(p) * lf_flig(ivt(p)) * wtcol(p) * leaf_prof(p,j)

                     ! fine root litter nitrogen fluxes
                     gap_mortality_n_to_litr_met_n(c,j) = gap_mortality_n_to_litr_met_n(c,j) + &
                          m_frootn_to_litter(p) * fr_flab(ivt(p)) * wtcol(p) * froot_prof(p,j)
                     gap_mortality_n_to_litr_cel_n(c,j) = gap_mortality_n_to_litr_cel_n(c,j) + &
                          m_frootn_to_litter(p) * fr_fcel(ivt(p)) * wtcol(p) * froot_prof(p,j)
                     gap_mortality_n_to_litr_lig_n(c,j) = gap_mortality_n_to_litr_lig_n(c,j) + &
                          m_frootn_to_litter(p) * fr_flig(ivt(p)) * wtcol(p) * froot_prof(p,j)

                     ! wood gap mortality nitrogen fluxes
                     gap_mortality_n_to_cwdn(c,j) = gap_mortality_n_to_cwdn(c,j)  + &
                          (m_livestemn_to_litter(p) + m_deadstemn_to_litter(p))  * wtcol(p) * stem_prof(p,j)
                     gap_mortality_n_to_cwdn(c,j) = gap_mortality_n_to_cwdn(c,j) + &
                          (m_livecrootn_to_litter(p) + m_deadcrootn_to_litter(p)) * wtcol(p) * croot_prof(p,j)

                     ! retranslocated N pool gap mortality fluxes
                     gap_mortality_n_to_litr_met_n(c,j) = gap_mortality_n_to_litr_met_n(c,j) + &
                          m_retransn_to_litter(p) * wtcol(p) * leaf_prof(p,j)

                     ! storage gap mortality nitrogen fluxes
                     gap_mortality_n_to_litr_met_n(c,j) = gap_mortality_n_to_litr_met_n(c,j) + &
                          m_leafn_storage_to_litter(p) * wtcol(p) * leaf_prof(p,j)
                     gap_mortality_n_to_litr_met_n(c,j) = gap_mortality_n_to_litr_met_n(c,j) + &
                          m_frootn_storage_to_litter(p) * wtcol(p) * froot_prof(p,j)
                     gap_mortality_n_to_litr_met_n(c,j)  = gap_mortality_n_to_litr_met_n(c,j) + &
                          (m_livestemn_storage_to_litter(p) + m_deadstemn_storage_to_litter(p))  * wtcol(p) * stem_prof(p,j)
                     gap_mortality_n_to_litr_met_n(c,j) = gap_mortality_n_to_litr_met_n(c,j) + &
                          (m_livecrootn_storage_to_litter(p) + m_deadcrootn_storage_to_litter(p)) * wtcol(p) * croot_prof(p,j)

                     ! transfer gap mortality nitrogen fluxes
                     gap_mortality_n_to_litr_met_n(c,j) = gap_mortality_n_to_litr_met_n(c,j) + &
                          m_leafn_xfer_to_litter(p) * wtcol(p) * leaf_prof(p,j)
                     gap_mortality_n_to_litr_met_n(c,j) = gap_mortality_n_to_litr_met_n(c,j) + &
                          m_frootn_xfer_to_litter(p) * wtcol(p) * froot_prof(p,j)
                     gap_mortality_n_to_litr_met_n(c,j) = gap_mortality_n_to_litr_met_n(c,j) + &
                          (m_livestemn_xfer_to_litter(p) + m_deadstemn_xfer_to_litter(p))  * wtcol(p) * stem_prof(p,j)
                     gap_mortality_n_to_litr_met_n(c,j) = gap_mortality_n_to_litr_met_n(c,j) + &
                          (m_livecrootn_xfer_to_litter(p) + m_deadcrootn_xfer_to_litter(p)) * wtcol(p) * croot_prof(p,j)


                  end if
               end if

            end do
         end do
      end do

    end associate

  end subroutine CNGap_PatchToColumn

end module CNGapMortalityMod
