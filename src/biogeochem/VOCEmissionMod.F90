module VOCEmissionMod

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Volatile organic compound emission
  !
  ! !USES:
  use shr_kind_mod       , only : r8 => shr_kind_r8
  use shr_log_mod        , only : errMsg => shr_log_errMsg
  use clm_varctl         , only : iulog
  use clm_varpar         , only : numpft, nlevcan
  use pftconMod          , only : ndllf_evr_tmp_tree,  ndllf_evr_brl_tree
  use pftconMod          , only : ndllf_dcd_brl_tree,  nbrdlf_evr_trp_tree
  use pftconMod          , only : nbrdlf_evr_tmp_tree, nbrdlf_dcd_brl_shrub
  use pftconMod          , only : nbrdlf_dcd_trp_tree, nbrdlf_dcd_tmp_tree
  use pftconMod          , only : nbrdlf_dcd_brl_tree, nbrdlf_evr_shrub
  use pftconMod          , only : nc3_arctic_grass   , nc3crop
  use pftconMod          , only : nc4_grass,           noveg
  use shr_megan_mod      , only : shr_megan_megcomps_n, shr_megan_megcomp_t, shr_megan_linkedlist
  use shr_megan_mod      , only : shr_megan_mechcomps_n, shr_megan_mechcomps, shr_megan_mapped_emisfctrs
  use MEGANFactorsMod    , only : Agro, Amat, Anew, Aold, betaT, ct1, ct2, LDF, Ceo
  use decompMod          , only : bounds_type
  use abortutils         , only : endrun
  use fileutils          , only : getfil
  use clm_varcon         , only : grlnd
  use atm2lndType        , only : atm2lnd_type
  use CanopyStateType    , only : canopystate_type
  use PhotosynthesisMod  , only : photosyns_type
  use SoilStateType      , only : soilstate_type
  use SolarAbsorbedType  , only : solarabs_type
  use TemperatureType    , only : temperature_type
  use PatchType          , only : patch                
  !
  implicit none
  private 
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  !
  ! !PUBLIC TYPES:
  type, public :: vocemis_type
     real(r8) , pointer, private :: Eopt_out_patch    (:)   ! Eopt coefficient
     real(r8) , pointer, private :: topt_out_patch    (:)   ! topt coefficient
     real(r8) , pointer, private :: alpha_out_patch   (:)   ! alpha coefficient
     real(r8) , pointer, private :: cp_out_patch      (:)   ! cp coefficient
     real(r8) , pointer, private :: paru_out_patch    (:)   ! 
     real(r8) , pointer, private :: par24u_out_patch  (:)   ! 
     real(r8) , pointer, private :: par240u_out_patch (:)   !  
     real(r8) , pointer, private :: para_out_patch    (:)   ! 
     real(r8) , pointer, private :: par24a_out_patch  (:)   ! 
     real(r8) , pointer, private :: par240a_out_patch (:)   ! 
     real(r8) , pointer, private :: gamma_out_patch   (:)   ! 
     real(r8) , pointer, private :: gammaL_out_patch  (:)   ! 
     real(r8) , pointer, private :: gammaT_out_patch  (:)   ! 
     real(r8) , pointer, private :: gammaP_out_patch  (:)   ! 
     real(r8) , pointer, private :: gammaA_out_patch  (:)   ! 
     real(r8) , pointer, private :: gammaS_out_patch  (:)   ! 
     real(r8) , pointer, private :: gammaC_out_patch  (:)   ! 
     real(r8) , pointer, private :: vocflx_tot_patch  (:)   ! total VOC flux into atmosphere [moles/m2/sec] 
     real(r8) , pointer, PUBLIC  :: vocflx_patch      (:,:) ! (num_mech_comps) MEGAN flux [moles/m2/sec] 
     real(r8) , pointer, private :: efisop_grc        (:,:) ! gridcell isoprene emission factors
  end type vocemis_type
  !
  ! !PRIVATE TYPES:
  type :: megan_out_type
     ! VOC fluxes structure for CLM history output
     real(r8), pointer, private  :: flux_out(:)   ! patch MEGAN flux [ug C m-2 h-1]
  end type megan_out_type
  type(megan_out_type), private, pointer :: meg_out(:) ! (n_megan_comps) points to output fluxes
  !
  logical, parameter :: debug = .false.

  character(len=*), parameter, private :: sourcefile = &
       __FILE__
  !------------------------------------------------------------------------

end module VOCEmissionMod


