module AerosolMod

#include "shr_assert.h"

  !-----------------------------------------------------------------------
  use shr_kind_mod     , only : r8 => shr_kind_r8
  use shr_log_mod      , only : errMsg => shr_log_errMsg
  use shr_infnan_mod   , only : nan => shr_infnan_nan, assignment(=)
  use decompMod        , only : bounds_type
  use clm_varpar       , only : nlevsno, nlevgrnd 
  use clm_time_manager , only : get_step_size
  use atm2lndType      , only : atm2lnd_type
  use WaterfluxType    , only : waterflux_type
  use WaterstateType   , only : waterstate_type
  use ColumnType       , only : col               
  use abortutils       , only : endrun
  !
  ! !PUBLIC TYPES:
  implicit none
  private
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  !
  ! !PUBLIC DATA MEMBERS:
  real(r8), public, parameter :: snw_rds_min = 54.526_r8          ! minimum allowed snow effective radius (also cold "fresh snow" value) [microns]
  real(r8), public            :: fresh_snw_rds_max = 204.526_r8   ! maximum warm fresh snow effective radius [microns]
  !
  type, public :: aerosol_type
     real(r8), pointer, public  :: mss_bcpho_col(:,:)      ! mass of hydrophobic BC in snow (col,lyr)     [kg]
     real(r8), pointer, public  :: mss_bcphi_col(:,:)      ! mass of hydrophillic BC in snow (col,lyr)    [kg]
     real(r8), pointer, public  :: mss_bctot_col(:,:)      ! total mass of BC in snow (pho+phi) (col,lyr) [kg]
     real(r8), pointer, public  :: mss_bc_col_col(:)       ! column-integrated mass of total BC           [kg]
     real(r8), pointer, public  :: mss_bc_top_col(:)       ! top-layer mass of total BC                   [kg]

     real(r8), pointer, public  :: mss_ocpho_col(:,:)      ! mass of hydrophobic OC in snow (col,lyr)     [kg]
     real(r8), pointer, public  :: mss_ocphi_col(:,:)      ! mass of hydrophillic OC in snow (col,lyr)    [kg]
     real(r8), pointer, public  :: mss_octot_col(:,:)      ! total mass of OC in snow (pho+phi) (col,lyr) [kg]
     real(r8), pointer, public  :: mss_oc_col_col(:)       ! column-integrated mass of total OC           [kg]
     real(r8), pointer, public  :: mss_oc_top_col(:)       ! top-layer mass of total OC                   [kg]

     real(r8), pointer, public  :: mss_dst1_col(:,:)       ! mass of dust species 1 in snow (col,lyr)     [kg]
     real(r8), pointer, public  :: mss_dst2_col(:,:)       ! mass of dust species 2 in snow (col,lyr)     [kg]
     real(r8), pointer, public  :: mss_dst3_col(:,:)       ! mass of dust species 3 in snow (col,lyr)     [kg]
     real(r8), pointer, public  :: mss_dst4_col(:,:)       ! mass of dust species 4 in snow (col,lyr)     [kg]
     real(r8), pointer, public  :: mss_dsttot_col(:,:)     ! total mass of dust in snow (col,lyr)         [kg]
     real(r8), pointer, public  :: mss_dst_col_col(:)      ! column-integrated mass of dust in snow       [kg]
     real(r8), pointer, public  :: mss_dst_top_col(:)      ! top-layer mass of dust in snow               [kg]

     real(r8), pointer, public  :: mss_cnc_bcphi_col(:,:)  ! mass concentration of hydrophilic BC in snow (col,lyr) [kg/kg]
     real(r8), pointer, public  :: mss_cnc_bcpho_col(:,:)  ! mass concentration of hydrophilic BC in snow (col,lyr) [kg/kg]
     real(r8), pointer, public  :: mss_cnc_ocphi_col(:,:)  ! mass concentration of hydrophilic OC in snow (col,lyr) [kg/kg]
     real(r8), pointer, public  :: mss_cnc_ocpho_col(:,:)  ! mass concentration of hydrophilic OC in snow (col,lyr) [kg/kg]
     real(r8), pointer, public  :: mss_cnc_dst1_col(:,:)   ! mass concentration of dust species 1 in snow (col,lyr) [kg/kg]
     real(r8), pointer, public  :: mss_cnc_dst2_col(:,:)   ! mass concentration of dust species 2 in snow (col,lyr) [kg/kg]
     real(r8), pointer, public  :: mss_cnc_dst3_col(:,:)   ! mass concentration of dust species 3 in snow (col,lyr) [kg/kg]
     real(r8), pointer, public  :: mss_cnc_dst4_col(:,:)   ! mass concentration of dust species 4 in snow (col,lyr) [kg/kg]

     real(r8), pointer, private :: flx_dst_dep_dry1_col(:) ! dust species 1 dry   deposition on ground (positive definite) [kg/s]
     real(r8), pointer, private :: flx_dst_dep_wet1_col(:) ! dust species 1 wet   deposition on ground (positive definite) [kg/s]
     real(r8), pointer, private :: flx_dst_dep_dry2_col(:) ! dust species 2 dry   deposition on ground (positive definite) [kg/s]
     real(r8), pointer, private :: flx_dst_dep_wet2_col(:) ! dust species 2 wet   deposition on ground (positive definite) [kg/s]
     real(r8), pointer, private :: flx_dst_dep_dry3_col(:) ! dust species 3 dry   deposition on ground (positive definite) [kg/s]
     real(r8), pointer, private :: flx_dst_dep_wet3_col(:) ! dust species 3 wet   deposition on ground (positive definite) [kg/s]
     real(r8), pointer, private :: flx_dst_dep_dry4_col(:) ! dust species 4 dry   deposition on ground (positive definite) [kg/s]
     real(r8), pointer, private :: flx_dst_dep_wet4_col(:) ! dust species 4 wet   deposition on ground (positive definite) [kg/s]
     real(r8), pointer, private :: flx_dst_dep_col(:)      ! total (dry+wet) dust deposition on ground (positive definite) [kg/s]

     real(r8), pointer, private :: flx_bc_dep_dry_col(:)   ! dry (BCPHO+BCPHI) BC deposition on ground (positive definite) [kg/s]
     real(r8), pointer, private :: flx_bc_dep_wet_col(:)   ! wet (BCPHI) BC       deposition on ground (positive definite) [kg/s]
     real(r8), pointer, private :: flx_bc_dep_pho_col(:)   ! hydrophobic BC       deposition on ground (positive definite) [kg/s]
     real(r8), pointer, private :: flx_bc_dep_phi_col(:)   ! hydrophillic BC      deposition on ground (positive definite) [kg/s]
     real(r8), pointer, private :: flx_bc_dep_col(:)       ! total (dry+wet) BC   deposition on ground (positive definite) [kg/s]

     real(r8), pointer, private :: flx_oc_dep_dry_col(:)   ! dry (OCPHO+OCPHI) OC deposition on ground (positive definite) [kg/s]
     real(r8), pointer, private :: flx_oc_dep_wet_col(:)   ! wet (OCPHI) OC       deposition on ground (positive definite) [kg/s]
     real(r8), pointer, private :: flx_oc_dep_pho_col(:)   ! hydrophobic OC       deposition on ground (positive definite) [kg/s]
     real(r8), pointer, private :: flx_oc_dep_phi_col(:)   ! hydrophillic OC      deposition on ground (positive definite) [kg/s]
     real(r8), pointer, private :: flx_oc_dep_col(:)       ! total (dry+wet) OC   deposition on ground (positive definite) [kg/s]

  end type aerosol_type

  character(len=*), parameter, private :: sourcefile = &
       __FILE__
  !-----------------------------------------------------------------------

end module AerosolMod
