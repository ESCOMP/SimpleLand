module CNFireEmissionsMod

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Gathers carbon emissions from fire sources to be sent to CAM-Chem via
  ! the coupler .... 
  ! Created by F. Vitt, and revised by F. Li
  ! !USES:
  use shr_kind_mod, only : r8 => shr_kind_r8
  use abortutils,   only : endrun
  use PatchType,    only : patch                
  use decompMod,    only : bounds_type
  use shr_fire_emis_mod,  only : shr_fire_emis_comps_n, shr_fire_emis_comp_t, shr_fire_emis_linkedlist
  use shr_fire_emis_mod,  only : shr_fire_emis_mechcomps_n, shr_fire_emis_mechcomps
  !
  implicit none
  private 
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  !
  ! !PRIVATE TYPES:
  type, private :: emis_t
     real(r8), pointer :: emis(:)
  end type emis_t
  !
  ! !PUBLIC TYPES:
  type, public :: fireemis_type
     real(r8),     pointer, public  :: fireflx_patch(:,:) ! carbon flux from fire sources (kg/m2/sec)
     real(r8),     pointer, public  :: ztop_patch(:)      ! height of the smoke plume (meters)
     type(emis_t), pointer, private :: comp(:)            ! fire emissions component (corresponds to emis factors table input file)
     type(emis_t), pointer, private :: mech(:)            ! cam-chem mechism species emissions
     type(emis_t),          private :: totfire            ! sum of all species emissions
  end type fireemis_type

end module CNFireEmissionsMod

