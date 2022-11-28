module clm_varpar

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Module containing CLM parameters
  !
  ! !USES:
  use shr_kind_mod , only: r8 => shr_kind_r8
  use spmdMod      , only: masterproc
  use clm_varctl   , only: iulog

  ! !PUBLIC TYPES:
  implicit none
  save

  ! Note - model resolution is read in from the surface dataset

  integer, parameter :: nlev_equalspace   = 15
  integer, parameter :: toplev_equalspace =  6
  integer            :: nlevsoi               ! number of hydrologically active soil layers
  integer            :: nlevsoifl             ! number of soil layers on input file
  integer            :: nlevgrnd              ! number of ground layers 
                                              ! (includes lower layers that are hydrologically inactive)
  integer            :: nlevurb               ! number of urban layers
  integer            :: nlevlak               ! number of lake layers
                                              ! (includes lower layers that are biogeochemically inactive)
  integer            :: nlevsno     =  -1     ! maximum number of snow layers
  integer, parameter :: nlevcan     =   1     ! number of leaf layers in canopy layer
  !ED variables
  integer, parameter :: numrad      =   2     ! number of solar radiation bands: vis, nir
  integer, parameter :: ndst        =   4     ! number of dust size classes (BGC only)
  integer, parameter :: mxpft       =  78     ! maximum number of PFT's for any mode;
  ! FIX(RF,032414) might we set some of these automatically from reading pft-physiology?
  integer, parameter :: numveg      =  16     ! number of veg types (without specific crop)
  integer, parameter :: nlayer      =   3     ! number of VIC soil layer --Added by AWang
  integer            :: nlayert               ! number of VIC soil layer + 3 lower thermal layers

  integer :: numpft      = mxpft   ! actual # of pfts (without bare)
  integer :: numcft      =  64     ! actual # of crops (includes unused CFTs that are merged into other CFTs)
  integer :: maxpatch_urb= 5       ! max number of urban patches (columns) in urban landunit

  ! Indices used in surface file read and set in clm_varpar_init

  integer :: natpft_lb          ! In PATCH arrays, lower bound of Patches on the natural veg landunit (i.e., bare ground index)
  integer :: natpft_ub          ! In PATCH arrays, upper bound of Patches on the natural veg landunit
  integer :: natpft_size        ! Number of Patches on natural veg landunit (including bare ground)

  ! The following variables pertain to arrays of all PFTs - e.g., those dimensioned (g,
  ! pft_index). These include unused CFTs that are merged into other CFTs. Thus, these
  ! variables do NOT give the actual number of CFTs on the crop landunit - that number
  ! will generally be less because CLM does not simulate all crop types (some crop types
  ! are merged into other types).
  integer :: cft_lb             ! In arrays of PFTs, lower bound of PFTs on the crop landunit
  integer :: cft_ub             ! In arrays of PFTs, upper bound of PFTs on the crop landunit
  integer :: cft_size           ! Number of PFTs on crop landunit in arrays of PFTs

  integer :: max_patch_per_col
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public clm_varpar_init          ! set parameters
  !
  !-----------------------------------------------------------------------

contains

  !------------------------------------------------------------------------------
  subroutine clm_varpar_init()
    !
    ! !DESCRIPTION:
    ! Initialize module variables 
    !
    ! !ARGUMENTS:
    implicit none
    !
    ! !LOCAL VARIABLES:
    !
    character(len=32) :: subname = 'clm_varpar_init'  ! subroutine name
    !------------------------------------------------------------------------------

    ! Crop settings and consistency checks

    numpft      = numveg  ! actual # of patches (without bare)
    numcft      =   2     ! actual # of crops

    ! For arrays containing all Patches (natural veg & crop), determine lower and upper bounds
    ! for (1) Patches on the natural vegetation landunit (includes bare ground, and includes
    ! crops if create_crop_landunit=false), and (2) CFTs on the crop landunit (no elements
    ! if create_crop_landunit=false)

    natpft_size = (numpft + 1) - numcft    ! note that numpft doesn't include bare ground -- thus we add 1
    cft_size    = numcft

    natpft_lb = 0
    natpft_ub = natpft_lb + natpft_size - 1
    cft_lb = natpft_ub + 1
    cft_ub = cft_lb + cft_size - 1

    max_patch_per_col= max(numpft+1, numcft, maxpatch_urb)

    nlevsoifl   =  10
    nlevurb     =  5
    nlevsoi     =  nlevsoifl
    nlevgrnd    =  15

    nlevlak     =  10     ! number of lake layers

    if ( masterproc )then
       write(iulog, *) 'CLM varpar subsurface discretization levels '
       write(iulog, '(a, i3)') '    nlevsoi = ', nlevsoi
       write(iulog, '(a, i3)') '    nlevgrnd = ', nlevgrnd
       write(iulog, '(a, i3)') '    nlevlak = ', nlevlak
       write(iulog, *)
    end if

  end subroutine clm_varpar_init

end module clm_varpar
