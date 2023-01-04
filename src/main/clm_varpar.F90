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

  integer            :: nlevsoi               ! number of hydrologically active soil layers
  integer            :: nlevsoifl             ! number of soil layers on input file
  integer            :: nlevgrnd              ! number of ground layers 
                                              ! (includes lower layers that are hydrologically inactive)
                                              ! (includes lower layers that are biogeochemically inactive)
  integer            :: nlevsno     =  -1     ! maximum number of snow layers
  !ED variables
  integer, parameter :: numrad      =   2     ! number of solar radiation bands: vis, nir
  integer, parameter :: ndst        =   4     ! number of dust size classes (BGC only)

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

    nlevsoifl   =  10
    nlevsoi     =  nlevsoifl
    nlevgrnd    =  15

    if ( masterproc )then
       write(iulog, *) 'CLM varpar subsurface discretization levels '
       write(iulog, '(a, i3)') '    nlevsoi = ', nlevsoi
       write(iulog, '(a, i3)') '    nlevgrnd = ', nlevgrnd
       write(iulog, *)
    end if

  end subroutine clm_varpar_init

end module clm_varpar
