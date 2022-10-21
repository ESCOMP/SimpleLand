module ndepStreamMod

  !----------------------------------------------------------------------- 
  ! !DESCRIPTION: 
  ! Contains methods for reading in nitrogen deposition data file
  ! Also includes functions for dynamic ndep file handling and 
  ! interpolation.
  !
  ! !USES
  use shr_kind_mod, only: r8 => shr_kind_r8, CL => shr_kind_cl 
  use mct_mod     , only: mct_ggrid
  use spmdMod     , only: mpicom, iam
  use clm_varctl  , only: iulog
  use abortutils  , only: endrun
  use decompMod   , only: bounds_type, gsmap_lnd_gdc2glo 
  use domainMod   , only: ldomain

  ! !PUBLIC TYPES:
  implicit none
  private
  save

  public :: clm_domain_mct ! Sets up MCT domain for this resolution

  ! ! PRIVATE TYPES

  character(len=*), parameter, private :: sourcefile = &
       __FILE__
  !==============================================================================

contains

 !==============================================================================
  subroutine clm_domain_mct(bounds, dom_clm)

    !-------------------------------------------------------------------
    ! Set domain data type for internal clm grid
    use clm_varcon  , only : re
    use domainMod   , only : ldomain
    use seq_flds_mod
    use mct_mod     , only : mct_ggrid, mct_gsMap_lsize, mct_gGrid_init
    use mct_mod     , only : mct_gsMap_orderedPoints, mct_gGrid_importIAttr
    use mct_mod     , only : mct_gGrid_importRAttr
    implicit none
    ! 
    ! arguments
    type(bounds_type), intent(in) :: bounds  
    type(mct_ggrid), intent(out)   :: dom_clm     ! Output domain information for land model
    !
    ! local variables
    integer :: g,i,j              ! index
    integer :: lsize              ! land model domain data size
    real(r8), pointer :: data(:)  ! temporary
    integer , pointer :: idata(:) ! temporary
    !-------------------------------------------------------------------
    !
    ! Initialize mct domain type
    ! lat/lon in degrees,  area in radians^2, mask is 1 (land), 0 (non-land)
    ! Note that in addition land carries around landfrac for the purposes of domain checking
    ! 
    lsize = mct_gsMap_lsize(gsmap_lnd_gdc2glo, mpicom)
    call mct_gGrid_init( GGrid=dom_clm, CoordChars=trim(seq_flds_dom_coord), &
                         OtherChars=trim(seq_flds_dom_other), lsize=lsize )
    !
    ! Allocate memory
    !
    allocate(data(lsize))
    !
    ! Determine global gridpoint number attribute, GlobGridNum, which is set automatically by MCT
    !
    call mct_gsMap_orderedPoints(gsmap_lnd_gdc2glo, iam, idata)
    call mct_gGrid_importIAttr(dom_clm,'GlobGridNum',idata,lsize)
    !
    ! Determine domain (numbering scheme is: West to East and South to North to South pole)
    ! Initialize attribute vector with special value
    !
    data(:) = -9999.0_R8 
    call mct_gGrid_importRAttr(dom_clm,"lat"  ,data,lsize) 
    call mct_gGrid_importRAttr(dom_clm,"lon"  ,data,lsize) 
    call mct_gGrid_importRAttr(dom_clm,"area" ,data,lsize) 
    call mct_gGrid_importRAttr(dom_clm,"aream",data,lsize) 
    data(:) = 0.0_R8     
    call mct_gGrid_importRAttr(dom_clm,"mask" ,data,lsize) 
    !
    ! Determine bounds
    !
    ! Fill in correct values for domain components
    ! Note aream will be filled in in the atm-lnd mapper
    !
    do g = bounds%begg,bounds%endg
       i = 1 + (g - bounds%begg)
       data(i) = ldomain%lonc(g)
    end do
    call mct_gGrid_importRattr(dom_clm,"lon",data,lsize) 

    do g = bounds%begg,bounds%endg
       i = 1 + (g - bounds%begg)
       data(i) = ldomain%latc(g)
    end do
    call mct_gGrid_importRattr(dom_clm,"lat",data,lsize) 

    do g = bounds%begg,bounds%endg
       i = 1 + (g - bounds%begg)
       data(i) = ldomain%area(g)/(re*re)
    end do
    call mct_gGrid_importRattr(dom_clm,"area",data,lsize) 

    do g = bounds%begg,bounds%endg
       i = 1 + (g - bounds%begg)
       data(i) = real(ldomain%mask(g), r8)
    end do
    call mct_gGrid_importRattr(dom_clm,"mask",data,lsize) 

    do g = bounds%begg,bounds%endg
       i = 1 + (g - bounds%begg)
       data(i) = real(ldomain%frac(g), r8)
    end do
    call mct_gGrid_importRattr(dom_clm,"frac",data,lsize) 

    deallocate(data)
    deallocate(idata)

  end subroutine clm_domain_mct
    
end module ndepStreamMod

