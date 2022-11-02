Module DryDepVelocity                                              

  !-----------------------------------------------------------------------  
  !  
  ! Purpose:  
  ! Deposition velocity (m/s) 
  !  
  ! Method:  
  ! This code simulates dry deposition velocities using the Wesely scheme.   
  ! Details of this method can be found in:  
  ! 
  ! M.L Wesely. Parameterization of surface resistances to gaseous dry deposition 
  ! in regional-scale numericl models. 1989. Atmospheric Environment vol.23 No.6  
  ! pp. 1293-1304. 
  ! 
  ! In Wesely (1998) "the magnitude of the dry deposition velocity can be found 
  ! as: 
  ! 
  !  |vd|=(ra+rb+rc)^-1 
  ! 
  ! where ra is the aerodynamic resistance (common to all gases) between a  
  ! specific height and the surface, rb is the quasilaminar sublayer resistance 
  ! (whose only dependence on the porperties of the gas of interest is its 
  ! molecular diffusivity in air), and rc is the bulk surface resistance". 
  ! 
  ! In this subroutine both ra and rb are calculated elsewhere in CLM.  
  ! 
  ! In Wesely (1989) rc is estimated for five seasonal categories and 11 landuse 
  ! types.  For each season and landuse type, Wesely compiled data into a  
  ! look-up-table for several parameters used to calculate rc. In this subroutine 
  ! the same values are used as found in wesely's look-up-tables, the only  
  ! difference is that this subroutine uses a CLM generated LAI to select values 
  ! from the look-up-table instead of seasonality.  Inaddition, Wesely(1989)  
  ! 
  ! Subroutine written to operate at the patch level. 
  ! 
  ! Output: 
  ! 
  ! vd(n_species) !Dry deposition velocity [m s-1] for each molecule or species 
  !  
  ! Author: Beth Holland and  James Sulzman 
  ! 
  ! Modified: Francis Vitt -- 30 Mar 2007
  ! Modified: Maria Val Martin -- 15 Jan 2014
  !  Corrected major bugs in the leaf and stomatal resitances. The code is now
  !  coupled to LAI and Rs uses the Ball-Berry Scheme. Also, corrected minor 
  !  bugs in rlu and rcl calculations. Added 
  !  no vegetation removal for CO. See README for details and 
  !     Val Martin et al., 2014 GRL for major corrections
  ! Modified: Louisa Emmons -- 30 November 2017
  !  Corrected the equation calculating stomatal resistance from rssun and rssha,
  !     and removed factor that scaled Rs to match observations
  !
  !----------------------------------------------------------------------- 

  use shr_kind_mod         , only : r8 => shr_kind_r8 
  use decompMod            , only : bounds_type
  !
  implicit none 
  private
  !
  type, public :: drydepvel_type

     real(r8), pointer, public  :: velocity_patch (:,:) ! Dry Deposition Velocity
   contains

     procedure , public  :: Init 
     procedure , private :: InitAllocate 

  end type drydepvel_type
  !----------------------------------------------------------------------- 

  character(len=*), parameter, private :: sourcefile = &
       __FILE__

CONTAINS 

  !------------------------------------------------------------------------
  subroutine Init(this, bounds)

    class(drydepvel_type) :: this
    type(bounds_type), intent(in) :: bounds  

    call this%InitAllocate(bounds)

  end subroutine Init

  !------------------------------------------------------------------------
  subroutine InitAllocate(this, bounds)
    !
    ! !USES:
    use shr_infnan_mod , only : nan => shr_infnan_nan, assignment(=)
    use seq_drydep_mod , only : n_drydep, drydep_method, DD_XLND
    !
    ! !ARGUMENTS:
    class(drydepvel_type) :: this
    type(bounds_type), intent(in) :: bounds  
    !
    ! !LOCAL VARIABLES:
    integer :: begp, endp
    !------------------------------------------------------------------------

    begp = bounds%begp; endp= bounds%endp

    ! Dry Deposition Velocity 
    if ( n_drydep > 0 .and. drydep_method == DD_XLND )then
       allocate(this%velocity_patch(begp:endp, n_drydep));  this%velocity_patch(:,:) = nan 
    end if

  end subroutine InitAllocate

end module DryDepVelocity
