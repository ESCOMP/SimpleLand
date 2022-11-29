module ActiveLayerMod
  
  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Module holding routines for calculation of active layer dynamics
  !
  ! !USES:
  use shr_kind_mod    , only : r8 => shr_kind_r8
  use CanopyStateType , only : canopystate_type
  !
  implicit none
  save
  private
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public:: alt_calc
  !-----------------------------------------------------------------------
  
contains

  !-----------------------------------------------------------------------
  subroutine alt_calc(num_soilc, filter_soilc, canopystate_inst)
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    integer                , intent(in)    :: num_soilc       ! number of soil columns in filter
    integer                , intent(in)    :: filter_soilc(:) ! filter for soil columns
    type(canopystate_type) , intent(inout) :: canopystate_inst
    !
    ! !LOCAL VARIABLES:
    integer :: c, fc  ! counters
    !-----------------------------------------------------------------------

    associate(                                                                & 
         alt                  =>    canopystate_inst%alt_col             ,    & ! Output:  [real(r8) (:)   ]  current depth of thaw                                                 
         altmax               =>    canopystate_inst%altmax_col          ,    & ! Output:  [real(r8) (:)   ]  maximum annual depth of thaw                                          
         altmax_lastyear      =>    canopystate_inst%altmax_lastyear_col ,    & ! Output:  [real(r8) (:)   ]  prior year maximum annual depth of thaw                               
         alt_indx             =>    canopystate_inst%alt_indx_col        ,    & ! Output:  [integer  (:)   ]  current depth of thaw                                                  
         altmax_indx          =>    canopystate_inst%altmax_indx_col     ,    & ! Output:  [integer  (:)   ]  maximum annual depth of thaw                                           
         altmax_lastyear_indx =>    canopystate_inst%altmax_lastyear_indx_col & ! Output:  [integer  (:)   ]  prior year maximum annual depth of thaw                                
         )

      do fc = 1,num_soilc
         c = filter_soilc(fc)
               alt(c)=0._r8
               alt_indx(c) = 0
               altmax(c) = alt(c)
               altmax_indx(c) = alt_indx(c)
               altmax_lastyear(c) = altmax(c)
               altmax_lastyear_indx(c) = altmax_indx(c)
      end do

    end associate 

  end subroutine alt_calc
  
end module ActiveLayerMod
