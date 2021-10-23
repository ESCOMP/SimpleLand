module FatesInterfaceMod

   ! ------------------------------------------------------------------------------------
   ! This is the FATES public API
   ! A host land model has defined and allocated a structure "fates" as
   ! defined by fates_interface_type
   !
   ! It is also likely/possible that this type is defined as a vector
   ! which is allocated by thread
   ! ------------------------------------------------------------------------------------

   use shr_kind_mod        , only : r8=> shr_kind_r8
   use abortutils          , only : endrun
   ! CIME Globals
   use shr_log_mod         , only : errMsg => shr_log_errMsg
   use shr_infnan_mod      , only : nan => shr_infnan_nan, assignment(=)

   implicit none

   public :: set_fates_global_elements

   character(len=*), parameter, private :: sourcefile = &
         __FILE__
   
   ! -------------------------------------------------------------------------------------
   ! Parameters that are dictated by the Host Land Model
   ! THESE ARE NOT DYNAMIC. SHOULD BE SET ONCE DURING INTIALIZATION.
   ! -------------------------------------------------------------------------------------

  
contains

    ! ===================================================================================
    
    subroutine set_fates_global_elements(use_fates)

       ! --------------------------------------------------------------------------------
       !
       ! This subroutine is called directly from the HLM, and is the first FATES routine
       ! that is called.
       !
       ! This subroutine MUST BE CALLED AFTER the FATES PFT parameter file has been read in,
       ! and the EDPftvarcon_inst structure has been made.
       ! This subroutine must ALSO BE CALLED BEFORE the history file dimensions
       ! are set.
       ! 
       ! This routine requires no information from the HLM. This routine is responsible
       ! for generating the globals that are required by the HLM that are entirely
       ! FATES derived.
       !
       ! --------------------------------------------------------------------------------

      !use CLMFatesParamInterfaceMod         , only : FatesReadParameters
      implicit none
      
      logical,intent(in) :: use_fates    ! Is fates turned on?
      
      integer :: i
      
      if (use_fates) then

            call endrun(msg=errMsg(sourcefile, __LINE__))
         
      else
         ! If we are not using FATES, the cohort dimension is still
         ! going to be initialized, lets set it to the smallest value
         ! possible so that the dimensioning info takes up little space

         !fates_maxElementsPerPatch = 1
      
         !fates_maxElementsPerSite = 1
         

      end if


    end subroutine set_fates_global_elements

    !==============================================================================================
    
end module FatesInterfaceMod
