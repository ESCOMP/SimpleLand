module SoilBiogeochemVerticalProfileMod

#include "shr_assert.h"

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  !  Calculate vertical profiles for distributing soil and litter C and N
  !
  ! !USES:
  use shr_kind_mod, only: r8 => shr_kind_r8
  !
  implicit none
  private
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public:: SoilBiogeochemVerticalProfile
  !
  real(r8), public :: surfprof_exp  = 10. ! how steep profile is for surface components (1/ e_folding depth) (1/m)      

  character(len=*), parameter, private :: sourcefile = &
       __FILE__
  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine SoilBiogeochemVerticalProfile(bounds, num_soilc,filter_soilc,num_soilp,filter_soilp, &
       canopystate_inst, soilstate_inst, soilbiogeochem_state_inst)
    !
    ! !DESCRIPTION:
    !  calculate vertical profiles for distributing soil and litter C and N
    !
    !  BUG(wjs, 2014-12-15, bugz 2107) 
    !  Because of this routine's placement in the driver sequence (it is
    !  called very early in each timestep, before weights are adjusted and filters are
    !  updated), it may be necessary for this routine to compute values over inactive as well
    !  as active points (since some inactive points may soon become active) - so that's what
    !  is done now. Currently, it seems to be okay to do this, because the variables computed
    !  here seem to only depend on quantities that are valid over inactive as well as active
    !  points. However, note that this routine is (mistakenly) called from two places
    !  currently - the above note applies to its call from the driver, but its call from
    !  CNDecompMod uses the standard filters that just apply over active points
    ! 
    ! !USES:
    use shr_log_mod             , only : errMsg => shr_log_errMsg
    use decompMod               , only : bounds_type
    use abortutils              , only : endrun
    use clm_varcon              , only : zsoi, dzsoi, zisoi, dzsoi_decomp, zmin_bedrock
    use clm_varpar              , only : nlevdecomp, nlevgrnd, nlevdecomp_full, maxpatch_pft
    use clm_varctl              , only : iulog
    use pftconMod               , only : noveg, pftcon
    use SoilBiogeochemStateType , only : soilbiogeochem_state_type
    use CanopyStateType         , only : canopystate_type
    use SoilStateType           , only : soilstate_type
    use ColumnType              , only : col                
    use PatchType               , only : patch                
    !
    ! !ARGUMENTS:
    type(bounds_type)               , intent(in)    :: bounds  
    integer                         , intent(in)    :: num_soilc                               ! number of soil columns in filter
    integer                         , intent(in)    :: filter_soilc(:)                         ! filter for soil columns
    integer                         , intent(in)    :: num_soilp                               ! number of soil patches in filter
    integer                         , intent(in)    :: filter_soilp(:)                         ! filter for soil patches
    type(canopystate_type)          , intent(in)    :: canopystate_inst
    type(soilstate_type)            , intent(in)    :: soilstate_inst				    
    type(soilbiogeochem_state_type) , intent(inout) :: soilbiogeochem_state_inst
    !
    ! !LOCAL VARIABLES:
    real(r8) :: surface_prof(1:nlevdecomp)
    real(r8) :: surface_prof_tot
    real(r8) :: rootfr_tot
    real(r8) :: cinput_rootfr(bounds%begp:bounds%endp, 1:nlevdecomp_full)      ! pft-native root fraction used for calculating inputs
    real(r8) :: col_cinput_rootfr(bounds%begc:bounds%endc, 1:nlevdecomp_full)  ! col-native root fraction used for calculating inputs
    integer  :: c, j, fc, p, fp, pi
    integer  :: alt_ind
    ! debugging temp variables
    real(r8) :: froot_prof_sum
    real(r8) :: croot_prof_sum
    real(r8) :: leaf_prof_sum
    real(r8) :: stem_prof_sum
    real(r8) :: ndep_prof_sum
    real(r8) :: nfixation_prof_sum
    real(r8) :: delta = 1.e-10
    integer  :: begp, endp
    integer  :: begc, endc
    character(len=32) :: subname = 'SoilBiogeochemVerticalProfile'
    !-----------------------------------------------------------------------

    begp = bounds%begp; endp= bounds%endp
    begc = bounds%begc; endc= bounds%endc

    associate(                                                                  & 
         altmax_lastyear_indx => canopystate_inst%altmax_lastyear_indx_col    , & ! Input:  [integer   (:)   ]  frost table depth (m)                              

         crootfr              => soilstate_inst%crootfr_patch                 , & ! Input:  [real(r8)  (:,:) ]  fraction of roots in each soil layer  (nlevgrnd)
         
         nfixation_prof       => soilbiogeochem_state_inst%nfixation_prof_col , & ! Input  :  [real(r8) (:,:) ]  (1/m) profile for N fixation additions          
         ndep_prof            => soilbiogeochem_state_inst%ndep_prof_col      , & ! Input  :  [real(r8) (:,:) ]  (1/m) profile for N fixation additions          
         leaf_prof            => soilbiogeochem_state_inst%leaf_prof_patch    , & ! Output :  [real(r8) (:,:) ]  (1/m) profile of leaves                         
         froot_prof           => soilbiogeochem_state_inst%froot_prof_patch   , & ! Output :  [real(r8) (:,:) ]  (1/m) profile of fine roots                     
         croot_prof           => soilbiogeochem_state_inst%croot_prof_patch   , & ! Output :  [real(r8) (:,:) ]  (1/m) profile of coarse roots                   
         stem_prof            => soilbiogeochem_state_inst%stem_prof_patch      & ! Output :  [real(r8) (:,:) ]  (1/m) profile of stems                          
         )

      ! for one layer decomposition model, set profiles to unity
      leaf_prof(begp:endp, :) = 1._r8
      froot_prof(begp:endp, :) = 1._r8
      croot_prof(begp:endp, :) = 1._r8
      stem_prof(begp:endp, :) = 1._r8
      nfixation_prof(begc:endc, :) = 1._r8
      ndep_prof(begc:endc, :) = 1._r8

      ! check to make sure integral of all profiles = 1.
      do fc = 1,num_soilc
         c = filter_soilc(fc)
         ndep_prof_sum = 0.
         nfixation_prof_sum = 0.
         do j = 1, nlevdecomp
            ndep_prof_sum = ndep_prof_sum + ndep_prof(c,j) *  dzsoi_decomp(j)
            nfixation_prof_sum = nfixation_prof_sum + nfixation_prof(c,j) *  dzsoi_decomp(j)
         end do
         if ( ( abs(ndep_prof_sum - 1._r8) > delta ) .or.  ( abs(nfixation_prof_sum - 1._r8) > delta ) ) then
            write(iulog, *) 'profile sums: ', ndep_prof_sum, nfixation_prof_sum
            write(iulog, *) 'c: ', c
            write(iulog, *) 'altmax_lastyear_indx: ', altmax_lastyear_indx(c)
            write(iulog, *) 'nfixation_prof: ', nfixation_prof(c,:)
            write(iulog, *) 'ndep_prof: ', ndep_prof(c,:)
            write(iulog, *) 'cinput_rootfr: ', cinput_rootfr(c,:)
            write(iulog, *) 'dzsoi_decomp: ', dzsoi_decomp(:)
            write(iulog, *) 'surface_prof: ', surface_prof(:)
            write(iulog, *) 'npfts(c): ', col%npatches(c)
            do p = col%patchi(c), col%patchi(c) + col%npatches(c) -1
               write(iulog, *) 'p, itype(p), wtcol(p): ', p, patch%itype(p), patch%wtcol(p)
               write(iulog, *) 'cinput_rootfr(p,:): ', cinput_rootfr(p,:)
            end do
            call endrun(msg=" ERROR: _prof_sum-1>delta"//errMsg(sourcefile, __LINE__))
         endif
      end do

      do fp = 1,num_soilp
         p = filter_soilp(fp)
         froot_prof_sum = 0.
         croot_prof_sum = 0.
         leaf_prof_sum = 0.
         stem_prof_sum = 0.
         do j = 1, nlevdecomp
            froot_prof_sum = froot_prof_sum + froot_prof(p,j) *  dzsoi_decomp(j)
            croot_prof_sum = croot_prof_sum + croot_prof(p,j) *  dzsoi_decomp(j)
            leaf_prof_sum = leaf_prof_sum + leaf_prof(p,j) *  dzsoi_decomp(j)
            stem_prof_sum = stem_prof_sum + stem_prof(p,j) *  dzsoi_decomp(j)
         end do
         if ( ( abs(froot_prof_sum - 1._r8) > delta ) .or.  ( abs(croot_prof_sum - 1._r8) > delta ) .or. &
              ( abs(stem_prof_sum - 1._r8) > delta ) .or.  ( abs(leaf_prof_sum - 1._r8) > delta ) ) then
            write(iulog, *) 'profile sums: ', froot_prof_sum, croot_prof_sum, leaf_prof_sum, stem_prof_sum
            call endrun(msg=' ERROR: sum-1 > delta'//errMsg(sourcefile, __LINE__))
         endif
      end do

    end associate

  end subroutine SoilBiogeochemVerticalProfile
  
end module SoilBiogeochemVerticalProfileMod
