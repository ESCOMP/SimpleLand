module SoilHydrologyInitTimeConstMod

  !------------------------------------------------------------------------------
  ! DESCRIPTION:
  ! Initialize time constant variables for SoilHydrologyType
  !
  ! !USES
  use shr_kind_mod      , only : r8 => shr_kind_r8
  use shr_log_mod       , only : errMsg => shr_log_errMsg
  use decompMod         , only : bounds_type
  use SoilHydrologyType , only : soilhydrology_type
  use LandunitType      , only : lun                
  use ColumnType        , only : col                
  !
  implicit none
  private
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: SoilHydrologyInitTimeConst
  !
  ! !PRIVATE MEMBER FUNCTIONS:
  private :: initSoilParVIC    ! Convert default CLM soil properties to VIC parameters
  private :: initCLMVICMap     ! Initialize map from VIC to CLM layers
  private :: linear_interp     ! function for linear interperation 

  character(len=*), parameter, private :: sourcefile = &
       __FILE__
  !-----------------------------------------------------------------------
  !
contains

  !-----------------------------------------------------------------------
  subroutine SoilHydrologyInitTimeConst(bounds, soilhydrology_inst) 
    !
    ! !USES:
    use shr_const_mod   , only : shr_const_pi
    use shr_spfn_mod    , only : shr_spfn_erf
    use clm_varpar      , only : nlevsoi
    use clm_varcon      , only : spval, pc
    use clm_varcon      , only : aquifer_water_baseline
    use column_varcon   , only : icol_road_perv
    !
    ! !ARGUMENTS:
    type(bounds_type)        , intent(in)    :: bounds                                    
    type(soilhydrology_type) , intent(inout) :: soilhydrology_inst
    !
    ! !LOCAL VARIABLES:
    integer            :: p,c,j,l,g,lev,nlevs 
    real(r8)           :: d, fd, dfdd
    real(r8)           :: clay,sand        ! temporaries
    real(r8), pointer  :: sandcol    (:,:) ! column level sand fraction for calculating VIC parameters
    real(r8), pointer  :: claycol    (:,:) ! column level clay fraction for calculating VIC parameters
    !-----------------------------------------------------------------------
    ! -----------------------------------------------------------------
    ! Initialize frost table
    ! -----------------------------------------------------------------

    soilhydrology_inst%wa_col(bounds%begc:bounds%endc)  = aquifer_water_baseline
    soilhydrology_inst%zwt_col(bounds%begc:bounds%endc) = 0._r8

    do c = bounds%begc,bounds%endc
       l = col%landunit(c)
       if (.not. lun%lakpoi(l)) then  !not lake
          if (lun%urbpoi(l)) then
             if (col%itype(c) == icol_road_perv) then
                ! Note that the following hard-coded constants (on the next two lines)
                ! seem implicitly related to aquifer_water_baseline
                soilhydrology_inst%wa_col(c)  = 4800._r8
                soilhydrology_inst%zwt_col(c) = (25._r8 + col%zi(c,nlevsoi)) - soilhydrology_inst%wa_col(c)/0.2_r8 /1000._r8  ! One meter below soil column
             else
                soilhydrology_inst%wa_col(c)  = spval
                soilhydrology_inst%zwt_col(c) = spval
             end if
             ! initialize frost_table, zwt_perched
             soilhydrology_inst%zwt_perched_col(c) = spval
             soilhydrology_inst%frost_table_col(c) = spval
          else
             ! Note that the following hard-coded constants (on the next two lines) seem
             ! implicitly related to aquifer_water_baseline
             soilhydrology_inst%wa_col(c)  = 4000._r8
             soilhydrology_inst%zwt_col(c) = (25._r8 + col%zi(c,nlevsoi)) - soilhydrology_inst%wa_col(c)/0.2_r8 /1000._r8  ! One meter below soil column
             ! initialize frost_table, zwt_perched to bottom of soil column
             soilhydrology_inst%zwt_perched_col(c) = col%zi(c,nlevsoi)
             soilhydrology_inst%frost_table_col(c) = col%zi(c,nlevsoi)
          end if
       end if
    end do

    associate(micro_sigma => col%micro_sigma)
      do c = bounds%begc, bounds%endc
         
         ! determine h2osfc threshold ("fill & spill" concept)
         ! set to zero for no h2osfc (w/frac_infclust =large)
         
         soilhydrology_inst%h2osfc_thresh_col(c) = 0._r8
         if (micro_sigma(c) > 1.e-6_r8 .and. (soilhydrology_inst%h2osfcflag /= 0)) then
            d = 0.0
            do p = 1,4
               fd   = 0.5*(1.0_r8+shr_spfn_erf(d/(micro_sigma(c)*sqrt(2.0)))) - pc
               dfdd = exp(-d**2/(2.0*micro_sigma(c)**2))/(micro_sigma(c)*sqrt(2.0*shr_const_pi))
               d    = d - fd/dfdd
            enddo
            soilhydrology_inst%h2osfc_thresh_col(c) = 0.5*d*(1.0_r8+shr_spfn_erf(d/(micro_sigma(c)*sqrt(2.0)))) + &
                 micro_sigma(c)/sqrt(2.0*shr_const_pi)*exp(-d**2/(2.0*micro_sigma(c)**2))         
            soilhydrology_inst%h2osfc_thresh_col(c) = 1.e3_r8 * soilhydrology_inst%h2osfc_thresh_col(c) !convert to mm from meters
         else
            soilhydrology_inst%h2osfc_thresh_col(c) = 0._r8
         endif

         if (soilhydrology_inst%h2osfcflag == 0) then 
            soilhydrology_inst%h2osfc_thresh_col(c) = 0._r8    ! set to zero for no h2osfc (w/frac_infclust =large)
         endif

         ! set decay factor
         soilhydrology_inst%hkdepth_col(c) = 1._r8/2.5_r8

      end do
    end associate

  end subroutine SoilhydrologyInitTimeConst

  !-----------------------------------------------------------------------
  subroutine initSoilParVIC(c, claycol, sandcol, om_fraccol, soilhydrology_inst)
    !
    ! !DESCRIPTION:
    ! Convert default CLM soil properties to VIC parameters
    ! to be used for runoff simulations (added by M. Huang)
    !
    ! !USES:
    use clm_varpar, only : nlevsoi, nlayert, nlayer
    !
    ! !ARGUMENTS:
    integer                  , intent(in)    :: c               ! column index
    real(r8)                 , pointer       :: sandcol(:,:)    ! read in - soil texture: percent sand
    real(r8)                 , pointer       :: claycol(:,:)    ! read in - soil texture: percent clay
    real(r8)                 , pointer       :: om_fraccol(:,:) ! read in - organic matter: kg/m3
    type(soilhydrology_type) , intent(inout) :: soilhydrology_inst

    ! !LOCAL VARIABLES:
    real(r8) :: om_watsat    = 0.9_r8             ! porosity of organic soil
    real(r8) :: om_hksat     = 0.1_r8             ! saturated hydraulic conductivity of organic soil [mm/s]
    real(r8) :: om_expt      = 3._r8+2._r8*2.7_r8 ! soil expt for VIC        
    real(r8) :: pc           = 0.5_r8             ! percolation threshold
    real(r8) :: pcbeta       = 0.139_r8           ! percolation exponent
    real(r8) :: xksat                             ! maximum hydraulic conductivity of soil [mm/s]
    real(r8) :: perc_frac                         ! "percolating" fraction of organic soil
    real(r8) :: perc_norm                         ! normalize to 1 when 100% organic soil
    real(r8) :: uncon_hksat                       ! series conductivity of mineral/organic soil
    real(r8) :: uncon_frac                        ! fraction of "unconnected" soil
    real(r8) :: temp_sum_frac                     ! sum of node fractions in each VIC layer
    real(r8) :: sandvic(1:nlayert)                ! temporary, weighted averaged sand% for VIC layers
    real(r8) :: clayvic(1:nlayert)                ! temporary, weighted averaged clay% for VIC layers
    real(r8) :: om_fracvic(1:nlayert)             ! temporary, weighted averaged organic matter fract for VIC layers
    integer  :: i, j                              ! indices
    !-------------------------------------------------------------------------------------------

    ! soilhydrology_inst%depth_col(:,:)           Output: layer depth of upper layer(m) 
    ! soilhydrology_inst%vic_clm_fract_col(:,:,:) Output: fraction of VIC layers in CLM layers
    ! soilhydrology_inst%c_param_col(:)           Output: baseflow exponent (Qb)
    ! soilhydrology_inst%expt_col(:,:)            Output: pore-size distribution related paramter(Q12)
    ! soilhydrology_inst%ksat_col(:,:)            Output: Saturated hydrologic conductivity (mm/s)
    ! soilhydrology_inst%phi_s_col(:,:)           Output: soil moisture dissusion parameter
    ! soilhydrology_inst%porosity_col(:,:)        Output: soil porosity
    ! soilhydrology_inst%max_moist_col(:,:)       Output: maximum soil moisture (ice + liq)

    !  map parameters between VIC layers and CLM layers
    soilhydrology_inst%c_param_col(c) = 2._r8

    ! map the CLM layers to VIC layers 
    do i = 1, nlayer      

       sandvic(i)    = 0._r8
       clayvic(i)    = 0._r8   
       om_fracvic(i) = 0._r8  
       temp_sum_frac = 0._r8     
       do j = 1, nlevsoi
          sandvic(i)    = sandvic(i)    + sandcol(c,j)    * soilhydrology_inst%vic_clm_fract_col(c,i,j)
          clayvic(i)    = clayvic(i)    + claycol(c,j)    * soilhydrology_inst%vic_clm_fract_col(c,i,j)
          om_fracvic(i) = om_fracvic(i) + om_fraccol(c,j) * soilhydrology_inst%vic_clm_fract_col(c,i,j) 
          temp_sum_frac = temp_sum_frac +                   soilhydrology_inst%vic_clm_fract_col(c,i,j)
       end do

       !average soil properties, M.Huang, 08/11/2010
       sandvic(i) = sandvic(i)/temp_sum_frac
       clayvic(i) = clayvic(i)/temp_sum_frac
       om_fracvic(i) = om_fracvic(i)/temp_sum_frac

       !make sure sand, clay and om fractions are between 0 and 100% 
       sandvic(i)    = min(100._r8 , sandvic(i))
       clayvic(i)    = min(100._r8 , clayvic(i))
       om_fracvic(i) = min(100._r8 , om_fracvic(i))
       sandvic(i)    = max(0._r8   , sandvic(i))
       clayvic(i)    = max(0._r8   , clayvic(i))
       om_fracvic(i) = max(0._r8   , om_fracvic(i))

       !calculate other parameters based on teh percentages
       soilhydrology_inst%porosity_col(c, i) = 0.489_r8 - 0.00126_r8*sandvic(i)
       soilhydrology_inst%expt_col(c, i)     = 3._r8+ 2._r8*(2.91_r8 + 0.159_r8*clayvic(i))
       xksat = 0.0070556 *( 10.**(-0.884+0.0153*sandvic(i)) )

       !consider organic matter, M.Huang 
       soilhydrology_inst%expt_col(c, i)    = &
            (1._r8 - om_fracvic(i))*soilhydrology_inst%expt_col(c, i)    + om_fracvic(i)*om_expt 
       soilhydrology_inst%porosity_col(c,i) = &
            (1._r8 - om_fracvic(i))*soilhydrology_inst%porosity_col(c,i) + om_watsat*om_fracvic(i) 

       ! perc_frac is zero unless perf_frac greater than percolation threshold
       if (om_fracvic(i) > pc) then
          perc_norm=(1._r8 - pc)**(-pcbeta)
          perc_frac=perc_norm*(om_fracvic(i) - pc)**pcbeta
       else
          perc_frac=0._r8
       endif
       ! uncon_frac is fraction of mineral soil plus fraction of "nonpercolating" organic soil
       uncon_frac=(1._r8-om_fracvic(i))+(1._r8-perc_frac)*om_fracvic(i)

       ! uncon_hksat is series addition of mineral/organic conductivites
       if (om_fracvic(i) < 1._r8) then
          uncon_hksat=uncon_frac/((1._r8-om_fracvic(i))/xksat &
               +((1._r8-perc_frac)*om_fracvic(i))/om_hksat)
       else
          uncon_hksat = 0._r8
       end if

       soilhydrology_inst%ksat_col(c,i)  = &
            uncon_frac*uncon_hksat + (perc_frac*om_fracvic(i))*om_hksat

       soilhydrology_inst%max_moist_col(c,i) = &
            soilhydrology_inst%porosity_col(c,i) * soilhydrology_inst%depth_col(c,i) * 1000._r8 !in mm!

       soilhydrology_inst%phi_s_col(c,i) = &
            -(exp((1.54_r8 - 0.0095_r8*sandvic(i) + &
            0.0063_r8*(100.0_r8-sandvic(i)-clayvic(i)))*log(10.0_r8))*9.8e-5_r8)

    end do ! end of loop over layers

  end subroutine initSoilParVIC

   !-----------------------------------------------------------------------
   subroutine initCLMVICMap(c, soilhydrology_inst)
     !
     ! !DESCRIPTION:
     ! Calculates mapping between CLM and VIC layers
     ! added by AWang, modified by M.Huang for CLM4 
     ! NOTE: in CLM h2osoil_liq unit is kg/m2, in VIC moist is mm
     ! h2osoi_ice is actually water equavlent ice content.
     !
     ! !USES:
     use clm_varpar  , only : nlevsoi, nlayer
     !
     ! !ARGUMENTS:
     integer , intent(in)  :: c
     type(soilhydrology_type), intent(inout) :: soilhydrology_inst
     !
     ! !REVISION HISTORY:
     ! Created by Maoyi Huang
     ! 11/13/2012, Maoyi Huang: rewrite the mapping modules in CLM4VIC 
     !
     ! !LOCAL VARIABLES
     real(r8) :: sum_frac(1:nlayer)                  ! sum of fraction for each layer
     real(r8) :: deltal(1:nlayer+1)                  ! temporary
     real(r8) :: zsum                                ! temporary
     real(r8) :: lsum                                ! temporary
     real(r8) :: temp                                ! temporary
     integer :: i, j, fc
     !-----------------------------------------------------------------------

     associate(                                                    & 
          dz            =>    col%dz    ,                          & ! Input:  [real(r8) (:,:)   ]  layer depth (m)                       
          zi            =>    col%zi    ,                          & ! Input:  [real(r8) (:,:)   ]  interface level below a "z" level (m) 
          z             =>    col%z     ,                          & ! Input:  [real(r8) (:,:)   ]  layer thickness (m)                   

          depth         =>    soilhydrology_inst%depth_col ,       & ! Input:  [real(r8) (:,:)   ]  layer depth of VIC (m)                
          vic_clm_fract =>    soilhydrology_inst%vic_clm_fract_col & ! Output: [real(r8) (:,:,:) ]  fraction of VIC layers in clm layers
          )

       !  set fraction of VIC layer in each CLM layer

       lsum = 0._r8
       do i = 1, nlayer
          deltal(i) = depth(c,i)
       end do
       do i = 1, nlayer
          zsum = 0._r8
          sum_frac(i) = 0._r8
          do j = 1, nlevsoi
             if( (zsum < lsum) .and. (zsum + dz(c,j) >= lsum ))  then
                call linear_interp(lsum, temp, zsum, zsum + dz(c,j), 0._r8, 1._r8)
                vic_clm_fract(c,i,j) = 1._r8 - temp
                if(lsum + deltal(i) < zsum + dz(c,j)) then
                   call linear_interp(lsum + deltal(i), temp, zsum, zsum + dz(c,j), 1._r8, 0._r8)
                   vic_clm_fract(c,i,j) = vic_clm_fract(c,i,j) - temp
                end if
             else if( (zsum < lsum + deltal(i)) .and. (zsum + dz(c,j) >= lsum + deltal(i)) ) then
                call linear_interp(lsum + deltal(i), temp, zsum, zsum + dz(c,j), 0._r8, 1._r8)
                vic_clm_fract(c,i,j) = temp
                if(zsum<=lsum) then
                   call linear_interp(lsum, temp, zsum, zsum + dz(c,j), 0._r8, 1._r8)
                   vic_clm_fract(c,i,j) = vic_clm_fract(c,i,j) - temp
                end if
             else if( (zsum >= lsum .and. zsum + dz(c,j) <= lsum + deltal(i)) )  then
                vic_clm_fract(c,i,j) = 1._r8
             else
                vic_clm_fract(c,i,j) = 0._r8
             end if
             zsum = zsum + dz(c,j)
             sum_frac(i) = sum_frac(i) + vic_clm_fract(c,i,j)
          end do                           ! end CLM layer calculation
          lsum = lsum + deltal(i)
       end do                             ! end VIC layer calcultion 

     end associate 

   end subroutine initCLMVICMap

   !-------------------------------------------------------------------
   subroutine linear_interp(x,y, x0, x1, y0, y1)
     !
     ! !DESCRIPTION:
     ! Provides linear interpolation
     !
     ! !ARGUMENTS:
     real(r8), intent(in)  :: x, x0, y0, x1, y1
     real(r8), intent(out) :: y
     !-------------------------------------------------------------------

     y = y0 + (x - x0) * (y1 - y0) / (x1 - x0)

   end subroutine linear_interp

end module SoilHydrologyInitTimeConstMod
