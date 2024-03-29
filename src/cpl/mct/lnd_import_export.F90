module lnd_import_export

  use shr_kind_mod , only: r8 => shr_kind_r8, cl=>shr_kind_cl
  use abortutils   , only: endrun
  use decompmod    , only: bounds_type
  use lnd2atmType  , only: lnd2atm_type
  use atm2lndType  , only: atm2lnd_type
  use clm_cpl_indices
  !
  implicit none
  !===============================================================================

contains

  !===============================================================================
  subroutine lnd_import( bounds, x2l, atm2lnd_inst)

    !---------------------------------------------------------------------------
    ! !DESCRIPTION:
    ! Convert the input data from the coupler to the land model 
    !
    ! !USES:
    use seq_flds_mod    , only: seq_flds_x2l_fields
    use clm_varctl      , only: iulog
    use clm_varcon      , only: rair, o2_molar_const
    use shr_const_mod   , only: SHR_CONST_TKFRZ
    use shr_string_mod  , only: shr_string_listGetName
    use domainMod       , only: ldomain
    use shr_infnan_mod  , only : isnan => shr_infnan_isnan
    !
    ! !ARGUMENTS:
    type(bounds_type)  , intent(in)    :: bounds   ! bounds
    real(r8)           , intent(in)    :: x2l(:,:) ! driver import state to land model
    type(atm2lnd_type) , intent(inout) :: atm2lnd_inst      ! clm internal input data type
    !
    ! !LOCAL VARIABLES:
    integer  :: g,i,k,nstep,ier      ! indices, number of steps, and error code
    real(r8) :: forc_rainc           ! rainxy Atm flux mm/s
    real(r8) :: e                    ! vapor pressure (Pa)
    real(r8) :: qsat                 ! saturation specific humidity (kg/kg)
    real(r8) :: forc_t               ! atmospheric temperature (Kelvin)
    real(r8) :: forc_q               ! atmospheric specific humidity (kg/kg)
    real(r8) :: forc_pbot            ! atmospheric pressure (Pa)
    real(r8) :: forc_rainl           ! rainxy Atm flux mm/s
    real(r8) :: forc_snowc           ! snowfxy Atm flux  mm/s
    real(r8) :: forc_snowl           ! snowfxl Atm flux  mm/s
    real(r8) :: esatw                ! saturation vapor pressure over water (Pa)
    real(r8) :: esati                ! saturation vapor pressure over ice (Pa)
    real(r8) :: a0,a1,a2,a3,a4,a5,a6 ! coefficients for esat over water
    real(r8) :: b0,b1,b2,b3,b4,b5,b6 ! coefficients for esat over ice
    real(r8) :: tdc, t               ! Kelvins to Celcius function and its input
    character(len=32) :: fname       ! name of field that is NaN
    character(len=32), parameter :: sub = 'lnd_import'

    ! Constants to compute vapor pressure
    parameter (a0=6.107799961_r8    , a1=4.436518521e-01_r8, &
         a2=1.428945805e-02_r8, a3=2.650648471e-04_r8, &
         a4=3.031240396e-06_r8, a5=2.034080948e-08_r8, &
         a6=6.136820929e-11_r8)

    parameter (b0=6.109177956_r8    , b1=5.034698970e-01_r8, &
         b2=1.886013408e-02_r8, b3=4.176223716e-04_r8, &
         b4=5.824720280e-06_r8, b5=4.838803174e-08_r8, &
         b6=1.838826904e-10_r8)
    !
    ! function declarations
    !
    tdc(t) = min( 50._r8, max(-50._r8,(t-SHR_CONST_TKFRZ)) )
    esatw(t) = 100._r8*(a0+t*(a1+t*(a2+t*(a3+t*(a4+t*(a5+t*a6))))))
    esati(t) = 100._r8*(b0+t*(b1+t*(b2+t*(b3+t*(b4+t*(b5+t*b6))))))
    !---------------------------------------------------------------------------

    ! Note that the precipitation fluxes received  from the coupler
    ! are in units of kg/s/m^2. To convert these precipitation rates
    ! in units of mm/sec, one must divide by 1000 kg/m^3 and multiply
    ! by 1000 mm/m resulting in an overall factor of unity.
    ! Below the units are therefore given in mm/s.


    do g = bounds%begg,bounds%endg
       i = 1 + (g - bounds%begg)

       ! Determine required receive fields

       atm2lnd_inst%forc_hgt_grc(g)                  = x2l(index_x2l_Sa_z,i)         ! zgcmxy  Atm state m
       atm2lnd_inst%forc_u_grc(g)                    = x2l(index_x2l_Sa_u,i)         ! forc_uxy  Atm state m/s
       atm2lnd_inst%forc_v_grc(g)                    = x2l(index_x2l_Sa_v,i)         ! forc_vxy  Atm state m/s
       atm2lnd_inst%forc_solad_grc(g,2)              = x2l(index_x2l_Faxa_swndr,i)   ! forc_sollxy  Atm flux  W/m^2
       atm2lnd_inst%forc_solad_grc(g,1)              = x2l(index_x2l_Faxa_swvdr,i)   ! forc_solsxy  Atm flux  W/m^2
       atm2lnd_inst%forc_solai_grc(g,2)              = x2l(index_x2l_Faxa_swndf,i)   ! forc_solldxy Atm flux  W/m^2
       atm2lnd_inst%forc_solai_grc(g,1)              = x2l(index_x2l_Faxa_swvdf,i)   ! forc_solsdxy Atm flux  W/m^2

       atm2lnd_inst%forc_q_not_downscaled_grc(g)     = x2l(index_x2l_Sa_shum,i)      ! forc_qxy  Atm state kg/kg
       atm2lnd_inst%forc_pbot_not_downscaled_grc(g)  = x2l(index_x2l_Sa_pbot,i)      ! ptcmxy  Atm state Pa
       atm2lnd_inst%forc_t_not_downscaled_grc(g)     = x2l(index_x2l_Sa_tbot,i)      ! forc_txy  Atm state K
       atm2lnd_inst%forc_lwrad_not_downscaled_grc(g) = x2l(index_x2l_Faxa_lwdn,i)    ! flwdsxy Atm flux  W/m^2

       forc_rainc                                    = x2l(index_x2l_Faxa_rainc,i)   ! mm/s
       forc_rainl                                    = x2l(index_x2l_Faxa_rainl,i)   ! mm/s
       forc_snowc                                    = x2l(index_x2l_Faxa_snowc,i)   ! mm/s
       forc_snowl                                    = x2l(index_x2l_Faxa_snowl,i)   ! mm/s

       ! Determine derived quantities for required fields

       forc_t = atm2lnd_inst%forc_t_not_downscaled_grc(g)
       forc_q = atm2lnd_inst%forc_q_not_downscaled_grc(g)
       forc_pbot = atm2lnd_inst%forc_pbot_not_downscaled_grc(g)
       
       atm2lnd_inst%forc_hgt_u_grc(g) = atm2lnd_inst%forc_hgt_grc(g)    !observational height of wind [m]
       atm2lnd_inst%forc_hgt_t_grc(g) = atm2lnd_inst%forc_hgt_grc(g)    !observational height of temperature [m]
       atm2lnd_inst%forc_hgt_q_grc(g) = atm2lnd_inst%forc_hgt_grc(g)    !observational height of humidity [m]
       atm2lnd_inst%forc_vp_grc(g)    = forc_q * forc_pbot  / (0.622_r8 + 0.378_r8 * forc_q)
       atm2lnd_inst%forc_rho_not_downscaled_grc(g) = &
            (forc_pbot - 0.378_r8 * atm2lnd_inst%forc_vp_grc(g)) / (rair * forc_t)
       atm2lnd_inst%forc_wind_grc(g)  = sqrt(atm2lnd_inst%forc_u_grc(g)**2 + atm2lnd_inst%forc_v_grc(g)**2)
       atm2lnd_inst%forc_solar_grc(g) = atm2lnd_inst%forc_solad_grc(g,1) + atm2lnd_inst%forc_solai_grc(g,1) + &
                                        atm2lnd_inst%forc_solad_grc(g,2) + atm2lnd_inst%forc_solai_grc(g,2)

       atm2lnd_inst%forc_rain_not_downscaled_grc(g)  = forc_rainc + forc_rainl
       atm2lnd_inst%forc_snow_not_downscaled_grc(g)  = forc_snowc + forc_snowl

       if (forc_t > SHR_CONST_TKFRZ) then
          e = esatw(tdc(forc_t))
       else
          e = esati(tdc(forc_t))
       end if
       qsat           = 0.622_r8*e / (forc_pbot - 0.378_r8*e)

       !modify specific humidity if precip occurs
       if(1==2) then
          if((forc_rainc+forc_rainl) > 0._r8) then
             forc_q = 0.95_r8*qsat
             !           forc_q = qsat
             atm2lnd_inst%forc_q_not_downscaled_grc(g) = forc_q
          endif
       endif

       ! Check that solar, specific-humidity and LW downward aren't negative
       if ( atm2lnd_inst%forc_lwrad_not_downscaled_grc(g) <= 0.0_r8 )then
          call endrun( sub//' ERROR: Longwave down sent from the atmosphere model is negative or zero' )
       end if
       if ( (atm2lnd_inst%forc_solad_grc(g,1) < 0.0_r8) .or.  (atm2lnd_inst%forc_solad_grc(g,2) < 0.0_r8) &
       .or. (atm2lnd_inst%forc_solai_grc(g,1) < 0.0_r8) .or.  (atm2lnd_inst%forc_solai_grc(g,2) < 0.0_r8) ) then
          call endrun( sub//' ERROR: One of the solar fields (indirect/diffuse, vis or near-IR)'// &
                       ' from the atmosphere model is negative or zero' )
       end if
       if ( atm2lnd_inst%forc_q_not_downscaled_grc(g) < 0.0_r8 )then
          call endrun( sub//' ERROR: Bottom layer specific humidty sent from the atmosphere model is less than zero' )
       end if

       ! Check if any input from the coupler is NaN
       if ( any(isnan(x2l(:,i))) )then
          write(iulog,*) '# of NaNs = ', count(isnan(x2l(:,i)))
          write(iulog,*) 'Which are NaNs = ', isnan(x2l(:,i))
          do k = 1, size(x2l(:,i))
             if ( isnan(x2l(k,i)) )then
                call shr_string_listGetName( seq_flds_x2l_fields, k, fname )
                write(iulog,*) trim(fname)
             end if
          end do
          write(iulog,*) 'gridcell index = ', g
          call endrun( sub//' ERROR: One or more of the input from the atmosphere model are NaN '// &
                       '(Not a Number from a bad floating point calculation)' )
       end if

    end do

  end subroutine lnd_import

  !===============================================================================

  subroutine lnd_export( bounds, lnd2atm_inst, l2x)

    !---------------------------------------------------------------------------
    ! !DESCRIPTION:
    ! Convert the data to be sent from the clm model to the coupler 
    ! 
    ! !USES:
    use shr_kind_mod       , only : r8 => shr_kind_r8
    use seq_flds_mod       , only : seq_flds_l2x_fields
    use clm_varctl         , only : iulog
    use clm_time_manager   , only : get_nstep, get_step_size  
    use domainMod          , only : ldomain
    use shr_string_mod     , only : shr_string_listGetName
    use shr_infnan_mod     , only : isnan => shr_infnan_isnan
    !
    ! !ARGUMENTS:
    implicit none
    type(bounds_type) , intent(in)    :: bounds  ! bounds
    type(lnd2atm_type), intent(inout) :: lnd2atm_inst ! clm land to atmosphere exchange data type
    real(r8)          , intent(out)   :: l2x(:,:)! land to coupler export state on land grid
    !
    ! !LOCAL VARIABLES:
    integer  :: g,i,k ! indices
    integer  :: ier   ! error status
    integer  :: nstep ! time step index
    integer  :: dtime ! time step   
    integer  :: num   ! counter
    character(len=32) :: fname       ! name of field that is NaN
    character(len=32), parameter :: sub = 'lnd_export'
    !---------------------------------------------------------------------------

    ! cesm sign convention is that fluxes are positive downward

    l2x(:,:) = 0.0_r8

    do g = bounds%begg,bounds%endg
       i = 1 + (g-bounds%begg)
       l2x(index_l2x_Sl_t,i)        =  lnd2atm_inst%t_rad_grc(g)
       l2x(index_l2x_Sl_snowh,i)    =  lnd2atm_inst%h2osno_grc(g)
       l2x(index_l2x_Sl_avsdr,i)    =  lnd2atm_inst%albd_grc(g,1)
       l2x(index_l2x_Sl_anidr,i)    =  lnd2atm_inst%albd_grc(g,2)
       l2x(index_l2x_Sl_avsdf,i)    =  lnd2atm_inst%albi_grc(g,1)
       l2x(index_l2x_Sl_anidf,i)    =  lnd2atm_inst%albi_grc(g,2)
       l2x(index_l2x_Sl_tref,i)     =  lnd2atm_inst%t_ref2m_grc(g)
       l2x(index_l2x_Sl_qref,i)     =  lnd2atm_inst%q_ref2m_grc(g)
       l2x(index_l2x_Sl_u10,i)      =  lnd2atm_inst%u_ref10m_grc(g)
       l2x(index_l2x_Fall_taux,i)   = -lnd2atm_inst%taux_grc(g)
       l2x(index_l2x_Fall_tauy,i)   = -lnd2atm_inst%tauy_grc(g)
       l2x(index_l2x_Fall_lat,i)    = -lnd2atm_inst%eflx_lh_tot_grc(g)
       l2x(index_l2x_Fall_sen,i)    = -lnd2atm_inst%eflx_sh_tot_grc(g)
       l2x(index_l2x_Fall_lwup,i)   = -lnd2atm_inst%eflx_lwrad_out_grc(g)
       l2x(index_l2x_Fall_evap,i)   = -lnd2atm_inst%qflx_evap_tot_grc(g)
       l2x(index_l2x_Fall_swnet,i)  =  lnd2atm_inst%fsa_grc(g)

       ! Additional fields for DUST, PROGSSLT, dry-deposition
       ! These are now standard fields, but the check on the index makes sure the driver handles them
       if (index_l2x_Sl_ram1      /= 0 )  l2x(index_l2x_Sl_ram1,i)     =  lnd2atm_inst%ram1_grc(g)
       if (index_l2x_Sl_fv        /= 0 )  l2x(index_l2x_Sl_fv,i)       =  lnd2atm_inst%fv_grc(g)
       if (index_l2x_Sl_soilw     /= 0 )  l2x(index_l2x_Sl_soilw,i)    =  0.5_r8
       if (index_l2x_Fall_flxdst1 /= 0 )  l2x(index_l2x_Fall_flxdst1,i)= -lnd2atm_inst%flxdst_grc(g,1)
       if (index_l2x_Fall_flxdst2 /= 0 )  l2x(index_l2x_Fall_flxdst2,i)= -lnd2atm_inst%flxdst_grc(g,2)
       if (index_l2x_Fall_flxdst3 /= 0 )  l2x(index_l2x_Fall_flxdst3,i)= -lnd2atm_inst%flxdst_grc(g,3)
       if (index_l2x_Fall_flxdst4 /= 0 )  l2x(index_l2x_Fall_flxdst4,i)= -lnd2atm_inst%flxdst_grc(g,4)

       ! sign convention is positive downward with 
       ! hierarchy of atm/glc/lnd/rof/ice/ocn.  
       ! I.e. water sent from land to rof is positive

       !  surface runoff is the sum of qflx_over, qflx_h2osfc_surf
       l2x(index_l2x_Flrl_rofsur,i) = lnd2atm_inst%qflx_rofliq_qsur_grc(g) &
            + lnd2atm_inst%qflx_rofliq_h2osfc_grc(g)

       !  subsurface runoff is the sum of qflx_drain and qflx_perched_drain
       l2x(index_l2x_Flrl_rofsub,i) = lnd2atm_inst%qflx_rofliq_qsub_grc(g) &
            + lnd2atm_inst%qflx_rofliq_drain_perched_grc(g)

       !  qgwl sent individually to coupler
       l2x(index_l2x_Flrl_rofgwl,i) = lnd2atm_inst%qflx_rofliq_qgwl_grc(g)

       ! ice  sent individually to coupler
       l2x(index_l2x_Flrl_rofi,i) = lnd2atm_inst%qflx_rofice_grc(g)

       ! Check if any output sent to the coupler is NaN
       if ( any(isnan(l2x(:,i))) )then
          write(iulog,*) '# of NaNs = ', count(isnan(l2x(:,i)))
          write(iulog,*) 'Which are NaNs = ', isnan(l2x(:,i))
          do k = 1, size(l2x(:,i))
             if ( isnan(l2x(k,i)) )then
                call shr_string_listGetName( seq_flds_l2x_fields, k, fname )
                write(iulog,*) trim(fname)
             end if
          end do
          write(iulog,*) 'gridcell index = ', g
          call endrun( sub//' ERROR: One or more of the output from CLM to the coupler are NaN ' )
       end if

    end do

  end subroutine lnd_export

end module lnd_import_export
