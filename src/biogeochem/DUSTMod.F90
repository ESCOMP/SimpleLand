module DUSTMod

  !----------------------------------------------------------------------- 
  ! !DESCRIPTION: 
  ! Routines in this module calculate Dust mobilization and dry deposition for dust.
  ! Simulates dust mobilization due to wind from the surface into the 
  ! lowest atmospheric layer. On output flx_mss_vrt_dst(ndst) is the surface dust 
  ! emission (kg/m**2/s) [ + = to atm].
  ! Calculates the turbulent component of dust dry deposition, (the turbulent deposition 
  ! velocity through the lowest atmospheric layer). CAM will calculate the settling 
  ! velocity through the whole atmospheric column. The two calculations will determine 
  ! the dust dry deposition flux to the surface.
  !                              
  ! !USES:
  use shr_kind_mod         , only : r8 => shr_kind_r8 
  use shr_log_mod          , only : errMsg => shr_log_errMsg
  use shr_infnan_mod       , only : nan => shr_infnan_nan, assignment(=)
  use clm_varpar           , only : dst_src_nbr, ndst, sz_nbr
  use clm_varcon           , only : grav, spval
  use clm_varctl           , only : iulog
  use abortutils           , only : endrun
  use decompMod            , only : bounds_type
  use LandunitType         , only : lun
  use ColumnType           , only : col
  use PatchType            , only : patch
  !  
  ! !PUBLIC TYPES
  implicit none
  private
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  !
  ! !PUBLIC DATA:
  !
  real(r8) , allocatable :: ovr_src_snk_mss(:,:)
  real(r8) , allocatable :: dmt_vwr(:) ![m] Mass-weighted mean diameter resolved
  real(r8) , allocatable :: stk_crc(:) ![frc] Correction to Stokes settling velocity
  real(r8) tmp1                        !Factor in saltation computation (named as in Charlie's code)
  real(r8) dns_aer                     ![kg m-3] Aerosol density
  !
  ! !PUBLIC DATA TYPES:
  !
  type, public :: dust_type

     real(r8), pointer, PUBLIC  :: flx_mss_vrt_dst_patch     (:,:) ! surface dust emission (kg/m**2/s) [ + = to atm] (ndst) 
     real(r8), pointer, private :: flx_mss_vrt_dst_tot_patch (:)   ! total dust flux into atmosphere
     real(r8), pointer, private :: vlc_trb_patch             (:,:) ! turbulent deposition velocity  (m/s) (ndst) 
     real(r8), pointer, private :: vlc_trb_1_patch           (:)   ! turbulent deposition velocity 1(m/s)
     real(r8), pointer, private :: vlc_trb_2_patch           (:)   ! turbulent deposition velocity 2(m/s)
     real(r8), pointer, private :: vlc_trb_3_patch           (:)   ! turbulent deposition velocity 3(m/s)
     real(r8), pointer, private :: vlc_trb_4_patch           (:)   ! turbulent deposition velocity 4(m/s)
     real(r8), pointer, private :: mbl_bsn_fct_col           (:)   ! basin factor

   contains

     procedure , public  :: Init
     procedure , private :: InitAllocate 
     procedure , private :: InitHistory  
     procedure , private :: InitCold     
     procedure , private :: InitDustVars ! Initialize variables used in subroutine Dust

  end type dust_type
  !------------------------------------------------------------------------

  character(len=*), parameter, private :: sourcefile = &
       __FILE__

contains

  !------------------------------------------------------------------------
  subroutine Init(this, bounds)

    class(dust_type) :: this
    type(bounds_type), intent(in) :: bounds  

    call this%InitAllocate (bounds)
    call this%InitHistory  (bounds)
    call this%InitCold     (bounds)
    call this%InitDustVars (bounds)

  end subroutine Init

  !------------------------------------------------------------------------
  subroutine InitAllocate(this, bounds)
    !
    ! !ARGUMENTS:
    class (dust_type) :: this
    type(bounds_type), intent(in) :: bounds  
    !
    ! !LOCAL VARIABLES:
    integer :: begp,endp
    integer :: begc,endc
    !------------------------------------------------------------------------

    begp = bounds%begp ; endp = bounds%endp
    begc = bounds%begc ; endc = bounds%endc

    allocate(this%flx_mss_vrt_dst_patch     (begp:endp,1:ndst)) ; this%flx_mss_vrt_dst_patch     (:,:) = nan
    allocate(this%flx_mss_vrt_dst_tot_patch (begp:endp))        ; this%flx_mss_vrt_dst_tot_patch (:)   = nan
    allocate(this%vlc_trb_patch             (begp:endp,1:ndst)) ; this%vlc_trb_patch             (:,:) = nan
    allocate(this%vlc_trb_1_patch           (begp:endp))        ; this%vlc_trb_1_patch           (:)   = nan
    allocate(this%vlc_trb_2_patch           (begp:endp))        ; this%vlc_trb_2_patch           (:)   = nan 
    allocate(this%vlc_trb_3_patch           (begp:endp))        ; this%vlc_trb_3_patch           (:)   = nan
    allocate(this%vlc_trb_4_patch           (begp:endp))        ; this%vlc_trb_4_patch           (:)   = nan
    allocate(this%mbl_bsn_fct_col           (begc:endc))        ; this%mbl_bsn_fct_col     (:)   = nan

  end subroutine InitAllocate

  !------------------------------------------------------------------------
  subroutine InitHistory(this, bounds)
    !
    ! !USES:
    use histFileMod, only : hist_addfld1d
    !
    !
    ! !ARGUMENTS:
    class (dust_type) :: this
    type(bounds_type), intent(in) :: bounds  
    !
    ! !LOCAL VARIABLES:
    integer :: begp,endp
    !------------------------------------------------------------------------

    begp = bounds%begp; endp = bounds%endp

    this%flx_mss_vrt_dst_tot_patch(begp:endp) = spval
    call hist_addfld1d (fname='DSTFLXT', units='kg/m2/s',  &
         avgflag='A', long_name='total surface dust emission', &
         ptr_patch=this%flx_mss_vrt_dst_tot_patch, set_lake=0._r8, set_urb=0._r8, default='inactive')

    this%vlc_trb_1_patch(begp:endp) = spval
    call hist_addfld1d (fname='DPVLTRB1', units='m/s',  &
         avgflag='A', long_name='turbulent deposition velocity 1', &
         ptr_patch=this%vlc_trb_1_patch, default='inactive')

    this%vlc_trb_2_patch(begp:endp) = spval
    call hist_addfld1d (fname='DPVLTRB2', units='m/s',  &
         avgflag='A', long_name='turbulent deposition velocity 2', &
         ptr_patch=this%vlc_trb_2_patch, default='inactive')

    this%vlc_trb_3_patch(begp:endp) = spval
    call hist_addfld1d (fname='DPVLTRB3', units='m/s',  &
         avgflag='A', long_name='turbulent deposition velocity 3', &
         ptr_patch=this%vlc_trb_3_patch, default='inactive')

    this%vlc_trb_4_patch(begp:endp) = spval
    call hist_addfld1d (fname='DPVLTRB4', units='m/s',  &
         avgflag='A', long_name='turbulent deposition velocity 4', &
         ptr_patch=this%vlc_trb_4_patch, default='inactive')

  end subroutine InitHistory

  !-----------------------------------------------------------------------
  subroutine InitCold(this, bounds)
    !
    ! !ARGUMENTS:
    class (dust_type) :: this
    type(bounds_type), intent(in) :: bounds  
    !
    ! !LOCAL VARIABLES:
    integer :: c,l
    !-----------------------------------------------------------------------

    ! Set basin factor to 1 for now

    do c = bounds%begc, bounds%endc
       l = col%landunit(c)

       if (.not.lun%lakpoi(l)) then
          this%mbl_bsn_fct_col(c) = 1.0_r8
       end if
    end do

  end subroutine InitCold

   !------------------------------------------------------------------------
   subroutine InitDustVars(this, bounds)
     !
     ! !DESCRIPTION: 
     !
     ! Compute source efficiency factor from topography
     ! Initialize other variables used in subroutine Dust:
     ! ovr_src_snk_mss(m,n) and tmp1.
     ! Define particle diameter and density needed by atm model
     ! as well as by dry dep model
     ! Source: Paul Ginoux (for source efficiency factor)
     ! Modifications by C. Zender and later by S. Levis
     ! Rest of subroutine from C. Zender's dust model
     !
     ! !USES
     use shr_const_mod , only: SHR_CONST_PI, SHR_CONST_RDAIR
     use shr_spfn_mod  , only: erf => shr_spfn_erf
     use decompMod     , only : get_proc_bounds
     !
     ! !ARGUMENTS:
     class(dust_type)  :: this
     type(bounds_type), intent(in) :: bounds  
     !
     ! !LOCAL VARIABLES
    integer  :: fc,c,l,m,n              ! indices
    real(r8) :: ovr_src_snk_frc
    real(r8) :: sqrt2lngsdi             ! [frc] Factor in erf argument
    real(r8) :: lndmaxjovrdmdni         ! [frc] Factor in erf argument
    real(r8) :: lndminjovrdmdni         ! [frc] Factor in erf argument
    real(r8) :: ryn_nbr_frc_thr_prx_opt ! [frc] Threshold friction Reynolds number approximation for optimal size
    real(r8) :: ryn_nbr_frc_thr_opt_fnc ! [frc] Threshold friction Reynolds factor for saltation calculation
    real(r8) :: icf_fct                 ! Interpartical cohesive forces factor for saltation calc
    real(r8) :: dns_fct                 ! Density ratio factor for saltation calculation
    real(r8) :: dmt_min(ndst)           ! [m] Size grid minimum
    real(r8) :: dmt_max(ndst)           ! [m] Size grid maximum
    real(r8) :: dmt_ctr(ndst)           ! [m] Diameter at bin center
    real(r8) :: dmt_dlt(ndst)           ! [m] Width of size bin
    real(r8) :: slp_crc(ndst)           ! [frc] Slip correction factor
    real(r8) :: vlm_rsl(ndst)           ! [m3 m-3] Volume concentration resolved
    real(r8) :: vlc_stk(ndst)           ! [m s-1] Stokes settling velocity
    real(r8) :: vlc_grv(ndst)           ! [m s-1] Settling velocity
    real(r8) :: ryn_nbr_grv(ndst)       ! [frc] Reynolds number at terminal velocity
    real(r8) :: cff_drg_grv(ndst)       ! [frc] Drag coefficient at terminal velocity
    real(r8) :: tmp                     ! temporary 
    real(r8) :: ln_gsd                  ! [frc] ln(gsd)
    real(r8) :: gsd_anl                 ! [frc] Geometric standard deviation
    real(r8) :: dmt_vma                 ! [m] Mass median diameter analytic She84 p.75 Tabl.1
    real(r8) :: dmt_nma                 ! [m] Number median particle diameter
    real(r8) :: lgn_dst                 ! Lognormal distribution at sz_ctr
    real(r8) :: eps_max                 ! [frc] Relative accuracy for convergence
    real(r8) :: eps_crr                 ! [frc] Current relative accuracy
    real(r8) :: itr_idx                 ! [idx] Counting index
    real(r8) :: dns_mdp                 ! [kg m-3] Midlayer density
    real(r8) :: mfp_atm                 ! [m] Mean free path of air
    real(r8) :: vsc_dyn_atm             ! [kg m-1 s-1] Dynamic viscosity of air
    real(r8) :: vsc_knm_atm             ! [kg m-1 s-1] Kinematic viscosity of air
    real(r8) :: vlc_grv_old             ! [m s-1] Previous gravitational settling velocity
    real(r8) :: series_ratio            ! Factor for logarithmic grid
    real(r8) :: lngsdsqrttwopi_rcp      ! Factor in lognormal distribution
    real(r8) :: sz_min(sz_nbr)          ! [m] Size Bin minima
    real(r8) :: sz_max(sz_nbr)          ! [m] Size Bin maxima
    real(r8) :: sz_ctr(sz_nbr)          ! [m] Size Bin centers
    real(r8) :: sz_dlt(sz_nbr)          ! [m] Size Bin widths
    
    ! constants
    real(r8), allocatable :: dmt_vma_src(:) ! [m] Mass median diameter       BSM96 p. 73 Table 2
    real(r8), allocatable :: gsd_anl_src(:) ! [frc] Geometric std deviation  BSM96 p. 73 Table 2
    real(r8), allocatable :: mss_frc_src(:) ! [frc] Mass fraction            BSM96 p. 73 Table 2

    real(r8) :: dmt_grd(5) =                  &     ! [m] Particle diameter grid
         (/ 0.1e-6_r8, 1.0e-6_r8, 2.5e-6_r8, 5.0e-6_r8, 10.0e-6_r8 /)
    real(r8), parameter :: dmt_slt_opt = 75.0e-6_r8    ! [m] Optim diam for saltation
    real(r8), parameter :: dns_slt = 2650.0_r8         ! [kg m-3] Density of optimal saltation particles
    !------------------------------------------------------------------------

    associate(& 
         mbl_bsn_fct  =>  this%mbl_bsn_fct_col   & ! Output:  [real(r8) (:)] basin factor                                       
         )

      ! allocate module variable
      allocate (ovr_src_snk_mss(dst_src_nbr,ndst))  
      allocate (dmt_vwr(ndst))
      allocate (stk_crc(ndst))

      ! allocate local variable
      allocate (dmt_vma_src(dst_src_nbr))  
      allocate (gsd_anl_src(dst_src_nbr))  
      allocate (mss_frc_src(dst_src_nbr))  

      dmt_vma_src(:) = (/ 0.832e-6_r8 , 4.82e-6_r8 , 19.38e-6_r8 /)        
      gsd_anl_src(:) = (/ 2.10_r8     , 1.90_r8    , 1.60_r8     /)        
      mss_frc_src(:) = (/ 0.036_r8    , 0.957_r8   , 0.007_r8 /)                  

      ! the following comes from (1) szdstlgn.F subroutine ovr_src_snk_frc_get
      !                      and (2) dstszdst.F subroutine dst_szdst_ini
      ! purpose(1): given one set (the "source") of lognormal distributions,
      !             and one set of bin boundaries (the "sink"), compute and return
      !             the overlap factors between the source and sink distributions
      ! purpose(2): set important statistics of size distributions

      do m = 1, dst_src_nbr
         sqrt2lngsdi = sqrt(2.0_r8) * log(gsd_anl_src(m))
         do n = 1, ndst
            lndmaxjovrdmdni = log(dmt_grd(n+1)/dmt_vma_src(m))
            lndminjovrdmdni = log(dmt_grd(n  )/dmt_vma_src(m))
            ovr_src_snk_frc = 0.5_r8 * (erf(lndmaxjovrdmdni/sqrt2lngsdi) - &
                 erf(lndminjovrdmdni/sqrt2lngsdi))
            ovr_src_snk_mss(m,n) = ovr_src_snk_frc * mss_frc_src(m)
         end do
      end do

      ! The following code from subroutine wnd_frc_thr_slt_get was placed 
      ! here because tmp1 needs to be defined just once

      ryn_nbr_frc_thr_prx_opt = 0.38_r8 + 1331.0_r8 * (100.0_r8*dmt_slt_opt)**1.56_r8

      if (ryn_nbr_frc_thr_prx_opt < 0.03_r8) then
         write(iulog,*) 'dstmbl: ryn_nbr_frc_thr_prx_opt < 0.03'
         call endrun(msg=errMsg(sourcefile, __LINE__))
      else if (ryn_nbr_frc_thr_prx_opt < 10.0_r8) then
         ryn_nbr_frc_thr_opt_fnc = -1.0_r8 + 1.928_r8 * (ryn_nbr_frc_thr_prx_opt**0.0922_r8)
         ryn_nbr_frc_thr_opt_fnc = 0.1291_r8 * 0.1291_r8 / ryn_nbr_frc_thr_opt_fnc
      else
         ryn_nbr_frc_thr_opt_fnc = 1.0_r8 - 0.0858_r8 * exp(-0.0617_r8*(ryn_nbr_frc_thr_prx_opt-10.0_r8))
         ryn_nbr_frc_thr_opt_fnc = 0.120_r8 * 0.120_r8 * ryn_nbr_frc_thr_opt_fnc * ryn_nbr_frc_thr_opt_fnc
      end if

      icf_fct = 1.0_r8 + 6.0e-07_r8 / (dns_slt * grav * (dmt_slt_opt**2.5_r8))
      dns_fct = dns_slt * grav * dmt_slt_opt
      tmp1 = sqrt(icf_fct * dns_fct * ryn_nbr_frc_thr_opt_fnc)

      ! Introducing particle diameter. Needed by atm model and by dry dep model.
      ! Taken from Charlie Zender's subroutines dst_psd_ini, dst_sz_rsl,
      ! grd_mk (dstpsd.F90) and subroutine lgn_evl (psdlgn.F90)

      ! Charlie allows logarithmic or linear option for size distribution
      ! however, he hardwires the distribution to logarithmic in his code
      ! therefore, I take his logarithmic code only
      ! furthermore, if dst_nbr == 4, he overrides the automatic grid calculation
      ! he currently works with dst_nbr = 4, so I only take the relevant code
      ! if ndst ever becomes different from 4, must add call grd_mk (dstpsd.F90)
      ! as done in subroutine dst_psd_ini
      ! note that here ndst = dst_nbr

      ! Override automatic grid with preset grid if available

      if (ndst == 4) then
         do n = 1, ndst
            dmt_min(n) = dmt_grd(n)                       ![m] Max diameter in bin
            dmt_max(n) = dmt_grd(n+1)                     ![m] Min diameter in bin
            dmt_ctr(n) = 0.5_r8 * (dmt_min(n)+dmt_max(n)) ![m] Diameter at bin ctr
            dmt_dlt(n) = dmt_max(n)-dmt_min(n)            ![m] Width of size bin
         end do
      else
         write(iulog,*) 'Dustini error: ndst must equal to 4 with current code'
         call endrun(msg=errMsg(sourcefile, __LINE__))
         !see more comments above end if ndst == 4
      end if

      ! Bin physical properties

      gsd_anl = 2.0_r8      ! [frc] Geometric std dev PaG77 p. 2080 Table1
      ln_gsd = log(gsd_anl)
      dns_aer = 2.5e+3_r8   ! [kg m-3] Aerosol density

      ! Set a fundamental statistic for each bin

      dmt_vma = 3.5000e-6_r8 ! [m] Mass median diameter analytic She84 p.75 Table1

      ! Compute analytic size statistics
      ! Convert mass median diameter to number median diameter (call vma2nma)

      dmt_nma = dmt_vma * exp(-3.0_r8*ln_gsd*ln_gsd) ! [m]

      ! Compute resolved size statistics for each size distribution
      ! In C. Zender's code call dst_sz_rsl

      do n = 1, ndst

         series_ratio = (dmt_max(n)/dmt_min(n))**(1.0_r8/sz_nbr)
         sz_min(1) = dmt_min(n)
         do m = 2, sz_nbr                            ! Loop starts at 2
            sz_min(m) = sz_min(m-1) * series_ratio
         end do

         ! Derived grid values
         do m = 1, sz_nbr-1                          ! Loop ends at sz_nbr-1
            sz_max(m) = sz_min(m+1)                  ! [m]
         end do
         sz_max(sz_nbr) = dmt_max(n)                 ! [m]

         ! Final derived grid values
         do m = 1, sz_nbr
            sz_ctr(m) = 0.5_r8 * (sz_min(m)+sz_max(m))
            sz_dlt(m) = sz_max(m)-sz_min(m)
         end do

         lngsdsqrttwopi_rcp = 1.0_r8 / (ln_gsd*sqrt(2.0_r8*SHR_CONST_PI))
         dmt_vwr(n) = 0.0_r8 ! [m] Mass wgted diameter resolved
         vlm_rsl(n) = 0.0_r8 ! [m3 m-3] Volume concentration resolved

         do m = 1, sz_nbr

            ! Evaluate lognormal distribution for these sizes (call lgn_evl)
            tmp = log(sz_ctr(m)/dmt_nma) / ln_gsd
            lgn_dst = lngsdsqrttwopi_rcp * exp(-0.5_r8*tmp*tmp) / sz_ctr(m)

            ! Integrate moments of size distribution
            dmt_vwr(n) = dmt_vwr(n) + sz_ctr(m) *                    &
                 SHR_CONST_PI / 6.0_r8 * (sz_ctr(m)**3.0_r8) * & ![m3] Volume
                 lgn_dst * sz_dlt(m)                ![# m-3] Number concentrn
            vlm_rsl(n) = vlm_rsl(n) +                                &
                 SHR_CONST_PI / 6.0_r8 * (sz_ctr(m)**3.0_r8) * & ![m3] Volume
                 lgn_dst * sz_dlt(m)                ![# m-3] Number concentrn

         end do

         dmt_vwr(n) = dmt_vwr(n) / vlm_rsl(n) ![m] Mass weighted diameter resolved

      end do

      ! calculate correction to Stokes' settling velocity (subroutine stk_crc_get)

      eps_max = 1.0e-4_r8
      dns_mdp = 100000._r8 / (295.0_r8*SHR_CONST_RDAIR) ![kg m-3] const prs_mdp & tpt_vrt

      ! Size-independent thermokinetic properties

      vsc_dyn_atm = 1.72e-5_r8 * ((295.0_r8/273.0_r8)**1.5_r8) * 393.0_r8 / &
           (295.0_r8+120.0_r8)      ![kg m-1 s-1] RoY94 p.102 tpt_mdp=295.0
      mfp_atm = 2.0_r8 * vsc_dyn_atm / &  !SeP97 p. 455 constant prs_mdp, tpt_mdp
           (100000._r8*sqrt(8.0_r8/(SHR_CONST_PI*SHR_CONST_RDAIR*295.0_r8)))
      vsc_knm_atm = vsc_dyn_atm / dns_mdp ![m2 s-1] Kinematic viscosity of air

      do m = 1, ndst
         slp_crc(m) = 1.0_r8 + 2.0_r8 * mfp_atm *                      &
              (1.257_r8+0.4_r8*exp(-1.1_r8*dmt_vwr(m)/(2.0_r8*mfp_atm))) / &
              dmt_vwr(m)                      ! [frc] Slip correction factor SeP97 p.464
         vlc_stk(m) = (1.0_r8/18.0_r8) * dmt_vwr(m) * dmt_vwr(m) * dns_aer * &
              grav * slp_crc(m) / vsc_dyn_atm ! [m s-1] SeP97 p.466
      end do

      ! For Reynolds number flows Re < 0.1 Stokes' velocity is valid for
      ! vlc_grv SeP97 p. 466 (8.42). For larger Re, inertial effects become
      ! important and empirical drag coefficients must be employed
      ! Implicit equation for Re, Cd, and Vt is SeP97 p. 467 (8.44)
      ! Using Stokes' velocity rather than iterative solution with empirical
      ! drag coefficient causes 60% errors for D = 200 um SeP97 p. 468

      ! Iterative solution for drag coefficient, Reynolds number, and terminal veloc
      do m = 1, ndst

         ! Initialize accuracy and counter
         eps_crr = eps_max + 1.0_r8  ![frc] Current relative accuracy
         itr_idx = 0                 ![idx] Counting index

         ! Initial guess for vlc_grv is exact for Re < 0.1
         vlc_grv(m) = vlc_stk(m)     ![m s-1]

         do while(eps_crr > eps_max)

            ! Save terminal velocity for convergence test
            vlc_grv_old = vlc_grv(m) ![m s-1]
            ryn_nbr_grv(m) = vlc_grv(m) * dmt_vwr(m) / vsc_knm_atm !SeP97 p.460

            ! Update drag coefficient based on new Reynolds number
            if (ryn_nbr_grv(m) < 0.1_r8) then
               cff_drg_grv(m) = 24.0_r8 / ryn_nbr_grv(m) !Stokes' law Sep97 p.463 (8.32)
            else if (ryn_nbr_grv(m) < 2.0_r8) then
               cff_drg_grv(m) = (24.0_r8/ryn_nbr_grv(m)) *    &
                    (1.0_r8 + 3.0_r8*ryn_nbr_grv(m)/16.0_r8 + &
                    9.0_r8*ryn_nbr_grv(m)*ryn_nbr_grv(m)*     &
                    log(2.0_r8*ryn_nbr_grv(m))/160.0_r8)        !Sep97 p.463 (8.32)
            else if (ryn_nbr_grv(m) < 500.0_r8) then
               cff_drg_grv(m) = (24.0_r8/ryn_nbr_grv(m)) * &
                    (1.0_r8 + 0.15_r8*ryn_nbr_grv(m)**0.687_r8) !Sep97 p.463 (8.32)
            else if (ryn_nbr_grv(m) < 2.0e5_r8) then
               cff_drg_grv(m) = 0.44_r8                         !Sep97 p.463 (8.32)
            else
               write(iulog,'(a,es9.2)') "ryn_nbr_grv(m) = ",ryn_nbr_grv(m)
               write(iulog,*)'Dustini error: Reynolds number too large in stk_crc_get()'
               call endrun(msg=errMsg(sourcefile, __LINE__))
            end if

            ! Update terminal velocity based on new Reynolds number and drag coeff
            ! [m s-1] Terminal veloc SeP97 p.467 (8.44)

            vlc_grv(m) = sqrt(4.0_r8 * grav * dmt_vwr(m) * slp_crc(m) * dns_aer / &
                 (3.0_r8*cff_drg_grv(m)*dns_mdp))   
            eps_crr = abs((vlc_grv(m)-vlc_grv_old)/vlc_grv(m)) !Relative convergence
            if (itr_idx == 12) then
               ! Numerical pingpong may occur when Re = 0.1, 2.0, or 500.0
               ! due to discontinuities in derivative of drag coefficient
               vlc_grv(m) = 0.5_r8 * (vlc_grv(m)+vlc_grv_old)  ! [m s-1]
            end if
            if (itr_idx > 20) then
               write(iulog,*) 'Dustini error: Terminal velocity not converging ',&
                    ' in stk_crc_get(), breaking loop...'
               goto 100                                        !to next iteration
            end if
            itr_idx = itr_idx + 1

         end do                                                !end while

100      continue   !Label to jump to when iteration does not converge
      end do   !end loop over size

      ! Compute factors to convert Stokes' settling velocities to
      ! actual settling velocities

      do m = 1, ndst
         stk_crc(m) = vlc_grv(m) / vlc_stk(m)
      end do

    end associate 

  end subroutine InitDustVars

end module DUSTMod
