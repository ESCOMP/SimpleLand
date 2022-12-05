module TemperatureType

#include "shr_assert.h"

  !------------------------------------------------------------------------------
  use shr_kind_mod    , only : r8 => shr_kind_r8
  use shr_log_mod     , only : errMsg => shr_log_errMsg
  use decompMod       , only : bounds_type
  use abortutils      , only : endrun
  use clm_varctl      , only : iulog
  use clm_varpar      , only : nlevsno, nlevgrnd, nlevlak, nlevlak, nlevurb
  use clm_varcon      , only : spval, ispval
  use GridcellType    , only : grc
  use LandunitType    , only : lun                
  use ColumnType      , only : col                
  use PatchType       , only : patch                
  !
  implicit none
  save
  private
  !
  type, public :: temperature_type

     ! Temperatures
     real(r8), pointer :: t_veg_patch              (:)   ! patch vegetation temperature (Kelvin)
     real(r8), pointer :: t_veg_day_patch          (:)   ! patch daytime  accumulative vegetation temperature (Kelvinx*nsteps), LUNA specific, from midnight to current step
     real(r8), pointer :: t_veg_night_patch        (:)   ! patch night-time accumulative vegetation temperature (Kelvin*nsteps), LUNA specific, from midnight to current step
     real(r8), pointer :: t_veg10_day_patch        (:)   ! 10 day running mean of patch daytime time vegetation temperature (Kelvin), LUNA specific, but can be reused
     real(r8), pointer :: t_veg10_night_patch      (:)   ! 10 day running mean of patch night time vegetation temperature (Kelvin), LUNA specific, but can be reused
     integer,  pointer :: ndaysteps_patch          (:)   ! number of daytime steps accumulated from mid-night, LUNA specific
     integer,  pointer :: nnightsteps_patch        (:)   ! number of nighttime steps accumulated from mid-night, LUNA specific
     real(r8), pointer :: t_h2osfc_col             (:)   ! col surface water temperature
     real(r8), pointer :: t_h2osfc_bef_col         (:)   ! col surface water temperature from time-step before  
     real(r8), pointer :: t_ssbef_col              (:,:) ! col soil/snow temperature before update (-nlevsno+1:nlevgrnd) 
     real(r8), pointer :: t_soisno_col             (:,:) ! col soil temperature (Kelvin)  (-nlevsno+1:nlevgrnd) 
     real(r8), pointer :: t_soi10cm_col            (:)   ! col soil temperature in top 10cm of soil (Kelvin)
     real(r8), pointer :: t_lake_col               (:,:) ! col lake temperature (Kelvin)  (1:nlevlak)          
     real(r8), pointer :: t_grnd_col               (:)   ! col ground temperature (Kelvin)
     real(r8), pointer :: snot_top_col             (:)   ! col temperature of top snow layer [K]
     real(r8), pointer :: dTdz_top_col             (:)   ! col temperature gradient in top layer  [K m-1]
     real(r8), pointer :: dt_veg_patch             (:)   ! patch change in t_veg, last iteration (Kelvin)

     real(r8), pointer :: dt_grnd_col              (:)   ! col change in t_grnd, last iteration (Kelvin)
     real(r8), pointer :: thv_col                  (:)   ! col virtual potential temperature (kelvin)
     real(r8), pointer :: thm_patch                (:)   ! patch intermediate variable (forc_t+0.0098*forc_hgt_t_patch)
     real(r8), pointer :: t_a10_patch              (:)   ! patch 10-day running mean of the 2 m temperature (K)
     real(r8), pointer :: t_a10min_patch           (:)   ! patch 10-day running mean of min 2-m temperature
     real(r8), pointer :: t_a5min_patch            (:)   ! patch 5-day running mean of min 2-m temperature

     real(r8), pointer :: taf_lun                  (:)   ! lun urban canopy air temperature (K)

     real(r8), pointer :: t_ref2m_patch            (:)   ! patch 2 m height surface air temperature (Kelvin)
     real(r8), pointer :: t_ref2m_min_patch        (:)   ! patch daily minimum of average 2 m height surface air temperature (K)
     real(r8), pointer :: t_ref2m_max_patch        (:)   ! patch daily maximum of average 2 m height surface air temperature (K)
     real(r8), pointer :: t_ref2m_min_inst_patch   (:)   ! patch instantaneous daily min of average 2 m height surface air temp (K)
     real(r8), pointer :: t_ref2m_max_inst_patch   (:)   ! patch instantaneous daily max of average 2 m height surface air temp (K)

     ! Accumulated quantities
     !
     ! TODO(wjs, 2014-08-05) Move these to the module(s) where they are used, to improve
     ! modularity. In cases where they are used by two completely different modules,
     ! which only use the same variable out of convenience, introduce a duplicate (point
     ! being: that way one parameterization is free to change the exact meaning of its
     ! accumulator without affecting the other).
     !
     real(r8), pointer :: gdd0_patch              (:)   ! patch growing degree-days base  0C from planting  (ddays)
     real(r8), pointer :: gdd8_patch              (:)   ! patch growing degree-days base  8C from planting  (ddays)
     real(r8), pointer :: gdd10_patch             (:)   ! patch growing degree-days base 10C from planting  (ddays)
     real(r8), pointer :: gdd020_patch            (:)   ! patch 20-year average of gdd0                     (ddays)
     real(r8), pointer :: gdd820_patch            (:)   ! patch 20-year average of gdd8                     (ddays)
     real(r8), pointer :: gdd1020_patch           (:)   ! patch 20-year average of gdd10                    (ddays)

     ! Heat content
     real(r8), pointer :: beta_col                 (:)   ! coefficient of convective velocity [-]
     real(r8), pointer :: heat1_grc                (:)   ! grc initial gridcell total heat content
     real(r8), pointer :: heat2_grc                (:)   ! grc post land cover change total heat content
     real(r8), pointer :: liquid_water_temp1_grc   (:)   ! grc initial weighted average liquid water temperature (K)
     real(r8), pointer :: liquid_water_temp2_grc   (:)   ! grc post land cover change weighted average liquid water temperature (K)

     ! Flags
     integer , pointer :: imelt_col                (:,:) ! flag for melting (=1), freezing (=2), Not=0 (-nlevsno+1:nlevgrnd) 

     ! Emissivities
     real(r8), pointer :: emv_patch                (:)   ! patch vegetation emissivity 

     ! Misc
     real(r8), pointer    :: xmf_col               (:)   ! total latent heat of phase change of ground water
     real(r8), pointer    :: xmf_h2osfc_col        (:)   ! latent heat of phase change of surface water
     real(r8), pointer    :: fact_col              (:,:) ! used in computing tridiagonal matrix
     real(r8), pointer    :: c_h2osfc_col          (:)   ! heat capacity of surface water

   contains

     procedure, public  :: Init
     procedure, public  :: Restart      
     procedure, private :: InitAllocate 
     procedure, private :: InitHistory  
     procedure, private :: InitCold     

  end type temperature_type

  character(len=*), parameter, private :: sourcefile = &
       __FILE__
  !------------------------------------------------------------------------

contains

  !------------------------------------------------------------------------
  subroutine Init(this, bounds)
    !
    ! !DESCRIPTION:
    !
    ! Initialization of the data type. Allocate data, setup variables
    ! for history output, and initialize values needed for a cold-start.
    !
    class(temperature_type)        :: this
    type(bounds_type) , intent(in) :: bounds  

    call this%InitAllocate ( bounds )
    call this%InitHistory ( bounds)
    call this%InitCold ( bounds)

  end subroutine Init

  !------------------------------------------------------------------------
  subroutine InitAllocate(this, bounds)
    !
    ! !DESCRIPTION:
    ! Initialize and allocate data structure
    !
    ! !USES:
    use shr_infnan_mod , only : nan => shr_infnan_nan, assignment(=)
    !
    ! !ARGUMENTS:
    class(temperature_type) :: this
    type(bounds_type), intent(in) :: bounds  
    !
    ! !LOCAL VARIABLES:
    integer :: begp, endp
    integer :: begc, endc
    integer :: begl, endl
    integer :: begg, endg
    !------------------------------------------------------------------------

    begp = bounds%begp; endp= bounds%endp
    begc = bounds%begc; endc= bounds%endc
    begl = bounds%begl; endl= bounds%endl
    begg = bounds%begg; endg= bounds%endg

    ! Temperatures
    allocate(this%t_veg_patch              (begp:endp))                      ; this%t_veg_patch              (:)   = nan
    allocate(this%t_h2osfc_col             (begc:endc))                      ; this%t_h2osfc_col             (:)   = nan
    allocate(this%t_h2osfc_bef_col         (begc:endc))                      ; this%t_h2osfc_bef_col         (:)   = nan
    allocate(this%t_ssbef_col              (begc:endc,-nlevsno+1:nlevgrnd))  ; this%t_ssbef_col              (:,:) = nan
    allocate(this%t_soisno_col             (begc:endc,-nlevsno+1:nlevgrnd))  ; this%t_soisno_col             (:,:) = nan
    allocate(this%t_lake_col               (begc:endc,1:nlevlak))            ; this%t_lake_col               (:,:) = nan
    allocate(this%t_grnd_col               (begc:endc))                      ; this%t_grnd_col               (:)   = nan
    allocate(this%snot_top_col             (begc:endc))                      ; this%snot_top_col             (:)   = nan
    allocate(this%dTdz_top_col             (begc:endc))                      ; this%dTdz_top_col             (:)   = nan
    allocate(this%dt_veg_patch             (begp:endp))                      ; this%dt_veg_patch             (:)   = nan

    allocate(this%t_soi10cm_col            (begc:endc))                      ; this%t_soi10cm_col            (:)   = nan
    allocate(this%dt_grnd_col              (begc:endc))                      ; this%dt_grnd_col              (:)   = nan
    allocate(this%thv_col                  (begc:endc))                      ; this%thv_col                  (:)   = nan
    allocate(this%thm_patch                (begp:endp))                      ; this%thm_patch                (:)   = nan
    allocate(this%t_a10_patch              (begp:endp))                      ; this%t_a10_patch              (:)   = nan
    allocate(this%t_a10min_patch           (begp:endp))                      ; this%t_a10min_patch           (:)   = nan
    allocate(this%t_a5min_patch            (begp:endp))                      ; this%t_a5min_patch            (:)   = nan

    allocate(this%taf_lun                  (begl:endl))                      ; this%taf_lun                  (:)   = nan

    allocate(this%t_ref2m_patch            (begp:endp))                      ; this%t_ref2m_patch            (:)   = nan
    allocate(this%t_ref2m_min_patch        (begp:endp))                      ; this%t_ref2m_min_patch        (:)   = nan
    allocate(this%t_ref2m_max_patch        (begp:endp))                      ; this%t_ref2m_max_patch        (:)   = nan
    allocate(this%t_ref2m_max_inst_patch   (begp:endp))                      ; this%t_ref2m_max_inst_patch   (:)   = nan
    allocate(this%t_ref2m_min_inst_patch   (begp:endp))                      ; this%t_ref2m_min_inst_patch   (:)   = nan

    ! Accumulated fields
    allocate(this%gdd0_patch               (begp:endp))                      ; this%gdd0_patch               (:)   = spval
    allocate(this%gdd8_patch               (begp:endp))                      ; this%gdd8_patch               (:)   = spval
    allocate(this%gdd10_patch              (begp:endp))                      ; this%gdd10_patch              (:)   = spval
    allocate(this%gdd020_patch             (begp:endp))                      ; this%gdd020_patch             (:)   = spval
    allocate(this%gdd820_patch             (begp:endp))                      ; this%gdd820_patch             (:)   = spval
    allocate(this%gdd1020_patch            (begp:endp))                      ; this%gdd1020_patch            (:)   = spval

    ! Heat content
    allocate(this%beta_col                 (begc:endc))                      ; this%beta_col                 (:)   = nan
    allocate(this%heat1_grc                (begg:endg))                      ; this%heat1_grc                (:)   = nan
    allocate(this%heat2_grc                (begg:endg))                      ; this%heat2_grc                (:)   = nan
    allocate(this%liquid_water_temp1_grc   (begg:endg))                      ; this%liquid_water_temp1_grc   (:)   = nan
    allocate(this%liquid_water_temp2_grc   (begg:endg))                      ; this%liquid_water_temp2_grc   (:)   = nan

    ! flags
    allocate(this%imelt_col                (begc:endc,-nlevsno+1:nlevgrnd))  ; this%imelt_col                (:,:) = huge(1)

    ! emissivities
    allocate(this%emv_patch                (begp:endp))                      ; this%emv_patch                (:)   = nan

    allocate(this%xmf_col                  (begc:endc))                      ; this%xmf_col                  (:)   = nan
    allocate(this%xmf_h2osfc_col           (begc:endc))                      ; this%xmf_h2osfc_col           (:)   = nan
    allocate(this%fact_col                 (begc:endc, -nlevsno+1:nlevgrnd)) ; this%fact_col                 (:,:) = nan
    allocate(this%c_h2osfc_col             (begc:endc))                      ; this%c_h2osfc_col             (:)   = nan

  end subroutine InitAllocate

  !------------------------------------------------------------------------
  subroutine InitHistory(this, bounds)
    !
    ! !DESCRIPTION:
    ! Setup the fields that can be output on history files.
    !
    ! !USES:
    use shr_infnan_mod , only : nan => shr_infnan_nan, assignment(=)
    use histFileMod    , only : hist_addfld1d, hist_addfld2d, no_snow_normal
    !
    ! !ARGUMENTS:
    class(temperature_type) :: this
    type(bounds_type), intent(in) :: bounds  
    !
    ! !LOCAL VARIABLES:
    integer           :: begp, endp
    integer           :: begc, endc
    integer           :: begl, endl
    integer           :: begg, endg
    character(10)     :: active
    character(100)    :: lname
    real(r8), pointer :: data2dptr(:,:), data1dptr(:) ! temp. pointers for slicing larger arrays
    !------------------------------------------------------------------------

    begp = bounds%begp; endp= bounds%endp
    begc = bounds%begc; endc= bounds%endc
    begl = bounds%begl; endl= bounds%endl
    begg = bounds%begg; endg= bounds%endg

    this%t_soisno_col(begc:endc,-nlevsno+1:0) = spval
    data2dptr => this%t_soisno_col(:,-nlevsno+1:0)
    call hist_addfld2d (fname='SNO_T', units='K', type2d='levsno',  &
         avgflag='A', long_name='Snow temperatures', &
         ptr_col=data2dptr, no_snow_behavior=no_snow_normal, default='inactive')

    this%t_ref2m_patch(begp:endp) = spval
    call hist_addfld1d (fname='TSA', units='K',  &
         avgflag='A', long_name='2m air temperature', &
         ptr_patch=this%t_ref2m_patch, default='inactive')

    this%t_ref2m_min_patch(begp:endp) = spval
    call hist_addfld1d (fname='TREFMNAV', units='K',  &
         avgflag='A', long_name='daily minimum of average 2-m temperature', &
         ptr_patch=this%t_ref2m_min_patch, default='inactive')

    this%t_ref2m_max_patch(begp:endp) = spval
    call hist_addfld1d (fname='TREFMXAV', units='K',  &
         avgflag='A', long_name='daily maximum of average 2-m temperature', &
         ptr_patch=this%t_ref2m_max_patch, default='inactive')

    this%t_grnd_col(begc:endc) = spval
    call hist_addfld1d (fname='TG',  units='K',  &
         avgflag='A', long_name='ground temperature', &
         ptr_col=this%t_grnd_col, c2l_scale_type='urbans', default='inactive')

    this%t_soisno_col(begc:endc,:) = spval
    call hist_addfld2d (fname='TSOI',  units='K', type2d='levgrnd', &
         avgflag='A', long_name='soil temperature (vegetated landunits only)', &
         ptr_col=this%t_soisno_col, l2g_scale_type='veg', default='inactive')

    this%t_soi10cm_col(begc:endc) = spval
    call hist_addfld1d (fname='TSOI_10CM',  units='K', &
         avgflag='A', long_name='soil temperature in top 10cm of soil', &
         ptr_col=this%t_soi10cm_col, set_urb=spval, default='inactive')

    active = "active"
    this%t_a10_patch(begp:endp) = spval
    call hist_addfld1d (fname='T10', units='K',  &
         avgflag='A', long_name='10-day running mean of 2-m temperature', &
         ptr_patch=this%t_a10_patch, default='inactive')

    this%heat1_grc(begg:endg) = spval
    call hist_addfld1d (fname='HEAT_CONTENT1',  units='J/m^2',  &
         avgflag='A', long_name='initial gridcell total heat content', &
         ptr_lnd=this%heat1_grc, default='inactive')
    call hist_addfld1d (fname='HEAT_CONTENT1_VEG',  units='J/m^2',  &
         avgflag='A', long_name='initial gridcell total heat content - vegetated landunits only', &
         ptr_lnd=this%heat1_grc, l2g_scale_type='veg', default='inactive')

    this%heat2_grc(begg:endg) = spval
    call hist_addfld1d (fname='HEAT_CONTENT2',  units='J/m^2',  &
         avgflag='A', long_name='post land cover change total heat content', &
         ptr_lnd=this%heat2_grc, default='inactive')  

  end subroutine InitHistory

  !-----------------------------------------------------------------------
  subroutine InitCold(this, bounds)
    !
    ! !DESCRIPTION:
    ! Initialize cold start conditions for module variables
    !
    ! !USES:
    use shr_kind_mod   , only : r8 => shr_kind_r8
    use shr_const_mod  , only : SHR_CONST_TKFRZ
    use clm_varcon     , only : denice, denh2o, sb
    use landunit_varcon, only : istwet, istsoil, istdlak, istice_mec
    use column_varcon  , only : icol_road_imperv, icol_roof, icol_sunwall
    use column_varcon  , only : icol_shadewall, icol_road_perv
    use clm_varctl     , only : iulog
    !
    ! !ARGUMENTS:
    class(temperature_type)        :: this
    type(bounds_type) , intent(in) :: bounds  
    !
    ! !LOCAL VARIABLES:
    integer  :: j,l,c,p ! indices
    integer  :: nlevs   ! number of levels
    real(r8) :: snowbd  ! temporary calculation of snow bulk density (kg/m3)
    real(r8) :: fmelt   ! snowbd/100
    integer  :: lev
    !-----------------------------------------------------------------------

    associate(snl => col%snl) ! Output: [integer (:)    ]  number of snow layers   
      ! Set snow/soil temperature
      ! t_lake only has valid values over non-lake   
      ! t_soisno, t_grnd and t_veg have valid values over all land 

      do c = bounds%begc,bounds%endc
         l = col%landunit(c)

         this%t_soisno_col(c,-nlevsno+1:nlevgrnd) = spval

         ! Snow level temperatures - all land points
         if (snl(c) < 0) then
            do j = snl(c)+1, 0
               this%t_soisno_col(c,j) = 250._r8
            end do
         end if

         ! Below snow temperatures - nonlake points (lake points are set below)
         if (.not. lun%lakpoi(l)) then 

            if (lun%itype(l)==istice_mec) then
               this%t_soisno_col(c,1:nlevgrnd) = 250._r8

            else if (lun%itype(l) == istwet) then
               this%t_soisno_col(c,1:nlevgrnd) = 277._r8

            else if (lun%urbpoi(l)) then
                  if (col%itype(c) == icol_road_perv .or. col%itype(c) == icol_road_imperv) then 
                     this%t_soisno_col(c,1:nlevgrnd) = 274._r8
                  else if (col%itype(c) == icol_sunwall .or. col%itype(c) == icol_shadewall &
                       .or. col%itype(c) == icol_roof) then
                     ! Set sunwall, shadewall, roof to fairly high temperature to avoid initialization
                     ! shock from large heating/air conditioning flux
                     this%t_soisno_col(c,1:nlevurb) = 292._r8
                  end if
            else
               this%t_soisno_col(c,1:nlevgrnd) = 274._r8

            endif
         endif
      end do

      ! Set Ground temperatures

      do c = bounds%begc,bounds%endc
         l = col%landunit(c)

         if (lun%lakpoi(l)) then 
            this%t_grnd_col(c) = 277._r8
         else
            this%t_grnd_col(c) = this%t_soisno_col(c,snl(c)+1)
         end if
      end do

      do c = bounds%begc,bounds%endc
         l = col%landunit(c)
         if (lun%lakpoi(l)) then ! lake
            this%t_lake_col(c,1:nlevlak) = this%t_grnd_col(c)
            this%t_soisno_col(c,1:nlevgrnd) = this%t_grnd_col(c)
         end if
      end do

      ! Set t_h2osfc_col

      this%t_h2osfc_col(bounds%begc:bounds%endc)  = 274._r8

      ! Set t_veg, t_ref2m

      do p = bounds%begp, bounds%endp
         c = patch%column(p)
         l = patch%landunit(p)

         this%t_veg_patch(p)   = 283._r8
         this%t_ref2m_patch(p) = 283._r8

      end do

    end associate

    do l = bounds%begl, bounds%endl 
       if (lun%urbpoi(l)) then
          this%taf_lun(l) = 283._r8
       end if
    end do

  end subroutine InitCold

  !------------------------------------------------------------------------
  subroutine Restart(this, bounds, ncid, flag)
    ! 
    ! !DESCRIPTION:
    ! Read/Write module information to/from restart file.
    !
    ! !USES:
    use shr_log_mod     , only : errMsg => shr_log_errMsg
    use spmdMod         , only : masterproc
    use abortutils      , only : endrun
    use ncdio_pio       , only : file_desc_t, ncd_double
    use restUtilMod
    !
    ! !ARGUMENTS:
    class(temperature_type) :: this
    type(bounds_type), intent(in)    :: bounds 
    type(file_desc_t), intent(inout) :: ncid   
    character(len=*) , intent(in)    :: flag   
    !
    ! !LOCAL VARIABLES:
    integer :: j,c       ! indices
    logical :: readvar   ! determine if variable is on initial file
    !-----------------------------------------------------------------------

    call restartvar(ncid=ncid, flag=flag, varname='T_SOISNO', xtype=ncd_double,   &
         dim1name='column', dim2name='levtot', switchdim=.true., &
         long_name='soil-snow temperature', units='K', &
         interpinic_flag='interp', readvar=readvar, data=this%t_soisno_col)

    call restartvar(ncid=ncid, flag=flag, varname='T_VEG', xtype=ncd_double,  &
         dim1name='pft', &
         long_name='vegetation temperature', units='K', &
         interpinic_flag='interp', readvar=readvar, data=this%t_veg_patch)

    call restartvar(ncid=ncid, flag=flag, varname='TH2OSFC', xtype=ncd_double,  &
         dim1name='column', &
         long_name='surface water temperature', units='K', &
         interpinic_flag='interp', readvar=readvar, data=this%t_h2osfc_col)
    if (flag=='read' .and. .not. readvar) then
       this%t_h2osfc_col(bounds%begc:bounds%endc) = 274.0_r8
    end if

    call restartvar(ncid=ncid, flag=flag, varname='T_LAKE', xtype=ncd_double,  &
         dim1name='column', dim2name='levlak', switchdim=.true., &
         long_name='lake temperature', units='K', &
         interpinic_flag='interp', readvar=readvar, data=this%t_lake_col)

    call restartvar(ncid=ncid, flag=flag, varname='T_GRND', xtype=ncd_double,  &
         dim1name='column', &
         long_name='ground temperature', units='K', &
         interpinic_flag='interp', readvar=readvar, data=this%t_grnd_col)

    call restartvar(ncid=ncid, flag=flag, varname='T_REF2M', xtype=ncd_double,  &
         dim1name='pft', &
         long_name='2m height surface air temperature', units='K', &
         interpinic_flag='interp', readvar=readvar, data=this%t_ref2m_patch)
    if (flag=='read' .and. .not. readvar) call endrun(msg=errMsg(sourcefile, __LINE__))

    call restartvar(ncid=ncid, flag=flag, varname='T_REF2M_MIN', xtype=ncd_double,  &
         dim1name='pft', &
         long_name='daily minimum of average 2 m height surface air temperature (K)', units='K', &
         interpinic_flag='interp', readvar=readvar, data=this%t_ref2m_min_patch)

    call restartvar(ncid=ncid, flag=flag, varname='T_REF2M_MAX', xtype=ncd_double,  &
         dim1name='pft', &
         long_name='daily maximum of average 2 m height surface air temperature (K)', units='K', &
         interpinic_flag='interp', readvar=readvar, data=this%t_ref2m_max_patch)

    call restartvar(ncid=ncid, flag=flag, varname='T_REF2M_MIN_INST', xtype=ncd_double,  &
         dim1name='pft', &
         long_name='instantaneous daily min of average 2 m height surface air temp (K)', units='K', &
         interpinic_flag='interp', readvar=readvar, data=this%t_ref2m_min_inst_patch)

    call restartvar(ncid=ncid, flag=flag, varname='T_REF2M_MAX_INST', xtype=ncd_double,  &
         dim1name='pft', &
         long_name='instantaneous daily max of average 2 m height surface air temp (K)', units='K', &
         interpinic_flag='interp', readvar=readvar, data=this%t_ref2m_max_inst_patch)

  end subroutine Restart

end module TemperatureType
