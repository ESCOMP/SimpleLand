module WaterstateType

#include "shr_assert.h"

  !------------------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Module variables for hydrology
  !
  ! !USES:
  use shr_kind_mod   , only : r8 => shr_kind_r8
  use shr_log_mod    , only : errMsg => shr_log_errMsg
  use decompMod      , only : bounds_type
  use clm_varpar     , only : nlevgrnd, nlevurb, nlevsno   
  use clm_varcon     , only : spval
  use LandunitType   , only : lun                
  use ColumnType     , only : col                
  !
  implicit none
  save
  private
  !
  ! !PUBLIC TYPES:
  type, public :: waterstate_type

     real(r8), pointer :: snow_depth_col         (:)   ! col snow height of snow covered area (m)

     real(r8), pointer :: h2osno_col             (:)   ! col snow water (mm H2O)
     real(r8), pointer :: h2osoi_liq_col         (:,:) ! col liquid water (kg/m2) (new) (-nlevsno+1:nlevgrnd)    
     real(r8), pointer :: h2osoi_ice_col         (:,:) ! col ice lens (kg/m2) (new) (-nlevsno+1:nlevgrnd)    
     real(r8), pointer :: h2osoi_liq_tot_col     (:)   ! vertically summed col liquid water (kg/m2) (new) (-nlevsno+1:nlevgrnd)    
     real(r8), pointer :: h2osoi_ice_tot_col     (:)   ! vertically summed col ice lens (kg/m2) (new) (-nlevsno+1:nlevgrnd)    
     real(r8), pointer :: h2osoi_vol_col         (:,:) ! col volumetric soil water (0<=h2osoi_vol<=watsat) [m3/m3]  (nlevgrnd)
     real(r8), pointer :: h2osoi_liqvol_col      (:,:) ! col volumetric liquid water content (v/v)
     real(r8), pointer :: tws_grc                (:)   ! grc total water storage (mm H2O)

     real(r8), pointer :: q_ref2m_patch          (:)   ! patch 2 m height surface specific humidity (kg/kg)

     ! Balance Checks
     real(r8), pointer :: endwb_col              (:)   ! water mass end of the time step

   contains

     procedure          :: Init         
     procedure          :: Restart      
     procedure, private :: InitAllocate 
     procedure, private :: InitHistory  
     procedure, private :: InitCold     

  end type waterstate_type

  character(len=*), parameter, private :: sourcefile = &
       __FILE__
 !------------------------------------------------------------------------

contains

  !------------------------------------------------------------------------
  subroutine Init(this, bounds, &
       h2osno_input_col, snow_depth_input_col, watsat_col, t_soisno_col)

    class(waterstate_type)            :: this
    type(bounds_type) , intent(in)    :: bounds  
    real(r8)          , intent(inout) :: h2osno_input_col(bounds%begc:)
    real(r8)          , intent(inout) :: snow_depth_input_col(bounds%begc:)
    real(r8)          , intent(inout) :: watsat_col(bounds%begc:, 1:)          ! volumetric soil water at saturation (porosity)
    real(r8)          , intent(inout) :: t_soisno_col(bounds%begc:, -nlevsno+1:) ! col soil temperature (Kelvin)

#ifdef __PGI
# if __PGIC__ == 14 && __PGIC_MINOR__ == 7
    ! COMPILER_BUG(bja, 2015-04, pgi 14.7-?) occurs at: call this%InitCold(...)
    ! PGF90-F-0000-Internal compiler error. normalize_forall_array: non-conformable
    ! not sure why this fixes things....
    real(r8), allocatable :: workaround_for_pgi_internal_compiler_error(:)
# endif
#endif

    call this%InitAllocate(bounds) 

    call this%InitHistory(bounds)

    call this%InitCold(bounds, &
       h2osno_input_col, snow_depth_input_col, watsat_col, t_soisno_col)

  end subroutine Init

  !------------------------------------------------------------------------
  subroutine InitAllocate(this, bounds)
    !
    ! !DESCRIPTION:
    ! Initialize module data structure
    !
    ! !USES:
    use shr_infnan_mod , only : nan => shr_infnan_nan, assignment(=)
    !
    ! !ARGUMENTS:
    class(waterstate_type) :: this
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

    allocate(this%snow_depth_col         (begc:endc))                     ; this%snow_depth_col         (:)   = nan
    allocate(this%h2osno_col             (begc:endc))                     ; this%h2osno_col             (:)   = nan   
    allocate(this%h2osoi_vol_col         (begc:endc, 1:nlevgrnd))         ; this%h2osoi_vol_col         (:,:) = nan
    allocate(this%h2osoi_liqvol_col      (begc:endc,-nlevsno+1:nlevgrnd)) ; this%h2osoi_liqvol_col      (:,:) = nan
    allocate(this%h2osoi_ice_col         (begc:endc,-nlevsno+1:nlevgrnd)) ; this%h2osoi_ice_col         (:,:) = nan
    allocate(this%h2osoi_liq_col         (begc:endc,-nlevsno+1:nlevgrnd)) ; this%h2osoi_liq_col         (:,:) = nan
    allocate(this%h2osoi_ice_tot_col     (begc:endc))                     ; this%h2osoi_ice_tot_col     (:)   = nan
    allocate(this%h2osoi_liq_tot_col     (begc:endc))                     ; this%h2osoi_liq_tot_col     (:)   = nan
    allocate(this%tws_grc                (begg:endg))                     ; this%tws_grc                (:)   = nan

    allocate(this%q_ref2m_patch          (begp:endp))                     ; this%q_ref2m_patch          (:)   = nan

    allocate(this%endwb_col              (begc:endc))                     ; this%endwb_col              (:)   = nan
  end subroutine InitAllocate

  !------------------------------------------------------------------------
  subroutine InitHistory(this, bounds)
    !
    ! !DESCRIPTION:
    ! Initialize module data structure
    !
    ! !USES:
    use shr_infnan_mod , only : nan => shr_infnan_nan, assignment(=)
    use clm_varpar     , only : nlevsoi
    use histFileMod    , only : hist_addfld1d, hist_addfld2d
    !
    ! !ARGUMENTS:
    class(waterstate_type) :: this
    type(bounds_type), intent(in) :: bounds  
    !
    ! !LOCAL VARIABLES:
    integer           :: begp, endp
    integer           :: begc, endc
    integer           :: begg, endg
    real(r8), pointer :: data2dptr(:,:)  ! temp. pointers for slicing larger arrays
    !------------------------------------------------------------------------

    begp = bounds%begp; endp= bounds%endp
    begc = bounds%begc; endc= bounds%endc
    begg = bounds%begg; endg= bounds%endg

    ! h2osno also includes snow that is part of the soil column (an 
    ! initial snow layer is only created if h2osno > 10mm). 

    data2dptr => this%h2osoi_vol_col(begc:endc,1:nlevsoi)
    call hist_addfld2d (fname='H2OSOI',  units='mm3/mm3', type2d='levsoi', &
         avgflag='A', long_name='volumetric soil water (vegetated landunits only)', &
         ptr_col=this%h2osoi_vol_col, l2g_scale_type='veg', default='inactive')

    call hist_addfld1d (fname='H2OSNO',  units='mm',  &
         avgflag='A', long_name='snow depth (liquid water)', &
         ptr_col=this%h2osno_col, c2l_scale_type='urbanf', default='inactive')

    this%tws_grc(begg:endg) = spval
    call hist_addfld1d (fname='TWS',  units='mm',  &
         avgflag='A', long_name='total water storage', &
         ptr_lnd=this%tws_grc, default='inactive')

    ! Humidity

    this%q_ref2m_patch(begp:endp) = spval
    call hist_addfld1d (fname='Q2M', units='kg/kg',  &
         avgflag='A', long_name='2m specific humidity', &
         ptr_patch=this%q_ref2m_patch, default='inactive')

  end subroutine InitHistory

  !-----------------------------------------------------------------------
  subroutine InitCold(this, bounds, &
       h2osno_input_col, snow_depth_input_col, watsat_col, t_soisno_col)
    !
    ! !DESCRIPTION:
    ! Initialize time constant variables and cold start conditions 
    !
    ! !USES:
    use shr_log_mod     , only : errMsg => shr_log_errMsg
    use shr_kind_mod    , only : r8 => shr_kind_r8
    use shr_const_mod   , only : SHR_CONST_TKFRZ
    use clm_varpar      , only : nlevsoi, nlevgrnd, nlevsno, nlevurb
    use landunit_varcon , only : istwet, istsoil, istcrop, istice_mec  
    use column_varcon   , only : icol_road_perv
    use column_varcon   , only : icol_road_imperv
    use clm_varcon      , only : denice, denh2o, spval, bdsno 
    use clm_varcon      , only : tfrz, spval
    use spmdMod         , only : masterproc
    use abortutils      , only : endrun
    use fileutils       , only : getfil
    use ncdio_pio       , only : file_desc_t, ncd_io
    !
    ! !ARGUMENTS:
    class(waterstate_type)                :: this
    type(bounds_type)     , intent(in)    :: bounds
    real(r8)              , intent(in)    :: h2osno_input_col(bounds%begc:)
    real(r8)              , intent(in)    :: snow_depth_input_col(bounds%begc:)
    real(r8)              , intent(in)    :: watsat_col(bounds%begc:, 1:)          ! volumetric soil water at saturation (porosity)
    real(r8)              , intent(in)    :: t_soisno_col(bounds%begc:, -nlevsno+1:) ! col soil temperature (Kelvin)
    !
    ! !LOCAL VARIABLES:
    integer            :: p,c,j,l,g,lev,nlevs 
    real(r8)           :: d
    type(file_desc_t)  :: ncid        
    !-----------------------------------------------------------------------

    SHR_ASSERT_ALL((ubound(h2osno_input_col)     == (/bounds%endc/))          , errMsg(sourcefile, __LINE__))
    SHR_ASSERT_ALL((ubound(snow_depth_input_col) == (/bounds%endc/))          , errMsg(sourcefile, __LINE__))
    SHR_ASSERT_ALL((ubound(watsat_col)           == (/bounds%endc,nlevgrnd/)) , errMsg(sourcefile, __LINE__))
    SHR_ASSERT_ALL((ubound(t_soisno_col)         == (/bounds%endc,nlevgrnd/)) , errMsg(sourcefile, __LINE__))

    ! The first three arrays are initialized from the input argument
    do c = bounds%begc,bounds%endc
       this%h2osno_col(c)             = h2osno_input_col(c) 
       this%snow_depth_col(c)         = snow_depth_input_col(c)
    end do

    associate(snl => col%snl) 

      !--------------------------------------------
      ! Set soil water
      !--------------------------------------------

      ! volumetric water is set first and liquid content and ice lens are obtained
      ! NOTE: h2osoi_vol, h2osoi_liq and h2osoi_ice only have valid values over soil
      ! and urban pervious road (other urban columns have zero soil water)

      this%h2osoi_vol_col(bounds%begc:bounds%endc,         1:) = spval
      this%h2osoi_liq_col(bounds%begc:bounds%endc,-nlevsno+1:) = spval
      this%h2osoi_ice_col(bounds%begc:bounds%endc,-nlevsno+1:) = spval
      do c = bounds%begc,bounds%endc
         l = col%landunit(c)
         if (.not. lun%lakpoi(l)) then  !not lake

            ! volumetric water
            if (lun%itype(l) == istsoil .or. lun%itype(l) == istcrop) then
               nlevs = nlevgrnd
               do j = 1, nlevs
                  if (j > nlevsoi) then
                     this%h2osoi_vol_col(c,j) = 0.0_r8
                  else
                     this%h2osoi_vol_col(c,j) = 0.15_r8
                  endif
               end do
            else if (lun%urbpoi(l)) then
               if (col%itype(c) == icol_road_perv) then
                  nlevs = nlevgrnd
                  do j = 1, nlevs
                     if (j <= nlevsoi) then
                        this%h2osoi_vol_col(c,j) = 0.3_r8
                     else
                        this%h2osoi_vol_col(c,j) = 0.0_r8
                     end if
                  end do
               else if (col%itype(c) == icol_road_imperv) then
                  nlevs = nlevgrnd
                  do j = 1, nlevs
                     this%h2osoi_vol_col(c,j) = 0.0_r8
                  end do
               else
                  nlevs = nlevurb
                  do j = 1, nlevs
                     this%h2osoi_vol_col(c,j) = 0.0_r8
                  end do
               end if
            else if (lun%itype(l) == istwet) then
               nlevs = nlevgrnd
               do j = 1, nlevs
                  if (j > nlevsoi) then
                     this%h2osoi_vol_col(c,j) = 0.0_r8
                  else
                     this%h2osoi_vol_col(c,j) = 1.0_r8
                  endif
               end do
            else if (lun%itype(l) == istice_mec) then
               nlevs = nlevgrnd 
               do j = 1, nlevs
                  this%h2osoi_vol_col(c,j) = 1.0_r8
               end do
            endif
            do j = 1, nlevs
               this%h2osoi_vol_col(c,j) = min(this%h2osoi_vol_col(c,j), watsat_col(c,j))
               if (t_soisno_col(c,j) <= SHR_CONST_TKFRZ) then
                  this%h2osoi_ice_col(c,j) = col%dz(c,j)*denice*this%h2osoi_vol_col(c,j)
                  this%h2osoi_liq_col(c,j) = 0._r8
               else
                  this%h2osoi_ice_col(c,j) = 0._r8
                  this%h2osoi_liq_col(c,j) = col%dz(c,j)*denh2o*this%h2osoi_vol_col(c,j)
               endif
            end do
            do j = -nlevsno+1, 0
               if (j > snl(c)) then
                  this%h2osoi_ice_col(c,j) = col%dz(c,j)*250._r8
                  this%h2osoi_liq_col(c,j) = 0._r8
               end if
            end do
         end if
      end do


      !--------------------------------------------
      ! Set Lake water
      !--------------------------------------------

      do c = bounds%begc, bounds%endc
         l = col%landunit(c)

         if (lun%lakpoi(l)) then
            do j = -nlevsno+1, 0
               if (j > snl(c)) then
                  this%h2osoi_ice_col(c,j) = col%dz(c,j)*bdsno
                  this%h2osoi_liq_col(c,j) = 0._r8
               end if
            end do
            do j = 1,nlevgrnd
               if (j <= nlevsoi) then ! soil
                  this%h2osoi_vol_col(c,j) = watsat_col(c,j)
                  this%h2osoi_liq_col(c,j) = spval
                  this%h2osoi_ice_col(c,j) = spval
               else                  ! bedrock
                  this%h2osoi_vol_col(c,j) = 0._r8
               end if
            end do
         end if
      end do

      !--------------------------------------------
      ! For frozen layers !TODO - does the following make sense ???? it seems to overwrite everything
      !--------------------------------------------

      do c = bounds%begc, bounds%endc
         do j = 1,nlevgrnd
            if (this%h2osoi_vol_col(c,j) /= spval) then
               if (t_soisno_col(c,j) <= tfrz) then
                  this%h2osoi_ice_col(c,j) = col%dz(c,j)*denice*this%h2osoi_vol_col(c,j)
                  this%h2osoi_liq_col(c,j) = 0._r8
               else
                  this%h2osoi_ice_col(c,j) = 0._r8
                  this%h2osoi_liq_col(c,j) = col%dz(c,j)*denh2o*this%h2osoi_vol_col(c,j)
               endif
            end if
         end do
      end do

    end associate

  end subroutine InitCold

  !------------------------------------------------------------------------
  subroutine Restart(this, bounds, ncid, flag, &
       watsat_col)
    ! 
    ! !DESCRIPTION:
    ! Read/Write module information to/from restart file.
    !
    ! !USES:
    use spmdMod          , only : masterproc
    use clm_varcon       , only : denice, denh2o, pondmx, watmin, spval, nameg
    use landunit_varcon  , only : istcrop, istdlak, istsoil  
    use column_varcon    , only : icol_roof, icol_sunwall, icol_shadewall
    use clm_time_manager , only : is_first_step
    use clm_varctl       , only : bound_h2osoi
    use ncdio_pio        , only : file_desc_t, ncd_io, ncd_double
    use restUtilMod
    !
    ! !ARGUMENTS:
    class(waterstate_type) :: this
    type(bounds_type), intent(in)    :: bounds 
    type(file_desc_t), intent(inout) :: ncid   ! netcdf id
    character(len=*) , intent(in)    :: flag   ! 'read' or 'write'
    real(r8)         , intent(in)    :: watsat_col (bounds%begc:, 1:)  ! volumetric soil water at saturation (porosity)
    !
    ! !LOCAL VARIABLES:
    integer  :: c,l,j,nlevs
    logical  :: readvar
    real(r8) :: maxwatsat    ! maximum porosity    
    real(r8) :: excess       ! excess volumetric soil water
    real(r8) :: totwat       ! total soil water (mm)
    !------------------------------------------------------------------------

    SHR_ASSERT_ALL((ubound(watsat_col) == (/bounds%endc,nlevgrnd/)) , errMsg(sourcefile, __LINE__))

    call restartvar(ncid=ncid, flag=flag, varname='H2OSNO', xtype=ncd_double,  &
         dim1name='column', &
         long_name='snow water', units='kg/m2', &
         interpinic_flag='interp', readvar=readvar, data=this%h2osno_col)

    call restartvar(ncid=ncid, flag=flag, varname='H2OSOI_LIQ', xtype=ncd_double,  &
         dim1name='column', dim2name='levtot', switchdim=.true., &
         long_name='liquid water', units='kg/m2', &
         interpinic_flag='interp', readvar=readvar, data=this%h2osoi_liq_col)

    call restartvar(ncid=ncid, flag=flag, varname='H2OSOI_ICE', xtype=ncd_double,   &
         dim1name='column', dim2name='levtot', switchdim=.true., &
         long_name='ice lens', units='kg/m2', &
         interpinic_flag='interp', readvar=readvar, data=this%h2osoi_ice_col)
         
    ! Determine volumetric soil water (for read only)
    if (flag == 'read' ) then
       do c = bounds%begc, bounds%endc
          l = col%landunit(c)
          if ( col%itype(c) == icol_sunwall   .or. &
               col%itype(c) == icol_shadewall .or. &
               col%itype(c) == icol_roof )then
             nlevs = nlevurb
          else
             nlevs = nlevgrnd
          end if
          if ( lun%itype(l) /= istdlak ) then ! This calculation is now done for lakes in initLake.
             do j = 1,nlevs
                this%h2osoi_vol_col(c,j) = this%h2osoi_liq_col(c,j)/(col%dz(c,j)*denh2o) &
                                         + this%h2osoi_ice_col(c,j)/(col%dz(c,j)*denice)
             end do
          end if
       end do
    end if

    ! If initial run -- ensure that water is properly bounded (read only)
    if (flag == 'read' ) then
       if ( is_first_step() .and. bound_h2osoi) then
          do c = bounds%begc, bounds%endc
             l = col%landunit(c)
             if ( col%itype(c) == icol_sunwall .or. col%itype(c) == icol_shadewall .or. &
                  col%itype(c) == icol_roof )then
                nlevs = nlevurb
             else
                nlevs = nlevgrnd
             end if
             do j = 1,nlevs
                l = col%landunit(c)
                if (lun%itype(l) == istsoil .or. lun%itype(l) == istcrop) then
                   this%h2osoi_liq_col(c,j) = max(0._r8,this%h2osoi_liq_col(c,j))
                   this%h2osoi_ice_col(c,j) = max(0._r8,this%h2osoi_ice_col(c,j))
                   this%h2osoi_vol_col(c,j) = this%h2osoi_liq_col(c,j)/(col%dz(c,j)*denh2o) &
                                       + this%h2osoi_ice_col(c,j)/(col%dz(c,j)*denice)
                   if (j == 1) then
                      maxwatsat = (watsat_col(c,j)*col%dz(c,j)*1000.0_r8 + pondmx) / (col%dz(c,j)*1000.0_r8)
                   else
                      maxwatsat =  watsat_col(c,j)
                   end if
                   if (this%h2osoi_vol_col(c,j) > maxwatsat) then 
                      excess = (this%h2osoi_vol_col(c,j) - maxwatsat)*col%dz(c,j)*1000.0_r8
                      totwat = this%h2osoi_liq_col(c,j) + this%h2osoi_ice_col(c,j)
                      this%h2osoi_liq_col(c,j) = this%h2osoi_liq_col(c,j) - &
                                           (this%h2osoi_liq_col(c,j)/totwat) * excess
                      this%h2osoi_ice_col(c,j) = this%h2osoi_ice_col(c,j) - &
                                           (this%h2osoi_ice_col(c,j)/totwat) * excess
                   end if
                   this%h2osoi_liq_col(c,j) = max(watmin,this%h2osoi_liq_col(c,j))
                   this%h2osoi_ice_col(c,j) = max(watmin,this%h2osoi_ice_col(c,j))
                   this%h2osoi_vol_col(c,j) = this%h2osoi_liq_col(c,j)/(col%dz(c,j)*denh2o) &
                                             + this%h2osoi_ice_col(c,j)/(col%dz(c,j)*denice)
                end if
             end do
          end do
       end if

    endif   ! end if if-read flag

    call restartvar(ncid=ncid, flag=flag, varname='SNOW_DEPTH', xtype=ncd_double,  & 
         dim1name='column', &
         long_name='snow depth', units='m', &
         interpinic_flag='interp', readvar=readvar, data=this%snow_depth_col) 

  end subroutine Restart

end module WaterstateType
