module UrbanTimeVarType

  !------------------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Urban Time Varying Data
  !
  ! !USES:
  use shr_kind_mod    , only : r8 => shr_kind_r8, CL => shr_kind_CL
  use shr_log_mod     , only : errMsg => shr_log_errMsg
  use abortutils      , only : endrun
  use decompMod       , only : bounds_type
  use clm_varctl      , only : iulog
  use landunit_varcon , only : isturb_MIN, isturb_MAX
  use clm_varcon      , only : spval
  use LandunitType    , only : lun                
  use GridcellType    , only : grc
  use mct_mod
  use shr_strdata_mod , only : shr_strdata_type
  !
  implicit none
  save
  private
  !
  !

  ! !PUBLIC TYPE
  type, public :: urbantv_type

     real(r8), public, pointer     :: t_building_max(:)    ! lun maximum internal building air temperature (K)
     type(shr_strdata_type)        :: sdat_urbantv         ! urban time varying input data stream
   contains

     ! !PUBLIC MEMBER FUNCTIONS:
     procedure, public :: Init              ! Allocate and initialize urbantv
     procedure, public :: urbantv_interp    ! Interpolate urban time varying stream
     
  end type urbantv_type

  !----------------------------------------------------------------------- 
  character(15), private :: stream_var_name(isturb_MIN:isturb_MAX)

  character(len=*), parameter, private :: sourcefile = &
       __FILE__

contains

  !-----------------------------------------------------------------------
  subroutine Init(this, bounds, NLFilename)
    !
    ! Allocate module variables and data structures
    !
    ! !USES:
    use shr_infnan_mod  , only : nan => shr_infnan_nan, assignment(=)
    use histFileMod     , only : hist_addfld1d
    !
    ! !ARGUMENTS:
    class(urbantv_type) :: this
    type(bounds_type) , intent(in) :: bounds  
    character(len=*)  , intent(in) :: NLFilename   ! Namelist filename
    !
    ! !LOCAL VARIABLES:
    integer		:: begl, endl
    !---------------------------------------------------------------------

    begl = bounds%begl; endl = bounds%endl

    ! Allocate urbantv data structure

    allocate(this%t_building_max      (begl:endl))          ; this%t_building_max      (:)   = nan

    call this%urbantv_interp(bounds)

    ! Add history fields
    call hist_addfld1d (fname='TBUILD_MAX', units='K',      &
          avgflag='A', long_name='prescribed maximum interior building temperature',   &
          ptr_lunit=this%t_building_max, default='inactive', set_nourb=spval, &
          l2g_scale_type='unity')


  end subroutine Init

  !-----------------------------------------------------------------------
  subroutine urbantv_interp(this, bounds)
  !
  ! !DESCRIPTION:
  ! Interpolate data stream information for urban time varying data.
  !
  ! !USES:
  use clm_time_manager, only : get_curr_date
  use spmdMod         , only : mpicom
  use shr_strdata_mod , only : shr_strdata_advance
  use clm_instur      , only : urban_valid
  !
  ! !ARGUMENTS:
  class(urbantv_type)           :: this
  type(bounds_type), intent(in) :: bounds
  !
  ! !LOCAL VARIABLES:
  logical :: found
  integer :: l, glun, ig, g, ip
  integer :: year    ! year (0, ...) for nstep+1
  integer :: mon     ! month (1, ..., 12) for nstep+1
  integer :: day     ! day of month (1, ..., 31) for nstep+1
  integer :: sec     ! seconds into current date for nstep+1
  integer :: mcdate  ! Current model date (yyyymmdd)
  integer :: lindx   ! landunit index
  integer :: gindx   ! gridcell index
  !-----------------------------------------------------------------------

   call get_curr_date(year, mon, day, sec)
   mcdate = year*10000 + mon*100 + day

   call shr_strdata_advance(this%sdat_urbantv, mcdate, sec, mpicom, 'urbantvdyn')

   do l = bounds%begl,bounds%endl
      if (lun%urbpoi(l)) then
         glun  = lun%gridcell(l)
         ip = mct_aVect_indexRA(this%sdat_urbantv%avs(1),trim(stream_var_name(lun%itype(l))))
         !
         ! Determine vector index corresponding to glun
         !
         ig = 0
         do g = bounds%begg,bounds%endg
            ig = ig+1
            if (g == glun) exit
         end do

         this%t_building_max(l) = this%sdat_urbantv%avs(1)%rAttr(ip,ig)
      else
         this%t_building_max(l) = spval
      end if
   end do

   found = .false.
   do l = bounds%begl,bounds%endl
      if (lun%urbpoi(l)) then
         glun  = lun%gridcell(l)
         !
         ! Determine vector index corresponding to glun
         !
         ig = 0
         do g = bounds%begg,bounds%endg
            ig = ig+1
            if (g == glun) exit
         end do

         if ( .not. urban_valid(g) .or. (this%t_building_max(l) <= 0._r8)) then
            found = .true.
            gindx = g
            lindx = l
            exit
         end if
      end if
   end do
   if ( found ) then
      write(iulog,*)'ERROR: no valid urban data for g= ',gindx
      write(iulog,*)'landunit type:   ',lun%itype(l)
      write(iulog,*)'urban_valid:     ',urban_valid(gindx)
      write(iulog,*)'t_building_max:  ',this%t_building_max(lindx)
      call endrun(msg=errmsg(sourcefile, __LINE__))
   end if


  end subroutine urbantv_interp

  !-----------------------------------------------------------------------

end module UrbanTimeVarType
