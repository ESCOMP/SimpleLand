module TemperatureType

#include "shr_assert.h"

  !------------------------------------------------------------------------------
  use shr_kind_mod    , only : r8 => shr_kind_r8
  use shr_log_mod     , only : errMsg => shr_log_errMsg
  use decompMod       , only : bounds_type
  use clm_varpar      , only : nlevsno, nlevgrnd, nlevurb
  use clm_varcon      , only : spval
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
     real(r8), pointer :: t_soisno_col (:,:) ! col soil temperature (Kelvin)  (-nlevsno+1:nlevgrnd) 
     real(r8), pointer :: t_grnd_col   (:)   ! col ground temperature (Kelvin)
     real(r8), pointer :: t_ref2m_patch(:)   ! patch 2 m height surface air temperature (Kelvin)

   contains

     procedure, public  :: Init
     procedure, public  :: Restart      
     procedure, private :: InitAllocate 
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
    allocate(this%t_soisno_col (begc:endc,-nlevsno+1:nlevgrnd))  ; this%t_soisno_col             (:,:) = nan
    allocate(this%t_grnd_col   (begc:endc))                      ; this%t_grnd_col               (:)   = nan
    allocate(this%t_ref2m_patch(begp:endp))                      ; this%t_ref2m_patch            (:)   = nan

  end subroutine InitAllocate

  !-----------------------------------------------------------------------
  subroutine InitCold(this, bounds)
    !
    ! !DESCRIPTION:
    ! Initialize cold start conditions for module variables
    !
    ! !USES:
    use shr_kind_mod   , only : r8 => shr_kind_r8
    use landunit_varcon, only : istwet, istice_mec
    use column_varcon  , only : icol_road_imperv, icol_roof, icol_sunwall
    use column_varcon  , only : icol_shadewall, icol_road_perv
    !
    ! !ARGUMENTS:
    class(temperature_type)        :: this
    type(bounds_type) , intent(in) :: bounds  
    !
    ! !LOCAL VARIABLES:
    integer  :: j,l,c,p ! indices
    !-----------------------------------------------------------------------

    associate(snl => col%snl) ! Output: [integer (:)    ]  number of snow layers   
      ! Set snow/soil temperature
      ! t_soisno, t_grnd have valid values over all land 

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
            this%t_soisno_col(c,1:nlevgrnd) = this%t_grnd_col(c)
         end if
      end do

      ! Set t_ref2m

      do p = bounds%begp, bounds%endp
         c = patch%column(p)
         l = patch%landunit(p)
         this%t_ref2m_patch(p) = 283._r8
      end do

    end associate

  end subroutine InitCold

  !------------------------------------------------------------------------
  subroutine Restart(this, bounds, ncid, flag)
    ! 
    ! !DESCRIPTION:
    ! Read/Write module information to/from restart file.
    !
    ! !USES:
    use shr_log_mod     , only : errMsg => shr_log_errMsg
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

    call restartvar(ncid=ncid, flag=flag, varname='T_GRND', xtype=ncd_double,  &
         dim1name='column', &
         long_name='ground temperature', units='K', &
         interpinic_flag='interp', readvar=readvar, data=this%t_grnd_col)

    call restartvar(ncid=ncid, flag=flag, varname='T_REF2M', xtype=ncd_double,  &
         dim1name='pft', &
         long_name='2m height surface air temperature', units='K', &
         interpinic_flag='interp', readvar=readvar, data=this%t_ref2m_patch)
    if (flag=='read' .and. .not. readvar) call endrun(msg=errMsg(sourcefile, __LINE__))

  end subroutine Restart

end module TemperatureType
