module lnd2glcMod

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Handle arrays used for exchanging data from land model to glc
  ! For now glc datais send and received on the lnd grid and decomposition.
  !
  ! The fields sent from the lnd component to the glc component via
  !  the coupler are labeled 's2x', or sno to coupler.
  ! The fields received by the lnd component from the glc component
  !  via the coupler are labeled 'x2s', or coupler to sno.
  ! 'Sno' is a misnomer in that the exchanged data are related to
  !  the ice beneath the snow, not the snow itself.  But by CESM convention,
  ! 'ice' refers to sea ice, not land ice.
  !
  ! !USES:
  use shr_kind_mod    , only : r8 => shr_kind_r8
  use shr_infnan_mod  , only : nan => shr_infnan_nan, assignment(=)
  use shr_log_mod     , only : errMsg => shr_log_errMsg
  use decompMod       , only : get_proc_bounds, bounds_type
  use domainMod       , only : ldomain
  use clm_varctl      , only : iulog
  use clm_varcon      , only : spval, tfrz, namec
  use column_varcon   , only : col_itype_to_icemec_class
  use landunit_varcon , only : istice_mec, istsoil
  use abortutils      , only : endrun
  use TemperatureType , only : temperature_type
  use LandunitType    , only : lun                
  use ColumnType      , only : col
  use TopoMod         , only : topo_type
  !
  ! !PUBLIC TYPES:
  implicit none
  private
  save

  ! land -> glc variables structure
  type, public :: lnd2glc_type
     real(r8), pointer :: tsrf_grc(:,:) => null()
     real(r8), pointer :: topo_grc(:,:) => null()
     real(r8), pointer :: qice_grc(:,:) => null()

   contains

     procedure, public  :: Init
     procedure, public  :: update_lnd2glc
     procedure, private :: InitAllocate

  end type lnd2glc_type

  ! !PUBLIC MEMBER FUNCTIONS:
  
  ! The following is public simply to support unit testing, and should not generally be
  ! called from outside this module.
  !
  ! Note that it is not a type-bound procedure, because it doesn't actually involve the
  ! lnd2glc_type. This suggests that perhaps it belongs in some other module.

  character(len=*), parameter, private :: sourcefile = &
       __FILE__
  !------------------------------------------------------------------------

contains

  !------------------------------------------------------------------------
  subroutine Init(this, bounds)

    class(lnd2glc_type) :: this
    type(bounds_type), intent(in) :: bounds  

    call this%InitAllocate(bounds)
    
  end subroutine Init

  !------------------------------------------------------------------------
  subroutine InitAllocate(this, bounds)
    !
    ! !DESCRIPTION:
    ! Initialize land variables required by glc
    !
    ! !USES:
    use clm_varcon , only : spval
    use histFileMod, only : hist_addfld1d
    !
    ! !ARGUMENTS:
    class(lnd2glc_type) :: this
    type(bounds_type), intent(in) :: bounds  
    !
    ! !LOCAL VARIABLES:
    integer :: begg,endg 
    !------------------------------------------------------------------------

    begg = bounds%begg; endg = bounds%endg

    allocate(this%tsrf_grc(begg:endg,0:10)) ; this%tsrf_grc(:,:)=0.0_r8
    allocate(this%topo_grc(begg:endg,0:10)) ; this%topo_grc(:,:)=0.0_r8
    allocate(this%qice_grc(begg:endg,0:10)) ; this%qice_grc(:,:)=0.0_r8

  end subroutine InitAllocate

  !------------------------------------------------------------------------------
  subroutine update_lnd2glc(this, bounds, num_do_smb_c, filter_do_smb_c, &
       temperature_inst, topo_inst)
    !
    ! !DESCRIPTION:
    ! Assign values to lnd2glc+
    !
    ! !ARGUMENTS:
    class(lnd2glc_type)    , intent(inout) :: this
    type(bounds_type)      , intent(in)    :: bounds  
    integer                , intent(in)    :: num_do_smb_c       ! number of columns in filter_do_smb_c
    integer                , intent(in)    :: filter_do_smb_c(:) ! column filter: columns where smb calculations are performed
    type(temperature_type) , intent(in)    :: temperature_inst
    type(topo_type)        , intent(in)    :: topo_inst
    !
    ! !LOCAL VARIABLES:
    integer  :: c, l, g, n, fc                   ! indices

    character(len=*), parameter :: subname = 'update_lnd2glc'
    !------------------------------------------------------------------------------

    ! Initialize to reasonable defaults

    this%qice_grc(bounds%begg : bounds%endg, :) = 0._r8
    this%tsrf_grc(bounds%begg : bounds%endg, :) = tfrz
    this%topo_grc(bounds%begg : bounds%endg, :) = 0._r8     
  
    ! Fill the lnd->glc data on the clm grid

    do fc = 1, num_do_smb_c
      c = filter_do_smb_c(fc)
      l = col%landunit(c)
      g = col%gridcell(c) 

      ! Set vertical index based on whether the column in question is glacier or vegetated.  
      if (lun%itype(l) == istice_mec) then
         n = col_itype_to_icemec_class(col%itype(c))
      else if (lun%itype(l) == istsoil) then
         n = 0  !0-level index (bareland information)
      else
         ! Other landunit types do not pass information in the lnd2glc fields.
         ! Note: for this to be acceptable, we need virtual vegetated columns in any grid
         ! cell that is made up solely of glacier plus some other special landunit (e.g.,
         ! glacier + lake) -- otherwise CISM wouldn't have any information for the non-
         ! glaciated portion of the grid cell.
         cycle
      end if

      ! Send surface temperature, topography, and SMB flux (qice) to coupler.
      ! t_soisno and topo_col are valid even in initialization, so tsrf and topo
      ! are set here regardless of the value of init.
      this%tsrf_grc(g,n) = temperature_inst%t_soisno_col(c,1)
      this%topo_grc(g,n) = topo_inst%topo_col(c)

    end do

  end subroutine update_lnd2glc

end module lnd2glcMod

