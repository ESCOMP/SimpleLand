module UrbanParamsType

  !------------------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Urban Constants
  !
  ! !USES:
  use shr_kind_mod , only : r8 => shr_kind_r8
  use shr_log_mod  , only : errMsg => shr_log_errMsg
  use abortutils   , only : endrun
  use decompMod    , only : bounds_type
  use clm_varctl   , only : iulog, fsurdat
  use clm_varcon   , only : grlnd
  use LandunitType , only : lun                
  !
  implicit none
  save
  private
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public  :: UrbanInput        ! Read in urban input data
  !
  ! !PRIVATE TYPE
  type urbinp_type
     real(r8), pointer :: canyon_hwr      (:,:)  
     real(r8), pointer :: wtlunit_roof    (:,:)  
     real(r8), pointer :: wtroad_perv     (:,:)  
     real(r8), pointer :: thick_wall      (:,:)
     real(r8), pointer :: thick_roof      (:,:)
  end type urbinp_type
  type (urbinp_type), public :: urbinp   ! urban input derived type

  ! !PUBLIC TYPE
  type, public :: urbanparams_type

     real(r8), pointer     :: thick_wall          (:)   ! lun total thickness of urban wall (m)
     real(r8), pointer     :: thick_roof          (:)   ! lun total thickness of urban roof (m)

   contains

     procedure, public :: Init 
     
  end type urbanparams_type
  !
  ! !Urban control variables
  character(len= *), parameter, public :: urban_wasteheat_on = 'ON_WASTEHEAT'  

  ! !PRIVATE MEMBER DATA:

  character(len=*), parameter, private :: sourcefile = &
       __FILE__
  !----------------------------------------------------------------------- 

contains

  !-----------------------------------------------------------------------
  subroutine Init(this, bounds)
    !
    ! Allocate module variables and data structures
    !
    ! !USES:
    use shr_infnan_mod  , only : nan => shr_infnan_nan, assignment(=)
    use landunit_varcon , only : isturb_MIN
    !
    ! !ARGUMENTS:
    class(urbanparams_type) :: this
    type(bounds_type)      , intent(in)    :: bounds  
    !
    ! !LOCAL VARIABLES:
    integer :: l,c,g  ! indices
    integer :: dindx  ! urban density type index
    integer :: begl, endl
    integer :: begg, endg
    !---------------------------------------------------------------------

    begl = bounds%begl; endl = bounds%endl
    begg = bounds%begg; endg = bounds%endg

    ! Allocate urbanparams data structure
    allocate(this%thick_wall          (begl:endl))          ; this%thick_wall          (:)   = nan
    allocate(this%thick_roof          (begl:endl))          ; this%thick_roof          (:)   = nan

    ! Initialize time constant urban variables

    do l = begl,endl

       ! "0" refers to urban wall/roof surface and "nlevsoi" refers to urban wall/roof bottom
       if (lun%urbpoi(l)) then

          g = lun%gridcell(l)
          dindx = lun%itype(l) - isturb_MIN + 1

          ! Landunit level initialization for urban wall and roof layers and interfaces

          lun%canyon_hwr(l)   = urbinp%canyon_hwr(g,dindx)
          lun%wtroad_perv(l)  = urbinp%wtroad_perv(g,dindx)
          lun%wtlunit_roof(l) = urbinp%wtlunit_roof(g,dindx)

          this%thick_wall(l)     = urbinp%thick_wall(g,dindx)
          this%thick_roof(l)     = urbinp%thick_roof(g,dindx)

       end if
    end do

    call UrbanInput(bounds%begg, bounds%endg, mode='finalize')

  end subroutine Init

  !-----------------------------------------------------------------------
  subroutine UrbanInput(begg, endg, mode)
    !
    ! !DESCRIPTION: 
    ! Allocate memory and read in urban input data
    !
    ! !USES:
    use clm_varpar      , only : numrad, nlevurb
    use landunit_varcon , only : numurbl
    use fileutils       , only : getavu, relavu, getfil, opnfil
    use spmdMod         , only : masterproc
    use domainMod       , only : ldomain
    use ncdio_pio       , only : file_desc_t, ncd_io, ncd_inqvdlen, ncd_inqfdims 
    use ncdio_pio       , only : ncd_pio_openfile, ncd_pio_closefile, ncd_inqdid, ncd_inqdlen
    !
    ! !ARGUMENTS:
    implicit none
    integer, intent(in) :: begg, endg
    character(len=*), intent(in) :: mode
    !
    ! !LOCAL VARIABLES:
    character(len=256) :: locfn      ! local file name
    type(file_desc_t)  :: ncid       ! netcdf id
    integer :: dimid                 ! netCDF id
    integer :: nw,n,k,i,j,ni,nj,ns   ! indices
    integer :: nlevurb_i             ! input grid: number of urban vertical levels
    integer :: numrad_i              ! input grid: number of solar bands (VIS/NIR)
    integer :: numurbl_i             ! input grid: number of urban landunits
    integer :: ier,ret               ! error status
    logical :: isgrid2d              ! true => file is 2d 
    logical :: readvar               ! true => variable is on dataset
    logical :: has_numurbl           ! true => numurbl dimension is on dataset
    character(len=32) :: subname = 'UrbanInput' ! subroutine name
    !-----------------------------------------------------------------------

    if ( nlevurb == 0 ) return

    if (mode == 'initialize') then

       ! Read urban data
       
       if (masterproc) then
          write(iulog,*)' Reading in urban input data from fsurdat file ...'
       end if
       
       call getfil (fsurdat, locfn, 0)
       call ncd_pio_openfile (ncid, locfn, 0)

       if (masterproc) then
          write(iulog,*) subname,trim(fsurdat)
       end if

       ! Check whether this file has new-format urban data
       call ncd_inqdid(ncid, 'numurbl', dimid, dimexist=has_numurbl)

       ! If file doesn't have numurbl, then it is old-format urban;
       ! in this case, set nlevurb to zero
       if (.not. has_numurbl) then
         nlevurb = 0
         if (masterproc) write(iulog,*)'PCT_URBAN is not multi-density, nlevurb set to 0'
       end if

       if ( nlevurb == 0 ) return

       ! Allocate dynamic memory
       allocate(urbinp%canyon_hwr(begg:endg, numurbl), &  
                urbinp%wtlunit_roof(begg:endg, numurbl), &  
                urbinp%wtroad_perv(begg:endg, numurbl), &
                urbinp%thick_wall(begg:endg, numurbl), &
                urbinp%thick_roof(begg:endg, numurbl), &
                stat=ier)
       if (ier /= 0) then
          call endrun(msg="Allocation error "//errmsg(sourcefile, __LINE__))
       endif

       call ncd_inqfdims (ncid, isgrid2d, ni, nj, ns)
       if (ldomain%ns /= ns .or. ldomain%ni /= ni .or. ldomain%nj /= nj) then
          write(iulog,*)trim(subname), 'ldomain and input file do not match dims '
          write(iulog,*)trim(subname), 'ldomain%ni,ni,= ',ldomain%ni,ni
          write(iulog,*)trim(subname), 'ldomain%nj,nj,= ',ldomain%nj,nj
          write(iulog,*)trim(subname), 'ldomain%ns,ns,= ',ldomain%ns,ns
          call endrun(msg=errmsg(sourcefile, __LINE__))
       end if

       call ncd_inqdid(ncid, 'nlevurb', dimid)
       call ncd_inqdlen(ncid, dimid, nlevurb_i)
       if (nlevurb_i /= nlevurb) then
          write(iulog,*)trim(subname)// ': parameter nlevurb= ',nlevurb, &
               'does not equal input dataset nlevurb= ',nlevurb_i
          call endrun(msg=errmsg(sourcefile, __LINE__))
       endif

       call ncd_inqdid(ncid, 'numrad', dimid)
       call ncd_inqdlen(ncid, dimid, numrad_i)
       if (numrad_i /= numrad) then
          write(iulog,*)trim(subname)// ': parameter numrad= ',numrad, &
               'does not equal input dataset numrad= ',numrad_i
          call endrun(msg=errmsg(sourcefile, __LINE__))
       endif
       call ncd_inqdid(ncid, 'numurbl', dimid)
       call ncd_inqdlen(ncid, dimid, numurbl_i)
       if (numurbl_i /= numurbl) then
          write(iulog,*)trim(subname)// ': parameter numurbl= ',numurbl, &
               'does not equal input dataset numurbl= ',numurbl_i
          call endrun(msg=errmsg(sourcefile, __LINE__))
       endif
       call ncd_io(ncid=ncid, varname='CANYON_HWR', flag='read', data=urbinp%canyon_hwr,&
            dim1name=grlnd, readvar=readvar)
       if (.not. readvar) then
          call endrun( msg='ERROR: CANYON_HWR NOT on fsurdat file '//errmsg(sourcefile, __LINE__))
       end if

       call ncd_io(ncid=ncid, varname='WTLUNIT_ROOF', flag='read', data=urbinp%wtlunit_roof, &
            dim1name=grlnd,  readvar=readvar)
       if (.not. readvar) then
          call endrun( msg=' ERROR: WTLUNIT_ROOF NOT on fsurdat file'//errmsg(sourcefile, __LINE__))
       end if

       call ncd_io(ncid=ncid, varname='WTROAD_PERV', flag='read', data=urbinp%wtroad_perv, &
            dim1name=grlnd, readvar=readvar)
       if (.not. readvar) then
          call endrun( msg=' ERROR: WTROAD_PERV NOT on fsurdat file'//errmsg(sourcefile, __LINE__))
       end if

       call ncd_io(ncid=ncid, varname='THICK_WALL', flag='read', data=urbinp%thick_wall, &
            dim1name=grlnd, readvar=readvar)
       if (.not. readvar) then
          call endrun( msg=' ERROR: THICK_WALL NOT on fsurdat file'//errmsg(sourcefile, __LINE__))
       end if

       call ncd_io(ncid=ncid, varname='THICK_ROOF', flag='read', data=urbinp%thick_roof, &
            dim1name=grlnd, readvar=readvar)
       if (.not. readvar) then
          call endrun( msg=' ERROR: THICK_ROOF NOT on fsurdat file'//errmsg(sourcefile, __LINE__))
       end if

       call ncd_pio_closefile(ncid)
       if (masterproc) then
          write(iulog,*)' Sucessfully read urban input data' 
          write(iulog,*)
       end if

    else if (mode == 'finalize') then

       if ( nlevurb == 0 ) return

       deallocate(urbinp%canyon_hwr, &
                  urbinp%wtlunit_roof, &
                  urbinp%wtroad_perv, &
                  urbinp%thick_wall, &
                  urbinp%thick_roof, &
                  stat=ier)
       if (ier /= 0) then
          call endrun(msg='initUrbanInput: deallocation error '//errmsg(sourcefile, __LINE__))
       end if
    else
       write(iulog,*)'initUrbanInput error: mode ',trim(mode),' not supported '
       call endrun(msg=errmsg(sourcefile, __LINE__))
    end if

  end subroutine UrbanInput

end module UrbanParamsType




