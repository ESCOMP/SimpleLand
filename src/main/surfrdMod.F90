module surfrdMod

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Contains methods for reading in surface data file and determining
  ! subgrid weights
  !
  ! !NOTES:
  ! TODO Currently reading domain file, although this is done in surfrd.
  ! In NUOPC version we will be reading ESMF mesh file. Until SLIM gets
  ! updated to NUOPC, we are leaving the calls to surfrd unchanged.
  !
  ! !USES:
#include "shr_assert.h"
  use shr_kind_mod    , only : r8 => shr_kind_r8
  use shr_log_mod     , only : errMsg => shr_log_errMsg
  use abortutils      , only : endrun
  use clm_varcon      , only : grlnd
  use clm_varctl      , only : iulog, scmlat, scmlon, single_column
  use ncdio_pio       , only : file_desc_t, var_desc_t, ncd_pio_openfile, ncd_pio_closefile
  use ncdio_pio       , only : ncd_io, check_var, ncd_inqfdims, ncd_inqdid
  use pio
  use spmdMod                         
  !
  ! !PUBLIC TYPES:
  implicit none
  save
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: surfrd_get_globmask  ! Reads global land mask (needed for setting domain decomp)
  public :: surfrd_get_grid      ! Read grid/ladnfrac data into domain (after domain decomp)
  !
  ! !PRIVATE DATA MEMBERS:

  character(len=*), parameter, private :: sourcefile = &
       __FILE__
  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine surfrd_get_globmask(filename, mask, ni, nj)
    !
    ! !DESCRIPTION:
    ! Read the surface dataset grid related information:
    ! This is the first routine called by clm_initialize 
    ! NO DOMAIN DECOMPOSITION  HAS BEEN SET YET
    !
    ! !USES:
    use fileutils , only : getfil
    !
    ! !ARGUMENTS:
    character(len=*), intent(in)    :: filename  ! grid filename
    integer         , pointer       :: mask(:)   ! grid mask 
    integer         , intent(out)   :: ni, nj    ! global grid sizes
    !
    ! !LOCAL VARIABLES:
    logical :: isgrid2d
    integer :: dimid,varid         ! netCDF id's
    integer :: ns                  ! size of grid on file
    integer :: n,i,j               ! index 
    integer :: ier                 ! error status
    type(file_desc_t)  :: ncid     ! netcdf id
    type(var_desc_t)   :: vardesc  ! variable descriptor
    character(len=256) :: varname  ! variable name
    character(len=256) :: locfn    ! local file name
    logical :: readvar             ! read variable in or not
    integer , allocatable :: idata2d(:,:)
    character(len=32) :: subname = 'surfrd_get_globmask' ! subroutine name
    !-----------------------------------------------------------------------

    if (filename == ' ') then
       mask(:) = 1
       RETURN
    end if

    if (masterproc) then
       if (filename == ' ') then
          write(iulog,*) trim(subname),' ERROR: filename must be specified '
          call endrun(msg=errMsg(sourcefile, __LINE__))
       endif
    end if

    call getfil( filename, locfn, 0 )
    call ncd_pio_openfile (ncid, trim(locfn), 0)

    ! Determine dimensions and if grid file is 2d or 1d

    call ncd_inqfdims(ncid, isgrid2d, ni, nj, ns)
    if (masterproc) then
       write(iulog,*)'lat/lon grid flag (isgrid2d) is ',isgrid2d
    end if

    allocate(mask(ns))
    mask(:) = 1

    if (isgrid2d) then
       allocate(idata2d(ni,nj))
       idata2d(:,:) = 1	
       call ncd_io(ncid=ncid, varname='LANDMASK', data=idata2d, flag='read', readvar=readvar)
       if (.not. readvar) then
          call ncd_io(ncid=ncid, varname='mask', data=idata2d, flag='read', readvar=readvar)
       end if
       if (readvar) then
          do j = 1,nj
          do i = 1,ni
             n = (j-1)*ni + i	
             mask(n) = idata2d(i,j)
          enddo
          enddo
       end if
       deallocate(idata2d)
    else
       call ncd_io(ncid=ncid, varname='LANDMASK', data=mask, flag='read', readvar=readvar)
       if (.not. readvar) then
          call ncd_io(ncid=ncid, varname='mask', data=mask, flag='read', readvar=readvar)
       end if
    end if
    if (.not. readvar) call endrun( msg=' ERROR: landmask not on fatmlndfrc file'//errMsg(sourcefile, __LINE__))

    call ncd_pio_closefile(ncid)

  end subroutine surfrd_get_globmask

  !-----------------------------------------------------------------------
  subroutine surfrd_get_grid(begg, endg, ldomain, filename)
    !
    ! !DESCRIPTION:
    ! THIS IS CALLED AFTER THE DOMAIN DECOMPOSITION HAS BEEN CREATED
    ! Read the surface dataset grid related information:
    ! o real latitude  of grid cell (degrees)
    ! o real longitude of grid cell (degrees)
    !
    ! !USES:
    use clm_varcon, only : spval, re
    use domainMod , only : domain_type, domain_init, domain_clean, lon1d, lat1d
    use fileutils , only : getfil
    !
    ! !ARGUMENTS:
    integer          ,intent(in)    :: begg, endg 
    type(domain_type),intent(inout) :: ldomain   ! domain to init
    character(len=*) ,intent(in)    :: filename  ! grid filename
    !
    ! !LOCAL VARIABLES:
    type(file_desc_t) :: ncid               ! netcdf id
    type(var_desc_t)  :: vardesc            ! variable descriptor
    integer :: beg                          ! local beg index
    integer :: end                          ! local end index
    integer :: ni,nj,ns                     ! size of grid on file
    integer :: dimid,varid                  ! netCDF id's
    integer :: start(1), count(1)           ! 1d lat/lon array sections
    integer :: ier,ret                      ! error status
    logical :: readvar                      ! true => variable is on input file 
    logical :: isgrid2d                     ! true => file is 2d lat/lon
    logical :: istype_domain                ! true => input file is of type domain
    real(r8), allocatable :: rdata2d(:,:)   ! temporary
    character(len=16) :: vname              ! temporary
    character(len=256):: locfn              ! local file name
    integer :: n                            ! indices
    real(r8):: eps = 1.0e-12_r8             ! lat/lon error tolerance
    character(len=32) :: subname = 'surfrd_get_grid'     ! subroutine name
!-----------------------------------------------------------------------

    if (masterproc) then
       if (filename == ' ') then
          write(iulog,*) trim(subname),' ERROR: filename must be specified '
          call endrun(msg=errMsg(sourcefile, __LINE__))
       endif
    end if

    call getfil( filename, locfn, 0 )
    call ncd_pio_openfile (ncid, trim(locfn), 0)

    ! Determine dimensions
    call ncd_inqfdims(ncid, isgrid2d, ni, nj, ns)

    ! Determine isgrid2d flag for domain
    call domain_init(ldomain, isgrid2d=isgrid2d, ni=ni, nj=nj, nbeg=begg, nend=endg)

    ! Determine type of file - old style grid file or new style domain file
    call check_var(ncid=ncid, varname='LONGXY', vardesc=vardesc, readvar=readvar) 
    if (readvar) istype_domain = .false.

    call check_var(ncid=ncid, varname='xc', vardesc=vardesc, readvar=readvar) 
    if (readvar) istype_domain = .true.

    ! Read in area, lon, lat

    if (istype_domain) then
       call ncd_io(ncid=ncid, varname= 'area', flag='read', data=ldomain%area, &
            dim1name=grlnd, readvar=readvar)
       ! convert from radians**2 to km**2
       ldomain%area = ldomain%area * (re**2)
       if (.not. readvar) call endrun( msg=' ERROR: area NOT on file'//errMsg(sourcefile, __LINE__))
       
       call ncd_io(ncid=ncid, varname= 'xc', flag='read', data=ldomain%lonc, &
            dim1name=grlnd, readvar=readvar)
       if (.not. readvar) call endrun( msg=' ERROR: xc NOT on file'//errMsg(sourcefile, __LINE__))
       
       call ncd_io(ncid=ncid, varname= 'yc', flag='read', data=ldomain%latc, &
            dim1name=grlnd, readvar=readvar)
       if (.not. readvar) call endrun( msg=' ERROR: yc NOT on file'//errMsg(sourcefile, __LINE__))
    else
       call ncd_io(ncid=ncid, varname= 'AREA', flag='read', data=ldomain%area, &
            dim1name=grlnd, readvar=readvar)
       if (.not. readvar) call endrun( msg=' ERROR: AREA NOT on file'//errMsg(sourcefile, __LINE__))
       
       call ncd_io(ncid=ncid, varname= 'LONGXY', flag='read', data=ldomain%lonc, &
            dim1name=grlnd, readvar=readvar)
       if (.not. readvar) call endrun( msg=' ERROR: LONGXY NOT on file'//errMsg(sourcefile, __LINE__))
       
       call ncd_io(ncid=ncid, varname= 'LATIXY', flag='read', data=ldomain%latc, &
            dim1name=grlnd, readvar=readvar)
       if (.not. readvar) call endrun( msg=' ERROR: LATIXY NOT on file'//errMsg(sourcefile, __LINE__))
    end if

    if (isgrid2d) then
       allocate(rdata2d(ni,nj), lon1d(ni), lat1d(nj))
       if (istype_domain) then
          vname = 'xc'
       else
          vname = 'LONGXY'
       end if
       call ncd_io(ncid=ncid, varname=trim(vname), data=rdata2d, flag='read', readvar=readvar)
       lon1d(:) = rdata2d(:,1)
       if (istype_domain) then
          vname = 'yc'
       else
          vname = 'LATIXY'
       end if
       call ncd_io(ncid=ncid, varname=trim(vname), data=rdata2d, flag='read', readvar=readvar)
       lat1d(:) = rdata2d(1,:)
       deallocate(rdata2d)
    end if

    ! Check lat limited to -90,90

    if (minval(ldomain%latc) < -90.0_r8 .or. &
        maxval(ldomain%latc) >  90.0_r8) then
       write(iulog,*) trim(subname),' WARNING: lat/lon min/max is ', &
            minval(ldomain%latc),maxval(ldomain%latc)
       ! call endrun( msg=' ERROR: lat is outside [-90,90]'//errMsg(sourcefile, __LINE__))
       ! write(iulog,*) trim(subname),' Limiting lat/lon to [-90/90] from ', &
       !     minval(domain%latc),maxval(domain%latc)
       ! where (ldomain%latc < -90.0_r8) ldomain%latc = -90.0_r8
       ! where (ldomain%latc >  90.0_r8) ldomain%latc =  90.0_r8
    endif

    call ncd_io(ncid=ncid, varname='LANDMASK', flag='read', data=ldomain%mask, &
         dim1name=grlnd, readvar=readvar)
    if (.not. readvar) then
       call ncd_io(ncid=ncid, varname='mask', flag='read', data=ldomain%mask, &
            dim1name=grlnd, readvar=readvar)
       if (.not. readvar) then
          call endrun( msg=' ERROR: LANDMASK NOT on fracdata file'//errMsg(sourcefile, __LINE__))
       end if
    end if

    call ncd_io(ncid=ncid, varname='LANDFRAC', flag='read', data=ldomain%frac, &
         dim1name=grlnd, readvar=readvar)
    if (.not. readvar) then
       call ncd_io(ncid=ncid, varname='frac', flag='read', data=ldomain%frac, &
            dim1name=grlnd, readvar=readvar)
       if (.not. readvar) then
          call endrun( msg=' ERROR: LANDFRAC NOT on fracdata file'//errMsg(sourcefile, __LINE__))
       end if
    end if

    call ncd_pio_closefile(ncid)

  end subroutine surfrd_get_grid

end module surfrdMod
