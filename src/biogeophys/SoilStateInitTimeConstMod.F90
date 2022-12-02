module SoilStateInitTimeConstMod

  !------------------------------------------------------------------------------
  ! DESCRIPTION:
  ! Set hydraulic and thermal properties 
  !
  ! !USES
  use SoilStateType , only : soilstate_type
  use LandunitType  , only : lun                
  use ColumnType    , only : col                
  use PatchType     , only : patch                
  !
  implicit none
  private
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public  :: SoilStateInitTimeConst
  !
  ! !PRIVATE DATA:
  ! Control variables (from namelist)

  character(len=*), parameter, private :: sourcefile = &
       __FILE__
  !-----------------------------------------------------------------------
  !
contains

  !-----------------------------------------------------------------------
  subroutine SoilStateInitTimeConst(bounds, soilstate_inst, nlfilename) 
    !
    ! !USES:
    use shr_kind_mod        , only : r8 => shr_kind_r8
    use shr_log_mod         , only : errMsg => shr_log_errMsg
    use shr_infnan_mod      , only : nan => shr_infnan_nan, assignment(=)
    use decompMod           , only : bounds_type
    use abortutils          , only : endrun
    use spmdMod             , only : masterproc
    use ncdio_pio           , only : file_desc_t, ncd_io, ncd_double, ncd_int, ncd_inqvdlen
    use ncdio_pio           , only : ncd_pio_openfile, ncd_pio_closefile, ncd_inqdlen
    use clm_varpar          , only : numpft, numrad
    use clm_varpar          , only : nlevsoi, nlevgrnd, nlevlak, nlevsoifl, nlayer, nlayert, nlevurb, nlevsno
    use clm_varcon          , only : zsoi, dzsoi, zisoi, spval
    use clm_varcon          , only : secspday, denh2o, denice, grlnd
    use clm_varctl          , only : iulog, fsurdat, paramfile
    use landunit_varcon     , only : istdlak, istwet, istsoil, istcrop, istice_mec
    use column_varcon       , only : icol_road_perv, icol_road_imperv 
    use fileutils           , only : getfil
    use GridcellType     , only : grc                
    !
    ! !ARGUMENTS:
    type(bounds_type)    , intent(in)    :: bounds  
    type(soilstate_type) , intent(inout) :: soilstate_inst
    character(len=*)     , intent(in)    :: nlfilename ! Namelist filename
    !
    ! !LOCAL VARIABLES:
    integer            :: p, lev, c, l, g, j            ! indices
    integer            :: dimid                         ! dimension id
    logical            :: readvar 
    type(file_desc_t)  :: ncid                          ! netcdf id
    real(r8) ,pointer  :: zsoifl (:)                    ! Output: [real(r8) (:)]  original soil midpoint 
    real(r8) ,pointer  :: zisoifl (:)                   ! Output: [real(r8) (:)]  original soil interface depth 
    real(r8) ,pointer  :: dzsoifl (:)                   ! Output: [real(r8) (:)]  original soil thickness 
    real(r8) ,pointer  :: gti (:)                       ! read in - fmax 
    character(len=256) :: locfn                         ! local filename
    integer            :: begp, endp
    integer            :: begc, endc
    integer            :: begg, endg
    !-----------------------------------------------------------------------

    begp = bounds%begp; endp= bounds%endp
    begc = bounds%begc; endc= bounds%endc
    begg = bounds%begg; endg= bounds%endg

    do c = begc,endc
       soilstate_inst%smpmin_col(c) = -1.e8_r8
    end do

    ! --------------------------------------------------------------------
    ! Initialize root fraction (computing from surface, d is depth in meter):
    ! --------------------------------------------------------------------

    do c = bounds%begc,bounds%endc
       l = col%landunit(c)
       if (lun%itype(l) == istsoil .or. lun%itype(l) == istcrop) then
          soilstate_inst%rootfr_col (c,nlevsoi+1:nlevgrnd) = 0._r8
       else
          ! Inactive CH4 columns 
          ! (Also includes (lun%itype(l)==istdlak .and.  allowlakeprod), which used to be
          ! in a separate branch of the conditional)
          soilstate_inst%rootfr_col (c,:) = spval
       end if
    end do

    ! --------------------------------------------------------------------
    ! Read surface dataset
    ! --------------------------------------------------------------------

    if (masterproc) then
       write(iulog,*) 'Attempting to read soil color boundary data .....'
    end if

    call getfil (fsurdat, locfn, 0)
    call ncd_pio_openfile (ncid, locfn, 0)

    ! Read fmax

    allocate(gti(begg:endg))
    call ncd_io(ncid=ncid, varname='FMAX', flag='read', data=gti, dim1name=grlnd, readvar=readvar)
    if (.not. readvar) then
       call endrun(msg=' ERROR: FMAX NOT on surfdata file'//errMsg(sourcefile, __LINE__)) 
    end if
    do c = begc, endc
       g = col%gridcell(c)
       soilstate_inst%wtfact_col(c) = gti(g)
    end do
    deallocate(gti)

    ! Close file

    call ncd_pio_closefile(ncid)

    ! --------------------------------------------------------------------
    ! get original soil depths to be used in interpolation of sand and clay
    ! --------------------------------------------------------------------

    allocate(zsoifl(1:nlevsoifl), zisoifl(0:nlevsoifl), dzsoifl(1:nlevsoifl))
    do j = 1, nlevsoifl
       zsoifl(j) = 0.025*(exp(0.5_r8*(j-0.5_r8))-1._r8)    !node depths
    enddo

    dzsoifl(1) = 0.5_r8*(zsoifl(1)+zsoifl(2))             !thickness b/n two interfaces
    do j = 2,nlevsoifl-1
       dzsoifl(j)= 0.5_r8*(zsoifl(j+1)-zsoifl(j-1))
    enddo
    dzsoifl(nlevsoifl) = zsoifl(nlevsoifl)-zsoifl(nlevsoifl-1)

    zisoifl(0) = 0._r8
    do j = 1, nlevsoifl-1
       zisoifl(j) = 0.5_r8*(zsoifl(j)+zsoifl(j+1))         !interface depths
    enddo
    zisoifl(nlevsoifl) = zsoifl(nlevsoifl) + 0.5_r8*dzsoifl(nlevsoifl)

    ! --------------------------------------------------------------------
    ! Set soil hydraulic and thermal properties: non-lake
    ! --------------------------------------------------------------------

    !   urban roof, sunwall and shadewall thermal properties used to 
    !   derive thermal conductivity and heat capacity are set to special 
    !   value because thermal conductivity and heat capacity for urban 
    !   roof, sunwall and shadewall are prescribed in SoilThermProp.F90 
    !   in SoilPhysicsMod.F90


    do c = begc, endc
       g = col%gridcell(c)
       l = col%landunit(c)

       if (lun%itype(l)==istwet .or. lun%itype(l)==istice_mec) then

          do lev = 1,nlevgrnd
             soilstate_inst%watsat_col(c,lev) = spval
          end do

       else if (lun%urbpoi(l) .and. (col%itype(c) /= icol_road_perv) .and. (col%itype(c) /= icol_road_imperv) )then

          ! Urban Roof, sunwall, shadewall properties set to special value
          do lev = 1,nlevgrnd
             soilstate_inst%watsat_col(c,lev) = spval
          end do

       else

          do lev = 1,nlevgrnd
             if (lun%itype(l) /= istdlak) then  ! soil columns of both urban and non-urban types
                ! TODO slevis: Temporary during dismantling for SLIM
                soilstate_inst%watsat_col(c,lev) = 1._r8
             end if
          end do
       end if
    end do

    ! TODO slevis: Temporary while dismantling for SLIM
    ! --------------------------------------------------------------------
    ! Set soil hydraulic and thermal properties: lake
    ! --------------------------------------------------------------------

    do c = begc, endc
       g = col%gridcell(c)
       l = col%landunit(c)
       if (lun%itype(l)==istdlak) then
          do lev = 1,nlevgrnd
             soilstate_inst%watsat_col(c,lev) = 0.489_r8
          end do
       endif
    end do

    ! --------------------------------------------------------------------
    ! Deallocate memory
    ! --------------------------------------------------------------------

    deallocate(zisoifl, zsoifl, dzsoifl)

  end subroutine SoilStateInitTimeConst

end module SoilStateInitTimeConstMod
