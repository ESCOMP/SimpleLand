module initInterpMod

  !----------------------------------------------------------------------- 
  ! Interpolate initial conditions file from one resolution and/or landmask
  ! to another resolution and/or landmask
  !----------------------------------------------------------------------- 

#include "shr_assert.h"
  use initInterpBounds, only : interp_bounds_type
  use initInterpMindist, only: set_mindist, subgrid_type
  use initInterp1dData, only : interp_1d_data
  use initInterp2dvar, only: interp_2dvar_type
  use initInterpMultilevelBase, only : interp_multilevel_type
  use initInterpMultilevelContainer, only : interp_multilevel_container_type
  use shr_kind_mod   , only: r8 => shr_kind_r8, r4 => shr_kind_r4
  use shr_const_mod  , only: SHR_CONST_PI, SHR_CONST_REARTH
  use shr_sys_mod    , only: shr_sys_flush
  use shr_log_mod    , only : errMsg => shr_log_errMsg
  use shr_string_mod , only: shr_string_listGetName
  use clm_varctl     , only: iulog
  use abortutils     , only: endrun
  use spmdMod        , only: masterproc
  use restUtilMod    , only: iflag_interp, iflag_copy, iflag_skip
  use ncdio_utils    , only: find_var_on_file
  use ncdio_pio
  use pio

  implicit none
  private
  save

  ! Public methods

  public :: initInterp

  ! Private methods

  private :: check_dim_subgrid
  private :: check_dim_level
  private :: findMinDist
  private :: set_subgrid_info
  private :: interp_0d_copy
  private :: interp_1d_double
  private :: interp_1d_int
  private :: interp_2d_double

  ! Private data
 
  character(len=*), parameter, private :: sourcefile = &
       __FILE__

contains

  !=======================================================================

  subroutine initInterp (filei, fileo, bounds)

    !----------------------------------------------------------------------- 
    ! Read initial data from netCDF instantaneous initial data history file 
    !-----------------------------------------------------------------------

    use decompMod, only: bounds_type

    ! --------------------------------------------------------------------
    ! arguments
    character(len=*)  , intent(in) :: filei     !input  initial dataset
    character(len=*)  , intent(in) :: fileo    !output initial dataset
    type(bounds_type) , intent(in) :: bounds
    !
    ! local variables
    integer            :: i,j,k,l,m,n     ! loop indices    
    integer            :: begi, endi      ! beginning/ending indices 
    integer            :: bego, endo      ! beginning/ending indices
    type(interp_bounds_type) :: bounds_i  ! input file bounds
    type(interp_bounds_type) :: bounds_o  ! output file bounds
    integer            :: nlevi,nlevo     ! input/output number of levels
    type(file_desc_t), target :: ncidi, ncido    ! input/output pio fileids 
    integer            :: dimleni,dimleno ! input/output dimension length       
    integer            :: nvars           ! number of variables
    character(len=256) :: varname         ! variable name
    character(len=256) :: varname_i_options ! possible variable names on input file
    character(len=256) :: varname_i       ! variable name on input file
    character(len=16)  :: vec_dimname     ! subgrid dimension name
    type(Var_desc_t)   :: vardesc         ! pio variable descriptor
    integer            :: xtypeo          ! netCDF variable type
    integer            :: varido          ! netCDF variable id
    integer            :: ndimso          ! netCDF number of dimensions
    integer            :: dimidso(3) = -1 ! netCDF dimension ids
    integer            :: status          ! return code
    integer            :: iflag_interpinic
    real(r8)           :: rvalue
    integer            :: ivalue 
    integer            :: spinup_state_i, spinup_state_o
    integer            :: decomp_cascade_state_i, decomp_cascade_state_o 
    integer            :: ngrcsi 
    integer            :: ngrcso 
    integer , pointer  :: grcindx(:) 
    logical , pointer  :: grc_activei(:), grc_activeo(:)
    integer , pointer  :: sgridindex(:)
    type(interp_multilevel_container_type) :: interp_multilevel_container
    type(interp_2dvar_type) :: var2d_i, var2d_o  ! holds metadata for 2-d variables
    !--------------------------------------------------------------------

    if (masterproc) then
       write (iulog,*) '**** Mapping clm initial data from input ',trim(filei),&
            '  to output ',trim(fileo),' ****'
    end if

    ! --------------------------------------------
    ! Open input and output initial conditions files (both just for reading now)
    ! --------------------------------------------

    call ncd_pio_openfile (ncidi, trim(filei) , 0)
    call ncd_pio_openfile (ncido, trim(fileo),  ncd_write)

    ! --------------------------------------------
    ! Determine dimensions and error checks on dimensions
    ! --------------------------------------------

    call check_dim_subgrid(ncidi, ncido, dimname ='gridcell', dimleni=ngrcsi, dimleno=ngrcso)

    if (masterproc) then
       write (iulog,*) 'input gridcells = ',ngrcsi,' output gridcells = ',ngrcso
    end if

    ! NOTE(wjs, 2015-10-31) The inclusion of must_be_same in these checks essentially
    ! duplicates checks done elsewhere - specifically, if the interpolator for a given
    ! level dimension is the copy interpolator, then that will do its own checks to
    ! ensure that the dimension sizes match. So we may want to remove the must_be_same
    ! argument, and make check_dim_level purely informational, in order to remove this
    ! maintenance problem - or maybe even remove check_dim_level entirely.
    call check_dim_level(ncidi, ncido, dimname='levgrnd', must_be_same=.false.)
    call check_dim_level(ncidi, ncido, dimname='numrad' , must_be_same=.true.)

    ! --------------------------------------------
    ! Find closest values for gridcells
    ! --------------------------------------------

    bounds_i = interp_bounds_type(begg = 1, endg = ngrcsi)

    bounds_o = interp_bounds_type(begg = bounds%begg, endg = bounds%endg)

    allocate(grc_activei(bounds_i%get_begg():bounds_i%get_endg()))

    allocate(grc_activeo(bounds_o%get_begg():bounds_o%get_endg()))

    allocate(grcindx(bounds_o%get_begg():bounds_o%get_endg()))

    ! For each output gridcell, find the input gridcell, grcindx, that is closest

    if (masterproc) then
       write(iulog,*)'finding minimum distance for gridcells'
    end if
    vec_dimname = 'gridcell'
    call findMinDist(vec_dimname, bounds_i%get_begg(), bounds_i%get_endg(), &
         bounds_o%get_begg(), bounds_o%get_endg(), ncidi, ncido, &
         grc_activei, grc_activeo, grcindx)

    ! ------------------------------------------------------------------------
    ! Set up interpolators for multi-level variables
    ! ------------------------------------------------------------------------

!   if (masterproc) then
!      write(iulog,*)'setting up interpolators for multi-level variables'
!   end if
!   interp_multilevel_container = interp_multilevel_container_type( &
!        ncid_source = ncidi, ncid_dest = ncido, &
!        bounds_source = bounds_i, bounds_dest = bounds_o)

    !------------------------------------------------------------------------          
    ! Read input initial data and write output initial data
    !------------------------------------------------------------------------          

    ! Only examing the snow interfaces above zi=0 => zisno and zsno have
    ! the same level dimension below

    ! Read input initial data and write output initial data
    ! Only examing the snow interfaces above zi=0 => zisno and zsno have
    ! the same level dimension below

    if (masterproc) then
       write(iulog,*)'reading in initial dataset'
    end if
    ! Get number of output variables and loop over them 
    status = pio_inquire(ncido, nVariables=nvars)
    do varido = 1, nvars

       !---------------------------------------------------          
       ! Given varido, get out variable data 
       !---------------------------------------------------          

       status = pio_inquire_variable(ncido, varid=varido, name=varname, &
            xtype=xtypeo, ndims=ndimso, dimids=dimidso)

       !---------------------------------------------------          
       ! If variable is zsoi, SKIP this variable
       !---------------------------------------------------          

       if  ( trim(varname) == 'zsoi' ) then
          if (masterproc) then
             write(iulog,*) 'Skipping     : ',trim(varname)
          end if
          CYCLE
       end if

       !---------------------------------------------------          
       ! If interpinic flag is set to skip on output file
       ! SKIP this variable
       !---------------------------------------------------          

       status = pio_inq_varid (ncido, trim(varname), vardesc)
       status = pio_get_att(ncido, vardesc, 'interpinic_flag', iflag_interpinic)
       if (iflag_interpinic == iflag_skip) then
          if (masterproc) then
             write (iulog,*) 'Skipping     : ', trim(varname)
          end if
          CYCLE
       end if

       !---------------------------------------------------          
       ! Read metadata telling us possible variable names on input file;
       ! determine which of these is present on the input file
       !---------------------------------------------------          

       ! 'varnames_on_old_files' is a colon-delimited list of possible variable names,
       ! enabling backwards compatibility with old input files
       status = pio_get_att(ncido, vardesc, 'varnames_on_old_files', varname_i_options)

       ! We expect the first name in the list to match the current variable name. If that
       ! isn't true, abort. This check is mainly to catch behavior changes in the code to
       ! write this attribute in restUtilMod: Significant changes in behavior there need
       ! to be matched by corresponding changes in this module. For example, if
       ! restUtilMod changes this attribute to exclude the first name (which, after all,
       ! is available via the variable name itself), then code in this module should
       ! change accordingly.
       call shr_string_listGetName(varname_i_options, 1, varname_i)
       if (varname_i /= varname) then
          if (masterproc) then
             write(iulog,*) 'ERROR: expect first element in varnames_on_old_files to match varname'
             write(iulog,*) 'varname = ', trim(varname)
             write(iulog,*) 'varnames_on_old_files = ', trim(varname_i_options)
             write(iulog,*) 'This likely indicates that the code to write the'
             write(iulog,*) 'varnames_on_old_files attribute list has changed behavior.'
          end if
          call endrun(msg='ERROR: expect first element in varnames_on_old_files to match varname'// &
               errMsg(sourcefile, __LINE__))
       end if

       ! Find which of the list of possible variables actually exists on the input file.
       call find_var_on_file(ncidi, varname_i_options, varname_i)

       ! Note that, if none of the options are found, varname_i will be set to the first
       ! variable in the list, in which case the following code will determine that we
       ! should skip this variable.

       !---------------------------------------------------
       ! If variable is on output file - but not on input file
       ! SKIP this variable
       !---------------------------------------------------          

       call pio_seterrorhandling(ncidi, PIO_BCAST_ERROR)
       status = pio_inq_varid(ncidi, name=varname_i, vardesc=vardesc)
       call pio_seterrorhandling(ncidi, PIO_INTERNAL_ERROR)
       if (status /= PIO_noerr) then
          if (masterproc) then
             write (iulog,*) 'Skipping     : ', trim(varname), ' variable is NOT on input file'
          end if
          CYCLE
       end if

       !---------------------------------------------------          
       ! For scalar outut variables 
       !---------------------------------------------------          

       if ( ndimso == 0 ) then

          if  ( trim(varname) .eq. 'spinup_state' ) then

             ! since we are copying soil variables, need to also copy spinup state 
             ! since otherwise if they are different then it will break the spinup procedure
             status = pio_inq_varid(ncidi, trim(varname_i), vardesc)
             status = pio_get_var(ncidi, vardesc, spinup_state_i)
             status = pio_inq_varid(ncido, trim(varname), vardesc)
             status = pio_get_var(ncido, vardesc, spinup_state_o)
             if ( spinup_state_i  /=  spinup_state_o ) then
                if (masterproc) then
                   write (iulog,*) 'Spinup states are different: Copying: ', &
                        trim(varname_i), ' => ', trim(varname)
                end if
                status = pio_put_var(ncido, vardesc, spinup_state_i)
             else
                if (masterproc) then
                   write (iulog,*) 'Spinup states match: Skipping: ', &
                        trim(varname_i), ' => ', trim(varname)
                end if
             endif
             
          else if  ( trim(varname) .eq. 'decomp_cascade_state' ) then

             ! ditto for the decomposition cascade
             status = pio_inq_varid(ncidi, trim(varname_i), vardesc)
             status = pio_get_var(ncidi, vardesc, decomp_cascade_state_i)
             status = pio_inq_varid(ncido, trim(varname), vardesc)
             status = pio_get_var(ncido, vardesc, decomp_cascade_state_o)
             if ( decomp_cascade_state_i  /=  decomp_cascade_state_o ) then
                call endrun(msg='ERROR: Decomposition cascade states are different'//errMsg(sourcefile, __LINE__))
             else
                if (masterproc) then
                   write (iulog,*) 'Decomposition cascade states match: Skipping: ', &
                        trim(varname_i), ' => ', trim(varname)
                end if
             endif

          else if (iflag_interpinic == iflag_copy) then

             call interp_0d_copy(varname, varname_i, xtypeo, ncidi, ncido)

          else if (iflag_interpinic == iflag_skip) then

             if (masterproc) then
                write(iulog,*) 'Skipping     : ',trim(varname)
             end if

          end if

       !---------------------------------------------------          
       ! For 1D output variables
       !---------------------------------------------------          

       else if ( ndimso == 1 ) then

          status = pio_inq_dimname(ncido, dimidso(1), vec_dimname)
          begi = bounds_i%get_beg(vec_dimname)
          endi = bounds_i%get_end(vec_dimname)
          bego = bounds_o%get_beg(vec_dimname)
          endo = bounds_o%get_end(vec_dimname)
          if ( vec_dimname == 'gridcell' )then
             sgridindex => grcindx
          else
             call endrun(msg='ERROR interpinic: 1D variable '//trim(varname)//&
                  'with unknown subgrid dimension: '//trim(vec_dimname)//&
                  errMsg(sourcefile, __LINE__))
          end if

          if ( xtypeo == pio_int )then
             call interp_1d_int ( varname, varname_i, vec_dimname, begi, endi, bego, endo, &
                  ncidi, ncido, sgridindex )
          else if ( xtypeo == pio_double )then                                 
             call interp_1d_double( varname, varname_i, vec_dimname, begi, endi, bego, endo, &
                  ncidi, ncido, sgridindex )
          else
             call endrun(msg='ERROR interpinic: 1D variable with unknown type: '//&
                  trim(varname)//errMsg(sourcefile, __LINE__))
          end if

       !---------------------------------------------------          
       ! For 2D output variables
       !---------------------------------------------------          

       else if ( ndimso == 2 )then

          if ( xtypeo /= pio_double )then
             call endrun(msg='ERROR interpinic: 2D variable with unknown type: '//&
                  trim(varname)//errMsg(sourcefile, __LINE__))
          end if

          var2d_i = interp_2dvar_type( &
               varname = varname_i, &
               ncid = ncidi, &
               file_is_dest = .false., &
               bounds = bounds_i)

          var2d_o = interp_2dvar_type( &
               varname = varname, &
               ncid = ncido, &
               file_is_dest = .true., &
               bounds = bounds_o)

          vec_dimname = var2d_o%get_vec_dimname()
          begi = bounds_i%get_beg(vec_dimname)
          endi = bounds_i%get_end(vec_dimname)
          bego = bounds_o%get_beg(vec_dimname)
          endo = bounds_o%get_end(vec_dimname)
          if ( vec_dimname == 'gridcell' )then
             sgridindex => grcindx
          else
             call endrun(msg='ERROR interpinic: 2D variable with unknown subgrid dimension: '//&
                  trim(varname)//errMsg(sourcefile, __LINE__))
          end if
          call interp_2d_double(var2d_i, var2d_o, &
               begi, endi, bego, endo, &
               sgridindex, &
               interp_multilevel_container)

       else

          call endrun(msg='ERROR interpinic: variable NOT scalar, 1D or 2D: '//&
               trim(varname)//errMsg(sourcefile, __LINE__))

       end if
       call shr_sys_flush(iulog)

    end do
    ! Close input file

    call pio_closefile(ncidi)

    ! Do some final cleanup of specific variables
    !
    ! NOTE(wjs, 2015-10-31) I really don't like having this variable-specific logic here,
    ! but I can't see a great way around this. (One alternative would be to do this
    ! cleanup when reading the restart file, e.g., in subgridRestMod. But (1) that feels
    ! like a hidden way to accomplish this, and (2) that would apply for *any* restart
    ! read, as opposed to just restart reads following init_interp, which could mask
    ! problems with other restart files.)

    ! Need to first sync the file so the previous writes complete before we try to re-read
    ! variables
    call pio_syncfile(ncido)

    ! BUG: (EBK, 2016-1-6, bugz: 2261) PIO2 in cime4.3.9 has a bug where you can't do 
    ! two writes of the same variable to a file. So we have to close and reopen the file. 
    ! When PIO2 is robust enough to handle that we can remove the following two lines.
    call pio_closefile(ncido)
    call ncd_pio_openfile (ncido, trim(fileo),  ncd_write)
    ! ENDBUG:

    if (masterproc) then
       write(iulog,*) 'Cleaning up / adjusting variables'
    end if

    ! Close output file

    call pio_closefile(ncido)

    if (masterproc) then
       write (iulog,*) ' Successfully created initial condition file mapped from input IC file'
    end if

  end subroutine initInterp

  !=======================================================================

  subroutine findMinDist( dimname, begi, endi, bego, endo, ncidi, ncido, &
       activei, activeo, minindx)

    ! --------------------------------------------------------------------
    !
    ! Find the PATCH distances based on the column distances already calculated
    !
    ! arguments
    character(len=*)  , intent(inout) :: dimname
    integer           , intent(in)    :: begi, endi
    integer           , intent(in)    :: bego, endo
    type(file_desc_t) , intent(inout) :: ncidi         
    type(file_desc_t) , intent(inout) :: ncido         
    logical           , intent(out)   :: activei(begi:endi)
    logical           , intent(out)   :: activeo(bego:endo)
    integer           , intent(out)   :: minindx(bego:endo)         
    !
    ! local variables 
    type(subgrid_type)   :: subgridi 
    type(subgrid_type)   :: subgrido 
    ! --------------------------------------------------------------------

    if (masterproc) then
       write(iulog,*)'calling set_subgrid_info for ',trim(dimname), ' for input'
    end if
    call set_subgrid_info(beg=begi, end=endi, dimname=dimname, use_glob=.true., &
         ncid=ncidi, active=activei, subgrid=subgridi)

    if (masterproc) then
       write(iulog,*)'calling set_subgrid_info for ',trim(dimname), ' for output'
    end if
    call set_subgrid_info(beg=bego, end=endo, dimname=dimname, use_glob=.false., &
         ncid=ncido, active=activeo, subgrid=subgrido)

    if (masterproc) then
       write(iulog,*)'calling set_mindist for ',trim(dimname)
    end if
    call set_mindist(begi, endi, bego, endo, activei, activeo, subgridi, subgrido, &
         minindx)

    deallocate(subgridi%lat, subgridi%lon, subgridi%coslat)
    deallocate(subgrido%lat, subgrido%lon, subgrido%coslat)
    
  end subroutine findMinDist

 !=======================================================================

  subroutine set_subgrid_info(beg, end, dimname, use_glob, ncid, active, subgrid)

    ! --------------------------------------------------------------------
    ! arguments
    integer            , intent(in)    :: beg, end
    type(file_desc_t)  , intent(inout) :: ncid
    character(len=*)   , intent(in)    :: dimname
    logical            , intent(in)    :: use_glob  ! if .true., use the 'glob' form of ncd_io
    logical            , intent(out)   :: active(beg:end)    
    type(subgrid_type) , intent(inout) :: subgrid
    !
    ! local variables
    integer              :: n
    integer, pointer     :: itemp(:) 
    real(r8), parameter  :: deg2rad  = SHR_CONST_PI/180._r8
    !-----------------------------------------------------------------------

    subgrid%name = dimname

    allocate(itemp(beg:end))
    allocate(subgrid%lat(beg:end), subgrid%lon(beg:end), subgrid%coslat(beg:end))
    do n = beg,end
       if (itemp(n) > 0) then
          active(n) = .true.
       else
          active(n) = .false.
       end if
       subgrid%lat(n) = subgrid%lat(n)*deg2rad
       subgrid%lon(n) = subgrid%lon(n)*deg2rad
       subgrid%coslat(n) = cos(subgrid%lat(n))
    end do

    deallocate(itemp)

  contains

    subroutine read_var_double(ncid, varname, data, dim1name, use_glob)
      ! Wraps the ncd_io call, providing logic related to whether we're using the 'glob'
      ! form of ncd_io
      type(file_desc_t)  , intent(inout) :: ncid
      character(len=*)   , intent(in)    :: varname
      real(r8), pointer  , intent(inout) :: data(:)
      character(len=*)   , intent(in)    :: dim1name
      logical            , intent(in)    :: use_glob  ! if .true., use the 'glob' form of ncd_io

      if (use_glob) then
         call ncd_io(ncid=ncid, varname=varname, flag='read', data=data)
      else
         call ncd_io(ncid=ncid, varname=varname, flag='read', data=data, dim1name=dim1name)
      end if
    end subroutine read_var_double

    subroutine read_var_int(ncid, varname, data, dim1name, use_glob)
      ! Wraps the ncd_io call, providing logic related to whether we're using the 'glob'
      ! form of ncd_io
      type(file_desc_t)  , intent(inout) :: ncid
      character(len=*)   , intent(in)    :: varname
      integer, pointer   , intent(inout) :: data(:)
      character(len=*)   , intent(in)    :: dim1name
      logical            , intent(in)    :: use_glob  ! if .true., use the 'glob' form of ncd_io

      if (use_glob) then
         call ncd_io(ncid=ncid, varname=varname, flag='read', data=data)
      else
         call ncd_io(ncid=ncid, varname=varname, flag='read', data=data, dim1name=dim1name)
      end if
    end subroutine read_var_int

  end subroutine set_subgrid_info

  !=======================================================================

  subroutine interp_0d_copy (varname, varname_i, xtype, ncidi, ncido)

    ! --------------------------------------------------------------------
    ! arguments
    character(len=*)   , intent(inout) :: varname   ! variable name on output file
    character(len=*)   , intent(in)    :: varname_i ! variable name on input file
    integer            , intent(in)    :: xtype
    type(file_desc_t)  , intent(inout) :: ncidi
    type(file_desc_t)  , intent(inout) :: ncido
    !
    ! local variables
    integer :: ivalue
    real(r8):: rvalue
    ! --------------------------------------------------------------------
 
    if (masterproc) then
       write(iulog,*) 'Copying      : ',trim(varname_i), ' => ', trim(varname)
    end if

    if (xtype == pio_int) then
       call ncd_io(ncid=ncidi, varname=trim(varname_i), flag='read' , data=ivalue)
       call ncd_io(ncid=ncido, varname=trim(varname), flag='write', data=ivalue)
    else if (xtype == pio_double) then
       call ncd_io(ncid=ncidi, varname=trim(varname_i), flag='read' , data=rvalue)
       call ncd_io(ncid=ncido, varname=trim(varname), flag='write', data=rvalue)
    else
       if (masterproc) then
          write(iulog,*)'ERROR interpinic: unhandled case for var ',trim(varname),' stopping' 
       end if
       call endrun(msg=errMsg(sourcefile, __LINE__))
    end if

  end subroutine interp_0d_copy

  !=======================================================================

  subroutine interp_1d_double (varname, varname_i, dimname, begi, endi, bego, endo, ncidi, ncido, &
       sgridindex)

    ! ------------------------ arguments ---------------------------------
    character(len=*)  , intent(inout) :: varname   ! variable name on output file
    character(len=*)  , intent(in)    :: varname_i ! variable name on input file
    character(len=*)  , intent(inout) :: dimname
    integer           , intent(in)    :: begi, endi
    integer           , intent(in)    :: bego, endo
    type(file_desc_t) , intent(inout) :: ncidi
    type(file_desc_t) , intent(inout) :: ncido
    integer           , intent(in)    :: sgridindex(bego:endo)
    ! --------------------------------------------------------------------

    ! ------------------------ local variables --------------------------
    real(r8), pointer :: rbufsli(:)     ! input array
    real(r8), pointer :: rbufslo(:)     ! output array
    ! --------------------------------------------------------------------

    if (masterproc) then
       write(iulog,*) 'Interpolating: ',trim(varname_i), ' => ', trim(varname)
    end if

    allocate (rbufsli(begi:endi), rbufslo(bego:endo))
    call ncd_io(ncid=ncidi, varname=trim(varname_i), flag='read', data=rbufsli)
    call ncd_io(ncid=ncido, varname=trim(varname), flag='read', data=rbufslo, &
         dim1name=dimname)

    call interp_1d_data(begi=begi, endi=endi, bego=bego, endo=endo, &
         sgridindex=sgridindex, keep_existing=.true., &
         data_in=rbufsli, data_out=rbufslo)

    call ncd_io(ncid=ncido, varname=trim(varname), flag='write', data=rbufslo, &
         dim1name=dimname)

    deallocate(rbufsli, rbufslo)

  end subroutine interp_1d_double

  !=======================================================================

  subroutine interp_1d_int (varname, varname_i, dimname, begi, endi, bego, endo, ncidi, ncido, &
       sgridindex)

    ! ------------------------ arguments ---------------------------------
    character(len=*)  , intent(inout) :: varname   ! variable name on output file
    character(len=*)  , intent(in)    :: varname_i ! variable name on input file
    character(len=*)  , intent(inout) :: dimname
    integer           , intent(in)    :: begi, endi
    integer           , intent(in)    :: bego, endo
    type(file_desc_t) , intent(inout) :: ncidi
    type(file_desc_t) , intent(inout) :: ncido
    integer           , intent(in)    :: sgridindex(bego:endo)
    ! --------------------------------------------------------------------

    ! ------------------------ local variables --------------------------
    integer , pointer :: ibufsli(:)     !input array
    integer , pointer :: ibufslo(:)     !output array
    ! --------------------------------------------------------------------

    if (masterproc) then
       write(iulog,*) 'Interpolating: ',trim(varname_i), ' => ', trim(varname)
    end if

    allocate (ibufsli(begi:endi), ibufslo(bego:endo))

    call ncd_io(ncid=ncidi, varname=trim(varname_i), flag='read', &
         data=ibufsli)
    call ncd_io(ncid=ncido, varname=trim(varname), flag='read', &
         data=ibufslo, dim1name=dimname)

    call interp_1d_data(begi=begi, endi=endi, bego=bego, endo=endo, &
         sgridindex=sgridindex, keep_existing=.true., &
         data_in=ibufsli, data_out=ibufslo)

    call ncd_io(ncid=ncido, varname=trim(varname), flag='write', &
         data=ibufslo, dim1name=dimname)

    deallocate (ibufsli, ibufslo)

  end subroutine interp_1d_int

  !=======================================================================

  subroutine interp_2d_double (var2di, var2do, &
       begi, endi, bego, endo, sgridindex, &
       interp_multilevel_container)

    ! --------------------------------------------------------------------
    ! arguments
    class(interp_2dvar_type), intent(inout) :: var2di  ! variable on input file
    class(interp_2dvar_type), intent(inout) :: var2do  ! variable on output file
    integer           , intent(in)    :: begi, endi
    integer           , intent(in)    :: bego, endo
    integer           , intent(in)    :: sgridindex(bego:)
    type(interp_multilevel_container_type), intent(in) :: interp_multilevel_container
    !
    ! local variables
    class(interp_multilevel_type), pointer :: multilevel_interpolator
    integer             :: no                   ! index
    integer             :: level                ! level index
    integer             :: nlevi                ! number of input levels
    real(r8), pointer   :: rbuf2do(:,:)         ! output array
    real(r8), pointer   :: rbuf1di(:)           ! one level of input array
    real(r8), pointer   :: rbuf2do_levelsi(:,:) ! array on output horiz grid, but input levels
    ! --------------------------------------------------------------------

    SHR_ASSERT_ALL((ubound(sgridindex) == (/endo/)), errMsg(sourcefile, __LINE__))
    SHR_ASSERT(var2di%get_vec_beg() == begi, errMsg(sourcefile, __LINE__))
    SHR_ASSERT(var2di%get_vec_end() == endi, errMsg(sourcefile, __LINE__))
    SHR_ASSERT(var2do%get_vec_beg() == bego, errMsg(sourcefile, __LINE__))
    SHR_ASSERT(var2do%get_vec_end() == endo, errMsg(sourcefile, __LINE__))

    multilevel_interpolator => interp_multilevel_container%find_interpolator( &
         lev_dimname = var2do%get_lev_dimname(), &
         vec_dimname = var2do%get_vec_dimname())

    if (masterproc) then
       write(iulog,*) 'Interpolating: ', &
            trim(var2di%get_varname()), ' => ', &
            trim(var2do%get_varname()), ': ', &
            multilevel_interpolator%get_description()
    end if

    call multilevel_interpolator%check_npts( &
         npts    = var2do%get_vec_npts(), &
         varname = var2do%get_varname())

    ! First do a horizontal interpolation. We need to separate the horizontal and vertical
    ! interpolation steps to avoid storing all levels of the source grid in memory at
    ! once: that is problematic in terms of memory requirements since the full source grid
    ! is stored in memory on every processor (in contrast to the destination grid, which
    ! is decomposed across processors).
    nlevi = var2di%get_nlev()
    allocate(rbuf2do_levelsi(bego:endo, nlevi))
    allocate(rbuf1di(begi:endi))
    do level = 1, nlevi
       ! COMPILER_BUG(wjs, 2015-11-25, cray8.4.0) The cray compiler has trouble
       ! resolving the generic reference here, giving the message: 'No specific
       ! match can be found for the generic subprogram call "READLEVEL"'. So we
       ! explicitly call the specific routine, rather than calling readlevel.
       call var2di%readlevel_double(rbuf1di, level)
       call interp_1d_data(begi=begi, endi=endi, bego=bego, endo=endo, &
            sgridindex=sgridindex, keep_existing=.false., &
            data_in=rbuf1di, data_out=rbuf2do_levelsi(:,level))
    end do
    deallocate(rbuf1di)

    ! Now do the vertical interpolation

    ! COMPILER_BUG(wjs, 2015-11-25, cray8.4.0) The cray compiler has trouble
    ! resolving the generic reference here, giving the message: 'No specific
    ! match can be found for the generic subprogram call "READVAR"'. So we
    ! explicitly call the specific routine, rather than calling readvar.
    call var2do%readvar_double(rbuf2do)
    do no = bego,endo
       ! Only do the interpolation for output points that have a corresponding input
       ! point. Other output points will remain at their original value.
       if (sgridindex(no) > 0) then
          call multilevel_interpolator%interp_multilevel( &
               data_dest    = rbuf2do(no,:), &
               data_source  = rbuf2do_levelsi(no,:), &
               index_dest   = no - bego + 1)
       end if
    end do

    ! COMPILER_BUG(wjs, 2015-11-25, cray8.4.0) The cray compiler has trouble
    ! resolving the generic reference here, giving the message: 'No specific
    ! match can be found for the generic subprogram call "WRITEVAR"'. So we
    ! explicitly call the specific routine, rather than calling writevar.
    call var2do%writevar_double(rbuf2do)
       
    deallocate(rbuf2do, rbuf2do_levelsi)

  end subroutine interp_2d_double

  !=======================================================================

  subroutine check_dim_subgrid(ncidi, ncido, dimname, dimleni, dimleno)

    ! --------------------------------------------------------------------
    ! arguments
    type(file_desc_t) , intent(inout) :: ncidi         
    type(file_desc_t) , intent(inout) :: ncido         
    character(len=*)  , intent(in)    :: dimname
    integer           , intent(out)   :: dimleni
    integer           , intent(out)   :: dimleno
    !
    ! local variables
    integer :: status
    integer :: dimid
    ! --------------------------------------------------------------------

    status = pio_inq_dimid (ncidi, dimname, dimid)
    status = pio_inq_dimlen(ncidi, dimid  , dimleni)
    status = pio_inq_dimid (ncido, dimname, dimid)
    status = pio_inq_dimlen(ncido, dimid  , dimleno)

  end subroutine check_dim_subgrid

  !=======================================================================

  subroutine check_dim_level(ncidi, ncido, dimname, must_be_same)
    ! Checks whether dimension size is the same for input and output. If 'must_be_same'
    ! is true, aborts if they disagree; otherwise, simply prints an informative message.

    ! --------------------------------------------------------------------
    ! arguments
    type(file_desc_t) , intent(inout) :: ncidi         
    type(file_desc_t) , intent(inout) :: ncido         
    character(len=*)  , intent(in)    :: dimname
    logical           , intent(in)    :: must_be_same
    !
    ! local variables
    integer :: status
    integer :: dimid
    integer :: dimleni, dimleno
    ! --------------------------------------------------------------------

    status = pio_inq_dimid (ncidi, dimname, dimid)
    status = pio_inq_dimlen(ncidi, dimid  , dimleni)
    status = pio_inq_dimid (ncido, dimname, dimid)
    status = pio_inq_dimlen(ncido, dimid  , dimleno)

    if (dimleni /= dimleno) then
       if (must_be_same) then
          write (iulog,*) 'ERROR interpinic: input and output ',trim(dimname),' values disagree'
          write (iulog,*) 'input dimlen = ',dimleni,' output dimlen = ',dimleno
          call endrun(msg=errMsg(sourcefile, __LINE__))
       else
          if (masterproc) then
             write (iulog,*) 'input and output ',trim(dimname),' values disagree'
             write (iulog,*) 'input nlevgrnd = ',dimleni,' output nlevgrnd = ',dimleno
             write (iulog,*) 'This is okay: vertical levels will be interpolated'
          end if
       end if
    end if

  end subroutine check_dim_level

end module initInterpMod
