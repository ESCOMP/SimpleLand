module lnd_import_export
  ! CTSM import and export fields exchanged with the coupler
  use ESMF                    , only : ESMF_GridComp, ESMF_State, ESMF_Mesh, ESMF_StateGet
  use ESMF                    , only : ESMF_KIND_R8, ESMF_SUCCESS, ESMF_MAXSTR, ESMF_LOGMSG_INFO
  use ESMF                    , only : ESMF_LogWrite, ESMF_LOGMSG_ERROR, ESMF_LogFoundError
  use ESMF                    , only : ESMF_STATEITEM_NOTFOUND, ESMF_StateItem_Flag
  use ESMF                    , only : operator(/=), operator(==)
  use NUOPC                   , only : NUOPC_CompAttributeGet, NUOPC_Advertise, NUOPC_IsConnected
  use NUOPC_Model             , only : NUOPC_ModelGet
  use shr_kind_mod            , only : r8 => shr_kind_r8, cx=>shr_kind_cx, cxx=>shr_kind_cxx, cs=>shr_kind_cs
  use shr_sys_mod             , only : shr_sys_abort
  use clm_varctl              , only : iulog
  use clm_time_manager        , only : get_nstep
  use decompmod               , only : bounds_type, get_proc_bounds
  use lnd2atmType             , only : lnd2atm_type
  use atm2lndType             , only : atm2lnd_type
  use domainMod               , only : ldomain
  use spmdMod                 , only : masterproc
  use shr_drydep_mod          , only : shr_drydep_readnl
  use shr_megan_mod           , only : shr_megan_readnl
  use nuopc_shr_methods       , only : chkerr
  use lnd_import_export_utils , only : check_for_errors, check_for_nans

  implicit none
  private ! except

  public  :: advertise_fields   ! Advertise the fields that can be sent/received
  public  :: realize_fields     ! Realize which fields will be sent and received
  public  :: import_fields      ! Import needed fields from mediator
  public  :: export_fields      ! Export fields from CTSM to mediator

  private :: fldlist_add
  private :: fldlist_realize
  private :: state_getimport_1d
  private :: state_getimport_2d
  private :: state_setexport_1d
  private :: state_setexport_2d
  private :: state_getfldptr
  private :: fldchk
  private :: ReadCapNamelist        ! Read in namelists governing import and export state

  type fld_list_type
     character(len=128) :: stdname
     integer :: ungridded_lbound = 0
     integer :: ungridded_ubound = 0
  end type fld_list_type

  integer, parameter     :: fldsMax = 100
  integer                :: fldsToLnd_num = 0
  integer                :: fldsFrLnd_num = 0
  type (fld_list_type)   :: fldsToLnd(fldsMax)
  type (fld_list_type)   :: fldsFrLnd(fldsMax)

  ! from lnd->atm
  character(len=cx)      :: carma_fields     ! List of CARMA fields from lnd->atm
  integer                :: drydep_nflds     ! number of dry deposition fields, ensure this is zero
  integer                :: megan_nflds      ! number of MEGAN voc fields from, ensure this is zero
  integer                :: emis_nflds       ! number of fire emission fields, ensure this is zero

  logical                :: flds_co2a        ! use case SLIM doesn't use CO2 so make sure isn't set
  logical                :: flds_co2b        ! use case (same as above)
  logical                :: flds_co2c        ! use case (same as above)
  logical                :: force_send_to_atm = .true.  ! Force sending export data to atmosphere even if ATM is not prognostic
  integer                :: glc_nec          ! number of glc elevation classes
  integer, parameter     :: debug = 0        ! internal debug level

  ! import fields
  character(*), parameter :: Sa_z                = 'Sa_z'
  character(*), parameter :: Sa_topo             = 'Sa_topo'
  character(*), parameter :: Sa_u                = 'Sa_u'
  character(*), parameter :: Sa_v                = 'Sa_v'
  character(*), parameter :: Sa_ptem             = 'Sa_ptem'
  character(*), parameter :: Sa_shum             = 'Sa_shum'
  character(*), parameter :: Sa_pbot             = 'Sa_pbot'
  character(*), parameter :: Sa_tbot             = 'Sa_tbot'
  character(*), parameter :: Faxa_rainc          = 'Faxa_rainc'
  character(*), parameter :: Faxa_rainl          = 'Faxa_rainl'
  character(*), parameter :: Faxa_snowc          = 'Faxa_snowc'
  character(*), parameter :: Faxa_snowl          = 'Faxa_snowl'
  character(*), parameter :: Faxa_lwdn           = 'Faxa_lwdn'
  character(*), parameter :: Faxa_swvdr          = 'Faxa_swvdr'
  character(*), parameter :: Faxa_swndr          = 'Faxa_swndr'
  character(*), parameter :: Faxa_swvdf          = 'Faxa_swvdf'
  character(*), parameter :: Faxa_swndf          = 'Faxa_swndf'

  ! export fields
  character(*), parameter :: Sl_lfrin       = 'Sl_lfrin'
  character(*), parameter :: Sl_t           = 'Sl_t'
  character(*), parameter :: Sl_snowh       = 'Sl_snowh'
  character(*), parameter :: Sl_avsdr       = 'Sl_avsdr'
  character(*), parameter :: Sl_anidr       = 'Sl_anidr'
  character(*), parameter :: Sl_avsdf       = 'Sl_avsdf'
  character(*), parameter :: Sl_anidf       = 'Sl_anidf'
  character(*), parameter :: Sl_tref        = 'Sl_tref'
  character(*), parameter :: Sl_qref        = 'Sl_qref'
  character(*), parameter :: Fall_taux      = 'Fall_taux'
  character(*), parameter :: Fall_tauy      = 'Fall_tauy'
  character(*), parameter :: Fall_lat       = 'Fall_lat'
  character(*), parameter :: Fall_sen       = 'Fall_sen'
  character(*), parameter :: Fall_lwup      = 'Fall_lwup'
  character(*), parameter :: Fall_evap      = 'Fall_evap'
  character(*), parameter :: Fall_swnet     = 'Fall_swnet'
  character(*), parameter :: Fall_flxdst    = 'Fall_flxdst'
  character(*), parameter :: Sl_u10         = 'Sl_u10'
  character(*), parameter :: Sl_ram1        = 'Sl_ram1'
  character(*), parameter :: Sl_fv          = 'Sl_fv'
  character(*), parameter :: Sl_soilw       = 'Sl_soilw'

  logical :: send_to_atm

  character(*),parameter :: F01 = "('(lnd_import_export) ',a,i5,2x,i5,2x,d21.14)"
  character(*),parameter :: u_FILE_u = &
       __FILE__

!===============================================================================
contains
!===============================================================================

  subroutine advertise_fields(gcomp, flds_scalar_name, atm_prognostic, rc)

    use shr_carma_mod     , only : shr_carma_readnl
    use shr_fire_emis_mod , only : shr_fire_emis_readnl
    use controlMod        , only : NLFilename

    ! input/output variables
    type(ESMF_GridComp)            :: gcomp
    character(len=*) , intent(in)  :: flds_scalar_name
    logical          , intent(in)  :: atm_prognostic
    integer          , intent(out) :: rc

    ! local variables
    type(ESMF_State)  :: importState
    type(ESMF_State)  :: exportState
    character(len=CS) :: cvalue
    logical           :: isPresent, isSet
    integer           :: n, num

    character(len=*), parameter :: subname='(lnd_import_export:advertise_fields)'
    !-------------------------------------------------------------------------------

    rc = ESMF_SUCCESS

    call NUOPC_ModelGet(gcomp, importState=importState, exportState=exportState, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! Determine  number of elevation classes
    call NUOPC_CompAttributeGet(gcomp, name='glc_nec', value=cvalue, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) glc_nec
    call ESMF_LogWrite('glc_nec = '// trim(cvalue), ESMF_LOGMSG_INFO)
    if (glc_nec < 1) then
       call shr_sys_abort('ERROR: In CLM4.5 and later, glc_nec must be at least 1.')
    end if

    !--------------------------------
    ! Advertise export fields
    !--------------------------------

    call ReadCapNamelist( NLFilename )

    ! Need to determine if there is no land for single column before the advertise call is done

    if (atm_prognostic .or. force_send_to_atm) then
       send_to_atm = .true.
    else
       send_to_atm = .false.
    end if

    if (send_to_atm) then
       call NUOPC_CompAttributeGet(gcomp, name='flds_co2a', value=cvalue, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       read(cvalue,*) flds_co2a
       call NUOPC_CompAttributeGet(gcomp, name='flds_co2b', value=cvalue, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       read(cvalue,*) flds_co2b
       call NUOPC_CompAttributeGet(gcomp, name='flds_co2c', value=cvalue, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       read(cvalue,*) flds_co2c
       if ( flds_co2a .or. flds_co2b .or. flds_co2c )then
          call shr_sys_abort(trim(subname)//": ERROR: SLIM can NOT pass CO2 fields to atmosphere")
       end if
    end if

    ! The following namelist reads should always be called regardless of the send_to_atm value

    ! (For dry-deposition, fire-emissions, and MEGAN just ensure they are NOT turned on)

    ! Dry Deposition velocities from land
    call shr_drydep_readnl("drv_flds_in", drydep_nflds)
    if ( drydep_nflds > 0 )then
       call shr_sys_abort(trim(subname)//": ERROR: dry deposition can NOT be turned on when running with SLIM")
    end if

    ! Fire emissions fluxes from land
    call shr_fire_emis_readnl('drv_flds_in', emis_nflds)
    if ( emis_nflds > 0 )then
       call shr_sys_abort(trim(subname)//": ERROR: fire emissions can NOT be turned on when running with SLIM")
    end if

    ! MEGAN VOC emissions fluxes from land
    call shr_megan_readnl('drv_flds_in', megan_nflds)
    if ( megan_nflds > 0 )then
       call shr_sys_abort(trim(subname)//": ERROR: MEGAN BVOC fields can NOT be turned on when running with SLIM")
    end if

    ! CARMA volumetric soil water from land
    ! TODO: is the following correct - the CARMA field exchange is very confusing in mct
    call shr_carma_readnl('drv_flds_in', carma_fields)

    ! export to atm
    call fldlist_add(fldsFrLnd_num, fldsFrlnd, trim(flds_scalar_name))
    call fldlist_add(fldsFrLnd_num, fldsFrlnd, 'Sl_lfrin')
    if (send_to_atm) then
       call fldlist_add(fldsFrLnd_num, fldsFrlnd, Sl_t          )
       call fldlist_add(fldsFrLnd_num, fldsFrlnd, Sl_tref       )
       call fldlist_add(fldsFrLnd_num, fldsFrlnd, Sl_qref       )
       call fldlist_add(fldsFrLnd_num, fldsFrlnd, Sl_avsdr      )
       call fldlist_add(fldsFrLnd_num, fldsFrlnd, Sl_anidr      )
       call fldlist_add(fldsFrLnd_num, fldsFrlnd, Sl_avsdf      )
       call fldlist_add(fldsFrLnd_num, fldsFrlnd, Sl_anidf      )
       call fldlist_add(fldsFrLnd_num, fldsFrlnd, Sl_snowh      )
       call fldlist_add(fldsFrLnd_num, fldsFrlnd, Sl_u10        )
       call fldlist_add(fldsFrLnd_num, fldsFrlnd, Sl_fv         )
       call fldlist_add(fldsFrLnd_num, fldsFrlnd, Sl_ram1       )
       call fldlist_add(fldsFrLnd_num, fldsFrlnd, Fall_taux     )
       call fldlist_add(fldsFrLnd_num, fldsFrlnd, Fall_tauy     )
       call fldlist_add(fldsFrLnd_num, fldsFrlnd, Fall_lat      )
       call fldlist_add(fldsFrLnd_num, fldsFrlnd, Fall_sen      )
       call fldlist_add(fldsFrLnd_num, fldsFrlnd, Fall_lwup     )
       call fldlist_add(fldsFrLnd_num, fldsFrlnd, Fall_evap     )
       call fldlist_add(fldsFrLnd_num, fldsFrlnd, Fall_swnet    )
       ! dust fluxes from land (4 sizes)
       call fldlist_add(fldsFrLnd_num, fldsFrLnd, Fall_flxdst, ungridded_lbound=1, ungridded_ubound=4)
       if (carma_fields /= ' ') then
          call fldlist_add(fldsFrLnd_num, fldsFrlnd, Sl_soilw) ! optional for carma
       end if
    end if

    call NUOPC_CompAttributeGet(gcomp, name="histaux_l2x1yrg", value=cvalue, &
         isPresent=isPresent, isSet=isSet, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! Now advertise above export fields
    do n = 1,fldsFrLnd_num
       call NUOPC_Advertise(exportState, standardName=fldsFrLnd(n)%stdname, &
            TransferOfferGeomObject='will provide', rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    enddo

    !--------------------------------
    ! Advertise import fields
    !--------------------------------

    call fldlist_add(fldsToLnd_num, fldsToLnd, trim(flds_scalar_name))

    ! from atm
    call fldlist_add(fldsToLnd_num, fldsToLnd, Sa_z         )
    call fldlist_add(fldsToLnd_num, fldsToLnd, Sa_topo      )
    call fldlist_add(fldsToLnd_num, fldsToLnd, Sa_u         )
    call fldlist_add(fldsToLnd_num, fldsToLnd, Sa_v         )
    call fldlist_add(fldsToLnd_num, fldsToLnd, Sa_ptem      )
    call fldlist_add(fldsToLnd_num, fldsToLnd, Sa_pbot      )
    call fldlist_add(fldsToLnd_num, fldsToLnd, Sa_tbot      )
    call fldlist_add(fldsToLnd_num, fldsToLnd, Sa_shum      )
    call fldlist_add(fldsToLnd_num, fldsToLnd, Faxa_lwdn    )
    call fldlist_add(fldsToLnd_num, fldsToLnd, Faxa_rainc   )
    call fldlist_add(fldsToLnd_num, fldsToLnd, Faxa_rainl   )
    call fldlist_add(fldsToLnd_num, fldsToLnd, Faxa_snowc   )
    call fldlist_add(fldsToLnd_num, fldsToLnd, Faxa_snowl   )
    call fldlist_add(fldsToLnd_num, fldsToLnd, Faxa_swndr   )
    call fldlist_add(fldsToLnd_num, fldsToLnd, Faxa_swvdr   )
    call fldlist_add(fldsToLnd_num, fldsToLnd, Faxa_swndf   )
    call fldlist_add(fldsToLnd_num, fldsToLnd, Faxa_swvdf   )

    ! Now advertise import fields
    do n = 1,fldsToLnd_num
       call NUOPC_Advertise(importState, standardName=fldsToLnd(n)%stdname, &
            TransferOfferGeomObject='will provide', rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    enddo

  end subroutine advertise_fields

  !===============================================================================
  subroutine realize_fields(importState, exportState, Emesh, flds_scalar_name, flds_scalar_num, rc)

    ! input/output variables
    type(ESMF_State) , intent(inout) :: importState
    type(ESMF_State) , intent(inout) :: exportState
    type(ESMF_Mesh)  , intent(in)    :: Emesh
    character(len=*) , intent(in)    :: flds_scalar_name
    integer          , intent(in)    :: flds_scalar_num
    integer          , intent(out)   :: rc

    ! local variables
    character(len=*), parameter :: subname='(lnd_import_export:realize_fields)'
    !---------------------------------------------------------------------------

    rc = ESMF_SUCCESS

    call fldlist_realize( &
         state=ExportState, &
         fldList=fldsFrLnd, &
         numflds=fldsFrLnd_num, &
         flds_scalar_name=flds_scalar_name, &
         flds_scalar_num=flds_scalar_num, &
         tag=subname//':clmExport',&
         mesh=Emesh, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call fldlist_realize( &
         state=importState, &
         fldList=fldsToLnd, &
         numflds=fldsToLnd_num, &
         flds_scalar_name=flds_scalar_name, &
         flds_scalar_num=flds_scalar_num, &
         tag=subname//':clmImport',&
         mesh=Emesh, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

  end subroutine realize_fields

  !===============================================================================
  subroutine import_fields( gcomp, bounds, &
       atm2lnd_inst, rc)

    !---------------------------------------------------------------------------
    ! Convert the input data from the mediator to the land model
    !---------------------------------------------------------------------------

    use clm_varcon              , only: rair, o2_molar_const, c13ratio
    use shr_const_mod           , only: SHR_CONST_TKFRZ
    use QSatMod                 , only: QSat
    use lnd_import_export_utils , only: derive_quantities, check_for_errors

    ! input/output variabes
    type(ESMF_GridComp)                         :: gcomp
    type(bounds_type)           , intent(in)    :: bounds         ! bounds
    type(atm2lnd_type)          , intent(inout) :: atm2lnd_inst   ! clm internal input data type
    integer                     , intent(out)   :: rc

    ! local variables
    type(ESMF_State)          :: importState
    type(ESMF_StateItem_Flag) :: itemFlag
    real(r8), pointer         :: dataPtr(:)
    real(r8), pointer         :: fldPtr1d(:)
    real(r8), pointer         :: fldPtr2d(:,:)
    character(len=CS)         :: fldname
    integer                   :: num
    integer                   :: begg, endg ! bounds
    integer                   :: g,i,k,n    ! indices
    real(r8)                  :: qsat_kg_kg ! saturation specific humidity (kg/kg)
    real(r8)                  :: forc_pbot  ! atmospheric pressure (Pa)
    real(r8)                  :: forc_rainc(bounds%begg:bounds%endg) ! rainxy Atm flux mm/s
    real(r8)                  :: forc_rainl(bounds%begg:bounds%endg) ! rainxy Atm flux mm/s
    real(r8)                  :: forc_snowc(bounds%begg:bounds%endg) ! snowfxy Atm flux  mm/s
    real(r8)                  :: forc_snowl(bounds%begg:bounds%endg) ! snowfxl Atm flux  mm/s
    real(r8)                  :: forc_noy(bounds%begg:bounds%endg)
    real(r8)                  :: forc_nhx(bounds%begg:bounds%endg)
    real(r8)                  :: frac_grc(bounds%begg:bounds%endg, 0:glc_nec)
    real(r8)                  :: topo_grc(bounds%begg:bounds%endg, 0:glc_nec)
    real(r8)                  :: hflx_grc(bounds%begg:bounds%endg, 0:glc_nec)
    real(r8)                  :: icemask_grc(bounds%begg:bounds%endg)
    real(r8)                  :: icemask_coupled_fluxes_grc(bounds%begg:bounds%endg)
    character(len=*), parameter :: subname='(lnd_import_export:import_fields)'
    !---------------------------------------------------------------------------

    rc = ESMF_SUCCESS
    call ESMF_LogWrite(subname//' called', ESMF_LOGMSG_INFO)

    ! Get import state
    call NUOPC_ModelGet(gcomp, importState=importState, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! Set bounds
    begg = bounds%begg; endg=bounds%endg

    ! Note: precipitation fluxes received  from the coupler
    ! are in units of kg/s/m^2. To convert these precipitation rates
    ! in units of mm/sec, one must divide by 1000 kg/m^3 and multiply
    ! by 1000 mm/m resulting in an overall factor of unity.
    ! Below the units are therefore given in mm/s.

    ! Required atm input fields
    call state_getimport_1d(importState, Sa_z      , atm2lnd_inst%forc_hgt_grc(begg:), rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_getimport_1d(importState, Sa_topo   , atm2lnd_inst%forc_topo_grc(begg:), rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_getimport_1d(importState, Sa_u      , atm2lnd_inst%forc_u_grc(begg:), rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_getimport_1d(importState, Sa_v      , atm2lnd_inst%forc_v_grc(begg:), rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_getimport_1d(importState, Sa_shum   , atm2lnd_inst%forc_q_not_downscaled_grc(begg:), rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_getimport_1d(importState, Sa_ptem   , atm2lnd_inst%forc_th_not_downscaled_grc(begg:), rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_getimport_1d(importState, Sa_pbot   , atm2lnd_inst%forc_pbot_not_downscaled_grc(begg:), rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_getimport_1d(importState, Sa_tbot   , atm2lnd_inst%forc_t_not_downscaled_grc(begg:), rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_getimport_1d(importState, Faxa_rainc, forc_rainc(begg:), rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_getimport_1d(importState, Faxa_rainl, forc_rainl(begg:), rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_getimport_1d(importState, Faxa_snowc, forc_snowc(begg:), rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_getimport_1d(importState, Faxa_snowl, forc_snowl(begg:), rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_getimport_1d(importState, Faxa_lwdn , atm2lnd_inst%forc_lwrad_not_downscaled_grc(begg:), rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_getimport_1d(importState, Faxa_swvdr, atm2lnd_inst%forc_solad_grc(begg:,1), rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_getimport_1d(importState, Faxa_swndr, atm2lnd_inst%forc_solad_grc(begg:,2), rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_getimport_1d(importState, Faxa_swvdf, atm2lnd_inst%forc_solai_grc(begg:,1), rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_getimport_1d(importState, Faxa_swndf, atm2lnd_inst%forc_solai_grc(begg:,2), rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    !--------------------------
    ! Derived quantities for required fields
    ! and corresponding error checks
    !--------------------------

    call derive_quantities(bounds, atm2lnd_inst, forc_rainc, forc_rainl, forc_snowc, forc_snowl)

    call check_for_errors(bounds, atm2lnd_inst )

  end subroutine import_fields

  !===============================================================================
  subroutine export_fields( gcomp, bounds, lnd2atm_inst, rc)

    !-------------------------------
    ! Pack the export state
    ! sign convention is positive downward with hierarchy of atm/glc/lnd/rof/ice/ocn.
    ! i.e. water sent from land to rof is positive
    !-------------------------------


    ! input/output variables
    type(ESMF_GridComp)                         :: gcomp
    type(bounds_type)           , intent(in)    :: bounds
    type(lnd2atm_type)          , intent(inout) :: lnd2atm_inst ! land to atmosphere exchange data type
    integer                     , intent(out)   :: rc

    ! local variables
    type(ESMF_State)  :: exportState
    real(r8), pointer :: fldPtr1d(:)
    real(r8), pointer :: fldPtr2d(:,:)
    character(len=CS) :: fldname
    integer           :: begg, endg
    integer           :: i, g, num
    real(r8)          :: data1d(bounds%begg:bounds%endg)
    character(len=*), parameter :: subname='(lnd_import_export:export_fields)'
    !---------------------------------------------------------------------------

    rc = ESMF_SUCCESS

    ! Get export state
    call NUOPC_ModelGet(gcomp, exportState=exportState, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! Set bounds
    begg = bounds%begg
    endg = bounds%endg

    ! -----------------------
    ! output to mediator
    ! -----------------------

    call state_setexport_1d(exportState, Sl_lfrin, ldomain%frac(begg:), init_spval=.false., rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! -----------------------
    ! output to atm
    ! -----------------------

    if (send_to_atm) then
       call state_setexport_1d(exportState, Sl_t      , lnd2atm_inst%t_rad_grc(begg:), &
            init_spval=.true., rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call state_setexport_1d(exportState, Sl_snowh  , lnd2atm_inst%h2osno_grc(begg:), &
            init_spval=.false., rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call state_setexport_1d(exportState, Sl_avsdr  , lnd2atm_inst%albd_grc(begg:,1), &
            init_spval=.true., rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call state_setexport_1d(exportState, Sl_anidr  , lnd2atm_inst%albd_grc(begg:,2), &
            init_spval=.true., rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call state_setexport_1d(exportState, Sl_avsdf  , lnd2atm_inst%albi_grc(begg:,1), &
            init_spval=.true., rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call state_setexport_1d(exportState, Sl_anidf  , lnd2atm_inst%albi_grc(begg:,2), &
            init_spval=.true., rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call state_setexport_1d(exportState, Sl_tref   , lnd2atm_inst%t_ref2m_grc(begg:), &
            init_spval=.true., rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call state_setexport_1d(exportState, Sl_qref   , lnd2atm_inst%q_ref2m_grc(begg:), &
            init_spval=.true., rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call state_setexport_1d(exportState, Fall_taux , lnd2atm_inst%taux_grc(begg:), &
            init_spval=.true., minus=.true., rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call state_setexport_1d(exportState, Fall_tauy , lnd2atm_inst%tauy_grc(begg:),  &
            init_spval=.true., minus=.true., rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call state_setexport_1d(exportState, Fall_lat  , lnd2atm_inst%eflx_lh_tot_grc(begg:), &
            init_spval=.true., minus=.true., rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call state_setexport_1d(exportState, Fall_sen  , lnd2atm_inst%eflx_sh_tot_grc(begg:), &
            init_spval=.true., minus=.true., rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call state_setexport_1d(exportState, Fall_lwup , lnd2atm_inst%eflx_lwrad_out_grc(begg:), &
            init_spval=.true., minus=.true., rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call state_setexport_1d(exportState, Fall_evap , lnd2atm_inst%qflx_evap_tot_grc(begg:), &
            init_spval=.true., minus=.true., rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call state_setexport_1d(exportState, Fall_swnet, lnd2atm_inst%fsa_grc(begg:), &
            init_spval=.true., rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

       ! optional fields
       call state_setexport_2d(exportState, Fall_flxdst, lnd2atm_inst%flxdst_grc(begg:,1:4), &
            init_spval=.true., minus= .true., rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call state_setexport_1d(exportState, Sl_u10, lnd2atm_inst%u_ref10m_grc(begg:), &
            init_spval=.false., rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call state_setexport_1d(exportState, Sl_ram1, lnd2atm_inst%ram1_grc(begg:), &
            init_spval=.true., rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call state_setexport_1d(exportState, Sl_fv, lnd2atm_inst%fv_grc(begg:), &
            init_spval=.true., rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       if (fldchk(exportState, Sl_soilw)) then
          call state_setexport_1d(exportState, Sl_soilw, lnd2atm_inst%h2osoi_vol_grc(begg:,1), &
               init_spval=.true., rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
       end if
    endif

  end subroutine export_fields

  !===============================================================================
  subroutine fldlist_add(num, fldlist, stdname, ungridded_lbound, ungridded_ubound)

    ! input/output variables
    integer,                    intent(inout) :: num
    type(fld_list_type),        intent(inout) :: fldlist(:)
    character(len=*),           intent(in)    :: stdname
    integer,          optional, intent(in)    :: ungridded_lbound
    integer,          optional, intent(in)    :: ungridded_ubound

    ! local variables
    integer :: rc
    character(len=*), parameter :: subname='(lnd_import_export:fldlist_add)'
    !-------------------------------------------------------------------------------

    ! Set up a list of field information

    num = num + 1
    if (num > fldsMax) then
       call ESMF_LogWrite(trim(subname)//": ERROR num > fldsMax "//trim(stdname), &
            ESMF_LOGMSG_ERROR, line=__LINE__, file=__FILE__)
       call shr_sys_abort(trim(subname)//": ERROR: num > fldsMax")
    endif
    fldlist(num)%stdname = trim(stdname)

    if (present(ungridded_lbound) .and. present(ungridded_ubound)) then
       fldlist(num)%ungridded_lbound = ungridded_lbound
       fldlist(num)%ungridded_ubound = ungridded_ubound
    end if

  end subroutine fldlist_add

  !===============================================================================
  subroutine fldlist_realize(state, fldList, numflds, flds_scalar_name, flds_scalar_num, mesh, tag, rc)

    use NUOPC , only : NUOPC_IsConnected, NUOPC_Realize
    use ESMF  , only : ESMF_MeshLoc_Element, ESMF_FieldCreate, ESMF_TYPEKIND_R8
    use ESMF  , only : ESMF_MAXSTR, ESMF_Field, ESMF_State, ESMF_Mesh, ESMF_StateRemove
    use ESMF  , only : ESMF_LogFoundError, ESMF_LOGMSG_INFO, ESMF_SUCCESS
    use ESMF  , only : ESMF_LogWrite, ESMF_LOGMSG_ERROR, ESMF_LOGERR_PASSTHRU

    ! input/output variables
    type(ESMF_State)    , intent(inout) :: state
    type(fld_list_type) , intent(in)    :: fldList(:)
    integer             , intent(in)    :: numflds
    character(len=*)    , intent(in)    :: flds_scalar_name
    integer             , intent(in)    :: flds_scalar_num
    character(len=*)    , intent(in)    :: tag
    type(ESMF_Mesh)     , intent(in)    :: mesh
    integer             , intent(inout) :: rc

    ! local variables
    integer                :: n
    type(ESMF_Field)       :: field
    character(len=80)      :: stdname
    character(len=*),parameter  :: subname='(lnd_import_export:fldlist_realize)'
    ! ----------------------------------------------

    rc = ESMF_SUCCESS

    do n = 1, numflds
       stdname = fldList(n)%stdname
       if (NUOPC_IsConnected(state, fieldName=stdname)) then
          if (stdname == trim(flds_scalar_name)) then
             call ESMF_LogWrite(trim(subname)//trim(tag)//" Field = "//trim(stdname)//" is connected on root pe", &
                  ESMF_LOGMSG_INFO)
             ! Create the scalar field
             call SetScalarField(field, flds_scalar_name, flds_scalar_num, rc=rc)
             if (ChkErr(rc,__LINE__,u_FILE_u)) return
          else
             ! Create the field
             if (fldlist(n)%ungridded_lbound > 0 .and. fldlist(n)%ungridded_ubound > 0) then
                field = ESMF_FieldCreate(mesh, ESMF_TYPEKIND_R8, name=stdname, meshloc=ESMF_MESHLOC_ELEMENT, &
                     ungriddedLbound=(/fldlist(n)%ungridded_lbound/), &
                     ungriddedUbound=(/fldlist(n)%ungridded_ubound/), &
                     gridToFieldMap=(/2/), rc=rc)
                if (ChkErr(rc,__LINE__,u_FILE_u)) return
             else
                field = ESMF_FieldCreate(mesh, ESMF_TYPEKIND_R8, name=stdname, meshloc=ESMF_MESHLOC_ELEMENT, rc=rc)
                if (ChkErr(rc,__LINE__,u_FILE_u)) return
             end if
             call ESMF_LogWrite(trim(subname)//trim(tag)//" Field = "//trim(stdname)//" is connected using mesh", &
                  ESMF_LOGMSG_INFO)
          endif

          ! NOW call NUOPC_Realize
          call NUOPC_Realize(state, field=field, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
       else
          if (stdname /= trim(flds_scalar_name)) then
             call ESMF_LogWrite(subname // trim(tag) // " Field = "// trim(stdname) // " is not connected.", &
                  ESMF_LOGMSG_INFO)
             call ESMF_StateRemove(state, (/stdname/), rc=rc)
             if (ChkErr(rc,__LINE__,u_FILE_u)) return
          end if
       end if
    end do

  contains  !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    subroutine SetScalarField(field, flds_scalar_name, flds_scalar_num, rc)
      ! ----------------------------------------------
      ! create a field with scalar data on the root pe
      ! ----------------------------------------------
      use ESMF, only : ESMF_Field, ESMF_DistGrid, ESMF_Grid
      use ESMF, only : ESMF_DistGridCreate, ESMF_GridCreate, ESMF_LogFoundError, ESMF_LOGERR_PASSTHRU
      use ESMF, only : ESMF_FieldCreate, ESMF_GridCreate, ESMF_TYPEKIND_R8

      type(ESMF_Field) , intent(inout) :: field
      character(len=*) , intent(in)    :: flds_scalar_name
      integer          , intent(in)    :: flds_scalar_num
      integer          , intent(inout) :: rc

      ! local variables
      type(ESMF_Distgrid) :: distgrid
      type(ESMF_Grid)     :: grid
      character(len=*), parameter :: subname='(lnd_import_export:SetScalarField)'
      ! ----------------------------------------------

      rc = ESMF_SUCCESS

      ! create a DistGrid with a single index space element, which gets mapped onto DE 0.
      distgrid = ESMF_DistGridCreate(minIndex=(/1/), maxIndex=(/1/), rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=u_FILE_u)) return

      grid = ESMF_GridCreate(distgrid, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=u_FILE_u)) return

      field = ESMF_FieldCreate(name=trim(flds_scalar_name), grid=grid, typekind=ESMF_TYPEKIND_R8, &
           ungriddedLBound=(/1/), ungriddedUBound=(/flds_scalar_num/), gridToFieldMap=(/2/), rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=u_FILE_u)) return

    end subroutine SetScalarField

  end subroutine fldlist_realize

  !===============================================================================
  subroutine state_getimport_1d(state, fldname, ctsmdata, rc)

    ! fill in ctsm import data for 1d field

    use ESMF, only : ESMF_LOGERR_PASSTHRU, ESMF_END_ABORT, ESMF_LogFoundError
    use ESMF, only : ESMF_Finalize

    ! input/output variabes
    type(ESMF_State) , intent(in)    :: state
    character(len=*) , intent(in)    :: fldname
    real(r8)         , intent(inout) :: ctsmdata(:)
    integer          , intent(out)   :: rc

    ! local variables
    real(r8), pointer :: fldPtr1d(:)
    integer           :: g
    character(len=*), parameter :: subname='(lnd_import_export:state_getimport_1d)'
    ! ----------------------------------------------

    rc = ESMF_SUCCESS

    call state_getfldptr(State, trim(fldname), fldptr1d=fldptr1d, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    do g = 1,size(ctsmdata)
       ctsmdata(g) = fldptr1d(g)
    end do
    call check_for_nans(ctsmdata, trim(fldname), 1)

  end subroutine state_getimport_1d

  !===============================================================================
  subroutine state_getimport_2d(state, fldname, ctsmdata, rc)

    ! fill in ctsm import data for 2d field

    use ESMF, only : ESMF_LOGERR_PASSTHRU, ESMF_END_ABORT, ESMF_LogFoundError
    use ESMF, only : ESMF_Finalize

    ! input/output variabes
    type(ESMF_State) , intent(in)    :: state
    character(len=*) , intent(in)    :: fldname
    real(r8)         , intent(inout) :: ctsmdata(:,:)
    integer          , intent(out)   :: rc

    ! local variables
    real(r8), pointer :: fldPtr2d(:,:)
    integer           :: g,n
    character(len=CS) :: cnum
    character(len=*), parameter :: subname='(lnd_import_export:state_getimport_1d)'
    ! ----------------------------------------------

    rc = ESMF_SUCCESS

    call state_getfldptr(state, trim(fldname), fldptr2d=fldptr2d, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    do n = 1,size(ctsmdata, dim=2)
       write(cnum,'(i0)') n
       do g = 1,size(ctsmdata,dim=1)
          ctsmdata(g,n) = fldptr2d(n,g)
       end do
       call check_for_nans(ctsmdata(:,n), trim(fldname)//trim(cnum), 1)
    end do

  end subroutine state_getimport_2d

  !===============================================================================
  subroutine state_setexport_1d(state, fldname, ctsmdata, init_spval, minus, rc)

    ! fill in ctsm export data for 1d field

    use ESMF          , only : ESMF_LOGERR_PASSTHRU, ESMF_END_ABORT, ESMF_LogFoundError
    use ESMF          , only : ESMF_Finalize
    use shr_const_mod , only : shr_const_spval

    ! input/output variabes
    type(ESMF_State) , intent(in) :: state
    character(len=*) , intent(in) :: fldname
    real(r8)         , intent(in) :: ctsmdata(:)
    logical          , intent(in) :: init_spval
    logical, optional, intent(in) :: minus
    integer          , intent(out):: rc

    ! local variables
    logical           :: l_minus ! local version of minus
    real(r8), pointer :: fldPtr1d(:)
    integer           :: g
    character(len=*), parameter :: subname='(lnd_export_export:state_setexport_1d)'
    ! ----------------------------------------------

    rc = ESMF_SUCCESS

    if (present(minus)) then
       l_minus = minus
    else
       l_minus = .false.
    end if

    call state_getfldptr(state, trim(fldname), fldptr1d=fldptr1d, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    if (init_spval) then
       fldptr1d(:) = shr_const_spval
    else
       fldptr1d(:) = 0._r8
    end if
    if (l_minus) then
       do g = 1,size(ctsmdata)
          fldptr1d(g) = -ctsmdata(g)
       end do
    else
       do g = 1,size(ctsmdata)
          fldptr1d(g) = ctsmdata(g)
       end do
    end if
    call check_for_nans(ctsmdata, trim(fldname), 1)

  end subroutine state_setexport_1d

  !===============================================================================
  subroutine state_setexport_2d(state, fldname, ctsmdata, init_spval, minus, rc)

    ! fill in ctsm export data for 2d field

    use ESMF          , only : ESMF_LOGERR_PASSTHRU, ESMF_END_ABORT, ESMF_LogFoundError
    use ESMF          , only : ESMF_Finalize
    use shr_const_mod , only : shr_const_spval

    ! input/output variabes
    type(ESMF_State) , intent(in) :: state
    character(len=*) , intent(in) :: fldname
    real(r8)         , intent(in) :: ctsmdata(:,:)
    logical,           intent(in) :: init_spval
    logical, optional, intent(in) :: minus
    integer          , intent(out):: rc

    ! local variables
    logical           :: l_minus         ! local version of minus
    real(r8), pointer :: fldPtr2d(:,:)
    integer           :: g, n
    character(len=CS) :: cnum
    character(len=*), parameter :: subname='(lnd_export_export:state_setexport_2d)'
    ! ----------------------------------------------

    rc = ESMF_SUCCESS

    if (present(minus)) then
       l_minus = minus
    else
       l_minus = .false.
    end if

    call state_getfldptr(state, trim(fldname), fldptr2d=fldptr2d, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    if (init_spval) then
       fldptr2d(:,:) = shr_const_spval
    else
       fldptr2d(:,:) = 0._r8
    end if
    do n = 1,size(ctsmdata, dim=2)
       write(cnum,'(i0)') n
       if (l_minus) then
          do g = 1,size(ctsmdata, dim=1)
             fldptr2d(n,g) = -ctsmdata(g,n)
          end do
       else
          do g = 1,size(ctsmdata, dim=1)
             fldptr2d(n,g) = ctsmdata(g,n)
          end do
       end if
       call check_for_nans(ctsmdata(:,n), trim(fldname)//trim(cnum), 1)
    end do

  end subroutine state_setexport_2d

  !===============================================================================
  subroutine state_getfldptr(State, fldname, fldptr1d, fldptr2d, rc)

    ! ----------------------------------------------
    ! Get pointer to a state field
    ! ----------------------------------------------

    use ESMF , only : ESMF_State, ESMF_Field, ESMF_Mesh, ESMF_FieldStatus_Flag
    use ESMF , only : ESMF_StateGet, ESMF_FieldGet, ESMF_MeshGet
    use ESMF , only : ESMF_FIELDSTATUS_COMPLETE, ESMF_FAILURE

    ! input/output variables
    type(ESMF_State),             intent(in)    :: State
    character(len=*),             intent(in)    :: fldname
    real(R8), pointer, optional , intent(out)   :: fldptr1d(:)
    real(R8), pointer, optional , intent(out)   :: fldptr2d(:,:)
    integer,                      intent(out)   :: rc

    ! local variables
    type(ESMF_FieldStatus_Flag) :: status
    type(ESMF_Field)            :: lfield
    character(len=*), parameter :: subname='(lnd_import_export:state_getfldptr)'
    ! ----------------------------------------------

    rc = ESMF_SUCCESS

    call ESMF_StateGet(State, itemName=trim(fldname), field=lfield, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    if (present(fldptr1d)) then
       call ESMF_FieldGet(lfield, farrayPtr=fldptr1d, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    else if (present(fldptr2d)) then
       call ESMF_FieldGet(lfield, farrayPtr=fldptr2d, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    else
       call shr_sys_abort("either fldptr1d or fldptr2d must be an input argument")
    end if

  end subroutine state_getfldptr

  !===============================================================================
  logical function fldchk(state, fldname)
    ! ----------------------------------------------
    ! Determine if field with fldname is in the input state
    ! ----------------------------------------------

    ! input/output variables
    type(ESMF_State), intent(in)  :: state
    character(len=*), intent(in)  :: fldname

    ! local variables
    type(ESMF_StateItem_Flag)   :: itemFlag
    ! ----------------------------------------------
    call ESMF_StateGet(state, trim(fldname), itemFlag)
    if (itemflag /= ESMF_STATEITEM_NOTFOUND) then
       fldchk = .true.
    else
       fldchk = .false.
    endif
  end function fldchk

  !===============================================================================
  subroutine ReadCapNamelist( NLFilename )

    ! ----------------------------------------------------
    ! Read in tne namelist for CTSM nuopc cap level items
    ! ----------------------------------------------------
    use clm_nlUtilsMod   , only : find_nlgroup_name
    use shr_mpi_mod      , only : shr_mpi_bcast
    use spmdMod          , only : mpicom
    use abortutils       , only : endrun
    use shr_log_mod      , only : errMsg => shr_log_errMsg
    ! !ARGUMENTS:
    character(len=*), intent(IN) :: NLFilename   ! Namelist filename
    ! !LOCAL VARIABLES:
    integer            :: nu_nml                           ! unit for namelist file
    integer            :: nml_error                        ! namelist i/o error flag
    character(*), parameter :: nml_name = "ctsm_nuopc_cap" ! MUST match with namelist name below
    namelist  /ctsm_nuopc_cap/ force_send_to_atm

    ! Read namelist
    if (masterproc) then
       open( newunit=nu_nml, file=trim(NLFilename), status='old', iostat=nml_error )
       call find_nlgroup_name(nu_nml, nml_name, status=nml_error)
       if (nml_error == 0) then
          read(nu_nml, nml=ctsm_nuopc_cap, iostat=nml_error)
          if (nml_error /= 0) then
             call endrun(msg='ERROR reading '//nml_name//' namelist'//errMsg(u_FILE_u, __LINE__))
          end if
       end if
       close(nu_nml)
    endif

    ! Broadcast namelist to all processors
    call shr_mpi_bcast(force_send_to_atm  , mpicom)

  end subroutine ReadCapNamelist

end module lnd_import_export
