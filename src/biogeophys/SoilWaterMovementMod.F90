module SoilWaterMovementMod

#include "shr_assert.h"

  !-----------------------------------------------------------------------
  ! DESCRIPTION
  ! module contains different subroutines to couple soil and root water interactions
  !
  ! created by Jinyun Tang, Mar 12, 2014
  use shr_log_mod    , only : errMsg => shr_log_errMsg
  use shr_kind_mod      , only : r8 => shr_kind_r8
  use shr_sys_mod         , only : shr_sys_flush
 
  !
  implicit none
  private
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: init_soilwater_movement
  !
  ! !PUBLIC DATA MEMBERS:

  ! !PRIVATE DATA MEMBERS:

  ! Solution method 
  integer, parameter :: zengdecker_2009 = 0
  integer, parameter :: moisture_form = 1
  integer, parameter :: mixed_form = 2
  integer, parameter :: head_form = 3

  ! Boundary conditions
  integer, parameter :: bc_head  = 0
  integer, parameter :: bc_flux  = 1
  integer, parameter :: bc_zero_flux  = 2
  integer, parameter :: bc_waterTable = 3
  integer, parameter :: bc_aquifer    = 4

  ! Soil hydraulic properties
  integer, parameter :: soil_hp_clapphornberg_1978=0
  integer, parameter :: soil_hp_vanGenuchten_1980=1

  real(r8),parameter :: m_to_mm = 1.e3_r8 !convert meters to mm

  integer :: soilwater_movement_method    ! method for solving richards equation
  integer :: upper_boundary_condition     ! named variable for the boundary condition
  integer :: lower_boundary_condition     ! named variable for the boundary condition

  ! Adaptive time stepping algorithmic control parameters
  real(r8) :: dtmin             ! minimum time step length (seconds)
  real(r8) :: verySmall         ! a very small number: used to check for sub step completion
  real(r8) :: xTolerUpper       ! tolerance to halve length of substep
  real(r8) :: xTolerLower       ! tolerance to double length of substep
  integer  :: expensive
  integer  :: inexpensive
  integer  :: flux_calculation

  character(len=*), parameter, private :: sourcefile = &
       __FILE__

  !-----------------------------------------------------------------------

contains

!#1
  !-----------------------------------------------------------------------
  subroutine init_soilwater_movement()
    !
    !DESCRIPTION
    !specify method for doing soil&root water interactions
    !
    ! !USES:
    use abortutils      , only : endrun   
    use fileutils       , only : getavu, relavu
    use spmdMod         , only : mpicom, masterproc
    use shr_mpi_mod     , only : shr_mpi_bcast
    use clm_varctl      , only : iulog
    use controlMod      , only : NLFilename
    use clm_nlUtilsMod  , only : find_nlgroup_name

    ! !ARGUMENTS:
    !------------------------------------------------------------------------------
    implicit none
    integer            :: nu_nml                     ! unit for namelist file
    integer            :: nml_error                  ! namelist i/o error flag
    character(*), parameter    :: subName = "('init_soilwater_movement')"

    !-----------------------------------------------------------------------

! MUST agree with name in namelist and read statement
    namelist /soilwater_movement_inparm/      &
         soilwater_movement_method,    &
         upper_boundary_condition,     &
         lower_boundary_condition,     &
         dtmin,                        &
         verySmall,                    &
         xTolerUpper,                  &
         xTolerLower,                  &
         expensive,                    &
         inexpensive,                  &
         flux_calculation

    ! Default values for namelist

    soilwater_movement_method = zengdecker_2009
    upper_boundary_condition = bc_flux
    lower_boundary_condition = bc_aquifer

    dtmin=60._r8          
    verySmall=1.e-8_r8    
    xTolerUpper=1.e-1_r8  
    xTolerLower=1.e-2_r8  
    expensive=42
    inexpensive=1
    flux_calculation=inexpensive  

    ! Read soilwater_movement namelist
    if (masterproc) then
       nu_nml = getavu()
       open( nu_nml, file=trim(NLFilename), status='old', iostat=nml_error )
       call find_nlgroup_name(nu_nml, 'soilwater_movement_inparm', status=nml_error)
       if (nml_error == 0) then
          read(nu_nml, nml=soilwater_movement_inparm,iostat=nml_error)
          if (nml_error /= 0) then
             call endrun(subname // ':: ERROR reading soilwater_movement namelist')
          end if
       else
          call endrun(subname // ':: ERROR reading soilwater_movement namelist')
       end if
       close(nu_nml)
       call relavu( nu_nml )

!  test for namelist consistency
       if((soilwater_movement_method == zengdecker_2009) .and. &
            (lower_boundary_condition /= bc_aquifer)) then
          call endrun(subname // ':: ERROR inconsistent soilwater_movement namelist: ZD09 must use bc_aquifer lbc')
       endif
    endif

    call shr_mpi_bcast(soilwater_movement_method, mpicom)
    call shr_mpi_bcast(upper_boundary_condition, mpicom)
    call shr_mpi_bcast(lower_boundary_condition, mpicom)
    call shr_mpi_bcast(dtmin, mpicom)
    call shr_mpi_bcast(verySmall, mpicom)
    call shr_mpi_bcast(xTolerUpper, mpicom)
    call shr_mpi_bcast(xTolerLower, mpicom)
    call shr_mpi_bcast(expensive, mpicom)
    call shr_mpi_bcast(inexpensive, mpicom)
    call shr_mpi_bcast(flux_calculation, mpicom)


    if (masterproc) then

       write(iulog,*) ' '
       write(iulog,*) 'soilwater_movement settings:'
       write(iulog,*) '  soilwater_movement_method  = ',soilwater_movement_method
       write(iulog,*) '  upper_boundary_condition   = ',upper_boundary_condition
       write(iulog,*) '  lower_boundary_condition   = ',lower_boundary_condition

       write(iulog,*) '  dtmin                      = ',dtmin
       write(iulog,*) '  verySmall                  = ',verySmall
       write(iulog,*) '  xTolerUpper                = ',xTolerUpper
       write(iulog,*) '  xTolerLower                = ',xTolerLower
       write(iulog,*) '  expensive                  = ',expensive
       write(iulog,*) '  inexpensive                = ',inexpensive
       write(iulog,*) '  flux_calculation           = ',flux_calculation
    endif

  end subroutine init_soilwater_movement
  

!#2
   !------------------------------------------------------------------------------   
   function use_aquifer_layer() result(lres)
     !
     !DESCRIPTION
     ! return true if an aquifer layer is used 
     ! otherwise false
     implicit none
     logical :: lres

     if(lower_boundary_condition == bc_aquifer .or. lower_boundary_condition == bc_watertable)then
        lres=.true.
     else
        lres=.false.
     endif
     return

   end function use_aquifer_layer

 end module SoilWaterMovementMod
