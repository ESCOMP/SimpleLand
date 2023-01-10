module initInterpMindist

  ! ------------------------------------------------------------------------
  ! This module contains most of the "interesting" logic of initInterp, in terms of
  ! finding the input column (or landunit, patch, etc.) to use as a template for each
  ! output column (etc.).
  !
  ! This is in a separate module to facilitate unit testing, since the full initInterp
  ! involves some awkward dependencies.
  ! ------------------------------------------------------------------------

  use shr_kind_mod   , only: r8 => shr_kind_r8
  use shr_log_mod    , only : errMsg => shr_log_errMsg
  use clm_varctl     , only: iulog
  use abortutils     , only: endrun
  use spmdMod        , only: masterproc
  use clm_varcon     , only: spval, re

  implicit none
  private
  save

  ! Public methods

  public :: set_mindist

  ! Public types

  type, public :: subgrid_type
     character(len=16) :: name               ! gridcell
     real(r8), pointer :: topoglc(:) => null()
     real(r8), pointer :: lat(:)
     real(r8), pointer :: lon(:)
     real(r8), pointer :: coslat(:)
  end type subgrid_type

  ! Private methods

  private :: is_sametype

  character(len=*), parameter, private :: sourcefile = &
       __FILE__

contains

  subroutine set_mindist(begi, endi, bego, endo, activei, activeo, subgridi, subgrido, &
       mindist_index)

    ! --------------------------------------------------------------------
    ! arguments
    integer            , intent(in)  :: begi, endi 
    integer            , intent(in)  :: bego, endo 
    logical            , intent(in)  :: activei(begi:endi) 
    logical            , intent(in)  :: activeo(bego:endo) 
    type(subgrid_type) , intent(in)  :: subgridi
    type(subgrid_type) , intent(in)  :: subgrido

    ! If false: if an output type cannot be found in the input, code aborts
    ! If true: if an output type cannot be found in the input, fill with closest natural
    ! veg column (using bare soil for patch-level variables)
    !
    ! NOTE: always treated as true for natural veg and crop landunits/columns/patches in
    ! the output - e.g., if we can't find the right column type to fill crop, we always
    ! use the closest natural veg column, regardless of the value of this flag.

    integer            , intent(out) :: mindist_index(bego:endo)
    !
    ! local variables
    real(r8) :: dx,dy
    real(r8) :: distmin,dist,hgtdiffmin,hgtdiff    
    integer  :: nsizei, nsizeo
    integer  :: ni,no,nmin,ier,n,noloc
    logical  :: closest
    ! --------------------------------------------------------------------

    mindist_index(bego:endo) = 0
    distmin = spval

!$OMP PARALLEL DO PRIVATE (ni,no,n,nmin,distmin,dx,dy,dist,closest,hgtdiffmin,hgtdiff)
    do no = bego,endo

       ! Only interpolate onto active points. Otherwise, the mere act of running
       ! init_interp (e.g., of a file onto itself) would lead to changes in a bunch of
       ! inactive points - e.g., going from their cold start initial conditions to some
       ! spunup initial conditions (from the closest active point of that type). This
       ! could potentially lead to different behavior in a transient run, if those points
       ! later became active; that's undesirable.
       if (activeo(no)) then 

          nmin    = 0
          distmin = spval
          hgtdiffmin = spval
          do ni = begi,endi
             if (activei(ni)) then
                if (is_sametype(ni, no, subgridi, subgrido)) then
                   dy = abs(subgrido%lat(no)-subgridi%lat(ni))*re
                   dx = abs(subgrido%lon(no)-subgridi%lon(ni))*re * &
                        0.5_r8*(subgrido%coslat(no)+subgridi%coslat(ni))
                   dist = dx*dx + dy*dy
                   if (associated(subgridi%topoglc) .and. associated(subgrido%topoglc)) then
                      hgtdiff = abs(subgridi%topoglc(ni) - subgrido%topoglc(no))
                   end if
                   closest = .false.
                   if ( dist < distmin ) then
                      closest = .true.
                      distmin = dist
                      nmin = ni
                      if (associated(subgridi%topoglc) .and. associated(subgrido%topoglc)) then
                         hgtdiffmin = hgtdiff
                      end if
                   end if
                   if (.not. closest) then
                      ! For glc_mec points, we first find the closest point in lat-lon
                      ! space (above). Then, within that closest point, we find the
                      ! closest column in topographic space; this second piece is done
                      ! here. Note that this ordering means that we could choose a column
                      ! with a very different topographic height from the target column,
                      ! if it is closer in lat-lon space than any similar-height columns.
                      if (associated(subgridi%topoglc) .and. associated(subgrido%topoglc)) then
                         hgtdiff = abs(subgridi%topoglc(ni) - subgrido%topoglc(no))
                         if ((dist == distmin) .and. (hgtdiff < hgtdiffmin)) then
                            closest = .true.
                            hgtdiffmin = hgtdiff
                            distmin = dist
                            nmin = ni
                         end if
                      end if
                   end if
                end if
             end if
          end do
          
          mindist_index(no) = nmin

       end if ! end if activeo block
    end do
!$OMP END PARALLEL DO
    
  end subroutine set_mindist

  !=======================================================================

  logical function is_sametype (ni, no, subgridi, subgrido)

    ! --------------------------------------------------------------------
    ! arguments
    integer           , intent(in)  :: ni 
    integer           , intent(in)  :: no 
    type(subgrid_type), intent(in)  :: subgridi
    type(subgrid_type), intent(in)  :: subgrido
    ! --------------------------------------------------------------------

    is_sametype = .false.

    if (trim(subgridi%name) == 'gridcell' .and. trim(subgrido%name) == 'gridcell') then
       is_sametype = .true.
    else 
       if (masterproc) then
          write(iulog,*)'ERROR interpinic: is_sametype check on input and output type not supported'
          write(iulog,*)'typei = ',trim(subgridi%name)
          write(iulog,*)'typeo = ',trim(subgrido%name)
       end if
       call endrun(msg=errMsg(sourcefile, __LINE__))
    end if

  end function is_sametype

end module initInterpMindist
