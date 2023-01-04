module decompInitMod

  !------------------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Module provides a descomposition into a clumped data structure which can
  ! be mapped back to atmosphere physics chunks.
  !
  ! !USES:
  use shr_kind_mod    , only : r8 => shr_kind_r8
  use shr_sys_mod     , only : shr_sys_flush
  use shr_log_mod     , only : errMsg => shr_log_errMsg
  use spmdMod         , only : masterproc, iam, npes, mpicom, comp_id
  use abortutils      , only : endrun
  use clm_varctl      , only : iulog
  use clm_varcon      , only : grlnd
  use GridcellType    , only : grc
  use decompMod
  use mct_mod         , only : mct_gsMap_init, mct_gsMap_ngseg, mct_gsMap_nlseg, mct_gsmap_gsize
  !
  ! !PUBLIC TYPES:
  implicit none
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public decompInit_lnd    ! initializes lnd grid decomposition into clumps and processors
  public decompInit_clumps ! initializes atm grid decomposition into clumps
  public decompInit_glcp   ! initializes g,l,c,p decomp info
  !
  ! !PRIVATE TYPES:
  private
  integer, pointer :: lcid(:)       ! temporary for setting ldecomp

  character(len=*), parameter, private :: sourcefile = &
       __FILE__
  !------------------------------------------------------------------------------

contains

  !------------------------------------------------------------------------------
  subroutine decompInit_lnd(lni,lnj,amask)
    !
    ! !DESCRIPTION:
    ! This subroutine initializes the land surface decomposition into a clump
    ! data structure.  This assumes each pe has the same number of clumps
    ! set by clump_pproc
    !
    ! !USES:
    use clm_varctl, only : nsegspc
    !
    ! !ARGUMENTS:
    implicit none
    integer , intent(in) :: amask(:)
    integer , intent(in) :: lni,lnj   ! domain global size
    !
    ! !LOCAL VARIABLES:
    integer :: lns                    ! global domain size
    integer :: ln,lj                  ! indices
    integer :: ag,an,ai,aj            ! indices
    integer :: numg                   ! number of land gridcells
    logical :: seglen1                ! is segment length one
    real(r8):: seglen                 ! average segment length
    real(r8):: rcid                   ! real value of cid
    integer :: cid,pid                ! indices
    integer :: n,m,ng                 ! indices
    integer :: ier                    ! error code
    integer :: beg,end,lsize,gsize    ! used for gsmap init
    integer, pointer :: gindex(:)     ! global index for gsmap init
    integer, pointer :: clumpcnt(:)   ! clump index counter
    !------------------------------------------------------------------------------

    lns = lni * lnj

    !--- set and verify nclumps ---
    if (clump_pproc > 0) then
       nclumps = clump_pproc * npes
       if (nclumps < npes) then
          write(iulog,*) 'decompInit_lnd(): Number of gridcell clumps= ',nclumps, &
               ' is less than the number of processes = ', npes
          call endrun(msg=errMsg(sourcefile, __LINE__))
       end if
    else
       write(iulog,*)'clump_pproc= ',clump_pproc,'  must be greater than 0'
       call endrun(msg=errMsg(sourcefile, __LINE__))
    end if

    ! allocate and initialize procinfo and clumps 
    ! beg and end indices initialized for simple addition of cells later 

    allocate(procinfo%cid(clump_pproc), stat=ier)
    if (ier /= 0) then
       write(iulog,*) 'decompInit_lnd(): allocation error for procinfo%cid'
       call endrun(msg=errMsg(sourcefile, __LINE__))
    endif
    procinfo%nclumps   = clump_pproc
    procinfo%cid(:)    = -1
    procinfo%ncells    = 0
    procinfo%begg      = 1
    procinfo%endg      = 0

    allocate(clumps(nclumps), stat=ier)
    if (ier /= 0) then
       write(iulog,*) 'decompInit_lnd(): allocation error for clumps'
       call endrun(msg=errMsg(sourcefile, __LINE__))
    end if
    clumps(:)%owner     = -1
    clumps(:)%ncells    = 0
    clumps(:)%begg      = 1
    clumps(:)%endg      = 0

    ! assign clumps to proc round robin 
    cid = 0
    do n = 1,nclumps
       pid = mod(n-1,npes)
       if (pid < 0 .or. pid > npes-1) then
          write(iulog,*) 'decompInit_lnd(): round robin pid error ',n,pid,npes
          call endrun(msg=errMsg(sourcefile, __LINE__))
       endif
       clumps(n)%owner = pid
       if (iam == pid) then
          cid = cid + 1
          if (cid < 1 .or. cid > clump_pproc) then
             write(iulog,*) 'decompInit_lnd(): round robin pid error ',n,pid,npes
             call endrun(msg=errMsg(sourcefile, __LINE__))
          endif
          procinfo%cid(cid) = n
       endif
    enddo

    ! count total land gridcells
    numg = 0
    do ln = 1,lns
       if (amask(ln) == 1) then
          numg = numg + 1
       endif
    enddo

    if (npes > numg) then
       write(iulog,*) 'decompInit_lnd(): Number of processes exceeds number ', &
            'of land grid cells',npes,numg
       call endrun(msg=errMsg(sourcefile, __LINE__))
    end if
    if (nclumps > numg) then
       write(iulog,*) 'decompInit_lnd(): Number of clumps exceeds number ', &
            'of land grid cells',nclumps,numg
       call endrun(msg=errMsg(sourcefile, __LINE__))
    end if

    if (float(numg)/float(nclumps) < float(nsegspc)) then
       seglen1 = .true.
       seglen = 1.0_r8
    else
       seglen1 = .false.
       seglen = dble(numg)/(dble(nsegspc)*dble(nclumps))
    endif

    if (masterproc) then
       write(iulog,*) ' decomp precompute numg,nclumps,seglen1,avg_seglen,nsegspc=', &
            numg,nclumps,seglen1,&
            sngl(seglen),sngl(dble(numg)/(seglen*dble(nclumps)))
    end if

    ! Assign gridcells to clumps (and thus pes) ---

    allocate(lcid(lns))
    lcid(:) = 0
    ng = 0
    do ln = 1,lns
       if (amask(ln) == 1) then
          ng = ng  + 1

          !--- give to clumps in order based on nsegspc
          if (seglen1) then
             cid = mod(ng-1,nclumps) + 1
          else
             rcid = (dble(ng-1)/dble(numg))*dble(nsegspc)*dble(nclumps)
             cid = mod(int(rcid),nclumps) + 1
          endif
          lcid(ln) = cid

          !--- give gridcell cell to pe that owns cid ---
          !--- this needs to be done to subsequently use function
          !--- get_proc_bounds(begg,endg) 
          if (iam == clumps(cid)%owner) then
             procinfo%ncells  = procinfo%ncells  + 1
          endif
          if (iam >  clumps(cid)%owner) then
             procinfo%begg = procinfo%begg + 1
          endif
          if (iam >= clumps(cid)%owner) then
             procinfo%endg = procinfo%endg + 1
          endif

          !--- give gridcell to cid ---
          !--- increment the beg and end indices ---
          clumps(cid)%ncells  = clumps(cid)%ncells  + 1
          do m = 1,nclumps
             if ((clumps(m)%owner >  clumps(cid)%owner) .or. &
                 (clumps(m)%owner == clumps(cid)%owner .and. m > cid)) then
                clumps(m)%begg = clumps(m)%begg + 1
             endif
             
             if ((clumps(m)%owner >  clumps(cid)%owner) .or. &
                 (clumps(m)%owner == clumps(cid)%owner .and. m >= cid)) then
                clumps(m)%endg = clumps(m)%endg + 1
             endif
          enddo

       end if
    enddo

    ! Set ldecomp

    allocate(ldecomp%gdc2glo(numg), stat=ier)
    if (ier /= 0) then
       write(iulog,*) 'decompInit_lnd(): allocation error1 for ldecomp, etc'
       call endrun(msg=errMsg(sourcefile, __LINE__))
    end if
    allocate(clumpcnt(nclumps),stat=ier)
    if (ier /= 0) then
       write(iulog,*) 'decompInit_lnd(): allocation error1 for clumpcnt'
       call endrun(msg=errMsg(sourcefile, __LINE__))
    end if

    ldecomp%gdc2glo(:) = 0
    ag = 0

    ! clumpcnt is the start gdc index of each clump

    clumpcnt = 0
    ag = 1
    do pid = 0,npes-1
    do cid = 1,nclumps
       if (clumps(cid)%owner == pid) then
         clumpcnt(cid) = ag
         ag = ag + clumps(cid)%ncells
       endif
    enddo
    enddo

    ! now go through gridcells one at a time and increment clumpcnt
    ! in order to set gdc2glo

    do aj = 1,lnj
    do ai = 1,lni
       an = (aj-1)*lni + ai
       cid = lcid(an)
       if (cid > 0) then
          ag = clumpcnt(cid)
          ldecomp%gdc2glo(ag) = an
          clumpcnt(cid) = clumpcnt(cid) + 1
       end if
    end do
    end do

    deallocate(clumpcnt)

    ! Set gsMap_lnd_gdc2glo (the global index here includes mask=0 or ocean points)

    call get_proc_bounds(beg, end)
    allocate(gindex(beg:end))
    do n = beg,end
       gindex(n) = ldecomp%gdc2glo(n)
    enddo
    lsize = end-beg+1
    gsize = lni * lnj
    call mct_gsMap_init(gsMap_lnd_gdc2glo, gindex, mpicom, comp_id, lsize, gsize)
    deallocate(gindex)

    ! Diagnostic output

    if (masterproc) then
       write(iulog,*)' Surface Grid Characteristics'
       write(iulog,*)'   longitude points               = ',lni
       write(iulog,*)'   latitude points                = ',lnj
       write(iulog,*)'   total number of land gridcells = ',numg
       write(iulog,*)' Decomposition Characteristics'
       write(iulog,*)'   clumps per process             = ',clump_pproc
       write(iulog,*)' gsMap Characteristics'
       write(iulog,*) '  lnd gsmap glo num of segs      = ',mct_gsMap_ngseg(gsMap_lnd_gdc2glo)
       write(iulog,*)
    end if

    call shr_sys_flush(iulog)

  end subroutine decompInit_lnd

  !------------------------------------------------------------------------------
  subroutine decompInit_clumps(lns,lni,lnj)
    !
    ! !DESCRIPTION:
    ! This subroutine initializes the land surface decomposition into a clump
    ! data structure.  This assumes each pe has the same number of clumps
    ! set by clump_pproc
    !
    ! !USES:
    use spmdMod
    !
    ! !ARGUMENTS:
    implicit none
    integer , intent(in) :: lns,lni,lnj ! land domain global size
    !
    ! !LOCAL VARIABLES:
    integer :: ln,an              ! indices
    integer :: i,g,l,k            ! indices
    integer :: cid,pid            ! indices
    integer :: n,m,np             ! indices
    integer :: anumg              ! lnd num gridcells
    integer :: icells             ! temporary
    integer :: begg, endg         ! temporary
    integer :: ier                ! error code
    integer, allocatable :: allvecg(:,:)  ! temporary vector "global"
    integer, allocatable :: allvecl(:,:)  ! temporary vector "local"
    integer :: ntest
    character(len=32), parameter :: subname = 'decompInit_clumps'
    !------------------------------------------------------------------------------

    !--- assign gridcells to clumps (and thus pes) ---
    call get_proc_bounds(begg, endg)

    allocate(allvecl(nclumps,5))   ! local  clumps [gcells]
    allocate(allvecg(nclumps,5))   ! global clumps [gcells]

    ! Determine the number of gridcells
    ! on this processor 

    allvecg= 0
    allvecl= 0
    do anumg = begg,endg
       an  = ldecomp%gdc2glo(anumg)
       cid = lcid(an)
       ln  = anumg
       allvecl(cid,1) = allvecl(cid,1) + 1
    enddo
    call mpi_allreduce(allvecl,allvecg,size(allvecg),MPI_INTEGER,MPI_SUM,mpicom,ier)

    ! Determine overall total gridcells and distribute
    ! gridcells over clumps

    numg = 0

    do cid = 1,nclumps
       icells   = allvecg(cid,1)  ! number of all clump cid gridcells (over all processors)

       !--- overall total ---
       numg = numg + icells             ! total number of gridcells

       !--- give gridcell to the proc that owns the cid ---
       !--- increment the beg and end indices ---
    end do

    do n = 1,nclumps
       if (clumps(n)%ncells /= allvecg(n,1)) then
          write(iulog ,*) 'decompInit_glcp(): allvecg error ncells ',iam,n,clumps(n)%ncells   ,allvecg(n,1)
          call endrun(msg=errMsg(sourcefile, __LINE__))
       endif
    enddo

    deallocate(allvecg,allvecl)
    deallocate(lcid)

  end subroutine decompInit_clumps

  !------------------------------------------------------------------------------
  subroutine decompInit_glcp(lns,lni,lnj)
    !
    ! !DESCRIPTION:
    ! Determine gsMaps for landunits, columns, patches and cohorts
    !
    ! !USES:
    use spmdMod
    use spmdGathScatMod
    !
    ! !ARGUMENTS:
    implicit none
    integer , intent(in) :: lns,lni,lnj ! land domain global size
    !
    ! !LOCAL VARIABLES:
    integer :: gi,li,ci,pi,coi    ! indices
    integer :: i,g,k,l,n,np       ! indices
    integer :: cid,pid            ! indices
    integer :: begg,endg          ! beg,end gridcells
    integer :: numg               ! total number of gridcells across all processors
    integer :: icells             ! temporary
    integer :: ier                ! error code
    integer :: npmin,npmax,npint  ! do loop values for printing
    integer :: clmin,clmax        ! do loop values for printing
    integer :: locsize,globsize   ! used for gsMap init
    integer :: ng                 ! number of gridcells in gsMap_lnd_gdc2glo
    integer :: val1, val2         ! temporaries
    integer, pointer :: gindex(:) ! global index for gsMap init
    integer, pointer :: arrayglob(:) ! temporaroy
    integer, pointer :: gstart(:),  gcount(:)
    integer, pointer :: ioff(:)
    integer, parameter :: dbug=1      ! 0 = min, 1=normal, 2=much, 3=max
    character(len=32), parameter :: subname = 'decompInit_glcp'
    !------------------------------------------------------------------------------

    !init 

    call get_proc_bounds(begg, endg)
    call get_proc_global(ng=numg)

    ! Determine global seg megs

    allocate(gstart(begg:endg))
    gstart(:) = 0
    allocate(gcount(begg:endg))
    gcount(:) = 0
    allocate(ioff(begg:endg)) 
    ioff(:) = 0

    ! Determine gcount

    do gi = begg,endg
       gcount(gi)  = 1         ! number of gridcells for local gridcell index gi
    enddo

    ! Determine gstart

    ! gather the gdc subgrid counts to masterproc in glo order
    ! compute glo ordered start indices from the counts
    ! scatter the subgrid start indices back out to the gdc gridcells
    ! set the local gindex array for the subgrid from the subgrid start and count arrays

    ng = mct_gsmap_gsize(gsmap_lnd_gdc2glo)
    allocate(arrayglob(ng))

    arrayglob(:) = 0
    call gather_data_to_master(gcount, arrayglob, grlnd)
    if (masterproc) then
       val1 = arrayglob(1)
       arrayglob(1) = 1
       do n = 2,ng
          val2 = arrayglob(n)
          arrayglob(n) = arrayglob(n-1) + val1
          val1 = val2
       enddo
    endif
    call scatter_data_from_master(gstart, arrayglob, grlnd)

    deallocate(arrayglob)

    ! Gridcell gsmap (compressed, no ocean points)

    allocate(gindex(begg:endg))
    i = begg-1
    do gi = begg,endg
       if (gcount(gi) <  1) then
          write(iulog,*) 'decompInit_glcp warning count g ',k,iam,g,gcount(g)
       endif
       do l = 1,gcount(gi)
          i = i + 1
          if (i < begg .or. i > endg) then
             write(iulog,*) 'decompInit_glcp error i ',i,begg,endg
             call endrun(msg=errMsg(sourcefile, __LINE__))
          endif
          gindex(i) = gstart(gi) + l - 1
       enddo
    enddo
    if (i /= endg) then
       write(iulog,*) 'decompInit_glcp error size ',i,begg,endg
       call endrun(msg=errMsg(sourcefile, __LINE__))
    endif
    locsize = endg-begg+1
    globsize = numg
    call mct_gsMap_init(gsmap_gce_gdc2glo, gindex, mpicom, comp_id, locsize, globsize)
    deallocate(gindex)

    ! Deallocate start/count arrays
    deallocate(gstart, gcount)
    deallocate(ioff)

    ! Diagnostic output

    if (masterproc) then
       write(iulog,*)' Surface Grid Characteristics'
       write(iulog,*)'   longitude points          = ',lni
       write(iulog,*)'   latitude points           = ',lnj
       write(iulog,*)'   total number of gridcells = ',numg
       write(iulog,*)' Decomposition Characteristics'
       write(iulog,*)'   clumps per process        = ',clump_pproc
       write(iulog,*)' gsMap Characteristics'
       write(iulog,*) '  lnd gsmap glo num of segs = ',mct_gsMap_ngseg(gsMap_lnd_gdc2glo)
       write(iulog,*) '  gce gsmap glo num of segs = ',mct_gsMap_ngseg(gsMap_gce_gdc2glo)
       write(iulog,*)
    end if

    ! Write out clump and proc info, one pe at a time, 
    ! barrier to control pes overwriting each other on stdout

    call shr_sys_flush(iulog)
    call mpi_barrier(mpicom,ier)
    npmin = 0
    npmax = npes-1
    npint = 1
    if (dbug == 0) then
       npmax = 0
    elseif (dbug == 1) then
       npmax = min(npes-1,4)
    elseif (dbug == 2) then
       npint = npes/8
    endif
    do np = npmin,npmax,npint
       pid = np
       if (dbug == 1) then
          if (np == 2) pid=npes/2-1
          if (np == 3) pid=npes-2
          if (np == 4) pid=npes-1
       endif
       pid = max(pid,0)
       pid = min(pid,npes-1)

       if (iam == pid) then
          write(iulog,*)
          write(iulog,*)'proc= ',pid,&
               ' beg gridcell= ',procinfo%begg, &
               ' end gridcell= ',procinfo%endg,                   &
               ' total gridcells per proc= ',procinfo%ncells
          write(iulog,*)'proc= ',pid,&
               ' lnd ngseg   = ',mct_gsMap_ngseg(gsMap_lnd_gdc2glo), &
               ' lnd nlseg   = ',mct_gsMap_nlseg(gsMap_lnd_gdc2glo,iam)
          write(iulog,*)'proc= ',pid,&
               ' gce ngseg   = ',mct_gsMap_ngseg(gsMap_gce_gdc2glo), &
               ' gce nlseg   = ',mct_gsMap_nlseg(gsMap_gce_gdc2glo,iam)
          write(iulog,*)'proc= ',pid,' nclumps = ',procinfo%nclumps

          clmin = 1
          clmax = procinfo%nclumps
          if (dbug == 1) then
            clmax = 1
          elseif (dbug == 0) then
            clmax = -1
          endif
          do n = clmin,clmax
             cid = procinfo%cid(n)
             write(iulog,*)'proc= ',pid,' clump no = ',n, &
                  ' clump id= ',procinfo%cid(n),    &
                  ' beg gridcell= ',clumps(cid)%begg, &
                  ' end gridcell= ',clumps(cid)%endg, &
                  ' total gridcells per clump= ',clumps(cid)%ncells
          end do
       end if
       call shr_sys_flush(iulog)
       call mpi_barrier(mpicom,ier)
    end do
    call shr_sys_flush(iulog)

  end subroutine decompInit_glcp

end module decompInitMod
