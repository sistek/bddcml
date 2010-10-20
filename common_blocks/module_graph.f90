module module_graph
! module for operations on graphs
! Jakub Sistek, Praha, 2009
implicit none 

contains 

subroutine graph_get_dual_mesh(nelem,nnod,inet,linet,nnet,lnnet,&
                               netn,lnetn,ietn,lietn,kietn,lkietn)
! get dual mesh - array of elements at nodes 

implicit none
! number of elements in mesh
integer, intent(in) :: nelem
! number of nodes in mesh
integer, intent(in) :: nnod

! PMD mesh description
integer, intent(in) :: linet
integer, intent(in) ::  inet(linet)
integer, intent(in) :: lnnet
integer, intent(in) ::  nnet(lnnet)

! dual mesh desctription 
integer, intent(in) :: lnetn
integer, intent(out) :: netn(lnetn)
integer, intent(in) :: lietn
integer, intent(out) :: ietn(lietn)
integer, intent(in) :: lkietn
integer, intent(out) :: kietn(lkietn)

! local variables
integer :: ie, indinet, indnode, ine, nne, indelmn, pointietn, inod

! Count elements touching particular node
! zero whole field of counts
      netn = 0
      indinet = 0
      do ie = 1,nelem
         nne = nnet(ie)
         do ine = 1,nne
            indinet = indinet + 1
            indnode = inet(indinet)
            netn(indnode) = netn(indnode) + 1
         end do
      end do

      ! Create array of pointers into ietn (transpose of inet)
      kietn(1) = 0
      do inod = 2,nnod
         kietn(inod) = kietn(inod-1) + netn(inod-1)
      end do

! zero whole field of counts and use it as secondary pointers
      netn = 0
      ietn = 0

      indinet = 0
      do ie = 1,nelem
         nne = nnet(ie)
         do ine = 1,nne
            indinet = indinet + 1
            indnode = inet(indinet)
            netn(indnode) = netn(indnode) + 1
            indelmn   = netn(indnode)
            pointietn = kietn(indnode)
            ietn(pointietn + indelmn) = ie
         end do
      end do
end subroutine graph_get_dual_mesh

subroutine graph_from_mesh_size(nelem,graphtype,neighbouring,inet,linet,nnet,lnnet,ietn,lietn,netn,lnetn,kietn,lkietn,&
                                xadj,lxadj, nedge, ladjncy, ladjwgt)
! find size for graph from mesh in PMD format 
use module_utils
implicit none
! number of elements in mesh
integer, intent(in) :: nelem
! type of output graph
integer, intent(in) :: graphtype
! prescribed value of number of shared nodes between two neighbours
integer, intent(in) :: neighbouring

! PMD mesh description
integer, intent(in) :: linet
integer, intent(in) ::  inet(linet)
integer, intent(in) :: lnnet
integer, intent(in) ::  nnet(lnnet)
! PMD dual mesh description
integer, intent(in) :: lietn
integer, intent(in) ::  ietn(lietn)
integer, intent(in) :: lnetn
integer, intent(in) ::  netn(lnetn)
integer, intent(in) :: lkietn
integer, intent(in) ::  kietn(lkietn)

! METIS graph description
integer, intent(in)  :: lxadj
integer, intent(out) ::  xadj(lxadj)
integer, intent(out) ::  nedge
integer, intent(out) ::  ladjncy
integer, intent(out) ::  ladjwgt

! local variables
integer,allocatable :: onerow(:), onerowweig(:)
integer :: nnetx, netnx, lonerow, lonerowweig, ie, indinet, indnode, ine, ionerow, nelmn, nne, pointietn, lorin, lorout

! prepare arrays for storing a row
      nnetx = maxval(nnet)
      netnx = maxval(netn)
      lonerow = nnetx * netnx
      lonerowweig = lonerow
      allocate(onerow(lonerow),onerowweig(lonerowweig))
      
      ! zero whole array
      xadj(1) = 1

      indinet = 0
      do ie = 1,nelem
         nne = nnet(ie)
         ! zero local row
         onerow     = 0
         onerowweig = 0
         ionerow = 0
         do ine = 1,nne
            indinet = indinet + 1
            indnode = inet(indinet)
            nelmn   = netn(indnode)
            pointietn = kietn(indnode)
            onerow(ionerow+1:ionerow + nelmn) = ietn(pointietn+1:pointietn+nelmn)
            ionerow = ionerow + nelmn
         end do
         lorin = ionerow
         ! parse onerow
         call graph_parse_onerow(ie,neighbouring,onerow,onerowweig,lorin,lorout)
         xadj(ie + 1) = xadj(ie) + count(onerowweig.ge.neighbouring)
      end do

      ladjncy = xadj(nelem+1)-1
      ! check the graph
      if (mod(ladjncy,2).ne.0) then
         write(*,*) 'GRAPH_FROM_MESH_SIZE: Number of nodes has to be even number!'
         call error_exit
      end if
      nedge = ladjncy / 2
      if (graphtype.eq.1) then
         ladjwgt = ladjncy
      else
         ladjwgt = 0
      end if

      deallocate(onerow,onerowweig)
end subroutine graph_from_mesh_size

subroutine graph_from_mesh(nelem,graphtype,neighbouring,inet,linet,nnet,lnnet,ietn,lietn,netn,lnetn,kietn,lkietn,&
                           xadj,lxadj, adjncy,ladjncy, adjwgt,ladjwgt)
! find size for graph from mesh in PMD format 
use module_utils
implicit none
! number of elements in mesh
integer, intent(in) :: nelem
! type of output graph
integer, intent(in) :: graphtype
! prescribed value of number of shared nodes between two neighbours
integer, intent(in) :: neighbouring

! PMD mesh description
integer, intent(in) :: linet
integer, intent(in) ::  inet(linet)
integer, intent(in) :: lnnet
integer, intent(in) ::  nnet(lnnet)
! PMD dual mesh description
integer, intent(in) :: lietn
integer, intent(in) ::  ietn(lietn)
integer, intent(in) :: lnetn
integer, intent(in) ::  netn(lnetn)
integer, intent(in) :: lkietn
integer, intent(in) ::  kietn(lkietn)

! METIS graph description
integer, intent(in) :: lxadj
integer, intent(in) ::  xadj(lxadj)
integer, intent(in) :: ladjncy
integer, intent(out) :: adjncy(ladjncy)
integer, intent(in) :: ladjwgt
integer, intent(out) :: adjwgt(ladjwgt)

! local variables
integer,allocatable :: onerow(:), onerowweig(:)
integer :: nnetx, netnx, lonerow, lonerowweig, ie, indinet, indnode, ine, ionerow, nelmn, nne, pointietn, lorin, lorout

! prepare arrays for storing a row
      nnetx = maxval(nnet)
      netnx = maxval(netn)
      lonerow = nnetx * netnx
      lonerowweig = lonerow
      allocate(onerow(lonerow),onerowweig(lonerowweig))
      
      indinet = 0
      do ie = 1,nelem
         onerow     = 0
         onerowweig = 0
         ionerow = 0
         nne = nnet(ie)
         do ine = 1,nne
            indinet = indinet + 1
            indnode = inet(indinet)
            nelmn = netn(indnode)
            pointietn = kietn(indnode)
            onerow(ionerow+1:ionerow + nelmn) = ietn(pointietn+1:pointietn+nelmn)
            ionerow = ionerow + nelmn
         end do
         lorin = ionerow
         ! parse onerow
         call graph_parse_onerow(ie,neighbouring,onerow,onerowweig,lorin,lorout)
         ! now only adjacencies above the level considered
         adjncy(xadj(ie):xadj(ie)+lorout-1) = onerow(1:lorout)
         if (graphtype.eq.1) then
            adjwgt(xadj(ie):xadj(ie+1)-1) = onerowweig(1:lorout)
         end if
      end do

      deallocate(onerow,onerowweig)
end subroutine graph_from_mesh

subroutine graph_check(nvertex,graphtype, xadj,lxadj, adjncy,ladjncy, adjwgt,ladjwgt)
! check the graph
use module_utils
implicit none
! number of vertices in graph
integer, intent(in) :: nvertex
! type of output graph
integer, intent(in) :: graphtype

! METIS graph description
integer, intent(in) :: lxadj
integer, intent(in) ::  xadj(lxadj)
integer, intent(in) :: ladjncy
integer, intent(out) :: adjncy(ladjncy)
integer, intent(in) :: ladjwgt
integer, intent(out) :: adjwgt(ladjwgt)

! local variables
integer :: iadjelm, indadjelm, iadjelmadj, indadjelmadj, ivertex, nadjelm, nadjelmadj
logical :: match

      ! Check the graph
      do ivertex = 1,nvertex
         nadjelm = xadj(ivertex+1) - xadj(ivertex)
         do iadjelm = 1,nadjelm
            indadjelm = adjncy(xadj(ivertex)-1 + iadjelm)
            match = .false. 
            indadjelmadj = 0
            nadjelmadj = xadj(indadjelm+1) - xadj(indadjelm)
            do iadjelmadj = 1,nadjelmadj
               if (adjncy(xadj(indadjelm)-1 + iadjelmadj).eq.ivertex) then
                  if (match) then
                     write(*,*) 'GRAPH_CHECK: Element ',ivertex,' multiply mentioned in the list of neighbours of element', &
                                indadjelm,'!'
                     call error_exit
                  end if
                  match = .true.
                  indadjelmadj = iadjelmadj
               end if
            end do
            if (.not.match) then
               write(*,*) 'GRAPH_CHECK: No match! Couple ',ivertex,indadjelm,' mentioned but not the opposite!'
               call error_exit
            end if
            if (graphtype.eq.1) then
               if (adjwgt(xadj(indadjelm)-1+indadjelmadj).ne.adjwgt(xadj(ivertex)-1+iadjelm)) then
                  write(*,*) 'GRAPH_CHECK: Non-matching edge weights between elements ', ivertex, indadjelm,'!'
!                  write(*,*) 'Indices of adjacent elements:'
!                  do ie2 = 1,nelem
!                     nadjelm = naetet(ie2)
!                     write(*,'(100i7)') ie2, iaetet(ie2,1:nadjelm)
!                  end do
!                  write(*,*) 'Weights of elements:'
!                  do ie2 = 1,nelem
!                     nadjelm = naetet(ie2)
!                     write(*,'(100i7)') ie2, edgeweights(ie2,1:nadjelm)
!                  end do
                  call error_exit
               end if
            end if
         end do
      end do
end subroutine graph_check

subroutine graph_write_to_file(idfile,nvertex,nedge,graphtype, xadj,lxadj, adjncy,ladjncy, adjwgt,ladjwgt)
! Write the list into file for METIS

use module_utils
implicit none
! unit for writing the graph
integer, intent(in) :: idfile
! number of vertices in graph
integer, intent(in) :: nvertex
! number of edges in graph
integer, intent(in) :: nedge
! type of output graph
integer, intent(in) :: graphtype

! METIS graph description
integer, intent(in) :: lxadj
integer, intent(in) ::  xadj(lxadj)
integer, intent(in) :: ladjncy
integer, intent(out) :: adjncy(ladjncy)
integer, intent(in) :: ladjwgt
integer, intent(out) :: adjwgt(ladjwgt)

! local variables
integer :: ie, nadje, j

      write(idfile,'(x,i10,x,i10,x,i10)') nvertex, nedge, graphtype
      do ie = 1,nvertex
         nadje = xadj(ie+1) - xadj(ie)
         if      (graphtype.eq.0) then
            ! no weights 
            write(idfile,'(600i9)') nadje, (adjncy(xadj(ie)-1+j), j = 1,nadje)
         else if (graphtype.eq.1) then
            ! weighted graph 
            write(idfile,'(600i9)') nadje, (adjncy(xadj(ie)-1+j), adjwgt(xadj(ie)-1+j), j = 1,nadje)
         else
            write(*,*) 'GRAPH_WRITE_TO_FILE: Graph type not supported: ',graphtype
            call error_exit
         end if
      end do
end subroutine graph_write_to_file

subroutine graph_parse_onerow(ie,neighbouring,onerow,onerowweig,lorin,lorout)

! contains quicksort routine for integers
use module_utils

implicit none

integer,intent(in) :: ie
integer,intent(in) :: neighbouring
integer,intent(in) ::    lorin
integer,intent(inout) :: onerow(lorin)
integer,intent(inout) :: onerowweig(lorin)
integer,intent(out) :: lorout

! local variables
integer :: valid, ivalid, nvalid, i, i2

         ! eliminate myself
         where(onerow .eq. ie) onerow = 0

         ! sort elements
         call iquick_sort(onerow,lorin)

         ! trim initial zeros
         !  - find where nonzeros start
         do i = 1,lorin
            if (onerow(i) .ne. 0) then
               i2 = i
               exit
            end if
         end do
         ! shift array
         onerow = cshift(onerow,i2-1)

         ! remove multiplicities
         ivalid = 1
         valid  = onerow(1)
         nvalid = 1
         do i = 2,lorin
            if (onerow(i).gt.valid) then
               onerowweig(ivalid) = nvalid
               nvalid = 1 ! reset count
               ivalid = ivalid + 1
               valid  = onerow(i)
               onerow(ivalid) = valid
            else if (onerow(i).lt.valid) then
               onerowweig(ivalid) = nvalid
               exit
            else               
               nvalid = nvalid + 1 ! add count
            end if
         end do
         do i = ivalid+1,lorin
            onerow(i)     = 0
            onerowweig(i) = 0
         end do

         lorout = ivalid
         ! no repeating indices in onerow

         ! eliminate down limit for adjacencies
         where(onerowweig(1:lorout).lt.neighbouring) onerow(1:lorout)     = 0
         where(onerowweig(1:lorout).lt.neighbouring) onerowweig(1:lorout) = 0

         lorout = count(onerow .ne. 0)
         onerow     = pack(onerow,onerow.ne.0)
         onerowweig = pack(onerowweig,onerowweig.ne.0)
         onerow(lorout+1:)     = 0
         onerowweig(lorout+1:) = 0

end subroutine graph_parse_onerow

subroutine graph_parse_onerow2(ie,neighbouring,onerow,onerowweig,lor)
! older (and slower) version of the routine
implicit none

integer,intent(in) :: ie
integer,intent(in) :: neighbouring
integer,intent(inout) :: onerow(:)
integer,intent(inout) :: onerowweig(:)
integer,intent(inout) :: lor

! local variables
integer :: io, indel

         ! eliminate myself
         where(onerow .eq. ie) onerow = 0
         ! eliminate multiplicities
         lor = count(onerow .ne. 0)
         onerow = pack(onerow,onerow.ne.0)
         onerow(lor+1:) = 0
         io = 0
         do 
            io = io + 1
            indel = onerow(io)
            if (indel.eq.0) exit
            onerowweig(io) = count(onerow.eq.indel)
            where (onerow(io+1:).eq.indel) onerow(io+1:) = 0
            lor = count(onerow .ne. 0)
            onerow = pack(onerow,onerow.ne.0)
            onerow(lor+1:) = 0
         end do
         ! no repeating indices in onerow
         ! eliminate down limit for adjacencies
         where(onerowweig(1:lor).lt.neighbouring) onerow(1:lor)     = 0
         where(onerowweig(1:lor).lt.neighbouring) onerowweig(1:lor) = 0
         lor = count(onerow .ne. 0)
         onerow     = pack(onerow,onerow.ne.0)
         onerowweig = pack(onerowweig,onerowweig.ne.0)
         onerow(lor+1:)     = 0
         onerowweig(lor+1:) = 0

end subroutine graph_parse_onerow2

subroutine graph_divide(graphtype,nvertex,xadj,lxadj,adjncy,ladjncy,vwgt,lvwgt,adjwgt,ladjwgt,nsub,edgecut,part,lpart)
! Divide graph by METIS
use module_utils
implicit none
! type of output graph
integer, intent(in) :: graphtype
integer, intent(in) :: nvertex
integer, intent(in) ::   lxadj
integer*4, intent(in) ::  xadj(lxadj)
integer, intent(in) ::   ladjncy
integer*4, intent(in) ::  adjncy(ladjncy)
integer, intent(in) ::   lvwgt
integer*4, intent(in) ::  vwgt(lvwgt)
integer, intent(in) ::   ladjwgt
integer*4, intent(in) ::  adjwgt(ladjwgt)
integer, intent(in) ::   nsub
integer, intent(out) ::  edgecut
integer, intent(in) ::    lpart
integer*4, intent(out) ::  part(lpart)

! local vars
integer*4 :: wgtflag 
integer*4 :: npart 
integer*4,parameter :: numflag = 1 !( 1 - Fortran-like arrays (0 - C-like arrays) 
integer*4 :: ec 
integer*4::               loptions
integer*4,allocatable::    options(:)

      ! weights
      if (graphtype.eq.1) then 
         ! use weights
         wgtflag = 1
      else
         ! no weights
         wgtflag = 0
      end if

      npart = nsub

      ! OPTIONS - use default options
      loptions = 5
      allocate(options(loptions))
      options = 0

      write(*,'(a,i6,a)') 'Calling METIS to divide into ',nsub,' subdomains...'
      if (nsub.eq.0) then
         write(*,'(a)') ' Illegal value of number of subdomains...'
         stop
      else if (nsub.eq.1) then
         edgecut = 0
         part = 1
      else if (nsub.gt.1 .and. nsub.le.8) then
         call METIS_PartGraphRecursive(nvertex,xadj,adjncy,vwgt,adjwgt,wgtflag,numflag,npart,options,ec,part)
      else
         call METIS_PartGraphKWay(nvertex,xadj,adjncy,vwgt,adjwgt,wgtflag,numflag,npart,options,ec,part)
      end if
      edgecut = ec

      deallocate(options)

end subroutine graph_divide

subroutine graph_pdivide_mesh(myid,nproc,comm,graphtype,neighbouring,nelem,nelem_loc,profile,&
                              inet_loc,linet_loc,nnet_loc,lnnet_loc,nsub,&
                              edgecut,part_loc,lpart_loc)
! parallel division of mesh using ParMETIS
use module_utils
implicit none
include "mpif.h"

! INPUT:
! parallel variables
integer, intent(in) :: myid, nproc, comm
! type of output graph
integer, intent(in) :: graphtype
! number of nodes to consider elements adjacent 
integer, intent(in) :: neighbouring
! number of all elements
integer, intent(in) :: nelem
! number of local elements assigned to my proc
integer, intent(in) :: nelem_loc
! info on times
logical, intent(in) :: profile 
! local INET array
integer, intent(in)   :: linet_loc
integer*4, intent(in) ::  inet_loc(linet_loc)
! local NNET array
integer, intent(in)   :: lnnet_loc
integer*4, intent(in) ::  nnet_loc(lnnet_loc)
! number of subdomains
integer, intent(in)   :: nsub
! OUTPUT:
integer,intent(out):: edgecut ! number of cut edges
integer, intent(in)    :: lpart_loc
integer*4, intent(out) ::  part_loc(lpart_loc) ! distribution

! local vars
integer :: i
integer ::             lnelempa
integer,allocatable ::  nelempa(:)

! ParMETIS vars
integer*4,parameter:: numflag = 1 ! (1 - Fortran-like arrays, 0 - C-like arrays)
integer*4::           wgtflag
integer::              lelmdist   
integer*4,allocatable:: elmdist(:)
integer ::             loptions
integer*4,allocatable::  options(:)
integer::            lwgt
real*4,allocatable :: wgt(:)
integer::              leptr   
integer*4,allocatable:: eptr(:)
integer*4:: ncon, ncommonnodes, nparts, ec
integer::            ltpwgts   
real*4,allocatable :: tpwgts(:)
integer::            lubvec   
real*4,allocatable :: ubvec(:)

! MPI vars 
integer :: ierr

      ! work according to number of subdomains
      if (nsub.eq.0) then
         ! zero subdomains is errorneous
         if (myid.eq.0) then
            write(*,'(a)') ' Illegal value of number of subdomains...'
         end if
         stop
      else if (nsub.eq.1) then
         ! one subdomain is trivial - no use of ParMETIS
         part_loc = 1
         edgecut = 0
      else
         ! use ParMETIS for more
         ! prepare arrays for ParMETIS
         ! get number of elements each processor works on
         lnelempa = nproc
         allocate(nelempa(lnelempa))
         call MPI_ALLGATHER(nelem_loc,1,MPI_INTEGER,nelempa,1,MPI_INTEGER,comm,ierr)
         ! check number of elements
         if (sum(nelempa).ne.nelem) then
            if (myid.eq.0) then
               write(*,'(a)') 'Error in number of elements for processors.'
               call flush(6)
               call error_exit
            end if
         end if

         ! ELMDIST
         lelmdist = nproc + 1
         allocate(elmdist(lelmdist))
         elmdist(1) = 1
         do i = 2,lelmdist
            elmdist(i) = elmdist(i-1) + nelempa(i-1)
         end do
         ! debug
         !write(*,*) 'elmdist',elmdist
         ! EPTR
         leptr = nelem_loc + 1
         allocate(eptr(leptr))
         eptr(1) = 1 ! In Fortran, start from 1
         do i = 2,leptr
            eptr(i) = eptr(i-1) + nnet_loc(i-1)
         end do
         ! EIND is the same as INET_LOC
         ! use weights (0 - no, 1 - yes)
         if (graphtype.eq.1) then 
            ! use weights
            wgtflag = 1
         else
            ! no weights
            wgtflag = 0
         end if
         ! WGT
         lwgt   = nelem_loc
         allocate(wgt(lwgt))
         wgt = 1.0
         ! NCON - number of constraints or weights for each vertex - influences TPWGTS
         ncon = 1
         ! NCOMMONNODES - number of nodes to call elements adjacent
         ncommonnodes = neighbouring
         ! NPARTS
         nparts = nsub
         ! TPWGTS
         ltpwgts = ncon*nparts
         allocate(tpwgts(ltpwgts))
         tpwgts = float(1)/nparts
         ! UBVEC - unbalance - refer to ParMETIS manual for proper behaviour
         lubvec = ncon
         allocate(ubvec(lubvec))
         ubvec = 1.05_4
         loptions = 3
         allocate(options(loptions))
         options(1) = 1
         if (profile) then
            options(2) = 1  ! timing info printed
         else
            options(2) = 0  ! timing info suppressed
         end if
         options(3) = 15 ! seed for random number generator

         if (myid.eq.0) then
            write(*,'(a,i6,a)') 'Calling ParMETIS to divide into ',nsub,' subdomains...'
            call flush(6)
         end if
         nparts = nsub
         call ParMETIS_V3_PartMeshKway(elmdist,eptr,inet_loc,wgt,wgtflag,numflag,ncon,ncommonnodes, nparts, tpwgts, ubvec, options,&
                                       ec, part_loc, comm)
         edgecut = ec
         if (myid.eq.0) then
            write(*,'(a)') ' ..done.'
            call flush(6)
         end if

         deallocate(options)
         deallocate(ubvec)
         deallocate(tpwgts)
         deallocate(wgt)
         deallocate(eptr)
         deallocate(elmdist)
         deallocate(nelempa)
      end if


end subroutine graph_pdivide_mesh

subroutine graph_pget_sub_neighbours(myid,nproc,comm,neighbouring,nelem,nelem_loc,nsub,nsub_loc,sub_start,&
                                     inet_loc,linet_loc,nnet_loc,lnnet_loc, iets,liets, debug, kadjsub,lkadjsub)
! parallel construction of graph of mesh and devising list of neigbouring subdomains
use module_utils
implicit none
include "mpif.h"

! INPUT:
! parallel variables
integer, intent(in) :: myid, nproc, comm
! number of nodes to consider elements adjacent 
integer, intent(in) :: neighbouring
! number of all elements
integer, intent(in) :: nelem
! number of local elements assigned to my proc
integer, intent(in) :: nelem_loc
! number of subdomains
integer, intent(in)   :: nsub
! number of subdomains on processor
integer, intent(in)   :: nsub_loc
! number of first subdomain on processor
integer, intent(in)   :: sub_start
! local INET array
integer, intent(in)   :: linet_loc
integer*4, intent(in) ::  inet_loc(linet_loc)
! local NNET array
integer, intent(in)   :: lnnet_loc
integer*4, intent(in) ::  nnet_loc(lnnet_loc)
! local IETS array
integer, intent(in)   :: liets
integer*4, intent(in) ::  iets(liets)
! debugging mode
logical,intent(in)    :: debug

! OUTPUT:
integer, intent(in)    :: lkadjsub ! nsub * nsub_loc
integer*4, intent(out) ::  kadjsub(lkadjsub) ! marked subdomains that share nodes

! local vars
integer :: i
integer ::             lnelempa
integer,allocatable ::  nelempa(:)

! ParMETIS vars
integer*4,parameter:: numflag = 1 ! (1 - Fortran-like arrays, 0 - C-like arrays)
integer::              lelmdist   
integer*4,allocatable:: elmdist(:)
integer::              leptr   
integer*4,allocatable:: eptr(:)
integer*4::     ncommonnodes
integer*4::     numdebug

! MPI vars 
integer :: ierr

      ! work according to number of subdomains
      if (nsub.eq.0) then
         ! zero subdomains is errorneous
         if (myid.eq.0) then
            write(*,'(a)') ' Illegal value of number of subdomains...'
         end if
         stop
      else if (nsub.eq.1) then
         ! one subdomain is trivial - no use of ParMETIS
         kadjsub = 0
      else
         ! use ParMETIS for more
         ! prepare arrays for ParMETIS
         ! get number of elements each processor works on
         lnelempa = nproc
         allocate(nelempa(lnelempa))
         call MPI_ALLGATHER(nelem_loc,1,MPI_INTEGER,nelempa,1,MPI_INTEGER,comm,ierr)
         ! check number of elements
         if (sum(nelempa).ne.nelem) then
            if (myid.eq.0) then
               write(*,'(a)') 'Error in number of elements for processors.'
               call flush(6)
               call error_exit
            end if
         end if

         ! ELMDIST
         lelmdist = nproc + 1
         allocate(elmdist(lelmdist))
         elmdist(1) = 1
         do i = 2,lelmdist
            elmdist(i) = elmdist(i-1) + nelempa(i-1)
         end do
         ! debug
         !write(*,*) 'elmdist',elmdist
         ! EPTR
         leptr = nelem_loc + 1
         ! debug
         !write(*,*) 'nelem_loc',nelem_loc
         !write(*,*) 'nnet_loc',nnet_loc
         allocate(eptr(leptr))
         eptr(1) = 1 ! In Fortran, start from 1
         do i = 2,leptr
            eptr(i) = eptr(i-1) + nnet_loc(i-1)
         end do
         ! EIND is the same as INET_LOC
         ! set debugging
         if (debug) then 
            numdebug = 1
         else
            numdebug = 0
         end if
         ! NCOMMONNODES - number of nodes to call elements adjacent
         ncommonnodes = neighbouring

         if (myid.eq.0) then
            write(*,'(a,i6,a)') 'Calling PGET_SUB_NEIGBOURS_C routine to find neigbours for ',nsub,' subdomains...'
            call flush(6)
         end if
         ! debug
         !write(*,*) 'elmdist',elmdist
         !write(*,*) 'eptr',eptr
         !write(*,*) 'inet_loc',inet_loc
         !write(*,*) 'numflag',numflag
         !write(*,*) 'ncommonnodes',ncommonnodes
         !call flush(6)
         call pget_sub_neighbours_c(elmdist,eptr,inet_loc,numflag,ncommonnodes,iets,liets, nsub, nsub_loc, sub_start,&
                                    kadjsub,lkadjsub, numdebug, comm)
         if (myid.eq.0) then
            write(*,'(a)') ' ..done.'
            call flush(6)
         end if

         deallocate(eptr)
         deallocate(elmdist)
         deallocate(nelempa)
      end if


end subroutine graph_pget_sub_neighbours

recursive subroutine graph_components(nvertex,xadj,lxadj,adjncy,ladjncy,components,lcomponents,ncomponents)

! find graph components in METIS-like graph

use module_utils
implicit none
integer, intent(in) :: nvertex
integer, intent(in) :: lxadj
integer, intent(in) ::  xadj(lxadj)
integer, intent(in) :: ladjncy
integer, intent(in) ::  adjncy(ladjncy)
integer, intent(in) ::    lcomponents
integer, intent(in out) :: components(lcomponents)
integer, intent(out) ::    ncomponents

! Local variable
integer :: i, icompo
integer,parameter :: sizelimit = 2000000

! initialize the components array
components = -1

! check size of graph - recursion fails for huge graphs
if (nvertex.gt.1000000) then
   write(*,*) 'GRAPH_COMPONENTS: I am recursive subroutine, use me only for graphs smaller than ',sizelimit,' vertices.'
   ncomponents = -1
   return
end if

icompo = 0
do i = 1,nvertex
   if (components(i) .le. 0) then
      icompo = icompo + 1
      call graph_components_1(i)
   end if
   ! debug
   ! print *,'components'
   ! print *,components
end do
ncomponents = icompo 
! check the components
if (any(components.eq.-1)) then
   write(*,*) 'GRAPH_COMPONENTS: Error - some vertices not associated.'
   call error_exit
end if
if (any(components.gt.ncomponents)) then
   write(*,*) 'GRAPH_COMPONENTS: Error - some component indices larger than allowed.'
   call error_exit
end if

contains

recursive subroutine graph_components_1(ivertex)

integer, intent(in) :: ivertex

!     Local variables
integer             :: ineib, indneib

! mark this node into the component
components(ivertex) = icompo

do ineib = xadj(ivertex),xadj(ivertex+1)-1
   indneib = adjncy(ineib)

   if (components(indneib).le.0) then
      call graph_components_1(indneib)
   end if
end do

end subroutine graph_components_1
end subroutine graph_components

subroutine graph_read_dimensions(idgraph,nvertex,nedge,graphtype,lxadj,ladjncy,ladjwgt)
use module_parsing
implicit none
! unit for graph file - supposed open
integer, intent(in) ::  idgraph
integer, intent(out) :: nvertex
integer, intent(out) :: nedge
integer, intent(out) :: graphtype
integer, intent(out) :: lxadj
integer, intent(out) :: ladjncy
integer, intent(out) :: ladjwgt

! read initial line in the file with graph
call rdline(idgraph)
read(line,*) nvertex, nedge, graphtype

! set dimensions for allocation
ladjncy = 2*nedge
if (graphtype.eq.1) then
   ladjwgt = ladjncy
else
   ladjwgt = 0
end if
lxadj   = nvertex + 1
end subroutine graph_read_dimensions

subroutine graph_read(idgraph,nvertex,nedge,graphtype,xadj,lxadj,adjncy,ladjncy,adjwgt,ladjwgt)
! import graph from file - without parsing (fast), no comments allowed
use module_utils
implicit none
! unit for graph file - supposed open
integer, intent(in) ::  idgraph
integer, intent(in) ::  nvertex
integer, intent(in) ::  nedge
integer, intent(in) ::  graphtype

integer, intent(in) ::  lxadj
integer, intent(out) ::  xadj(lxadj)
integer, intent(in) ::  ladjncy
integer, intent(out) ::  adjncy(ladjncy)
integer, intent(in) ::  ladjwgt
integer, intent(out) ::  adjwgt(ladjwgt)

! local variables
integer :: indadjncy, ivertex, j, nneib

indadjncy = 0
xadj(1)   = 1
do ivertex = 1,nvertex
   if      (graphtype.eq.1) then
      read(idgraph,*) nneib,(adjncy(indadjncy+j),adjwgt(indadjncy+j),j = 1,nneib)
   else if (graphtype.eq.0) then
      read(idgraph,*) nneib,(adjncy(indadjncy+j),j = 1,nneib)
   else
      write(*,*) 'GRAPH_READ: Type of graph not supported. graphtype =',graphtype
   end if
   indadjncy = indadjncy + nneib
   xadj(ivertex+1) = xadj(ivertex) + nneib
end do
! check the number of edges read correspond to the prescribed dimension
if (mod(indadjncy,2).ne.0) then
   write(*,*) 'GRAPH_READ: Error in graph - number of connections is odd but must be even.'
   call error_exit
end if
if (indadjncy/2 .ne. nedge) then
   write(*,*) 'GRAPH_READ: Error in graph - number of edges does not match.'
   call error_exit
end if
end subroutine graph_read

subroutine graph_read_parsing(idgraph,nvertex,nedge,graphtype,xadj,lxadj,adjncy,ladjncy,adjwgt,ladjwgt)
! import graph from file - with parsing (slow)
use module_parsing
use module_utils
implicit none
! unit for graph file - supposed open
integer, intent(in) ::  idgraph
integer, intent(in) ::  nvertex
integer, intent(in) ::  nedge
integer, intent(in) ::  graphtype

integer, intent(in) ::  lxadj
integer, intent(out) ::  xadj(lxadj)
integer, intent(in) ::  ladjncy
integer, intent(out) ::  adjncy(ladjncy)
integer, intent(in) ::  ladjwgt
integer, intent(out) ::  adjwgt(ladjwgt)

! local variables
integer :: indadjncy, indvertex, ivertex, nneib

indadjncy = 0
xadj(1)   = 1
do ivertex = 1,nvertex
   call rdline(idgraph)
   indvertex = 0
   call getstring
   read(string,*) nneib
   if (nneib.ne.0) then
      do
         call getstring
         indadjncy = indadjncy + 1
         indvertex = indvertex + 1
         read(string,*) adjncy(indadjncy)
         if (graphtype.eq.1) then
            call getstring
            read(string,*) adjwgt(indadjncy)
         end if
         if (kstring.eq.0) then
            exit
         end if
      end do
   end if
   ! check the line
   if (indvertex .ne. nneib) then
      write(*,*) 'GRAPH_READ_PARSING: Error in graph - number of connections mismatch.'
      write(*,*) 'GRAPH_READ_PARSING: ivertex =',ivertex,'nneib =',nneib,'indvertex =',indvertex
      call error_exit
   end if
   xadj(ivertex+1) = xadj(ivertex) + indvertex
end do
! check the number of edges read correspond to the prescribed dimension
if (mod(indadjncy,2).ne.0) then
   write(*,*) 'GRAPH_READ_PARSING: Error in graph - number of connections is odd but must be even.'
   call error_exit
end if
if (indadjncy/2 .ne. nedge) then
   write(*,*) 'GRAPH_READ_PARSING: Error in graph - number of edges does not match.'
   call error_exit
end if
end subroutine graph_read_parsing

end module module_graph
