! BDDCML - Multilevel BDDC
!
! This program is a free software.
! You can redistribute it and/or modify it under the terms of
! the GNU Lesser General Public License
! as published by the Free Software Foundation,
! either version 3 of the license,
! or (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU Lesser General Public License for more details
! <http://www.gnu.org/copyleft/lesser.html>.
!________________________________________________________________

module module_graph
! module for operations on graphs
! Jakub Sistek, Manchester, 2018
implicit none

contains

!*******************************************************************************
subroutine graph_get_dual_mesh(nelem,nnod,inet,linet,nnet,lnnet,&
                               netn,lnetn,ietn,lietn,kietn,lkietn)
!*******************************************************************************
! Get dual mesh from primal one - transpose of the sparse connectivity table.
! From Indices of Nodes on ElemenTs (INET) to Indices of ElemenTs at nodes
! (IETN).
!*******************************************************************************

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
      if (lkietn.gt.0) then
         kietn(1) = 0
      end if
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

!*******************************************************************************
subroutine graph_from_mesh(nelem,graphtype,neighbouring,inet,linet,nnet,lnnet,&
                           ietn,lietn,netn,lnetn,kietn,lkietn,&
                           nedge, xadj, adjncy, adjwgt)
!*******************************************************************************
! Construct a dual graph from the mesh and its dual mesh, corresponding to
! adjacency of elements. Corresponds to computing
! G = C C^T, where C is the connectivity matrix with rows representing elements
! and columns representing nodes.
!*******************************************************************************
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
integer, intent(out) ::  nedge
integer, allocatable, intent(out) :: xadj(:)
integer, allocatable, intent(out) :: adjncy(:)
integer, allocatable, intent(out) :: adjwgt(:)

! local variables
character(*),parameter:: routine_name = 'GRAPH_FROM_MESH'
integer,allocatable :: onerow(:), onerowweig(:)
integer :: nnetx, netnx, lonerow, lonerowweig, ie, indinet, indnode, ine
integer :: ionerow, nelmn, nne, pointietn, lorin, lorout
integer :: lxadj
integer :: ladjncy
integer :: ladjwgt

! prepare arrays for storing a row
      nnetx = maxval(nnet)
      netnx = maxval(netn)
      lonerow = nnetx * netnx
      lonerowweig = lonerow
      allocate(onerow(lonerow),onerowweig(lonerowweig))


! construct the array of pointers XADJ
      lxadj = nelem + 1
      allocate(xadj(lxadj))

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

! construct the array of edges ADJNCY and ADJWGT
      ladjncy = xadj(nelem+1)-1
      ! check the graph
      if (mod(ladjncy,2).ne.0) then
         call error(routine_name, 'Number of nodes has to be even number!')
      end if
      nedge = ladjncy / 2
      !if (graphtype.eq.1) then
      ladjwgt = ladjncy
      !else
      !   ladjwgt = 0
      !end if

      allocate(adjncy(ladjncy),adjwgt(ladjwgt))

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
            onerow(ionerow+1:ionerow + nelmn) = &
                ietn(pointietn+1:pointietn+nelmn)
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
      if (graphtype.eq.0) then
         adjwgt = 1
      end if

      deallocate(onerow,onerowweig)
end subroutine graph_from_mesh

!*******************************************************************************
subroutine graph_parse_onerow(ie,neighbouring,onerow,onerowweig,lorin,lorout)
!*******************************************************************************
! When constructing the graph, parse one row given the minimal neighbouring.
! Auxiliary routine for creating graph from mesh.
!*******************************************************************************
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
         i2 = 1
         do i = 1,lorin
            if (onerow(i) .ne. 0) then
               i2 = i
               exit
            end if
         end do
         ! shift array
         onerow = cshift(onerow,i2-1)

         ! remove multiplicities
         if (size(onerow).gt.0) then
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
         end if

         lorout = count(onerow .ne. 0)
         onerow     = pack(onerow,onerow.ne.0)
         onerowweig = pack(onerowweig,onerowweig.ne.0)
         onerow(lorout+1:)     = 0
         onerowweig(lorout+1:) = 0

end subroutine graph_parse_onerow

!*******************************************************************************
subroutine graph_check(nvertex,graphtype, xadj,lxadj, adjncy,ladjncy, &
                       adjwgt,ladjwgt)
!*******************************************************************************
! Check correctness of the graph.
!*******************************************************************************
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
integer, intent(in) ::  adjncy(ladjncy)
integer, intent(in) :: ladjwgt
integer, intent(in) ::  adjwgt(ladjwgt)

! local variables
character(*),parameter:: routine_name = 'GRAPH_CHECK'
integer :: iadjelm, indadjelm, iadjelmadj, indadjelmadj, ivertex
integer :: nadjelm, nadjelmadj
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
                     call error(routine_name, &
                                'Element multiply mentioned in the list of '// &
                                'neighbours of element', ivertex)
                  end if
                  match = .true.
                  indadjelmadj = iadjelmadj
               end if
            end do
            if (.not.match) then
               call error(routine_name, &
                          'No match! Couple mentioned but not the opposite!', &
                          ivertex)
            end if
            if (graphtype.eq.1) then
               if (adjwgt(xadj(indadjelm)-1+indadjelmadj) .ne. &
                   adjwgt(xadj(ivertex)-1+iadjelm)) then
                  call error(routine_name, &
                             'Non-matching edge weights'// &
                             'between elements ', ivertex)
                  call error_exit
               end if
            end if
         end do
      end do
end subroutine graph_check

!*******************************************************************************
subroutine graph_divide(graphtype, nvertex, &
                        xadj,lxadj, adjncy,ladjncy, &
                        vwgt,lvwgt, adjwgt,ladjwgt,&
                        nsub, contiguous_clusters, &
                        edgecut, part,lpart)
!*******************************************************************************
! Divide graph by METIS.
!*******************************************************************************

use module_utils
use iso_c_binding
implicit none
! type of output graph
integer, intent(in) :: graphtype
integer, intent(in) :: nvertex
integer, intent(in) :: lxadj
integer, intent(in) ::  xadj(lxadj)
integer, intent(in) :: ladjncy
integer, intent(in) ::  adjncy(ladjncy)
integer, intent(in) :: lvwgt
integer, intent(in) ::  vwgt(lvwgt)
integer, intent(in) :: ladjwgt
integer, intent(in) ::  adjwgt(ladjwgt)
integer, intent(in) ::   nsub
integer, intent(in) ::   contiguous_clusters
integer, intent(out) ::  edgecut
integer, intent(in) ::  lpart
integer, intent(out) ::  part(lpart)

! local vars
integer :: numflag = 1

interface
   subroutine graph_divide_c(numflag, graphtype, nvertex, &
                             xadj, lxadj, adjncy, ladjncy, &
                             vwgt, lvwgt, adjwgt, ladjwgt, &
                             nsub, contiguous_clusters, &
                             edgecut, part,lpart) &
              bind(c, name='graph_divide_c')
      use iso_c_binding
      implicit none
      integer(c_int) :: numflag
      integer(c_int) :: graphtype
      integer(c_int) :: nvertex
      integer(c_int) :: lxadj
      integer(c_int) :: xadj(lxadj)
      integer(c_int) :: ladjncy
      integer(c_int) :: adjncy(ladjncy)
      integer(c_int) :: lvwgt
      integer(c_int) :: vwgt(lvwgt)
      integer(c_int) :: ladjwgt
      integer(c_int) :: adjwgt(ladjwgt)
      integer(c_int) :: nsub
      integer(c_int) :: contiguous_clusters
      integer(c_int) :: edgecut
      integer(c_int) :: lpart
      integer(c_int) :: part(lpart)
   end subroutine graph_divide_c
end interface

      call graph_divide_c(numflag, graphtype, nvertex, &
                          xadj, lxadj, adjncy, ladjncy,&
                          vwgt, lvwgt, adjwgt, ladjwgt, &
                          nsub, contiguous_clusters, &
                          edgecut, part,lpart)

end subroutine graph_divide

!*******************************************************************************
recursive subroutine graph_components(nvertex, xadj,lxadj, adjncy,ladjncy, &
                                      components,lcomponents, ncomponents)
!*******************************************************************************
! Find graph components in a METIS-like graph.
!*******************************************************************************

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
character(*),parameter:: routine_name = 'GRAPH_COMPONENTS'
integer :: i, icompo
integer,parameter :: sizelimit = 1000000

! initialize the components array
      components = -1

! check size of graph - recursion fails for huge graphs
      if (nvertex.gt.sizelimit) then
         call warning(routine_name, &
                      'I am recursive subroutine, use me only for graphs '//&
                      'with number of vertices less than ',sizelimit)
         ncomponents = -1
         return
      end if

      icompo = 0
      do i = 1,nvertex
         if (components(i) .le. 0) then
            icompo = icompo + 1
            call graph_components_1(i)
         end if
      end do
      ncomponents = icompo

! check the components
      if (any(components.eq.-1)) then
         call error(routine_name, 'Some vertices not associated.')
      end if
      if (any(components.gt.ncomponents)) then
         call error(routine_name, ' Some component indices larger than allowed.')
      end if

contains

recursive subroutine graph_components_1(ivertex)

integer, intent(in) :: ivertex

! local vars
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

!*******************************************************************************
subroutine graph_write_to_file(idfile,nvertex,nedge,graphtype, xadj,lxadj, &
                               adjncy,ladjncy, adjwgt,ladjwgt)
!*******************************************************************************
! Write the list into file for METIS.
!*******************************************************************************

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
integer, intent(in) ::  adjncy(ladjncy)
integer, intent(in) :: ladjwgt
integer, intent(in) ::  adjwgt(ladjwgt)

! local variables
character(*),parameter:: routine_name = 'GRAPH_WRITE_TO_FILE'
integer :: ie, nadje, j

      write(idfile,'(x,i10,x,i10,x,i10)') nvertex, nedge, graphtype
      do ie = 1,nvertex
         nadje = xadj(ie+1) - xadj(ie)
         if      (graphtype.eq.0) then
            ! no weights
            write(idfile,'(600i9)') nadje, (adjncy(xadj(ie)-1+j), j = 1,nadje)
         else if (graphtype.eq.1) then
            ! weighted graph
            write(idfile,'(600i9)') nadje, &
                (adjncy(xadj(ie)-1+j), adjwgt(xadj(ie)-1+j), j = 1,nadje)
         else
            call error(routine_name, 'Graph type not supported: ', graphtype)
         end if
      end do
end subroutine graph_write_to_file

!*******************************************************************************
subroutine graph_read(idgraph,nvertex,nedge,graphtype, &
                      xadj,lxadj,adjncy,ladjncy,adjwgt,ladjwgt)
!*******************************************************************************
! Import graph from file - without parsing (fast), no comments allowed.
!*******************************************************************************
use module_parsing
use module_utils
implicit none
! unit for graph file - supposed open
integer, intent(in) ::  idgraph
integer, intent(out) ::  nvertex
integer, intent(out) ::  nedge
integer, intent(out) ::  graphtype

integer, intent(out) ::              lxadj
integer, intent(out), allocatable ::  xadj(:)
integer, intent(out) ::              ladjncy
integer, intent(out), allocatable ::  adjncy(:)
integer, intent(out) ::              ladjwgt
integer, intent(out), allocatable ::  adjwgt(:)

! local variables
character(*),parameter:: routine_name = 'GRAPH_READ'
integer :: indadjncy, ivertex, j, nneib

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

      allocate(xadj(lxadj), adjncy(ladjncy), adjwgt(ladjwgt))

      indadjncy = 0
      xadj(1)   = 1
      do ivertex = 1,nvertex
         if      (graphtype.eq.1) then
            read(idgraph,*) nneib,&
                            (adjncy(indadjncy+j),adjwgt(indadjncy+j),&
                             j = 1,nneib)
         else if (graphtype.eq.0) then
            read(idgraph,*) nneib,(adjncy(indadjncy+j),j = 1,nneib)
         else
            call error(routine_name, 'Type of graph not supported. graphtype:',&
                       graphtype)
         end if
         indadjncy = indadjncy + nneib
         xadj(ivertex+1) = xadj(ivertex) + nneib
      end do

! check the number of edges read correspond to the prescribed dimension
      if (mod(indadjncy,2).ne.0) then
         call error(routine_name, &
                    'Number of connections is odd but must be even.')
      end if
      if (indadjncy/2 .ne. nedge) then
         call error(routine_name, 'Number of edges does not match.')
      end if
end subroutine graph_read

!*******************************************************************************
subroutine graph_read_parsing(idgraph, nvertex, nedge, graphtype,&
                              xadj,lxadj, adjncy,ladjncy, adjwgt,ladjwgt)
!*******************************************************************************
! Import graph from file - with parsing. Potentially slow, but allows comments.
!*******************************************************************************
use module_parsing
use module_utils
implicit none
! unit for graph file - supposed open
integer, intent(in) ::  idgraph
integer, intent(out) ::  nvertex
integer, intent(out) ::  nedge
integer, intent(out) ::  graphtype

integer, intent(out) ::              lxadj
integer, intent(out), allocatable ::  xadj(:)
integer, intent(out) ::              ladjncy
integer, intent(out), allocatable ::  adjncy(:)
integer, intent(out) ::              ladjwgt
integer, intent(out), allocatable ::  adjwgt(:)

! local variables
character(*),parameter:: routine_name = 'GRAPH_READ_PARSING'
integer :: indadjncy, indvertex, ivertex, nneib

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

      allocate(xadj(lxadj), adjncy(ladjncy), adjwgt(ladjwgt))

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
            write(*,*) 'ivertex =',ivertex,&
                       'nneib =',nneib,&
                       'indvertex =',indvertex
            call error(routine_name, &
                       'Number of connections mismatch.')
         end if
         xadj(ivertex+1) = xadj(ivertex) + indvertex
      end do

! check the number of edges read correspond to the prescribed dimension
      if (mod(indadjncy,2).ne.0) then
         call error(routine_name, &
                    'Number of connections is odd but must be even.')
      end if
      if (indadjncy/2 .ne. nedge) then
         call error(routine_name, 'Number of edges does not match.')
      end if

end subroutine graph_read_parsing

end module module_graph
