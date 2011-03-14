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

program test_module_graph
! tester of module graph

use module_graph

implicit none

integer,parameter :: idgraph = 21

integer  :: nvertex, nedge, graphtype
integer  ::             lxadj
integer,allocatable  ::  xadj(:)
integer  ::             ladjncy
integer,allocatable  ::  adjncy(:)
integer  ::             ladjwgt
integer,allocatable  ::  adjwgt(:)
integer  ::            lcomponents
integer,allocatable  :: components(:)
integer  ::    ncomponents
integer  ::            lonerow
integer,allocatable  :: onerow(:)
integer,allocatable  :: onerowweig(:)
integer :: neighbouring, lorout

integer :: i, ie

      write(*,*) 'Test on connected graph with 1 component.'
! GRAPH - list of neighbours for elements
      open (unit=idgraph,file='test_module_graph_input_1.txt',status='old',form='formatted')
      rewind idgraph

! read initial line in the file with graph
      call graph_read_dimensions(idgraph,nvertex,nedge,graphtype,lxadj,ladjncy,ladjwgt)
      allocate(adjncy(ladjncy))
      allocate(adjwgt(ladjwgt))
      allocate(xadj(lxadj))
      call graph_read(idgraph,nvertex,nedge,graphtype,xadj,lxadj,adjncy,ladjncy,adjwgt,ladjwgt)
      close(idgraph)

      lcomponents = nvertex
      allocate(components(lcomponents))

      call graph_components(nvertex,xadj,lxadj,adjncy,ladjncy,components,lcomponents,ncomponents)
      write(*,*) 'ncomponents =', ncomponents
      write(*,*) 'vertex             component'
      do i = 1,nvertex
         write(*,*) i, components(i)
      end do
      if (all(components.eq.(/ 1, 1, 1, 1 /))) then
         write(*,*) 'Components seem O.K.'
      else   
         write(*,*) 'Components seem WRONG.'
      end if

      deallocate(components)
      deallocate(adjncy)
      deallocate(adjwgt)
      deallocate(xadj)

      write(*,*) 'Test on disconnected graph of 2 components.'
! GRAPH - list of neighbours for elements
      open (unit=idgraph,file='test_module_graph_input_2.txt',status='old',form='formatted')
      rewind idgraph

! read initial line in the file with graph
      call graph_read_dimensions(idgraph,nvertex,nedge,graphtype,lxadj,ladjncy,ladjwgt)
      allocate(adjncy(ladjncy))
      allocate(adjwgt(ladjwgt))
      allocate(xadj(lxadj))
      call graph_read(idgraph,nvertex,nedge,graphtype,xadj,lxadj,adjncy,ladjncy,adjwgt,ladjwgt)
      close(idgraph)

      lcomponents = nvertex
      allocate(components(lcomponents))

      call graph_components(nvertex,xadj,lxadj,adjncy,ladjncy,components,lcomponents,ncomponents)
      write(*,*) 'ncomponents =', ncomponents
      write(*,*) 'vertex             component'
      do i = 1,nvertex
         write(*,*) i, components(i)
      end do
      if (all(components.eq.(/ 1, 1, 2, 1 /))) then
         write(*,*) 'Components seem O.K.'
      else   
         write(*,*) 'Components seem WRONG.'
      end if

      deallocate(components)
      deallocate(adjncy)
      deallocate(adjwgt)
      deallocate(xadj)

      write(*,*) 'Test of parsing routine.'
      lonerow = 10
      allocate(onerow(lonerow),onerowweig(lonerow))
      onerow     = (/1, 4, 6, 5, 5, 3, 2, 2, 4, 5/)
      ie = 4
      write(*,*) 'Chosen element:', ie
      neighbouring = 2
      write(*,*) 'Chosen neighbouring:', neighbouring
      write(*,*) 'Array before parsing:'
      write(*,'(10i5)') onerow
      call graph_parse_onerow(ie,neighbouring,onerow,onerowweig,lonerow,lorout)
      write(*,*) 'Array after parsing:'
      write(*,'(10i5)') onerow
      write(*,*) 'Multiplicity:'
      write(*,'(10i5)') onerowweig
      write(*,*) 'Used length: ',lorout


      deallocate(onerow,onerowweig)

end program test_module_graph
