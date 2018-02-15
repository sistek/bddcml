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

module module_pp
! module for parallel preprocessing
! Jakub Sistek, Bologna, 2010
use, intrinsic :: iso_fortran_env
implicit none 

! adjustable parameters ############################
! type of real variables
      integer,parameter,private :: kr = REAL64
! debugging 
      logical,parameter,private :: debug = .false.
! profiling 
      logical,parameter,private :: profile = .false.
! adjustable parameters ############################

contains 



!**********************************************************************************
subroutine pp_pdivide_mesh(myid,nproc,comm,graphtype,neighbouring,nelem,nelem_loc,&
                           inet_loc,linet_loc,nnet_loc,lnnet_loc,nsub,&
                           edgecut,part_loc,lpart_loc)
!**********************************************************************************
! parallel division of mesh using ParMETIS
use module_utils
use iso_c_binding
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
integer*4,allocatable :: wgt(:)
integer::              leptr   
integer*4,allocatable:: eptr(:)
integer*4:: ncon, ncommonnodes, nparts, ec
integer::            ltpwgts   
real*4,allocatable :: tpwgts(:)
integer::            lubvec   
real*4,allocatable :: ubvec(:)

! MPI vars 
integer :: ierr
integer :: just_one

interface
   subroutine pdivide_mesh_c( elmdist, eptr, eind, elmwgt, wgtflag, numflag, &
                              ncon, ncommonnodes, nparts, tpwgts, ubvec, options, &
                              edgecut, part, commInt ) &
              bind(c, name='pdivide_mesh_c')
      use iso_c_binding
      implicit none
      integer(c_int) :: elmdist(*)
      integer(c_int) :: eptr(*)
      integer(c_int) :: eind(*) 
      integer(c_int) :: elmwgt(*)
      integer(c_int) :: wgtflag 
      integer(c_int) :: numflag 
      integer(c_int) :: ncon 
      integer(c_int) :: ncommonnodes 
      integer(c_int) :: nparts 
      real(c_float) :: tpwgts(ncon*nparts) 
      real(c_float) :: ubvec(ncon) 
      integer(c_int) :: options(*) 
      integer(c_int) :: edgecut 
      integer(c_int) :: part(*)
      integer(c_int) :: commInt
   end subroutine pdivide_mesh_c
end interface

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
         just_one = 1
         call MPI_ALLGATHER(nelem_loc,just_one,MPI_INTEGER,nelempa,just_one,MPI_INTEGER,comm,ierr)
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
         wgt = 1
         ! NCON - number of constraints or weights for each vertex - influences TPWGTS
         ncon = 1
         ! NCOMMONNODES - number of nodes to call elements adjacent
         ncommonnodes = neighbouring
         ! NPARTS
         nparts = nsub
         ! TPWGTS
         ltpwgts = ncon*nparts
         allocate(tpwgts(ltpwgts))
         tpwgts = 1. / float(nparts)
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

         if (myid.eq.0 .and. debug) then
            write(*,'(a,i6,a)') 'Calling ParMETIS to divide into ',nsub,' subdomains...'
            call flush(6)
         end if
         nparts = nsub
         ! portable call
         call pdivide_mesh_c(elmdist,eptr,inet_loc,wgt,wgtflag,numflag,ncon,ncommonnodes, nparts, &
                             tpwgts, ubvec, options, ec, part_loc, comm) 
         ! less portable call - works for MPI implementations using MPI_Comm int ( like mpich ), does not work with general types
         ! ( like in OpenMPI )
         !call ParMETIS_V3_PartMeshKway(elmdist,eptr,inet_loc,wgt,wgtflag,numflag,ncon,ncommonnodes, nparts, tpwgts, ubvec, options,&
         !                              ec, part_loc, comm)
         edgecut = ec
         if (myid.eq.0 .and. debug) then
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


end subroutine pp_pdivide_mesh

!**********************************************************************************
subroutine pp_pget_sub_neighbours(myid,nproc,comm,neighbouring,nelem,nelem_loc,nsub,nsub_loc,sub_start,&
                                  inet_loc,linet_loc,nnet_loc,lnnet_loc, iets,liets, debug, kadjsub,lkadjsub)
!**********************************************************************************
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

interface
   subroutine pget_sub_neighbours_c( elmdist, eptr, eind, numflag, &
                                     ncommonnodes, iets, liets, nsub, nsub_loc, sub_start, &
                                     kadjsub, lkadjsub, debug, &
                                     commInt ) &
              bind(c, name='pget_sub_neighbours_c')
      use iso_c_binding
      implicit none
      integer(c_int) :: elmdist(*)
      integer(c_int) :: eptr(*)
      integer(c_int) :: eind(*) 
      integer(c_int) :: numflag 
      integer(c_int) :: ncommonnodes 
      integer(c_int) :: liets 
      integer(c_int) :: iets(liets) 
      integer(c_int) :: nsub 
      integer(c_int) :: nsub_loc
      integer(c_int) :: sub_start
      integer(c_int) :: lkadjsub
      integer(c_int) :: kadjsub(lkadjsub)
      integer(c_int) :: debug
      integer(c_int) :: commInt
   end subroutine pget_sub_neighbours_c
end interface

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
               write(*,*) 'Error in number of elements for processors, sum(nelempa), nelem:', sum(nelempa), nelem
               write(*,*) 'nelempa:', nelempa
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
         leptr = max(leptr,2)
         ! debug
         !write(*,*) 'nelem_loc',nelem_loc
         !write(*,*) 'nnet_loc',nnet_loc
         allocate(eptr(leptr))
         eptr(1) = 1 ! In Fortran, start from 1
         if (nelem_loc .gt. 0) then
            do i = 2,leptr
               eptr(i) = eptr(i-1) + nnet_loc(i-1)
            end do
         end if
         if (nelem_loc.eq.0) then
            eptr(2) = 1
         end if

         ! EIND is the same as INET_LOC
         ! set debugging
         if (debug) then 
            numdebug = 1
         else
            numdebug = 0
         end if
         ! NCOMMONNODES - number of nodes to call elements adjacent
         ncommonnodes = neighbouring

         if (myid.eq.0 .and. debug) then
            write(*,'(a,i6,a)') 'Calling PGET_SUB_NEIGHBOURS_C routine to find neigbours for ',nsub,' subdomains...'
            call flush(6)
         end if
         ! debug
         !print *,'elmdist',elmdist
         !print *,'eptr',eptr
         !print *,'inet_loc',inet_loc
         
         !print *,'numflag',numflag
         !print *,'ncommonnodes',ncommonnodes
         !print *,'nelempa',nelempa
         !print *,'iets',iets
         !call flush(6)
         !call MPI_BARRIER(comm,ierr)

         call pget_sub_neighbours_c(elmdist,eptr,inet_loc,numflag,ncommonnodes,iets,liets, nsub, nsub_loc, sub_start,&
                                    kadjsub,lkadjsub, numdebug, comm)
         if (myid.eq.0 .and. debug) then
            write(*,'(a)') ' ..done.'
            call flush(6)
         end if

         deallocate(eptr)
         deallocate(elmdist)
         deallocate(nelempa)
      end if

end subroutine pp_pget_sub_neighbours

!************************************************************************************************
subroutine pp_get_sub_neighbours(neighbouring,nelem,nnod,nsub,inet,linet,nnet,lnnet,iets,liets,&
                                 kadjsub,lkadjsub)
!************************************************************************************************
!     Subroutine for getting the number and indices of neighbours of 
!     subdomain isub based on element graph.
!************************************************************************************************
use module_graph
use module_utils
implicit none

! INPUT:
! number of nodes to consider elements adjacent 
integer, intent(in) :: neighbouring
! number of all elements
integer, intent(in) :: nelem
! number of all nodes
integer, intent(in) :: nnod
! number of subdomains
integer, intent(in)   :: nsub

! INET array
integer, intent(in) :: linet
integer, intent(in) ::  inet(linet)
! NNET array
integer, intent(in) :: lnnet
integer, intent(in) ::  nnet(lnnet)
! IETS array
integer, intent(in) :: liets
integer, intent(in) ::  iets(liets)

! OUTPUT:
integer, intent(in)  :: lkadjsub ! nsub * nsub
integer, intent(out) ::  kadjsub(lkadjsub) ! marked subdomains that share nodes

! local vars
character(*),parameter:: routine_name = 'PP_GET_SUB_NEIGHBOURS'
! index of subdomain under 
integer :: isub

! number of elements in subdomain
integer :: nelems, nnods

! lists
integer :: lnetn, lietn, lkietn
integer,allocatable:: ietn(:), kietn(:), netn(:)

! global graph of mesh
integer :: lxadj
integer,allocatable ::  xadj(:)
integer :: ladjncy
integer,allocatable ::  adjncy(:)
integer :: ladjwgt
integer,allocatable ::  adjwgt(:)

! local subdomain mesh
integer ::            linets,   lnnets,   lisegns,   lisngns
integer,allocatable :: inets(:), nnets(:), isegns(:), isngns(:)
integer ::            ie

integer :: ie_loc, indel, indneibe, isubneib, iadje, n_graph_edge, pkadjsub

! graphtype
! 0 - graph without weights
! 1 - produce weighted graph
integer :: graphtype = 0

! Create graph of mesh
! Allocate proper sizes of two-dimensional field for list of touching elements
      lnetn  = nnod
      lietn  = linet
      lkietn = nnod
      allocate(netn(lnetn),ietn(lietn),kietn(lkietn))

! Create list of elements touching particular node
      call graph_get_dual_mesh(nelem,nnod,inet,linet,nnet,lnnet, &
                               netn,lnetn,ietn,lietn,kietn,lkietn)
! Create graph
      call graph_from_mesh(nelem,graphtype,neighbouring,inet,linet,nnet,lnnet,ietn,lietn,netn,lnetn,kietn,lkietn,&
                           n_graph_edge, xadj, adjncy, adjwgt)
      lxadj   = size(xadj)
      ladjncy = size(adjncy)
      ladjwgt = size(adjwgt)

! Check the graph
      call graph_check(nelem,graphtype, xadj,lxadj, adjncy,ladjncy, adjwgt,ladjwgt)

      deallocate(netn,ietn,kietn)

      ! loop over subdomains
      ! initialize kadjsub array
      kadjsub(:) = 0
      do isub = 1,nsub
      ! create submesh
         nelems = 0
         linets = 0
         do ie = 1,nelem
            if (iets(ie).eq.isub) then
               nelems = nelems + 1
               linets = linets + nnet(ie)
            end if
         end do
         lnnets  = nelems
         lisegns = nelems
         lisngns = linets
         allocate(inets(linets),nnets(lnnets),isegns(lisegns),isngns(lisngns))
         call pp_create_submesh(isub,nelem,inet,linet,&
                                nnet,lnnet,&
                                iets,liets,&
                                nnods,inets,linets,nnets,lnnets,&
                                isegns,lisegns,isngns,lisngns)

         ! first, use iadjs just for marking present subdomains, at the end, transform it to list
         pkadjsub = (isub-1)*nsub
         do ie_loc = 1,nelems
            indel = isegns(ie_loc)
         
            do iadje = xadj(indel),xadj(indel+1) - 1
               
               indneibe = adjncy(iadje)
               isubneib = iets(indneibe)
               if (isubneib.ne.isub) then
                  kadjsub(pkadjsub+isubneib) = 1
               end if
            end do
         end do
         deallocate(inets,nnets,isegns,isngns)
      end do

      deallocate(xadj)
      deallocate(adjncy,adjwgt)

end subroutine pp_get_sub_neighbours

!*******************************************************
subroutine pp_read_pmd_mesh_basic(idgmi,inet,linet,nnet,lnnet)
!*******************************************************
! routine that reads INET and NNET arrays of PMD mesh
use module_utils
implicit none

! INPUT:
! file with global mesh
integer, intent(in) :: idgmi
! OUTPUT:
! INET array
integer, intent(in)  :: linet
integer, intent(out) ::  inet(linet)
! NNET array
integer, intent(in)   :: lnnet
integer, intent(out) ::   nnet(lnnet)

! local vars
logical :: is_file_opened

! check that file is opened
inquire(idgmi,opened=is_file_opened)
if (.not.is_file_opened) then
   call error('PP_READ_PMD_MESH_BASIC','File idgmi not opened.')
end if

rewind idgmi
! read fields INET and NNET from GMIS file
read(idgmi,*) inet
read(idgmi,*) nnet

end subroutine pp_read_pmd_mesh_basic

!*********************************************************************
subroutine pp_get_problem_name(problemname,lproblemnamex,lproblemname)
!*********************************************************************
! routine that reads problemname from STDIN
use module_utils
implicit none

! INPUT:
! maximal langth of problem name
integer, intent(in)  :: lproblemnamex
! OUTPUT:
! name of the problem
character(lproblemnamex), intent(out) :: problemname
! actual length of the name
integer :: lproblemname

! local vars
integer :: i

      ! Name of the problem
      do i = 1,lproblemnamex
         problemname(i:i) = ' '
      end do

10    write(*,'(a)') 'Name of the problem: '
      call flush(6)
      read(*,*) problemname
      if (problemname.eq.' ') goto 10

      ! get length
      lproblemname = index(problemname,' ') - 1
      if (lproblemname.eq.-1) then
         lproblemname = lproblemnamex 
      end if
      ! pad the name with spaces
      do i = lproblemname+1,lproblemnamex
         problemname(i:i) = ' '
      end do

end subroutine pp_get_problem_name

!***************************************************************************
subroutine pp_pget_problem_name(comm,problemname,lproblemnamex,lproblemname)
!***************************************************************************
! routine that reads in parallel problemname from STDIN
use module_utils
implicit none
include "mpif.h"

! INPUT:
! communicator
integer, intent(in) :: comm
! maximal langth of problem name
integer, intent(in)  :: lproblemnamex
! OUTPUT:
! name of the problem
character(lproblemnamex), intent(out) :: problemname
! actual length of the name
integer :: lproblemname

! local vars
integer :: myid, ierr

! get my rank
call MPI_COMM_RANK(comm,myid,ierr)
! read mesh on root
if (myid.eq.0) then
   ! Name of the problem
   call pp_get_problem_name(problemname,lproblemnamex,lproblemname)
end if
! Broadcast of name of the problem      
!***************************************************************PARALLEL
call MPI_BCAST(lproblemname, 1,           MPI_INTEGER,   0, comm, ierr)
call MPI_BCAST(problemname, lproblemname, MPI_CHARACTER, 0, comm, ierr)
!***************************************************************PARALLEL

end subroutine pp_pget_problem_name

!******************************************************************************************************
subroutine pp_read_par_file(idpar, ndim, nsub, nelem, ndof, nnod, linet, tol, maxit, ndecrmax, meshdim)
!******************************************************************************************************
! routine that reads basic properties of problem from PAR file
use module_utils
implicit none

! INPUT:
! file with global boundary conditions
integer, intent(in) :: idpar
! OUTPUT:
! space dimensions
integer, intent(out) ::  ndim
! number of subdomains
integer, intent(out) ::  nsub
! number of elements
integer, intent(out) ::  nelem
! number of degrees of freedom
integer, intent(out) ::  ndof
! number of nodes
integer, intent(out) ::  nnod
! length of INET array
integer, intent(out) ::  linet
! tolerance for iterative solver
real(kr), intent(out) ::  tol
! maximal number of iterations
integer, intent(out) ::  maxit
! maximal number of iterations with nondecreasing residual
integer, intent(out) ::  ndecrmax
! dimension of mesh (NDIM is default, 2 for shells in 3D)
integer, intent(out) ::  meshdim

! local vars
logical :: is_file_opened
integer :: nnodc

! check that file is opened
   inquire(idpar,opened=is_file_opened)
   if (.not.is_file_opened) then
      call error('PP_READ_PAR_FILE','File idpar not opened.')
   end if

   rewind idpar
! read fields INET and NNET from GMIS file
   read(idpar,*) ndim, nsub, nelem, ndof, nnod, nnodc, linet, tol, maxit, ndecrmax
   read(idpar,*,err=18) meshdim
   goto 19
18 meshdim = ndim
19 continue

end subroutine pp_read_par_file

!*************************************************************************************************************
subroutine pp_pread_par_file(comm, idpar, ndim, nsub, nelem, ndof, nnod, linet, tol, maxit, ndecrmax, meshdim)
!*************************************************************************************************************
! routine that reads basic properties of problem from PAR file
use module_utils
implicit none
include "mpif.h"

! INPUT:
! communicator
integer, intent(in) :: comm
! file with global boundary conditions
integer, intent(in) :: idpar
! OUTPUT:
! space dimensions
integer, intent(out) ::  ndim
! number of subdomains
integer, intent(out) ::  nsub
! number of elements
integer, intent(out) ::  nelem
! number of degrees of freedom
integer, intent(out) ::  ndof
! number of nodes
integer, intent(out) ::  nnod
! length of INET array
integer, intent(out) ::  linet
! tolerance for iterative solver
real(kr), intent(out) ::  tol
! maximal number of iterations
integer, intent(out) ::  maxit
! maximal number of iterations with nondecreasing residual
integer, intent(out) ::  ndecrmax
! dimension of mesh (NDIM is default, 2 for shells in 3D)
integer, intent(out) ::  meshdim

! local vars
integer :: myid, ierr

! get my rank
call MPI_COMM_RANK(comm,myid,ierr)
! read mesh on root
if (myid.eq.0) then
   call pp_read_par_file(idpar, ndim, nsub, nelem, ndof, nnod, linet, tol, maxit, ndecrmax, meshdim)
end if
! Broadcast basic properties of the problem
!***************************************************************PARALLEL
call MPI_BCAST(ndim,     1,MPI_INTEGER,         0, comm, ierr)
call MPI_BCAST(nsub,     1,MPI_INTEGER,         0, comm, ierr)
call MPI_BCAST(nelem,    1,MPI_INTEGER,         0, comm, ierr)
call MPI_BCAST(ndof,     1,MPI_INTEGER,         0, comm, ierr)
call MPI_BCAST(nnod,     1,MPI_INTEGER,         0, comm, ierr)
call MPI_BCAST(linet,    1,MPI_INTEGER,         0, comm, ierr)
call MPI_BCAST(tol,      1,MPI_DOUBLE_PRECISION,0, comm, ierr)
call MPI_BCAST(maxit,    1,MPI_INTEGER,         0, comm, ierr)
call MPI_BCAST(ndecrmax, 1,MPI_INTEGER,         0, comm, ierr)
call MPI_BCAST(meshdim,  1,MPI_INTEGER,         0, comm, ierr)
!***************************************************************PARALLEL
end subroutine pp_pread_par_file

!************************************************
subroutine pp_read_pmd_bc_basic(idfvs,ifix,lifix)
!************************************************
! routine that reads IFIX array of PMD BC
use module_utils
implicit none

! INPUT:
! file with global boundary conditions
integer, intent(in) :: idfvs
! OUTPUT:
! IFIX array
integer, intent(in)  :: lifix
integer, intent(out) ::  ifix(lifix)

! local vars
logical :: is_file_opened

! check that file is opened
inquire(idfvs,opened=is_file_opened)
if (.not.is_file_opened) then
   call error('PP_READ_PMD_BC_BASIC','File idfvs not opened.')
end if

rewind idfvs
! read fields INET and NNET from GMIS file
read(idfvs,*) ifix

end subroutine pp_read_pmd_bc_basic

!*****************************************************
subroutine pp_read_pmd_bc(idfvs,ifix,lifix,fixv,lfixv)
!*****************************************************
! routine that reads IFIX and FIXV arrays of PMD BC
use module_utils
implicit none

! INPUT:
! file with global boundary conditions
integer, intent(in) :: idfvs
! OUTPUT:
! IFIX array
integer, intent(in)  :: lifix
integer, intent(out) ::  ifix(lifix)
! FIXV array
integer,  intent(in)  :: lfixv
real(kr), intent(out) ::  fixv(lfixv)

! read IFIX
call pp_read_pmd_bc_basic(idfvs,ifix,lifix)
! read FIXV
read(idfvs,*) fixv

end subroutine pp_read_pmd_bc

!******************************************************
subroutine pp_pread_pmd_bc_basic(comm,idfvs,ifix,lifix)
!******************************************************
! routine that reads in parallel IFIX array of PMD BC
use module_utils
implicit none
include "mpif.h"

! INPUT:
! communicator
integer, intent(in) :: comm
! file with global BC
integer, intent(in) :: idfvs
! OUTPUT:
! IFIX array
integer, intent(in)  :: lifix
integer, intent(out) ::  ifix(lifix)

! local vars
integer :: myid, ierr

! get my rank
call MPI_COMM_RANK(comm,myid,ierr)
! read mesh on root
if (myid.eq.0) then
   call pp_read_pmd_bc_basic(idfvs,ifix,lifix)
end if
!***************************************************************PARALLEL
call MPI_BCAST(ifix, lifix, MPI_INTEGER, 0, comm, ierr)
!***************************************************************PARALLEL
end subroutine pp_pread_pmd_bc_basic

!***********************************************************
subroutine pp_pread_pmd_bc(comm,idfvs,ifix,lifix,fixv,lfixv)
!***********************************************************
! routine that reads in parallel IFIX array of PMD BC
use module_utils
implicit none
include "mpif.h"

! INPUT:
! communicator
integer, intent(in) :: comm
! file with global BC
integer, intent(in) :: idfvs
! OUTPUT:
! IFIX array
integer, intent(in)  :: lifix
integer, intent(out) ::  ifix(lifix)
! FIXV array
integer, intent(in)   :: lfixv
real(kr), intent(out) ::  fixv(lfixv)

! local vars
integer :: myid, ierr

! get my rank
call MPI_COMM_RANK(comm,myid,ierr)
! read mesh on root
if (myid.eq.0) then
   ! read IFIX
   call pp_read_pmd_bc_basic(idfvs,ifix,lifix)
   ! read FIXV
   read(idfvs,*) fixv
end if
!***************************************************************PARALLEL
call MPI_BCAST(ifix, lifix, MPI_INTEGER, 0, comm, ierr)
call MPI_BCAST(fixv, lfixv, MPI_DOUBLE_PRECISION, 0, comm, ierr)
!***************************************************************PARALLEL
end subroutine pp_pread_pmd_bc

!*****************************************
subroutine pp_read_pmd_rhs(idrhs,rhs,lrhs)
!*****************************************
! routine that reads RHS array of PMD right-hand side
use module_utils
implicit none

! INPUT:
! file with global boundary conditions
integer, intent(in) :: idrhs
! OUTPUT:
! RHS array
integer,  intent(in)  :: lrhs
real(kr), intent(out) ::  rhs(lrhs)

! local vars
logical :: is_file_opened

! check that file is opened
inquire(idrhs,opened=is_file_opened)
if (.not.is_file_opened) then
   call error('PP_READ_PMD_RHS','File idrhs not opened.')
end if
! read RHS
read(idrhs) rhs

end subroutine pp_read_pmd_rhs

!***********************************************
subroutine pp_pread_pmd_rhs(comm,idrhs,rhs,lrhs)
!***********************************************
! routine that reads in parallel RHS array of PMD right-hand side
use module_utils
implicit none
include "mpif.h"

! INPUT:
! communicator
integer, intent(in) :: comm
! file with global BC
integer, intent(in) :: idrhs
! OUTPUT:
! RHS array
integer, intent(in)   :: lrhs
real(kr), intent(out) ::  rhs(lrhs)

! local vars
integer :: myid, ierr

! get my rank
call MPI_COMM_RANK(comm,myid,ierr)
! read mesh on root
if (myid.eq.0) then
   ! read IFIX
   call pp_read_pmd_rhs(idrhs,rhs,lrhs)
end if
!***************************************************************PARALLEL
call MPI_BCAST(rhs, lrhs, MPI_DOUBLE_PRECISION, 0, comm, ierr)
!***************************************************************PARALLEL
end subroutine pp_pread_pmd_rhs

!*******************************************************
subroutine pp_pread_pmd_mesh_basic(comm,idgmi,inet,linet,nnet,lnnet)
!*******************************************************
! routine that reads in parallel INET and NNET arrays of PMD mesh
use module_utils
implicit none
include "mpif.h"

! INPUT:
! communicator
integer, intent(in) :: comm
! file with global mesh
integer, intent(in) :: idgmi
! OUTPUT:
! INET array
integer, intent(in)  :: linet
integer, intent(out) ::  inet(linet)
! NNET array
integer, intent(in)   :: lnnet
integer, intent(out) ::   nnet(lnnet)

! local vars
integer :: myid, ierr

! get my rank
call MPI_COMM_RANK(comm,myid,ierr)
! read mesh on root
if (myid.eq.0) then
   call pp_read_pmd_mesh_basic(idgmi,inet,linet,nnet,lnnet)
end if
!***************************************************************PARALLEL
      call MPI_BCAST(inet, linet, MPI_INTEGER, 0, comm, ierr)
      call MPI_BCAST(nnet, lnnet, MPI_INTEGER, 0, comm, ierr)
!***************************************************************PARALLEL
end subroutine pp_pread_pmd_mesh_basic

!****************************************************************************
subroutine pp_read_pmd_mesh(idgmi,inet,linet,nnet,lnnet,nndf,lnndf,xyz,lxyz1,lxyz2)
!****************************************************************************
! routine to read complete PMD mesh data
use module_utils
implicit none
include "mpif.h"

! INPUT:
! file with global mesh
integer, intent(in) :: idgmi
! OUTPUT:
! INET array
integer, intent(in)  :: linet
integer, intent(out) ::  inet(linet)
! NNET array
integer, intent(in)   :: lnnet
integer, intent(out) ::   nnet(lnnet)
! NNDF array
integer, intent(in)   :: lnndf
integer, intent(out) ::   nndf(lnndf)
! XYZ array
integer, intent(in)   :: lxyz1, lxyz2
real(kr), intent(out) ::  xyz(lxyz1,lxyz2)

rewind idgmi
! read fields INET and NNET from GMIS file
call pp_read_pmd_mesh_basic(idgmi,inet,linet,nnet,lnnet)
read(idgmi,*) nndf
read(idgmi,*) xyz

end subroutine pp_read_pmd_mesh

!*******************************************************
subroutine pp_pread_pmd_mesh(comm,idgmi,inet,linet,nnet,lnnet,nndf,lnndf,xyz,lxyz1,lxyz2)
!*******************************************************
! routine that reads in parallel all PMD mesh
use module_utils
implicit none
include "mpif.h"

! INPUT:
! communicator
integer, intent(in) :: comm
! file with global mesh
integer, intent(in) :: idgmi
! OUTPUT:
! INET array
integer, intent(in)  :: linet
integer, intent(out) ::  inet(linet)
! NNET array
integer, intent(in)   :: lnnet
integer, intent(out) ::   nnet(lnnet)
! NNDF array
integer, intent(in)   :: lnndf
integer, intent(out) ::   nndf(lnndf)
! XYZ array
integer, intent(in)   :: lxyz1, lxyz2
real(kr), intent(out) ::  xyz(lxyz1,lxyz2)

! local vars
integer :: myid, ierr

! get my rank
call MPI_COMM_RANK(comm,myid,ierr)
! read mesh on root
if (myid.eq.0) then
   call pp_read_pmd_mesh(idgmi,inet,linet,nnet,lnnet,nndf,lnndf,xyz,lxyz1,lxyz2)
end if
!***************************************************************PARALLEL
      call MPI_BCAST(inet, linet, MPI_INTEGER, 0, comm, ierr)
      call MPI_BCAST(nnet, lnnet, MPI_INTEGER, 0, comm, ierr)
      call MPI_BCAST(nndf, lnndf, MPI_INTEGER, 0, comm, ierr)
      call MPI_BCAST(xyz, lxyz1*lxyz2, MPI_DOUBLE_PRECISION, 0, comm, ierr)
!***************************************************************PARALLEL
end subroutine pp_pread_pmd_mesh

!**********************************************************************************
subroutine pp_divide_mesh(graphtype,correct_division,neighbouring,nelem,nnod,&
                          inet,linet,nnet,lnnet,nsub,contiguous_subdomains,&
                          edgecut,part,lpart)
!**********************************************************************************
! serial division of mesh using METIS
use module_graph
use module_utils
implicit none

! INPUT:
! type of output graph
integer, intent(in) :: graphtype
! correct division to make subdomain continuous?
logical, intent(in) :: correct_division
! number of nodes to consider elements adjacent 
integer, intent(in) :: neighbouring
! number of elements
integer, intent(in) :: nelem
! number of nodes
integer, intent(in) :: nnod
! INET array
integer, intent(in)   :: linet
integer*4, intent(in) ::  inet(linet)
! NNET array
integer, intent(in)   :: lnnet
integer*4, intent(in) ::  nnet(lnnet)
! number of subdomains
integer, intent(in)   :: nsub
! should clusters be contiguous?
integer, intent(in)   :: contiguous_subdomains
! OUTPUT:
integer,intent(out):: edgecut ! number of cut edges
integer, intent(in)    :: lpart
integer*4, intent(out) ::  part(lpart) ! distribution

! local vars
character(*),parameter:: routine_name = 'PP_DIVIDE_MESH'

! type of algorithm for the division
! 0 - graph-based
! 1 - index-based
integer :: division_type = 0

select case (division_type)
   case (0)
      call pp_divide_mesh_graph(graphtype,correct_division,neighbouring,nelem,nnod,&
                                inet,linet,nnet,lnnet,nsub,contiguous_subdomains,&
                                edgecut,part,lpart)
   case (1)
      call pp_divide_mesh_chunks(nelem, nsub, part,lpart)
   case default
      call error(routine_name, 'Unknown division type.')
end select



end subroutine pp_divide_mesh

!**********************************************************************************
subroutine pp_divide_mesh_chunks(nelem, nsub, part,lpart)
!**********************************************************************************
! serial division of mesh using METIS
use module_graph
use module_utils
implicit none

! INPUT:
! number of elements
integer, intent(in) :: nelem
! number of subdomains
integer, intent(in)   :: nsub
! OUTPUT:
integer, intent(in)    :: lpart
integer, intent(out) ::  part(lpart) ! distribution

integer :: ie, isub
integer :: lel2sub
integer,allocatable :: el2sub(:)

      lel2sub = nsub + 1
      allocate(el2sub(lel2sub))

      call pp_distribute_linearly(nelem, nsub, el2sub,lel2sub)

      do isub = 1,nsub
         do ie = el2sub(isub), el2sub(isub+1)-1
            part(ie) = isub
         end do
      end do

      deallocate(el2sub)
      
end subroutine pp_divide_mesh_chunks

!**********************************************************************************
subroutine pp_divide_mesh_graph(graphtype,correct_division,neighbouring,nelem,nnod,&
                                inet,linet,nnet,lnnet,nsub,contiguous_subdomains,&
                                edgecut,part,lpart)
!**********************************************************************************
! serial division of mesh using METIS
use module_graph
use module_utils
implicit none

! INPUT:
! type of output graph
integer, intent(in) :: graphtype
! correct division to make subdomain continuous?
logical, intent(in) :: correct_division
! number of nodes to consider elements adjacent 
integer, intent(in) :: neighbouring
! number of elements
integer, intent(in) :: nelem
! number of nodes
integer, intent(in) :: nnod
! INET array
integer, intent(in)   :: linet
integer*4, intent(in) ::  inet(linet)
! NNET array
integer, intent(in)   :: lnnet
integer*4, intent(in) ::  nnet(lnnet)
! number of subdomains
integer, intent(in)   :: nsub
! should clusters be contiguous?
integer, intent(in)   :: contiguous_subdomains
! OUTPUT:
integer,intent(out):: edgecut ! number of cut edges
integer, intent(in)    :: lpart
integer*4, intent(out) ::  part(lpart) ! distribution

! local vars
character(*),parameter:: routine_name = 'PP_DIVIDE_MESH'
integer:: nedge, ncomponents
integer :: lnetn, lietn, lkietn
integer,allocatable:: ietn(:), kietn(:), netn(:)
integer*4::            ladjncy,   ladjwgt,   lxadj,   lvwgt
integer*4,allocatable:: adjncy(:), adjwgt(:), xadj(:), vwgt(:)
integer*4::            ladjncys,   ladjwgts,   lxadjs
integer*4,allocatable:: adjncys(:), adjwgts(:), xadjs(:)
integer ::           lcomponents
integer,allocatable:: components(:)
integer*4::             nelems, nnods
integer*4::             linets,   lnnets,   lietns,   lnetns,   lkietns,   lisegns,   lisngns
integer*4,allocatable::  inets(:), nnets(:), ietns(:), netns(:), kietns(:), isegns(:), isngns(:)
integer ::             lsubcomponents
integer,allocatable::   subcomponents(:)
integer ::             lsubpool
integer,allocatable::   subpool(:)
integer ::             lnelemsa
integer,allocatable::   nelemsa(:)
integer :: ie, isub, nedges, nsubcomponents, iadje, icomp, ies,&
           indcompx, indel, indneibe, indsubmin, isubneib, jsub,&
           necomp, necompx, nelemsadj, nesubmin, nneib
integer :: nelemsmax, nelemsmin
logical :: one_more_check_needed = .false.
real(kr) :: imbalance

      call info( routine_name, 'Building a graph ... ' )
! Allocate proper sizes of two-dimensional field for list of touching elements
      lnetn  = nnod
      lietn  = linet
      lkietn = nnod
      allocate(netn(lnetn),ietn(lietn),kietn(lkietn))

! Create list of elements touching particular node
      call graph_get_dual_mesh(nelem,nnod,inet,linet,nnet,lnnet, &
                               netn,lnetn,ietn,lietn,kietn,lkietn)
! Create graph
      call graph_from_mesh(nelem,graphtype,neighbouring,inet,linet,nnet,lnnet,ietn,lietn,netn,lnetn,kietn,lkietn,&
                           nedge, xadj, adjncy, adjwgt)
      lxadj   = size(xadj)
      ladjncy = size(adjncy)
      ladjwgt = size(adjwgt)
      call info( routine_name, 'done.' )
      ! free memory
      deallocate(netn,ietn,kietn)

      call info( routine_name, 'Check graph ... ' )
      !call info( routine_name, 'Neighbouring is: ', neighbouring )
! Check the graph
      call graph_check(nelem,graphtype, xadj,lxadj, adjncy,ladjncy, adjwgt,ladjwgt)
      call info( routine_name, 'done.' )
! Check components of the graph
      call info( routine_name, 'Check mesh components ... ' )
      lcomponents = nelem
      !print *, 'xadj', xadj
      !print *, 'adjncy', adjncy
      allocate(components(lcomponents))
      call graph_components(nelem, xadj,lxadj, adjncy,ladjncy, components,lcomponents, ncomponents)
      if (ncomponents.eq.1) then
         call info( routine_name, 'mesh seems continuous.' )
      else if (ncomponents.eq.-1) then
         call warning( routine_name, 'failed for size of mesh.' )
      else
         call warning( routine_name, 'mesh discontinuous - number of components: ',ncomponents )
      end if
      deallocate(components)


! Divide graph
      lvwgt = nelem
      allocate(vwgt(lvwgt))
      vwgt = 1
      call graph_divide(graphtype,nelem,xadj,lxadj,adjncy,ladjncy,vwgt,lvwgt,adjwgt,ladjwgt,&
                        nsub,contiguous_subdomains,edgecut,part,lpart)
      deallocate(vwgt)
      call info( routine_name, 'resulting number of cut edges:', edgecut )

      ! count numbers of elements in subdomains
      lnelemsa = nsub
      allocate(nelemsa(lnelemsa))
      do isub = 1,nsub
         nelemsa(isub) = count(part.eq.isub)
      end do

      call info( routine_name, 'CHECKING OF DIVISION: ' )

  123 continue

      nelemsmax = maxval(nelemsa)
      nelemsmin = minval(nelemsa)
      if (nelemsmin > 0) then
         imbalance = dble(nelemsmax-nelemsmin)/dble(nelemsmin)
      else
         imbalance = -1.
      endif
      call info( routine_name, 'final quality of division:' )
      call info( routine_name, 'largest subdomain contains elements :',nelemsmax )
      call info( routine_name, 'smallest subdomain contains elements :',nelemsmin )
      call info( routine_name, 'imbalance of division [%]:',imbalance*100 )

! check division into subdomains
      one_more_check_needed = .false.
      do isub = 1,nsub
         nelems = 0
         linets = 0
         do ie = 1,nelem
            if (part(ie).eq.isub) then
               nelems = nelems + 1
               linets = linets + nnet(ie)
            end if
         end do

         ! continue only if number of elements in the subdomain is larger than 0
         if (nelems == 0) then
            cycle
         end if

         lnnets  = nelems
         lisegns = nelems
         lisngns = linets
         allocate(inets(linets),nnets(lnnets),isegns(lisegns),isngns(lisngns))

         call pp_create_submesh(isub,nelem,inet,linet,nnet,lnnet,part,lpart,&
                                nnods,inets,linets,nnets,lnnets,&
                                isegns,lisegns,isngns,lisngns)
         lnetns  = nnods
         lietns  = linets
         lkietns = nnods
         allocate(netns(lnetns),ietns(lietns),kietns(lkietns))

! Create list of elements touching particular node
         call graph_get_dual_mesh(nelems,nnods,inets,linets,nnets,lnnets, &
                                  netns,lnetns,ietns,lietns,kietns,lkietns)
! Create graph
         call graph_from_mesh(nelems,graphtype,neighbouring,inets,linets,nnets,lnnets,&
                              ietns,lietns,netns,lnetns,kietns,lkietns,&
                              nedges, xadjs, adjncys, adjwgts)
         lxadjs   = size(xadjs)
         ladjncys = size(adjncys)
         ladjwgts = size(adjwgts)

! Check the graph
         call graph_check(nelems,graphtype, xadjs,lxadjs, adjncys,ladjncys, adjwgts,ladjwgts)
! Check components of the graph
         lsubcomponents = nelems
         allocate(subcomponents(lsubcomponents))
         call graph_components(nelems, xadjs,lxadjs, adjncys,ladjncys, subcomponents,lsubcomponents, nsubcomponents)
         if (nsubcomponents.gt.1) then
            call warning(routine_name, 'Subdomain discontinuous: ', isub )
            call warning(routine_name, '   - number of components: ', nsubcomponents )

            if (correct_division) then
               one_more_check_needed = .true.
               lsubpool = nsub
               allocate(subpool(lsubpool))
               ! glue discontinuous subdomains to neigbouring parts
               ! find largest subdomain - only this one will be assigned this number
               necompx  = 0
               indcompx = 0
               do icomp = 1,nsubcomponents
                  necomp = count(subcomponents .eq. icomp)
                  call info(routine_name, 'number of elements in subcomponent:',necomp )
                  if (necomp.gt.necompx) then
                     necompx  = necomp
                     indcompx = icomp
                  end if
               end do
               ! glue nondominant components to neighbouring subdomains
               do icomp = 1,nsubcomponents
                  if (icomp.ne.indcompx) then

                     ! find candidate subdomain
                     subpool = 0
                     do ies = 1,nelems
                        if (subcomponents(ies).eq.icomp) then
                           ! look into the graph to find neigbouring subdomains
                           indel = isegns(ies)
                           nneib = xadj(indel+1) - xadj(indel) 
                           ! mark subdomains
                           do iadje = xadj(indel),xadj(indel+1) - 1
                              indneibe = adjncy(iadje)
                              isubneib = part(indneibe)
                              subpool(isubneib) =  1
                           end do
                        end if
                     end do

                     nesubmin  = maxval(nelemsa)
                     indsubmin = 0
                     do jsub = 1,nsub
                        if (subpool(jsub).gt.0 .and. jsub.ne.isub) then
                           nelemsadj = nelemsa(jsub)
                           if (nelemsadj.le.nesubmin) then
                              indsubmin = jsub
                              nesubmin =  nelemsadj
                           end if
                        end if
                     end do
                     if (indsubmin.eq.0) then
                        call error( routine_name, ' Error in finding candidate subdomain for merging.' )
                     end if

                     ! join the component to the candidate subdomain
                     do ies = 1,nelems
                        if (subcomponents(ies).eq.icomp) then
                           indel = isegns(ies)

                           part(indel) = indsubmin
                           nelemsa(indsubmin) = nelemsa(indsubmin) + 1
                           nelemsa(isub)      = nelemsa(isub) - 1
                        end if
                     end do

                  end if
               end do
               deallocate(subpool)
            end if ! correct division

         end if

         deallocate(subcomponents)
         deallocate(adjncys,adjwgts)
         deallocate(xadjs)
         deallocate(netns,ietns,kietns)
         deallocate(inets,nnets,isegns,isngns)
      end do

      ! if some part of subdomains were glued to other subdomains, make one more check
      if (one_more_check_needed) then
         call info ( routine_name, 'Some parts were redistributed, perform one more check.' )
         goto 123
      end if

      deallocate(nelemsa)

      deallocate(adjncy,adjwgt)
      deallocate(xadj)

end subroutine pp_divide_mesh_graph

!************************************************************************
subroutine pp_create_submesh(isub,nelem,inet,linet,nnet,lnnet,iets,liets,&
                             nnods,inets,linets,nnets,lnnets,&
                             isegns,lisegns,isngns,lisngns)
!************************************************************************
!     Subroutine for creation of subdomain mesh description
!************************************************************************
use module_utils
implicit none

! number of elements and nodes
integer,intent(in) :: nelem
! subdomain number
integer,intent(in) :: isub
! global mesh description ala PMD
integer,intent(in) :: linet,       lnnet
integer,intent(in) ::  inet(linet), nnet(lnnet)
! division into subdomains
integer,intent(in) :: liets
integer,intent(in) ::  iets(liets)
! number of elements and nodes
integer,intent(out) :: nnods
! subdomain mesh description
integer,intent(in)  :: linets,        lnnets,        lisegns,         lisngns
integer,intent(out) ::  inets(linets), nnets(lnnets), isegns(lisegns), isngns(lisngns)

! local variables
integer :: ie, inod, inods, nne, pointinet, pointinets, indlocel

! Creation of field inets(linets) - still with global pointers
      pointinet  = 0
      pointinets = 0
      indlocel   = 0
      do ie = 1,nelem
         nne = nnet(ie)
         if (iets(ie).eq.isub) then
            inets(pointinets+1:pointinets+nne) = -inet(pointinet+1:pointinet+nne)

            indlocel = indlocel + 1

            nnets(indlocel)  = nne
            isegns(indlocel) = ie

            pointinets = pointinets + nne
         end if
         pointinet = pointinet + nne
      end do

! get imbedding of subdomain nodes into global nodes
      ! copy array
      isngns = -inets
      call iquick_sort(isngns,lisngns)
      call get_array_norepeat(isngns,lisngns,nnods)

      do inods = 1,nnods
         inod = isngns(inods)
         isngns(inods) = inod
         ! fix inet
         where (inets.eq.-inod) inets = inods
      end do

end subroutine pp_create_submesh

!**************************************************************
subroutine pp_get_globs(ndim,meshdim,nelem,nnod,nsub,&
                        inet,linet,nnet,lnnet,nndf,lnndf,xyz,lxyz1,lxyz2,&
                        remove_bc_nodes,ifix,lifix,iets,liets,ncornersmin,&
                        nnodi, ncorners,nedges,nfaces,&
                        kglobs,lkglobs, typeglobs,ltypeglobs)
!**************************************************************
!     Subroutine for selection of corners, edges and faces
!**************************************************************
use module_utils
implicit none

! INPUT:
! number of dimensions
integer,intent(in) :: ndim
! dimension of mesh ( 2 for shells, 3 for hexas)
integer,intent(in) :: meshdim
! number of elements
integer,intent(in) :: nelem
! number of nodes
integer,intent(in) :: nnod
! number of subdomains
integer,intent(in) :: nsub
! global mesh description ala PMD
integer,intent(in) :: linet,       lnnet,       lnndf
integer,intent(in) ::  inet(linet), nnet(lnnet), nndf(lnndf)
integer,intent(in) :: lxyz1, lxyz2
real(kr),intent(in) :: xyz(lxyz1,lxyz2)
! should nodes on Dirichlet BC be removed?
logical,intent(in) :: remove_bc_nodes
! indices of Dirichlet BC 
integer,intent(in) :: lifix
integer,intent(in) ::  ifix(lifix)
! division into subdomains
integer,intent(in) :: liets
integer,intent(in) ::  iets(liets)
! minimal number of selected corners - if resulting number is below this value, additional corners are generated randomly
integer,intent(in) :: ncornersmin
! OUTPUT:
! number of nodes at interface
integer,intent(out) :: nnodi
! number of corners - if value on input, this value is generated randomly and added to corners by CSA (corner selecting algorithm),
!                     otherwise, it is given a value on output
integer,intent(inout) :: ncorners
! number of edges
integer,intent(out) :: nedges
! number of faces
integer,intent(out) :: nfaces
! array with numbers of globs for each node
integer,intent(in) :: lkglobs
integer,intent(out) :: kglobs(lkglobs)
! array with type of globs for each node (3 - corner, 2 - edge, 1 - face)
integer,intent(in) :: ltypeglobs
integer,intent(out) :: typeglobs(ltypeglobs)

! local vars
integer:: inod, jnodi, jnod, isub, jsub, inodi, iglob, ing, iglobnode,&
          i, nmax, nng, nplaces, ie,&
          indsub, nnewvertex, ncornerss, nnodis, indisn, inodis, inodcs, inodcs1, &
          inodcs2, inodcs3, indi, inewnodes, indinode, ncorrections, &
          nglobn, ncornerspair, nshared, inodsh, ish, nnsub, iisub, &
          indjsub, inddof, ndofn, point, indnod, nremoved_corners_bc,&
          indnodeg, inodc, inc, nvertex
real(kr):: x1, y1, z1, x2, y2, z2, rnd, xish, yish, zish

integer:: indaux(1)

integer ::            lkinet
integer,allocatable :: kinet(:)
integer ::            lkmynodes,   lkinodes,   lkneibnodes
integer,allocatable :: kmynodes(:), kinodes(:), kneibnodes(:)
integer ::            ligingn
integer,allocatable :: igingn(:)
integer ::            lkdof
integer,allocatable :: kdof(:)
integer ::       lsublist1, lsublist2
integer,allocatable :: sublist(:,:)
integer ::            lnodeglob
integer,allocatable :: nodeglob(:)
integer ::            lnewvertex
integer,allocatable :: newvertex(:)
integer::             lsubneib
integer,allocatable::  subneib(:)
integer ::            licheck
integer,allocatable :: icheck(:)
integer::            lisingin
integer,allocatable:: isingin(:)
integer::            lplayground
integer,allocatable:: playground(:)
integer::             ldist,   larea
real(kr),allocatable:: dist(:), area(:)
integer::             lxyzsh1, lxyzsh2 
real(kr),allocatable:: xyzsh(:,:)
integer::             lxyzbase
real(kr),allocatable:: xyzbase(:)

type glob_type
   integer ::             itype
   logical ::             selected
   integer ::             nnod
   integer,allocatable :: nodes(:)
   integer ::             nsub
   integer,allocatable :: subdomains(:)
end type glob_type
integer  :: nglobs
integer::                      lglobs
type(glob_type), allocatable :: globs(:)

type neighbouring_type
   integer ::             nnsub
   integer,allocatable :: list(:)
end type neighbouring_type
integer::                      lneighbourings
type(neighbouring_type), allocatable :: neighbourings(:)

! Creation of field KINET(NELEM) with addresses before first global node of element IE in field inet
      lkinet = nelem
      allocate(kinet(lkinet))
      kinet(1) = 0
      do ie = 2,nelem
         kinet(ie) = kinet(ie-1) + nnet(ie-1)
      end do

! Recognize interface
      lkinodes  = nnod
      lkmynodes = nnod
      lkneibnodes = nnod
      allocate(kinodes(lkinodes),kmynodes(lkmynodes),kneibnodes(lkneibnodes))
! Begin loop over subdomains
      kinodes = 0
      do isub = 1,nsub
! Creation of field kmynodes(nnod) - local usage of global nodes
! if global node is in my subdomain, assigned 1, else remains 0
         call pp_mark_sub_nodes(isub,nelem,iets,liets,inet,linet,nnet,lnnet,&
                                kinet,lkinet,kmynodes,lkmynodes)
! Mark interface nodes
         where(kmynodes.eq.1) kinodes = kinodes + 1
      end do
! Count interface nodes
      nnodi = count(kinodes.gt.1)
      if (debug) then
         write(*,'(a,i8)') ' number of nodes on interface nnodi =',nnodi
      end if
      if (ncornersmin .gt. nnodi) then
            call error('PP_GET_GLOBS','demanded number of corners larger than number of interface nodes.')
      end if

! Create field of indices of global interface nodes in global numbering
      ligingn = nnodi
      allocate(igingn(ligingn))
      inodi = 0
      do inod = 1,nnod
         if (kinodes(inod).gt.1) then
            inodi = inodi + 1
            igingn(inodi) = inod
         end if
      end do

! What is the maximal number of subdomains sharing one node?
      nmax = maxval(kinodes)
! How many nodes are shared by the maximal number of subdomains?
      nplaces = count(kinodes.eq.nmax)
      if (debug) then
         write(*,'(a,i3,a,i7,a)') 'Maximum number of touching subdomains nmax =',nmax,' on ',nplaces,' places.'
      end if

! Recognize faces and edges
! Create list of subdomains that a node is included in
      lsublist1 = nnodi
      lsublist2 = nmax
      allocate(sublist(lsublist1,lsublist2))
      sublist = 0
      kinodes = 0
! Begin loop over subdomains
      do isub = 1,nsub
! Creation of field kmynodes(nnod) - local usage of global nodes
! if global node is in my subdomain, assigned 1, else remains 0
         call pp_mark_sub_nodes(isub,nelem,iets,liets,inet,linet,nnet,lnnet,&
                                kinet,lkinet,kmynodes,lkmynodes)
         do inodi = 1,nnodi
            inod = igingn(inodi)
            if (kmynodes(inod).eq.1) then
               kinodes(inod) = kinodes(inod) + 1
               indsub = kinodes(inod)
               sublist(inodi,indsub) = isub
            end if
         end do
      end do

!      do i = 1,nnod
!         write(*,*) (sublist(i,j), j = 1,nmax)
!      end do

! Identify topology of globs
      if (debug) then
         write(*,'(a)') '   identify globs by equivalence classes...'
      end if
      lnodeglob = nnodi
      allocate(nodeglob(lnodeglob))
      nodeglob = -1
      iglob = 0
      do inodi = 1,nnodi
         inod = igingn(inodi)
         if (count(sublist(inodi,:).ne.0).le.1) then
            call error('GETGLOBS','subdomain list - only one subdomain recognized at interface node')
         end if
         if(nodeglob(inodi).eq.-1) then
! Node was not yet assigned into a glob - so there is a new glob to define
            iglob = iglob + 1
            do jnodi = inodi,nnodi
               jnod = igingn(jnodi)
! Go through the remaining nodes
               if (all(sublist(jnodi,:).eq.sublist(inodi,:))) then
                  nodeglob(jnodi) = iglob
               end if
            end do
         end if
      end do
      nglobs = iglob
      if (debug) then
         write(*,'(a,i8)') 'Total number of globs is',nglobs
      end if
      if (count(nodeglob.eq.-1).gt.0) then
         call error('GETGLOBS','in nodeglob - unassociated node of interface with glob')
      end if
      
! glob type - 1 face, 2 edge, 3 vertex 
      lglobs = nglobs
      allocate(globs(lglobs))
      globs(:)%itype = 0
      globs(:)%selected = .true.
      do iglob = 1,nglobs
!        count nodes in each glob and allocate memory for nodes
         globs(iglob)%nnod = count(nodeglob.eq.iglob)
         allocate(globs(iglob)%nodes(globs(iglob)%nnod))
!        find nodes in the glob         
         iglobnode = 0
         do inodi = 1,nnodi
            if (nodeglob(inodi).eq.iglob) then
               iglobnode = iglobnode + 1
               globs(iglob)%nodes(iglobnode) = inodi
            end if
         end do
!        count subdomains in each glob and allocate memory for their list
         globs(iglob)%nsub = count(sublist(globs(iglob)%nodes(1),:).gt.0)
         allocate(globs(iglob)%subdomains(globs(iglob)%nsub))
!        find subdomains in the glob         
         ! CORE PART OF CLASSIFICATION FOLLOWS
         globs(iglob)%subdomains(:) = sublist(globs(iglob)%nodes(1),1:globs(iglob)%nsub)
         if (globs(iglob)%nsub.eq.2) then
!           such entity is a face
            globs(iglob)%itype = 1
         else if (globs(iglob)%nsub.gt.2.and.globs(iglob)%nnod.gt.1) then
!           such entity is an edge
            globs(iglob)%itype = 2
         else if (globs(iglob)%nsub.gt.2.and.globs(iglob)%nnod.eq.1) then
!           such entity is a vertex
            globs(iglob)%itype = 3
         else
            call warning('GETGLOBS','Unable to determine type of glob.')
         end if
      end do
      deallocate(sublist)
      deallocate(nodeglob)

      nfaces  = count(globs%itype.eq.1)
      nedges  = count(globs%itype.eq.2)
      nvertex = count(globs%itype.eq.3)
      if (nfaces + nedges + nvertex .ne. nglobs) then
         call error('GETGLOBS','Total number of globs does not match.')
      end if

      if (debug) then
         write(*,*) 'Situation in globs after equivalence class algorithm:'
         write(*,'(a,i8)') 'nfaces   = ',nfaces
         write(*,'(a,i8)') 'nedges   = ',nedges
         write(*,'(a,i8)') 'nvertex  = ',nvertex
         write(*,'(a)')    '--------------------'
         write(*,'(a,i8)') 'total    = ',nfaces + nedges + nvertex
      end if

      if (debug) then
         do iglob = 1,nglobs
            write(*,*) 'glob=============', iglob
            write(*,*) globs(iglob)%itype
            write(*,*) globs(iglob)%nnod
            write(*,*) globs(iglob)%nodes
            write(*,*) globs(iglob)%nsub
            write(*,*) globs(iglob)%subdomains
         end do
      end if


! Prepare space for new vertices
      lnewvertex = nnodi
      allocate(newvertex(lnewvertex))
      newvertex = 0

! mark current vertices into the newvertex array
      do iglob = 1,nglobs
         if (globs(iglob)%itype .eq. 3) then
            newvertex(globs(iglob)%nodes) = 1
         end if
      end do

! generate list of neighbouring subdomains for further algorithm
      lsubneib = nsub
      allocate(subneib(lsubneib))
      lneighbourings = nsub
      allocate(neighbourings(lneighbourings))
      do isub = 1,nsub
         subneib = 0
         do iglob = 1,nglobs
            if (meshdim.eq.2 .and. ndim .eq. 3) then
               ! faces and edges are interesting for shell elements
               if ((globs(iglob)%itype .eq. 1 .or. globs(iglob)%itype .eq. 2) .and. any(globs(iglob)%subdomains .eq. isub)) then
                  subneib(globs(iglob)%subdomains) = 1
               end if
            else 
               ! only faces are interesting in other cases
               if (globs(iglob)%itype .eq. 1 .and. any(globs(iglob)%subdomains .eq. isub)) then
                  subneib(globs(iglob)%subdomains) = 1
               end if
            end if
         end do
         ! my subdomain is also marked
         nnsub = count(subneib.eq.1) - 1
         neighbourings(isub)%nnsub = nnsub
         allocate(neighbourings(isub)%list(nnsub))
         iisub = 0
         do jsub = 1,nsub
            if (subneib(jsub).eq.1 .and. jsub.ne.isub) then
               iisub = iisub + 1
               neighbourings(isub)%list(iisub) = jsub
            end if
         end do
      end do
      if (debug) then
         do isub = 1,nsub
            write(*,*) 'isub',isub,'neighbours',neighbourings(isub)%list
         end do
      end if

      ! coordinates vector for nodes
      lxyzbase = ndim
      allocate (xyzbase(lxyzbase))

      if (debug) then
         write(*,*) '   correction of sufficient number of corners on each subdomain...'
      end if
      do isub = 1,nsub
! Creation of field kmynodes(nnod) - local usage of global nodes
! if global node is in my subdomain, assigned 1, else remains 0
         call pp_mark_sub_nodes(isub,nelem,iets,liets,inet,linet,nnet,lnnet,&
                                kinet,lkinet,kmynodes,lkmynodes)
         do jsub = 1,neighbourings(isub)%nnsub
            indjsub = neighbourings(isub)%list(jsub)
            kinodes = 0
! Creation of field kneibnodes(nnod) - local usage of global nodes
! if global node is in neighbouring subdomain, assigned 1, else remains 0
            call pp_mark_sub_nodes(indjsub,nelem,iets,liets,inet,linet,nnet,lnnet,&
                                   kinet,lkinet,kneibnodes,lkneibnodes)
            where( kmynodes.eq.1 .and. kneibnodes.eq.1 ) kinodes = 1
            nshared = count(kinodes.eq.1)
! number of corners already shared by the two subdomains
            ncornerspair = 0
            do inodi = 1,nnodi
               if (newvertex(inodi).gt.0 .and. kinodes(igingn(inodi)).eq.1) then
                  ncornerspair = ncornerspair + 1
               end if
            end do
            if (debug) then
               write (*,*) 'Subdomains',isub,' and ',indjsub,'already share ',ncornerspair,' corners.'
            end if
            ! generate at least 3 corners on the common face
            ! prepare field of local interface nodes
            lisingin = nshared
            allocate(isingin(lisingin))
            indisn = 0
            do inodi = 1,nnodi
               inod = igingn(inodi)
               if (kinodes(inod).eq.1) then
                  indisn = indisn + 1
                  isingin(indisn) = inodi
               end if
            end do
            ! Check you have filled the whole array
            if (indisn.ne.nshared) then
               write(*,*) 'Error in construction of subdomain interface.'
               stop
            end if
            ! here, there would be different size if the common interface is in fact multicomponent
            lplayground = nshared
            allocate(playground(lplayground))
            playground = 0
            ! Local coordinates
            lxyzsh1 = nshared
            lxyzsh2 = ndim
            ! localize coordinates of shared nodes
            allocate(xyzsh(lxyzsh1,lxyzsh2))
            xyzsh = xyz(igingn(isingin),:)

            ! mark already known corners in playground
            where (newvertex(isingin) .eq. 1) playground = 1
            if (nshared.gt.0) then
               ! no corner on subdomain yet
               ! make it the most remote node to the first interface node
               inodsh = 1

               xyzbase = xyzsh(inodsh,:)
               
               ! Find second corner by maximizing the distance of the first interface node
               ldist = nshared
               allocate(dist(ldist))
               do ish = 1,nshared
                  dist(ish) = sum((xyzsh(ish,:) - xyzbase(:))**2)
               end do
            !   dist = (xis-xbase)**2 + (yis-ybase)**2 + (zis-zbase)**2
               indaux  = maxloc(dist)
               inodcs1 = indaux(1)

               playground(inodcs1) = 1
               ncornerspair = 1

               deallocate(dist)
            end if
            if (ndim.gt.1 .and. meshdim.gt.1 .and. nshared.gt.1) then
               ! one corner is already set in playground, select the second
               inodcs = 0
               do inodis = 1,nshared
                  if (playground(inodis).eq.1) then
                     inodcs = inodis
                     exit
                  end if
               end do
               if (inodcs.eq.0) then
                  call error('PP_GETGLOBS','Problem finding already known corners in playground.')
               end if
               xyzbase = xyzsh(inodcs,:)
               
               ! Find second corner by maximizing the distance from the first one
               ldist = nshared
               allocate(dist(ldist))
               do ish = 1,nshared
                  dist(ish) = sum((xyzsh(ish,:) - xyzbase(:))**2)
               end do
!               dist = (xis-xbase)**2 + (yis-ybase)**2 + (zis-zbase)**2
               indaux  = maxloc(dist)
               inodcs2 = indaux(1)
               if (inodcs.eq.inodcs2) then
                  write(*,*) 'x,y,z'
                  write(*,*) xyzsh(ish,:)
                  write(*,*) 'Problem finding second corner on subdomain - same as first.'
                  write(*,*) 'dist'
                  write(*,*) dist
                  write(*,*) 'first corner', inodcs, ', proposed second corner',inodcs2
                  stop
               end if

               playground(inodcs2) = 1
               ncornerspair = 2

               deallocate(dist)
            end if
            ! it can not happen that there are only two nodes shared by two subdomains for quadratic elements
            !if (ndim.gt.2.and.(ncornerspair.eq.2.or..not.minimizecorners).and.meshdim.gt.2 .and. nshared.gt.2) then
            if (ndim.gt.2 .and. meshdim.gt.2 .and. nshared.gt.2) then
               ! two corners are already set in playground, select the third
               inodcs1 = 0
               do inodis = 1,nshared
                  if (playground(inodis).eq.1) then
                     inodcs1 = inodis
                     exit
                  end if
               end do
               if (inodcs1.eq.0) then
                  call error('GETGLOBS','Problem finding already known corners in playground.')
               end if
               inodcs2 = 0
               do inodis = inodcs1+1,nshared
                  if (playground(inodis).eq.1) then
                     inodcs2 = inodis
                     exit
                  end if
               end do
               if (inodcs2.eq.0) then
                  call error('GETGLOBS','Problem finding already known corners in playground.')
               end if
               x1 = xyzsh(inodcs1,1)
               y1 = xyzsh(inodcs1,2)
               z1 = xyzsh(inodcs1,3)
               x2 = xyzsh(inodcs2,1)
               y2 = xyzsh(inodcs2,2)
               z2 = xyzsh(inodcs2,3)
               
               ! Find third corner as the most orthogonal
               larea = nshared
               allocate(area(larea))
!               cosangle = (xis-x1)*(x2-x1) + (yis-y1)*(y2-y1) + (zis-z1)*(z2-z1)
               do ish = 1,nshared
                  xish = xyzsh(ish,1)
                  yish = xyzsh(ish,2)
                  zish = xyzsh(ish,3)
                  area(ish) = ((y2-y1)*(zish-z1) - (yish-y1)*(z2-z1))**2 &
                            + ((x2-x1)*(zish-z1) - (xish-x1)*(z2-z1))**2 &
                            + ((x2-x1)*(yish-y1) - (xish-x1)*(y2-y1))**2 
               end do
               indaux  = maxloc(abs(area))
               inodcs3 = indaux(1)
               if (inodcs1.eq.inodcs3.or.inodcs2.eq.inodcs3) then
                  write(*,*) 'WARNING: Problem finding third corner between subdomain ',isub,' and',jsub,&
                             ' - same as one of previous two.'
                  write(*,*) 'Already have:',inodcs1,inodcs2,' proposed third corner is', inodcs3
                  write(*,*) 'area'
                  do i = 1,larea
                     write(*,*) i, area(i)
                  end do
                  !stop
               else
                  playground(inodcs3) = 1
                  ncornerspair = 3
               end if

               deallocate(area)
            end if
            deallocate(xyzsh)

            ! now copy playground to global field of new vertices
            newvertex(isingin) = playground
            deallocate(playground)
            deallocate(isingin)
         end do
      end do

      ncorners = count(newvertex.eq.1)
      if (ncorners .lt. ncornersmin) then
! Set random corners form interface
      write(*,'(a,$)') '  Larger number of corners demanded than already selected - adding random corners...'
         do inc = ncorners+1,ncornersmin
  50        call get_random_number(rnd)
            indi = int(rnd*nnodi) + 1
            if (newvertex(indi).eq.1) goto 50
            newvertex(indi) = 1
         end do
         write(*,'(a)') '...done.'
      end if
      ncorners = count(newvertex.eq.1)

! Check subdomain by subdomain that the number of vertices is sufficient
      do isub = 1,nsub
! Creation of field kmynodes(nnod) - local usage of global nodes
! if global node is in my subdomain, assigned 1, else remains 0
         call pp_mark_sub_nodes(isub,nelem,iets,liets,inet,linet,nnet,lnnet,&
                                kinet,lkinet,kmynodes,lkmynodes)
         ncornerss = 0
         do inodi = 1,nnodi
            inod = igingn(inodi)
            if (kmynodes(inod).eq.1.and.newvertex(inodi).eq.1) then
               ncornerss = ncornerss + 1
            end if
         end do
         if (debug) then
            write(*,*) 'Number of corners on subdomain',isub,' is ',ncornerss
         end if
         ! if number of corners is lower than ndim, add some corners
         if (ncornerss.lt.ndim) then
            call warning('GETGLOBS','Number of corners on subdomain lower than dimension.')
            write(*,*) 'Number of corners on subdomain',isub,' is ',ncornerss
            nnodis = count(kmynodes.eq.1.and.kinodes.gt.1)
            write(*,*) 'Number of interface nodes on subdomain',isub,' is ',nnodis
            if (nnodis.eq.0) then
               if (any(iets.eq.isub)) then
                  write(*,*) 'Error in topology - subdomain ',isub,' has zero interface nodes.'
                  call error_exit
               end if
            end if
         end if
      end do
      deallocate(kinet)
      deallocate(kmynodes,kinodes,kneibnodes)

! remove vertices from list of globs - will be used more times
      do inodi = 1,nnodi
         if (newvertex(inodi).eq.1) then
            do iglob = 1,nglobs
               where (globs(iglob)%nodes.eq.inodi) globs(iglob)%nodes = 0
            end do
         end if
      end do

! Creation of field KDOF(NNOD) with addresses before first global
! dof of node
      lkdof = nnod
      allocate(kdof(lkdof))
      if (lkdof.gt.0) then
         kdof(1) = 0
      end if
      do inod = 2,nnod
         kdof(inod) = kdof(inod-1) + nndf(inod-1)
      end do
      nremoved_corners_bc = 0
      if (remove_bc_nodes) then
! remove Dirichlet boundary conditions from the list of globs
         do inodi = 1,nnodi
            indnod = igingn(inodi)
            inddof = kdof(indnod)
            ndofn = nndf(indnod)
            if (any(ifix(inddof+1:inddof+ndofn).gt.0)) then
               ! remove node from globs
               do iglob = 1,nglobs
                  where (globs(iglob)%nodes.eq.inodi) globs(iglob)%nodes = 0
               end do
               ! remove node from corners
   ! uncomment this if you allow corners on Dirichlet BC
   !            if (newvertex(inodi).gt.0) then
   !               newvertex(inodi) = 0
   !               nremoved_corners_bc = nremoved_corners_bc + 1
   !            end if
            end if
         end do

         if (debug) then
            write(*,*) 'Number of removals of corners for Dirichlet BC:', nremoved_corners_bc
         end if
      end if

! remove zeros from list of globs
      ncorrections = 0
      do iglob = 1,nglobs
         nng = globs(iglob)%nnod
         inewnodes = 0
         do ing = 1,nng
            indinode = globs(iglob)%nodes(ing)
            if (indinode.ne.0) then
               inewnodes = inewnodes + 1
               globs(iglob)%nodes(inewnodes) = indinode
            else
               ncorrections = ncorrections + 1
            end if
         end do
         globs(iglob)%nodes(inewnodes+1:) = 0
         globs(iglob)%nnod = inewnodes
      end do
      if (debug) then
         write(*,*) 'Number of removals of nodes from nodes of globs:',ncorrections
      end if

! Mark degenerated globs with negative integer in the field glob_type
      if (debug) then
         write(*,*) 'I am going to remove ', count(globs%nnod.eq.0.and.globs%itype.ne.-2), &
                    ' degenerated globs (including all vertices).'
      end if
      where (globs%nnod.eq.0) globs%itype = -2
      where (globs%nnod.eq.0) globs%selected = .false.

! All vertices should be removed by now - check that
      do iglob = 1,nglobs
         if (globs(iglob)%itype.eq.3.and.globs(iglob)%nnod.ne.0) then
            write(*,*) 'Strange vertex!', iglob
            write(*,*) 'Number of nodes =', globs(iglob)%nnod
            write(*,*) 'Nodes'
            inodi = globs(iglob)%nodes(1)
            write(*,*) inodi
            write(*,*) 'Number of interface nodes matching',count(igingn.eq.igingn(inodi))
         end if
      end do

      nfaces   = count(globs%itype.eq.1)
      nedges   = count(globs%itype.eq.2)
      nvertex  = count(globs%itype.eq.3)
      write(*,*) '   Situation in globs after corrections:'
      write(*,'(a,i8)') '    nfaces  = ',nfaces
      write(*,'(a,i8)') '    nedges  = ',nedges
      write(*,'(a,i8)') '    nvertex = ',nvertex
      write(*,'(a)')    '    --------------------'
      write(*,'(a,i8)') '    total   = ',nfaces + nedges + nvertex
      write(*,'(a)')    '    --------------------'
      write(*,'(a,i8)') '    corners = ',count(newvertex.eq.1)

! Perform check of interface coverage - union of faces, edges, corners and Dirichlet BC should form the whole interface 
! moreover, faces, edges and corners should be disjoint
      licheck = nnodi
      allocate(icheck(licheck))
      icheck = 0
      ! mark corners
      icheck = newvertex
      ! go through globs
      do iglob = 1,nglobs
         if (globs(iglob)%itype.eq.1.or.globs(iglob)%itype.eq.2) then
            nglobn = globs(iglob)%nnod
            icheck(globs(iglob)%nodes(1:nglobn)) = icheck(globs(iglob)%nodes(1:nglobn)) + 1
         end if
      end do
! mark Dirichlet boundary conditions from the list of globs
      if (remove_bc_nodes) then
         do inodi = 1,nnodi
            indnod = igingn(inodi)
            inddof = kdof(indnod)
            ndofn = nndf(indnod)
            if (any(ifix(inddof+1:inddof+ndofn).gt.0)) then
               icheck(inodi) = 1
            end if
         end do
      end if
      ! Every node should be present once
      if (any(icheck.ne.1)) then
         write(*,*) 'Error in interface coverage:'
         do inodi = 1,nnodi
            indnod = igingn(inodi)
            point = kdof(indnod)
            ndofn = nndf(indnod)
            if (icheck(inodi).ne.1) then
               write(*,*) 'Node ',igingn(inodi),' present ',icheck(inodi),' times.'
            end if
         end do
         stop
      end if
      deallocate(icheck)
      if (debug) then
         write(*,*) '     coverage of interface by globs successfully checked.'
      end if
      deallocate(kdof)

! Produce outputs
      nnewvertex = count(newvertex.eq.1)
      if (nnewvertex.ne.ncorners) then
         call error('PP_GET_GLOBS','Total number of corners does not match.')
      end if
      kglobs(:) = 0
      typeglobs(:) = 0

! order globs as corners, edges, faces
      inodc = 0 ! pseudo coarse nodes (includes corners, edges and faces)

      ! export corners
      do i = 1,nnodi
         if (newvertex(i).eq.1) then
            inodc = inodc + 1
            indnodeg = igingn(i)
            kglobs(indnodeg)    = inodc
            typeglobs(indnodeg) = 3
         end if
      end do
      ! export edges
      nedges = 0
      do iglob = 1,nglobs
         if (globs(iglob)%selected.and.globs(iglob)%itype.eq.2) then
            nedges = nedges + 1
            inodc = inodc + 1
            nng  = globs(iglob)%nnod
            do ing = 1,nng
               indnodeg = igingn(globs(iglob)%nodes(ing))
               kglobs(indnodeg)    = inodc
               typeglobs(indnodeg) = globs(iglob)%itype
            end do
         end if
      end do
      ! export faces
      nfaces = 0
      do iglob = 1,nglobs
         if (globs(iglob)%selected.and.globs(iglob)%itype.eq.1) then
            nfaces = nfaces + 1
            inodc = inodc + 1
            nng  = globs(iglob)%nnod
            do ing = 1,nng
               indnodeg = igingn(globs(iglob)%nodes(ing))
               kglobs(indnodeg)    = inodc
               typeglobs(indnodeg) = globs(iglob)%itype
            end do
         end if
      end do

! Clear memory
      do iglob = 1,nglobs
         deallocate(globs(iglob)%nodes)
         deallocate(globs(iglob)%subdomains)
      end do
      deallocate(globs)
      do isub = 1,nsub
         deallocate(neighbourings(isub)%list)
      end do
      deallocate(neighbourings)
      deallocate(newvertex)
      deallocate(subneib)
      deallocate(xyzbase)
      deallocate(igingn)

end subroutine pp_get_globs

!***********************************************************************
subroutine pp_mark_sub_nodes(isub,nelem,iets,liets,inet,linet,nnet,lnnet,&
                             kinet,lkinet,knodes,lknodes)
!***********************************************************************
!     Marks nodes of subdomain ISUB in array KNODES
!***********************************************************************
implicit none
      
! subdomain index
integer,intent(in):: isub 

! global number of elements
integer,intent(in):: nelem 

! indices of elements in subdomains
integer,intent(in):: liets 
integer,intent(in)::  iets(liets) 

! indices of nodes on elements
integer,intent(in):: linet 
integer,intent(in)::  inet(linet) 

! number of nodes on elements
integer,intent(in):: lnnet 
integer,intent(in)::  nnet(lnnet) 

! key to INET - pointers before element nodes
integer,intent(in):: lkinet 
integer,intent(in)::  kinet(lkinet) 

! key nodes - marked nodes for subdomain ISUB
integer,intent(in)::  lknodes 
integer,intent(out)::  knodes(lknodes) 

! local variables
integer:: ie, indsub, ine, nne, ipoint, indng

      knodes = 0
      do ie = 1,nelem
         indsub = iets(ie)
         if (indsub.eq.isub) then
            nne = nnet(ie)
            ipoint = kinet(ie)
            do ine = 1,nne
               indng = inet(ipoint + ine)
               knodes(indng) = 1
            end do
         end if
      end do

end subroutine pp_mark_sub_nodes

!****************************************************************
subroutine pp_distribute_linearly(nsub,nproc, sub2proc,lsub2proc)
!****************************************************************
!     Distribute linearly NSUB elements to NPROC parts
!****************************************************************
implicit none
      
! INPUT:
! subdomain number
integer,intent(in):: nsub 
! number of processors
integer,intent(in):: nproc 
! OUTPUT:
! sub2proc(nproc + 1) SUB2PROC array (ParMETIS-like)
integer,intent(in)::  lsub2proc 
integer,intent(out)::  sub2proc(lsub2proc) 

! local vars
integer :: nsub_locx, nsub_loc, iproc

nsub_locx = (nsub + nproc - 1)/nproc
sub2proc(1) = 1
do iproc = 1,nproc
   nsub_loc  = max(min(nsub - (iproc-1)*nsub_locx, nsub_locx),0)
   sub2proc(iproc+1) = sub2proc(iproc) + nsub_loc
end do

end subroutine pp_distribute_linearly

!*****************************************************************
subroutine pp_get_proc_for_sub(isub,comm,sub2proc,lsub2proc,iproc)
!*****************************************************************
!     Get processor id  for subdomain isub
!*****************************************************************
use module_utils
implicit none
include "mpif.h"
      
! INPUT:
! subdomain index
integer,intent(in):: isub 
! MPI communicator
integer,intent(in):: comm 
! sub2proc(nproc + 1) SUB2PROC array (ParMETIS-like)
integer,intent(in):: lsub2proc 
integer,intent(in)::  sub2proc(lsub2proc) 
! OUTPUT:
! index of processor who has the subdomain
integer,intent(out):: iproc 

! local vars
integer :: nproc, ierr, ip

call MPI_COMM_SIZE(comm,nproc,ierr)

! check data
if (nproc+1.ne.lsub2proc) then
   call error('PP_GET_PROC_FOR_SUB','lsub2proc length error')
end if
if (debug) then
   do ip = 2,nproc+1
      if (sub2proc(ip).lt.sub2proc(ip-1)) then
         call error('PP_GET_PROC_FOR_SUB','lsub2proc not monotone')
      end if
   end do
end if

iproc = 0
do ip = 2,lsub2proc
   if (isub.lt.sub2proc(ip)) then
      return
   else
      iproc = iproc + 1
   end if
end do

call error('PP_GET_PROC_FOR_SUB','processor not found for subdomain',isub)

end subroutine pp_get_proc_for_sub

!******************************************************************
subroutine pp_get_unique_tag(isub,jsub,comm,sub2proc,lsub2proc,tag)
!******************************************************************
! Subroutine for assigning unique tag
implicit none
! global numbers of subdomains
integer, intent(in) :: isub
integer, intent(in) :: jsub
! MPI communicator
integer,intent(in):: comm 
! sub2proc(nproc + 1) SUB2PROC array (ParMETIS-like)
integer,intent(in):: lsub2proc 
integer,intent(in)::  sub2proc(lsub2proc) 

! unique tag
integer, intent(out):: tag

! local variables
integer:: iproc, jproc, nsub_loc_j, isub_loc, jsub_loc

      call pp_get_proc_for_sub(isub,comm,sub2proc,lsub2proc,iproc)
      call pp_get_proc_for_sub(jsub,comm,sub2proc,lsub2proc,jproc)
      
      isub_loc = isub - sub2proc(iproc+1) + 1
      jsub_loc = jsub - sub2proc(jproc+1) + 1
      
      nsub_loc_j = sub2proc(jproc+2) - sub2proc(jproc+1)
      
      ! compute tag based on local numbers
      tag = isub_loc*nsub_loc_j + jsub_loc

end subroutine

!********************************************************************
subroutine pp_get_nevax(nelem,inet,linet,nnet,lnnet,nndf,lnndf,nevax)
!********************************************************************
implicit none

integer,intent(in) :: nelem
integer,intent(in) :: linet
integer,intent(in) ::  inet(linet)
integer,intent(in) :: lnnet
integer,intent(in) ::  nnet(lnnet)
integer,intent(in) :: lnndf
integer,intent(in) ::  nndf(lnndf)

integer,intent(out) ::  nevax ! number of variables on element

! Local variables
integer :: ie, nevab, indinet, ine, nne, indn

      nevax = 0
      indinet = 0
      do ie = 1, nelem
         nevab = 0
         nne = nnet(ie)
         do ine = 1,nne
            indinet = indinet + 1
            indn = inet(indinet)
            nevab = nevab + nndf(indn)
         end do
         if (nevab.gt.nevax) then
            nevax = nevab
         end if
     end do
end subroutine pp_get_nevax

end module module_pp
