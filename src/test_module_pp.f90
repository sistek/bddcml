! BDDCML - Multilevel BDDC
! Copyright (C) The BDDCML Team
! 
! This program is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
! 
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
! 
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <http://www.gnu.org/licenses/>.
! ____________________________________________________________________

program test_module_pp
! Tester of module_pp
      use module_pp
      use module_utils

      implicit none
      include "mpif.h"

      integer,parameter :: kr = kind(1.D0)

      ! parallel variables
      integer :: myid, comm_self, comm_all, nproc, ierr

      integer :: ndim, nsub, nelem, ndof, nnod, nnodc, linet
      integer :: idpar = 1, idgmi, idfvs, ides
      integer :: ncorners, ncornersmin, nfaces, nedges, meshdim, nnodi
      integer :: iproc, isub

      integer ::                    lnnet,   lnndf
      integer,allocatable:: inet(:), nnet(:), nndf(:)
      integer ::           liets
      integer,allocatable:: iets(:)
      integer ::           lifix
      integer,allocatable:: ifix(:)
      integer ::            lxyz1, lxyz2
      real(kr),allocatable:: xyz(:,:)
      integer ::           lkglobs
      integer,allocatable:: kglobs(:)
      integer ::           ltypeglobs
      integer,allocatable:: typeglobs(:)
      integer ::           lsub2proc
      integer,allocatable:: sub2proc(:)

      logical :: remove_bc_nodes

      character(90)  :: problemname 
      character(100) :: name
      character(100) :: filename

      real(kr) :: timeaux, timeaux1, timeaux2

      ! MPI initialization
!***************************************************************PARALLEL
      call MPI_INIT(ierr)
      ! Communicator
      comm_all  = MPI_COMM_WORLD
      comm_self = MPI_COMM_SELF
      call MPI_COMM_RANK(comm_all,myid,ierr)
      call MPI_COMM_SIZE(comm_all,nproc,ierr)
!***************************************************************PARALLEL

! Initial screen
      if (myid.eq.0) then
         write(*,'(a)') 'PP MODULE TESTER'
         write(*,'(a)') '================'

! Name of the problem
   10    write(*,'(a,$)') 'Name of the problem: '
         call flush(6)
         read(*,*) problemname
         if(problemname.eq.' ') goto 10

      end if
! Broadcast of name of the problem      
!***************************************************************PARALLEL
      call MPI_BCAST(problemname, 90, MPI_CHARACTER, 0, comm_all, ierr)
!***************************************************************PARALLEL

      if (myid.eq.0) then
         name = trim(problemname)//'.PAR'
         call allocate_unit(idpar)
         open (unit=idpar,file=name,status='old',form='formatted')
      end if

! Reading basic properties 
      if (myid.eq.0) then
         read(idpar,*) ndim, nsub, nelem, ndof, nnod, nnodc, linet
         read(idpar,*) meshdim
      end if
! Broadcast basic properties of the problem
!***************************************************************PARALLEL
      call MPI_BCAST(ndim,     1, MPI_INTEGER,         0, comm_all, ierr)
      call MPI_BCAST(nsub,     1, MPI_INTEGER,         0, comm_all, ierr)
      call MPI_BCAST(nelem,    1, MPI_INTEGER,         0, comm_all, ierr)
      call MPI_BCAST(ndof,     1, MPI_INTEGER,         0, comm_all, ierr)
      call MPI_BCAST(nnod,     1, MPI_INTEGER,         0, comm_all, ierr)
      call MPI_BCAST(nnodc,    1, MPI_INTEGER,         0, comm_all, ierr)
      call MPI_BCAST(linet,    1, MPI_INTEGER,         0, comm_all, ierr)
      call MPI_BCAST(meshdim,  1, MPI_INTEGER,         0, comm_all, ierr)
!***************************************************************PARALLEL

      ! open file with description of MESH
      if (myid.eq.0) then
         filename = trim(problemname)//'.GMIS'
         call allocate_unit(idgmi)
         open (unit=idgmi,file=filename,status='old',form='formatted')
         filename = trim(problemname)//'.FVS'
         call allocate_unit(idfvs)
         open (unit=idfvs,file=filename,status='old',form='formatted')
         filename = trim(problemname)//'.ES'
         call allocate_unit(ides)
         open (unit=ides,file=filename,status='old',form='formatted')
      end if
      linet = linet
      lnnet = nelem
      lnndf = nnod
      lxyz1 = nnod
      lxyz2 = ndim
      allocate(inet(linet),nnet(lnnet),nndf(lnndf),xyz(lxyz1,lxyz2))
      lifix = ndof
      allocate(ifix(lifix))
      liets = nelem
      allocate(iets(liets))
! read fields INET and NNET from GMIS file
      if (myid.eq.0) then
         read(idgmi,*) inet
         read(idgmi,*) nnet
         read(idgmi,*) nndf
         read(idgmi,*) xyz
         close(idgmi)
         read(idfvs,*) ifix
         close(idfvs)
         read(ides,*) iets
         close(ides)
      end if
!***************************************************************PARALLEL
      call MPI_BCAST(inet,  linet,MPI_INTEGER,              0, comm_all, ierr)
      call MPI_BCAST(nnet,  lnnet,MPI_INTEGER,              0, comm_all, ierr)
      call MPI_BCAST(nndf,  lnndf,MPI_INTEGER,              0, comm_all, ierr)
      call MPI_BCAST(xyz, lxyz1*lxyz2,MPI_DOUBLE_PRECISION, 0, comm_all, ierr)
      call MPI_BCAST(ifix,  lifix,MPI_INTEGER,              0, comm_all, ierr)
      call MPI_BCAST(iets,  liets,MPI_INTEGER,              0, comm_all, ierr)
!***************************************************************PARALLEL
   
!*******************************************AUX
! Measure time spent in DD module
      call MPI_BARRIER(comm_all,ierr)
      timeaux1 = MPI_WTIME()
!*******************************************AUX

      remove_bc_nodes = .true.
      ncornersmin = 0
      lkglobs    = nnod
      ltypeglobs = nnod
      allocate(kglobs(lkglobs),typeglobs(ltypeglobs))
      call pp_get_globs(ndim,meshdim,nelem,nnod,nsub,&
                        inet,linet,nnet,lnnet,nndf,lnndf,xyz,lxyz1,lxyz2,&
                        remove_bc_nodes,ifix,lifix,iets,liets,ncornersmin,&
                        nnodi, ncorners,nedges,nfaces,&
                        kglobs,lkglobs, typeglobs,ltypeglobs)
      write(*,*) 'myid = ',myid,'resulting vectors'
      write(*,*) 'myid = ',myid,'kglobs',kglobs
      write(*,*) 'myid = ',myid,'typeglobs',typeglobs

      deallocate(kglobs,typeglobs)
      deallocate(inet,nnet,nndf,xyz)
      deallocate(ifix)
      deallocate(iets)


      nsub = 13
      write(*,*) 'nsub',nsub
      lsub2proc = nproc + 1
      allocate(sub2proc(lsub2proc))
      call pp_distribute_linearly(nsub,nproc, sub2proc,lsub2proc)
      write(*,*) 'sub2proc',sub2proc

      do isub = 1,nsub
         call pp_get_proc_for_sub(isub,comm_all,sub2proc,lsub2proc,iproc)

         write(*,*) 'processor for subdomain ',isub,' is ',iproc
         call flush(6)
      end do
      
      deallocate(sub2proc)

!*******************************************AUX
! Measure time spent in DD module
      call MPI_BARRIER(comm_all,ierr)
      timeaux2 = MPI_WTIME()
      timeaux = timeaux2 - timeaux1
      if (myid.eq.0) then
         write(*,*) '***************************************'
         write(*,*) 'Time used for finding corners and globs  is ',timeaux,' s'
         write(*,*) '***************************************'
      end if

      ! MPI finalization
!***************************************************************PARALLEL
      call MPI_FINALIZE(ierr)
!***************************************************************PARALLEL


end program
