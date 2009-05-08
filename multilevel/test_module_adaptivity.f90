program test_module_adaptivity
! Tester of module_adaptivity
      use module_adaptivity
      use module_dd
      use module_utils

      implicit none
      include "mpif.h"

      integer,parameter :: kr = kind(1.D0)

      ! parallel variables
      integer :: myid, comm, nproc, ierr

      integer :: idpar, idpair 
      integer :: npair

      integer :: matrixtype
      ! how many pairs are assigned to a processor
      integer :: npair_locx
      integer :: ndim, nsub, nelem, ndof, nnod, nnodc, linet

      integer :: isub
      logical :: remove_original 

      integer :: i,j, ndofi, nnodi
      integer :: lschur
      real(kr),allocatable :: schur(:)

      character(90)  :: problemname 
      character(100) :: name
      character(100) :: filename

      logical :: debug = .true.

      ! MPI initialization
!***************************************************************PARALLEL
      call MPI_INIT(ierr)
      ! Communicator
      comm = MPI_COMM_WORLD
      call MPI_COMM_RANK(comm,myid,ierr)
      call MPI_COMM_SIZE(comm,nproc,ierr)
!***************************************************************PARALLEL

! Initial screen
      if (myid.eq.0) then
         write(*,'(a)') 'ADAPTIVITY TESTER'
         write(*,'(a)') '================='

! Name of the problem
   10    write(*,'(a,$)') 'Name of the problem: '
         read(*,*) problemname
         if(problemname.eq.' ') goto 10

      end if
! Broadcast of name of the problem      
!***************************************************************PARALLEL
      call MPI_BCAST(problemname, 90, MPI_CHARACTER, 0, comm, ierr)
!***************************************************************PARALLEL

      if (myid.eq.0) then
         name = trim(problemname)//'.PAR'
         call allocate_unit(idpar)
         open (unit=idpar,file=name,status='old',form='formatted')
      end if

! Reading basic properties 
      if (myid.eq.0) then
         read(idpar,*) ndim, nsub, nelem, ndof, nnod, nnodc, linet
      end if
! Broadcast basic properties of the problem
!***************************************************************PARALLEL
      call MPI_BCAST(ndim,     1, MPI_INTEGER,         0, comm, ierr)
      call MPI_BCAST(nsub,     1, MPI_INTEGER,         0, comm, ierr)
      call MPI_BCAST(nelem,    1, MPI_INTEGER,         0, comm, ierr)
      call MPI_BCAST(ndof,     1, MPI_INTEGER,         0, comm, ierr)
      call MPI_BCAST(nnod,     1, MPI_INTEGER,         0, comm, ierr)
      call MPI_BCAST(nnodc,    1, MPI_INTEGER,         0, comm, ierr)
      call MPI_BCAST(linet,    1, MPI_INTEGER,         0, comm, ierr)
!***************************************************************PARALLEL



      ! open file with description of pairs
      if (myid.eq.0) then
         filename = trim(problemname)//'.PAIR'
         call allocate_unit(idpair)
         open (unit=idpair,file=filename,status='old',form='formatted')
      end if
   
      print *, 'I am processor ',myid,': Hello nproc =',nproc
      ! SPD matrix
      matrixtype = 1

      call dd_init(nsub)
      call dd_distribute_subdomains(nsub,nproc)
      call dd_read_mesh_from_file(myid,trim(problemname))
      call dd_read_matrix_from_file(myid,trim(problemname),matrixtype)
      call dd_assembly_local_matrix(myid)
      remove_original = .false.
      call dd_matrix_tri2blocktri(myid,remove_original)
      do isub = 1,nsub
         call dd_prepare_schur(myid,comm,isub)
      end do
      call dd_prepare_adaptive_space(myid)

!      do isub = 1,nsub
!         call dd_get_interface_size(myid,isub,ndofi,nnodi)
!         lschur = ndofi*ndofi
!         allocate (schur(lschur))
!         call dd_get_schur(myid, isub, schur,lschur)
!         print *,'Schur complement matrix. isub = ',isub
!         do i = 1,ndofi
!            print '(100f10.6)',(schur((j-1)*ndofi + i),j = 1,ndofi)
!         end do
!         deallocate(schur)
!      end do
   
      call adaptivity_init(myid,comm,idpair,npair)

      print *, 'I am processor ',myid,': nproc = ',nproc, 'nsub = ',nsub
      call adaptivity_assign_pairs(npair,nproc,npair_locx)

      if (debug) then
         call adaptivity_print_pairs(myid)
      end if

      call adaptivity_solve_eigenvectors(myid,comm,npair_locx,npair,nproc)

      call adaptivity_finalize
   
      call dd_finalize
      print *, 'I am processor ',myid,': Hello 2'

      ! close file with description of pairs
      if (myid.eq.0) then
         close(idpair)
      end if

      ! MPI finalization
!***************************************************************PARALLEL
      call MPI_FINALIZE(ierr)
!***************************************************************PARALLEL

end program
