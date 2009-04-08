program test_module_adaptivity
! Tester of module_adaptivity
      use module_adaptivity
      use module_utils

      implicit none
      include "mpif.h"

      integer,parameter :: kr = kind(1.D0)

      !  parallel variables
      integer :: myid, comm, nproc, ierr

      integer :: idpair
      integer :: npair

      integer :: matrixtype, nsub
      ! how many pairs are assigned to a processor
      integer :: npair_locx

      character(7)  :: problemname = 'TESTLAP'
      character(20) :: filename

      logical :: debug = .true.

      ! MPI initialization
!***************************************************************PARALLEL
      call MPI_INIT(ierr)
      ! Communicator
      comm = MPI_COMM_WORLD
      call MPI_COMM_RANK(comm,myid,ierr)
      call MPI_COMM_SIZE(comm,nproc,ierr)
!***************************************************************PARALLEL

      ! open file with description of pairs
      if (myid.eq.0) then
         filename = problemname//'.PAIR'
         call allocate_unit(idpair)
         open (unit=idpair,file=filename,status='old',form='formatted')
      end if
   
      print *, 'I am processor ',myid,': Hello nproc =',nproc
      ! SPD matrix
      matrixtype = 1

      nsub = 2
      call adaptivity_init(myid,comm,idpair,npair)

      print *, 'I am processor ',myid,': nproc = ',nproc, 'nsub = ',nsub
      call adaptivity_assign_pairs(npair,nproc,npair_locx)

      if (debug) then
         call adaptivity_print_pairs(myid)
      end if

      call adaptivity_finalize
   
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
