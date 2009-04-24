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

      integer :: idpair
      integer :: npair

      integer :: matrixtype, nsub
      ! how many pairs are assigned to a processor
      integer :: npair_locx

      integer :: isub
      logical :: remove_original 

      character(7)  :: problemname = 'TESTGLB'
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
      call dd_init(nsub)
      call dd_distribute_subdomains(nsub,nproc)
      call dd_read_mesh_from_file(myid,problemname)
      call dd_read_matrix_from_file(myid,problemname,matrixtype)
      call dd_assembly_local_matrix(myid)
      remove_original = .false.
      call dd_matrix_tri2blocktri(myid,remove_original)
      do isub = 1,nsub
         call dd_prepare_schur(myid,comm,isub)
      end do

   
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
