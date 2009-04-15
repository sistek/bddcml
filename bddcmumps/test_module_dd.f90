program test_module_dd
! Tester of module_dd
      use module_dd
      use module_utils
      implicit none
      include "mpif.h"

      integer,parameter :: kr = kind(1.D0)

      integer :: matrixtype

      !  parallel variables
      integer :: myid, comm, nproc, ierr

      integer :: nsub, isub

      character(7) :: problemname = 'TESTLAP'

      integer,parameter :: lx = 2
      real(kr) :: x(lx)
      integer,parameter :: ly = 2
      real(kr) :: y(ly)

      integer idproc

      logical :: debug = .true.
      logical :: remove_original 

      ! MPI initialization
!***************************************************************PARALLEL
      call MPI_INIT(ierr)
      ! Communicator
      comm = MPI_COMM_WORLD
      call MPI_COMM_RANK(comm,myid,ierr)
      call MPI_COMM_SIZE(comm,nproc,ierr)
!***************************************************************PARALLEL

      print *, 'I am processor ',myid,': Hello nproc =',nproc
      ! SPD matrix
      matrixtype = 1

      nsub = 2
      call dd_init(nsub)

      print *, 'I am processor ',myid,': nproc = ',nproc, 'nsub = ',nsub
      call dd_distribute_subdomains(nsub,nproc)

      ! load mesh
      call dd_read_mesh_from_file(myid,problemname)

      ! load matrices into the structure
      call dd_read_matrix_from_file(myid,problemname,matrixtype)

      ! assembly matrices
      call dd_assembly_local_matrix(myid)

      ! convert matrix to blocks
      remove_original = .false.
      call dd_matrix_tri2blocktri(myid,remove_original)

      ! prepare Schur complements
      do isub = 1,nsub
         call dd_prepare_schur(myid,comm,isub)
      end do

      ! test who owns data for subdomain
      isub = 1
      call dd_where_is_subdomain(isub,idproc)
      write(*,*) 'Subdomain ',isub,' is held on processor ', idproc
      isub = 2
      call dd_where_is_subdomain(isub,idproc)
      write(*,*) 'Subdomain ',isub,' is held on processor ', idproc

      ! multiply vector of ones by subdomain # 1
      x = 1.
      y = 0.
      isub = 1
      call dd_multiply_by_schur(myid,isub,x,lx,y,ly)
      write(*,*) 'I am ',myid,'y =',y
      x = 1.
      y = 0.
      isub = 2
      call dd_multiply_by_schur(myid,isub,x,lx,y,ly)
      write(*,*) 'I am ',myid,'y =',y

      ! test who owns data for subdomain
      isub = 1
      call dd_prepare_c(myid,isub)
      isub = 2
      call dd_prepare_c(myid,isub)

      ! prepare augmented matrix for BDDC
      do isub = 1,nsub
         call dd_prepare_aug(myid,comm,isub)
      end do

      ! prepare coarse space basis functions for BDDC
      do isub = 1,nsub
         call dd_prepare_coarse(myid,isub)
      end do

      if (debug) then
         call dd_print_sub(myid)
      end if

      call dd_finalize
   
      print *, 'I am processor ',myid,': Hello 2'

      ! MPI finalization
!***************************************************************PARALLEL
      call MPI_FINALIZE(ierr)
!***************************************************************PARALLEL

end program
