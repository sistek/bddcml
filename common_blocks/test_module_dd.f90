program test_module_dd
! Tester of module_dd
      use module_pp
      use module_dd
      use module_utils
      implicit none
      include "mpif.h"

      integer,parameter :: kr = kind(1.D0)

      integer :: matrixtype

      !  parallel variables
      integer :: myid, comm_all, comm_self, nproc, ierr

      integer ::             lsub2proc
      integer,allocatable ::  sub2proc(:)

      integer :: nsub, isub

      ! only works for this problem
      character(7) :: problemname = 'TESTLAP'

      integer,parameter :: lx = 2
      real(kr) :: x(lx)
      integer,parameter :: ly = 2
      real(kr) :: y(ly)
      integer,parameter :: lz = 4
      real(kr) :: z(lz)

      integer idproc

      integer :: i

      logical :: debug = .true.
      logical :: remove_original 

      ! MPI initialization
!***************************************************************PARALLEL
      call MPI_INIT(ierr)
      ! Communicator
      comm_all  = MPI_COMM_WORLD
      comm_self = MPI_COMM_SELF
      call MPI_COMM_RANK(comm_all,myid,ierr)
      call MPI_COMM_SIZE(comm_all,nproc,ierr)
!***************************************************************PARALLEL

      print *, 'I am processor ',myid,': Hello nproc =',nproc
      ! SPD matrix
      matrixtype = 1

      nsub = 2
      call dd_init(nsub)

      print *, 'I am processor ',myid,': nproc = ',nproc, 'nsub = ',nsub
      lsub2proc = nproc + 1
      allocate(sub2proc(lsub2proc))
      call pp_distribute_subdomains(nsub,nproc,sub2proc,lsub2proc)
      call dd_distribute_subdomains(1,nsub,sub2proc,lsub2proc,nproc)
      deallocate(sub2proc)


      do isub = 1,nsub
         ! load mesh
         call dd_read_mesh_from_file(myid,isub,problemname)

         ! load matrices into the structure
         call dd_read_matrix_from_file(myid,isub,problemname,matrixtype)

         ! assembly matrices
         call dd_assembly_local_matrix(myid,isub)

         ! convert matrix to blocks
         remove_original = .false.
         call dd_matrix_tri2blocktri(myid,isub,remove_original)

         ! prepare Schur complements
         call dd_prepare_schur(myid,comm_self,isub)
      end do

      print *, 'I am here 1'
      call flush(6)

      ! create neigbouring 
      call dd_create_neighbouring(myid,nsub,comm_all)

      ! prepare weights
      call dd_weights_prepare(myid, nsub, comm_all)

      ! prepare reduced RHS
      call dd_prepare_reduced_rhs(myid,nsub,comm_all)

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

!      ! auxiliary routine, until reading directly the globs
!      isub = 1
!      call dd_get_cnodes(myid,isub)
!      isub = 2
!      call dd_get_cnodes(myid,isub)
!
!      isub = 1
!      call dd_prepare_c(myid,isub)
!      isub = 2
!      call dd_prepare_c(myid,isub)
!
!      ! prepare augmented matrix for BDDC
!      do isub = 1,nsub
!         call dd_prepare_aug(myid,comm_self,isub)
!      end do
!
!      ! prepare coarse space basis functions for BDDC
!      do isub = 1,nsub
!         call dd_prepare_coarse(myid,isub)
!      end do
!

      print *, 'I am here 2'
      call flush(6)

!      ! test selection to interface
      write(*,*) 'I am ',myid,'Mapping subdomain variables to subdomain interface vars.'
      do i = 1,4
         z(i) = float(i)
      end do
      isub = 1
      x = -1._kr
      call dd_map_sub_to_subi(myid,isub, z,lz, x,lx)
      write(*,*) 'I am ',myid,'z =',z
      write(*,*) 'I am ',myid,'x =',x
      do i = 1,4
         z(i) = float(i)
      end do
      isub = 2
      x = -1._kr
      call dd_map_sub_to_subi(myid,isub, z,lz, x,lx)
      write(*,*) 'I am ',myid,'z =',z
      write(*,*) 'I am ',myid,'x =',x


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
