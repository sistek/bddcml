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
      integer ::             lindexsub
      integer,allocatable ::  indexsub(:)

      integer :: nsub, isub, nsub_loc, isub_loc

      ! only works for this problem
      character(7) :: problemname = 'TESTLAP'

      integer,parameter :: lx = 2
      real(kr) :: x(lx)
      integer,parameter :: ly = 2
      real(kr) :: y(ly)
      integer,parameter :: lz = 4
      real(kr) :: z(lz)

      integer ::                         lsubdomains
      type(subdomain_type),allocatable :: subdomains(:)


      integer iproc

      integer :: i

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

      print *, 'I am processor ',myid,': nproc = ',nproc, 'nsub = ',nsub
      lsub2proc = nproc + 1
      allocate(sub2proc(lsub2proc))
      call pp_distribute_linearly(nsub,nproc,sub2proc,lsub2proc)

      nsub_loc = sub2proc(myid+2) - sub2proc(myid+1)
      lsubdomains = nsub_loc
      lindexsub   = nsub_loc
      allocate(subdomains(nsub_loc))
      allocate(indexsub(lindexsub))
      do isub_loc = 1,nsub_loc
         indexsub(isub_loc) = sub2proc(myid+1) + isub_loc - 1
      end do
      do isub_loc = 1,nsub_loc
         isub = indexsub(isub_loc)

         ! initialize subdomain
         call dd_init(subdomains(isub_loc),isub,nsub,comm_all)

         ! load mesh
         call dd_read_mesh_from_file(subdomains(isub_loc),trim(problemname))

         ! load matrices into the structure
         call dd_read_matrix_from_file(subdomains(isub_loc),matrixtype,trim(problemname))

         ! assembly matrices
         call dd_assembly_local_matrix(subdomains(isub_loc))
      end do

      print *, 'I am here 1'
      call flush(6)

      ! create neigbouring 
      call dd_create_neighbouring(subdomains,lsubdomains, sub2proc,lsub2proc, indexsub,lindexsub, comm_all)

      ! prepare weights
      call dd_weights_prepare(subdomains,lsubdomains, sub2proc,lsub2proc, indexsub,lindexsub, comm_all)

      do isub_loc = 1,nsub_loc
         isub = indexsub(isub_loc)

         ! convert matrix to blocks
         remove_original = .false.
         call dd_matrix_tri2blocktri(subdomains(isub_loc),remove_original)

         ! prepare Schur complements
         call dd_prepare_schur(subdomains(isub_loc),comm_self)
      end do

      ! prepare reduced RHS
      call dd_prepare_reduced_rhs_all(subdomains,lsubdomains, sub2proc,lsub2proc, indexsub,lindexsub, comm_all)

      ! multiply vector of ones by subdomain # 1
      x = 1.
      y = 0.
      do isub = 1,nsub
         call pp_get_proc_for_sub(isub,comm_all,sub2proc,lsub2proc,iproc)
         if (myid.eq.iproc) then
            call get_index(isub,indexsub,lindexsub,isub_loc)
            call dd_multiply_by_schur(subdomains(isub_loc),x,lx,y,ly)
            write(*,*) 'I am ',myid,'y =',y
         end if
      end do

      ! auxiliary routine, until reading directly the globs
      !do isub_loc = 1,nsub_loc
      !   ! construct coarse pseudo nodes
      !   call dd_construct_cnodes(subdomains(isub_loc))
      !   ! prepare matrix of constraints C
      !   call dd_prepare_c(subdomains(isub_loc))
      !   ! prepare augmented matrix for BDDC
      !   call dd_prepare_aug(subdomains(isub_loc),comm_self)
      !   ! prepare coarse space basis functions for BDDC
      !   call dd_prepare_coarse(subdomains(isub_loc),.false.)
      !end do

      print *, 'I am here 2'
      call flush(6)

      ! test selection to interface
      write(*,*) 'I am ',myid,'Mapping subdomain variables to subdomain interface vars.'
      do isub = 1,nsub
         do i = 1,4
            z(i) = float(i)
         end do
         x = -1._kr
         call pp_get_proc_for_sub(isub,comm_all,sub2proc,lsub2proc,iproc)
         if (myid.eq.iproc) then
            call get_index(isub,indexsub,lindexsub,isub_loc)

            call dd_map_sub_to_subi(subdomains(isub_loc), z,lz, x,lx)
            write(*,*) 'I am ',myid,'z =',z
            write(*,*) 'I am ',myid,'x =',x
         end if
      end do

      do isub_loc = 1,nsub_loc
         call dd_print_sub(subdomains(isub_loc))
      end do

      do isub_loc = 1,nsub_loc
         call dd_finalize(subdomains(isub_loc))
      end do
      deallocate(subdomains)
      deallocate(sub2proc)
      deallocate(indexsub)
   
      print *, 'I am processor ',myid,': Hello 2'

      ! MPI finalization
!***************************************************************PARALLEL
      call MPI_FINALIZE(ierr)
!***************************************************************PARALLEL

end program
