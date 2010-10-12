program test_module_dd2
! Tester of module_dd
      use module_dd
      use module_utils

      implicit none
      include "mpif.h"

      integer,parameter :: kr = kind(1.D0)

      ! parallel variables
      integer :: myid, comm_self, comm_all, nproc, ierr

      integer :: idpar

      integer :: matrixtype
      integer :: ndim, nsub, nelem, ndof, nnod, nnodc, linet

      integer :: isub
      logical :: remove_original 

      integer :: idlevel
      integer ::             lnndf_coarse
      integer,allocatable ::  nndf_coarse(:)
      integer :: iaux, ncnodes, nedge, nface


      integer :: glob_type

      character(90)  :: problemname 
      character(100) :: name

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
         write(*,'(a)') 'DD MODULE TESTER'
         write(*,'(a)') '================'

! Name of the problem
   10    write(*,'(a,$)') 'Name of the problem: '
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
!***************************************************************PARALLEL

      ! SPD matrix
      matrixtype = 1

!*******************************************AUX
! Measure time spent in DD module
      call MPI_BARRIER(comm_all,ierr)
      timeaux1 = MPI_WTIME()
!*******************************************AUX
      call dd_init(nsub)
      call dd_distribute_subdomains(nsub,nproc)
      call dd_read_mesh_from_file(myid,trim(problemname))
      call dd_read_matrix_from_file(myid,comm_all,trim(problemname),matrixtype)
      call dd_assembly_local_matrix(myid)
      remove_original = .false.
      call dd_matrix_tri2blocktri(myid,remove_original)
      do isub = 1,nsub
         call dd_prepare_schur(myid,comm_self,isub)
      end do


! create coarse mesh 
! read second level
      if (myid.eq.0) then
         name = trim(problemname)//'.L2'
         write(*,*) 'Reading data from file ',trim(name)
         call allocate_unit(idlevel)
         open (unit=idlevel,file=name,status='old',form='formatted')

         read(idlevel,*) iaux, iaux, iaux, iaux
         read(idlevel,*) nnodc, nedge, nface
      end if
!*****************************************************************MPI
      call MPI_BCAST(nnodc,1, MPI_INTEGER, 0, comm_all, ierr)
      call MPI_BCAST(nedge,1, MPI_INTEGER, 0, comm_all, ierr)
      call MPI_BCAST(nface,1, MPI_INTEGER, 0, comm_all, ierr)
!*****************************************************************MPI
       
      ncnodes = nnodc + nedge + nface

      ! prepare array of number of coarse dof in generalized coarse nodes
      lnndf_coarse = ncnodes
      allocate(nndf_coarse(lnndf_coarse))

      nndf_coarse = 0
      ! corners contain ndim coarse dof
      nndf_coarse(1:nnodc) = ndim
      ! edges contain ndim coarse dof
      nndf_coarse(nnodc+1:nnodc+nedge) = ndim


      ! auxiliary routine, until reading directly the globs
      do isub = 1,nsub
         call dd_get_cnodes(myid,isub)
      end do

      ! load arithmetic averages on edges
      glob_type = 2
      do isub = 1,nsub
         call dd_load_arithmetic_constraints(myid,isub,glob_type)
      end do

      ! prepare matrix C
      do isub = 1,nsub
         call dd_embed_cnodes(myid,isub,nndf_coarse,lnndf_coarse)
         call dd_prepare_c(myid,isub)
      end do

      ! prepare augmented matrix for BDDC
      do isub = 1,nsub
         call dd_prepare_aug(myid,comm_self,isub)
      end do

      ! prepare coarse space basis functions for BDDC
      do isub = 1,nsub
         call dd_prepare_coarse(myid,isub)
      end do

      call dd_finalize

!*******************************************AUX
! Measure time spent in DD module
      call MPI_BARRIER(comm_all,ierr)
      timeaux2 = MPI_WTIME()
      timeaux = timeaux2 - timeaux1
      if (myid.eq.0) then
         write(*,*) '***************************************'
         write(*,*) 'Time spent in DD setup is ',timeaux,' s'
         write(*,*) '***************************************'
      end if
!*******************************************AUX
   
      deallocate(nndf_coarse)

      ! MPI finalization
!***************************************************************PARALLEL
      call MPI_FINALIZE(ierr)
!***************************************************************PARALLEL

end program
