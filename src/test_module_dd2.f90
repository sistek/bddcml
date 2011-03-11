program test_module_dd2
! Tester of module_dd
      use module_dd
      use module_pp
      use module_utils

      implicit none
      include "mpif.h"

      integer,parameter :: kr = kind(1.D0)

      ! parallel variables
      integer :: myid, comm_self, comm_all, nproc, ierr

      integer :: idpar, idcn, idglb

      integer :: matrixtype
      integer :: ndim, nsub, nelem, ndof, nnod, nnodc, linet

      integer :: isub, isub_loc, nsub_loc
      logical :: remove_original, keep_global

      integer ::             lnndfc
      integer,allocatable ::  nndfc(:)
      integer :: ncnodes, nedge, nface, nglb

      integer ::             lsub2proc
      integer,allocatable ::  sub2proc(:)
      integer ::             lindexsub
      integer,allocatable ::  indexsub(:)

      integer ::             lnnglb
      integer,allocatable ::  nnglb(:)
      integer ::             linglb
      integer,allocatable ::  inglb(:)

      integer :: glob_type

      character(90)  :: problemname 
      character(100) :: name

      real(kr) :: timeaux, timeaux1, timeaux2

      integer ::                         lsubdomains
      type(subdomain_type),allocatable :: subdomains(:)

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

!*******************************************AUX
! Measure time spent in DD module
      call MPI_BARRIER(comm_all,ierr)
      timeaux1 = MPI_WTIME()
!*******************************************AUX
      lsub2proc = nproc + 1
      allocate(sub2proc(lsub2proc))
      call pp_distribute_linearly(nsub,nproc,sub2proc,lsub2proc)
      nsub_loc = sub2proc(myid+2) - sub2proc(myid+1)
      lsubdomains = nsub_loc
      lindexsub   = nsub_loc
      allocate(subdomains(nsub_loc),indexsub(lindexsub))
      do isub_loc = 1,nsub_loc
         indexsub(isub_loc) = sub2proc(myid+1) + isub_loc - 1
      end do
      do isub_loc = 1,nsub_loc
         isub = indexsub(isub_loc)

         call dd_init(subdomains(isub_loc),isub,nsub,comm_all)
         call dd_read_mesh_from_file(subdomains(isub_loc),trim(problemname))
         ! SPD matrix
         matrixtype = 1
         call dd_read_matrix_from_file(subdomains(isub_loc),matrixtype,trim(problemname))
         call dd_assembly_local_matrix(subdomains(isub_loc))
         remove_original = .false.
         call dd_matrix_tri2blocktri(subdomains(isub_loc),remove_original)
         call dd_prepare_schur(subdomains(isub_loc),comm_self)
      end do

! create coarse mesh 
! read second level
      if (myid.eq.0) then
         name = trim(problemname)//'.CN'
         write(*,*) 'Reading data from file ',trim(name)
         call allocate_unit(idcn)
         open (unit=idcn,file=name,status='old',form='formatted')
         read(idcn,*) nnodc
         close (idcn)
         name = trim(problemname)//'.GLB'
         write(*,*) 'Reading data from file ',trim(name)
         call allocate_unit(idglb)
         open (unit=idglb,file=name,status='old',form='formatted')
         read(idcn,*) nglb,linglb
         lnnglb = nglb
         allocate(nnglb(lnnglb),inglb(linglb))
         read(idglb,*) inglb
         read(idglb,*) nnglb
         read(idglb,*) nedge, nface
         close(idglb)
         deallocate(nnglb,inglb)
      end if
!*****************************************************************MPI
      call MPI_BCAST(nnodc,1, MPI_INTEGER, 0, comm_all, ierr)
      call MPI_BCAST(nedge,1, MPI_INTEGER, 0, comm_all, ierr)
      call MPI_BCAST(nface,1, MPI_INTEGER, 0, comm_all, ierr)
!*****************************************************************MPI
       
      ncnodes = nnodc + nedge + nface

      ! prepare array of number of coarse dof in generalized coarse nodes
      lnndfc = ncnodes
      allocate(nndfc(lnndfc))

      ! auxiliary routine, until reading directly the globs
      do isub_loc = 1,nsub_loc
         call dd_construct_cnodes(subdomains(isub_loc))
      end do

      ! load arithmetic averages on edges
      glob_type = 2
      do isub_loc = 1,nsub_loc
         call dd_load_arithmetic_constraints(subdomains(isub_loc),glob_type)
      end do

      ! embed coarse nodes and create NNDFC array
      call dd_embed_cnodes(subdomains,lsubdomains, &
                           indexsub,lindexsub,& 
                           comm_all, nndfc,lnndfc)

      ! prepare matrix C
      do isub_loc = 1,nsub_loc
         call dd_prepare_c(subdomains(isub_loc))
      end do

      ! prepare augmented matrix for BDDC
      do isub_loc = 1,nsub_loc
         call dd_prepare_aug(subdomains(isub_loc),comm_self)
      end do

      ! prepare coarse space basis functions for BDDC
      do isub_loc = 1,nsub_loc
         keep_global = .false.
         call dd_prepare_coarse(subdomains(isub_loc),keep_global)
      end do

      do isub_loc = 1,nsub_loc
         call dd_finalize(subdomains(isub_loc))
      end do
      deallocate(subdomains)
      deallocate(sub2proc,indexsub)

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
   
      deallocate(nndfc)

      ! MPI finalization
!***************************************************************PARALLEL
      call MPI_FINALIZE(ierr)
!***************************************************************PARALLEL

end program
