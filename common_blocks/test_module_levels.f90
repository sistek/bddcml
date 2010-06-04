program test_module_levels
! Tester of module_levels
      use module_levels
      use module_utils
      implicit none
      include "mpif.h"

      integer,parameter :: kr = kind(1.D0)

      !  parallel variables
      integer :: myid, comm_all, comm_self, nproc, ierr
      integer :: idpar

      ! number of levels
      integer :: nlevels

      integer :: ndim, nsub, nelem, ndof, nnod, nnodc, linet

      integer :: matrixtype

      character(90)  :: problemname 
      character(100) :: name
      !character(100) :: filename


      ! MPI initialization
!***************************************************************PARALLEL
      call MPI_INIT(ierr)
      ! Communicator
      comm_all  = MPI_COMM_WORLD
      comm_self = MPI_COMM_SELF
      call MPI_COMM_RANK(comm_all,myid,ierr)
      call MPI_COMM_SIZE(comm_all,nproc,ierr)
!***************************************************************PARALLEL

      write (*,*) 'I am processor ',myid,': Hello nproc =',nproc


! Initial screen
      if (myid.eq.0) then
         write(*,'(a)') 'MULTILEVEL TESTER'
         write(*,'(a)') '================='

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


      write (*,*) 'myid = ',myid,': Initializing LEVELS.'
      !============
      ! number of levels
      nlevels = 2
      !============
      call levels_init(nlevels,nsub)

      ! create first two levels
      ! open file with description of pairs
!      ilevel = 1
!      call levels_read_level_from_file(problemname,myid,comm_all,ndim,ilevel)
!
!      ilevel = 2
!      call levels_read_level_from_file(problemname,myid,comm_all,ndim,ilevel)
!
!      ! associate subdomains with first level
!      ilevel = 1
!      call levels_prepare_standard_level(ilevel,nsub,1,nsub)
!
!      call levels_prepare_last_level(myid,nproc,comm_all,comm_self,matrixtype,ndim,problemname)
      matrixtype = 1 ! SPD matrix
      call levels_pc_setup(problemname,myid,nproc,comm_all,comm_self,matrixtype,ndim,nsub)

   !   lvec = 15
   !   allocate(vec(lvec))
   !   vec = 1.0_kr
   !   call levels_pc_apply(vec,lvec, )
   !   if (myid.eq.0) then
   !      write(*,*) 'vec:'
   !      write(*,'(e18.9)') vec
   !   end if
   !   deallocate(vec)

      write (*,*) 'myid = ',myid,': Finalize LEVELS.'
      call levels_finalize

      ! MPI finalization
!***************************************************************PARALLEL
      call MPI_FINALIZE(ierr)
!***************************************************************PARALLEL

end program
