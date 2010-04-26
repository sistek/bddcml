module module_levels
!*******************
! Module for handling levels in multilevel BDDC 
! Jakub Sistek, Praha 2/2010

      implicit none

! type of real variables
      integer,parameter,private :: kr = kind(1.D0)
! numerical zero
      real(kr),parameter,private :: numerical_zero = 1.e-12_kr

! debugging 
      logical,parameter,private :: debug = .false.

! type for data about levels
      type levels_type
         integer ::             nelem    ! number of elements (subdomains) on level
         integer ::             nnod     ! number of nodes (corse nodes) on level
         integer ::             ndof     ! number of nodes (corse nodes) on level

         integer ::             nnodc    ! number of corners on level
         integer ::             nedge    ! number of edges on level
         integer ::             nface    ! number of faces on level

         ! description of subdomain mesh
         integer ::             linet    ! length of INET array 
         integer,allocatable ::  inet(:) ! INET array - indices of nodes on elements
         integer ::             lnnet    ! length of NNET array
         integer,allocatable ::  nnet(:) ! NNET array - number of nodes on elements
         integer ::             lnndf    ! length of NNDF array
         integer,allocatable ::  nndf(:) ! NNDF array - number of nodal degrees of freedom
         integer ::             liets     ! length of array IETS
         integer,allocatable ::  iets(:)  ! IETS array - indices of elements in numbering of elements in (level + 1)
         integer ::             lxyz1, lxyz2 ! length of array
         real(kr),allocatable :: xyz(:,:) ! coordinates of nodes of level
      end type levels_type

      integer, private ::                          nlevels
      integer, private ::                          llevels
      type(levels_type), allocatable, private ::    levels(:)

contains

!*************************
subroutine levels_init(nl)
!*************************
! Subroutine for initialization of levels data
      implicit none

! given number of levels
      integer,intent(in) :: nl

! initialize basic structure
      nlevels = nl
      llevels = nlevels
      allocate (levels(llevels))

end subroutine

!************************************************************************
subroutine levels_read_level_from_file(problemname,myid,comm,ndim,ilevel)
!************************************************************************
! Subroutine for loading first two levels from file

      use module_utils
      implicit none
      include "mpif.h"

! name of problem
      character(*),intent(in) :: problemname
! number of processor
      integer,intent(in) :: myid
! communicator
      integer,intent(in) :: comm
! dimension
      integer,intent(in) :: ndim
! index of level to import
      integer,intent(in) :: ilevel

! local variables
      integer :: idlevel
      integer :: ierr
      integer :: indlevel, nelem, nnod
      integer :: nnodc, nedge, nface

      integer ::             linet,   lnnet 
      integer, allocatable :: inet(:), nnet(:)

      integer ::              lxyz1, lxyz2
      real(kr), allocatable :: xyz(:,:)

      character(100) :: filename

      character(1) :: levelstring

      if (myid.eq.0) then
         if (ilevel.lt.10) then
            write(levelstring,'(i1)') ilevel
         else
            write(*,*) 'LEVELS_READ_LEVEL_FROM_FILE: Index of level too large...'
            call error_exit
         end if
         filename = trim(problemname)//'.L'//levelstring
         write(*,*) 'LEVELS_READ_LEVEL_FROM_FILE: Reading data from file ',trim(filename)
         call allocate_unit(idlevel)
         open (unit=idlevel,file=filename,status='old',form='formatted')
      end if

! read level header
      if (myid.eq.0) then
         read(idlevel,*) indlevel, nelem, nnod, linet
         read(idlevel,*) nnodc, nedge, nface
         if (indlevel.ne.ilevel) then
            write(*,*) 'LEVELS_READ_LEVEL_FROM_FILE: Data mismatch...'
            call error_exit
         end if
         if (nnod.ne.nnodc+nedge+nface) then
            write(*,*) 'LEVELS_READ_LEVEL_FROM_FILE: Number of pseudo nodes does not match.'
            call error_exit
         end if
      end if
!*****************************************************************MPI
      call MPI_BCAST(nelem,1, MPI_INTEGER, 0, comm, ierr)
      call MPI_BCAST(nnod, 1, MPI_INTEGER, 0, comm, ierr)
      call MPI_BCAST(linet,1, MPI_INTEGER, 0, comm, ierr)
      call MPI_BCAST(nnodc,1, MPI_INTEGER, 0, comm, ierr)
      call MPI_BCAST(nedge,1, MPI_INTEGER, 0, comm, ierr)
      call MPI_BCAST(nface,1, MPI_INTEGER, 0, comm, ierr)
!*****************************************************************MPI

! load data to structure
      levels(ilevel)%nelem = nelem
      levels(ilevel)%nnod  = nnod
      levels(ilevel)%linet = linet
      levels(ilevel)%nnodc = nnodc
      levels(ilevel)%nedge = nedge
      levels(ilevel)%nface = nface

      print *, 'I am here!'

! continue only for levels 2 and larger
      if (ilevel.ge.2) then

         lnnet = nelem
         lxyz1 = nnod
         lxyz2 = ndim

         allocate(inet(linet),nnet(lnnet),xyz(lxyz1,lxyz2))

         if (myid.eq.0) then
            read(idlevel,*) inet
            read(idlevel,*) nnet
            read(idlevel,*) xyz
         end if
!*****************************************************************MPI
         call MPI_BCAST(inet,linet,      MPI_INTEGER, 0, comm, ierr)
         call MPI_BCAST(nnet,lnnet,      MPI_INTEGER, 0, comm, ierr)
         call MPI_BCAST(xyz,lxyz1*lxyz2, MPI_INTEGER, 0, comm, ierr)
!*****************************************************************MPI

         ! save mesh into structure
         call levels_upload_level_mesh(ilevel,nelem,nnod,inet,linet,nnet,lnnet,xyz,lxyz1,lxyz2)

         deallocate(inet,nnet,xyz)
      end if
      if (myid.eq.0) then
         close(idlevel)
      end if

end subroutine

!***********************************************************************************************
subroutine levels_upload_level_mesh(ilevel, nelem,nnod, inet,linet, nnet,lnnet, xyz,lxyz1,lxyz2)
!***********************************************************************************************
! Subroutine for loading mesh of level into levels structure

      use module_utils
      implicit none

      integer,intent(in) :: ilevel, nelem, nnod
      integer,intent(in) :: linet
      integer,intent(in) ::  inet(linet)
      integer,intent(in) :: lnnet
      integer,intent(in) ::  nnet(lnnet)
      integer,intent(in) :: lxyz1, lxyz2
      real(kr),intent(in)::  xyz(lxyz1,lxyz2)

      ! local vars
      integer :: i, j

      ! check that array is allocated
      if (.not. allocated(levels) ) then
         write(*,*) 'LEVELS_UPLOAD_LEVEL_MESH: LEVELS not allocated for level  ',ilevel
         write(*,*) 'LEVELS_UPLOAD_LEVEL_MESH: Maybe missing call to LEVELS_INIT.'
         call error_exit
      end if
      ! check that enough space is available
      if (ilevel .gt. nlevels) then
         write(*,*) 'LEVELS_UPLOAD_LEVEL_MESH: ILEVEL exceeds size of LEVELS array.  ',ilevel
         call error_exit
      end if

      ! load data
      levels(ilevel)%nelem   = nelem
      levels(ilevel)%nnod    = nnod

      levels(ilevel)%linet   = linet
      allocate(levels(ilevel)%inet(linet))
      do i = 1,linet
         levels(ilevel)%inet(i) = inet(i)
      end do

      levels(ilevel)%lnnet   = lnnet
      allocate(levels(ilevel)%nnet(lnnet))
      do i = 1,lnnet
         levels(ilevel)%nnet(i) = nnet(i)
      end do

      levels(ilevel)%lxyz1 = lxyz1
      levels(ilevel)%lxyz2 = lxyz2
      allocate(levels(ilevel)%xyz(lxyz1,lxyz2))
      do j = 1,lxyz2
         do i = 1,lxyz1
            levels(ilevel)%xyz(i,j) = xyz(i,j)
         end do
      end do

end subroutine


!************************************
subroutine levels_clear_level(level)
!************************************
! Subroutine for deallocation data of one level
      implicit none

! local variables
      type(levels_type) :: level

! deallocate data
      if (allocated(level%inet)) then
         deallocate (level%inet)
      end if
      level%linet = 0
      if (allocated(level%nnet)) then
         deallocate (level%nnet)
      end if
      level%lnnet = 0
      if (allocated(level%nndf)) then
         deallocate (level%nndf)
      end if
      level%lnndf = 0
      if (allocated(level%iets)) then
         deallocate (level%iets)
      end if
      level%liets = 0
      if (allocated(level%xyz)) then
         deallocate (level%xyz)
      end if
      level%lxyz1 = 0
      level%lxyz2 = 0

      level%nelem = 0
      level%nnod  = 0
      level%ndof  = 0

      level%nnodc  = 0
      level%nedge  = 0
      level%nface  = 0

end subroutine

!*************************
subroutine levels_finalize
!*************************
! Subroutine for initialization of levels data
      implicit none

! local variables
      integer :: ilevel

! deallocate basic structure
      if (allocated(levels)) then
         ! clear data of particular subdomains
         do ilevel = 1,llevels
            call levels_clear_level(levels(ilevel))
         end do

         deallocate (levels)
      end if

end subroutine

end module module_levels

