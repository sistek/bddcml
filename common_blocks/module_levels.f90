module module_levels
!*******************
! Module for handling levels in multilevel BDDC 
! Jakub Sistek, Praha 2/2010

!     definition of MUMPS structure
      use dmumps_struc_def
      use module_dd

      implicit none

! adjustable parameters ############################
! type of real variables
      integer,parameter,private :: kr = kind(1.D0)
! numerical zero
      real(kr),parameter,private :: numerical_zero = 1.e-12_kr
! debugging 
      logical,parameter,private :: debug = .true.
! profiling 
      logical,parameter,private :: profile = .true.
! adjustable parameters ############################

! type for data about levels
      type levels_type

         integer :: comm_all  ! global communicator for the level
         integer :: comm_self ! local communicator for the level

         integer :: nelem    ! number of elements on level
         integer :: nnod     ! number of nodes on level
         integer :: ndof     ! number of dof on level

         ! description of mesh
         integer          ::   linet    ! length of INET array 
         integer,pointer  ::    inet(:) ! INET array - indices of nodes on elements
         integer          ::   lnnet    ! length of NNET array
         integer,pointer  ::    nnet(:) ! NNET array - number of nodes on elements
         integer          ::   lnndf    ! length of NNDF array
         integer,pointer  ::    nndf(:) ! NNDF array - number of nodal degrees of freedom
         integer          ::  lxyz1, lxyz2 ! length of array
         real(kr),pointer ::   xyz(:,:) ! coordinates of nodes of level
         integer          ::  lifix     ! length of IFIX array
         integer,pointer  ::   ifix(:)  ! indices of fixed variables
         integer          ::  lfixv     ! length of FIXV array
         real(kr),pointer ::   fixv(:)  ! values of fixed variables
         integer          ::  lrhs      ! length of RHS array
         real(kr),pointer ::   rhs(:)   ! values of right hand side
         integer          ::  lsol      ! length of SOL array
         real(kr),pointer ::   sol(:)   ! values of solution (initial at start, final at end)

         ! subdomain data
         integer ::              nsub     ! number of subdomains on level
         integer ::              nsub_loc ! number of subdomains on level assigned to the processor

         integer ::                          lsubdomains    ! lenth of array of subdomains (= nsub_loc)
         type(subdomain_type), allocatable :: subdomains(:) ! array of subdomains on level

         integer             :: lindexsub    
         integer,allocatable ::  indexsub(:) ! indices of elements in array subdomains
         integer             :: lsub2proc    ! nproc + 1 
         integer,allocatable ::  sub2proc(:) ! division of subdomains for processors 
                                             ! sub2proc(1) = 1, sub2proc(iproc+1) = sub2proc(iproc) + nsub

         integer ::             liets     ! length of array IETS
         integer,allocatable ::  iets(:)  ! IETS array - indices of elements in numbering of elements in (level + 1)

         ! descrtiption of coarse data
         integer ::             ncorner  ! number of corners on level
         integer ::             nedge    ! number of edges on level
         integer ::             nface    ! number of faces on level
         integer ::             nnodc    ! number of pseudo nodes on level = ncorner + nedge + nface
         integer ::             ndofc    ! number of coarse DOFs

         integer::             linetc    ! length of INET array 
         integer,allocatable :: inetc(:) ! INET array - indices of nodes on elements
         integer::             lnnetc    ! length of NNET array
         integer,allocatable :: nnetc(:) ! NNET array - number of nodes on elements
         integer::             lnndfc    ! length of NNDF array
         integer,allocatable :: nndfc(:) ! NNDF array - number of nodal degrees of freedom
         integer ::             lxyzc1, lxyzc2 ! length of array
         real(kr),allocatable :: xyzc(:,:) ! coordinates of nodes of level
         integer ::             lifixc   ! length of IFIXC array
         integer,allocatable ::  ifixc(:) ! indices of fixed variables
         integer ::             lfixvc   ! length of FIXVC array
         real(kr),allocatable :: fixvc(:) ! values of fixed variables
         integer ::             lrhsc    ! length of RHSC array
         real(kr),allocatable :: rhsc(:) ! values of right hand side
         integer ::             lsolc    ! length of SOLC array
         real(kr),allocatable :: solc(:) ! coarse residual/solution at level
         logical :: use_initial_solution = .false. ! should some initial solution be used for iterations?

         logical :: is_level_prepared = .false. ! level prepared

      end type levels_type

      integer, private ::                                  nlevels
      integer, private ::                                  iactive_level
      integer, private ::                                  llevels
      type(levels_type), allocatable, target, private ::    levels(:)

      type(DMUMPS_STRUC), private :: mumps_coarse  
      logical :: is_mumps_coarse_ready = .false.

contains

!************************************************************************
subroutine levels_init(nl,nsublev,lnsublev,nelem,nnod,ndof,&
                       inet,linet,nnet,lnnet,nndf,lnndf,xyz,lxyz1,lxyz2,&
                       ifix,lifix, fixv,lfixv, rhs,lrhs, sol,lsol)
!************************************************************************
! Subroutine for initialization of levels data
      use module_utils
      implicit none
      include "mpif.h"

! given number of levels
      integer,intent(in) :: nl
! number of subdomains in all levels
      integer,intent(in) :: lnsublev
      integer,intent(in) ::  nsublev(lnsublev)

      integer, intent(in):: nelem, nnod, ndof
      integer, intent(in):: linet
      integer, intent(in)::  inet(linet)
      integer, intent(in):: lnnet
      integer, intent(in)::  nnet(lnnet)
      integer, intent(in):: lnndf
      integer, intent(in)::  nndf(lnndf)
      integer, intent(in):: lxyz1, lxyz2
      real(kr), intent(in):: xyz(lxyz1,lxyz2)
      integer, intent(in):: lifix
      integer, intent(in)::  ifix(lifix)
      integer, intent(in):: lfixv
      real(kr), intent(in):: fixv(lfixv)
      integer, intent(in):: lrhs
      real(kr), intent(in):: rhs(lrhs)
      integer, intent(in):: lsol
      real(kr), intent(in):: sol(lsol)

      ! local vars 
      integer,parameter :: ilevel = 0
      integer :: i, j

! initial checks 
      if (nl.lt.2) then
         call error('LEVELS_INIT','Number of levels must be at least 2.')
      end if
      if (nsublev(nl).ne.1) then
         call error('LEVELS_INIT','Number of subdomains at last level must be 1.')
      end if

! initialize basic structure
      nlevels = nl
      llevels = nlevels
      allocate(levels(0:llevels))

! initialize zero level
      levels(ilevel)%nsub  = nelem
      levels(ilevel)%nnodc = nnod
      levels(ilevel)%ndofc = ndof

      levels(ilevel)%linetc = linet  
      allocate(levels(ilevel)%inetc(levels(ilevel)%linetc))
      do i = 1,levels(ilevel)%linetc
         levels(ilevel)%inetc(i) = inet(i)
      end do
      levels(ilevel)%lnnetc = lnnet  
      allocate(levels(ilevel)%nnetc(levels(ilevel)%lnnetc))
      do i = 1,levels(ilevel)%lnnetc
         levels(ilevel)%nnetc(i) = nnet(i)
      end do
      levels(ilevel)%lnndfc = lnndf  
      allocate(levels(ilevel)%nndfc(levels(ilevel)%lnndfc))
      do i = 1,levels(ilevel)%lnndfc
         levels(ilevel)%nndfc(i) = nndf(i)
      end do
      levels(ilevel)%lxyzc1 = lxyz1  
      levels(ilevel)%lxyzc2 = lxyz2  
      allocate(levels(ilevel)%xyzc(levels(ilevel)%lxyzc1,levels(ilevel)%lxyzc2))
      do j = 1,levels(ilevel)%lxyzc2
         do i = 1,levels(ilevel)%lxyzc1
            levels(ilevel)%xyzc(i,j) = xyz(i,j)
         end do
      end do
      levels(ilevel)%lifixc = lifix
      allocate(levels(ilevel)%ifixc(levels(ilevel)%lifixc))
      do i = 1,levels(ilevel)%lifixc
         levels(ilevel)%ifixc(i) = ifix(i)
      end do
      levels(ilevel)%lfixvc = lfixv
      allocate(levels(ilevel)%fixvc(levels(ilevel)%lfixvc))
      do i = 1,levels(ilevel)%lfixvc
         levels(ilevel)%fixvc(i) = fixv(i)
      end do
      levels(ilevel)%lrhsc = lrhs
      allocate(levels(ilevel)%rhsc(levels(ilevel)%lrhsc))
      do i = 1,levels(ilevel)%lrhsc
         levels(ilevel)%rhsc(i) = rhs(i)
      end do
      levels(ilevel)%lsolc = lsol
      allocate(levels(ilevel)%solc(levels(ilevel)%lsolc))
      do i = 1,levels(ilevel)%lsolc
         levels(ilevel)%solc(i) = sol(i)
      end do
      if (any(sol.ne.0._kr)) then
         levels(ilevel)%use_initial_solution = .true.
      else
         levels(ilevel)%use_initial_solution = .false.
      end if

      levels(ilevel)%is_level_prepared = .true.

end subroutine

!*************************************************************************************
subroutine levels_pc_setup(problemname,nsublev,lnsublev,load_division,load_globs,correct_division,&
                           neighbouring,matrixtype,ndim, meshdim, use_arithmetic, use_adaptive)
!*************************************************************************************
! subroutine for multilevel BDDC preconditioner setup
      use module_pp
      implicit none
      include "mpif.h"

! name of problem
      character(*),intent(in) :: problemname
! number of subdomains in all levels
      integer,intent(in) :: lnsublev
      integer,intent(in) ::  nsublev(lnsublev)
! use prepared division into subdomains?
      logical,intent(in) :: load_division
! use prepared selected corners and globs?
      logical,intent(in) :: load_globs
! should disconnected subdomains be connected?
      logical,intent(in) :: correct_division
! how many nodes are shared to called elements adjacent
      integer,intent(in) :: neighbouring
! type of matrix (0 - nosymetric, 1 - SPD, 2 - general symmetric)
      integer,intent(in) :: matrixtype
! dimension
      integer,intent(in) :: ndim
! dimension of mesh
      integer,intent(in) :: meshdim
! Use arithmetic averages on globs as constraints?
      logical,intent(in) :: use_arithmetic
! Use adaptive constraints on faces?
      logical,intent(in) :: use_adaptive

      ! local vars 
      integer :: isub, isub_loc, nsub, nsub_loc
      integer :: nproc, comm_self, comm_all, ierr, myid

      ! no debug
      !iactive_level = 1
      !call levels_read_level_from_file(problemname,comm_all,ndim,iactive_level)

      do iactive_level = 1,nlevels

         ! TODO: create proper communicators
         comm_all  = MPI_COMM_WORLD
         comm_self = MPI_COMM_SELF

         levels(iactive_level)%comm_all  = comm_all
         levels(iactive_level)%comm_self = comm_self

         ! orient in the communicator
         call MPI_COMM_RANK(comm_all,myid,ierr)
         call MPI_COMM_SIZE(comm_all,nproc,ierr)

         nsub = nsublev(iactive_level)
         levels(iactive_level)%nsub = nsub

         ! prepare distribution of subdomains for the level
         levels(iactive_level)%lsub2proc = nproc + 1
         allocate(levels(iactive_level)%sub2proc(levels(iactive_level)%lsub2proc))
         call pp_distribute_subdomains(nsub,nproc,levels(iactive_level)%sub2proc,levels(iactive_level)%lsub2proc)

         nsub_loc    = levels(iactive_level)%sub2proc(myid+2) - levels(iactive_level)%sub2proc(myid+1)
         levels(iactive_level)%nsub_loc = nsub_loc

         ! prepare array of indices of local subdomains for each level
         levels(iactive_level)%lindexsub = nsub_loc
         allocate(levels(iactive_level)%indexsub(levels(iactive_level)%lindexsub))
         do isub_loc = 1,nsub_loc
            isub = levels(iactive_level)%sub2proc(myid+1) + isub_loc - 1

            levels(iactive_level)%indexsub(isub_loc) = isub
         end do

         ! prepare space for subdomains
         levels(iactive_level)%lsubdomains = nsub_loc
         allocate(levels(iactive_level)%subdomains(levels(iactive_level)%lsubdomains))
         do isub_loc = 1,nsub_loc
            isub = levels(iactive_level)%indexsub(isub_loc)
            call dd_init(levels(iactive_level)%subdomains(isub_loc),isub,nsub,comm_all)
         end do
      end do

      ! associate subdomains with first level
      do iactive_level = 1,nlevels-1
         comm_all = levels(iactive_level)%comm_all
         ! orient in the communicator
         call MPI_COMM_RANK(comm_all,myid,ierr)
         call MPI_COMM_SIZE(comm_all,nproc,ierr)
         if (debug .and. myid.eq.0) then
            write(*,*) 'Preparing level',iactive_level
         end if
         call levels_prepare_standard_level(problemname,load_division,load_globs,&
                                            correct_division,neighbouring,matrixtype,ndim,meshdim,iactive_level,&
                                            use_arithmetic,use_adaptive)
      end do

      iactive_level = nlevels
      if (debug .and. myid.eq.0) then
         write(*,*) 'Preparing level',iactive_level
      end if
      comm_all = levels(iactive_level)%comm_all
      ! orient in the communicator
      call MPI_COMM_RANK(comm_all,myid,ierr)

      call levels_prepare_last_level(matrixtype)

end subroutine

!*******************************************************************
subroutine levels_read_level_from_file(problemname,comm,ndim,ilevel)
!*******************************************************************
! Subroutine for loading first two levels from file

      use module_utils
      implicit none
      include "mpif.h"

! name of problem
      character(*),intent(in) :: problemname
! communicator
      integer,intent(in) :: comm
! dimension
      integer,intent(in) :: ndim
! index of level to import
      integer,intent(in) :: ilevel

! local variables
! number of processor
      integer :: myid
      integer :: idlevel
      integer :: ierr
      integer :: indlevel, nelem, nnod, linet
      integer :: ncorner, nedge, nface, nnodc, ndofc
      integer :: nsub

      integer ::             linetc,   lnnetc ,  lnndfc
      integer, allocatable :: inetc(:), nnetc(:), nndfc(:)

      integer ::              lxyzc1, lxyzc2
      real(kr), allocatable :: xyzc(:,:)

      integer :: i,j

      character(500) :: filename

      character(1) :: levelstring

      ! find my id in the communicator
      call MPI_COMM_RANK(comm,myid,ierr)

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
         read(idlevel,*) ncorner, nedge, nface, linetc
         if (indlevel.ne.ilevel) then
            write(*,*) 'LEVELS_READ_LEVEL_FROM_FILE: Data mismatch...'
            call error_exit
         end if
      end if
!*****************************************************************MPI
      call MPI_BCAST(nelem, 1, MPI_INTEGER, 0, comm, ierr)
      call MPI_BCAST(nnod,  1, MPI_INTEGER, 0, comm, ierr)
      call MPI_BCAST(linet, 1, MPI_INTEGER, 0, comm, ierr)
      call MPI_BCAST(ncorner, 1, MPI_INTEGER, 0, comm, ierr)
      call MPI_BCAST(nedge, 1, MPI_INTEGER, 0, comm, ierr)
      call MPI_BCAST(nface, 1, MPI_INTEGER, 0, comm, ierr)
      call MPI_BCAST(linetc,1, MPI_INTEGER, 0, comm, ierr)
!*****************************************************************MPI

      nnodc = ncorner + nedge + nface

      nsub = levels(ilevel)%nsub
      lnnetc = nsub
      lnndfc = nnodc
      lxyzc1 = nnodc
      lxyzc2 = ndim

      allocate(inetc(linetc),nnetc(lnnetc),nndfc(lnndfc),xyzc(lxyzc1,lxyzc2))

      if (myid.eq.0) then
         read(idlevel,*) inetc
         read(idlevel,*) nnetc
         do i = 1,nnodc
            read(idlevel,*) (xyzc(i,j),j = 1,ndim)
         end do
      end if
!*****************************************************************MPI
      call MPI_BCAST(inetc,linetc,       MPI_INTEGER,          0, comm, ierr)
      call MPI_BCAST(nnetc,lnnetc,       MPI_INTEGER,          0, comm, ierr)
      call MPI_BCAST(xyzc,lxyzc1*lxyzc2, MPI_DOUBLE_PRECISION, 0, comm, ierr)
!*****************************************************************MPI
      ! when reading from L* file, array of numbers of degrees of freedom at nodes is not yet available
      ndofc = 0
      call zero(nndfc,lnndfc)

      ! save mesh into structure
      levels(ilevel)%nelem = nelem
      levels(ilevel)%nnod  = nnod
      levels(ilevel)%linet = linet
      call levels_upload_level_mesh(ilevel,&
                                    ncorner,nedge,nface,nnodc,ndofc,&
                                    inetc,linetc,nnetc,lnnetc,nndfc,lnndfc, xyzc,lxyzc1,lxyzc2)
      deallocate(inetc,nnetc,nndfc,xyzc)

      if (myid.eq.0) then
         close(idlevel)
      end if

end subroutine

!***********************************************************************************************
subroutine levels_upload_level_mesh(ilevel,ncorner,nedge,nface,nnodc,ndofc,&
                                    inetc,linetc,nnetc,lnnetc,nndfc,lnndfc,xyzc,lxyzc1,lxyzc2)
!***********************************************************************************************
! Subroutine for loading mesh of level into levels structure

      use module_utils
      implicit none

      integer,intent(in) :: ilevel
      integer,intent(in) :: ncorner, nedge, nface, nnodc, ndofc
      integer,intent(in) :: linetc
      integer,intent(in) ::  inetc(linetc)
      integer,intent(in) :: lnnetc
      integer,intent(in) ::  nnetc(lnnetc)
      integer,intent(in) :: lnndfc
      integer,intent(in) ::  nndfc(lnndfc)
      integer,intent(in) :: lxyzc1, lxyzc2
      real(kr),intent(in)::  xyzc(lxyzc1,lxyzc2)

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
      ! check sum
      if (nnodc .ne. ncorner + nedge + nface) then
         write(*,*) 'LEVELS_UPLOAD_LEVEL_MESH: coarse nodes number mismatch. nnodc = ',nnodc,&
                    'sum of globs',ncorner + nedge + nface
         call error_exit
      end if

! load data to structure
      levels(ilevel)%ncorner   = ncorner
      levels(ilevel)%nedge     = nedge
      levels(ilevel)%nface     = nface
      levels(ilevel)%nnodc     = nnodc
      levels(ilevel)%linetc    = linetc
      levels(ilevel)%ndofc     = ndofc

      allocate(levels(ilevel)%inetc(linetc))
      do i = 1,linetc
         levels(ilevel)%inetc(i) = inetc(i)
      end do

      levels(ilevel)%lnnetc  = lnnetc
      allocate(levels(ilevel)%nnetc(lnnetc))
      do i = 1,lnnetc
         levels(ilevel)%nnetc(i) = nnetc(i)
      end do

      levels(ilevel)%lnndfc = lnndfc
      allocate(levels(ilevel)%nndfc(lnndfc))
      do i = 1,lnndfc
         levels(ilevel)%nndfc(i) = nndfc(i)
      end do

      levels(ilevel)%lxyzc1 = lxyzc1
      levels(ilevel)%lxyzc2 = lxyzc2
      allocate(levels(ilevel)%xyzc(lxyzc1,lxyzc2))
      do j = 1,lxyzc2
         do i = 1,lxyzc1
            levels(ilevel)%xyzc(i,j) = xyzc(i,j)
         end do
      end do

end subroutine

!*****************************************************************************************
subroutine levels_prepare_standard_level(problemname,load_division,load_globs,&
                                         correct_division,neighbouring,matrixtype,ndim,meshdim,ilevel,&
                                         use_arithmetic,use_adaptive)
!*****************************************************************************************
! Subroutine for building the standard level
      use module_pp
      use module_adaptivity
      use module_sm
      use module_utils
      implicit none
      include "mpif.h"

      character(*),intent(in) :: problemname
      logical,intent(in) :: load_division
      logical,intent(in) :: load_globs
      logical,intent(in) :: correct_division
      integer,intent(in) :: neighbouring
      integer,intent(in) :: matrixtype
      integer,intent(in) :: ndim
      integer,intent(in) :: meshdim
      integer,intent(in) :: ilevel     ! index of level
      ! Use arithmetic averages on globs as constraints?
      logical,intent(in) :: use_arithmetic
      ! Use adaptive constraints on faces?
      logical,intent(in) :: use_adaptive

      ! local vars
      integer :: myid
      integer :: nproc
      integer :: comm_all, comm_self, ierr
      integer :: graphtype 
      integer :: ides, idcn, idglb

      integer :: ncorner, nedge, nface, isub, nnodc, ndofc, nelem, nnod, nnodi
      integer :: edgecut, ncornermin
      integer :: isub_loc, nsub, glbtype, nglb, nglbn

      integer ::             lifixs
      integer,allocatable ::  ifixs(:)
      integer :: ndofs, nelems, nnods, neighball
      integer :: pkadjsub

      integer ::           lkglobs
      integer,allocatable:: kglobs(:)
      integer ::           ltypeglobs
      integer,allocatable:: typeglobs(:)
      integer :: indglob, indinodc

      integer ::           lkadjsub
      integer,allocatable:: kadjsub(:)

      integer ::            linodc
      integer,allocatable :: inodc(:)
      integer ::            lnnglb
      integer,allocatable :: nnglb(:)
      integer ::            linglb
      integer,allocatable :: inglb(:)
      integer ::            lkglb
      integer,allocatable :: kglb(:)
      
      integer ::            lrhss,   lfixvs
      real(kr),allocatable:: rhss(:), fixvs(:)

      integer ::            linetc
      integer,allocatable :: inetc(:)
      integer ::            lnnetc
      integer,allocatable :: nnetc(:)
      integer ::            lnndfc
      integer,allocatable :: nndfc(:)
      integer ::            lxyzc1, lxyzc2
      real(kr),allocatable:: xyzc(:,:)
      integer ::            lkinetc
      integer,allocatable :: kinetc(:)

      integer,allocatable :: nnetcaux(:)
      integer,allocatable :: inetcaux(:)

      integer ::            lcnodenumbers
      integer,allocatable :: cnodenumbers(:)

      integer :: iglb, inc, indnc, inod, pointinglb
      integer :: indxyzc, ndofcs, nnodcs, pointinetc


      logical :: remove_original 
      logical :: remove_bc_nodes 
      logical :: keep_global 

      ! adaptivity variables
      integer :: idpair
      character(200) :: filename
      integer :: npair, npair_locx

      ! time variables
      real(kr) :: t_division, t_globs, t_matrix_import, t_adjacency, t_loc_mesh,&
                  t_loc_interface, t_loc_globs, t_loc_bc, t_loc_adjacency

      if (.not.levels(ilevel-1)%is_level_prepared) then
         call error('LEVELS_PREPARE_STANDARD_LEVEL', 'Previous level not ready.')
      end if

      ! orient in the communicator
      comm_all  = levels(ilevel)%comm_all
      comm_self = levels(ilevel)%comm_self
      call MPI_COMM_RANK(comm_all,myid,ierr)
      call MPI_COMM_SIZE(comm_all,nproc,ierr)

      ! make the connection with previous level
      levels(ilevel)%nelem = levels(ilevel-1)%nsub
      levels(ilevel)%nnod  = levels(ilevel-1)%nnodc
      levels(ilevel)%ndof  = levels(ilevel-1)%ndofc

      levels(ilevel)%linet = levels(ilevel-1)%linetc  
      levels(ilevel)%inet  => levels(ilevel-1)%inetc  
      levels(ilevel)%lnnet = levels(ilevel-1)%lnnetc  
      levels(ilevel)%nnet  => levels(ilevel-1)%nnetc  
      levels(ilevel)%lnndf = levels(ilevel-1)%lnndfc  
      levels(ilevel)%nndf  => levels(ilevel-1)%nndfc  
      levels(ilevel)%lxyz1 = levels(ilevel-1)%lxyzc1  
      levels(ilevel)%lxyz2 = levels(ilevel-1)%lxyzc2  
      levels(ilevel)%xyz   => levels(ilevel-1)%xyzc  
      levels(ilevel)%lifix = levels(ilevel-1)%lifixc  
      levels(ilevel)%ifix  => levels(ilevel-1)%ifixc  
      levels(ilevel)%lfixv = levels(ilevel-1)%lfixvc  
      levels(ilevel)%fixv  => levels(ilevel-1)%fixvc  
      if (ilevel.eq.1) then
         levels(ilevel)%lrhs = levels(ilevel-1)%lrhsc  
         levels(ilevel)%rhs  => levels(ilevel-1)%rhsc  
      end if
      levels(ilevel)%lsol = levels(ilevel-1)%lsolc  
      levels(ilevel)%sol  => levels(ilevel-1)%solc  

      ! initialize values
      nsub  = levels(ilevel)%nsub
      nnod  = levels(ilevel)%nnod
      nelem = levels(ilevel)%nelem

      ! make division into subdomains
      graphtype = 0 ! not weighted
      levels(ilevel)%liets = levels(ilevel)%nelem
      allocate(levels(ilevel)%iets(levels(ilevel)%liets))
      if (ilevel.eq.1 .and. load_division) then
         ! read the division from file *.ES
         if (myid.eq.0) then
            filename = trim(problemname)//'.ES'
            call allocate_unit(ides)
            open (unit=ides,file=filename,status='old',form='formatted')
            rewind ides

            read(ides,*) levels(ilevel)%iets
            close (ides)
         end if
!***************************************************************PARALLEL
         call MPI_BCAST(levels(ilevel)%iets,levels(ilevel)%liets, MPI_INTEGER, 0, comm_all, ierr)
!***************************************************************PARALLEL
         if (myid.eq.0) then
            write(*,*) ' Mesh division loaded from file ',trim(filename)
            call flush(6)
         end if
      else
         ! create division
         ! TODO: make this parallel 
!-----profile
         if (profile) then
            call MPI_BARRIER(comm_all,ierr)
            call time_start
         end if
!-----profile
         if (myid.eq.0) then
            call pp_divide_mesh(graphtype,correct_division,neighbouring,&
                                levels(ilevel)%nelem,levels(ilevel)%nnod,&
                                levels(ilevel)%inet,levels(ilevel)%linet,levels(ilevel)%nnet,levels(ilevel)%lnnet,nsub,&
                                edgecut,levels(ilevel)%iets,levels(ilevel)%liets)
         end if 
!***************************************************************PARALLEL
         call MPI_BCAST(levels(ilevel)%iets,levels(ilevel)%liets, MPI_INTEGER, 0, comm_all, ierr)
!***************************************************************PARALLEL
         if (myid.eq.0) then
            write(*,*) ' Mesh division created.'
            call flush(6)
         end if 
!-----profile
         if (profile) then
            call MPI_BARRIER(comm_all,ierr)
            call time_end(t_division)
            if (myid.eq.0) then
               write(*,*) '***********************************PROFILING'
               write(*,*) 'Time of creating division into subdomains: ',t_division,' s.'
               write(*,*) '***********************************PROFILING'
               call flush(6)
            end if
         end if
!-----profile
      end if

!-----profile 
      if (profile) then
         call MPI_BARRIER(comm_all,ierr)
         call time_start
      end if
!-----profile

      ! create subdomain mesh
      do isub_loc = 1,levels(ilevel)%nsub_loc
         isub = levels(ilevel)%indexsub(isub_loc)
         call dd_localize_mesh(levels(ilevel)%subdomains(isub_loc),isub,ndim,levels(ilevel)%nelem,levels(ilevel)%nnod,&
                               levels(ilevel)%inet,levels(ilevel)%linet,&
                               levels(ilevel)%nnet,levels(ilevel)%lnnet,&
                               levels(ilevel)%nndf,levels(ilevel)%lnndf,&
                               levels(ilevel)%xyz,levels(ilevel)%lxyz1,levels(ilevel)%lxyz2,&
                               levels(ilevel)%iets,levels(ilevel)%liets)
      end do

!-----profile
      if (profile) then
         call MPI_BARRIER(comm_all,ierr)
         call time_end(t_loc_mesh)
         if (myid.eq.0) then
            write(*,*) '***********************************PROFILING'
            write(*,*) 'Time of localizing mesh : ',t_loc_mesh,' s.'
            write(*,*) '***********************************PROFILING'
            call flush(6)
         end if
      end if
!-----profile

      ! for communication, any shared nodes are considered
!-----profile 
      if (profile) then
         call MPI_BARRIER(comm_all,ierr)
         call time_start
      end if
!-----profile
      ! TODO: make this parallel
      neighball = 1
      lkadjsub = nsub*nsub
      allocate(kadjsub(lkadjsub))
      if (myid.eq.0) then
         call pp_get_sub_neighbours(neighball,nelem,nnod,nsub,&
                                    levels(ilevel)%inet,levels(ilevel)%linet,&
                                    levels(ilevel)%nnet,levels(ilevel)%lnnet,&
                                    levels(ilevel)%iets,levels(ilevel)%liets,&
                                    kadjsub,lkadjsub)
      end if
!***************************************************************PARALLEL
      call MPI_BCAST(kadjsub,lkadjsub, MPI_INTEGER, 0, comm_all, ierr)
!***************************************************************PARALLEL
      ! debug
      !write(*,*) 'kadjsub'
      !do isub = 1,nsub
      !   pkadjsub = (isub-1)*nsub
      !   write(*,*) (kadjsub(pkadjsub+j),j = 1,nsub)
      !end do
!-----profile
      if (profile) then
         call MPI_BARRIER(comm_all,ierr)
         call time_end(t_adjacency)
         if (myid.eq.0) then
            write(*,*) '***********************************PROFILING'
            write(*,*) 'Time of finding neighbours : ',t_adjacency,' s.'
            write(*,*) '***********************************PROFILING'
            call flush(6)
         end if
      end if
!-----profile
!-----profile 
      if (profile) then
         call MPI_BARRIER(comm_all,ierr)
         call time_start
      end if
!-----profile

      ! create subdomain neighbours
      do isub_loc = 1,levels(ilevel)%nsub_loc
         isub = levels(ilevel)%indexsub(isub_loc)
         pkadjsub = (isub-1)*nsub
         call dd_localize_adj(levels(ilevel)%subdomains(isub_loc),nsub,kadjsub(pkadjsub+1),nsub)
      end do
      deallocate(kadjsub)

!-----profile
      if (profile) then
         call MPI_BARRIER(comm_all,ierr)
         call time_end(t_loc_adjacency)
         if (myid.eq.0) then
            write(*,*) '***********************************PROFILING'
            write(*,*) 'Time of localizing neighbours : ',t_loc_adjacency,' s.'
            write(*,*) '***********************************PROFILING'
            call flush(6)
         end if
      end if
!-----profile

!-----profile 
      if (profile) then
         call MPI_BARRIER(comm_all,ierr)
         call time_start
      end if
!-----profile

      ! generate interface data and shared nodes on level
      call dd_create_neighbouring(levels(ilevel)%subdomains,levels(ilevel)%lsubdomains,&
                                  levels(ilevel)%sub2proc,levels(ilevel)%lsub2proc,&
                                  levels(ilevel)%indexsub,levels(ilevel)%lindexsub, &
                                  comm_all)

!-----profile
      if (profile) then
         call MPI_BARRIER(comm_all,ierr)
         call time_end(t_loc_interface)
         if (myid.eq.0) then
            write(*,*) '***********************************PROFILING'
            write(*,*) 'Time of generating neighbouring data: ',t_loc_interface,' s.'
            write(*,*) '***********************************PROFILING'
            call flush(6)
         end if
      end if
!-----profile

      if (ilevel.eq.1 .and. load_globs) then
         ! read list of corners from *.CN file
         ! read list of globs from *.GLB file

         if (myid.eq.0) then
            filename = trim(problemname)//'.CN'
            call allocate_unit(idcn)
            open (unit=idcn,file=filename,status='old',form='formatted')
            rewind idcn

            read(idcn,*) ncorner
         end if
!***************************************************************PARALLEL
         call MPI_BCAST(ncorner,1, MPI_INTEGER, 0, comm_all, ierr)
!***************************************************************PARALLEL
         linodc = ncorner
         allocate(inodc(linodc))
         if (myid.eq.0) then
            read(idcn,*) inodc
            close(idcn)
         end if
!***************************************************************PARALLEL
         call MPI_BCAST(inodc,linodc, MPI_INTEGER, 0, comm_all, ierr)
!***************************************************************PARALLEL
         if (myid.eq.0) then
            write(*,*) ' Corners loaded from file ',trim(filename)
            call flush(6)
         end if

         if (myid.eq.0) then
            filename = trim(problemname)//'.GLB'
            call allocate_unit(idglb)
            open (unit=idglb,file=filename,status='old',form='formatted')
            rewind idglb

            read(idglb,*) nglb, linglb
         end if
!***************************************************************PARALLEL
         call MPI_BCAST(nglb,  1, MPI_INTEGER, 0, comm_all, ierr)
         call MPI_BCAST(linglb,1, MPI_INTEGER, 0, comm_all, ierr)
!***************************************************************PARALLEL
         lnnglb = nglb
         allocate (inglb(linglb),nnglb(lnnglb))
         if (myid.eq.0) then
            read(idglb,*) inglb
            read(idglb,*) nnglb
            read(idglb,*) nedge, nface
            close(idglb)
         end if
!***************************************************************PARALLEL
         call MPI_BCAST(inglb,linglb, MPI_INTEGER, 0, comm_all, ierr)
         call MPI_BCAST(nnglb,lnnglb, MPI_INTEGER, 0, comm_all, ierr)
         call MPI_BCAST(nedge, 1,     MPI_INTEGER, 0, comm_all, ierr)
         call MPI_BCAST(nface, 1,     MPI_INTEGER, 0, comm_all, ierr)
!***************************************************************PARALLEL
         nnodc = ncorner + nedge + nface
         if (myid.eq.0) then
            write(*,*) ' Globs loaded from file ',trim(filename)
            call flush(6)
         end if
      else
         ! TODO: make this parallel
!-----profile 
         if (profile) then
            call MPI_BARRIER(comm_all,ierr)
            call time_start
         end if
!-----profile
         ! find corners
         lkglobs    = nnod
         ltypeglobs = nnod
         ! on the first level, remove nodes at boundary consitions from globs
         if (ilevel.eq.1) then
            remove_bc_nodes = .true.
         else
            remove_bc_nodes = .false.
         end if
         ncornermin = 0 ! do not add anything randomly
         allocate(kglobs(lkglobs),typeglobs(ltypeglobs))
         if (myid.eq.0) then
            call pp_get_globs(ndim,meshdim,nelem,nnod,nsub,&
                              levels(ilevel)%inet,levels(ilevel)%linet,levels(ilevel)%nnet,levels(ilevel)%lnnet,&
                              levels(ilevel)%nndf,levels(ilevel)%lnndf,&
                              levels(ilevel)%xyz,levels(ilevel)%lxyz1,levels(ilevel)%lxyz2,&
                              remove_bc_nodes,levels(ilevel)%ifix,levels(ilevel)%lifix,&
                              levels(ilevel)%iets,levels(ilevel)%liets,ncornermin,&
                              nnodi, ncorner,nedge,nface,&
                              kglobs,lkglobs, typeglobs,ltypeglobs)
         end if
!***************************************************************PARALLEL
         call MPI_BCAST(nnodi,  1,       MPI_INTEGER, 0, comm_all, ierr)
         call MPI_BCAST(ncorner,1,       MPI_INTEGER, 0, comm_all, ierr)
         call MPI_BCAST(nedge,  1,       MPI_INTEGER, 0, comm_all, ierr)
         call MPI_BCAST(nface,  1,       MPI_INTEGER, 0, comm_all, ierr)
         call MPI_BCAST(kglobs,lkglobs,       MPI_INTEGER, 0, comm_all, ierr)
         call MPI_BCAST(typeglobs,ltypeglobs, MPI_INTEGER, 0, comm_all, ierr)
!***************************************************************PARALLEL
         nnodc = ncorner + nedge + nface

         !debug
         !print *,'kglobs',kglobs
         !print *,'typeglobs',typeglobs
         ! create array INODC
         linodc = ncorner
         allocate(inodc(linodc))
         indinodc = 0
         do inod = 1,nnod
            if (typeglobs(inod).eq.3) then
               indinodc = indinodc + 1
               inodc(indinodc) = inod
            end if
         end do
         ! arrays NNGLB and INGLB
         nglb = nedge + nface
         lnnglb = nglb
         linglb = count(typeglobs.eq.1 .or. typeglobs.eq.2)
         allocate(nnglb(lnnglb))
         allocate(inglb(linglb))
         call zero(nnglb,lnnglb)
         call zero(inglb,linglb)
         do inod = 1,nnod
            if (typeglobs(inod).eq.1 .or. typeglobs(inod).eq.2) then
               indglob = kglobs(inod) - ncorner
               ! now use array nnglb as pointers to inglb
               nnglb(indglob) = nnglb(indglob) + 1
            end if
         end do
         lkglb = nglb
         allocate(kglb(lkglb))
         if (nglb.gt.0) then
            kglb(1) = 0
            do iglb = 2,nglb
               kglb(iglb) = kglb(iglb-1) + nnglb(iglb-1)
            end do
         end if
         call zero(nnglb,lnnglb)
         do inod = 1,nnod
            if (typeglobs(inod).eq.1 .or. typeglobs(inod).eq.2) then
               indglob = kglobs(inod) - ncorner
               ! now use array nnglb as pointers to inglb
               nnglb(indglob) = nnglb(indglob) + 1
               inglb(kglb(indglob) + nnglb(indglob)) = inod
            end if
         end do
         deallocate(kglb)
         deallocate(kglobs,typeglobs)
         if (myid.eq.0) then
            write(*,*) ' Corners and globs identified.'
            call flush(6)
         end if
!-----profile
         if (profile) then
            call MPI_BARRIER(comm_all,ierr)
            call time_end(t_globs)
            if (myid.eq.0) then
               write(*,*) '***********************************PROFILING'
               write(*,*) 'Time of generating globs: ',t_globs,' s.'
               write(*,*) '***********************************PROFILING'
               call flush(6)
            end if
         end if
!-----profile
      end if

      ! debug
      !write(*,*) 'nnodc, ncorner, nedge, nface'
      !write(*,*) nnodc, ncorner, nedge, nface
      !write(*,*) 'inodc',inodc
      !write(*,*) 'nnglb',nnglb
      !write(*,*) 'inglb',inglb
      !call flush(6)

!-----profile 
      if (profile) then
         call MPI_BARRIER(comm_all,ierr)
         call time_start
      end if
!-----profile

      ! prepare coarse nnetc
      lnnetc = nsub
      allocate(nnetcaux(lnnetc))
      nnetcaux = 0
      ! create subdomain corners and globs
      do isub_loc = 1,levels(ilevel)%nsub_loc
         call dd_localize_cornersglobs(levels(ilevel)%subdomains(isub_loc),ncorner,inodc,linodc,nedge,nface,&
                                       nnglb,lnnglb,inglb,linglb,nnodcs)
         nnetcaux(isub) = nnodcs
      end do
      allocate(nnetc(lnnetc))
!*****************************************************************MPI
      call MPI_ALLREDUCE(nnetcaux,nnetc,lnnetc, MPI_INTEGER, MPI_SUM, comm_all, ierr) 
!*****************************************************************MPI
      deallocate(nnetcaux)
      ! debug
      !print *,'myid =',myid,'nnetc',nnetc

!-----profile
      if (profile) then
         call MPI_BARRIER(comm_all,ierr)
         call time_end(t_loc_globs)
         if (myid.eq.0) then
            write(*,*) '***********************************PROFILING'
            write(*,*) 'Time of localization of globs: ',t_loc_globs,' s.'
            write(*,*) '***********************************PROFILING'
            call flush(6)
         end if
      end if
!-----profile

!-----profile 
      if (profile) then
         call MPI_BARRIER(comm_all,ierr)
         call time_start
      end if
!-----profile

      ! create subdomain BC and RHS
      do isub_loc = 1,levels(ilevel)%nsub_loc

         ! find size of subdomain
         call dd_get_size(levels(ilevel)%subdomains(isub_loc), ndofs,nnods,nelems)

         ! make subdomain boundary conditions - IFIXS and FIXVS
         lifixs = ndofs
         lfixvs = ndofs
         allocate(ifixs(lifixs),fixvs(lfixvs))
         call dd_map_glob_to_sub(levels(ilevel)%subdomains(isub_loc), levels(ilevel)%ifix,levels(ilevel)%lifix, ifixs,lifixs)
         call dd_map_glob_to_sub(levels(ilevel)%subdomains(isub_loc), levels(ilevel)%fixv,levels(ilevel)%lfixv, fixvs,lfixvs)
         call dd_upload_bc(levels(ilevel)%subdomains(isub_loc), ifixs,lifixs, fixvs,lfixvs)
         deallocate(ifixs,fixvs)

         ! on first level, load also RHS
         if (ilevel.eq.1) then
            lrhss = ndofs
            allocate(rhss(lrhss))
            call dd_map_glob_to_sub(levels(ilevel)%subdomains(isub_loc), levels(ilevel)%rhs,levels(ilevel)%lrhs, rhss,lrhss)
            call dd_upload_rhs(levels(ilevel)%subdomains(isub_loc), rhss,lrhss)
            deallocate(rhss)
         end if
      end do

!-----profile
      if (profile) then
         call MPI_BARRIER(comm_all,ierr)
         call time_end(t_loc_bc)
         if (myid.eq.0) then
            write(*,*) '***********************************PROFILING'
            write(*,*) 'Time of localization of boundary conditions: ',t_loc_bc,' s.'
            write(*,*) '***********************************PROFILING'
            call flush(6)
         end if
      end if
!-----profile

      ! start preparing cnodes
      do isub_loc = 1,levels(ilevel)%nsub_loc
         call dd_construct_cnodes(levels(ilevel)%subdomains(isub_loc))
      end do

! Begin loop over subdomains
! prepare coarse mesh for next level
      linetc = sum(nnetc)
      allocate(inetcaux(linetc))
      inetcaux = 0

      lkinetc = nsub
      allocate(kinetc(lkinetc))
      if (nnodc.gt.0) then
         kinetc(1) = 0
         do inc = 2,nsub
            kinetc(inc) = kinetc(inc-1) + nnetc(inc-1)
         end do
      end if

      do isub_loc = 1,levels(ilevel)%nsub_loc

         call dd_get_coarse_size(levels(ilevel)%subdomains(isub_loc),ndofcs,nnodcs)

         lcnodenumbers = nnodcs
         allocate(cnodenumbers(lcnodenumbers))
         call dd_get_coarse_cnodes(levels(ilevel)%subdomains(isub_loc), cnodenumbers,lcnodenumbers)
         
! create subdomain description of coarse MESH
         pointinetc = kinetc(isub)
         inetcaux(pointinetc+1:pointinetc+nnodcs) = cnodenumbers
         pointinetc = pointinetc + nnodcs

         deallocate(cnodenumbers)

! Finish loop over subdomains
      end do
      deallocate(kinetc)
      allocate(inetc(linetc))
!*****************************************************************MPI
      call MPI_ALLREDUCE(inetcaux,inetc,linetc, MPI_INTEGER, MPI_SUM, comm_all, ierr) 
!*****************************************************************MPI
      deallocate(inetcaux)
      ! debug
      ! print *,'myid =',myid,'inetc',inetc
      ! check the inetc array
      if (any(inetc.eq.0)) then
         call error('LEVELS_PREPARE_STANDARD_LEVEL','Zeros in inetc array.')
      end if

! prepare array of nodes for coarse mesh
      lxyzc1 = nnodc
      lxyzc2 = ndim
      allocate(xyzc(lxyzc1,lxyzc2))

! copy coordinates of corners
      indxyzc = 0
      do inc = 1,ncorner
         indxyzc = indxyzc + 1

         indnc = inodc(inc)
         xyzc(indxyzc,:) = levels(ilevel)%xyz(indnc,:)
      end do
! create coordinates of globs as mean values
      pointinglb = 0
      do iglb = 1,nglb
         indxyzc = indxyzc + 1

         nglbn = nnglb(iglb)

         ! touch first node in glob
         xyzc(indxyzc,:) = sum(levels(ilevel)%xyz(inglb(pointinglb + 1:pointinglb + nglbn),:),dim=1)/nglbn

         pointinglb = pointinglb + nglbn
      end do
      deallocate(nnglb)
      deallocate(inglb)
      deallocate(inodc)
      
      ! Get matrix
!-----profile
      if (profile) then
         call MPI_BARRIER(comm_all,ierr)
         call time_start
      end if
!-----profile
      if (ilevel.eq.1) then
         if (load_division) then
            ! if division was loaded, subdomain files with element matrices should be ready, load them
            do isub_loc = 1,levels(ilevel)%nsub_loc
               call dd_read_matrix_from_file(levels(ilevel)%subdomains(isub_loc),matrixtype,trim(problemname))
            end do
         else
            ! if division of first level was created in the solver, use global file for input of element matrices
            call dd_read_matrix_by_root(levels(ilevel)%subdomains,levels(ilevel)%lsubdomains, comm_all,trim(problemname),&
                                        levels(ilevel)%nsub,levels(ilevel)%nelem,matrixtype,&
                                        levels(ilevel)%sub2proc,levels(ilevel)%lsub2proc,&
                                        levels(ilevel)%indexsub,levels(ilevel)%lindexsub,&
                                        levels(ilevel)%iets,levels(ilevel)%liets)
         end if
      else
         call dd_gather_matrix(levels(ilevel)%subdomains,levels(ilevel)%lsubdomains,&
                               levels(ilevel-1)%subdomains,levels(ilevel-1)%lsubdomains,&
                               comm_all,levels(ilevel)%nsub,levels(ilevel)%nelem,matrixtype,&
                               levels(ilevel)%sub2proc,levels(ilevel)%lsub2proc,&
                               levels(ilevel)%indexsub,levels(ilevel)%lindexsub,&
                               levels(ilevel-1)%sub2proc,levels(ilevel-1)%lsub2proc,&
                               levels(ilevel-1)%indexsub,levels(ilevel-1)%lindexsub,&
                               levels(ilevel)%iets,levels(ilevel)%liets)
      end if
!-----profile
      if (profile) then
         call MPI_BARRIER(comm_all,ierr)
         call time_end(t_matrix_import)
         if (myid.eq.0) then
            write(*,*) '***********************************PROFILING'
            write(*,*) 'Time of importing matrix: ',t_matrix_import,' s.'
            write(*,*) '***********************************PROFILING'
            call flush(6)
         end if
      end if
!-----profile

      ! Assembly matrix for each subdomain
      do isub_loc = 1,levels(ilevel)%nsub_loc
         call dd_assembly_local_matrix(levels(ilevel)%subdomains(isub_loc))
      end do

      ! For first level, prepare Schur complements
      if (ilevel.eq.1) then
         ! Schur only for first level
         remove_original = .false.
         do isub_loc = 1,levels(ilevel)%nsub_loc
            call dd_matrix_tri2blocktri(levels(ilevel)%subdomains(isub_loc),remove_original)
            call dd_prepare_schur(levels(ilevel)%subdomains(isub_loc),comm_self)
         end do
!         call dd_print_sub(myid)
          ! neigbouring should be already created
!         call dd_create_neighbouringi(myid,nsub,comm_all)
      end if

      ! weights
      call dd_weights_prepare(levels(ilevel)%subdomains,levels(ilevel)%lsubdomains, &
                              levels(ilevel)%sub2proc,levels(ilevel)%lsub2proc,&
                              levels(ilevel)%indexsub,levels(ilevel)%lindexsub,&
                              comm_all)

      if (ilevel.eq.1) then
         ! prepare reduced RHS
         call dd_prepare_reduced_rhs(levels(ilevel)%subdomains,levels(ilevel)%lsubdomains, &
                                     levels(ilevel)%sub2proc,levels(ilevel)%lsub2proc,&
                                     levels(ilevel)%indexsub,levels(ilevel)%lindexsub,& 
                                     comm_all)
      end if

! prepare nndfc
      ! in nndfc, nodes are ordered as corners - edges - faces
      lnndfc   = nnodc
      allocate (nndfc(lnndfc))

      nndfc = 0
      ! in corners, prescribe as many constraints as dimension
      nndfc(1:ncorner) = ndim
      if (use_arithmetic) then
         ! on edges
         nndfc(ncorner+1:ncorner+nedge) = ndim
         ! on faces
         if (use_adaptive) then
            nndfc(ncorner+nedge+1:ncorner+nedge+nface) = 0
         else
            nndfc(ncorner+nedge+1:ncorner+nedge+nface) = ndim
         end if
      else
         ! on edges
         nndfc(ncorner+1:ncorner+nedge) = 0
         ! on faces
         nndfc(ncorner+nedge+1:ncorner+nedge+nface) = 0
      end if

      ! BDDC data
      if (ilevel.eq.1) then
         keep_global = .false.
      else
         keep_global = .true.
      end if
      do isub_loc = 1,levels(ilevel)%nsub_loc
         ! load arithmetic averages on edges
         if (use_arithmetic) then
            glbtype = 2
            call dd_load_arithmetic_constraints(levels(ilevel)%subdomains(isub_loc),glbtype)
         end if
         ! load arithmetic averages on faces if adaptivity is not active
         if (.not.use_adaptive) then
            glbtype = 1
            call dd_load_arithmetic_constraints(levels(ilevel)%subdomains(isub_loc),glbtype)
         end if

         ! prepare matrix C for corners and arithmetic averages on edges
         call dd_embed_cnodes(levels(ilevel)%subdomains(isub_loc),nndfc,lnndfc)
         ! prepare matrix C for corners and arithmetic averages on edges
         call dd_prepare_c(levels(ilevel)%subdomains(isub_loc))
         ! prepare augmented matrix for BDDC
         call dd_prepare_aug(levels(ilevel)%subdomains(isub_loc),comm_self)
         ! prepare coarse space basis functions for BDDC
         call dd_prepare_coarse(levels(ilevel)%subdomains(isub_loc),keep_global)
      end do

      if (use_adaptive) then
         ! open file with description of pairs
         if (myid.eq.0) then
            filename = trim(problemname)//'.PAIR'
            call allocate_unit(idpair)
            open (unit=idpair,file=filename,status='old',form='formatted')
         end if

         call adaptivity_init(myid,comm_all,idpair,npair)
   
         call adaptivity_assign_pairs(npair,nproc,npair_locx)
   
         call adaptivity_solve_eigenvectors(levels(ilevel)%subdomains,levels(ilevel)%lsubdomains, &
                                            levels(ilevel)%sub2proc,levels(ilevel)%lsub2proc,&
                                            levels(ilevel)%indexsub,levels(ilevel)%lindexsub,&
                                            comm_all,npair_locx,npair)

         ! update nndf
         call adaptivity_update_ndof(nndfc,lnndfc,ncorner,nedge,nface)
   
         call adaptivity_finalize

         ! prepare AGAIN BDDC data
         do isub_loc = 1,levels(ilevel)%nsub_loc
            ! prepare matrix C for corners and arithmetic averages on edges
            call dd_embed_cnodes(levels(ilevel)%subdomains(isub_loc),nndfc,lnndfc)
            ! prepare matrix C for corners and arithmetic averages on edges
            call dd_prepare_c(levels(ilevel)%subdomains(isub_loc))
            ! prepare augmented matrix for BDDC
            call dd_prepare_aug(levels(ilevel)%subdomains(isub_loc),comm_self)
            ! prepare coarse space basis functions for BDDC
            call dd_prepare_coarse(levels(ilevel)%subdomains(isub_loc),keep_global)
         end do
      end if

      ! print the output
      !do isub_loc = 1,levels(ilevel)%nsub_loc
      !   call dd_print_sub(levels(ilevel)%subdomains(isub_loc))
      !end do

      ! upload coarse mesh
      ndofc = sum(nndfc)
      call levels_upload_level_mesh(ilevel, ncorner,nedge,nface,nnodc, ndofc, &
                                    inetc,linetc,nnetc,lnnetc,nndfc,lnndfc,xyzc,lxyzc1,lxyzc2)
      deallocate(inetc)
      deallocate(nnetc)
      deallocate(nndfc)
      deallocate(xyzc)

      ! prepare memory for coarse solutions
      levels(ilevel)%lsolc = ndofc
      allocate(levels(ilevel)%solc(levels(ilevel)%lsolc))

      levels(ilevel)%is_level_prepared = .true.

!      if (debug) then
!         write(*,*) 'ilevel = ',ilevel
!         write(*,*) 'LEVEL',ilevel-1,', indexsub = ',levels(ilevel-1)%indexsub
!      end if
end subroutine

!***********************************************
subroutine levels_prepare_last_level(matrixtype)
!***********************************************
! Subroutine for building the coarse problem on root process
      use module_mumps
      use module_sm
      use module_utils
      implicit none
      include "mpif.h"

      integer,intent(in) :: matrixtype

      ! local vars
      integer :: ilevel

      integer :: la, nnz
      integer,allocatable :: i_sparse(:),j_sparse(:)
      real(kr),allocatable :: a_sparse(:)

      integer :: ndof

      integer :: mumpsinfo

      ! MPI variables
      integer :: comm_all, myid, ierr

      ! time info
      real(kr) :: t_mumps_analysis, t_mumps_factorization

      ! last level has index of number of levels
      ilevel = nlevels

      ! orient in the communicator
      comm_all  = levels(ilevel)%comm_all
      call MPI_COMM_RANK(comm_all,myid,ierr)

      ! dof
      ndof = levels(ilevel-1)%ndofc

      !write(*,*) 'coarse ndof on myid',myid,'is',ndof

      if (.not.levels(ilevel-1)%is_level_prepared) then
         write(*,*) 'LEVELS_PREPARE_LAST_LEVEL: Previous level not ready.'
         call error_exit
      end if

      ! find length of coarse matrix
      call dd_get_my_coarsem_length(levels(ilevel-1)%subdomains,levels(ilevel-1)%lsubdomains,&
                                    levels(ilevel-1)%indexsub,levels(ilevel-1)%lindexsub,la)

      !write(*,*) 'myid =',myid,'la =',la

! Allocate proper size of matrix A on processor
      allocate(i_sparse(la), j_sparse(la), a_sparse(la))
      i_sparse = 0
      j_sparse = 0
      a_sparse = 0.0_kr

      ! load coarse matrix
      call dd_get_my_coarsem(levels(ilevel-1)%subdomains,levels(ilevel-1)%lsubdomains,matrixtype,&
                             levels(ilevel-1)%indexsub,levels(ilevel-1)%lindexsub, &
                             i_sparse, j_sparse, a_sparse, la)

! Assembly entries in matrix
      call sm_assembly(i_sparse,j_sparse,a_sparse,la,nnz)

      !write(*,*) 'myid =',myid,'la =',la
      !call sm_print(6, i_sparse, j_sparse, a_sparse, la, nnz)

! Initialize MUMPS
      call mumps_init(mumps_coarse,comm_all,matrixtype)

! Level of information from MUMPS
      mumpsinfo  = 1
      call mumps_set_info(mumps_coarse,mumpsinfo)

! Load matrix to MUMPS
      call mumps_load_triplet(mumps_coarse,ndof,nnz,i_sparse,j_sparse,a_sparse,la)

! Analyze matrix
!-----profile
      if (profile) then
         call MPI_BARRIER(comm_all,ierr)
         call time_start
      end if
!-----profile
      call mumps_analyze(mumps_coarse)
      if (myid.eq.0) then
         write(*,*)' Coarse matrix analyzed.'
         call flush(6)
      end if
!-----profile
      if (profile) then
         call MPI_BARRIER(comm_all,ierr)
         call time_end(t_mumps_analysis)
         if (myid.eq.0) then
            write(*,*) '***********************************PROFILING'
            write(*,*) 'Time of analysis of the coarse problem is: ',t_mumps_analysis,' s.'
            write(*,*) '***********************************PROFILING'
            call flush(6)
         end if
      end if
!-----profile

! Factorize matrix
!-----profile
      if (profile) then
         call MPI_BARRIER(comm_all,ierr)
         call time_start
      end if
!-----profile
      call mumps_factorize(mumps_coarse)
      if (myid.eq.0) then
         write(*,*)' Coarse matrix factorized.'
         call flush(6)
      end if
      if (profile) then
         call MPI_BARRIER(comm_all,ierr)
         call time_end(t_mumps_factorization)
         if (myid.eq.0) then
            write(*,*) '***********************************PROFILING'
            write(*,*) 'Time of factorization of the coarse problem is: ',t_mumps_factorization,' s.'
            call flush(6)
            write(*,*) '***********************************PROFILING'
         end if
      end if

      is_mumps_coarse_ready = .true.

! Clear memory
      deallocate(i_sparse, j_sparse, a_sparse)

end subroutine

!***************************************************
subroutine levels_pc_apply(krylov_data,lkrylov_data)
!***************************************************
! Apply last level coarse problem to array lvec
      use module_krylov_types_def
      use module_utils
      implicit none

      integer,intent(in)         ::       lkrylov_data
      type(pcg_data_type),intent(inout) :: krylov_data(lkrylov_data)

      ! local vars
      integer :: ilr
      integer :: comm_all, myid, ierr

! Apply subdomain corrections on first level 
      iactive_level = 1
      comm_all = levels(iactive_level)%comm_all
      ! orient in the communicator
      call MPI_COMM_RANK(comm_all,myid,ierr)
      if (debug .and. myid.eq.0) then
         write(*,*) 'Applying subdomain correction at level',iactive_level
         call flush(6)
      end if
      call levels_corsub_first_level(krylov_data,lkrylov_data)

! Sweep levels upwards - make coarse correction and produce coarse residual
      do iactive_level = 2,nlevels-1
         comm_all = levels(iactive_level)%comm_all
         ! orient in the communicator
         call MPI_COMM_RANK(comm_all,myid,ierr)
         if (debug .and. myid.eq.0) then
            write(*,*) 'Applying subdomain correction at level',iactive_level
            call flush(6)
         end if
         call levels_corsub_standard_level(iactive_level)
      end do

! Solve coarse porblem at last level
      iactive_level = nlevels
      comm_all = levels(iactive_level)%comm_all
      ! orient in the communicator
      call MPI_COMM_RANK(comm_all,myid,ierr)
      if (debug .and. myid.eq.0) then
         write(*,*) 'Solving problem at last level ',iactive_level
         call flush(6)
      end if
      call levels_solve_last_level

! Sweep levels downward - add coarse correction 
      do ilr = 2,nlevels-1
         ! make a reverse loop
         iactive_level = nlevels+1 - ilr

         comm_all = levels(iactive_level)%comm_all
         ! orient in the communicator
         call MPI_COMM_RANK(comm_all,myid,ierr)
         if (debug .and. myid.eq.0) then
            write(*,*) 'Adding coarse correction at level',iactive_level
            call flush(6)
         end if
         call levels_add_standard_level(iactive_level)
      end do

! Add coarse correction at the first level
      iactive_level = 1
      comm_all = levels(iactive_level)%comm_all
      ! orient in the communicator
      call MPI_COMM_RANK(comm_all,myid,ierr)
      if (debug .and. myid.eq.0) then
         write(*,*) 'Adding coarse correction at level',iactive_level
         call flush(6)
      end if
      call levels_add_first_level(krylov_data,lkrylov_data)

end subroutine


!*************************************************************
subroutine levels_corsub_first_level(krylov_data,lkrylov_data)
!*************************************************************
! Apply subdomain correction on first level and produce coarse resudual
      use module_krylov_types_def
      use module_utils
      implicit none
      include "mpif.h"

      integer,intent(in)    ::            lkrylov_data
      type(pcg_data_type),intent(inout) :: krylov_data(lkrylov_data)


      ! local vars
      integer,pointer ::  lresc
      real(kr),pointer ::  resc(:)

      real(kr),allocatable :: rescaux(:)
      integer ::             lrescs
      real(kr),allocatable :: rescs(:)
      integer ::             laux
      real(kr),allocatable :: aux(:)
      integer ::             laux2
      real(kr),allocatable :: aux2(:)

      integer :: lindexsub, ndofs, nnods, nelems, ndofaaugs, ndofcs, nnodcs
      integer :: isub_loc, i, nrhs
      integer :: ilevel
      logical :: transposed

      ! MPI vars
      integer :: comm_all, comm_self, myid, ierr


      ! on first level, resudial is collected from individual subdomains
      ilevel = 1

      ! check prerequisites
      if (.not.levels(ilevel)%is_level_prepared) then
         call error('LEVELS_CORSUB_FIRST_LEVEL','Level is not prepared.')
      end if

      ! set pointers
      lresc => levels(ilevel)%lsolc
      resc  => levels(ilevel)%solc

      ! get communicators
      comm_all  = levels(ilevel)%comm_all
      comm_self = levels(ilevel)%comm_self
      call MPI_COMM_RANK(comm_all,myid,ierr)

      lindexsub = levels(ilevel)%lindexsub
      ! check the length
      if (lkrylov_data .ne. lindexsub) then
         write(*,*) 'LEVELS_CORSUB_FIRST_LEVEL: Inconsistent length of data.'
         call error_exit
      end if

      ! prepare global residual 
      allocate(rescaux(lresc))
      call zero(rescaux,lresc)

      ! get local contribution to coarse residual
      do isub_loc = 1,levels(ilevel)%nsub_loc

         call dd_get_coarse_size(levels(ilevel)%subdomains(isub_loc),ndofcs,nnodcs)

         laux = krylov_data(isub_loc)%lresi
         allocate(aux(laux))
         do i = 1,laux
            aux(i) = krylov_data(isub_loc)%resi(i)
         end do

         ! aux = wi * resi
         call dd_weightsi_apply(levels(ilevel)%subdomains(isub_loc), aux,laux)

         lrescs = ndofcs
         allocate(rescs(lrescs))
         call zero(rescs,lrescs)

         ! rc = phis' * wi * resi
         transposed = .true.
         call dd_phisi_apply(levels(ilevel)%subdomains(isub_loc), transposed, aux,laux, rescs,lrescs)

         ! embed local resc to global one
         call dd_map_subc_to_globc(levels(ilevel)%subdomains(isub_loc), rescs,lrescs, rescaux,lresc)

         ! SUBDOMAIN CORRECTION
         ! prepare array of augmented size
         call dd_get_aug_size(levels(ilevel)%subdomains(isub_loc), ndofaaugs)
         laux2 = ndofaaugs
         allocate(aux2(laux2))
         call zero(aux2,laux2)
         call dd_get_size(levels(ilevel)%subdomains(isub_loc), ndofs,nnods,nelems)
         ! truncate the vector for embedding - zeros at the end
         call dd_map_subi_to_sub(levels(ilevel)%subdomains(isub_loc), aux,laux, aux2,ndofs)

         nrhs = 1
         call dd_solve_aug(levels(ilevel)%subdomains(isub_loc), aux2,laux2, nrhs)

         ! get interface part of the vector of preconditioned residual
         call zero(krylov_data(isub_loc)%z,krylov_data(isub_loc)%lz)
         call dd_map_sub_to_subi(levels(ilevel)%subdomains(isub_loc), aux2,ndofs, &
                                 krylov_data(isub_loc)%z,krylov_data(isub_loc)%lz)

         deallocate(aux2)
         deallocate(aux)
         deallocate(rescs)
      end do
!*****************************************************************MPI
      call MPI_ALLREDUCE(rescaux,resc,lresc, MPI_DOUBLE_PRECISION, MPI_SUM, comm_all, ierr) 
!*****************************************************************MPI
      deallocate(rescaux)
end subroutine

!**********************************************
subroutine levels_corsub_standard_level(ilevel)
!**********************************************
! Apply subdomain correction on standard level and produce coarse resudual
      use module_utils
      implicit none
      include "mpif.h"

      ! index of level
      integer,intent(in) ::  ilevel

      ! local vars
      ! residual on level
      integer,pointer :: lres
      real(kr),pointer :: res(:)

      ! coarse residual on level
      integer,pointer :: lresc
      real(kr),pointer :: resc(:)

      real(kr),allocatable :: resaux(:)
      real(kr),allocatable :: rescaux(:)
      integer ::             lrescs
      real(kr),allocatable :: rescs(:)
      integer ::             laux
      real(kr),allocatable :: aux(:)
      integer ::             laux2
      real(kr),allocatable :: aux2(:)

      integer :: ndofs, nnods, nelems, ndofaaugs, ndofcs, nnodcs
      integer :: isub_loc, i, nrhs
      logical :: transposed

      ! MPI vars
      integer :: comm_all, comm_self, myid, ierr

      ! check LEVEL
      if (ilevel.le.1 .or. ilevel.ge.nlevels) then
         call error('LEVELS_CORSUB_STANDARD_LEVEL','Illegal index of level ILEVEL')
      end if
      ! check prerequisites
      if (.not.levels(ilevel)%is_level_prepared) then
         call error('LEVELS_CORSUB_STANDARD_LEVEL','Current level is not prepared.')
      end if

      ! set pointers
      lres  => levels(ilevel-1)%lsolc
      res   => levels(ilevel-1)%solc
      lresc => levels(ilevel)%lsolc
      resc  => levels(ilevel)%solc

      ! get communicators
      comm_all  = levels(ilevel)%comm_all
      comm_self = levels(ilevel)%comm_self
      call MPI_COMM_RANK(comm_all,myid,ierr)

      ! prepare global coarse residual 
      allocate(rescaux(lresc))
      call zero(rescaux,lresc)

      ! prepare global residual 
      allocate(resaux(lres))
      call zero(resaux,lres)

      ! get local contribution to coarse residual
      do isub_loc = 1,levels(ilevel)%nsub_loc

         ! get dimensions
         call dd_get_size(levels(ilevel)%subdomains(isub_loc),ndofs,nnods,nelems)
         call dd_get_coarse_size(levels(ilevel)%subdomains(isub_loc),ndofcs,nnodcs)

         laux = ndofs
         allocate(aux(laux))
         call dd_map_glob_to_sub(levels(ilevel)%subdomains(isub_loc), res,lres, aux,laux)

         ! aux = wi * resi
         ! weigths are not on interface now
         call dd_weights_apply(levels(ilevel)%subdomains(isub_loc), aux,laux)

         lrescs = ndofcs
         allocate(rescs(lrescs))
         call zero(rescs,lrescs)

         ! rc = phis' * wi * resi
         transposed = .true.
         call dd_phis_apply(levels(ilevel)%subdomains(isub_loc), transposed, aux,laux, rescs,lrescs)

         ! embed local resc to global one
         call dd_map_subc_to_globc(levels(ilevel)%subdomains(isub_loc), rescs,lrescs, rescaux,lresc)

         ! SUBDOMAIN CORRECTION
         ! prepare array of augmented size
         call dd_get_aug_size(levels(ilevel)%subdomains(isub_loc), ndofaaugs)
         laux2 = ndofaaugs
         allocate(aux2(laux2))
         call zero(aux2,laux2)
         call dd_get_size(levels(ilevel)%subdomains(isub_loc), ndofs,nnods,nelems)
         ! truncate the vector for embedding - zeros at the end
         do i = 1,ndofs
            aux2(i) = aux(i)
         end do

         nrhs = 1
         call dd_solve_aug(levels(ilevel)%subdomains(isub_loc), aux2,laux2, nrhs)

         ! get global part of the vector of preconditioned residual
         call dd_map_sub_to_glob(levels(ilevel)%subdomains(isub_loc), aux2,ndofs, resaux,lres)

         deallocate(aux2)
         deallocate(aux)
         deallocate(rescs)
      end do
!*****************************************************************MPI
      call MPI_ALLREDUCE(rescaux,resc,lresc, MPI_DOUBLE_PRECISION, MPI_SUM, comm_all, ierr) 
      call MPI_ALLREDUCE(resaux, res, lres , MPI_DOUBLE_PRECISION, MPI_SUM, comm_all, ierr) 
!*****************************************************************MPI
      deallocate(rescaux)
      deallocate(resaux)
end subroutine

!*********************************
subroutine levels_solve_last_level
!*********************************
! Apply last level coarse problem to array lvec

      use module_mumps

      implicit none
      include "mpif.h"

      integer,pointer  ::  lsolc
      real(kr),pointer ::   solc(:)

      ! local vars
      integer :: ierr, myid, comm_all, comm_self
      integer :: ilevel

      ilevel = nlevels

      comm_all  = levels(ilevel)%comm_all
      comm_self = levels(ilevel)%comm_self
      call MPI_COMM_RANK(comm_all,myid,ierr)

      ! set pointers
      lsolc => levels(ilevel-1)%lsolc
       solc => levels(ilevel-1)%solc

! Solve problem matrix
      call mumps_resolve(mumps_coarse,solc,lsolc)
!***************************************************************PARALLEL
      call MPI_BCAST(solc, lsolc, MPI_DOUBLE_PRECISION, 0, comm_all, ierr)
!***************************************************************PARALLEL

end subroutine

!*******************************************
subroutine levels_add_standard_level(ilevel)
!*******************************************
! Collect residual on first level and produce coarse resudual
      use module_krylov_types_def
      use module_utils
      implicit none
      include "mpif.h"

      ! index of level
      integer,intent(in) ::  ilevel

      ! local vars
      ! residual on level
      integer,pointer :: lsol
      real(kr),pointer :: sol(:)

      ! coarse residual on level
      integer,pointer :: lsolc
      real(kr),pointer :: solc(:)

      integer ::             lsolcs
      real(kr),allocatable :: solcs(:)

      integer ::             lsols
      real(kr),allocatable :: sols(:)

      real(kr),allocatable :: solaux(:)
      real(kr),allocatable :: solaux2(:)

      integer :: ndofcs, nnodcs, nnods, nelems, ndofs
      integer :: isub_loc, i
      logical :: transposed

      ! MPI vars
      integer :: comm_all, comm_self, myid, ierr


      ! check LEVEL
      if (ilevel.le.1 .or. ilevel.ge.nlevels) then
         call error('LEVELS_CORSUB_STANDARD_LEVEL','Illegal index of level ILEVEL')
      end if
      ! check prerequisites
      if (.not.levels(ilevel)%is_level_prepared) then
         call error('LEVELS_APPLY_FIRST_LEVEL','Level is not prepared.')
      end if

      ! set pointers
      lsol => levels(ilevel-1)%lsolc
      sol  => levels(ilevel-1)%solc
      lsolc => levels(ilevel)%lsolc
      solc  => levels(ilevel)%solc

      ! get communicators
      comm_all  = levels(ilevel)%comm_all
      comm_self = levels(ilevel)%comm_self
      call MPI_COMM_RANK(comm_all,myid,ierr)

      ! prepare memory for coarse contribution
      allocate(solaux(lsol))
      allocate(solaux2(lsol))
      call zero(solaux,lsol)
      call zero(solaux2,lsol)

      do isub_loc = 1,levels(ilevel)%nsub_loc

         call dd_get_coarse_size(levels(ilevel)%subdomains(isub_loc), ndofcs,nnodcs)
         call dd_get_size(levels(ilevel)%subdomains(isub_loc), ndofs,nnods,nelems)

         lsolcs = ndofcs
         allocate(solcs(lsolcs))

         lsols  = ndofs
         allocate(sols(lsols))
         call zero(sols,lsols)

         ! restrict global solc to local solcs
         call dd_map_globc_to_subc(levels(ilevel)%subdomains(isub_loc), solc,lsolc, solcs,lsolcs)

         ! COARSE CORRECTION
         ! z_i = z_i + phis_i * uc_i
         transposed = .false.
         call dd_phis_apply(levels(ilevel)%subdomains(isub_loc), transposed, solcs,lsolcs, sols,lsols)
         ! apply weights
         ! z = wi * z
         call dd_weights_apply(levels(ilevel)%subdomains(isub_loc), sols,lsols)

         call dd_map_sub_to_glob(levels(ilevel)%subdomains(isub_loc), sols,lsols, solaux,lsol)

         deallocate(solcs)
         deallocate(sols)
      end do
!*****************************************************************MPI
      call MPI_ALLREDUCE(solaux, solaux2, lsol , MPI_DOUBLE_PRECISION, MPI_SUM, comm_all, ierr) 
!*****************************************************************MPI
      deallocate(solaux)
      
      ! add corrections - sol already contains subdomain contributions
      do i = 1,lsol
         sol(i) = sol(i) + solaux2(i)
      end do
      deallocate(solaux2)

end subroutine
!**********************************************************
subroutine levels_add_first_level(krylov_data,lkrylov_data)
!**********************************************************
! Collect residual on first level and produce coarse resudual
      use module_krylov_types_def
      use module_utils
      implicit none
      include "mpif.h"

      integer,intent(in)    ::            lkrylov_data
      type(pcg_data_type),intent(inout) :: krylov_data(lkrylov_data)

      ! local vars
      integer,pointer :: lsolc
      real(kr),pointer :: solc(:)
      integer ::             lsolcs
      real(kr),allocatable :: solcs(:)

      integer :: lindexsub, ndofcs, nnodcs
      integer :: isub_loc, i
      integer :: ilevel
      logical :: transposed

      ! MPI vars
      integer :: comm_all, comm_self, myid, ierr


      ! on first level, resudial is collected from individual subdomains
      ilevel = 1

      ! check prerequisites
      if (.not.levels(ilevel)%is_level_prepared) then
         call error('LEVELS_ADD_FIRST_LEVEL','Level is not prepared.')
      end if

      ! set pointers
      lsolc => levels(ilevel)%lsolc
      solc  => levels(ilevel)%solc

      ! get communicators
      comm_all  = levels(ilevel)%comm_all
      comm_self = levels(ilevel)%comm_self
      call MPI_COMM_RANK(comm_all,myid,ierr)

      lindexsub = levels(ilevel)%lindexsub
      ! check the length
      if (lkrylov_data .ne. lindexsub) then
         write(*,*) 'LEVELS_ADD_FIRST_LEVEL: Inconsistent length of data.'
         call error_exit
      end if

      do isub_loc = 1,levels(ilevel)%nsub_loc

         call dd_get_coarse_size(levels(ilevel)%subdomains(isub_loc), ndofcs,nnodcs)

         lsolcs = ndofcs
         allocate(solcs(lsolcs))

         ! restrict global solc to local solcs
         call dd_map_globc_to_subc(levels(ilevel)%subdomains(isub_loc), solc,lsolc, solcs,lsolcs)

         ! COARSE CORRECTION
         ! z_i = z_i + phis_i * uc_i
         transposed = .false.
         call dd_phisi_apply(levels(ilevel)%subdomains(isub_loc), transposed, solcs,lsolcs, &
                             krylov_data(isub_loc)%z,krylov_data(isub_loc)%lz)
         ! apply weights
         ! z = wi * z
         call dd_weightsi_apply(levels(ilevel)%subdomains(isub_loc), krylov_data(isub_loc)%z,krylov_data(isub_loc)%lz)

         ! load Z for communication
         call dd_comm_upload(levels(ilevel)%subdomains(isub_loc), krylov_data(isub_loc)%z,krylov_data(isub_loc)%lz)

         deallocate(solcs)
      end do
      ! communicate vector Z
      ! Interchange data
      call dd_comm_swapdata(levels(ilevel)%subdomains,levels(ilevel)%lsubdomains,&
                            levels(ilevel)%indexsub,levels(ilevel)%lindexsub,&    
                            levels(ilevel)%sub2proc,levels(ilevel)%lsub2proc,&
                            comm_all)
      ! Download data
      do isub_loc = 1,levels(ilevel)%nsub_loc

         call zero(krylov_data(isub_loc)%zadj,krylov_data(isub_loc)%lz)
         ! get contibution from neigbours
         call dd_comm_download(levels(ilevel)%subdomains(isub_loc), krylov_data(isub_loc)%zadj,krylov_data(isub_loc)%lz)
         ! join data
         do i = 1,krylov_data(isub_loc)%lz
            krylov_data(isub_loc)%z(i) = krylov_data(isub_loc)%z(i) + krylov_data(isub_loc)%zadj(i)
         end do

      end do

end subroutine

!***************************************************************
subroutine levels_get_number_of_subdomains(ilevel,nsub,nsub_loc)
!***************************************************************
! Subroutine for getting number of subdomains at levels
      use module_utils
      implicit none

! given index of level
      integer,intent(in) :: ilevel
! number of subdomains in that level
      integer,intent(out) :: nsub
      integer,intent(out) :: nsub_loc

      ! check prerequisites
      if (.not.levels(ilevel)%is_level_prepared) then
         call error('LEVELS_GET_NUMBER_OF_SUBDOMAINS','Level is not prepared.')
      end if

      nsub     = levels(ilevel)%nsub
      nsub_loc = levels(ilevel)%nsub_loc

end subroutine

!**************************************************************
subroutine levels_prepare_krylov_data(krylov_data,lkrylov_data)
!**************************************************************
! Subroutine for initialization of data for Krylov iterative method
! module for distributed Krylov data storage
      use module_krylov_types_def
! module with utility routines
      use module_utils

      implicit none

      integer,intent(in)               :: lkrylov_data
      type (pcg_data_type),intent(out) ::  krylov_data(lkrylov_data)

! local vars
      ! for Krylov data, index of level is 1
      integer,parameter :: ilevel = 1

      integer :: isub_loc
      integer :: ndofs, nnods, nelems, ndofis, nnodis
      integer :: lsolis, lresis, laps, lps, lzs

      integer ::             lsols
      real(kr),allocatable :: sols(:)

      ! check prerequisites
      if (.not.levels(ilevel)%is_level_prepared) then
         call error('LEVELS_PREPARE_KRYLOV_DATA','Level is not prepared.')
      end if
      ! check dimensions
      if (lkrylov_data.ne.levels(ilevel)%nsub_loc) then
         call error('LEVELS_PREPARE_KRYLOV_DATA','Dimension mismatch.')
      end if

      do isub_loc = 1,levels(ilevel)%nsub_loc
         ! determine size of subdomain
         call dd_get_size(levels(ilevel)%subdomains(isub_loc),ndofs,nnods,nelems)
         call dd_get_interface_size(levels(ilevel)%subdomains(isub_loc),ndofis,nnodis)

         ! allocate local subdomain solution
         lsols  = ndofs
         allocate(sols(lsols))

         if (levels(ilevel)%use_initial_solution) then
            ! set initial guess
            call dd_map_glob_to_sub(levels(ilevel)%subdomains(isub_loc),levels(ilevel)%sol,levels(ilevel)%lsol,sols,lsols)
         else
            ! set zero initial guess
            ! u_0 = 0
            call zero(sols,lsols)
         end if

         ! set initial guess to satisfy Dirichlet boundary conditions
         call dd_fix_bc(levels(ilevel)%subdomains(isub_loc), sols,lsols)

         ! allocate vectors for Krylov method 
         lsolis = ndofis
         krylov_data(isub_loc)%lsoli = lsolis
         allocate(krylov_data(isub_loc)%soli(lsolis))
         call zero(krylov_data(isub_loc)%soli,lsolis)
         lresis = ndofis
         krylov_data(isub_loc)%lresi = lresis
         allocate(krylov_data(isub_loc)%resi(lresis))
         call zero(krylov_data(isub_loc)%resi,lresis)

         ! restrict solution to interface
         call dd_map_sub_to_subi(levels(ilevel)%subdomains(isub_loc), sols,lsols, krylov_data(isub_loc)%soli,lsolis)
         deallocate(sols)

         ! set initial residual to RHS
         ! res = g
         call dd_get_reduced_rhs(levels(ilevel)%subdomains(isub_loc), krylov_data(isub_loc)%resi,lresis)

         laps = ndofis
         krylov_data(isub_loc)%lap = laps
         allocate(krylov_data(isub_loc)%ap(laps))
         allocate(krylov_data(isub_loc)%apadj(laps))

         lps = ndofis
         krylov_data(isub_loc)%lp = lps
         allocate(krylov_data(isub_loc)%p(lps))
         allocate(krylov_data(isub_loc)%padj(lps))

         lresis = krylov_data(isub_loc)%lresi
         allocate(krylov_data(isub_loc)%resiadj(lresis))

         lzs = ndofis
         krylov_data(isub_loc)%lz = lzs
         allocate(krylov_data(isub_loc)%z(lzs))
         allocate(krylov_data(isub_loc)%zadj(lzs))

      end do

end subroutine

!****************************************************************
subroutine levels_fix_bc_interface_dual(krylov_data,lkrylov_data)
!****************************************************************
! Subroutine for fixing boundary conditions in subdomain right hand sides 
! module for distributed Krylov data storage
      use module_krylov_types_def
! module with utility routines
      use module_utils

      implicit none

      integer,intent(in)                 :: lkrylov_data
      type (pcg_data_type),intent(inout) ::  krylov_data(lkrylov_data)

! local vars
      ! for Krylov data, index of level is 1
      integer,parameter :: ilevel = 1
      integer :: isub_loc

      ! check prerequisites
      if (.not.levels(ilevel)%is_level_prepared) then
         call error('LEVELS_FIX_BC_INTERFACE_DUAL','Level is not prepared.')
      end if
      ! check dimensions
      if (lkrylov_data.ne.levels(ilevel)%nsub_loc) then
         call error('LEVELS_FIX_BC_INTERFACE_DUAL','Dimension mismatch.')
      end if

      do isub_loc = 1,levels(ilevel)%nsub_loc
         ! determine size of subdomain
         call dd_fix_bc_interface_dual(levels(ilevel)%subdomains(isub_loc), &
                                       krylov_data(isub_loc)%resi,krylov_data(isub_loc)%lresi)
      end do

end subroutine

!***************************************************
subroutine levels_sm_apply(krylov_data,lkrylov_data)
!***************************************************
! Subroutine for postprocessong of data by Krylov iterative method
! module for distributed Krylov data storage
      use module_krylov_types_def
! module with utility routines
      use module_utils

      implicit none

      integer,intent(in)                 :: lkrylov_data
      type (pcg_data_type),intent(inout) ::  krylov_data(lkrylov_data)

! local vars
      ! for Krylov data, index of level is 1
      integer,parameter :: ilevel = 1
      integer :: i

      integer :: comm_all
      integer :: isub_loc

      ! check prerequisites
      if (.not.levels(ilevel)%is_level_prepared) then
         call error('LEVELS_POSTPROCESS_SOLUTION','Level is not prepared.')
      end if
      ! check dimensions
      if (lkrylov_data.ne.levels(ilevel)%nsub_loc) then
         call error('LEVELS_POSTPROCESS_SOLUTION','Dimension mismatch.')
      end if

      ! set communicator
      comm_all = levels(ilevel)%comm_all

      ! ap = A*p
      ! Upload data
      do isub_loc = 1,levels(ilevel)%nsub_loc

         call zero(krylov_data(isub_loc)%ap,krylov_data(isub_loc)%lap)

         call dd_multiply_by_schur(levels(ilevel)%subdomains(isub_loc),&
                                   krylov_data(isub_loc)%p,krylov_data(isub_loc)%lp, &
                                   krylov_data(isub_loc)%ap,krylov_data(isub_loc)%lap)

         call dd_comm_upload(levels(ilevel)%subdomains(isub_loc), krylov_data(isub_loc)%ap,krylov_data(isub_loc)%lap)
      end do
      ! Interchange data
      call dd_comm_swapdata(levels(ilevel)%subdomains,levels(ilevel)%lsubdomains,&
                            levels(ilevel)%indexsub,levels(ilevel)%lindexsub,&    
                            levels(ilevel)%sub2proc,levels(ilevel)%lsub2proc,&
                            comm_all)
      ! Download data
      do isub_loc = 1,levels(ilevel)%nsub_loc

         call zero(krylov_data(isub_loc)%apadj,krylov_data(isub_loc)%lap)
         ! get contibution from neigbours
         call dd_comm_download(levels(ilevel)%subdomains(isub_loc), krylov_data(isub_loc)%apadj,krylov_data(isub_loc)%lap)
         ! join data
         do i = 1,krylov_data(isub_loc)%lap
            krylov_data(isub_loc)%ap(i) = krylov_data(isub_loc)%ap(i) + krylov_data(isub_loc)%apadj(i)
         end do
      end do
end subroutine

!******************************************************************************************
subroutine levels_postprocess_solution(krylov_data,lkrylov_data,problemname,print_solution)
!******************************************************************************************
! Subroutine for postprocessong of data by Krylov iterative method
! module for distributed Krylov data storage
      use module_krylov_types_def
! module for handling subdomain data
      use module_dd
! module with utility routines
      use module_utils

      implicit none

      integer,intent(in)              :: lkrylov_data
      type (pcg_data_type),intent(in) ::  krylov_data(lkrylov_data)
      character(*),intent(in) :: problemname
      logical, intent(in) :: print_solution

! local vars
      ! for Krylov data, index of level is 1
      integer,parameter :: ilevel = 1
      integer :: isub_loc, isub
      integer :: ndofs, nnods, nelems

      integer ::             lsols
      real(kr),allocatable :: sols(:)

      ! check prerequisites
      if (.not.levels(ilevel)%is_level_prepared) then
         call error('LEVELS_POSTPROCESS_SOLUTION','Level is not prepared.')
      end if
      ! check dimensions
      if (lkrylov_data.ne.levels(ilevel)%nsub_loc) then
         call error('LEVELS_POSTPROCESS_SOLUTION','Dimension mismatch.')
      end if

      do isub_loc = 1,levels(ilevel)%nsub_loc
         isub = levels(ilevel)%indexsub(isub_loc)

         ! determine size of subdomain
         call dd_get_size(levels(ilevel)%subdomains(isub_loc),ndofs,nnods,nelems)

         ! allocate local subdomain solution
         lsols  = ndofs
         allocate(sols(lsols))

         ! set zero solution
         call zero(sols,lsols)

         ! resolve interior values
         call dd_resolve_interior(levels(ilevel)%subdomains(isub_loc), &
                                  krylov_data(isub_loc)%soli,krylov_data(isub_loc)%lsoli, &
                                  sols,lsols)

         ! write subdomain solution to independent disk files
         call dd_write_solution_to_file(trim(problemname),isub,sols,lsols,print_solution)

         deallocate(sols)
      end do

end subroutine

!**************************************************************
subroutine levels_destroy_krylov_data(krylov_data,lkrylov_data)
!**************************************************************
! Subroutine for clearing data for Krylov iterative method
! module for distributed Krylov data storage
      use module_krylov_types_def
! module with utility routines
      use module_utils

      implicit none

      integer,intent(in)                 :: lkrylov_data
      type (pcg_data_type),intent(inout) ::  krylov_data(lkrylov_data)

! local vars
      ! for Krylov data, index of level is 1
      integer,parameter :: ilevel = 1

      integer :: isub_loc

      ! check prerequisites
      if (.not.levels(ilevel)%is_level_prepared) then
         call error('LEVELS_DESTROY_KRYLOV_DATA','Level is not prepared.')
      end if
      ! check dimensions
      if (lkrylov_data.ne.levels(ilevel)%nsub_loc) then
         call error('LEVELS_DESTROY_KRYLOV_DATA','Dimension mismatch.')
      end if

      do isub_loc = 1,levels(ilevel)%nsub_loc
         deallocate(krylov_data(isub_loc)%ap)
         deallocate(krylov_data(isub_loc)%apadj)
         deallocate(krylov_data(isub_loc)%p)
         deallocate(krylov_data(isub_loc)%padj)
         deallocate(krylov_data(isub_loc)%resiadj)
         deallocate(krylov_data(isub_loc)%z)
         deallocate(krylov_data(isub_loc)%zadj)
         deallocate(krylov_data(isub_loc)%soli)
         deallocate(krylov_data(isub_loc)%resi)
      end do

end subroutine

!***********************************************************************************
subroutine levels_dd_dotprod_local(ilevel,isub_loc, vec1,lvec1, vec2,lvec2, dotprod)
!***********************************************************************************
! Subroutine used for indirect access to DD data
! multiplies vector 1 and vector 2 using weights in DD
      implicit none

      ! length of vector
      integer,intent(in) ::   ilevel 
      integer,intent(in) ::   isub_loc 
      ! vectors to multiply
      integer,intent(in) ::  lvec1
      real(kr), intent(in) :: vec1(lvec1)
      integer,intent(in) ::  lvec2
      real(kr), intent(in) :: vec2(lvec2)
      
      ! result
      real(kr), intent(out) :: dotprod

      ! add data from module and call function from adaptive module
      call dd_dotprod_local(levels(ilevel)%subdomains(isub_loc), vec1,lvec1, vec2,lvec2, dotprod)

end subroutine

!********************************************************
subroutine levels_adaptivity_mvecmult(n,x,lx,y,ly,idoper)
!********************************************************
! Subroutine used inside generation of adaptive constraints,
! helps eigensolver with multiplication of vector x by local Schur complement
! and storing result in y.
! Local matrix is accessed by module
      use module_adaptivity
      implicit none

      ! length of vector
      integer,intent(in) ::   n 
      integer,intent(in) ::  lx 
      real(kr),intent(in) ::  x(lx)
      integer,intent(in) ::  ly 
      real(kr),intent(out) :: y(ly)
      integer,intent(inout) ::  idoper 

      ! add data from module and call function from adaptive module
      call adaptivity_mvecmult(levels(iactive_level)%subdomains,levels(iactive_level)%lsubdomains,&
                               levels(iactive_level)%indexsub,levels(iactive_level)%lindexsub,&
                               n,x,lx,y,ly,idoper)
end subroutine

!************************************
subroutine levels_clear_level(level)
!************************************
! Subroutine for deallocation data of one level
      implicit none

! local variables
      type(levels_type) :: level
      integer :: isub_loc

! deallocate data
      if (allocated(level%inetc)) then
         deallocate (level%inetc)
      end if
      level%linet = 0
      if (allocated(level%nnetc)) then
         deallocate (level%nnetc)
      end if
      level%lnnet = 0
      if (allocated(level%nndfc)) then
         deallocate (level%nndfc)
      end if
      level%lnndf = 0
      if (allocated(level%iets)) then
         deallocate (level%iets)
      end if
      level%liets = 0
      if (allocated(level%xyzc)) then
         deallocate (level%xyzc)
      end if
      level%lxyz1 = 0
      level%lxyz2 = 0
      if (allocated(level%ifixc)) then
         deallocate (level%ifixc)
      end if
      level%lifixc = 0
      if (allocated(level%fixvc)) then
         deallocate (level%fixvc)
      end if
      level%lfixvc = 0
      if (allocated(level%rhsc)) then
         deallocate (level%rhsc)
      end if
      level%lrhsc = 0
      if (allocated(level%sub2proc)) then
         deallocate (level%sub2proc)
      end if
      level%lsub2proc = 0
      if (allocated(level%indexsub)) then
         deallocate (level%indexsub)
      end if
      level%lindexsub = 0
      if (allocated(level%solc)) then
         deallocate (level%solc)
      end if
      level%lsolc = 0
      if (allocated(level%subdomains)) then
         do isub_loc = 1,level%lsubdomains
            call dd_finalize(level%subdomains(isub_loc))
         end do
         deallocate (level%subdomains)
      end if
      level%lsubdomains = 0

      level%nelem = 0
      level%nnod  = 0
      level%ndof  = 0

      level%nsub  = 0

      level%ncorner  = 0
      level%nedge    = 0
      level%nface    = 0
      level%nnodc    = 0

      if (associated(level%inet)) then
         nullify(level%inet)
      end if
      level%linet = 0
      if (associated(level%nnet)) then
         nullify(level%nnet)
      end if
      level%lnnet = 0
      if (associated(level%nndf)) then
         nullify(level%nndf)
      end if
      level%lnndf = 0
      if (associated(level%xyz)) then
         nullify(level%xyz)
      end if
      level%lxyz1 = 0
      level%lxyz2 = 0

end subroutine

!*************************
subroutine levels_finalize
!*************************
! Subroutine for initialization of levels data
      use module_mumps
      implicit none

! local variables
      integer :: ilevel, start(1)
      
! destroy MUMPS structure of the last level
      if (is_mumps_coarse_ready) then
         call mumps_finalize(mumps_coarse)
      end if

! deallocate basic structure
      if (allocated(levels)) then
         ! clear data of particular subdomains
         start = lbound(levels)
         do ilevel = start(1),llevels
            call levels_clear_level(levels(ilevel))
         end do

         deallocate (levels)
      end if


end subroutine

end module module_levels

