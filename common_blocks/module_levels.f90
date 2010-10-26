module module_levels
!*******************
! Module for handling levels in multilevel BDDC 
! Jakub Sistek, Praha 2/2010

!     definition of MUMPS structure
      use dmumps_struc_def

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

         ! subdomain data
         integer ::              nsub     ! number of subdomains on level
         integer ::              nsub_loc ! number of subdomains on level assigned to the processor

         integer             :: lindexsub    
         integer,allocatable ::  indexsub(:) ! indices of elements in array sub
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

         logical :: is_level_prepared = .false. ! level prepared

      end type levels_type

      integer, private ::                          nlevels
      integer, private ::                          llevels
      type(levels_type), allocatable, target, private ::    levels(:)

      logical,private :: is_mumps_coarse_ready = .false.
      type(DMUMPS_STRUC), private :: mumps_coarse  

contains

!******************************
subroutine levels_init(nl,nsub)
!******************************
! Subroutine for initialization of levels data
      use module_dd
      implicit none

! given number of levels
      integer,intent(in) :: nl
! number of subdomains in all levels
      integer,intent(in) :: nsub

! initialize basic structure
      nlevels = nl
      llevels = nlevels
      allocate (levels(llevels))

! initialize DD module
      call dd_init(nsub)

end subroutine

!***************************************************************************************
subroutine levels_init_with_zero_level(nl,nsublev,lnsublev,nelem,nnod,ndof,&
                                       inet,linet,nnet,lnnet,nndf,lnndf,xyz,lxyz1,lxyz2,&
                                       ifix,lifix, fixv,lfixv, rhs,lrhs)
!***************************************************************************************
! Subroutine for initialization of levels data
      use module_dd
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

      ! local vars 
      integer :: ilevel, i, j, nsuball

! initialize DD module for all subdomains at all levels
      nsuball = sum(nsublev)
      call dd_init(nsuball)

! initialize basic structure
      nlevels = nl
      llevels = nlevels
      allocate (levels(0:llevels))

      ! debug
      ilevel = 0
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

      levels(ilevel)%is_level_prepared = .true.

end subroutine

!*************************************************************************************
subroutine levels_pc_setup(problemname,nsublev,lnsublev,use_preprocessor,correct_division,&
                           neighbouring,matrixtype,ndim, meshdim, use_arithmetic, use_adaptive)
!*************************************************************************************
! subroutine for multilevel BDDC preconditioner setup
      use module_pp
      use module_dd
      implicit none
      include "mpif.h"

! name of problem
      character(*),intent(in) :: problemname
! number of subdomains in all levels
      integer,intent(in) :: lnsublev
      integer,intent(in) ::  nsublev(lnsublev)
! was preprocessor run before the solver?
      logical,intent(in) :: use_preprocessor
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
      integer :: ilevel, indsub, isub, lindexsub, nsub, nsub_loc, nsuball_loc
      integer :: nproc, comm_all, ierr, myid

      ! no debug
      !ilevel = 1
      !call levels_read_level_from_file(problemname,comm_all,ndim,ilevel)

      ! load number of subdomains for each level
      indsub = 0
      do ilevel = 1,nlevels
         nsub = nsublev(ilevel)
         levels(ilevel)%nsub = nsub
      ! prepare array of indices of subdomains for each level
         lindexsub = nsub
         levels(ilevel)%lindexsub = lindexsub
         allocate(levels(ilevel)%indexsub(lindexsub))
         ! associate numbers of subdomains in array sub with this level
         do isub = 1,nsub
            indsub = indsub + 1

            levels(ilevel)%indexsub(isub) = indsub
         end do

      end do

      ! TODO: create proper communicators
      do ilevel = 1,nlevels
         levels(ilevel)%comm_all  = MPI_COMM_WORLD
         levels(ilevel)%comm_self = MPI_COMM_SELF
      end do

      nsuball_loc = 0
      do ilevel = 1,nlevels
         nsub = levels(ilevel)%nsub 

         comm_all = levels(ilevel)%comm_all
         ! orient in the communicator
         call MPI_COMM_RANK(comm_all,myid,ierr)
         call MPI_COMM_SIZE(comm_all,nproc,ierr)

         levels(ilevel)%lsub2proc = nproc + 1
         allocate(levels(ilevel)%sub2proc(levels(ilevel)%lsub2proc))
         call pp_distribute_subdomains(nsub,nproc,levels(ilevel)%sub2proc,levels(ilevel)%lsub2proc)

         call dd_distribute_subdomains(levels(ilevel)%indexsub(1),levels(ilevel)%indexsub(nsub),&
                                       levels(ilevel)%sub2proc,levels(ilevel)%lsub2proc, nproc)

         nsub_loc    = levels(ilevel)%sub2proc(myid+2) - levels(ilevel)%sub2proc(myid+1)
         levels(ilevel)%nsub_loc = nsub_loc

         nsuball_loc = nsuball_loc + nsub_loc
      end do
      ! TODO: initialize dd module only with subdomains local to processor and index them using lindexsub array

      ! associate subdomains with first level
      do ilevel = 1,nlevels-1
         call levels_prepare_standard_level(problemname,use_preprocessor,&
                                            correct_division,neighbouring,matrixtype,ndim,meshdim,ilevel,&
                                            use_arithmetic,use_adaptive)
      end do

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
      integer :: ncorner, nedge, nface, nnodc
      integer :: nsub

      integer ::             linetc,   lnnetc 
      integer, allocatable :: inetc(:), nnetc(:)

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
      lxyzc1 = nnodc
      lxyzc2 = ndim

      allocate(inetc(linetc),nnetc(lnnetc),xyzc(lxyzc1,lxyzc2))

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

      ! save mesh into structure
      levels(ilevel)%nelem = nelem
      levels(ilevel)%nnod  = nnod
      levels(ilevel)%linet = linet
      call levels_upload_level_mesh(ilevel,&
                                    ncorner,nedge,nface,nnodc,&
                                    inetc,linetc,nnetc,lnnetc,xyzc,lxyzc1,lxyzc2)
      deallocate(inetc,nnetc,xyzc)

      if (myid.eq.0) then
         close(idlevel)
      end if

end subroutine

!***********************************************************************************************
subroutine levels_upload_level_mesh(ilevel,ncorner,nedge,nface,nnodc,&
                                    inetc,linetc,nnetc,lnnetc,xyzc,lxyzc1,lxyzc2)
!***********************************************************************************************
! Subroutine for loading mesh of level into levels structure

      use module_utils
      implicit none

      integer,intent(in) :: ilevel
      integer,intent(in) :: ncorner, nedge, nface, nnodc
      integer,intent(in) :: linetc
      integer,intent(in) ::  inetc(linetc)
      integer,intent(in) :: lnnetc
      integer,intent(in) ::  nnetc(lnnetc)
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

      allocate(levels(ilevel)%inetc(linetc))
      do i = 1,linetc
         levels(ilevel)%inetc(i) = inetc(i)
      end do

      levels(ilevel)%lnnetc  = lnnetc
      allocate(levels(ilevel)%nnetc(lnnetc))
      do i = 1,lnnetc
         levels(ilevel)%nnetc(i) = nnetc(i)
      end do

      levels(ilevel)%lxyzc1 = lxyzc1
      levels(ilevel)%lxyzc2 = lxyzc2
      allocate(levels(ilevel)%xyzc(lxyzc1,lxyzc2))
      do j = 1,lxyzc2
         do i = 1,lxyzc1
            levels(ilevel)%xyzc(i,j) = xyzc(i,j)
         end do
      end do

      levels(ilevel)%is_level_prepared = .true.

end subroutine

!*****************************************************************************************
subroutine levels_prepare_standard_level(problemname,use_preprocessor,&
                                         correct_division,neighbouring,matrixtype,ndim,meshdim,ilevel,&
                                         use_arithmetic,use_adaptive)
!*****************************************************************************************
! Subroutine for building the standard level
      use module_pp
      use module_adaptivity
      use module_dd
      use module_sm
      use module_utils
      implicit none
      include "mpif.h"

      character(*),intent(in) :: problemname
      logical,intent(in) :: use_preprocessor
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
      integer :: i, comm_all, comm_self, ierr
      integer :: graphtype 
      integer :: ides, idcn, idglb

      integer :: ncorner, nedge, nface, isub, nnodc, ndofc, nelem, nnod, nnodi
      integer :: edgecut, ncornermin
      integer :: indsub, nsub, glbtype, nglb, nglbn, nglbv

      integer ::            linets,   lnnets,   lnndfs,   lifixs,   lkdofs,   lisegns,   lisngns
      integer,allocatable :: inets(:), nnets(:), nndfs(:), ifixs(:), kdofs(:), isegns(:), isngns(:)
      integer ::            ie, inods 
      integer :: ndofs, nelems, nnods, neighball
      integer :: pkadjsub, j, indadjs, isubneib, jsub, nadjs, ncorners, nglobs, ndofis, ndofn, &
                 ndofos

      integer ::           lkglobs
      integer,allocatable:: kglobs(:)
      integer ::           ltypeglobs
      integer,allocatable:: typeglobs(:)
      integer :: indg, indglob, indinodc, indins, indivs, indn, indn1

      integer ::           lkadjsub
      integer,allocatable:: kadjsub(:)

      integer ::           lkdof
      integer,allocatable:: kdof(:)

      integer ::            liins,   liivsvns,   liovsvns
      integer,allocatable :: iins(:), iivsvns(:), iovsvns(:)
      integer ::            lglobal_corner_numbers,   licnsins
      integer,allocatable :: global_corner_numbers(:), icnsins(:)
      integer::            lglobal_glob_numbers
      integer,allocatable:: global_glob_numbers(:)
      integer::            lglob_types
      integer,allocatable:: glob_types(:)
      integer::            lnglobnodess
      integer,allocatable:: nglobnodess(:)
      integer::            lnglobvars
      integer,allocatable:: nglobvars(:)
      integer::            liadjs
      integer,allocatable:: iadjs(:)
      integer::            lignsins1, lignsins2
      integer,allocatable:: ignsins(:,:)
      integer::            ligvsivns1, ligvsivns2
      integer,allocatable:: igvsivns(:,:)
      integer ::            nelemsadj, nnodsadj
      integer ::            linetsadj,   lnnetsadj,    lisegnsadj,    lisngnsadj
      integer,allocatable :: inetsadj(:), nnetsadj(:),  isegnsadj(:),  isngnsadj(:)
      integer ::            lishared
      integer,allocatable :: ishared(:)
      integer ::            lkinodes
      integer,allocatable :: kinodes(:)
      integer ::             nshared

      integer ::            linodc
      integer,allocatable :: inodc(:)
      integer ::            lnnglb
      integer,allocatable :: nnglb(:)
      integer ::            linglb
      integer,allocatable :: inglb(:)
      
      integer ::            lrhss,   lfixvs,   lxyzs1, lxyzs2
      real(kr),allocatable:: rhss(:), fixvs(:), xyzs(:,:)

      integer ::            linetc
      integer,allocatable :: inetc(:)
      integer ::            lnnetc
      integer,allocatable :: nnetc(:)
      integer ::            lnndfc
      integer,allocatable :: nndfc(:)
      integer ::            lxyzc1, lxyzc2
      real(kr),allocatable:: xyzc(:,:)

      integer ::            lcnodenumbers
      integer,allocatable :: cnodenumbers(:)

      integer :: idofn, idofis, idofos, iglb, iglbn, iglbv, iglobs, inc, indnc, &
                 indng, indns, indvs, inod, inodcs, inodis, nnodis, pointinglb
      integer :: indifix, indifixs
      integer :: indrhs, indrhss
      integer :: indxyzc, ndofcs, nnodcs, pointinetc


      logical :: remove_original 
      logical :: remove_bc_nodes 

      ! adaptivity variables
      integer :: idpair
      character(200) :: filename
      integer :: npair, npair_locx

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

      ! initialize values
      nsub  = levels(ilevel)%nsub
      nnod  = levels(ilevel)%nnod
      nelem = levels(ilevel)%nelem

      ! make division into subdomains
      graphtype = 0 ! not weighted
      levels(ilevel)%liets = levels(ilevel)%nelem
      allocate(levels(ilevel)%iets(levels(ilevel)%liets))
      if (ilevel.eq.1 .and. use_preprocessor) then
         ! read the division from file *.ES
         if (myid.eq.0) then
            filename = trim(problemname)//'.ES'
            call allocate_unit(ides)
            open (unit=ides,file=filename,status='old',form='formatted')
            rewind ides

            read(ides,*) levels(ilevel)%iets
            close (ides)
            write(*,*) 'myid = ',myid,'Mesh division loaded from file.'
            call flush(6)
         end if
      else
         ! create division
         ! TODO: make this parallel 
         if (myid.eq.0) then
            call pp_divide_mesh(graphtype,correct_division,neighbouring,&
                                levels(ilevel)%nelem,levels(ilevel)%nnod,&
                                levels(ilevel)%inet,levels(ilevel)%linet,levels(ilevel)%nnet,levels(ilevel)%lnnet,nsub,&
                                edgecut,levels(ilevel)%iets,levels(ilevel)%liets)
            write(*,*) 'myid = ',myid,'Mesh division created.'
            call flush(6)
         end if 
      end if
!***************************************************************PARALLEL
      call MPI_BCAST(levels(ilevel)%iets,levels(ilevel)%liets, MPI_INTEGER, 0, comm_all, ierr)
!***************************************************************PARALLEL

      if (ilevel.eq.1 .and. use_preprocessor) then
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
      else
         ! TODO: make this parallel
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
         write(*,*) 'myid = ',myid,'Corners and globs identified.'
         call flush(6)
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
         ! debug
         !write(*,*) 'inodc',inodc
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
               inglb(nnglb(indglob)) = inod
            end if
         end do
         deallocate(kglobs,typeglobs)
      end if
      !write(*,*) 'nnglb',nnglb
      !write(*,*) 'inglb',inglb

      ! for communication, any shared nodes are considered
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

      ! create array KDOF
      lkdof = nnod
      allocate(kdof(lkdof))
      kdof(1) = 0
      do inod = 2,nnod
         kdof(inod) = kdof(inod-1) + levels(ilevel)%nndf(inod-1)
      end do

      ! TODO: for local subdomains prepare local mesh
      ! prepare coarse nnetc
      lnnetc = nsub
      allocate(nnetc(lnnetc))
      nnetc = 0
      do isub = levels(ilevel)%sub2proc(myid+1), levels(ilevel)%sub2proc(myid+2)-1
         nelems = 0
         linets = 0
         do ie = 1,levels(ilevel)%nelem
            if (levels(ilevel)%iets(ie).eq.isub) then
               nelems = nelems + 1
               linets = linets + levels(ilevel)%nnet(ie)
            end if
         end do
         lnnets  = nelems
         lisegns = nelems
         lisngns = linets
         allocate(inets(linets),nnets(lnnets),isegns(lisegns),isngns(lisngns))

         call pp_create_submesh(isub,levels(ilevel)%nelem,levels(ilevel)%inet,levels(ilevel)%linet,&
                                levels(ilevel)%nnet,levels(ilevel)%lnnet,&
                                levels(ilevel)%iets,levels(ilevel)%liets,&
                                nnods,inets,linets,nnets,lnnets,&
                                isegns,lisegns,isngns,lisngns)
         lnndfs = nnods
         lkdofs = nnods
         allocate(nndfs(lnndfs),kdofs(lkdofs))

         ! get array nndfs
         do inods = 1,nnods
            nndfs(inods) = levels(ilevel)%nndf(isngns(inods))
         end do
! find local number of DOF on subdomain NDOFS
         ndofs = sum(nndfs)

! create array kdofs
         kdofs(1) = 0
         do inods = 2,nnods
            kdofs(inods) = kdofs(inods-1) + nndfs(inods-1)
         end do

         ! coordinates of subdomain nodes
         lxyzs1 = nnods
         lxyzs2 = ndim
         allocate(xyzs(lxyzs1,lxyzs2))
         do j = 1,ndim
            do i = 1,nnods
               indn = isngns(i)

               xyzs(i,j) = levels(ilevel)%xyz(indn,j)
            end do
         end do

         ! find number of adjacent subdomains NADJS
         pkadjsub = (isub-1)*nsub
         nadjs = 0
         do jsub = 1,nsub
            if (kadjsub(pkadjsub + jsub).eq.1) then
               nadjs = nadjs + 1
            end if
         end do
         ! debug
         !print *,'nadjs',nadjs
         ! create list of neighbouring subdomains
         liadjs = nadjs
         allocate(iadjs(liadjs))
         indadjs = 0
         do jsub = 1,nsub
            if (kadjsub(pkadjsub + jsub).eq.1) then
               indadjs = indadjs + 1
               iadjs(indadjs) = jsub
            end if
         end do
         ! debug
         !print *,'iadjs',iadjs

         lkinodes = nnods
         allocate(kinodes(lkinodes))
         kinodes = 0
         lishared = nnods
         allocate(ishared(lishared))
         ! this could be improved using graphs of mesh
         ! TODO: perform this in parallel 
         do isubneib = 1,nadjs
            jsub = iadjs(isubneib)

! create subdomain description of MESH
            nelemsadj = 0
            linetsadj = 0
            do ie = 1,nelem
               if (levels(ilevel)%iets(ie).eq.jsub) then
                  nelemsadj = nelemsadj + 1
                  linetsadj = linetsadj + levels(ilevel)%nnet(ie)
               end if
            end do
            lnnetsadj  = nelemsadj
            lisegnsadj = nelemsadj
            lisngnsadj = linetsadj
            allocate(inetsadj(linetsadj),nnetsadj(lnnetsadj),isegnsadj(lisegnsadj),isngnsadj(lisngnsadj))
   
            call pp_create_submesh(jsub,nelem,levels(ilevel)%inet,levels(ilevel)%linet,&
                                   levels(ilevel)%nnet,levels(ilevel)%lnnet,levels(ilevel)%iets,levels(ilevel)%liets,&
                                   nnodsadj,inetsadj,linetsadj,nnetsadj,lnnetsadj,&
                                   isegnsadj,lisegnsadj,isngnsadj,lisngnsadj)


            call get_array_intersection(isngns,nnods,isngnsadj,nnodsadj,ishared,lishared,nshared)
            ! check that some nodes are really shared
            if (nshared.eq.0) then
               write(*,*) 'Error in shared nodes -  seems as zero shared nodes with adjacent subdomain.'
               call error_exit
            end if

            ! mark interface nodes
            do i = 1,nshared
               indg = ishared(i)
               call get_index(indg,isngns,nnods,inods)
               kinodes(inods) = kinodes(inods) + 1
            end do
            deallocate(inetsadj,nnetsadj,isegnsadj,isngnsadj)
         end do
         deallocate(ishared)

! find number of interface nodes
         nnodis = 0 
         ndofis = 0
         ndofos = 0
         do inods = 1,nnods

            if (kinodes(inods).gt.0) then
               nnodis = nnodis + 1
               ndofis = ndofis + nndfs(inods)
            else
               ndofos = ndofos + nndfs(inods)
            end if
         end do
! generate mapping of interface nodes to subdomain nodes and the same for dofs 
         liins = nnodis
         allocate(iins(liins))
         liivsvns = ndofis
         allocate(iivsvns(liivsvns))
         liovsvns = ndofos
         allocate(iovsvns(liovsvns))

         inodis = 0
         idofis = 0
         idofos = 0
         do inods = 1,nnods
            ndofn  = nndfs(inods)

            if (kinodes(inods).gt.0) then
               inodis = inodis + 1

               iins(inodis) = inods
               do idofn = 1,ndofn 
                  idofis = idofis + 1

                  iivsvns(idofis) = kdofs(inods) + idofn
               end do
            else
               do idofn = 1,ndofn 
                  idofos = idofos + 1

                  iovsvns(idofos) = kdofs(inods) + idofn
               end do
            end if
         end do


         ! find number of coarse nodes on subdomain NCORNERS
         ncorners = 0
         do inc = 1,ncorner
            indnc = inodc(inc)
            if (any(isngns(1:nnods).eq.indnc)) then
               ncorners = ncorners + 1
            end if
         end do

         ! find mapping of corners
         lglobal_corner_numbers = ncorners
         allocate(global_corner_numbers(lglobal_corner_numbers))
         licnsins = ncorners
         allocate(icnsins(licnsins))

         inodcs = 0
         do inc = 1,ncorner
            indnc = inodc(inc)
            if (any(isngns(1:nnods).eq.indnc)) then
               inodcs = inodcs + 1

               ! mapping to global corner numbers
               global_corner_numbers(inodcs) = inc
               ! mapping to subdomain interface numbers
               call get_index(indnc,isngns,nnods,indns)
               if (indns .eq. -1) then
                  write(*,*) 'CREATE_SUB_FILES: Index of subdomain node not found.', indnc
                  call error_exit
               end if
               call get_index(indns,iins,liins,indins)
               if (indins .eq. -1) then
                  write(*,*) 'CREATE_SUB_FILES: Index of subdomain interface node not found.', indns
                  write(*,*) 'iins',iins
                  call error_exit
               end if
               icnsins(inodcs) = indins
            end if
         end do


         ! find local number of globs NGLOBS
         nglobs     = 0
         pointinglb = 0
         do iglb = 1,nglb
            nglbn = nnglb(iglb)

            ! touch first node in glob
            indn1 = inglb(pointinglb + 1)

            if (any(isngns(1:nnods).eq.indn1)) then
               nglobs = nglobs + 1
            end if

            pointinglb = pointinglb + nglbn
         end do

         ! mapping of globs
         lglobal_glob_numbers = nglobs
         allocate(global_glob_numbers(lglobal_glob_numbers))
         lnglobvars = nglobs
         allocate(nglobvars(lnglobvars))
         lnglobnodess = nglobs
         allocate(nglobnodess(lnglobnodess))
         lglob_types = nglobs
         allocate(glob_types(lglob_types))

         iglobs     = 0
         pointinglb = 0
         do iglb = 1,nglb
            nglbn = nnglb(iglb)

            ! touch first node in glob
            indn1 = inglb(pointinglb + 1)

            if (any(isngns(1:nnods).eq.indn1)) then

               iglobs = iglobs + 1

               nglbv = 0
               do iglbn = 1,nglbn
                  ndofn = levels(ilevel)%nndf(pointinglb + iglbn)

                  nglbv = nglbv + ndofn
               end do

               nglobvars(iglobs)  = nglbv
               nglobnodess(iglobs) = nglbn
               global_glob_numbers(iglobs) = iglb
            end if

            pointinglb = pointinglb + nglbn
         end do

         ! set type of glob
         glob_types = 1
         where (global_glob_numbers .le. nedge) glob_types = 2

         ! shift numbering behind corners
         global_glob_numbers = global_glob_numbers + ncorner

         ligvsivns1 = nglobs
         ligvsivns2 = maxval(nglobvars)
         allocate(igvsivns(ligvsivns1,ligvsivns2))
         lignsins1 = nglobs
         lignsins2 = maxval(nglobnodess)
         allocate(ignsins(lignsins1,lignsins2))
         iglobs     = 0
         pointinglb = 0
         do iglb = 1,nglb
            nglbn = nnglb(iglb)

            ! touch first node in glob
            indn1 = inglb(pointinglb + 1)

            if (any(isngns(1:nnods).eq.indn1)) then

               iglobs = iglobs + 1

               iglbv = 0
               do iglbn = 1,nglbn
                  
                  indng = inglb(pointinglb + iglbn)
                  ndofn = levels(ilevel)%nndf(indng)

                  call get_index(indng,isngns,nnods,indns)
                  if (indns .eq. -1) then
                     write(*,*) 'CREATE_SUB_FILES: Index of subdomain node not found.', indng
                     call error_exit
                  end if
                  call get_index(indns,iins,nnodis,indins)
                  if (indins .eq. -1) then
                     write(*,*) 'CREATE_SUB_FILES: Index of interface node not found.', indns
                     call error_exit
                  end if
                  ignsins(iglobs,iglbn) = indins

                  do idofn = 1,ndofn
                     iglbv = iglbv + 1

                     indvs = kdofs(indns) + idofn
                     call get_index(indvs,iivsvns,liivsvns,indivs)
                     if (indivs .eq. -1) then
                        write(*,*) 'CREATE_SUB_FILES: Index of subdomain interface dof not found.'
                        write(*,*) 'indng =',indng,'indns =',indns,'indvs = ',indvs,'indivs = ',indivs, 'isub = ',isub
                        write(*,*) 'iivsvns = ',iivsvns
                        stop
                     end if

                     igvsivns(iglobs,iglbv) = indivs
                  end do
               end do
            end if

            pointinglb = pointinglb + nglbn
         end do

         ! make subdomain boundary conditions - IFIXS and FIXVS
         lifixs = ndofs
         lfixvs = ndofs
         allocate(ifixs(lifixs),fixvs(lfixvs))
         ifixs = 0
         fixvs = 0._kr
         indifixs = 0
         do inods = 1,nnods
            inod = isngns(inods)

            ndofn   = nndfs(inods)
            indifix = kdof(inod)
            do idofn = 1,ndofn
               indifixs = indifixs + 1
               indifix  = indifix + 1
               if (levels(ilevel)%ifix(indifix).ne.0) then
                  ifixs(indifixs) = indifixs
                  fixvs(indifixs) = levels(ilevel)%fixv(indifix)
               end if
            end do
         end do

         ! on first level, load also RHS
         if (ilevel.eq.1) then
            lrhss = ndofs
            allocate(rhss(lrhss))
            indrhss = 0
            do inods = 1,nnods
               inod = isngns(inods)

               ndofn   = nndfs(inods)
               indrhs  = kdof(inod)
               do idofn = 1,ndofn
                  indrhss = indrhss + 1
                  indrhs  = indrhs  + 1
                  rhss(indrhss) = levels(ilevel)%rhs(indrhs)
               end do
            end do
         end if

         ! load data to structure
         call dd_upload_mesh(myid,isub, nelems, nnods, ndofs, ndim, nnodis, ndofis, ndofos, ncorners, nglobs, nadjs,&
                             nndfs,lnndfs, nnets,lnnets, inets,linets, isngns,lisngns, isegns,lisegns,&
                             xyzs,lxyzs1,lxyzs2, &
                             iins,liins, iivsvns,liivsvns, iovsvns,liovsvns,&
                             global_corner_numbers,lglobal_corner_numbers, icnsins,licnsins,&
                             global_glob_numbers,lglobal_glob_numbers,nglobnodess,lnglobnodess, nglobvars,lnglobvars,&
                             ignsins,lignsins1,lignsins2, igvsivns,ligvsivns1,ligvsivns2,glob_types,lglob_types, iadjs,liadjs)
         call dd_load_bc(myid,isub, ifixs,lifixs, fixvs,lfixvs)

         ! start preparing cnodes
         call dd_get_cnodes(myid,isub)

         ! on first level, upload also right hand side
         if (ilevel.eq.1) then
            call dd_load_rhs(myid,isub, rhss,lrhss)
            deallocate(rhss)
            ! for first level, load matrix from files
            ! TODO: parallel input
            call dd_read_matrix_from_file(myid,isub,trim(problemname),matrixtype)
         end if

         deallocate(ifixs,fixvs)
         deallocate(ignsins)
         deallocate(igvsivns)
         deallocate(glob_types)
         deallocate(nglobnodess)
         deallocate(nglobvars)
         deallocate(global_glob_numbers)
         deallocate(icnsins)
         deallocate(global_corner_numbers)
         deallocate(iovsvns)
         deallocate(iivsvns)
         deallocate(iins)
         deallocate(kinodes)
         deallocate(iadjs)
         deallocate(nndfs,kdofs)
         deallocate(xyzs)
         deallocate(inets,nnets,isegns,isngns)

         nnetc(isub) = ncorners + nglobs

      end do
      deallocate(kadjsub)

      ! Prepare file for second level
      ! nnetc  = number of coarse nodes that touch each subdomain

      linetc = sum(nnetc)

! Begin loop over subdomains
! prepare array for level 2
      allocate(inetc(linetc))
      inetc = 0

      pointinetc = 0
      do isub = levels(ilevel)%sub2proc(myid+1), levels(ilevel)%sub2proc(myid+2)-1
         indsub = levels(ilevel)%indexsub(isub)

         call dd_get_coarse_size(myid,indsub,ndofcs,nnodcs)
         lcnodenumbers = nnodcs
         allocate(cnodenumbers(lcnodenumbers))
         call dd_get_coarse_cnodes(myid,indsub,cnodenumbers,lcnodenumbers)
         
! create subdomain description of coarse MESH
         inetc(pointinetc+1:pointinetc+nnodcs) = cnodenumbers
         pointinetc = pointinetc + nnodcs

         deallocate(cnodenumbers)

! Finish loop over subdomains
      end do
      ! check the inetc array
      if (any(inetc.eq.0)) then
         write(*,*) 'Zeros in inetc array.'
         call error_exit
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
      deallocate(kdof)
      
      ! upload coarse mesh
      call levels_upload_level_mesh(ilevel, ncorner,nedge,nface,nnodc,&
                                    inetc,linetc,nnetc,lnnetc,xyzc,lxyzc1,lxyzc2)
      deallocate(inetc)
      deallocate(nnetc)
      deallocate(xyzc)

      do isub = 1,nsub
         indsub = levels(ilevel)%indexsub(isub)
         call dd_assembly_local_matrix(myid,indsub)
      end do
      if (ilevel.eq.1) then
         ! Schur only for first level
         remove_original = .false.
         do isub = 1,nsub
            indsub = levels(ilevel)%indexsub(isub)
            call dd_matrix_tri2blocktri(myid,indsub,remove_original)
            call dd_prepare_schur(myid,comm_self,indsub)
         end do
!         call dd_print_sub(myid)
         call dd_create_neighbouring(myid,nsub,comm_all)

         ! weights
         call dd_weights_prepare(myid, nsub, comm_all)

         ! reduced RHS
         call dd_prepare_reduced_rhs(myid,nsub, comm_all)

!         call dd_print_sub(myid)
      end if

! prepare nndfc
      ! in nndfc, nodes are ordered as corners - edges - faces
      ncorner  = levels(ilevel)%ncorner
      nedge    = levels(ilevel)%nedge
      nface    = levels(ilevel)%nface
      nnodc    = ncorner + nedge + nface
      ! check sum
      if (nnodc .ne. ncorner + nedge + nface) then
         write(*,*) 'LEVELS_PREPARE_STANDARD_LEVEL: coarse nodes number mismatch. nnodc = ',nnodc,&
                    'sum of globs',ncorner + nedge + nface
         call error_exit
      end if
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
      ! load arithmetic averages on edges
      if (use_arithmetic) then
         glbtype = 2
         do isub = 1,nsub
            indsub = levels(ilevel)%indexsub(isub)
            call dd_load_arithmetic_constraints(myid,indsub,glbtype)
         end do
         ! load arithmetic averages on faces
         if (.not.use_adaptive) then
            glbtype = 1
            do isub = 1,nsub
               indsub = levels(ilevel)%indexsub(isub)
               call dd_load_arithmetic_constraints(myid,indsub,glbtype)
            end do
         end if
      end if
      ! prepare matrix C for corners and arithmetic averages on edges
      do isub = 1,nsub
         indsub = levels(ilevel)%indexsub(isub)
         call dd_embed_cnodes(myid,indsub,nndfc,lnndfc)
         call dd_prepare_c(myid,indsub)
      end do

      ! prepare augmented matrix for BDDC
      do isub = 1,nsub
         indsub = levels(ilevel)%indexsub(isub)
         call dd_prepare_aug(myid,comm_self,indsub)
      end do

      ! prepare coarse space basis functions for BDDC
      do isub = 1,nsub
         indsub = levels(ilevel)%indexsub(isub)
         call dd_prepare_coarse(myid,indsub)
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
   
         call adaptivity_solve_eigenvectors(myid,comm_all,npair_locx,npair,nproc)

         ! update nndf
         call adaptivity_update_ndof(nndfc,lnndfc,ncorner,nedge,nface)
   
         call adaptivity_finalize

         ! prepare AGAIN matrix C, now for corners, arithmetic averages on edges and adaptive on faces
         do isub = 1,nsub
            indsub = levels(ilevel)%indexsub(isub)
            call dd_embed_cnodes(myid,indsub,nndfc,lnndfc)
            call dd_prepare_c(myid,indsub)
         end do

         ! prepare augmented matrix for BDDC
         do isub = 1,nsub
            indsub = levels(ilevel)%indexsub(isub)
            call dd_prepare_aug(myid,comm_self,indsub)
         end do

         ! prepare coarse space basis functions for BDDC
         do isub = 1,nsub
            indsub = levels(ilevel)%indexsub(isub)
            call dd_prepare_coarse(myid,indsub)
         end do

      end if

      ! print the output
!      call dd_print_sub(myid)

! load coarse NNDF array to structure
      ndofc = sum(nndfc)
      levels(ilevel)%ndofc = ndofc
      ! load nndf to next level
      levels(ilevel)%lnndfc = lnndfc
      allocate(levels(ilevel)%nndfc(lnndfc))
      do i = 1,lnndfc
         levels(ilevel)%nndfc(i) = nndfc(i)
      end do
      deallocate(nndfc)

!      if (debug) then
!         write(*,*) 'ilevel = ',ilevel
!         write(*,*) 'LEVEL',ilevel-1,', indexsub = ',levels(ilevel-1)%indexsub
!      end if
end subroutine

!***********************************************
subroutine levels_prepare_last_level(matrixtype)
!***********************************************
! Subroutine for building the coarse problem on root process

      use module_dd
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
      integer :: comm_all

      ! MPI variables
      integer :: myid, ierr

      ! time info
      real(kr) :: t_mumps_analysis, t_mumps_factorization

      ! last level has index of number of levels
      ilevel = nlevels

      ! dof
      ndof = levels(ilevel-1)%ndofc

      ! orient in the communicator
      comm_all = levels(ilevel)%comm_all
      call MPI_COMM_RANK(comm_all,myid,ierr)

      !write(*,*) 'coarse ndof on myid',myid,'is',ndof

      if (.not.levels(ilevel-1)%is_level_prepared) then
         write(*,*) 'LEVELS_PREPARE_LAST_LEVEL: Previous level not ready.'
         call error_exit
      end if

      ! find length of coarse matrix
      call dd_get_my_coarsem_length(myid,levels(ilevel-1)%indexsub,levels(ilevel-1)%lindexsub,la)

      !write(*,*) 'myid =',myid,'la =',la

! Allocate proper size of matrix A on processor
      allocate(i_sparse(la), j_sparse(la), a_sparse(la))
      i_sparse = 0
      j_sparse = 0
      a_sparse = 0.0_kr

      ! load coarse matrix
      call dd_get_my_coarsem(myid,matrixtype,levels(ilevel-1)%indexsub,levels(ilevel-1)%lindexsub, &
                             i_sparse, j_sparse, a_sparse, la)

! Assembly entries in matrix
      call sm_assembly(i_sparse,j_sparse,a_sparse,la,nnz)

      !write(*,*) 'myid =',myid,'la =',la
      !call sm_print(6, i_sparse, j_sparse, a_sparse, la, nnz)

! Initialize MUMPS
      call mumps_init(mumps_coarse,comm_all,matrixtype)
      write(*,*)'myid =',myid,': MUMPS Initialized'
      call flush(6)

! Level of information from MUMPS
      mumpsinfo  = 1
      call mumps_set_info(mumps_coarse,mumpsinfo)

! Load matrix to MUMPS
      call mumps_load_triplet(mumps_coarse,ndof,nnz,i_sparse,j_sparse,a_sparse,la)
      write(*,*)'myid =',myid,': Triplet loaded'
      call flush(6)

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

!******************************************************************
subroutine levels_pc_apply(myid,comm_all, krylov_data,lkrylov_data)
!******************************************************************
! Apply last level coarse problem to array lvec
      use module_krylov_types_def
      use module_dd
      use module_utils
      implicit none
      include "mpif.h"

      ! my ID
      integer,intent(in) :: myid
      ! MPI communicator
      integer,intent(in) :: comm_all

      integer,intent(in)    ::       lkrylov_data
      type(pcg_data_type),intent(inout) :: krylov_data(lkrylov_data)

      ! local vars
      integer ::             lsolc
      real(kr),allocatable :: solc(:)
      real(kr),allocatable :: solcaux(:)
      integer ::             lrescs
      real(kr),allocatable :: rescs(:)
      integer ::             laux
      real(kr),allocatable :: aux(:)
      integer ::             laux2
      real(kr),allocatable :: aux2(:)
      integer ::             lsolcs
      real(kr),allocatable :: solcs(:)

      integer :: lindexsub, ndofs, nnods, nelems, ndofaaugs, ndofcs, nnodcs
      integer :: is, i, isub, nrhs
      integer :: ilevel
      logical :: transposed

      ! MPI vars
      integer :: ierr

! Apply subdomain corrections on first level 
      ilevel = 1

      lindexsub = levels(ilevel)%lindexsub
      ! check the length
      if (lkrylov_data .ne. lindexsub) then
         write(*,*) 'LEVELS_PC_APPLY: Inconsistent length of data.'
         call error_exit
      end if

      ! prepare global residual 
      lsolc = levels(ilevel)%ndofc
      allocate(solcaux(lsolc))
      allocate(solc(lsolc))
      call zero(solcaux,lsolc)

      ! get local contribution to coarse residual
      do is = 1,lindexsub
         isub = levels(ilevel)%indexsub(is)

         call dd_get_coarse_size(myid,isub,ndofcs,nnodcs)
         if (ndofcs.le.0) then
            cycle
         end if

         laux = krylov_data(isub)%lresi
         allocate(aux(laux))
         do i = 1,laux
            aux(i) = krylov_data(isub)%resi(i)
         end do

         ! aux = wi * resi
         call dd_weights_apply(myid, isub, aux,laux)

         lrescs = ndofcs
         allocate(rescs(lrescs))
         call zero(rescs,lrescs)

         ! rc = phis' * wi * resi
         transposed = .true.
         call dd_phis_apply(myid,isub, transposed, aux,laux, rescs,lrescs)

         ! embed local resc to global one
         call dd_map_subc_to_globc(myid,isub, rescs,lrescs, solcaux,lsolc)

         ! SUBDOMAIN CORRECTION
         ! prepare array of augmented size
         call dd_get_aug_size(myid,isub, ndofaaugs)
         laux2 = ndofaaugs
         allocate(aux2(laux2))
         call zero(aux2,laux2)
         call dd_get_size(myid,isub, ndofs,nnods,nelems)
         ! truncate the vector for embedding - zeros at the end
         call dd_map_subi_to_sub(myid,isub, aux,laux, aux2,ndofs)

         nrhs = 1
         call dd_solve_aug(myid,isub, aux2,laux2, nrhs)

         ! get interface part of the vector of preconditioned residual
         call zero(krylov_data(isub)%z,krylov_data(isub)%lz)
         call dd_map_sub_to_subi(myid,isub, aux2,ndofs, krylov_data(isub)%z,krylov_data(isub)%lz)


         deallocate(aux2)
         deallocate(aux)
         deallocate(rescs)
      end do

      ! communicate coarse residual along processes
!***************************************************************PARALLEL
      call MPI_REDUCE(solcaux,solc,lsolc, MPI_DOUBLE_PRECISION, MPI_SUM, 0, comm_all, ierr) 
!***************************************************************PARALLEL
! Solve problem matrix
      call levels_apply_last_level(solc,lsolc)
!***************************************************************PARALLEL
      call MPI_BCAST(solc, lsolc, MPI_DOUBLE_PRECISION, 0, comm_all, ierr)
!***************************************************************PARALLEL

      do is = 1,lindexsub
         isub = levels(ilevel)%indexsub(is)

         call dd_get_coarse_size(myid,isub,ndofcs,nnodcs)
         if (ndofcs.le.0) then
            cycle
         end if

         lsolcs = ndofcs
         allocate(solcs(lsolcs))

         ! restrict global solc to local solcs
         call dd_map_globc_to_subc(myid,isub, solc,lsolc, solcs,lsolcs)

         ! COARSE CORRECTION
         ! z_i = z_i + phis_i * uc_i
         transposed = .false.
         call dd_phis_apply(myid,isub, transposed, solcs,lsolcs, krylov_data(isub)%z,krylov_data(isub)%lz)
         ! apply weights
         ! z = wi * z
         call dd_weights_apply(myid, isub, krylov_data(isub)%z,krylov_data(isub)%lz)

         ! load Z for communication
         call dd_comm_upload(myid, isub,  krylov_data(isub)%z,krylov_data(isub)%lz)

         deallocate(solcs)
      end do
      ! communicate vector Z
      ! Interchange data
      call dd_comm_swapdata(myid, lindexsub, comm_all)
      ! Download data
      do is = 1,lindexsub
         isub = levels(ilevel)%indexsub(is)

         if (krylov_data(isub)%is_mine) then

            call zero(krylov_data(isub)%zadj,krylov_data(isub)%lz)
            ! get contibution from neigbours
            call dd_comm_download(myid, isub,  krylov_data(isub)%zadj,krylov_data(isub)%lz)
            ! join data
            do i = 1,krylov_data(isub)%lz
               krylov_data(isub)%z(i) = krylov_data(isub)%z(i) + krylov_data(isub)%zadj(i)
            end do

         end if
      end do

      ! clear local memory
      deallocate(solc)
      deallocate(solcaux)


end subroutine



!*******************************************
subroutine levels_apply_last_level(vec,lvec)
!*******************************************
! Apply last level coarse problem to array lvec

      use module_mumps

      implicit none
      include "mpif.h"

      integer,intent(in)    ::  lvec
      real(kr),intent(inout) ::  vec(lvec)

! Solve problem matrix
      call mumps_resolve(mumps_coarse,vec,lvec)

end subroutine

!************************************
subroutine levels_clear_level(level)
!************************************
! Subroutine for deallocation data of one level
      implicit none

! local variables
      type(levels_type) :: level

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
      use module_dd
      use module_mumps
      implicit none

! local variables
      integer :: ilevel, start(1)
      
! destroy MUMPS structure of the last level
      if (is_mumps_coarse_ready) then
         call mumps_finalize(mumps_coarse)
      end if

! finalize DD module
      call dd_finalize

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

