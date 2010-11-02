module module_dd
!***************
! Module for handling domain decomposition structures
! Jakub Sistek, Denver, 4/2009, Praha 1/2010

!     definition of MUMPS structure
      use dmumps_struc_def
      implicit none

! type of real variables
      integer,parameter,private :: kr = kind(1.D0)
! numerical zero
      real(kr),parameter,private :: numerical_zero = 1.e-12_kr
! maximal allowed length of file names
      integer,parameter,private :: lfnamex = 130

! debugging 
      logical,parameter,private :: debug = .false.

      interface dd_map_glob_to_sub
         module procedure dd_map_glob_to_sub_int
         module procedure dd_map_glob_to_sub_real
      end interface dd_map_glob_to_sub

! type for storing adaptive constraints for selected globs
      type cnode_type
         integer :: itype                ! type of coarse node (1 - face, 2 - edge, 3 - corner)
         logical :: used       = .false. ! should this glob be applied in computation?
         logical :: adaptive   = .false. ! should adaptive constraints be applied at this glob?
         logical :: arithmetic = .false. ! should arithmetic constraints be applied at this glob?
         integer ::             lxyz
         real(kr),allocatable :: xyz(:) ! coordinates
         ! where the glob maps to
         integer :: global_cnode_number ! embedding into coarse nodes
         integer :: ncdof              ! number of coarse degrees of freedom
         integer,allocatable :: igcdof(:) ! indices of global coarse dof
         ! where it maps from
         integer :: nnod               ! number of nodes it contains
         integer,allocatable :: insin(:)  ! indices of glob nodes in subdomain interface numbering
         integer :: nvar               ! number of variables it contains
         integer,allocatable :: ivsivn(:)  ! indices of glob variables in subdomain interface numbering
         ! constraints on glob
         integer :: nnz      ! number of nonzeros the matrix contains
         integer :: lmatrix1 ! = ncdof
         integer :: lmatrix2 ! = nvar
         real(kr),allocatable :: matrix(:,:)
      end type cnode_type

! type for subdomain data
      type subdomain_type

         logical ::             is_initialized = .false.
         integer ::             isub    ! subdomain index
         integer ::             nsub    ! total number of subdomains
         integer ::             proc    ! index of processor taking care of this subdomain
         integer ::             comm    ! communicator to which proc refer

         logical ::             is_mesh_loaded = .false.
         integer ::             nelem   ! number of elements
         integer ::             nnod    ! number of nodes
         integer ::             ndof    ! number of degrees of freedom
         integer ::             ndim    ! dimension of the problem
         ! description of subdomain mesh
         integer ::             linet    ! length of INET array 
         integer,allocatable ::  inet(:) ! INET array - indices of nodes on elements
         integer ::             lnnet    ! length of NNET array
         integer,allocatable ::  nnet(:) ! NNET array - number of nodes on elements
         integer ::             lnndf    ! length of NNDF array
         integer,allocatable ::  nndf(:) ! NNDF array - number of nodal degrees of freedom
         integer ::             lisngn   ! length of array ISNGN
         integer,allocatable ::  isngn(:)! ISNGN array - indices of subdomain nodes in global numbering
         integer ::             lisvgvn   ! length of array ISVGVN
         integer,allocatable ::  isvgvn(:)! ISVGVN array - indices of subdomain variables in global variable numbering
         integer ::             lisegn   ! length of array ISEGN
         integer,allocatable ::  isegn(:)! ISEGN array - indices of subdomain elements in global numbering
         integer ::             lxyz1    ! length of array of coordinates
         integer ::             lxyz2    ! number of x,y,z vectors
         real(kr),allocatable :: xyz(:,:)! array of coordinates
      
         ! description of subdomain interface
         logical ::             is_interface_loaded = .false.
         integer ::             nnodi      ! number of nodes on interface
         integer ::             ndofi      ! number of dof on Iterface
         integer ::             ndofo      ! number of dof in iteriOr (ndof = ndofi + ndofo)
         integer ::             liin       ! length of IIN array 
         integer,allocatable ::  iin(:)    ! IIN array - indices of interface nodes
         integer ::             liivsvn    ! length of IIVSVN array 
         integer,allocatable ::  iivsvn(:) ! IIVSVN array - indices of Interface variables in subdomain variable numbering
         integer ::             liovsvn    ! length of IOVSVN array
         integer,allocatable ::  iovsvn(:) ! IOVSVN array - indices of interiOr variables in subdomain variable numbering

         ! boundary conditions
         logical ::             is_bc_present ! are some Dirichlet BC on subdomain?
         logical ::             is_bc_nonzero ! are some Dirichlet BC nonzero?
         logical ::             is_bc_loaded = .false. 
         integer ::             lifix         ! length of IFIX array
         integer,allocatable ::  ifix(:)      ! IFIX array - indices of fixed variables
         integer ::             lfixv         ! length of FIXV array
         real(kr),allocatable :: fixv(:)      ! FIXV array - fixed variables values
         integer ::             lbc           ! length of BC array
         real(kr),allocatable :: bc(:)        ! BC array - eliminated entries of stiffness matrix multiplied by values of fixed variables 

         ! description of corners
         logical ::             is_corners_loaded = .false.
         integer ::             ncorner                  ! number of corners on subdomain
         integer ::             lglobal_corner_number  ! length of array GLOBAL_CORNER_NUMBER
         integer, allocatable :: global_corner_number(:) ! global numbers of these corners - length NCORNER
         integer ::             licnsin                 ! length of array ICNSIN
         integer, allocatable :: icnsin(:)              ! ICNSIN array - indices of corner nodes in subdomain interface numbering

         ! description of globs
         logical ::             is_globs_loaded = .false.
         integer ::             nglob                  ! number of globs on subdomain
         integer ::             lglobal_glob_number    ! length of array GLOBAL_GLOB_NUMBER
         integer, allocatable :: global_glob_number(:) ! global numbers of these globs - lenght NGLOB
         integer ::             lnglobnodes           ! length of array NGLOBNODES
         integer, allocatable :: nglobnodes(:)        ! number of nodes in subdomain globs - lenght NGLOB
         integer ::             lignsin1              ! number of rows of IGNSIN array
         integer ::             lignsin2              ! number of cols of IGNSIN array
         integer, allocatable :: ignsin(:,:)          ! IGNSIN array - indices of glob nodes in subdomain interface numbering
         integer ::             lnglobvar              ! length of array NGLOBVAR
         integer, allocatable :: nglobvar(:)           ! number of variables in subdomain globs - lenght NGLOB
         integer ::             ligvsivn1              ! number of rows of IGVSIVN array
         integer ::             ligvsivn2              ! number of cols of IGVSIVN array
         integer, allocatable :: igvsivn(:,:)          ! IGVSIVN array - indices of glob variables in subdomain interface numbering
         integer ::             lglob_type             ! length of array GLOB_TYPE
         integer, allocatable :: glob_type(:)          ! type of globs ( 1 - face, 2 - edge)

         ! description of neighbouring of subdomain for data interchange
         logical ::             is_adj_loaded = .false. ! is adjacency of subdomains known?
         integer ::             nadj                   ! number of adjacent subdomains
         integer ::             liadj                  ! length of array IADJ
         integer, allocatable :: iadj(:)               ! indices of adjacent subdomains
         logical ::              is_neighbouring_ready =.false. ! are neighbouring arrays ready?
         integer ::             lnshnadj               ! length of array NSHNADJ
         integer, allocatable :: nshnadj(:)            ! number of nodes shared with adjacent subdomains
         integer ::             lkshvadj               ! length of array KSHVADJ
         integer, allocatable :: kshvadj(:)            ! pointers to COMMVEC array 
         integer ::             lishnadj               ! length of array ISHNADJ
         integer, allocatable :: ishnadj(:)            ! indices of nodes shared with adjacent subdomains
         integer ::              lcommvec
         real(kr),allocatable ::  commvec_out(:)       ! communication vector for sending data
         real(kr),allocatable ::  commvec_in(:)        ! communication vector for receiving data

         ! weights on interface
         integer ::              lwi
         real(kr),allocatable ::  wi(:)                ! weights at interface
         logical ::               is_weights_ready = .false. ! are weights ready?

         ! common description of joint coarse nodes/dofs
         logical ::                     is_cnodes_loaded   = .false. ! are coarse nodes activated?
         logical ::                     is_cnodes_embedded = .false. ! are coarse nodes embedded in global coarse problem?
         integer ::                     ncnodes ! number of coarse nodes (including corners and globs)
         type(cnode_type),allocatable :: cnodes(:) ! structure with all information about coarse nodes

      
         ! subdomain matrix 
         logical :: is_matrix_loaded = .false.
         ! Type of the matrix
         ! 0 - unsymmetric
         ! 1 - symmetric positive definite
         ! 2 - general symmetric
         integer :: matrixtype 
         !   type of storage
         !   0 - nothing stored
         !   1 - sparse single block full
         !   2 - sparse single block symmetric - only upper triangle
         !   3 - sparse blocks full
         !   4 - sparse blocks symmetric - only upper triangles
         !       | A11  A12 | interior variables
         !       | A21  A22 | interface variables
         integer :: istorage 
         logical :: is_assembled ! is the matrix assembled

         ! Matrix in IJA sparse format
         logical :: is_triplet = .false.
         integer :: nnza
         integer :: la
         integer,allocatable  :: i_a_sparse(:)
         integer,allocatable  :: j_a_sparse(:)
         real(kr),allocatable ::   a_sparse(:)

         ! Matrix blocks in IJA sparse format
         logical :: is_blocked = .false.
         integer :: nnza11
         integer :: la11
         integer,allocatable  :: i_a11_sparse(:)
         integer,allocatable  :: j_a11_sparse(:)
         real(kr),allocatable ::   a11_sparse(:)
         integer :: nnza12
         integer :: la12
         integer,allocatable  :: i_a12_sparse(:)
         integer,allocatable  :: j_a12_sparse(:)
         real(kr),allocatable ::   a12_sparse(:)
         integer :: nnza21
         integer :: la21
         integer,allocatable  :: i_a21_sparse(:)
         integer,allocatable  :: j_a21_sparse(:)
         real(kr),allocatable ::   a21_sparse(:)
         integer :: nnza22
         integer :: la22
         integer,allocatable  :: i_a22_sparse(:)
         integer,allocatable  :: j_a22_sparse(:)
         real(kr),allocatable ::   a22_sparse(:)
         
         logical :: is_mumps_interior_active = .false.
         logical :: is_interior_factorized   = .false.
         type(DMUMPS_STRUC) :: mumps_interior_block

         ! Matrices for BDDC 
         integer ::              ndofc            ! number of coarse degrees of freedom on subdomain
         !  matrix C with constraints on subdomain
         logical ::              is_c_loaded = .false.
         integer ::              nconstr ! number of constraints, i.e. no. of rows in C
         integer ::              nnzc
         integer ::              lc
         integer,allocatable  :: i_c_sparse(:)
         integer,allocatable  :: j_c_sparse(:)
         real(kr),allocatable ::   c_sparse(:)
         integer ::              lindrowc ! indices of rows in C matrix in global coarse dof
         integer,allocatable  ::  indrowc(:)

         !  matrix Aaug - augmented subdomain matrix
         integer ::              nnzaaug
         integer ::              laaug
         integer,allocatable  :: i_aaug_sparse(:)
         integer,allocatable  :: j_aaug_sparse(:)
         real(kr),allocatable ::   aaug_sparse(:)

         logical :: is_mumps_aug_active = .false.
         logical :: is_aug_factorized   = .false.
         type(DMUMPS_STRUC) :: mumps_aug

         ! coarse space basis functions on whole subdomain PHIS
         logical :: is_phis_prepared  = .false.
         integer ::               lphis1
         integer ::               lphis2
         real(kr),allocatable ::   phis(:,:)

         ! coarse space basis functions restricted to interface PHISI
         logical :: is_phisi_prepared  = .false.
         integer ::               lphisi1
         integer ::               lphisi2
         real(kr),allocatable ::   phisi(:,:)

         ! subdomain coarse matrix 
         logical :: is_coarse_prepared = .false.
         integer ::               coarse_matrixtype
         integer ::               lcoarsem
         real(kr),allocatable ::   coarsem(:)
         ! embedded into global coarse matrix by array INDROWC

         ! arrays connected to iterative methods
         logical :: is_rhs_loaded = .false.
         integer ::             lrhs         ! length of RHS array
         real(kr),allocatable :: rhs(:)      ! RHS array - right hand side restricted to subdomain
         logical :: is_reduced_rhs_loaded = .false.
         integer ::             lg           ! length of G array
         real(kr),allocatable :: g(:)        ! condensed right hand side on subdomain 

         integer ::             lsol         
         real(kr),allocatable :: sol(:)      ! array of solution restricted to subdomain
         integer ::             lsoli
         real(kr),allocatable :: soli(:)     ! array of solution at interface

      end type subdomain_type

contains

!***************************************
subroutine dd_init(sub,indsub,nsub,comm)
!***************************************
! Subroutine for initialization of subdomain data - assign a processor
      implicit none
! Subdomain structure
      type(subdomain_type),intent(inout) :: sub
! On which level is the subdomain
      integer,intent(in) :: indsub
! total number of subdomains
      integer,intent(in) :: nsub
! MPI communicator
      integer,intent(in) :: comm

! local vars
      integer :: myid, ierr

! orient in communicator
      call MPI_COMM_RANK(comm,myid,ierr)

! initialize basic structure
      sub%isub              = indsub
      sub%nsub              = nsub
      sub%proc              = myid
      sub%comm              = comm
      sub%is_initialized    = .true.

end subroutine

!*************************************************
subroutine dd_read_mesh_from_file(sub,problemname)
!*************************************************
! Subroutine for reading of subdomain data from separate SMD file
      use module_utils
      implicit none

! Subdomain structure
      type(subdomain_type),intent(inout) :: sub
! name of problem
      character(*),intent(in) :: problemname

! local variables
      integer :: isub
      integer :: idsmd, iglob, ivar, ignode
      character(lfnamex):: fname

      integer :: nelem, nnod, ndof, ndim, ncorner, nglob, nadj
      integer ::             lnndf,   lnnet,   linet 
      integer, allocatable :: nndf(:), nnet(:), inet(:)
      integer ::             lisngn
      integer, allocatable :: isngn(:)
      integer ::             lisvgvn
      integer, allocatable :: isvgvn(:)
      integer ::             lisegn
      integer, allocatable :: isegn(:)
      integer ::              lxyz1, lxyz2
      real(kr), allocatable :: xyz(:,:)
      integer :: nnodi, ndofi, ndofo
      integer ::             liin,   liivsvn,   liovsvn
      integer, allocatable :: iin(:), iivsvn(:), iovsvn(:)
      integer ::             lifix
      integer, allocatable :: ifix(:)
      integer ::              lfixv
      real(kr), allocatable :: fixv(:)
      integer ::             liadj
      integer, allocatable :: iadj(:)
      integer ::              lglobal_corner_number
      integer, allocatable ::  global_corner_number(:)
      integer ::              licnsin 
      integer, allocatable ::  icnsin(:)
      integer ::              lglobal_glob_number
      integer, allocatable ::  global_glob_number(:)
      integer ::              lglob_type
      integer, allocatable ::  glob_type(:)
      integer ::              lnglobnodes 
      integer, allocatable ::  nglobnodes(:)
      integer ::              lignsin1, lignsin2
      integer, allocatable ::  ignsin(:,:)
      integer ::              lnglobvar 
      integer, allocatable ::  nglobvar(:)
      integer ::              ligvsivn1, ligvsivn2
      integer, allocatable ::  igvsivn(:,:)
      integer ::              lrhss
      real(kr), allocatable :: rhss(:)

      if (.not.sub%is_initialized) then
         call error('DD_READ_MESH_FROM_FILE','subdomain not initialized')
      end if

      isub = sub%isub

      ! open subdomain SMD file with mesh description
      call getfname(problemname,isub,'SMD',fname)
      if (debug) then
         write(*,*) 'DD_READ_MESH_FROM_FILE: Opening file fname: ', trim(fname)
      end if
      call allocate_unit(idsmd)
      open (unit=idsmd,file=trim(fname),status='old',form='formatted')

      ! read data from file
      ! ---header
      read(idsmd,*) nelem, nnod, ndof, ndim

      ! ---NNDF array
      ! ---NNET array
      lnndf = nnod
      lnnet = nelem
      allocate (nndf(lnndf),nnet(lnnet))
      read(idsmd,*) nndf
      read(idsmd,*) nnet

      ! ---INET array
      linet = sum(nnet)
      allocate (inet(linet))
      read(idsmd,*) inet

      ! ---ISNGN array
      lisngn = nnod
      allocate (isngn(lisngn))
      read(idsmd,*) isngn

      ! ---ISVNGVN array
      lisvgvn = ndof
      allocate (isvgvn(lisvgvn))
      call zero(isvgvn,lisvgvn)

      ! ---ISEGN array
      lisegn = nelem
      allocate (isegn(lisegn))
      read(idsmd,*) isegn

      ! --- coordinates
      lxyz1 = nnod
      lxyz2 = ndim
      allocate(xyz(lxyz1,lxyz2))
      read(idsmd,*) xyz

      ! ---interface data
      read(idsmd,*) nnodi
      read(idsmd,*) ndofi
      read(idsmd,*) ndofo

      liin    = nnodi
      liivsvn = ndofi
      liovsvn = ndofo
      allocate (iin(liin), iivsvn(liivsvn), iovsvn(liovsvn))
      read(idsmd,*) iin
      read(idsmd,*) iivsvn
      read(idsmd,*) iovsvn

      ! --- boundary conditions
      lifix = ndof
      lfixv = ndof
      allocate(ifix(lifix),fixv(lfixv))
      read(idsmd,*) ifix
      read(idsmd,*) fixv

      ! --- neighbouring
      read(idsmd,*) nadj
      liadj = nadj
      allocate(iadj(liadj))
      read(idsmd,*) iadj

      ! --- corners 
      read(idsmd,*) ncorner
      lglobal_corner_number = ncorner
      licnsin               = ncorner
      allocate(global_corner_number(lglobal_corner_number),icnsin(licnsin))
      read(idsmd,*) global_corner_number
      read(idsmd,*) icnsin

      ! --- globs
      read(idsmd,*) nglob
      lglobal_glob_number = nglob
      lnglobvar           = nglob
      lnglobnodes         = nglob
      allocate(global_glob_number(lglobal_glob_number),nglobvar(lnglobvar),nglobnodes(lnglobnodes))
      read(idsmd,*) global_glob_number

      read(idsmd,*) nglobnodes
      lignsin1 = nglob
      lignsin2 = maxval(nglobnodes)
      allocate(ignsin(lignsin1,lignsin2))
      do iglob = 1,nglob
         ! read glob variables
         read(idsmd,*) (ignsin(iglob,ignode),ignode = 1,nglobnodes(iglob))
         ! pad the glob with zeros
         do ignode = nglobnodes(iglob) + 1,lignsin2
            ignsin(iglob,ignode) = 0
         end do
      end do
      read(idsmd,*) nglobvar
      ligvsivn1 = nglob
      ligvsivn2 = maxval(nglobvar)
      allocate(igvsivn(ligvsivn1,ligvsivn2))
      do iglob = 1,nglob
         ! read glob variables
         read(idsmd,*) (igvsivn(iglob,ivar),ivar = 1,nglobvar(iglob))
         ! pad the glob with zeros
         do ivar = nglobvar(iglob) + 1,ligvsivn2
            igvsivn(iglob,ivar) = 0
         end do
      end do
      lglob_type = nglob
      allocate(glob_type(lglob_type))
      read(idsmd,*) glob_type

      ! --- right hand side
      lrhss = ndof
      allocate(rhss(lrhss))
      read(idsmd,*) rhss

      close(idsmd)
      if (debug) then
         write(*,*) 'DD_READ_MESH_FROM_FILE: Data read successfully.'
      end if

      ! load data to structure
      call dd_upload_sub_mesh(sub, nelem, nnod, ndof, ndim, &
                              nndf,lnndf, nnet,lnnet, inet,linet, &
                              isngn,lisngn, isvgvn,lisvgvn, isegn,lisegn,&
                              xyz,lxyz1,lxyz2)
      call dd_upload_sub_adj(sub, nadj, iadj,liadj)
      !call dd_upload_sub_interface(sub, nnodi, ndofi, ndofo, &
      !                             iin,liin, iivsvn,liivsvn, iovsvn,liovsvn)
      call dd_upload_sub_corners(sub, ncorner,&
                                 global_corner_number,lglobal_corner_number, icnsin,licnsin)
      call dd_upload_sub_globs(sub, nglob, &
                               global_glob_number,lglobal_glob_number, &
                               nglobnodes,lnglobnodes, nglobvar,lnglobvar,&
                               ignsin,lignsin1,lignsin2, igvsivn,ligvsivn1,ligvsivn2,&
                               glob_type,lglob_type)
      call dd_upload_bc(sub, ifix,lifix, fixv,lfixv)
      call dd_upload_rhs(sub, rhss,lrhss)

      if (debug) then
         write(*,*) 'DD_READ_MESH_FROM_FILE: Data loaded successfully.'
      end if

      deallocate (nndf,nnet)
      deallocate (inet)
      deallocate (isngn)
      deallocate (isegn)
      deallocate (xyz)
      deallocate (iin, iivsvn, iovsvn)
      deallocate (iadj)
      deallocate (ifix,fixv)
      deallocate (glob_type)
      deallocate (global_glob_number,nglobvar,nglobnodes)
      deallocate (ignsin)
      deallocate (igvsivn)
      deallocate (global_corner_number,icnsin)
      deallocate (rhss)

end subroutine

!**************************************************************
subroutine dd_read_matrix_from_file(sub,matrixtype,problemname)
!**************************************************************
! Subroutine for initialization of subdomain data
      use module_utils
      use module_sm
      implicit none
      include "mpif.h"

! Subdomain structure
      type(subdomain_type),intent(inout) :: sub
! type of matrix
      integer,intent(in) :: matrixtype
! name of problem
      character(*),intent(in) :: problemname

! local variables
      character(lfnamex):: fname

      integer :: isub
      integer :: idelm, nelem, nnod, ndof

!     Matrix in IJA sparse format - triplet
      integer::  la, nnza
      integer,allocatable  :: i_sparse(:), j_sparse(:)
      real(kr),allocatable :: a_sparse(:)
      logical :: is_assembled


!     kdof array
      integer::              lkdof
      integer,allocatable  :: kdof(:)

!     BC array
      integer::              lbc
      real(kr),allocatable :: bc(:)

      integer :: inod
      integer :: ios

! check prerequisites
      if (.not.sub%is_initialized) then
         call error('DD_READ_MATRIX_FROM_FILE','Subdomain is not initialized.')
      end if

! get data from subdomain
      isub  = sub%isub
      nelem = sub%nelem
      nnod  = sub%nnod

! check prerequisites
      if (.not.sub%is_mesh_loaded) then
         write(*,*) 'DD_READ_MATRIX_FROM_FILE: Mesh is not loaded for subdomain:', isub
         call error_exit
      end if
      if (.not.sub%is_bc_loaded) then
         write(*,*) 'DD_READ_MATRIX_FROM_FILE: BC not loaded for subdomain:', isub
         call error_exit
      end if

      ! import PMD file to memory
      call sm_pmd_get_length(matrixtype,nelem,sub%inet,sub%linet,&
                             sub%nnet,sub%lnnet,sub%nndf,sub%lnndf,&
                             la)
      ! allocate memory for triplet
      allocate(i_sparse(la), j_sparse(la), a_sparse(la))

      ! Creation of field KDOF(NNOD) with addresses before first global dof of node
      lkdof = nnod
      allocate(kdof(lkdof))
      kdof(1) = 0
      do inod = 2,nnod
         kdof(inod) = kdof(inod-1) + sub%nndf(inod-1)
      end do

      ! try with subdomain ELM file with element matrices
      call allocate_unit(idelm)
      call getfname(problemname,isub,'ELM',fname)
      open (unit=idelm,file=trim(fname),status='old',form='unformatted',iostat=ios)
      if (ios.eq.0) then
         ! the subdomain file should exist, try reading it
         call sm_pmd_load(idelm,nelem,sub%inet,sub%linet,&
                          sub%nnet,sub%lnnet,sub%nndf,sub%lnndf,&
                          kdof,lkdof,&
                          i_sparse, j_sparse, a_sparse, la)
         close (idelm)
      else
         write(*,*) 'WARNING: File ',trim(fname), ' does not open. Trying with global ELM file...'
         ! (this could have a negative impact on performance)

         ! use global filename
         fname = trim(problemname)//'.ELM'
         open (unit=idelm,file=trim(fname),status='old',form='unformatted')
         call sm_pmd_load_masked(idelm,nelem,sub%inet,sub%linet,&
                                 sub%nnet,sub%lnnet,sub%nndf,sub%lnndf,&
                                 kdof,lkdof,sub%isegn,sub%lisegn,&
                                 i_sparse, j_sparse, a_sparse, la)
         close (idelm)
      end if
      if (debug) then
         write(*,*) 'DD_READ_MATRIX_FROM_FILE: File with element matrices read.'
      end if

      ! eliminate boundary conditions
      if (sub%is_bc_present) then
         ndof =  sub%ndof
         if (sub%is_bc_nonzero) then
            lbc = ndof
            allocate(bc(lbc))
         else
            lbc = 0

         end if

         ! eliminate natural BC
         call sm_apply_bc(sub%ifix,sub%lifix,sub%fixv,sub%lfixv,&
                          i_sparse,j_sparse,a_sparse,la, bc,lbc)
         if (sub%is_bc_nonzero) then
            call dd_load_eliminated_bc(sub, bc,lbc)
         end if
      end if

      ! load matrix to our structure
      nnza = la
      is_assembled = .false.
      call dd_load_matrix_triplet(sub, matrixtype, &
                                  i_sparse,j_sparse,a_sparse,la,nnza,is_assembled)


      if (debug) then
         write(*,*) 'DD_READ_MATRIX_FROM_FILE: isub =',isub,' matrix loaded.'
      end if

      deallocate(kdof)
      deallocate(i_sparse, j_sparse, a_sparse)
      if (allocated(bc)) then
         deallocate(bc)
      end if

end subroutine

!****************************************************************************************
subroutine dd_read_matrix_by_root(suba,lsuba, comm_all,problemname,nsub,nelem,matrixtype,&
                                  sub2proc,lsub2proc,indexsub,lindexsub,iets,liets)
!****************************************************************************************
! Subroutine for reading matrix from global *.ELM file and distribution of
! elements to subdomains
      use module_pp, only: pp_get_nevax
      use module_sm
      use module_utils
      implicit none
      include "mpif.h"

! array of sub structure
      integer,intent(in) ::                lsuba
      type(subdomain_type),intent(inout) :: suba(lsuba)
! communicator
      integer,intent(in) :: comm_all
! name of problem
      character(*),intent(in) :: problemname
! number of subdomains
      integer,intent(in) :: nsub
! number of elements
      integer,intent(in) :: nelem
! type of matrix
      integer,intent(in) :: matrixtype
! division of subdomains to processors
      integer,intent(in) :: lsub2proc
      integer,intent(in) ::  sub2proc(lsub2proc)
! global indices of local subdomains
      integer,intent(in) :: lindexsub
      integer,intent(in) ::  indexsub(lindexsub)
! indices of subdomains for elements
      integer,intent(in) :: liets
      integer,intent(in) ::  iets(liets)

! local variables
      integer :: myid, nproc, ierr
      integer :: stat(MPI_STATUS_SIZE)
      integer :: idelm, nelems, nnods, ndofs
      integer :: nevax_sub, nevax_loc, nevax, lelmx, lelm
      real(kr),allocatable :: elm(:)

!     Matrix in IJA sparse format - triplet
      integer::  la

!     kdof array
      integer::              lkdof
      integer,allocatable  :: kdof(:)

!     subp array of subdomains for processors
      integer::              lsubp
      integer,allocatable  :: subp(:)

!     locsubnumber array of local numbers subdomains for processors
      integer::              llocsubnumber
      integer,allocatable  :: locsubnumber(:)

!     BC array
      integer::              lbc
      real(kr),allocatable :: bc(:)

      integer :: i, ie, iproc, inods, isub, isub_loc, ios


! check dimension
      if (nelem.ne.liets) then
         call error('DD_READ_MATRIX_BY_ROOT','wrong length of array IETS')
      end if

! orient in communicator
      call MPI_COMM_RANK(comm_all,myid,ierr)
      call MPI_COMM_SIZE(comm_all,nproc,ierr)

      if (debug) then
         if (nproc+1.ne.lsub2proc) then
            call error('DD_READ_MATRIX_BY_ROOT','wrong length of array sub2proc')
         end if
         if (any(suba%comm .ne. comm_all)) then
            call error('DD_READ_MATRIX_BY_ROOT','communicator mismatch')
         end if
         if (any(suba%proc .ne. myid)) then
            call error('DD_READ_MATRIX_BY_ROOT','communicator mismatch')
         end if
         if (any(suba%nsub .ne. nsub)) then
            call error('DD_READ_MATRIX_BY_ROOT','number of subdomains mismatch')
         end if
      end if

      ! prepare arrays for quick mappings of processors for subdomains
      lsubp = nsub
      allocate(subp(lsubp))
      do iproc = 0,nproc-1
         do isub = sub2proc(iproc+1),sub2proc(iproc+2)-1
            subp(isub) = iproc
         end do
      end do
      llocsubnumber = nsub
      allocate(locsubnumber(llocsubnumber))
      locsubnumber = -1
      do isub_loc = 1,lindexsub
         isub = indexsub(isub_loc)

         locsubnumber(isub) = isub_loc
      end do

      !print *,'subp:',subp
      !call flush(6)
      !call MPI_BARRIER(comm_all,ierr)

      ! prepare memory for sparse matrix at each subdomain
      nevax_loc = 0
      do isub_loc = 1,lindexsub
         isub = indexsub(isub_loc)

         if (subp(isub).ne.myid) then
            cycle
         end if

         ! get data from subdomain
         nelems = suba(isub_loc)%nelem
         nnods  = suba(isub_loc)%nnod

         ! measure sparse matrix in unassembled format
         call sm_pmd_get_length(matrixtype,nelems,suba(isub_loc)%inet,suba(isub_loc)%linet,&
                                suba(isub_loc)%nnet,suba(isub_loc)%lnnet,suba(isub_loc)%nndf,suba(isub_loc)%lnndf,&
                                la)
         ! allocate memory for triplet
         suba(isub_loc)%matrixtype = matrixtype
         if (matrixtype.eq.0) then
            suba(isub_loc)%istorage   = 1
         else if (matrixtype.eq.1 .or. matrixtype.eq.2) then
            suba(isub_loc)%istorage   = 2
         else
            call error('DD_READ_MATRIX_BY_ROOT','strange type of matrix')
         end if
         suba(isub_loc)%la       = la
         allocate(suba(isub_loc)%i_a_sparse(la))
         allocate(suba(isub_loc)%j_a_sparse(la))
         !print *,'myid =',myid,'isub',isub,'la =',la
         !print *,'myid =',myid,'isub',isub,'isegn:',sub(isub_loc)%isegn
         !print *,'myid =',myid,'isub',isub,'inet:',sub(isub_loc)%inet
         !print *,'myid =',myid,'isub',isub,'nnet:',sub(isub_loc)%nnet
         !print *,'myid =',myid,'isub',isub,'nndf:',sub(isub_loc)%nndf
         !call flush(6)

         ! Creation of field KDOF(NNOD) with addresses before first global dof of node
         lkdof = nnods
         allocate(kdof(lkdof))
         if (lkdof.gt.0) then
            kdof(1) = 0
            do inods = 2,nnods
               kdof(inods) = kdof(inods-1) + suba(isub_loc)%nndf(inods-1)
            end do
         end if
         ! Prepare numbering of element matrices
         call sm_pmd_make_element_numbering(nelems,suba(isub_loc)%inet,suba(isub_loc)%linet,&
                                            suba(isub_loc)%nnet,suba(isub_loc)%lnnet,&
                                            suba(isub_loc)%nndf,suba(isub_loc)%lnndf,&
                                            kdof,lkdof,&
                                            suba(isub_loc)%i_a_sparse, suba(isub_loc)%j_a_sparse, la)
         deallocate(kdof)

         allocate(suba(isub_loc)%a_sparse(la))
         suba(isub_loc)%nnza = 0

         ! find nevax_sub
         call pp_get_nevax(nelems,&
                           suba(isub_loc)%inet,suba(isub_loc)%linet,&
                           suba(isub_loc)%nnet,suba(isub_loc)%lnnet,&
                           suba(isub_loc)%nndf,suba(isub_loc)%lnndf,nevax_sub)
         if (nevax_sub.gt.nevax_loc) then
            nevax_loc = nevax_sub
         end if
      end do
      ! get global maximum of nevax
!*****************************************************************MPI
      call MPI_ALLREDUCE(nevax_loc,nevax,1, MPI_INTEGER, MPI_MAX, comm_all, ierr) 
!*****************************************************************MPI

      ! prepare memory for one element matrix
      ! determine length by nevax
      if (matrixtype.eq.1 .or. matrixtype.eq.2) then
         lelmx = (nevax+1)*nevax / 2
      else
         lelmx = nevax*nevax
      end if
      allocate(elm(lelmx))

      ! open ELM file on root processor
      if (myid.eq.0) then
         ! ELM - element stiffness matrices - structure:
         call allocate_unit(idelm)
         open(unit=idelm,file=trim(problemname)//'.ELM',status='old',form='unformatted',iostat=ios)
         if (ios.ne.0) then
            write(*,*) 'DD_READ_MATRIX_BY_ROOT','Problem opening file ',trim(problemname)//'.ELM'
            call error_exit
         end if
         rewind idelm
      end if

      ! loop over elements
      do ie = 1,nelem
         isub = iets(ie)

         ! get processor taking care of this subdomain
         iproc = subp(isub)

         if (myid.eq.0) then
            ! root reads the matrix from file
            read(idelm) lelm,(elm(i),i = 1,lelm)
            if (iproc.eq.0) then
               ! if it is to be stored by 0, do not send messages, just store it
               isub_loc = locsubnumber(isub)
               do i = 1,lelm
                  suba(isub_loc)%a_sparse(suba(isub_loc)%nnza + i) = elm(i)
               end do
               suba(isub_loc)%nnza = suba(isub_loc)%nnza + lelm
            else
               ! send messages
               call MPI_SEND(lelm,1,MPI_INTEGER,iproc,isub,comm_all,ierr)
               call MPI_SEND(elm,lelm,MPI_DOUBLE_PRECISION,iproc,isub,comm_all,ierr)
            end if
         else 
            if (myid.eq.iproc) then
               call MPI_RECV(lelm,1,MPI_INTEGER,0,isub,comm_all,stat,ierr)
               call MPI_RECV(elm,lelm,MPI_DOUBLE_PRECISION,0,isub,comm_all,stat,ierr)

               isub_loc = locsubnumber(isub)
               do i = 1,lelm
                  suba(isub_loc)%a_sparse(suba(isub_loc)%nnza + i) = elm(i)
               end do
               suba(isub_loc)%nnza = suba(isub_loc)%nnza + lelm
            end if
         end if
      end do
      
      ! close *.ELM file
      if (myid.eq.0) then
         close(idelm)
      end if

      ! free memory
      deallocate(elm)

      ! finalize input of sparse matrix at each subdomain
      do isub_loc = 1,lindexsub
         isub = indexsub(isub_loc)

         if (subp(isub).ne.myid) then
            cycle
         end if

         la = suba(isub_loc)%la

         suba(isub_loc)%is_assembled = .false.

         suba(isub_loc)%is_matrix_loaded = .true.
         suba(isub_loc)%is_triplet       = .true.

      ! eliminate boundary conditions
         if (suba(isub_loc)%is_bc_present) then
            ndofs =  suba(isub_loc)%ndof
            if (suba(isub_loc)%is_bc_nonzero) then
               lbc = ndofs
               allocate(bc(lbc))
            else
               lbc = 0
            end if

            ! eliminate natural BC
            call sm_apply_bc(suba(isub_loc)%ifix,suba(isub_loc)%lifix,suba(isub_loc)%fixv,suba(isub_loc)%lfixv,&
                             suba(isub_loc)%i_a_sparse,suba(isub_loc)%j_a_sparse,suba(isub_loc)%a_sparse,la, bc,lbc)
            if (suba(isub_loc)%is_bc_nonzero) then
               call dd_load_eliminated_bc(suba(isub_loc), bc,lbc)
            end if
         end if

         if (allocated(bc)) then
            deallocate(bc)
         end if
         if (debug) then
            write(*,*) 'DD_READ_MATRIX_BY_ROOT: myid = ',myid, 'isub =',isub,' matrix loaded.'
         end if
      end do

      deallocate(locsubnumber)
      deallocate(subp)

end subroutine

!***********************************************************************************
subroutine dd_gather_matrix(suba,lsuba, elma,lelma, comm_all,nsub,nelem,matrixtype,&
                            sub2proc,lsub2proc,indexsub,lindexsub,&
                            elm2proc,lelm2proc,indexelm,lindexelm,&
                            iets,liets)
!***********************************************************************************
! Subroutine for reading matrix from global *.ELM file and distribution of
! elements to subdomains
      use module_pp, only: pp_get_nevax
      use module_sm
      use module_utils
      implicit none
      include "mpif.h"

! array of sub structure for actual subdomains
      integer,intent(in) ::                lsuba
      type(subdomain_type),intent(inout) :: suba(lsuba)
! array of sub structure for actual elements (subdomains on previous level)
      integer,intent(in) ::                lelma
      type(subdomain_type),intent(in) ::    elma(lelma)
! communicator including all processes involved in subdomain
      integer,intent(in) :: comm_all
! number of subdomains 
      integer,intent(in) :: nsub
! number of elements (subdomains on previous level)
      integer,intent(in) :: nelem
! type of matrix
      integer,intent(in) :: matrixtype
! division of subdomains to processors at previous level
      integer,intent(in) :: lsub2proc
      integer,intent(in) ::  sub2proc(lelm2proc)
! global indices of local subdomains
      integer,intent(in) :: lindexsub
      integer,intent(in) ::  indexsub(lindexsub)
! division of elements to processors 
      integer,intent(in) :: lelm2proc
      integer,intent(in) ::  elm2proc(lelm2proc)
! global indices of local elements
      integer,intent(in) :: lindexelm
      integer,intent(in) ::  indexelm(lindexelm)
! indices of subdomains for elements
      integer,intent(in) :: liets
      integer,intent(in) ::  iets(liets)

! local variables
      integer :: myid, nproc, ierr
      integer :: stat(MPI_STATUS_SIZE)
      integer :: nelems, nnods, ndofs
      integer :: nevax_sub, nevax_loc, nevax, lelmx, lelm
      real(kr),allocatable :: elm(:)

!     Matrix in IJA sparse format - triplet
      integer::  la

!     kdof array
      integer::              lkdof
      integer,allocatable  :: kdof(:)

!     subp array of subdomains for processors
      integer::              lsubp
      integer,allocatable  :: subp(:)

!     locsubnumber array of local numbers subdomains for processors
      integer::              llocsubnumber
      integer,allocatable  :: locsubnumber(:)

!     elmp array of elements for processors
      integer::              lelmp
      integer,allocatable  :: elmp(:)

!     locsubnumber array of local numbers subdomains for processors
      integer::              llocelmnumber
      integer,allocatable  :: locelmnumber(:)

!     BC array
      integer::              lbc
      real(kr),allocatable :: bc(:)

      integer :: i, ie, iproc, inods, isub, isub_loc, ie_loc, indel_loc,&
                 iprocelm, iprocsub


! check dimension
      if (nelem.ne.liets) then
         call error('DD_GATHER_MATRIX','wrong length of array IETS')
      end if

! orient in communicators
      call MPI_COMM_RANK(comm_all,myid,ierr)
      call MPI_COMM_SIZE(comm_all,nproc,ierr)

      if (nproc+1.ne.lsub2proc) then
         call error('DD_GATHER_MATRIX','wrong length of array sub2proc')
      end if
      if (nproc+1.ne.lelm2proc) then
         call error('DD_GATHER_MATRIX','wrong length of array elm2proc')
      end if

      ! prepare arrays for quick mappings of processors for subdomains
      lsubp = nsub
      allocate(subp(lsubp))
      do iproc = 0,nproc-1
         do isub = sub2proc(iproc+1),sub2proc(iproc+2)-1
            subp(isub) = iproc
         end do
      end do
      llocsubnumber = nsub
      allocate(locsubnumber(llocsubnumber))
      locsubnumber = -1
      do isub_loc = 1,lindexsub
         isub = indexsub(isub_loc)

         locsubnumber(isub) = isub_loc
      end do

      lelmp = nelem
      allocate(elmp(lelmp))
      do iproc = 0,nproc-1
         do ie = elm2proc(iproc+1),elm2proc(iproc+2)-1
            elmp(ie) = iproc
         end do
      end do
      llocelmnumber = nelem
      allocate(locelmnumber(llocelmnumber))
      locelmnumber = -1
      do ie_loc = 1,lindexelm
         ie = indexelm(ie_loc)

         locelmnumber(ie) = ie_loc
      end do

      ! prepare memory for sparse matrix at each subdomain
      nevax_loc = 0
      do isub_loc = 1,lindexsub
         isub = indexsub(isub_loc)

         if (subp(isub).ne.myid) then
            cycle
         end if

         ! get data from subdomain
         nelems = suba(isub_loc)%nelem
         nnods  = suba(isub_loc)%nnod

         ! measure sparse matrix in unassembled format
         call sm_pmd_get_length(matrixtype,nelems,suba(isub_loc)%inet,suba(isub_loc)%linet,&
                                suba(isub_loc)%nnet,suba(isub_loc)%lnnet,suba(isub_loc)%nndf,suba(isub_loc)%lnndf,&
                                la)
         ! allocate memory for triplet
         suba(isub_loc)%matrixtype = matrixtype
         if (matrixtype.eq.0) then
            suba(isub_loc)%istorage   = 1
         else if (matrixtype.eq.1 .or. matrixtype.eq.2) then
            suba(isub_loc)%istorage   = 2
         else
            call error('DD_READ_MATRIX_BY_ROOT','strange type of matrix')
         end if
         suba(isub_loc)%la       = la
         allocate(suba(isub_loc)%i_a_sparse(la))
         allocate(suba(isub_loc)%j_a_sparse(la))

         ! Creation of field KDOF(NNOD) with addresses before first global dof of node
         lkdof = nnods
         allocate(kdof(lkdof))
         kdof(1) = 0
         do inods = 2,nnods
            kdof(inods) = kdof(inods-1) + suba(isub_loc)%nndf(inods-1)
         end do
         ! Prepare numbering of element matrices
         call sm_pmd_make_element_numbering(nelems,suba(isub_loc)%inet,suba(isub_loc)%linet,&
                                            suba(isub_loc)%nnet,suba(isub_loc)%lnnet,&
                                            suba(isub_loc)%nndf,suba(isub_loc)%lnndf,&
                                            kdof,lkdof,&
                                            suba(isub_loc)%i_a_sparse, suba(isub_loc)%j_a_sparse, la)
         deallocate(kdof)

         allocate(suba(isub_loc)%a_sparse(la))
         suba(isub_loc)%nnza = 0

         ! find nevax_sub
         call pp_get_nevax(nelems,suba(isub_loc)%inet,suba(isub_loc)%linet,&
                           suba(isub_loc)%nnet,suba(isub_loc)%lnnet,&
                           suba(isub_loc)%nndf,suba(isub_loc)%lnndf,nevax_sub)
         if (nevax_sub.gt.nevax_loc) then
            nevax_loc = nevax_sub
         end if
      end do
      ! get global maximum of nevax
!*****************************************************************MPI
      call MPI_ALLREDUCE(nevax_loc,nevax,1, MPI_INTEGER, MPI_MAX, comm_all, ierr) 
!*****************************************************************MPI

      ! prepare memory for one element matrix
      ! determine length by nevax
      if (matrixtype.eq.1 .or. matrixtype.eq.2) then
         lelmx = (nevax+1)*nevax / 2
      else
         lelmx = nevax*nevax
      end if
      allocate(elm(lelmx))

      ! loop over elements
      do ie = 1,nelem
         
         isub = iets(ie)

         ! get processor taking care of this subdomain
         iprocsub = subp(isub)
         ! get processor taking care of this element
         iprocelm = elmp(ie)

         if (myid.eq.iprocelm) then
            indel_loc = locelmnumber(ie) 
            if (myid.ne.iprocelm) then
            ! send messages
               call MPI_SEND(elma(indel_loc)%lcoarsem,1,                      &
                             MPI_INTEGER,     iprocsub,ie,comm_all,ierr)
               call MPI_SEND(elma(indel_loc)%coarsem,elma(indel_loc)%lcoarsem,&
                             MPI_DOUBLE_PRECISION,iprocsub,ie,comm_all,ierr)
            else
               ! copy it to array elm
               lelm = elma(indel_loc)%lcoarsem
               do i = 1,lelm
                  elm(i) = elma(indel_loc)%coarsem(i)
               end do
            end if
         end if
         if (myid.eq.iprocsub) then
            ! storing the matrix - get it from owner by MPI
            if (myid.ne.iprocelm) then
               call MPI_RECV(lelm,1,  MPI_INTEGER,         iprocelm,ie,comm_all,stat,ierr)
               call MPI_RECV(elm,lelm,MPI_DOUBLE_PRECISION,iprocelm,ie,comm_all,stat,ierr)
            end if

            isub_loc = locsubnumber(isub)
            do i = 1,lelm
               suba(isub_loc)%a_sparse(suba(isub_loc)%nnza + i) = elm(i)
            end do
            suba(isub_loc)%nnza = suba(isub_loc)%nnza + lelm
         end if
      end do
      
      ! free memory
      deallocate(elm)

      ! finalize input of sparse matrix at each subdomain
      do isub_loc = 1,lindexsub
         isub = indexsub(isub_loc)

         if (subp(isub).ne.myid) then
            cycle
         end if

         la = suba(isub_loc)%la

         suba(isub_loc)%is_assembled = .false.
         suba(isub_loc)%is_matrix_loaded = .true.
         suba(isub_loc)%is_triplet       = .true.

      ! eliminate boundary conditions
         if (suba(isub_loc)%is_bc_present) then
            ndofs =  suba(isub_loc)%ndof
            if (suba(isub_loc)%is_bc_nonzero) then
               lbc = ndofs
               allocate(bc(lbc))
            else
               lbc = 0
            end if

            ! eliminate natural BC
            call sm_apply_bc(suba(isub_loc)%ifix,suba(isub_loc)%lifix,suba(isub_loc)%fixv,suba(isub_loc)%lfixv,&
                             suba(isub_loc)%i_a_sparse,suba(isub_loc)%j_a_sparse,suba(isub_loc)%a_sparse,la, bc,lbc)
            if (suba(isub_loc)%is_bc_nonzero) then
               call dd_load_eliminated_bc(suba(isub_loc), bc,lbc)
            end if
         end if

         if (allocated(bc)) then
            deallocate(bc)
         end if
         if (debug) then
            write(*,*) 'DD_GATHER_MATRIX: myid= ',myid, 'isub =',isub,': matrix loaded.'
         end if
      end do

      deallocate(locsubnumber)
      deallocate(subp)
      deallocate(locelmnumber)
      deallocate(elmp)

end subroutine

!********************************************************************************
subroutine dd_write_solution_to_file(problemname, isub, sol,lsol, print_solution)
!********************************************************************************
! Subroutine for writing subdomain solution into file
      use module_utils
      implicit none

      ! name of problem
      character(*),intent(in) :: problemname
      ! index of subdomain
      integer,intent(in):: isub

      ! solution
      integer,intent(in)  :: lsol
      real(kr),intent(in) ::  sol(lsol)

      ! print the solution to screen?
      logical,intent(in)  :: print_solution

      ! local vars
      integer :: idsols, i
      character(lfnamex):: fname

      ! open subdomain SOLS file for solution
      call getfname(problemname,isub,'SOLS',fname)
      if (debug) then
         write(*,*) 'DD_WRITE_SOLUTION_TO_FILE: Opening file fname: ', trim(fname)
      end if
      call allocate_unit(idsols)
      open (unit=idsols,file=trim(fname),status='replace',form='unformatted')

      rewind idsols
      write(idsols) (sol(i),i=1,lsol)
      
      if (print_solution) then
         write(*,*) 'isub =',isub,' solution: '
         write(*,'(e15.7)') sol
      end if

      close(idsols)

end subroutine

!**************************************************************************
subroutine dd_localize_mesh(sub,isub,ndim,nelem,nnod,&
                            inet,linet,nnet,lnnet,nndf,lnndf,xyz,lxyz1,lxyz2,&
                            iets,liets)
!**************************************************************************
! Subroutine for localization of global mesh to particular subdomain, 
! loads the data directly to the structure
      use module_pp, only: pp_create_submesh
      use module_utils
      implicit none

! Subdomain structure
      type(subdomain_type),intent(inout) :: sub
! number of subdomain
      integer,intent(in) :: isub
! space dimension
      integer,intent(in) :: ndim
! number of elements
      integer,intent(in) :: nelem
! number of nodes
      integer,intent(in) :: nnod
! global mesh data
      integer,intent(in) :: linet
      integer,intent(in) ::  inet(linet)
      integer,intent(in) :: lnnet
      integer,intent(in) ::  nnet(lnnet)
      integer,intent(in) :: lnndf
      integer,intent(in) ::  nndf(lnndf)
      integer,intent(in) :: lxyz1, lxyz2
      real(kr),intent(in) :: xyz(lxyz1,lxyz2)
! division into subdomains
      integer,intent(in) :: liets
      integer,intent(in) ::  iets(liets)

! local vars
      integer ::            lkdof
      integer,allocatable :: kdof(:)
      integer ::            linets,   lnnets,   lnndfs,   lkdofs,   lisegns,   lisngns,   lisvgvns
      integer,allocatable :: inets(:), nnets(:), nndfs(:), kdofs(:), isegns(:), isngns(:), isvgvns(:)
      integer ::             lxyzs1, lxyzs2
      real(kr),allocatable::  xyzs(:,:)

      integer ::             ndofs, nelems, nnods, ndofn
      integer ::             i, j, ie, inod, inods, idofn 
      integer ::             indn, indvg, indvs

      if (isub .ne. sub%isub) then
         write(*,*) 'DD_LOCALIZE_MESH: Subdomain index mismatch for subdomain:', isub
         call error_exit
      end if
! check prerequisites
      if (.not.sub%is_initialized) then
         write(*,*) 'DD_LOCALIZE_MESH: Subdomain ', isub,' not initialized.'
         call error_exit
      end if

      ! create array KDOF
      lkdof = nnod
      allocate(kdof(lkdof))
      if (lkdof.gt.0) then
         kdof(1) = 0
         do inod = 2,nnod
            kdof(inod) = kdof(inod-1) + nndf(inod-1)
         end do
      end if

      nelems = 0
      linets = 0
      do ie = 1,nelem
         if (iets(ie).eq.isub) then
            nelems = nelems + 1
            linets = linets + nnet(ie)
         end if
      end do
      lnnets  = nelems
      lisegns = nelems
      lisngns = linets
      allocate(inets(linets),nnets(lnnets),isegns(lisegns),isngns(lisngns))

      call pp_create_submesh(isub,nelem,inet,linet,&
                             nnet,lnnet,&
                             iets,liets,&
                             nnods,inets,linets,nnets,lnnets,&
                             isegns,lisegns,isngns,lisngns)
      ! correct lisngns
      lisngns = nnods
      lnndfs  = nnods
      lkdofs  = nnods
      allocate(nndfs(lnndfs),kdofs(lkdofs))

      ! get array nndfs
      do inods = 1,nnods
         nndfs(inods) = nndf(isngns(inods))
      end do
! find local number of DOF on subdomain NDOFS
      ndofs = sum(nndfs)

      ! prepare array ISVGVN
      lisvgvns = ndofs
      allocate(isvgvns(lisvgvns))
      indvs = 0
      do inods = 1,nnods
         inod = isngns(inods)

         ndofn   = nndfs(inods)
         indvg   = kdof(inod)
         ! in lack of memory for kdof array, this can be altered by indvg = sum(nndf(1:inod-1))
         do idofn = 1,ndofn
            indvs = indvs + 1
            indvg = indvg + 1

            isvgvns(indvs) = indvg
         end do
      end do

! create array kdofs
      if (nnods.gt.0) then
         kdofs(1) = 0
         do inods = 2,nnods
            kdofs(inods) = kdofs(inods-1) + nndfs(inods-1)
         end do
      end if

      ! coordinates of subdomain nodes
      lxyzs1 = nnods
      lxyzs2 = ndim
      allocate(xyzs(lxyzs1,lxyzs2))
      do j = 1,ndim
         do i = 1,nnods
            indn = isngns(i)

            xyzs(i,j) = xyz(indn,j)
         end do
      end do

      call dd_upload_sub_mesh(sub, nelems, nnods, ndofs, ndim, &
                              nndfs,lnndfs, nnets,lnnets, inets,linets, isngns,lisngns, &
                              isvgvns,lisvgvns, isegns,lisegns,&
                              xyzs,lxyzs1,lxyzs2) 

      deallocate(nndfs,kdofs)
      deallocate(xyzs)
      deallocate(inets,nnets,isegns,isngns,isvgvns)

      deallocate(kdof)

end subroutine

!****************************************************
subroutine dd_localize_adj(sub,nsub,kadjsub,lkadjsub)
!****************************************************
! Subroutine for localization of adjacency array to particular subdomain
! loads the data directly to the structure
      use module_utils
      implicit none

! Subdomain structure
      type(subdomain_type),intent(inout) :: sub
! number of subdomains
      integer,intent(in) :: nsub
! keys of neighbours of subdomain (length nsub): 
! 0 - subdomain is not my neighbour 
! 1 - subdomain is my neighbour
      integer,intent(in) :: lkadjsub
      integer,intent(in) ::  kadjsub(lkadjsub)

      ! local vars
      integer :: nadjs
      integer ::            liadjs
      integer,allocatable :: iadjs(:)

      integer :: isub, jsub, indadjs

      ! check dimension
      if (nsub.ne.lkadjsub) then
         call error('DD_LOCALIZE_ADJ','array dimension mismatch.')
      end if
      if (nsub.ne.sub%nsub) then
         call error('DD_LOCALIZE_ADJ','number of subdomains mismatch.')
      end if

      isub = sub%isub

      nadjs = 0
      do jsub = 1,nsub
         if (kadjsub(jsub).eq.1) then
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
         if (kadjsub(jsub).eq.1) then
            indadjs = indadjs + 1
            iadjs(indadjs) = jsub
         end if
      end do
      ! debug
      !print *,'iadjs',iadjs
      call dd_upload_sub_adj(sub, nadjs, iadjs,liadjs)
      deallocate(iadjs)

end subroutine

!*********************************************************************************
subroutine dd_localize_cornersglobs(sub,ncorner,inodc,linodc,&
                                    nedge,nface,nnglb,lnnglb,inglb,linglb, nnodcs)
!*********************************************************************************
! Subroutine for localization of corners and globs to particular subdomain
! loads the data directly to the structure
      use module_utils
      implicit none

! Subdomain structure
      type(subdomain_type),intent(inout) :: sub
! number of corners
      integer,intent(in) :: ncorner
! global indices of corners
      integer,intent(in) :: linodc
      integer,intent(in) ::  inodc(linodc)
! number of edges
      integer,intent(in) :: nedge
! number of faces
      integer,intent(in) :: nface
! number of nodes in globs
      integer,intent(in) :: lnnglb
      integer,intent(in) ::  nnglb(lnnglb)
! global indices of nodes in globs (lenght sum(nnglb))
      integer,intent(in) :: linglb
      integer,intent(in) ::  inglb(linglb)

! number of coarse pseudonodes (corners and globs) at subdomain
      integer,intent(out) ::  nnodcs

      ! local vars
      integer :: nglb

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
      integer::            lignsins1, lignsins2
      integer,allocatable:: ignsins(:,:)
      integer::            ligvsivns1, ligvsivns2
      integer,allocatable:: igvsivns(:,:)
      integer ::            liingns
      integer,allocatable :: iingns(:)
      integer ::            lkdofs
      integer,allocatable :: kdofs(:)

      integer :: idofn, inods, iglb, iglbn, iglbv, iglobs, inc, indnc, &
                 indng, indns, indvs, inodcs, pointinglb
      integer :: indins, indivs, indn1
      integer :: ncorners, nglobs
      integer :: nnods, ndofs, nelems, nnodis, ndofis, ndofn, nglbn, nglbv

      ! check dimension
      if (ncorner.ne.linodc) then
         call error('DD_LOCALIZE_CORNERSGLOBS','array dimension mismatch INODC.')
      end if
      nglb = nedge + nface
      if (nglb.ne.lnnglb) then
         call error('DD_LOCALIZE_CORNERSGLOBS','array dimension mismatch NNGLB.')
      end if
      if (.not.sub%is_interface_loaded) then
         call error('DD_LOCALIZE_CORNERSGLOBS','Interface not loaded yet.')
      end if

      ! prepare array for global indices of interface nodes 
      call dd_get_interface_size(sub,ndofis,nnodis)
      liingns = nnodis
      allocate(iingns(liingns))
      call dd_get_interface_global_numbers(sub, iingns,liingns)

! find number of coarse nodes on subdomain NCORNERS
      ncorners = 0
      do inc = 1,ncorner
         indnc = inodc(inc)
         if (any(iingns.eq.indnc)) then
            ncorners = ncorners + 1
         end if
      end do

! create array kdofs
      call dd_get_size(sub, ndofs,nnods,nelems)
      lkdofs = nnods
      allocate(kdofs(lkdofs))
      if (nnods.gt.0) then
         kdofs(1) = 0
         do inods = 2,nnods
            kdofs(inods) = kdofs(inods-1) + sub%nndf(inods-1)
         end do
      end if
      !print *,'kdofs:',kdofs
      !call flush(6)


      ! find mapping of corners
      lglobal_corner_numbers = ncorners
      allocate(global_corner_numbers(lglobal_corner_numbers))
      licnsins = ncorners
      allocate(icnsins(licnsins))

      inodcs = 0
      do inc = 1,ncorner
         indnc = inodc(inc)
         if (any(iingns.eq.indnc)) then
            inodcs = inodcs + 1

            ! mapping to global corner numbers
            global_corner_numbers(inodcs) = inc
            ! mapping to subdomain interface numbers
            call get_index(indnc,iingns,nnodis,indins)
            if (indins .eq. -1) then
               write(*,*) 'Index of subdomain interface node not found.', indnc
               write(*,*) 'iingns',iingns
               call error_exit
            end if
            icnsins(inodcs) = indins
         end if
      end do

      call dd_upload_sub_corners(sub, ncorners, global_corner_numbers,lglobal_corner_numbers, icnsins,licnsins)
      deallocate(icnsins)
      deallocate(global_corner_numbers)

      ! find local number of globs NGLOBS
      nglobs     = 0
      pointinglb = 0
      do iglb = 1,nglb
         nglbn = nnglb(iglb)

         ! touch first node in glob
         indn1 = inglb(pointinglb + 1)

         if (any(iingns.eq.indn1)) then
            nglobs = nglobs + 1
         end if

         pointinglb = pointinglb + nglbn
      end do

      ! get array of interface 

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

         if (any(iingns.eq.indn1)) then

            iglobs = iglobs + 1

            nglbv = 0
            do iglbn = 1,nglbn
               indng = inglb(pointinglb + iglbn)
               call get_index(indng,iingns,liingns,indins)
               if (indins .eq. -1) then
                  write(*,*) ' Index of interface node not found for global ', indng
                  call error_exit
               end if
               indns = sub%iin(indins)

               ndofn = sub%nndf(indns)

               nglbv = nglbv + ndofn
            end do

            nglobvars(iglobs)   = nglbv
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

         if (any(iingns.eq.indn1)) then

            iglobs = iglobs + 1

            iglbv = 0
            do iglbn = 1,nglbn
               
               indng = inglb(pointinglb + iglbn)
               call get_index(indng,iingns,nnodis,indins)
               if (indins .eq. -1) then
                  write(*,*) ' Index of interface node not found.', indng
                  call error_exit
               end if

               ignsins(iglobs,iglbn) = indins

               indns = sub%iin(indins)
               ndofn = sub%nndf(indns)

               do idofn = 1,ndofn
                  iglbv = iglbv + 1

                  indvs = kdofs(indns) + idofn
                  call get_index(indvs,sub%iivsvn,sub%liivsvn,indivs)
                  if (indivs .eq. -1) then
                     write(*,*) 'DD_LOCALIZE_CORNERSGLOBS: Index of subdomain interface dof not found.'
                     write(*,*) 'indng =',indng,'indns =',indns,'indvs = ',indvs,'indivs = ',indivs, 'isub = ',sub%isub
                     call error_exit
                  end if

                  igvsivns(iglobs,iglbv) = indivs
               end do
            end do
         end if

         pointinglb = pointinglb + nglbn
      end do
      deallocate(iingns)
      deallocate(kdofs)

      call dd_upload_sub_globs(sub, nglobs, global_glob_numbers,lglobal_glob_numbers,&
                               nglobnodess,lnglobnodess, nglobvars,lnglobvars,&
                               ignsins,lignsins1,lignsins2, igvsivns,ligvsivns1,ligvsivns2,&
                               glob_types,lglob_types)
      deallocate(ignsins)
      deallocate(igvsivns)
      deallocate(glob_types)
      deallocate(nglobnodess)
      deallocate(nglobvars)
      deallocate(global_glob_numbers)

      nnodcs = ncorners + nglobs

end subroutine

!*********************************************************************************
subroutine dd_load_matrix_triplet(sub, matrixtype, &
                                  i_sparse,j_sparse,a_sparse,la,nnza,is_assembled)
!*********************************************************************************
! Subroutine for loading sparse triplet to sub structure
      use module_utils
      implicit none

! Subdomain structure
      type(subdomain_type),intent(inout) :: sub
! Type of the matrix
! 0 - unsymmetric
! 1 - symmetric positive definite
! 2 - general symmetric
      integer,intent(in) :: matrixtype

      ! Matrix in IJA sparse format
      integer,intent(in)  :: la
      integer,intent(in)  :: nnza
      integer,intent(in)  :: i_sparse(la), j_sparse(la)
      real(kr),intent(in) :: a_sparse(la)
      logical,intent(in)  :: is_assembled

      ! local vars
      integer :: i

      ! check that matrix is not loaded
      if (sub%is_matrix_loaded) then
         write(*,*) 'DD_LOAD_MATRIX_TRIPLET: Matrix already loaded for subdomain: ',sub%isub
         call error_exit
      end if

      ! load data
      sub%matrixtype = matrixtype
      sub%istorage = 2
      sub%nnza     = nnza
      sub%la       = la
      sub%is_assembled = is_assembled
      allocate(sub%i_a_sparse(la))
      do i = 1,la
         sub%i_a_sparse(i) = i_sparse(i)
      end do
      allocate(sub%j_a_sparse(la))
      do i = 1,la
         sub%j_a_sparse(i) = j_sparse(i)
      end do
      allocate(sub%a_sparse(la))
      do i = 1,la
         sub%a_sparse(i)   = a_sparse(i)
      end do

      sub%is_matrix_loaded = .true.
      sub%is_triplet       = .true.

end subroutine

!********************************************
subroutine dd_load_eliminated_bc(sub, bc,lbc)
!********************************************
! Subroutine for loading eliminated entries of matrix multiplied by fixed
! variable values
      use module_utils
      implicit none

! Subdomain structure
      type(subdomain_type),intent(inout) :: sub
      ! eliminated values of fixed variables
      integer,intent(in) :: lbc
      real(kr),intent(in)::  bc(lbc)

      ! local vars
      integer :: i

      ! check if I store the subdomain
      if (sub%is_bc_loaded) then
         write(*,*) 'DD_LOAD_ELIMINATED_BC: BC already loaded for subdomain ',sub%isub
      end if
      if (.not.sub%is_bc_nonzero) then 
         write(*,*) 'DD_LOAD_ELIMINATED_BC: subdomain ',sub%isub,&
                    ': loading BC for homogenous conditions.'
         call error_exit
      end if

      ! load eliminated boundary conditions if they are present
      if (lbc .gt. 0) then
         sub%lbc = lbc
         allocate(sub%bc(lbc))
         do i = 1,lbc
            sub%bc(i) = bc(i)
         end do
      end if

end subroutine

!*********************************************************************************
subroutine dd_upload_sub_mesh(sub, nelem, nnod, ndof, ndim, &
                              nndf,lnndf, nnet,lnnet, inet,linet, &
                              isngn,lisngn, isvgvn,lisvgvn, isegn,lisegn,&
                              xyz,lxyz1,lxyz2)
!*********************************************************************************
! Subroutine for loading mesh data into sub structure
      use module_utils
      implicit none

! Subdomain structure
      type(subdomain_type),intent(inout) :: sub
      ! mesh
      integer,intent(in) :: nelem, nnod, ndof, ndim
      integer,intent(in) :: lnndf,       lnnet,       linet
      integer,intent(in) ::  nndf(lnndf), nnet(lnnet), inet(linet)
      integer,intent(in) :: lisngn
      integer,intent(in) ::  isngn(lisngn)
      integer,intent(in) :: lisvgvn
      integer,intent(in) ::  isvgvn(lisvgvn)
      integer,intent(in) :: lisegn
      integer,intent(in) ::  isegn(lisegn)
      integer,intent(in) :: lxyz1, lxyz2
      real(kr),intent(in)::  xyz(lxyz1,lxyz2)

      ! local vars
      integer :: i, j

      if (.not.sub%is_initialized) then
         write(*,*) 'DD_UPLOAD_SUB_MESH: Not initialized subdomain: ',sub%isub
         call error_exit
      end if

      ! load data
      sub%nelem   = nelem
      sub%nnod    = nnod
      sub%ndof    = ndof
      sub%ndim    = ndim

      sub%linet   = linet
      allocate(sub%inet(linet))
      do i = 1,linet
         sub%inet(i) = inet(i)
      end do

      sub%lnnet   = lnnet
      allocate(sub%nnet(lnnet))
      do i = 1,lnnet
         sub%nnet(i) = nnet(i)
      end do

      sub%lnndf   = lnndf
      allocate(sub%nndf(lnndf))
      do i = 1,lnndf
         sub%nndf(i) = nndf(i)
      end do

      sub%lisngn   = lisngn
      allocate(sub%isngn(lisngn))
      do i = 1,lisngn
         sub%isngn(i) = isngn(i)
      end do

      sub%lisvgvn   = lisvgvn
      allocate(sub%isvgvn(lisvgvn))
      do i = 1,lisvgvn
         sub%isvgvn(i) = isvgvn(i)
      end do

      sub%lisegn   = lisegn
      allocate(sub%isegn(lisegn))
      do i = 1,lisegn
         sub%isegn(i) = isegn(i)
      end do

      sub%lxyz1   = lxyz1
      sub%lxyz2   = lxyz2
      allocate(sub%xyz(lxyz1,lxyz2))
      do j = 1,lxyz2
         do i = 1,lxyz1
            sub%xyz(i,j) = xyz(i,j)
         end do
      end do

      sub%is_mesh_loaded = .true.

end subroutine

!**************************************************
subroutine dd_upload_sub_adj(sub, nadj, iadj,liadj)
!**************************************************
! Subroutine for loading adjacency data into sub structure
      implicit none

! Subdomain structure
      type(subdomain_type),intent(inout) :: sub
      integer,intent(in) :: nadj
      integer,intent(in) :: liadj
      integer,intent(in) ::  iadj(liadj)

      ! local vars
      integer :: i

      ! load data
      sub%nadj = nadj
      sub%liadj = liadj
      allocate(sub%iadj(liadj))
      do i = 1,liadj
         sub%iadj(i) = iadj(i)
      end do

      sub%is_adj_loaded = .true.

end subroutine

!***************************************************************************
subroutine dd_upload_sub_interface(sub, nnodi, ndofi, ndofo, &
                                   iin,liin, iivsvn,liivsvn, iovsvn,liovsvn)
!***************************************************************************
! Subroutine for loading interface data into sub structure
      implicit none

! Subdomain structure
      type(subdomain_type),intent(inout) :: sub
      integer,intent(in) :: nnodi, ndofi, ndofo
      integer,intent(in) :: liin,       liivsvn,         liovsvn
      integer,intent(in) ::  iin(liin),  iivsvn(liivsvn), iovsvn(liovsvn)

      ! local vars
      integer :: i

      ! load data
      sub%nnodi   = nnodi
      sub%ndofi   = ndofi
      sub%ndofo   = ndofo

      sub%liin = liin
      allocate(sub%iin(liin))
      do i = 1,liin
         sub%iin(i) = iin(i)
      end do

      sub%liivsvn = liivsvn
      allocate(sub%iivsvn(liivsvn))
      do i = 1,liivsvn
         sub%iivsvn(i) = iivsvn(i)
      end do

      sub%liovsvn = liovsvn
      allocate(sub%iovsvn(liovsvn))
      do i = 1,liovsvn
         sub%iovsvn(i) = iovsvn(i)
      end do

      sub%is_interface_loaded = .true.

end subroutine

!*******************************************************************************************
subroutine dd_upload_sub_corners(sub, ncorner,&
                                 global_corner_number,lglobal_corner_number, icnsin,licnsin)
!*******************************************************************************************
! Subroutine for loading corner data into sub structure
      implicit none

! Subdomain structure
      type(subdomain_type),intent(inout) :: sub
      integer,intent(in) :: ncorner
      integer,intent(in) :: lglobal_corner_number,                       licnsin
      integer,intent(in) ::  global_corner_number(lglobal_corner_number), icnsin(licnsin)

      ! local vars
      integer :: i

      ! corner data
      sub%ncorner = ncorner

      sub%lglobal_corner_number = lglobal_corner_number
      allocate(sub%global_corner_number(lglobal_corner_number))
      do i = 1,lglobal_corner_number
         sub%global_corner_number(i) = global_corner_number(i)
      end do

      sub%licnsin = licnsin
      allocate(sub%icnsin(licnsin))
      do i = 1,licnsin
         sub%icnsin(i) = icnsin(i)
      end do

      sub%is_corners_loaded = .true.

end subroutine

!*************************************************************************************
subroutine dd_upload_sub_globs(sub, nglob, &
                               global_glob_number,lglobal_glob_number, &
                               nglobnodes,lnglobnodes, nglobvar,lnglobvar,&
                               ignsin,lignsin1,lignsin2, igvsivn,ligvsivn1,ligvsivn2,&
                               glob_type,lglob_type)
!*************************************************************************************
! Subroutine for loading globs data into sub structure
      implicit none

! Subdomain structure
      type(subdomain_type),intent(inout) :: sub
      integer,intent(in) :: nglob
      integer,intent(in) :: lglobal_glob_number,                     lnglobnodes,             lnglobvar
      integer,intent(in) ::  global_glob_number(lglobal_glob_number), nglobnodes(lnglobnodes), nglobvar(lnglobvar)
      integer,intent(in) :: lignsin1, lignsin2
      integer,intent(in) ::  ignsin(lignsin1,lignsin2)
      integer,intent(in) :: ligvsivn1, ligvsivn2
      integer,intent(in) ::  igvsivn(ligvsivn1,ligvsivn2)
      integer,intent(in) :: lglob_type
      integer,intent(in) ::  glob_type(lglob_type)

      ! local vars
      integer :: i, j

      ! load glob data
      sub%nglob = nglob

      sub%lglobal_glob_number = lglobal_glob_number
      allocate(sub%global_glob_number(lglobal_glob_number))
      do i = 1,lglobal_glob_number
         sub%global_glob_number(i) = global_glob_number(i)
      end do

      sub%lnglobnodes = lnglobnodes
      allocate(sub%nglobnodes(lnglobnodes))
      do i = 1,lnglobnodes
         sub%nglobnodes(i) = nglobnodes(i)
      end do

      sub%lnglobvar = lnglobvar
      allocate(sub%nglobvar(lnglobvar))
      do i = 1,lnglobvar
         sub%nglobvar(i) = nglobvar(i)
      end do

      sub%lignsin1 = lignsin1
      sub%lignsin2 = lignsin2
      allocate(sub%ignsin(lignsin1,lignsin2))
      do i = 1,lignsin1
         do j = 1,lignsin2
            sub%ignsin(i,j) = ignsin(i,j)
         end do
      end do

      sub%ligvsivn1 = ligvsivn1
      sub%ligvsivn2 = ligvsivn2
      allocate(sub%igvsivn(ligvsivn1,ligvsivn2))
      do i = 1,ligvsivn1
         do j = 1,ligvsivn2
            sub%igvsivn(i,j) = igvsivn(i,j)
         end do
      end do

      sub%lglob_type = lglob_type
      allocate(sub%glob_type(lglob_type))
      do i = 1,lglob_type
         sub%glob_type(i) = glob_type(i)
      end do

      sub%is_globs_loaded = .true.

end subroutine

!**********************************************************
subroutine dd_upload_bc(sub, ifix,lifix, fixv,lfixv)
!**********************************************************
! Subroutine for loading boundary conditions on subdomain
      implicit none

! Subdomain structure
      type(subdomain_type),intent(inout) :: sub
      integer,intent(in) :: lifix
      integer,intent(in) ::  ifix(lifix)
      integer,intent(in) :: lfixv
      real(kr),intent(in)::  fixv(lfixv)

      ! local vars
      integer :: i
      logical :: import_bc

      if (any(ifix.ne.0)) then
         sub%is_bc_present = .true.
         import_bc = .true.
      else
         sub%is_bc_present = .false.
         import_bc = .false.
      end if
      if (any(fixv.ne.0.0_kr)) then
         sub%is_bc_nonzero = .true.
      else
         sub%is_bc_nonzero = .false.
      end if

      ! boundary conditions
      if (import_bc) then
         sub%lifix = lifix
         allocate(sub%ifix(lifix))
         do i = 1,lifix
            sub%ifix(i) = ifix(i)
         end do
   
         sub%lfixv = lfixv
         allocate(sub%fixv(lfixv))
         do i = 1,lfixv
            sub%fixv(i) = fixv(i)
         end do

      end if
      sub%is_bc_loaded = .true.

end subroutine

!**************************************
subroutine dd_upload_rhs(sub, rhs,lrhs)
!**************************************
! Subroutine for initialization of subdomain right hand side
      implicit none

! Subdomain structure
      type(subdomain_type),intent(inout) :: sub
      integer,intent(in) :: lrhs
      real(kr),intent(in)::  rhs(lrhs)

      ! local vars
      integer :: i

      ! boundary conditions
      sub%lrhs = lrhs
      allocate(sub%rhs(lrhs))
      do i = 1,lrhs
         sub%rhs(i) = rhs(i)
      end do

      sub%is_rhs_loaded = .true.
   
end subroutine

!**********************************
subroutine dd_fix_bc(sub, vec,lvec)
!**********************************
! Subroutine for enforcing Dirichlet boundary conditions in vector VEC
      use module_utils
      implicit none

! Subdomain structure
      type(subdomain_type),intent(in) :: sub
      integer,intent(in) ::    lvec
      real(kr),intent(inout)::  vec(lvec)

      ! check if there are nonzero BC at subdomain
      if (.not. sub%is_bc_present) then
         return
      end if
      if (.not. sub%is_bc_loaded) then
         call error('DD_FIX_BC','BC not yet loaded')
      end if
      ! check size
      if (sub%lifix .ne. lvec) then
         write(*,*) 'DD_FIX_BC: Vector size mismatch for subdomain ',sub%isub
         call error_exit
      end if

      ! enforce boundary conditions
      where (sub%ifix .ne. 0) vec = sub%fixv

end subroutine

!*************************************************
subroutine dd_fix_bc_interface_dual(sub, vec,lvec)
!*************************************************
! Subroutine for enforcing Dirichlet boundary conditions in interface residual vector
      use module_utils
      implicit none

! Subdomain structure
      type(subdomain_type),intent(in) :: sub
      integer,intent(in) ::    lvec
      real(kr),intent(inout)::  vec(lvec)

      ! local vars
      integer :: ndofi
      integer :: i, ind

      ! check if there are nonzero BC at subdomain
      if (.not. sub%is_bc_present) then
         return
      end if
      ! check size
      ndofi = sub%ndofi
      if (ndofi .ne. lvec) then
         write(*,*) 'DD_FIX_BC_INTERFACE_DUAL: Vector size mismatch: ',sub%isub
         call error_exit
      end if

      ! enforce boundary conditions
      do i = 1,ndofi
         ind = sub%iivsvn(i)
         if (sub%ifix(ind) .ne. 0) then
            vec(i) = 0._kr
         end if
      end do

end subroutine

!*******************************************************
subroutine dd_map_sub_to_subi(sub, vec,lvec, veci,lveci)
!*******************************************************
! Subroutine that maps subdomain vector to interface vector
      use module_utils
      implicit none

! Subdomain structure
      type(subdomain_type),intent(in) :: sub

      integer,intent(in)  ::  lvec
      real(kr),intent(in) ::   vec(lvec)

      integer,intent(in)  ::  lveci
      real(kr),intent(out) ::  veci(lveci)

      ! local vars
      integer :: i, ind
      integer :: ndofi, ndof

      ! check if mesh is loaded
      if (.not. sub%is_interface_loaded) then
         write(*,*) 'DD_MAP_SUB_TO_SUBI: Interface not loaded for subdomain: ',sub%isub
         call error_exit
      end if

      ! get dimensions
      ndof  = sub%ndof
      ndofi = sub%ndofi

      ! check dimensions
      if (lvec .ne. ndof .or. lveci .ne. ndofi) then
         write(*,*) 'DD_MAP_SUB_TO_SUBI: Vectors size mismatch, isub : ',sub%isub
         call error_exit
      end if

      do i = 1,ndofi
         ind = sub%iivsvn(i)

         veci(i) = vec(ind)
      end do

end subroutine

!***************************************************************
subroutine dd_map_glob_to_sub_int(sub, ivec,livec, ivecs,livecs)
!***************************************************************
! Subroutine that maps global vector of INTEGERSs to subdomain vector
      use module_utils
      implicit none

! Subdomain structure
      type(subdomain_type),intent(in) :: sub

      ! global array
      integer,intent(in) ::  livec
      integer,intent(in) ::   ivec(livec)

      ! subdomain array
      integer,intent(in)  ::  livecs
      integer,intent(out) ::   ivecs(livecs)

      ! local vars
      integer :: i, ind
      integer :: ndofs

      ! check if mesh is loaded
      if (.not. sub%is_mesh_loaded) then
         write(*,*) 'DD_MAP_GLOB_TO_SUB_INT: Mesh not loaded for subdomain: ',sub%isub
         call error_exit
      end if

      ! get dimensions
      ndofs = sub%ndof

      ! check dimensions
      if (livecs .ne. ndofs) then
         write(*,*) 'DD_MAP_GLOB_TO_SUB_INT: Vectors size mismatch, isub : ',sub%isub
         call error_exit
      end if

      do i = 1,ndofs
         ind = sub%isvgvn(i)

         ivecs(i) = ivec(ind)
      end do

end subroutine

!************************************************************
subroutine dd_map_glob_to_sub_real(sub, vec,lvec, vecs,lvecs)
!************************************************************
! Subroutine that maps global vector of REALs to subdomain vector
      use module_utils
      implicit none

! Subdomain structure
      type(subdomain_type),intent(in) :: sub

      ! global array
      integer,intent(in)  ::  lvec
      real(kr),intent(in) ::   vec(lvec)

      ! subdomain array
      integer,intent(in)  ::  lvecs
      real(kr),intent(out) ::  vecs(lvecs)

      ! local vars
      integer :: i, ind
      integer :: ndofs

      ! check if mesh is loaded
      if (.not. sub%is_mesh_loaded) then
         write(*,*) 'DD_MAP_GLOB_TO_SUB_REAL: Mesh not loaded for subdomain ',sub%isub
         call error_exit
      end if

      ! get dimensions
      ndofs = sub%ndof

      ! check dimensions
      if (lvecs .ne. ndofs) then
         write(*,*) 'DD_MAP_GLOB_TO_SUB_REAL: Vectors size mismatch for subdomain ',sub%isub
         call error_exit
      end if

      do i = 1,ndofs
         ind = sub%isvgvn(i)

         vecs(i) = vec(ind)
      end do

end subroutine

!*******************************************************
subroutine dd_map_sub_to_glob(sub, vecs,lvecs, vec,lvec)
!*******************************************************
! Subroutine that maps subdomain vector to global vector
      use module_utils
      implicit none

! Subdomain structure
      type(subdomain_type),intent(in) :: sub

      ! subdomain array
      integer,intent(in)  ::  lvecs
      real(kr),intent(in) ::   vecs(lvecs)

      ! global array
      integer,intent(in)  ::   lvec
      real(kr),intent(inout) :: vec(lvec)

      ! local vars
      integer :: i, ind
      integer :: ndofs

      ! check if mesh is loaded
      if (.not. sub%is_mesh_loaded) then
         write(*,*) 'DD_MAP_SUB_TO_GLOB: Mesh not loaded for subdomain: ',sub%isub
         call error_exit
      end if

      ! get dimensions
      ndofs = sub%ndof

      ! check dimensions
      if (lvecs .ne. ndofs) then
         write(*,*) 'DD_MAP_SUB_TO_GLOB: Vectors size mismatch, isub : ',sub%isub
         call error_exit
      end if

      do i = 1,ndofs
         ind = sub%isvgvn(i)

         vec(ind) = vec(ind) + vecs(i)
      end do

end subroutine

!*******************************************************
subroutine dd_map_subi_to_sub(sub, veci,lveci, vec,lvec)
!*******************************************************
! Subroutine that maps interface vector to subdomain vector
      use module_utils
      implicit none

! Subdomain structure
      type(subdomain_type),intent(in) :: sub

      integer,intent(in)  ::  lveci
      real(kr),intent(in) ::   veci(lveci)

      integer,intent(in)  ::  lvec
      real(kr),intent(out) ::  vec(lvec)

      ! local vars
      integer :: i, ind
      integer :: ndofi, ndof

      ! check pre
      if (.not. sub%is_interface_loaded) then
         write(*,*) 'DD_MAP_SUBI_TO_SUB: Interface not loaded for subdomain: ',sub%isub
         call error_exit
      end if

      ! get dimensions
      ndofi = sub%ndofi
      ndof  = sub%ndof

      ! check dimensions
      if (lvec .ne. ndof .or. lveci .ne. ndofi) then
         write(*,*) 'DD_MAP_SUBI_TO_SUB: Vectors size mismatch, isub : ',sub%isub
         call error_exit
      end if

      do i = 1,ndofi
         ind = sub%iivsvn(i)

         vec(ind) = vec(ind) + veci(i)
      end do

end subroutine

!*******************************************************
subroutine dd_map_sub_to_subo(sub, vec,lvec, veco,lveco)
!*******************************************************
! Subroutine that maps subdomain vector to interiOr vector
      use module_utils
      implicit none

! Subdomain structure
      type(subdomain_type),intent(in) :: sub

      integer,intent(in)  ::  lvec
      real(kr),intent(in) ::   vec(lvec)

      integer,intent(in)  ::  lveco
      real(kr),intent(out) ::  veco(lveco)

      ! local vars
      integer :: i, ind
      integer :: ndofo, ndof

      ! check pre
      if (.not. sub%is_interface_loaded) then
         write(*,*) 'DD_MAP_SUB_TO_SUBO: Mesh not loaded for subdomain: ',sub%isub
         call error_exit
      end if

      ! get dimensions
      ndof  = sub%ndof
      ndofo = sub%ndofo

      ! check dimensions
      if (lvec .ne. ndof .or. lveco .ne. ndofo) then
         write(*,*) 'DD_MAP_SUB_TO_SUB0: Vectors size mismatch, isub : ',sub%isub
         call error_exit
      end if

      do i = 1,ndofo
         ind = sub%iovsvn(i)

         veco(i) = vec(ind)
      end do

end subroutine

!***************************************************************
subroutine dd_map_subc_to_globc(sub, vecsc,lvecsc, vecgc,lvecgc)
!***************************************************************
! Subroutine that maps subdomain coarse vector to global coarse vector
      use module_utils
      implicit none

! Subdomain structure
      type(subdomain_type),intent(in) :: sub

      integer,intent(in)  ::  lvecsc
      real(kr),intent(in) ::   vecsc(lvecsc)

      integer,intent(in)  ::  lvecgc
      real(kr),intent(out) ::  vecgc(lvecgc)

      ! local vars
      integer :: i, ind
      integer :: ndofc

      ! check if mesh is loaded
      if (.not. sub%is_c_loaded) then
         write(*,*) 'DD_MAP_SUBC_TO_GLOBC: Embedding not loaded for subdomain: ',sub%isub
         call error_exit
      end if

      ! get dimensions
      ndofc = sub%ndofc

      ! check dimensions
      if (lvecsc .ne. ndofc .or. maxval(sub%indrowc) .gt. lvecgc) then
         write(*,*) 'DD_MAP_SUBC_TO_GLOBC: Vectors size mismatch, isub : ',sub%isub
         call error_exit
      end if

      ! embed the arrays
      ! vecgc = vecgc + Rc_i * vecsc
      do i = 1,ndofc
         ind = sub%indrowc(i)

         vecgc(ind) = vecgc(ind) + vecsc(i)
      end do

end subroutine

!***************************************************************
subroutine dd_map_globc_to_subc(sub, vecgc,lvecgc, vecsc,lvecsc)
!***************************************************************
! Subroutine that maps global coarse vector subdomain coarse vector to subdomain coarse vector

      use module_utils
      implicit none

! Subdomain structure
      type(subdomain_type),intent(in) :: sub

      integer,intent(in)  ::  lvecgc
      real(kr),intent(in) ::   vecgc(lvecgc)

      integer,intent(in)  ::  lvecsc
      real(kr),intent(out) ::  vecsc(lvecsc)

      ! local vars
      integer :: i, ind
      integer :: ndofc

      ! check if corse mesh is loaded
      if (.not. sub%is_c_loaded) then
         write(*,*) 'DD_MAP_SUBC_TO_GLOBC: Embedding not loaded for subdomain: ',sub%isub
         call error_exit
      end if

      ! get dimensions
      ndofc = sub%ndofc

      ! check dimensions
      if (lvecsc .ne. ndofc .or. maxval(sub%indrowc) .gt. lvecgc) then
         write(*,*) 'DD_MAP_SUBC_TO_GLOBC: Vectors size mismatch, isub : ',sub%isub
         call error_exit
      end if

      ! embed the arrays
      ! vecgc = vecgc + Rc_i * vecsc
      do i = 1,ndofc
         ind = sub%indrowc(i)

         vecsc(i) = vecgc(ind)
      end do

end subroutine

!***************************************
subroutine dd_assembly_local_matrix(sub)
!***************************************
! Subroutine for assemblage of matrix
      use module_utils
      use module_sm
      implicit none

! Subdomain structure
      type(subdomain_type),intent(inout) :: sub

      ! check the prerequisities
      if (.not.sub%is_matrix_loaded) then
         write(*,*) 'DD_ASSEMBLY_LOCAL_MATRIX: Matrix is not loaded for subdomain:', sub%isub
         call error_exit
      end if
      if (sub%is_assembled) then
         write(*,*) 'DD_ASSEMBLY_LOCAL_MATRIX: Matrix already assembled for subdomain:', sub%isub
         return 
      end if

      ! call matrix assemblage
      if (sub%istorage.eq.1 .or. sub%istorage.eq.2) then
         ! assembly triplet
         call sm_assembly(sub%i_a_sparse, sub%j_a_sparse, sub%a_sparse, sub%la, sub%nnza)
      else
         write(*,*) 'DD_ASSEMBLY_LOCAL_MATRIX: Matrix in unknown format type:',sub%istorage
         call error_exit
      end if
end subroutine

!*****************************************************
subroutine dd_matrix_tri2blocktri(sub,remove_original)
!*****************************************************
! Subroutine for conversion of sparse triplet into block format A11, A12, A21, A22
      use module_utils
      implicit none

! Subdomain structure
      type(subdomain_type),intent(inout) :: sub
      ! deallocate original matrix?
      logical,intent(in) :: remove_original

      ! local vars
      integer :: la, nnza, ndof, ndofi, ndofo
      logical :: is_symmetric_storage
      integer :: ia, inddofi, inddofo, idofi, idofo, irow, jcol
      integer :: la11, la12, la21, la22 
      integer :: ia11, ia12, ia21, ia22 
      real(kr) :: aval

      ! mask for interface entries
      integer ::            lmaski,   lmasko
      integer,allocatable :: maski(:), masko(:)

      ! check the prerequisities
      if (.not.sub%is_interface_loaded) then
         write(*,*) 'DD_MATRIX_TRI2BLOCKTRI: Interface is not loaded for subdomain:', sub%isub
         call error_exit
      end if
      if (.not.sub%is_matrix_loaded) then
         write(*,*) 'DD_MATRIX_TRI2BLOCKTRI: Matrix is not loaded for subdomain:', sub%isub
         call error_exit
      end if
      if (sub%istorage.eq.3 .or. sub%istorage.eq.4) then
         write(*,*) 'DD_MATRIX_TRI2BLOCKTRI: Matrix already in block triple format for subdomain:', sub%isub
         return
      end if
      if (.not.(sub%istorage.eq.1 .or. sub%istorage.eq.2)) then
         write(*,*) 'DD_MATRIX_TRI2BLOCKTRI: Matrix not in single block triple format, subdomain:', sub%isub,&
                    'type:',sub%istorage
         call error_exit
      end if

      if (sub%istorage .eq. 1) then
         is_symmetric_storage = .false.
      else
         is_symmetric_storage = .true.
      end if

      ! prepare mask of interface and interior nodes
      ndof   = sub%ndof
      lmaski = ndof
      lmasko = ndof
      allocate(maski(lmaski),masko(lmasko))
      call zero(maski,lmaski)
      call zero(masko,lmasko)

      ndofi  = sub%ndofi
      do idofi = 1,ndofi
         inddofi = sub%iivsvn(idofi)
         maski(inddofi) = idofi
      end do
      ndofo  = sub%ndofo
      do idofo = 1,ndofo
         inddofo = sub%iovsvn(idofo)
         masko(inddofo) = idofo 
      end do

      if (debug) then
         ! check consistency of masks
         if (any(maski.ne.0 .and. masko.ne.0)) then
            write(*,*) 'DD_MATRIX_TRI2BLOCKTRI: Check of mask consistency failed for subdomain:',sub%isub
            call error_exit
         end if
         if (any(maski.eq.0 .and. masko.eq.0)) then
            write(*,*) 'DD_MATRIX_TRI2BLOCKTRI: Check of mask coverage failed for subdomain:',sub%isub
            call error_exit
         end if
      end if

      ! count sizes of matrix blocks in question
      la   = sub%la
      nnza = sub%nnza
      la11 = 0
      la12 = 0
      la21 = 0
      la22 = 0
      do ia = 1,nnza
         irow = sub%i_a_sparse(ia)
         jcol = sub%j_a_sparse(ia)

         ! diagonal blocks
         if      (maski(irow).eq.0 .and. maski(jcol).eq.0) then
            la11 = la11 + 1
         else if (maski(irow).ne.0 .and. maski(jcol).ne.0) then
            la22 = la22 + 1
         ! offdiagonal blocks
         else if (maski(irow).eq.0 .and. maski(jcol).ne.0) then
            la12 = la12 + 1
         else if (maski(irow).ne.0 .and. maski(jcol).eq.0) then
            la21 = la21 + 1
         end if
      end do

      ! prepare space for blocks of matrix
      sub%la11   = la11
      sub%nnza11 = la11
      sub%la22   = la22
      sub%nnza22 = la22
      if (is_symmetric_storage) then
         la12 = la12 + la21
         sub%la12   = la12
         sub%nnza12 = la12
         sub%la21   = 0
         sub%nnza21 = 0
      else
         sub%la12   = la12
         sub%nnza12 = la12
         sub%la21   = la21
         sub%nnza21 = la21
      end if
      allocate(sub%i_a11_sparse(la11),sub%j_a11_sparse(la11),sub%a11_sparse(la11))
      allocate(sub%i_a22_sparse(la22),sub%j_a22_sparse(la22),sub%a22_sparse(la22))

      if (is_symmetric_storage) then
         allocate(sub%i_a12_sparse(la12),sub%j_a12_sparse(la12),sub%a12_sparse(la12))
      else
         allocate(sub%i_a12_sparse(la12),sub%j_a12_sparse(la12),sub%a12_sparse(la12))
         allocate(sub%i_a21_sparse(la21),sub%j_a21_sparse(la21),sub%a21_sparse(la21))
      end if

      ! convert matrix to blocks according to interface - denoted by 1 in MASKI array
      ia11 = 0
      ia12 = 0
      ia21 = 0
      ia22 = 0
      do ia = 1,nnza
         irow = sub%i_a_sparse(ia)
         jcol = sub%j_a_sparse(ia)
         aval = sub%a_sparse(ia)

         ! diagonal blocks
         if      (maski(irow).eq.0 .and. maski(jcol).eq.0) then
            ia11 = ia11 + 1
            sub%i_a11_sparse(ia11) = masko(irow)
            sub%j_a11_sparse(ia11) = masko(jcol)
            sub%a11_sparse(ia11)   = aval
         else if (maski(irow).ne.0 .and. maski(jcol).ne.0) then
            ia22 = ia22 + 1
            sub%i_a22_sparse(ia22) = maski(irow)
            sub%j_a22_sparse(ia22) = maski(jcol)
            sub%a22_sparse(ia22)   = aval
         ! offdiagonal blocks
         else if (maski(irow).eq.0 .and. maski(jcol).ne.0) then
            ia12 = ia12 + 1
            sub%i_a12_sparse(ia12) = masko(irow)
            sub%j_a12_sparse(ia12) = maski(jcol)
            sub%a12_sparse(ia12)   = aval
         else if (maski(irow).ne.0 .and. maski(jcol).eq.0) then
            if (is_symmetric_storage) then
               ia12 = ia12 + 1
               sub%i_a12_sparse(ia12) = masko(jcol)
               sub%j_a12_sparse(ia12) = maski(irow)
               sub%a12_sparse(ia12)   = aval
            else
               ia21 = ia21 + 1
               sub%i_a21_sparse(ia21) = maski(irow)
               sub%j_a21_sparse(ia21) = masko(jcol)
               sub%a21_sparse(ia21)   = aval
            end if
         end if
      end do

      if (is_symmetric_storage) then
         sub%istorage = 4
      else
         sub%istorage = 3
      end if

      sub%is_blocked = .true.

      deallocate(maski,masko)

      if (remove_original) then
         deallocate(sub%i_a_sparse, sub%j_a_sparse, sub%a_sparse)
         sub%la   = 0
         sub%nnza = 0
         sub%is_triplet = .false.
      end if
end subroutine

!*****************************************
subroutine dd_prepare_schur(sub,comm_self)
!*****************************************
! Subroutine for preparing data for computing with reduced problem
      use module_mumps
      use module_utils
      implicit none

! Subdomain structure
      type(subdomain_type),intent(inout) :: sub
      ! communicator
      integer,intent(in) :: comm_self

      ! local vars
      integer :: ndofo, la11, nnza11
      integer :: mumpsinfo
      logical :: parallel_analysis

      ! check the prerequisities
      if (.not. sub%is_blocked) then
         write(*,*) 'DD_PREPARE_SCHUR: Matrix is not in blocked format. isub = ',sub%isub
         call error_exit
      end if

      ! Initialize MUMPS
      call mumps_init(sub%mumps_interior_block,comm_self,sub%matrixtype)
      ! Level of information from MUMPS
      if (debug) then
         mumpsinfo = 2
      else
         mumpsinfo = 0
      end if
      call mumps_set_info(sub%mumps_interior_block,mumpsinfo)
      ! Load matrix to MUMPS
      ndofo  = sub%ndofo
      nnza11 = sub%nnza11
      la11   = sub%la11
      call mumps_load_triplet(sub%mumps_interior_block,ndofo,nnza11,&
                              sub%i_a11_sparse,sub%j_a11_sparse,sub%a11_sparse,la11)
      ! Analyze matrix
      parallel_analysis = .false.
      call mumps_analyze(sub%mumps_interior_block,parallel_analysis) 
      ! Factorize matrix 
      call mumps_factorize(sub%mumps_interior_block) 

      sub%is_mumps_interior_active = .true.
      sub%is_interior_factorized = .true.

end subroutine

!***************************************************
subroutine dd_load_arithmetic_constraints(sub,itype)
!***************************************************
! Subroutine for assemblage of matrix
      use module_utils
      implicit none
! Subdomain structure
      type(subdomain_type),intent(inout) :: sub

      integer,intent(in) :: itype ! type of globs (2 - edges, 1 - faces)

      ! local vars
      integer :: icnode, lmatrix1, lmatrix2, ncnodes, ncdof, nvar
      integer :: i, j

      ! check the prerequisities
      if (.not.allocated(sub%cnodes) .or. .not.sub%is_cnodes_loaded) then
         write(*,*) 'DD_LOAD_ARITHMETIC_CONSTRAINTS: Array for cnodes not ready.'
         call error_exit
      end if

      ! get number of coarse nodes
      ncnodes = sub%ncnodes

      ! get number of constraints on an arithmetic constraint
      ! ndim for arithmetic averages
      ncdof = sub%ndim

      ! generate arithmetic averages on coarse nodes of prescribed type (e.g. edges)
      do icnode = 1,ncnodes
         if (sub%cnodes(icnode)%itype .eq. itype) then
         
            nvar = sub%cnodes(icnode)%nvar

            lmatrix1 = ncdof
            lmatrix2 = nvar

            sub%cnodes(icnode)%ncdof    = ncdof
            sub%cnodes(icnode)%nnz      = nvar
            sub%cnodes(icnode)%lmatrix1 = lmatrix1
            sub%cnodes(icnode)%lmatrix2 = lmatrix2
            allocate(sub%cnodes(icnode)%matrix(lmatrix1,lmatrix2))

            call zero(sub%cnodes(icnode)%matrix,lmatrix1,lmatrix2)

            do i = 1,lmatrix1
               do j = i,lmatrix2,lmatrix1
                  sub%cnodes(icnode)%matrix(i,j) = 1._kr
               end do
            end do
            sub%cnodes(icnode)%arithmetic = .true.

            ! mark the coarse node as used
            sub%cnodes(icnode)%used = .true.
         end if
      end do

end subroutine

!*********************************************************************************
subroutine dd_load_adaptive_constraints(sub,gglob,cadapt,lcadapt1,lcadapt2,nvalid)
!*********************************************************************************
! Subroutine for assemblage of matrix
      use module_utils
      implicit none
! Subdomain structure
      type(subdomain_type),intent(inout) :: sub
      integer,intent(in) :: gglob

      integer,intent(in) :: lcadapt1, lcadapt2
      real(kr),intent(in) :: cadapt(lcadapt1,lcadapt2)

      integer,intent(out) :: nvalid ! number of valid constraints after regularization

! local variables

! Variables for regularization of constraints
      integer              :: lavg1, lavg2
      real(kr),allocatable ::  avg(:,:)

! LAPACK variables
      integer             :: lipiv
      integer,allocatable ::  ipiv(:)
      integer              :: lwork
      real(kr),allocatable ::  work(:)
      integer              :: ltau
      real(kr),allocatable ::  tau(:)
      integer             :: lapack_info, ldavg

      real(kr) :: normval, thresh_diag


      ! local vars
      integer :: ind_loc, lmatrix1, lmatrix2, nvarglb, ncdof
      integer :: i, j, indiv

      ! check the prerequisities
      if (.not.allocated(sub%cnodes) .or. .not.sub%is_cnodes_loaded) then
         write(*,*) 'DD_LOAD_ADAPTIVE_CONSTRAINTS: Array for cnodes not ready.'
         call error_exit
      end if

      ! find local (subdomain) index of the glob from its global number
      call get_index(gglob,sub%cnodes%global_cnode_number,sub%ncnodes,ind_loc)
      if (ind_loc.le.0) then
         write(*,*) 'DD_LOAD_ADAPTIVE_CONSTRAINTS: Index of glob not found!'
         call error_exit
      end if

      ! prepare space for these constraints in the structure
      nvarglb = sub%cnodes(ind_loc)%nvar 

      if (nvarglb.gt.lcadapt1) then
         write(*,*) 'DD_LOAD_ADAPTIVE_CONSTRAINTS: Dimension of matrix of averages is smaller than number of variables on glob.'
         write(*,*) 'DD_LOAD_ADAPTIVE_CONSTRAINTS: nvarglb =',nvarglb, 'lcadapt1 =',lcadapt1
         call error_exit
      end if

      ! regularize the matrix of constraints by QR decomposition
      ! copy part of constraints to array AVG
      lavg1 = nvarglb
      lavg2 = lcadapt2
      allocate(avg(lavg1,lavg2))
      
      do i = 1,nvarglb
         indiv = sub%cnodes(ind_loc)%ivsivn(i)

         do j = 1,lcadapt2
            avg(i,j) = cadapt(indiv,j)
         end do
      end do

      !write (*,*) 'AVG before QR'
      !do i = 1,lavg1
      !   write(*,*) (avg(i,j),j = 1,lavg2)
      !end do

      ! perform QR decomposition of AVG by LAPACK
      ! Prepare array for permutations
      lipiv = lavg2
      allocate(ipiv(lipiv))
      ipiv = 0
      ! prepare other LAPACK arrays
      ltau = lavg1
      allocate(tau(ltau))
      lwork = 3*lavg2 + 1
      allocate(work(lwork))

      ldavg = max(1,lavg1)
      ! QR decomposition
      call DGEQP3( lavg1, lavg2, avg, ldavg, ipiv, tau, work, lwork, lapack_info )

      !write (*,*) 'AVG after QR factorization'
      !do i = 1,lavg1
      !   write(*,*) (avg(i,j),j = 1,lavg2)
      !end do
      !write (*,*) 'IPIV after QR factorization'
      !write(*,*) (ipiv(j),j = 1,lavg2)

      ! determine number of columns to use
      ! threshold of 1% of maximal norm
      nvalid = 0
      if (lavg1.gt.0.and.lavg2.gt.0) then
         normval = abs(avg(1,1))
         if (normval.gt.numerical_zero) then
            thresh_diag = 0.01_kr * normval
            do i = 1,lavg2
               if (abs(avg(i,i)) .lt. thresh_diag) then
                  exit
               else
                  nvalid = i
               end if
            end do
         end if
      end if

      !write (*,*) 'Number of constraints to really use nvalid',nvalid

      ! construct Q in AVG array
      call DORGQR( lavg1, nvalid, nvalid, avg, ldavg, tau, work, lwork, lapack_info )

      deallocate(work)
      deallocate(tau)
      deallocate(ipiv)

      !write (*,*) 'AVG contains Q'
      !do i = 1,lavg1
      !   write(*,*) (avg(i,j),j = 1,nvalid)
      !end do

      if (debug) then
         if (nvalid.lt.lcadapt2) then
            write(*,*) 'DD_LOAD_ADAPTIVE_CONSTRAINTS: WARNING - almost linearly dependent constraints on glob ',gglob
            write(*,*) 'DD_LOAD_ADAPTIVE_CONSTRAINTS:           Number of constrains reduced from ',lcadapt2,' to ',nvalid
         end if
      end if

      ! copy transposed selected regularized constraints to the global structure
      ncdof = nvalid
      sub%cnodes(ind_loc)%ncdof = ncdof

      lmatrix1 = nvalid
      lmatrix2 = nvarglb
      sub%cnodes(ind_loc)%lmatrix1 = lmatrix1
      sub%cnodes(ind_loc)%lmatrix2 = lmatrix2
      allocate(sub%cnodes(ind_loc)%matrix(lmatrix1,lmatrix2))
      
      do i = 1,nvalid
         do j = 1,nvarglb

            sub%cnodes(ind_loc)%matrix(i,j) = avg(j,i)
         end do
      end do
      sub%cnodes(ind_loc)%nnz = lmatrix1*lmatrix2
      sub%cnodes(ind_loc)%adaptive   = .true.
      sub%cnodes(ind_loc)%arithmetic = .false.

      ! mark the coarse node as used
      sub%cnodes(ind_loc)%used = .true.

      deallocate(avg)

      if (debug) then
         write(*,*) 'DD_LOAD_ADAPTIVE_CONSTRAINTS: Loading matrix of subdomain ',sub%isub,&
                    ' local glob #',ind_loc
         do i = 1,lmatrix1
            print '(100f15.6)',(sub%cnodes(ind_loc)%matrix(i,j),j = 1,lmatrix2)
         end do
      end if
end subroutine

!****************************************************************
subroutine dd_get_adaptive_constraints_size(sub,iglb,lavg1,lavg2)
!****************************************************************
! Subroutine for inquiring sizes to allocate for number of adaptive averages
      use module_utils
      implicit none
! Subdomain structure
      type(subdomain_type),intent(in) :: sub
      ! local (subdomain) glob number
      integer,intent(in) :: iglb
      ! sizes of matrix of averages
      integer,intent(out) :: lavg1, lavg2

      ! check the prerequisities
      if (.not.allocated(sub%cnodes) .or. .not.sub%is_cnodes_loaded) then
         write(*,*) 'DD_GET_ADAPTIVE_CONSTRAINTS_SIZE: Array for coarse nodes constraints not ready.'
         call error_exit
      end if
      if (.not.allocated(sub%cnodes(iglb)%matrix)) then
         write(*,*) 'isub',sub%isub,'iglb',iglb,'type',sub%cnodes(iglb)%itype
         write(*,*) 'DD_GET_ADAPTIVE_CONSTRAINTS_SIZE: Matrix of constraints not allocated.'
         call error_exit
      end if

      lavg1 = sub%cnodes(iglb)%lmatrix1
      lavg2 = sub%cnodes(iglb)%lmatrix2
end subroutine

!***************************************************************
subroutine dd_get_adaptive_constraints(sub,iglb,avg,lavg1,lavg2)
!***************************************************************
! Subroutine for getting adaptive averages from the structure
      use module_utils
      implicit none
! Subdomain structure
      type(subdomain_type),intent(in) :: sub
      ! local (subdomain) glob number
      integer,intent(in) :: iglb

      ! matrix of averages
      integer,intent(in) ::  lavg1, lavg2
      real(kr),intent(out) :: avg(lavg1,lavg2)

      ! local variables
      integer :: i,j

      ! check the prerequisities
      if (.not.allocated(sub%cnodes) .or. .not.sub%is_cnodes_loaded) then
         write(*,*) 'DD_GET_ADAPTIVE_CONSTRAINTS: Array for coarse nodes constraints not ready.'
         call error_exit
      end if
      if (sub%cnodes(iglb)%lmatrix1.ne.lavg1 .or. &
          sub%cnodes(iglb)%lmatrix2.ne.lavg2) then
         write(*,*) 'DD_GET_ADAPTIVE_CONSTRAINTS: Matrix dimensions for averages do not match.'
         call error_exit
      end if

      do j = 1,lavg2
         do i = 1,lavg1
            avg(i,j) = sub%cnodes(iglb)%matrix(i,j)
         end do
      end do

end subroutine

!**********************************
subroutine dd_construct_cnodes(sub)
!**********************************
! Merging corners with globs - order first corners, then globs
      use module_utils
      implicit none
! Subdomain structure
      type(subdomain_type),intent(inout) :: sub

      ! local vars
      integer ::             ncorner
      integer ::             nglob
      integer ::             ncnodes, lcnodes
      integer ::             icnode, igcnode
      integer ::             inodc, indnode, indinode, i, nvar, iglob, nnodgl, indn, indin
      integer ::             lxyz

      ! check the prerequisities
      if (.not.sub%is_corners_loaded) then
         write(*,*) 'DD_GET_CNODES: Corners not loaded for subdomain:', sub%isub
         call error_exit
      end if
      if (.not.sub%is_globs_loaded) then
         write(*,*) 'DD_GET_CNODES: Globs not loaded for subdomain:', sub%isub
         call error_exit
      end if

      ! determine number of coarse nodes
      ncorner = sub%ncorner
      nglob   = sub%nglob
      ncnodes = ncorner + nglob
      sub%ncnodes = ncnodes

      lcnodes = ncnodes
      allocate(sub%cnodes(lcnodes))

      ! set counter
      icnode = 0

      ! copy corners
      do inodc = 1,ncorner
         icnode = icnode + 1

         indinode = sub%icnsin(inodc)
         indnode  = sub%iin(indinode)

         ! type of coarse node - corner
         sub%cnodes(icnode)%itype = 3
         ! is coarse node used?
         sub%cnodes(icnode)%used  = .true.
         ! global number
         igcnode =  sub%global_corner_number(inodc)
         sub%cnodes(icnode)%global_cnode_number  = igcnode
         ! number of nodes where it maps from
         sub%cnodes(icnode)%nnod = 1

         ! number of variables it maps from 
         nvar = sub%nndf(indnode)
         sub%cnodes(icnode)%nvar = nvar
         ! number of nonzeros it creates in matrix C
         sub%cnodes(icnode)%nnz = nvar

         ! fill coordinates
         lxyz = sub%ndim
         sub%cnodes(icnode)%lxyz = lxyz
         allocate(sub%cnodes(icnode)%xyz(lxyz))
         sub%cnodes(icnode)%xyz = sub%xyz(indnode,:)

         ! fill coarse node nodes
         allocate(sub%cnodes(icnode)%insin(1))
         sub%cnodes(icnode)%insin(1) = indinode

         ! fill coarse node variables
         allocate(sub%cnodes(icnode)%ivsivn(nvar))
         do i = 1,nvar
            sub%cnodes(icnode)%ivsivn(i) = (indinode-1)*nvar + i
         end do

      end do

      ! copy globs behind corners - DO NOT nullify the icnode counter to place them behind corners
      do iglob = 1,nglob
         icnode = icnode + 1

         ! type of coarse node - edge or face
         sub%cnodes(icnode)%itype = sub%glob_type(iglob)

         ! global number
         igcnode = sub%global_glob_number(iglob)
         sub%cnodes(icnode)%global_cnode_number = igcnode

         ! number of nodes where it maps from
         nnodgl = sub%nglobnodes(iglob)
         sub%cnodes(icnode)%nnod = nnodgl
         ! number of variables it maps from 
         nvar = sub%nglobvar(iglob)
         sub%cnodes(icnode)%nvar = nvar

         ! fill coordinates
         lxyz = sub%ndim
         sub%cnodes(icnode)%lxyz = lxyz
         allocate(sub%cnodes(icnode)%xyz(lxyz))
         sub%cnodes(icnode)%xyz(lxyz) = 0._kr
         ! create averaged coordinates
         do i = 1,nnodgl
            indin = sub%ignsin(iglob,i)
            indn  = sub%iin(indin)
            sub%cnodes(icnode)%xyz(:) = sub%cnodes(icnode)%xyz(:) + sub%xyz(indn,:)
         end do
         sub%cnodes(icnode)%xyz(:) = sub%cnodes(icnode)%xyz(:) / nnodgl

         ! fill coarse node nodes
         allocate(sub%cnodes(icnode)%insin(nnodgl))
         do i = 1,nnodgl
            sub%cnodes(icnode)%insin(i) = sub%ignsin(iglob,i)
         end do

         ! fill coarse node variables
         allocate(sub%cnodes(icnode)%ivsivn(nvar))
         do i = 1,nvar
            sub%cnodes(icnode)%ivsivn(i) = sub%igvsivn(iglob,i)
         end do
      end do

      sub%is_cnodes_loaded = .true.

      return
end subroutine

!*******************************************************
subroutine dd_embed_cnodes(sub,nndf_coarse,lnndf_coarse)
!*******************************************************
! Merging corners with globs - order first corners, then globs
      use module_utils
      implicit none

! Subdomain structure
      type(subdomain_type),intent(inout) :: sub
      ! coarse NNDF
      integer,intent(in) :: lnndf_coarse
      integer,intent(in) ::  nndf_coarse(lnndf_coarse)

      ! local vars
      integer ::             ncorner
      integer ::             nglob
      integer ::             ncnodes
      integer ::             icnode, igcnode
      integer ::             ncdof, indn, i, kcdof

      integer ::            lkdof_coarse
      integer,allocatable :: kdof_coarse(:)


      ! check the prerequisities
      if (.not.sub%is_cnodes_loaded) then
         write(*,*) 'DD_EMBED_CNODES: Coarse nodes not loaded for subdomain:', sub%isub
         call error_exit
      end if

      ! determine number of coarse nodes
      ncorner = sub%ncorner
      nglob = sub%nglob
      ncnodes = ncorner + nglob
      sub%ncnodes = ncnodes

      ! create array of global coarse dof KDOFC(ncorner) with addresses before first global dof
      lkdof_coarse = lnndf_coarse
      allocate(kdof_coarse(lkdof_coarse))
      kdof_coarse(1) = 0
      do indn = 2,lnndf_coarse
         kdof_coarse(indn) = kdof_coarse(indn-1) + nndf_coarse(indn-1)
      end do

      ! Update embedding arrays
      do icnode = 1,ncnodes
         ! number of coarse dof it contains
         igcnode = sub%cnodes(icnode)%global_cnode_number 
         ncdof = nndf_coarse(igcnode)
         if (allocated(sub%cnodes(icnode)%igcdof)) then
            deallocate(sub%cnodes(icnode)%igcdof)
         end if

         sub%cnodes(icnode)%ncdof = ncdof
         allocate(sub%cnodes(icnode)%igcdof(ncdof))
         ! fill coarse node dof
         kcdof = kdof_coarse(igcnode)
         do i = 1,ncdof
            sub%cnodes(icnode)%igcdof(i) = kcdof + i
         end do

      end do
      sub%is_cnodes_embedded = .true.

      deallocate(kdof_coarse)

      return
end subroutine

!***************************
subroutine dd_prepare_c(sub)
!***************************
! Subroutine for preparing matrix of constraints C in sparse triplet format
      use module_utils
      implicit none

! Subdomain structure
      type(subdomain_type),intent(inout) :: sub

      ! local vars
      integer ::             nnzc
      integer ::             lc
      integer,allocatable ::  i_c_sparse(:)
      integer,allocatable ::  j_c_sparse(:)
      real(kr),allocatable ::   c_sparse(:)

      integer ::             lindrowc
      integer,allocatable ::  indrowc(:)

      integer ::              lkdof
      real(kr),allocatable ::  kdof(:)

      integer :: nnod
      integer :: inod,& 
                 nconstr, icdof, icn, inzc, irowc, ivar, ncdof, &
                 ncnodes, nrowc, nvar, lmatrix2

      ! check the prerequisities
      if (.not.sub%is_cnodes_embedded) then
         write(*,*) 'DD_PREPARE_C: Coarse nodes not ready for subdomain:', sub%isub
         call error_exit
      end if

      ! find number of rows in C (constraints)
      nrowc = 0
      ! find number of nonzeros in C
      nnzc  = 0

      ncnodes = sub%ncnodes
      do icn = 1,ncnodes
         if (sub%cnodes(icn)%used) then
            nrowc = nrowc + sub%cnodes(icn)%ncdof
            nnzc  = nnzc  + sub%cnodes(icn)%nnz
         end if
      end do
      nconstr = nrowc


      ! Creation of field KDOF(NNOD) with addresses before first global
      ! dof of node
      nnod  = sub%nnod
      lkdof = nnod
      allocate(kdof(lkdof))
      kdof(1) = 0
      do inod = 2,nnod
         kdof(inod) = kdof(inod-1) + sub%nndf(inod-1)
      end do

      lc   = nnzc
      allocate (i_c_sparse(lc),j_c_sparse(lc),c_sparse(lc))

      lindrowc = nrowc
      allocate(indrowc(lindrowc))

      ! create constraints from corners and mapping from local coarse dof to
      ! global coarse dof
      irowc = 0
      inzc  = 0
      do icn = 1,ncnodes
         if (sub%cnodes(icn)%used) then

            if (sub%cnodes(icn)%itype .eq. 3) then
               ! for corners, do not read matrices but construct the sparse matrix directly
               nvar = sub%cnodes(icn)%nvar
               do ivar = 1,nvar
                  irowc = irowc + 1
                  inzc  = inzc  + 1

                  i_c_sparse(inzc) = irowc
                  j_c_sparse(inzc) = sub%cnodes(icn)%ivsivn(ivar)
                  c_sparse(inzc)   = 1._kr

                  indrowc(irowc) = sub%cnodes(icn)%igcdof(ivar)
               end do
            else if (((sub%cnodes(icn)%arithmetic .eqv. .true.) .and. (sub%cnodes(icn)%itype .eq. 2)) &
                     .or. &
                     ((sub%cnodes(icn)%arithmetic .eqv. .true.) .and. (sub%cnodes(icn)%adaptive .eqv. .false.) &
                      .and. (sub%cnodes(icn)%itype .eq. 1))) then
               ! use arithmetic averages for edges and faces if desired
               nvar  = sub%cnodes(icn)%nvar
               ncdof = sub%cnodes(icn)%ncdof
               do icdof = 1,ncdof
                  irowc = irowc + 1

                  do ivar = icdof,nvar,ncdof
                     inzc  = inzc  + 1

                     i_c_sparse(inzc) = irowc
                     j_c_sparse(inzc) = sub%cnodes(icn)%ivsivn(ivar)
                     c_sparse(inzc)   = 1._kr
                  end do

                  indrowc(irowc) = sub%cnodes(icn)%igcdof(icdof)
               end do
            else if ((sub%cnodes(icn)%itype .eq. 1) .and. &
                     (sub%cnodes(icn)%adaptive .eqv. .true.)) then
               ! use adaptive constraint on face

               ! copy the matrix of constraints on glob into sparse triplet of subdomain matrix C
               ! row by row
               nvar      = sub%cnodes(icn)%nvar
               lmatrix2  = sub%cnodes(icn)%lmatrix2
               if (nvar.ne.lmatrix2) then
                  write(*,*) 'DD_PREPARE_C: Matrix dimension does not match for subdomain', sub%isub
                  call error_exit
               end if

               ncdof = sub%cnodes(icn)%ncdof
               do icdof = 1,ncdof
                  irowc = irowc + 1

                  do ivar = 1,nvar
                     inzc  = inzc  + 1

                     i_c_sparse(inzc) = irowc
                     j_c_sparse(inzc) = sub%cnodes(icn)%ivsivn(ivar)
                     c_sparse(inzc)   = sub%cnodes(icn)%matrix(icdof,ivar)
                  end do

                  indrowc(irowc) = sub%cnodes(icn)%igcdof(icdof)
               end do
            else
               continue
            end if

            ! check matrix bounds
            if (inzc .gt. lc) then
               write(*,*) 'DD_PREPARE_C: Too many entries in matrix C for subdomain', sub%isub
               call error_exit
            end if
         end if
      end do

      ! check matrix bounds
      if (inzc .ne. lc) then
         write(*,*) 'DD_PREPARE_C: Dimension of matrix C mismatch for subdomain', sub%isub,'inzc =',inzc,'lc =',lc
         call error_exit
      end if
      
      ! load the new matrix directly to the structure
      call dd_load_c(sub,nconstr,i_c_sparse, j_c_sparse, c_sparse, lc, nnzc, indrowc,lindrowc)

      deallocate (indrowc)
      deallocate (i_c_sparse,j_c_sparse,c_sparse)
      deallocate (kdof)

end subroutine

!**********************************************************************************************
subroutine dd_load_c(sub, nconstr,i_c_sparse, j_c_sparse, c_sparse, lc, nnzc, indrowc,lindrowc)
!**********************************************************************************************
! Subroutine for loading sparse matrix C into the SUB structure
      use module_utils
      implicit none

! Subdomain structure
      type(subdomain_type),intent(inout) :: sub

      ! number of constraints
      integer,intent(in) :: nconstr
      ! matrix C in sparse format
      integer,intent(in) ::             nnzc
      integer,intent(in) ::             lc
      integer,intent(in) ::  i_c_sparse(lc)
      integer,intent(in) ::  j_c_sparse(lc)
      real(kr),intent(in) ::   c_sparse(lc)

      integer,intent(in) ::  lindrowc
      integer,intent(in) ::   indrowc(lindrowc)

      ! local variables
      integer :: i

      ! check the prerequisities
      if (sub%is_c_loaded) then
         ! replace the loaded C
         if (debug) then
            write(*,*) 'DD_LOAD_C: Subdomain', sub%isub,': Matrix C already allocated. Rewriting.'
         end if
         deallocate(sub%i_c_sparse,sub%j_c_sparse,sub%c_sparse)
         deallocate(sub%indrowc)
         sub%is_c_loaded = .false.
      end if

      ! load C now
      sub%nconstr = nconstr
      sub%nnzc = nnzc
      sub%lc   = lc
      allocate(sub%i_c_sparse(lc))
      do i = 1,lc
         sub%i_c_sparse(i) = i_c_sparse(i)
      end do
      allocate(sub%j_c_sparse(lc))
      do i = 1,lc
         sub%j_c_sparse(i) = j_c_sparse(i)
      end do
      allocate(sub%c_sparse(lc))
      do i = 1,lc
         sub%c_sparse(i)   = c_sparse(i)
      end do

      ! load indrowc
      sub%lindrowc = lindrowc
      allocate(sub%indrowc(lindrowc))
      do i = 1,lindrowc
         sub%indrowc(i)   = indrowc(i)
      end do
         
      sub%is_c_loaded = .true.

end subroutine

!***************************************
subroutine dd_prepare_aug(sub,comm_self)
!***************************************
! Subroutine for preparing augmented matrix for computing with BDDC preconditioner
! |A C^T|
! |C  0 |
! and its factorization
      use module_mumps
      use module_utils
      implicit none

! Subdomain structure
      type(subdomain_type),intent(inout) :: sub
      ! communicator
      integer,intent(in) :: comm_self

      ! local vars
      integer ::  nnzc, nnza, ndof, nconstr
      integer ::  nnzaaug, laaug, ndofaaug
      integer ::  i, iaaug
      integer ::  mumpsinfo, aaugmatrixtype
      integer ::  icol, icoli
      logical :: parallel_analysis

      ! check the prerequisities
      if (.not.sub%is_matrix_loaded) then
         write(*,*) 'DD_PREPARE_AUG: Matrix is not loaded for subdomain:', sub%isub
         call error_exit
      end if
      if (.not.sub%is_c_loaded) then
         write(*,*) 'DD_PREPARE_AUG: Matrix of constraints C not loaded for subdomain:', sub%isub
         call error_exit
      end if

      ! if augmented system is already allocated, clear it - this can happen for adaptivity
      if (sub%is_aug_factorized) then

         deallocate(sub%i_aaug_sparse,sub%j_aaug_sparse,sub%aaug_sparse)
         sub%nnzaaug = 0
         sub%laaug   = 0

         call mumps_finalize(sub%mumps_aug) 
         sub%is_mumps_aug_active = .false.
         sub%is_aug_factorized = .false.
      end if

      ! join the new matrix directly to the structure as
      ! in the unsymmetric case:
      ! A C^T
      ! C  0   
      ! in the symmetric case :
      ! \A C^T
      !     0   

      ndof     = sub%ndof
      nconstr  = sub%nconstr
      ndofaaug = ndof + nconstr
      nnza     = sub%nnza
      nnzc     = sub%nnzc

      if      (sub%matrixtype .eq. 0) then
         ! unsymmetric case:
         nnzaaug = nnza + 2*nnzc
      else if (sub%matrixtype .eq. 1 .or. sub%matrixtype .eq. 2) then
         ! symmetric case:
         nnzaaug = nnza + nnzc
      end if

      laaug   = nnzaaug
      allocate(sub%i_aaug_sparse(laaug),sub%j_aaug_sparse(laaug),sub%aaug_sparse(laaug))

      ! copy entries of A
      iaaug = 0
      do i = 1,nnza
         iaaug = iaaug + 1
         sub%i_aaug_sparse(iaaug) = sub%i_a_sparse(i) 
         sub%j_aaug_sparse(iaaug) = sub%j_a_sparse(i) 
         sub%aaug_sparse(iaaug)   = sub%a_sparse(i)   
      end do
      ! copy entries of right block of C^T with proper shift in columns
      do i = 1,nnzc
         iaaug = iaaug + 1
         icoli = sub%j_c_sparse(i)
         icol  = sub%iivsvn(icoli)
         sub%i_aaug_sparse(iaaug) = icol
         sub%j_aaug_sparse(iaaug) = sub%i_c_sparse(i) + ndof
         sub%aaug_sparse(iaaug)   = sub%c_sparse(i)   
      end do
      if      (sub%matrixtype .eq. 0) then
         ! unsymmetric case: apply lower block of C
         do i = 1,nnzc
            iaaug = iaaug + 1
            icoli = sub%j_c_sparse(i)
            icol  = sub%iivsvn(icoli)
            sub%i_aaug_sparse(iaaug) = sub%i_c_sparse(i) + ndof
            sub%j_aaug_sparse(iaaug) = icol
            sub%aaug_sparse(iaaug)   = sub%c_sparse(i)   
         end do
      end if
      if (iaaug.ne.nnzaaug) then
         write(*,*) 'DD_PREPARE_AUG: Actual length of augmented matrix does not match for subdomain:', sub%isub
         call error_exit
      end if
      sub%nnzaaug = nnzaaug
      sub%laaug   = laaug


      ! factorize matrix Aaug
      ! Set type of matrix
      if      (sub%matrixtype .eq. 0) then
         ! unsymmetric case:
         aaugmatrixtype = 0
      else if (sub%matrixtype .eq. 1 .or. sub%matrixtype .eq. 2) then
         ! in symmetric case, saddle point problem makes the augmented matrix indefinite,
         ! even if the original matrix is SPD:
         aaugmatrixtype = 2
      end if
      call mumps_init(sub%mumps_aug,comm_self,aaugmatrixtype)

      ! Verbosity level of MUMPS
      if (debug) then
         mumpsinfo = 2
      else
         mumpsinfo = 0
      end if
      call mumps_set_info(sub%mumps_aug,mumpsinfo)

      ! Load matrix to MUMPS
      ndof     = sub%ndof
      nconstr  = sub%nconstr
      ndofaaug = ndof + nconstr

      nnzaaug = sub%nnzaaug
      laaug   = sub%laaug
      call mumps_load_triplet(sub%mumps_aug,ndofaaug,nnzaaug,&
                              sub%i_aaug_sparse,sub%j_aaug_sparse,sub%aaug_sparse,laaug)
      ! Analyze matrix
      parallel_analysis = .false.
      call mumps_analyze(sub%mumps_aug,parallel_analysis) 
      ! Factorize matrix 
      call mumps_factorize(sub%mumps_aug) 

      sub%is_mumps_aug_active = .true.
      sub%is_aug_factorized = .true.

end subroutine

!********************************************
subroutine dd_prepare_coarse(sub,keep_global)
!********************************************
! Subroutine for solving of system 
! | A C^T|| phis | = | 0 |
! | C  0 ||lambda|   | I |
! phis are coarse space basis functions on subdomain
! Then the routine builds the local coarse matrix:
! Ac = phis^T * A * phis 

      use module_mumps
      use module_utils
      implicit none

! Subdomain structure
      type(subdomain_type),intent(inout) :: sub
      ! store global vector of PHIS instead of the one restricted to interface PHISI
      logical,intent(in) :: keep_global

      ! local vars
      integer ::  ndof, nconstr, ndofaaug, ndofi, matrixtype
      integer ::  ndofc
      integer ::  i, j, indphis, nrhs, &
                  indphisstart, indi, icoarsem, lcoarsem
      integer ::  lphisi1, lphisi2, lphis1, lphis2

      integer ::             lphis
      real(kr),allocatable :: phis(:)

      integer ::             lac1, lac2
      real(kr),allocatable :: ac(:,:)

      ! check the prerequisities
      if (.not.sub%is_aug_factorized) then
         write(*,*) 'DD_PREPARE_COARSE: Augmented matrix in not factorized for subdomain:', sub%isub
         call error_exit
      end if
      if (.not.sub%is_mumps_aug_active) then
         write(*,*) 'DD_PREPARE_COARSE: Augmented matrix solver in not ready for subdomain:', sub%isub
         call error_exit
      end if

      ! if coarse problem is ready, clear it and prepare it again - this can happen for adaptivity
      if (sub%is_coarse_prepared) then
         write(*,*) 'Reinitializing coarse problem for subdomain ', sub%isub

         if (allocated(sub%phis)) then
            deallocate(sub%phis)
            sub%lphis1 = 0
            sub%lphis2 = 0
            sub%is_phis_prepared   = .false.
         end if

         if (allocated(sub%phisi)) then
            deallocate(sub%phisi)
            sub%lphisi1 = 0
            sub%lphisi2 = 0
            sub%is_phisi_prepared   = .false.
         end if

         deallocate(sub%coarsem)
         sub%lcoarsem = 0
         sub%ndofc = 0
         sub%is_coarse_prepared = .false.
      end if

      ndof     = sub%ndof
      nconstr  = sub%nconstr
      ndofaaug = ndof + nconstr

      ! Prepare phis with multiple RHS as dense matrix - stored in an linear array for MUMPS
      lphis = ndofaaug*nconstr
      allocate(phis(lphis))
      ! zero all entries
      call zero(phis,lphis)

      ! put identity into the block of constraints
      do j = 1,nconstr
         indphis = (j-1)*ndofaaug + ndof + j
         phis(indphis) = 1._kr
      end do

      ! solve the system with multiple RHS
      nrhs = nconstr
      call dd_solve_aug(sub, phis,lphis, nrhs) 

      !debug
      !write(*,*) 'Subdomain ',sub%isub,' coarse basis functions phis:'
      !do i = 1,ndofaaug
      !   write(*,'(100f13.6)') (phis((j-1)*ndofaaug + i),j = 1,nconstr)
      !end do

      ! Build subdomain coarse matrix by the fact that phis^T*A*phis = -lambda 
      lac1 = nconstr
      lac2 = nconstr
      allocate(ac(lac1,lac2))
      do j = 1,nconstr
         indphisstart  = (j-1)*ndofaaug + ndof
         do i = 1,nconstr
            ac(i,j) = -phis(indphisstart + i)
         end do
      end do
      !debug
      !write(*,*) 'Subdomain ',sub%isub,' coarse matrix:'
      !do i = 1,nconstr
      !   write(*,'(100f13.6)') (ac(i,j),j = 1,nconstr)
      !end do

      if (keep_global) then
         ndof   = sub%ndof
         lphis1 = ndof
         lphis2 = nconstr
         allocate(sub%phis(lphis1,lphis2))
         sub%lphis1 = lphis1
         sub%lphis2 = lphis2
         do i = 1,ndof
            do j = 1,nconstr
               indphis = (j-1)*ndofaaug + i
               sub%phis(i,j) = phis(indphis)
            end do
         end do
         sub%is_phis_prepared   = .true.
      else
         ! restrict vector phis to interface unknowns and load it to the structure
         ndofi   = sub%ndofi
         lphisi1 = ndofi
         lphisi2 = nconstr
         allocate(sub%phisi(lphisi1,lphisi2))
         sub%lphisi1 = lphisi1
         sub%lphisi2 = lphisi2
         do i = 1,ndofi
            indi = sub%iivsvn(i)
            do j = 1,nconstr
               indphis = (j-1)*ndofaaug + indi

               sub%phisi(i,j) = phis(indphis)
            end do
         end do
         sub%is_phisi_prepared   = .true.
      end if

      ! load the coarse matrix to the structure in appropriate format
      matrixtype = sub%matrixtype
      if      (matrixtype.eq.0) then
         ! in unsymmetric case, load the whole coarse matrix columnwise
         lcoarsem = nconstr * nconstr
         allocate(sub%coarsem(lcoarsem))
         icoarsem = 0
         do j = 1,nconstr
            do i = 1,nconstr
               icoarsem = icoarsem + 1
               sub%coarsem(icoarsem) = ac(i,j)
            end do
         end do
         sub%lcoarsem = lcoarsem
      else if (matrixtype.eq.1 .or. matrixtype.eq.2) then
         ! in symmetric case, load the upper triangle columnwise
         lcoarsem = (nconstr+1)*nconstr/2
         allocate(sub%coarsem(lcoarsem))
         icoarsem = 0
         do j = 1,nconstr
            do i = 1,j
               icoarsem = icoarsem + 1
               sub%coarsem(icoarsem) = ac(i,j)
            end do
         end do
         sub%lcoarsem = lcoarsem
      end if
      ! check the counter
      if (icoarsem.ne.lcoarsem) then
         write(*,*) 'DD_PREPARE_COARSE: Check of coarse matrix length failed for subdomain: ',sub%isub
         call error_exit
      end if

      ! number of coarse degrees of freedom equals number of constraints in this implementation
      ndofc = nconstr
      sub%ndofc = ndofc

      sub%is_coarse_prepared = .true.

      deallocate(ac)
      deallocate(phis)
end subroutine

!****************************************
subroutine dd_get_aug_size(sub, ndofaaug)
!****************************************
! Subroutine for getting size of the system 
! | A C^T|| z_i | = | res_i |
! | C  0 || mu  |   |   0   |

      use module_utils
      implicit none

! Subdomain structure
      type(subdomain_type),intent(in) :: sub

      ! augmented problem size 
      integer,intent(out) :: ndofaaug

      ! local vars
      integer :: ndof, nconstr

      ! check the prerequisities
      if (.not.sub%is_c_loaded) then
         write(*,*) 'DD_GET_AUG_SIZE: C is not loaded for subdomain', sub%isub
         call error_exit
      end if
      if (.not.sub%is_mesh_loaded) then
         write(*,*) 'DD_GET_AUG_SIZE: Mesh is not loaded for subdomain', sub%isub
         call error_exit
      end if

      ndof     = sub%ndof
      nconstr  = sub%nconstr

      ndofaaug = ndof + nconstr

end subroutine

!*******************************************
subroutine dd_solve_aug(sub, vec,lvec, nrhs)
!*******************************************
! Subroutine for solving system 
! | A C^T|| z_i | = | res_i |
! | C  0 || mu  |   |   0   |
! on entry, vec contains res_i
! on exit, vec contains z_i

      use module_mumps
      use module_utils
      implicit none

! Subdomain structure
      type(subdomain_type),intent(inout) :: sub

      integer,intent(in) ::    lvec
      real(kr),intent(inout) :: vec(lvec)

      ! how many right hand sides are hidden in vec
      integer,intent(in) ::    nrhs

      ! local vars
      integer ::  ndofaaug

      ! check the prerequisities
      if (.not.sub%is_aug_factorized) then
         write(*,*) 'DD_SOLVE_AUG: Augmented matrix in not factorized for subdomain:', sub%isub
         call error_exit
      end if
      if (.not.sub%is_mumps_aug_active) then
         write(*,*) 'DD_SOLVE_AUG: Augmented matrix solver in not ready for subdomain:',sub%isub
         call error_exit
      end if
      if (mod(lvec,nrhs) .ne. 0) then
         write(*,*) 'DD_SOLVE_AUG: Unclear what the augmented size is:', sub%isub
         call error_exit
      end if

      call dd_get_aug_size(sub, ndofaaug)
      if (ndofaaug.ne.lvec/nrhs) then
         write(*,*) 'DD_SOLVE_AUG: Length of augmented system does not match:', sub%isub
         call error_exit
      end if

      ! solve the system with multiple RHS
      call mumps_resolve(sub%mumps_aug,vec,lvec,nrhs)

end subroutine

!*****************************************************************
subroutine dd_phisi_apply(sub, transposed, vec1,lvec1, vec2,lvec2)
!*****************************************************************
! Subroutine for multiplication of vector VEC1 by PHISI matrix
! vec2 = phisi * vec1 + vec2
! phisi are coarse space basis functions on subdomain restricted to interface

      use module_utils
      implicit none

! Subdomain structure
      type(subdomain_type),intent(in) :: sub

      ! is matrix transposed
      logical,intent(in) :: transposed

      ! input vector
      integer,intent(in) :: lvec1
      real(kr),intent(in) :: vec1(lvec1)
      ! output vector
      integer,intent(in) ::  lvec2
      real(kr),intent(out) :: vec2(lvec2)

      ! local vars
      logical :: wrong_dim

      ! BLAS vars
      character(1) :: TRANS
      integer :: M, N, LDA, INCX, INCY
      real(kr) :: alpha, beta

      ! check the prerequisities
      if (.not.sub%is_phisi_prepared) then
         write(*,*) 'DD_PHISI_APPLY: PHISI matrix not ready:', sub%isub
         call error_exit
      end if
      ! check dimensions
      wrong_dim = .false.
      if (transposed) then
         if (lvec1 .ne. sub%lphisi1 .or. lvec2 .ne. sub%lphisi2) then 
            wrong_dim = .true.
         end if
      else
         if (lvec1 .ne. sub%lphisi2 .or. lvec2 .ne. sub%lphisi1) then 
            wrong_dim = .true.
         end if
      end if
      if (wrong_dim) then
         write(*,*) 'DD_PHISI_APPLY: Dimensions mismatch:', sub%isub
         call error_exit
      end if

      ! checking done, perform multiply by BLAS
      if (transposed) then
         TRANS = 'T'
      else
         TRANS = 'N'
      end if
      M = sub%lphisi1
      N = sub%lphisi2
      ALPHA = 1._kr
      LDA = max(1,M)
      INCX = 1
      BETA = 1._kr ! sum second vector
      INCY = 1
      if (kr.eq.8) then
         ! double precision
         call DGEMV(TRANS,M,N,ALPHA,sub%phisi,LDA,vec1,INCX,BETA,vec2,INCY)
      else if (kr.eq.4) then
         ! single precision
         call SGEMV(TRANS,M,N,ALPHA,sub%phisi,LDA,vec1,INCX,BETA,vec2,INCY)
      end if

end subroutine

!****************************************************************
subroutine dd_phis_apply(sub, transposed, vec1,lvec1, vec2,lvec2)
!****************************************************************
! Subroutine for multiplication of vector VEC1 by PHIS matrix
! vec2 = phis * vec1 + vec2
! phis are coarse space basis functions on subdomain

      use module_utils
      implicit none

! Subdomain structure
      type(subdomain_type),intent(in) :: sub

      ! is matrix transposed
      logical,intent(in) :: transposed

      ! input vector
      integer,intent(in) :: lvec1
      real(kr),intent(in) :: vec1(lvec1)
      ! output vector
      integer,intent(in) ::  lvec2
      real(kr),intent(out) :: vec2(lvec2)

      ! local vars
      logical :: wrong_dim

      ! BLAS vars
      character(1) :: TRANS
      integer :: M, N, LDA, INCX, INCY
      real(kr) :: alpha, beta

      ! check the prerequisities
      if (.not.sub%is_phis_prepared) then
         write(*,*) 'DD_PHIS_APPLY: PHIS matrix not ready:', sub%isub
         call error_exit
      end if
      ! check dimensions
      wrong_dim = .false.
      if (transposed) then
         if (lvec1 .ne. sub%lphis1 .or. lvec2 .ne. sub%lphis2) then 
            wrong_dim = .true.
         end if
      else
         if (lvec1 .ne. sub%lphis2 .or. lvec2 .ne. sub%lphis1) then 
            wrong_dim = .true.
         end if
      end if
      if (wrong_dim) then
         write(*,*) 'DD_PHIS_APPLY: Dimensions mismatch:', sub%isub
         call error_exit
      end if

      ! checking done, perform multiply by BLAS
      if (transposed) then
         TRANS = 'T'
      else
         TRANS = 'N'
      end if
      M = sub%lphis1
      N = sub%lphis2
      ALPHA = 1._kr
      LDA = max(1,M)
      INCX = 1
      BETA = 1._kr ! sum second vector
      INCY = 1
      if (kr.eq.8) then
         ! double precision
         call DGEMV(TRANS,M,N,ALPHA,sub%phis,LDA,vec1,INCX,BETA,vec2,INCY)
      else if (kr.eq.4) then
         ! single precision
         call SGEMV(TRANS,M,N,ALPHA,sub%phis,LDA,vec1,INCX,BETA,vec2,INCY)
      end if

end subroutine

!********************************************************************
subroutine dd_get_my_coarsem_length(suba,lsuba,indexsub,lindexsub,la)
!********************************************************************
! Subroutine for obtaining length of coarse matrix
      implicit none

! array of sub structure for actual subdomains
      integer,intent(in) ::             lsuba
      type(subdomain_type),intent(in) :: suba(lsuba)
      ! subdomain number
      integer,intent(in) :: lindexsub
      integer,intent(in) :: indexsub(lindexsub)
      ! coarse matrix length
      integer,intent(out) :: la

      ! local vars
      integer :: isub, isub_loc, lcoarsem

      ! reset counter
      la = 0
      do isub_loc = 1,lindexsub
         isub = indexsub(isub_loc)

         lcoarsem = suba(isub_loc)%lcoarsem

         la = la + lcoarsem
      end do

end subroutine
      
!*****************************************************************************************************
subroutine dd_get_my_coarsem(suba,lsuba,matrixtype,indexsub,lindexsub,i_sparse, j_sparse, a_sparse,la)
!*****************************************************************************************************
! Subroutine for obtaining length of coarse matrix
      use module_utils
      implicit none

! array of sub structure for actual subdomains
      integer,intent(in) ::             lsuba
      type(subdomain_type),intent(in) :: suba(lsuba)
      ! matrixtype 
      integer,intent(in) :: matrixtype
      ! subdomain number
      integer,intent(in) :: lindexsub
      integer,intent(in) :: indexsub(lindexsub)
      ! coarse matrix length
      integer,intent(in) :: la
      integer,intent(out) ::  i_sparse(la)
      integer,intent(out) ::  j_sparse(la)
      real(kr),intent(out) :: a_sparse(la)

      ! local vars
      integer :: isub_loc, isub, lcoarsem, lindrowc
      integer :: irow, jcol, ielm, ia, i, j


      ia = 0
      do isub_loc = 1,lindexsub
         isub = indexsub(isub_loc)

         lindrowc = suba(isub_loc)%lindrowc
         lcoarsem = suba(isub_loc)%lcoarsem

         if (matrixtype.eq.0) then
            ! nonsymmetric matrix
            if (lcoarsem .ne. lindrowc**2) then
               write(*,*) 'DD_GET_MY_COARSEM: Error in coarse matrix length for nonsymmetric matrix.'
               call error_exit
            end if

            ielm = 0
            do j = 1,lindrowc
               jcol = suba(isub_loc)%indrowc(j)
               do i = 1,lindrowc
                  irow = suba(isub_loc)%indrowc(i)

                  ielm = ielm + 1
                  ia   = ia + 1

                  i_sparse(ia) = irow
                  j_sparse(ia) = jcol
                  a_sparse(ia) = suba(isub_loc)%coarsem(ielm)
               end do
            end do
         else if (matrixtype.eq.1 .or. matrixtype.eq.2) then
            ! symmetric matrix
            if (lcoarsem .ne. ((lindrowc+1) * lindrowc)/2) then
               write(*,*) 'DD_GET_MY_COARSEM: Error in coarse matrix length for symmetric matrix.'
               call error_exit
            end if

            ielm = 0
            do j = 1,lindrowc
               jcol = suba(isub_loc)%indrowc(j)
               do i = 1,j
                  irow = suba(isub_loc)%indrowc(i)

                  ielm = ielm + 1
                  ia   = ia + 1

                  i_sparse(ia) = irow
                  j_sparse(ia) = jcol
                  a_sparse(ia) = suba(isub_loc)%coarsem(ielm)
               end do
            end do
         else
            write(*,*) 'DD_GET_MY_COARSEM: Unknown matrixtype:',matrixtype
         end if
      end do
end subroutine
      

!*********************************************
subroutine dd_multiply_by_schur(sub,x,lx,y,ly)
!*********************************************
! Subroutine for multiplication of interface vector by Schur complement
      use module_utils
      use module_mumps
      use module_sm
      implicit none

! Subdomain structure
      type(subdomain_type),intent(inout) :: sub

      ! input vector
      integer,intent(in)  :: lx
      real(kr),intent(in) ::  x(lx)

      ! output vector
      integer,intent(in)   :: ly
      real(kr),intent(out) ::  y(ly)

      ! local vars
      integer ::              laux1
      real(kr),allocatable ::  aux1(:)
      integer ::              laux2
      real(kr),allocatable ::  aux2(:)

      integer :: ndofi, ndofo, nnza12, la12, nnza21, la21, nnza22, la22, &
                 matrixtype_aux, matrixtype
      integer :: i
      logical :: is_symmetric_storage

      ! check the prerequisities
      if (.not. (sub%is_interior_factorized)) then
         write(*,*) 'DD_PREPARE_SCHUR: Interior block not factorized yet.',sub%isub
         call error_exit
      end if
      if (.not. (sub%ndofi .eq. lx .or. .not. lx .eq. ly)) then
         write(*,*) 'DD_PREPARE_SCHUR: Inconsistent data size.'
         call error_exit
      end if
 
      ! prepare rhs vector for backsubstitution to problem A_11*aux1 = -A_12*x
      ndofo = sub%ndofo
      laux1 = ndofo
      allocate(aux1(laux1))

      ! prepare rhs vector for backsubstitution to problem A_11*aux1 = -A_12*x
      ! with offdiagonal blocks, use as nonsymmetric
      matrixtype_aux = 0
      nnza12     = sub%nnza12
      la12       = sub%la12
      call sm_vec_mult(matrixtype_aux, nnza12, &
                       sub%i_a12_sparse, sub%j_a12_sparse, sub%a12_sparse, la12, &
                       x,lx, aux1,laux1)

      ! resolve interior problem by MUMPS
      call mumps_resolve(sub%mumps_interior_block,aux1,laux1)

      if (sub%istorage .eq. 4) then
         is_symmetric_storage = .true.
      else
         is_symmetric_storage = .false.
      end if

      ! prepare auxiliary vector for multiplication
      ndofi = sub%ndofi
      laux2 = ndofi
      allocate(aux2(laux2))

      ! get aux2 = A_21*aux1, i.e. aux2 = A_21 * (A_11)^-1 * A_12 * x
      if (is_symmetric_storage) then
         matrixtype_aux = 0
         nnza12     = sub%nnza12
         la12       = sub%la12
         ! use the matrix with transposed indices in the call sm_vec_mult
         call sm_vec_mult(matrixtype_aux, nnza12, &
                          sub%j_a12_sparse, sub%i_a12_sparse, sub%a12_sparse, la12, &
                          aux1,laux1, aux2,laux2)
      else
         matrixtype_aux = 0
         nnza21     = sub%nnza21
         la21       = sub%la21
         call sm_vec_mult(matrixtype_aux, nnza21, &
                          sub%i_a21_sparse, sub%j_a21_sparse, sub%a21_sparse, la21, &
                          aux1,laux1, aux2,laux2)
      end if

      ! get y = A_22*x
      matrixtype = sub%matrixtype
      nnza22     = sub%nnza22
      la22       = sub%la22
      call sm_vec_mult(matrixtype, nnza22, &
                       sub%i_a22_sparse, sub%j_a22_sparse, sub%a22_sparse, la22, &
                       x,lx, y,ly)

      ! add results together to get y = y - aux2, i.e. y = A_22 * x - A_21 * (A_11)^-1 * A_12 * x, or y = (A_22 - A_21 * (A_11)^-1 * A_12) * x
      do i = 1,ly
         y(i) = y(i) - aux2(i)
      end do
 
      deallocate(aux1)
      deallocate(aux2)

end subroutine

!********************************************
subroutine dd_resolve_interior(sub,x,lx,y,ly)
!********************************************
! Subroutine for resolution of interior variables with interface fixed (aka solving Dirichlet problem)
      use module_mumps
      use module_sm
      use module_utils
      implicit none

! Subdomain structure
      type(subdomain_type),intent(inout) :: sub

      ! input vector
      integer,intent(in)  :: lx
      real(kr),intent(in) ::  x(lx)

      ! output vector
      integer,intent(in)   :: ly
      real(kr),intent(out) ::  y(ly)

      ! local vars
      integer ::   matrixtype_aux
      integer ::  nnza12, la12

      integer :: ndofo, ndof
      integer ::              laux1
      real(kr), allocatable :: aux1(:)
      integer ::              laux2
      real(kr), allocatable :: aux2(:)
      integer ::              laux3
      real(kr), allocatable :: aux3(:)

      integer :: io, ind, i

      ! check the prerequisities
      if (.not. (sub%is_blocked)) then
         write(*,*) 'DD_RESOLVE_INTERIOR: Matrix is not in blocked format. Call routine to do this.'
         call error_exit
      end if
      if (.not. (sub%is_interior_factorized)) then
         write(*,*) 'DD_RESOLVE_INTERIOR: Interior block not factorized yet.'
         call error_exit
      end if
      if (.not. (sub%ndofi .eq. lx)) then
         write(*,*) 'DD_RESOLVE_INTERIOR: Inconsistent data size on input.'
         call error_exit
      end if
      if (.not. (sub%ndof .eq. ly)) then
         write(*,*) 'DD_RESOLVE_INTERIOR: Inconsistent data size on output.'
         call error_exit
      end if

      ndofo = sub%ndofo
      laux1 = ndofo
      allocate(aux1(laux1))
 
      ndof = sub%ndof
      laux2 = ndof
      allocate(aux2(laux2))
      ! copy right hand side into aux2
      do i = 1,ndof
         aux2(i) = sub%rhs(i)
      end do
      ! fix BC in aux2
      if (sub%is_bc_present) then
         call sm_prepare_rhs(sub%ifix,sub%lifix,sub%bc,sub%lbc,aux2,laux2)
      end if

      laux3 = ndofo
      allocate(aux3(laux3))

      ! map aux2 (length of whole subdomain solution) to subdomain interior aux1
      call dd_map_sub_to_subo(sub, aux2,laux2, aux3,laux3)

      deallocate(aux2)

      ! prepare rhs vector for backsubstitution to problem A_11*aux1 = -A_12*x
      ! with offdiagonal blocks, use as nonsymmetric
      matrixtype_aux = 0
      nnza12     = sub%nnza12
      la12       = sub%la12
      call sm_vec_mult(matrixtype_aux, nnza12, &
                       sub%i_a12_sparse, sub%j_a12_sparse, sub%a12_sparse, la12, &
                       x,lx, aux1,laux1)

      ! prepare rhs into aux1
      ! aux = f - aux
      do io = 1,ndofo
         aux1(io) = aux3(io) - aux1(io)
      end do

      deallocate(aux3)

      ! resolve interior problem by MUMPS
      call mumps_resolve(sub%mumps_interior_block,aux1,laux1)

      ! embed vector of interior solution into subdomain solution y
      do io = 1,ndofo
         ind = sub%iovsvn(io)
         y(ind) = y(ind) + aux1(io)
      end do

      deallocate(aux1)

      ! embed interface solution x into y
      call dd_map_subi_to_sub(sub, x,lx, y,ly)
      
end subroutine

!*******************************************************************************************
subroutine dd_prepare_reduced_rhs(suba,lsuba,sub2proc,lsub2proc,indexsub,lindexsub,comm_all)
!*******************************************************************************************
! Subroutine for construction of reduced rhs
! g = f_2 - A_21 * A_11^-1 * f_1
      use module_mumps
      use module_sm
      use module_utils
      implicit none

! array of sub structure for actual subdomains
      integer,intent(in) ::                lsuba
      type(subdomain_type),intent(inout) :: suba(lsuba)
! division of subdomains to processors at previous level
      integer,intent(in) :: lsub2proc
      integer,intent(in) ::  sub2proc(lsub2proc)
      ! subdomain number of local subdomains
      integer,intent(in) :: lindexsub
      integer,intent(in) ::  indexsub(lindexsub)
      ! communicator
      integer,intent(in) :: comm_all

      ! local vars
      integer ::              laux1
      real(kr),allocatable ::  aux1(:)
      integer ::              laux2
      real(kr),allocatable ::  aux2(:)
      integer ::              laux3
      real(kr),allocatable ::  aux3(:)

      integer :: ndofi, ndofo, nnza12, la12, nnza21, la21, &
                 matrixtype_aux, lg, ndof
      integer :: i, isub_loc, isub
      logical :: is_symmetric_storage


      ! loop over subdomains
      do isub_loc = 1,lindexsub
         isub = indexsub(isub_loc)

         ! check the prerequisities
         if (.not. (suba(isub_loc)%is_interior_factorized)) then
            write(*,*) 'DD_PREPARE_REDUCED_RHS: Interior block not factorized yet.'
            call error_exit
         end if
         if (.not. (suba(isub_loc)%is_rhs_loaded)) then
            write(*,*) 'DD_PREPARE_REDUCED_RHS: RHS not loaded.'
            call error_exit
         end if
         if (.not. (suba(isub_loc)%is_weights_ready)) then
            write(*,*) 'DD_PREPARE_REDUCED_RHS: Weights not ready.'
            call error_exit
         end if
 
         ndof = suba(isub_loc)%ndof
         laux2 = ndof
         allocate(aux2(laux2))
         ! copy right hand side into aux2
         do i = 1,ndof
            aux2(i) = suba(isub_loc)%rhs(i)
         end do
         ! fix BC in aux2
         if (suba(isub_loc)%is_bc_present) then
            call sm_prepare_rhs(suba(isub_loc)%ifix,suba(isub_loc)%lifix,&
                                suba(isub_loc)%bc,suba(isub_loc)%lbc,aux2,laux2)
         end if

         ! prepare f_1
         ndofo = suba(isub_loc)%ndofo
         laux1 = ndofo
         allocate(aux1(laux1))

         ! map aux2 (length of whole subdomain solution) to subdomain interior aux1
         call dd_map_sub_to_subo(suba(isub_loc), aux2,laux2, aux1,laux1)

         deallocate(aux2)

         ! solve problem A_11*aux1 = f_1
         ! by MUMPS
         call mumps_resolve(suba(isub_loc)%mumps_interior_block,aux1,laux1)
         
         if (suba(isub_loc)%istorage .eq. 4) then
            is_symmetric_storage = .true.
         else
            is_symmetric_storage = .false.
         end if

         ! prepare auxiliary vector for multiplication
         ndofi = suba(isub_loc)%ndofi
         laux2 = ndofi
         allocate(aux2(laux2))

         ! get aux2 = A_21*aux1, i.e. aux2 = A_21 * (A_11)^-1 * f_1
         if (is_symmetric_storage) then
            matrixtype_aux = 0
            nnza12     = suba(isub_loc)%nnza12
            la12       = suba(isub_loc)%la12
            ! use the matrix with transposed indices in the call sm_vec_mult
            call sm_vec_mult(matrixtype_aux, nnza12, &
                             suba(isub_loc)%j_a12_sparse, suba(isub_loc)%i_a12_sparse, suba(isub_loc)%a12_sparse, la12,&
                             aux1,laux1, aux2,laux2)
         else
            matrixtype_aux = 0
            nnza21     = suba(isub_loc)%nnza21
            la21       = suba(isub_loc)%la21
            call sm_vec_mult(matrixtype_aux, nnza21, &
                             suba(isub_loc)%i_a21_sparse, suba(isub_loc)%j_a21_sparse, suba(isub_loc)%a21_sparse, la21, &
                             aux1,laux1, aux2,laux2)
         end if

         ! prepare f_2
         laux3 = ndofi
         allocate(aux3(laux3))
         ! copy proper part of RHS
         do i = 1,ndofi
            aux3(i) = suba(isub_loc)%rhs(suba(isub_loc)%iivsvn(i))
         end do
         ! weigh f_2
         do i = 1,ndofi
            aux3(i) = suba(isub_loc)%wi(i) * aux3(i)
         end do
         ! add eliminated BC to right hand side
         if (suba(isub_loc)%is_bc_nonzero) then
            do i = 1,ndofi
               aux3(i) = aux3(i) + suba(isub_loc)%bc(suba(isub_loc)%iivsvn(i))
            end do
         end if

         ! add results together to get aux3 = aux3 - aux2, i.e. aux3 = f_2 - A_21 * (A_11)^-1 * f_1
         do i = 1,ndofi
            aux3(i) = aux3(i) - aux2(i)
         end do

         ! load reduced RHS for communication structure
         call dd_comm_upload(suba(isub_loc), aux3,laux3) 
 
         ! save my own g into sub
         lg = ndofi
         suba(isub_loc)%lg = lg
         allocate(suba(isub_loc)%g(lg))
         do i = 1,lg
            suba(isub_loc)%g(i) = aux3(i)
         end do

         deallocate(aux1)
         deallocate(aux2)
         deallocate(aux3)
      end do

      ! communicate condensed right hand side
      call dd_comm_swapdata(suba,lsuba, indexsub,lindexsub, sub2proc,lsub2proc,comm_all)

      ! loop over subdomains
      do isub_loc = 1,lindexsub
         isub = indexsub(isub_loc)

         ndofi = suba(isub_loc)%ndofi
         laux3 = ndofi
         allocate(aux3(laux3))

         ! download contribution to g from my neighbours
         call dd_comm_download(suba(isub_loc), aux3,laux3) 

         ! sum the contributions
         lg = suba(isub_loc)%lg
         do i = 1,lg
            suba(isub_loc)%g(i) = suba(isub_loc)%g(i) + aux3(i)
         end do

         suba(isub_loc)%is_reduced_rhs_loaded = .true.
 
         deallocate(aux3)
      end do
end subroutine

!*******************************************
subroutine dd_get_reduced_rhs(sub, rhs,lrhs)
!*******************************************
! Subroutine for getting reduced rhs
      use module_utils
      implicit none

! Subdomain structure
      type(subdomain_type),intent(in) :: sub
      ! reduced rhs
      integer,intent(in) ::  lrhs
      real(kr),intent(out) :: rhs(lrhs)

      ! local vars
      integer :: i


      ! check the prerequisities
      if (.not.sub%is_reduced_rhs_loaded) then
         write(*,*) 'DD_GET_REDUCED_RHS: Reduced RHS is not loaded for subdomain:', sub%isub
         call error_exit
      end if
      ! check the size
      if (lrhs .ne. sub%lg) then
         write(*,*) 'DD_GET_REDUCED_RHS: RHS size mismatch for subdomain:', sub%isub
         call error_exit
      end if
 
      ! copy rhs
      do i = 1,lrhs
         rhs(i) = sub%g(i)
      end do

end subroutine

!************************************************
subroutine dd_get_interface_size(sub,ndofi,nnodi)
!************************************************
! Subroutine for finding size of subdomain data
      use module_utils
      implicit none
! Subdomain structure
      type(subdomain_type),intent(in) :: sub
! Length of vector of interface dof
      integer,intent(out) :: ndofi
! Length of vector of interface nodes
      integer,intent(out) :: nnodi

      ! check the prerequisities
      if (.not.sub%is_interface_loaded) then
         write(*,*) 'DD_GET_INTERFACE_SIZE: Interface not loaded for subdomain ',sub%isub
         call error_exit
      end if

      ! if all checks are OK, return subdomain interface size
      ndofi = sub%ndofi
      nnodi = sub%nnodi

end subroutine

!******************************************
subroutine dd_get_size(sub,ndof,nnod,nelem)
!******************************************
! Subroutine for finding size of subdomain data
      use module_utils
      implicit none
! Subdomain structure
      type(subdomain_type),intent(in) :: sub
! Number of subdomain dof
      integer,intent(out) :: ndof
! Number of subdomain nodes
      integer,intent(out) :: nnod
! Number of subdomain elements
      integer,intent(out) :: nelem

      if (.not.sub%is_mesh_loaded) then
         write(*,*) 'DD_GET_SIZE: Mesh not loaded for subdomain ',sub%isub
         call error_exit
      end if

      ! if all checks are OK, return subdomain size
      ndof  = sub%ndof
      nnod  = sub%nnod
      nelem = sub%nelem

end subroutine

!***********************************************
subroutine dd_get_coarse_size(sub,ndofc,ncnodes)
!***********************************************
! Subroutine for finding local size of subdomain coarse vector
      use module_utils
      implicit none
! Subdomain structure
      type(subdomain_type),intent(in) :: sub
! Length of vector of coarse dof on subdomain
      integer,intent(out) :: ndofc
! Length of vector of coarse nodes on subdomain
      integer,intent(out) :: ncnodes

      ! check prerequisites
      if (.not.sub%is_cnodes_loaded) then
         write(*,*) 'DD_GET_COARSE_SIZE: Coarse nodes not loaded for subdomain',sub%isub
         call error_exit
      end if

      ! if all checks are OK, return subdomain coarse size
      ndofc   = sub%ndofc
      ncnodes = sub%ncnodes

end subroutine

!****************************************************
subroutine dd_get_coarse_cnodes(sub,icnodes,licnodes)
!****************************************************
! Subroutine for getting subdomain coarse vector
      use module_utils
      implicit none
! Subdomain structure
      type(subdomain_type),intent(in) :: sub
! Global indices of coarse pseudo nodes
      integer,intent(in)  :: licnodes
      integer,intent(out) ::  icnodes(licnodes)

      ! local vars
      integer :: ncnodes, icn

      ! check prerequisites
      if (.not.sub%is_cnodes_loaded) then
         write(*,*) 'DD_GET_COARSE_CNODES: Coarse nodes not loaded for subdomain',sub%isub
         call error_exit
      end if

      ncnodes = sub%ncnodes

      ! check length of output array
      if (ncnodes .ne. licnodes) then
         write(*,*) 'DD_GET_COARSE_CNODES: requested array has wrong dimension'
         call error_exit
      end if

      ! if all checks are OK, return subdomain coarse nodes indices
      do icn = 1,ncnodes
         icnodes(icn) = sub%cnodes(icn)%global_cnode_number
      end do
end subroutine

!**********************************************
subroutine dd_get_number_of_crows(sub,lindrowc)
!**********************************************
! Subroutine for finding number of rows in C (constraints matrix) on subdomain
      use module_utils
      implicit none
! Subdomain structure
      type(subdomain_type),intent(in) :: sub
! Number of rows in C
      integer,intent(out) :: lindrowc

      if (.not.sub%is_c_loaded) then
         write(*,*) 'DD_GET_NUMBER_OF_CROWS: Coarse matrix not loaded for subdomain ',sub%isub
         call error_exit
      end if

      ! if all checks are OK, return number of rows in C
      lindrowc = sub%lindrowc

end subroutine

!*****************************************
subroutine dd_get_number_of_cnnz(sub,nnzc)
!*****************************************
! Subroutine for finding number of nonzeros in C (constraints matrix) on subdomain
      use module_utils
      implicit none
! Subdomain structure
      type(subdomain_type),intent(in) :: sub
! Number of nonzeros in C
      integer,intent(out) :: nnzc

      if (.not.sub%is_c_loaded) then
         write(*,*) 'DD_GET_NUMBER_OF_CNNZ: Matrix C not loaded yet for subdomain ',sub%isub
         call error_exit
      end if

      ! if all checks are OK, return number of nonzeros in C
      nnzc = sub%nnzc

end subroutine

!******************************************************
subroutine dd_get_subdomain_crows(sub,indrowc,lindrowc)
!******************************************************
! Subroutine for getting rows in C (constraints matrix) on subdomain
      use module_utils
      implicit none
! Subdomain structure
      type(subdomain_type),intent(in) :: sub
! Rows in C
      integer,intent(in) :: lindrowc
      integer,intent(out) :: indrowc(lindrowc)
! local vars
      integer :: i

      if (.not.sub%is_c_loaded) then
         write(*,*) 'DD_GET_SUBDOMAIN_CROWS: Matrix C not loaded yet for subdomain ',sub%isub
         call error_exit
      end if
      if (sub%lindrowc .ne. lindrowc) then
         write(*,*) 'DD_GET_SUBDOMAIN_CROWS: Size of array for rows of C not consistent.'
         write(*,*) 'lindrowc :', sub%lindrowc, 'array size: ',lindrowc
         call error_exit
      end if

      ! if all checks are OK, return rows of C
      do i = 1,lindrowc
         indrowc(i) = sub%indrowc(i)
      end do

end subroutine

!*****************************************************
subroutine dd_get_subdomain_corner_number(sub,ncorner)
!*****************************************************
! Subroutine for getting number of corners on subdomain
      use module_utils
      implicit none
! Subdomain structure
      type(subdomain_type),intent(in) :: sub
! Number of corners
      integer,intent(out) :: ncorner

      if (.not.sub%is_corners_loaded) then
         write(*,*) 'DD_GET_SUBDOMAIN_CONRER_NUMBER: Corners not loaded yet for subdomain ',sub%isub
         call error_exit
      end if

      ! if all checks are OK, return number of corners on subdomain
      ncorner = sub%ncorner

end subroutine

!**************************************************************************
subroutine dd_get_subdomain_c(sub,i_c_sparse,j_c_sparse,c_sparse,lc_sparse)
!**************************************************************************
! Subroutine for getting sparse matrix C
      use module_utils
      implicit none
! Subdomain structure
      type(subdomain_type),intent(in) :: sub
! Corners mapping
      integer,intent(in) ::  lc_sparse
      integer,intent(out) :: i_c_sparse(lc_sparse), j_c_sparse(lc_sparse)
      real(kr),intent(out) ::  c_sparse(lc_sparse)

! local vars
      integer :: nnzc, i

      if (.not.sub%is_c_loaded) then
         write(*,*) 'DD_GET_SUBDOMAIN_C: Matrix C is not loaded yet for subdomain ',sub%isub
         call error_exit
      end if
      if (sub%lc .ne. lc_sparse) then
         write(*,*) 'DD_GET_SUBDOMAIN_C: Size of array for matrix C mismatch.'
         write(*,*) 'lc :', sub%lc, 'array size: ',lc_sparse 
         call error_exit
      end if

      ! if all checks are OK, return subdomain C matrix
      nnzc = sub%nnzc
      do i = 1,nnzc
         i_c_sparse(i) = sub%i_c_sparse(i)
      end do
      do i = 1,nnzc
         j_c_sparse(i) = sub%j_c_sparse(i)
      end do
      do i = 1,nnzc
         c_sparse(i) = sub%c_sparse(i)
      end do

end subroutine

!***********************************************************
subroutine dd_get_subdomain_interface_nndf(sub,nndfi,lnndfi)
!***********************************************************
! Subroutine for getting subdomain interface nndfi array
      use module_utils
      implicit none
! Subdomain structure
      type(subdomain_type),intent(in) :: sub
! Array of numberf of DOF at interface nodes
      integer,intent(in) :: lnndfi
      integer,intent(out) :: nndfi(lnndfi)
! local vars
      integer :: nnodi, i, indnod

      if (.not.sub%is_interface_loaded) then
         write(*,*) 'DD_GET_SUBDOMAIN_INTERFACE_NNDF: Interface is not loaded yet for subdomain ',sub%isub
         call error_exit
      end if
      if (sub%nnodi .ne. lnndfi) then
         write(*,*) 'DD_GET_SUBDOMAIN_INTERFACE_NNDF: Size of array for interface not consistent.'
         write(*,*) 'nnodi :', sub%nnodi, 'array size: ',lnndfi 
         call error_exit
      end if

      ! if all checks are OK, construct array of numbers of DOF at interface nodes
      nnodi = sub%nnodi
      do i = 1,nnodi
         indnod   = sub%iin(i)
         nndfi(i) = sub%nndf(indnod)
      end do

end subroutine

!***********************************************
subroutine dd_get_subdomain_nndf(sub,nndf,lnndf)
!***********************************************
! Subroutine for getting subdomain nndf array
      use module_utils
      implicit none
! Subdomain structure
      type(subdomain_type),intent(in) :: sub
! Array of numberf of DOF at nodes
      integer,intent(in) :: lnndf
      integer,intent(out) :: nndf(lnndf)
! local vars
      integer :: nnod, i

      if (.not.sub%is_mesh_loaded) then
         write(*,*) 'DD_GET_SUBDOMAIN_NNDF: Mesh is not loaded yet for subdomain ',sub%isub
         call error_exit
      end if
      if (sub%nnod .ne. lnndf) then
         write(*,*) 'DD_GET_SUBDOMAIN_NNDF: Size of array not consistent.'
         write(*,*) 'nnod :', sub%nnod, 'array size: ',lnndf 
         call error_exit
      end if

      ! if all checks are OK, construct array of numbers of DOF 
      nnod = sub%nnod
      do i = 1,nnod
         nndf(i) = sub%nndf(i)
      end do

end subroutine

!**************************************************
subroutine dd_get_subdomain_isngn(sub,isngn,lisngn)
!**************************************************
! Subroutine for getting subdomain Indices of Subdomain Nodes in Global Numbering (ISNGN) array
      use module_utils
      implicit none
! Subdomain structure
      type(subdomain_type),intent(in) :: sub
! Array ISNGN
      integer,intent(in) :: lisngn
      integer,intent(out) :: isngn(lisngn)
! local vars
      integer :: nnod, i

      if (.not.sub%is_mesh_loaded) then
         write(*,*) 'DD_GET_SUBDOMAIN_ISNGN: Mesh is not loaded yet for subdomain ',sub%isub
         call error_exit
      end if
      if (sub%nnod .ne. lisngn) then
         write(*,*) 'DD_GET_SUBDOMAIN_ISNGN: Size of array not consistent.'
         write(*,*) 'nnod :', sub%nnod, 'array size: ',lisngn 
         call error_exit
      end if

      ! if all checks are OK, construct array of indices of subdomain nodes in global numbering
      nnod = sub%nnod
      do i = 1,nnod
         isngn(i) = sub%isngn(i)
      end do

end subroutine

!************************************************************
subroutine dd_get_interface_global_numbers(sub, iingn,liingn)
!************************************************************
! Subroutine for getting mapping of interface nodes into global node numbering
      use module_utils
      implicit none
! Subdomain structure
      type(subdomain_type),intent(in) :: sub
      integer,intent(in)  :: liingn
      integer,intent(out) ::  iingn(liingn)

      ! local vars
      integer :: nnodi
      integer :: i, ind

      ! check if mesh is loaded
      if (.not. sub%is_interface_loaded) then
         write(*,*) 'DD_GET_INTERFACE_GLOBAL_NUMBERS: Interface not loaded for subdomain: ',sub%isub
         call error_exit
      end if
      ! check dimensions
      if (sub%nnodi.ne.liingn) then
         write(*,*) 'DD_GET_INTERFACE_GLOBAL_NUMBERS: Interface dimensions mismatch for subdomain ',sub%isub
         call error_exit
      end if

      ! load data
      nnodi = sub%nnodi
      do i = 1,nnodi
         ! subdomain node number
         ind = sub%iin(i)
         ! global node number
         iingn(i) = sub%isngn(ind)
      end do
end subroutine

!***********************************************
subroutine dd_get_interface_nodes(sub, iin,liin)
!***********************************************
! Subroutine for getting mapping of interface nodes into subdomain node numbering
      use module_utils
      implicit none
! Subdomain structure
      type(subdomain_type),intent(in) :: sub
      integer,intent(in)  :: liin
      integer,intent(out) ::  iin(liin)

      ! local vars
      integer :: nnodi
      integer :: i

      ! check if mesh is loaded
      if (.not. sub%is_interface_loaded) then
         write(*,*) 'DD_GET_INTERFACE_NODES: Mesh not loaded for subdomain: ',sub%isub
         call error_exit
      end if
      ! check dimensions
      if (sub%nnodi.ne.liin) then
         write(*,*) 'DD_GET_INTERFACE_NODES: Interface dimensions mismatch for subdomain ',sub%isub
         call error_exit
      end if

      ! load data
      nnodi = sub%nnodi
      do i = 1,nnodi
         ! subdomain node number
         iin(i) = sub%iin(i)
      end do
end subroutine

!*********************************************************
subroutine dd_get_interface_variables(sub, iivsvn,liivsvn)
!*********************************************************
! Subroutine for getting mapping of interface variables into subdomain variables
      use module_utils
      implicit none
! Subdomain structure
      type(subdomain_type),intent(in) :: sub
      integer,intent(in)  :: liivsvn
      integer,intent(out) ::  iivsvn(liivsvn)

      ! local vars
      integer :: ndofi
      integer :: i

      ! check if interface is loaded
      if (.not. sub%is_interface_loaded) then
         write(*,*) 'DD_GET_INTERFACE_VARIABLES: Interface not loaded for subdomain: ',sub%isub
         call error_exit
      end if
      ! check dimensions
      if (sub%ndofi.ne.liivsvn) then
         write(*,*) 'DD_GET_INTERFACE_VARIABLES: Interface dimensions mismatch for subdomain ',sub%isub
         call error_exit
      end if

      ! load data
      ndofi = sub%ndofi
      do i = 1,ndofi
         iivsvn(i) = sub%iivsvn(i)
      end do
end subroutine

!*****************************************
subroutine dd_get_schur(sub, schur,lschur)
!*****************************************
! Subroutine for explicit construction of Schur complement
      use module_utils
      implicit none
! Subdomain structure
      type(subdomain_type),intent(inout) :: sub
      integer,intent(in)  ::  lschur
      real(kr),intent(out) ::  schur(lschur)

      ! local vars
      integer :: j, lmat, ndofi, pointschur

      integer              :: le
      real(kr),allocatable ::  e(:)

      ! check if interface is loaded
      if (.not. sub%is_interface_loaded) then
         write(*,*) 'DD_GET_SCHUR: Interface not loaded for subdomain: ',sub%isub
         call error_exit
      end if
      ! check dimensions
      ndofi = sub%ndofi
      lmat = ndofi*ndofi
      if (lmat.ne.lschur) then
         write(*,*) 'DD_GET_SCHUR: Interface dimensions mismatch for subdomain ',sub%isub
         call error_exit
      end if

      pointschur = 1
      le = ndofi
      allocate(e(le))
      do j = 1,ndofi
         ! construct vector of cartesian basis
         call zero(e,le)
         e(j) = 1._kr

         ! get column of S*I
         call dd_multiply_by_schur(sub,e,le,schur(pointschur),le)

         pointschur = pointschur + ndofi
      end do
      deallocate(e)

end subroutine

!*********************************************************************************************
subroutine dd_create_neighbouring(suba,lsuba, sub2proc,lsub2proc,indexsub,lindexsub, comm_all)
!*********************************************************************************************
! Subroutine for construction of array of interface and shared nodes with other subdomains
! based on subdomain data
      use module_utils
      use module_pp, only : pp_get_proc_for_sub 
      implicit none
      include "mpif.h"

! array of sub structure for actual subdomains
      integer,intent(in) ::                lsuba
      type(subdomain_type),intent(inout) :: suba(lsuba)

      integer,intent(in) :: lsub2proc
      integer,intent(in) ::  sub2proc(lsub2proc)
      integer,intent(in) :: lindexsub
      integer,intent(in) ::  indexsub(lindexsub)
      integer,intent(in) :: comm_all ! MPI communicator

      ! local vars
      integer :: isub_loc, isub, isubadj_loc, isubadj, nsub
      integer :: i, ia, inodis, inods, indg
      integer :: procadj
      integer :: kishnadj
      integer :: ierr, ireq
      integer :: nadj
      integer :: nnod
      integer :: nshared
      integer :: nsharedv, lcommvec, ndofn
      integer :: tag
      integer :: idofi, idofo, idofn, inod, inodi, &
                 kisngnadj, ndofi, ndofo, nnodi, nnodadj

      integer ::             lnshnadj
      integer,allocatable ::  nshnadj(:)
      integer ::             lkshvadj
      integer,allocatable ::  kshvadj(:)
      integer ::             lishnadj
      integer,allocatable ::  ishnadj(:)
      integer ::             lishared
      integer,allocatable ::  ishared(:)
      integer ::             lkinodes
      integer,allocatable ::  kinodes(:)
      integer ::             liin
      integer,allocatable ::  iin(:)
      integer ::             liivsvn
      integer,allocatable ::  iivsvn(:)
      integer ::             liovsvn
      integer,allocatable ::  iovsvn(:)
      integer ::             lkdof
      integer,allocatable ::  kdof(:)

      type sub_aux_type

         integer ::             lisngn
         integer,allocatable ::  isngn(:)
         integer ::             lnnodadj
         integer,allocatable ::  nnodadj(:)
         integer ::             lisngnadj
         integer,allocatable ::  isngnadj(:)

         ! MPI related arrays and variables
         integer ::            nreq
         integer ::            lrequest
         integer,allocatable :: request(:)
         integer             :: lstatarray1
         integer             :: lstatarray2
         integer,allocatable :: statarray(:,:)
      end type

      integer :: lsub_aux
      type(sub_aux_type),allocatable :: sub_aux(:)

      integer :: myid, nproc

! orient in communicators
      call MPI_COMM_RANK(comm_all,myid,ierr)
      call MPI_COMM_SIZE(comm_all,nproc,ierr)


      lsub_aux = lsuba
      allocate(sub_aux(lsub_aux))

      ! prepare subdomain data      
      do isub_loc = 1,lindexsub
         isub = indexsub(isub_loc)
         nsub = suba(isub_loc)%nsub

         if (.not. suba(isub_loc)%is_mesh_loaded) then
            write(*,*) 'DD_CREATE_NEIGHBOURING: myid =',myid,', Mesh not loaded for subdomain: ',isub
            call error_exit
         end if

         ! load data
         nadj  = suba(isub_loc)%nadj
         nnod  = suba(isub_loc)%nnod

         sub_aux(isub_loc)%lnnodadj = nadj
         allocate(sub_aux(isub_loc)%nnodadj(sub_aux(isub_loc)%lnnodadj))

         ! MPI arrays
         sub_aux(isub_loc)%lrequest = 2*nadj
         allocate(sub_aux(isub_loc)%request(sub_aux(isub_loc)%lrequest))
         sub_aux(isub_loc)%lstatarray1 = MPI_STATUS_SIZE
         sub_aux(isub_loc)%lstatarray2 = sub_aux(isub_loc)%lrequest
         allocate(sub_aux(isub_loc)%statarray(sub_aux(isub_loc)%lstatarray1,sub_aux(isub_loc)%lstatarray2))

! Determine sizes of interace of my neighbours
         ireq = 0
         do ia = 1,nadj
            ! get index of neighbour
            isubadj = suba(isub_loc)%iadj(ia)
   
            ! who owns this subdomain?
            call pp_get_proc_for_sub(isubadj,comm_all,sub2proc,lsub2proc,procadj)
   
            if (procadj .eq. myid) then
               ! I have subdomain data, simply copy necessary arrays
               call get_index(isubadj,indexsub,lindexsub,isubadj_loc)
               if (.not. suba(isubadj_loc)%is_mesh_loaded) then
                  write(*,*) 'DD_CREATE_NEIGHBOURING: myid =',myid,', Mesh not loaded for subdomain: ',isubadj
                  call error_exit
               end if
   
               sub_aux(isub_loc)%nnodadj(ia) = suba(isubadj_loc)%nnod
   
            else 
               ! interchange via MPI
   
               ! receive his data
               tag = isubadj*nsub + isub
               ireq = ireq + 1
               call MPI_IRECV(sub_aux(isub_loc)%nnodadj(ia),1,MPI_INTEGER,procadj,tag,comm_all,&
                              sub_aux(isub_loc)%request(ireq),ierr)
            !   print *, 'myid =',myid,'receiving', sub_aux(isub_loc)%nnodiadj(ia),'from',procadj,' tag',tag

               ! send him my data
               tag = isub*nsub + isubadj
               ireq = ireq + 1
               call MPI_ISEND(nnod,1,MPI_INTEGER,procadj,tag,comm_all,sub_aux(isub_loc)%request(ireq),ierr)
            !   print *, 'myid =',myid,'Sending', nnodi,'to ',procadj,' tag',tag
            end if
         end do
         sub_aux(isub_loc)%nreq = ireq
      end do

      ! waiting for communication to complete
      do isub_loc = 1,lindexsub
         isub = indexsub(isub_loc)
         call MPI_WAITALL(sub_aux(isub_loc)%nreq, sub_aux(isub_loc)%request, sub_aux(isub_loc)%statarray, ierr)
      end do

      do isub_loc = 1,lindexsub
         isub = indexsub(isub_loc)
         nsub = suba(isub_loc)%nsub

         ! load data
         nadj  = suba(isub_loc)%nadj
         nnod  = suba(isub_loc)%nnod

! Allocate array for global node numbers for neighbours
         sub_aux(isub_loc)%lisngnadj = sum(sub_aux(isub_loc)%nnodadj)
         allocate(sub_aux(isub_loc)%isngnadj(sub_aux(isub_loc)%lisngnadj))
         ! allocate data for myself
         sub_aux(isub_loc)%lisngn = nnod
         allocate(sub_aux(isub_loc)%isngn(sub_aux(isub_loc)%lisngn))
         call dd_get_subdomain_isngn(suba(isub_loc),sub_aux(isub_loc)%isngn,sub_aux(isub_loc)%lisngn)
   
! Interchange nodes in global numbering
         ireq = 0
         kisngnadj = 0
         do ia = 1,nadj
            ! get index of neighbour
            isubadj = suba(isub_loc)%iadj(ia)
            
            ! who owns this subdomain?
            call pp_get_proc_for_sub(isubadj,comm_all,sub2proc,lsub2proc,procadj)

            nnodadj = sub_aux(isub_loc)%nnodadj(ia)
   
            if (procadj .eq. myid) then
               ! I have subdomain data, simply copy necessary arrays
               call get_index(isubadj,indexsub,lindexsub,isubadj_loc)
               do i = 1,nnodadj
                  sub_aux(isub_loc)%isngnadj(kisngnadj + i) = suba(isubadj_loc)%isngn(i)
               end do
   
            else 
               ! interchange via MPI

               ! receive his data
               tag = isubadj*nsub + isub
               ireq = ireq + 1
               call MPI_IRECV(sub_aux(isub_loc)%isngnadj(kisngnadj + 1),nnodadj,MPI_INTEGER,&
                              procadj,tag,comm_all,sub_aux(isub_loc)%request(ireq),ierr)
   
               ! send him my data
               tag = isub*nsub + isubadj
               ireq = ireq + 1
               call MPI_ISEND(sub_aux(isub_loc)%isngn,sub_aux(isub_loc)%lisngn,MPI_INTEGER,&
                              procadj,tag,comm_all,sub_aux(isub_loc)%request(ireq),ierr)
            end if
            kisngnadj = kisngnadj + nnodadj
         end do
         sub_aux(isub_loc)%nreq = ireq
      end do

      ! waiting for communication to complete
      do isub_loc = 1,lindexsub
         isub = indexsub(isub_loc)
         call MPI_WAITALL(sub_aux(isub_loc)%nreq, sub_aux(isub_loc)%request, sub_aux(isub_loc)%statarray, ierr)
      end do

      ! Compare data with neighbours to detect interface
      do isub_loc = 1,lindexsub
         isub = indexsub(isub_loc)

         ! load data
         nadj  = suba(isub_loc)%nadj
         nnod  = suba(isub_loc)%nnod

         lkinodes = nnod
         allocate(kinodes(lkinodes))
         kinodes = 0
         lishared = nnod
         allocate(ishared(lishared))

         kisngnadj = 0
         do ia = 1,nadj
            ! get index of neighbour
            isubadj = suba(isub_loc)%iadj(ia)
   
            nnodadj = sub_aux(isub_loc)%nnodadj(ia)
   
            call get_array_intersection(sub_aux(isub_loc)%isngn,sub_aux(isub_loc)%lisngn,&
                                        sub_aux(isub_loc)%isngnadj(kisngnadj + 1),nnodadj,&
                                        ishared,lishared,nshared)
            if (nshared.le.0) then
               write(*,*) 'DD_CREATE_NEIGHBOURING: myid =',myid,', It seems as subdomain ',&
                           isub, ' does not share nodes with neighbour ',isubadj
               call error_exit
            end if
   
            ! mark interface
            nsharedv = 0
            do i = 1,nshared
               indg = ishared(i)
               call get_index(indg,suba(isub_loc)%isngn,nnod,inods)

               kinodes(inods) = kinodes(inods) + 1
            end do
   
            kisngnadj  = kisngnadj  + nnodadj
         end do

! find number of interface nodes
         nnodi = 0 
         ndofi = 0
         ndofo = 0
         do inod = 1,nnod
            if (kinodes(inod).gt.0) then
               nnodi = nnodi + 1
               ndofi = ndofi + suba(isub_loc)%nndf(inod)
            else
               ndofo = ndofo + suba(isub_loc)%nndf(inod)
            end if
         end do
! generate mapping of interface nodes to subdomain nodes and the same for dofs 
         liin = nnodi
         allocate(iin(liin))
         liivsvn = ndofi
         allocate(iivsvn(liivsvn))
         liovsvn = ndofo
         allocate(iovsvn(liovsvn))
         
         ! create array kdof
         lkdof = nnod
         allocate(kdof(lkdof))
         if (lkdof.gt.0) then
            kdof(1) = 0
            do inod = 2,nnod
               kdof(inod) = kdof(inod-1) + suba(isub_loc)%nndf(inod-1)
            end do
         end if


         inodi = 0
         idofi = 0
         idofo = 0
         do inod = 1,nnod
            ndofn  = suba(isub_loc)%nndf(inod)

            if (kinodes(inod).gt.0) then
               inodi = inodi + 1

               iin(inodi) = inod
               do idofn = 1,ndofn 
                  idofi = idofi + 1

                  iivsvn(idofi) = kdof(inod) + idofn
               end do
            else
               do idofn = 1,ndofn 
                  idofo = idofo + 1

                  iovsvn(idofo) = kdof(inod) + idofn
               end do
            end if
         end do
         deallocate(kdof)
         deallocate(kinodes)

         ! load data to structure
         call dd_upload_sub_interface(suba(isub_loc), nnodi, ndofi, ndofo, &
                                      iin,liin, iivsvn,liivsvn, iovsvn,liovsvn)
         ! debug
         !print *,'isub = ',isub,'iin',iin
         !print *,'isub = ',isub,'iivsvn',iivsvn
         !print *,'isub = ',isub,'iovsvn',iovsvn
         deallocate(iin)
         deallocate(iivsvn)
         deallocate(iovsvn)
   
         ! Compare data with neighbours to create neighbouring

         lnshnadj = nadj
         allocate(nshnadj(lnshnadj))
         lkshvadj = nadj + 1
         allocate(kshvadj(lkshvadj))
         lishnadj = nadj * nnodi
         allocate(ishnadj(lishnadj))

         kisngnadj = 0
         kishnadj  = 0
         lcommvec  = 0
         kshvadj(1) = 1
         do ia = 1,nadj
            ! get index of neighbour
            isubadj = suba(isub_loc)%iadj(ia)
   
            nnodadj = sub_aux(isub_loc)%nnodadj(ia)
   
            call get_array_intersection(sub_aux(isub_loc)%isngn,sub_aux(isub_loc)%lisngn,&
                                        sub_aux(isub_loc)%isngnadj(kisngnadj + 1),nnodadj,&
                                        ishared,lishared,nshared)
            if (nshared.le.0) then
               write(*,*) 'DD_CREATE_NEIGHBOURING: myid =',myid,', It seems as subdomain ',&
                           isub, ' does not share nodes with neighbour ',isubadj
               call error_exit
            end if
   
            nshnadj(ia) = nshared
   
            ! load shared nodes 
            nsharedv = 0
            do i = 1,nshared
               indg = ishared(i)
               call get_index(indg,suba(isub_loc)%isngn,suba(isub_loc)%lisngn,inods)
               call get_index(inods,suba(isub_loc)%iin,suba(isub_loc)%nnodi,inodis)

               ndofn = suba(isub_loc)%nndf(inods)
               nsharedv = nsharedv + ndofn

               ishnadj(kishnadj + i) = inodis
            end do
   
            kshvadj(ia + 1) = kshvadj(ia) + nsharedv
            lcommvec = lcommvec + nsharedv

            kisngnadj  = kisngnadj  + nnodadj
            kishnadj = kishnadj + nshared
         end do
         deallocate(ishared)
   
         ! load info about shared nodes into structure
         suba(isub_loc)%lnshnadj = lnshnadj
         allocate(suba(isub_loc)%nshnadj(lnshnadj))
         do i = 1,nadj
            suba(isub_loc)%nshnadj(i) = nshnadj(i)
         end do
   
         ! load info about number of shared dof into structure
         suba(isub_loc)%lkshvadj = lkshvadj
         allocate(suba(isub_loc)%kshvadj(lkshvadj))
         do i = 1,lkshvadj
            suba(isub_loc)%kshvadj(i) = kshvadj(i)
         end do

         ! truncate array ishnadj to really used indices
         suba(isub_loc)%lishnadj = sum(nshnadj)
         allocate(suba(isub_loc)%ishnadj(suba(isub_loc)%lishnadj))
         do i = 1,suba(isub_loc)%lishnadj
            suba(isub_loc)%ishnadj(i) = ishnadj(i)
         end do

         ! prepare communication arrays
         suba(isub_loc)%lcommvec = lcommvec
         
         ! activate the flag
         suba(isub_loc)%is_neighbouring_ready = .true.
   
         ! debug
         !print *,'isub = ',isub,'nshnadj',nshnadj
         !print *,'isub = ',isub,'kshvadj',kshvadj
         !print *,'isub = ',isub,'ishnadj',ishnadj
         deallocate(nshnadj)
         deallocate(kshvadj)
         deallocate(ishnadj)

         deallocate(sub_aux(isub_loc)%isngn)
         deallocate(sub_aux(isub_loc)%isngnadj)
         deallocate(sub_aux(isub_loc)%nnodadj)
         deallocate(sub_aux(isub_loc)%request)
         deallocate(sub_aux(isub_loc)%statarray)
      end do

      deallocate(sub_aux)

end subroutine

!!*******************************************************
!subroutine dd_create_neighbouringi(myid, nsub, comm_all)
!!*******************************************************
!! Subroutine for construction of array of shared nodes with other subdomains
!! based on interface data
!      use module_utils
!      implicit none
!      include "mpif.h"
!
!      integer,intent(in) :: myid, nsub
!      integer,intent(in) :: comm_all ! MPI communicator
!
!      ! local vars
!      integer :: i, ia, inds, inodis, inods, indg, isub
!      integer :: isubadj, procadj
!      integer :: kiinadj, kishnadj
!      integer :: ierr, ireq
!      integer :: nnadj
!      integer :: nadj
!      integer :: nnod
!      integer :: nnodi
!      integer :: nshared
!      integer :: nsharedv, lcommvec, ndofn
!      integer :: tag
!
!      integer ::             lnshnadj
!      integer,allocatable ::  nshnadj(:)
!      integer ::             lkshvadj
!      integer,allocatable ::  kshvadj(:)
!      integer ::             lishnadj
!      integer,allocatable ::  ishnadj(:)
!      integer ::             lishared
!      integer,allocatable ::  ishared(:)
!
!      type sub_aux_type
!
!         integer ::             lnnodiadj
!         integer,allocatable ::  nnodiadj(:)
!         integer ::             liinadj
!         integer,allocatable ::  iinadj(:)
!         integer ::             liingn
!         integer,allocatable ::  iingn(:)
!
!         ! MPI related arrays and variables
!         integer ::            nreq
!         integer ::            lrequest
!         integer,allocatable :: request(:)
!         integer             :: lstatarray1
!         integer             :: lstatarray2
!         integer,allocatable :: statarray(:,:)
!      end type
!
!      integer :: lsub_aux
!      type(sub_aux_type),allocatable :: sub_aux(:)
!
!      lsub_aux = nsub
!      allocate(sub_aux(lsub_aux))
!
!      ! check if I store the subdomain
!      do isub = 1,nsub
!
!         if (.not. sub(isub)%proc .eq. myid) then
!            if (debug) then
!               write(*,*) 'DD_CREATE_NEIGHBOURING: myid =',myid,', not my subdomain: ',isub
!            end if
!            cycle
!         end if
!         if (.not. sub(isub)%is_mesh_loaded) then
!            write(*,*) 'DD_CREATE_NEIGHBOURING: myid =',myid,', Mesh not loaded for subdomain: ',isub
!            call error_exit
!         end if
!
!         ! load data
!         nnod  = sub(isub)%nnod
!         nnodi = sub(isub)%nnodi
!         nadj  = sub(isub)%nadj
!
!         sub_aux(isub)%lnnodiadj = nadj
!         allocate(sub_aux(isub)%nnodiadj(sub_aux(isub)%lnnodiadj))
!
!         sub_aux(isub)%lrequest = 2*nadj
!         allocate(sub_aux(isub)%request(sub_aux(isub)%lrequest))
!         sub_aux(isub)%lstatarray1 = MPI_STATUS_SIZE
!         sub_aux(isub)%lstatarray2 = sub_aux(isub)%lrequest
!         allocate(sub_aux(isub)%statarray(sub_aux(isub)%lstatarray1,sub_aux(isub)%lstatarray2))
!
!! Determine sizes of interace of my neighbours
!         ireq = 0
!         do ia = 1,nadj
!            ! get index of neighbour
!            isubadj = sub(isub)%iadj(ia)
!   
!            ! who owns this subdomain?
!            procadj = sub(isubadj)%proc
!   
!            if (procadj .eq. myid) then
!               ! I have subdomain data, simply copy necessary arrays
!               if (.not. sub(isubadj)%is_mesh_loaded) then
!                  write(*,*) 'DD_CREATE_NEIGHBOURING: myid =',myid,', Mesh not loaded for subdomain: ',isubadj
!                  call error_exit
!               end if
!   
!               sub_aux(isub)%nnodiadj(ia) = sub(isubadj)%nnodi
!   
!            else 
!               ! interchange via MPI
!   
!               ! send him my data
!               tag = isub*nsub + isubadj
!               ireq = ireq + 1
!               call MPI_ISEND(nnodi,1,MPI_INTEGER,procadj,tag,comm_all,sub_aux(isub)%request(ireq),ierr)
!            !   print *, 'myid =',myid,'Sending', nnodi,'to ',procadj,' tag',tag
!               ! receive his data
!               tag = isubadj*nsub + isub
!               ireq = ireq + 1
!               call MPI_IRECV(sub_aux(isub)%nnodiadj(ia),1,MPI_INTEGER,procadj,tag,comm_all,sub_aux(isub)%request(ireq),ierr)
!            !   print *, 'myid =',myid,'receiving', sub_aux(isub)%nnodiadj(ia),'from',procadj,' tag',tag
!            end if
!         end do
!         sub_aux(isub)%nreq = ireq
!      end do
!
!      do isub = 1,nsub
!         if (.not. sub(isub)%proc .eq. myid) then
!            cycle
!         end if
!         call MPI_WAITALL(sub_aux(isub)%nreq, sub_aux(isub)%request, sub_aux(isub)%statarray, ierr)
!      end do
!
!      do isub = 1,nsub
!         if (.not. sub(isub)%proc .eq. myid) then
!            cycle
!         end if
!
!         ! load data
!         nnod  = sub(isub)%nnod
!         nnodi = sub(isub)%nnodi
!         nadj  = sub(isub)%nadj
!
!! Allocate array for global node numbers at interface of for neighbours
!         sub_aux(isub)%liinadj = sum(sub_aux(isub)%nnodiadj)
!         allocate(sub_aux(isub)%iinadj(sub_aux(isub)%liinadj))
!   
!
!! prepare array of interface nodes in global numbering
!         sub_aux(isub)%liingn = nnodi
!         allocate(sub_aux(isub)%iingn(sub_aux(isub)%liingn))
!         call dd_get_interface_global_numbers(myid, isub, sub_aux(isub)%iingn,sub_aux(isub)%liingn)
!
!!     Interchange interface nodes in global numbering
!         ireq = 0
!         kiinadj = 0
!         do ia = 1,nadj
!         ! get index of neighbour
!            isubadj = sub(isub)%iadj(ia)
!
!         ! who owns this subdomain?
!            procadj = sub(isubadj)%proc
!
!            nnadj = sub_aux(isub)%nnodiadj(ia)
!   
!            if (procadj .eq. myid) then
!               ! I have subdomain data, simply copy necessary arrays
!               do i = 1,nnadj
!                  inds     = sub(isubadj)%iin(i)
!                  sub_aux(isub)%iinadj(kiinadj + i) = sub(isubadj)%isngn(inds)
!               end do
!   
!            else 
!               ! interchange via MPI
!   
!               ! send him my data
!               tag = isub*nsub + isubadj
!               ireq = ireq + 1
!               call MPI_ISEND(sub_aux(isub)%iingn,sub_aux(isub)%liingn,MPI_INTEGER,&
!                              procadj,tag,comm_all,sub_aux(isub)%request(ireq),ierr)
!               ! receive his data
!               tag = isubadj*nsub + isub
!               ireq = ireq + 1
!               call MPI_IRECV(sub_aux(isub)%iinadj(kiinadj + 1),nnadj,MPI_INTEGER,&
!                              procadj,tag,comm_all,sub_aux(isub)%request(ireq),ierr)
!            end if
!            kiinadj = kiinadj + nnadj
!         end do
!         sub_aux(isub)%nreq = ireq
!      end do
!
!      do isub = 1,nsub
!         if (.not. sub(isub)%proc .eq. myid) then
!            cycle
!         end if
!         call MPI_WAITALL(sub_aux(isub)%nreq, sub_aux(isub)%request, sub_aux(isub)%statarray, ierr)
!      end do
!
!      ! Compare data at interface with neighbours 
!      do isub = 1,nsub
!         if (.not. sub(isub)%proc .eq. myid) then
!            cycle
!         end if
!
!         ! load data
!         nnod  = sub(isub)%nnod
!         nnodi = sub(isub)%nnodi
!         nadj  = sub(isub)%nadj
!
!         lishared = nnodi
!         allocate(ishared(lishared))
!         lnshnadj = nadj
!         allocate(nshnadj(lnshnadj))
!         lkshvadj = nadj + 1
!         allocate(kshvadj(lkshvadj))
!         lishnadj = nadj * nnodi
!         allocate(ishnadj(lishnadj))
!
!         kiinadj = 0
!         kishnadj = 0
!         lcommvec = 0
!         kshvadj(1) = 1
!         do ia = 1,nadj
!            ! get index of neighbour
!            isubadj = sub(isub)%iadj(ia)
!   
!            nnadj = sub_aux(isub)%nnodiadj(ia)
!   
!            call get_array_intersection(sub_aux(isub)%iingn,sub_aux(isub)%liingn,sub_aux(isub)%iinadj(kiinadj + 1),nnadj,&
!                                        ishared,lishared,nshared)
!            if (nshared.le.0) then
!               write(*,*) 'DD_CREATE_NEIGHBOURING: myid =',myid,', It seems as subdomain ',&
!                           isub, ' does not share nodes with neighbour ',isubadj
!               call error_exit
!            end if
!   
!            nshnadj(ia) = nshared
!   
!            ! load shared nodes 
!            nsharedv = 0
!            do i = 1,nshared
!               indg = ishared(i)
!               call get_index(indg,sub(isub)%isngn,nnod,inods)
!               call get_index(inods,sub(isub)%iin,nnodi,inodis)
!
!               ndofn = sub(isub)%nndf(inods)
!               nsharedv = nsharedv + ndofn
!
!               ishnadj(kishnadj + i) = inodis
!            end do
!   
!            kshvadj(ia + 1) = kshvadj(ia) + nsharedv
!            lcommvec = lcommvec + nsharedv
!
!            kiinadj  = kiinadj  + nnadj
!            kishnadj = kishnadj + nshared
!         end do
!   
!         ! load info about shared nodes into structure
!         sub(isub)%lnshnadj = lnshnadj
!         allocate(sub(isub)%nshnadj(lnshnadj))
!         do i = 1,nadj
!            sub(isub)%nshnadj(i) = nshnadj(i)
!         end do
!   
!         ! load info about number of shared dof into structure
!         sub(isub)%lkshvadj = lkshvadj
!         allocate(sub(isub)%kshvadj(lkshvadj))
!         do i = 1,lkshvadj
!            sub(isub)%kshvadj(i) = kshvadj(i)
!         end do
!
!         ! truncate array ishnadj to really used indices
!         sub(isub)%lishnadj = sum(nshnadj)
!         allocate(sub(isub)%ishnadj(sub(isub)%lishnadj))
!         do i = 1,sub(isub)%lishnadj
!            sub(isub)%ishnadj(i) = ishnadj(i)
!         end do
!
!         ! prepare communication arrays
!         sub(isub)%lcommvec = lcommvec
!         
!         ! activate the flag
!         sub(isub)%is_neighbouring_ready = .true.
!   
!         deallocate(ishared)
!         deallocate(nshnadj)
!         deallocate(kshvadj)
!         deallocate(ishnadj)
!
!         deallocate(sub_aux(isub)%iingn)
!         deallocate(sub_aux(isub)%iinadj)
!         deallocate(sub_aux(isub)%nnodiadj)
!         deallocate(sub_aux(isub)%request)
!         deallocate(sub_aux(isub)%statarray)
!      end do
!
!      deallocate(sub_aux)
!
!end subroutine

!***************************************
subroutine dd_comm_upload(sub, vec,lvec)
!***************************************
! Subroutine that loads data from subdomain vector VEC for communication
      use module_utils
      implicit none
! Subdomain structure
      type(subdomain_type),intent(inout) :: sub

      integer,intent(in) :: lvec
      real(kr),intent(in) :: vec(lvec)

      ! local vars
      integer ::            lkdofi
      integer,allocatable :: kdofi(:)

      integer :: ndofi, nnodi, lcommvec
      integer :: i, ia, inadj, indcommvec, indn, indshni, kishnadj, nadj, ndofn, nnadj

      ! check neighbouring
      if (.not. sub%is_neighbouring_ready) then
         write(*,*) 'DD_COMM_UPLOAD: Neighbouring is not ready for subdomain ',sub%isub
         call error_exit
      end if

      ndofi = sub%ndofi
      nnodi = sub%nnodi
      ! check size
      if (lvec .ne. ndofi) then
         write(*,*) 'DD_COMM_UPLOAD: Data size mismatch: lvec ',lvec, 'ndofi',ndofi
         call error_exit
      end if

      lcommvec = sub%lcommvec
      if (.not.allocated(sub%commvec_out)) then
         allocate(sub%commvec_out(lcommvec))
      end if
      if (.not.allocated(sub%commvec_in)) then
         allocate(sub%commvec_in(lcommvec))
      end if

      ! prepare array kdofi
      lkdofi = nnodi + 1
      allocate(kdofi(lkdofi))
      if (lkdofi.gt.0) then
         kdofi(1) = 1
         do i = 1,nnodi
            indn = sub%iin(i)
            ndofn = sub%nndf(indn)
            
            kdofi(i + 1) = kdofi(i) + ndofn
         end do
      end if

      ! load vector at interface into communication vector COMMVEC_OUT
      nadj  = sub%nadj
      kishnadj = 0
      indcommvec = 0
      do ia = 1,nadj
         nnadj = sub%nshnadj(ia)

         do inadj = 1,nnadj
            indshni = sub%ishnadj(kishnadj + inadj)

            ndofn = kdofi(indshni + 1) - kdofi(indshni)

            do i = 1,ndofn
               indcommvec = indcommvec + 1
               sub%commvec_out(indcommvec) = vec(kdofi(indshni)-1 + i)
            end do
         end do

         kishnadj = kishnadj + nnadj
      end do
         
      deallocate(kdofi)

end subroutine

!***************************************************************************************
subroutine dd_comm_swapdata(suba,lsuba, indexsub,lindexsub, sub2proc,lsub2proc,comm_all)
!***************************************************************************************
! Subroutine for interchange of data at subdomain interfaces using MPI
      use module_utils
      use module_pp, only : pp_get_proc_for_sub 
      implicit none
      include "mpif.h"

! array of sub structure for actual subdomains
      integer,intent(in) ::                lsuba
      type(subdomain_type),intent(inout) :: suba(lsuba)

      integer,intent(in) :: lindexsub
      integer,intent(in) ::  indexsub(lindexsub)
      integer,intent(in) :: lsub2proc
      integer,intent(in) ::  sub2proc(lsub2proc)
      integer,intent(in) :: comm_all ! MPI communicator

      ! local vars
      integer :: i, ia, isub, isub_loc, nsub
      integer :: isubadj, isubadj_loc, procadj

      integer :: myid, ierr, tag, ireq
      integer :: isubadj_ia
      integer :: kneib_isub
      integer :: kneib_isubadj
      integer :: lcommvec, nadj, nsharedv

      type sub_aux_type

         ! MPI related arrays and variables
         integer ::            nreq
         integer ::            lrequest
         integer,allocatable :: request(:)
         integer             :: lstatarray1
         integer             :: lstatarray2
         integer,allocatable :: statarray(:,:)
      end type

      integer :: lsub_aux
      type(sub_aux_type),allocatable :: sub_aux(:)

! orient in communicator
      call MPI_COMM_RANK(comm_all,myid,ierr)

      lsub_aux = lsuba
      allocate(sub_aux(lsub_aux))

      do isub_loc = 1,lindexsub
         isub = indexsub(isub_loc)
         nsub = suba(isub_loc)%nsub

         if (.not. suba(isub_loc)%is_neighbouring_ready) then
            write(*,*) 'DD_COMM_SWAPDATA: Neighbouring is not ready for subdomain ',isub
            call error_exit
         end if
         if (.not. allocated(suba(isub_loc)%commvec_out)) then
            write(*,*) 'DD_COMM_SWAPDATA: Array COMMVEC_OUT not allocated.  Perhaps missing call to dd_comm_upload. ',isub
            call error_exit
         end if

         ! load data
         nadj  = suba(isub_loc)%nadj
         lcommvec = suba(isub_loc)%lcommvec

         sub_aux(isub_loc)%lrequest = 2*nadj
         allocate(sub_aux(isub_loc)%request(sub_aux(isub_loc)%lrequest))
         sub_aux(isub_loc)%lstatarray1 = MPI_STATUS_SIZE
         sub_aux(isub_loc)%lstatarray2 = sub_aux(isub_loc)%lrequest
         allocate(sub_aux(isub_loc)%statarray(sub_aux(isub_loc)%lstatarray1,sub_aux(isub_loc)%lstatarray2))

!     Interchange interface nodes in global numbering
         ireq = 0
         do ia = 1,nadj
         ! get index of neighbour
            isubadj = suba(isub_loc)%iadj(ia)

            nsharedv = suba(isub_loc)%kshvadj(ia + 1) - suba(isub_loc)%kshvadj(ia)

            ! who owns this subdomain?
            call pp_get_proc_for_sub(isubadj,comm_all,sub2proc,lsub2proc,procadj)

            if (procadj .eq. myid) then
               ! I have subdomain data, simply copy necessary arrays
               call get_index(isubadj,indexsub,lindexsub,isubadj_loc)
               call get_index(isub,suba(isubadj_loc)%iadj,suba(isubadj_loc)%nadj,isubadj_ia)
               kneib_isub    = suba(isub_loc)%kshvadj(ia) - 1
               if ((suba(isubadj_loc)%kshvadj(isubadj_ia + 1) - suba(isubadj_loc)%kshvadj(isubadj_ia)) &
                  .ne. nsharedv) then
                  write(*,*) 'DD_COMM_SWAPDATA: Inconsistent lenght of copied array.'
                  call error_exit
               end if

               kneib_isubadj = suba(isubadj_loc)%kshvadj(isubadj_ia) - 1
               do i = 1,nsharedv
                  suba(isub_loc)%commvec_in(kneib_isub + i) = suba(isubadj_loc)%commvec_out(kneib_isubadj + i)
               end do
   
            else 
               ! interchange via MPI
   
               ! receive his data
               tag = isubadj*nsub + isub
               ireq = ireq + 1
               call MPI_IRECV(suba(isub_loc)%commvec_in(suba(isub_loc)%kshvadj(ia)),nsharedv,MPI_DOUBLE_PRECISION,&
                              procadj,tag,comm_all,sub_aux(isub_loc)%request(ireq),ierr)
               ! send him my data
               tag = isub*nsub + isubadj
               ireq = ireq + 1
               call MPI_ISEND(suba(isub_loc)%commvec_out(suba(isub_loc)%kshvadj(ia)),nsharedv,MPI_DOUBLE_PRECISION,&
                              procadj,tag,comm_all,sub_aux(isub_loc)%request(ireq),ierr)
            end if
         end do
         sub_aux(isub_loc)%nreq = ireq
      end do

      do isub_loc = 1,lindexsub
         call MPI_WAITALL(sub_aux(isub_loc)%nreq, sub_aux(isub_loc)%request, sub_aux(isub_loc)%statarray, ierr)
      end do

      ! clear memory
      do isub_loc = 1,lindexsub
         deallocate(sub_aux(isub_loc)%request)
         deallocate(sub_aux(isub_loc)%statarray)
      end do
      deallocate(sub_aux)

end subroutine

!*****************************************
subroutine dd_comm_download(sub, vec,lvec)
!*****************************************
! Subroutine that downloads data from communication to subdomain vector VEC
      use module_utils
      implicit none
! Subdomain structure
      type(subdomain_type),intent(inout) :: sub

      integer,intent(in) :: lvec
      real(kr),intent(out) :: vec(lvec)

      ! local vars
      integer ::            lkdofi
      integer,allocatable :: kdofi(:)

      integer :: ndofi, nnodi, lcommvec
      integer :: i, ia, inadj, indcommvec, indn, indshni, kishnadj, nadj, ndofn, nnadj

      ! check neighbouring
      if (.not. sub%is_neighbouring_ready) then
         write(*,*) 'DD_COMM_DOWNLOAD: Neighbouring is not ready for subdomain ',sub%isub
         call error_exit
      end if
      if (.not.allocated(sub%commvec_in)) then
         write(*,*) 'DD_COMM_DOWNLOAD: Communicated vector commvec_in not ready.'
         call error_exit
      end if

      lcommvec = sub%lcommvec
      ndofi = sub%ndofi
      nnodi = sub%nnodi

      ! check size
      if (lvec .ne. ndofi) then
         write(*,*) 'DD_COMM_DOWNLOAD: Data size mismatch: lvec ',lvec, 'ndofi',ndofi
         call error_exit
      end if

      ! prepare array kdofi
      lkdofi = nnodi + 1
      allocate(kdofi(lkdofi))
      kdofi(1) = 1
      do i = 1,nnodi
         indn = sub%iin(i)
         ndofn = sub%nndf(indn)
         
         kdofi(i + 1) = kdofi(i) + ndofn
      end do

      ! download vector at interface from communication vector COMMVEC_IN
      ! summ up repeated entries
      call zero(vec,lvec)
      nadj  = sub%nadj
      kishnadj = 0
      indcommvec = 0
      do ia = 1,nadj
         nnadj = sub%nshnadj(ia)

         do inadj = 1,nnadj
            indshni = sub%ishnadj(kishnadj + inadj)

            ndofn = kdofi(indshni + 1) - kdofi(indshni)

            do i = 1,ndofn
               indcommvec = indcommvec + 1
               vec(kdofi(indshni)-1 + i) = vec(kdofi(indshni)-1 + i) + sub%commvec_in(indcommvec)
            end do
         end do

         kishnadj = kishnadj + nnadj
      end do
         
      deallocate(kdofi)

      deallocate(sub%commvec_out)
      deallocate(sub%commvec_in)

end subroutine

!*****************************************************************************************
subroutine dd_weights_prepare(suba,lsuba, sub2proc,lsub2proc,indexsub,lindexsub, comm_all)
!*****************************************************************************************
! Subroutine for preparing weight at interface matrix
      use module_utils
      use module_pp, only : pp_get_proc_for_sub 
      implicit none
      include "mpif.h"

! array of sub structure for actual subdomains
      integer,intent(in) ::                lsuba
      type(subdomain_type),intent(inout) :: suba(lsuba)

      integer,intent(in) :: lsub2proc
      integer,intent(in) ::  sub2proc(lsub2proc)
      integer,intent(in) :: lindexsub
      integer,intent(in) ::  indexsub(lindexsub)
      integer,intent(in) :: comm_all ! MPI communicator

      ! local vars
      integer :: isub_loc, i
      integer :: ndofi

      integer ::             lrhoi
      real(kr),allocatable :: rhoi(:)
      real(kr),allocatable :: rhoiaux(:)
      integer ::             lwi
      real(kr),allocatable :: wi(:)

      ! what kind of weights should be used?
      ! choose only one
      logical :: weight_by_cardinality = .false.
      logical :: weight_by_stiffness   = .true.

      ! Prepare data for communication
      do isub_loc = 1,lindexsub
         ndofi = suba(isub_loc)%ndofi

         lrhoi = ndofi
         allocate(rhoi(lrhoi))

         if (weight_by_cardinality) then
            do i = 1,lrhoi
               rhoi(i) = 1._kr
            end do
         else if (weight_by_stiffness) then
            call dd_get_interface_diagonal(suba(isub_loc), rhoi,lrhoi)
         else
            write(*,*) 'DD_PREPARE_WEIGHTS: Type of weight not specified.'
            call error_exit
         end if

         call dd_comm_upload(suba(isub_loc), rhoi,lrhoi)

         deallocate(rhoi)
      end do

      ! Interchange data
      call dd_comm_swapdata(suba,lsuba, indexsub,lindexsub, sub2proc,lsub2proc,comm_all)

      ! Download communicated data 
      do isub_loc = 1,lindexsub
         ndofi = suba(isub_loc)%ndofi

         lrhoi = ndofi
         allocate(rhoi(lrhoi))
         allocate(rhoiaux(lrhoi))
         lwi = ndofi
         allocate(wi(lwi))

         call dd_comm_download(suba(isub_loc), rhoiaux,lrhoi)

         if (weight_by_cardinality) then
            do i = 1,lrhoi
               rhoi(i) = 1._kr
            end do
         else if (weight_by_stiffness) then
            call dd_get_interface_diagonal(suba(isub_loc), rhoi,lrhoi)
         end if

         ! compute weight
         do i = 1,ndofi
            wi(i) = rhoi(i) / (rhoi(i) + rhoiaux(i))
         end do

         ! load wi into structure
         suba(isub_loc)%lwi = lwi
         allocate(suba(isub_loc)%wi(lwi))
         do i = 1,lwi
            suba(isub_loc)%wi(i) = wi(i)
         end do
         suba(isub_loc)%is_weights_ready = .true.

         deallocate(wi)
         deallocate(rhoiaux)
         deallocate(rhoi)
      end do

end subroutine

!********************************************
subroutine dd_weightsi_apply(sub, veci,lveci)
!********************************************
! Subroutine for applying weight matrix on interface
      use module_utils
      implicit none
! Subdomain structure
      type(subdomain_type),intent(in) :: sub
      integer,intent(in) ::     lveci
      real(kr),intent(inout) ::  veci(lveci)

      ! local vars
      integer :: i
      integer :: ndofi, lwi

      if (.not.sub%is_weights_ready) then
         write(*,*) 'DD_WEIGHTSI_APPLY: Weights not ready for subdomain',sub%isub
         call error_exit
      end if
      ndofi = sub%ndofi
      ! check size
      if (lveci .ne. ndofi) then
         write(*,*) 'DD_WEIGHTSI_APPLY: Data size mismatch: lveci ',lveci, 'ndofi',ndofi
         call error_exit
      end if
      lwi = sub%lwi
      if (lveci .ne. lwi) then
         write(*,*) 'DD_WEIGHTSI_APPLY: Data size mismatch: lveci ',lveci, 'lwi',lwi
         call error_exit
      end if

      ! weight the vector
      do i = 1,ndofi
         veci(i) = sub%wi(i) * veci(i)
      end do

end subroutine

!*****************************************
subroutine dd_weights_apply(sub, vec,lvec)
!*****************************************
! Subroutine for applying weight matrix 

      use module_utils
      implicit none
! Subdomain structure
      type(subdomain_type),intent(in) :: sub
      integer,intent(in) ::     lvec
      real(kr),intent(inout) ::  vec(lvec)

      ! local vars
      integer :: i, indvs
      integer :: ndofi, lwi, ndof

      ! initial checking
      if (.not.sub%is_weights_ready) then
         write(*,*) 'DD_WEIGHTS_APPLY: Weights not ready for subdomain',sub%isub
         call error_exit
      end if
      if (.not. sub%is_mesh_loaded) then
         call error('DD_WEIGHTS_APPLY','mesh not loaded.')
      end if

      ndof  = sub%ndof
      ndofi = sub%ndofi

      ! check size
      if (lvec .ne. ndof) then
         write(*,*) 'DD_WEIGHTS_APPLY: Data size mismatch: lvec ',lvec, 'ndof',ndof
         call error_exit
      end if
      lwi = sub%lwi

      ! weight the vector at interface
      do i = 1,ndofi
         indvs = sub%iivsvn(i)
         vec(indvs) = sub%wi(i) * vec(indvs)
      end do

end subroutine

!****************************************************
subroutine dd_get_interface_diagonal(sub, rhoi,lrhoi)
!****************************************************
! Subroutine for getting subdomain diagonal from the structure
      use module_utils
      implicit none
! Subdomain structure
      type(subdomain_type),intent(in) :: sub
      integer,intent(in)  :: lrhoi
      real(kr),intent(out) ::  rhoi(lrhoi)

      ! local vars
      integer :: ndof, ndofi, la
      integer :: i, indiv, ia, ind

      integer ::              lrho
      real(kr),allocatable ::  rho(:)

      ! check if matrix is loaded
      if (.not. sub%is_matrix_loaded) then
         write(*,*) 'DD_GET_INTERFACE_DIAGONAL: Matrix not loaded for subdomain: ',sub%isub
         call error_exit
      end if
      ! check dimensions
      if (sub%ndofi.ne.lrhoi) then
         write(*,*) 'DD_GET_INTERFACE_DIAGONAL: Interface dimensions mismatch for subdomain ',sub%isub
         call error_exit
      end if

      ! load data
      ndof  = sub%ndof
      ndofi = sub%ndofi

      ! allocate vector for whole subdomain diagonal
      lrho = ndof
      allocate(rho(lrho))
      call zero(rho,lrho)

      la = sub%la
      do ia = 1,la
         if (sub%i_a_sparse(ia).eq.sub%j_a_sparse(ia)) then
            ind = sub%i_a_sparse(ia)

            rho(ind) = rho(ind) + sub%a_sparse(ia)
         end if
      end do

      do i = 1,ndofi
         indiv = sub%iivsvn(i)
         rhoi(i) = rho(indiv)
      end do

      deallocate(rho)
end subroutine

!****************************************************************
subroutine dd_dotprod_local(sub, vec1,lvec1, vec2,lvec2, dotprod)
!****************************************************************
! Subroutine for computing weighted dot product to be used in repeated entries with DD
! dotprod = vec1 * wi * vec2
      use module_utils
      implicit none
! Subdomain structure
      type(subdomain_type),intent(in) :: sub

      ! vectors to multiply
      integer,intent(in) ::  lvec1
      real(kr), intent(in) :: vec1(lvec1)
      integer,intent(in) ::  lvec2
      real(kr), intent(in) :: vec2(lvec2)
      
      ! result
      real(kr), intent(out) :: dotprod

      ! local vars
      integer :: i

      ! check the prerequisities
      if (.not. (sub%is_weights_ready)) then
         write(*,*) 'DD_DOTPROD_LOCAL: Weights not ready.'
         call error_exit
      end if

      ! check dimensions
      if (lvec1 .ne. lvec2) then
         write(*,*) 'DD_DOTPROD_LOCAL: Dimensions mismatch.'
         call error_exit
      end if
      if (lvec1 .ne. sub%lwi) then
         write(*,*) 'DD_DOTPROD_LOCAL: Dimensions mismatch with weights.'
         call error_exit
      end if

      dotprod = 0._kr
      do i = 1,lvec1
         dotprod = dotprod + vec1(i) * sub%wi(i) * vec2(i)
      end do

end subroutine
 
!***************************
subroutine dd_print_sub(sub)
!***************************
! Subroutine for printing the state of sub structure
      use module_sm
      implicit none
! Subdomain structure
      type(subdomain_type),intent(in) :: sub

! local variables
      integer :: i, j, ia, kishnadj

      write(*,*) '****** start subdomain export '
      write(*,*) '*** HEADER INFO :             '
      write(*,*) '     subdomain initialized:   ', sub%is_initialized
      write(*,*) '     global subdomain number: ', sub%isub
      write(*,*) '     number of subdomains:    ', sub%nsub
      write(*,*) '     processor number:        ', sub%proc
      write(*,*) '*** MESH INFO :               '
      write(*,*) '     mesh loaded:             ', sub%is_mesh_loaded
      write(*,*) '     number of elements:      ', sub%nelem
      write(*,*) '     number of nodes:         ', sub%nnod
      write(*,*) '     number of DOF:           ', sub%ndof
      write(*,*) '     number of dimensions:    ', sub%ndim
      write(*,*) '*** BOUNDARY CONDITIONS :     '
      write(*,*) '     is bc present:           ', sub%is_bc_present
      write(*,*) '     is bc nonzero:           ', sub%is_bc_nonzero
      write(*,*) '     is bc loaded:            ', sub%is_bc_loaded
      write(*,*) '*** CORNER INFO :             '
      write(*,*) '     are corners loaded:      ', sub%is_corners_loaded
      write(*,*) '     number of corners:       ', sub%ncorner
      write(*,*) '*** GLOB INFO :               '
      write(*,*) '     are globs loaded:        ', sub%is_globs_loaded
      write(*,*) '     number of globs:         ', sub%nglob
      write(*,*) '*** NEIGHBOURING INFO :       '
      write(*,*) '     interface ready:         ', sub%is_interface_loaded
      write(*,*) '     neighbouring ready:      ', sub%is_neighbouring_ready
      write(*,*) '     number of neighbours:    ', sub%nadj
      write(*,*) '     indices of neighbours:   ', sub%iadj
      if (sub%is_neighbouring_ready) then
         kishnadj = 0
         do ia = 1,sub%nadj
            write(*,*) '      number of nodes shared with subdomain ',sub%iadj(ia),'is: ', sub%nshnadj(ia)
            write(*,*) '      indices of nodes shared: ', sub%ishnadj(kishnadj+1:kishnadj+sub%nshnadj(ia))
            kishnadj = kishnadj + sub%nshnadj(ia)
         end do
      end if
      write(*,*) '*** WEIGHTS INFO :            '
      write(*,*) '     are weights ready?:      ', sub%is_weights_ready
      if (sub%is_weights_ready) then
         write(*,*) '     weights:                 ', sub%wi
      end if
      write(*,*) '*** COARSE NODES INFO :       '
      write(*,*) '     coarse nodes ready:      ', sub%is_cnodes_loaded
      write(*,*) '     number of coarse nodes:  ', sub%ncnodes
      do i = 1,sub%ncnodes
         call dd_print_cnode(sub,i)
      end do
      write(*,*) '*** MATRIX INFO :             '
      write(*,*) '     matrix loaded:           ', sub%is_matrix_loaded
      write(*,*) '     matrix blocked:          ', sub%is_blocked
      write(*,*) '     interior block factor.:  ', sub%is_interior_factorized
      if (debug) then
         if (sub%is_matrix_loaded.and.sub%is_triplet) then
            write(*,*) '     matrix data:           '
            call sm_print(6, sub%i_a_sparse, sub%j_a_sparse, sub%a_sparse, &
                          sub%la, sub%nnza)
         end if
         if (sub%is_matrix_loaded.and.sub%is_blocked) then
            write(*,*) '     matrix blocks:           '
            write(*,*) '     A_11:                     '
            call sm_print(6, sub%i_a11_sparse, sub%j_a11_sparse, sub%a11_sparse, &
                          sub%la11, sub%nnza11)
            write(*,*) '     A_12:                     '
            call sm_print(6, sub%i_a12_sparse, sub%j_a12_sparse, sub%a12_sparse, &
                          sub%la12, sub%nnza12)
            write(*,*) '     A_21:                     '
            call sm_print(6, sub%i_a21_sparse, sub%j_a21_sparse, sub%a21_sparse, &
                          sub%la21, sub%nnza21)
            write(*,*) '     A_22:                     '
            call sm_print(6, sub%i_a22_sparse, sub%j_a22_sparse, sub%a22_sparse, &
                          sub%la22, sub%nnza22)
         end if
      end if
      write(*,*) '*** BDDC INFO:                '
      write(*,*) '     matrix C loaded:         ', sub%is_c_loaded
      if (debug) then
         if (sub%is_c_loaded) then
            call sm_print(6, sub%i_c_sparse, sub%j_c_sparse, sub%c_sparse, &
                          sub%lc, sub%nnzc)
         end if
      end if
      write(*,*) '     matrix Kaug factorized:  ', sub%is_aug_factorized
      if (debug) then
         if (sub%is_matrix_loaded.and.sub%is_aug_factorized) then
            if (debug) then
               call sm_print(6, sub%i_aaug_sparse, sub%j_aaug_sparse, sub%aaug_sparse, &
                             sub%laaug, sub%nnzaaug)
            end if
         end if
      end if
      write(*,*) '     matrix PHIS prepared:    ', sub%is_phisi_prepared
      if (debug) then
         if (sub%is_coarse_prepared) then
            do i = 1,sub%lphisi1
               write(*,'(1000f13.6)') (sub%phisi(i,j),j = 1,sub%lphisi2)
            end do
         end if
      end if
      write(*,*) '     coarse matrix prepared:  ', sub%is_coarse_prepared
      if (debug) then
         if (sub%is_coarse_prepared) then
!      write(*,'(f13.6)') (sub%coarsem(j),j = 1,sub%lcoarsem)
            write(*,*) ' embedding of corse matrix :  '
            write(*,'(i8)') (sub%indrowc(j),j = 1,sub%lindrowc)
         end if
      end if
      ! PCG data
      write(*,*) '     reduced RHS loaded:  ', sub%is_reduced_rhs_loaded
      if (debug) then
         if (sub%is_reduced_rhs_loaded) then
!      write(*,'(f13.6)') (sub%coarsem(j),j = 1,sub%lcoarsem)
            write(*,*) ' reduced RHS :  '
            write(*,'(e15.5)') (sub%g(j),j = 1,sub%lg)
         end if
      end if
end subroutine

!************************************
subroutine dd_print_cnode(sub,icnode)
!************************************
! Subroutine for printing content of one coarse node
      implicit none
! Subdomain structure
      type(subdomain_type),intent(in) :: sub
      integer,intent(in) :: icnode

! basic structure
      write(*,*) '****** start coarse node export '
      write(*,*) '     coarse node number:      ', icnode
      write(*,*) '     type of coarse node:     ', sub%cnodes(icnode)%itype
      write(*,*) '     used for constraints?:   ', sub%cnodes(icnode)%used
      write(*,*) '     coordinates:             ', sub%cnodes(icnode)%xyz
      write(*,*) '****** where it maps to? '
      write(*,*) '     global coarse node number:', sub%cnodes(icnode)%global_cnode_number
      write(*,*) '     number of coarse degrees of freedom:', sub%cnodes(icnode)%ncdof
      write(*,*) '     indices of coarse dof:', sub%cnodes(icnode)%igcdof
!      write(*,*) '****** where it maps from? '
!      write(*,*) '     number of nodes it contains:', sub%cnodes(icnode)%nnod
!      write(*,*) '     indices of nodes on subdomain int:', sub%cnodes(icnode)%insin
      write(*,*) '     number of nodes it contains:', sub%cnodes(icnode)%nvar
      write(*,*) '     indices of variables on subdomain int:', sub%cnodes(icnode)%ivsivn
      write(*,*) '****** end coarse nodes export '
end subroutine

!**************************
subroutine dd_finalize(sub)
!**************************
! Subroutine for deallocation of data of subdomain
      use module_mumps
      implicit none
! Subdomain structure
      type(subdomain_type),intent(inout) :: sub

      ! local variables
      integer :: j

      if (allocated(sub%inet)) then
         deallocate(sub%inet)
      end if
      if (allocated(sub%nnet)) then
         deallocate(sub%nnet)
      end if
      if (allocated(sub%nndf)) then
         deallocate(sub%nndf)
      end if
      if (allocated(sub%isngn)) then
         deallocate(sub%isngn)
      end if
      if (allocated(sub%isvgvn)) then
         deallocate(sub%isvgvn)
      end if
      if (allocated(sub%isegn)) then
         deallocate(sub%isegn)
      end if
      if (allocated(sub%xyz)) then
         deallocate(sub%xyz)
      end if
      if (allocated(sub%iin)) then
         deallocate(sub%iin)
      end if
      if (allocated(sub%iivsvn)) then
         deallocate(sub%iivsvn)
      end if
      if (allocated(sub%iovsvn)) then
         deallocate(sub%iovsvn)
      end if
      if (allocated(sub%ifix)) then
         deallocate(sub%ifix)
      end if
      if (allocated(sub%fixv)) then
         deallocate(sub%fixv)
      end if
      if (allocated(sub%bc)) then
         deallocate(sub%bc)
      end if
      ! corners and globs
      if (allocated(sub%global_corner_number)) then
         deallocate(sub%global_corner_number)
      end if
      if (allocated(sub%icnsin)) then
         deallocate(sub%icnsin)
      end if
      if (allocated(sub%global_glob_number)) then
         deallocate(sub%global_glob_number)
      end if
      if (allocated(sub%nglobnodes)) then
         deallocate(sub%nglobnodes)
      end if
      if (allocated(sub%ignsin)) then
         deallocate(sub%ignsin)
      end if
      if (allocated(sub%nglobvar)) then
         deallocate(sub%nglobvar)
      end if
      if (allocated(sub%igvsivn)) then
         deallocate(sub%igvsivn)
      end if
      if (allocated(sub%glob_type)) then
         deallocate(sub%glob_type)
      end if
      if (allocated(sub%iadj)) then
         deallocate(sub%iadj)
      end if
      if (allocated(sub%nshnadj)) then
         deallocate(sub%nshnadj)
      end if
      if (allocated(sub%kshvadj)) then
         deallocate(sub%kshvadj)
      end if
      if (allocated(sub%ishnadj)) then
         deallocate(sub%ishnadj)
      end if
      if (allocated(sub%commvec_out)) then
         deallocate(sub%commvec_out)
      end if
      if (allocated(sub%commvec_in)) then
         deallocate(sub%commvec_in)
      end if
      if (allocated(sub%wi)) then
         deallocate(sub%wi)
      end if
      if (allocated(sub%cnodes)) then
         do j = 1,sub%ncnodes
            if (allocated(sub%cnodes(j)%xyz)) then
               deallocate(sub%cnodes(j)%xyz)
            end if
            if (allocated(sub%cnodes(j)%igcdof)) then
               deallocate(sub%cnodes(j)%igcdof)
            end if
            if (allocated(sub%cnodes(j)%insin)) then
               deallocate(sub%cnodes(j)%insin)
            end if
            if (allocated(sub%cnodes(j)%ivsivn)) then
               deallocate(sub%cnodes(j)%ivsivn)
            end if
            if (allocated(sub%cnodes(j)%matrix)) then
               deallocate(sub%cnodes(j)%matrix)
            end if
         end do
         deallocate(sub%cnodes)
      end if
      if (allocated(sub%i_a_sparse)) then
         deallocate(sub%i_a_sparse)
      end if
      if (allocated(sub%j_a_sparse)) then
         deallocate(sub%j_a_sparse)
      end if
      if (allocated(sub%a_sparse)) then
         deallocate(sub%a_sparse)
      end if
      if (allocated(sub%i_a11_sparse)) then
         deallocate(sub%i_a11_sparse)
      end if
      if (allocated(sub%j_a11_sparse)) then
         deallocate(sub%j_a11_sparse)
      end if
      if (allocated(sub%a11_sparse)) then
         deallocate(sub%a11_sparse)
      end if
      if (allocated(sub%i_a12_sparse)) then
         deallocate(sub%i_a12_sparse)
      end if
      if (allocated(sub%j_a12_sparse)) then
         deallocate(sub%j_a12_sparse)
      end if
      if (allocated(sub%a12_sparse)) then
         deallocate(sub%a12_sparse)
      end if
      if (allocated(sub%i_a21_sparse)) then
         deallocate(sub%i_a21_sparse)
      end if
      if (allocated(sub%j_a21_sparse)) then
         deallocate(sub%j_a21_sparse)
      end if
      if (allocated(sub%a21_sparse)) then
         deallocate(sub%a21_sparse)
      end if
      if (allocated(sub%i_a22_sparse)) then
         deallocate(sub%i_a22_sparse)
      end if
      if (allocated(sub%j_a22_sparse)) then
         deallocate(sub%j_a22_sparse)
      end if
      if (allocated(sub%a22_sparse)) then
         deallocate(sub%a22_sparse)
      end if
      if (sub%is_mumps_interior_active) then
         call mumps_finalize(sub%mumps_interior_block)
      end if

      ! BDDC matrices
      if (allocated(sub%i_c_sparse)) then
         deallocate(sub%i_c_sparse)
      end if
      if (allocated(sub%j_c_sparse)) then
         deallocate(sub%j_c_sparse)
      end if
      if (allocated(sub%c_sparse)) then
         deallocate(sub%c_sparse)
      end if
      if (allocated(sub%indrowc)) then
         deallocate(sub%indrowc)
      end if
      if (allocated(sub%i_aaug_sparse)) then
         deallocate(sub%i_aaug_sparse)
      end if
      if (allocated(sub%j_aaug_sparse)) then
         deallocate(sub%j_aaug_sparse)
      end if
      if (allocated(sub%aaug_sparse)) then
         deallocate(sub%aaug_sparse)
      end if
      if (allocated(sub%phis)) then
         deallocate(sub%phis)
      end if
      if (allocated(sub%phisi)) then
         deallocate(sub%phisi)
      end if
      if (allocated(sub%coarsem)) then
         deallocate(sub%coarsem)
      end if
      if (allocated(sub%rhs)) then
         deallocate(sub%rhs)
      end if
      if (allocated(sub%g)) then
         deallocate(sub%g)
      end if

      if (sub%is_mumps_aug_active) then
         call mumps_finalize(sub%mumps_aug)
      end if

!     Krylov vectors on subdomain
      if (allocated(sub%sol)) then
         deallocate(sub%sol)
      end if
      if (allocated(sub%soli)) then
         deallocate(sub%soli)
      end if

      sub%is_mesh_loaded          = .false.
      sub%is_interface_loaded     = .false.
      sub%is_adj_loaded           = .false.
      sub%is_corners_loaded       = .false.
      sub%is_globs_loaded         = .false.
      sub%is_neighbouring_ready   = .false.
      sub%is_matrix_loaded        = .false.
      sub%is_c_loaded             = .false.
      sub%is_phis_prepared        = .false.
      sub%is_phisi_prepared       = .false.
      sub%is_coarse_prepared      = .false.
      sub%is_aug_factorized       = .false.
      sub%is_interior_factorized  = .false.
      sub%is_weights_ready        = .false.

end subroutine

end module module_dd

