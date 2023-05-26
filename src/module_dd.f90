! BDDCML - Multilevel BDDC
! 
! This program is a free software.
! You can redistribute it and/or modify it under the terms of 
! the GNU Lesser General Public License 
! as published by the Free Software Foundation, 
! either version 3 of the license, 
! or (at your option) any later version.
! 
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU Lesser General Public License for more details
! <http://www.gnu.org/copyleft/lesser.html>.
!________________________________________________________________

module module_dd
!***************
! Module for handling domain decomposition structures
! Jakub Sistek, Denver, 4/2009, Praha 1/2010, Cambridge 5/2011

      use, intrinsic :: iso_fortran_env
      implicit none
!     definition of MUMPS structure
      include "dmumps_struc.h"

! type of real variables
      integer,parameter,private :: kr = REAL64
! numerical zero
      real(kr),parameter,private :: numerical_zero = 1.e-12_kr
! maximal allowed length of file names
      integer,parameter,private :: lfnamex = 130

! debugging mode
      logical,parameter,private :: debug = .false.

      interface dd_map_glob_to_sub
         module procedure dd_map_glob_to_sub_int
         module procedure dd_map_glob_to_sub_real
      end interface dd_map_glob_to_sub

! type for storing adaptive constraints for selected globs
      type cnode_type
         integer :: itype                ! type of coarse node (1 - face, 2 - edge, 3 - corner)
         logical :: used       = .false. ! should this glob be applied in computation?
         integer ::             lxyz
         real(kr),allocatable :: xyz(:) ! coordinates
         ! where the glob maps to
         integer :: global_cnode_number ! embedding into coarse nodes
         integer :: ncdof = 0              ! number of coarse degrees of freedom
         integer,allocatable :: igcdof(:) ! indices of global coarse dof
         ! where it maps from
         integer :: nnod  = 0             ! number of nodes it contains
         integer,allocatable :: insin(:)  ! indices of glob nodes in subdomain interface numbering
         integer :: nvar = 0              ! number of variables it contains
         integer,allocatable :: ivsivn(:)  ! indices of glob variables in subdomain interface numbering
         ! constraints on glob
         integer :: nnz      = 0 ! number of nonzeros the matrix contains
         integer :: lmatrix1 = 0 ! = ncdof
         integer :: lmatrix2 = 0 ! = nvar
         real(kr),allocatable :: matrix(:,:)
      end type cnode_type

! type for subdomain data
      type subdomain_type

         logical ::             is_initialized = .false.
         integer ::             isub    ! subdomain index
         integer ::             nsub    ! total number of subdomains
         integer ::             proc    ! index of processor taking care of this subdomain
         integer ::             comm    ! communicator to which proc refer

         logical ::             is_degenerated = .false. ! this flag serves for switching off degenerated subdomains 

         logical ::             is_mesh_loaded = .false.
         integer ::             nelem   ! number of elements
         integer ::             nnod    ! number of nodes
         integer ::             ndof    ! number of degrees of freedom
         integer ::             ndim    ! dimension of the problem
         integer ::             meshdim ! dimension of the mesh - e.g. shells in 3D have 2, beams have 1, etc.
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
         ! components of the nodal graph
         logical ::             find_components     = .false.  
         logical ::             use_dual_mesh_graph = .false. ! use dual mesh graph for detecting components?
         integer ::             neighbouring = 1              ! how many nodes need element to share to create an edge in graph?
         logical ::             is_nodal_components_loaded = .false.
         integer ::             nnodal_components    ! number of local subcomponents of the mesh - 1 for contiguous domain
         integer ::             lnodal_components    ! lenght of NODAL_COMPONENTS array - number of nodes
         integer,allocatable ::  nodal_components(:) ! array of indices of local subcomponents for disconnected subdomains
      
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
         logical ::             is_bc_present = .false. ! are some Dirichlet BC on subdomain?
         logical ::             is_bc_nonzero = .false. ! are some Dirichlet BC nonzero?
         logical ::             is_bc_loaded = .false. 
         integer ::             lifix         ! length of IFIX array
         integer,allocatable ::  ifix(:)      ! IFIX array - indices of fixed variables
         integer ::         lfixv         ! length of FIXV array
         real(kr),allocatable :: fixv(:)      ! FIXV array - fixed variables values
         integer ::         lbc = 0       ! length of BC array
         real(kr),allocatable :: bc(:)        ! BC array - eliminated entries of stiffness matrix multiplied by values of fixed variables 

         ! matrix of fixed rows in IJA sparse format
         logical :: is_fixed_rows_matrix_loaded = .false.
         integer :: nnza_fixed
         integer :: la_fixed
         integer,allocatable  :: i_a_fixed_sparse(:)
         integer,allocatable  :: j_a_fixed_sparse(:)
         real(kr),allocatable ::   a_fixed_sparse(:)

         ! matrix of user defined constraints
         logical :: is_user_constraints_loaded = .false.
         integer ::              luser_constraints1 ! number of constraints
         integer ::              luser_constraints2 ! number of columns (=NNODS)
         real(kr),allocatable ::  user_constraints(:,:) ! array of user defined constraints

         ! matrix of element data
         logical :: is_element_data_loaded = .false.
         integer ::              lelement_data1 ! number of additional data per element 
         integer ::              lelement_data2 ! number of columns (=NELEMS)
         real(kr),allocatable ::  element_data(:,:) ! array of element data (e.g. constants on elements, 
                                                    ! possibly used for construction of weights)

         ! array of dof data
         logical :: is_dof_data_loaded = .false.
         integer ::              ldof_data ! length (=NDOFS)
         real(kr),allocatable ::  dof_data(:) ! array of dof data (e.g. weights for the weight matrix )

         ! reactions at fixed variables
         logical :: is_reactions_ready = .false.
         integer ::         lrea = 0       ! length of REA array
         real(kr),allocatable :: rea(:)    ! REA array - reactions = A*sol - rhs

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
         integer ::             lnsubglobs             ! length of array NSUBGLOBS
         integer, allocatable :: nsubglobs(:)          ! number of subdomains in glob (only neighbours)
         integer ::             lglob_subs1            ! number of rows of GLOB_SUBS array
         integer ::             lglob_subs2            ! number of cols of GLOB_SUBS array
         integer, allocatable :: glob_subs(:,:)        ! indices of subdomains in glob (only neighbours)

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
         integer ::             lishnncadj             ! length of array ISHNNCADJ
         integer, allocatable :: ishnncadj(:)          ! indices of nodal components of nodes shared with adjacent subdomains 
                                                       ! (in numbering of the neighbour)
         integer ::              lcommvec
         real(kr),allocatable ::  commvec_out(:)       ! communication vector for sending data
         real(kr),allocatable ::  commvec_in(:)        ! communication vector for receiving data

         ! weights on interface
         integer ::              lwi
         real(kr),allocatable ::  wi(:)                ! weights at interface
         integer ::               weights_type         ! type of weights:
                                                       ! 0 - weights by cardinality
                                                       ! 1 - weights by diagonal stiffness
                                                       ! 2 - weights based on first row of element data
                                                       ! 3 - weights based on dof data
                                                       ! 4 - weights by Marta Certikova - unit load
                                                       ! 5 - weights by Marta Certikova - unit jump
                                                       ! 6 - weights by Schur row sums for whole subdomain
                                                       ! 7 - weights by Schur row sums computed face by face
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
         integer ::            mumps_interior_block_factor_size

         logical :: is_subdomain_matrix_factorized = .false.
         type(DMUMPS_STRUC) :: mumps_subdomain_matrix

         ! explicit Schur complement if desired
         integer ::             lschur1
         integer ::             lschur2
         real(kr), allocatable :: schur(:,:)
         logical :: is_explicit_schur_prepared = .false.

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

         integer ::             laaug_dense1
         integer ::             laaug_dense2
         real(kr), allocatable :: aaug_dense(:,:)
         integer ::              laaug_ipiv
         integer, allocatable ::  aaug_ipiv(:)

         logical :: is_aug_dense_active = .false.
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

         ! dual coarse space basis functions on whole subdomain PHIS
         logical :: is_phis_dual_prepared  = .false.
         integer ::               lphis_dual1
         integer ::               lphis_dual2
         real(kr),allocatable ::   phis_dual(:,:)

         ! dual coarse space basis functions restricted to interface PHISI
         logical :: is_phisi_dual_prepared  = .false.
         integer ::               lphisi_dual1
         integer ::               lphisi_dual2
         real(kr),allocatable ::   phisi_dual(:,:)

         ! subdomain coarse matrix 
         logical :: is_coarse_prepared = .false.
         integer ::               coarse_matrixtype
         integer ::               lcoarsem = 0
         real(kr),allocatable ::   coarsem(:)
         ! embedded into global coarse matrix by array INDROWC

         ! arrays connected to iterative methods
         logical :: is_rhs_loaded = .false.
         logical :: is_rhs_complete = .true. ! if TRUE, the solver weights the RHS before using it to manage repeated entries
         integer ::             lrhs         ! length of RHS array
         real(kr),allocatable :: rhs(:)      ! RHS array - right hand side restricted to subdomain
         logical :: is_reduced_rhs_loaded = .false.
         integer ::             lg           ! length of G array
         real(kr),allocatable :: g(:)        ! condensed right hand side on subdomain 
         logical :: is_interior_solution_loaded = .false.
         integer ::             lsolo
         real(kr),allocatable :: solo(:)     ! array of solution at interior of subdomain

         logical :: is_solution_loaded = .false.
         integer ::             lsol         
         real(kr),allocatable :: sol(:)      ! array of solution restricted to subdomain

      end type subdomain_type

      type sub_aux_type

         integer ::         lisngn
         integer,pointer ::  isngn(:)
         integer,pointer ::  isnc(:)
         integer ::             lnnodadj
         integer,allocatable ::  nnodadj(:)
         integer ::             lisngnadj
         integer,allocatable ::  isngnadj(:)
         integer,allocatable ::  isncadj(:)

         integer ::             nglobs
         integer ::             ncorners
         integer ::             nedges
         integer ::             nfaces
         integer ::             lnsubnode
         integer,allocatable ::  nsubnode(:)
         integer ::             lglobsubs1
         integer ::             lglobsubs2
         integer ::             lglobsubs3
         integer,allocatable ::  globsubs(:,:,:)
         integer ::             lglobtypes
         integer,allocatable ::  globtypes(:)
         integer ::             lkglobs
         integer,allocatable ::  kglobs(:)
         integer ::             licgn
         integer,allocatable ::  icgn(:)
         integer ::             licgcn
         integer,allocatable ::  icgcn(:)
         integer ::             liegn
         integer,allocatable ::  iegn(:)
         integer ::             liegcn
         integer,allocatable ::  iegcn(:)
         integer ::             lifgn
         integer,allocatable ::  ifgn(:)
         integer ::             lifgcn
         integer,allocatable ::  ifgcn(:)

         integer ::            lineibcorners
         integer,allocatable :: ineibcorners(:)

         integer,pointer     :: global_indices(:)
         integer,pointer     :: global_coarse_indices(:)

         integer,pointer     :: comm_array_numbers(:)
         integer,allocatable :: comm_gathered_array1(:)
         integer,allocatable :: comm_gathered_array2(:)
         integer,pointer     :: comm_array_result1(:)
         integer,pointer     :: comm_array_result2(:)

         integer,pointer     :: comm_array_in(:)
         integer,pointer     :: comm_array_out(:)

      end type

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

!!*************************************************
!subroutine dd_read_mesh_from_file(sub,problemname)
!!*************************************************
!! Subroutine for reading of subdomain data from separate SMD file
!      use module_utils
!      implicit none
!
!! Subdomain structure
!      type(subdomain_type),intent(inout) :: sub
!! name of problem
!      character(*),intent(in) :: problemname
!
!! local variables
!      integer :: isub
!      integer :: idsmd, iglob, ivar, ignode
!      character(lfnamex):: fname
!
!      integer :: nelem, nnod, ndof, ndim, ncorner, nglob, nadj
!      integer ::             lnndf,   lnnet,   linet 
!      integer, allocatable :: nndf(:), nnet(:), inet(:)
!      integer ::             lisngn
!      integer, allocatable :: isngn(:)
!      integer ::             lisvgvn
!      integer, allocatable :: isvgvn(:)
!      integer ::             lisegn
!      integer, allocatable :: isegn(:)
!      integer ::              lxyz1, lxyz2
!      real(kr), allocatable :: xyz(:,:)
!      integer :: nnodi, ndofi, ndofo
!      integer ::             liin,   liivsvn,   liovsvn
!      integer, allocatable :: iin(:), iivsvn(:), iovsvn(:)
!      integer ::             lifix
!      integer, allocatable :: ifix(:)
!      integer ::              lfixv
!      real(kr), allocatable :: fixv(:)
!      integer ::             liadj
!      integer, allocatable :: iadj(:)
!      integer ::              lglobal_corner_number
!      integer, allocatable ::  global_corner_number(:)
!      integer ::              licnsin 
!      integer, allocatable ::  icnsin(:)
!      integer ::              lglobal_glob_number
!      integer, allocatable ::  global_glob_number(:)
!      integer ::              lglob_type
!      integer, allocatable ::  glob_type(:)
!      integer ::              lnglobnodes 
!      integer, allocatable ::  nglobnodes(:)
!      integer ::              lignsin1, lignsin2
!      integer, allocatable ::  ignsin(:,:)
!      integer ::              lnglobvar 
!      integer, allocatable ::  nglobvar(:)
!      integer ::              ligvsivn1, ligvsivn2
!      integer, allocatable ::  igvsivn(:,:)
!      integer ::              lrhss
!      real(kr), allocatable :: rhss(:)
!
!      if (.not.sub%is_initialized) then
!         call error('DD_READ_MESH_FROM_FILE','subdomain not initialized')
!      end if
!
!      isub = sub%isub
!
!      ! open subdomain SMD file with mesh description
!      call getfname(problemname,isub,'SMD',fname)
!      if (debug) then
!         write(*,*) 'DD_READ_MESH_FROM_FILE: Opening file fname: ', trim(fname)
!      end if
!      call allocate_unit(idsmd)
!      open (unit=idsmd,file=trim(fname),status='old',form='formatted')
!
!      ! read data from file
!      ! ---header
!      read(idsmd,*) nelem, nnod, ndof, ndim
!
!      ! ---NNDF array
!      ! ---NNET array
!      lnndf = nnod
!      lnnet = nelem
!      allocate (nndf(lnndf),nnet(lnnet))
!      read(idsmd,*) nndf
!      read(idsmd,*) nnet
!
!      ! ---INET array
!      linet = sum(nnet)
!      allocate (inet(linet))
!      read(idsmd,*) inet
!
!      ! ---ISNGN array
!      lisngn = nnod
!      allocate (isngn(lisngn))
!      read(idsmd,*) isngn
!
!      ! ---ISVNGVN array
!      lisvgvn = ndof
!      allocate (isvgvn(lisvgvn))
!      call zero(isvgvn,lisvgvn)
!
!      ! ---ISEGN array
!      lisegn = nelem
!      allocate (isegn(lisegn))
!      read(idsmd,*) isegn
!
!      ! --- coordinates
!      lxyz1 = nnod
!      lxyz2 = ndim
!      allocate(xyz(lxyz1,lxyz2))
!      read(idsmd,*) xyz
!
!      ! ---interface data
!      read(idsmd,*) nnodi
!      read(idsmd,*) ndofi
!      read(idsmd,*) ndofo
!
!      liin    = nnodi
!      liivsvn = ndofi
!      liovsvn = ndofo
!      allocate (iin(liin), iivsvn(liivsvn), iovsvn(liovsvn))
!      read(idsmd,*) iin
!      read(idsmd,*) iivsvn
!      read(idsmd,*) iovsvn
!
!      ! --- boundary conditions
!      lifix = ndof
!      lfixv = ndof
!      allocate(ifix(lifix),fixv(lfixv))
!      read(idsmd,*) ifix
!      read(idsmd,*) fixv
!
!      ! --- neighbouring
!      read(idsmd,*) nadj
!      liadj = nadj
!      allocate(iadj(liadj))
!      read(idsmd,*) iadj
!
!      ! --- corners 
!      read(idsmd,*) ncorner
!      lglobal_corner_number = ncorner
!      licnsin               = ncorner
!      allocate(global_corner_number(lglobal_corner_number),icnsin(licnsin))
!      read(idsmd,*) global_corner_number
!      read(idsmd,*) icnsin
!
!      ! --- globs
!      read(idsmd,*) nglob
!      lglobal_glob_number = nglob
!      lnglobvar           = nglob
!      lnglobnodes         = nglob
!      allocate(global_glob_number(lglobal_glob_number),nglobvar(lnglobvar),nglobnodes(lnglobnodes))
!      read(idsmd,*) global_glob_number
!
!      read(idsmd,*) nglobnodes
!      lignsin1 = nglob
!      lignsin2 = maxval(nglobnodes)
!      allocate(ignsin(lignsin1,lignsin2))
!      do iglob = 1,nglob
!         ! read glob variables
!         read(idsmd,*) (ignsin(iglob,ignode),ignode = 1,nglobnodes(iglob))
!         ! pad the glob with zeros
!         do ignode = nglobnodes(iglob) + 1,lignsin2
!            ignsin(iglob,ignode) = 0
!         end do
!      end do
!      read(idsmd,*) nglobvar
!      ligvsivn1 = nglob
!      ligvsivn2 = maxval(nglobvar)
!      allocate(igvsivn(ligvsivn1,ligvsivn2))
!      do iglob = 1,nglob
!         ! read glob variables
!         read(idsmd,*) (igvsivn(iglob,ivar),ivar = 1,nglobvar(iglob))
!         ! pad the glob with zeros
!         do ivar = nglobvar(iglob) + 1,ligvsivn2
!            igvsivn(iglob,ivar) = 0
!         end do
!      end do
!      lglob_type = nglob
!      allocate(glob_type(lglob_type))
!      read(idsmd,*) glob_type
!
!      ! --- right hand side
!      lrhss = ndof
!      allocate(rhss(lrhss))
!      read(idsmd,*) rhss
!
!      close(idsmd)
!      if (debug) then
!         write(*,*) 'DD_READ_MESH_FROM_FILE: Data read successfully.'
!      end if
!
!      ! load data to structure
!      call dd_upload_sub_mesh(sub, nelem, nnod, ndof, ndim, &
!                              nndf,lnndf, nnet,lnnet, 0, inet,linet, &
!                              isngn,lisngn, isvgvn,lisvgvn, isegn,lisegn,&
!                              xyz,lxyz1,lxyz2)
!      call dd_upload_sub_adj(sub, nadj, iadj,liadj)
!      !call dd_upload_sub_interface(sub, nnodi, ndofi, ndofo, &
!      !                             iin,liin, iivsvn,liivsvn, iovsvn,liovsvn)
!      call dd_upload_sub_corners(sub, ncorner,&
!                                 global_corner_number,lglobal_corner_number, icnsin,licnsin)
!      call dd_upload_sub_globs(sub, nglob, &
!                               global_glob_number,lglobal_glob_number, &
!                               nglobnodes,lnglobnodes, nglobvar,lnglobvar,&
!                               ignsin,lignsin1,lignsin2, igvsivn,ligvsivn1,ligvsivn2,&
!                               glob_type,lglob_type)
!      call dd_upload_bc(sub, ifix,lifix, fixv,lfixv)
!      call dd_upload_rhs(sub, rhss,lrhss, .true.)
!
!      if (debug) then
!         write(*,*) 'DD_READ_MESH_FROM_FILE: Data loaded successfully.'
!      end if
!
!      deallocate (nndf,nnet)
!      deallocate (inet)
!      deallocate (isngn)
!      deallocate (isegn)
!      deallocate (xyz)
!      deallocate (iin, iivsvn, iovsvn)
!      deallocate (iadj)
!      deallocate (ifix,fixv)
!      deallocate (glob_type)
!      deallocate (global_glob_number,nglobvar,nglobnodes)
!      deallocate (ignsin)
!      deallocate (igvsivn)
!      deallocate (global_corner_number,icnsin)
!      deallocate (rhss)
!
!end subroutine

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
         call sm_pmd_load(matrixtype,idelm,nelem,sub%inet,sub%linet,&
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
         call sm_pmd_load_masked(matrixtype,idelm,nelem,sub%inet,sub%linet,&
                                 sub%nnet,sub%lnnet,sub%nndf,sub%lnndf,&
                                 kdof,lkdof,sub%isegn,sub%lisegn,&
                                 i_sparse, j_sparse, a_sparse, la)
         close (idelm)
      end if
      if (debug) then
         write(*,*) 'DD_READ_MATRIX_FROM_FILE: File with element matrices read.'
      end if

      nnza = la

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
         call sm_apply_bc(matrixtype,sub%ifix,sub%lifix,sub%fixv,sub%lfixv,&
                          nnza, i_sparse,j_sparse,a_sparse,la, bc,lbc, .false., 0)
         if (sub%is_bc_nonzero) then
            call dd_load_eliminated_bc(sub, bc,lbc)
         else
            sub%lbc = 0
         end if
      end if

      ! load matrix to our structure
      is_assembled = .false.
      call dd_load_matrix_triplet(sub, matrixtype, 0, &
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
subroutine dd_read_matrix_by_root(suba,lsuba, comm_all,idelm,nsub,nelem,matrixtype,&
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
! Fortran unit with attached file with element matrices
      integer,intent(in) :: idelm
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
      character(*),parameter:: routine_name = 'DD_READ_MATRIX_BY_ROOT'
      integer :: myid, nproc, ierr
      integer :: stat(MPI_STATUS_SIZE)
      integer :: nelems, nnods
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

      integer :: i, ie, iproc, inods, isub, isub_loc 
      logical        :: is_opened 
      character :: is_readable


! check dimension
      if (nelem.ne.liets) then
         call error(routine_name,'wrong length of array IETS')
      end if

! orient in communicator
      call MPI_COMM_RANK(comm_all,myid,ierr)
      call MPI_COMM_SIZE(comm_all,nproc,ierr)

      if (debug) then
         if (nproc+1.ne.lsub2proc) then
            call error(routine_name,'wrong length of array sub2proc')
         end if
         if (any(suba%comm .ne. comm_all)) then
            call error(routine_name,'communicator mismatch')
         end if
         if (any(suba%proc .ne. myid)) then
            call error(routine_name,'communicator mismatch')
         end if
         if (any(suba%nsub .ne. nsub)) then
            call error(routine_name,'number of subdomains mismatch')
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
            call error(routine_name,'strange type of matrix:',matrixtype)
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
         call sm_pmd_make_element_numbering(matrixtype,nelems,suba(isub_loc)%inet,suba(isub_loc)%linet,&
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

      ! check that ELM file is opened on root processor
      if (myid.eq.0) then
         ! ELM - element stiffness matrices - structure:
         inquire (idelm, opened=is_opened, read=is_readable )
         if (.not. is_opened .or. trim(is_readable) .eq. 'N' ) then
            write (*,*) is_opened, is_readable
            call error(routine_name,'Cannot access file with element matrices, unit IDELM.')
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
      
      ! set proper flags after matrix is loaded
      do isub_loc = 1,lindexsub

         suba(isub_loc)%is_matrix_loaded = .true.
         suba(isub_loc)%is_triplet       = .true.
         suba(isub_loc)%is_assembled     = .false.

      end do

      ! free memory
      deallocate(elm)
      deallocate(locsubnumber)
      deallocate(subp)

end subroutine

!*****************************************************************************************
subroutine dd_fix_constraints(suba,lsuba, comm_all, sub2proc,lsub2proc,indexsub,lindexsub)
!*****************************************************************************************
! Subroutine for eliminating constraints from RHS and fixing them in the matrix
      use module_sm
      use module_utils
      implicit none
      include "mpif.h"

! array of sub structure
      integer,intent(in) ::                lsuba
      type(subdomain_type),intent(inout) :: suba(lsuba)
! communicator
      integer,intent(in) :: comm_all
! division of subdomains to processors
      integer,intent(in) :: lsub2proc
      integer,intent(in) ::  sub2proc(lsub2proc)
! global indices of local subdomains
      integer,intent(in) :: lindexsub
      integer,intent(in) ::  indexsub(lindexsub)

! local variables
      character(*),parameter:: routine_name = 'DD_FIX_CONSTRAINTS'
      integer :: ierr
      logical :: is_bc_present_loc, is_bc_present

!     BC arrays
      integer::              lbc
      real(kr),allocatable :: bc(:)
      integer::              lbci
      real(kr),allocatable :: bci(:)

      integer :: isub, isub_loc, nnza, la, ndofs, ndofis, ia, irow
      real(kr) :: aval, fixedval

      logical,parameter :: store_fixed_rows = .true.
      integer :: la_fixed

      ! find if any nonzero BC is present
      is_bc_present_loc = .false.
      do isub_loc = 1,lindexsub
         if (suba(isub_loc)%is_bc_present) then
            is_bc_present_loc = .true.
         end if
      end do
!*****************************************************************MPI
      call MPI_ALLREDUCE(is_bc_present_loc,is_bc_present,1, MPI_LOGICAL, MPI_LOR, comm_all, ierr) 
!*****************************************************************MPI

      ! eliminate constraints if they are present
      do isub_loc = 1,lindexsub
         isub = indexsub(isub_loc)

         nnza = suba(isub_loc)%nnza
         la   = suba(isub_loc)%la

         suba(isub_loc)%is_assembled = .false.

         ! eliminate boundary conditions
         ndofs =  suba(isub_loc)%ndof

         if (is_bc_present) then
            lbc = ndofs
            allocate(bc(lbc))
            bc(:) = 0._kr
         end if

         if (suba(isub_loc)%is_fixed_rows_matrix_loaded) then
            ! already eliminated rows of the matrix, just recreate bc by
            ! multiplication of such matrix entries with prescribed values of solution
            if (suba(isub_loc)%is_bc_present) then
               call sm_vec_mult(suba(isub_loc)%matrixtype, suba(isub_loc)%nnza_fixed, &
                                suba(isub_loc)%i_a_fixed_sparse, &
                                suba(isub_loc)%j_a_fixed_sparse, &
                                suba(isub_loc)%a_fixed_sparse, & 
                                suba(isub_loc)%la_fixed, &
                                suba(isub_loc)%fixv,suba(isub_loc)%lfixv, &
                                bc,lbc)
               ! flip sign
               bc = -bc

               ! correct entries on diagonal
               do ia = 1,suba(isub_loc)%nnza
                  if (suba(isub_loc)%i_a_sparse(ia).eq.suba(isub_loc)%j_a_sparse(ia)) then
                     irow  = suba(isub_loc)%i_a_sparse(ia)
                     if (size(suba(isub_loc)%ifix).gt.0) then
                        if (suba(isub_loc)%ifix(irow).gt.0) then
                           aval     = suba(isub_loc)%a_sparse(ia)
                           fixedval = suba(isub_loc)%fixv(irow)

                           bc(irow) = aval*fixedval
                        end if
                     end if
                  end if
               end do
            end if
         else
            ! first elimination of the matrix, store the eliminated entries of the matrix
            ! eliminate natural BC
            la_fixed = 0
            if (suba(isub_loc)%is_bc_present) then
               la_fixed = 0
               if ( store_fixed_rows ) then
                   call sm_count_fixed_rows_matrix_size(suba(isub_loc)%matrixtype,suba(isub_loc)%ifix,suba(isub_loc)%lifix,nnza,&
                                                        suba(isub_loc)%i_a_sparse,suba(isub_loc)%j_a_sparse, la, la_fixed)
               end if
            end if
            suba(isub_loc)%la_fixed   = la_fixed
            suba(isub_loc)%nnza_fixed = la_fixed
            allocate( suba(isub_loc)%i_a_fixed_sparse(la_fixed), &
                      suba(isub_loc)%j_a_fixed_sparse(la_fixed), &
                      suba(isub_loc)%a_fixed_sparse(la_fixed) )
            if (suba(isub_loc)%is_bc_present) then
               call sm_apply_bc(suba(isub_loc)%matrixtype,&
                                suba(isub_loc)%ifix,suba(isub_loc)%lifix,suba(isub_loc)%fixv,suba(isub_loc)%lfixv,&
                                nnza, suba(isub_loc)%i_a_sparse,suba(isub_loc)%j_a_sparse,suba(isub_loc)%a_sparse,la, bc,lbc, &
                                store_fixed_rows, suba(isub_loc)%la_fixed,&
                                suba(isub_loc)%i_a_fixed_sparse,&
                                suba(isub_loc)%j_a_fixed_sparse,&
                                suba(isub_loc)%a_fixed_sparse)
               if ( store_fixed_rows ) then
                   !print *, 'la_fixed, nnza_fixed',suba(isub_loc)%la_fixed,suba(isub_loc)%nnza_fixed
                   !call sm_print(6, suba(isub_loc)%i_a_fixed_sparse, suba(isub_loc)%j_a_fixed_sparse, suba(isub_loc)%a_fixed_sparse, &
                   !              suba(isub_loc)%la_fixed, suba(isub_loc)%la_fixed)
                   !call flush(6)
                   call sm_assembly(suba(isub_loc)%i_a_fixed_sparse, &
                                    suba(isub_loc)%j_a_fixed_sparse, &
                                    suba(isub_loc)%a_fixed_sparse,& 
                                    suba(isub_loc)%la_fixed, suba(isub_loc)%nnza_fixed)
               end if
            end if
            suba(isub_loc)%is_fixed_rows_matrix_loaded = .true.
         end if

         if (is_bc_present) then
            ndofis = suba(isub_loc)%ndofi
            lbci = ndofis
            allocate(bci(lbci))
            bci(:) = 0._kr
            call dd_map_sub_to_subi(suba(isub_loc), bc,lbc, bci,lbci) 
            call dd_comm_upload(suba(isub_loc), bci,lbci) 
            call dd_load_eliminated_bc(suba(isub_loc), bc,lbc)
            deallocate(bci)
            deallocate(bc)
         end if
         if (debug) then
            call info(routine_name,' matrix loaded for subdomain',isub)
         end if
      end do
      if (is_bc_present) then
         call dd_comm_swapdata(suba,lsuba, indexsub,lindexsub, sub2proc,lsub2proc,comm_all)
         ! finalize input of boundary conditions
         do isub_loc = 1,lindexsub
            isub = indexsub(isub_loc)

            ndofis = suba(isub_loc)%ndofi
            lbci = ndofis
            allocate(bci(lbci))
            bci(:) = 0._kr

            call dd_comm_download(suba(isub_loc), bci,lbci) 
            ! add correction from neighbours
            call dd_map_subi_to_sub(suba(isub_loc), bci,lbci, suba(isub_loc)%bc,suba(isub_loc)%lbc) 

            deallocate(bci)

            suba(isub_loc)%is_bc_present = .true.

            if (debug) then
               call info(routine_name,' bc loaded for subdomain',isub)
            end if
         end do
      end if


end subroutine

!*******************************************************************************************
subroutine dd_compute_reactions(suba,lsuba, comm_all, sub2proc,lsub2proc,indexsub,lindexsub)
!*******************************************************************************************
! Subroutine for eliminating constraints from RHS and fixing them in the matrix
      use module_sm
      use module_utils
      implicit none
      include "mpif.h"

! array of sub structure
      integer,intent(in) ::                lsuba
      type(subdomain_type),intent(inout) :: suba(lsuba)
! communicator
      integer,intent(in) :: comm_all
! division of subdomains to processors
      integer,intent(in) :: lsub2proc
      integer,intent(in) ::  sub2proc(lsub2proc)
! global indices of local subdomains
      integer,intent(in) :: lindexsub
      integer,intent(in) ::  indexsub(lindexsub)

! local variables
      character(*),parameter:: routine_name = 'DD_COMPUTE_REACTIONS'
      integer :: ierr
      logical :: is_bc_present_loc, is_bc_present

!     BC arrays
      integer::              lreai
      real(kr),allocatable :: reai(:)

      integer :: isub, isub_loc, ndofs, ndofis
      integer :: lrea

      ! find if any nonzero BC is present
      is_bc_present_loc = .false.
      do isub_loc = 1,lindexsub
         if (suba(isub_loc)%is_bc_present) then
            is_bc_present_loc = .true.
         end if
      end do
!*****************************************************************MPI
      call MPI_ALLREDUCE(is_bc_present_loc,is_bc_present,1, MPI_LOGICAL, MPI_LOR, comm_all, ierr) 
!*****************************************************************MPI

      ! eliminate constraints if they are present
      do isub_loc = 1,lindexsub
         isub = indexsub(isub_loc)

         if (suba(isub_loc)%is_bc_present .and. .not. suba(isub_loc)%is_fixed_rows_matrix_loaded) then
            print *, suba(isub_loc)%is_bc_present, suba(isub_loc)%is_fixed_rows_matrix_loaded
            call error( routine_name, &
                        'Matrix of fixed rows not prepared,'//&
                        'not asked for reactions in DD_FIX_CONSTRAINTS for subdomain',&
                        suba(isub_loc)%isub )
         end if

         ! eliminate boundary conditions
         ndofs =  suba(isub_loc)%ndof

         lrea = ndofs
         suba(isub_loc)%lrea = lrea
         if (.not.allocated(suba(isub_loc)%rea)) then
            allocate(suba(isub_loc)%rea(lrea))
         end if
         suba(isub_loc)%rea(:) = 0._kr

         ! eliminate natural BC
         if (suba(isub_loc)%is_bc_present) then
            call sm_vec_mult(suba(isub_loc)%matrixtype, suba(isub_loc)%nnza_fixed, &
                             suba(isub_loc)%i_a_fixed_sparse, suba(isub_loc)%j_a_fixed_sparse, suba(isub_loc)%a_fixed_sparse, &
                             suba(isub_loc)%la_fixed, &
                             suba(isub_loc)%sol,suba(isub_loc)%lsol, suba(isub_loc)%rea,suba(isub_loc)%lrea )
            ! at the moment, rea = A_fixed*sol
            ! change it to rea = A_fixed * sol - rhs
            if (.not. suba(isub_loc)%is_rhs_complete) then
               if (suba(isub_loc)%is_bc_present) then
                  where (suba(isub_loc)%ifix .ne. 0 ) suba(isub_loc)%rea = suba(isub_loc)%rea - suba(isub_loc)%rhs
               end if
            end if

            ndofis = suba(isub_loc)%ndofi
            lreai = ndofis
            allocate(reai(lreai))
            reai(:) = 0._kr
            call dd_map_sub_to_subi(suba(isub_loc), suba(isub_loc)%rea,suba(isub_loc)%lrea, reai,lreai) 
            call dd_comm_upload(suba(isub_loc), reai,lreai) 
            deallocate(reai)
         end if
      end do
      if (is_bc_present) then
         call dd_comm_swapdata(suba,lsuba, indexsub,lindexsub, sub2proc,lsub2proc,comm_all)
         ! finalize input of boundary conditions
         do isub_loc = 1,lindexsub
            isub = indexsub(isub_loc)

            ndofis = suba(isub_loc)%ndofi
            lreai = ndofis
            allocate(reai(lreai))
            reai(:) = 0._kr

            call dd_comm_download(suba(isub_loc), reai,lreai) 
            ! add correction from neighbours
            call dd_map_subi_to_sub(suba(isub_loc), reai,lreai, suba(isub_loc)%rea,suba(isub_loc)%lrea) 

            deallocate(reai)

            ! at the moment, rea = A_fixed*sol
            ! change it to rea = A_fixed * sol - rhs
            if (suba(isub_loc)%is_rhs_complete) then
               if (suba(isub_loc)%is_bc_present .and. suba(isub_loc)%lifix.gt.0 ) then
                  if ( suba(isub_loc)%lifix .ne. suba(isub_loc)%lrea .or. suba(isub_loc)%lifix .ne. suba(isub_loc)%lrhs ) then
                     print *, suba(isub_loc)%lifix, suba(isub_loc)%lrea, suba(isub_loc)%lrhs
                     call error( routine_name, 'Dimension of arrays mismatch.')
                  end if
                  where (suba(isub_loc)%ifix .ne. 0 ) suba(isub_loc)%rea = suba(isub_loc)%rea - suba(isub_loc)%rhs
               end if
            end if

            suba(isub_loc)%is_reactions_ready = .true.
         end do
      end if


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
      character(*),parameter:: routine_name = 'DD_GATHER_MATRIX'

      ! MPI related arrays and variables
      integer :: myid, nproc, ierr, ireq, nreq, nreqx
      integer :: nsendsx, nrecvsx
      integer ::            lrequest
      integer,allocatable :: request(:)
      integer             :: lstatarray1
      integer             :: lstatarray2
      integer,allocatable :: statarray(:,:)

      integer :: nelems, nnods, ndofs
      !integer :: nevax_sub, nevax_loc, nevax, lelmx, lelm
      integer :: lelm
      !real(kr),allocatable :: elm(:)

!     Matrix in IJA sparse format - triplet
      integer::  nnza, la

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

      integer :: ie, iproc, inods, isub, isub_loc, ie_loc, indel_loc,&
                 iprocelm, iprocsub
      integer :: nne, ndofe
      integer,allocatable  :: nsubelem(:)
      integer,allocatable  :: pinet(:)


! check dimension
      if (nelem.ne.liets) then
         call error(routine_name,'wrong length of array IETS')
      end if

! orient in communicators
      call MPI_COMM_RANK(comm_all,myid,ierr)
      call MPI_COMM_SIZE(comm_all,nproc,ierr)

      if (nproc+1.ne.lsub2proc) then
         call error(routine_name,'wrong length of array sub2proc')
      end if
      if (nproc+1.ne.lelm2proc) then
         call error(routine_name,'wrong length of array elm2proc')
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
      !nevax_loc = 0
      nrecvsx = 0
      do isub_loc = 1,lindexsub
         isub = indexsub(isub_loc)

         if (subp(isub).ne.myid) then
            cycle
         end if

         ! get data from subdomain
         nelems = suba(isub_loc)%nelem
         nnods  = suba(isub_loc)%nnod

         nrecvsx = nrecvsx + nelems

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
            call error(routine_name,'strange type of matrix',matrixtype)
         end if
         suba(isub_loc)%la       = la
         allocate(suba(isub_loc)%i_a_sparse(la))
         allocate(suba(isub_loc)%j_a_sparse(la))

         ! Creation of field KDOF(NNOD) with addresses before first global dof of node
         lkdof = nnods
         allocate(kdof(lkdof))
         if (lkdof.gt.0) then
            kdof(1) = 0
         end if
         do inods = 2,nnods
            kdof(inods) = kdof(inods-1) + suba(isub_loc)%nndf(inods-1)
         end do
         ! Prepare numbering of element matrices
         call sm_pmd_make_element_numbering(matrixtype,nelems,suba(isub_loc)%inet,suba(isub_loc)%linet,&
                                            suba(isub_loc)%nnet,suba(isub_loc)%lnnet,&
                                            suba(isub_loc)%nndf,suba(isub_loc)%lnndf,&
                                            kdof,lkdof,&
                                            suba(isub_loc)%i_a_sparse, suba(isub_loc)%j_a_sparse, la)
         deallocate(kdof)

         allocate(suba(isub_loc)%a_sparse(la))
         suba(isub_loc)%nnza = 0

         ! find nevax_sub
         !call pp_get_nevax(nelems,suba(isub_loc)%inet,suba(isub_loc)%linet,&
         !                  suba(isub_loc)%nnet,suba(isub_loc)%lnnet,&
         !                  suba(isub_loc)%nndf,suba(isub_loc)%lnndf,nevax_sub)
         !if (nevax_sub.gt.nevax_loc) then
         !   nevax_loc = nevax_sub
         !end if
      end do
      ! get global maximum of nevax
!*****************************************************************MPI
      !call MPI_ALLREDUCE(nevax_loc,nevax,1, MPI_INTEGER, MPI_MAX, comm_all, ierr) 
!*****************************************************************MPI

      ! prepare memory for one element matrix
      ! determine length by nevax
      !if (matrixtype.eq.1 .or. matrixtype.eq.2) then
      !   lelmx = (nevax+1)*nevax / 2
      !else
      !   lelmx = nevax*nevax
      !end if
      !allocate(elm(lelmx))

      allocate(pinet(lindexsub))
      pinet = 0
      allocate(nsubelem(lindexsub))
      nsubelem = 0

      nsendsx = elm2proc(myid+2) - elm2proc(myid+1) ! not excluding copies in memory
      nreqx = nsendsx + nrecvsx

      ! prepare MPI data for processor
      ! double the size for two sided non-blocking communication
      lrequest = nreqx
      allocate(request(lrequest))
      lstatarray1 = MPI_STATUS_SIZE
      lstatarray2 = nreqx
      allocate(statarray(lstatarray1,lstatarray2))

      ! counter of requests
      ireq = 0

      ! loop over elements
      do ie = 1,nelem
         
         isub = iets(ie)

         ! get processor taking care of this subdomain
         iprocsub = subp(isub)
         ! get processor taking care of this element
         iprocelm = elmp(ie)

         if (myid.eq.iprocelm) then
            indel_loc = locelmnumber(ie) 
            if (myid.ne.iprocsub) then
            ! send messages
               !call MPI_SEND(elma(indel_loc)%lcoarsem,1,                      &
               !              MPI_INTEGER,     iprocsub,ie,comm_all,ierr)
               ireq = ireq + 1
               call MPI_ISEND(elma(indel_loc)%coarsem,elma(indel_loc)%lcoarsem,&
                              MPI_DOUBLE_PRECISION,iprocsub,ie,comm_all,request(ireq),ierr)
            !else
               ! copy it to array elm
               !lelm = elma(indel_loc)%lcoarsem
               !do i = 1,lelm
               !   elm(i) = elma(indel_loc)%coarsem(i)
               !end do
            end if
         end if
         if (myid.eq.iprocsub) then

            isub_loc = locsubnumber(isub)
            ! increase counter of local elements in subdomain
            nsubelem(isub_loc) = nsubelem(isub_loc) + 1
            ie_loc = nsubelem(isub_loc) ! index of element in local numbering

            nne = suba(isub_loc)%nnet(ie_loc)
            ndofe = sum(suba(isub_loc)%nndf(suba(isub_loc)%inet(pinet(isub_loc)+1:pinet(isub_loc)+nne)))
            if (matrixtype.eq.1 .or. matrixtype.eq.2) then
               lelm = (ndofe+1)*ndofe / 2
            else
               lelm = ndofe*ndofe
            end if

            ! storing the matrix - get it from owner by MPI
            if (suba(isub_loc)%la > 0) then
               if (myid.ne.iprocelm) then
                  !call MPI_RECV(lelm,1,  MPI_INTEGER,         iprocelm,ie,comm_all,stat,ierr)
                  ireq = ireq + 1
                  call MPI_IRECV(suba(isub_loc)%a_sparse(suba(isub_loc)%nnza + 1),lelm,&
                                 MPI_DOUBLE_PRECISION,iprocelm,ie,comm_all,request(ireq),ierr)
               else
                  indel_loc = locelmnumber(ie) 
                  suba(isub_loc)%a_sparse(suba(isub_loc)%nnza + 1:suba(isub_loc)%nnza + lelm) = elma(indel_loc)%coarsem
               end if
            end if
            suba(isub_loc)%nnza = suba(isub_loc)%nnza + lelm

            !if (lelm /= lelm_expected) then
            !   call error(routine_name, "Mismatch in the expected element length.")
            !end if
            !call info(routine_name, "I am successfully here!")

            !do i = 1,lelm
            !   suba(isub_loc)%a_sparse(suba(isub_loc)%nnza + i) = elm(i)
            !end do

            pinet(isub_loc) = pinet(isub_loc) + nne

         end if
      end do

      nreq = ireq
      if (nreq.gt.0) then
         call MPI_WAITALL(nreq, request, statarray, ierr)
      end if

      !if (nsends /= ireqs) then
      !   print *, "nsends, ireq", nsends, ireqs
      !   call error(routine_name, "Number of send does not match.")
      !end if
      
      ! free memory
      !deallocate(elm)

      ! finalize input of sparse matrix at each subdomain
      do isub_loc = 1,lindexsub
         isub = indexsub(isub_loc)

         if (subp(isub).ne.myid) then
            cycle
         end if

         la   = suba(isub_loc)%la
         nnza = suba(isub_loc)%nnza

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
            call sm_apply_bc(matrixtype,suba(isub_loc)%ifix,suba(isub_loc)%lifix,suba(isub_loc)%fixv,suba(isub_loc)%lfixv,&
                             nnza,suba(isub_loc)%i_a_sparse,suba(isub_loc)%j_a_sparse,suba(isub_loc)%a_sparse,la, bc,lbc, &
                             .false., 0 )
            if (suba(isub_loc)%is_bc_nonzero) then
               call dd_load_eliminated_bc(suba(isub_loc), bc,lbc)
            end if
         end if

         if (allocated(bc)) then
            deallocate(bc)
         end if
         if (debug) then
            call info(routine_name,'matrix loaded for subdomain',isub)
         end if
      end do

      deallocate(locsubnumber)
      deallocate(subp)
      deallocate(locelmnumber)
      deallocate(elmp)

      deallocate(pinet)
      deallocate(nsubelem)

end subroutine

!!*************************************************
!subroutine dd_set_mirror_subdomain(sub_in,sub_out)
!!*************************************************
!! Subroutine for mirroring data of sub_in in sub_out
!! mirroring = integers are copied, arrays are set as pointers
!      use module_utils
!      implicit none
!
!      type(subdomain_type),intent(in)    :: sub_in
!      type(subdomain_type),intent(inout) :: sub_out
!
!      ! local vars
!      character(*),parameter:: routine_name = 'DD_SET_MIRROR_SUBDOMAIN'
!
!      ! check the prerequisities
!      if (.not.sub_in%is_initialized) then
!         call error(routine_name,'Subdomain not initialized:', sub_in%isub)
!      end if
!      if (.not.sub_in%is_interface_loaded) then
!         call error(routine_name,'Interface is not loaded for subdomain:', sub_in%isub)
!      end if
!      if (.not.sub_in%is_matrix_loaded) then
!         call error(routine_name,'Matrix is not loaded for subdomain:', sub_in%isub)
!      end if
!
!      ! if checks are OK, copy data and set pointers
!      sub_out%is_degenerated = sub_in%is_degenerated 
!
!      sub_out%nelem = sub_in%nelem   ! number of elements
!      sub_out%nnod  = sub_in%nnod    ! number of nodes
!      sub_out%ndof  = sub_in%ndof    ! number of degrees of freedom
!      sub_out%ndim  = sub_in%ndim    ! dimension of the problem
!      ! description of subdomain mesh
!      sub_out%linet   =  sub_in%linet   ! length of INET array 
!      sub_out%inet    => sub_in%inet    ! INET array - indices of nodes on elements
!      sub_out%lnnet   =  sub_in%lnnet   ! length of NNET array
!      sub_out%nnet    => sub_in%nnet    ! NNET array - number of nodes on elements
!      sub_out%lnndf   =  sub_in%lnndf   ! length of NNDF array
!      sub_out%nndf    => sub_in%nndf    ! NNDF array - number of nodal degrees of freedom
!      sub_out%lisngn  =  sub_in%lisngn  ! length of array ISNGN
!      sub_out%isngn   => sub_in%isngn   ! ISNGN array - indices of subdomain nodes in global numbering
!      sub_out%lisvgvn =  sub_in%lisvgvn ! length of array ISVGVN
!      sub_out%isvgvn  => sub_in%isvgvn  ! ISVGVN array - indices of subdomain variables in global variable numbering
!      sub_out%lisegn  =  sub_in%lisegn  ! length of array ISEGN
!      sub_out%isegn   => sub_in%isegn   ! ISEGN array - indices of subdomain elements in global numbering
!      sub_out%lxyz1   =  sub_in%lxyz1   ! length of array of coordinates
!      sub_out%lxyz2   =  sub_in%lxyz2   ! number of x,y,z vectors
!      sub_out%xyz     => sub_in%xyz     ! array of coordinates
!      sub_out%is_mesh_loaded = sub_in%is_mesh_loaded
!      
!      ! description of subdomain interface
!      sub_out%nnodi   =  sub_in%nnodi     ! number of nodes on interface
!      sub_out%ndofi   =  sub_in%ndofi     ! number of dof on Iterface
!      sub_out%ndofo   =  sub_in%ndofo     ! number of dof in iteriOr (ndof = ndofi + ndofo)
!      sub_out%liin    =  sub_in%liin      ! length of IIN array 
!      sub_out%iin     => sub_in%iin       ! IIN array - indices of interface nodes
!      sub_out%liivsvn =  sub_in%liivsvn   ! length of IIVSVN array 
!      sub_out%iivsvn  => sub_in%iivsvn    ! IIVSVN array - indices of Interface variables in subdomain variable numbering
!      sub_out%liovsvn =  sub_in%liovsvn   ! length of IOVSVN array
!      sub_out%iovsvn  => sub_in%iovsvn    ! IOVSVN array - indices of interiOr variables in subdomain variable numbering
!      sub_out%is_interface_loaded = sub_in%is_interface_loaded
!
!      ! boundary conditions
!      sub_out%lifix         =  sub_in%lifix         ! length of IFIX array
!      sub_out%ifix          => sub_in%ifix          ! IFIX array - indices of fixed variables
!      sub_out%lfixv         =  sub_in%lfixv         ! length of FIXV array
!      sub_out%fixv          => sub_in%fixv          ! FIXV array - fixed variables values
!      sub_out%lbc           =  sub_in%lbc           ! length of BC array
!      sub_out%bc            => sub_in%bc            ! BC array - eliminated entries of stiffness matrix multiplied by values of fixed variables 
!      sub_out%is_bc_present = sub_in%is_bc_present  ! are some Dirichlet BC on subdomain?
!      sub_out%is_bc_nonzero = sub_in%is_bc_nonzero  ! are some Dirichlet BC nonzero?
!      sub_out%is_bc_loaded  = sub_in%is_bc_loaded  
!
!      ! description of subdomain matrix
!      sub_out%matrixtype   =  sub_in%matrixtype  
!      sub_out%istorage     =  sub_in%istorage    
!      sub_out%is_assembled =  sub_in%is_assembled
!      sub_out%is_triplet   =  sub_in%is_triplet 
!      sub_out%nnza         =  sub_in%nnza       
!      sub_out%la           =  sub_in%la         
!      sub_out%i_a_sparse   => sub_in%i_a_sparse
!      sub_out%j_a_sparse   => sub_in%j_a_sparse
!      sub_out%a_sparse     => sub_in%a_sparse
!      sub_out%is_matrix_loaded = sub_in%is_matrix_loaded
!
!end subroutine


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
subroutine dd_localize_mesh(sub,isub,ndim,meshdim,nelem,nnod,&
                            inet,linet,nnet,lnnet,nndf,lnndf,xyz,lxyz1,lxyz2,&
                            iets,liets,&
                            find_components, use_dual_mesh_graph, neighbouring)
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
! mesh dimension
      integer,intent(in) :: meshdim
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

! detection of components
      logical,intent(in) :: find_components      
      logical,intent(in) :: use_dual_mesh_graph  
      integer,intent(in) :: neighbouring ! element_neighbouring

! local vars
      character(*),parameter:: routine_name = 'DD_LOCALIZE_MESH'
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
         call error(routine_name,'subdomain index mismatch for subdomain:', isub)
      end if
! check prerequisites
      if (.not.sub%is_initialized) then
         call error(routine_name,'subdomain not initialized.',isub)
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

      call dd_upload_sub_mesh(sub, nelems, nnods, ndofs, ndim, meshdim, &
                              nndfs,lnndfs, nnets,lnnets, 0, inets,linets, isngns,lisngns, &
                              isvgvns,lisvgvns, isegns,lisegns,&
                              xyzs,lxyzs1,lxyzs2, &
                              find_components, use_dual_mesh_graph, neighbouring) 

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
         call warning('DD_LOCALIZE_ADJ','zero kadjsub length.')
      end if
      if (nsub.ne.sub%nsub) then
         call error('DD_LOCALIZE_ADJ','number of subdomains mismatch.')
      end if

      isub = sub%isub

      nadjs = 0
      if (lkadjsub.gt.0) then
         do jsub = 1,nsub
            if (kadjsub(jsub).eq.1) then
               nadjs = nadjs + 1
            end if
         end do
      end if
      ! debug
      !print *,'nadjs for sub',isub,':',nadjs

      ! create list of neighbouring subdomains
      liadjs = nadjs
      allocate(iadjs(liadjs))
      indadjs = 0
      if (lkadjsub.gt.0) then
         do jsub = 1,nsub
            if (kadjsub(jsub).eq.1) then
               indadjs = indadjs + 1
               iadjs(indadjs) = jsub
            end if
         end do
      end if
      ! debug
      !print *,'iadjs for sub',isub,':',iadjs
      call flush(6)
      call dd_upload_sub_adj(sub, nadjs, iadjs,liadjs)
      deallocate(iadjs)

end subroutine

!!*********************************************************************************
!subroutine dd_localize_cornersglobs(sub,ncorner,inodc,linodc,&
!                                    nedge,nface,nnglb,lnnglb,inglb,linglb, nnodcs)
!!*********************************************************************************
!! Subroutine for localization of corners and globs to particular subdomain
!! loads the data directly to the structure
!      use module_utils
!      implicit none
!
!! Subdomain structure
!      type(subdomain_type),intent(inout) :: sub
!! number of corners
!      integer,intent(in) :: ncorner
!! global indices of corners
!      integer,intent(in) :: linodc
!      integer,intent(in) ::  inodc(linodc)
!! number of edges
!      integer,intent(in) :: nedge
!! number of faces
!      integer,intent(in) :: nface
!! number of nodes in globs
!      integer,intent(in) :: lnnglb
!      integer,intent(in) ::  nnglb(lnnglb)
!! global indices of nodes in globs (lenght sum(nnglb))
!      integer,intent(in) :: linglb
!      integer,intent(in) ::  inglb(linglb)
!
!! number of coarse pseudonodes (corners and globs) at subdomain
!      integer,intent(out) ::  nnodcs
!
!      ! local vars
!      integer :: nglb
!
!      integer ::            lglobal_corner_numbers,   licnsins
!      integer,allocatable :: global_corner_numbers(:), icnsins(:)
!      integer::            lglobal_glob_numbers
!      integer,allocatable:: global_glob_numbers(:)
!      integer::            lglob_types
!      integer,allocatable:: glob_types(:)
!      integer::            lnglobnodess
!      integer,allocatable:: nglobnodess(:)
!      integer::            lnglobvars
!      integer,allocatable:: nglobvars(:)
!      integer::            lignsins1, lignsins2
!      integer,allocatable:: ignsins(:,:)
!      integer::            ligvsivns1, ligvsivns2
!      integer,allocatable:: igvsivns(:,:)
!      integer ::            liingns
!      integer,allocatable :: iingns(:)
!      integer ::            lkdofs
!      integer,allocatable :: kdofs(:)
!
!      integer :: idofn, inods, iglb, iglbn, iglbv, iglobs, inc, indnc, &
!                 indng, indns, indvs, inodcs, pointinglb
!      integer :: indins, indivs, indn1
!      integer :: ncorners, nglobs
!      integer :: nnods, ndofs, nelems, nnodis, ndofis, ndofn, nglbn, nglbv
!
!      ! check dimension
!      if (ncorner.ne.linodc) then
!         call error('DD_LOCALIZE_CORNERSGLOBS','array dimension mismatch INODC.')
!      end if
!      nglb = nedge + nface
!      if (nglb.ne.lnnglb) then
!         call error('DD_LOCALIZE_CORNERSGLOBS','array dimension mismatch NNGLB.')
!      end if
!      if (.not.sub%is_interface_loaded) then
!         call error('DD_LOCALIZE_CORNERSGLOBS','Interface not loaded yet.')
!      end if
!
!      ! prepare array for global indices of interface nodes 
!      call dd_get_interface_size(sub,ndofis,nnodis)
!      liingns = nnodis
!      allocate(iingns(liingns))
!      call dd_get_interface_global_numbers(sub, iingns,liingns)
!
!! find number of coarse nodes on subdomain NCORNERS
!      ncorners = 0
!      do inc = 1,ncorner
!         indnc = inodc(inc)
!         if (any(iingns.eq.indnc)) then
!            ncorners = ncorners + 1
!         end if
!      end do
!
!! create array kdofs
!      call dd_get_size(sub, ndofs,nnods,nelems)
!      lkdofs = nnods
!      allocate(kdofs(lkdofs))
!      if (nnods.gt.0) then
!         kdofs(1) = 0
!         do inods = 2,nnods
!            kdofs(inods) = kdofs(inods-1) + sub%nndf(inods-1)
!         end do
!      end if
!      !print *,'kdofs:',kdofs
!      !call flush(6)
!
!
!      ! find mapping of corners
!      lglobal_corner_numbers = ncorners
!      allocate(global_corner_numbers(lglobal_corner_numbers))
!      licnsins = ncorners
!      allocate(icnsins(licnsins))
!
!      inodcs = 0
!      do inc = 1,ncorner
!         indnc = inodc(inc)
!         if (any(iingns.eq.indnc)) then
!            inodcs = inodcs + 1
!
!            ! mapping to global corner numbers
!            global_corner_numbers(inodcs) = inc
!            ! mapping to subdomain interface numbers
!            call get_index(indnc,iingns,nnodis,indins)
!            if (indins .eq. -1) then
!               write(*,*) 'Index of subdomain interface node not found.', indnc
!               write(*,*) 'iingns',iingns
!               call error_exit
!            end if
!            icnsins(inodcs) = indins
!         end if
!      end do
!
!      call dd_upload_sub_corners(sub, ncorners, global_corner_numbers,lglobal_corner_numbers, icnsins,licnsins)
!      deallocate(icnsins)
!      deallocate(global_corner_numbers)
!
!      ! find local number of globs NGLOBS
!      nglobs     = 0
!      pointinglb = 0
!      do iglb = 1,nglb
!         nglbn = nnglb(iglb)
!
!         ! touch first node in glob
!         indn1 = inglb(pointinglb + 1)
!
!         if (any(iingns.eq.indn1)) then
!            nglobs = nglobs + 1
!         end if
!
!         pointinglb = pointinglb + nglbn
!      end do
!
!      ! get array of interface 
!
!      ! mapping of globs
!      lglobal_glob_numbers = nglobs
!      allocate(global_glob_numbers(lglobal_glob_numbers))
!      lnglobvars = nglobs
!      allocate(nglobvars(lnglobvars))
!      lnglobnodess = nglobs
!      allocate(nglobnodess(lnglobnodess))
!      lglob_types = nglobs
!      allocate(glob_types(lglob_types))
!
!      iglobs     = 0
!      pointinglb = 0
!      do iglb = 1,nglb
!         nglbn = nnglb(iglb)
!
!         ! touch first node in glob
!         indn1 = inglb(pointinglb + 1)
!
!         if (any(iingns.eq.indn1)) then
!
!            iglobs = iglobs + 1
!
!            nglbv = 0
!            do iglbn = 1,nglbn
!               indng = inglb(pointinglb + iglbn)
!               call get_index(indng,iingns,liingns,indins)
!               if (indins .eq. -1) then
!                  write(*,*) ' Index of interface node not found for global ', indng
!                  call error_exit
!               end if
!               indns = sub%iin(indins)
!
!               ndofn = sub%nndf(indns)
!
!               nglbv = nglbv + ndofn
!            end do
!
!            nglobvars(iglobs)   = nglbv
!            nglobnodess(iglobs) = nglbn
!            global_glob_numbers(iglobs) = iglb
!         end if
!
!         pointinglb = pointinglb + nglbn
!      end do
!
!      ! set type of glob
!      glob_types = 1
!      where (global_glob_numbers .le. nedge) glob_types = 2
!      ! shift numbering behind corners
!      global_glob_numbers = global_glob_numbers + ncorner
!
!      ligvsivns1 = nglobs
!      ligvsivns2 = maxval(nglobvars)
!      allocate(igvsivns(ligvsivns1,ligvsivns2))
!      lignsins1 = nglobs
!      lignsins2 = maxval(nglobnodess)
!      allocate(ignsins(lignsins1,lignsins2))
!      iglobs     = 0
!      pointinglb = 0
!      do iglb = 1,nglb
!         nglbn = nnglb(iglb)
!
!         ! touch first node in glob
!         indn1 = inglb(pointinglb + 1)
!
!         if (any(iingns.eq.indn1)) then
!
!            iglobs = iglobs + 1
!
!            iglbv = 0
!            do iglbn = 1,nglbn
!               
!               indng = inglb(pointinglb + iglbn)
!               call get_index(indng,iingns,nnodis,indins)
!               if (indins .eq. -1) then
!                  write(*,*) ' Index of interface node not found.', indng
!                  call error_exit
!               end if
!
!               ignsins(iglobs,iglbn) = indins
!
!               indns = sub%iin(indins)
!               ndofn = sub%nndf(indns)
!
!               do idofn = 1,ndofn
!                  iglbv = iglbv + 1
!
!                  indvs = kdofs(indns) + idofn
!                  call get_index(indvs,sub%iivsvn,sub%liivsvn,indivs)
!                  if (indivs .eq. -1) then
!                     write(*,*) 'DD_LOCALIZE_CORNERSGLOBS: Index of subdomain interface dof not found.'
!                     write(*,*) 'indng =',indng,'indns =',indns,'indvs = ',indvs,'indivs = ',indivs, 'isub = ',sub%isub
!                     call error_exit
!                  end if
!
!                  igvsivns(iglobs,iglbv) = indivs
!               end do
!            end do
!         end if
!
!         pointinglb = pointinglb + nglbn
!      end do
!      deallocate(iingns)
!      deallocate(kdofs)
!
!      call dd_upload_sub_globs(sub, nglobs, global_glob_numbers,lglobal_glob_numbers,&
!                               nglobnodess,lnglobnodess, nglobvars,lnglobvars,&
!                               ignsins,lignsins1,lignsins2, igvsivns,ligvsivns1,ligvsivns2,&
!                               glob_types,lglob_types)
!      deallocate(ignsins)
!      deallocate(igvsivns)
!      deallocate(glob_types)
!      deallocate(nglobnodess)
!      deallocate(nglobvars)
!      deallocate(global_glob_numbers)
!
!      nnodcs = ncorners + nglobs
!
!end subroutine

!*********************************************************************************
subroutine dd_load_matrix_triplet(sub, matrixtype, numshift,&
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

! shift of entries ( 1 for C arrays, 0 for Fortran )
      integer,intent(in) :: numshift
      ! Matrix in IJA sparse format
      integer,intent(in)  :: la
      integer,intent(in)  :: nnza
      integer,intent(in)  :: i_sparse(la), j_sparse(la)
      real(kr),intent(in) :: a_sparse(la)
      logical,intent(in)  :: is_assembled

      ! local vars
      character(*),parameter:: routine_name = 'DD_LOAD_MATRIX_TRIPLET'
      integer :: i
      integer :: la_corr

      ! check that matrix is not loaded
      if (sub%is_matrix_loaded) then
         call error(routine_name,' Matrix already loaded for subdomain: ',sub%isub )
      end if

      ! load data
      sub%matrixtype = matrixtype
      if (matrixtype.eq.0) then
         sub%istorage   = 1
      else if (matrixtype.eq.1 .or. matrixtype.eq.2) then
         sub%istorage   = 2
      else
         call error(routine_name,'strange type of matrix:', matrixtype)
      end if
      if (matrixtype.eq.0 .or. matrixtype.eq.2) then
         ! add some space to avoid problems with e.g. Navier-Stokes equations with the zero block (2,2)
         la_corr = la + 1
      else
         la_corr = la
      end if
      sub%nnza     = nnza
      sub%la       = la_corr
      sub%is_assembled = is_assembled
      allocate(sub%i_a_sparse(la_corr))
      do i = 1,la
         sub%i_a_sparse(i) = i_sparse(i) + numshift
      end do
      do i = la+1,la_corr
         sub%i_a_sparse(i) = 1
      end do
      allocate(sub%j_a_sparse(la_corr))
      do i = 1,la
         sub%j_a_sparse(i) = j_sparse(i) + numshift
      end do
      do i = la+1,la_corr
         sub%j_a_sparse(i) = 1
      end do
      allocate(sub%a_sparse(la_corr))
      do i = 1,la
         sub%a_sparse(i)   = a_sparse(i)
      end do
      do i = la+1,la_corr
         sub%a_sparse(i) = 0._kr
      end do

      sub%is_matrix_loaded = .true.
      sub%is_triplet       = .true.

end subroutine

!*********************************************************************************
subroutine dd_check_local_matrix(sub)
!*********************************************************************************
! Subroutine for checking a sparse triplet, the array is assumed sorted
      use module_sm
      implicit none

! Subdomain structure
      type(subdomain_type),intent(inout) :: sub

      ! local vars
      character(*),parameter:: routine_name = 'DD_CHECK_LOCAL_MATRIX'

      if (.not. sub%is_degenerated) then
         call sm_check_matrix(sub%matrixtype,sub%ndof,sub%ndof,sub%i_a_sparse,sub%j_a_sparse,sub%la,sub%nnza)
      end if

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
      character(*),parameter:: routine_name = 'DD_LOAD_ELIMINATED_BC'
      integer :: i

      ! check if I store the subdomain
      if (.not.sub%is_bc_loaded) then
         call warning(routine_name,'Boundary Conditions not loaded for subdomain ',sub%isub)
      end if

      ! load eliminated boundary conditions if they are present
      if (allocated(sub%bc)) then
         ! just check the size
         if (sub%lbc.ne.lbc) then
            call error(routine_name,'Array size mismatch for subdomain ',sub%isub)
         end if
      else
         sub%lbc = lbc
         allocate(sub%bc(lbc))
      end if
      if (lbc .gt. 0) then
         do i = 1,lbc
            sub%bc(i) = bc(i)
         end do
      end if

end subroutine

!*********************************************************************************
subroutine dd_upload_sub_mesh(sub, nelem, nnod, ndof, ndim, meshdim, &
                              nndf,lnndf, nnet,lnnet, numshift, inet,linet, &
                              isngn,lisngn, isvgvn,lisvgvn, isegn,lisegn,&
                              xyz,lxyz1,lxyz2, &
                              find_components, use_dual_mesh_graph, neighbouring)
!*********************************************************************************
! Subroutine for loading mesh data into sub structure
      use module_utils
      use module_graph
      implicit none

! Subdomain structure
      type(subdomain_type),intent(inout) :: sub
      ! mesh
      integer,intent(in) :: nelem, nnod, ndof, ndim, meshdim
      integer,intent(in) :: numshift
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
      logical,intent(in) :: find_components      
      logical,intent(in) :: use_dual_mesh_graph  
      integer,intent(in) :: neighbouring ! element_neighbouring

      ! local vars
      character(*),parameter:: routine_name = 'DD_UPLOAD_SUB_MESH'
      integer,parameter :: node_neighbouring = 1  ! any element connecting nodes counts as a graph edge (default)
      integer :: i, j
      integer :: ie, indinet, indnode, ine, nne, indcomponent, nelcomponents

      ! check continuity of nodal graph
      integer             :: lnetn,  lietn,  lkietn,   lkinet
      integer, allocatable :: netn(:),ietn(:),kietn(:), kinet(:)
      integer             :: lelcomponents
      integer, allocatable :: elcomponents(:)

      integer ::              graphtype
      integer ::              ngraph_edge
      integer ::              ngraph_vertex
      integer ::              lxadj,   ladjncy,   ladjwgt
      integer,allocatable ::   xadj(:), adjncy(:), adjwgt(:)

      ! check data consistency
      if (.not.sub%is_initialized) then
         call error(routine_name, 'Not initialized subdomain: ',sub%isub)
      end if
      if ( lnnet .ne. nelem ) then
         call error(routine_name, 'LNNET size not equal to NELEM for subdomain: ',sub%isub)
      end if
      if ( lnndf .ne. nnod ) then
         call error(routine_name, 'LNNDF size not equal to NNOD for subdomain: ',sub%isub)
      end if
      if ( lisngn .ne. nnod ) then
         call error(routine_name, 'LISNGN size not equal to NNOD for subdomain: ',sub%isub)
      end if
      if ( lisngn .ne. nnod ) then
         call error(routine_name, 'LISNGN size not equal to NNOD for subdomain: ',sub%isub)
      end if
      if ( lisegn .ne. nelem ) then
         call error(routine_name, 'LISEGN size not equal to NELEM for subdomain: ',sub%isub)
      end if
      if ( lxyz1 .ne. nnod ) then
         call error(routine_name, 'LXYZ1 size not equal to NNOD for subdomain: ',sub%isub)
      end if
      if ( lxyz2 .ne. ndim ) then
         call error(routine_name, 'LXYZ2 size not equal to NDIM for subdomain: ',sub%isub)
      end if
      if (sum(nnet) .ne. linet) then
         call error(routine_name, 'NNET does not sum up to LINET for subdomain: ',sub%isub)
      end if
      if (sum(nndf) .ne. ndof) then
         call error(routine_name, 'NNDF does not sum up to NDOF for subdomain: ',sub%isub)
      end if

      ! load data
      sub%nelem   = nelem
      sub%nnod    = nnod
      sub%ndof    = ndof
      sub%ndim    = ndim
      sub%meshdim = meshdim

      !write (*,*) 'subdomain statistics:'
      !write (*,*) ' nelem:   ', sub%nelem
      !write (*,*) ' nnod:    ', sub%nnod
      !write (*,*) ' ndof:    ', sub%ndof
      !write (*,*) ' ndim:    ', sub%ndim
      !write (*,*) ' meshdim: ', sub%meshdim

      sub%linet   = linet
      allocate(sub%inet(linet))
      do i = 1,linet
         sub%inet(i) = inet(i) + numshift
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
         sub%isngn(i) = isngn(i) + numshift
      end do

      sub%lisvgvn   = lisvgvn
      allocate(sub%isvgvn(lisvgvn))
      do i = 1,lisvgvn
         sub%isvgvn(i) = isvgvn(i) + numshift
      end do

      sub%lisegn   = lisegn
      allocate(sub%isegn(lisegn))
      do i = 1,lisegn
         sub%isegn(i) = isegn(i) + numshift
      end do

      sub%lxyz1   = lxyz1
      sub%lxyz2   = lxyz2
      allocate(sub%xyz(lxyz1,lxyz2))
      do j = 1,lxyz2
         do i = 1,lxyz1
            sub%xyz(i,j) = xyz(i,j)
         end do
      end do

      if (ndof.eq.0) then
         sub%is_degenerated = .true.
      else
         sub%is_degenerated = .false.
      end if
      sub%is_mesh_loaded = .true.

      ! determine nodal components of the mesh
      sub%lnodal_components = nnod
      allocate(sub%nodal_components(sub%lnodal_components))
      if (find_components .and. .not.sub%is_degenerated) then
         ! prepare kinet
         lkinet = nelem
         allocate(kinet(lkinet))
         if (lkinet.gt.0) then
            kinet(1) = 0
         end if
         do i = 2,lkinet
            kinet(i) = kinet(i-1) + nnet(i-1)
         end do

         ! prepare dual mesh
         lnetn  = nnod
         lietn  = linet
         lkietn = nnod
         allocate(netn(lnetn),ietn(lietn),kietn(lkietn))
         call graph_get_dual_mesh(sub%nelem,sub%nnod,&
                                  sub%inet,sub%linet,sub%nnet,sub%lnnet,&
                                  netn,lnetn,ietn,lietn,kietn,lkietn)

         if (use_dual_mesh_graph) then
            ! prepare graph of subdomain elements - i.e. dual graph
            graphtype    = 1 ! weighted
            ngraph_vertex = nelem
            ! how many nodes have to share two elements to call them adjacent in a graph
            ! set properly !!!!!!!!!!!!!!!!!
            !element_neighbouring = meshdim
            !element_neighbouring = 1
            !!!!!!!!!!!!!!!!!!
            call graph_from_mesh(ngraph_vertex,graphtype,neighbouring,&
                                 sub%inet,sub%linet,&
                                 sub%nnet,sub%lnnet,&
                                 ietn,lietn, netn,lnetn,&
                                 kietn,lkietn,&
                                 ngraph_edge, xadj, adjncy, adjwgt)

            lxadj   = size(xadj)
            ladjncy = size(adjncy)
            ladjwgt = size(adjwgt)
            !call graph_check(ngraph_vertex,graphtype, xadj,lxadj, adjncy,ladjncy, adjwgt,ladjwgt)
   
            ! determine continuity of ELEMENT components
            lelcomponents = nelem
            allocate(elcomponents(lelcomponents))
            call graph_components(ngraph_vertex,xadj,lxadj,adjncy,ladjncy,&
                                  elcomponents,lelcomponents, nelcomponents)
            
            ! copy element components into its nodes 
            sub%nnodal_components = nelcomponents
            indinet = 0
            do ie = 1,nelem
               nne = sub%nnet(ie)

               indcomponent = elcomponents(ie)

               do ine = 1,nne
                  indinet = indinet + 1
                  indnode = sub%inet(indinet)
                  sub%nodal_components(indnode) = indcomponent
               end do
            end do

            deallocate(elcomponents)

         else
            ! prepare graph of subdomain nodes - i.e. primal graph
            graphtype    = 0 ! unweighted
            ngraph_vertex = nnod
            call graph_from_mesh(ngraph_vertex,graphtype,node_neighbouring,&
                                 ietn,lietn, netn,lnetn,&
                                 sub%inet,sub%linet,&
                                 sub%nnet,sub%lnnet,&
                                 kinet,lkinet,&
                                 ngraph_edge, xadj, adjncy, adjwgt)
            lxadj   = size(xadj)
            ladjncy = size(adjncy)
            ladjwgt = size(adjwgt)
            call graph_check(ngraph_vertex,graphtype, xadj,lxadj, adjncy,ladjncy, adjwgt,ladjwgt)


            ! determine continuity of components
            call graph_components(ngraph_vertex,xadj,lxadj,adjncy,ladjncy,&
                                  sub%nodal_components,sub%lnodal_components,sub%nnodal_components)
         end if

         ! data for creating components
         sub%find_components = find_components
         sub%use_dual_mesh_graph = use_dual_mesh_graph
         sub%neighbouring = neighbouring

         deallocate(kinet)
         deallocate(netn,ietn,kietn)
         deallocate(xadj)
         deallocate(adjncy,adjwgt)
      else
         ! assume only one component
         sub%nodal_components  = 1
         sub%nnodal_components = 1
      end if
      sub%is_nodal_components_loaded = .true.

      ! debug
      !call info(routine_name,'number of components:',sub%nnodal_components)

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

      ! the routine can be called multiple times to update the information
      if (allocated(sub%iadj)) then
         deallocate(sub%iadj)
      end if

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
                               glob_type,lglob_type, &
                               nsubglobs,lnsubglobs, &
                               glob_subs,lglob_subs1,lglob_subs2)
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
      integer,intent(in) :: lnsubglobs
      integer,intent(in) ::  nsubglobs(lnsubglobs)
      integer,intent(in) :: lglob_subs1, lglob_subs2
      integer,intent(in) ::  glob_subs(lglob_subs1,lglob_subs2)

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

      sub%lnsubglobs = lnsubglobs
      allocate(sub%nsubglobs(lnsubglobs))
      do i = 1,lnsubglobs
         sub%nsubglobs(i) = nsubglobs(i)
      end do

      sub%lglob_subs1 = lglob_subs1
      sub%lglob_subs2 = lglob_subs2
      allocate(sub%glob_subs(lglob_subs1,lglob_subs2))
      do i = 1,lglob_subs1
         do j = 1,lglob_subs2
            sub%glob_subs(i,j) = glob_subs(i,j)
         end do
      end do

      sub%is_globs_loaded = .true.

end subroutine

!***************************************************
subroutine dd_upload_bc(sub, ifix,lifix, fixv,lfixv)
!***************************************************
! Subroutine for loading boundary conditions on subdomain
      use module_utils
      implicit none

! Subdomain structure
      type(subdomain_type),intent(inout) :: sub
      integer,intent(in) :: lifix
      integer,intent(in) ::  ifix(lifix)
      integer,intent(in) :: lfixv
      real(kr),intent(in)::  fixv(lfixv)

      ! local vars
      character(*),parameter:: routine_name = 'DD_UPLOAD_BC'
      integer :: i
      logical :: import_bc

      if (any(ifix.ne.0)) then
         if (.not. sub%is_bc_loaded) then
            sub%is_bc_present = .true.
         end if
         import_bc = .true.
         if (any(fixv.ne.0.0_kr)) then
            sub%is_bc_nonzero = .true.
         else
            sub%is_bc_nonzero = .false.
         end if
      else
         if (.not. sub%is_bc_loaded) then
            sub%is_bc_present = .false.
         end if
         sub%is_bc_nonzero = .false.
         import_bc = .false.
      end if

      ! boundary conditions
      if (import_bc) then
         ! check sizes
         if (sub%is_bc_loaded) then
            if ( sub%lifix .ne. lifix .or.  sub%lfixv .ne. lfixv ) then
               call error( routine_name, 'Size mismatch for subdomain', sub%isub )
            end if
            if (any(sub%ifix .ne. ifix)) then
               call error( routine_name, 'Indices of fixed variables IFIX cannot be changed for subdomain', sub%isub )
            end if
         else
            sub%lifix = lifix
            sub%lfixv = lfixv

            allocate(sub%ifix(sub%lifix))
            allocate(sub%fixv(sub%lfixv))
         end if
      else
         sub%lifix = 0
         sub%lfixv = 0
         if (.not. sub%is_bc_loaded) then
            allocate(sub%ifix(sub%lifix))
            allocate(sub%fixv(sub%lfixv))
         end if
         sub%ifix = 0
         sub%fixv = 0._kr
      end if

      if (import_bc) then
         do i = 1,lifix
            sub%ifix(i) = ifix(i)
         end do
         do i = 1,lfixv
            sub%fixv(i) = fixv(i)
         end do
      end if

      sub%is_bc_loaded = .true.

end subroutine

!*******************************************************
subroutine dd_upload_rhs(sub, rhs,lrhs, is_rhs_complete)
!*******************************************************
! Subroutine for initialization of subdomain right hand side
      use module_utils
      implicit none

! Subdomain structure
      type(subdomain_type),intent(inout) :: sub
      integer,intent(in) :: lrhs
      real(kr),intent(in)::  rhs(lrhs)
      logical :: is_rhs_complete

      ! local vars
      character(*),parameter:: routine_name = 'DD_UPLOAD_RHS'
      integer :: i

      ! is array allocated?
      if (.not.allocated(sub%rhs)) then
         ! if not, allocate it
         sub%lrhs = lrhs
         allocate(sub%rhs(lrhs))
      else
         ! if yes, check that it is allocated to correct dimension
         if (lrhs.ne.sub%lrhs) then
            call error(routine_name,'Dimension of subdomain RHS mismatch for subdomain:',sub%isub)
         end if
      end if
      do i = 1,lrhs
         sub%rhs(i) = rhs(i)
      end do

      sub%is_rhs_loaded = .true.
      sub%is_rhs_complete = is_rhs_complete

      sub%is_reduced_rhs_loaded = .false.

end subroutine

!******************************************************
subroutine dd_upload_interior_solution(sub, solo,lsolo)
!******************************************************
! Subroutine for uploading interior solution to subdomain structure
      use module_utils
      implicit none

! Subdomain structure
      type(subdomain_type),intent(inout) :: sub
      integer,intent(in) :: lsolo
      real(kr),intent(in)::  solo(lsolo)

      ! local vars
      character(*),parameter:: routine_name = 'DD_UPLOAD_INTERIOR_SOLUTION'
      integer :: i
      integer :: ndofo

      ! check prerequisites
      if (.not. sub%is_interface_loaded) then
         call error(routine_name,'Interface not yet loaded.')
      end if

      ndofo = sub%ndofo
      ! check dimension
      if (lsolo.ne.ndofo) then
         call error(routine_name,'Array SOLO dimension mismatch.')
      end if

      ! is array allocated?
      if (.not.allocated(sub%solo)) then
         ! if not, allocate it
         sub%lsolo = lsolo
         allocate(sub%solo(lsolo))
      else
         ! if yes, check that it is allocated to correct dimension
         if (lsolo.ne.sub%lsolo) then
            call error(routine_name,'Dimension of subdomain SOLO mismatch for subdomain:',sub%isub)
         end if
      end if
      do i = 1,lsolo
         sub%solo(i) = solo(i)
      end do

      sub%is_interior_solution_loaded = .true.

end subroutine

!*******************************************
subroutine dd_upload_solution(sub, sol,lsol)
!*******************************************
! Subroutine for uploading solution to subdomain structure
      use module_utils
      implicit none

! Subdomain structure
      type(subdomain_type),intent(inout) :: sub
      integer,intent(in) :: lsol
      real(kr),intent(in)::  sol(lsol)

      ! local vars
      character(*),parameter:: routine_name = 'DD_UPLOAD_SOLUTION'
      integer :: i
      integer :: ndof

      ! check prerequisites
      if (.not. sub%is_initialized) then
         call error( routine_name,'Not initialized subdomain: ',sub%isub)
      end if

      ndof = sub%ndof
      ! check dimension
      if (lsol.ne.ndof) then
         call error(routine_name,'Array SOL dimension mismatch.')
      end if

      ! is array allocated?
      if (.not.allocated(sub%sol)) then
         ! if not, allocate it
         sub%lsol = lsol
         allocate(sub%sol(lsol))
      else
         ! if yes, check that it is allocated to correct dimension
         if (lsol.ne.sub%lsol) then
            call error(routine_name,'Dimension of subdomain SOL mismatch for subdomain:',sub%isub)
         end if
      end if
      do i = 1,lsol
         sub%sol(i) = sol(i)
      end do

      sub%is_solution_loaded = .true.

end subroutine

!*****************************************************************************************************
subroutine dd_upload_sub_user_constraints(sub, user_constraints,luser_constraints1,luser_constraints2)
!*****************************************************************************************************
! Subroutine for uploading solution to subdomain structure
      use module_utils
      implicit none

! Subdomain structure
      type(subdomain_type),intent(inout) :: sub
      integer,intent(in) :: luser_constraints1
      integer,intent(in) :: luser_constraints2
      real(kr),intent(in)::  user_constraints(luser_constraints1*luser_constraints2)

      ! local vars
      character(*),parameter:: routine_name = 'DD_UPLOAD_SUB_USER_CONSTRAINTS'
      integer :: i,j

      ! check prerequisites
      if (.not. sub%is_initialized) then
         call error( routine_name,'Not initialized subdomain: ',sub%isub)
      end if

      ! check dimension
      if (luser_constraints2.ne.sub%nnod .and. luser_constraints2.gt.0) then
         call error(routine_name,'Array USER_CONSTRAINTS dimension mismatch.')
      end if

      ! is array allocated?
      if (.not.allocated(sub%user_constraints)) then
         ! if not, allocate it
         sub%luser_constraints1 = luser_constraints1
         sub%luser_constraints2 = luser_constraints2
         allocate(sub%user_constraints(luser_constraints1,luser_constraints2))
      else
         ! if yes, check that it is allocated to correct dimension
         if (luser_constraints1.ne.sub%luser_constraints1 .or. luser_constraints2.ne.sub%luser_constraints2) then
            call error(routine_name,'Dimension of subdomain USER_CONSTRAINTS mismatch for subdomain:',sub%isub)
         end if
      end if
      ! copy constraints and convert it from linear row-wise stored array to two-dimensional matrix
      do i = 1,luser_constraints1
         do j = 1,luser_constraints2
            sub%user_constraints(i,j) = user_constraints((i-1)*luser_constraints2 + j)
         end do
      end do

      sub%is_user_constraints_loaded = .true.

end subroutine

!*****************************************************************************************************
subroutine dd_upload_sub_element_data(sub, element_data,lelement_data1,lelement_data2)
!*****************************************************************************************************
! Subroutine for uploading solution to subdomain structure
      use module_utils
      implicit none

! Subdomain structure
      type(subdomain_type),intent(inout) :: sub
      integer,intent(in) :: lelement_data1
      integer,intent(in) :: lelement_data2
      real(kr),intent(in)::  element_data(lelement_data1*lelement_data2)

      ! local vars
      character(*),parameter:: routine_name = 'DD_UPLOAD_SUB_ELEMENT_DATA'
      integer :: i,j

      ! check prerequisites
      if (.not. sub%is_initialized) then
         call error( routine_name,'Not initialized subdomain: ',sub%isub)
      end if

      ! check dimension
      if (lelement_data2.ne.sub%nelem .and. lelement_data2.gt.0) then
         call error(routine_name,'Array ELEMENT_DATA dimension mismatch.')
      end if

      ! is array allocated?
      if (.not.allocated(sub%element_data)) then
         ! if not, allocate it
         sub%lelement_data1 = lelement_data1
         sub%lelement_data2 = lelement_data2
         allocate(sub%element_data(lelement_data1,lelement_data2))
      else
         ! if yes, check that it is allocated to correct dimension
         if (lelement_data1.ne.sub%lelement_data1 .or. lelement_data2.ne.sub%lelement_data2) then
            call error(routine_name,'Dimension of subdomain ELEMENT_DATA mismatch for subdomain:',sub%isub)
         end if
      end if
      ! copy constraints and convert it from linear row-wise stored array to two-dimensional matrix
      do i = 1,lelement_data1
         do j = 1,lelement_data2
            sub%element_data(i,j) = element_data((i-1)*lelement_data2 + j)
         end do
      end do

      sub%is_element_data_loaded = .true.

end subroutine

!*****************************************************************************************************
subroutine dd_upload_sub_dof_data(sub, dof_data,ldof_data)
!*****************************************************************************************************
! Subroutine for uploading solution to subdomain structure
      use module_utils
      implicit none

! Subdomain structure
      type(subdomain_type),intent(inout) :: sub
      integer,intent(in) :: ldof_data
      real(kr),intent(in)::  dof_data(ldof_data)

      ! local vars
      character(*),parameter:: routine_name = 'DD_UPLOAD_SUB_DOF_DATA'
      integer :: i

      ! check prerequisites
      if (.not. sub%is_initialized) then
         call error( routine_name,'Not initialized subdomain: ',sub%isub)
      end if

      ! is array allocated?
      if (.not.allocated(sub%dof_data)) then
         ! if not, allocate it
         sub%ldof_data = ldof_data
         allocate(sub%dof_data(ldof_data))
      else
         ! if yes, check that it is allocated to correct dimension
         if (ldof_data.ne.sub%ldof_data ) then
            call error(routine_name,'Dimension of subdomain DOF_DATA mismatch for subdomain:',sub%isub)
         end if
      end if
      ! copy constraints and convert it from linear row-wise stored array to two-dimensional matrix
      do i = 1,ldof_data
         sub%dof_data(i) = dof_data(i)
      end do

      sub%is_dof_data_loaded = .true.

end subroutine

!*********************************************
subroutine dd_download_solution(sub, sol,lsol)
!*********************************************
! Subroutine for downloading solution from subdomain structure
      use module_utils
      implicit none

! Subdomain structure
      type(subdomain_type),intent(inout) :: sub
! Subdomain solution
      integer,intent(in) ::  lsol
      real(kr),intent(out)::  sol(lsol)

      ! local vars
      character(*),parameter:: routine_name = 'DD_DOWNLOAD_SOLUTION'
      integer :: i
      integer :: ndof

      ! check prerequisites
      if (.not. sub%is_solution_loaded) then
         call error( routine_name,'Solution is not loaded for subdomain: ',sub%isub)
      end if

      ndof = sub%ndof
      ! check dimension
      if (lsol.ne.ndof) then
         call error(routine_name,'Array SOL dimension mismatch.')
      end if

      do i = 1,lsol
         sol(i) = sub%sol(i)
      end do

end subroutine

!**********************************************
subroutine dd_download_reactions(sub, rea,lrea)
!**********************************************
! Subroutine for downloading reactions from subdomain structure
      use module_utils
      implicit none

! Subdomain structure
      type(subdomain_type),intent(inout) :: sub
! Subdomain reactions
      integer,intent(in) ::  lrea
      real(kr),intent(out)::  rea(lrea)

      ! local vars
      character(*),parameter:: routine_name = 'DD_DOWNLOAD_REACTIONS'
      integer :: i
      integer :: ndof

      ! check prerequisites
      if (.not. sub%is_reactions_ready) then
         call error( routine_name,'Reactions are not computed for subdomain: ',sub%isub)
      end if

      ndof = sub%ndof
      ! check dimension
      if (lrea.ne.ndof) then
         call error(routine_name,'Array REA dimension mismatch.')
      end if

      do i = 1,lrea
         rea(i) = sub%rea(i)
      end do

end subroutine

!*************************************************
subroutine dd_add_interior_solution(sub, sol,lsol)
!*************************************************
! Subroutine for adding interior solution (stored inside structure) to subdomain solution
      use module_utils
      implicit none

! Subdomain structure
      type(subdomain_type),intent(inout) :: sub
      integer,intent(in) ::    lsol
      real(kr),intent(inout)::  sol(lsol)

      ! local vars
      character(*),parameter:: routine_name = 'DD_ADD_INTERIOR_SOLUTION'
      integer :: i, ind
      integer :: ndof, ndofo

      ! check prerequisites
      if (.not. sub%is_interior_solution_loaded) then
         call error(routine_name,'Interior solution not loaded.')
      end if

      ndof  = sub%ndof

      ! check dimension
      if (lsol.ne.ndof) then
         call error(routine_name,'Array SOL dimension mismatch.')
      end if

      ndofo = sub%ndofo
      do i = 1,ndofo
         ind = sub%iovsvn(i)
         sol(ind) = sol(ind) + sub%solo(i)
      end do

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
      if (sub%lifix .gt. 0 .and. sub%lifix .ne. lvec) then
         write(*,*) 'DD_FIX_BC: Vector size mismatch for subdomain ',sub%isub
         call error_exit
      end if

      ! enforce boundary conditions
      if (sub%lifix .gt. 0) then
         where (sub%ifix .ne. 0) vec = sub%fixv
      end if

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
      if (sub%lifix.gt.0) then
         do i = 1,ndofi
            ind = sub%iivsvn(i)
            if (sub%ifix(ind) .ne. 0) then
               vec(i) = 0._kr
            end if
         end do
      end if

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

!***********************************************************
subroutine dd_map_sub_to_subi_int(sub, vec,lvec, veci,lveci)
!***********************************************************
! Subroutine that maps subdomain vector to interface vector
      use module_utils
      implicit none

! Subdomain structure
      type(subdomain_type),intent(in) :: sub

      integer,intent(in)  ::  lvec
      integer,intent(in) ::   vec(lvec)

      integer,intent(in)  ::  lveci
      integer,intent(out) ::  veci(lveci)

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
         write(*,*) 'DD_MAP_SUB_TO_SUBO: Interface not loaded for subdomain: ',sub%isub
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

!*******************************************************
subroutine dd_map_subo_to_sub(sub, veco,lveco, vec,lvec)
!*******************************************************
! Subroutine that maps subdomain interiOr vector to subdomain vector
      use module_utils
      implicit none

! Subdomain structure
      type(subdomain_type),intent(in) :: sub

      integer,intent(in)  ::  lveco
      real(kr),intent(in) ::   veco(lveco)

      integer,intent(in)  ::  lvec
      real(kr),intent(out) ::  vec(lvec)

      ! local vars
      character(*),parameter:: routine_name = 'DD_MAP_SUBO_TO_SUB'
      integer :: i, ind
      integer :: ndofo, ndof

      ! check pre
      if (.not. sub%is_interface_loaded) then
         call error(routine_name,'Interface not loaded for subdomain: ',sub%isub)
      end if

      ! get dimensions
      ndof  = sub%ndof
      ndofo = sub%ndofo

      ! check dimensions
      if (lvec .ne. ndof .or. lveco .ne. ndofo) then
         call error(routine_name,'Vectors size mismatch, isub : ',sub%isub)
      end if

      do i = 1,ndofo
         ind = sub%iovsvn(i)

         vec(ind) = vec(ind) + veco(i)
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
      use module_sm
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
      maski(:) = 0
      masko(:) = 0

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

!*************************************************************
subroutine dd_prepare_schur(sub,comm_self,use_explicit_schurs)
!*************************************************************
! Subroutine for preparing data for computing with reduced problem
      use module_mumps
      use module_utils
      implicit none

! Subdomain structure
      type(subdomain_type),intent(inout) :: sub
      ! communicator
      integer,intent(in) :: comm_self
      ! use explicit Schur?
      logical,intent(in) :: use_explicit_schurs

      ! local vars
      integer :: ndofo, la11, nnza11
      integer :: la, ndof, ndofi, nnza
      integer :: schur_size, lschur
      integer :: mumpsinfo
      integer :: iparallel
      logical :: remove_original

      ! check the prerequisities
      if (.not.sub%is_interface_loaded) then
         write(*,*) 'DD_PREPARE_SCHUR: Interface is not loaded for subdomain:', sub%isub
         call error_exit
      end if
      if (.not.sub%is_matrix_loaded) then
         write(*,*) 'DD_PREPARE_SCHUR: Matrix is not loaded for subdomain:', sub%isub
         call error_exit
      end if
      if (sub%istorage.eq.3 .or. sub%istorage.eq.4) then
         write(*,*) 'DD_PREPARE_SCHUR: Matrix already in block triple format for subdomain:', sub%isub
         return
      end if
      if (.not.(sub%istorage.eq.1 .or. sub%istorage.eq.2)) then
         write(*,*) 'DD_PREPARE_SCHUR: Matrix not in single block triple format, subdomain:', sub%isub,&
                    'type:',sub%istorage
         call error_exit
      end if

      ! Load matrix to MUMPS
      if (use_explicit_schurs) then
         ndof  = sub%ndof
         ndofi = sub%ndofi
         nnza  = sub%nnza
         la    = sub%la

         ! Initialize MUMPS
         call mumps_init(sub%mumps_subdomain_matrix,comm_self,sub%matrixtype)
         ! Level of information from MUMPS
         if (debug) then
            mumpsinfo = 2
         else
            mumpsinfo = 0
         end if
         call mumps_set_info(sub%mumps_subdomain_matrix,mumpsinfo)

         call mumps_load_triplet_centralized(sub%mumps_subdomain_matrix, ndof, nnza, &
                                             sub%i_a_sparse, sub%j_a_sparse, sub%a_sparse, la)

         schur_size = ndofi
         call mumps_set_schur(sub%mumps_subdomain_matrix, sub%iivsvn, schur_size)

         sub%lschur1 = schur_size
         sub%lschur2 = schur_size
         allocate(sub%schur(schur_size,schur_size))

         lschur = schur_size*schur_size

         ! Associate Schur complement
         call mumps_assoc_schur(sub%mumps_subdomain_matrix, sub%schur, &
                                lschur, schur_size)

         ! Analyze matrix
         iparallel = 1 ! force serial analysis
         call mumps_analyze(sub%mumps_subdomain_matrix,iparallel) 

         ! Factorize matrix 
         call mumps_factorize(sub%mumps_subdomain_matrix) 

         sub%is_explicit_schur_prepared = .true.

      else
         ! block the matrix into four blocks
         remove_original = .false.
         call dd_matrix_tri2blocktri(sub,remove_original)

         ndofo  = sub%ndofo
         nnza11 = sub%nnza11
         la11   = sub%la11
         if (ndofo.gt.0) then
            ! Initialize MUMPS
            call mumps_init(sub%mumps_interior_block,comm_self,sub%matrixtype)
            ! Level of information from MUMPS
            if (debug) then
               mumpsinfo = 2
            else
               mumpsinfo = 0
            end if
            call mumps_set_info(sub%mumps_interior_block,mumpsinfo)

            call mumps_load_triplet_centralized(sub%mumps_interior_block,ndofo,nnza11,&
                                                sub%i_a11_sparse,sub%j_a11_sparse,sub%a11_sparse,nnza11)
            ! Analyze matrix
            iparallel = 1 ! force serial analysis
            call mumps_analyze(sub%mumps_interior_block,iparallel) 
            ! Factorize matrix 
            call mumps_factorize(sub%mumps_interior_block) 
            ! find the size of the factors in interior block
            call mumps_get_factor_size(sub%mumps_interior_block, sub%mumps_interior_block_factor_size )

            sub%is_mumps_interior_active = .true.
         else
            sub%is_mumps_interior_active = .false.
         end if
         sub%is_interior_factorized = .true.
      end if


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
      character(*),parameter:: routine_name = 'DD_LOAD_ARITHMETIC_CONSTRAINTS'
      integer :: icnode, ncnodes, ncdof, nvar
      integer :: pointdof, nnod, inod, ndofn, idofn, indi, ind, ind1, ind2
      real(kr) :: val
      integer ::             lmatrix1, lmatrix2
      real(kr),allocatable :: matrix(:,:)
      integer :: nnz_new

      ! check the prerequisities
      if (.not.allocated(sub%cnodes) .or. .not.sub%is_cnodes_loaded) then
         call error(routine_name, 'Array for cnodes not ready for subdomain ',sub%isub)
      end if

      ! get number of coarse nodes
      ncnodes = sub%ncnodes

      ! generate arithmetic averages on coarse nodes of prescribed type (e.g. edges)
      do icnode = 1,ncnodes
         if (sub%cnodes(icnode)%itype .eq. itype) then

            ! get number of constraints on an arithmetic constraint
            nvar = sub%cnodes(icnode)%nvar

            ! maximal number of dofs at a node
            ncdof = 0
            nnod = sub%cnodes(icnode)%nnod
            do inod = 1,nnod
               indi = sub%cnodes(icnode)%insin(inod)
               ind  = sub%iin(indi)

               ndofn = sub%nndf(ind)
               if (ndofn.gt.ncdof) then
                  ncdof = ndofn
               end if
            end do

            lmatrix1 = ncdof
            lmatrix2 = nvar
            allocate(matrix(lmatrix1,lmatrix2))
            matrix(:,:) = 0._kr

            pointdof = 0
            do inod = 1,nnod
               indi = sub%cnodes(icnode)%insin(inod)
               ind  = sub%iin(indi)

               val = 1._kr

               ndofn = sub%nndf(ind)
               do idofn = 1,ndofn
                  ind1 = idofn
                  ind2 = pointdof + idofn
                  ! check indices
                  if (ind1.gt.lmatrix1) then
                     call error(routine_name,'row index out of bounds for matrix of constraints')
                  end if
                  if (ind2.gt.lmatrix2) then
                     call error(routine_name,'column index out of bounds for matrix of constraints')
                  end if

                  matrix(ind1,ind2) = val
               end do
               pointdof = pointdof + ndofn
            end do

            nnz_new = nvar ! number of new nonzeroes equal NVAR thanks to the structure
            call dd_append_cnode_constraits(sub%cnodes(icnode),matrix,lmatrix1,lmatrix2,nnz_new)
            deallocate(matrix)

            ! mark the coarse node as used
            sub%cnodes(icnode)%used = .true.
         end if
      end do

end subroutine

!**************************************************************************
subroutine dd_load_adaptive_constraints(sub,gglob,cadapt,lcadapt1,lcadapt2)
!**************************************************************************
! Subroutine for assemblage of matrix
      use module_utils
      implicit none
! Subdomain structure
      type(subdomain_type),intent(inout) :: sub
      integer,intent(in) :: gglob

      integer,intent(in) :: lcadapt1, lcadapt2
      real(kr),intent(in) :: cadapt(lcadapt1,lcadapt2)

! local variables
      character(*),parameter:: routine_name = 'DD_LOAD_ADAPTIVE_CONSTRAINTS'

      ! local vars
      integer :: ind_loc, nvarglb
      integer :: i, j, indiv
      integer ::             lmatrix1, lmatrix2
      real(kr),allocatable :: matrix(:,:)
      integer :: nnz_new

      ! check the prerequisities
      if (.not.allocated(sub%cnodes) .or. .not.sub%is_cnodes_loaded) then
         call error(routine_name,'Array for cnodes not ready.')
      end if

      ! find local (subdomain) index of the glob from its global number
      call get_index(gglob,sub%cnodes%global_cnode_number,sub%ncnodes,ind_loc)
      if (ind_loc.le.0) then
         call error(routine_name,' Index of local glob not found for global ',gglob)
      end if

      ! prepare space for these constraints in the structure
      nvarglb = sub%cnodes(ind_loc)%nvar 

      if (nvarglb.gt.lcadapt1) then
         call error(routine_name,' Number of variables at glob seems larger than interface size of subdomain for glob',gglob)
      end if

      ! append selected constraints to the global structure
      lmatrix1 = lcadapt2
      lmatrix2 = nvarglb
      allocate(matrix(lmatrix1,lmatrix2))
      matrix(:,:) = 0._kr

      ! copy transposed constraints
      do i = 1,nvarglb
         indiv = sub%cnodes(ind_loc)%ivsivn(i)

         do j = 1,lcadapt2
            matrix(j,i) = cadapt(indiv,j)
         end do
      end do

      if (debug) then
         call info(routine_name,' Loading adaptive matrix of globs of subdomain ',sub%isub)
         call info(routine_name,' local glob #',ind_loc)
         do i = 1,lmatrix1
            write(*,'(100f15.6)') (matrix(i,j),j = 1,lmatrix2)
         end do
      end if

      nnz_new = lcadapt2 * nvarglb
      call dd_append_cnode_constraits(sub%cnodes(ind_loc),matrix,lmatrix1,lmatrix2,nnz_new)
      deallocate(matrix)

      ! mark the coarse node as used
      sub%cnodes(ind_loc)%used = .true.

end subroutine

!*********************************************
subroutine dd_load_user_constraints(sub,itype)
!*********************************************
! Subroutine for creating constraints based on user supplied data
      use module_utils
      implicit none
! Subdomain structure
      type(subdomain_type),intent(inout) :: sub

      integer,intent(in) :: itype ! type of globs (2 - edges, 1 - faces)

      ! local vars
      character(*),parameter:: routine_name = 'DD_LOAD_USER_CONSTRAINTS'
      integer :: icnode, ncnodes, ncdof, nvar
      integer :: pointdof, nnod, inod, ndofn, idofn, indi, ind, ind1, ind2, &
                 iconstr, nuser_constraints
      integer :: nnz_new
      real(kr) :: val
      integer ::             lmatrix1, lmatrix2
      real(kr),allocatable :: matrix(:,:)

      ! check the prerequisities
      if (.not.sub%is_user_constraints_loaded) then
         call error(routine_name, 'User constraints not loaded for subdomain ',sub%isub)
      end if
      if (.not.allocated(sub%cnodes) .or. .not.sub%is_cnodes_loaded) then
         call error(routine_name, 'Array for cnodes not ready for subdomain ',sub%isub)
      end if

      ! get number of coarse nodes
      ncnodes = sub%ncnodes

      ! generate user constraints on coarse nodes of prescribed type (e.g. edges)
      do icnode = 1,ncnodes
         if (sub%cnodes(icnode)%itype .eq. itype) then

            ! get number of constraints on an arithmetic constraint
            nvar = sub%cnodes(icnode)%nvar

            ! maximal number of dofs at a node per constraint
            ncdof = 0
            nnod = sub%cnodes(icnode)%nnod
            do inod = 1,nnod
               indi = sub%cnodes(icnode)%insin(inod)
               ind  = sub%iin(indi)

               ndofn = sub%nndf(ind)
               if (ndofn.gt.ncdof) then
                  ncdof = ndofn
               end if
            end do
            nuser_constraints = sub%luser_constraints1

            lmatrix1 = ncdof * nuser_constraints
            lmatrix2 = nvar

            allocate(matrix(lmatrix1,lmatrix2))
            matrix(:,:) = 0._kr

            do iconstr = 1,nuser_constraints
               pointdof = 0
               do inod = 1,nnod
                  indi = sub%cnodes(icnode)%insin(inod)
                  ind  = sub%iin(indi)

                  val = sub%user_constraints(iconstr,ind)

                  ndofn = sub%nndf(ind)
                  do idofn = 1,ndofn
                     ind1 = (iconstr-1)*ncdof + idofn
                     ind2 = pointdof + idofn
                     ! check indices
                     if (ind1.gt.lmatrix1) then
                        call error(routine_name,'row index out of bounds for matrix of constraints')
                     end if
                     if (ind2.gt.lmatrix2) then
                        call error(routine_name,'column index out of bounds for matrix of constraints')
                     end if

                     matrix(ind1,ind2) = val
                  end do
                  pointdof = pointdof + ndofn
               end do
            end do

            !print *,'indices:', sub%iin(sub%cnodes(icnode)%insin)
            !print *,'matrix:'  
            !do i = 1,lmatrix1
            !   print *, matrix(i,:)
            !end do

            nnz_new = count(matrix.ne.0._kr)   ! number of new nonzeroes equal NVAR thanks to the structure
            if (nnz_new.gt.0) then
               call dd_append_cnode_constraits(sub%cnodes(icnode),matrix,lmatrix1,lmatrix2,nnz_new)
            else
               call warning(routine_name, 'User constraints are plain zeros on coarse node ',&
                            sub%cnodes(icnode)%global_cnode_number)
            end if
            deallocate(matrix)

            ! mark the coarse node as used
            sub%cnodes(icnode)%used = .true.
         end if

      end do

end subroutine

!*************************************************************************
subroutine dd_append_cnode_constraits(cnode,matrix,lmatrix1,lmatrix2, nnz)
!*************************************************************************
! Subroutine for inquiring sizes to allocate for number of adaptive averages
      use module_utils
      implicit none
! Coarse node structure
      type(cnode_type),intent(inout) :: cnode

      ! matrix of appended constraints
      integer,intent(in) :: lmatrix1, lmatrix2
      real(kr),intent(in) :: matrix(lmatrix1,lmatrix2)

      ! number of nonzeros in matrix
      integer,intent(in) :: nnz

      ! local vars
      character(*),parameter:: routine_name = 'DD_APPEND_CNODE_CONSTRAITS'
      integer :: nnz_old, ncdof_old 
      integer::              lmatrix_old1, lmatrix_old2
      real(kr),allocatable :: matrix_old(:,:)
      integer::              lmatrix_new1, lmatrix_new2

      ! check the prerequisities
      if (allocated(cnode%matrix) .and. cnode%lmatrix1 .ne. cnode%ncdof) then
         call error(routine_name, 'Matrix size mismatch for existing coarse node',cnode%global_cnode_number)
      end if
      if (allocated(cnode%matrix) .and. cnode%lmatrix2 .ne. lmatrix2) then
         call error(routine_name, 'Cannot append matrix to existing coarse node',cnode%global_cnode_number)
      end if
      if (cnode%nvar .ne. lmatrix2) then
         call error(routine_name, 'Matrix size mismatch for coarse node',cnode%global_cnode_number)
      end if
      !if (count(matrix.ne.0._kr) .ne. nnz) then
      !   call info( routine_name, 'nnz = ',nnz )
      !   call info( routine_name, 'count(matrix.ne.0._kr) = ',count(matrix.ne.0._kr) )
      !   call error(routine_name, 'Matrix number of nonzeros mismatch for coarse node',cnode%global_cnode_number)
      !end if

      ncdof_old = cnode%ncdof
      nnz_old   = cnode%nnz 

      ! store existing matrix
      lmatrix_old1 = cnode%lmatrix1
      lmatrix_old2 = cnode%lmatrix2
      if (allocated(cnode%matrix)) then
         allocate(matrix_old(lmatrix_old1,lmatrix_old2))
         matrix_old = cnode%matrix
         deallocate(cnode%matrix)
      end if

      lmatrix_new1 = lmatrix_old1 + lmatrix1
      lmatrix_new2 = lmatrix2
      allocate(cnode%matrix(lmatrix_new1,lmatrix_new2))

      ! copy back existing matrix
      if (allocated(matrix_old)) then
         cnode%matrix(1:lmatrix_old1,:) = matrix_old
         deallocate(matrix_old)
      end if
      ! append new matrix to the bottom
      cnode%matrix(lmatrix_old1+1:lmatrix_new1,:) = matrix
      cnode%lmatrix1 = lmatrix_new1
      cnode%lmatrix2 = lmatrix_new2

      ! update properties
      cnode%ncdof = ncdof_old + lmatrix1
      cnode%nnz   = nnz_old   + nnz
end subroutine

!*************************************************
subroutine dd_orthogonalize_constraints(sub,itype)
!*************************************************
! Subroutine for assemblage of matrix
      use module_utils
      implicit none
! Subdomain structure
      type(subdomain_type),intent(inout) :: sub

      integer,intent(in) :: itype ! type of globs (2 - edges, 1 - faces)

! local variables
      character(*),parameter:: routine_name = 'DD_ORTHOGONALIZE_CONSTRAINTS'

      integer::              lconstraints1, lconstraints2
      real(kr),allocatable :: constraints(:,:)

! LAPACK variables
      integer             :: lipiv
      integer,allocatable ::  ipiv(:)
      integer              :: lwork
      real(kr),allocatable ::  work(:)
      integer              :: ltau
      real(kr),allocatable ::  tau(:)
      integer             :: lapack_info, ldconstraints

      real(kr) :: normval, thresh_diag

      integer :: lmatrix2, nvar, icnode
      integer :: lmatrix1_old, ncdof_old, ncnodes
      integer :: i, j
      integer :: limit_loop
      integer ::  nnz, nnz_old, nvalid
      real(kr) :: val

      ! check the prerequisities
      if (.not.allocated(sub%cnodes) .or. .not.sub%is_cnodes_loaded) then
         call error(routine_name,'Array for cnodes not ready.')
      end if

      ! get number of coarse nodes
      ncnodes = sub%ncnodes

      ! generate arithmetic averages on coarse nodes of prescribed type (e.g. edges)
      do icnode = 1,ncnodes
         if (sub%cnodes(icnode)%itype .eq. itype) then

            nvar         = sub%cnodes(icnode)%nvar
            ncdof_old    = sub%cnodes(icnode)%ncdof
            nnz_old      = sub%cnodes(icnode)%nnz

            lmatrix1_old = sub%cnodes(icnode)%lmatrix1
            lmatrix2     = sub%cnodes(icnode)%lmatrix2

            ! perform checks
            if (.not. allocated(sub%cnodes(icnode)%matrix)) then
               call error(routine_name, 'Matrix not allocated, nothing to regularize.',sub%cnodes(icnode)%global_cnode_number)
            end if
            if (lmatrix1_old .ne. ncdof_old) then
               call error(routine_name, 'Constraints size mismatch for glob',sub%cnodes(icnode)%global_cnode_number)
            end if
            if (lmatrix2 .ne. nvar) then
               call error(routine_name, 'Constraints size mismatch for glob',sub%cnodes(icnode)%global_cnode_number)
            end if

            ! space for transposed constraints
            lconstraints1 = lmatrix2
            lconstraints2 = lmatrix1_old
            allocate(constraints(lconstraints1,lconstraints2))
            constraints = transpose(sub%cnodes(icnode)%matrix)

            !write (*,*) 'CONSTRAINTS before QR'
            !do i = 1,lconstraints1
            !   write(*,*) (constraints(i,j),j = 1,lconstraints2)
            !end do

            ! perform QR decomposition of AVG by LAPACK
            ! Prepare array for permutations
            lipiv = lconstraints2
            allocate(ipiv(lipiv))
            ipiv = 0
            ! prepare other LAPACK arrays
            ltau = lconstraints1
            allocate(tau(ltau))
            lwork = 3*lconstraints2 + 1
            allocate(work(lwork))

            ldconstraints = max(1,lconstraints1)
            ! QR decomposition
            call DGEQP3( lconstraints1, lconstraints2, constraints, ldconstraints, ipiv, tau, work, lwork, lapack_info )

            !write (*,*) 'constraints after QR factorization'
            !do i = 1,lconstraints1
            !   write(*,*) (constraints(i,j),j = 1,lconstraints2)
            !end do
            !write (*,*) 'IPIV after QR factorization'
            !write(*,*) (ipiv(j),j = 1,lconstraints2)

            ! determine number of columns to use
            ! threshold of 1% of maximal norm
            nvalid = 0
            if (lconstraints1.gt.0.and.lconstraints2.gt.0) then
               normval = abs(constraints(1,1))
               if (normval.gt.numerical_zero) then
                  thresh_diag = 0.01_kr * normval
                  limit_loop = min(lconstraints1,lconstraints2)
                  do i = 1,limit_loop
                     if (abs(constraints(i,i)) .lt. thresh_diag) then
                        exit
                     else
                        nvalid = i
                     end if
                  end do
               end if
            end if

            !write (*,*) 'Number of constraints to really use nvalid',nvalid

            ! construct Q in constraints array
            call DORGQR( lconstraints1, nvalid, nvalid, constraints, ldconstraints, tau, work, lwork, lapack_info )

            deallocate(work)
            deallocate(tau)
            deallocate(ipiv)

            !write (*,*) 'constraints contains Q'
            !do i = 1,lconstraints1
            !   write(*,*) (constraints(i,j),j = 1,nvalid)
            !end do

            if (debug) then
               if (nvalid.lt.lconstraints2) then
                  call warning(routine_name,' Almost linearly dependent constraints on glob ',&
                               sub%cnodes(icnode)%global_cnode_number)
                  call warning(routine_name,' Number of constrains reduced to ',nvalid)
               end if
            end if

            ! copy transposed selected regularized constraints to the global structure
            sub%cnodes(icnode)%ncdof = nvalid

            sub%cnodes(icnode)%lmatrix1 = nvalid
            sub%cnodes(icnode)%lmatrix2 = lmatrix2
            deallocate(sub%cnodes(icnode)%matrix)
            allocate(sub%cnodes(icnode)%matrix(sub%cnodes(icnode)%lmatrix1,sub%cnodes(icnode)%lmatrix2))
            sub%cnodes(icnode)%matrix(:,:) = 0._kr

            nnz = 0
            do i = 1,nvalid
               do j = 1,nvar

                  val = constraints(j,i)

                  if (abs(val).gt.0._kr) then
                     sub%cnodes(icnode)%matrix(i,j) = val
                     nnz = nnz + 1
                  end if
               end do
            end do
            sub%cnodes(icnode)%nnz = nnz

            deallocate(constraints)

         end if
      end do
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
      character(*),parameter:: routine_name = 'DD_CONSTRUCT_CNODES'
      integer ::             ncorner
      integer ::             nglob
      integer ::             ncnodes, lcnodes
      integer ::             nnodi, ndofn
      integer ::             icnode, igcnode
      integer ::             inodc, indnode, indinode, i, nvar, iglob, nnodgl, indn, indin
      integer ::             lxyz

      integer ::             lkdofi
      integer,allocatable ::  kdofi(:)

      ! check the prerequisities
      if (.not.sub%is_corners_loaded) then
         call error(routine_name,'Corners not loaded for subdomain:', sub%isub)
      end if
      if (.not.sub%is_globs_loaded) then
         call error(routine_name,'Globs not loaded for subdomain:', sub%isub)
      end if

      ! determine number of coarse nodes
      ncorner = sub%ncorner
      nglob   = sub%nglob
      nnodi   = sub%nnodi
      ncnodes = ncorner + nglob
      sub%ncnodes = ncnodes

      lcnodes = ncnodes
      allocate(sub%cnodes(lcnodes))

      ! prepare KDOFI
      lkdofi = nnodi
      allocate(kdofi(lkdofi))
      if (lkdofi.gt.0) then
         kdofi(1) = 0
         do i = 2,nnodi
            indn = sub%iin(i-1)
            ndofn = sub%nndf(indn)

            kdofi(i) = kdofi(i-1) + ndofn
         end do
      end if

      ! set counter
      icnode = 0

      ! copy corners
      do inodc = 1,ncorner
         icnode = icnode + 1

         indinode = sub%icnsin(inodc)
         indnode  = sub%iin(indinode)

         ! type of coarse node - corner
         sub%cnodes(icnode)%itype = 3

         ! global number
         igcnode =  sub%global_corner_number(inodc)
         sub%cnodes(icnode)%global_cnode_number  = igcnode
         ! number of nodes where it maps from
         sub%cnodes(icnode)%nnod = 1

         ! number of variables it maps from 
         nvar = sub%nndf(indnode)
         sub%cnodes(icnode)%nvar = nvar

         ! fill coordinates
         lxyz = sub%ndim
         sub%cnodes(icnode)%lxyz = lxyz
         allocate(sub%cnodes(icnode)%xyz(lxyz))
         sub%cnodes(icnode)%xyz = sub%xyz(indnode,1:lxyz)

         ! fill coarse node nodes
         allocate(sub%cnodes(icnode)%insin(1))
         sub%cnodes(icnode)%insin(1) = indinode

         ! fill coarse node variables
         allocate(sub%cnodes(icnode)%ivsivn(nvar))
         do i = 1,nvar
            sub%cnodes(icnode)%ivsivn(i) = kdofi(indinode) + i
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
         sub%cnodes(icnode)%xyz(:) = 0._kr
         ! create averaged coordinates
         do i = 1,nnodgl
            indin = sub%ignsin(iglob,i)
            indn  = sub%iin(indin)
            sub%cnodes(icnode)%xyz(1:lxyz) = sub%cnodes(icnode)%xyz(1:lxyz) + sub%xyz(indn,1:lxyz)
         end do
         sub%cnodes(icnode)%xyz(1:lxyz) = sub%cnodes(icnode)%xyz(1:lxyz) / nnodgl
         ! debug
         !print *,'coords of glob'
         !print *,sub%cnodes(icnode)%xyz(:)
         !call flush(6)

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

      deallocate(kdofi)

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
      character(*),parameter:: routine_name = 'DD_PREPARE_C'
      integer ::             nnzc
      integer ::             lc
      integer,allocatable ::  i_c_sparse(:)
      integer,allocatable ::  j_c_sparse(:)
      real(kr),allocatable ::   c_sparse(:)

      integer ::             lindrowc
      integer,allocatable ::  indrowc(:)

      integer ::             lkdof
      integer,allocatable ::  kdof(:)

      integer :: nnod
      integer :: inod,& 
                 nconstr, icdof, icn, inzc, irowc, ivar, ncdof, &
                 ncnodes, nrowc, nvar, lmatrix1, lmatrix2
      real(kr) :: val

      ! check the prerequisities
      if (.not.sub%is_cnodes_embedded) then
         call error(routine_name, 'Coarse nodes not ready for subdomain:', sub%isub)
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
      if (lkdof.gt.0) then
         kdof(1) = 0
      end if
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
            ! copy the matrix of constraints on glob into sparse triplet of subdomain matrix C
            ! row by row
            nvar      = sub%cnodes(icn)%nvar
            ncdof     = sub%cnodes(icn)%ncdof
            lmatrix1  = sub%cnodes(icn)%lmatrix1
            lmatrix2  = sub%cnodes(icn)%lmatrix2
            if (nvar.ne.lmatrix2) then
               call info(routine_name, 'nvar', nvar)
               call info(routine_name, 'lmatrix2', lmatrix2)
               call error(routine_name, 'Second matrix dimension does not match for subdomain', sub%isub)
            end if
            if (ncdof.ne.lmatrix1) then
               call error(routine_name, 'First matrix dimension does not match for subdomain', sub%isub)
            end if
            do icdof = 1,ncdof
               irowc = irowc + 1

               do ivar = 1,nvar
                  val = sub%cnodes(icn)%matrix(icdof,ivar)
                  if (abs(val).gt.0._kr) then
                     inzc  = inzc  + 1
                     ! check length
                     if (inzc.gt.lc) then
                        call error(routine_name,'out of bounds for matrix C', sub%isub)
                     end if

                     i_c_sparse(inzc) = irowc
                     j_c_sparse(inzc) = sub%cnodes(icn)%ivsivn(ivar)
                     c_sparse(inzc)   = val
                  end if
               end do

               indrowc(irowc) = sub%cnodes(icn)%igcdof(icdof)
            end do
         end if
      end do
      if (inzc.ne.nnzc) then
         call info(routine_name,'nnzc =',nnzc)
         call info(routine_name,'inzc =',inzc)
         call error(routine_name,'some entries from matrix C are missing', sub%isub)
      end if

      ! check matrix bounds
      if (inzc .ne. lc) then
         call error(routine_name, 'Dimension of matrix C mismatch for subdomain', sub%isub)
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
      use module_sm
      use module_utils
      implicit none

! Subdomain structure
      type(subdomain_type),intent(inout) :: sub
      ! communicator
      integer,intent(in) :: comm_self

      ! local vars
      character(*),parameter:: routine_name = 'DD_PREPARE_AUG'
      integer ::  nnzc, nnza, ndof, nconstr
      integer ::  nnzaaug, laaug, ndofaaug
      integer ::  i, iaaug
      integer ::  mumpsinfo, aaugmatrixtype
      integer ::  icol, icoli
      integer ::  jcoli, iconstr, jconstr
      integer ::  ndofi, ndofiaug
      integer :: iparallel

      ! LAPACK related variables
      integer, external :: ILAENV

      integer ::  nb
      integer :: lwork
      real(kr), allocatable :: work(:)
      integer :: lapack_info

      ! check the prerequisities
      !if (sub%is_degenerated) then
      !   return
      !end if

      if (.not.sub%is_interface_loaded) then
         call error(routine_name,'Interface not loaded for subdomain:', sub%isub)
      end if
      if (.not.sub%is_matrix_loaded) then
         call error(routine_name,'Matrix is not loaded for subdomain:', sub%isub)
      end if
      if (.not.sub%is_c_loaded) then
         call error(routine_name,'Matrix of constraints C not loaded for subdomain:', sub%isub)
      end if


      if (sub%is_explicit_schur_prepared) then
         ! if augmented system is already allocated, clear it - this can happen for adaptivity
         ! we are in the explicit world, allocate the large dense matrix
         if (sub%is_aug_factorized) then

            deallocate(sub%aaug_dense)
            sub%laaug_dense1 = 0
            sub%laaug_dense2 = 0

            sub%is_mumps_aug_active = .false.
            sub%is_aug_factorized = .false.
         end if

         ! join the new matrix directly to the structure as
         ! in the unsymmetric case:
         ! S C^T
         ! C  0
         ! in the symmetric case :
         ! \S C^T
         !     0
         ndof      = sub%ndof
         ndofi     = sub%ndofi
         nconstr   = sub%nconstr
         ndofiaug  = ndofi + nconstr
         nnzc      = sub%nnzc

         sub%laaug_dense1 = ndofiaug
         sub%laaug_dense2 = ndofiaug

         allocate(sub%aaug_dense(sub%laaug_dense1,sub%laaug_dense2))
         sub%aaug_dense = 0._kr

         ! copy Schur complement as the (1,1) block of the matrix
         sub%aaug_dense(1:ndofi,1:ndofi) = sub%schur

         ! copy entries of right block of C^T with proper shift in columns
         do i = 1,nnzc
            jconstr = sub%i_c_sparse(i)
            icoli   = sub%j_c_sparse(i)
            if (icoli.gt.sub%ndofi) then
               !print *,'ndofi',sub%ndofi
               !print *,'icoli',icoli
               !print *,'j_c_sparse',sub%j_c_sparse
               call error(routine_name,'out of bounds of interface for subdomain',sub%isub)
            end if

            sub%aaug_dense(icoli,ndofi + jconstr) = sub%c_sparse(i)
         end do
         if      (sub%matrixtype .eq. 0) then
            ! unsymmetric case: apply lower block of C
            do i = 1,nnzc
               iconstr = sub%i_c_sparse(i)
               jcoli   = sub%j_c_sparse(i)

               sub%aaug_dense(ndofi + iconstr, jcoli) = sub%c_sparse(i)
            end do
         end if

         if (sub%is_degenerated) then
            goto 133
         end if

         ! factorize matrix Aaug by LAPACK
         ! Set type of matrix
         if      (sub%matrixtype .eq. 0) then
            ! unsymmetric case, use LU
            sub%laaug_ipiv = sub%laaug_dense1
            allocate(sub%aaug_ipiv(sub%laaug_ipiv))
            call DGETRF(sub%laaug_dense1, sub%laaug_dense2, sub%aaug_dense, sub%laaug_dense1, sub%aaug_ipiv, lapack_info)

         else if (sub%matrixtype .eq. 1 .or. sub%matrixtype .eq. 2) then
            ! in symmetric case, saddle point problem makes the augmented matrix indefinite,
            ! even if the original matrix is SPD, use LDLT
            !SUBROUTINE DSYTRF( UPLO, N, A, LDA, IPIV, WORK, LWORK, INFO )
            sub%laaug_ipiv = sub%laaug_dense1
            allocate(sub%aaug_ipiv(sub%laaug_ipiv))
            nb = ILAENV(1, 'DSYTRF', 'U',  sub%laaug_dense1, 0, 0, 0)
            lwork = sub%laaug_dense1*nb
            allocate(work(lwork))
            call DSYTRF('U', sub%laaug_dense1, sub%aaug_dense, sub%laaug_dense1, sub%aaug_ipiv, work, lwork, lapack_info)
            deallocate(work)
         else 
            call error(routine_name,'Matrixtype not set for subdomain:', sub%isub)
         end if

 133     sub%is_aug_factorized = .true.
         sub%is_aug_dense_active = .true.

      else
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
            if (icoli.gt.sub%ndofi) then
               print *,'ndofi',sub%ndofi
               print *,'icoli',icoli
               print *,'j_c_sparse',sub%j_c_sparse
               call error(routine_name,'out of bounds of interface for subdomain',sub%isub)
            end if
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
            call error(routine_name,'Actual length of augmented matrix does not match for subdomain:', sub%isub)
         end if
         sub%nnzaaug = nnzaaug
         sub%laaug   = laaug

!         call sm_print(100+sub%isub, &
!                       sub%i_aaug_sparse, sub%j_aaug_sparse, sub%aaug_sparse, &
!                       sub%laaug, sub%nnzaaug)

         if (.not. sub%is_degenerated) then

            ! factorize matrix Aaug
            ! Set type of matrix
            if      (sub%matrixtype .eq. 0) then
               ! unsymmetric case:
               aaugmatrixtype = 0
            else if (sub%matrixtype .eq. 1 .or. sub%matrixtype .eq. 2) then
               ! in symmetric case, saddle point problem makes the augmented matrix indefinite,
               ! even if the original matrix is SPD:
               aaugmatrixtype = 2
            else 
               call error(routine_name,'Matrixtype not set for subdomain:', sub%isub)
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

            if (ndofaaug.eq.0) then
               !print *, 'ndof = ', ndof
               !print *, 'nconstr = ', nconstr
               !print *, 'nelem = ', sub%nelem
               !print *, 'nnod = ', sub%nnod
               call error(routine_name, 'This is strange - subdomain not degenerated but still ndofaaug = 0,',&
                          sub%isub)
            end if

            nnzaaug = sub%nnzaaug
            laaug   = sub%laaug
            call mumps_load_triplet_centralized(sub%mumps_aug,ndofaaug,nnzaaug,&
                                                sub%i_aaug_sparse,sub%j_aaug_sparse,sub%aaug_sparse,nnzaaug)
            ! Analyze matrix
            iparallel = 1 ! force serial analysis
            call mumps_analyze(sub%mumps_aug,iparallel) 
            ! Factorize matrix 
            call mumps_factorize(sub%mumps_aug) 
         end if

         sub%is_mumps_aug_active = .true.
         sub%is_aug_factorized = .true.
      end if

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
! for unsymmetric problems, it also solves
! | A^T C^T|| phis_dual | = | 0 |
! | C   0  ||lambda_dual|   | I |
! and 
! Ac = phis_dual^T * A * phis 

      use module_mumps
      use module_utils
      implicit none

! Subdomain structure
      type(subdomain_type),intent(inout) :: sub
      ! store global vector of PHIS instead of the one restricted to interface PHISI
      logical,intent(in) :: keep_global

      ! local vars
      character(*),parameter:: routine_name = 'DD_PREPARE_COARSE'
      integer ::  ndof, nconstr, ndofi, matrixtype
      integer ::  ndofaaug, shift
      integer ::  ndofc
      integer ::  i, j, indphis, nrhs, &
                  indphisstart, indi, icoarsem, lcoarsem
      integer ::  lphisi1, lphisi2, lphis1, lphis2
      integer ::  lphisi_dual1, lphisi_dual2, lphis_dual1, lphis_dual2
      logical :: solve_adjoint

      integer ::             lphis
      real(kr),allocatable :: phis(:)

      integer ::             lac1, lac2
      real(kr),allocatable :: ac(:,:)

      ! check the prerequisities
      !if (sub%is_degenerated) then
      !   return
      !end if
      if (.not.sub%is_aug_factorized) then
         call error(routine_name, 'Augmented matrix is not factorized for subdomain:', sub%isub)
      end if
      if (.not. sub%is_aug_dense_active .and. .not.sub%is_mumps_aug_active) then
         call error(routine_name, 'Augmented matrix solver in not ready for subdomain:', sub%isub)
      end if

      ! if coarse problem is ready, clear it and prepare it again - this can happen for adaptivity
      if (sub%is_coarse_prepared) then
         if (debug) then
            call info(routine_name,'Reinitializing coarse problem for subdomain ', sub%isub)
         end if

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

         if (allocated(sub%phis_dual)) then
            deallocate(sub%phis_dual)
            sub%lphis_dual1 = 0
            sub%lphis_dual2 = 0
            sub%is_phis_dual_prepared   = .false.
         end if

         if (allocated(sub%phisi_dual)) then
            deallocate(sub%phisi_dual)
            sub%lphisi_dual1 = 0
            sub%lphisi_dual2 = 0
            sub%is_phisi_dual_prepared   = .false.
         end if

         deallocate(sub%coarsem)
         sub%lcoarsem = 0
         sub%ndofc = 0
         sub%is_coarse_prepared = .false.
      end if

      ndof     = sub%ndof
      ndofi    = sub%ndofi
      nconstr  = sub%nconstr

      if (sub%is_aug_dense_active) then
         shift = ndofi
         ndofaaug = ndofi + nconstr
      else
         shift = ndof
         ndofaaug  = ndof + nconstr
      end if

      lphis = ndofaaug*nconstr
      allocate(phis(lphis))
      ! zero all entries
      phis(:) = 0._kr

      ! put identity into the block of constraints
      do j = 1,nconstr
         indphis = (j-1)*ndofaaug + shift + j
         phis(indphis) = 1._kr
      end do

      ! solve the system with multiple RHS
      nrhs = nconstr
      solve_adjoint = .false.
      !print *, 'lphis, nrhs, nconstr, ndofi, ndof', lphis, nrhs, nconstr, ndofi, ndof
      call dd_solve_aug(sub, phis,lphis, nrhs, solve_adjoint) 

      !debug
      !write(*,*) 'Subdomain ',sub%isub,' coarse basis functions phis:'
      !do i = 1,ndofaaug
      !   write(*,'(100f13.6)') (phis((j-1)*ndofaaug + i),j = 1,nconstr)
      !end do

      ! Build subdomain coarse matrix by the fact that phis^T*A*phis = -lambda 
      ! this is correct also for unsymmetric problems where Ac_i = phis_dual^T*A*phis = -lambda
      lac1 = nconstr
      lac2 = nconstr
      allocate(ac(lac1,lac2))
      do j = 1,nconstr
         indphisstart  = (j-1)*ndofaaug + shift
         do i = 1,nconstr
            ac(i,j) = -phis(indphisstart + i)
         end do
      end do
      !debug
      !write(*,*) 'Subdomain ',sub%isub,' coarse matrix:'
      !do i = 1,nconstr
      !   write(*,'(100f13.6)') (ac(i,j),j = 1,nconstr)
      !end do

      ! Prepare phis with multiple RHS as dense matrix - stored in an linear array for MUMPS
      if (sub%is_aug_dense_active) then
         if (keep_global) then
             call error(routine_name, 'Cannot give global phis when using dense Schur complement.')
         end if

         ! copy vector phis to interface unknowns and load it to the structure
         lphisi1 = ndofi
         lphisi2 = nconstr
         allocate(sub%phisi(lphisi1,lphisi2))
         sub%lphisi1 = lphisi1
         sub%lphisi2 = lphisi2

         do i = 1,ndofi
            do j = 1,nconstr
               indphis = (j-1)*ndofaaug + i

               sub%phisi(i,j) = phis(indphis)
            end do
         end do
         sub%is_phisi_prepared   = .true.
      else
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
         end if
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

      ! compute dual basis functions from the adjoint problem
      ! A^T
      if (sub%matrixtype .eq. 0) then
         ! zero all entries
         phis(:) = 0._kr

         ! put identity into the block of constraints
         do j = 1,nconstr
            indphis = (j-1)*ndofaaug + shift + j
            phis(indphis) = 1._kr
         end do

         ! solve the system with multiple RHS
         nrhs = nconstr
         solve_adjoint = .true.
         call dd_solve_aug(sub, phis,lphis, nrhs, solve_adjoint) 

         if (sub%is_aug_dense_active) then
            if (keep_global) then
                call error(routine_name, 'Cannot give global phis when using dense Schur complement.')
            end if

            ! copy vector phis to interface unknowns and load it to the structure
            lphisi_dual1 = ndofi
            lphisi_dual2 = nconstr
            allocate(sub%phisi_dual(lphisi_dual1,lphisi_dual2))
            sub%lphisi_dual1 = lphisi_dual1
            sub%lphisi_dual2 = lphisi_dual2

            do i = 1,ndofi
               do j = 1,nconstr
                  indphis = (j-1)*ndofaaug + i

                  sub%phisi_dual(i,j) = phis(indphis)
               end do
            end do
            sub%is_phisi_dual_prepared   = .true.
         else
            if (keep_global) then
               ndof   = sub%ndof
               lphis_dual1 = ndof
               lphis_dual2 = nconstr
               allocate(sub%phis_dual(lphis_dual1,lphis_dual2))
               sub%lphis_dual1 = lphis_dual1
               sub%lphis_dual2 = lphis_dual2
               do i = 1,ndof
                  do j = 1,nconstr
                     indphis = (j-1)*ndofaaug + i
                     sub%phis_dual(i,j) = phis(indphis)
                  end do
               end do
               sub%is_phis_dual_prepared   = .true.
            end if
            ! restrict vector phis to interface unknowns and load it to the structure
            ndofi   = sub%ndofi
            lphisi_dual1 = ndofi
            lphisi_dual2 = nconstr
            allocate(sub%phisi_dual(lphisi_dual1,lphisi_dual2))
            sub%lphisi_dual1 = lphisi_dual1
            sub%lphisi_dual2 = lphisi_dual2
            do i = 1,ndofi
               indi = sub%iivsvn(i)
               do j = 1,nconstr
                  indphis = (j-1)*ndofaaug + indi

                  sub%phisi_dual(i,j) = phis(indphis)
               end do
            end do
            sub%is_phisi_dual_prepared   = .true.
         end if
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
         call error(routine_name, 'Check of coarse matrix length failed for subdomain: ',sub%isub)
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
      integer :: ndof, ndofi, nconstr

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
      ndofi    = sub%ndofi
      nconstr  = sub%nconstr

      if (sub%is_aug_dense_active) then
         ndofaaug = ndofi + nconstr
      else
         ndofaaug = ndof + nconstr
      end if

end subroutine

!**********************************************************
subroutine dd_solve_aug(sub, vec,lvec, nrhs, solve_adjoint)
!**********************************************************
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

      ! should the adjoint system be solved? This solves transposed system. 
      logical,intent(in) ::    solve_adjoint

      ! local vars
      integer ::  ldb
      character(1) :: transa
      integer ::  lapack_info

      if (sub%is_degenerated) then
         return
      end if

      ! check the prerequisities
      if (.not.sub%is_aug_factorized) then
         write(*,*) 'DD_SOLVE_AUG: Augmented matrix in not factorized for subdomain:', sub%isub
         call error_exit
      end if
      if (.not. sub%is_aug_dense_active .and. .not.sub%is_mumps_aug_active) then
         write(*,*) 'DD_SOLVE_AUG: Augmented matrix solver in not ready for subdomain:',sub%isub
         call error_exit
      end if

      if (nrhs > 0) then
         if (mod(lvec,nrhs) .ne. 0) then
            write(*,*) 'DD_SOLVE_AUG: Unclear what the augmented size is:', sub%isub
            call error_exit
         end if
      end if

      if (sub%is_aug_dense_active) then
         ldb = lvec/nrhs
         if (ldb /= sub%ndofi+sub%nconstr) then
            write(*,*) 'DD_SOLVE_AUG: Length of augmented system does not match:', sub%isub
            print *, 'ldb, expected ', ldb, sub%ndofi + sub%nconstr
            call error_exit
         end if
         ! solve the problem using LU or LDLT
         ! Set type of matrix
         if      (sub%matrixtype .eq. 0) then
            ! unsymmetric case, use LU
            if (solve_adjoint) then
                transa = 'T'
            else
                transa = 'N'
            end if
            call DGETRS(transa, sub%laaug_dense1, nrhs, sub%aaug_dense, sub%laaug_dense1, sub%aaug_ipiv, vec, ldb, lapack_info)
         else if (sub%matrixtype .eq. 1 .or. sub%matrixtype .eq. 2) then
            ! in symmetric case, saddle point problem makes the augmented matrix indefinite,
            ! even if the original matrix is SPD, use LDLT
            call DSYTRS('U', sub%laaug_dense1, nrhs, sub%aaug_dense, sub%laaug_dense1, sub%aaug_ipiv, vec, ldb, lapack_info)
         else 
            write(*,*) 'DD_SOLVE_AUG: Matrixtype not set for subdomain:', sub%isub
            call error_exit
         end if
      else
         if (nrhs > 0) then
            if (lvec/nrhs /= sub%ndof+sub%nconstr) then          ! intentional integer divide
               write(*,*) 'DD_SOLVE_AUG: Length of augmented system does not match:', sub%isub
               call error_exit
            end if
         end if

         ! solve the system with multiple RHS
         if (nrhs > 0) then
            call mumps_resolve(sub%mumps_aug,vec,lvec,nrhs,solve_adjoint)
         end if
      end if

end subroutine

!**************************************************
subroutine dd_get_phisi_size(sub, lphisi1, lphisi2)
!**************************************************
! Subroutine for getting size of PHISI matrix
! phisi are coarse space basis functions on subdomain restricted to interface

      use module_utils
      implicit none

! Subdomain structure
      type(subdomain_type),intent(in) :: sub

      ! output size
      integer,intent(out) :: lphisi1
      integer,intent(out) :: lphisi2

      ! local vars
      character(*),parameter:: routine_name = 'DD_GET_PHISI_SIZE'

      ! check the prerequisities
      if (.not.sub%is_phisi_prepared) then
         call error(routine_name, 'PHISI matrix not ready:', sub%isub)
      end if

      ! copy sizes
      lphisi1 = sub%lphisi1
      lphisi2 = sub%lphisi2

end subroutine

!***************************************************
subroutine dd_get_phisi(sub, phisi, lphisi1,lphisi2)
!***************************************************
! Subroutine for getting PHISI matrix
! phisi are coarse space basis functions on subdomain restricted to interface

      use module_utils
      implicit none

! Subdomain structure
      type(subdomain_type),intent(in) :: sub

      ! phisi size
      integer,intent(in) :: lphisi1
      integer,intent(in) :: lphisi2

      ! output - dense matrix of coarse basis functions
      real(kr),intent(out) :: phisi(lphisi1,lphisi2)

      ! local vars
      character(*),parameter:: routine_name = 'DD_GET_PHISI'
      integer :: i,j

      ! check the prerequisities
      if (.not.sub%is_phisi_prepared) then
         call error(routine_name, 'PHISI matrix not ready:', sub%isub)
      end if

      ! check dimensions
      if (lphisi1 .ne. sub%lphisi1) then
         call error(routine_name, 'PHISI matrix first dimension mismatch:', sub%isub)
      end if
      if (lphisi2 .ne. sub%lphisi2) then
         call error(routine_name, 'PHISI matrix second dimension mismatch:', sub%isub)
      end if

      ! checks OK, copy matrix
      do j = 1,lphisi2
         do i = 1,lphisi1
            
            phisi(i,j) = sub%phisi(i,j)
         end do
      end do

end subroutine

!*****************************************************
subroutine dd_phisi_apply(sub, vec1,lvec1, vec2,lvec2)
!*****************************************************
! Subroutine for multiplication of vector VEC1 by PHISI matrix
! vec2 = phisi * vec1 + vec2
! phisi are coarse space basis functions on subdomain restricted to interface

      use module_utils
      implicit none

! Subdomain structure
      type(subdomain_type),intent(in) :: sub

      ! input vector
      integer,intent(in) :: lvec1
      real(kr),intent(in) :: vec1(lvec1)
      ! output vector
      integer,intent(in) ::  lvec2
      real(kr),intent(out) :: vec2(lvec2)

      ! local vars
      character(*),parameter:: routine_name = 'DD_PHISI_APPLY'

      ! BLAS vars
      character(1) :: TRANS
      integer :: M, N, LDA, INCX, INCY
      real(kr) :: alpha, beta

      ! check the prerequisities
      if (.not.sub%is_phisi_prepared) then
         call error(routine_name, 'PHISI matrix not ready for sub: ', sub%isub)
      end if
      ! check dimensions
      if (lvec1 .ne. sub%lphisi2 .or. lvec2 .ne. sub%lphisi1) then 
         call error(routine_name, 'Dimensions mismatch for sub:', sub%isub)
      end if

      ! checking done, perform multiply by BLAS
      TRANS = 'N'
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

!****************************************************
subroutine dd_phis_apply(sub, vec1,lvec1, vec2,lvec2)
!****************************************************
! Subroutine for multiplication of vector VEC1 by PHIS matrix
! vec2 = phis * vec1 + vec2
! phis are coarse space basis functions on subdomain

      use module_utils
      implicit none

! Subdomain structure
      type(subdomain_type),intent(in) :: sub

      ! input vector
      integer,intent(in) :: lvec1
      real(kr),intent(in) :: vec1(lvec1)
      ! output vector
      integer,intent(in) ::  lvec2
      real(kr),intent(out) :: vec2(lvec2)

      ! local vars
      character(*),parameter:: routine_name = 'DD_PHIS_APPLY'

      ! BLAS vars
      character(1) :: TRANS
      integer :: M, N, LDA, INCX, INCY
      real(kr) :: alpha, beta

      ! check the prerequisities
      if (sub%is_degenerated) then
         return
      end if
      if (.not.sub%is_phis_prepared) then
         call error(routine_name, 'PHIS matrix not ready for sub: ', sub%isub)
      end if
      ! check dimensions
      if (lvec1 .ne. sub%lphis2 .or. lvec2 .ne. sub%lphis1) then 
         call error(routine_name, 'Dimensions mismatch for sub:', sub%isub)
      end if

      ! checking done, perform multiply by BLAS
      TRANS = 'N'
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

!**********************************************************
subroutine dd_phisi_dual_apply(sub, vec1,lvec1, vec2,lvec2)
!**********************************************************
! Subroutine for multiplication of vector VEC1 by PHISI matrix
! vec2 = phisi_dual^T * vec1 + vec2
! phisi_dual are dual coarse space basis functions on subdomain restricted to interface

      use module_utils
      implicit none

! Subdomain structure
      type(subdomain_type),intent(in) :: sub

      ! input vector
      integer,intent(in) :: lvec1
      real(kr),intent(in) :: vec1(lvec1)
      ! output vector
      integer,intent(in) ::  lvec2
      real(kr),intent(out) :: vec2(lvec2)

      ! local vars
      character(*),parameter:: routine_name = 'DD_PHISI_DUAL_APPLY'

      ! BLAS vars
      character(1) :: TRANS
      integer :: M, N, LDA, INCX, INCY
      real(kr) :: alpha, beta

      ! check the prerequisities
      if (sub%is_degenerated) then
         return
      end if
      if ((sub%matrixtype .eq. 1 .or. sub%matrixtype .eq. 2) ) then
          
         if ( .not.sub%is_phisi_prepared) then
            call error(routine_name, 'PHISI matrix not ready for sub: ', sub%isub)
         end if
         ! check dimensions
         if (lvec1 .ne. sub%lphisi1 .or. lvec2 .ne. sub%lphisi2) then 
            call error(routine_name, 'Dimensions mismatch for sub:', sub%isub)
         end if
      else
         if (.not.sub%is_phisi_dual_prepared) then
            call error(routine_name, 'PHISI_DUAL matrix not ready for sub: ', sub%isub)
         end if
         ! check dimensions
         if (lvec1 .ne. sub%lphisi_dual1 .or. lvec2 .ne. sub%lphisi_dual2) then 
            call error(routine_name, 'Dimensions mismatch for sub:', sub%isub)
         end if
      end if

      ! checking done, perform multiply by BLAS
      TRANS = 'T'
      ALPHA = 1._kr
      INCX = 1
      BETA = 1._kr ! sum second vector
      INCY = 1
      if ((sub%matrixtype .eq. 1 .or. sub%matrixtype .eq. 2) ) then
         M = sub%lphisi1
         N = sub%lphisi2
         LDA = max(1,M)
         if (kr.eq.8) then
            ! double precision
            call DGEMV(TRANS,M,N,ALPHA,sub%phisi,LDA,vec1,INCX,BETA,vec2,INCY)
         else if (kr.eq.4) then
            ! single precision
            call SGEMV(TRANS,M,N,ALPHA,sub%phisi,LDA,vec1,INCX,BETA,vec2,INCY)
         end if
      else
         M = sub%lphisi_dual1
         N = sub%lphisi_dual2
         LDA = max(1,M)
         if (kr.eq.8) then
            ! double precision
            call DGEMV(TRANS,M,N,ALPHA,sub%phisi_dual,LDA,vec1,INCX,BETA,vec2,INCY)
         else if (kr.eq.4) then
            ! single precision
            call SGEMV(TRANS,M,N,ALPHA,sub%phisi_dual,LDA,vec1,INCX,BETA,vec2,INCY)
         end if
      end if

end subroutine

!*********************************************************
subroutine dd_phis_dual_apply(sub, vec1,lvec1, vec2,lvec2)
!*********************************************************
! Subroutine for multiplication of vector VEC1 by PHIS matrix
! vec2 = phis * vec1 + vec2
! phis are coarse space basis functions on subdomain

      use module_utils
      implicit none

! Subdomain structure
      type(subdomain_type),intent(in) :: sub

      ! input vector
      integer,intent(in) :: lvec1
      real(kr),intent(in) :: vec1(lvec1)
      ! output vector
      integer,intent(in) ::  lvec2
      real(kr),intent(out) :: vec2(lvec2)

      ! local vars
      character(*),parameter:: routine_name = 'DD_PHIS_DUAL_APPLY'

      ! BLAS vars
      character(1) :: TRANS
      integer :: M, N, LDA, INCX, INCY
      real(kr) :: alpha, beta

      ! check the prerequisities
      if (sub%is_degenerated) then
         return
      end if
      if ((sub%matrixtype .eq. 1 .or. sub%matrixtype .eq. 2) ) then
         if (.not.sub%is_phis_prepared) then
            call error(routine_name, 'PHIS matrix not ready for sub: ', sub%isub)
         end if
      ! check dimensions
         if (lvec1 .ne. sub%lphis1 .or. lvec2 .ne. sub%lphis2) then 
            call error(routine_name, 'Dimensions mismatch for sub:', sub%isub)
         end if
      else
         if (.not.sub%is_phis_dual_prepared) then
            call error(routine_name, 'PHIS_DUAL matrix not ready for sub: ', sub%isub)
         end if
      ! check dimensions
         if (lvec1 .ne. sub%lphis_dual1 .or. lvec2 .ne. sub%lphis_dual2) then 
            call error(routine_name, 'Dimensions mismatch for sub:', sub%isub)
         end if

      end if

      ! checking done, perform multiply by BLAS
      TRANS = 'T'
      ALPHA = 1._kr
      INCX = 1
      BETA = 1._kr ! sum second vector
      INCY = 1
      if ((sub%matrixtype .eq. 1 .or. sub%matrixtype .eq. 2) ) then
         M = sub%lphis1
         N = sub%lphis2
         LDA = max(1,M)
         if (kr.eq.8) then
            ! double precision
            call DGEMV(TRANS,M,N,ALPHA,sub%phis,LDA,vec1,INCX,BETA,vec2,INCY)
         else if (kr.eq.4) then
            ! single precision
            call SGEMV(TRANS,M,N,ALPHA,sub%phis,LDA,vec1,INCX,BETA,vec2,INCY)
         end if
      else
         M = sub%lphis_dual1
         N = sub%lphis_dual2
         LDA = max(1,M)
         if (kr.eq.8) then
            ! double precision
            call DGEMV(TRANS,M,N,ALPHA,sub%phis_dual,LDA,vec1,INCX,BETA,vec2,INCY)
         else if (kr.eq.4) then
            ! single precision
            call SGEMV(TRANS,M,N,ALPHA,sub%phis_dual,LDA,vec1,INCX,BETA,vec2,INCY)
         end if
      end if

end subroutine

!*******************************************
subroutine dd_get_coarsem_size(sub,lcoarsem)
!*******************************************
! Subroutine for obtaining SIZE of coarse matrix of subdomain
      use module_utils
      implicit none

! sub structure for actual subdomain
      type(subdomain_type),intent(in) :: sub

      ! coarse matrix length
      integer,intent(out)  :: lcoarsem

      ! local vars
      character(*),parameter:: routine_name = 'DD_GET_COARSEM_SIZE'

      ! check the prerequisities
      if (.not.sub%is_coarse_prepared) then
         call error(routine_name, 'Coarse matrix not ready:', sub%isub)
      end if

      ! if checks are OK, copy data
      lcoarsem = sub%lcoarsem

end subroutine
      
!**********************************************
subroutine dd_get_coarsem(sub,coarsem,lcoarsem)
!**********************************************
! Subroutine for obtaining coarse matrix of subdomain
      use module_utils
      implicit none

! sub structure for actual subdomain
      type(subdomain_type),intent(in) :: sub

      ! coarse matrix length
      integer,intent(in)  :: lcoarsem
      ! coarse matrix 
      real(kr),intent(out) :: coarsem(lcoarsem)

      ! local vars
      character(*),parameter:: routine_name = 'DD_GET_COARSEM'
      integer :: i

      ! check the prerequisities
      if (.not.sub%is_coarse_prepared) then
         call error(routine_name, 'Coarse matrix not ready:', sub%isub)
      end if

      ! check dimensions
      if (lcoarsem .ne. sub%lcoarsem) then
         call error(routine_name, 'coarse matrix dimension mismatch:', sub%isub)
      end if

      ! if checks are OK, copy data
      do i = 1,lcoarsem
         coarsem(i) = sub%coarsem(i)
      end do

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
               write(*,*) 'lcoarsem:', lcoarsem, 'lindrowc', lindrowc
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
      
!*****************************************************************
subroutine dd_construct_interior_residual(sub,sol,lsol,reso,lreso)
!*****************************************************************
! Subroutine for construction of interior residual
      use module_utils
      use module_sm
      implicit none

! Subdomain structure
      type(subdomain_type),intent(inout) :: sub

      ! input vector
      integer,intent(in)  :: lsol
      real(kr),intent(in) ::  sol(lsol)

      ! output vector
      integer,intent(in)   :: lreso
      real(kr),intent(out) ::  reso(lreso)

      ! local vars
      character(*),parameter:: routine_name = 'DD_CONSTRUCT_INTERIOR_RESIDUAL'
      integer ::              laux
      real(kr),allocatable ::  aux(:)
      integer ::              lsoli
      real(kr),allocatable ::  soli(:)
      integer ::              lsolo
      real(kr),allocatable ::  solo(:)

      integer :: ndofi, ndofo, ndof, nnza12, la12, nnza11, la11, &
                 matrixtype_aux, matrixtype
      integer :: i
      logical :: remove_original

      ! check the prerequisities
      if (.not. (sub%is_blocked)) then
         remove_original = .false.
         call dd_matrix_tri2blocktri(sub,remove_original)
      end if
      
      ndof  = sub%ndof
      ndofo = sub%ndofo

      ! check dimensions
      if (lsol.ne.ndof) then
         call error(routine_name,'Vector SOL does not have correct dimension.',sub%isub)
      end if
      if (lreso.ne.ndofo) then
         call error(routine_name,'Vector RESO does not have correct dimension.',sub%isub)
      end if
 
      ! prepare rhs vector for backsubstitution to problem A_11*aux1 = -A_12*x
      if (ndofo.eq.0) then
         return
      end if

      ndofi = sub%ndofi

      ! prepare interface part of solution
      lsoli = ndofi
      allocate(soli(lsoli))
      call dd_map_sub_to_subi(sub, sol,lsol, soli,lsoli)

      ! prepare interiOr part of solution
      lsolo = ndofo
      allocate(solo(lsolo))
      call dd_map_sub_to_subo(sub, sol,lsol, solo,lsolo)
   
      ! get reso = A_11*solo
      matrixtype = sub%matrixtype
      nnza11     = sub%nnza11
      la11       = sub%la11
      call sm_vec_mult(matrixtype, nnza11, &
                       sub%i_a11_sparse, sub%j_a11_sparse, sub%a11_sparse, la11, &
                       solo,lsolo, reso,lreso)

      laux = ndofo
      allocate(aux(laux))
      ! prepare rhs vector 
      ! with offdiagonal blocks, use as nonsymmetric
      ! aux1 = A_12 * soli  
      matrixtype_aux = 0
      nnza12     = sub%nnza12
      la12       = sub%la12
      call sm_vec_mult(matrixtype_aux, nnza12, &
                       sub%i_a12_sparse, sub%j_a12_sparse, sub%a12_sparse, la12, &
                       soli,lsoli, aux,laux)

      ! add results together to get reso = reso + aux, i.e. reso = A_11 * solo + A_12 * soli
      do i = 1,lreso
         reso(i) = reso(i) + aux(i)
      end do
 
      deallocate(aux)
      deallocate(solo)
      deallocate(soli)

end subroutine

!**************************************************
subroutine dd_multiply_by_schur(sub,x,lx,y,ly,ncol)
!**************************************************
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

      ! number of columns in x and y
      integer,intent(in)   :: ncol

      ! local vars
      character(*),parameter:: routine_name = 'DD_MULTIPLY_BY_SCHUR'
      integer ::              laux1
      real(kr),allocatable ::  aux1(:)
      integer ::              laux2
      real(kr),allocatable ::  aux2(:)

      integer :: ndofi, ndofo, nnza12, la12, nnza21, la21, nnza22, la22, &
                 matrixtype_aux, matrixtype
      integer :: j
      logical :: is_symmetric_storage

      ! check the prerequisities
      if (.not. (sub%is_explicit_schur_prepared) .and. .not. (sub%is_interior_factorized)) then
         call error( routine_name, 'Interior block not factorized yet.', sub%isub)
      end if
      if (.not. (sub%ndofi .eq. lx/ncol .or. .not. lx .eq. ly)) then
         call error( routine_name, 'Inconsistent data size.', sub%isub)
      end if
 
      ! use explicit Schur complement if it is prepared
      if (sub%is_explicit_schur_prepared) then
         ! copy x to y
         y = x
         if      (sub%istorage == 2) then
              call dsymv('U', sub%lschur1, 1._kr, sub%schur, sub%lschur1, x, 1, 0._kr, y, 1)
         else if (sub%istorage == 1) then
              call dgemv('N', sub%lschur1, sub%lschur2, 1.0_kr, sub%schur, sub%lschur1, x, 1, 0._kr, y, 1)
         else
            call error( routine_name, 'Illegal storage type.', sub%isub)
         end if
      else
         ! prepare rhs vector for backsubstitution to problem A_11*aux1 = -A_12*x
         ndofi = sub%ndofi
         ndofo = sub%ndofo
         if (ndofo.gt.0) then
            laux1 = ndofo * ncol
            allocate(aux1(laux1))
   
            ! prepare rhs vector for backsubstitution to problem A_11*aux1 = -A_12*x
            ! with offdiagonal blocks, use as nonsymmetric
            matrixtype_aux = 0
            nnza12     = sub%nnza12
            la12       = sub%la12

            do j = 1,ncol
               if (ndofi > 0) then
                  call sm_vec_mult(matrixtype_aux, nnza12, &
                                   sub%i_a12_sparse, sub%j_a12_sparse, sub%a12_sparse, la12, &
                                   x((j-1)*ndofi + 1),ndofi, aux1((j-1)*ndofo + 1),ndofo)
               end if
            end do
   
            ! resolve interior problem by MUMPS
            call mumps_resolve(sub%mumps_interior_block,aux1,laux1,ncol)
   
            if (sub%istorage .eq. 4) then
               is_symmetric_storage = .true.
            else
               is_symmetric_storage = .false.
            end if
   
            ! prepare auxiliary vector for multiplication
            ndofi = sub%ndofi
            laux2 = ndofi * ncol
            allocate(aux2(laux2))
   
            ! get aux2 = A_21*aux1, i.e. aux2 = A_21 * (A_11)^-1 * A_12 * x
            do j = 1,ncol
               if (is_symmetric_storage) then
                  matrixtype_aux = 0
                  nnza12     = sub%nnza12
                  la12       = sub%la12
                  ! use the matrix with transposed indices in the call sm_vec_mult
                  if (laux2 > 0) then
                     call sm_vec_mult(matrixtype_aux, nnza12, &
                                      sub%j_a12_sparse, sub%i_a12_sparse, sub%a12_sparse, la12, &
                                      aux1((j-1)*ndofo + 1),ndofo, aux2((j-1)*ndofi + 1),ndofi)
                  end if
               else
                  matrixtype_aux = 0
                  nnza21     = sub%nnza21
                  la21       = sub%la21
                  if (laux2 > 0) then
                     call sm_vec_mult(matrixtype_aux, nnza21, &
                                      sub%i_a21_sparse, sub%j_a21_sparse, sub%a21_sparse, la21, &
                                      aux1((j-1)*ndofo + 1),ndofo, aux2((j-1)*ndofi + 1),ndofi)
                  end if
               end if
            end do
         end if

         ! get y = A_22*x
         matrixtype = sub%matrixtype
         nnza22     = sub%nnza22
         la22       = sub%la22
         do j = 1,ncol
            if (ndofi > 0) then
               call sm_vec_mult(matrixtype, nnza22, &
                                sub%i_a22_sparse, sub%j_a22_sparse, sub%a22_sparse, la22, &
                                x((j-1)*ndofi + 1),ndofi, y((j-1)*ndofi + 1),ndofi)
            end if
         end do

         ! add results together to get y = y - aux2, i.e. y = A_22 * x - A_21 * (A_11)^-1 * A_12 * x, or y = (A_22 - A_21 * (A_11)^-1 * A_12) * x
         if (ndofo.gt.0) then
            y = y - aux2
 
            deallocate(aux1)
            deallocate(aux2)
         end if
      end if

end subroutine

!!*********************************************
!subroutine dd_multiply_by_schur(sub,x,lx,y,ly)
!!*********************************************
!! Subroutine for multiplication of interface vector by Schur complement
!      use module_utils
!      use module_mumps
!      use module_sm
!      implicit none
!
!! Subdomain structure
!      type(subdomain_type),intent(inout) :: sub
!
!      ! input vector
!      integer,intent(in)  :: lx
!      real(kr),intent(in) ::  x(lx)
!
!      ! output vector
!      integer,intent(in)   :: ly
!      real(kr),intent(out) ::  y(ly)
!
!      ! local vars
!      integer ::              laux1
!      real(kr),allocatable ::  aux1(:)
!      integer ::              laux2
!      real(kr),allocatable ::  aux2(:)
!
!      integer :: ndofi, ndofo, nnza12, la12, nnza21, la21, nnza22, la22, &
!                 matrixtype_aux, matrixtype
!      integer :: i
!      logical :: is_symmetric_storage
!
!      ! check the prerequisities
!      if (.not. (sub%is_interior_factorized)) then
!         write(*,*) 'DD_PREPARE_SCHUR: Interior block not factorized yet.',sub%isub
!         call error_exit
!      end if
!      if (.not. (sub%ndofi .eq. lx .or. .not. lx .eq. ly)) then
!         write(*,*) 'DD_PREPARE_SCHUR: Inconsistent data size.'
!         call error_exit
!      end if
! 
!      ! prepare rhs vector for backsubstitution to problem A_11*aux1 = -A_12*x
!      ndofo = sub%ndofo
!      if (ndofo.gt.0) then
!         laux1 = ndofo
!         allocate(aux1(laux1))
!   
!         ! prepare rhs vector for backsubstitution to problem A_11*aux1 = -A_12*x
!         ! with offdiagonal blocks, use as nonsymmetric
!         matrixtype_aux = 0
!         nnza12     = sub%nnza12
!         la12       = sub%la12
!         call sm_vec_mult(matrixtype_aux, nnza12, &
!                          sub%i_a12_sparse, sub%j_a12_sparse, sub%a12_sparse, la12, &
!                          x,lx, aux1,laux1)
!   
!         ! resolve interior problem by MUMPS
!         call mumps_resolve(sub%mumps_interior_block,aux1,laux1)
!   
!         if (sub%istorage .eq. 4) then
!            is_symmetric_storage = .true.
!         else
!            is_symmetric_storage = .false.
!         end if
!   
!         ! prepare auxiliary vector for multiplication
!         ndofi = sub%ndofi
!         laux2 = ndofi
!         allocate(aux2(laux2))
!   
!         ! get aux2 = A_21*aux1, i.e. aux2 = A_21 * (A_11)^-1 * A_12 * x
!         if (is_symmetric_storage) then
!            matrixtype_aux = 0
!            nnza12     = sub%nnza12
!            la12       = sub%la12
!            ! use the matrix with transposed indices in the call sm_vec_mult
!            call sm_vec_mult(matrixtype_aux, nnza12, &
!                             sub%j_a12_sparse, sub%i_a12_sparse, sub%a12_sparse, la12, &
!                             aux1,laux1, aux2,laux2)
!         else
!            matrixtype_aux = 0
!            nnza21     = sub%nnza21
!            la21       = sub%la21
!            call sm_vec_mult(matrixtype_aux, nnza21, &
!                             sub%i_a21_sparse, sub%j_a21_sparse, sub%a21_sparse, la21, &
!                             aux1,laux1, aux2,laux2)
!         end if
!      end if
!
!      ! get y = A_22*x
!      matrixtype = sub%matrixtype
!      nnza22     = sub%nnza22
!      la22       = sub%la22
!      call sm_vec_mult(matrixtype, nnza22, &
!                       sub%i_a22_sparse, sub%j_a22_sparse, sub%a22_sparse, la22, &
!                       x,lx, y,ly)
!
!      ! add results together to get y = y - aux2, i.e. y = A_22 * x - A_21 * (A_11)^-1 * A_12 * x, or y = (A_22 - A_21 * (A_11)^-1 * A_12) * x
!      if (ndofo.gt.0) then
!         do i = 1,ly
!            y(i) = y(i) - aux2(i)
!         end do
! 
!         deallocate(aux1)
!         deallocate(aux2)
!      end if
!
!end subroutine

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
      integer :: nrhs

      ! check the prerequisities
      if (.not. sub%is_explicit_schur_prepared .and. .not. sub%is_blocked) then
         write(*,*) 'DD_RESOLVE_INTERIOR: Matrix is not in blocked format. Call routine to do this.'
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

      if (sub%is_explicit_schur_prepared) then
         ! get the expanded solution from MUMPS

         nrhs = 1
         call mumps_get_expanded_solution(sub%mumps_subdomain_matrix, nrhs, x,lx, y,ly)
      else
         if (.not. (sub%is_blocked)) then
            write(*,*) 'DD_RESOLVE_INTERIOR: Matrix is not in blocked format. Call routine to do this.'
            call error_exit
         end if
         if (.not. (sub%is_interior_factorized)) then
            write(*,*) 'DD_RESOLVE_INTERIOR: Interior block not factorized yet.'
            call error_exit
         end if
         ndofo = sub%ndofo
         ! continue only for nontrivial interior block
         if (ndofo.gt.0) then
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
            call dd_solve_interior_problem(sub,aux1,laux1)

            ! embed vector of interior solution into subdomain solution y
            do io = 1,ndofo
               ind = sub%iovsvn(io)
               y(ind) = y(ind) + aux1(io)
            end do

            deallocate(aux1)
         end if

         ! embed interface solution x into y
         call dd_map_subi_to_sub(sub, x,lx, y,ly)
      end if
      
end subroutine

!*************************************************
subroutine dd_solve_interior_problem(sub,vec,lvec)
!*************************************************
! Subroutine for resolution of interior variables with prepared right-hand side
      use module_mumps
      use module_utils
      implicit none

! Subdomain structure
      type(subdomain_type),intent(inout) :: sub

      ! input vector
      integer,intent(in)     :: lvec
      real(kr),intent(inout) ::  vec(lvec)

      ! local vars
      character(*),parameter:: routine_name = 'DD_SOLVE_INTERIOR_PROBLEM'
      integer :: ndofo
      integer :: ndof
      real(kr),allocatable ::  vecs(:)
      integer :: lvecs
      integer :: nrhs

      ! check dimensions
      ndofo = sub%ndofo
      if (lvec .ne. ndofo) then
         call error(routine_name,'Data size mismatch for subdomain:',sub%isub)
      end if

      ! continue only for nontrivial interior block
      if (ndofo.eq.0) then
         return
      end if

      if (sub%is_explicit_schur_prepared) then
         ndof = sub%ndof
         lvecs = ndof
         allocate(vecs(lvecs))
         vecs = 0._kr
         !expand solution
         call dd_map_subo_to_sub(sub, vec,lvec, vecs,lvecs)

         nrhs = 1
         call mumps_get_interior_solution(sub%mumps_subdomain_matrix, nrhs, vecs,lvecs)

         vec = 0._kr
         call dd_map_sub_to_subo(sub, vecs,lvecs, vec,lvec)
         deallocate(vecs)
      else
         ! check the prerequisities
         if (.not. (sub%is_interior_factorized)) then
            call error(routine_name,'Interior block not factorized yet for subdomain:',sub%isub)
         end if

         ! resolve interior problem by MUMPS
         call mumps_resolve(sub%mumps_interior_block,vec,lvec)
      end if

end subroutine

!***********************************************************************************************
subroutine dd_prepare_reduced_rhs_all(suba,lsuba,sub2proc,lsub2proc,indexsub,lindexsub,comm_all)
!***********************************************************************************************
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
      character(*),parameter:: routine_name = 'DD_PREPARE_REDUCED_RHS_ALL'
      integer ::              lrhs
      real(kr),allocatable ::  rhs(:)
      integer ::              lbc
      real(kr),allocatable ::  bc(:)
      integer ::              lsolo
      real(kr),allocatable ::  solo(:)
      integer ::              lg
      real(kr),allocatable ::  g(:)

      integer :: ndofi, ndofo, ndof
      integer :: i, isub_loc, isub

      ! loop over subdomains
      do isub_loc = 1,lindexsub
         isub = indexsub(isub_loc)

         ! check the prerequisities
         if (.not. suba(isub_loc)%is_explicit_schur_prepared .and. &
             .not. suba(isub_loc)%is_interior_factorized    ) then
            call error(routine_name,'Interior block not factorized yet.')
         end if
         if (.not. (suba(isub_loc)%is_rhs_loaded)) then
            call error(routine_name,'RHS not loaded.')
         end if
         if (.not. (suba(isub_loc)%is_weights_ready)) then
            call error(routine_name,'Weights not ready.')
         end if
 
         ndof = suba(isub_loc)%ndof
         ndofo = suba(isub_loc)%ndofo
         ndofi = suba(isub_loc)%ndofi

         ! prepare subdomain rhs array
         lrhs = ndof
         allocate(rhs(lrhs))

         ! prepare interior solution array
         lsolo = ndofo
         allocate(solo(lsolo))

         ! prepare reduced rhs array
         lg = ndofi
         allocate(g(lg))

         ! copy right hand side into aux2
         do i = 1,ndof
            rhs(i) = suba(isub_loc)%rhs(i)
         end do

         ! only if restriction of RHS is loaded, i.e. not assembled subdomain RHS,
         if ( suba(isub_loc)%is_rhs_complete ) then
            ! apply weights on interface
            call dd_weights_apply(suba(isub_loc), rhs,lrhs)
         end if

         ! fix BC in aux2
         if (suba(isub_loc)%is_bc_present) then
            ! prepare subdomain BC array
            lbc = ndof
            allocate(bc(lbc))
            ! copy BC
            do i = 1,ndof
               bc(i) = suba(isub_loc)%bc(i)
            end do

            ! apply weights on interface
            call dd_weights_apply(suba(isub_loc), bc,lbc)
            call sm_prepare_rhs(suba(isub_loc)%ifix,suba(isub_loc)%lifix,&
                                bc,lbc,rhs,lrhs)
            deallocate(bc)
         end if

         ! prepare interior portion of solution
         call dd_prepare_reduced_rhs(suba(isub_loc),rhs,lrhs,solo,lsolo,g,lg)

         ! store interior solution in the structure
         call dd_upload_interior_solution(suba(isub_loc),solo,lsolo)

         ! load reduced RHS for communication structure
         call dd_comm_upload(suba(isub_loc), g,lg) 
 
         ! save my own g into sub
         ! is it already allocated?
         if (.not.allocated(suba(isub_loc)%g)) then
            ! if not, allocate it
            suba(isub_loc)%lg = lg
            allocate(suba(isub_loc)%g(lg))
         else
            ! if yes, check that it is allocated to correct dimension
            if (lg.ne.suba(isub_loc)%lg) then
               call error(routine_name,'Dimension of subdomain G mismatch for subdomain:',isub)
            end if
         end if
         do i = 1,lg
            suba(isub_loc)%g(i) = g(i)
         end do

         deallocate(rhs)
         deallocate(solo)
         deallocate(g)
      end do

      ! communicate condensed right hand side
      call dd_comm_swapdata(suba,lsuba, indexsub,lindexsub, sub2proc,lsub2proc,comm_all)

      ! loop over subdomains
      do isub_loc = 1,lindexsub
         isub = indexsub(isub_loc)

         ndofi = suba(isub_loc)%ndofi

         ! download contribution to g from my neighbours
         call dd_comm_download(suba(isub_loc), suba(isub_loc)%g,suba(isub_loc)%lg) 

         suba(isub_loc)%is_reduced_rhs_loaded = .true.
 
      end do
end subroutine

!****************************************************************
subroutine dd_prepare_reduced_rhs(sub,rhs,lrhs, solo,lsolo, g,lg)
!****************************************************************
! Subroutine for construction of interior solution SOLO and reduced rhs gi
! |A_11 A_12| u_1 = rhs_1
! |A_21 A_22| u_2   rhs_2
! solo = (A_11)^-1 * rhs_1
! g   = rhs_2 - A_21 * solo
      use module_mumps
      use module_sm
      use module_utils
      implicit none

! array of sub structure for actual subdomains
      type(subdomain_type),intent(inout) :: sub
! subdomain RHS
      integer,intent(in) ::  lrhs
      real(kr),intent(in) ::  rhs(lrhs)
! subdomain interiOr solution
      integer,intent(in) ::   lsolo
      real(kr),intent(out) ::  solo(lsolo)
! subdomain reduced RHS
      integer,intent(in) ::   lg
      real(kr),intent(out) ::  g(lg)

      ! local vars
      character(*),parameter:: routine_name = 'DD_PREPARE_REDUCED_RHS'
      integer :: ndofi, ndofo, ndof, nnza12, la12, nnza21, la21, &
                 matrixtype_aux
      integer :: i
      logical :: is_symmetric_storage
      integer :: nrhs
      real(kr),allocatable ::  rhs_aux(:)


      ndof  = sub%ndof
      ndofo = sub%ndofo
      ndofi = sub%ndofi

      ! check dimensions 
      if (lrhs.ne.ndof) then
         call error(routine_name,'Dimension of RHS array mismatch.')
      end if
      if (lsolo.ne.ndofo) then
         call error(routine_name,'Dimension of SOLO array mismatch.')
      end if
      if (lg.ne.ndofi) then
         call error(routine_name,'Dimension of G array mismatch.')
      end if

      if (sub%is_explicit_schur_prepared) then
         ! for now, do not return the interior solution

         allocate(rhs_aux(lrhs))
         rhs_aux = rhs
         solo = 0._kr
         nrhs = 1
         call mumps_get_interior_solution(sub%mumps_subdomain_matrix, nrhs, rhs_aux,lrhs)
         call dd_map_sub_to_subo(sub, rhs_aux,lrhs, solo,lsolo)

         nrhs = 1
         rhs_aux = rhs
         call mumps_get_reduced_rhs(sub%mumps_subdomain_matrix, nrhs, rhs_aux,lrhs, g,lg)

         deallocate(rhs_aux)
      else
         ! check the prerequisities
         if (.not. (sub%is_interior_factorized)) then
            call error(routine_name,'Interior block not factorized yet.')
         end if

         ! initialize solo and g
         solo(:) = 0._kr
         g(:)    = 0._kr

         ! continue for nontrivial interior block
         if (ndofo.gt.0) then

            ! prepare f_1
            ! map RHS (length of whole subdomain solution) to subdomain interior SOLO
            call dd_map_sub_to_subo(sub, rhs,lrhs, solo,lsolo)
            ! SOLO now contains interior RHS

            ! solve problem A_11*solo = f_1
            ! by MUMPS
            call mumps_resolve(sub%mumps_interior_block,solo,lsolo)
            ! SOLO now contains interior solution
            
            if (sub%istorage .eq. 4) then
               is_symmetric_storage = .true.
            else
               is_symmetric_storage = .false.
            end if

            ! get g = A_21*sol1, i.e. g = A_21 * (A_11)^-1 * rhs_1
            if (is_symmetric_storage) then
               matrixtype_aux = 0
               nnza12     = sub%nnza12
               la12       = sub%la12
               ! use the matrix with transposed indices in the call sm_vec_mult
               call sm_vec_mult(matrixtype_aux, nnza12, &
                                sub%j_a12_sparse, sub%i_a12_sparse, sub%a12_sparse, la12,&
                                solo,lsolo, g,lg)
            else
               matrixtype_aux = 0
               nnza21     = sub%nnza21
               la21       = sub%la21
               call sm_vec_mult(matrixtype_aux, nnza21, &
                                sub%i_a21_sparse, sub%j_a21_sparse, sub%a21_sparse, la21, &
                                solo,lsolo, g,lg)
            end if
         end if

         ! add rhs_2
         ! copy proper part of RHS
         do i = 1,ndofi
            g(i) = rhs(sub%iivsvn(i)) - g(i) 
         end do
      end if

end subroutine

!***************************************
subroutine dd_get_reduced_rhs(sub, g,lg)
!***************************************
! Subroutine for getting reduced rhs
      use module_utils
      implicit none

! Subdomain structure
      type(subdomain_type),intent(in) :: sub
      ! reduced rhs
      integer,intent(in) ::  lg
      real(kr),intent(out) :: g(lg)

      ! local vars
      character(*),parameter:: routine_name = 'DD_GET_REDUCED_RHS'
      integer :: i


      ! check the prerequisities
      if (.not.sub%is_reduced_rhs_loaded) then
         call error(routine_name,'Reduced RHS is not loaded, perhaps missing call to BDDCML_SETUP_NEW_DATA', sub%isub)
      end if
      ! check the size
      if (lg .ne. sub%lg) then
         call error(routine_name,'RHS size mismatch for subdomain:', sub%isub)
      end if
 
      ! copy g
      do i = 1,lg
         g(i) = sub%g(i)
      end do

end subroutine

!***************************************
subroutine dd_fix_reduced_rhs(sub, g,lg)
!***************************************
! Subroutine for fixing boundary conditions in reduced rhs
      use module_sm
      use module_utils
      implicit none

! Subdomain structure
      type(subdomain_type),intent(in) :: sub
      ! reduced rhs
      integer,intent(in) ::    lg
      real(kr),intent(inout) :: g(lg)

      ! local vars
      character(*),parameter:: routine_name = 'DD_FIX_REDUCED_RHS'

      integer :: ndof
      integer ::              lrhs
      real(kr), allocatable :: rhs(:)


      ! check the prerequisities
      if (.not.sub%is_bc_loaded) then
         call error(routine_name,'BC is not loaded for subdomain:', sub%isub)
      end if
      ! check the size
      if (lg .ne. sub%ndofi) then
         call error(routine_name,'RHS size mismatch for subdomain:', sub%isub)
      end if
 
      if (sub%is_bc_present) then
         ! fix nonhomogenous RHS
         ndof = sub%ndof

         lrhs = ndof
         allocate(rhs(lrhs))
         rhs(:) = 0._kr

         call dd_map_subi_to_sub(sub, g,lg, rhs,lrhs)

         call sm_prepare_rhs(sub%ifix,sub%lifix,sub%bc,sub%lbc,rhs,lrhs)

         call dd_map_sub_to_subi(sub, rhs,lrhs, g,lg)
         deallocate(rhs)
      end if

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

!***********************************************
subroutine dd_is_mesh_loaded(sub,is_mesh_loaded)
!***********************************************
! Subroutine for querying if mesh is loaded for the subdomain
      use module_utils
      implicit none
! Subdomain structure
      type(subdomain_type),intent(in) :: sub
! Is mesh loaded?
      logical,intent(out) :: is_mesh_loaded

      is_mesh_loaded = sub%is_mesh_loaded
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

!*************************************************
subroutine dd_get_mesh_basic_size(sub,linet,lnnet)
!*************************************************
! Subroutine for finding size of subdomain mesh
      use module_utils
      implicit none
! Subdomain structure
      type(subdomain_type),intent(in) :: sub
! Length of subdomain INET array
      integer,intent(out) :: linet
! Length of subdomain NNET array
      integer,intent(out) :: lnnet

! local vars
      character(*),parameter:: routine_name = 'DD_GET_MESH_BASIC_SIZE'


      if (.not.sub%is_mesh_loaded) then
         call error(routine_name,'Mesh not loaded for subdomain ',sub%isub)
      end if

      ! if all checks are OK, return data
      linet = sub%linet
      lnnet = sub%lnnet

end subroutine

!*************************************************************************
subroutine dd_get_mesh_basic(sub,use_global_indices,inet,linet,nnet,lnnet)
!*************************************************************************
! Subroutine for finding size of subdomain mesh
      use module_utils
      implicit none
! Subdomain structure
      type(subdomain_type),intent(in) :: sub
! Should INET array contain global indices rather than local to subdomain (native to storage)
      logical,intent(in) :: use_global_indices
! Length of subdomain INET array
      integer,intent(in) :: linet
! Local INET array
      integer,intent(out) :: inet(linet)
! Length of subdomain NNET array
      integer,intent(in) :: lnnet
! Local NNET array
      integer,intent(out) :: nnet(lnnet)

! local vars
      character(*),parameter:: routine_name = 'DD_GET_MESH_BASIC'
      integer :: i, indnode, indnodeg

      ! check prerequisites
      if (.not.sub%is_mesh_loaded) then
         call error(routine_name,'Mesh not loaded for subdomain ',sub%isub)
      end if

      ! check dimensions
      if (linet.ne.sub%linet) then
         call error(routine_name,'Dimension of INET mismatch for subdomain ',sub%isub)
      end if
      if (lnnet.ne.sub%lnnet) then
         call error(routine_name,'Dimension of NNET mismatch for subdomain ',sub%isub)
      end if

      ! if all checks are OK, copy INET and NNET arrays
      do i = 1,linet
         inet(i) = sub%inet(i)
      end do
      do i = 1,lnnet
         nnet(i) = sub%nnet(i)
      end do

      ! embed the INET array into global numbering if required
      if (use_global_indices) then
         do i = 1,linet
            indnode  = inet(i)
            indnodeg = sub%isngn(indnode)

            inet(i) = indnodeg
         end do
      end if
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

!***********************************************************
subroutine dd_get_coarse_coordinates(sub,xyzc,lxyzc1,lxyzc2)
!***********************************************************
! Subroutine for getting coordinates of subdomain coarse nodes
      use module_utils
      implicit none
! Subdomain structure
      type(subdomain_type),intent(in) :: sub
! Coordinates of subdomain coarse nodes
      integer,intent(in)  :: lxyzc1
      integer,intent(in)  :: lxyzc2
      real(kr),intent(out) :: xyzc(lxyzc1,lxyzc2)

      ! local vars
      character(*),parameter:: routine_name = 'DD_GET_COARSE_COORDINATES'
      integer :: ncnodes, icn, ndim, idm 

      ! check prerequisites
      if (.not.sub%is_cnodes_loaded) then
         call error(routine_name, 'Coarse nodes not loaded for subdomain',sub%isub)
      end if

      ncnodes = sub%ncnodes
      ndim    = sub%ndim

      ! check length of output array
      if (ncnodes .ne. lxyzc1) then
         call error(routine_name, 'First dimension of array for coordinates mismatch for subdomain',sub%isub)
      end if
      if (ndim .ne. lxyzc2) then
         call error(routine_name, 'Second dimension of array for coordinates mismatch for subdomain',sub%isub)
      end if

      ! if all checks are OK, return subdomain coarse nodes coordinates
      do icn = 1,ncnodes
         do idm = 1,ndim
            xyzc(icn,idm) = sub%cnodes(icn)%xyz(idm)
         end do
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

!*****************************************************
subroutine dd_get_subdomain_isvgvn(sub,isvgvn,lisvgvn)
!*****************************************************
! Subroutine for getting subdomain Indices of Subdomain Variables in Global Variable Numbering (ISVGVN) array
      use module_utils
      implicit none
! Subdomain structure
      type(subdomain_type),intent(in) :: sub
! Array ISNGN
      integer,intent(in) :: lisvgvn
      integer,intent(out) :: isvgvn(lisvgvn)
! local vars
      character(*),parameter:: routine_name = 'DD_GET_SUBDOMAIN_ISVGVN'
      integer :: ndof, i

      if (.not.sub%is_mesh_loaded) then
         call error(routine_name,'Mesh is not loaded yet for subdomain ',sub%isub)
      end if
      if (sub%ndof .ne. lisvgvn) then
         call error(routine_name,'Size of array ISVGVN not consistent.',sub%isub)
      end if

      ! if all checks are OK, construct array of indices of subdomain nodes in global numbering
      ndof = sub%ndof
      do i = 1,ndof
         isvgvn(i) = sub%isvgvn(i)
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

!****************************************
subroutine dd_prepare_explicit_schur(sub)
!****************************************
! Subroutine for explicit construction of Schur complement
      use module_utils
      implicit none
! Subdomain structure
      type(subdomain_type),intent(inout) :: sub


      ! local vars
      character(*),parameter:: routine_name = 'DD_PREPARE_EXPLICIT_SCHUR'
      integer :: j, ndofi

      integer  ::  lschur1, lschur2

      integer              :: le
      real(kr),allocatable ::  e(:)

      ! check if interface is loaded
      if (.not. sub%is_interface_loaded) then
         call error(routine_name,'Interface not loaded for subdomain: ',sub%isub)
      end if

      ! get dimensions
      ndofi = sub%ndofi
      lschur1 = ndofi
      lschur2 = ndofi
      allocate(sub%schur(lschur1,lschur2))
      sub%lschur1 = lschur1
      sub%lschur2 = lschur2

      le = ndofi
      allocate(e(le))
      do j = 1,ndofi
         ! construct vector of cartesian basis
         e(:) = 0._kr
         e(j) = 1._kr

         ! get column of S*I
         call dd_multiply_by_schur(sub,e,le,sub%schur(1:lschur1,j),le,1)

      end do
      deallocate(e)

      sub%is_explicit_schur_prepared = .true.

end subroutine

!*************************************************
subroutine dd_get_schur(sub,schur,lschur1,lschur2)
!*************************************************
! Subroutine for getting Schur complement from the structure
      use module_utils
      implicit none
! Subdomain structure
      type(subdomain_type),intent(inout) :: sub

      integer,intent(in)  ::  lschur1, lschur2
      real(kr),intent(out) ::  schur(lschur1,lschur2)

      ! local vars
      character(*),parameter:: routine_name = 'DD_GET_SCHUR'
      integer :: i, j

      ! check if explicit Schur is loaded
      if (.not. sub%is_explicit_schur_prepared) then
         call error(routine_name,'Schur complement not ready for subdomain: ',sub%isub)
      end if
      ! check dimensions
      if (lschur1 .ne. sub%lschur1) then
         call error(routine_name,'Schur complement first dimension mismatch for subdomain: ',sub%isub)
      end if
      if (lschur2 .ne. sub%lschur2) then
         call error(routine_name,'Schur complement second dimension mismatch for subdomain: ',sub%isub)
      end if

      ! if all checks are OK, copy Schur complement

      do j = 1,lschur2
         do i = 1,lschur1
            schur(i,j) = sub%schur(i,j)
         end do
      end do

end subroutine

!****************************************
subroutine dd_destroy_explicit_schur(sub)
!****************************************
! Subroutine for deallocating Schur complement from the structure
      use module_utils
      implicit none
! Subdomain structure
      type(subdomain_type),intent(inout) :: sub

      ! local vars
      character(*),parameter:: routine_name = 'DD_DESTROY_EXPLICIT_SCHUR'

      ! check if explicit Schur is loaded
      if (.not. sub%is_explicit_schur_prepared) then
         call error(routine_name,'Schur complement not ready for subdomain: ',sub%isub)
      end if
      if (.not. allocated(sub%schur)) then
         call error(routine_name,'Schur complement seems not allocated for subdomain: ',sub%isub)
      end if

      deallocate(sub%schur)
      sub%lschur1 = 0
      sub%lschur2 = 0

      sub%is_explicit_schur_prepared = .false.

end subroutine

!*********************************************************************************************
subroutine dd_create_neighbouring(suba,lsuba, sub2proc,lsub2proc,indexsub,lindexsub, comm_all)
!*********************************************************************************************
! Subroutine for construction of array of interface and shared nodes with other subdomains
! based on subdomain data
      use module_utils
      use module_pp, only : pp_get_proc_for_sub, pp_get_unique_tag
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
      character(*),parameter:: routine_name = 'DD_CREATE_NEIGHBOURING'
      integer :: isub_loc, isub, isubadj
      integer :: i, ia, inodis, inods, indg
      integer :: kishnadj
      integer :: nsub
      integer :: nadj
      integer :: nnod
      integer :: nshared
      integer :: nsharedv, lcommvec, ndofn
      integer :: idofi, idofo, idofn, inod, inodi, &
                 kisngnadj, ndofi, ndofo, nnodi, nnodadj, &
                 indcomponentadj, inodsadj

      integer ::             lkadjsub
      integer,allocatable ::  kadjsub(:)

      integer ::             iadjcorr

      integer ::             lnshnadj
      integer,allocatable ::  nshnadj(:)
      integer ::             lkshvadj
      integer,allocatable ::  kshvadj(:)
      integer ::             lishnadj
      integer,allocatable ::  ishnadj(:)
      integer ::             lishnncadj
      integer,allocatable ::  ishnncadj(:)
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

      integer :: lsub_aux
      type(sub_aux_type),allocatable :: sub_aux(:)

      ! MPI related arrays and variables
      integer :: myid, nproc, procadj
      integer :: isubadj_loc
      integer :: ierr
      integer :: tag
      integer :: ireq, nreq
      integer ::            lrequest
      integer,allocatable :: request(:)
      integer             :: lstatarray1
      integer             :: lstatarray2
      integer,allocatable :: statarray(:,:)

! orient in communicators
      call MPI_COMM_RANK(comm_all,myid,ierr)
      call MPI_COMM_SIZE(comm_all,nproc,ierr)

      nsub = sub2proc(nproc+1) - 1

      lsub_aux = lsuba
      allocate(sub_aux(lsub_aux))


      ! prepare memory for each subdomain
      ireq = 0
      do isub_loc = 1,lindexsub
         isub = indexsub(isub_loc)
   
         if (.not. suba(isub_loc)%is_mesh_loaded) then
            call error(routine_name, 'Mesh not loaded for subdomain: ',isub)
         end if
   
         ! load data
         nadj  = suba(isub_loc)%nadj
         nnod  = suba(isub_loc)%nnod
   
         sub_aux(isub_loc)%lnnodadj = nadj
         allocate(sub_aux(isub_loc)%nnodadj(sub_aux(isub_loc)%lnnodadj))
   
         ireq = ireq + nadj
      end do
      nreq = ireq

      ! MPI arrays
      lrequest = nreq * 4
      allocate(request(lrequest))
      lstatarray1 = MPI_STATUS_SIZE
      lstatarray2 = nreq * 4
      allocate(statarray(lstatarray1,lstatarray2))

      ! prepare subdomain data      
      ireq = 0
      do isub_loc = 1,lindexsub
         isub = indexsub(isub_loc)
   
         ! load data
         nadj  = suba(isub_loc)%nadj
         nnod  = suba(isub_loc)%nnod
   
   ! Determine sizes of interace of my neighbours
         do ia = 1,nadj
            ! get index of neighbour
            isubadj = suba(isub_loc)%iadj(ia)
      
            ! who owns this subdomain?
            call pp_get_proc_for_sub(isubadj,comm_all,sub2proc,lsub2proc,procadj)
            if (procadj.gt.nproc-1.or.procadj.lt.0) then
               call error(routine_name,'out of range of processors')
            end if
      
            if (procadj .ne. myid) then
               ! interchange via MPI
      
               ! make a receiving request for his data
               call pp_get_unique_tag(isubadj,isub,comm_all,sub2proc,lsub2proc,tag)
               ireq = ireq + 1
               call MPI_IRECV(sub_aux(isub_loc)%nnodadj(ia),1,MPI_INTEGER,procadj,tag,comm_all,&
                              request(ireq),ierr)
               !print *, 'myid =',myid,'receiving one integer from',procadj,' tag',tag
            else 
               ! I have subdomain data, simply copy necessary arrays
               call get_index(isubadj,indexsub,lindexsub,isubadj_loc)
               if (.not. suba(isubadj_loc)%is_mesh_loaded) then
                  call error(routine_name, 'Mesh not loaded for subdomain: ',isubadj)
               end if
      
               sub_aux(isub_loc)%nnodadj(ia) = suba(isubadj_loc)%nnod
            end if
         end do
      end do
      ! prepare subdomain data      
      do isub_loc = 1,lindexsub
         isub = indexsub(isub_loc)
   
         ! load data
         nadj  = suba(isub_loc)%nadj
         nnod  = suba(isub_loc)%nnod
   
   ! Determine sizes of interace of my neighbours
         do ia = 1,nadj
            ! get index of neighbour
            isubadj = suba(isub_loc)%iadj(ia)
      
            ! who owns this subdomain?
            call pp_get_proc_for_sub(isubadj,comm_all,sub2proc,lsub2proc,procadj)
            if (procadj.gt.nproc-1.or.procadj.lt.0) then
               call error(routine_name,'out of range of processors')
            end if
      
            if (procadj .ne. myid) then
               ! interchange via MPI
      
               ! send him my data
               call pp_get_unique_tag(isub,isubadj,comm_all,sub2proc,lsub2proc,tag)
               ireq = ireq + 1
               call MPI_ISEND(nnod,1,MPI_INTEGER,procadj,tag,comm_all,request(ireq),ierr)
               !print *, 'myid =',myid,'Sending', nnod,'to ',procadj,' tag',tag
            end if
         end do
      end do
      nreq = ireq
      ! waiting for communication to complete in this round
      if (nreq.gt.0) then
         call MPI_WAITALL(nreq,request,statarray,ierr)
      end if

      ! Prepare memory to interchange subdomain indices

      do isub_loc = 1,lindexsub
         isub = indexsub(isub_loc)
   
         ! load data
         nadj  = suba(isub_loc)%nadj
         nnod  = suba(isub_loc)%nnod
   
   ! Allocate array for global node numbers for neighbours
         sub_aux(isub_loc)%lisngnadj = sum(sub_aux(isub_loc)%nnodadj)
         allocate(sub_aux(isub_loc)%isngnadj(sub_aux(isub_loc)%lisngnadj))
         allocate(sub_aux(isub_loc)%isncadj( sub_aux(isub_loc)%lisngnadj))
         ! allocate data for myself
         sub_aux(isub_loc)%lisngn = nnod
         allocate(sub_aux(isub_loc)%isngn(sub_aux(isub_loc)%lisngn))
         allocate(sub_aux(isub_loc)%isnc( sub_aux(isub_loc)%lisngn))

         ! prepare local ISNGN array
         do i = 1,nnod
            sub_aux(isub_loc)%isngn(i) = suba(isub_loc)%isngn(i)
            sub_aux(isub_loc)%isnc(i)  = suba(isub_loc)%nodal_components(i)
         end do
      end do
!!
!!      ! now we need to exchange the arrays among all neighbours 
!!      do isub_loc = 1,lindexsub
!!         sub_aux(isub_loc)%comm_array_in => sub_aux(isub_loc)%isngn
!!      end do
!!
!!      call dd_interchange_integer_arrays(suba,lsuba, &
!!                                         sub2proc,lsub2proc, indexsub,lindexsub, &
!!                                         comm_all, &
!!                                         sub_aux,lsub_aux)
!!
!!      do isub_loc = 1,lindexsub
!!
!!         sub_aux(isub_loc)%lisngnadj = size(sub_aux(isub_loc)%comm_array_out)
!!         allocate(sub_aux(isub_loc)%isngnadj(sub_aux(isub_loc)%lisngnadj))
!!         sub_aux(isub_loc)%isngnadj = sub_aux(isub_loc)%comm_array_out
!!         deallocate(sub_aux(isub_loc)%comm_array_out)
!!         deallocate(sub_aux(isub_loc)%comm_array_numbers)
!!      end do
!!
!!      do isub_loc = 1,lindexsub
!!         sub_aux(isub_loc)%comm_array_in => sub_aux(isub_loc)%isnc
!!      end do
!!
!!      call dd_interchange_integer_arrays(suba,lsuba, &
!!                                         sub2proc,lsub2proc, indexsub,lindexsub, &
!!                                         comm_all, &
!!                                         sub_aux,lsub_aux)
!!
!!      do isub_loc = 1,lindexsub
!!
!!         allocate(sub_aux(isub_loc)%isncadj(sub_aux(isub_loc)%lisngnadj))
!!         sub_aux(isub_loc)%isncadj = sub_aux(isub_loc)%comm_array_out
!!         deallocate(sub_aux(isub_loc)%comm_array_out)
!!         sub_aux(isub_loc)%lnnodadj = size(sub_aux(isub_loc)%comm_array_numbers)
!!         allocate(sub_aux(isub_loc)%nnodadj(sub_aux(isub_loc)%lnnodadj))
!!         sub_aux(isub_loc)%nnodadj = sub_aux(isub_loc)%comm_array_numbers
!!         deallocate(sub_aux(isub_loc)%comm_array_numbers)
!!      end do

      ireq = 0
      do isub_loc = 1,lindexsub
         isub = indexsub(isub_loc)
         nsub = suba(isub_loc)%nsub
   
         ! load data
         nadj  = suba(isub_loc)%nadj
         nnod  = suba(isub_loc)%nnod
   
   ! Interchange nodes in global numbering
         kisngnadj = 0
         do ia = 1,nadj
            ! get index of neighbour
            isubadj = suba(isub_loc)%iadj(ia)
            
            ! who owns this subdomain?
            call pp_get_proc_for_sub(isubadj,comm_all,sub2proc,lsub2proc,procadj)
   
            nnodadj = sub_aux(isub_loc)%nnodadj(ia)
            if (nnodadj.eq.0) then
               call error(routine_name,'adjacent subdomain with zero shared nodes?',isubadj)
            end if
      
            if (procadj .ne. myid) then
               ! interchange via MPI
   
               ! make a receiving request for his data
               call pp_get_unique_tag(isubadj,isub,comm_all,sub2proc,lsub2proc,tag)
               ireq = ireq + 1
               call MPI_IRECV(sub_aux(isub_loc)%isngnadj(kisngnadj + 1),nnodadj,MPI_INTEGER,&
                              procadj,tag,comm_all,request(ireq),ierr)
               ireq = ireq + 1
               call MPI_IRECV(sub_aux(isub_loc)%isncadj(kisngnadj + 1),nnodadj,MPI_INTEGER,&
                              procadj,tag,comm_all,request(ireq),ierr)
            else 
               ! I have subdomain data, simply copy necessary arrays
               call get_index(isubadj,indexsub,lindexsub,isubadj_loc)
               do i = 1,nnodadj
                  sub_aux(isub_loc)%isngnadj(kisngnadj + i) = suba(isubadj_loc)%isngn(i)
                  sub_aux(isub_loc)%isncadj (kisngnadj + i) = suba(isubadj_loc)%nodal_components(i)
               end do
            end if
            kisngnadj = kisngnadj + nnodadj
         end do
      end do
      do isub_loc = 1,lindexsub
         isub = indexsub(isub_loc)
         nsub = suba(isub_loc)%nsub
   
         ! load data
         nadj  = suba(isub_loc)%nadj
         nnod  = suba(isub_loc)%nnod
   
   ! Interchange nodes in global numbering
         kisngnadj = 0
         do ia = 1,nadj
            ! get index of neighbour
            isubadj = suba(isub_loc)%iadj(ia)
            
            ! who owns this subdomain?
            call pp_get_proc_for_sub(isubadj,comm_all,sub2proc,lsub2proc,procadj)
   
            nnodadj = sub_aux(isub_loc)%nnodadj(ia)
            if (procadj .ne. myid) then
               ! send him my data
               call pp_get_unique_tag(isub,isubadj,comm_all,sub2proc,lsub2proc,tag)
               ireq = ireq + 1
               call MPI_ISEND(sub_aux(isub_loc)%isngn(1),nnod,MPI_INTEGER, procadj,tag,comm_all,request(ireq),ierr)
               ireq = ireq + 1
               call MPI_ISEND(sub_aux(isub_loc)%isnc(1), nnod,MPI_INTEGER, procadj,tag,comm_all,request(ireq),ierr)
            end if
            kisngnadj = kisngnadj + nnodadj
         end do
      end do
      nreq = ireq
      ! waiting for communication to complete
      if (nreq.gt.0) then
         call MPI_WAITALL(nreq, request, statarray, ierr)
      end if

      deallocate(request)
      deallocate(statarray)

      ! Compare data with neighbours to detect interface
      lkadjsub = nsub
      allocate(kadjsub(lkadjsub))
      do isub_loc = 1,lindexsub
         isub = indexsub(isub_loc)

         ! load data
         nadj  = suba(isub_loc)%nadj
         nnod  = suba(isub_loc)%nnod

         lkinodes = nnod
         allocate(kinodes(lkinodes))
         kinodes = 0

         kisngnadj = 0
         do ia = 1,nadj
            ! get index of neighbour
            isubadj = suba(isub_loc)%iadj(ia)
   
            nnodadj = sub_aux(isub_loc)%nnodadj(ia)
   
            lishared = sub_aux(isub_loc)%lisngn + nnodadj
            allocate(ishared(lishared))
            call get_array_intersection(sub_aux(isub_loc)%isngn,sub_aux(isub_loc)%lisngn,&
                                        sub_aux(isub_loc)%isngnadj(kisngnadj + 1),nnodadj,&
                                        ishared,lishared,nshared)
            !if (nshared.gt.0) then
            !   write(*,*) routine_name,': myid =',myid,', It seems as subdomain ',&
            !               isub, ' does not share nodes with neighbour ',isubadj
            !   call error_exit
            !end if
   
            ! mark interface
            nsharedv = 0
            do i = 1,nshared
               indg = ishared(i)
               call get_index(indg,suba(isub_loc)%isngn,nnod,inods)

               kinodes(inods) = kinodes(inods) + 1
            end do
            deallocate(ishared)
   
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
         deallocate(kinodes)
         deallocate(kdof)

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
         nadj  = suba(isub_loc)%nadj

         lnshnadj = nadj
         allocate(nshnadj(lnshnadj))
         nshnadj = 0
         lkshvadj = nadj + 1
         allocate(kshvadj(lkshvadj))
         lishnadj = nadj * nnodi
         allocate(ishnadj(lishnadj))
         lishnncadj = nadj * nnodi
         allocate(ishnncadj(lishnncadj))

         kadjsub = 0
         kisngnadj = 0
         kishnadj  = 0
         lcommvec  = 0
         kshvadj(1) = 1
         iadjcorr = 0
         do ia = 1,nadj
            ! get index of neighbour
            isubadj = suba(isub_loc)%iadj(ia)
   
            nnodadj = sub_aux(isub_loc)%nnodadj(ia)
   
            lishared = sub_aux(isub_loc)%lisngn + nnodadj
            allocate(ishared(lishared))
            call get_array_intersection(sub_aux(isub_loc)%isngn,sub_aux(isub_loc)%lisngn,&
                                        sub_aux(isub_loc)%isngnadj(kisngnadj + 1),nnodadj,&
                                        ishared,lishared,nshared)
            if (nshared.gt.0) then
               !write(*,*) routine_name,': myid =',myid,', It seems as subdomain ',&
               !            isub, ' does not share nodes with neighbour ',isubadj
               !print *,'isub',isub,'nnodadj',sub_aux(isub_loc)%nnodadj
               !print *,'isub',isub,'nadj',suba(isub_loc)%nadj
               !call error_exit
               iadjcorr = iadjcorr + 1

               kadjsub(isubadj) = 1

               nshnadj(iadjcorr) = nshared
   
               ! load shared nodes 
               nsharedv = 0
               do i = 1,nshared
                  indg = ishared(i)
                  call get_index(indg,suba(isub_loc)%isngn,suba(isub_loc)%lisngn,inods)
                  call get_index(inods,suba(isub_loc)%iin,suba(isub_loc)%nnodi,inodis)

                  ndofn = suba(isub_loc)%nndf(inods)
                  nsharedv = nsharedv + ndofn

                  ishnadj(kishnadj + i) = inodis

                  call get_index(indg,sub_aux(isub_loc)%isngnadj(kisngnadj + 1),nnodadj,inodsadj)
                  indcomponentadj = sub_aux(isub_loc)%isncadj(kisngnadj + inodsadj)

                  ishnncadj(kishnadj + i) = indcomponentadj
               end do

               kshvadj(iadjcorr + 1) = kshvadj(iadjcorr) + nsharedv
               lcommvec = lcommvec + nsharedv

               kishnadj = kishnadj + nshared
            end if

            deallocate(ishared)

            kisngnadj  = kisngnadj  + nnodadj
         end do
   
         ! load info about shared nodes into structure
         suba(isub_loc)%lnshnadj = iadjcorr
         allocate(suba(isub_loc)%nshnadj(lnshnadj))
         do i = 1,suba(isub_loc)%lnshnadj
            suba(isub_loc)%nshnadj(i) = nshnadj(i)
         end do
   
         ! load info about number of shared dof into structure
         suba(isub_loc)%lkshvadj = iadjcorr+1
         allocate(suba(isub_loc)%kshvadj(lkshvadj))
         do i = 1,suba(isub_loc)%lkshvadj
            suba(isub_loc)%kshvadj(i) = kshvadj(i)
         end do

         ! truncate array ishnadj to really used indices
         suba(isub_loc)%lishnadj = sum(nshnadj)
         allocate(suba(isub_loc)%ishnadj(suba(isub_loc)%lishnadj))
         do i = 1,suba(isub_loc)%lishnadj
            suba(isub_loc)%ishnadj(i) = ishnadj(i)
         end do

         ! truncate array ishnncadj to really used indices
         suba(isub_loc)%lishnncadj = sum(nshnadj)
         allocate(suba(isub_loc)%ishnncadj(suba(isub_loc)%lishnncadj))
         do i = 1,suba(isub_loc)%lishnncadj
            suba(isub_loc)%ishnncadj(i) = ishnncadj(i)
         end do

         ! prepare communication arrays
         suba(isub_loc)%lcommvec = lcommvec

         ! update adjacency info
         call dd_localize_adj(suba(isub_loc),nsub,kadjsub,lkadjsub)
         
         ! activate the flag
         suba(isub_loc)%is_neighbouring_ready = .true.
   
         ! debug
         !print *,'isub = ',isub,'nshnadj',nshnadj
         !print *,'isub = ',isub,'kshvadj',kshvadj
         !print *,'isub = ',isub,'ishnadj',ishnadj
         !print *,'isub = ',isub,'ishnncadj',ishnncadj
         deallocate(nshnadj)
         deallocate(kshvadj)
         deallocate(ishnadj)
         deallocate(ishnncadj)

         deallocate(sub_aux(isub_loc)%nnodadj)
         deallocate(sub_aux(isub_loc)%isngn)
         deallocate(sub_aux(isub_loc)%isngnadj)
         deallocate(sub_aux(isub_loc)%isnc)
         deallocate(sub_aux(isub_loc)%isncadj)
      end do
      deallocate(sub_aux)

      deallocate(kadjsub)

end subroutine

!********************************************************************************************************
subroutine dd_guess_neighbouring_by_nodes(suba,lsuba, sub2proc,lsub2proc,indexsub,lindexsub, comm_all,&
                                          kadjsub,lkadjsub)
!********************************************************************************************************
! Subroutine for estimating which subdomain is neighbour to which subdomain based on bounds on global node indices

      use module_utils
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
      integer ::            lkadjsub
      integer,intent(out) :: kadjsub(lkadjsub)

      ! local vars
      character(*),parameter:: routine_name = 'DD_ESTIMATE_NEIGHBOURING_BY_NODES'

      ! MPI related arrays and variables
      integer :: myid, nproc
      integer :: ierr

      integer ::            lnode_id_bounds
      integer,allocatable :: node_id_upper_bounds(:)
      integer,allocatable :: node_id_upper_bounds_aux(:)
      integer,allocatable :: node_id_lower_bounds(:)
      integer,allocatable :: node_id_lower_bounds_aux(:)


      integer :: isub, isub_loc, nsub
      integer :: my_lower_bound
      integer :: my_upper_bound

! orient in communicators
      call MPI_COMM_RANK(comm_all,myid,ierr)
      call MPI_COMM_SIZE(comm_all,nproc,ierr)

      nsub = sub2proc(nproc+1) - 1

! allocate arrays for lower and upper bounds
      lnode_id_bounds = nsub
      allocate(node_id_lower_bounds(lnode_id_bounds))
      allocate(node_id_lower_bounds_aux(lnode_id_bounds))
      node_id_lower_bounds_aux = 0
      allocate(node_id_upper_bounds(lnode_id_bounds))
      allocate(node_id_upper_bounds_aux(lnode_id_bounds))
      node_id_upper_bounds_aux = 0

      ! Extract bounds for each local subdomain
      do isub_loc = 1,lindexsub
         isub = indexsub(isub_loc)
         if ( .not. suba(isub_loc)%is_mesh_loaded ) then
            call error(routine_name, 'Mesh data not available for subdomain',isub)
         end if

         node_id_lower_bounds_aux(isub) = minval(suba(isub_loc)%isngn)
         node_id_upper_bounds_aux(isub) = maxval(suba(isub_loc)%isngn)
      end do

! communicate the data
!*****************************************************************MPI
      call MPI_ALLREDUCE(node_id_lower_bounds_aux,node_id_lower_bounds,lnode_id_bounds, &
                         MPI_INTEGER, MPI_SUM, comm_all, ierr) 
      call MPI_ALLREDUCE(node_id_upper_bounds_aux,node_id_upper_bounds,lnode_id_bounds, &
                         MPI_INTEGER, MPI_SUM, comm_all, ierr) 
!*****************************************************************MPI
      deallocate(node_id_lower_bounds_aux)
      deallocate(node_id_upper_bounds_aux)

      ! check data
      if (any(node_id_lower_bounds.eq.0)) then
         call error(routine_name, 'Error in lower bounds of node indices.')
      end if
      if (any(node_id_upper_bounds.eq.0)) then
         call error(routine_name, 'Error in upper bounds of node indices.')
      end if

      ! find possible intersections
      kadjsub = 0
      do isub_loc = 1,lindexsub
         isub = indexsub(isub_loc)

         my_lower_bound = minval(suba(isub_loc)%isngn)
         my_upper_bound = maxval(suba(isub_loc)%isngn)

         where ( .not. ( my_lower_bound.gt.node_id_upper_bounds .or. my_upper_bound.lt.node_id_lower_bounds) ) &
                 kadjsub((isub_loc-1)*nsub+1:isub_loc*nsub) = 1

         ! keep zero on diagonal
         kadjsub((isub_loc-1)*nsub + isub) = 0
      end do

      deallocate(node_id_lower_bounds)
      deallocate(node_id_upper_bounds)

end subroutine

!********************************************************************************************************
subroutine dd_guess_neighbouring_by_bb(ndim, suba,lsuba, sub2proc,lsub2proc,indexsub,lindexsub, comm_all,&
                                       kadjsub,lkadjsub)
!********************************************************************************************************
! Subroutine for estimating which subdomain is neighbour to which subdomain based on geometrical 
! bounding boxes

      use module_utils
      implicit none
      include "mpif.h"

      integer,intent(in) :: ndim

! array of sub structure for actual subdomains
      integer,intent(in) ::                lsuba
      type(subdomain_type),intent(inout) :: suba(lsuba)

      integer,intent(in) :: lsub2proc
      integer,intent(in) ::  sub2proc(lsub2proc)
      integer,intent(in) :: lindexsub
      integer,intent(in) ::  indexsub(lindexsub)
      integer,intent(in) :: comm_all ! MPI communicator
      integer ::            lkadjsub
      integer,intent(out) :: kadjsub(lkadjsub)

      ! local vars
      character(*),parameter:: routine_name = 'DD_ESTIMATE_NEIGHBOURING_BY_BB'

      integer ::            d

      integer ::            lbb_bounds1
      integer ::            lbb_bounds2
      real(kr),allocatable :: bb_upper_bounds(:,:)
      real(kr),allocatable :: bb_upper_bounds_aux(:,:)
      real(kr),allocatable :: bb_lower_bounds(:,:)
      real(kr),allocatable :: bb_lower_bounds_aux(:,:)

      integer :: isub, isub_loc, nsub
      real(kr) :: my_lower_bound
      real(kr) :: my_upper_bound
      ! MPI related arrays and variables
      integer :: myid, nproc
      integer :: ierr


! orient in communicators
      call MPI_COMM_RANK(comm_all,myid,ierr)
      call MPI_COMM_SIZE(comm_all,nproc,ierr)

      nsub = sub2proc(nproc+1) - 1

! allocate arrays for lower and upper bounds
      lbb_bounds1 = nsub
      lbb_bounds2 = ndim
      allocate(bb_lower_bounds(lbb_bounds1,lbb_bounds2))
      allocate(bb_lower_bounds_aux(lbb_bounds1,lbb_bounds2))
      bb_lower_bounds_aux = 0._kr
      allocate(bb_upper_bounds(lbb_bounds1,lbb_bounds2))
      allocate(bb_upper_bounds_aux(lbb_bounds1,lbb_bounds2))
      bb_upper_bounds_aux = 0._kr

      ! Extract bounds for each local subdomain
      do isub_loc = 1,lindexsub
         isub = indexsub(isub_loc)
         if ( .not. suba(isub_loc)%is_mesh_loaded ) then
            call error(routine_name, 'Mesh data not available for subdomain',isub)
         end if

         do d = 1,ndim
            bb_lower_bounds_aux(isub,d) = minval(suba(isub_loc)%xyz(:,d))
            bb_upper_bounds_aux(isub,d) = maxval(suba(isub_loc)%xyz(:,d))
         end do
      end do

! communicate the data
!*****************************************************************MPI
      call MPI_ALLREDUCE(bb_lower_bounds_aux,bb_lower_bounds,lbb_bounds1*lbb_bounds2, &
                         MPI_DOUBLE_PRECISION, MPI_SUM, comm_all, ierr) 
      call MPI_ALLREDUCE(bb_upper_bounds_aux,bb_upper_bounds,lbb_bounds1*lbb_bounds2, &
                         MPI_DOUBLE_PRECISION, MPI_SUM, comm_all, ierr) 
!*****************************************************************MPI
      deallocate(bb_lower_bounds_aux)
      deallocate(bb_upper_bounds_aux)

      if (any( fuzzyLessThan( bb_upper_bounds, bb_lower_bounds ))) then
         call warning(routine_name, &
                      'Suspicious subdomain with zero size of bounding box at some dimension.')
      end if

      ! debug
      !write(*,*) 'bounding box'
      !do isub = 1,nsub
      !   write(*,'(f5.2,a,f5.2,f5.2,a,f5.2,f5.2,a,f5.2)') ( bb_lower_bounds(isub,d),'-',bb_upper_bounds(isub,d), d = 1,ndim)
      !end do

      ! find possible intersections
      kadjsub = 1
      do isub_loc = 1,lindexsub
         isub = indexsub(isub_loc)

         do d = 1,ndim
            my_lower_bound = bb_lower_bounds(isub,d)
            my_upper_bound = bb_upper_bounds(isub,d)

            where ( ( fuzzyLessThan(bb_upper_bounds(:,d), my_lower_bound) .or. &
                      fuzzyLessThan(my_upper_bound,bb_lower_bounds(:,d)) ) ) &
                    kadjsub((isub_loc-1)*nsub+1:isub_loc*nsub) = 0
         end do

         ! keep zero on diagonal
         kadjsub((isub_loc-1)*nsub + isub) = 0
      end do
      if (all(kadjsub.eq.0) .and. lindexsub.gt.0 .and. nproc.gt.1) then
         call warning(routine_name,'This is strange - seems as if none of my subdomains has any neighbour, myid = ',myid)
      end if

      deallocate(bb_lower_bounds)
      deallocate(bb_upper_bounds)

end subroutine

!********************************************************************************************************
subroutine dd_guess_neighbouring(ndim, suba,lsuba, sub2proc,lsub2proc,indexsub,lindexsub, comm_all,&
                                 kadjsub,lkadjsub)
!********************************************************************************************************
! Subroutine for estimating which subdomain is neighbour to which subdomain 

      use module_utils
      implicit none
      include "mpif.h"

      integer,intent(in) :: ndim

! array of sub structure for actual subdomains
      integer,intent(in) ::                lsuba
      type(subdomain_type),intent(inout) :: suba(lsuba)

      integer,intent(in) :: lsub2proc
      integer,intent(in) ::  sub2proc(lsub2proc)
      integer,intent(in) :: lindexsub
      integer,intent(in) ::  indexsub(lindexsub)
      integer,intent(in) :: comm_all ! MPI communicator
      integer ::            lkadjsub
      integer,intent(out) :: kadjsub(lkadjsub)

      ! local vars
      character(*),parameter:: routine_name = 'DD_ESTIMATE_NEIGHBOURING'

      integer,allocatable :: kadjsub_bb(:)
      integer,allocatable :: kadjsub_node_id(:)

! allocate arrays for lower and upper bounds
      allocate(kadjsub_bb(lkadjsub))
      allocate(kadjsub_node_id(lkadjsub))

      ! perform guess based on node ids
      call dd_guess_neighbouring_by_nodes(suba,lsuba, sub2proc,lsub2proc, indexsub,lindexsub, comm_all, &
                                          kadjsub_bb,lkadjsub )

      ! perform guess based on bounding boxes of subdomains
      call dd_guess_neighbouring_by_bb(ndim,suba,lsuba, sub2proc,lsub2proc, indexsub,lindexsub, comm_all, &
                                       kadjsub_node_id,lkadjsub )

      ! final guess is the intersection of both guesses
      kadjsub = 0
      where (kadjsub_bb.eq.1 .and. kadjsub_node_id.eq.1) kadjsub = 1

      deallocate(kadjsub_bb)
      deallocate(kadjsub_node_id)

end subroutine

!*******************************************************************************************************
subroutine dd_create_globs(suba,lsuba, sub2proc,lsub2proc,indexsub,lindexsub, comm_all, remove_bc_nodes,&
                           meshdim, use_corner_constraints, &
                           ncorner, nedge, nface)
!*******************************************************************************************************
! Subroutine for finding corners in BDDC
! based on subdomain data
      use module_utils
      use module_pp
      use module_graph
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
! should nodes on Dirichlet BC be removed?
      logical,intent(in) :: remove_bc_nodes

! what is the dimension of mesh
      integer,intent(in) :: meshdim

      ! additional search for corners
      ! F - only vertices are considered as corners
      ! T - additional face-based selection
      logical,intent(in) :: use_corner_constraints

! basic properties of the coarse problem
      integer,intent(out) :: ncorner
      integer,intent(out) :: nedge
      integer,intent(out) :: nface

      ! local vars
      character(*),parameter:: routine_name = 'DD_CREATE_GLOBS'
! ######## parameters to set
      ! check for disconnected faces
      ! F - face not analysed for continuity
      ! T - search based algorithm run on each component separately
      logical, parameter :: check_face_components = .true.
! ######## parameters to set

      integer :: isub_loc, isub, isubadj
      integer :: ia, i, inodi, indni, indn, inodshaux
      integer :: nadj, nnodi, ndofi, nshn, ishn, indshn, nsubnx
      integer :: kishnadj
      integer :: iglobs, nglobs, jnodi, nsubn, nglobn, ndim, nsub_loc
      integer :: icorners, iedges, ifaces, ncorners, nedges, nfaces
      integer :: inds, indg
      integer :: start, finish, step, maxl
      real(kr):: maxv
      integer :: indedges, indfaces
      real(kr):: x1, y1, z1, x2, y2, z2, xish, yish, zish
      logical :: we_share_a_face
      integer :: iold, indnew

      ! check continuity
      integer ::            linterface_components
      integer,allocatable :: interface_components(:)
      integer ::             ninterface_components
      integer ::            lmy_interface_components
      integer,allocatable :: my_interface_components(:)
      integer ::             nmy_components
      integer ::            ladj_interface_components
      integer,allocatable :: adj_interface_components(:)
      integer ::             nadj_components
      integer ::            lcompind
      integer,allocatable :: compind(:)
      integer ::             incomp, indcomp
      integer ::             interface_componentsize, iinterface_component

      integer :: indshsub

      ! MPI related variables 
      integer :: ierr, nproc, myid

      integer ::             lnsubnode
      integer,pointer ::  nsubnode(:)
      integer ::             lglobsubs1, lglobsubs2, lglobsubs3
      integer,pointer ::  globsubs(:,:,:)
      integer ::             lkglobs
      integer,pointer ::      kglobs(:)
      integer ::             lglobtypes
      integer,pointer ::      globtypes(:)
      integer ::              indcomponentadj
      integer ::              indcomponent

      integer ::             indglb
      integer ::             index_old, index_new, iedg, ifac

      integer ::             liingns
      integer,allocatable ::  iingns(:)
      integer ::             licnsins
      integer,allocatable ::  icnsins(:)
      integer :: indi
      integer :: icnodes
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
      integer ::            lkdofi
      integer,allocatable :: kdofi(:)
      integer ::            lifixi
      integer,allocatable :: ifixi(:)
      integer :: idofn, ndofn, indivs 
      integer :: pifixi

      integer ::            lnsubglobs
      integer,allocatable :: nsubglobs(:)
      integer ::            lglob_subs1, lglob_subs2
      integer,allocatable :: glob_subs(:,:)
      integer ::             iglob

      integer::            ndofi_loc, ndofi_sub, ndofig
      integer::            nnodi_loc, nnodi_sub, nnodig

      integer::             ldist,   larea
      real(kr),allocatable:: dist(:), area(:)
      integer::             lxyzsh1, lxyzsh2 
      real(kr),allocatable:: xyzsh(:,:)
      integer::             lxyzbase
      real(kr),allocatable:: xyzbase(:)
      integer ::             inodcf(4)
      integer ::             indaux(1)
      integer ::             end_bound
      integer ::             iround

      integer ::             nintersection
      integer ::             lintersection
      integer,allocatable ::  intersection(:)

      integer :: lsub_aux
      type(sub_aux_type),allocatable,target :: sub_aux(:)

      integer             :: ncs
      integer             :: lcornermerge
      integer,allocatable ::  cornermerge(:)
      integer             :: indcg, icscn
      integer             :: minindex


      ! orient in communicator
      call MPI_COMM_RANK(comm_all,myid,ierr)
      call MPI_COMM_SIZE(comm_all,nproc,ierr)

      nsub_loc = sub2proc(myid+2) - sub2proc(myid+1) 
      lsub_aux = nsub_loc
      allocate(sub_aux(lsub_aux))

      ! counter of interface unknowns
      ndofi_loc = 0
      nnodi_loc = 0

      ! Compare data with neighbours to detect interface
      do isub_loc = 1,lindexsub
         isub = indexsub(isub_loc)

         ! check prerequisites
         if (.not. suba(isub_loc)%is_neighbouring_ready) then
            call error(routine_name,'Neighbouring not ready for subdomain',isub)
         end if
         if (remove_bc_nodes .and. .not. suba(isub_loc)%is_bc_loaded) then
            call error(routine_name,'Boundary conditions not ready for subdomain',isub)
         end if

         ! load data
         nadj   = suba(isub_loc)%nadj
         nnodi  = suba(isub_loc)%nnodi
         ndim   = suba(isub_loc)%ndim

         lnsubnode = nnodi
         sub_aux(isub_loc)%lnsubnode = lnsubnode
         allocate(sub_aux(isub_loc)%nsubnode(lnsubnode))
         nsubnode => sub_aux(isub_loc)%nsubnode
         nsubnode(:) = 0

         kishnadj  = 0
         do ia = 1,nadj
            ! get index of neighbour
            isubadj = suba(isub_loc)%iadj(ia)

            ! number of nodes shared with this subdomain
            nshn = suba(isub_loc)%nshnadj(ia)

            do ishn = 1,nshn
               ! index of shared node in interface numbering
               indshn = suba(isub_loc)%ishnadj(kishnadj + ishn)

               ! mark the node to array nsubnode
               nsubnode(indshn) = nsubnode(indshn) + 1

            end do
   
            kishnadj = kishnadj + nshn
         end do

         ! check coverage of interface
         if (any(nsubnode.eq.0)) then
            call error(routine_name,'Coverage of interface by globs failed for subdomain',isub)
         end if

         ! find maximal number of subdomains sharing a node
         nsubnx = maxval(nsubnode)
         ! create list of subdomains at interface nodes
         lglobsubs1 = nnodi
         lglobsubs2 = nsubnx
         lglobsubs3 = 3 ! neighbouring subdomain | component in neighbouring subdomain | my own component
         sub_aux(isub_loc)%lglobsubs1 = lglobsubs1
         sub_aux(isub_loc)%lglobsubs2 = lglobsubs2
         sub_aux(isub_loc)%lglobsubs3 = lglobsubs3
         allocate(sub_aux(isub_loc)%globsubs(lglobsubs1,lglobsubs2,lglobsubs3))
         globsubs => sub_aux(isub_loc)%globsubs
         globsubs = 0
         ! use now nsubnode as counters
         nsubnode(:) = 0
         kishnadj  = 0
         do ia = 1,nadj
            ! get index of neighbour
            isubadj = suba(isub_loc)%iadj(ia)

            ! number of nodes shared with this subdomain
            nshn = suba(isub_loc)%nshnadj(ia)

            do ishn = 1,nshn
               ! index of shared node in interface numbering
               indshn = suba(isub_loc)%ishnadj(kishnadj + ishn)

               ! index of component of the adjacent subdomain 
               indcomponentadj = suba(isub_loc)%ishnncadj(kishnadj + ishn)

               ! determine my own component
               ! index of shared node in subdomain numbering
               indshsub = suba(isub_loc)%iin(indshn)

               if (.not. suba(isub_loc)%is_nodal_components_loaded ) then
                  call error( routine_name, 'Cannot determine my own component, isub: ', suba(isub_loc)%isub )
               end if
               indcomponent = suba(isub_loc)%nodal_components(indshsub)

               ! mark the node to array nsubnode
               nsubnode(indshn) = nsubnode(indshn) + 1

               globsubs(indshn,nsubnode(indshn),1) = isubadj
               ! component of adjacent subdomain is stored in third dimension of the array
               globsubs(indshn,nsubnode(indshn),2) = indcomponentadj
               ! component of actual subdomain is stored in third dimension of the array
               globsubs(indshn,nsubnode(indshn),3) = indcomponent
            end do
   
            kishnadj = kishnadj + nshn
         end do
         ! debug
         !do inodi = 1,nnodi
         !   print *,'isub = ',isub,'inodi',inodi,'nsubnode',(globsubs(inodi,i),i = 1,nsubnode(inodi))
         !end do

         ndofi_sub = 0
         nnodi_sub = 0
         do inodi = 1,nnodi
            if ( max( isub, maxval(globsubs(inodi,1:nsubnode(inodi),1)) ) .eq.  isub ) then
               inds = suba(isub_loc)%iin(inodi)
               ndofn = suba(isub_loc)%nndf(inds)

               ndofi_sub = ndofi_sub + ndofn
               nnodi_sub = nnodi_sub + 1
            end if
         end do

         ! local contribution to global number of interface unknowns (without overlaps)
         ndofi_loc = ndofi_loc + ndofi_sub
         nnodi_loc = nnodi_loc + nnodi_sub

         ! Identify topology of globs
         lglobtypes = nnodi
         sub_aux(isub_loc)%lglobtypes  = lglobtypes
         allocate(sub_aux(isub_loc)%globtypes(lglobtypes))
         globtypes => sub_aux(isub_loc)%globtypes
         lkglobs = nnodi
         sub_aux(isub_loc)%lkglobs     = lkglobs   
         allocate(sub_aux(isub_loc)%kglobs(lkglobs))
         kglobs => sub_aux(isub_loc)%kglobs

         ! Initialize KGLOBS
         kglobs = -1
         globtypes(:) = 0

         iglobs   = 0
         icorners = 0
         iedges   = 0
         ifaces   = 0
         do inodi = 1,nnodi
            if (kglobs(inodi) .eq. -1) then
               ! node is not assigned to glob yet
               iglobs = iglobs + 1
               kglobs(inodi) = iglobs
               ! search for nodes shared by the same set of subdomains and same components
               nsubn = nsubnode(inodi)
               nglobn = 1
               do jnodi = inodi+1,nnodi
                  ! check that the node is not yet assigned to glob and it is not a corner
                  if (kglobs(jnodi).eq.-1 .and. globtypes(jnodi).ne.3) then
                     ! check that number of subdomains sharing the node matches
                     if (nsubnode(jnodi).eq.nsubn) then
                        ! check that all subdomain indices and all components are the same
                        if (all(globsubs(jnodi,1:nsubn,:).eq.globsubs(inodi,1:nsubn,:))) then
                           ! the node belongs to the same glob
                           kglobs(jnodi) = iglobs
                           nglobn = nglobn + 1
                           if (nsubn.eq.1) then
                              ! it is a face
                              globtypes(jnodi) = 1
                           else
                              ! it is an edge
                              globtypes(jnodi) = 2
                           end if
                        end if
                     end if
                  end if
               end do
               ! determine type of glob
               if (nsubn.eq.1) then
                  ! it is a face
                  ifaces = ifaces + 1
                  globtypes(inodi) = 1
               else
                  if (nglobn.gt.1) then
                     ! it is an edge
                     iedges = iedges + 1
                     globtypes(inodi) = 2
                  else
                     ! it is a vertex, add it as a corner
                     icorners = icorners + 1
                     globtypes(inodi) = 3
                  end if
               end if
            end if
         end do
         nglobs   = iglobs
         ncorners = icorners
         nfaces   = ifaces
         nedges   = iedges
         !print *,'initial isub = ',isub,'nglobs',nglobs
         !print *,'initial isub = ',isub,'ncorners',ncorners
         !print *,'initial isub = ',isub,'nedges',nedges
         !print *,'initial isub = ',isub,'nfaces',nfaces
         !print *,'isub = ',isub,'kglobs',kglobs
         !print *,'isub = ',isub,'globtypes',globtypes
         !call flush(6)

         ! where nsubnode is 1, there is a face
         where (nsubnode.eq.1) globtypes = 1

         if (use_corner_constraints) then
            call dd_create_globs_select_corners
         end if

         ! number of all subdomain corners
         ncorners = count(globtypes.eq.3)
         sub_aux(isub_loc)%ncorners = ncorners

         ! make local lists of corners
         sub_aux(isub_loc)%licgn = ncorners
         ! Indices of Corners in Global Numbering
         allocate(sub_aux(isub_loc)%icgn(sub_aux(isub_loc)%licgn))

         icorners = 0
         do inodi = 1,nnodi
            nnodi  = suba(isub_loc)%nnodi

            if (globtypes(inodi).eq.3) then
               icorners = icorners + 1

               inds = suba(isub_loc)%iin(inodi)
               indg = suba(isub_loc)%isngn(inds)

               sub_aux(isub_loc)%icgn(icorners) = indg
            end if
         end do

         nullify(globtypes)
         nullify(kglobs)
         nullify(globsubs)
         nullify(nsubnode)
      end do

      ! now we need to exchange the arrays of corners among all neighbours - first exchange their numbers
      do isub_loc = 1,lindexsub
         sub_aux(isub_loc)%comm_array_in => sub_aux(isub_loc)%icgn
      end do

      call dd_interchange_integer_arrays(suba,lsuba, &
                                         sub2proc,lsub2proc, indexsub,lindexsub, &
                                         comm_all, &
                                         sub_aux,lsub_aux)

      do isub_loc = 1,lindexsub

         sub_aux(isub_loc)%lineibcorners = size(sub_aux(isub_loc)%comm_array_out)
         allocate(sub_aux(isub_loc)%ineibcorners(sub_aux(isub_loc)%lineibcorners))
         sub_aux(isub_loc)%ineibcorners = sub_aux(isub_loc)%comm_array_out
         deallocate(sub_aux(isub_loc)%comm_array_out)
         deallocate(sub_aux(isub_loc)%comm_array_numbers)
      end do

      ! process array with corners from neigbours
      do isub_loc = 1,lindexsub
         isub = indexsub(isub_loc)

         nadj  = suba(isub_loc)%nadj
         nnodi = suba(isub_loc)%nnodi

         ! merge array of local corners and corners by others
         lcornermerge = sub_aux(isub_loc)%ncorners + sub_aux(isub_loc)%lineibcorners
         allocate(cornermerge(lcornermerge))

         cornermerge(1:sub_aux(isub_loc)%ncorners)  = sub_aux(isub_loc)%icgn
         cornermerge(sub_aux(isub_loc)%ncorners+1:) = sub_aux(isub_loc)%ineibcorners
         deallocate(sub_aux(isub_loc)%ineibcorners)

         ! sort local array of corners before its send
         call iquick_sort(cornermerge,lcornermerge)
         ! remove repeated indices
         call get_array_norepeat(cornermerge,lcornermerge,ncs)

         ! prepare array for global indices of interface nodes 
         liingns = nnodi
         allocate(iingns(liingns))
         call dd_get_interface_global_numbers(suba(isub_loc), iingns,liingns)
         ! now, it contains all corners from the current and neighbouring subdomains
         ! perform intersection with local nodes

         lintersection = liingns + ncs
         allocate(intersection(lintersection))
         call get_array_intersection(iingns,liingns,cornermerge,ncs,&
                                     intersection,lintersection,nintersection)

         ! update local corners to union of all corners including neighbours
         deallocate(sub_aux(isub_loc)%icgn)
         sub_aux(isub_loc)%ncorners = nintersection
         sub_aux(isub_loc)%licgn    = nintersection
         allocate(sub_aux(isub_loc)%icgn(sub_aux(isub_loc)%licgn))
         sub_aux(isub_loc)%icgn = intersection(1:nintersection)

         ! debug
         !print *,'isub = ', isub, 'ncorners',sub_aux(isub_loc)%ncorners, 'icgn', sub_aux(isub_loc)%icgn
         !call flush(6)

         deallocate(intersection)
         deallocate(cornermerge)
         deallocate(iingns)
      end do

      do isub_loc = 1,lindexsub
         isub = indexsub(isub_loc)
         nnodi  = suba(isub_loc)%nnodi

         globtypes => sub_aux(isub_loc)%globtypes
         globsubs  => sub_aux(isub_loc)%globsubs
         nsubnode  => sub_aux(isub_loc)%nsubnode

         liingns = nnodi
         allocate(iingns(liingns))
         call dd_get_interface_global_numbers(suba(isub_loc), iingns,liingns)

         ! now, it contains all corners from the current and neighbouring subdomains
         ! number of corners the subdomain is responsible for numbering of by the rule:
         ! a corner is given a number by the subdomain with the lowest index
         do icorners = 1,sub_aux(isub_loc)%ncorners

            indg = sub_aux(isub_loc)%icgn(icorners)

            ! get interface index of first node of glob
            call get_index(indg,iingns,liingns,inodi)
            if (inodi .eq. -1) then
               call error(routine_name,'Index of corner not found ', indedges)
            end if

            globtypes(inodi) = 3

         end do

         !print *,'isub = ',isub,'globtypes',globtypes
         !call flush(6)

         deallocate(iingns)

         nullify(globtypes)
         nullify(globsubs)
         nullify(nsubnode)
      end do

      ! number coarse corners based on icgn
      do isub_loc = 1,lindexsub

         sub_aux(isub_loc)%licgcn = sub_aux(isub_loc)%ncorners
         ! Indices of Corners in Global Coarse Numbering
         allocate(sub_aux(isub_loc)%icgcn(sub_aux(isub_loc)%licgcn))

         sub_aux(isub_loc)%global_indices        => sub_aux(isub_loc)%icgn
         sub_aux(isub_loc)%global_coarse_indices => sub_aux(isub_loc)%icgcn
      end do

      call dd_create_globs_number_coarse_quantity_dofs(0, ncorner)

      ! now finish numbering of non-local corners 
      do isub_loc = 1,lindexsub
         isub = indexsub(isub_loc)

         nnodi  = suba(isub_loc)%nnodi

         ! debug
         !print *, 'isub',isub,'ICGN ', sub_aux(isub_loc)%icgn
         !print *, 'isub',isub,'ICGCN', sub_aux(isub_loc)%icgcn
         !call flush(6)

         if (any(sub_aux(isub_loc)%icgcn.eq.0)) then
            call error(routine_name,'Incomplete information in ICGCN array for isub ', isub)
         end if
         
         licnsins = sub_aux(isub_loc)%ncorners
         allocate(icnsins(licnsins))
         icnsins = 0

         do inodi = 1,nnodi
            inds = suba(isub_loc)%iin(inodi)
            indg = suba(isub_loc)%isngn(inds)

            if (any(sub_aux(isub_loc)%icgn.eq.indg)) then

               call get_index_sorted(indg,sub_aux(isub_loc)%icgn,size(sub_aux(isub_loc)%icgn),icscn)
               if (icscn .eq. -1) then
                  write(*,*) 'icgn:',sub_aux(isub_loc)%icgn
                  call error(routine_name,'Global index not found in ICGN', indg)
               end if

               icnsins(icscn) = inodi
            end if
         end do
         ! check the array
         if (any(icnsins.eq.0)) then
            write(*,*) 'isub',isub,'IIN', suba(isub_loc)%isngn(suba(isub_loc)%iin)
            write(*,*) 'isub',isub,'ICNSINS', icnsins
            call error(routine_name, 'array of local corners ICNSINS incomplete')
         end if

         ! upload subdomain corners 
         call dd_upload_sub_corners(suba(isub_loc), sub_aux(isub_loc)%ncorners, &
                                    sub_aux(isub_loc)%icgcn, sub_aux(isub_loc)%ncorners, &
                                    icnsins,licnsins)

         deallocate(icnsins)
      end do

      do isub_loc = 1,lindexsub
         isub = indexsub(isub_loc)

         kglobs    => sub_aux(isub_loc)%kglobs
         globtypes => sub_aux(isub_loc)%globtypes
         globsubs  => sub_aux(isub_loc)%globsubs
         nsubnode  => sub_aux(isub_loc)%nsubnode

         nnodi  = suba(isub_loc)%nnodi
         ncorners  = sub_aux(isub_loc)%ncorners

         ! reinitialize kglobs
         kglobs    = -1
         globtypes = 0

         icnodes   = 0

         ! first mark corners in subdomain interface
         do icorners = 1,ncorners
            indi = suba(isub_loc)%icnsin(icorners)

            kglobs(indi)    = icorners
            globtypes(indi) = 3
         end do

         !print *,'isub = ',isub,'globtypes',globtypes
         !call flush(6)

         ! at this point, local numbering of corners is already sorted before globs

         ! next, find edges on each subdomain
         iedges   = 0
         do inodi = 1,nnodi
            if (kglobs(inodi) .eq. -1) then
               ! node is not assigned to glob yet
               nsubn = nsubnode(inodi)
               if (nsubn.gt.1) then ! it is an edge 
                  iedges = iedges + 1
                  kglobs(inodi)    = ncorners + iedges
                  globtypes(inodi) = 2

                  ! search for nodes shared by the same set of subdomains
                  do jnodi = inodi+1,nnodi
                     ! check that the node is not yet assigned to glob and it is not a corner
                     if (kglobs(jnodi).eq.-1) then
                        ! check that number of subdomains sharing the node matches
                        if (nsubnode(jnodi).eq.nsubn) then
                           ! check that all indices are the same
                           if (all(globsubs(jnodi,1:nsubn,:).eq.globsubs(inodi,1:nsubn,:))) then
                              ! the node belongs to the same glob
                              kglobs(jnodi) = ncorners + iedges
                              ! it is an edge
                              globtypes(jnodi) = 2
                           end if
                        end if
                     end if
                  end do
               end if
            end if
         end do
         nedges   = iedges
         ! find faces each subdomain
         ifaces   = 0
         do inodi = 1,nnodi
            if (kglobs(inodi) .eq. -1) then
               ! node is not assigned to glob yet
               nsubn = nsubnode(inodi)
               if (nsubn.ne.1) then ! it is a face, check it!
                  call error(routine_name,'The only unassigned nodes should belong to a face!')
               end if

               ifaces = ifaces + 1
               icnodes = icnodes + 1
               kglobs(inodi)    = ncorners + nedges + ifaces
               globtypes(inodi) = 1

               ! search for nodes with shared by the same set of subdomains
               nglobn = 1
               do jnodi = inodi+1,nnodi
                  ! check that the node is not yet assigned to glob and it is not a corner
                  if (kglobs(jnodi).eq.-1) then
                     ! check that number of subdomains sharing the node matches
                     if (nsubnode(jnodi).eq.nsubn) then
                        ! check that all indices are the same
                        if (all(globsubs(jnodi,1:nsubn,:).eq.globsubs(inodi,1:nsubn,:))) then
                           ! the node belongs to the same glob
                           kglobs(jnodi) = ncorners + nedges + ifaces
                           ! it is a face
                           globtypes(jnodi) = 1
                        end if
                     end if
                  end if
               end do
            end if
         end do
         nfaces = ifaces

         ! debug
         !print *,'isub = ',isub,'ncorners',ncorners
         !print *,'isub = ',isub,'nedges',nedges
         !print *,'isub = ',isub,'nfaces',nfaces
         !print *,'isub = ',isub,'kglobs',kglobs
         !print *,'isub = ',isub,'globtypes',globtypes
         !call flush(6)
 
         ! take care of Dirichlet BC
         if (remove_bc_nodes) then
            if (suba(isub_loc)%is_bc_present) then
               ! prepare kdofi
               lkdofi = nnodi
               allocate(kdofi(lkdofi))
               if (lkdofi.gt.0) then
                  kdofi(1) = 0
               end if
               do inodi = 2,nnodi
                  inds = suba(isub_loc)%iin(inodi-1)
                  ndofn = suba(isub_loc)%nndf(inds)

                  kdofi(inodi) = kdofi(inodi-1) + ndofn
               end do
               
               ndofi = suba(isub_loc)%ndofi
               lifixi = ndofi
               allocate(ifixi(lifixi))

               call dd_map_sub_to_subi_int(suba(isub_loc),suba(isub_loc)%ifix,suba(isub_loc)%lifix,ifixi,lifixi)

               do inodi = 1,nnodi
                  if (globtypes(inodi) .ne. 3) then
                     inds = suba(isub_loc)%iin(inodi)
                     ndofn = suba(isub_loc)%nndf(inds)
                     pifixi = kdofi(inodi)
                     if (any(ifixi(pifixi+1:pifixi+ndofn).gt.0)) then
                        kglobs(inodi)    = -2
                        globtypes(inodi) = -2
                     end if
                  end if
               end do
               !print *,'isub =',isub,'kglobs:',kglobs
               !print *,'isub =',isub,'globtypes:',globtypes

               ! update number of edges and faces
               iedges = 0
               do iedg = 1,nedges
                  index_old = ncorners + iedg
                  if (any(kglobs.eq.index_old)) then
                     iedges = iedges + 1
                     index_new = ncorners + iedges
                     where (kglobs.eq.index_old) kglobs = index_new
                  end if
               end do
               ifaces = 0
               do ifac = 1,nfaces
                  index_old = ncorners + nedges + ifac
                  if (any(kglobs.eq.index_old)) then
                     ifaces = ifaces + 1
                     index_new = ncorners + iedges + ifaces
                     where(kglobs.eq.index_old) kglobs = index_new
                  end if
               end do
               ! now update numbers
               nedges = iedges
               nfaces = ifaces

               deallocate(ifixi)
               deallocate(kdofi)
            end if
         end if

         ! debug
         !print *,'after removing BC isub = ',isub,'ncorners',ncorners
         !print *,'after removing BC isub = ',isub,'nedges',nedges
         !print *,'after removing BC isub = ',isub,'nfaces',nfaces

         ! make a note on numbers to the structure
         sub_aux(isub_loc)%ncorners = ncorners
         sub_aux(isub_loc)%nedges   = nedges
         sub_aux(isub_loc)%nfaces   = nfaces

         nullify(globtypes)
         nullify(kglobs)
         nullify(globsubs)
         nullify(nsubnode)
      end do

      ! number coarse corners based on icgn
      do isub_loc = 1,lindexsub

         sub_aux(isub_loc)%licgcn = sub_aux(isub_loc)%nedges + sub_aux(isub_loc)%nfaces
         deallocate(sub_aux(isub_loc)%icgcn)
         allocate(sub_aux(isub_loc)%icgcn(sub_aux(isub_loc)%licgcn))

         sub_aux(isub_loc)%global_indices        => sub_aux(isub_loc)%icgn
         sub_aux(isub_loc)%global_coarse_indices => sub_aux(isub_loc)%icgcn
      end do

      ! Prepare lists of global indices of first nodes at edges and faces
      do isub_loc = 1,lindexsub
         isub = indexsub(isub_loc)

         ! load data
         nadj   = suba(isub_loc)%nadj
         nnodi  = suba(isub_loc)%nnodi
         ndim   = suba(isub_loc)%ndim

         ! set pointers
         lnsubnode = sub_aux(isub_loc)%lnsubnode 
         nsubnode => sub_aux(isub_loc)%nsubnode

         lglobsubs1 = sub_aux(isub_loc)%lglobsubs1
         lglobsubs2 = sub_aux(isub_loc)%lglobsubs2
         globsubs => sub_aux(isub_loc)%globsubs

         lglobtypes = sub_aux(isub_loc)%lglobtypes
         globtypes => sub_aux(isub_loc)%globtypes

         lkglobs = sub_aux(isub_loc)%lkglobs
         kglobs => sub_aux(isub_loc)%kglobs

         ncorners = sub_aux(isub_loc)%ncorners 
         nedges   = sub_aux(isub_loc)%nedges   
         nfaces   = sub_aux(isub_loc)%nfaces   

         ! prepare list of edges and faces at each subdomain
         sub_aux(isub_loc)%liegn  = nedges
         allocate(sub_aux(isub_loc)%iegn(sub_aux(isub_loc)%liegn))
         do iedges = 1,nedges
            indedges = ncorners + iedges

            ! get interface index of first node of glob
            !call get_index(indedges,kglobs,lkglobs,inodi)
            !if (inodi .eq. -1) then
            !   call error(routine_name,'Index of edge not found in kglobs ', indedges)
            !end if

            !inds = suba(isub_loc)%iin(inodi)
            !indg = suba(isub_loc)%isngn(inds)

            ! add this edge index to the list
            minindex = minval(suba(isub_loc)%isngn(suba(isub_loc)%iin),kglobs == indedges)
            sub_aux(isub_loc)%iegn(iedges) = minindex
         end do
         sub_aux(isub_loc)%lifgn  = nfaces
         allocate(sub_aux(isub_loc)%ifgn(sub_aux(isub_loc)%lifgn))
         do ifaces = 1,nfaces
            indfaces = ncorners + nedges + ifaces

            ! get interface index of first node of glob
            !call get_index(indfaces,kglobs,lkglobs,inodi)
            !if (inodi .eq. -1) then
            !   call error(routine_name,'Index of face not found in kglobs ', indfaces)
            !end if

            !inds = suba(isub_loc)%iin(inodi)
            !indg = suba(isub_loc)%isngn(inds)

            !print *, 'isub',isub,'faces position', minindex, 'vs. ',indg

            ! add this face index to the list
            minindex = minval(suba(isub_loc)%isngn(suba(isub_loc)%iin),kglobs == indfaces)
            sub_aux(isub_loc)%ifgn(ifaces) = minindex
         end do

         nullify(globtypes)
         nullify(kglobs)
         nullify(globsubs)
         nullify(nsubnode)
      end do

      ! number coarse edges based on iegn
      do isub_loc = 1,lindexsub

         sub_aux(isub_loc)%liegcn = sub_aux(isub_loc)%nedges
         ! Indices of Edges in Global Coarse Numbering
         allocate(sub_aux(isub_loc)%iegcn(sub_aux(isub_loc)%liegcn))
         sub_aux(isub_loc)%iegcn = 0

         sub_aux(isub_loc)%global_indices        => sub_aux(isub_loc)%iegn
         sub_aux(isub_loc)%global_coarse_indices => sub_aux(isub_loc)%iegcn

         ! debug
         !print *, 'isub',isub,'IEGN ', sub_aux(isub_loc)%iegn
         !print *, 'isub',isub,'IEGCN', sub_aux(isub_loc)%iegcn
         !call flush(6)
      end do

      call dd_create_globs_number_coarse_quantity_dofs(ncorner, nedge)

      ! number coarse faces based on ifgn
      do isub_loc = 1,lindexsub

         ! debug
         !print *, 'isub',isub,'IEGN ', sub_aux(isub_loc)%iegn
         !print *, 'isub',isub,'IEGCN', sub_aux(isub_loc)%iegcn
         !call flush(6)

         sub_aux(isub_loc)%lifgcn = sub_aux(isub_loc)%nfaces
         ! Indices of Faces in Global Coarse Numbering
         allocate(sub_aux(isub_loc)%ifgcn(sub_aux(isub_loc)%lifgcn))
         sub_aux(isub_loc)%ifgcn = 0

         sub_aux(isub_loc)%global_indices        => sub_aux(isub_loc)%ifgn
         sub_aux(isub_loc)%global_coarse_indices => sub_aux(isub_loc)%ifgcn
      end do

      call dd_create_globs_number_coarse_quantity_dofs(ncorner + nedge, nface)

      ! debug
      !do isub_loc = 1,lindexsub
      !   print *, 'isub',isub,'IFGN ', sub_aux(isub_loc)%ifgn
      !   print *, 'isub',isub,'IFGCN', sub_aux(isub_loc)%ifgcn
      !   call flush(6)
      !end do

      ! now all processes has global array of first indices of edges and faces available, perform localization of globs
      do isub_loc = 1,lindexsub
         isub = indexsub(isub_loc)

         ! load data
         nadj   = suba(isub_loc)%nadj
         nnodi  = suba(isub_loc)%nnodi
         ndim   = suba(isub_loc)%ndim

         ! set pointers
         lnsubnode = sub_aux(isub_loc)%lnsubnode 
         nsubnode => sub_aux(isub_loc)%nsubnode

         lglobsubs1 = sub_aux(isub_loc)%lglobsubs1
         lglobsubs2 = sub_aux(isub_loc)%lglobsubs2
         globsubs => sub_aux(isub_loc)%globsubs

         lglobtypes = sub_aux(isub_loc)%lglobtypes
         globtypes => sub_aux(isub_loc)%globtypes

         lkglobs = sub_aux(isub_loc)%lkglobs
         kglobs => sub_aux(isub_loc)%kglobs

         ncorners = sub_aux(isub_loc)%ncorners  
         nedges   = sub_aux(isub_loc)%nedges  
         nfaces   = sub_aux(isub_loc)%nfaces  

         ! prepare array for global indices of interface nodes 
         liingns = nnodi
         allocate(iingns(liingns))
         call dd_get_interface_global_numbers(suba(isub_loc), iingns,liingns)

         ! mapping of globs
         nglobs = nedges + nfaces
         lglobal_glob_numbers = nglobs
         allocate(global_glob_numbers(lglobal_glob_numbers))
         lnglobvars = nglobs
         allocate(nglobvars(lnglobvars))
         nglobvars   = 0
         lnglobnodess = nglobs
         allocate(nglobnodess(lnglobnodess))
         nglobnodess = 0
         lglob_types = nglobs
         allocate(glob_types(lglob_types))

         ! put edges
         global_glob_numbers(:nedges)   = sub_aux(isub_loc)%iegcn
         glob_types(:nedges)   = 2
         ! put faces
         global_glob_numbers(nedges+1:) = sub_aux(isub_loc)%ifgcn
         glob_types(nedges+1:) = 1

         ! compute number of nodes and variables at individual globs
         do inodi = 1,nnodi
            if (globtypes(inodi) .eq. 1 .or. globtypes(inodi) .eq. 2 ) then
               ! it is a face or edge

               indglb = kglobs(inodi) - ncorners

               inds = suba(isub_loc)%iin(inodi)
               ndofn = suba(isub_loc)%nndf(inds)

               nglobnodess(indglb) = nglobnodess(indglb) + 1
               nglobvars(indglb)   = nglobvars(indglb)  + ndofn
            end if
         end do

         ! debug
         !print *,'globtypes',globtypes
         !print *,'nglobnodess',nglobnodess
         !print *,'nglobvars',nglobvars
         ! debug
         !print *,'global_glob_numbers',global_glob_numbers

         ligvsivns1 = nglobs
         ligvsivns2 = maxval(nglobvars)
         allocate(igvsivns(ligvsivns1,ligvsivns2))
         igvsivns(:,:) = 0
         lignsins1 = nglobs
         lignsins2 = maxval(nglobnodess)
         allocate(ignsins(lignsins1,lignsins2))
         ignsins(:,:) = 0

         lkdofi = nnodi
         allocate(kdofi(lkdofi))
         if (lkdofi.gt.0) then
            kdofi(1) = 0
         end if
         do inodi = 2,nnodi
            inds = suba(isub_loc)%iin(inodi-1)
            ndofn = suba(isub_loc)%nndf(inds)

            kdofi(inodi) = kdofi(inodi-1) + ndofn
         end do
         ! debug
         !print *,'kdofi',kdofi
         !call flush(6)

         ! use nglobnodess as counter
         nglobnodess(:) = 0
         nglobvars(:) = 0
         do inodi = 1,nnodi
            if (globtypes(inodi) .eq. 1 .or. globtypes(inodi) .eq. 2 ) then
               ! it is a face or edge
               indglb = kglobs(inodi) - ncorners

               inds = suba(isub_loc)%iin(inodi)
               ndofn = suba(isub_loc)%nndf(inds)


               nglobnodess(indglb) = nglobnodess(indglb) + 1
               ignsins(indglb,nglobnodess(indglb)) = inodi

               do idofn = 1,ndofn
                  nglobvars(indglb) = nglobvars(indglb)  + 1

                  indivs = kdofi(inodi) + idofn

                  igvsivns(indglb,nglobvars(indglb)) = indivs
               end do
            end if
         end do
         deallocate(kdofi)
         ! debug
         !print *,'igvsivns'
         !do i = 1,ligvsivns1
         !   print *,(igvsivns(i,j),j = 1,nglobvars(i))
         !end do

         lnsubglobs = nglobs
         allocate(nsubglobs(lnsubglobs))
         do iglob = 1,nglobs
            if (nglobnodess(iglob) .eq. 0) then
               call error( routine_name, 'Glob with no nodes?' )
            end if
            indi = ignsins(iglob,1)

            nsubglobs(iglob) = nsubnode(indi)
         end do

         lglob_subs1 = nglobs
         lglob_subs2 = maxval(nsubglobs)
         allocate(glob_subs(lglob_subs1,lglob_subs2))
         glob_subs = 0
         do iglob = 1,nglobs
            indi = ignsins(iglob,1)

            do i = 1,nsubnode(indi)
               glob_subs(iglob,i) = globsubs(indi,i,1)
            end do
         end do

         !print *, 'nsubglobs', nsubglobs
         !print *, 'glob_subs', glob_subs

         call dd_upload_sub_globs(suba(isub_loc), nglobs, global_glob_numbers,lglobal_glob_numbers,&
                                  nglobnodess,lnglobnodess, nglobvars,lnglobvars,&
                                  ignsins,lignsins1,lignsins2, igvsivns,ligvsivns1,ligvsivns2,&
                                  glob_types,lglob_types, &
                                  nsubglobs,lnsubglobs,&
                                  glob_subs,lglob_subs1,lglob_subs2)
         deallocate(nsubglobs)
         deallocate(glob_subs)
         deallocate(ignsins)
         deallocate(igvsivns)
         deallocate(glob_types)
         deallocate(nglobnodess)
         deallocate(nglobvars)
         deallocate(global_glob_numbers)

         deallocate(iingns)

         nullify(globtypes)
         nullify(kglobs)
         nullify(globsubs)
         nullify(nsubnode)
      end do

      ! get global size of interface problem
!*****************************************************************MPI
      call MPI_ALLREDUCE(nnodi_loc,nnodig,1, MPI_INTEGER, MPI_SUM, comm_all, ierr) 
      call MPI_ALLREDUCE(ndofi_loc,ndofig,1, MPI_INTEGER, MPI_SUM, comm_all, ierr) 
!*****************************************************************MPI

      ! root print the summary of selection of globs
      if (myid.eq.0) then
         call info( routine_name, '   Total size of interface problem: ' )
         call info( routine_name, '    number of unknowns on interface  = ',ndofig )
         call info( routine_name, '    number of nodes on interface     = ',nnodig )
         call info( routine_name, '   The following globs were recognized: ' )
         call info( routine_name, '    nface   = ',nface )
         call info( routine_name, '    nedge   = ',nedge )
         call info( routine_name, '    ncorner = ',ncorner )
      end if

      ! clear auxiliary structure
      do isub_loc = 1,lindexsub
         deallocate(sub_aux(isub_loc)%nsubnode)
         deallocate(sub_aux(isub_loc)%globsubs)
         deallocate(sub_aux(isub_loc)%globtypes)
         deallocate(sub_aux(isub_loc)%kglobs)
         deallocate(sub_aux(isub_loc)%icgn)
         deallocate(sub_aux(isub_loc)%icgcn)
         deallocate(sub_aux(isub_loc)%iegn)
         deallocate(sub_aux(isub_loc)%iegcn)
         deallocate(sub_aux(isub_loc)%ifgn)
         deallocate(sub_aux(isub_loc)%ifgcn)
      end do
      deallocate(sub_aux)

contains

!**********************************************
      subroutine dd_create_globs_select_corners
!**********************************************
      ! routine for selection of corners for one subdomain
      ! This is a local subroutine with access to global data of DD_CREATE_GLOBS

      implicit none

      ! now with subdomains sharing a face, run the corner selecting algorithm
      kishnadj  = 0
      do ia = 1,nadj
         ! get index of neighbour
         isubadj = suba(isub_loc)%iadj(ia)

         ! number of nodes shared with this subdomain
         nshn = suba(isub_loc)%nshnadj(ia)

         we_share_a_face = .false.
         do ishn = 1,nshn
            ! index of shared node in interface numbering
            indshn = suba(isub_loc)%ishnadj(kishnadj + ishn)

            if (globtypes(indshn).eq.1) then
               we_share_a_face = .true.
               exit
            end if
         end do

         if (we_share_a_face) then

            ! check components of graph
            ! combined interface components - may be multiple of the above
            linterface_components = nshn
            allocate(interface_components(linterface_components))


            ! determine continuity of components
            if (check_face_components) then
               ! my own components
               lmy_interface_components = nshn
               allocate(my_interface_components(lmy_interface_components))
               ! interface components of adjacent subdomain
               ladj_interface_components = nshn
               allocate(adj_interface_components(ladj_interface_components))

               ! create my own interface components
               do ishn = 1,nshn
                  ! index of shared node in interface numbering
                  indshn   = suba(isub_loc)%ishnadj(kishnadj + ishn)

                  ! index of shared node in subdomain numbering
                  indshsub = suba(isub_loc)%iin(indshn)

                  ! index of interface component of adjacent subdomain
                  indcomponentadj = suba(isub_loc)%ishnncadj(kishnadj + ishn)

                  my_interface_components(ishn)  = suba(isub_loc)%nodal_components(indshsub)
                  adj_interface_components(ishn) = indcomponentadj
               end do
               ! both arrays can have holes - not all components of
               ! adjacent subdomain have common interface with actual subdomain
               ! renumber components from 1 to number_of_different_indices
               indnew = 0
               do iold = 1,maxval(my_interface_components)
                  if (any(my_interface_components.eq.iold)) then
                     ! increase new index
                     indnew = indnew + 1
                     where (my_interface_components.eq.iold) my_interface_components = indnew
                  end if
               end do
               nmy_components = indnew
               indnew = 0
               do iold = 1,maxval(adj_interface_components)
                  if (any(adj_interface_components.eq.iold)) then
                     ! increase new index
                     indnew = indnew + 1
                     where (adj_interface_components.eq.iold) adj_interface_components = indnew
                  end if
               end do
               nadj_components = indnew

               ! combine the components
               interface_components = (my_interface_components - 1) * nadj_components + adj_interface_components
               ! exclude repetitions
               indnew = 0
               do iold = 1,maxval(interface_components)
                  if (any(interface_components.eq.iold)) then
                     ! increase new index
                     indnew = indnew + 1
                     where (interface_components.eq.iold) interface_components = indnew
                  end if
               end do
               ninterface_components = indnew

               if (ninterface_components.gt.1) then
                   call info(routine_name,'disconnected interface - number of interface components:',ninterface_components)
               end if

               deallocate(my_interface_components)
               deallocate(adj_interface_components)
            else
               interface_components  = 1
               ninterface_components = 1
            end if


            ! perform the triangle check for each component
            do iinterface_component = 1,ninterface_components

               ! initialize array
               inodcf = 0

               interface_componentsize = count(interface_components .eq. iinterface_component)

               lcompind = interface_componentsize
               allocate(compind(lcompind))
               
               indcomp = 0
               do ishn = 1,nshn
                  if (interface_components(ishn) .eq. iinterface_component) then
                     indcomp = indcomp + 1

                     compind(indcomp) = ishn
                  end if
               end do
               if (indcomp.ne.lcompind) then
                  call error(routine_name,'dimension mismatch in component size',indcomp)
               end if

               ! prepare coordinates of common interface
               lxyzsh1 = interface_componentsize
               lxyzsh2 = ndim
               ! localize coordinates of shared nodes
               allocate(xyzsh(lxyzsh1,lxyzsh2))
               do incomp = 1,interface_componentsize
                  ishn  = compind(incomp)
                  indni = suba(isub_loc)%ishnadj(kishnadj + ishn)
                  indn  = suba(isub_loc)%iin(indni)
                  xyzsh(incomp,1:ndim) = suba(isub_loc)%xyz(indn,1:ndim)
               end do

               ! search an optimal triangle
               inodcf = 0
               if (interface_componentsize.eq.0) then
                  call error(routine_name,'There appears that there are zero shared nodes for subdomain:',isub)
               end if
               if (interface_componentsize.gt.0) then
                  lxyzbase = ndim
                  allocate(xyzbase(lxyzbase))
                  ! make it the most remote node to the first interface node
                  inodshaux = 1
                  xyzbase(1:ndim) = xyzsh(inodshaux,1:ndim)
               
                  ! Find second corner by maximizing the distance of the first shared node
                  ldist = interface_componentsize
                  allocate(dist(ldist))
                  do incomp = 1,interface_componentsize
                     dist(incomp) = sum((xyzsh(incomp,1:ndim) - xyzbase(1:ndim))**2)
                  end do
                  indaux  = maxloc(dist)

                  ! set index of first new corner
                  inodcf(1) = indaux(1)

                  deallocate(xyzbase)
                  deallocate(dist)
               end if
               if (meshdim.gt.1.and.interface_componentsize.gt.1) then
                  lxyzbase = ndim
                  allocate(xyzbase(lxyzbase))
                  ! one corner is already selected, select the second
                  xyzbase(1:ndim) = xyzsh(inodcf(1),1:ndim)
               
                  ! Find second corner by maximizing the distance from the first one
                  ldist = interface_componentsize
                  allocate(dist(ldist))
                  do incomp = 1,interface_componentsize
                     dist(incomp) = sum((xyzsh(incomp,1:ndim) - xyzbase(1:ndim))**2)
                  end do
                  indaux  = maxloc(dist)

                  ! set index of second new corner
                  inodcf(2) = indaux(1)

                  ! if geometry fails, select any node different from the first, regardless of the coordinates
                  if (inodcf(2).eq.inodcf(1)) then
                     !write(*,*) 'dist:',dist
                     !write(*,*) 'coords:'
                     !do incomp = 1,interface_componentsize
                     !   write(*,*) xyzsh(incomp,1:ndim)
                     !end do
                     !write(*,*) 'base:'
                     !write(*,*) xyzbase(1:ndim)
                     !call flush(6)
                     call warning(routine_name, 'Problem finding second corner on subdomain - same as first. ',&
                                  inodcf(2))
                     ! perform search for different corner
                     do incomp = 1,interface_componentsize
                        if ( incomp .ne. inodcf(1) ) then
                           inodcf(2) = incomp
                           exit
                        end if
                     end do
                  end if

                  deallocate(xyzbase)
                  deallocate(dist)
               end if
               if (meshdim.gt.2.and.interface_componentsize.gt.2) then
               ! two corners are already set, select the third
                  x1 = xyzsh(inodcf(1),1)
                  y1 = xyzsh(inodcf(1),2)
                  z1 = xyzsh(inodcf(1),3)
                  x2 = xyzsh(inodcf(2),1)
                  y2 = xyzsh(inodcf(2),2)
                  z2 = xyzsh(inodcf(2),3)
               
                  ! Find third corner as the one maximizing area of triangle
                  larea = interface_componentsize
                  allocate(area(larea))
                  do incomp = 1,interface_componentsize
                     xish = xyzsh(incomp,1)
                     yish = xyzsh(incomp,2)
                     zish = xyzsh(incomp,3)
                     area(incomp) = ((y2-y1)*(zish-z1) - (yish-y1)*(z2-z1))**2 &
                                  + ((x2-x1)*(zish-z1) - (xish-x1)*(z2-z1))**2 &
                                  + ((x2-x1)*(yish-y1) - (xish-x1)*(y2-y1))**2 
                  end do
                  ! find maximum in the area array - perform two checks in reverse order
                  do iround = 1,2
                     if (iround.eq.1) then
                        start  = 1
                        finish = larea
                        step   = 1
                     else 
                        start  = larea
                        finish = 1
                        step   = -1
                     end if
                     ! search maximum
                     maxv = 0._kr
                     maxl = 0
                     do i = start,finish,step
                        if (area(i).gt.maxv) then
                           maxv = area(i)
                           maxl = i
                        end if
                     end do
                     indaux(1) = maxl
                     !indaux  = maxloc(area)

                     ! set index of the third new corner
                     inodcf(2+iround) = indaux(1)
                  end do

                  ! if geometry fails, select any node different from the first two, regardless of the coordinates
                  if ((inodcf(3).eq.inodcf(1).and.inodcf(4).eq.inodcf(1)).or.&
                      (inodcf(3).eq.inodcf(2).and.inodcf(4).eq.inodcf(2))) then
                     !write(*,*) 'area:',area
                     !call flush(6)
                     call warning(routine_name, &
                                  'Problem finding third corner on subdomain - same as first or second.', &
                                  inodcf(3))
                     ! perform search for different corners
                     do incomp = 1,interface_componentsize
                        if ( incomp .ne. inodcf(1) .and. incomp .ne. inodcf(2) ) then
                           inodcf(3) = incomp
                           exit
                        end if
                     end do
                     do incomp = 1,interface_componentsize
                        if ( incomp .ne. inodcf(1) .and. &
                             incomp .ne. inodcf(2) .and. &
                             incomp .ne. inodcf(3) ) then
                           inodcf(4) = incomp
                           exit
                        end if
                     end do
                  end if

                  deallocate(area)
               end if

               ! mark corner in subdomain interface
               !print *,'inodcf',suba(isub_loc)%isngn(suba(isub_loc)%iin(suba(isub_loc)%ishnadj(kishnadj + inodcf)))
               if (meshdim.le.2) then
                  end_bound = 2
               else
                  end_bound = 4
               end if
               do i = 1,end_bound
                  if (inodcf(i).ne.0) then

                     indni = suba(isub_loc)%ishnadj(kishnadj + compind(inodcf(i)))

                     globtypes(indni) = 3
                  end if
               end do

               deallocate(xyzsh)
               deallocate(compind)
            end do


            deallocate(interface_components)

         end if

         kishnadj = kishnadj + nshn
      end do
      end subroutine

!***************************************************************************
      subroutine dd_create_globs_number_coarse_quantity_dofs(shift, ncoarse)
!***************************************************************************
! routine for distributed numbering of coarse degrees of freedom
! This is a local subroutine with access to global data of DD_CREATE_GLOBS

      implicit none
      ! shift of coarse numbering 
      ! shift = 0                for corners 
      ! shift = ncorner          for edges 
      ! shift = ncorner + nedge  for edges 
      integer, intent(in)  :: shift
      ! resulting number of coarse entities (corners, edges, or faces)
      integer, intent(out) :: ncoarse
      ! local vars
      integer :: nsub, icoarses

      integer ::            array_length
      integer ::            indnc

      integer ::             ncoarse_unique
      integer ::             lncoarse_uniquea
      integer,allocatable ::  ncoarse_uniquea(:)
      integer,allocatable ::  ncoarse_uniquea_aux(:)


      nsub     = sub2proc(nproc+1) - 1
      lncoarse_uniquea = nsub + 1

      allocate(ncoarse_uniquea(lncoarse_uniquea))
      allocate(ncoarse_uniquea_aux(lncoarse_uniquea))

      ncoarse_uniquea_aux = 0

      do isub_loc = 1,lindexsub
         isub = indexsub(isub_loc)
         nnodi  = suba(isub_loc)%nnodi

         globtypes => sub_aux(isub_loc)%globtypes
         globsubs  => sub_aux(isub_loc)%globsubs
         nsubnode  => sub_aux(isub_loc)%nsubnode

         liingns = nnodi
         allocate(iingns(liingns))
         call dd_get_interface_global_numbers(suba(isub_loc), iingns,liingns)

         ! now, it contains all coarse entities from the current and neighbouring subdomains
         ! number of coarse entities the subdomain is responsible for numbering of by the rule:
         ! a coarse entities is given a number by the subdomain with the lowest index
         ncoarse_unique = 0
         do icoarses = 1,size(sub_aux(isub_loc)%global_indices)

            indg = sub_aux(isub_loc)%global_indices(icoarses)

            ! get interface index of first node of glob
            call get_index(indg,iingns,liingns,inodi)
            if (inodi .eq. -1) then
               call error(routine_name,'Index of coarse entity not found ', indedges)
            end if

            if ( min(isub, minval(globsubs(inodi,1:nsubnode(inodi),1))) .eq. isub ) then
               ! I am the subdomain with the smallest number that shares this node
               ncoarse_unique = ncoarse_unique + 1
            end if
         end do

         ncoarse_uniquea_aux(isub) = ncoarse_unique

         deallocate(iingns)

         nullify(globtypes)
         nullify(globsubs)
      end do
!*****************************************************************MPI
      call MPI_ALLREDUCE(ncoarse_uniquea_aux,ncoarse_uniquea,lncoarse_uniquea, &
                         MPI_INTEGER, MPI_SUM, comm_all, ierr) 
!*****************************************************************MPI
      deallocate(ncoarse_uniquea_aux)

      ! now array NCOARSE_UNIQUEA contains number of unique coarse dof for each
      ! subdomain, change it to starts
      ncoarse_uniquea(2:) = ncoarse_uniquea(:nsub)
      ncoarse_uniquea(1)  = 1
      do i = 2,lncoarse_uniquea
         ncoarse_uniquea(i) = ncoarse_uniquea(i-1) + ncoarse_uniquea(i)
      end do
      ! global number of coarse entities
      ncoarse = ncoarse_uniquea(lncoarse_uniquea) - 1

      ! now number local unique coarse entities 
      do isub_loc = 1,lindexsub
         isub = indexsub(isub_loc)
         nnodi  = suba(isub_loc)%nnodi

         globtypes => sub_aux(isub_loc)%globtypes
         globsubs  => sub_aux(isub_loc)%globsubs
         nsubnode  => sub_aux(isub_loc)%nsubnode

         sub_aux(isub_loc)%global_coarse_indices   = 0

         ! number of coarse entities the subdomain is responsible for numbering of by the rule:
         ! a coarse entities is given a number by the subdomain with the lowest index
         indcg = ncoarse_uniquea(isub)
         do inodi = 1,nnodi
            inds = suba(isub_loc)%iin(inodi)
            indg = suba(isub_loc)%isngn(inds)

            if (any(sub_aux(isub_loc)%global_indices.eq.indg)) then

               if ( min(isub, minval(globsubs(inodi,1:nsubnode(inodi),1))) .eq. isub ) then
                  ! I am the subdomain with the smallest number that shares this node

                  !call get_index_sorted(indg,sub_aux(isub_loc)%global_indices,size(sub_aux(isub_loc)%global_indices),icscn)
                  ! the array may not be sorted
                  call get_index(indg,sub_aux(isub_loc)%global_indices,size(sub_aux(isub_loc)%global_indices),icscn)
                  if (icscn .eq. -1) then
                     write(*,*) 'INDG',indg,'GLOBAL_INDICES',sub_aux(isub_loc)%global_indices
                     call error(routine_name,'Global index not found in first call to GLOBAL_INDICES ', indg)
                  end if

                  sub_aux(isub_loc)%global_coarse_indices(icscn) = shift + indcg

                  indcg = indcg + 1
               end if
            end if
         end do

         nullify(globtypes)
         nullify(globsubs)
      end do

      deallocate(ncoarse_uniquea)

      ! now we need to exchange the arrays of corners among all neighbours - first exchange their numbers
      do isub_loc = 1,lindexsub
         isub = indexsub(isub_loc)

         sub_aux(isub_loc)%comm_array_in => sub_aux(isub_loc)%global_indices
      end do

      call dd_interchange_integer_arrays(suba,lsuba, &
                                         sub2proc,lsub2proc, indexsub,lindexsub, &
                                         comm_all, &
                                         sub_aux,lsub_aux)

      do isub_loc = 1,lindexsub
         isub = indexsub(isub_loc)

         array_length = size(sub_aux(isub_loc)%comm_array_out)
         allocate(sub_aux(isub_loc)%comm_gathered_array1(array_length))
         sub_aux(isub_loc)%comm_gathered_array1 = sub_aux(isub_loc)%comm_array_out
         deallocate(sub_aux(isub_loc)%comm_array_out)
         deallocate(sub_aux(isub_loc)%comm_array_numbers)
      end do

      do isub_loc = 1,lindexsub
         isub = indexsub(isub_loc)
         sub_aux(isub_loc)%comm_array_in => sub_aux(isub_loc)%global_coarse_indices
      end do

      call dd_interchange_integer_arrays(suba,lsuba, &
                                         sub2proc,lsub2proc, indexsub,lindexsub, &
                                         comm_all, &
                                         sub_aux,lsub_aux)

      do isub_loc = 1,lindexsub
         isub = indexsub(isub_loc)

         array_length = size(sub_aux(isub_loc)%comm_array_out)
         allocate(sub_aux(isub_loc)%comm_gathered_array2(array_length))
         sub_aux(isub_loc)%comm_gathered_array2 = sub_aux(isub_loc)%comm_array_out
         deallocate(sub_aux(isub_loc)%comm_array_out)
         deallocate(sub_aux(isub_loc)%comm_array_numbers)
      end do

      ! process the array
      do isub_loc = 1,lindexsub
         isub = indexsub(isub_loc)

         !print *, 'comm_gathered_array1', sub_aux(isub_loc)%comm_gathered_array1
         !print *, 'comm_gathered_array2', sub_aux(isub_loc)%comm_gathered_array2
         !call flush(6)

         array_length = count(sub_aux(isub_loc)%comm_gathered_array2 .gt. 0)
         allocate(sub_aux(isub_loc)%comm_array_result1(array_length))
         allocate(sub_aux(isub_loc)%comm_array_result2(array_length))

         indnc = 0
         do i = 1,size(sub_aux(isub_loc)%comm_gathered_array1)

            if (sub_aux(isub_loc)%comm_gathered_array2(i).gt.0) then
               indnc = indnc + 1

               sub_aux(isub_loc)%comm_array_result1(indnc) = sub_aux(isub_loc)%comm_gathered_array1(i)
               sub_aux(isub_loc)%comm_array_result2(indnc) = sub_aux(isub_loc)%comm_gathered_array2(i)
            end if
         end do

         ! sort the arrays
         call iquick_sort_simultaneous(sub_aux(isub_loc)%comm_array_result1, size(sub_aux(isub_loc)%comm_array_result1), &
                                       sub_aux(isub_loc)%comm_array_result2, size(sub_aux(isub_loc)%comm_array_result2))

         ! now get through my array and fill the holes
         do i = 1,size(sub_aux(isub_loc)%global_indices)
            indg  = sub_aux(isub_loc)%global_indices(i)

            if (sub_aux(isub_loc)%global_coarse_indices(i).eq.0) then
               ! this point was numbered by a neighbour
               call get_index_sorted(indg,sub_aux(isub_loc)%comm_array_result1, size(sub_aux(isub_loc)%comm_array_result1),icscn)
               if (icscn .eq. -1) then
                  call error(routine_name,'Global index not found in second call to COMM_ARRAY_RESULT1 ', indg)
               end if

               sub_aux(isub_loc)%global_coarse_indices(i) = sub_aux(isub_loc)%comm_array_result2(icscn)
            end if
         end do

         deallocate(sub_aux(isub_loc)%comm_gathered_array1)
         deallocate(sub_aux(isub_loc)%comm_gathered_array2)
         deallocate(sub_aux(isub_loc)%comm_array_result1)
         deallocate(sub_aux(isub_loc)%comm_array_result2)
      end do

      end subroutine

end subroutine

!*******************************************************************************
subroutine dd_interchange_integer_arrays(suba,lsuba, &
                                         sub2proc,lsub2proc, indexsub,lindexsub, &
                                         comm_all, &
                                         sub_aux,lsub_aux)
!*******************************************************************************
! Procedure for distributed interchange of integer data of subdomain and its neighbours. 
! This is a local subroutine with access to global data of DD_CREATE_GLOBS

      use module_utils
      use module_pp, only : pp_get_proc_for_sub, pp_get_unique_tag
      implicit none
      include "mpif.h"

      integer,intent(in) ::             lsuba
      type(subdomain_type),intent(in) :: suba(lsuba)

      integer,intent(in) :: lsub2proc
      integer,intent(in) ::  sub2proc(lsub2proc)

      integer,intent(in) :: lindexsub
      integer,intent(in) ::  indexsub(lindexsub)

      integer,intent(in) :: comm_all ! MPI communicator

      integer,intent(in) ::              lsub_aux
      type(sub_aux_type),intent(inout) :: sub_aux(lsub_aux)

      ! local vars
      character(*),parameter:: routine_name = 'DD_INTERCHANGE_INTEGER_ARRAYS'
      integer ::            number_to_send, number_to_receive
      integer ::            array_length, kineibarray

      integer             :: nadj, nadjx, ia
      integer             :: isub, isubadj, isubadj_loc
      integer             :: isub_loc, nsub_loc

      ! MPI related arrays and variables
      integer ::            myid, nproc
      integer ::            ierr
      integer ::            tag
      integer ::            nreq, ireq
      integer ::            lrequest
      integer,allocatable :: request(:)
      integer             :: lstatarray1
      integer             :: lstatarray2
      integer,allocatable :: statarray(:,:)
      integer ::            neibproc

      ! orient in communicator
      call MPI_COMM_RANK(comm_all,myid,ierr)
      call MPI_COMM_SIZE(comm_all,nproc,ierr)

      nsub_loc = sub2proc(myid+2) - sub2proc(myid+1) 

      nadjx = maxval(suba%nadj)
      lrequest = 2*nsub_loc*nadjx
      allocate(request(lrequest))
      lstatarray1 = MPI_STATUS_SIZE
      lstatarray2 = 2*nsub_loc*nadjx
      allocate(statarray(lstatarray1,lstatarray2))

      ! now we need to exchange the arrays of corners among all neighbours - first exchange their numbers
      ireq = 0
      do isub_loc = 1,lindexsub
         isub = indexsub(isub_loc)

         nadj = suba(isub_loc)%nadj

         ! send/receive data from all neighbours
         allocate(sub_aux(isub_loc)%comm_array_numbers(nadj))
         ! initialize the array
         sub_aux(isub_loc)%comm_array_numbers = -1

         do ia = 1,nadj

            ! get index of neighbour
            isubadj = suba(isub_loc)%iadj(ia)
            call pp_get_proc_for_sub(isubadj,comm_all,sub2proc,lsub2proc,neibproc)
   
            if (neibproc /= myid) then
               ! pass messages
               !ireq = ireq + 1
               !call pp_get_unique_tag(isub,isubadj,comm_all,sub2proc,lsub2proc, tag)
               !call MPI_ISEND(number_to_send,   1,MPI_INTEGER,neibproc,tag, comm_all,request(ireq),ierr)
               call pp_get_unique_tag(isubadj,isub,comm_all,sub2proc,lsub2proc, tag)
               ireq = ireq + 1
               call MPI_IRECV(sub_aux(isub_loc)%comm_array_numbers(ia),1,MPI_INTEGER,neibproc,tag, comm_all,request(ireq),ierr)
            else
               ! simple copy in memory
               ! find index of the subdomain
               call get_index(isubadj,indexsub,lindexsub,isubadj_loc)
               if (isubadj_loc .eq. -1) then
                  write(*,*) 'indexsub',indexsub
                  call error(routine_name, 'Index of subdomain not found.',isubadj)
               end if 
               sub_aux(isub_loc)%comm_array_numbers(ia) = size(sub_aux(isubadj_loc)%comm_array_in)
            end if

         end do
      end do
      do isub_loc = 1,lindexsub
         isub = indexsub(isub_loc)

         nadj = suba(isub_loc)%nadj

         number_to_send = size(sub_aux(isub_loc)%comm_array_in)

         do ia = 1,nadj

            ! get index of neighbour
            isubadj = suba(isub_loc)%iadj(ia)
            call pp_get_proc_for_sub(isubadj,comm_all,sub2proc,lsub2proc,neibproc)
   
            if (neibproc /= myid) then
               ! pass messages
               call pp_get_unique_tag(isub,isubadj,comm_all,sub2proc,lsub2proc, tag)
               ireq = ireq + 1
               call MPI_ISEND(number_to_send,   1,MPI_INTEGER,neibproc,tag, comm_all,request(ireq),ierr)
            end if

         end do
      end do
      nreq = ireq
      if (nreq.gt.0) then
         call MPI_WAITALL(nreq, request, statarray, ierr)
      end if

      ! now we need to exchange the arrays of corners among all neighbours - now exchange indices
      ireq = 0
      do isub_loc = 1,lindexsub
         isub = indexsub(isub_loc)

         nadj = suba(isub_loc)%nadj

         !number_to_send = size(sub_aux(isub_loc)%comm_array_in)

         array_length = sum(sub_aux(isub_loc)%comm_array_numbers)
         allocate(sub_aux(isub_loc)%comm_array_out(array_length))

         ! send/receive data from all neighbours
         kineibarray  = 1
         do ia = 1,nadj

            ! get index of neighbour
            isubadj = suba(isub_loc)%iadj(ia)
            call pp_get_proc_for_sub(isubadj,comm_all,sub2proc,lsub2proc,neibproc)

            number_to_receive = sub_aux(isub_loc)%comm_array_numbers(ia)

            if (neibproc /= myid) then
               ! pass messages
               !if (number_to_send.gt.0) then
               !   ireq = ireq + 1
               !   call pp_get_unique_tag(isub,isubadj,comm_all,sub2proc,lsub2proc, tag)
               !   call MPI_ISEND(sub_aux(isub_loc)%comm_array_in, number_to_send, &
               !                  MPI_INTEGER,neibproc,tag, comm_all,request(ireq),ierr)
               !end if
               if (number_to_receive.gt.0) then
                  call pp_get_unique_tag(isubadj,isub,comm_all,sub2proc,lsub2proc, tag)
                  ireq = ireq + 1
                  call MPI_IRECV(sub_aux(isub_loc)%comm_array_out(kineibarray),number_to_receive,&
                                 MPI_INTEGER,neibproc,tag, comm_all,request(ireq),ierr)
               end if
            else
               ! simple copy in memory
               ! find index of the subdomain
               call get_index(isubadj,indexsub,lindexsub,isubadj_loc)
               if (isubadj_loc .eq. -1) then
                  write(*,*) 'indexsub',indexsub
                  call error(routine_name, 'Index of subdomain not found.',isubadj)
               end if
               if (number_to_receive /= size(sub_aux(isubadj_loc)%comm_array_in)) then
                  write (*,*) number_to_receive, size(sub_aux(isubadj_loc)%comm_array_in)
                  write (*,*) 'subdomains in play:', isub, isubadj
                  write (*,*) 'local index of adjacent subdomain:', isubadj_loc
                  call error(routine_name, 'Array size mismatch for subdomain:', isub)
               end if
               sub_aux(isub_loc)%comm_array_out(kineibarray:kineibarray+number_to_receive-1) = sub_aux(isubadj_loc)%comm_array_in
            end if

            kineibarray = kineibarray + number_to_receive
         end do
      end do
      do isub_loc = 1,lindexsub
         isub = indexsub(isub_loc)

         nadj = suba(isub_loc)%nadj

         number_to_send = size(sub_aux(isub_loc)%comm_array_in)

         ! send/receive data from all neighbours
         do ia = 1,nadj

            ! get index of neighbour
            isubadj = suba(isub_loc)%iadj(ia)
            call pp_get_proc_for_sub(isubadj,comm_all,sub2proc,lsub2proc,neibproc)

            if (neibproc /= myid) then
               ! pass message
               if (number_to_send.gt.0) then
                  call pp_get_unique_tag(isub,isubadj,comm_all,sub2proc,lsub2proc, tag)
                  ireq = ireq + 1
                  call MPI_ISEND(sub_aux(isub_loc)%comm_array_in, number_to_send, &
                                 MPI_INTEGER,neibproc,tag, comm_all,request(ireq),ierr)
               end if
            end if

         end do
      end do
      nreq = ireq
      if (nreq.gt.0) then
         call MPI_WAITALL(nreq, request, statarray, ierr)
      end if

      !print *, 'All messages arrived  in the routine.'
      !call flush(6)

      !do isub_loc = 1,lindexsub
      !   deallocate(sub_aux(isub_loc)%comm_array_numbers)
      !end do

      deallocate(request)
      deallocate(statarray)

end subroutine



!**************************************************************************
subroutine dd_create_pairs(nsub, suba,lsuba, indexsub,lindexsub, comm_all,&
                           pairs,lpairs1,lpairs2, npair)
!**************************************************************************
! Subroutine for finding pairs for adaptivity
! based on subdomain data
      use module_utils
      implicit none
      include "mpif.h"

! number of subdomains
      integer, intent(in) :: nsub

! array of sub structure for actual subdomains
      integer,intent(in) ::                lsuba
      type(subdomain_type),intent(inout) :: suba(lsuba)

      integer,intent(in) :: lindexsub
      integer,intent(in) ::  indexsub(lindexsub)
      integer,intent(in) :: comm_all ! MPI communicator
      ! list of pairs for adaptivity
      integer,intent(in) :: lpairs1,lpairs2
      integer,intent(out) :: pairs(lpairs1,lpairs2)
      ! actual number of pairs (number of really used rows in pairs - may be lower that dimension of pairs)
      integer,intent(out) :: npair

      ! local vars
      character(*),parameter:: routine_name = 'DD_CREATE_PAIRS'
      integer :: isub_loc, isub, jsub
      integer :: ncnodes
      integer :: iproc
      integer :: nface_loc
      integer :: npairp, npairp_proc, kpair
      integer ::            lnpairpa
      integer,allocatable :: npairpa(:)
      integer :: npairdoubled
      integer :: i, icnodes, ipairp
      integer :: ign, jgn, ipair, idp

      ! MPI related variables 
      integer :: ierr, nproc, myid
      integer :: stat(MPI_STATUS_SIZE)

      integer ::             lglobal_cnode_number_loc
      integer,allocatable ::  global_cnode_number_loc(:)
      integer ::             lglobal_cnode_number
      integer,allocatable ::  global_cnode_number(:)
      integer ::             lpair_subdomain_loc
      integer,allocatable ::  pair_subdomain_loc(:)
      integer ::             lpair_subdomain
      integer,allocatable ::  pair_subdomain(:)

      integer ::             lsubdomain_costs
      integer,allocatable ::  subdomain_costs_aux(:)
      integer,allocatable ::  subdomain_costs(:)

! orient in communicator
      call MPI_COMM_RANK(comm_all,myid,ierr)
      call MPI_COMM_SIZE(comm_all,nproc,ierr)

      ! find maximal size of local pairs
      nface_loc = 0
      do isub_loc = 1,lindexsub
         ! check prerequisites
         if (.not. suba(isub_loc)%is_cnodes_loaded) then
            call error(routine_name,'Coarse nodes not ready for subdomain',isub)
         end if

         nface_loc = nface_loc + suba(isub_loc)%ncnodes
      end do

      ! prepare local array for finding pairs
      lglobal_cnode_number_loc = nface_loc
      allocate(global_cnode_number_loc(lglobal_cnode_number_loc))
      lpair_subdomain_loc      = nface_loc
      allocate(pair_subdomain_loc(lpair_subdomain_loc))

      ipairp = 0
      do isub_loc = 1,lindexsub
         isub = indexsub(isub_loc)


         ! load data
         ncnodes   = suba(isub_loc)%ncnodes

         ! loop over coarse nodes
         do icnodes = 1,ncnodes
            if (suba(isub_loc)%cnodes(icnodes)%itype.eq.1) then
               ! submit this face as a pair
               ipairp = ipairp + 1

               ! check bounds
               if (ipairp.gt.nface_loc) then
                  call error(routine_name,'not enough space for local pairs')
               end if

               global_cnode_number_loc(ipairp) = suba(isub_loc)%cnodes(icnodes)%global_cnode_number
               pair_subdomain_loc(ipairp)      = isub
            end if
         end do
      end do
      npairp = ipairp

      ! Create global array of pairs
      lnpairpa = nproc
      allocate(npairpa(lnpairpa))
!*****************************************************************MPI
      call MPI_ALLGATHER(npairp,1, MPI_INTEGER,npairpa,1, MPI_INTEGER,  comm_all, ierr) 
!*****************************************************************MPI

      npairdoubled = sum(npairpa)
      ! check that it is even number
      if (mod(npairdoubled,2).ne.0) then
         call error(routine_name,'Twice the number of pairs should be even number.')
      end if
      npair = npairdoubled / 2


      ! determine the global size of array of pair subdomains
      lglobal_cnode_number = npairdoubled
      allocate(global_cnode_number(lglobal_cnode_number))
      lpair_subdomain = npairdoubled
      allocate(pair_subdomain(lpair_subdomain))

      lsubdomain_costs = nsub
      allocate( subdomain_costs_aux(lsubdomain_costs), subdomain_costs(lsubdomain_costs) )
      subdomain_costs_aux = 0

      do isub_loc = 1,lindexsub
         isub = indexsub(isub_loc)

         subdomain_costs_aux(isub) = suba(isub_loc)%mumps_interior_block_factor_size
      end do
!*****************************************************************MPI
      call MPI_ALLREDUCE(subdomain_costs_aux,subdomain_costs,lsubdomain_costs, &
                         MPI_INTEGER, MPI_MAX, comm_all, ierr) 
!*****************************************************************MPI
      deallocate(subdomain_costs_aux)


      if (myid.eq.0) then
         ! first copy my own inodc_loc into global array
         kpair = 0
         do i = 1,npairp
            global_cnode_number(kpair + i) = global_cnode_number_loc(i)
            pair_subdomain(kpair + i)      = pair_subdomain_loc(i)
         end do
         kpair = kpair + npairp
            
         ! then receive data from others
         do iproc = 1,nproc-1
            npairp_proc = npairpa(iproc+1)

            if (npairp_proc .gt.0) then
               call MPI_RECV(global_cnode_number(kpair+1),npairp_proc,MPI_INTEGER,iproc,iproc,comm_all,stat,ierr)
               call MPI_RECV(pair_subdomain(kpair+1),npairp_proc,MPI_INTEGER,iproc,iproc,comm_all,stat,ierr)
            end if

            kpair = kpair + npairp_proc
         end do

      else
         ! send my inodc_loc to root
         if (npairp.gt.0) then
            call MPI_SEND(global_cnode_number_loc,npairp,MPI_INTEGER,0,myid,comm_all,ierr)
            call MPI_SEND(pair_subdomain_loc,npairp,MPI_INTEGER,0,myid,comm_all,ierr)
         end if
      end if

      deallocate(npairpa)
      deallocate(global_cnode_number_loc)
      deallocate(pair_subdomain_loc)

      ! now root sorts corners and remove duplicities
      if (myid.eq.0) then

         ! sort local array of corners before its send
         call iquick_sort_simultaneous(global_cnode_number,lglobal_cnode_number,pair_subdomain,lpair_subdomain)
      end if
      ! root already has global array of pairs, populate it along communicator
!*****************************************************************MPI
      call MPI_BCAST(global_cnode_number,lglobal_cnode_number, MPI_INTEGER, 0, comm_all, ierr) 
      call MPI_BCAST(pair_subdomain,lpair_subdomain, MPI_INTEGER, 0, comm_all, ierr) 
!*****************************************************************MPI
      !print *,'myid =',myid,'ncorner',ncorner
      !call flush(6)

      ! now parse array of pairs to get global pairs

      idp = 0
      do ipair = 1,npair
         idp = idp + 1
         isub = pair_subdomain(idp)
         ign  = global_cnode_number(idp)

         idp = idp + 1
         jsub = pair_subdomain(idp)
         jgn  = global_cnode_number(idp)

         ! check we have pair
         if (ign.ne.jgn) then
            call error(routine_name,'global indices not in pairs',ign)
         end if
         ! mark pair's global cnode number
         pairs(ipair,1) = ign
         ! mark pair's first sharing subdomain
         pairs(ipair,2) = isub
         ! mark pair's second sharing subdomain
         pairs(ipair,3) = jsub
         ! mark pair's first subdomain factorization cost
         pairs(ipair,4) = subdomain_costs(isub)
         ! mark pair's second subdomain factorization cost
         pairs(ipair,5) = subdomain_costs(jsub)
      end do

      deallocate(global_cnode_number)
      deallocate(pair_subdomain)
      deallocate(subdomain_costs)

      ! root print the summary of selection of globs
      if (myid.eq.0) then
         call info ( routine_name, 'Number of pairs for adaptivity:', npair )
      end if

end subroutine

!*****************************************************************
subroutine dd_embed_cnodes(suba,lsuba, indexsub,lindexsub, comm_all,&
                           nndfc,lnndfc)
!*****************************************************************
! Subroutine for embedding of local coarse nodes into global array while creating actual number of degrees of freedom at global nodes 
! based on subdomain data
      use module_utils
      implicit none
      include "mpif.h"

! array of sub structure for actual subdomains
      integer,intent(in) ::                lsuba
      type(subdomain_type),intent(inout) :: suba(lsuba)

      integer,intent(in) :: lindexsub
      integer,intent(in) ::  indexsub(lindexsub)
      integer,intent(in) :: comm_all ! MPI communicator
      ! number of DOF at nodes
      integer,intent(in) :: lnndfc
      integer,intent(out) :: nndfc(lnndfc)

      ! local vars
      character(*),parameter:: routine_name = 'DD_GET_NNDFC'
      integer :: isub_loc, isub
      integer :: ncnodes
      integer :: ncnodesw
      integer :: iproc
      integer :: ncnodes_loc, nvalues
      integer :: ncnodesp, ncnodesp_proc, kcnodes
      integer ::            lncnodespa
      integer,allocatable :: ncnodespa(:)
      integer :: i, icnodes, icnodesp

      ! MPI related variables 
      integer :: ierr, nproc, myid
      integer :: stat(MPI_STATUS_SIZE)

      integer ::             lglobal_cnode_number_loc
      integer,allocatable ::  global_cnode_number_loc(:)
      integer ::             lglobal_cnode_numberw
      integer,allocatable ::  global_cnode_numberw(:)
      integer ::             lnndfc_loc
      integer,allocatable ::  nndfc_loc(:)
      integer ::             lnndfcw
      integer,allocatable ::  nndfcw(:)

      integer ::             lkdofc
      integer,allocatable ::  kdofc(:)
      integer :: icnode, igcnode, indn, kcdof, ncdof

! orient in communicator
      call MPI_COMM_RANK(comm_all,myid,ierr)
      call MPI_COMM_SIZE(comm_all,nproc,ierr)

      ! find maximal size of local pairs
      ncnodes_loc = 0
      do isub_loc = 1,lindexsub
         ! check prerequisites
         if (.not. suba(isub_loc)%is_cnodes_loaded) then
            call error(routine_name,'Coarse nodes not ready for subdomain',isub)
         end if

         ncnodes_loc = ncnodes_loc + suba(isub_loc)%ncnodes
      end do

      ! prepare local array for finding pairs
      lglobal_cnode_number_loc = ncnodes_loc
      allocate(global_cnode_number_loc(lglobal_cnode_number_loc))
      lnndfc_loc               = ncnodes_loc
      allocate(nndfc_loc(lnndfc_loc))

      icnodesp = 0
      do isub_loc = 1,lindexsub
         isub = indexsub(isub_loc)

         ! load data
         ncnodes   = suba(isub_loc)%ncnodes

         ! loop over coarse nodes
         do icnodes = 1,ncnodes
            icnodesp = icnodesp + 1

            ! check bounds
            if (icnodesp.gt.ncnodes_loc) then
               call error(routine_name,'not enough space for local nodes')
            end if

            global_cnode_number_loc(icnodesp) = suba(isub_loc)%cnodes(icnodes)%global_cnode_number
            nndfc_loc(icnodesp)               = suba(isub_loc)%cnodes(icnodes)%ncdof
         end do
      end do
      ncnodesp = icnodesp

      ! Create global array of pairs
      lncnodespa = nproc
      allocate(ncnodespa(lncnodespa))
!*****************************************************************MPI
      call MPI_ALLGATHER(ncnodesp,1, MPI_INTEGER,ncnodespa,1, MPI_INTEGER, comm_all, ierr) 
!*****************************************************************MPI

      ncnodesw = sum(ncnodespa)

      if (myid.eq.0) then

         ! determine the global size of array of pair subdomains
         lglobal_cnode_numberw = ncnodesw
         allocate(global_cnode_numberw(lglobal_cnode_numberw))
         lnndfcw = ncnodesw
         allocate(nndfcw(lnndfcw))

         ! first copy my own nndfc_loc into global array
         kcnodes = 0
         do i = 1,ncnodesp
            global_cnode_numberw(kcnodes + i) = global_cnode_number_loc(i)
            nndfcw(kcnodes + i)               = nndfc_loc(i)
         end do
         kcnodes = kcnodes + ncnodesp
            
         ! then receive data from others
         do iproc = 1,nproc-1
            ncnodesp_proc = ncnodespa(iproc+1)

            if (ncnodesp_proc .gt.0) then
               call MPI_RECV(global_cnode_numberw(kcnodes+1),ncnodesp_proc,MPI_INTEGER,iproc,iproc,comm_all,stat,ierr)
               call MPI_RECV(nndfcw(kcnodes+1),ncnodesp_proc,MPI_INTEGER,iproc,iproc,comm_all,stat,ierr)
            end if

            kcnodes = kcnodes + ncnodesp_proc
         end do

      else
         ! send my arrays to root
         if (ncnodesp.gt.0) then
            call MPI_SEND(global_cnode_number_loc,ncnodesp,MPI_INTEGER,0,myid,comm_all,ierr)
            call MPI_SEND(nndfc_loc,ncnodesp,MPI_INTEGER,0,myid,comm_all,ierr)
         end if
      end if

      deallocate(ncnodespa)
      deallocate(global_cnode_number_loc)
      deallocate(nndfc_loc)


      ! now root sorts corners and remove duplicities
      if (myid.eq.0) then

         ! sort local array of corners before its send
         call iquick_sort_simultaneous(global_cnode_numberw,lglobal_cnode_numberw,nndfcw,lnndfcw)
         ! check what I got
         do i = 1,ncnodesw-1
            if (global_cnode_numberw(i).eq.global_cnode_numberw(i+1)) then
               if (nndfcw(i) .ne. nndfcw(i+1)) then
                  call error(routine_name,'different number of dof at a coarse node from different subdomains')
               end if
            end if
         end do
         call get_array_norepeat_simultaneous(global_cnode_numberw,lglobal_cnode_numberw,nndfcw,lnndfcw,nvalues)

         ! check what I got
         if (nvalues.ne.lnndfc) then
            call error(routine_name,'number of global coarse nodes does not match')
         end if
         do i = 1,nvalues
            if (global_cnode_numberw(i) .ne. i) then
               call error(routine_name,'there is mismatch in global indices of coarse nodes')
            end if
         end do

         ! if checks are OK, copy resulting array
         do i = 1,nvalues
            nndfc(i) = nndfcw(i)
         end do

         deallocate(global_cnode_numberw)
         deallocate(nndfcw)
      end if
      ! root already has global array of pairs, populate it along communicator
!*****************************************************************MPI
      call MPI_BCAST(nndfc,lnndfc, MPI_INTEGER, 0, comm_all, ierr) 
!*****************************************************************MPI
      !print *,'myid =',myid,'nndfc',nndfc
      !call flush(6)


      ! create array of global coarse dof KDOFC(ncorner) with addresses before first global dof
      lkdofc = lnndfc
      allocate(kdofc(lkdofc))
      if ( lkdofc .gt. 0 ) then
         kdofc(1) = 0
      end if
      do indn = 2,lnndfc
         kdofc(indn) = kdofc(indn-1) + nndfc(indn-1)
      end do

      ! Update embedding arrays
      do isub_loc = 1,lindexsub
         isub = indexsub(isub_loc)

         ncnodes =  suba(isub_loc)%ncnodes

         do icnode = 1,ncnodes
            ! number of coarse dof it contains
            igcnode = suba(isub_loc)%cnodes(icnode)%global_cnode_number 
            ncdof = nndfc(igcnode)
            if (allocated(suba(isub_loc)%cnodes(icnode)%igcdof)) then
               deallocate(suba(isub_loc)%cnodes(icnode)%igcdof)
            end if

            suba(isub_loc)%cnodes(icnode)%ncdof = ncdof
            allocate(suba(isub_loc)%cnodes(icnode)%igcdof(ncdof))
            ! fill coarse node dof
            kcdof = kdofc(igcnode)
            do i = 1,ncdof
               suba(isub_loc)%cnodes(icnode)%igcdof(i) = kcdof + i
            end do

         end do
         suba(isub_loc)%is_cnodes_embedded = .true.
      end do

      deallocate(kdofc)

end subroutine


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
      use module_pp, only : pp_get_proc_for_sub, pp_get_unique_tag
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
      character(*),parameter:: routine_name = 'DD_COMM_SWAPDATA'
      integer :: i, ia, isub, isub_loc, nsub
      integer :: isubadj, isubadj_loc, procadj

      integer :: isubadj_ia
      integer :: kneib_isub
      integer :: kneib_isubadj
      integer :: nadj, nsharedv
      integer :: nsubd

      ! MPI related arrays and variables
      integer :: myid, nproc, ierr, tag, ireq, nreq
      integer ::            lrequest
      integer,allocatable :: request(:)
      integer             :: lstatarray1
      integer             :: lstatarray2
      integer,allocatable :: statarray(:,:)

! orient in communicator
      call MPI_COMM_RANK(comm_all,myid,ierr)
      call MPI_COMM_SIZE(comm_all,nproc,ierr)

      nsub = sub2proc(nproc+1) - 1

      ! get number of requests
      ireq = 0
      do isub_loc = 1,lindexsub
         isub = indexsub(isub_loc)
         nsubd = suba(isub_loc)%nsub

         if (nsubd.ne.nsub) then
            call error(routine_name,'number of subdomains mismatch from DD module and from arguments for sub',isub)
         end if
         if (.not. suba(isub_loc)%is_neighbouring_ready) then
            call error(routine_name,'Neighbouring is not ready for subdomain ',isub)
         end if
         if (.not. allocated(suba(isub_loc)%commvec_out)) then
            call error(routine_name,'Array COMMVEC_OUT not allocated. Perhaps missing call to dd_comm_upload. ',isub)
         end if

         ! load data
         nadj  = suba(isub_loc)%nadj

         ireq = ireq + nadj
      end do

      ! prepare MPI data for processor
      ! double the size for two sided non-blocking communication
      nreq = 2 * ireq
      lrequest = nreq
      allocate(request(lrequest))
      lstatarray1 = MPI_STATUS_SIZE
      lstatarray2 = nreq
      allocate(statarray(lstatarray1,lstatarray2))

!     Interchange interface nodes in global numbering
      ireq = 0
      do isub_loc = 1,lindexsub
         isub = indexsub(isub_loc)

         nadj = suba(isub_loc)%nadj

         do ia = 1,nadj
         ! get index of neighbour
            isubadj = suba(isub_loc)%iadj(ia)

            nsharedv = suba(isub_loc)%kshvadj(ia + 1) - suba(isub_loc)%kshvadj(ia)
            if (nsharedv.gt.0) then

               ! who owns this subdomain?
               call pp_get_proc_for_sub(isubadj,comm_all,sub2proc,lsub2proc,procadj)

               if (procadj .ne. myid) then
                  ! interchange via MPI
   
                  ! raise a request to receive data
                  call pp_get_unique_tag(isubadj,isub,comm_all,sub2proc,lsub2proc,tag)
                  ireq = ireq + 1
                  call MPI_IRECV(suba(isub_loc)%commvec_in(suba(isub_loc)%kshvadj(ia)),nsharedv,MPI_DOUBLE_PRECISION,&
                                 procadj,tag,comm_all,request(ireq),ierr)
               else 
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
               end if
            end if
   
         end do
      end do
      do isub_loc = 1,lindexsub
         isub = indexsub(isub_loc)

         nadj = suba(isub_loc)%nadj

!     Interchange interface nodes in global numbering
         do ia = 1,nadj
         ! get index of neighbour
            isubadj = suba(isub_loc)%iadj(ia)

            nsharedv = suba(isub_loc)%kshvadj(ia + 1) - suba(isub_loc)%kshvadj(ia)
            if (nsharedv.gt.0) then

               ! who owns this subdomain?
               call pp_get_proc_for_sub(isubadj,comm_all,sub2proc,lsub2proc,procadj)

               if (procadj .ne. myid) then
                  ! interchange via MPI
   
                  ! send him my data
                  call pp_get_unique_tag(isub,isubadj,comm_all,sub2proc,lsub2proc,tag)
                  ireq = ireq + 1
                  call MPI_ISEND(suba(isub_loc)%commvec_out(suba(isub_loc)%kshvadj(ia)),nsharedv,MPI_DOUBLE_PRECISION,&
                                 procadj,tag,comm_all,request(ireq),ierr)
               end if
            end if
   
         end do
      end do
      nreq = ireq
      if (nreq.gt.0) then
         call MPI_WAITALL(nreq, request, statarray, ierr)
      end if

      ! clear memory
      deallocate(request)
      deallocate(statarray)

end subroutine

!*****************************************
subroutine dd_comm_download(sub, vec,lvec)
!*****************************************
! Subroutine that downloads data from communication and add it to subdomain vector VEC
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
      ! sum up repeated entries
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

!*******************************************************************************************************
subroutine dd_weights_prepare(suba,lsuba, sub2proc,lsub2proc,indexsub,lindexsub, comm_all, weights_type)
!*******************************************************************************************************
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
      ! type of weights:
      ! 0 - weights by cardinality
      ! 1 - weights by diagonal stiffness
      ! 2 - weights based on first row of element data
      ! 3 - weights based on dof data
      ! 4 - weights by Marta Certikova - unit load
      ! 5 - weights by Marta Certikova - unit jump
      ! 6 - weights by Schur row sums for whole subdomain
      ! 7 - weights by Schur row sums computed face by face
      ! 8 - weights by Schur diagonal - available for explicit Schur complements
      integer,intent(in) :: weights_type 

      ! local vars
      character(*),parameter:: routine_name = 'DD_WEIGHTS_PREPARE'
      integer :: isub_loc, i
      integer :: ndofi

      integer ::             lrhoi
      real(kr),allocatable :: rhoi(:)
      real(kr),allocatable :: rhoiaux(:)
      integer ::             lwi
      real(kr),allocatable :: wi(:)

      integer :: myid, ierr

      call MPI_COMM_RANK(comm_all,myid,ierr)
      if (myid.eq.0) then
          call info(routine_name,'Using weights type: ',weights_type)
      end if

      ! Prepare data for communication
      do isub_loc = 1,lindexsub
         ndofi = suba(isub_loc)%ndofi

         lrhoi = ndofi
         allocate(rhoi(lrhoi))
         call dd_get_my_coefficients_for_weights(suba(isub_loc), weights_type, rhoi,lrhoi)
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

         rhoiaux(:) = 0._kr
         call dd_comm_download(suba(isub_loc), rhoiaux,lrhoi)

         call dd_get_my_coefficients_for_weights(suba(isub_loc), weights_type, rhoi,lrhoi)

         ! check that denominator is nonzero
         if (any(abs(rhoi + rhoiaux) .eq. 0._kr)) then
             write(*,*) 'rhoi' 
             write(*,*) rhoi 
             write(*,*) 'rhoiaux' 
             write(*,*) rhoiaux 
             call error( routine_name, 'Zero in denominator in computation of weights for subdomain:', suba(isub_loc)%isub )
         end if

         ! compute weight
         wi = rhoi / (rhoi + rhoiaux)

         ! load wi into structure
         suba(isub_loc)%lwi = lwi
         allocate(suba(isub_loc)%wi(lwi))
         do i = 1,lwi
            suba(isub_loc)%wi(i) = wi(i)
         end do
         suba(isub_loc)%weights_type     = weights_type
         suba(isub_loc)%is_weights_ready = .true.

         deallocate(wi)
         deallocate(rhoiaux)
         deallocate(rhoi)
      end do

end subroutine

!***************************************************************************
subroutine dd_get_my_coefficients_for_weights(sub, weights_type, rhoi,lrhoi)
!***************************************************************************
! Subroutine for counting neighbours at dofs at the interface
      use module_utils
      implicit none
! Subdomain structure
      type(subdomain_type),intent(inout) :: sub
! type of weights:
! 0 - weights by cardinality
! 1 - weights by diagonal stiffness
! 2 - weights based on first row of element data
! 3 - weights based on dof data
! 4 - weights by Marta Certikova - unit load
! 5 - weights by Marta Certikova - unit jump
! 6 - weights by Schur row sums for whole subdomain
! 7 - weights by Schur row sums computed face by face
! 8 - weights by Schur diagonal - available for explicit Schur complements
      integer,intent(in) :: weights_type 

! interface vector of my coefficients
      integer,intent(in) ::  lrhoi
      real(kr),intent(out)::  rhoi(lrhoi)

      ! local vars
      character(*),parameter:: routine_name = 'DD_GET_MY_COEFFICIENTS_FOR_WEIGHTS'

      ! check input data
      if (lrhoi .ne. sub%ndofi) then
         call error(routine_name, 'Array lenght mismatch for subdomain',sub%isub)
      end if

      if (weights_type .eq. 0) then
         rhoi = 1._kr
      else if (weights_type .eq. 1) then
         call dd_get_interface_diagonal(sub, rhoi,lrhoi)
         ! if zeros can be on diagonal, correct it
         if (count(rhoi == 0.) > 0) then
            call warning(routine_name, 'Replacing zeros in weights by diagonal stiffness', &
                         count(rhoi == 0.))
            where (rhoi == 0.) rhoi = 1._kr
         end if
      else if (weights_type .eq. 2) then
         call dd_get_interface_element_data(sub, rhoi,lrhoi)
      else if (weights_type .eq. 3) then
         call dd_get_interface_dof_data(sub, rhoi,lrhoi)
      else if (weights_type .eq. 4) then
         call dd_generate_interface_unit_load(sub, rhoi,lrhoi)
      else if (weights_type .eq. 5) then
         call dd_generate_interface_unit_jump(sub, rhoi,lrhoi)
      else if (weights_type .eq. 6) then
         call dd_generate_interface_schur_row_sums(sub, rhoi,lrhoi)
      else if (weights_type .eq. 7) then
         call dd_generate_interface_face_schur_row_sums(sub, rhoi,lrhoi)
      else if (weights_type .eq. 8) then
         if (sub%is_explicit_schur_prepared) then
            call dd_get_interface_schur_diagonal(sub, rhoi,lrhoi)
         else
            call warning(routine_name, 'Explicit Schur is not available, using stiffness scaling for ', sub%isub)
            call dd_get_interface_diagonal(sub, rhoi,lrhoi)
         end if
      else
         call error(routine_name,'Type of weights not supported:',weights_type)
      end if

end subroutine

!*************************************************************************
subroutine dd_count_neighbours_of_interface_dof(sub, nneighbiv,lnneighbiv)
!*************************************************************************
! Subroutine for counting neighbours at dofs at the interface
      use module_utils
      implicit none
! Subdomain structure
      type(subdomain_type),intent(in) :: sub
      integer,intent(in) :: lnneighbiv
      integer,intent(out)::  nneighbiv(lnneighbiv)

      ! local vars
      character(*),parameter:: routine_name = 'DD_COUNT_NEIGHBOURS_OF_INTERFACE_DOF'
      integer :: i, isubadj, ia, ishn, idofn, indn, indni
      integer :: nnodi, nshn, ndofn
      integer :: kishnadj

      integer ::            lkdofi
      integer,allocatable :: kdofi(:)

      ! check input data
      if (lnneighbiv .ne. sub%ndofi) then
         call error(routine_name, 'Array lenght mismatch for subdomain',sub%isub)
      end if

      ! zero the array
      nneighbiv = 0

      ! prepare array kdofi
      nnodi = sub%nnodi
      lkdofi = nnodi + 1
      allocate(kdofi(lkdofi))
      kdofi(1) = 1
      do i = 1,nnodi
         indn = sub%iin(i)
         ndofn = sub%nndf(indn)
         
         kdofi(i + 1) = kdofi(i) + ndofn
      end do

      ! loop over neighbours
      kishnadj  = 0
      do ia = 1,sub%nadj
         ! get index of neighbour
         isubadj = sub%iadj(ia)

         ! number of nodes shared with this subdomain
         nshn = sub%nshnadj(ia)

         do ishn = 1,nshn
            ! index of shared node in interface numbering
            indni  = sub%ishnadj(kishnadj + ishn)
            indn   = sub%iin(indni)

            ndofn = sub%nndf(indn)

            do idofn = 1,ndofn
               nneighbiv(kdofi(indni)-1 + idofn) = nneighbiv(kdofi(indni)-1 + idofn) + 1
            end do
         end do
   
         kishnadj = kishnadj + nshn
      end do

      ! check coverage of interface
      if (any(nneighbiv.eq.0)) then
         call error(routine_name,'Coverage of interface failed for subdomain',sub%isub)
      end if

      deallocate(kdofi)

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
      rho(:) = 0._kr

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

!********************************************************
subroutine dd_get_interface_element_data(sub, rhoi,lrhoi)
!********************************************************
! Subroutine for getting subdomain diagonal from the structure
      use module_utils
      use module_graph
      implicit none
! Subdomain structure
      type(subdomain_type),intent(in) :: sub
      integer,intent(in)  :: lrhoi
      real(kr),intent(out) ::  rhoi(lrhoi)

      ! local vars
      character(*),parameter:: routine_name = 'DD_GET_INTERFACE_ELEMENT_DATA'
      integer :: nnod, nnodi, linet 
      integer :: i, ind, nen, ien, ndofn, idofn, inde, indn, pointen

      real(kr) :: average_coef, sum_coef, coef

      integer             :: lnetn,  lietn,  lkietn
      integer, allocatable :: netn(:),ietn(:),kietn(:)

      ! check if mesh is loaded
      if (.not. sub%is_mesh_loaded) then
         call error(routine_name, 'Mesh not loaded for subdomain: ',sub%isub)
      end if
      ! check if matrix is loaded
      if (.not. sub%is_element_data_loaded) then
         call error(routine_name, 'Element data not loaded for subdomain: ',sub%isub)
      end if
      ! check dimensions
      if (sub%ndofi.ne.lrhoi) then
         call error(routine_name, 'Interface dimensions mismatch for subdomain ',sub%isub)
      end if
      if (sub%lelement_data1 .lt. 1) then
         call error(routine_name, 'At least one row of element data needed for subdomain',sub%isub)
      end if
      if (sub%nelem.ne.sub%lelement_data2) then
         call error(routine_name, 'Dimensions mismatch for subdomain ',sub%isub)
      end if

      ! load data
      nnod  = sub%nnod
      nnodi = sub%nnodi
      linet = sub%linet

      ! prepare dual mesh
      lnetn  = nnod
      lietn  = linet
      lkietn = nnod
      allocate(netn(lnetn),ietn(lietn),kietn(lkietn))
      call graph_get_dual_mesh(sub%nelem,sub%nnod,&
                               sub%inet,sub%linet,sub%nnet,sub%lnnet,&
                               netn,lnetn,ietn,lietn,kietn,lkietn)

      ind = 0
      do i = 1,nnodi
         indn = sub%iin(i)
         ndofn = sub%nndf(indn)

         pointen = kietn(indn)
         nen     = netn(indn)
         sum_coef = 0._kr 
         do ien = 1,nen
            inde = ietn(pointen + ien)
            coef = sub%element_data(1,inde)

            sum_coef = sum_coef + coef
         end do
         average_coef = sum_coef / nen

         do idofn = 1,ndofn
            ind = ind + 1
            if (ind .gt. lrhoi) then
               call error(routine_name, 'Index out of bounds for RHOI', ind)
            end if
            rhoi(ind) = average_coef
         end do
      end do
      !write(*,*) 'rhoi', rhoi

      deallocate(netn,ietn,kietn)
end subroutine

!****************************************************
subroutine dd_get_interface_dof_data(sub, rhoi,lrhoi)
!****************************************************
! Subroutine for getting weights from data in subdomain dof_data 
      use module_utils
      implicit none
! Subdomain structure
      type(subdomain_type),intent(in) :: sub
      integer,intent(in)  :: lrhoi
      real(kr),intent(out) ::  rhoi(lrhoi)

      ! local vars
      character(*),parameter:: routine_name = 'DD_GET_INTERFACE_DOF_DATA'
      integer :: i, indv, ndofi

      ! check if mesh is loaded
      if (.not. sub%is_mesh_loaded) then
         call error(routine_name, 'Mesh not loaded for subdomain: ',sub%isub)
      end if
      ! check if dof data are loaded
      if (.not. sub%is_dof_data_loaded) then
          call error(routine_name, 'DOF data not loaded for subdomain: ',sub%isub)
      end if
      ! check dimensions
      if (sub%ndofi.ne.lrhoi) then
         call error(routine_name, 'Interface dimensions mismatch for subdomain ',sub%isub)
      end if
      if (sub%ndof.ne.sub%ldof_data) then
          call error(routine_name, 'DOF data dimensions mismatch for subdomain ',sub%isub)
      end if

      ! load data
      ndofi = sub%ndofi
      do i = 1,ndofi
         indv = sub%iivsvn(i)

         if (indv .gt. sub%ldof_data) then
            call error(routine_name, 'Index out of bounds for DOF_DATA', indv)
         end if
         rhoi(i) = sub%dof_data(indv)
      end do
      !write(*,*) 'rhoi', rhoi
end subroutine

!******************************************************
subroutine dd_generate_interface_unit_load(sub, vi,lvi)
!******************************************************
! Subroutine for generating weights as a solution to local problems with unit
! load.
! Based on the paper
! M. Certikova, J. Sistek, P. Burda: On selection of interface weights in domain decomposition methods,
! Proceedings of Programs and Algorithms of Numerical Mathematics 16, Institute
! of Mathematics AS CR, pp. 35-44, 2013
      use module_utils
      implicit none
! Subdomain structure
      type(subdomain_type),intent(inout) :: sub
      integer,intent(in)  :: lvi
      real(kr),intent(out) :: vi(lvi)

      ! local vars
      character(*),parameter:: routine_name = 'DD_GENERATE_INTERFACE_UNIT_LOAD'
      integer :: i, j, ndof, nnodi, ndofi, ndofaug, ndofn, nglob, nnadj, ncorner
      integer :: ndofnx 
      integer :: lri1, lr1, lvaug1
      integer :: lri2, lr2, lvaug2
      integer :: nrhs
      logical :: solve_adjoint

      integer :: isubadj, jsub, ia, idofn, iglob, inadj, indn, indshni, kishnadj, nadj, icorner, &
                 indi, ind, inodi

      real(kr), allocatable :: ri(:,:)
      real(kr), allocatable :: r(:,:)
      real(kr), allocatable :: vaug(:,:)

      integer ::            lkdofi
      integer,allocatable :: kdofi(:)

      ! check if mesh is loaded
      if (.not. sub%is_mesh_loaded) then
         call error(routine_name, 'Mesh not loaded for subdomain: ',sub%isub)
      end if
      ! check if dof data are loaded
      if (.not. sub%is_coarse_prepared) then
         call error(routine_name, 'Coarse problem is not factorized for subdomain: ',sub%isub)
      end if
      ! check dimensions
      if (sub%ndofi.ne.lvi) then
         call error(routine_name, 'Interface dimensions mismatch for subdomain ',sub%isub)
      end if

      ! prepare array kdofi
      nnodi = sub%nnodi
      lkdofi = nnodi + 1
      allocate(kdofi(lkdofi))
      kdofi(1) = 1
      do i = 1,nnodi
         indn = sub%iin(i)
         ndofn = sub%nndf(indn)
         
         kdofi(i + 1) = kdofi(i) + ndofn
      end do

      ndofnx = minval(sub%nndf)
      !if (ndofnx .eq. 0) then
      !   call error(routine_name, 'Some nodes have zero degrees of freedom for subdomain ',sub%isub)
      !end if

      ! prepare unit load vector at subdomain interface
      ndofi = sub%ndofi
      lri1 = ndofi
      lri2 = ndofnx
      allocate(ri(lri1,lri2))

      ndof = sub%ndof
      lr1 = ndof
      lr2 = ndofnx
      allocate(r(lr1,lr2))

      ! prepare array of augmented size
      call dd_get_aug_size(sub, ndofaug)
      lvaug1 = ndofaug
      lvaug2 = ndofnx
      allocate(vaug(lvaug1,lvaug2))

      vi(:) = 0._kr

      ! only continue for non-degenerate subdomains
      if (ndofnx .gt. 0 ) then

         ! loop over globs and identify faces
         nglob = sub%nglob
         nadj  = sub%nadj
         do iglob = 1,nglob

            ! only select faces 
            if (sub%glob_type(iglob) .eq. 1) then

               ! the second subdomain
               if (sub%nsubglobs(iglob).lt.1) then
                  call error(routine_name,'Face appears to have no neighbour for subdomain:',sub%isub)
               end if
               if (sub%lglob_subs2.lt.1) then
                  call error(routine_name,'Wrong dimension of array glob_subs for subdomain:',sub%isub)
               end if
               jsub = sub%glob_subs(iglob,1)

               ! find all nodes shared with this subdomain
               kishnadj = 0
               do ia = 1,nadj
                  nnadj = sub%nshnadj(ia)

                  ! get index of neighbour
                  isubadj = sub%iadj(ia)
                  ! perform the solves only if the subdomain matches the face neighbour
                  if (isubadj .eq. jsub) then

                     ! fill in ones in selected direction
                     ri = 0._kr
                     do inadj = 1,nnadj
                        indshni = sub%ishnadj(kishnadj + inadj)

                        do idofn = 1,ndofnx
                           ri(kdofi(indshni)-1 + idofn,idofn) = 1._kr
                        end do
                     end do

                     !print *, 'ri before'
                     !do i = 1,lri1
                     !   write (*,'(i2,2x,10e7.1)') i, ri(i,:)
                     !end do

                     r = 0._kr
                     do idofn = 1,ndofnx
                        call dd_map_subi_to_sub(sub, ri(1,idofn),lri1, r(1,idofn),lr1)
                     end do

                     ! subdomain data are the first part of the augmented vector
                     vaug = 0._kr
                     do j = 1,ndofnx
                        do i = 1,ndof
                           vaug(i,j) = r(i,j)
                        end do
                     end do

                     ! solve augmented problem
                     nrhs = ndofnx
                     solve_adjoint = .false.
                     call dd_solve_aug(sub, vaug,lvaug1*lvaug2, nrhs, solve_adjoint)

                     ! extract interface values
                     ri = 0._kr
                     do idofn = 1,ndofnx
                        call dd_map_sub_to_subi(sub, vaug(1,idofn),ndof, ri(1,idofn),lri1)
                     end do

                     ! copy result to the vector of weights
                     do inadj = 1,nnadj
                        indshni = sub%ishnadj(kishnadj + inadj)

                        do idofn = 1,ndofnx
                           vi(kdofi(indshni)-1 + idofn) = vi(kdofi(indshni)-1 + idofn) + ri(kdofi(indshni)-1 + idofn,idofn)
                        end do
                     end do

                     !print *, 'ri after'
                     !do i = 1,lri1
                     !   write (*,'(i2,2x,10e7.1)') i, ri(i,:), vi(i)
                     !end do

                  end if
                  kishnadj = kishnadj + nnadj
               end do
            end if
         end do
      end if

      ! add ones to corners
      ncorner = sub%ncorner  
      do icorner = 1,ncorner
         indi = sub%icnsin(icorner) 

         do idofn = 1,ndofnx
            vi(kdofi(indi)-1+idofn) = 1._kr
         end do
      end do

      !print *, 'vi with corners'
      !do i = 1,lvi
      !   write (*,'(i2,2x,10e7.1)') i, vi(i)
      !end do

      ! add ones to dofs above the ndofnx
      do inodi = 1,nnodi
         ind = sub%iin(inodi)
         ndofn = sub%nndf(ind)
         do idofn = ndofnx+1,ndofn
            vi(kdofi(inodi)-1+idofn) = 1._kr
         end do
      end do

      ! other weights may be related to average constraints
      !if (any(vi.le.0._kr)) then
      !   call warning(routine_name,'zeros in weights for subdomain',sub%isub)
      !end if
      where (vi.eq.0._kr) vi = 1._kr

      ! avoid negative weights
      vi = abs(vi)

      ! make results compatible with scheme of their construction 
      vi = 1._kr / vi

      deallocate(ri)
      deallocate(r)
      deallocate(vaug)
      deallocate(kdofi)
end subroutine

!***********************************************************
subroutine dd_generate_interface_schur_row_sums(sub, vi,lvi)
!***********************************************************
! Subroutine for generating weights as a solution to local problems with unit jump, but using 
! different value at each dof - w_i = S*ones(ndofi,1)
      use module_utils
      implicit none
! Subdomain structure
      type(subdomain_type),intent(inout) :: sub
      integer,intent(in)  :: lvi
      real(kr),intent(out) :: vi(lvi)

      ! local vars
      character(*),parameter:: routine_name = 'DD_GENERATE_INTERFACE_SCHUR_ROW_SUMS'
      integer :: ndofi
      integer :: lri
      real(kr), allocatable :: ri(:)
      integer :: ncol

      ! check if mesh is loaded
      if (.not. sub%is_mesh_loaded) then
         call error(routine_name, 'Mesh not loaded for subdomain: ',sub%isub)
      end if
      ! check if dof data are loaded
      if (.not. sub%is_interior_factorized) then
         call error(routine_name, 'Interior problem is not factorized for subdomain: ',sub%isub)
      end if
      ! check dimensions
      if (sub%ndofi.ne.lvi) then
         call error(routine_name, 'Interface dimensions mismatch for subdomain ',sub%isub)
      end if


      ! prepare unit load vector at subdomain interface
      ndofi = sub%ndofi
      vi(:) = 0._kr

      lri = ndofi
      allocate(ri(lri))
      ri = 1._kr

      ! multiply the vectors by Schur complement
      ncol = 1
      call dd_multiply_by_schur(sub, ri,lri, vi,lvi, ncol)

      ! avoid negative and zero weights
      vi = abs(vi)
      where (vi.lt.numerical_zero) vi = numerical_zero

      deallocate(ri)

end subroutine

!****************************************************************
subroutine dd_generate_interface_face_schur_row_sums(sub, vi,lvi)
!****************************************************************
! Subroutine for generating weights as a solution to local problems with unit jump, but using 
! different value at each dof - w_i = S*ones(ndofi,1)
! Compared to dd_generate_interface_schur_row_sums, it is constructed locally face by face

      use module_utils
      implicit none
! Subdomain structure
      type(subdomain_type),intent(inout) :: sub
      integer,intent(in)  :: lvi
      real(kr),intent(out) :: vi(lvi)

      ! local vars
      character(*),parameter:: routine_name = 'DD_GENERATE_INTERFACE_FACE_SCHUR_ROW_SUMS'
      integer :: i, nnodi, ndofi, ndofn, nglob, nnadj
      integer :: ndofnx 
      integer :: lri1, lri2 
      real(kr), allocatable :: ri(:)
      integer :: lsi1, lsi2 
      real(kr), allocatable :: si(:)
      integer :: ncol

      integer :: isubadj, jsub, ia, idofn, iglob, inadj, indn, indshni, kishnadj, nadj, ind, inodi

      integer ::            lkdofi
      integer,allocatable :: kdofi(:)

      ! check if mesh is loaded
      if (.not. sub%is_mesh_loaded) then
         call error(routine_name, 'Mesh not loaded for subdomain: ',sub%isub)
      end if
      ! check if dof data are loaded
      if (.not. sub%is_interior_factorized) then
         call error(routine_name, 'Interior problem is not factorized for subdomain: ',sub%isub)
      end if
      ! check dimensions
      if (sub%ndofi.ne.lvi) then
         call error(routine_name, 'Interface dimensions mismatch for subdomain ',sub%isub)
      end if

      ! prepare array kdofi
      nnodi = sub%nnodi
      lkdofi = nnodi + 1
      allocate(kdofi(lkdofi))
      kdofi(1) = 1
      do i = 1,nnodi
         indn = sub%iin(i)
         ndofn = sub%nndf(indn)
         
         kdofi(i + 1) = kdofi(i) + ndofn
      end do

      ndofnx = minval(sub%nndf)

      ! prepare unit load vector at subdomain interface
      ndofi = sub%ndofi
      lri1 = ndofi
      lri2 = ndofnx
      allocate(ri(lri1*lri2))
      lsi1 = ndofi
      lsi2 = ndofnx
      allocate(si(lsi1*lsi2))

      vi(:) = 0._kr

      if (ndofnx.gt.0) then

         ! loop over globs and identify faces
         nglob = sub%nglob
         nadj  = sub%nadj
         do iglob = 1,nglob

            ! only select faces 
            if (sub%glob_type(iglob) .eq. 1) then

               ! the second subdomain
               if (sub%nsubglobs(iglob).lt.1) then
                  call error(routine_name,'Face appears to have no neighbour for subdomain:',sub%isub)
               end if
               if (sub%lglob_subs2.lt.1) then
                  call error(routine_name,'Wrong dimension of array glob_subs for subdomain:',sub%isub)
               end if
               jsub = sub%glob_subs(iglob,1)

               ! find all nodes shared with this subdomain
               kishnadj = 0
               do ia = 1,nadj
                  nnadj = sub%nshnadj(ia)

                  ! get index of neighbour
                  isubadj = sub%iadj(ia)
                  ! perform the solves only if the subdomain matches the face neighbour
                  if (isubadj .eq. jsub) then

                     ! fill in ones in selected direction
                     ri = 0._kr
                     do inadj = 1,nnadj
                        indshni = sub%ishnadj(kishnadj + inadj)

                        do idofn = 1,ndofnx
                           ri((idofn-1)*lri1 + (kdofi(indshni)-1 + idofn)) = 1._kr
                        end do
                     end do

                     !print *, 'ri before'
                     !do i = 1,lri1
                     !   write (*,'(i2,2x,10e7.1)') i, ri(i,:)
                     !end do

                     ! multiply the vectors by Schur complement
                     ncol = ndofnx
                     call dd_multiply_by_schur(sub, ri,lri1*lri2, si,lsi1*lsi2, ncol)

                     ! fill in ones in selected direction
                     do inadj = 1,nnadj
                        indshni = sub%ishnadj(kishnadj + inadj)

                        do idofn = 1,ndofnx
                           vi(kdofi(indshni)-1 + idofn) = vi(kdofi(indshni)-1 + idofn) &
                                                        + si((idofn-1)*lsi1 + kdofi(indshni)-1 + idofn)
                        end do
                     end do

                     !print *, 'ri after'
                     !do i = 1,lri1
                     !   write (*,'(i2,2x,10e7.1)') i, ri(i,:), vi(i)
                     !end do

                  end if
                  kishnadj = kishnadj + nnadj
               end do
            end if
         end do
      end if

      ! add ones to corners
      !ncorner = sub%ncorner  
      !do icorner = 1,ncorner
      !   indi = sub%icnsin(icorner) 

      !   do idofn = 1,ndofnx
      !      vi(kdofi(indi)-1+idofn) = 1._kr
      !   end do
      !end do

      !print *, 'vi with corners'
      !do i = 1,lvi
      !   write (*,'(i2,2x,10e7.1)') i, vi(i)
      !end do

      ! add ones to dofs above the ndofnx
      do inodi = 1,nnodi
         ind = sub%iin(inodi)
         ndofn = sub%nndf(ind)
         do idofn = ndofnx+1,ndofn
            vi(kdofi(inodi)-1+idofn) = 1._kr
         end do
      end do

      !if (any(vi.le.0._kr)) then
      !   call warning(routine_name,'zeros in weights for subdomain',sub%isub)
      !   print *, 'vi', vi
      !end if

      ! avoid negative and zero weights
      vi = abs(vi)
      where (vi.lt.numerical_zero) vi = numerical_zero

      deallocate(ri)
      deallocate(si)
      deallocate(kdofi)
end subroutine

!******************************************************
subroutine dd_generate_interface_unit_jump(sub, vi,lvi)
!******************************************************
! Subroutine for generating weights as a solution to local problems with unit jump.
! alpha = dT * S_i * d / (dT * (S_i + S_j) * d)
! Based on the paper
! M. Certikova, J. Sistek, P. Burda: On selection of interface weights in domain decomposition methods,
! Proceedings of Programs and Algorithms of Numerical Mathematics 16, Institute
! of Mathematics AS CR, pp. 35-44, 2013
! Constructed locally face by face.

      use module_utils
      implicit none
! Subdomain structure
      type(subdomain_type),intent(inout) :: sub
      integer,intent(in)  :: lvi
      real(kr),intent(out) :: vi(lvi)

      ! local vars
      character(*),parameter:: routine_name = 'DD_GENERATE_INTERFACE_UNIT_JUMP'
      integer :: i, nnodi, ndofi, ndofn, nglob, nnadj
      integer :: ndofnx 
      integer :: lri1, lri2 
      real(kr), allocatable :: ri(:)
      integer :: lsi1, lsi2 
      real(kr), allocatable :: si(:)
      integer :: ncol

      integer :: isubadj, jsub, ia, idofn, iglob, inadj, indn, indshni, kishnadj, nadj, ind, inodi
      integer :: range_start, range_end

      integer ::            lkdofi
      integer,allocatable :: kdofi(:)
      integer ::             lalphas
      real(kr),allocatable :: alphas(:)

      ! check if mesh is loaded
      if (.not. sub%is_mesh_loaded) then
         call error(routine_name, 'Mesh not loaded for subdomain: ',sub%isub)
      end if
      ! check if dof data are loaded
      if (.not. sub%is_interior_factorized) then
         call error(routine_name, 'Interior problem is not factorized for subdomain: ',sub%isub)
      end if
      ! check dimensions
      if (sub%ndofi.ne.lvi) then
         call error(routine_name, 'Interface dimensions mismatch for subdomain ',sub%isub)
      end if

      ! prepare array kdofi
      nnodi = sub%nnodi
      lkdofi = nnodi + 1
      allocate(kdofi(lkdofi))
      kdofi(1) = 1
      do i = 1,nnodi
         indn = sub%iin(i)
         ndofn = sub%nndf(indn)
         
         kdofi(i + 1) = kdofi(i) + ndofn
      end do

      ndofnx = minval(sub%nndf)

      ! prepare unit load vector at subdomain interface
      ndofi = sub%ndofi
      lri1 = ndofi
      lri2 = ndofnx
      allocate(ri(lri1*lri2))
      lsi1 = ndofi
      lsi2 = ndofnx
      allocate(si(lsi1*lsi2))

      lalphas = ndofnx
      allocate(alphas(lalphas))

      vi(:) = 0._kr

      if (ndofnx.gt.0) then

         ! loop over globs and identify faces
         nglob = sub%nglob
         nadj  = sub%nadj
         do iglob = 1,nglob

            ! only select faces 
            if (sub%glob_type(iglob) .eq. 1) then

               ! the second subdomain
               if (sub%nsubglobs(iglob).lt.1) then
                  call error(routine_name,'Face appears to have no neighbour for subdomain:',sub%isub)
               end if
               if (sub%lglob_subs2.lt.1) then
                  call error(routine_name,'Wrong dimension of array glob_subs for subdomain:',sub%isub)
               end if
               jsub = sub%glob_subs(iglob,1)

               ! find all nodes shared with this subdomain
               kishnadj = 0
               do ia = 1,nadj
                  nnadj = sub%nshnadj(ia)

                  ! get index of neighbour
                  isubadj = sub%iadj(ia)
                  ! perform the solves only if the subdomain matches the face neighbour
                  if (isubadj .eq. jsub) then

                     ! fill in ones in selected direction
                     ri = 0._kr
                     do inadj = 1,nnadj
                        indshni = sub%ishnadj(kishnadj + inadj)

                        do idofn = 1,ndofnx
                           ri((idofn-1)*lri1 + (kdofi(indshni)-1 + idofn)) = 1._kr
                        end do
                     end do

                     !print *, 'ri before'
                     !do i = 1,lri1
                     !   write (*,'(i2,2x,10e7.1)') i, ri(i,:)
                     !end do

                     ! multiply the vectors by Schur complement
                     ncol = ndofnx
                     call dd_multiply_by_schur(sub, ri,lri1*lri2, si,lsi1*lsi2, ncol)

                     ! complete alpha = ri^T * (S_i * ri)
                     alphas = 0._kr
                     do idofn = 1,ndofnx
                        range_start = (idofn-1)*ndofi + 1
                        range_end   = idofn*ndofi
                        alphas(idofn) = dot_product(ri(range_start:range_end),si(range_start:range_end))
                     end do

                     ! fill in ones in selected direction
                     do inadj = 1,nnadj
                        indshni = sub%ishnadj(kishnadj + inadj)

                        do idofn = 1,ndofnx
                           vi(kdofi(indshni)-1 + idofn) = vi(kdofi(indshni)-1 + idofn) &
                                                        + alphas(idofn)
                        end do
                     end do

                     !print *, 'ri after'
                     !do i = 1,lri1
                     !   write (*,'(i2,2x,10e7.1)') i, ri(i,:), vi(i)
                     !end do

                  end if
                  kishnadj = kishnadj + nnadj
               end do
            end if
         end do
      end if

      ! add ones to corners
      !ncorner = sub%ncorner  
      !do icorner = 1,ncorner
      !   indi = sub%icnsin(icorner) 

      !   do idofn = 1,ndofnx
      !      vi(kdofi(indi)-1+idofn) = 1._kr
      !   end do
      !end do

      !print *, 'vi with corners'
      !do i = 1,lvi
      !   write (*,'(i2,2x,10e7.1)') i, vi(i)
      !end do

      ! add ones to dofs above the ndofnx
      do inodi = 1,nnodi
         ind = sub%iin(inodi)
         ndofn = sub%nndf(ind)
         do idofn = ndofnx+1,ndofn
            vi(kdofi(inodi)-1+idofn) = 1._kr
         end do
      end do

      !if (any(vi.le.0._kr)) then
      !   call warning(routine_name,'zeros in weights for subdomain',sub%isub)
      !   print *, 'vi', vi
      !end if

      ! avoid negative and zero weights
      vi = abs(vi)
      where (vi.lt.numerical_zero) vi = numerical_zero

      deallocate(alphas)
      deallocate(ri)
      deallocate(si)
      deallocate(kdofi)
end subroutine

!**********************************************************
subroutine dd_get_interface_schur_diagonal(sub, rhoi,lrhoi)
!**********************************************************
! Subroutine for getting subdomain diagonal from the structure
      use module_utils
      implicit none
! Subdomain structure
      type(subdomain_type),intent(in) :: sub
      integer,intent(in)  :: lrhoi
      real(kr),intent(out) ::  rhoi(lrhoi)

      ! local vars
      integer :: ndofi
      integer :: i

      ! check dimensions
      if (sub%ndofi.ne.lrhoi) then
         write(*,*) 'DD_GET_INTERFACE_SCHUR_DIAGONAL: Interface dimensions mismatch for subdomain ',sub%isub
         call error_exit
      end if

      ! check if matrix is loaded
      if (.not. sub%is_explicit_schur_prepared) then
         write(*,*) 'DD_GET_INTERFACE_SCHUR_DIAGONAL: Explicit Schur complement is not prepared for subdomain: ',sub%isub
         call error_exit
      end if

      ! load data
      ndofi = sub%ndofi
      do i = 1,ndofi
         rhoi(i) = sub%schur(i,i)
      end do
end subroutine

!**************************************************************************
subroutine dd_dotprod_subdomain_local(sub, vec1,lvec1, vec2,lvec2, dotprod)
!**************************************************************************
! Subroutine for computing weighted dot product to be used in repeated entries with DD
! dotprod = vec1 * wi * vec2, assumes vectors of subdomain length
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
      character(*),parameter:: routine_name = 'DD_DOTPROD_SUBDOMAIN_LOCAL'
      integer :: i, indvs
      integer :: ndofi, ndofo

      ! check the prerequisities
      if (.not. (sub%is_weights_ready)) then
         write(*,*) 'DD_DOTPROD_LOCAL: Weights not ready.'
         call error_exit
      end if

      ! check dimensions
      if (lvec1 .ne. lvec2) then
         call error( routine_name, 'Dimensions mismatch.' )
      end if
      if (lvec1 .ne. sub%ndof) then
         call error( routine_name, 'Dimensions mismatch with local subdomain size.' )
      end if

      ndofi = sub%ndofi
      ndofo = sub%ndofo

      dotprod = 0._kr
      ! weight the vector at interface
      do i = 1,ndofi
         indvs = sub%iivsvn(i)
         dotprod = dotprod + vec1(indvs) * sub%wi(i) * vec2(indvs)
      end do
      ! unweighted contribution at interior
      do i = 1,ndofo
         indvs = sub%iovsvn(i)
         dotprod = dotprod + vec1(indvs) * vec2(indvs)
      end do

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
 
!*********************************
subroutine dd_print_sub(iunit,sub)
!*********************************
! Subroutine for printing the state of sub structure
      use module_sm
      implicit none
! Subdomain structure
      type(subdomain_type),intent(in) :: sub
      integer,intent(in) :: iunit


! local variables
      integer :: i, j, ia, kishnadj

      write(iunit,*) '****** start subdomain export '
      write(iunit,*) '*** HEADER INFO :             '
      write(iunit,*) '     subdomain initialized:   ', sub%is_initialized
      write(iunit,*) '     global subdomain number: ', sub%isub
      write(iunit,*) '     number of subdomains:    ', sub%nsub
      write(iunit,*) '     processor number:        ', sub%proc
      write(iunit,*) '*** MESH INFO :               '
      write(iunit,*) '     mesh loaded:             ', sub%is_mesh_loaded
      write(iunit,*) '     number of elements:      ', sub%nelem
      write(iunit,*) '     number of nodes:         ', sub%nnod
      write(iunit,*) '     number of DOF:           ', sub%ndof
      write(iunit,*) '     number of dimensions:    ', sub%ndim
      write(iunit,*) '*** BOUNDARY CONDITIONS :     '
      write(iunit,*) '     is bc present:           ', sub%is_bc_present
      write(iunit,*) '     is bc nonzero:           ', sub%is_bc_nonzero
      write(iunit,*) '     is bc loaded:            ', sub%is_bc_loaded
      write(iunit,*) '*** CORNER INFO :             '
      write(iunit,*) '     are corners loaded:      ', sub%is_corners_loaded
      write(iunit,*) '     number of corners:       ', sub%ncorner
      write(iunit,*) '*** GLOB INFO :               '
      write(iunit,*) '     are globs loaded:        ', sub%is_globs_loaded
      write(iunit,*) '     number of globs:         ', sub%nglob
      write(iunit,*) '*** NEIGHBOURING INFO :       '
      write(iunit,*) '     interface ready:         ', sub%is_interface_loaded
      write(iunit,*) '     neighbouring ready:      ', sub%is_neighbouring_ready
      write(iunit,*) '     number of neighbours:    ', sub%nadj
      write(iunit,*) '     indices of neighbours:   ', sub%iadj
      if (sub%is_neighbouring_ready) then
         kishnadj = 0
         do ia = 1,sub%nadj
            write(iunit,*) '      number of nodes shared with subdomain ',sub%iadj(ia),'is: ', sub%nshnadj(ia)
            write(iunit,*) '      indices of nodes shared: ', sub%ishnadj(kishnadj+1:kishnadj+sub%nshnadj(ia))
            kishnadj = kishnadj + sub%nshnadj(ia)
         end do
      end if
      write(iunit,*) '     indices of interface nodes:   ', sub%iin
      write(iunit,*) '*** WEIGHTS INFO :            '
      write(iunit,*) '     are weights ready?:      ', sub%is_weights_ready
      if (sub%is_weights_ready) then
         write(iunit,*) '     weights:                 ', sub%wi
      end if
      write(iunit,*) '*** COARSE NODES INFO :       '
      write(iunit,*) '     coarse nodes ready:      ', sub%is_cnodes_loaded
      write(iunit,*) '     number of coarse nodes:  ', sub%ncnodes
      if (sub%is_cnodes_loaded) then
         do i = 1,sub%ncnodes
            call dd_print_cnode(iunit,sub,i)
         end do
      end if
      write(iunit,*) '*** MATRIX INFO :             '
      write(iunit,*) '     matrix loaded:           ', sub%is_matrix_loaded
      write(iunit,*) '     matrix blocked:          ', sub%is_blocked
      write(iunit,*) '     interior block factor.:  ', sub%is_interior_factorized
      if (debug) then
         if (sub%is_matrix_loaded.and.sub%is_triplet) then
            write(iunit,*) '     matrix data:           '
            call sm_print(iunit, sub%i_a_sparse, sub%j_a_sparse, sub%a_sparse, &
                          sub%la, sub%nnza)
         end if
         if (sub%is_matrix_loaded.and.sub%is_blocked) then
            write(iunit,*) '     matrix blocks:           '
            write(iunit,*) '     A_11:                     '
            call sm_print(iunit, sub%i_a11_sparse, sub%j_a11_sparse, sub%a11_sparse, &
                          sub%la11, sub%nnza11)
            write(iunit,*) '     A_12:                     '
            call sm_print(iunit, sub%i_a12_sparse, sub%j_a12_sparse, sub%a12_sparse, &
                          sub%la12, sub%nnza12)
            write(iunit,*) '     A_21:                     '
            call sm_print(iunit, sub%i_a21_sparse, sub%j_a21_sparse, sub%a21_sparse, &
                          sub%la21, sub%nnza21)
            write(iunit,*) '     A_22:                     '
            call sm_print(iunit, sub%i_a22_sparse, sub%j_a22_sparse, sub%a22_sparse, &
                          sub%la22, sub%nnza22)
         end if
      end if
      write(iunit,*) '*** BDDC INFO:                '
      write(iunit,*) '     matrix C loaded:         ', sub%is_c_loaded
      if (debug) then
         if (sub%is_c_loaded) then
            call sm_print(iunit, sub%i_c_sparse, sub%j_c_sparse, sub%c_sparse, &
                          sub%lc, sub%nnzc)
         end if
      end if
      write(iunit,*) '     matrix Kaug factorized:  ', sub%is_aug_factorized
      if (debug) then
         if (sub%is_matrix_loaded.and.sub%is_aug_factorized) then
            if (debug) then
               call sm_print(iunit, sub%i_aaug_sparse, sub%j_aaug_sparse, sub%aaug_sparse, &
                             sub%laaug, sub%nnzaaug)
            end if
         end if
      end if
      write(iunit,*) '     matrix PHIS prepared:    ', sub%is_phisi_prepared
      if (debug) then
         if (sub%is_coarse_prepared) then
            do i = 1,sub%lphisi1
               write(iunit,'(1000f13.6)') (sub%phisi(i,j),j = 1,sub%lphisi2)
            end do
         end if
      end if
      write(iunit,*) '     matrix PHIS_DUAL prepared:    ', sub%is_phisi_dual_prepared
      if (debug) then
         if (sub%is_coarse_prepared) then
            do i = 1,sub%lphisi_dual1
               write(iunit,'(1000f13.6)') (sub%phisi_dual(i,j),j = 1,sub%lphisi_dual2)
            end do
         end if
      end if
      write(iunit,*) '     coarse matrix prepared:  ', sub%is_coarse_prepared
      if (debug) then
         if (sub%is_coarse_prepared) then
!      write(*,'(f13.6)') (sub%coarsem(j),j = 1,sub%lcoarsem)
            write(iunit,*) ' embedding of corse matrix :  '
            write(iunit,'(i8)') (sub%indrowc(j),j = 1,sub%lindrowc)
         end if
      end if
      ! PCG data
      write(iunit,*) '     reduced RHS loaded:  ', sub%is_reduced_rhs_loaded
      if (debug) then
         if (sub%is_reduced_rhs_loaded) then
!      write(*,'(f13.6)') (sub%coarsem(j),j = 1,sub%lcoarsem)
            write(iunit,*) ' reduced RHS :  '
            write(iunit,'(e15.5)') (sub%g(j),j = 1,sub%lg)
         end if
      end if
end subroutine

!******************************************
subroutine dd_print_cnode(iunit,sub,icnode)
!******************************************
! Subroutine for printing content of one coarse node
      implicit none
! Subdomain structure
      type(subdomain_type),intent(in) :: sub
      integer,intent(in) :: icnode
      integer,intent(in) :: iunit

! basic structure
      write(iunit,*) '****** start coarse node export '
      write(iunit,*) '     coarse node number:      ', icnode
      write(iunit,*) '     type of coarse node:     ', sub%cnodes(icnode)%itype
      write(iunit,*) '     used for constraints?:   ', sub%cnodes(icnode)%used
      write(iunit,*) '     coordinates:             ', sub%cnodes(icnode)%xyz
      write(iunit,*) '****** where it maps to? '
      write(iunit,*) '     global coarse node number:', sub%cnodes(icnode)%global_cnode_number
      write(iunit,*) '     number of coarse degrees of freedom:', sub%cnodes(icnode)%ncdof
!      write(iunit,*) '     indices of coarse dof:', sub%cnodes(icnode)%igcdof
!      write(iunit,*) '****** where it maps from? '
!      write(iunit,*) '     number of nodes it contains:', sub%cnodes(icnode)%nnod
!      write(iunit,*) '     indices of nodes on subdomain int:', sub%cnodes(icnode)%insin
      write(iunit,*) '     number of nodes it contains:', sub%cnodes(icnode)%nvar
      write(iunit,*) '     indices of variables on subdomain int:', sub%cnodes(icnode)%ivsivn
      write(iunit,*) '****** end coarse nodes export '
end subroutine

!************************************************
subroutine dd_plot_subdomain_data_vtu(prefix,sub)
!************************************************
! Subroutine for plotting the subdomain data in VTU format 
      use module_paraview

      implicit none
      character(*), intent(in) :: prefix           ! basename of vtu files
! Subdomain structure
      type(subdomain_type),intent(in) :: sub

      ! local vars
      integer ::  idvtu
      integer ::             lsubdomain
      integer, allocatable :: subdomain(:)

      ! write solution to a separate VTU file
      call paraview_open_subdomain_file(prefix,sub%isub,idvtu)

      ! write header of VTU file
      call paraview_write_mesh(idvtu, sub%nelem,sub%nnod, sub%ndim, sub%inet,sub%linet, sub%nnet,sub%lnnet, &
                               sub%xyz,sub%lxyz1,sub%lxyz2)

      ! write cell data
      call paraview_open_celldata(idvtu)

      lsubdomain = sub%nelem
      allocate( subdomain(lsubdomain) )
      subdomain = sub%isub
      call paraview_write_dataarray(idvtu,1,'subdomain',subdomain,lsubdomain)
      deallocate( subdomain )
      call paraview_close_celldata(idvtu)

      ! finalize the file
      call paraview_finalize_file(idvtu)
    
      ! close file
      call paraview_close_subdomain_file(idvtu)

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
      if (allocated(sub%nodal_components)) then
         deallocate(sub%nodal_components)
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
      sub%lbc = 0
      if (allocated(sub%i_a_fixed_sparse)) then
         deallocate(sub%i_a_fixed_sparse)
      end if
      if (allocated(sub%j_a_fixed_sparse)) then
         deallocate(sub%j_a_fixed_sparse)
      end if
      if (allocated(sub%a_fixed_sparse)) then
         deallocate(sub%a_fixed_sparse)
      end if
      sub%nnza_fixed = 0
      sub%la_fixed   = 0

      if (allocated(sub%user_constraints)) then
         deallocate(sub%user_constraints)
      end if
      sub%luser_constraints1 = 0
      sub%luser_constraints2 = 0

      if (allocated(sub%element_data)) then
         deallocate(sub%element_data)
      end if
      sub%lelement_data1 = 0
      sub%lelement_data2 = 0

      if (allocated(sub%dof_data)) then
         deallocate(sub%dof_data)
      end if
      sub%ldof_data = 0

      if (allocated(sub%rea)) then
         deallocate(sub%rea)
      end if
      sub%lrea = 0
      sub%is_reactions_ready = .false.
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
      if (allocated(sub%nsubglobs)) then
         deallocate(sub%nsubglobs)
      end if
      if (allocated(sub%glob_subs)) then
         deallocate(sub%glob_subs)
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
      if (allocated(sub%ishnncadj)) then
         deallocate(sub%ishnncadj)
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
         sub%is_mumps_interior_active = .false.
         sub%mumps_interior_block_factor_size = 0
      end if
      if (allocated(sub%schur)) then
         deallocate(sub%schur)
      end if
      if (allocated(sub%aaug_dense)) then
         deallocate(sub%aaug_dense)
      end if
      if (allocated(sub%aaug_ipiv)) then
         deallocate(sub%aaug_ipiv)
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
      if (allocated(sub%phis_dual)) then
         deallocate(sub%phis_dual)
      end if
      if (allocated(sub%phisi_dual)) then
         deallocate(sub%phisi_dual)
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
      if (allocated(sub%solo)) then
         deallocate(sub%solo)
      end if

      if (sub%is_mumps_aug_active) then
         if (.not. sub%is_degenerated) then
            call mumps_finalize(sub%mumps_aug)
         end if
         sub%is_mumps_aug_active = .false.
      end if

      if (allocated(sub%sol)) then
         deallocate(sub%sol)
      end if

      sub%is_mesh_loaded             = .false.
      sub%is_nodal_components_loaded = .false.
      sub%is_interface_loaded        = .false.
      sub%is_adj_loaded              = .false.
      sub%is_corners_loaded          = .false.
      sub%is_globs_loaded            = .false.
      sub%is_neighbouring_ready      = .false.
      sub%is_matrix_loaded           = .false.
      sub%is_user_constraints_loaded = .false.
      sub%is_element_data_loaded     = .false.
      sub%is_dof_data_loaded         = .false.
      sub%is_c_loaded                = .false.
      sub%is_phis_prepared           = .false.
      sub%is_phisi_prepared          = .false.
      sub%is_phis_dual_prepared      = .false.
      sub%is_phisi_dual_prepared     = .false.
      sub%is_coarse_prepared         = .false.
      sub%is_aug_factorized          = .false.
      sub%is_interior_factorized     = .false.
      sub%is_weights_ready           = .false.

end subroutine

end module module_dd

