module module_dd
!***************
! Module for handling domain decomposition structures
! Jakub Sistek, Denver, 4/2009, Praha 1/2010

!     definition of MUMPS srtucture
      use dmumps_struc_def
      use module_mumps
      implicit none

! type of real variables
      integer,parameter,private :: kr = kind(1.D0)
! numerical zero
      real(kr),parameter,private :: numerical_zero = 1.e-12_kr

! debugging 
      logical,parameter,private :: debug = .false.

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
         logical ::             is_sub_identified = .false.
         integer ::             isub    ! index of subdomain
         integer ::             nsub    ! number of subdomains

         logical ::             is_proc_assigned = .false.
         integer ::             proc    ! index of processor cantaining this subdomain

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

         integer ::             lxyz1    ! length of array of coordinates
         integer ::             lxyz2    ! number of x,y,z vectors
         real(kr),allocatable :: xyz(:,:)! array of coordinates
      
         ! description of subdomain interface
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
         integer ::             lifix         ! length of IFIX array
         integer,allocatable ::  ifix(:)      ! IFIX array - indices of fixed variables
         integer ::             lfixv         ! length of FIXV array
         real(kr),allocatable :: fixv(:)      ! FIXV array - fixed variables values
         integer ::             lbc           ! length of BC array
         real(kr),allocatable :: bc(:)        ! BC array - eliminated entries of stiffness matrix multiplied by values of fixed variables 

         ! description of corners
         integer ::             nnodc                  ! number of corners on subdomain
         integer ::             lglobal_corner_number  ! length of array GLOBAL_CORNER_NUMBER
         integer, allocatable :: global_corner_number(:) ! global numbers of these corners - length NNODC
         integer ::             licnsin                 ! length of array ICNSIN
         integer, allocatable :: icnsin(:)              ! ICNSIN array - indices of corner nodes in subdomain interface numbering
         integer ::             lncdf                   ! length of array NCDF
         integer, allocatable :: ncdf(:)                ! NCDF array - numbers of corner dof - length NNODC

         ! description of globs
         integer ::             nglob                  ! number of globs on subdomain
         integer ::             lglobal_glob_number    ! length of array GLOBAL_GLOB_NUMBER
         integer, allocatable :: global_glob_number(:) ! global numbers of these globs - lenght NGLOB
         integer ::             lnglobvar              ! length of array NGLOBVAR
         integer, allocatable :: nglobvar(:)           ! number of variables in subdomain globs - lenght NGLOB
         integer ::             ligvsivn1              ! number of rows of IGVSIVN array
         integer ::             ligvsivn2              ! number of cols of IGVSIVN array
         integer, allocatable :: igvsivn(:,:)          ! IGVSIVN array - indices of glob variables in subdomain interface numbering
         integer ::             lglob_type             ! length of array GLOB_TYPE
         integer, allocatable :: glob_type(:)          ! type of globs ( 1 - face, 2 - edge)
                                                       ! data are stored by rows
         integer ::             lngdf                  ! length of array NGDF
         integer, allocatable :: ngdf(:)               ! number of degrees of freedom associated with a glob (e.g. number of averages on glob) - lenght NGLOB

         ! common description of joint coarse nodes/dofs
         logical ::                     is_cnodes_loaded = .false. ! are coarse nodes activated?
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

         ! coarse space basis functions on interface PHISI
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

      end type subdomain_type

      integer, private ::                          lsub
      type(subdomain_type), allocatable, private :: sub(:)

contains

!***********************
subroutine dd_init(nsub)
!***********************
! Subroutine for initialization of subdomain data
      implicit none
! Global number of subdomains
      integer,intent(in) :: nsub

! local variables
      integer :: isub

! initialize basic structure
      lsub = nsub
      allocate (sub(lsub))
      do isub = 1,nsub
         sub(isub)%isub = isub
         sub(isub)%nsub = nsub
         sub(isub)%is_sub_identified = .true.
      end do

end subroutine

!**********************************************
subroutine dd_distribute_subdomains(nsub,nproc)
!**********************************************
! Subroutine for distribution of load to processors
      implicit none
! Global number of subdomains
      integer,intent(in) :: nsub
! Number of processors
      integer,intent(in) :: nproc

! local variables
      integer :: iproc, isub_loc, isub, nsub_locx

      nsub_locx = (nsub + nproc - 1)/nproc
      isub = 0
      do iproc = 0,nproc-1
         do isub_loc = 1,nsub_locx
            isub = isub + 1
            if (isub.le.nsub) then
               sub(isub)%proc = iproc
               sub(isub)%is_proc_assigned = .true.
            end if
         end do
      end do
end subroutine

!**************************************************
subroutine dd_read_mesh_from_file(myid,problemname)
!**************************************************
! Subroutine for initialization of subdomain data
      use module_utils
      implicit none

! number of processor
      integer,intent(in) :: myid
! name of problem
      character(*),intent(in) :: problemname

! local variables
      integer :: isub, idsmd, iglob, ivar
      character(100):: fname

      integer :: nelem, nnod, ndof, ndim, nnodc, nglob
      integer ::             lnndf,   lnnet,   linet 
      integer, allocatable :: nndf(:), nnet(:), inet(:)
      integer ::             lisngn
      integer, allocatable :: isngn(:)
      integer ::              lxyz1, lxyz2
      real(kr), allocatable :: xyz(:,:)
      integer :: nnodi, ndofi, ndofo
      integer ::             liin,   liivsvn,   liovsvn
      integer, allocatable :: iin(:), iivsvn(:), iovsvn(:)
      integer ::             lifix
      integer, allocatable :: ifix(:)
      integer ::              lfixv
      real(kr), allocatable :: fixv(:)
      integer ::              lglobal_corner_number
      integer, allocatable ::  global_corner_number(:)
      integer ::              licnsin 
      integer, allocatable ::  icnsin(:)
      integer ::              lglobal_glob_number
      integer, allocatable ::  global_glob_number(:)
      integer ::              lglob_type
      integer, allocatable ::  glob_type(:)
      integer ::              lnglobvar 
      integer, allocatable ::  nglobvar(:)
      integer ::              ligvsivn1, ligvsivn2
      integer, allocatable ::  igvsivn(:,:)

      do isub = 1,lsub
         if (sub(isub)%proc .eq. myid) then

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

            ! --- corners 
            read(idsmd,*) nnodc
            lglobal_corner_number = nnodc
            licnsin               = nnodc
            allocate(global_corner_number(lglobal_corner_number),icnsin(licnsin))
            read(idsmd,*) global_corner_number
            read(idsmd,*) icnsin

            ! --- globs
            read(idsmd,*) nglob
            lglobal_glob_number = nglob
            lnglobvar           = nglob
            allocate(global_glob_number(lglobal_glob_number),nglobvar(lnglobvar))
            read(idsmd,*) global_glob_number
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

            close(idsmd)
            if (debug) then
               write(*,*) 'DD_READ_MESH_FROM_FILE: Data read successfully.'
            end if

            ! load data to structure
            call dd_upload_mesh(myid,isub, nelem, nnod, ndof, ndim, nnodi, ndofi, ndofo, nnodc, nglob,&
                                nndf,lnndf, nnet,lnnet, inet,linet, isngn,lisngn,&
                                xyz,lxyz1,lxyz2, &
                                iin,liin, iivsvn,liivsvn, iovsvn,liovsvn,&
                                global_corner_number,lglobal_corner_number, icnsin,licnsin,&
                                global_glob_number,lglobal_glob_number,nglobvar,lnglobvar,&
                                igvsivn,ligvsivn1,ligvsivn2,glob_type,lglob_type)
            call dd_load_bc(myid,isub, ifix,lifix, fixv,lfixv)

            if (debug) then
               write(*,*) 'DD_READ_MESH_FROM_FILE: Data loaded successfully.'
            end if

            deallocate (nndf,nnet)
            deallocate (inet)
            deallocate (isngn)
            deallocate (xyz)
            deallocate (iin, iivsvn, iovsvn)
            deallocate (ifix,fixv)
            deallocate (glob_type)
            deallocate (global_glob_number,nglobvar)
            deallocate (igvsivn)
            deallocate (global_corner_number,icnsin)
         end if
      end do

end subroutine

!***************************************************************
subroutine dd_read_matrix_from_file(myid,problemname,matrixtype)
!***************************************************************
! Subroutine for initialization of subdomain data
      use module_utils
      use module_sm
      implicit none

! number of processor
      integer,intent(in) :: myid
! name of problem
      character(*),intent(in) :: problemname
! type of matrix
      integer,intent(in) :: matrixtype

! local variables
      character(100):: fname

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

      integer :: inod, isub

      do isub = 1,lsub
         if (sub(isub)%proc .eq. myid) then

            ! read data from file
            if (.not.sub(isub)%is_mesh_loaded) then
               write(*,*) 'DD_READ_MATRIX_FROM_FILE: Mesh is not loaded for subdomain:', isub
               call error_exit
            end if

            ! get data from subdomain
            nelem = sub(isub)%nelem
            nnod  = sub(isub)%nnod

            ! import PMD file to memory
            call sm_pmd_get_length(matrixtype,nelem,sub(isub)%inet,sub(isub)%linet,&
                                   sub(isub)%nnet,sub(isub)%lnnet,sub(isub)%nndf,sub(isub)%lnndf,&
                                   la)

            ! open subdomain ELM file with mesh description
            call getfname(problemname,isub,'ELM',fname)
            if (debug) then
               write(*,*) 'DD_READ_MATRIX_FROM_FILE: Opening file fname: ', trim(fname)
            end if
            call allocate_unit(idelm)
            open (unit=idelm,file=trim(fname),status='old',form='unformatted')

            ! Creation of field KDOF(NNOD) with addresses before first global
            ! dof of node
            lkdof = nnod
            allocate(kdof(lkdof))
            kdof(1) = 0
            do inod = 2,nnod
               kdof(inod) = kdof(inod-1) + sub(isub)%nndf(inod-1)
            end do

            ! allocate memory for triplet
            allocate(i_sparse(la), j_sparse(la), a_sparse(la))

            call sm_pmd_load(idelm,nelem,sub(isub)%inet,sub(isub)%linet,&
                             sub(isub)%nnet,sub(isub)%lnnet,sub(isub)%nndf,sub(isub)%lnndf,&
                             kdof,lkdof,&
                             i_sparse, j_sparse, a_sparse, la)
            close (idelm)
            if (debug) then
               write(*,*) 'DD_READ_MATRIX_FROM_FILE: File with element matrices read.'
            end if

            ! eliminate boundary conditions
            if (sub(isub)%is_bc_present) then
               ndof =  sub(isub)%ndof
               if (sub(isub)%is_bc_nonzero) then
                  lbc = ndof
                  allocate(bc(lbc))
               else
                  lbc = 0

               end if

               ! eliminate natural BC
               call sm_apply_bc(sub(isub)%ifix,sub(isub)%lifix,sub(isub)%fixv,sub(isub)%lfixv,&
                                i_sparse,j_sparse,a_sparse,la, bc,lbc)
               if (sub(isub)%is_bc_nonzero) then
                  call dd_load_eliminated_bc(myid, isub, bc,lbc)
               end if
            end if

            ! load matrix to our structure
            nnza = la
            is_assembled = .false.
            call dd_load_matrix_triplet(myid, isub, matrixtype, &
                                        i_sparse,j_sparse,a_sparse,la,nnza,is_assembled)


            if (debug) then
               write(*,*) 'DD_READ_MATRIX_FROM_FILE: isub =',isub,' matrix loaded.'
            end if

            deallocate(kdof)
            deallocate(i_sparse, j_sparse, a_sparse)
            if (allocated(bc)) then
               deallocate(bc)
            end if
         end if
      end do

end subroutine

!*********************************************************************************
subroutine dd_load_matrix_triplet(myid, isub, matrixtype, &
                                  i_sparse,j_sparse,a_sparse,la,nnza,is_assembled)
!*********************************************************************************
! Subroutine for loading sparse triplet to sub structure
      implicit none

      integer,intent(in) :: myid, isub
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

      ! check if I store the subdomain
      if (.not. sub(isub)%proc .eq. myid) then
         if (debug) then
            write(*,*) 'DD_LOAD_MATRIX_TRIPLET: myid =',myid,', not my subdomain: ',isub
         end if
         return
      end if

      ! load data
      sub(isub)%matrixtype = matrixtype
      sub(isub)%istorage = 2
      sub(isub)%nnza     = nnza
      sub(isub)%la       = la
      sub(isub)%is_assembled = is_assembled
      allocate(sub(isub)%i_a_sparse(la))
      do i = 1,la
         sub(isub)%i_a_sparse(i) = i_sparse(i)
      end do
      allocate(sub(isub)%j_a_sparse(la))
      do i = 1,la
         sub(isub)%j_a_sparse(i) = j_sparse(i)
      end do
      allocate(sub(isub)%a_sparse(la))
      do i = 1,la
         sub(isub)%a_sparse(i)   = a_sparse(i)
      end do

      sub(isub)%is_matrix_loaded = .true.
      sub(isub)%is_triplet       = .true.

end subroutine

!***************************************************
subroutine dd_load_eliminated_bc(myid, isub, bc,lbc)
!***************************************************
! Subroutine for loading eliminated entries of matrix multiplied by fixed
! variable values
      use module_utils
      implicit none

      integer,intent(in) :: myid, isub

      ! eliminated values of fixed variables
      integer,intent(in) :: lbc
      real(kr),intent(in)::  bc(lbc)

      ! local vars
      integer :: i

      ! check if I store the subdomain
      if (.not. sub(isub)%proc .eq. myid) then
         if (debug) then
            write(*,*) 'DD_LOAD_ELIMINATED_BC: myid =',myid,', not my subdomain: ',isub
         end if
         return
      end if

      if (.not.sub(isub)%is_bc_nonzero) then 
         write(*,*) 'DD_LOAD_ELIMINATED_BC: myid =',myid,', subdomain: ',isub,&
                    'loading BC for homogenous conditions.'
         call error_exit
      end if

      ! load eliminated boundary conditions if they are present
      if (lbc .gt. 0) then
         sub(isub)%lbc = lbc
         allocate(sub(isub)%bc(lbc))
         do i = 1,lbc
            sub(isub)%bc(i) = bc(i)
         end do
      end if

end subroutine

!*********************************************************************************
subroutine dd_upload_mesh(myid, isub, nelem, nnod, ndof, ndim, nnodi, ndofi, ndofo, nnodc, nglob,&
                          nndf,lnndf, nnet,lnnet, inet,linet, isngn,lisngn,&
                          xyz,lxyz1,lxyz2, &
                          iin,liin, iivsvn,liivsvn, iovsvn,liovsvn, &
                          global_corner_number,lglobal_corner_number, icnsin,licnsin,&
                          global_glob_number,lglobal_glob_number, nglobvar,lnglobvar,&
                          igvsivn,ligvsivn1,ligvsivn2,glob_type,lglob_type)
!*********************************************************************************
! Subroutine for loading mesh data into sub structure
      implicit none

      integer,intent(in) :: myid, isub, nelem, nnod, ndof, ndim, nnodi, ndofi, ndofo, nnodc, nglob
      integer,intent(in) :: lnndf,       lnnet,       linet
      integer,intent(in) ::  nndf(lnndf), nnet(lnnet), inet(linet)
      integer,intent(in) :: lisngn
      integer,intent(in) ::  isngn(lisngn)
      integer,intent(in) :: lxyz1, lxyz2
      real(kr),intent(in)::  xyz(lxyz1,lxyz2)
      integer,intent(in) :: liin,       liivsvn,         liovsvn
      integer,intent(in) ::  iin(liin),  iivsvn(liivsvn), iovsvn(liovsvn)
      integer,intent(in) :: lglobal_corner_number,                       licnsin
      integer,intent(in) ::  global_corner_number(lglobal_corner_number), icnsin(licnsin)
      integer,intent(in) :: lglobal_glob_number,                     lnglobvar
      integer,intent(in) ::  global_glob_number(lglobal_glob_number), nglobvar(lnglobvar)
      integer,intent(in) :: ligvsivn1, ligvsivn2
      integer,intent(in) ::  igvsivn(ligvsivn1,ligvsivn2)
      integer,intent(in) :: lglob_type
      integer,intent(in) ::  glob_type(lglob_type)

      ! local vars
      integer :: i, j

      ! check if I store the subdomain
      if (.not. sub(isub)%proc .eq. myid) then
         if (debug) then
            write(*,*) 'DD_UPLOAD_MESH: myid =',myid,', not my subdomain: ',isub
         end if
         return
      end if

      ! load data
      sub(isub)%nelem   = nelem
      sub(isub)%nnod    = nnod
      sub(isub)%ndof    = ndof
      sub(isub)%ndim    = ndim

      sub(isub)%linet   = linet
      allocate(sub(isub)%inet(linet))
      do i = 1,linet
         sub(isub)%inet(i) = inet(i)
      end do

      sub(isub)%lnnet   = lnnet
      allocate(sub(isub)%nnet(lnnet))
      do i = 1,lnnet
         sub(isub)%nnet(i) = nnet(i)
      end do

      sub(isub)%lnndf   = lnndf
      allocate(sub(isub)%nndf(lnndf))
      do i = 1,lnndf
         sub(isub)%nndf(i) = nndf(i)
      end do

      sub(isub)%lisngn   = lisngn
      allocate(sub(isub)%isngn(lisngn))
      do i = 1,lisngn
         sub(isub)%isngn(i) = isngn(i)
      end do

      sub(isub)%lxyz1   = lxyz1
      sub(isub)%lxyz2   = lxyz2
      allocate(sub(isub)%xyz(lxyz1,lxyz2))
      do j = 1,lxyz2
         do i = 1,lxyz1
            sub(isub)%xyz(i,j) = xyz(i,j)
         end do
      end do

      ! interface data
      sub(isub)%nnodi   = nnodi
      sub(isub)%ndofi   = ndofi
      sub(isub)%ndofo   = ndofo

      sub(isub)%liin = liin
      allocate(sub(isub)%iin(liin))
      do i = 1,liin
         sub(isub)%iin(i) = iin(i)
      end do

      sub(isub)%liivsvn = liivsvn
      allocate(sub(isub)%iivsvn(liivsvn))
      do i = 1,liivsvn
         sub(isub)%iivsvn(i) = iivsvn(i)
      end do

      sub(isub)%liovsvn = liovsvn
      allocate(sub(isub)%iovsvn(liovsvn))
      do i = 1,liovsvn
         sub(isub)%iovsvn(i) = iovsvn(i)
      end do

      ! corner data
      sub(isub)%nnodc = nnodc

      sub(isub)%lglobal_corner_number = lglobal_corner_number
      allocate(sub(isub)%global_corner_number(lglobal_corner_number))
      do i = 1,lglobal_corner_number
         sub(isub)%global_corner_number(i) = global_corner_number(i)
      end do

      sub(isub)%licnsin = licnsin
      allocate(sub(isub)%icnsin(licnsin))
      do i = 1,licnsin
         sub(isub)%icnsin(i) = icnsin(i)
      end do

      ! glob data
      sub(isub)%nglob = nglob

      sub(isub)%lglobal_glob_number = lglobal_glob_number
      allocate(sub(isub)%global_glob_number(lglobal_glob_number))
      do i = 1,lglobal_glob_number
         sub(isub)%global_glob_number(i) = global_glob_number(i)
      end do

      sub(isub)%lnglobvar = lnglobvar
      allocate(sub(isub)%nglobvar(lnglobvar))
      do i = 1,lnglobvar
         sub(isub)%nglobvar(i) = nglobvar(i)
      end do

      sub(isub)%ligvsivn1 = ligvsivn1
      sub(isub)%ligvsivn2 = ligvsivn2
      allocate(sub(isub)%igvsivn(ligvsivn1,ligvsivn2))
      do i = 1,ligvsivn1
         do j = 1,ligvsivn2
            sub(isub)%igvsivn(i,j) = igvsivn(i,j)
         end do
      end do

      sub(isub)%lglob_type = lglob_type
      allocate(sub(isub)%glob_type(lglob_type))
      do i = 1,lglob_type
         sub(isub)%glob_type(i) = glob_type(i)
      end do

      sub(isub)%is_mesh_loaded = .true.

end subroutine

!********************************************************
subroutine dd_load_bc(myid, isub, ifix,lifix, fixv,lfixv)
!********************************************************
! Subroutine for initialization of subdomain data
      implicit none

      integer,intent(in) :: myid, isub
      integer,intent(in) :: lifix
      integer,intent(in) ::  ifix(lifix)
      integer,intent(in) :: lfixv
      real(kr),intent(in)::  fixv(lfixv)

      ! local vars
      integer :: i
      logical :: import_bc


      ! check if I store the subdomain
      if (.not. sub(isub)%proc .eq. myid) then
         if (debug) then
            write(*,*) 'DD_LOAD_BC: myid =',myid,', not my subdomain: ',isub
         end if
         return
      end if

      if (any(ifix.ne.0)) then
         sub(isub)%is_bc_present = .true.
         import_bc = .true.
      else
         sub(isub)%is_bc_present = .false.
         import_bc = .false.
      end if
      if (any(fixv.ne.0.0_kr)) then
         sub(isub)%is_bc_nonzero = .true.
      else
         sub(isub)%is_bc_nonzero = .false.
      end if

      ! boundary conditions
      if (import_bc) then
         sub(isub)%lifix = lifix
         allocate(sub(isub)%ifix(lifix))
         do i = 1,lifix
            sub(isub)%ifix(i) = ifix(i)
         end do
   
         sub(isub)%lfixv = lfixv
         allocate(sub(isub)%fixv(lfixv))
         do i = 1,lfixv
            sub(isub)%fixv(i) = fixv(i)
         end do
      end if

end subroutine

!****************************************
subroutine dd_assembly_local_matrix(myid)
!****************************************
! Subroutine for assemblage of matrix
      use module_utils
      use module_sm
      implicit none

      integer,intent(in) :: myid

      ! local vars
      integer :: isub

      do isub = 1,lsub
         if (sub(isub)%proc .eq. myid) then

            ! check the prerequisities
            if (.not.sub(isub)%is_matrix_loaded) then
               write(*,*) 'DD_ASSEMBLY_LOCAL_MATRIX: Matrix is not loaded for subdomain:', isub
               call error_exit
            end if

            if (debug) then
               write(*,*) 'DD_ASSEMBLY_LOCAL_MATRIX: Starting assembly of sparse matrix for subdomain:',isub
            end if

            ! call matrix assemblage
            if (sub(isub)%istorage.eq.1 .or. sub(isub)%istorage.eq.2) then
               ! assembly triplet
               call sm_assembly(sub(isub)%i_a_sparse, sub(isub)%j_a_sparse, sub(isub)%a_sparse, sub(isub)%la, sub(isub)%nnza)
            else
               write(*,*) 'DD_ASSEMBLY_LOCAL_MATRIX: Matrix in unknown format type:',sub(isub)%istorage
               call error_exit
            end if

            if (debug) then
               write(*,*) 'DD_ASSEMBLY_LOCAL_MATRIX: ... done.'
            end if

         end if
      end do
end subroutine

!******************************************************
subroutine dd_matrix_tri2blocktri(myid,remove_original)
!******************************************************
! Subroutine for conversion of sparse triplet into block format A11, A12, A21, A22
      use module_utils
      implicit none

      integer,intent(in) :: myid
      ! deallocate original matrix?
      logical,intent(in) :: remove_original

      ! local vars
      integer :: isub, la, nnza, ndof, ndofi, ndofo
      logical :: is_symmetric_storage
      integer :: ia, inddofi, inddofo, idofi, idofo, irow, jcol
      integer :: la11, la12, la21, la22 
      integer :: ia11, ia12, ia21, ia22 
      real(kr) :: aval

      ! mask for interface entries
      integer ::            lmaski,   lmasko
      integer,allocatable :: maski(:), masko(:)

      do isub = 1,lsub
         if (sub(isub)%proc .eq. myid) then

            ! check the prerequisities
            if (.not.sub(isub)%is_mesh_loaded) then
               write(*,*) 'DD_MATRIX_TRI2BLOCKTRI: Mesh is not loaded for subdomain:', isub
               call error_exit
            end if
            if (.not.sub(isub)%is_matrix_loaded) then
               write(*,*) 'DD_MATRIX_TRI2BLOCKTRI: Matrix is not loaded for subdomain:', isub
               call error_exit
            end if
            if (sub(isub)%istorage.eq.3 .or. sub(isub)%istorage.eq.4) then
               write(*,*) 'DD_MATRIX_TRI2BLOCKTRI: Matrix already in block triple format, subdomain.'
               cycle
            end if
            if (.not.(sub(isub)%istorage.eq.1 .or. sub(isub)%istorage.eq.2)) then
               write(*,*) 'DD_MATRIX_TRI2BLOCKTRI: Matrix not in single block triple format, subdomain:', isub,&
                          'type:',sub(isub)%istorage
               call error_exit
            end if

            if (sub(isub)%istorage .eq. 1) then
               is_symmetric_storage = .false.
            else
               is_symmetric_storage = .true.
            end if

            ! prepare mask of interface and interior nodes
            ndof   = sub(isub)%ndof
            lmaski = ndof
            lmasko = ndof
            allocate(maski(lmaski),masko(lmasko))
            call zero(maski,lmaski)
            call zero(masko,lmasko)

            ndofi  = sub(isub)%ndofi
            do idofi = 1,ndofi
               inddofi = sub(isub)%iivsvn(idofi)
               maski(inddofi) = idofi
            end do
            ndofo  = sub(isub)%ndofo
            do idofo = 1,ndofo
               inddofo = sub(isub)%iovsvn(idofo)
               masko(inddofo) = idofo 
            end do

            if (debug) then
               ! check consistency of masks
               if (any(maski.ne.0 .and. masko.ne.0)) then
                  write(*,*) 'DD_MATRIX_TRI2BLOCKTRI: Check of mask consistency failed for subdomain:',isub
                  call error_exit
               end if
               if (any(maski.eq.0 .and. masko.eq.0)) then
                  write(*,*) 'DD_MATRIX_TRI2BLOCKTRI: Check of mask coverage failed for subdomain:',isub
                  call error_exit
               end if
            end if

            ! count sizes of matrix blocks in question
            la   = sub(isub)%la
            nnza = sub(isub)%nnza
            la11 = 0
            la12 = 0
            la21 = 0
            la22 = 0
            do ia = 1,nnza
               irow = sub(isub)%i_a_sparse(ia)
               jcol = sub(isub)%j_a_sparse(ia)

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
            sub(isub)%la11   = la11
            sub(isub)%nnza11 = la11
            sub(isub)%la22   = la22
            sub(isub)%nnza22 = la22
            if (is_symmetric_storage) then
               la12 = la12 + la21
               sub(isub)%la12   = la12
               sub(isub)%nnza12 = la12
               sub(isub)%la21   = 0
               sub(isub)%nnza21 = 0
            else
               sub(isub)%la12   = la12
               sub(isub)%nnza12 = la12
               sub(isub)%la21   = la21
               sub(isub)%nnza21 = la21
            end if
            allocate(sub(isub)%i_a11_sparse(la11),sub(isub)%j_a11_sparse(la11),sub(isub)%a11_sparse(la11))
            allocate(sub(isub)%i_a22_sparse(la22),sub(isub)%j_a22_sparse(la22),sub(isub)%a22_sparse(la22))

            if (is_symmetric_storage) then
               allocate(sub(isub)%i_a12_sparse(la12),sub(isub)%j_a12_sparse(la12),sub(isub)%a12_sparse(la12))
            else
               allocate(sub(isub)%i_a12_sparse(la12),sub(isub)%j_a12_sparse(la12),sub(isub)%a12_sparse(la12))
               allocate(sub(isub)%i_a21_sparse(la21),sub(isub)%j_a21_sparse(la21),sub(isub)%a21_sparse(la21))
            end if

            ! convert matrix to blocks according to interface - denoted by 1 in MASKI array
            ia11 = 0
            ia12 = 0
            ia21 = 0
            ia22 = 0
            do ia = 1,nnza
               irow = sub(isub)%i_a_sparse(ia)
               jcol = sub(isub)%j_a_sparse(ia)
               aval = sub(isub)%a_sparse(ia)

               ! diagonal blocks
               if      (maski(irow).eq.0 .and. maski(jcol).eq.0) then
                  ia11 = ia11 + 1
                  sub(isub)%i_a11_sparse(ia11) = masko(irow)
                  sub(isub)%j_a11_sparse(ia11) = masko(jcol)
                  sub(isub)%a11_sparse(ia11)   = aval
               else if (maski(irow).ne.0 .and. maski(jcol).ne.0) then
                  ia22 = ia22 + 1
                  sub(isub)%i_a22_sparse(ia22) = maski(irow)
                  sub(isub)%j_a22_sparse(ia22) = maski(jcol)
                  sub(isub)%a22_sparse(ia22)   = aval
               ! offdiagonal blocks
               else if (maski(irow).eq.0 .and. maski(jcol).ne.0) then
                  ia12 = ia12 + 1
                  sub(isub)%i_a12_sparse(ia12) = masko(irow)
                  sub(isub)%j_a12_sparse(ia12) = maski(jcol)
                  sub(isub)%a12_sparse(ia12)   = aval
               else if (maski(irow).ne.0 .and. maski(jcol).eq.0) then
                  if (is_symmetric_storage) then
                     ia12 = ia12 + 1
                     sub(isub)%i_a12_sparse(ia12) = masko(jcol)
                     sub(isub)%j_a12_sparse(ia12) = maski(irow)
                     sub(isub)%a12_sparse(ia12)   = aval
                  else
                     ia21 = ia21 + 1
                     sub(isub)%i_a21_sparse(ia21) = maski(irow)
                     sub(isub)%j_a21_sparse(ia21) = masko(jcol)
                     sub(isub)%a21_sparse(ia21)   = aval
                  end if
               end if
            end do

            if (is_symmetric_storage) then
               sub(isub)%istorage = 4
            else
               sub(isub)%istorage = 3
            end if

            sub(isub)%is_blocked = .true.

            deallocate(maski,masko)

            if (remove_original) then
               deallocate(sub(isub)%i_a_sparse, sub(isub)%j_a_sparse, sub(isub)%a_sparse)
               sub(isub)%la   = 0
               sub(isub)%nnza = 0
               sub(isub)%is_triplet = .false.
            end if

         end if
      end do

end subroutine

!******************************************
subroutine dd_prepare_schur(myid,comm,isub)
!******************************************
! Subroutine for preparing data for computing with reduced problem
      use module_utils
      use module_mumps
      implicit none

      ! processor ID
      integer,intent(in) :: myid
      ! communicator
      integer,intent(in) :: comm
      ! subdomain
      integer,intent(in) :: isub

      ! local vars
      integer :: ndofo, la11, nnza11
      integer :: mumpsinfo

      ! check the prerequisities
      if (.not.sub(isub)%is_proc_assigned) then
         write(*,*) 'DD_PREPARE_SCHUR: Processor not assigned for subdomain:', isub
         call error_exit
      end if
      if (sub(isub)%proc .ne. myid) then
         if (debug) then
            write(*,*) 'DD_PREPARE_SCHUR: myid =',myid,'Subdomain', isub,' is not mine.'
         end if
         return
      end if
      if (.not.sub(isub)%is_matrix_loaded) then
         write(*,*) 'DD_PREPARE_SCHUR: Matrix is not loaded for subdomain:', isub
         call error_exit
      end if
      if (.not. (sub(isub)%is_blocked)) then
         write(*,*) 'DD_PREPARE_SCHUR: Matrix is not in blocked format. Call routine to do this.'
         call error_exit
      end if

      sub(isub)%is_mumps_interior_active = .true.

      ! Initialize MUMPS
      call mumps_init(sub(isub)%mumps_interior_block,comm,sub(isub)%matrixtype)
      ! Level of information from MUMPS
      if (debug) then
         mumpsinfo = 2
      else
         mumpsinfo = 0
      end if
      call mumps_set_info(sub(isub)%mumps_interior_block,mumpsinfo)
      ! Load matrix to MUMPS
      ndofo  = sub(isub)%ndofo
      nnza11 = sub(isub)%nnza11
      la11   = sub(isub)%la11
      call mumps_load_triplet(sub(isub)%mumps_interior_block,ndofo,nnza11,&
                              sub(isub)%i_a11_sparse,sub(isub)%j_a11_sparse,sub(isub)%a11_sparse,la11)
      ! Analyze matrix
      call mumps_analyze(sub(isub)%mumps_interior_block) 
      ! Factorize matrix 
      call mumps_factorize(sub(isub)%mumps_interior_block) 

      sub(isub)%is_interior_factorized = .true.

end subroutine

!*********************************************************
subroutine dd_load_arithmetic_constraints(myid,isub,itype)
!*********************************************************
! Subroutine for assemblage of matrix
      use module_utils
      implicit none
      integer,intent(in) :: myid
      integer,intent(in) :: isub
      integer,intent(in) :: itype ! type of globs (2 - edges, 1 - faces)

      ! local vars
      integer :: icnode, lmatrix1, lmatrix2, ncnodes, ncdof, nvar
      integer :: i, j

      ! check the prerequisities
      if (sub(isub)%proc .ne. myid) then
         if (debug) then
            write(*,*) 'DD_LOAD_ARITHMETIC_CONSTRAINTS: myid =',myid,'Subdomain', isub,' is not mine.'
         end if
         return
      end if
      if (.not.allocated(sub(isub)%cnodes)) then
         write(*,*) 'DD_LOAD_ARITHMETIC_CONSTRAINTS: Array for cnodes not ready.'
         call error_exit
      end if

      ! get number of coarse nodes
      ncnodes = sub(isub)%ncnodes

      ! get number of constraints on an arithmetic constraint
      ! ndim for arithmetic averages
      ncdof = sub(isub)%ndim

      ! generate arithmetic averages on coarse nodes of prescribed type (e.g. edges)
      do icnode = 1,ncnodes
         if (sub(isub)%cnodes(icnode)%itype .eq. itype) then
         
            nvar = sub(isub)%cnodes(icnode)%nvar

            lmatrix1 = ncdof
            lmatrix2 = nvar

            sub(isub)%cnodes(icnode)%ncdof    = ncdof
            sub(isub)%cnodes(icnode)%nnz      = nvar
            sub(isub)%cnodes(icnode)%lmatrix1 = lmatrix1
            sub(isub)%cnodes(icnode)%lmatrix2 = lmatrix2
            allocate(sub(isub)%cnodes(icnode)%matrix(lmatrix1,lmatrix2))

            call zero(sub(isub)%cnodes(icnode)%matrix,lmatrix1,lmatrix2)

            do i = 1,lmatrix1
               do j = i,lmatrix2,lmatrix1
                  sub(isub)%cnodes(icnode)%matrix(i,j) = 1._kr
               end do
            end do
            !print *, 'Matrix C for an edge'
            !do i = 1,lmatrix1
            !   print '(100f5.1)', (sub(isub)%cnodes(icnode)%matrix(i,j),j=1,lmatrix2)
            !end do
            sub(isub)%cnodes(icnode)%arithmetic = .true.

            ! mark the coarse node as used
            sub(isub)%cnodes(icnode)%used = .true.
         end if
      end do

end subroutine

!***************************************************************************
subroutine dd_load_adaptive_constraints(isub,gglob,cadapt,lcadapt1,lcadapt2)
!***************************************************************************
! Subroutine for assemblage of matrix
      use module_utils
      implicit none
      integer,intent(in) :: isub
      integer,intent(in) :: gglob

      integer,intent(in) :: lcadapt1, lcadapt2
      real(kr),intent(in) :: cadapt(lcadapt1,lcadapt2)

      ! local vars
      integer :: ind_loc, lmatrix1, lmatrix2, nvarglb, ncdof
      integer :: i, j, indiv

      ! check the prerequisities
      if (.not.allocated(sub(isub)%cnodes)) then
         write(*,*) 'DD_LOAD_ADAPTIVE_CONSTRAINTS: Array for cnodes not ready.'
         call error_exit
      end if

      ! find local (subdomain) index of the glob from its global number
      call get_index(gglob,sub(isub)%cnodes%global_cnode_number,sub(isub)%ncnodes,ind_loc)
      if (ind_loc.le.0) then
         write(*,*) 'DD_LOAD_ADAPTIVE_CONSTRAINTS: Index of glob not found!'
         call error_exit
      end if

      ! prepare space for these constraints in the structure
      nvarglb = sub(isub)%cnodes(ind_loc)%nvar 

      if (nvarglb.gt.lcadapt1) then
         write(*,*) 'DD_LOAD_ADAPTIVE_CONSTRAINTS: Dimension of matrix of averages is smaller than number of variables on glob.'
         write(*,*) 'DD_LOAD_ADAPTIVE_CONSTRAINTS: nvarglb =',nvarglb, 'lcadapt1 =',lcadapt1
         call error_exit
      end if

      ! regularize the matrix of constraints by QR decomposition
      ! TODO

      ! copy transposed selected variables to the global structure
      ncdof = lcadapt2
      sub(isub)%cnodes(ind_loc)%ncdof = ncdof

      lmatrix1 = lcadapt2
      lmatrix2 = nvarglb
      sub(isub)%cnodes(ind_loc)%lmatrix1 = lmatrix1
      sub(isub)%cnodes(ind_loc)%lmatrix2 = lmatrix2
      allocate(sub(isub)%cnodes(ind_loc)%matrix(lmatrix1,lmatrix2))
      
      do i = 1,lmatrix1
         do j = 1,nvarglb
            indiv = sub(isub)%cnodes(ind_loc)%ivsivn(j)

            sub(isub)%cnodes(ind_loc)%matrix(i,j) = cadapt(indiv,i)
         end do
      end do
      sub(isub)%cnodes(ind_loc)%nnz = lmatrix1*lmatrix2
      sub(isub)%cnodes(ind_loc)%adaptive = .true.

      ! mark the coarse node as used
      sub(isub)%cnodes(ind_loc)%used = .true.

      if (debug) then
         write(*,*) 'DD_LOAD_ADAPTIVE_CONSTRAINTS: Loading matrix of subdomain ',isub,' local glob #',ind_loc
         do i = 1,lmatrix1
            print '(100f15.6)',(sub(isub)%cnodes(ind_loc)%matrix(i,j),j = 1,lmatrix2)
         end do
      end if
end subroutine

!**********************************************************************
subroutine dd_get_adaptive_constraints_size(myid,isub,iglb,lavg1,lavg2)
!**********************************************************************
! Subroutine for inquiring sizes to allocate for number of adaptive averages
      use module_utils
      implicit none
      ! processor ID
      integer,intent(in) :: myid
      ! subdomain number
      integer,intent(in) :: isub
      ! local (subdomain) glob number
      integer,intent(in) :: iglb

      ! sizes of matrix of averages
      integer,intent(out) :: lavg1, lavg2

      ! check the prerequisities
      if (.not.allocated(sub)) then
         write(*,*) 'DD_GET_ADAPTIVE_CONSTRAINTS_SIZE: Main DD structure is not ready.'
         call error_exit
      end if
      if (sub(isub)%proc .ne. myid) then
         write(*,*) 'DD_GET_ADAPTIVE_CONSTRAINTS_SIZE: Not my subdomain.'
         call error_exit
      end if
      if (.not.allocated(sub(isub)%cnodes)) then
         write(*,*) 'DD_GET_ADAPTIVE_CONSTRAINTS_SIZE: Array for coarse nodes constraints not ready.'
         call error_exit
      end if
      if (.not.allocated(sub(isub)%cnodes(iglb)%matrix)) then
         write(*,*) 'isub',isub,'iglb',iglb,'type',sub(isub)%cnodes(iglb)%itype
         write(*,*) 'DD_GET_ADAPTIVE_CONSTRAINTS_SIZE: Matrix constraints are not allocated.'
         call error_exit
      end if

      lavg1 = sub(isub)%cnodes(iglb)%lmatrix1
      lavg2 = sub(isub)%cnodes(iglb)%lmatrix2
end subroutine

!*********************************************************************
subroutine dd_get_adaptive_constraints(myid,isub,iglb,avg,lavg1,lavg2)
!*********************************************************************
! Subroutine for getting adaptive averages from the structure
      use module_utils
      implicit none
      ! processor ID
      integer,intent(in) :: myid
      ! subdomain number
      integer,intent(in) :: isub
      ! local (subdomain) glob number
      integer,intent(in) :: iglb

      ! matrix of averages
      integer,intent(in) ::  lavg1, lavg2
      real(kr),intent(out) :: avg(lavg1,lavg2)

      ! local variables
      integer :: i,j

      ! check the prerequisities
      if (.not.allocated(sub)) then
         write(*,*) 'DD_GET_ADAPTIVE_CONSTRAINTS_SIZE: Main DD structure is not ready.'
         call error_exit
      end if
      if (sub(isub)%proc .ne. myid) then
         write(*,*) 'DD_GET_ADAPTIVE_CONSTRAINTS_SIZE: Not my subdomain ',isub
         return
      end if
      if (sub(isub)%proc .ne. myid) then
         write(*,*) 'DD_GET_ADAPTIVE_CONSTRAINTS_SIZE: Not my subdomain.'
         call error_exit
      end if
      if (.not.allocated(sub(isub)%cnodes)) then
         write(*,*) 'DD_GET_ADAPTIVE_CONSTRAINTS_SIZE: Array for coarse nodes constraints not ready.'
         call error_exit
      end if
      if (sub(isub)%cnodes(iglb)%lmatrix1.ne.lavg1 .or. &
          sub(isub)%cnodes(iglb)%lmatrix2.ne.lavg2) then
         write(*,*) 'DD_GET_ADAPTIVE_CONSTRAINTS_SIZE: Matrix dimensions for averages do not match.'
         call error_exit
      end if

      do j = 1,lavg2
         do i = 1,lavg1
            avg(i,j) = sub(isub)%cnodes(iglb)%matrix(i,j)
         end do
      end do

end subroutine

!**********************************
subroutine dd_get_cnodes(myid,isub)
!**********************************
! Merging corners with globs - order first corners, then globs
      use module_sm
      use module_utils
      implicit none

      ! processor ID
      integer,intent(in) :: myid
      ! subdomain
      integer,intent(in) :: isub

      ! local vars
      integer ::             nnodc
      integer ::             nglob
      integer ::             ncnodes, lcnodes
      integer ::             icnode
      integer ::             inodc, indnode, indinode, i, kcdof, ncdof, nvar, iglob, nnodgl
      integer ::             lxyz


      ! check the prerequisities
      if (sub(isub)%proc .ne. myid) then
         if (debug) then
            write(*,*) 'DD_GET_CNODES: myid =',myid,'Subdomain', isub,' is not mine.'
         end if
         return
      end if
      if (.not.sub(isub)%is_mesh_loaded) then
         write(*,*) 'DD_GET_CNODES: Mesh is not loaded for subdomain:', isub
         call error_exit
      end if

      ! determine number of coarse nodes
      nnodc = sub(isub)%nnodc
      nglob = sub(isub)%nglob
      ncnodes = nnodc + nglob
      sub(isub)%ncnodes = ncnodes

      lcnodes = ncnodes
      allocate(sub(isub)%cnodes(lcnodes))

      ! set counter
      icnode = 0

      ! copy corners
      do inodc = 1,nnodc
         icnode = icnode + 1

         ! type of coarse node - corner
         sub(isub)%cnodes(icnode)%itype = 3
         ! is coarse node used?
         sub(isub)%cnodes(icnode)%used  = .true.
         ! global number
         sub(isub)%cnodes(icnode)%global_cnode_number  = sub(isub)%global_corner_number(inodc)
         ! number of coarse dof it contains
         ncdof = sub(isub)%ndim
         sub(isub)%cnodes(icnode)%ncdof = ncdof
         ! number of nodes where it maps from
         sub(isub)%cnodes(icnode)%nnod = 1
         ! number of variables it maps from 
         nvar = sub(isub)%ndim
         sub(isub)%cnodes(icnode)%nvar = nvar
         ! number of nonzeros it creates in matrix C
         sub(isub)%cnodes(icnode)%nnz = nvar

         ! fill coordinates
         lxyz = sub(isub)%ndim
         sub(isub)%cnodes(icnode)%lxyz = lxyz
         allocate(sub(isub)%cnodes(icnode)%xyz(lxyz))
         indinode = sub(isub)%icnsin(inodc)
         indnode  = sub(isub)%iin(indinode)
         sub(isub)%cnodes(icnode)%xyz = sub(isub)%xyz(indnode,:)

         ! fill coarse node dof
         allocate(sub(isub)%cnodes(icnode)%igcdof(ncdof))
         kcdof = (sub(isub)%global_corner_number(inodc)-1)*nvar 
         do i = 1,ncdof
            sub(isub)%cnodes(icnode)%igcdof(i) = kcdof + i
         end do

         ! fill coarse node nodes
         allocate(sub(isub)%cnodes(icnode)%insin(1))
         sub(isub)%cnodes(icnode)%insin(1) = indinode

         ! fill coarse nodes variables
         allocate(sub(isub)%cnodes(icnode)%ivsivn(nvar))
         do i = 1,nvar
            sub(isub)%cnodes(icnode)%ivsivn(i) = (indinode-1)*nvar + i
         end do

      end do

      ! copy globs behind corners - DO NOT nullify the icnode counter to place them behind corners
      do iglob = 1,nglob
         icnode = icnode + 1

         ! type of coarse node - edge or face
         sub(isub)%cnodes(icnode)%itype = sub(isub)%glob_type(iglob)
         ! is coarse node used?
         ! this is set by routines dd_load_arithmetic_constraints and
         !                         dd_load_adaptive_constraints
         ! default behaviour for globs 
         !if      (sub(isub)%glob_type(iglob) .eq. 2 .and. use_edges) then
         !   sub(isub)%cnodes(icnode)%used  = .true.
         !else if (sub(isub)%glob_type(iglob) .eq. 1 .and. use_faces_arithmetic) then
         !   sub(isub)%cnodes(icnode)%used  = .true.
         !   sub(isub)%cnodes(icnode)%arithmetic = .true.
         !else if (sub(isub)%glob_type(iglob) .eq. 1 .and. use_faces_adaptive) then
         !   sub(isub)%cnodes(icnode)%used  = .true.
         !   sub(isub)%cnodes(icnode)%adaptive = .true.
         !end if

         ! global number
         sub(isub)%cnodes(icnode)%global_cnode_number = sub(isub)%global_glob_number(iglob)
         ! number of coarse dof it contains
         ! ndim for arithmetic averages
         !ncdof = sub(isub)%ndim
         !sub(isub)%cnodes(icnode)%ncdof = ncdof
         ! number of nodes where it maps from
         nnodgl = sub(isub)%nglobvar(iglob)/sub(isub)%ndim
         sub(isub)%cnodes(icnode)%nnod = nnodgl
         ! number of variables it maps from 
         nvar = sub(isub)%nglobvar(iglob)
         sub(isub)%cnodes(icnode)%nvar = nvar

         ! fill coordinates
         lxyz = sub(isub)%ndim
         sub(isub)%cnodes(icnode)%lxyz = lxyz
         allocate(sub(isub)%cnodes(icnode)%xyz(lxyz))
         !TODO: create averaged coordinates
         sub(isub)%cnodes(icnode)%xyz = 0.

         ! fill coarse node dof
         allocate(sub(isub)%cnodes(icnode)%igcdof(ncdof))
         !CONTINU HERE! - does not work for adaptive constraints, only for arithmetic averages
         kcdof = (sub(isub)%global_glob_number(iglob)-1)*nvar
         do i = 1,ncdof
            sub(isub)%cnodes(icnode)%igcdof(i) = kcdof + i
         end do

         ! fill coarse node nodes
         ! TODO: not for globs

         ! fill coarse nodes variables
         allocate(sub(isub)%cnodes(icnode)%ivsivn(nvar))
         do i = 1,nvar
            sub(isub)%cnodes(icnode)%ivsivn(i) = sub(isub)%igvsivn(iglob,i)
         end do

         sub(isub)%is_cnodes_loaded = .true.

      end do

      return
end subroutine

!*********************************
subroutine dd_prepare_c(myid,isub)
!*********************************
! Subroutine for preparing matrix of constraints C in sparse triplet format
      use module_sm
      use module_utils
      implicit none

      ! processor ID
      integer,intent(in) :: myid
      ! subdomain
      integer,intent(in) :: isub

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
      if (sub(isub)%proc .ne. myid) then
         if (debug) then
            write(*,*) 'DD_PREPARE_C: myid =',myid,'Subdomain', isub,' is not mine.'
         end if
         return
      end if
      if (.not.sub(isub)%is_mesh_loaded) then
         write(*,*) 'DD_PREPARE_C: Mesh is not loaded for subdomain:', isub
         call error_exit
      end if
      if (.not.sub(isub)%is_matrix_loaded) then
         write(*,*) 'DD_PREPARE_C: Matrix is not loaded for subdomain:', isub
         call error_exit
      end if

      ! find number of rows in C (constraints)
      nrowc = 0
      ! find number of nonzeros in C
      nnzc  = 0

      ncnodes = sub(isub)%ncnodes
      do icn = 1,ncnodes
         if (sub(isub)%cnodes(icn)%used) then
            nrowc = nrowc + sub(isub)%cnodes(icn)%ncdof
            nnzc  = nnzc + sub(isub)%cnodes(icn)%nnz
         end if
      end do
      nconstr = nrowc


      ! Creation of field KDOF(NNOD) with addresses before first global
      ! dof of node
      nnod  = sub(isub)%nnod
      lkdof = nnod
      allocate(kdof(lkdof))
      kdof(1) = 0
      do inod = 2,nnod
         kdof(inod) = kdof(inod-1) + sub(isub)%nndf(inod-1)
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
         if (sub(isub)%cnodes(icn)%used) then

            if (sub(isub)%cnodes(icn)%itype .eq. 3) then
               ! for corners, do not read matrices but construct the sparse matrix directly
               nvar = sub(isub)%cnodes(icn)%nvar
               do ivar = 1,nvar
                  irowc = irowc + 1
                  inzc  = inzc  + 1

                  i_c_sparse(inzc) = irowc
                  j_c_sparse(inzc) = sub(isub)%cnodes(icn)%ivsivn(ivar)
                  c_sparse(inzc)   = 1._kr

                  indrowc(irowc) = sub(isub)%cnodes(icn)%igcdof(ivar)
               end do
            else if (sub(isub)%cnodes(icn)%itype .eq. 2 .or. &
                     (sub(isub)%cnodes(icn)%itype .eq. 1 .and.&
                      sub(isub)%cnodes(icn)%arithmetic .eqv. .true.)) then
               ! use arithmetic averages for edges and faces if desired
               nvar  = sub(isub)%cnodes(icn)%nvar
               ncdof = sub(isub)%cnodes(icn)%ncdof
               do icdof = 1,ncdof
                  irowc = irowc + 1

                  do ivar = icdof,nvar,ncdof
                     inzc  = inzc  + 1

                     i_c_sparse(inzc) = irowc
                     j_c_sparse(inzc) = sub(isub)%cnodes(icn)%ivsivn(ivar)
                     c_sparse(inzc)   = 1._kr
                  end do

                  indrowc(irowc) = sub(isub)%cnodes(icn)%igcdof(icdof)
               end do
            else if ((sub(isub)%cnodes(icn)%itype .eq. 1 .and.&
                      sub(isub)%cnodes(icn)%adaptive .eqv. .true.)) then
               ! use adaptive constraint on face

               ! copy the matrix of constraints on glob into sparse triplet of subdomain matrix C
               ! row by row
               nvar      = sub(isub)%cnodes(icn)%nvar
               lmatrix2  = sub(isub)%cnodes(icn)%lmatrix2
               if (nvar.ne.lmatrix2) then
                  write(*,*) 'DD_PREPARE_C: Matrix dimension does not match for subdomain', isub
                  call error_exit
               end if

               ncdof = sub(isub)%cnodes(icn)%ncdof
               do icdof = 1,ncdof
                  irowc = irowc + 1

                  do ivar = 1,nvar
                     inzc  = inzc  + 1

                     i_c_sparse(inzc) = irowc
                     j_c_sparse(inzc) = sub(isub)%cnodes(icn)%ivsivn(ivar)
                     c_sparse(inzc)   = sub(isub)%cnodes(icn)%matrix(icdof,ivar)
                  end do

                  indrowc(irowc) = sub(isub)%cnodes(icn)%igcdof(icdof)
               end do
            else
               continue
            end if

            ! check matrix bounds
            if (inzc .gt. lc) then
               write(*,*) 'DD_PREPARE_C: Too many entries in matrix C for subdomain', isub
               call error_exit
            end if
         end if
      end do

      ! check matrix bounds
      if (inzc .ne. lc) then
         write(*,*) 'DD_PREPARE_C: Dimension of matrix Cmismatch for subdomain', isub
         call error_exit
      end if
      
      ! load the new matrix directly to the structure
      call dd_load_c(myid,isub,nconstr,i_c_sparse, j_c_sparse, c_sparse, lc, nnzc, indrowc,lindrowc)

      deallocate (indrowc)
      deallocate (i_c_sparse,j_c_sparse,c_sparse)
      deallocate (kdof)

end subroutine

!***************************************************************************************************
subroutine dd_load_c(myid,isub,nconstr,i_c_sparse, j_c_sparse, c_sparse, lc, nnzc, indrowc,lindrowc)
!***************************************************************************************************
! Subroutine for loading sparse matrix C into the SUB structure
      use module_utils
      implicit none

      ! processor ID
      integer,intent(in) :: myid
      ! subdomain
      integer,intent(in) :: isub
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
      if (sub(isub)%proc .ne. myid) then
         if (debug) then
            write(*,*) 'DD_LOAD_C: myid =',myid,'Subdomain', isub,' is not mine.'
         end if
         return
      end if
      if (sub(isub)%is_c_loaded) then
         ! replace the loaded C
         if (debug) then
            write(*,*) 'DD_LOAD_C: myid =',myid,'Subdomain', isub,': Matrix C already allocated. Rewriting.'
         end if
         deallocate(sub(isub)%i_c_sparse,sub(isub)%j_c_sparse,sub(isub)%c_sparse)
         deallocate(sub(isub)%indrowc)
         sub(isub)%is_c_loaded = .false.
      end if

      ! load C now
      sub(isub)%nconstr = nconstr
      sub(isub)%nnzc = nnzc
      sub(isub)%lc   = lc
      allocate(sub(isub)%i_c_sparse(lc))
      do i = 1,lc
         sub(isub)%i_c_sparse(i) = i_c_sparse(i)
      end do
      allocate(sub(isub)%j_c_sparse(lc))
      do i = 1,lc
         sub(isub)%j_c_sparse(i) = j_c_sparse(i)
      end do
      allocate(sub(isub)%c_sparse(lc))
      do i = 1,lc
         sub(isub)%c_sparse(i)   = c_sparse(i)
      end do

      ! load indrowc
      sub(isub)%lindrowc = lindrowc
      allocate(sub(isub)%indrowc(lindrowc))
      do i = 1,lindrowc
         sub(isub)%indrowc(i)   = indrowc(i)
      end do
         
      sub(isub)%is_c_loaded = .true.

end subroutine

!****************************************
subroutine dd_prepare_aug(myid,comm,isub)
!****************************************
! Subroutine for preparing augmented matrix for computing with BDDC preconditioned
! |A C^T|
! |C  0 |
! and its factorization
      use module_sm
      use module_mumps
      use module_utils
      implicit none

      ! processor ID
      integer,intent(in) :: myid
      ! communicator
      integer,intent(in) :: comm
      ! subdomain
      integer,intent(in) :: isub

      ! local vars
      integer ::  nnzc, nnza, ndof, nconstr
      integer ::  nnzaaug, laaug, ndofaaug
      integer ::  i, iaaug
      integer ::  mumpsinfo, aaugmatrixtype
      integer ::  icol, icoli

      ! check the prerequisities
      if (sub(isub)%proc .ne. myid) then
         if (debug) then
            write(*,*) 'DD_PREPARE_AUG: myid =',myid,'Subdomain', isub,' is not mine.'
         end if
         return
      end if
      if (.not.sub(isub)%is_matrix_loaded) then
         write(*,*) 'DD_PREPARE_AUG: Matrix is not loaded for subdomain:', isub
         call error_exit
      end if
      if (.not.sub(isub)%is_c_loaded) then
         write(*,*) 'DD_PREPARE_AUG: Matrix of constraints C not loaded for subdomain:', isub
         call error_exit
      end if

      ! join the new matrix directly to the structure as
      ! in the unsymmetric case:
      ! A C^T
      ! C  0   
      ! in the symmetric case :
      ! \A C^T
      !     0   

      ndof     = sub(isub)%ndof
      nconstr  = sub(isub)%nconstr
      ndofaaug = ndof + nconstr
      nnza     = sub(isub)%nnza
      nnzc     = sub(isub)%nnzc

      if      (sub(isub)%matrixtype .eq. 0) then
         ! unsymmetric case:
         nnzaaug = nnza + 2*nnzc
      else if (sub(isub)%matrixtype .eq. 1 .or. sub(isub)%matrixtype .eq. 2) then
         ! symmetric case:
         nnzaaug = nnza + nnzc
      end if

      laaug   = nnzaaug
      allocate(sub(isub)%i_aaug_sparse(laaug),sub(isub)%j_aaug_sparse(laaug),sub(isub)%aaug_sparse(laaug))

      ! copy entries of A
      iaaug = 0
      do i = 1,nnza
         iaaug = iaaug + 1
         sub(isub)%i_aaug_sparse(iaaug) = sub(isub)%i_a_sparse(i) 
         sub(isub)%j_aaug_sparse(iaaug) = sub(isub)%j_a_sparse(i) 
         sub(isub)%aaug_sparse(iaaug)   = sub(isub)%a_sparse(i)   
      end do
      ! copy entries of right block of C^T with proper shift in columns
      do i = 1,nnzc
         iaaug = iaaug + 1
         icoli = sub(isub)%j_c_sparse(i)
         icol  = sub(isub)%iivsvn(icoli)
         sub(isub)%i_aaug_sparse(iaaug) = icol
         sub(isub)%j_aaug_sparse(iaaug) = sub(isub)%i_c_sparse(i) + ndof
         sub(isub)%aaug_sparse(iaaug)   = sub(isub)%c_sparse(i)   
      end do
      if      (sub(isub)%matrixtype .eq. 0) then
         ! unsymmetric case: apply lower block of C
         do i = 1,nnzc
            iaaug = iaaug + 1
            icoli = sub(isub)%j_c_sparse(i)
            icol  = sub(isub)%iivsvn(icoli)
            sub(isub)%i_aaug_sparse(iaaug) = sub(isub)%i_c_sparse(i) + ndof
            sub(isub)%j_aaug_sparse(iaaug) = icol
            sub(isub)%aaug_sparse(iaaug)   = sub(isub)%c_sparse(i)   
         end do
      end if
      if (iaaug.ne.nnzaaug) then
         write(*,*) 'DD_PREPARE_AUG: Actual length of augmented matrix does not match for subdomain:', isub
         call error_exit
      end if
      sub(isub)%nnzaaug = nnzaaug
      sub(isub)%laaug   = laaug


      ! factorize matrix Aaug
      sub(isub)%is_mumps_aug_active = .true.

      ! Set type of matrix
      if      (sub(isub)%matrixtype .eq. 0) then
         ! unsymmetric case:
         aaugmatrixtype = 0
      else if (sub(isub)%matrixtype .eq. 1 .or. sub(isub)%matrixtype .eq. 2) then
         ! in symmetric case, saddle point problem makes the augmented matrix indefinite,
         ! even if the original matrix is SPD:
         aaugmatrixtype = 2
      end if
      call mumps_init(sub(isub)%mumps_aug,comm,aaugmatrixtype)

      ! Verbosity level of MUMPS
      if (debug) then
         mumpsinfo = 2
      else
         mumpsinfo = 0
      end if
      call mumps_set_info(sub(isub)%mumps_aug,mumpsinfo)

      ! Load matrix to MUMPS
      ndof     = sub(isub)%ndof
      nconstr  = sub(isub)%nconstr
      ndofaaug = ndof + nconstr

      nnzaaug = sub(isub)%nnzaaug
      laaug   = sub(isub)%laaug
      call mumps_load_triplet(sub(isub)%mumps_aug,ndofaaug,nnzaaug,&
                                  sub(isub)%i_aaug_sparse,sub(isub)%j_aaug_sparse,sub(isub)%aaug_sparse,laaug)
      ! Analyze matrix
      call mumps_analyze(sub(isub)%mumps_aug) 
      ! Factorize matrix 
      call mumps_factorize(sub(isub)%mumps_aug) 

      sub(isub)%is_aug_factorized = .true.

end subroutine

!**************************************
subroutine dd_prepare_coarse(myid,isub)
!**************************************
! Subroutine for solving of system 
! | A C^T|| phis | = | 0 |
! ! C  0 ||lambda|   | I |
! phis are coarse space basis functions on subdomain
! Then the routine builds the local coarse matrix:
! Ac = phis^T * A * phis 

      use module_sm
      use module_mumps
      use module_utils
      implicit none

      ! processor ID
      integer,intent(in) :: myid
      ! subdomain
      integer,intent(in) :: isub

      ! local vars
      integer ::  ndof, nconstr, ndofaaug, ndofi, matrixtype
      integer ::  ndofc
      integer ::  i, j, indphis, nrhs, &
                  indphisstart, indi, icoarsem, lcoarsem
      integer ::  lphisi1, lphisi2

      integer ::             lphis
      real(kr),allocatable :: phis(:)

      integer ::             lac1, lac2
      real(kr),allocatable :: ac(:,:)

      ! check the prerequisities
      if (sub(isub)%proc .ne. myid) then
         if (debug) then
            write(*,*) 'DD_PREPARE_COARSE: myid =',myid,'Subdomain', isub,' is not mine.'
         end if
         return
      end if
      if (.not.sub(isub)%is_aug_factorized) then
         write(*,*) 'DD_PREPARE_COARSE: Augmented matrix in not factorized for subdomain:', isub
         call error_exit
      end if
      if (.not.sub(isub)%is_mumps_aug_active) then
         write(*,*) 'DD_PREPARE_COARSE: Augmented matrix solver in not ready for subdomain:', isub
         call error_exit
      end if

      ndof     = sub(isub)%ndof
      nconstr  = sub(isub)%nconstr
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
      call mumps_resolve(sub(isub)%mumps_aug,phis,lphis,nrhs)

      if (debug) then
         write(*,*) 'Subdomain ',isub,' coarse basis functions phis:'
         do i = 1,ndofaaug
            write(*,'(100f13.6)') (phis((j-1)*ndofaaug + i),j = 1,nconstr)
         end do
      end if

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
      if (debug) then
         write(*,*) 'Subdomain ',isub,' coarse matrix:'
         do i = 1,nconstr
            write(*,'(100f13.6)') (ac(i,j),j = 1,nconstr)
         end do
      end if

      ! restrict vector phis to interface unknowns and load it to the structure
      ndofi   = sub(isub)%ndofi
      lphisi1 = ndofi
      lphisi2 = nconstr
      allocate(sub(isub)%phisi(lphisi1,lphisi2))
      sub(isub)%lphisi1 = lphisi1
      sub(isub)%lphisi2 = lphisi2
      do i = 1,ndofi
         indi = sub(isub)%iivsvn(i)
         do j = 1,nconstr
            indphis = (j-1)*ndofaaug + indi

            sub(isub)%phisi(i,j) = phis(indphis)
         end do
      end do
      sub(isub)%is_phisi_prepared   = .true.

      ! load the coarse matrix to the structure in appropriate format
      matrixtype = sub(isub)%matrixtype
      if      (matrixtype.eq.0) then
         ! in unsymmetric case, load the whole coarse matrix columnwise
         lcoarsem = nconstr * nconstr
         allocate(sub(isub)%coarsem(lcoarsem))
         icoarsem = 0
         do j = 1,nconstr
            do i = 1,nconstr
               icoarsem = icoarsem + 1
               sub(isub)%coarsem(icoarsem) = ac(i,j)
            end do
         end do
         sub(isub)%lcoarsem = lcoarsem
      else if (matrixtype.eq.1 .or. matrixtype.eq.2) then
         ! in symmetric case, load the upper triangle columnwise
         lcoarsem = (nconstr+1)*nconstr/2
         allocate(sub(isub)%coarsem(lcoarsem))
         icoarsem = 0
         do j = 1,nconstr
            do i = 1,j
               icoarsem = icoarsem + 1
               sub(isub)%coarsem(icoarsem) = ac(i,j)
            end do
         end do
         sub(isub)%lcoarsem = lcoarsem
      end if
      ! check the counter
      if (icoarsem.ne.lcoarsem) then
         write(*,*) 'DD_PREPARE_COARSE: Check of coarse matrix length failed for subdomain: ',isub
         call error_exit
      end if

      ! number of coarse degrees of freedom equals number of constraints in this implementation
      ndofc = nconstr
      sub(isub)%ndofc = ndofc

      sub(isub)%is_coarse_prepared = .true.

      deallocate(ac)
      deallocate(phis)
end subroutine

!***************************************************
subroutine dd_multiply_by_schur(myid,isub,x,lx,y,ly)
!***************************************************
! Subroutine for multiplication of interface vector by Schur complement
      use module_utils
      use module_mumps
      use module_sm
      implicit none

      ! processor ID
      integer,intent(in) :: myid
      ! subdomain number
      integer,intent(in) :: isub

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

      if (sub(isub)%proc .eq. myid) then

         ! check the prerequisities
         if (.not.sub(isub)%is_matrix_loaded) then
            write(*,*) 'DD_PREPARE_SCHUR: Matrix is not loaded for subdomain:', isub
            call error_exit
         end if
         if (.not. (sub(isub)%is_blocked)) then
            write(*,*) 'DD_PREPARE_SCHUR: Matrix is not in blocked format. Call routine to do this.'
            call error_exit
         end if
         if (.not. (sub(isub)%is_interior_factorized)) then
            write(*,*) 'DD_PREPARE_SCHUR: Interior block not factorized yet.'
            call error_exit
         end if
 
         ! prepare rhs vector for backsubstitution to problem A_11*aux1 = -A_12*x
         ndofo = sub(isub)%ndofo
         laux1 = ndofo
         allocate(aux1(laux1))

         ! with offdiagonal blocks, use as nonsymmetric
         matrixtype_aux = 0
         nnza12     = sub(isub)%nnza12
         la12       = sub(isub)%la12
         call sm_vec_mult(matrixtype_aux, nnza12, &
                          sub(isub)%i_a12_sparse, sub(isub)%j_a12_sparse, sub(isub)%a12_sparse, la12, &
                          x,lx, aux1,laux1)

         ! resolve interior problem by MUMPS
         call mumps_resolve(sub(isub)%mumps_interior_block,aux1,laux1)
         
         if (sub(isub)%istorage .eq. 4) then
            is_symmetric_storage = .true.
         else
            is_symmetric_storage = .false.
         end if

         ! prepare auxiliary vector for multiplication
         ndofi = sub(isub)%ndofi
         laux2 = ndofi
         allocate(aux2(laux2))

         ! get aux2 = A_21*aux1, i.e. aux2 = A_21 * (A_11)^-1 * A_12 * x
         if (is_symmetric_storage) then
            matrixtype_aux = 0
            nnza12     = sub(isub)%nnza12
            la12       = sub(isub)%la12
            ! use the matrix with transposed indices in the call sm_vec_mult
            call sm_vec_mult(matrixtype_aux, nnza12, &
                             sub(isub)%j_a12_sparse, sub(isub)%i_a12_sparse, sub(isub)%a12_sparse, la12, &
                             aux1,laux1, aux2,laux2)
         else
            matrixtype_aux = 0
            nnza21     = sub(isub)%nnza21
            la21       = sub(isub)%la21
            call sm_vec_mult(matrixtype_aux, nnza21, &
                             sub(isub)%i_a21_sparse, sub(isub)%j_a21_sparse, sub(isub)%a21_sparse, la21, &
                             aux1,laux1, aux2,laux2)
         end if

         ! get y = A_22*x
         matrixtype = sub(isub)%matrixtype
         nnza22     = sub(isub)%nnza22
         la22       = sub(isub)%la22
         call sm_vec_mult(matrixtype, nnza22, &
                          sub(isub)%i_a22_sparse, sub(isub)%j_a22_sparse, sub(isub)%a22_sparse, la22, &
                          x,lx, y,ly)

         ! add results together to get y = y - aux2, i.e. y = A_22 * x - A_21 * (A_11)^-1 * A_12 * x, or y = (A_22 - A_21 * (A_11)^-1 * A_12) * x
         do i = 1,ly
            y(i) = y(i) - aux2(i)
         end do
 
         deallocate(aux1)
         deallocate(aux2)
      end if

end subroutine

!********************************************
subroutine dd_where_is_subdomain(isub,idproc)
!********************************************
! Subroutine for locating subdomain data
      implicit none
! Index of subdomain whose data I want to locate
      integer,intent(in) :: isub
! Index of processor where the data are located (if errors, idproc < 0)
      integer,intent(out) :: idproc

      if (.not.allocated(sub).or.isub.gt.lsub) then
         write(*,*) 'DD_WHERE_IS_SUBDOMAIN: Trying to localize nonexistent subdomain.'
         idproc = -1
         return
      end if
      if (.not.sub(isub)%is_sub_identified) then
         write(*,*) 'DD_WHERE_IS_SUBDOMAIN: Subdomain is not identified.'
         idproc = -2
         return
      end if
      if (sub(isub)%isub .ne. isub) then
         write(*,*) 'DD_WHERE_IS_SUBDOMAIN: Subdomain has strange number.'
         idproc = -3
         return
      end if
      if (.not.sub(isub)%is_proc_assigned) then
         write(*,*) 'DD_WHERE_IS_SUBDOMAIN: Processor is not assigned.'
         idproc = -4
         return
      end if

      ! if all checks are OK, return processor number
      idproc = sub(isub)%proc

end subroutine

!******************************************************
subroutine dd_get_interface_size(myid,isub,ndofi,nnodi)
!******************************************************
! Subroutine for finding size of subdomain data
      implicit none
! processor ID
      integer,intent(in) :: myid
! Index of subdomain whose data I want to get
      integer,intent(in) :: isub
! Length of vector of interface dof
      integer,intent(out) :: ndofi
! Length of vector of interface nodes
      integer,intent(out) :: nnodi

      if (.not.allocated(sub).or.isub.gt.lsub) then
         write(*,*) 'DD_GET_INTERFACE_SIZE: Trying to localize nonexistent subdomain.'
         ndofi = -1
         return
      end if
      if (.not.sub(isub)%is_sub_identified) then
         write(*,*) 'DD_GET_INTERFACE_SIZE: Subdomain is not identified.'
         ndofi = -2
         return
      end if
      if (sub(isub)%isub .ne. isub) then
         write(*,*) 'DD_GET_INTERFACE_SIZE: Subdomain has strange number.'
         ndofi = -3
         return
      end if
      if (.not.sub(isub)%is_proc_assigned) then
         write(*,*) 'DD_GET_INTERFACE_SIZE: Processor is not assigned.'
         ndofi = -4
         return
      end if
      if (sub(isub)%proc .ne. myid) then
         write(*,*) 'DD_GET_INTERFACE_SIZE: Subdomain ',isub,' is not mine, myid = ',myid
         ndofi = -5
         return
      end if
      if (.not.sub(isub)%is_mesh_loaded) then
         write(*,*) 'DD_GET_INTERFACE_SIZE: Mesh is not loaded yet for subdomain ',isub
         ndofi = -6
         return
      end if

      ! if all checks are OK, return subdomain interface size
      ndofi = sub(isub)%ndofi
      nnodi = sub(isub)%nnodi

end subroutine

!****************************************************
subroutine dd_get_number_of_crows(myid,isub,lindrowc)
!****************************************************
! Subroutine for finding size of subdomain data
      implicit none
! processor ID
      integer,intent(in) :: myid
! Index of subdomain whose data I want to get
      integer,intent(in) :: isub
! Number of corners
      integer,intent(out) :: lindrowc

      if (.not.allocated(sub).or.isub.gt.lsub) then
         write(*,*) 'DD_GET_NUMBER_OF_CROWS: Trying to localize nonexistent subdomain.'
         lindrowc = -1
         return
      end if
      if (.not.sub(isub)%is_sub_identified) then
         write(*,*) 'DD_GET_NUMBER_OF_CROWS: Subdomain is not identified.'
         lindrowc = -2
         return
      end if
      if (sub(isub)%isub .ne. isub) then
         write(*,*) 'DD_GET_NUMBER_OF_CROWS: Subdomain has strange number.'
         lindrowc = -3
         return
      end if
      if (.not.sub(isub)%is_proc_assigned) then
         write(*,*) 'DD_GET_NUMBER_OF_CROWS: Processor is not assigned.'
         lindrowc = -4
         return
      end if
      if (sub(isub)%proc .ne. myid) then
         write(*,*) 'DD_GET_NUMBER_OF_CROWS: Subdomain ',isub,' is not mine, myid = ',myid
         lindrowc = -5
         return
      end if
      if (.not.sub(isub)%is_mesh_loaded) then
         write(*,*) 'DD_GET_NUMBER_OF_CROWS: Mesh is not loaded yet for subdomain ',isub
         lindrowc = -6
         return
      end if
      if (.not.sub(isub)%is_cnodes_loaded) then
         write(*,*) 'DD_GET_NUMBER_OF_CROWS: Coarse nodes not loaded yet for subdomain ',isub
         lindrowc = -6
         return
      end if

      ! if all checks are OK, return subdomain interface size
      lindrowc = sub(isub)%lindrowc

end subroutine

!***********************************************
subroutine dd_get_number_of_cnnz(myid,isub,nnzc)
!***********************************************
! Subroutine for finding size of subdomain data
      implicit none
! processor ID
      integer,intent(in) :: myid
! Index of subdomain whose data I want to get
      integer,intent(in) :: isub
! Number of corners
      integer,intent(out) :: nnzc

      if (.not.allocated(sub).or.isub.gt.lsub) then
         write(*,*) 'DD_GET_NUMBER_OF_CNNZ: Trying to localize nonexistent subdomain.'
         nnzc = -1
         return
      end if
      if (.not.sub(isub)%is_sub_identified) then
         write(*,*) 'DD_GET_NUMBER_OF_CNNZ: Subdomain is not identified.'
         nnzc = -2
         return
      end if
      if (sub(isub)%isub .ne. isub) then
         write(*,*) 'DD_GET_NUMBER_OF_CNNZ: Subdomain has strange number.'
         nnzc = -3
         return
      end if
      if (.not.sub(isub)%is_proc_assigned) then
         write(*,*) 'DD_GET_NUMBER_OF_CNNZ: Processor is not assigned.'
         nnzc = -4
         return
      end if
      if (sub(isub)%proc .ne. myid) then
         write(*,*) 'DD_GET_NUMBER_OF_CNNZ: Subdomain ',isub,' is not mine, myid = ',myid
         nnzc = -5
         return
      end if
      if (.not.sub(isub)%is_mesh_loaded) then
         write(*,*) 'DD_GET_NUMBER_OF_CNNZ: Mesh is not loaded yet for subdomain ',isub
         nnzc = -6
         return
      end if
      if (.not.sub(isub)%is_c_loaded) then
         write(*,*) 'DD_GET_NUMBER_OF_CNNZ: Matrix C not loaded yet for subdomain ',isub
         nnzc = -6
         return
      end if

      ! if all checks are OK, return subdomain interface size
      nnzc = sub(isub)%nnzc

end subroutine

!************************************************************
subroutine dd_get_subdomain_crows(myid,isub,indrowc,lindrowc)
!************************************************************
! Subroutine for getting corner nodes
      implicit none
! processor ID
      integer,intent(in) :: myid
! Index of subdomain whose data I want to get
      integer,intent(in) :: isub
! Corners
      integer,intent(in) :: lindrowc
      integer,intent(out) :: indrowc(lindrowc)
! local vars
      integer :: i

      if (.not.allocated(sub).or.isub.gt.lsub) then
         write(*,*) 'DD_GET_SUBDOMAIN_CROWS: Trying to localize nonexistent subdomain.'
         return
      end if
      if (.not.sub(isub)%is_sub_identified) then
         write(*,*) 'DD_GET_SUBDOMAIN_CROWS: Subdomain is not identified.'
         return
      end if
      if (sub(isub)%isub .ne. isub) then
         write(*,*) 'DD_GET_SUBDOMAIN_CROWS: Subdomain has strange number.'
         return
      end if
      if (.not.sub(isub)%is_proc_assigned) then
         write(*,*) 'DD_GET_SUBDOMAIN_CROWS: Processor is not assigned.'
         return
      end if
      if (sub(isub)%proc .ne. myid) then
         write(*,*) 'DD_GET_SUBDOMAIN_CROWS: Subdomain ',isub,' is not mine, myid = ',myid
         return
      end if
      if (.not.sub(isub)%is_mesh_loaded) then
         write(*,*) 'DD_GET_SUBDOMAIN_CROWS: Mesh is not loaded yet for subdomain ',isub
         return
      end if
      if (sub(isub)%lindrowc .ne. lindrowc) then
         write(*,*) 'DD_GET_SUBDOMAIN_CROWS: Size of array for corners not consistent.'
         write(*,*) 'lindrowc :', sub(isub)%lindrowc, 'array size: ',lindrowc
         return
      end if

      ! if all checks are OK, return number of cnodes on subdomain
      do i = 1,lindrowc
         indrowc(i) = sub(isub)%indrowc(i)
      end do

end subroutine

!*********************************************************
subroutine dd_get_subdomain_corner_number(myid,isub,nnodc)
!*********************************************************
! Subroutine for getting corner nodes
      implicit none
! processor ID
      integer,intent(in) :: myid
! Index of subdomain whose data I want to get
      integer,intent(in) :: isub
! Corners
      integer,intent(out) :: nnodc

      if (.not.allocated(sub).or.isub.gt.lsub) then
         write(*,*) 'DD_GET_SUBDOMAIN_CONRER_NUMBER: Trying to localize nonexistent subdomain.'
         return
      end if
      if (.not.sub(isub)%is_sub_identified) then
         write(*,*) 'DD_GET_SUBDOMAIN_CONRER_NUMBER: Subdomain is not identified.'
         return
      end if
      if (sub(isub)%isub .ne. isub) then
         write(*,*) 'DD_GET_SUBDOMAIN_CONRER_NUMBER: Subdomain has strange number.'
         return
      end if
      if (.not.sub(isub)%is_proc_assigned) then
         write(*,*) 'DD_GET_SUBDOMAIN_CONRER_NUMBER: Processor is not assigned.'
         return
      end if
      if (sub(isub)%proc .ne. myid) then
         write(*,*) 'DD_GET_SUBDOMAIN_CONRER_NUMBER: Subdomain ',isub,' is not mine, myid = ',myid
         return
      end if
      if (.not.sub(isub)%is_mesh_loaded) then
         write(*,*) 'DD_GET_SUBDOMAIN_CONRER_NUMBER: Mesh is not loaded yet for subdomain ',isub
         return
      end if

      ! if all checks are OK, return number of corners on subdomain
      nnodc = sub(isub)%nnodc

end subroutine

!**************************************************************************
subroutine dd_get_subdomain_cmap(myid,isub,i_c_sparse,j_c_sparse,lc_sparse)
!**************************************************************************
! Subroutine for getting corner nodes embedding into interface
      implicit none
! processor ID
      integer,intent(in) :: myid
! Index of subdomain whose data I want to get
      integer,intent(in) :: isub
! Corners mapping
      integer,intent(in) ::  lc_sparse
      integer,intent(out) :: i_c_sparse(lc_sparse), j_c_sparse(lc_sparse)

! local vars
      integer :: nnzc, i

      if (.not.allocated(sub).or.isub.gt.lsub) then
         write(*,*) 'DD_GET_SUBDOMAIN_CMAP: Trying to localize nonexistent subdomain.'
         return
      end if
      if (.not.sub(isub)%is_sub_identified) then
         write(*,*) 'DD_GET_SUBDOMAIN_CMAP: Subdomain is not identified.'
         return
      end if
      if (sub(isub)%isub .ne. isub) then
         write(*,*) 'DD_GET_SUBDOMAIN_CMAP: Subdomain has strange number.'
         return
      end if
      if (.not.sub(isub)%is_proc_assigned) then
         write(*,*) 'DD_GET_SUBDOMAIN_CMAP: Processor is not assigned.'
         return
      end if
      if (sub(isub)%proc .ne. myid) then
         write(*,*) 'DD_GET_SUBDOMAIN_CMAP: Subdomain ',isub,' is not mine, myid = ',myid
         return
      end if
      if (.not.sub(isub)%is_mesh_loaded) then
         write(*,*) 'DD_GET_SUBDOMAIN_CMAP: Mesh is not loaded yet for subdomain ',isub
         return
      end if
      if (.not.sub(isub)%is_c_loaded) then
         write(*,*) 'DD_GET_SUBDOMAIN_CMAP: Matrix C is not loaded yet for subdomain ',isub
         return
      end if
      if (sub(isub)%lc .ne. lc_sparse) then
         write(*,*) 'DD_GET_SUBDOMAIN_CMAP: Size of array for corners not consistent.'
         write(*,*) 'nnodc :', sub(isub)%lc, 'array size: ',lc_sparse 
         return
      end if

      ! if all checks are OK, return subdomain interface size
      nnzc = sub(isub)%nnzc
      do i = 1,nnzc
         i_c_sparse(i) = sub(isub)%i_c_sparse(i)
      end do
      do i = 1,nnzc
         j_c_sparse(i) = sub(isub)%j_c_sparse(i)
      end do

end subroutine

!*****************************************************************
subroutine dd_get_subdomain_interface_nndf(myid,isub,nndfi,lnndfi)
!*****************************************************************
! Subroutine for getting subdomain nndfi array
      implicit none
! processor ID
      integer,intent(in) :: myid
! Index of subdomain whose data I want to get
      integer,intent(in) :: isub
! Array of numberf of DOF at interface nodes
      integer,intent(in) :: lnndfi
      integer,intent(out) :: nndfi(lnndfi)
! local vars
      integer :: nnodi, i, indnod

      if (.not.allocated(sub).or.isub.gt.lsub) then
         write(*,*) 'DD_GET_SUBDOMAIN_INTERFACE_NNDF: Trying to localize nonexistent subdomain.'
         return
      end if
      if (.not.sub(isub)%is_sub_identified) then
         write(*,*) 'DD_GET_SUBDOMAIN_INTERFACE_NNDF: Subdomain is not identified.'
         return
      end if
      if (sub(isub)%isub .ne. isub) then
         write(*,*) 'DD_GET_SUBDOMAIN_INTERFACE_NNDF: Subdomain has strange number.'
         return
      end if
      if (.not.sub(isub)%is_proc_assigned) then
         write(*,*) 'DD_GET_SUBDOMAIN_INTERFACE_NNDF: Processor is not assigned.'
         return
      end if
      if (sub(isub)%proc .ne. myid) then
         write(*,*) 'DD_GET_SUBDOMAIN_INTERFACE_NNDF: Subdomain ',isub,' is not mine, myid = ',myid
         return
      end if
      if (.not.sub(isub)%is_mesh_loaded) then
         write(*,*) 'DD_GET_SUBDOMAIN_INTERFACE_NNDF: Mesh is not loaded yet for subdomain ',isub
         return
      end if
      if (.not.sub(isub)%is_mesh_loaded) then
         write(*,*) 'DD_GET_SUBDOMAIN_INTERFACE_NNDF: Mesh is not loaded yet for subdomain ',isub
         return
      end if
      if (sub(isub)%nnodi .ne. lnndfi) then
         write(*,*) 'DD_GET_SUBDOMAIN_INTERFACE_NNDF: Size of array for interface not consistent.'
         write(*,*) 'nnodc :', sub(isub)%nnodi, 'array size: ',lnndfi 
         return
      end if

      ! if all checks are OK, construct array of numbers of DOF at interface nodes
      nnodi = sub(isub)%nnodi
      do i = 1,nnodi
         indnod   = sub(isub)%iin(i)
         nndfi(i) = sub(isub)%nndf(indnod)
      end do

end subroutine

!*******************************************************************
subroutine dd_get_interface_global_numbers(myid, isub, iingn,liingn)
!*******************************************************************
! Subroutine for getting mapping of interface nodes into global node numbering
      use module_utils
      implicit none

      integer,intent(in) :: myid, isub
      integer,intent(in)  :: liingn
      integer,intent(out) ::  iingn(liingn)

      ! local vars
      integer :: nnodi
      integer :: i, ind

      ! check if I store the subdomain
      if (.not. sub(isub)%proc .eq. myid) then
         if (debug) then
            write(*,*) 'DD_GET_INTERFACE_GLOBAL_NUMBERS: myid =',myid,', not my subdomain: ',isub
         end if
         return
      end if
      ! check if mesh is loaded
      if (.not. sub(isub)%is_mesh_loaded) then
         if (debug) then
            write(*,*) 'DD_GET_INTERFACE_GLOBAL_NUMBERS: myid =',myid,', Mesh not loaded for subdomain: ',isub
         end if
         call error_exit
      end if
      ! check dimensions
      if (sub(isub)%nnodi.ne.liingn) then
         if (debug) then
            write(*,*) 'DD_GET_INTERFACE_GLOBAL_NUMBERS: myid =',myid,', Interface dimensions mismatch for subdomain ',isub
         end if
         call error_exit
      end if

      ! load data
      nnodi = sub(isub)%nnodi
      do i = 1,nnodi
         ! subdomain node number
         ind = sub(isub)%iin(i)
         ! global node number
         iingn(i) = sub(isub)%isngn(ind)
      end do
end subroutine

!************************************************
subroutine dd_get_schur(myid, isub, schur,lschur)
!************************************************
! Subroutine for explicit construction of Schur complement
      use module_utils
      implicit none

      integer,intent(in) :: myid, isub
      integer,intent(in)  ::  lschur
      real(kr),intent(out) ::  schur(lschur)

      ! local vars
      integer :: j, lmat, ndofi, pointschur

      integer              :: le
      real(kr),allocatable ::  e(:)

      ! check if I store the subdomain
      if (.not. sub(isub)%proc .eq. myid) then
         if (debug) then
            write(*,*) 'DD_GET_SCHUR: myid =',myid,', not my subdomain: ',isub
         end if
         return
      end if
      ! check if mesh is loaded
      if (.not. sub(isub)%is_mesh_loaded) then
         if (debug) then
            write(*,*) 'DD_GET_SCHUR: myid =',myid,', Mesh not loaded for subdomain: ',isub
         end if
         call error_exit
      end if
      ! check dimensions
      ndofi = sub(isub)%ndofi
      lmat = ndofi*ndofi
      if (lmat.ne.lschur) then
         if (debug) then
            write(*,*) 'DD_GET_SCHUR: myid =',myid,', Interface dimensions mismatch for subdomain ',isub
         end if
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
         call dd_multiply_by_schur(myid,isub,e,le,schur(pointschur),le)

         pointschur = pointschur + ndofi
      end do
      deallocate(e)

end subroutine

!***********************************************************
subroutine dd_get_interface_diagonal(myid, isub, rhoi,lrhoi)
!***********************************************************
! Subroutine for getting subdomain diagonal from the structure
      use module_utils
      implicit none

      integer,intent(in) :: myid, isub
      integer,intent(in)  :: lrhoi
      real(kr),intent(out) ::  rhoi(lrhoi)

      ! local vars
      integer :: ndof, ndofi, la
      integer :: i, indiv, ia, ind

      integer ::              lrho
      real(kr),allocatable ::  rho(:)

      ! check if I store the subdomain
      if (.not. sub(isub)%proc .eq. myid) then
         if (debug) then
            write(*,*) 'DD_GET_INTERFACE_DIAGONAL: myid =',myid,', not my subdomain: ',isub
         end if
         return
      end if
      ! check if matrix is loaded
      if (.not. sub(isub)%is_matrix_loaded) then
         if (debug) then
            write(*,*) 'DD_GET_INTERFACE_DIAGONAL: myid =',myid,', Matrix not loaded for subdomain: ',isub
         end if
         call error_exit
      end if
      ! check dimensions
      if (sub(isub)%ndofi.ne.lrhoi) then
         if (debug) then
            write(*,*) 'DD_GET_INTERFACE_DIAGONAL: myid =',myid,', Interface dimensions mismatch for subdomain ',isub
         end if
         call error_exit
      end if

      ! load data
      ndof  = sub(isub)%ndof
      ndofi = sub(isub)%ndofi

      ! allocate vector for whole subdomain diagonal
      lrho = ndof
      allocate(rho(lrho))
      call zero(rho,lrho)

      la = sub(isub)%la
      do ia = 1,la
         if (sub(isub)%i_a_sparse(ia).eq.sub(isub)%j_a_sparse(ia)) then
            ind = sub(isub)%i_a_sparse(ia)

            rho(ind) = rho(ind) + sub(isub)%a_sparse(ia)
         end if
      end do

      do i = 1,ndofi
         indiv = sub(isub)%iivsvn(i)
         rhoi(i) = rho(indiv)
      end do

      deallocate(rho)
end subroutine

!**********************************
subroutine dd_clear_subdomain(isub)
!**********************************
! Subroutine for deallocation of data of subdomain isub
      implicit none
! Index of subdomain whose data I want to deallocate
      integer,intent(in) :: isub

      ! local variables
      integer :: j

      if (.not.allocated(sub).or.isub.gt.lsub) then
         write(*,*) 'DD_CLEAR_SUBDOMAIN: Trying to clear nonexistent subdomain.'
         return
      end if

      if (allocated(sub(isub)%inet)) then
         deallocate(sub(isub)%inet)
      end if
      if (allocated(sub(isub)%nnet)) then
         deallocate(sub(isub)%nnet)
      end if
      if (allocated(sub(isub)%nndf)) then
         deallocate(sub(isub)%nndf)
      end if
      if (allocated(sub(isub)%isngn)) then
         deallocate(sub(isub)%isngn)
      end if
      if (allocated(sub(isub)%xyz)) then
         deallocate(sub(isub)%xyz)
      end if
      if (allocated(sub(isub)%iin)) then
         deallocate(sub(isub)%iin)
      end if
      if (allocated(sub(isub)%iivsvn)) then
         deallocate(sub(isub)%iivsvn)
      end if
      if (allocated(sub(isub)%iovsvn)) then
         deallocate(sub(isub)%iovsvn)
      end if
      if (allocated(sub(isub)%ifix)) then
         deallocate(sub(isub)%ifix)
      end if
      if (allocated(sub(isub)%fixv)) then
         deallocate(sub(isub)%fixv)
      end if
      if (allocated(sub(isub)%bc)) then
         deallocate(sub(isub)%bc)
      end if
      ! corners and globs
      if (allocated(sub(isub)%global_corner_number)) then
         deallocate(sub(isub)%global_corner_number)
      end if
      if (allocated(sub(isub)%icnsin)) then
         deallocate(sub(isub)%icnsin)
      end if
      if (allocated(sub(isub)%ncdf)) then
         deallocate(sub(isub)%ncdf)
      end if
      if (allocated(sub(isub)%global_glob_number)) then
         deallocate(sub(isub)%global_glob_number)
      end if
      if (allocated(sub(isub)%nglobvar)) then
         deallocate(sub(isub)%nglobvar)
      end if
      if (allocated(sub(isub)%igvsivn)) then
         deallocate(sub(isub)%igvsivn)
      end if
      if (allocated(sub(isub)%glob_type)) then
         deallocate(sub(isub)%glob_type)
      end if
      if (allocated(sub(isub)%ngdf)) then
         deallocate(sub(isub)%ngdf)
      end if
      if (allocated(sub(isub)%cnodes)) then
         do j = 1,sub(isub)%ncnodes
            if (allocated(sub(isub)%cnodes(j)%xyz)) then
               deallocate(sub(isub)%cnodes(j)%xyz)
            end if
            if (allocated(sub(isub)%cnodes(j)%igcdof)) then
               deallocate(sub(isub)%cnodes(j)%igcdof)
            end if
            if (allocated(sub(isub)%cnodes(j)%insin)) then
               deallocate(sub(isub)%cnodes(j)%insin)
            end if
            if (allocated(sub(isub)%cnodes(j)%ivsivn)) then
               deallocate(sub(isub)%cnodes(j)%ivsivn)
            end if
            if (allocated(sub(isub)%cnodes(j)%matrix)) then
               deallocate(sub(isub)%cnodes(j)%matrix)
            end if
         end do
         deallocate(sub(isub)%cnodes)
      end if
      if (allocated(sub(isub)%i_a_sparse)) then
         deallocate(sub(isub)%i_a_sparse)
      end if
      if (allocated(sub(isub)%j_a_sparse)) then
         deallocate(sub(isub)%j_a_sparse)
      end if
      if (allocated(sub(isub)%a_sparse)) then
         deallocate(sub(isub)%a_sparse)
      end if
      if (allocated(sub(isub)%i_a11_sparse)) then
         deallocate(sub(isub)%i_a11_sparse)
      end if
      if (allocated(sub(isub)%j_a11_sparse)) then
         deallocate(sub(isub)%j_a11_sparse)
      end if
      if (allocated(sub(isub)%a11_sparse)) then
         deallocate(sub(isub)%a11_sparse)
      end if
      if (allocated(sub(isub)%i_a12_sparse)) then
         deallocate(sub(isub)%i_a12_sparse)
      end if
      if (allocated(sub(isub)%j_a12_sparse)) then
         deallocate(sub(isub)%j_a12_sparse)
      end if
      if (allocated(sub(isub)%a12_sparse)) then
         deallocate(sub(isub)%a12_sparse)
      end if
      if (allocated(sub(isub)%i_a21_sparse)) then
         deallocate(sub(isub)%i_a21_sparse)
      end if
      if (allocated(sub(isub)%j_a21_sparse)) then
         deallocate(sub(isub)%j_a21_sparse)
      end if
      if (allocated(sub(isub)%a21_sparse)) then
         deallocate(sub(isub)%a21_sparse)
      end if
      if (allocated(sub(isub)%i_a22_sparse)) then
         deallocate(sub(isub)%i_a22_sparse)
      end if
      if (allocated(sub(isub)%j_a22_sparse)) then
         deallocate(sub(isub)%j_a22_sparse)
      end if
      if (allocated(sub(isub)%a22_sparse)) then
         deallocate(sub(isub)%a22_sparse)
      end if
      if (sub(isub)%is_mumps_interior_active) then
         call mumps_finalize(sub(isub)%mumps_interior_block)
      end if

      ! BDDC matrices
      if (allocated(sub(isub)%i_c_sparse)) then
         deallocate(sub(isub)%i_c_sparse)
      end if
      if (allocated(sub(isub)%j_c_sparse)) then
         deallocate(sub(isub)%j_c_sparse)
      end if
      if (allocated(sub(isub)%c_sparse)) then
         deallocate(sub(isub)%c_sparse)
      end if
      if (allocated(sub(isub)%indrowc)) then
         deallocate(sub(isub)%indrowc)
      end if
      if (allocated(sub(isub)%i_aaug_sparse)) then
         deallocate(sub(isub)%i_aaug_sparse)
      end if
      if (allocated(sub(isub)%j_aaug_sparse)) then
         deallocate(sub(isub)%j_aaug_sparse)
      end if
      if (allocated(sub(isub)%aaug_sparse)) then
         deallocate(sub(isub)%aaug_sparse)
      end if
      if (allocated(sub(isub)%phisi)) then
         deallocate(sub(isub)%phisi)
      end if
      if (allocated(sub(isub)%coarsem)) then
         deallocate(sub(isub)%coarsem)
      end if

      if (sub(isub)%is_mumps_aug_active) then
         call mumps_finalize(sub(isub)%mumps_aug)
      end if

      sub(isub)%is_mesh_loaded          = .false.
      sub(isub)%is_matrix_loaded        = .false.
      sub(isub)%is_c_loaded             = .false.
      sub(isub)%is_phisi_prepared       = .false.
      sub(isub)%is_coarse_prepared      = .false.
      sub(isub)%is_aug_factorized       = .false.
      sub(isub)%is_interior_factorized  = .false.

end subroutine

!****************************
subroutine dd_print_sub(myid)
!****************************
! Subroutine for printing the state of sub array
      use module_sm
      implicit none

      integer,intent(in) :: myid

! local variables
      integer :: isub, i, j

! basic structure
      write(*,*)    '******************************'
      write(*,*)    'I am processor                ', myid
      write(*,*)    '   lsub:                      ', lsub
      write(*,*)    '   Is sub allocated?          ', allocated(sub)
      do isub = 1,lsub
         write(*,*) '****** start subdomain export '
         write(*,*) ' local subdomain number :     ', isub
         write(*,*) '*** HEADER INFO :             '
         write(*,*) '     subdomain identified:    ', sub(isub)%is_sub_identified
         write(*,*) '     global subdomain number: ', sub(isub)%isub
         write(*,*) '     number of subdomains:    ', sub(isub)%nsub
         write(*,*) '*** PROCESSOR INFO :          '
         write(*,*) '     processor assigned:      ', sub(isub)%is_proc_assigned
         write(*,*) '     processor number:        ', sub(isub)%proc
         if (sub(isub)%proc .ne. myid) then
            cycle
         end if
         write(*,*) '*** MESH INFO :               '
         write(*,*) '     mesh loaded:             ', sub(isub)%is_mesh_loaded
         write(*,*) '     number of elements:      ', sub(isub)%nelem
         write(*,*) '     number of nodes:         ', sub(isub)%nnod
         write(*,*) '     number of DOF:           ', sub(isub)%ndof
         write(*,*) '     number of dimensions:    ', sub(isub)%ndim
         write(*,*) '*** BOUNDARY CONDITIONS :     '
         write(*,*) '     is bc present:           ', sub(isub)%is_bc_present
         write(*,*) '     is bc nonzero:           ', sub(isub)%is_bc_nonzero
         write(*,*) '*** CORNER INFO :             '
         write(*,*) '     number of corners:       ', sub(isub)%nnodc
         write(*,*) '*** GLOB INFO :               '
         write(*,*) '     number of globs:         ', sub(isub)%nglob
         write(*,*) '*** COARSE NODES INFO :       '
         write(*,*) '     number of coarse nodes:  ', sub(isub)%ncnodes
         do i = 1,sub(isub)%ncnodes
            call dd_print_cnode(isub,i)
         end do
         write(*,*) '*** MATRIX INFO :             '
         write(*,*) '     matrix loaded:           ', sub(isub)%is_matrix_loaded
         write(*,*) '     matrix blocked:          ', sub(isub)%is_blocked
         write(*,*) '     interior block factor.:  ', sub(isub)%is_interior_factorized
         if (debug) then
            if (sub(isub)%is_matrix_loaded.and.sub(isub)%is_triplet) then
               write(*,*) '     matrix data:           '
               call sm_print(6, sub(isub)%i_a_sparse, sub(isub)%j_a_sparse, sub(isub)%a_sparse, &
                             sub(isub)%la, sub(isub)%nnza)
            end if
            if (sub(isub)%is_matrix_loaded.and.sub(isub)%is_blocked) then
               write(*,*) '     matrix blocks:           '
               write(*,*) '     A_11:                     '
               call sm_print(6, sub(isub)%i_a11_sparse, sub(isub)%j_a11_sparse, sub(isub)%a11_sparse, &
                             sub(isub)%la11, sub(isub)%nnza11)
               write(*,*) '     A_12:                     '
               call sm_print(6, sub(isub)%i_a12_sparse, sub(isub)%j_a12_sparse, sub(isub)%a12_sparse, &
                             sub(isub)%la12, sub(isub)%nnza12)
               write(*,*) '     A_21:                     '
               call sm_print(6, sub(isub)%i_a21_sparse, sub(isub)%j_a21_sparse, sub(isub)%a21_sparse, &
                             sub(isub)%la21, sub(isub)%nnza21)
               write(*,*) '     A_22:                     '
               call sm_print(6, sub(isub)%i_a22_sparse, sub(isub)%j_a22_sparse, sub(isub)%a22_sparse, &
                             sub(isub)%la22, sub(isub)%nnza22)
            end if
         end if
         write(*,*) '*** BDDC INFO:                '
         write(*,*) '     matrix C loaded:         ', sub(isub)%is_c_loaded
         if (debug) then
            if (sub(isub)%is_c_loaded) then
               call sm_print(6, sub(isub)%i_c_sparse, sub(isub)%j_c_sparse, sub(isub)%c_sparse, &
                             sub(isub)%lc, sub(isub)%nnzc)
            end if
         end if
         write(*,*) '     matrix Kaug factorized:  ', sub(isub)%is_aug_factorized
         if (debug) then
            if (sub(isub)%is_matrix_loaded.and.sub(isub)%is_aug_factorized) then
               if (debug) then
                  call sm_print(6, sub(isub)%i_aaug_sparse, sub(isub)%j_aaug_sparse, sub(isub)%aaug_sparse, &
                                sub(isub)%laaug, sub(isub)%nnzaaug)
               end if
            end if
         end if
         write(*,*) '     matrix PHIS prepared:    ', sub(isub)%is_phisi_prepared
         if (debug) then
            do i = 1,sub(isub)%lphisi1
               write(*,'(1000f13.6)') (sub(isub)%phisi(i,j),j = 1,sub(isub)%lphisi2)
            end do
         end if
         write(*,*) '     coarse matrix prepared:  ', sub(isub)%is_coarse_prepared
!         write(*,'(f13.6)') (sub(isub)%coarsem(j),j = 1,sub(isub)%lcoarsem)
         write(*,*) ' embedding of corse matrix :  '
         write(*,'(i8)') (sub(isub)%indrowc(j),j = 1,sub(isub)%lindrowc)
      end do
      write(*,*) '******************************'

end subroutine

!*************************************
subroutine dd_print_cnode(isub,icnode)
!*************************************
! Subroutine for printing content of one coarse node
      implicit none

      integer,intent(in) :: isub, icnode

! basic structure
      write(*,*) '****** start coarse node export '
      write(*,*) '     coarse node number:      ', icnode
      write(*,*) '     type of coarse node:     ', sub(isub)%cnodes(icnode)%itype
      write(*,*) '     used for constraints?:   ', sub(isub)%cnodes(icnode)%used
      write(*,*) '     coordinates:             ', sub(isub)%cnodes(icnode)%xyz
      write(*,*) '****** where it maps to? '
      write(*,*) '     global coarse node number:', sub(isub)%cnodes(icnode)%global_cnode_number
      write(*,*) '     number of coarse degrees of freedom:', sub(isub)%cnodes(icnode)%ncdof
      write(*,*) '     indices of coarse dof:', sub(isub)%cnodes(icnode)%igcdof
!      write(*,*) '****** where it maps from? '
!      write(*,*) '     number of nodes it contains:', sub(isub)%cnodes(icnode)%nnod
!      write(*,*) '     indices of nodes on subdomain int:', sub(isub)%cnodes(icnode)%insin
      write(*,*) '     number of nodes it contains:', sub(isub)%cnodes(icnode)%nvar
      write(*,*) '     indices of variables on subdomain int:', sub(isub)%cnodes(icnode)%ivsivn
      write(*,*) '****** end coarse nodes export '
end subroutine

!*********************
subroutine dd_finalize
!*********************
! Subroutine for initialization of subdomain data
      implicit none

! local variables
      integer :: isub

! deallocate basic structure
      if (allocated(sub)) then
         ! clear data of particular subdomains
         do isub = 1,lsub
            call dd_clear_subdomain(isub)
         end do

         deallocate (sub)
      end if

end subroutine

end module module_dd
