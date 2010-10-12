module module_levels
!*******************
! Module for handling levels in multilevel BDDC 
! Jakub Sistek, Praha 2/2010

!     definition of MUMPS structure
      use dmumps_struc_def

      implicit none

! type of real variables
      integer,parameter,private :: kr = kind(1.D0)
! numerical zero
      real(kr),parameter,private :: numerical_zero = 1.e-12_kr

! debugging 
      logical,parameter,private :: debug = .true.

! type for data about levels
      type levels_type
         integer ::             nelem    ! number of elements (subdomains) on level
         integer ::             nnod     ! number of nodes (corse nodes) on level
         integer ::             ndof     ! number of dof (corse dof) on level

         integer ::             nsub     ! number of subdomains on level

         integer ::             nnodc    ! number of corners on level
         integer ::             nedge    ! number of edges on level
         integer ::             nface    ! number of faces on level

         ! description of mesh
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

         ! subdomain data
         integer             :: lindexsub    
         integer,allocatable ::  indexsub(:) ! indices of elements in array sub

      end type levels_type

      integer, private ::                          nlevels
      integer, private ::                          llevels
      type(levels_type), allocatable, private ::    levels(:)

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

!*****************************************************************************************
subroutine levels_pc_setup(problemname,myid,nproc,comm_all,comm_self,matrixtype,ndim,nsub,&
                           use_arithmetic, use_adaptive)
!*****************************************************************************************
! subroutine for multilevel BDDC preconditioner setup
      implicit none
      include "mpif.h"

! name of problem
      character(*),intent(in) :: problemname
! number of processor
      integer,intent(in) :: myid
! number of all processors
      integer,intent(in) :: nproc
! communicator for all processors
      integer,intent(in) :: comm_all
! communicator for single processor
      integer,intent(in) :: comm_self
! type of matrix (0 - nosymetric, 1 - SPD, 2 - general symmetric)
      integer,intent(in) :: matrixtype
! dimension
      integer,intent(in) :: ndim
! number of subdomains
      integer,intent(in) :: nsub
! Use arithmetic averages on globs as constraints?
      logical,intent(in) :: use_arithmetic
! Use adaptive constraints on faces?
      logical,intent(in) :: use_adaptive


      ! local vars
      integer :: ilevel

      ilevel = 1
      call levels_read_level_from_file(problemname,myid,comm_all,ndim,ilevel)

      ilevel = 2
      call levels_read_level_from_file(problemname,myid,comm_all,ndim,ilevel)

      ! associate subdomains with first level
      do ilevel = 1,nlevels-1
         call levels_prepare_standard_level(ilevel,nsub,1,nsub)
      end do
      call levels_prepare_last_level(myid,nproc,comm_all,comm_self,matrixtype,ndim,problemname,&
                                     use_arithmetic, use_adaptive)

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

      integer :: i,j

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

! continue only for levels 2 and larger
      if (ilevel.ge.2) then

         lnnet = nelem
         lxyz1 = nnod
         lxyz2 = ndim

         allocate(inet(linet),nnet(lnnet),xyz(lxyz1,lxyz2))

         if (myid.eq.0) then
            read(idlevel,*) inet
            read(idlevel,*) nnet
            do i = 1,nnod
               read(idlevel,*) (xyz(i,j),j = 1,ndim)
            end do
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

!*************************************************************************
subroutine levels_prepare_standard_level(ilevel,nsub,isubstart,isubfinish)
!*************************************************************************
! Subroutine for building the standard level
      use module_utils
      implicit none
      include "mpif.h"

      integer,intent(in) :: ilevel     ! index of level
      integer,intent(in) :: nsub       ! number of subdomains
      integer,intent(in) :: isubstart  ! index of first subdomain of this level
      integer,intent(in) :: isubfinish ! index of last subdomain of this level

      ! local vars
      integer :: lindexsub, i, ind

      lindexsub = nsub
      levels(ilevel)%lindexsub = lindexsub
      allocate(levels(ilevel)%indexsub(lindexsub))

      ! associate numbers of subdomains in array sub with this level
      ind = 0
      do i = isubstart,isubfinish
         ind = ind + 1

         levels(ilevel)%indexsub(ind) = i
      end do

end subroutine

!**********************************************************************************************
subroutine levels_prepare_last_level(myid,nproc,comm_all,comm_self,matrixtype,ndim,problemname,&
                                     use_arithmetic,use_adaptive)
!**********************************************************************************************
! Subroutine for building the coarse problem on root process

      use module_adaptivity
      use module_dd
      use module_sm
      use module_utils
      implicit none
      include "mpif.h"

      integer,intent(in) :: myid
      integer,intent(in) :: nproc
      integer,intent(in) :: comm_all
      integer,intent(in) :: comm_self
      integer,intent(in) :: matrixtype
      integer,intent(in) :: ndim
      character(90),intent(in) :: problemname

      ! Use arithmetic averages on globs as constraints?
      logical,intent(in) :: use_arithmetic
      ! Use adaptive constraints on faces?
      logical,intent(in) :: use_adaptive

      ! local vars
      integer :: ilevel
      integer :: glob_type
      integer :: nnod, nnodc, nedge, nface, isub, nsub, ndof
      integer :: i

      integer :: la, nnz
      integer,allocatable :: i_sparse(:),j_sparse(:)
      real(kr),allocatable :: a_sparse(:)

      integer ::            lnndf
      integer,allocatable :: nndf(:)

      integer :: mumpsinfo
      logical :: remove_original 

      ! adaptivity variables
      integer :: idpair
      character(100) :: filename
      integer :: npair, npair_locx

      ! last level has index of number of levels
      ilevel = nlevels

! number of elements
      nsub = levels(ilevel)%nelem

! prepare nndf
      nnodc = levels(ilevel)%nnodc
      nedge = levels(ilevel)%nedge
      nface = levels(ilevel)%nface

      ! in nndf, nodes are ordered as corners - edges - faces
      nnod  = nnodc + nedge + nface
      lnndf = nnod
      allocate (nndf(lnndf))

      nndf = 0
      ! in corners, prescribe as many constraints as dimension
      nndf(1:nnodc) = ndim
      if (use_arithmetic) then
         ! on edges
         nndf(nnodc+1:nnodc+nedge) = ndim
         ! on faces
         if (use_adaptive) then
            nndf(nnodc+nedge+1:nnodc+nedge+nface) = 0
         else
            nndf(nnodc+nedge+1:nnodc+nedge+nface) = ndim
         end if
      else
         ! on edges
         nndf(nnodc+1:nnodc+nedge) = 0
         ! on faces
         nndf(nnodc+nedge+1:nnodc+nedge+nface) = 0
      end if

      call dd_distribute_subdomains(nsub,nproc)
      call dd_read_mesh_from_file(myid,trim(problemname))
      call dd_read_matrix_from_file(myid,comm_all,trim(problemname),matrixtype)
      call dd_assembly_local_matrix(myid)
      remove_original = .false.
      call dd_matrix_tri2blocktri(myid,remove_original)

      do isub = 1,nsub
         call dd_prepare_schur(myid,comm_self,isub)
      end do

!      call dd_print_sub(myid)
      call dd_create_neighbouring(myid,nsub,comm_all)

      ! weights
      call dd_weights_prepare(myid, nsub, comm_all)

      ! reduced RHS
      call dd_prepare_reduced_rhs(myid,nsub, comm_all)

!      call dd_print_sub(myid)

      ! BDDC data
      ! start preparing cnodes
      do isub = 1,nsub
         call dd_get_cnodes(myid,isub)
      end do
      ! load arithmetic averages on edges
      if (use_arithmetic) then
         glob_type = 2
         do isub = 1,nsub
            call dd_load_arithmetic_constraints(myid,isub,glob_type)
         end do
         ! load arithmetic averages on faces
         if (.not.use_adaptive) then
            glob_type = 1
            do isub = 1,nsub
               call dd_load_arithmetic_constraints(myid,isub,glob_type)
            end do
         end if
      end if
      ! prepare matrix C for corners and arithmetic averages on edges
      do isub = 1,nsub
         call dd_embed_cnodes(myid,isub,nndf,lnndf)
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
         call adaptivity_update_ndof(nndf,lnndf,nnodc,nedge,nface)
   
         call adaptivity_finalize

         ! prepare AGAIN matrix C, now for corners, arithmetic averages on edges and adaptive on faces
         do isub = 1,nsub
            call dd_embed_cnodes(myid,isub,nndf,lnndf)
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

      end if

      ! print the output
!      call dd_print_sub(myid)

      ndof = sum(nndf)
      ! load nndf to level
      levels(ilevel)%ndof = ndof
      levels(ilevel)%lnndf = lnndf
      allocate(levels(ilevel)%nndf(lnndf))
      do i = 1,lnndf
         levels(ilevel)%nndf(i) = nndf(i)
      end do

!      if (debug) then
!         write(*,*) 'ilevel = ',ilevel
!         write(*,*) 'LEVEL',ilevel-1,', indexsub = ',levels(ilevel-1)%indexsub
!      end if

      if (levels(ilevel)%nelem .ne. levels(ilevel-1)%lindexsub) then
         write(*,*) 'LEVELS_PREPARE_LAST_LEVEL: Error in levels consistency.'
         call error_exit
      end if

      ! find length of coarse matrix
      call dd_get_my_coarsem_length(myid,levels(ilevel-1)%indexsub,levels(ilevel-1)%lindexsub,la)

!      write(*,*) 'myid =',myid,'la =',la

! Allocate proper size of matrix A on processor
      allocate(i_sparse(la), j_sparse(la), a_sparse(la))
      i_sparse = 0
      j_sparse = 0
      a_sparse = 0.0D0

      ! load coarse matrix
      call dd_get_my_coarsem(myid,matrixtype,levels(ilevel-1)%indexsub,levels(ilevel-1)%lindexsub, &
                             i_sparse, j_sparse, a_sparse, la)

! Assembly entries in matrix
      call sm_assembly(i_sparse,j_sparse,a_sparse,la,nnz)

      nnz = la
!      write(*,*) 'myid =',myid,'la =',la
!      call sm_print(6, i_sparse, j_sparse, a_sparse, la, nnz)

! Initialize MUMPS
      call mumps_init(mumps_coarse,comm_all,matrixtype)
      write(*,*)'myid =',myid,': MUMPS Initialized'
      call flush(6)

! Level of information from MUMPS
      mumpsinfo  = 1
      call mumps_set_info(mumps_coarse,mumpsinfo)

! Load matrix to MUMPS
      ndof = levels(ilevel)%ndof
      !write(*,*) 'coarse ndof',ndof
      call mumps_load_triplet(mumps_coarse,ndof,nnz,i_sparse,j_sparse,a_sparse,la)
      write(*,*)'myid =',myid,': Triplet loaded'
      call flush(6)

! Analyze matrix
      call mumps_analyze(mumps_coarse)
      write(*,*)'myid =',myid,': Matrix analyzed'
      call flush(6)

! Factorize matrix
      call mumps_factorize(mumps_coarse)
      write(*,*)'myid =',myid,': Matrix factorized'
      call flush(6)

      is_mumps_coarse_ready = .true.

! Clear memory
      deallocate(i_sparse, j_sparse, a_sparse)
      deallocate(nndf)

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

      integer :: lindexsub, ndofs, nnods, nelems, ndofaaugs, ndofcs, ncnodess
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
      lsolc = levels(ilevel + 1)%ndof
      allocate(solcaux(lsolc))
      allocate(solc(lsolc))
      call zero(solcaux,lsolc)

      ! get local contribution to coarse residual
      do is = 1,lindexsub
         isub = levels(ilevel)%indexsub(is)

         call dd_get_coarse_size(myid,isub,ndofcs,ncnodess)
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

         call dd_get_coarse_size(myid,isub,ndofcs,ncnodess)
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
      if (allocated(level%indexsub)) then
         deallocate (level%indexsub)
      end if
      level%lxyz1 = 0
      level%lxyz2 = 0

      level%nelem = 0
      level%nnod  = 0
      level%ndof  = 0

      level%nnodc  = 0
      level%nedge  = 0
      level%nface  = 0

      level%lindexsub = 0

end subroutine

!*************************
subroutine levels_finalize
!*************************
! Subroutine for initialization of levels data
      use module_dd
      use module_mumps
      implicit none

! local variables
      integer :: ilevel
      
! destroy MUMPS structure of the last level
      if (is_mumps_coarse_ready) then
         call mumps_finalize(mumps_coarse)
      end if

! finalize DD module
      call dd_finalize

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

