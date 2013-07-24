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

module module_levels
!*******************
! Module for handling levels in multilevel BDDC 
! Jakub Sistek, Bologna 11/2010

!     definition of MUMPS structure
      use module_dd

      implicit none

! adjustable parameters ############################
! type of real variables
      integer,parameter,private :: kr = kind(1.D0)
! numerical zero
      real(kr),parameter,private :: numerical_zero = 1.e-12_kr
! debugging 
      logical,parameter,private :: debug = .false.
! plot input data in ParaView - useful for debugging
      logical,parameter,private :: plot_inputs = .true.
! profiling 
      logical,private :: profile = .true.
! damping division
      logical,parameter,private :: damp_division = .false.
! maximal allowed length of file names
      integer,parameter,private :: lfnamex = 130
! should parallel search of globs be used? (some corrections on globs may not be available)
      logical,parameter,private :: levels_parallel_globs = .true.
! should parallel search of neighbours be used? (distributed graph rather than serial graph)
      logical,parameter,private :: levels_parallel_neighbouring = .true.
! use ParMETIS for determination of neighbouring? Alternatively, guessing based on index ranges will be used
      logical,parameter,private :: levels_parmetis_neighbouring = .false.
! use globs from *.GLB and *.CN file?
      logical,parameter,private :: levels_load_globs = .false.
! load PAIRS for adaptivity from file?
      logical,parameter,private :: levels_load_pairs = .false.
! should disconnected subdomains be connected? (not applicable for parallel division)
      logical,parameter,private :: levels_correct_division = .false.
! adjustable parameters ############################

! shift of indexing between C and Fortran
      integer,private :: levels_numshift

! type for data about levels
      type levels_type

         logical :: i_am_active_in_this_level = .false.! this is a flag for switching processors at levels
         logical :: is_new_comm_created       = .false.! this is a flag to remember, if a new communicator was created for this level
         integer :: comm_all  ! global communicator for the level
         integer :: comm_self ! local communicator for the level
         
         
         logical :: global_loading = .true. ! global data are normally expected - this is switched off by local loaded subdomains

         integer :: nelem    ! number of elements on level
         integer :: nnod     ! number of nodes on level
         integer :: ndof     ! number of dof on level
         integer :: ndim     ! number of spacial dimensions
         integer :: meshdim  ! topological dimension of mesh: for 3D problems, = ndim, for 3D shells = 2, for 3D beams  = 1

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

         integer ::             idelm    ! unit with opened Fortran unformatted file with element matrices

         integer ::             neighbouring = 0 ! how many nodes must be shared to call two elements neighbours
         logical ::      load_division = .false. ! should division into subdomains be read from file?

         logical :: compute_reactions = .false. ! should reactions be computed on the level?
         

         real(kr) :: adaptivity_estimate ! estimate of condition number on given level

         logical :: is_level_prepared = .false. ! whole level prepared

      end type levels_type

      integer, private ::                                  nlevels
      integer, private ::                                  iactive_level
      integer, private ::                                  llevels
      type(levels_type), allocatable, target, private ::    levels(:)

      type(DMUMPS_STRUC), private :: mumps_coarse  
      logical :: is_mumps_coarse_ready = .false.

contains

!*************************************************************************************
subroutine levels_init(nl,nsublev,lnsublev,nsub_loc_1,comm_init,verbose_level,numbase)
!*************************************************************************************
! Subroutine for initialization of levels data and creating communicators
      ! NOTE:
      ! although the code supports reducing communicators, some routines require
      ! them to be in triangular shape, meaning following IDs:
      ! L(n)   :   0 1 2
      ! L(n-1) :   0 1 2 3 4
      ! ...
      ! L(3)   :   0 1 2 3 4 5 6 7 8
      ! L(2)   :   0 1 2 3 4 5 6 7 8 9 10
      ! L(1)   :   0 1 2 3 4 5 6 7 8 9 10
      ! L(0)   :   0 1 2 3 4 5 6 7 8 9 10
      use module_utils
      use module_pp
      implicit none
      include "mpif.h"

! given number of levels
      integer,intent(in) :: nl
! global number of subdomains in all levels - i.e. sum of all subdomains at first level
      integer,intent(in) :: lnsublev
      integer,intent(in) ::  nsublev(lnsublev)

      integer,intent(inout) ::  nsub_loc_1    ! number of subdomains on first level dedicated to the process - if -1, 
                                              ! the solver decides and returns the value, to which it was allocated
                                              ! note that you can either determine all or none, i.e. all cores should have -1
                                              ! if this value should be computed

      integer, intent(in)::     comm_init     ! initial global communicator (possibly MPI_COMM_WORLD)
      integer, intent(in)::     verbose_level ! verbosity level
      integer, intent(in)::     numbase       ! first index of arrays ( 0 for C, 1 for Fortran )

      ! local vars 
      character(*),parameter:: routine_name = 'LEVELS_INIT'
      integer,parameter :: ilevel = 0
      integer :: isub, isub_loc, nsub, nsub_loc
      integer :: myid, nproc, comm_self, comm_all, ierr
      integer ::            lranks
      integer,allocatable :: ranks(:)
      integer :: mpi_group_new, mpi_group_old, comm_new, comm_old
      integer :: myid_old, nproc_old
      integer :: nsub_1
      integer :: i, ir

! initial checks 
      if (nl.lt.2) then
         call error(routine_name,'Number of levels must be at least 2.')
      end if
      if (nsublev(nl).ne.1) then
         call error(routine_name,'Number of subdomains at last level must be 1.')
      end if
      if ( nsub_loc_1 .ne. -1 ) then ! check that the sum applies
!***************************************************************PARALLEL
         call MPI_ALLREDUCE( nsub_loc_1, nsub_1, 1, MPI_INTEGER, MPI_SUM, comm_init, ierr) 
!***************************************************************PARALLEL
         ! the result should compare to number of levels at first level
         if ( nsub_1 .ne. nsublev(1) ) then
            call error(routine_name,'Number of local subdomains on first level must sum up to given global number:', nsublev(1))
         end if
      end if
      if (any(nsublev(2:nl).ge.nsublev(1:nl-1))) then
         call error(routine_name,'Number of subdomains must decrease monotonically with increasing level.')
      end if

! set verbosity
      if (verbose_level .eq. 2 ) then
         call levels_set_profile_on
      end if

! set index shift for C/Fortran
      levels_numshift = 1 - numbase

! initialize basic structure
      nlevels = nl
      llevels = nlevels
      allocate(levels(0:llevels))

      ! initialize zero level
      iactive_level = 0
      levels(iactive_level)%i_am_active_in_this_level = .true.
      levels(iactive_level)%is_new_comm_created = .false.
      levels(iactive_level)%comm_all  = comm_init
      levels(iactive_level)%comm_self = MPI_COMM_SELF

      ! prepare communicators for levels
      do iactive_level = 1,nlevels-1
         nsub = nsublev(iactive_level)
         levels(iactive_level)%nsub = nsub

         ! if I did not compute on previous level, do not compute on this one as well
         if (.not.levels(iactive_level-1)%i_am_active_in_this_level) then
            levels(iactive_level)%i_am_active_in_this_level = .false.
            levels(iactive_level)%is_new_comm_created = .false.
            return
         end if

         comm_self   = MPI_COMM_SELF
         if (iactive_level.eq.1) then
            comm_old = comm_init
         else
            comm_old = levels(iactive_level-1)%comm_all
         end if
         ! get absolute number of processes
         ! find ID in old communicator
         call MPI_COMM_SIZE(comm_old,nproc_old,ierr)
         call MPI_COMM_RANK(comm_old,myid_old,ierr)

         if (nsub.ge.nproc_old .or. iactive_level.eq.1) then
            comm_all = comm_old
            levels(iactive_level)%is_new_comm_created       = .false.
            levels(iactive_level)%i_am_active_in_this_level = .true.
         else
            if (debug) then
               if (myid_old.eq.0) then
                  call info(routine_name,'reducing communicator for level',iactive_level)
               end if
            end if
            lranks = nsub
            allocate(ranks(lranks))
            ! TODO: build new communicators with better data locality (derived from division on previous level)
            do isub = 1,nsub
               ranks(isub) = isub - 1 
            end do
            ! extract group of processors of previous level
            call MPI_COMM_GROUP(comm_old,mpi_group_old,ierr)
            ! prepare new group of processes
            call MPI_GROUP_INCL(mpi_group_old,nsub,ranks,mpi_group_new,ierr)
            ! create new communicator
            call MPI_COMM_CREATE(comm_all,mpi_group_new,comm_new,ierr)
            ! free groups
            call MPI_GROUP_FREE(mpi_group_new,ierr)
            call MPI_GROUP_FREE(mpi_group_old,ierr)

            comm_all = comm_new
            levels(iactive_level)%is_new_comm_created  = .true.

            ! I won't work at this level 
            if (.not.any(ranks.eq.myid_old)) then
               ! I have not been included to the new communicator
               levels(iactive_level)%i_am_active_in_this_level = .false.
            else
               ! I have been included to the new communicator
               levels(iactive_level)%i_am_active_in_this_level = .true.
            end if
            deallocate(ranks)
         end if

         levels(iactive_level)%comm_self = comm_self
         levels(iactive_level)%comm_all  = comm_all

         if (levels(iactive_level)%i_am_active_in_this_level) then

            ! orient in the communicator
            call MPI_COMM_SIZE(comm_all,nproc,ierr)
            call MPI_COMM_RANK(comm_all,myid,ierr)

            ! prepare distribution of subdomains for the level
            levels(iactive_level)%lsub2proc = nproc + 1
            allocate(levels(iactive_level)%sub2proc(levels(iactive_level)%lsub2proc))

            if ( iactive_level .eq. 1 .and. nsub_loc_1 .ne. -1 ) then ! user defined subdomain distribution on first level
!***************************************************************PARALLEL
               call MPI_ALLGATHER( nsub_loc_1, 1, MPI_INTEGER, levels(iactive_level)%sub2proc, 1, MPI_INTEGER, &
                                   comm_all,ierr)
!***************************************************************PARALLEL
               ! the array now contains counts, change it to starts
               do i = 2,nproc
                  levels(iactive_level)%sub2proc(i) = levels(iactive_level)%sub2proc(i-1) + levels(iactive_level)%sub2proc(i)
               end do
               ! shift it one back and add one 
               do ir = 0, nproc - 1 ! reverse index
                  i = nproc + 1 - ir

                  levels(iactive_level)%sub2proc(i) = levels(iactive_level)%sub2proc(i-1) + 1
               end do
               ! put one in the beginning
               levels(iactive_level)%sub2proc(1) = 1

               nsub_loc = nsub_loc_1

            else
               call pp_distribute_linearly(nsub,nproc,levels(iactive_level)%sub2proc,levels(iactive_level)%lsub2proc)
               nsub_loc    = levels(iactive_level)%sub2proc(myid+2) - levels(iactive_level)%sub2proc(myid+1)
               ! return to user the number of assigned subdomains
               if ( iactive_level .eq. 1 ) then
                  nsub_loc_1 = nsub_loc
               end if
            end if

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

         else
            levels(iactive_level)%lsub2proc = 0
            allocate(levels(iactive_level)%sub2proc(levels(iactive_level)%lsub2proc))
            levels(iactive_level)%nsub_loc  = 0
            levels(iactive_level)%lindexsub = 0
            allocate(levels(iactive_level)%indexsub(levels(iactive_level)%lindexsub))
            levels(iactive_level)%lsubdomains = 0
            allocate(levels(iactive_level)%subdomains(levels(iactive_level)%lsubdomains))
         end if

      end do

      ! prepare communicators for last level
      iactive_level = nlevels
      levels(iactive_level)%comm_self = comm_self
      levels(iactive_level)%comm_all  = levels(iactive_level-1)%comm_all
      levels(iactive_level)%i_am_active_in_this_level = levels(iactive_level-1)%i_am_active_in_this_level
      levels(iactive_level)%is_new_comm_created       = .false.

end subroutine

!***********************************************************************************
subroutine levels_upload_global_data(nelem,nnod,ndof,ndim,meshdim,&
                                     inet,linet,nnet,lnnet,nndf,lnndf,xyz,lxyz1,lxyz2,&
                                     ifix,lifix, fixv,lfixv, rhs,lrhs, sol,lsol, idelm, &
                                     neighbouring, load_division)
!***********************************************************************************
! Subroutine for loading global data as zero level
      use module_utils
      implicit none

      integer, intent(in):: nelem, nnod, ndof, ndim, meshdim
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
      integer,intent(in) :: idelm  
! how many nodes are shared to called elements adjacent
      integer,intent(in) :: neighbouring
! use prepared division into subdomains?
      logical,intent(in) :: load_division

      ! local vars 
      character(*),parameter:: routine_name = 'LEVELS_UPLOAD_GLOBAL_DATA'
      integer :: i,j

      ! store unit of file with element matrices at first level
      levels(1)%idelm = idelm

      ! set active level to zero
      iactive_level = 0

      ! check number of elements
      if (levels(1)%nsub.gt.nelem) then
         call error(routine_name,'Number of subdomains at first level must not be larger than number of elements.')
      end if

! initialize zero level
      levels(iactive_level)%nsub    = nelem
      levels(iactive_level)%nnodc   = nnod
      levels(iactive_level)%ndofc   = ndof
      levels(iactive_level)%ndim    = ndim
      levels(iactive_level)%meshdim = meshdim

      levels(iactive_level)%linetc = linet  
      if ( .not. allocated(levels(iactive_level)%inetc)) then
         allocate(levels(iactive_level)%inetc(levels(iactive_level)%linetc))
      end if
      do i = 1,levels(iactive_level)%linetc
         levels(iactive_level)%inetc(i) = inet(i) + levels_numshift
      end do
      levels(iactive_level)%lnnetc = lnnet  
      if ( .not. allocated(levels(iactive_level)%nnetc)) then
         allocate(levels(iactive_level)%nnetc(levels(iactive_level)%lnnetc))
      end if
      do i = 1,levels(iactive_level)%lnnetc
         levels(iactive_level)%nnetc(i) = nnet(i)
      end do
      levels(iactive_level)%lnndfc = lnndf  
      if ( .not. allocated(levels(iactive_level)%nndfc)) then
         allocate(levels(iactive_level)%nndfc(levels(iactive_level)%lnndfc))
      end if
      do i = 1,levels(iactive_level)%lnndfc
         levels(iactive_level)%nndfc(i) = nndf(i)
      end do
      levels(iactive_level)%lxyzc1 = lxyz1  
      levels(iactive_level)%lxyzc2 = lxyz2  
      if ( .not. allocated(levels(iactive_level)%xyzc)) then
         allocate(levels(iactive_level)%xyzc(levels(iactive_level)%lxyzc1,levels(iactive_level)%lxyzc2))
      end if
      do j = 1,levels(iactive_level)%lxyzc2
         do i = 1,levels(iactive_level)%lxyzc1
            levels(iactive_level)%xyzc(i,j) = xyz(i,j)
         end do
      end do
      levels(iactive_level)%lifixc = lifix
      if ( .not. allocated(levels(iactive_level)%ifixc)) then
         allocate(levels(iactive_level)%ifixc(levels(iactive_level)%lifixc))
      end if
      do i = 1,levels(iactive_level)%lifixc
         levels(iactive_level)%ifixc(i) = ifix(i)
      end do
      levels(iactive_level)%lfixvc = lfixv
      if ( .not. allocated(levels(iactive_level)%fixvc)) then
         allocate(levels(iactive_level)%fixvc(levels(iactive_level)%lfixvc))
      end if
      do i = 1,levels(iactive_level)%lfixvc
         levels(iactive_level)%fixvc(i) = fixv(i)
      end do
      levels(iactive_level)%lrhsc = lrhs
      if ( .not. allocated(levels(iactive_level)%rhsc)) then
         allocate(levels(iactive_level)%rhsc(levels(iactive_level)%lrhsc))
      end if
      do i = 1,levels(iactive_level)%lrhsc
         levels(iactive_level)%rhsc(i) = rhs(i)
      end do
      levels(iactive_level)%lsolc = lsol
      if ( .not. allocated(levels(iactive_level)%solc)) then
         allocate(levels(iactive_level)%solc(levels(iactive_level)%lsolc))
      end if
      do i = 1,levels(iactive_level)%lsolc
         levels(iactive_level)%solc(i) = sol(i)
      end do
      levels(iactive_level)%use_initial_solution = .true.


      levels(iactive_level)%is_level_prepared = .true.

      iactive_level = 1
      levels(iactive_level)%neighbouring  = neighbouring
      levels(iactive_level)%load_division = load_division

      ! set computing of reactions for the first level
      levels(iactive_level)%compute_reactions = .true.


end subroutine

!*******************************************************************************
subroutine levels_change_global_data(ifix,lifix, fixv,lfixv, rhs,lrhs, sol,lsol)
!*******************************************************************************
! Subroutine for changing global data as zero level
      use module_utils
      implicit none

      integer, intent(in):: lifix
      integer, intent(in)::  ifix(lifix)
      integer, intent(in):: lfixv
      real(kr), intent(in):: fixv(lfixv)
      integer, intent(in):: lrhs
      real(kr), intent(in):: rhs(lrhs)
      integer, intent(in):: lsol
      real(kr), intent(in):: sol(lsol)

      ! local vars 
      character(*),parameter:: routine_name = 'LEVELS_CHANGE_GLOBAL_DATA'
      integer :: i

      ! set active level to zero
      iactive_level = 0

      ! verify that IFIX did not change
      if (levels(iactive_level)%lifixc .ne. lifix) then
         call error( routine_name, 'Dimension mismatch in IFIX')
      end if
      if ( .not. allocated(levels(iactive_level)%ifixc)) then
         call error( routine_name, 'IFIX not allocated')
      end if
      if ( any(levels(iactive_level)%ifixc.ne.ifix)) then
         call error( routine_name, 'IFIX cannot be changed when loading new data')
      end if
      if ( .not. allocated(levels(iactive_level)%fixvc)) then
         allocate(levels(iactive_level)%fixvc(levels(iactive_level)%lfixvc))
      end if
      do i = 1,levels(iactive_level)%lfixvc
         levels(iactive_level)%fixvc(i) = fixv(i)
      end do
      levels(iactive_level)%lrhsc = lrhs
      if ( .not. allocated(levels(iactive_level)%rhsc)) then
         allocate(levels(iactive_level)%rhsc(levels(iactive_level)%lrhsc))
      end if
      do i = 1,levels(iactive_level)%lrhsc
         levels(iactive_level)%rhsc(i) = rhs(i)
      end do
      levels(iactive_level)%lsolc = lsol
      if ( .not. allocated(levels(iactive_level)%solc)) then
         allocate(levels(iactive_level)%solc(levels(iactive_level)%lsolc))
      end if
      do i = 1,levels(iactive_level)%lsolc
         levels(iactive_level)%solc(i) = sol(i)
      end do
      levels(iactive_level)%use_initial_solution = .true.

      ! new global data loaded, now localize them to subdomains
      call levels_localize_new_data

end subroutine

!***********************************************************************************
subroutine levels_upload_subdomain_data(nelem, nnod, ndof, ndim, meshdim, &
                                        isub, nelems, nnods, ndofs, &
                                        inet,linet, nnet,lnnet, nndf,lnndf, &
                                        isngn,lisngn, isvgvn,lisvgvn, isegn,lisegn, &
                                        xyz,lxyz1,lxyz2, &
                                        ifix,lifix, fixv,lfixv, &
                                        rhs,lrhs, is_rhs_complete, &
                                        sol,lsol, &
                                        matrixtype, i_sparse, j_sparse, a_sparse, la, is_assembled, &
                                        user_constraints,luser_constraints1,luser_constraints2, &
                                        element_data,lelement_data1,lelement_data2, &
                                        dof_data,ldof_data )
!***********************************************************************************
! Subroutine for loading LOCAL data of one subdomain at first level
      use module_utils
      implicit none

      integer, intent(in):: nelem, nnod, ndof, ndim, meshdim
      integer, intent(in):: isub, nelems, nnods, ndofs
      integer, intent(in):: linet
      integer, intent(in)::  inet(linet)
      integer, intent(in):: lnnet
      integer, intent(in)::  nnet(lnnet)
      integer, intent(in):: lnndf
      integer, intent(in)::  nndf(lnndf)
      integer, intent(in):: lisngn
      integer, intent(in)::  isngn(lisngn)
      integer, intent(in):: lisvgvn
      integer, intent(in)::  isvgvn(lisvgvn)
      integer, intent(in):: lisegn
      integer, intent(in)::  isegn(lisegn)
      integer, intent(in):: lxyz1, lxyz2
      real(kr), intent(in):: xyz(lxyz1,lxyz2)
      integer, intent(in):: lifix
      integer, intent(in)::  ifix(lifix)
      integer, intent(in):: lfixv
      real(kr), intent(in):: fixv(lfixv)
      integer, intent(in):: lrhs
      real(kr), intent(in):: rhs(lrhs)
      logical, intent(in):: is_rhs_complete
      integer, intent(in):: lsol
      real(kr), intent(in):: sol(lsol)
      ! LOCAL matrix triplet i, j, a(i,j) 
      integer, intent(in)::  matrixtype    ! type of matrix (MUMPS-like)
      integer, intent(in)::  la            ! length of previous arrays (= number of nonzeros for assembled matrix)
      integer, intent(in)::  i_sparse(la)  ! array of row indices
      integer, intent(in)::  j_sparse(la)  ! array of column indices
      real(kr), intent(in):: a_sparse(la)  ! array of values
      logical, intent(in)::  is_assembled  ! is the array assembled? 
                                           !  FALSE = no, it can contain repeated entries
                                           !  TRUE  = yes, it is sorted and doesn't contain repeated index pairs
      integer, intent(in)::  luser_constraints1 ! number of rows in matrix of constraints
      integer, intent(in)::  luser_constraints2 ! number of columns in matrix of constraints
      real(kr), intent(in):: user_constraints(luser_constraints1*luser_constraints2) ! array for additional constraints
      integer, intent(in)::  lelement_data1 ! number of rows in matrix of element data
      integer, intent(in)::  lelement_data2 ! number of columns in matrix of element data
      real(kr), intent(in):: element_data(lelement_data1*lelement_data2) ! array for additional data about elements
      integer, intent(in)::  ldof_data           ! length of array DOF_DATA
      real(kr), intent(in):: dof_data(ldof_data) ! array for additional data on dofs      

      ! local vars 
      character(*),parameter:: routine_name = 'LEVELS_UPLOAD_SUBDOMAIN_DATA'

      logical :: is_mesh_loaded 
      integer :: isub_loc

      ! set active level to zero
      iactive_level = 0

      if ( .not. levels(iactive_level)%is_level_prepared ) then
         ! initialize zero level
         levels(iactive_level)%nsub    = nelem
         levels(iactive_level)%nnodc   = nnod
         levels(iactive_level)%ndofc   = ndof
         levels(iactive_level)%ndim    = ndim
         levels(iactive_level)%meshdim = meshdim

         levels(iactive_level)%lnnetc = nelem  
         levels(iactive_level)%lnndfc = nnod 
         levels(iactive_level)%lxyzc1 = nnod  
         levels(iactive_level)%lxyzc2 = ndim  
         levels(iactive_level)%lifixc = ndof
         levels(iactive_level)%lfixvc = ndof
         levels(iactive_level)%lrhsc  = ndof

         levels(iactive_level)%lsolc  = ndof

         ! mark 0th level as ready
         levels(iactive_level)%is_level_prepared = .true.
      end if

      ! set active level to one
      iactive_level = 1

      ! note that no global data are loaded
      levels(iactive_level)%global_loading = .false.

      ! local index of subdomain
      call get_index(isub+levels_numshift, levels(iactive_level)%indexsub,levels(iactive_level)%lindexsub, isub_loc)
      if ( isub_loc .eq. -1 ) then
         call error( routine_name, 'Index of subdomain not found among local subdomains: ', isub )
      end if

      if (.not. allocated(levels(iactive_level)%subdomains)) then
         call error( routine_name,'memory not prepared for subdomain', isub )
      end if

      ! check that the subdomain is empty
      call dd_is_mesh_loaded( levels(iactive_level)%subdomains(isub_loc), is_mesh_loaded )
      if ( is_mesh_loaded ) then
         call error( routine_name, 'It appears that mesh was already loaded for subdomain', isub )
      end if

      call dd_upload_sub_mesh(levels(iactive_level)%subdomains(isub_loc), nelems, nnods, ndofs, ndim, &
                              nndf,lnndf, nnet,lnnet, levels_numshift, inet,linet, &
                              isngn,lisngn, isvgvn,lisvgvn, isegn,lisegn,&
                              xyz,lxyz1,lxyz2)
      call dd_upload_bc(levels(iactive_level)%subdomains(isub_loc), ifix,lifix, fixv,lfixv)
      call dd_upload_rhs(levels(iactive_level)%subdomains(isub_loc), rhs,lrhs, is_rhs_complete)

      ! upload initial solution to DD structure
      call dd_upload_solution(levels(iactive_level)%subdomains(isub_loc), sol,lsol)
      levels(iactive_level)%use_initial_solution = .true.

      ! check matrix format
      ! type of matrix (0 - nosymetric, 1 - SPD, 2 - general symmetric)
      if ( ( matrixtype.eq.1 .or. matrixtype.eq.2 ) .and. any( i_sparse .gt. j_sparse ) ) then
         call error( routine_name, &
                     'Matrix in wrong format - only upper triangle is expected for symmetric problems, matrixtype = ', &
                     matrixtype )
      end if
      ! load matrix to our structure
      call dd_load_matrix_triplet(levels(iactive_level)%subdomains(isub_loc), matrixtype, levels_numshift, &
                                  i_sparse,j_sparse,a_sparse,la,la,is_assembled)

      ! assembly it if needed
      if (.not. is_assembled) then
         call dd_assembly_local_matrix(levels(iactive_level)%subdomains(isub_loc))
      end if

      ! check matrix
      call dd_check_local_matrix(levels(iactive_level)%subdomains(isub_loc))

      ! load user's constraints
      call dd_upload_sub_user_constraints(levels(iactive_level)%subdomains(isub_loc), &
                                          user_constraints,luser_constraints1,luser_constraints2)

      ! load element data 
      call dd_upload_sub_element_data(levels(iactive_level)%subdomains(isub_loc), &
                                      element_data,lelement_data1,lelement_data2)

      ! load dof data 
      call dd_upload_sub_dof_data(levels(iactive_level)%subdomains(isub_loc), dof_data,ldof_data)


      ! initialize first level if all subdomains were loaded
      levels(iactive_level)%nelem   = nelem
      levels(iactive_level)%nnodc   = nnod
      levels(iactive_level)%ndofc   = ndof
      levels(iactive_level)%ndim    = ndim
      levels(iactive_level)%meshdim = meshdim


end subroutine

!***********************************************************************************
subroutine levels_change_subdomain_data(isub, &
                                        ifix,lifix, fixv,lfixv, &
                                        rhs,lrhs, is_rhs_complete, &
                                        sol,lsol)
!***********************************************************************************
! Subroutine for loading new LOCAL data - right-hand side, new values of Dirichlet
! boundary conditions LOCAL data of one subdomain at first level
      use module_utils
      implicit none

      integer, intent(in):: isub
      integer, intent(in):: lifix
      integer, intent(in)::  ifix(lifix)
      integer, intent(in):: lfixv
      real(kr), intent(in):: fixv(lfixv)
      integer, intent(in):: lrhs
      real(kr), intent(in):: rhs(lrhs)
      logical, intent(in):: is_rhs_complete
      integer, intent(in):: lsol
      real(kr), intent(in):: sol(lsol)

      ! local vars 
      character(*),parameter:: routine_name = 'LEVELS_CHANGE_SUBDOMAIN_DATA'
      logical :: is_mesh_loaded
      integer :: isub_loc

      ! set active level to one
      iactive_level = 1

      ! check prerequisites
      if ( .not. levels(iactive_level)%is_level_prepared ) then
         call error(routine_name, 'Cannot change data at this moment, level is not ready.')
      end if
      ! note that no global data are loaded
      if (levels(iactive_level)%global_loading) then
         call error(routine_name, 'Global loading should be used.')
      end if

      ! local index of subdomain
      call get_index(isub+levels_numshift, levels(iactive_level)%indexsub,levels(iactive_level)%lindexsub, isub_loc)
      if ( isub_loc .eq. -1 ) then
         call error( routine_name, 'Index of subdomain not found among local subdomains: ', isub )
      end if

      if (.not. allocated(levels(iactive_level)%subdomains)) then
         call error( routine_name,'memory not prepared for subdomain', isub )
      end if

      ! check that the subdomain is loaded
      call dd_is_mesh_loaded( levels(iactive_level)%subdomains(isub_loc), is_mesh_loaded )
      if ( .not. is_mesh_loaded ) then
         call error( routine_name, 'It appears that mesh is not loaded for subdomain', isub )
      end if

      ! set new boundary conditions
      call dd_upload_bc(levels(iactive_level)%subdomains(isub_loc), ifix,lifix, fixv,lfixv)

      ! set new right hand side 
      call dd_upload_rhs(levels(iactive_level)%subdomains(isub_loc), rhs,lrhs, is_rhs_complete)

      ! set new initial solution 
      call dd_upload_solution(levels(iactive_level)%subdomains(isub_loc), sol,lsol)
      levels(iactive_level)%use_initial_solution = .true.
      
end subroutine

!*******************************
subroutine levels_setup_new_data
!*******************************
! Subroutine for building new reduced right-hand side and boundary conditions
      use module_utils
      implicit none


      ! local vars 
      character(*),parameter:: routine_name = 'LEVELS_SETUP_NEW_DATA'
      integer :: comm_all, myid, ierr

      ! time variables
      real(kr) :: t_reduced_rhs_prepare

      ! set active level to one
      iactive_level = 1

      ! orient in the communicator
      comm_all  = levels(iactive_level)%comm_all
      call MPI_COMM_RANK(comm_all,myid,ierr)

      ! check prerequisites
      if ( .not. levels(iactive_level)%is_level_prepared ) then
         call error(routine_name, 'Cannot change data at this moment, level is not ready.')
      end if

      call dd_fix_constraints(levels(iactive_level)%subdomains,levels(iactive_level)%lsubdomains, comm_all,&
                              levels(iactive_level)%sub2proc,levels(iactive_level)%lsub2proc,&
                              levels(iactive_level)%indexsub,levels(iactive_level)%lindexsub)

      ! prepare reduced RHS
!-----profile
      if (profile) then
         call MPI_BARRIER(comm_all,ierr)
         call time_start
      end if
!-----profile
      call dd_prepare_reduced_rhs_all(levels(iactive_level)%subdomains,levels(iactive_level)%lsubdomains, &
                                      levels(iactive_level)%sub2proc,levels(iactive_level)%lsub2proc,&
                                      levels(iactive_level)%indexsub,levels(iactive_level)%lindexsub,& 
                                      comm_all)
!-----profile
      if (profile) then
         call MPI_BARRIER(comm_all,ierr)
         call time_end(t_reduced_rhs_prepare)
         if (myid.eq.0) then
            call time_print('preparing new reduced right-hand side',t_reduced_rhs_prepare)
         end if
      end if
!-----profile
      
end subroutine
!*************************************************************************************
subroutine levels_pc_setup( parallel_division,&
                            matrixtype, &
                            use_arithmetic_constraints,&
                            use_adaptive_constraints,&
                            use_user_constraints, &
                            weights_type )
!*************************************************************************************
! subroutine for multilevel BDDC preconditioner setup
      use module_pp
      use module_utils
      implicit none
      include "mpif.h"

! should parallel division be used (ParMETIS instead of METIS)?
      logical,intent(in) :: parallel_division
! type of matrix (0 - nosymetric, 1 - SPD, 2 - general symmetric)
      integer,intent(in) :: matrixtype
! Use arithmetic averages on globs as constraints?
      logical,intent(in) :: use_arithmetic_constraints
! Use adaptive constraints on faces?
      logical,intent(in) :: use_adaptive_constraints
! Use user defined constraints on faces?
      logical,intent(in) :: use_user_constraints
! what type of weights should be used on the interface ?
      ! type of weights:
      ! 0 - weights by cardinality
      ! 1 - weights by diagonal stiffness
      ! 2 - weights based on first row of element data
      ! 3 - weights based on dof data
      ! 4 - weights by Marta Certikova - unit load
      ! 5 - weights by Marta Certikova - unit jump
      ! 6 - weights by Schur row sums for whole subdomain
      ! 7 - weights by Schur row sums computed face by face
      integer,intent(in) :: weights_type


      ! local vars 
      character(*),parameter:: routine_name = 'LEVELS_PC_SETUP'
      integer :: comm_all, ierr, myid, nproc
      real(kr) :: cond_est

      real(kr) :: t_pc_setup

      ! no debug
      !call levels_read_level_from_file(problemname,comm_all,iactive_level)

!-----profile
      if (profile) then
         comm_all = levels(0)%comm_all
         call MPI_BARRIER(comm_all,ierr)
         call time_start
      end if
!-----profile

      ! welcome screen
      iactive_level = 1
      comm_all = levels(iactive_level-1)%comm_all
      call MPI_COMM_RANK(comm_all,myid,ierr)
      call MPI_COMM_SIZE(comm_all,nproc,ierr)
      if (myid.eq.0) then
         call info( routine_name, '==========================================' )
         call info( routine_name, 'Entering BDDCML library                   ' )
         call info( routine_name, '==========================================' )
         call info( routine_name, 'problem properties:                       ' )
         call info( routine_name, 'number of unknowns:  ',levels(iactive_level-1)%ndofc )
         call info( routine_name, 'number of nodes:     ',levels(iactive_level-1)%nnodc )
         call info( routine_name, 'number of elements:  ',levels(iactive_level-1)%nsub )
         call info( routine_name, 'number of subdomains:',levels(iactive_level)%nsub )
         call info( routine_name, 'number of processors:',nproc )
         if (levels(iactive_level)%global_loading) then
            call info( routine_name, 'using mode with GLOBAL data loading' )
         else
            call info( routine_name, 'using mode with PER-SUBDOMAIN data loading' )
         end if
         call info( routine_name, '==========================================' )
      end if

      ! check input
      if (.not. (levels_load_globs .eqv. levels_load_pairs).and. use_adaptive_constraints) then
         call warning(routine_name,'loading globs while not loading pairs or vice versa may result in wrong constraints')
      end if

      ! prepare standard levels 
      do iactive_level = 1,nlevels-1
         if (iactive_level.eq.1 .or. (iactive_level.gt.1.and.levels(iactive_level - 1)%i_am_active_in_this_level)) then
            comm_all = levels(iactive_level-1)%comm_all
            ! orient in the communicator
            call MPI_COMM_RANK(comm_all,myid,ierr)
            if (debug .and. myid.eq.0) then
               call info(routine_name,'Preparing level',iactive_level)
            end if
            call levels_prepare_standard_level(parallel_division,&
                                               matrixtype,iactive_level,&
                                               use_arithmetic_constraints,&
                                               use_adaptive_constraints,&
                                               use_user_constraints,&
                                               weights_type)
         end if
      end do

      ! prediction of condition number
      if (use_adaptive_constraints) then
         comm_all = levels(1)%comm_all
         ! orient in the communicator
         call MPI_COMM_RANK(comm_all,myid,ierr)
         if (myid.eq.0) then
            cond_est = 1._kr
            do iactive_level = 1,nlevels-1
               cond_est = cond_est * levels(iactive_level)%adaptivity_estimate
            end do
            call info(routine_name,'Expected estimated multilevel condition number: ',cond_est)
         end if
      end if

      ! prepare last level
      iactive_level = nlevels
      if (levels(iactive_level - 1)%i_am_active_in_this_level) then
         comm_all = levels(iactive_level-1)%comm_all
         ! orient in the communicator
         call MPI_COMM_RANK(comm_all,myid,ierr)
         if (debug .and. myid.eq.0) then
            call info(routine_name,'Preparing level',iactive_level)
         end if
         call levels_prepare_last_level(matrixtype)
      end if

!-----profile
      if (profile) then
         comm_all = levels(0)%comm_all
         call MPI_BARRIER(comm_all,ierr)
         call time_end(t_pc_setup)
         if (myid.eq.0) then
            call time_print('set-up of the preconditioner',t_pc_setup)
         end if
      end if
!-----profile
end subroutine

!**************************************************************
subroutine levels_read_level_from_file(problemname,comm,ilevel)
!**************************************************************
! Subroutine for loading first two levels from file

      use module_utils
      implicit none
      include "mpif.h"

! name of problem
      character(*),intent(in) :: problemname
! communicator
      integer,intent(in) :: comm
! index of level to import
      integer,intent(in) :: ilevel

! local variables
      character(*),parameter:: routine_name = 'LEVELS_READ_LEVEL_FROM_FILE'
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

      character(lfnamex) :: filename

      character(1) :: levelstring

      ! find my id in the communicator
      call MPI_COMM_RANK(comm,myid,ierr)

      if (myid.eq.0) then
         if (ilevel.lt.10) then
            write(levelstring,'(i1)') ilevel
         else
            call error(routine_name,'Index of level too large for reading from file:',ilevel)
         end if
         filename = trim(problemname)//'.L'//levelstring
         if (debug) then
            call info(routine_name,' Reading data from file '//trim(filename))
         end if
         call allocate_unit(idlevel)
         open (unit=idlevel,file=filename,status='old',form='formatted')
      end if

! read level header
      if (myid.eq.0) then
         read(idlevel,*) indlevel, nelem, nnod, linet
         read(idlevel,*) ncorner, nedge, nface, linetc
         if (indlevel.ne.ilevel) then
            call error(routine_name,'Level number mismatch...')
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
      lxyzc2 = levels(ilevel)%ndim

      allocate(inetc(linetc),nnetc(lnnetc),nndfc(lnndfc),xyzc(lxyzc1,lxyzc2))

      if (myid.eq.0) then
         read(idlevel,*) inetc
         read(idlevel,*) nnetc
         do i = 1,nnodc
            read(idlevel,*) (xyzc(i,j),j = 1,levels(ilevel)%ndim)
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

!*************************************************
subroutine levels_damp_division(ilevel,iets,liets)
!*************************************************
! Subroutine for damping division into subdomains at level to file

      use module_utils
      implicit none

! index of level to import
      integer,intent(in) :: ilevel
! division into subdomains
      integer,intent(in) :: liets
      integer,intent(in) ::  iets(liets)

! local variables
      character(*),parameter:: routine_name = 'LEVELS_DAMP_DIVISION'
      character(lfnamex) :: filename
      character(1) :: levelstring
      integer :: idlevel

      if (ilevel.lt.10) then
         write(levelstring,'(i1)') ilevel
      else
         call error(routine_name,'Index of level too large for reading from file:',ilevel)
      end if
      filename = 'partition_l'//levelstring//'.ES'
      if (debug) then
         call info(routine_name,' Damping division to file '//trim(filename))
      end if
      call allocate_unit(idlevel)
      open (unit=idlevel,file=filename,status='replace',form='formatted')

! write division into file 
      write(idlevel,*) iets
      close(idlevel)
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
      character(*),parameter:: routine_name = 'LEVELS_UPLOAD_LEVEL_MESH'
      integer :: i, j


      ! check that array is allocated
      if (.not. allocated(levels) ) then
         call error(routine_name,'LEVELS not allocated.')
      end if
      ! check that level is within the allowed range
      if (ilevel .lt. lbound(levels,1) .or. ilevel .gt. ubound(levels,1)) then
         call error(routine_name,'ILEVEL out of range of LEVELS array: ',ilevel)
      end if
      ! check sum
      if (nnodc .ne. ncorner + nedge + nface) then
         call error(routine_name,'Coarse nodes number mismatch',nnodc)
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

!*******************************************************************************
subroutine levels_prepare_standard_level(parallel_division,&
                                         matrixtype,ilevel,&
                                         use_arithmetic_constraints,&
                                         use_adaptive_constraints,&
                                         use_user_constraints,&
                                         weights_type)
!*******************************************************************************
! Subroutine for building the standard level
      use module_pp
      use module_adaptivity
      use module_sm
      use module_paraview
      use module_utils
      implicit none
      include "mpif.h"

      logical,intent(in) :: parallel_division ! should parallel division be used?
      integer,intent(in) :: matrixtype
      integer,intent(in) :: ilevel     ! index of level
      ! Use arithmetic averages on globs as constraints?
      logical,intent(in) :: use_arithmetic_constraints
      ! Use adaptive constraints on faces?
      logical,intent(in) :: use_adaptive_constraints
      ! Use user defined constraints on faces?
      logical,intent(in) :: use_user_constraints
      ! type of weights
      integer,intent(in) :: weights_type

      ! local vars
      character(*),parameter:: routine_name = 'LEVELS_PREPARE_STANDARD_LEVEL'

      integer :: myid
      integer :: nproc
      integer :: comm_all, comm_self, ierr
      integer :: graphtype 
      integer :: ides

      integer :: ncorner, nedge, nface, isub, nnodc, ndofc, nelem, nnod
      integer :: edgecut
      integer :: nsub_loc, isub_loc, nsub, glbtype, sub_start

      !  division variables
      integer :: stat(MPI_STATUS_SIZE)
      integer:: nelem_loc, nelem_locx
      ! pointers until mesh is local 
      integer,pointer :: linet
      integer,pointer ::  inet(:)
      integer,pointer :: lnnet
      integer,pointer ::  nnet(:)
      integer,pointer :: liets
      integer,pointer ::  iets(:)
      ! linear partitioning of elements
      integer ::           lelm2proc
      integer,allocatable:: elm2proc(:)
      ! local mesh
      integer ::           lnnet_loc
      integer,allocatable:: nnet_loc(:)
      integer ::             linet_loc
      integer*4,allocatable:: inet_loc(:)
      ! division
      integer ::           lelm2sub
      integer,allocatable:: elm2sub(:)
      integer ::           lpart_loc
      integer,allocatable:: part_loc(:)
      integer ::           lpart_loc_proc
      integer,allocatable:: part_loc_proc(:)
      integer:: el_start
      integer:: el_start_send, el_finish_send, ie, indinet, indproc
      integer:: length_send, length
      integer:: ie_loc, indel, indsub
      integer:: el_start_proc, ie_loc_proc, nelem_loc_proc
      integer:: i,j

      ! data for neighbouring
      integer :: ndofs, nelems, nnods, neighball
      integer :: linets, lnnets
      integer :: pkadjsub
      integer ::           lkadjsub
      integer,allocatable:: kadjsub(:)
      integer :: kinet_loc, knnet_loc
      logical :: use_global_indices
      integer ::           lnelemsa
      integer,allocatable:: nelemsa(:)
      integer,allocatable:: nelemsaaux(:)
      integer ::           liets_linear
      integer,allocatable:: iets_linear(:)
      integer ::            indiets

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
      real(kr),allocatable :: xyzcaux(:,:)

      integer ::            lcnodenumbers
      integer,allocatable :: cnodenumbers(:)
      integer ::             lxyzcs1, lxyzcs2
      real(kr),allocatable :: xyzcs(:,:)
      real(kr) :: init_value

      integer :: inc, indc, indcs
      integer :: ndofcs, nnodcs, pointinetc

      ! data concerning pairs for adaptivity
      integer :: idpair
      character(lfnamex) :: filename
      integer :: npair
      integer ::            lpairs1, lpairs2
      integer,allocatable :: pairs(:,:)
      integer :: ipair, ldata
      integer ::            lpair2proc
      integer,allocatable :: pair2proc(:)

      logical :: remove_original 
      logical :: remove_bc_nodes 
      logical :: keep_global 
      logical :: is_mesh_loaded 

      ! variables for new Parmetis Communicator
      integer ::            lranks
      integer,allocatable :: ranks(:)
      integer ::            lnelempa
      integer,allocatable :: nelempa(:)
      integer ::             comm_parmetis, myid_parmetis, nproc_parmetis, mpi_group_all, mpi_group_parmetis
      logical :: is_comm_parmetis_created
      integer ::             nproc_used, kranks, iproc
      integer ::             nelem_loc_min
      integer ::             lkadjsub_row

      ! variables for matrix distribution among levels
      integer ::            comm_prev, myid_prev, nproc_prev
      integer ::            liets_aux
      integer,allocatable :: iets_aux(:)
      integer ::            lsub2proc_aux
      integer,allocatable :: sub2proc_aux(:)

      logical,parameter :: use_explicit_schurs = .false.

      character(1) :: levelstring
      character(256) :: prefix
      character(300) :: command


      ! time variables
      real(kr) :: t_division, t_matrix_import, t_adjacency, t_loc_mesh,&
                  t_loc_interface, t_loc_bc, t_loc_adjacency,&
                  t_matrix_assembly, t_schur_prepare, t_weights_prepare,&
                  t_reduced_rhs_prepare, t_prepare_c, t_prepare_aug,&
                  t_prepare_coarse, t_standard_coarse_prepare, t_adaptive_coarse_prepare,&
                  t_par_globs_search, t_construct_cnodes 

      if (.not.levels(ilevel-1)%is_level_prepared .and.  levels(ilevel)%nsub_loc .gt. 0 ) then
         call error(routine_name, 'Previous level not ready:', ilevel-1)
      end if

      ! orient in the communicator
      if (levels(ilevel)%i_am_active_in_this_level) then
         comm_all  = levels(ilevel)%comm_all
         comm_self = levels(ilevel)%comm_self
         call MPI_COMM_RANK(comm_all,myid,ierr)
         call MPI_COMM_SIZE(comm_all,nproc,ierr)
      end if

      ! make the connection with previous level
      levels(ilevel)%nelem   = levels(ilevel-1)%nsub
      levels(ilevel)%nnod    = levels(ilevel-1)%nnodc
      levels(ilevel)%ndof    = levels(ilevel-1)%ndofc
      levels(ilevel)%ndim    = levels(ilevel-1)%ndim
      levels(ilevel)%meshdim = levels(ilevel-1)%meshdim
      ! set neighbouring for higher levels than 1 - it is never given
      if ( ilevel.gt.1 ) then
         if      ( levels(ilevel)%meshdim .eq. 3 ) then
            levels(ilevel)%neighbouring = 4
         else if ( levels(ilevel)%meshdim .eq. 2 ) then
            levels(ilevel)%neighbouring = 3
         else if ( levels(ilevel)%meshdim .eq. 1 ) then
            levels(ilevel)%neighbouring = 1
         else
            call error ( routine_name, 'Unsupported mesh dimension :', levels(ilevel)%meshdim )
         end if
      end if
      ! for debugging purposes
      !if (ilevel.gt.1) then
      !   levels(ilevel)%load_division = levels(ilevel-1)%load_division
      !end if

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

      ! debug
      !print *,'nnod',levels(ilevel)%nnod
      !print *,'nelem',levels(ilevel)%nelem
      !print *,'nndf',levels(ilevel)%nndf

      ! initialize values
      nsub      = levels(ilevel)%nsub
      nsub_loc  = levels(ilevel)%nsub_loc
      nnod      = levels(ilevel)%nnod
      nelem     = levels(ilevel)%nelem

      ! only perform division of into subdomains for global loading
      call MPI_BCAST(levels(ilevel)%global_loading, 1, MPI_LOGICAL, 0, levels(ilevel-1)%comm_all, ierr)
      if (.not. levels(ilevel)%global_loading .and. ilevel .eq. 1) then
         goto 1234
      end if

      ! make division into subdomains
      graphtype = 0 ! not weighted
      levels(ilevel)%liets = levels(ilevel)%nelem
      allocate(levels(ilevel)%iets(levels(ilevel)%liets))

      ! use previous communicator for distributing resulting IETS
      if (.not.levels(ilevel)%i_am_active_in_this_level) then
         ! jump to obtain resulting division
         goto 111
      end if

      ! this could be generalized to more levels if needed
      !if (ilevel.eq.1 .and. levels(ilevel)%load_division) then
      if (levels(ilevel)%load_division) then
         ! read the division from file *.ES
         if (myid.eq.0) then
            if (ilevel.lt.10) then
               write(levelstring,'(i1)') ilevel
            else
               call error(routine_name,'Index of level too large for reading from file:',ilevel)
            end if
            filename = 'partition_l'//levelstring//'.ES'
            call allocate_unit(ides)
            open (unit=ides,file=filename,status='old',form='formatted')
            rewind ides

            read(ides,*) levels(ilevel)%iets
            close (ides)
         end if
         if (debug) then
            if (myid.eq.0) then
               call info(routine_name, 'Mesh division loaded from file '//trim(filename))
            end if
         end if
      else ! create partition and not load it from file
         ! distribute the mesh among processors
!-----profile
         if (profile) then
            call MPI_BARRIER(comm_all,ierr)
            call time_start
         end if
!-----profile
         !if (parallel_division.and.nproc.gt.1) then
         if (parallel_division.and.nproc.gt.1.and.ilevel.eq.1) then
            ! TODO: use INET only on root
            ! at the moment, all processors have access to global INET array
            linet => levels(ilevel)%linet
            inet  => levels(ilevel)%inet
            lnnet => levels(ilevel)%lnnet
            nnet  => levels(ilevel)%nnet
            liets => levels(ilevel)%liets
            iets  => levels(ilevel)%iets

            if (nproc.gt.nsub) then
               if (myid.eq.0) then
                  call info(routine_name, 'creating communicator for ParMETIS')
               end if
               call MPI_COMM_GROUP(comm_all,mpi_group_all,ierr)
               lranks = nsub
               allocate(ranks(lranks))
               do isub = 1,nsub
                  ranks(isub) = isub - 1 
               end do
               ! prepare new group of processes
               call MPI_GROUP_INCL(mpi_group_all,nsub,ranks,mpi_group_parmetis,ierr)
               ! create new communicator
               call MPI_COMM_CREATE(comm_all,mpi_group_parmetis,comm_parmetis,ierr)
               is_comm_parmetis_created = .true.
               ! free groups
               call MPI_GROUP_FREE(mpi_group_parmetis,ierr)
               call MPI_GROUP_FREE(mpi_group_all,ierr)

               ! I won't work on partitioning, skip to getting IETS array
               if (.not.any(ranks.eq.myid)) then
                  deallocate(ranks)
                  goto 123
               end if
               deallocate(ranks)
            else
               comm_parmetis = comm_all
               is_comm_parmetis_created = .false.
            end if

            call MPI_COMM_RANK(comm_parmetis,myid_parmetis,ierr)
            call MPI_COMM_SIZE(comm_parmetis,nproc_parmetis,ierr)

            ! first create linear distribution of elements among processors
            lelm2proc = nproc_parmetis + 1
            allocate(elm2proc(lelm2proc))
            call pp_distribute_linearly(nelem,nproc_parmetis,elm2proc,lelm2proc)

            ! find local number of elements in linear distribution and maximal value
            nelem_loc  = elm2proc(myid_parmetis+2) - elm2proc(myid_parmetis+1)
!************************************************************PARALLEL
            call MPI_ALLREDUCE(nelem_loc,nelem_locx,1, MPI_INTEGER, MPI_MAX, comm_parmetis, ierr) 
!************************************************************PARALLEL
            ! debug
            !print *,'myid = ',myid_parmetis,'nelem_loc',nelem_loc,'nelem_locx',nelem_locx
            
            ! distribute NNET array
            lnnet_loc = nelem_loc
            allocate(nnet_loc(lnnet_loc))
            el_start = elm2proc(myid_parmetis+1)
            do ie_loc = 1,nelem_loc
               indel = el_start + ie_loc - 1
               nnet_loc(ie_loc) = nnet(indel)
            end do
      
            ! debug
            !print *, 'myid =',myid_parmetis,'nnet_loc = ',nnet_loc
      
            linet_loc = sum(nnet_loc(1:lnnet_loc))
            allocate(inet_loc(linet_loc))
      
            ! now distribute INET array
            if (myid_parmetis.eq.0) then
               ! root process copies its data and distributes the inet array by messages
               length = sum(nnet_loc)
               inet_loc = inet(1:length)
      
               ! now send messages to others
               indinet = length + 1
               do indproc = 1,nproc_parmetis - 1
                  el_start_send  = indproc * nelem_locx + 1
                  el_finish_send = min((indproc+1) * nelem_locx,nelem)
                  length_send = 0
                  do ie = el_start_send,el_finish_send
                     length_send = length_send + nnet(ie)
                  end do
                  if (length_send.gt.0) then
!************************************************************PARALLEL
                     call MPI_SEND(inet(indinet),length_send,MPI_INTEGER,indproc,indproc,comm_parmetis,ierr)
!************************************************************PARALLEL
                  end if
                  indinet = indinet + length_send
               end do
               ! free memory on root
               nullify(inet)
               nullify(nnet)
            else
               if (linet_loc.gt.0) then
!************************************************************PARALLEL
                  call MPI_RECV(inet_loc,linet_loc,MPI_INTEGER,0,myid_parmetis,comm_parmetis,stat,ierr)
!************************************************************PARALLEL
               end if
            end if
      
            ! debug
            !write(*,*) 'myid =',myid_parmetis,'inet_loc = ',inet_loc

            ! prepare initial iets
            lelm2sub = nsub + 1
            allocate(elm2sub(lelm2sub))
            call pp_distribute_linearly(nelem,nsub,elm2sub,lelm2sub)
            do isub = 1,nsub
               do ie = elm2sub(isub),elm2sub(isub+1)-1
                  iets(ie) = isub
               end do
            end do
            deallocate(elm2sub)
      
            lpart_loc = nelem_loc
            allocate(part_loc(lpart_loc))
      
            ! prepare initial linear distribution of subdomains
            do ie_loc = 1,nelem_loc
               indel  = el_start + ie_loc - 1
               indsub = iets(indel)
      
               part_loc(ie_loc) = indsub
            end do
               
            ! debug
            !print *, 'myid =',myid_parmetis,'part_loc = ',part_loc
            !call flush(6)
      
      ! divide mesh into subdomains by ParMetis
            graphtype = 0 ! no weights
            call pp_pdivide_mesh(myid_parmetis,nproc_parmetis,comm_parmetis,graphtype,&
                                 levels(ilevel)%neighbouring,nelem,nelem_loc,&
                                 inet_loc,linet_loc,nnet_loc,lnnet_loc,nsub,&
                                 edgecut,part_loc,lpart_loc)
            if (myid.eq.0) then
               call info( routine_name, 'Mesh divided. Resulting number of cut edges:',edgecut )
            end if
      
      ! prepare memory for the global array on root
            if (myid_parmetis.eq.0) then
               ! store my data
               do ie_loc = 1,nelem_loc
                  indel  = el_start + ie_loc - 1
                  iets(indel) = part_loc(ie_loc)
               end do

               ! receive data from others processors
               do iproc = 1,nproc_parmetis-1
                  nelem_loc_proc = elm2proc(iproc+2) - elm2proc(iproc+1)
                  el_start_proc  = elm2proc(iproc+1)
                  lpart_loc_proc = nelem_loc_proc
                  allocate(part_loc_proc(lpart_loc_proc))
                  call MPI_RECV(part_loc_proc,lpart_loc_proc,MPI_INTEGER,iproc,iproc,comm_parmetis,stat,ierr)
                  ! store data from processor
                  do ie_loc_proc = 1,nelem_loc_proc
                     indel  = el_start_proc + ie_loc_proc - 1
                     iets(indel) = part_loc_proc(ie_loc_proc)
                  end do
                  deallocate(part_loc_proc)
               end do
               ! debug
               !print *,'iets',iets
            else
               ! send my part_loc to root
               call MPI_SEND(part_loc,lpart_loc,MPI_INTEGER,0,myid_parmetis,comm_parmetis,ierr)
            end if

      ! clear communicator
            if (is_comm_parmetis_created) then
               call MPI_COMM_FREE(comm_parmetis,ierr)
               is_comm_parmetis_created = .false.
            end if
      ! free some memory
            deallocate(part_loc)
            deallocate(inet_loc)
            deallocate(nnet_loc)
            deallocate(elm2proc)
            
      ! free memory
            nullify(inet)
            nullify(linet)
            nullify(nnet)
            nullify(lnnet)
            nullify(iets)
            nullify(liets)
  
  123       continue     
         else ! divide mesh serially on root process
            if (myid.eq.0) then
               call pp_divide_mesh(graphtype,levels_correct_division,levels(ilevel)%neighbouring,&
                                   levels(ilevel)%nelem,levels(ilevel)%nnod,&
                                   levels(ilevel)%inet,levels(ilevel)%linet,levels(ilevel)%nnet,levels(ilevel)%lnnet,nsub,&
                                   edgecut,levels(ilevel)%iets,levels(ilevel)%liets)
            end if 
         end if
         if (debug) then
            if (myid.eq.0) then
               call info(routine_name, 'Mesh division created.')
            end if 
         end if 
         if (damp_division) then
            if (myid.eq.0) then
               call levels_damp_division(ilevel,levels(ilevel)%iets,levels(ilevel)%liets)
            end if
         end if
!-----profile
         if (profile) then
            call MPI_BARRIER(comm_all,ierr)
            call time_end(t_division)
            if (myid.eq.0) then
               call time_print('creating division into subdomains',t_division)
            end if
         end if
!-----profile
      end if
      ! populate IETS along previous communicator from root processor
!************************************************************PARALLEL
      call MPI_BCAST(levels(ilevel)%iets,levels(ilevel)%liets, MPI_INTEGER, 0, comm_all, ierr)
!************************************************************PARALLEL
      ! check division - that number of subdomains equal largest entry in partition array IETS
      !if (maxval(levels(ilevel)%iets).ne.nsub) then
      !   !write(*,*) 'IETS:',levels(ilevel)%iets
      !   call warning(routine_name,'Partition does not contain all subdomains.')
      !end if

!-----profile 
      if (profile) then
         call MPI_BARRIER(comm_all,ierr)
         call time_start
      end if
!-----profile

      ! create subdomain mesh
      do isub_loc = 1,nsub_loc
         isub = levels(ilevel)%indexsub(isub_loc)
         ! check division - that number of subdomains equal largest entry in partition array IETS
         if (.not.any(levels(ilevel)%iets.eq.isub)) then
            !write(*,*) 'IETS:',levels(ilevel)%iets
            call warning(routine_name,'Partition does not contain subdomain',isub)
         end if
         call dd_localize_mesh(levels(ilevel)%subdomains(isub_loc),isub,levels(ilevel)%ndim,&
                               levels(ilevel)%nelem,levels(ilevel)%nnod,&
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
            call time_print('localizing mesh',t_loc_mesh)
         end if
      end if
!-----profile

!-----profile 
      if (profile) then
         call MPI_BARRIER(comm_all,ierr)
         call time_start
      end if
!-----profile

      ! create subdomain BC, RHS, and SOL on first level
      if (ilevel.eq.1) then
         call levels_localize_new_data
      end if

!-----profile
      if (profile) then
         call MPI_BARRIER(comm_all,ierr)
         call time_end(t_loc_bc)
         if (myid.eq.0) then
            call time_print('localizing boundary conditions and right-hand side',t_loc_bc)
         end if
      end if
!-----profile
      
! only for global loading
 1234 continue

      ! export data to ParaView for debugging purposes
      if (plot_inputs) then
         if (ilevel.lt.10) then
            write(levelstring,'(i1)') ilevel
         else
            call error(routine_name,'Index of level too large for writing to string:',ilevel)
         end if
         ! set the names for files
         prefix = 'bddcml_data_on_level_'//levelstring

         ! make sure the directory exists
         if (myid.eq.0) then
            command = 'mkdir -p '//trim(prefix)
            call system(trim(command))
         end if
         call MPI_BARRIER(comm_all,ierr)

         ! write local files
         do i = 1,levels(ilevel)%nsub_loc
            call dd_plot_subdomain_data_vtu( prefix, levels(ilevel)%subdomains(i) )
         end do

         ! export the umbrella PVD file
         if (myid.eq.0) then
            call paraview_export_pvd_file(prefix, nsub)
         end if
      end if

      ! check that subdomains are properly loaded
      do i = 1,levels(ilevel)%nsub_loc
         call dd_is_mesh_loaded( levels(ilevel)%subdomains(i), is_mesh_loaded )
         if ( .not. is_mesh_loaded ) then
            call error( routine_name, 'It appears that not all expected subdomains were loaded on proc: ', myid )
         end if
      end do
      ! assure that all processes has same value of ndim as root
!**************************************************************MPI
      call MPI_BCAST(levels(ilevel)%ndim,     1,  MPI_INTEGER, 0, comm_all, ierr)
      call MPI_BCAST(levels(ilevel)%meshdim,  1,  MPI_INTEGER, 0, comm_all, ierr)
!**************************************************************MPI

      ! for communication, any shared nodes are considered
!-----profile 
      if (profile) then
         call MPI_BARRIER(comm_all,ierr)
         call time_start
      end if
!-----profile
      neighball = 1
      if (levels_parallel_neighbouring)then 
         if (levels_parmetis_neighbouring) then
            nelem_loc  = 0
            do isub_loc = 1,nsub_loc
               ! find size of subdomain and mesh
               call dd_get_size(levels(ilevel)%subdomains(isub_loc), ndofs,nnods,nelems)

               ! add counters
               nelem_loc = nelem_loc + nelems
            end do
            lnelempa = nproc
            allocate(nelempa(lnelempa))
!***************************************************************PARALLEL
            call MPI_ALLGATHER(nelem_loc,1,MPI_INTEGER,nelempa,1,MPI_INTEGER,comm_all,ierr)
!***************************************************************PARALLEL

            nelem_loc_min = minval(nelempa)
            if (nelem_loc_min.eq.0) then
               if (myid.eq.0) then
                  call info(routine_name, 'creating communicator for ParMETIS')
               end if
               call MPI_COMM_GROUP(comm_all,mpi_group_all,ierr)
               lranks = nproc
               allocate(ranks(lranks))
               kranks = 0
               do iproc = 0,nproc-1
                  if (nelempa(iproc+1).gt.0 ) then
                     kranks = kranks + 1
                     ranks(kranks) = iproc
                  end if
               end do
               nproc_used = kranks
               !print *,'ranks',ranks(1:kranks)
               !print *,'kranks',kranks
               !call flush(6)
               ! prepare new group of processes
               call MPI_GROUP_INCL(mpi_group_all,nproc_used,ranks(1:nproc_used),mpi_group_parmetis,ierr)
               ! create new communicator
               call MPI_COMM_CREATE(comm_all,mpi_group_parmetis,comm_parmetis,ierr)
               is_comm_parmetis_created = .true.
               ! free groups
               call MPI_GROUP_FREE(mpi_group_parmetis,ierr)
               call MPI_GROUP_FREE(mpi_group_all,ierr)

               ! I won't work on partitioning, skip to getting IETS array
               if (.not.any(ranks(1:nproc_used).eq.myid)) then
                  deallocate(nelempa)
                  deallocate(ranks)
                  goto 124
               end if
               deallocate(ranks)
            else
               comm_parmetis = comm_all
               is_comm_parmetis_created = .false.
            end if
            deallocate(nelempa)

            call MPI_COMM_RANK(comm_parmetis,myid_parmetis,ierr)
            call MPI_COMM_SIZE(comm_parmetis,nproc_parmetis,ierr)

            lkadjsub = nsub * nsub_loc
            allocate(kadjsub(lkadjsub))
            kadjsub = 0

            lnelemsa = nsub ! number of elements in subdomains
            allocate(nelemsaaux(lnelemsa),nelemsa(lnelemsa))
            nelemsaaux = 0 ! number of elements in subdomains

            linet_loc  = 0 ! length of local INET
            lnnet_loc  = 0 ! length of local NNET
            nelem_loc  = 0
            do isub_loc = 1,nsub_loc
               isub = levels(ilevel)%indexsub(isub_loc)

               ! find size of subdomain and mesh
               call dd_get_size(levels(ilevel)%subdomains(isub_loc), ndofs,nnods,nelems)
               call dd_get_mesh_basic_size(levels(ilevel)%subdomains(isub_loc),linets,lnnets)

               ! add counters
               linet_loc = linet_loc + linets
               lnnet_loc = lnnet_loc + lnnets
               nelem_loc = nelem_loc + nelems

               nelemsaaux(isub) = nelems
            end do
!***************************************************************PARALLEL
            call MPI_ALLREDUCE(nelemsaaux, nelemsa, lnelemsa,MPI_INTEGER, MPI_SUM, comm_parmetis, ierr)
!***************************************************************PARALLEL
            deallocate(nelemsaaux)
            ! debug
            !write(*,*) 'nelemsaaux',nelemsaaux
            !write(*,*) 'nelemsa',nelemsa

            ! get local INET and NNET
            ! take care of zero arrays
            linet_loc = max(1,linet_loc)
            lnnet_loc = max(1,lnnet_loc)
            allocate(nnet_loc(lnnet_loc),inet_loc(linet_loc))
            nnet_loc = 0
            inet_loc = 0
            kinet_loc = 1
            knnet_loc = 1
            do isub_loc = 1,nsub_loc
               isub = levels(ilevel)%indexsub(isub_loc)

               ! find size of subdomain mesh
               call dd_get_mesh_basic_size(levels(ilevel)%subdomains(isub_loc),linets,lnnets)

               ! debug
               !print *,'myid =',myid
               !print *,'nsub_loc =',nsub_loc
               !print *,'linet_loc =',linet_loc
               !print *,'kinet_loc =',kinet_loc
               !print *,'linets =',linets
               !print *,'lnnets =',lnnets
               call flush(6)

               ! get portion of INET and NNET arrays
               use_global_indices = .true.
               if (linets.gt.0) then
                  call dd_get_mesh_basic(levels(ilevel)%subdomains(isub_loc),use_global_indices,&
                                         inet_loc(kinet_loc),linets,nnet_loc(knnet_loc),lnnets)
               end if

               kinet_loc = kinet_loc + linets
               knnet_loc = knnet_loc + lnnets
            end do
            if (nelem_loc.eq.0) then
               inet_loc = 0
               inet_loc = 0
            end if
            ! create fake IETS array with elements ordered subdomain by subdomain
            liets_linear = nelem
            allocate(iets_linear(liets_linear))
            indiets = 0
            do isub = 1,nsub
               do ie = 1,nelemsa(isub)
                  indiets = indiets + 1
                  iets_linear(indiets) = isub
               end do
            end do
            !print *, 'nelemsa',nelemsa
            !print *, 'iets_linear',iets_linear
            deallocate(nelemsa)

            sub_start = levels(ilevel)%sub2proc(myid+1)
            call pp_pget_sub_neighbours(myid_parmetis,nproc_parmetis,comm_parmetis,neighball,nelem,nelem_loc,nsub,&
                                        nsub_loc,sub_start,&
                                        inet_loc,linet_loc,nnet_loc,lnnet_loc, &
                                        iets_linear, liets_linear, &
                                        debug, kadjsub,lkadjsub)
         ! clear communicator
            if (is_comm_parmetis_created) then
               call MPI_COMM_FREE(comm_parmetis,ierr)
               is_comm_parmetis_created = .false.
            end if
            deallocate(nnet_loc,inet_loc)
            deallocate(iets_linear)
 124        continue
         else
            ! use guess on integer ranges
            lkadjsub = nsub * nsub_loc
            allocate(kadjsub(lkadjsub))
            kadjsub = 0
            call dd_guess_neighbouring(levels(ilevel)%ndim,levels(ilevel)%subdomains,levels(ilevel)%lsubdomains, &
                                       levels(ilevel)%sub2proc,levels(ilevel)%lsub2proc, &
                                       levels(ilevel)%indexsub,levels(ilevel)%lindexsub, comm_all, &
                                       kadjsub,lkadjsub )
            !print *,'myid',myid,'kadjsub',kadjsub
         end if
      else
         lkadjsub = max(nsub*nsub,1)
         allocate(kadjsub(lkadjsub))
         kadjsub = 0
         if (myid.eq.0) then
            call pp_get_sub_neighbours(neighball,nelem,nnod,nsub,&
                                       levels(ilevel)%inet,levels(ilevel)%linet,&
                                       levels(ilevel)%nnet,levels(ilevel)%lnnet,&
                                       levels(ilevel)%iets,levels(ilevel)%liets,&
                                       kadjsub,lkadjsub)
         end if
         ! check zeros on diagonal
         do i = 1,nsub
            if (kadjsub((i-1)*nsub + i) .ne. 0) then
               call error(routine_name,'nonzero on diagonal of adjacency matrix')
            end if
         end do
         ! check symmetry
         do i = 1,nsub
            do j = i+1,nsub
               if (kadjsub((i-1)*nsub + j) .ne. kadjsub((j-1)*nsub + i)) then
                  call error(routine_name,'adjacency matrix not symmetric')
               end if
            end do
         end do
!***************************************************************PARALLEL
         call MPI_BCAST(kadjsub,lkadjsub, MPI_INTEGER, 0, comm_all, ierr)
!***************************************************************PARALLEL
      end if

      ! debug
      !write(*,*) 'kadjsub'
      !do isub_loc = 1,nsub_loc
      !   if (levels_parallel_neighbouring) then
      !      pkadjsub = (isub_loc-1)*nsub
      !   else
      !      isub = levels(ilevel)%indexsub(isub_loc)
      !      pkadjsub = (isub-1)*nsub
      !   end if
      !   write(*,*) (kadjsub(pkadjsub+j),j = 1,nsub)
      !end do
!-----profile
      if (profile) then
         call MPI_BARRIER(comm_all,ierr)
         call time_end(t_adjacency)
         if (myid.eq.0) then
            call time_print('finding neighbours',t_adjacency)
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
      do isub_loc = 1,nsub_loc
         if (levels_parallel_neighbouring) then
            pkadjsub = (isub_loc-1)*nsub
         else
            isub = levels(ilevel)%indexsub(isub_loc)
            pkadjsub = (isub-1)*nsub
         end if
         if (.not.allocated(kadjsub)) then
            lkadjsub_row = 0
         else
            lkadjsub_row = nsub
         end if
         call dd_localize_adj(levels(ilevel)%subdomains(isub_loc),nsub,kadjsub(pkadjsub+1),lkadjsub_row)
      end do
      if (allocated(kadjsub)) then
         deallocate(kadjsub)
      end if

!-----profile
      if (profile) then
         call MPI_BARRIER(comm_all,ierr)
         call time_end(t_loc_adjacency)
         if (myid.eq.0) then
            call time_print('localizing neighbours',t_loc_adjacency)
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
            call time_print('generating neighbouring data',t_loc_interface)
         end if
      end if
!-----profile

      ! PARALLEL IDENTIFICATION OF GLOBS
!-----profile 
      if (profile) then
         call MPI_BARRIER(comm_all,ierr)
         call time_start
      end if
!-----profile

      ! generate globs on level
      ! on the first level, remove nodes at boundary consitions from globs
      if (ilevel.eq.1) then
         remove_bc_nodes = .true.
      else
         remove_bc_nodes = .false.
      end if
      call dd_create_globs(levels(ilevel)%subdomains,levels(ilevel)%lsubdomains,&
                           levels(ilevel)%sub2proc,levels(ilevel)%lsub2proc,&
                           levels(ilevel)%indexsub,levels(ilevel)%lindexsub, &
                           comm_all,remove_bc_nodes, &
                           levels(ilevel)%meshdim, &
                           ncorner,nedge,nface)
!      call dd_create_globs2(levels(ilevel)%subdomains,levels(ilevel)%lsubdomains,&
!                           levels(ilevel)%sub2proc,levels(ilevel)%lsub2proc,&
!                           levels(ilevel)%indexsub,levels(ilevel)%lindexsub, &
!                           comm_all,remove_bc_nodes, &
!                           .false., ilevel, &
!                           levels(ilevel)%meshdim, &
!                           ncorner,nedge,nface)
      nnodc = ncorner + nedge + nface
!-----profile
      if (profile) then
         call MPI_BARRIER(comm_all,ierr)
         call time_end(t_par_globs_search)
         if (myid.eq.0) then
            call time_print('generating globs',t_par_globs_search)
         end if
      end if
!-----profile

      ! debug
      !write(*,*) 'nnodc, ncorner, nedge, nface'
      !write(*,*) nnodc, ncorner, nedge, nface
      !write(*,*) 'inodc',inodc
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
      do isub_loc = 1,nsub_loc
         isub = levels(ilevel)%indexsub(isub_loc)
         
         call dd_construct_cnodes(levels(ilevel)%subdomains(isub_loc))
         call dd_get_coarse_size(levels(ilevel)%subdomains(isub_loc),ndofc,nnodcs)

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
         call time_end(t_construct_cnodes)
         if (myid.eq.0) then
            call time_print('constructing coarse nodes',t_construct_cnodes)
         end if
      end if
!-----profile

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

! prepare coordinates for coarse mesh
      lxyzc1 = nnodc
      lxyzc2 = levels(ilevel)%ndim
      allocate(xyzc(lxyzc1,lxyzc2))
      allocate(xyzcaux(lxyzc1,lxyzc2))
      ! initialize array
      init_value = -Huge(init_value)
      xyzcaux = init_value

      do isub_loc = 1,nsub_loc
         isub = levels(ilevel)%indexsub(isub_loc)

         call dd_get_coarse_size(levels(ilevel)%subdomains(isub_loc),ndofcs,nnodcs)

         lcnodenumbers = nnodcs
         allocate(cnodenumbers(lcnodenumbers))
         call dd_get_coarse_cnodes(levels(ilevel)%subdomains(isub_loc), cnodenumbers,lcnodenumbers)
         
! create subdomain description of coarse MESH
         pointinetc = kinetc(isub)
         inetcaux(pointinetc+1:pointinetc+nnodcs) = cnodenumbers
         pointinetc = pointinetc + nnodcs

         lxyzcs1 = nnodcs
         lxyzcs2 = levels(ilevel)%ndim
         allocate(xyzcs(lxyzcs1,lxyzcs2))

         call dd_get_coarse_coordinates(levels(ilevel)%subdomains(isub_loc), xyzcs,lxyzcs1,lxyzcs2)

         ! copy local coordinates into global array
         do indcs = 1,nnodcs
            indc = cnodenumbers(indcs)
            if (xyzcaux(indc,1) .eq. init_value) then
               ! I am the first one to write local coordinates to the global
               xyzcaux(indc,1:levels(ilevel)%ndim) = xyzcs(indcs,1:levels(ilevel)%ndim)
            else if (any(abs(xyzcaux(indc,1:levels(ilevel)%ndim) - xyzcs(indcs,1:levels(ilevel)%ndim)).gt.numerical_zero)) then

               ! debug
               !print *,'myid',myid,'isub',isub
               !print *,'cnodenumbers',cnodenumbers
               !print *,'xyzcs'
               !do i = 1,nnodcs
               !   print *,xyzcs(i,1:levels(ilevel)%ndim)
               !end do
               !call flush(6)
               !print *,xyzcaux(indc,1:levels(ilevel)%ndim) , 'vs.', xyzcs(indcs,1:levels(ilevel)%ndim)
               call error(routine_name,'Different coarse coordinates from two subdomains.')
            end if
         end do

         deallocate(xyzcs)
         deallocate(cnodenumbers)

! Finish loop over subdomains
      end do
      deallocate(kinetc)
      allocate(inetc(linetc))
!*****************************************************************MPI
      call MPI_ALLREDUCE(inetcaux,inetc,linetc, MPI_INTEGER, MPI_SUM, comm_all, ierr) 
      call MPI_ALLREDUCE(xyzcaux,xyzc,lxyzc1*lxyzc2, MPI_DOUBLE_PRECISION, MPI_MAX, comm_all, ierr) 
!*****************************************************************MPI
      deallocate(inetcaux)
      deallocate(xyzcaux)
      ! debug
      ! print *,'myid =',myid,'inetc',inetc
      ! check the inetc array
      if (any(inetc.eq.0)) then
         call error(routine_name,'Zeros in inetc array.')
      end if

      ! Get matrix
      if (.not. levels(ilevel)%global_loading .and. ilevel .eq. 1) then
         goto 1235
      end if

!-----profile
      if (profile) then
         call MPI_BARRIER(comm_all,ierr)
         call time_start
      end if
!-----profile
 111  continue
      ! IMPORTANT: at levels > 1, this routine must be run by all members of previous level
      if (ilevel.eq.1) then
         ! uncomment this to use independent subdomain files
         !do isub_loc = 1,nsub_loc
         !   call dd_read_matrix_from_file(levels(ilevel)%subdomains(isub_loc),matrixtype,trim(problemname))
         !end do
         ! use global file for input of element matrices
         call dd_read_matrix_by_root(levels(ilevel)%subdomains,levels(ilevel)%lsubdomains, comm_all,levels(ilevel)%idelm,&
                                     levels(ilevel)%nsub,levels(ilevel)%nelem,matrixtype,&
                                     levels(ilevel)%sub2proc,levels(ilevel)%lsub2proc,&
                                     levels(ilevel)%indexsub,levels(ilevel)%lindexsub,&
                                     levels(ilevel)%iets,levels(ilevel)%liets)
      else
         ! this is a tricky routine - it needs to be run at all processes at PREVIOUS level to properly distribute matrices
         ! prepare auxiliary arrays necessary by the routine
         comm_prev = levels(ilevel-1)%comm_all
         call MPI_COMM_RANK(comm_prev,myid_prev,ierr)
         call MPI_COMM_SIZE(comm_prev,nproc_prev,ierr)
         ! distribute IETS along all processes at previous level
         if (myid_prev.eq.0) then 
            liets_aux = levels(ilevel)%liets
         end if
!**************************************************************MPI
         call MPI_BCAST(liets_aux,  1,  MPI_INTEGER, 0, comm_prev, ierr)
!**************************************************************MPI
         allocate(iets_aux(liets_aux))
         if (myid_prev.eq.0) then 
            ! copy IETS at root
            do i = 1,liets_aux
               iets_aux(i) = levels(ilevel)%iets(i)
            end do
         end if
!**************************************************************MPI
         call MPI_BCAST(iets_aux,  liets_aux,  MPI_INTEGER, 0, comm_prev, ierr)
!**************************************************************MPI
         ! distribute SUB2PROC along all processes at previous level
         if (myid_prev.eq.0) then 
            lsub2proc_aux = nproc_prev+1
         end if
!**************************************************************MPI
         call MPI_BCAST(lsub2proc_aux,  1,  MPI_INTEGER, 0, comm_prev, ierr)
!**************************************************************MPI
         allocate(sub2proc_aux(lsub2proc_aux))
         if (myid_prev.eq.0) then 
            ! copy SUB2PROC at root
            do i = 1,levels(ilevel)%lsub2proc
               sub2proc_aux(i) = levels(ilevel)%sub2proc(i)
            end do
            ! make auxiliary record to last parts of sub2proc
            do i = levels(ilevel)%lsub2proc,levels(ilevel-1)%lsub2proc
               sub2proc_aux(i) = nsub+1
            end do
         end if
!**************************************************************MPI
         call MPI_BCAST(sub2proc_aux,  lsub2proc_aux,  MPI_INTEGER, 0, comm_prev, ierr)
!**************************************************************MPI
         call dd_gather_matrix(levels(ilevel)%subdomains,levels(ilevel)%lsubdomains,&
                               levels(ilevel-1)%subdomains,levels(ilevel-1)%lsubdomains,&
                               comm_prev,levels(ilevel)%nsub,levels(ilevel-1)%nsub,matrixtype,&
                               sub2proc_aux,lsub2proc_aux,&
                               levels(ilevel)%indexsub,levels(ilevel)%lindexsub,&
                               levels(ilevel-1)%sub2proc,levels(ilevel-1)%lsub2proc,&
                               levels(ilevel-1)%indexsub,levels(ilevel-1)%lindexsub,&
                               iets_aux,liets_aux)
         deallocate(iets_aux)
         deallocate(sub2proc_aux)
      end if
      ! IMPORTANT: at levels > 1, previous routine must be run by all members of previous level, 
      !            processors not active at current level can now exit
      if (.not.levels(ilevel)%i_am_active_in_this_level) then
         return
      end if

!-----profile
      if (profile) then
         call MPI_BARRIER(comm_all,ierr)
         call time_end(t_matrix_import)
         if (myid.eq.0) then
            call time_print('importing matrix',t_matrix_import)
         end if
      end if
!-----profile

!-----for local data input, continue here (matrices are loaded by user)
 1235 continue

      if (ilevel.eq.1) then
         call dd_fix_constraints(levels(ilevel)%subdomains,levels(ilevel)%lsubdomains, comm_all,&
                                 levels(ilevel)%sub2proc,levels(ilevel)%lsub2proc,&
                                 levels(ilevel)%indexsub,levels(ilevel)%lindexsub)
      end if

      ! Assembly matrix for each subdomain
!-----profile
      if (profile) then
         call MPI_BARRIER(comm_all,ierr)
         call time_start
      end if
!-----profile
      do isub_loc = 1,nsub_loc
         call dd_assembly_local_matrix(levels(ilevel)%subdomains(isub_loc))
      end do
!-----profile
      if (profile) then
         call MPI_BARRIER(comm_all,ierr)
         call time_end(t_matrix_assembly)
         if (myid.eq.0) then
            call time_print('assembling matrices',t_matrix_assembly)
         end if
      end if
!-----profile

      ! For first level, prepare Schur complements
!-----profile
      if (profile) then
         call MPI_BARRIER(comm_all,ierr)
         call time_start
      end if
!-----profile
      ! Schur only for first level
      remove_original = .false.
      do isub_loc = 1,nsub_loc
         call dd_matrix_tri2blocktri(levels(ilevel)%subdomains(isub_loc),remove_original)
         call dd_prepare_schur(levels(ilevel)%subdomains(isub_loc),comm_self)
      end do
!-----profile
      if (profile) then
         call MPI_BARRIER(comm_all,ierr)
         call time_end(t_schur_prepare)
         if (myid.eq.0) then
            call time_print('preparing Schur complement matrices',t_schur_prepare)
         end if
      end if
!-----profile

      ! BDDC data
      do isub_loc = 1,nsub_loc
         ! load arithmetic averages on corners
         glbtype = 3
         call dd_load_arithmetic_constraints(levels(ilevel)%subdomains(isub_loc),glbtype)
         ! load arithmetic averages on edges
         if (use_arithmetic_constraints) then
            glbtype = 2
            call dd_load_arithmetic_constraints(levels(ilevel)%subdomains(isub_loc),glbtype)
         end if
         ! load arithmetic averages on faces if adaptivity is not active
         if (.not.use_adaptive_constraints) then
            glbtype = 1
            call dd_load_arithmetic_constraints(levels(ilevel)%subdomains(isub_loc),glbtype)
         end if
         ! prepare user constraints on faces on the first level
         if (use_user_constraints .and. ilevel.eq.1) then
            glbtype = 1
            call dd_load_user_constraints(levels(ilevel)%subdomains(isub_loc),glbtype)
            call dd_orthogonalize_constraints(levels(ilevel)%subdomains(isub_loc),glbtype)
         end if
      end do

! prepare nndfc and embed cnodes and get array with numbers of constraints
      ! in nndfc, nodes are ordered as corners - edges - faces
      lnndfc   = nnodc
      allocate (nndfc(lnndfc))
      call zero(nndfc,lnndfc)
      call dd_embed_cnodes(levels(ilevel)%subdomains,levels(ilevel)%lsubdomains, &
                           levels(ilevel)%indexsub,levels(ilevel)%lindexsub,& 
                           comm_all, nndfc,lnndfc)

      ! prepare matrix C for corners and arithmetic averages on edges
!-----profile
      if (profile) then
         call MPI_BARRIER(comm_all,ierr)
         call time_start
      end if
!-----profile
      do isub_loc = 1,nsub_loc
         call dd_prepare_c(levels(ilevel)%subdomains(isub_loc))
      end do
!-----profile
      if (profile) then
         call MPI_BARRIER(comm_all,ierr)
         call time_end(t_prepare_c)
         if (myid.eq.0) then
            call time_print('preparing C',t_prepare_c)
         end if
      end if
!-----profile
      ! prepare augmented matrix for BDDC
!-----profile
      if (profile) then
         call MPI_BARRIER(comm_all,ierr)
         call time_start
      end if
!-----profile
      do isub_loc = 1,nsub_loc
         call dd_prepare_aug(levels(ilevel)%subdomains(isub_loc),comm_self)
      end do
!-----profile
      if (profile) then
         call MPI_BARRIER(comm_all,ierr)
         call time_end(t_prepare_aug)
         if (myid.eq.0) then
            call time_print('preparing augmented problem',t_prepare_aug)
         end if
      end if
!-----profile
      ! prepare coarse space basis functions for BDDC
!-----profile
      if (profile) then
         call MPI_BARRIER(comm_all,ierr)
         call time_start
      end if
      if (ilevel.eq.1) then
         keep_global = .false.
      else
         keep_global = .true.
      end if
!-----profile
      do isub_loc = 1,nsub_loc
         call dd_prepare_coarse(levels(ilevel)%subdomains(isub_loc),keep_global)
      end do
!-----profile
      if (profile) then
         call MPI_BARRIER(comm_all,ierr)
         call time_end(t_prepare_coarse)
         if (myid.eq.0) then
            call time_print('preparing coarse probem matrices and shape functions ',t_prepare_coarse)
         end if
      end if
!-----profile

!-----profile
      if (profile) then
         call MPI_BARRIER(comm_all,ierr)
         call time_end(t_standard_coarse_prepare)
         if (myid.eq.0) then
            call time_print('preparing standard coarse problem',t_standard_coarse_prepare)
         end if
      end if
!-----profile

      if (use_adaptive_constraints) then
!-----profile
         if (profile) then
            call MPI_BARRIER(comm_all,ierr)
            call time_start
         end if
!-----profile
         if (ilevel.eq.1 .and. levels_load_pairs) then
            ! read pairs from *.PAIR file
            if (myid.eq.0) then
               filename = 'pairs_l1.PAIR'
               call allocate_unit(idpair)
               open (unit=idpair,file=filename,status='old',form='formatted')
            end if
            if (myid.eq.0) then
               read(idpair,*) npair
            end if
!*****************************************************************MPI
            call MPI_BCAST(npair,1, MPI_INTEGER, 0, comm_all, ierr)
!*****************************************************************MPI
            lpairs1 = npair
            lpairs2 = 3
            allocate(pairs(lpairs1,lpairs2))
            call zero(pairs,lpairs1,lpairs2)
            if (myid.eq.0) then
               do ipair = 1,npair
                  ! first column is associated with processors - initialize it to -1 = no processor assigned
                  read(idpair,*) (pairs(ipair,j), j = 1,3)
               end do
               close(idpair)
            end if
            ldata = lpairs1*lpairs2
!*****************************************************************MPI
            call MPI_BCAST(pairs,lpairs1*lpairs2, MPI_INTEGER, 0, comm_all, ierr)
!*****************************************************************MPI
         else
            ! generate pairs
            lpairs1 = nface
            lpairs2 = 3
            allocate(pairs(lpairs1,lpairs2))
            call dd_create_pairs(levels(ilevel)%subdomains,levels(ilevel)%lsubdomains, &
                                 levels(ilevel)%indexsub,levels(ilevel)%lindexsub,& 
                                 comm_all,&
                                 pairs,lpairs1,lpairs2, npair)
         end if

         ! check number of pairs
         if (npair.ne.nface) then
            call warning(routine_name,'number of found pairs different from number of faces')
         end if

         call adaptivity_init(comm_all,pairs,lpairs1,lpairs2, npair)
         deallocate(pairs)

         lpair2proc = nproc + 1
         allocate(pair2proc(lpair2proc))
         call pp_distribute_linearly(npair,nproc,pair2proc,lpair2proc)

         if (use_explicit_schurs) then
            do isub_loc = 1,nsub_loc
               call dd_prepare_explicit_schur(levels(ilevel)%subdomains(isub_loc))
            end do
         end if
         call adaptivity_solve_eigenvectors(levels(ilevel)%subdomains,levels(ilevel)%lsubdomains, &
                                            levels(ilevel)%sub2proc,levels(ilevel)%lsub2proc,&
                                            levels(ilevel)%indexsub,levels(ilevel)%lindexsub,&
                                            pair2proc,lpair2proc, comm_all, use_explicit_schurs,&
                                            matrixtype, levels(ilevel)%adaptivity_estimate)

         if (use_explicit_schurs) then
            do isub_loc = 1,nsub_loc
               call dd_destroy_explicit_schur(levels(ilevel)%subdomains(isub_loc))
            end do
         end if

         call adaptivity_finalize
         deallocate(pair2proc)

         do isub_loc = 1,nsub_loc
            glbtype = 1
            call dd_orthogonalize_constraints(levels(ilevel)%subdomains(isub_loc),glbtype)
         end do

         call dd_embed_cnodes(levels(ilevel)%subdomains,levels(ilevel)%lsubdomains, &
                              levels(ilevel)%indexsub,levels(ilevel)%lindexsub,& 
                              comm_all, nndfc,lnndfc)

         ! prepare AGAIN BDDC data
         do isub_loc = 1,nsub_loc
            ! prepare matrix C for and updated constraints on faces
            call dd_prepare_c(levels(ilevel)%subdomains(isub_loc))
            ! prepare augmented matrix for BDDC
            call dd_prepare_aug(levels(ilevel)%subdomains(isub_loc),comm_self)
            ! prepare coarse space basis functions for BDDC
            call dd_prepare_coarse(levels(ilevel)%subdomains(isub_loc),keep_global)
         end do
!-----profile
         if (profile) then
            call MPI_BARRIER(comm_all,ierr)
            call time_end(t_adaptive_coarse_prepare)
            if (myid.eq.0) then
               call time_print('preparing adaptive coarse problem',t_adaptive_coarse_prepare)
            end if
         end if
!-----profile
      end if

      ! print the output
      !do isub_loc = 1,nsub_loc
      !   call dd_print_sub(levels(ilevel)%subdomains(isub_loc))
      !end do

      ! weights
!-----profile
      if (profile) then
         call MPI_BARRIER(comm_all,ierr)
         call time_start
      end if
!-----profile

      ! set type of weights
      !if (matrixtype.eq.1) then
      !   ! use diagonal stiffness for SPD problems
      !   !weights_type = 1
      !   weights_type = 4
      !else
      !   if      (all(levels(ilevel)%subdomains%ldof_data .gt. 0)) then
      !      ! use dof data if they are available
      !      weights_type = 3
      !   else if (all(levels(ilevel)%subdomains%lelement_data1 .gt. 0)) then
      !      ! use element data if they are available
      !      weights_type = 2
      !   else
      !      ! use simple cardinality otherwise
      !      weights_type = 0
      !   end if
      !end if
      call dd_weights_prepare(levels(ilevel)%subdomains,levels(ilevel)%lsubdomains, &
                              levels(ilevel)%sub2proc,levels(ilevel)%lsub2proc,&
                              levels(ilevel)%indexsub,levels(ilevel)%lindexsub,&
                              comm_all, weights_type)
!-----profile
      if (profile) then
         call MPI_BARRIER(comm_all,ierr)
         call time_end(t_weights_prepare)
         if (myid.eq.0) then
            call time_print('preparing weights',t_weights_prepare)
         end if
      end if
!-----profile

      if (ilevel.eq.1) then
         ! prepare reduced RHS
!-----profile
         if (profile) then
            call MPI_BARRIER(comm_all,ierr)
            call time_start
         end if
!-----profile
         call dd_prepare_reduced_rhs_all(levels(ilevel)%subdomains,levels(ilevel)%lsubdomains, &
                                         levels(ilevel)%sub2proc,levels(ilevel)%lsub2proc,&
                                         levels(ilevel)%indexsub,levels(ilevel)%lindexsub,& 
                                         comm_all)
!-----profile
         if (profile) then
            call MPI_BARRIER(comm_all,ierr)
            call time_end(t_reduced_rhs_prepare)
            if (myid.eq.0) then
               call time_print('preparing reduced right-hand side',t_reduced_rhs_prepare)
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
      character(*),parameter:: routine_name = 'LEVELS_PREPARE_LAST_LEVEL'

      integer :: ilevel

      integer :: la, nnz
      integer,allocatable :: i_sparse(:),j_sparse(:)
      real(kr),allocatable :: a_sparse(:)

      integer :: ndof

      integer :: mumpsinfo, iparallel

      ! MPI variables
      integer :: comm_all, comm_self, myid, nproc, ierr

      ! time info
      real(kr) :: t_mumps_analysis, t_mumps_factorization

      ! last level has index of number of levels
      ilevel = nlevels

      if (.not.levels(ilevel-1)%is_level_prepared) then
         call error(routine_name, 'Previous level not ready.')
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
      ndof  = levels(ilevel)%ndof

      !write(*,*) 'coarse ndof on myid',myid,'is',ndof

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

      !write(*,*) 'myid =',myid,'la =',la,'nnz =',nnz
      !call sm_print(6, i_sparse, j_sparse, a_sparse, la, nnz)
      !call flush(6)

! Initialize MUMPS
      call mumps_init(mumps_coarse,comm_all,matrixtype)

! Level of information from MUMPS
      mumpsinfo  = 1
      call mumps_set_info(mumps_coarse,mumpsinfo)

! Load matrix to MUMPS
      call mumps_load_triplet_distributed(mumps_coarse,ndof,nnz,i_sparse,j_sparse,a_sparse,nnz)

! Analyze matrix
!-----profile
      if (profile) then
         call MPI_BARRIER(comm_all,ierr)
         call time_start
      end if
!-----profile
      iparallel = 0 ! let MUMPS decide
      call mumps_analyze(mumps_coarse,iparallel)
      if (debug) then
         if (myid.eq.0) then
            call info(routine_name, 'Coarse matrix analyzed.')
         end if
      end if
!-----profile
      if (profile) then
         call MPI_BARRIER(comm_all,ierr)
         call time_end(t_mumps_analysis)
         if (myid.eq.0) then
            call time_print('analysis of the coarse problem',t_mumps_analysis)
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
      if (debug) then
         if (myid.eq.0) then
            call info(routine_name, 'Coarse matrix factorized.')
         end if
      end if
!-----profile
      if (profile) then
         call MPI_BARRIER(comm_all,ierr)
         call time_end(t_mumps_factorization)
         if (myid.eq.0) then
            call time_print('factorization of the coarse problem',t_mumps_factorization)
         end if
      end if
!-----profile

      is_mumps_coarse_ready = .true.

! Clear memory
      deallocate(i_sparse, j_sparse, a_sparse)

end subroutine

!**********************************
subroutine levels_localize_new_data
!**********************************
      use module_utils
      implicit none
      include "mpif.h"

      ! local vars
      character(*),parameter:: routine_name = 'LEVELS_LOCALIZE_NEW_DATA'
      integer :: ilevel
      integer :: nsub_loc, isub_loc
      integer :: nelems, nnods, ndofs

      integer ::             lifixs
      integer,allocatable ::  ifixs(:)

      integer ::            lrhss,   lfixvs,   lsols
      real(kr),allocatable:: rhss(:), fixvs(:), sols(:)

      ilevel = 1
      nsub_loc = levels(ilevel)%nsub_loc
      do isub_loc = 1,nsub_loc
   
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
   
         ! load subdomain RHS
         lrhss = ndofs
         allocate(rhss(lrhss))
         call dd_map_glob_to_sub(levels(ilevel)%subdomains(isub_loc), levels(ilevel)%rhs,levels(ilevel)%lrhs, rhss,lrhss)
         call dd_upload_rhs(levels(ilevel)%subdomains(isub_loc), rhss,lrhss, .true. )
         deallocate(rhss)

         ! load subdomain SOLS
         lsols = ndofs
         allocate(sols(lsols))
         call dd_map_glob_to_sub(levels(ilevel)%subdomains(isub_loc), levels(ilevel)%sol,levels(ilevel)%lsol, sols,lsols)
         call dd_upload_solution(levels(ilevel)%subdomains(isub_loc), sols,lsols)
         deallocate(sols)

      end do

end subroutine

!*****************************************************************
subroutine levels_pc_apply(common_krylov_data,lcommon_krylov_data)
!*****************************************************************
! Apply last level coarse problem to array lvec
      use module_krylov_types_def
      use module_utils
      implicit none

      integer,intent(in)         ::                 lcommon_krylov_data
      type(common_krylov_data_type),intent(inout) :: common_krylov_data(lcommon_krylov_data)

      ! local vars
      character(*),parameter:: routine_name = 'LEVELS_PC_APPLY'
      integer :: ilr
      integer :: comm_all, myid, ierr

! Apply subdomain corrections on first level 
      iactive_level = 1
      comm_all = levels(iactive_level)%comm_all
      ! orient in the communicator
      call MPI_COMM_RANK(comm_all,myid,ierr)
      if (debug .and. myid.eq.0) then
         call info(routine_name,'Applying subdomain correction at level',iactive_level)
      end if
      call levels_corsub_first_level(common_krylov_data,lcommon_krylov_data)

! Sweep levels upwards - make coarse correction and produce coarse residual
      do iactive_level = 2,nlevels-1
         if (levels(iactive_level)%i_am_active_in_this_level) then
            comm_all = levels(iactive_level)%comm_all
            ! orient in the communicator
            call MPI_COMM_RANK(comm_all,myid,ierr)
            if (debug .and. myid.eq.0) then
               call info(routine_name,'Applying subdomain correction at level',iactive_level)
            end if
            call levels_corsub_standard_level(iactive_level)
         end if
      end do

! Solve coarse porblem at last level
      iactive_level = nlevels
      if (levels(iactive_level)%i_am_active_in_this_level) then
         comm_all = levels(iactive_level)%comm_all
         ! orient in the communicator
         call MPI_COMM_RANK(comm_all,myid,ierr)
         if (debug .and. myid.eq.0) then
            call info(routine_name,'Solving problem at last level ',iactive_level)
         end if
         call levels_solve_last_level
      end if

! Sweep levels downward - add coarse correction 
      do ilr = 2,nlevels-1
         ! make a reverse loop
         iactive_level = nlevels+1 - ilr

         if (levels(iactive_level)%i_am_active_in_this_level) then
            comm_all = levels(iactive_level)%comm_all
            ! orient in the communicator
            call MPI_COMM_RANK(comm_all,myid,ierr)
            if (debug .and. myid.eq.0) then
               call info(routine_name,'Adding coarse correction at level',iactive_level)
            end if
            call levels_add_standard_level(iactive_level)
         end if
      end do

! Add coarse correction at the first level
      iactive_level = 1
      comm_all = levels(iactive_level)%comm_all
      ! orient in the communicator
      call MPI_COMM_RANK(comm_all,myid,ierr)
      if (debug .and. myid.eq.0) then
         call info(routine_name,'Adding coarse correction at level',iactive_level)
      end if
      call levels_add_first_level(common_krylov_data,lcommon_krylov_data)

end subroutine


!***************************************************************************
subroutine levels_corsub_first_level(common_krylov_data,lcommon_krylov_data)
!***************************************************************************
! Apply subdomain correction on first level and produce coarse resudual
      use module_krylov_types_def
      use module_utils
      implicit none
      include "mpif.h"

      integer,intent(in)    ::                      lcommon_krylov_data
      type(common_krylov_data_type),intent(inout) :: common_krylov_data(lcommon_krylov_data)


      ! local vars
      character(*),parameter:: routine_name = 'LEVELS_CORSUB_FIRST_LEVEL'
      integer,pointer ::  lresc
      real(kr),pointer ::  resc(:)

      real(kr),allocatable :: rescaux(:)
      integer ::             lrescs
      real(kr),allocatable :: rescs(:)
      integer ::             laux
      real(kr),allocatable :: aux(:)
      integer ::             laux2
      real(kr),allocatable :: aux2(:)

      integer :: ndofs, nnods, nelems, ndofaaugs, ndofcs, nnodcs
      integer :: nsub_loc, isub_loc, i, nrhs
      integer :: ilevel
      logical :: solve_adjoint

      ! MPI vars
      integer :: comm_all, comm_self, myid, ierr


      ! on first level, resudial is collected from individual subdomains
      ilevel = 1

      ! check prerequisites
      if (.not.levels(ilevel)%is_level_prepared) then
         call error(routine_name,'Level is not prepared.')
      end if

      ! set pointers
      lresc => levels(ilevel)%lsolc
      resc  => levels(ilevel)%solc

      ! get communicators
      comm_all  = levels(ilevel)%comm_all
      comm_self = levels(ilevel)%comm_self
      call MPI_COMM_RANK(comm_all,myid,ierr)

      ! get local number of subdomains
      nsub_loc = levels(ilevel)%nsub_loc
      ! check the length
      if (lcommon_krylov_data .ne. nsub_loc) then
         call error(routine_name,'Inconsistent length of data.')
      end if

      ! prepare global residual 
      allocate(rescaux(lresc))
      call zero(rescaux,lresc)

      ! get local contribution to coarse residual
      do isub_loc = 1,nsub_loc

         call dd_get_coarse_size(levels(ilevel)%subdomains(isub_loc),ndofcs,nnodcs)

         laux = common_krylov_data(isub_loc)%lvec_in
         allocate(aux(laux))
         do i = 1,laux
            aux(i) = common_krylov_data(isub_loc)%vec_in(i)
         end do

         ! aux = wi * resi
         call dd_weightsi_apply(levels(ilevel)%subdomains(isub_loc), aux,laux)

         lrescs = ndofcs
         allocate(rescs(lrescs))
         call zero(rescs,lrescs)

         ! rc = phis_dual' * wi * resi
         call dd_phisi_dual_apply(levels(ilevel)%subdomains(isub_loc), aux,laux, rescs,lrescs)

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
         solve_adjoint = .false.
         call dd_solve_aug(levels(ilevel)%subdomains(isub_loc), aux2,laux2, nrhs, solve_adjoint)

         ! get interface part of the vector of preconditioned residual
         call zero(common_krylov_data(isub_loc)%vec_out,common_krylov_data(isub_loc)%lvec_out)
         call dd_map_sub_to_subi(levels(ilevel)%subdomains(isub_loc), aux2,ndofs, &
                                 common_krylov_data(isub_loc)%vec_out,common_krylov_data(isub_loc)%lvec_out)

         deallocate(aux2)
         deallocate(aux)
         deallocate(rescs)
      end do
!*****************************************************************MPI
      call MPI_REDUCE(rescaux,resc,lresc, MPI_DOUBLE_PRECISION, MPI_SUM, 0, comm_all, ierr) 
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
      character(*),parameter:: routine_name = 'LEVELS_CORSUB_STANDARD_LEVEL'
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
      integer ::             lress
      real(kr),allocatable :: ress(:)
      integer ::             lgs
      real(kr),allocatable :: gs(:)
      integer ::             lresaugs
      real(kr),allocatable :: resaugs(:)

      integer :: ndofs, nnods, nelems, ndofaaugs, ndofcs, nnodcs, ndofis, nnodis
      integer :: nsub_loc, isub_loc, i, nrhs, isub
      logical :: solve_adjoint

      ! MPI vars
      integer :: comm_all, comm_self, myid, ierr

      ! check LEVEL
      if (ilevel.le.1 .or. ilevel.ge.nlevels) then
         call error(routine_name,'Illegal index of level ILEVEL')
      end if
      ! check prerequisites
      if (.not.levels(ilevel)%is_level_prepared) then
         call error(routine_name,'Current level is not prepared.')
      end if

      ! orient in communicators
      comm_all  = levels(ilevel)%comm_all
      comm_self = levels(ilevel)%comm_self
      call MPI_COMM_RANK(comm_all,myid,ierr)

      ! set pointers
      lres  => levels(ilevel-1)%lsolc
      res   => levels(ilevel-1)%solc
      lresc => levels(ilevel)%lsolc
      resc  => levels(ilevel)%solc

      ! assure distribution of previous data from root
!*****************************************************************MPI
      call MPI_BCAST(res, lres , MPI_DOUBLE_PRECISION, 0, comm_all, ierr) 
!*****************************************************************MPI

      ! local number of subdomains
      nsub_loc = levels(ilevel)%nsub_loc

      ! get local contribution to coarse residual and load it into all subdomains of the level
      do isub_loc = 1,nsub_loc
         isub = levels(ilevel)%indexsub(isub_loc)

         ! get dimensions
         call dd_get_size(levels(ilevel)%subdomains(isub_loc),ndofs,nnods,nelems)

         lress = ndofs
         allocate(ress(lress))
         call dd_map_glob_to_sub(levels(ilevel)%subdomains(isub_loc), res,lres, ress,lress)

         ! upload residual on subdomain to DD structure as RHS
         call dd_upload_rhs(levels(ilevel)%subdomains(isub_loc), ress,lress, .true.)
         deallocate(ress)
      end do

      ! prepare reduced RHS, i.e. G for all subdomains of the level
      call dd_prepare_reduced_rhs_all(levels(ilevel)%subdomains,levels(ilevel)%lsubdomains, &
                                      levels(ilevel)%sub2proc,levels(ilevel)%lsub2proc,&
                                      levels(ilevel)%indexsub,levels(ilevel)%lindexsub,& 
                                      comm_all)

      ! prepare global coarse residual 
      allocate(rescaux(lresc))
      call zero(rescaux,lresc)

      ! prepare global residual 
      allocate(resaux(lres))
      call zero(resaux,lres)

      do isub_loc = 1,nsub_loc
         isub = levels(ilevel)%indexsub(isub_loc)

         ! get dimensions
         call dd_get_size(levels(ilevel)%subdomains(isub_loc),ndofs,nnods,nelems)
         call dd_get_interface_size(levels(ilevel)%subdomains(isub_loc),ndofis,nnodis)
         call dd_get_coarse_size(levels(ilevel)%subdomains(isub_loc),ndofcs,nnodcs)

         ! load reduced rhs on subdomain
         lgs = ndofis
         allocate(gs(lgs))
         call dd_get_reduced_rhs(levels(ilevel)%subdomains(isub_loc), gs,lgs)

         ! gs = w * gs
         call dd_weightsi_apply(levels(ilevel)%subdomains(isub_loc), gs,lgs)

         ! make new residual on whole subdomain from gs
         ! prepare local residual from condensed RHS
         lress = ndofs
         allocate(ress(lress))
         call zero(ress,lress)

         call dd_map_subi_to_sub(levels(ilevel)%subdomains(isub_loc),gs,lgs,ress,lress)

         lrescs = ndofcs
         allocate(rescs(lrescs))
         call zero(rescs,lrescs)

         ! rc = phis_dual' * wi * resi
         call dd_phis_dual_apply(levels(ilevel)%subdomains(isub_loc), ress,lress, rescs,lrescs)

         ! embed local resc to global one
         call dd_map_subc_to_globc(levels(ilevel)%subdomains(isub_loc), rescs,lrescs, rescaux,lresc)

         ! SUBDOMAIN CORRECTION
         ! prepare array of augmented size
         call dd_get_aug_size(levels(ilevel)%subdomains(isub_loc), ndofaaugs)
         lresaugs = ndofaaugs
         allocate(resaugs(lresaugs))
         call zero(resaugs,lresaugs)
         ! truncate the vector for embedding - zeros at the end
         do i = 1,ndofs
            resaugs(i) = ress(i)
         end do

         nrhs = 1
         solve_adjoint = .false.
         call dd_solve_aug(levels(ilevel)%subdomains(isub_loc), resaugs,lresaugs, nrhs, solve_adjoint)

         ! apply weights
         ! z = wi * z
         call dd_weights_apply(levels(ilevel)%subdomains(isub_loc), resaugs,ndofs)

         ! get global part of the vector of preconditioned residual
         call dd_map_sub_to_glob(levels(ilevel)%subdomains(isub_loc), resaugs,ndofs, resaux,lres)

         deallocate(resaugs)
         deallocate(ress)
         deallocate(gs)
         deallocate(rescs)
      end do
      ! store data only on root
!*****************************************************************MPI
      call MPI_REDUCE(rescaux,resc,lresc, MPI_DOUBLE_PRECISION, MPI_SUM, 0, comm_all, ierr) 
      call MPI_REDUCE(resaux, res, lres , MPI_DOUBLE_PRECISION, MPI_SUM, 0, comm_all, ierr) 
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
      ! right hand side and solution are correct only on root
      call mumps_resolve(mumps_coarse,solc,lsolc)

end subroutine

!*******************************************
subroutine levels_add_standard_level(ilevel)
!*******************************************
! Collect residual on first level and produce coarse resudual
      use module_utils
      implicit none
      include "mpif.h"

      ! index of level
      integer,intent(in) ::  ilevel

      ! local vars
      character(*),parameter:: routine_name = 'LEVELS_ADD_STANDARD_LEVEL'
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

      integer ::             lresos
      real(kr),allocatable :: resos(:)

      real(kr),allocatable :: solaux(:)
      real(kr),allocatable :: solaux2(:)

      integer :: ndofcs, nnodcs, nnods, nelems, ndofs, ndofis, nnodis, ndofos
      integer :: nsub_loc, isub_loc, i, isub

      ! MPI vars
      integer :: comm_all, comm_self, myid, ierr


      ! check LEVEL
      if (ilevel.le.1 .or. ilevel.ge.nlevels) then
         call error(routine_name,'Illegal index of level ILEVEL',ilevel)
      end if
      ! check prerequisites
      if (.not.levels(ilevel)%is_level_prepared) then
         call error(routine_name,'Level is not prepared:',ilevel)
      end if

      ! get communicators
      comm_all  = levels(ilevel)%comm_all
      comm_self = levels(ilevel)%comm_self
      call MPI_COMM_RANK(comm_all,myid,ierr)

      ! set pointers
      lsol => levels(ilevel-1)%lsolc
      sol  => levels(ilevel-1)%solc
      lsolc => levels(ilevel)%lsolc
      solc  => levels(ilevel)%solc

      ! assure correct distribution of coarse solution
!***************************************************************PARALLEL
      call MPI_BCAST(solc, lsolc, MPI_DOUBLE_PRECISION, 0, comm_all, ierr)
      call MPI_BCAST(sol,  lsol, MPI_DOUBLE_PRECISION, 0, comm_all, ierr)
!***************************************************************PARALLEL

      ! local number of subdomains
      nsub_loc = levels(ilevel)%nsub_loc

      ! prepare memory for coarse contribution
      allocate(solaux(lsol))
      allocate(solaux2(lsol))
      call zero(solaux,lsol)
      call zero(solaux2,lsol)

      do isub_loc = 1,nsub_loc

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
         call dd_phis_apply(levels(ilevel)%subdomains(isub_loc), solcs,lsolcs, sols,lsols)
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
      
      ! add corrections - sol already contains subdomain contributions
      ! TODO: this could be done only among neighbouring subdomains exactly as on the first level
      !       using commvec instead of global allreduce

      do i = 1,lsol
         sol(i) = sol(i) + solaux2(i)
      end do
      ! now SOL contains coarse and subdomain corrections

      ! interior POST-CORRECTION
      call zero(solaux,lsol) 
      do isub_loc = 1,nsub_loc
         isub = levels(ilevel)%indexsub(isub_loc)

         call dd_get_size(levels(ilevel)%subdomains(isub_loc), ndofs,nnods,nelems)
         call dd_get_interface_size(levels(ilevel)%subdomains(isub_loc),ndofis,nnodis)
         ndofos = ndofs - ndofis

         lsols  = ndofs
         allocate(sols(lsols))
         call zero(sols,lsols)

         ! restrict global corrections to local sols
         call dd_map_glob_to_sub(levels(ilevel)%subdomains(isub_loc), sol,lsol, sols,lsols)

         ! construct subdomain residual
         lresos  = ndofos
         allocate(resos(lresos))
         call dd_construct_interior_residual(levels(ilevel)%subdomains(isub_loc),&
                                             sols,lsols, resos,lresos)
         ! solve interior problem
         call dd_solve_interior_problem(levels(ilevel)%subdomains(isub_loc), &
                                        resos,lresos)

         ! map interior solution to subdomain
         call zero(sols,lsols)
         call dd_map_subo_to_sub(levels(ilevel)%subdomains(isub_loc), resos,lresos, sols,lsols)
         deallocate(resos)

         ! prolong subdomain sols to global sol
         call dd_map_sub_to_glob(levels(ilevel)%subdomains(isub_loc), sols,lsols, solaux,lsol)
         deallocate(sols)
      end do
!*****************************************************************MPI
      call MPI_ALLREDUCE(solaux, solaux2, lsol , MPI_DOUBLE_PRECISION, MPI_SUM, comm_all, ierr) 
!*****************************************************************MPI

      ! subtract interior correction
      do i = 1,lsol
         sol(i) = sol(i) - solaux2(i)
      end do
      deallocate(solaux2)

      ! add stored solution at interior
      call zero(solaux,lsol) 
      do isub_loc = 1,nsub_loc

         call dd_get_size(levels(ilevel)%subdomains(isub_loc), ndofs,nnods,nelems)

         lsols  = ndofs
         allocate(sols(lsols))

         ! restrict global corrections to local sols
         call dd_map_glob_to_sub(levels(ilevel)%subdomains(isub_loc), sol,lsol, sols,lsols)

         ! add subdomain interior correction
         call dd_add_interior_solution(levels(ilevel)%subdomains(isub_loc), sols,lsols)

         ! apply weights
         call dd_weights_apply(levels(ilevel)%subdomains(isub_loc), sols,lsols)

         ! restrict prolong subdomain sols to global sol
         call dd_map_sub_to_glob(levels(ilevel)%subdomains(isub_loc), sols,lsols, solaux,lsol)
         deallocate(sols)
      end do
!*****************************************************************MPI
      call MPI_REDUCE(solaux, sol, lsol , MPI_DOUBLE_PRECISION, MPI_SUM, 0, comm_all, ierr) 
!*****************************************************************MPI
      deallocate(solaux)

end subroutine

!************************************************************************
subroutine levels_add_first_level(common_krylov_data,lcommon_krylov_data)
!************************************************************************
! Collect residual on first level and produce coarse resudual
      use module_krylov_types_def
      use module_utils
      implicit none
      include "mpif.h"

      integer,intent(in)    ::                      lcommon_krylov_data
      type(common_krylov_data_type),intent(inout) :: common_krylov_data(lcommon_krylov_data)

      ! local vars
      character(*),parameter:: routine_name = 'LEVELS_ADD_FIRST_LEVEL'
      integer,pointer :: lsolc
      real(kr),pointer :: solc(:)
      integer ::             lsolcs
      real(kr),allocatable :: solcs(:)

      integer :: ndofcs, nnodcs
      integer :: nsub_loc, isub_loc
      integer :: ilevel

      ! MPI vars
      integer :: comm_all, comm_self, myid, ierr


      ! on first level, resudial is collected from individual subdomains
      ilevel = 1

      ! check prerequisites
      if (.not.levels(ilevel)%is_level_prepared) then
         call error(routine_name,'Level is not prepared.')
      end if

      ! get communicators
      comm_all  = levels(ilevel)%comm_all
      comm_self = levels(ilevel)%comm_self
      call MPI_COMM_RANK(comm_all,myid,ierr)

      ! set pointers
      lsolc => levels(ilevel)%lsolc
      solc  => levels(ilevel)%solc

      ! assure correct distribution of coarse solution
!***************************************************************PARALLEL
      call MPI_BCAST(solc, lsolc, MPI_DOUBLE_PRECISION, 0, comm_all, ierr)
!***************************************************************PARALLEL

      ! get local number of subdomains
      nsub_loc = levels(ilevel)%nsub_loc
      ! check the length
      if (lcommon_krylov_data .ne. nsub_loc) then
         call error(routine_name,'Inconsistent length of data.')
      end if

      do isub_loc = 1,nsub_loc

         call dd_get_coarse_size(levels(ilevel)%subdomains(isub_loc), ndofcs,nnodcs)

         lsolcs = ndofcs
         allocate(solcs(lsolcs))

         ! restrict global solc to local solcs
         call dd_map_globc_to_subc(levels(ilevel)%subdomains(isub_loc), solc,lsolc, solcs,lsolcs)

         ! COARSE CORRECTION
         ! z_i = z_i + phis_i * uc_i
         call dd_phisi_apply(levels(ilevel)%subdomains(isub_loc), solcs,lsolcs, &
                             common_krylov_data(isub_loc)%vec_out,common_krylov_data(isub_loc)%lvec_out)
         ! apply weights
         ! z = wi * z
         call dd_weightsi_apply(levels(ilevel)%subdomains(isub_loc), &
                                common_krylov_data(isub_loc)%vec_out,common_krylov_data(isub_loc)%lvec_out)

         ! load Z for communication
         call dd_comm_upload(levels(ilevel)%subdomains(isub_loc), &
                             common_krylov_data(isub_loc)%vec_out,common_krylov_data(isub_loc)%lvec_out)

         deallocate(solcs)
      end do
      ! communicate vector Z
      ! Interchange data
      call dd_comm_swapdata(levels(ilevel)%subdomains,levels(ilevel)%lsubdomains,&
                            levels(ilevel)%indexsub,levels(ilevel)%lindexsub,&    
                            levels(ilevel)%sub2proc,levels(ilevel)%lsub2proc,&
                            comm_all)
      ! Download data
      do isub_loc = 1,nsub_loc

         ! get contibution from neighbours
         call dd_comm_download(levels(ilevel)%subdomains(isub_loc), &
                               common_krylov_data(isub_loc)%vec_out,common_krylov_data(isub_loc)%lvec_out)
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
      ! local vars
      character(*),parameter:: routine_name = 'LEVELS_GET_NUMBER_OF_SUBDOMAINS'

      ! check prerequisites
      if (.not.levels(ilevel)%is_level_prepared) then
         call error(routine_name,'Level is not prepared.')
      end if

      nsub     = levels(ilevel)%nsub
      nsub_loc = levels(ilevel)%nsub_loc

end subroutine

!***********************************************************************************
subroutine levels_prepare_interface_initial_data(isub_loc,solis,lsolis,resis,lresis)
!***********************************************************************************
! Subroutine for initialization of solution and residual at interface of subdomain
! module with utility routines
      use module_utils

      implicit none

      integer,intent(in)   :: lsolis
      real(kr),intent(out) ::  solis(lsolis)
      integer,intent(in)   :: lresis
      real(kr),intent(out) ::  resis(lresis)

! local vars
      character(*),parameter:: routine_name = 'LEVELS_PREPARE_INTERFACE_INITIAL_DATA'
      ! for Krylov data, index of level is 1
      integer,parameter :: ilevel = 1

      integer :: isub_loc
      integer :: ndofs, nnods, nelems, ndofis, nnodis

      integer ::             lsols
      real(kr),allocatable :: sols(:)

      ! check prerequisites
      if (.not.levels(ilevel)%is_level_prepared) then
         call error(routine_name,'Level is not prepared.')
      end if
      ! check dimensions
      if (isub_loc.gt.levels(ilevel)%nsub_loc) then
         call error(routine_name,'Dimension mismatch.')
      end if

      ! determine size of subdomain
      call dd_get_size(levels(ilevel)%subdomains(isub_loc),ndofs,nnods,nelems)
      call dd_get_interface_size(levels(ilevel)%subdomains(isub_loc),ndofis,nnodis)

      ! allocate local subdomain solution
      lsols  = ndofs
      allocate(sols(lsols))

      if (levels(ilevel)%use_initial_solution) then
         ! set initial guess
         if ( levels(ilevel)%global_loading ) then
            call dd_map_glob_to_sub(levels(ilevel)%subdomains(isub_loc),levels(ilevel)%sol,levels(ilevel)%lsol,sols,lsols)
         else
            call dd_download_solution(levels(ilevel)%subdomains(isub_loc), sols,lsols)
         end if
      else
         ! set zero initial guess
         ! u_0 = 0
         call zero(sols,lsols)
      end if

      ! set initial guess to satisfy Dirichlet boundary conditions
      call dd_fix_bc(levels(ilevel)%subdomains(isub_loc), sols,lsols)

      ! restrict solution to interface
      call zero(solis,lsolis)
      call dd_map_sub_to_subi(levels(ilevel)%subdomains(isub_loc), sols,lsols, solis,lsolis)
      deallocate(sols)

      ! set initial residual to RHS
      ! res = g
      call dd_get_reduced_rhs(levels(ilevel)%subdomains(isub_loc), resis,lresis)

end subroutine

!*****************************************************************
subroutine levels_sm_apply(common_krylov_data,lcommon_krylov_data)
!*****************************************************************
! Subroutine for multiplication of distributed vector by System Matrix
! module for distributed Krylov data storage
      use module_krylov_types_def
! module with utility routines
      use module_utils

      implicit none

      integer,intent(in)                           :: lcommon_krylov_data
      type (common_krylov_data_type),intent(inout) ::  common_krylov_data(lcommon_krylov_data)

! local vars
      character(*),parameter:: routine_name = 'LEVELS_SM_APPLY'
      ! for Krylov data, index of level is 1
      integer,parameter :: ilevel = 1

      integer :: comm_all
      integer :: ierr, myid
      integer :: isub_loc

      integer :: ncol

      ! check prerequisites
      if (.not.levels(ilevel)%is_level_prepared) then
         call error(routine_name,'Level is not prepared:',ilevel)
      end if
      ! check dimensions
      if (lcommon_krylov_data.ne.levels(ilevel)%nsub_loc) then
         call error(routine_name,'Dimension mismatch.')
      end if

      ! set communicator
      comm_all = levels(ilevel)%comm_all
      call MPI_COMM_RANK(comm_all,myid,ierr)

      ! vec_out = A*vec_in
      ! Upload data
      do isub_loc = 1,levels(ilevel)%nsub_loc

         call zero(common_krylov_data(isub_loc)%vec_out,common_krylov_data(isub_loc)%lvec_out)

         ncol = 1
         call dd_multiply_by_schur(levels(ilevel)%subdomains(isub_loc),&
                                   common_krylov_data(isub_loc)%vec_in,common_krylov_data(isub_loc)%lvec_in, &
                                   common_krylov_data(isub_loc)%vec_out,common_krylov_data(isub_loc)%lvec_out,ncol)

         call dd_comm_upload(levels(ilevel)%subdomains(isub_loc), &
                             common_krylov_data(isub_loc)%vec_out,common_krylov_data(isub_loc)%lvec_out)
      end do

      ! Interchange data
      call dd_comm_swapdata(levels(ilevel)%subdomains,levels(ilevel)%lsubdomains,&
                            levels(ilevel)%indexsub,levels(ilevel)%lindexsub,&    
                            levels(ilevel)%sub2proc,levels(ilevel)%lsub2proc,&
                            comm_all)
      ! Download data
      do isub_loc = 1,levels(ilevel)%nsub_loc

         ! add contibution from neighbours
         call dd_comm_download(levels(ilevel)%subdomains(isub_loc), &
                               common_krylov_data(isub_loc)%vec_out,common_krylov_data(isub_loc)%lvec_out)
      end do
end subroutine

!***************************************************************
subroutine levels_postprocess_solution(krylov_data,lkrylov_data)
!***************************************************************
! Subroutine for postprocessing of data by Krylov iterative method and storing
! solution into DD structure
! module for distributed Krylov data storage
      use module_krylov_types_def
! module for handling subdomain data
      use module_dd
! module with utility routines
      use module_utils
      implicit none
      include "mpif.h"

      integer,intent(in)                        :: lkrylov_data
      type (common_krylov_data_type),intent(in) ::  krylov_data(lkrylov_data)

! local vars
      character(*),parameter:: routine_name = 'LEVELS_POSTPROCESS_SOLUTION'
      ! for Krylov data, index of level is 1
      integer,parameter :: ilevel = 1
      integer :: isub_loc, isub
      integer :: ndofs, nnods, nelems

      integer ::             lsols
      real(kr),allocatable :: sols(:)


      ! check prerequisites
      if (.not.levels(ilevel)%is_level_prepared) then
         call error(routine_name,'Level is not prepared:',ilevel)
      end if
      ! check dimensions
      if (lkrylov_data.ne.levels(ilevel)%nsub_loc) then
         call error(routine_name,'Dimension mismatch.')
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
                                  krylov_data(isub_loc)%vec_in,krylov_data(isub_loc)%lvec_in, &
                                  sols,lsols)

         ! upload solution to DD structure
         call dd_upload_solution(levels(ilevel)%subdomains(isub_loc), sols,lsols)

         deallocate(sols)
      end do

      if (levels(ilevel)%compute_reactions) then
         call dd_compute_reactions(levels(ilevel)%subdomains,levels(ilevel)%lsubdomains, levels(ilevel)%comm_all, &
                                   levels(ilevel)%sub2proc,levels(ilevel)%lsub2proc,&
                                   levels(ilevel)%indexsub,levels(ilevel)%lindexsub)
      end if

end subroutine

!*******************************************************
subroutine levels_dd_download_solution(isub, sols,lsols)
!*******************************************************
! Subroutine for obtaining LOCAL solution from DD structure.
! Only calls the function from DD module.
! module for handling subdomain data
      use module_dd
! module with utility routines
      use module_utils

      implicit none
      include "mpif.h"

      ! global index of subdomain
      integer, intent(in) ::  isub

      ! local solution 
      integer, intent(in) ::  lsols
      real(kr), intent(out) :: sols(lsols)

! local vars
      character(*),parameter:: routine_name = 'LEVELS_DD_DOWNLOAD_SOLUTION'
      ! for Krylov data, index of level is 1
      integer,parameter :: ilevel = 1

      integer :: isub_loc

      ! local index of subdomain
      call get_index(isub+levels_numshift, levels(ilevel)%indexsub,levels(ilevel)%lindexsub, isub_loc)
      if ( isub_loc .eq. -1 ) then
         call error( routine_name, 'Index of subdomain not found among local subdomains: ', isub )
      end if

      ! obtain solution from DD structure
      call dd_download_solution(levels(ilevel)%subdomains(isub_loc), sols,lsols)

end subroutine

!**********************************************
subroutine levels_get_global_solution(sol,lsol)
!**********************************************
! Subroutine for obtaining local solution from LEVELS and gathering global
! solution on root.
! It downloads it from the DD structures
! module for handling subdomain data
      use module_dd
! module with utility routines
      use module_utils

      implicit none
      include "mpif.h"

      ! global solution - allocated only at root
      integer, intent(in) ::  lsol
      real(kr), intent(out) :: sol(lsol)

! local vars
      character(*),parameter:: routine_name = 'LEVELS_GET_GLOBAL_SOLUTION'
      ! for Krylov data, index of level is 1
      integer,parameter :: ilevel = 1

      ! subdomain solution
      integer ::              lsols
      real(kr), allocatable :: sols(:)

      integer ::             lisvgvn
      integer,allocatable ::  isvgvn(:)

      integer ::             lacceptedp
      integer,allocatable ::  acceptedp(:)

      integer :: ndofs, nnods, nelems, isub, isub_loc

      ! MPI variables
      integer :: ierr, myid, nproc, iproc, comm_all
      integer :: i, indv, nsub_locp
      integer :: stat(MPI_STATUS_SIZE)

      ! orient in communicators
      comm_all  = levels(ilevel)%comm_all
      call MPI_COMM_RANK(comm_all,myid,ierr)
      call MPI_COMM_SIZE(comm_all,nproc,ierr)

      ! check that array is properly allocated at root
      if (myid.eq.0) then
         if (lsol.ne.levels(ilevel)%ndof) then
            call error( routine_name, 'Space for global solution not properly allocated on root.')
         end if
         call zero(sol,lsol)
      end if

      do isub_loc = 1,levels(ilevel)%nsub_loc
         isub = levels(ilevel)%indexsub(isub_loc)

         ! determine size of subdomain
         call dd_get_size(levels(ilevel)%subdomains(isub_loc),ndofs,nnods,nelems)

         ! allocate local subdomain solution
         lsols  = ndofs
         allocate(sols(lsols))

         ! obtain solution from DD structure
         call dd_download_solution(levels(ilevel)%subdomains(isub_loc), sols,lsols)

         ! send solution to root
         lisvgvn = ndofs
         allocate(isvgvn(lisvgvn))
         call dd_get_subdomain_isvgvn(levels(ilevel)%subdomains(isub_loc),isvgvn,lisvgvn)
         if (myid.ne.0) then
            ! now send the data to root
            call MPI_SEND(ndofs,1,MPI_INTEGER,0,isub,comm_all,ierr)
            call MPI_SEND(isvgvn,lisvgvn,MPI_INTEGER,0,isub,comm_all,ierr)
            call MPI_SEND(sols,lsols,MPI_DOUBLE_PRECISION,0,isub,comm_all,ierr)
         else
            ! as root, first write my solution, then process others
            do i = 1,ndofs
               indv = isvgvn(i)
               sol(indv) = sols(i)
            end do
         end if
         deallocate(isvgvn)
         deallocate(sols)
      end do

      ! construct solution at root
      if (myid.eq.0) then
         lacceptedp = nproc
         allocate(acceptedp(lacceptedp))
         call zero(acceptedp,lacceptedp)
         acceptedp(1) = levels(ilevel)%nsub_loc
         ! get solutions from others 
         do 
            do iproc = 1,nproc-1
               nsub_locp = levels(ilevel)%sub2proc(iproc+2) - levels(ilevel)%sub2proc(iproc+1)
               if (acceptedp(iproc+1).lt.nsub_locp) then

                  ! I haven't got all processors subdomains yet, get one
                  isub = levels(ilevel)%sub2proc(iproc+1) + acceptedp(iproc+1)
                  call MPI_RECV(ndofs,1,MPI_INTEGER,iproc,isub,comm_all,stat,ierr)
                  lisvgvn = ndofs
                  allocate(isvgvn(lisvgvn))
                  call MPI_RECV(isvgvn,lisvgvn,MPI_INTEGER,iproc,isub,comm_all,stat,ierr)
                  lsols = ndofs
                  allocate(sols(lsols))
                  call MPI_RECV(sols,lsols,MPI_DOUBLE_PRECISION,iproc,isub,comm_all,stat,ierr)

                  ! as root, write received solution
                  do i = 1,ndofs
                     indv = isvgvn(i)

                     sol(indv) = sols(i)
                  end do

                  deallocate(sols)
                  deallocate(isvgvn)
                  ! add counter
                  acceptedp(iproc+1) = acceptedp(iproc+1) + 1
               end if
            end do
            if (all(acceptedp.eq.levels(ilevel)%sub2proc(2:nproc+1) - levels(ilevel)%sub2proc(1:nproc))) then
               exit
            end if
         end do
         deallocate(acceptedp)
      end if

end subroutine

!********************************************************
subroutine levels_dd_download_reactions(isub, reas,lreas)
!********************************************************
! Subroutine for obtaining LOCAL reactions from DD structure.
! Only calls the function from DD module.
! module for handling subdomain data
      use module_dd
! module with utility routines
      use module_utils

      implicit none
      include "mpif.h"

      ! global index of subdomain
      integer, intent(in) ::  isub

      ! local reactions 
      integer, intent(in) ::  lreas
      real(kr), intent(out) :: reas(lreas)

! local vars
      character(*),parameter:: routine_name = 'LEVELS_DD_DOWNLOAD_REACTIONS'
      ! for Krylov data, index of level is 1
      integer,parameter :: ilevel = 1

      integer :: isub_loc

      ! local index of subdomain
      call get_index(isub+levels_numshift, levels(ilevel)%indexsub,levels(ilevel)%lindexsub, isub_loc)
      if ( isub_loc .eq. -1 ) then
         call error( routine_name, 'Index of subdomain not found among local subdomains: ', isub )
      end if

      ! obtain solution from DD structure
      call dd_download_reactions(levels(ilevel)%subdomains(isub_loc), reas,lreas)

end subroutine

!***********************************************
subroutine levels_get_global_reactions(rea,lrea)
!***********************************************
! Subroutine for obtaining local reactions from LEVELS and gathering global
! vector of reactions on root.
! It downloads it from the DD structures
! module for handling subdomain data
      use module_dd
! module with utility routines
      use module_utils

      implicit none
      include "mpif.h"

      ! global reactions - allocated only at root
      integer, intent(in) ::  lrea
      real(kr), intent(out) :: rea(lrea)

! local vars
      character(*),parameter:: routine_name = 'LEVELS_GET_GLOBAL_REACTIONS'
      ! for Krylov data, index of level is 1
      integer,parameter :: ilevel = 1

      ! subdomain reactions
      integer ::              lreas
      real(kr), allocatable :: reas(:)

      integer ::             lisvgvn
      integer,allocatable ::  isvgvn(:)

      integer ::             lacceptedp
      integer,allocatable ::  acceptedp(:)

      integer :: ndofs, nnods, nelems, isub, isub_loc

      ! MPI variables
      integer :: ierr, myid, nproc, iproc, comm_all
      integer :: i, indv, nsub_locp
      integer :: stat(MPI_STATUS_SIZE)

      ! orient in communicators
      comm_all  = levels(ilevel)%comm_all
      call MPI_COMM_RANK(comm_all,myid,ierr)
      call MPI_COMM_SIZE(comm_all,nproc,ierr)

      ! check that array is properly allocated at root
      if (myid.eq.0) then
         if (lrea.ne.levels(ilevel)%ndof) then
            call error( routine_name, 'Space for global reactions not properly allocated on root.')
         end if
         call zero(rea,lrea)
      end if

      do isub_loc = 1,levels(ilevel)%nsub_loc
         isub = levels(ilevel)%indexsub(isub_loc)

         ! determine size of reactions
         call dd_get_size(levels(ilevel)%subdomains(isub_loc),ndofs,nnods,nelems)

         ! allocate local reactions reactions
         lreas  = ndofs
         allocate(reas(lreas))

         ! obtain reactions from DD structure
         call dd_download_reactions(levels(ilevel)%subdomains(isub_loc), reas,lreas)

         ! send reactions to root
         lisvgvn = ndofs
         allocate(isvgvn(lisvgvn))
         call dd_get_subdomain_isvgvn(levels(ilevel)%subdomains(isub_loc),isvgvn,lisvgvn)
         if (myid.ne.0) then
            ! now send the data to root
            call MPI_SEND(ndofs,1,MPI_INTEGER,0,isub,comm_all,ierr)
            call MPI_SEND(isvgvn,lisvgvn,MPI_INTEGER,0,isub,comm_all,ierr)
            call MPI_SEND(reas,lreas,MPI_DOUBLE_PRECISION,0,isub,comm_all,ierr)
         else
            ! as root, first write my reactions, then process others
            do i = 1,ndofs
               indv = isvgvn(i)
               rea(indv) = reas(i)
            end do
         end if
         deallocate(isvgvn)
         deallocate(reas)
      end do

      ! construct reactions at root
      if (myid.eq.0) then
         lacceptedp = nproc
         allocate(acceptedp(lacceptedp))
         call zero(acceptedp,lacceptedp)
         acceptedp(1) = levels(ilevel)%nsub_loc
         ! get reactions from others 
         do 
            do iproc = 1,nproc-1
               nsub_locp = levels(ilevel)%sub2proc(iproc+2) - levels(ilevel)%sub2proc(iproc+1)
               if (acceptedp(iproc+1).lt.nsub_locp) then

                  ! I haven't got all processors subdomains yet, get one
                  isub = levels(ilevel)%sub2proc(iproc+1) + acceptedp(iproc+1)
                  call MPI_RECV(ndofs,1,MPI_INTEGER,iproc,isub,comm_all,stat,ierr)
                  lisvgvn = ndofs
                  allocate(isvgvn(lisvgvn))
                  call MPI_RECV(isvgvn,lisvgvn,MPI_INTEGER,iproc,isub,comm_all,stat,ierr)
                  lreas = ndofs
                  allocate(reas(lreas))
                  call MPI_RECV(reas,lreas,MPI_DOUBLE_PRECISION,iproc,isub,comm_all,stat,ierr)

                  ! as root, write received reactions
                  do i = 1,ndofs
                     indv = isvgvn(i)

                     rea(indv) = reas(i)
                  end do

                  deallocate(reas)
                  deallocate(isvgvn)
                  ! add counter
                  acceptedp(iproc+1) = acceptedp(iproc+1) + 1
               end if
            end do
            if (all(acceptedp.eq.levels(ilevel)%sub2proc(2:nproc+1) - levels(ilevel)%sub2proc(1:nproc))) then
               exit
            end if
         end do
         deallocate(acceptedp)
      end if

end subroutine

!*****************************************************************************************
subroutine levels_dd_dotprod_subdomain_local(ilevel,isub, vec1,lvec1, vec2,lvec2, dotprod)
!*****************************************************************************************
! Subroutine used for indirect access to DD data
! multiplies vector 1 and vector 2 using weights in DD at interface, but uses
! subdomain arrays (i.e. including intrior )
      use module_utils
      implicit none

      ! length of vector
      integer,intent(in) ::   ilevel 
      integer,intent(in) ::   isub
      ! vectors to multiply
      integer,intent(in) ::  lvec1
      real(kr), intent(in) :: vec1(lvec1)
      integer,intent(in) ::  lvec2
      real(kr), intent(in) :: vec2(lvec2)
      
      ! result
      real(kr), intent(out) :: dotprod

      ! local vars
      character(*),parameter:: routine_name = 'LEVELS_DD_DOTPROD_SUBDOMAIN_LOCAL'
      integer :: isub_loc

      ! local index of subdomain
      call get_index(isub+levels_numshift, levels(ilevel)%indexsub,levels(ilevel)%lindexsub, isub_loc)
      if ( isub_loc .eq. -1 ) then
         call error( routine_name, 'Index of subdomain not found among local subdomains: ', isub )
      end if

      ! add data from module and call function from adaptive module
      call dd_dotprod_subdomain_local(levels(ilevel)%subdomains(isub_loc), vec1,lvec1, vec2,lvec2, dotprod)

end subroutine

!***********************************************************************************
subroutine levels_dd_dotprod_local(ilevel,isub_loc, vec1,lvec1, vec2,lvec2, dotprod)
!***********************************************************************************
! Subroutine used for indirect access to DD data
! multiplies vector 1 and vector 2 using weights in DD, only uses vectors at
! interface
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

!***********************************************************************
subroutine levels_dd_get_interface_size(ilevel,isub_loc, ndofis, nnodis)
!***********************************************************************
! Subroutine used for indirect access to DD data
! obtain interface size of subdomain isub_loc at level ilevel
      use module_utils
      implicit none

      ! length of vector
      integer,intent(in) ::   ilevel 
      integer,intent(in) ::   isub_loc 
      
      ! size of interface
      integer, intent(out) :: ndofis
      integer, intent(out) :: nnodis

      ! local vars
      character(*),parameter:: routine_name = 'LEVELS_DD_GET_INTERFACE_SIZE'

      ! check prerequisites
      if (.not.levels(ilevel)%is_level_prepared) then
         call error(routine_name,'Level is not prepared:',ilevel)
      end if

      ! call DD module function with data from module
      call dd_get_interface_size(levels(ilevel)%subdomains(isub_loc),ndofis,nnodis)

end subroutine

!*********************************************************************
subroutine levels_dd_fix_bc_interface_dual(ilevel,isub_loc,resi,lresi)
!*********************************************************************
! Subroutine used for indirect access to DD data
! module for distributed Krylov data storage
! module with utility routines
      use module_utils

      implicit none

      integer,intent(in)     :: ilevel
      integer,intent(in)     :: isub_loc
      integer,intent(in)     :: lresi
      real(kr),intent(inout) ::  resi(lresi)

! local vars
      character(*),parameter:: routine_name = 'LEVELS_DD_FIX_BC_INTERFACE_DUAL'

      ! check prerequisites
      if (.not.levels(ilevel)%is_level_prepared) then
         call error(routine_name,'Level is not prepared:',ilevel)
      end if

      ! fix BC on subdomain
      call dd_fix_bc_interface_dual(levels(ilevel)%subdomains(isub_loc), resi,lresi)

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
      if (allocated(level%iets)) then
         deallocate (level%iets)
      end if
      level%liets = 0
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

      level%nelem   = 0
      level%nnod    = 0
      level%ndof    = 0
      level%ndim    = 0
      level%meshdim = 0

      level%nsub  = 0

      level%ncorner  = 0
      level%nedge    = 0
      level%nface    = 0
      level%nnodc    = 0

      level%linet = 0
      level%lnnet = 0
      level%lnndf = 0
      level%lxyz1 = 0
      level%lxyz2 = 0
      level%lsol = 0

end subroutine

!*************************
subroutine levels_finalize
!*************************
! Subroutine for initialization of levels data
      use module_mumps
      implicit none
      include "mpif.h"

! local variables
      integer :: ilevel, start(1), irl
      integer :: ierr
      
! destroy MUMPS structure of the last level
      if (is_mumps_coarse_ready) then
         call mumps_finalize(mumps_coarse)
         is_mumps_coarse_ready = .false.
      end if

! deallocate basic structure
      if (allocated(levels)) then
         ! clear data of particular subdomains
         start = lbound(levels)
         do irl = start(1),llevels
            ! reverse levels
            ilevel = llevels + start(1) - irl

            call levels_clear_level(levels(ilevel))

            if (levels(ilevel)%is_new_comm_created .and. levels(ilevel)%i_am_active_in_this_level) then
               call MPI_COMM_FREE(levels(ilevel)%comm_all,ierr)
            end if
         end do

         deallocate (levels)
      end if
      llevels = 0
      nlevels = 0
      iactive_level = 0


end subroutine

!*******************************
subroutine levels_set_profile_on
!*******************************
! auxiliary routine to switch profiling flag on
profile = .true.
end subroutine

!********************************
subroutine levels_set_profile_off
!********************************
! auxiliary routine to switch profiling flag off
profile = .false.
end subroutine

end module module_levels

