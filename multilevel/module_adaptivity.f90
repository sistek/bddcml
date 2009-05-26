module module_adaptivity
!***********************
! Module for adaptive search of constraints for BDDC preconditioner
! Jakub Sistek, Denver, 3/2009

! type of reals
integer,parameter,private  :: kr = kind(1.D0)
! numerical zero
real(kr),parameter,private :: numerical_zero = 1.e-12_kr

! treshold on eigenvalues to define an adaptive constraint
real(kr),parameter,private :: treshold_eigval = 1.5_kr

! debugging 
logical,parameter,private :: debug = .false.

! table of pairs of eigenproblems to compute
! structure:
!  PROC | IGLOB | ISUB | JSUB 
integer,private            :: lpair_subdomains1 
integer,parameter,private  :: lpair_subdomains2 = 4
integer,allocatable,private :: pair_subdomains(:,:)


integer,private :: comm_calls = 0

integer,private :: comm_myid

integer,private :: neigvec, problemsize

integer,private :: ndofi_i, ndofi_j
integer,private :: nnodci, nnodcj
integer,private :: nnodi_i,nnodi_j

integer,private :: comm_myplace1, comm_myplace2 
integer,private :: comm_myisub, comm_myjsub 
integer,private :: comm_mygglob

integer,private ::             lbufsend_i,    lbufsend_j
real(kr),allocatable,private :: bufsend_i(:),  bufsend_j(:)
integer,private ::             lbufsend   , lbufrecv
real(kr),allocatable,private :: bufsend(:),  bufrecv(:)
integer,private ::            lkbufsend   
integer,allocatable,private :: kbufsend(:)
integer,private ::            llbufa   
integer,allocatable,private :: lbufa(:)

! array for serving to eigensolvers
integer,private ::            ninstructions 
integer,private ::            linstructions1
integer,parameter,private ::  linstructions2 = 4
integer,allocatable,private :: instructions(:,:)

! weigth matrix in interesting variables
integer,private ::             lweight
real(kr),allocatable,private :: weight(:)

! dependence of interface variables
integer,private ::            lpairslavery
integer,allocatable,private :: pairslavery(:)

! matrix of local constraints D_ij
integer,private ::             lddij ! leading dimension
integer,private ::              ldij1,ldij2
real(kr),allocatable,private ::  dij(:,:)

! MPI related arrays and variables
integer,private ::            lrequest
integer,allocatable,private :: request(:)
integer,private             :: lstatarray1
integer,private             :: lstatarray2
integer,allocatable,private :: statarray(:,:)
integer,private             :: comm_comm

! LAPACK QR related variables
integer,private ::             ltau
real(kr),allocatable,private :: tau(:)
integer,private ::             lwork
real(kr),allocatable,private :: work(:)

contains

!*************************************************
subroutine adaptivity_init(myid,comm,idpair,npair)
!*************************************************
! Subroutine for initialization of adaptive search of constraints
      use module_utils
      implicit none
      include "mpif.h"

! number of processor
      integer,intent(in) :: myid
! communicator
      integer,intent(in) :: comm
! file unit with opened list of pairs
      integer,intent(in) :: idpair

! local variables
      integer :: npair
      integer :: ierr
      integer :: ldata

      integer :: ipair, j

! read pairs to be computed data
      if (myid.eq.0) then
         read(idpair,*) npair
      end if
!*****************************************************************MPI
      call MPI_BCAST(npair,1, MPI_INTEGER, 0, comm, ierr)
!*****************************************************************MPI
      lpair_subdomains1 = npair
      allocate(pair_subdomains(lpair_subdomains1,lpair_subdomains2))
      if (myid.eq.0) then
         do ipair = 1,npair
            ! first column is associated with processors - initialize it to -1 = no processor assigned
            pair_subdomains(ipair,1) = -1
            read(idpair,*) (pair_subdomains(ipair,j), j = 2,lpair_subdomains2)
         end do
      end if
      close(idpair)
      ldata = lpair_subdomains1*lpair_subdomains2
!*****************************************************************MPI
      call MPI_BCAST(pair_subdomains,ldata, MPI_INTEGER, 0, comm, ierr)
!*****************************************************************MPI
 
      ! print what is loaded
      if (debug) then
         call adaptivity_print_pairs(myid)
      end if

      return
end subroutine

!*********************************************************
subroutine adaptivity_assign_pairs(npair,nproc,npair_locx)
!*********************************************************
! Subroutine for distribution of load to processors
      implicit none

! Global number of pairs to compute eigenproblems
      integer,intent(in) :: npair
! Number of processors
      integer,intent(in) :: nproc

! local variables
      integer :: iproc, ipair_loc, ipair, npair_locx

      npair_locx = (npair + nproc - 1)/nproc
      ipair = 0
      do iproc = 0,nproc-1
         do ipair_loc = 1,npair_locx
            ipair = ipair + 1
            if (ipair.le.npair) then
               pair_subdomains(ipair,1) = iproc
            end if
         end do
      end do
end subroutine

!******************************************************************************************************************************
subroutine adaptivity_get_active_pairs(iround,nproc,npair,npair_locx,active_pairs,lactive_pairs,nactive_pairs,all_pairs_solved)
!******************************************************************************************************************************
! Subroutine for activating and deactivating pairs
      implicit none

! number of round 
      integer,intent(in) :: iround
! number of pairs
      integer,intent(in) :: npair
! maximal local number of eigenproblems at one processor
      integer,intent(in) :: npair_locx
! number of processors
      integer,intent(in) :: nproc
! indices of active pairs
      integer,intent(in) :: lactive_pairs
      integer,intent(out) :: active_pairs(lactive_pairs)
! number of active pairs
      integer, intent(out) :: nactive_pairs
! set this to true if all pairs are solved
      logical, intent(out) :: all_pairs_solved

! local variables
      integer :: ipair, indpair, iactive_pair, i

      if (iround.gt.npair_locx) then
         all_pairs_solved = .true.
         return
      else
         all_pairs_solved = .false.
      end if

      indpair      = iround
      ipair        = 0
      iactive_pair = 0
      do i = 1,nproc
         ipair = ipair + 1
         if (indpair.le.npair) then
            iactive_pair = iactive_pair + 1
            active_pairs(ipair) = indpair
         else
            active_pairs(ipair) = 0
         end if

         indpair = indpair + npair_locx
      end do
      nactive_pairs = iactive_pair
end subroutine

!**********************************************************************
subroutine adaptivity_get_my_pair(iround,myid,npair_locx,npair,my_pair)
!**********************************************************************
! Subroutine for getting number of pair to solve
      implicit none

! number of round 
      integer,intent(in) :: iround
! processor ID
      integer,intent(in) :: myid
! Maximal local number of eigenproblems at one processor
      integer,intent(in) :: npair_locx
! Global number of pairs to compute eigenproblems
      integer,intent(in) :: npair
! number of pair to solve
      integer,intent(out) :: my_pair 

      my_pair = (myid*npair_locx) + iround
      if (my_pair.gt.npair) then
         my_pair = -1
      end if
end subroutine

!*************************************************************************
subroutine adaptivity_solve_eigenvectors(myid,comm,npair_locx,npair,nproc)
!*************************************************************************
! Subroutine for parallel solution of distributed eigenproblems
      use module_dd
      use module_utils
      implicit none
      include "mpif.h"

! number of processor
      integer,intent(in) :: myid
! communicator
      integer,intent(in) :: comm

! Maximal local number of eigenproblems at one processor
      integer,intent(in) :: npair_locx

! Global number of pairs to compute eigenproblems
      integer,intent(in) :: npair

! Number of processors
      integer,intent(in) :: nproc

! Maximal number of eigenvectors per problem
      integer,parameter :: neigvecx = 20 

! local variables
      integer :: isub, jsub, ipair, iactive_pair, iround
      integer :: gglob, my_pair, &
                 nactive_pairs, owner,&
                 place1, place2, pointbuf, i, j, iinstr, indcorner, icommon,&
                 indc_i, indc_j, indi_i, indi_j, ndofn, nconstr, iconstr, inodi,&
                 pointv_i, pointv_j, shift, indcommon, ndofcomm, idofn, irhoicomm,&
                 point_i, point_j, indiv, nadaptive, ioper, nadaptive_rcv, ind

      integer :: ndofi
      integer :: nnodc
      integer :: nnodi


      integer ::            lpair_data
      integer,allocatable :: pair_data(:)

      ! numbers of active pairs
      integer ::            lactive_pairs
      integer,allocatable :: active_pairs(:)

      logical :: all_pairs_solved

      ! eigenvectors and eigenvalues
      integer ::             leigvec
      real(kr),allocatable :: eigvec(:)
      integer ::             leigval
      real(kr),allocatable :: eigval(:)

      ! corner information
      integer ::            ncommon_corners
      integer ::            lglobal_corner_number_i,   lglobal_corner_number_j
      integer,allocatable :: global_corner_number_i(:), global_corner_number_j(:)
      integer ::            lglobal_corner_number
      integer,allocatable :: global_corner_number(:)
      integer ::            lcommon_corners
      integer,allocatable :: common_corners(:)
      integer ::            licnsin_i,   licnsin_j
      integer,allocatable :: icnsin_i(:), icnsin_j(:)
      integer ::            licnsin
      integer,allocatable :: icnsin(:)
      integer ::            lnndfi_i,   lnndfi_j
      integer,allocatable :: nndfi_i(:), nndfi_j(:)
      integer ::            lnndfi
      integer,allocatable :: nndfi(:)
      integer ::            lkdofi_i,   lkdofi_j
      integer,allocatable :: kdofi_i(:), kdofi_j(:)
      integer ::             lrhoi_i,   lrhoi_j
      real(kr),allocatable :: rhoi_i(:), rhoi_j(:)
      integer ::             lrhoi
      real(kr),allocatable :: rhoi(:)
      integer ::            liingn_i,   liingn_j
      integer,allocatable :: iingn_i(:), iingn_j(:)
      integer ::            liingn
      integer,allocatable :: iingn(:)
      integer ::             lrhoicomm
      real(kr),allocatable :: rhoicomm(:)
      integer ::             lconstraints1, lconstraints2
      real(kr),allocatable :: constraints(:,:)

      integer ::             lcadapt1, lcadapt2
      real(kr),allocatable :: cadapt(:,:)

      integer ::            ncommon_interface
      integer ::            lcommon_interface
      integer,allocatable :: common_interface(:)

      real(kr),external :: ddot

      ! LAPACK QR related variables
      integer :: lapack_info

      ! LOBPCG related variables
      integer ::  lobpcg_maxit, lobpcg_verbosity, use_vec_values, lobpcg_iter
      real(kr) :: lobpcg_tol

      ! MPI related variables
      integer :: ierr, ireq, nreq
      
      ! allocate table for work instructions - the worst case is that in each
      ! round, I have to compute all the subdomains, i.e. 2 for each pair
      linstructions1 = 2*nproc
      allocate(instructions(linstructions1,linstructions2))
      ! prepare MPI arrays 
      lrequest = 2*nproc + 2
      allocate(request(lrequest))
      lstatarray1 = MPI_STATUS_SIZE
      lstatarray2 = lrequest
      allocate(statarray(lstatarray1,lstatarray2))

      ! loop over number of rounds for solution of eigenproblems
      lactive_pairs = nproc
      allocate(active_pairs(lactive_pairs))

      ! prepare space for pair_data
      lpair_data = lpair_subdomains2
      allocate(pair_data(lpair_data))

      iround = 0
      ! Loop over rounds of eigenvalue solves
      do 
         iround = iround + 1

         ! each round of eigenproblems has its structure - determine active pairs
         call adaptivity_get_active_pairs(iround,nproc,npair,npair_locx,active_pairs,lactive_pairs,nactive_pairs,all_pairs_solved)
         if (all_pairs_solved) then
            exit
         end if

         ! determine which pair I compute
         call adaptivity_get_my_pair(iround,myid,npair_locx,npair,my_pair)

         if (my_pair.ge.0) then

            call adaptivity_get_pair_data(my_pair,pair_data,lpair_data)
   
            comm_mygglob = pair_data(2)
            comm_myisub  = pair_data(3)
            comm_myjsub  = pair_data(4)

            ! where are these subdomains ?
            call dd_where_is_subdomain(comm_myisub,comm_myplace1)
            call dd_where_is_subdomain(comm_myjsub,comm_myplace2)
         end if
         !write(90+myid,*) 'myid =',myid, 'pair_data:'
         !write(90+myid,*) pair_data


         ! determine working instructions for sending subdomain matrices
         ! go through pairs that are active in this round
         instructions  = 0
         ninstructions = 0
         pointbuf      = 1
         do ipair = 1,nactive_pairs
            iactive_pair = active_pairs(ipair)
            
            call adaptivity_get_pair_data(iactive_pair,pair_data,lpair_data)
            owner  = pair_data(1)
            gglob  = pair_data(2)
            isub   = pair_data(3)
            jsub   = pair_data(4)

            ! where are these subdomains ?
            call dd_where_is_subdomain(isub,place1)
            call dd_where_is_subdomain(jsub,place2)
            !write(90+myid,*) 'myid =',myid, 'place1:',place1,'place2:',place2


            if (myid.eq.place1) then
               ! add instruction
               ninstructions = ninstructions + 1

               ! who I will send data to
               instructions(ninstructions,1) = owner
               ! subdomain number
               instructions(ninstructions,2) = isub
               ! global glob number
               instructions(ninstructions,3) = gglob
            end if

            if (myid.eq.place2) then
               ! add instruction
               ninstructions = ninstructions + 1

               ! who I will send data to
               instructions(ninstructions,1) = owner
               ! subdomain number
               instructions(ninstructions,2) = jsub
               ! global glob number
               instructions(ninstructions,3) = gglob
            end if
         end do

         ! the scheme for communication is ready
         !write(90+myid,*) 'myid =',myid, 'instructions:'
         !do i = 1,ninstructions
         !   write(90+myid,*) instructions(i,:)
         !end do
         !call flush(90+myid)

         ! build the local matrix of projection on common globs for active pair
         !  get sizes of interface of subdomains in my problem
         ireq = 0
         if (my_pair.ge.0) then
            ! receive sizes of interfaces of subdomains involved in my problem

            ireq = ireq + 1
            call MPI_IRECV(ndofi_i,1,MPI_INTEGER,comm_myplace1,comm_myisub,comm,request(ireq),ierr)

            ireq = ireq + 1
            call MPI_IRECV(ndofi_j,1,MPI_INTEGER,comm_myplace2,comm_myjsub,comm,request(ireq),ierr)
         end if
         ! send sizes of subdomains involved in problems
         do iinstr = 1,ninstructions
            owner = instructions(iinstr,1)
            isub  = instructions(iinstr,2)
            call dd_get_interface_size(myid,isub,ndofi,nnodi)

            ireq = ireq + 1
            call MPI_ISEND(ndofi,1,MPI_INTEGER,owner,isub,comm,request(ireq),ierr)
         end do
         nreq = ireq
         if (nreq.gt.size(request)) then
            write(*,*) 'ADAPTIVITY_SOLVE_EIGENVECTORS: Not enough space for MPI requests.'
            call error_exit
         end if
         call MPI_WAITALL(nreq, request, statarray, ierr)
         if (debug) then
            write(*,*) 'I am ',myid, 'All messages in pack 1 received, MPI is fun!.'
            call flush(6)
         end if

         !if (my_pair.ge.0) then
         !   write(90+myid,*) 'ndofi_i',ndofi_i,'ndofi_j',ndofi_j
         !end if

         !  get number of corners to find common subset
         ireq = 0
         if (my_pair.ge.0) then
            ! receive numbers of corners

            ireq = ireq + 1
            call MPI_IRECV(nnodci,1,MPI_INTEGER,comm_myplace1,comm_myisub,comm,request(ireq),ierr)

            ireq = ireq + 1
            call MPI_IRECV(nnodcj,1,MPI_INTEGER,comm_myplace2,comm_myjsub,comm,request(ireq),ierr)
         end if
         ! send sizes of subdomains involved in problems
         do iinstr = 1,ninstructions
            owner = instructions(iinstr,1)
            isub  = instructions(iinstr,2)
            call dd_get_number_of_corners(myid,isub,nnodc)

            ireq = ireq + 1
            call MPI_ISEND(nnodc,1,MPI_INTEGER,owner,isub,comm,request(ireq),ierr)
         end do
         nreq = ireq
         call MPI_WAITALL(nreq, request, statarray, ierr)
         if (debug) then
            write(*,*) 'I am ',myid, 'All messages in pack 2 received, MPI is fun!.'
            call flush(6)
         end if

         !if (my_pair.ge.0) then
         !   write(90+myid,*) 'nnodci',nnodci,'nnodcj',nnodcj
         !end if

         !  get number of nodes on interface
         ireq = 0
         if (my_pair.ge.0) then
            ! receive numbers of corners

            ireq = ireq + 1
            call MPI_IRECV(nnodi_i,1,MPI_INTEGER,comm_myplace1,comm_myisub,comm,request(ireq),ierr)

            ireq = ireq + 1
            call MPI_IRECV(nnodi_j,1,MPI_INTEGER,comm_myplace2,comm_myjsub,comm,request(ireq),ierr)
         end if
         ! send sizes of subdomains involved in problems
         do iinstr = 1,ninstructions
            owner = instructions(iinstr,1)
            isub  = instructions(iinstr,2)
            call dd_get_interface_size(myid,isub,ndofi,nnodi)

            ireq = ireq + 1
            call MPI_ISEND(nnodi,1,MPI_INTEGER,owner,isub,comm,request(ireq),ierr)
         end do
         nreq = ireq
         call MPI_WAITALL(nreq, request, statarray, ierr)
         if (debug) then
            write(*,*) 'I am ',myid, 'All messages in pack 2.5 received, MPI is fun!.'
            call flush(6)
         end if
         !if (my_pair.ge.0) then
         !   write(90+myid,*) 'nnodi_i',nnodi_i,'nnodi_j',nnodi_j
         !end if

         ! allocate space for corners
         if (my_pair.ge.0) then
            problemsize = ndofi_i + ndofi_j
            lglobal_corner_number_i = nnodci
            lglobal_corner_number_j = nnodcj
            allocate(global_corner_number_i(lglobal_corner_number_i),global_corner_number_j(lglobal_corner_number_j))
            licnsin_i = nnodci
            licnsin_j = nnodcj
            allocate(icnsin_i(licnsin_i),icnsin_j(licnsin_j))
            lnndfi_i = nnodi_i
            lnndfi_j = nnodi_j
            allocate(nndfi_i(lnndfi_i),nndfi_j(lnndfi_j))
            lrhoi_i = ndofi_i
            lrhoi_j = ndofi_j
            allocate(rhoi_i(lrhoi_i),rhoi_j(lrhoi_j))
            liingn_i = nnodi_i
            liingn_j = nnodi_j
            allocate(iingn_i(liingn_i),iingn_j(liingn_j))
         end if
         ! get data about corners
         ireq = 0
         if (my_pair.ge.0) then
            ! receive global numbers of corners

            ireq = ireq + 1
            call MPI_IRECV(global_corner_number_i,nnodci,MPI_INTEGER,comm_myplace1,comm_myisub,comm,request(ireq),ierr)

            ireq = ireq + 1
            call MPI_IRECV(global_corner_number_j,nnodcj,MPI_INTEGER,comm_myplace2,comm_myjsub,comm,request(ireq),ierr)
         end if
         ! send sizes of subdomains involved in problems
         do iinstr = 1,ninstructions
            owner = instructions(iinstr,1)
            isub  = instructions(iinstr,2)
            call dd_get_number_of_corners(myid,isub,nnodc)

            lglobal_corner_number = nnodc
            allocate(global_corner_number(lglobal_corner_number))
            call dd_get_subdomain_corners(myid,isub,global_corner_number,lglobal_corner_number)

            ireq = ireq + 1
            call MPI_ISEND(global_corner_number,nnodc,MPI_INTEGER,owner,isub,comm,request(ireq),ierr)
            deallocate(global_corner_number)
         end do
         nreq = ireq
         call MPI_WAITALL(nreq, request, statarray, ierr)
         if (debug) then
            write(*,*) 'I am ',myid, 'All messages in pack 3 received, MPI is fun!.'
            call flush(6)
         end if
         !if (my_pair.ge.0) then
         !   write(90+myid,*) 'global_corner_number_i'
         !   write(90+myid,*)  global_corner_number_i
         !   write(90+myid,*) 'global_corner_number_j'
         !   write(90+myid,*)  global_corner_number_j
         !end if
         ! get data about corners
         ireq = 0
         if (my_pair.ge.0) then
            ! receive global numbers of corners

            ireq = ireq + 1
            call MPI_IRECV(icnsin_i,nnodci,MPI_INTEGER,comm_myplace1,comm_myisub,comm,request(ireq),ierr)

            ireq = ireq + 1
            call MPI_IRECV(icnsin_j,nnodcj,MPI_INTEGER,comm_myplace2,comm_myjsub,comm,request(ireq),ierr)
         end if
         ! send sizes of subdomains involved in problems
         do iinstr = 1,ninstructions
            owner = instructions(iinstr,1)
            isub  = instructions(iinstr,2)
            call dd_get_number_of_corners(myid,isub,nnodc)

            licnsin = nnodc
            allocate(icnsin(licnsin))
            call dd_get_subdomain_corners_map(myid,isub,icnsin,licnsin)

            ireq = ireq + 1
            call MPI_ISEND(icnsin,nnodc,MPI_INTEGER,owner,isub,comm,request(ireq),ierr)
            deallocate(icnsin)
         end do
         nreq = ireq
         call MPI_WAITALL(nreq, request, statarray, ierr)
         if (debug) then
            write(*,*) 'I am ',myid, 'All messages in pack 4 received, MPI is fun!.'
            call flush(6)
         end if

         ! get numbers of dof on interface
         ireq = 0
         if (my_pair.ge.0) then
            ! receive global numbers of corners

            ireq = ireq + 1
            call MPI_IRECV(nndfi_i,nnodi_i,MPI_INTEGER,comm_myplace1,comm_myisub,comm,request(ireq),ierr)

            ireq = ireq + 1
            call MPI_IRECV(nndfi_j,nnodi_j,MPI_INTEGER,comm_myplace2,comm_myjsub,comm,request(ireq),ierr)
         end if
         ! send sizes of subdomains involved in problems
         do iinstr = 1,ninstructions
            owner = instructions(iinstr,1)
            isub  = instructions(iinstr,2)
            call dd_get_interface_size(myid,isub,ndofi,nnodi)

            lnndfi = nnodi
            allocate(nndfi(lnndfi))
            call dd_get_subdomain_interface_nndf(myid,isub,nndfi,lnndfi)

            ireq = ireq + 1
            call MPI_ISEND(nndfi,nnodi,MPI_INTEGER,owner,isub,comm,request(ireq),ierr)
            deallocate(nndfi)
         end do
         nreq = ireq
         call MPI_WAITALL(nreq, request, statarray, ierr)
         if (debug) then
            write(*,*) 'I am ',myid, 'All messages in pack 5 received, MPI is fun!.'
            call flush(6)
         end if
         !if (my_pair.ge.0) then
         !   write(90+myid,*) 'nndfi_i'
         !   write(90+myid,*)  nndfi_i
         !   write(90+myid,*) 'nndfi_j'
         !   write(90+myid,*)  nndfi_j
         !   write(90+myid,*) 'icnsin_i'
         !   write(90+myid,*)  icnsin_i
         !   write(90+myid,*) 'icnsin_j'
         !   write(90+myid,*)  icnsin_j
         !end if

         ! find intersection of corners
         if (my_pair.ge.0) then
            lcommon_corners = min(lglobal_corner_number_i,lglobal_corner_number_j)
            allocate(common_corners(lcommon_corners))
            call get_array_intersection(global_corner_number_i,lglobal_corner_number_i,&
                                        global_corner_number_j,lglobal_corner_number_j,&
                                        common_corners,lcommon_corners,ncommon_corners)
         end if

         ! build matrix D_ij^T for the pair of subdomains
         if (my_pair.ge.0) then
            ! first determine the size of D_ij
            nconstr = 0
            do icommon = 1,ncommon_corners
               indcorner = common_corners(icommon)
               ! find index of this corner
               call get_index(indcorner,global_corner_number_i,lglobal_corner_number_i,indc_i) 
               if (indc_i.gt.0 ) then
                  indi_i = icnsin_i(indc_i)
                  ndofn  = nndfi_i(indi_i)

                  nconstr = nconstr + ndofn
               else
                  call get_index(indcorner,global_corner_number_j,lglobal_corner_number_j,indc_j) 
                  indi_j = icnsin_j(indc_j)
                  ndofn  = nndfi_j(indi_j)

                  nconstr = nconstr + ndofn
               end if
            end do
            ! prepare space for dense matrix D_ij
            ldij1 = ndofi_i + ndofi_j
            ldij2 = nconstr
            allocate(dij(ldij1,ldij2))
            call zero(dij,ldij1,ldij2)

            ! prepare arrays kdofi_i and kdofi_j
            lkdofi_i = nnodi_i
            lkdofi_j = nnodi_j
            allocate(kdofi_i(lkdofi_i),kdofi_j(lkdofi_j))
            kdofi_i(1) = 0
            do inodi = 2,nnodi_i
               kdofi_i(inodi) = kdofi_i(inodi-1) + nndfi_i(inodi-1)
            end do
            kdofi_j(1) = 0
            do inodi = 2,nnodi_j
               kdofi_j(inodi) = kdofi_j(inodi-1) + nndfi_j(inodi-1)
            end do

            iconstr = 0
            do icommon = 1,ncommon_corners
               indcorner = common_corners(icommon)
               ! find index of this corner
               call get_index(indcorner,global_corner_number_i,lglobal_corner_number_i,indc_i) 
               if (indc_i.gt.0) then
                  indi_i = icnsin_i(indc_i)
                  ndofn  = nndfi_i(indi_i)
               end if
               call get_index(indcorner,global_corner_number_j,lglobal_corner_number_j,indc_j) 
               if (indc_j.gt.0) then
                  indi_j = icnsin_j(indc_j)
                  ndofn  = nndfi_i(indi_i)
               end if

               shift = ndofi_i
               if (indc_i.gt.0) then
                  pointv_i = kdofi_i(indi_i)
               end if
               if (indc_j.gt.0) then
                  pointv_j = kdofi_j(indi_j)
               end if
               do i = 1,ndofn
                  iconstr = iconstr + 1

                  ! contribution of subdomain i
                  if (indc_i.gt.0) then
                     dij(pointv_i+i,iconstr)       = 1._kr
                  end if
                  ! contribution of subdomain j
                  if (indc_j.gt.0) then
                     dij(shift+pointv_j+i,iconstr) = -1._kr
                  end if
               end do
            end do
         end if

!         if (my_pair.ge.0) then
!            do i = 1,ldij1
!               print *,(dij(i,j),j = 1,ldij2)
!            end do
!         end if

         ! Prepare projection onto null D_ij
         if (my_pair.ge.0) then
            ! QR decomposition of matrix Dij^T
            ! LAPACK arrays
            ltau = ldij2
            allocate(tau(ltau))
            lwork = ldij2
            allocate(work(lwork))
            ! leading dimension
            lddij = max(1,ldij1)
            call DGEQRF( ldij1, ldij2, dij, lddij, tau, work, lwork, lapack_info)
            if (lapack_info.ne.0) then
               write(*,*) 'ADAPTIVITY_SOLVE_EIGENVECTORS: Error in LAPACK QR factorization of matrix D_ij^T: ', lapack_info
               call error_exit
            end if
            ! in space of D_ij are now stored factors R and Householder reflectors v
         end if

         ! prepare operator (I-R_ijE_ij)
         ! get data about interface weigths
         ireq = 0
         if (my_pair.ge.0) then
            ! receive diagonal entries

            ireq = ireq + 1
            call MPI_IRECV(rhoi_i,ndofi_i,MPI_DOUBLE_PRECISION,comm_myplace1,comm_myisub,comm,request(ireq),ierr)

            ireq = ireq + 1
            call MPI_IRECV(rhoi_j,ndofi_j,MPI_DOUBLE_PRECISION,comm_myplace2,comm_myjsub,comm,request(ireq),ierr)
         end if
         do iinstr = 1,ninstructions
            owner = instructions(iinstr,1)
            isub  = instructions(iinstr,2)
            call dd_get_interface_size(myid,isub,ndofi,nnodi)

            lrhoi = ndofi
            allocate(rhoi(lrhoi))
            call dd_get_interface_diagonal(myid,isub, rhoi,lrhoi)

            ireq = ireq + 1
            call MPI_ISEND(rhoi,ndofi,MPI_DOUBLE_PRECISION,owner,isub,comm,request(ireq),ierr)
            deallocate(rhoi)
         end do
         nreq = ireq
         call MPI_WAITALL(nreq, request, statarray, ierr)
         if (debug) then
            write(*,*) 'I am ',myid, 'All messages in pack 7 received, MPI is fun!.'
            call flush(6)
         end if
         !if (my_pair.ge.0) then
         !   write(90+myid,*)'rhoi_i'
         !   write(90+myid,*) rhoi_i
         !   write(90+myid,*)'rhoi_j'
         !   write(90+myid,*) rhoi_j
         !end if

         ! get data about common interface
         ireq = 0
         if (my_pair.ge.0) then
            ! receive mapping of interface nodes into global nodes

            ireq = ireq + 1
            call MPI_IRECV(iingn_i,nnodi_i,MPI_INTEGER,comm_myplace1,comm_myisub,comm,request(ireq),ierr)

            ireq = ireq + 1
            call MPI_IRECV(iingn_j,nnodi_j,MPI_INTEGER,comm_myplace2,comm_myjsub,comm,request(ireq),ierr)
         end if
         ! send sizes of subdomains involved in problems
         do iinstr = 1,ninstructions
            owner = instructions(iinstr,1)
            isub  = instructions(iinstr,2)
            call dd_get_interface_size(myid,isub,ndofi,nnodi)

            liingn = nnodi
            allocate(iingn(liingn))
            call dd_get_interface_global_numbers(myid,isub,iingn,liingn)

            ireq = ireq + 1
            call MPI_ISEND(iingn,nnodi,MPI_INTEGER,owner,isub,comm,request(ireq),ierr)
            deallocate(iingn)
         end do
         nreq = ireq
         call MPI_WAITALL(nreq, request, statarray, ierr)
         if (debug) then
            write(*,*)  'I am ',myid,'All messages in pack 8 received, MPI is fun!.'
            call flush(6)
         end if
         !if (my_pair.ge.0) then
         !   write(90+myid,*)'iingn_i'
         !   write(90+myid,*) iingn_i
         !   write(90+myid,*)'iingn_j'
         !   write(90+myid,*) iingn_j
         !end if

         ! find common intersection of interface nodes
         if (my_pair.ge.0) then
            lcommon_interface = min(nnodi_i,nnodi_j)
            allocate(common_interface(lcommon_interface))
            if (my_pair.ge.0) then
               call get_array_intersection(iingn_i,liingn_i, iingn_j,liingn_j,&
                                           common_interface,lcommon_interface,ncommon_interface)
            end if
            ! determine size of common interface
            ndofcomm = 0
            do icommon = 1,ncommon_interface
               indcommon = common_interface(icommon)
               call get_index(indcommon,iingn_i,liingn_i,indi_i)
               ndofn = nndfi_i(indi_i)

               ndofcomm = ndofcomm + ndofn
            end do
            ! prepare space for sum of interface entries
            lrhoicomm = ndofcomm
            allocate(rhoicomm(lrhoicomm))
            irhoicomm = 0
            ! sum diagonal entries on interface
            do icommon = 1,ncommon_interface
               indcommon = common_interface(icommon)
               call get_index(indcommon,iingn_i,liingn_i,indi_i)
               point_i = kdofi_i(indi_i)
               call get_index(indcommon,iingn_j,liingn_j,indi_j)
               point_j = kdofi_j(indi_j)

               ndofn = nndfi_i(indi_i)

               do idofn = 1,ndofn
                  irhoicomm = irhoicomm + 1

                  rhoicomm(irhoicomm) = rhoi_i(point_i + idofn) + rhoi_j(point_j + idofn)
               end do
            end do
            ! prepare vector of mask and weigths representing (I - R_ij E_ij)
            lweight = ndofi_i + ndofi_j
            allocate(weight(lweight))
            do i = 1,problemsize
               weight(i) = 1._kr
            end do

            lpairslavery = ndofi_i + ndofi_j
            allocate(pairslavery(lpairslavery))
            call zero(pairslavery,lpairslavery)

            irhoicomm = 0
            shift = ndofi_i
            do icommon = 1,ncommon_interface
               indcommon = common_interface(icommon)
               call get_index(indcommon,iingn_i,liingn_i,indi_i)
               point_i = kdofi_i(indi_i)
               call get_index(indcommon,iingn_j,liingn_j,indi_j)
               point_j = kdofi_j(indi_j)

               ndofn = nndfi_i(indi_i)

               do idofn = 1,ndofn
                  irhoicomm = irhoicomm + 1

                  weight(        point_i + idofn) = rhoi_i(point_i + idofn) / rhoicomm(irhoicomm)
                  weight(shift + point_j + idofn) = rhoi_j(point_j + idofn) / rhoicomm(irhoicomm)

                  ! prepare slavery array - first of pair is the master, the second slave
                  indiv = point_i + idofn
                  pairslavery(        point_i + idofn) = indiv
                  pairslavery(shift + point_j + idofn) = indiv
               end do
            end do
            !write(90+myid,*) 'Pair slavery:', pairslavery

            deallocate(rhoicomm)
            !write(90+myid,*) 'weight',weight
         end if


         ! prepare space for eigenvectors
         ireq = 0
         if (my_pair.ge.0) then
            neigvec     = min(neigvecx,ndofcomm)
            leigvec = neigvec * problemsize
            leigval = neigvec
            allocate(eigvec(leigvec),eigval(leigval))

            ! prepare space for buffers
            lbufsend_i = ndofi_i
            lbufsend_j = ndofi_j
            allocate(bufsend_i(lbufsend_i),bufsend_j(lbufsend_j))

            ! distribute sizes of chunks of eigenvectors
            ireq = ireq + 1
            call MPI_ISEND(lbufsend_i,1,MPI_INTEGER,comm_myplace1,comm_myisub,comm,request(ireq),ierr)

            ireq = ireq + 1
            call MPI_ISEND(lbufsend_j,1,MPI_INTEGER,comm_myplace2,comm_myjsub,comm,request(ireq),ierr)
         end if
         ! prepare pointers to buffers with chunks of eigenvectors and their size
         llbufa = ninstructions
         allocate(lbufa(llbufa))
         do iinstr = 1,ninstructions
            owner = instructions(iinstr,1)
            isub  = instructions(iinstr,2)

            ireq = ireq + 1
            call MPI_IRECV(lbufa(iinstr),1,MPI_INTEGER,owner,isub,comm,request(ireq),ierr)
         end do
         nreq = ireq
         call MPI_WAITALL(nreq, request, statarray, ierr)
         ! prepare arrays kbufsend and 
         lkbufsend = ninstructions
         allocate(kbufsend(lkbufsend))
         if (lkbufsend.gt.0) then
            kbufsend(1)  = 1
            do i = 2,ninstructions
               kbufsend(i) = kbufsend(i-1) + lbufa(i-1)
            end do
         end if
         if (debug) then
            write(*,*)  'I am ',myid,'All messages in pack 9 received, MPI is fun!.'
            call flush(6)
         end if
         lbufsend = 0
         do i = 1,ninstructions
            lbufsend = lbufsend + lbufa(i)
         end do
         lbufrecv = lbufsend
         allocate(bufrecv(lbufrecv),bufsend(lbufsend))

         ! All arrays are ready for iterational process
         ! call LOBCPG
         comm_comm = comm
         comm_myid = myid
         if (my_pair.ge.0) then
            lobpcg_tol   = 1.e-5_kr
            lobpcg_maxit = 1000
            lobpcg_verbosity = 0
            use_vec_values = 0
            if (use_vec_values .eq. 1) then
              ! read initial guess of vectors
               write(*,*) 'neigvec =',neigvec,'problemsize =',problemsize
               open(unit = 99,file='start_guess.txt')
               do i = 1,problemsize
                  read(99,*) (eigvec((j-1)*problemsize + i),j = 1,neigvec)
               end do
               close(99)
            end if
               
            if (debug) then
               write(*,*) 'myid =',myid,', I am calling eigensolver for pair ',my_pair
               call flush(6)
            end if
            call lobpcg_driver(problemsize,neigvec,lobpcg_tol,lobpcg_maxit,lobpcg_verbosity,use_vec_values,&
                               eigval,eigvec,lobpcg_iter,ierr)
            if (ierr.ne.0) then
               write(*,'(a,i4,a,i6)') 'ADAPTIVITY_SOLVE_EIGENVECTORS: WARNING - LOBPCG exited with nonzero code ',ierr, &
                                      ' for pair',my_pair
            end if
            if (debug) then
               write(*,*) 'myid =',myid,', LOBPCG converged in ',lobpcg_iter,' iterations.'
               call flush(6)
            end if
            ! turn around the eigenvalues to be the largest
            eigval = -eigval
            write(*,*) 'eigval for pair:',my_pair,' between subdomains ',comm_myisub,' and',comm_myjsub
            write(*,'(f20.10)') eigval
            !write(90+myid,*) 'x ='
            !do i = 1,problemsize
            !   write(90+myid,'(30f13.7)') (eigvec((j-1)*problemsize + i),j = 1,neigvec)
            !end do
            !do i = 1,problemsize
            !   write(88, '(30e15.5)'), (eigvec((j-1)*problemsize + i),j = 1,neigvec)
            !end do
         end if
         call adaptivity_fake_lobpcg_driver

         ! select eigenvectors of eigenvalues exceeding treshold
         if (my_pair.ge.0) then

            nadaptive = count(eigval.ge.treshold_eigval)
            if (debug) then
               write(*,*) 'ADAPTIVITY_SOLVE_EIGENVECTORS: I am going to add ',nadaptive,' constraints for pair ',my_pair
               call flush(6)
            end if
            lconstraints1 = problemsize
            lconstraints2 = nadaptive
            allocate(constraints(lconstraints1,lconstraints2))
            ! construct constraints out of these eigenvectors
            ioper = 1 ! multiply by A
            do j = 1,nadaptive
               call adaptivity_mvecmult(problemsize,eigvec((j-1)*problemsize + 1),&
                    problemsize,constraints(1,j),problemsize,ioper)
            end do
         end if
         call adaptivity_fake_lobpcg_driver

         if (my_pair.ge.0) then
            !write(90+myid,*) 'Constraints to be added on pair ',my_pair
            !do i = 1,problemsize
            !   write(90+myid,'(100f15.3)') (constraints(i,j),j = 1,nadaptive)
            !end do
            !call flush(90+myid)
         end if

         ! REALLOCATE BUFFERS
         deallocate(bufrecv,bufsend)
         if (my_pair.ge.0) then
            deallocate(bufsend_i,bufsend_j)
         end if


         ! distribute adaptive constraints to slaves
         ! prepare space for these constraints
         ireq = 0
         if (my_pair.ge.0) then
            ! prepare space for buffers
            lbufsend_i = ndofi_i * nadaptive
            lbufsend_j = ndofi_j * nadaptive
            allocate(bufsend_i(lbufsend_i),bufsend_j(lbufsend_j))

            ! distribute sizes of chunks of eigenvectors
            ireq = ireq + 1
            call MPI_ISEND(lbufsend_i,1,MPI_INTEGER,comm_myplace1,comm_myisub,comm,request(ireq),ierr)

            ireq = ireq + 1
            call MPI_ISEND(lbufsend_j,1,MPI_INTEGER,comm_myplace2,comm_myjsub,comm,request(ireq),ierr)
         end if
         ! prepare pointers to buffers with chunks of eigenvectors and their size
         do iinstr = 1,ninstructions
            owner = instructions(iinstr,1)
            isub  = instructions(iinstr,2)

            ireq = ireq + 1
            call MPI_IRECV(lbufa(iinstr),1,MPI_INTEGER,owner,isub,comm,request(ireq),ierr)
         end do
         nreq = ireq
         call MPI_WAITALL(nreq, request, statarray, ierr)
         ! prepare arrays kbufsend and 
         if (lkbufsend .gt. 0) then
            kbufsend(1)  = 1
            do i = 2,ninstructions
               kbufsend(i) = kbufsend(i-1) + lbufa(i-1)
            end do
         end if
         if (debug) then
            write(*,*)  'myid =',myid,'All messages in pack 10 received, MPI is fun!.'
            call flush(6)
         end if
         lbufsend = 0
         do i = 1,ninstructions
            lbufsend = lbufsend + lbufa(i)
         end do
         lbufrecv = lbufsend
         allocate(bufrecv(lbufrecv),bufsend(lbufsend))

         ! distribute adaptive constraints to slaves
         ireq = 0
         if (my_pair.ge.0) then
            ind = 0
            do j = 1,nadaptive
               do i = 1,ndofi_i
                  ind = ind + 1
                  bufsend_i(ind) = constraints(i,j)
               end do
            end do

            shift = ndofi_i
            ind = 0
            do j = 1,nadaptive
               do i = 1,ndofi_j
                  ind = ind + 1
                  bufsend_j(ind) = constraints(shift + i,j)
                  ! revert sign so that subdomains i and j are getting the same input
                  bufsend_j(ind) = -bufsend_j(ind)
               end do
            end do

            ! distribute chunks of eigenvectors
            if (lbufsend_i .gt. 0) then
               ireq = ireq + 1
               call MPI_ISEND(bufsend_i,lbufsend_i,MPI_DOUBLE_PRECISION,comm_myplace1,comm_myisub,comm,request(ireq),ierr)
            end if

            if (lbufsend_j .gt. 0) then
               ireq = ireq + 1
               call MPI_ISEND(bufsend_j,lbufsend_j,MPI_DOUBLE_PRECISION,comm_myplace2,comm_myjsub,comm,request(ireq),ierr)
            end if
         end if
         do iinstr = 1,ninstructions
            owner = instructions(iinstr,1)
            isub  = instructions(iinstr,2)
            gglob = instructions(iinstr,3)
            call dd_get_interface_size(myid,isub,ndofi,nnodi)

            if (lbufa(iinstr) .gt. 0) then
               ireq = ireq + 1
               call MPI_IRECV(bufrecv(kbufsend(iinstr)),lbufa(iinstr),MPI_DOUBLE_PRECISION,owner,isub,comm,request(ireq),ierr)
            end if
         end do
         nreq = ireq
         call MPI_WAITALL(nreq, request, statarray, ierr)
         if (debug) then
            write(*,*)  'myid =',myid,'All messages in pack 11 received, MPI is fun!.'
            call flush(6)
         end if

         ! processors now own data of adaptively found constraints on their subdomains, they have to filter globs and load them into the structure
         do iinstr = 1,ninstructions
            owner = instructions(iinstr,1)
            isub  = instructions(iinstr,2)
            gglob = instructions(iinstr,3)
            call dd_get_interface_size(myid,isub,ndofi,nnodi)

            nadaptive_rcv = lbufa(iinstr)/ndofi
            lcadapt1 = ndofi
            lcadapt2 = nadaptive_rcv
            allocate(cadapt(lcadapt1,lcadapt2))

            ! copy matrix of constraints from MPI buffer to temporary matrix
            pointbuf = kbufsend(iinstr) - 1
            ind = 0
            do j = 1,nadaptive_rcv
               do i = 1,ndofi
                  ind = ind + 1
                  cadapt(i,j) = bufrecv(pointbuf + ind)
               end do
            end do

            !write(90+myid,*) 'myid = ',myid, 'cadapt:'
            !do i = 1,ndofi
            !   write(90+myid,'(10f16.7)') (cadapt(i,j),j = 1,lcadapt2)
            !end do

            ! load constraints into DD structure 
            call dd_load_adaptive_constraints(isub,gglob,cadapt,lcadapt1,lcadapt2)

            deallocate(cadapt)
         end do

         deallocate(lbufa)
         deallocate(kbufsend) 
         deallocate(bufrecv,bufsend)
         if (my_pair.ge.0) then
            deallocate(constraints)
            deallocate(eigvec,eigval)
            deallocate(bufsend_i,bufsend_j)
            deallocate(weight)
            deallocate(pairslavery)
            deallocate(work)
            deallocate(common_interface)
            deallocate(iingn_i,iingn_j)
            deallocate(rhoi_i,rhoi_j)
            deallocate(tau)
            deallocate(kdofi_i,kdofi_j)
            deallocate(dij)
            deallocate(common_corners)
            deallocate(nndfi_i,nndfi_j)
            deallocate(icnsin_i,icnsin_j)
            deallocate(global_corner_number_i,global_corner_number_j)
         end if
      end do ! loop over rounds of eigenvalue pairs

      deallocate(statarray)
      deallocate(request)
      deallocate(pair_data)
      deallocate(active_pairs)
      deallocate(instructions)

end subroutine

!***************************************
subroutine adaptivity_fake_lobpcg_driver
!***************************************
! Subroutine for looping to compute all Schur complement multiplications in eigenproblems
      implicit none

! local variables
      integer :: idoper

      ! these have no meaning and are present only for matching the arguments
      integer :: n, lx = 0, ly = 0
      real(kr) :: x , y


      ! repeat loop until idmat not equal 3
      ! here the only important argument is idoper
      ! all other necessary data are obtained through modules
      idoper = 3
      do 
         call lobpcg_mvecmult_f(n,x,lx,y,ly,idoper)
         ! stopping criterion
         if (idoper .eq. -3) then
            exit
         end if
      end do

end subroutine

!*************************************************
subroutine adaptivity_mvecmult(n,x,lx,y,ly,idoper)
!*************************************************
! realizes the multiplications with the strange matrices in the local eigenvalue problems for adaptive BDDC

use module_dd
use module_utils
implicit none
include "mpif.h"

real(kr),external :: ddot

! length of vector
integer,intent(in) ::   n 
integer,intent(in) ::  lx 
real(kr),intent(in) ::  x(lx)
integer,intent(in) ::  ly 
real(kr),intent(out) :: y(ly)
! determine which matrix should be multiplied
!  1 - A = P(I-RE)'S(I_RE)P
!  2 - B = PSP
!  3 - not called from LOBPCG by from fake looper - just perform demanded multiplications by S
!  -3 - set on exit if iterational process should be stopped now
integer,intent(inout) ::  idoper 

! local
integer :: i
integer :: iinstr, isub, owner, point1, point2, point, length, do_i_compute, is_active

real(kr) :: xaux(lx)
real(kr) :: xaux2(lx)
real(kr) :: xaux3(lx)

logical :: i_am_slave, all_slaves

integer :: ireq, nreq, ierr

comm_calls = comm_calls + 1

! Check if all processors simply called the fake routine - if so, finalize
if (idoper.eq.3) then
   i_am_slave = .true.
else
   i_am_slave = .false.
end if
call MPI_ALLREDUCE(i_am_slave,all_slaves,1,MPI_LOGICAL,MPI_LAND,comm_comm,ierr)
if (all_slaves) then
   idoper = -3
   return
end if

! Determine active pairs - regardless where I was called from
ireq = 0
do iinstr = 1,ninstructions
   owner = instructions(iinstr,1)
   isub  = instructions(iinstr,2)

   ireq = ireq + 1
   call MPI_IRECV(instructions(iinstr,4),1,MPI_INTEGER,owner,isub,comm_comm,request(ireq),ierr)
end do
if (idoper.eq.3) then
   do_i_compute = 0
else
   do_i_compute = 1
end if
ireq = ireq + 1
call MPI_ISEND(do_i_compute,1,MPI_INTEGER,comm_myplace1,comm_myisub,comm_comm,request(ireq),ierr)
ireq = ireq + 1
call MPI_ISEND(do_i_compute,1,MPI_INTEGER,comm_myplace2,comm_myjsub,comm_comm,request(ireq),ierr)
nreq = ireq
call MPI_WAITALL(nreq, request, statarray, ierr)

ireq = 0
! common operations for A and B - checks and projection to null(D_ij)
if (idoper.eq.1 .or. idoper.eq.2) then
   ! check the dimensions
   if (n .ne. problemsize) then
       write(*,*) 'ADAPTIVITY_LOBPCG_MVECMULT: Vector size mismatch.'
       call error_exit
   end if
   if (lx .ne. ly) then
       write(*,*) 'ADAPTIVITY_LOBPCG_MVECMULT: Data size mismatch: lx,ly',lx,ly
       call error_exit
   end if
   if (lx .ne. n) then
       write(*,*) 'ADAPTIVITY_LOBPCG_MVECMULT: Data size mismatch. lx, n', lx, n
       call error_exit
   end if

   ! make temporary copy
   do i = 1,lx
      xaux(i) = x(i)
   end do
   
   ! xaux = P xaux
   call adaptivity_apply_null_projection(xaux,lx)

   if (idoper.eq.1) then
      ! apply (I-RE) xaux = (I-RR'D_P) xaux
      !  - apply weights
      call adaptivity_apply_weights(weight,lweight,xaux,lx,xaux2,lx)
      !  - apply the R^T operator
      call adaptivity_apply_RT(pairslavery,lpairslavery,xaux2,lx)
      !  - apply the R operator
      call adaptivity_apply_R(pairslavery,lpairslavery,xaux2,lx)
      ! Ix - REx
      do i = 1,problemsize
         xaux(i) = xaux(i) - xaux2(i)
      end do
   end if


   ! distribute sizes of chunks of eigenvectors
   ireq = ireq + 1
   point1 = 1
   call MPI_ISEND(xaux(point1),ndofi_i,MPI_DOUBLE_PRECISION,comm_myplace1,comm_myisub,comm_comm,request(ireq),ierr)

   ireq = ireq + 1
   point2 = ndofi_i + 1
   call MPI_ISEND(xaux(point2),ndofi_j,MPI_DOUBLE_PRECISION,comm_myplace2,comm_myjsub,comm_comm,request(ireq),ierr)
end if

! What follows is performed independently on from where I was called
! obtain chunks of eigenvectors to multiply them
do iinstr = 1,ninstructions
   owner     = instructions(iinstr,1)
   isub      = instructions(iinstr,2)
   is_active = instructions(iinstr,4)

   if (is_active .eq. 1) then
      ireq = ireq + 1
      call MPI_IRECV(bufrecv(kbufsend(iinstr)),lbufa(iinstr),MPI_DOUBLE_PRECISION,owner,isub,comm_comm,request(ireq),ierr)
   end if
end do
nreq = ireq
call MPI_WAITALL(nreq, request, statarray, ierr)

! Multiply subdomain vectors by Schur complement
do iinstr = 1,ninstructions
   owner     = instructions(iinstr,1)
   isub      = instructions(iinstr,2)
   is_active = instructions(iinstr,4)

   if (is_active .eq. 1) then
      point  = kbufsend(iinstr)
      length = lbufa(iinstr)
      call dd_multiply_by_schur(comm_myid,isub,bufrecv(point),length,bufsend(point),length)
   end if
end do

! distribute multiplied vectors
ireq = 0
do iinstr = 1,ninstructions
   owner     = instructions(iinstr,1)
   isub      = instructions(iinstr,2)
   is_active = instructions(iinstr,4)

   if (is_active .eq. 1) then
      ireq = ireq + 1
      call MPI_ISEND(bufsend(kbufsend(iinstr)),lbufa(iinstr),MPI_DOUBLE_PRECISION,owner,isub,comm_comm,request(ireq),ierr)
   end if
end do

! Continue only of I was called from LOBPCG
if (idoper.eq.1 .or. idoper.eq.2) then
   ireq = ireq + 1
   point1 = 1
   call MPI_IRECV(xaux(point1),ndofi_i,MPI_DOUBLE_PRECISION,comm_myplace1,comm_myisub,comm_comm,request(ireq),ierr)

   ireq = ireq + 1
   point2 = ndofi_i + 1
   call MPI_IRECV(xaux(point2),ndofi_j,MPI_DOUBLE_PRECISION,comm_myplace2,comm_myjsub,comm_comm,request(ireq),ierr)
end if

! Wait for all vectors reach their place
nreq = ireq
call MPI_WAITALL(nreq, request, statarray, ierr)

! Continue only of I was called from LOBPCG
if (idoper.eq.1 .or. idoper.eq.2) then

   ! reverse sign of vector if this is the A (required for LOBPCG)
   if (idoper.eq.1) then
      do i = 1,problemsize
         xaux(i) = -xaux(i)
      end do
   end if

   if (idoper.eq.1) then
      ! apply (I-RE)^T = (I - E^T R^T) = (I - D_P^T * R * R^T)
      ! copy the array
      do i = 1,problemsize
         xaux2(i) = xaux(i)
      end do
      !  - apply the R^T operator
      call adaptivity_apply_RT(pairslavery,lpairslavery,xaux2,lx)
      !  - apply the R operator
      call adaptivity_apply_R(pairslavery,lpairslavery,xaux2,lx)
      !  - apply weights
      call adaptivity_apply_weights(weight,lweight,xaux2,lx,xaux3,lx)
      ! Ix - E'R'x
      do i = 1,problemsize
         xaux(i) = xaux(i) - xaux3(i)
      end do
   end if

   ! xaux = P xaux
   call adaptivity_apply_null_projection(xaux,lx)

   ! copy result to y
   do i = 1,lx
      y(i) = xaux(i)
   end do

end if

end subroutine

!****************************************************
subroutine adaptivity_apply_null_projection(vec,lvec)
!****************************************************
! Subroutine for application of projection onto null(D_ij)
! P = I - D' (DD')^-1 D = I - Q_1Q_1'
! where D' = QR = [ Q_1 | Q_2 ] R, is the FULL QR decomposition of D', Q_1 has
! n columns, which corresponds to number of rows in D
      use module_utils, only: zero

      implicit none
      integer, intent(in) ::    lvec
      real(kr), intent(inout) :: vec(lvec)

! local variables
      integer :: lddij, ldvec, lapack_info

      ! apply prepared projection using LAPACK as P = (I-Q_1Q_1') as P = Q_2Q_2',
      ! where Q
      lddij = ldij1
      ldvec = problemsize
      ! xaux = Q_2' * xaux
      call DORMQR( 'Left', 'Transpose',     ldij1, 1, ldij2, dij, lddij, &
                   tau, vec, ldvec, &
                   work,lwork, lapack_info)
      ! put zeros in first N positions of vec
      call zero(vec,ldij2)
      ! xaux = Q_2 * xaux
      call DORMQR( 'Left', 'Non-Transpose', ldij1, 1, ldij2, dij, lddij, &
                   tau, vec, ldvec, &
                   work,lwork, lapack_info)
end subroutine

!**************************************************************************
subroutine adaptivity_apply_weights(dp,ldp,vec_in,lvec_in,vec_out,lvec_out)
!**************************************************************************
! Subroutine for application of weight matrix D_P on vector VEC 
! D_P - stored as diagonal
! vec_out = D_P * vec_in
      use module_utils, only : error_exit
      implicit none
      integer, intent(in)  :: ldp
      real(kr), intent(in) ::  dp(ldp)
      integer, intent(in)  :: lvec_in
      real(kr), intent(in) ::  vec_in(lvec_in)
      integer, intent(in)  :: lvec_out
      real(kr), intent(out) :: vec_out(lvec_out)

! local variables
      integer :: i

      ! check the length of vector for data
      if (ldp .ne. lvec_in .or. lvec_out .ne. lvec_in) then
         write(*,*) 'ADAPTIVITY_APPLY_WEIGHTS: Data size does not match.'
         call error_exit
      end if

      do i = 1,lvec_in
         vec_out(i) = dp(i) * vec_in(i) 
      end do
end subroutine

!*******************************************************
subroutine adaptivity_apply_R(slavery,lslavery,vec,lvec)
!*******************************************************
! Subroutine for application of operator R from domain decomposition
! R : W_hat -> W
! copies values of master variables to slave variables
! vec = R * vec
      use module_utils, only : error_exit
      implicit none
      integer, intent(in) :: lslavery
      integer, intent(in) ::  slavery(lslavery)
      integer, intent(in)  ::   lvec
      real(kr), intent(inout) :: vec(lvec)

! local variables
      integer :: i, indi

      ! check the length of vector for data
      if (lslavery .ne. lvec ) then
         write(*,*) 'ADAPTIVITY_APPLY_R: Data size does not match.'
         call error_exit
      end if

      do i = 1,lvec
         if (slavery(i).ne.0.and.slavery(i).ne.i) then
            indi = pairslavery(i)
            vec(i) = vec(indi)
         end if
      end do
end subroutine

!********************************************************
subroutine adaptivity_apply_RT(slavery,lslavery,vec,lvec)
!********************************************************
! Subroutine for application of operator R^T from domain decomposition
! R^T : W -> W_hat
! works as sum of disconnected entries on master variables, 
! slave variables are meaningless after the action
! vec = R' * vec
      use module_utils, only : error_exit
      implicit none
      integer, intent(in) :: lslavery
      integer, intent(in) ::  slavery(lslavery)
      integer, intent(in)  ::   lvec
      real(kr), intent(inout) :: vec(lvec)

! local variables
      integer :: i, indi

      ! check the length of vector for data
      if (lslavery .ne. lvec ) then
         write(*,*) 'ADAPTIVITY_APPLY_RT: Data size does not match.'
         call error_exit
      end if

      do i = 1,lvec
         if (slavery(i).ne.0.and.slavery(i).ne.i) then
            indi = slavery(i)
            vec(indi) = vec(indi) + vec(i) 
         end if
      end do
end subroutine

!***************************************************************
subroutine adaptivity_get_pair_data(idpair,pair_data,lpair_data)
!***************************************************************
! Subroutine for getting info about pairs to the global structure
      use module_utils
      implicit none

! pair number
      integer,intent(in) :: idpair
! length of vector for data
      integer,intent(in) :: lpair_data
! vector of data for pair IDPAIR
      integer,intent(out) :: pair_data(lpair_data)

! local variables
      integer :: i

      ! check the length of vector for data
      if (lpair_data .ne. lpair_subdomains2) then
         write(*,*) 'ADAPTIVITY_GET_PAIR_DATA: Size not sufficient for getting info about pair.'
         call error_exit
      end if
      ! check that the info about pair is available
      if (.not.allocated(pair_subdomains)) then
         write(*,*) 'ADAPTIVITY_GET_PAIR_DATA: Structure with global pair data is not allocated.'
         call error_exit
      end if
      if (pair_subdomains(idpair,1).eq.-1) then
         write(*,*) 'ADAPTIVITY_GET_PAIR_DATA: Incomplete information about pair - processor not assigned.'
         call error_exit
      end if
      if (any(pair_subdomains(idpair,2:).eq.0)) then
         write(*,*) 'ADAPTIVITY_GET_PAIR_DATA: Incomplete information about pair - zeros in subdomain data.'
         call error_exit
      end if

      ! after checking, get info about pair from the global structure and load it
      do i = 1,lpair_data
         pair_data(i) = pair_subdomains(idpair,i)
      end do

end subroutine


!**************************************
subroutine adaptivity_print_pairs(myid)
!**************************************
! Subroutine for printing data about pairs to screen
      implicit none

! number of processor
      integer,intent(in) :: myid

! local variables
      integer :: ipair, j

      write(*,*) 'Info about loaded pairs on processor ',myid,',',lpair_subdomains1,' pairs loaded:'
      if (allocated(pair_subdomains)) then
         do ipair = 1,lpair_subdomains1
            write(*,'(6i10)') (pair_subdomains(ipair,j),j = 1,lpair_subdomains2)
         end do
      else 
         write(*,*) 'ADAPTIVITY_PRINT_PAIRS: Array of pairs is not allocated.'
      end if
end subroutine

!*****************************
subroutine adaptivity_finalize
!*****************************
! Subroutine for finalization of adaptivity
      implicit none

! clear memory
      if (allocated(pair_subdomains)) then
         deallocate(pair_subdomains)
      end if

      return
end subroutine

end module module_adaptivity


