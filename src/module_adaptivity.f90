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

module module_adaptivity
!***********************
! Module for adaptive search of constraints for BDDC preconditioner
! Jakub Sistek, Bologna, 11/2010, Praha 1/2011

use module_dd

!=================================================
! basic parameters related to adaptivity
! type of reals
integer,parameter,private  :: kr = kind(1.D0)
! numerical zero
real(kr),parameter,private :: numerical_zero = 1.e-12_kr

! threshold on eigenvalues to define an adaptive constraint
! eigenvectors for eigenvalues above this value are used for creation of constraints 
real(kr),parameter,private :: threshold_eigval_default = 1.5_kr
logical,parameter,private  :: read_threshold_from_file = .false.
! LOBPCG related variables
! maximal number of LOBPCG iterations
integer,parameter,private ::  lobpcg_maxit = 15
! precision of LOBPCG solver - worst residual
real(kr),parameter,private :: lobpcg_rel_tol   = 1.e-9_kr
! maximal number of eigenvectors per problem
! this number is used for sufficient size of problems 
! for small problems, size of glob is used
integer,parameter,private ::  neigvecx     = 10  
! verbosity of LOBPCG solver
! 0 - no output
! 1 - some output
! 2 - maximal output
integer,parameter,private ::  lobpcg_verbosity = 0
! loading old values of initial guess of eigenvectors
! 0 - generate new random vectors in LOBCPCG - may have poor performance
! 1 - generate new random vectors in Fortran and pass these to LOBPCG - better behaviour
integer,parameter,private ::  use_vec_values = 1
! using preconditioner for LOBPCG
! 0 - no preconditioning (default)
! 1 - local BDDC preconditioner (!!! EXPERIMENTAL !!!)
integer,parameter,private ::  lobpcg_preconditioner = 1
! using nullspace projection for eigenproblems - may improve robustness but can increase time dramatically
logical,parameter,private ::  apply_null_projection = .false.
! if using preconditioner and it fails, try again without preconditioner?
logical,parameter,private ::  try_harder = .false.
! loading computed eigenvectors
! T - compute new vectors
! F - load eigenvectors from file
logical,parameter,private :: recompute_vectors = .true.
! debugging 
logical,parameter,private :: debug = .false.
logical,parameter,private :: profile = .true.

! maximal allowed length of file names
integer,parameter,private :: lfnamex = 130

real(kr),private :: threshold_eigval
integer,parameter,private :: idbase = 100 ! basic unit to add myid for independent units for procs

!=================================================

! table of pairs of eigenproblems to compute
! structure:
!  IGLOB | ISUB | JSUB 
integer,private            :: lpair_subdomains1 
integer,parameter,private  :: lpair_subdomains2 = 3
integer,allocatable,private :: pair_subdomains(:,:)

logical,private :: i_compute_pair
logical,private :: i_compute_multiplications

integer,private :: comm_calls

integer,private :: comm_myid

integer,private :: neigvec, problemsize

integer,private :: ndofi_i, ndofi_j
integer,private :: lindrowc_i, lindrowc_j
integer,private :: nnodi_i,nnodi_j
integer ::            lindrowc_adapt_i,   lindrowc_adapt_j
integer,allocatable :: indrowc_adapt_i(:), indrowc_adapt_j(:)

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

! matrix of local nullspace basis nullB
logical ::                     is_nullB_ready = .false.
integer,private ::             ldnullB ! leading dimension
integer,private ::              lnullB1,lnullB2
real(kr),allocatable,private ::  nullB(:,:)

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
integer,private ::             ltau3
real(kr),allocatable,private :: tau3(:)
integer,private ::             lwork3
real(kr),allocatable,private :: work3(:)

! local coarse matrix
integer,private ::             lcoarsem_adapt1
integer,private ::             lcoarsem_adapt2
real(kr),allocatable,private :: coarsem_adapt(:,:)
integer,private ::             lkceigval
real(kr),allocatable,private :: kceigval(:)
integer,private ::              comm_lresc
real(kr),allocatable,private :: comm_resc(:)
real(kr),allocatable,private :: comm_resc_i(:)
real(kr),allocatable,private :: comm_resc_j(:)

! auxiliary arrays (necessary in each LOBPCG iteration)
integer ::             comm_lxaux
real(kr),allocatable :: comm_xaux(:)
integer ::             comm_lxaux2
real(kr),allocatable :: comm_xaux2(:)

contains

!************************************************************
subroutine adaptivity_init(comm,pairs,lpairs1,lpairs2, npair)
!************************************************************
! Subroutine for initialization of adaptive search of constraints
      use module_utils
      implicit none
      include "mpif.h"

! communicator
      integer,intent(in) :: comm
! pairs of subdomains
      integer,intent(in) :: lpairs1
      integer,intent(in) :: lpairs2
      integer,intent(in) :: pairs(lpairs1,lpairs2)
! number of pairs
      integer,intent(in) :: npair

! local variables
! MPI variables
      integer :: myid
      integer :: ierr

      integer :: ipair, j
      integer :: idthresh
      real(kr) :: threshold_file

      ! orient in the communicator
      call MPI_COMM_RANK(comm,myid,ierr)

! copy pairs to module array
      lpair_subdomains1 = npair
      allocate(pair_subdomains(lpair_subdomains1,lpair_subdomains2))
      call zero(pair_subdomains,lpair_subdomains1,lpair_subdomains2)
      do ipair = 1,npair
         ! first column is associated with processors - initialize it to -1 = no processor assigned
         do j = 1,3
            pair_subdomains(ipair,j) = pairs(ipair,j)
         end do
      end do
 
! set threshold
      if (read_threshold_from_file) then
         ! root reads the value
         if (myid.eq.0) then
            call allocate_unit(idthresh)
            open (unit=idthresh,file='threshold.txt',status='old',form='formatted')

            read(idthresh,*) threshold_file
            close(idthresh)
         end if
         ! broadcast it among processors
!*****************************************************************MPI
         call MPI_BCAST(threshold_file,1, MPI_DOUBLE_PRECISION, 0, comm, ierr)
!*****************************************************************MPI
         ! set the threshold
         threshold_eigval = threshold_file
      else
         threshold_eigval = threshold_eigval_default
      end if

      if (myid.eq.0) then
         write(*,*) 'ADAPTIVITY SETUP====================='
         write(*,*) 'threshold for selection: ', threshold_eigval
         write(*,*) 'max LOBPCG iterations: ',  lobpcg_maxit
         write(*,*) 'LOBPCG relative tolerance: ',  lobpcg_rel_tol
         write(*,*) 'max number of computed eigenvectors: ', neigvecx
         write(*,*) 'LOBPCG preconditioner (0 = no,1 = BDDC): ', lobpcg_preconditioner
         write(*,*) 'END ADAPTIVITY SETUP================='
      end if

end subroutine

!******************************************************************************************************************
subroutine adaptivity_get_active_pairs(iround,nproc,pair2proc,lpair2proc, active_pairs,lactive_pairs,nactive_pairs)
!******************************************************************************************************************
! Subroutine for activating and deactivating pairs
      use module_utils
      implicit none

! number of round 
      integer,intent(in) :: iround
! number of processors
      integer,intent(in) :: nproc
! distribution of pairs
      integer,intent(in) :: lpair2proc
      integer,intent(in) ::  pair2proc(lpair2proc)
! indices of active pairs
      integer,intent(in) :: lactive_pairs
      integer,intent(out) :: active_pairs(lactive_pairs)
! number of active pairs
      integer, intent(out) :: nactive_pairs

! local variables
      character(*),parameter:: routine_name = 'ADAPTIVITY_GET_ACTIVE_PAIRS'
      integer :: iproc, indpair

      ! check length
      if (lpair2proc.ne.nproc+1) then
         call error(routine_name,'array dimension mismatch for pairs')
      end if
      if (lactive_pairs.ne.nproc) then
         call error(routine_name,'array dimension mismatch for active pairs')
      end if

      nactive_pairs = 0
      do iproc = 0,nproc-1
         indpair = pair2proc(iproc+1) + iround - 1
         if (indpair.lt.pair2proc(iproc+2)) then
            nactive_pairs = nactive_pairs + 1
            active_pairs(iproc+1) = indpair
         else
            active_pairs(iproc+1) = -1 
         end if
      end do

end subroutine

!******************************************************************************************
subroutine adaptivity_solve_eigenvectors(suba,lsuba,sub2proc,lsub2proc,indexsub,lindexsub,&
                                         pair2proc,lpair2proc,comm_all,&
                                         use_explicit_schurs,matrixtype, est)
!******************************************************************************************
! Subroutine for parallel solution of distributed eigenproblems
      use module_dd
      use module_pp
      use module_utils
      implicit none
      include "mpif.h"

! array of sub structure
      integer,intent(in) ::                lsuba
      type(subdomain_type),intent(inout) :: suba(lsuba)

! division of subdomains to processors
      integer,intent(in) :: lsub2proc
      integer,intent(in) ::  sub2proc(lsub2proc)
! global indices of local subdomains
      integer,intent(in) :: lindexsub
      integer,intent(in) ::  indexsub(lindexsub)
! division of pairs to processors (independent of subdomains)
      integer,intent(in) :: lpair2proc
      integer,intent(in) ::  pair2proc(lpair2proc)

! communicator global
      integer,intent(in) :: comm_all
! should explicit Schur complements be used and sent?
      logical,intent(in) :: use_explicit_schurs
! type of matrix
      integer,intent(in) :: matrixtype

! prediction of condition number
      real(kr),intent(out) :: est


! local variables
      character(*),parameter:: routine_name = 'ADAPTIVITY_SOLVE_EIGENVECTORS'
      integer :: isub, jsub, iround, isub_loc
! Maximal local number of eigenproblems at one processor
      integer :: npair_locx
! Global number of pairs to compute eigenproblems
      integer :: npair

      integer :: indpair, iproc
      integer :: gglob, my_pair, &
                 nactive_pairs, owner,&
                 place1, place2, pointbuf, i, j, iinstr, icommon,&
                 indi_i, indi_j, ndofn, inodi,&
                 shift, indcommon, ndofcomm, idofn, irhoicomm,&
                 point_i, point_j, indiv, nadaptive, ioper, nadaptive_rcv, ind, &
                 nnzc_i, nnzc_j, inddrow, ic, indirow, indirow_loc, indjcol_loc, &
                 neigvecf, problemsizef
      integer :: nvalid_i, nvalid_j, nvalid
      integer :: idmyunit

      real(kr) :: trA, lobpcg_tol

      integer :: ndofi
      integer :: nnodi
      integer :: nnzc

      character(lfnamex) :: filename

      integer ::            lpair_data
      integer,allocatable :: pair_data(:)

      ! numbers of active pairs
      integer ::            lactive_pairs
      integer,allocatable :: active_pairs(:)

      ! eigenvectors and eigenvalues
      integer ::             leigvec
      real(kr),allocatable :: eigvec(:)
      integer ::             leigval
      real(kr),allocatable :: eigval(:)

      ! coarse nodes information
      integer ::            ncommon_crows
      integer ::            lindrowc_i,   lindrowc_j
      integer,allocatable :: indrowc_i(:), indrowc_j(:)
      integer ::            lindrowc
      integer,allocatable :: indrowc(:)
      integer ::            lc_sparse_i
      integer,allocatable :: i_c_sparse_i(:),j_c_sparse_i(:)
      real(kr),allocatable ::  c_sparse_i(:)
      integer ::            lc_sparse_j
      integer,allocatable :: i_c_sparse_j(:),j_c_sparse_j(:)
      real(kr),allocatable ::  c_sparse_j(:)
      integer ::            lc_sparse
      integer,allocatable :: i_c_sparse(:),j_c_sparse(:)
      real(kr),allocatable ::  c_sparse(:)
      integer ::            lcommon_crows
      integer,allocatable :: common_crows(:)
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

      ! regularization of the problem
      integer ::              lphisi1_i, lphisi2_i
      real(kr),allocatable ::  phisi_i(:,:)
      integer ::              lphisi1_j, lphisi2_j
      real(kr),allocatable ::  phisi_j(:,:)
      integer ::              lphisi1, lphisi2
      real(kr),allocatable ::  phisi(:,:)
      integer ::              lz1, lz2
      real(kr),allocatable ::  z(:,:)
      integer ::              shifti, shiftj
      integer ::              lbz1, lbz2
      real(kr),allocatable ::  bz(:,:)
      integer ::              lzbz1, lzbz2
      real(kr),allocatable ::  zbz(:,:)
      integer ::              lzbzeigval
      real(kr),allocatable ::  zbzeigval(:)
      real(kr) ::              null_tol
      integer ::               null_dim

      ! local BDDC
      integer :: irow, jcol
      integer ::              lcoarsem
      real(kr),allocatable ::  coarsem(:)
      integer ::              lcoarsem_i
      real(kr),allocatable ::  coarsem_i(:)
      integer ::              lcoarsem_j
      real(kr),allocatable ::  coarsem_j(:)
      integer ::              icoarsem
      integer :: i_loc, indg
      integer ::            lcoarsem_embed
      integer,allocatable :: coarsem_embed(:)
      integer :: i_upp, j_upp, indc_loc
      integer :: no_prec


      integer ::             lcadapt1, lcadapt2
      real(kr),allocatable :: cadapt(:,:)

      integer ::            ncommon_interface
      integer ::            lcommon_interface
      integer,allocatable :: common_interface(:)

      ! data for computing with explicit Schur complements
      integer :: lschur1,lschur2
      real(kr),allocatable :: schur(:,:)
      integer :: lschur1_i,lschur2_i
      real(kr),allocatable :: schur_i(:,:)
      integer :: lschur1_j,lschur2_j
      real(kr),allocatable :: schur_j(:,:)
      integer ::             lmata1,lmata2
      real(kr),allocatable :: mata(:,:)
      integer ::             lmatb1,lmatb2
      real(kr),allocatable :: matb(:,:)
      integer ::             ldoubleschur1,ldoubleschur2
      real(kr),allocatable :: doubleschur(:,:)
      integer ::             lxaux
      real(kr),allocatable :: xaux(:)
      integer ::             lxaux2
      real(kr),allocatable :: xaux2(:)
      integer :: jrcol

      real(kr),external :: ddot

      ! LAPACK QR related variables
      integer :: lapack_info
      ! LAPACK eigenproblems
      integer::              lwork2
      real(kr),allocatable :: work2(:)
      integer::              leiglap
      real(kr),allocatable :: eiglap(:)
      real(kr) :: stab

      ! LOBPCG related variables
      integer ::  lobpcg_iter ! final number of iterations
      real(kr)::   est_round, est_loc

      ! MPI related variables
      integer :: myid, nproc, ierr, ireq, nreq

      ! time vars
      real(kr) :: time, &
                  time_obtain_accu = 0._kr, &
                  time_null_prep_accu = 0._kr, &
                  time_null_mult_accu = 0._kr, &
                  time_null_comp_accu = 0._kr, &
                  time_precond_accu = 0._kr, &
                  time_solve_accu = 0._kr, &
                  time_postp_accu = 0._kr

      ! orient in the communicator
      call MPI_COMM_RANK(comm_all,myid,ierr)
      call MPI_COMM_SIZE(comm_all,nproc,ierr)
      
      ! find maximal number of pairs per proc
      npair_locx = maxval(pair2proc(2:nproc+1) - pair2proc(1:nproc))
      npair = pair2proc(nproc+1)-1

      ! allocate table for work instructions - the worst case is that in each
      ! round, I have to compute all the subdomains, i.e. 2 for each pair
      linstructions1 = 2*nproc
      allocate(instructions(linstructions1,linstructions2))
      ! prepare MPI arrays 
      lrequest = 2*nproc + 2*3
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
      ! unit for local reading/writing from/to file
      idmyunit = idbase + myid

      ! Loop over rounds of eigenvalue solves
      do iround = 1,npair_locx
         comm_calls = 0

         ! each round of eigenproblems has its structure - determine active pairs
         call adaptivity_get_active_pairs(iround,nproc,pair2proc,lpair2proc, active_pairs,lactive_pairs,nactive_pairs)

         ! determine which pair I compute
         my_pair = active_pairs(myid+1)
         if (my_pair.gt.0) then
            i_compute_pair = .true.
         else
            i_compute_pair = .false.
         end if

         call zero(pair_data,lpair_data)
         if (i_compute_pair) then

            call adaptivity_get_pair_data(my_pair,pair_data,lpair_data)
   
            comm_mygglob = pair_data(1)
            comm_myisub  = pair_data(2)
            comm_myjsub  = pair_data(3)

            ! where are these subdomains ?
            call pp_get_proc_for_sub(comm_myisub,comm_all,sub2proc,lsub2proc,comm_myplace1)
            call pp_get_proc_for_sub(comm_myjsub,comm_all,sub2proc,lsub2proc,comm_myplace2)
         else
            comm_myisub   = -1
            comm_myjsub   = -1
            comm_myplace1 = -1
            comm_myplace2 = -1
         end if

         ! determine working instructions for sending subdomain matrices
         ! go through pairs that are active in this round
         instructions  = 0
         ninstructions = 0
         do iproc = 0,nproc-1
            indpair = active_pairs(iproc+1)
            
            if (indpair.gt.0) then
               call adaptivity_get_pair_data(indpair,pair_data,lpair_data)
               gglob  = pair_data(1)
               isub   = pair_data(2)
               jsub   = pair_data(3)

               ! where are these subdomains ?
               call pp_get_proc_for_sub(isub,comm_all,sub2proc,lsub2proc,place1)
               call pp_get_proc_for_sub(jsub,comm_all,sub2proc,lsub2proc,place2)

               if (myid.eq.place1) then
                  ! add instruction
                  ninstructions = ninstructions + 1

                  ! who I will send data to
                  instructions(ninstructions,1) = iproc
                  ! subdomain number
                  instructions(ninstructions,2) = isub
                  ! global glob number
                  instructions(ninstructions,3) = gglob
               end if

               if (myid.eq.place2) then
                  ! add instruction
                  ninstructions = ninstructions + 1

                  ! who I will send data to
                  instructions(ninstructions,1) = iproc
                  ! subdomain number
                  instructions(ninstructions,2) = jsub
                  ! global glob number
                  instructions(ninstructions,3) = gglob
               end if
            end if
         end do
         if (ninstructions.gt.0) then
            i_compute_multiplications = .true.
         else
            i_compute_multiplications = .false.
         end if

         ! the scheme for communication is ready
         ! debug
         !write(idmyunit,*) 'myid =',myid, 'instructions:'
         !do i = 1,ninstructions
         !   write(idmyunit,*) instructions(i,:)
         !end do
         !call flush(idmyunit)

!-----profile
         if (profile) then
            call MPI_BARRIER(comm_all,ierr)
            call time_start
         end if
!-----profile
         ! build the local matrix of projection on common globs for active pair
         !  get sizes of interface of subdomains in my problem
         ireq = 0
         if (i_compute_pair) then
            ! receive sizes of interfaces of subdomains involved in my problem

            ireq = ireq + 1
            call MPI_IRECV(ndofi_i,1,MPI_INTEGER,comm_myplace1,comm_myisub,comm_all,request(ireq),ierr)

            ireq = ireq + 1
            call MPI_IRECV(ndofi_j,1,MPI_INTEGER,comm_myplace2,comm_myjsub,comm_all,request(ireq),ierr)
         end if
         ! send sizes of subdomains involved in problems
         do iinstr = 1,ninstructions
            owner = instructions(iinstr,1)
            isub  = instructions(iinstr,2)
            call get_index(isub,indexsub,lindexsub,isub_loc)

            call dd_get_interface_size(suba(isub_loc),ndofi,nnodi)

            call MPI_SEND(ndofi,1,MPI_INTEGER,owner,isub,comm_all,ierr)
         end do
         nreq = ireq
         if (nreq.gt.0) then
            call MPI_WAITALL(nreq, request, statarray, ierr)
         end if
         ! debug
         !call info(routine_name,'All messages in pack 1 received on proc',myid)
         !if (i_compute_pair) then
         !   write(idmyunit,*) 'ndofi_i',ndofi_i,'ndofi_j',ndofi_j
         !end if

         !  get number of rows in C and find common subset
         ireq = 0
         if (i_compute_pair) then
            ! receive numbers of coarse nodes

            ireq = ireq + 1
            call MPI_IRECV(lindrowc_i,1,MPI_INTEGER,comm_myplace1,comm_myisub,comm_all,request(ireq),ierr)

            ireq = ireq + 1
            call MPI_IRECV(lindrowc_j,1,MPI_INTEGER,comm_myplace2,comm_myjsub,comm_all,request(ireq),ierr)
         end if
         ! send sizes of subdomains involved in problems
         do iinstr = 1,ninstructions
            owner = instructions(iinstr,1)
            isub  = instructions(iinstr,2)
            call get_index(isub,indexsub,lindexsub,isub_loc)

            call dd_get_number_of_crows(suba(isub_loc),lindrowc)

            call MPI_SEND(lindrowc,1,MPI_INTEGER,owner,isub,comm_all,ierr)
         end do
         nreq = ireq
         if (nreq.gt.0) then
            call MPI_WAITALL(nreq, request, statarray, ierr)
         end if
         ! debug
         !call info(routine_name,'All messages in pack 2 received on proc',myid)
         !if (i_compute_pair) then
         !   write(idmyunit,*) 'lindrowc_i',lindrowc_i,'lindrowc_j',lindrowc_j
         !   call flush(idmyunit)
         !end if

         !  get number of nonzeros in C 
         ireq = 0
         if (i_compute_pair) then
            ! receive numbers of coarse nodes

            ireq = ireq + 1
            call MPI_IRECV(nnzc_i,1,MPI_INTEGER,comm_myplace1,comm_myisub,comm_all,request(ireq),ierr)

            ireq = ireq + 1
            call MPI_IRECV(nnzc_j,1,MPI_INTEGER,comm_myplace2,comm_myjsub,comm_all,request(ireq),ierr)
         end if
         ! send sizes of subdomains involved in problems
         do iinstr = 1,ninstructions
            owner = instructions(iinstr,1)
            isub  = instructions(iinstr,2)
            call get_index(isub,indexsub,lindexsub,isub_loc)

            call dd_get_number_of_cnnz(suba(isub_loc),nnzc)

            call MPI_SEND(nnzc,1,MPI_INTEGER,owner,isub,comm_all,ierr)
         end do
         nreq = ireq
         if (nreq.gt.0) then
            call MPI_WAITALL(nreq, request, statarray, ierr)
         end if
         ! debug
         !call info(routine_name,'All messages in pack 3 received on proc',myid)
         !if (i_compute_pair) then
         !   write(idmyunit,*) 'nnzc_i',nnzc_i,'nnzc_j',nnzc_j
         !   call flush(idmyunit)
         !end if

         !  get number of nodes on interface
         ireq = 0
         if (i_compute_pair) then
            ! receive numbers of interface nodes

            ireq = ireq + 1
            call MPI_IRECV(nnodi_i,1,MPI_INTEGER,comm_myplace1,comm_myisub,comm_all,request(ireq),ierr)

            ireq = ireq + 1
            call MPI_IRECV(nnodi_j,1,MPI_INTEGER,comm_myplace2,comm_myjsub,comm_all,request(ireq),ierr)
         end if
         ! send sizes of subdomains involved in problems
         do iinstr = 1,ninstructions
            owner = instructions(iinstr,1)
            isub  = instructions(iinstr,2)
            call get_index(isub,indexsub,lindexsub,isub_loc)

            call dd_get_interface_size(suba(isub_loc),ndofi,nnodi)

            call MPI_SEND(nnodi,1,MPI_INTEGER,owner,isub,comm_all,ierr)
         end do
         nreq = ireq
         if (nreq.gt.0) then
            call MPI_WAITALL(nreq, request, statarray, ierr)
         end if
         ! debug
         !call info(routine_name,'All messages in pack 4 received on proc',myid)
         !if (i_compute_pair) then
         !   write(idmyunit,*) 'nnodi_i',nnodi_i,'nnodi_j',nnodi_j
         !   call flush(idmyunit)
         !end if



         ! allocate space for coarse nodes
         if (i_compute_pair) then
            ! determine problem size
            problemsize = ndofi_i + ndofi_j

            allocate(indrowc_i(lindrowc_i),indrowc_j(lindrowc_j))
            lnndfi_i = nnodi_i
            lnndfi_j = nnodi_j
            allocate(nndfi_i(lnndfi_i),nndfi_j(lnndfi_j))
            lrhoi_i = ndofi_i
            lrhoi_j = ndofi_j
            allocate(rhoi_i(lrhoi_i),rhoi_j(lrhoi_j))
            liingn_i = nnodi_i
            liingn_j = nnodi_j
            allocate(iingn_i(liingn_i),iingn_j(liingn_j))
    ! communicate matrix C
            lc_sparse_i = nnzc_i
            lc_sparse_j = nnzc_j
            allocate(i_c_sparse_i(lc_sparse_i),i_c_sparse_j(lc_sparse_j))
            allocate(j_c_sparse_i(lc_sparse_i),j_c_sparse_j(lc_sparse_j))
            allocate(  c_sparse_i(lc_sparse_i),  c_sparse_j(lc_sparse_j))
         end if

         ! get data about coarse nodes
         ireq = 0
         if (i_compute_pair) then
            ! receive global numbers of coarse nodes

            ireq = ireq + 1
            call MPI_IRECV(indrowc_i,lindrowc_i,MPI_INTEGER,comm_myplace1,comm_myisub,comm_all,request(ireq),ierr)

            ireq = ireq + 1
            call MPI_IRECV(indrowc_j,lindrowc_j,MPI_INTEGER,comm_myplace2,comm_myjsub,comm_all,request(ireq),ierr)
         end if
         ! send sizes of subdomains involved in problems
         do iinstr = 1,ninstructions
            owner = instructions(iinstr,1)
            isub  = instructions(iinstr,2)
            call get_index(isub,indexsub,lindexsub,isub_loc)

            call dd_get_number_of_crows(suba(isub_loc),lindrowc)

            allocate(indrowc(lindrowc))
            call dd_get_subdomain_crows(suba(isub_loc),indrowc,lindrowc)

            call MPI_SEND(indrowc,lindrowc,MPI_INTEGER,owner,isub,comm_all,ierr)
            deallocate(indrowc)
         end do
         nreq = ireq
         if (nreq.gt.0) then
            call MPI_WAITALL(nreq, request, statarray, ierr)
         end if
         ! debug
         !call info(routine_name,'All messages in pack 5 received on proc',myid)
         !if (i_compute_pair) then
         !   write(idmyunit,*) 'indrowc_i'
         !   write(idmyunit,*)  indrowc_i
         !   write(idmyunit,*) 'indrowc_j'
         !   write(idmyunit,*)  indrowc_j
         !   call flush(idmyunit)
         !end if


         ! get data for matrix C
         ireq = 0
         if (my_pair.ge.0) then
            ireq = ireq + 1
            call MPI_IRECV(i_c_sparse_i,lc_sparse_i,MPI_INTEGER,comm_myplace1,comm_myisub,comm_all,request(ireq),ierr)
            ireq = ireq + 1
            call MPI_IRECV(j_c_sparse_i,lc_sparse_i,MPI_INTEGER,comm_myplace1,comm_myisub,comm_all,request(ireq),ierr)
            ireq = ireq + 1
            call MPI_IRECV(c_sparse_i,lc_sparse_i,MPI_DOUBLE_PRECISION,comm_myplace1,comm_myisub,comm_all,request(ireq),ierr)

            ireq = ireq + 1
            call MPI_IRECV(i_c_sparse_j,lc_sparse_j,MPI_INTEGER,comm_myplace2,comm_myjsub,comm_all,request(ireq),ierr)
            ireq = ireq + 1
            call MPI_IRECV(j_c_sparse_j,lc_sparse_j,MPI_INTEGER,comm_myplace2,comm_myjsub,comm_all,request(ireq),ierr)
            ireq = ireq + 1
            call MPI_IRECV(c_sparse_j,lc_sparse_j,MPI_DOUBLE_PRECISION,comm_myplace2,comm_myjsub,comm_all,request(ireq),ierr)
         end if
         do iinstr = 1,ninstructions
            owner = instructions(iinstr,1)
            isub  = instructions(iinstr,2)
            call get_index(isub,indexsub,lindexsub,isub_loc)

            call dd_get_number_of_cnnz(suba(isub_loc),nnzc)

            lc_sparse = nnzc
            allocate(i_c_sparse(lc_sparse),j_c_sparse(lc_sparse),c_sparse(lc_sparse))
            call dd_get_subdomain_c(suba(isub_loc),i_c_sparse,j_c_sparse,c_sparse,lc_sparse)

            call MPI_SEND(i_c_sparse,lc_sparse,MPI_INTEGER,owner,isub,comm_all,ierr)
            call MPI_SEND(j_c_sparse,lc_sparse,MPI_INTEGER,owner,isub,comm_all,ierr)
            call MPI_SEND(c_sparse,lc_sparse,MPI_DOUBLE_PRECISION,owner,isub,comm_all,ierr)
            deallocate(i_c_sparse,j_c_sparse,c_sparse)
         end do
         nreq = ireq
         if (nreq.gt.0) then
            call MPI_WAITALL(nreq, request, statarray, ierr)
         end if
         ! debug
         !call info(routine_name,'All messages in pack 6 received on proc',myid)
         !if (i_compute_pair) then
         !   write(idmyunit,*) 'i_c_sparse_i'
         !   write(idmyunit,*)  i_c_sparse_i
         !   write(idmyunit,*) 'j_c_sparse_i'
         !   write(idmyunit,*)  j_c_sparse_i
         !   write(idmyunit,*) 'c_sparse_i'
         !   write(idmyunit,*)  c_sparse_i
         !   write(idmyunit,*) 'i_c_sparse_j'
         !   write(idmyunit,*)  i_c_sparse_j
         !   write(idmyunit,*) 'j_c_sparse_j'
         !   write(idmyunit,*)  j_c_sparse_j
         !   write(idmyunit,*) 'c_sparse_j'
         !   write(idmyunit,*)  c_sparse_j
         !   call flush(idmyunit)
         !end if

         ! get numbers of dof on interface
         ireq = 0
         if (i_compute_pair) then
            ! receive global numbers of coarse nodes

            ireq = ireq + 1
            call MPI_IRECV(nndfi_i,nnodi_i,MPI_INTEGER,comm_myplace1,comm_myisub,comm_all,request(ireq),ierr)

            ireq = ireq + 1
            call MPI_IRECV(nndfi_j,nnodi_j,MPI_INTEGER,comm_myplace2,comm_myjsub,comm_all,request(ireq),ierr)
         end if
         ! send sizes of subdomains involved in problems
         do iinstr = 1,ninstructions
            owner = instructions(iinstr,1)
            isub  = instructions(iinstr,2)
            call get_index(isub,indexsub,lindexsub,isub_loc)

            call dd_get_interface_size(suba(isub_loc),ndofi,nnodi)

            lnndfi = nnodi
            allocate(nndfi(lnndfi))
            call dd_get_subdomain_interface_nndf(suba(isub_loc),nndfi,lnndfi)

            call MPI_SEND(nndfi,nnodi,MPI_INTEGER,owner,isub,comm_all,ierr)
            deallocate(nndfi)
         end do
         nreq = ireq
         if (nreq.gt.0) then
            call MPI_WAITALL(nreq, request, statarray, ierr)
         end if
         ! debug
         !call info(routine_name,'All messages in pack 7 received on proc',myid)
         !if (my_pair.ge.0) then
         !   write(idmyunit,*) 'nndfi_i'
         !   write(idmyunit,*)  nndfi_i
         !   write(idmyunit,*) 'nndfi_j'
         !   write(idmyunit,*)  nndfi_j
         !   call flush(idmyunit)
         !end if

         ! find intersection of coarse nodes
         if (i_compute_pair) then
            lcommon_crows = min(lindrowc_i,lindrowc_j)
            allocate(common_crows(lcommon_crows))
            call get_array_intersection(indrowc_i,lindrowc_i,&
                                        indrowc_j,lindrowc_j,&
                                        common_crows,lcommon_crows,ncommon_crows)
            ! debug
            !write(idmyunit,*) 'common_crows', common_crows(1:ncommon_crows)
            !call flush(idmyunit)
         end if


         ! build matrix D_ij^T for the pair of subdomains
         if (i_compute_pair) then
            ! prepare space for dense matrix D_ij
            ldij1 = ndofi_i + ndofi_j
            ldij2 = ncommon_crows
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

            ! join C_i and C_j together to make the matrix D_ij = [C_i  -C_j] on common rows
            ! C_i
            do ic = 1,nnzc_i
               indirow_loc = i_c_sparse_i(ic)
               indirow     = indrowc_i(indirow_loc)

               if (any(indirow.eq.common_crows)) then
                  indjcol_loc = j_c_sparse_i(ic)

                  call get_index(indirow,common_crows,lcommon_crows,inddrow) 

                  dij(indjcol_loc,inddrow) = c_sparse_i(ic)
               end if
            end do
            ! - C_j
            do ic = 1,nnzc_j
               indirow_loc = i_c_sparse_j(ic)
               indirow     = indrowc_j(indirow_loc)

               if (any(indirow.eq.common_crows)) then
                  ! shift the index
                  indjcol_loc = j_c_sparse_j(ic) + ndofi_i

                  call get_index(indirow,common_crows,lcommon_crows,inddrow) 

                  dij(indjcol_loc,inddrow) = -c_sparse_j(ic)
               end if
            end do

            ! debug
            !   do i = 1,ldij1
            !      write(*,'(100f4.1)') (dij(i,j),j = 1,ldij2)
            !   end do

        ! Prepare projection onto null D_ij
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
               call error(routine_name,' in LAPACK QR factorization of matrix D_ij^T: ', lapack_info)
            end if
            ! in space of D_ij are now stored factors R and Householder reflectors v
         end if

         ! prepare operator (I-R_ijE_ij)
         ! get data about interface weigths
         ireq = 0
         if (i_compute_pair) then
            ! receive diagonal entries
            ireq = ireq + 1
            call MPI_IRECV(rhoi_i,ndofi_i,MPI_DOUBLE_PRECISION,comm_myplace1,comm_myisub,comm_all,request(ireq),ierr)

            ireq = ireq + 1
            call MPI_IRECV(rhoi_j,ndofi_j,MPI_DOUBLE_PRECISION,comm_myplace2,comm_myjsub,comm_all,request(ireq),ierr)
         end if
         do iinstr = 1,ninstructions
            owner = instructions(iinstr,1)
            isub  = instructions(iinstr,2)
            call get_index(isub,indexsub,lindexsub,isub_loc)

            call dd_get_interface_size(suba(isub_loc),ndofi,nnodi)

            lrhoi = ndofi
            allocate(rhoi(lrhoi))
            call dd_get_interface_diagonal(suba(isub_loc), rhoi,lrhoi)

            call MPI_SEND(rhoi,ndofi,MPI_DOUBLE_PRECISION,owner,isub,comm_all,ierr)

            deallocate(rhoi)
         end do
         nreq = ireq
         if (nreq.gt.0) then
            call MPI_WAITALL(nreq, request, statarray, ierr)
         end if
         ! debug
         !call info(routine_name,'All messages in pack 8 received on proc',myid)
         !if (i_compute_pair) then
         !   write(idmyunit,*)'rhoi_i'
         !   write(idmyunit,*) rhoi_i
         !   write(idmyunit,*)'rhoi_j'
         !   write(idmyunit,*) rhoi_j
         !   call flush(idmyunit)
         !end if

         ! get data about common interface
         ireq = 0
         if (i_compute_pair) then
            ! receive mapping of interface nodes into global nodes

            ireq = ireq + 1
            call MPI_IRECV(iingn_i,nnodi_i,MPI_INTEGER,comm_myplace1,comm_myisub,comm_all,request(ireq),ierr)

            ireq = ireq + 1
            call MPI_IRECV(iingn_j,nnodi_j,MPI_INTEGER,comm_myplace2,comm_myjsub,comm_all,request(ireq),ierr)
         end if
         ! send sizes of subdomains involved in problems
         do iinstr = 1,ninstructions
            owner = instructions(iinstr,1)
            isub  = instructions(iinstr,2)
            call get_index(isub,indexsub,lindexsub,isub_loc)

            call dd_get_interface_size(suba(isub_loc),ndofi,nnodi)

            liingn = nnodi
            allocate(iingn(liingn))
            call dd_get_interface_global_numbers(suba(isub_loc),iingn,liingn)

            call MPI_SEND(iingn,nnodi,MPI_INTEGER,owner,isub,comm_all,ierr)
            deallocate(iingn)
         end do
         nreq = ireq
         if (nreq.gt.0) then
            call MPI_WAITALL(nreq, request, statarray, ierr)
         end if
         ! debug
         !call info(routine_name,'All messages in pack 9 received on proc',myid)
         !if (i_compute_pair) then
         !   write(idmyunit,*)'iingn_i'
         !   write(idmyunit,*) iingn_i
         !   write(idmyunit,*)'iingn_j'
         !   write(idmyunit,*) iingn_j
         !end if

         ! compute trace of A
         if (i_compute_pair) then
            trA = sum(rhoi_i) + sum(rhoi_j)
         end if


         ! find common intersection of interface nodes
         if (i_compute_pair) then
            lcommon_interface = min(nnodi_i,nnodi_j)
            allocate(common_interface(lcommon_interface))
            call get_array_intersection(iingn_i,liingn_i, iingn_j,liingn_j,&
                                        common_interface,lcommon_interface,ncommon_interface)
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
            !write(idmyunit,*) 'Pair slavery:', pairslavery

            deallocate(rhoicomm)
            !write(idmyunit,*) 'weight',weight
         end if

         if (apply_null_projection) then
            ! get data about coarse basis functions to regularize the eigenproblem
            ireq = 0
            if (i_compute_pair) then
               ! receive global numbers of coarse nodes

               ireq = ireq + 1
               call MPI_IRECV(lphisi1_i,1,MPI_INTEGER,comm_myplace1,comm_myisub,comm_all,request(ireq),ierr)
               ireq = ireq + 1
               call MPI_IRECV(lphisi2_i,1,MPI_INTEGER,comm_myplace1,comm_myisub,comm_all,request(ireq),ierr)

               ireq = ireq + 1
               call MPI_IRECV(lphisi1_j,1,MPI_INTEGER,comm_myplace2,comm_myjsub,comm_all,request(ireq),ierr)
               ireq = ireq + 1
               call MPI_IRECV(lphisi2_j,1,MPI_INTEGER,comm_myplace2,comm_myjsub,comm_all,request(ireq),ierr)
            end if
            ! send sizes of subdomains involved in problems
            do iinstr = 1,ninstructions
               owner = instructions(iinstr,1)
               isub  = instructions(iinstr,2)
               call get_index(isub,indexsub,lindexsub,isub_loc)

               call dd_get_phisi_size(suba(isub_loc),lphisi1,lphisi2)

               call MPI_SEND(lphisi1,1,MPI_INTEGER,owner,isub,comm_all,ierr)
               call MPI_SEND(lphisi2,1,MPI_INTEGER,owner,isub,comm_all,ierr)
            end do
            nreq = ireq
            if (nreq.gt.0) then
               call MPI_WAITALL(nreq, request, statarray, ierr)
            end if
            ! debug
            !call info(routine_name,'All messages in pack 10 received on proc',myid)
            !if (i_compute_pair) then
            !   write(*,*) 'lphisi1_i'
            !   write(*,*)  lphisi1_i 
            !   write(*,*) 'lphisi2_i'
            !   write(*,*)  lphisi2_i 
            !   write(*,*) 'lphisi1_j'
            !   write(*,*)  lphisi1_j 
            !   write(*,*) 'lphisi2_j'
            !   write(*,*)  lphisi2_j 
            !   call flush(6)
            !end if

            ! get data for matrices phisi of subdomains in pair
            ireq = 0
            if (i_compute_pair) then

               ! prepare matrix Z as a superset of rigid body modes - based on
               ! coarse basis functions of the subdomains
               ! matriz Z = [ phisi_i     0    ]
               !            [   0      phisi_j ]
               allocate(phisi_i(lphisi1_i,lphisi2_i))
               allocate(phisi_j(lphisi1_j,lphisi2_j))

               ireq = ireq + 1
               call MPI_IRECV(phisi_i,lphisi1_i*lphisi2_i,MPI_DOUBLE_PRECISION,comm_myplace1,comm_myisub,&
                              comm_all,request(ireq),ierr)
               ireq = ireq + 1
               call MPI_IRECV(phisi_j,lphisi1_j*lphisi2_j,MPI_DOUBLE_PRECISION,comm_myplace2,comm_myjsub,&
                              comm_all,request(ireq),ierr)

            end if
            do iinstr = 1,ninstructions
               owner = instructions(iinstr,1)
               isub  = instructions(iinstr,2)
               call get_index(isub,indexsub,lindexsub,isub_loc)

               call dd_get_phisi_size(suba(isub_loc),lphisi1,lphisi2)

               allocate(phisi(lphisi1,lphisi2))
               call dd_get_phisi(suba(isub_loc),phisi, lphisi1,lphisi2)

               call MPI_SEND(phisi,lphisi1*lphisi2,MPI_DOUBLE_PRECISION,owner,isub,comm_all,ierr)
               deallocate(phisi)
            end do
            nreq = ireq
            if (nreq.gt.0) then
               call MPI_WAITALL(nreq, request, statarray, ierr)
            end if
            ! debug
            !call info(routine_name,'All messages in pack 11 received on proc',myid)
         end if
!-----profile
         if (profile) then
            call MPI_BARRIER(comm_all,ierr)
            call time_end(time)
            time_obtain_accu = time_obtain_accu + time
         end if
!-----profile

         ! get data for matrices phisi of subdomains in pair
         ! find null(B) as Z*null(Z'BZ), where
         ! null(B) subset of Range(Z)
!-----profile
         if (profile) then
            call MPI_BARRIER(comm_all,ierr)
            call time_start
         end if
!-----profile
         if (apply_null_projection) then
            if (i_compute_pair) then

               ! prepare matrix Z as a superset of rigid body modes - based on
               ! coarse basis functions of the subdomains
               ! matriz Z = [ phisi_i     0    ]
               !            [   0      phisi_j ]
               lz1 = lphisi1_i + lphisi1_j
               lz2 = lphisi2_i + lphisi2_j

               allocate(z(lz1,lz2))
               call zero(z,lz1,lz2)

               ! copy phisi_i into Z
               do j = 1,lphisi2_i
                  do i = 1,lphisi1_i
                     z(i,j) = phisi_i(i,j)
                  end do
               end do
               ! copy phisi_j into Z as block diagonal
               shifti = lphisi1_i
               shiftj = lphisi2_i
               do j = 1,lphisi2_j
                  do i = 1,lphisi1_j
                     z(shifti + i,shiftj + j) = phisi_j(i,j)
                  end do
               end do

               deallocate(phisi_i)
               deallocate(phisi_j)

               ! debug
               ! print dense matrix Z
               !call write_matrix(6,z,'e8.2')
            end if
         end if

         ! prepare space for eigenvectors
         ireq = 0
         ! prepare pointers to buffers with chunks of eigenvectors and their size
         llbufa = ninstructions
         allocate(lbufa(llbufa))
         do iinstr = 1,ninstructions
            owner = instructions(iinstr,1)
            isub  = instructions(iinstr,2)

            ireq = ireq + 1
            call MPI_IRECV(lbufa(iinstr),1,MPI_INTEGER,owner,isub,comm_all,request(ireq),ierr)
         end do
         if (i_compute_pair) then
            neigvec     = min(neigvecx,ndofcomm-ncommon_crows)
            if (debug) then
               if (neigvec.ne.neigvecx) then
                  write (*,*) 'ADAPTIVITY_SOLVE_EIGENVECTORS: Number of vectors reduced to',neigvec,' for pair',my_pair
               end if
            end if
            leigvec = neigvec * problemsize
            leigval = neigvec
            allocate(eigvec(leigvec),eigval(leigval))

            ! prepare space for buffers
            lbufsend_i = ndofi_i
            lbufsend_j = ndofi_j
            allocate(bufsend_i(lbufsend_i),bufsend_j(lbufsend_j))

            ! distribute sizes of chunks of eigenvectors
            call MPI_SEND(lbufsend_i,1,MPI_INTEGER,comm_myplace1,comm_myisub,comm_all,ierr)

            call MPI_SEND(lbufsend_j,1,MPI_INTEGER,comm_myplace2,comm_myjsub,comm_all,ierr)
         end if
         nreq = ireq
         if (nreq.gt.0) then
            call MPI_WAITALL(nreq, request, statarray, ierr)
         end if
         ! debug
         !call info(routine_name,'All messages in pack 12 received on proc',myid)

         ! prepare arrays kbufsend 
         lkbufsend = ninstructions
         allocate(kbufsend(lkbufsend))
         if (lkbufsend.gt.0) then
            kbufsend(1)  = 1
            do i = 2,ninstructions
               kbufsend(i) = kbufsend(i-1) + lbufa(i-1)
            end do
         end if
         lbufsend = 0
         do i = 1,ninstructions
            lbufsend = lbufsend + lbufa(i)
         end do
         lbufrecv = lbufsend
         allocate(bufrecv(lbufrecv),bufsend(lbufsend))

         ! set communicators
         comm_comm = comm_all
         comm_myid = myid

         if (use_explicit_schurs) then
            continue
         else
            if (i_compute_pair) then
               ! prepare auxiliary space used in each iteration of eigensolver
               comm_lxaux  = problemsize
               comm_lxaux2 = problemsize
               allocate(comm_xaux(comm_lxaux),comm_xaux2(comm_lxaux2))
            end if
         end if


!-----profile
         if (profile) then
            call MPI_BARRIER(comm_all,ierr)
            call time_end(time)
            time_null_prep_accu = time_null_prep_accu + time
         end if
!-----profile

!-----profile
         if (profile) then
            call MPI_BARRIER(comm_all,ierr)
            call time_start
         end if
!-----profile
         ! now construct matrix B*Z
         if (apply_null_projection) then
            if (use_explicit_schurs) then
               continue
            else
               if (i_compute_pair) then

                  lbz1 = lz1
                  lbz2 = lz2
                  allocate(bz(lbz1,lbz2))

                  ioper = 2 ! multiply by B
                  do j = 1,lz2
                     call adaptivity_mvecmult(suba,lsuba,indexsub,lindexsub,&
                                              problemsize,z(1,j),&
                                              problemsize,bz(1,j),problemsize,ioper)
                  end do
               end if
               call adaptivity_fake_lobpcg_driver
            end if
         end if
!-----profile
         if (profile) then
            call MPI_BARRIER(comm_all,ierr)
            call time_end(time)
            time_null_mult_accu = time_null_mult_accu + time
         end if
!-----profile

!-----profile
         if (profile) then
            call MPI_BARRIER(comm_all,ierr)
            call time_start
         end if
!-----profile
         if (apply_null_projection) then
            if (i_compute_pair) then
               ! debug
               ! print dense matrix BZ
               !call write_matrix(6,bz,'e8.2')
               ! construct matrix zbz = Z'*B*Z = Z'*BZ
               lzbz1 = lz2
               lzbz2 = lz2
               allocate(zbz(lzbz1,lzbz2))
               ! use BLAS for the multiply
               ! DGEMM - perform one of the matrix-matrix operations   C := alpha*op( A )*op( B ) + beta*C
               call DGEMM('T', 'N', lz2, lz2, lz1, 1.0_kr, z(:,:), lz1, bz(:,:), lbz1, 0._kr, zbz(:,:), lzbz1)
               !zbz = matmul(transpose(z),bz)
               ! debug
               ! print dense matrix ZBZ
               !call write_matrix(6,zbz,'e8.2')

               ! determine nullspace of ZBZ by eigenvalue decomposition
               lzbzeigval = lzbz1
               allocate(zbzeigval(lzbzeigval))
               ! determine size of work array
               lwork2 = 1
               allocate(work2(lwork2))
               lwork2 = -1
               ! first call routine just to find optimal size of WORK2
               call DSYEV( 'Vectors', 'Upper', lzbz1, zbz(:,:), lzbz1, zbzeigval, work2, lwork2, lapack_info )
               if (lapack_info.ne.0) then
                  call error(routine_name,'in LAPACK during finding size for nullspace eigenproblem solution:',lapack_info)
               end if
               lwork2 = int(work2(1))
               deallocate(work2)
               allocate(work2(lwork2))
               ! now call LAPACK to solve the eigenproblem
               call DSYEV( 'Vectors', 'Upper', lzbz1, zbz(:,:), lzbz1, zbzeigval, work2, lwork2, lapack_info )
               deallocate(work2)
               if (lapack_info.ne.0) then
                  call error(routine_name,'in LAPACK during solving nullspace eigenproblems:',lapack_info)
               end if

               ! debug
               ! print eigenvalues
               !write(6,*) 'eigenvalues by LAPACK:'
               !write(6,'(e8.2)') zbzeigval

               ! set tolerance - inspired by Octave null() function, relaxing epsilon a bit: 
               ! max (size (A)) * max (svd (A)) * eps
               null_tol = lzbz1 * zbzeigval(lzbzeigval) * 4 * epsilon(zbzeigval) 
               ! own way
               !null_tol = 1.e-10_kr * (zbzeigval(1) + zbzeigval(lzbzeigval)) / 2._kr 

               ! determine dimension of nullspace of Z'BZ
               null_dim = count(zbzeigval .lt. null_tol)

               ! if the largest eigenvalue is below tolerance, all values are in null_space
               if (zbzeigval(lzbzeigval).lt.numerical_zero) then
                  null_dim = lzbzeigval
               end if

               if (debug) then
                  call info(routine_name,'dimension of nullspace',null_dim)
               end if

               ! find null(B) = Z*null(Z'BZ) store it in BZ
               lnullB1 = lz1
               lnullB2 = null_dim
               ! leading dimension
               ldnullB = lnullB1
               allocate(nullB(lnullB1,lnullB2))
               call DGEMM('N', 'N', lz1, null_dim, lz2, 1.0_kr, z(:,:), lz1, zbz(:,1:null_dim), lzbz1, &
                          0._kr, nullB(:,:), ldnullB)
               
            ! Prepare projection onto complement of null(B) - orthogonalize null(B) -> get Q -> P = I-QQ'
               ! QR decomposition of matrix Z*null(Z'BZ) - stored in BZ 
               ! LAPACK arrays
               ltau3 = null_dim
               allocate(tau3(ltau3))
               lwork3 = max(1,null_dim)
               allocate(work3(lwork3))
               call DGEQRF( lnullB1, lnullB2, nullB(:,:), ldnullB, tau3, work3, lwork3, lapack_info)
               if (lapack_info.ne.0) then
                  call error(routine_name,' in LAPACK QR factorization of matrix null(B): ', lapack_info)
               end if
               ! in space of nullB are now stored factors R and Householder reflectors v
               is_nullB_ready = .true.

               deallocate(zbzeigval)
               deallocate(z)
               deallocate(bz)
               deallocate(zbz)
            end if
         end if
!-----profile
         if (profile) then
            call MPI_BARRIER(comm_all,ierr)
            call time_end(time)
            time_null_comp_accu = time_null_comp_accu + time
         end if
!-----profile
         ! debug
         ! check the nullspace
         !if (i_compute_pair) then
         !   ! prepare auxiliary space used in each iteration of eigensolver

         !   ioper = 2 ! multiply by B
         !   do j = 1,null_dim
         !      call adaptivity_mvecmult(suba,lsuba,indexsub,lindexsub,&
         !                               problemsize,bz(1,j),&
         !                               problemsize,z(1,j),problemsize,ioper)
         !   end do
         !end if
         !call adaptivity_fake_lobpcg_driver
         !if (i_compute_pair) then
         !   if (null_dim.gt.0) then
         !      call write_matrix(6,z(:,1:null_dim),'e8.2')
         !      call flush(6)
         !   end if
         !end if

         ! All arrays are ready for solving eigenproblems
         if (use_explicit_schurs) then
            ireq = 0
            if (i_compute_pair) then
               ! prepare space for local Schur complements
               lschur1_i = ndofi_i
               lschur2_i = ndofi_i
               allocate(schur_i(lschur1_i,lschur2_i))
               lschur1_j = ndofi_j
               lschur2_j = ndofi_j
               allocate(schur_j(lschur1_j,lschur2_j))

               ! receive Schur complements
               ireq = ireq + 1
               call MPI_IRECV(schur_i,lschur1_i*lschur2_i,MPI_DOUBLE_PRECISION,comm_myplace1,comm_myisub,comm_all,&
                              request(ireq),ierr)

               ireq = ireq + 1
               call MPI_IRECV(schur_j,lschur1_j*lschur2_j,MPI_DOUBLE_PRECISION,comm_myplace2,comm_myjsub,comm_all,&
                              request(ireq),ierr)
            end if
            ! send local Schur complements
            do iinstr = 1,ninstructions
               owner = instructions(iinstr,1)
               isub  = instructions(iinstr,2)
               call get_index(isub,indexsub,lindexsub,isub_loc)

               call dd_get_interface_size(suba(isub_loc),ndofi,nnodi)

               lschur1 = ndofi
               lschur2 = ndofi
               allocate(schur(lschur1,lschur2))
               call dd_get_schur(suba(isub_loc),schur,lschur1,lschur2)

               call MPI_SEND(schur,lschur1*lschur2,MPI_DOUBLE_PRECISION,owner,isub,comm_all,ierr)
               deallocate(schur)
            end do
            nreq = ireq
            if (nreq.gt.0) then
               call MPI_WAITALL(nreq, request, statarray, ierr)
            end if
            ! debug
            !call info(routine_name,'All messages in pack 2.1 received on proc',myid)

            if (i_compute_pair) then
               ! I have Schur complements, build matrices A and B for the
               ! generalized eigenvalue problem

               ! embed local Schur complements to matrices of the pair
               ldoubleschur1 = problemsize
               ldoubleschur2 = problemsize
               allocate(doubleschur(ldoubleschur1,ldoubleschur2))
               call zero(doubleschur,ldoubleschur1,ldoubleschur2)

               ! initialize the matrices with local Schur complements
               do j = 1,lschur2_i
                  do i = 1,lschur1_i
                     doubleschur(i,j) = schur_i(i,j)
                  end do
               end do
               shift = lschur1_i
               do j = 1,lschur2_j
                  do i = 1,lschur1_j
                     doubleschur(shift+i,shift+j) = schur_j(i,j)
                  end do
               end do
               deallocate(schur_i)
               deallocate(schur_j)

               ! TODO: build the matrix using matrix-matrix BLAS operations

               ! build matrices of generalized eigenproblem A and B
               lmata1 = problemsize
               lmata2 = problemsize
               allocate(mata(lmata1,lmata2))
               call zero(mata,lmata1,lmata2)
               lmatb1 = problemsize
               lmatb2 = problemsize
               allocate(matb(lmatb1,lmatb2))
               call zero(matb,lmatb1,lmatb2)

               lxaux = problemsize
               allocate(xaux(lxaux))
               lxaux2 = problemsize
               allocate(xaux2(lxaux2))

               ! build matrix A = P(I-RE)^T S (I-RE)P 
               do jcol = 1,problemsize
                  ! prepare auxiliary vector - a column of identity
                  call zero(xaux,lxaux)
                  xaux(jcol) = 1._kr

                  ! xaux = P xaux
                  call adaptivity_apply_null_projection(xaux,lxaux)

                  do i = 1,problemsize
                     xaux2(i) = xaux(i)
                  end do
                  ! apply (I-RE) xaux = (I-RR'D_P) xaux
                  !  - apply weights
                  call adaptivity_apply_weights(weight,lweight,xaux2,lxaux2)
                  !  - apply the R^T operator
                  call adaptivity_apply_RT(pairslavery,lpairslavery,xaux2,lxaux2)
                  !  - apply the R operator
                  call adaptivity_apply_R(pairslavery,lpairslavery,xaux2,lxaux2)
                  ! Ix - REx
                  do i = 1,problemsize
                     xaux(i) = xaux(i) - xaux2(i)
                  end do
                  ! apply double Schur
                  ! y = Sx
                  call DGEMV('N',problemsize,problemsize,1._kr,doubleschur,problemsize,xaux,1,0._kr,xaux2,1)
                  ! make a copy
                  do i = 1,problemsize
                     xaux(i) = xaux2(i)
                  end do
                  !  - apply the R^T operator
                  call adaptivity_apply_RT(pairslavery,lpairslavery,xaux2,lxaux2)
                  !  - apply the R operator
                  call adaptivity_apply_R(pairslavery,lpairslavery,xaux2,lxaux2)
                  !  - apply weights
                  call adaptivity_apply_weights(weight,lweight,xaux2,lxaux2)
                  ! Ix - E'R'x
                  do i = 1,problemsize
                     xaux(i) = xaux(i) - xaux2(i)
                  end do
                  ! xaux = P xaux
                  call adaptivity_apply_null_projection(xaux,lxaux)

                  ! make a copy
                  do i = 1,problemsize
                     mata(i,jcol) = xaux(i)
                  end do
               end do
               ! build matrix B = PSP + t*(I-P)
               ! find stabilization parameter
               stab = 0._kr 
               do i = 1,problemsize
                  if (abs(doubleschur(i,i)).gt.stab) then
                     stab = abs(doubleschur(i,i))
                  end if
               end do
               ! debug
               !print *,'stabilization parameter:',stab

               do jcol = 1,problemsize
                  ! prepare auxiliary vector - a column of identity
                  call zero(xaux,lxaux)
                  xaux(jcol) = 1._kr

                  ! xaux = P xaux
                  call adaptivity_apply_null_projection(xaux,lxaux)
                  ! add a note of this vector t*(I-P)
                  do i = 1,problemsize
                     matb(i,jcol) = -stab*xaux(i)
                  end do
                  matb(jcol,jcol) = matb(jcol,jcol) + stab

                  ! apply double Schur
                  ! y = Sx
                  call DGEMV('N',problemsize,problemsize,1._kr,doubleschur,problemsize,xaux,1,0._kr,xaux2,1)
                  ! make a copy
                  do i = 1,problemsize
                     xaux(i) = xaux2(i)
                  end do
                  ! xaux = P xaux
                  call adaptivity_apply_null_projection(xaux,lxaux)

                  ! make a copy
                  do i = 1,problemsize
                     matb(i,jcol) = matb(i,jcol) + xaux(i)
                  end do
               end do

               ! matrices A and B are ready, call LAPACK eigensolver
               leiglap = problemsize
               allocate(eiglap(leiglap))
               ! determine size of work array
               lwork2 = 1
               allocate(work2(lwork2))
               lwork2 = -1
               ! first call routine just to find optimal size of WORK2
               call DSYGV( 1,'V','U', problemsize, mata, problemsize, matb, problemsize, eiglap, work2,lwork2, lapack_info)
               if (lapack_info.ne.0) then
                  call error(routine_name,'in LAPACK during finding size for eigenproblem solution:',lapack_info)
               end if
               lwork2 = int(work2(1))
               deallocate(work2)
               print *,'I am here, LAPACK OK, lwork2:',lwork2
               allocate(work2(lwork2))
               ! now call LAPACK to solve the eigenproblem
               call DSYGV( 1,'V','U', problemsize, mata, problemsize, matb, problemsize, eiglap, work2,lwork2, lapack_info)
               deallocate(work2)
               print *,'I am here 2, LAPACK OK, lwork2:',lwork2
               if (lapack_info.ne.0) then
                  call error(routine_name,'in LAPACK during solving eigenproblems:',lapack_info)
               end if

               ! revert eigval to contain maximum
               !print *,'eiglap'
               !print '(f15.7)',eiglap
               !print *,'mata'
               !do i = 1,lmata1
               !   print '(f16.3)',mata(i,:)
               !end do
               !print *,'matb'
               !do i = 1,lmatb1
               !   print '(100f16.3)',matb(i,:)
               !end do
               !call flush(6)

               ! copy eigenvalues and eigenvectors
               do jcol = 1,neigvec
                  jrcol = problemsize - jcol + 1

                  ! copy eigenvalue
                  eigval(jcol) = eiglap(jrcol)

                  ! copy eigenvector
                  do i = 1,problemsize
                     eigvec((jcol-1)*problemsize + i) = mata(i,jrcol)
                  end do
               end do
               deallocate(eiglap)
               deallocate(mata)
               deallocate(matb)
            end if
         else

!-----profile
            if (profile) then
               call MPI_BARRIER(comm_all,ierr)
               call time_start
            end if
!-----profile
            if (lobpcg_preconditioner .eq. 1) then
               ! set up local BDDC preconditioner for LOBPCG from existing data

               ireq = 0
               if (i_compute_pair) then

                  ! check type of matrix - for LOBPCG, any other than SPD is forbidden
                  if (matrixtype.ne.1) then
                     call error(routine_name,'adaptivity without explicit Schur complements not supported for non SPD matrices')
                  end if

                  ! prepare space for local coarse matrices of slaves
                  ! since LOBPCG does not work for indefinite nor general matrices, symmetric storage is assumed
                  lcoarsem_i = ((lindrowc_i+1) * lindrowc_i) /2
                  lcoarsem_j = ((lindrowc_j+1) * lindrowc_j) /2

                  allocate(coarsem_i(lcoarsem_i))
                  allocate(coarsem_j(lcoarsem_j))

                  ! distribute sizes of common_crows (with global indices of rows of C_i and C_j)
                  ireq = ireq + 1
                  call MPI_IRECV(coarsem_i,lcoarsem_i,MPI_DOUBLE_PRECISION,comm_myplace1,comm_myisub,comm_all,&
                                 request(ireq),ierr)
                  ireq = ireq + 1
                  call MPI_IRECV(coarsem_j,lcoarsem_j,MPI_DOUBLE_PRECISION,comm_myplace2,comm_myjsub,comm_all,&
                                 request(ireq),ierr)
               end if
               nreq = ireq

               do iinstr = 1,ninstructions
                  owner = instructions(iinstr,1)
                  isub  = instructions(iinstr,2)
                  
                  call get_index(isub,indexsub,lindexsub,isub_loc)

                  ! extract the local coarse matrix and send it to master of pair
                  call dd_get_coarsem_size(suba(isub_loc),lcoarsem)

                  allocate(coarsem(lcoarsem))
                  call dd_get_coarsem(suba(isub_loc),coarsem,lcoarsem)

                  ! send the matrix
                  call MPI_SEND(coarsem,lcoarsem,MPI_DOUBLE_PRECISION,owner,isub,comm_all,ierr)

                  deallocate(coarsem)
               end do

               ! wait for local coarse matrices
               if (nreq.gt.0) then
                  call MPI_WAITALL(nreq, request, statarray, ierr)
               end if
               ! debug
               !call info(routine_name,'All messages in pack 30.3 received on proc',myid)

               if (i_compute_pair) then

                  ! assemble global coarse matrix 
                  ! determine length
                  lcoarsem_adapt1 = lindrowc_i + lindrowc_j - ncommon_crows
                  lcoarsem_adapt2 = lcoarsem_adapt1
                  allocate(coarsem_adapt(lcoarsem_adapt1,lcoarsem_adapt2))
                  call zero(coarsem_adapt,lcoarsem_adapt1,lcoarsem_adapt2)

                  ! construct embedding of local coarse matrices to coarse matrix of the pair
                  ! order unknowns as:  unique_i   |   unique_j   | common_crows
                  lcoarsem_embed = lcoarsem_adapt1
                  allocate(coarsem_embed(lcoarsem_embed))

                  ! debug
                  !print *, 'indrowc_i',indrowc_i
                  !print *, 'indrowc_j',indrowc_j
                  !print *, 'ncommon_crows',ncommon_crows

                  i_loc = 0
                  ! global coarse dof unique to i
                  do i = 1,lindrowc_i
                     indg = indrowc_i(i)

                     if (.not. any(common_crows(1:ncommon_crows) .eq. indg)) then
                        i_loc = i_loc + 1

                        coarsem_embed(i_loc) = indg
                     end if
                  end do
                  ! global coarse dof unique to j
                  do i = 1,lindrowc_j
                     indg = indrowc_j(i)

                     if (.not. any(common_crows(1:ncommon_crows) .eq. indg)) then
                        i_loc = i_loc + 1

                        coarsem_embed(i_loc) = indg
                     end if
                  end do
                  ! global coarse dof common to i and j
                  do i = 1,ncommon_crows
                     indg = common_crows(i)

                     i_loc = i_loc + 1
                     coarsem_embed(i_loc) = indg
                  end do
                  ! debug
                  !print *, 'coarsem_embed',coarsem_embed

                  ! prepare local versions of indrowc_ij arrays
                  lindrowc_adapt_i = lindrowc_i
                  lindrowc_adapt_j = lindrowc_j
                  allocate(indrowc_adapt_i(lindrowc_adapt_i))
                  allocate(indrowc_adapt_j(lindrowc_adapt_j))

                  do i = 1,lindrowc_i
                     indg = indrowc_i(i)
                     
                     call get_index(indg,coarsem_embed,lcoarsem_embed,indc_loc)
                     if (indc_loc .eq. -1) then
                        call error(routine_name,'Index of coarse dof not found.', indg)
                     end if

                     indrowc_adapt_i(i) = indc_loc
                  end do
                  do i = 1,lindrowc_j
                     indg = indrowc_j(i)
                     
                     call get_index(indg,coarsem_embed,lcoarsem_embed,indc_loc)
                     if (indc_loc .eq. -1) then
                        call error(routine_name,'Index of coarse dof not found.', indg)
                     end if

                     indrowc_adapt_j(i) = indc_loc
                  end do
                  ! debug
                  !print *, 'indrowc_adapt_i',indrowc_adapt_i
                  !print *, 'indrowc_adapt_j',indrowc_adapt_j

                  ! assemble matrices to coarsem_adapt while storing them in upper triangle of new 2D array
                  ! K_C = [ K_C_i  0     x        ] with overlapping common global coarse dof 
                  !       [  0    K_C_j  x        ]
                  !       [  x     x   K_C_common ]

                  ! add K_C_i
                  icoarsem = 0
                  do j = 1,lindrowc_i
                     jcol = indrowc_adapt_i(j)
                     do i = 1,j
                        irow = indrowc_adapt_i(i)

                        icoarsem = icoarsem + 1

                        if (irow.gt.jcol) then
                           i_upp = jcol
                           j_upp = irow
                        else
                           i_upp = irow
                           j_upp = jcol
                        end if

                        coarsem_adapt(i_upp,j_upp) = coarsem_adapt(i_upp,j_upp) + coarsem_i(icoarsem)
                     end do
                  end do
                  ! add K_C_j
                  icoarsem = 0
                  do j = 1,lindrowc_j
                     jcol = indrowc_adapt_j(j)
                     do i = 1,j
                        irow = indrowc_adapt_j(i)

                        icoarsem = icoarsem + 1

                        if (irow.gt.jcol) then
                           i_upp = jcol
                           j_upp = irow
                        else
                           i_upp = irow
                           j_upp = jcol
                        end if

                        coarsem_adapt(i_upp,j_upp) = coarsem_adapt(i_upp,j_upp) + coarsem_j(icoarsem)
                     end do
                  end do

                  deallocate(coarsem_embed)
                  deallocate(coarsem_i)
                  deallocate(coarsem_j)

                  ! debug
                  ! print matrix
                  !call write_matrix(6,coarsem_adapt,'e8.2')
                  !call flush(6)

                  ! prepare eigendecomposition of the coarse matrix of the pair by LAPACK
                  ! eigenvalue decomposition of local coarse matrix K_C_loc
                  lkceigval = lcoarsem_adapt1
                  allocate(kceigval(lkceigval))

                  ! determine size of work array
                  lwork2 = 1
                  allocate(work2(lwork2))
                  lwork2 = -1
                  ! first call routine just to find optimal size of WORK2
                  call DSYEV( 'Vectors', 'Upper', lcoarsem_adapt1, coarsem_adapt(:,:), lcoarsem_adapt1, &
                              kceigval, work2, lwork2, lapack_info )
                  if (lapack_info.ne.0) then
                     call error(routine_name,'in LAPACK during finding size for eigendecomposition of local coarse matrix:',&
                                lapack_info)
                  end if
                  lwork2 = int(work2(1))
                  deallocate(work2)
                  allocate(work2(lwork2))
                  ! now call LAPACK to solve the eigenproblem
                  call DSYEV( 'Vectors', 'Upper', lcoarsem_adapt1, coarsem_adapt(:,:), lcoarsem_adapt1, &
                              kceigval, work2, lwork2, lapack_info )
                  deallocate(work2)
                  if (lapack_info.ne.0) then
                     call error(routine_name,'in LAPACK during eigendecomposition of local coarse matrix:',lapack_info)
                  end if

                  ! debug
                  ! print eigenvalues
                  !write(6,*) 'eigenvalues by LAPACK:'
                  !write(6,'(e8.2)') kceigval

                  ! set tolerance - inspired by Octave null() function, relaxing epsilon a bit: 
                  ! max (size (A)) * max (svd (A)) * eps
                  null_tol = lcoarsem_adapt1 * kceigval(lkceigval) * 4 * epsilon(kceigval) 
                  ! own way
                  !null_tol = 1.e-10_kr * (zbzeigval(1) + zbzeigval(lzbzeigval)) / 2._kr 

                  ! determine dimension of nullspace of K_C
                  null_dim = count(kceigval .lt. null_tol)
                  if (debug) then
                     call warning(routine_name,'nullspace dimension of local coarse matrix',null_dim)
                  end if

                  ! compute inverse of Lambda for nonsingular eigenvalues
                  do i = 1,null_dim
                     kceigval(i) = 0._kr
                  end do
                  do i = null_dim + 1,lcoarsem_adapt1
                     kceigval(i) = 1._kr / kceigval(i)
                  end do

                  ! prepare memory for iterations
                  comm_lresc = lcoarsem_adapt1
                  allocate(comm_resc(comm_lresc))
                  allocate(comm_resc_i(comm_lresc))
                  allocate(comm_resc_j(comm_lresc))


                  if (debug) then
                     call info(routine_name,'local coarse matrix factorized for pair ',my_pair)
                  end if

               end if

               ! end set-up local BDDC preconditioner
            end if
!-----profile
            if (profile) then
               call MPI_BARRIER(comm_all,ierr)
               call time_end(time)
               time_precond_accu = time_precond_accu + time
            end if
!-----profile


!-----profile
            if (profile) then
               call MPI_BARRIER(comm_all,ierr)
               call time_start
            end if
!-----profile
            ! call LOBCPG
            if (i_compute_pair) then
               if (use_vec_values .eq. 1) then

                  ! Initialize the array with own random number generator
                  !   reinitialize the random number sequence for each pair for reproducibility 
                  !   (LOBPCG is sensitive for initial guess)
                  call initialize_random_number
                  do i = 1,neigvec*problemsize
                     call get_random_number(eigvec(i))
                  end do

                  ! debug
                  ! read initial guess of vectors
                  !write(*,*) 'neigvec =',neigvec,'problemsize =',problemsize
                  !open(unit = idmyunit,file='start_guess.txt')
                  !do i = 1,problemsize
                  !   read(idmyunit,*) (eigvec((j-1)*problemsize + i),j = 1,neigvec)
                  !end do
                  !close(idmyunit)
               end if

               ! set name of individual file for glob
               call getfname('pair',my_pair,'EIG',filename)
                  
               if (recompute_vectors) then
                  ! set LOBPCG absolute tolerance based on some values of diagonal of A
                  lobpcg_tol = lobpcg_rel_tol * trA / problemsize

                  ! compute eigenvectors and store them into file
                  if (debug) then
                     write(*,*) 'myid =',myid,', I am calling eigensolver for pair ',my_pair
                     write(*,*) 'myid =',myid,', LOBPCG tolerance: ',lobpcg_tol
                  end if
                  call lobpcg_driver(problemsize,neigvec,lobpcg_tol,lobpcg_maxit,lobpcg_verbosity,use_vec_values,&
                                     lobpcg_preconditioner,&
                                     eigval,eigvec,lobpcg_iter,ierr)
                  if (ierr.ne.0) then
                     call warning(routine_name,'LOBPCG exited with nonzero code for pair',my_pair)
                     if (debug) then
                        call warning(routine_name,'LOBPCG error code: ',ierr)
                     end if
                     ! if problems with preconditioner are detected, try it without preconditioner
                     if (ierr.eq.-1 .and. lobpcg_preconditioner .ne. 0 .and. try_harder) then
                        call warning(routine_name,'trying without preconditioner for pair',my_pair)
                        no_prec = 0
                        call lobpcg_driver(problemsize,neigvec,lobpcg_tol,lobpcg_maxit,lobpcg_verbosity,use_vec_values,&
                                           no_prec,&
                                           eigval,eigvec,lobpcg_iter,ierr)
                     end if
                  end if
                  write(*,*) 'myid =',myid,', LOBPCG finished in ',lobpcg_iter,' iterations for pair ',my_pair

                  ! turn around the eigenvalues to be the largest
                  eigval = -eigval

                  if (debug) then
                     open(unit = idmyunit,file=filename,status='replace',form='formatted')
                     write(idmyunit,*) neigvec, problemsize
                     do i = 1,neigvec
                        write(idmyunit,*) eigval(i)
                     end do
                     do i = 1,problemsize
                        write(idmyunit,*) (eigvec((j-1)*problemsize + i),j = 1,neigvec)
                     end do
                     close(idmyunit)
                     write(*,*) 'myid =',myid,', Eigenvalues stored in file ',trim(filename)
                  end if
                  !write(90+myid,*) 'x ='
                  !do i = 1,problemsize
                  !   write(90+myid,'(30f13.7)') (eigvec((j-1)*problemsize + i),j = 1,neigvec)
                  !end do
                  !do i = 1,problemsize
                  !   write(88, '(30e15.5)'), (eigvec((j-1)*problemsize + i),j = 1,neigvec)
                  !end do
               else
                  ! read eigenvectors from file
                  if (debug) then
                     write(*,*) 'myid =',myid,', I am reading eigenvalues for pair ',my_pair,' from file ',trim(filename)
                     open(unit = idmyunit,file=filename,status='old',form='formatted')
                     read(idmyunit,*) neigvecf, problemsizef
                     if (neigvecf.ne.neigvec .or. problemsizef.ne.problemsize) then
                        call error(routine_name,'Wrong dimensions of input eigenvectors for glob', gglob)
                     end if
                     do i = 1,neigvec
                        read(idmyunit,*) eigval(i)
                     end do
                     do i = 1,problemsize
                        read(idmyunit,*) (eigvec((j-1)*problemsize + i),j = 1,neigvec)
                     end do
                     close(idmyunit)
                  end if
               end if
            end if
            call adaptivity_fake_lobpcg_driver
!-----profile
            if (profile) then
               call MPI_BARRIER(comm_all,ierr)
               call time_end(time)
               time_solve_accu = time_solve_accu + time
            end if
!-----profile
         end if


!-----profile
         if (profile) then
            call MPI_BARRIER(comm_all,ierr)
            call time_start
         end if
!-----profile
         ! select eigenvectors of eigenvalues exceeding threshold
         est_loc = 0._kr
         if (i_compute_pair) then

            if (debug) then
               write(*,*) 'eigval for pair:',my_pair,' between subdomains ',comm_myisub,' and',comm_myjsub
               write(*,'(f20.10)') eigval
            end if

            nadaptive = count(eigval.ge.threshold_eigval)

            if (debug) then
               write(*,*) 'ADAPTIVITY_SOLVE_EIGENVECTORS: I am going to add ',nadaptive,' constraints for pair ',my_pair
            end if

            ! find estimator of condition number
            if (nadaptive.lt.neigvec) then
               est_loc = eigval(nadaptive + 1)
            else
               est_loc = eigval(nadaptive)
            end if

            lconstraints1 = problemsize
            lconstraints2 = nadaptive
            allocate(constraints(lconstraints1,lconstraints2))
         end if

         if (use_explicit_schurs) then
            if (i_compute_pair) then
               ! multiply eigenvectors locally
               do jcol = 1,nadaptive

                  ! copy eigvec to xaux
                  do i = 1,problemsize
                     xaux(i) = eigvec((jcol-1)*problemsize + i)
                  end do

                  ! xaux = P xaux
                  call adaptivity_apply_null_projection(xaux,lxaux)

                  do i = 1,problemsize
                     xaux2(i) = xaux(i)
                  end do
                  ! apply (I-RE) xaux = (I-RR'D_P) xaux
                  !  - apply weights
                  call adaptivity_apply_weights(weight,lweight,xaux2,lxaux2)
                  !  - apply the R^T operator
                  call adaptivity_apply_RT(pairslavery,lpairslavery,xaux2,lxaux2)
                  !  - apply the R operator
                  call adaptivity_apply_R(pairslavery,lpairslavery,xaux2,lxaux2)
                  ! Ix - REx
                  do i = 1,problemsize
                     xaux(i) = xaux(i) - xaux2(i)
                  end do
                  ! apply double Schur
                  ! y = Sx
                  call DGEMV('N',problemsize,problemsize,1._kr,doubleschur,problemsize,xaux,1,0._kr,xaux2,1)
                  ! make a copy
                  do i = 1,problemsize
                     xaux(i) = xaux2(i)
                  end do
                  !  - apply the R^T operator
                  call adaptivity_apply_RT(pairslavery,lpairslavery,xaux2,lxaux2)
                  !  - apply the R operator
                  call adaptivity_apply_R(pairslavery,lpairslavery,xaux2,lxaux2)
                  !  - apply weights
                  call adaptivity_apply_weights(weight,lweight,xaux2,lxaux2)
                  ! Ix - E'R'x
                  do i = 1,problemsize
                     xaux(i) = xaux(i) - xaux2(i)
                  end do
                  ! xaux = P xaux
                  call adaptivity_apply_null_projection(xaux,lxaux)

                  ! make a copy
                  do i = 1,problemsize
                     constraints(i,jcol) = xaux(i)
                  end do
               end do

               deallocate(xaux)
               deallocate(xaux2)
               deallocate(doubleschur)
            end if
         else
            if (i_compute_pair) then
               ! construct constraints out of these eigenvectors
               ioper = 1 ! multiply by A
               do j = 1,nadaptive
                  call adaptivity_mvecmult(suba,lsuba,indexsub,lindexsub,&
                                           problemsize,eigvec((j-1)*problemsize + 1),&
                                           problemsize,constraints(1,j),problemsize,ioper)
               end do
               ! clean auxiliary arrays
               deallocate(comm_xaux,comm_xaux2)
            end if
            call adaptivity_fake_lobpcg_driver
         end if

         ! global communication
         call MPI_ALLREDUCE(est_loc,est_round,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm_all,ierr)
         !if (my_pair.ge.0) then
         !   write(*,*) 'Constraints to be added on pair ',my_pair
         !   do i = 1,problemsize
         !      write(*,'(100f15.3)') (constraints(i,j),j = 1,nadaptive)
         !   end do
         !end if

         ! REALLOCATE BUFFERS
         deallocate(bufrecv,bufsend)
         if (i_compute_pair) then
            deallocate(bufsend_i,bufsend_j)
         end if

         ! distribute adaptive constraints to slaves
         ! prepare space for these constraints
         ireq = 0
         ! prepare pointers to buffers with chunks of eigenvectors and their size
         do iinstr = 1,ninstructions
            owner = instructions(iinstr,1)
            isub  = instructions(iinstr,2)

            ireq = ireq + 1
            call MPI_IRECV(lbufa(iinstr),1,MPI_INTEGER,owner,isub,comm_all,request(ireq),ierr)
         end do
         if (i_compute_pair) then
            ! prepare space for buffers
            lbufsend_i = ndofi_i * nadaptive
            lbufsend_j = ndofi_j * nadaptive
            allocate(bufsend_i(lbufsend_i),bufsend_j(lbufsend_j))

            ! distribute sizes of chunks of eigenvectors
            call MPI_SEND(lbufsend_i,1,MPI_INTEGER,comm_myplace1,comm_myisub,comm_all,ierr)
            call MPI_SEND(lbufsend_j,1,MPI_INTEGER,comm_myplace2,comm_myjsub,comm_all,ierr)
         end if
         nreq = ireq
         if (nreq.gt.0) then
            call MPI_WAITALL(nreq, request, statarray, ierr)
         end if
         ! debug
         !call info(routine_name,'All messages in pack 13 received on proc',myid)

         ! prepare arrays kbufsend and 
         if (lkbufsend .gt. 0) then
            kbufsend(1)  = 1
            do i = 2,ninstructions
               kbufsend(i) = kbufsend(i-1) + lbufa(i-1)
            end do
         end if
         lbufsend = 0
         do i = 1,ninstructions
            lbufsend = lbufsend + lbufa(i)
         end do
         lbufrecv = lbufsend
         allocate(bufrecv(lbufrecv),bufsend(lbufsend))

         ! distribute adaptive constraints to slaves
         ireq = 0
         do iinstr = 1,ninstructions
            owner = instructions(iinstr,1)
            isub  = instructions(iinstr,2)
            gglob = instructions(iinstr,3)
            call get_index(isub,indexsub,lindexsub,isub_loc)

            call dd_get_interface_size(suba(isub_loc),ndofi,nnodi)

            if (lbufa(iinstr) .gt. 0) then
               ireq = ireq + 1
               call MPI_IRECV(bufrecv(kbufsend(iinstr)),lbufa(iinstr),MPI_DOUBLE_PRECISION,owner,isub,comm_all,request(ireq),ierr)
            end if
         end do
         if (i_compute_pair) then
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
               call MPI_SEND(bufsend_i,lbufsend_i,MPI_DOUBLE_PRECISION,comm_myplace1,comm_myisub,comm_all,ierr)
            end if

            if (lbufsend_j .gt. 0) then
               call MPI_SEND(bufsend_j,lbufsend_j,MPI_DOUBLE_PRECISION,comm_myplace2,comm_myjsub,comm_all,ierr)
            end if
         end if
         nreq = ireq
         if (nreq.gt.0) then
            call MPI_WAITALL(nreq, request, statarray, ierr)
         end if
         ! debug
         !call info(routine_name,'All messages in pack 14 received on proc',myid)

         ! processors now own data of adaptively found constraints on their subdomains, they have to filter globs and load them into the structure
         ! gather back data about really loaded constraints
         ireq = 0
         if (i_compute_pair) then
            ! receive mapping of interface nodes into global nodes

            ireq = ireq + 1
            call MPI_IRECV(nvalid_i,1,MPI_INTEGER,comm_myplace1,comm_myisub,comm_all,request(ireq),ierr)

            ireq = ireq + 1
            call MPI_IRECV(nvalid_j,1,MPI_INTEGER,comm_myplace2,comm_myjsub,comm_all,request(ireq),ierr)
         end if
         do iinstr = 1,ninstructions
            owner = instructions(iinstr,1)
            isub  = instructions(iinstr,2)
            gglob = instructions(iinstr,3)
            call get_index(isub,indexsub,lindexsub,isub_loc)

            call dd_get_interface_size(suba(isub_loc),ndofi,nnodi)

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
            call dd_load_adaptive_constraints(suba(isub_loc),gglob,cadapt,lcadapt1,lcadapt2, nvalid)

            call MPI_SEND(nvalid,1,MPI_INTEGER,owner,isub,comm_all,ierr)

            deallocate(cadapt)
         end do
         nreq = ireq
         if (nreq.gt.0) then
            call MPI_WAITALL(nreq, request, statarray, ierr)
         end if
         ! debug
         !call info(routine_name,'All messages in pack 15 received on proc',myid)

         if (i_compute_pair) then
            ! check that number of loaded constraints left and write match
            if (nvalid_i .ne. nvalid_j) then
               call error(routine_name,'Number of valid constraints does not match for pair ',my_pair)
            end if

         end if

!-----profile
         if (profile) then
            call MPI_BARRIER(comm_all,ierr)
            call time_end(time)
            time_postp_accu = time_postp_accu + time
         end if
!-----profile
         deallocate(lbufa)
         deallocate(kbufsend) 
         deallocate(bufrecv,bufsend)
         if (i_compute_pair) then
            deallocate(constraints)
            deallocate(eigvec,eigval)
            deallocate(bufsend_i,bufsend_j)
            deallocate(weight)
            deallocate(pairslavery)
            deallocate(common_interface)
            deallocate(iingn_i,iingn_j)
            deallocate(rhoi_i,rhoi_j)
            deallocate(work)
            deallocate(tau)
            if (apply_null_projection) then
               deallocate(work3)
               deallocate(tau3)
               deallocate(nullB)
               is_nullB_ready = .false.
            end if
            if (lobpcg_preconditioner .eq. 1) then
               deallocate(indrowc_adapt_i)
               deallocate(indrowc_adapt_j)
               deallocate(kceigval)
               deallocate(coarsem_adapt)
               deallocate(comm_resc)
               deallocate(comm_resc_i)
               deallocate(comm_resc_j)
            end if
            deallocate(kdofi_i,kdofi_j)
            deallocate(dij)
            deallocate(common_crows)
            deallocate(nndfi_i,nndfi_j)
            deallocate(i_c_sparse_i,i_c_sparse_j)
            deallocate(j_c_sparse_i,j_c_sparse_j)
            deallocate(c_sparse_i,c_sparse_j)
            deallocate(indrowc_i,indrowc_j)
         end if

         ! find estimate of condition number
         est = max(est,est_round)

         if (debug) then
            if (myid.eq.0) then
               call info(routine_name, 'completed round',iround)
               call info(routine_name, 'number of calls to matrix multiply',comm_calls)
            end if
         end if
      end do ! loop over rounds of eigenvalue pairs

      if (profile) then
         if (myid.eq.0) then
            call time_print('obtaining data',time_obtain_accu)
            call time_print('computing preparations in nullspace',time_null_prep_accu)
            call time_print('computing multiplies in nullspace',time_null_mult_accu)
            call time_print('computing nullspace by eig decomposition',time_null_comp_accu)
            call time_print('preconditioner setup',time_precond_accu)
            call time_print('eigenproblems solution',time_solve_accu)
            call time_print('postprocessing constraints',time_postp_accu)
         end if
      end if

      if (debug) then
         if (myid.eq.0) then
            write(*,*) 'Expected estimated condition number: ',est
         end if
      end if

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
      integer :: n = 0, lx = 0, ly = 0
      real(kr) :: x, y


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

!*******************************************************************************
subroutine adaptivity_mvecmult(suba,lsuba,indexsub,lindexsub,n,x,lx,y,ly,idoper)
!*******************************************************************************
! realizes the multiplications with the strange matrices in the local eigenvalue problems for adaptive BDDC

use module_dd
use module_utils
implicit none
include "mpif.h"

real(kr),external :: ddot

! array of sub structure
      integer,intent(in) ::                lsuba
      type(subdomain_type),intent(inout) :: suba(lsuba)
! global indices of local subdomains
      integer,intent(in) :: lindexsub
      integer,intent(in) ::  indexsub(lindexsub)

! length of vector
integer,intent(in) ::   n 
integer,intent(in) ::  lx 
real(kr),intent(in) ::  x(lx)
integer,intent(in) ::  ly 
real(kr),intent(out) :: y(ly)
! determine which matrix should be multiplied
!  1 - A = P(I-RE)'S(I-RE)P
!  2 - B = P_barPSPP_bar
!  3 - not called from LOBPCG, called from fake looper - just perform demanded multiplications by S
!  5 - preconditioning by local BDDC
!  -3 - set on exit if iterational process should be stopped now
integer,intent(inout) ::  idoper 

! local vars
character(*),parameter:: routine_name = 'ADAPTIVITY_MVECMULT'
integer :: i
integer :: iinstr, isub, isub_loc, owner, point1, point2, point, length, do_i_compute, is_active

! small BDDC related vars
integer ::             laux2
real(kr),allocatable :: aux2(:)
integer ::             lrescs
real(kr),allocatable :: rescs(:)
integer ::             lsolis
real(kr),allocatable :: solis(:)
integer ::             nrhs, nnods, nelems, ndofs, ndofaaugs, lindrowc, ndofi, nnodi
logical :: transposed
integer :: stat(MPI_STATUS_SIZE)

integer :: indc

logical :: i_am_slave, all_slaves

logical :: suppress_preconditioning = .false.

integer :: ireq, nreq, ierr

comm_calls = comm_calls + 1


! debug
!print *,'myid = ',comm_myid,'idoper = ',idoper
!call flush(6)

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
if (i_compute_pair) then
   if (idoper.eq.3) then
      ! called from outside of LOBPCG
      do_i_compute = 0
   else if (idoper.eq.5) then
      ! use as preconditioner
      do_i_compute = 2
   else
      ! use as matrix multiply
      do_i_compute = 1
   end if

   call MPI_SEND(do_i_compute,1,MPI_INTEGER,comm_myplace1,comm_myisub,comm_comm,ierr)
   call MPI_SEND(do_i_compute,1,MPI_INTEGER,comm_myplace2,comm_myjsub,comm_comm,ierr)
end if
nreq = ireq
if (nreq.gt.0) then
   call MPI_WAITALL(nreq, request, statarray, ierr)
end if

ireq = 0
! What follows is performed independently on from where I was called
! obtain chunks of eigenvectors to multiply them
do iinstr = 1,ninstructions
   owner     = instructions(iinstr,1)
   isub      = instructions(iinstr,2)
   is_active = instructions(iinstr,4)

   if (is_active .gt. 0) then
      ireq = ireq + 1
      call MPI_IRECV(bufrecv(kbufsend(iinstr)),lbufa(iinstr),MPI_DOUBLE_PRECISION,owner,isub,comm_comm,request(ireq),ierr)
   end if
end do
! common operations for A and B - checks and projection to null(D_ij)
if (idoper.eq.1 .or. idoper.eq.2 .or. idoper.eq.5 ) then
   ! check the dimensions
   if (n .ne. problemsize) then
       call error(routine_name, 'Vector size mismatch.', n)
   end if
   if (lx .ne. ly) then
       call error(routine_name, 'Data size mismatch: lx =',lx)
   end if
   if (lx .ne. n) then
       call error(routine_name, 'Data size mismatch. lx = ', lx)
   end if
   ! are arrays allocated ?
   if (.not.allocated(comm_xaux) .or. .not.allocated(comm_xaux2)) then
       call error(routine_name, 'Auxiliary arrays not allocated. idoper =',idoper)
   end if
   ! correct size ?
   if (comm_lxaux .ne. comm_lxaux2) then
       call error(routine_name, 'Array size not match.',comm_lxaux2)
   end if
   ! correct size ?
   if (comm_lxaux .ne. lx) then
       call error(routine_name, 'Array size not match.',comm_lxaux)
   end if


   ! make temporary copy
   do i = 1,lx
      comm_xaux(i) = x(i)
   end do
   
   ! debug
   !if (idoper.eq.5) then
   !   print *, 'comm_xaux before'
   !   print '(f9.3)',comm_xaux
   !end if

   ! xaux = P_bar xaux
   if (idoper.eq.2) then
      if (is_nullB_ready) then
         call adaptivity_apply_complementary_projection(comm_xaux,comm_lxaux)
      end if
   end if

   ! apply projection Pi - to fulfill common constraints
   ! xaux = P xaux (in all regimes)
   call adaptivity_apply_null_projection(comm_xaux,comm_lxaux)

   if (idoper.eq.1) then
      ! make temporary copy
      do i = 1,lx
         comm_xaux2(i) = comm_xaux(i)
      end do
      ! apply (I-RE) xaux = (I-RR'D_P) xaux
      !  - apply weights
      call adaptivity_apply_weights(weight,lweight,comm_xaux2,comm_lxaux2)
      !  - apply the R^T operator
      call adaptivity_apply_RT(pairslavery,lpairslavery,comm_xaux2,comm_lxaux2)
      !  - apply the R operator
      call adaptivity_apply_R(pairslavery,lpairslavery,comm_xaux2,comm_lxaux2)
      ! Ix - REx
      do i = 1,problemsize
         comm_xaux(i) = comm_xaux(i) - comm_xaux2(i)
      end do
   end if

   ! distribute sizes of chunks of eigenvectors
   point1 = 1
   call MPI_SEND(comm_xaux(point1),ndofi_i,MPI_DOUBLE_PRECISION,comm_myplace1,comm_myisub,comm_comm,ierr)

   point2 = ndofi_i + 1
   call MPI_SEND(comm_xaux(point2),ndofi_j,MPI_DOUBLE_PRECISION,comm_myplace2,comm_myjsub,comm_comm,ierr)
end if
nreq = ireq
if (nreq.gt.0) then
   call MPI_WAITALL(nreq, request, statarray, ierr)
end if

! If I am preconditioning, get coarse residuals
ireq = 0
! Continue only of I was called from LOBPCG
if (idoper.eq.5) then
   ireq = ireq + 1
   call MPI_IRECV(comm_resc_i,comm_lresc,MPI_DOUBLE_PRECISION,comm_myplace1,comm_myisub,comm_comm,request(ireq),ierr)

   ireq = ireq + 1
   call MPI_IRECV(comm_resc_j,comm_lresc,MPI_DOUBLE_PRECISION,comm_myplace2,comm_myjsub,comm_comm,request(ireq),ierr)
end if
! Multiply subdomain vectors by Schur complement
do iinstr = 1,ninstructions
   owner     = instructions(iinstr,1)
   isub      = instructions(iinstr,2)
   is_active = instructions(iinstr,4)

   if (is_active .eq. 1) then
      point  = kbufsend(iinstr)
      length = lbufa(iinstr)

      ! Multiply by Schur complement
      call get_index(isub,indexsub,lindexsub,isub_loc)
      call dd_multiply_by_schur(suba(isub_loc),bufrecv(point),length,bufsend(point),length)
   end if
   if (is_active .eq. 2) then
      point  = kbufsend(iinstr)
      length = lbufa(iinstr)

      ! SUBDOMAIN CORRECTION
      ! prepare array of augmented size
      call get_index(isub,indexsub,lindexsub,isub_loc)
      call dd_get_aug_size(suba(isub_loc), ndofaaugs)
      laux2 = ndofaaugs
      allocate(aux2(laux2))
      call zero(aux2,laux2)
      call dd_get_size(suba(isub_loc), ndofs,nnods,nelems)
      ! truncate the vector for embedding - zeros at the end
      call dd_map_subi_to_sub(suba(isub_loc), bufrecv(point),length, aux2,ndofs)

      nrhs = 1
      call dd_solve_aug(suba(isub_loc), aux2,laux2, nrhs)

      ! get interface part of the vector of preconditioned residual
      call dd_get_number_of_crows(suba(isub_loc),lindrowc)

      lrescs = lindrowc
      allocate(rescs(lrescs))
      call zero(rescs,lrescs)

      ! rc = phis' * x
      transposed = .true.
      call dd_phisi_apply(suba(isub_loc), transposed, bufrecv(point),length, rescs,lrescs)

      call MPI_SEND(rescs,lrescs,MPI_DOUBLE_PRECISION,owner,isub,comm_comm,ierr)
      deallocate(rescs)

      call dd_map_sub_to_subi(suba(isub_loc), aux2,ndofs, bufsend(point),length)
      deallocate(aux2)

   end if
end do
! Wait for all vectors reach their place
nreq = ireq
if (nreq.gt.0) then
   call MPI_WAITALL(nreq, request, statarray, ierr)
end if

! if I am in preconditioning and I own a pair, apply coarse correction of BDDC
ireq = 0
if (idoper.eq.5) then

   ! assemble local residuals to single residual
   call zero(comm_resc,comm_lresc)
   do i = 1,lindrowc_adapt_i
      indc = indrowc_adapt_i(i)

      comm_resc(indc) = comm_resc(indc) + comm_resc_i(i)
   end do
   do i = 1,lindrowc_adapt_j
      indc = indrowc_adapt_j(i)

      comm_resc(indc) = comm_resc(indc) + comm_resc_j(i)
   end do

   ! debug
   !print *, 'resc before'
   !comm_resc = 1._kr
   !print '(f9.3)',comm_resc

   ! solve the local coarse problem by BLAS
   ! u = V * Lambda^-1 * V' f
   laux2 = comm_lresc
   allocate(aux2(laux2))
   ! apply V'
   call DGEMV('Transpose',lcoarsem_adapt1,lcoarsem_adapt2,1._kr,coarsem_adapt,lcoarsem_adapt1,comm_resc,1,0._kr,aux2,1)
   ! apply Lambda^-1
   do i = 1,laux2
      aux2(i) = kceigval(i) * aux2(i)
   end do
   ! apply V
   call DGEMV('Non-transpose',lcoarsem_adapt1,lcoarsem_adapt2,1._kr,coarsem_adapt,lcoarsem_adapt1,aux2,1,0._kr,comm_resc,1)
   deallocate(aux2)

   ! debug
   !print *, 'resc after'
   !print '(e15.6)',comm_resc

   ! pick local solutions
   do i = 1,lindrowc_adapt_i
      indc = indrowc_adapt_i(i)

      comm_resc_i(i) = comm_resc(indc)
   end do
   do i = 1,lindrowc_adapt_j
      indc = indrowc_adapt_j(i)

      comm_resc_j(i) = comm_resc(indc)
   end do

   ! distribute coarse solution
   ireq = ireq + 1
   call MPI_ISEND(comm_resc_i,lindrowc_adapt_i,MPI_DOUBLE_PRECISION,comm_myplace1,comm_myisub,comm_comm,request(ireq),ierr)
   ireq = ireq + 1
   call MPI_ISEND(comm_resc_j,lindrowc_adapt_j,MPI_DOUBLE_PRECISION,comm_myplace2,comm_myjsub,comm_comm,request(ireq),ierr)

end if

! Apply corrections at subdomain level
do iinstr = 1,ninstructions
   owner     = instructions(iinstr,1)
   isub      = instructions(iinstr,2)
   is_active = instructions(iinstr,4)

   if (is_active .eq. 2) then
      point  = kbufsend(iinstr)
      length = lbufa(iinstr)

      call get_index(isub,indexsub,lindexsub,isub_loc)

      ! receive coarse correction
      ! get interface part of the vector of preconditioned residual
      call dd_get_number_of_crows(suba(isub_loc),lindrowc)

      lrescs = lindrowc
      allocate(rescs(lrescs))
      call MPI_RECV(rescs,lrescs,MPI_DOUBLE_PRECISION,owner,isub,comm_comm,stat,ierr)

! COARSE CORRECTION
      call dd_get_interface_size(suba(isub_loc),ndofi,nnodi)
      if (ndofi.ne.length) then
         call error(routine_name,'Interface size mismatch')
      end if

      lsolis = ndofi
      allocate(solis(lsolis))
      call zero(solis,lsolis)

      ! z_i = z_i + phis_i * uc_i
      transposed = .false.
      call dd_phisi_apply(suba(isub_loc), transposed, rescs,lrescs, solis,lsolis)
      deallocate(rescs)

      ! add coarse correction to already prepared subdomain correction
      do i = 1,length
         bufsend(point + i - 1) = bufsend(point + i - 1) + solis(i)
      end do

      deallocate(solis)

      if (suppress_preconditioning) then
         ! plain copy
         do i = 1,length
            bufsend(point + i - 1) = bufrecv(point + i - 1)
         end do
      end if

   end if
end do

! Wait for all vectors reach their place
nreq = ireq
if (nreq.gt.0) then
   call MPI_WAITALL(nreq, request, statarray, ierr)
end if

! distribute multiplied vectors
ireq = 0
! Continue only if I was called from LOBPCG
if (idoper.eq.1 .or. idoper.eq.2 .or. idoper.eq.5) then
   ireq = ireq + 1
   point1 = 1
   call MPI_IRECV(comm_xaux(point1),ndofi_i,MPI_DOUBLE_PRECISION,comm_myplace1,comm_myisub,comm_comm,request(ireq),ierr)

   ireq = ireq + 1
   point2 = ndofi_i + 1
   call MPI_IRECV(comm_xaux(point2),ndofi_j,MPI_DOUBLE_PRECISION,comm_myplace2,comm_myjsub,comm_comm,request(ireq),ierr)
end if
do iinstr = 1,ninstructions
   owner     = instructions(iinstr,1)
   isub      = instructions(iinstr,2)
   is_active = instructions(iinstr,4)

   if (is_active .gt. 0) then
      call MPI_SEND(bufsend(kbufsend(iinstr)),lbufa(iinstr),MPI_DOUBLE_PRECISION,owner,isub,comm_comm,ierr)
   end if
end do

! Wait for all vectors reach their place
nreq = ireq
if (nreq.gt.0) then
   call MPI_WAITALL(nreq, request, statarray, ierr)
end if


! Continue only of I was called from LOBPCG
if (idoper.eq.1 .or. idoper.eq.2 .or. idoper.eq.5) then

   ! reverse sign of vector if this is the A (required for LOBPCG)
   if (idoper.eq.1) then
      do i = 1,problemsize
         comm_xaux(i) = -comm_xaux(i)
      end do
   end if

   if (idoper.eq.1) then
      ! apply (I-RE)^T = (I - E^T R^T) = (I - D_P^T * R * R^T)
      ! copy the array
      do i = 1,problemsize
         comm_xaux2(i) = comm_xaux(i)
      end do
      !  - apply the R^T operator
      call adaptivity_apply_RT(pairslavery,lpairslavery,comm_xaux2,comm_lxaux2)
      !  - apply the R operator
      call adaptivity_apply_R(pairslavery,lpairslavery,comm_xaux2,comm_lxaux2)
      !  - apply weights
      call adaptivity_apply_weights(weight,lweight,comm_xaux2,comm_lxaux2)
      ! Ix - E'R'x
      do i = 1,problemsize
         comm_xaux(i) = comm_xaux(i) - comm_xaux2(i)
      end do
   end if

   ! apply projection Pi - to fulfill common constraints
   ! xaux = P xaux (in all regimes)
   call adaptivity_apply_null_projection(comm_xaux,comm_lxaux)

   ! xaux = P_bar xaux
   if (idoper.eq.2) then
      if (is_nullB_ready) then
         call adaptivity_apply_complementary_projection(comm_xaux,comm_lxaux)
      end if
   end if

   ! copy result to y
   do i = 1,lx
      y(i) = comm_xaux(i)
   end do

   ! debug
   !if (idoper.eq.5) then
   !   print *, 'ly',ly
   !   print *, 'y after'
   !   print '(e16.3)',y
   !end if

end if

end subroutine

!*************************************************************
subroutine adaptivity_apply_complementary_projection(vec,lvec)
!*************************************************************
! Subroutine for application of projection onto complement of null(B)
! P = I - nullB (nullB' nullB) nullB' = I - Q_1Q_1'
! where nullB = QR = [ Q_1 | Q_2 ] R, is the FULL QR decomposition of nullB, Q_1 has
! n columns, which corresponds to number of columns in null(B) 
      use module_utils, only: zero, error

      implicit none
      integer, intent(in) ::    lvec
      real(kr), intent(inout) :: vec(lvec)

! local variables
      character(*),parameter:: routine_name = 'ADAPTIVITY_APPLY_COMPLEMENTARY_PROJECTION'
      integer :: lapack_info

      ! check dimension
      if (lvec.ne.lnullB1 .or. lvec.ne.problemsize) then
         call error(routine_name,'input dimension mismatch')
      end if

      ! apply prepared projection using LAPACK as P = (I-Q_1Q_1') as P = Q_2Q_2',
      ! where Q
      ! xaux = Q_2' * xaux
      call DORMQR( 'Left', 'Transpose',     lvec, 1, lnullB2, nullB, ldnullB, &
                   tau3, vec, lvec, &
                   work3,lwork3, lapack_info)
      if (lapack_info.ne.0) then
         call error(routine_name,'in LAPACK during first application of Q',lapack_info)
      end if
      ! put zeros in first N positions of vec
      call zero(vec,lnullB2)
      ! xaux = Q_2 * xaux
      call DORMQR( 'Left', 'Non-Transpose', lvec, 1, lnullB2, nullB, ldnullB, &
                   tau3, vec, lvec, &
                   work3,lwork3, lapack_info)
      if (lapack_info.ne.0) then
         call error(routine_name,'in LAPACK during second application of Q',lapack_info)
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
      integer :: ldvec, lapack_info

      ! apply prepared projection using LAPACK as P = (I-Q_1Q_1') as P = Q_2Q_2',
      ! where Q
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

!***************************************************
subroutine adaptivity_apply_weights(dp,ldp,vec,lvec)
!***************************************************
! Subroutine for application of weight matrix D_P on vector VEC 
! D_P - stored as diagonal
! vec_out = D_P * vec_in
      use module_utils, only : error_exit
      implicit none
      integer, intent(in)  :: ldp
      real(kr), intent(in) ::  dp(ldp)
      integer, intent(in)  ::    lvec
      real(kr), intent(inout) ::  vec(lvec)

! local variables
      integer :: i

      ! check the length of vector for data
      if (ldp .ne. lvec) then
         write(*,*) 'ADAPTIVITY_APPLY_WEIGHTS: Data size does not match.'
         call error_exit
      end if

      do i = 1,lvec
         vec(i) = dp(i) * vec(i) 
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
      if (any(pair_subdomains(idpair,1:lpair_subdomains2).eq.0)) then
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


