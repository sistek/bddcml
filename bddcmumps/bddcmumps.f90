!***********************************************************************
program bddcmumps
!***********************************************************************
! Program for solving systems of linear equations with Balancing Domain
! Decomposition based on Constraints (BDDC) preconditioning.
!
! Programmed by Jakub Sistek                       Denver,        2/2009
!***********************************************************************

! module for using sparse matrices 
use module_sm
! module for using BDDC preconditioner
use module_bddc1
! module for using MUMPS
use module_mumps
! module with auxiliary utilities
use module_utils
      
implicit none
include "mpif.h"

! Setting of real type kr
integer,parameter:: kr = kind(1.D0)

! Type of matrix - determine its storage, important for MUMPS package
! 0 - unsymmetric                 -> full element matrices
! 1 - symmetric positive definite -> only upper triangle of element matrix
! 2 - symmetric general           -> only upper triangle of element matrix
!  for elasticity, use 1 - matrix is SPD - PCG
!  for Stokes, use 2 - matrix is symmetric indefinite - MINRES
integer,parameter:: matrixtype = 1

! Approach to enforce Gw = 0, ie. averages
! 0 - do not apply averages
! 1 - apply averages via Lagrange multipliers G (S_c)^-1 G^T * mu = G * E^T res
! 2 - apply averages via projection,P = I-G^T*(GG^T)^-1*G
! 3 - change of variables P*T^T*A*T*P & projection
integer,parameter:: averages_approach = 3

! What approach should be used for averaging operator D_P
! 1 - aritmetic averaging (cardinality)
! 2 - weighted averaging - according to diagonal stiffness of original matrix
!  for elasticity, use 2 - robust with respect to jumps in coefficients
!  for Stokes, use 1 - matrix entries may be 0 on diagonal !
integer,parameter:: weight_approach = 1

! Level of information from MUMPS package
! 0 - suppressed printing of all information
! 1 - only errors printed
! 2 - errors and global information printed
! 3 - errors, diagnostics, and global information printed
integer,parameter:: mumpsinfo = 2

! Iterate on staticaly condensed problem on interface?
logical,parameter:: iterate_on_reduced = .true.

! Iterate on transformed problem ? (This option does not work in combination
! with iterations on reduced problem above.)
logical,parameter:: iterate_on_transformed = .false.

! Level of messages of timing
! 0 - no output about times
! 1 - basic times printed
! 2 - rigorous timing output for profiling
integer,parameter:: timeinfo = 1

! Print solution to the screen at the end ?
logical,parameter:: print_solution = .false.

! Use preconditioner?
logical,parameter:: use_preconditioner = .true.

! Use this structure of MUMPS for routines from mumps
type(DMUMPS_STRUC) :: schur_mumps

! Disk units
integer,parameter:: idpar = 1, idfvs = 4, &
                    idelm = 10, idrhs = 11, idint = 12, &
                    idrhss = 15, idfvss = 16, idsols = 17, idsol = 18, idglb = 19, &
                    idgmist = 20, idtr = 21

! Lenght of names 
integer,parameter:: lname1x = 8, lnamex = 15, lfnamex = 20
! Name of problem
character(lname1x)::  name1
character(lnamex) ::  name
character(lfnamex)::  fname

! Global parameters of problem
integer:: lname1, ndim, nsub, nelem, ndof, nnod, nnodc, linet, maxit, ndecrmax, &
          nglb
real(kr):: tol, pival

! Variables for communication
integer:: comm, myid, nproc, ierr
integer:: start, finish, nsub_loc, nsub_locx

! Mesh description in W_hat
integer::            lifix
integer,allocatable:: ifix(:)

! Mesh description in W_tilde - basic GLOBAL dimensions
integer:: nnodt, ndoft
! Mesh description in W_tilde - basic LOCAL dimensions
integer:: nelems , linets

! Mesh description in W_tilde - local arrays to subdomain
integer::            linetst,   lnnetst
integer,allocatable:: inetst(:), nnetst(:)

! Mesh description in W_tilde - global arrays
integer ::            lnndft,   lifixt,   lkdoft,   lslavery,   lihntn
integer,allocatable :: nndft(:), ifixt(:), kdoft(:), slavery(:), ihntn(:)
 
! Mesh description in W and W_tilde - interface
integer ::            ligingn,   liintt
integer,allocatable :: igingn(:), iintt(:)

! Real arrays in W_hat
integer ::            lrhs,   lfixv,   lsol
real(kr),allocatable:: rhs(:), fixv(:), sol(:)

! Real arrays in W_tilde
integer::             lrhst,   lfixvt,   lsolt,   lbct = 0, lsolintt
real(kr),allocatable:: rhst(:), fixvt(:), solt(:), bct(:),   solintt(:)

! Matrix in IJA sparse format - triplet
integer::  la, nnz
integer,allocatable  :: i_sparse(:), j_sparse(:)
real(kr),allocatable :: a_sparse(:)

! Description of globs
integer ::            linglb,   lnnglb
integer,allocatable :: inglb(:), nnglb(:)

! Variables for measuring time
real(kr):: time

! Coefficient of sparsity of the matrix
real(kr):: sparsity

! Local variables
integer :: i, j, inodt, iaux, isub, nnz_proj_est, nnz_tr_proj_est, la_pure, la_proj, la_transform,&
           ndofs, inod, isol, ndofn, nnodi, ini, indnt, ndofnt

! auxiliary array
integer ::             laux
real(kr),allocatable :: aux(:)

logical :: nonzero_bc_loc = .false. , nonzero_bc
logical :: is_this_the_first_run

logical :: match_mask(2)

logical :: parallel_analysis

character(100) :: filename, problemname

! MPI initialization
!***************************************************************PARALLEL
! Communicator
      comm = MPI_COMM_WORLD
      call MPI_INIT(ierr)
      call MPI_COMM_RANK(comm,myid,ierr)
      call MPI_COMM_SIZE(comm,nproc,ierr)
!***************************************************************PARALLEL

! Beginning of mesuring time
      call bddc_time_start(comm)
      
! Initial screen
      if (myid.eq.0) then
         write(*,'(a)') 'BDDCMUMPS - BDDC solver using MUMPS'
         write(*,'(a)') '==================================='

! Name of the problem
   10    write(*,'(a,$)') 'Name of the problem: '
         read(*,*) name1
         if(name1.eq.' ') goto 10
         lname1   = index(name1,' ') - 1
         if(lname1.eq.-1) then
           lname1 = lname1x
         end if
         call flush(6)
      end if
! Broadcast of name of the problem      
!***************************************************************PARALLEL
      call MPI_BCAST(name1, 8, MPI_CHARACTER, 0, comm, ierr)
      call MPI_BCAST(lname1,  1, MPI_INTEGER, 0, comm, ierr)
!***************************************************************PARALLEL

! Open disk files on root processor
      if (myid.eq.0) then
! PAR - basic properties of the problem    
         name = name1(1:lname1)//'.PAR'
         open (unit=idpar,file=name,status='old',form='formatted')

! FVS - fixed variables
!  * IFIX(LIFIX) * FIXV(LFIXV) * 
         name = name1(1:lname1)//'.FVS'
         open (unit=idfvs,file=name,status='old',form='formatted')

! INT - interface variables
!  * IGINGN(LIGINGN) *
         if (iterate_on_reduced) then
            name = name1(1:lname1)//'.INT'
            open (unit=idint,file=name,status='old',form='formatted')
         end if

! RHS - global right hand side 
!  * RHS(LRHS) * 
         name = name1(1:lname1)//'.RHS'
         open (unit=idrhs,file=name,status='old',form='unformatted')
        
! SOL - global solution
         name = name1(1:lname1)//'.SOL'
         open (unit=idsol,file=name,status='replace',form='unformatted')
      end if

! Reading basic properties 
      if (myid.eq.0) then
         read(idpar,*) ndim, nsub, nelem, ndof, nnod, nnodc, linet, &
                       tol, maxit, ndecrmax, iaux, pival

         write(*,*)'Characteristics of the problem ',name1(1:lname1),':'
         write(*,*)'  number of processors            nproc =',nproc
         write(*,*)'  number of dimensions             ndim =',ndim
         write(*,*)'  number of subdomains             nsub =',nsub
         write(*,*)'  number of elements global       nelem =',nelem
         write(*,*)'  number of DOF                    ndof =',ndof
         write(*,*)'  number of nodes global           nnod =',nnod
         write(*,*)'  number of constrained nodes     nnodc =',nnodc
         write(*,*)'  lenght of field INET            linet =',linet
         write(*,*)'Characteristics of iterational process:'
         write(*,*)'  tolerance of error                tol =',tol
         write(*,*)'  maximum number of iterations    maxit =',maxit
         write(*,*)'  number of incresing residual ndecrmax =',ndecrmax
         write(*,*)'Characteristics of frontal solution:'
         write(*,*)'  minimal pivot                   pival =',pival
         call flush(6)
      end if
! Broadcast basic properties of the problem
!***************************************************************PARALLEL
      call MPI_BCAST(ndim,     1, MPI_INTEGER,         0, comm, ierr)
      call MPI_BCAST(nsub,     1, MPI_INTEGER,         0, comm, ierr)
      call MPI_BCAST(nelem,    1, MPI_INTEGER,         0, comm, ierr)
      call MPI_BCAST(ndof,     1, MPI_INTEGER,         0, comm, ierr)
      call MPI_BCAST(nnod,     1, MPI_INTEGER,         0, comm, ierr)
      call MPI_BCAST(nnodc,    1, MPI_INTEGER,         0, comm, ierr)
      call MPI_BCAST(linet,    1, MPI_INTEGER,         0, comm, ierr)
      call MPI_BCAST(tol,      1, MPI_DOUBLE_PRECISION,0, comm, ierr)
      call MPI_BCAST(maxit,    1, MPI_INTEGER,         0, comm, ierr)
      call MPI_BCAST(ndecrmax, 1, MPI_INTEGER,         0, comm, ierr)
      call MPI_BCAST(pival,    1, MPI_DOUBLE_PRECISION,0, comm, ierr)
!***************************************************************PARALLEL

! Check the demanded number of processors
      if (nproc.ne.nsub) then
         if (myid.eq.0) then
            write(*,*) 'ERROR: This program has to run on nproc = ', nsub, &
                       'processors exclusively.'
         end if
         stop
      end if

! Distribution of subdomains to processors
      nsub_locx = (nsub + nproc - 1)/nproc
      nsub_loc = min(nsub - myid*nsub_locx, nsub_locx)
      start = nsub_locx * myid + 1
      finish = start + nsub_loc - 1
      write(*,*) 'myid =',myid,'start = ',start,'finish = ',finish
! change for more subdomains per processor 
      isub = start

! GMISTS - basic mesh data for subdomain in W_tilde space
!  * INETTS(LINETS) * NNETTS(NELEMS) * NNDFT(NNODT) * SLAVERY(NNODT) * IHNTN(NNOD)
      call bddc_getfname(name1,lname1,isub,'GMISTS',fname)
      open (unit=idgmist,file=fname,status='old',form='formatted')
      read(idgmist,*)  nnodt, ndoft
      read(idgmist,*)  nelems, linets, ndofs
      linetst = linets
      lnnetst = nelems
      lnndft  = nnodt
      lkdoft  = nnodt
      lslavery = nnodt
      lihntn  = nnod
      allocate(inetst(linetst),nnetst(lnnetst),nndft(lnndft),slavery(lslavery), &
               ihntn(lihntn),kdoft(lkdoft))
      read(idgmist,*) inetst
      read(idgmist,*) nnetst
      read(idgmist,*) nndft
      read(idgmist,*) slavery
      read(idgmist,*) ihntn

! Creation of field KDOFT(NNODT) with addresses before first global dof of node
      kdoft(1) = 0
      do inodt = 2,nnodt
         kdoft(inodt) = kdoft(inodt-1) + nndft(inodt-1)
      end do

! Read fixed variables
      lifix = ndof
      lfixv = ndof
      allocate(ifix(lifix),fixv(lfixv))
      if (myid.eq.0) then
         read(idfvs,*) ifix
         read(idfvs,*) fixv
      end if
!*****************************************************************MPI
      call MPI_BCAST(ifix,lifix, MPI_INTEGER,          0, comm, ierr)
      call MPI_BCAST(fixv,lfixv, MPI_DOUBLE_PRECISION, 0, comm, ierr)
!*****************************************************************MPI
      ! convert IFIX and FIXV to W_tilde
      lifixt = ndoft
      lfixvt = ndoft
      allocate(ifixt(lifixt),fixvt(lfixvt))
      call bddc_convert_ht(nnod,nnodt,nndft,lnndft,ihntn,lihntn,slavery,lslavery,kdoft,lkdoft,fixv,lfixv,fixvt,lfixvt)
      call bddc_convert_ht_int(nnod,nnodt,nndft,lnndft,ihntn,lihntn,slavery,lslavery,kdoft,lkdoft,ifix,lifix,ifixt,lifixt)
      ! keep only arrays in W_tilde
      deallocate(ifix,fixv)

! Read interface nodes
      if (iterate_on_reduced) then
         if (myid.eq.0) then
            read(idint,*) nnodi
            ligingn = nnodi
            allocate(igingn(ligingn))
            read(idint,*) igingn
         end if
         liintt = ndoft
         allocate(iintt(liintt))
         iintt = 0
         if (myid.eq.0) then
            do ini = 1,nnodi
               indnt  = ihntn(igingn(ini))
               ndofnt = nndft(indnt)
   
               iintt(kdoft(indnt)+1:kdoft(indnt)+ndofnt) = 1
            end do
            deallocate(igingn)
         end if
!*****************************************************************MPI
         call MPI_BCAST(iintt,liintt, MPI_INTEGER,        0, comm, ierr)
!*****************************************************************MPI
         call bddc_R_int(nnodt,nndft,lnndft,slavery,lslavery,kdoft,lkdoft,iintt,liintt)
      end if


! Approach to averages
      select case (averages_approach)
      case(0)
         if (myid.eq.0) then
            write(*,*) 'Approach to averages: NO AVERAGES APPLIED'
            call flush(6)
         end if
      case(1)
         if (myid.eq.0) then
            write(*,*) 'Approach to averages: LAGRANGE MULTIPLIERS'
            call flush(6)
         end if
      case(2)
         if (myid.eq.0) then
            write(*,*) 'Approach to averages: PROJECTION ON NULLSPACE OF G'
            call flush(6)
         end if
      case(3)
         if (myid.eq.0) then
            write(*,*) 'Approach to averages: CHANGE OF VARIABLES + PROJECTION'
            call flush(6)
         end if
      case default
         if (myid.eq.0) then
            write(*,*) 'bddcmumps: Illegal value of AVERAGES_APPROACH.', averages_approach
            call flush(6)
         end if
         stop
      end select
! Glob data
      ! GLB - glob information
      !  * nglb * linglb
      !  * inglb(linglb) *
      !  * nnglb(nglb) *
      select case (averages_approach)
      case(0)
         ! using no globs requested
         nglb   = 0
         linglb = 0
      case(1)
         ! For Lagrange multipliers, GLOBAL arrays of glob description are necessary
         ! Read number of globs
         if (myid.eq.0) then
            name = name1(1:lname1)//'.GLB'
            open (unit=idglb,file=name,status='old',form='formatted')
            read(idglb,*) nglb, linglb
         end if
         !***************************************************************PARALLEL
         call MPI_BCAST(nglb, 1, MPI_INTEGER, 0, comm, ierr)
         call MPI_BCAST(linglb, 1, MPI_INTEGER, 0, comm, ierr)
         !***************************************************************PARALLEL
         linglb = linglb
         lnnglb = nglb
         allocate(inglb(linglb),nnglb(lnnglb))
         if (myid.eq.0) then
            read(idglb,*) inglb
            read(idglb,*) nnglb
            close(idglb)
         end if
         !***************************************************************PARALLEL
         call MPI_BCAST(inglb,linglb, MPI_INTEGER, 0, comm, ierr)
         call MPI_BCAST(nnglb,lnnglb, MPI_INTEGER, 0, comm, ierr)
         !***************************************************************PARALLEL
      case(2,3)
         ! For projection and transformation, LOCAL glob description is used
         call bddc_getfname(name1,lname1,isub,'GLB',fname)
         open (unit=idglb,file=fname,status='old',form='formatted')
         read(idglb,*) nglb, linglb
         linglb = linglb
         lnnglb = nglb
         allocate(inglb(linglb),nnglb(lnnglb))
         read(idglb,*) inglb
         read(idglb,*) nnglb
         close(idglb)
      end select

! Load sparse matrix
      ! find the length for matrix entries
      call sm_pmd_get_length(matrixtype,nelems,inetst,linetst,nnetst,lnnetst,nndft,lnndft,la_pure)
      ! estimate sparsity coefficient
      sparsity = real(la_pure,kr)/(ndofs*ndofs)
      write(*,*) 'myid =',myid,': Sparsity =',sparsity
      ! add space for entries generated by change of variables and/or projection
      la_proj      = 0
      la_transform = 0
      if      (averages_approach.eq.2) then
         call bddc_P_nnz_est(matrixtype,sparsity,ndofs,ndim,nglb,inglb,linglb,nnglb,lnnglb,&
                             slavery,lslavery, nnz_proj_est)
         la_proj = nnz_proj_est
         write(*,*) 'myid =',myid,': Space estimated for projection =',la_proj
      else if (averages_approach.eq.3) then
         call bddc_T_nnz_est(myid,matrixtype,sparsity,ndofs,ndim,nglb,inglb,linglb,nnglb,lnnglb,&
                             slavery,lslavery,nnz_tr_proj_est)
         la_transform = nnz_tr_proj_est
         write(*,*) 'myid =',myid,': Space estimated for change of basis =',la_transform
      end if
      la = la_pure + la_proj + la_transform
      ! prepare memory for distributed sparse matrix IJA
      write(*,*) 'myid =',myid,': Trying to allocate for sparse matrix',2*sm_showsize(la,8),' [MB].'
      allocate(i_sparse(la), j_sparse(la), a_sparse(la))

! Run this twice - in the first run, prepare interior solution u_int
! (if iteration on reduced problem is requested)
! A_11 u_int = f_int
! In the second run, multipy it by A_21 and eliminate this from rhs
      is_this_the_first_run = .true.
! jump here after the interior solution was found
      lrhst = ndoft
      allocate(rhst(lrhst))
      if (iterate_on_reduced) then
         lsolintt = ndoft
         allocate(solintt(lsolintt))
      end if
      lsolt = ndoft
      allocate(solt(lsolt))
 20   continue
      
      if (is_this_the_first_run) then
         if (myid.eq.0) then
            write (*,*) 'Building the system for the FIRST time to find the internal solution.'
            call flush(6)
         end if
      else
         if (myid.eq.0) then
            write (*,*) 'Building the system for the SECOND time.'
            call flush(6)
         end if
      end if

! Read right hand side
      lrhs = ndof
      allocate(rhs(lrhs))
      if (myid.eq.0) then
         rewind idrhs
         read(idrhs) rhs
      end if
!***************************************************************PARALLEL
      call MPI_BCAST(rhs,lrhs, MPI_DOUBLE_PRECISION, 0, comm, ierr)
!***************************************************************PARALLEL
      ! convert rhs to W_tilde
      call bddc_convert_ht(nnod,nnodt,nndft,lnndft,ihntn,lihntn,slavery,lslavery,kdoft,lkdoft,rhs,lrhs,rhst,lrhst)
      ! keep only field in tilde
      deallocate(rhs)

! ELMS - element stiffness matrices of subdomain - structure:
      problemname = name1(1:lname1)
      call getfname(trim(problemname),isub,'ELM',filename)
      open(unit=idelm,file=filename,status='old',form='unformatted')
      call sm_pmd_load(idelm,nelems,inetst,linetst,nnetst,lnnetst,nndft,lnndft,kdoft,lkdoft,&
                       i_sparse, j_sparse, a_sparse, la_pure)
      close(idelm)

! debug      
!      call sm_print(6,i_sparse, j_sparse, a_sparse,la_pure, la_pure)
!      ltestm1 = ndofs
!      ltestm2 = ndofs
!      allocate(testm(ltestm1,ltestm2))
!      call sm_to_dm(matrixtype,i_sparse, j_sparse, a_sparse,la_pure, testm,ltestm1,ltestm2)
!      write(88,'(i8)') ndofs
!      do i = 1,ltestm1
!         write(88,'(1000f10.4)') (testm(i,j),j = 1,ltestm2)
!      end do
!      write(88,'(f10.4)') rhst


! Apply boundary conditions and prepare vector BC with effect of non-homogenous boudary conditions
      if (any(ifixt.ne.0.and.fixvt.ne.0.0_kr)) then
         nonzero_bc_loc = .true.
      end if
!*****************************************************************MPI
      call MPI_ALLREDUCE(nonzero_bc_loc,nonzero_bc,1,MPI_LOGICAL,MPI_LOR,comm,ierr)
!*****************************************************************MPI
      ! If there are nonzero Dirichlet BC, we will need array BCT to eliminate them
      if (nonzero_bc) then
         lbct   = ndoft
         allocate(bct(lbct))
      end if
      ! eliminate natural BC
      call sm_apply_bc(ifixt,lifixt,fixvt,lfixvt,i_sparse,j_sparse,a_sparse,la_pure, bct,lbct)
      if (nonzero_bc) then
         call bddc_RRT(comm,nnodt,nndft,lnndft,slavery,lslavery,kdoft,lkdoft,bct,lbct)
      end if
! Prepare RHS
      call sm_prepare_rhs(ifixt,lifixt,bct,lbct,rhst,lrhst)
      if (iterate_on_reduced) then
         ! eliminate in the first run also the interface variables - fixed to zeros
         if (is_this_the_first_run) then
            laux = lfixvt
            allocate(aux(laux))
            aux = 0.0_kr
            call sm_apply_bc(iintt,liintt,aux,laux,i_sparse,j_sparse,a_sparse,la_pure, bct,lbct)
            deallocate(aux)
            if (nonzero_bc) then
               call bddc_RRT(comm,nnodt,nndft,lnndft,slavery,lslavery,kdoft,lkdoft,bct,lbct)
            end if
            call sm_prepare_rhs(iintt,liintt,bct,lbct,rhst,lrhst)
         end if
      end if

! Assembly entries in matrix
      call sm_assembly(i_sparse,j_sparse,a_sparse,la_pure,nnz)
      write(*,*) 'myid =',myid,': Matrix loaded and assembled'
      call flush(6)

      ! remove BCT array if it was used
      if (allocated(bct)) then
         deallocate(bct)
      end if

! debug
!      call sm_to_dm(matrixtype,i_sparse, j_sparse, a_sparse,nnz, testm,ltestm1,ltestm2)
!      write(89,'(i8)') ndofs
!      do i = 1,ltestm1
!         write(89,'(1000f10.4)') (testm(i,j),j = 1,ltestm2)
!      end do
!      write(89,'(f10.4)') rhst
!      deallocate(testm)

      if (iterate_on_reduced) then
         if (is_this_the_first_run) then
            ! Initialize MUMPS
            call mumps_init(schur_mumps,comm,matrixtype)
            ! Level of information from MUMPS
            call mumps_set_info(schur_mumps,mumpsinfo)
            ! Load matrix to MUMPS
            call mumps_load_triplet(schur_mumps,ndoft,nnz,i_sparse,j_sparse,a_sparse,la)
            ! Analyze matrix
            parallel_analysis = .true.
            call mumps_analyze(schur_mumps,parallel_analysis) 
            ! Factorize matrix 
            call mumps_factorize(schur_mumps)
            ! Solve the system for given rhs
            call mumps_resolve(schur_mumps,rhst,lrhst)
!*****************************************************************MPI
            call MPI_BCAST(rhst,lrhst, MPI_DOUBLE_PRECISION, 0, comm, ierr)
!*****************************************************************MPI

            ! Now rhst contains solution in interior nodes u_0
            solintt = rhst
   
            is_this_the_first_run = .false.
            goto 20
   
         else
      ! Construct g - reduced RHS

      ! Multiply the block of the matrix A21 and rhst
            match_mask(1) = .true.
            match_mask(2) = .false.
            call sm_vec_mult_mask(matrixtype, nnz, i_sparse, j_sparse, a_sparse, la, &
                                  solintt,lsolintt, solt,lsolt, &
                                  iintt,liintt,match_mask)
            call bddc_RRT(comm,nnodt,nndft,lnndft,slavery,lslavery,kdoft,lkdoft,solt,lsolt)

            where (iintt.eq.0) rhst = 0.0_kr

      ! Check that A12*u01 has zeros in interiors
            if (any(iintt.eq.0.and.solt.ne.0.0_kr)) then
               write(*,*) 'Solution should have zeros in interiors.'
               stop
            end if

            rhst = rhst - solt

      ! Now rhst should have the structure (0 g)^T
         end if
      end if


! Initial vector of solution - prescribed fixed variables, zero elsewhere
      solt = 0.0_kr
      if (iterate_on_reduced) then
         where(iintt.eq.1) solt = fixvt
      else
         where(ifixt.ne.0) solt = fixvt
      end if
      if (iterate_on_reduced) then
         if (any(solt.ne.0.0_kr)) then
! Project solt onto the space of energy minimal functions over interiors
            where(iintt.eq.0) solt = 0.0_kr
            ! Multiply the block of the matrix A12 and p
            laux = lsolt
            allocate(aux(laux))
            ! A12 * solt
            match_mask(1) = .false.
            match_mask(2) = .true.
            call sm_vec_mult_mask(matrixtype, nnz, i_sparse, j_sparse, a_sparse, la, &
                                  solt,lsolt, aux,laux, &
                                  iintt,liintt, match_mask)
            call bddc_RRT(comm,nnodt,nndft,lnndft,slavery,lslavery,kdoft,lkdoft,aux,laux)
            aux = -aux
            ! Solve the system for given rhs
            call mumps_resolve(schur_mumps,aux,laux)
!*****************************************************************MPI
            call MPI_BCAST(aux,laux, MPI_DOUBLE_PRECISION, 0, comm, ierr)
!*****************************************************************MPI
            where(iintt.eq.0) solt = aux
            deallocate(aux)
      !debug
! Check that the vector is really energy - minimal
!         call sm_vec_mult(matrixtype, nnz, i_sparse, j_sparse, a_sparse, la, &
!                          p,lp, ap,lap)
!      if (any (abs(ap).gt.10e-16 * normrhs .and. iintt.eq.0)) then
!         write(*,*) 'Energy minimality check failed!'
!         write(*,*) 'Sum of error entries',sum(abs(ap),mask=iintt.eq.0)
!         stop
!      end if
!      write(*,*) 'myid = ',myid,'p after preconditioner'
!      write(*,'(f15.9)') p(ihntn)
            end if
      end if

! Open file for storing transformation matrices
      if (averages_approach.eq.3) then
         call bddc_getfname(name1,lname1,isub,'TR',fname)
         open(unit=idtr,file=fname,status='replace',form='unformatted')
      end if

! Call Krylov method for solution of the system
!      if      (matrixtype.eq.1) then
         if (myid.eq.0) then
            write (*,*) 'Calling PCG method for solution.'
            call flush(6)
         end if
         call bddcpcg(myid,comm,iterate_on_reduced,iterate_on_transformed,matrixtype,nnz,a_sparse,i_sparse,j_sparse,la, &
                      mumpsinfo,timeinfo,use_preconditioner,weight_approach,averages_approach,&
                      ndim,nglb,inglb,linglb,nnglb,lnnglb,&
                      nnodt,ndoft,ihntn,lihntn,slavery,lslavery,kdoft,lkdoft, nndft,lnndft, iintt,liintt, &
                      idtr, schur_mumps, rhst,lrhst,tol,maxit,ndecrmax,solt,lsolt)
!      else if (matrixtype.eq.2) then
!         if (myid.eq.0) then
!            write (*,*) 'Calling MINRES method for solution.'
!            call flush(6)
!         end if
!         call bddcminres(myid,comm,iterate_on_reduced,iterate_on_transformed,matrixtype,nnz,a_sparse,i_sparse,j_sparse,la, &
!                         mumpsinfo,timeinfo,use_preconditioner,weight_approach,averages_approach,&
!                         ndim,nglb,inglb,linglb,nnglb,lnnglb,&
!                         nnodt,ndoft,ihntn,lihntn,slavery,lslavery,kdoft,lkdoft, nndft,lnndft, iintt,liintt, &
!                         idtr, schur_mumps, rhst,lrhst,tol,maxit,ndecrmax,solt,lsolt)
!      else
!         if (myid.eq.0) then
!            write (*,*) 'Matrixtype not supported by iterations.', matrixtype
!         end if
!         call error_exit
!      end if


! Close file for transformation matrices
      if (averages_approach.eq.3) then
         close (idtr)
      end if

! Now add interface (energy minimal on interiors) solution SOLT to interior solution SOLINTT
      if (iterate_on_reduced) then
!         where (iintt.eq.0) solt = solt + solintt
         solt = solt + solintt
      end if

! Convert solution to W_hat
      lsol = ndof
      allocate(sol(lsol))
      call bddc_convert_th(nnod,nndft,lnndft,ihntn,lihntn,kdoft,lkdoft,solt,lsolt,sol,lsol)
      deallocate(solt)
      if (allocated(solintt)) then
         deallocate(solintt)
      end if

! Write solution
      if (myid.eq.0) then
         if (print_solution) then
            write(*,*) 'Solution is ...'
            isol = 0
            do inod = 1,nnod
               ndofn = nndft(ihntn(inod))
               write(*,'(i5, 5f15.10)') inod, (sol(isol + j), j = 1,ndofn)
               isol = isol + ndofn
            end do
         end if

         write(*,*) 'Writing solution on disk ...'
         call flush(6)
         write(idsol) (sol(i), i = 1,lsol)
! Nonsense writing in place where postprocessor expect reactions
         write(idsol) (sol(i), i = 1,lsol), 0.0_kr, 0.0_kr, 0.0_kr
         write(*,*) '...done'
         write(*,*) 'Solution has been written into file ',name1(1:lname1),'.SOL'
         write(*,*) 'Warning: At the moment solver does not ',&
                    'resolve reaction forces. Record of these does not',&
                    ' make sense and is present only to make ',&
                    'postprocessor str3 happy with input data.'
         call flush(6)
         
         close (idsol)
      end if

! Clear memory
      ! Finalize the instance of MUMPS for Schur complement
      if (iterate_on_reduced) then
         call mumps_finalize(schur_mumps)
      end if

      if (averages_approach.ne.0) then
         deallocate(inglb,nnglb)
      end if
      if (allocated(iintt)) then
         deallocate(iintt)
      end if
      deallocate(inetst,nnetst,nndft,slavery,ihntn,kdoft)
      deallocate(rhst)
      deallocate(ifixt,fixvt)
      deallocate(a_sparse,i_sparse,j_sparse)
      deallocate(sol)

! End of mesuring time
      call bddc_time_end(comm,time)
      if (myid.eq.0.and.timeinfo.ge.1) then
         write(*,*) '=========================='
         write(*,*) 'Time of whole run = ',time
         write(*,*) '=========================='
         call flush(6)
      end if

! Finalize MPI
!***************************************************************PARALLEL
      call MPI_FINALIZE(ierr)
!***************************************************************PARALLEL

      write(*,*) 'myid =',myid,': O.K.'
      
end program

!******************************************************************************************************
subroutine bddcpcg(myid,comm,iterate_on_reduced,iterate_on_transformed,matrixtype,nnz,a_sparse,i_sparse,j_sparse,la, &
                   mumpsinfo,timeinfo,use_preconditioner,weight_approach,averages_approach,ndim,nglb,inglb,linglb,nnglb,lnnglb,&
                   nnodt,ndoft,ihntn,lihntn,slavery,lslavery,kdoft,lkdoft,nndft,lnndft, iintt, liintt,&
                   idtr, schur_mumps, res,lres,tol,maxit,ndecrmax,sol,lsol)
!******************************************************************************************************
! Preconditioned conjugate gradient solver based on BDDC
!******************************************************************************************************
! Use module for BDDC
      use module_bddc1

! Use module for sparse matrices
      use module_sm

! Use module MUMPS
      use module_mumps

! Use module with utilities
      use module_utils

      implicit none
      include "mpif.h"

! Setting of real type kr
      integer,parameter:: kr = kind(1.D0)

! Iterate on staticaly condensed problem ?
      logical,intent(in) :: iterate_on_reduced

! Iterate on transformed problem
      logical,intent(in) :: iterate_on_transformed

! Matrix in sparse IJA format
      integer,intent(in) :: matrixtype
      integer,intent(in) :: nnz, la
      integer,intent(inout) :: i_sparse(la),j_sparse(la)
      real(kr),intent(inout):: a_sparse(la)

! Verbose level of MUMPS
      integer,intent(in) :: mumpsinfo

! Verbose level of times
      integer,intent(in) :: timeinfo

! Use preconditioner?
      logical,intent(in) :: use_preconditioner

! Approach to averaging operator D_P
      integer,intent(in) :: weight_approach

! Description of globs
      integer,intent(in) :: averages_approach
      integer,intent(in) :: ndim
      integer,intent(in) :: nglb
      integer,intent(in) :: linglb,        lnnglb
      integer,intent(in) ::  inglb(linglb), nnglb(lnnglb)

! Space W_tilde
      integer,intent(in) :: nnodt, ndoft
      integer,intent(in) :: lihntn,        lslavery,          lkdoft,        lnndft,        liintt
      integer,intent(in) ::  ihntn(lihntn), slavery(lslavery), kdoft(lkdoft), nndft(lnndft), iintt(liintt)

! Variables for communication
      integer,intent(in):: comm, myid

! Disk unit for transformation
      integer,intent(in):: idtr

! Use this structure of MUMPS for routines from mumps
      type(DMUMPS_STRUC) :: schur_mumps

! Right hand side
      integer,intent(in) ::    lres
      real(kr),intent(inout)::  res(lres)
! Solution
      integer,intent(in)::    lsol
      real(kr),intent(inout):: sol(lsol)

! Tolerance
      real(kr),intent(in) :: tol
! Maximum number of iterations
      integer,intent(in) :: maxit
! Maximum number of iterations without residual decrease
      integer,intent(in) :: ndecrmax
      
! Local variables of PCG algorithm
      real(kr):: alpha, beta, normres, normrhs, relres, lastres, rmr, rmrold, pap, pap_loc
      integer::              lp,   lap,   lh
      real(kr), allocatable:: p(:), ap(:), h(:)

      integer:: ierr, iter, ndecr, nnz_mult
      real(kr):: time
      logical:: nonzero_initial_solution_loc = .false. , nonzero_initial_solution
      integer:: nnz_transform

! Variables for condition number estimation
      integer ::             lw,   lwchol,   lx,   ly
      real(kr),allocatable :: w(:), wchol(:), x(:), y(:)
      integer :: nw, nwx  
      real(kr) :: cond

      logical :: match_mask(2)

      !integer :: i

! Beginning of mesuring time
      call bddc_time_start(comm)

! Prepare fields for PCG - residual will be associated with RHS
      lp   = ndoft
      lap  = ndoft
      lh   = ndoft
      allocate(p(lp),ap(lap),h(lh))

!***************************************************CONDITION NUMBER ESTIMATION
! Prepare matrix for condition number estimation
      lw     = (maxit+1)*(maxit+1)
      lwchol = lw
      lx     = maxit + 1
      ly     = maxit + 1
      allocate(w(lw),wchol(lw),x(lx),y(ly))
      w     = 0.0_kr
      wchol = 0.0_kr
      x     = 0.0_kr
      y     = 0.0_kr
!***************************************************CONDITION NUMBER ESTIMATION

!****************************************************************************BDDC
! Beginning of mesuring time
      call bddc_time_start(comm)
! Initialize an instance of the MUMPS package
      call bddc_init(myid,comm,matrixtype,mumpsinfo,timeinfo,use_preconditioner,iterate_on_transformed,&
                     ndoft,nnz,i_sparse,j_sparse,a_sparse,la, &
                     weight_approach,averages_approach,ndim,nglb,inglb,linglb,nnglb,lnnglb,&
                     nnodt,nndft,lnndft,ihntn,lihntn,slavery,lslavery,kdoft,lkdoft, &
                     sol,lsol, res,lres, nnz_transform,idtr)
      call bddc_time_end(comm,time)
      if (myid.eq.0.and.timeinfo.ge.1) then
         write(*,*) 'Initialized BDDC'
         write(*,*) '===================================='
         write(*,*) 'Time of initializing of BDDC = ',time
         write(*,*) '===================================='
         call flush(6)
      end if
!****************************************************************************BDDC

! Iterate on transformed problem?
      if (iterate_on_transformed) then
         nnz_mult = nnz+nnz_transform
      else
         nnz_mult = nnz
      end if

! Prepare initial residual
! res_0 = rhs - A*sol_0
      if (any(sol.ne.0_kr)) then
         if (myid.eq.0) then
            write(*,*) 'Nonzero initial solution -> create initial vector of residual.'
            nonzero_initial_solution_loc = .true.
         end if
      end if
!*****************************************************************MPI
      call MPI_ALLREDUCE(nonzero_initial_solution_loc,nonzero_initial_solution,1, &
                         MPI_LOGICAL,MPI_LOR,comm,ierr)
!*****************************************************************MPI
      if (nonzero_initial_solution) then
         if (iterate_on_reduced) then
            ! A21 * sol
            match_mask(1) = .true.
            match_mask(2) = .false.
            call sm_vec_mult_mask(matrixtype,nnz_mult, i_sparse, j_sparse, a_sparse, la, sol,lsol, ap,lap, &
                                  iintt,liintt,match_mask)
            ! A22 * sol
            match_mask(1) = .true.
            match_mask(2) = .true.
            call sm_vec_mult_mask(matrixtype,nnz_mult, i_sparse, j_sparse, a_sparse, la, sol,lsol, h,lh, &
                                  iintt,liintt,match_mask)
            ap = ap + h
         else
            call sm_vec_mult(matrixtype,nnz_mult, i_sparse, j_sparse, a_sparse, la, sol,lsol, ap,lap)
         end if
         call bddc_RRT(comm,nnodt,nndft,lnndft,slavery,lslavery,kdoft,lkdoft,ap,lap)
         res = res - ap
      end if
      if (iterate_on_reduced) then
         if (any(res.ne.0.0_kr .and. iintt.eq.0)) then
            write(*,*) 'WARNING: Initial check of energy minimal function failed.'
         end if
         where(iintt.eq.0) res = 0.0_kr
      end if
     
! determine norm of initial rhs
! ||rhs||
      normrhs = bddc_normvec(comm,res,lres)

! Check of zero right hand side => all zero solution
      if (.not.any(res.ne.0.0_kr)) then
         if (myid.eq.0) then
            write(*,*) 'all zero RHS => solution without change'
         end if
         goto 77
      end if
     
! Initial action of the preconditioner
! p = M_BDDC*res 
      if (myid.eq.0) then
         write(*,*) ' Initial action of preconditioner'
         call flush(6)
      end if
      if (iterate_on_reduced) then
      ! check that residual is energy minimal
         if (any(abs(res).gt.1e-15_kr .and. iintt .eq. 0)) then
            write(*,*) 'Residual is not energy minimal!!!!!!!!!!!!!!!!'
            stop
         end if
      end if
!      write(*,*) 'myid = ',myid,'res before preconditioner'
!      write(*,'(f15.9)') res(ihntn)
      if (use_preconditioner) then
         call bddc_M(myid,comm,iterate_on_transformed,nnodt,nndft,lnndft,slavery,lslavery,kdoft,lkdoft,res,lres,p,lp,rmr)
      else
         call bddc_M_fake(comm,res,lres,p,lp,rmr)
      end if
      if (iterate_on_reduced) then
! Project p onto the space of energy minimal functions over interiors
         ! Multiply the block of the matrix A12 and p
         match_mask(1) = .false.
         match_mask(2) = .true.
         call sm_vec_mult_mask(matrixtype, nnz, i_sparse, j_sparse, a_sparse, la, &
                               p,lp, ap,lap, &
                               iintt,liintt, match_mask)
         call bddc_RRT(comm,nnodt,nndft,lnndft,slavery,lslavery,kdoft,lkdoft,ap,lap)
         ap = -ap
         ! Solve the system for given rhs
         call mumps_resolve(schur_mumps,ap,lap)
!*****************************************************************MPI
         call MPI_BCAST(ap,lap, MPI_DOUBLE_PRECISION, 0, comm, ierr)
!*****************************************************************MPI
         where(iintt.eq.0) p = ap
! debug
! Check that the vector is really energy - minimal
!         call sm_vec_mult(matrixtype, nnz, i_sparse, j_sparse, a_sparse, la, &
!                          p,lp, ap,lap)
!      if (any (abs(ap).gt.10e-16 * normrhs .and. iintt.eq.0)) then
!         write(*,*) 'Energy minimality check failed!'
!         write(*,*) 'Sum of error entries',sum(abs(ap),mask=iintt.eq.0)
!         stop
!      end if
!      write(*,*) 'myid = ',myid,'p after preconditioner'
!      write(*,'(f15.9)') p(ihntn)
      end if


! Control of positive definiteness of preconditioner matrix
      if (rmr.le.0.0_kr) then
         if (myid.eq.0) then
            write(*,*) 'WARNING: Preconditioner not positive definite! rmr =',rmr
         end if
      end if


! Setting up the properties for decreasing residual
      ndecr   = 0
      lastres = 1.0_kr

! Measure time of all iterations
      call bddc_time_start(comm)
!***********************************************************************
!*************************MAIN LOOP OVER ITERATIONS*********************
!***********************************************************************
      do iter = 1,maxit

! Multiplication of P by local system matrix 
! ap = A*p
         call bddc_time_start(comm)
         if (myid.eq.0) then
            write(*,*) ' Multiplication by system matrix'
            call flush(6)
         end if
         if (iterate_on_reduced) then
            ! A21 * p
            match_mask(1) = .true.
            match_mask(2) = .false.
            call sm_vec_mult_mask(matrixtype,nnz_mult, i_sparse, j_sparse, a_sparse, la, p,lp, ap,lap, &
                                  iintt,liintt, match_mask)
            ! A22 * p
            match_mask(1) = .true.
            match_mask(2) = .true.
            call sm_vec_mult_mask(matrixtype,nnz_mult, i_sparse, j_sparse, a_sparse, la, p,lp, h,lh, &
                                  iintt,liintt, match_mask)
            ap = ap + h
         else
            call sm_vec_mult(matrixtype,nnz_mult, i_sparse, j_sparse, a_sparse, la, p,lp, ap,lap)
         end if

         call bddc_time_end(comm,time)
         if (myid.eq.0.and.timeinfo.ge.2) then
            write(*,*) '===================================='
            write(*,*) 'Time of matrix multiplication = ',time
            write(*,*) '===================================='
         end if

!***************************************************************PARALLEL
! Scalar product of vectors of old search direction and ap
! pap = p*ap 
         pap_loc = dot_product(p,ap)
!***************************************************************************MPI
         call MPI_ALLREDUCE(pap_loc,pap,1,MPI_DOUBLE_PRECISION,MPI_SUM,comm,ierr)
!***************************************************************************MPI

! Make the AP global (AFTER PAP computation!)
         call bddc_RRT(comm,nnodt,nndft,lnndft,slavery,lslavery,kdoft,lkdoft,ap,lap)

! Control of positive definitenes of system matrix
         if (pap.le.0.0_kr) then
            if (myid.eq.0) then
               write(*,*) ' System matrix not positive definite!, pap =',pap
            end if
!            call MPI_ABORT(comm, 78, ierr)
         end if

! Determination of step lenght ALPHA
         alpha = rmr/pap

! Correction of solution and residual
         sol = sol + alpha*p
         res = res - alpha*ap
         if (iterate_on_reduced) then
            where (iintt.eq.0) res = 0.0_kr
         end if

! Determine global norm of RES
! ||res||
         normres  = bddc_normvec(comm,res,lres)

! Evaluation of relative residual
         relres = normres/normrhs
            
! Write residual to screen
         if (myid.eq.0) then
            write(* ,5001) iter, relres
 5001       format(1X,'iteration iter = ',I4,2X,'relres = ',F25.18)
         end if

! Check convergence
!  relres < tol
         if (relres.lt.tol) then
!***************************************************CONDITION NUMBER ESTIMATION
            nw = iter-1
!***************************************************CONDITION NUMBER ESTIMATION
            if (myid.eq.0) then
               write(*,*)'Number of PCG iterations:',iter
            end if
            exit
         end if

! Check of decreasing of residual
! relres < lastres
         if (relres.lt.lastres) then
            ndecr = 0
         else
            ndecr = ndecr + 1
            if (ndecr.ge.ndecrmax) then
               if (myid.eq.0) then
                  write(*,*)'Residual did not decrease for',ndecrmax,' iterations'
               end if
               stop
            end if
         end if
         lastres = relres

! Shift rmr
         rmrold = rmr

! Action of preconditioner M on residual vector RES 
! h = M*res
         if (myid.eq.0) then
            write(*,*) ' Action of preconditioner'
            call flush(6)
         end if
         call bddc_time_start(comm)
         if (use_preconditioner) then
            call bddc_M(myid,comm,iterate_on_transformed,nnodt,nndft,lnndft,slavery,lslavery,kdoft,lkdoft,res,lres,h,lh,rmr)
         else
            call bddc_M_fake(comm,res,lres,h,lh,rmr)
         end if
         if (iterate_on_reduced) then
! Project p onto the space of energy minimal functions over interiors
            ! Multiply the block of the matrix A12 and p
            match_mask(1) = .false.
            match_mask(2) = .true.
            call sm_vec_mult_mask(matrixtype, nnz, i_sparse, j_sparse, a_sparse, la, &
                                h,lh, ap,lap, &
                                iintt,liintt, match_mask)
            call bddc_RRT(comm,nnodt,nndft,lnndft,slavery,lslavery,kdoft,lkdoft,ap,lap)
            ap = -ap
            ! Solve the system for given rhs
            call mumps_resolve(schur_mumps,ap,lap)
!*****************************************************************MPI
            call MPI_BCAST(ap,lap, MPI_DOUBLE_PRECISION, 0, comm, ierr)
!*****************************************************************MPI
            where(iintt.eq.0) h = ap
         end if

         call bddc_time_end(comm,time)
         if (myid.eq.0.and.timeinfo.ge.2) then
            write(*,*) '===================================='
            write(*,*) 'Time of preconditioning = ',time
            write(*,*) '===================================='
         end if

! Control of positive definiteness of preconditioner matrix
         if (rmr.le.0.0_kr) then
            if (myid.eq.0) then
               write(*,*) 'WARNING: Preconditioner not positive definite! rmr =',rmr
            end if
         end if

! Determination of parameter BETA
         beta = rmr/rmrold

! Determination of new step direction
         p = h + beta*p

!***************************************************CONDITION NUMBER ESTIMATION
! Filling matrix for the Lanczos method
         w((iter-1)*(maxit+1) + iter) = w((iter-1)*(maxit+1) + iter) + 1/alpha
         w((iter)*(maxit+1) + iter + 1) = beta/alpha
         w((iter-1)*(maxit+1) + iter + 1) = -sqrt(beta)/alpha
         w((iter)*(maxit+1) + iter) = w((iter-1)*(maxit+1) + iter + 1)
!***************************************************CONDITION NUMBER ESTIMATION

! Check if iterations reached the limit count
         if (iter.eq.maxit) then
            if (myid.eq.0) then
               write(*,*) 'bddcpcg: Iterations reached prescribed limit without reaching precision.'
            end if
            stop
         end if

      end do
!*************************END OF MAIN CYCLE OVER ITERATIONS*************
! End of mesuring time
      call bddc_time_end(comm,time)
      if (myid.eq.0.and.timeinfo.ge.1) then
         write(*,*) '===================================='
         write(*,*) 'Time of all iterations = ',time
         write(*,*) 'Time per 1 iteration = ',time/float(iter)
         write(*,*) '===================================='
         call flush(6)
      end if

!***************************************************CONDITION NUMBER ESTIMATION
! Condition number estimation on root processor
      if (myid.eq.0) then
         write(*,*) '================================================'
         write(*,*) 'ESTIMATION OF CONDITION NUMBER BY LANCZOS METHOD'
         nwx = maxit + 1
         call condtri(nw,nwx,w,lw, cond)
         write(*,*) 'Condition number cond = ',cond
         write(*,*) '================================================'
      end if
!***************************************************CONDITION NUMBER ESTIMATION

! Transform solution back to original variables
      if (iterate_on_transformed) then
         call bddc_T_apply(comm,sol,lsol)
      end if


! Zero initial solution
   77 continue

! End of mesuring time
      call bddc_time_end(comm,time)
      if (myid.eq.0.and.timeinfo.ge.1) then
         write(*,*) '===================================='
         write(*,*) 'Time of PCG routine = ',time
         write(*,*) '===================================='
      end if

! Clean the memory
      call bddc_finalize
!***************************************************CONDITION NUMBER ESTIMATION
      deallocate(w,wchol,x,y)
!***************************************************CONDITION NUMBER ESTIMATION
      deallocate(p,ap,h)

      return
      end subroutine


 !******************************************************************************************************
subroutine bddcminres(myid,comm,iterate_on_reduced,iterate_on_transformed,matrixtype,nnz,a_sparse,i_sparse,j_sparse,la, &
                      mumpsinfo,timeinfo,use_preconditioner,weight_approach,averages_approach,ndim,nglb,inglb,linglb,nnglb,lnnglb,&
                      nnodt,ndoft,ihntn,lihntn,slavery,lslavery,kdoft,lkdoft,nndft,lnndft, iintt, liintt,&
                      idtr, schur_mumps, res,lres,tol,maxit,ndecrmax,sol,lsol)
!******************************************************************************************************
! Preconditioned minimal residual solver based on BDDC
!******************************************************************************************************
! Use module for BDDC
      use module_bddc1

! Use module for sparse matrices
      use module_sm

! Use module MUMPS
      use module_mumps

! Use module with utilities
      use module_utils

      implicit none
      include "mpif.h"

! Setting of real type kr
      integer,parameter:: kr = kind(1.D0)

! Iterate on staticaly condensed problem ?
      logical,intent(in) :: iterate_on_reduced

! Iterate on transformed problem
      logical,intent(in) :: iterate_on_transformed

! Matrix in sparse IJA format
      integer,intent(in) :: matrixtype
      integer,intent(in) :: nnz, la
      integer,intent(inout) :: i_sparse(la),j_sparse(la)
      real(kr),intent(inout):: a_sparse(la)

! Verbose level of MUMPS
      integer,intent(in) :: mumpsinfo

! Verbose level of times
      integer,intent(in) :: timeinfo

! Use preconditioner?
      logical,intent(in) :: use_preconditioner

! Approach to averaging operator D_P
      integer,intent(in) :: weight_approach

! Description of globs
      integer,intent(in) :: averages_approach
      integer,intent(in) :: ndim
      integer,intent(in) :: nglb
      integer,intent(in) :: linglb,        lnnglb
      integer,intent(in) ::  inglb(linglb), nnglb(lnnglb)

! Space W_tilde
      integer,intent(in) :: nnodt, ndoft
      integer,intent(in) :: lihntn,        lslavery,          lkdoft,        lnndft,        liintt
      integer,intent(in) ::  ihntn(lihntn), slavery(lslavery), kdoft(lkdoft), nndft(lnndft), iintt(liintt)

! Variables for communication
      integer,intent(in):: comm, myid

! Disk unit for transformation
      integer,intent(in):: idtr

! Use this structure of MUMPS for routines from mumps
      type(DMUMPS_STRUC),intent(inout) :: schur_mumps
      
! Right hand side
      integer,intent(in) ::    lres
      real(kr),intent(inout)::  res(lres)
! Solution
      integer,intent(in)::    lsol
      real(kr),intent(inout):: sol(lsol)

! Tolerance
      real(kr),intent(in) :: tol
! Maximum number of iterations
      integer,intent(in) :: maxit
! Maximum number of iterations without residual decrease
      integer,intent(in) :: ndecrmax
      
! Local variables of PCG algorithm
      real(kr):: normres, normrhs, relres, lastres, rmr, cold, c, cnew, alpha0, alpha1, alpha2, alpha3, eta,&
                 gammanew, gamma, gammaold, delta_loc, delta, snew, s, sold
      integer::              lvnew,   lv,   lvold,   lwnew,   lw,   lwold,   lznew,   lz,   lap
      real(kr), allocatable:: vnew(:), v(:), vold(:), wnew(:), w(:), wold(:), znew(:), z(:), ap(:)

      integer:: ierr, iter, ndecr, nnz_mult
      real(kr):: time
      logical:: nonzero_initial_solution_loc = .false. , nonzero_initial_solution
      integer:: nnz_transform

      logical :: match_mask(2)

! Beginning of mesuring time
      call bddc_time_start(comm)
      
! Prepare fields for PMINRES - residual will be associated with RHS
      lvnew = ndoft
      lv    = ndoft
      lvold = ndoft
      lwnew = ndoft
      lw    = ndoft
      lwold = ndoft
      lz    = ndoft
      lznew = ndoft
      lap   = ndoft
      allocate(vnew(lvnew),v(lv),vold(lvold),wnew(lwnew),w(lw),wold(lwold),z(lz),znew(lznew),ap(lap))
      vold = 0._kr
      w    = 0._kr
      wold = 0._kr


!****************************************************************************BDDC
! Beginning of mesuring time
      call bddc_time_start(comm)
! Initialize an instance of the MUMPS package
      call bddc_init(myid,comm,matrixtype,mumpsinfo,timeinfo,use_preconditioner,iterate_on_transformed,&
                     ndoft,nnz,i_sparse,j_sparse,a_sparse,la, &
                     weight_approach,averages_approach,ndim,nglb,inglb,linglb,nnglb,lnnglb,&
                     nnodt,nndft,lnndft,ihntn,lihntn,slavery,lslavery,kdoft,lkdoft, &
                     sol,lsol, res,lres, nnz_transform,idtr)
      call bddc_time_end(comm,time)
      if (myid.eq.0.and.timeinfo.ge.1) then
         write(*,*) 'Initialized BDDC'
         write(*,*) '===================================='
         write(*,*) 'Time of initializing of BDDC = ',time
         write(*,*) '===================================='
      end if
!****************************************************************************BDDC

! Iterate on transformed problem?
      if (iterate_on_transformed) then
         nnz_mult = nnz+nnz_transform
      else
         nnz_mult = nnz
      end if

! Prepare initial residual
! res_0 = rhs - A*sol_0
      if (any(sol.ne.0_kr)) then
         if (myid.eq.0) then
            write(*,*) 'Nonzero initial solution -> create initial vector of residual.'
            nonzero_initial_solution_loc = .true.
         end if
      end if
!*****************************************************************MPI
      call MPI_ALLREDUCE(nonzero_initial_solution_loc,nonzero_initial_solution,1, &
                         MPI_LOGICAL,MPI_LOR,comm,ierr)
!*****************************************************************MPI
      if (nonzero_initial_solution) then
         if (iterate_on_reduced) then
            ! A21 * sol
            match_mask(1) = .true.
            match_mask(2) = .false.
            call sm_vec_mult_mask(matrixtype,nnz_mult, i_sparse, j_sparse, a_sparse, la, sol,lsol, ap,lap, &
                                  iintt,liintt,match_mask)
            ! A22 * sol
            match_mask(1) = .true.
            match_mask(2) = .true.
            call sm_vec_mult_mask(matrixtype,nnz_mult, i_sparse, j_sparse, a_sparse, la, sol,lsol, v,lv, &
                                  iintt,liintt,match_mask)
            ap = ap + v
         else
            call sm_vec_mult(matrixtype,nnz_mult, i_sparse, j_sparse, a_sparse, la, sol,lsol, ap,lap)
         end if
         call bddc_RRT(comm,nnodt,nndft,lnndft,slavery,lslavery,kdoft,lkdoft,ap,lap)
         v = res - ap
      else
         v = res
      end if
      if (iterate_on_reduced) then
         if (any(res.ne.0.0_kr .and. iintt.eq.0)) then
            write(*,*) 'WARNING: Initial check of energy minimal function failed.'
         end if
         where(iintt.eq.0) res = 0.0_kr
      end if

      !if (myid.eq.0) then
      !   write(*,*) 'debug: initial v ='
      !   do i = 1,lv
      !      write(*,'(i8, f14.6)') i, v(i)
      !   end do
      !end if


! Check of zero right hand side => all zero solution
      if (.not.any(v.ne.0.0_kr)) then
         if (myid.eq.0) then
            write(*,*) 'all zero RHS => solution without change'
         end if
         goto 77
      end if
     
! Initial action of the preconditioner
! p = M_BDDC*res 
      if (myid.eq.0) then
         write(*,*) ' Initial action of preconditioner'
      end if
      if (iterate_on_reduced) then
      ! check that residual is energy minimal
         if (any(abs(v).gt.1e-15_kr .and. iintt .eq. 0)) then
            write(*,*) 'Residual is not energy minimal!!!!!!!!!!!!!!!!'
            stop
         end if
      end if
      if (use_preconditioner) then
         call bddc_M(myid,comm,iterate_on_transformed,nnodt,nndft,lnndft,slavery,lslavery,kdoft,lkdoft,v,lv,z,lz,rmr)
! if preconditioner is not positive definite, use gamma 0.
         if (rmr .lt. 0._kr) then
            if (myid.eq.0) then
               write(*,*) 'Ignoring preconditioner for indefinitness.'
            end if
            call bddc_M_fake(comm,v,lv,z,lz,rmr)
         end if
      else
         call bddc_M_fake(comm,v,lv,z,lz,rmr)
      end if
      if (iterate_on_reduced) then
! Project p onto the space of energy minimal functions over interiors
         ! Multiply the block of the matrix A12 and p
         match_mask(1) = .false.
         match_mask(2) = .true.
         call sm_vec_mult_mask(matrixtype, nnz, i_sparse, j_sparse, a_sparse, la, &
                               z,lz, ap,lap, &
                               iintt,liintt, match_mask)
         call bddc_RRT(comm,nnodt,nndft,lnndft,slavery,lslavery,kdoft,lkdoft,ap,lap)
         ap = -ap
         ! Solve the system for given rhs
         call mumps_resolve(schur_mumps,ap,lap)
!*****************************************************************MPI
         call MPI_BCAST(ap,lap, MPI_DOUBLE_PRECISION, 0, comm, ierr)
!*****************************************************************MPI
         where(iintt.eq.0) z = ap
! debug
! Check that the vector is really energy - minimal
!         call sm_vec_mult(matrixtype, nnz, i_sparse, j_sparse, a_sparse, la, &
!                          p,lp, ap,lap)
!      if (any (abs(ap).gt.10e-16 * normrhs .and. iintt.eq.0)) then
!         write(*,*) 'Energy minimality check failed!'
!         write(*,*) 'Sum of error entries',sum(abs(ap),mask=iintt.eq.0)
!         stop
!      end if
!      write(*,*) 'myid = ',myid,'p after preconditioner'
!      write(*,'(f15.9)') p(ihntn)
      end if

! determine norm of initial rhs
! ||rhs||_{M^-1}
      !normrhs = bddc_normvec(comm,v,lv)
      normrhs = rmr
      if (myid.eq.0) then
         write(*,*) 'debug: normrhs =', normrhs
      end if


! determine gamma
      gamma = sqrt(rmr)
      if (myid.eq.0) then
         write(*,*) 'debug: gamma  =', gamma
      end if

      eta  = gamma
      sold = 0._kr
      s    = 0._kr
      cold = 1._kr
      c    = 1._kr
      gammaold = 1._kr

! Setting up the properties for decreasing residual
      ndecr   = 0
      lastres = 1.0_kr

! Measure time of all iterations
      call bddc_time_start(comm)
!***********************************************************************
!*************************MAIN LOOP OVER ITERATIONS*********************
!***********************************************************************
      do iter = 1,maxit

         z = z/gamma
         !if (myid.eq.0) then
         !   write(*,*) 'debug: z ='
         !   do i = 1,lz
         !      write(*,'(i8, f14.6)') i, z(i)
         !   end do
         !end if


! Multiplication of P by local system matrix 
! ap = A*p
         call bddc_time_start(comm)
         if (myid.eq.0) then
            write(*,*) ' Multiplication by system matrix'
         end if
         if (iterate_on_reduced) then
            ! A21 * z
            match_mask(1) = .true.
            match_mask(2) = .false.
            call sm_vec_mult_mask(matrixtype,nnz_mult, i_sparse, j_sparse, a_sparse, la, z,lz, ap,lap, &
                                  iintt,liintt, match_mask)
            ! A22 * z
            match_mask(1) = .true.
            match_mask(2) = .true.
            call sm_vec_mult_mask(matrixtype,nnz_mult, i_sparse, j_sparse, a_sparse, la, z,lz, vnew,lvnew, &
                                  iintt,liintt, match_mask)
            ap = ap + vnew 
         else
            call sm_vec_mult(matrixtype,nnz_mult, i_sparse, j_sparse, a_sparse, la, z,lz, ap,lap)
         end if
         call bddc_time_end(comm,time)
         if (myid.eq.0.and.timeinfo.ge.2) then
            write(*,*) '===================================='
            write(*,*) 'Time of matrix multiplication = ',time
            write(*,*) '===================================='
         end if

!***************************************************************PARALLEL

! Scalar product of vectors of old search direction and ap
! pap = p*ap 
         delta_loc = dot_product(z,ap)
!***************************************************************************MPI
         call MPI_ALLREDUCE(delta_loc,delta,1,MPI_DOUBLE_PRECISION,MPI_SUM,comm,ierr)
!***************************************************************************MPI

! Make the AP global (AFTER PAP computation!)
         call bddc_RRT(comm,nnodt,nndft,lnndft,slavery,lslavery,kdoft,lkdoft,ap,lap)
         !if (myid.eq.0) then
         !   write(*,*) 'debug: ap ='
         !   do i = 1,lap
         !      write(*,'(i8, f14.6)') i, ap(i)
         !   end do
         !end if

! Control of positive definitenes of system matrix
!         if (delta.le.0.0D0) then
!            if (myid.eq.0) then
!               write(*,*) ' System matrix not positive definite!, delta =',delta
!            end if
!!            call MPI_ABORT(comm, 78, ierr)
!         end if
         if (myid.eq.0) then
            write(*,*) 'debug: delta  =', delta
         end if

         vnew = ap - delta/gamma*v - gamma/gammaold*vold
         if (iterate_on_reduced) then
            where (iintt.eq.0) vnew = 0.0_kr
         end if
         !if (myid.eq.0) then
         !   write(*,*) 'debug: vnew ='
         !   do i = 1,lvnew
         !      write(*,'(i8, f14.6)') i, vnew(i)
         !   end do
         !end if
         !


! Action of preconditioner M on residual vector VNEW 
! znew = M*vnew
         if (myid.eq.0) then
            write(*,*) ' Action of preconditioner'
         end if
         call bddc_time_start(comm)
         if (use_preconditioner) then
            call bddc_M(myid,comm,iterate_on_transformed,nnodt,nndft,lnndft,slavery,lslavery,kdoft,lkdoft,vnew,lvnew,znew,lznew,rmr)
            if (rmr .lt. 0._kr) then
               if (myid.eq.0) then
                  write(*,*) 'Ignoring preconditioner for indefinitness.'
               end if
               call bddc_M_fake(comm,vnew,lvnew,znew,lznew,rmr)
            end if
         else
            call bddc_M_fake(comm,vnew,lvnew,znew,lznew,rmr)
         end if
         if (iterate_on_reduced) then
! Project p onto the space of energy minimal functions over interiors
            ! Multiply the block of the matrix A12 and p
            match_mask(1) = .false.
            match_mask(2) = .true.
            call sm_vec_mult_mask(matrixtype, nnz, i_sparse, j_sparse, a_sparse, la, &
                                znew,lznew, ap,lap, &
                                iintt,liintt, match_mask)
            call bddc_RRT(comm,nnodt,nndft,lnndft,slavery,lslavery,kdoft,lkdoft,ap,lap)
            ap = -ap
            ! Solve the system for given rhs
            call mumps_resolve(schur_mumps,ap,lap)
!*****************************************************************MPI
            call MPI_BCAST(ap,lap, MPI_DOUBLE_PRECISION, 0, comm, ierr)
!*****************************************************************MPI
            where(iintt.eq.0) znew = ap
         end if

         call bddc_time_end(comm,time)
         if (myid.eq.0.and.timeinfo.ge.2) then
            write(*,*) '===================================='
            write(*,*) 'Time of preconditioning = ',time
            write(*,*) '===================================='
         end if


         gammanew = sqrt(rmr)

         alpha0 = c*delta - cold*s*gamma
         alpha1 = sqrt(alpha0**2 + gammanew**2)
         alpha2 = s*delta + cold*c*gamma
         alpha3 = sold*gamma
         cnew   = alpha0/alpha1
         snew   = gammanew/alpha1

         wnew   = (z - alpha3*wold - alpha2*w)/alpha1
         !if (myid.eq.0) then
         !   write(*,*) 'debug: wnew ='
         !   do i = 1,lwnew
         !      write(*,'(i8, f14.6)') i, wnew(i)
         !   end do
         !end if
         sol    = sol + cnew*eta*wnew

         eta    = -snew*eta

! Determine global norm of RES
! ||res||_{M^-1}
         !normres  = bddc_normvec(comm,vnew,lvnew)
         normres  = rmr
         if (myid.eq.0) then
            write(*,*) 'debug: normres =', normres
         end if

! Evaluation of relative residual
         relres = normres/normrhs
            
! Print residual to screen
         if (myid.eq.0) then
            write(* ,5001) iter, relres
 5001       format(1X,'iteration iter = ',I4,2X,'relres = ',F25.18)
         end if

! Check convergence
!  relres < tol
         if (relres.lt.tol) then
            if (myid.eq.0) then
               write(*,*)'Number of MINRES iterations:',iter
            end if
            exit
         end if

! Check of decreasing of residual
! relres < lastres
         if (relres.lt.lastres) then
            ndecr = 0
         else
            ndecr = ndecr + 1
            if (ndecr.ge.ndecrmax) then
               if (myid.eq.0) then
                  write(*,*)'Residual did not decrease for',ndecrmax,' iterations'
               end if
               stop
            end if
         end if
         lastres = relres

! Shifts
         vold = v
         v    = vnew

         z    = znew

         wold = w
         w    = wnew

         cold = c
         c    = cnew

         gammaold = gamma
         gamma    = gammanew

         sold     = s
         s        = snew


! Check if iterations reached the limit count
         if (iter.eq.maxit) then
            if (myid.eq.0) then
               write(*,*) 'bddcminres: Iterations reached prescribed limit without reaching precision.'
            end if
            stop
         end if

      end do
!*************************END OF MAIN CYCLE OVER ITERATIONS*************
! End of mesuring time
      call bddc_time_end(comm,time)
      if (myid.eq.0.and.timeinfo.ge.1) then
         write(*,*) '===================================='
         write(*,*) 'Time of all iterations = ',time
         write(*,*) '===================================='
      end if


! Zero initial solution
   77 continue

! End of mesuring time
      call bddc_time_end(comm,time)
      if (myid.eq.0.and.timeinfo.ge.1) then
         write(*,*) '===================================='
         write(*,*) 'Time of MINRES routine = ',time
         write(*,*) '===================================='
      end if

! Clean the memory
      call bddc_finalize

      deallocate(vnew,v,vold,wnew,w,wold,z,znew,ap)

      return
 end subroutine


