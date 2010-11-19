!***********************************************************************
program bddcml
!***********************************************************************
! module for distributed Krylov data storage
      use module_krylov_types_def
! module for preprocessing
      use module_pp
! module for preconditioner
      use module_levels
! Program name
      use module_utils

      implicit none
      
      include "mpif.h"

!######### PARAMETERS TO SET
! precision of floats
      integer,parameter :: kr = kind(1.D0)

! use arithmetic constraints?
      logical,parameter :: use_arithmetic = .true.

! use adaptive constraints?
      logical,parameter :: use_adaptive = .false.

! these options set following type of constraints
!----------------------------------------------------- 
!   \ use_arithmetic |      TRUE     |     FALSE     |
! use_adaptive \     |               |               |
!----------------------------------------------------|
!    TRUE            | edges: arith. | edges: -      |
!                    | faces: adapt. | faces: adapt. |
!----------------------------------------------------|
!    FALSE           | edges: arith. | edges: -      |
!                    | faces: arith. | faces: -      |
!----------------------------------------------------- 

! use prepared division into subdomains on first level in file *.ES?
      logical,parameter :: load_division = .true.
! use prepared selection of corners in file *.CN and description of globs for first level in file *.GLB?
      logical,parameter :: load_globs = .false.
! use prepared file with pairs for adaptivity (*.PAIR) on first level?
      logical,parameter :: load_pairs = .false.
! should parallel division be used (ParMETIS instead of METIS)?
      logical,parameter :: parallel_division = .true.
! correct disconnected subdomains to make them continuous (not allowed for parallel divisions and loaded divisions)
      logical,parameter :: correct_division = .true. 
! should parallel search of neighbours be used? (distributed graph rather than serial graph)
      logical,parameter :: parallel_neighbouring = .true.
! should parallel search of globs be used? (some corrections on globs may not be available)
      logical,parameter :: parallel_globs = .true.
! are you debugging the code?
      logical,parameter :: debug = .true.
! maximal length of problemname
      integer,parameter:: lproblemnamex = 100
! maximal length of any used file - should be reasonably larger than length of problem to allow suffices
      integer,parameter:: lfilenamex = 130
! print solution on screen?
      logical,parameter :: print_solution = .false.
! write solution to a single file instead of distributed files?
      logical,parameter :: write_solution_by_root = .true.

!######### END OF PARAMETERS TO SET
      character(*),parameter:: routine_name = 'BDDCML'

      !  parallel variables
      integer :: myid, comm_all, comm_self, nproc, ierr
      integer :: idpar, idml, idgmi, idfvs, idrhs

      integer :: ndim, nsub, nelem, ndof, nnod, linet
      integer :: maxit, ndecrmax
      real(kr) :: tol

! number of levels
      integer :: nlevels
! subdomains in levels
      integer ::            lnsublev
      integer,allocatable :: nsublev(:)
      integer :: nsub_loc ! number of locally stored subdomains
      integer :: ilevel
      integer :: iaux

      integer ::                    lnnet,   lnndf
      integer,allocatable:: inet(:), nnet(:), nndf(:)
      integer ::           lxyz1,   lxyz2
      real(kr),allocatable:: xyz(:,:)
      integer ::            lifix
      integer,allocatable::  ifix(:)
      integer ::            lfixv
      real(kr),allocatable:: fixv(:)
      integer ::            lrhs
      real(kr),allocatable:: rhs(:)
      integer ::            lsol
      real(kr),allocatable:: sol(:)

      integer :: matrixtype

      integer :: lproblemname
      integer :: meshdim
      integer :: neighbouring

      character(lproblemnamex) :: problemname 
      character(lfilenamex)    :: filename

      integer ::                          lpcg_data
      type (pcg_data_type), allocatable :: pcg_data(:)

      ! time variables
      real(kr) :: t_total, t_import, t_distribute, t_init, t_pc_setup, t_pcg, t_postproc



      ! MPI initialization
!***************************************************************PARALLEL
      call MPI_INIT(ierr)
      ! Communicator
      comm_all  = MPI_COMM_WORLD
      comm_self = MPI_COMM_SELF
      call MPI_COMM_RANK(comm_all,myid,ierr)
      call MPI_COMM_SIZE(comm_all,nproc,ierr)
!***************************************************************PARALLEL

! Initial screen
      if (myid.eq.0) then
         write(*,'(a)') ' _____  ____   ____    ____  __   __ _      '
         write(*,'(a)') '|  _  \|  _  \|  _  \ / __ \|  \ /  | |     '
         write(*,'(a)') '| |_|  | | \  | | \  | /  \_|   V   | |     '
         write(*,'(a)') '|  ___/| |  | | |  | | |    | |\ /| | |     '
         write(*,'(a)') '|  _  \| |  | | |  | | |   _| | V | | |     '
         write(*,'(a)') '| |_|  | |_/  | |_/  | \__/ | |   | | |____ '
         write(*,'(a)') '|_____/|_____/|_____/ \____/|_|   |_|______|'
         write(*,'(a)') '===========multilevel BDDC solver==========='
      end if

! Name of the problem
      call pp_pget_problem_name(comm_all,problemname,lproblemnamex,lproblemname)

! measuring time
      call MPI_BARRIER(comm_all,ierr)
      call time_start

      if (myid.eq.0) then
         filename = problemname(1:lproblemname)//'.PAR'
         call allocate_unit(idpar)
         open (unit=idpar,file=filename,status='old',form='formatted')
      end if
      call pp_pread_par_file(comm_all,idpar, ndim, nsub, nelem, ndof, nnod, linet, tol, maxit, ndecrmax, meshdim)
      if (myid.eq.0) then
         close (idpar)
      end if

! Reading basic properties 
      if (myid.eq.0) then
         write(*,*)'Characteristics of the problem ',problemname(1:lproblemname), ':'
         write(*,*)'  number of processors            nproc =',nproc
         write(*,*)'  number of dimensions             ndim =',ndim
         write(*,*)'  number of elements global       nelem =',nelem
         write(*,*)'  number of DOF                    ndof =',ndof
         write(*,*)'  number of nodes global           nnod =',nnod
         write(*,*)'  lenght of field INET            linet =',linet
         write(*,*)'  mesh dimension                meshdim =',meshdim
         write(*,*)'Characteristics of iterational process:'
         write(*,*)'  tolerance of error                tol =',tol
         write(*,*)'  maximum number of iterations    maxit =',maxit
         write(*,*)'  number of incresing residual ndecrmax =',ndecrmax
         call flush(6)
      end if

      if (myid.eq.0) then
         filename = problemname(1:lproblemname)//'.ML'
         call allocate_unit(idml)
         open (unit=idml,file=filename,status='old',form='formatted')
         rewind idml
         read(idml,*) nlevels
      end if
!***************************************************************PARALLEL
      call MPI_BCAST(nlevels,1,MPI_INTEGER, 0, comm_all, ierr)
!***************************************************************PARALLEL
      lnsublev = nlevels
      allocate(nsublev(lnsublev))
      if (myid.eq.0) then
         read(idml,*) nsublev
      end if
!***************************************************************PARALLEL
      call MPI_BCAST(nsublev,lnsublev,MPI_INTEGER, 0, comm_all, ierr)
!***************************************************************PARALLEL
      if (myid.eq.0) then
         close (idml)
      end if
      if (myid.eq.0) then
         write(*,*)'  number of levels              nlevels =',nlevels
         write(*,*)'  number of subdomains in levels        =',nsublev
         call flush(6)
      end if
      if (nsub.ne.nsublev(1)) then
         if (myid.eq.0) then
            call error(routine_name,'number of subdomains at first level mismatch')
         end if
      end if

      call time_start
      lnnet = nelem
      lnndf = nnod
      lxyz1 = nnod
      lxyz2 = ndim
      allocate(inet(linet),nnet(lnnet),nndf(lnndf),xyz(lxyz1,lxyz2))
      lifix = ndof
      lfixv = ndof
      allocate(ifix(lifix),fixv(lfixv))
      lrhs = ndof
      allocate(rhs(lrhs))
      if (myid.eq.0) then
         write (*,'(a,$)') 'Reading data files ...'
         call flush(6)
      ! read PMD mesh 
         filename = problemname(1:lproblemname)//'.GMIS'
         call allocate_unit(idgmi)
         open (unit=idgmi,file=filename,status='old',form='formatted')
         call pp_read_pmd_mesh(idgmi,inet,linet,nnet,lnnet,nndf,lnndf,xyz,lxyz1,lxyz2)
         close(idgmi)

      ! read PMD boundary conditions 
         filename = problemname(1:lproblemname)//'.FVS'
         call allocate_unit(idfvs)
         open (unit=idfvs,file=filename,status='old',form='formatted')
         call pp_read_pmd_bc(idfvs,ifix,lifix,fixv,lfixv)
         close(idfvs)

      ! read PMD right-hand side
         filename = problemname(1:lproblemname)//'.RHS'
         call allocate_unit(idrhs)
         open (unit=idrhs,file=filename,status='old',form='unformatted')
         call pp_read_pmd_rhs(idrhs,rhs,lrhs)
         close(idrhs)
         write (*,'(a)') 'done.'
         call flush(6)
      end if
      call MPI_BARRIER(comm_all,ierr)
      call time_end(t_import)

      call time_start
      if (myid.eq.0) then
         write (*,'(a,$)') 'Distributing initial data ...'
         call flush(6)
      end if
!***************************************************************PARALLEL
      call MPI_BCAST(inet, linet, MPI_INTEGER, 0, comm_all, ierr)
      call MPI_BCAST(nnet, lnnet, MPI_INTEGER, 0, comm_all, ierr)
      call MPI_BCAST(nndf, lnndf, MPI_INTEGER, 0, comm_all, ierr)
      call MPI_BCAST(xyz, lxyz1*lxyz2, MPI_DOUBLE_PRECISION, 0, comm_all, ierr)
      call MPI_BCAST(ifix, lifix, MPI_INTEGER, 0, comm_all, ierr)
      call MPI_BCAST(fixv, lfixv, MPI_DOUBLE_PRECISION, 0, comm_all, ierr)
      call MPI_BCAST(rhs, lrhs, MPI_DOUBLE_PRECISION, 0, comm_all, ierr)
!***************************************************************PARALLEL

      if (myid.eq.0) then
         write (*,'(a)') 'done.'
         call flush(6)
      end if
      call MPI_BARRIER(comm_all,ierr)
      call time_end(t_distribute)

      ! prepare initial solution
      lsol = ndof
      allocate(sol(lsol))
      sol = 0._kr

      if (myid.eq.0) then
         write (*,'(a)') 'Initializing LEVELS ...'
         call flush(6)
      end if

      call time_start
      call levels_init(nlevels,nsublev,lnsublev,comm_all)
      call levels_load_global_data(nelem,nnod,ndof,&
                                   inet,linet,nnet,lnnet,nndf,lnndf,xyz,lxyz1,lxyz2,&
                                   ifix,lifix,fixv,lfixv,rhs,lrhs,sol,lsol)
      deallocate(inet,nnet,nndf,xyz)
      deallocate(ifix,fixv)
      deallocate(rhs)
      deallocate(sol)
      call MPI_BARRIER(comm_all,ierr)
      call time_end(t_init)
      if (myid.eq.0) then
         write (*,'(a)') 'Initializing LEVELS done.'
         call flush(6)
      end if

      if (myid.eq.0) then
         write (*,'(a)') 'Minimal number of shared nodes to call elements adjacent: '
         call flush(6)
         read (*,*) neighbouring
      end if
! Broadcast basic properties of the problem
!***************************************************************PARALLEL
      call MPI_BCAST(neighbouring,1,MPI_INTEGER,      0, comm_all, ierr)
!***************************************************************PARALLEL

      if (myid.eq.0) then
         write (*,'(a)') 'Preconditioner SETUP ...'
         call flush(6)
      end if
! PRECONDITIONER SETUP
      call MPI_BARRIER(comm_all,ierr)
      call time_start
      matrixtype = 1 ! SPD matrix
      call levels_pc_setup(problemname(1:lproblemname),load_division,load_globs,load_pairs,&
                           parallel_division,correct_division,parallel_neighbouring,neighbouring,&
                           parallel_globs,matrixtype,ndim,meshdim,use_arithmetic,use_adaptive)

      call MPI_BARRIER(comm_all,ierr)
      call time_end(t_pc_setup)

      if (myid.eq.0) then
         write (*,'(a)') 'Preconditioner SETUP done.'
         call flush(6)
      end if

      ilevel = 1
      call levels_get_number_of_subdomains(ilevel,iaux,nsub_loc)
      lpcg_data = nsub_loc
      allocate(pcg_data(lpcg_data))
      ! prepare data and memory for PCG
      call levels_prepare_krylov_data(pcg_data,lpcg_data)

      ! call PCG method
      call MPI_BARRIER(comm_all,ierr)
      call time_start
      call bddcpcg(pcg_data,lpcg_data, nsub_loc, myid,comm_all,tol,maxit,ndecrmax,debug)
      call MPI_BARRIER(comm_all,ierr)
      call time_end(t_pcg)

      ! Postprocessing of solution - computing interior values
      call time_start
      call levels_postprocess_solution(pcg_data,lpcg_data,problemname(1:lproblemname),&
                                       print_solution,write_solution_by_root)
      call time_end(t_postproc)

      ! Clear memory of PCG
      call levels_destroy_krylov_data(pcg_data,lpcg_data)
      deallocate(pcg_data)

      deallocate(nsublev)

      if (myid.eq.0) then
         write (*,'(a)') 'Finalizing LEVELS ...'
         call flush(6)
      end if
      call levels_finalize
      if (myid.eq.0) then
         write (*,'(a)') 'Finalizing LEVELS done.'
         call flush(6)
      end if

      call MPI_BARRIER(comm_all,ierr)
      call time_end(t_total)

      ! Information about times
      if (myid.eq.0) then
         write(*,'(a)')         ' TIMES OF RUN OF ROUTINES:'
         write(*,'(a,f11.3,a)') '  reading data      ',t_import, ' s'
         write(*,'(a,f11.3,a)') '  distributing data ',t_distribute, ' s'
         write(*,'(a,f11.3,a)') '  initialization    ',t_init, ' s'
         write(*,'(a,f11.3,a)') '  pc_setup          ',t_pc_setup, ' s'
         write(*,'(a,f11.3,a)') '  PCG               ',t_pcg, ' s'
         write(*,'(a,f11.3,a)') '  postprocessing    ',t_postproc, ' s'
         write(*,'(a)')         '  ______________________________'
         write(*,'(a,f11.3,a)') '  total             ',t_total,    ' s'
      end if

      ! MPI finalization
!***************************************************************PARALLEL
      call MPI_FINALIZE(ierr)
!***************************************************************PARALLEL
      end program

!*********************************************************************************************
      subroutine bddcpcg(pcg_data,lpcg_data, nsub_loc, myid,comm_all,tol,maxit,ndecrmax,debug)
!*********************************************************************************************
! subroutine realizing PCG algorithm with vectors distributed by subdomains

! module for distributed Krylov data storage
      use module_krylov_types_def
! module for preconditioner
      use module_levels
! Program name
      use module_utils

      implicit none
      
      include "mpif.h"

      integer,parameter :: kr = kind(1.D0)

      integer,intent(in) ::                 lpcg_data
      type (pcg_data_type), intent(inout) :: pcg_data(lpcg_data)

      ! number of locally stored subdomains on first level
      integer,intent(in) :: nsub_loc

      ! parallel variables
      integer,intent(in) :: myid, comm_all 

      ! limit on iterations
      integer,intent(in) :: maxit

      ! limit on iterations with increasing residual
      integer,intent(in) :: ndecrmax

      ! desired accuracy of relative residual
      real(kr),intent(in) :: tol

      ! Are you debugging the code?
      logical, intent(in) :: debug

      ! local vars
      character(*),parameter:: routine_name = 'BDDCPCG'
      integer,parameter :: ilevel = 1
      integer :: isub_loc, i
      integer :: iter, ndecr
      integer :: lz, lsoli, lp

      ! PCG vars
      real(kr) :: normrhs, normres2, normres, normres2_loc, normres2_sub
      real(kr) :: rmp, rmp_loc, rmp_sub
      real(kr) :: pap, pap_loc, pap_sub
      real(kr) :: rmpold
      real(kr) :: alpha, beta
      real(kr) :: relres, lastres

      ! MPI vars
      integer :: ierr

      ! Condition number estimation
      real(kr),allocatable :: diag(:)
      real(kr),allocatable :: subdiag(:)
      integer :: nw, ldiag, lsubdiag
      real(kr) :: cond

      ! Prepare data for Lanczos estimation
      ldiag    = maxit + 1
      lsubdiag = maxit 
      allocate(diag(ldiag))
      allocate(subdiag(lsubdiag))
      call zero(diag,ldiag)
      call zero(subdiag,lsubdiag)

      ! get initial residual
      ! r_0 = g - A*u_0
      ! ap = A*u_0
      ! first copy solution to p
      ! p = soli
      do isub_loc = 1,nsub_loc
         lsoli = pcg_data(isub_loc)%lsoli
         do i = 1,lsoli
            pcg_data(isub_loc)%p(i) = pcg_data(isub_loc)%soli(i)
         end do
      end do

      ! ap = A*u_0
      call levels_sm_apply(pcg_data,lpcg_data)

      ! update residual
      ! r_0 = g - A*u_0
      do isub_loc = 1,nsub_loc
         do i = 1,pcg_data(isub_loc)%lresi
            pcg_data(isub_loc)%resi(i) = pcg_data(isub_loc)%resi(i) - pcg_data(isub_loc)%ap(i)
         end do
      end do
      call levels_fix_bc_interface_dual(pcg_data,lpcg_data)

      ! compute norm of right hand side
      normres2_loc = 0._kr
      do isub_loc = 1,nsub_loc
         call levels_dd_dotprod_local(ilevel,isub_loc,pcg_data(isub_loc)%resi,pcg_data(isub_loc)%lresi, &
                                      pcg_data(isub_loc)%resi,pcg_data(isub_loc)%lresi, &
                                      normres2_sub)
         normres2_loc = normres2_loc + normres2_sub
      end do
!***************************************************************PARALLEL
      call MPI_ALLREDUCE(normres2_loc,normres2, 1, MPI_DOUBLE_PRECISION,&
                         MPI_SUM, comm_all, ierr) 
!***************************************************************PARALLEL
      normrhs = sqrt(normres2)
      if (debug) then
         if (myid.eq.0) then
            call info(routine_name,'Norm of the right hand side =',normrhs)
         end if
      end if

      ! Check of zero right hand side => all zero solution
      if (normrhs.eq.0.0D0) then
         if (myid.eq.0) then
            call warning(routine_name,'initial residual zero => initial solution exact')
         end if
         return 
      end if


! Initial action of the preconditioner M on residual vector RESI
! M*resi => p
      if (debug) then
         if (myid.eq.0) then
            call info(routine_name,' Initial action of preconditioner')
         end if
      end if
      call levels_pc_apply(pcg_data,lpcg_data)
      ! produced new z

      ! write z
      !do isub_loc = 1,nsub_loc
      !   write(*,*) 'myid',myid,'z', pcg_data(isub_loc)%z(1:pcg_data(isub_loc)%lz) 
      !end do

      ! compute rmp = res'*M*res
      ! ||f||
      rmp_loc = 0._kr
      do isub_loc = 1,nsub_loc
         call levels_dd_dotprod_local(ilevel,isub_loc, &
                                      pcg_data(isub_loc)%resi,pcg_data(isub_loc)%lresi, &
                                      pcg_data(isub_loc)%z,pcg_data(isub_loc)%lz, &
                                      rmp_sub)
         rmp_loc = rmp_loc + rmp_sub
      end do
!***************************************************************PARALLEL
      call MPI_ALLREDUCE(rmp_loc,rmp, 1, MPI_DOUBLE_PRECISION,          &
                         MPI_SUM, comm_all, ierr) 
!***************************************************************PARALLEL

! Control of positive definiteness of preconditioner matrix
      if (rmp.le.0._kr) then
         if (myid.eq.0) then
            call error(routine_name,'Preconditioner not positive definite!')
         end if
      end if

      if (debug) then
         if (myid.eq.0) then
            call info(routine_name,'rmp initial =',rmp)
         end if
      end if

! copy z to p
      ! p = z
      do isub_loc = 1,nsub_loc
         lz = pcg_data(isub_loc)%lz
         do i = 1,lz
            pcg_data(isub_loc)%p(i) = pcg_data(isub_loc)%z(i)
         end do
      end do


! Setting up the properties for decreasing residual
      ndecr   = 0
      lastres = 1.0D0
 
!***********************************************************************
!*************************MAIN LOOP OVER ITERATIONS*********************
!***********************************************************************
      do iter = 1,maxit

         ! multiply by system matrix
         ! ap = A * p
         if (debug) then
            if (myid.eq.0) then
               call info(routine_name,' Action of system matrix')
            end if
         end if
         call levels_sm_apply(pcg_data,lpcg_data)

         ! write ap
         !do isub_loc = 1,nsub_loc
         !   write(*,*) 'myid',myid,'ap', pcg_data(isub_loc)%ap(1:pcg_data(isub_loc)%lap) 
         !end do

         ! Scalar product of vectors of old search direction and ap - p*ap => pap
         pap_loc = 0._kr
         do isub_loc = 1,nsub_loc
            call levels_dd_dotprod_local(ilevel,isub_loc, &
                                         pcg_data(isub_loc)%p,pcg_data(isub_loc)%lp, &
                                         pcg_data(isub_loc)%ap,pcg_data(isub_loc)%lap, &
                                         pap_sub)
            pap_loc = pap_loc + pap_sub
         end do
!***************************************************************PARALLEL
         call MPI_ALLREDUCE(pap_loc,pap, 1, MPI_DOUBLE_PRECISION,          &
                            MPI_SUM, comm_all, ierr) 
!***************************************************************PARALLEL

! Control of positive definiteness of system matrix
         if (pap.le.0._kr) then
            if (myid.eq.0) then
               call error(routine_name,'System matrix not positive definite!')
            end if
         end if

         if (debug) then
            if (myid.eq.0) then
               call info(routine_name,'pap =',pap)
            end if
         end if

         ! Determination of step lenght ALPHA
         alpha = rmp/pap

         ! Correction of solution vector SOLI and residual vector RES
         ! u   = u   + alpha*p
         ! res = res - alpha*ap
         do isub_loc = 1,nsub_loc
            lsoli = pcg_data(isub_loc)%lsoli
            do i = 1,lsoli
               pcg_data(isub_loc)%soli(i) = pcg_data(isub_loc)%soli(i) + alpha * pcg_data(isub_loc)%p(i)
               pcg_data(isub_loc)%resi(i) = pcg_data(isub_loc)%resi(i) - alpha * pcg_data(isub_loc)%ap(i)
            end do
         end do

         ! determine norm of residual 
         ! normres = ||resi||
         normres2_loc = 0._kr
         do isub_loc = 1,nsub_loc
            call levels_dd_dotprod_local(ilevel,isub_loc, &
                                         pcg_data(isub_loc)%resi,pcg_data(isub_loc)%lresi, &
                                         pcg_data(isub_loc)%resi,pcg_data(isub_loc)%lresi, &
                                         normres2_sub)
            normres2_loc = normres2_loc + normres2_sub
         end do
!***************************************************************PARALLEL
         call MPI_ALLREDUCE(normres2_loc,normres2, 1, MPI_DOUBLE_PRECISION, &
                            MPI_SUM, comm_all, ierr) 
!***************************************************************PARALLEL
         normres = sqrt(normres2)
         if (debug) then
            if (myid.eq.0) then
               call info(routine_name,'normres =',normres)
            end if
         end if

         ! Evaluation of stopping criterion
         relres = normres/normrhs
            
         if (myid.eq.0) then
            write(* ,5001) iter, relres
 5001       format(1X,'iteration iter = ',I4,2X,'relres = ',F25.18)
         end if

         if (relres.lt.tol) then
            nw = iter-1
            if (myid.eq.0) then
               call info(routine_name,'Number of PCG iterations:',iter)
            end if
            exit
         end if

         ! Check number of iterations
         if (iter.eq.maxit) then
            nw = iter-1
            if (myid.eq.0) then
               call warning(routine_name,'Maximal number of iterations reached, precision not achieved.')
            end if
            exit
         end if

         ! Check of decreasing of residual
         if (relres.lt.lastres) then
            ndecr = 0
         else
            ndecr = ndecr + 1
            if (ndecr.ge.ndecrmax) then
               if (myid.eq.0) then
                  call error(routine_name,'Residual did not decrease for maximal number of iterations:',ndecrmax)
               end if
            end if
         end if
         lastres = relres

! Action of the preconditioner M on residual vector RES 
! M*resi => z
         if (debug) then
            if (myid.eq.0) then
               call info(routine_name,' Action of preconditioner')
            end if
         end if
         call levels_pc_apply(pcg_data,lpcg_data)
         ! produced new z

         ! write z
         !do isub = 1,nsub
         !   if (pcg_data(isub)%is_mine) then
         !      write(*,*) 'myid',myid,'z', pcg_data(isub)%z(1:pcg_data(isub)%lz) 
         !   end if
         !end do

         ! shift generation of res'*M*res
         rmpold = rmp

         ! compute rmp = res'*M*res
         ! ||f||
         rmp_loc = 0._kr
         do isub_loc = 1,nsub_loc
            call levels_dd_dotprod_local(ilevel,isub_loc, &
                                         pcg_data(isub_loc)%resi,pcg_data(isub_loc)%lresi, &
                                         pcg_data(isub_loc)%z,pcg_data(isub_loc)%lz, &
                                         rmp_sub)
            rmp_loc = rmp_loc + rmp_sub
         end do
!***************************************************************PARALLEL
         call MPI_ALLREDUCE(rmp_loc,rmp, 1, MPI_DOUBLE_PRECISION,          &
                            MPI_SUM, comm_all, ierr) 
!***************************************************************PARALLEL

         ! Check of positive definiteness of preconditioner matrix
         if (rmp.le.0._kr) then
            if (myid.eq.0) then
               call error(routine_name,'Preconditioner not positive definite!')
            end if
         end if

         if (debug) then
            if (myid.eq.0) then
               call info(routine_name,'rmp =',rmp)
            end if
         end if

         ! Determination of parameter BETA
         beta = rmp/rmpold

         ! Determination of new step direction P
         ! p = z + beta*p
         do isub_loc = 1,nsub_loc
            lp = pcg_data(isub_loc)%lp
            do i = 1,lp
               pcg_data(isub_loc)%p(i) = pcg_data(isub_loc)%z(i) + beta * pcg_data(isub_loc)%p(i)
            end do
         end do

         ! Filling matrix for the Lanczos method
         diag(iter) = diag(iter) + 1/alpha
         diag(iter+1) = beta/alpha
         subdiag(iter) = -sqrt(beta)/alpha

      end do
!*************************END OF MAIN LOOP OVER ITERATIONS**************

! Condition number estimation on root processor
      call condsparse(nw,diag,nw,subdiag,nw-1, cond)
      if (myid.eq.0) then
         write(*,*) '================================================'
         write(*,*) 'ESTIMATION OF CONDITION NUMBER BY LANCZOS METHOD'
         write(*,*) 'Condition number cond = ',cond
         write(*,*) '================================================'
      end if
      deallocate(diag)
      deallocate(subdiag)

      end subroutine

