!***********************************************************************
program bddcml
!***********************************************************************
! module for distributed Krylov data storage
      use module_krylov_types_def
! module for preconditioner
      use module_levels
! module for domain decomposition data
      use module_dd
! Program name
      use module_utils

      implicit none
      
      include "mpif.h"

!######### PARAMETERS TO SET
! precision of floats
      integer,parameter :: kr = kind(1.D0)

! number of levels
      integer,parameter :: nlevels = 2

! use arithmetic constraints?
      logical,parameter :: use_arithmetic = .true.

! use adaptive constraints?
      logical,parameter :: use_adaptive = .false.

!######### END OF PARAMETERS TO SET

      !  parallel variables
      integer :: myid, comm_all, comm_self, nproc, ierr
      integer :: idpar

      integer :: ndim, nsub, nelem, ndof, nnod, nnodc, linet
      integer :: maxit, ndecrmax
      real(kr) :: tol

      integer :: isub

      integer :: matrixtype

      integer :: nnods, nelems, ndofs, ndofis, nnodis
      integer :: lsolis, lresis
      integer :: lps, laps, lzs

      integer ::             lsols
      real(kr),allocatable :: sols(:)


      character(90)  :: problemname 
      character(100) :: filename

      integer ::                          lpcg_data
      type (pcg_data_type), allocatable :: pcg_data(:)

      logical :: print_solution 

      ! time variables
      real(kr) :: t_total, t_pc_setup, t_pcg


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
         write(*,'(a)') 'MULTILEVEL BDDC solver'
         write(*,'(a)') '======================'

! Name of the problem
   10    write(*,'(a,$)') 'Name of the problem: '
         read(*,*) problemname
         if(problemname.eq.' ') goto 10

      end if
! Broadcast of name of the problem      
!***************************************************************PARALLEL
      call MPI_BCAST(problemname, 90, MPI_CHARACTER, 0, comm_all, ierr)
!***************************************************************PARALLEL

! measuring time
      call time_start(comm_all)

      if (myid.eq.0) then
         filename = trim(problemname)//'.PAR'
         call allocate_unit(idpar)
         open (unit=idpar,file=filename,status='old',form='formatted')
      end if

! Reading basic properties 
      if (myid.eq.0) then
         read(idpar,*) ndim, nsub, nelem, ndof, nnod, nnodc, linet,     &
                       tol, maxit, ndecrmax
         write(*,*)'Characteristics of the problem ',trim(problemname), ':'
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
         call flush(6)
      end if
! Broadcast basic properties of the problem
!***************************************************************PARALLEL
      call MPI_BCAST(ndim,     1,MPI_INTEGER,         0, comm_all, ierr)
      call MPI_BCAST(nsub,     1,MPI_INTEGER,         0, comm_all, ierr)
      call MPI_BCAST(nelem,    1,MPI_INTEGER,         0, comm_all, ierr)
      call MPI_BCAST(ndof,     1,MPI_INTEGER,         0, comm_all, ierr)
      call MPI_BCAST(nnod,     1,MPI_INTEGER,         0, comm_all, ierr)
      call MPI_BCAST(nnodc,    1,MPI_INTEGER,         0, comm_all, ierr)
      call MPI_BCAST(linet,    1,MPI_INTEGER,         0, comm_all, ierr)
      call MPI_BCAST(tol,      1,MPI_DOUBLE_PRECISION,0, comm_all, ierr)
      call MPI_BCAST(maxit,    1,MPI_INTEGER,         0, comm_all, ierr)
      call MPI_BCAST(ndecrmax, 1,MPI_INTEGER,         0, comm_all, ierr)
!***************************************************************PARALLEL

      write (*,*) 'myid = ',myid,': Initializing LEVELS.'
      call levels_init(nlevels,nsub)

! PRECONDITIONER SETUP
      call time_start(comm_all)
      matrixtype = 1 ! SPD matrix
      call levels_pc_setup(problemname,myid,nproc,comm_all,comm_self,matrixtype,ndim,nsub,&
                           use_arithmetic,use_adaptive)
      call time_end(comm_all,t_pc_setup)

! Input zeros as initial vector of solution
      lpcg_data = nsub
      allocate(pcg_data(lpcg_data))
      do isub = 1,nsub

         ! deterine size of subdomain
         call dd_get_size(myid,isub,ndofs,nnods,nelems)
         if (ndofs.ge.0) then
            pcg_data(isub)%is_mine = .true.
         else
            pcg_data(isub)%is_mine = .false.
            cycle
         end if
         call dd_get_interface_size(myid,isub,ndofis,nnodis)

         ! allocate local subdomain solution
         lsols  = ndofs
         allocate(sols(lsols))

         ! set zero initial guess
         ! u_0 = 0
         call zero(sols,lsols)

         ! set initial guess to satisfy Dirichlet boundary conditions
         call dd_fix_bc(myid,isub, sols,lsols)

         ! allocate vectors for Krylov method 
         lsolis = ndofis
         lresis = ndofis
         pcg_data(isub)%lsoli = lsolis
         allocate(pcg_data(isub)%soli(lsolis))
         call zero(pcg_data(isub)%soli,lsolis)
         pcg_data(isub)%lresi = lresis
         allocate(pcg_data(isub)%resi(lresis))
         call zero(pcg_data(isub)%resi,lresis)

         ! restrict solution to interface
         call dd_map_sub_to_subi(myid,isub, sols,lsols, pcg_data(isub)%soli,lsolis)
         deallocate(sols)

         ! set initial residual to RHS
         ! res = g
         call dd_get_reduced_rhs(myid,isub, pcg_data(isub)%resi,lresis)
      end do

      ! Prepare memory for PCG
      do isub = 1,nsub
         if (pcg_data(isub)%is_mine) then

            call dd_get_interface_size(myid,isub,ndofis,nnodis)

            laps = ndofis
            pcg_data(isub)%lap = laps
            allocate(pcg_data(isub)%ap(laps))
            allocate(pcg_data(isub)%apadj(laps))

            lps = ndofis
            pcg_data(isub)%lp = lps
            allocate(pcg_data(isub)%p(lps))
            allocate(pcg_data(isub)%padj(lps))

            lresis = pcg_data(isub)%lresi
            allocate(pcg_data(isub)%resiadj(lresis))

            lzs = ndofis
            pcg_data(isub)%lz = lzs
            allocate(pcg_data(isub)%z(lzs))
            allocate(pcg_data(isub)%zadj(lzs))

         end if
      end do

      ! call PCG method
      call time_start(comm_all)
      call bddcpcg(pcg_data,lpcg_data, nsub, myid,comm_all,tol,maxit,ndecrmax)
      call time_end(comm_all,t_pcg)

      ! Clear memory of PCG
      do isub = 1,nsub
         if (pcg_data(isub)%is_mine) then
            deallocate(pcg_data(isub)%ap)
            deallocate(pcg_data(isub)%apadj)
            deallocate(pcg_data(isub)%p)
            deallocate(pcg_data(isub)%padj)
            deallocate(pcg_data(isub)%resiadj)
            deallocate(pcg_data(isub)%z)
            deallocate(pcg_data(isub)%zadj)
         end if
      end do

      ! Postprocessing of solution - computing interior values
      do isub = 1,nsub
         if (pcg_data(isub)%is_mine) then

            ! determine size of subdomain
            call dd_get_size(myid,isub,ndofs,nnods,nelems)
            call dd_get_interface_size(myid,isub,ndofis,nnodis)
            lsolis = ndofis

            ! allocate local subdomain solution
            lsols  = ndofs
            allocate(sols(lsols))

            ! set zero solution
            call zero(sols,lsols)

            ! resolve interior values
            call dd_resolve_interior(myid,isub, pcg_data(isub)%soli,lsolis, sols,lsols)

            ! write subdomain solution to disk file
            print_solution = .false.
            call dd_write_solution_to_file(myid,problemname,isub,sols,lsols,print_solution)

            deallocate(sols)
         end if
      end do

      ! Clear memory
      do isub = 1,nsub
         if (pcg_data(isub)%is_mine) then
            deallocate(pcg_data(isub)%soli)
            deallocate(pcg_data(isub)%resi)
         end if
      end do
      deallocate(pcg_data)

      write (*,*) 'myid = ',myid,': Finalize LEVELS.'
      call levels_finalize

      call time_end(comm_all,t_total)

      ! Information about times
      if (myid.eq.0) then
         write(*,'(a)')         ' TIMES OF RUN OF ROUTINES:'
         write(*,'(a,f11.3,a)') '  pc_setup  ',t_pc_setup, ' s'
         write(*,'(a,f11.3,a)') '  PCG       ',t_pcg, ' s'
         write(*,'(a)')         '  ______________________'
         write(*,'(a,f11.3,a)') '  total     ',t_total,    ' s'
      end if

      ! MPI finalization
!***************************************************************PARALLEL
      call MPI_FINALIZE(ierr)
!***************************************************************PARALLEL
      end program

!***********************************************************************************
      subroutine bddcpcg(pcg_data,lpcg_data, nsub, myid,comm_all,tol,maxit,ndecrmax)
!***********************************************************************************
! subroutine realizing PCG algorithm with vectors distributed by subdomains

! module for distributed Krylov data storage
      use module_krylov_types_def
! module for preconditioner
      use module_levels
! module for domain decomposition data
      use module_dd
! Program name
      use module_utils

      implicit none
      
      include "mpif.h"

      integer,parameter :: kr = kind(1.D0)

      integer,intent(in) ::                 lpcg_data
      type (pcg_data_type), intent(inout) :: pcg_data(lpcg_data)

      ! number of subdomains
      integer,intent(in) :: nsub

      ! parallel variables
      integer,intent(in) :: myid, comm_all 

      ! limit on iterations
      integer,intent(in) :: maxit

      ! limit on iterations with increasing residual
      integer,intent(in) :: ndecrmax

      ! desired accuracy of relative residual
      real(kr),intent(in) :: tol

      ! local vars
      integer :: isub, i
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

      ! Are you debugging the code?
      logical, parameter :: debug = .true.

      ! Prepare data for Lanczos estimation
      ldiag    = maxit
      lsubdiag = maxit - 1
      allocate(diag(ldiag))
      allocate(subdiag(lsubdiag))
      call zero(diag,ldiag)
      call zero(subdiag,lsubdiag)

      ! get initial residual
      ! r_0 = g - A*u_0
      ! ap = A*u_0
      ! first copy solution to p
      ! p = soli
      do isub = 1,nsub
         if (pcg_data(isub)%is_mine) then
            lsoli = pcg_data(isub)%lsoli
            do i = 1,lsoli
               pcg_data(isub)%p(i) = pcg_data(isub)%soli(i)
            end do  
         end if
      end do

      ! ap = A*u_0
      call bddc_sm_apply(pcg_data,lpcg_data, nsub, myid, comm_all)

      ! update residual
      ! r_0 = g - A*u_0
      do isub = 1,nsub
         if (pcg_data(isub)%is_mine) then

            do i = 1,pcg_data(isub)%lresi
               pcg_data(isub)%resi(i) = pcg_data(isub)%resi(i) - pcg_data(isub)%ap(i)
            end do
            call dd_fix_bc_interface_dual(myid,isub, pcg_data(isub)%resi,pcg_data(isub)%lresi)
         end if
      end do

      ! compute norm of right hand side
      normres2_loc = 0._kr
      do isub = 1,nsub
         if (pcg_data(isub)%is_mine) then

            call dd_dotprod_local(myid,isub,pcg_data(isub)%resi,pcg_data(isub)%lresi, &
                                            pcg_data(isub)%resi,pcg_data(isub)%lresi, &
                                            normres2_sub)

            normres2_loc = normres2_loc + normres2_sub
         end if
      end do
!***************************************************************PARALLEL
      call MPI_ALLREDUCE(normres2_loc,normres2, 1, MPI_DOUBLE_PRECISION,&
                         MPI_SUM, comm_all, ierr) 
!***************************************************************PARALLEL
      normrhs = sqrt(normres2)
      if (debug) then
         if (myid.eq.0) then
            write(*,*) 'Norm of the right hand side: ',normrhs
         end if
      end if

      ! Check of zero right hand side => all zero solution
      if (normrhs.eq.0.0D0) then
         if (myid.eq.0) then
            write(*,*) 'initial residual zero => initial solution exact'
         end if
         return 
      end if


! Initial action of the preconditioner M on residual vector RESI
! M*resi => p
      if (myid.eq.0) then
         write(*,*) ' Initial action of preconditioner'
      end if
      call levels_pc_apply(myid,comm_all, pcg_data,lpcg_data)
      ! produced new z

      ! write z
      !do isub = 1,nsub
      !   if (pcg_data(isub)%is_mine) then
      !      write(*,*) 'myid',myid,'z', pcg_data(isub)%z(1:pcg_data(isub)%lz) 
      !   end if
      !end do

      ! compute rmp = res'*M*res
      ! ||f||
      rmp_loc = 0._kr
      do isub = 1,nsub
         if (pcg_data(isub)%is_mine) then

            call dd_dotprod_local(myid,isub,pcg_data(isub)%resi,pcg_data(isub)%lresi, &
                                            pcg_data(isub)%z,pcg_data(isub)%lz, &
                                            rmp_sub)

            rmp_loc = rmp_loc + rmp_sub
         end if
      end do
!***************************************************************PARALLEL
      call MPI_ALLREDUCE(rmp_loc,rmp, 1, MPI_DOUBLE_PRECISION,          &
                         MPI_SUM, comm_all, ierr) 
!***************************************************************PARALLEL

! Control of positive definiteness of preconditioner matrix
      if (rmp.le.0._kr) then
         if (myid.eq.0) then
            write(*,*) 'Preconditioner not positive definite!'
            call error_exit
         end if
      end if

      if (debug) then
         if (myid.eq.0) then
            write(*,*) 'rmp initial: ',rmp
         end if
      end if

! copy z to p
      ! p = z
      do isub = 1,nsub
         if (pcg_data(isub)%is_mine) then
            lz = pcg_data(isub)%lz
            do i = 1,lz
               pcg_data(isub)%p(i) = pcg_data(isub)%z(i)
            end do
         end if
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
         call bddc_sm_apply(pcg_data,lpcg_data, nsub, myid, comm_all)

         ! write ap
         !do isub = 1,nsub
         !   if (pcg_data(isub)%is_mine) then
         !      write(*,*) 'myid',myid,'ap', pcg_data(isub)%ap(1:pcg_data(isub)%lap) 
         !   end if
         !end do

         ! Scalar product of vectors of old search direction and ap - p*ap => pap
         pap_loc = 0._kr
         do isub = 1,nsub
            if (pcg_data(isub)%is_mine) then

               call dd_dotprod_local(myid,isub,pcg_data(isub)%p,pcg_data(isub)%lp, &
                                               pcg_data(isub)%ap,pcg_data(isub)%lap, &
                                               pap_sub)

               pap_loc = pap_loc + pap_sub
            end if
         end do
!***************************************************************PARALLEL
         call MPI_ALLREDUCE(pap_loc,pap, 1, MPI_DOUBLE_PRECISION,          &
                            MPI_SUM, comm_all, ierr) 
!***************************************************************PARALLEL

! Control of positive definiteness of system matrix
         if (pap.le.0._kr) then
            if (myid.eq.0) then
               write(*,*) ' System matrix not positive definite!'
               call error_exit
            end if
         end if

         if (debug) then
            if (myid.eq.0) then
               write(*,*) 'iter',iter,'pap : ',pap
            end if
         end if

         ! Determination of step lenght ALPHA
         alpha = rmp/pap

         ! Correction of solution vector SOLI and residual vector RES
         ! u   = u   + alpha*p
         ! res = res - alpha*ap
         do isub = 1,nsub
            if (pcg_data(isub)%is_mine) then
               lsoli = pcg_data(isub)%lsoli
               do i = 1,lsoli
                  pcg_data(isub)%soli(i) = pcg_data(isub)%soli(i) + alpha * pcg_data(isub)%p(i)
                  pcg_data(isub)%resi(i) = pcg_data(isub)%resi(i) - alpha * pcg_data(isub)%ap(i)
               end do
            end if
         end do

         ! determine norm of residual 
         ! normres = ||resi||
         normres2_loc = 0._kr
         do isub = 1,nsub
            if (pcg_data(isub)%is_mine) then

               call dd_dotprod_local(myid,isub,pcg_data(isub)%resi,pcg_data(isub)%lresi, &
                                               pcg_data(isub)%resi,pcg_data(isub)%lresi, &
                                               normres2_sub)

               normres2_loc = normres2_loc + normres2_sub
            end if
         end do
!***************************************************************PARALLEL
         call MPI_ALLREDUCE(normres2_loc,normres2, 1, MPI_DOUBLE_PRECISION, &
                            MPI_SUM, comm_all, ierr) 
!***************************************************************PARALLEL
         normres = sqrt(normres2)
         if (debug) then
            if (myid.eq.0) then
               write(*,*) 'iter',iter,'normres : ',normres
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
               write(*,*)'Number of PCG iterations:',iter
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
                  write(*,*)'Residual did not decrease for',ndecrmax,   ' iterations'
                  call error_exit
               end if
            end if
         end if
         lastres = relres

! Action of the preconditioner M on residual vector RES 
! M*resi => z
         if (myid.eq.0) then
            write(*,*) ' Action of preconditioner'
         end if

         call levels_pc_apply(myid,comm_all, pcg_data,lpcg_data)
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
         do isub = 1,nsub
            if (pcg_data(isub)%is_mine) then

               call dd_dotprod_local(myid,isub,pcg_data(isub)%resi,pcg_data(isub)%lresi, &
                                               pcg_data(isub)%z,pcg_data(isub)%lz, &
                                               rmp_sub)

               rmp_loc = rmp_loc + rmp_sub
            end if
         end do
!***************************************************************PARALLEL
         call MPI_ALLREDUCE(rmp_loc,rmp, 1, MPI_DOUBLE_PRECISION,          &
                            MPI_SUM, comm_all, ierr) 
!***************************************************************PARALLEL

         ! Check of positive definiteness of preconditioner matrix
         if (rmp.le.0._kr) then
            if (myid.eq.0) then
               write(*,*) 'Preconditioner not positive definite!'
               call error_exit
            end if
         end if

         if (debug) then
            if (myid.eq.0) then
               write(*,*) 'iter',iter,'rmp : ',rmp
            end if
         end if

         ! Determination of parameter BETA
         beta = rmp/rmpold

         ! Determination of new step direction P
         ! p = z + beta*p
         do isub = 1,nsub
            if (pcg_data(isub)%is_mine) then
               lp = pcg_data(isub)%lp
               do i = 1,lp
                  pcg_data(isub)%p(i) = pcg_data(isub)%z(i) + beta * pcg_data(isub)%p(i)
               end do
            end if
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

!*****************************************************************************
      subroutine bddc_sm_apply(krylov_data,lkrylov_data, nsub, myid, comm_all)
!*****************************************************************************
! subroutine for multiplication of vector krylov_data(isub)%p by system
! matrix S to produce vector krylov_data(isub)%ap
! ap = S*p

! module for distributed Krylov data storage
      use module_krylov_types_def
! module for domain decomposition data
      use module_dd
! Program name
      use module_utils

      implicit none
      
      include "mpif.h"

      integer,parameter :: kr = kind(1.D0)

      integer,intent(in) ::                 lkrylov_data
      type (pcg_data_type), intent(inout) :: krylov_data(lkrylov_data)

      ! number of subdomains
      integer,intent(in) :: nsub

      ! parallel variables
      integer,intent(in) :: myid, comm_all 

      ! local vars
      integer :: isub, i

      ! ap = A*p
      ! Upload data
      do isub = 1,nsub
         if (krylov_data(isub)%is_mine) then

            call zero(krylov_data(isub)%ap,krylov_data(isub)%lap)

            call dd_multiply_by_schur(myid,isub,&
                                      krylov_data(isub)%p,krylov_data(isub)%lp, &
                                      krylov_data(isub)%ap,krylov_data(isub)%lap)

            call dd_comm_upload(myid, isub,  krylov_data(isub)%ap,krylov_data(isub)%lap)
         end if
      end do
      ! Interchange data
      call dd_comm_swapdata(myid, nsub, comm_all)
      ! Download data
      do isub = 1,nsub
         if (krylov_data(isub)%is_mine) then

            call zero(krylov_data(isub)%apadj,krylov_data(isub)%lap)
            ! get contibution from neigbours
            call dd_comm_download(myid, isub,  krylov_data(isub)%apadj,krylov_data(isub)%lap)
            ! join data
            do i = 1,krylov_data(isub)%lap
               krylov_data(isub)%ap(i) = krylov_data(isub)%ap(i) + krylov_data(isub)%apadj(i)
            end do
         end if
      end do

      end subroutine

