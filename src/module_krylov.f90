! BDDCML - Multilevel BDDC
! Copyright (C) The BDDCML Team
! 
! This program is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
! 
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
! 
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <http://www.gnu.org/licenses/>.
! ____________________________________________________________________

module module_krylov
! module for some Krylov subspace iterative methods suitable for DD implementation

      implicit none
! adjustable parameters ############################
! type of real variables
      integer,parameter,private :: kr = kind(1.D0)
! numerical zero
      real(kr),parameter,private :: numerical_zero = 1.e-12_kr
! debugging 
      logical,parameter,private :: debug = .true.
! profiling 
      logical,parameter,private :: profile = .true.
! adjustable parameters ############################

      contains

!*****************************************************************************
      subroutine krylov_bddcpcg(problemname, comm_all,tol,maxit,ndecrmax,&
                                print_solution, write_solution_by_root)
!*****************************************************************************
! subroutine realizing PCG algorithm with vectors distributed by subdomains

! module for distributed Krylov data storage
      use module_krylov_types_def
! module for preconditioner
      use module_levels
! Program name
      use module_utils

      implicit none
      
      include "mpif.h"

      ! name of the problem
      character(*),intent(in) :: problemname

      ! parallel variables
      integer,intent(in) :: comm_all 

      ! limit on iterations
      integer,intent(in) :: maxit

      ! limit on iterations with increasing residual
      integer,intent(in) :: ndecrmax

      ! desired accuracy of relative residual
      real(kr),intent(in) :: tol

      ! print solution on screen?
      logical, intent(in) :: print_solution

      ! write solution to a single file instead of distributed files?
      logical, intent(in) :: write_solution_by_root

      ! local vars
      character(*),parameter:: routine_name = 'KRYLOV_BDDCPCG'
      integer,parameter :: ilevel = 1

      ! data for storing actual PCG data
      integer ::                                  lpcg_data
      type (pcg_data_type), allocatable, target :: pcg_data(:)

      ! data for auxiliary manipulation with preconditioner and system matrix 
      integer ::                                     lcommon_krylov_data
      type (common_krylov_data_type), allocatable ::  common_krylov_data(:)

      integer :: myid
      integer :: nsub, nsub_loc
      integer :: isub_loc, i
      integer :: iter, ndecr
      integer :: lsoli, lp
      integer :: ndofis, nnodis

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

      ! time variables
      real(kr) :: t_sm_apply, t_pc_apply
      real(kr) :: t_postproc

      ! orient in the communicator
      call MPI_COMM_RANK(comm_all,myid,ierr)

      ! Prepare data for Lanczos estimation
      ldiag    = maxit + 1
      lsubdiag = maxit 
      allocate(diag(ldiag))
      allocate(subdiag(lsubdiag))
      call zero(diag,ldiag)
      call zero(subdiag,lsubdiag)

      ! prepare data and memory for PCG
      call levels_get_number_of_subdomains(ilevel,nsub,nsub_loc)
      lcommon_krylov_data = nsub_loc
      allocate(common_krylov_data(lcommon_krylov_data))
      lpcg_data = nsub_loc
      allocate(pcg_data(lpcg_data))
      do isub_loc = 1,nsub_loc
         call levels_dd_get_interface_size(ilevel,isub_loc, ndofis, nnodis)
         pcg_data(isub_loc)%lsoli = ndofis
         allocate(pcg_data(isub_loc)%soli(pcg_data(isub_loc)%lsoli))
         pcg_data(isub_loc)%lresi = ndofis
         allocate(pcg_data(isub_loc)%resi(pcg_data(isub_loc)%lresi))
         pcg_data(isub_loc)%lap   = ndofis
         allocate(pcg_data(isub_loc)%ap(pcg_data(isub_loc)%lap))
         pcg_data(isub_loc)%lp    = ndofis
         allocate(pcg_data(isub_loc)%p(pcg_data(isub_loc)%lp))
         pcg_data(isub_loc)%lz    = ndofis
         allocate(pcg_data(isub_loc)%z(pcg_data(isub_loc)%lz))
      end do

      ! prepare initial solution and right-hand side
      do isub_loc = 1,nsub_loc
         call levels_prepare_interface_initial_data(isub_loc,pcg_data(isub_loc)%soli,pcg_data(isub_loc)%lsoli,&
                                                             pcg_data(isub_loc)%resi,pcg_data(isub_loc)%lresi)
         ! fix boundary conditions in residual to zero
         call levels_dd_fix_bc_interface_dual(ilevel,isub_loc,pcg_data(isub_loc)%resi,pcg_data(isub_loc)%lresi)
      end do

      ! get initial residual
      ! r_0 = g - A*u_0
      ! ap = A*u_0
      ! first set pointers to soli and ap
      do isub_loc = 1,nsub_loc
         common_krylov_data(isub_loc)%lvec_in  = pcg_data(isub_loc)%lsoli
         common_krylov_data(isub_loc)%vec_in  => pcg_data(isub_loc)%soli
         common_krylov_data(isub_loc)%lvec_out = pcg_data(isub_loc)%lap
         common_krylov_data(isub_loc)%vec_out => pcg_data(isub_loc)%ap
      end do
      call MPI_BARRIER(comm_all,ierr)
      call time_start
      call levels_sm_apply(common_krylov_data,lcommon_krylov_data)
      call MPI_BARRIER(comm_all,ierr)
      call time_end(t_sm_apply)
      if (myid.eq.0) then
         call time_print('application of system matrix',t_sm_apply)
      end if

      ! update residual
      ! r_0 = g - A*u_0
      do isub_loc = 1,nsub_loc
         do i = 1,pcg_data(isub_loc)%lresi
            pcg_data(isub_loc)%resi(i) = pcg_data(isub_loc)%resi(i) - pcg_data(isub_loc)%ap(i)
         end do
      end do
      ! fix boundary conditions in residual to zero
      do isub_loc = 1,nsub_loc
         call levels_dd_fix_bc_interface_dual(ilevel,isub_loc,pcg_data(isub_loc)%resi,pcg_data(isub_loc)%lresi)
      end do

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
      ! first set pointers to resi and p
      do isub_loc = 1,nsub_loc
         common_krylov_data(isub_loc)%lvec_in  = pcg_data(isub_loc)%lresi
         common_krylov_data(isub_loc)%vec_in  => pcg_data(isub_loc)%resi
         common_krylov_data(isub_loc)%lvec_out = pcg_data(isub_loc)%lp
         common_krylov_data(isub_loc)%vec_out => pcg_data(isub_loc)%p
      end do
      call MPI_BARRIER(comm_all,ierr)
      call time_start
      call levels_pc_apply(common_krylov_data,lcommon_krylov_data)
      call MPI_BARRIER(comm_all,ierr)
      call time_end(t_pc_apply)
      if (myid.eq.0) then
         call time_print('application of preconditioner',t_pc_apply)
      end if
      ! produced new p

      ! compute rmp = res'*M*res
      ! ||f||
      rmp_loc = 0._kr
      do isub_loc = 1,nsub_loc
         call levels_dd_dotprod_local(ilevel,isub_loc, &
                                      pcg_data(isub_loc)%resi,pcg_data(isub_loc)%lresi, &
                                      pcg_data(isub_loc)%p,pcg_data(isub_loc)%lp, &
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
            call warning(routine_name,'Preconditioner not positive definite!')
         end if
      end if

      if (debug) then
         if (myid.eq.0) then
            call info(routine_name,'rmp initial =',rmp)
         end if
      end if

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
         ! first set pointers to soli 
         do isub_loc = 1,nsub_loc
            common_krylov_data(isub_loc)%lvec_in  = pcg_data(isub_loc)%lp
            common_krylov_data(isub_loc)%vec_in  => pcg_data(isub_loc)%p
            common_krylov_data(isub_loc)%lvec_out = pcg_data(isub_loc)%lap
            common_krylov_data(isub_loc)%vec_out => pcg_data(isub_loc)%ap
         end do
         call levels_sm_apply(common_krylov_data,lcommon_krylov_data)

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
               call warning(routine_name,'System matrix not positive definite!')
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
            write(*,'(a,i5,a,f25.18)') 'iteration: ',iter,', relative residual: ',relres
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
         ! first set pointers to resi and p
         do isub_loc = 1,nsub_loc
            common_krylov_data(isub_loc)%lvec_in  = pcg_data(isub_loc)%lresi
            common_krylov_data(isub_loc)%vec_in  => pcg_data(isub_loc)%resi
            common_krylov_data(isub_loc)%lvec_out = pcg_data(isub_loc)%lz
            common_krylov_data(isub_loc)%vec_out => pcg_data(isub_loc)%z
         end do
         call levels_pc_apply(common_krylov_data,lcommon_krylov_data)
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
               call warning(routine_name,'Preconditioner not positive definite!')
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

      ! Postprocessing of solution - computing interior values
      call time_start
      ! first set pointers to soli
      do isub_loc = 1,nsub_loc
         common_krylov_data(isub_loc)%lvec_in  = pcg_data(isub_loc)%lsoli
         common_krylov_data(isub_loc)%vec_in  => pcg_data(isub_loc)%soli
      end do
      call levels_postprocess_solution(common_krylov_data,lcommon_krylov_data,problemname,&
                                       print_solution,write_solution_by_root)
      call time_end(t_postproc)
      if (myid.eq.0) then
         call time_print('postprocessing of solution',t_postproc)
      end if

      ! Clear memory of PCG
      do isub_loc = 1,nsub_loc
         nullify(common_krylov_data(isub_loc)%vec_in)
         nullify(common_krylov_data(isub_loc)%vec_out)
      end do
      deallocate(common_krylov_data)
      do isub_loc = 1,nsub_loc
         deallocate(pcg_data(isub_loc)%soli)
         deallocate(pcg_data(isub_loc)%resi)
         deallocate(pcg_data(isub_loc)%ap)
         deallocate(pcg_data(isub_loc)%p)
         deallocate(pcg_data(isub_loc)%z)
      end do
      deallocate(pcg_data)

      end subroutine

!****************************************************************************
      subroutine krylov_bddcbicgstab(problemname, comm_all,tol,maxit,ndecrmax,&
                                     print_solution, write_solution_by_root)
!****************************************************************************
! subroutine realizing BICGSTAB algorithm with vectors distributed by subdomains

! module for distributed Krylov data storage
      use module_krylov_types_def
! module for preconditioner
      use module_levels
! Program name
      use module_utils

      implicit none
      
      include "mpif.h"

      ! name of the problem
      character(*),intent(in) :: problemname

      ! parallel variables
      integer,intent(in) :: comm_all 

      ! limit on iterations
      integer,intent(in) :: maxit

      ! limit on iterations with increasing residual
      integer,intent(in) :: ndecrmax

      ! desired accuracy of relative residual
      real(kr),intent(in) :: tol

      ! print solution on screen?
      logical, intent(in) :: print_solution

      ! write solution to a single file instead of distributed files?
      logical, intent(in) :: write_solution_by_root

      ! local vars
      character(*),parameter:: routine_name = 'KRYLOV_BDDCBICGSTAB'
      integer,parameter :: ilevel = 1

      ! data for storing actual BICGSTAB data
      integer ::                                       lbicgstab_data
      type (bicgstab_data_type), allocatable, target :: bicgstab_data(:)

      ! data for auxiliary manipulation with preconditioner and system matrix 
      integer ::                                     lcommon_krylov_data
      type (common_krylov_data_type), allocatable ::  common_krylov_data(:)

      integer :: myid
      integer :: nsub, nsub_loc
      integer :: isub_loc, i
      integer :: iter, ndecr
      integer :: ndofis, nnodis
      integer :: lsoli, lp

      ! BICGSTAB vars
      real(kr) :: normrhs, normres2, normres, normres2_loc, normres2_sub
      real(kr) :: tt, tt_loc, tt_sub
      real(kr) :: ts, ts_loc, ts_sub
      real(kr) :: vrstab, vrstab_loc, vrstab_sub
      real(kr) :: rho, rhoold, rho_loc, rho_sub
      real(kr) :: alpha, beta
      real(kr) :: omega
      real(kr) :: relres, lastres

      ! MPI vars
      integer :: ierr

      ! time variables
      real(kr) :: t_postproc

      ! orient in the communicator
      call MPI_COMM_RANK(comm_all,myid,ierr)

      ! find number of subdomains
      call levels_get_number_of_subdomains(ilevel,nsub,nsub_loc)

      ! prepare data and memory for BICGSTAB
      lcommon_krylov_data = nsub_loc
      allocate(common_krylov_data(lcommon_krylov_data))
      lbicgstab_data = nsub_loc
      allocate(bicgstab_data(lbicgstab_data))
      do isub_loc = 1,nsub_loc
         call levels_dd_get_interface_size(ilevel,isub_loc, ndofis, nnodis)
         bicgstab_data(isub_loc)%lsoli = ndofis
         allocate(bicgstab_data(isub_loc)%soli(bicgstab_data(isub_loc)%lsoli))
         bicgstab_data(isub_loc)%lresi = ndofis
         allocate(bicgstab_data(isub_loc)%resi(bicgstab_data(isub_loc)%lresi))
         bicgstab_data(isub_loc)%lresistab = ndofis
         allocate(bicgstab_data(isub_loc)%resistab(bicgstab_data(isub_loc)%lresistab))
         bicgstab_data(isub_loc)%lv = ndofis
         allocate(bicgstab_data(isub_loc)%v(bicgstab_data(isub_loc)%lv))
         bicgstab_data(isub_loc)%lp    = ndofis
         allocate(bicgstab_data(isub_loc)%p(bicgstab_data(isub_loc)%lp))
         bicgstab_data(isub_loc)%ly    = ndofis
         allocate(bicgstab_data(isub_loc)%y(bicgstab_data(isub_loc)%ly))
         bicgstab_data(isub_loc)%lz    = ndofis
         allocate(bicgstab_data(isub_loc)%z(bicgstab_data(isub_loc)%lz))
         bicgstab_data(isub_loc)%ls    = ndofis
         allocate(bicgstab_data(isub_loc)%s(bicgstab_data(isub_loc)%ls))
         bicgstab_data(isub_loc)%lt    = ndofis
         allocate(bicgstab_data(isub_loc)%t(bicgstab_data(isub_loc)%lt))
      end do

      do isub_loc = 1,nsub_loc
         ! prepare initial solution and right-hand side
         call levels_prepare_interface_initial_data(isub_loc,bicgstab_data(isub_loc)%soli,bicgstab_data(isub_loc)%lsoli,&
                                                             bicgstab_data(isub_loc)%resi,bicgstab_data(isub_loc)%lresi)
         ! fix boundary conditions in residual to zero
         call levels_dd_fix_bc_interface_dual(ilevel,isub_loc,bicgstab_data(isub_loc)%resi,bicgstab_data(isub_loc)%lresi)
      end do

      ! get initial residual
      ! r_0 = g - A*u_0
      ! ap = A*u_0
      ! first set pointers to soli and ap
      do isub_loc = 1,nsub_loc
         common_krylov_data(isub_loc)%lvec_in  = bicgstab_data(isub_loc)%lsoli
         common_krylov_data(isub_loc)%vec_in  => bicgstab_data(isub_loc)%soli
         common_krylov_data(isub_loc)%lvec_out = bicgstab_data(isub_loc)%lv
         common_krylov_data(isub_loc)%vec_out => bicgstab_data(isub_loc)%v
      end do
      call levels_sm_apply(common_krylov_data,lcommon_krylov_data)

      ! update residual
      ! r_0 = g - A*u_0
      do isub_loc = 1,nsub_loc
         do i = 1,bicgstab_data(isub_loc)%lresi
            bicgstab_data(isub_loc)%resi(i) = bicgstab_data(isub_loc)%resi(i) - bicgstab_data(isub_loc)%v(i)
         end do
      end do
      ! fix boundary conditions in residual to zero
      do isub_loc = 1,nsub_loc
         call levels_dd_fix_bc_interface_dual(ilevel,isub_loc,bicgstab_data(isub_loc)%resi,bicgstab_data(isub_loc)%lresi)
      end do

      ! compute norm of right hand side
      normres2_loc = 0._kr
      do isub_loc = 1,nsub_loc
         call levels_dd_dotprod_local(ilevel,isub_loc,bicgstab_data(isub_loc)%resi,bicgstab_data(isub_loc)%lresi, &
                                      bicgstab_data(isub_loc)%resi,bicgstab_data(isub_loc)%lresi, &
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

      ! BICGSTAB initialization
      rhoold = 1._kr
      alpha  = 1._kr
      omega  = 1._kr
      ! v = 0
      do isub_loc = 1,nsub_loc
         call zero(bicgstab_data(isub_loc)%v,bicgstab_data(isub_loc)%lv)
      end do
      ! p = 0
      do isub_loc = 1,nsub_loc
         call zero(bicgstab_data(isub_loc)%p,bicgstab_data(isub_loc)%lp)
      end do
      ! shadow residual
      do isub_loc = 1,nsub_loc
         do i = 1,bicgstab_data(isub_loc)%lresi
            bicgstab_data(isub_loc)%resistab(i) = bicgstab_data(isub_loc)%resi(i)
         end do
      end do

! Setting up the properties for decreasing residual
      ndecr   = 0
      lastres = 1.0D0
 
!***********************************************************************
!*************************MAIN LOOP OVER ITERATIONS*********************
!***********************************************************************
      do iter = 1,maxit

         ! Scalar product of vectors of res and resstab
         ! rho = res*resstab
         rho_loc = 0._kr
         do isub_loc = 1,nsub_loc
            call levels_dd_dotprod_local(ilevel,isub_loc,bicgstab_data(isub_loc)%resi,bicgstab_data(isub_loc)%lresi, &
                                         bicgstab_data(isub_loc)%resistab,bicgstab_data(isub_loc)%lresistab, &
                                         rho_sub)
            rho_loc = rho_loc + rho_sub
         end do
!***************************************************************PARALLEL
         call MPI_ALLREDUCE(rho_loc,rho, 1, MPI_DOUBLE_PRECISION,&
                            MPI_SUM, comm_all, ierr) 
!***************************************************************PARALLEL
         if (debug) then
            if (myid.eq.0) then
               call info(routine_name,'rho =',rho)
            end if
         end if

         beta = rho*alpha/(rhoold*omega)
         if (debug) then
            if (myid.eq.0) then
               call info(routine_name,'beta =',beta)
            end if
         end if

         !p = res + beta*(p - omega*v)
         do isub_loc = 1,nsub_loc
            lp = bicgstab_data(isub_loc)%lp
            do i = 1,lp
               bicgstab_data(isub_loc)%p(i) = bicgstab_data(isub_loc)%resi(i) &
                                            + beta * (bicgstab_data(isub_loc)%p(i) - omega * bicgstab_data(isub_loc)%v(i))
            end do
         end do

         ! Action of preconditioner M on vector P 
         ! y = M*p
         if (debug) then
            if (myid.eq.0) then
               call info(routine_name,' Action of preconditioner')
            end if
         end if
         ! first set properly pointers
         do isub_loc = 1,nsub_loc
            common_krylov_data(isub_loc)%lvec_in  = bicgstab_data(isub_loc)%lp
            common_krylov_data(isub_loc)%vec_in  => bicgstab_data(isub_loc)%p
            common_krylov_data(isub_loc)%lvec_out = bicgstab_data(isub_loc)%ly
            common_krylov_data(isub_loc)%vec_out => bicgstab_data(isub_loc)%y
         end do
         call levels_pc_apply(common_krylov_data,lcommon_krylov_data)

         ! Multiplication of Y by local system matrix 
         ! v = A*y
         if (debug) then
            if (myid.eq.0) then
               call info(routine_name,' Action of system matrix')
            end if
         end if
         ! first set properly pointers
         do isub_loc = 1,nsub_loc
            common_krylov_data(isub_loc)%lvec_in  = bicgstab_data(isub_loc)%ly
            common_krylov_data(isub_loc)%vec_in  => bicgstab_data(isub_loc)%y
            common_krylov_data(isub_loc)%lvec_out = bicgstab_data(isub_loc)%lv
            common_krylov_data(isub_loc)%vec_out => bicgstab_data(isub_loc)%v
         end do
         call levels_sm_apply(common_krylov_data,lcommon_krylov_data)

         ! Scalar product of vectors of v and resstab
         ! vrstab = v*resstab
         vrstab_loc = 0._kr
         do isub_loc = 1,nsub_loc
            call levels_dd_dotprod_local(ilevel,isub_loc,bicgstab_data(isub_loc)%v,bicgstab_data(isub_loc)%lv, &
                                         bicgstab_data(isub_loc)%resistab,bicgstab_data(isub_loc)%lresistab, &
                                         vrstab_sub)
            vrstab_loc = vrstab_loc + vrstab_sub
         end do
!***************************************************************PARALLEL
         call MPI_ALLREDUCE(vrstab_loc,vrstab, 1, MPI_DOUBLE_PRECISION,&
                            MPI_SUM, comm_all, ierr) 
!***************************************************************PARALLEL
         if (debug) then
            if (myid.eq.0) then
               call info(routine_name,'vrstab =',vrstab)
            end if
         end if

         alpha = rho/vrstab
         if (debug) then
            if (myid.eq.0) then
               call info(routine_name,'alpha =',alpha)
            end if
         end if

         ! build half step
         ! soli = soli + alpha * y
         ! s = res - alpha*v
         do isub_loc = 1,nsub_loc
            lsoli = bicgstab_data(isub_loc)%lsoli
            do i = 1,lsoli
               bicgstab_data(isub_loc)%soli(i) = bicgstab_data(isub_loc)%soli(i) &
                                               + alpha * bicgstab_data(isub_loc)%y(i) 
               bicgstab_data(isub_loc)%s(i) = bicgstab_data(isub_loc)%resi(i) &
                                            - alpha * bicgstab_data(isub_loc)%v(i)
            end do
         end do
         ! determine norm of residual 
         ! normres = ||resi||
         normres2_loc = 0._kr
         do isub_loc = 1,nsub_loc
            call levels_dd_dotprod_local(ilevel,isub_loc, &
                                         bicgstab_data(isub_loc)%s,bicgstab_data(isub_loc)%ls, &
                                         bicgstab_data(isub_loc)%s,bicgstab_data(isub_loc)%ls, &
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
               call info(routine_name,'normres half =',normres)
            end if
         end if

         ! Evaluation of stopping criterion
         relres = normres/normrhs

! Print residual to screen
         if (myid.eq.0) then
            write(*,'(a,i5,a,f25.18)') 'iteration: ',iter-1,'.5, relative residual: ',relres
         end if

! Check convergence in the half step
!  relres < tol
         if (relres.lt.tol) then
            if (myid.eq.0) then
               write(*,'(a,a,i5,a)') routine_name,': Number of BICGSTAB iterations: ',iter-1,'.5'
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

         ! Action of preconditioner M on vector S 
         ! z = M*s
         if (debug) then
            if (myid.eq.0) then
               call info(routine_name,' Action of preconditioner')
            end if
         end if
         ! first set properly pointers
         do isub_loc = 1,nsub_loc
            common_krylov_data(isub_loc)%lvec_in  = bicgstab_data(isub_loc)%ls
            common_krylov_data(isub_loc)%vec_in  => bicgstab_data(isub_loc)%s
            common_krylov_data(isub_loc)%lvec_out = bicgstab_data(isub_loc)%lz
            common_krylov_data(isub_loc)%vec_out => bicgstab_data(isub_loc)%z
         end do
         call levels_pc_apply(common_krylov_data,lcommon_krylov_data)

         ! Multiplication of Z by local system matrix 
         ! t = A*z
         if (debug) then
            if (myid.eq.0) then
               call info(routine_name,' Action of system matrix')
            end if
         end if
         ! first set properly pointers
         do isub_loc = 1,nsub_loc
            common_krylov_data(isub_loc)%lvec_in  = bicgstab_data(isub_loc)%lz
            common_krylov_data(isub_loc)%vec_in  => bicgstab_data(isub_loc)%z
            common_krylov_data(isub_loc)%lvec_out = bicgstab_data(isub_loc)%lt
            common_krylov_data(isub_loc)%vec_out => bicgstab_data(isub_loc)%t
         end do
         call levels_sm_apply(common_krylov_data,lcommon_krylov_data)

         ! Scalar product of vectors s and t
         ! ts = s * t
         ts_loc = 0._kr
         do isub_loc = 1,nsub_loc
            call levels_dd_dotprod_local(ilevel,isub_loc,bicgstab_data(isub_loc)%s,bicgstab_data(isub_loc)%ls, &
                                         bicgstab_data(isub_loc)%t,bicgstab_data(isub_loc)%lt, &
                                         ts_sub)
            ts_loc = ts_loc + ts_sub
         end do
!***************************************************************PARALLEL
         call MPI_ALLREDUCE(ts_loc,ts, 1, MPI_DOUBLE_PRECISION,&
                            MPI_SUM, comm_all, ierr) 
!***************************************************************PARALLEL

         ! Scalar product of vectors t and t
         ! tt = t * t
         tt_loc = 0._kr
         do isub_loc = 1,nsub_loc
            call levels_dd_dotprod_local(ilevel,isub_loc,bicgstab_data(isub_loc)%t,bicgstab_data(isub_loc)%lt, &
                                         bicgstab_data(isub_loc)%t,bicgstab_data(isub_loc)%lt, &
                                         tt_sub)
            tt_loc = tt_loc + tt_sub
         end do
!***************************************************************PARALLEL
         call MPI_ALLREDUCE(tt_loc,tt, 1, MPI_DOUBLE_PRECISION,&
                            MPI_SUM, comm_all, ierr) 
!***************************************************************PARALLEL
         if (debug) then
            if (myid.eq.0) then
               call info(routine_name,'tt =',tt)
            end if
         end if

         omega = ts/tt
         if (debug) then
            if (myid.eq.0) then
               call info(routine_name,'omega =',omega)
            end if
         end if

         ! Final correction of solution vector SOLI and residual vector RES
         !soli = soli + omega*z
         !res  = s - omega*t
         do isub_loc = 1,nsub_loc
            lsoli = bicgstab_data(isub_loc)%lsoli
            do i = 1,lsoli
               bicgstab_data(isub_loc)%soli(i) = bicgstab_data(isub_loc)%soli(i) &
                                               + omega * bicgstab_data(isub_loc)%z(i) 
               bicgstab_data(isub_loc)%resi(i) = bicgstab_data(isub_loc)%s(i) - omega * bicgstab_data(isub_loc)%t(i)
            end do
         end do

         ! determine norm of residual 
         ! normres = ||resi||
         normres2_loc = 0._kr
         do isub_loc = 1,nsub_loc
            call levels_dd_dotprod_local(ilevel,isub_loc, &
                                         bicgstab_data(isub_loc)%resi,bicgstab_data(isub_loc)%lresi, &
                                         bicgstab_data(isub_loc)%resi,bicgstab_data(isub_loc)%lresi, &
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
            
! Print residual to screen
         if (myid.eq.0) then
            write(*,'(a,i5,a,f25.18)') 'iteration: ',iter,', relative residual: ',relres
         end if

! Check convergence
!  relres < tol
         if (relres.lt.tol) then
            if (myid.eq.0) then
               write(*,'(a,a,i5)') routine_name,': Number of BICGSTAB iterations: ',iter
            end if
            exit
         end if

         ! Check number of iterations
         if (iter.eq.maxit) then
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

! Shift rho
         rhoold = rho

      end do
!*************************END OF MAIN LOOP OVER ITERATIONS**************

      ! Postprocessing of solution - computing interior values
      call time_start
      ! first set pointers to soli
      do isub_loc = 1,nsub_loc
         common_krylov_data(isub_loc)%lvec_in  = bicgstab_data(isub_loc)%lsoli
         common_krylov_data(isub_loc)%vec_in  => bicgstab_data(isub_loc)%soli
      end do
      call levels_postprocess_solution(common_krylov_data,lcommon_krylov_data,problemname,&
                                       print_solution,write_solution_by_root)
      call time_end(t_postproc)
      if (myid.eq.0) then
         call time_print('postprocessing of solution',t_postproc)
      end if

      ! Clear memory of BICGSTAB
      do isub_loc = 1,nsub_loc
         nullify(common_krylov_data(isub_loc)%vec_in)
         nullify(common_krylov_data(isub_loc)%vec_out)
      end do
      deallocate(common_krylov_data)
      do isub_loc = 1,nsub_loc
         deallocate(bicgstab_data(isub_loc)%soli)
         deallocate(bicgstab_data(isub_loc)%resi)
         deallocate(bicgstab_data(isub_loc)%resistab)
         deallocate(bicgstab_data(isub_loc)%v)
         deallocate(bicgstab_data(isub_loc)%p)
         deallocate(bicgstab_data(isub_loc)%y)
         deallocate(bicgstab_data(isub_loc)%z)
         deallocate(bicgstab_data(isub_loc)%s)
         deallocate(bicgstab_data(isub_loc)%t)
      end do
      deallocate(bicgstab_data)

      end subroutine

end module
