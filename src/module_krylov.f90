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

    module module_krylov
    ! module for some Krylov subspace iterative methods suitable for DD implementation
          use, intrinsic :: iso_fortran_env
    ! module for distributed Krylov data storage
          use module_krylov_types_def

          implicit none
    ! adjustable parameters ############################
    ! type of real variables
          integer,parameter,private :: kr = REAL64
    ! debugging 
          logical,parameter,private :: debug = .false.
    ! profiling 
          logical,private ::           profile = .false.
    ! adjustable parameters ############################

    ! data necessary for recycling of Krylov subspace
          logical,private :: is_recycling_prepared = .false.
          integer,private                                        :: nactive_cols_recycling_basis = 0
          integer,private                                        :: lrecycling_basis = 0
          type(krylov_recycling_data_type), allocatable, private ::  recycling_basis(:) ! one per each local subdomain
          !real(kr),allocatable :: recycling_vtw(:,:) ! matrix of all p'*A*p at interface or its LU factors
          !integer,allocatable :: recycling_ipiv(:)   ! permutation from LAPACK LU
          ! the above matrix should be diagonal, so for a well-conditioned basis, only its diagonal is needed
          integer ::             lrecycling_idvtw
          real(kr),allocatable :: recycling_idvtw(:) ! array of the inverse of the diagonal of p'*A*p at interface 
          !logical,private :: recycling_is_inverse_prepared = .false.

          contains

    !******************************************************************************************************
          subroutine krylov_bddcpcg(comm_all,tol,maxit,ndecrmax, recycling, max_number_of_stored_vectors, &
                                    num_iter, converged_reason, cond)
    !******************************************************************************************************
    ! subroutine realizing PCG algorithm with vectors distributed by subdomains

    ! module for preconditioner
          use module_levels
    ! Program name
          use module_utils

          implicit none
          
          include "mpif.h"

          ! parallel variables
          integer,intent(in) :: comm_all 

          ! limit on iterations
          integer,intent(in) :: maxit

          ! limit on iterations with increasing residual
          integer,intent(in) :: ndecrmax

          ! desired accuracy of relative residual
          real(kr),intent(in) :: tol

          ! should recycling of Krylov space be used?
          logical,intent(in) :: recycling 

          ! if recycling should be used, how many vectors of the Krylov basis do you want to store?
          integer,intent(in) :: max_number_of_stored_vectors

          ! resulting number of iterations
          integer,intent(out) :: num_iter

          ! convergence reason
          !  =  0 - converged relative residual
          !  = -1 - reached limit on number of iterations
          !  = -2 - reached limit on number of iterations with nondecreasing residual
          integer,intent(out) :: converged_reason

          ! estimated condition number
          real(kr),intent(out) :: cond

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
          integer :: isub_loc, i, j
          integer :: iter, ndecr
          integer :: lsoli, lp
          integer :: ndofis, nnodis

          ! PCG vars
          real(kr) :: normrhs, normrhs2, normrhs2_loc, normrhs2_sub
          real(kr) :: normres, normres0, normres2, normres2_loc, normres2_sub
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

          ! Recycling of Krylov spaces
          integer :: jcol, ibasis, jbasis
          integer :: lstored
          integer ::              lvtb
          real(kr), allocatable :: vtb(:)
          real(kr), allocatable :: vtb_loc(:)
          real(kr) :: vtb_sub
          integer ::              lvtw
          real(kr), allocatable :: vtw(:,:)
          real(kr), allocatable :: vtw_loc(:,:)
          real(kr) :: vtw_sub
          real(kr) :: check_orthogonality
          real(kr) :: scaling_of_basis
          !integer :: recycling_info

          ! time variables
          real(kr) :: t_sm_apply, t_pc_apply
          real(kr) :: t_postproc
          real(kr) :: t_krylov_solve
          real(kr) :: t_recycling_projection

    !-----profile
          if (profile) then
             call MPI_BARRIER(comm_all,ierr)
             call time_start
          end if
    !-----profile

          ! orient in the communicator
          call MPI_COMM_RANK(comm_all,myid,ierr)

          ! Prepare data for Lanczos estimation
          ldiag    = maxit + 1
          lsubdiag = maxit 
          allocate(diag(ldiag))
          allocate(subdiag(lsubdiag))
          diag(:) = 0._kr
          subdiag(:) = 0._kr

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

          ! edits for recycling of Krylov space
          if (recycling) then
             ! allocate when called first
             if (.not. is_recycling_prepared) then

                lrecycling_idvtw = max_number_of_stored_vectors
                allocate(recycling_idvtw(lrecycling_idvtw))

                lrecycling_basis = nsub_loc
                allocate(recycling_basis(lrecycling_basis))

                ! prepare arrays to default value
                do isub_loc = 1,nsub_loc
                   call levels_dd_get_interface_size(ilevel,isub_loc, ndofis, nnodis)

                   recycling_basis(isub_loc)%lv1 = ndofis
                   recycling_basis(isub_loc)%lv2 = max_number_of_stored_vectors
                   allocate(recycling_basis(isub_loc)%v(recycling_basis(isub_loc)%lv1,recycling_basis(isub_loc)%lv2))
                   recycling_basis(isub_loc)%lw1 = ndofis
                   recycling_basis(isub_loc)%lw2 = max_number_of_stored_vectors
                   allocate(recycling_basis(isub_loc)%w(recycling_basis(isub_loc)%lw1,recycling_basis(isub_loc)%lw2))
                end do

                is_recycling_prepared = .true.
             else
                ! check the matching size
                if (lrecycling_basis.ne.nsub_loc) then
                   call error( routine_name, 'Size mismatch of array RECYCLING_BASIS' )
                end if
             end if
          end if

          ! prepare initial solution and right-hand side
          do isub_loc = 1,nsub_loc
             call levels_prepare_interface_initial_data(isub_loc,pcg_data(isub_loc)%soli,pcg_data(isub_loc)%lsoli,&
                                                                 pcg_data(isub_loc)%resi,pcg_data(isub_loc)%lresi)
             ! fix boundary conditions in residual to zero
             call levels_dd_fix_bc_interface_dual(ilevel,isub_loc,pcg_data(isub_loc)%resi,pcg_data(isub_loc)%lresi)
          end do

          ! compute norm of right-hand side
          normrhs2_loc = 0._kr
          do isub_loc = 1,nsub_loc
             call levels_dd_dotprod_local(ilevel,isub_loc,pcg_data(isub_loc)%resi,pcg_data(isub_loc)%lresi, &
                                          pcg_data(isub_loc)%resi,pcg_data(isub_loc)%lresi, &
                                          normrhs2_sub)
             normrhs2_loc = normrhs2_loc + normrhs2_sub
          end do
    !***************************************************************PARALLEL
          call MPI_ALLREDUCE(normrhs2_loc,normrhs2, 1, MPI_DOUBLE_PRECISION,&
                             MPI_SUM, comm_all, ierr) 
    !***************************************************************PARALLEL
          normrhs = sqrt(normrhs2)
          if (debug) then
             if (myid.eq.0) then
                call info(routine_name,'Norm of the right-hand side =',normrhs)
             end if
          end if

          ! Check of zero right-hand side => all zero solution
          if (normrhs.eq.0.0D0) then
             if (myid.eq.0) then
                call warning(routine_name,'initial right-hand side zero => zero solution')
             end if
             return 
          end if

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
          if (myid.eq.0 .and. profile) then
             call time_print('application of system matrix',t_sm_apply)
          end if

          do isub_loc = 1,nsub_loc
             ! fix boundary conditions in residual to zero
             call levels_dd_fix_bc_interface_dual(ilevel,isub_loc,pcg_data(isub_loc)%ap,pcg_data(isub_loc)%lap)
          end do

          ! update residual
          ! r_0 = g - A*u_0
          do isub_loc = 1,nsub_loc
             ! fix boundary conditions in residual to zero
             do i = 1,pcg_data(isub_loc)%lresi
                pcg_data(isub_loc)%resi(i) = pcg_data(isub_loc)%resi(i) - pcg_data(isub_loc)%ap(i)
             end do
          end do

          ! compute norm of the initial residual
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
          normres0 = sqrt(normres2)
          if (debug) then
             if (myid.eq.0) then
                call info(routine_name,'Norm of the initial residual =',normres0)
             end if
          end if

          ! Check of zero right-hand side => all zero solution
          if (normres0.eq.0.0D0) then
             if (myid.eq.0) then
                call warning(routine_name,'initial residual zero => initial solution exact')
             end if
             return 
          end if

          if (recycling .and. nactive_cols_recycling_basis .gt. 0) then 

             if (debug) then
                ! check orthogonality of basis
                ! V'*W
                lvtw     = nactive_cols_recycling_basis
                allocate(vtw_loc(lvtw,lvtw))
                vtw_loc = 0._kr
                allocate(vtw(lvtw,lvtw))
                do ibasis = 1,nactive_cols_recycling_basis
                   do jbasis = 1,nactive_cols_recycling_basis
                      do isub_loc = 1,nsub_loc
                         call levels_dd_dotprod_local(ilevel,isub_loc,&
                                                      recycling_basis(isub_loc)%v(:,ibasis),recycling_basis(isub_loc)%lv1, &
                                                      recycling_basis(isub_loc)%w(:,jbasis),recycling_basis(isub_loc)%lw1, &
                                                      vtw_sub)
                         vtw_loc(ibasis,jbasis) = vtw_loc(ibasis,jbasis) + vtw_sub
                      end do
                   end do
                end do
    !***************************************************************PARALLEL
                call MPI_ALLREDUCE(vtw_loc,vtw, lvtw*lvtw, MPI_DOUBLE_PRECISION,&
                                   MPI_SUM, comm_all, ierr) 
    !***************************************************************PARALLEL
                deallocate(vtw_loc)

                !if (myid.eq.0) then
                !   write(*,*) 'V^T*W'
                !   do i = 1,lvtw
                !      write(*,'(1000f12.5)') (vtw(i,j), j = 1,lvtw)
                !   end do
                !end if

                ! inv(D)*V'*W
                do ibasis = 1,nactive_cols_recycling_basis
                   !vtw(:,ibasis) = recycling_idvtw(1:nactive_cols_recycling_basis) * vtw(:,ibasis) 
                   call recycling_apply_diagonal_inverse(vtw(1,ibasis),lvtw)
                end do

                ! subtract identity matrix
                ! inv(D)*V'*W - I
                do ibasis = 1,nactive_cols_recycling_basis
                   vtw(ibasis,ibasis) = vtw(ibasis,ibasis) - 1._kr
                end do
                check_orthogonality = sqrt(sum(vtw**2))
                if (myid.eq.0) then
                   write(*,*) 'D^{-1}*V^T*W - I'
                   do i = 1,lvtw
                      write(*,'(1000e12.5)') (vtw(i,j), j = 1,lvtw)
                   end do
                   if (check_orthogonality .gt. 1.e-4_kr) then
                      call warning( routine_name, 'Krylov basis has problems with orthogonality:', check_orthogonality )
                   end if
                end if
                deallocate(vtw)
             end if 

             ! Initial projection of the right-hand side onto the Krylov basis
             call MPI_BARRIER(comm_all,ierr)
             call time_start

             ! project the right-hand side onto the stored Krylov space
             lvtb     = nactive_cols_recycling_basis
             allocate(vtb_loc(lvtb))
             vtb_loc = 0._kr
             allocate(vtb(lvtb))

             ! V'*b
             do ibasis = 1,nactive_cols_recycling_basis
                do isub_loc = 1,nsub_loc
                   call levels_dd_dotprod_local(ilevel,isub_loc,&
                                                recycling_basis(isub_loc)%v(:,ibasis),recycling_basis(isub_loc)%lv1, &
                                                pcg_data(isub_loc)%resi,pcg_data(isub_loc)%lresi, &
                                                vtb_sub)
                   vtb_loc(ibasis) = vtb_loc(ibasis) + vtb_sub
                end do
             end do
    !***************************************************************PARALLEL
             call MPI_ALLREDUCE(vtb_loc,vtb, lvtb, MPI_DOUBLE_PRECISION,&
                                MPI_SUM, comm_all, ierr) 
    !***************************************************************PARALLEL
             deallocate(vtb_loc)

             ! inv(D)*V'*b
             !vtb(1:nactive_cols_recycling_basis) = vtb(1:nactive_cols_recycling_basis) &
             !                                    * recycling_idvtw(1:nactive_cols_recycling_basis)
             call recycling_apply_diagonal_inverse(vtb,lvtb)

             ! update solution by projection
             ! u <- u + Pb = u + V*inv(D)*V'*b = u + V*IDIAG*V'
             do isub_loc = 1,nsub_loc
                ! u_P - projected part of solution
                pcg_data(isub_loc)%z = matmul(recycling_basis(isub_loc)%v(:,1:nactive_cols_recycling_basis), &
                                              vtb(1:nactive_cols_recycling_basis) )
                ! u_I <- u_I + u_P - projected part of solution
                pcg_data(isub_loc)%soli = pcg_data(isub_loc)%soli + pcg_data(isub_loc)%z
             end do

             deallocate(vtb)

             ! update residual
             ! r_0 = r_0 - A*u_P
             ! ap = A*u_0
             ! first set pointers to soli and ap
             do isub_loc = 1,nsub_loc
                common_krylov_data(isub_loc)%lvec_in  = pcg_data(isub_loc)%lz
                common_krylov_data(isub_loc)%vec_in  => pcg_data(isub_loc)%z
                common_krylov_data(isub_loc)%lvec_out = pcg_data(isub_loc)%lap
                common_krylov_data(isub_loc)%vec_out => pcg_data(isub_loc)%ap
             end do
             call MPI_BARRIER(comm_all,ierr)
             call time_start
             call levels_sm_apply(common_krylov_data,lcommon_krylov_data)
             call MPI_BARRIER(comm_all,ierr)
             call time_end(t_sm_apply)
             if (myid.eq.0 .and. profile) then
                call time_print('application of system matrix',t_sm_apply)
             end if

             ! update residual
             ! r_0 = g - A*u_P
             do isub_loc = 1,nsub_loc
                do i = 1,pcg_data(isub_loc)%lresi
                   pcg_data(isub_loc)%resi(i) = pcg_data(isub_loc)%resi(i) - pcg_data(isub_loc)%ap(i)
                end do
             end do
             ! fix boundary conditions in residual to zero
             do isub_loc = 1,nsub_loc
                call levels_dd_fix_bc_interface_dual(ilevel,isub_loc,pcg_data(isub_loc)%resi,pcg_data(isub_loc)%lresi)
             end do

             ! compute norm of right-hand side
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
             normres = sqrt(normres2)
             if (debug) then
                if (myid.eq.0) then
                   call info(routine_name,'Norm of the first residual =',normres)
                end if
             end if

             call MPI_BARRIER(comm_all,ierr)
             call time_end(t_recycling_projection)
             if (myid.eq.0 .and. profile) then
                call time_print('RHS projection onto existing Krylov basis',t_recycling_projection)
             end if

             ! Evaluation of stopping criterion
             ! compute it as relative residual with respect to initial residual
             relres = normres/normres0
             ! compute it as relative residual with respect to right-hand side
             !relres = normres/normrhs
             if (relres.lt.tol) then
                iter = 0
                if (myid.eq.0) then
                   call info(routine_name,'Number of PCG iterations:',iter)
                end if
                num_iter = iter
                converged_reason = 0
                nw = 0
                goto 123
             end if
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
          if (myid.eq.0 .and. profile) then
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
          !if (rmp.le.0._kr) then
          !   if (myid.eq.0) then
          !      call warning(routine_name,'Preconditioner not positive definite!')
          !   end if
          !end if

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

             if (recycling .and. nactive_cols_recycling_basis.gt.0) then
                ! orthogonalize vector P with respect to the columns the stored Krylov space directions
                do isub_loc = 1,nsub_loc
                   common_krylov_data(isub_loc)%vec_in  => pcg_data(isub_loc)%p
                   common_krylov_data(isub_loc)%lvec_in =  pcg_data(isub_loc)%lp
                end do
                call krylov_orthogonalize_icgs( comm_all, common_krylov_data, lcommon_krylov_data)
                !call krylov_orthogonalize_cgs( comm_all, common_krylov_data, lcommon_krylov_data)
                !call krylov_orthogonalize_mgs( comm_all, common_krylov_data, lcommon_krylov_data)


                if (debug) then
                   ! check orthogonality of basis
                   ! V'*W
                   lvtw     = nactive_cols_recycling_basis
                   allocate(vtw_loc(lvtw,lvtw))
                   vtw_loc = 0._kr
                   allocate(vtw(lvtw,lvtw))
                   do ibasis = 1,nactive_cols_recycling_basis
                      do jbasis = 1,nactive_cols_recycling_basis
                         do isub_loc = 1,nsub_loc
                            call levels_dd_dotprod_local(ilevel,isub_loc,&
                                                         recycling_basis(isub_loc)%v(:,ibasis),recycling_basis(isub_loc)%lv1, &
                                                         recycling_basis(isub_loc)%w(:,jbasis),recycling_basis(isub_loc)%lw1, &
                                                         vtw_sub)
                            vtw_loc(ibasis,jbasis) = vtw_loc(ibasis,jbasis) + vtw_sub
                         end do
                      end do
                   end do
    !***************************************************************PARALLEL
                   call MPI_ALLREDUCE(vtw_loc,vtw, lvtw*lvtw, MPI_DOUBLE_PRECISION,&
                                      MPI_SUM, comm_all, ierr) 
    !***************************************************************PARALLEL
                   deallocate(vtw_loc)

                   if (myid.eq.0) then
                      write(*,*) 'V^T*W'
                      do i = 1,lvtw
                         write(*,'(1000f20.13)') (vtw(i,j), j = 1,lvtw)
                      end do
                   end if

                   ! inv(D)*V'*W
                   do ibasis = 1,nactive_cols_recycling_basis
                      !vtw(:,ibasis) = recycling_idvtw(1:nactive_cols_recycling_basis) * vtw(:,ibasis) 
                      call recycling_apply_diagonal_inverse(vtw(1,ibasis),lvtw)
                   end do

                   ! subtract identity matrix
                   ! inv(D)*V'*W - I
                   do ibasis = 1,nactive_cols_recycling_basis
                      vtw(ibasis,ibasis) = vtw(ibasis,ibasis) - 1._kr
                   end do
                   check_orthogonality = sqrt(sum(vtw**2))
                   if (myid.eq.0) then
                      write(*,*) 'D^{-1}*V^T*W - I'
                      do i = 1,lvtw
                         write(*,'(1000e12.5)') (vtw(i,j), j = 1,lvtw)
                      end do
                      if (check_orthogonality .gt. 1.e-4_kr) then
                         call warning( routine_name, 'Krylov basis has problems with orthogonality:', check_orthogonality )
                      end if
                   end if
                   deallocate(vtw)
                end if 

             end if

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

             if (recycling) then
                !if (nactive_cols_recycling_basis .le.  max_number_of_stored_vectors ) then

                ! if there is no room for storing a new vector, move the basis forward to fit the latest vector
                if (nactive_cols_recycling_basis .eq. max_number_of_stored_vectors ) then
                   lstored = max_number_of_stored_vectors - 1
                   do isub_loc = 1,nsub_loc
                      recycling_basis(isub_loc)%v(:,1:lstored) = recycling_basis(isub_loc)%v(:,2:lstored+1)
                      recycling_basis(isub_loc)%w(:,1:lstored) = recycling_basis(isub_loc)%w(:,2:lstored+1)
                   end do
                   recycling_idvtw(1:lstored) = recycling_idvtw(2:lstored+1)
                   
                   nactive_cols_recycling_basis = nactive_cols_recycling_basis - 1
                end if

                ! scaling factor - make the basis vectors normalized
                scaling_of_basis = 1._kr / sqrt(pap)

                jcol = nactive_cols_recycling_basis + 1
                do isub_loc = 1,nsub_loc
                   ! V <- [ V p ]
                   recycling_basis(isub_loc)%v(:,jcol) = scaling_of_basis * pcg_data(isub_loc)%p(:)
                   ! W <- [ W ap ]
                   recycling_basis(isub_loc)%w(:,jcol) = scaling_of_basis * pcg_data(isub_loc)%ap(:)
                end do
                ! IDIAG <- [ IDIAG 1./p'Ap ]
                !recycling_idvtw(jcol) = 1._kr / pap
                recycling_idvtw(jcol) = 1._kr 

                nactive_cols_recycling_basis = nactive_cols_recycling_basis + 1

                ! update the V'*W matrix and its factors
                !lvtw     = nactive_cols_recycling_basis
                !allocate(vtw_loc(lvtw,lvtw))
                !vtw_loc = 0._kr
                !do ibasis = 1,nactive_cols_recycling_basis
                !   do jbasis = 1,nactive_cols_recycling_basis
                !      do isub_loc = 1,nsub_loc
                !         call levels_dd_dotprod_local(ilevel,isub_loc,&
                !                                      recycling_basis(isub_loc)%v(:,ibasis),recycling_basis(isub_loc)%lv1, &
                !                                      recycling_basis(isub_loc)%w(:,jbasis),recycling_basis(isub_loc)%lw1, &
                !                                      vtw_sub)
                !         vtw_loc(ibasis,jbasis) = vtw_loc(ibasis,jbasis) + vtw_sub
                !      end do
                !   end do
                !end do
                !if (allocated(recycling_vtw)) then
                !   deallocate(recycling_vtw)
                !end if
                !allocate(recycling_vtw(lvtw,lvtw))
    !************************************************************PARALLEL
                !call MPI_ALLREDUCE(vtw_loc,recycling_vtw, lvtw*lvtw, MPI_DOUBLE_PRECISION,&
                !                   MPI_SUM, comm_all, ierr) 
    !************************************************************PARALLEL
                !deallocate(vtw_loc)

                !if (myid.eq.0) then
                !   write(*,*) 'V^T*W'
                !   do i = 1,lvtw
                !      write(*,'(1000f20.13)') (recycling_vtw(i,j), j = 1,lvtw)
                !   end do
                !end if

                ! prepare factors
                !if (allocated(recycling_ipiv)) then
                !   deallocate(recycling_ipiv)
                !end if
                !allocate(recycling_ipiv(lvtw))
                !call DGETRF( lvtw, lvtw, recycling_vtw, lvtw, recycling_ipiv, recycling_info )
                !recycling_is_inverse_prepared = .true.

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
             ! compute it as relative residual with respect to initial residual
             relres = normres/normres0
             ! compute it as relative residual with respect to right-hand side
             !relres = normres/normrhs
                
             if (myid.eq.0) then
                call info (routine_name, 'iteration: ',iter)
                call info (routine_name, '          relative residual: ',relres)
             end if

             if (relres.lt.tol) then
                nw = iter-1
                if (myid.eq.0) then
                   call info(routine_name,'Number of PCG iterations:',iter)
                end if
                num_iter = iter
                converged_reason = 0
                exit
             end if

             ! Check number of iterations
             if (iter.eq.maxit) then
                nw = iter-1
                if (myid.eq.0) then
                   call warning(routine_name,'Maximal number of iterations reached, precision not achieved.')
                end if
                num_iter = iter
                converged_reason = -1
                exit
             end if

             ! Check of decreasing of residual
             if (relres.lt.lastres) then
                ndecr = 0
             else
                ndecr = ndecr + 1
                if (ndecr.ge.ndecrmax) then
                   nw = iter-1
                   if (myid.eq.0) then
                      call warning(routine_name,'Residual did not decrease for maximal number of iterations:',ndecrmax)
                   end if
                   num_iter = iter
                   converged_reason = -2
                   exit
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
             diag(iter) = diag(iter) + 1._kr/alpha
             diag(iter+1) = beta/alpha
             if (beta.ge.0._kr) then
                subdiag(iter) = -sqrt(beta)/alpha
             else
                subdiag(iter) = 0._kr
             end if

          end do
    !*************************END OF MAIN LOOP OVER ITERATIONS**************
    123   continue

    ! Condition number estimation on root processor, if there are no NaNs
          if (nw.gt.0) then
             call krylov_condsparse(nw,diag,nw,subdiag,nw-1, cond)
          else
             cond = 1._kr
          end if
          if (myid.eq.0) then
             call info(routine_name, '================================================')
             call info(routine_name, 'ESTIMATION OF CONDITION NUMBER BY LANCZOS METHOD')
             call info(routine_name, 'Condition number cond = ',cond                   )
             call info(routine_name, '================================================')
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
          call levels_postprocess_solution(common_krylov_data,lcommon_krylov_data)
          call time_end(t_postproc)
          if (myid.eq.0 .and. profile) then
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

    !-----profile
          if (profile) then
             call MPI_BARRIER(comm_all,ierr)
             call time_end(t_krylov_solve)
             if (myid.eq.0) then
                call time_print('solution by Krylov method',t_krylov_solve)
             end if
          end if
    !-----profile

          end subroutine

    !*******************************************************************************************
          subroutine krylov_bddcbicgstab(comm_all,tol,maxit,ndecrmax, num_iter,converged_reason)
    !*******************************************************************************************
    ! subroutine realizing BICGSTAB algorithm with vectors distributed by subdomains

    ! module for preconditioner
          use module_levels
    ! Program name
          use module_utils

          implicit none
          
          include "mpif.h"

          ! parallel variables
          integer,intent(in) :: comm_all 

          ! limit on iterations
          integer,intent(in) :: maxit

          ! limit on iterations with increasing residual
          integer,intent(in) :: ndecrmax

          ! desired accuracy of relative residual
          real(kr),intent(in) :: tol

          ! resulting number of iterations
          integer,intent(out) :: num_iter

          ! convergence reason
          !  =  0 - converged relative residual
          !  = -1 - reached limit on number of iterations
          !  = -2 - reached limit on number of iterations with nondecreasing residual
          integer,intent(out) :: converged_reason

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
          real(kr) :: t_postproc, t_krylov_solve

    !-----profile
          if (profile) then
             call MPI_BARRIER(comm_all,ierr)
             call time_start
          end if
    !-----profile

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

          ! compute norm of right-hand side
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
                call info(routine_name,'Norm of the right-hand side =',normrhs)
             end if
          end if

          ! Check of zero right-hand side => all zero solution
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
             bicgstab_data(isub_loc)%v(:) = 0._kr
          end do
          ! p = 0
          do isub_loc = 1,nsub_loc
             bicgstab_data(isub_loc)%p(:) = 0._kr
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
                call info (routine_name, 'iteration: ',dble(iter-0.5) )
                call info (routine_name, '          relative residual: ',relres)
             end if

    ! Check convergence in the half step
    !  relres < tol
             if (relres.lt.tol) then
                if (myid.eq.0) then
                   call info (routine_name, 'Number of BICGSTAB iterations: ',dble(iter-0.5) )
                end if
                num_iter = iter
                converged_reason = 0
                exit
             end if

             ! Check of decreasing of residual
             if (relres.lt.lastres) then
                ndecr = 0
             else
                ndecr = ndecr + 1
                if (ndecr.ge.ndecrmax) then
                   if (myid.eq.0) then
                      call warning(routine_name,'Residual did not decrease for maximal number of iterations:',ndecrmax)
                   end if
                   num_iter = iter
                   converged_reason = -2
                   exit
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
                call info (routine_name, 'iteration: ',dble(iter))
                call info (routine_name, '          relative residual: ',relres)
             end if

    ! Check convergence
    !  relres < tol
             if (relres.lt.tol) then
                if (myid.eq.0) then
                   call info (routine_name, 'Number of BICGSTAB iterations: ',dble(iter) )
                end if
                num_iter = iter
                converged_reason = 0
                exit
             end if

             ! Check number of iterations
             if (iter.eq.maxit) then
                if (myid.eq.0) then
                   call warning(routine_name,'Maximal number of iterations reached, precision not achieved.')
                end if
                num_iter = iter
                converged_reason = -1
                exit
             end if

             ! Check of decreasing of residual
             if (relres.lt.lastres) then
                ndecr = 0
             else
                ndecr = ndecr + 1
                if (ndecr.ge.ndecrmax) then
                   if (myid.eq.0) then
                      call warning(routine_name,'Residual did not decrease for maximal number of iterations:',ndecrmax)
                   end if
                   num_iter = iter
                   converged_reason = -2
                   exit
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
      call levels_postprocess_solution(common_krylov_data,lcommon_krylov_data)
      call time_end(t_postproc)
      if (myid.eq.0 .and. profile) then
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

!-----profile
      if (profile) then
         call MPI_BARRIER(comm_all,ierr)
         call time_end(t_krylov_solve)
         if (myid.eq.0) then
            call time_print('solution by Krylov method',t_krylov_solve)
         end if
      end if
!-----profile

      end subroutine

    !**************************************************************************************************
          subroutine krylov_bddcsteepestdescent(comm_all,tol,maxit,ndecrmax, num_iter,converged_reason)
    !**************************************************************************************************
    ! subroutine realizing steepest descent iteration algorithm with vectors distributed by subdomains

    ! module for preconditioner
          use module_levels
    ! Program name
          use module_utils

          implicit none
          
          include "mpif.h"

          ! parallel variables
          integer,intent(in) :: comm_all 

          ! limit on iterations
          integer,intent(in) :: maxit

          ! limit on iterations with increasing residual
          integer,intent(in) :: ndecrmax

          ! desired accuracy of relative residual
          real(kr),intent(in) :: tol

          ! resulting number of iterations
          integer,intent(out) :: num_iter

          ! convergence reason
          !  =  0 - converged relative residual
          !  = -1 - reached limit on number of iterations
          !  = -2 - reached limit on number of iterations with nondecreasing residual
          integer,intent(out) :: converged_reason

          ! local vars
          character(*),parameter:: routine_name = 'KRYLOV_BDDCSTEEPESTDESCENT'
          integer,parameter :: ilevel = 1

          ! data for storing actual data
          integer ::                                              lsteepestdescent_data
          type (steepestdescent_data_type), allocatable, target :: steepestdescent_data(:)

          ! data for auxiliary manipulation with preconditioner and system matrix 
          integer ::                                     lcommon_krylov_data
          type (common_krylov_data_type), allocatable ::  common_krylov_data(:)

          integer :: myid
          integer :: nsub, nsub_loc
          integer :: isub_loc
          integer :: iter, ndecr
          integer :: ndofis, nnodis

          ! steepest descent vars
          real(kr) :: normrhs, normres2, normres, normres2_loc, normres2_sub
          real(kr) :: rz, rz_loc, rz_sub
          real(kr) :: zaz, zaz_loc, zaz_sub
          real(kr) :: alpha 
          real(kr) :: relres, lastres

          ! MPI vars
          integer :: ierr

          ! time variables
          real(kr) :: t_postproc, t_solve

    !-----profile
          if (profile) then
             call MPI_BARRIER(comm_all,ierr)
             call time_start
          end if
    !-----profile

          ! orient in the communicator
          call MPI_COMM_RANK(comm_all,myid,ierr)

          ! find number of subdomains
          call levels_get_number_of_subdomains(ilevel,nsub,nsub_loc)

          ! prepare data and memory for steepest descent
          lcommon_krylov_data = nsub_loc
          allocate(common_krylov_data(lcommon_krylov_data))
          lsteepestdescent_data = nsub_loc
          allocate(steepestdescent_data(lsteepestdescent_data))
          do isub_loc = 1,nsub_loc
             call levels_dd_get_interface_size(ilevel,isub_loc, ndofis, nnodis)
             steepestdescent_data(isub_loc)%lsoli = ndofis
             allocate(steepestdescent_data(isub_loc)%soli(steepestdescent_data(isub_loc)%lsoli))
             steepestdescent_data(isub_loc)%lresi = ndofis
             allocate(steepestdescent_data(isub_loc)%resi(steepestdescent_data(isub_loc)%lresi))

             steepestdescent_data(isub_loc)%lz = ndofis
             allocate(steepestdescent_data(isub_loc)%z(steepestdescent_data(isub_loc)%lz))
             steepestdescent_data(isub_loc)%laz = ndofis
             allocate(steepestdescent_data(isub_loc)%az(steepestdescent_data(isub_loc)%laz))
          end do

          do isub_loc = 1,nsub_loc
             ! prepare initial solution and right-hand side
             call levels_prepare_interface_initial_data(isub_loc,&
                       steepestdescent_data(isub_loc)%soli,steepestdescent_data(isub_loc)%lsoli,&
                       steepestdescent_data(isub_loc)%resi,steepestdescent_data(isub_loc)%lresi)
             ! fix boundary conditions in residual to zero
             call levels_dd_fix_bc_interface_dual(ilevel,isub_loc,&
                       steepestdescent_data(isub_loc)%resi,steepestdescent_data(isub_loc)%lresi)
          end do

          ! get initial residual
          ! r_0 = g - A*u_0
          ! au = A*u_0
          ! first set pointers to soli and ap
          do isub_loc = 1,nsub_loc
             common_krylov_data(isub_loc)%lvec_in  = steepestdescent_data(isub_loc)%lsoli
             common_krylov_data(isub_loc)%vec_in  => steepestdescent_data(isub_loc)%soli
             common_krylov_data(isub_loc)%lvec_out = steepestdescent_data(isub_loc)%laz
             common_krylov_data(isub_loc)%vec_out => steepestdescent_data(isub_loc)%az
          end do
          call levels_sm_apply(common_krylov_data,lcommon_krylov_data)

          ! update residual
          ! r_0 = g - A*u_0
          do isub_loc = 1,nsub_loc
             steepestdescent_data(isub_loc)%resi = steepestdescent_data(isub_loc)%resi - steepestdescent_data(isub_loc)%az
          end do
          ! fix boundary conditions in residual to zero
          do isub_loc = 1,nsub_loc
             call levels_dd_fix_bc_interface_dual(ilevel,isub_loc,&
                     steepestdescent_data(isub_loc)%resi,steepestdescent_data(isub_loc)%lresi)
          end do

          ! compute norm of right-hand side
          normres2_loc = 0._kr
          do isub_loc = 1,nsub_loc
             call levels_dd_dotprod_local(ilevel,isub_loc,&
                                          steepestdescent_data(isub_loc)%resi,steepestdescent_data(isub_loc)%lresi, &
                                          steepestdescent_data(isub_loc)%resi,steepestdescent_data(isub_loc)%lresi, &
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
                call info(routine_name,'Norm of the right-hand side =',normrhs)
             end if
          end if

          ! Check of zero right-hand side => all zero solution
          if (normrhs.eq.0.0D0) then
             if (myid.eq.0) then
                call warning(routine_name,'initial residual zero => initial solution exact')
             end if
             return 
          end if

    ! Setting up the properties for decreasing residual
          ndecr   = 0
          lastres = 1.0_kr

     
    !***********************************************************************
    !*************************MAIN LOOP OVER ITERATIONS*********************
    !***********************************************************************
          do iter = 1,maxit

             ! Action of preconditioner M on vector resi 
             ! z = M*resi
             if (debug) then
                if (myid.eq.0) then
                   call info(routine_name,' Action of preconditioner')
                end if
             end if
             ! first set properly pointers
             do isub_loc = 1,nsub_loc
                common_krylov_data(isub_loc)%lvec_in  = steepestdescent_data(isub_loc)%lresi
                common_krylov_data(isub_loc)%vec_in  => steepestdescent_data(isub_loc)%resi
                common_krylov_data(isub_loc)%lvec_out = steepestdescent_data(isub_loc)%lz
                common_krylov_data(isub_loc)%vec_out => steepestdescent_data(isub_loc)%z
             end do
             call levels_pc_apply(common_krylov_data,lcommon_krylov_data)

             ! find A*z
             ! az = A*z
             ! first set pointers to soli and ap
             do isub_loc = 1,nsub_loc
                common_krylov_data(isub_loc)%lvec_in  = steepestdescent_data(isub_loc)%lz
                common_krylov_data(isub_loc)%vec_in  => steepestdescent_data(isub_loc)%z
                common_krylov_data(isub_loc)%lvec_out = steepestdescent_data(isub_loc)%laz
                common_krylov_data(isub_loc)%vec_out => steepestdescent_data(isub_loc)%az
             end do
             call levels_sm_apply(common_krylov_data,lcommon_krylov_data)

             ! determine ( r' * z ) 
             ! normsoli = ||soli||
             rz_loc = 0._kr
             do isub_loc = 1,nsub_loc
                call levels_dd_dotprod_local(ilevel,isub_loc, &
                                             steepestdescent_data(isub_loc)%resi,steepestdescent_data(isub_loc)%lresi, &
                                             steepestdescent_data(isub_loc)%z,steepestdescent_data(isub_loc)%lz, &
                                             rz_sub)
                rz_loc = rz_loc + rz_sub
             end do
    !***************************************************************PARALLEL
             call MPI_ALLREDUCE(rz_loc,rz, 1, MPI_DOUBLE_PRECISION, &
                                MPI_SUM, comm_all, ierr) 
    !***************************************************************PARALLEL

             ! determine ( z' * A * z ) 
             zaz_loc = 0._kr
             do isub_loc = 1,nsub_loc
                call levels_dd_dotprod_local(ilevel,isub_loc, &
                                             steepestdescent_data(isub_loc)%z,steepestdescent_data(isub_loc)%lz, &
                                             steepestdescent_data(isub_loc)%az,steepestdescent_data(isub_loc)%laz, &
                                             zaz_sub)
                zaz_loc = zaz_loc + zaz_sub
             end do
    !***************************************************************PARALLEL
             call MPI_ALLREDUCE(zaz_loc,zaz, 1, MPI_DOUBLE_PRECISION, &
                                MPI_SUM, comm_all, ierr) 
    !***************************************************************PARALLEL

             ! determine alpha
             ! alpha = ( r' * z ) / ( z' * A * z )
             alpha = rz / zaz

             ! update solution
             ! soli = soli + alpha * z 
             do isub_loc = 1,nsub_loc
                steepestdescent_data(isub_loc)%soli = steepestdescent_data(isub_loc)%soli &
                                                    + alpha * steepestdescent_data(isub_loc)%z 
             end do

             ! update residual
             ! resi = resi - alpha * Az 
             do isub_loc = 1,nsub_loc
                steepestdescent_data(isub_loc)%resi = steepestdescent_data(isub_loc)%resi &
                                                    - alpha * steepestdescent_data(isub_loc)%az 
                call levels_dd_fix_bc_interface_dual(ilevel,isub_loc,&
                                           steepestdescent_data(isub_loc)%resi,steepestdescent_data(isub_loc)%lresi)
             end do

             ! determine norm of residual 
             ! normres = ||resi||
             normres2_loc = 0._kr
             do isub_loc = 1,nsub_loc
                call levels_dd_dotprod_local(ilevel,isub_loc, &
                                             steepestdescent_data(isub_loc)%resi,steepestdescent_data(isub_loc)%lresi, &
                                             steepestdescent_data(isub_loc)%resi,steepestdescent_data(isub_loc)%lresi, &
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
                call info (routine_name, 'iteration: ',iter)
                call info (routine_name, '          relative residual: ',relres)
             end if

    ! Check convergence in the half step
    !  relres < tol
             if (relres.lt.tol) then
                if (myid.eq.0) then
                   call info (routine_name, ': Number of steepest descent iterations: ',iter )
                end if
                num_iter = iter
                converged_reason = 0
                exit
             end if

             ! Check number of iterations
             if (iter.eq.maxit) then
                if (myid.eq.0) then
                   call warning(routine_name,'Maximal number of iterations reached, precision not achieved.')
                end if
                num_iter = iter
                converged_reason = -1
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

          end do
    !*************************END OF MAIN LOOP OVER ITERATIONS**************

          ! Postprocessing of solution - computing interior values
          call time_start
          ! first set pointers to soli
          do isub_loc = 1,nsub_loc
             common_krylov_data(isub_loc)%lvec_in  = steepestdescent_data(isub_loc)%lsoli
             common_krylov_data(isub_loc)%vec_in  => steepestdescent_data(isub_loc)%soli
          end do
      call levels_postprocess_solution(common_krylov_data,lcommon_krylov_data)
      call time_end(t_postproc)
      if (myid.eq.0 .and. profile) then
         call time_print('postprocessing of solution',t_postproc)
      end if

      ! Clear memory of BICGSTAB
      do isub_loc = 1,nsub_loc
         nullify(common_krylov_data(isub_loc)%vec_in)
         nullify(common_krylov_data(isub_loc)%vec_out)
      end do
      deallocate(common_krylov_data)
      do isub_loc = 1,nsub_loc
         deallocate(steepestdescent_data(isub_loc)%soli)
         deallocate(steepestdescent_data(isub_loc)%resi)

         deallocate(steepestdescent_data(isub_loc)%z)
         deallocate(steepestdescent_data(isub_loc)%az)
      end do
      deallocate(steepestdescent_data)

!-----profile
      if (profile) then
         call MPI_BARRIER(comm_all,ierr)
         call time_end(t_solve)
         if (myid.eq.0) then
            call time_print('solution by steepest descent method',t_solve)
         end if
      end if
!-----profile

      end subroutine

  !************************************************************************
      subroutine krylov_orthogonalize_cgs(comm_all,krylov_data,lkrylov_data)
  !************************************************************************
      ! project vector P onto the stored Krylov space A*p
      ! using classical Gram-Schmidt orthogonalization
      use module_levels
      use module_utils

      implicit none
      
      include "mpif.h"

      ! parallel variables
      integer,intent(in) :: comm_all 

      ! data for PCG
      integer,intent(in) ::                         lkrylov_data
      type(common_krylov_data_type),intent(inout) :: krylov_data(lkrylov_data) 

      ! local vars
      character(*),parameter:: routine_name = 'KRYLOV_ORTHOGONALIZE_CGS'
      integer,parameter :: ilevel = 1

      integer ::              lwtp
      real(kr), allocatable :: wtp(:)
      real(kr), allocatable :: wtp_loc(:)
      real(kr) :: wtp_sub

      integer :: isub_loc, ibasis, nsub, nsub_loc 
      integer :: ierr

      ! find number of subdomains
      call levels_get_number_of_subdomains(ilevel,nsub,nsub_loc)

      lwtp     = nactive_cols_recycling_basis 
      allocate(wtp_loc(lwtp))
      wtp_loc = 0._kr
      allocate(wtp(lwtp))

      ! W'*p
      do ibasis = 1,nactive_cols_recycling_basis
         do isub_loc = 1,nsub_loc
            call levels_dd_dotprod_local(ilevel,isub_loc,&
                                         recycling_basis(isub_loc)%w(:,ibasis),recycling_basis(isub_loc)%lw1, &
                                         krylov_data(isub_loc)%vec_in,krylov_data(isub_loc)%lvec_in, &
                                         wtp_sub)
            wtp_loc(ibasis) = wtp_loc(ibasis) + wtp_sub
         end do
      end do
!***************************************************************PARALLEL
      call MPI_ALLREDUCE(wtp_loc,wtp, lwtp, MPI_DOUBLE_PRECISION,&
                         MPI_SUM, comm_all, ierr) 
!***************************************************************PARALLEL
      deallocate(wtp_loc)

      ! inv(D) * W' * p
      call recycling_apply_diagonal_inverse(wtp,lwtp)
      !wtp = wtp * recycling_idvtw(1:nactive_cols_recycling_basis)

      ! correct search direction vector
      ! p <- p - V*inv(D)*W'p = p - V*IDIAG*W'p
      do isub_loc = 1,nsub_loc
         krylov_data(isub_loc)%vec_in = krylov_data(isub_loc)%vec_in &
                                      - matmul(recycling_basis(isub_loc)%v(:,1:nactive_cols_recycling_basis), &
                                               wtp(1:nactive_cols_recycling_basis) )
      end do

      deallocate(wtp)

      end subroutine

  !**************************************************************************
      subroutine krylov_orthogonalize_icgs(comm_all,krylov_data,lkrylov_data)
  !**************************************************************************
      ! project vector P onto the stored Krylov space A*p
      ! using iterated classical Gram-Schmidt orthogonalization
      ! based on Algorithm A3 from 
      ! Tebbens, Hnetynkova, Plesinger, Strakos, Tichy
      ! Analyza metod pro maticove vypocty
      ! Matfyzpress, Praha, 2012

      use module_levels
      use module_utils

      implicit none
      
      include "mpif.h"

      ! parallel variables
      integer,intent(in) :: comm_all 

      ! data for PCG
      integer,intent(in) ::                         lkrylov_data
      type(common_krylov_data_type),intent(inout) :: krylov_data(lkrylov_data) 

      ! local vars
      integer :: i

      ! perform two iterations of classical Gram-Schmidt
      do i = 1,2
         call krylov_orthogonalize_cgs( comm_all, krylov_data, lkrylov_data)
      end do

      end subroutine

  !*************************************************************************
      subroutine krylov_orthogonalize_mgs(comm_all,krylov_data,lkrylov_data)
  !*************************************************************************
      ! project vector P onto the stored Krylov space A*p
      ! using modified Gram-Schmidt orthogonalization
      ! this algorithm needs a lot of synchronization by global MPI functions
      use module_levels
      use module_utils

      implicit none
      
      include "mpif.h"

      ! parallel variables
      integer,intent(in) :: comm_all 

      ! data for PCG
      integer,intent(in) ::                         lkrylov_data
      type(common_krylov_data_type),intent(inout) :: krylov_data(lkrylov_data) 

      ! local vars
      character(*),parameter:: routine_name = 'KRYLOV_ORTHOGONALIZE_MGS'
      integer,parameter :: ilevel = 1

      real(kr) :: wtp
      real(kr) :: wtp_loc
      real(kr) :: wtp_sub

      integer :: isub_loc, ibasis, nsub, nsub_loc 
      integer :: ierr

      ! find number of subdomains
      call levels_get_number_of_subdomains(ilevel,nsub,nsub_loc)

      ! W'*p
      do ibasis = 1,nactive_cols_recycling_basis
         wtp_loc = 0._kr
         do isub_loc = 1,nsub_loc
            call levels_dd_dotprod_local(ilevel,isub_loc,&
                                         recycling_basis(isub_loc)%w(:,ibasis),recycling_basis(isub_loc)%lw1, &
                                         krylov_data(isub_loc)%vec_in,krylov_data(isub_loc)%lvec_in, &
                                         wtp_sub)
            wtp_loc = wtp_loc + wtp_sub
         end do
!***************************************************************PARALLEL
         call MPI_ALLREDUCE(wtp_loc,wtp, 1, MPI_DOUBLE_PRECISION,&
                            MPI_SUM, comm_all, ierr) 
!***************************************************************PARALLEL

         ! inv(D) * W' * p
         wtp = wtp * recycling_idvtw(ibasis)

         ! correct search direction vector
         ! p <- p - V*inv(D)+W'p = p - V*IDIAG*W'p
         do isub_loc = 1,nsub_loc
            krylov_data(isub_loc)%vec_in = krylov_data(isub_loc)%vec_in &
                                         - wtp * recycling_basis(isub_loc)%v(:,ibasis) 
         end do
      end do

      end subroutine

      !****************************************************
      subroutine recycling_apply_diagonal_inverse(vec,lvec)
      !****************************************************
      ! application of the inverse of the diagonal matrix in recycling
      use module_utils

      implicit none

      integer,intent(in)     :: lvec
      real(kr),intent(inout) ::  vec(lvec) 

      ! local vars
      character(*),parameter:: routine_name = 'RECYCLING_APPLY_DIAGONAL_INVERSE'
      !integer :: nrhs 
      !integer :: lapack_info 

      if (lvec /= nactive_cols_recycling_basis) then
         call error( routine_name, 'size mismatch', lvec )
      end if
      ! only diagonal
      vec(1:lvec) = vec(1:lvec) * recycling_idvtw(1:nactive_cols_recycling_basis)

      ! actually apply the inverse of the Grammian
      !if ( .not. recycling_is_inverse_prepared ) then
      !   call error( routine_name, 'Gramian of the basis is not factorized', lvec )
      !end if

      !nrhs = 1
      !call DGETRS('N', lvec, nrhs, recycling_vtw, lvec, recycling_ipiv, vec, lvec, lapack_info)
      !if (lapack_info /= 0) then
      !   call error( routine_name, 'error in LAPACK, info ', lapack_info )
      !end if

      end subroutine

      !***********************************************
      subroutine krylov_condsparse(nw,d,ld,e,le, cond)
      !***********************************************
      ! Routine that estimates condition number of a real symmetric tridiagonal matrix 
      ! using LAPACK
            use module_utils
            implicit none
      
      ! used dimension of matrix
            integer,intent(in) :: nw
      
      ! diagonal of matrix
            integer,intent(in) :: ld
            real(kr),intent(in) :: d(ld)
      ! subdiagonal of matrix
            integer,intent(in) :: le
            real(kr),intent(in) :: e(le)
      
      ! estimate of the condition number
            real(kr),intent(out) :: cond
      
            ! auxiliary variables
            character(*),parameter:: routine_name = 'KRYLOV_CONDSPARSE'
            integer ::  iaux
            real(kr) :: raux(1)
      
            real(kr) :: eigmax
      
            ! LAPACK
            integer :: lapack_info
      
            ! checks
            if (ld.ne.nw .or. le .ne. ld-1) then
               call error(routine_name,'Dimensions mismatch.')
            end if
      
            ! return silently if only 1 or less iterations were performed
            if (nw.le.1) then
               cond = 1._kr
               return
            end if
      
            ! LAPACK routine for searching eigenvalues (returned in ascending order)
            iaux = 1
            call DSTEV( 'N', nw, d, e, raux, iaux, raux, lapack_info)
            if (debug) then
               write(*,*) 'eigenvalues = '
               write(*,*) d(1:nw)
            end if
      
            ! compute condition number
            eigmax = d(nw)
            ! do not get the lowest eigenvalue from the Lanczos sequence - the Ritz value may not converge (cf Treffethen, Bau)
            if (debug) then
               write(*,*) 'eigmax = ',eigmax
            end if
            cond = eigmax
      
      end subroutine



      !*******************************
      subroutine krylov_set_profile_on
      !*******************************
      ! auxiliary routine to switch profiling flag on
      profile = .true.
      end subroutine
      
      !********************************
      subroutine krylov_set_profile_off
      !********************************
      ! auxiliary routine to switch profiling flag off
      profile = .false.
      end subroutine

      !*************************
      subroutine krylov_finalize
      !*************************
      ! deallocating data from recycling
      !local vars
      integer :: i

      do i = 1, lrecycling_basis
         if (allocated(recycling_basis(i)%v)) then
            deallocate(recycling_basis(i)%v)
            recycling_basis(i)%lv1 = 0
            recycling_basis(i)%lv2 = 0
         end if
         if (allocated(recycling_basis(i)%w)) then
            deallocate(recycling_basis(i)%w)
            recycling_basis(i)%lw1 = 0
            recycling_basis(i)%lw2 = 0
         end if
      end do
      if (allocated(recycling_basis)) then
         deallocate(recycling_basis)
         lrecycling_basis = 0
         nactive_cols_recycling_basis = 0
      end if
      if (allocated(recycling_idvtw)) then
         deallocate(recycling_idvtw)
         lrecycling_idvtw = 0
      end if
      is_recycling_prepared = .false.

      end subroutine

end module
