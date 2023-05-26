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


! Wrapper module for generating C functions. 
! Jakub Sistek, Praha 12/2010

module module_bddcml_c_wrapper
!*****************************
! Module for BDDCML interface
! Jakub Sistek, Praha 2015

use iso_c_binding, only: c_double, c_int
use module_bddcml

implicit none

integer,parameter,private :: c_real_type    = c_double
integer,parameter,private :: c_integer_type = c_int

contains

!*********************************************************************************************************************
subroutine bddcml_init_c(nl, nsublev,lnsublev, nsub_loc_1, comm_init, verbose_level, numbase, just_direct_solve_int) &
           bind(c)
!*********************************************************************************************************************
! initialization of LEVELS module

! given number of levels
      integer(c_integer_type),intent(in) :: nl

! GLOBAL numbers of subdomains in all levels
      integer(c_integer_type),intent(in) :: lnsublev
      integer(c_integer_type),intent(in) ::  nsublev(lnsublev)

! LOCAL number of subdomains assigned to the process
!     >= 0 - number of local subdomains - sum up across processes to nsublev(1)
!     -1   - let solver decide, then the value is returned ( determining linear partition )
      integer(c_integer_type), intent(inout):: nsub_loc_1 

! initial global communicator (possibly MPI_COMM_WORLD)
      integer(c_integer_type), intent(in):: comm_init 

! level of verbosity
!     0 - only errors printed
!     1 - some output
!     2 - detailed output
      integer(c_integer_type), intent(in):: verbose_level 

! first index of arrays ( 0 for C, 1 for Fortran )
      integer(c_integer_type), intent(in) :: numbase

! if >0, only perform parallel solve by a direct solver, e.g. parallel MUMPS
! if ==0 perform regular solve by BDDC
      integer(c_integer_type), intent(in)::     just_direct_solve_int 

      call bddcml_init(nl, nsublev,lnsublev, nsub_loc_1, comm_init, verbose_level, numbase, just_direct_solve_int)

end subroutine

!**************************************************************************************
subroutine bddcml_upload_global_data_c(nelem,nnod,ndof,ndim,meshdim,&
                                       inet,linet,nnet,lnnet,nndf,lnndf,xyz,lxyz1,lxyz2,&
                                       ifix,lifix, fixv,lfixv, rhs,lrhs, sol,lsol, idelm, &
                                       neighbouring, load_division_int) &
           bind(c)
!**************************************************************************************
! Subroutine for loading global data as zero level
      ! GLOBAL number of elements
      integer(c_integer_type), intent(in):: nelem
      ! GLOBAL number of nodes
      integer(c_integer_type), intent(in):: nnod
      ! GLOBAL number of degrees of freedom
      integer(c_integer_type), intent(in):: ndof

      ! space dimensions
      integer(c_integer_type),intent(in) :: ndim 

      ! mesh dimension
      ! for 3D problems, = ndim
      ! for 3D shells = 2
      ! for 3D beams  = 1
      integer(c_integer_type),intent(in) :: meshdim


      ! GLOBAL Indices of Nodes on ElemenTs
      integer(c_integer_type), intent(in):: linet
      integer(c_integer_type), intent(in)::  inet(linet)

      ! GLOBAL Number of Nodes on ElemenTs
      integer(c_integer_type), intent(in):: lnnet
      integer(c_integer_type), intent(in)::  nnet(lnnet)

      ! GLOBAL Number of Nodal Degrees of Freedom
      integer(c_integer_type), intent(in):: lnndf
      integer(c_integer_type), intent(in)::  nndf(lnndf)

      ! GLOBAL Coordinates of nodes (X | Y | Z) or as one array (all X, all Y, all Z)
      integer(c_integer_type), intent(in):: lxyz1, lxyz2
      real(c_real_type), intent(in):: xyz(lxyz1,lxyz2)

      ! GLOBAL Indices of FIXed variables - all dof with Dirichlet BC are marked with its number
      integer(c_integer_type), intent(in):: lifix
      integer(c_integer_type), intent(in)::  ifix(lifix)

      ! GLOBAL FIXed Variables - where IFIX is nonzero, value of Dirichlet BC
      integer(c_integer_type), intent(in):: lfixv
      real(c_real_type), intent(in):: fixv(lfixv)

      ! GLOBAL Right-Hand Side
      integer(c_integer_type), intent(in):: lrhs
      real(c_real_type), intent(in):: rhs(lrhs)

      ! GLOBAL initial SOLution - initial approximation for iterative method
      integer(c_integer_type), intent(in):: lsol
      real(c_real_type), intent(in):: sol(lsol)

      ! opened unit with unformatted file with element matrices
      integer(c_integer_type),intent(in) :: idelm  

      ! how many nodes should be shared for calling element adjacent in graph?
      integer(c_integer_type),intent(in) :: neighbouring 

      ! use division from 'partition_l1.ES' file?
      integer(c_integer_type),intent(in) :: load_division_int

      call bddcml_upload_global_data(nelem,nnod,ndof,ndim,meshdim,&
                                     inet,linet,nnet,lnnet,nndf,lnndf,xyz,lxyz1,lxyz2,&
                                     ifix,lifix, fixv,lfixv, rhs,lrhs, sol,lsol, idelm, &
                                     neighbouring, load_division_int)

end subroutine

!**************************************************************************************
subroutine bddcml_upload_subdomain_data_c(nelem, nnod, ndof, ndim, meshdim, &
                                          isub, nelems, nnods, ndofs, &
                                          inet,linet, nnet,lnnet, nndf,lnndf, &
                                          isngn,lisngn, isvgvn,lisvgvn, isegn,lisegn, &
                                          xyz,lxyz1,lxyz2, &
                                          ifix,lifix, fixv,lfixv, &
                                          rhs,lrhs, is_rhs_complete_int, &
                                          sol,lsol, &
                                          matrixtype, i_sparse, j_sparse, a_sparse, la, is_assembled_int,&
                                          user_constraints,luser_constraints1,luser_constraints2,&
                                          element_data,lelement_data1,lelement_data2, &
                                          dof_data,ldof_data, &
                                          find_components_int, use_dual_mesh_graph_int, neighbouring) &
           bind(c)
!**************************************************************************************
! Subroutine for loading global data as zero level

      ! GLOBAL number of elements
      integer(c_integer_type), intent(in):: nelem
      ! GLOBAL number of nodes
      integer(c_integer_type), intent(in):: nnod
      ! GLOBAL number of degrees of freedom
      integer(c_integer_type), intent(in):: ndof
      ! GLOBAL number of spacial dimensions
      integer(c_integer_type), intent(in):: ndim
      ! mesh dimension
      ! for 3D problems, = ndim
      ! for 3D shells = 2
      ! for 3D beams  = 1
      integer(c_integer_type),intent(in) :: meshdim

      ! GLOBAL index of subdomain
      integer(c_integer_type), intent(in):: isub
      ! LOCAL number of elements
      integer(c_integer_type), intent(in):: nelems
      ! LOCAL number of nodes
      integer(c_integer_type), intent(in):: nnods
      ! LOCAL number of degrees of freedom
      integer(c_integer_type), intent(in):: ndofs

      ! LOCAL Indices of Nodes on ElemenTs
      integer(c_integer_type), intent(in):: linet
      integer(c_integer_type), intent(in)::  inet(linet)

      ! LOCAL Number of Nodes on ElemenTs
      integer(c_integer_type), intent(in):: lnnet
      integer(c_integer_type), intent(in)::  nnet(lnnet)

      ! LOCAL Number of Nodal Degrees of Freedom
      integer(c_integer_type), intent(in):: lnndf
      integer(c_integer_type), intent(in)::  nndf(lnndf)

      ! Indices of Subdomain Nodes in Global Numbering (local to global map of nodes)
      integer(c_integer_type), intent(in):: lisngn
      integer(c_integer_type), intent(in)::  isngn(lisngn)

      ! Indices of Subdomain Variables in Global Variable Numbering (local to global map of variables)
      integer(c_integer_type), intent(in):: lisvgvn
      integer(c_integer_type), intent(in)::  isvgvn(lisvgvn)

      ! Indices of Subdomain Elements in Global Numbering (local to global map of elements)
      integer(c_integer_type), intent(in):: lisegn
      integer(c_integer_type), intent(in)::  isegn(lisegn)

      ! LOCAL Coordinates of nodes (X | Y | Z (if present) ) or as one array (all X, all Y, all Z) 
      ! lxyz1 = NNODS, lxyz2 = NDIM
      integer(c_integer_type), intent(in):: lxyz1, lxyz2
      real(c_real_type), intent(in):: xyz(lxyz1,lxyz2)

      ! LOCAL Indices of FIXed variables - all dof with Dirichlet BC are marked with its local number
      integer(c_integer_type), intent(in):: lifix
      integer(c_integer_type), intent(in)::  ifix(lifix)

      ! LOCAL FIXed Variables - where IFIX is nonzero, value of Dirichlet BC
      integer(c_integer_type), intent(in):: lfixv
      real(c_real_type), intent(in):: fixv(lfixv)

      ! LOCAL Right-Hand Side (where nodes of subdomains coincide, values are repeated)
      integer(c_integer_type), intent(in):: lrhs
      real(c_real_type), intent(in):: rhs(lrhs)

      ! LOCAL initial SOLution - initial approximation for iterative method
      integer(c_integer_type), intent(in):: lsol
      real(c_real_type), intent(in):: sol(lsol)
      integer(c_integer_type), intent(in)::  is_rhs_complete_int  ! is the right-hand side complete?  
                                                  !  0 = no, e.g. if it is assembled subdomain-by-subdomain 
                                                  !      and the interface entries are not interchanged
                                                  !  1 = yes, e.g. if it was created as a restriction of
                                                  !      the global RHS array and so entries on interface are
                                                  !      contained in several subdomains - weights will be applied

      ! LOCAL matrix triplet i, j, a(i,j) 
      ! Type of the matrix
      ! 0 - unsymmetric
      ! 1 - symmetric positive definite
      ! 2 - general symmetric
      integer(c_integer_type), intent(in)::  matrixtype 
      integer(c_integer_type), intent(in)::  la            ! length of previous arrays (= number of nonzeros for assembled matrix)
      integer(c_integer_type), intent(in)::  i_sparse(la)  ! array of row indices
      integer(c_integer_type), intent(in)::  j_sparse(la)  ! array of column indices
      real(c_real_type), intent(in):: a_sparse(la)  ! array of values
      integer(c_integer_type), intent(in)::  is_assembled_int  ! is the array assembled? 
                                               !  0 = no, it can contain repeated entries
                                               !  1 = yes, it is sorted and doesn't contain repeated index pairs
      ! LOCAL array of constraints defined by user
      integer(c_integer_type), intent(in)::  luser_constraints1 ! number of rows in matrix of constraints, i.e. number of user constraints
      integer(c_integer_type), intent(in)::  luser_constraints2 ! number of columns in matrix of constraints, ( = NNODS)
      real(c_real_type), intent(in):: user_constraints(luser_constraints1*luser_constraints2) ! array for additional constraints

      ! LOCAL array of additional element data - can be used e.g. for construction of averaging operator E
      integer(c_integer_type), intent(in)::  lelement_data1 ! number of rows in matrix of element data, i.e. number of data per element
      integer(c_integer_type), intent(in)::  lelement_data2 ! number of columns in matrix of constraints, ( = NELEMS)
      real(c_real_type), intent(in):: element_data(lelement_data1*lelement_data2) ! array for data on elements

      ! LOCAL array of additional data at degrees of freedom 
      integer(c_integer_type), intent(in)::  ldof_data  ! number of entries in dof_data ( = NDOFS)
      real(c_real_type), intent(in):: dof_data(ldof_data) ! array for data on degrees of freedom

      ! should the mesh components be detected ? 
      ! Should be the same for all subdomains.
      integer(c_integer_type), intent(in)::  find_components_int  ! 0 = mesh components will not be detected and the subdomain is considered to be connected
                                                                  ! 1 = components will be detected based on the dual graph of the mesh 
                                                                  ! This is recommended for unstructured meshes where a subdomains might be of disconnected parts. 
                                                                  ! However it may become expensive for high order elements.

      ! if find_components_int = 1, should the dual graph of mesh be used for detecting components? 
      ! Not accessed if detection of subdomain components is switched off by find_components_int = 0.
      ! Should be the same for all subdomains.
      integer(c_integer_type), intent(in)::  use_dual_mesh_graph_int  ! 0 = mesh components will be detected from primal graph of the mesh, 
                                                                      !     where element nodes are graph vertices and a graph edge is introduced 
                                                                      !     if they are connected by the same element
                                                                      ! 1 = dual graph of mesh will be used for detecting components.
                                                                      !     Dual graph of mesh contains elements as graph vertices and a graph edge is
                                                                      !     introduced if two elements share at least the number of nodes prescribed by
                                                                      !     the neighbouring parameter.

      ! How many nodes need to be shared by two elements to declare a graph edge between them? 
      ! Accessed only if find_components_int = 1 and use_dual_mesh_graph_int = 1. A graph edge is introduced between two elements if
      ! they share number of nodes >= neighbouring. Should be the same for all subdomains.
      integer(c_integer_type), intent(in)::  neighbouring         


      call bddcml_upload_subdomain_data(nelem, nnod, ndof, ndim, meshdim, &
                                        isub, nelems, nnods, ndofs, &
                                        inet,linet, nnet,lnnet, nndf,lnndf, &
                                        isngn,lisngn, isvgvn,lisvgvn, isegn,lisegn, &
                                        xyz,lxyz1,lxyz2, &
                                        ifix,lifix, fixv,lfixv, &
                                        rhs,lrhs, is_rhs_complete_int, &
                                        sol,lsol, &
                                        matrixtype, i_sparse, j_sparse, a_sparse, la, is_assembled_int,&
                                        user_constraints,luser_constraints1,luser_constraints2,&
                                        element_data,lelement_data1,lelement_data2, &
                                        dof_data,ldof_data, &
                                        find_components_int, use_dual_mesh_graph_int, neighbouring) 
end subroutine

!************************************************************************************************
subroutine bddcml_setup_preconditioner_c(matrixtype, use_defaults_int, &
                                         parallel_division_int, &
                                         use_corner_constraints_int, &
                                         use_arithmetic_constraints_int, &
                                         use_adaptive_constraints_int, &
                                         use_user_constraints_int, &
                                         weights_type) &
           bind(c)
!************************************************************************************************
      ! numerical properties of the matrix (MUMPS-like notation)
      !  0 - general
      !  1 - symmetric positive definite
      !  2 - symmetric general
      integer(c_integer_type),intent(in) :: matrixtype 

      ! use default values for all other settings?
      integer(c_integer_type),intent(in) :: use_defaults_int

      ! should parallel division be used (ParMETIS instead of METIS)?
      integer(c_integer_type),intent(in) :: parallel_division_int

      ! use continuity at corners as constraints?
      integer(c_integer_type),intent(in) :: use_corner_constraints_int

      ! use arithmetic constraints?
      integer(c_integer_type),intent(in) :: use_arithmetic_constraints_int

      ! use adaptive constraints?
      integer(c_integer_type),intent(in) :: use_adaptive_constraints_int

      ! use user constraints?
      integer(c_integer_type),intent(in) :: use_user_constraints_int

      ! what type of weights should be used on the interface ?
      ! 0 - weights by cardinality
      ! 1 - weights by diagonal stiffness
      ! 2 - weights based on first row of element data
      ! 3 - weights based on dof data
      ! 4 - weights by Marta Certikova - unit load
      ! 5 - weights by Marta Certikova - unit jump
      ! 6 - weights by Schur row sums for whole subdomain
      ! 7 - weights by Schur row sums computed face by face
      integer(c_integer_type),intent(in) :: weights_type

      call bddcml_setup_preconditioner(matrixtype, use_defaults_int, &
                                       parallel_division_int, &
                                       use_corner_constraints_int, &
                                       use_arithmetic_constraints_int, &
                                       use_adaptive_constraints_int, &
                                       use_user_constraints_int, &
                                       weights_type) 
end subroutine


!******************************************************************
subroutine bddcml_solve_c(comm_all,method,tol,maxit,ndecrmax, &
                          recycling_int, max_number_of_stored_vectors, &
                          num_iter,converged_reason,condition_number) &
           bind(c)
!******************************************************************
      ! parallel variables
      integer(c_integer_type),intent(in) :: comm_all 

      ! optional paramaters - use either none or all of them
      ! preferable Krylov method
      ! -1 - use defaults (tol, maxit, and ndecrmax not accessed)
      ! 0 - PCG
      ! 1 - BICGSTAB
      ! 2 - steepest descent method
      ! 5 - direct solve by MUMPS
      integer(c_integer_type), intent(in) :: method

      ! desired accuracy of relative residual
      real(c_real_type), intent(in) :: tol

      ! limit on iterations
      integer(c_integer_type), intent(in) :: maxit

      ! limit on iterations with increasing residual
      integer(c_integer_type), intent(in) :: ndecrmax

      ! Should the Krylov subspace be recycled? If yes, basis will be stored and
      ! used for subsequent right hand sides. It also launches re-orthogonalization 
      ! of the Krylov basis.
      integer(c_integer_type), intent(in) :: recycling_int

      ! If the Krylov subspace is recycled, what is the upper bound for the number of stored vectors?
      ! With increasing number, the memory consumption also increases.
      integer(c_integer_type), intent(in) :: max_number_of_stored_vectors

      ! resulting number of iterations
      integer(c_integer_type),intent(out) :: num_iter

      ! convergence reason
      !  =  0 - converged relative residual
      !  = -1 - reached limit on number of iterations
      !  = -2 - reached limit on number of iterations with nondecreasing residual
      integer(c_integer_type),intent(out) :: converged_reason

      ! estimated condition number ( for PCG only )
      real(c_real_type),intent(out) :: condition_number

      call bddcml_solve(comm_all,method,tol,maxit,ndecrmax, &
                        recycling_int, max_number_of_stored_vectors, &
                        num_iter,converged_reason,condition_number) 
end subroutine

!**************************************************************
subroutine bddcml_download_local_solution_c(isub, sols,lsols) &
           bind(c)
!**************************************************************
! Subroutine for getting local solution,
! i.e. restriction of solution vector to subdomain (no weights are applied)

      ! GLOBAL index of subdomain
      integer(c_integer_type), intent(in)::  isub

      ! LOCAL solution 
      integer(c_integer_type), intent(in)::  lsols
      real(c_real_type), intent(out):: sols(lsols)

      call bddcml_download_local_solution(isub, sols,lsols)
end subroutine

!********************************************************
subroutine bddcml_download_global_solution_c(sol, lsol) &
           bind(c)
!********************************************************
! Subroutine for getting global solution at root process
      ! GLOBAL solution - required to be allocated only at root
      integer(c_integer_type), intent(in)::  lsol
      real(c_real_type), intent(out):: sol(lsol)

      call bddcml_download_global_solution(sol, lsol) 
end subroutine

!***************************************************************
subroutine bddcml_download_local_reactions_c(isub, reas,lreas) &
           bind(c)
!***************************************************************
! Subroutine for getting local reactions at Dirichlet BC,
! i.e. restriction of vector of reactions to subdomain (no weights are applied)

      ! GLOBAL index of subdomain
      integer(c_integer_type), intent(in)::  isub

      ! LOCAL reactions 
      integer(c_integer_type), intent(in)::  lreas
      real(c_real_type), intent(out):: reas(lreas)

      call bddcml_download_local_reactions(isub, reas,lreas)
end subroutine

!*********************************************************
subroutine bddcml_download_global_reactions_c(rea, lrea) &
           bind(c)
!*********************************************************
! Subroutine for getting global vector of reactions at Dirichlet BC at root process

      ! GLOBAL  - required to be allocated only at root
      integer(c_integer_type), intent(in)::  lrea
      real(c_real_type), intent(out):: rea(lrea)

      call bddcml_download_global_reactions(rea, lrea)
end subroutine

!***********************************************************************************
subroutine bddcml_change_global_data_c(ifix,lifix, fixv,lfixv, rhs,lrhs, sol,lsol) &
           bind(c)
!***********************************************************************************
! Subroutine for changing global RHS, and values of Dirichlet boundary conditions
      ! GLOBAL Indices of FIXed variables - all dof with Dirichlet BC are marked with its number
      integer(c_integer_type), intent(in):: lifix
      integer(c_integer_type), intent(in)::  ifix(lifix)

      ! GLOBAL FIXed Variables - where IFIX is nonzero, value of Dirichlet BC
      integer(c_integer_type), intent(in):: lfixv
      real(c_real_type), intent(in):: fixv(lfixv)

      ! GLOBAL Right-Hand Side
      integer(c_integer_type), intent(in):: lrhs
      real(c_real_type), intent(in):: rhs(lrhs)

      ! GLOBAL initial SOLution - initial approximation for iterative method
      integer(c_integer_type), intent(in):: lsol
      real(c_real_type), intent(in):: sol(lsol)

      call bddcml_change_global_data(ifix,lifix, fixv,lfixv, rhs,lrhs, sol,lsol)
end subroutine

!***********************************************************************
subroutine bddcml_change_subdomain_data_c(isub, &
                                          ifix,lifix, fixv,lfixv, &
                                          rhs,lrhs, is_rhs_complete_int, &
                                          sol,lsol) &
           bind(c)
!***********************************************************************
! Subroutine for loading global data as zero level

      ! GLOBAL index of subdomain
      integer(c_integer_type), intent(in):: isub

      ! LOCAL Indices of FIXed variables - all dof with Dirichlet BC are marked with its local number
      ! cannot change - just for reference
      integer(c_integer_type), intent(in):: lifix
      integer(c_integer_type), intent(in)::  ifix(lifix)

      ! LOCAL FIXed Variables - where IFIX is nonzero, value of Dirichlet BC
      integer(c_integer_type), intent(in):: lfixv
      real(c_real_type), intent(in):: fixv(lfixv)

      ! LOCAL Right-Hand Side (where nodes of subdomains coincide, values are repeated)
      integer(c_integer_type), intent(in):: lrhs
      real(c_real_type), intent(in):: rhs(lrhs)

      ! LOCAL initial SOLution - initial approximation for iterative method
      integer(c_integer_type), intent(in):: lsol
      real(c_real_type), intent(in):: sol(lsol)
      integer(c_integer_type), intent(in)::  is_rhs_complete_int  ! is the right-hand side complete?  
                                                  !  0 = no, e.g. if it is assembled subdomain-by-subdomain 
                                                  !      and the interface entries are not interchanged
                                                  !  1 = yes, e.g. if it was created as a restriction of
                                                  !      the global RHS array and so entries on interface are
                                                  !      contained in several subdomains - weights will be applied
      call bddcml_change_subdomain_data(isub, &
                                        ifix,lifix, fixv,lfixv, &
                                        rhs,lrhs, is_rhs_complete_int, &
                                        sol,lsol)
end subroutine

!*************************************
subroutine bddcml_setup_new_data_c() &
           bind(c)
!*************************************
      call bddcml_setup_new_data 
end subroutine

!*******************************************************************************
subroutine bddcml_dotprod_subdomain_c( isub, vec1,lvec1, vec2,lvec2, dotprod ) &
           bind(c)
!*******************************************************************************
! Auxiliary subroutine to compute scalar product of two vectors of lenght of
! subdomain exploiting interface weights from the solver. This routine is useful 
! if we want to compute global norm or dot product based on vectors restricted to 
! subdomains. Since interface values are contained in several vectors for
! several subdomains, this dot product or norm cannot be determined without
! weights.

      ! GLOBAL index of subdomain
      integer(c_integer_type),intent(in) ::   isub 
      ! vectors to multiply
      integer(c_integer_type),intent(in) ::  lvec1        ! length of the first vector
      real(c_real_type), intent(in) :: vec1(lvec1) ! first vector
      integer(c_integer_type),intent(in) ::  lvec2        ! length of the second vector
      real(c_real_type), intent(in) :: vec2(lvec2) ! second vector - may be same as first
      
      ! result = vec1' * weights * vec2
      real(c_real_type), intent(out) :: dotprod

      ! local vars
      integer(c_integer_type),parameter :: ilevel = 1

      call bddcml_dotprod_subdomain( isub, vec1,lvec1, vec2,lvec2, dotprod )
end subroutine

!*******************************
subroutine bddcml_finalize_c() &
           bind(c)
!*******************************
! finalization of LEVELS module
      call bddcml_finalize
end subroutine

end module module_bddcml_c_wrapper
