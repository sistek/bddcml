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


! Interface providing Fortran functions that may be called from applications
! to exploit multilevel adaptive BDDC solver 
! Jakub Sistek, Praha 12/2010

!******************************************************************************************
subroutine bddcml_init(nl, nsublev,lnsublev, nsub_loc_1, comm_init, verbose_level, numbase)
!******************************************************************************************
! initialization of LEVELS module
      use module_levels
      use module_utils , only : suppress_output_on
      use module_krylov , only : krylov_set_profile_on
      implicit none

! given number of levels
      integer,intent(in) :: nl

! GLOBAL numbers of subdomains in all levels
      integer,intent(in) :: lnsublev
      integer,intent(in) ::  nsublev(lnsublev)

! LOCAL number of subdomains assigned to the process
!     >= 0 - number of local subdomains - sum up across processes to nsublev(1)
!     -1   - let solver decide, then the value is returned ( determining linear partition )
      integer, intent(inout):: nsub_loc_1 

! initial global communicator (possibly MPI_COMM_WORLD)
      integer, intent(in):: comm_init 

! level of verbosity
!     0 - only errors printed
!     1 - some output
!     2 - detailed output
      integer, intent(in):: verbose_level 

! first index of arrays ( 0 for C, 1 for Fortran )
      integer, intent(in) :: numbase

      call levels_init(nl,nsublev,lnsublev,nsub_loc_1,comm_init,verbose_level,numbase)

      select case ( verbose_level )
      case (0)
          call suppress_output_on
      case (1)
          continue ! default behaviour
      case (2)
          call krylov_set_profile_on
      end select


end subroutine

!**************************************************************************************
subroutine bddcml_upload_global_data(nelem,nnod,ndof,ndim,meshdim,&
                                     inet,linet,nnet,lnnet,nndf,lnndf,xyz,lxyz1,lxyz2,&
                                     ifix,lifix, fixv,lfixv, rhs,lrhs, sol,lsol, idelm, &
                                     neighbouring, load_division_int)
!**************************************************************************************
! Subroutine for loading global data as zero level
      use module_levels
      use module_utils
      implicit none
      integer,parameter :: kr = kind(1.D0)

      ! GLOBAL number of elements
      integer, intent(in):: nelem
      ! GLOBAL number of nodes
      integer, intent(in):: nnod
      ! GLOBAL number of degrees of freedom
      integer, intent(in):: ndof

      ! space dimensions
      integer,intent(in) :: ndim 

      ! mesh dimension
      ! for 3D problems, = ndim
      ! for 3D shells = 2
      ! for 3D beams  = 1
      integer,intent(in) :: meshdim


      ! GLOBAL Indices of Nodes on ElemenTs
      integer, intent(in):: linet
      integer, intent(in)::  inet(linet)

      ! GLOBAL Number of Nodes on ElemenTs
      integer, intent(in):: lnnet
      integer, intent(in)::  nnet(lnnet)

      ! GLOBAL Number of Nodal Degrees of Freedom
      integer, intent(in):: lnndf
      integer, intent(in)::  nndf(lnndf)

      ! GLOBAL Coordinates of nodes (X | Y | Z) or as one array (all X, all Y, all Z)
      integer, intent(in):: lxyz1, lxyz2
      real(kr), intent(in):: xyz(lxyz1,lxyz2)

      ! GLOBAL Indices of FIXed variables - all dof with Dirichlet BC are marked with its number
      integer, intent(in):: lifix
      integer, intent(in)::  ifix(lifix)

      ! GLOBAL FIXed Variables - where IFIX is nonzero, value of Dirichlet BC
      integer, intent(in):: lfixv
      real(kr), intent(in):: fixv(lfixv)

      ! GLOBAL Right-Hand Side
      integer, intent(in):: lrhs
      real(kr), intent(in):: rhs(lrhs)

      ! GLOBAL initial SOLution - initial approximation for iterative method
      integer, intent(in):: lsol
      real(kr), intent(in):: sol(lsol)

      ! opened unit with unformatted file with element matrices
      integer,intent(in) :: idelm  

      ! how many nodes should be shared for calling element adjacent in graph?
      integer,intent(in) :: neighbouring 

      ! use division from 'partition_l1.ES' file?
      integer,intent(in) :: load_division_int


      ! local variables
      logical :: load_division          

      call logical2integer(load_division_int,load_division)

      call levels_upload_global_data(nelem,nnod,ndof,ndim,meshdim,&
                                     inet,linet,nnet,lnnet,nndf,lnndf,xyz,lxyz1,lxyz2,&
                                     ifix,lifix,fixv,lfixv,rhs,lrhs, sol,lsol, idelm, &
                                     neighbouring,load_division)

end subroutine

!**************************************************************************************
subroutine bddcml_upload_subdomain_data(nelem, nnod, ndof, ndim, meshdim, &
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
                                        dof_data,ldof_data)
!**************************************************************************************
! Subroutine for loading global data as zero level
      use module_levels
      use module_utils
      implicit none
      integer,parameter :: kr = kind(1.D0)

      ! GLOBAL number of elements
      integer, intent(in):: nelem
      ! GLOBAL number of nodes
      integer, intent(in):: nnod
      ! GLOBAL number of degrees of freedom
      integer, intent(in):: ndof
      ! GLOBAL number of spacial dimensions
      integer, intent(in):: ndim
      ! mesh dimension
      ! for 3D problems, = ndim
      ! for 3D shells = 2
      ! for 3D beams  = 1
      integer,intent(in) :: meshdim

      ! GLOBAL index of subdomain
      integer, intent(in):: isub
      ! LOCAL number of elements
      integer, intent(in):: nelems
      ! LOCAL number of nodes
      integer, intent(in):: nnods
      ! LOCAL number of degrees of freedom
      integer, intent(in):: ndofs

      ! LOCAL Indices of Nodes on ElemenTs
      integer, intent(in):: linet
      integer, intent(in)::  inet(linet)

      ! LOCAL Number of Nodes on ElemenTs
      integer, intent(in):: lnnet
      integer, intent(in)::  nnet(lnnet)

      ! LOCAL Number of Nodal Degrees of Freedom
      integer, intent(in):: lnndf
      integer, intent(in)::  nndf(lnndf)

      ! Indices of Subdomain Nodes in Global Numbering (local to global map of nodes)
      integer, intent(in):: lisngn
      integer, intent(in)::  isngn(lisngn)

      ! Indices of Subdomain Variables in Global Variable Numbering (local to global map of variables)
      integer, intent(in):: lisvgvn
      integer, intent(in)::  isvgvn(lisvgvn)

      ! Indices of Subdomain Elements in Global Numbering (local to global map of elements)
      integer, intent(in):: lisegn
      integer, intent(in)::  isegn(lisegn)

      ! LOCAL Coordinates of nodes (X | Y | Z (if present) ) or as one array (all X, all Y, all Z) 
      ! lxyz1 = NNODS, lxyz2 = NDIM
      integer, intent(in):: lxyz1, lxyz2
      real(kr), intent(in):: xyz(lxyz1,lxyz2)

      ! LOCAL Indices of FIXed variables - all dof with Dirichlet BC are marked with its local number
      integer, intent(in):: lifix
      integer, intent(in)::  ifix(lifix)

      ! LOCAL FIXed Variables - where IFIX is nonzero, value of Dirichlet BC
      integer, intent(in):: lfixv
      real(kr), intent(in):: fixv(lfixv)

      ! LOCAL Right-Hand Side (where nodes of subdomains coincide, values are repeated)
      integer, intent(in):: lrhs
      real(kr), intent(in):: rhs(lrhs)

      ! LOCAL initial SOLution - initial approximation for iterative method
      integer, intent(in):: lsol
      real(kr), intent(in):: sol(lsol)
      integer, intent(in)::  is_rhs_complete_int  ! is the right-hand side complete?  
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
      integer, intent(in)::  matrixtype 
      integer, intent(in)::  la            ! length of previous arrays (= number of nonzeros for assembled matrix)
      integer, intent(in)::  i_sparse(la)  ! array of row indices
      integer, intent(in)::  j_sparse(la)  ! array of column indices
      real(kr), intent(in):: a_sparse(la)  ! array of values
      integer, intent(in)::  is_assembled_int  ! is the array assembled? 
                                               !  0 = no, it can contain repeated entries
                                               !  1 = yes, it is sorted and doesn't contain repeated index pairs
      ! LOCAL array of constraints defined by user
      integer, intent(in)::  luser_constraints1 ! number of rows in matrix of constraints, i.e. number of user constraints
      integer, intent(in)::  luser_constraints2 ! number of columns in matrix of constraints, ( = NNODS)
      real(kr), intent(in):: user_constraints(luser_constraints1*luser_constraints2) ! array for additional constraints

      ! LOCAL array of additional element data - can be used e.g. for construction of averaging operator E
      integer, intent(in)::  lelement_data1 ! number of rows in matrix of element data, i.e. number of data per element
      integer, intent(in)::  lelement_data2 ! number of columns in matrix of constraints, ( = NELEMS)
      real(kr), intent(in):: element_data(lelement_data1*lelement_data2) ! array for data on elements

      ! LOCAL array of additional data at degrees of freedom 
      integer, intent(in)::  ldof_data  ! number of entries in dof_data ( = NDOFS)
      real(kr), intent(in):: dof_data(ldof_data) ! array for data on degrees of freedom

      ! local vars
      logical :: is_assembled
      logical :: is_rhs_complete

      ! translate integer to logical
      call logical2integer(is_assembled_int,    is_assembled)
      call logical2integer(is_rhs_complete_int, is_rhs_complete)

      call levels_upload_subdomain_data(nelem, nnod, ndof, ndim, meshdim, &
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

end subroutine

!************************************************************************************************
subroutine bddcml_setup_preconditioner(matrixtype, use_defaults_int, &
                                       parallel_division_int, &
                                       use_arithmetic_constraints_int, &
                                       use_adaptive_constraints_int, &
                                       use_user_constraints_int, &
                                       weights_type)
!************************************************************************************************
! setup multilevel preconditioner
      use module_levels
      use module_utils

      implicit none

      ! numerical properties of the matrix (MUMPS-like notation)
      !  0 - general
      !  1 - symmetric positive definite
      !  2 - symmetric general
      integer,intent(in) :: matrixtype 

      ! use default values for all other settings?
      integer,intent(in) :: use_defaults_int

      ! should parallel division be used (ParMETIS instead of METIS)?
      integer,intent(in) :: parallel_division_int

      ! use arithmetic constraints?
      integer,intent(in) :: use_arithmetic_constraints_int

      ! use adaptive constraints?
      integer,intent(in) :: use_adaptive_constraints_int

      ! use user constraints?
      integer,intent(in) :: use_user_constraints_int

      ! what type of weights should be used on the interface ?
      ! 0 - weights by cardinality
      ! 1 - weights by diagonal stiffness
      ! 2 - weights based on first row of element data
      ! 3 - weights based on dof data
      ! 4 - weights by Marta Certikova - unit load
      ! 5 - weights by Marta Certikova - unit jump
      ! 6 - weights by Schur row sums for whole subdomain
      ! 7 - weights by Schur row sums computed face by face
      integer,intent(in) :: weights_type

      ! local variables
      ! ######################################
      ! default values of levels parameters
      logical :: levels_parallel_division          = .true.
      logical :: levels_use_arithmetic_constraints = .true.
      logical :: levels_use_adaptive_constraints   = .false. ! experimental - not used by default
      logical :: levels_use_user_constraints       = .false. ! experimental - not used by default
      integer :: levels_weights_type               = 0       ! arithmetic averaging - should always work
      ! ######################################
      logical :: use_defaults
      logical :: parallel_division
      logical :: use_arithmetic_constraints
      logical :: use_adaptive_constraints
      logical :: use_user_constraints

      character(*),parameter:: routine_name = 'BDDCML_SETUP_PRECONDITIONER'

      ! convert logical parameters (entered as integers for C compatibility)
      call logical2integer(use_defaults_int,use_defaults)

      if (use_defaults) then
         call info(routine_name,'using default preconditioner parameters')
      else
         call logical2integer(parallel_division_int,parallel_division)
         call logical2integer(use_arithmetic_constraints_int,use_arithmetic_constraints)
         call logical2integer(use_adaptive_constraints_int,use_adaptive_constraints)
         call logical2integer(use_user_constraints_int,use_user_constraints)
         levels_parallel_division      =  parallel_division    
         levels_use_arithmetic_constraints = use_arithmetic_constraints
         levels_use_adaptive_constraints   = use_adaptive_constraints
         levels_use_user_constraints       = use_user_constraints         
         levels_weights_type               = weights_type
      end if

      ! setup preconditioner
      call levels_pc_setup(levels_parallel_division,&
                           matrixtype,&
                           levels_use_arithmetic_constraints,&
                           levels_use_adaptive_constraints,&
                           levels_use_user_constraints,&
                           levels_weights_type)

end subroutine


!******************************************************************
subroutine bddcml_solve(comm_all,method,tol,maxit,ndecrmax, &
                        recycling_int, max_number_of_stored_vectors, &
                        num_iter,converged_reason,condition_number)
!******************************************************************
! solution of problem by a Krylov subspace iterative method
      use module_krylov
      use module_utils
      implicit none
      integer,parameter :: kr = kind(1.D0)

      ! parallel variables
      integer,intent(in) :: comm_all 

      ! optional paramaters - use either none or all of them
      ! preferable Krylov method
      ! -1 - use defaults (tol, maxit, and ndecrmax not accessed)
      ! 0 - PCG
      ! 1 - BICGSTAB
      integer, intent(in) :: method

      ! desired accuracy of relative residual
      real(kr), intent(in) :: tol

      ! limit on iterations
      integer, intent(in) :: maxit

      ! limit on iterations with increasing residual
      integer, intent(in) :: ndecrmax

      ! Should the Krylov subspace be recycled? If yes, basis will be stored and
      ! used for subsequent right hand sides. It also launches re-orthogonalization 
      ! of the Krylov basis.
      integer, intent(in) :: recycling_int

      ! If the Krylov subspace is recycled, what is the upper bound for the number of stored vectors?
      ! With increasing number, the memory consumption also increases.
      integer, intent(in) :: max_number_of_stored_vectors

      ! resulting number of iterations
      integer,intent(out) :: num_iter

      ! convergence reason
      !  =  0 - converged relative residual
      !  = -1 - reached limit on number of iterations
      !  = -2 - reached limit on number of iterations with nondecreasing residual
      integer,intent(out) :: converged_reason

      ! estimated condition number ( for PCG only )
      real(kr),intent(out) :: condition_number

      ! local variables
      ! ######################################
      ! default values of krylov parameters

      integer ::  krylov_method                       = 1
      real(kr) :: krylov_tol                          = 1.e-6_kr
      integer ::  krylov_maxit                        = 1000
      integer ::  krylov_ndecrmax                     = 30
      logical ::  krylov_recycling                    = .false.
      integer ::  krylov_max_number_of_stored_vectors = 100
      ! ######################################
      character(*),parameter:: routine_name = 'BDDCML_SOLVE'
      logical :: recycling

      call logical2integer(recycling_int,recycling)

      ! determine Krylov method and parameters
      if (method .ne. -1) then
         krylov_method    = method
         krylov_tol       = tol
         krylov_maxit     = maxit
         krylov_ndecrmax  = ndecrmax
         krylov_recycling = recycling
         krylov_max_number_of_stored_vectors = max_number_of_stored_vectors
      else
         call info(routine_name,'using default Krylov solver parameters')
      end if

      if (krylov_method.eq.0) then 
         ! use PCG 
         call krylov_bddcpcg(comm_all,krylov_tol,krylov_maxit,krylov_ndecrmax, &
                             krylov_recycling, krylov_max_number_of_stored_vectors, &
                             num_iter, converged_reason, condition_number)
      else if (krylov_method.eq.1) then 
         ! use BICGSTAB 
         call krylov_bddcbicgstab(comm_all,krylov_tol,krylov_maxit,krylov_ndecrmax, &
                                  num_iter, converged_reason)
         condition_number = -1._kr ! condition number is not computed for BICGSTAB
      else
         call error(routine_name,'unknown iterative method',krylov_method)
      end if

end subroutine

!**********************************************************
subroutine bddcml_download_local_solution(isub, sols,lsols)
!**********************************************************
! Subroutine for getting local solution,
! i.e. restriction of solution vector to subdomain (no weights are applied)
      use module_levels
      use module_utils
      implicit none
      integer,parameter :: kr = kind(1.D0)

      ! GLOBAL index of subdomain
      integer, intent(in)::  isub

      ! LOCAL solution 
      integer, intent(in)::  lsols
      real(kr), intent(out):: sols(lsols)

      call levels_dd_download_solution(isub, sols,lsols)

end subroutine

!****************************************************
subroutine bddcml_download_global_solution(sol, lsol)
!****************************************************
! Subroutine for getting global solution at root process
      use module_levels
      use module_utils
      implicit none
      integer,parameter :: kr = kind(1.D0)

      ! GLOBAL solution - required to be allocated only at root
      integer, intent(in)::  lsol
      real(kr), intent(out):: sol(lsol)

      call levels_get_global_solution(sol,lsol)

end subroutine

!***********************************************************
subroutine bddcml_download_local_reactions(isub, reas,lreas)
!***********************************************************
! Subroutine for getting local reactions at Dirichlet BC,
! i.e. restriction of vector of reactions to subdomain (no weights are applied)
      use module_levels
      use module_utils
      implicit none
      integer,parameter :: kr = kind(1.D0)

      ! GLOBAL index of subdomain
      integer, intent(in)::  isub

      ! LOCAL reactions 
      integer, intent(in)::  lreas
      real(kr), intent(out):: reas(lreas)

      call levels_dd_download_reactions(isub, reas,lreas)

end subroutine

!*****************************************************
subroutine bddcml_download_global_reactions(rea, lrea)
!*****************************************************
! Subroutine for getting global vector of reactions at Dirichlet BC at root process
      use module_levels
      use module_utils
      implicit none
      integer,parameter :: kr = kind(1.D0)

      ! GLOBAL  - required to be allocated only at root
      integer, intent(in)::  lrea
      real(kr), intent(out):: rea(lrea)

      call levels_get_global_reactions(rea,lrea)

end subroutine

!*******************************************************************************
subroutine bddcml_change_global_data(ifix,lifix, fixv,lfixv, rhs,lrhs, sol,lsol)
!*******************************************************************************
! Subroutine for changing global RHS, and values of Dirichlet boundary conditions
      use module_levels
      use module_utils
      implicit none
      integer,parameter :: kr = kind(1.D0)

      ! GLOBAL Indices of FIXed variables - all dof with Dirichlet BC are marked with its number
      integer, intent(in):: lifix
      integer, intent(in)::  ifix(lifix)

      ! GLOBAL FIXed Variables - where IFIX is nonzero, value of Dirichlet BC
      integer, intent(in):: lfixv
      real(kr), intent(in):: fixv(lfixv)

      ! GLOBAL Right-Hand Side
      integer, intent(in):: lrhs
      real(kr), intent(in):: rhs(lrhs)

      ! GLOBAL initial SOLution - initial approximation for iterative method
      integer, intent(in):: lsol
      real(kr), intent(in):: sol(lsol)

      call levels_change_global_data(ifix,lifix,fixv,lfixv,rhs,lrhs, sol,lsol)

end subroutine

!***********************************************************************
subroutine bddcml_change_subdomain_data(isub, &
                                        ifix,lifix, fixv,lfixv, &
                                        rhs,lrhs, is_rhs_complete_int, &
                                        sol,lsol)
!***********************************************************************
! Subroutine for loading global data as zero level
      use module_levels
      use module_utils
      implicit none
      integer,parameter :: kr = kind(1.D0)

      ! GLOBAL index of subdomain
      integer, intent(in):: isub

      ! LOCAL Indices of FIXed variables - all dof with Dirichlet BC are marked with its local number
      ! cannot change - just for reference
      integer, intent(in):: lifix
      integer, intent(in)::  ifix(lifix)

      ! LOCAL FIXed Variables - where IFIX is nonzero, value of Dirichlet BC
      integer, intent(in):: lfixv
      real(kr), intent(in):: fixv(lfixv)

      ! LOCAL Right-Hand Side (where nodes of subdomains coincide, values are repeated)
      integer, intent(in):: lrhs
      real(kr), intent(in):: rhs(lrhs)

      ! LOCAL initial SOLution - initial approximation for iterative method
      integer, intent(in):: lsol
      real(kr), intent(in):: sol(lsol)
      integer, intent(in)::  is_rhs_complete_int  ! is the right-hand side complete?  
                                                  !  0 = no, e.g. if it is assembled subdomain-by-subdomain 
                                                  !      and the interface entries are not interchanged
                                                  !  1 = yes, e.g. if it was created as a restriction of
                                                  !      the global RHS array and so entries on interface are
                                                  !      contained in several subdomains - weights will be applied

      ! local vars
      logical :: is_rhs_complete

      ! translate integer to logical
      call logical2integer(is_rhs_complete_int, is_rhs_complete)

      call levels_change_subdomain_data(isub, &
                                        ifix,lifix, fixv,lfixv, &
                                        rhs,lrhs, is_rhs_complete, &
                                        sol,lsol)
end subroutine

!*******************************
subroutine bddcml_setup_new_data
!*******************************
! finalization of LEVELS module
      use module_levels
      implicit none

      call levels_setup_new_data
end subroutine

!*******************************************************************************
subroutine bddcml_dotprod_subdomain( isub, vec1,lvec1, vec2,lvec2, dotprod )
!*******************************************************************************
! Auxiliary subroutine to compute scalar product of two vectors of lenght of
! subdomain exploiting interface weights from the solver. This routine is useful 
! if we want to compute global norm or dot product based on vectors restricted to 
! subdomains. Since interface values are contained in several vectors for
! several subdomains, this dot product or norm cannot be determined without
! weights.
      use module_levels
      implicit none
      integer,parameter :: kr = kind(1.D0)

      ! GLOBAL index of subdomain
      integer,intent(in) ::   isub 
      ! vectors to multiply
      integer,intent(in) ::  lvec1        ! length of the first vector
      real(kr), intent(in) :: vec1(lvec1) ! first vector
      integer,intent(in) ::  lvec2        ! length of the second vector
      real(kr), intent(in) :: vec2(lvec2) ! second vector - may be same as first
      
      ! result = vec1' * weights * vec2
      real(kr), intent(out) :: dotprod

      ! local vars
      integer,parameter :: ilevel = 1

      call levels_dd_dotprod_subdomain_local(ilevel, isub, vec1,lvec1, vec2,lvec2, dotprod)

end subroutine

!*************************
subroutine bddcml_finalize
!*************************
! finalization of LEVELS module
      use module_levels
      use module_krylov
      implicit none

      call levels_finalize
      call krylov_finalize
end subroutine

