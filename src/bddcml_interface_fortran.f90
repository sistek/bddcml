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

!****************************************************
subroutine bddcml_init(nl,nsublev,lnsublev,comm_init)
!****************************************************
! initialization of LEVELS module
      use module_levels
      implicit none

! given number of levels
      integer,intent(in) :: nl
! number of subdomains in all levels
      integer,intent(in) :: lnsublev
      integer,intent(in) ::  nsublev(lnsublev)
! initial global communicator (possibly MPI_COMM_WORLD)
      integer, intent(in):: comm_init 

      call levels_init(nl,nsublev,lnsublev,comm_init)

end subroutine

!**************************************************************************************
subroutine bddcml_upload_global_data(nelem,nnod,ndof,&
                                     numbase, inet,linet,nnet,lnnet,nndf,lnndf,xyz,lxyz1,lxyz2,&
                                     ifix,lifix, fixv,lfixv, rhs,lrhs, sol,lsol, idelm)
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

      ! beginning index of nodes ( 0 for C, 1 for Fortran )
      integer, intent(in) :: numbase

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

      call levels_upload_global_data(nelem,nnod,ndof,&
                                     numbase, inet,linet,nnet,lnnet,nndf,lnndf,xyz,lxyz1,lxyz2,&
                                     ifix,lifix,fixv,lfixv,rhs,lrhs, sol,lsol, idelm)

end subroutine

!**************************************************************************************
subroutine bddcml_upload_local_data(nelem, nnod, ndof, ndim, &
                                    isub, nelems, nnods, ndofs, &
                                    numbase, inet,linet, nnet,lnnet, nndf,lnndf, &
                                    isngn,lisngn, isvgvn,lisvgvn, isegn,lisegn, &
                                    xyz,lxyz1,lxyz2, &
                                    ifix,lifix, fixv,lfixv, &
                                    rhs,lrhs, &
                                    matrixtype, i_sparse, j_sparse, a_sparse, la, is_assembled_int)
!**************************************************************************************
! Subroutine for loading global data as zero level
! Only one subdomain is allowed to be loaded at each processor
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

      ! index of subdomain
      integer, intent(in):: isub
      ! LOCAL number of elements
      integer, intent(in):: nelems
      ! LOCAL number of nodes
      integer, intent(in):: nnods
      ! LOCAL number of degrees of freedom
      integer, intent(in):: ndofs

      ! beginning index of nodes ( 0 for C, 1 for Fortran )
      integer, intent(in) :: numbase

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
      ! not supported in present version
      !integer, intent(in):: lsol
      !real(kr), intent(in):: sol(lsol)

      ! LOCAL matrix triplet i, j, a(i,j) 
      ! Type of the matrix
      ! 0 - unsymmetric
      ! 1 - symmetric positive definite
      ! 2 - general symmetric
      integer, intent(in)::  matrixtype 
      integer, intent(in)::  i_sparse(la)  ! array of row indices
      integer, intent(in)::  j_sparse(la)  ! array of column indices
      real(kr), intent(in):: a_sparse(la)  ! array of values
      integer, intent(in)::  la            ! length of previous arrays (= number of nonzeros for assembled matrix)
      integer, intent(in)::  is_assembled_int  ! is the array assembled? 
                                               !  0 = no, it can contain repeated entries
                                               !  1 = yes, it is sorted and doesn't contain repeated index pairs


      ! local vars
      logical :: is_assembled

      ! translate integer to logical
      call logical2integer(is_assembled_int,is_assembled)

      call levels_upload_local_data(nelem, nnod, ndof, ndim, &
                                    isub, nelems, nnods, ndofs, &
                                    numbase, inet,linet, nnet,lnnet, nndf,lnndf, &
                                    isngn,lisngn, isvgvn,lisvgvn, isegn,lisegn, &
                                    xyz,lxyz1,lxyz2, &
                                    ifix,lifix, fixv,lfixv, &
                                    rhs,lrhs, &
                                    matrixtype, i_sparse, j_sparse, a_sparse, la, is_assembled)

end subroutine

!************************************************************************************************
subroutine bddcml_setup_preconditioner(matrixtype, ndim, meshdim, neighbouring, &
                                       use_defaults_int, load_division_int,&
                                       parallel_division_int,correct_division_int,parallel_neighbouring_int,&
                                       parallel_globs_int,use_arithmetic_int,use_adaptive_int)
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

      ! space dimensions
      integer,intent(in) :: ndim 

      ! mesh dimension
      ! for 3D problems, = ndim
      ! for 3D shells = 2
      ! for 3D beams  = 1
      integer,intent(in) :: meshdim

      ! how many nodes should be shared for calling element adjacent in graph?
      integer,intent(in) :: neighbouring 

      ! use default values for all other variables?
      integer,intent(in) :: use_defaults_int

      ! use division from *.ES file?
      integer,intent(in) :: load_division_int

      ! should parallel division be used (ParMETIS instead of METIS)?
      integer,intent(in) :: parallel_division_int

      ! correct disconnected subdomains to make them continuous (not allowed for parallel divisions and loaded divisions)
      integer,intent(in) :: correct_division_int

      ! should parallel search of neighbours be used? (distributed graph using ParMETIS rather than own serial graph)
      integer,intent(in) :: parallel_neighbouring_int

      ! should parallel search of globs be used? (some corrections on globs may not be available)
      integer,intent(in) :: parallel_globs_int

      ! use arithmetic constraints?
      integer,intent(in) :: use_arithmetic_int

      ! use adaptive constraints?
      integer,intent(in) :: use_adaptive_int

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

      ! local variables
      ! ######################################
      ! default values of levels parameters
      logical :: levels_load_division          = .false.
      logical :: levels_parallel_division      = .true.
      logical :: levels_correct_division       = .false.
      logical :: levels_parallel_neighbouring  = .true.
      logical :: levels_parallel_globs         = .true.
      logical :: levels_use_arithmetic         = .true.
      logical :: levels_use_adaptive           = .false. ! experimental - not used by default
      ! ######################################
      logical :: use_defaults
      logical :: load_division
      logical :: parallel_division
      logical :: correct_division
      logical :: parallel_neighbouring
      logical :: parallel_globs
      logical :: use_arithmetic
      logical :: use_adaptive

      character(*),parameter:: routine_name = 'BDDCML_SETUP_PRECONDITIONER'

      ! use globs from *.GLB and *.CN file?
      logical :: levels_load_globs = .false.
      ! load PAIRS from file?
      logical :: levels_load_pairs = .false.

      ! convert logical parameters (entered as integers for C compatibility)
      call logical2integer(use_defaults_int,use_defaults)
      call logical2integer(load_division_int,load_division)
      call logical2integer(parallel_division_int,parallel_division)
      call logical2integer(correct_division_int,correct_division)
      call logical2integer(parallel_neighbouring_int,parallel_neighbouring)
      call logical2integer(parallel_globs_int,parallel_globs)
      call logical2integer(use_arithmetic_int,use_arithmetic)
      call logical2integer(use_adaptive_int,use_adaptive)

      if (use_defaults) then
         call info(routine_name,'using default preconditioner parameters')
      else
         levels_load_division          =  load_division        
         levels_parallel_division      =  parallel_division    
         levels_correct_division       =  correct_division     
         levels_parallel_neighbouring  =  parallel_neighbouring
         levels_parallel_globs         =  parallel_globs       
         levels_use_arithmetic         =  use_arithmetic       
         levels_use_adaptive           =  use_adaptive         
      end if

      ! setup preconditioner
      call levels_pc_setup(levels_load_division,levels_load_globs,levels_load_pairs,&
                           levels_parallel_division,levels_correct_division,levels_parallel_neighbouring,neighbouring,&
                           levels_parallel_globs,matrixtype,ndim,meshdim,levels_use_arithmetic,levels_use_adaptive)

end subroutine


!******************************************************************
subroutine bddcml_solve(comm_all,method,tol,maxit,ndecrmax, &
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

      integer ::  krylov_method   = 1
      real(kr) :: krylov_tol      = 1.e-6_kr
      integer ::  krylov_maxit    = 1000
      integer ::  krylov_ndecrmax = 30
      ! ######################################
      character(*),parameter:: routine_name = 'BDDCML_SOLVE'

      ! determine Krylov method and parameters
      if (method .ne. -1) then
         krylov_method   = method
         krylov_tol      = tol
         krylov_maxit    = maxit
         krylov_ndecrmax = ndecrmax
      else
         call info(routine_name,'using default Krylov solver parameters')
      end if

      if (krylov_method.eq.0) then 
         ! use PCG 
         call krylov_bddcpcg(comm_all,krylov_tol,krylov_maxit,krylov_ndecrmax, &
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

!***************************************************************
subroutine bddcml_download_local_solution(sols, lsols, norm_sol)
!***************************************************************
! Subroutine for getting local solution,
! i.e. restriction of solution vector to subdomain (no weights are applied)
! Only one subdomain is allowed to be at each processor
      use module_levels
      use module_utils
      implicit none
      integer,parameter :: kr = kind(1.D0)

      ! LOCAL solution 
      integer, intent(in)::  lsols
      real(kr), intent(out):: sols(lsols)
      real(kr), intent(out):: norm_sol

      call levels_dd_download_solution(sols,lsols, norm_sol)

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

!*************************
subroutine bddcml_finalize
!*************************
! finalization of LEVELS module
      use module_levels
      implicit none

      call levels_finalize
end subroutine
