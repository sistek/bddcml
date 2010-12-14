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
                                     inet,linet,nnet,lnnet,nndf,lnndf,xyz,lxyz1,lxyz2,&
                                     ifix,lifix, fixv,lfixv, rhs,lrhs, sol,lsol)
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

      call levels_upload_global_data(nelem,nnod,ndof,&
                                     inet,linet,nnet,lnnet,nndf,lnndf,xyz,lxyz1,lxyz2,&
                                     ifix,lifix,fixv,lfixv,rhs,lrhs,sol,lsol)

end subroutine

!************************************************************************************************
subroutine bddcml_setup_preconditioner(problemname, matrixtype, ndim, meshdim, neighbouring, &
                                       use_defaults, load_division,&
                                       parallel_division,correct_division,parallel_neighbouring,&
                                       parallel_globs,use_arithmetic,use_adaptive)
!************************************************************************************************
! setup multilevel preconditioner
      use module_levels
      use module_utils

      implicit none

      ! name of the problem
      character(*),intent(in) :: problemname

      ! numerical properties of the matrix (MUMPS-like notation)
      !  0 - general
      !  1 - symmetric positive definite
      !  2 - symmetric general
      integer,intent(in) :: matrixtype 

      ! space dimensions
      integer,intent(in) :: ndim 

      ! mesh dimension
      ! for 3D problems, = ndim
      ! shells = 2
      ! beams  = 1
      integer,intent(in) :: meshdim

      ! how many nodes should be shared for calling element adjacent in graph?
      integer,intent(in) :: neighbouring 

      ! use default values for all other variables?
      logical,intent(in) :: use_defaults

      ! use division from *.ES file?
      logical,intent(in) :: load_division

      ! should parallel division be used (ParMETIS instead of METIS)?
      logical,intent(in) :: parallel_division

      ! correct disconnected subdomains to make them continuous (not allowed for parallel divisions and loaded divisions)
      logical,intent(in) :: correct_division 

      ! should parallel search of neighbours be used? (distributed graph using ParMETIS rather than own serial graph)
      logical,intent(in) :: parallel_neighbouring

      ! should parallel search of globs be used? (some corrections on globs may not be available)
      logical,intent(in) :: parallel_globs

      ! use arithmetic constraints?
      logical,intent(in) :: use_arithmetic

      ! use adaptive constraints?
      logical,intent(in) :: use_adaptive

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
      character(*),parameter:: routine_name = 'BDDCML_SETUP_PRECONDITIONER'

      ! use globs from *.GLB and *.CN file?
      logical :: levels_load_globs = .false.
      ! load PAIRS from file?
      logical :: levels_load_pairs = .false.

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
      call levels_pc_setup(trim(problemname),levels_load_division,levels_load_globs,levels_load_pairs,&
                           levels_parallel_division,levels_correct_division,levels_parallel_neighbouring,neighbouring,&
                           levels_parallel_globs,matrixtype,ndim,meshdim,levels_use_arithmetic,levels_use_adaptive)

end subroutine


!*****************************************************************
subroutine bddcml_solve(problemname, comm_all,&
                        print_solution, write_solution_by_root, &
                        method,tol,maxit,ndecrmax)
!*****************************************************************
! solution of problem by a Krylov subspace iterative method
      use module_krylov
      use module_utils
      implicit none
      integer,parameter :: kr = kind(1.D0)

      ! name of the problem
      character(*),intent(in) :: problemname

      ! parallel variables
      integer,intent(in) :: comm_all 

      ! print solution on screen?
      logical, intent(in) :: print_solution

      ! write solution to a single file instead of distributed files?
      logical, intent(in) :: write_solution_by_root

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
         call krylov_bddcpcg(trim(problemname), comm_all,krylov_tol,krylov_maxit,krylov_ndecrmax,&
                             print_solution, write_solution_by_root)
      else if (krylov_method.eq.1) then 
         ! use BICGSTAB 
         call krylov_bddcbicgstab(trim(problemname), comm_all,krylov_tol,krylov_maxit,krylov_ndecrmax,&
                                  print_solution, write_solution_by_root)
      else
         call error(routine_name,'unknown iterative method',krylov_method)
      end if

end subroutine

!*************************
subroutine bddcml_finalize
!*************************
! finalization of LEVELS module
      use module_levels
      implicit none

      call levels_finalize
end subroutine

