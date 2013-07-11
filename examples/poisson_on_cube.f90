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

!**********************
program poisson_on_cube
!**********************
! '===========Possion on cube solver==========='
! ' Solves problem                             '
! '         -/\u = 1 on D = [0,1]^3,           '
! '            u = 0 on dD,                    '
! ' using FEM and the BDDCML solver.           '
! '============================================'
 
! utilities for error prints, etc.
      use module_utils
! functions for exporting to ParaView format
      use module_paraview

      implicit none
      
      include "mpif.h"

! precision of floats - depends on precision chosen for compiling the BDDCML library, do not change
      integer,parameter :: kr = kind(1.D0)

! **************************
! GENERAL BDDCML PARAMETERS:
! **************************

! beginning index of arrays ( 0 for C, 1 for Fortran )
      integer, parameter :: numbase = 1

! verbosity of BDDCML ( 0 - only fatal errors, 1 - mild output, 2 - detailed output )
      integer,parameter:: verbose_level = 0

! export solution to VTU files?
      logical,parameter :: export_solution = .true.

! what is the name of that file (resp. collection of files)
      character(*),parameter :: output_file_prefix = 'poisson_solution'

! *************************
! KRYLOV METHOD PARAMETERS:
! *************************

! Krylov subspace iterative method to be used
!     -1 - use solver defaults
!     0 - PCG
!     1 - BICGSTAB (choose for general symmetric and general matrices)
      integer,parameter :: krylov_method = 0  

! use recycling of Krylov subspace
!     0 - no recycling used
!     1 - basis of the Krylov subspace will be orthogonalized and also used for new right hand sides
      integer,parameter :: recycling_int = 1
! size of the Krylov subspace basis to store
      integer,parameter :: max_number_of_stored_vectors = 50

! maximum number of iterations of a Krylov subspace method
      integer,parameter :: maxit = 500

! maximum number of iterations of a Krylov subspace method with non-decreasing residual
      integer,parameter :: ndecrmax = 50

! relative precision of the Krylov subspace method ||residual||/||right-hand side||
      real(kr),parameter :: tol = 1.e-6_kr

! *******************************
! BDDC PRECONDITIONER PARAMETERS:
! *******************************

! use default values in preconditioner? In such case, all other parameters are ignored
      integer,parameter :: use_preconditioner_defaults = 0

! use arithmetic constraints on edges and faces?
      integer,parameter :: use_arithmetic_constraints = 1

! use adaptive constraints on faces?
      integer,parameter :: use_adaptive_constraints = 0

! use user constraints? - not used in this example
      integer,parameter :: use_user_constraints = 0

! what type of weights use on interface?
! 0 - weights by cardinality
! 1 - weights by diagonal stiffness
! 2 - weights based on first row of element data
! 3 - weights based on dof data
! 4 - weights by Marta Certikova - unit load
! 5 - weights by Marta Certikova - unit jump
      integer,parameter :: weights_type = 0

! should parallel division be used (ParMETIS instead of METIS) on the first level?
      integer,parameter :: parallel_division = 1

! *******************
! PROBLEM PARAMETERS:
! *******************
! numerical properties of the matrix (MUMPS-like notation)
!     0 - general (full storage)
!     1 - symmetric positive definite (only triangle stored)
!     2 - symmetric general (only triangle stored)
      integer,parameter :: matrixtype = 1  

! assuming tri-linear hexahedral finite elements
!     z
!   ^
!   |                                                                        
!   |                                                                        
!   5----------8            
!   |\         |\           
!   | \        | \          
!   |  \       |  \         
!   |   6----------7        
!   |   |      |   |        
!   1---|------4---|--> y   
!    \  |       \  |        
!     \ |        \ |        
!      \|         \|        
!       2----------3        
!        \
!         \
!         `'
!           x

! number of degrees of freedom on element
      integer,parameter :: ndof_per_element = 8

! element matrix on reference element [0,1]^3
      real(kr) :: element_matrix_ref(ndof_per_element * ndof_per_element) = &
      (/ 1._kr/3._kr , 0._kr       ,-1._kr/12._kr, 0._kr       , 0._kr       ,-1._kr/12._kr,-1._kr/12._kr,-1._kr/12._kr,&
         0._kr       , 1._kr/3._kr , 0._kr       ,-1._kr/12._kr,-1._kr/12._kr, 0._kr       ,-1._kr/12._kr,-1._kr/12._kr,&
        -1._kr/12._kr, 0._kr       , 1._kr/3._kr , 0._kr       ,-1._kr/12._kr,-1._kr/12._kr, 0._kr       ,-1._kr/12._kr,&
         0._kr       ,-1._kr/12._kr, 0._kr       , 1._kr/3._kr ,-1._kr/12._kr,-1._kr/12._kr,-1._kr/12._kr, 0._kr       ,&
         0._kr       ,-1._kr/12._kr,-1._kr/12._kr,-1._kr/12._kr, 1._kr/3._kr , 0._kr       ,-1._kr/12._kr, 0._kr       ,&
        -1._kr/12._kr, 0._kr       ,-1._kr/12._kr,-1._kr/12._kr, 0._kr       , 1._kr/3._kr , 0._kr       ,-1._kr/12._kr,&
        -1._kr/12._kr,-1._kr/12._kr, 0._kr       ,-1._kr/12._kr,-1._kr/12._kr, 0._kr       , 1._kr/3._kr , 0._kr       ,&
        -1._kr/12._kr,-1._kr/12._kr,-1._kr/12._kr, 0._kr       , 0._kr       ,-1._kr/12._kr, 0._kr       , 1._kr/3._kr /)

! spacial dimension
      integer,parameter :: ndim    = 3 

! topological dimension of elements elements, would be lower for shells or beams
      integer,parameter :: meshdim = 3 

      character(*),parameter:: routine_name = 'POISSON_ON_CUBE'

      ! input read from command line
      integer  :: num_el_per_sub_edge, num_sub_per_cube_edge ! basic properties of the cubes-in-cubes problem

      integer  :: num_el_per_cube_edge                        
      real(kr) :: hsize, el_vol                              ! elements size and volume

      !  parallel variables
      integer :: myid, comm_all, nproc, ierr

      integer :: nsub  ! number of subdomains on the first level 
      integer :: nelem ! number of elements 
      integer :: ndof  ! number of degrees of freedom 
      integer :: nnod  ! number of nodes

      integer :: nlevels ! number of levels
      ! subdomains in levels
      integer ::            lnsublev
      integer,allocatable :: nsublev(:)
      integer ::             nsub_loc_1
      ! bounds of subdomain numbers for each process
      integer ::            lsub2proc
      integer,allocatable::  sub2proc(:)

      ! scaled element matrix
      real(kr) :: element_matrix(ndof_per_element * ndof_per_element)

      ! local subdomain data
      integer :: nelems  ! subdomain number of elements
      integer :: ndofs   ! subdomain number on degrees of freedom
      integer :: nnods   ! subdomain number of nodes
      integer ::           linets,   lnnets,   lnndfs
      integer,allocatable:: inets(:), nnets(:), nndfs(:)
      integer ::           lxyzs1,   lxyzs2
      real(kr),allocatable:: xyzs(:,:)
      integer ::            lifixs
      integer,allocatable::  ifixs(:)
      integer ::            lfixvs
      real(kr),allocatable:: fixvs(:)
      integer ::            lrhss
      real(kr),allocatable:: rhss(:)
      integer ::            lsols
      real(kr),allocatable:: sols(:)
      integer ::           lisegns,   lisngns,   lisvgvns
      integer,allocatable:: isegns(:), isngns(:), isvgvns(:)

      ! matrix in coordinate format - triplets (i,j,a_ij)
      integer ::            la
      integer,allocatable::  i_sparse(:)
      integer,allocatable::  j_sparse(:)
      real(kr),allocatable:: a_sparse(:)

      ! user constraints - not really used here
      integer ::              luser_constraints1
      integer ::              luser_constraints2
      real(kr),allocatable ::  user_constraints(:)

      ! data for elements - not really used here
      integer ::              lelement_data1
      integer ::              lelement_data2
      real(kr),allocatable ::  element_data(:)

      ! data for dofs - not really used here
      integer ::              ldof_data
      real(kr),allocatable ::  dof_data(:)

      ! data about resulting convergence
      integer ::  num_iter, converged_reason 
      real(kr) :: condition_number
      real(kr) :: normRn_sol, normRn2, normRn2_loc, normRn2_sub 
      real(kr) :: normL2_sol, normL2_loc, normL2_sub 
      real(kr) :: normLinf_sol, normLinf_loc

      ! time variables
      real(kr) :: t_total, t_load, t_pc_setup, t_krylov

      ! small variables - indices, etc.
      integer :: ia, indinets, nne, idof, jdof, lelm
      integer :: ie, i, isub, j, ir
      integer :: is_rhs_complete_int
      integer :: is_assembled_int
      character(len=32) :: aux
      real(kr) :: coarsening
      character(len=256) :: command

      ! MPI initialization
!***************************************************************PARALLEL
      call MPI_INIT(ierr)
      ! Communicator
      comm_all  = MPI_COMM_WORLD
      call MPI_COMM_RANK(comm_all,myid,ierr)
      call MPI_COMM_SIZE(comm_all,nproc,ierr)
!***************************************************************PARALLEL

! Initial screen
      if (myid.eq.0) then
         write(*,'(a)') ' ===========Possion on cube solver=========== '
         write(*,'(a)') '| Solves problem                             |'
         write(*,'(a)') '|        -/\u = 1 on D = [0,1]^3,            |'
         write(*,'(a)') '|           u = 0 on dD,                     |'
         write(*,'(a)') '| using FEM and the BDDCML solver.           |'
         write(*,'(a)') ' ============================================ '
      end if

! Number of elements in an edge of a subdomain and number of subdomains in an edge of the unit cube
      if ( myid .eq. 0 ) then
         if(iargc().eq.3) then
            call getarg(1,aux)
            read(aux,*) num_el_per_sub_edge
            call getarg(2,aux)
            read(aux,*) num_sub_per_cube_edge
            call getarg(3,aux)
            read(aux,*) nlevels
         else
            write (*,'(a)') ' Usage: mpirun -np X ./poisson_on_cube NUM_EL_PER_SUB_EDGE NUM_SUB_PER_CUBE_EDGE NLEVELS'
            call error(routine_name,'trouble getting problem sizes')
         end if
      end if
! Broadcast of name of the problem      
!***************************************************************PARALLEL
      call MPI_BCAST(num_el_per_sub_edge,   1, MPI_INTEGER,   0, comm_all, ierr)
      call MPI_BCAST(num_sub_per_cube_edge, 1, MPI_INTEGER,   0, comm_all, ierr)
      call MPI_BCAST(nlevels,               1, MPI_INTEGER,   0, comm_all, ierr)
!***************************************************************PARALLEL

! measuring time
      call MPI_BARRIER(comm_all,ierr)
      call time_start

! number of elements on an edge of the unit cube
      num_el_per_cube_edge = num_el_per_sub_edge * num_sub_per_cube_edge
! element size
      hsize = 1._kr/ num_el_per_cube_edge
! element volume
      el_vol = hsize**ndim
! total number of elements
      nelem = num_el_per_cube_edge**ndim
! total number of subdomains
      nsub  = num_sub_per_cube_edge**ndim
! total number of nodes
      nnod  = (num_el_per_cube_edge+1)**ndim
! total number of degrees of freedom - equal to number of nodes for scalar problem
      ndof  = nnod

! initialize levels
      lnsublev = nlevels
      allocate(nsublev(lnsublev))
      if (nlevels.eq.2) then
         nsublev(1) = nsub
         nsublev(2) = 1
      else if (nlevels.gt.2) then
         ! determine coarsening factor
         coarsening = nsub**(1._kr/(nlevels-1))
         ! prescribe number of subdomains on levels so that coarsening is fixed between levels
         nsublev(1) = nsub
         do i = 2,nlevels-1
            ir = nlevels - i + 1
            nsublev(i) = int(coarsening**(ir-1))
            if (mod(nsublev(i),2).ne.0) then
               nsublev(i) = nsublev(i) + 1
            end if
         end do
         nsublev(nlevels) = 1
      else
         call error(routine_name,'Unsupported number of levels:',nlevels)
      end if

! Basic properties 
      if (myid.eq.0) then
         write(*,*)'Characteristics of the problem :'
         write(*,*)'  number of processors            nproc =',nproc
         write(*,*)'  number of dimensions             ndim =',ndim
         write(*,*)'  mesh dimension                meshdim =',meshdim
         write(*,*)'  number of elements global       nelem =',nelem
         write(*,*)'  number of subdomains             nsub =',nsub
         write(*,*)'  number of nodes global           nnod =',nnod
         write(*,*)'  number of DOF                    ndof =',ndof
         write(*,*)'  number of levels              nlevels =',nlevels
         write(*,*)'  number of subdomains in levels        =',nsublev
         write(*,*)'Characteristics of iterational process:'
         write(*,*)'  tolerance of error                tol =',tol
         write(*,*)'  maximum number of iterations    maxit =',maxit
         write(*,*)'  number of incresing residual ndecrmax =',ndecrmax
         write(*,*)'  using recycling of Krylov method ?     ',recycling_int
         call flush(6)
      end if

      if (myid.eq.0) then
         write (*,'(a)') 'Initializing BDDCML ...'
         call flush(6)
      end if
      ! tell me how much subdomains should I load
      nsub_loc_1 = -1
      call bddcml_init(nlevels, nsublev,lnsublev, nsub_loc_1, comm_all, verbose_level, numbase)
      if (myid.eq.0) then
         write (*,'(a)') 'Initializing BDDCML done.'
         call flush(6)
      end if

      lsub2proc = nproc + 1
      allocate(sub2proc(lsub2proc))
!***************************************************************PARALLEL
      call MPI_ALLGATHER( nsub_loc_1, 1, MPI_INTEGER, sub2proc(1), 1, MPI_INTEGER, comm_all, ierr)
!***************************************************************PARALLEL
      ! the array now contains counts, change it to starts
      do i = 2,nproc
         sub2proc(i) = sub2proc(i-1) + sub2proc(i)
      end do

      ! shift it one back and add one 
      do ir = 0, nproc - 1 ! reverse index
         i = nproc + 1 - ir

         sub2proc(i) = sub2proc(i-1) + 1
      end do
      ! put one in the beginning
      sub2proc(1) = 1

      ! create and load subdomains
      if (myid.eq.0) then
         write (*,'(a)') 'Loading data ...'
         call flush(6)
      end if
      call time_start
      ! Loop over subdomains and load them to BDDCML
      do isub = sub2proc(myid+1), sub2proc(myid+2) - 1

         ! create mesh for subdomain
         nelems   = num_el_per_sub_edge**ndim
         nnods    = (num_el_per_sub_edge+1)**ndim
         ndofs    = nnods
         linets   = nelems * 8
         lnnets   = nelems
         lnndfs   = nnods
         lisegns  = nelems
         lisngns  = nnods
         lisvgvns = ndofs
         allocate(inets(linets),nnets(lnnets),nndfs(lnndfs),isegns(lisegns),isngns(lisngns),isvgvns(lisvgvns))
         lxyzs1   = nnods
         lxyzs2   = ndim
         allocate(xyzs(lxyzs1,lxyzs2))
         lifixs   = nnods
         lfixvs   = nnods
         allocate(ifixs(lifixs),fixvs(lfixvs))

         ! create subdomain mesh and boundary conditions
         call prepare_subdomain_data(isub, num_sub_per_cube_edge, num_el_per_sub_edge, hsize, &
                                     inets,linets, nnets,lnnets, nndfs,lnndfs, &
                                     isegns,lisegns, isngns,lisngns, isvgvns,lisvgvns, &
                                     xyzs,lxyzs1,lxyzs2, &
                                     ifixs,lifixs, fixvs,lfixvs)

         ! create local right hand side
         lrhss = ndofs
         allocate(rhss(lrhss))
         rhss = 1._kr * el_vol
         ! erase right-hand side at the boundary
         where (ifixs.gt.0) rhss = 0._kr
         is_rhs_complete_int = 1

         ! create local initial solution
         lsols = ndofs
         allocate(sols(lsols))
         sols = 0._kr

         ! create local subdomain matrix for each subdomain
         ! full element matrix scaled based on element size
         element_matrix = hsize * element_matrix_ref 

         ! how much space the upper triangle of the element matrix occupies
         lelm = ndof_per_element * (ndof_per_element + 1) / 2
         ! space for all upper triangles of element matrics
         la   = nelems*lelm
         allocate(i_sparse(la), j_sparse(la), a_sparse(la))

         ! copy the upper triangle of the element matrix to the sparse triplet
         ia = 0
         indinets = 0
         do ie = 1,nelems
            nne = nnets(ie)
            do j = 1,ndof_per_element
               jdof = inets(indinets+j)
               do i = 1,j
                  idof = inets(indinets+i)
                  ia = ia + 1

                  if (idof.le.jdof) then
                     i_sparse(ia) = idof
                     j_sparse(ia) = jdof
                  else
                     ! transpose the entry
                     i_sparse(ia) = jdof
                     j_sparse(ia) = idof
                  end if
                  a_sparse(ia) = element_matrix((j-1)*ndof_per_element + i)

               end do
            end do
            indinets = indinets + nne
         end do
         is_assembled_int = 0

         ! prepare user constraints - not really used here
         luser_constraints1 = 0
         luser_constraints2 = 0
         allocate(user_constraints(luser_constraints1*luser_constraints2))

         ! prepare element data - not really used
         lelement_data1 = 0
         lelement_data2 = 0
         allocate(element_data(lelement_data1*lelement_data2))

         ! prepare dof data - not really used
         ldof_data = 0 
         allocate(dof_data(ldof_data))

         call bddcml_upload_subdomain_data(nelem, nnod, ndof, ndim, meshdim, &
                                           isub, nelems, nnods, ndofs, &
                                           inets,linets, nnets,lnnets, nndfs,lnndfs, &
                                           isngns,lisngns, isvgvns,lisvgvns, isegns,lisegns, &
                                           xyzs,lxyzs1,lxyzs2, &
                                           ifixs,lifixs, fixvs,lfixvs, &
                                           rhss,lrhss, is_rhs_complete_int, &
                                           sols,lsols, &
                                           matrixtype, i_sparse, j_sparse, a_sparse, la, is_assembled_int, &
                                           user_constraints,luser_constraints1,luser_constraints2, &
                                           element_data,lelement_data1,lelement_data2,&
                                           dof_data,ldof_data)

         deallocate(inets,nnets,nndfs,isegns,isngns,isvgvns)
         deallocate(xyzs)
         deallocate(ifixs,fixvs)
         deallocate(rhss)
         deallocate(sols)
         deallocate(i_sparse, j_sparse, a_sparse)
         deallocate(user_constraints)
         deallocate(element_data)
         deallocate(dof_data)
      end do
      call MPI_BARRIER(comm_all,ierr)
      call time_end(t_load)
      if (myid.eq.0) then
         write (*,'(a)') 'Loading data done.'
         call flush(6)
      end if

      if (myid.eq.0) then
         write (*,'(a)') 'Preconditioner set-up ...'
         call flush(6)
      end if
! PRECONDITIONER SETUP
      call MPI_BARRIER(comm_all,ierr)
      call time_start
      call bddcml_setup_preconditioner(matrixtype,&
                                       use_preconditioner_defaults, &
                                       parallel_division,&
                                       use_arithmetic_constraints,&
                                       use_adaptive_constraints,&
                                       use_user_constraints,&
                                       weights_type)
      call MPI_BARRIER(comm_all,ierr)
      call time_end(t_pc_setup)

      if (myid.eq.0) then
         write (*,'(a)') 'Preconditioner set-up done.'
         call flush(6)
      end if

      ! call Krylov method
      if (myid.eq.0) then
         write (*,'(a)') 'Calling Krylov method ...'
         call flush(6)
      end if
      call MPI_BARRIER(comm_all,ierr)
      call time_start
      ! call with setting of iterative properties
      call bddcml_solve(comm_all, krylov_method, tol,maxit,ndecrmax, recycling_int, max_number_of_stored_vectors, &
                        num_iter, converged_reason, condition_number)
      call MPI_BARRIER(comm_all,ierr)
      call time_end(t_krylov)
      if (myid.eq.0) then
         write (*,'(a)') 'Krylov method done.'
         call flush(6)
      end if
      if (myid.eq.0) then
         write(*,'(a)')          ' Output of PCG: =============='
         write(*,'(a,i7)')       ' Number of iterations: ', num_iter
         write(*,'(a,i1)')       ' Convergence reason:   ', converged_reason
         if ( condition_number .ge. 0._kr ) then 
            write(*,'(a,f11.3)') ' Condition number: ', condition_number
         end if
         write(*,'(a)')          ' ============================='
         call flush(6)
      end if

      normRn2_loc  = 0._kr
      normL2_loc   = 0._kr
      normLinf_loc = 0._kr
      if (export_solution) then
         ! make sure the directory exists
         if (myid.eq.0) then
            command = 'mkdir -p '//trim(output_file_prefix)
            call system(trim(command))
         end if
         call MPI_BARRIER(comm_all,ierr)
      end if
      do isub = sub2proc(myid+1), sub2proc(myid+2) - 1

         ! download local solution
         nnods    = (num_el_per_sub_edge+1)**ndim
         ndofs    = nnods
         lsols = ndofs
         allocate(sols(lsols))
         call bddcml_download_local_solution(isub, sols,lsols)

         ! compute norm of local solution
         call bddcml_dotprod_subdomain( isub, sols,lsols, sols,lsols, normRn2_sub )
         normRn2_loc = normRn2_loc + normRn2_sub

         ! re-create mesh for subdomain
         nelems   = num_el_per_sub_edge**ndim
         nnods    = (num_el_per_sub_edge+1)**ndim
         ndofs    = nnods
         linets   = nelems * 8
         lnnets   = nelems
         lnndfs   = nnods
         lisegns  = nelems
         lisngns  = nnods
         lisvgvns = ndofs
         allocate(inets(linets),nnets(lnnets),nndfs(lnndfs),isegns(lisegns),isngns(lisngns),isvgvns(lisvgvns))
         lxyzs1   = nnods
         lxyzs2   = ndim
         allocate(xyzs(lxyzs1,lxyzs2))
         lifixs   = nnods
         lfixvs   = nnods
         allocate(ifixs(lifixs),fixvs(lfixvs))

         call prepare_subdomain_data(isub, num_sub_per_cube_edge, num_el_per_sub_edge, hsize, &
                                     inets,linets, nnets,lnnets, nndfs,lnndfs, &
                                     isegns,lisegns, isngns,lisngns, isvgvns,lisvgvns, &
                                     xyzs,lxyzs1,lxyzs2, &
                                     ifixs,lifixs, fixvs,lfixvs)

         ! compute L_2 norm of the solution
         normL2_sub = 0._kr
         indinets = 0
         do ie = 1,nelems
            ! number of nodes on element
            nne = nnets(ie)
            normL2_sub = normL2_sub + el_vol * sum(sols(inets(indinets+1:indinets+nne))) / nne
            indinets = indinets + nne
         end do
         normL2_loc = normL2_loc + normL2_sub

         ! find maximum of the solution
         normLinf_loc = max(normLinf_loc,maxval(sols))

         if (export_solution) then
            call export_vtu_file_with_solution(output_file_prefix, &
                                               isub, nelems, nnods, &
                                               inets,linets, nnets,lnnets, &
                                               xyzs,lxyzs1,lxyzs2, &
                                               sols,lsols)
         end if
         deallocate(inets,nnets,nndfs,isegns,isngns,isvgvns)
         deallocate(xyzs)
         deallocate(ifixs,fixvs)
         deallocate(sols)
      end do
      if (export_solution) then
         ! export the umbrella PVD file
         if (myid.eq.0) then
            call paraview_export_pvd_file('poisson_solution', nsub)
         end if
      end if

      if (myid.eq.0) then
         write (*,'(a)') 'Finalizing BDDCML ...'
         call flush(6)
      end if
      call bddcml_finalize
      if (myid.eq.0) then
         write (*,'(a)') 'Finalizing BDDCML done.'
         call flush(6)
      end if

      ! find global R^n norm of solution
      call MPI_ALLREDUCE(normRn2_loc, normRn2, 1, MPI_DOUBLE_PRECISION, MPI_SUM, comm_all, ierr)
      normRn_sol = sqrt( normRn2 )
      if (myid.eq.0) then
      end if

      ! find global L_2 norm of solution
      call MPI_ALLREDUCE(normL2_loc, normL2_sol, 1, MPI_DOUBLE_PRECISION, MPI_SUM, comm_all, ierr)
      if (myid.eq.0) then
      end if

      ! find global L_inf norm of solution
      call MPI_ALLREDUCE(normLinf_loc, normLinf_sol, 1, MPI_DOUBLE_PRECISION, MPI_MAX, comm_all, ierr)
      if (myid.eq.0) then
      end if
      if (myid.eq.0) then
         write(*,'(a)') ' Solution properties========'
         write(*,'(a,f15.10)') ' L_2 norm:   ', normL2_sol
         write(*,'(a,f15.10)') ' L_inf norm: ', normLinf_sol
         write(*,'(a,f15.10)') ' R^n norm:   ', normRn_sol
         write(*,'(a)') ' ==========================='
      end if

      call MPI_BARRIER(comm_all,ierr)
      call time_end(t_total)

      ! Information about times
      if (myid.eq.0) then
         write(*,'(a)')         ' Profiling infirmation==========='
         write(*,'(a)')         ' TIMES OF RUN OF ROUTINES:'
         write(*,'(a,f11.3,a)') '  loading           ',t_load, ' s'
         write(*,'(a,f11.3,a)') '  pc_setup          ',t_pc_setup, ' s'
         write(*,'(a,f11.3,a)') '  Krylov method     ',t_krylov, ' s'
         write(*,'(a)')         '  _______________________________'
         write(*,'(a,f11.3,a)') '  total             ',t_total,    ' s'
         write(*,'(a)')         ' ================================'
      end if

      ! MPI finalization
!***************************************************************PARALLEL
      call MPI_FINALIZE(ierr)
!***************************************************************PARALLEL

      deallocate(nsublev)
      deallocate(sub2proc)

      end program

!************************************************************************************************
      subroutine prepare_subdomain_data(isub,num_sub_per_cube_edge, num_el_per_sub_edge, hsize, &
                                        inets,linets, nnets,lnnets, nndfs,lnndfs, &
                                        isegns,lisegns, isngns,lisngns, isvgvns,lisvgvns, &
                                        xyzs,lxyzs1,lxyzs2, &
                                        ifixs,lifixs, fixvs,lfixvs)
!************************************************************************************************
! subroutine creating data for one subdomain
      use module_utils

      implicit none
! precision of floats
      integer,parameter :: kr = kind(1.D0)

      integer, intent(in) :: isub                  ! global subdomain index
      integer, intent(in) :: num_sub_per_cube_edge ! number of subdomains in one edge of a cube 
      integer, intent(in) :: num_el_per_sub_edge   ! number of elements in one edge of a subdomain 
      real(kr), intent(in) :: hsize                ! element size

      integer, intent(in) :: linets 
      integer, intent(out) :: inets(linets)
      integer, intent(in) :: lnnets 
      integer, intent(out) :: nnets(lnnets)
      integer, intent(in) :: lnndfs 
      integer, intent(out) :: nndfs(lnndfs)
      integer, intent(in) :: lisegns
      integer, intent(out) :: isegns(lisegns)
      integer, intent(in) :: lisngns
      integer, intent(out) :: isngns(lisngns)
      integer, intent(in) :: lisvgvns
      integer, intent(out) :: isvgvns(lisvgvns)
      integer, intent(in) ::  lxyzs1,lxyzs2
      real(kr), intent(out) :: xyzs(lxyzs1,lxyzs2)
      integer, intent(in) :: lifixs
      integer, intent(out) :: ifixs(lifixs)
      integer, intent(in) ::  lfixvs
      real(kr), intent(out) :: fixvs(lfixvs)

      ! local vars
      character(*),parameter:: routine_name = 'PREPARE_SUBDOMAIN_DATA'
      integer ::  num_sub_xy
      integer ::  ind_sub_x, ind_sub_y, ind_sub_z
      integer ::  num_el_per_cube_edge
      integer ::  num_nodes_per_sub_edge, num_nodes_per_cube_edge
      integer ::  i, j, k
      integer ::  ig, jg, kg
      integer ::  indng, indns
      integer ::  indelg, indels
      integer ::  indinets
      integer ::  nne
      integer ::  n1, n2, n3, n4, n5, n6, n7, n8

      ! number of elements on one edge of the cube
      num_el_per_cube_edge = num_el_per_sub_edge*num_sub_per_cube_edge

      ! determine subdomain indices along coordinate axes
      ! intentional integer divisons
      num_sub_xy = num_sub_per_cube_edge**2

      ind_sub_z = (isub-1)/num_sub_xy + 1
      ind_sub_y = (isub - ((ind_sub_z-1) * num_sub_xy) - 1)/num_sub_per_cube_edge + 1
      ind_sub_x =  isub - ((ind_sub_z-1) * num_sub_xy) - ((ind_sub_y-1) * num_sub_per_cube_edge) 

      ! debug
      !write(*,*) 'subdomain index and coord indices:',isub, ind_sub_x, ind_sub_y, ind_sub_z

      ! initialize boundary conditions
      ifixs = 0
      fixvs = 0._kr

      ! number nodes and degrees of freedom
      num_nodes_per_sub_edge  = num_el_per_sub_edge + 1
      num_nodes_per_cube_edge = num_el_per_cube_edge + 1

      indns = 0
      do k = 1,num_nodes_per_sub_edge
         kg = (ind_sub_z-1)*(num_nodes_per_sub_edge-1) + k
         do j = 1,num_nodes_per_sub_edge
            jg = (ind_sub_y-1)*(num_nodes_per_sub_edge-1) + j
            do i = 1,num_nodes_per_sub_edge
               ig = (ind_sub_x-1)*(num_nodes_per_sub_edge-1) + i

               ! increase counter of local nodes
               indns = indns + 1
               ! compute global node index
               indng = ig + (jg-1)*num_nodes_per_cube_edge + (kg-1)*num_nodes_per_cube_edge**2

               isngns(indns) = indng

               ! compute coordinates 
               xyzs(indns,1) = (ig-1)*hsize
               xyzs(indns,2) = (jg-1)*hsize
               xyzs(indns,3) = (kg-1)*hsize

               ! for Poisson problem, there is only one dof per node, 
               nndfs(indns) = 1
               !and thus the numbering of nodes and dofs is the same,
               isvgvns(indns) = indng

               ! if node is on the boundary, fix boundary conditions
               if ( ig.eq.1 .or. ig.eq.num_nodes_per_cube_edge .or. &
                    jg.eq.1 .or. jg.eq.num_nodes_per_cube_edge .or. &
                    kg.eq.1 .or. kg.eq.num_nodes_per_cube_edge ) then

                    ifixs(indns) = 1
                    fixvs(indns) = 0._kr
               end if
            end do
         end do
      end do
      if (indns.ne.lisngns) then
         call error( routine_name, 'Some bug in node index computing for sub', isub )
      end if
      ! debug
      !write(*,*) 'isub',isub,'isngns',isngns


      ! create element connectivity
      indels = 0
      indinets = 0
      nne = 8
      do k = 1,num_el_per_sub_edge
         kg = (ind_sub_z-1)*num_el_per_sub_edge + k
         do j = 1,num_el_per_sub_edge
            jg = (ind_sub_y-1)*num_el_per_sub_edge + j
            do i = 1,num_el_per_sub_edge
               ig = (ind_sub_x-1)*num_el_per_sub_edge + i

               ! increase counter of local elements
               indels = indels + 1
               ! compute global element index
               indelg = ig + (jg-1)*num_el_per_cube_edge + (kg-1)*num_el_per_cube_edge**2

               ! compute local node index of the first node
               indns  = i + (j-1)*num_nodes_per_sub_edge + (k-1)*num_nodes_per_sub_edge**2

               ! compute indices of the eight nodes of each element
               n1 = indns
               n2 = indns + 1
               n3 = n2 + num_nodes_per_sub_edge
               n4 = n3 - 1
               n5 = n1 + num_nodes_per_sub_edge**2
               n6 = n2 + num_nodes_per_sub_edge**2
               n7 = n3 + num_nodes_per_sub_edge**2
               n8 = n4 + num_nodes_per_sub_edge**2

               inets(indinets + 1) = n1
               inets(indinets + 2) = n2
               inets(indinets + 3) = n3
               inets(indinets + 4) = n4
               inets(indinets + 5) = n5
               inets(indinets + 6) = n6
               inets(indinets + 7) = n7
               inets(indinets + 8) = n8

               indinets = indinets + nne

               ! number of nodes on element is constant for all elements
               nnets(indels) = nne

               ! embedding of local elements into global numbers
               isegns(indels) = indelg
            end do
         end do
      end do
      ! debug
      !write(*,*) 'isub',isub,'isegns',isegns
      !write(*,*) 'isub',isub,'inets',inets
      !write(*,*) 'isub',isub,'xyzs',xyzs
      !write(*,*) 'isub',isub,'ifixs',ifixs
      !write(*,*) 'isub',isub,'fixvs',fixvs

      end subroutine

!************************************************************************************************
      subroutine export_vtu_file_with_solution(prefix, &
                                               isub, nelems, nnods, &
                                               inets,linets, nnets,lnnets, &
                                               xyzs,lxyzs1,lxyzs2, &
                                               sols,lsols)
!************************************************************************************************
! subroutine for exporting a VTU file
      use module_utils
      use module_paraview

      implicit none
! precision of floats
      integer,parameter :: kr = kind(1.D0)

      character(*), intent(in) :: prefix           ! basename of vtu files
      integer, intent(in) :: isub                  ! global subdomain index
      integer, intent(in) :: nelems                ! subdomain number of elements
      integer, intent(in) :: nnods                 ! subdomain number of nodes

      integer, intent(in) :: linets 
      integer, intent(in) :: inets(linets)
      integer, intent(in) :: lnnets 
      integer, intent(in) :: nnets(lnnets)
      integer, intent(in) :: lxyzs1,lxyzs2
      real(kr), intent(in) :: xyzs(lxyzs1,lxyzs2)
      integer, intent(in) :: lsols
      real(kr), intent(in) :: sols(lsols)

      ! local vars
      character(*),parameter:: routine_name = 'EXPORT_VTU_FILE_WITH_SOLUTION'
      integer ::  idvtu
      integer ::             lsubdomain
      integer, allocatable :: subdomain(:)

      ! write solution to a separate VTU file
      call paraview_open_subdomain_file(prefix,isub,idvtu)

      ! write header of VTU file
      call paraview_write_mesh(idvtu, nelems,nnods, inets,linets, nnets,lnnets, xyzs,lxyzs1,lxyzs2)

      ! write cell data
      call paraview_open_celldata(idvtu)

      lsubdomain = nelems
      allocate( subdomain(lsubdomain) )
      subdomain = isub
      call paraview_write_dataarray(idvtu,1,'subdomain',subdomain,lsubdomain)
      deallocate( subdomain )
      call paraview_close_celldata(idvtu)

      ! write point data
      call paraview_open_pointdata(idvtu)
      ! export solution
      call paraview_write_dataarray(idvtu,1,'Solution',sols,lsols)
      call paraview_close_pointdata(idvtu)

      ! finalize the file
      call paraview_finalize_file(idvtu)
    
      ! close file
      call paraview_close_subdomain_file(idvtu)

      end subroutine
