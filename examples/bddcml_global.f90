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

!***********************************************************************
program bddcml_global
!***********************************************************************
! Example solver which uses multilevel BDDC library with GLOBAL data loading
!
! basic bddcml module
      use module_bddcml
! module for preprocessing
      use module_pp
! Program name
      use module_utils

      implicit none
      
      include "mpif.h"

!######### PARAMETERS TO SET
! precision of floats
      integer,parameter :: kr = kind(1.D0)

! numerical properties of the matrix (MUMPS-like notation)
!     0 - general (full storage)
!     1 - symmetric positive definite (only triangle stored)
!     2 - symmetric general (only triangle stored)
      integer :: matrixtype = 1 
! Krylov subspace iterative method to be used
!     -1 - use solver defaults
!     0 - PCG
!     1 - BICGSTAB (choose for general symmetric and general matrices)
!     2 - steepest descent method
!     5 - direct solve by MUMPS
      integer :: krylov_method = 0  

! use default values in preconditioner? In such case, all other parameters are ignored
      integer,parameter :: use_preconditioner_defaults = 0

! use arithmetic constraints?
      integer,parameter :: use_arithmetic_constraints = 1

! use adaptive constraints?
      integer,parameter :: use_adaptive_constraints = 0

! use user constraints?
      integer,parameter :: use_user_constraints = 0

! what type of weights use on interface?
      ! 0 - weights by cardinality
      ! 1 - weights by diagonal stiffness
      ! 2 - weights based on first row of element data
      ! 3 - weights based on dof data
      ! 4 - weights by Marta Certikova - unit load
      ! 5 - weights by Marta Certikova - unit jump
      ! 6 - weights by Schur row sums for whole subdomain
      ! 7 - weights by Schur row sums computed face by face
      integer,parameter :: weights_type = 4

! beginning index of arrays ( 0 for C, 1 for Fortran )
      integer, parameter :: numbase = 1

! Just a direct solve by MUMPS?
      integer, parameter :: just_direct_solve_int = 0

! use recycling of Krylov subspace
      integer :: recycling_int = 1
      integer :: max_number_of_stored_vectors = 500


! use prepared division into subdomains on first level in file *.ES?
      integer,parameter :: load_division = 0
! should parallel division be used (ParMETIS instead of METIS)?
      integer,parameter :: parallel_division = 1
! maximal length of problemname
      integer,parameter:: lproblemnamex = 100
! maximal length of any used file - should be reasonably larger than length of problem to allow suffices
      integer,parameter:: lfilenamex = 130
! verbosity ( 0 - only fatal errors, 1 - mild output, 2 - detailed output )
      integer,parameter:: verbose_level = 1
! print solution on screen?
      logical,parameter :: print_solution = .false.

!######### END OF PARAMETERS TO SET
      character(*),parameter:: routine_name = 'BDDCML_GLOBAL'

      !  parallel variables
      integer :: myid, comm_all, comm_self, nproc, ierr
      integer :: idpar, idml, idgmi, idfvs, idrhs, idsol

      integer :: i

      integer(4) :: iargc

      integer :: idelm, ios

      integer :: ndim, nsub, nelem, ndof, nnod, linet
      integer :: maxit, ndecrmax
      real(kr) :: tol

! number of levels
      integer :: nlevels
! subdomains in levels
      integer ::            lnsublev
      integer,allocatable :: nsublev(:)

      ! number of local subdomains on level 1
      integer :: nsub_loc_1

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
      integer ::            linisol
      real(kr),allocatable:: inisol(:)
      integer ::            lsol
      real(kr),allocatable:: sol(:)
      integer ::            lrea
      real(kr),allocatable:: rea(:)
      real(kr) :: sum_rea(3)

      integer :: idofn, ndofn, inddof, inod


      integer :: lproblemname
      integer :: meshdim
      character(2) :: aux
      integer :: neighbouring

      character(lproblemnamex) :: problemname 
      character(lfilenamex)    :: filename

      ! data about resulting convergence
      integer :: num_iter, converged_reason 
      real(kr) :: condition_number

      ! time variables
      real(kr) :: t_total, t_import, t_distribute, t_init, t_load, t_pc_setup, t_pcg


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
         write(*,'(a)') '==================GLOBAL===================='
      end if

! Name of the problem as the first argument, 
! number of nodes shared by two elements to call them adjacent as the second argument
      problemname = ' '
      if ( myid .eq. 0 ) then
         if(iargc().eq.2) then
            call getarg(1,problemname)
            call getarg(2,aux)
            read ( aux, * ) neighbouring
         else
            write (*,'(a)') ' Usage: mpirun -np X ./bddcml_local PROBLEMNAME M              '
            write (*,'(a)') '  M - minimal number of shared nodes to call elements adjacent '
            call error(routine_name,'trouble getting problemname and neighbouring')
         end if
         ! get length
         lproblemname = index(problemname,' ') - 1
         if (lproblemname.eq.-1) then
            lproblemname = lproblemnamex 
         end if
         ! pad the name with spaces
         do i = lproblemname+1,lproblemnamex
            problemname(i:i) = ' '
         end do
      end if
! Broadcast of name of the problem      
!***************************************************************PARALLEL
      call MPI_BCAST(lproblemname, 1,           MPI_INTEGER,   0, comm_all, ierr)
      call MPI_BCAST(problemname, lproblemname, MPI_CHARACTER, 0, comm_all, ierr)
      call MPI_BCAST(neighbouring, 1,           MPI_INTEGER,   0, comm_all, ierr)
!***************************************************************PARALLEL

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
      if (matrixtype.ne.1.and.krylov_method .eq. 0) then
         if (myid.eq.0) then
            call warning(routine_name,'PCG method is not considered robust for non SPD matrices')
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
      linisol = ndof
      allocate(inisol(linisol))
      inisol = 0._kr

      if (myid.eq.0) then
         write (*,'(a)') 'Initializing LEVELS ...'
         call flush(6)
      end if

      ! attach file with element matrices
      if (myid.eq.0) then
         ! ELM - element stiffness matrices - structure:
         call allocate_unit(idelm)
         open(unit=idelm,file=problemname(1:lproblemname)//'.ELM',status='old',form='unformatted',iostat=ios)
         if (ios.ne.0) then
            call error(routine_name,'Problem opening file '//trim(problemname)//'.ELM')
         end if
         rewind idelm
      end if

      call time_start
      nsub_loc_1 = -1
      call bddcml_init(nlevels, nsublev,lnsublev, nsub_loc_1, comm_all, verbose_level, numbase, &
                       just_direct_solve_int)
      call MPI_BARRIER(comm_all,ierr)
      call time_end(t_init)
      if (myid.eq.0) then
         write (*,'(a)') 'Initializing LEVELS done.'
         call flush(6)
      end if
      if (myid.eq.0) then
         write (*,'(a)') 'Loading data ...'
         call flush(6)
      end if
      call time_start
      call bddcml_upload_global_data(nelem,nnod,ndof,ndim, meshdim, &
                                     inet,linet,nnet,lnnet,nndf,lnndf,xyz,lxyz1,lxyz2,&
                                     ifix,lifix,fixv,lfixv,rhs,lrhs,inisol,linisol, idelm, neighbouring,load_division)
      call MPI_BARRIER(comm_all,ierr)
      call time_end(t_load)
      if (myid.eq.0) then
         write (*,'(a)') 'Loading data done.'
         call flush(6)
      end if
      deallocate(inet,nnet,xyz)

      if (myid.eq.0) then
         write (*,'(a)') 'Preconditioner SETUP ...'
         call flush(6)
      end if

! PRECONDITIONER SETUP
      call MPI_BARRIER(comm_all,ierr)
      call time_start
      call bddcml_setup_preconditioner(matrixtype, &
                                       use_preconditioner_defaults,&
                                       parallel_division,&
                                       use_arithmetic_constraints,&
                                       use_adaptive_constraints,&
                                       use_user_constraints,&
                                       weights_type)
      call MPI_BARRIER(comm_all,ierr)
      call time_end(t_pc_setup)

      if (myid.eq.0) then
         write (*,'(a)') 'Preconditioner SETUP done.'
         call flush(6)
      end if

      ! call PCG method
      call MPI_BARRIER(comm_all,ierr)
      call time_start
      ! call with setting of iterative properties
      call bddcml_solve(comm_all, krylov_method, tol,maxit,ndecrmax, recycling_int, max_number_of_stored_vectors, &
                        num_iter, converged_reason, condition_number)
      if (myid.eq.0) then
          write(*,*) 'Number of iterations: ', num_iter
          write(*,*) 'Convergence reason: ', converged_reason
          if ( condition_number .ge. 0._kr ) then 
             write(*,*) 'Condition number: ', condition_number
          end if
      end if
      call MPI_BARRIER(comm_all,ierr)
      call time_end(t_pcg)

      deallocate(nsublev)

      if (myid.eq.0) then
         lsol = ndof
         lrea = ndof
      else 
         lsol = 0
         lrea = 0
      end if
      allocate(rea(lrea))
      allocate(sol(lsol))

      ! download global solution - all processors have to call this
      call bddcml_download_global_solution(sol,lsol)

      ! download global reactions - all processors have to call this
      call bddcml_download_global_reactions(rea,lrea)

      if (myid.eq.0) then
         ! compute sums of reactions into coordinate axes
         inddof = 0
         sum_rea = 0._kr
         do inod = 1,nnod
            ndofn = nndf(inod)
            do idofn = 1,ndofn
               inddof = inddof + 1
               sum_rea(idofn) = sum_rea(idofn) + rea(inddof)
            end do
         end do

         ! solution is ready, write it to SOL file
         call allocate_unit(idsol)
         open (unit=idsol,file=trim(problemname)//'.SOL',status='replace',form='unformatted')
         rewind idsol
      
         ! solution
         write(idsol) (sol(i), i = 1,lsol)
         ! reactions
         write(idsol) (rea(i), i = 1,lrea), (sum_rea(i), i = 1,3)
         write(*,*) 'Solution has been written into file ',trim(problemname)//'.SOL'
         close(idsol)

         if (print_solution) then
            write(*,*) ' solution | reactions'
            write(*,'(2e15.7)') (sol(i), rea(i), i = 1,lsol)
            write(*,*) ' sums of reactions into coordinate axes: ',sum_rea
         end if

      end if
      deallocate(nndf)
      deallocate(sol)
      deallocate(rea)

      ! load the data again
      call bddcml_change_global_data(ifix,lifix,fixv,lfixv,rhs,lrhs,inisol,linisol)
      call bddcml_setup_new_data

      ! call PCG method
      call MPI_BARRIER(comm_all,ierr)
      call time_start
      ! call with setting of iterative properties
      call bddcml_solve(comm_all, krylov_method, tol,maxit,ndecrmax, recycling_int, max_number_of_stored_vectors, &
                        num_iter, converged_reason, condition_number)
      if (myid.eq.0) then
          write(*,*) 'Number of iterations: ', num_iter
          write(*,*) 'Convergence reason: ', converged_reason
          if ( condition_number .ge. 0._kr ) then 
             write(*,*) 'Condition number: ', condition_number
          end if
      end if
      call MPI_BARRIER(comm_all,ierr)
      call time_end(t_pcg)


      if (myid.eq.0) then
         write (*,'(a)') 'Finalizing LEVELS ...'
         call flush(6)
      end if
      call bddcml_finalize
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
         write(*,'(a,f11.3,a)') '  loading           ',t_load, ' s'
         write(*,'(a,f11.3,a)') '  pc_setup          ',t_pc_setup, ' s'
         write(*,'(a,f11.3,a)') '  Krylov method     ',t_pcg, ' s'
         write(*,'(a)')         '  ______________________________'
         write(*,'(a,f11.3,a)') '  total             ',t_total,    ' s'
      end if

      ! close *.ELM file
      if (myid.eq.0) then
         close(idelm)
      end if

      deallocate(ifix,fixv)
      deallocate(rhs)
      deallocate(inisol)

      ! MPI finalization
!***************************************************************PARALLEL
      call MPI_FINALIZE(ierr)
!***************************************************************PARALLEL
      end program

