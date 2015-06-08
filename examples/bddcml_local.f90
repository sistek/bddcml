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
program bddcml_local
!***********************************************************************
! Example solver which uses multilevel BDDC library with LOCAL data loading
!
! basic bddcml module
      use module_bddcml
! module for preprocessing
      use module_pp
! module for sparse matrices
      use module_sm
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
      integer,parameter :: weights_type = 1

! beginning index of arrays ( 0 for C, 1 for Fortran )
      integer, parameter :: numbase = 1

! Just a direct solve by MUMPS?
      integer, parameter :: just_direct_solve_int = 0

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

! use recycling of Krylov subspace
      integer :: recycling_int = 1
      integer :: max_number_of_stored_vectors = 100

!######### END OF PARAMETERS TO SET
      character(*),parameter:: routine_name = 'BDDCML_LOCAL'

      !  parallel variables
      integer :: myid, comm_all, comm_self, nproc, ierr
      integer :: idpar, idml, idgmi, idfvs, idrhs, idelm, ides, idsols

      integer :: ndim, nsub, nelem, ndof, nnod, linet
      integer :: maxit, ndecrmax
      real(kr) :: tol

      integer(4) :: iargc

! number of levels
      integer :: nlevels
! subdomains in levels
      integer ::            lnsublev
      integer,allocatable :: nsublev(:)
      integer ::             nsub_loc_1

      integer ::            lsub2proc
      integer,allocatable::  sub2proc(:)
      integer ::            lndofsa
      integer,allocatable::  ndofsa(:)

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
      integer ::           liets
      integer,allocatable:: iets(:)
      integer ::            lkdof
      integer,allocatable::  kdof(:)

      ! local data
      integer :: nelems, ndofs, nnods
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
      integer ::            lkdofs
      integer,allocatable::  kdofs(:)


      ! matrix triplet
      integer ::            la
      integer,allocatable::  i_sparse(:)
      integer,allocatable::  j_sparse(:)
      real(kr),allocatable:: a_sparse(:)

      ! data not used here
      integer ::              luser_constraints1
      integer ::              luser_constraints2
      real(kr),allocatable ::  user_constraints(:)

      ! data for elements
      integer ::              lelement_data1
      integer ::              lelement_data2
      real(kr),allocatable ::  element_data(:)

      ! data for dofs
      integer ::              ldof_data
      real(kr),allocatable ::  dof_data(:)

      integer :: lproblemname
      integer :: meshdim

      character(lproblemnamex) :: problemname 
      character(lfilenamex)    :: filename

      ! small variables - indices, etc.
      integer :: ie, idofn, ind, indn, indvg, i, indvs, inod, inods, isub, j, ndofn, ir
      integer :: is_rhs_complete_int
      integer :: is_assembled_int

      ! data about resulting convergence
      integer :: num_iter, converged_reason 
      real(kr) :: condition_number
      real(kr) :: norm_sol, norm2, norm2_loc, norm2_sub 

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
         write(*,'(a)') '==================LOCAL====================='
      end if

! Name of the problem as the first argument
      problemname = ' '
      if ( myid .eq. 0 ) then
         if(iargc().eq.1) then
            call getarg(1,problemname)
         else
            write (*,'(a)') ' Usage: mpirun -np X ./bddcml_local PROBLEMNAME'
            call error(routine_name,'trouble getting problemname')
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
      lsol = ndof
      allocate(sol(lsol))
      sol = 0._kr

      ! read division into subdomains
      liets = nelem
      allocate(iets(liets))
      if (myid.eq.0) then
         filename = 'partition_l1.ES'
         call allocate_unit(ides)
         open (unit=ides,file=filename,status='old',form='formatted')
         rewind ides

         read(ides,*) iets
         close (ides)
      end if
      ! populate IETS along previous communicator
!***************************************************************PARALLEL
      call MPI_BCAST(iets,liets, MPI_INTEGER, 0, comm_all, ierr)
!***************************************************************PARALLEL

! create array kdof
      lkdof = nnod
      allocate(kdof(lkdof))
      if (nnod.gt.0) then
         kdof(1) = 0
         do inod = 2,nnod
            kdof(inod) = kdof(inod-1) + nndf(inod-1)
         end do
      end if

      if (myid.eq.0) then
         write (*,'(a)') 'Initializing LEVELS ...'
         call flush(6)
      end if
      call time_start
      ! tell me how much subdomains should I load
      nsub_loc_1 = -1
      call bddcml_init(nlevels, nsublev,lnsublev, nsub_loc_1, comm_all, verbose_level, numbase, &
                       just_direct_solve_int)
      call MPI_BARRIER(comm_all,ierr)
      call time_end(t_init)
      write (*,*) 'Initializing LEVELS done, locally owned subdomains: ', nsub_loc_1
      call flush(6)

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

      if (myid.eq.0) then
         write (*,'(a)') 'Loading data ...'
         call flush(6)
      end if
      call time_start
      lndofsa = nsub_loc_1
      allocate (ndofsa(lndofsa)) ! note number of degrees of freedom for each subdomain
      do isub = sub2proc(myid+1), sub2proc(myid+2) - 1
         ! localize mesh

         nelems = 0
         linets = 0
         do ie = 1,nelem
            if (iets(ie).eq.isub) then
               nelems = nelems + 1
               linets = linets + nnet(ie)
            end if
         end do
         lnnets  = nelems
         lisegns = nelems
         lisngns = linets
         allocate(inets(linets),nnets(lnnets),isegns(lisegns),isngns(lisngns))

         call pp_create_submesh(isub,nelem,inet,linet,nnet,lnnet,iets,liets,&
                                nnods,inets,linets,nnets,lnnets,&
                                isegns,lisegns,isngns,lisngns)

         ! correct lisngns
         lisngns = nnods
         lnndfs  = nnods
         allocate(nndfs(lnndfs))

         ! get array nndfs
         do inods = 1,nnods
            nndfs(inods) = nndf(isngns(inods))
         end do
! find local number of DOF on subdomain NDOFS
         ndofs = sum(nndfs)
         ndofsa(isub - sub2proc(myid+1) + 1) = ndofs

! create array kdofs
         lkdofs = nnods
         allocate(kdofs(lkdofs))
         if (nnods.gt.0) then
            kdofs(1) = 0
            do inods = 2,nnods
               kdofs(inods) = kdofs(inods-1) + nndfs(inods-1)
            end do
         end if

         ! prepare array ISVGVN
         lisvgvns = ndofs
         allocate(isvgvns(lisvgvns))
         indvs = 0
         do inods = 1,nnods
            inod = isngns(inods)

            ndofn   = nndfs(inods)
            indvg   = kdof(inod)
            ! in lack of memory for kdof array, this can be altered by indvg = sum(nndf(1:inod-1))
            do idofn = 1,ndofn
               indvs = indvs + 1
               indvg = indvg + 1

               isvgvns(indvs) = indvg
            end do
         end do

         ! coordinates of subdomain nodes
         lxyzs1 = nnods
         lxyzs2 = ndim
         allocate(xyzs(lxyzs1,lxyzs2))
         do j = 1,ndim
            do i = 1,nnods
               indn = isngns(i)

               xyzs(i,j) = xyz(indn,j)
            end do
         end do

         lrhss = ndofs
         allocate(rhss(lrhss))
         lifixs = ndofs
         allocate(ifixs(lifixs))
         lfixvs = ndofs
         allocate(fixvs(lfixvs))
         lsols = ndofs
         allocate(sols(lsols))

         ! localize RHS, IFIX and FIXV
         do i = 1,ndofs
            ind = isvgvns(i)

            rhss(i)  = rhs(ind)

            ifixs(i) = ifix(ind)
            fixvs(i) = fixv(ind)

            sols(i)  = sol(i)
         end do
         is_rhs_complete_int = 1

! prepare matrix
! Load sparse matrix
         ! find the length for matrix entries
         call sm_pmd_get_length(matrixtype,nelems,inets,linets,nnets,lnnets,nndfs,lnndfs,la)
         ! prepare memory
         allocate(i_sparse(la), j_sparse(la), a_sparse(la))
! ELMS - element stiffness matrices of subdomain - structure:
         call getfname(problemname(1:lproblemname),isub,'ELM',filename)
         call allocate_unit(idelm)
         open(unit=idelm,file=filename,status='old',form='unformatted')
         call sm_pmd_load(matrixtype,idelm,nelems,inets,linets,nnets,lnnets,nndfs,lnndfs,kdofs,lkdofs,&
                          i_sparse, j_sparse, a_sparse, la)


         is_assembled_int = 0

         ! prepare user constraints
         luser_constraints1 = 0
         luser_constraints2 = 0
         allocate(user_constraints(luser_constraints1*luser_constraints2))

         ! prepare element data
         lelement_data1 = 1
         lelement_data2 = nelems
         allocate(element_data(lelement_data1*lelement_data2))
         element_data = 2._kr

         ! prepare dof data
         ldof_data = ndofs 
         allocate(dof_data(ldof_data))
         dof_data = 3._kr


         ! experiment a bit
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
         deallocate(inets,nnets,nndfs,xyzs)
         deallocate(kdofs)
         deallocate(rhss)
         deallocate(ifixs,fixvs)
         deallocate(sols)
         deallocate(isvgvns)
         deallocate(isngns)
         deallocate(isegns)
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
         write (*,'(a)') 'Preconditioner SETUP ...'
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
         write (*,'(a)') 'Preconditioner SETUP done.'
         call flush(6)
      end if

      ! call PCG method
      call MPI_BARRIER(comm_all,ierr)
      call time_start
      ! call with setting of iterative properties
      if (krylov_method /= 5 .and. just_direct_solve_int /= 0) then
         call error( routine_name, 'You have to set just_direct_solve_int to nonzero if you want to use MUMPS')
      end if
      if (krylov_method == 5 .and. just_direct_solve_int == 0) then
         call error( routine_name, 'You have to set just_direct_solve_int to nonzero if you want to use MUMPS')
      end if
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

      norm2_loc = 0._kr
      do isub = sub2proc(myid+1), sub2proc(myid+2) - 1
         ! download local solution
         ndofs = ndofsa(isub - sub2proc(myid+1) + 1)
         lsols = ndofs
         allocate(sols(lsols))
         call bddcml_download_local_solution(isub, sols,lsols)

         ! write solution to separate file
         ! open subdomain SOLS file for solution
         call getfname(problemname(1:lproblemname),isub,'SOLS',filename)
         call info(routine_name,' Opening file fname: '//trim(filename))
         call allocate_unit(idsols)
         open (unit=idsols,file=trim(filename),status='replace',form='unformatted')

         rewind idsols
         write(idsols) (sols(i),i=1,lsols)
         
         if (print_solution) then
            write(*,*) 'isub =',isub,' solution: '
            write(*,'(e15.7)') sols
         end if
         close(idsols)

         if (krylov_method.ne.5) then
            call bddcml_dotprod_subdomain( isub, sols,lsols, sols,lsols, norm2_sub )
         end if

         norm2_loc = norm2_loc + norm2_sub

         deallocate(sols)
      end do

      ! find global norm of solution
      call MPI_ALLREDUCE(norm2_loc, norm2, 1, MPI_DOUBLE_PRECISION, MPI_SUM, comm_all, ierr)
      norm_sol = sqrt( norm2 )

      if (krylov_method.ne.5) then
         if (myid.eq.0) then
             write(*,*) 'Norm of solution is: ', norm_sol
         end if
      end if

      do isub = sub2proc(myid+1), sub2proc(myid+2) - 1
         ! localize mesh

         nelems = 0
         linets = 0
         do ie = 1,nelem
            if (iets(ie).eq.isub) then
               nelems = nelems + 1
               linets = linets + nnet(ie)
            end if
         end do
         lnnets  = nelems
         lisegns = nelems
         lisngns = linets
         allocate(inets(linets),nnets(lnnets),isegns(lisegns),isngns(lisngns))

         call pp_create_submesh(isub,nelem,inet,linet,nnet,lnnet,iets,liets,&
                                nnods,inets,linets,nnets,lnnets,&
                                isegns,lisegns,isngns,lisngns)

         ! correct lisngns
         lisngns = nnods
         lnndfs  = nnods
         allocate(nndfs(lnndfs))

         ! get array nndfs
         do inods = 1,nnods
            nndfs(inods) = nndf(isngns(inods))
         end do
! find local number of DOF on subdomain NDOFS
         ndofs = sum(nndfs)

! create array kdofs
         lkdofs = nnods
         allocate(kdofs(lkdofs))
         if (nnods.gt.0) then
            kdofs(1) = 0
            do inods = 2,nnods
               kdofs(inods) = kdofs(inods-1) + nndfs(inods-1)
            end do
         end if

         ! prepare array ISVGVN
         lisvgvns = ndofs
         allocate(isvgvns(lisvgvns))
         indvs = 0
         do inods = 1,nnods
            inod = isngns(inods)

            ndofn   = nndfs(inods)
            indvg   = kdof(inod)
            ! in lack of memory for kdof array, this can be altered by indvg = sum(nndf(1:inod-1))
            do idofn = 1,ndofn
               indvs = indvs + 1
               indvg = indvg + 1

               isvgvns(indvs) = indvg
            end do
         end do
         lrhss = ndofs
         allocate(rhss(lrhss))
         lifixs = ndofs
         allocate(ifixs(lifixs))
         lfixvs = ndofs
         allocate(fixvs(lfixvs))
         lsols = ndofs
         allocate(sols(lsols))

         ! localize RHS, IFIX and FIXV
         do i = 1,ndofs
            ind = isvgvns(i)

            rhss(i)  = rhs(ind)

            ifixs(i) = ifix(ind)
            fixvs(i) = fixv(ind)

            sols(i)  = sol(i)
         end do
         is_rhs_complete_int = 1

         ! experiment a bit
         call bddcml_change_subdomain_data(isub, &
                                           ifixs,lifixs, fixvs,lfixvs, &
                                           rhss,lrhss, is_rhs_complete_int, &
                                           sols,lsols)
         deallocate(inets,nnets,nndfs)
         deallocate(kdofs)
         deallocate(rhss)
         deallocate(ifixs,fixvs)
         deallocate(sols)
         deallocate(isvgvns)
         deallocate(isngns)
         deallocate(isegns)
      end do

      ! call PCG method
      call MPI_BARRIER(comm_all,ierr)
      call time_start
      ! call with setting of iterative properties
      if (krylov_method.ne.5) then
         call bddcml_setup_new_data
      end if
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

      norm2_loc = 0._kr
      do isub = sub2proc(myid+1), sub2proc(myid+2) - 1
         ! download local solution
         ndofs = ndofsa(isub - sub2proc(myid+1) + 1)
         lsols = ndofs
         allocate(sols(lsols))
         call bddcml_download_local_solution(isub, sols,lsols)

         ! write solution to separate file
         ! open subdomain SOLS file for solution
         call getfname(problemname(1:lproblemname),isub,'SOLS',filename)
         call info(routine_name,' Opening file fname: '//trim(filename))
         call allocate_unit(idsols)
         open (unit=idsols,file=trim(filename),status='replace',form='unformatted')

         rewind idsols
         write(idsols) (sols(i),i=1,lsols)
         
         if (print_solution) then
            write(*,*) 'isub =',isub,' solution: '
            write(*,'(e15.7)') sols
         end if
         close(idsols)

         if (krylov_method.ne.5) then
            call bddcml_dotprod_subdomain( isub, sols,lsols, sols,lsols, norm2_sub )
         end if

         norm2_loc = norm2_loc + norm2_sub

         deallocate(sols)
      end do

      ! find global norm of solution
      call MPI_ALLREDUCE(norm2_loc, norm2, 1, MPI_DOUBLE_PRECISION, MPI_SUM, comm_all, ierr)
      norm_sol = sqrt( norm2 )

      if (krylov_method.ne.5) then
         if (myid.eq.0) then
             write(*,*) 'Norm of solution is: ', norm_sol
         end if
      end if

      deallocate(nsublev)
      deallocate (ndofsa)
      deallocate(sub2proc)

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

      deallocate(inet,nnet,nndf,xyz)
      deallocate(iets)
      deallocate(kdof)
      deallocate(ifix,fixv)
      deallocate(rhs)
      deallocate(sol)

      ! MPI finalization
!***************************************************************PARALLEL
      call MPI_FINALIZE(ierr)
!***************************************************************PARALLEL
      end program

