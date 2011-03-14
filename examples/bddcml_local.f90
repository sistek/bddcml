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
! Example solver which uses multilevel BDDC library with GLOBAL data loading
!
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
      integer :: krylov_method = 0  

! use default values in preconditioner? In such case, all other parameters are ignored
      integer,parameter :: use_preconditioner_defaults = 0

! use arithmetic constraints?
      integer,parameter :: use_arithmetic = 1

! use adaptive constraints?
      integer,parameter :: use_adaptive = 0

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

! use prepared division into subdomains on first level in file *.ES?
      integer,parameter :: load_division = 1
! use prepared selection of corners in file *.CN and description of globs for first level in file *.GLB?
      integer,parameter :: load_globs = 0
! use prepared file with pairs for adaptivity (*.PAIR) on first level?
      integer,parameter :: load_pairs = 0
! should parallel division be used (ParMETIS instead of METIS)?
      integer,parameter :: parallel_division = 1
! correct disconnected subdomains to make them continuous (not allowed for parallel divisions and loaded divisions)
      integer,parameter :: correct_division = 1
! should parallel search of neighbours be used? (distributed graph rather than serial graph)
      integer,parameter :: parallel_neighbouring = 1
! should parallel search of globs be used? (some corrections on globs may not be available)
      integer,parameter :: parallel_globs = 1
! are you debugging the code?
      integer,parameter :: debug = 1
! maximal length of problemname
      integer,parameter:: lproblemnamex = 100
! maximal length of any used file - should be reasonably larger than length of problem to allow suffices
      integer,parameter:: lfilenamex = 130
! print solution on screen?
      integer,parameter :: print_solution = 0
! write solution to a single file instead of distributed files?
      integer,parameter :: write_solution_by_root = 1

!######### END OF PARAMETERS TO SET
      character(*),parameter:: routine_name = 'BDDCML'

      !  parallel variables
      integer :: myid, comm_all, comm_self, nproc, ierr
      integer :: idpar, idml, idgmi, idfvs, idrhs, idelm, ides

      integer :: ndim, nsub, nelem, ndof, nnod, linet
      integer :: maxit, ndecrmax
      real(kr) :: tol

! number of levels
      integer :: nlevels
! subdomains in levels
      integer ::            lnsublev
      integer,allocatable :: nsublev(:)

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
      integer :: nelems, ndofs, nnods, nadjs
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
      integer ::            lkadjsub
      integer,allocatable::  kadjsub(:)
      integer ::            liadjs
      integer,allocatable::  iadjs(:)
      integer ::            lkdofs
      integer,allocatable::  kdofs(:)


      ! matrix triplet
      integer ::            la
      integer,allocatable::  i_sparse(:)
      integer,allocatable::  j_sparse(:)
      real(kr),allocatable:: a_sparse(:)

      integer :: lproblemname
      integer :: meshdim
      integer :: neighbouring

      character(lproblemnamex) :: problemname 
      character(lfilenamex)    :: filename

      ! small variables - indices, etc.
      integer :: ie, idofn, ind, indn, indvg, ia, i, indvs, inod, inods, isub, j, ndofn
      integer :: is_assembled_int

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
      end if

! Name of the problem
      call pp_pget_problem_name(comm_all,problemname,lproblemnamex,lproblemname)

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
      if (nsub.ne.nproc) then
         if (myid.eq.0) then
            call error(routine_name,'This solver can run solely on number of subdomains'&
                                    //'equal to the number of processors.')
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
         filename = trim(problemname)//'.ES'
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

      ! localize mesh
      isub = myid + 1

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

! prepare matrix
! Load sparse matrix
      ! find the length for matrix entries
      call sm_pmd_get_length(matrixtype,nelems,inets,linets,nnets,lnnets,nndfs,lnndfs,la)
      ! prepare memory
      allocate(i_sparse(la), j_sparse(la), a_sparse(la))
! ELMS - element stiffness matrices of subdomain - structure:
      call getfname(problemname(1:lproblemname),isub,'ELM',filename)
      open(unit=idelm,file=filename,status='old',form='unformatted')
      call sm_pmd_load(idelm,nelems,inets,linets,nnets,lnnets,nndfs,lnndfs,kdofs,lkdofs,&
                       i_sparse, j_sparse, a_sparse, la)


      if (myid.eq.0) then
         write (*,'(a)') 'Initializing LEVELS ...'
         call flush(6)
      end if

      call time_start
      call bddcml_init(nlevels,nsublev,lnsublev,comm_all)
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
      is_assembled_int = 0

      if (myid.eq.0) then
         write (*,'(a)') 'Minimal number of shared nodes to call elements adjacent: '
         call flush(6)
         read (*,*) neighbouring
      end if

      lkadjsub = nsub * nsub
      allocate(kadjsub(lkadjsub))
      kadjsub = 0
      call pp_get_sub_neighbours(1,nelem,nnod,nsub, &
                                 inet,linet,nnet,lnnet,iets,liets,&
                                 kadjsub,lkadjsub)
      ! parse kadjsub
      nadjs = count(kadjsub((isub-1)*nsub + 1 : isub*nsub).eq.1)
      liadjs = nadjs
      allocate(iadjs(liadjs))
      ia = 0
      do i = 1,nsub
         if (kadjsub((isub-1)*nsub + i).eq.1) then
            ia = ia + 1

            iadjs(ia) = i
         end if
      end do
      deallocate(kadjsub)

      call bddcml_upload_local_data(nelem, nnod, ndof, ndim, &
                                    isub, nelems, nnods, ndofs, &
                                    inets,linets, nnets,lnnets, nndfs,lnndfs, &
                                    isngns,lisngns, isvgvns,lisvgvns, isegns,lisegns, &
                                    xyzs,lxyzs1,lxyzs2, &
                                    ifixs,lifixs, fixvs,lfixvs, &
                                    rhss,lrhss, &
                                    nadjs, iadjs,liadjs, &
                                    matrixtype, i_sparse, j_sparse, a_sparse, la, is_assembled_int)
      call MPI_BARRIER(comm_all,ierr)
      call time_end(t_load)
      if (myid.eq.0) then
         write (*,'(a)') 'Loading data done.'
         call flush(6)
      end if
      deallocate(inet,nnet,nndf,xyz)
      deallocate(kdof)
      deallocate(ifix,fixv)
      deallocate(rhs)
      deallocate(sol)
      deallocate(inets,nnets,nndfs,xyzs)
      deallocate(kdofs)
      deallocate(rhss)
      deallocate(ifixs,fixvs)
      deallocate(sols)
      deallocate(isvgvns)
      deallocate(isngns)
      deallocate(isegns)
      deallocate(iadjs)
      deallocate(i_sparse, j_sparse, a_sparse)

! Broadcast basic properties of the problem
!***************************************************************PARALLEL
      call MPI_BCAST(neighbouring,1,MPI_INTEGER,      0, comm_all, ierr)
!***************************************************************PARALLEL

      if (myid.eq.0) then
         write (*,'(a)') 'Preconditioner SETUP ...'
         call flush(6)
      end if
! PRECONDITIONER SETUP
      call MPI_BARRIER(comm_all,ierr)
      call time_start
      call bddcml_setup_preconditioner(problemname(1:lproblemname), matrixtype, ndim, meshdim, neighbouring, &
                                       use_preconditioner_defaults, load_division,&
                                       parallel_division,correct_division,parallel_neighbouring,&
                                       parallel_globs,use_arithmetic,use_adaptive)
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
      call bddcml_solve(problemname(1:lproblemname), comm_all,&
                        print_solution, write_solution_by_root, &
                        krylov_method, tol,maxit,ndecrmax)
      call MPI_BARRIER(comm_all,ierr)
      call time_end(t_pcg)


      deallocate(nsublev)

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

      ! MPI finalization
!***************************************************************PARALLEL
      call MPI_FINALIZE(ierr)
!***************************************************************PARALLEL
      end program

