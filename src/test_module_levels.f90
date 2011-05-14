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

program test_module_levels
! Tester of module_levels
      use module_pp
      use module_levels
      use module_utils
      implicit none
      include "mpif.h"

      integer,parameter :: kr = kind(1.D0)

      ! parameters to set ######################
      ! maximal length of problemname
      integer,parameter:: lproblemnamex = 100
      ! maximal length of any used file - should be reasonably larger than length of problem to allow suffices
      integer,parameter:: lfilenamex = 130
      ! use prepared division into subdomains in file *.ES?
      logical,parameter :: load_division = .false.
      ! type of matrix (0 - general, 1 - SPD, 2 - general symmetric)
      integer,parameter :: matrixtype = 1 ! SPD matrix
      ! should arithmetic averages be used?
      logical,parameter :: use_arithmetic = .true.
      ! should adaptive averages be used?
      logical,parameter :: use_adaptive   = .false.
      ! should parallel division be used (e.g. ParMETIS instead of METIS)?
      logical,parameter :: parallel_division = .true.
      ! how many nodes has to be shared by two elements to consider them adjacent?
      integer,parameter :: neighbouring = 4
      ! what is the base for array indexing ( 1 in Fortran )
      integer,parameter :: numbase = 1
      ! end of parameters to set ######################


      ! number of levels
      integer :: nlevels
      ! subdomains in levels
      integer ::            lnsublev
      integer,allocatable :: nsublev(:)

      integer :: nsub_loc_1

      !  parallel variables
      integer :: myid, comm_all, comm_self, nproc, ierr
      integer :: idpar, idgmi, idml, idfvs, idrhs

      integer :: ndim, nsub, nelem, ndof, nnod, linet
      integer :: lproblemname
      integer :: maxit, ndecrmax, meshdim
      real(kr):: tol

      integer :: idelm, ios

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

      character(lproblemnamex) :: problemname 
      character(lfilenamex) :: name
      character(lfilenamex) :: filename

      ! MPI initialization
!***************************************************************PARALLEL
      call MPI_INIT(ierr)
      ! Communicator
      comm_all  = MPI_COMM_WORLD
      comm_self = MPI_COMM_SELF
      call MPI_COMM_RANK(comm_all,myid,ierr)
      call MPI_COMM_SIZE(comm_all,nproc,ierr)
!***************************************************************PARALLEL

      write (*,*) 'I am processor ',myid,': Hello nproc =',nproc


! Initial screen
      if (myid.eq.0) then
         write(*,'(a)') 'MULTILEVEL TESTER'
         write(*,'(a)') '================='
      end if
      ! name of the problem
      call pp_pget_problem_name(comm_all,problemname,lproblemnamex,lproblemname)

      if (myid.eq.0) then
         name = problemname(1:lproblemname)//'.PAR'
         call allocate_unit(idpar)
         open (unit=idpar,file=name,status='old',form='formatted')
      end if
      call pp_pread_par_file(comm_all,idpar, ndim, nsub, nelem, ndof, nnod, linet, tol, maxit, ndecrmax, meshdim)
      if (myid.eq.0) then
         close (idpar)
      end if


      write (*,*) 'myid = ',myid,': Initializing LEVELS.'
       ! read PMD mesh and input it as zero level coarse problem
      lnnet = nelem
      lnndf = nnod
      lxyz1 = nnod
      lxyz2 = ndim
      allocate(inet(linet),nnet(lnnet),nndf(lnndf),xyz(lxyz1,lxyz2))
      if (myid.eq.0) then
         filename = problemname(1:lproblemname)//'.GMIS'
         call allocate_unit(idgmi)
         open (unit=idgmi,file=filename,status='old',form='formatted')
      end if
      call pp_pread_pmd_mesh(comm_all,idgmi,inet,linet,nnet,lnnet,nndf,lnndf,xyz,lxyz1,lxyz2)
      if (myid.eq.0) then
         close(idgmi)
      end if

      ! read PMD boundary conditions 
      lifix = ndof
      lfixv = ndof
      allocate(ifix(lifix),fixv(lfixv))
      if (myid.eq.0) then
         filename = problemname(1:lproblemname)//'.FVS'
         call allocate_unit(idfvs)
         open (unit=idfvs,file=filename,status='old',form='formatted')
      end if
      call pp_pread_pmd_bc(comm_all,idfvs,ifix,lifix,fixv,lfixv)
      if (myid.eq.0) then
         close(idfvs)
      end if

      ! read PMD right-hand side
      lrhs = ndof
      allocate(rhs(lrhs))
      if (myid.eq.0) then
         filename = problemname(1:lproblemname)//'.RHS'
         call allocate_unit(idrhs)
         open (unit=idrhs,file=filename,status='old',form='unformatted')
      end if
      call pp_pread_pmd_rhs(comm_all,idrhs,rhs,lrhs)
      if (myid.eq.0) then
         close(idrhs)
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
         write(*,*)'  number of subdomains in levels: ',nsublev
         call flush(6)
      end if

      ! prepare initial solution
      lsol = ndof
      allocate(sol(lsol))
      sol = 0._kr

      ! attach file with element matrices
      if (myid.eq.0) then
         ! ELM - element stiffness matrices - structure:
         call allocate_unit(idelm)
         open(unit=idelm,file=trim(problemname)//'.ELM',status='old',form='unformatted',iostat=ios)
         if (ios.ne.0) then
            call error('TEST_MODULE_LEVELS','Problem opening file '//trim(problemname)//'.ELM')
         end if
         rewind idelm
      end if


      write (*,*) 'myid = ',myid,': Initializing LEVELS.'
      nsub_loc_1 = -1
      call levels_init(nlevels,nsublev,lnsublev,nsub_loc_1,comm_all,1,numbase)
      deallocate(nsublev)
      write (*,*) 'myid = ',myid,': Number of assigned subdomains:',nsub_loc_1
      call levels_upload_global_data(nelem,nnod,ndof,ndim,meshdim,&
                                     inet,linet,nnet,lnnet,nndf,lnndf,xyz,lxyz1,lxyz2,&
                                     ifix,lifix,fixv,lfixv,rhs,lrhs,sol,lsol,idelm, neighbouring, load_division)
      deallocate(inet,nnet,nndf,xyz)
      deallocate(ifix,fixv)
      deallocate(rhs)
      deallocate(sol)

!      call levels_init(nlevels,nsub)

      ! create first two levels
      ! open file with description of pairs
!      ilevel = 1
!      call levels_read_level_from_file(problemname,myid,comm_all,ndim,ilevel)
!
!      ilevel = 2
!      call levels_read_level_from_file(problemname,myid,comm_all,ndim,ilevel)
!
!      ! associate subdomains with first level
!      ilevel = 1
!      call levels_prepare_standard_level(ilevel,nsub,1,nsub)
!
!      call levels_prepare_last_level(myid,nproc,comm_all,comm_self,matrixtype,ndim,problemname)
      call levels_pc_setup(parallel_division, &
                           matrixtype,use_arithmetic, use_adaptive)

   !   lvec = 15
   !   allocate(vec(lvec))
   !   vec = 1.0_kr
   !   call levels_pc_apply(vec,lvec, )
   !   if (myid.eq.0) then
   !      write(*,*) 'vec:'
   !      write(*,'(e18.9)') vec
   !   end if
   !   deallocate(vec)

      write (*,*) 'myid = ',myid,': Finalize LEVELS.'
      call levels_finalize

      ! MPI finalization
!***************************************************************PARALLEL
      call MPI_FINALIZE(ierr)
!***************************************************************PARALLEL

end program
