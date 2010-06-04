!***********************************************************************
      program bddcml
!***********************************************************************
! module for distributed Krylov data storage
      use module_krylov_types_def
! module for preconditioner
      use module_levels
! module for domain decomposition data
      use module_dd
! Program name
      use module_utils

      implicit none
      
      include "mpif.h"

      integer,parameter :: kr = kind(1.D0)

      !  parallel variables
      integer :: myid, comm_all, comm_self, nproc, ierr
      integer :: idpar

      ! number of levels
      integer :: nlevels

      integer :: ndim, nsub, nelem, ndof, nnod, nnodc, linet
      integer :: maxit, ndecrmax
      real(kr) :: tol

      integer :: isub

      integer :: matrixtype

      integer :: nnods, nelems, ndofs, ndofis, nnodis
      integer :: lsolis, lresis

      integer :: lsols
      real(kr),allocatable :: sols(:)


      character(90)  :: problemname 
      character(100) :: filename

      integer ::                          lpcg_data
      type (pcg_data_type), allocatable :: pcg_data(:)


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
         write(*,'(a)') 'MULTILEVEL BDDC solver'
         write(*,'(a)') '======================'

! Name of the problem
   10    write(*,'(a,$)') 'Name of the problem: '
         read(*,*) problemname
         if(problemname.eq.' ') goto 10

      end if
! Broadcast of name of the problem      
!***************************************************************PARALLEL
      call MPI_BCAST(problemname, 90, MPI_CHARACTER, 0, comm_all, ierr)
!***************************************************************PARALLEL

      if (myid.eq.0) then
         filename = trim(problemname)//'.PAR'
         call allocate_unit(idpar)
         open (unit=idpar,file=filename,status='old',form='formatted')
      end if

! Reading basic properties 
      if (myid.eq.0) then
         read(idpar,*) ndim, nsub, nelem, ndof, nnod, nnodc, linet,     &
                       tol, maxit, ndecrmax
         write(*,*)'Characteristics of the problem ',trim(problemname), ':'
         write(*,*)'  number of processors            nproc =',nproc
         write(*,*)'  number of dimensions             ndim =',ndim
         write(*,*)'  number of subdomains             nsub =',nsub
         write(*,*)'  number of elements global       nelem =',nelem
         write(*,*)'  number of DOF                    ndof =',ndof
         write(*,*)'  number of nodes global           nnod =',nnod
         write(*,*)'  number of constrained nodes     nnodc =',nnodc
         write(*,*)'  lenght of field INET            linet =',linet
         write(*,*)'Characteristics of iterational process:'
         write(*,*)'  tolerance of error                tol =',tol
         write(*,*)'  maximum number of iterations    maxit =',maxit
         write(*,*)'  number of incresing residual ndecrmax =',ndecrmax
         call flush(6)
      end if
! Broadcast basic properties of the problem
!***************************************************************PARALLEL
      call MPI_BCAST(ndim,     1,MPI_INTEGER,         0, comm_all, ierr)
      call MPI_BCAST(nsub,     1,MPI_INTEGER,         0, comm_all, ierr)
      call MPI_BCAST(nelem,    1,MPI_INTEGER,         0, comm_all, ierr)
      call MPI_BCAST(ndof,     1,MPI_INTEGER,         0, comm_all, ierr)
      call MPI_BCAST(nnod,     1,MPI_INTEGER,         0, comm_all, ierr)
      call MPI_BCAST(nnodc,    1,MPI_INTEGER,         0, comm_all, ierr)
      call MPI_BCAST(linet,    1,MPI_INTEGER,         0, comm_all, ierr)
      call MPI_BCAST(tol,      1,MPI_DOUBLE_PRECISION,0, comm_all, ierr)
      call MPI_BCAST(maxit,    1,MPI_INTEGER,         0, comm_all, ierr)
      call MPI_BCAST(ndecrmax, 1,MPI_INTEGER,         0, comm_all, ierr)
!***************************************************************PARALLEL

      write (*,*) 'myid = ',myid,': Initializing LEVELS.'
      !============
      ! number of levels
      nlevels = 2
      !============
      call levels_init(nlevels,nsub)

      matrixtype = 1 ! SPD matrix
      call levels_pc_setup(problemname,myid,nproc,comm_all,comm_self,matrixtype,ndim,nsub)


! Input zeros as initial vector of solution
      lpcg_data = nsub
      allocate(pcg_data(lpcg_data))
      do isub = 1,nsub

         ! deterine size of subdomain
         call dd_get_size(myid,isub,ndofs,nnods,nelems)
         if (ndofs.ge.0) then
            pcg_data(isub)%is_mine = .true.
         else
            pcg_data(isub)%is_mine = .false.
            cycle
         end if
         call dd_get_interface_size(myid,isub,ndofis,nnodis)

         ! allocate local subdomain solution
         lsols  = ndofs
         allocate(sols(lsols))

         ! set zero initial guess
         ! u_0 = 0
         call zero(sols,lsols)

         ! set initial guess to satisfy Dirichlet boundary conditions
         call dd_fix_bc(myid,isub, sols,lsols)

         ! allocate vectors for Krylov method 
         lsolis = ndofis
         lresis = ndofis
         pcg_data(isub)%lsoli = lsolis
         allocate(pcg_data(isub)%soli(lsolis))
         call zero(pcg_data(isub)%soli,lsolis)
         pcg_data(isub)%lresi = lresis
         allocate(pcg_data(isub)%resi(lresis))
         call zero(pcg_data(isub)%resi,lresis)

         ! restrict solution to interface
         call dd_map_sub_to_subi(myid,isub, sols,lsols, pcg_data(isub)%soli,lsolis)
         deallocate(sols)

         ! set initial residual to RHS
         ! res = g
         call dd_get_reduced_rhs(myid,isub, pcg_data(isub)%resi,lresis)
         write(*,*) 'myid =',myid, 'g = ',pcg_data(isub)%resi

      end do

      call bddcpcg(pcg_data,lpcg_data, nsub, myid,comm_all,tol,maxit,ndecrmax)




      ! Clear memory
      do isub = 1,nsub
         if (pcg_data(isub)%is_mine) then
            deallocate(pcg_data(isub)%soli)
            deallocate(pcg_data(isub)%resi)
         end if
      end do
      deallocate(pcg_data)

      write (*,*) 'myid = ',myid,': Finalize LEVELS.'
      call levels_finalize

      ! MPI finalization
!***************************************************************PARALLEL
      call MPI_FINALIZE(ierr)
!***************************************************************PARALLEL
      end program

!***********************************************************************************
      subroutine bddcpcg(pcg_data,lpcg_data, nsub, myid,comm_all,tol,maxit,ndecrmax)
!***********************************************************************************
! subroutine realizing PCG algorithm with vectors distributed by subdomains

! module for distributed Krylov data storage
      use module_krylov_types_def
! module for preconditioner
      use module_levels
! module for domain decomposition data
      use module_dd
! Program name
      use module_utils

      implicit none
      
      include "mpif.h"

      integer,parameter :: kr = kind(1.D0)

      integer,intent(in) ::                 lpcg_data
      type (pcg_data_type), intent(inout) :: pcg_data(lpcg_data)

      ! number of subdomains
      integer,intent(in) :: nsub

      ! parallel variables
      integer,intent(in) :: myid, comm_all 

      ! limit on iterations
      integer,intent(in) :: maxit

      ! limit on iterations with increasing residual
      integer,intent(in) :: ndecrmax

      ! desired accuracy of relative residual
      real(kr),intent(in) :: tol

      ! local vars
      integer :: ndofis, nnodis
      integer :: lp, lap, lresi, lsoli, lz
      integer :: isub, i

      ! PCG vars
      real(kr) :: normrhs, normres2, normres
      real(kr) :: normres2_loc
      real(kr) :: normres2_sub

      ! MPI vars
      integer :: ierr

      ! prepare memory for PCG
      do isub = 1,nsub
         if (pcg_data(isub)%is_mine) then

            call dd_get_interface_size(myid,isub,ndofis,nnodis)

            lap = ndofis
            pcg_data(isub)%lap = lap
            allocate(pcg_data(isub)%ap(lap))
            allocate(pcg_data(isub)%apadj(lap))

            lp  = ndofis
            pcg_data(isub)%lp = lp
            allocate(pcg_data(isub)%p(lp))
            allocate(pcg_data(isub)%padj(lp))

            lresi = pcg_data(isub)%lresi
            allocate(pcg_data(isub)%resiadj(lresi))

            lz = ndofis
            pcg_data(isub)%lz = lz
            allocate(pcg_data(isub)%z(lz))
            allocate(pcg_data(isub)%zadj(lz))

         end if
      end do


      ! get initial residual
      ! r_0 = g - A*u_0
      ! ap = A*u_0
      do isub = 1,nsub
         if (pcg_data(isub)%is_mine) then

            call zero(pcg_data(isub)%ap,pcg_data(isub)%lap)

            call dd_multiply_by_schur(myid,isub,&
                                      pcg_data(isub)%soli,pcg_data(isub)%lsoli, &
                                      pcg_data(isub)%ap,pcg_data(isub)%lap)

            call dd_comm_upload(myid, isub,  pcg_data(isub)%ap,pcg_data(isub)%lap)
         end if
      end do

      ! Interchange data
      call dd_comm_swapdata(myid, nsub, comm_all)
      ! Download data
      do isub = 1,nsub
         if (pcg_data(isub)%is_mine) then

            call zero(pcg_data(isub)%apadj,pcg_data(isub)%lap)
            ! get contibution from neigbours
            call dd_comm_download(myid, isub,  pcg_data(isub)%apadj,pcg_data(isub)%lap)
            ! join data
            do i = 1,pcg_data(isub)%lap
               pcg_data(isub)%ap(i) = pcg_data(isub)%ap(i) + pcg_data(isub)%apadj(i)
            end do

            ! update residual
            ! r_0 = g - A*u_0
            do i = 1,pcg_data(isub)%lresi
               pcg_data(isub)%resi(i) = pcg_data(isub)%resi(i) - pcg_data(isub)%ap(i)
            end do
         end if
      end do

      ! compute norm of right hand side
      normres2_loc = 0._kr
      do isub = 1,nsub
         if (pcg_data(isub)%is_mine) then

            call dd_dotprod_local(myid,isub,pcg_data(isub)%resi,pcg_data(isub)%lresi, &
                                            pcg_data(isub)%resi,pcg_data(isub)%lresi, &
                                            normres2_sub)

            normres2_loc = normres2_loc + normres2_sub
         end if
      end do
!***************************************************************PARALLEL
      call MPI_ALLREDUCE(normres2_loc,normres2, 1, MPI_DOUBLE_PRECISION,&
                         MPI_SUM, comm_all, ierr) 
!***************************************************************PARALLEL
      normrhs = sqrt(normres2)
      if (myid.eq.1) then
         write(*,*) 'Norm of the right hand side: ',normrhs
      end if

      ! Check of zero right hand side => all zero solution
      if (normrhs.eq.0.0D0) then
         if (myid.eq.0) then
            write(*,*) 'initial residual zero => initial solution exact'
         end if
         return 
      end if


! Initial action of the preconditioner M on residual vector RESI
! M*resi => p
      if (myid.eq.0) then
         write(*,*) ' Initial action of preconditioner'
      end if
      call levels_pc_apply(myid,comm_all, pcg_data,lpcg_data)
      ! produced new z

     
 


      ! Clear memory of PCG
      do isub = 1,nsub
         if (pcg_data(isub)%is_mine) then
            deallocate(pcg_data(isub)%ap)
            deallocate(pcg_data(isub)%apadj)
            deallocate(pcg_data(isub)%p)
            deallocate(pcg_data(isub)%padj)
            deallocate(pcg_data(isub)%resiadj)
            deallocate(pcg_data(isub)%z)
            deallocate(pcg_data(isub)%zadj)
         end if
      end do

      end subroutine

! Input zeros as initial vector of solution
!      call rzero(soli,lsoli)
!
!C put Dirichlet BC into apropriate positions of initial solution
!
!      do 1005 isub_loc = 1,nsub_loc
!         isub = start + isub_loc - 1
!
!         jnndfs    = (isub_loc-1)*nnodsx
!         jiinsn    = (isub_loc-1)*nnodisx
!         jifixns   = (isub_loc-1)*lsolsx
!         jsoli     = (isub_loc-1)*lsolisx
!         ji        = (isub_loc-1)*lsolisx
!
!C Creation of field kdofs
!         kdofs(1) = 0
!         nnods = nnodsa(isub)
!         do in = 2,nnods
!            kdofs(in) = kdofs(in-1) + nndfs(jnndfs + in-1)
!         end do
!
!C Prepare initial solution
!         call getfname(name1,isub,'FVSS',fname)
!         open(unit=idfvss,file=fname,status='old',
!     *        form='unformatted')
!         rewind idfvss
!C GFIXV is vector of fixed variables - here only first column is used
!         call rzero(gfixv,lgfixv)
!         lsols  = lsolsa(isub)
!         do i = 1,lsols
!            read(idfvss) gfixv(i)
!         end do
!         close(idfvss)
!
!         nnodis = nnodisa(isub)
!         indvi = 0
!         do i = 1,nnodis
!            indns = iinsn(jiinsn + i)
!            ndofn = nndfs(jnndfs + indns)
!            point1 = kdofs(indns)
!            do idofn = 1,ndofn
!               indvs = point1 + idofn
!               indvi = indvi + 1
!               if(ifixns(jifixns + indvs).gt.0) then
!                  soli(jsoli + indvi) = gfixv(indvs)
!               end if
!            end do
!         end do
!
!C Find staticly condensed RHS : g = rhsi - Kio*(Koo^-1)*rhso
!         call getfname(name1,isub,'RHSS',fname)
!         open(unit=idrhss,file=fname,status='old',
!     *        form='unformatted')
!         rewind idrhss
!C Read ASLOD - subdomain right hand side - not weighted!
!         call rzero(aslod,laslod)
!         lsols = lsolsa(isub)
!         do i = 1,lsols
!            read(idrhss) aslod(i)
!         end do
!         close(idrhss)
      
!C Weight the rhs in interface nodes

!         nnodis = nnodisa(isub)
!         indi = 0
!         do ii = 1,nnodis
!            ins = iinsn(jiinsn + ii)
!            ndofn = nndfs(jnndfs + ins)
!            ind  = kdofs(ins)
!            do idofn = 1,ndofn
!               ind = ind + 1
!               indi = indi + 1
!               aslod(ind) = wi(ji + indi)*aslod(ind)
!            end do
!         end do
!
!         call getfname(name1,isub,'EQ2',fname)
!         open(unit=ideq2, file=fname, status='old',
!     *        form='unformatted', access='direct', recl=lrec)
!         nrhs = 1
!         nreq2 = nreq2a(isub)
!         lsols  = lsolsa(isub)
!         call frors(isub_loc, mfronx, ideq2,
!     *              lsolsx, lsols, nrhs,       
!     *              ibuf,libuf, rbuf,lrbuf, gload,lgload,
!     *              aslod,laslod, rea,lrea, ifixins,lifixins, 
!     *              gfixv,lgfixv,nreq2)
!         close(ideq2)
!
!         nnodis = nnodisa(isub)
!         indi = 0
!         do ii = 1,nnodis
!            ins = iinsn(jiinsn + ii)
!            ndofn = nndfs(jnndfs + ins)
!            ind = kdofs(ins)
!            do idofn = 1,ndofn
!               ind = ind + 1
!               indi = indi + 1
!               resaux(ji + indi) = rea(ind)
!            end do
!         end do
!
! 1005 end do
!
!C Introducing vector of local residual RESAUX
!      do 1010 isub_loc = 1,nsub_loc
!         isub = start + isub_loc - 1
!
!         ji     = (isub_loc-1)*lsolisx
!         jiinsn = (isub_loc-1)*nnodisx
!         jnndfs = (isub_loc-1)*nnodsx
!         jifixns = (isub_loc-1)*lsolsx
!
!C Creation of field kdofs
!         kdofs(1) = 0
!         nnods = nnodsa(isub)
!         do in = 2,nnods
!            kdofs(in) = kdofs(in-1) + nndfs(jnndfs + in-1)
!         end do
!C In the vector of residual, natural fixed variables are zeros
!         nnodis = nnodisa(isub)
!         indi = 0
!         do i = 1,nnodis
!            indns = iinsn(jiinsn + i)
!            inds = kdofs(indns)
!            ndofn = nndfs(jnndfs + indns)
!            do idofn = 1,ndofn
!               indi = indi + 1
!               inds = inds + 1
!               if (ifixns(jifixns + inds).gt.0) then
!                  resaux(ji + indi) = 0.0D0
!               end if
!            end do
!         end do
! 1010 end do
!
!C Create vector of global residual RES
!C***************************************************************PARALLEL
!      call swapdata(myid, nsub_loc, nsub_locx, start, finish, comm,
!     *              nsub, nadjsx, nnodsx, nnodisx, lsolisx, 
!     *              nsnodsx, nsdofsx,
!     *              nadjsa,lnadjsa, nnodisa,lnnodisa, 
!     *              lsolisa,llsolisa, nsnodsa,lnsnodsa,
!     *              iadjs,liadjs, nsdofadj,lnsdofadj,
!     *              ishnin,lishnin, iinsn,liinsn,
!     *              commveco,commveci,lcommvec, nndfs,lnndfs,
!     *              kdofi,lkdofi,
!     *              res,resaux,lres)
!C***************************************************************PARALLEL
!
!C Scalar product of RES and RESAUX to determine global norm of RES
!      normres2_loc = 0.0D0
!      do 1020 isub_loc = 1,nsub_loc
!         isub = start + isub_loc - 1
!         
!         call scalpvv(lsolisa(isub), lsolisx, isub_loc, 
!     *                res,lres, resaux,lres, normres2_loc)
! 1020 end do
!C***************************************************************PARALLEL
!      call MPI_ALLREDUCE(normres2_loc,normres2, 1, MPI_DOUBLE_PRECISION,
!     *                   MPI_SUM, comm, ierr) 
!C***************************************************************PARALLEL
!      normrhs = sqrt(normres2)
!
!C Check of zero right hand side => all zero solution
!      if (normrhs.eq.0.0D0) then
!         if (myid.eq.0) then
!            write(*,*) 'all zero RHS => all zero solution'
!         end if
!         goto 58
!      end if
!      ! debug
!!      print *, 'norm rhs initial *************************',normrhs
!
!C Initial action of the preconditioner M on residual vector RES
!C M*res => p
!C Collect residual into vector SOLC
!      if (myid.eq.0) then
!         write(*,*) ' Initial action of preconditioner'
!         write(*,*) '  collect residual'
!      end if
!      call rzero(solcaux,lsolc)
!      do 1030 isub_loc = 1,nsub_loc
!         isub = start + isub_loc - 1
!        
!         lsolis = lsolisa(isub)
!         nnodcs = nnodcsa(isub)
!         nglbs  = nglbsa(isub)
!         nrhss  = nrhssa(isub)
!         call colres(isub_loc, nnodc,
!     *               lphisisx, nnodsx, lsolisx, nnodcsx,
!     *               nglbsx, linglbsx, 
!     *               lsolis, nnodcs, nglbs, nrhss,
!     *               icngcn,licngcn, icnsn,licnsn, nndfs,lnndfs,
!     *               isglbgglbn,lisglbgglbn, inglbs,linglbs,
!     *               nnglbs,lnnglbs,
!     *               kdofc,lkdofc,
!     *               phisi,lphisi, wi,lwi,
!     *               reshelp1,lreshelp1, reshelp2,lreshelp2,
!     *               res,lres, solcaux,lsolc)
! 1030 end do
!C***************************************************************PARALLEL
!      call MPI_REDUCE(solcaux,solc,lsolc, MPI_DOUBLE_PRECISION, MPI_SUM,
!     *                                    0, comm, ierr) 
!C***************************************************************PARALLEL
!
!C******************************************************************MUMPS
!C Solution of the coarse problem KCG*UC = RESC (here stored in UC)
!      ! debug
!      !if (myid.eq.0) then
!      !   print *,'Ini norm of c. residual',sqrt(dot_product(solc,solc))
!      !end if
!      call mumps_resolve(mumps_struc,solc,lsolc)
!      ! debug
!      !if (myid.eq.0) then
!      !   print *,'Ini norm of c. solution',sqrt(dot_product(solc,solc))
!      !end if
!C******************************************************************MUMPS
!
!C Substructure correction into PAUX
!      if (myid.eq.0) then
!         write(*,*) '  substructure correction'
!      end if
!      call rzero(paux,lp)
!      do 1040 isub_loc = 1,nsub_loc
!         isub = start + isub_loc - 1
!
!         nnods  = nnodsa(isub)
!         lsols  = lsolsa(isub)
!         nnodis = nnodisa(isub)
!         nglbv  = nrhssa(isub) - lsolcsa(isub)
!         nreq1  = nreq1a(isub)
!         call cors(isub_loc, start, ideq1, lrec, name1,lname1,
!     *             nnodsx, lsolsx, nnodisx, lsolisx, lcfsx, nrhssx,
!     *             nnods, lsols, nnodis, nglbv,
!     *             nreq1, mfronx,
!     *             iinsn,liinsn,
!     *             nndfs,lnndfs, kdofs,lkdofs, ipiv,lipiv,
!     *             wi,lwi, cf,lcf, dualm,ldualm, scaling,lscaling,
!     *             dualrhs,ldualrhs, ifixcns,lifixcns, ibuf,libuf, 
!     *             gfixv,lgfixv, gload,lgload, rbuf,lrbuf, rea,lrea,
!     *             aslod,laslod, res,lres, paux,lp)
! 1040 end do
!
!C***************************************************************PARALLEL
!      call MPI_BCAST(solc, lsolc, MPI_DOUBLE_PRECISION, 0, comm, ierr)
!C***************************************************************PARALLEL
!
!C Coarse grid correction and substructure correction into PAUX
!      if (myid.eq.0) then
!         write(*,*) '  coarse grid correction'
!      end if
!      do 1045 isub_loc = 1,nsub_loc
!         isub = start + isub_loc - 1
!
!         lsolis = lsolisa(isub)
!         nnodcs = nnodcsa(isub)
!         nglbs  = nglbsa(isub)
!         nrhss  = nrhssa(isub)
!         call corc(isub_loc, start, nnodc,
!     *             nnodsx, lsolisx, nnodcsx,lphisisx, nglbsx,linglbsx,
!     *             lsolis, nnodcs, nglbs, nrhss, 
!     *             icngcn,licngcn, icnsn,licnsn,
!     *             nndfs,lnndfs, 
!     *             isglbgglbn,lisglbgglbn, inglbs,linglbs,
!     *             nnglbs,lnnglbs,
!     *             kdofc,lkdofc,
!     *             phisi,lphisi, wi,lwi,
!     *             solhelp1,lsolhelp1, solhelp2,lsolhelp2,
!     *             solc,lsolc, paux,lp)
! 1045 end do
!
!C Construction of global vector of search direction P
!C***************************************************************PARALLEL
!      call swapdata(myid, nsub_loc, nsub_locx, start, finish, comm,
!     *              nsub, nadjsx, nnodsx, nnodisx, lsolisx, 
!     *              nsnodsx, nsdofsx,
!     *              nadjsa,lnadjsa, nnodisa,lnnodisa, 
!     *              lsolisa,llsolisa, nsnodsa,lnsnodsa,
!     *              iadjs,liadjs, nsdofadj,lnsdofadj,
!     *              ishnin,lishnin, iinsn,liinsn,
!     *              commveco,commveci,lcommvec, nndfs,lnndfs,
!     *              kdofi,lkdofi,
!     *              p,paux,lp)
!C***************************************************************PARALLEL
!
!C Scalar product of vectors of residual and direction - res*p => rmp
!      rmp_loc = 0.0D0
!      do 1050 isub_loc = 1,nsub_loc
!         isub = start + isub_loc - 1
!         
!         call scalpvv(lsolisa(isub), lsolisx, isub_loc, 
!     *                res,lres, paux,lp, rmp_loc)
! 1050 end do
!C***************************************************************PARALLEL
!      call MPI_ALLREDUCE(rmp_loc,rmp, 1, MPI_DOUBLE_PRECISION,
!     *                                MPI_SUM, comm, ierr) 
!C***************************************************************PARALLEL
!
!C Control of positive definiteness of preconditioner matrix
!      if (rmp.le.0.0D0) then
!         if (myid.eq.0) then
!            write(*,*) 'Preconditioner not positive definite!'
!            call MPI_ABORT(comm, 78, ierr)
!         end if
!      end if
!      ! debug
!!      print *, 'rmp initial *************************',rmp
!
!C Setting up the properties for decreasing residual
!      ndecr   = 0
!      lastres = 1.0D0
!      
!
!C***********************************************************************
!C*************************MAIN LOOP OVER ITERATIONS*********************
!C***********************************************************************
!
!      do 1060 iter = 1,maxit
!
!C Multiplication of P by local system matrix - A*p => ap
!         if (myid.eq.0) then
!            write(*,*) ' Multiplication by system matrix'
!         end if
!         call rzero(apaux,lap)
!         do 2000 isub_loc = 1,nsub_loc
!            isub = start + isub_loc - 1
!
!            nnods  = nnodsa(isub)
!            lsols  = lsolsa(isub)
!            nnodis = nnodisa(isub)
!            nreq2  = nreq2a(isub)
!            call smmult(isub_loc, start, ideq2,lrec,name1,lname1,
!     *                  nnodsx, lsolsx, nnodisx, lsolisx,
!     *                  nnods, lsols, nnodis, nreq2, mfronx,
!     *                  iinsn,liinsn, nndfs,lnndfs,
!     *                  ifixins,lifixins, kdofs,lkdofs, ibuf,libuf,
!     *                  gfixv,lgfixv, gload,lgload, 
!     *                  rbuf,lrbuf, rea,lrea, aslod,laslod, 
!     *                  p,lp, apaux,lap)
! 2000    end do
!
!
!C Construction of global vector AP
!C***************************************************************PARALLEL
!      call swapdata(myid, nsub_loc, nsub_locx, start, finish, comm,
!     *              nsub, nadjsx, nnodsx, nnodisx, lsolisx, 
!     *              nsnodsx, nsdofsx,
!     *              nadjsa,lnadjsa, nnodisa,lnnodisa, 
!     *              lsolisa,llsolisa, nsnodsa,lnsnodsa,
!     *              iadjs,liadjs, nsdofadj,lnsdofadj,
!     *              ishnin,lishnin, iinsn,liinsn,
!     *              commveco,commveci,lcommvec, nndfs,lnndfs,
!     *              kdofi,lkdofi,
!     *              ap,apaux,lap)
!C***************************************************************PARALLEL
!
!C Scalar product of vectors of old search direction and ap - p*ap => pap
!         pap_loc = 0.0D0
!         do 2010 isub_loc = 1,nsub_loc
!            isub = start + isub_loc - 1
!            
!            call scalpvv(lsolisa(isub), lsolisx, isub_loc, 
!     *                   p,lp, apaux,lap, pap_loc)
! 2010    end do
!C***************************************************************PARALLEL
!         call MPI_ALLREDUCE(pap_loc,pap, 1, MPI_DOUBLE_PRECISION,
!     *                      MPI_SUM, comm, ierr) 
!C***************************************************************PARALLEL
!
!C Control of positive definitenes of system matrix
!         if (pap.le.0.0D0) then
!            if (myid.eq.0) then
!               write(*,*) ' System matrix not positive definite!'
!               call MPI_ABORT(comm, 78, ierr)
!            end if
!         end if
!         ! debug
!         !print *, 'iter',iter,'pap  *************************',pap
!
!C Determination of step lenght ALPHA
!         alpha = rmp/pap
!
!C Correction of solution vector SOLI and residual vector RES
!         do 2020 isub_loc = 1,nsub_loc
!            isub = start + isub_loc - 1
!            
!            ji      = (isub_loc-1)*lsolisx
!            jiinsn  = (isub_loc-1)*nnodisx
!            jnndfs  = (isub_loc-1)*nnodsx
!            jifixns = (isub_loc-1)*lsolsx
!
!C Creation of field kdofs
!            kdofs(1) = 0
!            nnods = nnodsa(isub)
!            do in = 2,nnods
!               kdofs(in) = kdofs(in-1) + nndfs(jnndfs + in-1)
!            end do
!
!C In the vector of residual, natural fixed variables are zeros
!            nnodis = nnodisa(isub)
!            indi = 0
!            do i = 1,nnodis
!               indns = iinsn(jiinsn + i)
!               inds = kdofs(indns)
!               ndofn = nndfs(jnndfs + indns)
!               do idofn = 1,ndofn
!                  indi = indi + 1
!                  inds = inds + 1
!                  if ((ifixns(jifixns + inds).gt.0).and.
!     *                (p(ji + indi).ne.0.0D0)) then
!                     write(*,*) 'Correction of solution not consistent'
!                  end if
!
!                  soli(ji + indi) = soli(ji + indi) + alpha*p(ji + indi)
!                  res(ji + indi) = res(ji + indi) - alpha*ap(ji + indi)
!C Residual is zero in fixed variables
!
!                  if (ifixns(jifixns + inds).gt.0) then
!                     res(ji + indi) = 0.0D0
!                  end if
!                  resaux(ji + indi) = wi(ji + indi)*res(ji + indi)
!               end do
!            end do
! 2020    end do
!
!C Scalar product of RES and RESAUX to determine global norm of RES
!         normres2_loc = 0.0D0
!         do 2030 isub_loc = 1,nsub_loc
!            isub = start + isub_loc - 1
!            
!            call scalpvv(lsolisa(isub), lsolisx, isub_loc, 
!     *                   res,lres, resaux,lres, normres2_loc)
! 2030    end do
!C***************************************************************PARALLEL
!         call MPI_ALLREDUCE(normres2_loc,normres2, 1, 
!     *                      MPI_DOUBLE_PRECISION, MPI_SUM, comm, ierr)
!C***************************************************************PARALLEL
!         normres = sqrt(normres2)
!         ! debug
!!         print *, 'iter',iter,'normres *************************',
!!     *            normres
!         ! debug
!!         print *, 'iter',iter,'normrhs *************************',
!!     *            normrhs
!
!C Evaluation of stop criterion
!         relres = normres/normrhs
!            
!         if (myid.eq.0) then
!            write(* ,5001) iter, relres
! 5001       format(1X,'iteration iter = ',I4,2X,'relres = ',F25.18)
!         end if
!
!         if (relres.lt.tol) then
!            nw = iter-1
!            if (myid.eq.0) then
!               write(*,*)'Number of PCG iterations:',iter
!            end if
!            goto 58
!         end if
!
!C Check of decreasing of residual
!         if (relres.lt.lastres) then
!            ndecr = 0
!         else
!            ndecr = ndecr + 1
!            if (ndecr.ge.ndecrmax) then
!               if (myid.eq.0) then
!                  write(*,*)'Residual did not decrease for',ndecrmax,
!     *                      ' iterations'
!                  call MPI_ABORT(comm, 78, ierr)
!               end if
!            end if
!         end if
!         lastres = relres
!
!C Action of the preconditioner M on residual vector RES - M*res => h
!C Collect residual into vector SOLC
!         if (myid.eq.0) then
!            write(*,*) ' Action of preconditioner'
!            write(*,*) '  collect residual'
!         end if
!         call rzero(solcaux,lsolc)
!         do 2040 isub_loc = 1,nsub_loc
!            isub = start + isub_loc - 1
!            
!            lsolis = lsolisa(isub)
!            nnodcs = nnodcsa(isub)
!            nglbs  = nglbsa(isub)
!            nrhss  = nrhssa(isub)
!            call colres(isub_loc, nnodc,
!     *                  lphisisx, nnodsx, lsolisx, nnodcsx,
!     *                  nglbsx, linglbsx, 
!     *                  lsolis, nnodcs, nglbs, nrhss,
!     *                  icngcn,licngcn, icnsn,licnsn, nndfs,lnndfs,
!     *                  isglbgglbn,lisglbgglbn, inglbs,linglbs,
!     *                  nnglbs,lnnglbs,
!     *                  kdofc,lkdofc,
!     *                  phisi,lphisi, wi,lwi,
!     *                  reshelp1,lreshelp1, reshelp2,lreshelp2,
!     *                  res,lres, solcaux,lsolc)
! 2040    end do
!C***************************************************************PARALLEL
!         call MPI_REDUCE(solcaux,solc,lsolc, MPI_DOUBLE_PRECISION, 
!     *                                       MPI_SUM, 0, comm, ierr) 
!C***************************************************************PARALLEL
!
!C******************************************************************MUMPS
!C Solution of the coarse problem KCG*UC = RESC (here stored in UC)
!         call mumps_resolve(mumps_struc,solc,lsolc)
!C******************************************************************MUMPS
!
!C Substructure correction into HAUX
!         if (myid.eq.0) then
!            write(*,*) '  substructure correction'
!         end if
!         call rzero(haux,lh)
!         do 2050 isub_loc = 1,nsub_loc
!            isub = start + isub_loc - 1
!
!            nnods  = nnodsa(isub)
!            lsols  = lsolsa(isub)
!            nnodis = nnodisa(isub)
!            nglbv  = nrhssa(isub) - lsolcsa(isub)
!            nreq1  = nreq1a(isub)
!            call cors(isub_loc, start, ideq1, lrec, name1,lname1,
!     *                nnodsx, lsolsx, nnodisx, lsolisx, lcfsx, nrhssx,
!     *                nnods, lsols, nnodis, nglbv,
!     *                nreq1, mfronx,
!     *                iinsn,liinsn,
!     *                nndfs,lnndfs, kdofs,lkdofs, ipiv,lipiv,
!     *                wi,lwi, cf,lcf, dualm,ldualm, scaling,lscaling,
!     *                dualrhs,ldualrhs, ifixcns,lifixcns, ibuf,libuf, 
!     *                gfixv,lgfixv, gload,lgload, rbuf,lrbuf, rea,lrea,
!     *                aslod,laslod, res,lres, haux,lh)
! 2050    end do
!
!C***************************************************************PARALLEL
!         call MPI_BCAST(solc, lsolc, MPI_DOUBLE_PRECISION, 0,comm,ierr)
!C***************************************************************PARALLEL
!
!C Coarse grid correction into HAUX
!         if (myid.eq.0) then
!            write(*,*) '  coarse grid correction'
!         end if
!         do 2052 isub_loc = 1,nsub_loc
!            isub = start + isub_loc - 1
!
!            lsolis = lsolisa(isub)
!            nnodcs = nnodcsa(isub)
!            nglbs  = nglbsa(isub)
!            nrhss  = nrhssa(isub)
!            call corc(isub_loc, start, nnodc,
!     *                nnodsx, lsolisx, nnodcsx,lphisisx,nglbsx,linglbsx,
!     *                lsolis, nnodcs, nglbs, nrhss, 
!     *                icngcn,licngcn, icnsn,licnsn,
!     *                nndfs,lnndfs, 
!     *                isglbgglbn,lisglbgglbn, inglbs,linglbs,
!     *                nnglbs,lnnglbs,
!     *                kdofc,lkdofc,
!     *                phisi,lphisi, wi,lwi,
!     *                solhelp1,lsolhelp1, solhelp2,lsolhelp2,
!     *                solc,lsolc, haux,lh)
! 2052    end do
!
!C Construction of global vector H
!C***************************************************************PARALLEL
!      call swapdata(myid, nsub_loc, nsub_locx, start, finish, comm,
!     *              nsub, nadjsx, nnodsx, nnodisx, lsolisx, 
!     *              nsnodsx, nsdofsx,
!     *              nadjsa,lnadjsa, nnodisa,lnnodisa, 
!     *              lsolisa,llsolisa, nsnodsa,lnsnodsa,
!     *              iadjs,liadjs, nsdofadj,lnsdofadj,
!     *              ishnin,lishnin, iinsn,liinsn,
!     *              commveco,commveci,lcommvec, nndfs,lnndfs,
!     *              kdofi,lkdofi,
!     *              h,haux,lh)
!C***************************************************************PARALLEL
!
!         rmpold = rmp
!
!C Scalar product of vectors RES * H - res*h => rmp
!         rmp_loc = 0.0D0
!         do 2060 isub_loc = 1,nsub_loc
!            isub = start + isub_loc - 1
!            
!            call scalpvv(lsolisa(isub), lsolisx, isub_loc, 
!     *                   res,lres, haux,lh, rmp_loc)
! 2060    end do
!C***************************************************************PARALLEL
!         call MPI_ALLREDUCE(rmp_loc,rmp, 1, MPI_DOUBLE_PRECISION,
!     *                      MPI_SUM, comm, ierr) 
!C***************************************************************PARALLEL
!
!C Determination of parameter BETA
!         beta = rmp/rmpold
!         ! debug
!!         print *, 'iter',iter,'rmp *************************',
!!     *              rmp
!         ! debug
!!         print *, 'iter',iter,'beta *************************',
!!     *              beta
!         
!C Determination of new step direction P
!         do 2070 isub_loc = 1,nsub_loc
!            isub = start + isub_loc - 1
!
!            ji = (isub_loc-1)*lsolisx
!           
!            lsolis = lsolisa(isub)
!            do i = 1,lsolis
!               p(ji + i) = h(ji + i) + beta*p(ji + i)
!            end do
! 2070    end do
!C Filling matrix for the Lanczos method
!         w((iter-1)*(maxit+1) + iter) = w((iter-1)*(maxit+1) + iter)
!     1                                + 1/alpha
!         w((iter)*(maxit+1) + iter + 1) = beta/alpha
!         w((iter-1)*(maxit+1) + iter + 1) = -sqrt(beta)/alpha
!         w((iter)*(maxit+1) + iter) = w((iter-1)*(maxit+1) + iter + 1)
! 1060 end do
!C*************************END OF MAIN CYCLE OVER ITERATIONS*************
!
!58    continue
!
!C Postprocessing of solution - computing interior values
!      if (myid.eq.0) then
!         write(*,*) 'Postprocessing of solution'
!      end if
!      do 1070 isub_loc = 1,nsub_loc
!         isub = start + isub_loc - 1
!         
!         nnods  = nnodsa(isub)
!         lsols  = lsolsa(isub)
!         nnodis = nnodisa(isub)
!         nnodos = nnodosa(isub)
!         nreq2  = nreq2a(isub)
!         call post(isub_loc, start, ideq2,lrec, idfvss, idrhss, idsols, 
!     *             name1,lname1,
!     *             nnodsx, lsolsx, nnodisx, lsolisx, 
!     *             nnodosx, lsolosx,
!     *             nnods, lsols, nnodis,
!     *             nreq2,mfronx,
!     *             iinsn,liinsn,
!     *             nndfs,lnndfs, kdofs,lkdofs,
!     *             wi,lwi, 
!     *             aslod,laslod,
!     *             ifixins,lifixins, ibuf,libuf,
!     *             gfixv,lgfixv, gload,lgload, 
!     *             rbuf,lrbuf, rea,lrea,
!     *             soli,lsoli)
! 1070 end do
!
!C 59   continue
!
!C End of mesuring time
!C***************************************************************PARALLEL
!      call MPI_BARRIER(comm,ierr)
!      time2 = MPI_WTIME()
!      tpcg = time2 - time1
!C***************************************************************PARALLEL
!
!C Condition number estimation on root processor
!      if (myid.eq.0) then
!         write(*,*) '================================================'
!         write(*,*) 'ESTIMATION OF CONDITION NUMBER BY LANCZOS METHOD'
!         nwx = maxit + 1
!         call condtri(nw,nwx,w,lw, cond)
!C         call condest(nw,nwx,w,wchol,lw,x,y, cond)
!         write(*,*) 'Condition number cond = ',cond
!         write(*,*) '================================================'
!      end if
!
!
!
!C***********************************************************************
!      subroutine colres(isub_loc, nnodc,
!     *                  lphisisx, nnodsx, lsolisx, nnodcsx,
!     *                  nglbsx, linglbsx, 
!     *                  lsolis, nnodcs, nglbs, nrhss,
!     *                  icngcn,licngcn, icnsn,licnsn, nndfs,lnndfs,
!     *                  isglbgglbn,lisglbgglbn, inglbs,linglbs,
!     *                  nnglbs,lnnglbs,
!     *                  kdofc,lkdofc,
!     *                  phisi,lphisi, wi,lwi,
!     *                  reshelp1,lreshelp1, reshelp2,lreshelp2,
!     *                  res,lres, solcaux,lsolc)
!C***********************************************************************
!C Subroutine for collecting residual in BDDC
!C***********************************************************************
!      
!      implicit none
!      
!      integer*4 isub_loc, nnodc, lphisisx, nnodsx, lsolisx, nnodcsx,
!     *          nglbsx, linglbsx, lsolis, nnodcs, nglbs, nrhss,
!     *          licngcn, licnsn, lnndfs, 
!     *          lisglbgglbn, linglbs, lnnglbs,
!     *          lkdofc,
!     *          lres, lsolc, lphisi, lwi, lreshelp1, lreshelp2, 
!     *          i, j, ji, jphisi, jrow, jicnsn, jnndfs,
!     *          jicngcn, jisglbgglbn, jinglbs, jnnglbs, indcs, indns, 
!     *          ndofn, indnc, indc, idofn,
!     *          iglbg, indinglbs, nglbn
!      integer*4 icngcn(licngcn), icnsn(licnsn), nndfs(lnndfs),
!     *          isglbgglbn(lisglbgglbn), inglbs(linglbs),
!     *          nnglbs(lnnglbs),
!     *          kdofc(lkdofc)
!      
!      real*8 res(lres), solcaux(lsolc), phisi(lphisi), wi(lwi), 
!     *       reshelp1(lreshelp1), reshelp2(lreshelp2)
!     
!      ji     = (isub_loc-1)*lsolisx
!      jphisi = (isub_loc-1)*lphisisx
!      jicnsn = (isub_loc-1)*nnodcsx
!      jnndfs = (isub_loc-1)*nnodsx
!      jicngcn = (isub_loc-1)*nnodcsx
!      jisglbgglbn = (isub_loc-1)*nglbsx
!      jinglbs = (isub_loc-1)*linglbsx
!      jnnglbs = (isub_loc-1)*nglbsx
!
!      do i = 1,lsolis
!         reshelp1(i) = wi(ji + i)*res(ji + i)
!      end do
!
!      call rzero(reshelp2,lreshelp2)
!
!      do i = 1,nrhss
!         jrow = (i-1)*lsolis
!         do j = 1,lsolis
!            if ((phisi(jphisi + jrow + j).ne.0.0D0).and.
!     *          (reshelp1(j).ne.0.0D0)) then
!               reshelp2(i)=reshelp2(i)+phisi(jphisi+jrow+j)*reshelp1(j)
!            end if   
!         end do
!      end do
!
!      indcs = 0
!      do i = 1,nnodcs
!         indns = icnsn(jicnsn + i)
!         ndofn = nndfs(jnndfs + indns)
!         indnc = icngcn(jicngcn + i)
!         indc  = kdofc(indnc)
!         
!         do idofn = 1,ndofn
!            indcs = indcs + 1
!            indc  = indc + 1
!            
!            solcaux(indc) = solcaux(indc) + reshelp2(indcs)
!         end do
!      end do   
!      indinglbs = 0
!      do i = 1,nglbs
!         iglbg = isglbgglbn(jisglbgglbn + i)
!         indns = inglbs(jinglbs + indinglbs + 1)
!         ndofn = nndfs(jnndfs + indns)
!         nglbn = nnglbs(jnnglbs + i)
!         indc  = kdofc(nnodc + iglbg)
!         
!         do idofn = 1,ndofn
!            indcs = indcs + 1
!            indc  = indc + 1
!            
!            solcaux(indc) = solcaux(indc) + reshelp2(indcs)
!         end do
!         indinglbs = indinglbs + nglbn
!      end do   
!
!      return
!      end
!
!C***********************************************************************
!      subroutine smmult(isub_loc, start, ideq2,lrec,name1,lname1,
!     *                  nnodsx, lsolsx, nnodisx, lsolisx,
!     *                  nnods, lsols, nnodis, nreq2, mfronx,
!     *                  iinsn,liinsn, nndfs,lnndfs,
!     *                  ifixins,lifixins, kdofs,lkdofs, ibuf,libuf,
!     *                  gfixv,lgfixv, gload,lgload, 
!     *                  rbuf,lrbuf, rea,lrea, aslod,laslod, 
!     *                  p,lp, ap,lap)
!C***********************************************************************
!C Subroutine for multiplication by system matrix in BDDC
!C***********************************************************************
!      use module_utils
!
!      implicit none
!      
!      integer*4 isub_loc, start, ideq2,lrec, lname1, 
!     *          nnodsx, lsolsx, nnodisx, lsolisx, 
!     *          nnods, lsols, nnodis, nreq2, mfronx,
!     *          liinsn, lnndfs, lifixins, lkdofs,
!     *          lp, lap,
!     *          libuf, lgfixv, lgload, lrbuf,
!     *          laslod, lrea,
!     *          ji, jiinsn, jnndfs,
!     *          isub, ind, nrhs, indi, ii, ins, ndofn,
!     *          idofn, in
!      integer*4 iinsn(liinsn), nndfs(lnndfs),
!     *          ifixins(lifixins), kdofs(lkdofs), ibuf(libuf)
!
!      real*8 p(lp), ap(lap), 
!     *       aslod(laslod),
!     *       gfixv(lgfixv), gload(lgload), rbuf(lrbuf),
!     *       rea(lrea)
!      
!      character name1*100, fname*130
!
!      isub = start + isub_loc - 1
!
!      ji     = (isub_loc-1)*lsolisx
!      jiinsn = (isub_loc-1)*nnodisx
!      jnndfs = (isub_loc-1)*nnodsx
!
!C Creation of field kdofs
!      kdofs(1) = 0
!      do in = 2,nnods
!         kdofs(in) = kdofs(in-1) + nndfs(jnndfs + in-1)
!      end do
!
!C Prepare fixed variables      
!      call rzero(gfixv,lgfixv)
!      indi = 0
!      do ii = 1,nnodis
!         ins = iinsn(jiinsn + ii)
!         ndofn = nndfs(jnndfs + ins)
!         ind = kdofs(ins)
!         do idofn = 1,ndofn
!            ind = ind + 1
!            indi = indi + 1
!            gfixv(ind) = p(ji + indi)
!         end do
!      end do
!
!C Prepare RHS and reactions
!      call rzero(aslod,laslod)
!      call rzero(rea,lrea)
!
!      call getfname(name1,isub,'EQ2',fname)
!      open(unit=ideq2, file=fname, status='old',
!     *     form='unformatted', access='direct', recl=lrec)
!
!C Solve subdomain problem
!      nrhs = 1
!      call frors(isub_loc, mfronx, ideq2,
!     *           lsolsx, lsols, nrhs,       
!     *           ibuf,libuf, rbuf,lrbuf, gload,lgload,
!     *           aslod,laslod, rea,lrea, ifixins,lifixins, 
!     *           gfixv,lgfixv,nreq2)
!      close(ideq2)
!
!      indi = 0
!      do ii = 1,nnodis
!         ins = iinsn(jiinsn + ii)
!         ndofn = nndfs(jnndfs + ins)
!         ind = kdofs(ins)
!         do idofn = 1,ndofn
!            ind = ind + 1
!            indi = indi + 1
!            ap(ji + indi) = -rea(ind)
!         end do
!      end do      
!      
!      return
!      end
!
!C***********************************************************************
!      subroutine post(isub_loc, start, ideq2,lrec, idfvss,idrhss,idsols,
!     *                name1,lname1,
!     *                nnodsx, lsolsx, nnodisx, lsolisx, 
!     *                nnodosx, lsolosx,
!     *                nnods, lsols, nnodis,
!     *                nreq2,mfronx,
!     *                iinsn,liinsn,
!     *                nndfs,lnndfs, kdofs,lkdofs,
!     *                wi,lwi, 
!     *                aslod,laslod,
!     *                ifixins,lifixins, ibuf,libuf,
!     *                gfixv,lgfixv, gload,lgload, 
!     *                rbuf,lrbuf, rea,lrea,
!     *                soli,lsoli)
!C***********************************************************************
!C Subroutine for postprocessing after solution by BDDC
!C Computes interior values of staticly condensed solution
!C***********************************************************************
!      use module_utils
!      
!      implicit none
!      
!      integer*4 isub_loc, start, ideq2,lrec, idfvss, idrhss, idsols, 
!     *          lname1,
!     *          nnodsx, lsolsx, nnodisx, lsolisx, nnodosx, lsolosx,
!     *          nnods, lsols, nnodis,
!     *          nreq2, mfronx,
!     *          liinsn, lnndfs, lkdofs, 
!     *          lsoli,
!     *          laslod, lifixins, libuf, 
!     *          lgfixv, lgload, 
!     *          lrbuf,
!     *          lwi, lrea,
!     *          ji, jo, jisngn, jiinsn, jnndfs, jionsn, 
!     *          jifixns, 
!     *          ind, ii,
!     *          indi, inds,
!     *          nrhs, isub, ins, ndofn, idofn, in, i
!      integer*4 iinsn(liinsn),
!     *          nndfs(lnndfs), kdofs(lkdofs),
!     *          ifixins(lifixins), 
!     *          ibuf(libuf)
!
!      real*8 soli(lsoli), wi(lwi), 
!     *       aslod(laslod),
!     *       gfixv(lgfixv), gload(lgload), rbuf(lrbuf),
!     *       rea(lrea)
!     
!      character name1*100, fname*130
!
!      isub = start + isub_loc - 1
!
!      ji       = (isub_loc-1)*lsolisx
!      jo       = (isub_loc-1)*lsolosx
!      jisngn   = (isub_loc-1)*nnodsx
!      jiinsn   = (isub_loc-1)*nnodisx
!      jionsn   = (isub_loc-1)*nnodosx
!      jnndfs   = (isub_loc-1)*nnodsx
!      jifixns  = (isub_loc-1)*lsolsx
!
!
!C Creation of field kdofs
!      kdofs(1) = 0
!      do in = 2,nnods
!         kdofs(in) = kdofs(in-1) + nndfs(jnndfs + in-1)
!      end do
!
!C Prepare fixed variables
!      call getfname(name1,isub,'FVSS',fname)
!      open(unit=idfvss,file=fname,status='old',form='unformatted')
!      rewind idfvss
!C GFIXV is vector of fixed variables - here only first column is used
!      call rzero(gfixv,lgfixv)
!      do i = 1,lsols
!         read(idfvss) gfixv(i)
!      end do
!      close(idfvss)
!      indi = 0
!      do ii = 1,nnodis
!         ins = iinsn(jiinsn + ii)
!         ndofn = nndfs(jnndfs + ins)
!         ind = kdofs(ins)
!         do idofn = 1,ndofn
!            ind = ind + 1
!            indi = indi + 1
!            gfixv(ind) = soli(ji+indi)
!         end do
!      end do
!
!C Prepare RHS
!      call getfname(name1,isub,'RHSS',fname)
!      open(unit=idrhss,file=fname,status='old',form='unformatted')
!      rewind idrhss
!      call rzero(aslod,laslod)
!      do i = 1,lsols
!         read(idrhss) aslod(i)
!      end do
!      close(idrhss)
!
!      call getfname(name1,isub,'EQ2',fname)
!      open(unit=ideq2, file=fname, status='old',
!     *     form='unformatted', access='direct', recl=lrec)
!
!C Solve subdomain problem
!      nrhs = 1
!      call frors(isub_loc, mfronx, ideq2,
!     *           lsolsx, lsols, nrhs,       
!     *           ibuf,libuf, rbuf,lrbuf, gload,lgload,
!     *           aslod,laslod, rea,lrea, ifixins,lifixins, 
!     *           gfixv,lgfixv,nreq2)
!      close(ideq2)
!
!C Write solution in interiors into subdomain solution file
!C interface solution is weighted
!      indi = 0
!      do ii = 1,nnodis
!         ins   = iinsn(jiinsn + ii)
!         ndofn = nndfs(jnndfs + ins)
!         inds  = kdofs(ins)
!         do idofn = 1,ndofn
!            indi = indi + 1
!            inds = inds + 1
!            aslod(inds) = wi(ji + indi)*aslod(inds)
!         end do
!      end do
!      call getfname(name1,isub,'SOLS',fname)
!      open(unit=idsols,file=fname,status='unknown',
!     *     form='unformatted')
!      rewind idsols
!      write(idsols) (aslod(i),i=1,lsols)
!      close(idsols)
!
!      return
!      end
!
!C***********************************************************************
!      subroutine swapdata(myid, nsub_loc, nsub_locx, start, finish,comm,
!     *                    nsub, nadjsx, nnodsx, nnodisx, lsolisx, 
!     *                    nsnodsx, nsdofsx,
!     *                    nadjsa,lnadjsa, nnodisa,lnnodisa, 
!     *                    lsolisa,llsolisa, nsnodsa,lnsnodsa,
!     *                    iadjs,liadjs, nsdofadj,lnsdofadj,
!     *                    ishnin,lishnin, iinsn,liinsn,
!     *                    commveco,commveci,lcommvec, nndfs,lnndfs,
!     *                    kdofi,lkdofi,
!     *                    vec,vecaux,lvec)
!C***********************************************************************
!C Subroutine for SEND/RECV communication
!C***********************************************************************
!      
!      implicit none
!      
!      include "mpif.h"
!
!      integer*4 myid, ierr, comm, tag
!      integer*4 nsub_loc, nsub_locx, start, finish, 
!     *          nsub, nadjsx, nnodsx, nnodisx, lsolisx, nsnodsx,nsdofsx,
!     *          lnadjsa, lnnodisa, llsolisa, lnsnodsa, lvec, 
!     *          liadjs, lnsdofadj,
!     *          lishnin, liinsn, lcommvec, lnndfs, lkdofi,
!     *          jiinsn, jnndfs, jishnin, jcommvec, jvec,
!     *          nnodis, ndofn, nsnods, nadjs, nadjsneib, nsdofadjs,
!     *          lsolis,
!     *          i, j, k, isub_loc, iproc, isub, isubneib, isubneib_loc,
!     *          indn, indcommvec, indni, idofn, indi,  
!     *          ind, ind2, indneib, ireq, nreq
!      integer*4 request(2*nsub_locx*nadjsx)
!      integer*4 statarray(MPI_STATUS_SIZE, 2*nsub_locx*nadjsx)
!      integer*4 nadjsa(lnadjsa), nnodisa(lnnodisa), lsolisa(llsolisa), 
!     *          nsnodsa(lnsnodsa),
!     *          iadjs(liadjs), nsdofadj(lnsdofadj),
!     *          ishnin(lishnin),
!     *          iinsn(liinsn), nndfs(lnndfs), kdofi(lkdofi)
!
!      real*8 vec(lvec), vecaux(lvec),
!     *       commveco(lcommvec), commveci(lcommvec)
!      
!C Prepare vector for communication
!      do 1000 isub_loc = 1,nsub_loc
!         isub = start + isub_loc - 1
!
!         jiinsn    = (isub_loc-1)*nnodisx
!         jnndfs    = (isub_loc-1)*nnodsx
!         jishnin   = (isub_loc-1)*nsnodsx
!         jcommvec  = (isub_loc-1)*nsdofsx
!         jvec      = (isub_loc-1)*lsolisx
!
!C Creation of field kdofi
!         kdofi(1) = 0
!         nnodis = nnodisa(isub)
!         do i = 2,nnodis
!            indn  = iinsn(jiinsn + i)
!            ndofn = nndfs(jnndfs + indn)
!            kdofi(i) = kdofi(i-1) + ndofn
!         end do
!
!         nsnods = nsnodsa(isub)
!         indcommvec = jcommvec
!         do i = 1,nsnods
!            indni = ishnin(jishnin + i)
!            indn  = iinsn(jiinsn + indni)
!            ndofn = nndfs(jnndfs + indn)
!            do idofn = 1,ndofn
!               indcommvec = indcommvec + 1
!               indi = jvec + kdofi(indni) + idofn
!               commveco(indcommvec) = vecaux(indi)
!            end do
!         end do
! 1000 end do
!
!C Swap data
!      ireq = 1
!      do 1010 isub_loc = 1,nsub_loc
!         isub = start + isub_loc - 1
!         
!         jiinsn    = (isub_loc-1)*nnodisx
!         jnndfs    = (isub_loc-1)*nnodsx
!         jishnin   = (isub_loc-1)*nsnodsx
!         jcommvec  = (isub_loc-1)*nsdofsx
!         jvec      = (isub_loc-1)*lsolisx
!
!         ind = (isub_loc-1)*nsdofsx
!         nadjs = nadjsa(isub)
!         do i = 1,nadjs
!            isubneib = iadjs((isub_loc-1)*nadjsx + i)
!            if ((isubneib.ge.start).and.(isubneib.le.finish)) then
!               isubneib_loc = isubneib - start + 1
!               indneib = (isubneib_loc-1)*nsdofsx
!               nadjsneib = nadjsa(isubneib)
!               do j = 1,nadjsneib
!                  if (iadjs((isubneib_loc-1)*nadjsx + j).eq.isub) then
!                     nsdofadjs = nsdofadj((isub_loc-1)*nadjsx + i)
!                     do k = 1,nsdofadjs
!                        commveci(indneib + k) = commveco(ind + k)
!                     end do
!                     if (nsdofadj((isub_loc-1)*nadjsx + i).ne.
!     *                   nsdofadj((isubneib_loc-1)*nadjsx + j)) then
!                        write(*,*)'myid =',myid,': nsdofadj does not
!     *                             match!'
!                        call MPI_ABORT(comm, 78, ierr)
!                     end if   
!                  end if
!                  indneib = indneib+nsdofadj((isubneib_loc-1)*nadjsx+j)
!               end do
!            else
!               iproc = (isubneib-1)/nsub_locx
!
!c               write(*,*) 'I am',myid,' working on sub',isub,
!c     *                   ' sending to',iproc,'about sub',
!c     *                    isubneib,' with id ',isubneib*isub
!
!C Determination of tag for communication
!               tag = isub*nsub + isubneib
!  
!               call MPI_ISEND(commveco(ind + 1),
!     *                        nsdofadj((isub_loc-1)*nadjsx + i),
!     *                        MPI_DOUBLE_PRECISION, iproc,
!     *                        tag, comm,
!     *                        request(ireq), ierr)
!
!               ireq = ireq + 1
!
!c               write(*,*) 'I am',myid,' working on sub',isub,
!c     *               ' recieving from',iproc,'about sub',
!c     *                 isub,' with id',isubneib*isub
!   
!
!C Determination of tag for communication
!               tag = isubneib*nsub + isub
!   
!               call MPI_IRECV(commveci(ind + 1),
!     *                        nsdofadj((isub_loc-1)*nadjsx + i),
!     *                        MPI_DOUBLE_PRECISION, iproc, 
!     *                        tag, comm,
!     *                        request(ireq), ierr)
!               
!               ireq = ireq + 1
!            end if
!            ind = ind + nsdofadj((isub_loc-1)*nadjsx + i)
!         end do
! 1010 end do
!      
!      nreq = ireq-1
!       
!      call MPI_WAITALL(nreq, request, statarray, ierr)
! 
!C Copy obtained data to vector VEC
!      call rzero(vec,lvec)
!      ind  = 0
!      ind2 = 0
!      do 1020 isub_loc = 1,nsub_loc
!         isub = start + isub_loc - 1
!         
!         jiinsn    = (isub_loc-1)*nnodisx
!         jnndfs    = (isub_loc-1)*nnodsx
!         jishnin   = (isub_loc-1)*nsnodsx
!         jcommvec  = (isub_loc-1)*nsdofsx
!         jvec      = (isub_loc-1)*lsolisx
!
!C Creation of field kdofi
!         kdofi(1) = 0
!         nnodis = nnodisa(isub)
!         do i = 2,nnodis
!            indn  = iinsn(jiinsn + i)
!            ndofn = nndfs(jnndfs + indn)
!            kdofi(i) = kdofi(i-1) + ndofn
!         end do
!
!         nsnods = nsnodsa(isub)
!         indcommvec = jcommvec
!         do i = 1,nsnods
!            indni = ishnin(jishnin + i)
!            indn  = iinsn(jiinsn + indni)
!            ndofn = nndfs(jnndfs + indn)
!            do idofn = 1,ndofn
!               indcommvec = indcommvec + 1
!               indi = jvec + kdofi(indni) + idofn
!               vec(indi) = vec(indi) + commveci(indcommvec)
!            end do
!         end do
!         lsolis = lsolisa(isub)
!         do i = 1,lsolis
!            vec(jvec + i) = vec(jvec + i) + vecaux(jvec + i)
!         end do
! 1020 end do
!
!      return
!      end
!
!C***********************************************************************
!      subroutine scalpvv(n, nx, isub_loc, a,la, b,lb, c)
!C***********************************************************************
!C     Computes scalar product of vectors A(NX) and B(NX) and saves    
!C     the result into C
!C***********************************************************************
!
!      implicit none
!      
!      integer*4 n, nx, isub_loc, la, lb, i, ind
!
!      real*8 c
!      real*8 a(la),b(lb)
!
!      ind = (isub_loc-1)*nx
!      do i = 1,n
!         if ((a(ind + i).ne.0.0D0).and.(b(ind + i).ne.0.0D0)) then
!            c = c + a(ind + i)*b(ind + i)
!         end if   
!      end do
!
!      return
!      end
!
!      SUBROUTINE FRODS(ISUB_LOC, MW, ICL, NDOFN, MFRON, PIVAL, 
!     *                 IDELM, IDEQ,
!     *                 NELEMX, LSOLX, LINETX,
!     *                 NELEM, LSOL, NNOD, LLINET,
!     *                 INET,LINET, NNET,LNNET, IFIX,LIFIX,
!     *                 LOCEL,NDEST,ILOC,NEVAX, 
!     *                 NACVA,LNACVA, ELM,LELMX,
!     *                 GSTIF,LGSTIF, IBUF,LIBUF, RBUF,LRBUF, NREQ)
!C*********************************************************************
!C     FRODS (FRONTAL - DIRECT SOLUTION) - PRIMY  CHOD FRONTALNIHO    *
!C ALGORITMU PRO STANDARTNI OKRAJOVOU ULOHU S NENULOVYMI              *
!C PREDEPSANYMI POSUVY                                                *
!C     INPUT: MW,ICL,MFRON,PIVAL,IDGMI,IDP,IDELM,IDEQ1,NELEM,NNOD,    *
!C INET(LINET),NNET(LNNET),IFIX(LIFIX)                                *
!C     OUTPUT: -                                                      *
!C*********************************************************************
!C
!      IMPLICIT  NONE
!      
!      integer*4 ISUB_LOC, MW, ICL, NDOFN, MFRON, IDELM, IDEQ,
!     *          NELEM, LSOL, NNOD, LLINET,
!     *          NELEMX, LSOLX, LINETX,
!     *          LINET, LNNET, LIFIX, NEVAX, LNACVA, LELMX, LGSTIF,
!     *          LIBUF, LRBUF, NREQ,
!     *          pointifix, pointinet, pointnnet,
!     *          truelifix, truelinet, truelnnet
!      integer*4 INET(LINET),NNET(LNNET),IFIX(LIFIX),LOCEL(NEVAX),
!     *          NDEST(NEVAX),ILOC(NEVAX),NACVA(LNACVA), IBUF(LIBUF)
!     
!       real*8 PIVAL
!       real*8 ELM(LELMX), GSTIF(LGSTIF), RBUF(LRBUF)
!C
!C
!      pointifix = (ISUB_LOC-1)*LSOLX + 1
!      pointinet = (ISUB_LOC-1)*LINETX + 1
!      pointnnet = (ISUB_LOC-1)*NELEMX + 1
!
!      truelifix = LSOL
!      truelinet = LLINET
!      truelnnet = NELEM
!
!C
!C     ** PRIMY CHOD FRONTALNIHO ALGORITMU **
!      CALL DIRF(MW,ICL,NDOFN,MFRON,PIVAL,IDELM,IDEQ,NELEM,NNOD,
!     *          INET(pointinet),truelinet, 
!     *          NNET(pointnnet),truelnnet, 
!     *          IFIX(pointifix),truelifix,
!     *          LOCEL,NDEST,ILOC,NEVAX,
!     *          NACVA,LNACVA, ELM,LELMX, GSTIF,LGSTIF,
!     *          IBUF,LIBUF ,RBUF,LRBUF, NREQ)
!C
!      RETURN
!      END
!
!      SUBROUTINE FRORS(ISUB_LOC, MFRON, IDEQ,
!     *                 LSOLSX, LSOL, NRHS, 
!     *                 IBUF,LIBUF, RBUF,LRBUF, GLOAD,LGLOAD,
!     *                 ASLOD,LASLOD, REA,LREA, IFIXS,LIFIXS, 
!     *                 GFIXV,LGFIXV, NREQ)
!C*********************************************************************
!C     FRORS (FRONTAL - RESOLUTION AND BACKSUBSTITUTION) ZPETNY CHOD  *
!C FRONTALNIHO ALGORITMU PRO STANDARTNI OKRAJOVOU ULOHU S NENULOVYMI  *
!C FIXOVANYMI PROMENNYMI, RUZNYMI PRO KAZDOU PRAVOU STRANU            *
!C     INPUT: MW,MFRON,NRHS,IDP,IDFV,IDEQ1,NREQ1                      *
!C     OUTPUT: ASLOD(NLSOL)                                           *
!C*********************************************************************
!C
!      IMPLICIT  NONE
!      
!      integer*4 ISUB_LOC, MFRON, LFIXV, IDEQ,
!     *          LSOLSX, LSOL, NRHS, 
!     *          LIBUF, LRBUF, LIFIXS, LGFIXV, LGLOAD, LASLOD, LREA,NREQ,
!     *          pointifixs, truelifixs, truelgload, truelaslod, 
!     *          truelgfixv
!      integer*4 IBUF(LIBUF), IFIXS(LIFIXS)
!      
!      real*8 RBUF(LRBUF), GFIXV(LGFIXV), GLOAD(LGLOAD), ASLOD(LASLOD),
!     *       REA(LREA) 
!
!      pointifixs = (ISUB_LOC-1)*LSOLSX + 1
!      truelifixs = LSOL
!
!      LFIXV = LSOL
!
!      truelgload = NRHS * MFRON
!      truelaslod = NRHS * LSOL
!      truelgfixv = NRHS * LFIXV
!
!C     ** ZPETNY CHOD FRONTALNIHO ALGORITMU **
!      CALL RESF(MFRON,NRHS,IDEQ, 
!     *          IBUF,LIBUF, RBUF,LRBUF,
!     *          GLOAD,truelgload, ASLOD,REA,truelaslod,
!     *          IFIXS(pointifixs),truelifixs, 
!     *          GFIXV,truelgfixv, LFIXV, NREQ)
!
!      RETURN
!      END
!
!C***********************************************************************
!      subroutine izero(l,ll)
!C***********************************************************************
!C Subroutine which zero integer field L of lenght LL
!C***********************************************************************
!
!      implicit none
!      
!      integer*4 ll, i
!      integer*4 l(ll)
!
!      do i = 1,ll
!         l(i) = 0
!      end do
!
!      return
!      end
!
!C***********************************************************************
!      subroutine rzero(r,lr)
!C***********************************************************************
!C Subroutine which zero doble precision field R of lenght LR
!C***********************************************************************
!
!      implicit none
!      
!      integer*4 lr, i
!      
!      real*8 r(lr)
!
!      do i = 1,lr
!         r(i) = 0.0D0
!      end do
!
!      return
!      end
!
!C***********************************************************************
!      subroutine cors(isub_loc, start, ideq1, lrec, name1,lname1,
!     *               nnodsx, lsolsx, nnodisx, lsolisx, lcfsx, nrhssx,
!     *               nnods, lsols, nnodis, nglbv,
!     *               nreq1, mfronx,
!     *               iinsn,liinsn,
!     *               nndfs,lnndfs, kdofs,lkdofs, ipiv,lipiv,
!     *               wi,lwi, cf,lcf, dualm,ldualm, scaling,lscaling,
!     *               dualrhs,ldualrhs, ifixcns,lifixcns, ibuf,libuf, 
!     *               gfixv,lgfixv, gload,lgload, rbuf,lrbuf, rea,lrea,
!     *               aslod,laslod,res,lres, v,lv)
!C***********************************************************************
!C Subroutine for substructure correction in BDDC
!C***********************************************************************
!
!      use module_utils
!
!      implicit none
!      
!      integer*4 isub_loc, start, ideq1, lrec, lname1, 
!     *          nnodsx, lsolsx, nnodisx, lsolisx, lcfsx, nrhssx,
!     *          nnods, lsols, nnodis, nglbv,
!     *          nreq1, mfronx,
!     *          liinsn, lnndfs, lkdofs, lipiv,
!     *          lifixcns, libuf, lres, lv, lwi, lcf, ldualm, lscaling,
!     *          ldualrhs, lgfixv, lgload, lrbuf, lrea, laslod,
!     *          ind, isub, nrhs, 
!     *          ji, jiinsn, jnndfs, jcf, jdualm, jipiv,
!     *          ndofn, idofn, indi, ii, ins, in, i, j, iglbv, info, ldim
!      integer*4 iinsn(liinsn), nndfs(lnndfs), kdofs(lkdofs),
!     *          ipiv(lipiv),
!     *          ifixcns(lifixcns), ibuf(libuf)
!      
!      real*8 value1, value2, sum
!      real*8 res(lres), v(lv), wi(lwi), cf(lcf), dualm(ldualm),
!     *       scaling(lscaling), dualrhs(ldualrhs),
!     *       gfixv(lgfixv), gload(lgload), rbuf(lrbuf),
!     *       rea(lrea), aslod(laslod)
!
!      character name1*100, fname*130
!
!      isub = start + isub_loc - 1
!
!      ji      = (isub_loc-1)*lsolisx
!      jiinsn  = (isub_loc-1)*nnodisx
!      jnndfs  = (isub_loc-1)*nnodsx
!      jcf     = (isub_loc-1)*lcfsx
!      jdualm  = (isub_loc-1)*nrhssx*nrhssx
!      jipiv   = (isub_loc-1)*nrhssx
!
!C Creation of field kdofs
!      kdofs(1) = 0
!      do in = 2,nnods
!         kdofs(in) = kdofs(in-1) + nndfs(jnndfs + in-1)
!      end do
!
!C*************************************************************************
!C Step 1 of the algorithm
!C Backward step of frontal solver with Res as RHS and zero fixed variables
!C to find Kff^-1*Res
!      call rzero(aslod,laslod)
!      indi = 0
!      do ii = 1,nnodis
!         ins = iinsn(jiinsn + ii)
!         ndofn = nndfs(jnndfs + ins)
!         ind = kdofs(ins)
!         do idofn = 1,ndofn
!            ind = ind + 1
!            indi = indi + 1
!            aslod(ind) = wi(ji + indi)*res(ji + indi)
!         end do
!      end do
!
!C Create vector of fixed variables GFIXV
!      call rzero(gfixv,lgfixv)
!
!      call getfname(name1,isub,'EQ1',fname)
!      open(unit=ideq1, file=fname, status='old',
!     *     form='unformatted', access='direct', recl=lrec)
!
!      nrhs = 1
!      call frors(isub_loc, mfronx, ideq1,
!     *           lsolsx, lsols, nrhs,       
!     *           ibuf,libuf, rbuf,lrbuf, gload,lgload,
!     *           aslod,laslod, rea,lrea, ifixcns,lifixcns, 
!     *           gfixv,lgfixv,nreq1)
!C In aslod is now Kff^-1*Res
!     
!C Step 2 of the algorithm
!C Multiply by Cf to find Cf*Kff^-1Res
!      do iglbv = 1,nglbv
!         sum = 0.0D0
!         do j = 1,lsols
!            value1 = cf(jcf + (iglbv-1)*lsols + j)
!            value2 = aslod(j)
!            if (value1.ne.0.0D0.and.value2.ne.0.0D0) then
!               sum = sum + value1*value2
!            end if
!         end do
!         dualrhs(iglbv) = sum
!      end do
!
!C Step 3 of the algorithm
!C Solve the dual problem using LAPACK routine
!      do iglbv = 1,nglbv
!         dualrhs(iglbv) = scaling(isub_loc) * dualrhs(iglbv)
!      end do
!      nrhs = 1
!C LAPACK patch - LAPACK does not support zero leading edge which
!C                is the case for sudomain without globs
!      LDIM = MAX(1,NGLBV)
!      call DGETRS( 'N', nglbv, nrhs, dualm(jdualm+1), LDIM,
!     *             ipiv(jipiv+1), dualrhs, LDIM, info)
!      if (info.ne.0) then
!         write(*,*) 'Error in Lapack solution of dual problem'
!         stop
!      end if
!
!C Step 4 of the algorithm
!C Multiply by -CfT to get -CfT*mu and add Res
!      call rzero(aslod,laslod)
!      do i = 1,lsols
!         sum = 0.0D0
!         do iglbv = 1,nglbv
!            value1 = cf(jcf + (iglbv-1)*lsols + i)
!            value2 = dualrhs(iglbv)
!            if (value1.ne.0.0D0.and.value2.ne.0.0D0) then
!               sum = sum + value1*value2
!            end if
!         end do
!         aslod(i) = -sum
!      end do
!
!      indi = 0
!      do ii = 1,nnodis
!         ins = iinsn(jiinsn + ii)
!         ndofn = nndfs(jnndfs + ins)
!         ind = kdofs(ins)
!         do idofn = 1,ndofn
!            ind = ind + 1
!            indi = indi + 1
!            aslod(ind) = aslod(ind) + wi(ji + indi)*res(ji + indi)
!         end do
!      end do
!
!C Step 5 of the algorithm
!C Backward step of frontal solver with RHS = -CfT*mu + Res
!      nrhs = 1
!      call frors(isub_loc, mfronx, ideq1,
!     *           lsolsx, lsols, nrhs,       
!     *           ibuf,libuf, rbuf,lrbuf, gload,lgload,
!     *           aslod,laslod, rea,lrea, ifixcns,lifixcns, 
!     *           gfixv,lgfixv,nreq1)
!      close(ideq1)
!C*************************************************************************
!
!      indi = 0
!      do ii = 1,nnodis
!         ins = iinsn(jiinsn + ii)
!         ndofn = nndfs(jnndfs + ins)
!         ind = kdofs(ins)
!         do idofn = 1,ndofn
!            ind = ind + 1
!            indi = indi + 1
!            v(ji + indi) = v(ji + indi) + wi(ji + indi)*aslod(ind)
!         end do
!      end do
!
!      return
!      end
!
!C***********************************************************************
!      subroutine corc(isub_loc, start, nnodc,
!     *               nnodsx, lsolisx, nnodcsx,lphisisx, nglbsx,linglbsx,
!     *               lsolis, nnodcs, nglbs, nrhss, 
!     *               icngcn,licngcn, icnsn,licnsn,
!     *               nndfs,lnndfs, 
!     *               isglbgglbn,lisglbgglbn, inglbs,linglbs,
!     *               nnglbs,lnnglbs,
!     *               kdofc,lkdofc,
!     *               phisi,lphisi, wi,lwi,
!     *               solhelp1,lsolhelp1, solhelp2,lsolhelp2,
!     *               solc,lsolc, v,lv)
!C***********************************************************************
!C Subroutine for coarse grid correction in BDDC
!C***********************************************************************
!
!      implicit none
!      
!      integer*4 isub_loc, start, nnodc,
!     *          nnodsx, lsolisx, nnodcsx, lphisisx, nglbsx, linglbsx,
!     *          lsolis, nnodcs, nglbs, nrhss,
!     *          licngcn, licnsn, lnndfs, 
!     *          lisglbgglbn, linglbs, lnnglbs,
!     *          lkdofc, 
!     *          lsolc, lv,
!     *          lphisi, lwi,
!     *          lsolhelp1, lsolhelp2,
!     *          i, j, isub, 
!     *          jphisi, ji, jicngcn, jicnsn, jnndfs, 
!     *          jisglbgglbn, jinglbs, jnnglbs, jrow,
!     *          indcs, indns, ndofn, indnc, indc, idofn,
!     *          iglbg, indinglbs, nglbn
!      integer*4 icngcn(licngcn), icnsn(licnsn),
!     *          nndfs(lnndfs),
!     *          isglbgglbn(lisglbgglbn), inglbs(linglbs),
!     *          nnglbs(lnnglbs),
!     *          kdofc(lkdofc)
!      
!      real*8 solc(lsolc), v(lv), 
!     *       phisi(lphisi),  wi(lwi), 
!     *       solhelp1(lsolhelp1), solhelp2(lsolhelp2)
!
!      isub = start + isub_loc - 1
!
!      jphisi  = (isub_loc-1)*lphisisx
!      ji      = (isub_loc-1)*lsolisx
!      jicnsn  = (isub_loc-1)*nnodcsx
!      jnndfs  = (isub_loc-1)*nnodsx
!      jicngcn = (isub_loc-1)*nnodcsx
!      jisglbgglbn = (isub_loc-1)*nglbsx
!      jinglbs = (isub_loc-1)*linglbsx
!      jnnglbs = (isub_loc-1)*nglbsx
!
!C Coarse grid correction
!      indcs = 0
!      do i = 1,nnodcs
!         indns = icnsn(jicnsn + i)
!         ndofn = nndfs(jnndfs + indns)
!         indnc = icngcn(jicngcn + i)
!         indc  = kdofc(indnc)
!         
!         do idofn = 1,ndofn
!            indcs = indcs + 1
!            indc  = indc + 1
!            
!            solhelp1(indcs) = solc(indc)
!         end do
!      end do
!      indinglbs = 0
!      do i = 1,nglbs
!         iglbg = isglbgglbn(jisglbgglbn + i)
!         indns = inglbs(jinglbs + indinglbs + 1)
!         ndofn = nndfs(jnndfs + indns)
!         nglbn = nnglbs(jnnglbs + i)
!         indc  = kdofc(nnodc + iglbg)
!         
!         do idofn = 1,ndofn
!            indcs = indcs + 1
!            indc  = indc + 1
!            
!            solhelp1(indcs) = solc(indc)
!         end do
!         indinglbs = indinglbs + nglbn
!      end do   
!
!      call rzero(solhelp2,lsolhelp2)
!
!      do i = 1,lsolis
!         do j = 1,nrhss
!            jrow = (j-1)*lsolis
!            if ((phisi(jphisi + jrow + i).ne.0.0D0).and.
!     *          (solhelp1(j).ne.0.0D0)) then
!               solhelp2(i) = solhelp2(i) 
!     *                     + phisi(jphisi + jrow + i)*solhelp1(j)
!            end if   
!         end do
!      end do
!
!      do i = 1,lsolis
!         v(ji + i) = v(ji + i) + wi(ji + i)*solhelp2(i)
!      end do
!
!      return
!      end
!
!C***********************************************************************
!      subroutine showlenght(namevar,lenght,nbytes)
!C***********************************************************************
!C Computes size of field in Megabytes and print the value at stdout
!C***********************************************************************
!
!      implicit none
!
!      integer*4 lenght,nbytes
!      real*8 lenghtmb
!      
!      character namevar*(*), realint*7
!
!C set type of array
!      if(nbytes.eq.4) then
!         realint = 'INTEGER'
!      else if(nbytes.eq.8) then
!         realint = 'REAL   '
!      else
!         realint = 'STRANGE'
!      end if
!
!C determine size in megabytes
!      lenghtmb = lenght*nbytes/1048576.0d0
!C
!      write(*,7000) trim(realint),namevar,lenght,lenghtmb
! 7000 format(1X,'Lenght of ',a,' array ',a,' is ',i10,',i.e.',
!     *       f10.3,' MB.')
!C
!      end
!
!
