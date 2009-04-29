module module_bddc
!*****************
! Module for realization of the BDDC preconditioner
! Jakub Sistek, Denver, 12/2007

use dmumps_struc_def
implicit none

! type of reals
integer,parameter,private  :: kr = kind(1.D0)
! numerical zero
real(kr),parameter,private :: numerical_zero = 1.e-12_kr

! Use this structure of MUMPS for routines from mumps
type(DMUMPS_STRUC),private :: bddc_mumps

! Weight matrix D_P shared among routines in this module
integer,private :: ldp
real(kr),allocatable,private :: dp(:)

! BUF is shared among routines of the module and serves for global MPI communication
! prevents too many allocations
integer,private :: lbuf
real(kr),allocatable,private :: buf(:)

! Setting if dual problem (i.e. Lagrange multipliers) should be used to apply averages
logical,private :: use_dual_problem = .false.
! Setting if projection on nullspace of G should be used to apply averages
logical,private :: use_projection   = .false.
! Setting if change of variables should be used to apply averages
logical,private :: use_transform    = .true.
logical,private :: use_adaptive_averages   = .true.
logical,private :: use_arithmetic_averages = .false.

! Data for DUAL PROBLEM approach
! matrix G with constraints
integer,private              :: lg1 = 0, lg2 = 0
real(kr),allocatable,private :: g(:,:)
! dual matrix for averages
integer,private ::             ldualm1 = 0, ldualm2 = 0
real(kr),allocatable,private :: dualm(:,:)
integer,private ::            lipiv = 0
integer,allocatable,private :: ipiv(:)
! dual RHS for averages
integer,private ::             ldualrhs = 0
real(kr),allocatable,private :: dualrhs(:)

! Data for PROJECTION approach (used also for CHANGE OF VARIABLES)
! matrix F for projections, R^T*F = G where G^T = QR
integer,private              :: nnz_f = 0, lf_sparse, f_type = 0
integer,private              :: nrowf = 0, ncolf = 0
integer,allocatable,private  :: i_f_sparse(:), j_f_sparse(:)
real(kr),allocatable,private :: f_sparse(:)

! Data for CHANGE OF VARIABLES approach
! disk unit for transformation
integer,private ::             idtr = 0
! sparse matrix of transformation
integer,private              :: nnz_t = 0, lt_sparse, t_type = 0
integer,allocatable,private  :: i_t_sparse(:), j_t_sparse(:)
real(kr),allocatable,private :: t_sparse(:)

! Time measurements
integer,parameter,private :: level_time_max = 20
real(kr),private ::          times(level_time_max) = 0._kr
integer,private  ::          level_time = 0
integer,private  ::          time_verbose = 1

contains

!*********************************************************************************************
subroutine bddc_init(myid,comm,matrixtype,mumpsinfo,timeinfo,iterate_on_transformed,&
                     ndoft,nnz,i_sparse,j_sparse,a_sparse,la, &
                     weight_approach,averages_approach,ndim,nglb,inglb,linglb,nnglb,lnnglb,&
                     nnodt,nndft,lnndft,ihntn,lihntn,slavery,lslavery,kdoft,lkdoft, &
                     solt,lsolt, rest,lrest, nnz_transform,idtrmain)
!*********************************************************************************************
! Subroutine for initialization of BDDC preconditioner
! module for using MUMPS package
      use module_mumps
      use module_sm
       
      implicit none
      include "mpif.h"

! Process number within communicator COMM
      integer,intent(in) :: myid

! MPI communicator
      integer,intent(in) :: comm

! Type of matrix
! 0 - unsymmetric                 
! 1 - symmetric positive definite 
! 2 - symmetric general          
      integer,intent(in) :: matrixtype

! What approach should be used for averaging D_P
! 1 - aritmetic averaging (cardinality)
! 2 - weighted averaging - according to diagonal stiffness of original matrix
      integer,intent(in):: weight_approach

! What approach do you want to use to handle matrix G of averages
! 0 - do not apply averages
! 1 - apply averages via Lagrange multipliers G (S_c)^-1 G^T * mu = G * E^T res
! 2 - apply averages via projection,P = I-G^T*(GG^T)^-1*G
! 3 - change of variables P*T^T*A*T*P & projection
      integer,intent(in):: averages_approach

! Level of information from MUMPS
      integer,intent(in) :: mumpsinfo

! Level of messages of timing
! 0 - no output about times
! 1 - basic times printed
! 2 - rigorous timing output for profiling
      integer,intent(in) :: timeinfo

! Iterate on transformed problem?
      logical,intent(in) :: iterate_on_transformed

! Length of global solution in W_tilde
      integer,intent(in) :: ndoft

! Number of nonzeros in matrix
      integer,intent(in) :: nnz

! Matrix in sparse format IJA of local triplets
      integer,intent(in)     :: la
      integer,intent(inout)  :: i_sparse(la), j_sparse(la)
      real(kr),intent(inout) :: a_sparse(la)

! Number of dimensions
      integer,intent(in) :: ndim
! Description of globs
      integer,intent(in) :: nglb
      integer,intent(in) :: linglb,        lnnglb
      integer,intent(in) ::  inglb(linglb), nnglb(lnnglb)

! Description of W_tilde
      integer, intent(in) :: nnodt, lnndft, lihntn, lslavery, lkdoft
      integer, intent(in) :: nndft(lnndft), ihntn(lihntn), slavery(lslavery), kdoft(lkdoft)

! Initial solution and residual
      integer,intent(in)     :: lsolt,       lrest
      real(kr),intent(inout) ::  solt(lsolt), rest(lrest)

! nonzeros from change of variables
      integer,intent(out) :: nnz_transform

! Disk unit for transformation
      integer, optional, intent(in) :: idtrmain

! Local variables
      ! nonzeros from projection
      integer :: nnz_proj      = 0
      integer :: nglbg
      integer :: ierr

      real(kr) :: time

! Optional arguments
      if (present(idtrmain)) then
         idtr = idtrmain
      end if

! setting flags of approach to averages application
      select case (averages_approach)
      case (0)
         use_dual_problem = .false.
         use_projection   = .false.
         use_transform    = .false.
      case (1)
         use_dual_problem = .true.
         use_projection   = .false.
         use_transform    = .false.
      case (2)
         use_dual_problem = .false.
         use_projection   = .true.
         use_transform    = .false.
      case (3)
         use_dual_problem = .false.
         use_projection   = .false.
         use_transform    = .true.
      case default
         write(*,*) 'bddc_init: Unknown type of averages approach. Maybe AVERAGES_APPROACH not set.'
      end select

! Switch off any application of averages if no averages are applied
      call MPI_ALLREDUCE(nglb,nglbg,1,MPI_INTEGER,MPI_SUM,comm,ierr)
      if (nglbg.eq.0) then
         use_dual_problem = .false.
         use_projection   = .false.
         use_transform    = .false.
      end if

! Set time verbose level
      time_verbose = timeinfo

! Initialize MUMPS
      call mumps_init(bddc_mumps,comm,matrixtype)

! Level of information from MUMPS
      call mumps_set_info(bddc_mumps,mumpsinfo)

! Allocate buffer for MPI communication
      lbuf = ndoft
      allocate(buf(lbuf))

! Prepare projection onto null G
      if (use_projection) then
         call bddc_time_start(comm)
         call bddc_P_init(myid,comm,ndim,nglb,inglb,linglb,nnglb,lnnglb,slavery,lslavery,nnodt,ndoft,kdoft,lkdoft,&
                          matrixtype,nnz,i_sparse,j_sparse,a_sparse,la, nnz_proj)
         call bddc_time_end(comm,time)
         if (myid.eq.0.and.time_verbose.ge.1) then
            write(*,*) '==========================================='
            write(*,*) 'Time of construction of projection = ',time
            write(*,*) '==========================================='
            call flush(6)
         end if
      end if

! Prepare transformation of coordinates
      nnz_transform = 0
      if (use_transform) then
         call bddc_time_start(comm)
         call bddc_T_init(myid,comm,ndim,nglb,inglb,linglb,nnglb,lnnglb,slavery,lslavery,nnodt,ndoft,kdoft,lkdoft,&
                          iterate_on_transformed,matrixtype,nnz,i_sparse,j_sparse,a_sparse,la, nnz_transform, &
                          solt,lsolt, nnz_proj)
         call bddc_time_end(comm,time)
         if (myid.eq.0.and.time_verbose.ge.1) then
            write(*,*) '==========================================='
            write(*,*) 'Time of transforming the matrix = ',time
            write(*,*) '==========================================='
            call flush(6)
         end if
      end if

! Load matrix to MUMPS
!      ltestm1 = ndoft
!      ltestm2 = ndoft
!      allocate(testm(ltestm1,ltestm2))
!      call sm_to_dm(matrixtype,i_sparse, j_sparse, a_sparse,nnz+nnz_transform+nnz_proj, testm,ltestm1,ltestm2)
!      write(98,'(i8)') ndoft
!      do i = 1,ltestm1
!         write(98,'(1000f10.4)') (testm(i,j),j = 1,ltestm2)
!      end do
!      deallocate(testm)

      call mumps_load_triplet(bddc_mumps,ndoft,nnz+nnz_transform+nnz_proj,i_sparse,j_sparse,a_sparse,la)

! Analyze matrix
      call bddc_time_start(comm)
      call mumps_analyze(bddc_mumps) 
      call bddc_time_end(comm,time)
      if (myid.eq.0.and.time_verbose.ge.1) then
         write(*,*) '==========================================='
         write(*,*) 'Time of MUMPS analysis = ',time
         write(*,*) '==========================================='
         call flush(6)
      end if

! Factorize matrix
      call bddc_time_start(comm)
      call mumps_factorize(bddc_mumps)
      call bddc_time_end(comm,time)
      if (myid.eq.0.and.time_verbose.ge.1) then
         write(*,*) '==========================================='
         write(*,*) 'Time of MUMPS factorization = ',time
         write(*,*) '==========================================='
         call flush(6)
      end if

! Create weigth matrix DP
      ldp = ndoft
      allocate(dp(ldp))
      select case (weight_approach)
      case (1)
         ! according to cardinality
         call bddc_DP_build_card(comm,nnodt,nndft,lnndft,slavery,lslavery,kdoft,lkdoft,dp,ldp)
      case (2)
         ! according to diagonal stiffness
         call bddc_DP_build_stiff(comm,nnz,i_sparse,j_sparse,a_sparse,la, &
                                  nnodt,nndft,lnndft,slavery,lslavery,kdoft,lkdoft,dp,ldp)
      case default
         if (myid.eq.0) then
            write(*,*) 'bddc_init: Illegal value of WEIGHT_APPROACH:', weight_approach
         end if
         stop
      end select

! Correct initial solution and residual
      if (iterate_on_transformed) then
         solt = dp*solt
         call MPI_ALLREDUCE(solt,buf,lsolt,MPI_DOUBLE_PRECISION,MPI_SUM,comm,ierr)
         solt = buf
         if (use_transform) then
            call bddc_time_start(comm)
            call bddc_T_apply(comm,rest,lrest,transposed=.true.)
            call bddc_time_end(comm,time)
            if (myid.eq.0.and.time_verbose.ge.2) then
               write(*,*) '=============================================='
               write(*,*) 'Time of application of initial transformation = ',time
               write(*,*) '=============================================='
               call flush(6)
            end if
         end if
      end if

! Prepare dual problem for solution using Lagrange multipliers
      if (use_dual_problem) then
         ! Prepare space for matrix G
         call bddc_G_aritmetic_size(nglb,ndim,ndoft,inglb,linglb,nnglb,lnnglb,ihntn,lihntn,slavery,lslavery, &
                                    lg1,lg2)
         if (myid.eq.0) then
            write(*,*) 'number of rows in G is', lg1
         end if
         allocate(g(lg1,lg2))
         ! Build G
         call bddc_G_aritmetic_build(nglb,ndim,nnodt,inglb,linglb,nnglb,lnnglb,ihntn,lihntn,slavery,lslavery, &
                                     kdoft,lkdoft, g,lg1,lg2)
         call bddc_dual_init(myid,g,lg1,lg2)
      end if

      return
end subroutine

!******************************************************************************************
subroutine bddc_DP_build_card(comm,nnodt,nndft,lnndft,slavery,lslavery,kdoft,lkdoft,dp,ldp)
!******************************************************************************************
! Subroutine for creation of weight matrix D_P based on CARDINALITY, 
! i.e. number of subdomains sharing the node
       
      implicit none

! MPI communicator
      integer,intent(in) :: comm

! Description of W_tilde
      integer, intent(in) :: nnodt
      integer, intent(in) :: lnndft, lslavery, lkdoft
      integer, intent(in) :: nndft(lnndft), slavery(lslavery), kdoft(lkdoft)

! Weight matrix
      integer, intent(in) :: ldp
      real(kr),intent(out)::  dp(ldp)

      dp = 1.0_kr
      call bddc_RRT(comm,nnodt,nndft,lnndft,slavery,lslavery,kdoft,lkdoft,dp,ldp)
      where(dp.ne.1.0_kr) dp = 1.0_kr/dp

      return
end subroutine

!**************************************************************************************
subroutine bddc_DP_build_stiff(comm,nnz,i_sparse,j_sparse,a_sparse,la, &
                               nnodt,nndft,lnndft,slavery,lslavery,kdoft,lkdoft,dp,ldp)
!**************************************************************************************
! Subroutine for creation of weight matrix D_P based on DIAGONAL STIFFNESS OF MATRIX
       
      implicit none

! MPI communicator
      integer,intent(in) :: comm

! Number of nonzeros in matrix
      integer,intent(in) :: nnz

! Matrix in sparse format AIJ of local triplets
      integer,intent(in)  :: la
      integer,intent(in)  :: i_sparse(la), j_sparse(la)
      real(kr),intent(in) :: a_sparse(la)

! Description of W_tilde
      integer, intent(in) :: nnodt, lnndft, lslavery, lkdoft
      integer, intent(in) :: nndft(lnndft), slavery(lslavery), kdoft(lkdoft)

! Weight matrix
      integer, intent(in) :: ldp
      real(kr),intent(out)::  dp(ldp)

! Local variables
      integer  :: ia, indi, indj
      integer  :: laux
      real(kr),allocatable :: aux(:)

      dp = 0._kr
      do ia = 1,nnz
         indi = i_sparse(ia)
         indj = j_sparse(ia)
         if (indi.eq.indj) then
            dp(indi) = a_sparse(ia)
         end if
      end do
      laux = ldp
      allocate(aux(laux))
      aux = dp
      call bddc_RRT(comm,nnodt,nndft,lnndft,slavery,lslavery,kdoft,lkdoft,dp,ldp)
      dp = aux/dp
      deallocate(aux)

      return
end subroutine

!**************************************************************************************************
subroutine bddc_M(myid,comm,iterate_on_transformed,nnodt,nndft,lnndft,slavery,lslavery,kdoft,lkdoft,&
                  rest,lrest,ht,lht,rmr)
!**************************************************************************************************
! Subroutine for realization of BDDC preconditioner M_BDDC
! module for using MUMPS
      use module_mumps

      implicit none
      include "mpif.h"

! MPI rank
      integer, intent(in) :: myid

! MPI communicator
      integer, intent(in) :: comm

! Iterate on transformed problem?
      logical,intent(in) :: iterate_on_transformed

! Description of space W_tilde
      integer, intent(in) :: nnodt, lnndft, lslavery, lkdoft
      integer, intent(in) :: nndft(lnndft), slavery(lslavery), kdoft(lkdoft)

! Residual
      integer,  intent(in) :: lrest
      real(kr), intent(in) :: rest(lrest)
! Preconditioned residual
      integer,  intent(in)  :: lht
      real(kr), intent(out) :: ht(lht)

! rMr parameter in PCG
      real(kr), intent(out) :: rmr

! Local variables
      real(kr) :: rmr_loc
      integer ::  ierr

      real(kr) :: time

! Make copy of residual to initial ht
      ht = rest

! Apply weight matrix to residual
! D_P * ht => ht
      ht = dp * ht

! Reduce weighted residual
!***************************************************************************MPI
      if (lbuf.ne.ldp.or..not.allocated(buf)) then
         write(*,*) 'bddc_M: MPI Buffer not ready. Maybe missing BDDC_INIT call.'
         stop
      end if
      call MPI_ALLREDUCE(ht,buf,lht,MPI_DOUBLE_PRECISION,MPI_SUM,comm,ierr)
      ht = buf
!***************************************************************************MPI

! Add correction from Lagrange multipliers
! ht - G^T*(G*(A^c)^-1*G^T)^-1*G*(A^c)^-1*ht => ht
      if (use_dual_problem) then
         call bddc_time_start(comm)
         ! assumes same values in buf and ht
         call bddc_dual_solve(myid,g,lg1,lg2,buf,lbuf, ht,lht)
         call bddc_time_end(comm,time)
         if (myid.eq.0.and.time_verbose.ge.2) then
            write(*,*) '=============================================='
            write(*,*) 'Time of dual problem solution = ',time
            write(*,*) '=============================================='
            call flush(6)
         end if
      end if

! Transform residual to new variables
! T^T * ht => ht
      if (.not.iterate_on_transformed) then
         if (use_transform) then
            call bddc_time_start(comm)
            call bddc_T_apply(comm,ht,lht,transposed=.true.)
            call bddc_time_end(comm,time)
            if (myid.eq.0.and.time_verbose.ge.2) then
               write(*,*) '=============================================='
               write(*,*) 'Time of application of transformation on residual = ',time
               write(*,*) '=============================================='
               call flush(6)
            end if
         end if
      end if

! Project residual onto null G
! P*ht => ht, where P = I - G^T*(G*G^T)^-1*G = I - F^T*F, R^T*F = G, G^T = Q*R
      if (use_projection.or.use_transform) then
         ! matrix F is different in both cases
         call bddc_time_start(comm)
         call bddc_P_apply(comm,ht,lht)
         call bddc_time_end(comm,time)
         if (myid.eq.0.and.time_verbose.ge.2) then
            write(*,*) '==========================================='
            write(*,*) 'Time of application of projection = ',time
            write(*,*) '==========================================='
            call flush(6)
         end if
      end if

! Solve the system with weighted residual as right hand side
! Atilde^-1 * ht => ht
      call bddc_time_start(comm)
      call mumps_resolve(bddc_mumps,ht,lht)
      ! distribute solution to all processors
!***************************************************************************MPI
      call MPI_BCAST(ht,lht, MPI_DOUBLE_PRECISION, 0, comm, ierr)
!***************************************************************************MPI
      call bddc_time_end(comm,time)
      if (myid.eq.0.and.time_verbose.ge.2) then
         write(*,*) '==========================================='
         write(*,*) 'Time of MUMPS backsubstitution = ',time
         write(*,*) '==========================================='
      end if

! Transform preconditioned residual into old variables
! T*ht => ht
      if (.not.iterate_on_transformed) then
         if (use_transform) then
            call bddc_time_start(comm)
            call bddc_T_apply(comm,ht,lht)
            call bddc_time_end(comm,time)
            if (myid.eq.0.and.time_verbose.ge.2) then
               write(*,*) '=============================================='
               write(*,*) 'Time of application of transformation on residual = ',time
               write(*,*) '=============================================='
            end if
         end if
      end if

! Apply weight matrix
! D_P * ht => ht
      ht = dp * ht

! determine the coefficient rMr
! rMr = rest * D_P * T * Atilde^-1 * P * T^T * D_P * rest
      rmr_loc = dot_product(rest,ht)
!***************************************************************************MPI
      call MPI_ALLREDUCE(rmr_loc,rmr,1,MPI_DOUBLE_PRECISION,MPI_SUM,comm,ierr)
!***************************************************************************MPI
 
! apply operator RR^T to preconditioned residual
! R * R^T * ht => ht
      call bddc_RRT(comm,nnodt,nndft,lnndft,slavery,lslavery,kdoft,lkdoft,ht,lht)
      
      return
end subroutine

!**************************************************************************************************
subroutine bddc_M_fake(comm,rest,lrest,ht,lht,rmr)
!**************************************************************************************************
! Subroutine for void preconditioning - multiplying by identity
! module for using MUMPS

      implicit none
      include "mpif.h"

! MPI communicator
      integer, intent(in) :: comm

! Residual
      integer,  intent(in) :: lrest
      real(kr), intent(in) :: rest(lrest)
! Preconditioned residual
      integer,  intent(in)  :: lht
      real(kr), intent(out) :: ht(lht)

! rMr parameter in PCG
      real(kr), intent(out) :: rmr

! Local variables
      real(kr) :: rmr_loc
      integer ::  ierr

! Make copy of residual to initial ht
      ht = rest

! determine the coefficient rMr
! rMr = rest * D_P * T * Atilde^-1 * P * T^T * D_P * rest
      rmr_loc = bddc_dot_product(comm,rest,lrest,ht,lht)
!***************************************************************************MPI
      call MPI_ALLREDUCE(rmr_loc,rmr,1,MPI_DOUBLE_PRECISION,MPI_SUM,comm,ierr)
!***************************************************************************MPI
 
      return
end subroutine

!****************************************
subroutine bddc_dual_init(myid,g,lg1,lg2)
!****************************************
! Subroutine for creation of dual problem from matrix G - problem for Lagrange multipliers
! G (Ac)^-1 G^T and its factorization by LAPACK
! For using this routine, G is constructed for ALL globs and is same on all processors.

! Use module for MUMPS package
      use module_mumps
       
      implicit none

! Process number
      integer,intent(in) :: myid

! Matrix G of averages
      integer,intent(in)  :: lg1, lg2
      real(kr),intent(in) :: g(lg1,lg2)

! Local variables
      integer  :: lh1, lh2
      real(kr),allocatable :: h(:,:)

      integer :: nrhs, lapack_info, ldim, lh

! prepare space for H where Ac*H = G^T
      lh1 = lg2
      lh2 = lg1
      allocate(h(lh1,lh2)) 
      h  = transpose(g)
      lh = lh1*lh2
      nrhs = lh2
      call mumps_resolve(bddc_mumps,h,lh,nrhs)

! find G*H => dualm
      ldualm1 = lg1
      ldualm2 = lg1
      allocate(dualm(ldualm1,ldualm2))
      if (myid.eq.0) then
         dualm = matmul(g,h)
      end if
      deallocate(h)

! Factorize dual problem matrix DUALM with LAPACK routine
      ldim  = max(1,ldualm1)
      lipiv = ldualm1
      allocate(ipiv(lipiv))
      if (myid.eq.0) then
         call DGETF2(ldualm1, ldualm2, dualm, ldim, ipiv, lapack_info)
         if (lapack_info.ne.0) then
            write(*,*) 'Error in Lapack factorization of dual problem'
            stop
         end if
         write(*,*) 'Dual matrix factorized'
      end if

! Prepare dual RHS for backsubstitutions
      ldualrhs = ldualm1
      allocate(dualrhs(ldualrhs))

      return
end subroutine

!**********************************************************
subroutine bddc_dual_solve(myid,g,lg1,lg2,aux,laux, ht,lht)
!**********************************************************
! Subroutine for addition of correction from dual problem - Lagrange multipliers
! ht - G^T*(G*(A^c)^-1*G^T)^-1*G*(A^c)^-1*ht => ht

! module for using MUMPS
      use module_mumps

      implicit none

! MPI rank
      integer, intent(in) :: myid

! Matrix G of averages
      integer,intent(in)  :: lg1, lg2
      real(kr),intent(in) :: g(lg1,lg2)

! initial value of ht
      integer,  intent(in) :: laux
      real(kr), intent(in) :: aux(laux)

! preconditioned residual
      integer,  intent(in)    :: lht
      real(kr), intent(inout) :: ht(lht)

! Local variables
      integer :: lapack_info, ldim

! Solve the system with weighted residual as right hand side on master processor
! Ac^-1 * ht => ht
      call mumps_resolve(bddc_mumps,ht,lht)
      
! Only root can continue
      if (myid.eq.0) then
         ! construct dual RHS 
         ! G*ht => dualrhs
         if (.not.allocated(dualrhs)) then
            write(*,*) 'bddc_dual_solve: DUALRHS not ready, maybe missing call for bddc_init.'
            stop
         end if
         dualrhs = matmul(g,ht)
         ! solve by LAPACK 
         ! (G*(Ac)^-1*G^T)^-1*dualrhs => dualrhs
         ldim = max(1,ldualm1)
         call DGETRS( 'N', ldualm1, 1, dualm, ldim, ipiv, dualrhs, ldim, lapack_info)
         if (lapack_info.ne.0) then
            write(*,*) 'bddc_dual_solve: Error in Lapack factorization of dual problem'
            stop
         end if
         ! add the correction from dual problem 
         ! ht - G^T * dualrhs => ht
         ht = aux - matmul(transpose(g),dualrhs)
      end if
      
      return
end subroutine

!****************************************************************************************
subroutine bddc_P_nnz_est(matrixtype,sparsity,ndofs,ndim,nglb,inglb,linglb,nnglb,lnnglb,&
                          slavery,lslavery,nnz_proj_est)
!****************************************************************************************
! Subroutine for estimation of space necessary for NEW entries of the projection onto null G
       
      implicit none

! Type of matrix
! 0 - unsymmetric                 
! 1 - symmetric positive definite 
! 2 - symmetric general          
      integer,intent(in) :: matrixtype

! Sparsity of the matrix
      real(kr),intent(in) :: sparsity

! Number of variables on subdomain
      integer,intent(in) :: ndofs

! Number of dimensions
      integer,intent(in) :: ndim
! Description of globs
      integer,intent(in) :: nglb
      integer,intent(in) :: linglb,        lnnglb
      integer,intent(in) ::  inglb(linglb), nnglb(lnnglb)

! Description of spcae W_tilde
      integer,intent(in) :: lslavery      
      integer,intent(in) ::  slavery(lslavery)

! Number of values added into matrix
      integer,intent(out) :: nnz_proj_est

! local variables
      integer :: koef, indinglb, nglbn, indnglb, indnglb_m, nsubglb, ndofn, iglb
      integer :: block1, block2, ncols, nrows
      real(kr) :: sparsity1D

! Estimate 1D sparsity
      sparsity1D = sqrt(sparsity)

      indinglb = 0
      ncols    = 0
      nrows    = 0
      do iglb = 1,nglb
         nglbn = nnglb(iglb)
         indnglb   = inglb(indinglb + 1)
         indnglb_m = slavery(indnglb)
         nsubglb   = count(slavery.eq.indnglb_m)
         ndofn     = ndim

         ncols = ncols + nsubglb * ndofn * nglbn
         nrows = nrows + ndofn * (nsubglb-1)

         indinglb = indinglb + nglbn
      end do

      ! nonzero blocks
      block1 = nrows * sparsity1D * ndofs * ncols
      block2 = ncols * ncols

      select case (matrixtype)
         case(0)
            koef = 2
         case(1,2)
            koef = 1.1
         case default
            write(*,*) 'Unknown type of matrix. Maybe MATRIXTYPE not set.'
            stop
      end select
      nnz_proj_est = koef*block1 + koef*block2

      return
end subroutine

!***************************************************************************************************************
subroutine bddc_P_init(myid,comm,ndim,nglb,inglb,linglb,nnglb,lnnglb,slavery,lslavery,nnodt,ndoft,kdoft,lkdoft,&
                       matrixtype,nnz,i_sparse,j_sparse,a_sparse,la, nnz_proj)
!***************************************************************************************************************
! Subroutine for setting up the projection onto null G
!!!!! Works only for SYMMETRIC matrices!!!!!
! A_avg = P*A*P + t(I-P), where P = I - G^T*(G*G^T)^-1*G = I - F^T*F, R^T*F = G, G^T = Q*R, t stabilization parameter
! Storage in A after projection:
! A | t*F^T*F | -A*F^T*F | -F^T*F*A | F^T*F*A*F^T*F |
!nnz|                   nnz_proj                    |
! For symmetric case, block -F^T*F*A  is obtained by flipping block -A*F^T*F
! along diagonal and double entries on diagonal.
! Blocks of projection are assemblaged.

! Use sparse matrices module
      use module_sm
       
      implicit none
      include "mpif.h"

! Process number
      integer,intent(in) :: myid

! MPI communicator
      integer, intent(in) :: comm

! Number of dimensions
      integer,intent(in) :: ndim
! Description of globs
      integer,intent(in) :: nglb
      integer,intent(in) :: linglb,        lnnglb
      integer,intent(in) ::  inglb(linglb), nnglb(lnnglb)

! Description of spcae W_tilde
      integer,intent(in) :: nnodt, ndoft
      integer,intent(in) :: lslavery,           lkdoft
      integer,intent(in) ::  slavery(lslavery),  kdoft(lkdoft)

! Type of matrix
! 0 - unsymmetric                 
! 1 - symmetric positive definite 
! 2 - symmetric general          
      integer,intent(in) :: matrixtype

! Matrix in sparse format AIJ of local triplets
      integer,intent(in)     :: nnz, la
      integer,intent(inout)  :: i_sparse(la), j_sparse(la)
      real(kr),intent(inout) :: a_sparse(la)

! Number of values added into matrix
      integer,intent(out) :: nnz_proj

! Local variables
      integer  :: lg1, lg2
      real(kr),allocatable :: g(:,:)
      integer  :: lr1, lr2
      real(kr),allocatable :: r(:,:)
      integer :: ltau
      real(kr),allocatable :: tau(:)
      integer :: lwork
      real(kr),allocatable :: work(:)
      integer ::             lf1 = 0, lf2 = 0
      real(kr),allocatable :: f(:,:)
      integer  :: lh1, lh2
      real(kr),allocatable :: h(:,:)
      integer  :: liglbntn
      integer,allocatable :: iglbntn(:)
      integer  :: laux1, laux2
      real(kr),allocatable :: aux(:,:)

      integer :: iconstr, idofn, lapack_info = 0, indnglb_m, &
                 indnglb_s, iglb, indfline, indnglb, jfline, inodglb, ldim, &
                 nconstr_loc, nrhs, ndofn, nglbn, nsubglb, &
                 indinglb, inodt, ig, jg, igtg, jgtg, jnodt, ipoint, jpoint, &
                 point_a, start_j, i, store_type
      real(kr):: value

      integer :: ia, ierr, nnz_add, nnz_new, point, space_left, iaux
      real(kr):: t_stab_loc, t_stab, scalar

      real(kr):: time

      logical :: i_am_master_of_this_glob

! Check the applicability of the routine - symmetric matrices
      if (matrixtype.ne.1.and.matrixtype.ne.2) then
         if (myid.eq.0) then
            write(*,*) 'bddc_P_init: This routine currently works only with symmetric matrices.'
         end if
         stop
      end if

! count number of lines in F
      indinglb = 0
      nconstr_loc = 0
      do iglb = 1,nglb
         nglbn = nnglb(iglb)
         indnglb   = inglb(indinglb + 1)
         indnglb_m = slavery(indnglb)
         nsubglb   = count(slavery.eq.indnglb_m)
         ndofn     = ndim

         nconstr_loc = nconstr_loc + ndofn*(nsubglb-1)

         indinglb = indinglb + nglbn
      end do

! find stabilization factor t for projection
      t_stab_loc = 0
      write(*,*) 'nnz in proj start= ',nnz
      do ia = 1,nnz
         if (i_sparse(ia).eq.j_sparse(ia)) then
            t_stab_loc = max(abs(a_sparse(ia)),t_stab_loc)
         end if
      end do
      call MPI_ALLREDUCE(t_stab_loc,t_stab,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr)
      if (myid.eq.0) then
         write(*,*) 'Stabilization parameter in projection t =',t_stab
      end if

! Construct F in globs on the subdomain
      call bddc_time_start(comm)

      lf1 = nconstr_loc
      lf2 = ndoft
      allocate(f(lf1,lf2))
      f = 0._kr
      indinglb = 0
      indfline = 0
      nnz_new  = 0
      do iglb = 1,nglb
         nglbn = nnglb(iglb)
         indnglb   = inglb(indinglb + 1)
         indnglb_m = slavery(indnglb)
         if (indnglb_m.eq.indnglb) then
            i_am_master_of_this_glob = .true.
         else
            i_am_master_of_this_glob = .false.
         end if

         nsubglb   = count(slavery.eq.indnglb_m)
         ndofn     = ndim

         ! Find mapping of glob unknows into W_tilde
         liglbntn = nsubglb * nglbn
         allocate(iglbntn(liglbntn))
         do inodglb = 1,nglbn
            indnglb   = inglb(indinglb + inodglb)
            indnglb_m = slavery(indnglb)
            iglbntn(inodglb) = indnglb_m
            indnglb_s =  0
            ! loop over columns of G
            do iconstr = 1,nsubglb-1
               ! find slaves
               do inodt = indnglb_s + 1 ,nnodt
                  if (slavery(inodt).eq.indnglb_m.and.inodt.ne.indnglb_m) then
                     indnglb_s = inodt
                     exit
                  end if
               end do
               if (indnglb_s.eq.0) then
                  write(*,*) 'bddc_P_init: Error in searching slave index.'
                  call flush(6)
                  stop
               end if
               iglbntn(iconstr*nglbn + inodglb) = indnglb_s
            end do
         end do

         ! prepare matrix G - local for glob
         lg1       = nsubglb - 1
         lg2       = nsubglb * nglbn
         allocate(g(lg1,lg2))
         g = 0._kr
         
         do iconstr = 1,nsubglb-1
            g(iconstr,1:nglbn)                           =  1._kr
            g(iconstr,iconstr*nglbn+1:(iconstr+1)*nglbn) = -1._kr
         end do
!         write(*,*) 'Matrix G'
!         do i = 1,lg1
!            write(*,'(30f10.5)') g(i,:)
!         end do
!         call flush(6)

!        QR decomposition of matrix G^T
         lr1 = lg2
         lr2 = lg1
         allocate(r(lr1,lr2))
         r = transpose(g)
         ! LAPACK arrays
         ltau = lr2
         allocate(tau(ltau))
         lwork = lr2
         allocate(work(lwork))
         ! leading dimension
         ldim = max(1,lr1)
         call DGEQR2( lr1, lr2, r, ldim, tau, work, lapack_info)
         if (lapack_info.ne.0) then
            write(*,*) 'bddc_P_init: Error in LAPACK QR factorization of G: ', lapack_info
            call flush(6)
            stop
         end if
         deallocate(tau,work)
!         write(*,*) 'Matrix R'
!         do i = 1,lr1
!            write(*,'(30f10.5)') r(i,:)
!         end do
!         call flush(6)

!        find local block of F, such that R^T F = G, store F in G
         ldim = max(1,lr2)
         nrhs = lg2
         call DTRTRS( 'L',   'N',   'N', lr2 , nrhs, transpose(r), ldim, g, ldim, lapack_info)
         if (lapack_info.ne.0) then
            write(*,*) 'Error in LAPACK solution of R^T * F = G :', lapack_info
            call flush(6)
            stop
         end if
         deallocate(r)

         ! copy block of F stored in G into global F
         do jg = 1,lg2
            indnglb = iglbntn(jg)
            point   = kdoft(indnglb)
            ! loop over rows of G
            do ig = 1,lg1
               ! loop over diagonal of degrees of freedom at the node
               jfline = (ig-1)*ndofn
               do idofn = 1,ndofn
                  f(indfline + jfline + idofn,point + idofn) = g(ig,jg)
               end do
            end do
         end do

         ! add entries from t_stab*F^T*F - for this glob
         if (i_am_master_of_this_glob) then
            point_a = nnz + nnz_new
            ia = 0
            do igtg = 1,lg2
               inodt = iglbntn(igtg)
               ipoint = kdoft(inodt)
               select case (matrixtype)
                  case(0)
                     start_j = 1
                  case(1,2)
                     start_j = igtg
                  case default
                     write(*,*) 'bddc_P_init: Unknown value of MATRIXTYPE.'
                     stop
               end select
               do jgtg = start_j,lg2
                  jnodt = iglbntn(jgtg)
                  jpoint = kdoft(jnodt)

                  value = dot_product(g(:,igtg),g(:,jgtg))
                  do idofn = 1,ndofn
                     ia = ia + 1
                     i_sparse(point_a + ia) = ipoint + idofn
                     j_sparse(point_a + ia) = jpoint + idofn
                     a_sparse(point_a + ia) = t_stab*value
                  end do
               end do
            end do
            nnz_new = nnz_new + ia
         end if

         deallocate(g)
         deallocate(iglbntn)
         indfline = indfline + ndofn*(nsubglb-1)
         indinglb = indinglb + nglbn
      end do
!      write(*,*) 'myid:',myid,'Matrix F'
!      do i = 1,lf1
!         write(*,'(30f10.5)') f(i,:)
!      end do
!      call flush(6)
      call bddc_time_end(comm,time)
      if (myid.eq.0.and.time_verbose.ge.2) then
         write(*,*) '==========================================='
         write(*,*) 'Time of F construction = ',time
         write(*,*) '==========================================='
      end if


      call bddc_time_start(comm)
      ! find Ac * F^T => H
      lh1 = lf2
      lh2 = lf1
      allocate(h(lh1,lh2))
      call sm_mat_mult(matrixtype, nnz, i_sparse, j_sparse, a_sparse, la, &
                       transpose(f),lf2,lf1, h,lh1,lh2)
      call bddc_time_end(comm,time)
      if (myid.eq.0.and.time_verbose.ge.2) then
         write(*,*) '==========================================='
         write(*,*) 'Time of H construction = ',time
         write(*,*) '==========================================='
      end if
!      write(*,*) 'Matrix H'
!      do i = 1,lh1
!         write(*,'(30f10.5)') h(i,:)
!      end do
!      call flush(6)
      
! Convert F to sparse matrix
      lf_sparse = count(f.ne.0._kr) 
      write(*,*) 'myid =',myid,': F allocated to nnz_f =',lf_sparse
      nrowf     = lf1
      ncolf     = lf2
      allocate(i_f_sparse(lf_sparse),j_f_sparse(lf_sparse),f_sparse(lf_sparse))
      call sm_from_dm(f_type, f,lf1,lf2, i_f_sparse, j_f_sparse, f_sparse, lf_sparse, nnz_f)
      ! clear memory
      deallocate(f)
      write(*,*) 'myid =',myid,': F converted to sparse, nnz_f =',nnz_f
      call flush(6)

! Add entries from projection to the matrix
      call bddc_time_start(comm)

      ! add entries from -Ac*F^T*F - F^T*F*Ac = -H*F - (H*F)^T
      call bddc_time_start(comm)
      point      = nnz + nnz_new + 1
      space_left = la - nnz - nnz_new
      scalar     = -1._kr
      store_type = 0
!     F^T * H^T
      call sm_from_sm_mat_mult(f_type, nnz_f, j_f_sparse, i_f_sparse, f_sparse, lf_sparse, &
                               transpose(h),lh2,lh1, store_type,ndoft, &
                               i_sparse(point),j_sparse(point),a_sparse(point),space_left, nnz_add, scalar)
      ! double entries on diagonal
      where (i_sparse(point:point+nnz_add).eq.j_sparse(point:point+nnz_add)) &
            a_sparse(point:point+nnz_add) = 2._kr*a_sparse(point:point+nnz_add)
      ! reverse the lower triangle into upper and double entries on diagonal
      do i = point, point+nnz_add
         if (i_sparse(i).gt.j_sparse(i)) then
            iaux = i_sparse(i)
            i_sparse(i) = j_sparse(i)
            j_sparse(i) = iaux
         end if
      end do
      nnz_new    = nnz_new + nnz_add
      call bddc_time_end(comm,time)
      if (myid.eq.0.and.time_verbose.ge.2) then
         write(*,*) '=============================================='
         write(*,*) 'Time of 1st matrix-matrix multiplication = ',time
         write(*,*) '=============================================='
      end if

      ! add entries from F^T*F*Ac*F^T*F = F^T*((F*H)*F)
      call bddc_time_start(comm)
      laux1 = nrowf
      laux2 = nrowf
      allocate(aux(laux1,laux2))
      ! Multiply F*H => aux
      call sm_mat_mult(f_type, nnz_f, i_f_sparse, j_f_sparse, f_sparse, lf_sparse, &
                       h,lh1,lh2, aux,laux1,laux2) 
      ! Multiply (F*H)*F = [F^T*(F*H)^T]^T => H^T
      call sm_mat_mult(f_type, nnz_f, j_f_sparse, i_f_sparse, f_sparse, lf_sparse, &
                       transpose(aux),laux2,laux1, h,lh1,lh2)
      deallocate(aux)
      call bddc_time_end(comm,time)
      if (myid.eq.0.and.time_verbose.ge.2) then
         write(*,*) '=============================================='
         write(*,*) 'Time of F*H multiplications = ',time
         write(*,*) '=============================================='
      end if

      call bddc_time_start(comm)
      point      = nnz + nnz_new + 1
      space_left = la - nnz - nnz_new
      scalar     = 1._kr
      store_type = 1
      ! F^T * AUX2
      call sm_from_sm_mat_mult(f_type, nnz_f, j_f_sparse, i_f_sparse, f_sparse, lf_sparse, &
                               transpose(h),lh2,lh1, store_type, ndoft, &
                               j_sparse(point),i_sparse(point),a_sparse(point),space_left, nnz_add, scalar)
      nnz_new    = nnz_new + nnz_add
      call bddc_time_end(comm,time)
      if (myid.eq.0.and.time_verbose.ge.2) then
         write(*,*) '=============================================='
         write(*,*) 'Time of 3rd matrix-matrix multiplication = ',time
         write(*,*) '=============================================='
      end if
      deallocate(h)

      call bddc_time_end(comm,time)
      if (myid.eq.0.and.time_verbose.ge.2) then
         write(*,*) '=============================================='
         write(*,*) 'Time of all matrix-matrix multiplications = ',time
         write(*,*) '=============================================='
      end if

! Assembly new entries in matrix
      write(*,*) 'myid =',myid,': Space really needed for projection =',nnz_new
      call flush(6)
      call sm_assembly(i_sparse(nnz+1),j_sparse(nnz+1),a_sparse(nnz+1),nnz_new, nnz_proj)
      write(*,*) 'myid =',myid,': Number of nonzeros from projection nnz_proj =',nnz_proj
      write(*,*) 'myid =',myid,': Local fill-in by projection nnz_proj/nnz =',float(nnz_proj)/nnz
      call flush(6)

      return
end subroutine

!***********************************
subroutine bddc_P_apply(comm,ht,lht)
!***********************************
! Subroutine for projection of residual onto null G
! Module for using sparse matrices
      use module_sm

      implicit none
      include "mpif.h"

! MPI communicator
      integer, intent(in) :: comm

! preconditioned residual
      integer,  intent(in)    :: lht
      real(kr), intent(inout) :: ht(lht)

! Local variables
      integer :: ierr

      integer ::             laux1
      real(kr),allocatable :: aux1(:)
      integer ::             laux2
      real(kr),allocatable :: aux2(:)

! construct  F*ht => aux1
      laux1 = nrowf
      allocate(aux1(laux1))
      call sm_vec_mult(f_type, nnz_f, i_f_sparse, j_f_sparse, f_sparse, lf_sparse, &
                       ht,lht, aux1,laux1)

! construct  F^T*aux1 => aux2
      laux2 = ncolf
      allocate(aux2(laux2))
      call sm_vec_mult(f_type, nnz_f, j_f_sparse, i_f_sparse, f_sparse, lf_sparse, &
                       aux1,laux1, aux2,laux2)
      deallocate(aux1)

! apply projection
! ht - F^T*F*ht => ht
      call MPI_ALLREDUCE(aux2,buf,lbuf,MPI_DOUBLE_PRECISION,MPI_SUM,comm,ierr)
      ht = ht - buf

      deallocate(aux2)
      
      return
end subroutine

!*********************************************************************************************
subroutine bddc_T_nnz_est(myid,matrixtype,sparsity,ndofs,ndim,nglb,inglb,linglb,nnglb,lnnglb,&
                          slavery,lslavery, nnz_tr_proj_est)
!*********************************************************************************************
! Subroutine for estimation of space necessary for entries from the change of variables.
       
      implicit none

! Process number
      integer,intent(in) :: myid

! Type of matrix
! 0 - unsymmetric                 
! 1 - symmetric positive definite 
! 2 - symmetric general          
      integer,intent(in) :: matrixtype

! Sparsity of the matrix
      real(kr),intent(in) :: sparsity

! Number of unknowns on subdomain
      integer,intent(in) :: ndofs

! Number of dimensions
      integer,intent(in) :: ndim
! Description of globs
      integer,intent(in) :: nglb
      integer,intent(in) :: linglb,        lnnglb
      integer,intent(in) ::  inglb(linglb), nnglb(lnnglb)

! Description of spcae W_tilde
      integer,intent(in) :: lslavery
      integer,intent(in) ::  slavery(lslavery)

! Number of values added into matrix
      integer,intent(out) :: nnz_tr_proj_est

! local variables
      integer :: koef, block, indinglb, nglbn, nglbv, indnglb, ndofn, iglb, nglbtr, &
                 indinglbtr, iavgn, navgn
      integer :: nnz_transform_est, nnz_proj_est
      integer ::             linglbtr
      integer,allocatable ::  inglbtr(:)
      integer ::             lnnglbtr
      integer,allocatable ::  nnglbtr(:)
      real(kr) :: sparsity1D

! Estimate coefficient of 1D sparsity
      sparsity1D = sqrt(sparsity)

      indinglb = 0
      nnz_transform_est = 0
      nglbtr   = 0
      do iglb = 1,nglb
         nglbn   = nnglb(iglb)
         indnglb = inglb(indinglb + 1)
         ndofn   = ndim

         ! for aritmetic averages - each glob leads to one single node constraint
         navgn   = 1
         nglbtr  = nglbtr + navgn

! count nonzeros in transformation on subdomain
         ! number of variables in glob
         nglbv = nglbn * ndofn
         ! rectangular block that gets dense by transformation
         block = sparsity1D * ndofs * nglbv

         nnz_transform_est = nnz_transform_est + block

         indinglb = indinglb + nglbn
      end do
      select case (matrixtype)
         case(0)
            koef = 3
         case(1,2)
            koef = 2
         case default
            write(*,*) 'Unknown type of matrix. Maybe MATRIXTYPE not set.'
            stop
      end select
      nnz_transform_est = koef*nnz_transform_est

! Virtual array NNGLB for problem in new variables
      lnnglbtr = nglbtr
      allocate(nnglbtr(lnnglbtr))
      nnglbtr = 1
      linglbtr = nglbtr
      allocate(inglbtr(linglbtr))
      indinglb   = 0
      indinglbtr = 0
      do iglb = 1,nglb
         nglbn = nnglb(iglb)

         ! number of averages on glob
         navgn = 1
         do iavgn = 1,navgn
            inglbtr(indinglbtr + iavgn) = inglb(indinglb + iavgn)
         end do

         indinglb   = indinglb   + nglbn
         indinglbtr = indinglbtr + navgn
      end do

! Check of array construction
      if (any(inglbtr.eq.0)) then
         write(*,*) 'bddc_T_nnz_est: Zeros in INGLBTR.'
         stop
      end if

! Call estimation of space for projection after transformation
      call bddc_P_nnz_est(matrixtype,sparsity,ndofs,ndim,nglbtr,inglbtr,linglbtr,nnglbtr,lnnglbtr,&
                          slavery,lslavery, nnz_proj_est)

! Sum the memory demands
      nnz_tr_proj_est = nnz_transform_est + nnz_proj_est

      write(*,*) 'myid =',myid,': Space estimated for transformation =',nnz_transform_est
      call flush(6)
      write(*,*) 'myid =',myid,': Space estimated for projection after transformation =',nnz_proj_est
      call flush(6)

      deallocate(nnglbtr)
      deallocate(inglbtr)

      return
end subroutine

!***************************************************************************************************************
subroutine bddc_T_init(myid,comm,ndim,nglb,inglb,linglb,nnglb,lnnglb,slavery,lslavery,nnodt,ndoft,kdoft,lkdoft,&
                       iterate_on_transformed,matrixtype,nnz,i_sparse,j_sparse,a_sparse,la, nnz_transform, &
                       solt,lsolt, nnz_proj)
!***************************************************************************************************************
! Subroutine for setting up the change of variables
!!!!! Works only for SYMMETRIC matrices!!!!!
! Storage in A after transformation:
! A | AT | TA | TAT | projected_matrix (cf bddc_P_init)|
!nnz| nnz_transform |            nnz_proj              |
! For symmetric case, block T*A  is obtained by flipping block A*T
! along diagonal and double entries on diagonal.
! Blocks of transformation are assemblaged.

! Use sparse matrices module
      use module_sm
      use module_dd
      use module_utils
       
      implicit none
      include "mpif.h"

! Process number
      integer,intent(in) :: myid

! MPI communicator
      integer, intent(in) :: comm

! Number of dimensions
      integer,intent(in) :: ndim
! Description of globs
      integer,intent(in) :: nglb
      integer,intent(in) :: linglb,        lnnglb
      integer,intent(in) ::  inglb(linglb), nnglb(lnnglb)

! Description of spcae W_tilde
      integer,intent(in) :: nnodt, ndoft
      integer,intent(in) :: lslavery,          lkdoft
      integer,intent(in) ::  slavery(lslavery), kdoft(lkdoft)

! Iterate on transformed problem?
      logical,intent(in) :: iterate_on_transformed

! Type of matrix
! 0 - unsymmetric                 
! 1 - symmetric positive definite 
! 2 - symmetric general          
      integer,intent(in) :: matrixtype

! Matrix in sparse format IJA of local triplets
      integer,intent(in)     :: nnz, la
      integer,intent(inout)  :: i_sparse(la), j_sparse(la)
      real(kr),intent(inout) :: a_sparse(la)

! Solution for transformation
      integer,intent(in)     :: lsolt
      real(kr),intent(inout) ::  solt(lsolt)

! Number of values added into matrix from change of variables
      integer,intent(out) :: nnz_transform
! Number of values added into matrix from projection after transformation
      integer,intent(out) :: nnz_proj

! Local variables
      integer  :: lavg1, lavg2
      real(kr),allocatable :: avg(:,:)
      integer  :: lt1, lt2
      real(kr),allocatable :: t(:,:)
      integer  :: ltdof1, ltdof2
      real(kr),allocatable :: tdof(:,:)
      integer  :: liglbvgvn
      integer,allocatable :: iglbvgvn(:)
      integer  ::           lnnglbtr
      integer,allocatable :: nnglbtr(:)
      integer  ::           linglbtr
      integer,allocatable :: inglbtr(:)

      integer :: idofn, iglb, indnglb, ndofn, nglbn, indinglb, &
                 iglbn, jglbn, nglbv, iv, jv, &
                 ibuf, ia, &
                 point, space_left, store_type, &
                 iglbntr, nglbntr, indinglbtr, nglbtr, navgn, indi, indin, &
                 indj, indjn, jdofn, pointi, pointj, isub
      integer :: nnz_right, nnz_left, nnz_new, nnz_add, nnz_transform_as, &
                 nnz_t_est, nnz_t_est_loc, point_t, space_t_left, nnz_t_add
      logical :: correct_sol = .false.
      real(kr):: val, valnew

      real(kr):: time

! Check the applicability of the routine - symmetric matrices
      if (matrixtype.ne.1.and.matrixtype.ne.2) then
         if (myid.eq.0) then
            write(*,*) 'bddc_P_init: This routine currently works only with symmetric matrices.'
         end if
         stop
      end if

! Application of transformation matrix from the right A*T
      call bddc_time_start(comm)
      indinglb  = 0
      nnz_new   = 0
      if (iterate_on_transformed) then
! transform initial solution if nonzero
         if (any(solt.ne.0._kr)) then
            correct_sol = .true.
         else
            correct_sol = .false.
         end if
      end if
      ! count new number of constraints
      nglbtr    = 0
      ! nonzeros in T for subdomain
      nnz_t_est = 0
      do iglb = 1,nglb
         nglbn   = nnglb(iglb)
         indnglb = inglb(indinglb + 1)
         ndofn   = ndim


         if (use_arithmetic_averages) then
            ! number of averages on glob
            navgn = 1
   
            ! prepare matrix AVG - local for glob, only averages => rectangular
            lavg1 = navgn
            lavg2 = nglbn
            allocate(avg(lavg1,lavg2))
            ! case of single aritmetic average
            avg(1,:) = 1._kr
   
            ! Number of constraints become the number of nodes in constraint
            nglbtr = nglbtr + navgn
   
   !         write(*,*) 'Matrix Tinv'
   !         do i = 1,lt1
   !            write(*,'(30f10.5)') t(i,:)
   !         end do
   !         call flush(6)
   
            ! correct solution if desired
            ! u^hat = T u
            if (iterate_on_transformed.and.correct_sol) then
               do iglbn = 1,navgn
                  indin = inglb(indinglb + iglbn)
                  pointi = kdoft(indin)
                  do idofn = 1,ndofn
                     indi = pointi + idofn
   
                     valnew = 0._kr
                     do jglbn = 1,nglbn
                        indjn = inglb(indinglb + jglbn)
                        pointj = kdoft(indjn)
   
                        val = avg(iglbn,jglbn)
                        do jdofn = 1,ndofn
                           indj = pointj + jdofn
   
                           if (solt(indj).ne.0._kr.and.val.ne.0._kr) then
                              valnew = valnew + val*solt(indj)
                           end if
                        end do
                     end do
                     solt(indi) = valnew
                  end do
               end do
            end if
   
   ! Find inverse of the matrix T
            lt1 = nglbn
            lt2 = nglbn
            allocate(t(lt1,lt2))
            call bddc_T_inverse(avg,lavg1,lavg2,t,lt1,lt2)
            deallocate(avg)
   
   !         write(*,*) 'Matrix Tinv after inversion'
   !         do i = 1,lt1
   !            write(*,'(30f10.5)') t(i,:)
   !         end do
   !         call flush(6)
   
   ! Inflate the inversion of matrix from nodes to globs
            nglbv = ndofn*nglbn
            ltdof1 = nglbv
            ltdof2 = nglbv
            allocate(tdof(ltdof1,ltdof2))
            tdof = 0._kr
            iv = 0
            do iglbn = 1,nglbn
               jv = 0
               do jglbn = 1,nglbn
                  do idofn = 1,ndofn
                     tdof(iv + idofn, jv + idofn) = t(iglbn,jglbn)
                  end do
                  jv = jv + ndofn
               end do
               iv = iv + ndofn
            end do
            deallocate(t)
         else if (use_adaptive_averages) then
            ! number of averages on glob
            ! same number of subdomains and processors is enforced
            isub = myid + 1
            call dd_get_adaptive_constraints_size(myid,isub,iglb,lavg1,lavg2)
            print *,'I am here: lavg1 =',lavg1,'lavg2 =',lavg2

            ! prepare matrix AVG - local for glob, only averages => rectangular
            allocate(avg(lavg1,lavg2))
            ! case of single aritmetic average
            avg(1,:) = 1._kr
   
            ! Number of constraints become the number of nodes in constraint
            nglbtr = nglbtr + navgn
   
   !         write(*,*) 'Matrix Tinv'
   !         do i = 1,lt1
   !            write(*,'(30f10.5)') t(i,:)
   !         end do
   !         call flush(6)
   
            ! correct solution if desired
            ! u^hat = T u
            if (iterate_on_transformed.and.correct_sol) then
               do iglbn = 1,navgn
                  indin = inglb(indinglb + iglbn)
                  pointi = kdoft(indin)
                  do idofn = 1,ndofn
                     indi = pointi + idofn
   
                     valnew = 0._kr
                     do jglbn = 1,nglbn
                        indjn = inglb(indinglb + jglbn)
                        pointj = kdoft(indjn)
   
                        val = avg(iglbn,jglbn)
                        do jdofn = 1,ndofn
                           indj = pointj + jdofn
   
                           if (solt(indj).ne.0._kr.and.val.ne.0._kr) then
                              valnew = valnew + val*solt(indj)
                           end if
                        end do
                     end do
                     solt(indi) = valnew
                  end do
               end do
            end if
   
   ! Find inverse of the matrix T
            lt1 = nglbn
            lt2 = nglbn
            allocate(t(lt1,lt2))
            call bddc_T_inverse(avg,lavg1,lavg2,t,lt1,lt2)
            deallocate(avg)
   
   !         write(*,*) 'Matrix Tinv after inversion'
   !         do i = 1,lt1
   !            write(*,'(30f10.5)') t(i,:)
   !         end do
   !         call flush(6)
   
   ! Inflate the inversion of matrix from nodes to globs
            nglbv = ndofn*nglbn
            ltdof1 = nglbv
            ltdof2 = nglbv
            allocate(tdof(ltdof1,ltdof2))
            tdof = 0._kr
            iv = 0
            do iglbn = 1,nglbn
               jv = 0
               do jglbn = 1,nglbn
                  do idofn = 1,ndofn
                     tdof(iv + idofn, jv + idofn) = t(iglbn,jglbn)
                  end do
                  jv = jv + ndofn
               end do
               iv = iv + ndofn
            end do
            deallocate(t)
         else 
            write(*,*) 'BDDC_T_INIT: No averages implied.'
            call error_exit
         end if


!         write(*,*) 'Matrix Tdof after inflating'
!         do i = 1,ltdof1
!            write(*,'(30f10.5)') tdof(i,:)
!         end do
!         call flush(6)

! Subtract unity from matrix TDOF
! This small trick allows the matrix to only generate new values without
! changing values of the matrix before transformation.
         do iv = 1,nglbv
            tdof(iv,iv) = tdof(iv,iv) - 1._kr
         end do

!         write(*,*) 'Matrix Tdof after subtracting identity'
!         do i = 1,ltdof1
!            write(*,'(30f10.5)') tdof(i,:)
!         end do
!         call flush(6)

! Find mapping of glob variables into global variables
         liglbvgvn = nglbv
         allocate(iglbvgvn(liglbvgvn))
         iv = 0
         do iglbn = 1,nglbn
            indnglb = inglb(indinglb + iglbn)
            point   = kdoft(indnglb)
            do idofn = 1,ndofn
               iglbvgvn(iv+idofn) = point+idofn
            end do
            iv = iv + ndofn
         end do

! Find new entries in the matrix from Ac*T
         point = nnz + nnz_new
         space_left = la - point 
         store_type = 0
         call sm_from_sm_mat_mult_emb(matrixtype, nnz, i_sparse, j_sparse, a_sparse, nnz, ndoft, ndoft, &
                                      tdof,ltdof1,ltdof2, iglbvgvn,liglbvgvn, iglbvgvn,liglbvgvn, &
                                      store_type,i_sparse(point+1),j_sparse(point+1),a_sparse(point+1),space_left,nnz_add)
         nnz_new = nnz_new + nnz_add

         write(idtr) ltdof1, ltdof2
         write(idtr) iglbvgvn
         write(idtr) tdof

         deallocate(tdof)
         deallocate(iglbvgvn)

         ! add nonzeros in T to counter
         nnz_t_est_loc = ndofn*navgn*nglbn + 2*(nglbv - navgn*ndofn)
         nnz_t_est     = nnz_t_est + nnz_t_est_loc

         indinglb = indinglb + nglbn
      end do
      nnz_right = nnz_new

!      write(*,*) 'Matrix without transformed entries'
!      call sm_print(6, i_sparse, j_sparse, a_sparse, la, nnz)
!      write(*,*) 'Matrix with transformed entries','nnz_right =',nnz_right
!      call sm_print(6, i_sparse, j_sparse, a_sparse, la, nnz + nnz_right)

      call bddc_time_end(comm,time)
      if (myid.eq.0.and.time_verbose.ge.2) then
         write(*,*) '======================================================='
         write(*,*) 'Time of right multiplication by transformation = ',time
         write(*,*) '======================================================='
         call flush(6)
      end if

      ! Set file with tranformation matrices to beginning
      rewind idtr
      ! Prepare sparse matrix T
      write(*,*) 'myid =',myid,': Estimated number of nonzeros in transformation T =',nnz_t_est
      call flush(6)
      lt_sparse = nnz_t_est
      allocate(i_t_sparse(lt_sparse),j_t_sparse(lt_sparse),t_sparse(lt_sparse))

! Application of transformation matrix from the left T^T*(Ac*T)
      call bddc_time_start(comm)
      indinglb  = 0
      nnz_new   = 0
      nnz_t     = 0
      do iglb = 1,nglb
         nglbn   = nnglb(iglb)
         indnglb = inglb(indinglb + 1)
         ndofn   = ndim
 
         ! Read transforming matrix and its embedding from file
         read(idtr) ltdof1, ltdof2
         liglbvgvn = ltdof1
         allocate(tdof(ltdof1,ltdof2),iglbvgvn(liglbvgvn))
         read(idtr) iglbvgvn
         read(idtr) tdof

! Find new entries in the matrix from T^T*(Ac*T)
         point = nnz + nnz_right + nnz_new
         space_left = la - point 
         store_type = 0
         ! only multiplication by (Ac*T) part of the sparse matrix
         call sm_from_sm_mat_mult_emb(store_type,nnz_right,j_sparse(nnz+1),i_sparse(nnz+1),a_sparse(nnz+1),nnz_right,ndoft,ndoft,&
                                      tdof,ltdof1,ltdof2, iglbvgvn,liglbvgvn, iglbvgvn,liglbvgvn, &
                                      matrixtype,i_sparse(point+1),j_sparse(point+1),a_sparse(point+1),space_left,nnz_add)
         nnz_new = nnz_new + nnz_add

         ! convert tdof to sparse T
         point_t      = nnz_t
         space_t_left = lt_sparse - point_t 
         t_type = 0
         call sm_from_dm_emb(t_type, tdof,ltdof1,ltdof2, iglbvgvn,liglbvgvn, iglbvgvn,liglbvgvn,&
                             i_t_sparse(point_t+1), j_t_sparse(point_t+1), t_sparse(point_t+1), space_t_left, nnz_t_add)
         nnz_t = nnz_t + nnz_t_add

         deallocate(tdof)
         deallocate(iglbvgvn)

         indinglb = indinglb + nglbn
      end do
      nnz_left = nnz_new
      
      ! Transpose the lower triangle in entries of (Ac*T) and multiply diagonal by 2
      do ia = nnz+1, nnz+nnz_right
         if (i_sparse(ia).gt.j_sparse(ia)) then
            ibuf = i_sparse(ia)
            i_sparse(ia) = j_sparse(ia)
            j_sparse(ia) = ibuf
         end if
         if (i_sparse(ia).eq.j_sparse(ia)) then
            a_sparse(ia) = 2._kr * a_sparse(ia)
         end if
      end do

      ! Add number of new entries
      write(*,*) 'right',nnz_right,'left',nnz_left
      nnz_transform = nnz_right + nnz_left

!      write(*,*) 'Matrix with transformed entries','nnz_right =',nnz_right,' nnz_left=',nnz_left
!      call sm_print(6, i_sparse, j_sparse, a_sparse, la, nnz + nnz_transform)

      write(*,*) 'myid =',myid,': Final number of nonzeros in transformation T =',nnz_t
      call flush(6)

      call bddc_time_end(comm,time)
      if (myid.eq.0.and.time_verbose.ge.2) then
         write(*,*) '======================================================='
         write(*,*) 'Time of left multiplication by transformation = ',time
         write(*,*) '======================================================='
         call flush(6)
      end if

      ! Prepare array INGLB and NNGLB for changed variables
      lnnglbtr = nglbtr
      allocate(nnglbtr(lnnglbtr))
      ! only one constraint per glob
      nnglbtr = 1

      ! Construct array INGLB - assumption of averages in FIRST nnglbtr(iglb) lines
      linglbtr = nglbtr
      allocate(inglbtr(linglbtr))
      indinglb   = 0
      indinglbtr = 0
      do iglb = 1,nglb
         nglbn   = nnglb(iglb)
         nglbntr = nnglbtr(iglb)

         ! copy first NGLBNTR entries of INGLB into INGLBTR
         do iglbntr = 1,nglbntr
            inglbtr(indinglbtr + iglbntr) = inglb(indinglb + iglbntr)
         end do

         indinglb   = indinglb   + nglbn
         indinglbtr = indinglbtr + nglbntr
      end do

      ! check of non-zeros in NNGLBTR
      if (any(nnglbtr.eq.0)) then
         write(*,*) 'bddc_T_init: Zeros in transformed number of nodes in glob.'
         stop
      end if
      ! check of non-zeros in INGLBTR
      if (any(inglbtr.eq.0)) then
         write(*,*) 'bddc_T_init: Zeros in transformed glob description.'
         stop
      end if

! Assembly new entries in matrix
      write(*,*) 'myid =',myid,': Space really needed by transformation =',nnz_transform
      call flush(6)
      call sm_assembly(i_sparse(nnz+1),j_sparse(nnz+1),a_sparse(nnz+1),nnz_transform, nnz_transform_as)
      nnz_transform = nnz_transform_as
      write(*,*) 'myid =',myid,': Number of nonzeros from transformation nnz_transform =',nnz_transform
      write(*,*) 'myid =',myid,': Local fill-in by transformation nnz_transform/nnz =',float(nnz_transform)/nnz
      call flush(6)


!      call sm_print(6,i_sparse,j_sparse,a_sparse,la,nnz+nnz_transform)
! Set-up projection after transformation
      call bddc_time_start(comm)
      call bddc_P_init(myid,comm,ndim,nglbtr,inglbtr,linglbtr,nnglbtr,lnnglbtr,&
                       slavery,lslavery,nnodt,ndoft,kdoft,lkdoft,&
                       matrixtype,nnz+nnz_transform,i_sparse,j_sparse,a_sparse,la, nnz_proj)
      call bddc_time_end(comm,time)
      if (myid.eq.0.and.time_verbose.ge.1) then
         write(*,*) '================================================================'
         write(*,*) 'Time of construction of projection after transformation = ',time
         write(*,*) '================================================================'
         call flush(6)
      end if

      deallocate(nnglbtr)
      deallocate(inglbtr)

      return
end subroutine

!************************************************
subroutine bddc_T_apply(comm,vec,lvec,transposed)
!************************************************
! Subroutine for application of transformation into new variables and its transpose.
! Module for using sparse matrices.
      use module_sm
       
      implicit none
      include "mpif.h"

! MPI communicator
      integer, intent(in) :: comm

! Vector to be transformed
      integer,intent(in)     :: lvec
      real(kr),intent(inout) ::  vec(lvec)

! Request of transposed operator
      logical, optional, intent(in) :: transposed

! Local variables
      integer  :: lbuf2
      real(kr),allocatable :: buf2(:)

      integer  :: ierr

      logical :: use_transpose

! setting-up transposition if desired
      if (present(transposed)) then
         use_transpose = transposed
      else
         use_transpose = .false.
      end if

! prepare buffer
      lbuf2 = lvec
      allocate(buf2(lbuf2))

! apply transformation
      if (use_transpose) then
         call sm_vec_mult(t_type, nnz_t, j_t_sparse, i_t_sparse, t_sparse, lt_sparse, &
                          vec,lvec, buf2,lbuf2)
      else
         call sm_vec_mult(t_type, nnz_t, i_t_sparse, j_t_sparse, t_sparse, lt_sparse, &
                          vec,lvec, buf2,lbuf2)
      end if

! communicate buffer with correction to vec
      call MPI_ALLREDUCE(buf2,buf,lbuf,MPI_DOUBLE_PRECISION,MPI_SUM,comm,ierr)

! correct the vector
      vec = vec + buf

! clear memory
      deallocate(buf2)

      return
end subroutine

!***************************************************
subroutine bddc_T_inverse(avg,lavg1,lavg2,t,lt1,lt2)
!***************************************************
! Subroutine for generation of transformation matrix T based on averages on glob
! stored in AVG.
! Performs LU factorization with column permutation 
! of a rectangular matrix AVG(LAVG1,LAVG2), LAVG1 <= LAVG2
! in first LAVG1 columns appears an invertible upper triangular block
! ang lower block of L without identity on diagonal
       
      implicit none

! rectangular matrix
      integer, intent(in)    :: lavg1, lavg2
      real(kr),intent(inout) ::  avg(lavg1,lavg2)

! permutations
      integer,intent(in)   :: lt1, lt2
      real(kr),intent(out) ::  t(lt1,lt2)

! local variables
      integer             :: irow, indrow
      integer             :: lipiv
      integer,allocatable ::  ipiv(:)
      integer              :: laux
      real(kr),allocatable ::  aux(:)
       
! Check dimensions of arrays
      if (lavg1.gt.lavg2) then
         write(*,*) 'bddc_T_inverse: Wrong shape - matrix has more rows than columns.', lavg1, lavg2
         stop
      end if
      if (lt2.ne.lavg2) then
         write(*,*) 'bddc_T_inverse: Array T has different number of columns than AVG.', lt2, lavg2
         stop
      end if

      ! Prepare array for permutations
      lipiv = lavg2
      allocate(ipiv(lipiv))

      call bddc_T_LU(avg,lavg1,lavg2,ipiv,lipiv)
      call bddc_T_blinv(avg,lavg1,lavg2)

      ! construct new array T
      t = 0._kr
      t(1:lavg1,:) = avg
      do irow = lavg1 + 1,lt1
         t(irow,irow) = 1._kr
      end do

      ! permute the matrix
      laux = lt2
      allocate(aux(laux))
      do irow = 1,lt1
         indrow = ipiv(irow)
         if (indrow.gt.irow) then
            ! exchange rows
            aux         = t(irow,:)
            t(irow,:)   = t(indrow,:)
            t(indrow,:) = aux
         end if
      end do

      deallocate(aux)
      deallocate(ipiv)

      return
end subroutine

!***********************************************
subroutine bddc_T_LU(avg,lavg1,lavg2,ipiv,lipiv)
!***********************************************
! Subroutine for LU factorization with column permutation 
! of a rectangular matrix AVG(LAVG1,LAVG2), LAVG1 <= LAVG2
! in first LAVG1 columns appears an invertible upper triangular block
! ang lower block of L without identity on diagonal
       
      implicit none

! rectangular matrix
      integer, intent(in)    :: lavg1, lavg2
      real(kr),intent(inout) ::  avg(lavg1,lavg2)

! permutations
      integer,intent(in)  :: lipiv
      integer,intent(out) ::  ipiv(lipiv)

! Local variables
      integer  :: locaux(1), pivot, ibuf, irow, jcol, jrow
      real(kr) :: coef
      
      integer              :: laux
      real(kr),allocatable ::  aux(:)

! Check dimensions of arrays
      if (lavg1.gt.lavg2) then
         write(*,*) 'bddc_T_perm: Wrong shape - matrix has more rows than columns.', lavg1, lavg2
         stop
      end if
      if (lipiv.ne.lavg2) then
         write(*,*) 'bddc_T_perm: Array IPIV has different length than number of columns.', lipiv, lavg2
         stop
      end if

! Initial array of permutations
      do jcol = 1,lavg2
         ipiv(jcol) = jcol
      end do

! buffer for column exchange
      laux = lavg1
      allocate(aux(laux))

! Gauss elimination
      ! loop over rows
      do irow = 1,lavg1
         ! find pivot on the line
         locaux = maxloc(abs(avg(irow,irow:)))
         pivot  = irow - 1 + locaux(1)
         if (abs(avg(irow,pivot)).lt.numerical_zero) then
            write(*,*) 'bddc_T_perm: Pivot less than numerical zero => dependent rows.'
            stop
         end if

         if (pivot.ne.irow) then
            ! column exchange
            aux          = avg(:,irow)
            avg(:,irow)  = avg(:,pivot)
            avg(:,pivot) = aux
            ibuf         = ipiv(irow)
            ipiv(irow)   = ipiv(pivot)
            ipiv(pivot)  = ibuf
         end if

         ! loop in the column
         do jrow = irow + 1,lavg1
            ! eliminate the line
            coef = avg(jrow,irow) / avg(irow,irow)

            ! generate line of matrix U
            avg(jrow,irow+1:) = avg(jrow,irow+1:) - coef * avg(irow,irow+1:)

            ! find the entry of L
            avg(jrow,irow) = coef
         end do
      end do

      deallocate(aux)

      return
end subroutine

!***************************************
subroutine bddc_T_blinv(avg,lavg1,lavg2)
!***************************************
! Subroutine for block inversion of the first square block
! of a rectangular matrix AVG(LAVG1,LAVG2), LAVG1 <= LAVG2
! in first LAVG1 columns assumes LU factors
! At the end, inverse of this block is put there
! initial matrix (   A      B   )
! virtually      (   0      I   )
!
! result         ( A^-1  -A^-1*B)
! virtually      (   0      I   )
       
      implicit none

! rectangular matrix
      integer, intent(in)    :: lavg1, lavg2
      real(kr),intent(inout) ::  avg(lavg1,lavg2)

! Local variables
      integer  :: lapack_info, jcol, ldim

      integer             :: lipiv
      integer,allocatable ::  ipiv(:)

      integer              :: lwork
      real(kr),allocatable ::  work(:)
      
! Check dimensions of arrays
      if (lavg1.gt.lavg2) then
         write(*,*) 'bddc_T_blinv: Wrong shape - matrix has more rows than columns.', lavg1, lavg2
         stop
      end if

! Prepare array of row permutations - trivial
      lipiv = lavg1
      allocate(ipiv(lipiv))
      do jcol = 1,lavg1
         ipiv(jcol) = jcol
      end do

! Invert the first block of AVG using LAPACK
      lwork = max(1,lavg1)
      allocate(work(lwork))
      ldim  = max(1,lavg1)
      call DGETRI(lavg1, avg, ldim, ipiv, work, lwork, lapack_info)
      if (lapack_info.ne.0) then
         write(*,*) 'Error in Lapack factorization of transformation matrix'
         stop
      end if
      deallocate(work)
      deallocate(ipiv)

! Calculate -A^-1*B
      if (lavg2.gt.lavg1) then
         avg(:,lavg1+1:) = matmul(avg(:,1:lavg1),-1._kr*avg(:,lavg1+1:))
      end if

      return
end subroutine

!************************************************************************************
subroutine bddc_RRT(comm,nnodt,nndft,lnndft,slavery,lslavery,kdoft,lkdoft,vect,lvect)
!************************************************************************************
! Subroutine for realization of operator R * R^T: W_tilde -> W_tilde on array of
! dimension of W_tilde 
! Interchange of data across interface.
      implicit none
      include "mpif.h"

! MPI communicator
      integer, intent(in) :: comm
! Description of space W_tilde
      integer, intent(in) :: nnodt, lnndft, lslavery, lkdoft, lvect
      integer, intent(in) :: nndft(lnndft), slavery(lslavery), kdoft(lkdoft)
! Vector to be communicated
      real(kr), intent(inout) :: vect(lvect)

! Local variables
      integer :: ierr
      logical :: sweep = .false.

! R^T * vect => vect
      call bddc_RT(nnodt,nndft,lnndft,slavery,lslavery,kdoft,lkdoft,vect,lvect)
! R * R^T * vect => vect
      call bddc_R(nnodt,nndft,lnndft,slavery,lslavery,kdoft,lkdoft,vect,lvect)

! Interchange data across interfaces.
!***************************************************************************MPI
      ! if buffer is not prepared, allocate it, use it and destroy it
      if (.not.allocated(buf)) then
         lbuf = lvect
         allocate(buf(lbuf))
         sweep = .true.
      end if
      if (lbuf.ne.lvect) then
         write(*,*) 'bddc_RRT: MPI Buffer size do not match vector size.'
         stop
      end if
      call MPI_ALLREDUCE(vect,buf,lvect,MPI_DOUBLE_PRECISION,MPI_SUM,comm,ierr)
      vect = buf
      ! clean array if it was not allocated before the routine
      if (sweep) then
         lbuf = 0
         deallocate(buf)
         sweep = .false.
      end if
!***************************************************************************MPI

      return
end subroutine

!******************************************************************************
subroutine bddc_RT(nnodt,nndft,lnndft,slavery,lslavery,kdoft,lkdoft,vect,lvect)
!******************************************************************************
! Subroutine for realization of operator R^T: W_tilde -> W_hat on array of
! dimension W_tilde - slave entries are summed up to master, then have no meaning
      implicit none

      integer,intent(in)     :: nnodt, lnndft, lslavery, lkdoft, lvect
      integer,intent(in)     :: nndft(lnndft), slavery(lslavery), kdoft(lkdoft)
      real(kr),intent(inout) :: vect(lvect)

! Local variables
      integer:: indvt, jndvt, inodt, ndofn

      do inodt = 1,nnodt
         if (slavery(inodt).ne.0.and.slavery(inodt).ne.inodt) then
            jndvt = kdoft(inodt)
            indvt = kdoft(slavery(inodt))
            ndofn = nndft(inodt)
            vect(indvt+1 : indvt+ndofn) = vect(indvt+1 : indvt+ndofn) + vect(jndvt+1 : jndvt+ndofn) 
         end if
      end do

      return
end subroutine

!*****************************************************************************
subroutine bddc_R(nnodt,nndft,lnndft,slavery,lslavery,kdoft,lkdoft,vect,lvect)
!*****************************************************************************
! Subroutine for realization of operator R: W_hat -> W_tilde on array of
! dimension W_tilde - slave entries have no meaning before the routine
      implicit none

      integer,intent(in) :: nnodt, lnndft, lslavery, lkdoft, lvect
      integer,intent(in) :: nndft(lnndft), slavery(lslavery), kdoft(lkdoft)
      real(kr),intent(inout) :: vect(lvect)

! Local variables
      integer:: indvt, jndvt, inodt, ndofn

      do inodt = 1,nnodt
         if (slavery(inodt).ne.0) then
            jndvt = kdoft(inodt)
            indvt = kdoft(slavery(inodt))
            ndofn = nndft(inodt)
            vect(jndvt+1 : jndvt+ndofn) = vect(indvt+1 : indvt+ndofn)
         end if
      end do

      return
end subroutine

!*********************************************************************************
subroutine bddc_R_int(nnodt,nndft,lnndft,slavery,lslavery,kdoft,lkdoft,vect,lvect)
!*********************************************************************************
! Subroutine for realization of operator R: W_hat -> W_tilde on INTEGER array of
! dimension W_tilde - slave entries have no meaning before the routine
      implicit none

      integer,intent(in)    :: nnodt, lnndft, lslavery, lkdoft, lvect
      integer,intent(in)    :: nndft(lnndft), slavery(lslavery), kdoft(lkdoft)
      integer,intent(inout) :: vect(lvect)

! Local variables
      integer:: indvt, jndvt, inodt, ndofn

      do inodt = 1,nnodt
         if (slavery(inodt).ne.0) then
            jndvt = kdoft(inodt)
            indvt = kdoft(slavery(inodt))
            ndofn = nndft(inodt)
            vect(jndvt+1 : jndvt+ndofn) = vect(indvt+1 : indvt+ndofn)
         end if
      end do

      return
end subroutine

!*************************************
function bddc_normvec(comm,vect,lvect)
!*************************************
! Subroutine for calculation of norm of vector in W_tilde with entries copied to
! more places
! ||vec|| = sqrt(vec_tilde * D_P * vec_tilde)
      implicit none
      include "mpif.h"

! MPI communicator
      integer, intent(in) :: comm

! Norm of vector weighted by D_P
      real(kr) :: bddc_normvec

! Vector for evaluation of the weighted norm
      integer,intent(in)  :: lvect
      real(kr),intent(in) :: vect(lvect)

! Local variables
      real(kr) :: normvec2
      

      normvec2 = bddc_dot_product(comm,vect,lvect,vect,lvect)
      bddc_normvec = sqrt(normvec2)

      return
end function

!********************************************************
function bddc_dot_product(comm,vect1,lvect1,vect2,lvect2)
!********************************************************
! Subroutine for calculation of dot product of two vectors in W_tilde with entries copied to
! more places
! (vec1,vec2) = vec1 * D_P * vec2
      implicit none
      include "mpif.h"

! MPI communicator
      integer, intent(in) :: comm

! Resulting scalar product weighted by D_P
      real(kr) :: bddc_dot_product

! Vector for evaluation of the weighted scalar product
      integer,intent(in)  :: lvect1,        lvect2
      real(kr),intent(in) ::  vect1(lvect1), vect2(lvect2)

! Local variables
      real(kr) :: bddc_dot_product_loc
      integer :: ierr
      

      if (allocated(dp)) then
         bddc_dot_product_loc = sum(vect1 * dp * vect2)
         call MPI_ALLREDUCE(bddc_dot_product_loc,bddc_dot_product,1,MPI_DOUBLE_PRECISION,MPI_SUM,comm,ierr)
      else
         write (*,*) 'BDDC_DOT_PRODUCT: Weigth matrix D_P not allocated. Maybe a call to BDDC_INIT missing.'
         stop
      end if

      return
end function

!**********************************************************************************************************
subroutine bddc_G_aritmetic_size(nglb,ndim,ndoft,inglb,linglb,nnglb,lnnglb,ihntn,lihntn,slavery,lslavery, &
                                 lg1,lg2)
!**********************************************************************************************************
! Subroutine for finding dimensions of matrix G with constraints based on aritmetic averages

      implicit none

! Number of globs
      integer,intent(in) :: nglb
! Number of dimensions
      integer,intent(in) :: ndim
! Description of globs
      integer,intent(in) :: linglb, lnnglb
      integer,intent(in) :: inglb(linglb), nnglb(lnnglb)

! Description of space W_tilde 
      integer,intent(in)  :: ndoft
      integer,intent(in)  :: lihntn, lslavery
      integer,intent(in)  :: ihntn(lihntn), slavery(lslavery)

! Dimensions of matrix G
      integer,intent(out)  :: lg1, lg2

! Local variables
      integer:: nglbv, nglbvn, indinglb, indnglb, indntglb, nsubn, navgn, ndofn, iglb

! find number of lines in G
      nglbv = 0
      indinglb = 0
      do iglb = 1,nglb
      ! examine the first node of glob - in hat
         indnglb  = inglb(indinglb + 1)
         indntglb = ihntn(indnglb)
         ! number of subdomains having this node
         nsubn  = count(slavery.eq.indntglb)
         ! number of constraints of this node is one less
         navgn = nsubn - 1
         ndofn = ndim
         ! number of lines that the glob adds to the G matrix
         nglbvn = navgn * ndofn
         ! add to the total number
         nglbv = nglbv + nglbvn

         indinglb = indinglb + nnglb(iglb)
      end do
! set the sizes of matrix G
      lg1 = nglbv 
      lg2 = ndoft

      return
end subroutine

!***********************************************************************************************************
subroutine bddc_G_aritmetic_build(nglb,ndim,nnodt,inglb,linglb,nnglb,lnnglb,ihntn,lihntn,slavery,lslavery, &
                                  kdoft,lkdoft, g,lg1,lg2)
!***********************************************************************************************************
! Subroutine for creation of matrix G with constraints based on aritmetic averages

      implicit none

! Number of globs
      integer,intent(in) :: nglb
! Number of dimensions
      integer,intent(in) :: ndim
! Description of globs
      integer,intent(in) :: linglb, lnnglb
      integer,intent(in) :: inglb(linglb), nnglb(lnnglb)

! Description of space W_tilde 
      integer,intent(in)  :: nnodt
      integer,intent(in)  :: lihntn, lslavery, lkdoft
      integer,intent(in)  :: ihntn(lihntn), slavery(lslavery), kdoft(lkdoft)

! matrix G with constraints
      integer,intent(in)    :: lg1, lg2
      real(kr),intent(out)  :: g(lg1,lg2)

! Local variables
      integer:: indinglb, indnglb, indntglb, indntglb_m, point_m, indntglb_s, point_s, &
                nsubn, navgn, ndofn, iglb, nng, jg, javgn, iavgn, idofn, ing, inodt

      real(kr) :: avgval

! initialize G
      g = 0._kr

! go through globs to generate matrix G
      indinglb = 0
      jg       = 0
      do iglb = 1,nglb
         ! explore the first node of glob - in hat
         indnglb  = inglb(indinglb + 1)
         indntglb = ihntn(indnglb)
         ! number of subdomains having this node
         nsubn  = count(slavery.eq.indntglb)
         ! number of constraints if this node is one less
         navgn  = nsubn - 1
         ndofn  = ndim
         nng    = nnglb(iglb)
         ! value of aritmetic mean on average
         avgval = 1._kr/nng

         ! loop over nodes of glob
         do ing = 1,nng
            ! find master
            indnglb    = inglb(indinglb + ing)
            indntglb_m = ihntn(indnglb)
            point_m    = kdoft(indntglb_m)
            ! find slaves
            indntglb_s =  0
            javgn      = 0
            do iavgn   = 1,navgn
               do inodt = indntglb_s + 1 ,nnodt
                  if (slavery(inodt).eq.indntglb_m.and.inodt.ne.indntglb_m) then
                     indntglb_s = inodt
                     exit
                  end if
               end do
               point_s = kdoft(indntglb_s)

               do idofn = 1,ndofn
                  g(jg + javgn + idofn,point_m + idofn) =  avgval
                  g(jg + javgn + idofn,point_s + idofn) = -avgval
               end do

               javgn = javgn + ndofn
            end do
         end do

         jg = jg + ndofn*navgn
         indinglb = indinglb + nng
      end do

      return
end subroutine

!*****************************************************************************************************************
subroutine bddc_convert_ht(nnod,nnodt,nndft,lnndft,ihntn,lihntn,slavery,lslavery,kdoft,lkdoft,vec,lvec,vect,lvect)
!*****************************************************************************************************************
! Subroutine for conversion of REAL array in W_hat to W_tilde of corresponding lengths
      implicit none

      integer,intent(in)  :: nnod, nnodt, lnndft, lihntn, lslavery, lkdoft, lvec, lvect
      integer,intent(in)  :: nndft(lnndft), ihntn(lihntn), slavery(lslavery), kdoft(lkdoft)
      real(kr),intent(in) :: vec(lvec)
      real(kr),intent(out):: vect(lvect)

! Local variables
      integer:: indv, indvt, inod, ndofn, indnt

      indv = 0
      do inod = 1,nnod
         indnt = ihntn(inod)
         ndofn = nndft(indnt)
         indvt = kdoft(indnt)
         vect(indvt+1 : indvt+ndofn) = vec(indv+1 : indv+ndofn)
         indv = indv + ndofn
      end do

      call bddc_R(nnodt,nndft,lnndft,slavery,lslavery,kdoft,lkdoft,vect,lvect)

      return
end subroutine

!*********************************************************************************************************************
subroutine bddc_convert_ht_int(nnod,nnodt,nndft,lnndft,ihntn,lihntn,slavery,lslavery,kdoft,lkdoft,vec,lvec,vect,lvect)
!*********************************************************************************************************************
! Subroutine for conversion of INTEGER array in W_hat to W_tilde of corresponding lengths
      implicit none

      integer,intent(in) :: nnod, nnodt, lnndft, lihntn, lslavery, lkdoft, lvec, lvect
      integer,intent(in) :: nndft(lnndft), ihntn(lihntn), slavery(lslavery), kdoft(lkdoft)
      integer,intent(in) :: vec(lvec)
      integer,intent(out):: vect(lvect)

! Local variables
      integer:: indv, indvt, inod, ndofn, indnt

      indv = 0
      do inod = 1,nnod
         indnt = ihntn(inod)
         ndofn = nndft(indnt)
         indvt = kdoft(indnt)
         vect(indvt+1 : indvt+ndofn) = vec(indv+1 : indv+ndofn)
         indv = indv + ndofn
      end do

      call bddc_R_int(nnodt,nndft,lnndft,slavery,lslavery,kdoft,lkdoft,vect,lvect)

      return
end subroutine

!******************************************************************************************
subroutine bddc_convert_th(nnod,nndft,lnndft,ihntn,lihntn,kdoft,lkdoft,vect,lvect,vec,lvec)
!******************************************************************************************
! Subroutine for conversion of REAL array in W_tilde to W_hat of corresponding lengths
      implicit none

      integer,intent(in)  :: nnod, lnndft, lihntn, lkdoft, lvect, lvec
      integer,intent(in)  :: nndft(lnndft), ihntn(lihntn), kdoft(lkdoft)
      real(kr),intent(in) :: vect(lvect)
      real(kr),intent(out):: vec(lvec) 

! Local variables
      integer:: indv, indvt, inod, ndofn, indnt

      indv = 0
      do inod = 1,nnod
         indnt = ihntn(inod)
         ndofn = nndft(indnt)
         indvt = kdoft(indnt)
         vec(indv+1 : indv+ndofn) = vect(indvt+1 : indvt+ndofn) 
         indv = indv + ndofn
      end do

      return
end subroutine

!***********************
subroutine bddc_finalize
!***********************
! Subroutine for finalization of BDDC preconditioner
! Use module for MUMPS package
      use module_mumps
       
      implicit none

! Finalize MUMPS
      call mumps_finalize(bddc_mumps)

! Deallocate weigth matrix
      deallocate(dp)

! Deallocate MPI buffer
      deallocate(buf)

! Clean memory after using dual problem
      if (use_dual_problem) then
         ! deallocate matrix G
         deallocate(g)
         ! deallocate dual matrix
         deallocate(dualm)
         ! deallocate vector of permutations
         deallocate(ipiv)
         ! deallocate dual RHS
         deallocate(dualrhs)
      end if

! Clean memory after using projection onto null G
      if (use_projection.or.use_transform) then
         ! deallocate matrix F = R^-T * G^T
         deallocate(i_f_sparse, j_f_sparse, f_sparse)
         ! close file with transformation matrices of globs
      end if

! Clean memory after using transformation
      if (use_transform) then
         ! deallocate matrix F = R^-T * G^T
         deallocate(i_t_sparse, j_t_sparse, t_sparse)
         ! close file with transformation matrices of globs
      end if

      return
end subroutine


!*******************************************************
subroutine bddc_getfname(name1,lname1,isub,suffix,fname)
!*******************************************************
!     Prepares name of file for subdomain ISUB with structure:
!     name of problem / number of subdomain + . + suffix

      implicit none
      
! subdomain number
      integer,intent(in) :: isub

! name of the problem
      integer,intent(in) ::     lname1
      character(*),intent(in) :: name1 

! suffix of the generated file
      character(*),intent(in) :: suffix

! generated name
      character(*),intent(out) :: fname

! local variables
      integer lsuffix, lfname, zerostart, zeroend, i

      fname = ' '
      lsuffix = len(suffix)
      lfname = len(fname)
      fname(1:lname1) = name1(1:lname1)
      fname(lname1+1:lname1+1) = '/'
      zerostart = lname1+2
      zeroend = lfname-lsuffix-1
      do i = zerostart,zeroend
         fname(i:i) = '0'
      end do
      fname(zeroend+1:zeroend+1) = '.'
      fname(zeroend+2:zeroend+lsuffix+1) = suffix

      if (isub.lt.10) then
         write(fname(zeroend:zeroend),'(i1)') isub
      else if (isub.lt.100) then
         write(fname(zeroend-1:zeroend),'(i2)') isub
      else if (isub.lt.1000) then
         write(fname(zeroend-2:zeroend),'(i3)') isub
      else if (isub.lt.10000) then
         write(fname(zeroend-3:zeroend),'(i4)') isub
      else
         write(*,*) 'isub = ',isub,': Out of range for file name!'
      end if

end subroutine

!*******************************
subroutine bddc_time_start(comm)
!*******************************
! Routine that starts new time measurement
! Levels of timing work like opening and closing brackets,
! measuring time elapsed between a pair of them.
! This routine is like opening a bracket, 
! see routine BDDC_TIME_END for the opposite.

      implicit none
      include "mpif.h"
      
! MPI communicator
      integer,intent(in) :: comm
      
! Local variables
      integer  :: ierr

! add new level of timing
      level_time = level_time + 1

! check if it is not too many
      if (level_time.gt.level_time_max) then
         write(*,*) 'bddc_time: Maximal number of time levels reached.'
         stop
      end if

! measure the time and add it to the times array
!***************************************************************PARALLEL
      call MPI_BARRIER(comm,ierr)
      times(level_time) = MPI_WTIME()
!***************************************************************PARALLEL


      return
end subroutine

!**********************************
subroutine bddc_time_end(comm,time)
!**********************************
! Routine that finish time measurement of a level
! Levels of timing work like opening and closing brackets,
! measuring time elapsed between a pair of them.
! This routine is like closing a bracket, 
! see routine BDDC_TIME_START for the opposite.

      implicit none
      include "mpif.h"

! MPI communicator
      integer,intent(in) :: comm
      
! Time of run
      real(kr),intent(out) :: time
      
! Local variables
      integer  :: ierr
      real(kr) :: current_time

! measure the time
!***************************************************************PARALLEL
      call MPI_BARRIER(comm,ierr)
      current_time = MPI_WTIME()
!***************************************************************PARALLEL

! check if it is not too few
      if (level_time.le.0) then
         write(*,*) 'bddc_time: All time levels already finished.'
         stop
      end if

! find the elapsed time
      time =  current_time - times(level_time)

! subtract one level of timing
      level_time = level_time - 1

      return
end subroutine

end module module_bddc

