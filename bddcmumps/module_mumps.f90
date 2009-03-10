module module_mumps
!******************
! Module for more convenient using of MUMPS
! Jakub Sistek, Denver, 12/2007

integer,parameter,private:: kr = kind(1.D0)

contains

!*************************************************
      subroutine mumps_init(mumps,comm,matrixtype)
!*************************************************
! Initializes run of MUMPS
      implicit none
      include "dmumps_struc.h"
      type(DMUMPS_STRUC),intent(inout) :: mumps
      integer,intent(in):: comm, matrixtype

! Define a communicator for this instance of MUMPS
      mumps%COMM = comm

! Set matrix type
! 0 - unsymmetric
! 1 - symmetric positive definite
! 2 - symmetric general
      mumps%SYM = matrixtype

! Set type of parallelism
! 0 - master does not work on factorization
! 1 - master works on factorization
      mumps%PAR = 1

! Initialize an instance of the package
      mumps%JOB = -1

      CALL DMUMPS(mumps)

      return
      end subroutine

!***********************************************
      subroutine mumps_set_info(mumps,mumpsinfo)
!***********************************************
! Defines level of information from MUMPS printed to screen (unit 6)
      implicit none
      include "dmumps_struc.h"
      type(DMUMPS_STRUC),intent(inout) :: mumps
      integer,intent(in):: mumpsinfo

      select case (mumpsinfo)
      case (0)
         ! errors suppressed
         mumps%ICNTL(1) = 0
         ! diagnostics suppressed
         mumps%ICNTL(2) = 0
         ! global information from host suppressed
         mumps%ICNTL(3) = 0
      case (1)
         ! errors printed on unit 6
         mumps%ICNTL(1) = 6
         ! diagnostics suppressed
         mumps%ICNTL(2) = 0
         ! global information from host suppressed
         mumps%ICNTL(3) = 0
      case (2)
         ! errors printed on unit 6
         mumps%ICNTL(1) = 6
         ! diagnostics suppressed
         mumps%ICNTL(2) = 0
         ! global information from host printed on unit 6
         mumps%ICNTL(3) = 6
      case (3)
         ! errors printed on unit 6
         mumps%ICNTL(1) = 6
         ! diagnostics printed on unit 6
         mumps%ICNTL(2) = 6
         ! global information from host printed on unit 6
         mumps%ICNTL(3) = 6
      case default
         write(*,*) 'mumps_set_info: Illegal value of MUMPSINFO (0-3):',mumpsinfo
         stop
      end select

      return
      end subroutine

!*****************************************************************************
      subroutine mumps_load_triplet(mumps,n,nnz,i_sparse,j_sparse,a_sparse,la)
!*****************************************************************************
! Associates parts of MUMPS structure with program data.
! Subroutine for loading local triplets of distributed sparse matrix in IJA format
! i_sparse, j_sparse, a_sparse are non-zero entries in IJA
! nnz =< la

      implicit none
      include "dmumps_struc.h"
      type(DMUMPS_STRUC),intent(inout) :: mumps
      integer,intent(in):: n, nnz, la
      integer,intent(in),target:: i_sparse(la), j_sparse(la)
      real(kr),intent(in),target:: a_sparse(la)

! Specify distributed matrix
      mumps%ICNTL(5) = 0
! Specify local triplet entry
      mumps%ICNTL(18) = 3

! Fill in the MUMPS structure
      ! problem dimension
      mumps%N       = n
      ! number of nonzeros on processor
      mumps%NZ_loc  = nnz
      ! row indices of matrix entries
      mumps%IRN_loc => i_sparse
      ! column indices of matrix entries
      mumps%JCN_loc => j_sparse
      ! values of matrix entries
      mumps%A_loc   => a_sparse

      return
      end subroutine

!*******************************************************************
      subroutine mumps_set_schur(mumps,listvar_schur,llistvar_schur)
!*******************************************************************
! Associates parts of MUMPS structure with program data.
! Subroutine for setting Schur complement dimension and variables

      implicit none
      include "dmumps_struc.h"
      type(DMUMPS_STRUC),intent(inout) :: mumps
      integer,intent(in)::       llistvar_schur
      integer,intent(in),target:: listvar_schur(llistvar_schur)

! Specify using Schur 
!  1 - complement centralized at host
!  2 - complement distributed - only lower triangle for symmetric case, whole matrix otherwise
!  3 - complement distributed - whole matrix for both sym and unsym case
      mumps%ICNTL(19) = 2

! Fill in the MUMPS structure
      ! Schur dimension
      mumps%SIZE_SCHUR = llistvar_schur 
      ! indices of Schur elements
      mumps%LISTVAR_SCHUR => listvar_schur

      return
      end subroutine

!************************************
      subroutine mumps_analyze(mumps)
!************************************
! Performs the analysis of matrix by MUMPS.
      implicit none
      include "dmumps_struc.h"
      type(DMUMPS_STRUC),intent(inout) :: mumps

! Job type = 1 for matrix analysis
      mumps%JOB = 1
      CALL DMUMPS(mumps)

      return
      end subroutine

!*****************************************************
      subroutine mumps_get_schur_size(mumps,mloc,nloc)
!*****************************************************
! Return size of local block of Schur complement

      implicit none
      include "dmumps_struc.h"
      type(DMUMPS_STRUC),intent(inout) :: mumps
      integer, intent(out) :: mloc, nloc

! Set the dimensions
      mloc = mumps%SCHUR_MLOC
      nloc = mumps%SCHUR_NLOC

      return
      end subroutine

!**********************************************************
      subroutine mumps_assoc_schur(mumps,schur,lschur,mloc)
!**********************************************************
! Associates parts of MUMPS structure with program data.
! Subroutine for setting Schur complement dimension and variables

      implicit none
      include "dmumps_struc.h"
      type(DMUMPS_STRUC),intent(inout) :: mumps
! Local block of Schur complement
      integer,intent(in)::        lschur
      real(kr),intent(in),target:: schur(lschur)
! Local number of rows in Schur complement
      integer,intent(in):: mloc

! Set the leading dimension for Schur complement
      mumps%SCHUR_LLD = mloc
! Associate memory for Schur complement
      mumps%SCHUR => schur

      return
      end subroutine

!**************************************
      subroutine mumps_factorize(mumps)
!**************************************
! Performs factorization of matrix by MUMPS
      implicit none
      include "dmumps_struc.h"
      type(DMUMPS_STRUC),intent(inout) :: mumps

! Job type = 2 for factorization
      mumps%JOB = 2
      CALL DMUMPS(mumps)

      return
      end subroutine

!********************************************
subroutine mumps_resolve(mumps,rhs,lrhs,nrhs)
!********************************************
! Performs backward step of multifrontal algorithm by MUMPS
! in the beginning, RHS contains RHS
! in the end, RHS contains solution

      implicit none
      include "dmumps_struc.h"
      type(DMUMPS_STRUC),intent(inout) :: mumps
      integer,intent(in):: lrhs
      real(kr),intent(inout),target:: rhs(lrhs)
      integer,intent(in),optional:: nrhs

! Associate RHS to the memory location of RHS
      mumps%RHS => rhs
      ! if number of right hand sides is specified, use it
      if (present(nrhs)) then
         mumps%NRHS = nrhs
         ! leading dimension
         mumps%LRHS = lrhs/nrhs
      end if

! Job type = 3 for backsubstitution
      mumps%JOB = 3
      CALL DMUMPS(mumps)

      ! get back to defaults
      if (present(nrhs)) then
         mumps%NRHS = 1
      end if

      return
      end subroutine
         
!*************************************
      subroutine mumps_finalize(mumps)
!*************************************
! Destroy the instance (deallocate internal data structures)
      implicit none
      include "dmumps_struc.h"
      type(DMUMPS_STRUC),intent(inout) :: mumps

      mumps%JOB = -2
      CALL DMUMPS(mumps)

      return
      end subroutine

end module module_mumps
