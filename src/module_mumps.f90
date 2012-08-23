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

module module_mumps
!******************
! Module for more convenient using of MUMPS
! Jakub Sistek, Denver, 12/2007

integer,parameter,private:: kr = kind(1.D0)

! may relaxing memory be tried in factorization?
logical,parameter,private:: allow_memory_relaxation = .true.
logical,parameter,private:: debug = .true.

contains

!*************************************************
      subroutine mumps_init(mumps,comm,matrixtype)
!*************************************************
! Initializes run of MUMPS
      use module_utils
      implicit none
      include "dmumps_struc.h"
      type(DMUMPS_STRUC),intent(inout) :: mumps
      integer,intent(in):: comm, matrixtype
      ! local vars
      character(*),parameter:: routine_name = 'MUMPS_INIT'

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

! Check output
      if (mumps%INFOG(1) .ne. 0) then
         call error(routine_name,'MUMPS error',mumps%INFOG(1))
      end if

      end subroutine

!***********************************************
      subroutine mumps_set_info(mumps,mumpsinfo)
!***********************************************
! Defines level of information from MUMPS printed to screen (unit 6)
      use module_utils
      implicit none
      include "dmumps_struc.h"
      type(DMUMPS_STRUC),intent(inout) :: mumps
      integer,intent(in):: mumpsinfo
      ! local vars
      character(*),parameter:: routine_name = 'MUMPS_SET_INFO'

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
         call error(routine_name,'Illegal value of MUMPSINFO (0-3):',mumpsinfo)
      end select

      end subroutine

!*****************************************************************************************
      subroutine mumps_load_triplet_centralized(mumps,n,nnz,i_sparse,j_sparse,a_sparse,la)
!*****************************************************************************************
! Associates parts of MUMPS structure with program data.
! Subroutine for loading local triplets of distributed sparse matrix in IJA format
! i_sparse, j_sparse, a_sparse are non-zero entries in IJA
! nnz =< la

      implicit none
      include "dmumps_struc.h"
      include "mpif.h"
      type(DMUMPS_STRUC),intent(inout) :: mumps
      integer,intent(in):: n, nnz, la
      integer,intent(in),target:: i_sparse(la), j_sparse(la)
      real(kr),intent(in),target:: a_sparse(la)

! Specify assembled matrix
      mumps%ICNTL(5) = 0
! Specify centralized triplet entry
      mumps%ICNTL(18) = 0

! Fill in the MUMPS structure
      ! problem dimension
      mumps%N       = n
      ! number of nonzeros on processor
      mumps%NZ      = nnz
      ! row indices of matrix entries
      mumps%IRN     => i_sparse
      ! column indices of matrix entries
      mumps%JCN     => j_sparse
      ! values of matrix entries
      mumps%A       => a_sparse

      end subroutine

!*****************************************************************************************
      subroutine mumps_load_triplet_distributed(mumps,n,nnz,i_sparse,j_sparse,a_sparse,la)
!*****************************************************************************************
! Associates parts of MUMPS structure with program data.
! Subroutine for loading local triplets of distributed sparse matrix in IJA format
! i_sparse, j_sparse, a_sparse are non-zero entries in IJA
! nnz =< la

      implicit none
      include "dmumps_struc.h"
      include "mpif.h"
      type(DMUMPS_STRUC),intent(inout) :: mumps
      integer,intent(in):: n, nnz, la
      integer,intent(in),target:: i_sparse(la), j_sparse(la)
      real(kr),intent(in),target:: a_sparse(la)

! Specify assembled matrix
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
      mumps%ICNTL(19) = 1

! Fill in the MUMPS structure
      ! Schur dimension
      mumps%SIZE_SCHUR = llistvar_schur 
      ! indices of Schur elements
      mumps%LISTVAR_SCHUR => listvar_schur

      end subroutine

!**********************************************
      subroutine mumps_analyze(mumps,iparallel)
!**********************************************
! Performs the analysis of matrix by MUMPS.
      use module_utils
      implicit none
      include "dmumps_struc.h"
      type(DMUMPS_STRUC),intent(inout) :: mumps
      ! should analysis be parallel?
      ! 0 - let MUMPS decide
      ! 1 - force serial analysis
      ! 2 - force parallel analysis
      integer, intent(in) :: iparallel 
      ! local vars
      character(*),parameter:: routine_name = 'MUMPS_ANALYZE'

      ! analyze parameters
      if (iparallel .eq. 0) then
         ! let MUMPS decide
         mumps%ICNTL(28) = 0
      else if (iparallel .eq. 1) then
         ! serial analysis
         mumps%ICNTL(28) = 1
      else if (iparallel .eq. 2) then
         mumps%ICNTL(28) = 2
         ! Parmetis or PT-SCOTCH?
         ! 0 - automatic choice
         ! 1 - PT-SCOTCH
         ! 2 - Parmetis
         mumps%ICNTL(29) = 0
      else
         call error(routine_name,'Illegal value of iparallel:',iparallel)
      end if

! Job type = 1 for matrix analysis
      mumps%JOB = 1
      CALL DMUMPS(mumps)

! Check output
      if (mumps%INFOG(1) .ne. 0) then
         call error(routine_name,'MUMPS error',mumps%INFOG(1))
      end if

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

      end subroutine

!**********************************************************
      subroutine mumps_assoc_schur(mumps,schur,lschur,mloc)
!**********************************************************
! Associates parts of MUMPS structure with program data.
! Subroutine for setting Schur complement dimension and variables

      use module_utils
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

      end subroutine

!**************************************
      subroutine mumps_factorize(mumps)
!**************************************
! Performs factorization of matrix by MUMPS
      use module_utils
      implicit none
      include "dmumps_struc.h"
      type(DMUMPS_STRUC),intent(inout) :: mumps
      include "mpif.h"

      ! local vars
      character(*),parameter:: routine_name = 'MUMPS_FACTORIZE'
      integer :: myid, errproc, ierr, comm, errcode
      integer,parameter :: mem_relax_lim = 1000 ! limit of allowed memory relaxation (in percent)


! Job type = 2 for factorization
 13   mumps%JOB = 2
      CALL DMUMPS(mumps)

! Check output
      if (mumps%INFOG(1) .ne. 0) then
         if (allow_memory_relaxation) then
            ! find communicator
            comm = mumps%COMM
            ! orient in the communicator
            call MPI_COMM_RANK(comm,myid,ierr)
            ! if problem can be solved by relaxing memory, do it, otherwise exit on error
            ! set error process on root
            if (myid.eq.0) then
               if (mumps%INFO(1) .eq. -1) then
                  errproc = mumps%INFO(2)
               else
                  errproc = myid
               end if
            end if
            ! broadcast the ID of errorneous process
            call MPI_BCAST(errproc, 1, MPI_INTEGER, 0, comm, ierr)
            if (myid.eq.errproc) then
               errcode = mumps%INFO(1)
            end if
            ! broadcast the error code from errorneous process
            call MPI_BCAST(errcode, 1, MPI_INTEGER, errproc, comm, ierr)
            ! if the error is with local memory, try to rerun the factorization
            if (errcode.eq.-9) then
               if (myid.eq.errproc) then 
                  if (mumps%ICNTL(14) .lt. mem_relax_lim) then
                     ! add 20% to relaxation parameter up to limit
                     mumps%ICNTL(14) = mumps%ICNTL(14) + 50
                     mumps%ICNTL(14) = min(mumps%ICNTL(14),mem_relax_lim)
                     if (debug) then
                        call warning(routine_name,'Reruning factorization with memory relaxation % ',mumps%ICNTL(14))
                     end if
                  else
                     call error(routine_name,'unable to increase memory relaxation parameter for MUMPS factorization from ',&
                                mumps%ICNTL(14))
                  end if
               end if
               ! rerun the factorization
               goto 13
            end if
         end if
         call error(routine_name,'MUMPS error',mumps%INFOG(1))
      end if

      end subroutine

!*************************************************************
subroutine mumps_resolve(mumps,rhs,lrhs,nrhs,solve_transposed)
!*************************************************************
! Performs backward step of multifrontal algorithm by MUMPS
! in the beginning, RHS contains RHS
! in the end, RHS contains solution

      use module_utils
      implicit none
      include "dmumps_struc.h"
      type(DMUMPS_STRUC),intent(inout) :: mumps
      integer,intent(in):: lrhs
      real(kr),intent(inout),target:: rhs(lrhs)
      integer,intent(in),optional:: nrhs
      logical,intent(in),optional:: solve_transposed

      ! local vars
      character(*),parameter:: routine_name = 'MUMPS_RESOLVE'

! Associate RHS to the memory location of RHS
      mumps%RHS => rhs
      ! if number of right hand sides is specified, use it
      if (present(nrhs)) then
         mumps%NRHS = nrhs
         ! leading dimension
         mumps%LRHS = lrhs/nrhs
      end if

      ! default is 1 which means A x = b
      mumps%ICNTL(9) = 1
      if (present(solve_transposed)) then
         if (solve_transposed ) then
            ! should a transposed system A^T x = b be solved?
            mumps%ICNTL(9) = 0
         end if
      end if

! Job type = 3 for backsubstitution
      mumps%JOB = 3
      CALL DMUMPS(mumps)

! Check output
      if (mumps%INFOG(1) .ne. 0) then
         call error(routine_name,'MUMPS error',mumps%INFOG(1))
      end if

      ! get back to defaults
      if (present(nrhs)) then
         mumps%NRHS = 1
      end if

      end subroutine
         
!*************************************
      subroutine mumps_finalize(mumps)
!*************************************
! Destroy the instance (deallocate internal data structures)
      use module_utils
      implicit none
      include "dmumps_struc.h"
      type(DMUMPS_STRUC),intent(inout) :: mumps

      ! local vars
      character(*),parameter:: routine_name = 'MUMPS_FINALIZE'

      mumps%JOB = -2
      CALL DMUMPS(mumps)

! Check output
      if (mumps%INFOG(1) .ne. 0) then
         call error(routine_name,'MUMPS error',mumps%INFOG(1))
      end if

      end subroutine

end module module_mumps
