program test_module_mumps
! Tester of module_mumps
      use dmumps_struc_def
      use module_mumps
      implicit none
      include "mpif.h"
! Use this structure of MUMPS for routines from mumps
      type(DMUMPS_STRUC) :: test_mumps

      integer,parameter :: kr = kind(1.D0)

! Problem dimension
      integer,parameter :: n = 6
! Length of sparse matrix
      integer,parameter :: la = 13

! Matrix
      real(kr) :: a_sparse(la) = (/ 4._kr, 1._kr, 1._kr, 8._kr, 1._kr, 2._kr, &
                                    4._kr, 1._kr, 4._kr, 1._kr, 8._kr, 1._kr, 4._kr /)
      integer ::  i_sparse(la) = (/ 1, 1, 1, 2, 2, 2, 3, 3, 4, 4, 5, 5, 6 /)
      integer ::  j_sparse(la) = (/ 1, 2, 4, 2, 3, 5, 3, 6, 4, 5, 5, 6, 6 /)
! Number of non-zeros
      integer :: nnz = la

! Right hand side
      real(kr) :: rhs(n) = (/ 1._kr, 2._kr, 1._kr, 1._kr, 2._kr, 1._kr /)
! Reference solution
      integer:: i
      real(kr) :: solution_ref(n) = (/(1._kr/6._kr ,i = 1,n)/)
! Vector of solution
      integer:: lsolution = n
      real(kr) :: solution(n) 

!  parallel variables
      integer :: comm, ierr

!  local variables
      integer :: matrixtype, mumpsinfo
   

! MPI initialization
      call MPI_INIT(ierr)

! Communicator
      comm = MPI_COMM_WORLD
! SPD matrix
      matrixtype = 1
! MUMPS initialization
      call mumps_init(test_mumps,comm,matrixtype)

! Level of information from MUMPS
      mumpsinfo = 0
      call mumps_set_info(test_mumps,mumpsinfo)

! Load matrix to MUMPS
      call mumps_load_triplet(test_mumps,n,nnz,i_sparse,j_sparse,a_sparse,la)

! Analyze matrix
      call mumps_analyze(test_mumps)

! Analyze matrix
      call mumps_factorize(test_mumps)

! Solve the problem
      solution = rhs
      call mumps_resolve(test_mumps,solution,lsolution)

! Visual check the solution
      write(*,*) 'Position | Solution by MUMPS | Reference solution '
      write(*,'(i6,8x, f12.7,8x, f12.7)') ( i, solution(i), solution_ref(i), i = 1,n )

! Repeat the problem backward step
      solution = rhs
      call mumps_resolve(test_mumps,solution,lsolution)

! Visual check the solution
      write(*,*) 'Another backward step...'
      write(*,*) 'Position | Solution by MUMPS | Reference solution '
      write(*,'(i6,8x, f12.7,8x, f12.7)') ( i, solution(i), solution_ref(i), i = 1,n )

! Finalize MUMPS
      call mumps_finalize(test_mumps)
   
! Finalize MPI
      call MPI_FINALIZE(ierr)

end program
