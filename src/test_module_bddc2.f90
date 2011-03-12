! BDDCML - Multilevel BDDC
! Copyright (C) The BDDCML Team
! 
! This program is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
! 
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
! 
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <http://www.gnu.org/licenses/>.
! ____________________________________________________________________

program test_module_bddc2
! Tester of module_mumps
      use module_bddc

      implicit none
      include "mpif.h"
      integer,parameter :: kr = kind(1.D0)

! Problem dimension
      integer,parameter :: n = 6
! Length of sparse matrix
      integer,parameter :: la = 13

! Description of mesh
      integer :: lsolt = n, nnodt = n
      integer :: lnndft = n
      integer :: nndft(n) = (/ 1, 1, 1, 1, 1, 1 /)
      integer :: lihntn = n
      integer :: ihntn(n) = (/ 1, 2, 3, 4, 5, 6 /)
      integer :: lslavery = n
      integer :: slavery(n) = (/ 0, 0, 0, 0, 0, 0 /)
      integer :: lkdoft = n
      integer :: kdoft(n) = (/ 0, 1, 2, 3, 4, 5 /)

! Matrix
      real(kr) :: a_sparse(la) = (/ 4._kr, 1._kr, 1._kr, 8._kr, 1._kr, 2._kr, &
                                    4._kr, 1._kr, 4._kr, 1._kr, 8._kr, 1._kr, 4._kr /)
      integer ::  i_sparse(la) = (/ 1, 1, 1, 2, 2, 2, 3, 3, 4, 4, 5, 5, 6 /)
      integer ::  j_sparse(la) = (/ 1, 2, 4, 2, 3, 5, 3, 6, 4, 5, 5, 6, 6 /)
! Number of non-zeros
      integer :: nnz = la

! Weight approach
      integer :: weight_approach = 1

! Matrix G
      integer :: averages_approach = 0
      integer :: ndim = 2
      integer :: nglb = 0
      integer :: linglb = 0, lnnglb = 0
      integer,allocatable :: inglb(:), nnglb(:)

! Residual
      integer  :: lrest = n
      real(kr) :: rest(n) = (/ 1._kr, 2._kr, 1._kr, 1._kr, 2._kr, 1._kr /)

! Preconditioned residual
      integer  :: lht = n
      real(kr) :: ht(n) 

! Reference solution
      integer:: i
      real(kr) :: solution_ref(n)  = (/(1._kr/12._kr ,i = 1,n)/)
      real(kr) :: solution_ref2(n) = (/(1._kr/6._kr  ,i = 1,n)/)

!  parallel variables
      integer :: myid, comm, ierr

! Iterate on transformed problem?
      logical :: iterate_on_transformed = .false.

!  local variables
      integer :: matrixtype, mumpsinfo, nnz_transform = 0, timeinfo = 1
      real(kr) :: rmr
   

! MPI initialization
      call MPI_INIT(ierr)

! Communicator
      comm = MPI_COMM_WORLD
      call MPI_COMM_RANK(comm,myid,ierr)
! SPD matrix
      matrixtype = 1
! Level of information from MUMPS
      mumpsinfo = 0
! Initialize BDDC
      call bddc_init(myid,comm,matrixtype,mumpsinfo,timeinfo,iterate_on_transformed,lsolt,nnz,i_sparse,j_sparse,a_sparse,la, &
                     weight_approach,averages_approach,ndim,nglb,inglb,linglb,nnglb,lnnglb,&
                     nnodt,nndft,lnndft,ihntn,lihntn,slavery,lslavery,kdoft,lkdoft,&
                     rest,lrest, rest,lrest, nnz_transform)

! Call BDDC preconditioner
      call bddc_M(myid,comm,iterate_on_transformed,nnodt,nndft,lnndft,slavery,lslavery,kdoft,lkdoft,rest,lrest,ht,lht,rmr)

! Visual check the solution
      if (myid.eq.0) then
         write(*,*) 'Position | Solution by BDDC | Reference solution '
         write(*,'(i6,8x, f12.7,8x, f12.7)') ( i, ht(i), solution_ref(i), i = 1,lsolt )
      end if

! Interchange variables on interface
      call bddc_RRT(comm,nnodt,nndft,lnndft,slavery,lslavery,kdoft,lkdoft,ht,lht)

! Visual check the solution
      if (myid.eq.0) then
         write(*,*) 'Position | Solution by BDDC | Reference solution '
         write(*,'(i6,8x, f12.7,8x, f12.7)') ( i, ht(i), solution_ref2(i), i = 1,lsolt )
      end if

! Finalize BDDC
      call bddc_finalize
   
! Finalize MPI
      call MPI_FINALIZE(ierr)

end program
