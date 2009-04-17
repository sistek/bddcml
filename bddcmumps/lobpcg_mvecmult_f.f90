subroutine lobpcg_mvecmult_f(n,x,lx,y,ly,idmatrix)
! realizes the multiplications with the strange matrices in the local eigenvalue problems
! for adaptive BDDC

use module_dd
use module_adaptivity_comm
use module_utils
implicit none
include "mpif.h"
integer,parameter :: kr = kind(1.D0)

! length of vector
integer,intent(in) ::   n 
integer,intent(in) ::  lx 
real(kr),intent(in) ::  x(lx)
integer,intent(in) ::  ly 
real(kr),intent(out) :: y(ly)
! determine which matrix should be multiplied
!  1 - A = P(I-RE)'S(I_RE)P
!  2 - B = PSP
!  3 - not called from LOBPCG by from fake looper - just perform demanded multiplications by S
!  -3 - set on exit if iterational process should be stopped now
integer,intent(inout) ::  idmatrix 

! local
integer :: ivec

logical :: i_am_slave = .false., all_slaves

! Check if all processors simply called the fake routine - if so, finalize
if (idmatrix.eq.3) then
   i_am_slave = .true.
end if
call MPI_ALLREDUCE(i_am_slave,all_slaves,1,MPI_LOGICAL,MPI_LAND,comm_comm,ierr)
if (all_slaves) then
   idmatrix = -3
   return
end if

if (idmatrix.eq.1) then
! multiply by matrix A

   call zero(y,ly)
   ! check the dimensions
   if (n .ne. problemsize) then
       write(*,*) 'LOBPCG_MVECMULT_F: Vector size mismatch.'
       call error_exit
   end if
   if (lx .ne. ly) then
       write(*,*) 'LOBPCG_MVECMULT_F: Data size mismatch.'
       call error_exit
   end if
   if (lx .ne. n*neigvec) then
       write(*,*) 'LOBPCG_MVECMULT_F: Data size mismatch.'
       call error_exit
   end if


else if (idmatrix.eq.2) then
! multiply by matrix B

else if (idmatrix.eq.3) then
! just work as a slave - multiply by Schur complements

end if

end subroutine
