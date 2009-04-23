subroutine lobpcg_mvecmult_f(n,x,lx,y,ly,idoper)
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
integer,intent(inout) ::  idoper 

! local
integer :: i, indi
integer :: iinstr, isub, owner, point1, point2, point, length

! LAPACK related arrays 
integer :: lapack_info 

integer :: ldxaux
real(kr) :: xaux(lx)
real(kr) :: xaux2(lx)

logical :: i_am_slave, all_slaves

comm_calls = comm_calls + 1
print *,'I am in multiply for LOBPCG!, myid = ',comm_myid,'comm_calls =',comm_calls,'idoper =', idoper

! Check if all processors simply called the fake routine - if so, finalize
if (idoper.eq.3) then
   i_am_slave = .true.
else
   i_am_slave = .false.
end if
call MPI_ALLREDUCE(i_am_slave,all_slaves,1,MPI_LOGICAL,MPI_LAND,comm_comm,ierr)
if (all_slaves) then
   idoper = -3
   return
end if

ireq = 0
! common operations - checks and projection to null(D_ij)
if (idoper.eq.1 .or. idoper.eq.2) then
   ! check the dimensions
   if (n .ne. problemsize) then
       write(*,*) 'LOBPCG_MVECMULT_F: Vector size mismatch.'
       call error_exit
   end if
   if (lx .ne. ly) then
       write(*,*) 'LOBPCG_MVECMULT_F: Data size mismatch: lx,ly',lx,ly
       call error_exit
   end if
   if (lx .ne. n) then
       write(*,*) 'LOBPCG_MVECMULT_F: Data size mismatch. lx, n', lx, n
       call error_exit
   end if

   ! make temporary copy
   do i = 1,lx
      xaux(i) = x(i)
   end do
   
   print *, 'xaux initial = ',xaux

   ! apply prepared projection using LAPACK as (I-QQ^T)xaux
   lddij  = ldij1
   ldxaux = problemsize
   call DORMQR( 'Left', 'Transpose', ldij1, 1, ldij2, dij, lddij, &
                tau, xaux, ldxaux, &
                work,lwork, lapack_info)
   call DORMQR( 'Left', 'Non-Transpose', ldij1, 1, ldij2, dij, lddij, &
                tau, xaux, ldxaux, &
                work,lwork, lapack_info)
   do i = 1,lx
      xaux(i) = x(i) - xaux(i)
   end do

   print *, 'xaux after projection = ',xaux

   if (idoper.eq.1) then
      ! multiply by matrix A - apply (I-RE) or I - R*R^T*D_P
      !  - apply weigths
      do i = 1,problemsize
         xaux2(i) = weight(i) * xaux(i) 
      end do
      !  - apply the R^T operator
      do i = 1,problemsize
         if (pairslavery(i).ne.0.and.pairslavery(i).ne.i) then
            indi = pairslavery(i)
            xaux2(indi) = xaux2(indi) + xaux2(i) 
         end if
      end do
      !  - apply the R operator
      do i = 1,problemsize
         if (pairslavery(i).ne.0.and.pairslavery(i).ne.i) then
            indi = pairslavery(i)
            xaux2(i) = xaux2(indi)
         end if
      end do
      ! Ix - REx
      do i = 1,problemsize
         xaux(i) = xaux(i) - xaux2(i)
      end do
   end if

   print *, 'xaux after I -RE = ',xaux

   ! distribute sizes of chunks of eigenvectors
   ireq = ireq + 1
   point1 = 1
   call MPI_ISEND(xaux(point1),ndofi_i,MPI_DOUBLE_PRECISION,comm_myplace1,comm_myisub,comm_comm,request(ireq),ierr)

   ireq = ireq + 1
   point2 = ndofi_i + 1
   call MPI_ISEND(xaux(point2),ndofi_j,MPI_DOUBLE_PRECISION,comm_myplace2,comm_myjsub,comm_comm,request(ireq),ierr)
end if

! What follows is performed independently on from where I was called
   ! obtain chunks of eigenvectors to multiply them
do iinstr = 1,ninstructions
   owner = instructions(iinstr,1)
   isub  = instructions(iinstr,2)

   ireq = ireq + 1
   call MPI_IRECV(bufrecv(kbufsend(iinstr)),lbufa(iinstr),MPI_DOUBLE_PRECISION,owner,isub,comm_comm,request(ireq),ierr)
end do
nreq = ireq

call MPI_WAITALL(nreq, request, statarray, ierr)
ireq = 0

! Multiply subdomain vectors by Schur complement
do iinstr = 1,ninstructions
   owner = instructions(iinstr,1)
   isub  = instructions(iinstr,2)

   point  = kbufsend(iinstr)
   length = lbufa(iinstr)
   call dd_multiply_by_schur(comm_myid,isub,bufrecv(point),length,bufsend(point),length)
end do

! distribute multiplied vectors
do iinstr = 1,ninstructions
   owner = instructions(iinstr,1)
   isub  = instructions(iinstr,2)

   ireq = ireq + 1
   call MPI_ISEND(bufsend(kbufsend(iinstr)),lbufa(iinstr),MPI_DOUBLE_PRECISION,owner,isub,comm_comm,request(ireq),ierr)
end do

! Continue only of I was called from LOBPCG
if (idoper.eq.1 .or. idoper.eq.2) then
   ireq = ireq + 1
   point1 = 1
   call MPI_IRECV(xaux(point1),ndofi_i,MPI_DOUBLE_PRECISION,comm_myplace1,comm_myisub,comm_comm,request(ireq),ierr)

   ireq = ireq + 1
   point2 = ndofi_i + 1
   call MPI_IRECV(xaux(point2),ndofi_j,MPI_DOUBLE_PRECISION,comm_myplace2,comm_myjsub,comm_comm,request(ireq),ierr)
end if
! Wait for all vectors reach their place
call MPI_WAITALL(nreq, request, statarray, ierr)
! Continue only of I was called from LOBPCG
if (idoper.eq.1 .or. idoper.eq.2) then

   ! reverse sign of vector
   do i = 1,problemsize
      xaux(i) = -xaux(i)
   end do

   print *, 'xaux after Schur = ',xaux

   if (idoper.eq.1) then
      ! apply (I-RE)^T = (I - E^T R^T) = (I - D_P^T * R * R^T)
      ! copy the array
      do i = 1,problemsize
         xaux2(i) = xaux(i)
      end do
      !  - apply the R^T operator
      do i = 1,problemsize
         if (pairslavery(i).ne.0.and.pairslavery(i).ne.i) then
            indi = pairslavery(i)
            xaux2(indi) = xaux2(indi) + xaux2(i) 
         end if
      end do
      !  - apply the R operator
      do i = 1,problemsize
         if (pairslavery(i).ne.0.and.pairslavery(i).ne.i) then
            indi = pairslavery(i)
            xaux2(i) = xaux2(indi)
         end if
      end do
      !  - apply weigths D_P
      do i = 1,problemsize
         xaux2(i) = weight(i) * xaux2(i) 
      end do
      ! Ix - REx
      do i = 1,problemsize
         xaux(i) = xaux(i) - xaux2(i)
      end do
   end if

   print *, 'xaux after (I-RE)^T = ',xaux

   ! apply prepared projection using LAPACK as Pxaux = (I-QQ^T)xaux
   lddij  = ldij1
   ldxaux = problemsize
   call DORMQR( 'Left', 'Transpose', ldij1, 1, ldij2, dij, lddij, &
                tau, xaux, ldxaux, &
                work,lwork, lapack_info)
   call DORMQR( 'Left', 'Non-Transpose', ldij1, 1, ldij2, dij, lddij, &
                tau, xaux, ldxaux, &
                work,lwork, lapack_info)
   do i = 1,lx
      xaux(i) = x(i) - xaux(i)
   end do

   print *, 'xaux after I-QQ^T = ',xaux

   ! copy result to y
   do i = 1,lx
      y(i) = xaux(i)
   end do

end if

end subroutine
