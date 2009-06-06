module module_utils
! Module for various auxiliary utilities

implicit none
integer,parameter,private :: kr = kind(1.D0)

real(kr),parameter :: numerical_zero = 1.e-14

interface zero
   module procedure zero_int_1
   module procedure zero_int_2
   module procedure zero_real_1
   module procedure zero_real_2
end interface zero

contains

subroutine zero_int_1(ia,lia)
! zero integer array

integer,intent(in)  :: lia
integer,intent(out) ::  ia(lia)
! local
integer :: i

do i = 1,lia
   ia(i) = 0
end do
end subroutine zero_int_1

subroutine zero_int_2(ia,lia1,lia2)
! zero integer array

integer,intent(in)  :: lia1, lia2
integer,intent(out) ::  ia(lia1,lia2)
! local
integer :: i,j

do j = 1,lia2
   do i = 1,lia1
      ia(i,j) = 0
   end do
end do
end subroutine zero_int_2

subroutine zero_real_1(x,lx)
! zero real array

integer,intent(in)   :: lx
real(kr),intent(out) ::  x(lx)
! local
integer :: i

do i = 1,lx
   x(i) = 0.0_kr
end do
end subroutine zero_real_1

subroutine zero_real_2(x,lx1,lx2)
! zero real array

integer,intent(in)   :: lx1, lx2
real(kr),intent(out) ::  x(lx1,lx2)
! local
integer :: i,j

do j = 1,lx2
   do i = 1,lx1
      x(i,j) = 0.0_kr
   end do
end do
end subroutine zero_real_2

!*************************************************
subroutine getfname(problemname,isub,suffix,fname)
!*************************************************
!     Prepares name of file for subdomain ISUB with structure:
!     name of problem / number of subdomain + . + suffix

      implicit none
      
! subdomain number
      integer,intent(in) :: isub

! name of the problem
      character(*),intent(in) :: problemname 

! suffix of the generated file
      character(*),intent(in) :: suffix

! generated name
      character(*),intent(out) :: fname

! local variables
      character(4) :: numberstr

! zero filename
      fname = ' '

      ! generate name for 4 digits
      numberstr = '0000'
      if (isub.lt.10) then
         write(numberstr(4:4),'(i1)') isub
      else if (isub.lt.100) then
         write(numberstr(3:4),'(i2)') isub
      else if (isub.lt.1000) then
         write(numberstr(2:4),'(i3)') isub
      else if (isub.lt.10000) then
         write(numberstr(1:4),'(i4)') isub
      else
         write(*,*) 'isub = ',isub,': Out of range for file name!'
         fname = ' '
         return
      end if

      fname = trim(problemname)//'/'//trim(problemname)//'_'//trim(numberstr)//'.'//trim(suffix)

end subroutine

! find an empty IO unit (originally written by Jaroslav Hajek)
subroutine allocate_unit(iunit)
integer,parameter:: iunit_start = 7
! 
! purpose:      find an unused i/o unit
! arguments:
! iunit         the unit allocate
integer,intent(out):: iunit
logical:: op
iunit = iunit_start
do
  inquire(unit=iunit,opened=op)
  if (.not.op) exit
  iunit = iunit + 1
  if (iunit > 1e5) then
    iunit = 0
    exit
  end if
end do
end subroutine

subroutine error(mname,msg)
character(*),intent(in):: mname,msg
integer,parameter:: err_unit = 0
character(*),parameter:: err_fmt = '("ERROR in ",a,": ",a)'
write(err_unit,err_fmt) mname,msg
call error_exit
end subroutine

subroutine warning(mname,msg)
character(*),intent(in):: mname,msg
integer,parameter:: err_unit = 0
character(*),parameter:: err_fmt = '("WARNING in ",a,": ",a)'
write(err_unit,err_fmt) mname,msg
end subroutine

subroutine error_exit
call flush(6)
call abort()
end subroutine

!****************************************************
subroutine build_operator(op_routine,mat,lmat1,lmat2)
!****************************************************
! Subroutine for building explicit matrix corresponding to an implicit operator
      implicit none
      external :: op_routine
      integer,intent(in) :: lmat1, lmat2
      real(kr),intent(out) ::  mat(lmat1,lmat2)

! local vars
      integer  ::             lvec
      real(kr),allocatable ::  vec(:)

      integer  :: i,j

      lvec = lmat1
      allocate (vec(lvec))

      ! allocate columns of identity to apply the operator
      do j = 1,lmat1
         call zero(vec,lvec)

         vec(j) = 1._kr

         call op_routine(vec,lvec)
         ! copy result to matrix
         do i = 1,lmat1
            mat(i,j) = vec(i)
         end do
      end do

      deallocate(vec)
end subroutine

!****************************************************************************
subroutine get_array_union(array1,larray1,array2,larray2,union,lunion,nunion)
!****************************************************************************
! Subroutine for getting union of two integer arrays ordered
      implicit none
! input arrays
      integer,intent(in) :: larray1,         larray2
      integer,intent(in) ::  array1(larray1), array2(larray2)
! output union of arrays 
      integer,intent(in)  :: lunion
      integer,intent(out) ::  union(lunion)
! number of entries in union
      integer,intent(out) ::  nunion

! local vars
      integer :: i, ivalid
      real(kr) :: valid

      ! join arrays to large array
      do i = 1,larray1
         union(i) = array1(i)
      end do
      do i = 1,larray2
         union(larray1 + i) = array2(i)
      end do

      ! sort input arrays
      call iquick_sort(union,lunion)

      ! remove multiplicities
      ivalid = 1
      valid  = union(1)
      do i = 2,lunion
         if (union(i).ne.valid) then
            ivalid = ivalid + 1
            valid = union(i)
            union(ivalid) = valid
         end if
      end do
      nunion = ivalid

      ! pad the rest with zeros
      do i = nunion + 1,lunion
         union(i) = 0
      end do
end subroutine
!********************************************************************************************************
subroutine get_array_intersection(array1,larray1,array2,larray2,intersection,lintersection,nintersection)
!********************************************************************************************************
! Subroutine for getting intersection of two integer arrays ordered
      implicit none
! input arrays
      integer,intent(in) :: larray1,         larray2
      integer,intent(in) ::  array1(larray1), array2(larray2)
! output intersection of arrays 
      integer,intent(in)  :: lintersection
      integer,intent(out) ::  intersection(lintersection)
! number of entries in intersection
      integer,intent(out) ::  nintersection

! local vars
      integer :: workarray1(larray1),workarray2(larray2)
      integer :: i, j, ind1, ind2, iinter

      ! copy arrays to local memory
      do i = 1,larray1
         workarray1(i) = array1(i)
      end do
      do i = 1,larray2
         workarray2(i) = array2(i)
      end do

      ! sort input arrays
      call iquick_sort(workarray1,larray1)
      call iquick_sort(workarray2,larray2)

      ! find common intersection in one loop through the arrays
      i = 1
      j = 1
      iinter = 0
      ind1 = workarray1(i)
      ind2 = workarray2(j)
      do 
         if      (ind1.eq.ind2) then
            ! they match - add to intersection
            iinter = iinter + 1
            intersection(iinter) = ind1
            i = i + 1
            j = j + 1
         else if (ind1.gt.ind2) then
            j = j + 1
         else 
            i = i + 1
         end if
         ! set stopping criterion
         if (i.gt.larray1 .or. j.gt.larray2) then
            exit
         end if
         ind1 = workarray1(i)
         ind2 = workarray2(j)
      end do
      nintersection = iinter
      do iinter = nintersection + 1,lintersection
         intersection(iinter) = 0
      end do
end subroutine

RECURSIVE SUBROUTINE iquick_sort(list,llist)

! Quick sort routine from:
! Brainerd, W.S., Goldberg, C.H. & Adams, J.C. (1990) "Programmer's Guide to
! Fortran 90", McGraw-Hill  ISBN 0-07-000248-7, pages 149-150.
! modified by Jakub Sistek

IMPLICIT NONE
INTEGER,INTENT(IN) :: llist
INTEGER,INTENT(IN OUT)  :: list(llist)

CALL iquick_sort_1(1, llist)

CONTAINS

RECURSIVE SUBROUTINE iquick_sort_1(left_end, right_end)

INTEGER, INTENT(IN) :: left_end, right_end

!     Local variables
INTEGER             :: i, j
INTEGER             :: reference, temp
INTEGER, PARAMETER  :: max_simple_sort_size = 6

IF (right_end < left_end + max_simple_sort_size) THEN
  ! Use interchange sort for small lists
  CALL interchange_sort(left_end, right_end)

ELSE
  ! Use partition ("quick") sort
  reference = list((left_end + right_end)/2)
  i = left_end - 1; j = right_end + 1

  DO
    ! Scan list from left end until element >= reference is found
    DO
      i = i + 1
      IF (list(i) >= reference) EXIT
    END DO
    ! Scan list from right end until element <= reference is found
    DO
      j = j - 1
      IF (list(j) <= reference) EXIT
    END DO


    IF (i < j) THEN
      ! Swap two out-of-order elements
      temp = list(i); list(i) = list(j); list(j) = temp
    ELSE IF (i == j) THEN
      i = i + 1
      EXIT
    ELSE
      EXIT
    END IF
  END DO

  IF (left_end < j) CALL iquick_sort_1(left_end, j)
  IF (i < right_end) CALL iquick_sort_1(i, right_end)
END IF

END SUBROUTINE iquick_sort_1


SUBROUTINE interchange_sort(left_end, right_end)

INTEGER, INTENT(IN) :: left_end, right_end

!     Local variables
INTEGER             :: i, j
INTEGER             :: temp

DO i = left_end, right_end - 1
  DO j = i+1, right_end
    IF (list(i) > list(j)) THEN
      temp = list(i); list(i) = list(j); list(j) = temp
    END IF
  END DO
END DO

END SUBROUTINE interchange_sort

END SUBROUTINE iquick_sort

!*************************************************
subroutine get_index(ivalue,iarray,liarray,iindex)
!*************************************************
! Subroutine for returning first index of element in an array that equals prescribed value
      implicit none
! prescribed value
      integer,intent(in) :: ivalue
! input array
      integer,intent(in) :: liarray
      integer,intent(in) ::  iarray(liarray)
! output index
      integer,intent(out) ::  iindex
! local vars
      integer :: i

      ! find first index where value of element equals prescribed
      iindex = -1
      do i = 1,liarray
         if (iarray(i) .eq. ivalue) then
            iindex = i
            exit
         end if
      end do
end subroutine

!************************************
subroutine condtri(nw,nwx,w,lw, cond)
!************************************
! Routine that estimates condition number of a real symmetric tridiagonal matrix 
! using LAPACK

! used dimension of matrix
      integer,intent(in) :: nw
! allocated size of matrix (leading dimension)
      integer,intent(in) :: nwx

! tridiagonal matrix
      integer,intent(in) :: lw
      real(kr),intent(in) :: w(lw)

! estimate of the condition number
      real(kr),intent(out) :: cond

! local variables

! arrays on stack
      ! diagonal of the matrix
      real(kr) :: d(nw)
      real(kr) :: e(nw-1)

      ! auxiliary variables
      integer ::  iaux
      real(kr) :: raux(1)

      real(kr) :: eigmax, eigmin
      integer ::  i

      ! LAPACK
      integer :: lapack_info

      ! prepare data for LAPACK
      ! diagonal
      do i = 1,nw
         d(i) = w((i-1)*nwx + i)
      end do
      ! subdiagonal
      do i = 1,nw-1
         e(i) = w((i-1)*nwx + i + 1)
      end do

      ! LAPACK routine for searching eigenvalues (returned in ascending order)
      iaux = 1
      call DSTEV( 'N', nw, d, e, raux, iaux, raux, lapack_info)

      ! compute condition number
      eigmax = d(nw)
      eigmin = d(1)
      if (eigmin.lt.1._kr) then
         write(*,*) 'CONDTRI: WARNING - Lowest eigenvalue is less than 1.'
         eigmin = 1._kr
      end if
      if (abs(eigmin).gt.numerical_zero) then
         cond = eigmax/eigmin
      else
         cond = huge(cond)
      end if

end subroutine

end module module_utils

