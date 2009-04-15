module module_utils
! Module for various auxiliary utilities

implicit none
integer,parameter,private :: kr = kind(1.D0)

interface zero
   module procedure zero_int
   module procedure zero_real
end interface zero

contains

subroutine zero_int(ia,lia)
! zero integer array

integer,intent(in)  :: lia
integer,intent(out) ::  ia(lia)
! local
integer :: i

do i = 1,lia
   ia(i) = 0
end do
end subroutine zero_int

subroutine zero_real(x,lx)
! zero real array

integer,intent(in)   :: lx
real(kr),intent(out) ::  x(lx)
! local
integer :: i

do i = 1,lx
   x(i) = 0.0_kr
end do
end subroutine zero_real

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
call abort()
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

end module module_utils

