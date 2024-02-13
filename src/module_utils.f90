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

module module_utils
! Module for various auxiliary utilities

use, intrinsic :: iso_fortran_env
implicit none
integer,parameter,private :: kr = REAL64

logical,parameter,private :: debug = .false.

! standard UNIX units
logical,private :: suppress_output = .false.
integer,parameter,private:: unit_err    = 0 ! stderr
integer,parameter,private:: unit_stdout = 6 ! stdout

interface info
   module procedure info_plain
   module procedure info_with_number_int
   module procedure info_with_number_rp
end interface info

interface warning
   module procedure warning_plain
   module procedure warning_with_number_int
   module procedure warning_with_number_rp
end interface warning

interface error
   module procedure error_plain
   module procedure error_with_number_int
   module procedure error_with_number_rp
end interface error

interface time_print
   module procedure time_print_plain
   module procedure time_print_with_number_int
end interface time_print

interface rquick_sort
   module procedure rquick_sort_sp
   module procedure rquick_sort_dp
end interface rquick_sort

! Time measurements
integer,parameter,private :: level_time_max = 30
real(kr),private ::          times(level_time_max) = 0._kr
logical,private ::           is_wall_time(level_time_max) 
integer,private  ::          level_time = 0

! data for pseudo-random numbers
logical,private :: initialized_random_numbers = .false.

contains

!***********************************
subroutine integer2string(i,stringi)
!***********************************
!     converts integer i to string
      implicit none
      
! subdomain number
      integer,intent(in) :: i

! name of the problem
      character(*),intent(out) :: stringi 

! local vars
      character(*),parameter:: routine_name = 'INTEGER2STRING'
      integer :: l, j
      character(2) :: fmtstr
      logical :: assigned
      
      ! generate name for
      l = len(stringi)
      do j = 1,l
         stringi(j:j) = '0'
      end do
      assigned = .false.
      do j = 1,l
         write(fmtstr,'(i2)') j
         if (i.lt.10**j) then
            write(stringi(l-j+1:l),'(i'//fmtstr//')') i
            assigned = .true.
            exit
         end if
      end do
      if (.not.assigned) then
         call error( routine_name, 'Index too large for assignement.', i)
      end if

end subroutine

!************************************
subroutine integer2logical(tf_int,tf)
!************************************
! translates C-like logical value as integer 0 = false, 1 = true into Fortran logical type
implicit none
integer, intent(in) :: tf_int
logical, intent(out) :: tf

! local vars
character(*),parameter:: routine_name = 'INTEGER2LOGICAL'

if      (tf_int.eq.0) then
   tf = .false.
else if (tf_int.eq.1) then
   tf = .true.
else
   call error(routine_name, 'Illegal value of integer true/false type - not 0 nor 1 but ', tf_int)
end if

end subroutine

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
      character(5) :: numberstr

! zero filename
      fname = ' '

      ! generate name for 4 digits
      call integer2string(isub, numberstr)

      fname = trim(problemname)//'/'//trim(problemname)//'_'//trim(numberstr)//'.'//trim(suffix)

end subroutine

!******************************
subroutine allocate_unit(iunit)
!******************************
! find an empty IO unit (originally written by Jaroslav Hajek)
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

!****************************
subroutine suppress_output_on
!****************************
! set output suppressing flag on
suppress_output = .true.
end subroutine suppress_output_on

!*****************************
subroutine suppress_output_off
!*****************************
! set output suppressing flag off
suppress_output = .false.
end subroutine suppress_output_off

!*******************************
subroutine info_plain(mname,msg)
!*******************************
character(*),intent(in):: mname,msg
character(*),parameter:: info_fmt = '(a,": ",a)'
if ( .not. suppress_output ) then
   write(unit_stdout,info_fmt) mname,msg
   call flush(unit_stdout)
end if
end subroutine

!*************************************************
subroutine info_with_number_int(mname,msg,inumber)
!*************************************************
character(*),intent(in):: mname,msg
integer,intent(in) :: inumber
character(*),parameter:: info_fmt = '(a,": ",a,X,i10)'
if ( .not. suppress_output ) then
   write(unit_stdout,info_fmt) mname,msg,inumber
   call flush(unit_stdout)
end if
end subroutine

!************************************************
subroutine info_with_number_rp(mname,msg,rnumber)
!************************************************
character(*),intent(in):: mname,msg
real(kr),intent(in) :: rnumber
character(*),parameter:: info_fmt = '(a,": ",a,X,e17.9)'
if ( .not. suppress_output ) then
   write(unit_stdout,info_fmt) mname,msg,rnumber
   call flush(unit_stdout)
end if
end subroutine

!**********************************
subroutine warning_plain(mname,msg)
!**********************************
character(*),intent(in):: mname,msg
character(*),parameter:: wrn_fmt = '("WARNING in ",a,": ",a)'
if ( .not. suppress_output ) then
   write(unit_stdout,wrn_fmt) mname,msg
   call flush(unit_stdout)
end if
end subroutine

!****************************************************
subroutine warning_with_number_int(mname,msg,inumber)
!****************************************************
character(*),intent(in):: mname,msg
integer,intent(in) :: inumber
character(*),parameter:: wrn_fmt = '("WARNING in ",a,": ",a,X,i10)'
if ( .not. suppress_output ) then
   write(unit_stdout,wrn_fmt) mname,msg,inumber
   call flush(unit_stdout)
end if
end subroutine

!***************************************************
subroutine warning_with_number_rp(mname,msg,rnumber)
!***************************************************
character(*),intent(in):: mname,msg
real(kr),intent(in) :: rnumber
character(*),parameter:: wrn_fmt = '("WARNING in ",a,": ",a,X,e17.9)'
if ( .not. suppress_output ) then
   write(unit_stdout,wrn_fmt) mname,msg,rnumber
   call flush(unit_stdout)
end if
end subroutine

!********************************
subroutine error_plain(mname,msg)
!********************************
character(*),intent(in):: mname,msg
character(*),parameter:: err_fmt = '("ERROR in ",a,": ",a)'
write(unit_err,err_fmt) mname,msg
call flush(unit_err)
call error_exit
end subroutine

!**************************************************
subroutine error_with_number_int(mname,msg,inumber)
!**************************************************
character(*),intent(in):: mname,msg
integer,intent(in) :: inumber
character(*),parameter:: err_fmt = '("ERROR in ",a,": ",a,X,i10)'
write(unit_err,err_fmt) mname,msg,inumber
call flush(unit_err)
call error_exit
end subroutine

!*************************************************
subroutine error_with_number_rp(mname,msg,rnumber)
!*************************************************
character(*),intent(in):: mname,msg
real(kr),intent(in) :: rnumber
character(*),parameter:: err_fmt = '("ERROR in ",a,": ",a,X,e17.9)'
write(unit_err,err_fmt) mname,msg,rnumber
call flush(unit_err)
call error_exit
end subroutine

subroutine error_exit
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
         vec(:) = 0._kr

         vec(j) = 1._kr

         call op_routine(vec,lvec)
         ! copy result to matrix
         do i = 1,lmat1
            mat(i,j) = vec(i)
         end do
      end do

      deallocate(vec)
end subroutine

!***********************************************************
subroutine write_matrix(idunit,mat,format_string,info,trans)
!***********************************************************
! write matrix to unit
! arguments:
! idunit        the unit to write to
! mat           the matrix to write
! format_string (optional) the format to use, e.g. 'E16.4'. default is E23.16
! info          (optional) the exit status
! trans         (optional) write transposed matrix
implicit none
integer,intent(in):: idunit
real(kr),intent(in):: mat(:,:)
character(*),intent(in),optional:: format_string
integer,intent(out),optional:: info
logical,intent(in),optional:: trans
character(50):: tfmt
logical:: ttrans
integer:: ios
integer:: i

! setup format
if (present(format_string)) then
  tfmt = '(    (' // format_string //',1x))'
else
  tfmt = '(    (' // 'E23.16' //',1x))'
end if
ttrans = .false.
if (present(trans)) ttrans = trans
if (ttrans) then
  write(tfmt(2:5),'(I4)') size(mat,1)
  do i = 1,size(mat,1)
    write(idunit,tfmt,iostat=ios) mat(:,i)
    if (ios /= 0) exit
  end do
else
  write(tfmt(2:5),'(I4)') size(mat,2)
  do i = 1,size(mat,1)
    write(idunit,tfmt,iostat=ios) mat(i,:)
    if (ios /= 0) exit
  end do
end if
if (present(info)) info = ios

end subroutine


!*******************************************
recursive subroutine iquick_sort(list,llist)
!*******************************************
! Intiger quick sort routine from:
! Brainerd, W.S., Goldberg, C.H. & Adams, J.C. (1990) "Programmer's Guide to
! Fortran 90", McGraw-Hill  ISBN 0-07-000248-7, pages 149-150.
! modified by Jakub Sistek

implicit none
integer,intent(in) :: llist
integer,intent(inout)  :: list(llist)

call iquick_sort_1(1, llist)

contains

recursive subroutine iquick_sort_1(left_end, right_end)

integer, intent(in) :: left_end, right_end

!     Local variables
integer             :: i, j
integer             :: reference, temp
integer, parameter  :: max_simple_sort_size = 6

if (right_end < left_end + max_simple_sort_size) then
  ! use interchange sort for small lists
  call iinterchange_sort(left_end, right_end)

else
  ! use partition ("quick") sort
  reference = list((left_end + right_end)/2)
  i = left_end - 1; j = right_end + 1

  do
    ! scan list from left end until element >= reference is found
    do
      i = i + 1
      if (list(i) >= reference) exit
    end do
    ! scan list from right end until element <= reference is found
    do
      j = j - 1
      if (list(j) <= reference) exit
    end do


    if (i < j) then
      ! swap two out-of-order elements
      temp = list(i); list(i) = list(j); list(j) = temp
    else if (i == j) then
      i = i + 1
      exit
    else
      exit
    end if
  end do

  if (left_end < j) call iquick_sort_1(left_end, j)
  if (i < right_end) call iquick_sort_1(i, right_end)
end if

end subroutine iquick_sort_1


subroutine iinterchange_sort(left_end, right_end)

integer, intent(in) :: left_end, right_end

!     local variables
integer             :: i, j
integer             :: temp

do i = left_end, right_end - 1
  do j = i+1, right_end
    if (list(i) > list(j)) then
      temp = list(i); list(i) = list(j); list(j) = temp
    end if
  end do
end do

end subroutine iinterchange_sort
end subroutine iquick_sort

!***********************************************************************
recursive subroutine iquick_sort_simultaneous(list1,llist1,list2,llist2)
!***********************************************************************
! quick sort that sorts two arrays based on sorting of first array
! Quick sort routine from:
! Brainerd, W.S., Goldberg, C.H. & Adams, J.C. (1990) "Programmer's Guide to
! Fortran 90", McGraw-Hill  ISBN 0-07-000248-7, pages 149-150.
! modified by Jakub Sistek

implicit none
integer,intent(in) :: llist1
integer,intent(in out)  :: list1(llist1)
integer,intent(in) :: llist2
integer,intent(in out)  :: list2(llist2)

! local vars
character(*),parameter:: routine_name = 'IQUICK_SORT_SIMULTANEOUS'

! check dimension
if (llist1.ne.llist2) then
   call error(routine_name,'dimension of arrays mismatch')
end if

call iquick_sort_simultaneous_1(1, llist1)

contains

recursive subroutine iquick_sort_simultaneous_1(left_end, right_end)

integer, intent(in) :: left_end, right_end

!     local variables
integer             :: i, j
integer             :: reference, temp
integer, parameter  :: max_simple_sort_size = 6

if (right_end < left_end + max_simple_sort_size) then
  ! use interchange sort for small lists
  call iinterchange_sort_simultaneous(left_end, right_end)

else
  ! use partition ("quick") sort
  reference = list1((left_end + right_end)/2)
  i = left_end  - 1 
  j = right_end + 1

  do
    ! scan list from left end until element >= reference is found
    do
      i = i + 1
      if (list1(i) >= reference) exit
    end do
    ! scan list from right end until element <= reference is found
    do
      j = j - 1
      if (list1(j) <= reference) exit
    end do


    if (i < j) then
      ! swap two out-of-order elements
      temp = list1(i) ; list1(i) = list1(j) ; list1(j) = temp
      temp = list2(i) ; list2(i) = list2(j) ; list2(j) = temp
    else if (i == j) then
      i = i + 1
      exit
    else
      exit
    end if
  end do

  if (left_end < j) call iquick_sort_simultaneous_1(left_end, j)
  if (i < right_end) call iquick_sort_simultaneous_1(i, right_end)
end if

end subroutine iquick_sort_simultaneous_1

subroutine iinterchange_sort_simultaneous(left_end, right_end)

integer, intent(in) :: left_end, right_end

!     local variables
integer             :: i, j
integer             :: temp

do i = left_end, right_end - 1
  do j = i+1, right_end
    if (list1(i) > list1(j)) then
      temp = list1(i); list1(i) = list1(j); list1(j) = temp
      temp = list2(i); list2(i) = list2(j); list2(j) = temp
    end if
  end do
end do

end subroutine iinterchange_sort_simultaneous

end subroutine iquick_sort_simultaneous

!***************************************************************************************
recursive subroutine iquick_array_sort_simultaneous(array1,larray1,larray2,list2,llist2)
!***************************************************************************************
! quick sort that sorts two arrays based on sorting of first array
! Quick sort routine from:
! Brainerd, W.S., Goldberg, C.H. & Adams, J.C. (1990) "Programmer's Guide to
! Fortran 90", McGraw-Hill  ISBN 0-07-000248-7, pages 149-150.
! modified by Jakub Sistek

implicit none
integer,intent(in) :: larray1
integer,intent(in) :: larray2
integer,intent(inout)  :: array1(larray1,larray2)
integer,intent(in) :: llist2
integer,intent(inout)  :: list2(llist2)

! local vars
character(*),parameter:: routine_name = 'IQUICK_ARRAY_SORT_SIMULTANEOUS'

! check dimension
if (larray1.ne.llist2) then
   call error(routine_name,'dimension of arrays mismatch')
end if

CALL iquick_array_sort_simultaneous_1(1, larray1)

contains

recursive subroutine iquick_array_sort_simultaneous_1(left_end, right_end)

integer, intent(in) :: left_end, right_end

!     local variables
integer             :: i, j
integer             :: temp1
integer             :: reference(larray2), temp(larray2)
integer, parameter  :: max_simple_sort_size = 6

if (right_end < left_end + max_simple_sort_size) then
  ! use interchange sort for small lists
  call iinterchange_array_sort_simultaneous(left_end, right_end)

else
  ! Use partition ("quick") sort
  reference = array1((left_end + right_end)/2,:)
  i = left_end  - 1 
  j = right_end + 1

  do
    ! scan list from left end until element >= reference is found
    do
      i = i + 1
      if ( .not. is_less_than(array1(i,:),reference,larray2) ) exit
    end do
    ! scan list from right end until element <= reference is found
    do
      j = j - 1
      if ( .not. is_less_than(reference,array1(j,:),larray2) ) exit
    end do

    if (i < j) then
      ! swap two out-of-order elements
      temp  = array1(i,:) ; array1(i,:) = array1(j,:) ; array1(j,:) = temp
      temp1 = list2(i) ;    list2(i) = list2(j) ;       list2(j)   = temp1
    else if (i == j) then
      i = i + 1
      exit
    else
      exit
    end if
  end do

  if (left_end < j) call iquick_array_sort_simultaneous_1(left_end, j)
  if (i < right_end) call iquick_array_sort_simultaneous_1(i, right_end)
end if

end subroutine iquick_array_sort_simultaneous_1

subroutine iinterchange_array_sort_simultaneous(left_end, right_end)

integer, intent(in) :: left_end, right_end

!     local variables
integer             :: i, j
integer             :: temp1
integer             :: temp(larray2)

do i = left_end, right_end - 1
  do j = i+1, right_end
    if ( is_less_than(array1(j,:),array1(i,:),larray2) ) then
      temp  = array1(i,:); array1(i,:) = array1(j,:); array1(j,:) = temp
      temp1 = list2(i); list2(i) = list2(j); list2(j) = temp1
    end if
  end do
end do

end subroutine iinterchange_array_sort_simultaneous

function is_less_than(array1,array2,length) result(ilt)
implicit none
integer,intent(in) :: length
integer,intent(in) :: array1(length)
integer,intent(in) :: array2(length)
logical :: ilt

integer i

ilt = .false.
do i = 1,length
   if   (array1(i) < array2(i)) then
      ilt = .true.
      exit
   else if (array2(i) < array1(i)) then
      ilt = .false.
      exit
   else 
      ! element is equal at this point
      cycle
   end if
end do

end function
end subroutine iquick_array_sort_simultaneous


!**************************************************
subroutine iquick_sort_2(list1,llist1,list2,llist2)
!**************************************************
! routine to sort array of two columns first by first index and then by second index
implicit none
integer,intent(in) ::   llist1
integer,intent(inout) :: list1(llist1)
integer,intent(in) ::   llist2
integer,intent(inout) :: list2(llist2)
! local vars
character(*),parameter:: routine_name = 'IQUICK_SORT_2'
integer :: ibegin, iend, indbegin, indend
integer :: nvals

! check dimension
if (llist1.ne.llist2) then
   call error(routine_name,'dimension of arrays mismatch')
end if

! sort the array first by first index
call iquick_sort_simultaneous(list1,llist1,list2,llist2)

! Sort by second index
ibegin = 1
iend   = ibegin
do while (iend.le.llist1)
   indbegin = list1(ibegin)
   indend   = list1(iend)
   ! find the index of last same value
   do while (indend.eq.indbegin)
      iend = iend + 1
      ! correction at the end of the array
      if (iend.gt.llist1) exit
      indend = list1(iend)
   end do
   iend = iend - 1
      
   nvals = iend-ibegin + 1
   call iquick_sort(list2(ibegin:iend),nvals)
      
   ibegin = iend + 1
   iend   = ibegin
end do

end subroutine iquick_sort_2

!**************************************************************************
recursive subroutine rquick_sort_sp(list,llist1,llist2,indsort,istart,iend)
!**************************************************************************
! Quick sort routine from:
! Brainerd, W.S., Goldberg, C.H. & Adams, J.C. (1990) "Programmer's Guide to
! Fortran 90", McGraw-Hill  ISBN 0-07-000248-7, pages 149-150.
! modified by Jakub Sistek
! sort two-dimensional real array according to indsort-th column

implicit none
integer,intent(in) ::      llist1, llist2
real(4),intent(inout) ::  list(llist1,llist2)
integer,intent(in) ::      indsort
integer,intent(in) ::      istart, iend

call rquick_sort_1_sp(istart, iend)

contains

recursive subroutine rquick_sort_1_sp(left_end, right_end)

integer, intent(in) :: left_end, right_end

!     local variables
integer             :: i, j, k
real(4)            :: reference, temp
integer, parameter  :: max_simple_sort_size = 6

if (right_end < left_end + max_simple_sort_size) then
  ! use interchange sort for small lists
  call rinterchange_sort_sp(left_end, right_end)

else
  ! use partition ("quick") sort
  reference = list((left_end + right_end)/2,indsort)
  i = left_end  - 1 
  j = right_end + 1

  do
    ! scan list from left end until element >= reference is found
    do
      i = i + 1
      if (list(i,indsort) >= reference) exit
    end do
    ! scan list from right end until element <= reference is found
    do
      j = j - 1
      if (list(j,indsort) <= reference) exit
    end do


    if (i < j) then
      ! swap two out-of-order elements
      do k = 1,llist2
         temp = list(i,k) 
         list(i,k) = list(j,k) 
         list(j,k) = temp
      end do
    else if (i == j) then
      i = i + 1
      exit
    else
      exit
    end if
  end do

  if (left_end < j) call rquick_sort_1_sp(left_end, j)
  if (i < right_end) call rquick_sort_1_sp(i, right_end)
end if

end subroutine rquick_sort_1_sp

subroutine rinterchange_sort_sp(left_end, right_end)

integer, intent(in) :: left_end, right_end

!     local variables
integer             :: i, j, k
real(4)            :: temp

do i = left_end, right_end - 1
  do j = i+1, right_end
    if (list(i,indsort) > list(j,indsort)) then
      ! swap two out-of-order rows
      do k = 1,llist2
         temp = list(i,k) 
         list(i,k) = list(j,k) 
         list(j,k) = temp
      end do
    end if
  end do
end do

end subroutine rinterchange_sort_sp
end subroutine rquick_sort_sp

!**************************************************************************
recursive subroutine rquick_sort_dp(list,llist1,llist2,indsort,istart,iend)
!**************************************************************************
! Quick sort routine from:
! Brainerd, W.S., Goldberg, C.H. & Adams, J.C. (1990) "Programmer's Guide to
! Fortran 90", McGraw-Hill  ISBN 0-07-000248-7, pages 149-150.
! modified by Jakub Sistek
! sort two-dimensional real array according to indsort-th column

implicit none
integer,intent(in) ::      llist1, llist2
real(8),intent(inout) ::  list(llist1,llist2)
integer,intent(in) ::      indsort
integer,intent(in) ::      istart, iend

call rquick_sort_1_dp(istart, iend)

contains

recursive subroutine rquick_sort_1_dp(left_end, right_end)

integer, intent(in) :: left_end, right_end

!     local variables
integer             :: i, j, k
real(8)            :: reference, temp
integer, parameter  :: max_simple_sort_size = 6

if (right_end < left_end + max_simple_sort_size) then
  ! use interchange sort for small lists
  call rinterchange_sort_dp(left_end, right_end)

else
  ! use partition ("quick") sort
  reference = list((left_end + right_end)/2,indsort)
  i = left_end  - 1 
  j = right_end + 1

  do
    ! scan list from left end until element >= reference is found
    do
      i = i + 1
      if (list(i,indsort) >= reference) exit
    end do
    ! scan list from right end until element <= reference is found
    do
      j = j - 1
      if (list(j,indsort) <= reference) exit
    end do


    if (i < j) then
      ! swap two out-of-order elements
      do k = 1,llist2
         temp = list(i,k) 
         list(i,k) = list(j,k) 
         list(j,k) = temp
      end do
    else if (i == j) then
      i = i + 1
      exit
    else
      exit
    end if
  end do

  if (left_end < j) call rquick_sort_1_dp(left_end, j)
  if (i < right_end) call rquick_sort_1_dp(i, right_end)
end if

end subroutine rquick_sort_1_dp

subroutine rinterchange_sort_dp(left_end, right_end)

integer, intent(in) :: left_end, right_end

!     local variables
integer             :: i, j, k
real(8)            :: temp

do i = left_end, right_end - 1
  do j = i+1, right_end
    if (list(i,indsort) > list(j,indsort)) then
      ! swap two out-of-order rows
      do k = 1,llist2
         temp = list(i,k) 
         list(i,k) = list(j,k) 
         list(j,k) = temp
      end do
    end if
  end do
end do

end subroutine rinterchange_sort_dp
end subroutine rquick_sort_dp

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
      integer :: valid

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

      ! shortcut for degenerate data
      if (larray1.eq.0 .or. larray2.eq.0) then
         intersection = 0
         nintersection = 0
         return
      end if

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

!***************************************************
subroutine get_array_norepeat(array,larray,nentries)
!***************************************************
! Subroutine for removing repeated entries from SORTED array
! shifting entries at beginning, pad with zeros
! Example
! input:  1 1 1 3 3 5 9
! output: 1 3 5 9 0 0 0
! nentries  = 4
      implicit none
! input arrays
      integer,intent(in) ::    larray
      integer,intent(inout) ::  array(larray)
! number of entries
      integer,intent(out) ::   nentries

! local vars
      integer :: i, j, ind

! return if array has zero length
      if (larray .eq. 0) then
         nentries = 0
         return
      end if

      ! find common intersection in one loop through the arrays
      i = 1
      j = 1
      ind = array(i)
      if (ind.eq.0) then
         write(*,*) 'GET_ARRAY_NOREPEAT: Error - array starts with zero!'
         call error_exit
      end if
      do j = 2,larray
         if      (array(j).gt.ind) then
            i   = i + 1
            ind = array(j)
            array(i) = ind
         else 
            if (debug .and. array(j).lt.ind) then
               write(*,*) 'GET_ARRAY_NOREPEAT: Error - array supposed sorted!'
               call error_exit
            end if
         end if
      end do
      ! pad the rest
      nentries = i
      do i = nentries + 1,larray
         array(i) = 0
      end do
end subroutine

!*********************************************************************************
subroutine get_array_norepeat_simultaneous(array1,larray1,array2,larray2,nentries)
!*********************************************************************************
! Subroutine for removing repeated entries from SORTED array ARRAY1
! shifting entries in both ARRAYS at beginning, pad with zeros
! Example
! input:  1 1 1 3 3 5 9
! output: 1 3 5 9 0 0 0
! nentries  = 4
      implicit none
! input arrays
      integer,intent(in) ::    larray1
      integer,intent(inout) ::  array1(larray1)
      integer,intent(in) ::    larray2
      integer,intent(inout) ::  array2(larray2)
! number of entries
      integer,intent(out) ::   nentries

! local vars
      character(*),parameter:: routine_name = 'GET_ARRAY_NOREPEAT_SIMULTANEOUS'
      integer :: i, j, ind

! return if array has zero length
      if (larray1 .eq. 0) then
         nentries = 0
         return
      end if
      ! check
      if (larray1 .ne. larray2) then
         call error(routine_name,'dimension mismatch')
      end if

      ! pack both arrays
      i = 1
      j = 1
      ind = array1(i)
      if (ind.eq.0) then
         call error(routine_name, 'array1 starts with zero')
      end if
      do j = 2,larray1
         if      (array1(j).gt.ind) then
            i   = i + 1
            ind = array1(j)
            array1(i) = ind
            array2(i) = array2(j)
         else 
            if (debug .and. array1(j).lt.ind) then
               call error(routine_name, 'array1 supposed to be sorted')
            end if
         end if
      end do
      ! pad the rest
      nentries = i
      do i = nentries + 1,larray1
         array1(i) = 0
         array2(i) = 0
      end do
end subroutine

!**********************************************************************
subroutine get_array_norepeat_2(array1,larray1,array2,larray2,nentries)
!**********************************************************************
! Subroutine for removing repeated pairs entries from SORTED arrays ARRAY1 and ARRAY2
! shifting entries in both ARRAYS at beginning, pad with zeros
! Example
! input1:  1 1 1 3 3 5 9 9
! input2:  1 2 2 1 2 1 1 1
! output1: 1 1 3 3 5 9 0 0 
! output2: 1 2 1 2 1 1 0 0 
! nentries  = 6
      implicit none
! input arrays
      integer,intent(in) ::    larray1
      integer,intent(inout) ::  array1(larray1)
      integer,intent(in) ::    larray2
      integer,intent(inout) ::  array2(larray2)
! number of entries
      integer,intent(out) ::   nentries

! local vars
      character(*),parameter:: routine_name = 'GET_ARRAY_NOREPEAT_2'
      integer :: i, j, ind1, ind2

! return if array has zero length
      if (larray1 .eq. 0) then
         nentries = 0
         return
      end if
      ! check
      if (larray1 .ne. larray2) then
         call error(routine_name,'dimension mismatch')
      end if

      ! pack both arrays
      i = 1
      j = 1
      ind1 = array1(i)
      ind2 = array2(i)
      if (ind1.eq.0) then
         call error(routine_name, 'array1 starts with zero')
      end if
      if (ind2.eq.0) then
         call error(routine_name, 'array2 starts with zero')
      end if
      do j = 2,larray1
         ! core part
         if      (array1(j).gt.ind1 .or. array2(j).gt.ind2) then
            i   = i + 1
            ind1 = array1(j)
            ind2 = array2(j)
            array1(i) = ind1
            array2(i) = ind2
         else 
            if (debug .and. array1(j).lt.ind1) then
               call error(routine_name, 'array1 supposed to be sorted')
            end if
         end if
      end do
      ! pad the rest
      nentries = i
      do i = nentries + 1,larray1
         array1(i) = 0
         array2(i) = 0
      end do
end subroutine

!*********************************
subroutine push_back(iarray, ival)
!*********************************
! an analog to C++ push_back member function of std::vector
! makes array iarray larger by one and adds value ival at the end of the array
implicit none
integer,allocatable,intent(inout) :: iarray(:)
integer,intent(in) :: ival

! local vars
integer,allocatable :: tmp_arr(:)

allocate(tmp_arr(size(iarray)+1))
tmp_arr(1:size(iarray)) = iarray
tmp_arr(size(iarray)+1) = ival
deallocate(iarray)
allocate(iarray(size(tmp_arr)))
call move_alloc(tmp_arr, iarray)
end subroutine

!*************************************
subroutine counts2starts(array,larray)
!*************************************
! routine that converts counts to starts
! e.g. 
! input:  | 3 4 2  2  2  X | 
! output: | 1 4 8 10 12 14 |
implicit none
integer, intent(in) ::   larray
integer, intent(inout) :: array(larray)
! local vars
integer :: i

do i = 2,larray
   array(i) = array(i-1) + array(i)
end do
! shift it one back and add one 
do i = larray, 2, - 1 
   array(i) = array(i-1) + 1
end do
! put one in the beginning
array(1) = 1
end subroutine

!********************************************************
subroutine get_index_sorted(ivalue,iarray,liarray,iindex)
!********************************************************
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
      integer :: left_end, right_end

      iindex = -1
      left_end  = 1
      right_end = liarray

      call get_index_sorted_1

      contains

      recursive subroutine get_index_sorted_1
      implicit none
    
      !     Local variables
      character(*),parameter:: routine_name = 'GET_INDEX_SORTED_1'
      integer             :: indnew, valref
      integer, parameter  :: max_simple_search_size = 10 
      integer :: iindex_loc
    
      if (right_end < left_end + max_simple_search_size) then
        ! Use interchange sort for small lists
        call get_index(ivalue,iarray(left_end:right_end),right_end - left_end + 1,iindex_loc)
        if (iindex_loc.lt.0) then
           !print *,'left_end, right_end', left_end, right_end
           !print *,'iarray', iarray(left_end:right_end)
           call warning(routine_name,'index does not appear to be in list')
           iindex = -1
        else
           iindex = left_end + iindex_loc - 1
        end if
      else
        ! Set new limits
        indnew = (left_end + right_end)/2
        valref = iarray(indnew)
        if (ivalue.gt.valref) then
           left_end = indnew
        else if (ivalue.le.valref) then
           right_end = indnew
        else 
           call error(routine_name,'unable to split interval')
        end if
      
        call get_index_sorted_1
      end if
      end subroutine
end subroutine

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

!*************************************************
logical elemental pure function fuzzyLessThan(x,y)
!*************************************************
! function for comparison of two reals with some tolerance
! x - eps < y

real(kr),intent(in) :: x,y
! local vars
real(kr),parameter :: overlap = 1.e-8

fuzzyLessThan = x+overlap .lt. y 

end function

!**********************************
subroutine initialize_random_number
!**********************************
! subroutine for initialization of random number generation
implicit none
integer :: lseed
integer,allocatable :: seed(:)
integer :: i

      ! find size of seed
      call random_seed(size = lseed)

      ! generate new seed
      allocate(seed(lseed))
      ! simply make a sequence 3,3,3,...
      do i = 1,lseed
         seed(i) = 3
      end do

      call random_seed(put = seed(1:lseed))
      deallocate(seed)

      initialized_random_numbers = .true.

end subroutine

!******************************
subroutine get_random_number(x)
!******************************
! subroutine for random number generation
implicit none
real(kr), intent(out) :: x

      ! in the first call of the routine, initialize the random number generator, then use the sequence
      if (.not.initialized_random_numbers) then
         if (debug) then
            write (*,*) 'Initializing random number generator.'
         end if
         call initialize_random_number
      else
         if (debug) then
            write (*,*) 'Using initialized random number generator.'
         end if
      end if

      call random_number(x)

end subroutine

!**********************************
subroutine time_start(use_cpu_time)
!**********************************
! Routine that starts new time measurement
! Levels of timing work like opening and closing brackets,
! measuring time elapsed between a pair of them.
! This routine is like opening a bracket, 
! see routine TIME_END for the opposite.
      implicit none
      logical, intent(in), optional :: use_cpu_time
      
! Local variables
      character(*),parameter:: routine_name = 'TIME_START'
      logical :: uct = .false.
      
      uct = .false.
      if (present(use_cpu_time)) then
         if (use_cpu_time) then
            uct = .true.
         end if
      end if

! add new level of timing
      level_time = level_time + 1
      if (debug) then
         call info(routine_name,'Starting time level ',level_time)
      end if

! check if it is not too many
      if (level_time.gt.level_time_max) then
         write(*,*) 'TIME_START: Maximal number of time levels reached.'
         stop
      end if

! measure the time and add it to the times array
      if (uct) then
         ! use CPU time
         call processor_time(times(level_time))
         is_wall_time(level_time) = .false.
      else
         call wall_time(times(level_time))
         is_wall_time(level_time) = .true.
      end if

      return
end subroutine

!************************
subroutine time_end(time)
!************************
! Routine that finish time measurement of a level
! Levels of timing work like opening and closing brackets,
! measuring time elapsed between a pair of them.
! This routine is like closing a bracket, 
! see routine TIME_START for the opposite.

      implicit none

! Time of run
      real(kr),intent(out) :: time
      
! Local variables
      character(*),parameter:: routine_name = 'TIME_END'
      real(kr) :: current_time

! check if it is not too few
      if (level_time.le.0) then
         write(*,*) 'TIME_END: All time levels already finished.'
         stop
      end if

! measure the time
      if (is_wall_time(level_time)) then
         call wall_time(current_time)
      else
         call processor_time(current_time)
      end if

! find the elapsed time
      time =  current_time - times(level_time)

! subtract one level of timing
      if (debug) then
         call info(routine_name,'Closing time level ',level_time)
      end if
      level_time = level_time - 1

      return
end subroutine

!************************************
subroutine time_print_plain(msg,time)
!************************************
! Routine that prints info about time
implicit none
character(*),intent(in):: msg
real(kr),intent(in)    :: time
character(*),parameter:: info_fmt = '("Time of ",a,": ",f13.3," s")'
if ( .not. suppress_output ) then
   write(unit_stdout,'(a)') '****PROFILING*****************************'
   write(unit_stdout,info_fmt) msg,time
   write(unit_stdout,'(a)') '****PROFILING*****************************'
   call flush(unit_stdout)
end if
end subroutine

!**************************************************
subroutine time_print_with_number_int(msg,num,time)
!**************************************************
! Routine that prints info about time
implicit none
character(*),intent(in):: msg
integer,intent(in)     :: num
real(kr),intent(in)    :: time
character(*),parameter:: info_fmt = '("Time of ",a," ",i5,": ",f13.3," s")'
if ( .not. suppress_output ) then
   write(unit_stdout,'(a)') '****PROFILING*****************************'
   write(unit_stdout,info_fmt) msg,num,time
   write(unit_stdout,'(a)') '****PROFILING*****************************'
   call flush(unit_stdout)
end if
end subroutine

!**********************
subroutine wall_time(t)
!**********************
! return wall clock time in s
implicit none
real(kr) :: t
integer, dimension(8):: v
integer(8) :: timems
!integer count,rate,cmax
!call system_clock(count,rate,cmax)
!t=dble(count)/dble(rate)
call date_and_time(values=v) 
timems=(v(8)+1000*(v(7)+60*(v(6)+60*(v(5)+24*(v(3))))))
t= timems / 1000._kr
end subroutine wall_time

!***************************
subroutine processor_time(t)
!***************************
! return CPU time in s
implicit none
real(kr) :: t
call cpu_time(t) 
end subroutine processor_time


end module module_utils

