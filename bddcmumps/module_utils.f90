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


end module module_utils

