module module_errors
! this module is based on Jaroslav Hajek's auxiliary routines
implicit none
private
public assert,err_unit,error,warning
integer,parameter:: err_unit = 0
logical,save:: interrupted = .false.

interface error
  module procedure error1
  module procedure error2
  module procedure error3
end interface

interface warning
  module procedure warning1
  module procedure warning2
  module procedure warning3
end interface

character(*),parameter:: err_fmt = '("Error in ",a,": ",a)'
character(*),parameter:: warn_fmt = '("Warning from ",a,": ",a)'

contains

subroutine assert(cond,modname,i)
logical,intent(in):: cond
character(*),intent(in):: modname
integer,intent(in):: i
if (.not. cond) then
  write(err_unit,'("assertion #",I1," failed in module ",a)') i,modname
end if
end subroutine

subroutine wait_input
integer:: ios
read(*,*,iostat=ios) 
end subroutine

subroutine error1(mname,msg)
character(*),intent(in):: mname,msg
write(err_unit,err_fmt) mname,msg
call error_exit
end subroutine

subroutine error2(mname,msg,arg1)
character(*),intent(in):: mname,msg,arg1
write(err_unit,err_fmt) mname,msg//' '//arg1
call error_exit
end subroutine

subroutine error3(mname,msg,arg1,arg2)
character(*),intent(in):: mname,msg,arg1,arg2
write(err_unit,err_fmt) mname,msg//' '//arg1//' '//arg2
call error_exit
end subroutine

subroutine warning1(mname,msg)
character(*),intent(in):: mname,msg
write(err_unit,warn_fmt) mname,msg
end subroutine

subroutine warning2(mname,msg,arg1)
character(*),intent(in):: mname,msg,arg1
write(err_unit,warn_fmt) mname,msg//' '//arg1
end subroutine

subroutine warning3(mname,msg,arg1,arg2)
character(*),intent(in):: mname,msg,arg1,arg2
write(err_unit,warn_fmt) mname,msg//' '//arg1//' '//arg2
end subroutine

subroutine error_exit
stop
end subroutine

end module

