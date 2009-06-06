subroutine lobpcg_mvecmult_f(n,x,lx,y,ly,idoper)
! this is an externally callable Fortran routine that is called from C driver
! and just calls the module function for adaptive BDDC
use module_adaptivity 
implicit none
! type of reals
integer,parameter :: kr = kind(1.D0)

! length of vector
integer,intent(in) ::   n 
integer,intent(in) ::  lx 
real(kr),intent(in) ::  x(lx)
integer,intent(in) ::  ly 
real(kr),intent(out) :: y(ly)
integer,intent(inout) ::  idoper 

call adaptivity_mvecmult(n,x,lx,y,ly,idoper)
end subroutine
