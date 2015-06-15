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

subroutine lobpcg_mvecmult_f(n,x,lx,y,ly,idoper) bind(c)
! this is an externally callable Fortran routine that is called from C driver
! and just calls the module function for adaptive BDDC
use iso_c_binding
use module_levels 
implicit none
! type of reals
integer,parameter :: kr = kind(1.D0)

! length of vector
integer(c_int),intent(in) ::   n 
integer(c_int),intent(in) ::  lx 
real(c_double),intent(in) ::   x(lx)
integer(c_int),intent(in) ::  ly 
real(c_double),intent(out) ::  y(ly)
integer(c_int),intent(inout) ::  idoper 

call levels_adaptivity_mvecmult(n,x,lx,y,ly,idoper)

end subroutine
