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

subroutine lobpcg_mvecmult_f(x,lx1,lx2, y,ly1,ly2, idoper)
! this is an externally callable Fortran routine that is called from C driver
! and just calls the module function for adaptive BDDC
use module_levels 
implicit none
! type of reals
integer,parameter :: kr = kind(1.D0)

integer,intent(in) ::  lx1,lx2 
real(kr),intent(in) ::  x(lx1,lx2)
integer,intent(in) ::  ly1,ly2 
real(kr),intent(out) :: y(ly1,ly2)
integer,intent(inout) ::  idoper 

call levels_adaptivity_mvecmult(x,lx1,lx2, y,ly1,ly2, idoper)

end subroutine
