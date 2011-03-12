! BDDCML - Multilevel BDDC
! Copyright (C) The BDDCML Team
! 
! This program is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
! 
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
! 
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <http://www.gnu.org/licenses/>.
! ____________________________________________________________________

subroutine lobpcg_mvecmult_f(n,x,lx,y,ly,idoper)
! this is an externally callable Fortran routine that is called from C driver
! and just calls the module function for adaptive BDDC
use module_levels 
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

call levels_adaptivity_mvecmult(n,x,lx,y,ly,idoper)

end subroutine
