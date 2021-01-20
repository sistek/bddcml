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
use, intrinsic :: iso_fortran_env
use module_levels
implicit none
! type of reals
integer,parameter :: kr = REAL64

! length of vector
integer(c_int),intent(in)             :: n
integer(c_int),intent(in)             :: lx
complex(c_double_complex),intent(in)  :: x(lx)
integer(c_int),intent(in)             :: ly
complex(c_double_complex),intent(out) :: y(ly)
integer(c_int),intent(inout)          :: idoper

call levels_adaptivity_mvecmult(n,x,lx,y,ly,idoper)

end subroutine

subroutine lobpcg_zhegv(itype, jobz, uplo, n, a, lda, b, ldb, w, work, lwork, rwork, info) &
           bind(c,name='lobpcg_zhegv')
! Fortran wrapper of LAPACK ZHEGV for unique symbol name callable from C
use iso_c_binding
implicit none

integer(c_int),intent(in)               :: itype
character(c_char),intent(in)            :: jobz
character(c_char),intent(in)            :: uplo
integer(c_int),intent(in)               :: n
complex(c_double_complex),intent(inout) :: a(lda*n)
integer(c_int),intent(in)               :: lda
complex(c_double_complex),intent(inout) :: b(n)
integer(c_int),intent(in)               :: ldb
real(c_double),intent(out)              :: w
complex(c_double_complex),intent(out)   :: work(lwork)
integer(c_int),intent(in)               :: lwork
real(c_double),intent(in)               :: rwork(3*n-2)
integer(c_int),intent(out)              :: info

call zhegv(itype, jobz, uplo, n, a, lda, b, ldb, w, work, lwork, rwork, info)

end subroutine

subroutine lobpcg_zpotrf(uplo, n, a, lda, info) &
           bind(c,name='lobpcg_zpotrf')
! Fortran wrapper of LAPACK DPOTRF for unique symbol name callable from C
use iso_c_binding
implicit none

character(c_char),intent(in)            :: uplo
integer(c_int),intent(in)               :: n
complex(c_double_complex),intent(inout) :: a(lda*n)
integer(c_int),intent(in)               :: lda
integer(c_int),intent(out)              :: info

call dpotrf(uplo, n, a, lda, info)

end subroutine
