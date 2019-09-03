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

module module_densela
!********************
! Module for interfacing dense numerical libraries, like LAPACK,
! PLASMA, MAGMA etc.

      use, intrinsic :: iso_fortran_env
#if defined(BDDCML_WITH_MAGMA)
      use magma
#endif
      implicit none

! type of real variables
      integer,parameter,private :: kr = REAL64

! debugging mode
      logical,parameter,private :: debug = .false.

! library to use
      integer,parameter :: DENSELA_LAPACK = 1
      integer,parameter :: DENSELA_MAGMA  = 2

      !integer,parameter,private :: library = DENSELA_MAGMA
      !integer,parameter,private :: library = DENSELA_LAPACK

contains

!****************************************************
subroutine densela_getrf(library, m, n, A, lda, ipiv)
!****************************************************
! LU factorization of a matrix.
      use module_utils
      implicit none
! Numerical library to use
      integer,intent(in) :: library
! Number of rows in the matrix
      integer,intent(in) :: m
! Number of columns in the matrix
      integer,intent(in) :: n
! Matrix A on input, LU factors on output
      real(kr),intent(inout) :: A(lda,*)
! Leading dimension of A
      integer,intent(in) :: lda
! Sequence of pivots
      integer,intent(out),allocatable :: ipiv(:)

! local vars
      character(*),parameter:: routine_name = 'densela_getrf'
      integer :: lapack_info

      if (.not. allocated(ipiv)) then
         allocate(ipiv(min(m,n)))
      end if
      if (size(ipiv) < min(m,n)) then
         call error(routine_name, 'Size of ipiv not sufficient.')
      end if

      select case (library)
         case (DENSELA_LAPACK)
            ! LAPACK version
            if      (kr == REAL64) then
               ! double precision
               call dgetrf(m, n, A, lda, ipiv, lapack_info)
            else if (kr == REAL32) then
               ! single precision
               call sgetrf(m, n, A, lda, ipiv, lapack_info)
            end if
#if defined(BDDCML_WITH_MAGMA)
         case (DENSELA_MAGMA)
            !if      (kr == REAL64) then
               ! double precision
            call magmaf_dgetrf(m, n, A, lda, ipiv, lapack_info)
            !else if (kr == REAL32) then
               ! single precision
            !   call magmaf_sgetrf(m, n, A, lda, ipiv, lapack_info)
            !end if
#endif
         case default
            call error(routine_name, "Illegal library.")
      end select

end subroutine

!**********************************************************************
subroutine densela_getrs(library, trans, n, nrhs, A, lda, ipiv, B, ldb)
!**********************************************************************
! Solve by LU factorization
      use module_utils
      implicit none
! Numerical library to use
      integer,intent(in) :: library
! Solve transposed problem
      character, intent(in) :: trans
! Number of rows and columns in the matrix
      integer,intent(in) :: n
! Number of right-hand sides
      integer,intent(in) :: nrhs
! Matrix A with LU factors
      real(kr),intent(in) :: A(lda,*)
! Leading dimension of A
      integer,intent(in) :: lda
! Sequence of pivots
      integer,intent(in) :: ipiv(:)
! On entry, the right hand side matrix B.
! On exit, the solution matrix X.
      real(kr),intent(inout) :: B(ldb,*)
! Leading dimension of B
      integer,intent(in) :: ldb
	 		
! local vars
      character(*),parameter:: routine_name = 'densela_getrs'
      integer :: lapack_info

      select case (library)
         case (DENSELA_LAPACK)
            ! LAPACK version
            if      (kr == REAL64) then
               ! double precision
               call dgetrs(trans, n, nrhs, A, lda, ipiv, B, ldb, lapack_info)
            else if (kr == REAL32) then
               ! single precision
               call sgetrs(trans, n, nrhs, A, lda, ipiv, B, ldb, lapack_info)
            end if
#if defined(BDDCML_WITH_MAGMA)
         case (DENSELA_MAGMA)
            if      (kr == REAL64) then
               ! double precision
               call dgetrs(trans, n, nrhs, A, lda, ipiv, B, ldb, lapack_info)
            else if (kr == REAL32) then
               ! single precision
               call sgetrs(trans, n, nrhs, A, lda, ipiv, B, ldb, lapack_info)
            end if
#endif
         case default
            call error(routine_name, "Illegal library.")
      end select

end subroutine

!*******************************************************
subroutine densela_sytrf(library, uplo, n, A, lda, ipiv)
!*******************************************************
! LDLT factorization of a matrix.
      use module_utils
      implicit none
! Numerical library to use
      integer,intent(in) :: library
! Upper or lower part of the matrix
      character, intent(in) :: uplo
! Number of rows and columns in the matrix
      integer,intent(in) :: n
! Matrix A on input, LDLT factors on output
      real(kr),intent(inout) :: A(lda,*)
! Leading dimension of A
      integer,intent(in) :: lda
! Sequence of pivots
      integer,intent(out),allocatable :: ipiv(:)

      integer, external :: ilaenv

! local vars
      character(*),parameter:: routine_name = 'densela_sytrf'
      real(kr), allocatable :: work(:)
      integer :: nb
      integer :: lwork
      integer :: lapack_info
      character(6) :: function_name

      if (.not. allocated(ipiv)) then
         allocate(ipiv(n))
      end if
      if (size(ipiv) < n) then
         call error(routine_name, 'Size of ipiv not sufficient.')
      end if

      select case (library)
         case (DENSELA_LAPACK)
            ! LAPACK version
            if      (kr == REAL64) then
               ! double precision
               function_name = 'DSYTRF'
            else if (kr == REAL32) then
               ! single precision
               function_name = 'SSYTRF'
            end if

            nb = ilaenv(1, function_name, uplo,  n, 0, 0, 0)
            lwork = lda*nb
            allocate(work(lwork))

            if      (kr == REAL64) then
               ! double precision
               call dsytrf(uplo, n, A, lda, ipiv, work, lwork, lapack_info)
            else if (kr == REAL32) then
               ! single precision
               call ssytrf(uplo, n, A, lda, ipiv, work, lwork, lapack_info)
            end if

            deallocate(work)
#if defined(BDDCML_WITH_MAGMA)
         case (DENSELA_MAGMA)
            !if      (kr == REAL64) then
               ! double precision
            call magmaf_dsytrf(uplo, n, A, lda, ipiv, lapack_info)
            !else if (kr == REAL32) then
               ! single precision
            !   call magmaf_ssytrf(uplo, n, A, lda, ipiv, lapack_info)
            !end if
#endif
         case default
            call error(routine_name, "Illegal library.")
      end select

end subroutine

!*********************************************************************
subroutine densela_sytrs(library, uplo, n, nrhs, A, lda, ipiv, B, ldb)
!*********************************************************************
! LDLT factorization of a matrix.
      use module_utils
      implicit none
! Numerical library to use
      integer,intent(in) :: library
! Solve transposed problem
      character, intent(in) :: uplo
! Number of rows and columns in the matrix
      integer,intent(in) :: n
! Number of right-hand sides
      integer,intent(in) :: nrhs
! Matrix A with LDLT factors
      real(kr),intent(in) :: A(lda,*)
! Leading dimension of A
      integer,intent(in) :: lda
! Sequence of pivots
      integer,intent(in) :: ipiv(:)
! On entry, the right hand side matrix B.
! On exit, the solution matrix X.
      real(kr),intent(inout) :: B(ldb,*)
! Leading dimension of B
      integer,intent(in) :: ldb

! local vars
      character(*),parameter:: routine_name = 'densela_sytrs'
      integer :: lapack_info

      select case (library)
         case (DENSELA_LAPACK)
            ! LAPACK version
            if      (kr == REAL64) then
               ! double precision
               call dsytrs(uplo, n, nrhs, A, lda, ipiv, B, ldb, lapack_info)
            else if (kr == REAL32) then
               ! single precision
               call ssytrs(uplo, n, nrhs, A, lda, ipiv, B, ldb, lapack_info)
            end if
#if defined(BDDCML_WITH_MAGMA)
         case (DENSELA_MAGMA)
            if      (kr == REAL64) then
               ! double precision
               call dsytrs(uplo, n, nrhs, A, lda, ipiv, B, ldb, lapack_info)
            else if (kr == REAL32) then
               ! single precision
               call ssytrs(uplo, n, nrhs, A, lda, ipiv, B, ldb, lapack_info)
            end if
#endif
         case default
            call error(routine_name, "Illegal library.")
      end select

end subroutine

!***********************************************************************************
subroutine densela_gemv(library, trans, m, n, alpha, A, lda, x, incx, beta, y, incy)
!***********************************************************************************
! Matrix-vector product in the form y := alpha*A*x + beta*y, for general A
      use module_utils
      use iso_c_binding
      implicit none
! Numerical library to use
      integer,intent(in) :: library
! Solve transposed problem
      character, intent(in) :: trans
! Number of rows in the matrix
      integer,intent(in) :: m
! Number of columns in the matrix
      integer,intent(in) :: n
! Coefficient alpha 
      real(kr),intent(in) :: alpha
! Matrix A
      real(kr),intent(in) :: A(lda,*)
! Leading dimension of A
      integer,intent(in) :: lda
! Input array X
      real(kr),intent(in) :: x(*)
! Increment of X
      integer,intent(in) :: incx
! Coefficient beta 
      real(kr),intent(in) :: beta
! Input/output array Y
      real(kr),intent(inout) :: y(*)
! Increment of Y
      integer,intent(in) :: incy

! local vars
      character(*),parameter:: routine_name = 'densela_gemv'
#if defined(BDDCML_WITH_MAGMA)
      integer(8) :: queue
      integer(c_int) :: device
      integer(8) :: dA
      integer(8) :: dx
      integer(8) :: dy
      integer :: ldda
      integer :: ierr
      integer :: lx, ly
#endif

      select case (library)
         case (DENSELA_LAPACK)
            ! LAPACK version
            if      (kr == REAL64) then
               ! double precision
               call dgemv(trans, m, n, alpha, A, lda, x, incx, beta, y, incy)
            else if (kr == REAL32) then
               ! single precision
               call sgemv(trans, m, n, alpha, A, lda, x, incx, beta, y, incy)
            end if
#if defined(BDDCML_WITH_MAGMA)
         case (DENSELA_MAGMA)
            ! simly call LAPACK anyway
            !call dgemv(trans, m, n, alpha, A, lda, x, incx, beta, y, incy)

            call magmaf_getdevice(device)
            call magmaf_queue_create(device, queue)

            ! allocate memory on GPU
            ierr = magmaf_dmalloc(dA, m*n)

            if (trans == 'T' .or. trans == 't') then
               lx = m
               ly = n
            else
               lx = n
               ly = m
            end if
            ierr = magmaf_dmalloc(dx, lx)
            ierr = magmaf_dmalloc(dy, ly)

            ! copy the matrix from CPU to GPU
            lddA = m
            call magmaf_dsetmatrix( m, n, A, lda, dA, lddA, queue )

            ! copy vectors from CPU to GPU
            call magmaf_dsetvector(lx, x, incx, dx, 1, queue)
            call magmaf_dsetvector(ly, y, incy, dy, 1, queue)

            !if      (kr == REAL64) then
            !   ! double precision
            call magmaf_dgemv(trans, m, n, alpha, dA, lddA, dx, 1, beta, dy, 1, queue)
            !else if (kr == REAL32) then
            !   ! single precision
            !   call magmaf_sgemv(trans, m, n, alpha, A, lda, x, incx, beta, y, incy, queue)
            !end if

            ! copy data from GPU to CPU
            call magmaf_dgetvector(ly, dy, 1, y, incy, queue)

            call magmaf_queue_destroy(queue)

            ! free memory on GPU
            ierr = magmaf_free(dA)
            ierr = magmaf_free(dx)
            ierr = magmaf_free(dy)
#endif
         case default
            call error(routine_name, "Illegal library.")
      end select

end subroutine

!*******************************************************************************
subroutine densela_symv(library, uplo, n, alpha, A, lda, x, incx, beta, y, incy)
!*******************************************************************************
! Matrix-vector product in the form y := alpha*A*x + beta*y, for symmetric A
      use module_utils
      use iso_c_binding
      implicit none
! Numerical library to use
      integer,intent(in) :: library
! Upper or lower part of the matrix
      character, intent(in) :: uplo
! Size of the matrix
      integer,intent(in) :: n
! Coefficient alpha 
      real(kr),intent(in) :: alpha
! Matrix A
      real(kr),intent(in) :: A(lda,*)
! Leading dimension of A
      integer,intent(in) :: lda
! Input array X
      real(kr),intent(in) :: x(*)
! Increment of X
      integer,intent(in) :: incx
! Coefficient beta 
      real(kr),intent(in) :: beta
! Input/output array Y
      real(kr),intent(inout) :: y(*)
! Increment of Y
      integer,intent(in) :: incy

! local vars
      character(*),parameter:: routine_name = 'densela_symv'
#if defined(BDDCML_WITH_MAGMA)
      integer(8) :: queue
      integer(c_int) :: device
      integer(8) :: dA
      integer(8) :: dx
      integer(8) :: dy
      integer :: ldda
      integer :: ierr
#endif

      select case (library)
         case (DENSELA_LAPACK)
            ! LAPACK version
            if      (kr == REAL64) then
               ! double precision
               call dsymv(uplo, n, alpha, A, lda, x, incx, beta, y, incy)
            else if (kr == REAL32) then
               ! single precision
               call ssymv(uplo, n, alpha, A, lda, x, incx, beta, y, incy)
            end if
#if defined(BDDCML_WITH_MAGMA)
         case (DENSELA_MAGMA)
            ! simply call LAPACK anyway
            !call dsymv(uplo, n, alpha, A, lda, x, incx, beta, y, incy)

            call magmaf_getdevice(device)
            call magmaf_queue_create(device, queue)

            ! allocate memory on GPU
            ierr = magmaf_dmalloc(dA, n*n)
            ierr = magmaf_dmalloc(dx, n)
            ierr = magmaf_dmalloc(dy, n)

            ! copy the matrix from CPU to GPU
            lddA = n
            call magmaf_dsetmatrix( n, n, A, lda, dA, lddA, queue )

            ! copy vectors from CPU to GPU
            call magmaf_dsetvector(n, x, incx, dx, 1, queue)
            call magmaf_dsetvector(n, y, incx, dy, 1, queue)

            !if      (kr == REAL64) then
            !   ! double precision
            call magmaf_dsymv(uplo, n, alpha, dA, lda, dx, incx, beta, dy, incy, queue)
            !else if (kr == REAL32) then
            !   ! single precision
            !   call magmaf_ssymv(uplo, n, alpha, A, lda, x, incx, beta, y, incy, queue)
            !end if
 
            ! copy data from GPU to CPU
            call magmaf_dgetvector(n, dy, 1, y, incy, queue)

            call magmaf_queue_destroy(queue)

            ! free memory on GPU
            ierr = magmaf_free(dA)
            ierr = magmaf_free(dx)
            ierr = magmaf_free(dy)
#endif
         case default
            call error(routine_name, "Illegal library.")
      end select

end subroutine

end module module_densela