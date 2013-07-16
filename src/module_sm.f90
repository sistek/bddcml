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

module module_sm
!***************
! Module for operations with sparse matrices in triplets I,J,A
! Jakub Sistek, Denver, 12/2007

! type of real variables
      integer,parameter,private :: kr = kind(1.D0)
! numerical zero
      real(kr),parameter,private :: numerical_zero = 1.e-16_kr

contains

!*********************************************************************************
subroutine sm_pmd_get_length(matrixtype,nelem,inet,linet,nnet,lnnet,nndf,lnndf,la)
!*********************************************************************************
! Subroutine for getting size of sparse matrix to be allocated from PMD
! description of mesh.
! Assumes storing all elements in memory without assemblage.

      implicit none

! Type of matrix - determine its storage
! 0 - unsymmetric                 -> full element matrices
! 1 - symmetric positive definite -> only upper triangle of element matrix
! 2 - symmetric general           -> only upper triangle of element matrix
      integer,intent(in) :: matrixtype

! Problem size - number of elements
      integer,intent(in) :: nelem

! description of PMD mesh
      integer,intent(in) :: linet,       lnnet,       lnndf
      integer,intent(in) ::  inet(linet), nnet(lnnet), nndf(lnndf)

! Length of matrix in IJA sparse format
      integer,intent(out)  :: la

! Local variables
      integer :: ie, nevab, ndofn, iinet, nne, inod, inodg, lelm
                
! Find the lenght of sparse matrix to be allocated LA 
      iinet = 0
      la = 0
      do ie = 1,nelem
         nevab = 0
         nne = nnet(ie)
         do inod = 1,nne
            iinet = iinet + 1
            inodg = inet(iinet)
            ndofn = nndf(inodg)
            nevab = nevab + ndofn
         end do
         if (matrixtype.eq.1.or.matrixtype.eq.2) then
            lelm = (nevab*(nevab+1))/2
         else if (matrixtype.eq.0) then
            lelm = nevab*nevab
         else
            write(*,*) 'sm_pmd_get_length: Unknown type of matrix.'
            stop
         end if
         la = la + lelm
      end do

      return
end subroutine

!***************************************************************************************************
subroutine sm_pmd_get_length_quart(matrixtype,nelem,inet,linet,nnet,lnnet,nndf,lnndf,iingn,liingn, &
                                   la11,la12,la21,la22)
!***************************************************************************************************
! Subroutine for getting size of 4 parts of a sparse matrix to be allocated from PMD
! description of mesh.
! Assumes storing all elements in memory without assemblage.
! The blocking is denoted, all blocks are sparse
! | A11  A12 |
! | A21  A22 |
! Variables to nodes in 22 block are determined by array iinode(ninode)

      implicit none

! Type of matrix - determine its storage
! 0 - unsymmetric                 -> full element matrices
! 1 - symmetric positive definite -> only upper triangle of element matrix
! 2 - symmetric general           -> only upper triangle of element matrix
      integer,intent(in) :: matrixtype

! Problem size - number of elements
      integer,intent(in) :: nelem

! description of PMD mesh
      integer,intent(in) :: linet,       lnnet,       lnndf
      integer,intent(in) ::  inet(linet), nnet(lnnet), nndf(lnndf)

! indices of interface nodes
      integer,intent(in) :: liingn
      integer,intent(in) ::  iingn(liingn)

! Length of matrix in IJA sparse format
      integer,intent(out)  :: la11, la12, la21, la22

! Local variables
      integer :: ie, nevab1, nevab2, ndofn, iinet, nne, inod, inodg, lelm11, lelm12, lelm21, lelm22 
                
! Find the lenght of sparse matrix to be allocated LA 
      iinet = 0
      la11  = 0
      la12  = 0
      la21  = 0
      la22  = 0
      do ie = 1,nelem
         nevab1 = 0
         nevab2 = 0
         nne = nnet(ie)
         do inod = 1,nne
            iinet = iinet + 1
            inodg = inet(iinet)
            ndofn = nndf(inodg)
            if (any(iingn.eq.inodg)) then
               ! node is in block 2
               nevab2 = nevab2 + ndofn
            else
               ! node is in block 1
               nevab1 = nevab1 + ndofn
            end if
         end do
         if (matrixtype.eq.1.or.matrixtype.eq.2) then
            ! symmetric case
            lelm11 = (nevab1*(nevab1+1))/2
            lelm22 = (nevab2*(nevab2+1))/2
            lelm12 = nevab1*nevab2
            lelm21 = 0
         else if (matrixtype.eq.0) then
            ! unsymmetric case
            lelm11 = nevab1*nevab1
            lelm22 = nevab2*nevab2
            lelm12 = nevab1*nevab2
            lelm21 = nevab2*nevab1
         else
            write(*,*) 'sm_pmd_get_length: Unknown type of matrix.'
            stop
         end if
         la11 = la11 + lelm11
         la12 = la12 + lelm12
         la21 = la21 + lelm21
         la22 = la22 + lelm22
      end do

      return
end subroutine

!*******************************************************************************
subroutine sm_pmd_load(matrixtype, idelm,nelem,inet,linet,nnet,lnnet,nndf,lnndf,&
                       kdof,lkdof, i_sparse, j_sparse, a_sparse, la)
!*******************************************************************************
! Subroutine for reading element matrices in PMD format and storing them in IJA
! sparse format with repeated entries (without assembly), stores element by element
! element matrices are read column-by-column, for symmetric storage only up to
! the diagonal

      implicit none

! Type of matrix - determine its storage
! 0 - unsymmetric                 -> full element matrices
! 1 - symmetric positive definite -> only upper triangle of element matrix
! 2 - symmetric general           -> only upper triangle of element matrix
      integer,intent(in) :: matrixtype

! unit associated with file with element matrices, number of elements
      integer,intent(in) :: idelm, nelem

! description of mesh
      integer,intent(in) :: linet, lnnet, lnndf, lkdof
      integer,intent(in) :: inet(linet), nnet(lnnet), nndf(lnndf), kdof(lkdof)

! Matrix in IJA sparse format
      integer,intent(in)  :: la
      integer,intent(out)  :: i_sparse(la), j_sparse(la)
      real(kr),intent(out) :: a_sparse(la)

! Local variables
      integer:: ie, i, inda, lelm

      ! Read element matrices
      inda  = 0
      do ie = 1,nelem
         read(idelm) lelm,(a_sparse(inda + i), i = 1,lelm)
         inda = inda + lelm
      end do

      ! Prepare numbering of element matrices
      call sm_pmd_make_element_numbering(matrixtype,nelem,inet,linet,nnet,lnnet,nndf,lnndf,kdof,lkdof,&
                                         i_sparse, j_sparse, la)

      return
end subroutine

!**************************************************************************************
subroutine sm_pmd_load_masked(matrixtype, idelm,nelem,inet,linet,nnet,lnnet,nndf,lnndf,&
                              kdof,lkdof, isegn,lisegn, i_sparse, j_sparse, a_sparse, la)
!**************************************************************************************
! Subroutine for reading element matrices in PMD format and storing them in IJA
! sparse format with repeated entries (without assembly), stores element by element
! selects elements from global file according to ISEGN array

      implicit none

! Type of matrix - determine its storage
! 0 - unsymmetric                 -> full element matrices
! 1 - symmetric positive definite -> only upper triangle of element matrix
! 2 - symmetric general           -> only upper triangle of element matrix
      integer,intent(in) :: matrixtype

! unit associated with file with element matrices, number of elements
      integer,intent(in) :: idelm, nelem

! description of mesh
      integer,intent(in) :: linet, lnnet, lnndf, lkdof
      integer,intent(in) :: inet(linet), nnet(lnnet), nndf(lnndf), kdof(lkdof)
      integer,intent(in) :: lisegn
      integer,intent(in) ::  isegn(lisegn)

! Matrix in IJA sparse format
      integer,intent(in)  :: la
      integer,intent(out)  :: i_sparse(la), j_sparse(la)
      real(kr),intent(out) :: a_sparse(la)

! Local variables
      integer:: ie, i, inda, lelm, indposition, inde

      ! Read element matrices from a global file
      inda  = 0
      indposition = 1
      do ie = 1,nelem
         inde = isegn(ie)
         ! move in file to the proper element position
         do i = 1,inde-indposition
            read(idelm) 
         end do
         indposition = inde
         read(idelm) lelm,(a_sparse(inda + i), i = 1,lelm)
         indposition = indposition + 1

         inda = inda + lelm
      end do

      ! Prepare numbering of element matrices
      call sm_pmd_make_element_numbering(matrixtype,nelem,inet,linet,nnet,lnnet,nndf,lnndf,kdof,lkdof,&
                                         i_sparse, j_sparse, la)

      return
end subroutine

!*******************************************************************************************
subroutine sm_pmd_make_element_numbering(matrixtype,nelem,inet,linet,nnet,lnnet,nndf,lnndf,&
                                         kdof,lkdof, i_sparse, j_sparse, la)
!*******************************************************************************************
! Subroutine for reading element matrices in PMD format and storing them in IJA
! sparse format with repeated entries (without assembly), stores element by element

      use module_utils, only: error
      implicit none

! Type of matrix - determine its storage
! 0 - unsymmetric                 -> full element matrices
! 1 - symmetric positive definite -> only upper triangle of element matrix
! 2 - symmetric general           -> only upper triangle of element matrix
      integer,intent(in) :: matrixtype

! number of elements
      integer,intent(in) :: nelem

! description of mesh
      integer,intent(in) :: linet, lnnet, lnndf, lkdof
      integer,intent(in) :: inet(linet), nnet(lnnet), nndf(lnndf), kdof(lkdof)

! Matrix in IJA sparse format - only rows and columns
      integer,intent(in)  :: la
      integer,intent(out)  :: i_sparse(la), j_sparse(la)

! Local variables
      character(*),parameter:: routine_name = 'SM_PMD_MAKE_ELEMENT_NUMBERING'
      integer,allocatable :: kdofe(:)
      integer:: inddof, idofn, ndofn, ive, jve, inda, nevab, nevax, ie, iinet, &
                inod, inodg, nve, nne
      integer:: uppbound

! Find the maximum dimension of element matrix NEVAX
      iinet = 0
      nevax = 0
      do ie = 1,nelem
         nevab = 0
         nne = nnet(ie)
         do inod = 1,nne
            iinet = iinet + 1
            inodg = inet(iinet)
            ndofn = nndf(inodg)
            nevab = nevab + ndofn
         end do
         nevax = max(nevax,nevab)
      end do

! Allocate the index array on element KDOFE
      allocate(kdofe(nevax))

! Read element matrices
      iinet = 0
      inda  = 0
      do ie = 1,nelem
         nne = nnet(ie)
         ! Construct array of pointers kdofe with global dofs on element
         ive = 0
         kdofe = 0
         do inod = 1,nne
            iinet = iinet + 1
            inodg = inet(iinet)
            ndofn = nndf(inodg)
            inddof = kdof(inodg)
            do idofn = 1,ndofn
               ive = ive + 1
               inddof = inddof + 1
               kdofe(ive) = inddof
            end do
         end do
         nve = ive

!         write(*,*) ie, kdofe

         ! Distribute the indices of entries
         do jve = 1,nve
            ! set proper end for the row loop
            if (matrixtype.eq.1.or.matrixtype.eq.2) then
               ! symmetric case
               uppbound = jve
            else if (matrixtype.eq.0) then
               ! unsymmetric case
               uppbound = nve
            else
               call error(routine_name, 'Unknown type of matrix: ',matrixtype)
            end if
            do ive = 1,uppbound 

               ! Move in the main matrix array
               inda = inda + 1

               ! Fill in the indices of entries
               if (matrixtype .eq.0 .or. ( kdofe(ive).le.kdofe(jve) ) ) then
                  ! in nonsymmetric storage or case of upper triangular entry, standard position
                  i_sparse(inda) = kdofe(ive)
                  j_sparse(inda) = kdofe(jve)
               else
                  ! in symmetric storage and case of entry bellow diagonal, mirror index along diagonal
                  i_sparse(inda) = kdofe(jve)
                  j_sparse(inda) = kdofe(ive)
               end if
            end do
         end do
      end do

! Free memory when element matrices are read
      deallocate(kdofe)

      return
end subroutine

!************************************************************************************************
subroutine sm_count_fixed_rows_matrix_size(matrixtype,ifix,lifix,nnz,i_sparse,j_sparse,la,&
                                           la_fixed)
!************************************************************************************************
! Subroutine for calculating the size of the sparse matrix for storing fixed rows of the original matrix 
! using PMD arrays
      implicit none

! Type of matrix - determine its storage
! 0 - unsymmetric                 -> full element matrices
! 1 - symmetric positive definite -> only upper triangle of element matrix
! 2 - symmetric general           -> only upper triangle of element matrix
      integer,intent(in) :: matrixtype

! indices of fixed variables
      integer,intent(in) :: lifix
      integer,intent(in) :: ifix(lifix)

! Matrix in IJA sparse format
      integer,intent(in) :: nnz
      integer,intent(in) :: la
      integer,intent(inout) :: i_sparse(la), j_sparse(la)

      integer,intent(out) ::    la_fixed


! Local variables
      character(*),parameter:: routine_name = 'SM_COUNT_FIXED_ROWS_MATRIX_SIZE'
      integer :: ia, indi, indj
      integer :: ifixed 

! Eliminate boundary conditions
      ifixed = 0
      do ia = 1,nnz
         indi = i_sparse(ia)
         indj = j_sparse(ia)

         if (ifix(indi).ne.0 .or. ( ( matrixtype.eq.1 .or. matrixtype.eq.2 ) .and. ifix(indj).ne.0 ) ) then
            ifixed = ifixed + 1
         end if

      end do

      la_fixed = ifixed

end subroutine

!************************************************************************************************
subroutine sm_apply_bc(matrixtype,ifix,lifix,fixv,lfixv,nnz,i_sparse,j_sparse,a_sparse,la,bc,lbc,&
                       store_fixed_rows, la_fixed, i_fixed_sparse, j_fixed_sparse, a_fixed_sparse)
!************************************************************************************************
! Subroutine for application of Dirichlet boudary conditions on a sparse matrix
! using PMD arrays
! Eliminate boundary conditions from right hand side, put zeros to the row and
! column and one on the diagonal. 
! The value of fixed variable put in RHS.

      use module_utils
      implicit none

! Type of matrix - determine its storage
! 0 - unsymmetric                 -> full element matrices
! 1 - symmetric positive definite -> only upper triangle of element matrix
! 2 - symmetric general           -> only upper triangle of element matrix
      integer,intent(in) :: matrixtype

! indices of fixed variables
      integer,intent(in) :: lifix
      integer,intent(in) :: ifix(lifix)

! values of fixed variables
      integer,intent(in) :: lfixv
      real(kr),intent(in):: fixv(lfixv)

! Matrix in IJA sparse format
      integer,intent(in) :: nnz
      integer,intent(in) :: la
      integer,intent(inout) :: i_sparse(la), j_sparse(la)
      real(kr),intent(inout):: a_sparse(la)

! Should the fixed rows be stored during elimination?
      logical,intent(in) ::            store_fixed_rows

! Matrix of fixed rows in IJA sparse format
      integer,intent(in) ::             la_fixed
      integer,intent(inout),optional :: i_fixed_sparse(la_fixed), j_fixed_sparse(la_fixed)
      real(kr),intent(inout),optional:: a_fixed_sparse(la_fixed)

! values of fixed variables
      integer,intent(in)   :: lbc
      real(kr),intent(out) :: bc(lbc)

! Local variables
      character(*),parameter:: routine_name = 'SM_APPLY_BC'
      integer :: ia, indi, indj
      real(kr):: fixval, aval, value, auxval
      integer :: i, newplace
      integer ::            lis_diag_fixed
      logical,allocatable :: is_diag_fixed(:)
      integer :: ifixed

! Zero BC
      if (lbc.gt.0) then
         bc = 0.0_kr
      end if

! Auxiliary value to put on diagonal after eliminating BC
      auxval = 1.0_kr

! Auxiliary array checking that no diagonal entry was omitted - even for saddle point systems
      lis_diag_fixed = lifix
      allocate( is_diag_fixed(lis_diag_fixed) )
      is_diag_fixed = .false.

! Eliminate boundary conditions
      ifixed = 0
      do ia = 1,nnz
         indi = i_sparse(ia)
         indj = j_sparse(ia)

         ! store the row if desired
         if (store_fixed_rows) then
            if (ifix(indi).ne.0 .or. ( ( matrixtype.eq.1 .or. matrixtype.eq.2 ) .and. ifix(indj).ne.0 ) ) then
               ifixed = ifixed + 1
               if (ifixed.gt.la_fixed) then
                  call error ( routine_name, ' No space left in matrix of fixed rows.' )
               end if
               i_fixed_sparse(ifixed) = indi
               j_fixed_sparse(ifixed) = indj
               a_fixed_sparse(ifixed) = a_sparse(ia)
            end if
         end if

         ! off-diagonal entries
         if (indi.ne.indj) then
            ! do this only in the case of symmetric storage

            if (matrixtype.ne.0) then
               ! fixed row
               if (ifix(indi).ne.0) then
                  fixval = fixv(indi)
                  aval   = a_sparse(ia)
                  if (fixval.ne.0.0_kr.and.aval.ne.0.0_kr) then
                     value = aval*fixval
                     if (ifix(indj).eq.0) then
                        bc(indj) = bc(indj) - value
                     end if
                  end if
               end if
            end if

            ! fixed column
            if (ifix(indj).ne.0) then
               fixval = fixv(indj)
               aval   = a_sparse(ia)
               if (fixval.ne.0.0_kr.and.aval.ne.0.0_kr) then
                  value = aval*fixval
                  if (ifix(indi).eq.0) then
                     bc(indi) = bc(indi) - value
                  end if
               end if
            end if

            ! zero entry
            if (ifix(indi).ne.0.or.ifix(indj).ne.0) then
               a_sparse(ia) = 0.0_kr
            end if
         ! diagonal entries
         else 
            ! fixed row
            if (ifix(indi).ne.0) then
               fixval = fixv(indi)
               aval   = a_sparse(ia)
               if (abs(aval).lt.numerical_zero) then
                  a_sparse(ia) = auxval
                  aval         = auxval
               end if
               if (fixval.ne.0.0_kr.and.aval.ne.0.0_kr) then
                  value = aval*fixval
                  bc(indi) = bc(indi) + value
               end if
               is_diag_fixed(indi) = .true.
            end if
         end if

      end do

      ! make sure that all diagonal entries were fixed - problem with saddle point systems
      ! eventually write them behind existing entries
      newplace = nnz + 1
      do i = 1,lifix
         if ( ifix(i) .gt. 0 .and. .not. is_diag_fixed(i) ) then
            ! value was not fixed, add it

            ! check space
            if ( newplace .gt. la ) then
                call error ( routine_name, ' No space left in matrix for adding new entries induced by BC.' )
            end if

            fixval = fixv(i)

            i_sparse(newplace) = i
            j_sparse(newplace) = i

            a_sparse(newplace) = auxval
            value = auxval*fixval
            bc(i) = bc(i) + value

            is_diag_fixed(i) = .true.

            newplace = newplace + 1
         end if
      end do

      deallocate( is_diag_fixed )

end subroutine

!****************************************************
subroutine sm_prepare_rhs(ifix,lifix,bc,lbc,rhs,lrhs)
!****************************************************
! Subroutine for preparation of RHS so that it satisfy Dirichlet BC 
! For nonhomogeneous Dirichlet BC using array BC from subroutine SM_APPLY_BC.

      implicit none

      integer,intent(in)     :: lifix, lbc, lrhs
      integer,intent(in)     :: ifix(lifix)
      real(kr),intent(in)    :: bc(lbc)
      real(kr),intent(inout) :: rhs(lrhs)

! Natural RHS is zero in constraints
      if (lifix.gt.0) then
         where(ifix.ne.0) rhs = 0.0_kr
      end if

! Add eliminated variables from BC (if it has sense)
      if (lbc.eq.lrhs) then
         rhs = rhs + bc
      end if

! Right hand side is now consistent with matrix

      return
end subroutine

!************************************************************
subroutine sm_assembly(i_sparse, j_sparse, a_sparse, la, nnz)
!************************************************************
! Subroutine for assemblage of sparse matrix in IJA format 
! Sorts the matrix and adds entries with same indices

      implicit none

! Matrix in IJA sparse format
      integer,intent(in) :: la
      integer,intent(inout) :: i_sparse(la), j_sparse(la)
      real(kr),intent(inout):: a_sparse(la)
! Number of nonzeros
      integer,intent(out) :: nnz
      
! Local variables
      integer:: inz, indi, indj, indi_old, indj_old, ia

! In case of zero length skip to the end
      if (la.lt.1) then
         nnz = 0
         goto 77
      end if

! Sort entries by indices I and J
      call sm_sort(i_sparse,j_sparse,a_sparse,la)

! Assembly entries
      inz = 1
      indi_old = i_sparse(1)
      indj_old = j_sparse(1)
      do ia = 2,la
         indi = i_sparse(ia)
         indj = j_sparse(ia)
         !if (abs(a_sparse(ia)).gt.numerical_zero) then
         if (indi.ne.indi_old.or.indj.ne.indj_old) then
            inz = inz + 1
            a_sparse(inz) = a_sparse(ia)
            i_sparse(inz) = i_sparse(ia)
            j_sparse(inz) = j_sparse(ia)
            indi_old = indi
            indj_old = indj
         else
            a_sparse(inz) = a_sparse(inz) + a_sparse(ia)
         end if
         !end if
      end do
      nnz = inz

! Zero the tail
      a_sparse(nnz+1:la) = 0._kr

  77  return
end subroutine

!**********************************************************
RECURSIVE SUBROUTINE sm_sort(i_sparse,j_sparse,a_sparse,la)
!**********************************************************
! Subroutine for sorting sparse matrix with repeated entries.
! Sorts entries by I (row) index and then by J (column) index.
! Based on Quicksort algorithm by
! Brainerd, W.S., Goldberg, C.H. & Adams, J.C. (1990) "Programmer's Guide to
! Fortran 90", McGraw-Hill  ISBN 0-07-000248-7, pages 149-150.
! Modified by Alan Miller.
! Edited by Jakub Sistek for purpose of sorting entries of sparse matrices IJA.

implicit none

! Matrix in IJA sparse format
integer, intent(in) :: la
integer,intent(inout),target :: i_sparse(la), j_sparse(la)
real(kr),intent(inout),target:: a_sparse(la)

integer, pointer  :: p_i_sparse(:), p_j_sparse(:)
real(kr), pointer :: p_a_sparse(:)

! Local variables
integer :: itemp, jtemp, indend, iend, indbegin, ibegin
real(kr):: atemp

! Associate pointers with arrays
p_i_sparse => i_sparse
p_j_sparse => j_sparse
p_a_sparse => a_sparse

! Sort elements by row index
call quick_sort(1,la)

! Interchange rows and columns for next step
p_i_sparse => j_sparse
p_j_sparse => i_sparse

! Sort elements by column index
ibegin = 1
iend   = ibegin
do while (iend.le.la)
   indbegin = p_j_sparse(ibegin)
   indend   = p_j_sparse(iend)
   ! find the index of last same value
   do while (indend.eq.indbegin)
      iend = iend + 1
      ! correction at the end of the array
      if (iend.gt.la) exit
      indend = p_j_sparse(iend)
   end do
   iend = iend - 1

!   write(*,*) 'ibegin',ibegin,'iend',iend
   call quick_sort(ibegin,iend)

   ibegin = iend + 1
   iend   = ibegin
end do

! Nullify pointers
nullify(p_i_sparse)
nullify(p_j_sparse)
nullify(p_a_sparse)

contains

RECURSIVE SUBROUTINE quick_sort(left_end, right_end)

INTEGER, INTENT(IN) :: left_end, right_end

!     Local variables
INTEGER             :: i, j
REAL                :: reference
INTEGER, PARAMETER  :: max_simple_sort_size = 6

IF (right_end < left_end + max_simple_sort_size) THEN
  ! Use interchange sort for small lists
  CALL interchange_sort(left_end, right_end)

ELSE
  ! Use partition ("quick") sort
  reference = p_i_sparse((left_end + right_end)/2)
  i = left_end - 1; j = right_end + 1

  DO
    ! Scan p_i_sparse from left end until element >= reference is found
    DO
      i = i + 1
      IF (p_i_sparse(i) >= reference) EXIT
    END DO
    ! Scan p_i_sparse from right end until element <= reference is found
    DO
      j = j - 1
      IF (p_i_sparse(j) <= reference) EXIT
    END DO


    IF (i < j) THEN
      ! Swap two out-of-order elements
      jtemp = p_j_sparse(i); p_j_sparse(i) = p_j_sparse(j); p_j_sparse(j) = jtemp
      atemp = p_a_sparse(i); p_a_sparse(i) = p_a_sparse(j); p_a_sparse(j) = atemp
      itemp = p_i_sparse(i); p_i_sparse(i) = p_i_sparse(j); p_i_sparse(j) = itemp
    ELSE IF (i == j) THEN
      i = i + 1
      EXIT
    ELSE
      EXIT
    END IF
  END DO

  IF (left_end < j) CALL quick_sort(left_end, j)
  IF (i < right_end) CALL quick_sort(i, right_end)
END IF

END SUBROUTINE quick_sort

SUBROUTINE interchange_sort(left_end, right_end)

INTEGER, INTENT(IN) :: left_end, right_end

!     Local variables
INTEGER             :: i, j

DO i = left_end, right_end - 1
  DO j = i+1, right_end
    IF (p_i_sparse(i) > p_i_sparse(j)) THEN
      jtemp = p_j_sparse(i); p_j_sparse(i) = p_j_sparse(j); p_j_sparse(j) = jtemp
      atemp = p_a_sparse(i); p_a_sparse(i) = p_a_sparse(j); p_a_sparse(j) = atemp
      itemp = p_i_sparse(i); p_i_sparse(i) = p_i_sparse(j); p_i_sparse(j) = itemp
    END IF
  END DO
END DO

END SUBROUTINE interchange_sort

END SUBROUTINE sm_sort

!**************************************************************************
subroutine sm_vec_mult(matrixtype, nnz, i_sparse, j_sparse, a_sparse, la, &
                       vec_in,lvec_in, vec_out,lvec_out)
!**************************************************************************
! Subroutine for multiplication of sparse matrix in IJA format 
! with a vector vec_in, result in vec_out

      implicit none

! Type of the matrix
! 0 - unsymmetric
! 1 - symmetric positive definite
! 2 - general symmetric
      integer :: matrixtype

! Matrix in IJA sparse format
      integer,intent(in) :: nnz, la
      integer,intent(in) :: i_sparse(la), j_sparse(la)
      real(kr),intent(in):: a_sparse(la)

! Vectors
      integer,intent(in)  :: lvec_in, lvec_out
      real(kr),intent(in) :: vec_in(lvec_in)
      real(kr),intent(out):: vec_out(lvec_out)

! Zero the output vector
      vec_out = 0.0_kr

! Decide according to value of matrixtype which multiplication is appropriate
      select case (matrixtype)
      case (0)
         ! matrix is stored whole
         call sm_vec_mult_full(nnz, i_sparse, j_sparse, a_sparse, la, &
                               vec_in,lvec_in, vec_out,lvec_out)
      case (1,2)
         ! only one triangle of matrix is stored
         call sm_vec_mult_sym(nnz, i_sparse, j_sparse, a_sparse, la, &
                              vec_in,lvec_in, vec_out,lvec_out)
      case default
         write(*,*) 'sm_vec_mult: Unknown type of matrix. Maybe MATRIXTYPE not set.'
         stop
      end select

      return
end subroutine

!**************************************************************************
subroutine sm_vec_mult_mask(matrixtype, nnz, i_sparse, j_sparse, a_sparse, la, &
                            vec_in,lvec_in, vec_out,lvec_out, &
                            mask,lmask,match)
!**************************************************************************
! Subroutine for multiplication of sparse matrix in IJA format 
! with a vector vec_in, result in vec_out

      implicit none

! Type of the matrix
! 0 - unsymmetric
! 1 - symmetric positive definite
! 2 - general symmetric
      integer :: matrixtype

! Matrix in IJA sparse format
      integer,intent(in) :: nnz, la
      integer,intent(in) :: i_sparse(la), j_sparse(la)
      real(kr),intent(in):: a_sparse(la)

! Vectors
      integer,intent(in)  :: lvec_in, lvec_out
      real(kr),intent(in) :: vec_in(lvec_in)
      real(kr),intent(out):: vec_out(lvec_out)

! Blocks of matrix
      ! selects a block of matrix, e.g. corresponding to interface
      integer,intent(in)  :: lmask
      integer,intent(in)  ::  mask(lmask)
      ! describes, if the row and column index should be selected in mask
      ! true true   -> interface block    A_22
      ! true false  -> intermediate block A_21
      ! false true  -> intermediate block A_12
      ! false false -> interior block     A_11
      logical,intent(in)  ::  match(2)

! Zero the output vector
      vec_out = 0.0_kr

! Decide according to value of matrixtype which multiplication is appropriate
      select case (matrixtype)
      case (0)
         ! matrix is stored whole
         call sm_vec_mult_full_mask(nnz, i_sparse, j_sparse, a_sparse, la, &
                                    vec_in,lvec_in, vec_out,lvec_out, &
                                    mask,lmask,match)
      case (1,2)
         ! only one triangle of matrix is stored
         call sm_vec_mult_sym_mask(nnz, i_sparse, j_sparse, a_sparse, la, &
                                   vec_in,lvec_in, vec_out,lvec_out, &
                                   mask,lmask,match)
      case default
         write(*,*) 'sm_vec_mult_mask: Unknown type of matrix. Maybe MATRIXTYPE not set.'
         stop
      end select

      return
end subroutine

!*******************************************************************
subroutine sm_vec_mult_full(nnz, i_sparse, j_sparse, a_sparse, la, &
                            vec_in,lvec_in, vec_out,lvec_out)
!*******************************************************************
! Subroutine for multiplication of FULL sparse matrix in IJA format 
! with a vector vec_in, result in vec_out
      use module_utils, only: error_exit

      implicit none
! Matrix in IJA sparse format
      integer,intent(in) :: nnz, la
      integer,intent(in) :: i_sparse(la), j_sparse(la)
      real(kr),intent(in):: a_sparse(la)

! Vectors
      integer,intent(in)  :: lvec_in, lvec_out
      real(kr),intent(in) :: vec_in(lvec_in)
      real(kr),intent(out):: vec_out(lvec_out)

! local variables
      integer:: ia, indi, indj

! Do the multiplication for each non-zero entry
      do ia = 1,nnz
         indi  = i_sparse(ia)
         indj  = j_sparse(ia)

         vec_out(indi) = vec_out(indi) + a_sparse(ia)*vec_in(indj)
      end do

      return
end subroutine

!*******************************************************************
subroutine sm_vec_mult_full_mask(nnz, i_sparse, j_sparse, a_sparse, la, &
                                 vec_in,lvec_in, vec_out,lvec_out, &
                                 mask,lmask,match)
!*******************************************************************
! Subroutine for multiplication of FULL sparse matrix in IJA format 
! with a vector vec_in, result in vec_out
      use module_utils, only: error_exit

      implicit none
! Matrix in IJA sparse format
      integer,intent(in) :: nnz, la
      integer,intent(in) :: i_sparse(la), j_sparse(la)
      real(kr),intent(in):: a_sparse(la)

! Vectors
      integer,intent(in)  :: lvec_in, lvec_out
      real(kr),intent(in) :: vec_in(lvec_in)
      real(kr),intent(out):: vec_out(lvec_out)

! Blocks of matrix
      ! selects a block of matrix, e.g. corresponding to interface
      integer,intent(in)  :: lmask
      integer,intent(in)  ::  mask(lmask)
      ! describes, if the row and column index should be selected in mask
      ! true true   -> interface block    A_22
      ! true false  -> intermediate block A_21
      ! false true  -> intermediate block A_12
      ! false false -> interior block     A_11
      logical,intent(in)  ::  match(2)

! local variables
      integer:: ia, indi, indj
      logical:: rowmatch, colmatch

! Do the multiplication for each non-zero entry
      do ia = 1,nnz
         indi  = i_sparse(ia)
         indj  = j_sparse(ia)

         ! only if the mask is satisfied
         if (mask(indi) .ne. 0) then
            rowmatch = .true.
         else
            rowmatch = .false.
         end if
         if (mask(indj) .ne. 0) then
            colmatch = .true.
         else
            colmatch = .false.
         end if
         if (.not. ((rowmatch.eqv.match(1)) .and. (colmatch.eqv.match(2)))) then
            cycle
         end if

         vec_out(indi) = vec_out(indi) + a_sparse(ia)*vec_in(indj)
      end do

      return
end subroutine

!******************************************************************
subroutine sm_vec_mult_sym(nnz, i_sparse, j_sparse, a_sparse, la, &
                           vec_in,lvec_in, vec_out,lvec_out, &
                           mask,lmask,match)
!******************************************************************
! Subroutine for multiplication of SYMMETRIC sparse matrix in IJA format 
! with a vector vec_in, result in vec_out
! Only upper triangle of the matrix is stored
      use module_utils, only: error_exit

      implicit none

! Matrix in IJA sparse format
      integer,intent(in) :: nnz, la
      integer,intent(in) :: i_sparse(la), j_sparse(la)
      real(kr),intent(in):: a_sparse(la)

! Vectors
      integer,intent(in)  :: lvec_in, lvec_out
      real(kr),intent(in) :: vec_in(lvec_in)
      real(kr),intent(out):: vec_out(lvec_out)

! Blocks of matrix
      ! selects a block of matrix, e.g. corresponding to interface
      integer,optional,intent(in)  :: lmask
      integer,optional,intent(in)  ::  mask(:)
      ! describes, if the row and column index should be selected in mask
      ! true true   -> interface block    A_22
      ! true false  -> intermediate block A_21
      ! false true  -> intermediate block A_12
      ! false false -> interior block     A_11
      logical,optional,intent(in)  ::  match(2)

! local variables
      integer:: ia, indi, indj
      logical:: use_mask, rowmatch, colmatch, multnormal = .false., multtranspose = .false., offdiagblock=.false.

      if (present(mask)) then
         use_mask = .true.
         if (.not.(match(1).eqv.match(2))) then
            offdiagblock = .true.
         else
            offdiagblock = .false.
         end if
         if (size(mask).ne.lmask) then
            write(*,*) 'SM_VEC_MULT_FULL: Inproper length of mask.'
            call error_exit
         end if
      else
         use_mask = .false.
      end if

! Do the multiplication for each non-zero entry and its transpose
      do ia = 1,nnz
         indi  = i_sparse(ia)
         indj  = j_sparse(ia)

         ! only if the mask is satisfied
         if (use_mask) then
            if (mask(indi) .ne. 0) then
               rowmatch  = .true.
            else
               rowmatch  = .false.
            end if
            if (mask(indj) .ne. 0) then
               colmatch  = .true.
            else
               colmatch  = .false.
            end if
            multnormal    = (rowmatch.eqv.match(1)) .and. (colmatch.eqv.match(2))
            multtranspose = (colmatch.eqv.match(1)) .and. (rowmatch.eqv.match(2))
            if (.not.(multnormal .or. multtranspose)) then
               cycle
            end if
         end if

         if (.not.offdiagblock) then
            vec_out(indi) = vec_out(indi) + a_sparse(ia)*vec_in(indj)
            if (indi.ne.indj) then
               vec_out(indj) = vec_out(indj) + a_sparse(ia)*vec_in(indi)
            end if
         else
            if (multnormal) then
               vec_out(indi) = vec_out(indi) + a_sparse(ia)*vec_in(indj)
            else if (multtranspose) then
               vec_out(indj) = vec_out(indj) + a_sparse(ia)*vec_in(indi)
            else
               write(*,*) 'sm_vec_mult_sym:',' Strange combination of logical values.'
               stop
            end if
         end if
      end do

! default optional parameters
      offdiagblock = .false.

      return
end subroutine

!******************************************************************
subroutine sm_vec_mult_sym_mask(nnz, i_sparse, j_sparse, a_sparse, la, &
                                vec_in,lvec_in, vec_out,lvec_out, &
                                mask,lmask,match)
!******************************************************************
! Subroutine for multiplication of SYMMETRIC sparse matrix in IJA format 
! with a vector vec_in, result in vec_out
! Only upper triangle of the matrix is stored
      use module_utils, only: error_exit

      implicit none

! Matrix in IJA sparse format
      integer,intent(in) :: nnz, la
      integer,intent(in) :: i_sparse(la), j_sparse(la)
      real(kr),intent(in):: a_sparse(la)

! Vectors
      integer,intent(in)  :: lvec_in, lvec_out
      real(kr),intent(in) :: vec_in(lvec_in)
      real(kr),intent(out):: vec_out(lvec_out)

! Blocks of matrix
      ! selects a block of matrix, e.g. corresponding to interface
      integer,intent(in)  :: lmask
      integer,intent(in)  ::  mask(:)
      ! describes, if the row and column index should be selected in mask
      ! true true   -> interface block    A_22
      ! true false  -> intermediate block A_21
      ! false true  -> intermediate block A_12
      ! false false -> interior block     A_11
      logical,intent(in)  ::  match(2)

! local variables
      integer:: ia, indi, indj
      logical:: rowmatch, colmatch, multnormal = .false., multtranspose = .false., offdiagblock=.false.

      if (.not.(match(1).eqv.match(2))) then
         offdiagblock = .true.
      else
         offdiagblock = .false.
      end if
      if (size(mask).ne.lmask) then
         write(*,*) 'SM_VEC_MULT_FULL: Inproper length of mask.'
         call error_exit
      end if

! Do the multiplication for each non-zero entry and its transpose
      do ia = 1,nnz
         indi  = i_sparse(ia)
         indj  = j_sparse(ia)

         ! only if the mask is satisfied
         if (mask(indi) .ne. 0) then
            rowmatch  = .true.
         else
            rowmatch  = .false.
         end if
         if (mask(indj) .ne. 0) then
            colmatch  = .true.
         else
            colmatch  = .false.
         end if
         multnormal    = (rowmatch.eqv.match(1)) .and. (colmatch.eqv.match(2))
         multtranspose = (colmatch.eqv.match(1)) .and. (rowmatch.eqv.match(2))
         if (.not.(multnormal .or. multtranspose)) then
            cycle
         end if

         if (.not.offdiagblock) then
            vec_out(indi) = vec_out(indi) + a_sparse(ia)*vec_in(indj)
            if (indi.ne.indj) then
               vec_out(indj) = vec_out(indj) + a_sparse(ia)*vec_in(indi)
            end if
         else
            if (multnormal) then
               vec_out(indi) = vec_out(indi) + a_sparse(ia)*vec_in(indj)
            else if (multtranspose) then
               vec_out(indj) = vec_out(indj) + a_sparse(ia)*vec_in(indi)
            else
               write(*,*) 'sm_vec_mult_sym:',' Strange combination of logical values.'
               stop
            end if
         end if
      end do

! default optional parameters
      offdiagblock = .false.

      return
end subroutine

!****************************************************************************
subroutine sm_mat_mult(matrixtype, nnz, i_sparse, j_sparse, a_sparse, la, &
                       mat_in,lmat_in1,lmat_in2, mat_out,lmat_out1,lmat_out2)
!****************************************************************************
! Subroutine for multiplication of sparse matrix in IJA format 
! with a dense matrix mat_in, result in mat_out

      implicit none

! Type of the matrix
! 0 - unsymmetric
! 1 - symmetric positive definite
! 2 - general symmetric
      integer,intent(in) :: matrixtype

! Matrix in IJA sparse format
      integer,intent(in) :: nnz, la
      integer,intent(in) :: i_sparse(la), j_sparse(la)
      real(kr),intent(in):: a_sparse(la)

! Dense matrices
      integer,intent(in)   :: lmat_in1, lmat_in2
      real(kr),intent(in)  :: mat_in(lmat_in1, lmat_in2)
      integer,intent(in)   :: lmat_out1, lmat_out2
      real(kr),intent(out) :: mat_out(lmat_out1, lmat_out2)

! Local variables
      integer :: icol

! Check the consistency of matrices
      if (lmat_in2.ne.lmat_out2) then
         write(*,*) 'sm_mat_mult: Error - inconsistent matrices dimensions'
         stop
      end if

! Loop over columns of dense matrices
      do icol = 1,lmat_in2
         if (any(abs(mat_in(:,icol)).gt.numerical_zero)) then
            call sm_vec_mult(matrixtype, nnz, i_sparse, j_sparse, a_sparse, la, &
                             mat_in(1,icol),lmat_in1, mat_out(1,icol),lmat_out1)
         else
            mat_out(:,icol) = 0._kr
         end if
      end do

      return
end subroutine

!**************************************************************************************************************
subroutine sm_from_sm_mat_mult(matrixtype, nnz, i_sparse, j_sparse, a_sparse, la, &
                               mat_in,lmat_in1,lmat_in2, &
                               store_type,mrow,i_mout_sparse,j_mout_sparse,mout_sparse,lmout_sparse,nnz_mout, &
                               scalar)
!**************************************************************************************************************
! Subroutine for multiplication of sparse matrix in IJA format 
! with a dense matrix mat_in.
! Result is stored directly in SPARSE format.

      implicit none

! Type of the matrix
! 0 - unsymmetric
! 1 - symmetric positive definite
! 2 - general symmetric
      integer :: matrixtype

! Matrix in IJA sparse format
      integer,intent(in) :: nnz, la
      integer,intent(in) :: i_sparse(la), j_sparse(la)
      real(kr),intent(in):: a_sparse(la)

! Dense matrix
      integer,intent(in)   :: lmat_in1, lmat_in2
      real(kr),intent(in)  :: mat_in(lmat_in1, lmat_in2)

! Type of storage of result
! 0   - store all entries
! 1,2 - store upper triangle
      integer :: store_type

! Resulting matrix in IJA sparse format
      integer,intent(in)   :: mrow
      integer,intent(in)   :: lmout_sparse
      integer,intent(out)  :: i_mout_sparse(lmout_sparse), j_mout_sparse(lmout_sparse)
      real(kr),intent(out) :: mout_sparse(lmout_sparse)
      integer,intent(out)  :: nnz_mout

! Scalar coefficient that multiplies all entries of resulting matrix
      real(kr),optional,intent(in):: scalar

! Local variables
      integer ::             lvec
      real(kr),allocatable :: vec(:)
      integer :: icol, nnz_add, space_left, length
      logical :: use_scalar

! Handle optional arguments
      if (present(scalar)) then
         use_scalar = .true.
      else
         use_scalar = .false.
      end if

! Prepare vector for dense result of multiplication
      lvec = mrow
      allocate(vec(lvec))

! Loop over columns of dense matrix
      nnz_mout = 0
      do icol = 1,lmat_in2
         if (any(mat_in(:,icol).ne.0._kr)) then
            ! Find the resulting vector for one
            call sm_vec_mult(matrixtype, nnz, i_sparse, j_sparse, a_sparse, la, &
                             mat_in(1,icol),lmat_in1, vec,lvec)

            ! Convert thr vector into a sparse matrix
            select case (store_type)
               case(0)
                  length = lvec
               case(1,2)
                  length = icol
               case default
                  write(*,*) 'sm_from_sm_mat_mult: Unknown type of storage. Maybe STORE_TYPE not set.'
                  stop
            end select
            space_left = lmout_sparse - nnz_mout
            call sm_from_vec(vec,length, i_mout_sparse(nnz_mout+1), mout_sparse(nnz_mout+1), space_left, nnz_add)
            ! correct the column index
            j_mout_sparse(nnz_mout+1:nnz_mout+nnz_add) = icol
            ! Apply the scalar coefficient
            if (use_scalar) then
               mout_sparse(nnz_mout+1:nnz_mout+nnz_add) = scalar*mout_sparse(nnz_mout+1:nnz_mout+nnz_add)
            end if
            nnz_mout = nnz_mout + nnz_add
         end if
      end do

! Clear memory
      deallocate(vec)

      return
end subroutine

!*************************************************************************************************************
subroutine sm_from_sm_mat_mult_emb(matrixtype, nnz, i_sparse, j_sparse, a_sparse, la, narow, nacol, &
                                   mat_in,lmat_in1,lmat_in2, mask_row,lmask_row,mask_col,lmask_col, &
                                   store_type,i_mout_sparse,j_mout_sparse,mout_sparse,lmout_sparse,nnz_mout, &
                                   scalar)
!*************************************************************************************************************
! Subroutine for multiplication of sparse matrix in IJA format 
! with a dense matrix mat_in embedded by mask into some bigger matrix.
! Result is stored directly in sparse format.

      implicit none

! Type of the matrix
! 0 - unsymmetric
! 1 - symmetric positive definite
! 2 - general symmetric
      integer :: matrixtype

! Matrix in IJA sparse format
      integer,intent(in) :: nnz, la
      integer,intent(in) :: i_sparse(la), j_sparse(la)
      real(kr),intent(in):: a_sparse(la)
! number of rows and columns in the sparse matrix
      integer,intent(in) :: narow, nacol

! Dense matrix
      integer,intent(in)   :: lmat_in1, lmat_in2
      real(kr),intent(in)  :: mat_in(lmat_in1, lmat_in2)

! Embedding arrays
      integer,intent(in)  :: lmask_row,           lmask_col
      integer,intent(in)  ::  mask_row(lmask_row), mask_col(lmask_col)

! Type of storage of result
! 0 - store all entries
! 1,2 - store upper triangle
      integer :: store_type

! Resulting matrix in IJA sparse format
      integer,intent(in)   :: lmout_sparse
      integer,intent(out)  :: i_mout_sparse(lmout_sparse), j_mout_sparse(lmout_sparse)
      real(kr),intent(out) :: mout_sparse(lmout_sparse)
      integer,intent(out)  :: nnz_mout

! Scalar coefficient
      real(kr),optional,intent(in):: scalar

! Local variables
      integer ::             lvecin
      real(kr),allocatable :: vecin(:)
      integer ::             lvecout
      real(kr),allocatable :: vecout(:)
      integer :: icol, nnz_add, space_left, length, i
      logical :: use_scalar

! Handle optional arguments
      if (present(scalar)) then
         use_scalar = .true.
      else
         use_scalar = .false.
      end if

! Check consistency of dimensions of mask arrays
      if (lmask_row.ne.lmat_in1) then
         write(*,*) 'sm_from_sm_mat_mult_emb: Number of rows in input matrix does not match with mask array.'
         stop
      end if
      if (lmask_col.ne.lmat_in2) then
         write(*,*) 'sm_from_sm_mat_mult_emb: Number of columns in input matrix does not match with mask array.'
         stop
      end if
      if (any(mask_row.gt.nacol)) then
         write(*,*) 'sm_from_sm_mat_mult_emb: Mask of rows contains indices larger than dimension of sparse matrix.'
         stop
      end if

! Prepare vector for dense input of matrix multiplication
      lvecin = nacol
      allocate(vecin(lvecin))

! Prepare vector for dense result of matrix multiplication
      lvecout = narow
      allocate(vecout(lvecout))

! Loop over columns of dense matrices
      nnz_mout = 0
      do icol = 1,lmat_in2
         if (any(mat_in(:,icol).ne.0._kr)) then
! Find the resulting vector for one column
            ! Prepare the input vector
            vecin = 0._kr
            do i = 1,lmask_row
               vecin(mask_row(i)) = mat_in(i,icol)
            end do

            ! Do the matrix-vector multiplication
            call sm_vec_mult(matrixtype, nnz, i_sparse, j_sparse, a_sparse, la, &
                             vecin,lvecin, vecout,lvecout)

            ! Convert thr vector into a sparse matrix
            select case (store_type)
               case(0)
                  length = lvecout
               case(1,2)
                  length = mask_col(icol)
               case default
                  write(*,*) 'sm_from_sm_mat_mult_emb: Unknown type os storage. Maybe STORE_TYPE not set.'
                  stop
            end select
            space_left = lmout_sparse - nnz_mout
            call sm_from_vec(vecout,length, i_mout_sparse(nnz_mout+1), mout_sparse(nnz_mout+1), space_left, nnz_add)
            ! correct the column index
            j_mout_sparse(nnz_mout+1:nnz_mout+nnz_add) = mask_col(icol)

            ! Apply the scalar coefficient
            if (use_scalar) then
               mout_sparse(nnz_mout+1:nnz_mout+nnz_add) = scalar*mout_sparse(nnz_mout+1:nnz_mout+nnz_add)
            end if
            nnz_mout = nnz_mout + nnz_add
         end if
      end do

! Clear memory
      deallocate(vecin)
      deallocate(vecout)

      return
end subroutine

!*********************************************************************
subroutine sm_check_matrix(ibound,jbound, i_sparse, j_sparse, la, nnz)
!*********************************************************************
! Subroutine for checking that the matrix. It checks:
! - matrix is sorted
! - no index exceeds prescribed limits
! - does not contain zero rows or columns

      use module_utils
      implicit none

! Input matrix
      integer,intent(in)  :: ibound, jbound
      integer,intent(in)  :: la
      integer,intent(in) :: i_sparse(la)
      integer,intent(in) :: j_sparse(la)
      integer,intent(out) :: nnz

! Local variables
      character(*),parameter:: routine_name = 'SM_CHECK_MATRIX'
      integer :: ia, irow, jcol, indrow
      integer ::             lcol_counts
      integer, allocatable :: col_counts(:)

      ! bounds
      if (minval(i_sparse) .lt. 1 .or. maxval(i_sparse) .gt. ibound) then
         call error(routine_name,' Some row index entries out of range.')
      end if
      if (minval(j_sparse) .lt. 1 .or. maxval(j_sparse) .gt. jbound) then
         call error(routine_name,' Some column index entries out of range.')
      end if

      ! check that the array is sorted and does not contain zero rows or columns
      if (nnz.gt.0) then
         lcol_counts = jbound
         allocate(col_counts(lcol_counts))
         col_counts = 0
         indrow = 0
         do ia = 1,nnz
            irow = i_sparse(ia)
            jcol = j_sparse(ia)

            if      (irow.eq.indrow + 1) then
               ! new line starting
               indrow = irow 
            else if (irow .eq. indrow) then
               ! keeping at the same line
               continue
            else if (irow.gt.indrow+1) then
               call error(routine_name,'There seems to be a zero row in the matrix after row ',indrow)
            else
               call error(routine_name,'Matrix appears to be unsorted. Sort it prior the call to '//trim(routine_name))
            end if

            col_counts(jcol) = col_counts(jcol) + 1
         end do

         if (any(col_counts.eq.0)) then
            call error(routine_name,'There seems to be a zero column in the matrix around column ', minloc(col_counts,1))
         end if

         deallocate(col_counts)
      end if

      return
end subroutine

!************************************************************
subroutine sm_from_vec(vec,lvec, i_sparse, a_sparse, la, nnz)
!************************************************************
! Subroutine for conversion of a dense vector 
! into a sparse pair of rows and values of nonzeros.

      implicit none

! Input dense vector
      integer,intent(in) :: lvec
      real(kr),intent(in)::  vec(lvec)

! Output pair of rows and values with nonzeros
      integer,intent(in)  :: la
      integer,intent(out) :: i_sparse(la)
      real(kr),intent(out):: a_sparse(la)
      integer,intent(out) :: nnz

! Local variables
      integer :: ia, i
      real(kr):: val

      ia = 0
      do i = 1,lvec
         val = vec(i)
         if (abs(val).gt.numerical_zero) then
            ia = ia + 1
            if (ia.gt.la) then
               write(*,*) 'sm_from_vec: No space left for conversion to sparse matrix.'
               stop
            end if
            i_sparse(ia) = i
            a_sparse(ia) = val
         end if
      end do
      nnz = ia

      return
end subroutine

!**********************************************************************************
subroutine sm_from_dm(matrixtype, m,lm1,lm2, i_sparse, j_sparse, a_sparse, la, nnz)
!**********************************************************************************
! Subroutine for conversion of a dense matrix M into a sparse triplet in IJA form.

      implicit none

! Type of the matrix
! 0 - unsymmetric
! 1 - symmetric positive definite
! 2 - general symmetric
      integer :: matrixtype

! Input dense matrix
      integer,intent(in) :: lm1, lm2
      real(kr),intent(in):: m(lm1,lm2)

! Matrix in AIJ sparse format
      integer,intent(in)  :: la
      integer,intent(out) :: i_sparse(la), j_sparse(la)
      real(kr),intent(out):: a_sparse(la)
      integer,intent(out) :: nnz

! Local variables
      integer :: jcol, nnz_add, lvec, pointa, space_left

      nnz = 0
      do jcol = 1,lm2
         if (any(abs(m(:,jcol)).gt.numerical_zero)) then
            select case (matrixtype)
            case (0)
               ! matrix is stored whole
               lvec = lm1
            case (1,2)
               ! only one triangle of matrix is stored
               lvec = jcol
            case default
               write(*,*) 'sm_from_dm: Unknown type of matrix. Maybe MATRIXTYPE not set.'
               stop
            end select
            ! call routine for conversion of vector column-wise
            pointa     = nnz
            space_left = la - nnz
            call sm_from_vec(m(1,jcol),lvec, i_sparse(pointa+1), a_sparse(pointa+1), space_left, nnz_add)
            ! make column index
            j_sparse(pointa+1:pointa+nnz_add) = jcol

            ! add number of new entries
            nnz = nnz + nnz_add
         end if
      end do

      return
end subroutine

!**************************************************************************************
subroutine sm_from_dm_emb(matrixtype, m,lm1,lm2, mask_row,lmask_row,mask_col,lmask_col,&
                          i_sparse, j_sparse, a_sparse, la, nnz)
!**************************************************************************************
! Subroutine for conversion of a dense matrix M 
! embedded into larger matrix by MASK arrays into a sparse triplet IJA.

      implicit none

! Type of the matrix
! 0 - unsymmetric
! 1 - symmetric positive definite
! 2 - general symmetric
      integer :: matrixtype

! Input dense matrix
      integer,intent(in) :: lm1, lm2
      real(kr),intent(in):: m(lm1,lm2)

! Embedding matrices
      integer,intent(in)  :: lmask_row,           lmask_col
      integer,intent(in)  ::  mask_row(lmask_row), mask_col(lmask_col)

! Matrix in IJA sparse format
      integer,intent(in)  :: la
      integer,intent(out) :: i_sparse(la), j_sparse(la)
      real(kr),intent(out):: a_sparse(la)
      integer,intent(out) :: nnz

! Local variables
      integer :: irow, jcol, nnz_add, lvec, pointa, space_left

      nnz = 0
      do jcol = 1,lm2
         if (any(abs(m(:,jcol)).gt.numerical_zero)) then
            select case (matrixtype)
            case (0)
               ! matrix is stored whole
               lvec = lm1
            case (1,2)
               ! only one triangle of matrix is stored
               lvec = jcol
            case default
               write(*,*) 'sm_from_dm_emb: Unknown type of matrix. Maybe MATRIXTYPE not set.'
               stop
            end select
            ! call routine for conversion of vector column-wise
            pointa     = nnz
            space_left = la - nnz
            call sm_from_vec(m(1,jcol),lvec, i_sparse(pointa+1), a_sparse(pointa+1), space_left, nnz_add)
            ! make column index
            j_sparse(pointa+1:pointa+nnz_add) = mask_col(jcol)

            ! correct the row indices by embedding
            ! first negative indices 
            do irow = 1,lm1
               where (i_sparse(pointa+1:pointa+nnz_add).eq.irow) i_sparse(pointa+1:pointa+nnz_add) = - mask_row(irow)
            end do
            ! check if all is done well
            if (any(i_sparse(pointa+1:pointa+nnz_add).ge.0)) then
               write(*,*) 'sm_from_dm_emb: Positive entries in I_SPARSE - wrong embedding.'
               stop
            end if
            ! swap the sign
            i_sparse(pointa+1:pointa+nnz_add) = - i_sparse(pointa+1:pointa+nnz_add)

            ! add number of new entries
            nnz = nnz + nnz_add
         end if
      end do

      return
end subroutine

!********************************************************************************************************
subroutine sm_from_dmatmult(matrixtype, m,lm1,lm2, n,ln1,ln2, i_sparse,j_sparse,a_sparse,la, nnz, scalar)
!********************************************************************************************************
! Subroutine for multiplication of two dense matrices plus a scalar    scalar * M * N
! and immediate conversion of the product into a sparse triplet IJA

      implicit none

! Type of the matrix
! 0 - unsymmetric
! 1 - symmetric positive definite
! 2 - general symmetric
      integer :: matrixtype

! Input dense matrices
      integer,intent(in) :: lm1, lm2
      real(kr),intent(in):: m(lm1,lm2)
      integer,intent(in) :: ln1, ln2
      real(kr),intent(in):: n(ln1,ln2)

! Resulting matrix in sparse format IJA
      integer,intent(in)  :: la
      integer,intent(out) :: i_sparse(la), j_sparse(la)
      real(kr),intent(out):: a_sparse(la)
      integer,intent(out) :: nnz

! Scalar coefficient
      real(kr),optional,intent(in):: scalar

! Local variables
      integer :: irow, jcol, ia, start_col
      real(kr):: aval
      logical :: use_scalar

! Handle optional scalar argument
      if (present(scalar)) then
         use_scalar = .true.
      else
         use_scalar = .false.
      end if

! check the consistency of matrix dimensions
      if (lm2.ne.ln1) then
         write(*,*) 'sm_from_dmatmult: Inconsistent matrix dimensions for multiplication.'
         stop
      end if

! perform the multiplication
      ia = 0
      do irow = 1,lm1
         if (any(m(irow,:).ne.0._kr)) then
            select case (matrixtype)
               case (0)
                  ! matrix is stored whole
                  start_col = 1
               case (1,2)
                  ! only one triangle of matrix is stored
                  start_col = irow
               case default
                  write(*,*) 'sm_from_dmatmult: Unknown type of matrix. Maybe MATRIXTYPE not set.'
                  stop
            end select
            do jcol = start_col,ln2
               ! determine the entry in resulting matrix
               if (any(n(:,jcol).ne.0._kr)) then
                  aval = dot_product(m(irow,:),n(:,jcol))

                  ! apply scalar 
                  if (use_scalar) then
                     aval = scalar * aval
                  end if
                  ! if it is nonzero, add it to the triplet
                  if (abs(aval).gt.numerical_zero) then
                     ia = ia + 1
                     if (ia.gt.la) then
                        write(*,*) 'sm_from_dmatmult: No space left for conversion to sparse matrix.'
                        stop
                     end if
                     i_sparse(ia) = irow
                     j_sparse(ia) = jcol
                     a_sparse(ia) = aval
                  end if
               end if
            end do
         end if
      end do
      nnz = ia

      return
end subroutine

!************************************************************************
subroutine sm_to_dm(matrixtype, i_sparse,j_sparse,a_sparse,la, m,lm1,lm2)
!************************************************************************
! Subroutine for conversion of a matrix A stored in a sparse triplet IJA
! to dense matrix M

      implicit none

! Type of the matrix
! 0 - unsymmetric
! 1 - symmetric positive definite
! 2 - general symmetric
      integer :: matrixtype

! Input matrix in sparse format IJA
      integer,intent(in)  :: la
      integer,intent(in) :: i_sparse(la), j_sparse(la)
      real(kr),intent(in):: a_sparse(la)

! Output dense matrix
      integer,intent(in)  :: lm1, lm2
      real(kr),intent(out):: m(lm1,lm2)

! Local variables
      integer :: irow, jcol, ia

! check the consistency of matrix dimensions
      if (maxval(i_sparse).gt.lm1) then
         write(*,*) 'sm_to_dm: Inconsistent matrix dimensions - too many rows in sparse matrix.'
         stop
      end if
      if (maxval(j_sparse).gt.lm2) then
         write(*,*) 'sm_to_dm: Inconsistent matrix dimensions - too many columns in sparse matrix.'
         stop
      end if

! perform the conversion
      m = 0.0_kr
      do ia = 1,la
         irow = i_sparse(ia)
         jcol = j_sparse(ia)
         m(irow,jcol) = m(irow,jcol) + a_sparse(ia)
         if (matrixtype.eq.1.or.matrixtype.eq.2) then
            if (irow.ne.jcol) then
               m(jcol,irow) = m(jcol,irow) + a_sparse(ia)
            end if
         end if
      end do

      return
end subroutine

!****************************************************************
subroutine sm_print(iunit, i_sparse, j_sparse, a_sparse, la, nnz)
!****************************************************************
! Subroutine for printing of sparse matrix in IJA format 
! in ASCII into unit IUNIT

      implicit none

! Number of unit to print to
      integer,intent(in) :: iunit

! Matrix in IJA sparse format
      integer,intent(in) :: la
      integer,intent(in) :: i_sparse(la), j_sparse(la)
      real(kr),intent(in):: a_sparse(la)
!     number of nonzeros
      integer,intent(in),optional :: nnz

! Local variables
      integer:: nz, ia

! Number of non-zeros is optional, if it is not present, use whole lenght of A
      if (present(nnz)) then
         nz = nnz
      else
         nz = la
      end if

! Write the matrix on the unit
      write(iunit,'(i8,i8,e28.8)') (i_sparse(ia),j_sparse(ia),a_sparse(ia),ia = 1,nz)

      return
end subroutine

!**********************************
function sm_showsize(lenght,nbytes)
!**********************************
! Function returns size in megabytes of LENGHT values that have NBYTES size each

      implicit none

      real(kr):: sm_showsize

      integer,intent(in) :: lenght, nbytes
      
! determine size in megabytes
      sm_showsize = lenght*nbytes/1048576.0d0

      return
end function

end module module_sm

