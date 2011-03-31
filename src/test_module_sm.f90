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

program test_module_sm
! Tester of module_sm
      use module_sm

      implicit none
      integer,parameter :: kr = kind(1.D0)

! Problem dimension
      integer,parameter :: lsol = 6, nnod = 6, linet = 8, nelem = 2, nnodi = 2

! Description of mesh
      integer :: inet(linet) = (/ 1, 2, 5, 4, 2, 3, 6, 5 /)
      integer :: lnnet = nelem
      integer :: nnet(nelem) = (/ 4, 4 /)
      integer :: lnndf = nnod
      integer :: nndf(nnod) = (/ 1, 1, 1, 1, 1, 1 /)
      integer :: lkdof = nnod
      integer :: kdof(nnod) = (/ 0, 1, 2, 3, 4, 5 /)

! Description of interface
      integer :: liingn = nnodi
      integer ::  iingn(nnodi) = (/ 2, 5 /)
      

! Right hand side
      integer :: lrhs = lsol
      real(kr) :: rhs(lsol) = (/ 1._kr, 2._kr, 1._kr, 1._kr, 2._kr, 1._kr /)

! Fixed variables
      integer :: lifix = lsol
      integer :: ifix(lsol) = (/ 1, 0, 0, 1, 0, 0 /)
      integer :: lfixv = lsol
      real(kr) :: fixv(lsol) = (/ 1._kr, 0._kr, 0._kr, 1._kr, 0._kr, 0._kr /)

! Element matrix
      integer,parameter :: lelm = 10
      real(kr) :: elm(lelm) = (/ 4._kr, 1._kr, 4._kr, 0._kr, 1._kr, 4._kr, 1._kr, 0._kr, 1._kr, 4._kr  /)

! Mask
      integer :: lmask = lsol
      integer ::  mask(lsol) = (/ 0, 1, 0, 0, 1, 0 /)

! Unit for disk file with element matrices
      integer :: idelm = 1

! Reference solution
      integer:: i
      real(kr) :: sol_ref_1(lsol) = (/ 6._kr, 12._kr, 6._kr, 6._kr, 12._kr, 6._kr /)
      real(kr) :: sol_ref_2(lsol) = (/ 6._kr, 11._kr, 5._kr, 5._kr,  9._kr, 4._kr /)
      real(kr) :: sol_ref_3(lsol) = (/ 4._kr,  1._kr, 1._kr, 4._kr,  1._kr, 1._kr /)
      real(kr) :: sol_ref_4(lsol) = (/ 5._kr,  0._kr, 5._kr, 5._kr,  0._kr, 5._kr /)
      real(kr) :: sol_ref_5(lsol) = (/ 0._kr,  2._kr, 0._kr, 0._kr,  2._kr, 0._kr /)
      real(kr) :: sol_ref_6(lsol) = (/ 1._kr,  0._kr, 1._kr, 1._kr,  0._kr, 1._kr /)
      real(kr) :: sol_ref_7(lsol) = (/ 0._kr, 10._kr, 0._kr, 0._kr, 10._kr, 0._kr /)

! Full matrix to be converted to sparse
      integer,parameter :: lm1 = 5, lm2 = 5
      real(kr) ::          m(lm1,lm2)

! Matrix AIJ in sparse format
      integer :: la
      integer :: nnz
      real(kr),allocatable :: a_sparse(:)
      integer,allocatable  :: i_sparse(:)
      integer,allocatable  :: j_sparse(:)

! Matrix AIJ in sparse format - quarters
      integer :: la11,  la12,  la21,  la22

! Vector for matrix vector testing
      integer :: lones = lsol
      real(kr) :: ones(lsol) = 1.0

! Vector of boudary conditions
      integer :: lbc = lsol
      real(kr) :: bc(lsol)

! Vector of solution
      real(kr) :: sol(lsol) 

! Local variables
      integer :: matrixtype = 1
      integer :: im

      
! Prepare virtual file with element matrices
      open(unit=idelm, form='unformatted', status='scratch')
      write(idelm) lelm, (elm(i),i = 1,lelm)
      write(idelm) lelm, (elm(i),i = 1,lelm)
      rewind idelm
   
! Load sparse matrix - using matrix quarters
      call sm_pmd_get_length_quart(matrixtype,nelem,inet,linet,nnet,lnnet,nndf,lnndf,iingn,liingn, la11,la12,la21,la22)
! Print the matrix lengths
      write(*,*) 'Matrix lengths ...'
      write(*,*) ' la11 = ', la11
      write(*,*) ' la12 = ', la12
      write(*,*) ' la21 = ', la21
      write(*,*) ' la22 = ', la22
      

! Load sparse matrix - usual routine
      call sm_pmd_get_length(matrixtype,nelem,inet,linet,nnet,lnnet,nndf,lnndf,la)
! Prepare memory for sparse matrix
      allocate(i_sparse(la), j_sparse(la), a_sparse(la))
! Load sparse matrix
      call sm_pmd_load(matrixtype,idelm,nelem,inet,linet,nnet,lnnet,nndf,lnndf,kdof,lkdof,&
                       i_sparse, j_sparse, a_sparse, la)
      close(idelm)

! Print the matrix
      write(*,*) 'Matrix after loading...'
      call sm_print(6, i_sparse, j_sparse, a_sparse, la)

! Sort the matrix
      call sm_sort(i_sparse, j_sparse, a_sparse, la)
      write(*,*) 'Matrix after plain sorting...'
      call sm_print(6, i_sparse, j_sparse, a_sparse, la)

! Assembly the matrix
      call sm_assembly(i_sparse, j_sparse, a_sparse, la, nnz)

! Print the matrix
      write(*,*) 'Matrix after assembly...'
      call sm_print(6, i_sparse, j_sparse, a_sparse, la, nnz)

! Apply symmetric matrix vector multiplication
      call sm_vec_mult(matrixtype, nnz, i_sparse, j_sparse, a_sparse, la, &
                       ones,lones, sol,lsol)

! Visual check of result
      write(*,*) 'Position | Solution by module | Reference solution '
      write(*,'(i6,8x, f12.7,8x, f12.7)') ( i, sol(i), sol_ref_1(i), i = 1, lsol)

! Apply full matrix vector multiplication
      matrixtype = 0
      call sm_vec_mult(matrixtype, nnz, i_sparse, j_sparse, a_sparse, la, &
                       ones,lones, sol,lsol)

! Visual check of result
      write(*,*) 'Position | Solution by module | Reference solution '
      write(*,'(i6,8x, f12.7,8x, f12.7)') ( i, sol(i), sol_ref_2(i), i = 1, lsol)

      write(*,*) 'Check the masked array multiplication'
      matrixtype = 1
      call sm_vec_mult_mask(matrixtype, nnz, i_sparse, j_sparse, a_sparse, la, &
                            ones,lones, sol,lsol, mask,lmask, (/.false.,.false./))
! Visual check of result
      write(*,*) 'Position | Solution by module | Reference solution '
      write(*,'(i6,8x, f12.7,8x, f12.7)') ( i, sol(i), sol_ref_4(i), i = 1, lsol)

      call sm_vec_mult_mask(matrixtype, nnz, i_sparse, j_sparse, a_sparse, la, &
                            ones,lones, sol,lsol, mask,lmask, (/.true.,.false./))
! Visual check of result
      write(*,*) 'Position | Solution by module | Reference solution '
      write(*,'(i6,8x, f12.7,8x, f12.7)') ( i, sol(i), sol_ref_5(i), i = 1, lsol)

      call sm_vec_mult_mask(matrixtype, nnz, i_sparse, j_sparse, a_sparse, la, &
                            ones,lones, sol,lsol, mask,lmask, (/.false.,.true./))
! Visual check of result
      write(*,*) 'Position | Solution by module | Reference solution '
      write(*,'(i6,8x, f12.7,8x, f12.7)') ( i, sol(i), sol_ref_6(i), i = 1, lsol)

      call sm_vec_mult_mask(matrixtype, nnz, i_sparse, j_sparse, a_sparse, la, &
                            ones,lones, sol,lsol, mask,lmask, (/.true.,.true./))
! Visual check of result
      write(*,*) 'Position | Solution by module | Reference solution '
      write(*,'(i6,8x, f12.7,8x, f12.7)') ( i, sol(i), sol_ref_7(i), i = 1, lsol)

! Apply boundary conditions 
      fixv = 1.0_kr
      call sm_apply_bc(matrixtype,ifix,lifix,fixv,lfixv,i_sparse,j_sparse,a_sparse,nnz,bc,lbc) 
      write(*,*) 'Matrix after application of BC...'
      call sm_print(6, i_sparse, j_sparse, a_sparse, la, nnz)

! Prepare RHS
      call sm_prepare_rhs(ifix,lifix,bc,lbc,rhs,lrhs)

! Visual check of result
      write(*,*) 'Position | Solution by module | Reference solution '
      write(*,'(i6,8x, f12.7,8x, f12.7)') ( i, rhs(i), sol_ref_3(i), i = 1, lsol)

! Filter out the zeros from matrix
      call sm_assembly(i_sparse, j_sparse, a_sparse, la, nnz)
      write(*,*) 'Matrix after filtering out zeros...'
      call sm_print(6, i_sparse, j_sparse, a_sparse, la, nnz)

! Clear memory
      deallocate(i_sparse, j_sparse, a_sparse)

      write(*,*) 'Test of conversion of dense matrix into sparse...'
      m(1,:) = (/ 0._kr, 1._kr, 0._kr, 3._kr, 0._kr /)
      m(2,:) = (/ 0._kr, 0._kr, 1._kr, 0._kr, 0._kr /)
      m(3,:) = (/ 0._kr, 0._kr, 0._kr, 1._kr, 0._kr /)
      m(4,:) = (/ 0._kr, 0._kr, 0._kr, 0._kr, 1._kr /)
      m(5,:) = (/ 1._kr, 0._kr, 0._kr, 0._kr, 1._kr /)
      write(*,*) 'Dense matrix...'
      do im = 1,lm1
         write(*,'(40f13.7)') m(im,:)
      end do

! Convert matrix to sparse
      la = count(m.ne.0._kr)
      allocate(i_sparse(la), j_sparse(la), a_sparse(la))
      matrixtype = 0
      call sm_from_dm(matrixtype, m,lm1,lm2, i_sparse, j_sparse, a_sparse, la, nnz)
      write(*,*) 'And its full sparse representation...' 
      call sm_print(6, i_sparse, j_sparse, a_sparse, la, nnz)
      i_sparse = 0
      j_sparse = 0
      a_sparse = 0._kr
      matrixtype = 1
      call sm_from_dm(matrixtype, m,lm1,lm2, i_sparse, j_sparse, a_sparse, la, nnz)
      write(*,*) 'Full sparse representation of upper triangle...' 
      call sm_print(6, i_sparse, j_sparse, a_sparse, la, nnz)
      deallocate(i_sparse, j_sparse, a_sparse)

end program
