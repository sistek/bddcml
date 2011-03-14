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

program test_module_bddc
! Tester of module_bddc
      use module_bddc

      implicit none
      include "mpif.h"
      integer,parameter :: kr = kind(1.D0)

! Problem dimension
      integer,parameter :: lsolt = 28, nnodt = 14, lsol = 24, nnod = 12

! Description of mesh
      integer :: lnndft = nnodt
      integer :: nndft(nnodt) = (/ 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2 /)
      integer :: lkdoft = nnodt
      integer :: kdoft(nnodt) = (/ 0, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26 /)
      integer :: lslavery = nnodt
      integer :: slavery(nnodt) = (/ 0, 0, 3, 0, 5, 0, 0, 3, 0, 5, 0, 0, 0, 0 /)
      integer :: lihntn = nnod
      integer :: ihntn(nnod) = (/ 1, 13, 7, 2, 3, 9, 4, 5, 11, 6, 14, 12 /)

! Right hand side
      integer :: lrhs = lsol
      real(kr) :: rhs(lsol) = (/ 0._kr, 0._kr, 0._kr, 0._kr, 0._kr, 0._kr, 0._kr, &
                                 0._kr, 1._kr, 1._kr, 0._kr, 0._kr, 0._kr, 0._kr, &
                                 2._kr, 2._kr, 0._kr, 0._kr, 0._kr, 0._kr, 0._kr, &
                                 0._kr, 0._kr, 0._kr /)
      integer :: lrhst = lsolt
      real(kr) :: rhst(lsolt)

! Reference solution
      real(kr) :: sol_ref_0(lsolt) = (/ 0._kr, 0._kr, 0._kr, 0._kr, 1._kr, 1._kr, 0._kr, &
                                        0._kr, 2._kr, 2._kr, 0._kr, 0._kr, 0._kr, 0._kr, &
                                        1._kr, 1._kr, 0._kr, 0._kr, 2._kr, 2._kr, 0._kr, &
                                        0._kr, 0._kr, 0._kr, 0._kr, 0._kr, 0._kr, 0._kr  /)
      real(kr) :: sol_ref_1(lsolt) = (/ 0._kr, 0._kr, 0._kr, 0._kr, 2._kr, 2._kr, 0._kr, &
                                        0._kr, 4._kr, 4._kr, 0._kr, 0._kr, 0._kr, 0._kr, &
                                        1._kr, 1._kr, 0._kr, 0._kr, 2._kr, 2._kr, 0._kr, &
                                        0._kr, 0._kr, 0._kr, 0._kr, 0._kr, 0._kr, 0._kr  /)
      real(kr) :: sol_ref_2(lsolt) = (/ 0._kr, 0._kr, 0._kr, 0._kr, 2._kr, 2._kr, 0._kr, &
                                        0._kr, 4._kr, 4._kr, 0._kr, 0._kr, 0._kr, 0._kr, &
                                        2._kr, 2._kr, 0._kr, 0._kr, 4._kr, 4._kr, 0._kr, &
                                        0._kr, 0._kr, 0._kr, 0._kr, 0._kr, 0._kr, 0._kr  /)

! rectangular matrix for testing of averages
      integer :: lavgref1, lavgref2
      real(kr), allocatable :: avgref(:,:)
! rectangular matrix for testing of averages
      integer :: lavg1, lavg2
      real(kr), allocatable :: avg(:,:)
! array of permutations
      integer :: lipiv
      integer, allocatable :: ipiv(:)
! Explicit factor L
      integer :: lfl1, lfl2
      real(kr), allocatable :: fl(:,:)
! Explicit factor U
      integer :: lfu1, lfu2
      real(kr), allocatable :: fu(:,:)
! Space for check of identity
      integer :: lidcheck1, lidcheck2
      real(kr), allocatable :: idcheck(:,:)
! Space for check of zero
      integer :: lzerocheck1, lzerocheck2
      real(kr), allocatable :: zerocheck(:,:)

! Local variables
      integer :: i
      

! Test the converting routine
      call bddc_convert_ht(nnod,nnodt,nndft,lnndft,ihntn,lihntn,slavery,lslavery,kdoft,lkdoft,rhs,lrhs,rhst,lrhst)
      write(*,*) 'Position | Solution by module | Reference solution '
      write(*,'(i6,8x, f12.7,8x, f12.7)') ( i, rhst(i), sol_ref_0(i), i = 1, lsolt)


! Test the interchange routine
      call bddc_RT(nnodt,nndft,lnndft,slavery,lslavery,kdoft,lkdoft,rhst,lrhst)

! Visual check of result
      write(*,*) 'Position | Solution by module | Reference solution '
      write(*,'(i6,8x, f12.7,8x, f12.7)') ( i, rhst(i), sol_ref_1(i), i = 1, lsolt)

! Test the interchange routine
      call bddc_R(nnodt,nndft,lnndft,slavery,lslavery,kdoft,lkdoft,rhst,lrhst)

! Visual check of result
      write(*,*) 'Position | Solution by module | Reference solution '
      write(*,'(i6,8x, f12.7,8x, f12.7)') ( i, rhst(i), sol_ref_2(i), i = 1, lsolt)

! Test the routine for LU factorization of a rectangular block
      ! reference matrix
      lavgref1 = 3
      lavgref2 = 7
      allocate(avgref(lavgref1,lavgref2))
      avgref(1,:) = (/ 0.1_kr, 0.2_kr, 0.3_kr, 1._kr, 0._kr, 0._kr, 9._kr /)
      avgref(2,:) = (/ 0.1_kr, 0.2_kr, 0.3_kr, 2._kr, 6._kr, 0._kr, 0._kr /)
      avgref(3,:) = (/ 0.1_kr, 0.2_kr, 0.3_kr, 3._kr, 4._kr, 5._kr, 0._kr /)

      lavg1 = lavgref1
      lavg2 = lavgref2
      allocate(avg(lavg1,lavg2))
      avg = avgref
      lipiv = lavg2
      allocate(ipiv(lipiv))
      write(*,*) 'Matrix AVG before run'
      do i = 1,lavg1
         write(*,'(40f13.5)') avg(i,:)
      end do
      call bddc_T_LU(avg,lavg1,lavg2,ipiv,lipiv)
      ! visual check
      write(*,*) 'Matrix AVG after run, should have LU factors in first square block'
      do i = 1,lavg1
         write(*,'(40f13.5)') avg(i,:)
      end do
      write(*,*) 'Vector of permutations'
      write(*,'(40i13)') ipiv
      lfl1 = lavg1
      lfl2 = lavg1
      allocate(fl(lfl1,lfl2))
      fl = 0._kr
      do i = 1,lavg1
        fl(i,1:i-1) = avg(i,1:i-1)
        fl(i,i)     = 1._kr
      end do
      write(*,*) 'Factor L'
      do i = 1,lavg1
         write(*,'(40f13.5)') fl(i,:)
      end do
      lfu1 = lavg1
      lfu2 = lavg2
      allocate(fu(lfu1,lfu2))
      fu = 0._kr
      do i = 1,lavg1
        fu(i,i:) = avg(i,i:)
      end do
      write(*,*) 'Factor U'
      do i = 1,lavg1
         write(*,'(40f13.5)') fu(i,:)
      end do

      call bddc_T_blinv(avg,lavg1,lavg2)
      write(*,*) 'Inversion of the first block A^-1 and B'
      do i = 1,lavg1
         write(*,'(40f13.5)') avg(i,:)
      end do

      lidcheck1 = lavg1
      lidcheck2 = lavg1
      allocate(idcheck(lidcheck1,lidcheck2))

      ! I = L*U*A^-1
      idcheck = matmul(fl,matmul(fu(:,1:lavg1),avg(:,1:lavg1)))
      write(*,*) 'Identity obtained after L*U*A^-1'
      do i = 1,lavg1
         write(*,'(40f13.5)') idcheck(i,:)
      end do

      lzerocheck1 = lavg1
      lzerocheck2 = lavg2 - lavg1
      allocate(zerocheck(lzerocheck1,lzerocheck2))

      ! 0 = L*U*-A^-1*B + B
      zerocheck = matmul(fl,matmul(fu(:,1:lavg1),avg(:,lavg1+1:))) + fu(:,lavg1+1:)
      write(*,*) 'Zero obtained after A*-A^-1*B + B'
      do i = 1,lavg1
         write(*,'(40f13.5)') zerocheck(i,:)
      end do

      deallocate(zerocheck)
      deallocate(idcheck)
      deallocate(fu)
      deallocate(fl)
      deallocate(avg)
      deallocate(avgref)
      deallocate(ipiv)

end program
