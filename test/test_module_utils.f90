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

! Tester of module_utils
program test_module_utils
      use module_utils
      use module_test
      implicit none

      integer,parameter :: kr = kind(1.D0)

      integer :: i
      integer :: ivalue, iindex
      integer :: nentries

      integer,parameter :: larray1 = 5
      integer           :: array1(larray1) = [9,8,3,5,6]
      integer           :: iindex5 = 4

      integer,parameter :: larray2 = 3
      integer           :: array2(larray2) = [8,16,6]

      integer,parameter :: lintersection = 3
      integer           ::  intersection(lintersection) 
      integer           ::  intersection_check(2) = [6,8]

      integer,parameter :: lunion = 8
      integer           ::  union(lunion) 
      integer           ::  union_check(6) = [3,5,6,8,9,16]

      integer           ::  nintersection, nunion

      integer,parameter :: larray3 = 7
      integer           :: array3(larray3) = [1,1,1,3,3,5,9]
      integer           :: array3_check(4) = [1,3,5,9]

      integer,parameter :: larray4 = 30
      integer           :: array4(larray4)        = [1,1,1,3,3,5,9,5,90,2,4,4,6,2,15,12,91,67,89,45,32,41,44,45,46,32,12,24,39,40]
      integer           :: array4_sorted(larray4) = [1,1,1,2,2,3,3,4,4,5,5,6,9,12,12,15,24,32,32,39,40,41,44,45,45,46,67,89,90,91]
      integer           :: iindex4_5 = 10

      integer,parameter :: larray5 = 6
      integer           :: array5(larray5)        = [90,9,5,1,3,4]
      integer           :: array5_sorted(larray5) = [1,3,4,5,9,90]
      integer,parameter :: larray6 = 6
      integer           :: array6(larray6)           = [1,2,3,4,5,6]
      integer           :: array6_reordered(larray6) = [4,5,6,3,2,1]

      integer           :: array7(6,3)        = reshape([90, 9, 9, 1, 3, 4, &
                                                         2,  3, 1, 2, 3, 4, &
                                                         1,  5, 3, 4, 8, 8], shape(array7))
      integer           :: array7_sorted(6,3) = reshape([1, 3, 4, 9, 9, 90, &
                                                         2, 3, 4, 1, 3, 2, &
                                                         4, 8, 8, 3, 5, 1], shape(array7))

      integer,parameter :: larray8 = 6
      integer           :: array8(larray8)           = [1,2,3,4,5,6]
      integer           :: array8_reordered(larray8) = [4,5,6,3,2,1]

      integer           :: array9(6)         = [90, 9, 9, 1, 3, 4]
      integer           :: array9_sorted(6)  = [1, 3, 4, 9, 9, 90]
      integer           :: array10(6)        = [ 2,  3, 1, 2, 3, 4]
      integer           :: array10_sorted(6) = [ 2,  3, 4, 1, 3, 2]

      real(kr)          :: array11(6,3)        = reshape([90, 9, 9, 1, 3, 4, &
                                                          6,  4, 1, 2, 3, 5, &
                                                          1,  5, 3, 4, 8, 8], shape(array7))
      real(kr)          :: array11_sorted(6,3) = reshape([ 9, 1,  3, 9, 4, 90, &
                                                           1, 2,  3, 4, 5, 6, &
                                                           3, 4,  8, 5, 8, 1], shape(array7))

      integer,parameter :: larray12 = 7
      integer           :: array12(larray12) = [1,1,1,3,3,5,9]
      integer           :: array13(larray12) = [2,3,4,1,3,2,15]
      integer           :: array12_check(4)  = [1,3,5,9]
      integer           :: array13_check(4)  = [2,1,2,15]

      integer,parameter :: larray14 = 7
      integer           :: array14(larray14) = [1,1,1,3,3,5,9]
      integer           :: array15(larray14) = [1,1,2,3,3,2,15]
      integer           :: array14_check(5)  = [1,1,3,5,9]
      integer           :: array15_check(5)  = [1,2,3,2,15]

      integer,allocatable :: array16(:) 
      integer             :: array16_check(6) = [5,4,3,2,1,0]

      integer,parameter :: larray17 = 7
      integer           :: array17(larray17)       = [4,2,1,2,3,4,-1]
      integer           :: array17_check(larray17) = [1,5,7,8,10,13,17]

      integer,parameter :: larray18 = 10 
      real(kr)          :: array18(larray18)       = &
                           [2._kr,3._kr,         4._kr,        5._kr,6._kr,             5.999999999_kr,8._kr,9._kr,10._kr,11._kr]
      real(kr)          :: array19(larray18)       = &
                           [1._kr,2.999999995_kr,4.00000002_kr,5._kr,5.9999999999999_kr,6._kr,         7._kr,8._kr,9._kr, 10._kr]
      logical           :: array1819_check(larray18) = &
                           [.false.,.false.,.true.,.false.,.false.,.false.,.false.,.false.,.false.,.false.]

      integer,parameter :: larray20 = 30 
      real(kr)          :: array20(larray20)
      real(kr)          :: array21(larray20)

      real(kr)          :: time_of_execution_wall
      real(kr)          :: time_of_execution_cpu

      ! start timing
      call time_start
      call time_start(.true.) ! measure CPU time rather than wall time

      !==============
      ! test sortings
      !==============
      ! integer array
      call iquick_sort(array4,larray4)
      call print_test_result(all(array4 == array4_sorted), 'integer quick sort')

      ! integer array with simultaneous renumbering
      call iquick_sort_simultaneous(array5,larray5, array6,larray6)
      call print_test_result(all(array5 == array5_sorted) .and. all(array6 == array6_reordered), &
                             'integer quick sort simultaneous')

      ! integer 2-dimensional array with simultaneous renumbering
      call iquick_array_sort_simultaneous(array7,size(array7,1), size(array7,2), array8,larray8)
      call print_test_result(all(array7 == array7_sorted) .and. all(array8 == array8_reordered), &
                             'integer quick sort array simultaneous')

      ! sorting two arrays first by array 1 and then by array 2
      call iquick_sort_2(array9,size(array9), array10,size(array10))
      call print_test_result(all(array9 == array9_sorted) .and. all(array10 == array10_sorted), &
                             'integer quick sort two arrays simultaneous')

      ! sorting a real array by quick sort on a given column
      call rquick_sort(array11,size(array11,1),size(array11,2), 2, 1, size(array11,1))
      call print_test_result(all(array11 == array11_sorted), &
                             'real quick sort arrays simultaneous')

      !======================
      ! test array operations
      !======================
      ! testing union of two arrays
      call get_array_union(array1,larray1,array2,larray2,union,lunion,nunion)
      call print_test_result(all(union(1:nunion) == union_check), 'array union')

      ! testing intersection of two arrays
      call get_array_intersection(array1,larray1,array2,larray2,intersection,lintersection,nintersection)
      call print_test_result(all(intersection(1:nintersection) == intersection_check), 'array intersection')

      ! test nonrepeated indices
      call get_array_norepeat(array3,larray3,nentries)
      call print_test_result(all(array3(1:nentries) == array3_check), 'nonrepeating values in array')

      ! test nonrepeated indices with the second array for permutation
      call get_array_norepeat_simultaneous(array12,size(array12),array13,size(array13),nentries)
      call print_test_result(all(array12(1:nentries) == array12_check) .and. &
                             all(array13(1:nentries) == array13_check), &
                             'nonrepeating values simultaneous')

      ! test nonrepeated pairs of indices
      call get_array_norepeat_2(array14,size(array14),array15,size(array15),nentries)
      call print_test_result(all(array14(1:nentries) == array14_check) .and. &
                             all(array15(1:nentries) == array15_check), &
                             'nonrepeating values in two arrays')

      ! push back function for integer vector like C++
      allocate(array16(0))
      call push_back(array16, 5)
      call push_back(array16, 4)
      call push_back(array16, 3)
      call push_back(array16, 2)
      call push_back(array16, 1)
      call push_back(array16, 0)
      call print_test_result(all(array16 == array16_check), 'push_back test')
      deallocate(array16)

      ! counts to starts conversion
      call counts2starts(array17,size(array17))
      call print_test_result(all(array17 == array17_check), 'counts2starts test')

      ! fuzzy less than
      call print_test_result(all(fuzzyLessThan(array18,array19) .eqv. array1819_check), 'fuzzyLessThan test')

      !==================
      ! test index search
      !==================
      ! searching index in a sorted array by binary search
      ivalue = 5
      call get_index_sorted(ivalue,array4,larray4,iindex)
      call print_test_result(iindex == iindex4_5, 'index search in sorted array')

      ! searching index of a closest value in a sorted array by binary search
      ivalue = 66
      call suppress_output_on
      call get_index_sorted(ivalue,array4,larray4,iindex)
      call print_test_result(iindex == -1, 'index search in sorted array without the value')
      call suppress_output_off

      ! searching index in general array
      ivalue = 5
      call get_index(ivalue,array1,larray1,iindex)
      call print_test_result(iindex == iindex5, 'index search') 

      !====================
      ! test random numbers
      !====================
      call initialize_random_number
      do i = 1,size(array20)
          call get_random_number(array20(i))
      end do
      call initialize_random_number
      do i = 1,size(array21)
          call get_random_number(array21(i))
      end do
      call print_test_result(all(array20 == array21), 'random number test')


      ! finish timing
      call time_end(time_of_execution_cpu)
      call time_print('the execution (CPU) was:', time_of_execution_cpu)
      call time_end(time_of_execution_wall)
      call time_print('the execution (WALL) was on proc', 0, time_of_execution_wall)

end program test_module_utils
