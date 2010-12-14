program test_module_utils
! Tester of module_utils
      use module_utils
      implicit none

      integer,parameter :: kr = kind(1.D0)

      integer :: i, ivalue, iindex
      integer :: nentries

      integer,parameter :: larray1 = 5
      integer           :: array1(larray1) = (/9,8,3,5,6/)
      integer,parameter :: larray2 = 3
      integer           :: array2(larray2) = (/8,16,6/)
      integer,parameter :: lintersection = 3
      integer           ::  intersection(lintersection) 
      integer,parameter :: lunion = 8
      integer           ::  union(lunion) 
      integer           ::  nintersection, nunion
      integer,parameter :: larray3 = 7
      integer           :: array3(larray3) = (/1,1,1,3,3,5,9/)
      integer,parameter :: larray4 = 30
      integer           :: array4(larray4) = (/1,1,1,3,3,5,9,5,90,2,4,4,6,2,15,12,91,67,89,45,32,41,44,45,46,32,12,24,39,40/)

      real(kr) :: rnd

      character(100) :: fname
      character(100) :: name1

      ! test naming routine
      name1 = 'TESTNAME'
      do i = 1,13000,301
         call getfname(name1,i,'TEST',fname)
         write(*,*) trim(fname)
      end do

      ! sorting
      call get_array_intersection(array1,larray1,array2,larray2,intersection,lintersection,nintersection)
      print *, 'Test of computing intersection of 2 integer arrays.'
      print *, 'array1'
      print *,  array1
      print *, 'array2'
      print *,  array2
      print *, 'intersection'
      print *,  intersection
      print *, 'valid number'
      print *,  nintersection
      call get_array_union(array1,larray1,array2,larray2,union,lunion,nunion)
      print *, 'Test of computing union of 2 integer arrays.'
      print *, 'array1'
      print *,  array1
      print *, 'array2'
      print *,  array2
      print *, 'union'
      print *,  union
      print *, 'valid number'
      print *,  nunion

      ! searching index in array
      ivalue = 5
      call get_index(ivalue,array1,larray1,iindex)
      print *,'Index of value ',ivalue,' in array' 
      print *,  array1
      print *,'is ',iindex

      call iquick_sort(array4,larray4)
      call get_index_sorted(ivalue,array4,larray4,iindex)
      print *,'Index of value ',ivalue,' in sorted array'
      print *,  array4
      print *,'is ',iindex
      ivalue = 66
      call get_index_sorted(ivalue,array4,larray4,iindex)
      print *,'Index of value ',ivalue,' in sorted array'
      print *,  array4
      print *,'is ',iindex

      ! test nonrepeated indices
      print *, 'Test of removing repeated entries in array:'
      print *, 'array before'
      print *,  array3
      call get_array_norepeat(array3,larray3,nentries)
      print *, 'array after'
      print *,  array3
      print *, 'nentries =', nentries

      ! random numbers
      call get_random_number(rnd)
      write(*,*) 'rnd = ',rnd
      call get_random_number(rnd)
      write(*,*) 'rnd = ',rnd
      call get_random_number(rnd)
      write(*,*) 'rnd = ',rnd

end program
