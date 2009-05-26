program test_module_utils
! Tester of module_utils
      use module_utils
      implicit none

      integer,parameter :: kr = kind(1.D0)

      integer :: i, ivalue, iindex

      integer,parameter :: larray1 = 5
      integer           :: array1(larray1) = (/9,8,3,5,6/)
      integer,parameter :: larray2 = 3
      integer           :: array2(larray2) = (/8,16,6/)
      integer,parameter :: lintersection = 3
      integer           ::  intersection(lintersection) 
      integer,parameter :: lunion = 8
      integer           ::  union(lunion) 
      integer           ::  nintersection, nunion

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

end program
