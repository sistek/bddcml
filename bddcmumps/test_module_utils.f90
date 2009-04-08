program test_module_dd
! Tester of module_dd
      use module_utils
      implicit none

      integer,parameter :: kr = kind(1.D0)

      integer :: i

      character(100) :: fname
      character(100) :: name1

      ! test naming routine
      name1 = 'TESTNAME'
      do i = 1,13000,301
         call getfname(name1,i,'TEST',fname)
         write(*,*) trim(fname)
      end do

end program
