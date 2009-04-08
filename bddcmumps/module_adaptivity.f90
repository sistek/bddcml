module module_adaptivity
!***********************
! Module for adaptive search of constraints for BDDC preconditioner
! Jakub Sistek, Denver, 3/2009

! type of reals
integer,parameter,private  :: kr = kind(1.D0)
! numerical zero
real(kr),parameter,private :: numerical_zero = 1.e-12_kr

! debugging 
logical,parameter,private :: debug = .true.

! table of pairs of eigenproblems to compute
! structure:
!  PROC | ISUB | IGLBISUB | JSUB | IGLBJSUB | NVAR
integer,private            :: lpair_subdomains1 
integer,parameter,private  :: lpair_subdomains2 = 6
integer,allocatable,private :: pair_subdomains(:,:)

contains

!*************************************************
subroutine adaptivity_init(myid,comm,idpair,npair)
!*************************************************
! Subroutine for initialization of adaptive search of constraints
      use module_utils
      implicit none
      include "mpif.h"

! number of processor
      integer,intent(in) :: myid
! communicator
      integer,intent(in) :: comm
! file unit with opened list of pairs
      integer,intent(in) :: idpair

! local variables
      integer :: npair
      integer :: ierr
      integer :: ldata

      integer :: ipair, j

! read pairs to be computed data
      if (myid.eq.0) then
         read(idpair,*) npair
      end if
!*****************************************************************MPI
      call MPI_BCAST(npair,1, MPI_INTEGER, 0, comm, ierr)
!*****************************************************************MPI
      lpair_subdomains1 = npair
      allocate(pair_subdomains(lpair_subdomains1,lpair_subdomains2))
      if (myid.eq.0) then
         do ipair = 1,npair
            ! first column is associated with processors - initialize it to -1 = no processor assigned
            pair_subdomains(ipair,1) = -1
            read(idpair,*) (pair_subdomains(ipair,j), j = 2,lpair_subdomains2)
         end do
      end if
      close(idpair)
      ldata = lpair_subdomains1*lpair_subdomains2
!*****************************************************************MPI
      call MPI_BCAST(pair_subdomains,ldata, MPI_INTEGER, 0, comm, ierr)
!*****************************************************************MPI
 
      ! print what is loaded
      if (debug) then
         call adaptivity_print_pairs(myid)
      end if

      return
end subroutine

!*********************************************************
subroutine adaptivity_assign_pairs(npair,nproc,npair_locx)
!*********************************************************
! Subroutine for distribution of load to processors
      implicit none

! Global number of pairs to compute eigenproblems
      integer,intent(in) :: npair
! Number of processors
      integer,intent(in) :: nproc

! local variables
      integer :: iproc, ipair_loc, ipair, npair_locx

      npair_locx = (npair + nproc - 1)/nproc
      ipair = 0
      do iproc = 0,nproc-1
         do ipair_loc = 1,npair_locx
            ipair = ipair + 1
            if (ipair.le.npair) then
               pair_subdomains(ipair,1) = iproc
            end if
         end do
      end do
end subroutine

!******************************************************************************
subroutine adaptivity_solve_eigenvectors(myid,comm,npair_locx,npair,nsub,nproc)
!******************************************************************************
! Subroutine for parallel solution of distributed eigenproblems
      implicit none

! number of processor
      integer,intent(in) :: myid
! communicator
      integer,intent(in) :: comm

! Maximal local number of eigenproblems at one processor
      integer,intent(in) :: npair_locx

! Global number of pairs to compute eigenproblems
      integer,intent(in) :: npair

! Number of subdomains
      integer,intent(in) :: nsub

! Number of processors
      integer,intent(in) :: nproc

! local variables
      integer :: isub, jsub, ipair
      integer :: lbuf

      integer :: pair_data(lpair_subdomains2)

      ! allocate table for work instructions - the worst case is that in each
      ! round, I have to compute all the subdomains, i.e. 2 for each pair
      lintructions1 = (nsub + nproc)/nproc
      lintructions2 = 5
      allocate(instructions(lintructions1,lintructions2))

      ! loop over number of rounds for solution of eigenproblems
      do ipair_round = 1,npair_locx
         ! each round of eigenproblems has its structure - determine active pairs

         ! determine who I compute
         ! I am sending this pair
         my_pair = (myid*npair_locx) + ipair_round

         if (my_pair.le.npair) then

            call adaptivity_get_pair_data(my_pair,pair_data,lpair_data)
   
            myisub     = pair_data(2)
            myisubgl1  = pair_data(3)
            myjsub     = pair_data(4)
            myisubgl2  = pair_data(5)
            mylvec     = pair_data(6)

         end if

         ! determine working instructions for sending subdomain matrices
         ! go through pairs that are active in this round
         instructions  = 0
         ninstructions = 0
         pointbuf      = 1
         lbufrecv      = 0
         do iproc = 1,nproc
            iactive_pair = ipair_round + (iproc-1)*npair_locx
            if (iactive_pair .gt. npair) then
               cycle
            end if
            
            call adaptivity_get_pair_data(iactive_pair,pair_data,lpair_data)
            owner  = pair_data(1)
            isub   = pair_data(2)
            isubgl = pair_data(3)
            jsub   = pair_data(4)
            jsubgl = pair_data(5)
            lvec   = pair_data(6)

            ! where are these subdomains ?
            call dd_where_is_subdomain(isub,place1)
            call dd_where_is_subdomain(jsub,place2)

            if (myid.eq.place1) then
               ! add instruction
               ninstructions = ninstructions + 1

               ! who I will send data to
               instructions(ninstructions,1) = owner
               ! subdomain number
               instructions(ninstructions,2) = isub
               ! glob number
               instructions(ninstructions,3) = isubgl
            end if

            if (myid.eq.place2) then
               ! add instruction
               ninstructions = ninstructions + 1

               ! who I will send data to
               instructions(ninstructions,1) = owner
               ! subdomain number
               instructions(ninstructions,2) = jsub
               ! glob number
               instructions(ninstructions,3) = jsubgl
            end if
         end do
               
         ! scheme for this round are ready, start the game
         ! distribute cards, i.e.








      end do

      deallocate(instructions)

end subroutine

!***************************************************************
subroutine adaptivity_get_pair_data(idpair,pair_data,lpair_data)
!***************************************************************
! Subroutine for getting info about pairs to the global structure
      use utils
      implicit none

! pair number
      integer,intent(in) :: idpair
! length of vector for data
      integer,intent(in) :: lpair_data
! vector of data for pair IDPAIR
      integer,intent(in) :: pair_data(lpair_data)

! local variables
      integer :: i

      ! check the length of vector for data
      if (lpair_data .ne. lpair_subdomains2) then
         write(*,*) 'ADAPTIVITY_GET_PAIR_DATA: Size not sufficient for getting info about pair.'
         call error_exit
      end if
      ! check that the info about pair is available
      if (.not.allocated(pair_subdomains)) then
         write(*,*) 'ADAPTIVITY_GET_PAIR_DATA: Structure with global pair data is not allocated.'
         call error_exit
      end if
      if (pair_subdomains(idpair,1).eq.-1) then
         write(*,*) 'ADAPTIVITY_GET_PAIR_DATA: Incomplete information about pair - processor not assigned.'
         call error_exit
      end if
      if (any(pair_subdomains(idpair,2:).eq.0) then
         write(*,*) 'ADAPTIVITY_GET_PAIR_DATA: Incomplete information about pair - zeros in subdomain data.'
         call error_exit
      end if

      ! after checking, get info about pair from the global structure and load it
      do i = 1,lpair_data
         pair_data(i) = pair_subdomains(idpair,i)
      end do

end subroutine


!**************************************
subroutine adaptivity_print_pairs(myid)
!**************************************
! Subroutine for printing data about pairs to screen
      implicit none

! number of processor
      integer,intent(in) :: myid

! local variables
      integer :: ipair, j

      write(*,*) 'Info about loaded pairs on processor ',myid,',',lpair_subdomains1,' pairs loaded:'
      if (allocated(pair_subdomains)) then
         do ipair = 1,lpair_subdomains1
            write(*,'(6i10)') (pair_subdomains(ipair,j),j = 1,lpair_subdomains2)
         end do
      else 
         write(*,*) 'ADAPTIVITY_PRINT_PAIRS: Array of pairs is not allocated.'
      end if
end subroutine

!*****************************
subroutine adaptivity_finalize
!*****************************
! Subroutine for finalization of adaptivity
      implicit none

! clear memory
      if (allocated(pair_subdomains)) then
         deallocate(pair_subdomains)
      end if

      return
end subroutine

end module module_adaptivity

