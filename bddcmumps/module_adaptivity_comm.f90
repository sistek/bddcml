module module_adaptivity_comm
!****************************
! Module for communication within LOBPCG solution of eigenproblems
! Jakub Sistek, Denver, 4/2009

! type of reals
integer,parameter,private  :: kr = kind(1.D0)
! numerical zero
real(kr),parameter,private :: numerical_zero = 1.e-12_kr


      integer :: neigvec, problemsize

      integer :: ndofi_i, ndofi_j
      integer :: nnodci, nnodcj
      integer :: nnodi_i,nnodi_j


      integer ::             lbufsend_i,    lbufsend_j
      real(kr),allocatable :: bufsend_i(:),  bufsend_j(:)
      integer ::             lbufsend   
      real(kr),allocatable :: bufsend(:)
      integer ::            lkbufsend   
      integer,allocatable :: kbufsend(:)
      integer ::            llbufa   
      integer,allocatable :: lbufa(:)

      ! array for serving to eigensolvers
      integer ::            linstructions1
      integer,parameter ::  linstructions2 = 5
      integer,allocatable :: instructions(:,:)

      ! corner information
      integer ::             limre
      real(kr),allocatable :: imre(:)

      ! matrix of local constraints D_ij
      integer ::              ldij1,ldij2
      real(kr),allocatable ::  dij(:,:)

      ! MPI related arrays and variables
      integer ::            lrequest
      integer,allocatable :: request(:)
      integer             :: lstatarray1
      integer             :: lstatarray2
      integer,allocatable :: statarray(:,:)
      integer             :: ireq, nreq, ierr, comm_comm

      ! LAPACK QR related variables
      integer ::             ltau
      real(kr),allocatable :: tau(:)

end module module_adaptivity_comm

