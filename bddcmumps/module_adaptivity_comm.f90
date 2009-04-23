module module_adaptivity_comm
!****************************
! Module for communication within LOBPCG solution of eigenproblems
! Jakub Sistek, Denver, 4/2009

! type of reals
integer,parameter,private  :: kr = kind(1.D0)
! numerical zero
real(kr),parameter,private :: numerical_zero = 1.e-12_kr

      integer :: comm_calls = 0

      integer :: comm_myid

      integer :: neigvec, problemsize

      integer :: ndofi_i, ndofi_j
      integer :: nnodci, nnodcj
      integer :: nnodi_i,nnodi_j

      integer :: comm_myplace1, comm_myplace2 
      integer :: comm_myisub, comm_myjsub 

      integer ::             lbufsend_i,    lbufsend_j
      real(kr),allocatable :: bufsend_i(:),  bufsend_j(:)
      integer ::             lbufsend   , lbufrecv
      real(kr),allocatable :: bufsend(:),  bufrecv(:)
      integer ::            lkbufsend   
      integer,allocatable :: kbufsend(:)
      integer ::            llbufa   
      integer,allocatable :: lbufa(:)

      ! array for serving to eigensolvers
      integer ::            ninstructions 
      integer ::            linstructions1
      integer,parameter ::  linstructions2 = 5
      integer,allocatable :: instructions(:,:)

      ! weigth matrix in interesting variables
      integer ::             lweight
      real(kr),allocatable :: weight(:)

      ! dependence of interface variables
      integer ::            lpairslavery
      integer,allocatable :: pairslavery(:)

      ! matrix of local constraints D_ij
      integer ::             lddij ! leading dimension
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
      integer ::             lwork
      real(kr),allocatable :: work(:)

end module module_adaptivity_comm

