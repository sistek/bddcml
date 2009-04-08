module module_adaptivity_comm
!****************************
! Module for communication instructions
! Jakub Sistek, Denver, 4/2009

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

end module module_adaptivity_comm

