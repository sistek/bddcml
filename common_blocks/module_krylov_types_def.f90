module module_krylov_types_def
!*****************************
! Module for definition of vectors for Krylov method
! Jakub Sistek, Praha 6/2010

      implicit none

! type of real variables
      integer,parameter,private :: kr = kind(1.D0)

! type for storing PCG data in distributed manner
      type pcg_data_type
         integer ::             lsoli
         real(kr),allocatable :: soli(:)     ! array of solution at interface
         integer ::             lresi
         real(kr),allocatable :: resi(:)     ! array of residual at interface
         real(kr),allocatable :: resiadj(:)  ! array of residual at interface from adjacent subdomains
         integer ::             lp
         real(kr),allocatable :: p(:)        ! array for search direction p
         real(kr),allocatable :: padj(:)     ! array for search direction p from adjacent subdomains
         integer ::             lap
         real(kr),allocatable :: ap(:)       ! array for A*p
         real(kr),allocatable :: apadj(:)    ! array for A*p from adjacent subdomains
         integer ::             lz
         real(kr),allocatable :: z(:)        ! array for preconditioned residual M*res
         real(kr),allocatable :: zadj(:)     ! array for preconditioned residual M*res from adjacent subdomains
      end type pcg_data_type

end module
