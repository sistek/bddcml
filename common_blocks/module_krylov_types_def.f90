module module_krylov_types_def
!*****************************
! Module for definition of vectors for Krylov method
! Jakub Sistek, Praha 6/2010

      implicit none

! type of real variables
      integer,parameter,private :: kr = kind(1.D0)

! type for auxiliary handling of data for Krylov methods in distributed manner
      type common_krylov_data_type
         integer ::            lvec_in
         real(kr),pointer ::    vec_in(:)  ! vector that inputs SM or PC applications
         integer ::            lvec_out
         real(kr),pointer ::    vec_out(:) ! vector that outputs SM or PC applications
      end type common_krylov_data_type

! type for storing data for PCG method in distributed manner
      type pcg_data_type
         integer ::             lsoli
         real(kr),allocatable :: soli(:)     ! array of solution at interface
         integer ::             lresi
         real(kr),allocatable :: resi(:)     ! array of residual at interface
         integer ::             lp
         real(kr),allocatable :: p(:)        ! array for search direction p
         integer ::             lap
         real(kr),allocatable :: ap(:)       ! array for A*p
         integer ::             lz
         real(kr),allocatable :: z(:)        ! array for preconditioned residual M*res
      end type pcg_data_type

! type for storing data for BICGSTAB method in distributed manner
      type bicgstab_data_type
         integer ::             lsoli
         real(kr),allocatable :: soli(:)     ! array of solution at interface
         integer ::             lresi
         real(kr),allocatable :: resi(:)     ! array of residual at interface
         integer ::             lresistab
         real(kr),allocatable :: resistab(:) ! array of stabilizing residual at interface
         integer ::             lv
         real(kr),allocatable :: v(:)     
         integer ::             lp
         real(kr),allocatable :: p(:)     
         integer ::             ly
         real(kr),allocatable :: y(:)        
         integer ::             lz
         real(kr),allocatable :: z(:)        
         integer ::             ls
         real(kr),allocatable :: s(:)       
         integer ::             lt
         real(kr),allocatable :: t(:)      
      end type bicgstab_data_type

end module
