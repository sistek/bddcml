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

module module_krylov_types_def
!*****************************
! Module for definition of vectors for Krylov method
! Jakub Sistek, Praha 6/2010

      use, intrinsic :: iso_fortran_env
      implicit none

! type of real variables
      integer,parameter,private :: kr = REAL64

! type for auxiliary handling of data for Krylov methods in distributed manner
      type common_krylov_data_type
         integer ::               lvec_in
         complex(kr),pointer ::    vec_in(:)  ! vector that inputs SM or PC applications
         integer ::               lvec_out
         complex(kr),pointer ::    vec_out(:) ! vector that outputs SM or PC applications
      end type common_krylov_data_type

! type for storing data for PCG method in distributed manner
      type pcg_data_type
         integer ::                lsoli
         complex(kr),allocatable :: soli(:)     ! array of solution at interface
         integer ::                lresi
         complex(kr),allocatable :: resi(:)     ! array of residual at interface
         integer ::                lp
         complex(kr),allocatable :: p(:)        ! array for search direction p
         integer ::                lap
         complex(kr),allocatable :: ap(:)       ! array for A*p
         integer ::                lz
         complex(kr),allocatable :: z(:)        ! array for preconditioned residual M*res
      end type pcg_data_type

! type for storing data for BICGSTAB method in distributed manner
      type bicgstab_data_type
         integer ::                lsoli
         complex(kr),allocatable :: soli(:)     ! array of solution at interface
         integer ::                lresi
         complex(kr),allocatable :: resi(:)     ! array of residual at interface
         integer ::                lresistab
         complex(kr),allocatable :: resistab(:) ! array of stabilizing residual at interface
         integer ::                lv
         complex(kr),allocatable :: v(:)
         integer ::                lp
         complex(kr),allocatable :: p(:)
         integer ::                ly
         complex(kr),allocatable :: y(:)
         integer ::                lz
         complex(kr),allocatable :: z(:)
         integer ::                ls
         complex(kr),allocatable :: s(:)
         integer ::                lt
         complex(kr),allocatable :: t(:)
      end type bicgstab_data_type

! type for storing data for steepest descent iterative method in distributed manner
      type steepestdescent_data_type
         integer ::                lsoli
         complex(kr),allocatable :: soli(:)     ! array of solution at interface
         integer ::                lresi
         complex(kr),allocatable :: resi(:)     ! array of residual at interface
         integer ::                lz
         complex(kr),allocatable :: z(:)        ! M * residual
         integer ::                laz
         complex(kr),allocatable :: az(:)       ! array for A*u
      end type steepestdescent_data_type


! type for storing data for recycling of Krylov subspaces
      type krylov_recycling_data_type
         integer ::                lv1
         integer ::                lv2
         complex(kr),allocatable :: v(:,:)   ! array of all search directions p on interface
         integer ::                lw1
         integer ::                lw2
         complex(kr),allocatable :: w(:,:)   ! array of all A*p at interface
      end type krylov_recycling_data_type

end module

