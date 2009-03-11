module module_adaptivity
!***********************
! Module for adaptive search of constraints for BDDC preconditioner
! Jakub Sistek, Denver, 3/2009

! type of reals
integer,parameter,private  :: kr = kind(1.D0)
! numerical zero
real(kr),parameter,private :: numerical_zero = 1.e-12_kr

! description of pairs for eigenproblems
integer,private :: npair ! number of pairs
! table of pairs of eigenproblems to compute
integer,allocatable,private :: pair_subdomains(:,:)

contains

!*****************************************************************************************************************
subroutine adaptivity_init(myid,comm,idpair,nnod,ndof,nglb,inglb,linglb,nnglb,lnnglb,nndf,lnndf,npair,selectedglobs,lselectedglobs)
!*****************************************************************************************************************
! Subroutine for initialization of adaptive search of constraints
      implicit none
      include "mpif.h"

! local variables
      integer ::            lkglb,   lkdof
      integer,allocatable :: kglb(:), kdof(:)

      integer :: ierr

! read glob data
      if (myid.eq.0) then
         read(idpair,*) npair
      end if
!*****************************************************************MPI
      call MPI_BCAST(npair,1, MPI_INTEGER, 0, comm, ierr)
!*****************************************************************MPI
      allocate(pair_subdomains(npair,2)
      if (myid.eq.0) then
         do ipair = 1,npair
            read(idpair,*) pair_subdomains(ipair,:)
         end do
      end if
!*****************************************************************MPI
      call MPI_BCAST(pair_subdomains,2*npair, MPI_INTEGER, 0, comm, ierr)
!*****************************************************************MPI

      

! Creation of field KGLB(NGLB) with addresses before first global node of glob
      allocate(kglb(lkglb))
      kglb(1) = 0
      do iglb = 2,nglb
         kglb(iglb) = kglb(iglb-1) + nnglb(iglb-1)
      end do

! Creation of field KDOF(NNOD) with addresses before first global dof of node
      allocate(kdof(lkdof))
      kdof(1) = 0
      do inod = 2,nnod
         kdof(inod) = kdof(inod-1) + nndf(inod-1)
      end do

      allocate(globmask(lglobmask))
      do ipair = 1,npair
         indglb = selectedglobs(ipair)

         ! create mask for variables of this glob
         globmask = 0
         nglbn = nnglb(indglb)
         do iglbn = 1,nglbn
            indglbn = inglb(kglb+iglbn)
            ndofn = nndf(indglbn)
            do idofn = 1,ndofn
               globmask(kdof(indglbn) + idofn) = 1
            end do
         end do

         ! multiply by system matrix
         call sm_vec_mult(matrixtype, nnz, i_sparse, j_sparse, a_sparse, la, &
                          solintt,lsolintt, solt,lsolt, &
                          globmask,lglobmask,(/.true.,.true./))

         
      end do

      deallocate(globmask)
      deallocate(kglb)
      deallocate(kdoft)

      return
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

