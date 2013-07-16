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

module module_paraview
! Module for export to VTK ASCII XML data file loadable by ParaView

      implicit none
! adjustable parameters ############################
! type of real variables
      integer,parameter,private :: kr = kind(1.D0)
! numerical zero
      real(kr),parameter,private :: numerical_zero = 1.e-12_kr
! debugging 
      logical,parameter,private :: debug = .false.
! adjustable parameters ############################

interface paraview_write_dataarray
   module procedure paraview_write_dataarray_int
   module procedure paraview_write_dataarray_dp
   module procedure paraview_write_dataarray2d_int
   module procedure paraview_write_dataarray2d_dp
end interface paraview_write_dataarray

contains

!*********************************************************
subroutine paraview_open_subdomain_file(prefix,isub,idvtu)
!*********************************************************
! Subroutine for opening a VTU file for subdomain data
!*********************************************************
      use module_utils

      implicit none
      
      ! basename of vtu files
      character(*), intent(in) :: prefix           
      ! global subdomain index
      integer, intent(in) :: isub                  
      ! unit with opened file
      integer, intent(out) :: idvtu

      ! local vars
      character(*),parameter:: routine_name = 'PARAVIEW_OPEN_SUBDOMAIN_FILE'
      character(len=256) :: filename

      filename = ' '
      call getfname(trim(prefix),isub,'vtu',filename)
      if (debug) then
         call info(routine_name,' Opening file: '//trim(filename))
      end if
      call allocate_unit(idvtu)

      open (unit=idvtu,file=trim(filename),status='replace',form='formatted')

end subroutine paraview_open_subdomain_file

!**********************************************
subroutine paraview_close_subdomain_file(idvtu)
!**********************************************
! Subroutine for opening a VTU file for subdomain data
!**********************************************
      implicit none
      
      ! unit with opened file
      integer, intent(in) :: idvtu

      close (idvtu)

end subroutine paraview_close_subdomain_file

!************************************************
subroutine paraview_export_pvd_file(prefix, nsub)
!************************************************
! subroutine for exporting an umbrella PVD file for a set of VTU files
! for convenient handling by ParaView
      use module_utils

      implicit none

      character(*), intent(in) :: prefix           ! basename of vtu files
      integer, intent(in) :: nsub                  ! total number of subdomains

      ! local vars
      character(*),parameter:: routine_name = 'EXPORT_PVD_FILE'
      integer ::  idpvd
      character(len=256) :: filename
      integer :: isub

      filename = trim(prefix)//'.pvd'
      call info(routine_name,' Opening file: '//trim(filename))
      call allocate_unit(idpvd)

      open (unit=idpvd,file=trim(filename),status='replace',form='formatted')

      ! write header of PVD file
      write(idpvd,'(a)')             '<?xml version="1.0"?>' 
      write(idpvd,'(a)')             '<VTKFile type="Collection" version="0.1" byte_order="LittleEndian">'
      write(idpvd,'(a)')             '  <Collection>'
      ! write entries with individual VTU files
      do isub = 1,nsub
         call getfname(trim(prefix),isub,'vtu',filename)
         write(idpvd,'(a,i10,a)')    '    <DataSet part="',isub-1,'" file="./'//trim(filename)//'"/>'
      end do
      write(idpvd,'(a)')             '  </Collection>'
      write(idpvd,'(a)')             '</VTKFile>'

      close(idpvd)

end subroutine paraview_export_pvd_file

!*****************************************************************************************
subroutine paraview_write_mesh(idvtu, nelem,nnod, inet,linet, nnet,lnnet, xyz,lxyz1,lxyz2)
!*****************************************************************************************
! Subroutine for starting header and writing the mesh into it
!*****************************************************************************************
      use module_utils
      implicit none
      
      ! disk unit
      integer, intent(in) :: idvtu
      ! number of elements
      integer, intent(in) :: nelem
      ! number of nodes
      integer, intent(in) :: nnod
      ! array with mesh description
      integer, intent(in) :: linet,       lnnet
      integer, intent(in) ::  inet(linet), nnet(lnnet)
      ! array with coordinates
      integer, intent(in) :: lxyz1, lxyz2 
      real(kr), intent(in) :: xyz(lxyz1,lxyz2)
      
      ! local vars
      character(*),parameter:: routine_name = 'PARAVIEW_WRITE_MESH'
      integer :: i, ie, indinet, nne
      integer :: offset
      integer :: VTKtype



! write header of VTU file
      write(idvtu,'(a)')             '<?xml version="1.0"?>' 
      write(idvtu,'(a)')             '<VTKFile type="UnstructuredGrid" byte_order="LittleEndian">'
      write(idvtu,'(a)')             '  <UnstructuredGrid>'
      write(idvtu,'(a,i10,a,i10,a)') '  <Piece NumberOfPoints="', nnod,'" NumberOfCells="', nelem,'">'
      
      ! write nodal coordinates
      write(idvtu,'(a)')             '    <Points>'
      write(idvtu,'(a,i1,a)')        '      <DataArray type="Float64" NumberOfComponents="',lxyz2,&
                                     '" Name="Coordinates" format="ascii">'

      do i = 1,lxyz1
         write(idvtu,'(3e15.7)')       xyz(i,:)
      end do
      write(idvtu,'(a)')             '      </DataArray>'
      write(idvtu,'(a)')             '    </Points>'

      ! write cells
      write(idvtu,'(a)')             '    <Cells>'
    
      ! connectivity
      write(idvtu,'(a)')             '      <DataArray type="Int32" NumberOfComponents="1" Name="connectivity" format="ascii">'

      indinet = 0
      do ie = 1,nelem
         nne = nnet(ie) 
         write(idvtu,'(10i15)')         inet(indinet+1:indinet+nne) - 1 ! remove one to start from 0
         indinet = indinet + nne
      end do
      write(idvtu,'(a)')             '      </DataArray>'

      ! offset
      write(idvtu,'(a)')             '      <DataArray type="Int32" NumberOfComponents="1" Name="offsets" format="ascii">'
      offset = 0
      do ie = 1,nelem
         nne = nnet(ie) 
         offset = offset + nne
         write(idvtu,'(10i15)')         offset
      end do
      write(idvtu,'(a)')             '      </DataArray>'

      ! VTKtype
      write(idvtu,'(a)')             '      <DataArray type="Int32" NumberOfComponents="1" Name="types" format="ascii">'
      do ie = 1,nelem
         nne = nnet(ie)
         select case (nne) 
            case (8)
               ! hexahedron
               VTKtype = 12
            case (4)
               ! tertahedron
               VTKtype = 10
            case default
               ! default to convex point set
               VTKtype = 41
         end select
         write(idvtu,'(10i15)')         VTKtype
      end do
      write(idvtu,'(a)')             '      </DataArray>'
      write(idvtu,'(a)')             '    </Cells>'
                   
end subroutine paraview_write_mesh

!***************************************
subroutine paraview_open_celldata(idvtu)
!***************************************
! Subroutine for opening celldata in VTU file
!***************************************
      implicit none
      
      ! unit with opened file
      integer, intent(in) :: idvtu

      write(idvtu,'(a)') '    <CellData>'

end subroutine paraview_open_celldata

!***************************************
subroutine paraview_close_celldata(idvtu)
!***************************************
! Subroutine for closing celldata in VTU file
!***************************************
      implicit none
      
      ! unit with opened file
      integer, intent(in) :: idvtu

      write(idvtu,'(a)') '    </CellData>'

end subroutine paraview_close_celldata

!****************************************
subroutine paraview_open_pointdata(idvtu)
!****************************************
! Subroutine for opening pointdata in VTU file
!****************************************
      implicit none
      
      ! unit with opened file
      integer, intent(in) :: idvtu

      write(idvtu,'(a)') '    <PointData>'

end subroutine paraview_open_pointdata

!*****************************************
subroutine paraview_close_pointdata(idvtu)
!*****************************************
! Subroutine for closing pointdata in VTU file
!*****************************************
      implicit none
      
      ! unit with opened file
      integer, intent(in) :: idvtu

      write(idvtu,'(a)') '    </PointData>'

end subroutine paraview_close_pointdata

!******************************************************************************************
subroutine paraview_write_dataarray_int(idvtu,number_of_components,array_name,array,larray)
!******************************************************************************************
! Subroutine for writing a one-dimensional integer array into a VTU file
!******************************************************************************************
      use module_utils
      implicit none
      
      ! unit with opened file
      integer, intent(in) :: idvtu
      ! how many components are in the array (e.g. 1 for scalar, 3 for vectors in 3D)
      integer, intent(in) :: number_of_components
      ! name of the array
      character(*), intent(in) :: array_name           
      ! the array
      integer, intent(in) :: larray
      integer, intent(in) :: array(larray)

      ! local vars
      real(kr) :: realaux(1)

      call paraview_write_dataarray_generic(idvtu,number_of_components,array_name,array,realaux,larray,.true.)

end subroutine paraview_write_dataarray_int

!*****************************************************************************************
subroutine paraview_write_dataarray_dp(idvtu,number_of_components,array_name,array,larray)
!*****************************************************************************************
! Subroutine for writing a one-dimensional double precision array into a VTU file
!*****************************************************************************************
      use module_utils
      implicit none
      
      ! unit with opened file
      integer, intent(in) :: idvtu
      ! how many components are in the array (e.g. 1 for scalar, 3 for vectors in 3D)
      integer, intent(in) :: number_of_components
      ! name of the array
      character(*), intent(in) :: array_name           
      ! the array
      integer, intent(in) :: larray
      real(kr), intent(in) :: array(larray)

      ! local vars
      integer :: intaux(1)

      call paraview_write_dataarray_generic(idvtu,number_of_components,array_name,intaux,array,larray,.false.)

end subroutine paraview_write_dataarray_dp

!*****************************************************************************************************
subroutine paraview_write_dataarray2d_int(idvtu,number_of_components,array_name,array,larray1,larray2)
!*****************************************************************************************************
! Subroutine for writing a two-dimensional integer array into a VTU file
!*****************************************************************************************************
      use module_utils
      implicit none
      
      ! unit with opened file
      integer, intent(in) :: idvtu
      ! how many components are in the array (e.g. 1 for scalar, 3 for vectors in 3D)
      integer, intent(in) :: number_of_components
      ! name of the array
      character(*), intent(in) :: array_name           
      ! the array
      integer, intent(in) :: larray1, larray2
      integer, intent(in) :: array(larray1,larray2)

      ! local vars
      integer  :: larray
      real(kr) :: realaux(1)
      integer  :: new_shape(2)

      ! set sizes of reshaped array to make the array one-dimensional
      larray = larray1 * larray2
      new_shape(1) = larray
      new_shape(2) = 1

      call paraview_write_dataarray_generic(idvtu,number_of_components,array_name,reshape(array,new_shape),realaux,larray,.true.)

end subroutine paraview_write_dataarray2d_int

!****************************************************************************************************
subroutine paraview_write_dataarray2d_dp(idvtu,number_of_components,array_name,array,larray1,larray2)
!****************************************************************************************************
! Subroutine for writing a two-dimensional double precision array into a VTU file
!****************************************************************************************************
      use module_utils
      implicit none
      
      ! unit with opened file
      integer, intent(in) :: idvtu
      ! how many components are in the array (e.g. 1 for scalar, 3 for vectors in 3D)
      integer, intent(in) :: number_of_components
      ! name of the array
      character(*), intent(in) :: array_name           
      ! the array
      integer, intent(in) :: larray1, larray2
      real(kr), intent(in) :: array(larray1,larray2)

      ! local vars
      integer  :: larray
      integer :: intaux(1)
      integer :: new_shape(2)

      ! set sizes of reshaped array to make the array one-dimensional
      larray = larray1 * larray2
      new_shape(1) = larray
      new_shape(2) = 1

      call paraview_write_dataarray_generic(idvtu,number_of_components,array_name,intaux,reshape(array,new_shape),larray,.false.)

end subroutine paraview_write_dataarray2d_dp

!**********************************************************************************************************************
subroutine paraview_write_dataarray_generic(idvtu,number_of_components,array_name,arrayint,arrayreal,larray,is_integer)
!**********************************************************************************************************************
! Generic subroutine for writing an array into a VTU file
!**********************************************************************************************************************
      use module_utils
      implicit none
      
      ! unit with opened file
      integer, intent(in) :: idvtu
      ! how many components are in the array (e.g. 1 for scalar, 3 for vectors in 3D)
      integer, intent(in) :: number_of_components
      ! name of the array
      character(*), intent(in) :: array_name           
      ! the array
      integer, intent(in) :: larray
      integer, intent(in) :: arrayint(larray)
      real(kr),intent(in) :: arrayreal(larray)
      ! if true, integer array will be used, otherwise the real array. The other array is not accessed.
      logical, intent(in) :: is_integer

      ! local vars
      character(*),parameter:: routine_name = 'PARAVIEW_WRITE_DATAARRAY'
      character(10) :: num_type
      integer :: icol, i
      integer :: npoint

      ! set type of numbers
      num_type = ''
      if (is_integer) then
         num_type = 'Int32'
      else
         num_type = 'Float64'
      end if

      write(idvtu,'(a,i1,a)')  '      <DataArray type="'//trim(num_type)//'" NumberOfComponents="',number_of_components,&
                                       '" Name="'//trim(array_name)//'" format="ascii">'
      npoint = larray/number_of_components
      if (mod(larray,npoint).ne.0) then
         call error( routine_name, 'Length of dataarray mismatch. Expected length: ', npoint * number_of_components )
      end if
      do i = 1,npoint
         if (is_integer) then
            write(idvtu,'(3i15)')   ( arrayint( (icol-1)*npoint + i ), icol = 1,number_of_components )
         else
            write(idvtu,'(3e15.7)') ( arrayreal( (icol-1)*npoint + i ), icol = 1,number_of_components )
         end if
      end do
      write(idvtu,'(a)')       '      </DataArray>'

end subroutine paraview_write_dataarray_generic

!***************************************
subroutine paraview_finalize_file(idvtu)
!***************************************
! Subroutine for finalizing a VTU file
!***************************************
      implicit none
      
      ! unit with opened file
      integer, intent(in) :: idvtu

      write(idvtu,'(a)')             '  </Piece>'
      write(idvtu,'(a)')             ' </UnstructuredGrid>'
      write(idvtu,'(a)')             '</VTKFile>'
    
end subroutine paraview_finalize_file

end module module_paraview

