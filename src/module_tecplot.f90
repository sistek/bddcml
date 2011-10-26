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

module module_tecplot
! Module for export to Tecplot ASCII data file
implicit none
integer,parameter,private :: kr = kind(1.D0)

interface tecplot_export_block_variable
  module procedure tecplot_export_block_variable_sp
  module procedure tecplot_export_block_variable_dp
end interface tecplot_export_block_variable

contains

!***************************************************
subroutine tecplot_header(iddat,ndim,name1,varnames)
!***************************************************
! Subroutine for export header to TECPLOT file
!***************************************************
implicit none

! disk unit
integer, intent(in) :: iddat
! number of dimensions
integer, intent(in) :: ndim
! problem name
character(*),intent(in)    :: name1
! variable names (excluding coordinates)
character(*),optional,intent(in) :: varnames

! local variables
character(13) :: coordstring

! write title
write(iddat,*) 'TITLE = "',trim(name1),'"'

! write names of variables
if      (ndim.eq.3) then
        coordstring = '"X", "Y", "Z"'
else if (ndim.eq.2) then
        coordstring = '"X", "Y"'
else 
        write(*,*) 'TECPLOT_HEADER: Strange number of dimensions. ndim =',ndim
        stop
end if
if (present(varnames) .and. varnames.ne.' ') then
   write(iddat,*) 'VARIABLES = '//trim(coordstring)//', '//trim(varnames)
else
   write(iddat,*) 'VARIABLES = '//trim(coordstring)
end if

return
end subroutine tecplot_header

!************************************************************************************
subroutine tecplot_start_fe_zone(iddat,ndim,nnod,nelem,nvar,cellcentered,datapacking)
!************************************************************************************
! Subroutine for starting a new finite element zone in TECPLOT data file
!************************************************************************************
implicit none

! disk unit
integer, intent(in) :: iddat
! space dimensions
integer, intent(in) :: ndim
! number of nodes
integer, intent(in) :: nnod
! number of elements
integer, intent(in) :: nelem
! number of variables to export
integer, intent(in) :: nvar
! is this variable cellcentered?
logical :: cellcentered(nvar)
! type of data packing
!  0 - block
!  1 - point
integer, intent(in) :: datapacking

! local variables
character(15)  :: zonetype
character(5)   :: datapackingst
character(1000):: varlocation
character(2)   :: ivarst

integer :: ivar, i

! Type of elements
if      (ndim.eq.3) then
   zonetype = 'FEBRICK'
else if (ndim.eq.2) then
   zonetype = 'FEQUADRILATERAL'
else
   write(*,*) 'TECPLOT_START_FE_ZONE: Strange number of dimensions, ndim =',ndim
   stop
end if

! Type of packing data
if      (datapacking.eq.0) then
   datapackingst = 'BLOCK'
else if (datapacking.eq.1) then
   datapackingst = 'POINT'
else
   write(*,*) 'TECPLOT_START_FE_ZONE: Strange datapacking, datapacking =',datapacking
   stop
end if

! Type of storing variables
if (any(cellcentered)) then
   varlocation = ' VARLOCATION = ('
   do i = 1,nvar
      ivar =  ndim + i
      if (cellcentered(i)) then
         ivarst = ' '
         if      (ivar.lt.10) then
            write(ivarst(2:2),'(i1)') ivar
         else if (ivar.lt.100) then
            write(ivarst(1:2),'(i2)') ivar
         else
            write(*,*) 'TECPLOT_START_FE_ZONE: Out of range ivar =',ivar
            stop
         end if
         varlocation = trim(varlocation)//' ['//trim(ivarst)//'] = CELLCENTERED'
      end if
   end do
   varlocation = trim(varlocation)//')'
else 
   varlocation = ' '
end if

if (nelem.gt.0) then
   ! start finite element zone
   write(iddat,*) ' ZONE N = ',nnod,' E = ',nelem, &
                  ' ZONETYPE = '//trim(zonetype), &
                  ' DATAPACKING = '//trim(datapackingst), &
                  trim(varlocation)
else if (nelem.eq.0) then
   ! start point zone (ordered)
   write(iddat,*) ' ZONE I = ',nnod, &
                  ' DATAPACKING = '//trim(datapackingst), &
                  trim(varlocation)
else
   write(*,*) 'TECPLOT_START_FE_ZONE: Unknown number of elements:',nelem
   stop
end if

                   
return
end subroutine tecplot_start_fe_zone

!*******************************************************************************
subroutine tecplot_start_ordered_zone(iddat,npointx,npointy,npointz,datapacking)
!*******************************************************************************
! Subroutine for starting a new ordered zone in TECPLOT data file
!*******************************************************************************
implicit none

! disk unit
integer, intent(in) :: iddat
! number of points in x,y,z
integer, intent(in) :: npointx, npointy, npointz
! type of data packing
!  0 - block
!  1 - point
integer, intent(in) :: datapacking

! local variables
character(15)  :: zonetype
character(5)   :: datapackingst

! Type of elements
zonetype = 'ORDERED'

! Type of packing data
if      (datapacking.eq.0) then
   datapackingst = 'BLOCK'
else if (datapacking.eq.1) then
   datapackingst = 'POINT'
else
   write(*,*) 'TECPLOT_START_ORDERED_ZONE: Strange datapacking, datapacking =',datapacking
   stop
end if

if (npointx .gt. 0 ) then
   write(iddat,'(a,i7.7,a)')       ' ZONE'
   write(iddat,'(a,i7.7,a)')       '  I=',npointx,','
   if (npointy .gt. 0 ) then
      write(iddat,'(a,i7.7,a)')    '  J=',npointy,','
      if (npointz .gt. 0 ) then
         write(iddat,'(a,i7.7,a)') '  K=',npointz,','
      end if
   end if

   write(iddat,*) ' ZONETYPE = '//trim(zonetype), &
                  ' DATAPACKING = '//trim(datapackingst)
end if
                   
return
end subroutine tecplot_start_ordered_zone

!**********************************************************
subroutine tecplot_export_block_variable_sp(iddat,var,lvar)
!**********************************************************
! Subroutine for exporting a variable in SINGLE PRECISION block format to the TECPLOT data file
!**********************************************************
implicit none

! disk unit
integer, intent(in) :: iddat
! array with variable
integer, intent(in)  :: lvar
real*4, intent(in) ::  var(lvar)

write(iddat,'(5e16.6)') var
                   
return
end subroutine tecplot_export_block_variable_sp

!**********************************************************
subroutine tecplot_export_block_variable_dp(iddat,var,lvar)
!**********************************************************
! Subroutine for exporting a variable in DOUBLE PRECISION block format to the TECPLOT data file
!**********************************************************
implicit none

! disk unit
integer, intent(in) :: iddat
! array with variable
integer, intent(in)  :: lvar
double precision, intent(in) ::  var(lvar)

write(iddat,'(5e16.6)') var
                   
return
end subroutine tecplot_export_block_variable_dp

!*******************************************************************
subroutine tecplot_export_point_data_dp(iddat,array,larray1,larray2)
!*******************************************************************
! Subroutine for exporting a data table in double precision POINT format 
! to the TECPLOT data file
! each row of array correspond to a line in Tecplot file
!*******************************************************************
implicit none

! disk unit
integer, intent(in) :: iddat
! array with variable
integer, intent(in)  :: larray1, larray2
double precision, intent(in) ::  array(larray1,larray2)

! local
integer :: i

do i = 1,larray1
   write(iddat,'(8e18.9)') array(i,:)
end do
                   
return
end subroutine tecplot_export_point_data_dp

!***********************************************************************************
subroutine tecplot_connectivity_table(iddat,ndim,nelem,inet,linet,nnet,lnnet,shells)
!***********************************************************************************
! Subroutine for exporting a the finite element conectivity table
!***********************************************************************************
implicit none

! disk unit
integer, intent(in) :: iddat
! space dimensions
integer, intent(in) :: ndim
! number of elements
integer, intent(in) :: nelem
! array with mesh description
integer, intent(in) :: linet,       lnnet
integer, intent(in) ::  inet(linet), nnet(lnnet)
! are these elements shells?
logical, optional,intent(in) :: shells

! local variables
integer :: ie, nne, i, indinet
logical :: shell_elements = .false.

if (present(shells)) then
   shell_elements = shells
else
   ! default for using shell elements is false
   shell_elements = .false.
end if

! connectivity table
      if (ndim.eq.3) then
      ! 3D
         indinet = 0
         do ie = 1,nelem
            nne = nnet(ie)
            if      (nne.eq.20) then
!     quadratic cubes
               write(iddat,7002) (inet(indinet+i), i = 1,8)
            else if (nne.eq.15) then
!     quadratic prisms
               write(iddat,7002) (inet(indinet+i), i = 1,3), inet(indinet+3), &
                                 (inet(indinet+i), i = 4,6), inet(indinet+6)
            else if (nne.eq.10) then
!     quadratic tetra
               write(iddat,7002) (inet(indinet+i), i = 1,3), inet(indinet+3), &
                                 (inet(indinet+4), i = 1,4)
            else if (nne.eq.8) then
               if (shell_elements) then
!     quadratic semiloofs quadrilaterals
                  write(iddat,7002) (inet(indinet+i), i = 1,4),&
                                    (inet(indinet+i), i = 1,4)
               else
!     linear cubes
                  write(iddat,7002) (inet(indinet+i), i = 1,8)
               end if
            else if (nne.eq.6) then
               if (shell_elements) then
!     quadratic semiloofs triangles
                  write(iddat,7002) (inet(indinet+i), i = 1,3),inet(indinet+3),&
                                    (inet(indinet+i), i = 1,3),inet(indinet+3)
               else
!     linear prisms
                  write(iddat,7002) (inet(indinet+i), i = 1,3),inet(indinet+3),&
                                    (inet(indinet+i), i = 4,6),inet(indinet+6)
               end if
            else if (nne.eq.4) then
!     linear tetra
               write(iddat,7002) (inet(indinet+i), i = 1,3), inet(indinet+3), (inet(indinet+4),i = 1,4)
            else 
               write(*,*) 'TECPLOT_CONNECTIVITY_TABLE: Strange number of nodes for 3D finite element, nne = ', nne
               stop
            end if
            indinet = indinet + nne
         end do
      else if (ndim.eq.2) then
      ! 2D
         indinet = 0
         do ie = 1,nelem
            nne = nnet(ie)
            if      (nne.eq.8) then
               write(iddat,7002) (inet(indinet+i), i = 1,4)
            else if (nne.eq.6) then
               write(iddat,7002) (inet(indinet+i), i = 1,3), inet(indinet+3)
            else if (nne.eq.4) then
               write(iddat,7002) (inet(indinet+i), i = 1,4)
            else
               write(*,*) 'TECPLOT_CONNECTIVITY_TABLE: Strange number of nodes for 2D finite element, nne = ', nne
               stop
            end if
            indinet = indinet + nne
         end do
      else 
         write(*,*) 'TECPLOT_CONNECTIVITY_TABLE: Strange number of dimensions, ndim =',ndim
         stop
      end if

 7002 format(8i15)

return
end subroutine tecplot_connectivity_table

!*********************************************************************************************
subroutine tecplot_get_point_data_header(iddat,nvar,ordered,nx,ny,nz,orderxyz,title,varstring)
!*********************************************************************************************
! Subroutine for reading Tecplot header in datapacking=point format
!*********************************************************************************************
use module_parsing
implicit none

! disk unit
integer, intent(in) :: iddat
! number of variables in the file
integer, intent(out) :: nvar
! is the set ordered
logical, intent(out) :: ordered
! number of points in x, y and z
integer, intent(out) :: nx, ny, nz
! coordinates x, y and z in file
integer, intent(out) :: orderxyz(3)
! string of variable names 
character(*),intent(out) :: title
! string of variable names 
character(*),intent(out) :: varstring

! local vars
integer :: npoint(3)


! initialize counter of lines
fileline = 0

title = ' '
do
   call rdline(iddat)
!   print *,'line',trim(line), 'kstring',kstring
   call getstring

!   print *,'string:',trim(string)
   if (trim(string).eq.'TITLE') then 
      if (kstring.eq.0) then
         call rdline(iddat)
      end if
      call getstring
      if (kstring.eq.0) then
         call rdline(iddat)
      end if
      call getstring
      title = string(2:len(trim(string))-1)
   end if

!   print *,'string:',trim(string)
   if (trim(string).eq.'VARIABLES') then 
      exit
   end if
end do

call getstring
if (kstring.eq.0) then
   call rdline(iddat)
end if

nvar = 0
orderxyz = 0
varstring = ' '
do
   call getstring
   if (trim(string).eq.'ZONE') then
      exit
   end if
!   print *,'string:',trim(string)

   nvar = nvar + 1
   if      (string(2:2) .eq. 'X') then
      orderxyz(1) = nvar
   else if (string(2:2) .eq. 'Y') then
      orderxyz(2) = nvar
   else if (string(2:2) .eq. 'Z') then
      orderxyz(3) = nvar
   else
! create the string of variables
      !if (varstring .ne. ' ') then
      !   varstring = trim(varstring)//','
      !end if
      varstring = trim(varstring)//' '//trim(string)
   end if

   if (kstring.eq.0) then
      call rdline(iddat)
   end if
end do

do 
   if (kstring.eq.0) then
      call rdline(iddat)
   end if
   call getstring

   if (string(1:2).eq.'I=') then
      read(string(3:),*) npoint(1)
      if (kstring.eq.0) then
         call rdline(iddat)
      end if
      call getstring
      read(string(3:),*) npoint(2)
      if (kstring.eq.0) then
         call rdline(iddat)
      end if
      call getstring
      read(string(3:),*) npoint(3)
      if (kstring.eq.0) then
         call rdline(iddat)
      end if
      call getstring
      if (string.eq.'ZONETYPE=Ordered') then
         ordered = .true.
      else
         ordered = .false.
      end if
      exit
   end if
end do

nx = npoint(orderxyz(1))
ny = npoint(orderxyz(2))
nz = npoint(orderxyz(3))

! pad the rest of the file up to values
do 
   call rdline(iddat)
   if (line(1:3).eq.'DT=') then
      do while (kstring.ne.0)
         call getstring
      end do
      exit
   end if
end do

! check orderxyz
if (any(orderxyz.eq.0)) then
   write(*,*) 'TECPLOT_GET_POINT_DATA_HEADER: Missing coordinates in file.'
   stop
end if

return
end subroutine tecplot_get_point_data_header

end module module_tecplot

