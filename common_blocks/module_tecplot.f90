module module_tecplot
! Module for export to Tecplot ASCII data file
implicit none
integer,parameter,private :: kr = kind(1.D0)

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

!*********************************************************************************
subroutine tecplot_start_zone(iddat,ndim,nnod,nelem,nvar,cellcentered,datapacking)
!*********************************************************************************
! Subroutine for starting a new zone in TECPLOT data file
!*********************************************************************************
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
   write(*,*) 'TECPLOT_START_ZONE: Strange number of dimensions, ndim =',ndim
   stop
end if


! Type of packing data
if      (datapacking.eq.0) then
   datapackingst = 'BLOCK'
else if (datapacking.eq.1) then
   datapackingst = 'POINT'
else
   write(*,*) 'TECPLOT_START_ZONE: Strange datapacking, datapacking =',datapacking
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
            write(*,*) 'TECPLOT_START_ZONE: Out of range ivar =',ivar
            stop
         end if
         varlocation = trim(varlocation)//'['//trim(ivarst)//'] = CELLCENTERED'
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
   write(*,*) 'TECPLOT_START_ZONE: Unknown number of elements:',nelem
   stop
end if

                   
return
end subroutine tecplot_start_zone

!*******************************************************
subroutine tecplot_export_block_variable(iddat,var,lvar)
!*******************************************************
! Subroutine for exporting a variable in block format to the TECPLOT data file
!*******************************************************
implicit none

! disk unit
integer, intent(in) :: iddat
! array with variable
integer, intent(in)  :: lvar
real(kr), intent(in) ::  var(lvar)

write(iddat,'(5e15.6)') var
                   
return
end subroutine tecplot_export_block_variable

!****************************************************************************
subroutine tecplot_connectivity_table(iddat,ndim,nelem,inet,linet,nnet,lnnet)
!****************************************************************************
! Subroutine for exporting a the finite element conectivity table
!****************************************************************************
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

! local variables
integer :: ie, nne, i, indinet

! connectivity table
      if (ndim.eq.3) then
      ! 3D
         indinet = 0
         do ie = 1,nelem
            nne = nnet(ie)
            if      (nne.eq.20) then
               write(iddat,7002) (inet(indinet+i), i = 1,8)
            else if (nne.eq.15) then
               write(iddat,7002) (inet(indinet+i), i = 1,3), inet(indinet+3), &
                                 (inet(indinet+i), i = 4,6), inet(indinet+6)
            else if (nne.eq.10) then
               write(iddat,7002) (inet(indinet+i), i = 1,3), inet(indinet+3), &
                                 (inet(indinet+4), i = 1,4)
!     linear cubes
            else if (nne.eq.8) then
               write(iddat,7002) (inet(indinet+i), i = 1,8)
!     linear tetra
            else if (nne.eq.4) then
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

!********************************************************************
subroutine tecplot_get_point_data_header(iddat,nvar,ordered,nx,ny,nz)
!********************************************************************
! Subroutine for exporting a the finite element conectivity table
!********************************************************************
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

! initialize counter of lines
fileline = 0

do
   call rdline(iddat)
!   print *,'line',trim(line), 'kstring',kstring
   call getstring

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
do
   call getstring
   if (trim(string).eq.'ZONE') then
      exit
   end if
!   print *,'string:',trim(string)

   nvar = nvar + 1

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
      read(string(3:),*) nx
      if (kstring.eq.0) then
         call rdline(iddat)
      end if
      call getstring
      read(string(3:),*) ny
      if (kstring.eq.0) then
         call rdline(iddat)
      end if
      call getstring
      read(string(3:),*) nz
      exit
   end if
end do

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

return
end subroutine tecplot_get_point_data_header

end module module_tecplot

