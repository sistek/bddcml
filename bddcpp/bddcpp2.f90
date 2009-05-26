!***********************************************************************
program bddcpp
!***********************************************************************
! BDDC preprocessor
! Program for preparing subdomain data from global PMD data
!
! *.GMIS - ASCII file with mesh data
! *.FVS  - ASCII file with fixed variables
! *.ELM  - binary file with element matrices
! *.RHS  - binary file with global right hand side vector
!
! plus requires extra files
!
! *.PAR  - ASCII file with basic parameters of problem
! either of the two: (if ES does not exist, it is created from ES0)
! *.ES0  - ASCII file with list of subdomains's numbers of elements,
!          numbered from 0 to nsub-1
! *.ES   - ASCII file with list of subd. numbers of elements, from 1 
! *.CN   - ASCII file with list of "corner nodes"
!
! programmed by Jakub Sistek            Denver                  5.3.2009 
!***********************************************************************

implicit none
      
logical :: debug = .false.

integer,parameter:: kr = kind(1.D0)
integer,parameter:: idpar = 1, idbase = 100, idgmi = 2, ides = 4, idfvs = 3, &
                    idelm = 10, idrhs = 11, iddat = 12, idcn = 13, ides0 = 14, &
                    idrhss = 15, idfvss = 16, idsols = 17, idsol = 18, idglb = 19, &
                    idgmist = 20, idgmis = 21, idgraph = 22, idint = 23, idpair = 24

integer,parameter:: lname1x = 8, lnamex = 100, lfnamex = 20

integer:: lname1
integer:: ndim, nsub, nelem, ndof, nnod, nnodc

integer ::           linet,   lnnet,   lnndf,   lkinet,   liets,   lkmynodes,   lkinodes,   lkneibnodes,   lkdof,   lkadjsnodes
integer,allocatable:: inet(:), nnet(:), nndf(:), kinet(:), iets(:), kmynodes(:), kinodes(:), kneibnodes(:), kdof(:), kadjsnodes(:)
integer ::           linodc
integer,allocatable:: inodc(:)
 
integer ::           lifix
integer,allocatable:: ifix(:)

integer ::            lelm,   lx,   ly,   lz,   lrhs,   lfixv,   lsol   
real(kr),allocatable:: elm(:), x(:), y(:), z(:), rhs(:), fixv(:), sol(:)
integer ::             lxyz1, lxyz2
real(kr),allocatable::  xyz(:,:)

character(lname1x):: name1
character(lnamex)::  name
character(lfnamex)::  fname

type glob_type
   integer ::             itype
   logical ::             selected
   integer ::             nnod
   integer,allocatable :: nodes(:)
   integer ::             nsub
   integer,allocatable :: subdomains(:)
end type glob_type
integer::                      lglobs
type(glob_type), allocatable :: globs(:)

integer :: nc, ios


! Initial screen
      write(*,'(a)') ' _____  _____  _____   ____  _____  _____  '
      write(*,'(a)') '|  _  \|  _  \|  _  \ / __ \|  _  \|  _  \ '
      write(*,'(a)') '| |_|  | | \  | | \  | /  \_| |_|  | |_|  |'
      write(*,'(a)') '|  ___/| |  | | |  | | |    |  ___/|  ___/ '
      write(*,'(a)') '|  _  \| |  | | |  | | |   _| |    | |     '
      write(*,'(a)') '| |_|  | |_/  | |_/  | \__/ | |    | |     '
      write(*,'(a)') '|_____/|_____/|_____/ \____/|_|    |_|     '
      write(*,'(a)') 'BDDCPP - Pre/Post-processor for BDDC method'
      write(*,'(a)') '==========================================='

! Name of the problem
   10 write(*,'(a,$)') 'Name of the problem: '
      read(*,*) name1
      if(name1.eq.' ') goto 10
      lname1   = index(name1,' ') - 1
      if(lname1.eq.-1) then
        lname1 = lname1x
      end if

! Open files
! PAR - basic properties of the problem    
      name = name1(1:lname1)//'.PAR'
      open (unit=idpar,file=name,status='old',form='formatted')

! Read basic parameters for allocation
      read(idpar,*) ndim, nsub, nelem, ndof, nnod, nnodc, linet

!*****************
! Select operation
 20   print '("===================")'
      print '("And what we do now?")'
      print '("1 - generate element graph for the mesh")'
      print '("2 - divide mesh into subdomains by METIS")'
      print '("3 - create subdomain files with element matrices")'
      print '("4 - select globs")'
      print '("5 - generate space W_tilde")'
      print '("6 - export division of subdomains for TECPLOT")'
      print '("7 - assembly global solution from subdomain files")'
      print '("8 - export solution to TECPLOT")'
      print '("9 - select corner nodes (random selection from interface)")'
      print '("10- create file ES from ES0")'
      print '("11- create subdomain mesh files SMD")'
      print '("0 - quit")'
      print '("your choice: ",$)'; read *,nc
      print '("===================")'
      select case (nc)
         case (1)
            call graphpmd
            goto 20
         case (2)
            call meshdivide
            goto 20
         case (3)
            call createem(name1(1:lname1))
            goto 20
         case (4)
            call getglobs
            goto 20
         case (5)
            call createtilde
            goto 20
         case (6)
            call exptecplotonesubzones
            goto 20
         case (7)
            call solassembly
            goto 20
         case (8)
            call exptecplotsolution
            goto 20
         case (9)
            call getcorners
            goto 20
         case (10)
            call es0toes
            goto 20
         case (11)
            call create_sub_files(name1(1:lname1))
            goto 20
         case (0)
            write(*,*) 'O.K.'
            stop
         case default
            goto 20
      end select

contains

!***********************************************************************
subroutine exptecplotonesubzones
!***********************************************************************
!     Subroutine for export to TECPLOT file of 3D data, 
!     one zone per subdomain
!***********************************************************************
use module_tecplot
use module_errors
implicit none
      
! Local variables

! number of variables
integer,parameter :: nvar = 1
! is this variable cellcentered?
logical :: cellcentered(nvar) = (/.true./) 
! type of data packing
! type of data packing
!  0 - block
!  1 - point
integer :: datapacking = 0

integer ::             linets,   lnnets,   lkinet
integer, allocatable :: inets(:), nnets(:), kinet(:)
integer::              linodc
integer,allocatable::   inodc(:)
integer ::             lmapping
integer, allocatable :: mapping(:), subind(:)
integer ::             lkmynodes
integer, allocatable :: kmynodes(:)
integer ::              lsubaux
real(kr), allocatable :: subaux(:)
integer :: nne, indnnets, indinet, indinets, i, imapping, inod, &
           ie, isub, inods
integer :: nnods, nelems
character(1) :: yn

logical :: exportcorners

! Import of basic geometry
! GMIS - basic mesh data - structure:
!  * INET(LINET) * NNET(LNNET) * NNDF(LNNDF) * XYF(LXYF) *
      name = name1(1:lname1)//'.GMIS'
      open (unit=idgmi,file=name,status='old',form='formatted')
      rewind idgmi
      linet = linet
      lnnet = nelem
      lnndf = nnod
      lx    = nnod
      ly    = nnod
      lz    = nnod
      allocate(inet(linet),nnet(lnnet),nndf(lnndf))
      if      (ndim.eq.3) then
         allocate(x(lx),y(ly),z(lz))
      else if (ndim.eq.2) then
         allocate(x(lx),y(ly))
      else
         write(*,*) 'EXPTECPLOTSUBZONES: Strange number of space dimensions, ndim = ',ndim
      end if
! read fields INET, NNET, NNDF from file IDGMI
      read(idgmi,*) inet
      read(idgmi,*) nnet
      read(idgmi,*) nndf
      if      (ndim.eq.3) then
         read(idgmi,*) x,y,z
      else if (ndim.eq.2) then
         read(idgmi,*) x,y
      end if
      close(idgmi)

! ES - list of global element numbers in subdomains - structure:
      name = name1(1:lname1)//'.ES'
      open (unit=ides,file=name,status='old',form='formatted')
      liets = nelem
      allocate(iets(liets))
      read(ides,*) iets
      close(ides)

  55  write(*,'(a,$)') 'Do you want to export corners (y/n)? '
      read(*,*) yn
      select case (yn)
         case ('Y','y')
            exportcorners = .true.
         case ('N','n')
            exportcorners = .false.
         case default
            print *, 'Unknown answer, try again ...'
            goto 55
      end select

! CN - list of corners
      if (exportcorners) then
         name = name1(1:lname1)//'.CN'
         open (unit=idcn,file=name,status='old',form='formatted', iostat=ios)
         if(ios.ne.0) then
            call error('exptecplotonesubzones','File '//trim(name)//' does not exist, corners will not be exported.')
         else
            ! read number of corners
            read(idcn,*) nnodc
            ! allocate array for indices of corners
            linodc = nnodc
            allocate(inodc(linodc))
            ! read indices of corners in global node numbering in hat
            read(idcn,*) inodc
            close(idcn)
         end if
      end if

! DAT - data for TECPLOT:
      name = name1(1:lname1)//'_sub.dat'
      open (unit=iddat,file=name,status='replace',form='formatted')

! write header
      call tecplot_header(iddat,ndim,name1,lname1,'"SUBDOMAIN"')

! Creation of field KINET(NELEM) with addresses before first global node of element IE in field inet
      lkinet = nelem
      allocate(kinet(lkinet))
      kinet(1) = 0
      do ie = 2,nelem
         kinet(ie) = kinet(ie-1) + nnet(ie-1)
      end do

      lkmynodes = nnod
      allocate(kmynodes(lkmynodes))
! loop over subdomains
      do isub = 1,nsub
! prepare subdomain data
         nelems = count(iets == isub)
! Begin loop over subdomains
! Creation of field kmynodes(nnod) - local usage of global nodes
! if global node is in my subdomain, assigned 1, else remains 0
         call pp_mark_sub_nodes(isub,nelem,iets,liets,inet,linet,nnet,lnnet,&
                                kinet,lkinet,kmynodes,lkmynodes)
         nnods = count(kmynodes == 1)
         lmapping = nnods
         allocate(mapping(lmapping))
         imapping = 0
         do inod = 1,nnod
            if (kmynodes(inod) == 1) then
               imapping = imapping + 1
               mapping(imapping) = inod
            end if
         end do

         linets = 0
         do ie = 1,nelem
            if (iets(ie) == isub) then
               linets = linets + nnet(ie)
            end if
         end do
         lnnets = nelems
         allocate (inets(linets),nnets(lnnets))
         indinet  = 0
         indinets = 0
         indnnets = 0
         do ie = 1,nelem
            nne = nnet(ie)
            if (iets(ie) == isub) then
               indnnets = indnnets + 1
               nnets(indnnets) = nne
               inets(indinets+1:indinets+nne) = inet(indinet+1:indinet+nne)
               indinets = indinets + nne
            end if
            indinet = indinet + nne
         end do

         do inods = 1,nnods
            kmynodes(mapping(inods)) = inods
         end do
         do i = 1,linets
            inets(i) = kmynodes(inets(i))
         end do


! start zone
         call tecplot_start_zone(iddat,ndim,nnods,nelems,nvar,cellcentered,datapacking)
! export coordinates
         call tecplot_export_block_variable(iddat,x(mapping),nnods)
         call tecplot_export_block_variable(iddat,y(mapping),nnods)
         if (ndim.eq.3) then
            call tecplot_export_block_variable(iddat,z(mapping),nnods)
         end if
! export variables
         allocate(subind(nelems))
         subind = isub
         call tecplot_export_block_variable(iddat,dble(subind),nelems)
         deallocate(subind)

         call tecplot_connectivity_table(iddat,ndim,nelems,inets,linets,nnets,lnnets)

         deallocate(mapping)
         deallocate(inets,nnets)

      end do
      deallocate(kinet,kmynodes)

! export corners as a next zone
      if (exportcorners) then
         datapacking  = 0
         cellcentered = .false.
         call tecplot_start_zone(iddat,ndim,nnodc,0,nvar,cellcentered,datapacking)
         call tecplot_export_block_variable(iddat,x(inodc),nnodc)
         call tecplot_export_block_variable(iddat,y(inodc),nnodc)
         if (ndim.eq.3) then
            call tecplot_export_block_variable(iddat,z(inodc),nnodc)
         end if
         lsubaux = nnodc
         allocate(subaux(lsubaux))
         subaux = nsub + 1
         call tecplot_export_block_variable(iddat,subaux,nnodc)
         deallocate(subaux)
      end if

      print *, 'Division exported into TECPLOT file ', trim(name)
      close(iddat)

! Clear memory
      deallocate(iets)
      deallocate(inet,nnet,nndf)
      if      (ndim.eq.3) then
         deallocate(x,y,z)
      else if (ndim.eq.2) then
         deallocate(x,y)
      end if
      if (allocated(inodc)) then
         deallocate(inodc)
      end if

end subroutine exptecplotonesubzones

!***********************************************************************
subroutine exptecplotsolution
!***********************************************************************
!     Subroutine for export to TECPLOT file of solution
!***********************************************************************
use module_tecplot
use module_errors
implicit none
      
! Local variables
integer:: indinet, nne, ie, indn, ine, inod, node1, node2, ncne, ndofn, &
          idofn, indsol, inddim

integer ::             lvx,   lvy,   lvz,   lp
real(kr),allocatable :: vx(:), vy(:), vz(:), p(:)

integer ::             ldisplacement1, ldisplacement2
real(kr),allocatable :: displacement(:,:)

! number of variables
integer :: nvar
! is this variable cellcentered?
logical,allocatable :: cellcentered(:)
! type of data packing
! type of data packing
!  0 - block
!  1 - point
integer :: datapacking = 0
! what kind of problem
integer :: problemid
logical :: flow = .false.
logical :: elasticity = .false.

! Select problem
20    write(*,'(a)') "Kind of problem ?"
      write(*,'(a)') "1 - elasticity"
      write(*,'(a)') "2 - incompressible flow"
      write(*,'(a,$)') "problem ID: "
      read (*,*) problemid
      select case (problemid)
         case (1)
            elasticity = .true.
         case (2)
            flow = .true.
         case default
            write(*,'(a)') "Unknown problem ID, try again ..."
            goto 20
      end select

! Import of basic geometry
! GMIS - basic mesh data - structure:
!  * INET(LINET) * NNET(LNNET) * NNDF(LNNDF) * XYF(LXYF) *
      name = name1(1:lname1)//'.GMIS'
      open (unit=idgmi,file=name,status='old',form='formatted')
      rewind idgmi
      linet = linet
      lnnet = nelem
      lnndf = nnod
      lxyz1 = nnod
      lxyz2 = ndim
      allocate(inet(linet),nnet(lnnet),nndf(lnndf))
      allocate(xyz(lxyz1,lxyz2))
! read fields INET, NNET, NNDF from file IDGMI
      read(idgmi,*) inet
      read(idgmi,*) nnet
      read(idgmi,*) nndf
      read(idgmi,*) xyz
      close(idgmi)

! ES - list of global element numbers in subdomains - structure:
      name = name1(1:lname1)//'.SOL'
      open (unit=idsol,file=name,status='old',form='unformatted')
      lsol = ndof
      allocate(sol(lsol))
      read(idsol) sol
      close(idsol)

! Creation of field KDOF(NNOD) with addresses before first global
! dof of node
      lkdof = nnod
      allocate(kdof(lkdof))
      kdof(1) = 0
      do inod = 2,nnod
         kdof(inod) = kdof(inod-1) + nndf(inod-1)
      end do

      if (flow) then
! Convert solution to velocities
         lvx = nnod
         lvy = nnod
         lp  = nnod
         allocate(vx(lvx),vy(lvy),p(lp))

! Interpolate pressure to midsides
         indinet = 0
         do ie = 1,nelem
            nne = nnet(ie)
            ncne = nne/2
            do ine = 1,nne
               indn = inet(indinet + ine)
               ndofn = nndf(indn)
               vx(indn) = sol(kdof(indn)+1)
               vy(indn) = sol(kdof(indn)+2)
               if (ndofn.eq.ndim+1) then
                  p(indn) = sol(kdof(indn)+3)
               else
                  node1 = inet(indinet+ine-ncne)
                  if (ine.ne.nne) then
                     node2 = inet(indinet+ine-ncne+1)
                  else
                     node2 = inet(indinet+1)
                  end if
                  p(indn) = (p(node1) + p(node2))/2
               end if
            end do
            indinet = indinet + nne
         end do
      endif

      if (elasticity) then
         ldisplacement1 = nnod
         ldisplacement2 = ndim
         allocate (displacement(ldisplacement1,ldisplacement2))
         displacement = 0

         indsol = 0
         do inod = 1,nnod
            ndofn = nndf(inod)
            do idofn = 1,ndofn
               indsol = indsol + 1
               displacement(inod,idofn) = sol(indsol)
            end do
         end do
         ! check length
         if (indsol.ne.lsol) then
            call error('EXPTECPLOTSOLUTION','Length of vector of solution does not match.')
         end if
      end if

! Visual check
!      do inod = 1,nnod
!         write(*,'(i5,5f15.10)') inod, vx(inod), vy(inod), p(inod)
!      end do

! DAT - data for TECPLOT:
      name = name1(1:lname1)//'_sol.dat'
      open (unit=iddat,file=name,status='replace',form='formatted')

! write header
      if (flow) then
         if (ndim.eq.2) then
            nvar = 3
            allocate(cellcentered(nvar))
            cellcentered = .false.
            call tecplot_header(iddat,ndim,name1,lname1,'"VX", "VY", "P"')
         else if (ndim.eq.3) then
            nvar = 4
            allocate(cellcentered(nvar))
            cellcentered = .false.
            call tecplot_header(iddat,ndim,name1,lname1,'"VX", "VY", "VZ", "P"')
         end if
      else if (elasticity) then
         if (ndim.eq.2) then
            nvar = 3
            allocate(cellcentered(nvar))
            cellcentered = .false.
            call tecplot_header(iddat,ndim,name1,lname1,'"DX", "DY", "DMAGNITUDE", "XDEF", "YDEF" ')
         else if (ndim.eq.3) then
            call tecplot_header(iddat,ndim,name1,lname1,'"DX", "DY", "DZ", "DMAGNITUDE", "XDEF", "YDEF", "ZDEF" ')
            nvar = 4
            allocate(cellcentered(nvar))
            cellcentered = .false.
         end if
      else
         call warning('EXPTECPLOTSOLUTION','No file will be exported.')
         return
      end if

! start zone
      call tecplot_start_zone(iddat,ndim,nnod,nelem,nvar,cellcentered,datapacking)

! export coordinates
      do inddim = 1,ndim
         call tecplot_export_block_variable(iddat,xyz(:,inddim),lxyz1)
      end do

! export variables
      if (flow) then
         call tecplot_export_block_variable(iddat,vx,lvx)
         call tecplot_export_block_variable(iddat,vy,lvy)
         if (ndim.eq.3) then
            call tecplot_export_block_variable(iddat,vz,lvz)
         end if
         call tecplot_export_block_variable(iddat,p, lp)
      else if (elasticity) then
         do inddim = 1,ndim
            call tecplot_export_block_variable(iddat,displacement(:,inddim),ldisplacement1)
         end do
         ! export displacement magnitude
         call tecplot_export_block_variable(iddat,sqrt(sum(displacement(:,:)**2,DIM=2)),ldisplacement1)
         do inddim = 1,ndim
            call tecplot_export_block_variable(iddat,xyz(:,inddim)+displacement(:,inddim),lxyz1)
         end do
      end if
   
      call tecplot_connectivity_table(iddat,ndim,nelem,inet,linet,nnet,lnnet)

      print *, 'Solution exported into TECPLOT file ', trim(name)
      close(iddat)

! Clear memory
      deallocate(cellcentered)
      deallocate(kdof)
      if (flow) then
         deallocate(vx,vy,p)
      end if
      if (elasticity) then
         deallocate(displacement)
      end if
      deallocate(sol)
      deallocate(inet,nnet,nndf)
      deallocate(xyz)

end subroutine exptecplotsolution

!***********************************************************************
subroutine createem(problemname)
!***********************************************************************
!     Subroutine for creation files with element matrices for subdomains
!***********************************************************************
use module_utils
implicit none

character(*),intent(in) :: problemname

! Local variables
integer :: idelms
integer :: i, inod, indn, idofn, ine, isub, nne, lelmx, nevab, ndofn, nevax, point, ie, &
           indsub

character(100) :: filename

      write(*,*) 'Generating subdomain files with element matrices...'

! Import of basic geometry
! GMIS - basic mesh data - structure:
!  * INET(LINET) * NNET(LNNET) * NNDF(LNNDF) * XYF(LXYF) *
      name = name1(1:lname1)//'.GMIS'
      open (unit=idgmi,file=name,status='old',form='formatted')
      rewind idgmi
      linet = linet
      lnnet = nelem
      lnndf = nnod
      allocate(inet(linet),nnet(lnnet),nndf(lnndf))
! read fields INET, NNET, NNDF from file IDGMI
      read(idgmi,*) inet
      read(idgmi,*) nnet
      read(idgmi,*) nndf
      close(idgmi)

! ES - list of global element numbers in subdomains - structure:
      name = name1(1:lname1)//'.ES'
      open (unit=ides,file=name,status='old',form='formatted')
      liets = nelem
      allocate(iets(liets))
      read(ides,*) iets
      close(ides)

! ELM - element stiffness matrices - structure:
      name = name1(1:lname1)//'.ELM'
      open(unit=idelm,file=name,status='old',form='unformatted')
      rewind idelm
        
! Creation of field KINET(NELEM) with addresses before first global node of element IE in field inet
      lkinet = nelem
      allocate(kinet(lkinet))
      kinet(1) = 0
      do ie = 2,nelem
         kinet(ie) = kinet(ie-1) + nnet(ie-1)
      end do

! Find nevax
      nevax = 0
      do ie = 1, nelem
         nevab = 0
         nne = nnet(ie)
         do ine = 1,nne
            indn = inet(kinet(ie) + ine)
            nevab = nevab + nndf(indn)
         end do
         if (nevab.gt.nevax) then
            nevax = nevab
         end if
      end do
      
      lelmx = (nevax+1)*nevax/2
      
      allocate(elm(lelmx))

! If number of subdomains is not extensive,
      if (nsub.le.1024) then
! open a file for each subdomain and distribute element like cards in one pass through element matrices.
         do isub = 1,nsub
            call getfname(problemname,isub,'ELM',filename)
            idelms = idbase + isub
            open(unit=idelms, file=filename, status='replace', form='unformatted')
         end do
!        Write the element matrix to a proper file
         do ie = 1,nelem
            isub = iets(ie)
            read(idelm)   lelm,(elm(i), i = 1,lelm)
            idelms = idbase + isub
            write(idelms) lelm,(elm(i), i = 1,lelm)
         end do
! Close all element files
         do isub = 1,nsub
            idelms = idbase + isub
            close(idelms)
         end do
      else
!     Create element files one by one in NSUB passes through element matrices -
!     it is slow but robust
         do isub = 1,nsub
            call getfname(problemname,isub,'ELM',filename)
            idelms = idbase + 1
            open(unit=idelms, file=filename, status='replace', form='unformatted')
            rewind idelm
            do ie = 1,nelem
               indsub = iets(ie)
               ! Is the element in current subdomain?
               if (indsub.eq.isub) then
                  read(idelm)   lelm,(elm(i), i = 1,lelm)
                  write(idelms) lelm,(elm(i), i = 1,lelm)
               else
                  ! move one record
                  read(idelm)  
               end if
            end do
            close(idelms)
         end do
      end if

      write(*,*) '...done'

      deallocate(elm)
      close(idelm)

! Creation of field KDOF(NNOD) with addresses before first global
! dof of node
      lkdof = nnod
      allocate(kdof(lkdof))
      kdof(1) = 0
      do inod = 2,nnod
         kdof(inod) = kdof(inod-1) + nndf(inod-1)
      end do

! RHS - global right hand side - structure:
      name = name1(1:lname1)//'.RHS'
      open (unit=idrhs,file=name,status='old',form='unformatted')
      rewind idrhs
      lrhs  = ndof
      allocate(rhs(lrhs))
! Read RHS
      read(idrhs) rhs
      close(idrhs)

! Distribute RHS and FIXV
      write(*,*) 'Generating subdomain files with right hand sides...'
      lkmynodes = nnod
      allocate(kmynodes(lkmynodes))
      do isub = 1,nsub
         call bddc_getfname(name1,lname1,isub,'RHSS',fname)
         open(unit=idrhss, file=fname, status='replace', form='unformatted')
! Begin loop over subdomains
! Creation of field kmynodes(nnod) - local usage of global nodes
! if global node is in my subdomain, assigned 1, else remains 0
         call pp_mark_sub_nodes(isub,nelem,iets,liets,inet,linet,nnet,lnnet,&
                                kinet,lkinet,kmynodes,lkmynodes)
         do inod = 1,nnod
            if (kmynodes(inod).eq.1) then
               ndofn = nndf(inod)
               point = kdof(inod)
               do idofn = 1,ndofn
                  write(idrhss) rhs(point + idofn)
               end do
            end if
         end do
! Close all RHSS files
         close(idrhss)
      end do
      deallocate(kmynodes)
      deallocate(rhs)
      write(*,*) '...done'
      
! FVS - fixed variables
!  * IFIX(LIFIX) * FIXV(LFIXV) * 
      name = name1(1:lname1)//'.FVS'
      open (unit=idfvs,file=name,status='old',form='formatted')
      rewind idfvs
      lifix = ndof
      lfixv = ndof
      allocate(ifix(lifix),fixv(lfixv))
! Read fixed variables
      read(idfvs,*) ifix
      read(idfvs,*) fixv
      close(idfvs)

! Distribute RHS and FIXV
      write(*,*) 'Generating subdomain files with fixed variables...'
      lkmynodes = nnod
      allocate(kmynodes(lkmynodes))
      do isub = 1,nsub
         call bddc_getfname(name1,lname1,isub,'FVSS',fname)
         open(unit=idfvss, file=fname, status='replace', form='unformatted')
! Begin loop over subdomains
! Creation of field kmynodes(nnod) - local usage of global nodes
! if global node is in my subdomain, assigned 1, else remains 0
         call pp_mark_sub_nodes(isub,nelem,iets,liets,inet,linet,nnet,lnnet,&
                                kinet,lkinet,kmynodes,lkmynodes)
         do inod = 1,nnod
            if (kmynodes(inod).eq.1) then
               ndofn = nndf(inod)
               point = kdof(inod)
               do idofn = 1,ndofn
                  write(idfvss) fixv(point + idofn)
               end do
            end if
         end do
! Close all RHSS files
         close(idfvss)
      end do
      deallocate(kmynodes)
      deallocate(ifix,fixv)
      write(*,*) '...done'

      deallocate(kinet)
      deallocate(kdof)
      deallocate(inet,nnet,nndf)
      deallocate(iets)

end subroutine createem

!***********************************************************************
subroutine getglobs
!***********************************************************************
!     Subroutine for glob recognition and selection.
!     Also creates set of corners such that each pair of subdomains 
!     sharing a face is joint by at least three corners.
!***********************************************************************
use module_errors
implicit none

! Local variables
integer::            lsublist1, lsublist2
integer,allocatable:: sublist(:,:)
integer::            lnodeglob
integer,allocatable:: nodeglob(:)
integer::            ligingn
integer,allocatable:: igingn(:)

integer::            linglb,   lnnglb
integer,allocatable:: inglb(:), nnglb(:)
integer::            lnewvertex
integer,allocatable:: newvertex(:)
integer::            lisingin
integer,allocatable:: isingin(:)
integer::            lplayground
integer,allocatable:: playground(:)
integer::            linodc
integer,allocatable:: inodc(:)
integer::            licheck
integer,allocatable:: icheck(:)
integer::             ldist,   larea
real(kr),allocatable:: dist(:), area(:)
integer::             lxyz1, lxyz2 
real(kr),allocatable:: xyz(:,:)
integer::             lxyzsh1, lxyzsh2 
real(kr),allocatable:: xyzsh(:,:)
integer::             lxyzbase
real(kr),allocatable:: xyzbase(:)
integer::             lsubneib
integer,allocatable::  subneib(:)

integer:: kselection
integer:: inod, jnodi, jnod, isub, jsub, inodi, iglob, iglb, ing, iinglb, iglobnode,&
          i, nnodi, nmax, nng, nglb, nglob, nplaces, ie,&
          indsub, nnewvertex, nnodcs, nnodis, indisn, inodis, inodcs, inodcs1, &
          inodcs2, inodcs3, indi, inewnodes, indinode, ncorrections, &
          indn, nglobn, nnodcpair, nshared, inodsh, ish, nnsub, iisub, &
          indjsub, inddof, ndofn, point, indnod, nremoved_corners_bc,&
          ipair, npair
integer:: indaux(1)
integer:: nface, nedge, nvertex
real(kr):: x1, y1, z1, x2, y2, z2, rnd, xish, yish, zish
integer:: nnodcold, inc

logical:: minimizecorners
character(1) :: yn

type neighbouring_type
   integer ::             nnsub
   integer,allocatable :: list(:)
end type neighbouring_type
integer::                      lneighbourings
type(neighbouring_type), allocatable :: neighbourings(:)

! Import of basic geometry
! GMIS - basic mesh data - structure:
!  * INET(LINET) * NNET(LNNET) * NNDF(LNNDF) * XYF(LXYF) *
      name = name1(1:lname1)//'.GMIS'
      open (unit=idgmi,file=name,status='old',form='formatted')
      rewind idgmi
      linet = linet
      lnnet = nelem
      lnndf = nnod
      lxyz1 = nnod
      lxyz2 = ndim
      allocate(inet(linet),nnet(lnnet),nndf(lnndf))
      allocate(xyz(lxyz1,lxyz2))
! read fields INET, NNET, NNDF from file IDGMI
      read(idgmi,*) inet
      read(idgmi,*) nnet
      read(idgmi,*) nndf
      read(idgmi,*) xyz
      close(idgmi)

! ES - list of global element numbers in subdomains - structure:
      name = name1(1:lname1)//'.ES'
      open (unit=ides,file=name,status='old',form='formatted')
      liets = nelem
      allocate(iets(liets))
      read(ides,*) iets
      close(ides)

! Creation of field KINET(NELEM) with addresses before first global node of element IE in field inet
      lkinet = nelem
      allocate(kinet(lkinet))
      kinet(1) = 0
      do ie = 2,nelem
         kinet(ie) = kinet(ie-1) + nnet(ie-1)
      end do

! Minimize corners?
  55  write(*,'(a,$)') 'Do you want to minimize amount of corners (y/n)? ("no" guarantee better distribution of corners) '
      read(*,*) yn
      select case (yn)
         case ('Y','y')
      print *, ndim, nsub, nelem, ndof, nnod, nnodc, linet
            minimizecorners = .true.
         case ('N','n')
            minimizecorners = .false.
         case default
            print *, 'Unknown answer, try again ...'
            goto 55
      end select

! Recognize interface
      lkinodes  = nnod
      lkmynodes = nnod
      lkneibnodes = nnod
      allocate(kinodes(lkinodes),kmynodes(lkmynodes),kneibnodes(lkneibnodes))
! Begin loop over subdomains
      kinodes = 0
      do isub = 1,nsub
! Creation of field kmynodes(nnod) - local usage of global nodes
! if global node is in my subdomain, assigned 1, else remains 0
         call pp_mark_sub_nodes(isub,nelem,iets,liets,inet,linet,nnet,lnnet,&
                                kinet,lkinet,kmynodes,lkmynodes)
! Mark interface nodes
         where(kmynodes.eq.1) kinodes = kinodes + 1
      end do
! Count interface nodes
      nnodi = count(kinodes.gt.1)
      write(*,'(a,i8)') 'Number of nodes on interface nnodi =',nnodi

! Create field indices of global interface nodes in global numbering
      ligingn = nnodi
      allocate(igingn(ligingn))
      inodi = 0
      do inod = 1,nnod
         if (kinodes(inod).gt.1) then
            inodi = inodi + 1
            igingn(inodi) = inod
         end if
      end do

! What is the maximal number of subdomains sharing one node?
      nmax = maxval(kinodes)
! How many nodes are shared by the maximal number of subdomains?
      nplaces = count(kinodes.eq.nmax)
      write(*,'(a,i3,a,i7,a)') 'Maximum number of touching subdomains nmax =',nmax,' on ',nplaces,' places.'

! Recognize faces and edges
! Create list of subdomains that a node is included in
      lsublist1 = nnodi
      lsublist2 = nmax
      allocate(sublist(lsublist1,lsublist2))
      sublist = 0
      kinodes = 0
! Begin loop over subdomains
      do isub = 1,nsub
! Creation of field kmynodes(nnod) - local usage of global nodes
! if global node is in my subdomain, assigned 1, else remains 0
         call pp_mark_sub_nodes(isub,nelem,iets,liets,inet,linet,nnet,lnnet,&
                                kinet,lkinet,kmynodes,lkmynodes)
         do inodi = 1,nnodi
            inod = igingn(inodi)
            if (kmynodes(inod).eq.1) then
               kinodes(inod) = kinodes(inod) + 1
               indsub = kinodes(inod)
               sublist(inodi,indsub) = isub
            end if
         end do
      end do

!      do i = 1,nnod
!         write(*,*) (sublist(i,j), j = 1,nmax)
!      end do

! Identify topology of globs
      write(*,'(a)') 'Identify globs by equivalence classes...'
      lnodeglob = nnodi
      allocate(nodeglob(lnodeglob))
      nodeglob = -1
      iglob = 0
      do inodi = 1,nnodi
         inod = igingn(inodi)
         if (count(sublist(inodi,:).ne.0).le.1) then
            call error('GETGLOBS','subdomain list - only one subdomain recognized at interface node')
         end if
         if(nodeglob(inodi).eq.-1) then
! Node was not yet assigned into a glob - so there is a new glob to define
            iglob = iglob + 1
            do jnodi = inodi,nnodi
               jnod = igingn(jnodi)
! Go through the remaining nodes
               if (all(sublist(jnodi,:).eq.sublist(inodi,:))) then
                  nodeglob(jnodi) = iglob
               end if
            end do
         end if
      end do
      nglob = iglob
      write(*,'(a,i8)') 'Total number of globs is',nglob
      if (count(nodeglob.eq.-1).gt.0) then
         call error('GETGLOBS','in nodeglob - unassociated node of interface with glob')
      end if
      
! glob type - 1 face, 2 edge, 3 vertex 
      lglobs = nglob
      allocate(globs(lglobs))
      globs(:)%itype = 0
      globs(:)%selected = .false.
      do iglob = 1,nglob
!        count nodes in each glob and allocate memory for nodes
         globs(iglob)%nnod = count(nodeglob.eq.iglob)
         allocate(globs(iglob)%nodes(globs(iglob)%nnod))
!        find nodes in the glob         
         iglobnode = 0
         do inodi = 1,nnodi
            if (nodeglob(inodi).eq.iglob) then
               iglobnode = iglobnode + 1
               globs(iglob)%nodes(iglobnode) = inodi
            end if
         end do
!        count subdomains in each glob and allocate memory for their list
         globs(iglob)%nsub = count(sublist(globs(iglob)%nodes(1),:).gt.0)
         allocate(globs(iglob)%subdomains(globs(iglob)%nsub))
!        find subdomains in the glob         
         globs(iglob)%subdomains(:) = sublist(globs(iglob)%nodes(1),1:globs(iglob)%nsub)
         if (globs(iglob)%nsub.eq.2) then
!           such entity is a face
            globs(iglob)%itype = 1
         else if (globs(iglob)%nsub.gt.2.and.globs(iglob)%nnod.gt.1) then
!           such entity is an edge
            globs(iglob)%itype = 2
         else if (globs(iglob)%nsub.gt.2.and.globs(iglob)%nnod.eq.1) then
!           such entity is a vertex
            globs(iglob)%itype = 3
         else
            call warning('GETGLOBS','Unable to determine type of glob.')
         end if
      end do
      deallocate(sublist)
      deallocate(nodeglob)

      nface   = count(globs%itype.eq.1)
      nedge   = count(globs%itype.eq.2)
      nvertex = count(globs%itype.eq.3)
      if (nface + nedge + nvertex .ne. nglob) then
         call error('GETGLOBS','Total number of globs does not match.')
      end if

      write(*,*) 'Situation in globs after equivalence class algorithm:'
      write(*,'(a,i8)') 'nface   = ',nface
      write(*,'(a,i8)') 'nedge   = ',nedge
      write(*,'(a,i8)') 'nvertex = ',nvertex
      write(*,'(a)')    '--------------------'
      write(*,'(a,i8)') 'total   = ',nface + nedge + nvertex

      if (debug) then
         do iglob = 1,nglob
            write(*,*) 'glob=============', iglob
            write(*,*) globs(iglob)%itype
            write(*,*) globs(iglob)%nnod
            write(*,*) globs(iglob)%nodes
            write(*,*) globs(iglob)%nsub
            write(*,*) globs(iglob)%subdomains
         end do
      end if

! Prepare space for new vertices
      lnewvertex = nnodi
      allocate(newvertex(lnewvertex))
      newvertex = 0

! mark current vertices into the newvertex array
      do iglob = 1,nglob
         if (globs(iglob)%itype .eq. 3) then
            newvertex(globs(iglob)%nodes) = 1
         end if
      end do

! generate list of neighbouring subdomains for further algorithm
      lsubneib = nsub
      allocate(subneib(lsubneib))
      lneighbourings = nsub
      allocate(neighbourings(lneighbourings))
      do isub = 1,nsub
         subneib = 0
         do iglob = 1,nglob
            ! only faces are interesting
            if (globs(iglob)%itype .eq. 1 .and. any(globs(iglob)%subdomains .eq. isub)) then
               subneib(globs(iglob)%subdomains) = 1
            end if
         end do
         ! my subdomain is also marked
         nnsub = count(subneib.eq.1) - 1
         neighbourings(isub)%nnsub = nnsub
         allocate(neighbourings(isub)%list(nnsub))
         iisub = 0
         do jsub = 1,nsub
            if (subneib(jsub).eq.1 .and. jsub.ne.isub) then
               iisub = iisub + 1
               neighbourings(isub)%list(iisub) = jsub
            end if
         end do
      end do
      if (debug) then
         do isub = 1,nsub
            write(*,*) 'isub',isub,'neighbours',neighbourings(isub)%list
         end do
      end if

      ! coordinates vector for nodes
      lxyzbase = ndim
      allocate (xyzbase(lxyzbase))

      write(*,*) 'Correction of sufficient number of corners on each subdomain...'
      do isub = 1,nsub
! Creation of field kmynodes(nnod) - local usage of global nodes
! if global node is in my subdomain, assigned 1, else remains 0
         call pp_mark_sub_nodes(isub,nelem,iets,liets,inet,linet,nnet,lnnet,&
                                kinet,lkinet,kmynodes,lkmynodes)
         do jsub = 1,neighbourings(isub)%nnsub
            indjsub = neighbourings(isub)%list(jsub)
            kinodes = 0
! Creation of field kneibnodes(nnod) - local usage of global nodes
! if global node is in neighbouring subdomain, assigned 1, else remains 0
            call pp_mark_sub_nodes(indjsub,nelem,iets,liets,inet,linet,nnet,lnnet,&
                                   kinet,lkinet,kneibnodes,lkneibnodes)
            where( kmynodes.eq.1 .and. kneibnodes.eq.1 ) kinodes = 1
            nshared = count(kinodes.eq.1)
! number of corners already shared by the two subdomains
            nnodcpair = 0
            do inodi = 1,nnodi
               if (newvertex(inodi).gt.0 .and. kinodes(igingn(inodi)).eq.1) then
                  nnodcpair = nnodcpair + 1
               end if
            end do
            if (debug) then
               write (*,*) 'Subdomains',isub,' and ',indjsub,'already share ',nnodcpair,' corners.'
            end if
            ! generate at least 3 corners on the common face
            ! prepare field of local interface nodes
            lisingin = nshared
            allocate(isingin(lisingin))
            indisn = 0
            do inodi = 1,nnodi
               inod = igingn(inodi)
               if (kinodes(inod).eq.1) then
                  indisn = indisn + 1
                  isingin(indisn) = inodi
               end if
            end do
            ! Check you have filled the whole array
            if (indisn.ne.nshared) then
               write(*,*) 'Error in construction of subdomain interface.'
               stop
            end if
            lplayground = nshared
            allocate(playground(lplayground))
            playground = 0
            ! Local coordinates
            lxyzsh1 = nshared
            lxyzsh2 = ndim
            ! localize coordinates of shared nodes
            allocate(xyzsh(lxyzsh1,lxyzsh2))
            xyzsh = xyz(igingn(isingin),:)

            ! mark already known corners in playground
            where (newvertex(isingin) .eq. 1) playground = 1
            if (nnodcpair.eq.0.or..not.minimizecorners) then
               ! no corner on subdomain yet
               ! make it the most remote node to the first interface node
               inodsh = 1

               xyzbase = xyzsh(inodsh,:)
               
               ! Find second corner by maximizing the distance of the first interface node
               ldist = nshared
               allocate(dist(ldist))
               do ish = 1,nshared
                  dist(ish) = sum((xyzsh(ish,:) - xyzbase(:))**2)
               end do
            !   dist = (xis-xbase)**2 + (yis-ybase)**2 + (zis-zbase)**2
               indaux  = maxloc(dist)
               inodcs1 = indaux(1)

               playground(inodcs1) = 1
               nnodcpair = 1

               deallocate(dist)
            end if
            if (nnodcpair.eq.1.or..not.minimizecorners) then
               ! one corner is already set in playground, select the second
               inodcs = 0
               do inodis = 1,nshared
                  if (playground(inodis).eq.1) then
                     inodcs = inodis
                     exit
                  end if
               end do
               if (inodcs.eq.0) then
                  write(*,*) 'Problem finding already known corners in playground.'
                  stop
               end if
               xyzbase = xyzsh(inodcs,:)
               
               ! Find second corner by maximizing the distance from the first one
               ldist = nshared
               allocate(dist(ldist))
               do ish = 1,nshared
                  dist(ish) = sum((xyzsh(ish,:) - xyzbase(:))**2)
               end do
!               dist = (xis-xbase)**2 + (yis-ybase)**2 + (zis-zbase)**2
               indaux  = maxloc(dist)
               inodcs2 = indaux(1)
               if (inodcs.eq.inodcs2) then
                  write(*,*) 'x,y,z'
                  write(*,*) xyzsh(ish,:)
                  write(*,*) 'Problem finding second corner on subdomain - same as first.'
                  write(*,*) 'dist'
                  write(*,*) dist
                  write(*,*) 'first corner', inodcs, ', proposed second corner',inodcs2
                  stop
               end if

               playground(inodcs2) = 1
               nnodcpair = 2

               deallocate(dist)
            end if
            if (ndim.gt.2.and.(nnodcpair.eq.2.or..not.minimizecorners)) then
               ! two corners are already set in playground, select the third
               inodcs1 = 0
               do inodis = 1,nshared
                  if (playground(inodis).eq.1) then
                     inodcs1 = inodis
                     exit
                  end if
               end do
               if (inodcs1.eq.0) then
                  write(*,*) 'Problem finding already known corners in playground.'
                  stop
               end if
               inodcs2 = 0
               do inodis = inodcs1+1,nshared
                  if (playground(inodis).eq.1) then
                     inodcs2 = inodis
                     exit
                  end if
               end do
               if (inodcs2.eq.0) then
                  write(*,*) 'Problem finding already known corners in playground.'
                  stop
               end if
               x1 = xyzsh(inodcs1,1)
               y1 = xyzsh(inodcs1,2)
               z1 = xyzsh(inodcs1,3)
               x2 = xyzsh(inodcs2,1)
               y2 = xyzsh(inodcs2,2)
               z2 = xyzsh(inodcs2,3)
               
               ! Find third corner as the most orthogonal
               larea = nshared
               allocate(area(larea))
!               cosangle = (xis-x1)*(x2-x1) + (yis-y1)*(y2-y1) + (zis-z1)*(z2-z1)
               do ish = 1,nshared
                  xish = xyzsh(ish,1)
                  yish = xyzsh(ish,2)
                  zish = xyzsh(ish,3)
                  area(ish) = ((y2-y1)*(zish-z1) - (yish-y1)*(z2-z1))**2 &
                            + ((x2-x1)*(zish-z1) - (xish-x1)*(z2-z1))**2 &
                            + ((x2-x1)*(yish-y1) - (xish-x1)*(y2-y1))**2 
               end do
               indaux  = maxloc(abs(area))
               inodcs3 = indaux(1)
               if (inodcs1.eq.inodcs3.or.inodcs2.eq.inodcs3) then
                  write(*,*) 'Problem finding third corner on subdomain - same as one of previous two.'
                  write(*,*) 'Already have:',inodcs1,inodcs2,' proposed third corner is', inodcs3
                  write(*,*) 'area'
                  do i = 1,larea
                     write(*,*) i, area(i)
                  end do
                  stop
               end if

               playground(inodcs3) = 1
               nnodcpair = 3

               deallocate(area)
            end if
            deallocate(xyzsh)

            ! now copy playground to global field of new vertices
            newvertex(isingin) = playground
            deallocate(playground)
            deallocate(isingin)
         end do
      end do


! Check subdomain by subdomain that the number of vertices is sufficient
      do isub = 1,nsub
! Creation of field kmynodes(nnod) - local usage of global nodes
! if global node is in my subdomain, assigned 1, else remains 0
         call pp_mark_sub_nodes(isub,nelem,iets,liets,inet,linet,nnet,lnnet,&
                                kinet,lkinet,kmynodes,lkmynodes)
         nnodcs = 0
         do inodi = 1,nnodi
            inod = igingn(inodi)
            if (kmynodes(inod).eq.1.and.newvertex(inodi).eq.1) then
               nnodcs = nnodcs + 1
            end if
         end do
         if (debug) then
            write(*,*) 'Number of corners on subdomain',isub,' is ',nnodcs
         end if
         ! if number of corners is lower than ndim, add some corners
         if (nnodcs.lt.ndim) then
            call warning('getglobs','Number of corners on subdomain lower than dimension.')
            write(*,*) 'Number of corners on subdomain',isub,' is ',nnodcs
            nnodis = count(kmynodes.eq.1.and.kinodes.gt.1)
            write(*,*) 'Number of interface nodes on subdomain',isub,' is ',nnodis
            if (nnodis.eq.0) then
               write(*,*) 'Error in topology - subdomain ',isub,' has zero interface nodes.'
               stop
            end if
         end if
      end do
      deallocate(xyz)
      deallocate(iets)
      deallocate(kinet)
      deallocate(kmynodes,kinodes,kneibnodes)

! remove vertices from list of globs - will be used more times
 70   do inodi = 1,nnodi
         if (newvertex(inodi).eq.1) then
            do iglob = 1,nglob
               where (globs(iglob)%nodes.eq.inodi) globs(iglob)%nodes = 0
            end do
         end if
      end do

! FVS - fixed variables
!  * IFIX(LIFIX) * FIXV(LFIXV) * 
      name = name1(1:lname1)//'.FVS'
      open (unit=idfvs,file=name,status='old',form='formatted')
      rewind idfvs
      lifix = ndof
      allocate(ifix(lifix))
! Read fixed variables
      read(idfvs,*) ifix
      close(idfvs)
! Creation of field KDOF(NNOD) with addresses before first global
! dof of node
      lkdof = nnod
      allocate(kdof(lkdof))
      kdof(1) = 0
      do inod = 2,nnod
         kdof(inod) = kdof(inod-1) + nndf(inod-1)
      end do
! remove Dirichlet boundary conditions from the list of globs
      nremoved_corners_bc = 0
      do inodi = 1,nnodi
         indnod = igingn(inodi)
         inddof = kdof(indnod)
         ndofn = nndf(indnod)
         if (any(ifix(inddof+1:inddof+ndofn).gt.0)) then
            ! remove node from globs
            do iglob = 1,nglob
               where (globs(iglob)%nodes.eq.inodi) globs(iglob)%nodes = 0
            end do
            ! remove node from corners
! uncomment this if you allow corners on Dirichlet BC
!            if (newvertex(inodi).gt.0) then
!               newvertex(inodi) = 0
!               nremoved_corners_bc = nremoved_corners_bc + 1
!            end if
         end if
      end do
      write(*,*) 'Number of removals of corners for Dirichlet BC:', nremoved_corners_bc
! remove zeros from list of globs
      ncorrections = 0
      do iglob = 1,nglob
         nng = globs(iglob)%nnod
         inewnodes = 0
         do ing = 1,nng
            indinode = globs(iglob)%nodes(ing)
            if (indinode.ne.0) then
               inewnodes = inewnodes + 1
               globs(iglob)%nodes(inewnodes) = indinode
            else
               ncorrections = ncorrections + 1
            end if
         end do
         globs(iglob)%nodes(inewnodes+1:) = 0
         globs(iglob)%nnod = inewnodes
      end do
      write(*,*) 'Number of removals of nodes from nodes of globs:',ncorrections

! Mark degenerated globs with negative integer in the field glob_type
      write(*,*) 'I am going to remove ', count(globs%nnod.eq.0.and.globs%itype.ne.-2), &
                 ' degenerated globs (including all vertices).'
      where (globs%nnod.eq.0) globs%itype = -2

! All vertices should be removed by now - check that
      do iglob = 1,nglob
         if (globs(iglob)%itype.eq.3.and.globs(iglob)%nnod.ne.0) then
            write(*,*) 'Strange vertex!', iglob
            write(*,*) 'Number of nodes =', globs(iglob)%nnod
            write(*,*) 'Nodes'
            inodi = globs(iglob)%nodes(1)
            write(*,*) inodi
            write(*,*) 'Number of interface nodes matching',count(igingn.eq.igingn(inodi))
         end if
      end do

      nface   = count(globs%itype.eq.1)
      nedge   = count(globs%itype.eq.2)
      nvertex = count(globs%itype.eq.3)
      write(*,*) 'Situation in globs after corrections:'
      write(*,'(a,i8)') 'nface   = ',nface
      write(*,'(a,i8)') 'nedge   = ',nedge
      write(*,'(a,i8)') 'nvertex = ',nvertex
      write(*,'(a)')    '--------------------'
      write(*,'(a,i8)') 'total   = ',nface + nedge + nvertex
      write(*,'(a)')    '--------------------'
      write(*,'(a,i8)') 'corners = ',count(newvertex.eq.1)

! Perform check of interface coverage - union of faces, edges, corners and Dirichlet BC should form the whole interface 
! moreover, faces, edges and corners should be disjoint
      licheck = nnodi
      allocate(icheck(licheck))
      icheck = 0
      ! mark corners
      icheck = newvertex
      ! go through globs
      do iglob = 1,nglob
         if (globs(iglob)%itype.eq.1.or.globs(iglob)%itype.eq.2) then
            nglobn = globs(iglob)%nnod
            icheck(globs(iglob)%nodes(1:nglobn)) = icheck(globs(iglob)%nodes(1:nglobn)) + 1
         end if
      end do
! mark Dirichlet boundary conditions from the list of globs
      do inodi = 1,nnodi
         indnod = igingn(inodi)
         inddof = kdof(indnod)
         ndofn = nndf(indnod)
         if (any(ifix(inddof+1:inddof+ndofn).gt.0)) then
            icheck(inodi) = 1
         end if
      end do
      ! Every node should be present once
      if (any(icheck.ne.1)) then
         write(*,*) 'Error in interface coverage:'
         do inodi = 1,nnodi
            indnod = igingn(inodi)
            point = kdof(indnod)
            ndofn = nndf(indnod)
            if (icheck(inodi).ne.1) then
               write(*,*) 'Node ',igingn(inodi),' present ',icheck(inodi),' times.'
            end if
         end do
         stop
      end if
      deallocate(icheck)
      write(*,*) 'Coverage of interface by globs successfully checked.'

      deallocate(ifix)
      deallocate(kdof)

! Select globs
! Choosing 
3     print *, 'Already selected: ',count(globs%selected),' globs.'
      print *, 'What do you want to select?'
      print *, '1 - specify faces'
      print *, '2 - specify edges'
      print *, '3 - all faces'
      print *, '4 - all edges'
      print *, '5 - all faces and edges'
      print *, '6 - zero selection'
      print *, '7 - random selection of more corners'
      print *, '8 - import existing file with corners *.CN'
      print *, '0 - export globs and back to main menu'
      print *, '9 - back to main menu (no saving)'
      write (*, '(a,$)') 'Your choice: '
      read *, kselection

      select case(kselection)
      case(1)
         call globselection(1,nglob)
         goto 3
      case(2)
         call globselection(2,nglob)
         goto 3
      case(3)
         where (globs%itype.eq.1) globs%selected = .true.
         goto 3
      case(4)
         where (globs%itype.eq.2) globs%selected = .true.
         goto 3
      case(5)
         where (globs%itype.eq.1) globs%selected = .true.
         where (globs%itype.eq.2) globs%selected = .true.
         goto 3
      case(6)
         globs%selected = .false.
         goto 3
      case(7)
! Random selection if desired ! Check number of corners
         nnodcold = count(newvertex.eq.1)
 2       write (*, '(a,$)') 'Get target number of corners: '
         read *, nnodc
         if (nnodc.gt.nnodi) then
            write(*,*) 'Chosen number of corners higher than number of interface nodes!(',nnodi,')'
            write(*,*) 'Unable to proceed, try again...'
            goto 2
         else if (nnodc.lt.nnodcold) then
            write(*,*) 'Chosen number is lower than number of already selected corners!(',nnodcold,')'
            write(*,*) 'Unable to proceed, try again...'
            goto 2
         end if

! Set random corners form interface
         print *, 'Selecting random corners...'
         do inc = nnodcold+1,nnodc
  50        call RANDOM_NUMBER(rnd)
            indi = int(rnd*nnodi) + 1
            if (newvertex(indi).eq.1) goto 50
            newvertex(indi) = 1
         end do
         print *, '...done.'
         goto 70

      case(8)

! CN - list of corners
         name = name1(1:lname1)//'.CN'
         open (unit=idcn,file=name,status='old',form='formatted', iostat=ios)
         if( ios.ne.0) then
            write(*,*) 'Error opening file ',trim(name),'. Nothing changed.'
            goto 3
         else
            ! read number of corners
            read(idcn,*) nnodc
            ! allocate array for indices of corners
            linodc = nnodc
            allocate(inodc(linodc))
            ! read indices of corners in global node numbering in hat
            read(idcn,*) inodc
            close(idcn)

            ! Mark these corners
            newvertex = 0
            do inc = 1,nnodc
               indn = inodc(inc)
               where(igingn.eq.indn) newvertex = 1
            end do
            deallocate(inodc)
            goto 70
         end if

      case(0)

! GLB - list of globs
         name = name1(1:lname1)//'.GLB'
         open (unit=idglb,file=name,status='replace',form='formatted')

! transform the globs into pmd-style arrays
         nglb = count(globs%selected)
! determine length of array inglb
         linglb = 0
         do iglob = 1,nglob
            if (globs(iglob)%selected) then
               nng  = globs(iglob)%nnod
               linglb = linglb + nng
            end if
         end do
! prepare arrays
         lnnglb = nglb
         linglb = linglb
         allocate(nnglb(lnnglb),inglb(linglb))
! fill arrays 
         iglb   = 0
         iinglb = 0
         do iglob = 1,nglob
            if (globs(iglob)%selected) then
               iglb = iglb + 1
               nng  = globs(iglob)%nnod
               nnglb(iglb) = nng
               do ing = 1,nng
                  iinglb = iinglb + 1
                  inglb(iinglb) = igingn(globs(iglob)%nodes(ing))
               end do
            end if
         end do

! write into the GLB file
         write(idglb,'(2i15)') nglb, linglb
         write(idglb,'(4i15)') inglb
         write(idglb,'(4i15)') nnglb
         
         close(idglb)
         print *, 'Created file ', trim(name),' with a list of selected ',nglb,' globs.'

         deallocate(nnglb,inglb)

      case(9)
         continue
      case default
         write(*,*) 'Unknown answer, try again.'
         goto 3
      end select

! CN - list of corner nodes
      name = name1(1:lname1)//'.CN'
      open (unit=idcn,file=name,status='replace',form='formatted')

      nnewvertex = count(newvertex.eq.1)
      write(idcn,*) nnewvertex
      do i = 1,nnodi
         if (newvertex(i).eq.1) then
            write(idcn,*) igingn(i)
         end if
      end do
      if (nnewvertex.eq.0) then
         write(idcn,*) ' 0'
      end if
      print *, 'Created file ', trim(name), ' with a list of ',nnewvertex,' corners.'
      close(idcn)


! List interface nodes to file *.INT
! INT - interface
      name = name1(1:lname1)//'.INT'
      open (unit=idint,file=name,status='replace',form='formatted')
      write(idint,*) nnodi
      write(idint,*) igingn
      print *, 'Created file ', trim(name), ' with list of',nnodi,' interface nodes.'
      close(idint)

! Pairs
      ! prepare data about pairs
      npair = count(globs%itype.eq.1.and. globs%selected)

! PAIR - pairs of globs to compute by adaptivity
      name = name1(1:lname1)//'.PAIR'
      open (unit=idpair,file=name,status='replace',form='formatted')
      write(idpair,*) npair
      ipair = 0
      iglb  = 0
      do iglob = 1,nglob
         if (globs(iglob)%selected) then
            iglb = iglb + 1
         end if
         if (globs(iglob)%itype.eq.1 .and. globs(iglob)%selected) then
            ipair = ipair + 1

            write (idpair,*) iglb, globs(iglob)%subdomains(1),globs(iglob)%subdomains(2)
         end if
      end do
      ! check the number
      if (ipair .ne. npair) then
         write(*,*) 'GETGLOBS: Check of number of pairs to export failed.'
         stop
      end if
      print *, 'Created file ', trim(name), ' with list of',npair,' pairs for adaptivity.'
      close(idpair)
      

! Clear memory
      do iglob = 1,nglob
         deallocate(globs(iglob)%nodes)
         deallocate(globs(iglob)%subdomains)
      end do
      deallocate(globs)
      do isub = 1,nsub
         deallocate(neighbourings(isub)%list)
      end do
      deallocate(neighbourings)
      deallocate(newvertex)
      deallocate(subneib)
      deallocate(xyzbase)
      deallocate(igingn)
      deallocate(inet,nnet,nndf)

      return

end subroutine getglobs

!***************************************************************************************
subroutine globselection(iglobtype,nglob)
! Subroutine that helps in selecting globs
!***************************************************************************************
implicit none

integer,intent(in) :: iglobtype
integer,intent(in) :: nglob

! Local variables
integer :: indglob, iglob

character(10) :: string
character(1) :: yn


  55  write(*,'(a,$)') 'Do you want to see the list (y/n)?'
      read(*,*) yn
      select case (yn)
         case ('Y','y')
            write(*,*) 'selected, glob number, number of nodes, subdomains'
            do iglob = 1,nglob
               if (globs(iglob)%itype.eq.iglobtype) then
                  write(*,*) globs(iglob)%selected, iglob, globs(iglob)%nnod, globs(iglob)%subdomains
               end if
            end do
         case ('N','n')
            continue
         case default
            print *, 'Unknown answer, try again ...'
            goto 55
      end select

! Choosing 
      string = ' '
      print *, 'Select:'
      print *, 'All - all globs of the type'
      print *, 'No - no globs of the type'
      print *, 'number of a glob'
      print *, '0 - back to menu'
      write (*, '(a,$)') 'Your choice: '
      read *, string

      if (string(1:1).eq.'a'.or.string(1:1).eq.'A') then
! select all globs of the type
         where(globs%itype.eq.iglobtype) globs%selected = .true.
         goto 55
      else if (string(1:1).eq.'n'.or.string(1:1).eq.'N') then
! deselect all globs of the type
         where(globs%itype.eq.iglobtype) globs%selected = .false.
         goto 55
      else if (string(1:1).eq.'0') then
         return
      else
         read(string,*) indglob
         globs(indglob)%selected = .true.
         goto 55
      end if
end subroutine globselection

!***********************************************************************
subroutine create_sub_files(problemname)
!***********************************************************************
!     Subroutine for creation of subdomain files
!***********************************************************************
use module_utils
implicit none

! name of problem
character(*),intent(in) :: problemname

integer ::            linets,   lnnets,   lnndfs,   lifixs,   lkdofs,   lisngns
integer,allocatable :: inets(:), nnets(:), nndfs(:), ifixs(:), kdofs(:), isngns(:)
integer ::            liins,   liivsvns,   liovsvns
integer,allocatable :: iins(:), iivsvns(:), iovsvns(:)
integer ::            lglobal_corner_numbers,   licnsins
integer,allocatable :: global_corner_numbers(:), icnsins(:)
integer::            linglb,   lnnglb
integer,allocatable:: inglb(:), nnglb(:)
integer::            lglobal_glob_numbers
integer,allocatable:: global_glob_numbers(:)
integer::            lnglobvars
integer,allocatable:: nglobvars(:)
integer::            ligvsivns1, ligvsivns2
integer,allocatable:: igvsivns(:,:)

integer ::            lrhss,   lfixvs,   lxyzs1, lxyzs2
real(kr),allocatable:: rhss(:), fixvs(:), xyzs(:,:)

character(100) :: filename

! local variables
integer :: nglb
integer:: isub, ie, inc, ndofn, nne, &
          inod, idofn, indifixs, indifix, indnc, indnnets, indrhs, indrhss, &
          inods, jsub, pointinet, pointinets, i, j, indn, idofis, idofos, inodis, &
          indns, indins, inodcs, iglbn, nglbn, iglb, iglobs, indn1,&
          nglbv, pointinglb, indivs, iglbv, ivar, indng, indvs
integer:: nadjs, nnodcs, nnods, nelems, ndofs, nnodis, ndofis, ndofos, nglobs
integer:: idsmd

! GMIS - basic mesh data - structure:
!  * INET(LINET) * NNET(LNNET) * NNDF(LNNDF) * XYF(LXYF) *
      name = name1(1:lname1)//'.GMIS'
      open (unit=idgmi,file=name,status='old',form='formatted')
      rewind idgmis

! FVS - fixed variables
!  * IFIX(LIFIX) * FIXV(LFIXV) * 
      name = name1(1:lname1)//'.FVS'
      open (unit=idfvs,file=name,status='old',form='formatted')

! RHS - Right Hand Side
      name = name1(1:lname1)//'.RHS'
      open (unit=idrhs,file=name,status='old',form='unformatted')

! ES - list of global element numbers in subdomains - structure:
      name = name1(1:lname1)//'.ES'
      open (unit=ides,file=name,status='old',form='formatted')

! CN - list of global corner nodes - structure:
      name = name1(1:lname1)//'.CN'
      open (unit=idcn,file=name,status='old',form='formatted')

! Import global data******************************
! Mesh data
! read fields INET, NNET, NNDF from file IDGMI
      lnnet = nelem
      lnndf = nnod
      lkinet= nelem
      lkdof = nnod
      allocate(inet(linet),nnet(lnnet),nndf(lnndf),x(lx),y(ly),z(lz),&
               kinet(lkinet),kdof(lkdof))
      lxyz1 = nnod
      lxyz2 = ndim
      allocate(xyz(lxyz1,lxyz2))
      read(idgmi,*) inet
      read(idgmi,*) nnet
      read(idgmi,*) nndf
! read fields INET, NNET, NNDF from file IDGMI
      read(idgmi,*) xyz
      close(idgmi)

! Creation of field KINET(NELEM) with addresses before first global node of element IE in field inet
      kinet(1) = 0
      do ie = 2,nelem
         kinet(ie) = kinet(ie-1) + nnet(ie-1)
      end do
! create array KDOF
      kdof(1) = 0
      do inod = 2,nnod
         kdof(inod) = kdof(inod-1) + nndf(inod-1)
      end do

! Boundary conditions
! read fields IFIX
      lifix = ndof
      lfixv = ndof
      allocate(ifix(lifix),fixv(lfixv))
      read(idfvs,*) ifix
      read(idfvs,*) fixv

! Right hand side
! read fields RHS
      lrhs = ndof
      allocate(rhs(lrhs))
      read(idrhs) rhs

! Division of elements into subdomains
! read field IETS from file IDES
      liets = nelem
      allocate(iets(liets))
      read(ides,*) iets

! Corners
! read file CN
      read(idcn,*) nnodc
      linodc = nnodc
      allocate(inodc(linodc))
      read(idcn,*) inodc

! Globs
! GLB - list of globs
      name = name1(1:lname1)//'.GLB'
      open (unit=idglb,file=name,status='old',form='formatted')
      read(idglb,*) nglb, linglb
      linglb = linglb
      lnnglb = nglb
      allocate (inglb(linglb),nnglb(lnnglb))
      read(idglb,*) inglb
      read(idglb,*) nnglb
      close(idglb)

! End of reading global data**********************

! Prepare subdomain data
      lkmynodes   = nnod
      lkadjsnodes = nnod
      lkinodes    = nnod
      allocate(kmynodes(lkmynodes),kinodes(lkinodes),kadjsnodes(lkadjsnodes))

! Recognize interface
      kinodes = 0
      do isub = 1,nsub
! Creation of field kmynodes(nnod) - local usage of global nodes
! if global node is in my subdomain, assigned 1, else remains 0
         call pp_mark_sub_nodes(isub,nelem,iets,liets,inet,linet,nnet,lnnet,&
                                kinet,lkinet,kmynodes,lkmynodes)
! Mark interface nodes
         where(kmynodes.eq.1) kinodes = kinodes + 1
      end do

! Begin loop over subdomains
      do isub = 1,nsub

! Creation of field kmynodes(nnod) - local usage of global nodes
         call pp_mark_sub_nodes(isub,nelem,iets,liets,inet,linet,nnet,lnnet,kinet,lkinet,kmynodes,lkmynodes)

! find local number of nodes on subdomain NNODS
         nnods = count(kmynodes.ne.0)

! find local number of DOF on subdomain NDOFS
         ndofs = sum(kmynodes*nndf)


! find number of adjacent subdomains NADJS
         nadjs = 0
         do jsub = 1,nsub
            if (jsub.ne.isub) then
               call pp_mark_sub_nodes(jsub,nelem,iets,liets,inet,linet,nnet,lnnet,&
                                      kinet,lkinet,kadjsnodes,lkadjsnodes)
               if (any(kmynodes.eq.1.and.kadjsnodes.eq.1)) then
                  nadjs = nadjs + 1
               end if
            end if
         end do

! find local number of elements on subdomain NELEMS
         nelems = count(iets.eq.isub)

! create subdomain description of MESH
         linets = 0
         do ie = 1,nelem
            if (iets(ie).eq.isub) then
               linets = linets + nnet(ie)
            end if
         end do
         lnnets = nelems
         lnndfs = nnods
         lkdofs = nnods
         allocate(inets(linets),nnets(lnnets),nndfs(lnndfs),kdofs(lkdofs))
         lisngns = nnods
         allocate(isngns(lisngns))
         pointinet  = 0
         pointinets = 0
         indnnets   = 0
         do ie = 1,nelem
            nne = nnet(ie)
            if (iets(ie).eq.isub) then
               inets(pointinets+1:pointinets+nne) = -inet(pointinet+1:pointinet+nne)
               indnnets = indnnets + 1
               nnets(indnnets) = nne
               pointinets = pointinets + nne
            end if
            pointinet = pointinet + nne
         end do
         inods = 0
         do inod = 1,nnod
            if (kmynodes(inod).eq.1) then
               inods = inods + 1
               nndfs(inods) = nndf(inod)
               isngns(inods) = inod
               ! fix inet
               where (inets.eq.-inod) inets = inods
            end if
         end do
! create array kdofs
         kdofs(1) = 0
         do inods = 2,nnods
            kdofs(inods) = kdofs(inods-1) + nndfs(inods-1)
         end do

! find number of interface nodes
         nnodis = 0 
         ndofis = 0
         ndofos = 0
         do inods = 1,nnods
            indn = isngns(inods)

            if (kinodes(indn).gt.1) then
               nnodis = nnodis + 1
               ndofis = ndofis + nndfs(inods)
            else
               ndofos = ndofos + nndfs(inods)
            end if
         end do
! generate mapping of interface nodes to subdomain nodes and the same for dofs 
         liins = nnodis
         allocate(iins(liins))
         liivsvns = ndofis
         allocate(iivsvns(liivsvns))
         liovsvns = ndofos
         allocate(iovsvns(liovsvns))

         inodis = 0
         idofis = 0
         idofos = 0
         do inods = 1,nnods
            indn = isngns(inods)
            ndofn  = nndfs(inods)

            if (kinodes(indn).gt.1) then
               inodis = inodis + 1

               iins(inodis) = inods
               do idofn = 1,ndofn 
                  idofis = idofis + 1

                  iivsvns(idofis) = kdofs(inods) + idofn
               end do
            else
               do idofn = 1,ndofn 
                  idofos = idofos + 1

                  iovsvns(idofos) = kdofs(inods) + idofn
               end do
            end if
         end do

! find number of coarse nodes on subdomain NNODCS
         nnodcs = 0
         do inc = 1,nnodc
            indnc = inodc(inc)
            if (kmynodes(indnc).eq.1) then
               nnodcs = nnodcs + 1
            end if
         end do

! find mapping of corners
         lglobal_corner_numbers = nnodcs
         allocate(global_corner_numbers(lglobal_corner_numbers))
         licnsins = nnodcs
         allocate(icnsins(licnsins))

         inodcs = 0
         do inc = 1,nnodc
            indnc = inodc(inc)
            if (kmynodes(indnc).eq.1) then
               inodcs = inodcs + 1

               ! mapping to global corner numbers
               global_corner_numbers(inodcs) = inc
               ! mapping to subdomain interface numbers
               call get_index(indnc,isngns,lisngns,indns)
               if (indns .eq. -1) then
                  write(*,*) 'CREATE_SUB_FILES: Index of subdomain node not found.'
                  stop
               end if
               call get_index(indns,iins,liins,indins)
               if (indins .eq. -1) then
                  write(*,*) 'CREATE_SUB_FILES: Index of subdomain interface node not found.'
                  stop
               end if
               icnsins(inodcs) = indins
            end if
         end do

! find local number of globs NGLOBS
         nglobs     = 0
         pointinglb = 0
         do iglb = 1,nglb
            nglbn = nnglb(iglb)

            ! touch first node in glob
            indn1 = inglb(pointinglb + 1)

            if (kmynodes(indn1).eq.1) then
               nglobs = nglobs + 1
            end if

            pointinglb = pointinglb + nglbn
         end do

! mapping of globs
         lglobal_glob_numbers = nglobs
         allocate(global_glob_numbers(lglobal_glob_numbers))
         lnglobvars = nglobs
         allocate(nglobvars(lnglobvars))

         iglobs     = 0
         pointinglb = 0
         do iglb = 1,nglb
            nglbn = nnglb(iglb)

            ! touch first node in glob
            indn1 = inglb(pointinglb + 1)

            if (kmynodes(indn1).eq.1) then

               iglobs = iglobs + 1

               nglbv = 0
               do iglbn = 1,nglbn
                  ndofn = nndf(pointinglb + iglbn)

                  nglbv = nglbv + ndofn
               end do

               nglobvars(iglobs) = nglbv
               global_glob_numbers(iglobs) = iglb
            end if

            pointinglb = pointinglb + nglbn
         end do

         ligvsivns1 = nglobs
         ligvsivns2 = maxval(nglobvars)
         allocate(igvsivns(ligvsivns1,ligvsivns2))
         iglobs     = 0
         pointinglb = 0
         do iglb = 1,nglb
            nglbn = nnglb(iglb)

            ! touch first node in glob
            indn1 = inglb(pointinglb + 1)

            if (kmynodes(indn1).eq.1) then

               iglobs = iglobs + 1

               iglbv = 0
               do iglbn = 1,nglbn
                  
                  indng  = inglb(pointinglb + iglbn)
                  ndofn = nndf(indn)

                  call get_index(indng,isngns,lisngns,indns)
                  if (indns .eq. -1) then
                     write(*,*) 'CREATE_SUB_FILES: Index of subdomain node not found.'
                     stop
                  end if

                  do idofn = 1,ndofn
                     iglbv = iglbv + 1

                     indvs = kdofs(indns) + idofn
                     call get_index(indvs,iivsvns,liivsvns,indivs)
                     if (indivs .eq. -1) then
                        write(*,*) 'CREATE_SUB_FILES: Index of subdomain interface dof not found.'
                        write(*,*) 'indng =',indng,'indns =',indns,'indvs = ',indvs,'indivs = ',indivs, 'isub = ',isub
                        write(*,*) 'iivsvns = ',iivsvns
                        stop
                     end if

                     igvsivns(iglobs,iglbv) = indivs
                  end do
               end do
            end if

            pointinglb = pointinglb + nglbn
         end do



! coordinates of subdomain nodes
         lxyzs1 = nnods
         lxyzs2 = ndim
         allocate(xyzs(lxyzs1,lxyzs2))
         do j = 1,ndim
            do i = 1,nnods
               indn = isngns(i)

               xyzs(i,j) = xyz(indn,j)
            end do
         end do

! make subdomain boundary conditions - IFIXS and FIXVS
         lifixs = ndofs
         lfixvs = ndofs
         allocate(ifixs(lifixs),fixvs(lfixvs))
         ifixs = 0
         fixvs = 0._kr
         indifixs = 0
         do inod = 1,nnod
            if (kmynodes(inod).eq.1) then
               ndofn   = nndf(inod)
               indifix = kdof(inod)
               do idofn = 1,ndofn
                  indifixs = indifixs + 1
                  indifix  = indifix + 1
                  if (ifix(indifix).ne.0) then
                     ifixs(indifixs) = indifixs
                     fixvs(indifixs) = fixv(indifix)
                  end if
               end do
            end if
         end do

         lrhss = ndofs
         allocate(rhss(lrhss))
         indrhss = 0
         do inod = 1,nnod
            if (kmynodes(inod).eq.1) then
               ndofn   = nndf(inod)
               indrhs  = kdof(inod)
               do idofn = 1,ndofn
                  indrhss = indrhss + 1
                  indrhs  = indrhs  + 1
                  rhss(indrhss) = rhs(indrhs)
               end do
            end if
         end do

! ---export subdomain file
         ! open subdomain SMD file with mesh description
         call getfname(problemname,isub,'SMD',filename)
         write(*,*) 'CREATE_SUB_FILES: Creating file: ', trim(filename)
         call allocate_unit(idsmd)
         open (unit=idsmd,file=trim(filename),status='replace',form='formatted')

         ! ---header
         write(idsmd,*) nelems, nnods, ndofs, ndim

         ! ---NNDF array
         ! ---NNET array
         write(idsmd,*) nndfs
         write(idsmd,*) nnets

         ! ---INET array
         write(idsmd,*) inets

         ! ---ISNGN array
         write(idsmd,*) isngns

         ! --- coordinates
         write(idsmd,*) xyzs

         ! ---interface data
         write(idsmd,*) nnodis
         write(idsmd,*) ndofis
         write(idsmd,*) ndofos
         write(idsmd,*) iins
         write(idsmd,*) iivsvns
         write(idsmd,*) iovsvns

         ! --- boundary conditions
         write(idsmd,*) ifixs
         write(idsmd,*) fixvs

         ! --- corners 
         write(idsmd,*) nnodcs
         write(idsmd,*) global_corner_numbers
         write(idsmd,*) icnsins

         ! --- globs
         write(idsmd,*) nglobs
         write(idsmd,*) global_glob_numbers
         write(idsmd,*) nglobvars
         do iglobs = 1,nglobs
            ! glob variables
            write(idsmd,*) (igvsivns(iglobs,ivar),ivar = 1,nglobvars(iglobs))
         end do

         close(idsmd)

         deallocate(igvsivns)
         deallocate(nglobvars)
         deallocate(global_glob_numbers)
         deallocate(icnsins)
         deallocate(global_corner_numbers)
         deallocate(iovsvns)
         deallocate(iivsvns)
         deallocate(iins)
         deallocate(xyzs)
         deallocate(rhss)
         deallocate(ifixs,fixvs)
         deallocate(inets,nnets,nndfs,kdofs)
         deallocate(isngns)

! Finish loop over subdomains
      end do

! Clear memory
      deallocate(inglb,nnglb)
      deallocate(inodc)
      deallocate(ifix,fixv)
      deallocate(kinet)
      deallocate(kdof)
      deallocate(kmynodes,kinodes,kadjsnodes)
      deallocate(inet,nnet,nndf)
      deallocate(iets)
      deallocate(xyz)
      deallocate(rhs)

      return
end subroutine create_sub_files
         
!***********************************************************************
subroutine createtilde
!***********************************************************************
!     Subroutine for creation of space W_tilde
!***********************************************************************
implicit none

integer ::            lihntn,   lslavery,   linett,   lnndft
integer,allocatable :: ihntn(:), slavery(:), inett(:), nndft(:)

integer ::           linodc
integer,allocatable:: inodc(:)

integer:: linglb, lnnglb
integer,allocatable:: inglb(:), nnglb(:)

integer:: linglb_loc, lnnglb_loc
integer,allocatable:: inglb_loc(:), nnglb_loc(:)
integer:: lkmynodest
integer,allocatable:: kmynodest(:)
integer::            lndofsa
integer,allocatable:: ndofsa(:)

! local variables
integer:: isub, indc, nnodt, lsolt, inodt, in, indtilde, nelems, linets, pinet,&
          nglb, nglb_loc, indinglb, indinglb_loc, indnnglb_loc, indnodet, &
          nglbn, ind, iglb, i, ie, ipoint, indng, inc, indsub, ine, ndofn, nne, &
          ndofs, inod

! Import of basic geometry
! GMIS - basic mesh data - structure:
!  * INET(LINET) * NNET(LNNET) * NNDF(LNNDF) * XYF(LXYF) *
      name = name1(1:lname1)//'.GMIS'
      open (unit=idgmi,file=name,status='old',form='formatted')
      rewind idgmi
      linet = linet
      lnnet = nelem
      lnndf = nnod
      allocate(inet(linet),nnet(lnnet),nndf(lnndf))
! read fields INET, NNET, NNDF from file IDGMI
      read(idgmi,*) inet
      read(idgmi,*) nnet
      read(idgmi,*) nndf

! ES - list of global element numbers in subdomains - structure:
      name = name1(1:lname1)//'.ES'
      open (unit=ides,file=name,status='old',form='formatted')
      liets = nelem
      allocate(iets(liets))
      read(ides,*) iets
      close(ides)

! Creation of field KINET(NELEM) with addresses before first global node of element IE in field inet
      lkinet = nelem
      allocate(kinet(lkinet))
      kinet(1) = 0
      do ie = 2,nelem
         kinet(ie) = kinet(ie-1) + nnet(ie-1)
      end do

! Insert list of corner nodes
      ! CN - list of corner nodes
      name = name1(1:lname1)//'.CN'
      open (unit=idcn,file=name,status='old',form='formatted')
      ! read number of corners
      read(idcn,*) nnodc
      ! allocate array for indices of corners
      linodc = nnodc
      allocate(inodc(linodc))
      ! read indices of corners in global node numbering in hat
      read(idcn,*) inodc
      close(idcn)
      
      lkmynodes = nnod
      lkinodes  = nnod
      allocate(kmynodes(lkmynodes),kinodes(lkinodes))
! Prepare field with interface nodes
      kinodes = 0
! Array of subdomain variables
      lndofsa = nsub
      allocate(ndofsa(lndofsa))
      do isub = 1,nsub
! Creation of field kmynodes(nnod) - local usage of global nodes
! if global node is in my subdomain, assigned 1, else remains 0
         call pp_mark_sub_nodes(isub,nelem,iets,liets,inet,linet,nnet,lnnet,&
                                kinet,lkinet,kmynodes,lkmynodes)
! Mark interface nodes
         where(kmynodes.eq.1) kinodes = kinodes + 1

         ! count subdomain degrees of freedom
         ndofs = 0
         do inod = 1,nnod
            if (kmynodes(inod).eq.1) then
               ndofs = ndofs + nndf(inod)
            end if
         end do
         ndofsa(isub) = ndofs
      end do
! Mark corners in kinodes as zeros
      do i = 1,nnodc
         indc = inodc(i)
         kinodes(indc) = 0
      end do
      
! Determine size of the W_tilde in nodes
      nnodt = nnod + sum(kinodes) - count(kinodes.ne.0)
!      write(*,*) 'nnodt =', nnodt

! Prepare fields for W_hat to W_tilde embedding
      ! ihntn(nnod)    - indices of hat nodes in tilde numbering
      ! slavery(nnodt) - field of dependencies among nodes in W_tilde
      !                - 0 if there is no dependece
      !                - index of master node in tilde numbering otherwise
      lihntn   = nnod
      lslavery = nnodt
      allocate(ihntn(lihntn),slavery(lslavery))
      ihntn   = 0
      slavery = 0
! Create field INETT(LINET) - INET in tilde numbering
!              NNDFT(NNODT) - NNDF in tilde dimension
      linett = linet
      lnndft = nnodt
      allocate(inett(linett),nndft(lnndft))
      inett = 0
      nndft = 0

! Generate W_tilde
      inodt = 0
      do isub = 1,nsub
! Creation of field kmynodes(nnod) - local usage of global nodes
! if global node is in my subdomain, assigned 1, else remains 0
         call pp_mark_sub_nodes(isub,nelem,iets,liets,inet,linet,nnet,lnnet,&
                                kinet,lkinet,kmynodes,lkmynodes)
! Generate the mapping for everything except corners
         do in = 1,nnod
            if (kmynodes(in).gt.0.and.count(inodc.eq.in).eq.0) then
            ! node is in subdomain and is not a corner
            ! it means it is a unique in W_tilde
               inodt = inodt + 1
               if (ihntn(in).eq.0) then
                  ! it has not been assigned to mapping of hat into tilde yet, do it
                  ihntn(in) = inodt
               else
                  ! the node has already been assigned in mapping - it is an
                  ! interface node - create dependence in slavery field
                  slavery(inodt)     = ihntn(in)
                  slavery(ihntn(in)) = ihntn(in)
               end if
               ! make a note 
               kmynodes(in) = -inodt
            end if
         end do
! Generate INETT for nodes that are unique to subdomain
         do ie = 1,nelem
            indsub = iets(ie)
            if (indsub.eq.isub) then
               nne = nnet(ie)
               ipoint = kinet(ie)
               do ine = 1,nne
                  indng = inet(ipoint + ine)
                  inett(ipoint + ine) = kmynodes(indng)
               end do
            end if
         end do
      end do
! Generate mapping of hat into tilde for corners - put them at the end of new list
      do inc = 1,nnodc
         inodt = inodt + 1
         indc = inodc(inc)
         ihntn(indc) = inodt
         slavery(inodt) = inodt
         ! Finish array INETT for corners
         where (inet.eq.indc) inett = -inodt
      end do
! Array INETT now has negative entries, reverse it
      inett = -inett
! Check that the numbering reached the absolute value of nodes in tilde
      if (inodt.ne.nnodt) then
         write(*,*) 'Error in creation of W_tilde mapping'
         stop
      end if
! Check that the all nodes in hat are associated with nodes in tilde
      if (count(ihntn.eq.0).ne.0) then
         write(*,*) 'Some nodes in W_hat not associated with W_tilde'
         stop
      end if
! Check that whole INETT is nonzero
      if (count(inett.eq.0).ne.0) then
         write(*,*) 'Some nodes in INETT are zero'
         stop
      end if

! Create field nndft
      do i = 1,nnod
         ndofn = nndf(i)
         indtilde = ihntn(i)
         nndft(indtilde) = ndofn
         where (slavery.eq.indtilde) nndft = ndofn
      end do
!     Check that whole NNDFT is nonzero
      if (count(nndft.eq.0).ne.0) then
         write(*,*) 'Some nodes in NNDFT are zero'
         stop
      end if

! Determine lenght of solution in W_tilde
      lsolt = sum(nndft)

! Generate file with description of W_tilde
      ! GMIST 
      do isub = 1,nsub
         ! Find nelems
         nelems = count(iets.eq.isub)
         ! Find linets
         linets = 0
         do ie = 1,nelem
            if (iets(ie).eq.isub) then
               linets = linets + nnet(ie)
            end if
         end do

         ! open GMISTS
         call bddc_getfname(name1,lname1,isub,'GMISTS',fname)
         open (unit=idgmist,file=fname,status='replace',form='formatted')

         ! open GMISS
         call bddc_getfname(name1,lname1,isub,'GMISS',fname)
         open (unit=idgmis,file=fname,status='replace',form='formatted')

         ! write header
         write(idgmist,*)  nnodt, lsolt
         write(idgmist,*)  nelems, linets, ndofsa(isub)

         write(idgmis,*)   nnod, ndof
         write(idgmis,*)   nelems, linets

         ! write inetst
         pinet = 0
         do ie = 1,nelem
            nne = nnet(ie)
            if (iets(ie).eq.isub) then
               write(idgmist,*) inett(pinet+1:pinet+nne)
               write(idgmis,*)  inet(pinet+1:pinet+nne)
            end if
            pinet = pinet + nne
         end do

         ! write nnetst
         do ie = 1,nelem
            if (iets(ie).eq.isub) then
               write(idgmist,*) nnet(ie)
               write(idgmis,*)  nnet(ie)
            end if
         end do

         write(idgmist,*) nndft

         write(idgmis,*)  nndf

         write(idgmist,*) slavery
         write(idgmist,*) ihntn
         close(idgmist)
         close(idgmis)
      end do
      write(*,*) 'Generated files with description of W_tilde and W_hat spaces.'

! Localize globs to subdomains
! GLB - list of globs
      name = name1(1:lname1)//'.GLB'
      open (unit=idglb,file=name,status='old',form='formatted')
      read(idglb,*) nglb, linglb
      linglb = linglb
      lnnglb = nglb
      allocate (inglb(linglb),nnglb(lnnglb))
      read(idglb,*) inglb
      read(idglb,*) nnglb
      close(idglb)

      lkmynodest = nnodt
      allocate(kmynodest(lkmynodest))

      do isub = 1,nsub
         ! mark subdomain nodes
         call pp_mark_sub_nodes(isub,nelem,iets,liets,inett,linet,nnet,lnnet,&
                                kinet,lkinet,kmynodest,lkmynodest)
         ! count subdomain globs
         nglb_loc   = 0
         linglb_loc = 0
         indinglb   = 0
         do iglb = 1,nglb
            nglbn = nnglb(iglb)
            indnodet = ihntn(inglb(indinglb+1))
            if (any(kmynodest.eq.1.and.slavery.eq.indnodet)) then
               ! subdomain contains this node
               
               nglb_loc   = nglb_loc + 1
               linglb_loc = linglb_loc + nglbn
            end if
            indinglb = indinglb + nglbn
         end do

         ! prepare space for local arrays
         linglb_loc = linglb_loc
         lnnglb_loc = nglb_loc
         allocate(inglb_loc(linglb_loc),nnglb_loc(lnnglb_loc))

         ! generate local arrays INGLB_LOC and NNGLB_LOC
         indinglb     = 0
         indinglb_loc = 0
         indnnglb_loc = 0
         do iglb = 1,nglb
            nglbn = nnglb(iglb)
            indnodet = ihntn(inglb(indinglb+1))
            if (any(kmynodest.eq.1.and.slavery.eq.indnodet)) then
               ! subdomain contains this node
               
               indnnglb_loc = indnnglb_loc + 1
               nnglb_loc(indnnglb_loc) = nglbn
               do i = 1,nglbn
                  indnodet = ihntn(inglb(indinglb + i))
                  ind = 0
                  do inodt = 1,nnodt
                     if (slavery(inodt).eq.indnodet.and.kmynodest(inodt).eq.1) then
                        ind = inodt
                        exit
                     end if
                  end do
                  if (ind.eq.0) then
                     write(*,*) 'Error in glob localization.'
                     stop
                  end if
                  indinglb_loc = indinglb_loc + 1
                  inglb_loc(indinglb_loc) = ind
               end do
            end if
            indinglb = indinglb + nglbn
         end do

         ! write globs into subdomain file
         call bddc_getfname(name1,lname1,isub,'GLB',fname)
         open (unit=idglb,file=fname,status='replace',form='formatted')

         write(idglb,'(2i15)') nglb_loc, linglb_loc
         write(idglb,'(4i15)') inglb_loc
         write(idglb,'(4i15)') nnglb_loc

         close(idglb)
         write(*,'(a,a,a,i6,a,i6)') 'Generated file ',trim(fname), &
                                    ' with local globs in W_tilde space. Subdomain ',isub,',nglb =',nglb_loc

         deallocate(inglb_loc,nnglb_loc)
      end do

! Clear memory
      deallocate(ndofsa)
      deallocate(kinet)
      deallocate(kmynodes,kinodes)
      deallocate(inet,nnet,nndf)
      deallocate(iets)
      deallocate(kmynodest)
      deallocate(inglb,nnglb)
      deallocate(inodc)
      deallocate(ihntn,slavery)
      deallocate(inett,nndft)

      return
end subroutine createtilde
         
!***********************************************************************
subroutine getcorners
!***********************************************************************
!     Subroutine for corner nodes selection
!***********************************************************************
implicit none

integer:: indmax(1), indmin(1), indni, indnc, kselection, &
          indi, indg, inodi, buffer, isub, i, ie,&
          inc, inod, indnodi, j, nmax, nnodi

integer::            ligingn
integer,allocatable:: igingn(:)
integer::            linodc
integer,allocatable:: inodc(:)
integer::            licngin
integer,allocatable:: icngin(:)


real(kr):: xbase, ybase, zbase, scalarp, cangle, distance, rnd, meas
real(kr):: limcangle = 0.99995_kr
real(kr), allocatable:: norm(:), measure(:), xishift(:), yishift(:), zishift(:)

logical:: badangle

! Import of basic geometry
! GMIS - basic mesh data - structure:
!  * INET(LINET) * NNET(LNNET) * NNDF(LNNDF) * XYF(LXYF) *
      name = name1(1:lname1)//'.GMIS'
      open (unit=idgmi,file=name,status='old',form='formatted')
      rewind idgmi
      linet = linet
      lnnet = nelem
      lnndf = nnod
      allocate(inet(linet),nnet(lnnet),nndf(lnndf))
      lx = nnod
      ly = nnod
      lz = nnod
      allocate(x(lx),y(ly),z(lz))
! read fields INET, NNET, NNDF from file IDGMI
      read(idgmi,*) inet
      read(idgmi,*) nnet
      read(idgmi,*) nndf
      read(idgmi,*) x,y,z

! ES - list of global element numbers in subdomains - structure:
      name = name1(1:lname1)//'.ES'
      open (unit=ides,file=name,status='old',form='formatted')
      liets = nelem
      allocate(iets(liets))
      read(ides,*) iets
      close(ides)

! Creation of field KINET(NELEM) with addresses before first global node of element IE in field inet
      lkinet = nelem
      allocate(kinet(lkinet))
      kinet(1) = 0
      do ie = 2,nelem
         kinet(ie) = kinet(ie-1) + nnet(ie-1)
      end do

! Recognize interface
      lkinodes  = nnod
      lkmynodes = nnod
      allocate(kinodes(lkinodes),kmynodes(lkmynodes))
! Begin loop over subdomains
      kinodes = 0
      do isub = 1,nsub
! Creation of field kmynodes(nnod) - local usage of global nodes
! if global node is in my subdomain, assigned 1, else remains 0
         call pp_mark_sub_nodes(isub,nelem,iets,liets,inet,linet,nnet,lnnet,&
                                kinet,lkinet,kmynodes,lkmynodes)
! Mark interface nodes
         where(kmynodes.eq.1) kinodes = kinodes + 1
      end do
! Count interface nodes
      nnodi = count(kinodes.gt.1)
      write(*,'(a,i8)') 'Number of nodes on interface nnodi =',nnodi
! Create field indices of global interface nodes in global numbering
      ligingn = nnodi
      allocate(igingn(ligingn))
      inodi = 0
      do inod = 1,nnod
         if (kinodes(inod).gt.1) then
            inodi = inodi + 1
            igingn(inodi) = inod
         end if
      end do


! Check number of corners
 2    write (*, '(a,$)') 'Get number of corners: '
      read *, nnodc
      if (nnodc.gt.nnodi) then
         write(*,*) 'Chosen number of corners higher than number of interface nodes!'
         write(*,*) 'Unable to proceed, try again...'
         goto 2
      end if

! Choosing the method
 3    print *, 'Method of selecting corner nodes:'
      print *, '1 - random selection from interface nodes'
      print *, '2 - maximizing distance from previously selected nodes'
      print *, '3 - maximizing area of triangles among previously selected nodes'
      write (*, '(a,$)') 'Your choice: '
      read *, kselection

      linodc  = nnodc
      licngin = nnodc
      allocate(inodc(linodc),icngin(licngin))

      select case(kselection)
      case(1)
! Set random corners form interface
         print *, 'Selecting random corners...'
         do inc = 1,nnodc
  50        call RANDOM_NUMBER(rnd)
            indi = int(rnd*nnodi) + 1
            indg = igingn(indi)
            if (count(inodc.eq.indg).gt.0) goto 50
            inodc(inc) = indg
         end do

      case(2)
! Select coordinates of interfaces
         print *, 'Selecting corners by maximizing distance among them...'
         allocate(xishift(nnodi),yishift(nnodi),zishift(nnodi),norm(nnodi),measure(nnodi))
         do i = 1,nnodi
            indnodi = igingn(i)
            xishift(i) = x(indnodi)
            yishift(i) = y(indnodi)
            zishift(i) = z(indnodi)
         end do

! Set first corner
         indmax = maxloc(kinodes)
         nmax   = maxval(kinodes)
         inodc(1) = indmax(1)
         indmin = minloc(abs(igingn-inodc(1)))
         icngin(1) = indmin(1)
         xbase = x(inodc(1))
         ybase = y(inodc(1))
         zbase = z(inodc(1))

! Shift coordinates of interface nodes
         xishift = xishift - xbase
         yishift = yishift - ybase
         zishift = zishift - zbase

! Set second corner
         norm = sqrt(xishift**2 + yishift**2 + zishift**2)
         indmax = maxloc(norm)
         inodc(2) = igingn(indmax(1))
         icngin(2) = indmax(1)
      
! Set the rest of corners
         do inc = 3,nnodc
! Measure triangles
            measure = 0._kr
            do i = 1,nnodi
               indnodi = igingn(i)
               if (count(inodc.eq.indnodi).eq.0) then
                  meas = 0._kr
                  badangle = .false.
                  do j = 2,inc-1
                     indni = icngin(j)
                     distance = sqrt((xishift(i)-xishift(indni))**2 &
                                   + (yishift(i)-yishift(indni))**2 &
                                   + (zishift(i)-zishift(indni))**2 )
                     scalarp = xishift(i)*xishift(indni) + yishift(i)*yishift(indni) + zishift(i)*zishift(indni)
                     cangle = scalarp/(norm(i) * norm(indni))
                     if (abs(cangle).gt.limcangle) badangle = .true.
                     meas = meas + distance
                  end do
                  if (badangle) meas = 0.0_kr
                  measure(i) = meas
               end if
            end do
            if (maxval(measure).eq.0.0_kr) then 
               stop 'Limit for badness probably too low, algorithm for selecting corners fail.'
            end if
           
            indmax = maxloc(measure)
            inodc(inc) = igingn(indmax(1))
            icngin(inc) = indmax(1)
         end do

         deallocate(xishift,yishift,zishift,norm,measure)

      case(3)

         print *, 'Selecting corners by maximizing area of triangle among them...'
! Select coordinates of interfaces
         allocate(xishift(nnodi),yishift(nnodi),zishift(nnodi),norm(nnodi),measure(nnodi))
         do i = 1,nnodi
            indnodi = igingn(i)
            xishift(i) = x(indnodi)
            yishift(i) = y(indnodi)
            zishift(i) = z(indnodi)
         end do

! Set first corner
         indmax = maxloc(kinodes)
         nmax = maxval(kinodes)
         inodc(1) = indmax(1)
         indmin = minloc(abs(igingn-inodc(1)))
         icngin(1) = indmin(1)
         xbase = x(inodc(1))
         ybase = y(inodc(1))
         zbase = z(inodc(1))

! Shift coordinates of interface nodes
         xishift = xishift - xbase
         yishift = yishift - ybase
         zishift = zishift - zbase

! Set second corner
         norm = sqrt(xishift**2 + yishift**2 + zishift**2)
         indmax = maxloc(norm)
         inodc(2) = igingn(indmax(1))
         icngin(2) = indmax(1)
      
! Set the rest of corners
         do inc = 3,nnodc
! Measure triangles
            measure = 0._kr
            do i = 1,nnodi
               indnodi = igingn(i)
               if (count(inodc.eq.indnodi).eq.0) then
                  meas = 0._kr
                  badangle = .false.
                  do j = 2,inc-1
                     indni = icngin(j)
                     scalarp = xishift(i)*xishift(indni) + yishift(i)*yishift(indni) + zishift(i)*zishift(indni)
                     cangle = scalarp/(norm(i) * norm(indni))
                     if (abs(cangle).gt.limcangle) badangle = .true.
                     meas = meas + scalarp
                  end do
                  if (badangle) meas = 0.0_kr
                  measure(i) = meas
               end if
            end do
            if (maxval(measure).eq.0.0_kr) then 
               stop 'Limit for badness probably too strict, algorithm for selecting corners fail.'
            end if
           
            indmax = maxloc(measure)
            inodc(inc) = igingn(indmax(1))
            icngin(inc) = indmax(1)
         end do

         deallocate(xishift,yishift,zishift,norm,measure)

      case default
         print *, 'Unknown selection...'
         goto 3
      end select

! Sort vector inodc
      do inc = nnodc,2,-1
         indmax = maxloc(inodc(1:inc))
         if (inc.ne.indmax(1)) then
            buffer = inodc(inc)
            inodc(inc) = inodc(indmax(1))
            inodc(indmax(1)) = buffer
         end if
      end do

! Check that all corner nodes are interface nodes
      do inc = 1,nnodc
         indnc = inodc(inc)
         if (count(igingn.eq.indnc).eq.0) then
            write(*,*) 'Warning: Node ',indnc,' selected as a corner is not on interface.'
         end if
      end do

! CN - list of corner nodes
      name = name1(1:lname1)//'.CN'
      open (unit=idcn,file=name,status='replace',form='formatted')
      write(idcn,*) nnodc
      write(idcn,'(8i15)') inodc
      close(idcn)

      print *, 'Created file with corner nodes ', trim(name)

      deallocate(x,y,z)
      deallocate(igingn)
      deallocate(inodc,icngin)
      deallocate(kinet)
      deallocate(kmynodes,kinodes)
      deallocate(inet,nnet,nndf)
      deallocate(iets)

end subroutine getcorners

!***********************************************************************
subroutine solassembly
!***********************************************************************
!     Subroutine for global solution assembly
!***********************************************************************
implicit none

! Local variables
integer ::             lsols
real(kr),allocatable :: sols(:)

integer :: idofn, i, isol, inod, isub, isols, point, ndofn, ie

! Import of basic geometry
! GMIS - basic mesh data - structure:
!  * INET(LINET) * NNET(LNNET) * NNDF(LNNDF) * XYF(LXYF) *
      name = name1(1:lname1)//'.GMIS'
      open (unit=idgmi,file=name,status='old',form='formatted')
      rewind idgmi
      linet = linet
      lnnet = nelem
      lnndf = nnod
      allocate(inet(linet),nnet(lnnet),nndf(lnndf))
! read fields INET, NNET, NNDF from file IDGMI
      read(idgmi,*) inet
      read(idgmi,*) nnet
      read(idgmi,*) nndf
      close(idgmi)

! ES - list of global element numbers in subdomains - structure:
      name = name1(1:lname1)//'.ES'
      open (unit=ides,file=name,status='old',form='formatted')
      liets = nelem
      allocate(iets(liets))
      read(ides,*) iets
      close(ides)

! Creation of field KINET(NELEM) with addresses before first global node of element IE in field inet
      lkinet = nelem
      allocate(kinet(lkinet))
      kinet(1) = 0
      do ie = 2,nelem
         kinet(ie) = kinet(ie-1) + nnet(ie-1)
      end do

! Creation of field KDOF(NNOD) with addresses before first global
! dof of node
      lkdof = nnod
      allocate(kdof(lkdof))
      kdof(1) = 0
      do inod = 2,nnod
         kdof(inod) = kdof(inod-1) + nndf(inod-1)
      end do

! Prepare space for global solution 
      lsol = ndof
      allocate(sol(lsol))
      sol = 0._kr

! Assembly solution
      write(*,*) 'Assembling global solution from subdomain files...'
      lkmynodes = nnod
      allocate(kmynodes(lkmynodes))
      do isub = 1,nsub
         call bddc_getfname(name1,lname1,isub,'SOLS',fname)
         open(unit=idsols, file=fname, status='old', form='unformatted')
! Begin loop over subdomains
! Creation of field kmynodes(nnod) - local usage of global nodes
! if global node is in my subdomain, assigned 1, else remains 0
         call pp_mark_sub_nodes(isub,nelem,iets,liets,inet,linet,nnet,lnnet,&
                                kinet,lkinet,kmynodes,lkmynodes)

         lsols = 0
         do inod = 1,nnod
            if (kmynodes(inod).eq.1) then
               ndofn = nndf(inod)
               lsols = lsols + ndofn
            end if
         end do
         
         allocate(sols(lsols))

         read(idsols) (sols(i),i = 1,lsols)

         isols = 0
         do inod = 1,nnod
            if (kmynodes(inod).eq.1) then
               ndofn = nndf(inod)
               point = kdof(inod)
               do idofn = 1,ndofn
                  isols = isols + 1
                  isol = point + idofn
                  sol(isol) = sol(isol) + sols(isols)
               end do
            end if
         end do

         deallocate(sols)

! Close SOLS files
         close(idsols)
      end do
      
! SOL - global solution
      name = name1(1:lname1)//'.SOL'
      open(unit=idsol,file=name,status='replace',form='unformatted')
      rewind idsol
      
      write(idsol) (sol(i), i = 1,lsol)
! Nonsense writing in place where postprocessor expect reactions
      write(idsol) (sol(i), i = 1,lsol), 0.0D0, 0.0D0, 0.0D0
      write(*,*) 'Solution has been written into file ',name1(1:lname1),'.SOL'
      write(*,*) 'Warning: At the moment solver does not ',&
                 'resolve reaction forces. Record of these does not',&
                 ' make sense and is present only to make ',&
                 'postprocessor str3 happy with input data.'

      write(*,*) '...done'
      close(idsol)

      deallocate(kmynodes)
      deallocate(sol)
      deallocate(kdof)
      deallocate(kinet)
      deallocate(iets)
      deallocate(inet,nnet,nndf)

end subroutine solassembly
         
!***********************************************************************
subroutine graphpmd
!***********************************************************************
!     Generate weighted graph for the mesh for METIS
!***********************************************************************
implicit none

integer:: ie, nne, ine, indinet, indnode, netnx, indelmn, nelmn, nedges, &
          nadjelmx, nadje, nadjelm, iadjelm, indadjelm, nadjelmadj, &
          iadjelmadj, indadjelmadj, j, nnetx, lor, liaetetx, &
          ionerow
integer,allocatable:: ietn(:,:), netn(:), &
                      iaetet(:,:), naetet(:), &
                      edgeweights(:,:), &
                      onerow(:), onerowweig(:)

! Set limit for considering neighbours
! 1 - all adjacencies included
! 2 - corner adjacencies excluded
! n - do not consider adjacencies with number of shared nodes lower or equal to n
integer :: neighbouring = 1
 
logical :: match

      write (*,'(a,$)') 'Minimal number of shared nodes to call elements adjacent: '
      read (*,*) neighbouring

      write (*,'(a,$)') 'Input PMD mesh data ... '
! prepare arrays for mesh data
      allocate(inet(linet),nnet(nelem),nndf(nnod),netn(nnod))

! GMIF - basic mesh data - structure:
!  * INETF(LINETF) * NNETF(LNNETF) * NNDFF(LNNDFF) * XYF(LXYF) *
      name = name1(1:lname1)//'.GMIS'
      open (unit=idgmi,file=name,status='old',form='formatted')
! Import of basic fields
! read fields INET, NNET, NNDF from file IDGMI
      rewind idgmi
      read(idgmi,*) inet
      read(idgmi,*) nnet
      read(idgmi,*) nndf

      write (*,'(a)') 'done.'

      nnetx = maxval(nnet)

      write (*,'(a,$)') 'Create list of elements at a node ... '
! Count elements touching particular node
! zero whole field of counts
      netn = 0
      indinet = 0
      do ie = 1,nelem
         nne = nnet(ie)
         do ine = 1,nne
            indinet = indinet + 1
            indnode = inet(indinet)
            netn(indnode) = netn(indnode) + 1
         end do
      end do

! Allocate proper sizes of two-dimensional field for list of touching elements
      netnx = maxval(netn)
      liaetetx = nnetx*netnx
      allocate(ietn(nnod,netnx))
      allocate(onerow(liaetetx),onerowweig(liaetetx))

! Create list of elements touching particular node
! zero whole field of counts
      netn = 0
      ietn = 0

      indinet = 0
      do ie = 1,nelem
         nne = nnet(ie)
         do ine = 1,nne
            indinet = indinet + 1
            indnode = inet(indinet)
            netn(indnode) = netn(indnode) + 1
            indelmn = netn(indnode)
            ietn(indnode,indelmn) = ie
         end do
      end do
      write (*,'(a)') 'done.'

      write (*,'(a,$)') 'Count list of neigbours of elements ... '
! Count elements adjacent to element
      allocate(naetet(nelem))
      naetet = 0

      indinet = 0
      do ie = 1,nelem
         nne = nnet(ie)
         onerow     = 0
         onerowweig = 0
         ionerow = 0
         do ine = 1,nne
            indinet = indinet + 1
            indnode = inet(indinet)
            nelmn = netn(indnode)
            onerow(ionerow+1:ionerow + nelmn) = ietn(indnode,1:nelmn)
            ionerow = ionerow + nelmn
         end do
         ! parse onerow
         call parse_onerow(ie,neighbouring,onerow,onerowweig,lor)
         naetet(ie) = count(onerowweig.ge.neighbouring)
      end do
      nedges = sum(naetet) / 2

      write (*,'(a)') 'done.'

      write (*,'(a,$)') 'Create list of neigbours of elements ... '
! Allocate proper field for elements
      nadjelmx = maxval(naetet)
      allocate(iaetet(nelem,nadjelmx),edgeweights(nelem,nadjelmx))

! Create list of neighbouring elements
      indinet = 0
      do ie = 1,nelem
         onerow     = 0
         onerowweig = 0
         ionerow = 0
         nne = nnet(ie)
         do ine = 1,nne
            indinet = indinet + 1
            indnode = inet(indinet)
            nelmn = netn(indnode)
            onerow(ionerow+1:ionerow + nelmn) = ietn(indnode,1:nelmn)
            ionerow = ionerow + nelmn
         end do
         ! parse onerow
         call parse_onerow(ie,neighbouring,onerow,onerowweig,lor)
         ! now only adjacencies above the level considered
         iaetet(ie,1:lor) = onerow(1:lor)
         edgeweights(ie,1:lor) = onerowweig(1:lor)

      end do
      write (*,'(a)') 'done.'

      write (*,'(a,$)') 'Check graph ... '
! Check the graph
      do ie = 1,nelem
         nadjelm = naetet(ie)
         do iadjelm = 1,nadjelm
            indadjelm = iaetet(ie,iadjelm)
            match = .false. 
            indadjelmadj = 0
            nadjelmadj = naetet(indadjelm)
            do iadjelmadj = 1,nadjelmadj
               if (iaetet(indadjelm,iadjelmadj).eq.ie) then
                  if (match) then
                     write(*,*) 'Element ',ie,' multiply mentioned in the list of neighbours of element', &
                                indadjelm,'!'
                     stop
                  end if
                  match = .true.
                  indadjelmadj = iadjelmadj
               end if
            end do
            if (.not.match) then
               write(*,*) 'No match! Couple ',ie,indadjelm,' mentioned but not the opposite!'
               stop
            end if
            if (edgeweights(indadjelm,indadjelmadj).ne.edgeweights(ie,iadjelm)) then
               write(*,*) 'Non-matching edge weights between elements ', ie, indadjelm,'!'
!               write(*,*) 'Indices of adjacent elements:'
!               do ie2 = 1,nelem
!                  nadjelm = naetet(ie2)
!                  write(*,'(100i7)') ie2, iaetet(ie2,1:nadjelm)
!               end do
!               write(*,*) 'Weights of elements:'
!               do ie2 = 1,nelem
!                  nadjelm = naetet(ie2)
!                  write(*,'(100i7)') ie2, edgeweights(ie2,1:nadjelm)
!               end do
               stop
            end if
         end do
      end do
      write (*,'(a)') 'done.'

      write (*,'(a,$)') 'Export graph to file ... '
! GRAPH - list of neighbours for elements
      name = name1(1:lname1)//'.GRAPH'
      open (unit=idgraph,file=name,status='replace',form='formatted')
      rewind idgraph

! Write the list into file for Chaco
      write(idgraph,'(x,i10,x,i10,a)') nelem, nedges, '  1'
      do ie = 1,nelem
         nadje = naetet(ie)
         write(idgraph,'(200i9)') (iaetet(ie,j), edgeweights(ie,j), j = 1,nadje)
      end do
      write (*,'(a)') 'done.'
      write(*,'(a,a)') 'Graph exported into file ',trim(name)

      deallocate(inet,nnet,nndf,netn,ietn,iaetet,edgeweights,naetet,onerow,onerowweig)

end subroutine graphpmd

subroutine parse_onerow(ie,neighbouring,onerow,onerowweig,lor)
implicit none

integer,intent(in) :: ie
integer,intent(in) :: neighbouring
integer,intent(inout) :: onerow(:)
integer,intent(inout) :: onerowweig(:)
integer,intent(inout) :: lor

! local variables
integer :: io, indel

         ! eliminate myself
         where(onerow .eq. ie) onerow = 0
         ! eliminate multiplicities
         lor = count(onerow .ne. 0)
         onerow = pack(onerow,onerow.ne.0)
         onerow(lor+1:) = 0
         io = 0
         do 
            io = io + 1
            indel = onerow(io)
            if (indel.eq.0) exit
            onerowweig(io) = count(onerow.eq.indel)
            where (onerow(io+1:).eq.indel) onerow(io+1:) = 0
            lor = count(onerow .ne. 0)
            onerow = pack(onerow,onerow.ne.0)
            onerow(lor+1:) = 0
         end do
         ! no repeating indices in onerow
         ! eliminate down limit for adjacencies
         where(onerowweig(1:lor).lt.neighbouring) onerow(1:lor)     = 0
         where(onerowweig(1:lor).lt.neighbouring) onerowweig(1:lor) = 0
         lor = count(onerow .ne. 0)
         onerow     = pack(onerow,onerow.ne.0)
         onerowweig = pack(onerowweig,onerowweig.ne.0)
         onerow(lor+1:)     = 0
         onerowweig(lor+1:) = 0

end subroutine
      

!***********************************************************************
subroutine meshdivide
!***********************************************************************
!     Divide the mesh into subdomains using METIS
!***********************************************************************
use module_parsing
implicit none

integer*4:: iaux, indvertex, indadjncy, ivertex
integer*4:: wgtflag = 1, numflag = 1, edgecut, nedge, nvertex
integer*4::            ladjncy,   ladjwgt,   lxadj,   lvwgt,   loptions,   liets
integer*4,allocatable:: adjncy(:), adjwgt(:), xadj(:), vwgt(:), options(:), iets(:)
 
! GRAPH - list of neighbours for elements
      name = name1(1:lname1)//'.GRAPH'
      open (unit=idgraph,file=name,status='old',form='formatted')
      rewind idgraph

! read initial line in the file with graph
      call rdline(idgraph)
      read(line,*) nvertex, nedge, iaux

      ladjncy = 2*nedge
      allocate(adjncy(ladjncy))
      ladjwgt = ladjncy
      allocate(adjwgt(ladjwgt))
      lxadj   = nvertex + 1
      allocate(xadj(lxadj))

      indadjncy = 0
      xadj(1)   = 1
      do ivertex = 1,nvertex
         call rdline(idgraph)
         indvertex = 0
         do
            indadjncy = indadjncy + 1
            indvertex = indvertex + 1
            call getstring
            read(string,*) adjncy(indadjncy)
            call getstring
            read(string,*) adjwgt(indadjncy)
            if (kstring.eq.0) then
               exit
            end if
         end do
         xadj(ivertex+1) = xadj(ivertex) + indvertex
      end do
      close(idgraph)
      write(*,'(a,a,a,i8,a)') 'File ',trim(name), ' imported. It has', fileline,' lines.'

      loptions = 5
      allocate(options(loptions))
      options = 0
      lvwgt   = 1
      allocate(vwgt(lvwgt))
      liets = nelem
      allocate(iets(liets))

      write(*,'(a,i6,a)') 'Calling METIS to divide into ',nsub,' subdomains...'
      if (nsub.le.8) then
         call METIS_PartGraphRecursive(nvertex,xadj,adjncy,vwgt,adjwgt,wgtflag,numflag,nsub,options,edgecut,iets)
      else
         call METIS_PartGraphKWay(nvertex,xadj,adjncy,vwgt,adjwgt,wgtflag,numflag,nsub,options,edgecut,iets)
      end if
      write(*,'(a,i9)') 'resulting number of cut edges:',edgecut

! ES - list of global element numbers in subdomains - structure:
      name = name1(1:lname1)//'.ES'
      open (unit=ides,file=name,status='replace',form='formatted')
      write(ides,'(1i5)') iets
      close(ides)

      write(*,'(a,a)') 'Division exported into file ',trim(name)

      deallocate(vwgt)
      deallocate(iets)
      deallocate(options)
      deallocate(adjncy)
      deallocate(adjwgt)
      deallocate(xadj)

      return

end subroutine meshdivide

!***********************************************************************
subroutine pp_mark_sub_nodes(isub,nelem,iets,liets,inet,linet,nnet,lnnet,&
                             kinet,lkinet,knodes,lknodes)
!***********************************************************************
!     Marks nodes of subdomain ISUB in array KNODES
!***********************************************************************
implicit none
      
! subdomain index
integer,intent(in):: isub 

! global number of elements
integer,intent(in):: nelem 

! indices of elements in subdomains
integer,intent(in):: liets 
integer,intent(in)::  iets(liets) 

! indices of nodes on elements
integer,intent(in):: linet 
integer,intent(in)::  inet(linet) 

! number of nodes on elements
integer,intent(in):: lnnet 
integer,intent(in)::  nnet(lnnet) 

! key to INET - pointers before element nodes
integer,intent(in):: lkinet 
integer,intent(in)::  kinet(lkinet) 

! key nodes - marked nodes for subdomain ISUB
integer,intent(in)::  lknodes 
integer,intent(out)::  knodes(lknodes) 

! local variables
integer:: ie, indsub, ine, nne, ipoint, indng

      knodes = 0
      do ie = 1,nelem
         indsub = iets(ie)
         if (indsub.eq.isub) then
            nne = nnet(ie)
            ipoint = kinet(ie)
            do ine = 1,nne
               indng = inet(ipoint + ine)
               knodes(indng) = 1
            end do
         end if
      end do

end subroutine pp_mark_sub_nodes

!*****************
subroutine es0toes
!*****************
implicit none

! Local variables
integer :: ie, ios
integer ::            liaux
integer,allocatable :: iaux(:)

      write(*,*) 'Trying to create .ES from .ES0 ...'
      name = name1(1:lname1)//'.ES0'
      open (unit=ides0,file=name,status='old',form='formatted', iostat=ios)
      if( ios/=0) then
         write(*,*) 'Error opening file .ES0: ', ios
         stop
      else
         name = name1(1:lname1)//'.ES'
         open (unit=ides,file=name,status='replace',form='formatted')
         ! in a loop over number of elements, read the old list and increase it by 1
         liaux = nelem
         allocate (iaux(liaux))
         read(ides0,*) iaux
         do ie = 1, nelem
            write(ides,*) iaux(ie) + 1
         end do
         deallocate (iaux)
!         do ie = 1, nelem
!            read(ides0,*) indsub
!            write(ides,*) indsub + 1
!         end do
         close(ides0)
         close(ides)
         open (unit=ides,file=name,status='old',form='formatted')
      end if

      return
end subroutine es0toes

!*******************************************************
subroutine bddc_getfname(name1,lname1,isub,suffix,fname)
!*******************************************************
!     Prepares name of file for subdomain ISUB with structure:
!     name of problem / number of subdomain + . + suffix

      implicit none
      
! subdomain number
      integer,intent(in) :: isub

! name of the problem
      integer,intent(in) ::     lname1
      character(*),intent(in) :: name1 

! suffix of the generated file
      character(*),intent(in) :: suffix

! generated name
      character(*),intent(out) :: fname

! local variables
      integer lsuffix, lfname, zerostart, zeroend, i

      fname = ' '
      lsuffix = len(suffix)
      lfname = len(fname)
      fname(1:lname1) = name1(1:lname1)
      fname(lname1+1:lname1+1) = '/'
      zerostart = lname1+2
      zeroend = lfname-lsuffix-1
      do i = zerostart,zeroend
         fname(i:i) = '0'
      end do
      fname(zeroend+1:zeroend+1) = '.'
      fname(zeroend+2:zeroend+lsuffix+1) = suffix

      if (isub.lt.10) then
         write(fname(zeroend:zeroend),'(i1)') isub
      else if (isub.lt.100) then
         write(fname(zeroend-1:zeroend),'(i2)') isub
      else if (isub.lt.1000) then
         write(fname(zeroend-2:zeroend),'(i3)') isub
      else if (isub.lt.10000) then
         write(fname(zeroend-3:zeroend),'(i4)') isub
      else
         write(*,*) 'isub = ',isub,': Out of range for file name!'
      end if

end subroutine

end
