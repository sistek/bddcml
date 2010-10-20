!***********************************************************************
program ppbddcml
!***********************************************************************
! parallel preprocessor for BDDC in multilevel setting
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
! *.ES   - ASCII file with list of subd. numbers of elements, from 1 
! *.CN   - ASCII file with list of "corner nodes"
!
! programmed by Jakub Sistek            Bologna               14.10.2010 
!***********************************************************************

use module_graph
use module_utils

implicit none
include "mpif.h"


! ######### PARAMETERS TO SET:
! precision
integer,parameter:: kr = kind(1.D0)
! debugging mode - a lot of output
logical,parameter :: debug = .true.
! profiling mode - timing info
logical,parameter :: profile = .true.
! limit on number of concurrently opened files on filesystem
integer,parameter :: num_files_lim = 400
! values assumed as zero 
real(kr),parameter :: numerical_zero = 1e-15_kr
! maximal number of entries in element matrix
integer,parameter :: lelmx = 7921
!######### END OF PARAMETERS TO SET


integer,parameter:: idpar = 1, idgmi = 2, ides = 4, &
                    idbase = 100, idelm = 10

!  parallel variables
integer :: myid, comm_all, comm_self, nproc, ierr 
integer :: stat(MPI_STATUS_SIZE)

integer:: ndim, nsub, nelem, ndof, nnod, nnodc, meshdim

integer ::           linet,   lnnet
integer,allocatable:: inet(:), nnet(:)
integer ::           lpart
integer,allocatable:: part(:)
integer ::           liets
integer,allocatable:: iets(:)


integer ::           lnnet_loc
integer,allocatable:: nnet_loc(:)
integer ::             linet_loc
integer*4,allocatable:: inet_loc(:)

integer ::           lpart_loc
integer,allocatable:: part_loc(:)

integer,parameter:: lproblemnamex = 100, lfilenamex = 130
character(lproblemnamex):: problemname
character(lfilenamex)   :: filename

integer:: nelem_locx, nelem_loc, el_start, el_finish
integer:: el_start_send, el_finish_send, i, ie, indinet, indproc
integer:: length_send, length
integer:: ie_loc, indel, indsub, nelem_sub

integer :: neighbouring, graphtype, edgecut

integer ::            lnelemsa
integer,allocatable :: nelemsa(:), nelemsaaux(:)

real(kr),allocatable :: elm(:)
integer :: idelms, iproc, isub, lelm, isub_loc, nelems
integer :: sub_start, sub_finish, nsub_locx, nsub_loc

integer ::            lkadjsub
integer,allocatable :: kadjsub(:)

integer ::            lnadj
integer,allocatable :: nadj(:)
integer ::            ladj
integer,allocatable :: adj(:)
integer :: point_kadjsub, indadj, nadjs


integer :: indinet_loc, indnnet_loc
integer :: ine, nne, indiets

! time variables
real(kr) :: t_part, t_element_files, t_neighbourings, t_total

      ! MPI initialization
!***************************************************************PARALLEL
      call MPI_INIT(ierr)
      ! Communicator
      comm_all  = MPI_COMM_WORLD
      comm_self = MPI_COMM_SELF
      call MPI_COMM_RANK(comm_all,myid,ierr)
      call MPI_COMM_SIZE(comm_all,nproc,ierr)
!***************************************************************PARALLEL

! Initial screen
      if (myid.eq.0) then
         write(*,'(a)') ' _____  _____  _____  ____   ____    ____  __   __ _      '
         write(*,'(a)') '|  _  \|  _  \|  _  \|  _  \|  _  \ / __ \|  \ /  | |     '
         write(*,'(a)') '| |_|  | |_|  | |_|  | | \  | | \  | /  \_|   V   | |     '
         write(*,'(a)') '|  ___/|  ___/|  ___/| |  | | |  | | |    | |\ /| | |     '
         write(*,'(a)') '| |    | |    |  _  \| |  | | |  | | |   _| | V | | |     '
         write(*,'(a)') '| |    | |    | |_|  | |_/  | |_/  | \__/ | |   | | |____ '
         write(*,'(a)') '|_|    |_|    |_____/|_____/|_____/ \____/|_|   |_|______|'
         write(*,'(a)') 'PPBDDCML - Parallel Pre/Post-processor for BDDC multilevel'
         write(*,'(a)') '=========================================================='


         ! Name of the problem
   10    write(*,'(a)') 'Name of the problem: '
         call flush(6)
         read(*,*) problemname
         if(problemname.eq.' ') goto 10

      end if
! Broadcast of name of the problem      
!***************************************************************PARALLEL
      call MPI_BCAST(problemname, lproblemnamex, MPI_CHARACTER, 0, comm_all, ierr)
!***************************************************************PARALLEL

      if (myid.eq.0) then
         write (*,'(a)') 'Minimal number of shared nodes to call elements adjacent: '
         call flush(6)
         read (*,*) neighbouring
      end if
! Broadcast basic properties of the problem
!***************************************************************PARALLEL
      call MPI_BCAST(neighbouring,1,MPI_INTEGER,      0, comm_all, ierr)
!***************************************************************PARALLEL

!*****************PROFILING
      if (profile) then
         call MPI_BARRIER(comm_all,ierr)
         call time_start
      end if
!*****************PROFILING

      if (myid.eq.0) then
         filename = trim(problemname)//'.PAR'
         open (unit=idpar,file=filename,status='old',form='formatted')
      end if

! Reading basic properties 
      if (myid.eq.0) then
         read(idpar,*) ndim, nsub, nelem, ndof, nnod, nnodc, linet
         read(idpar,*,err=18) meshdim
         goto 19
 18      meshdim = ndim
         write (*,*) 'Default dimension of the mesh: ',meshdim
 19      continue
         write(*,*)'Characteristics of the problem ',trim(problemname), ':'
         write(*,*)'  number of processors            nproc =',nproc
         write(*,*)'  number of dimensions             ndim =',ndim
         write(*,*)'  number of subdomains             nsub =',nsub
         write(*,*)'  number of elements global       nelem =',nelem
         write(*,*)'  number of DOF                    ndof =',ndof
         write(*,*)'  number of nodes global           nnod =',nnod
         write(*,*)'  number of constrained nodes     nnodc =',nnodc
         write(*,*)'  length of field INET            linet =',linet
         write (*,*) 'Mesh dimension read from the PAR file: ',meshdim
         if (meshdim .eq. 0) then
            meshdim = ndim
            write (*,*) 'Correction of mesh dimension to make sense: ',meshdim
         end if
         call flush(6)
      end if
! Broadcast basic properties of the problem
!***************************************************************PARALLEL
      call MPI_BCAST(ndim,     1,MPI_INTEGER,         0, comm_all, ierr)
      call MPI_BCAST(nsub,     1,MPI_INTEGER,         0, comm_all, ierr)
      call MPI_BCAST(nelem,    1,MPI_INTEGER,         0, comm_all, ierr)
      call MPI_BCAST(ndof,     1,MPI_INTEGER,         0, comm_all, ierr)
      call MPI_BCAST(nnod,     1,MPI_INTEGER,         0, comm_all, ierr)
      call MPI_BCAST(nnodc,    1,MPI_INTEGER,         0, comm_all, ierr)
      call MPI_BCAST(linet,    1,MPI_INTEGER,         0, comm_all, ierr)
      call MPI_BCAST(meshdim,  1,MPI_INTEGER,         0, comm_all, ierr)
!***************************************************************PARALLEL

!*****************PROFILING
      if (profile) then
         call MPI_BARRIER(comm_all,ierr)
         call time_start
      end if
!*****************PROFILING

! read mesh on root proces 
      if (myid.eq.0) then
         filename = trim(problemname)//'.GMIS'
         open (unit=idgmi,file=filename,status='old',form='formatted')
         rewind idgmi
         linet = linet
         lnnet = nelem
         allocate(inet(linet),nnet(lnnet))
! read fields INET and NNET from GMIS file
         read(idgmi,*) inet
         read(idgmi,*) nnet
         close(idgmi)
      end if

! determine initial element ranges for each subdomain 
      nelem_locx = (nelem + nproc - 1)/nproc
      nelem_loc  = min(nelem - myid*nelem_locx, nelem_locx)
      el_start   = (myid*nelem_locx) + 1
      el_finish  = el_start + nelem_loc - 1
      if (debug) then
         write(*,*) 'myid =',myid,'el_start = ',el_start,'el_finish = ',el_finish
         call flush(6)
      end if

! DISTRIBUTE THE MESH AMONG PROCESSORS
      ! fisrt distribute NNET array
      lnnet_loc = nelem_locx
      allocate(nnet_loc(lnnet_loc))
      call MPI_SCATTER(nnet,lnnet_loc,MPI_INTEGER,nnet_loc,lnnet_loc,MPI_INTEGER,0,comm_all,ierr)

      ! pad the nnet_loc at the last processor
      if (lnnet_loc .ne. nelem_loc) then
         do i = nelem_loc+1,lnnet_loc
            nnet_loc(i) = 0
         end do
      end if

      ! debug
      !write(*,*) 'myid =',myid,'nnet_loc = ',nnet_loc

      linet_loc = sum(nnet_loc)
      allocate(inet_loc(linet_loc))

      ! now distribute INET array
      if (myid.eq.0) then
         ! root process copies its data and distributes the inet array by messages
         length = sum(nnet_loc)
         inet_loc = inet(1:length)

         ! now send messages to others
         indinet = length + 1
         do indproc = 1,nproc - 1
            el_start_send  = indproc * nelem_locx + 1
            el_finish_send = min((indproc+1) * nelem_locx,nelem)
            length_send = 0
            do ie = el_start_send,el_finish_send
               length_send = length_send + nnet(ie)
            end do
            if (length_send.gt.0) then
               call MPI_SEND(inet(indinet),length_send,MPI_INTEGER,indproc,indproc,comm_all,ierr)
            end if
            indinet = indinet + length_send
         end do
         ! free memory on root
         deallocate(inet,nnet)
      else
         if (linet_loc.gt.0) then
            call MPI_RECV(inet_loc,linet_loc,MPI_INTEGER,0,myid,comm_all,stat,ierr)
         end if
      end if

      ! debug
      !write(*,*) 'myid =',myid,'inet_loc = ',inet_loc

      lpart_loc = nelem_locx
      allocate(part_loc(lpart_loc))

      ! prepare initial distribution of subdomains
      nelem_sub = (nelem + nsub - 1) / nsub
      do ie_loc = 1,nelem_loc
         indel  = el_start + ie_loc - 1
         indsub = int((indel-1) / nelem_sub) + 1

         part_loc(ie_loc) = indsub
      end do
      ! pad the nnet_loc at the last processor
      if (lpart_loc .ne. nelem_loc) then
         do i = nelem_loc+1,lpart_loc
            part_loc(i) = 0
         end do
      end if
         
      ! debug
      !write(*,*) 'myid =',myid,'part_loc = ',part_loc

! divide mesh into subdomains by ParMetis
      graphtype = 0 ! no weights
      call graph_pdivide_mesh(myid,nproc,comm_all,graphtype,neighbouring,nelem,nelem_loc,profile,&
                              inet_loc,linet_loc,nnet_loc,lnnet_loc,nsub,&
                              edgecut,part_loc,lpart_loc)
      if (myid.eq.0) then
         write(*,'(a,i9)') 'Mesh divided. Resulting number of cut edges:',edgecut
         call flush(6)
      end if


! prepare memory for the global array on root
      lpart = nproc*nelem_locx
      allocate(part(lpart))
! get the global part array
      call MPI_ALLGATHER(part_loc,lpart_loc,MPI_INTEGER,part,lpart_loc,MPI_INTEGER,comm_all,ierr)

! free some memory
      deallocate(part_loc)
      deallocate(inet_loc)
      deallocate(nnet_loc)
      
! create array IETS and copy proper part of PART array to it
      liets = nelem
      allocate(iets(liets))
      do i = 1,nelem
         iets(i) = part(i)
      end do
! free memory
      deallocate(part)

      if (myid.eq.0) then
         ! debug
         ! write(*,*) 'iets',iets

         ! write the obtained global array to disk
         filename = trim(problemname)//'.ES'
         open (unit=ides,file=filename,status='replace',form='formatted')
         write(ides,'(1i5)') iets
         close(ides)

         write(*,'(a,a)') 'Division exported into file ',trim(filename)
      end if

!*****************PROFILING
      if (profile) then
         call MPI_BARRIER(comm_all,ierr)
         call time_end(t_part)
      end if
!*****************PROFILING

! CREATE SUBDOMAIN FILES WITH ELEMENT MATRICES
!*****************PROFILING
      if (profile) then
         call MPI_BARRIER(comm_all,ierr)
         call time_start
      end if
!*****************PROFILING


      ! prepare memory for one element matrix
      allocate(elm(lelmx))

! Distribute subdomains
      nsub_locx = (nsub + nproc - 1)/nproc
      nsub_loc  = min(nsub - myid*nsub_locx, nsub_locx)
      sub_start   = (myid*nsub_locx) + 1
      sub_finish  = sub_start + nsub_loc - 1
      if (debug) then
         write(*,*) 'myid =',myid,'sub_start = ',sub_start,'sub_finish = ',sub_finish
         call flush(6)
      end if
      if (nsub_loc.gt.num_files_lim) then
         write(*,*) 'myid =',myid,':ERROR - Number of subdomain files greater than the limit.'
         call error_exit
      end if


      ! open ELM file on root processor
      if (myid.eq.0) then
         ! ELM - element stiffness matrices - structure:
         filename = trim(problemname)//'.ELM'
         open(unit=idelm,file=filename,status='old',form='unformatted')
         rewind idelm
      end if

      ! open subdomain files at each processor
      do isub = sub_start,sub_finish
         idelms = idbase + isub - sub_start + 1
         call getfname(problemname,isub,'ELM',filename)
         if (debug) then
            write(*,*) 'myid =',myid,': Creating file ',trim(filename)
         end if
         open(unit=idelms, file=filename, status='replace', form='unformatted')
      end do

      ! loop over elements
      do ie = 1,nelem
         isub = iets(ie)

         ! get processor taking care of this subdomain
         iproc = int((isub-1) / nsub_locx) 

         if (myid.eq.0) then
            ! root reads the matrix from file
            read(idelm) lelm,(elm(i),i = 1,lelm)
            if (iproc.eq.0) then
               ! if it is to be stored by 0, do not send messages, just store it
               idelms = idbase + isub - sub_start + 1
               write(idelms) lelm,(elm(i),i = 1,lelm)
            else
               ! send messages
               call MPI_SEND(lelm,1,MPI_INTEGER,iproc,isub,comm_all,ierr)
               call MPI_SEND(elm,lelm,MPI_DOUBLE_PRECISION,iproc,isub,comm_all,ierr)
            end if
         else 
            if (myid.eq.iproc) then
               call MPI_RECV(lelm,1,MPI_INTEGER,0,isub,comm_all,stat,ierr)
               call MPI_RECV(elm,lelm,MPI_DOUBLE_PRECISION,0,isub,comm_all,stat,ierr)

               idelms = idbase + isub - sub_start + 1
               write(idelms) lelm,(elm(i),i = 1,lelm)
            end if
         end if
      end do
      
      ! close files
      do isub = sub_start,sub_finish
         idelms = idbase + isub - sub_start + 1
         close(idelms)
      end do
      if (myid.eq.0) then
         close(idelm)
      end if

      ! free memory
      deallocate(elm)

!*****************PROFILING
      if (profile) then
         call MPI_BARRIER(comm_all,ierr)
         call time_end(t_element_files)
      end if
!*****************PROFILING

! FIND NEIGBHBOURINGS
!*****************PROFILING
      if (profile) then
         call MPI_BARRIER(comm_all,ierr)
         call time_start
      end if
!*****************PROFILING

! read mesh on root proces 
      if (myid.eq.0) then
         filename = trim(problemname)//'.GMIS'
         open (unit=idgmi,file=filename,status='old',form='formatted')
         rewind idgmi
      end if
      linet = linet
      lnnet = nelem
      allocate(inet(linet),nnet(lnnet))
! read fields INET and NNET from GMIS file
      if (myid.eq.0) then
         read(idgmi,*) inet
         read(idgmi,*) nnet
         close(idgmi)
      end if
!***************************************************************PARALLEL
      call MPI_BCAST(inet,  linet,MPI_INTEGER,         0, comm_all, ierr)
      call MPI_BCAST(nnet,  lnnet,MPI_INTEGER,         0, comm_all, ierr)
!***************************************************************PARALLEL

      lnelemsa = nsub ! number of elements in subdomains
      allocate(nelemsaaux(lnelemsa),nelemsa(lnelemsa))

      nelemsaaux = 0 ! number of elements in subdomains
      nelem_loc  = 0 ! number of elements on processors
      linet_loc  = 0 ! length of local linet
      do isub_loc = 1,nsub_loc
         isub = sub_start + isub_loc - 1

         ! pick my subdomains
         ! find numbers of elements in subdomains
         nelems = 0
         do ie = 1,nelem
            if (iets(ie).eq.isub) then
               nelems = nelems + 1
               linet_loc = linet_loc + nnet(ie)
            end if
         end do
         nelemsaaux(isub) = nelems

         nelem_loc = nelem_loc + nelems
      end do
!***************************************************************PARALLEL
      call MPI_ALLREDUCE(nelemsaaux, nelemsa, lnelemsa,MPI_INTEGER, MPI_SUM, comm_all, ierr)
!***************************************************************PARALLEL
      deallocate(nelemsaaux)
      ! debug
      !write(*,*) 'nelemsaaux',nelemsaaux
      !write(*,*) 'nelemsa',nelemsa

      ! get local INET and NNET
      lnnet_loc = nelem_loc
      allocate(nnet_loc(lnnet_loc),inet_loc(linet_loc))
      indnnet_loc = 0
      indinet_loc = 0
      do isub_loc = 1,nsub_loc
         isub = sub_start + isub_loc - 1

         ! pick my subdomains
         ! find numbers of elements in subdomains
         indinet     = 0
         do ie = 1,nelem
            nne = nnet(ie)
            if (iets(ie).eq.isub) then
               indnnet_loc = indnnet_loc + 1
               nnet_loc(indnnet_loc) = nne
               do ine = 1,nne
                  inet_loc(indinet_loc+ine) = inet(indinet+ine)
               end do
               indinet_loc = indinet_loc + nne
            end if
            indinet = indinet + nne
         end do
      end do
      ! create new IETS with elements ordered subdomain by subdomain
      indiets = 0
      do isub = 1,nsub
         do ie = 1,nelemsa(isub)
            indiets = indiets + 1
            iets(indiets) = isub
         end do
      end do
      deallocate(inet,nnet)

      ! debug
      !write(*,*) 'nnet_loc',nnet_loc
      !write(*,*) 'inet_loc',inet_loc

! find adjacent subdomains from parallel graph
      lkadjsub = nsub * nsub_loc
      allocate(kadjsub(lkadjsub))
      kadjsub = 0
      call graph_pget_sub_neighbours(myid,nproc,comm_all,1,nelem,nelem_loc,nsub,nsub_loc,sub_start,&
                                     inet_loc,linet_loc,nnet_loc,lnnet_loc, iets,liets,debug, kadjsub,lkadjsub)

      ! decode kadjsub
      lnadj = nsub_loc
      allocate(nadj(lnadj))
      do isub_loc = 1,nsub_loc

         point_kadjsub = (isub_loc-1)*nsub

         nadjs  = 0
         do i = 1,nsub
            if (kadjsub(point_kadjsub + i) .eq. 1) then
               nadjs = nadjs + 1
            end if
         end do
         nadj(isub_loc) = nadjs
      end do
      ladj = sum(nadj)
      allocate(adj(ladj))
      indadj = 0
      do isub_loc = 1,nsub_loc

         point_kadjsub = (isub_loc-1)*nsub

         do i = 1,nsub
            if (kadjsub(point_kadjsub + i) .eq. 1) then
               indadj = indadj + 1
               adj(indadj) = i
            end if
         end do
      end do
      ! check the length 
      if (ladj.ne.indadj) then
         write(*,*) 'Error: Length mismatch.'
         call error_exit
      end if
      deallocate(kadjsub)

      ! debug
      write(*,*) 'nadj',nadj
      write(*,*) 'adj',adj


      deallocate(adj)
      deallocate(nadj)
      deallocate(nnet_loc,inet_loc)

      deallocate(nelemsa)

!*****************PROFILING
      if (profile) then
         call MPI_BARRIER(comm_all,ierr)
         call time_end(t_neighbourings)
      end if
!*****************PROFILING

! free memory
      deallocate(iets)

!*****************PROFILING
      if (profile) then
         call MPI_BARRIER(comm_all,ierr)
         call time_end(t_total)

         ! Information about times
         if (myid.eq.0) then
            write(*,'(a)')         ' TIMES OF RUN:    '
            write(*,'(a,f11.3,a)') '  partitioning    ',t_part,          ' s'
            write(*,'(a,f11.3,a)') '  subdomain files ',t_element_files, ' s'
            write(*,'(a,f11.3,a)') '  neighbours      ',t_neighbourings, ' s'
            write(*,'(a)')         '  ___________________________'
            write(*,'(a,f11.3,a)') '  total           ',t_total,         ' s'
         end if
      end if
!*****************PROFILING

      ! MPI finalization
!***************************************************************PARALLEL
      call MPI_FINALIZE(ierr)
!***************************************************************PARALLEL

      end program

