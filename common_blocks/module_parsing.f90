module module_parsing
! Module for parsing ASCII files
! programmed by Jakub Sistek <sistek@vzlu.cz>, May 2007
      implicit none
      integer,parameter,private:: llinex = 2000, lstringx = 255
      character(llinex),private:: lineraw
      integer,private:: linerest = 0

! Variables passed to calling program
      character(llinex):: line      ! the rest of line from file
      character(lstringx):: string  ! stripped string
      integer:: fileline            ! number of line to be read
      integer:: lline               ! lenght of the rest of line till continuous non-blank characters
      integer:: lstring             ! lenght of string
      integer:: kstring             ! is there another string in the line?

contains

subroutine rdline(idfile)
!************************
!     Reads next non-comment line from unit IDFILE
!     Strips off leading blanks
!     Ignores everything after and including "#"
!     If ";" is found, behaves like a new line
!
!     LINE returns the line
!     LLINE returns the number of characters in non-blank portion
!     FILELINE is number of lines
!
!     If e.o.f. is reached, LINE returns 'EOF'
!     If read error occurs, LINE returns 'ERR'
      integer:: idfile, kexl, i, endline, lrest
      
      if (linerest.eq.0) then
 20      read (idfile,'(a)',end=80,err=90) lineraw
         fileline = fileline + 1
!---- skip comment line
         if (lineraw(1:1).eq.'#') goto 20
!---- change tabulators into spaces
         do i = 1,len(lineraw)
            if (lineraw(i:i).eq.char(9)) lineraw(i:i) = ' '
         end do
!---- skip blank line
         if (lineraw.eq.' ') goto 20
      end if
      endline = index(lineraw,';') - 1
      if (endline.lt.1) then 
         line = lineraw
         linerest = 0
      else
         line(1:endline) = lineraw(1:endline)
         do i = endline+1,len(line)
            line(i:i) = ' '
         end do
         lrest = len(lineraw)-endline-1
         lineraw(1:lrest) = lineraw(endline+2:len(lineraw))
         do i = lrest+1,len(lineraw)
            lineraw(i:i) = ' '
         end do
         linerest = 1
         if (lineraw.eq.' ') linerest = 0
      end if

!---- strip off leading blanks and do normal return after signif. line
      call stripline
      kexl = index(line(1:lline),'#')
      if (kexl.gt.1) lline = kexl-1
      return

   80 line = 'EOF '
      return
   90 line = 'ERR '
      call error('reading of line failed')
end subroutine rdline

subroutine stripline
!*******************
!     Strips leading blanks off line
!     and returns length of non-blank part.
!     changes LINE
!     changes LLINE

      integer:: n, k1, k2, k

      n = len(line)

!---- find last non-blank character
      do k2 = n, 1, -1
        if (line(k2:k2).ne.' ') goto 11
      end do
      k2 = 0
   11 continue

!---- find first non-blank character
      do k1 = 1, k2
        if (line(k1:k1).ne.' ') goto 21
      end do
   21 continue

!---- number of non-blank characters
      lline = k2 - k1 + 1
      if (lline.eq.0) return

!---- shift line so first character is non-blank
      line(1:lline) = line(k1:k2)

!---- pad tail of line with blanks
      do k = lline+1, n
        line(k:k) = ' '
      end do
end subroutine stripline

subroutine getstring
!*******************
!     Returns first string on line separately in STRING and its length LSTRING
!     moves next string in LINE to the first position in LINE, consequently changes LLINE
!     sets key KSTRING to 1 if another string left in LINE, sets KSTRING to 0 if LINE remains empty
      
      integer:: nextstring, i, llinenew, lspace, lbracket

!     find first blank space or bracket
      lspace     = index(line(1:lline),' ') - 1
      lbracket   = index(line(1:lline),'}') - 1
      if (lspace.eq.-1.and.lbracket.eq.-1) then
         lstring = lline
      else if (lspace.eq.-1) then
         lstring = lbracket
      else if (lbracket.eq.-1) then
         lstring = lspace
      else if (lbracket.gt.-1.and.lspace.gt.-1) then
         lstring = min(lbracket,lspace)
      else
         call error('Error in separating strings in',line)
      end if
      
      if (line(1:1).eq.'{') lstring = 1
      if (line(1:1).eq.'}') lstring = 1

!---- shift STRING so first character is non-blank
      string(1:lstring) = line(1:lstring)

!---- pad tail of STRING with blanks
      do i = lstring + 1,len(string) 
        string(i:i) = ' '
      end do

!---- shift LINE so first character is non-blank
      nextstring = 0
      do i = lstring + 1,lline
         if (line(i:i).ne.' ') then 
            nextstring = i
            exit
         end if
      end do

      if (nextstring.ne.0) then
!---- there is another string on the line so shift it to the first pos.
         llinenew = lline - nextstring + 1
         line(1:llinenew) = line(nextstring:lline)
!---- pad tail of LINE with blanks
         do i = llinenew + 1, len(line)
           line(i:i) = ' '
         end do
         lline = llinenew
         kstring = 1
      else
!---- there in no other string on the line
         kstring = 0
      end if

end subroutine getstring

subroutine lc2uc(input)
!**********************
! Routine for lower to upper case transformations
      character*(*) input
      integer:: n, i, k
      character(26) lcase, ucase
      data lcase / 'abcdefghijklmnopqrstuvwxyz' /
      data ucase / 'ABCDEFGHIJKLMNOPQRSTUVWXYZ' /
 
      n = len(input)
      do i=1, n
        k = index( lcase , input(i:i) )
        if (k.gt.0) input(i:i) = ucase(k:k)
      end do
end subroutine lc2uc

subroutine error(message,problem)
!********************************
! Terminates the code and prints the reason from MESSAGE and prints FILELINE
      character(*):: message
      character(*),optional:: problem
      if (present(problem)) then
         write(*,'(a,a,2x,a,a,i6)') 'Error in input file: ',message, problem,', line ',fileline
      else
         write(*,'(a,a,a,i6)') 'Error in input file: ',message,' , line ',fileline
      end if
      stop
end subroutine error

end module module_parsing

