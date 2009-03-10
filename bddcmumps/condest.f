C***********************************************************************
      subroutine condest(nw,nwx,w,wchol,lw,x,y, cond)
C***********************************************************************
C Subroutine for calculation of condition number to matrix W
C***********************************************************************

      implicit none

      integer*4 nw, nwx, lw, i, j, maxit

      real*8 tol, eigmax, eigmin, cond, pival, shift
      real*8 w(lw), wchol(lw), x(nwx), y(nwx)
      logical singularity

C***********************************************************************
C User specified parameters of inverse iteration and power method
      tol   = 1.0E-9
      maxit = 1000
C Minimal pivot for Choleski factorization
      pival = 1.0E-9
C***********************************************************************
      
      if (nw.eq.0) then
         cond = 1.0
         goto 77
      end if

C Find largest eigenvalue of matrix W
      call powmet(nwx,nw,w,x,y,tol,maxit, eigmax)
      write(*,*) 'eigmax:',eigmax

C Find lowest eigenvalue of matrix W      
      shift = 1.0D0
      call invit(nwx,nw,w,wchol,x,y,tol,maxit,pival,shift, 
     *           eigmin,singularity)
      if (singularity) then
         write(*,*) 'Shift',shift,' is eigenvalue!'
         cond = -1.0D0
         goto 77
      end if

C Print the matrix into file
      open(unit=99,file='matice.txt',form='formatted')
      do i = 1,nw
         write(99,1001) (wchol((i-1)*nwx+j),j = 1,nw)
      end do
      close(99)
 1001 format(100f12.5)

      write(*,*) 'eigmin:',eigmin

C Evaluate condition number to matrix W
      cond = eigmax/eigmin
      
 77   return
      end

C***********************************************************************
      subroutine invit(ng,n,a,achol,x,y,tol,maxit,pival,shift, 
     *                 eigval,singularity)
C***********************************************************************
C Subroutine for inverse iteration power method
C Find the eigenvalue of matrix A that is closest to prescribed SHIFT
C modified for tridiagonal matrices 
C Jakub Sistek 12.2.2007
C***********************************************************************
      implicit none

      integer*4 ng,n,maxit, i, j, iter

      real*8 tol, shift, eigval, normy, theta, error, pival, sum
      real*8 a(ng*ng), achol(ng*ng), x(ng), y(ng)

      logical singularity

C Initial guess
      call rzero(y,n)
      y(1) = 1.0D0

C Apply shift
      do i = 1,n
         a((i-1)*ng+i) = a((i-1)*ng+i) - shift
      end do

C Print the matrix
C      do i = 1,n
C         write(*,1001) (a((i-1)*ng+j),j = 1,n)
C      end do
C 1001 format(100f12.5)
      
C Decomposition of matrix A by Choleski : (A-shift*I) = LL^T
      call chol(ng,n,pival,a,achol,singularity)
      if(singularity) return
      
C Loop over iterations of inverse iteration 
      do 1000 iter = 1,maxit
C        write(*,*) 'iter = ',iter

C ||y||_2
         normy = 0.0D0
         do i = 1,n
            normy = normy + y(i)*y(i)
         end do
         normy = sqrt(normy)

C Normalising of vector Y: y = y/||y||_2
         do i = 1,n
            y(i) = y(i)/normy
         end do

C  x = y
         do i = 1,n
            x(i) = y(i)
         end do

C Solving y = (A-shift*I)^-1 x
C only tridiagonal matrix
c         call rzero(y,n)
C First backsubstitution
         do i = 1,n
            sum = 0.0D0
            do j = 1,i-1
               sum = sum + achol((i-1)*ng + j)*y(j)
            end do
            y(i) = (y(i) - sum)/achol((i-1)*ng + i)
         end do
C Second backsubstitution
         do i = n,1,-1
            sum = 0.0D0
            do j = i+1,n
               sum = sum + achol((j-1)*ng + i)*y(j)
            end do
            y(i) = (y(i) - sum)/achol((i-1)*ng + i)
         end do

C theta = xT*y
         theta = 0.0D0
         do i = 1,n
            theta = theta + x(i)*y(i)
         end do

C error = ||y - theta*x||_2
         error = 0.0D0
         do i = 1,n
            error = error + (y(i)-theta*x(i))**2
         end do
         error = sqrt(error)

C Test prints
C         write(*,*) iter,eigvaln,nas,y

C Stop criterion
         if(error.lt.tol*abs(theta)) then
            write(*,'(a,i5)') ' number of inverse iterations:',iter
            goto 51
         end if
 1000 end do
      write(*,*) ' Number of inverse iterations exceeded limit!'

   51 continue

C Return eigenvalue
      eigval = shift + 1.0D0/theta
C eigenvector
      do i = 1,n
         y(i) = y(i)/theta
      end do

      return
      end

C***********************************************************************
      subroutine powmet(ng,n,a,x,y,tol,maxit, eigval)
C***********************************************************************
C Subroutine for power method
C Find the largest eigenvalue of matrix A modified for tridiag matrices 
C Jakub Sistek 12.2.2007
C***********************************************************************
      implicit none

      integer*4 ng,n,maxit, i, iter

      real*8 tol, eigval, normy, theta, error
      real*8 a(ng*ng), x(ng), y(ng)

C Initial guess
      call rzero(y,n)
      y(1) = 1.0D0

C Loop over iterations of inverse iteration 
      do 1000 iter = 1,maxit
C        write(*,*) 'iter = ',iter

C ||y||_2
         normy = 0.0D0
         do i = 1,n
            normy = normy + y(i)*y(i)
         end do
         normy = sqrt(normy)

C Normalising of vector Y: x = y/||y||_2
         do i = 1,n
            x(i) = y(i)/normy
         end do

C Multiplication of matrix A and vector X
C only tridiagonal matrix
         call rzero(y,n)
C Multiply diagonal
         do i = 1,n
            y(i) = y(i) + a((i-1)*ng + i)*x(i)
         end do
C Multiply subdiagonal
         do i = 2,n
            y(i) = y(i) + a((i-1)*ng + i - 1)*x(i-1)
         end do
C Multiply superdiagonal
         do i = 1,n-1
            y(i) = y(i) + a((i-1)*ng + i + 1)*x(i+1)
         end do

C theta = xT*y
         theta = 0.0D0
         do i = 1,n
            theta = theta + x(i)*y(i)
         end do

C error = ||y - theta*x||_2
         error = 0.0D0
         do i = 1,n
            error = error + (y(i)-theta*x(i))**2
         end do
         error = sqrt(error)

C Stop criterion
         if(error.lt.tol*abs(theta)) then
            write(*,'(a,i5)') ' number of iterations of power method:',
     *                        iter
            goto 51
         end if
 1000 end do
      write(*,*) ' Number of iterations of power method exceeded limit!'

   51 continue

C Return eigenvalue
      eigval = theta
C eigenvector
      do i = 1,n
         y(i) = y(i)
      end do

      return
      end

C********************************************************************
      SUBROUTINE CHOL(NG,N,PIVAL,B,D,SINGULARITY)
C********************************************************************
C Podprogram pro Choleskeho rozklad tridiagonalni matice B = D*DT   *
C input:  PIVAL - nejmensi pivot                                    *
C         NG - alokovany rozmer matic                               *
C         N - rozmer matic                                          *
C         B - matice pro rozklad tridiagonalni                      *
C output  D - dolni trojuhelnikova matice vznikla Choleskeho roz.   *
C         SINGULARITY - singularnost matice (logicka)               * 
C                                                                   *
C programoval: Jakub Sistek 4.9.2006                                *
C********************************************************************
C
      IMPLICIT NONE
      INTEGER NG, N, K
      REAL*8 PIVAL, BSDIAG
      REAL*8 B(NG*NG),D(NG*NG)

      LOGICAL SINGULARITY

      SINGULARITY = .FALSE.
C
      CALL RZERO(D,NG*NG)
C
      DO 1010 K = 1,N

C        UPLNY VYBER HLAVNIHO PRVKU JAKO MAXIMA Z PRVKU
C        K-TEHO SLOUPCE MATICE NA DIAGONALE
C         BM = B((K-1)*NG + K)
C         KM = K
C         DO L = K+1,N
C            IF (B((L-1)*NG + L).GT.BM) THEN
C               BM = B((L-1)*NG + L)
C               KM = L
C            END IF
C         END DO

C        PREHOZENI I-TEHO A KM-TEHO RADKU A SLOUPCE MATICE
C         IF(KM.NE.K) THEN
C            write(*,*) 'Pivotuji'
C           VYMENA RADKU
C            DO J = 1,N
C               BUFF = B((K-1)*NG + J)
C               B((K-1)*NG+J)  = B((KM-1)*NG+J)
C               B((KM-1)*NG+J) = BUFF
C            END DO
C           VYMENA SLOUPCE
C            DO J = 1,N
C               BUFF = B((J-1)*NG + K)
C               B((J-1)*NG+K)  = B((J-1)*NG+KM)
C               B((J-1)*NG+KM) = BUFF
C            END DO
C         END IF

         IF (K.GT.1) THEN
            BSDIAG = B((K-1)*NG + K) 
     *             - D((K-1)*NG + K-1)*D((K-1)*NG + K-1)
         ELSE
            BSDIAG = B((K-1)*NG + K) 
         END IF

C         WRITE(*,*) 'chol: BSDIAG = ',BSDIAG
         IF (BSDIAG.LT.PIVAL) THEN
            write(*,*) 'Choleski: Diagonal entry:',BSDIAG,',pivot',PIVAL
            WRITE(*,*) 'Choleski: Pivot less than limit!'
            SINGULARITY = .TRUE.
            RETURN
         END IF
         D((K-1)*NG + K) = SQRT(BSDIAG)

         D(K*NG + K) = B(K*NG + K) / D((K-1)*NG + K)
 1010 ENDDO
c     WRITE(*,*) 'Chol. O.K.'
C
      RETURN
      END

C***********************************************************************
      subroutine rzero(r,lr)
C***********************************************************************
C Subroutine which zero doble precision field R of lenght LR
C***********************************************************************

      implicit none

      integer*4 lr, i

      real*8 r(lr)

      do i = 1,lr
         r(i) = 0.0D0
      end do

      return
      end

