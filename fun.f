      function fun(n,x)
      use phase
      use constant
      implicit none
      integer n,i,j
      real x(n)
      real :: fun
      real*8 chisquare
      real*8 xd(n)
      do i=1,n
        xd(i)=xbound(1,i)+dble(x(i))
     +   *(xbound(2,i)-xbound(1,i)) 
      end do
      chisquare=0.0d0
      call calphase(xd,cutoff,n)
      do j=1,phasenum
        do i=1,melabnum
        chisquare=(phasetheory(i,j+1)-phaseexp(i,j+1))**2+chisquare
        end do
      end do
      chisquare=chisquare/(phasenum*melabnum)
      chisquare=dsqrt(chisquare)
c the fortran program can only be used in maxum problems      
      chisquare=1.0d0/(chisquare-0.5d0)
      fun=real(chisquare)
      RETURN
      END
c this subroutine readin the standard file (like SP07-100.txt)
c and save the experemental data in module phase
