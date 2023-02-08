      program ga
c======================================================================
c     Sample driver program for pikaia.f
c======================================================================
      use phase
      use constant      
      implicit none
      integer n, seed, i, status
      parameter (n=24)
      real ctrl(12), x(n), f, fun
c     note this xs' range ([0,1])  
      external fun
c
c        (twod is an example fitness function, a smooth 2-d landscape)
c
c     First, initialize the random-number generator
c
      write(*,'(/A$)') ' Random number seed (I*4)? '
      read(*,*) seed
      call rninit(seed)
c
c     Set control variables (use defaults)
      do 10 i=1,12
         ctrl(i) = -1
   10 continue
      ctrl(2)=1
      
      call ini_phase
      call constant_int

c     Now call pikaia
      call pikaia(fun,n,ctrl,x,f,status)
c
c     Print the results
      write(*,*) ' status: ',status
      write(*,*) '      x: ',x
      write(*,*) '      f: ',100.0d0-f

      call releaseval
      call constant_rel
      stop
      end