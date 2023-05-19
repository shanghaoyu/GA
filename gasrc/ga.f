      program main
c======================================================================
c     Sample driver program for pikaia.f
c======================================================================
      use phase
      use constant     
      use omp_lib
      implicit none
      integer n, i, status
      real ctrl(12), f, fun
      real,allocatable :: x(:)
      
c     note this xs' range ([0,1])  
      external fun
      external pikaia
c
c        (twod is an example fitness function, a smooth 2-d landscape)
c
c     First, initialize the random-number generator
c
      call omp_set_num_threads(80)
c
c     Set control variables (use defaults)
      do 10 i=1,12
         ctrl(i) = -1
   10 continue
      ctrl(1)=2000
      ctrl(2)=1000
      ctrl(10)=1
      
      call ini_phase
      call constant_int
      n=num_par
      allocate(x(num_par))
      call rninit()
c     Now call pikaia
      call pikaia(fun,n,ctrl,x,f,status)
c
c     Print the results
      write(*,*) ' status: ',status
      do i=1,n
      write(*,"(f10.6)") x(i)
      end do
      write(*,*) '      f: ',0.5+1.0/f

      call releaseval
      call constant_rel
      stop
      end