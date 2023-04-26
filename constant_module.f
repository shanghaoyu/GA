      module constant
        implicit none
        real*8,save  :: cutoff
        real*8,save, allocatable :: xbound(:,:)
        integer,save :: num_par 
        contains
          subroutine constant_int
          implicit none
          character(len=20) a
          integer n,i
          open(unit=10,file='constant.d')
          read(10,"(a6,2x,f10.4)") a,cutoff
          read(10,"(A10,2X,I5)") a,n
          num_par=n
          allocate(xbound(2,n))
          read(10,*)        
          do i=1,n
          read(10,*) xbound(:,i)
          end do
          close(10)
          end subroutine
          subroutine constant_rel
            deallocate(xbound)
          end subroutine
      end module