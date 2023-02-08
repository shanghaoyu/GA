        
      module phase
c this module include globle varables will be used in 
c all the subroutines,the initial subroutine set the 
c initial values and allocate the varables.
        real*8,save,allocatable :: phaseexp(:,:),phasetheory(:,:)
        integer,save :: phasenum,melabnum
        integer,save :: phasename(40,2)
        contains
        subroutine ini_phase
c this subroutine need a standard input file "expdata.txt"
c read in the phasenum phasename,allocate phaseexp phasetheory
c read in phaseexp

c this subroutine should be called once and only once
            integer  :: i,a,b
            open (unit=12,file='expdata.d')
            read (12,*)
            read (12,*)
            i=0
            a=0
            b=0
            phasename=0
            do while(a >= 0)
              read(12,*) a,b
              i=i+1
              phasename(i,1)=a
              phasename(i,2)=b
            end do
            phasenum=i-1
            allocate(phaseexp(41,phasenum+1))
            allocate(phasetheory(41,phasenum+1))
            phaseexp=0.0d0
            phasetheory=0.0d0
            read(12,*)
            do i=1,41
                read (12,*,iostat=ios) phaseexp(i,:)
                if (ios /= 0) exit
            end do
            close(12)
            melabnum=i-1
        return          
        end subroutine

        subroutine calphase(contactvalue,cutoff,n)
c array contact contains the constants in contact terms        
c n is the dimension
c cutoff is the cutoff of the nuclear force
c this subroutine will use two files 'input.d'
            integer :: n 
            real*8 :: contactvalue(n),cutoff
            open(UNIT=10,file='input.d') 
            write(10,*) phasenum
            do i=1,phasenum
            write(10,*) phasename(i,:)
            end do
            write(10,*) melabnum
            do i=1,melabnum
            write(10,*) phaseexp(i,1)
            end do
1002        format(f10.6)
            do i=1,n
            write(10,1002) contactvalue(i)
            end do
            write(10,*) cutoff
            close(10)

c  calculate theoretic phases

            CALL system('phase.exe')
     
c  read in the theoretic phases
            open(UNIT=11,file='output.d')
            do i=1,melabnum
            read(11,*) phasetheory(i,:)   
            end do
            close(11)
        return
        end subroutine

        subroutine  releaseval
            deallocate(phaseexp)
            deallocate(phasetheory)
        return
        end subroutine

      end module