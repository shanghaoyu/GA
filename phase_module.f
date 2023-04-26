        
c this version is for openmp 
      module phase
c this module include globle varables will be used in 
c all the subroutines,the initial subroutine set the 
c initial values and allocate the varables.
        use omp_lib
        real*8,save,allocatable :: phaseexp(:,:)
        integer,save :: phasenum,melabnum
        integer,save :: phasename(40,2),phasetheory(50,50)
!$omp threadprivate(phasetheory)
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
            integer :: n ,id
            character*2 :: filenum
            real*8 :: contactvalue(n),cutoff
            id=omp_get_thread_num()
            if(id < 10)then
            write(filenum,'(I1)') id 
            else
              write(filenum,'(I2)') id
            end if
            id=10+id
            open(unit=id,file='input'//filenum) 
            write(id,*) phasenum
            do i=1,phasenum
            write(id,*) phasename(i,:)
            end do
            write(id,*) melabnum
            do i=1,melabnum
            write(id,*) phaseexp(i,1)
            end do
1002        format(f10.6)
            do i=1,n
            write(id,1002) contactvalue(i)
            end do
            write(id,*) cutoff
            close(id)

c  calculate theoretic phases

            CALL system('./phase.o '//filenum)
     
c  read in the theoretic phases
            open(id,file='output'//filenum)
            do i=1,melabnum
            read(id,*) phasetheory(i,1:phasenum)   
            end do
            close(id)
        return
        end subroutine

        subroutine  releaseval
            deallocate(phaseexp)

        return
        end subroutine

      end module