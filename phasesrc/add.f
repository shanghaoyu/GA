c***********************************************************************
c    version one :March 2022
c    version two :June 2022 (change the representation of theparameters
c    in lsj representation)
c
c***********************************************************************
c    
c    in this code you can add any terms you want in the chiral force
c    using module addterms,note that you have to use mass of nucleon as
c    the unit of energy
c    
c***********************************************************************  
c
c    author: shang haoyu 
c            physics school
c            Peking university
c    email:  shy@stu.pku.edu.cn
c
c
c***********************************************************************   
c    module kqxy contains the variable we use mostly in this code

c    global varialbe xlab,yalb to deliver xmev and ymev,
c    mass is the mass of n p or np
c    pottype to deliver inn 

c    function initialize evaluate the variables we use in this code

c    function normk, normq means the length of vector k and q
c    thier variable z means cos(theta)
      module const 
         real*8 ::pi=3.141592653589793d0
         real*8 ::ga=1.29d0
         real*8 ::mpi=138.0390d0
         real*8 ::fpi=92.4d0
      end module
      module kqxy
        use const
        real*8,save :: xlab,ylab,mass,x,y,dwn,wnq,dwnq,x2,y2,c(9)
        integer,save :: pottype

        contains

          subroutine initialize 
            implicit real*8 (a-h,o-z)
            real*8 wnn(3)
            real*8 matrix(9,9)
            common /cnn/ inn
            common /cpot/ v(6),xmev,ymev
            real*8 conta(9)            
            common /con/ conta
            data wnn /938.272d0,938.9182d0,939.5653d0/
            data((matrix(j,i),i=1,9),j=1,9)/
     1 0.0198943d0,       0.0d0,       0.0d0,       0.0d0,        0.0d0,
     * 0.0596831d0,       0.0d0,       0.0d0,       0.0d0,
     2       0.0d0, 0.0099471d0,-0.0049735d0,-0.0149207d0, -0.0149207d0,
     *       0.0d0, 0.0298415d0,       0.0d0, -0.024867d0,
     3       0.0d0, 0.0397887d0, 0.0198943d0, 0.0596831d0,  0.0596831d0,
     *       0.0d0, 0.1193662d0,       0.0d0, 0.0994718d0,
     4 -0.019894d0,       0.0d0,       0.0d0,       0.0d0,        0.0d0,
     * 0.0198943d0,       0.0d0,       0.0d0,       0.0d0,
     5       0.0d0,-0.0099471d0,-0.0049735d0, 0.0149207d0,        0.0d0,
     *       0.0d0, 0.0099471d0, 0.0140674d0,-0.0099471d0,
     6       0.0d0, -0.039788d0, 0.0198943d0,-0.0596831d0,        0.0d0,
     *       0.0d0, 0.0397887d0, 0.0562697d0, 0.0397887d0,
     7       0.0d0,       0.0d0,-0.0397887d0,       0.0d0, -0.0596831d0,
     *       0.0d0,       0.0d0,       0.0d0, 0.0994718d0,
     8       0.0d0,       0.0d0, 0.0099471d0,       0.0d0, -0.0149207d0,
     *       0.0d0,       0.0d0,-0.0422023d0, 0.0049735d0,
     9       0.0d0,       0.0d0,-0.0397887d0,       0.0d0,  0.0596831d0,
     *       0.0d0,       0.0d0,-0.1688093d0,-0.0198943d0/
c   here we evaluate the number we use
                        
            xlab=xmev
            ylab=ymev
            pottype=inn
            mass=wnn(pottype)
            c=0.0d0
            do i=1,9
               do k=1,9
                   c(i)=c(i)+matrix(i,k)*conta(k)
               end do
            end do
            dwn=1.0d0/mass
            wnq=mass*mass
            dwnq=dwn*dwn

            x=xlab*dwn
            y=ylab*dwn
            x2=x*x
            y2=y*y

          end subroutine
        
          real*8 function normk(z)
            implicit real*8 (a-h,o-z)
            real*8 z
            normk=dsqrt(x*x+y*y+2.0d0*x*y*z)/2.0d0
            return
          end function
        
          real*8 function normq(z)
            implicit real*8 (a-h,o-z)
            real*8 z    
            normq=dsqrt(x*x+y*y-2.0d0*x*y*z)
            return
          end function
c   function w(q) in (20) PRC 66,014002(2002)  
         
          real*8 function wfunc(z)
             implicit real*8 (a-h,o-z)
             real*8 z    
             wfunc=dsqrt(4.0d0*(mpi*dwn)**2+normq(z)**2)
             return
          end function
c   function L(q) in (19) PRC 66,014002(2002)         
          real*8 function lfunc(z)
            implicit real*8 (a-h,o-z)
            real*8 z    
            lfunc=wfunc(z)/normq(z)*dlog((wfunc(z)+normq(z))
     1       /(2.0d0*mpi*dwn))
            return
          end function
         
        
        
      end module

c    in module add-terms we can add potential terms arbitrary
      
      module addterms
        use kqxy
        
c    we set some logical variable to control whether the terms 
c    are contained(private),true means they are contained
        logical :: vc=.true.
        logical :: vss=.true.
        logical :: vso=.true.
        logical :: vsigl=.true.
        logical :: vt=.true.
        logical :: vsigk=.true.

        private vc,vt,vss,vso,vsigl,vsigk
c    in the following terms w means terms with tau1*tau2
c    z is cos(theta) 

!!!! notice here we use mass of nucleon as the uint of energy

        contains
          real*8 function vcentral(z)
            real*8 z
            integer ::type=1
            if(vc) then
               select case(type)
               case(1)
                  vcentral=0.01d0*wnq*c(1)+(c(2)*wnq*wnq
     1 *1.0d-8*normq(z)**2+c(3)*wnq*wnq*1.0d-8*normk(z)**2)
               case(2)
            vcentral=0.01d0*wnq*c(1)/(1.0d0
     1      -(c(2)*wnq*wnq*1.0d-8*normq(z)**2
     2      +c(3)*wnq*wnq*1.0d-8*normk(z)**2)
     3      /0.01d0*wnq*c(1))
               case(3)
            vcentral=0.01d0*wnq*c(1)/(1.0d0
     1      -(c(2)*wnq*wnq*1.0d-8*normq(z)**2)
     2      /0.01d0*wnq*c(1))+c(3)*wnq*wnq
     3      *1.0d-8*normk(z)**2
               end select


            else
            vcentral=0.0d0
            end if
          return
          end function

          real*8 function vspinspin(z)
          real*8 z
          if(vss) then
          vspinspin=(0.01d0*wnq*c(4)+c(5)*wnq*wnq*1.0d-8*normq(z)**2
     1      +c(6)*wnq*wnq*1.0d-8*normk(z)**2) 
          else
          vspinspin=0.0d0
          end if
          return
          end function

          real*8 function vspinobit(z)
          real*8 z
          if(vso) then
          vspinobit=c(7)*wnq*wnq*1.0d-8
          else
          vspinobit=0.0d0
          end if
          return
          end function

          real*8 function vsigmaL(z)
          real*8 z
          if(vsigl) then
          vsigmaL=0.0d0
          else
          vsigmaL=0.0d0
          end if
          return
          end function

          real*8 function vtensor(z)
            real*8 z
            if(vt) then
            vtensor=c(8)*wnq*wnq*1.0d-8 
            else
            vtensor=0.0d0
            end if
          return
          end function

          real*8 function vsigmak(z)
          real*8 z
          if(vsigk) then
          vsigmak=c(9)*wnq*wnq*1.0d-8
          else
          vsigmak=0.0d0
          end if
          return
          end function
      end module   


c   we write another legendre 
        function legendre(x,j)
c
c
c   x is the independent variable
c
        real*8  x,legendrem1,a,b
        real*8  legendre
        integer j
c
c
c
c        compute legendre polynom for j equals zero
c
c
        if (j.gt.0) go to 1
        legendre=1.d0
        legendrem1=0.d0
        if (j.lt.0) legendre=0.d0
        return
c
c
c
c        compute legendre polynoms for j equals one
c
c
c
    1 legendre=x
        legendrem1=1.d0
        if (j.eq.1) return
c
c
c
c        compute legendre polynom for j greater or equal two
c
c
c
        do 2 i=2,j
        a=x*legendre
        b=a-legendrem1
        legendrem1=legendre
    2 legendre=-b/dfloat(i)+b+a
c

        return
        end

        function alj(func,l,j)
c   function alj calculate the integral a^j(l), l is obit-angular
c   momentum and j is total angular momentum,func is the function 
c   being integrated

c   ALJ^J(q',q)=pi*INTEGRATE[fuc(q',q,z) z^l P_J(z)，-1,1]

        implicit real*8 (a-h,o-z)
        real*8,external :: func
        real*8,external :: legendre
        external gset
        integer l,j
        real*8 ct(96),wt(96),pj(96),pi,x(273),a(273)
        real*8 alj
        

c   now we just set the ct and wt 
        data x(226)/0.999689503883230d0/, a(226)/0.000796792065552d0/
        data x(227)/0.998364375863181d0/, a(227)/0.001853960788946d0/
        data x(228)/0.995981842987209d0/, a(228)/0.002910731817934d0/
        data x(229)/0.992543900323762d0/, a(229)/0.003964554338444d0/
        data x(230)/0.988054126329623d0/, a(230)/0.005014202742927d0/
        data x(231)/0.982517263563014d0/, a(231)/0.006058545504235d0/
        data x(232)/0.975939174585136d0/, a(232)/0.007096470791153d0/
        data x(233)/0.968326828463264d0/, a(233)/0.008126876925698d0/
        data x(234)/0.959688291448742d0/, a(234)/0.009148671230783d0/
        data x(235)/0.950032717784437d0/, a(235)/0.010160770535008d0/
        data x(236)/0.939370339752755d0/, a(236)/0.011162102099838d0/
        data x(237)/0.927712456722308d0/, a(237)/0.012151604671088d0/
        data x(238)/0.915071423120898d0/, a(238)/0.013128229566961d0/
        data x(239)/0.901460635315852d0/, a(239)/0.014090941772314d0/
        data x(240)/0.886894517402420d0/, a(240)/0.015038721026994d0/
        data x(241)/0.871388505909296d0/, a(241)/0.015970562902562d0/
        data x(242)/0.854959033434601d0/, a(242)/0.016885479864245d0/
        data x(243)/0.837623511228187d0/, a(243)/0.017782502316045d0/
        data x(244)/0.819400310737931d0/, a(244)/0.018660679627411d0/
        data x(245)/0.800308744139140d0/, a(245)/0.019519081140145d0/
        data x(246)/0.780369043867433d0/, a(246)/0.020356797154333d0/
        data x(247)/0.759602341176647d0/, a(247)/0.021172939892191d0/
        data x(248)/0.738030643744400d0/, a(248)/0.021966644438744d0/
        data x(249)/0.715676812348967d0/, a(249)/0.022737069658329d0/
        data x(250)/0.692564536642171d0/, a(250)/0.023483399085926d0/
        data x(251)/0.668718310043916d0/, a(251)/0.024204841792364d0/
        data x(252)/0.644163403784967d0/, a(252)/0.024900633222483d0/
        data x(253)/0.618925840125468d0/, a(253)/0.025570036005349d0/
        data x(254)/0.593032364777572d0/, a(254)/0.026212340735672d0/
        data x(255)/0.566510418561397d0/, a(255)/0.026826866725591d0/
        data x(256)/0.539388108324357d0/, a(256)/0.027412962726029d0/
        data x(257)/0.511694177154667d0/, a(257)/0.027970007616848d0/
        data x(258)/0.483457973920596d0/, a(258)/0.028497411065085d0/
        data x(259)/0.454709422167743d0/, a(259)/0.028994614150555d0/
        data x(260)/0.425478988407300d0/, a(260)/0.029461089958167d0/ 
        data x(261)/0.395797649828908d0/, a(261)/0.029896344136328d0/ 
        data x(262)/0.365696861472313d0/, a(262)/0.030299915420827d0/ 
        data x(263)/0.335208522892625d0/, a(263)/0.030671376123669d0/ 
        data x(264)/0.304364944354496d0/, a(264)/0.031010332586313d0/ 
        data x(265)/0.273198812591049d0/, a(265)/0.031316425596861d0/ 
        data x(266)/0.241743156163840d0/, a(266)/0.031589330770727d0/ 
        data x(267)/0.210031310460567d0/, a(267)/0.031828758894411d0/ 
        data x(268)/0.178096882367618d0/, a(268)/0.032034456231992d0/ 
        data x(269)/0.145973714654896d0/, a(269)/0.032206204794030d0/
        data x(270)/0.113695850110665d0/, a(270)/0.032343822568575d0/
        data x(271)/0.081297495464425d0/, a(271)/0.032447163714064d0/
        data x(272)/0.048812985136049d0/, a(272)/0.032516118713868d0/
        data x(273)/0.016276744849602d0/, a(273)/0.032550614492363d0/
        wt=0.0d0 
        ct=0.0d0
        do 3 i=1,48
        wt(i)=a(225+i)   
    3   ct(i)=-x(225+i)
        do 4 i=96,49,-1
        ct(i)=-ct(97-i)
    4   wt(i)=wt(97-i)
        pj=0.0d0
        alj=0.0d0
        pi=3.141592653589793d0
c   determine the number of gauss points we use 
        
c   now we set 96
!        call gset(-1.0d0,1.0d0,16,ct,wt)
c   ct,wt  gauss zero and weight
            
c   set the jth legendre value 
        do 1 i=1,96
    1   pj(i)=legendre(ct(i),j)
        
c   calculate the integration
        do 2 i=1,96
    2   alj=alj+wt(i)*pj(i)*(ct(i)**l)*func(ct(i))*pi

        return
        end

c   this function calculate the cutoff
c   lambda is cutoff energy ,n adjust the sharp degree 

        function cutoff(lambda,n)
        implicit real*8 (a-h,o-z)
        real*8 lambda,t,expo
        real*8 cutoff
        integer n
        common /cpot/ v(6),xmev,ymev 
        t=dfloat(n)
        expo=(xmev/lambda)**(2.0d0*t)+(ymev/lambda)**(2.0d0*t)
        cutoff=dexp(-expo)  
        end

c   this subroutine calculate lsj decomposition
c   we use pot(6) to record them
c   0v(singlet), 1v(uncoupled triplet), v++, v--, v+-, v-+ (coupled)



        subroutine lsjdecomposition(pot,j)
        use kqxy
        use addterms
        implicit real*8 (a-h,o-z)
        real*8 pot(6),jd,jdp1,jd2p1
        real*8,external :: cutoff
        real*8,external :: alj
        integer j
        real*8 lambda
        common /cut/ lambda
        
        call initialize
        pot=0.0d0
        jd=dfloat(j)
        jdp1=jd+1.0d0
        jd2p1=2.0d0*jd+1.0d0

c   central force part
c   when j=0,there is only two terms do not equal to zero
        
        if (j .eq. 0) then
        pot(1)=2.0d0*alj(vcentral,0,j)
        pot(3)=2.0d0*alj(vcentral,0,j+1)
        else 
        pot(1)=2.0d0*alj(vcentral,0,j)
        pot(2)=2.0d0*alj(vcentral,0,j)
        pot(3)=2.0d0*alj(vcentral,0,j+1)
        pot(4)=2.0d0*alj(vcentral,0,j-1)
        pot(5)=0.0d0
        pot(6)=0.0d0
        end if

c   spin-spin force part

        if (j .eq. 0)then
        pot(1)=pot(1)-6.0d0*alj(vspinspin,0,j)
        pot(3)=pot(3)+2.0d0*alj(vspinspin,0,j+1)
        else 
        pot(1)=pot(1)-6.0d0*alj(vspinspin,0,j)
        pot(2)=pot(2)+2.0d0*alj(vspinspin,0,j)
        pot(3)=pot(3)+2.0d0*alj(vspinspin,0,j+1)
        pot(4)=pot(4)+2.0d0*alj(vspinspin,0,j-1)
        pot(5)=pot(5)+0.0d0
        pot(6)=pot(6)+0.0d0
        end if        

c   spin-obit force part
     
        if (j .eq. 0) then
        pot(1)=pot(1)+0.0d0
        pot(3)=pot(3)+2.0d0*x*y*(jd+2.0d0)/(2.0d0*jd+3.0d0)
     1  *(alj(vspinobit,0,j+2)-alj(vspinobit,0,j))
        else
        pot(1)=pot(1)+0.0d0
        pot(2)=pot(2)+2.0d0*x*y/jd2p1*(alj(vspinobit,0,j+1)
     1  -alj(vspinobit,0,j-1))
        pot(3)=pot(3)+2.0d0*x*y*(jd+2.0d0)/(2.0d0*jd+3.0d0)
     1  *(alj(vspinobit,0,j+2)-alj(vspinobit,0,j))
        pot(4)=pot(4)+2.0d0*x*y*(jd-1.0d0)/(2.0d0*jd-1.0d0)
     1  *(alj(vspinobit,0,j-2)-alj(vspinobit,0,j))
        pot(5)=pot(5)+0.0d0
        pot(6)=pot(6)+0.0d0
        end if

c   sigmaL force part
        
        if (j .eq. 0) then
        pot(1)=pot(1)+2.0d0*x2*y2*(alj(vsigmaL,2,j)-alj(vsigmaL,0,j))
        pot(3)=pot(3)+2.0d0*x2*y2*((2.0d0*jd+3.0d0)/jd2p1
     1  *alj(vsigmaL,0,j+1)-2.0d0/jd2p1*alj(vsigmaL,1,j)
     2  -alj(vsigmaL,2,j+1))
        else
        pot(1)=pot(1)+2.0d0*x2*y2*(alj(vsigmaL,2,j)-alj(vsigmaL,0,j))
        pot(2)=pot(2)+2.0d0*x2*y2*(-alj(vsigmaL,0,j)+(jd-1.0d0)/jd2p1
     1  *alj(vsigmaL,1,j+1)+(jd+2.0d0)/jd2p1*alj(vsigmaL,1,j-1))
        pot(3)=pot(3)+2.0d0*x2*y2*((2.0d0*jd+3.0d0)/jd2p1
     1  *alj(vsigmaL,0,j+1)-2.0d0/jd2p1*alj(vsigmaL,1,j)
     2  -alj(vsigmaL,2,j+1))
        pot(4)=pot(4)+2.0d0*x2*y2*((2*jd-1.0d0)/jd2p1
     1  *alj(vsigmaL,0,j-1)+2.0d0/jd2p1*alj(vsigmaL,1,j)
     2  -alj(vsigmaL,2,j-1))
        pot(5)=pot(5)-4.0d0*x2*y2*dsqrt(jd*jdp1)/(jd2p1)**2*
     1  (alj(vsigmaL,0,j+1)-alj(vsigmaL,0,j-1))
        pot(6)=pot(5)
        end if 
        
c   tensor force part

        if (j .eq. 0) then
        pot(1)=pot(1)+2.0d0*(-(x2+y2)*alj(vtensor,0,j)
     1  +2.0d0*x*y*alj(vtensor,1,j))
        pot(3)=pot(3)+2.0d0/jd2p1*(-(x2+y2)*alj(vtensor,0,j+1)
     1  +2.0d0*x*y*alj(vtensor,0,j))
        else
        pot(1)=pot(1)+2.0d0*(-(x2+y2)*alj(vtensor,0,j)
     1  +2.0d0*x*y*alj(vtensor,1,j))
        pot(2)=pot(2)+2.0d0*((x2+y2)*alj(vtensor,0,j)
     1  -2.0d0*x*y/jd2p1*(jd*alj(vtensor,0,j+1)
     2  +jdp1*alj(vtensor,0,j-1)))
        pot(3)=pot(3)+2.0d0/jd2p1*(-(x2+y2)*alj(vtensor,0,j+1)
     1  +2.0d0*x*y*alj(vtensor,0,j))
        pot(4)=pot(4)+2.0d0/jd2p1*((x2+y2)*alj(vtensor,0,j-1)
     1  -2.0d0*x*y*alj(vtensor,0,j))
        pot(5)=pot(5)-4.0d0*dsqrt(jd*jdp1)/jd2p1*(y2*alj(vtensor,0,j+1)
     1  +x2*alj(vtensor,0,j-1)-2.0d0*x*y*alj(vtensor,0,j) )   
        pot(6)=pot(6)-4.0d0*dsqrt(jd*jdp1)/jd2p1*(y2*alj(vtensor,0,j-1)
     1  +x2*alj(vtensor,0,j+1)-2.0d0*x*y*alj(vtensor,0,j) )
        end if 

c   sigmak force part

      if (j .eq. 0) then
        pot(1)=pot(1)+0.5d0*(-(x2+y2)*alj(vsigmak,0,j)
     1  -2.0d0*x*y*alj(vsigmak,1,j))
        pot(3)=pot(3)+0.5d0/jd2p1*(-(x2+y2)*alj(vsigmak,0,j+1)
     1  -2.0d0*x*y*alj(vsigmak,0,j))
        else
        pot(1)=pot(1)+0.5d0*(-(x2+y2)*alj(vsigmak,0,j)
     1  -2.0d0*x*y*alj(vsigmak,1,j))
        pot(2)=pot(2)+0.5d0*((x2+y2)*alj(vsigmak,0,j)
     1  +2.0d0*x*y/jd2p1*(jd*alj(vsigmak,0,j+1)
     2  +jdp1*alj(vsigmak,0,j-1)))
        pot(3)=pot(3)+0.5d0/jd2p1*(-(x2+y2)*alj(vsigmak,0,j+1)
     1  -2.0d0*x*y*alj(vsigmak,0,j))
        pot(4)=pot(4)+0.5d0/jd2p1*((x2+y2)*alj(vsigmak,0,j-1)
     1  +2.0d0*x*y*alj(vsigmak,0,j))
        pot(5)=pot(5)-1.0d0*dsqrt(jd*jdp1)/jd2p1*(y2*alj(vsigmak,0,j+1)
     1  +x2*alj(vsigmak,0,j-1)+2.0d0*x*y*alj(vsigmak,0,j) )   
        pot(6)=pot(6)-1.0d0*dsqrt(jd*jdp1)/jd2p1*(y2*alj(vsigmak,0,j-1)
     1  +x2*alj(vsigmak,0,j+1)+2.0d0*x*y*alj(vsigmak,0,j) )
        end if 

c   multiply the cutoff and the overall factor 1/(2*pi)^3

c   there is another factor sqrt(m/E)

        ex=dsqrt(1.0d0+x*x)
        ey=dsqrt(1.0d0+y*y)
        ree=dsqrt(ex*ey)
        expexp=cutoff(lambda,2) 
        do 1 i=1,6
    1   pot(i)=pot(i)/(2.0d0*pi)**3*dwnq/ree*expexp
        end
        