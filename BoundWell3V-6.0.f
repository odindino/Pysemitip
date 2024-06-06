c   *********************** BoundEner3V ***********************
c
C   COMPUTE BASIS STATES FOR PERIODIC SYSTEM CONSISTING OF SEMICONDUCTOR
C   SLAB (A QUANTUM WELL) AND VACUUM REGION, USING EFFECTIVE MASS IN 
C   THE SEMICONDUCTOR.
c   SLAB THICKNESS IS EL1, AND SLAB+VACUUM THICKNESS IS EL2.
c   HANDLES VALENCE BAND CASE (BOUND STATES ONLY).
c
      subroutine BoundEner3V(effm,el1,el2,v0,ener,amp1,amp2,n)
c
c   ener is energy of bound states (zero is bottom of well, with depth v0)
c   amp1 is amplitude OUTSIDE of well
c   amp2 is amplitude INSIDE of well
c
      IMPLICIT REAL*8(A-H,O-Z)
      dimension ener(100),amp1(100),amp2(100),tmp(100)
      logical err
c
      effm1=effm
      pi=4.*datan(1.d0)
c   'a' is half the well width in this routine
      a=el1/2
      el=el2/2.
      alpha=(el/a)-1.d0
      hbar22m=197.3**2/(2.*0.511e6)
      c=dsqrt(v0*a**2/hbar22m)
c
c   regular state
c
      do 150 i=1,n
         xsi=i*pi/2.
            do 100 iter=1,1000
               xsiold=xsi
               eta=dsqrt(c**2+xsiold**2/effm)
               if (mod(i,2).eq.1) then
               xsi=datan(effm*eta*dtanh(eta*alpha)/xsiold)+(i-1)*pi/2
               else
               xsi=datan(effm*eta/(xsiold*dtanh(eta*alpha)))+(i-1)*pi/2
               end if
               if (dabs(xsi-xsiold).lt.1.d-6) go to 110
100         continue
            write(6,*) '*** error - reached max iterations'
110         continue
         ener(i)=hbar22m*(xsi/a)**2/effm
c
c   check
c
         xsi=a*sqrt(effm*ener(i)/hbar22m)
         eta=a*sqrt(abs(v0+ener(i))/hbar22m)
         err=.false.
         if (mod(i,2).eq.1) then
            if ((abs((xsi*tan(xsi)*dtanh(eta*alpha)/effm)/eta)-1)
     &         .gt.1.d-4) err=.true.
         else
            if ((abs(xsi*dtanh(eta*alpha)/(tan(xsi)*effm*eta))-1)
     &         .gt.1.d-4) err=.true.
         end if
C         if (err) write(6,*) '*** error in solution, i=',i
c
c   get amplitudes
c
         wk=dsqrt(ener(i)*effm/hbar22m)
         wkappa=dsqrt((v0+ener(i))/hbar22m)
c         write(6,*) '1 i,wk,wkappa =',i,wk,wkappa
c         wk=xsi/a
c         wkappa=eta/a
c         write(6,*) '2 i,wk,wkappa =',i,wk,wkappa
         if (mod(i,2).eq.1) then
            atmp=(a+sin(2.*wk*a)/(2.*wk))+
     &            (cos(wk*a)/cosh(wkappa*(el-a)))**2*
     &            (el-a+sinh(2.*wkappa*(el-a))/(2.*wkappa))
            amp2(i)=1./dsqrt(atmp)
            amp1(i)=amp2(i)*cos(wk*a)/cosh(wkappa*(el-a))
         else
            atmp=(a-sin(2.*wk*a)/(2.*wk))+
     &            (sin(wk*a)/sinh(wkappa*(el-a)))**2*
     &            (-el+a+sinh(2.*wkappa*(el-a))/(2.*wkappa))
            amp2(i)=1./dsqrt(atmp)
            amp1(i)=amp2(i)*sin(wk*a)/sinh(wkappa*(el-a))
         end if
c         tmp(i)=0.
c         do 140 ix=0,nint(el*100)
c           x=ix/100.
c           if (x.le.a) then
c              if (mod(i,2).eq.1) then
c                 psi=amp2(i)*cos(wk*x)
c              else
c                 psi=amp2(i)*sin(wk*x)
c              end if
c           else
c              if (mod(i,2).eq.1) then
c              psi=amp1(i)*cosh(wkappa*(el-x))
c              else
c              psi=amp1(i)*sinh(wkappa*(el-x))
c              end if
c           end if
c           write(10+i,*) x,psi
c           tmp(i)=tmp(i)+psi**2*0.01
c140      continue
c         write(6,*) 'i,sum =',i,2.*tmp(i)
c         if (mod(i,10).eq.0) read(5,*)
150   continue
c
      return
      end
