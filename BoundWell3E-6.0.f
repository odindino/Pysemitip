c   *********************** BoundEner3E.f ********************
c
C   COMPUTE BASIS STATES FOR PERIODIC SYSTEM CONSISTING OF SEMICONDUCTOR
C   SLAB (A QUANTUM WELL) AND VACUUM REGION, USING EFFECTIVE MASS IN 
C   THE SEMICONDUCTOR.
c   SLAB THICKNESS IS EL1, AND SLAB+VACUUM THICKNESS IS EL2.
C   HANDLES CONDUCTION BAND CASE (BOTH BOUND AND EXTENDED STATES).
c
      subroutine BoundEner3E(effm,el1,el2,v0,ener,amp1,amp2,n,nstate)
c
c   ener is energy of bound states (zero is bottom of well, with depth v0)
c   amp1 is amplitude OUTSIDE of well
c   amp2 is amplitude INSIDE of well
c   n is the number of bound states, nstate is the number of states desired
c
      IMPLICIT REAL*8(A-H,O-Z)
      dimension ener(100),amp1(100),amp2(100),tmp(100)
      logical err,specfinal
      external feven,fodd,feven2,fodd2
      common/qw/c,effm1,alpha,v00,func
c
      v00=v0
      effm1=effm
      pi=4.*datan(1.d0)
c   'a' is half the well width in this routine
      a=el1/2
      el=el2/2.
      alpha=(el/a)-1.d0
      hbar22m=197.3**2/(2.*0.511e6)
      c=dsqrt(abs(v0)*a**2/hbar22m)
      n=0
c
c   skip bound states if necessary
c
      if (v0.lt.0.) go to 170
      n=1+int(2.*c*sqrt(effm)/pi)
      specfinal=.true.
c
c   *********** bound states *********
c
c   check if final state is special one (usually is, except in
c   unusual circumstance for odd parity, i.e. n even, in which case
c   it doesn't exist)
c
      if (mod(n,2).eq.0) then
         delxsi0=datan(effm/(alpha*(n-1)*pi/2.))
         do 50 i=1,100
            delxsi1=datan(effm/((alpha*(n-1)*pi/2.)+delxsi0))
            if (abs(delxsi1-delxsi0).lt.1.e-6) go to 60
            delxsi0=delxsi1
50       continue
         write(6,*) '*** error - BoundWell, delxsi, too many iterations'
c         read(5,*)
60       delxsi=delxsi1    
         if ((c*sqrt(effm)-(n-1)*pi/2.).le.delxsi) then
            specfinal=.false.
            n=n-1
         end if
      end if
c      write(6,*) 'number of bound states =',n
c
c   non-final state
c
      do 150 i=1,n
         xsi=i*pi/2.
         if (i.ne.n.or.(.not.specfinal)) then
            do 100 iter=1,1000
               xsiold=xsi
               eta=dsqrt(c**2-xsiold**2/effm)
               if (mod(i,2).eq.1) then
               xsi=datan(effm*eta*dtanh(eta*alpha)/xsiold)+(i-1)*pi/2
               else
               xsi=datan(effm*eta/(xsiold*dtanh(eta*alpha)))+(i-1)*pi/2
               end if
               if (dabs(xsi-xsiold).lt.1.d-6) go to 110
100         continue
            write(6,*) '*** error - reached max iterations'
c            read(5,*)
110         continue
         else
c
c   final state
c
            xmin=(i-1)*pi/2.
            if (mod(n,2).eq.0) xmin=xmin+delxsi
            xmax=c*sqrt(effm)
            ep=1.d-6
            if (mod(i,2).eq.1) then
               call dgsect(xmin,xmax,ep,feven)
            else
               call dgsect(xmin,xmax,ep,fodd)
            end if
            xsi=(xmin+xmax)/2.
         end if
         ener(i)=hbar22m*(xsi/a)**2/effm
c
c   check
c
         xsi=a*sqrt(effm*ener(i)/hbar22m)
         eta=a*sqrt(abs(v0-ener(i))/hbar22m)
c         write(6,*) i,xsi,eta
         err=.false.
         if (mod(i,2).eq.1) then
            if ((abs((xsi*tan(xsi)/effm)/(eta*dtanh(eta*alpha)))-1)
     &            .gt.1.d-4) err=.true.
         else
            if ((abs(xsi*dtanh(eta*alpha)/(tan(xsi)*effm*eta))-1)
     &            .gt.1.d-4) err=.true.
         end if
         if (err) THEN
C            write(6,*) '*** error in solution, i=',i
C            write(16,*) '*** error in solution, i=',i
         END IF
c
c   get amplitudes
c
         wk=dsqrt(ener(i)*effm/hbar22m)
         wkappa=dsqrt((v0-ener(i))/hbar22m)
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
c   get extended states, and order them appropriately
c
c   ********* extended even states **********
c
170   if (v0.ge.0.) then
         xsi=c*sqrt(effm)
      else
         xsi=0.
      end if
      i=n
180   if (i.ge.nstate) go to 300
         i=i+1
185      tmp1=0.
         iter=0
190      if (tmp1*tmp1sav.lt.0.) go to 200
            iter=iter+1
            tmp1sav=tmp1
            xsiold=xsi
            xsi=xsi+0.01
            if (v0.ge.0.) then
               eta=dsqrt(abs((xsi**2/effm)-c**2))
            else
               eta=dsqrt((xsi**2/effm)+c**2)
            end if
            tmp1=xsi*tan(xsi)/effm+eta*tan(eta*alpha)
c            write(6,*) i,iter,xsi,eta
c            read(5,*)
            go to 190
200      continue
c         write(6,*) 'found even soln !',i,iter,xsi,eta
         ep=1.d-6
         xmin=xsiold
         xmax=xsi
         call dgsect(xmin,xmax,ep,feven2)
         xsi=(xmin+xmax)/2.
c         write(6,*) 'found finer even soln !',i,iter,xsi,eta
c         read(5,*)
         if (func.gt.1.d-6) go to 185
         ii=n+2*(i-n-1)+1+mod(n,2)
c         write(6,*) 'i,ii =',i,ii
c         read(5,*)
         ener(ii)=hbar22m*(xsi/a)**2/effm
c
c   check
c
         xsi=a*sqrt(effm*ener(ii)/hbar22m)
         eta=a*sqrt(abs(ener(ii)-v0)/hbar22m)
c         write(6,*) ii,xsi,eta,sngl((ii-1)*pi/2.)
         err=.false.
         if ((abs((xsi*tan(xsi)/effm)/(eta*dtan(eta*alpha)))-1)
     &         .gt.1.d-4) err=.true.
         if (err) THEN
C            write(6,*) '*** error in even solution, i=',i
C            write(16,*) '*** error in even solution, i=',i
         END IF
c
c   get amplitudes
c
         wk=dsqrt(ener(ii)*effm/hbar22m)
         wkappa=dsqrt((ener(ii)-v0)/hbar22m)
c         write(6,*) '1 ii,wk,wkappa =',ii,wk,wkappa
c         wk=xsi/a
c         wkappa=eta/a
c         write(6,*) '2 ii,wk,wkappa =',ii,wk,wkappa
         atmp=(a+sin(2.*wk*a)/(2.*wk))+
     &         (cos(wk*a)/cos(wkappa*(el-a)))**2*
     &         (el-a+sin(2.*wkappa*(el-a))/(2.*wkappa))
         amp2(ii)=1./dsqrt(atmp)
         amp1(ii)=amp2(ii)*cos(wk*a)/cos(wkappa*(el-a))
c         tmp(ii)=0.
c         do 240 ix=0,nint(el*100)
c           x=ix/100.
c           if (x.le.a) then
c              psi=amp2(ii)*cos(wk*x)
c           else
c           psi=amp1(ii)*cos(wkappa*(el-x))
c           end if
c           write(10+ii,*) x,psi
c           tmp(ii)=tmp(ii)+psi**2*0.01
c240      continue
c         write(6,*) 'ii,sum =',ii,2.*tmp(ii)
c         if (mod(i,10).eq.0) read(5,*)
250   continue
      go to 180
300   continue
c
c   ********* extended odd states **********
c
370   if (v0.ge.0.) then
         xsi=c*sqrt(effm)
      else
         xsi=0.
      end if
      i=n
380   if (i.ge.nstate) go to 500
         i=i+1
385      tmp1=0.
         iter=0
390      if (tmp1*tmp1sav.lt.0.) go to 400
            iter=iter+1
            tmp1sav=tmp1
            xsiold=xsi
            xsi=xsi+0.01
            if (v0.ge.0.) then
               eta=dsqrt(abs((xsi**2/effm)-c**2))
            else
               eta=dsqrt((xsi**2/effm)+c**2)
            end if
            tmp1=tan(xsi)*effm/xsi+tan(eta*alpha)/eta
c            write(6,*) i,iter,xsi,eta
c            read(5,*)
            go to 390
400      continue
c         write(6,*) 'found odd soln !',i,iter,xsi,eta
         ep=1.d-6
         xmin=xsiold
         xmax=xsi
         call dgsect(xmin,xmax,ep,fodd2)
         xsi=(xmin+xmax)/2.
c         write(6,*) 'found finer odd soln !',i,iter,xsi,eta
c         read(5,*)
         if (func.gt.1.d-6) go to 385
         ii=n+2*(i-n-1)+1+mod(n+1,2)
c         write(6,*) 'i,ii =',i,ii
c         read(5,*)
         ener(ii)=hbar22m*(xsi/a)**2/effm
c
c   check
c
         xsi=a*sqrt(effm*ener(ii)/hbar22m)
         eta=a*sqrt(abs(ener(ii)-v0)/hbar22m)
c         write(6,*) ii,xsi,eta,sngl((ii-1)*pi/2.)
         err=.false.
         if ((abs((tan(xsi)*effm*eta)/(xsi*dtan(eta*alpha)))-1)
     &         .gt.1.d-4) err=.true.
         if (err) THEN
C            write(6,*) '*** error in odd solution, i=',i
C            write(16,*) '*** error in odd solution, i=',i
         END IF
c
c   get amplitudes
c
         wk=dsqrt(ener(ii)*effm/hbar22m)
         wkappa=dsqrt((ener(ii)-v0)/hbar22m)
c         write(6,*) '1 ii,wk,wkappa =',ii,wk,wkappa
c         wk=xsi/a
c         wkappa=eta/a
c         write(6,*) '2 ii,wk,wkappa =',ii,wk,wkappa
         atmp=(a-sin(2.*wk*a)/(2.*wk))+
     &         (sin(wk*a)/sin(wkappa*(el-a)))**2*
     &         (el-a-sin(2.*wkappa*(el-a))/(2.*wkappa))
         amp2(ii)=1./dsqrt(atmp)
         amp1(ii)=amp2(ii)*sin(wk*a)/sin(wkappa*(el-a))
c         tmp(ii)=0.
c         do 440 ix=0,nint(el*100)
c           x=ix/100.
c           if (x.le.a) then
c              psi=amp2(ii)*sin(wk*x)
c           else
c           psi=amp1(ii)*sin(wkappa*(el-x))
c           end if
c           write(10+ii,*) x,psi
c           tmp(ii)=tmp(ii)+psi**2*0.01
c440      continue
c         write(6,*) 'ii,sum =',ii,2.*tmp(ii)
c         if (mod(i,10).eq.0) read(5,*)
450   continue
      go to 380
500   continue
      return
      end
c
c   function for even parity bound states
c
      function feven(xsi)
      IMPLICIT REAL*8(A-H,O-Z)
      common/qw/c,effm,alpha,v0,func
      eta=dsqrt(dabs(c**2-xsi**2/effm))
      feven=(xsi*dtan(xsi)/(effm*dtanh(eta*alpha))-eta)**2
      return
      end
c
c   function for odd parity bound state
c
      function fodd(xsi)
      IMPLICIT REAL*8(A-H,O-Z)
      common/qw/c,effm,alpha,v0,func
      eta=dsqrt(dabs(c**2-xsi**2/effm))
      fodd=(xsi*dtanh(eta*alpha)/(dtan(xsi)*effm)+eta)**2
      return
      end
c
c   function for even parity bound states
c
      function feven2(xsi)
      IMPLICIT REAL*8(A-H,O-Z)
      common/qw/c,effm,alpha,v0,func
      if (v0.ge.0.) then
         eta=dsqrt(abs((xsi**2/effm)-c**2))
      else
         eta=dsqrt((xsi**2/effm)+c**2)
      end if
      feven2=(xsi*tan(xsi)/effm+eta*tan(eta*alpha))**2
c      write(6,*) 'xsi,feven2=',xsi,feven2
      func=feven2
      return
      end
c
c   function for odd parity bound state
c
      function fodd2(xsi)
      IMPLICIT REAL*8(A-H,O-Z)
      common/qw/c,effm,alpha,v0,func
      if (v0.ge.0.) then
         eta=dsqrt(abs((xsi**2/effm)-c**2))
      else
         eta=dsqrt((xsi**2/effm)+c**2)
      end if
      fodd2=(tan(xsi)*effm/xsi+tan(eta*alpha)/eta)**2
c      write(6,*) 'xsi,fodd2=',xsi,fodd2
      func=fodd2
      return
      end
