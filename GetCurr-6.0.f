C
C   SUM STATES TO GET TUNNEL CURRENT
C
      SUBROUTINE GETcurr(IBAND,BIAS,SEP,E0,BARR,TK1,TK2,EF,EPSI2,
     &NSTATES,EFTIP,NBARR2,NVDIM2,BARRPROF,NEIGENDIM,CURR)
C
      DIMENSION EPSI2(4,NEIGENDIM),BARRPROF(NVDIM2)
      DOUBLE PRECISION SUM,SUM2
      REAL KAPPA
c   value of C below is 2m/hbar^2 in units of 1/(eV nm^2)
      DATA C/26.254/RQUANT/12900./
      PI=4.*ATAN(1.)
      WKFTIP=sqrt(C*EFTIP)
C
      SUM=0.
      do 120 I=1,NSTATES
         ener=E0+EPSI2(1,I)
         occ=fd(ener-bias,ef,tk2)-fd(ener,ef,tk1)
C         trans=EPSI2(2,I)*exp(-2.*sep*EPSI2(3,I))
C         KAPPA=sqrt(C*AMAX1(0.,BARR-ENER)+EPSI2(4,I))
C         TRANS=EPSI2(2,I)*exp(-2.*SEP*KAPPA)
         SUM2=0.D0
         DELVAC=sep/float(NBARR2-1)
         do 100 J=1,NBARR2
            TERM=SQRT(C*amax1(0.,BARRPROF(J)-ENER)+EPSI2(4,I))
            IF (J.EQ.1.OR.J.EQ.NBARR2) TERM=TERM/2.
            SUM2=SUM2+TERM
100      CONTINUE
         SKAPPA=SUM2*DELVAC
c         write(6,*) 'sep*kappa =',SEP*KAPPA,skappa
c         read(5,*)
c         stop
         TRANS=EPSI2(2,I)*exp(-2.*SKAPPA)
         SUM=SUM+TRANS*OCC
120   continue
      CURR=16.*PI*WKFTIP*SUM/(C*RQUANT)
      return
      end
