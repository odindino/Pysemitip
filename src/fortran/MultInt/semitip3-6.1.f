C   ******************** SEMITIP3 ************************
C
C   E-FIELD COMPUTATIONS FOR TIP IN PROXIMITY TO SEMICONDUCTOR, IN 3D
C
C   VERSION 1.0 - WRITTEN BY R. M. FEENSTRA, DEC 2004
C                 BASED ON SEMITIP VERSION 2.1
C           1.1 - MAR/05, MODIFY GSECT2,3 STEP SIZE IF BIAS=0
C           1.2 - DEC/05, ADD IINV PARAMETER
C           1.3 - APR/06, CHANGE INTERPOLATION FOR Pot0 (WITH DP)
C           1.4 - APR/06, ADD BOUNDARY CONDITION AT HETEROINTERFACE
C           2.0 - APR/07, CHANGE COORDINATE SYSTEM TO EXACTLY MATCH TIP
C                         (INTRODUCE C,Z0; MODIFY DELV,A,ETAT)
C           3.0 - NOV/10, INTRODUCED CALL TO RHOSURF FOR SURFACE CHARGE
C           6.0 - FEB/11, TOTAL REWRITE, COMBINE WITH SEMITIP2_6.0; 
C                         MODIFIED LOOPS ON DOUBLING AND REMOVE HETEROINTERFACE
C           6.1 - FEB/11. REVISED CALLS TO SEMMIN AND SURFMIN, CORRECTED 
C                         INTERPOLATION FOR DOUBLING POTENTIAL ARRAYS
C
C   CALLS ROUTINES RHOSURF AND RHOBULK SUPPLIED BY USER, FOR EVALUATING CHARGE DENSITIES.
C
C   INPUT AND OUTPUT PARAMETERS:
C      SEP=SEPARATION BETWEEN SEMICONDUCTOR AND END OF TIP (NM)
C      RAD=RADIUS OF HYPERBOLOID WHICH FORMS TIP (NM)
C      SLOPE=SLOPE OF TIP SHANK
C      ETAT=ETA PARAMETER FOR THE HYPERBOLIC COORDINATE SYSTEM
C      A=LOCATION OF FOCUS FOR THE HYPERBOLIC COORDINATE SYSTEM
C      Z0=LOCATION OF ORIGIN FOR THE HYPERBOLIC COORDINATE SYSTEM
C      C=ORIGIN PARAMETER FOR FOR THE HYPERBOLIC COORDINATE SYSTEM
C      VAC=ARRAY OF POTENTIAL VALUES IN VACUUM
C      SEM=ARRAY OF POTENTIAL VALUES IN SEMICONDUCTOR
C      TIP=ARRAY WITH VALUES OF .TRUE. IF INSIDE TIP, .FALSE. OTHERWISE
C      VSINT=ARRAY OF POTENTIAL VALUES AT SEMICONDUCTOR SURFACE
C      R=ARRAY OF R VALUES
C      S=ARRAY OF Z VALUES IN SEMICONDUCTOR
C      DELV=ARRAY OF Z-SPACING VALUES IN VACUUM
C      DELR=R-SPACING (FOR SMALL R) IN VACUUM AND SEMICONDUCTOR
C      DELS=Z-SPACING (FOR SMALL Z) IN SEMICONDUCTOR
C      DELP=SPACING OF ANGLES (RADIANS)
C      NRDIM=R-DIMENSION FOR VAC, SEM, AND VSINT ARRAYS
C      NVDIM=Z-DIMENSION FOR VAC ARRAY
C      NSDIM=Z-DIMENSION FOR SEM ARRAY
C      NPDIM=ANGULAR DIMENSION FOR VAC, SEM, AND VSINT ARRAYS
C      NR=CURRENT VALUE OF R-DIMENSION FOR VAC, SEM, AND VSINT ARRAYS
C      NV=CURRENT VALUE OF Z-DIMENSION FOR VAC ARRAY
C      NS=CURRENT VALUE OF Z-DIMENSION FOR SEM ARRAY
C      NP=CURRENT VALUE OF ANGULAR DIMENSION FOR VAC, SEM, AND VSINT ARRAYS
C      BIAS=BIAS OF TIP RELATIVE TO A POINT FAR INSIDE THE SEMICONDUCTOR
C      IWRIT >= 1 IF WANT WRITTEN OUTPUT, 0 OTHERWISE
C      ITMAX=ARRAY OF MAXIMUM ITERATION SPECIFICATIONS
C      EP=ARRAY OF STOPPING CRITERIONS
C      IPMAX=MAXIMUM NUMBER OF SCALING LOOPS FOR GRID DOUBLING
C      Pot0=MAX SURFACE BAND BENDING
C      IERR=ERROR INDICATOR
C      IINIT=1 TO INITIALIZE THE POTENTIAL WITH FIRST GUESS, OR 0 NOT TO
C      EPSIL=DIELECTRIC CONSTANT OF SEMICONDUCTOR
C
C
C   SEMMIN USED TO FIND OPTIMAL UPDATE EQUATION FOR SOLN IN SEMICONDUCTOR
C
      FUNCTION SEMMIN(Pot)
      COMMON/SEMSURFMIN/EPSIL,X,Y,S,STEMP,DENOM,I,J,K,NR,NS,NP
      DATA EEP/1.80943E-20/
      RHO=RHOBULK(Pot,X,Y,S,I,J,K,NR,NS,NP)
      TEMP=STEMP-RHO*EEP/EPSIL
      SEMMIN=ABS(Pot-TEMP/DENOM)
      RETURN
      END
C
C   SURFMIN USED TO FIND OPTIMAL UPDATE EQUATION FOR SOLN ON SURFACE
C
      FUNCTION SURFMIN(Pot)
      COMMON/SEMSURFMIN/EPSIL,X,Y,S,STEMP,DENOM,I,J,K,NR,NS,NP
      DATA EEP/1.80943E-20/
      RHO=RHOSURF(Pot,X,Y,I,K,NR,NP)
      TEMP=STEMP-RHO*EEP*1.E7
      SURFMIN=ABS(Pot-TEMP/DENOM)
      RETURN
      END
C
C   MAIN PROGRAM
C
      SUBROUTINE SEMITIP3(SEP,RAD,SLOPE,ETAT,A,Z0,C,VAC,TIP,SEM,
     &VSINT,R,S,DELV,DELR0,DELS0,DELP,NRDIM,NVDIM,NSDIM,NPDIM,NR,NV,
     &NS,NP,BIAS,IWRIT,ITMAX,EP,IPMAX,Pot0,IERR,IINIT,MIRROR,EPSIL)
C
      DIMENSION VAC(2,NRDIM,NVDIM,NPDIM),SEM(2,NRDIM,NSDIM,NPDIM),
     &VSINT(2,NRDIM,NPDIM),R(NRDIM),DELR(NRDIM),DELV(NRDIM),S(NSDIM),
     &DELS(NSDIM),ITMAX(10),EP(10),DELXSI(NRDIM)
      LOGICAL TIP(NRDIM,NVDIM,NPDIM)
      PI=4.*ATAN(1.)
C   BOUNDARY CONDITION (0=DIRICHLET, 1=VON NEUMANN)
      IBC=0
      IERR=0
C
C   CONSTRUCT THE TIP AND VACUUM GRID AND INITIALIZE POTENTIAL
C
      ETAT=1./SQRT(1.+1./SLOPE**2)
      A=RAD*SLOPE**2/ETAT
      SPRIME=A*ETAT
      Z0=SEP-SPRIME
      C=Z0/SPRIME
      DELETA=ETAT/FLOAT(NV)
      CETAT=ALOG((1.+ETAT)/(1.-ETAT))
      IF (IWRIT.EQ.0) GO TO 100
         WRITE(6,*) 'ETAT, A, Z0, C =',ETAT,A,Z0,C
         WRITE(16,*) 'ETAT, A, Z0, C =',ETAT,A,Z0,C
         WRITE(6,*) 'NR,NS,NV,NP =',NR,NS,NV,NP
         WRITE(16,*) 'NR,NS,NV,NP =',NR,NS,NV,NP
         WRITE(6,90) DELR0,DELS0,(1.+C)*A*DELETA,DELP
         WRITE(16,90) DELR0,DELS0,(1.+C)*A*DELETA,DELP
90       FORMAT(' DELR,DELS,DELV,DELP =',4G13.5)
100   CONTINUE
      IF (NR.GT.NRDIM.OR.NV.GT.NVDIM.OR.NS.GT.NSDIM.OR.NP.GT.NPDIM) 
     &GO TO 600
      DO 150 I=1,NR
         R(I)=(2*NR*DELR0/PI)*TAN(PI*(I-0.5)/(2.*NR))
         X2M1=(R(I)/A)**2
         IF (I.NE.1) XSISAV=XSI
         XSI=DSQRT(1.D0+X2M1)
         IF (I.EQ.1) THEN
            DELR(I)=R(I)
            DELXSI(I)=DSQRT(1.D0+X2M1)-1.D0
         ELSE
            DELR(I)=R(I)-R(I-1)
            DELXSI(I)=XSI-XSISAV
         END IF
         DELV(I)=(SQRT(A**2+R(I)**2)+C*A)*ETAT/FLOAT(NV)
         DO 120 J=1,NV-1
            ETA=J*ETAT/FLOAT(NV)
            Z=A*ETA*(XSI+C)
            ZP=Z*(J+0.5)/FLOAT(J)
            RP=A*SQRT(X2M1*(1.-ETA**2))
            IF (ZP.GT.(A*ETAT*(SQRT(1.+RP**2/((1.-ETAT**2)*A**2))+C)
     &           -P(RP))) GO TO 110
            DO 105 K=1,NP
               IF (IINIT.EQ.1) THEN
                  VAC(1,I,J,K)=0.
                  VAC(2,I,J,K)=0.
               END IF
               TIP(I,J,K)=.FALSE.
105         CONTINUE
            GO TO 120
110         DO 115 K=1,NP
               VAC(1,I,J,K)=BIAS
               VAC(2,I,J,K)=BIAS
               TIP(I,J,K)=.TRUE.
115         CONTINUE
120      CONTINUE
         DO 130 K=1,NP
            VAC(1,I,NV,K)=BIAS
            VAC(2,I,NV,K)=BIAS
            TIP(I,NV,K)=.TRUE.
130      CONTINUE
150   CONTINUE
      IF (IINIT.EQ.1) THEN
         DO 160 I=1,NR
            DO 155 K=1,NP
               VSINT(1,I,K)=0.
               VSINT(2,I,K)=0.
155         CONTINUE
160      CONTINUE
      END IF
      DO 180 J=1,NS
         S(J)=(2*NS*DELS0/PI)*TAN(PI*(J-0.5)/(2.*NS))
         IF (J.EQ.1) THEN
            DELS(J)=S(J)
         ELSE
            DELS(J)=S(J)-S(J-1)
         END IF
         IF (IINIT.EQ.1) THEN
            DO 170 I=1,NR
               DO 165 K=1,NP
                  SEM(1,I,J,K)=0.
                  SEM(2,I,J,K)=0.
165            CONTINUE
170         CONTINUE
         END IF
180   CONTINUE
      IF (.NOT.TIP(1,1,1)) GO TO 190
         IERR=1
         WRITE(6,*) '*** ERROR - VACUUM GRID SPACING TOO LARGE'
         WRITE(16,*) '*** ERROR - VACUUM GRID SPACING TOO LARGE'
190   CONTINUE
      IF (IWRIT.NE.0) THEN
      WRITE(6,*) 'LARGEST RADIUS, DEPTH =',R(NR),S(NS)
      WRITE(16,*) 'LARGEST RADIUS, DEPTH =',R(NR),S(NS)
      END IF
      IF (IERR.EQ.1) GO TO 600
C
C   INITIAL GUESS
C
      DO 260 K=1,NP
         DO 250 I=1,NR
            DO 210 J=NV,1,-1
               IF (.NOT.TIP(I,J,K)) GO TO 220
210         CONTINUE
220         JSAV=J
            IF (IINIT.EQ.1) THEN
               DO 230 J=JSAV,1,-1
                  ETA=J*DELETA*FLOAT(NV)/FLOAT(JSAV+1)
                  VAC(1,I,J,K)=BIAS*ALOG((1.+ETA)/(1.-ETA))/CETAT
                  VAC(2,I,J,K)=VAC(1,I,J,K)
230            CONTINUE
            END IF
250      CONTINUE
260   CONTINUE
      Pot0=PCENT(0,VAC,SEM,VSINT,NRDIM,NVDIM,NSDIM,NPDIM,NP)
C
C   ITERATE POISSON'S EQN
C
      DO 590 IP=1,IPMAX
C
         IF (IWRIT.NE.0) WRITE(6,*) 'SOLUTION #',IP
         IF (IWRIT.NE.0) WRITE(16,*) 'SOLUTION #',IP
         ITM=ITMAX(IP)
         EPI=EP(IP)
         CALL ITER3(VAC,TIP,SEM,VSINT,R,DELR,DELV,DELXSI,S,DELS,BIAS,
     &       DELR0,DELS0,DELP,DELETA,A,NRDIM,NVDIM,NSDIM,NPDIM,NR,NV,
     &       NS,NP,EPI,ITM,Pot0,IWRIT,ETAT,C,MIRROR,IERR,EPSIL,IBC)
         IF (IP.EQ.IPMAX) GO TO 600
C
         IF (NR*2.GT.NRDIM.OR.NV*2.GT.NVDIM.OR.NS*2.GT.NSDIM.OR.
     &       NP*2.GT.NPDIM) GO TO 500
C
C   DOUBLE ALL ARRAYS
C
         NR=NR*2
         NS=NS*2
         NV=NV*2
         NP=NP*2
         DELR0=DELR0/2.
         DELS0=DELS0/2.
         DELETA=DELETA/2.
         DELP=DELP/2.
         IF (IWRIT.EQ.0) GO TO 300
            WRITE(6,*) 'NR,NS,NV,NP =',NR,NS,NV,NP
            WRITE(16,*) 'NR,NS,NV,NP =',NR,NS,NV,NP
            WRITE(6,90) DELR0,DELS0,(1.+C)*A*DELETA,DELP
            WRITE(16,90) DELR0,DELS0,(1.+C)*A*DELETA,DELP
300      CONTINUE
         DO 450 I=1,NR
            R(I)=(2*NR*DELR0/PI)*TAN(PI*(I-0.5)/(2.*NR))
            X2M1=(R(I)/A)**2
            IF (I.NE.1) XSISAV=XSI
            XSI=DSQRT(1.D0+X2M1)
            IF (I.EQ.1) THEN
               DELR(I)=R(I)
               DELXSI(I)=DSQRT(1.D0+X2M1)-1.D0
            ELSE
               DELR(I)=R(I)-R(I-1)
               DELXSI(I)=XSI-XSISAV
            END IF
            DELV(I)=(SQRT(A**2+R(I)**2)+C*A)*ETAT/FLOAT(NV)
            DO 430 J=1,NV-1
               ETA=J*ETAT/FLOAT(NV)
               Z=A*ETA*(XSI+C)
               ZP=Z*(J+0.5)/FLOAT(J)
               RP=A*SQRT(X2M1*(1.-ETA**2))
               IF (ZP.LE.(A*ETAT*(SQRT(1.+RP**2/((1.-ETAT**2)*A**2))
     &            +C)-P(RP))) THEN
                  DO 405 K=1,NP
                     IF (J.EQ.1) THEN
                        TEMP=(VAC(2,(I-1)/2+1,1,(K-1)/2+1)+
     &                     VSINT(2,(I-1)/2+1,(K-1)/2+1))/2.
                     ELSE
                     IF (MOD(J,2).EQ.0) THEN
                        TEMP=VAC(2,(I-1)/2+1,J/2,(K-1)/2+1)
                     ELSE
                        TEMP=(VAC(2,(I-1)/2+1,J/2,(K-1)/2+1)+
     &                     VAC(2,(I-1)/2+1,(J+1)/2,(K-1)/2+1))/2.
                     END IF
                     END IF
                     VAC(1,I,J,K)=TEMP
                     TIP(I,J,K)=.FALSE.
405               CONTINUE
               ELSE
                  DO 420 K=1,NP   
                     VAC(1,I,J,K)=BIAS
                     TIP(I,J,K)=.TRUE.
420               CONTINUE
               END IF
430         CONTINUE
            DO 440 K=1,NP
               VAC(1,I,NV,K)=BIAS
               TIP(I,NV,K)=.TRUE.
440         CONTINUE               
450      CONTINUE
         DO 460 I=1,NR
            DO 455 K=1,NP
               VSINT(1,I,K)=VSINT(2,(I-1)/2+1,(K-1)/2+1)
455         CONTINUE               
460      CONTINUE
         DO 463 J=1,NS
            S(J)=(2*NS*DELS0/PI)*TAN(PI*(J-0.5)/(2.*NS))
            IF (J.EQ.1) THEN
               DELS(J)=S(J)
            ELSE
               DELS(J)=S(J)-S(J-1)
            END IF
463      CONTINUE
         DO 478 K=1,NP
            DO 476 I=1,NR
               SEM(1,I,1,K)=(SEM(2,(I-1)/2+1,1,(K-1)/2+1)+
     &                  VSINT(2,(I-1)/2+1,(K-1)/2+1))/2.
               DO 470 J=2,NS-2,2
                   SEM(1,I,J,K)=(3.*SEM(2,(I-1)/2+1,J/2,(K-1)/2+1)+
     &                     SEM(2,(I-1)/2+1,J/2+1,(K-1)/2+1))/4.                   
470            CONTINUE
               DO 475 J=3,NS-1,2
                   SEM(1,I,J,K)=(SEM(2,(I-1)/2+1,J/2,(K-1)/2+1)+
     &                     3.*SEM(2,(I-1)/2+1,J/2+1,(K-1)/2+1))/4.
475            CONTINUE
               SEM(1,I,NS,K)=SEM(2,(I-1)/2+1,NS/2,(K-1)/2+1)
476         CONTINUE
478      CONTINUE
         DO 499 K=1,NP
            DO 492 J=1,NV
               DO 490 I=1,NR
                  VAC(2,I,J,K)=VAC(1,I,J,K)
490            CONTINUE
492         CONTINUE
            DO 494 I=1,NR
               VSINT(2,I,K)=VSINT(1,I,K)
494         CONTINUE
            DO 498 J=1,NS
               DO 496 I=1,NR
                  SEM(2,I,J,K)=SEM(1,I,J,K)
496            CONTINUE
498         CONTINUE
499      CONTINUE
         WRITE(6,*) 'LARGEST RADIUS, DEPTH =',R(NR),S(NS)
         WRITE(16,*) 'LARGEST RADIUS, DEPTH =',R(NR),S(NS)
         GO TO 590
C
500      IF (NS*2.GT.NSDIM) GO TO 600
C
C   DOUBLE ONLY SEM ARRAY ELEMENTS
C
         NS=NS*2
         DELS0=DELS0/2.
         IF (IWRIT.EQ.0) GO TO 505
            WRITE(6,*) 'NS =',NS
            WRITE(6,*) 'DELS =',DELS0
            WRITE(16,*) 'NS =',NS
            WRITE(16,*) 'DELS =',DELS0
505      CONTINUE
         DO 508 J=1,NS
            S(J)=(2*NS*DELS0/PI)*TAN(PI*(J-0.5)/(2.*NS))
            IF (J.EQ.1) THEN
               DELS(J)=S(J)
            ELSE
               DELS(J)=S(J)-S(J-1)
            END IF
508      CONTINUE
         DO 540 K=1,NP
            DO 532 I=1,NR
               SEM(1,I,1,K)=(SEM(2,I,1,K)+VSINT(2,I,K))/2.
               DO 520 J=2,NS-2,2
                  SEM(1,I,J,K)=(3.*SEM(2,I,J/2,K)+SEM(2,I,J/2+1,K))/4.
520            CONTINUE
               DO 530 J=3,NS-1,2
                  SEM(1,I,J,K)=(SEM(2,I,J/2,K)+3.*SEM(2,I,J/2+1,K))/4.
530            CONTINUE
               SEM(1,I,NS,K)=SEM(2,I,NS/2,K)
532         CONTINUE
540      CONTINUE
         DO 560 K=1,NP
            DO 550 J=1,NS
               DO 545 I=1,NR
                  SEM(2,I,J,K)=SEM(1,I,J,K)
545            CONTINUE
550         CONTINUE
560      CONTINUE
         WRITE(6,*) 'LARGEST DEPTH =',S(NS)
         WRITE(16,*) 'LARGEST DEPTH =',S(NS)
590   CONTINUE
C
600   CONTINUE
      IF (IWRIT.NE.0) WRITE(6,*) 'NUMBER OF ITERATIONS = ',ITM
      IF (IWRIT.NE.0) WRITE(6,*) 'BAND BENDING AT MIDPOINT = ',Pot0
      IF (IWRIT.NE.0) WRITE(16,*) 'NUMBER OF ITERATIONS = ',ITM
      IF (IWRIT.NE.0) WRITE(16,*) 'BAND BENDING AT MIDPOINT = ',Pot0
C
      RETURN
      END
C
C   FIND POTENTIAL ON CENTRAL AXIS
C
      FUNCTION PCENT(JJ,VAC,SEM,VSINT,NRDIM,NVDIM,NSDIM,NPDIM,NP)
      DIMENSION VAC(2,NRDIM,NVDIM,NPDIM),SEM(2,NRDIM,NSDIM,NPDIM),
     &VSINT(2,NRDIM,NPDIM)
      J=ABS(JJ)
      I=1
      SUM=0.
      IF (JJ.EQ.0) THEN
         DO 100 K=1,NP
            SUM=SUM+(9.*VSINT(1,I,K)-VSINT(1,I+1,K))/8.
100      CONTINUE
      ELSE
      IF (JJ.GT.0) THEN
         DO 200 K=1,NP
            SUM=SUM+(9.*VAC(1,I,J,K)-VAC(1,I+1,J,K))/8.
200      CONTINUE
      ELSE
         DO 300 K=1,NP
            SUM=SUM+(9.*SEM(1,I,J,K)-SEM(1,I+1,J,K))/8.
300      CONTINUE
      END IF
      END IF
      PCENT=SUM/FLOAT(NP)
      RETURN
      END
C
C   3-D SOLUTION (AZIMUTHAL SYMMETRY) FOR POTENTIAL
C
      SUBROUTINE ITER3(VAC,TIP,SEM,VSINT,R,DELR,DELV,DELXSI,S,DELS,BIAS,
     &DELR0,DELS0,DELP,DELETA,A,NRDIM,NVDIM,NSDIM,NPDIM,NR,NV,NS,NP,
     &EP,ITMAX,Pot0,IWRIT,ETAT,C,MIRROR,IERR,EPSIL,IBC)
C
      DIMENSION VAC(2,NRDIM,NVDIM,NPDIM),SEM(2,NRDIM,NSDIM,NPDIM),
     &VSINT(2,NRDIM,NPDIM),R(NRDIM),DELR(NRDIM),DELV(NRDIM),S(NSDIM),
     &DELS(NSDIM),DELXSI(NRDIM)
      LOGICAL TIP(NRDIM,NVDIM,NPDIM)
      COMMON/SEMSURFMIN/EPSIL1,X,Y,SJ,STEMP,DENOM,I,J,K,NR1,NS1,NP1
      EXTERNAL SEMMIN,SURFMIN
c   constant below is e/epsilon_0 in units of V cm, times 1.e-14 cm^2/nm^2
      DATA EEP/1.80943E-20/
C   PARAMETERS BELOW PASSED TO SEMMIN AND SURFMIN THROUGH SEMSURFMIN COMMON BLOCK
      NR1=NR
      NS1=NS
      NP1=NP
      EPSIL1=EPSIL
C
      PotSAV=0.
      C2=C*C
      C3=C2*C
      C2M1=C2-1.
      C2P1=C2+1.
      C2P2=C2+2.
      C2P6=C2+6.
      TC2P2=3.*C2+2.
      DO 500 ITER=1,ITMAX
C
C   DO FINITE-DIFFERENCE SOLN IN VACUUM
C
         DO 365 K=1,NP
            DO 360 I=1,NR
               X2M1=(R(I)/A)**2
               XSI=DSQRT(1.D0+X2M1)
               XSI2=XSI*XSI
               XSI3=XSI2*XSI
               XSI4=XSI3*XSI
               XSI5=XSI4*XSI
               DELXSI=R(I)*DELR(I)/(XSI*A**2)
               DO 350 J=1,NV-1
                  IF (TIP(I,J,K)) GO TO 360
                  ETA=J*DELETA
                  ETA2=ETA*ETA
                  ETA3=ETA2*ETA
                  ETA4=ETA3*ETA
                  ETA5=ETA4*ETA
                  OME2=1.-ETA**2
                  X2ME2=XSI**2-ETA**2
                  X2ME2C=XSI*(XSI+C)-ETA**2*(C*XSI+1.)
                  X2ME2C2=X2ME2C*X2ME2C
                  T1=X2M1*((XSI+C)**2-ETA**2*(XSI*C*2.+C2+1.))/X2ME2C
                  T2=OME2*X2ME2/X2ME2C
                  T3=X2ME2C/(X2M1*OME2)
                  T4=-C*ETA*X2M1*OME2/X2ME2C
                  T5=(C3+3.*C2*XSI+C*C2P2*XSI2+3.*C2*XSI3+4.*C*XSI4+
     &2.*XSI5+ETA4*(C3+TC2P2*XSI+C*C2P6*XSI2+3.*C2*XSI3)-
     &2.*ETA2*(C*C2M1+3.*C2*XSI+C*C2P6*XSI2+TC2P2*XSI3+C*XSI4))/X2ME2C2
                  T6=-ETA*(C2+4.*C*XSI+C2*XSI2+2.*XSI4+
     &ETA4*(2.+C2+4.*C*XSI+C2*XSI2)-
     &2.*ETA2*(C2+4.*C*XSI+C2P2*XSI2))/X2ME2C2
C
                  IF (J.EQ.1) THEN
                     IF (I.EQ.1) THEN
                        VACIM1JM1K=PCENT(J,VAC,SEM,VSINT,
     &                                      NRDIM,NVDIM,NSDIM,NPDIM,NP)
                     ELSE
                        VACIM1JM1K=VSINT(1,I-1,K)
                     END IF
                     VACIJM1K=VSINT(1,I,K)
                     IF (I.NE.NR) THEN
                        VACIP1JM1K=VSINT(1,I+1,K)
                     ELSE
                        VACIP1JM1K=IBC*VSINT(1,I,K)
                     END IF
                  ELSE
                     IF (I.EQ.1) THEN
                        VACIM1JM1K=PCENT(J-1,VAC,SEM,VSINT,
     &                                      NRDIM,NVDIM,NSDIM,NPDIM,NP)
                     ELSE
                        VACIM1JM1K=VAC(1,I-1,J-1,K)
                     END IF
                     VACIJM1K=VAC(1,I,J-1,K)
                     IF (I.NE.NR) THEN
                        VACIP1JM1K=VAC(1,I+1,J-1,K)
                     ELSE
                        VACIP1JM1K=IBC*VAC(1,I,J-1,K)
                     END IF
                  END IF
                  IF (I.EQ.1) THEN
                     VACIP1JK=VAC(1,I+1,J,K)
                     VACIM1JK=PCENT(J,VAC,SEM,VSINT,
     &                              NRDIM,NVDIM,NSDIM,NPDIM,NP)
                     VACIP1JP1K=VAC(1,I+1,J+1,K)
                     VACIM1JP1K=PCENT(J+1,VAC,SEM,VSINT,
     &                                NRDIM,NVDIM,NSDIM,NPDIM,NP)
                     DELXSII=DELXSI(I)
                     DELXSIIP1=DELXSI(I+1)
                     DELXSI2=DELXSI(I+1)+DELXSI(I)
                  ELSE
                  IF (I.EQ.NR) THEN
                     VACIP1JK=IBC*VAC(1,I,J,K)
                     VACIM1JK=VAC(1,I-1,J,K)
                     VACIP1JP1K=IBC*VAC(1,I,J+1,K)
                     VACIM1JP1K=VAC(1,I-1,J+1,K)
                     DELXSII=DELXSI(I)
                     DELXSIIP1=DELXSI(I)
                     DELXSI2=DELXSI(I)+DELXSI(I)
                  ELSE
                     VACIP1JK=VAC(1,I+1,J,K)
                     VACIM1JK=VAC(1,I-1,J,K)
                     VACIP1JP1K=VAC(1,I+1,J+1,K)
                     VACIM1JP1K=VAC(1,I-1,J+1,K)
                     DELXSII=DELXSI(I)
                     DELXSIIP1=DELXSI(I+1)
                     DELXSI2=DELXSI(I+1)+DELXSI(I)
                  END IF
                  END IF
                  IF (K.EQ.1) THEN
                     VACIJKP1=VAC(1,I,J,K+1)
                     IF (MIRROR.EQ.1) THEN
                        VACIJKM1=VAC(1,I,J,1)
                     ELSE
                        VACIJKM1=VAC(1,I,J,NP)
                     END IF
                  ELSE
                  IF (K.EQ.NP) THEN
                     IF (MIRROR.EQ.1) THEN
                        VACIJKP1=VAC(1,I,J,NP)
                     ELSE
                        VACIJKP1=VAC(1,I,J,1)
                     END IF
                     VACIJKM1=VAC(1,I,J,K-1)
                  ELSE
                     VACIJKP1=VAC(1,I,J,K+1)
                     VACIJKM1=VAC(1,I,J,K-1)
                  END IF
                  END IF
                  TEMP=
     &              T1*2.*(VACIP1JK/DELXSIIP1+VACIM1JK/DELXSII)/DELXSI2+
     &               T2*(VAC(1,I,J+1,K)+VACIJM1K)/DELETA**2+
     &               T3*(VACIJKP1+VACIJKM1)/DELP**2+
     &               T4*(VACIP1JP1K-VACIM1JP1K-
     &                   VACIP1JM1K+VACIM1JM1K)/(DELXSI2*DELETA)+
     &               T5*(VACIP1JK-VACIM1JK)/(DELXSI2)+
     &               T6*(VAC(1,I,J+1,K)-VACIJM1K)/(2.*DELETA)
C     
                  VAC(2,I,J,K)=
     &               TEMP/(2.*T1*(1/DELXSIIP1+1./DELXSII)/DELXSI2+
     &                     2.*T2/DELETA**2+2.*T3/DELP**2)
350            CONTINUE
360         CONTINUE
365      CONTINUE
         DO 385 K=1,NP
            DO 380 J=1,NV
               DO 370 I=1,NR
                  VAC(1,I,J,K)=VAC(2,I,J,K)
370            CONTINUE
380         CONTINUE
385      CONTINUE
C
C   DO FINITE-DIFFERENCE SOLN AT SURFACE
C
         DO 405 K=1,NP
            DO 400 I=1,NR
               X=R(I)*COS((K-0.5)*DELP)
               Y=R(I)*SIN((K-0.5)*DELP)
               SURFOLD=VSINT(1,I,K)
               IF (TIP(I,3,K)) GO TO 390
C   3RD ORDER SOLUTION IN SEMICONDUCTOR & VACUUM
               STEMP=
     &(3.*VAC(1,I,1,K)-(9./6.)*VAC(1,I,2,K)+
     &                   (1./3.)*VAC(1,I,3,K))/DELV(I)+
     &EPSIL*(3.75*SEM(1,I,1,K)-(5./6.)*SEM(1,I,2,K)+
     &                    0.15*SEM(1,I,3,K))/DELS0
               DENOM=((11./6.)/DELV(I)+(46./15.)*EPSIL/DELS0)
               GO TO 396
C   3RD ORDER SOLUTION IN SEMICONDUCTOR, 2ND ORDER IN VACUUM
390            IF (TIP(I,2,K)) GO TO 395
               STEMP=
     &(2.*VAC(1,I,1,K)-0.5*VAC(1,I,2,K))/DELV(I)+
     &EPSIL*(3.75*SEM(1,I,1,K)-(5./6.)*SEM(1,I,2,K)+
     &                    0.15*SEM(1,I,3,K))/DELS0
               DENOM=(1.5/DELV(I)+(46./15.)*EPSIL/DELS0)
               GO TO 396
C   3RD ORDER SOLUTION IN SEMICONDUCTOR, 1ST ORDER IN VACUUM
395            STEMP=
     &VAC(1,I,1,K)/DELV(I)+
     &EPSIL*(3.75*SEM(1,I,1,K)-(5./6.)*SEM(1,I,2,K)+
     &                    0.15*SEM(1,I,3,K))/DELS0
               DENOM=(1./DELV(I)+(46./15.)*EPSIL/DELS0)
396            CONTINUE
C
C      SOLVE NONLINEAR PROBLEM
C
               RHO=RHOSURF(VSINT(1,I,K),X,Y,I,K,NR,NP)
               TEMP=STEMP-RHO*EEP*1.E7
               SURFNEW=TEMP/DENOM
               DELSURF=AMAX1(1.E-6,ABS(BIAS)/1.E6)
               CALL GSECT(SURFMIN,SURFOLD,SURFNEW,DELSURF)
               VSINT(2,I,K)=(SURFOLD+SURFNEW)/2.
400         CONTINUE
405      CONTINUE
         DO 415 K=1,NP
            DO 410 I=1,NR
               VSINT(1,I,K)=VSINT(2,I,K)
410         CONTINUE
415      CONTINUE
         IF (MOD(ITER,100).EQ.0) THEN
            PotSAV2=PotSAV
            PotSAV=Pot0
            Pot0=PCENT(0,VAC,SEM,VSINT,NRDIM,NVDIM,NSDIM,NPDIM,NP)
            IF (IWRIT.NE.0) WRITE(6,*) 'ITER,Pot0 =',ITER,Pot0
            IF (IWRIT.NE.0) WRITE(16,*) 'ITER,Pot0 =',ITER,Pot0
         END IF
C
C   DO FINITE-DIFFERENCE SOLN IN BULK SEMICONDUCTOR
C
         DO 460 K=1,NP
            DO 450 J=1,NS
               DO 440 I=1,NR
                  SEMOLD=SEM(1,I,J,K)
                  ENER=EF-SEM(1,I,J,K)
                  RSAV=R(I)
                  X=R(I)*COS((K-0.5)*DELP)
                  Y=R(I)*SIN((K-0.5)*DELP)
                  IF (I.EQ.1) THEN
                     SEMIM1JK=PCENT(-J,VAC,SEM,VSINT,
     &                               NRDIM,NVDIM,NSDIM,NPDIM,NP)
                     SEMIP1JK=SEM(1,I+1,J,K)
                     DELR2=DELR(I+1)+DELR(I)
                     DELRIP1=DELR(I+1)
                     DELRI=DELR(I)
                  ELSE
                  IF (I.EQ.NR) THEN
                     SEMIM1JK=SEM(1,I-1,J,K)
                     SEMIP1JK=IBC*SEM(1,I,J,K)
                     DELR2=DELR(I)+DELR(I)
                     DELRIP1=DELR(I)
                     DELRI=DELR(I)
                  ELSE
                     SEMIM1JK=SEM(1,I-1,J,K)
                     SEMIP1JK=SEM(1,I+1,J,K)
                     DELR2=DELR(I+1)+DELR(I)
                     DELRIP1=DELR(I+1)
                     DELRI=DELR(I)
                  END IF
                  END IF
                  IF (J.EQ.1) THEN
                     SEMIJP1K=SEM(1,I,J+1,K)
                     SEMIJM1K=VSINT(1,I,K)
                     DELS2=DELS(J+1)+DELS(J)
                     DELSJP1=DELS(J+1)
                     DELSJ=DELS(J)
                  ELSE
                  IF (J.EQ.NS) THEN
                     SEMIJP1K=IBC*SEM(1,I,J,K)
                     SEMIJM1K=SEM(1,I,J-1,K)
                     DELS2=DELS(J)+DELS(J)
                     DELSJP1=DELS(J)
                     DELSJ=DELS(J)
                  ELSE
                     SEMIJP1K=SEM(1,I,J+1,K)
                     SEMIJM1K=SEM(1,I,J-1,K)
                     DELS2=DELS(J+1)+DELS(J)
                     DELSJP1=DELS(J+1)
                     DELSJ=DELS(J)
                  END IF
                  END IF
                  IF (K.EQ.1) THEN
                     SEMIJKP1=SEM(1,I,J,K+1)
                     IF (MIRROR.EQ.1) THEN
                        SEMIJKM1=SEM(1,I,J,1)
                     ELSE
                        SEMIJKM1=SEM(1,I,J,NP)
                     END IF
                  ELSE
                  IF (K.EQ.NP) THEN
                     IF (MIRROR.EQ.1) THEN
                        SEMIJKP1=SEM(1,I,J,NP)
                     ELSE
                        SEMIJKP1=SEM(1,I,J,1)
                     END IF
                     SEMIJKM1=SEM(1,I,J,K-1)
                  ELSE
                     SEMIJKP1=SEM(1,I,J,K+1)
                     SEMIJKM1=SEM(1,I,J,K-1)
                  END IF
                  END IF
                  STEMP=
     &               2.*(SEMIP1JK/DELRIP1+SEMIM1JK/DELRI)/DELR2+
     &               2.*(SEMIJP1K/DELSJP1+SEMIJM1K/DELSJ)/DELS2+
     &               (SEMIP1JK-SEMIM1JK)/(R(I)*DELR2)+
     &               (SEMIJKP1+SEMIJKM1)/(R(I)**2*DELP**2)
C
C      SOLVE NONLINEAR PROBLEM
C
                  RHO=RHOBULK(SEM(1,I,J,K),X,Y,S(J),I,J,K,NR,NS,NP)
                  TEMP=STEMP-RHO*EEP/EPSIL
                  DENOM=2.*(1./DELRIP1+1./DELRI)/DELR2+
     &               2.*(1./DELSJP1+1./DELSJ)/DELS2+2./(R(I)**2*DELP**2)
                  SEMNEW=TEMP/DENOM
                  DELSEM=AMAX1(1.E-6,ABS(BIAS)/1.E6)
C   PARAMETER BELOW PASSED TO SEMMIN THROUGH SEMSURFMIN COMMON BLOCK
                  SJ=S(J)
                  CALL GSECT(SEMMIN,SEMOLD,SEMNEW,DELSEM)
                  SEM(2,I,J,K)=(SEMOLD+SEMNEW)/2.
440            CONTINUE
450         CONTINUE
460      CONTINUE
         DO 485 K=1,NP
            DO 480 J=1,NS
               DO 470 I=1,NR
                  SEM(1,I,J,K)=SEM(2,I,J,K)
470            CONTINUE
480         CONTINUE   
485      CONTINUE
C
C   OUTPUT PROFILE (EACH ITERATION; FOR DEBUGGING)
C
         IF (IWRIT.GE.3) THEN
            DO 487 J=NV,1,-1
              WRITE(11,*) -J*DELV(1),VAC(1,1,J,1),VAC(1,NR,J,1)
487         CONTINUE
            WRITE(11,*) 0.,VSINT(1,1,1),VSINT(1,NR,1)
            DO 490 J=1,NS
                WRITE(11,*) S(J),SEM(1,1,J,1),SEM(1,NR,J,1)
490         CONTINUE
            WRITE(12,*) -R(1),VSINT(1,1,1)
            DO 495 I=1,NR
                WRITE(12,*) R(I),VSINT(1,I,1)
495         CONTINUE
            CLOSE(11)
            CLOSE(12)
         END IF
         IF ((MOD(ITER,100).EQ.0.AND.ABS(Pot0-PotSAV).LT.EP.AND.
     &        ABS(PotSAV-PotSAV2).LT.2.*EP)) GO TO 510
500   CONTINUE
510   ITMAX=ITER
      Pot0=PCENT(0,VAC,SEM,VSINT,NRDIM,NVDIM,NSDIM,NPDIM,NP)
      RETURN
      END

