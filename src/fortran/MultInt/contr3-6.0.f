C   ******************** CONTR3 *********************
C
C   DRAW CONTOURS OF THE POTENTIAL DISTRIBUTIONS (3D CASE)
C
C   VERSION 6.0 - FEB/11 - SAME AS PRIOR QWCONTR (WITHOUT NUMBERS)
C
      SUBROUTINE CONTR3(ETA1,VAC,TIP,SEM,VSINT,R,S,DELV,NRDIM,NVDIM,
     &NSDIM,NPDIM,NR,NV,NS,NP,NUMC,DELPOT,MIRROR,KPLOT1,KPLOT2)
C
C   NUMC=NUMBER OF CONTOURS
C   DELPOT=SPACING OF POTENTIAL CONTOURS (IF 0 USE MAX_P-MIN_P/(NUMC+1))
C
      DIMENSION VAC(2,NRDIM,NVDIM,NPDIM),SEM(2,NRDIM,NSDIM,NPDIM),
     &VSINT(2,NRDIM,NPDIM),R(NRDIM),S(NSDIM),DELV(NRDIM)
      LOGICAL TIP (NRDIM,NVDIM,NPDIM),KDONE(1000)
C
C   DRAW TIP
C
      NRSKIP=NR/500
      IF (NRSKIP.EQ.0) NRSKIP=1
      DO 300 I=-NR,NR,NRSKIP
         IF (I.EQ.0) GO TO 300
         IF (I.GT.0) THEN
            II=I
            RSAV=R(II)
         ELSE
            II=-I
            RSAV=-R(II)
         END IF
         DO 250 J=1,NV
            IF (.NOT.TIP(II,J,1)) GO TO 250
            WRITE(20,*) RSAV*SQRT(1.-(ETA1*J/NV)**2),-J*DELV(II)
            GO TO 300
250      CONTINUE
300   CONTINUE
C
C   SEARCH FOR MIN, MAX POINTS IN POTENTIAL
C
      PMIN=1.E10
      PMAX=-1.E10
      DO 400 I=1,NR
         DO 390 K=1,NP,NP-1
         DO 330 J=1,NV
            IF (PMIN.LT.VAC(1,I,J,K)) GO TO 320
               PMIN=VAC(1,I,J,K)
320         IF (PMAX.GT.VAC(1,I,J,K)) GO TO 330
               PMAX=VAC(1,I,J,K)
330      CONTINUE
         IF (PMIN.LT.VSINT(1,I,K)) GO TO 340
            PMIN=VSINT(1,I,K)
340      IF (PMAX.GT.VSINT(1,I,K)) GO TO 350
            PMAX=VSINT(1,I,K)
350      DO 370 J=1,NS
            IF (PMIN.LT.SEM(1,I,J,K)) GO TO 360
               PMIN=SEM(1,I,J,K)
360         IF (PMAX.GT.SEM(1,I,J,K)) GO TO 370
               PMAX=SEM(1,I,J,K)
370      CONTINUE
390      CONTINUE
400   CONTINUE
      WRITE(6,*) 'MIN, MAX POTENTIAL VALUES =',PMIN,PMAX
      WRITE(16,*) 'MIN, MAX POTENTIAL VALUES =',PMIN,PMAX
C
C   DRAW CONTOURS
C
410   IF (DELPOT.NE.0.) GO TO 450
      DELPOT=(PMAX-PMIN)/(NUMC+1)
      WRITE(6,*) 'CONTOUR SPACING =',DELPOT
      WRITE(16,*) 'CONTOUR SPACING =',DELPOT
C
450   DO 600 I=-NR,NR,NRSKIP
         IF (I.EQ.0) GO TO 600
         IF (I.GT.0) THEN
            II=I
            RSAV=R(II)
            KP=KPLOT1
         ELSE
            II=-I
            RSAV=-R(II)
            KP=KPLOT2
         END IF
         DO 500 K=1,NUMC
            KDONE(K)=.FALSE.
500      CONTINUE
         DO 530 K=1,NUMC
         DO 520 J=NS,1,-1
            P=K*DELPOT+PMIN
            IF (J.EQ.1) GO TO 510
            IF ((SEM(1,II,J,KP).GE.P.AND.SEM(1,II,J-1,KP).LE.P).OR.
     &          (SEM(1,II,J,KP).LE.P.AND.SEM(1,II,J-1,KP).GE.P))GOTO 515
            GO TO 520
510         IF ((SEM(1,II,J,KP).GE.P.AND.VSINT(1,II,KP).LE.P).OR.
     &          (SEM(1,II,J,KP).LE.P.AND.VSINT(1,II,KP).GE.P)) GO TO 515
            GO TO 520
515         WRITE(20+K,*) RSAV,S(J)
            KDONE(K)=.TRUE.
            GO TO 530
520      CONTINUE
530      CONTINUE
         DO 540 K=1,NUMC
         P=K*DELPOT+PMIN
         IF ((VSINT(1,II,KP).GE.P.AND.VAC(1,II,1,KP).LE.P).OR.
     &       (VSINT(1,II,KP).LE.P.AND.VAC(1,II,1,KP).GE.P)) GO TO 535
         GO TO 540
535      IF (.NOT.KDONE(K)) WRITE(20+K,*) RSAV,0.
         KDONE(K)=.TRUE.
540      CONTINUE
         DO 570 K=1,NUMC
         DO 560 J=1,NV-1
            IF (TIP(II,J,KP)) GO TO 560
            P=K*DELPOT+PMIN
            IF ((VAC(1,II,J,KP).GE.P.AND.VAC(1,II,J+1,KP).LE.P).OR.
     &          (VAC(1,II,J,KP).LE.P.AND.VAC(1,II,J+1,KP).GE.P))GOTO 550
            GO TO 560
550         IF (.NOT.KDONE(K))
     &         WRITE(20+K,*) RSAV*SQRT(1.-(ETA1*J/NV)**2),-J*DELV(II)
            KDONE(K)=.TRUE.
            GO TO 570
560      CONTINUE
570      CONTINUE
600   CONTINUE
C
C   NOTE POSSIBLE ERRORS BELOW:  GO TO -NR AND ADD STATEMENTS BELOW THAT !!
C
      DO 700 I=NR,-NR,-NRSKIP
         IF (I.EQ.0) GO TO 700
         IF (I.GT.0) THEN
            II=I
            RSAV=R(II)
            KP=KPLOT1
         ELSE
            II=-I
            RSAV=-R(II)
            KP=KPLOT2
         END IF
         DO 605 K=1,NUMC
            KDONE(K)=.FALSE.
605      CONTINUE
         DO 630 K=1,NUMC
         DO 620 J=NV-1,1,-1
            IF (TIP(II,J,KP)) GO TO 620
            P=K*DELPOT+PMIN
            IF ((VAC(1,II,J,KP).GE.P.AND.VAC(1,II,J+1,KP).LE.P).OR.
     &          (VAC(1,II,J,KP).LE.P.AND.VAC(1,II,J+1,KP).GE.P))GOTO 610
            GO TO 620
610         IF (.NOT.KDONE(K)) 
     &         WRITE(20+K,*) RSAV*SQRT(1.-(ETA1*J/NV)**2),-J*DELV(II)
           KDONE(K)=.TRUE.
            GO TO 630
620      CONTINUE
630      CONTINUE
         DO 640 K=1,NUMC
         P=K*DELPOT+PMIN
         IF ((VSINT(1,II,KP).GE.P.AND.VAC(1,II,1,KP).LE.P).OR.
     &       (VSINT(1,II,KP).LE.P.AND.VAC(1,II,1,KP).GE.P)) GO TO 635
         GO TO 640
635      IF (.NOT.KDONE(K)) WRITE(20+K,*) RSAV,0.
         KDONE(K)=.TRUE.
640      CONTINUE
         DO 670 K=1,NUMC
         DO 660 J=1,NS
            P=K*DELPOT+PMIN
            IF (J.EQ.1) GO TO 645
            IF ((SEM(1,II,J,KP).GE.P.AND.SEM(1,II,J-1,KP).LE.P).OR.
     &          (SEM(1,II,J,KP).LE.P.AND.SEM(1,II,J-1,KP).GE.P))GOTO 650
            GO TO 660
645         IF ((SEM(1,II,J,KP).GE.P.AND.VSINT(1,II,KP).LE.P).OR.
     &          (SEM(1,II,J,KP).LE.P.AND.VSINT(1,II,KP).GE.P)) GO TO 650
            GO TO 660
650         IF (.NOT.KDONE(K)) WRITE(20+K,*) RSAV,S(J)
            KDONE(K)=.TRUE.
            GO TO 670
660      CONTINUE
670      CONTINUE
700   CONTINUE
C
      RETURN
      END