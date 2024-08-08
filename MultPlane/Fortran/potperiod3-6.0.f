C   ******************** POTPERIOD ************************
C
C   CONSTRUCT A PERIODIC POTENTIAL WITH EQUAL GRID SPACING
C   FOR USE IN PLANE WAVE COMPUTATIONS, STARTING WITH CIRCULARLY
C   SYMMETRIC COORDINATES
C   ASSUME MIRROR PLANE IN X, Y AND Z DIRECTIONS
C
      SUBROUTINE POTPERIOD(SEP,NX,NY,NZ,XDEL,YDEL,ZDEL,
     &VAC,TIP,SEM,VSINT,PVAC,PSEM,R,S,DELV,DELR,DELS,DELP,NRDIM,
     &NVDIM,NSDIM,NPDIM,NXDIM,NXDIM2,NYDIM,NZDIM,NR,NV,NS,NP,IWRIT,IERR)
C
      DIMENSION VAC(2,NRDIM,NVDIM,NPDIM),SEM(2,NRDIM,NSDIM,NPDIM),
     &VSINT(2,NRDIM,NPDIM),R(NRDIM),DELV(NRDIM),S(NSDIM),
     &PVAC(NXDIM2,NYDIM,NVDIM),PSEM(NXDIM2,NYDIM,NZDIM)
      LOGICAL TIP(NRDIM,NVDIM,NPDIM)
      DATA PI/3.141592654/
C
      DELV0=SEP/NV
c      WRITE(6,*) 'NUMBER OF Z POINTS IN VACUUM, SPACING =',NV,DELV0
c      WRITE(6,*) 'NUMBER OF X POINTS IN SEMICOND, SPACING =',NX,XDEL
c      WRITE(6,*) 'NUMBER OF Y POINTS IN SEMICOND, SPACING =',NY,YDEL
c      WRITE(6,*) 'NUMBER OF Z POINTS IN SEMICOND, SPACING =',NZ,ZDEL
      IF (NX.GT.NXDIM.OR.NY.GT.NYDIM.OR.NZ.GT.NZDIM) THEN
         WRITE(6,*) '*** ERROR - INSUFFICIENT STORAGE FOR PERIODIC',
     &              ' POTENTIAL'
         RETURN
      END IF
C
C   POTENTIAL INTO SEMICONDUCTOR (AND ON SURFACE)
C
      DO 280 K=1,NZ
         Z=(K-1)*ZDEL
         IF (Z.LT.S(1)) THEN
            FZ=Z/S(1)
            KK=0
         ELSE
            DO 210 KK=1,NS
               IF (S(KK).LE.Z.AND.S(KK+1).GT.Z) GO TO 220
210         CONTINUE
            WRITE(6,*) '*** ERROR - INTERPOLATED S NOT FOUND'
            write(6,*) z,ns,kk
            WRITE(6,*) 'press ENTER to continue'
            READ(5,*)
220         FZ=(Z-S(KK))/(S(KK+1)-S(KK))
         END IF
         DO 270 J=1,NY
            Y=(J-1)*YDEL
            DO 260 I=1,2*NX-2
               X=(I-NX+1)*XDEL
               PHI=ATAN2(Y,X)
               IP=NINT((PHI/DELP)+0.5)
               IF (IP.LE.0) IP=1
               IF (IP.GT.NP) IP=NP
               R0=SQRT(X**2+Y**2)
               IF (I.EQ.(NX-1).AND.J.EQ.1) THEN
                  IF (KK.EQ.0) THEN
                     PSEM(NX-1,1,K)=((9.*(VSINT(1,1,1)+VSINT(1,1,NP))
     &                 -(VSINT(1,2,1)+VSINT(1,2,NP)))/16.)*(1.-FZ)
     &                 +((9.*(SEM(1,1,1,1)+SEM(1,1,1,NP))
     &                 -(SEM(1,2,1,1)+SEM(1,2,1,NP)))/16.)*FZ
                  ELSE
                     PSEM(NX-1,1,K)=((9.*(SEM(1,1,KK,1)+SEM(1,1,KK,NP))
     &                 -(SEM(1,2,KK,1)+SEM(1,2,KK,NP)))/16.)*(1.-FZ)
     &                 +((9.*(SEM(1,1,KK+1,1)+SEM(1,1,KK+1,NP))
     &                 -(SEM(1,2,KK+1,1)+SEM(1,2,KK+1,NP)))/16.)*FZ
                  END IF
               ELSE
                 IF (R0.LE.R(1)) THEN
                  F=R0/R(1)
                  IF (KK.EQ.0) THEN
            PSEM(I,J,K)=(((9.*(VSINT(1,1,1)+VSINT(1,1,NP))
     &                 -(VSINT(1,2,1)+VSINT(1,2,NP)))/16.)*(1.-F)
     &                  +VSINT(1,1,IP)*F)*(1.-FZ)
     &                +(((9.*(SEM(1,1,1,1)+SEM(1,1,1,NP))
     &                 -(SEM(1,2,1,1)+SEM(1,2,1,NP)))/16.)*(1.-F)
     &                  +SEM(1,1,1,IP)*F)*FZ
                  ELSE
            PSEM(I,J,K)=(((9.*(SEM(1,1,KK,1)+SEM(1,1,KK,NP))
     &                 -(SEM(1,2,KK,1)+SEM(1,2,KK,NP)))/16.)*(1.-F)
     &                  +SEM(1,1,KK,IP)*F)*(1.-FZ)
     &                +(((9.*(SEM(1,1,KK+1,1)+SEM(1,1,KK+1,NP))
     &                 -(SEM(1,2,KK+1,1)+SEM(1,2,KK+1,NP)))/16.)*(1.-F)
     &                  +SEM(1,1,KK+1,IP)*F)*FZ
                  END IF
                 ELSE
                  DO 240 II=1,NR
                     IF (R(II).LE.R0.AND.R(II+1).GT.R0) GO TO 250
240               CONTINUE
                  WRITE(6,*) '*** ERROR - INTERPOLATED R NOT FOUND'
                  WRITE(6,*) 'II,NR,R0,R=',II,NR,R0,(R(III),III=1,10)
                  WRITE(6,*) 'press ENTER to continue'
                  READ(5,*)
250               F=(R0-R(II))/(R(II+1)-R(II))
                  IF (KK.EQ.0) THEN
               PSEM(I,J,K)=(VSINT(1,II,IP)*(1.-F)+VSINT(1,II+1,IP)*F)*
     &           (1.-FZ)+(SEM(1,II,1,IP)*(1.-F)+SEM(1,II+1,1,IP)*F)*FZ
                  ELSE
               PSEM(I,J,K)=(SEM(1,II,KK,IP)*(1.-F)+SEM(1,II+1,KK,IP)*F)*
     &       (1.-FZ)+(SEM(1,II,KK+1,IP)*(1.-F)+SEM(1,II+1,KK+1,IP)*F)*FZ
                  END IF
                 END IF
               END IF
260         CONTINUE
270      CONTINUE
280   CONTINUE
C
C   POTENTIAL INTO VACUUM
C
      DO 380 K=1,NV
         Z=K*DELV0
         DO 370 J=1,NY
            Y=(J-1)*YDEL
            DO 360 I=1,2*NX-2
               X=(I-NX+1)*XDEL
               PHI=ATAN2(Y,X)
               IP=NINT((PHI/DELP)+0.5)
               IF (IP.LE.0) IP=1
               IF (IP.GT.NP) IP=NP
               R0=SQRT(X**2+Y**2)
               IF (I.EQ.(NX-1).AND.J.EQ.1) THEN
                  PVAC(NX-1,1,K)=(9.*(VAC(1,1,K,1)+VAC(1,1,K,NP))
     &                          -(VAC(1,2,K,1)+VAC(1,2,K,NP)))/16.
               ELSE
                IF (R0.LE.R(1)) THEN
                  F=R0/R(1)
                  PVAC(I,J,K)=((9.*(VAC(1,1,K,1)+VAC(1,1,K,NP))
     &                      -(VAC(1,2,K,1)+VAC(1,2,K,NP)))/16.)*(1.-F)
     &                      +VAC(1,1,K,IP)*F
                ELSE
                  DO 310 II=1,NR
                     IF (R(II).LE.R0.AND.R(II+1).GT.R0) GO TO 320
310               CONTINUE
                  WRITE(6,*) '*** ERROR - INTERPOLATED R NOT FOUND'
                  WRITE(6,*) 'II,NR,R0,R=',II,NR,R0,(R(III),III=1,10)
                  WRITE(6,*) 'press ENTER to continue'
                  READ(5,*)
320               F=(R0-R(II))/(R(II+1)-R(II))
                  DELV1=DELV(II)*(1.-F)+DELV(II+1)*F
                  IF (Z.LT.DELV1) THEN
                     FZ=Z/DELV1
                     KK=0
                  ELSE
                     DO 330 KK=1,NV
                      IF (KK*DELV1.LE.Z.AND.(KK+1)*DELV1.GT.Z) GO TO 340
330                  CONTINUE
                     WRITE(6,*) '*** ERROR - INTERPOLATED Z NOT FOUND'
                     WRITE(6,*) 'press ENTER to continue'
                     READ(5,*)
340                  FZ=(Z-KK*DELV1)/DELV1
                  END IF
                  IF (KK.EQ.0) THEN
              PVAC(I,J,K)=(VSINT(1,II,IP)*(1.-F)+VSINT(1,II+1,IP)*F)*
     &           (1.-FZ)+(VAC(1,II,1,IP)*(1.-F)+VAC(1,II+1,1,IP)*F)*FZ
                  ELSE
              PVAC(I,J,K)=(VAC(1,II,KK,IP)*(1.-F)+VAC(1,II+1,KK,IP)*F)*
     &       (1.-FZ)+(VAC(1,II,KK+1,IP)*(1.-F)+VAC(1,II+1,KK+1,IP)*F)*FZ
                  END IF
                END IF
               END IF
360         CONTINUE
370      CONTINUE
380   CONTINUE
C
      RETURN
      END
