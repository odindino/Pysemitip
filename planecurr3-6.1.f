C   ******************** PLANECURR3 ************************
C
C   COMPUTE TUNNEL CURRENTS USING PLANE WAVES FOR QUANTUM DOT
C   ASSUME MIRROR PLANE IN X DIRECTION FOR POTENTIAL
C
C   VERSION 6.0 - MAY/11
C           6.1 - NOV/12, RENAME CBEDGE, VBEDGE TO DELCBEDGE, DELVBEDGE
C
      SUBROUTINE planecurr(SEP,VAC,TIP,SEM,VSINT,R,S,DELV,DELR,DELS,
     &DELP,NRDIM,NVDIM,NSDIM,NPDIM,NR,NV,NS,NP,NXDIM,NXDIM2,NYDIM,
     &NZDIM,PVAC,PSEM,PSURF,VACWID,CHI,EFTIP,EGAP,AVBL,AVBH,AVBSO,ACB,
     &ESO,BIAS,DELPHI,PHI0,TK,EF,EMAX,CURR,CURRV,CURRC,IWRIT,IERR,
     &NBARR2,NVDIM2,BARRPROF,NEIGENDIM,ZVACDEL,ICOMP)
C
      DIMENSION VAC(2,NRDIM,NVDIM,NPDIM),SEM(2,NRDIM,NSDIM,NPDIM),
     &VSINT(2,NRDIM,NPDIM),R(NRDIM),DELV(NRDIM),S(NSDIM),AVB(3),
     &PVAC(NXDIM2,NYDIM,NVDIM),PSEM(NXDIM2,NYDIM,NZDIM),EMAX(4),
     &PSURF(NXDIM,NYDIM),EPSI2(4,NEIGENDIM),BARRPROF(NVDIM2)
      LOGICAL TIP(NRDIM,NVDIM,NPDIM)
c   value of C below is 2m/hbar^2 in units of 1/(eV nm^2)
      DATA C/26.254/RQUANT/12900./
      PI=4.*ATAN(1.)
      CURRL=0.
      CURRH=0.
      CURRSO=0.
      CURRC=0.
C
C   BARR BELOW IS AVERAGE BARRIER; BARR1 IS ADJUSTED VALUE TO USE IN 
C   FOR BASIS SET (PHI0 IS ADDED IN FULL SOLN OF SCHRODINGER EQN);
C   BARR2 HAS ZERO LEVEL SHIFTED, TO ACCOUNT FOR ZERO IN ENERGIES
C
C   VALENCE BAND; LIGHT HOLES
C
      IF (ICOMP.EQ.0) GO TO 400
      IF (IWRIT.GT.0) THEN
      WRITE(6,*) '*********** LIGHT-HOLE VALENCE BAND *************'
      WRITE(16,*) '*********** LIGHT-HOLE VALENCE BAND *************'
      END IF
      EFFM=AVBL
      BARR=CHI+EGAP+DELPHI/2.+BIAS/2.+PHI0/2.
      BARR1=BARR-PHI0
      CALL GETstates(-1,SEP,VAC,TIP,SEM,VSINT,R,S,DELV,DELR,DELS,DELP,
     &NRDIM,NVDIM,NSDIM,NPDIM,NR,NV,NS,NP,NXDIM,NXDIM2,NYDIM,NZDIM,
     &PVAC,PSEM,PSURF,VACWID,BARR1,EFFM,EMAX(1),IWRIT,IERR,EPSI2,
     &NSTATES,NEIGENDIM,ZVACDEL)
      E0=0.
      BARR2=BARR+E0
      CALL GETcurr(-1,BIAS,SEP,E0,BARR2,TK,TK,EF,EPSI2,NSTATES,EFTIP,
     &NBARR2,NVDIM2,BARRPROF,NEIGENDIM,CURRL)
      IF (IWRIT.GT.0) THEN
      WRITE(6,*) 'LIGHT-HOLE CURRENT =',CURRL
      WRITE(16,*) 'LIGHT-HOLE CURRENT =',CURRL
      END IF
      IF (IWRIT.EQ.3) write(51,*) bias,CURRL
c      go to 400
C
C   VALENCE BAND; HEAVY HOLES
C
      IF (IWRIT.GT.0) THEN
      WRITE(6,*) '*********** HEAVY-HOLE VALENCE BAND *************'
      WRITE(16,*) '*********** HEAVY-HOLE VALENCE BAND *************'
      END IF
      EFFM=AVBH
      BARR=CHI+EGAP+DELPHI/2.+BIAS/2.+PHI0/2.
      BARR1=BARR-PHI0
      CALL GETstates(-1,SEP,VAC,TIP,SEM,VSINT,R,S,DELV,DELR,DELS,DELP,
     &NRDIM,NVDIM,NSDIM,NPDIM,NR,NV,NS,NP,NXDIM,NXDIM2,NYDIM,NZDIM,
     &PVAC,PSEM,PSURF,VACWID,BARR1,EFFM,EMAX(2),IWRIT,IERR,EPSI2,
     &NSTATES,NEIGENDIM,ZVACDEL)
      E0=0.
      BARR2=BARR+E0
      CALL GETcurr(-1,BIAS,SEP,E0,BARR2,TK,TK,EF,EPSI2,NSTATES,EFTIP,
     &NBARR2,NVDIM2,BARRPROF,NEIGENDIM,CURRH)
      IF (IWRIT.GT.0) THEN
      WRITE(6,*) 'HEAVY-HOLE CURRENT =',CURRH
      WRITE(16,*) 'HEAVY-HOLE CURRENT =',CURRH
      END IF
      IF (IWRIT.EQ.3) write(52,*) bias,CURRH
C
C   VALENCE BAND; SPLIT-OFF BAND
C
      IF (IWRIT.GT.0) THEN
      WRITE(6,*) '*********** SPLIT-OFF VALENCE BAND *************'
      WRITE(16,*) '*********** SPLIT-OFF VALENCE BAND *************'
      END IF
      EFFM=AVBSO
      BARR=CHI+EGAP+DELPHI/2.+BIAS/2.+PHI0/2.
      BARR1=BARR-PHI0
      CALL GETstates(-1,SEP,VAC,TIP,SEM,VSINT,R,S,DELV,DELR,DELS,DELP,
     &NRDIM,NVDIM,NSDIM,NPDIM,NR,NV,NS,NP,NXDIM,NXDIM2,NYDIM,NZDIM,
     &PVAC,PSEM,PSURF,VACWID,BARR1,EFFM,EMAX(3),IWRIT,IERR,EPSI2,
     &NSTATES,NEIGENDIM,ZVACDEL)
      E0=-ESO
      BARR2=BARR+E0
      CALL GETcurr(-1,BIAS,SEP,E0,BARR2,TK,TK,EF,EPSI2,NSTATES,EFTIP,
     &NBARR2,NVDIM2,BARRPROF,NEIGENDIM,CURRSO)
      IF (IWRIT.GT.0) THEN
      WRITE(6,*) 'SPLIT-OFF CURRENT =',CURRSO
      WRITE(16,*) 'SPLIT-OFF CURRENT =',CURRSO
      CURRV=CURRL+CURRH+CURRSO
      WRITE(6,*) 'TOTAL HOLE CURRENT =',CURRV
      WRITE(16,*) 'TOTAL HOLE CURRENT =',CURRV
      END IF
      IF (IWRIT.EQ.3) write(53,*) bias,CURRSO
C
C   CONDUCTION BAND
C
400   IF (ICOMP.EQ.-1) GO TO 500    
      IF (IWRIT.GT.0) THEN
      WRITE(6,*) '*********** CONDUCTION BAND *************'
      WRITE(16,*) '*********** CONDUCTION BAND *************'
      END IF
      EFFM=ACB
      BARR=CHI+DELPHI/2.+BIAS/2.+PHI0/2.
      BARR1=BARR-PHI0
      CALL GETstates(1,SEP,VAC,TIP,SEM,VSINT,R,S,DELV,DELR,DELS,DELP,
     &NRDIM,NVDIM,NSDIM,NPDIM,NR,NV,NS,NP,NXDIM,NXDIM2,NYDIM,NZDIM,
     &PVAC,PSEM,PSURF,VACWID,BARR1,EFFM,EMAX(4),IWRIT,IERR,EPSI2,
     &NSTATES,NEIGENDIM,ZVACDEL)
      E0=EGAP
      BARR2=BARR+E0
      CALL GETcurr(1,BIAS,SEP,E0,BARR2,TK,TK,EF,EPSI2,NSTATES,EFTIP,
     &NBARR2,NVDIM2,BARRPROF,NEIGENDIM,CURRC)
      IF (IWRIT.GT.0) THEN
      WRITE(6,*) 'ELECTRON CURRENT =',CURRC
      WRITE(16,*) 'ELECTRON CURRENT =',CURRC
      END IF
C      
500   CURR=CURRV+CURRC
      RETURN
      END
C
C  ************** GET EFFECTIVE MASS STATES ******************
C
      SUBROUTINE GETstates(IBAND,SEP,VAC,TIP,SEM,VSINT,R,S,DELV,
     &DELR,DELS,DELP,NRDIM,NVDIM,NSDIM,NPDIM,NR,NV,NS,NP,NXDIM,NXDIM2,
     &NYDIM,NZDIM,PVAC,PSEM,PSURF,VACWID,BARR,EFFM,EMAX,IWRIT,IERR,
     &EPSI2,NSTATES,NEIGENDIM,ZVACDEL)
C
C   IBAND=-1 FOR VB, 1 FOR CB
C
C   max value of nkz is nz-2 (can use smaller value)
C   max value of nkr is nx-1 (can use smaller value)
C
      parameter(nkx=6,nky=6,nkz=6,nkxm1=nkx-1,nkx2m1=2*nkx-1,
     &nk3=nkx2m1*nky*nkz,npack=nk3*(nk3+1)/2,nauxmx=4*nk3,NZVACDIM=100)
C
      real*8 cmat(NK3,NK3),w(NK3),zmat(NK3,NK3),aux1(nk3),aux2(nk3),
     &ener(1000),amp1(1000),amp2(1000),DEFFM,ELZ1,ELZ2,V0
      DIMENSION VAC(2,NRDIM,NVDIM,NPDIM),SEM(2,NRDIM,NSDIM,NPDIM),
     &VSINT(2,NRDIM,NPDIM),R(NRDIM),DELV(NRDIM),S(NSDIM),
     &PVAC(NXDIM2,NYDIM,NVDIM),PSEM(NXDIM2,NYDIM,NZDIM),
     &COSX(NKX,NXDIM2),SINX(NKXM1,NXDIM2),COSY(NKY,NYDIM),
     &HARM(NKX,NKY,NKZ,NZDIM),TUNN(NKX,NKY,NKZ,NZVACDIM),
     &wq(NKX,NKY,NKZ),wkap(NKX,NKY,NKZ),PSURF(NXDIM,NYDIM),
     &EPSI2(4,NEIGENDIM),AM1(NKX,NKY,NKZ),AM2(NKX,NKY,NKZ),
     &plot1(nxdim2,nydim),plot2(nxdim2,nydim),
     &plot3(nxdim2,nydim),plot4(nxdim2,nydim),plot5(nxdim2,nydim),
     &plot6(nxdim2,nydim),plot7(nxdim2,nydim),plot8(nxdim2,nydim),
     &plot9(nxdim2,nydim),img(512,512)
      LOGICAL TIP(NRDIM,NVDIM,NPDIM)
      DOUBLE PRECISION SUM1,SUM2,SUM3,SUM4,SUM5,WKX(NKX),
     &WKY(NKY)
c   value of C below is 2m/hbar^2 in units of 1/(eV nm^2)
      DATA C/26.254/RQUANT/12900./
      PI=4.*ATAN(1.)
C      
C   set IPOT=1 to include computed potential, 0 to neglect it
C
      IPOT=1
C
C   set IVAC=1 to include vacuum, 0 to neglect it (and thus have full dot)
C
      IVAC=1
C      
      DEFFM=EFFM
      IF (IWRIT.GT.0) THEN
      WRITE(6,*) 'ENTERING GETstates'
      WRITE(6,*) 'effm =',effm
      WRITE(16,*) 'ENTERING GETstates'
      WRITE(16,*) 'effm =',effm
      IF (IPOT.EQ.0) then
         WRITE(6,*) 'ZEROING ELECTROSTATIC POTENTIAL'
         WRITE(16,*) 'ZEROING ELECTROSTATIC POTENTIAL'
      end if
      IF (IVAC.EQ.0) then
         WRITE(6,*) 'NO VACUUM REGION; FULL QUANTUM DOT'
         WRITE(16,*) 'NO VACUUM REGION; FULL QUANTUM DOT'
      end if
      END IF
c
c   construct periodic cartesian potential
c
      NX=2*NKX
      NY=2*NKY
      NZ=2*NKZ
      IF (IWRIT.GT.0) THEN
      write(6,*) 'energy cutoff =',emax
      write(16,*) 'energy cutoff =',emax
      write(6,*) 'number of k-points =',nkx,nky,nkz
      write(16,*) 'number of k-points =',nkx,nky,nkz
      END IF
      XMAX=NKX*PI/SQRT(C*EFFM*EMAX)
      YMAX=NKY*PI/SQRT(C*EFFM*EMAX)
      ZMAX=NKZ*PI/SQRT(C*EFFM*EMAX)
      EL1=ZMAX*2
      EL2=EL1+VACWID
      IF (IWRIT.GT.0) THEN
      WRITE(6,*) 'slab size X,Y,Z =',XMAX*2,YMAX*2,ZMAX*2
      WRITE(16,*) 'slab size X,Y,Z =',XMAX*2,YMAX*2,ZMAX*2
      END IF
      XDEL=XMAX/(NX-1)
      YDEL=YMAX/(NY-1)
      ZDEL=ZMAX/(NZ-1)
      IF (IWRIT.GT.0) THEN
      write(6,*) 'grid size X,Y,Z =',xdel,ydel,zdel
      write(16,*) 'grid size X,Y,Z =',xdel,ydel,zdel
      END IF
      CALL potperiod(SEP,NX,NY,NZ,XDEL,YDEL,ZDEL,
     &VAC,TIP,SEM,VSINT,PVAC,PSEM,R,S,DELV,DELR,DELS,DELP,NRDIM,
     &NVDIM,NSDIM,NPDIM,NXDIM,NXDIM2,NYDIM,NZDIM,NR,NV,NS,NP,IWRIT,IERR)
C
C   ADD THE OFFSET TO THE POTENTIAL
C
      IF (NX.GT.NXDIM.OR.NY.GT.NYDIM.OR.NZ.GT.NZDIM) THEN
         WRITE(6,*) '*** ERROR - INTERPOLATED POTENTIAL ARRAY TOO SMALL'
         WRITE(16,*)'*** ERROR - INTERPOLATED POTENTIAL ARRAY TOO SMALL'
         write(6,*) nx,ny,nz,nxdim,nydim,nzdim
         write(16,*) nx,ny,nz,nxdim,nydim,nzdim
         STOP
      END IF
      IF (NKX.GT.(NX-1).OR.NKY.GT.(NY-1).OR.NKZ.GT.(NZ-2)) THEN
         WRITE(6,*) '*** ERROR - NUMBER OF K-POINTS TOO LARGE'
         WRITE(16,*) '*** ERROR - NUMBER OF K-POINTS TOO LARGE'
         STOP
      END IF
      ELX=XMAX*2
      ELY=YMAX*2
      ELZ1=EL1
      ELZ2=EL2
      EL=EL2/2.
      DO 140 K=1,NZ
         Z=(K-1)*ZDEL
         DO 135 J=1,NY
            Y=(J-1)*YDEL
            DO 130 I=1,2*NX-2
               X=(I-NX+1)*XDEL
               IF (IPOT.EQ.0) PSEM(I,J,K)=0.
               IF (K.EQ.1) PSURF(I,J)=PSEM(I,J,K)
               IF (IBAND.EQ.-1) THEN
                  PSEM(I,J,K)=PSEM(I,J,K)+DELVBEDGE(X,Y,Z,I,J,K)
               ELSE
                  PSEM(I,J,K)=PSEM(I,J,K)+DELCBEDGE(X,Y,Z,I,J,K)
               END IF
               if (iwrit.ge.7) then
               if (k.eq.1) plot1(i,j)=psem(i,j,k)
               if (k.eq.2) plot2(i,j)=psem(i,j,k)
               if (k.eq.3) plot3(i,j)=psem(i,j,k)
               if (k.eq.4) plot4(i,j)=psem(i,j,k)
               if (k.eq.5) plot5(i,j)=psem(i,j,k)
               if (k.eq.6) plot6(i,j)=psem(i,j,k)
               if (k.eq.7) plot7(i,j)=psem(i,j,k)
               if (k.eq.8) plot8(i,j)=psem(i,j,k)
               if (k.eq.9) plot9(i,j)=psem(i,j,k)
               end if
130         CONTINUE
135      CONTINUE
140   CONTINUE
      if (iwrit.ge.7) then
         ilog=0
         iexpan=3
         shade=0.
         igraycut=1
        call PlotGray(plot1,nxdim2,nydim,nxdim2,nydim,ilog,iexpan,shade,
     &   ibin,ixbin,iybin,px,py,pz,ipow,ycorrec,nfilt,ixfilt,iyfilt,
     &   irfilt,disc,nscal,zsmin,zsmax,igraycut,grang,smin,smax,img)
         call dispir(1,img)
        call PlotGray(plot2,nxdim2,nydim,nxdim2,nydim,ilog,iexpan,shade,
     &   ibin,ixbin,iybin,px,py,pz,ipow,ycorrec,nfilt,ixfilt,iyfilt,
     &   irfilt,disc,nscal,zsmin,zsmax,igraycut,grang,smin,smax,img)
         call dispir(2,img)
        call PlotGray(plot3,nxdim2,nydim,nxdim2,nydim,ilog,iexpan,shade,
     &   ibin,ixbin,iybin,px,py,pz,ipow,ycorrec,nfilt,ixfilt,iyfilt,
     &   irfilt,disc,nscal,zsmin,zsmax,igraycut,grang,smin,smax,img)
         call dispir(3,img)
        call PlotGray(plot4,nxdim2,nydim,nxdim2,nydim,ilog,iexpan,shade,
     &   ibin,ixbin,iybin,px,py,pz,ipow,ycorrec,nfilt,ixfilt,iyfilt,
     &   irfilt,disc,nscal,zsmin,zsmax,igraycut,grang,smin,smax,img)
         call dispir(4,img)
        call PlotGray(plot5,nxdim2,nydim,nxdim2,nydim,ilog,iexpan,shade,
     &   ibin,ixbin,iybin,px,py,pz,ipow,ycorrec,nfilt,ixfilt,iyfilt,
     &   irfilt,disc,nscal,zsmin,zsmax,igraycut,grang,smin,smax,img)
         call dispir(5,img)
        call PlotGray(plot6,nxdim2,nydim,nxdim2,nydim,ilog,iexpan,shade,
     &   ibin,ixbin,iybin,px,py,pz,ipow,ycorrec,nfilt,ixfilt,iyfilt,
     &   irfilt,disc,nscal,zsmin,zsmax,igraycut,grang,smin,smax,img)
         call dispir(6,img)
        call PlotGray(plot7,nxdim2,nydim,nxdim2,nydim,ilog,iexpan,shade,
     &   ibin,ixbin,iybin,px,py,pz,ipow,ycorrec,nfilt,ixfilt,iyfilt,
     &   irfilt,disc,nscal,zsmin,zsmax,igraycut,grang,smin,smax,img)
         call dispir(7,img)
        call PlotGray(plot8,nxdim2,nydim,nxdim2,nydim,ilog,iexpan,shade,
     &   ibin,ixbin,iybin,px,py,pz,ipow,ycorrec,nfilt,ixfilt,iyfilt,
     &   irfilt,disc,nscal,zsmin,zsmax,igraycut,grang,smin,smax,img)
         call dispir(8,img)
        call PlotGray(plot9,nxdim2,nydim,nxdim2,nydim,ilog,iexpan,shade,
     &   ibin,ixbin,iybin,px,py,pz,ipow,ycorrec,nfilt,ixfilt,iyfilt,
     &   irfilt,disc,nscal,zsmin,zsmax,igraycut,grang,smin,smax,img)
         call dispir(9,img)
      end if
C
C   GET K VALUES IN X DIRECTION
C
C     CALL SLAB2(NKX,ELX,WKX)
      DO 160 KX=1,NKX
         WKX(KX)=(KX-1)*PI/XMAX
160   CONTINUE
      DO 165 KY=1,NKY
         WKY(KY)=(KY-1)*PI/YMAX
165   CONTINUE
C
C   LOOP OVER ODD AND EVEN STATES
C
      DO 990 IPARITY=0,1
      IF (IWRIT.GT.0) THEN
      IF (IPARITY.EQ.0) THEN
         WRITE(6,*) '*** ODD Z PARITY STATES'
         WRITE(16,*) '*** ODD Z PARITY STATES'
      ELSE
         WRITE(6,*) '*** EVEN Z PARITY STATES'
         WRITE(16,*) '*** EVEN Z PARITY STATES'
      END IF
      END IF
C
C   GET K VALUES IN Z-DIRECTION
C
      nkz2=nkz*2
      DO 175 KX=1,NKX
         DO 173 KY=1,NKY
            IF (IVAC.NE.0) THEN
            WKPARR2=WKX(KX)**2+WKY(KY)**2
            V0=BARR+WKPARR2*((-IBAND/EFFM)+1.)/C
C            WRITE(6,*) V0,ELZ1,ELZ2
c            READ(5,*)
            IF (IBAND.EQ.-1) THEN
               CALL BoundEner3V(deffm,ELZ1,ELZ2,v0,ener,amp1,amp2,nkz2)
            ELSE
c               write(6,*) 'V0,NKZ2 =',v0,nkz2
               CALL BoundEner3E(deffm,ELZ1,ELZ2,v0,ener,amp1,amp2,
     &                                                      nbound,nkz2)
            END IF
            do 170 kz=1,nkz
               wq(KX,KY,KZ)=dsqrt(C*ener(2*kz-IPARITY)*effm)
               ediff=v0-IBAND*ener(2*kz-IPARITY)
c
c   wkap>0 is indicator of sinh,cosh wavefucntion;
c   wkap<0 is indicator of sin,cos wavefunction
c
               if (ediff.gt.0) then
               wkap(KX,KY,KZ)=dsqrt(C*(v0-IBAND*ener(2*kz-IPARITY)))
               else
               wkap(KX,KY,KZ)=-dsqrt(-C*(v0-IBAND*ener(2*kz-IPARITY)))
               end if
C               write(6,*) v0-IBAND*ener(2*kz-IPARITY),wkap(KX,KY,KZ)
C               read(5,*)
               am2(KX,KY,KZ)=amp2(2*kz-IPARITY)
               am1(KX,KY,KZ)=amp1(2*kz-IPARITY)
170         continue
c      WRITE(6,*) '# KZ POINTS, MAX E =',NKZ,WQ(KX,KY,NKZ)**2/(C*EFFM)
         ELSE
         DO 172 KZ=1,NKZ
            IF (IPARITY.EQ.1) THEN
               WQ(KX,KY,KZ)=(KZ-1)*PI/ZMAX
            ELSE
               WQ(KX,KY,KZ)=(KZ-0.5)*PI/ZMAX
            END IF
            WKAP(KX,KY,KZ)=0.
            IF (IPARITY.EQ.1.AND.KZ.EQ.1) THEN
               am2(KX,KY,KZ)=1./SQRT(2.*ZMAX)
            ELSE
               am2(KX,KY,KZ)=1./SQRT(ZMAX)
            END IF
            am1(KX,KY,KZ)=0.
172      CONTINUE
         END IF      
173      CONTINUE
175   CONTINUE
      V0=BARR
C
C   COMPUTE WAVEFUNCTION VALUES
C
      DO 185 KX=1,NKX
         DO 180 I=1,2*NX-2
            IF (KX.EQ.1) THEN
               COSX(KX,I)=1./SQRT(2.*XMAX)
            ELSE
               COSX(KX,I)=COS((KX-1)*PI*(I-NX+1)/(NX-1))/SQRT(XMAX)
               SINX(KX-1,I)=SIN((KX-1)*PI*(I-NX+1)/(NX-1))/SQRT(XMAX)
            END IF
180      CONTINUE
185   CONTINUE
      DO 195 KY=1,NKY
         DO 190 J=1,NY
            IF (KY.EQ.1) THEN
               COSY(KY,J)=1./SQRT(2.*YMAX)
            ELSE
               COSY(KY,J)=COS((KY-1)*PI*(J-1)/(NY-1))/SQRT(YMAX)
            END IF
190      CONTINUE
195   CONTINUE
      ZVACDEL=0.05
      NZVAC=NINT((ELZ2-ELZ1)/(2.*ZVACDEL))
      IF (NZVAC.GT.NZVACDIM) THEN
         WRITE(6,*) '*** ERROR - NZVACDIM TOO SMALL, NZVAC =',NZVAC
         WRITE(16,*) '*** ERROR - NZVACDIM TOO SMALL, NZVAC =',NZVAC
         WRITE(6,*) 'PRESS ENTER TO EXIT PROGRAM'
         WRITE(16,*) 'PRESS ENTER TO EXIT PROGRAM'
         READ(5,*)
         STOP
      END IF
      ZVACDEL=(ELZ2-ELZ1)/(2.*(NZVAC-1))
      IF (IWRIT.NE.0) THEN
      WRITE(6,*) 'VACUUM NZ,DELZ =',NZVAC,ZVACDEL
      WRITE(16,*) 'VACUUM NZ,DELZ =',NZVAC,ZVACDEL
      END IF
      DO 250 KX=1,NKX
      DO 240 KY=1,NKY
      DO 210 KZ=1,NKZ
         DO 200 J=1,NZ
         IF (IPARITY.EQ.0) THEN
            HARM(KX,KY,KZ,J)=AM2(KX,KY,KZ)*SIN(WQ(KX,KY,KZ)*
     &                                                    (NZ-J)*ZDEL)
         ELSE
            HARM(KX,KY,KZ,J)=AM2(KX,KY,KZ)*COS(WQ(KX,KY,KZ)*
     &                                                    (NZ-J)*ZDEL)
         END IF
200      CONTINUE
         DO 205 J=1,NZVAC
            Z=(J-1)*ZVACDEL+(ELZ1/2.)
            if (wkap(kx,ky,kz).gt.0.) then
            IF (IPARITY.EQ.0) THEN
              TUNN(KX,KY,KZ,J)=AM1(KX,KY,KZ)*SINH(WKAP(KX,KY,KZ)*(EL-Z))
            ELSE
              TUNN(KX,KY,KZ,J)=AM1(KX,KY,KZ)*COSH(WKAP(KX,KY,KZ)*(EL-Z))
            END IF
            else
            IF (IPARITY.EQ.0) THEN
              TUNN(KX,KY,KZ,J)=AM1(KX,KY,KZ)*SIN(-WKAP(KX,KY,KZ)*(EL-Z))
            ELSE
              TUNN(KX,KY,KZ,J)=AM1(KX,KY,KZ)*COS(-WKAP(KX,KY,KZ)*(EL-Z))
            END IF
            end if            
205      CONTINUE
210   CONTINUE
240   CONTINUE
250   CONTINUE      
C
C   COMPUTE MATRIX ELEMENTS
C
      IF (IWRIT.GT.0) THEN
      write(6,*) 'COMPUTING MATRIX ELEMENTS'
      write(16,*) 'COMPUTING MATRIX ELEMENTS'
      END IF
      DO 500 KZP=1,NKZ
       DO 470 KYP=1,NKY
        DO 450 KXP=1,NKX
         KKP=KXP+(KYP-1)*NKX2M1+(KZP-1)*NKX2M1*NKY
         DO 420 KZ=1,NKZ
          DO 410 KY=1,NKY
           DO 400 KX=1,NKX
            KK=KX+(KY-1)*NKX2M1+(KZ-1)*NKX2M1*NKY
            if (kk.gt.kkp) go to 400
            SUM1=0.D0
            SUM2=0.D0
            SUM3=0.D0
            SUM4=0.D0
            DO 380 K=1,NZ
             Z=(J-1)*ZDEL
             DO 360 J=1,NY
              Y=(J-1)*YDEL
              DO 350 I=1,2*NX-2
               X=(I-NX+1)*XDEL
               TMP1=COSX(KXP,I)*COSY(KYP,J)*HARM(KXP,KYP,KZP,K)*
     &           PSEM(I,J,K)*COSX(KX,I)*COSY(KY,J)*HARM(KX,KY,KZ,K)
               IF (KXP.NE.1)
     &           TMP2=SINX(KXP-1,I)*COSY(KYP,J)*HARM(KXP,KYP,KZP,K)*
     &           PSEM(I,J,K)*COSX(KX,I)*COSY(KY,J)*HARM(KX,KY,KZ,K)
               IF (KX.NE.1)     
     &           TMP3=COSX(KXP,I)*COSY(KYP,J)*HARM(KXP,KYP,KZP,K)*
     &           PSEM(I,J,K)*SINX(KX-1,I)*COSY(KY,J)*HARM(KX,KY,KZ,K)
               IF (KXP.NE.1.AND.KX.NE.1)     
     &           TMP4=SINX(KXP-1,I)*COSY(KYP,J)*HARM(KXP,KYP,KZP,K)*
     &           PSEM(I,J,K)*SINX(KX-1,I)*COSY(KY,J)*HARM(KX,KY,KZ,K)
C               TMP1=COSX(KXP,I)*COSY(KYP,J)*HARM(KZP,K)*
C     &           COSX(KX,I)*COSY(KY,J)*HARM(KZ,K)
C               IF (KXP.NE.1)
C     &           TMP2=SINX(KXP-1,I)*COSY(KYP,J)*HARM(KZP,K)*
C     &           COSX(KX,I)*COSY(KY,J)*HARM(KZ,K)
C               IF (KX.NE.1)     
C     &           TMP3=COSX(KXP,I)*COSY(KYP,J)*HARM(KZP,K)*
C     &         SINX(KX-1,I)*COSY(KY,J)*HARM(KZ,K)
C               IF (KXP.NE.1.AND.KX.NE.1)     
C     &           TMP4=SINX(KXP-1,I)*COSY(KYP,J)*HARM(KZP,K)*
C     &         SINX(KX-1,I)*COSY(KY,J)*HARM(KZ,K)
               IF (J.NE.1.AND.J.NE.NY) THEN
                  TMP1=TMP1*2
                  TMP2=TMP2*2
                  TMP3=TMP3*2
                  TMP4=TMP4*2
               END IF
               IF (K.NE.1.AND.K.NE.NZ) THEN
                  TMP1=TMP1*2
                  TMP2=TMP2*2
                  TMP3=TMP3*2
                  TMP4=TMP4*2
               END IF
               SUM1=SUM1+TMP1
               IF (KXP.NE.1) SUM2=SUM2+TMP2
               IF (KX.NE.1) SUM3=SUM3+TMP3
               IF (KXP.NE.1.AND.KX.NE.1) SUM4=SUM4+TMP4
350           CONTINUE
360          CONTINUE
380         CONTINUE
            CMAT(KKP,KK)=SUM1*XDEL*YDEL*ZDEL
            IF (KXP.NE.1) CMAT(KKP+NKX-1,KK)=SUM2*XDEL*YDEL*ZDEL
            IF (KX.NE.1) CMAT(KKP,KK+NKX-1)=SUM3*XDEL*YDEL*ZDEL
            IF (KX.NE.1.AND.KXP.NE.1) 
     &            CMAT(KKP+NKX-1,KK+NKX-1)=SUM4*XDEL*YDEL*ZDEL
C      WRITE(6,*) KKP,KK,CMAT(KKP,KK)
C      IF (KXP.NE.1) WRITE(6,*) KKP+NKX-1,KK,CMAT(KKP+NKX-1,KK)
C      IF (KX.NE.1) WRITE(6,*) KKP,KK+NKX-1,CMAT(KKP,KK+NKX-1)
C      IF (KX.NE.1.AND.KXP.NE.1) WRITE(6,*) KKP+NKX-1,KK+NKX-1,
C     &   CMAT(KKP+NKX-1,KK+NKX-1)
C      READ(5,*)
            IF (KKP.EQ.KK) THEN
               CMAT(KKP,KK)=CMAT(KKP,KK)+IBAND*
     &            (WKX(KX)**2+WKY(KY)**2+WQ(KX,KY,KZ)**2)/(C*EFFM)
               IF (KX.NE.1.AND.KXP.NE.1)
     &            CMAT(KKP+NKX-1,KK+NKX-1)=CMAT(KKP+NKX-1,KK+NKX-1)+
     &            IBAND*(WKX(KX)**2+WKY(KY)**2+WQ(KX,KY,KZ)**2)/(C*EFFM)
            END IF
400        CONTINUE
410       CONTINUE
420      CONTINUE
450     CONTINUE
470    CONTINUE
500   CONTINUE
      IF (IVAC.NE.0) THEN
      DO 700 KZP=1,NKZ
       DO 670 KYP=1,NKY
        DO 650 KXP=1,NKX
         KKP=KXP+(KYP-1)*NKX2M1+(KZP-1)*NKX2M1*NKY
         DO 620 KZ=1,NKZ
          DO 610 KY=1,NKY
           DO 600 KX=1,NKX
            if (kk.gt.kkp) go to 600
            KK=KX+(KY-1)*NKX2M1+(KZ-1)*NKX2M1*NKY
            SUM1=0.D0
            SUM2=0.D0
            SUM3=0.D0
            SUM4=0.D0
C            JMAX1=MAX0(20,MIN0(100,NINT((10./WKAP(KX,KY,KZ))/ZVACDEL)))
C            JMAX2=MAX0(20,MIN0(100,NINT((10./WKAP(KXP,KYP,KZP))/ZVACDEL)))
C            JMAX=MIN0(JMAX1,JMAX2)
            DO 580 K=1,NZVAC
             Z=(K-1)*ZVACDEL
             DO 570 J=1,NY
              DO 550 I=1,2*NX-2
              TMP1=COSX(KXP,I)*COSY(KYP,J)*TUNN(KXP,KYP,KZP,K)*
     &           PSURF(I,J)*COSX(KX,I)*COSY(KY,J)*TUNN(KX,KY,KZ,K)
               IF (KXP.NE.1)
     &           TMP2=SINX(KXP-1,I)*COSY(KYP,J)*TUNN(KXP,KYP,KZP,K)*
     &           PSURF(I,J)*COSX(KX,I)*COSY(KY,J)*TUNN(KX,KY,KZ,K)
               IF (KX.NE.1)
     &           TMP3=COSX(KXP,I)*COSY(KYP,J)*TUNN(KXP,KYP,KZP,K)*
     &           PSURF(I,J)*SINX(KX-1,I)*COSY(KY,J)*TUNN(KX,KY,KZ,K)
               IF (KXP.NE.1.AND.KX.NE.1)
     &           TMP4=SINX(KXP-1,I)*COSY(KYP,J)*TUNN(KXP,KYP,KZP,K)*
     &           PSURF(I,J)*SINX(KX-1,I)*COSY(KY,J)*TUNN(KX,KY,KZ,K)
C              TMP1=COSX(KXP,I)*COSY(KYP,J)*TUNN(KZP,K)*
C     &           COSX(KX,I)*COSY(KY,J)*TUNN(KZ,K)
C               IF (KXP.NE.1)
C     &      TMP2=SINX(KXP-1,I)*COSY(KYP,J)*TUNN(KZP,K)*
C     &           COSX(KX,I)*COSY(KY,J)*TUNN(KZ,K)
C               IF (KX.NE.1)
C     &        TMP3=COSX(KXP,I)*COSY(KYP,J)*TUNN(KZP,K)*
C     &         SINX(KX-1,I)*COSY(KY,J)*TUNN(KZ,K)
C               IF (KXP.NE.1.AND.KX.NE.1)
C     &      TMP4=SINX(KXP-1,I)*COSY(KYP,J)*TUNN(KZP,K)*
C     &         SINX(KX-1,I)*COSY(KY,J)*TUNN(KZ,K)
               IF (J.NE.1.AND.J.NE.NY) THEN
                  TMP1=TMP1*2
                  TMP2=TMP2*2
                  TMP3=TMP3*2
                  TMP4=TMP4*2
               END IF
               IF (K.NE.1.AND.K.NE.NZ) THEN
                  TMP1=TMP1*2
                  TMP2=TMP2*2
                  TMP3=TMP3*2
                  TMP4=TMP4*2
               END IF
               SUM1=SUM1+TMP1
               IF (KXP.NE.1) SUM2=SUM2+TMP2
               IF (KX.NE.1) SUM3=SUM3+TMP3
               IF (KXP.NE.1.AND.KX.NE.1) SUM4=SUM4+TMP4
550           CONTINUE
570          CONTINUE
580         CONTINUE
            CMAT(KKP,KK)=CMAT(KKP,KK)+SUM1*XDEL*YDEL*ZVACDEL
            IF (KXP.NE.1)
     &      CMAT(KKP+NKX-1,KK)=CMAT(KKP+NKX-1,KK)+SUM2*XDEL*YDEL*ZVACDEL
            IF (KX.NE.1)
     &      CMAT(KKP,KK+NKX-1)=CMAT(KKP,KK+NKX-1)+SUM3*XDEL*YDEL*ZVACDEL
            IF (KX.NE.1.AND.KXP.NE.1)
     &           CMAT(KKP+NKX-1,KK+NKX-1)=CMAT(KKP+NKX-1,KK+NKX-1)+
     &           SUM4*XDEL*YDEL*ZVACDEL
C     WRITE(6,*) KKP,KK,CMAT(KKP,KK)
C     READ(5,*)
600        CONTINUE
610       CONTINUE
620      CONTINUE
650     CONTINUE
670    CONTINUE
700   CONTINUE
      END IF
      do 720 kkp=1,nk3
         do 710 kk=1,nk3
            if (kk.gt.kkp) cmat(kkp,kk)=cmat(kk,kkp)
710      continue
720   continue
C
C   SOLVE THE EIGENVALUE PROBLEM
C
      n    = nk3                                                        
C
      IF (IWRIT.GT.0) THEN
      WRITE(6,*) 'SOLVING EIGENVALUE PROBLEM'
      WRITE(16,*) 'SOLVING EIGENVALUE PROBLEM'
      END IF
C
      call  rs(n,n,cmat,w,1,zmat,aux1,aux2,ierr)
c
      IF (IBAND.EQ.-1) THEN
         IF (IWRIT.GT.0) THEN
         write(6,*) 'highest eigenvalues ='
         write(16,*) 'highest eigenvalues ='
         END IF
         do 780 i=NK3,MAX0(1,NK3-6),-1
            IF (IWRIT.GT.0) THEN
            write(6,*) w(i)
            write(16,*) w(i)
            END IF
780      continue
      ELSE
         IF (IWRIT.GT.0) THEN
         write(6,*) 'lowest eigenvalues ='
         write(16,*) 'lowest eigenvalues ='
         END IF
         do 785 i=MIN0(5,NK3),1,-1
            IF (IWRIT.GT.0) THEN
            write(6,*) w(i)
            write(16,*) w(i)
            END IF
785      continue
      END IF
      do 788 i=1,NK3
         IF (IWRIT.GE.4) write(35+IPARITY*10,*) w(i)
788   continue
C      READ(5,*)
C
C   EPSI2(1,*) HOLDS ENERGIES, AND EPSI2(2,*) HOLDS WAVEFCN^2;
C   EPSI2(3,*) HOLDS DECAY CONSTANT (ONLY VALID FOR DECAYING STATES)
C   EPSI2(4,*) HOLDS THE EXPECTATION VALUE FOR K-PARALLEL
C
      IF (IPARITY.EQ.0) THEN
      DO 795 I=1,NK3
         EPSI2(1,I)=W(I)
         sum1=0.d0
         sum2=0.d0
         sum3=0.d0
         sum4=0.d0
         do 790 k=1,nk3
            kz=((k-1)/(nkx2m1*nky))+1
            ky=((k-1-(kz-1)*nkx2m1*nky)/nkx2m1)+1
            kx=k-(kz-1)*nkx2m1*nky-(ky-1)*nkx2m1
            if (kx.le.nkx) then
               sum1=sum1+zmat(k,i)*COSX(KX,NX-1)*COSY(KY,1)*
     &            HARM(KX,KY,KZ,1)
               sum2=sum2+zmat(k,i)*COSX(KX,NX-1)*COSY(KY,1)* 
     &                                                 TUNN(KX,KY,KZ,1)
               sum3=sum3+zmat(k,i)*COSX(KX,NX-1)*COSY(KY,1)*
     &                                                 TUNN(KX,KY,KZ,2)
               sum4=sum4+zmat(k,i)**2*(WKX(KX)**2+WKY(KY)**2)
            else
               sum1=sum1+zmat(k,i)*SINX(KX-NKX,NX-1)*COSY(KY,1)*
     &            HARM(KX-NKX+1,KY,KZ,1)
               sum2=sum2+zmat(k,i)*SINX(KX-NKX,NX-1)*COSY(KY,1)* 
     &                                            TUNN(KX-NKX+1,KY,KZ,1)
               sum3=sum3+zmat(k,i)*SINX(KX-NKX,NX-1)*COSY(KY,1)*
     &                                            TUNN(KX-NKX+1,KY,KZ,2)
               sum4=sum4+zmat(k,i)**2*(WKX(KX-NKX+1)**2+WKY(KY)**2)
            end if
790      continue
         EPSI2(2,I)=sum1**2
         EPSI2(3,I)=(DLOG(ABS(SUM2))-DLOG(ABS(SUM3)))/ZVACDEL
         EPSI2(3,I)=amax1(0.,EPSI2(3,I))
         EPSI2(4,I)=sum4
795   CONTINUE
      ELSE
      DO 820 I=1,NK3
         DO 800 J=1,NK3+I-1
            IF (W(I).LT.EPSI2(1,J)) GO TO 810
800      CONTINUE
810      DO 812 K=NK3+I-1,J,-1
            EPSI2(1,K+1)=EPSI2(1,K)
            EPSI2(2,K+1)=EPSI2(2,K)
            EPSI2(3,K+1)=EPSI2(3,K)
            EPSI2(4,K+1)=EPSI2(4,K)
812      CONTINUE
         EPSI2(1,J)=W(I)
         sum1=0.d0
         sum2=0.d0
         sum3=0.d0
         sum4=0.d0
         do 815 k=1,NK3
            kz=((k-1)/(nkx2m1*nky))+1
            ky=((k-1-(kz-1)*nkx2m1*nky)/nkx2m1)+1
            kx=k-(kz-1)*nkx2m1*nky-(ky-1)*nkx2m1
            if (kx.le.nkx) then
               sum1=sum1+zmat(k,i)*COSX(KX,NX-1)*COSY(KY,1)*
     &            HARM(KX,KY,KZ,1)
               sum2=sum2+zmat(k,i)*COSX(KX,NX-1)*COSY(KY,1)*
     &                                                  TUNN(KX,KY,KZ,1)
               sum3=sum3+zmat(k,i)*COSX(KX,NX-1)*COSY(KY,1)*
     &                                                  TUNN(KX,KY,KZ,2)
               sum4=sum4+zmat(k,i)**2*(WKX(KX)**2+WKY(KY)**2)
            else
               sum1=sum1+zmat(k,i)*SINX(KX-NKX,NX-1)*COSY(KY,1)*
     &            HARM(KX-NKX+1,KY,KZ,1)
               sum2=sum2+zmat(k,i)*SINX(KX-NKX,NX-1)*COSY(KY,1)*
     &                                            TUNN(KX-NKX+1,KY,KZ,1)
               sum3=sum3+zmat(k,i)*SINX(KX-NKX,NX-1)*COSY(KY,1)*
     &                                            TUNN(KX-NKX+1,KY,KZ,2)
               sum4=sum4+zmat(k,i)**2*(WKX(KX-NKX+1)**2+WKY(KY)**2)
            end if
815      continue
         EPSI2(2,J)=sum1**2
         EPSI2(3,J)=(DLOG(ABS(SUM2))-DLOG(ABS(SUM3)))/ZVACDEL
         EPSI2(3,J)=amax1(0.,EPSI2(3,J))
         EPSI2(4,J)=sum4
820   CONTINUE
      do 825 I=2*NK3,1,-1
         IF (IWRIT.GE.4) write(40,*) EPSI2(1,I),EPSI2(2,I),EPSI2(3,I)
825   continue
      END IF
      NSTATES=2*NK3
C
C   OUTPUT WAVEFUNCTION VALUES
C
      IF (IBAND.EQ.-1) THEN
         IWF1=NK3
         IWF2=NK3-1
         IWF3=NK3-2
         IWF4=NK3-3
         IWF5=NK3-4
      ELSE
         IWF1=1
         IWF2=2
         IWF3=3
         IWF4=4
         IWF5=5
      END IF
      IF (IWRIT.GE.5) THEN
      DO 850 j=nz,1,-1
         Z=(j-1)*ZDEL
         write(31+IPARITY*10,*) (j-1)*ZDEL,PSEM(1,1,j)
         sum1=0.d0
         sum2=0.d0
         sum3=0.d0
         sum4=0.d0
         sum5=0.d0
         do 830 k=1,nk3
            kz=((k-1)/(nkx2m1*nky))+1
            ky=((k-1-(kz-1)*nkx2m1*nky)/nkx2m1)+1
            kx=k-(kz-1)*nkx2m1*nky-(ky-1)*nkx2m1
            if (kx.le.nkx) then
               sum1=sum1+zmat(k,IWF1)*COSX(KX,NX-1)*COSY(KY,1)*
     &            HARM(kx,ky,kz,j)
               sum2=sum2+zmat(k,IWF2)*COSX(KX,NX-1)*COSY(KY,1)*
     &            HARM(kx,ky,kz,j)
               sum3=sum3+zmat(k,IWF3)*COSX(KX,NX-1)*COSY(KY,1)*
     &            HARM(kx,ky,kz,j)
               sum4=sum4+zmat(k,IWF4)*COSX(KX,NX-1)*COSY(KY,1)*
     &            HARM(kx,ky,kz,j)
               sum5=sum5+zmat(k,IWF5)*COSX(KX,NX-1)*COSY(KY,1)*
     &            HARM(kx,ky,kz,j)
            else
               sum1=sum1+zmat(k,IWF1)*SINX(KX-NKX,NX-1)*
     &            COSY(KY,1)*HARM(kx-NKX+1,ky,kz,j)
               sum2=sum2+zmat(k,IWF2)*SINX(KX-NKX,NX-1)*
     &            COSY(KY,1)*HARM(kx-NKX+1,ky,kz,j)
               sum3=sum3+zmat(k,IWF3)*SINX(KX-NKX,NX-1)*
     &            COSY(KY,1)*HARM(kx-NKX+1,ky,kz,j)
               sum4=sum4+zmat(k,IWF4)*SINX(KX-NKX,NX-1)*
     &            COSY(KY,1)*HARM(kx-NKX+1,ky,kz,j)
               sum5=sum5+zmat(k,IWF5)*SINX(KX-NKX,NX-1)*
     &            COSY(KY,1)*HARM(kx-NKX+1,ky,kz,j)
            end if
830      continue
         write(32+IPARITY*10,840) z,sngl(sum1),sngl(sum2),sngl(sum3),
     &      sngl(sum4),sngl(sum5)
840      format(6g12.4)         
850   CONTINUE
      DO 900 j=1,20
         Z=(J-1)*ZVACDEL
         write(31+IPARITY*10,*) z,sngl(v0)
         sum1=0.d0
         sum2=0.d0
         sum3=0.d0
         sum4=0.d0
         sum5=0.d0
         do 880 k=1,nk3
            kz=((k-1)/(nkx2m1*nky))+1
            ky=((k-1-(kz-1)*nkx2m1*nky)/nkx2m1)+1
            kx=k-(kz-1)*nkx2m1*nky-(ky-1)*nkx2m1
            if (kx.le.nkx) then
               sum1=sum1+zmat(k,IWF1)*COSX(KX,NX-1)*COSY(KY,1)*
     &            TUNN(kx,ky,kz,j)
               sum2=sum2+zmat(k,IWF2)*COSX(KX,NX-1)*COSY(KY,1)*
     &            TUNN(kx,ky,kz,j)
               sum3=sum3+zmat(k,IWF3)*COSX(KX,NX-1)*COSY(KY,1)*
     &            TUNN(kx,ky,kz,j)
               sum4=sum4+zmat(k,IWF4)*COSX(KX,NX-1)*COSY(KY,1)*
     &            TUNN(kx,ky,kz,j)
               sum5=sum5+zmat(k,IWF5)*COSX(KX,NX-1)*COSY(KY,1)*
     &            TUNN(kx,ky,kz,j)
            else
               sum1=sum1+zmat(k,IWF1)*SINX(KX-NKX,NX-1)*
     &            COSY(KY,1)*TUNN(kx-NKX+1,ky,kz,j)
               sum2=sum2+zmat(k,IWF2)*SINX(KX-NKX,NX-1)*
     &            COSY(KY,1)*TUNN(kx-NKX+1,ky,kz,j)
               sum3=sum3+zmat(k,IWF3)*SINX(KX-NKX,NX-1)*
     &            COSY(KY,1)*TUNN(kx-NKX+1,ky,kz,j)
               sum4=sum4+zmat(k,IWF4)*SINX(KX-NKX,NX-1)*
     &            COSY(KY,1)*TUNN(kx-NKX+1,ky,kz,j)
               sum5=sum5+zmat(k,IWF5)*SINX(KX-NKX,NX-1)*
     &            COSY(KY,1)*TUNN(kx-NKX+1,ky,kz,j)
            end if
880      continue
         write(32+IPARITY*10,890) -Z,sngl(sum1),sngl(sum2),sngl(sum3),
     &      sngl(sum4),sngl(sum5)
890      format(6g12.4)         
900   CONTINUE
      DO 950 i=1,2*NX-2
         write(33+IPARITY*10,*) (i-NX+1)*XDEL,PSEM(i,1,1)
         sum1=0.d0
         sum2=0.d0
         sum3=0.d0
         sum4=0.d0
         sum5=0.d0
         do 930 k=1,NK3
            kz=((k-1)/(nkx2m1*nky))+1
            ky=((k-1-(kz-1)*nkx2m1*nky)/nkx2m1)+1
            kx=k-(kz-1)*nkx2m1*nky-(ky-1)*nkx2m1
            if (kx.le.nkx) then
               sum1=sum1+zmat(k,IWF1)*COSX(KX,I)*COSY(KY,1)*
     &            HARM(kx,ky,kz,1)
               sum2=sum2+zmat(k,IWF2)*COSX(KX,I)*COSY(KY,1)*
     &            HARM(kx,ky,kz,1)
               sum3=sum3+zmat(k,IWF3)*COSX(KX,I)*COSY(KY,1)*
     &            HARM(kx,ky,kz,1)
               sum4=sum4+zmat(k,IWF4)*COSX(KX,I)*COSY(KY,1)*
     &            HARM(kx,ky,kz,1)
               sum5=sum5+zmat(k,IWF5)*COSX(KX,I)*COSY(KY,1)*
     &            HARM(kx,ky,kz,1)
            else
               sum1=sum1+zmat(k,IWF1)*SINX(KX-NKX,I)*COSY(KY,1)*
     &            HARM(kx-NKX+1,ky,kz,1)
               sum2=sum2+zmat(k,IWF2)*SINX(KX-NKX,I)*COSY(KY,1)*
     &            HARM(kx-NKX+1,ky,kz,1)
               sum3=sum3+zmat(k,IWF3)*SINX(KX-NKX,I)*COSY(KY,1)*
     &            HARM(kx-NKX+1,ky,kz,1)
               sum4=sum4+zmat(k,IWF4)*SINX(KX-NKX,I)*COSY(KY,1)*
     &            HARM(kx-NKX+1,ky,kz,1)
               sum5=sum5+zmat(k,IWF5)*SINX(KX-NKX,I)*COSY(KY,1)*
     &            HARM(kx-NKX+1,ky,kz,1)
            end if
930      continue
         write(34+IPARITY*10,940) (I-NX+1)*XDEL,sngl(sum1),sngl(sum2),
     &      sngl(sum3),sngl(sum4),sngl(sum5)
940      format(6g12.4)
950   CONTINUE
      END IF
C
C   OUTPUT LDOS
C
      IF (IWRIT.GE.1) THEN
      DO 970 IWF1=1,NK3
         sum1=0.d0
         do 960 k=1,nk3
            kz=((k-1)/(nkx2m1*nky))+1
            ky=((k-1-(kz-1)*nkx2m1*nky)/nkx2m1)+1
            kx=k-(kz-1)*nkx2m1*nky-(ky-1)*nkx2m1
            if (kx.le.nkx) then
               sum1=sum1+zmat(k,IWF1)*COSX(KX,NX-1)*COSY(KY,1)*
     &            HARM(kx,ky,kz,1)
            else
               sum1=sum1+zmat(k,IWF1)*SINX(KX-NKX,NX-1)*
     &            COSY(KY,1)*HARM(kx-NKX+1,ky,kz,1)
            end if
960      continue
         IF (IWRIT.GE.6) write(80+IPARITY,*) sngl(w(iwf1)),
     &                                                  abs(sngl(sum1))
970   CONTINUE
      END IF
990   CONTINUE
c
      RETURN
      END
