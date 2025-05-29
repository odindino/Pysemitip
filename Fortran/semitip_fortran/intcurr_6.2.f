C   ***************** INTCURR ******************
C
C   ROUTINES FOR COMPUTING TUNNEL CURRENT FOR EFFECTIVE MASS BULK 
C   BANDS, USING BARDEEN FORMALISM AND T&H APPROXIMATION, BY SOLVING
C   SCHRODINGER EQN BY 1D NUMERICAL INTEGRATION OF POTENTIAL.
C
C   IF INPUT WITH NWK=1, THEN ONLY COMPUTE QUANTUM CHARGE DENSITIES
C
C   VERSION 1.0 - WRITTEN BY R. M. FEENSTRA, 2007-2009
C           2.0 - AUG/09, INCLUDES ARRAYS CDESEM,CDESURF,CDLSEM,CDLSURF,
C                 CDEVAC,CDLVAC FOR SELF-CONSISTENT CHARGE DENSITIES;
C                 HENCE COMPATIBLE WITH SEMISTIP2 VERSION 4 
C           2.1 - DEC/10, ELIMINATE USAGE OF 2-BAND MODEL, AND MODIFY
C                         OUTPUT OF CURRENT VALUES
C           6.0 - FEB/11, REPACKAGE ROUTINES TO BE COMPATIBLE WITH SEMITIP V6
C           6.1 - FEB/11, ELIMINATE USE OF /SEMI/ COMMON BLOCK
C           6.2 - MAY/13, ADD IMAGE POTENTIAL OPTION
C
C   CALLS POTEXPAND VERSION 6.1 FOR EXPANDING POTENTIAL CURVE INTO
C   A FORM SUITABLE FOR NUMERICAL INTEGRATION.
C
C   EV,EVSO,EC = BAND EDGE POSITIONS DEEP INSIDE SEMICONDUCTOR
C   PROF=ELECTROSTATIC POTENTIAL IN SEMICONDUCTOR
C   VBPROF=VB EDGE PROFILE, INCLUDING ELECTROSTATIC POTENTIAL
C   CBPROF=CB EDGE PROFILE, INCLUDING ELECTROSTATIC POTENTIAL
C   PMAX=MAXIMUM ENERGY IN VB PROFILE
C   PMIN=MINIMUM ENERGY IN CB PROFILE
C
C   OUTPUT INTO LOGICAL UNIT NUMBER LUN+I:
C      I=0, QUANTUM NO. AND ENERGY OF LOCALIZED STATES, VS BIAS
C      I=1, WAVEFUNCTIONS OF LOCALIZED STATES  (AT FINAL BIAS)
C      I=2, CHARGE DENSITIES OF LOCALIZED STATES  (AT FINAL BIAS)
C      I=3, EXTENDED STATE WAVEFCN AT KPARR=0 (AT FINAL BIAS)
C
      SUBROUTINE INTCURR(IMPOT,BARR,PROF,NBARR1,NV,NS,NSP,
     &   NVDIM,NSDIM,S,SEP,BIAS,EF,CHI,EFTIP,CPot,EGAP,TK,AVBH,AVBL,
     &   AVBSO,ACB,ESO,E2HH,E2SO,nee,nwk,Pot0,NVDIM1,NVDIM2,NSDIM2,
     &   EXPANI,NLOC,CURRV,CURRV0,CURRC,CURRC0,CURR,CURR0,IWRIT,
     &   ICOMP,CDESEM,CDESURF,CDLSEM,CDLSURF,CDEVAC,CDLVAC)
C
      DIMENSION S(NSDIM),BARR(NVDIM1),PROF(NSDIM),NLOC(4),
     &CDESEM(NSDIM),CDLSEM(NSDIM),CDEVAC(NVDIM1),CDLVAC(NVDIM1),
     &VBPROF(NSDIM),CBPROF(NSDIM)
      real kappa,lambda
      DOUBLE PRECISION SUM1,SUM2,SUM1S,SUM2S,SUM1P,SUM2P
      DATA RQUANT/12900./
      PI=4.*ATAN(1.)
      tk1=tk
      tk2=tk
      IF (ICOMP.EQ.1) THEN
      CDESURF=0.
      CDLSURF=0.
      DO 65 J=1,NSP
         CDESEM(J)=0.
         CDLSEM(J)=0.
65    CONTINUE            
      DO 75 J=1,NV+1
         CDEVAC(J)=0.
         CDLVAC(J)=0.
75    CONTINUE
      END IF
c
c   valence band
c
      DO 100 J=1,NSP
         SZ=S(J)
         VBPROF(J)=PROF(J)+VBEDGE(SZ)
         IF (J.EQ.1) THEN
            PMAX=VBPROF(J)
         ELSE
            PMAX=AMAX1(PMAX,VBPROF(J))
         END IF
100   CONTINUE
      EV=PROF(NS)+VBEDGE(S(NS))
      IF (iwrit.eq.4) THEN
         write(6,*) 'MAXIMUM ENERGY IN VB PROFILE =',PMAX
         write(16,*) 'MAXIMUM ENERGY IN VB PROFILE =',PMAX
         WRITE(6,*) 'VB EDGE ENERGY DEEP INSIDE SEMICONDUCTOR =',EV
         WRITE(16,*) 'VB EDGE ENERGY DEEP INSIDE SEMICONDUCTOR =',EV
      END IF
      IF (IWRIT.GE.1) THEN
         DO 110 J=1,NSP
            WRITE(96,*) S(J),VBPROF(J)
110      CONTINUE
         CLOSE(96)
      END IF
      LUN=30
      CALL VBCURR1(IMPOT,BARR,NBARR1,S,PROF,NSP,NSDIM,NVDIM,NV,NVDIM1,
     &NSDIM2,NVDIM2,SEP,Pot0,EXPANI,PMAX,AVBL,TK1,TK2,BIAS,
     &EF,NEE,NWK,NL,EV,EFTIP,EGAP,LUN,CDESEM,CDESURF,CDLSEM,
     &CDLSURF,CDEVAC,CDLVAC,CURRVL,CURRV0L,IWRIT,ICOMP,VBPROF)
      IF (IWRIT.NE.0) THEN
      write(6,*) 'number of VB light-hole localized states =',nl
      write(16,*) 'number of VB light-hole localized states =',nl
      END IF
      NLOC(1)=NL
      IF (IWRIT.GE.2) CLOSE(LUN+1)
      IF (IWRIT.GE.2) CLOSE(LUN+2)
      IF (IWRIT.GE.4) CLOSE(LUN+3)
      LUN=40
      CALL VBCURR1(IMPOT,BARR,NBARR1,S,PROF,NSP,NSDIM,NVDIM,NV,NVDIM1,
     &NSDIM2,NVDIM2,SEP,Pot0,EXPANI,PMAX,AVBH,TK1,TK2,BIAS,
     &EF,NEE,NWK,NL,EV,EFTIP,E2HH,LUN,CDESEM,CDESURF,CDLSEM,
     &CDLSURF,CDEVAC,CDLVAC,CURRVH,CURRV0H,IWRIT,ICOMP,VBPROF)
      IF (IWRIT.NE.0) THEN
      write(6,*) 'number of VB heavy-hole localized states =',nl
      write(16,*) 'number of VB heavy-hole localized states =',nl
      END IF
      NLOC(2)=NL
      IF (IWRIT.GE.2) CLOSE(LUN+1)
      IF (IWRIT.GE.2) CLOSE(LUN+2)
      IF (IWRIT.GE.4) CLOSE(LUN+3)
      EVSO=EV-ESO
      PMAXSO=PMAX-ESO
      DO 200 J=1,NSP
         VBPROF(J)=VBPROF(J)-ESO
200   CONTINUE         
      LUN=50
      CALL VBCURR1(IMPOT,BARR,NBARR1,S,PROF,NSP,NSDIM,NVDIM,NV,NVDIM1,
     &NSDIM2,NVDIM2,SEP,Pot0,EXPANI,PMAXSO,AVBSO,TK1,TK2,BIAS,
     &EF,NEE,NWK,NL,EVSO,EFTIP,E2SO,LUN,CDESEM,CDESURF,CDLSEM,
     &CDLSURF,CDEVAC,CDLVAC,CURRVSO,CURRV0SO,IWRIT,ICOMP,VBPROF)
      IF (IWRIT.NE.0) THEN
      write(6,*) 'number of VB split-off localized states =',nl
      write(16,*) 'number of VB split-off localized states =',nl
      END IF
      NLOC(3)=NL
      IF (IWRIT.GE.2) CLOSE(LUN+1)
      IF (IWRIT.GE.2) CLOSE(LUN+2)
      IF (IWRIT.GE.4) CLOSE(LUN+3)
      currv=currvL+currvH+currvSO
      currv0=currv0L+currv0H+currv0SO
c
c   conduction band
c
      DO 300 J=1,NSP
         SZ=S(J)
         CBPROF(J)=PROF(J)+CBEDGE(SZ)
         IF (J.EQ.1) THEN
            PMIN=CBPROF(J)
         ELSE
            PMIN=AMIN1(PMIN,CBPROF(J))
         END IF
300   CONTINUE         
      EC=PROF(NS)+CBEDGE(S(NS))
      IF (iwrit.eq.4) THEN
         write(6,*) 'MINIMUM ENERGY IN CB PROFILE =',PMIN
         write(16,*) 'MINIMUM ENERGY IN CB PROFILE =',PMIN
         WRITE(6,*) 'CB EDGE ENERGY DEEP INSIDE SEMICONDUCTOR =',EC
         WRITE(16,*) 'CB EDGE ENERGY DEEP INSIDE SEMICONDUCTOR =',EC
      END IF
      IF (IWRIT.GE.1) THEN
         DO 310 J=1,NSP
            WRITE(97,*) S(J),CBPROF(J)
310      CONTINUE
         CLOSE(97)
      END IF
      LUN=60
      CALL CBCURR1(IMPOT,BARR,NBARR1,S,PROF,NSP,NSDIM,NVDIM,NV,NVDIM1,
     &NSDIM2,NVDIM2,SEP,Pot0,EXPANI,PMIN,ACB,TK1,TK2,BIAS,
     &EF,NEE,NWK,NL,EC,EFTIP,EGAP,LUN,CDESEM,CDESURF,CDLSEM,
     &CDLSURF,CDEVAC,CDLVAC,CURRC,CURRC0,IWRIT,ICOMP,CBPROF)
      IF (IWRIT.NE.0) THEN
      write(6,*) 'number of CB localized states =',nl
      write(16,*) 'number of CB localized states =',nl
      END IF
      NLOC(4)=NL
      IF (IWRIT.GE.2) CLOSE(LUN+1)
      IF (IWRIT.GE.2) CLOSE(LUN+2)
      IF (IWRIT.GE.4) CLOSE(LUN+3)
      curr=currv+currc
      curr0=currv0+currc0
      return
      end
C
C   VB TUNNEL CURRENT
C
      SUBROUTINE VBCURR1(IMPOT,BARR,NBARR1,S,PROF,NSP,NSDIM,NVDIM,NV,
     &NVDIM1,NSDIM2,NVDIM2,SEP,Pot0,EXPANI,PMAX,EFFM,TK1,TK2,BIAS,
     &EF,NE,NWK,NLOC,EV,EFTIP,E2BAND,LUN,CDESEM,CDESURF,CDLSEM,
     &CDLSURF,CDEVAC,CDLVAC,CURRV,CURRV0,IWRIT,ICOMP,VBPROF)
C
      REAL KAPPA
      DIMENSION BARR(NVDIM1),PROF(NSDIM),S(NSDIM),BARR2(NVDIM2),
     &PROF2(NSDIM2),S2(NSDIM2),JSEM(NSDIM2),NEXSEM(NSDIM),
     &CDESEM(NSDIM),CDLSEM(NSDIM),CDEVAC(NVDIM1),CDLVAC(NVDIM1),
     &VBPROF(NSDIM)
      DOUBLE PRECISION SUM,SUM2,SUM3,PSIVAC(NVDIM2),PSISEM(NSDIM2)
c   value of C below is 2m/hbar^2 in units of 1/(eV nm^2)
      DATA C/26.254/RQUANT/12900./
      PI=4.*ATAN(1.)
      NLOC=0
      CURRV=0.
      CURRV0=0.
      WKFTIP=sqrt(C*EFTIP)
C
C   EXTENDED STATES
C
      emax=EV
      IF (NWK.NE.1) THEN      
         emin=amin1(ef-10.*tk1,ef+bias-10.*tk2)
      ELSE
         emin=ef-10.*tk1
      END IF
      dele=(emax-emin)/ne
      sum=0.
      if (dele.le.0.) go to 200
      wkmax=sqrt(C*EFFM*(emax-emin))
      SEMSTEP=(2.*pi/wkmax)/EXPANI
      KAPPA=SQRT(C*(amax1(BARR(NBARR1),BARR(1))-EMIN))
      VACSTEP=(2.*pi/KAPPA)/EXPANI
      Pot0p=Pot0+VBEDGE(0.)
      CALL POTEXPAND(IMPOT,SEP,NV,Pot0p,S,NSP,NSDIM,BARR,NBARR1,BARR2,
     &NBARR2,NVDIM1,NVDIM2,VBPROF,PROF2,NSDIM2,S2,NS2,VACSTEP,SEMSTEP,
     &JSEM,NEXSEM,NEXVAC,IWRIT)
      DELVAC=sep/float(NBARR2-1)
      delwk=wkmax/nwk
      do 120 iwky=0,nwk-1
         wky=iwky*delwk
         do 115 iwkx=0,nwk-1
         wkx=iwkx*delwk
         wkparr=sqrt(wkx**2+wky**2)
         eparr=wkparr**2/(EFFM*C)
         nwkdeg=8
         if (iwkx.eq.0) nwkdeg=nwkdeg/2
         if (iwky.eq.0) nwkdeg=nwkdeg/2
         if (iwkx.eq.iwky) nwkdeg=nwkdeg/2
         if (iwky.gt.iwkx) go to 115
         do 110 ie=1,ne
            ener=emax-(ie-0.5)*dele
            if (eparr.ge.(emax-ener)) go to 110
            occtip=fd(ener-bias,ef,tk2)
            occsem=fd(ener,ef,tk1)
            occdiff=occtip-occsem
            call VBwf(IMPOT,wf,wfderiv,WKSEM,ener,wkparr,sep,bias,
     &        BARR2,NVDIM2,NBARR2,PROF2,NS2,NSDIM2,S2,EFFM,EV,E2BAND,
     &        PSIVAC,PSISEM)
            IF (IWKX.EQ.0.AND.IWKY.EQ.0) THEN

C   WRITE OUTPUT IF DESIRED
C
               IF (IWRIT.GE.4) THEN
                 do 50 J=NBARR2,1,-1
                    write(LUN+3,*) -(J-1)*DELVAC,PSIVAC(J)
50               CONTINUE
                 do 60 J=1,NS2
                    write(LUN+3,*) S2(J),PSISEM(J)
60               CONTINUE
               END IF
C
C   INCLUDE IN CHARGE DENSITY
C
               IF (ICOMP.EQ.1) THEN
               IF (IE.EQ.1) THEN
                  SUM3=0.
                  do 70 ieocc=ie,ne
                     enerocc=emax-ieocc*dele
                     OCCSEM2=1.-FD(ENEROCC,EF,TK1)
                     SUM3=SUM3+OCCSEM2
70                continue
               ELSE
                  enerocc=emax-(ie-1)*dele
                  OCCSEM2=1.-FD(ENEROCC,EF,TK1)
                  SUM3=SUM3-OCCSEM2
               END IF
               OCCINT=SUM3*DELE
               CMULT=OCCINT*(EFFM*C/(2.*PI))*(DELE*EFFM*C/(2.*PI*WKSEM))
               SUM2=0.
               do 80 j=NBARR2,1,-1
                  TMP=(PSIVAC(J))**2*CMULT
                  IF (MOD(J-1,NEXVAC).EQ.0) THEN
                     JTMP=((J-1)/NEXVAC)+1
                     CDEVAC(JTMP)=CDEVAC(JTMP)+TMP
                  END IF
                  IF (J.NE.1) THEN
                     SUM2=SUM2+TMP*DELVAC
                  ELSE
                     SUM2=SUM2+TMP*0.5*DELVAC
                  END IF
80             CONTINUE
               CDESURF=CDESURF+SUM2
               do 90 j=1,NS2
                  TMP=(PSISEM(J))**2*CMULT
                  IF (JSEM(J).NE.0) THEN
                     CDESEM(JSEM(J))=CDESEM(JSEM(J))+
     &                                          TMP/NEXSEM(JSEM(J))
                  ELSE
                     WRITE(6,*) 'ERROR , INTCURR - J=0'
                     READ(5,*)
                  END IF
90             CONTINUE
            END IF
            END IF
C
C   EVALUATE EVERYTHING; eperp, kappa in vacuum
C
            EPERP=ener-WKPARR**2/C
            KAPPA=SQRT(C*(BARR2(NBARR2)-EPERP))
            trans=2.*nwkdeg*(2.*wf)**2*WKFTIP/(WKSEM/EFFM)
            sum=sum+trans*occdiff
110      continue
115      continue
120   continue
      currv=sum*dele*delwk**2/(4.*PI**2*RQUANT)
C
C   LOCALIZED STATES
C
200   IF (PMAX.LE.EV) RETURN
      emax=PMAX
c   could delete statements below!
      IF (NWK.NE.1) THEN      
         emin=amax1(EV,amin1(ef-10.*tk1,ef+bias-10.*tk2))
      ELSE
         emin=amax1(EV,ef-10.*tk1)
      END IF
      dele=(emax-emin)/ne
      emin2=ef-10.*tk1
      dele2=(emax-emin2)/ne
      sum=0.
      if (dele.le.0.) go to 460
      wkmax=sqrt(C*EFFM*(emax-emin))
      SEMSTEP=(2.*pi/wkmax)/EXPANI
      KAPPA=SQRT(C*(amax1(BARR(NBARR1),BARR(1))-EMIN))
      VACSTEP=(2.*pi/KAPPA)/EXPANI
      Pot0p=Pot0+VBEDGE(0.)
      CALL POTEXPAND(IMPOT,SEP,NV,Pot0p,S,NSP,NSDIM,BARR,NBARR1,BARR2,
     &NBARR2,NVDIM1,NVDIM2,VBPROF,PROF2,NSDIM2,S2,NS2,VACSTEP,SEMSTEP,
     &JSEM,NEXSEM,NEXVAC,IWRIT)
      DELVAC=sep/float(NBARR2-1)
      delwk=wkmax/nwk
      do 320 iwky=0,nwk-1
         wky=iwky*delwk
         do 315 iwkx=0,nwk-1
         wkx=iwkx*delwk
         wkparr=sqrt(wkx**2+wky**2)
         eparr=wkparr**2/(EFFM*C)
         n=0
         nsav=0
         nwkdeg=8
         if (iwkx.eq.0) nwkdeg=nwkdeg/2
         if (iwky.eq.0) nwkdeg=nwkdeg/2
         if (iwkx.eq.iwky) nwkdeg=nwkdeg/2
         if (iwky.gt.iwkx) go to 315
         do 310 ie=1,ne
            ener=emax-(ie-0.5)*dele
            if (eparr.ge.(emax-ener)) go to 310
            occtip=fd(ener-bias,ef,tk2)
            occsem=fd(ener,ef,tk1)
            occdiff=occtip-occsem
            call VBloc(IMPOT,n,wf,wfderiv,ener,wkparr,sep,bias,BARR2,
     &          NVDIM2,NBARR2,PROF2,NS2,NSDIM2,S2,EFFM,EV,E2BAND,PSIVAC,
     &          PSISEM)
            if (n.eq.nsav) go to 310
            IF (IWKX.EQ.0.AND.IWKY.EQ.0) THEN
               IF (PSISEM(1).NE.0) THEN
                  NLOC=NLOC+1
C
C   FOUND LOCALIZED STATE; WRITE OUTPUT IF DESIRED
C
                  IF (IWRIT.GT.1) THEN
                     write(6,*) 'VB localized state at energy ',ener
                     write(16,*)'VB localized state at energy ',ener
                  END IF
                  IF (IWRIT.GE.1) THEN
                     WRITE(LUN,*) BIAS,N,ENER,OCCDIFF
                  END IF
                  IF (IWRIT.GE.2) THEN
                     do 250 J=NBARR2,1,-1
                        write(LUN+1,*) -(J-1)*DELVAC,PSIVAC(J)
250                  CONTINUE
                     do 260 J=1,NS2
                        write(LUN+1,*) S2(J),PSISEM(J)
260                  CONTINUE
                  END IF
C
C   INCLUDE IN CHARGE DENSITY
C
                  IF (ICOMP.EQ.1) THEN
                  OCCINT=AMAX1(0.,(ENER-EF))
                  IF (TK1.NE.0.) THEN
                     SUM2=0.D0
                     DO 265 IEOCC=1,NE
                        ENEROCC=ENER-(IEOCC-0.5)*DELE2
                        IF (ENEROCC.LT.EMIN2) GO TO 270
                        OCCSEM2=1.-FD(ENEROCC,EF,TK1)
                        SUM2=SUM2+OCCSEM2
265                  CONTINUE
270                  CONTINUE      
                     SUM2=SUM2*DELE2
                     OCCINT=SUM2
                  END IF
                  CMULT=OCCINT*(EFFM*C/(2.*PI))
                  SUM2=0.
                  do 275 j=NBARR2,1,-1
                     TMP=(PSIVAC(J))**2*CMULT
                     IF (MOD(J-1,NEXVAC).EQ.0) THEN
                        JTMP=((J-1)/NEXVAC)+1
                        CDLVAC(JTMP)=CDLVAC(JTMP)+TMP
                     END IF
                     IF (IWRIT.GE.2) WRITE(LUN+2,*) -(J-1)*DELVAC,TMP
                     IF (J.NE.1) THEN
                        SUM2=SUM2+TMP*DELVAC
                     ELSE
                        SUM2=SUM2+TMP*0.5*DELVAC
                     END IF
275                 CONTINUE
                  CDLSURF=CDLSURF+SUM2
                  do 290 j=1,NS2
                     TMP=(PSISEM(J))**2*CMULT
                     IF (IWRIT.GE.2) WRITE(LUN+2,*) S2(J),TMP
                     IF (JSEM(J).NE.0) THEN
                       CDLSEM(JSEM(J))=CDLSEM(JSEM(J))+
     &                                            TMP/NEXSEM(JSEM(J))
                     ELSE
                     WRITE(6,*) 'ERROR , INTCURR - J=0'
                     READ(5,*)
                     END IF
290               CONTINUE
                  END IF
               END IF
            END IF
C
C   EVALUATE EVERYTHING; eperp, kappa in vacuum
C
            EPERP=ener-WKPARR**2/C
            KAPPA=SQRT(C*(BARR2(NBARR2)-EPERP))
            trans=nwkdeg*(n-nsav)*2.*(2.*wf)**2*WKFTIP
            sum=sum+trans*occdiff
            nsav=n
310      continue
315      continue
320   continue
      currv0=sum*delwk**2/(C*2.*PI*RQUANT)
460   RETURN
      END
C
C   CB TUNNEL CURRENT
C
      SUBROUTINE CBCURR1(IMPOT,BARR,NBARR1,S,PROF,NSP,NSDIM,NVDIM,NV,
     &NVDIM1,NSDIM2,NVDIM2,SEP,Pot0,EXPANI,PMIN,EFFM,TK1,TK2,BIAS,
     &EF,NE,NWK,NLOC,EC,EFTIP,E2BAND,LUN,CDESEM,CDESURF,CDLSEM,
     &CDLSURF,CDEVAC,CDLVAC,CURRC,CURRC0,IWRIT,ICOMP,CBPROF)
C
      REAL KAPPA
      DIMENSION BARR(NVDIM1),PROF(NSDIM),S(NSDIM),BARR2(NVDIM2),
     &PROF2(NSDIM2),S2(NSDIM2),JSEM(NSDIM2),NEXSEM(NSDIM),
     &CDESEM(NSDIM),CDLSEM(NSDIM),CDEVAC(NVDIM1),CDLVAC(NVDIM1),
     &CBPROF(NSDIM)
      DOUBLE PRECISION SUM,SUM2,SUM3,PSIVAC(NVDIM2),PSISEM(NSDIM2)
c   value of C below is 2m/hbar^2 in units of 1/(eV nm^2)
      DATA C/26.254/RQUANT/12900./
      PI=4.*ATAN(1.)
      NLOC=0
      CURRC=0.
      CURRC0=0.
      WKFTIP=sqrt(C*EFTIP)
C
C   EXTENDED STATES
C
      emin=EC
      IF (NWK.NE.1) THEN
         emax=amax1(ef+10.*tk1,ef+bias+10.*tk2)
      ELSE
         emax=ef+10.*tk1
      END IF
      dele=(emax-emin)/ne
      sum=0.
      if (dele.le.0.) go to 200
      wkmax=sqrt(C*EFFM*(emax-emin))
      SEMSTEP=(2.*pi/wkmax)/EXPANI
      KAPPA=SQRT(C*(amax1(BARR(NBARR1),BARR(1))-EMIN))
      VACSTEP=(2.*pi/KAPPA)/EXPANI
      Pot0p=Pot0+CBEDGE(0.)
      CALL POTEXPAND(IMPOT,SEP,NV,Pot0p,S,NSP,NSDIM,BARR,NBARR1,BARR2,
     &NBARR2,NVDIM1,NVDIM2,CBPROF,PROF2,NSDIM2,S2,NS2,VACSTEP,SEMSTEP,
     &JSEM,NEXSEM,NEXVAC,IWRIT)
      DELVAC=sep/float(NBARR2-1)
      delwk=wkmax/nwk
      do 120 iwky=0,nwk-1
         wky=iwky*delwk
         do 115 iwkx=0,nwk-1
         wkx=iwkx*delwk
         wkparr=sqrt(wkx**2+wky**2)
         eparr=wkparr**2/(EFFM*C)
         nwkdeg=8
         if (iwkx.eq.0) nwkdeg=nwkdeg/2
         if (iwky.eq.0) nwkdeg=nwkdeg/2
         if (iwkx.eq.iwky) nwkdeg=nwkdeg/2
         if (iwky.gt.iwkx) go to 115
         do 110 ie=1,ne
            ener=emin+(ie-0.5)*dele
            if (eparr.ge.(ener-emin)) go to 110
            occtip=fd(ener-bias,ef,tk2)
            occsem=fd(ener,ef,tk1)
            occdiff=occtip-occsem
            call CBwf(IMPOT,wf,wfderiv,WKSEM,ener,wkparr,sep,bias,
     &        BARR2,NVDIM2,NBARR2,PROF2,NS2,NSDIM2,S2,EFFM,EC,E2BAND,
     &        PSIVAC,PSISEM)
            IF (IWKX.EQ.0.AND.IWKY.EQ.0) THEN
C
C   WRITE OUTPUT IF DESIRED
C
               IF (IWRIT.GE.4) THEN
                 do 50 J=NBARR2,1,-1
                    write(LUN+3,*) -(J-1)*DELVAC,PSIVAC(J)
50               CONTINUE
                 do 60 J=1,NS2
                    write(LUN+3,*) S2(J),PSISEM(J)
60               CONTINUE
               END IF
C
C   INCLUDE IN CHARGE DENSITY
C
               IF (ICOMP.EQ.1) THEN
               IF (IE.EQ.1) THEN
                  SUM3=0.
                  do 70 ieocc=ie,ne
                     enerocc=emin+ieocc*dele
                     OCCSEM2=FD(ENEROCC,EF,TK1)
                     SUM3=SUM3+OCCSEM2
70                continue
               ELSE
                  enerocc=emin+(ie-1)*dele
                  OCCSEM2=FD(ENEROCC,EF,TK1)
                  SUM3=SUM3-OCCSEM2
               END IF
               OCCINT=-SUM3*DELE
               CMULT=OCCINT*(EFFM*C/(2.*PI))*(DELE*EFFM*C/(2.*PI*WKSEM))
               SUM2=0.
               do 80 j=NBARR2,1,-1
                  TMP=(PSIVAC(J))**2*CMULT
                  IF (MOD(J-1,NEXVAC).EQ.0) THEN
                     JTMP=((J-1)/NEXVAC)+1
                     CDEVAC(JTMP)=CDEVAC(JTMP)+TMP
                  END IF
                  IF (J.NE.1) THEN
                     SUM2=SUM2+TMP*DELVAC
                  ELSE
                     SUM2=SUM2+TMP*0.5*DELVAC
                  END IF
80             CONTINUE
               CDESURF=CDESURF+SUM2
               do 90 j=1,NS2
                  TMP=(PSISEM(J))**2*CMULT
                  IF (JSEM(J).NE.0) THEN
                    CDESEM(JSEM(J))=CDESEM(JSEM(J))+
     &                                         TMP/NEXSEM(JSEM(J))
                  ELSE
                     WRITE(6,*) 'ERROR , INTCURR - J=0'
                     READ(5,*)
                  END IF
90             CONTINUE
               END IF
            END IF
C
C   EVALUATE EVERYTHING; eperp, kappa in vacuum
C
            EPERP=ener-WKPARR**2/C
            KAPPA=SQRT(C*(BARR2(NBARR2)-EPERP))
            trans=2.*nwkdeg*(2.*wf)**2*WKFTIP/(WKSEM/EFFM)
            sum=sum+trans*occdiff
110      continue
115      continue
120   continue
      currc=sum*dele*delwk**2/(4.*PI**2*RQUANT)
C
C   LOCALIZED STATES
C
200   IF (PMIN.GE.EC) RETURN
      emin=PMIN
c   could delete statements below!
      IF (NWK.NE.1) THEN
         emax=amin1(EC,amax1(ef+10.*tk1,ef+bias+10.*tk2))
      ELSE
         emax=amin1(EC,ef+10.*tk1)
      END IF
      dele=(emax-emin)/ne
      emax2=ef+10.*tk1
      dele2=(emax2-emin)/ne
      sum=0.
      if (dele.le.0.) go to 460
      wkmax=sqrt(C*EFFM*(emax-emin))
      SEMSTEP=(2.*pi/wkmax)/EXPANI
      KAPPA=SQRT(C*(amax1(BARR(NBARR1),BARR(1))-EMIN))
      VACSTEP=(2.*pi/KAPPA)/EXPANI
      Pot0p=Pot0+CBEDGE(0.)
      CALL POTEXPAND(IMPOT,SEP,NV,Pot0p,S,NSP,NSDIM,BARR,NBARR1,BARR2,
     &NBARR2,NVDIM1,NVDIM2,CBPROF,PROF2,NSDIM2,S2,NS2,VACSTEP,SEMSTEP,
     &JSEM,NEXSEM,NEXVAC,IWRIT)
      DELVAC=sep/float(NBARR2-1)
      delwk=wkmax/nwk
      do 320 iwky=0,nwk-1
         wky=iwky*delwk
         do 315 iwkx=0,nwk-1
         wkx=iwkx*delwk
         wkparr=sqrt(wkx**2+wky**2)
         eparr=wkparr**2/(EFFM*C)
         n=0
         nsav=0
         nwkdeg=8
         if (iwkx.eq.0) nwkdeg=nwkdeg/2
         if (iwky.eq.0) nwkdeg=nwkdeg/2
         if (iwkx.eq.iwky) nwkdeg=nwkdeg/2
         if (iwky.gt.iwkx) go to 315
         do 310 ie=1,ne
            ener=emin+(ie-0.5)*dele
            if (eparr.ge.(ener-emin)) go to 310
            occtip=fd(ener-bias,ef,tk2)
            occsem=fd(ener,ef,tk1)
            occdiff=occtip-occsem
            call CBloc(IMPOT,n,wf,wfderiv,ener,wkparr,sep,bias,BARR2,
     &          NVDIM2,NBARR2,PROF2,NS2,NSDIM2,S2,EFFM,EC,E2BAND,PSIVAC,
     &          PSISEM)
            if (n.eq.nsav) go to 310
            IF (IWKX.EQ.0.AND.IWKY.EQ.0) THEN
               IF (PSISEM(1).NE.0) THEN
                  NLOC=NLOC+1
C
C   FOUND LOCALIZED STATE; WRITE OUTPUT IF DESIRED
C
                  IF (IWRIT.GT.1) THEN
                     write(6,*) 'CB localized state at energy ',ener
                     write(16,*)'CB localized state at energy ',ener
                  END IF
                  IF (IWRIT.GE.1) THEN
                     WRITE(LUN,*) BIAS,N,ENER,OCCDIFF
                  END IF
                  IF (IWRIT.GE.2) THEN
                     do 250 J=NBARR2,1,-1
                        write(LUN+1,*) -(J-1)*DELVAC,PSIVAC(J)
250                  CONTINUE
                     do 260 J=1,NS2
                        write(LUN+1,*) S2(J),PSISEM(J)
260                  CONTINUE
                  END IF
C
C   INCLUDE IN CHARGE DENSITY
C
                  IF (ICOMP.EQ.1) THEN
                  OCCINT=-AMAX1(0.,(EF-ENER))
                  IF (TK1.NE.0.) THEN
                     SUM2=0.D0
                     DO 265 IEOCC=1,NE
                        ENEROCC=ENER+(IEOCC-0.5)*DELE2
                        IF (ENEROCC.GT.EMAX2) GO TO 270
                        OCCSEM2=FD(ENEROCC,EF,TK1)
                        SUM2=SUM2+OCCSEM2
265                  CONTINUE
270                  CONTINUE      
                     SUM2=SUM2*DELE2
                     OCCINT=-SUM2
                  END IF
                  CMULT=OCCINT*(EFFM*C/(2.*PI))
                  SUM2=0.
                  do 275 j=NBARR2,1,-1
                     TMP=(PSIVAC(J))**2*CMULT
                     IF (MOD(J-1,NEXVAC).EQ.0) THEN
                        JTMP=((J-1)/NEXVAC)+1
                        CDLVAC(JTMP)=CDLVAC(JTMP)+TMP
                     END IF
                     IF (IWRIT.GE.2) WRITE(LUN+2,*) -(J-1)*DELVAC,TMP
                     IF (J.NE.1) THEN
                        SUM2=SUM2+TMP*DELVAC
                     ELSE
                        SUM2=SUM2+TMP*0.5*DELVAC
                     END IF
275                 CONTINUE
                  CDLSURF=CDLSURF+SUM2
                  do 290 j=1,NS2
                     TMP=(PSISEM(J))**2*CMULT
                     IF (IWRIT.GE.2) WRITE(LUN+2,*) S2(J),TMP
                     IF (JSEM(J).NE.0) THEN
                       CDLSEM(JSEM(J))=CDLSEM(JSEM(J))+
     &                                            TMP/NEXSEM(JSEM(J))
                     ELSE
                     WRITE(6,*) 'ERROR , INTCURR - J=0'
                     READ(5,*)
                     END IF
290               CONTINUE
                  END IF
               END IF
            END IF
C
C   EVALUATE EVERYTHING; eperp, kappa in vacuum
C
            EPERP=ener-WKPARR**2/C
            KAPPA=SQRT(C*(BARR2(NBARR2)-EPERP))
            trans=nwkdeg*(n-nsav)*2.*(2.*wf)**2*WKFTIP
            sum=sum+trans*occdiff
            nsav=n
310      continue
315      continue
320   continue
      currc0=sum*delwk**2/(C*2.*PI*RQUANT)
460   RETURN
      END
c
c   integrate VB wavefunction from tip across to sample, return
c   values of wf and derivative of wf at tip surface
c
      subroutine VBwf(IMPOT,wf,wfderiv,WKSEM,e,wkparr,sep,bias,BARR2,
     &NVDIM2,NBARR2,PROF2,NS2,NSDIM2,S2,effm,EV,E2BAND,PSIVAC,PSISEM)
      dimension PROF2(NSDIM2),BARR2(NVDIM2),S2(NSDIM2)
      DOUBLE PRECISION PSI,DPSI,PSIVAC(NVDIM2),PSISEM(NSDIM2)
c   value of C below is 2m/hbar^2 in units of 1/(eV nm^2)
      DATA C/26.254/
      PI=4.*ATAN(1.)
c
      wf=0.
      wfderiv=0.
      WKSEM=1.
      eperp=e+wkparr**2/(C*effm)
      if (eperp.ge.EV) return
c
c   determine initial conditions for wavefunction
c
      eperp=e-wkparr**2/(C)
      PSI=1.D0
      PSIVAC(NBARR2)=PSI
      imin=1
      imax=NBARR2
      if (eperp.lt.BARR2(1).and.eperp.lt.BARR2(NBARR2)) go to 80
c   eperp above vacuum barrier; find crossing points
      do 30 i=1,NBARR2
         if (eperp.lt.BARR2(i)) go to 40
30    continue
40    imin=i
      do 50 i=NBARR2,1,-1
         PSIVAC(I)=PSI
         if (eperp.lt.BARR2(i)) go to 60
50    continue
60    imax=i
      if (imax.gt.imin) go to 80
      write(6,*) '*** error - eperp above vacuum barrier'
      write(16,*) '*** error - eperp above vacuum barrier'
      return
80    DPSI=PSI*sqrt(c*(BARR2(imax)-eperp))
      wf=psi
      wfderiv=dpsi
c
c   integrate through vacuum
c
      DELVAC=sep/float(NBARR2-1)
      do 100 i=imax-1,1,-1
         if (IMPOT.EQ.1.AND.(BARR2(i)-eperp).le.0.) go to 100
         PSI=PSI+DPSI*DELVAC
         PSIVAC(I)=PSI
         DPSI=DPSI+C*(BARR2(i)-eperp)*PSI*DELVAC
100   CONTINUE
c
c   match across vacuum-semiconductor interface
c
      PSI=PSI
      DPSI=DPSI*effm
c
c   integrate through semiconductor
c
      eperp=e+wkparr**2/(C*effm)
      PSI=PSI+DPSI*S2(1)
      PSISEM(1)=PSI
      ebarr=eperp-(PROF2(1))
C      if (ebarr.gt.0.) ebarr=ebarr-ebarr**2/E2BAND
      DPSI=DPSI+C*effm*(ebarr)*PSI*S2(1)
      do 200 i=2,NS2
         dels=(S2(I)-S2(I-1))
         PSI=PSI+DPSI*dels
         if (.not.(dabs(psi).lt.1.d100)) go to 500
         PSISEM(I)=PSI
         ebarr=eperp-(PROF2(i))
C         if (ebarr.gt.0.) ebarr=ebarr-ebarr**2/E2BAND
         DPSI=DPSI+C*effm*(ebarr)*PSI*dels
200   CONTINUE
c
c   determine amplitude, and evaluate across barrier
c
      WKSEM=sqrt(C*effm*(EV-eperp))
      phase=atan(psi*wksem/dpsi)
      amp=sqrt(2.)*sin(phase)/psi
      wf=wf*amp
      wfderiv=wfderiv*amp
C
C   NORMALIZE WAVEFUNCTION
C
      DO 250 I=NBARR2,1,-1
         PSIVAC(I)=PSIVAC(I)*AMP
250   CONTINUE         
      DO 300 I=1,NS2
         PSISEM(I)=PSISEM(I)*AMP
300   CONTINUE         
      delsmax=S2(NS2)-S2(NS2-1)
      if (delsmax/(2.*pi/wksem).GT.0.25) then
         write(6,*) '*** CAUTION *** RATIO OF SEMICOND. STEP SIZE ',
     &      'TO WAVELENGTH = ',delsmax/(2.*pi/wksem)
         write(16,*) '*** CAUTION *** RATIO OF SEMICOND. STEP SIZE ',
     &      'TO WAVELENGTH = ',delsmax/(2.*pi/wksem)
      end if
      return
C
C   OVERFLOW OF WAVEFUNCTION
C
500   wf=0.
      wfderiv=0.
      return
      end
c
c   integrate CB wavefunction from tip across to sample, return
c   values of wf and derivative of wf at tip surface
c
      subroutine CBwf(IMPOT,wf,wfderiv,WKSEM,e,wkparr,sep,bias,BARR2,
     &NVDIM2,NBARR2,PROF2,NS2,NSDIM2,S2,effm,EC,E2BAND,PSIVAC,PSISEM)
      dimension PROF2(NSDIM2),BARR2(NVDIM2),S2(NSDIM2)
      DOUBLE PRECISION PSI,DPSI,PSIVAC(NVDIM2),PSISEM(NSDIM2)
c   value of C below is 2m/hbar^2 in units of 1/(eV nm^2)
      DATA C/26.254/
      PI=4.*ATAN(1.)
c
      wf=0.
      wfderiv=0.
      WKSEM=1.
      eperp=e-wkparr**2/(C*effm)
      if (eperp.le.EC) return
c
c   determine initial conditions for wavefunction
c
      eperp=e-wkparr**2/(C)
      PSI=1.D0
      PSIVAC(NBARR2)=PSI
      imin=1
      imax=NBARR2
      if (eperp.lt.BARR2(1).and.eperp.lt.BARR2(NBARR2)) go to 80
c   eperp above vacuum barrier; find crossing points
      do 30 i=1,NBARR2
         if (eperp.lt.BARR2(i)) go to 40
30    continue
40    imin=i
      do 50 i=NBARR2,1,-1
         PSIVAC(I)=PSI
         if (eperp.lt.BARR2(i)) go to 60
50    continue
60    imax=i
      if (imax.gt.imin) go to 80
      write(6,*) '*** error - eperp above vacuum barrier'
      write(16,*) '*** error - eperp above vacuum barrier'
      return
80    DPSI=PSI*sqrt(c*(BARR2(imax)-eperp))
      wf=psi
      wfderiv=dpsi
c
c   integrate through vacuum
c
      DELVAC=sep/float(NBARR2-1)
      do 100 i=imax-1,1,-1
         if (IMPOT.EQ.1.AND.(BARR2(i)-eperp).le.0.) go to 100
         PSI=PSI+DPSI*DELVAC
         PSIVAC(I)=PSI
         DPSI=DPSI+C*(BARR2(i)-eperp)*PSI*DELVAC
100   CONTINUE
c
c   match across vacuum-semiconductor interface
c
      PSI=PSI
      DPSI=DPSI*effm
c
c   integrate through semiconductor
c
      eperp=e-wkparr**2/(C*effm)
      PSI=PSI+DPSI*S2(1)
      PSISEM(1)=PSI
      ebarr=(PROF2(1))-eperp
C      if (ebarr.gt.0.) ebarr=ebarr-ebarr**2/E2BAND
      DPSI=DPSI+C*effm*(ebarr)*PSI*S2(1)
      do 200 i=2,NS2
         dels=(S2(I)-S2(I-1))
         PSI=PSI+DPSI*dels
         if (.not.(dabs(psi).lt.1.d100)) go to 500
         PSISEM(I)=PSI
         ebarr=(PROF2(i))-eperp
C         if (ebarr.gt.0.) ebarr=ebarr-ebarr**2/E2BAND
         DPSI=DPSI+C*effm*(ebarr)*PSI*dels
200   CONTINUE
c
c   determine amplitude, and evaluate across barrier
c
      WKSEM=sqrt(C*effm*(eperp-EC))
      phase=atan(psi*wksem/dpsi)
      amp=sqrt(2.)*sin(phase)/psi
      wf=wf*amp
      wfderiv=wfderiv*amp
C
C   NORMALIZE WAVEFUNCTION
C
      DO 250 I=NBARR2,1,-1
         PSIVAC(I)=PSIVAC(I)*AMP
250   CONTINUE
      DO 300 I=1,NS2
         PSISEM(I)=PSISEM(I)*AMP
300   CONTINUE         
      delsmax=S2(NS2)-S2(NS2-1)
      if (delsmax/(2.*pi/wksem).GT.0.25) then
         write(6,*) '*** CAUTION *** RATIO OF SEMICOND. STEP SIZE ',
     &      'TO WAVELENGTH = ',delsmax/(2.*pi/wksem)
         write(16,*) '*** CAUTION *** RATIO OF SEMICOND. STEP SIZE ',
     &      'TO WAVELENGTH = ',delsmax/(2.*pi/wksem)
      end if
      return
C
C   OVERFLOW OF WAVEFUNCTION
C
500   wf=0.
      wfderiv=0.
      return
      end
c
c   integrate localized VB wavefunction from tip across to sample,
c   looking for switches in sign of wf to enumerate localized states
c
      subroutine VBloc(IMPOT,nsign,wf,wfderiv,e,wkparr,sep,bias,BARR2,
     &NVDIM2,NBARR2,PROF2,NS2,NSDIM2,S2,effm,EV,E2BAND,PSIVAC,PSISEM)
      dimension PROF2(NSDIM2),BARR2(NVDIM2),S2(NSDIM2)
      DOUBLE PRECISION PSI,DPSI,SUM1,SUM2,PSIVAC(NVDIM2),PSISEM(NSDIM2)
c   value of C below is 2m/hbar^2 in units of 1/(eV nm^2)
      DATA C/26.254/
c
      wf=0.
      wfderiv=0.
      isav=NS2
      eperp=e+wkparr**2/(C*effm)
      if (eperp.le.EV) return
      nsign=0
c
c   determine initial conditions for wavefunction
c
      eperp=e-wkparr**2/(C)
      SUM1=0.D0
      SUM2=0.D0
      PSI=1.D0
      PSIVAC(NBARR2)=PSI
      imin=1
      imax=NBARR2
      if (eperp.lt.BARR2(1).and.eperp.lt.BARR2(NBARR2)) go to 80
c   eperp above vacuum barrier; find crossing points
      do 30 i=1,NBARR2
         if (eperp.lt.BARR2(i)) go to 40
30    continue
40    imin=i
      do 50 i=NBARR2,1,-1
         PSIVAC(I)=PSI
         if (eperp.lt.BARR2(i)) go to 60
50    continue
60    imax=i
      if (imax.gt.imin) go to 80
      write(6,*) '*** error - eperp above vacuum barrier'
      write(16,*) '*** error - eperp above vacuum barrier'
      return
80    DPSI=PSI*sqrt(c*(BARR2(imax)-eperp))
      wf=psi
      wfderiv=dpsi
c
c   integrate through vacuum
c
      DELVAC=sep/float(NBARR2-1)
      do 100 i=imax-1,1,-1
         if (IMPOT.EQ.1.AND.(BARR2(i)-eperp).le.0.) go to 100
         PSI=PSI+DPSI*DELVAC
         PSIVAC(I)=PSI
         SUM1=SUM1+PSI**2*DELVAC
         DPSI=DPSI+C*(BARR2(i)-eperp)*PSI*DELVAC
100   CONTINUE
c
c   match across vacuum-semiconductor interface
c
      PSI=PSI
      DPSI=DPSI*effm
c
c   integrate through semiconductor
c
      eperp=e+wkparr**2/(C*effm)
      PSI=PSI+DPSI*S2(1)
      PSISEM(1)=PSI
      SUM1=SUM1+PSI**2*S2(1)
      ebarr=eperp-(PROF2(1))
C      if (ebarr.gt.0.) ebarr=ebarr-ebarr**2/E2BAND
      DPSI=DPSI+C*effm*(ebarr)*PSI*S2(1)
      do 200 i=2,NS2-1
         dels=(S2(I)-S2(I-1))
         PSISAV=PSI
         PSI=PSI+DPSI*dels
         PSISEM(I)=PSI
         if (PSISAV*PSI.lt.0.) then
            nsign=nsign+1
            isav=i
            SUM2=SUM2+SUM1
            SUM1=0.D0
         end if
         SUM1=SUM1+PSI**2*dels
         ebarr=eperp-(PROF2(i))
C         if (ebarr.gt.0.) ebarr=ebarr-ebarr**2/E2BAND
         DPSI=DPSI+C*effm*(ebarr)*PSI*dels
200   CONTINUE
C
C   NORMALIZE WAVEFUNCTION
C
      IF (SUM2.NE.0.) THEN
         AMP=1./SQRT(SUM2)
      ELSE
         AMP=1./SQRT(SUM1)
      END IF
      wf=wf*AMP
      wfderiv=wfderiv*AMP
      DO 250 I=NBARR2,1,-1
         PSIVAC(I)=PSIVAC(I)*AMP
250   CONTINUE         
      DO 300 I=1,ISAV-1
         PSISEM(I)=PSISEM(I)*AMP
300   CONTINUE         
      do 350 i=isav,NS2
         psisem(i)=0.
350   continue
c
      return
      end
c
c   integrate localized CB wavefunction from tip across to sample,
c   looking for switches in sign of wf to enumerate localized states
c
      subroutine CBloc(IMPOT,nsign,wf,wfderiv,e,wkparr,sep,bias,BARR2,
     &NVDIM2,NBARR2,PROF2,NS2,NSDIM2,S2,effm,EC,E2BAND,PSIVAC,PSISEM)
      dimension PROF2(NSDIM2),BARR2(NVDIM2),S2(NSDIM2)
      DOUBLE PRECISION PSI,DPSI,SUM1,SUM2,PSIVAC(NVDIM2),PSISEM(NSDIM2)
c   value of C below is 2m/hbar^2 in units of 1/(eV nm^2)
      DATA C/26.254/
c
      wf=0.
      wfderiv=0.
      isav=NS2
      eperp=e-wkparr**2/(C*effm)
      if (eperp.ge.EC) return
      nsign=0
c
c   determine initial conditions for wavefunction
c
      eperp=e-wkparr**2/(C)
      SUM1=0.D0
      SUM2=0.D0
      PSI=1.D0
      PSIVAC(NBARR2)=PSI
      imin=1
      imax=NBARR2
      if (eperp.lt.BARR2(1).and.eperp.lt.BARR2(NBARR2)) go to 80
c   eperp above vacuum barrier; find crossing points
      do 30 i=1,NBARR2
         if (eperp.lt.BARR2(i)) go to 40
30    continue
40    imin=i
      do 50 i=NBARR2,1,-1
         PSIVAC(I)=PSI
         if (eperp.lt.BARR2(i)) go to 60
50    continue
60    imax=i
      if (imax.gt.imin) go to 80
      write(6,*) '*** error - eperp above vacuum barrier'
      write(16,*) '*** error - eperp above vacuum barrier'
      return
80    DPSI=PSI*sqrt(c*(BARR2(imax)-eperp))
      wf=psi
      wfderiv=dpsi
c
c   integrate through vacuum
c
      DELVAC=sep/float(NBARR2-1)
      do 100 i=imax-1,1,-1
         if (IMPOT.EQ.1.AND.(BARR2(i)-eperp).le.0.) go to 100
         PSI=PSI+DPSI*DELVAC
         PSIVAC(I)=PSI
         SUM1=SUM1+PSI**2*DELVAC
         DPSI=DPSI+C*(BARR2(i)-eperp)*PSI*DELVAC
100   CONTINUE
c
c   match across vacuum-semiconductor interface
c
      PSI=PSI
      DPSI=DPSI*effm
c
c   integrate through semiconductor
c
      eperp=e-wkparr**2/(C*effm)
      PSI=PSI+DPSI*S2(1)
      PSISEM(1)=PSI
      SUM1=SUM1+PSI**2*S2(1)
      ebarr=(PROF2(1))-eperp
C      if (ebarr.gt.0.) ebarr=ebarr-ebarr**2/E2BAND
      DPSI=DPSI+C*effm*(ebarr)*PSI*S2(1)
      do 200 i=2,NS2-1
         dels=(S2(I)-S2(I-1))
         PSISAV=PSI
         PSI=PSI+DPSI*dels
         PSISEM(I)=PSI
         if (PSISAV*PSI.lt.0.) then
            nsign=nsign+1
            isav=i
            SUM2=SUM2+SUM1
            SUM1=0.D0
         end if
         SUM1=SUM1+PSI**2*dels
         ebarr=(PROF2(i))-eperp
C         if (ebarr.gt.0.) ebarr=ebarr-ebarr**2/E2BAND
         DPSI=DPSI+C*effm*(ebarr)*PSI*dels
200   CONTINUE
C
C   NORMALIZE WAVEFUNCTION
C
      IF (SUM2.NE.0.) THEN
         AMP=1./SQRT(SUM2)
      ELSE
         AMP=1./SQRT(SUM1)
      END IF
      wf=wf*AMP
      wfderiv=wfderiv*AMP
      DO 250 I=NBARR2,1,-1
         PSIVAC(I)=PSIVAC(I)*AMP
250   CONTINUE         
      DO 300 I=1,ISAV-1
         PSISEM(I)=PSISEM(I)*AMP
300   CONTINUE         
      do 350 i=isav,NS2
         psisem(i)=0.
350   continue
c
      return
      end
