import numpy as np
from potperiod3 import potperiod
from getcurr import getcurr

def planecurr(SEP, VAC, TIP, SEM, VSINT, R, S, DELV, DELR, DELS, DELP,
              NRDIM, NVDIM, NSDIM, NPDIM, NR, NV, NS, NP, NXDIM, NXDIM2, NYDIM,
              NZDIM, PVAC, PSEM, PSURF, VACWID, CHI, EFTIP, EGAP, AVBL, AVBH, AVBSO, ACB,
              ESO, BIAS, DELPHI, PHI0, TK, EF, EMAX, CURR, CURRV, CURRC, IWRIT, IERR,
              NBARR2, NVDIM2, BARRPROF, NEIGENDIM, ZVACDEL, ICOMP):
    # Constants
    C = 26.254
    PI = np.pi

    CURRL = CURRH = CURRSO = CURRC = CURRV = 0.0

    if ICOMP != 0:
        if IWRIT > 0:
            print('*********** LIGHT-HOLE VALENCE BAND *************')
        
        EFFM = AVBL
        BARR = CHI + EGAP + DELPHI / 2. + BIAS / 2. + PHI0 / 2.
        BARR1 = BARR - PHI0
        
        EPSI2, NSTATES = getstates(-1, SEP, VAC, TIP, SEM, VSINT, R, S, DELV, DELR, DELS, DELP,
                                   NRDIM, NVDIM, NSDIM, NPDIM, NR, NV, NS, NP, NXDIM, NXDIM2, NYDIM, NZDIM,
                                   PVAC, PSEM, PSURF, VACWID, BARR1, EFFM, EMAX[0], IWRIT, IERR, NEIGENDIM, ZVACDEL)
        
        E0 = 0.0
        BARR2 = BARR + E0
        CURRL = getcurr(-1, BIAS, SEP, E0, BARR2, TK, TK, EF, EPSI2, NSTATES, EFTIP,
                        NBARR2, NVDIM2, BARRPROF, NEIGENDIM)
        
        if IWRIT > 0:
            print('LIGHT-HOLE CURRENT =', CURRL)

    if IWRIT > 0:
        print('*********** HEAVY-HOLE VALENCE BAND *************')
    
    EFFM = AVBH
    BARR = CHI + EGAP + DELPHI / 2. + BIAS / 2. + PHI0 / 2.
    BARR1 = BARR - PHI0
    
    EPSI2, NSTATES = getstates(-1, SEP, VAC, TIP, SEM, VSINT, R, S, DELV, DELR, DELS, DELP,
                               NRDIM, NVDIM, NSDIM, NPDIM, NR, NV, NS, NP, NXDIM, NXDIM2, NYDIM, NZDIM,
                               PVAC, PSEM, PSURF, VACWID, BARR1, EFFM, EMAX[1], IWRIT, IERR, NEIGENDIM, ZVACDEL)
    
    E0 = 0.0
    BARR2 = BARR + E0
    CURRH = getcurr(-1, BIAS, SEP, E0, BARR2, TK, TK, EF, EPSI2, NSTATES, EFTIP,
                    NBARR2, NVDIM2, BARRPROF, NEIGENDIM)
    
    if IWRIT > 0:
        print('HEAVY-HOLE CURRENT =', CURRH)

    if IWRIT > 0:
        print('*********** SPLIT-OFF VALENCE BAND *************')
    
    EFFM = AVBSO
    BARR = CHI + EGAP + DELPHI / 2. + BIAS / 2. + PHI0 / 2.
    BARR1 = BARR - PHI0
    
    EPSI2, NSTATES = getstates(-1, SEP, VAC, TIP, SEM, VSINT, R, S, DELV, DELR, DELS, DELP,
                               NRDIM, NVDIM, NSDIM, NPDIM, NR, NV, NS, NP, NXDIM, NXDIM2, NYDIM, NZDIM,
                               PVAC, PSEM, PSURF, VACWID, BARR1, EFFM, EMAX[2], IWRIT, IERR, NEIGENDIM, ZVACDEL)
    
    E0 = -ESO
    BARR2 = BARR + E0
    CURRSO = getcurr(-1, BIAS, SEP, E0, BARR2, TK, TK, EF, EPSI2, NSTATES, EFTIP,
                     NBARR2, NVDIM2, BARRPROF, NEIGENDIM)
    
    if IWRIT > 0:
        print('SPLIT-OFF CURRENT =', CURRSO)
        CURRV = CURRL + CURRH + CURRSO
        print('TOTAL HOLE CURRENT =', CURRV)

    if ICOMP != -1:
        if IWRIT > 0:
            print('*********** CONDUCTION BAND *************')
        
        EFFM = ACB
        BARR = CHI + DELPHI / 2. + BIAS / 2. + PHI0 / 2.
        BARR1 = BARR - PHI0
        
        EPSI2, NSTATES = getstates(1, SEP, VAC, TIP, SEM, VSINT, R, S, DELV, DELR, DELS, DELP,
                                   NRDIM, NVDIM, NSDIM, NPDIM, NR, NV, NS, NP, NXDIM, NXDIM2, NYDIM, NZDIM,
                                   PVAC, PSEM, PSURF, VACWID, BARR1, EFFM, EMAX[3], IWRIT, IERR, NEIGENDIM, ZVACDEL)
        
        E0 = EGAP
        BARR2 = BARR + E0
        CURRC = getcurr(1, BIAS, SEP, E0, BARR2, TK, TK, EF, EPSI2, NSTATES, EFTIP,
                        NBARR2, NVDIM2, BARRPROF, NEIGENDIM)
        
        if IWRIT > 0:
            print('ELECTRON CURRENT =', CURRC)
    
    CURR = CURRV + CURRC
    return CURR


def getstates(IBAND, SEP, VAC, TIP, SEM, VSINT, R, S, DELV, DELR, DELS, DELP,
              NRDIM, NVDIM, NSDIM, NPDIM, NR, NV, NS, NP, NXDIM, NXDIM2, NYDIM, NZDIM,
              PVAC, PSEM, PSURF, VACWID, BARR, EFFM, EMAX, IWRIT, IERR, NEIGENDIM, ZVACDEL):
    # Constants
    NKX = 6
    NKY = 6
    NKZ = 6
    NKXM1 = NKX - 1
    NKX2M1 = 2 * NKX - 1
    NK3 = NKX2M1 * NKY * NKZ
    NZVACDIM = 100
    C = 26.254
    PI = np.pi
    
    plot1 = np.zeros((NXDIM2, NYDIM))
    plot2 = np.zeros((NXDIM2, NYDIM))
    plot3 = np.zeros((NXDIM2, NYDIM))
    plot4 = np.zeros((NXDIM2, NYDIM))
    plot5 = np.zeros((NXDIM2, NYDIM))
    plot6 = np.zeros((NXDIM2, NYDIM))
    plot7 = np.zeros((NXDIM2, NYDIM))
    plot8 = np.zeros((NXDIM2, NYDIM))
    plot9 = np.zeros((NXDIM2, NYDIM))
    
    # Initializing matrices and arrays using NumPy
    CMAT = np.zeros((NK3, NK3))
    W = np.zeros(NK3)
    ZMAT = np.zeros((NK3, NK3))
    ENER = np.zeros(1000)
    AMP1 = np.zeros(1000)
    AMP2 = np.zeros(1000)
    COSX = np.zeros((NKX, NXDIM2))
    SINX = np.zeros((NKXM1, NXDIM2))
    COSY = np.zeros((NKY, NYDIM))
    HARM = np.zeros((NKX, NKY, NKZ, NZDIM))
    TUNN = np.zeros((NKX, NKY, NKZ, NZVACDIM))
    WQ = np.zeros((NKX, NKY, NKZ))
    WKAP = np.zeros((NKX, NKY, NKZ))
    AM1 = np.zeros((NKX, NKY, NKZ))
    AM2 = np.zeros((NKX, NKY, NKZ))

    # Parameters
    IPOT = 1
    IVAC = 1
    DEFFM = EFFM

    if IWRIT > 0:
        print('ENTERING GETstates')
        print('EFFM =', EFFM)
        if IPOT == 0:
            print('ZEROING ELECTROSTATIC POTENTIAL')
        if IVAC == 0:
            print('NO VACUUM REGION; FULL QUANTUM DOT')

    NX = 2 * NKX
    NY = 2 * NKY
    NZ = 2 * NKZ
    if IWRIT > 0:
        print('energy cutoff =', EMAX)
        print('number of k-points =', NKX, NKY, NKZ)

    XMAX = NKX * PI / np.sqrt(C * EFFM * EMAX)
    YMAX = NKY * PI / np.sqrt(C * EFFM * EMAX)
    ZMAX = NKZ * PI / np.sqrt(C * EFFM * EMAX)
    EL1 = ZMAX * 2
    EL2 = EL1 + VACWID
    if IWRIT > 0:
        print(f'slab size X,Y,Z =   {XMAX * 2:.9f} {YMAX * 2:.9f} {ZMAX * 2:.9f}')
        print(f'grid size X,Y,Z =  {XMAX / (NX - 1):.9f} {YMAX / (NY - 1):.9f} {ZMAX / (NZ - 1):.9f}')
        print('*** ODD Z PARITY STATES')

    XDEL = XMAX / (NX - 1)
    YDEL = YMAX / (NY - 1)
    ZDEL = ZMAX / (NZ - 1)
    
    potperiod(SEP, NX, NY, NZ, XDEL, YDEL, ZDEL,
              VAC, TIP, SEM, VSINT, PVAC, PSEM, R, S, DELV, DELR, DELS, DELP, NRDIM,
              NVDIM, NSDIM, NPDIM, NXDIM, NXDIM2, NYDIM, NZDIM, NR, NV, NS, NP, IWRIT, IERR)

    if NX > NXDIM or NY > NYDIM or NZ > NZDIM:
        print('*** ERROR - INTERPOLATED POTENTIAL ARRAY TOO SMALL')
        print(NX, NY, NZ, NXDIM, NYDIM, NZDIM)
        raise ValueError('Interpolated potential array too small')

    if NKX > (NX - 1) or NKY > (NY - 1) or NKZ > (NZ - 2):
        print('*** ERROR - NUMBER OF K-POINTS TOO LARGE')
        raise ValueError('Number of k-points too large')

    ELX = XMAX * 2
    ELY = YMAX * 2
    ELZ1 = EL1
    ELZ2 = EL2
    EL = EL2 / 2.0
    DELVBEDGE = np.zeros((NXDIM, NYDIM, NZDIM))
    DELCBEDGE = np.zeros((NXDIM, NYDIM, NZDIM))

    for K in range(1, NZ + 1):
        for J in range(1, NY + 1):
            for I in range(1, 2 * NX - 1):
                if IPOT == 0:
                    PSEM[I][J][K] = 0.0
                if K == 1:
                    PSURF[I][J] = PSEM[I][J][K]
                if IBAND == -1:
                    PSEM[I][J][K] += DELVBEDGE[I][J][K]
                else:
                    PSEM[I][J][K] += DELCBEDGE[I][J][K]

                if IWRIT >= 7:
                    if K == 1: plot1[I][J] = PSEM[I][J][K]
                    if K == 2: plot2[I][J] = PSEM[I][J][K]
                    if K == 3: plot3[I][J] = PSEM[I][J][K]
                    if K == 4: plot4[I][J] = PSEM[I][J][K]
                    if K == 5: plot5[I][J] = PSEM[I][J][K]
                    if K == 6: plot6[I][J] = PSEM[I][J][K]
                    if K == 7: plot7[I][J] = PSEM[I][J][K]
                    if K == 8: plot8[I][J] = PSEM[I][J][K]
                    if K == 9: plot9[I][J] = PSEM[I][J][K]

    WKX = np.zeros(NKX)
    WKY = np.zeros(NKY)
    for KX in range(1, NKX + 1):
        for I in range(1, 2 * NX - 1):
            if KX == 1:
                COSX[KX - 1][I - 1] = 1.0 / np.sqrt(2.0 * XMAX)
            else:
                COSX[KX - 1][I - 1] = np.cos((KX - 1) * PI * (I - NX + 1) / (NX - 1)) / np.sqrt(XMAX)
                SINX[KX - 2][I - 1] = np.sin((KX - 1) * PI * (I - NX + 1) / (NX - 1)) / np.sqrt(XMAX)
        WKX[KX - 1] = (KX - 1) * PI / XMAX

    for KY in range(1, NKY + 1):
        for J in range(1, NY + 1):
            if KY == 1:
                COSY[KY - 1][J - 1] = 1.0 / np.sqrt(2.0 * YMAX)
            else:
                COSY[KY - 1][J - 1] = np.cos((KY - 1) * PI * (J - 1) / (NY - 1)) / np.sqrt(YMAX)
        WKY[KY - 1] = (KY - 1) * PI / YMAX

    ZVACDEL = 0.05
    NZVAC = int((ELZ2 - ELZ1) / (2.0 * ZVACDEL))
    if NZVAC > NZVACDIM:
        print('*** ERROR - NZVACDIM TOO SMALL, NZVAC =', NZVAC)
        raise ValueError('NZVACDIM TOO SMALL, NZVAC =', NZVAC)
    ZVACDEL = (ELZ2 - ELZ1) / (2.0 * (NZVAC - 1))
    if IWRIT > 0:
        print(f'VACUUM NZ,DELZ =  {NZVAC} {ZVACDEL:.9f}')

    IPARITY = 0
    for KX in range(1, NKX + 1):
        for KY in range(1, NKY + 1):
            for KZ in range(1, NKZ + 1):
                idx = KX + (KY - 1) * NKX + (KZ - 1) * NKX * NKY - 1
                for J in range(1, NZ + 1):
                    if IPARITY == 0:
                        HARM[KX - 1][KY - 1][KZ - 1][J - 1] = AMP2[idx] * np.sin(WQ[KX - 1][KY - 1][KZ - 1] * (NZ - J) * ZDEL)
                    else:
                        HARM[KX - 1][KY - 1][KZ - 1][J - 1] = AMP2[idx] * np.cos(WQ[KX - 1][KY - 1][KZ - 1] * (NZ - J) * ZDEL)
                for J in range(1, NZVAC + 1):
                    Z = (J - 1) * ZVACDEL + (ELZ1 / 2.0)
                    if WKAP[KX - 1][KY - 1][KZ - 1] > 0.0:
                        if IPARITY == 0:
                            TUNN[KX - 1][KY - 1][KZ - 1][J - 1] = AMP1[idx] * np.sinh(WKAP[KX - 1][KY - 1][KZ - 1] * (EL - Z))
                        else:
                            TUNN[KX - 1][KY - 1][KZ - 1][J - 1] = AMP1[idx] * np.cosh(WKAP[KX - 1][KY - 1][KZ - 1] * (EL - Z))
                    else:
                        if IPARITY == 0:
                            TUNN[KX - 1][KY - 1][KZ - 1][J - 1] = AMP1[idx] * np.sin(-WKAP[KX - 1][KY - 1][KZ - 1] * (EL - Z))
                        else:
                            TUNN[KX - 1][KY - 1][KZ - 1][J - 1] = AMP1[idx] * np.cos(-WKAP[KX - 1][KY - 1][KZ - 1] * (EL - Z))

    if IWRIT > 0:
        print('COMPUTING MATRIX ELEMENTS')

    CMAT = calculate_matrix_elements(COSX, SINX, COSY, HARM, PSEM, NKX, NKY, NKZ, NX, NY, NZ, NXDIM, NXDIM2, NZDIM, IBAND, WKX, WKY, WQ, C, EFFM)

    if IWRIT > 0:
        print('SOLVING EIGENVALUE PROBLEM')

    # 使用 NumPy 的特徵值計算
    W, ZMAT = np.linalg.eigh(CMAT)

    if IBAND == -1:
        if IWRIT > 0:
            print('highest eigenvalues =')
            for eigenvalue in W[-6:]:
                print(f'{eigenvalue:.15f}')
    else:
        if IWRIT > 0:
            print('lowest eigenvalues =')
            for eigenvalue in W[:6]:
                print(f'{eigenvalue:.15f}')

    EPSI2 = np.zeros((4, NEIGENDIM))
    NSTATES = 0

    if IWRIT >= 5:
        IWF1 = 1
        IWF2 = 2
        IWF3 = 3
        IWF4 = 4
        IWF5 = 5
        for J in range(NZ - 1, -1, -1):
            Z = J * ZDEL
            print(Z, PSEM[0][0][J])
            SUM1 = SUM2 = SUM3 = SUM4 = SUM5 = 0.0
            for K in range(NK3):
                KZ = (K // (NKX2M1 * NKY)) + 1
                KY = ((K % (NKX2M1 * NKY)) // NKX2M1) + 1
                KX = (K % NKX2M1) + 1
                if KX <= NKX:
                    SUM1 += ZMAT[K][IWF1 - 1] * COSX[KX - 1][NX - 1] * COSY[KY - 1][1] * HARM[KX - 1][KY - 1][KZ - 1][J]
                    SUM2 += ZMAT[K][IWF2 - 1] * COSX[KX - 1][NX - 1] * COSY[KY - 1][1] * HARM[KX - 1][KY - 1][KZ - 1][J]
                    SUM3 += ZMAT[K][IWF3 - 1] * COSX[KX - 1][NX - 1] * COSY[KY - 1][1] * HARM[KX - 1][KY - 1][KZ - 1][J]
                    SUM4 += ZMAT[K][IWF4 - 1] * COSX[KX - 1][NX - 1] * COSY[KY - 1][1] * HARM[KX - 1][KY - 1][KZ - 1][J]
                    SUM5 += ZMAT[K][IWF5 - 1] * COSX[KX - 1][NX - 1] * COSY[KY - 1][1] * HARM[KX - 1][KY - 1][KZ - 1][J]
                else:
                    SUM1 += ZMAT[K][IWF1 - 1] * SINX[KX - NKX - 1][NX - 1] * COSY[KY - 1][1] * HARM[KX - 1][KY - 1][KZ - 1][J]
                    SUM2 += ZMAT[K][IWF2 - 1] * SINX[KX - NKX - 1][NX - 1] * COSY[KY - 1][1] * HARM[KX - 1][KY - 1][KZ - 1][J]
                    SUM3 += ZMAT[K][IWF3 - 1] * SINX[KX - NKX - 1][NX - 1] * COSY[KY - 1][1] * HARM[KX - 1][KY - 1][KZ - 1][J]
                    SUM4 += ZMAT[K][IWF4 - 1] * SINX[KX - NKX - 1][NX - 1] * COSY[KY - 1][1] * HARM[KX - 1][KY - 1][KZ - 1][J]
                    SUM5 += ZMAT[K][IWF5 - 1] * SINX[KX - NKX - 1][NX - 1] * COSY[KY - 1][1] * HARM[KX - 1][KY - 1][KZ - 1][J]
            print(SUM1, SUM2, SUM3, SUM4, SUM5)
    return EPSI2, NSTATES
    """
    for K in range(NK3):
        if W[K] > EMAX:
            break
        NSTATES += 1
        if NSTATES > 4 * NEIGENDIM:
            print(f'*** ERROR - TOO MANY STATES, NSTATES = {NSTATES}')
            raise ValueError(f'Too many states, NSTATES = {NSTATES}')
        EPSI2[0][K] = W[K]
    if IWRIT > 0:
        print(f'NSTATES = {NSTATES}')
    """
     

def calculate_matrix_elements(COSX, SINX, COSY, HARM, PSEM, NKX, NKY, NKZ, NX, NY, NZ, NXDIM, NXDIM2, NZDIM, IBAND, WKX, WKY, WQ, C, EFFM):
    NK3 = NKX * NKY * NKZ
    CMAT = np.zeros((NK3, NK3))
    XDEL = 1.0  # Example value, replace with actual calculation
    YDEL = 1.0  # Example value, replace with actual calculation
    ZDEL = 1.0  # Example value, replace with actual calculation

    for KZP in range(1, NKZ + 1):
        for KYP in range(1, NKY + 1):
            for KXP in range(1, NKX + 1):
                KKP = KXP + (KYP - 1) * (2 * NKX - 1) + (KZP - 1) * (2 * NKX - 1) * NKY
                if KKP > NK3:
                    continue
                for KZ in range(1, NKZ + 1):
                    for KY in range(1, NKY + 1):
                        for KX in range(1, NKX + 1):
                            KK = KX + (KY - 1) * (2 * NKX - 1) + (KZ - 1) * (2 * NKX - 1) * NKY
                            if KK > NK3:
                                continue
                            if KK > KKP:
                                continue
                            SUM1 = SUM2 = SUM3 = SUM4 = 0.0
                            for K in range(1, NZ + 1):
                                for J in range(1, NY + 1):
                                    for I in range(1, 2 * NX - 1):
                                        TMP1 = COSX[KXP - 1, I - 1] * COSY[KYP - 1, J - 1] * HARM[KXP - 1, KYP - 1, KZP - 1, K - 1] * PSEM[I - 1, J - 1, K - 1] * COSX[KX - 1, I - 1] * COSY[KY - 1, J - 1] * HARM[KX - 1, KY - 1, KZ - 1, K - 1]
                                        TMP2 = TMP3 = TMP4 = 0.0
                                        if KXP != 1:
                                            TMP2 = SINX[KXP - 2, I - 1] * COSY[KYP - 1, J - 1] * HARM[KXP - 1, KYP - 1, KZP - 1, K - 1] * PSEM[I - 1, J - 1, K - 1] * COSX[KX - 1, I - 1] * COSY[KY - 1, J - 1] * HARM[KX - 1, KY - 1, KZ - 1, K - 1]
                                        if KX != 1:
                                            TMP3 = COSX[KXP - 1, I - 1] * COSY[KYP - 1, J - 1] * HARM[KXP - 1, KYP - 1, KZP - 1, K - 1] * PSEM[I - 1, J - 1, K - 1] * SINX[KX - 2, I - 1] * COSY[KY - 1, J - 1] * HARM[KX - 1, KY - 1, KZ - 1, K - 1]
                                        if KXP != 1 and KX != 1:
                                            TMP4 = SINX[KXP - 2, I - 1] * COSY[KYP - 1, J - 1] * HARM[KXP - 1, KYP - 1, KZP - 1, K - 1] * PSEM[I - 1, J - 1, K - 1] * SINX[KX - 2, I - 1] * COSY[KY - 1, J - 1] * HARM[KX - 1, KY - 1, KZ - 1, K - 1]
                                        if J != 1 and J != NY:
                                            TMP1 *= 2
                                            TMP2 *= 2
                                            TMP3 *= 2
                                            TMP4 *= 2
                                        if K != 1 and K != NZ:
                                            TMP1 *= 2
                                            TMP2 *= 2
                                            TMP3 *= 2
                                            TMP4 *= 2
                                        SUM1 += TMP1
                                        if KXP != 1:
                                            SUM2 += TMP2
                                        if KX != 1:
                                            SUM3 += TMP3
                                        if KXP != 1 and KX != 1:
                                            SUM4 += TMP4
                            CMAT[KKP - 1, KK - 1] = SUM1 * XDEL * YDEL * ZDEL
                            if KXP != 1:
                                if KKP - 1 + NKX - 1 < NK3 and KK - 1 < NK3:
                                    CMAT[KKP - 1 + NKX - 1, KK - 1] = SUM2 * XDEL * YDEL * ZDEL
                            if KX != 1:
                                if KKP - 1 < NK3 and KK - 1 + NKX - 1 < NK3:
                                    CMAT[KKP - 1, KK - 1 + NKX - 1] = SUM3 * XDEL * YDEL * ZDEL
                            if KX != 1 and KXP != 1:
                                if KKP - 1 + NKX - 1 < NK3 and KK - 1 + NKX - 1 < NK3:
                                    CMAT[KKP - 1 + NKX - 1, KK - 1 + NKX - 1] = SUM4 * XDEL * YDEL * ZDEL
                            if KKP == KK:
                                CMAT[KKP - 1, KK - 1] += IBAND * (WKX[KX - 1]**2 + WKY[KY - 1]**2 + WQ[KX - 1, KY - 1, KZ - 1]**2) / (C * EFFM)
                                if KX != 1 and KXP != 1:
                                    if KKP - 1 + NKX - 1 < NK3 and KK - 1 + NKX - 1 < NK3:
                                        CMAT[KKP - 1 + NKX - 1, KK - 1 + NKX - 1] += IBAND * (WKX[KX - 1]**2 + WKY[KY - 1]**2 + WQ[KX - 1, KY - 1, KZ - 1]**2) / (C * EFFM)

    return CMAT
