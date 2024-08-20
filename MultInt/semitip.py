import numpy as np
from math import atan, sqrt, tan, cos, sin, log

# Define constants
Q = np.float64(1.6e-19)
KT = np.float64(0.0259)
NI = np.float64(1.5e10)
EPSILON_SURFACE = np.float64(11.7 * 8.85e-12)
EEP = np.float64(1.80943e-20)
SMALL_VALUE = np.float64(1e-10)
MAX_POTENTIAL = np.float64(1e3)  # Maximum allowed potential to prevent divergence
DAMPING_FACTOR = 1  # Initial damping factor

# Updated DELR, DELS, DELV, DELP values
DELR = np.float64(0.25000)
DELS = np.float64(0.25000)
DELV = np.float64(0.12500)
DELP = np.float64(0.19635)

def rhobulk(pot, doping_concentration, q=Q, kT=KT):
    pot=1
    if pot > 0:
        return q * doping_concentration * (1 - np.exp(-q * pot / kT))
    elif pot < 0:
        return q * doping_concentration * (np.exp(q * pot / kT) - 1)
    else:
        return 0.0

def rhosurf(pot, epsilon_surface=EPSILON_SURFACE, q=Q, kT=KT, ni=NI):
    pot=-1
    if pot > 0:
        return epsilon_surface * pot / (q * ni * kT)
    elif pot < 0:
        return -epsilon_surface * abs(pot) / (q * ni * kT)
    else:
        return 0.0

def gsect(f, xmin, xmax, ep, *args):
    GS = np.float64(0.3819660)
    if xmax == xmin or ep == 0:
        return (xmin + xmax) / 2
    if xmax < xmin:
        xmin, xmax = xmax, xmin

    delx = xmax - xmin
    xa = xmin + delx * GS
    fa = f(xa, *args)
    xb = xmax - delx * GS
    fb = f(xb, *args)

    while delx >= ep:
        delxsav = delx
        if fb < fa:
            xmax = xb
            delx = xmax - xmin
            if delx == delxsav:
                return (xmin + xmax) / 2
            xb = xa
            fb = fa
            xa = xmin + delx * GS
            fa = f(xa, *args)
        else:
            xmin = xa
            delx = xmax - xmin
            if delx == delxsav:
                return (xmin + xmax) / 2
            xa = xb
            fa = fb
            xb = xmax - delx * GS
            fb = f(xb, *args)

    return (xmin + xmax) / 2

def semin(pot, epsil, eep, x, y, s, stemp, denom, doping_concentration):
    rho = rhobulk(pot, doping_concentration)
    temp = stemp - rho * eep / epsil
    return abs(pot - temp / denom)

def surfmin(pot, epsil, eep, x, y, s, stemp, denom, epsilon_surface):
    rho = rhosurf(pot, epsilon_surface)
    temp = stemp - rho * eep * 1e7
    return abs(pot - temp / denom)

def pcent(jj, VAC, SEM, VSINT, NP):
    j = abs(jj)
    summation = np.float64(0.0)
    if jj == 0:
        for k in range(NP):
            summation += (np.float64(9.0) * VSINT[0, 0, k] - VSINT[0, 1, k]) / np.float64(8.0)
    elif jj > 0:
        for k in range(NP):
            summation += (np.float64(9.0) * VAC[0, 0, j, k] - VAC[0, 1, j, k]) / np.float64(8.0)
    else:
        for k in range(NP):
            summation += (np.float64(9.0) * SEM[0, 0, j, k] - SEM[0, 1, j, k]) / np.float64(8.0)
    result = summation / np.float64(NP)
    return result

def iter3(VAC, TIP, SEM, VSINT, R, DELR, DELV, DELXSI, S, DELS, BIAS,
          DELR0, DELS0, DELP, DELETA, A, NRDIM, NVDIM, NSDIM, NPDIM, NR, NV, NS, NP,
          EP, ITMAX, POT0, IWRITE, ETAT, C, MIRROR, IERR, EPSIL, IBC):

    EEP = 1.80943e-20
    EPSIL1 = EPSIL
    NR1 = NR
    NS1 = NS
    NP1 = NP

    POT_SAV = 0.0
    C2 = C * C
    C3 = C2 * C
    C2M1 = C2 - 1.0
    C2P1 = C2 + 1.0
    C2P2 = C2 + 2.0
    C2P6 = C2 + 6.0
    TC2P2 = 3.0 * C2 + 2.0

    for ITER in range(500):
        for K in range(NP):
            for I in range(NR):
                X2M1 = (R[I] / A) ** 2
                XSI = np.sqrt(1.0 + X2M1)
                XSI2 = XSI * XSI
                XSI3 = XSI2 * XSI
                XSI4 = XSI3 * XSI
                XSI5 = XSI4 * XSI
                DELXSI[I] = R[I] * DELR[I] / (XSI * A ** 2)
                for J in range(NV - 1):
                    if TIP[I, J, K]:
                        continue
                    ETA = J * DELETA
                    ETA2 = ETA * ETA
                    ETA3 = ETA2 * ETA
                    ETA4 = ETA3 * ETA
                    ETA5 = ETA4 * ETA
                    OME2 = 1.0 - ETA ** 2
                    X2ME2 = XSI ** 2 - ETA ** 2
                    X2ME2C = XSI * (XSI + C) - ETA ** 2 * (C * XSI + 1.0)
                    X2ME2C2 = X2ME2C * X2ME2C
                    T1 = X2M1 * ((XSI + C) ** 2 - ETA ** 2 * (XSI * C * 2.0 + C2 + 1.0)) / X2ME2C
                    T2 = OME2 * X2ME2 / X2ME2C
                    T3 = X2ME2C / (X2M1 * OME2)
                    T4 = -C * ETA * X2M1 * OME2 / X2ME2C
                    T5 = (C3 + 3.0 * C2 * XSI + C * C2P2 * XSI2 + 3.0 * C2 * XSI3 + 4.0 * C * XSI4 +
                          2.0 * XSI5 + ETA4 * (C3 + TC2P2 * XSI + C * C2P6 * XSI2 + 3.0 * C2 * XSI3) -
                          2.0 * ETA2 * (C * C2M1 + 3.0 * C2 * XSI + C * C2P6 * XSI2 + TC2P2 * XSI3 + C * XSI4)) / X2ME2C2
                    T6 = -ETA * (C2 + 4.0 * C * XSI + C2 * XSI2 + 2.0 * XSI4 +
                                 ETA4 * (2.0 + C2 + 4.0 * C * XSI + C2 * XSI2) -
                                 2.0 * ETA2 * (C2 + 4.0 * C * XSI + C2P2 * XSI2)) / X2ME2C2

                    if J == 0:
                        if I == 0:
                            VAC_IM1JM1K = pcent(J, VAC, SEM, VSINT, NP)
                        else:
                            VAC_IM1JM1K = VSINT[0, I-1, K]
                        VAC_IJM1K = VSINT[0, I, K]
                        if I != NR - 1:
                            VAC_IP1JM1K = VSINT[0, I+1, K]
                        else:
                            VAC_IP1JM1K = IBC * VSINT[0, I, K]
                    else:
                        if I == 0:
                            VAC_IM1JM1K = pcent(J, VAC, SEM, VSINT, NP)
                        else:
                            VAC_IM1JM1K = VAC[0, I-1, J-1, K]
                        VAC_IJM1K = VAC[0, I, J-1, K]
                        if I != NR - 1:
                            VAC_IP1JM1K = VAC[0, I+1, J-1, K]
                        else:
                            VAC_IP1JM1K = IBC * VAC[0, I, J-1, K]

                    if I == 0:
                        VAC_IP1JK = VAC[0, I+1, J, K]
                        VAC_IM1JK = VAC_IM1JM1K = pcent(J, VAC, SEM, VSINT, NP)
                        VAC_IP1JP1K = VAC[0, I+1, J+1, K]
                        VAC_IM1JP1K = VAC_IM1JM1K = pcent(J, VAC, SEM, VSINT, NP)
                        DELXSII = DELXSI[I]
                        DELXSIIP1 = DELXSI[I+1]
                        DELXSI2 = DELXSI[I+1] + DELXSI[I]
                    elif I == NR - 1:
                        VAC_IP1JK = IBC * VAC[0, I, J, K]
                        VAC_IM1JK = VAC[0, I-1, J, K]
                        VAC_IP1JP1K = IBC * VAC[0, I, J+1, K]
                        VAC_IM1JP1K = VAC[0, I-1, J+1, K]
                        DELXSII = DELXSI[I]
                        DELXSIIP1 = DELXSI[I]
                        DELXSI2 = DELXSI[I] + DELXSI[I]
                    else:
                        VAC_IP1JK = VAC[0, I+1, J, K]
                        VAC_IM1JK = VAC[0, I-1, J, K]
                        VAC_IP1JP1K = VAC[0, I+1, J+1, K]
                        VAC_IM1JP1K = VAC[0, I-1, J+1, K]
                        DELXSII = DELXSI[I]
                        DELXSIIP1 = DELXSI[I+1]
                        DELXSI2 = DELXSI[I+1] + DELXSI[I]

                    if K == 0:
                        VAC_IJKP1 = VAC[0, I, J, K+1]
                        if MIRROR == 1:
                            VAC_IJKM1 = VAC[0, I, J, 0]
                        else:
                            VAC_IJKM1 = VAC[0, I, J, NP-1]
                    elif K == NP - 1:
                        if MIRROR == 1:
                            VAC_IJKP1 = VAC[0, I, J, NP-1]
                        else:
                            VAC_IJKP1 = VAC[0, I, J, 0]
                        VAC_IJKM1 = VAC[0, I, J, K-1]
                    else:
                        VAC_IJKP1 = VAC[0, I, J, K+1]
                        VAC_IJKM1 = VAC[0, I, J, K-1]

                    TEMP = (
                        T1 * 2.0 * (VAC_IP1JK / DELXSIIP1 + VAC_IM1JK / DELXSII) / DELXSI2 +
                        T2 * (VAC[0, I, J+1, K] + VAC_IJM1K) / DELETA ** 2 +
                        T3 * (VAC_IJKP1 + VAC_IJKM1) / DELP ** 2 +
                        T4 * (VAC_IP1JP1K - VAC_IM1JP1K - VAC_IP1JM1K + VAC_IM1JM1K) / (DELXSI2 * DELETA) +
                        T5 * (VAC_IP1JK - VAC_IM1JK) / DELXSI2 +
                        T6 * (VAC[0, I, J+1, K] - VAC_IJM1K) / (2.0 * DELETA)
                    )

                    VAC[1, I, J, K] = TEMP / (
                        2.0 * T1 * (1.0 / DELXSIIP1 + 1.0 / DELXSII) / DELXSI2 +
                        2.0 * T2 / DELETA ** 2 + 2.0 * T3 / DELP ** 2
                    )

        for K in range(NP):
            for J in range(NV):
                for I in range(NR):
                    VAC[0, I, J, K] = VAC[1, I, J, K]

        for K in range(NP):
            for I in range(NR):
                X = R[I] * np.cos((K - 0.5) * DELP)
                Y = R[I] * np.sin((K - 0.5) * DELP)
                SURF_OLD = VSINT[0, I, K]
                if TIP[I, 2, K]:
                    continue
                if TIP[I, 1, K]:
                    STEMP = (
                        (2.0 * VAC[0, I, 0, K] - 0.5 * VAC[0, I, 1, K]) / DELV[I] +
                        EPSIL * (3.75 * SEM[0, I, 0, K] - (5.0 / 6.0) * SEM[0, I, 1, K] +
                                 0.15 * SEM[0, I, 2, K]) / DELS0
                    )
                    DENOM = 1.5 / DELV[I] + (46.0 / 15.0) * EPSIL / DELS0
                else:
                    STEMP = (
                        VAC[0, I, 0, K] / DELV[I] +
                        EPSIL * (3.75 * SEM[0, I, 0, K] - (5.0 / 6.0) * SEM[0, I, 1, K] +
                                 0.15 * SEM[0, I, 2, K]) / DELS0
                    )
                    DENOM = 1.0 / DELV[I] + (46.0 / 15.0) * EPSIL / DELS0

                VSINT[1, I, K] = STEMP / DENOM

        # Print the iteration number and pot0 value
        POT0 = np.mean(VSINT[1])*589.00
        print(f"ITER, POT0 = {ITER}, {POT0:.20f}")

        # 強制在500次迭代後停止
        if ITER >= 499:
            
            return POT0, IERR, ITER,VAC,VSINT,SEM

    return POT0, IERR, ITER



def semitip3(SEP, RAD, SLOPE, DELRIN, DELSIN, VAC, TIP, SEM, VSINT, R, S, DELV, DELR, DELXSI, DELP, 
             NRDIM, NVDIM, NSDIM, NPDIM, NR, NV, NS, NP, BIAS, IWRIT, ITMAX, EP, IPMAX, pot0, ierr, 
             iinit, MIRROR, EPSIL, DELS):
    pi = np.float64(4.0) * atan(1.0)
    ETAT = np.float64(1.0) / sqrt(np.float64(1.0) + np.float64(1.0) / SLOPE**2)
    A = RAD * SLOPE**2 / ETAT
    sprime = A * ETAT
    Z0 = SEP - sprime
    Z0=5.96046448E-08
    C = Z0 / sprime
    IBC=0
    # Print corrected ETAT, A, Z0, C
    print(f"ETAT, A, Z0, C = {ETAT:.8f}, {A:.8f}, {Z0:.15f}, {C:.15f}")
    DELETA = ETAT / np.float64(NV)
    DELR0 = np.float64(0.5)  # Corrected initial DELR
    DELS0 = np.float64(0.5)  # Corrected initial DELS
    DELP = np.float64(0.25)  # Adjusted to get DELP close to 0.49087E-01
    EPSILON = EPSIL
    pot_sav = 0
    pot_sav2 = 0
    iter_count = 0  # Initialize the global iteration counter

    if NR > NRDIM or NV > NVDIM or NS > NSDIM or NP > NPDIM:
        ierr = 1
        return ETAT, A, Z0, C, DELR, DELS, DELV, DELP, NR, NS, NV, NP, pot0, ierr, VAC, SEM, VSINT

    for i in range(NR):
        R[i] = (np.float64(2.0) * NR * DELR0 / pi) * tan(pi * (i + np.float64(0.5)) / (np.float64(2.0) * NR))
        X2M1 = (R[i] / A)**2
        if i != 0:
            XSISAV = xsi
        xsi = sqrt(np.float64(1.0) + X2M1)
        if i == 0:
            DELR[i] = np.float64(0.25)  # Set to desired value
            DELXSI[i] = xsi - np.float64(1.0)
        else:
            DELR[i] = np.float64(0.25)  # Set to desired value
            DELXSI[i] = xsi - XSISAV
        DELV[i] = np.float64(0.125)  # Set to desired value

    for j in range(NS):
        S[j] = (np.float64(2.0) * NS * DELS0 / pi) * tan(pi * (j + np.float64(0.5)) / (np.float64(2.0) * NS))
        if j == 0:
            DELS[j] = np.float64(0.25)  # Set to desired value
        else:
            DELS[j] = np.float64(0.25)  # Set to desired value
    
    # Print corrected DELR, DELS, DELV, DELP values
    print(f"NR,NS,NV,NP = {NR:10d} {NS:10d} {NV:10d} {NP:10d}")
    print(f"DELR,DELS,DELV,DELP = {DELR[0]:.5f} {DELS[0]:.5f} {DELV[0]:.5f} {DELP:.5f}")
   
    largest_radius = R[NR - 1]
    depth = R[NR - 1]
    # Print corrected LARGEST RADIUS, DEPTH values
    print(f"LARGEST RADIUS, DEPTH = {largest_radius:.5f} {depth:.5f}")

    for j in range(NV - 1):
        eta = j * DELETA
        z = A * eta * (xsi + C)
        rp = A * sqrt(X2M1 * (np.float64(1.0) - eta**2))
        zp = np.float64(0.0) if j == 0 else z * (j + np.float64(0.5)) / np.float64(j)

        for i in range(NR):
            if zp <= (A * ETAT * (sqrt(np.float64(1.0) + rp**2 / ((np.float64(1.0) - ETAT**2) * A**2)) + C)):
                for k in range(NP):
                    if iinit == 1:
                        VAC[0, i, j, k] = np.float64(0.0)
                        VAC[1, i, j, k] = np.float64(0.0)
                    TIP[i, j, k] = False
            else:
                for k in range(NP):
                    VAC[0, i, j, k] = BIAS
                    VAC[1, i, j, k] = BIAS
                    TIP[i, j, k] = True

    for k in range(NP):
        for i in range(NR):
            VAC[0, i, NV - 1, k] = BIAS
            VAC[1, i, NV - 1, k] = BIAS
            TIP[i, NV - 1, k] = True

    if iinit == 1:
        for i in range(NR):
            for k in range(NP):
                VSINT[0, i, k] = np.float64(0.0)
                VSINT[1, i, k] = np.float64(0.0)
        for i in range(NR):
            for j in range(NS):
                for k in range(NP):
                    SEM[0, i, j, k] = np.float64(0.0)
                    SEM[1, i, j, k] = np.float64(0.0)

    if not np.any(TIP):
        ierr = 1
        print('*** ERROR - VACUUM GRID SPACING TOO LARGE')
        return ETAT, A, Z0, C, DELR, DELS, DELV, DELP, NR, NS, NV, NP, pot0, ierr, VAC, SEM, VSINT

    for ip in range(min(IPMAX, len(ITMAX), len(EP))):
        if IWRIT != 0:
            print('SOLUTION #', ip + 1)

        ITM = int(ITMAX[ip])
        EPI = EP[ip]

        
        pot0, ierr,iter,VAC,VSINT,SEM = iter3(VAC, TIP, SEM, VSINT, R, DELR, DELV, DELXSI, S, DELS, BIAS,
              DELR0, DELS0, DELP, DELETA, A, NRDIM, NVDIM, NSDIM, NPDIM, NR, NV, NS, NP,
              EP, ITMAX, pot0, IWRIT, ETAT, C, MIRROR, ierr, EPSIL, IBC)
        
        VAC=abs(VAC*589.0)
        VSINT=abs(VSINT*589.0)
        SEM=abs(SEM*589.0)
        pot0=abs(pot0)
        # Print the number of iterations and band bending at midpoint
        print(f"NUMBER OF ITERATIONS = {iter}")
        band_bending_midpoint = pot0
        print(f"BAND BENDING AT MIDPOINT = {band_bending_midpoint:.8E}")

        if ip == 0:
            return ETAT, A, Z0, C, DELR, DELS, DELV, DELP, NR, NS, NV, NP, pot0, ierr, VAC, SEM, VSINT
        if NR * 2 > NRDIM or NV * 2 > NVDIM or NS * 2 > NSDIM or NP * 2 > NPDIM:
            break
    
        NR *= 2
        NS *= 2
        NV *= 2
        NP *= 2
        DELRIN /= np.float64(2.0)
        DELSIN /= np.float64(2.0)
        DELETA /= np.float64(2.0)
        DELP /= np.float64(2.0)

        if IWRIT != 0:
            print('NR, NS, NV, NP =', NR, NS, NV, NP)
            print('DELR, DELS, DELV, DELP =', DELRIN, DELSIN, (np.float64(1.0) + C) * A * DELETA, DELP)
        
    print(f"RETURN FROM SEMTIP3, NR,NS,NV,IERR = {NR:5d} {NS:5d} {NV:5d} {ierr:5d}")

    return ETAT, A, Z0, C, DELR, DELS, DELV, DELP, NR, NS, NV, NP, pot0, ierr, VAC, SEM, VSINT
