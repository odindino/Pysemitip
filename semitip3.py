import numpy as np
from math import atan, sqrt, log, pi, tan, cos, sin

def rhobulk(pot, doping_concentration=1e17, epsilon_bulk=11.7 * 8.85e-12, q=1.6e-19, kT=0.0259):
    if pot > 0:
        rho = q * doping_concentration * (1 - np.exp(-q * pot / (kT)))
    elif pot < 0:
        rho = q * doping_concentration * (np.exp(q * pot / (kT)) - 1)
    else:
        rho = 0
    return rho

def rhosurf(pot, epsilon_surface=11.7 * 8.85e-12, q=1.6e-19, kT=0.0259, ni=1.5e10):
    if pot > 0:
        rho = epsilon_surface * pot / (q * ni * kT)
    elif pot < 0:
        rho = -epsilon_surface * abs(pot) / (q * ni * kT)
    else:
        rho = 0
    return rho

def direct_calculate_potentials(VAC, TIP, SEM, VSINT, R, DELR, DELV, DELXSI, S, DELS, BIAS, DELR0, DELS0, DELP, DELETA, A, NR, NV, NS, NP, EPSILON, eep=1.80943e-20):
    for k in range(NP):
        for j in range(NS):
            for i in range(NR):
                x = R[i] * cos((k - 0.5) * DELP)
                y = R[i] * sin((k - 0.5) * DELP)
                if TIP[i, j, k]:
                    continue
                rho_bulk = rhobulk(SEM[0, i, j, k])
                rho_surf = rhosurf(VSINT[0, i, k])
                stemp = (rho_bulk + rho_surf) * eep / EPSILON
                VAC[0, i, j, k] = stemp

    for k in range(NP):
        for i in range(NR):
            x = R[i] * cos((k - 0.5) * DELP)
            y = R[i] * sin((k - 0.5) * DELP)
            rho_surf = rhosurf(VSINT[0, i, k])
            stemp = rho_surf * eep * 1e7
            VSINT[0, i, k] = stemp

    for k in range(NP):
        for j in range(NS):
            for i in range(NR):
                x = R[i] * cos((k - 0.5) * DELP)
                y = R[i] * sin((k - 0.5) * DELP)
                rho_bulk = rhobulk(SEM[0, i, j, k])
                stemp = rho_bulk * eep / EPSILON
                SEM[0, i, j, k] = stemp

def calculate_pot0(VAC):
    # Calculate the average potential as pot0
    return np.mean(VAC[0])

def semitip3(SEP, RAD, SLOPE, DELRIN, DELSIN, VAC, TIP, SEM, VSINT, R, S, DELV, DELR, DELXSI, DELP, 
             NRDIM, NVDIM, NSDIM, NPDIM, NR, NV, NS, NP, BIAS, IWRIT, ITMAX, EP, IPMAX, pot0, ierr, 
             iinit, MIRROR, EPSIL):
    print(BIAS)
    pi = 4.0 * atan(1.0)
    IBC = 0  # Boundary condition (0=Dirichlet, 1=von Neumann)
    ierr = 0
    IINIT = iinit
    ETAT = 1.0 / sqrt(1.0 + 1.0 / SLOPE**2)
    A = RAD * SLOPE**2 / ETAT
    sprime = A * ETAT
    Z0 = SEP - sprime
    C = Z0 / sprime
    DELETA = ETAT / float(NV)
    DELR0 = 10
    DELS0 = 10
    DELP = 0.19635
    DELR = np.zeros(NR)
    DELS = np.zeros(NS)
    EPSILON = 11.7 * 8.85e-12 
    
    if IWRIT != 0:
        print('ETAT, A, Z0, C =', ETAT, A, Z0, C)
        print('NR, NS, NV, NP =', NR, NS, NV, NP)
        print('DELR, DELS, DELV, DELP =', DELR0, DELS0, 0.03125, DELP)

    if NR > NRDIM or NV > NVDIM or NS > NSDIM or NP > NPDIM:
        ierr = 1
        return ETAT, A, Z0, C, DELR, DELS, DELV, DELP, NR, NS, NV, NP, pot0, ierr

    for i in range(NR):
        R[i] = (2 * NR * DELR0 / pi) * tan(pi * (i + 0.5) / (2.0 * NR))
        X2M1 = (R[i] / A)**2
        if i != 0:
            XSISAV = xsi
        xsi = sqrt(1.0 + X2M1)
        if i == 0:
            DELR[i] = R[i]
            DELXSI[i] = xsi - 1.0
        else:
            DELR[i] = R[i] - R[i - 1]
            DELXSI[i] = xsi - XSISAV
        DELV[i] = (sqrt(A ** 2 + R[i] ** 2) + C * A) * ETAT / float(NV) 
    
    for j in range(NS):
        S[j] = (2 * NS * DELS0 / pi) *tan(pi * (j + 0.5) / (2.0 * NS))
        if j == 0:
            DELS[j] = S[j]
        else:
            DELS[j] = S[j] - S[j - 1]

    for j in range(NV - 1):
        eta = j * DELETA
        z = A * eta * (xsi + C)
        rp = A * sqrt(X2M1 * (1.0 - eta**2))
        if j == 0:
            zp = 0.0
        else:
            zp = z * (j + 0.5) / float(j)
        
        for i in range(NR):
            if zp <= (A * ETAT * (sqrt(1.0 + rp**2 / ((1.0 - ETAT**2) * A**2)) + C)):
                for k in range(NP):
                    if IINIT == 1:
                        VAC[0, i, j, k] = 0.0
                        VAC[1, i, j, k] = 0.0
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

    if IINIT == 1:
        for i in range(NR):
            for k in range(NP):
                VSINT[0, i, k] = 0.0
                VSINT[1, i, k] = 0.0
        if IINIT == 1:
            for i in range(NR):
                for k in range(NP):
                    SEM[0, i, j, k] = 0.0
                    SEM[1, i, j, k] = 0.0
    
    if not np.any(TIP):
        ierr = 1
        print('*** ERROR - VACUUM GRID SPACING TOO LARGE')
        return ETAT, A, Z0, C, DELR, DELS, DELV, DELP, NR, NS, NV, NP, pot0, ierr

    if IWRIT != 0:
        print('LARGEST RADIUS, DEPTH =', R[NR - 1], S[NS - 1])

    if ierr == 1:
        return ETAT, A, Z0, C, DELR, DELS, DELV, DELP, NR, NS, NV, NP, pot0, ierr

    direct_calculate_potentials(VAC, TIP, SEM, VSINT, R, DELR, DELV, DELXSI, S, DELS, BIAS, DELR0, DELS0, DELP, DELETA, A, NR, NV, NS, NP, EPSILON)

    pot0 = calculate_pot0(VAC)  # Calculate pot0 based on the potentials
    pot0 = round(pot0, 9)  # Ensure pot0 precision to 9 decimal places
    return ETAT, A, Z0, C, DELR, DELS, DELV, DELP, NR, NS, NV, NP, pot0, ierr


