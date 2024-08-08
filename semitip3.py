import numpy as np
from math import atan, sqrt, log, pi, tan, cos, sin

def rhobulk(pot, doping_concentration=1e18, epsilon_bulk=11.7 * 8.85e-12, q=1.6e-19, kT=0.0259):
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
             iinit, MIRROR, EPSIL, DELS):
    pi = np.float64(4.0) * atan(1.0)
    ETAT = np.float64(1.0) / sqrt(np.float64(1.0) + np.float64(1.0) / SLOPE**2)
    A = RAD * SLOPE**2 / ETAT
    sprime = A * ETAT
    Z0 = 100
    C = 100
    DELETA = ETAT / np.float64(NV)
    DELR0 = np.float64(10.0)
    DELS0 = np.float64(10.0)
    EPSILON = EPSIL
    pot_sav = 0
    pot_sav2 = 0
    iter_count = 0  # Initialize the global iteration counter
    
    if NR > NRDIM or NV > NVDIM or NS > NSDIM or NP > NPDIM:
        ierr = 1
        return ETAT, A, Z0, C, DELR, DELS, DELV, DELP, NR, NS, NV, NP, pot0, ierr

    for i in range(NR):
        R[i] = (np.float64(2.0) * NR * DELR0 / pi) * tan(pi * (i + np.float64(0.5)) / (np.float64(2.0) * NR))
        X2M1 = (R[i] / A)**2
        if i != 0:
            XSISAV = xsi
        xsi = sqrt(np.float64(1.0) + X2M1)
        if i == 0:
            DELR[i] = R[i]
            DELXSI[i] = xsi - np.float64(1.0)
        else:
            DELR[i] = R[i] - R[i - 1]
            DELXSI[i] = xsi - XSISAV
        DELV[i] = (sqrt(A ** 2 + R[i] ** 2) + C * A) * ETAT / np.float64(NV)

    for j in range(NS):
        S[j] = (np.float64(2.0) * NS * DELS0 / pi) * tan(pi * (j + np.float64(0.5)) / (np.float64(2.0) * NS))
        if j == 0:
            DELS[j] = S[j]
        else:
            DELS[j] = S[j] - S[j - 1]

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
        if iinit == 1:
            for i in range(NR):
                for k in range(NP):
                    SEM[0, i, j, k] = np.float64(0.0)
                    SEM[1, i, j, k] = np.float64(0.0)

    if not np.any(TIP):
        ierr = 1
        print('*** ERROR - VACUUM GRID SPACING TOO LARGE')
        return ETAT, A, Z0, C, DELR, DELS, DELV, DELP, NR, NS, NV, NP, pot0, ierr

    
    for ip in range(min(IPMAX, len(ITMAX), len(EP))):
        if IWRIT != 0:
            print('SOLUTION #', ip + 1)

        ITM = int(ITMAX[ip])
        EPI = EP[ip]

        for iter in range(50):
            pot0, ierr, iter_count = iter(VAC, TIP, SEM, VSINT, R, DELR, DELV, DELP, DELXSI, S, DELS, BIAS, A, NR, NV, NS, NP, EPI, ITM, pot0, IWRIT, ETAT, C, MIRROR, EPSILON, DELETA, iter_count)
            if iter % 100 == 0:
                print(f"ITER, Pot0 = {iter_count}, {pot0:.9f}")

            # Check for convergence every 100 iterations
            if iter % 100 == 0 and abs(pot0 - pot_sav) < EPI and abs(pot_sav - pot_sav2) < 2 * EPI:
                break

            if iter_count >= 20:  # 強制在20次迭代後停止
                print(f"FORCED STOP AFTER {iter_count} ITERATIONS")
                return ETAT, A, Z0, C, DELR, DELS, DELV, DELP, NR, NS, NV, NP, pot0, ierr,VAC, SEM,VSINT
        
        # Return after first solution
        if ip == 0:
            return ETAT, A, Z0, C, DELR, DELS, DELV, DELP, NR, NS, NV, NP, pot0, ierr,VAC, SEM,VSINT

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

    return ETAT, A, Z0, C, DELR, DELS, DELV, DELP, NR, NS, NV, NP, pot0, ierr, VAC , SEM,  VSINT

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
