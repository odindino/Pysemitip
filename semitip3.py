import numpy as np
from math import atan, sqrt, log, pi, tan, cos, sin
from gsect import*

def rhobulk(pot, x, y, s, i, j, k, NR, NS, NP):
    epsilon_bulk = 11.7 * 8.85e-12
    q = 1.6e-19
    kT = 0.0259
    ni = 1.5e10
    
    doping_concentration = 1e17
    
    if pot > 0:
        rho = q * doping_concentration * (1 - np.exp(-q * pot / (kT)))
    elif pot < 0:
        rho = q * doping_concentration * (np.exp(q * pot / (kT)) - 1)
    else:
        rho = 0
    
    return rho

def rhosurf(pot, x, y, i, k, NR, NP):
    epsilon_surface = 11.7 * 8.85e-12
    q = 1.6e-19
    kT = 0.0259
    ni = 1.5e10
    
    if pot > 0:
        rho = epsilon_surface * pot / (q * ni * kT)
    elif pot < 0:
        rho = -epsilon_surface * abs(pot) / (q * ni * kT)
    else:
        rho = 0
    
    return rho

def semin(pot, epsil, eep, x, y, s, stemp, denom, i, j, k, NR, NS, NP):
    rho = rhobulk(pot, x, y, s, i, j, k, NR, NS, NP)
    temp = stemp - rho * eep / epsil
    return abs(pot - temp / denom)

def surfmin(pot, epsil, eep, x, y, s, stemp, denom, i, k, NR, NP):
    rho = rhosurf(pot, x, y, i, k, NR, NP)
    temp = stemp - rho * eep * 1e7
    return abs(pot - temp / denom)

def pcent(jj, VAC, SEM, VSINT, NRDIM, NVDIM, NSDIM, NPDIM, NP):
    j = abs(jj)
    i = 0
    summation = 0.0
    if jj == 0:
        for k in range(NP):
            summation += (9.0 * VSINT[0, i, k] - VSINT[0, i + 1, k]) / 8.0
    elif jj > 0:
        for k in range(NP):
            summation += (9.0 * VAC[0, i, j, k] - VAC[0, i + 1, j, k]) / 8.0
    else:
        for k in range(NP):
            summation += (9.0 * SEM[0, i, j, k] - SEM[0, i + 1, j, k]) / 8.0
    return summation / float(NP)

def iter3(VAC, TIP, SEM, VSINT, R, DELR, DELV, DELXSI, S, DELS, BIAS, DELR0, DELS0, DELP, DELETA, A, NRDIM, NVDIM, NSDIM, NPDIM, NR, NV, NS, NP, EP, ITMAX, pot0, IWRIT, ETAT, C, MIRROR, ierr, EPSILON, IBC):
    eep = 1.80943e-20
    c2 = C * C
    c3 = c2 * C
    c2m1 = c2 - 1.0
    c2p1 = c2 + 1.0
    c2p2 = c2 + 2.0
    c2p6 = c2 + 6.0
    tc2p2 = 3.0 * c2 + 2.0
    pot_sav = pot0
    pot_sav2 = pot0
    

    damping_factor = 1  # 阻尼系数
    step_size = 1  # 步长大小

    for iter in range(ITMAX):
        if iter != 0 :
            pot_sav2 = pot_sav
            pot_sav = pot0
            pot0 = pcent(0, VAC, SEM, VSINT, NRDIM, NVDIM, NSDIM, NPDIM, NP)
            print(f"ITER,Pot0 = {iter} {pot0:.9f}")
        for k in range(NP):
            for i in range(NR):
                x2m1 = (R[i] / A) ** 2
                xsi = sqrt(1.0 + x2m1)
                xsi2 = xsi * xsi
                xsi3 = xsi2 * xsi
                xsi4 = xsi3 * xsi
                xsi5 = xsi4 * xsi
                delxsi_val = R[i] * DELR[i] / (xsi * A ** 2)
                for j in range(NV - 1):
                    if TIP[i, j, k]:
                        continue
                    eta = j * DELETA
                    eta2 = eta * eta
                    eta3 = eta2 * eta
                    eta4 = eta3 * eta
                    eta5 = eta4 * eta
                    ome2 = 1.0 - eta ** 2
                    x2me2 = xsi ** 2 - eta ** 2
                    x2me2c = xsi * (xsi + C) - eta ** 2 * (C * xsi + 1.0)
                    x2me2c2 = x2me2c * x2me2c
                    t1 = x2m1 * ((xsi + C) ** 2 - eta ** 2 * (xsi * C * 2.0 + c2 + 1.0)) / x2me2c
                    t2 = ome2 * x2me2 / x2me2c
                    t3 = x2me2c / (x2m1 * ome2)
                    t4 = -C * eta * x2m1 * ome2 / x2me2c
                    t5 = (c3 + 3.0 * c2 * xsi + C * c2p2 * xsi2 + 3.0 * c2 * xsi3 + 4.0 * xsi4 + 2.0 * xsi5 +
                          eta4 * (c3 + tc2p2 * xsi + C * c2p6 * xsi2 + 3.0 * c2 * xsi3) -
                          2.0 * eta2 * (C * c2m1 + 3.0 * c2 * xsi + C * c2p6 * xsi2 + tc2p2 * xsi3 + C * xsi4)) / x2me2c2
                    t6 = -eta * (c2 + 4.0 * C * xsi + c2 * xsi2 + 2.0 * xsi4 +
                                 eta4 * (2.0 + c2 + 4.0 * C * xsi + c2 * xsi2) -
                                 2.0 * eta2 * (c2 + 4.0 * C * xsi + c2p2 * xsi2)) / x2me2c2

                    if j == 0:
                        if i == 0:
                            vac_im1jm1k = pcent(j, VAC, SEM, VSINT, NRDIM, NVDIM, NSDIM, NPDIM, NP)
                        else:
                            vac_im1jm1k = VSINT[0, i - 1, k]
                        vac_ijm1k = VSINT[0, i, k]
                        if i != NR - 1:
                            vac_ip1jm1k = VSINT[0, i + 1, k]
                        else:
                            vac_ip1jm1k = IBC * VSINT[0, i, k]
                    else:
                        if i == 0:
                            vac_im1jm1k = pcent(j - 1, VAC, SEM, VSINT, NRDIM, NVDIM, NSDIM, NPDIM, NP)
                        else:
                            vac_im1jm1k = VAC[0, i - 1, j - 1, k]
                        vac_ijm1k = VAC[0, i, j - 1, k]
                        if i != NR - 1:
                            vac_ip1jm1k = VAC[0, i + 1, j - 1, k]
                        else:
                            vac_ip1jm1k = IBC * VAC[0, i, j - 1, k]

                    if i == 0:
                        vac_ip1jk = VAC[0, i + 1, j, k]
                        vac_im1jk = pcent(j, VAC, SEM, VSINT, NRDIM, NVDIM, NSDIM, NPDIM, NP)
                        vac_ip1jp1k = VAC[0, i + 1, j + 1, k]
                        vac_im1jp1k = pcent(j + 1, VAC, SEM, VSINT, NRDIM, NVDIM, NSDIM, NPDIM, NP)
                        delxsi_i = DELXSI[i]
                        delxsi_ip1 = DELXSI[i + 1]
                        delxsi2 = DELXSI[i + 1] + DELXSI[i]
                    elif i == NR - 1:
                        vac_ip1jk = IBC * VAC[0, i, j, k]
                        vac_im1jk = VAC[0, i - 1, j, k]
                        vac_ip1jp1k = IBC * VAC[0, i, j + 1, k]
                        vac_im1jp1k = VAC[0, i - 1, j + 1, k]
                        delxsi_i = DELXSI[i]
                        delxsi_ip1 = DELXSI[i]
                        delxsi2 = DELXSI[i] + DELXSI[i]
                    else:
                        vac_ip1jk = VAC[0, i + 1, j, k]
                        vac_im1jk = VAC[0, i - 1, j, k]
                        vac_ip1jp1k = VAC[0, i + 1, j + 1, k]
                        vac_im1jp1k = VAC[0, i - 1, j + 1, k]
                        delxsi_i = DELXSI[i]
                        delxsi_ip1 = DELXSI[i + 1]
                        delxsi2 = DELXSI[i + 1] + DELXSI[i]

                    if k == 0:
                        vac_ijkp1 = VAC[0, i, j, k + 1]
                        if MIRROR == 1:
                            vac_ijkm1 = VAC[0, i, j, 0]
                        else:
                            vac_ijkm1 = VAC[0, i, j, NP - 1]
                    elif k == NP - 1:
                        if MIRROR == 1:
                            vac_ijkp1 = VAC[0, i, j, NP - 1]
                        else:
                            vac_ijkp1 = VAC[0, i, j, 0]
                        vac_ijkm1 = VAC[0, i, j, k - 1]
                    else:
                        vac_ijkp1 = VAC[0, i, j, k + 1]
                        vac_ijkm1 = VAC[0, i, j, k - 1]

                    temp = (t1 * 2.0 * (vac_ip1jk / delxsi_ip1 + vac_im1jk / delxsi_i) / delxsi2 +
                            t2 * (VAC[0, i, j + 1, k] + vac_ijm1k) / DELETA ** 2 +
                            t3 * (vac_ijkp1 + vac_ijkm1) / DELP ** 2 +
                            t4 * (vac_ip1jp1k - vac_im1jp1k - vac_ip1jm1k + vac_im1jm1k) / (delxsi2 * DELETA) +
                            t5 * (vac_ip1jk - vac_im1jk) / delxsi2 +
                            t6 * (VAC[0, i, j + 1, k] - vac_ijm1k) / (2.0 * DELETA))

                    VAC[1, i, j, k] = temp / (2.0 * t1 * (1 / delxsi_ip1 + 1.0 / delxsi_i) / delxsi2 +
                                              2.0 * t2 / DELETA ** 2 + 2.0 * t3 / DELP ** 2)

        for k in range(NP):

            for j in range(NV):
                for i in range(NR):
                    VAC[0, i, j, k] = damping_factor * VAC[1, i, j, k] + (1 - damping_factor) * VAC[0, i, j, k]
                    
        for k in range(NP):
            for i in range(NR):
                x = R[i] * cos((k - 0.5) * DELP)
                y = R[i] * sin((k - 0.5) * DELP)
                surf_old = VSINT[0, i, k]
                if TIP[i, 3, k]:
                    continue
                stemp = ((3.0 * VAC[0, i, 0, k] - (9.0 / 6.0) * VAC[0, i, 1, k] + (1.0 / 3.0) * VAC[0, i, 2, k]) / DELV[i] +
                         EPSILON * (3.75 * SEM[0, i, 0, k] - (5.0 / 6.0) * SEM[0, i, 1, k] + 0.15 * SEM[0, i, 2, k]) / DELS0)
                denom = ((11.0 / 6.0) / DELV[i] + (46.0 / 15.0) * EPSILON / DELS0)

                rho = rhosurf(VSINT[0, i, k], x, y, i, k, NR, NP)
                temp = stemp - rho * eep * 1e7
                surf_new = temp / denom
                del_surf = max(1e-6, abs(BIAS) / 1e6)

                surf_new = gsect(surfmin, surf_old, surf_new, del_surf, EPSILON, eep, x, y, S[i], stemp, denom, i, k, NR, NP)

                VSINT[1, i, k] = damping_factor * surf_new + (1 - damping_factor) * surf_old

        for k in range(NP):
            
            for i in range(NR):
                VSINT[0, i, k] = VSINT[1, i, k]
            
        

        for k in range(NP):
            for j in range(NS):
                for i in range(NR):
                    sem_old = SEM[0, i, j, k]
                    ener = 0.0 - SEM[0, i, j, k]
                    rsav = R[i]
                    x = R[i] * cos((k - 0.5) * DELP)
                    y = R[i] * sin((k - 0.5) * DELP)
                    if i == 0:
                        sem_im1jk = pcent(-j, VAC, SEM, VSINT, NRDIM, NVDIM, NSDIM, NPDIM, NP)
                        sem_ip1jk = SEM[0, i + 1, j, k]
                        delr2 = DELR[i + 1] + DELR[i]
                        delr_ip1 = DELR[i + 1]
                        delr_i = DELR[i]
                    elif i == NR - 1:
                        sem_im1jk = SEM[0, i - 1, j, k]
                        sem_ip1jk = IBC * SEM[0, i, j, k]
                        delr2 = DELR[i] + DELR[i]
                        delr_ip1 = DELR[i]
                        delr_i = DELR[i]
                    else:
                        sem_im1jk = SEM[0, i - 1, j, k]
                        sem_ip1jk = SEM[0, i + 1, j, k]
                        delr2 = DELR[i + 1] + DELR[i]
                        delr_ip1 = DELR[i + 1]
                        delr_i = DELR[i]

                    if j == 0:
                        sem_ijp1k = SEM[0, i, j + 1, k]
                        sem_ijm1k = VSINT[0, i, k]
                        dels2 = DELS[j + 1] + DELS[j]
                        dels_jp1 = DELS[j + 1]
                        dels_j = DELS[j]
                    elif j == NS - 1:
                        sem_ijp1k = IBC * SEM[0, i, j, k]
                        sem_ijm1k = SEM[0, i, j - 1, k]
                        dels2 = DELS[j] + DELS[j]
                        dels_jp1 = DELS[j]
                        dels_j = DELS[j]
                    else:
                        sem_ijp1k = SEM[0, i, j + 1, k]
                        sem_ijm1k = SEM[0, i, j - 1, k]
                        dels2 = DELS[j + 1] + DELS[j]
                        dels_jp1 = DELS[j + 1]
                        dels_j = DELS[j]

                    if k == 0:
                        sem_ijkp1 = SEM[0, i, j, k + 1]
                        if MIRROR == 1:
                            sem_ijkm1 = SEM[0, i, j, 0]
                        else:
                            sem_ijkm1 = SEM[0, i, j, NP - 1]
                    elif k == NP - 1:
                        if MIRROR == 1:
                            sem_ijkp1 = SEM[0, i, j, NP - 1]
                        else:
                            sem_ijkp1 = SEM[0, i, j, 0]
                        sem_ijkm1 = SEM[0, i, j, k - 1]
                    else:
                        sem_ijkp1 = SEM[0, i, j, k + 1]
                        sem_ijkm1 = SEM[0, i, j, k - 1]

                    stemp = (2.0 * (sem_ip1jk / delr_ip1 + sem_im1jk / delr_i) / delr2 +
                             2.0 * (sem_ijp1k / dels_jp1 + sem_ijm1k / dels_j) / dels2 +
                             (sem_ip1jk - sem_im1jk) / (R[i] * delr2) +
                             (sem_ijkp1 + sem_ijkm1) / (R[i] ** 2 * DELP ** 2))

                    rho = rhobulk(SEM[0, i, j, k], x, y, S[j], i, j, k, NR, NS, NP)
                    temp = stemp - rho * eep / EPSILON
                    denom = (2.0 * (1.0 / delr_ip1 + 1.0 / delr_i) / delr2 +
                             2.0 * (1.0 / dels_jp1 + 1.0 / dels_j) / dels2 +
                             2.0 / (R[i] ** 2 * DELP ** 2))
                    sem_new = temp / denom
                    del_sem = max(1e-6, abs(BIAS) / 1e6)

                    sem_new = gsect(semin, sem_old, sem_new, del_sem, EPSILON, eep, x, y, S[j], stemp, denom, i, j, k, NR, NS, NP)

                    SEM[1, i, j, k] = damping_factor * sem_new + (1 - damping_factor) * sem_old

        for k in range(NP):
            for j in range(NS):
                for i in range(NR):
                    SEM[0, i, j, k] = SEM[1, i, j, k]
            
        if iter % 100 == 0 and iter != 0 and abs(pot0 - pot_sav) < EP and abs(pot_sav - pot_sav2) < 2.0 * EP:
            break

    
    pot0 = round(pot0, 9)  # Ensure pot0 precision to 9 decimal places
    return pot0, ierr
# Main semitip3 function
def semitip3(SEP, RAD, SLOPE, DELRIN, DELSIN, VAC, TIP, SEM, VSINT, R, S, DELV, DELR, DELXSI, DELP, 
             NRDIM, NVDIM, NSDIM, NPDIM, NR, NV, NS, NP, BIAS, IWRIT, ITMAX, EP, IPMAX, pot0, ierr, 
             iinit, MIRROR, EPSIL):
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
    #NR, NS, NV, NP = 32 32 16 16
    # Initialize R and DELR arrays
    
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
    
    # Initialize S and DELS arrays
    for j in range(NS):
        S[j] = (2 * NS * DELS0 / pi) *tan(pi * (j + 0.5) / (2.0 * NS))
        if j == 0:
            DELS[j] = S[j]
        else:
            DELS[j] = S[j] - S[j - 1]
    print(DELS)
    # Initialize VAC and TIP arrays
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
    

    pot0 = pcent(0, VAC, SEM, VSINT, NRDIM, NVDIM, NSDIM, NPDIM, NP)
    
    for ip in range(min(IPMAX, len(ITMAX), len(EP))):
        if IWRIT != 0:
            print('SOLUTION #', ip + 1)

        ITM = int(ITMAX[ip])
        EPI = EP[ip]

        for iter in range(0,401,100):
            if iter != 0:
                pot0, ierr = iter3(VAC, TIP, SEM, VSINT, R, DELR, DELV, DELXSI, S, DELS, BIAS, DELRIN, DELSIN, DELP, DELETA, A, NRDIM, NVDIM, NSDIM, NPDIM, NR, NV, NS, NP, EPI, 100, 0.5, IWRIT, ETAT, C, MIRROR, ierr, EPSILON, IBC)
                if iter != 0 and iter % 100 == 0:
                    print(f"ITER,Pot0 = {iter} {pot0:.9f}")

        if ip == IPMAX - 1:
            return ETAT, A, Z0, C, DELR, DELS, DELV, DELP, NR, NS, NV, NP, pot0, ierr

        if NR * 2 > NRDIM or NV * 2 > NVDIM or NS * 2 > NSDIM or NP * 2 > NPDIM:
            break

        NR *= 2
        NS *= 2
        NV *= 2
        NP *= 2
        DELRIN /= 2.0
        DELSIN /= 2.0
        DELETA /= 2.0
        DELP /= 2.0

        if IWRIT != 0:
            print('NR, NS, NV, NP =', NR, NS, NV, NP)
            print('DELR, DELS, DELV, DELP =', DELRIN, DELSIN, (1.0 + C) * A * DELETA, DELP)

    return ETAT, A, Z0, C, DELR, DELS, DELV, DELP, NR, NS, NV, NP, pot0, ierr

def gsect(f, xmin, xmax, ep, *args):
    GS = 0.3819660
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
                result = (xmin + xmax) / 2
                return result
            xb = xa
            fb = fa
            xa = xmin + delx * GS
            fa = f(xa, *args)
        else:
            xmin = xa
            delx = xmax - xmin
            if delx == delxsav:
                result = (xmin + xmax) / 2
                return result
            xa = xb
            fa = fb
            xb = xmax - delx * GS
            fb = f(xb, *args)

    result = (xmin + xmax) / 2
    return result
