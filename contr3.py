import numpy as np

def contr3(ETA1, VAC, TIP, SEM, VSINT, R, S, DELV, NRDIM, NVDIM, NSDIM, NPDIM, NR, NV, NS, NP, NUMC, DELPOT, MIRROR, KPLOT1, KPLOT2):
    
    NRS = NR // 500
    if NRS == 0:
        NRS = 1

    # Draw TIP
    for i in range(-NR + 1, NR + 1, NRS):
        if i == 0:
            continue
        if i > 0:
            II = i
            RSAV = R[II - 1]
        else:
            II = -i
            RSAV = -R[II - 1]
        
        for j in range(1, NV + 1):
            if not TIP[II - 1, j - 1, 0]:
                continue
            print(RSAV * np.sqrt(1 - (ETA1 * j / NV) ** 2), -j * DELV[II - 1])
            break

    # Search for min, max points in potential
    PMIN = 1e10
    PMAX = -1e10
    for i in range(1, NR + 1):
        for k in range(1, NP + 1, NP - 1):
            for j in range(1, NV + 1):
                if PMIN > VAC[0, i - 1, j - 1, k - 1]:
                    PMIN = VAC[0, i - 1, j - 1, k - 1]
                if PMAX < VAC[0, i - 1, j - 1, k - 1]:
                    PMAX = VAC[0, i - 1, j - 1, k - 1]
            if PMIN > VSINT[0, i - 1, k - 1]:
                PMIN = VSINT[0, i - 1, k - 1]
            if PMAX < VSINT[0, i - 1, k - 1]:
                PMAX = VSINT[0, i - 1, k - 1]
            for j in range(1, NS + 1):
                if PMIN > SEM[0, i - 1, j - 1, k - 1]:
                    PMIN = SEM[0, i - 1, j - 1, k - 1]
                if PMAX < SEM[0, i - 1, j - 1, k - 1]:
                    PMAX = SEM[0, i - 1, j - 1, k - 1]

    print('MIN, MAX POTENTIAL VALUES =', PMIN, PMAX)

    # Draw contours
    if DELPOT == 0:
        DELPOT = (PMAX - PMIN) / (NUMC + 1)
        print('CONTOUR SPACING =', DELPOT)

    for i in range(-NR + 1, NR + 1, NRS):
        if i == 0:
            continue
        if i > 0:
            II = i
            RSAV = R[II - 1]
            KP = KPLOT1
        else:
            II = -i
            RSAV = -R[II - 1]
            KP = KPLOT2

        KDONE = np.zeros(NUMC, dtype=bool)

        for k in range(1, NUMC + 1):
            for j in range(NS, 0, -1):
                P = k * DELPOT + PMIN
                if j == 1:
                    if (SEM[0, II - 1, j - 1, KP - 1] >= P and VSINT[0, II - 1, KP - 1] <= P) or (SEM[0, II - 1, j - 1, KP - 1] <= P and VSINT[0, II - 1, KP - 1] >= P):
                        print(RSAV, S[j - 1])
                        KDONE[k - 1] = True
                        break
                else:
                    if (SEM[0, II - 1, j - 1, KP - 1] >= P and SEM[0, II - 1, j - 2, KP - 1] <= P) or (SEM[0, II - 1, j - 1, KP - 1] <= P and SEM[0, II - 1, j - 2, KP - 1] >= P):
                        print(RSAV, S[j - 1])
                        KDONE[k - 1] = True
                        break

            for k in range(1, NUMC + 1):
                P = k * DELPOT + PMIN
                if (VSINT[0, II - 1, KP - 1] >= P and VAC[0, II - 1, 0, KP - 1] <= P) or (VSINT[0, II - 1, KP - 1] <= P and VAC[0, II - 1, 0, KP - 1] >= P):
                    if not KDONE[k - 1]:
                        print(RSAV, 0)
                        KDONE[k - 1] = True

            for k in range(1, NUMC + 1):
                for j in range(1, NV):
                    if TIP[II - 1, j - 1, KP - 1]:
                        continue
                    P = k * DELPOT + PMIN
                    if (VAC[0, II - 1, j - 1, KP - 1] >= P and VAC[0, II - 1, j, KP - 1] <= P) or (VAC[0, II - 1, j - 1, KP - 1] <= P and VAC[0, II - 1, j, KP - 1] >= P):
                        if not KDONE[k - 1]:
                            print(RSAV * np.sqrt(1 - (ETA1 * j / NV) ** 2), -j * DELV[II - 1])
                            KDONE[k - 1] = True
                            break
