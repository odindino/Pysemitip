import numpy as np
def potexpand(IMPOT, SEP, NV, POT0P, S, NS, NSDIM, BARR, NBARR1, BARR2, NBARR2, NVDIM1, NVDIM2, PROF, PROF2, NSDIM2, S2, NS2, VACSTEP, SEMSTEP, JSEM, NEXSEM, NEXVAC, IWRIT):
    # Initialize common variables
    KAPPA = 0.0
    LAMBDA_VAL = 0.0
    EGAP, ED, EA, ACB, AVB, CD, CA, EPSIL, TK, IDEG, IINV = 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0

    # Expand vacuum barrier
    NEXPAN = max(1, int(round((SEP / NV) / VACSTEP)))
    if IMPOT == 1:
        NEXPAN *= 10
    NEXVAC = NEXPAN
    if IWRIT > 1:
        print(f'expansion factor for barrier = {NEXPAN}')
    
    NBARR2 = NEXPAN * (NBARR1 - 1) + 1
    BARR2[NBARR2 - 1] = BARR[NBARR1 - 1]
    for J in range(NBARR1 - 1, 0, -1):
        B2 = BARR[J]
        B1 = BARR[J - 1]
        for K in range(NEXPAN - 1, -1, -1):
            BARR2[(J - 1) * NEXPAN + K] = (B2 * K + B1 * (NEXPAN - K)) / NEXPAN
    
    if IWRIT >= 3:
        for I in range(NBARR2):
            print(-(I) * SEP / float(NBARR2 - 1), BARR2[I])

    if IWRIT > 1:
        print(f'number of expanded points in vacuum = {NBARR2}')

    LAMBDA_VAL = 3.81 ** 2 * 0.1 * np.log(2) / (2 * 2 * SEP)
    if IMPOT == 1:
        for J in range(1, NBARR2 - 1):
            BARR2[J] = BARR2[J] - 1.15 * LAMBDA_VAL * (NBARR2 - 1) ** 2 / ((J) * (float(NBARR2) - J))
    
    if IWRIT >= 3:
        for I in range(NBARR2):
            print(-(I) * SEP / float(NBARR2 - 1), BARR2[I])

    # Expand the potential profile in semiconductor
    NEXSEM.fill(0)
    NEXPAN = max(1, int(round(2 * S[0] / SEMSTEP)))
    if IWRIT > 1:
        print(f'initial expansion factor for semiconductor = {NEXPAN}')

    KK = 0
    for J in range(NS):
        if J == 0:
            NEXPAN = max(1, int(round(S[0] / SEMSTEP)))
        else:
            NEXPAN = max(1, int(round((S[J] - S[J - 1]) / SEMSTEP)))
        
        if NEXPAN % 2 == 0:
            NEXPAN += 1
        
        for K in range(1, NEXPAN + 1):
            KK += 1
            if J == 0:
                JSEM[KK - 1] = J + 1
            else:
                if K <= (NEXPAN / 2):
                    JSEM[KK - 1] = J
                else:
                    JSEM[KK - 1] = J + 1
            NEXSEM[JSEM[KK - 1] - 1] += 1
            if KK > NSDIM2:
                print('*** ERROR - NSDIM2 TOO SMALL ', KK, NSDIM2)
                print('PRESS THE ENTER KEY TO EXIT')
                input()
                exit()

            if J == 0:
                PROF2[KK - 1] = ((NEXPAN - K) * POT0P + K * PROF[J]) / NEXPAN
            else:
                PROF2[KK - 1] = ((NEXPAN - K) * PROF[J - 1] + K * PROF[J]) / NEXPAN
            
            if J == 0:
                S2[KK - 1] = ((NEXPAN - K) * 0.0 + K * S[J]) / NEXPAN
            else:
                S2[KK - 1] = ((NEXPAN - K) * S[J - 1] + K * S[J]) / NEXPAN
    
    NS2 = KK
    if IWRIT > 1:
        print(f'number of expanded points in semiconductor = {NS2}')

    return NEXVAC, NBARR2, BARR2, NS2, PROF2, S2