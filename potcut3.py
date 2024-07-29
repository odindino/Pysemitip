def PCENT(JJ, VAC, SEM, VSINT, NRDIM, NVDIM, NSDIM, NPDIM, NP):
    J = abs(JJ)
    I = 1
    SUM = 0.0
    if JJ == 0:
        for K in range(1, NP + 1):
            SUM += (9.0 * VSINT[0, I, K-1] - VSINT[0, I+1, K-1]) / 8.0
    elif JJ > 0:
        for K in range(1, NP + 1):
            SUM += (9.0 * VAC[0, I, J, K-1] - VAC[0, I+1, J, K-1]) / 8.0
    else:
        for K in range(1, NP + 1):
            SUM += (9.0 * SEM[0, I, J, K-1] - SEM[0, I+1, J, K-1]) / 8.0
    return SUM / float(NP)


def potcut3(ICUT, VAC, TIP, SEM, VSINT, NRDIM, NVDIM, NSDIM, NPDIM, NV, NS, NP, SEP, S, DELV, POT0, BIAS, CHI, CPot, EGAP, BARR, PROF, NBARR1, NVDIM1, NVDIM2, IWRIT):
    if ICUT == 0:
        POT0 = PCENT(0, VAC, SEM, VSINT, NRDIM, NVDIM, NSDIM, NPDIM, NP)
    else:
        POT0 = VSINT[1, ICUT, 1]

    NBARR1 = 0
    BARR[0] = CHI + EGAP + POT0

    for J in range(NV):
        if TIP[0, J, 0]:
            break

        if ICUT == 0:
            BARR[J + 1] = CHI + EGAP + PCENT(J, VAC, SEM, VSINT, NRDIM, NVDIM, NSDIM, NPDIM, NP)
        else:
            Z = SEP * J / NV
            JP = int(Z / DELV[ICUT])
            F = (Z - JP * DELV[ICUT]) / DELV[ICUT]
            if JP == 0:
                BARR[J + 1] = CHI + EGAP + VSINT[1, ICUT, 0] * (1. - F) + VAC[1, ICUT, JP + 1, 0] * F
            else:
                BARR[J + 1] = CHI + EGAP + VAC[1, ICUT, JP, 0] * (1. - F) + VAC[1, ICUT, JP + 1, 0] * F

    BARR[J + 1] = CHI + EGAP + (BIAS + CPot)
    NBARR1 = J + 1

    for J in range(NS):
        if ICUT == 0:
            PROF[J] = PCENT(-J, VAC, SEM, VSINT, NRDIM, NVDIM, NSDIM, NPDIM, NP)
        else:
            PROF[J] = SEM[1, ICUT, J, 0]

    return BARR, PROF, NBARR1