def pcent(JJ, VAC, SEM, VSINT, NRDIM, NVDIM, NSDIM, NPDIM, NP):
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


def potcut3(icut, vac, tip, sem, vsint, nrdim, nvdim, nsdim, npdim, nv, ns, np, sep, s, delv, pot0, bias, chi, cpot, egap, barr, prof, nbarr1, nvdim1, nvdim2, iwrit):
    if icut == 0:
        pot0 = pcent(0, vac, sem, vsint, nrdim, nvdim, nsdim, npdim, np)
    else:
        pot0 = vsint[1, icut, 1]

    nbarr1 = 0
    barr[0] = chi + egap + pot0

    for j in range(nv):
        if tip[0, j, 0]:
            break

        if icut == 0:
            barr[j + 1] = chi + egap + \
                pcent(j, vac, sem, vsint, nrdim, nvdim, nsdim, npdim, np)
        else:
            z = sep * j / nv
            jp = int(z / delv[icut])
            f = (z - jp * delv[icut]) / delv[icut]
            if jp == 0:
                barr[j + 1] = chi + egap + vsint[1, icut, 0] * \
                    (1. - f) + vac[1, icut, jp + 1, 0] * f
            else:
                barr[j + 1] = chi + egap + vac[1, icut, jp, 0] * \
                    (1. - f) + vac[1, icut, jp + 1, 0] * f

    barr[j + 1] = chi + egap + (bias + cpot)
    nbarr1 = j + 1

    for j in range(ns):
        if icut == 0:
            prof[j] = pcent(-j, vac, sem, vsint, nrdim,
                            nvdim, nsdim, npdim, np)
        else:
            prof[j] = sem[1, icut, j, 0]

    return barr, prof, nbarr1
