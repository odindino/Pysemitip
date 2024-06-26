import numpy as np
from semirhomult import fd


def getcurr(iband, bias, sep, e0, barr, tk1, tk2, ef, epsi2, nstates, eftip, nbarr2, nvdim2, barrprof, neigendim):
    c = 26.254  # 2m/hbar^2 in units of 1/(eV nm^2)
    rquant = 12900.0
    pi = np.pi
    wkftip = np.sqrt(c * eftip)

    total_sum = 0.0

    for i in range(nstates):
        ener = e0 + epsi2[0, i]
        occ = fd(ener - bias, ef, tk2) - fd(ener, ef, tk1)
        sum2 = 0.0
        delvac = sep / float(nbarr2 - 1)

        for j in range(nbarr2):
            term = np.sqrt(c * max(0.0, barrprof[j] - ener) + epsi2[3, i])
            if j == 0 or j == nbarr2 - 1:
                term = term / 2.0
            sum2 += term

        skapka = sum2 * delvac
        trans = epsi2[1, i] * np.exp(-2.0 * skapka)
        total_sum += trans * occ

    curr = 16.0 * pi * wkftip * total_sum / (c * rquant)
    return curr
