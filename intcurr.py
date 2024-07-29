import numpy as np
from potexpand import *

# Example parameters, you need to adjust these based on your specific needs
impot = 1
nbarr1 = 100
sep = 1.0
pot0 = 1.0
expani = 1.0
pmax = 1.0
effm = 1.0
tk1 = 1.0
tk2 = 1.0
bias = 1.0
ef = 1.0
ne = 10
nwk = 10
ev = 1.0
eftip = 1.0
e2band = 1.0
lun = 10
iwrit = 4
icomp = 1

nvdim1 = 100
nvdim2 = 200
nsdim2 = 200
nbarr2 = 100
ns2 = 200
nsp = 200
vacstep = 0.1
jsem = np.zeros(nsdim2, dtype=int)
nexsem = np.zeros(nsp, dtype=int)
prof2 = np.zeros(nsdim2)
s2 = np.zeros(nsdim2)
psivac = np.zeros(nvdim2)
psisam = np.zeros(ns2)
wksem = np.zeros(nvdim2)
nexvac = np.zeros(nvdim2, dtype=int)
psisem = np.zeros(ns2)

barr = np.zeros(nvdim1)
s = np.zeros(nsp)
prof = np.zeros(nsp)
cdesem = np.zeros(nsp)
cdlsem = np.zeros(nsp)
cdevac = np.zeros(nvdim1)
cdlvac = np.zeros(nvdim1)

vbprof = np.zeros(nsp)
cbprof = np.zeros(nsp)


def fd(energy, ef, tk):
    """Fermi energy function"""
    return 1.0 / (1.0 + np.exp((energy - ef) / tk))


def vbwf(impot, e, wkparr, sep, bias, barr2,
         prof2, s2, effm, ev, e2band):
    """Integrate VB wavefunction from tip to sample"""
    C = 26.254
    pi = 4.0 * np.arctan(1.0)

    wf = 0.0
    wfderiv = 0.0
    wksem = 1.0
    eperp = e + wkparr**2 / (C * effm)

    if eperp >= ev:
        return wf, wfderiv

    psi = 1.0
    psivac[nbarr2 - 1] = psi
    imin = 0
    imax = nbarr2 - 1

    if eperp < barr2[0] and eperp < barr2[nbarr2 - 1]:
        for i in range(nbarr2):
            if eperp < barr2[i]:
                imin = i
                break

        for i in range(nbarr2 - 1, -1, -1):
            psivac[i] = psi
            if eperp < barr2[i]:
                imax = i
                break

        if imax <= imin:
            print('*** error - eperp above vacuum barrier')
            return wf, wfderiv

    dpsi = psi * np.sqrt(C * (barr2[imax] - eperp))
    wf = psi
    wfderiv = dpsi

    delvac = sep / float(nbarr2 - 1)
    for i in range(imax - 1, -1, -1):
        if impot == 1 and (barr2[i] - eperp) <= 0:
            continue
        psi += dpsi * delvac
        psivac[i] = psi
        dpsi += C * (barr2[i] - eperp) * psi * delvac

    psi = psi
    dpsi *= effm

    eperp = e + wkparr**2 / (C * effm)
    psi += dpsi * s2[0]
    psisem[0] = psi
    ebarr = eperp - prof2[0]
    e2band = 1.0
    if ebarr > 0:
        ebarr -= ebarr**2 / e2band

    dpsi += C * effm * ebarr * psi * s2[0]

    for i in range(1, ns2):
        dels = s2[i] - s2[i - 1]
        psi += dpsi * dels

        if np.abs(psi) >= 1e100:
            wf = 0.0
            wfderiv = 0.0
            return wf, wfderiv

        psisem[i] = psi
        ebarr = eperp - prof2[i]
        e2band = 1.0
        if ebarr > 0:
            ebarr -= ebarr**2 / e2band

        dpsi += C * effm * ebarr * psi * dels

    wksem = np.sqrt(C * effm * (ev - eperp))
    phase = np.arctan(psi * wksem / dpsi)
    amp = np.sqrt(2.0) * np.sin(phase) / psi
    wf *= amp
    wfderiv *= amp

    for i in range(nbarr2 - 1, -1, -1):
        psivac[i] *= amp

    for i in range(ns2):
        psisem[i] *= amp

    delsmax = s2[ns2 - 1] - s2[ns2 - 2]
    if delsmax / (2.0 * pi / wksem) > 0.25:
        print('*** CAUTION *** RATIO OF SEMICOND. STEP SIZE TO WAVELENGTH = ',
              delsmax / (2.0 * pi / wksem))

    return wf, wfderiv


def cbwf(impot, e, wkparr, sep, bias, barr2,
         prof2, s2, effm, ec, e2band):
    """Integrate CB wavefunction from tip to sample"""
    C = 26.254
    pi = 4.0 * np.arctan(1.0)

    wf = 0.0
    wfderiv = 0.0
    wksem = 1.0
    eperp = e - wkparr**2 / (C * effm+ 1e-10)

    if eperp <= ec:
        return wf, wfderiv

    psi = 1.0
    psivac[nbarr2 - 1] = psi
    imin = 0
    imax = nbarr2 - 1

    if eperp < barr2[0] and eperp < barr2[nbarr2 - 1]:
        for i in range(nbarr2):
            if eperp < barr2[i]:
                imin = i
                break

        for i in range(nbarr2 - 1, -1, -1):
            psivac[i] = psi
            if eperp < barr2[i]:
                imax = i
                break

        if imax <= imin:
            print('*** error - eperp above vacuum barrier')
            return wf, wfderiv

    dpsi = psi * np.sqrt(C * (barr2[imax] - eperp))
    wf = psi
    wfderiv = dpsi

    delvac = sep / float(nbarr2 - 1)
    for i in range(imax - 1, -1, -1):
        if impot == 1 and (barr2[i] - eperp) <= 0:
            continue
        psi += dpsi * delvac
        psivac[i] = psi
        dpsi += C * (barr2[i] - eperp) * psi * delvac

    psi = psi
    dpsi *= effm

    eperp = e - wkparr**2 / (C * effm)
    psi += dpsi * s2[0]
    psisem[0] = psi
    ebarr = prof2[0] - eperp

    if ebarr > 0:
        ebarr -= ebarr**2 / e2band

    dpsi += C * effm * ebarr * psi * s2[0]

    for i in range(1, ns2):
        dels = s2[i] - s2[i - 1]
        psi += dpsi * dels

        if np.abs(psi) >= 1e100:
            wf = 0.0
            wfderiv = 0.0
            return wf, wfderiv

        psisem[i] = psi
        ebarr = prof2[i] - eperp

        if ebarr > 0:
            ebarr -= ebarr**2 / e2band

        dpsi += C * effm * ebarr * psi * dels

    wksem = np.sqrt(C * effm * (eperp - ec))
    phase = np.arctan(psi * wksem / dpsi)
    amp = np.sqrt(2.0) * np.sin(phase) / psi
    wf *= amp
    wfderiv *= amp

    for i in range(nbarr2 - 1, -1, -1):
        psivac[i] *= amp

    for i in range(ns2):
        psisem[i] *= amp

    delsmax = s2[ns2 - 1] - s2[ns2 - 2]
    if delsmax / (2.0 * pi / wksem) > 0.25:
        print('*** CAUTION *** RATIO OF SEMICOND. STEP SIZE TO WAVELENGTH = ',
              delsmax / (2.0 * pi / wksem))

    return wf, wfderiv


def vbloc(impot, e, wkparr, sep, bias, barr2,
          prof2, s2, effm, ev, e2band):
    """Integrate localized VB wavefunction from tip to sample, looking for switches in sign of wf to enumerate localized states"""
    C = 26.254

    wf = 0.0
    wfderiv = 0.0
    isav = ns2
    eperp = e + wkparr**2 / (C * effm)
    if eperp <= ev:
        return wf, wfderiv, 0

    nsign = 0

    psi = 1.0
    psivac[nbarr2 - 1] = psi
    imin = 0
    imax = nbarr2 - 1

    if eperp < barr2[0] and eperp < barr2[nbarr2 - 1]:
        for i in range(nbarr2):
            if eperp < barr2[i]:
                imin = i
                break

        for i in range(nbarr2 - 1, -1, -1):
            psivac[i] = psi
            if eperp < barr2[i]:
                imax = i
                break

        if imax <= imin:
            print('*** error - eperp above vacuum barrier')
            return wf, wfderiv, nsign

    dpsi = psi * np.sqrt(C * (barr2[imax] - eperp))
    wf = psi
    wfderiv = dpsi

    delvac = sep / float(nbarr2 - 1)
    for i in range(imax - 1, -1, -1):
        if impot == 1 and (barr2[i] - eperp) <= 0:
            continue
        psi += dpsi * delvac
        psivac[i] = psi
        dpsi += C * (barr2[i] - eperp) * psi * delvac

    psi = psi
    dpsi *= effm

    eperp = e + wkparr**2 / (C * effm)
    psi += dpsi * s2[0]
    psisem[0] = psi
    ebarr = eperp - prof2[0]

    if ebarr > 0:
        ebarr -= ebarr**2 / e2band

    dpsi += C * effm * ebarr * psi * s2[0]

    for i in range(1, ns2 - 1):
        dels = s2[i] - s2[i - 1]
        psisav = psi
        psi += dpsi * dels
        psisem[i] = psi

        if psisav * psi < 0:
            nsign += 1
            isav = i

        ebarr = eperp - prof2[i]

        if ebarr > 0:
            ebarr -= ebarr**2 / e2band

        dpsi += C * effm * ebarr * psi * dels

    return wf, wfderiv, nsign


def cbloc(impot, e, wkparr, sep, bias, barr2,
          prof2, s2, effm, ec, e2band):
    """Integrate localized CB wavefunction from tip to sample, looking for switches in sign of wf to enumerate localized states"""
    C = 26.254

    wf = 0.0
    wfderiv = 0.0
    isav = ns2
    eperp = e - wkparr**2 / (C * effm)
    if eperp >= ec:
        return wf, wfderiv, 0

    nsign = 0

    psi = 1.0
    psivac[nbarr2 - 1] = psi
    imin = 0
    imax = nbarr2 - 1

    if eperp < barr2[0] and eperp < barr2[nbarr2 - 1]:
        for i in range(nbarr2):
            if eperp < barr2[i]:
                imin = i
                break

        for i in range(nbarr2 - 1, -1, -1):
            psivac[i] = psi
            if eperp < barr2[i]:
                imax = i
                break

        if imax <= imin:
            print('*** error - eperp above vacuum barrier')
            return wf, wfderiv, nsign

    dpsi = psi * np.sqrt(C * (barr2[imax] - eperp))
    wf = psi
    wfderiv = dpsi

    delvac = sep / float(nbarr2 - 1)
    for i in range(imax - 1, -1, -1):
        if impot == 1 and (barr2[i] - eperp) <= 0:
            continue
        psi += dpsi * delvac
        psivac[i] = psi
        dpsi += C * (barr2[i] - eperp) * psi * delvac

    psi = psi
    dpsi *= effm

    eperp = e - wkparr**2 / (C * effm)
    psi += dpsi * s2[0]
    psisem[0] = psi
    ebarr = prof2[0] - eperp

    if ebarr > 0:
        ebarr -= ebarr**2 / e2band

    dpsi += C * effm * ebarr * psi * s2[0]

    for i in range(1, ns2 - 1):
        dels = s2[i] - s2[i - 1]
        psisav = psi
        psi += dpsi * dels
        psisem[i] = psi

        if psisav * psi < 0:
            nsign += 1
            isav = i

        ebarr = prof2[i] - eperp

        if ebarr > 0:
            ebarr -= ebarr**2 / e2band

        dpsi += C * effm * ebarr * psi * dels

    return wf, wfderiv, nsign


def intcurr(barr, prof, nbarr1, nv, ns, nsp,
            nvdim, nsdim, s, sep, bias, ef, chi, eftip, cpot, egap, tk, avbh, avbl,
            avbso, acb, eso, e2hh, e2so, nee, nwk, pot0, nvdim1, nvdim2, nsdim2,
            expani, nloc, currv, currv0, currc, currc0, curr, curr0, iwrit,
            icomp, cdesem, cdesurf, cdlsem, cdlsurf, cdevac, cdlvac):
    """Main subroutine translated from FORTRAN to Python"""
    pi = 4.0 * np.arctan(1.0)
    rquant = 12900.0
    tk1 = tk
    tk2 = tk
    icomp = np.isscalar(icomp)
    if icomp == 1:
        cdesurf.fill(0.0)
        cdlsurf.fill(0.0)
        cdesem.fill(0.0)
        cdlsem.fill(0.0)
        cdevac.fill(0.0)
        cdlvac.fill(0.0)

    s = np.asarray(s)
    if s.size == 0:
        raise ValueError("s array must have a size greater than 0")
    
    # 確保 nsdim 是一個整數且有效
    nsdim = int(nsdim) if isinstance(nsdim, (int, float)) and nsdim > 0 else s.size

    vbprof = prof[:nsdim] + np.array([vbedge(1.0) for _ in s[:nsdim]])
    pmax = np.max(vbprof)
    ev = prof[ns] + vbedge(s[ns])

    if iwrit == 4:
        print(f'MAXIMUM ENERGY IN VB PROFILE = {pmax}')
        print(f'VB EDGE ENERGY DEEP INSIDE SEMICONDUCTOR = {ev}')

    if iwrit >= 1:
        with open("vbprof_output.txt", "w") as f:
            for j in range(nsp):
                f.write(f'{s[j]} {vbprof[j]}\n')
    ns = 0  # 初始化 ns，根據需要調整
    barr2 = np.zeros(nvdim2)  # 初始
    lun = 30
    currv = np.zeros(3)  # 假設需要 3 個元素
    currv0 = np.zeros(3)  # 假設需要 3 個元素
    currc = np.zeros(1)  # 假設需要 1 個元素
    currc0 = np.zeros(1)  # 假設需要 1 個元素
    curr = np.zeros(1)  # 假設需要 1 個元素
    curr0 = np.zeros(1)
    vbcurr1(impot, barr, nbarr1, s, prof, nsp, nsdim, nvdim, nv, nvdim1,
            nsdim2, nvdim2, sep, pot0, expani, pmax, avbl, tk1, tk2, bias,
            ef, nee, nwk, nloc[0], ev, eftip, e2hh, lun, cdesem, cdesurf, cdlsem,
            cdlsurf, cdevac, cdlvac, currv[0], currv0[0], iwrit, icomp, vbprof)

    if iwrit != 0:
        print(f'number of VB light-hole localized states = {nloc[0]}')

    lun = 40
    vbcurr1(impot, barr, nbarr1, s, prof, nsp, nsdim, nvdim, nv, nvdim1,
            nsdim2, nvdim2, sep, pot0, expani, pmax, avbh, tk1, tk2, bias,
            ef, nee, nwk, nloc[1], ev, eftip, e2hh, lun, cdesem, cdesurf, cdlsem,
            cdlsurf, cdevac, cdlvac, currv[1], currv0[1], iwrit, icomp, vbprof)

    if iwrit != 0:
        print(f'number of VB heavy-hole localized states = {nloc[1]}')

    evso = ev - eso
    pmaxso = pmax - eso
    vbprof -= eso

    lun = 50
    vbcurr1(impot, barr, nbarr1, s, prof, nsp, nsdim, nvdim, nv, nvdim1,
            nsdim2, nvdim2, sep, pot0, expani, pmaxso, avbso, tk1, tk2, bias,
            ef, nee, nwk, nloc[2], evso, eftip, e2so, lun, cdesem, cdesurf, cdlsem,
            cdlsurf, cdevac, cdlvac, currv[2], currv0[2], iwrit, icomp, vbprof)

    if iwrit != 0:
        print(f'number of VB split-off localized states = {nloc[2]}')

    currv_sum = sum(currv)
    currv0_sum = sum(currv0)

    cbprof = prof[:nsdim] + np.array([cbedge(x) for x in s[:nsdim]])
    pmin = np.min(cbprof)
    ec = prof[ns] + cbedge(s[ns])

    if iwrit == 4:
        print(f'MINIMUM ENERGY IN CB PROFILE = {pmin}')
        print(f'CB EDGE ENERGY DEEP INSIDE SEMICONDUCTOR = {ec}')

    if iwrit >= 1:
        with open("cbprof_output.txt", "w") as f:
            for j in range(nsp):
                f.write(f'{s[j]} {cbprof[j]}\n')

    lun = 60
    cbcurr1(impot, barr, nbarr1, s, prof, nsp, nsdim, nvdim, nv, nvdim1,
            nsdim2, nvdim2, sep, pot0, expani, pmin, acb, tk1, tk2, bias,
            ef, nee, nwk, nloc[3], ec, eftip, egap, lun, cdesem, cdesurf, cdlsem,
            cdlsurf, cdevac, cdlvac, currc[0], currc0[0], iwrit, icomp, cbprof)

    if iwrit != 0:
        print(f'number of CB localized states = {nloc[3]}')

    currc_sum = sum(currc)
    currc0_sum = sum(currc0)
    curr[0] = currv_sum + currc_sum
    curr0[0] = currv0_sum + currc0_sum
    return


def vbcurr1(impot, barr, nbarr1, s, prof, nsp, nsdim, nvdim, nv,
            nvdim1, nsdim2, nvdim2, sep, pot0, expani, pmax, effm, tk1, tk2, bias,
            ef, ne, nwk, nloc, ev, eftip, e2band, lun, cdesem, cdesurf, cdlsem,
            cdlsurf, cdevac, cdlvac, currv, currv0, iwrit, icomp, vbprof):
    """Calculate the tunnel current for valence band (VB) using the Bardeen formalism and T&H approximation."""
    C = 26.254
    RQUANT = 12900.0
    PI = 4.0 * np.arctan(1.0)

    nloc = 0
    currv = 0.0
    currv0 = 0.0
    wkftip = np.sqrt(C * eftip)

    emax = ev
    if nwk != 1:
        emin = min(ef - 10.0 * tk1, ef + bias - 10.0 * tk2)
    else:
        emin = ef - 10.0 * tk1
    ne = 10    
    dele = (emax - emin) / ne
    sum_val = 0.0
    if dele <= 0.0:
        return currv, currv0, nloc
    expani = 1.0
    wkmax = np.sqrt(C * effm * (emax - emin))
    semstep = (2.0 * PI / wkmax) / expani
    kappa = np.sqrt(C * (max(barr[:nbarr1]) - emin))
    vacstep = (2.0 * PI / kappa) / expani
    pot0p = pot0 + vbedge(0)

    # Initialize barr2
    barr2 = np.zeros(nvdim2)
    
    potexpand(impot, sep, nv, pot0p, s, nsp, nsdim, barr, nbarr1, nvdim1, nvdim2, prof, prof2, nsdim2, s2, vacstep, semstep, jsem, nexsem, iwrit)
    nwk = 10
    delvac = sep / float(nvdim2 - 1)
    delwk = wkmax / nwk

    for iwky in range(nwk):
        wky = iwky * delwk
        for iwkx in range(nwk):
            wkx = iwkx * delwk
            wkparr = np.sqrt(wkx**2 + wky**2)
            eparr = wkparr**2 / (effm * C)
            nwkdeg = 8
            if iwkx == 0:
                nwkdeg //= 2
            if iwky == 0:
                nwkdeg //= 2
            if iwkx == iwky:
                nwkdeg //= 2
            if iwky > iwkx:
                continue

            for ie in range(1, ne + 1):
                ener = emax - (ie - 0.5) * dele
                if eparr >= (emax - ener):
                    continue
                occtip = fd(ener - bias, ef, tk2)
                occsem = fd(ener, ef, tk1)
                occdiff = occtip - occsem
                wf, wfderiv = vbwf(impot, ener, wkparr, sep, bias, vbprof,
                                   prof2, s2, effm, ev, e2band)
                if iwkx == 0 and iwky == 0:
                    if iwrit >= 4:
                        with open(f'{lun + 3}.txt', "w") as f:
                            for j in range(nvdim2 - 1, -1, -1):
                                f.write(f"{-(j - 1) * delvac} {psivac[j]}\n")
                            for j in range(ns2):
                                f.write(f"{s2[j]} {psisam[j]}\n")
                    if icomp == 1:
                        if ie == 1:
                            sum3 = 0.0
                            for ieocc in range(ie, ne + 1):
                                enerocc = emax - ieocc * dele
                                occsem2 = 1.0 - fd(enerocc, ef, tk1)
                                sum3 += occsem2
                        else:
                            enerocc = emax - (ie - 1) * dele
                            occsem2 = 1.0 - fd(enerocc, ef, tk1)
                            sum3 -= occsem2
                        occint = sum3 * dele
                        wksem = 1.0
                        if np.all(wksem == 0):
                            cmult = 0
                        else:
                            cmult = occint * (effm * C / (2.0 * PI)) * (dele * effm * C / (2.0 * PI * wksem))
                        sum2 = 0.0
                        for j in range(200 - 1, -1, -1):
                            if j < nvdim2:
                                tmp = (psivac[j])**2 * cmult
                                if nexvac[j] != 0 and (j - 1) % nexvac[j] == 0:
                                    jtmp = ((j - 1) // nexvac[j]) + 1
                                    if jtmp < len(cdevac):
                                        cdevac[jtmp] += tmp
                                if j != 1:
                                    sum2 += tmp * delvac
                                else:
                                    sum2 += tmp * 0.5 * delvac
                        cdesurf += sum2
                        for j in range(ns2):
                            if j < len(psisem):
                                tmp = (psisam[j])**2 * cmult
                                if jsem[j] != 0:
                                    cdesem[jsem[j]] += tmp / nexsem[jsem[j]]
                                '''
                                else:
                                    print("ERROR , INTCURR - J=0")
                                    input()
                                '''
            eperp = ener - wkparr**2 / C
            kappa = np.sqrt(C * (barr[64] - eperp))
            wksem = 1.0
            trans = 2.0 * nwkdeg * (2.0 * wf)**2 * wkftip / (wksem / effm)
            sum_val += trans * occdiff

    currv = sum_val * dele * delwk**2 / (4.0 * PI**2 * RQUANT)

    if pmax <= ev:
        return currv, currv0, nloc

    emax = pmax
    if nwk != 1:
        emin = max(ev, min(ef - 10.0 * tk1, ef + bias - 10.0 * tk2))
    else:
        emin = max(ev, ef - 10.0 * tk1)
    dele = (emax - emin) / ne
    emin2 = ef - 10.0 * tk1
    dele2 = (emax - emin2) / ne
    sum_val = 0.0
    if dele <= 0.0:
        return currv, currv0, nloc
    expani = 1.0
    wkmax = np.sqrt(C * effm * (emax - emin))
    semstep = (2.0 * PI / wkmax) / expani
    kappa = np.sqrt(C * (max(barr[:nbarr1]) - emin))
    vacstep = (2.0 * PI / kappa) / expani
    pot0p = pot0 + vbedge(0)

    # Re-initialize barr2
    barr2 = np.zeros(nvdim2)

    potexpand(impot, sep, nv, pot0p, s, nsp, nsdim, barr, nbarr1, nvdim1, nvdim2, prof, prof2, nsdim2, s2, vacstep, semstep, jsem, nexsem, iwrit)
    nwk = 10
    delvac = sep / float(nvdim2 - 1)
    delwk = wkmax / nwk

    for iwky in range(nwk):
        wky = iwky * delwk
        for iwkx in range(nwk):
            wkx = iwkx * delwk
            wkparr = np.sqrt(wkx**2 + wky**2)
            eparr = wkparr**2 / (effm * C)
            n = 0
            nsav = 0
            nwkdeg = 8
            if iwkx == 0:
                nwkdeg //= 2
            if iwky == 0:
                nwkdeg //= 2
            if iwkx == iwky:
                nwkdeg //= 2
            if iwky > iwkx:
                continue

            for ie in range(1, ne + 1):
                ener = emax - (ie - 0.5) * dele
                if eparr >= (emax - ener):
                    continue
                occtip = fd(ener - bias, ef, tk2)
                occsem = fd(ener, ef, tk1)
                occdiff = occtip - occsem
                wf, wfderiv, n = vbloc(impot, ener, wkparr, sep, bias, vbprof,
                                       prof2, s2, effm, ev, e2band)
                if n == nsav:
                    continue
                if iwkx == 0 and iwky == 0:
                    if psisem[0] != 0:
                        nloc += 1
                        if iwrit > 1:
                            print(f"VB localized state at energy {ener}")
                        if iwrit >= 1:
                            with open(f'{lun}.txt', "w") as f:
                                f.write(f"{bias} {n} {ener} {occdiff}\n")
                        if iwrit >= 2:
                            with open(f'{lun + 1}.txt', "w") as f:
                                for j in range(nvdim2 - 1, -1, -1):
                                    f.write(f"{-(j - 1) * delvac} {psivac[j]}\n")
                                for j in range(ns2):
                                    f.write(f"{s2[j]} {psisam[j]}\n")

                        if icomp == 1:
                            occint = max(0.0, ener - ef)
                            if tk1 != 0:
                                sum2 = 0.0
                                for ieocc in range(1, ne + 1):
                                    enerocc = ener - (ieocc - 0.5) * dele2
                                    if enerocc < emin2:
                                        break
                                    occsem2 = 1.0 - fd(enerocc, ef, tk1)
                                    sum2 += occsem2
                                sum2 *= dele2
                                occint = sum2

                            cmult = occint * (effm * C / (2.0 * PI))
                            sum2 = 0.0
                            for j in range(nvdim2 - 1, -1, -1):
                                if j < nvdim2:
                                    tmp = (psivac[j])**2 * cmult
                                    if nexvac[j] != 0 and (j - 1) % nexvac[j] == 0:
                                        jtmp = ((j - 1) // nexvac[j]) + 1
                                        if jtmp < len(cdlvac):
                                            cdlvac[jtmp] += tmp
                                    if iwrit >= 2:
                                        with open(f'{lun + 2}.txt', "w") as f:
                                            f.write(f"{-(j - 1) * delvac} {tmp}\n")
                                    if j != 1:
                                        sum2 += tmp * delvac
                                    else:
                                        sum2 += tmp * 0.5 * delvac
                            cdlsurf += sum2
                            for j in range(ns2):
                                if j < len(psisem):
                                    tmp = (psisam[j])**2 * cmult
                                    if iwrit >= 2:
                                        with open(f'{lun + 2}.txt', "w") as f:
                                            f.write(f"{s2[j]} {tmp}\n")
                                    if jsem[j] != 0:
                                        cdlsem[jsem[j]] += tmp / nexsem[jsem[j]]
                                    else:
                                        print("ERROR , INTCURR - J=0")
                                        input()

            eperp = ener - wkparr**2 / C
            kappa = np.sqrt(C * (barr[nvdim2] - eperp))
            trans = nwkdeg * (n - nsav) * 2.0 * (2.0 * wf)**2 * wkftip
            sum_val += trans * occdiff
            nsav = n

    currv0 = sum_val * delwk**2 / (C * 2.0 * PI * RQUANT)
    return currv, currv0, nloc






def cbcurr1(impot, barr, nbarr1, s, prof, nsp, nsdim, nvdim, nv,
            nvdim1, nsdim2, nvdim2, sep, pot0, expani, pmin, effm, tk1, tk2, bias,
            ef, ne, nwk, nloc, ec, eftip, e2band, lun, cdesem, cdesurf, cdlsem,
            cdlsurf, cdevac, cdlvac, currc, currc0, iwrit, icomp, cbprof):
    """Calculate the tunnel current for conduction band (CB) using the Bardeen formalism and T&H approximation."""
    C = 26.254
    RQUANT = 12900.0
    PI = 4.0 * np.arctan(1.0)

    nloc = 0
    currv = 0.0
    currv0 = 0.0
    wkftip = np.sqrt(C * eftip)

    emax = ef  # Assuming emax should be defined as ef
    if nwk != 1:
        emin = min(ef - 10.0 * tk1, ef + bias - 10.0 * tk2)
    else:
        emin = ef - 10.0 * tk1
    ne = 10    
    dele = (emax - emin) / ne
    sum_val = 0.0
    if dele <= 0.0:
        return currv, currv0, nloc
    expani = 1.0
    wkmax = np.sqrt(C * effm * (emax - emin))
    wkmax = max(wkmax, 1e-10)  # Prevent divide by zero
    semstep = (2.0 * PI / wkmax) / expani
    kappa = np.sqrt(C * (max(barr[:nbarr1]) - emin))
    kappa = max(kappa, 1e-10)  # Prevent divide by zero
    vacstep = (2.0 * PI / kappa) / expani
    pot0p = pot0 + cbedge(0)

    # Initialize barr2
    barr2 = np.zeros(nvdim2)
    nwk = 10
    potexpand(impot, sep, nv, pot0p, s, nsp, nsdim, barr, nbarr1, nvdim1, nvdim2, prof, prof2, nsdim2, s2, vacstep, semstep, jsem, nexsem, iwrit)
    delvac = sep / float(nbarr2 - 1)
    delwk = wkmax / nwk

    for iwky in range(nwk):
        wky = iwky * delwk
        for iwkx in range(nwk):
            wkx = iwkx * delwk
            wkparr = np.sqrt(wkx**2 + wky**2)
            eparr = wkparr**2 / (effm * C + 1e-10)  # Prevent divide by zero
            nwkdeg = 8
            if iwkx == 0:
                nwkdeg //= 2
            if iwky == 0:
                nwkdeg //= 2
            if iwkx == iwky:
                nwkdeg //= 2
            if iwky > iwkx:
                continue

            for ie in range(1, ne + 1):
                ener = emin + (ie - 0.5) * dele
                if eparr >= (ener - emin):
                    continue
                occtip = fd(ener - bias, ef, tk2)
                occsem = fd(ener, ef, tk1)
                occdiff = occtip - occsem
                wf, wfderiv = cbwf(impot, ener, wkparr, sep, bias, cbprof,
                                   prof2, s2, effm, ec, e2band)
                if iwkx == 0 and iwky == 0:
                    if iwrit >= 4:
                        with open(f'{lun + 3}.txt', "w") as f:
                            for j in range(nbarr2 - 1, -1, -1):
                                f.write(f"{-(j - 1) * delvac} {psivac[j]}\n")
                            for j in range(ns2):
                                f.write(f"{s2[j]} {psisam[j]}\n")
                    if icomp == 1:
                        if ie == 1:
                            sum3 = 0.0
                            for ieocc in range(ie, ne + 1):
                                enerocc = emin + ieocc * dele
                                occsem2 = fd(enerocc, ef, tk1)
                                sum3 += occsem2
                        else:
                            enerocc = emin + (ie - 1) * dele
                            occsem2 = fd(enerocc, ef, tk1)
                            sum3 -= occsem2
                        occint = -sum3 * dele
                        wksem = max(np.sqrt(C * effm * (emax - emin)), 1e-10)  # Prevent divide by zero
                        cmult = 0 if np.all(wksem == 0) else occint * (effm * C / (2.0 * PI)) * (dele * effm * C / (2.0 * PI * wksem))
                        sum2 = 0.0
                        for j in range(nbarr2 - 1, -1, -1):
                            tmp = (psivac[j])**2 * cmult
                            if nexvac[j] != 0 and (j - 1) % nexvac[j] == 0:
                                jtmp = ((j - 1) // nexvac) + 1
                                cdevac[jtmp] += tmp
                            if j != 1:
                                sum2 += tmp * delvac
                            else:
                                sum2 += tmp * 0.5 * delvac
                        cdesurf += sum2
                        for j in range(ns2):
                            tmp = (psisam[j])**2 * cmult
                            if jsem[j] != 0:
                                cdesem[jsem[j]] += tmp / nexsem[jsem[j]]
                            """
                            else:
                                print("ERROR , INTCURR - J=0")
                                input()
                            """

            eperp = ener - wkparr**2 / C
            kappa = np.sqrt(C * (barr[64] - eperp))
            kappa = max(kappa, 1e-10)  # Prevent invalid value in sqrt
            wksem = max(np.sqrt(C * effm * (emax - emin)), 1e-10)  # Prevent divide by zero
            trans = 2.0 * nwkdeg * (2.0 * wf)**2 * wkftip / (wksem / effm)
            sum_val += trans * occdiff

    currc = sum_val * dele * delwk**2 / (4.0 * PI**2 * RQUANT)

    if pmin >= ec:
        return currc, currc0, nloc

    emin = pmin
    if nwk != 1:
        emax = min(ec, max(ef + 10.0 * tk1, ef + bias + 10.0 * tk2))
    else:
        emax = min(ec, ef + 10.0 * tk1)
    dele = (emax - emin) / ne
    emax2 = ef + 10.0 * tk1
    dele2 = (emax2 - emin) / ne
    sum_val = 0.0
    if dele <= 0.0:
        return currc, currc0, nloc

    wkmax = np.sqrt(C * effm * (emax - emin))
    wkmax = max(wkmax, 1e-10)  # Prevent divide by zero
    semstep = (2.0 * PI / wkmax) / expani
    kappa = np.sqrt(C * (max(barr[:nbarr1]) - emin))
    kappa = max(kappa, 1e-10)  # Prevent divide by zero
    vacstep = (2.0 * PI / kappa) / expani
    pot0p = pot0 + cbedge(0)
    potexpand(impot, sep, nv, pot0p, s, nsp, nsdim, barr, nbarr1, nvdim1, nvdim2, prof, prof2, nsdim2, s2, vacstep, semstep, jsem, nexsem, iwrit)
    nwk = 10
    delvac = sep / float(nbarr2 - 1)
    delwk = wkmax / nwk

    for iwky in range(nwk):
        wky = iwky * delwk
        for iwkx in range(nwk):
            wkx = iwkx * delwk
            wkparr = np.sqrt(wkx**2 + wky**2)
            eparr = wkparr**2 / (effm * C + 1e-10)  # Prevent divide by zero
            n = 0
            nsav = 0
            nwkdeg = 8
            if iwkx == 0:
                nwkdeg //= 2
            if iwky == 0:
                nwkdeg //= 2
            if iwkx == iwky:
                nwkdeg //= 2
            if iwky > iwkx:
                continue

            for ie in range(1, ne + 1):
                ener = emin + (ie - 0.5) * dele
                if eparr >= (ener - emin):
                    continue
                occtip = fd(ener - bias, ef, tk2)
                occsem = fd(ener, ef, tk1)
                occdiff = occtip - occsem
                wf, wfderiv, n = cbloc(impot, ener, wkparr, sep, bias, cbprof,
                                       prof2, s2, effm, ec, e2band)
                if n == nsav:
                    continue
                if iwkx == 0 and iwky == 0:
                    if psisem[0] != 0:
                        nloc += 1
                        if iwrit > 1:
                            print(f"CB localized state at energy {ener}")
                        if iwrit >= 1:
                            with open(f'{lun}.txt', "w") as f:
                                f.write(f"{bias} {n} {ener} {occdiff}\n")
                        if iwrit >= 2:
                            with open(f'{lun + 1}.txt', "w") as f:
                                for j in range(nbarr2 - 1, -1, -1):
                                    f.write(f"{-(j - 1) * delvac} {psivac[j]}\n")
                                for j in range(ns2):
                                    f.write(f"{s2[j]} {psisam[j]}\n")

                        if icomp == 1:
                            occint = -max(0.0, ef - ener)
                            if tk1 != 0:
                                sum2 = 0.0
                                for ieocc in range(1, ne + 1):
                                    enerocc = ener + (ieocc - 0.5) * dele2
                                    if enerocc > emax2:
                                        break
                                    occsem2 = fd(enerocc, ef, tk1)
                                    sum2 += occsem2
                                sum2 *= dele2
                                occint = -sum2

                            cmult = occint * (effm * C / (2.0 * PI))
                            sum2 = 0.0
                            for j in range(nbarr2 - 1, -1, -1):
                                tmp = (psivac[j])**2 * cmult
                                if (j - 1) % nexvac == 0:
                                    jtmp = ((j - 1) // nexvac) + 1
                                    cdlvac[jtmp] += tmp
                                if iwrit >= 2:
                                    with open(f'{lun + 2}.txt', "w") as f:
                                        f.write(f"{-(j - 1) * delvac} {tmp}\n")
                                if j != 1:
                                    sum2 += tmp * delvac
                                else:
                                    sum2 += tmp * 0.5 * delvac
                            cdlsurf += sum2
                            for j in range(ns2):
                                tmp = (psisam[j])**2 * cmult
                                if iwrit >= 2:
                                    with open(f'{lun + 2}.txt', "w") as f:
                                        f.write(f"{s2[j]} {tmp}\n")
                                if jsem[j] != 0:
                                    cdlsem[jsem[j]] += tmp / nexsem[jsem[j]]
                                else:
                                    print("ERROR , INTCURR - J=0")
                                    input()

            eperp = ener - wkparr**2 / C
            kappa = np.sqrt(C * (barr[64] - eperp))
            kappa = max(kappa, 1e-10)  # Prevent invalid value in sqrt
            wksem = max(np.sqrt(C * effm * (emax - emin)), 1e-10)  # Prevent divide by zero
            trans = nwkdeg * (n - nsav) * 2.0 * (2.0 * wf)**2 * wkftip
            sum_val += trans * occdiff
            nsav = n

    currc0 = sum_val * delwk**2 / (C * 2.0 * PI * RQUANT)
    return currc, currc0, nloc





def vbedge(z):
    """Function to calculate the vacuum band edge."""
    if z < 0:
        return 0.0
    return 1.0


def cbedge(z):
    """Function to calculate the conduction band edge."""
    return np.where(z > 0, 1.0, 0.0)

# Example call to intcurr (you need to provide actual values for the parameters)