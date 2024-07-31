import numpy as np
from potexpand import potexpand

# Example parameters, adjust these based on your specific needs
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

def intcurr(impot, e, wkparr, sep, barr2, s2, effm, ev, e2band, is_cb=False):
    """Integrate wavefunction from tip to sample"""
    C = 26.254
    pi = 4.0 * np.arctan(1.0)
    wf, wfderiv = 0.0, 0.0
    eperp = e - wkparr**2 / (C * effm) if is_cb else e + wkparr**2 / (C * effm)
    if (is_cb and eperp <= ev) or (not is_cb and eperp >= ev):
        return wf, wfderiv

    psi = 1.0
    dpsi = psi * np.sqrt(C * (barr2[-1] - eperp))
    wf, wfderiv = psi, dpsi
    delvac = sep / float(len(barr2) - 1)

    for i in range(len(barr2) - 2, -1, -1):
        if impot == 1 and (barr2[i] - eperp) <= 0:
            continue
        psi += dpsi * delvac
        psivac[i] = psi
        dpsi += C * (barr2[i] - eperp) * psi * delvac

    dpsi *= effm
    eperp = e - wkparr**2 / (C * effm) if is_cb else e + wkparr**2 / (C * effm)
    psi += dpsi * s2[0]
    psisem[0] = psi
    ebarr = prof2[0] - eperp if is_cb else eperp - prof2[0]

    if ebarr > 0:
        ebarr -= ebarr**2 / e2band

    dpsi += C * effm * ebarr * psi * s2[0]

    for i in range(1, len(s2)):
        dels = s2[i] - s2[i - 1]
        psi += dpsi * dels

        if np.abs(psi) >= 1e100:
            return 0.0, 0.0

        psisem[i] = psi
        ebarr = prof2[i] - eperp if is_cb else eperp - prof2[i]

        if ebarr > 0:
            ebarr -= ebarr**2 / e2band

        dpsi += C * effm * ebarr * psi * dels

    wf *= np.sqrt(2.0) / psi
    wfderiv *= np.sqrt(2.0) / psi

    for i in range(len(barr2)):
        psivac[i] *= np.sqrt(2.0) / psi

    for i in range(len(s2)):
        psisem[i] *= np.sqrt(2.0) / psi

    return wf, wfderiv

def calculate_current(impot, barr, prof, nbarr1, nv, ns, nsp, nvdim, nsdim, s, sep, bias, ef, chi, eftip, cpot, egap, tk, avbh, avbl, avbso, acb, eso, e2hh, e2so, nee, nwk, pot0, nvdim1, nvdim2, nsdim2, expani, nloc, currv, currv0, currc, currc0, curr, curr0, iwrit, icomp, cdesem, cdesurf, cdlsem, cdlsurf, cdevac, cdlvac):
    """Main subroutine translated from FORTRAN to Python"""
    pi = 4.0 * np.arctan(1.0)
    rquant = 12900.0
    tk1, tk2 = tk, tk
    cdesurf.fill(0.0)
    cdlsurf.fill(0.0)
    cdesem.fill(0.0)
    cdlsem.fill(0.0)
    cdevac.fill(0.0)
    cdlvac.fill(0.0)

    vbprof = prof[:nsdim] + np.array([vbedge(1.0) for _ in s[:nsdim]])
    pmax = np.max(vbprof)
    ev = prof[ns-1] + vbedge(s[ns-1])

    if iwrit == 4:
        print(f'MAXIMUM ENERGY IN VB PROFILE = {pmax}')
        print(f'VB EDGE ENERGY DEEP INSIDE SEMICONDUCTOR = {ev}')

    if iwrit >= 1:
        with open("vbprof_output.txt", "w") as f:
            for j in range(nsp):
                f.write(f'{s[j]} {vbprof[j]}\n')

    ns = 0  # Initialize ns as needed
    barr2 = np.zeros(nvdim2)
    currv = np.zeros(3)
    currv0 = np.zeros(3)
    currc = np.zeros(1)
    currc0 = np.zeros(1)
    curr = np.zeros(1)
    curr0 = np.zeros(1)

    def run_vbcurr(lun, av, ev, e2, nloc_idx):
        nloc[nloc_idx] = vbcurr1(impot, barr, nbarr1, s, prof, nsp, nsdim, nvdim, nv, nvdim1, nsdim2, nvdim2, sep, pot0, expani, pmax, av, tk1, tk2, bias, ef, nee, nwk, nloc[nloc_idx], ev, eftip, e2, lun, cdesem, cdesurf, cdlsem, cdlsurf, cdevac, cdlvac, currv[nloc_idx], currv0[nloc_idx], iwrit, icomp, vbprof)

    run_vbcurr(30, avbl, ev, e2hh, 0)
    run_vbcurr(40, avbh, ev, e2hh, 1)
    run_vbcurr(50, avbso, ev - eso, e2so, 2)

    print(f'number of VB light-hole localized states = {nloc[0]}')
    print(f'number of VB heavy-hole localized states = {nloc[1]}')
    print(f'number of VB split-off localized states = {nloc[2]}')

    currv_sum = np.sum(currv)
    currv0_sum = np.sum(currv0)

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

    nloc[3] = cbcurr1(impot, barr, nbarr1, s, prof, nsp, nsdim, nvdim, nv, nvdim1, nsdim2, nvdim2, sep, pot0, expani, pmin, acb, tk1, tk2, bias, ef, nee, nwk, nloc[3], ec, eftip, egap, lun, cdesem, cdesurf, cdlsem, cdlsurf, cdevac, cdlvac, currc[0], currc0[0], iwrit, icomp, cbprof)

    print(f'number of CB localized states = {nloc[3]}')

    currc_sum = np.sum(currc)
    currc0_sum = np.sum(currc0)
    curr[0] = currv_sum + currc_sum
    curr0[0] = currv0_sum + currc0_sum
    return

def vbcurr1(impot, barr, nbarr1, s, prof, nsp, nsdim, nvdim, nv, nvdim1, nsdim2, nvdim2, sep, pot0, expani, pmax, effm, tk1, tk2, bias, ef, ne, nwk, nloc, ev, eftip, e2band, lun, cdesem, cdesurf, cdlsem, cdlsurf, cdevac, cdlvac, currv, currv0, iwrit, icomp, vbprof):
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
    if dele <= 0.0:
        return currv, currv0, nloc

    expani = 1.0
    wkmax = np.sqrt(C * effm * (emax - emin))
    semstep = (2.0 * PI / wkmax) / expani
    kappa = np.sqrt(C * (np.max(barr[:nbarr1]) - emin))
    vacstep = (2.0 * PI / kappa) / expani
    pot0p = pot0 + vbedge(0)

    potexpand(impot, sep, nv, pot0p, s, nsp, nsdim, barr, nbarr1, nvdim1, nvdim2, prof, prof2, nsdim2, s2, vacstep, semstep, jsem, nexsem, iwrit)
    delvac = sep / float(nvdim2 - 1)
    delwk = wkmax / nwk

    for iwky in range(nwk):
        wky = iwky * delwk
        for iwkx in range(nwk):
            wkx = iwkx * delwk
            wkparr = np.sqrt(wkx**2 + wky**2)
            eparr = wkparr**2 / (effm * C)
            nwkdeg = 8
            nwkdeg //= 2 if iwkx == 0 else nwkdeg
            nwkdeg //= 2 if iwky == 0 else nwkdeg
            nwkdeg //= 2 if iwkx == iwky else nwkdeg
            if iwky > iwkx:
                continue

            for ie in range(1, ne + 1):
                ener = emax - (ie - 0.5) * dele
                if eparr >= (emax - ener):
                    continue
                occtip = fd(ener - bias, ef, tk2)
                occsem = fd(ener, ef, tk1)
                occdiff = occtip - occsem
                wf, wfderiv = intcurr(impot, ener, wkparr, sep, vbprof, s2, effm, ev, e2band)
                if iwkx == 0 and iwky == 0 and ie == 1 and iwrit >= 4:
                    with open(f'{lun + 3}.txt', "w") as f:
                        for j in range(nvdim2 - 1, -1, -1):
                            f.write(f"{-(j - 1) * delvac} {psivac[j]}\n")
                        for j in range(ns2):
                            f.write(f"{s2[j]} {psisam[j]}\n")
                    if icomp == 1:
                        sum3 = np.sum([1.0 - fd(emax - (ieocc * dele), ef, tk1) for ieocc in range(ie, ne + 1)]) if ie == 1 else sum3 - (1.0 - fd(emax - ((ie - 1) * dele), ef, tk1))
                        occint = sum3 * dele
                        wksem = 1.0
                        cmult = 0 if np.all(wksem == 0) else occint * (effm * C / (2.0 * PI)) * (dele * effm * C / (2.0 * PI * wksem))
                        sum2 = np.sum([(psivac[j]**2 * cmult) if (nexvac[j] == 0 or (j - 1) % nexvac[j] != 0) else (cdevac[((j - 1) // nexvac[j]) + 1] + (psivac[j]**2 * cmult)) for j in range(nvdim2 - 1, -1, -1)]) * delvac
                        cdesurf += sum2
                        sum2 = np.sum([(psisam[j]**2) * cmult / nexsem[jsem[j]] if nexsem[jsem[j]] != 0 else 0 for j in range(ns2)])
                        cdesem += sum2

            eperp = ener - wkparr**2 / C
            kappa = np.sqrt(C * (barr[64] - eperp))
            trans = 2.0 * nwkdeg * (2.0 * wf)**2 * wkftip / (wksem / effm)
            currv += trans * occdiff

    currv = currv * dele * delwk**2 / (4.0 * PI**2 * RQUANT)

    if pmax <= ev:
        return currv, currv0, nloc

    emax = pmax
    emin = max(ev, min(ef - 10.0 * tk1, ef + bias - 10.0 * tk2)) if nwk != 1 else max(ev, ef - 10.0 * tk1)
    dele = (emax - emin) / ne
    sum_val = 0.0

    if dele <= 0.0:
        return currv, currv0, nloc

    pot0p = pot0 + vbedge(0)
    potexpand(impot, sep, nv, pot0p, s, nsp, nsdim, barr, nbarr1, nvdim1, nvdim2, prof, prof2, nsdim2, s2, vacstep, semstep, jsem, nexsem, iwrit)
    delvac = sep / float(nvdim2 - 1)
    delwk = wkmax / nwk

    for iwky in range(nwk):
        wky = iwky * delwk
        for iwkx in range(nwk):
            wkx = iwkx * delwk
            wkparr = np.sqrt(wkx**2 + wky**2)
            eparr = wkparr**2 / (effm * C)
            nsav = 0
            nwkdeg = 8
            nwkdeg //= 2 if iwkx == 0 else nwkdeg
            nwkdeg //= 2 if iwky == 0 else nwkdeg
            nwkdeg //= 2 if iwkx == iwky else nwkdeg
            if iwky > iwkx:
                continue

            for ie in range(1, ne + 1):
                ener = emax - (ie - 0.5) * dele
                if eparr >= (emax - ener):
                    continue
                occtip = fd(ener - bias, ef, tk2)
                occsem = fd(ener, ef, tk1)
                occdiff = occtip - occsem
                wf, wfderiv, n = intcurr(impot, ener, wkparr, sep, vbprof, s2, effm, ev, e2band, is_cb=False)
                if n == nsav:
                    continue
                if iwkx == 0 and iwky == 0 and psisem[0] != 0:
                    nloc += 1
                    if iwrit > 1:
                        print(f"VB localized state at energy {ener}")
                    if iwrit >= 1:
                        with open(f'{lun}.txt', "w") as f:
                            f.write(f"{bias} {n} {ener} {occdiff}\n")
                    if icomp == 1:
                        occint = max(0.0, ener - ef)
                        if tk1 != 0:
                            sum2 = np.sum([1.0 - fd(ener - (ieocc - 0.5) * dele, ef, tk1) for ieocc in range(1, ne + 1)]) * dele
                            occint = sum2

                        cmult = occint * (effm * C / (2.0 * PI))
                        sum2 = np.sum([(psivac[j]**2 * cmult) if (nexvac[j] == 0 or (j - 1) % nexvac[j] != 0) else (cdevac[((j - 1) // nexvac[j]) + 1] + (psivac[j]**2 * cmult)) for j in range(nvdim2 - 1, -1, -1)]) * delvac
                        cdlsurf += sum2
                        sum2 = np.sum([(psisam[j]**2) * cmult / nexsem[jsem[j]] if nexsem[jsem[j]] != 0 else 0 for j in range(ns2)])
                        cdlsem += sum2

            eperp = ener - wkparr**2 / C
            kappa = np.sqrt(C * (barr[64] - eperp))
            trans = nwkdeg * (n - nsav) * 2.0 * (2.0 * wf)**2 * wkftip
            currv += trans * occdiff
            nsav = n

    currv0 = currv * delwk**2 / (C * 2.0 * PI * RQUANT)
    return currv, currv0, nloc


def cbcurr1(impot, barr, nbarr1, s, prof, nsp, nsdim, nvdim, nv, nvdim1, nsdim2, nvdim2, sep, pot0, expani, pmin, effm, tk1, tk2, bias, ef, ne, nwk, nloc, ec, eftip, e2band, lun, cdesem, cdesurf, cdlsem, cdlsurf, cdevac, cdlvac, currc, currc0, iwrit, icomp, cbprof):
    """Calculate the tunnel current for conduction band (CB) using the Bardeen formalism and T&H approximation."""
    C = 26.254
    RQUANT = 12900.0
    PI = 4.0 * np.arctan(1.0)

    nloc = 0
    currc = 0.0
    currc0 = 0.0
    wkftip = np.sqrt(C * eftip)

    emax = ef
    if nwk != 1:
        emin = min(ef - 10.0 * tk1, ef + bias - 10.0 * tk2)
    else:
        emin = ef - 10.0 * tk1

    ne = 10    
    dele = (emax - emin) / ne
    if dele <= 0.0:
        return currc, currc0, nloc

    expani = 1.0
    wkmax = np.sqrt(C * effm * (emax - emin))
    semstep = (2.0 * PI / wkmax) / expani
    kappa = np.sqrt(C * (np.max(barr[:nbarr1]) - emin))
    vacstep = (2.0 * PI / kappa) / expani
    pot0p = pot0 + cbedge(0)

    potexpand(impot, sep, nv, pot0p, s, nsp, nsdim, barr, nbarr1, nvdim1, nvdim2, prof, prof2, nsdim2, s2, vacstep, semstep, jsem, nexsem, iwrit)
    delvac = sep / float(nvdim2 - 1)
    delwk = wkmax / nwk

    for iwky in range(nwk):
        wky = iwky * delwk
        for iwkx in range(nwk):
            wkx = iwkx * delwk
            wkparr = np.sqrt(wkx**2 + wky**2)
            eparr = wkparr**2 / (effm * C)
            nwkdeg = 8
            nwkdeg //= 2 if iwkx == 0 else nwkdeg
            nwkdeg //= 2 if iwky == 0 else nwkdeg
            nwkdeg //= 2 if iwkx == iwky else nwkdeg
            if iwky > iwkx:
                continue

            for ie in range(1, ne + 1):
                ener = emin + (ie - 0.5) * dele
                if eparr >= (ener - emin):
                    continue
                occtip = fd(ener - bias, ef, tk2)
                occsem = fd(ener, ef, tk1)
                occdiff = occtip - occsem
                wf, wfderiv = intcurr(impot, ener, wkparr, sep, cbprof, s2, effm, ec, e2band, is_cb=True)
                if iwkx == 0 and iwky == 0 and ie == 1 and iwrit >= 4:
                    with open(f'{lun + 3}.txt', "w") as f:
                        for j in range(nvdim2 - 1, -1, -1):
                            f.write(f"{-(j - 1) * delvac} {psivac[j]}\n")
                        for j in range(ns2):
                            f.write(f"{s2[j]} {psisam[j]}\n")
                    if icomp == 1:
                        sum3 = np.sum([fd(emin + (ieocc * dele), ef, tk1) for ieocc in range(ie, ne + 1)]) if ie == 1 else sum3 - fd(emin + ((ie - 1) * dele), ef, tk1)
                        occint = -sum3 * dele
                        wksem = 1.0
                        cmult = 0 if np.all(wksem == 0) else occint * (effm * C / (2.0 * PI)) * (dele * effm * C / (2.0 * PI * wksem))
                        sum2 = np.sum([(psivac[j]**2 * cmult) if (nexvac[j] == 0 or (j - 1) % nexvac[j] != 0) else (cdevac[((j - 1) // nexvac[j]) + 1] + (psivac[j]**2 * cmult)) for j in range(nvdim2 - 1, -1, -1)]) * delvac
                        cdesurf += sum2
                        sum2 = np.sum([(psisam[j]**2) * cmult / nexsem[jsem[j]] if nexsem[jsem[j]] != 0 else 0 for j in range(ns2)])
                        cdesem += sum2

            eperp = ener - wkparr**2 / C
            kappa = np.sqrt(C * (barr[64] - eperp))
            trans = 2.0 * nwkdeg * (2.0 * wf)**2 * wkftip / (wksem / effm)
            currc += trans * occdiff

    currc = currc * dele * delwk**2 / (4.0 * PI**2 * RQUANT)

    if pmin >= ec:
        return currc, currc0, nloc

    emin = pmin
    emax = min(ec, max(ef + 10.0 * tk1, ef + bias + 10.0 * tk2)) if nwk != 1 else min(ec, ef + 10.0 * tk1)
    dele = (emax - emin) / ne
    sum_val = 0.0

    if dele <= 0.0:
        return currc, currc0, nloc

    pot0p = pot0 + cbedge(0)
    potexpand(impot, sep, nv, pot0p, s, nsp, nsdim, barr, nbarr1, nvdim1, nvdim2, prof, prof2, nsdim2, s2, vacstep, semstep, jsem, nexsem, iwrit)
    delvac = sep / float(nvdim2 - 1)
    delwk = wkmax / nwk

    for iwky in range(nwk):
        wky = iwky * delwk
        for iwkx in range(nwk):
            wkx = iwkx * delwk
            wkparr = np.sqrt(wkx**2 + wky**2)
            eparr = wkparr**2 / (effm * C)
            nsav = 0
            nwkdeg = 8
            nwkdeg //= 2 if iwkx == 0 else nwkdeg
            nwkdeg //= 2 if iwky == 0 else nwkdeg
            nwkdeg //= 2 if iwkx == iwky else nwkdeg
            if iwky > iwkx:
                continue

            for ie in range(1, ne + 1):
                ener = emin + (ie - 0.5) * dele
                if eparr >= (ener - emin):
                    continue
                occtip = fd(ener - bias, ef, tk2)
                occsem = fd(ener, ef, tk1)
                occdiff = occtip - occsem
                wf, wfderiv, n = intcurr(impot, ener, wkparr, sep, cbprof, s2, effm, ec, e2band, is_cb=True)
                if n == nsav:
                    continue
                if iwkx == 0 and iwky == 0 and psisem[0] != 0:
                    nloc += 1
                    if iwrit > 1:
                        print(f"CB localized state at energy {ener}")
                    if iwrit >= 1:
                        with open(f'{lun}.txt', "w") as f:
                            f.write(f"{bias} {n} {ener} {occdiff}\n")
                    if icomp == 1:
                        occint = max(0.0, ef - ener)
                        if tk1 != 0:
                            sum2 = np.sum([fd(emin + (ieocc - 0.5) * dele, ef, tk1) for ieocc in range(1, ne + 1)]) * dele
                            occint = sum2

                        cmult = occint * (effm * C / (2.0 * PI))
                        sum2 = np.sum([(psivac[j]**2 * cmult) if (nexvac[j] == 0 or (j - 1) % nexvac[j] != 0) else (cdevac[((j - 1) // nexvac[j]) + 1] + (psivac[j]**2 * cmult)) for j in range(nvdim2 - 1, -1, -1)]) * delvac
                        cdlsurf += sum2
                        sum2 = np.sum([(psisam[j]**2) * cmult / nexsem[jsem[j]] if nexsem[jsem[j]] != 0 else 0 for j in range(ns2)])
                        cdlsem += sum2

            eperp = ener - wkparr**2 / C
            kappa = np.sqrt(C * (barr[64] - eperp))
            trans = nwkdeg * (n - nsav) * 2.0 * (2.0 * wf)**2 * wkftip
            currc += trans * occdiff
            nsav = n

    currc0 = currc * delwk**2 / (C * 2.0 * PI * RQUANT)
    return currc, currc0, nloc

def vbedge(s):
    return 0.0  # Implement this function based on your requirements

def cbedge(s):
    return 0.0  # Implement this function based on your requirements
