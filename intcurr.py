import numpy as np
from potexpand import*

# Constants
C = 26.254  # 1/(eV nm^2)
RQUANT = 12900.0
PI = 4.0 * np.arctan(1.0)
nbarr2 = 100  # Example value, adjust as necessary
ns2 = 50      # Example value, adjust as necessary
nexsem = np.zeros(10, dtype=int)
nexvac = np.zeros(10, dtype=int)
b1 = 1.0  # Example coefficient for valence band edge
b2 = 0.1  # Example coefficient for valence band edge
a1 = 1.5  # Example coefficient for conduction band edge
a2 = 0.2  # Example coefficient for conduction band edge
e2hh = 1.0  # Heavy-hole band energy offset
e2so = 0.5  # Split-off band energy offset
barr2 = np.zeros(nbarr2)
prof2 = np.zeros(ns2)
s2 = np.linspace(0, 1, ns2)
jsem = np.zeros(ns2, dtype=int)
psivac = np.zeros(nbarr2)
psisem = np.zeros(ns2)

def atand(x):
    return np.degrees(np.arctan(x))

def fd(e, ef, tk):
    """Fermi-Dirac distribution."""
    k_boltz = 8.617333262145e-5  # Boltzmann constant in eV/K
    return 1.0 / (1.0 + np.exp((e - ef) / (k_boltz * tk)))

def vbcurr(impot, barr, nbarr1, s, prof, nsp, nsdim, nvdim, nv, nvdim1, nsdim2, nvdim2, sep, pot0, expani, pmax, effm, tk1, tk2, bias, ef, nee, nwk, nloc, ev, eftip, e2band, lun, cdesem, cdesurf, cdlsem, cdlsurf, cdevac, cdlvac, currv, currv0, iwrit, icomp, vbprof):
    """VB Tunnel Current"""
    nwk=1
    nee=1
    expani=1
    kappa = 0.0
    wkftip = np.sqrt(C * eftip)
    sum1 = sum2 = sum1s = sum2s = sum1p = sum2p = 0.0
    currv = 0.0
    currv0 = 0.0
    nwk=1
    nee=1
    expani=1
    emax = ev
    if nwk != 1:
        emin = min(ef - 10.0 * tk1, ef + bias - 10.0 * tk2)
    else:
        emin = ef - 10.0 * tk1
    dele = (emax - emin) / nee
    if dele <= 0:
        return currv, currv0

    wkmax = np.sqrt(C * effm * (emax - emin))
    if wkmax == 0:
        return currv, currv0  # Avoid division by zero

    semstep = (2.0 * PI / wkmax) / expani
    kappa = np.sqrt(C * max(barr[nbarr1 - 1], barr[0]) - emin)
    if kappa == 0:
        return currv, currv0  # Avoid division by zero

    vacstep = (2.0 * PI / kappa) / expani
    pot0p = pot0 + vbprof[0]
    
    potexpand(impot, sep, nv, pot0p, s, nsp, nsdim, barr, nbarr1, barr2, nbarr2, nvdim1, nvdim2, vbprof, prof2, nsdim2, s2, ns2, vacstep, semstep, jsem, nexsem, nexvac, iwrit)
    
    delvac = sep / float(nbarr2 - 1)
    delwk = wkmax / nwk

    for iwky in range(nwk):
        wky = iwky * delwk
        for iwkx in range(nwk):
            wkx = iwkx * delwk
            wkparr = np.sqrt(wkx ** 2 + wky ** 2)
            eparr = wkparr ** 2 / (effm * C)
            nwkdeg = 8
            if iwkx == 0:
                nwkdeg /= 2
            if iwky == 0:
                nwkdeg /= 2
            if iwkx == iwky:
                nwkdeg /= 2
            if iwky > iwkx:
                continue
            for ie in range(1, nee + 1):
                ener = emax - (ie - 0.5) * dele
                if eparr >= (emax - ener):
                    continue
                occtip = fd(ener - bias, ef, tk2)
                occsem = fd(ener, ef, tk1)
                occdiff = occtip - occsem
                wf, wfderiv, wksem = vbwf(impot, 0.0, 0.0, 0.0, ener, wkparr, sep, bias, barr2, nvdim2, nbarr2, prof2, ns2, nsdim2, s2, effm, ev, e2band, np.zeros(nvdim2), np.zeros(nsdim2))
                if iwkx == 0 and iwky == 0 and iwrit >= 4:
                    for j in range(nbarr2 - 1, -1, -1):
                        print(-(j - 1) * delvac, psivac[j])
                    for j in range(ns2):
                        print(s2[j], psisem[j])

                if icomp == 1:
                    if ie == 1:
                        sum3 = 0.0
                        for ieocc in range(ie, nee + 1):
                            enerocc = emax - ieocc * dele
                            occsem2 = 1.0 - fd(enerocc, ef, tk1)
                            sum3 += occsem2
                    else:
                        enerocc = emax - (ie - 1) * dele
                        occsem2 = 1.0 - fd(enerocc, ef, tk1)
                        sum3 -= occsem2
                    occint = sum3 * dele
                    cmult = occint * (effm * C / (2.0 * PI)) * (dele * effm * C / (2.0 * PI * wksem))
                    sum2 = 0.0
                    for j in range(nbarr2 - 1, -1, -1):
                        tmp = (psivac[j]) ** 2 * cmult
                        if (j - 1) % nexvac == 0:
                            jtmp = (j - 1) // nexvac + 1
                            cdevac[jtmp - 1] += tmp
                        if j != 1:
                            sum2 += tmp * delvac
                        else:
                            sum2 += tmp * 0.5 * delvac
                    cdesurf += sum2
                    for j in range(ns2):
                        tmp = (psisem[j]) ** 2 * cmult
                        if jsem[j] != 0:
                            cdesem[jsem[j] - 1] += tmp / nexsem[jsem[j] - 1]
                        else:
                            print('ERROR , INTCURR - J=0')
                            input()

                eperp = ener - wkparr ** 2 / C
                kappa = np.sqrt(C * (barr2[nbarr2 - 1] - eperp))
                trans = 2.0 * nwkdeg * (2.0 * wf) ** 2 * wkftip / (wksem / effm)
                sum1 += trans * occdiff

    currv = sum1 * dele * delwk ** 2 / (4.0 * PI ** 2 * RQUANT)

    # Localized States
    if pmax <= ev:
        return currv, currv0
    emax = pmax
    if nwk != 1:
        emin = max(ev, min(ef - 10.0 * tk1, ef + bias - 10.0 * tk2))
    else:
        emin = max(ev, ef - 10.0 * tk1)
    dele = (emax - emin) / nee
    emin2 = ef - 10.0 * tk1
    dele2 = (emax - emin2) / nee
    sum1 = 0.0
    if dele <= 0:
        return currv, currv0

    wkmax = np.sqrt(C * effm * (emax - emin))
    if wkmax == 0:
        return currv, currv0  # Avoid division by zero

    semstep = (2.0 * PI / wkmax) / expani
    kappa = np.sqrt(C * max(barr[nbarr1 - 1], barr[0]) - emin)
    if kappa == 0:
        return currv, currv0  # Avoid division by zero

    vacstep = (2.0 * PI / kappa) / expani
    pot0p = pot0 + vbprof[0]
    
    potexpand(impot, sep, nv, pot0p, s, nsp, nsdim, barr, nbarr1, barr2, nbarr2, nvdim1, nvdim2, vbprof, prof2, nsdim2, s2, ns2, vacstep, semstep, jsem, nexsem, nexvac, iwrit)

    delvac = sep / float(nbarr2 - 1)
    delwk = wkmax / nwk

    for iwky in range(nwk):
        wky = iwky * delwk
        for iwkx in range(nwk):
            wkx = iwkx * delwk
            wkparr = np.sqrt(wkx ** 2 + wky ** 2)
            eparr = wkparr ** 2 / (effm * C)
            n = 0
            nsav = 0
            nwkdeg = 8
            if iwkx == 0:
                nwkdeg /= 2
            if iwky == 0:
                nwkdeg /= 2
            if iwkx == iwky:
                nwkdeg /= 2
            if iwky > iwkx:
                continue
            for ie in range(1, nee + 1):
                ener = emax - (ie - 0.5) * dele
                if eparr >= (emax - ener):
                    continue
                occtip = fd(ener - bias, ef, tk2)
                occsem = fd(ener, ef, tk1)
                occdiff = occtip - occsem
                n, wf, wfderiv = vbloc(impot, n, 0.0, 0.0, ener, wkparr, sep, bias, barr2, nvdim2, nbarr2, prof2, ns2, nsdim2, s2, effm, ev, e2band, np.zeros(nvdim2), np.zeros(nsdim2))
                if n == nsav:
                    continue
                if iwkx == 0 and iwky == 0:
                    if psisem[0] != 0.0:
                        nloc[0] = nloc[0] + 1
                        if iwrit > 1:
                            print('VB localized state at energy ', ener)
                        if iwrit >= 1:
                            print(bias, n, ener, occdiff)
                        if iwrit >= 2:
                            for j in range(nbarr2 - 1, -1, -1):
                                print(-(j - 1) * delvac, psivac[j])
                            for j in range(ns2):
                                print(s2[j], psisem[j])

                        if icomp == 1:
                            occint = max(0.0, (ener - ef))
                            if tk1 != 0.0:
                                sum2 = 0.0
                                for ieocc in range(1, nee + 1):
                                    enerocc = ener - (ieocc - 0.5) * dele2
                                    if enerocc < emin2:
                                        break
                                    occsem2 = 1.0 - fd(enerocc, ef, tk1)
                                    sum2 += occsem2
                                sum2 *= dele2
                                occint = sum2
                            cmult = occint * (effm * C / (2.0 * PI))
                            sum2 = 0.0
                            for j in range(nbarr2 - 1, -1, -1):
                                tmp = (psivac[j]) ** 2 * cmult
                                if (j - 1) % nexvac == 0:
                                    jtmp = (j - 1) // nexvac + 1
                                    cdlvac[jtmp - 1] += tmp
                                if iwrit >= 2:
                                    print(-(j - 1) * delvac, tmp)
                                if j != 1:
                                    sum2 += tmp * delvac
                                else:
                                    sum2 += tmp * 0.5 * delvac
                            cdlsurf += sum2
                            for j in range(ns2):
                                tmp = (psisem[j]) ** 2 * cmult
                                if iwrit >= 2:
                                    print(s2[j], tmp)
                                if jsem[j] != 0:
                                    cdlsem[jsem[j] - 1] += tmp / nexsem[jsem[j] - 1]
                                else:
                                    print('ERROR , INTCURR - J=0')
                                    input()

                eperp = ener - wkparr ** 2 / C
                kappa = np.sqrt(C * (barr2[nbarr2 - 1] - eperp))
                trans = nwkdeg * (n - nsav) * 2.0 * (2.0 * wf) ** 2 * wkftip
                sum1 += trans * occdiff
                nsav = n

    currv0 = sum1 * delwk ** 2 / (C * 2.0 * PI * RQUANT)
    
    return currv, currv0


def vbwf(impot, wf, wfderiv, wksem, e, wkparr, sep, bias, barr2, nvdim2, nbarr2, prof2, ns2, nsdim2, s2, effm, ev, e2band, psivac, psisem):
    """Integrate valence band wavefunction from tip across to sample, return values of wf and derivative of wf at tip surface."""
    wf = 0.0
    wfderiv = 0.0
    wksem = 1.0
    eperp = e + wkparr ** 2 / (C * effm)
    if eperp >= ev:
        return wf, wfderiv, wksem

    # Determine initial conditions for wavefunction
    eperp = e - wkparr ** 2 / C
    psi = 1.0
    psivac[nbarr2 - 1] = psi
    imin = 1
    imax = nbarr2
    if eperp < barr2[0] and eperp < barr2[nbarr2 - 1]:
        return wf, wfderiv, wksem

    for i in range(nbarr2):
        if eperp < barr2[i]:
            imin = i
            break

    for i in range(nbarr2 - 1, -1, -1):
        psivac[i] = psi
        if eperp < barr2[i]:
            imax = i
            break

    if imax > imin:
        return wf, wfderiv, wksem

    dpsi = psi * np.sqrt(C * (barr2[imax - 1] - eperp))
    wf = psi
    wfderiv = dpsi

    # Integrate through vacuum
    delvac = sep / float(nbarr2 - 1)
    for i in range(imax - 2, -1, -1):
        if impot == 1 and (barr2[i] - eperp) <= 0:
            continue
        psi = psi + dpsi * delvac
        psivac[i] = psi
        dpsi = dpsi + C * (barr2[i] - eperp) * psi * delvac

    # Match across vacuum-semiconductor interface
    psi = psi
    dpsi = dpsi * effm

    # Integrate through semiconductor
    eperp = e + wkparr ** 2 / (C * effm)
    psi = psi + dpsi * s2[0]
    psisem[0] = psi
    ebarr = eperp - prof2[0]
    dpsi = dpsi + C * effm * ebarr * psi * s2[0]
    for i in range(1, ns2):
        dels = s2[i] - s2[i - 1]
        psi = psi + dpsi * dels
        if abs(psi) >= 1e100:
            return 0.0, 0.0, wksem
        psisem[i] = psi
        ebarr = eperp - prof2[i]
        dpsi = dpsi + C * effm * ebarr * psi * dels

    # Determine amplitude, and evaluate across barrier
    wksem = np.sqrt(C * effm * (ev - eperp))
    phase = np.arctan(psi * wksem / dpsi)
    amp = np.sqrt(2.0) * np.sin(phase) / psi
    wf = wf * amp
    wfderiv = wfderiv * amp

    # Normalize wavefunction
    for i in range(nbarr2 - 1, -1, -1):
        psivac[i] = psivac[i] * amp
    for i in range(ns2):
        psisem[i] = psisem[i] * amp

    delsmax = s2[ns2 - 1] - s2[ns2 - 2]
    if delsmax / (2.0 * PI / wksem) > 0.25:
        print("*** CAUTION *** RATIO OF SEMICOND. STEP SIZE TO WAVELENGTH =", delsmax / (2.0 * PI / wksem))

    return wf, wfderiv, wksem

def cbwf(impot, wf, wfderiv, wksem, e, wkparr, sep, bias, barr2, nvdim2, nbarr2, prof2, ns2, nsdim2, s2, effm, ec, e2band, psivac, psisem):
    """Integrate conduction band wavefunction from tip across to sample, return values of wf and derivative of wf at tip surface."""
    c = 26.254
    pi = 4.0 * np.arctan(1.0)

    wf = 0.0
    wfderiv = 0.0
    wksem = 1.0
    eperp = e - wkparr ** 2 / (c * effm)
    if eperp <= ec:
        return wf, wfderiv, wksem

    # Determine initial conditions for wavefunction
    eperp = e - wkparr ** 2 / c
    psi = 1.0
    psivac[nbarr2 - 1] = psi
    imin = 1
    imax = nbarr2
    if eperp < barr2[0] and eperp < barr2[nbarr2 - 1]:
        return wf, wfderiv, wksem

    for i in range(nbarr2):
        if eperp < barr2[i]:
            imin = i
            break

    for i in range(nbarr2 - 1, -1, -1):
        psivac[i] = psi
        if eperp < barr2[i]:
            imax = i
            break

    if imax > imin:
        return wf, wfderiv, wksem

    dpsi = psi * np.sqrt(c * (barr2[imax] - eperp))
    wf = psi
    wfderiv = dpsi

    # Integrate through vacuum
    delvac = sep / float(nbarr2 - 1)
    for i in range(imax - 1, 0, -1):
        if impot == 1 and (barr2[i] - eperp) <= 0:
            continue
        psi = psi + dpsi * delvac
        psivac[i] = psi
        dpsi = dpsi + c * (barr2[i] - eperp) * psi * delvac

    # Match across vacuum-semiconductor interface
    psi = psi
    dpsi = dpsi * effm

    # Integrate through semiconductor
    eperp = e - wkparr ** 2 / (c * effm)
    psi = psi + dpsi * s2[0]
    psisem[0] = psi
    ebarr = prof2[0] - eperp
    dpsi = dpsi + c * effm * ebarr * psi * s2[0]
    for i in range(1, ns2):
        dels = s2[i] - s2[i - 1]
        psi = psi + dpsi * dels
        if abs(psi) >= 1e100:
            return 0.0, 0.0, wksem
        psisem[i] = psi
        ebarr = prof2[i] - eperp
        dpsi = dpsi + c * effm * ebarr * psi * dels

    # Determine amplitude, and evaluate across barrier
    wksem = np.sqrt(c * effm * (eperp - ec))
    phase = np.arctan(psi * wksem / dpsi)
    amp = np.sqrt(2.0) * np.sin(phase) / psi
    wf = wf * amp
    wfderiv = wfderiv * amp

    # Normalize wavefunction
    for i in range(nbarr2 - 1, -1, -1):
        psivac[i] = psivac[i] * amp
    for i in range(ns2):
        psisem[i] = psisem[i] * amp

    delsmax = s2[ns2 - 1] - s2[ns2 - 2]
    if delsmax / (2.0 * pi / wksem) > 0.25:
        print("*** CAUTION *** RATIO OF SEMICOND. STEP SIZE TO WAVELENGTH =", delsmax / (2.0 * pi / wksem))

    return wf, wfderiv, wksem

def vbloc(impot, nsign, wf, wfderiv, e, wkparr, sep, bias, barr2, nvdim2, nbarr2, prof2, ns2, nsdim2, s2, effm, ev, e2band, psivac, psisem):
    """Integrate localized VB wavefunction from tip across to sample, looking for switches in sign of wf to enumerate localized states."""
    c = 26.254

    wf = 0.0
    wfderiv = 0.0
    isav = ns2
    eperp = e + wkparr ** 2 / (c * effm)
    if eperp <= ev:
        return nsign, wf, wfderiv

    nsign = 0

    # Determine initial conditions for wavefunction
    eperp = e - wkparr ** 2 / c
    sum1 = 0.0
    sum2 = 0.0
    psi = 1.0
    psivac[nbarr2 - 1] = psi
    imin = 1
    imax = nbarr2
    if eperp < barr2[0] and eperp < barr2[nbarr2 - 1]:
        return nsign, wf, wfderiv

    for i in range(nbarr2):
        if eperp < barr2[i]:
            imin = i
            break

    for i in range(nbarr2 - 1, -1, -1):
        psivac[i] = psi
        if eperp < barr2[i]:
            imax = i
            break

    if imax > imin:
        return nsign, wf, wfderiv

    dpsi = psi * np.sqrt(c * (barr2[imax] - eperp))
    wf = psi
    wfderiv = dpsi

    # Integrate through vacuum
    delvac = sep / float(nbarr2 - 1)
    for i in range(imax - 1, 0, -1):
        if impot == 1 and (barr2[i] - eperp) <= 0:
            continue
        psi = psi + dpsi * delvac
        psivac[i] = psi
        sum1 += psi ** 2 * delvac
        dpsi = dpsi + c * (barr2[i] - eperp) * psi * delvac

    # Match across vacuum-semiconductor interface
    psi = psi
    dpsi = dpsi * effm

    # Integrate through semiconductor
    eperp = e + wkparr ** 2 / (c * effm)
    psi = psi + dpsi * s2[0]
    psisem[0] = psi
    sum1 += psi ** 2 * s2[0]
    ebarr = eperp - prof2[0]
    dpsi = dpsi + c * effm * ebarr * psi * s2[0]
    for i in range(1, ns2 - 1):
        dels = s2[i] - s2[i - 1]
        psisav = psi
        psi = psi + dpsi * dels
        psisem[i] = psi
        if psisav * psi < 0.0:
            nsign += 1
            isav = i
            sum2 += sum1
            sum1 = 0.0
        sum1 += psi ** 2 * dels
        ebarr = eperp - prof2[i]
        dpsi = dpsi + c * effm * ebarr * psi * dels

    # Normalize wavefunction
    if sum2 != 0.0:
        amp = 1.0 / np.sqrt(sum2)
    else:
        amp = 1.0 / np.sqrt(sum1)
    wf = wf * amp
    wfderiv = wfderiv * amp
    for i in range(nbarr2 - 1, -1, -1):
        psivac[i] = psivac[i] * amp
    for i in range(isav):
        psisem[i] = psisem[i] * amp
    for i in range(isav, ns2):
        psisem[i] = 0.0

    return nsign, wf, wfderiv

def cbloc(impot, nsign, wf, wfderiv, e, wkparr, sep, bias, barr2, nvdim2, nbarr2, prof2, ns2, nsdim2, s2, effm, ec, e2band, psivac, psisem):
    """Integrate localized CB wavefunction from tip across to sample, looking for switches in sign of wf to enumerate localized states."""
    c = 26.254

    wf = 0.0
    wfderiv = 0.0
    isav = ns2
    eperp = e - wkparr ** 2 / (c * effm)
    if eperp >= ec:
        return nsign, wf, wfderiv

    nsign = 0

    # Determine initial conditions for wavefunction
    eperp = e - wkparr ** 2 / c
    sum1 = 0.0
    sum2 = 0.0
    psi = 1.0
    psivac[nbarr2 - 1] = psi
    imin = 1
    imax = nbarr2
    if eperp < barr2[0] and eperp < barr2[nbarr2 - 1]:
        return nsign, wf, wfderiv

    for i in range(nbarr2):
        if eperp < barr2[i]:
            imin = i
            break

    for i in range(nbarr2 - 1, -1, -1):
        psivac[i] = psi
        if eperp < barr2[i]:
            imax = i
            break

    if imax > imin:
        return nsign, wf, wfderiv

    dpsi = psi * np.sqrt(c * (barr2[imax] - eperp))
    wf = psi
    wfderiv = dpsi

    # Integrate through vacuum
    delvac = sep / float(nbarr2 - 1)
    for i in range(imax - 1, 0, -1):
        if impot == 1 and (barr2[i] - eperp) <= 0:
            continue
        psi = psi + dpsi * delvac
        psivac[i] = psi
        sum1 += psi ** 2 * delvac
        dpsi = dpsi + c * (barr2[i] - eperp) * psi * delvac

    # Match across vacuum-semiconductor interface
    psi = psi
    dpsi = dpsi * effm

    # Integrate through semiconductor
    eperp = e - wkparr ** 2 / (c * effm)
    psi = psi + dpsi * s2[0]
    psisem[0] = psi
    sum1 += psi ** 2 * s2[0]
    ebarr = prof2[0] - eperp
    dpsi = dpsi + c * effm * ebarr * psi * s2[0]
    for i in range(1, ns2 - 1):
        dels = s2[i] - s2[i - 1]
        psisav = psi
        psi = psi + dpsi * dels
        psisem[i] = psi
        if psisav * psi < 0.0:
            nsign += 1
            isav = i
            sum2 += sum1
            sum1 = 0.0
        sum1 += psi ** 2 * dels
        ebarr = prof2[i] - eperp
        dpsi = dpsi + c * effm * ebarr * psi * dels

    # Normalize wavefunction
    if sum2 != 0.0:
        amp = 1.0 / np.sqrt(sum2)
    else:
        amp = 1.0 / np.sqrt(sum1)
    wf = wf * amp
    wfderiv = wfderiv * amp
    for i in range(nbarr2 - 1, -1, -1):
        psivac[i] = psivac[i] * amp
    for i in range(isav):
        psisem[i] = psisem[i] * amp
    for i in range(isav, ns2):
        psisem[i] = 0.0

    return nsign, wf, wfderiv

def cbcurr(impot, barr, nbarr1, s, prof, nsp, nsdim, nvdim, nv, nvdim1, nsdim2, nvdim2, sep, pot0, expani, pmin, effm, tk1, tk2, bias, ef, nee, nwk, nloc, ec, eftip, e2band, lun, cdesem, cdesurf, cdlsem, cdlsurf, cdevac, cdlvac, currc, currc0, iwrit, icomp, cbprof):
    """CB Tunnel Current"""
    nwk=1
    nee=1
    expani=1
    kappa = 0.0
    wkftip = np.sqrt(C * eftip)
    sum1 = sum2 = sum1s = sum2s = sum1p = sum2p = 0.0
    currc = 0.0
    currc0 = 0.0
    emin = ec
    if nwk != 1:
        emax = max(ef + 10.0 * tk1, ef + bias + 10.0 * tk2)
    else:
        emax = ef + 10.0 * tk1
    dele = (emax - emin) / nee
    if dele <= 0:
        return currc, currc0

    wkmax = np.sqrt(C * effm * (emax - emin))
    if wkmax == 0:
        return currc, currc0  # Avoid division by zero

    semstep = (2.0 * PI / wkmax) / expani
    kappa = np.sqrt(C * max(barr[nbarr1 - 1], barr[0]) - emin)
    if kappa == 0:
        return currc, currc0  # Avoid division by zero

    vacstep = (2.0 * PI / kappa) / expani
    pot0p = pot0 + cbprof[0]
    
    potexpand(impot, sep, nv, pot0p, s, nsp, nsdim, barr, nbarr1, barr2, nbarr2, nvdim1, nvdim2, cbprof, prof2, nsdim2, s2, ns2, vacstep, semstep, jsem, nexsem, nexvac, iwrit)
    
    delvac = sep / float(nbarr2 - 1)
    delwk = wkmax / nwk

    for iwky in range(nwk):
        wky = iwky * delwk
        for iwkx in range(nwk):
            wkx = iwkx * delwk
            wkparr = np.sqrt(wkx ** 2 + wky ** 2)
            eparr = wkparr ** 2 / (effm * C)
            nwkdeg = 8
            if iwkx == 0:
                nwkdeg /= 2
            if iwky == 0:
                nwkdeg /= 2
            if iwkx == iwky:
                nwkdeg /= 2
            if iwky > iwkx:
                continue
            for ie in range(1, nee + 1):
                ener = emin + (ie - 0.5) * dele
                if eparr >= (ener - emin):
                    continue
                occtip = fd(ener - bias, ef, tk2)
                occsem = fd(ener, ef, tk1)
                occdiff = occtip - occsem
                wf, wfderiv, wksem = cbwf(impot, 0.0, 0.0, 0.0, ener, wkparr, sep, bias, barr2, nvdim2, nbarr2, prof2, ns2, nsdim2, s2, effm, ec, e2band, np.zeros(nvdim2), np.zeros(nsdim2))
                if iwkx == 0 and iwky == 0 and iwrit >= 4:
                    for j in range(nbarr2 - 1, -1, -1):
                        print(-(j - 1) * delvac, psivac[j])
                    for j in range(ns2):
                        print(s2[j], psisem[j])

                if icomp == 1:
                    if ie == 1:
                        sum3 = 0.0
                        for ieocc in range(ie, nee + 1):
                            enerocc = emin + ieocc * dele
                            occsem2 = fd(enerocc, ef, tk1)
                            sum3 += occsem2
                    else:
                        enerocc = emin + (ie - 1) * dele
                        occsem2 = fd(enerocc, ef, tk1)
                        sum3 -= occsem2
                    occint = -sum3 * dele
                    cmult = occint * (effm * C / (2.0 * PI)) * (dele * effm * C / (2.0 * PI * wksem))
                    sum2 = 0.0
                    for j in range(nbarr2 - 1, -1, -1):
                        tmp = (psivac[j]) ** 2 * cmult
                        if (j - 1) % nexvac == 0:
                            jtmp = (j - 1) // nexvac + 1
                            cdevac[jtmp - 1] += tmp
                        if j != 1:
                            sum2 += tmp * delvac
                        else:
                            sum2 += tmp * 0.5 * delvac
                    cdesurf += sum2
                    for j in range(ns2):
                        tmp = (psisem[j]) ** 2 * cmult
                        if jsem[j] != 0:
                            cdesem[jsem[j] - 1] += tmp / nexsem[jsem[j] - 1]
                        else:
                            print('ERROR , INTCURR - J=0')
                            input()

                eperp = ener - wkparr ** 2 / C
                kappa = np.sqrt(C * (barr2[nbarr2 - 1] - eperp))
                trans = 2.0 * nwkdeg * (2.0 * wf) ** 2 * wkftip / (wksem / effm)
                sum1 += trans * occdiff

    currc = sum1 * dele * delwk ** 2 / (4.0 * PI ** 2 * RQUANT)

    # Localized States
    if pmin >= ec:
        return currc, currc0
    emin = pmin
    if nwk != 1:
        emax = min(ec, max(ef + 10.0 * tk1, ef + bias + 10.0 * tk2))
    else:
        emax = min(ec, ef + 10.0 * tk1)
    dele = (emax - emin) / nee
    emax2 = ef + 10.0 * tk1
    dele2 = (emax2 - emin) / nee
    sum1 = 0.0
    if dele <= 0:
        return currc, currc0

    wkmax = np.sqrt(C * effm * (emax - emin))
    if wkmax == 0:
        return currc, currc0  # Avoid division by zero

    semstep = (2.0 * PI / wkmax) / expani
    kappa = np.sqrt(C * max(barr[nbarr1 - 1], barr[0]) - emin)
    if kappa == 0:
        return currc, currc0  # Avoid division by zero

    vacstep = (2.0 * PI / kappa) / expani
    pot0p = pot0 + cbprof[0]
    
    potexpand(impot, sep, nv, pot0p, s, nsp, nsdim, barr, nbarr1, barr2, nbarr2, nvdim1, nvdim2, cbprof, prof2, nsdim2, s2, ns2, vacstep, semstep, jsem, nexsem, nexvac, iwrit)

    delvac = sep / float(nbarr2 - 1)
    delwk = wkmax / nwk

    for iwky in range(nwk):
        wky = iwky * delwk
        for iwkx in range(nwk):
            wkx = iwkx * delwk
            wkparr = np.sqrt(wkx ** 2 + wky ** 2)
            eparr = wkparr ** 2 / (effm * C)
            n = 0
            nsav = 0
            nwkdeg = 8
            if iwkx == 0:
                nwkdeg /= 2
            if iwky == 0:
                nwkdeg /= 2
            if iwkx == iwky:
                nwkdeg /= 2
            if iwky > iwkx:
                continue
            for ie in range(1, nee + 1):
                ener = emin + (ie - 0.5) * dele
                if eparr >= (ener - emin):
                    continue
                occtip = fd(ener - bias, ef, tk2)
                occsem = fd(ener, ef, tk1)
                occdiff = occtip - occsem
                n, wf, wfderiv = cbloc(impot, n, 0.0, 0.0, ener, wkparr, sep, bias, barr2, nvdim2, nbarr2, prof2, ns2, nsdim2, s2, effm, ec, e2band, np.zeros(nvdim2), np.zeros(nsdim2))
                if n == nsav:
                    continue
                if iwkx == 0 and iwky == 0:
                    if psisem[0] != 0.0:
                        nloc[1] = nloc[1] + 1
                        if iwrit > 1:
                            print('CB localized state at energy ', ener)
                        if iwrit >= 1:
                            print(bias, n, ener, occdiff)
                        if iwrit >= 2:
                            for j in range(nbarr2 - 1, -1, -1):
                                print(-(j - 1) * delvac, psivac[j])
                            for j in range(ns2):
                                print(s2[j], psisem[j])

                        if icomp == 1:
                            occint = -max(0.0, (ef - ener))
                            if tk1 != 0.0:
                                sum2 = 0.0
                                for ieocc in range(1, nee + 1):
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
                                tmp = (psivac[j]) ** 2 * cmult
                                if (j - 1) % nexvac == 0:
                                    jtmp = (j - 1) // nexvac + 1
                                    cdlvac[jtmp - 1] += tmp
                                if iwrit >= 2:
                                    print(-(j - 1) * delvac, tmp)
                                if j != 1:
                                    sum2 += tmp * delvac
                                else:
                                    sum2 += tmp * 0.5 * delvac
                            cdlsurf += sum2
                            for j in range(ns2):
                                tmp = (psisem[j]) ** 2 * cmult
                                if iwrit >= 2:
                                    print(s2[j], tmp)
                                if jsem[j] != 0:
                                    cdlsem[jsem[j] - 1] += tmp / nexsem[jsem[j] - 1]
                                else:
                                    print('ERROR , INTCURR - J=0')
                                    input()

                eperp = ener - wkparr ** 2 / C
                kappa = np.sqrt(C * (barr2[nbarr2 - 1] - eperp))
                trans = nwkdeg * (n - nsav) * 2.0 * (2.0 * wf) ** 2 * wkftip
                sum1 += trans * occdiff
                nsav = n

    currc0 = sum1 * delwk ** 2 / (C * 2.0 * PI * RQUANT)
    
    return currc, currc0


# 主程式
def intcurr(impot, barr, prof, nbarr1, nv, ns, nsp, nvdim, nsdim, s, sep, bias, ef, chi, 
             eftip, cpot, egap, tk, avbh, avbl, avbso, acb, eso, nee, nwk, expani ,pot0, nvdim1, nvdim2, nsdim2, nloc, currv, currv0, currc, currc0, 
             curr, curr0, iwrit, icomp, cdesem, cdesurf, cdlsem, cdlsurf, cdevac, cdlvac):
    nwk=1
    nee=1
    expani=1     
    cdesem = np.zeros(nsdim)
    cdlsem = np.zeros(nsdim)
    cdevac = np.zeros(nvdim1)
    cdlvac = np.zeros(nvdim1)
    vbprof = np.zeros(nsdim)
    cbprof = np.zeros(nsdim)
    nexsem = np.zeros(nvdim2, dtype=int)

    rquant = 12900.0
    pi = 4.0 * np.arctan(1.0)
    tk1 = tk
    tk2 = tk

    if icomp == 1:
        cdesurf = 0.0
        cdlsurf = 0.0
        cdesem.fill(0.0)
        cdlsem.fill(0.0)
        cdevac.fill(0.0)
        cdlvac.fill(0.0)

    # Initialize pmax to a very small value
    pmax = -np.inf

    # Valence Band
    for j in range(nsp):
        sz = s[j]
        vbprof[j] = prof[j] + vbedge(sz)
        if j == 0:
            pmax = vbprof[j]
        else:
            pmax = max(pmax, vbprof[j])

    ev = prof[ns - 1] + vbedge(s[ns - 1])
    if iwrit == 4:
        print(f'MAXIMUM ENERGY IN VB PROFILE = {pmax}')
        print(f'VB EDGE ENERGY DEEP INSIDE SEMICONDUCTOR = {ev}')
    if iwrit >= 1:
        for j in range(nsp):
            print(s[j], vbprof[j])

    lun = 30
    currvl, currv0l = vbcurr(impot, barr, nbarr1, s, prof, nsp, nsdim, nvdim, nv, nvdim1, nsdim2, nvdim2, sep, pot0, expani, pmax, avbl, tk1, tk2, bias, ef, nee, nwk, nloc, ev, eftip, egap, lun, cdesem, cdesurf, cdlsem, cdlsurf, cdevac, cdlvac, currv, currv0, iwrit, icomp, vbprof)
    if iwrit != 0:
        print(f'number of VB light-hole localized states = {nloc}')

    lun = 40
    currvh, currv0h = vbcurr(impot, barr, nbarr1, s, prof, nsp, nsdim, nvdim, nv, nvdim1, nsdim2, nvdim2, sep, pot0, expani, pmax, avbh, tk1, tk2, bias, ef, nee, nwk, nloc, ev, eftip, e2hh, lun, cdesem, cdesurf, cdlsem, cdlsurf, cdevac, cdlvac, currv, currv0, iwrit, icomp, vbprof)
    if iwrit != 0:
        print(f'number of VB heavy-hole localized states = {nloc}')

    evso = ev - eso
    pmaxso = pmax - eso
    for j in range(nsp):
        vbprof[j] = vbprof[j] - eso

    lun = 50
    currvso, currv0so = vbcurr(impot, barr, nbarr1, s, prof, nsp, nsdim, nvdim, nv, nvdim1, nsdim2, nvdim2, sep, pot0, expani, pmaxso, avbso, tk1, tk2, bias, ef, nee, nwk, nloc, evso, eftip, e2so, lun, cdesem, cdesurf, cdlsem, cdlsurf, cdevac, cdlvac, currv, currv0, iwrit, icomp, vbprof)
    if iwrit != 0:
        print(f'number of VB split-off localized states = {nloc}')

    currv = currvl + currvh + currvso
    currv0 = currv0l + currv0h + currv0so

    # Initialize pmin to a very large value
    pmin = np.inf

    # Conduction Band
    for j in range(nsp):
        sz = s[j]
        cbprof[j] = prof[j] + cbedge(sz)
        if j == 0:
            pmin = cbprof[j]
        else:
            pmin = min(pmin, cbprof[j])

    ec = prof[ns - 1] + cbedge(s[ns - 1])
    if iwrit == 4:
        print(f'MINIMUM ENERGY IN CB PROFILE = {pmin}')
        print(f'CB EDGE ENERGY DEEP INSIDE SEMICONDUCTOR = {ec}')
    if iwrit >= 1:
        for j in range(nsp):
            print(s[j], cbprof[j])

    lun = 60
    currc, currc0 = cbcurr(impot, barr, nbarr1, s, prof, nsp, nsdim, nvdim, nv, nvdim1, nsdim2, nvdim2, sep, pot0, expani, pmin, acb, tk1, tk2, bias, ef, nee, nwk, nloc, ec, eftip, egap, lun, cdesem, cdesurf, cdlsem, cdlsurf, cdevac, cdlvac, currc, currc0, iwrit, icomp, cbprof)
    if iwrit != 0:
        print(f'number of CB localized states = {nloc}')

    curr = currv + currc
    curr0 = currv0 + currc0

    return currv, currv0, currc, currc0, curr, curr0

def cbedge(sz):
    a1 = 1.0  # 根據擬合或計算得到的線性項係數
    a2 = 0.1  # 根據擬合或計算得到的二次項係數

    """
    Calculate the conduction band edge.
    
    Parameters:
    sz (float): Input parameter.
    a1 (float): Coefficient for the linear term.
    a2 (float): Coefficient for the quadratic term.
    
    Returns:
    float: Conduction band edge.
    """
    return a1 * sz + a2 * sz ** 2
def vbedge(sz):
    b1 = 1.0  # 根據擬合或計算得到的線性項係數
    b2 = 0.1  # 根據擬合或計算得到的二次項係數

    """
    Calculate the valence band edge.
    
    Parameters:
    sz (float): Input parameter.
    b1 (float): Coefficient for the linear term.
    b2 (float): Coefficient for the quadratic term.
    
    Returns:
    float: Valence band edge.
    """
    return b1 * sz + b2 * sz ** 2
