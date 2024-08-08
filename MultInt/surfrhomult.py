import numpy as np
from semirhomult import fd


def SURFRHO(IAR, DELE, ESTART, NE, NEDIM, RHOSTAB, ISTK, TK, EN0, EN, DENS, FWHM, ECENT, surf):
    if NE > NEDIM:
        raise ValueError("*** ERROR - NE > NEDIM; PROGRAM HALTED")

    if ISTK == 1:
        for I in range(1, NE + 1):
            EF1 = (I - 1) * DELE + ESTART
            if DENS[IAR, 1] == 0:
                RHOSTAB[IAR, I - 1] = rhos1(IAR, EF1, DELE, surf)
            elif DENS[IAR, 0] == 0:
                RHOSTAB[IAR, I - 1] = rhos2(IAR, EF1, DELE, surf)
            else:
                RHOSTAB[IAR, I - 1] = rhos(IAR, EF1, DELE, surf)
    else:
        if DENS[IAR, 0] == 0 or DENS[IAR, 1] == 0:
            NEN = int(np.round((EN0[IAR] - ESTART) / DELE)) + 1
            RHOSTAB[IAR, NEN - 1] = 0
            SUM = 0
            for I in range(NEN, NE):
                EF1 = (I - 1) * DELE + ESTART
                SUM += sigsum(IAR, EF1)
                RHOSTAB[IAR, I] = SUM * DELE
            SUM = 0
            for I in range(NEN - 2, -1, -1):
                EF1 = (I - 1) * DELE + ESTART
                SUM += sigsum(IAR, EF1)
                RHOSTAB[IAR, I] = SUM * DELE
        else:
            NEN = int(np.round((EN0[IAR] - ESTART) / DELE)) + 1
            RHOSTAB[IAR, NEN - 1] = 0
            for I in range(NEN, NE):
                EF1 = (I - 1) * DELE + ESTART
                RHOSTAB[IAR, I] = rhos(IAR, EF1, DELE, surf)
            for I in range(NEN - 2, -1, -1):
                EF1 = (I - 1) * DELE + ESTART
                RHOSTAB[IAR, I] = rhos(IAR, EF1, DELE, surf)
    return RHOSTAB


def sig(IAR, ID, ENER, surf):
    PI = np.pi
    """
    Parameters:
    IAR: Index of the area
    ID: Index of the state (1 or 2)
    ENER: Energy value
    surf: Dictionary containing surface state parameters

    Returns:
    Surface state spectrum value
    """
    fwhm = surf["FWHM"]
    EN = surf["EN"]
    dens = surf["DENS"]
    ecent = surf["ECENT"]

    if fwhm[IAR, ID] == 0.0:
        sig_value = dens[IAR, ID]
        if ENER > EN[IAR, ID]:
            sig_value = -sig_value
        return sig_value

    width = fwhm[IAR, ID] / (2.0 * np.sqrt(2.0 * np.log(2.0)))
    sig_value = -np.exp(-((ENER - (EN[IAR, ID] + ecent[IAR, ID])) ** 2) /
                        (2.0 * width ** 2)) + np.exp(-((ENER - (EN[IAR, ID] - ecent[IAR, ID])) ** 2) / (2.0 * width ** 2))
    sig_value = sig_value * dens[IAR, ID] / (np.sqrt(2.0 * PI) * width)
    return sig_value


def sigsum(iar, ener):
    return sig(iar, 0, ener) + sig(iar, 1, ener)


def rhos1(iar, ef1, dele, surf):
    sum_val = 0.0
    e = surf["EN"][iar, 0]

    if ef1 == surf["EN"][iar, 0]:
        return sum_val

    if surf["ISTK"] == 0:
        if ef1 < surf["EN"][iar, 0]:
            while e < ef1:
                e += dele
                sum_val += sig(iar, 0, e, surf) * dele
        else:
            while e > ef1:
                e -= dele
                sum_val += sig(iar, 0, e, surf) * dele
        return sum_val

    if ef1 < surf["EN"][iar, 0]:
        while e < ef1 + 10.0 * surf["TK"]:
            e += dele
            sum_val += sig(iar, 0, e, surf) * fd(e, ef1, surf["TK"]) * dele
    else:
        while e > ef1 - 10.0 * surf["TK"]:
            e -= dele
            sum_val += sig(iar, 0, e, surf) * \
                (1.0 - fd(e, ef1, surf["TK"])) * dele

    return sum_val


def rhos2(iar, ef1, dele, surf):
    sum_val = 0.0
    e = surf["EN"][iar, 1]

    if ef1 == surf["EN"][iar, 1]:
        return sum_val

    if surf["ISTK"] == 0:
        if ef1 < surf["EN"][iar, 1]:
            while e < ef1:
                e += dele
                sum_val += sig(iar, 1, e, surf) * dele
        else:
            while e > ef1:
                e -= dele
                sum_val += sig(iar, 1, e, surf) * dele
        return sum_val

    if ef1 < surf["EN"][iar, 1]:
        while e < ef1 + 10.0 * surf["TK"]:
            e += dele
            sum_val += sig(iar, 1, e, surf) * fd(e, ef1, surf["TK"]) * dele
    else:
        while e > ef1 - 10.0 * surf["TK"]:
            e -= dele
            sum_val += sig(iar, 1, e, surf) * \
                (1.0 - fd(e, ef1, surf["TK"])) * dele

    return sum_val


def rhos(iar, ef1, dele, surf):
    return rhos1(iar, ef1, dele, surf) + rhos2(iar, ef1, dele, surf)


def enfind(iar, en1, en2, ne, surf):
    en0 = en1
    estart = en1
    dele = abs(en1 - en2) / float(ne)
    if dele == 0:
        return en0
    if en2 < en1:
        dele = -dele

    for ie in range(ne + 1):
        ef1 = estart + ie * dele
        sigtmp = rhos(iar, ef1, abs(dele), surf)
        if dele > 0 and sigtmp <= 0:
            en0 = ef1
            break
        elif dele < 0 and sigtmp >= 0:
            en0 = ef1
            break

    return en0
