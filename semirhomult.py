import numpy as np


def fjint(J, ETA):
    PI = 3.141592654
    if ETA > 40:
        return np.sqrt(ETA**(J+2)) / (J/2.0 + 1)
    if ETA < -8:
        result = np.sqrt(PI) * np.exp(max(-40., ETA)) / 2.0
        if J == 3:
            result *= 3.0 / 2.0
        return result

    def fj(X):
        return np.sqrt(X**J) / (1.0 + np.exp(X - ETA))

    def trap(F, XMIN, XMAX, NSTEP):
        DELX = (XMAX - XMIN) / float(NSTEP)
        SUM = (F(XMIN) + F(XMAX)) / 2.0
        if NSTEP >= 2:
            for I in range(1, NSTEP):
                X = XMIN + I * DELX
                SUM += F(X)
        return SUM * DELX

    return trap(fj, 0.0, 20.0 + ETA, 1000)


def fd(e, ef, tk):
    if tk == 0:
        if e == ef:
            return 0.5
        if e < ef:
            return 1.0
        return 0.0
    ratio = (e - ef) / tk
    if ratio > 40:
        return 0.0
    if ratio < -40:
        return 1.0
    return 1.0 / (1.0 + np.exp(ratio))


def rhocb(IREG, EF, Pot, semi, fjint_function):
    C = 6.815e21
    if semi["IINV"][IREG] == 2 or semi["IINV"][IREG] == 3:
        return 0.0
    if semi["TK"] != 0.0:
        return C * np.sqrt((semi["ACB"][IREG] * semi["TK"])**3) * fjint_function(1, (EF - semi["EGAP"][IREG] - semi["DELVB"][IREG] - Pot) / semi["TK"])
    if (EF - semi["EGAP"][IREG] - semi["DELVB"][IREG] - Pot) <= 0.0:
        return 0.0
    return (2.0 * C / 3.0) * np.sqrt((semi["ACB"][IREG] * (EF - semi["EGAP"][IREG] - semi["DELVB"][IREG] - Pot))**3)


def rhovb(IREG, EF, Pot, semi, fjint_function):
    C = 6.815e21
    if semi["IINV"][IREG] == 1 or semi["IINV"][IREG] == 3:
        return 0.0
    if semi["TK"] != 0.0:
        return C * np.sqrt((semi["AVB"][IREG] * semi["TK"])**3) * fjint_function(1, (-EF + semi["DELVB"][IREG] + Pot) / semi["TK"])
    if (-EF + semi["DELVB"][IREG] + Pot) <= 0.0:
        return 0.0
    return (2.0 * C / 3.0) * np.sqrt((semi["AVB"][IREG] * (-EF + semi["DELVB"][IREG] + Pot))**3)


def rhod(IREG, EF, Pot, semi):
    RHOD = semi["CD"][IREG]
    if semi["IDEG"][IREG] == 1:
        return RHOD
    EXPO = EF - semi["EGAP"][IREG] - \
        semi["DELVB"][IREG] + semi["ED"][IREG] - Pot
    if semi["TK"] != 0.0:
        EXPO /= semi["TK"]
        if EXPO < -40:
            return semi["CD"][IREG]
        if EXPO > 40:
            return 0.0
        return semi["CD"][IREG] / (1.0 + 2.0 * np.exp(EXPO))
    if EXPO <= 0:
        return semi["CD"][IREG]
    return 0.0


def rhoa(IREG, EF, Pot, semi):
    RHOA = semi["CA"][IREG]
    if semi["IDEG"][IREG] == 1:
        return RHOA
    EXPO = semi["EA"][IREG] - EF + semi["DELVB"][IREG] + Pot
    if semi["TK"] != 0.0:
        EXPO /= semi["TK"]
        if EXPO < -40:
            return semi["CA"][IREG]
        if EXPO > 40:
            return 0.0
        return semi["CA"][IREG] / (1.0 + 4.0 * np.exp(EXPO))
    if EXPO <= 0:
        return semi["CA"][IREG]
    return 0.0


def rhob(IREG, EF, Pot, rhocb_function, rhoa_function, rhovb_function, rhod_function, semi, fjint_function):
    RHOE = rhocb_function(IREG, EF, Pot, semi, fjint_function) + \
        rhoa_function(IREG, EF, Pot, semi)
    RHOH = rhovb_function(IREG, EF, Pot, semi, fjint_function) + \
        rhod_function(IREG, EF, Pot, semi)
    return -RHOE + RHOH


def arho(EF, IREG, rhob_function, rhocb_function, rhoa_function, rhovb_function, rhod_function, semi, fjint_function):
    return abs(rhob_function(IREG, EF, 0, rhocb_function, rhoa_function, rhovb_function, rhod_function, semi, fjint_function))


def effind(IREG, semi, arho_function, gsect_function, rhob_function,
           rhocb_function, rhoa_function, rhovb_function, rhod_function, fjint_function):
    IINVSAV = semi["IINV"][IREG]
    semi["IINV"][IREG] = 0

    if semi["CD"][IREG] == 0 and semi["CA"][IREG] == 0:
        EF = semi["EGAP"][IREG] / 2 + 0.75 * semi["TK"] * \
            np.log(semi["AVB"][IREG] / semi["ACB"][IREG])
    elif semi["TK"] == 0:
        if semi["CA"][IREG] > semi["CD"][IREG]:
            EF = semi["EA"][IREG] / 2
        elif semi["CA"][IREG] < semi["CD"][IREG]:
            EF = semi["EGAP"][IREG] - semi["ED"][IREG] / 2
        else:
            EF = (semi["EGAP"][IREG] - semi["ED"][IREG] + semi["EA"][IREG]) / 2
    else:
        ESTART = -0.1 + semi["DELVB"][IREG]
        NE = 1000
        DELE = (semi["EGAP"][IREG] + semi["DELVB"][IREG] + 0.2) / float(NE)
        RMIN = arho_function(ESTART, IREG, rhob_function, rhocb_function,
                             rhoa_function, rhovb_function, rhod_function, semi, fjint_function)
        IESAV = 1

        for IE in range(2, NE + 1):
            ENER = ESTART + (IE - 1) * DELE
            RTMP = arho_function(ENER, IREG, rhob_function, rhocb_function,
                                 rhoa_function, rhovb_function, rhod_function, semi, fjint_function)
            if RTMP <= RMIN:
                IESAV = IE
                RMIN = RTMP

        EFMIN = ESTART + (IESAV - 2) * DELE
        EFMAX = ESTART + IESAV * DELE
        EF_OPT = gsect_function(arho_function, EFMIN, EFMAX, 1.e-6, IREG, rhob_function, rhocb_function,
                                rhoa_function, rhovb_function, rhod_function, semi, fjint_function)
        EF = (EFMIN + EFMAX) / 2

    semi["IINV"][IREG] = IINVSAV
    return EF


def semirho(IREG, DELE, ESTART, NE, NEDIM, RHOBTAB, ICOMP, RHOCBTAB, RHOVBTAB, rhocb_function, rhoa_function, rhovb_function, rhod_function, semi, fjint_function):
    if NE > NEDIM:
        raise ValueError('*** ERROR - NE > NEDIM; PROGRAM HALTED')
    for I in range(1, NE + 1):
        EF1 = (I - 1) * DELE + ESTART
        RHOCBSAV = rhocb_function(IREG, EF1, 0, semi, fjint_function)
        RHOVBSAV = rhovb_function(IREG, EF1, 0, semi, fjint_function)
        RHOBTAB[IREG, I] = -RHOCBSAV - \
            rhoa_function(IREG, EF1, 0, semi) + RHOVBSAV + \
            rhod_function(IREG, EF1, 0, semi)
        if ICOMP == 1:
            RHOCBTAB[IREG, I] = RHOCBSAV
            RHOVBTAB[IREG, I] = RHOVBSAV