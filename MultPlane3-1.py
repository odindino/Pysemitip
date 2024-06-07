import numpy as np

# Constants
EPSIL0 = 8.854185E-12
E = 1.60210E-19
PI = 4.0 * np.arctan(1.0)

# Parameters (example values)
NRDIM, NVDIM, NSDIM, NPDIM = 512, 64, 512, 64
NVDIM1, NVDIM2 = NVDIM + 1, 2048
NSDIM2, NEDIM = 20000, 50000
NXDIM, NXDIM2 = 64, 126
NYDIM, NZDIM, NEIGENDIM, NREGDIM, NARDIM = 64, 64, 8200, 2, 2

# Arrays initialization
VAC = np.zeros((2, NRDIM, NVDIM, NPDIM))
SEM = np.zeros((2, NRDIM, NSDIM, NPDIM))
VSINT = np.zeros((2, NRDIM, NPDIM))
R = np.zeros(NRDIM)
S = np.zeros(NSDIM)
DELV = np.zeros(NRDIM)
ITMAX = np.zeros(10)
EP = np.zeros(10)
BBIAS = np.zeros(1000)
BARR = np.zeros(NVDIM1)
PROF = np.zeros(NSDIM)
EMAX = np.zeros(4)
PVAC = np.zeros((NXDIM2, NYDIM, NVDIM))
PSEM = np.zeros((NXDIM2, NYDIM, NZDIM))
PSURF = np.zeros((NXDIM2, NYDIM))
BARR2 = np.zeros(NVDIM2)
AVBL = np.zeros(NREGDIM)
AVBH = np.zeros(NREGDIM)
AVBSO = np.zeros(NREGDIM)
ESO = np.zeros(NREGDIM)
TIP = np.zeros((NRDIM, NVDIM, NPDIM), dtype=bool)

# Common blocks
SEMI = {
    "TK": 0,
    "EGAP": np.zeros(NREGDIM),
    "ED": np.zeros(NREGDIM),
    "EA": np.zeros(NREGDIM),
    "ACB": np.zeros(NREGDIM),
    "AVB": np.zeros(NREGDIM),
    "CD": np.zeros(NREGDIM),
    "CA": np.zeros(NREGDIM),
    "IDEG": np.zeros(NREGDIM),
    "IINV": np.zeros(NREGDIM),
    "DELVB": np.zeros(NREGDIM)
}

PROTRU = {
    "RAD2": 0
}

SURF = {
    "ISTK": 0,
    "TK1": 0,
    "EN0": np.zeros(NARDIM),
    "EN": np.zeros((NARDIM, 2)),
    "DENS": np.zeros((NARDIM, 2)),
    "FWHM": np.zeros((NARDIM, 2)),
    "ECENT": np.zeros((NARDIM, 2))
}

CD_BLOCK = {
    "EF": 0,
    "ESTART": 0,
    "DELE": 0,
    "NE": 0,
    "RHOBTAB": np.zeros((NREGDIM, NEDIM)),
    "RHOSTAB": np.zeros((NARDIM, NEDIM)),
    "XSTEP1": 0,
    "XSTEP2": 0
}

TIPPOS = {
    "X0": 0,
    "Y0": 0
}


def gsect(f, xmin, xmax, ep, *args):
    GS = 0.3819660
    if xmax == xmin or ep == 0:
        return (xmin + xmax) / 2
    if xmax < xmin:
        xmin, xmax = xmax, xmin

    delx = xmax - xmin
    xa = xmin + delx * GS
    fa = f(xa, *args)
    xb = xmax - delx * GS
    fb = f(xb, *args)

    while delx >= ep:
        delxsav = delx
        if fb < fa:
            xmax = xb
            delx = xmax - xmin
            if delx == delxsav:
                return (xmin + xmax) / 2
            xb = xa
            fb = fa
            xa = xmin + delx * GS
            fa = f(xa, *args)
        else:
            xmin = xa
            delx = xmax - xmin
            if delx == delxsav:
                return (xmin + xmax) / 2
            xa = xb
            fa = fb
            xb = xmax - delx * GS
            fb = f(xb, *args)
    return (xmin + xmax) / 2

# # fermi-dirac integral


# def fjint(J, ETA, trap_function):
#     PI = 3.141592654
#     if ETA > 40:
#         return np.sqrt(ETA ** (J + 2)) / (J / 2.0 + 1)
#     if ETA < -8:
#         result = np.sqrt(PI) * np.exp(max(-40.0, ETA)) / 2.0
#         if J == 3:
#             result *= 3.0 / 2.0
#         return result
#     return trap_function(fj, 0.0, 20.0 + ETA, 1000)

# # trapezoidal integration routine


# def trap(f, XMIN, XMAX, NSTEP):
#     DELX = (XMAX - XMIN) / float(NSTEP)
#     SUM = (f(XMIN) + f(XMAX)) / 2.0
#     if NSTEP >= 2:
#         for i in range(1, NSTEP):
#             X = XMIN + i * DELX
#             SUM += f(X)
#     return SUM * DELX


# def fj(X, J, ETA):
#     return np.sqrt(X ** J) / (1.0 + np.exp(X - ETA))

# # fermi-dirac occupation factor


# def fd(e, ef, tk):
#     if tk == 0:
#         if e == ef:
#             return 0.5
#         return 1.0 if e < ef else 0.0
#     ratio = (e - ef) / tk
#     if ratio > 40:
#         return 0.0
#     if ratio < -40:
#         return 1.0
#     return 1.0 / (1.0 + np.exp(ratio))

# # electron density in conduction band


# def rhocb(IREG, EF, Pot, semi, fjint_function):
#     C = 6.815e21
#     if semi["IINV"][IREG] == 2 or semi["IINV"][IREG] == 3:
#         return 0.0
#     if semi["TK"] != 0.0:
#         return C * np.sqrt((semi["ACB"][IREG] * semi["TK"]) ** 3) * fjint_function(1, (EF - semi["EGAP"][IREG] - semi["DELVB"][IREG] - Pot) / semi["TK"])
#     if (EF - semi["EGAP"][IREG] - semi["DELVB"][IREG] - Pot) <= 0.0:
#         return 0.0
#     return (2.0 * C / 3.0) * np.sqrt((semi["ACB"][IREG] * (EF - semi["EGAP"][IREG] - semi["DELVB"][IREG] - Pot)) ** 3)

# # hole density in valence band


# def rhovb(IREG, EF, Pot, semi, fjint_function):
#     C = 6.815e21
#     if semi["IINV"][IREG] == 1 or semi["IINV"][IREG] == 3:
#         return 0.0
#     if semi["TK"] != 0.0:
#         return C * np.sqrt((semi["AVB"][IREG] * semi["TK"]) ** 3) * fjint_function(1, (-EF + semi["DELVB"][IREG] + Pot) / semi["TK"])
#     if (-EF + semi["DELVB"][IREG] + Pot) <= 0.0:
#         return 0.0
#     return (2.0 * C / 3.0) * np.sqrt((semi["AVB"][IREG] * (-EF + semi["DELVB"][IREG] + Pot)) ** 3)

# # ionized donor density


# def rhod(IREG, EF, Pot, semi):
#     RHOD = semi["CD"][IREG]
#     if semi["IDEG"][IREG] == 1:
#         return RHOD

#     EXPO = EF - semi["EGAP"][IREG] - \
#         semi["DELVB"][IREG] + semi["ED"][IREG] - Pot
#     if semi["TK"] != 0.0:
#         EXPO /= semi["TK"]
#         if EXPO < -40:
#             return RHOD
#         if EXPO > 40:
#             return 0.0
#         return RHOD / (1.0 + 2.0 * np.exp(EXPO))

#     if EXPO <= 0.0:
#         return RHOD
#     return 0.0

# # ionized acceptor density


# def rhoa(IREG, EF, Pot, semi):
#     RHOA = semi["CA"][IREG]
#     if semi["IDEG"][IREG] == 1:
#         return RHOA

#     EXPO = semi["EA"][IREG] - EF + semi["DELVB"][IREG] + Pot
#     if semi["TK"] != 0.0:
#         EXPO /= semi["TK"]
#         if EXPO < -40:
#             return RHOA
#         if EXPO > 40:
#             return 0.0
#         return RHOA / (1.0 + 4.0 * np.exp(EXPO))

#     if EXPO <= 0.0:
#         return RHOA
#     return 0.0

# # Total density of bulk charges


# def rhob(IREG, EF, Pot, rhocb_function, rhoa_function, rhovb_function, rhod_function, semi, fjint_function):
#     RHOE = rhocb_function(IREG, EF, Pot, semi, fjint_function) + \
#         rhoa_function(IREG, EF, Pot, semi)
#     RHOH = rhovb_function(IREG, EF, Pot, semi, fjint_function) + \
#         rhod_function(IREG, EF, Pot, semi)
#     return -RHOE + RHOH


# def rhos(IAR, EF1, DELE, rhos1_function, rhos2_function):
#     return rhos1_function(IAR, EF1, DELE) + rhos2_function(IAR, EF1, DELE)


# def arho(EF, IREG, rhob_function, rhocb_function, rhoa_function, rhovb_function, rhod_function, semi, fjint_function):
#     return abs(rhob_function(IREG, EF, 0, rhocb_function, rhoa_function, rhovb_function, rhod_function, semi, fjint_function))


# def enfind(IAR, EN1, EN2, NE, rhos_function):
#     EN0 = EN1
#     ESTART = EN1
#     DELE = abs(EN1 - EN2) / float(NE)
#     if DELE == 0:
#         return EN0
#     if EN2 < EN1:
#         DELE = -DELE

#     for IE in range(NE + 1):
#         EF1 = ESTART + IE * DELE
#         SIGTMP = rhos_function(IAR, EF1, abs(DELE))
#         if (DELE > 0 and SIGTMP <= 0) or (DELE <= 0 and SIGTMP >= 0):
#             EN0 = EF1
#             break
#     return EN0

def sig(iar, index, e):
    # Placeholder for the SIG function implementation
    # This function needs to be implemented based on your requirements
    return 1.0  # example value


def rhos1(iar, ef1, dele, surf):
    sum_val = 0.0
    e = surf["EN"][iar, 0]

    if ef1 == surf["EN"][iar, 0]:
        return sum_val

    if surf["ISTK"] == 0:
        if ef1 < surf["EN"][iar, 0]:
            while e < ef1:
                e += dele
                sum_val += sig(iar, 0, e) * dele
        else:
            while e > ef1:
                e -= dele
                sum_val += sig(iar, 0, e) * dele
        return sum_val

    if ef1 < surf["EN"][iar, 0]:
        while e < ef1 + 10.0 * surf["TK"]:
            e += dele
            sum_val += sig(iar, 0, e) * fd(e, ef1, surf["TK"]) * dele
    else:
        while e > ef1 - 10.0 * surf["TK"]:
            e -= dele
            sum_val += sig(iar, 0, e) * (1.0 - fd(e, ef1, surf["TK"])) * dele

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
                sum_val += sig(iar, 1, e) * dele
        else:
            while e > ef1:
                e -= dele
                sum_val += sig(iar, 1, e) * dele
        return sum_val

    if ef1 < surf["EN"][iar, 1]:
        while e < ef1 + 10.0 * surf["TK"]:
            e += dele
            sum_val += sig(iar, 1, e) * fd(e, ef1, surf["TK"]) * dele
    else:
        while e > ef1 - 10.0 * surf["TK"]:
            e -= dele
            sum_val += sig(iar, 1, e) * (1.0 - fd(e, ef1, surf["TK"])) * dele

    return sum_val


def rhos(iar, ef1, dele, surf):
    return rhos1(iar, ef1, dele, surf) + rhos2(iar, ef1, dele, surf)


def sigsum(iar, ener):
    return sig(iar, 0, ener) + sig(iar, 1, ener)


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

# # Find Fermi-level


# def effind(IREG, semi, arho_function, gsect_function, rhob_function, rhocb_function, rhoa_function, rhovb_function, rhod_function, fjint_function):
#     IINVSAV = semi["IINV"][IREG]
#     semi["IINV"][IREG] = 0

#     if semi["CD"][IREG] == 0 and semi["CA"][IREG] == 0:
#         EF = semi["EGAP"][IREG] / 2 + 0.75 * semi["TK"] * \
#             np.log(semi["AVB"][IREG] / semi["ACB"][IREG])
#     elif semi["TK"] == 0:
#         if semi["CA"][IREG] > semi["CD"][IREG]:
#             EF = semi["EA"][IREG] / 2
#         elif semi["CA"][IREG] < semi["CD"][IREG]:
#             EF = semi["EGAP"][IREG] - semi["ED"][IREG] / 2
#         else:
#             EF = (semi["EGAP"][IREG] - semi["ED"][IREG] + semi["EA"][IREG]) / 2
#     else:
#         ESTART = -0.1 + semi["DELVB"][IREG]
#         NE = 1000
#         DELE = (semi["EGAP"][IREG] + semi["DELVB"][IREG] + 0.2) / float(NE)
#         RMIN = arho_function(ESTART, IREG, rhob_function, rhocb_function,
#                              rhoa_function, rhovb_function, rhod_function, semi, fjint_function)
#         IESAV = 1

#         for IE in range(2, NE + 1):
#             ENER = ESTART + (IE - 1) * DELE
#             RTMP = arho_function(ENER, IREG, rhob_function, rhocb_function,
#                                  rhoa_function, rhovb_function, rhod_function, semi, fjint_function)
#             if RTMP <= RMIN:
#                 IESAV = IE
#                 RMIN = RTMP

#         EFMIN = ESTART + (IESAV - 2) * DELE
#         EFMAX = ESTART + IESAV * DELE
#         gsect_function(arho_function, EFMIN, EFMAX, 1.e-6, IREG, rhob_function, rhocb_function,
#                        rhoa_function, rhovb_function, rhod_function, semi, fjint_function)
#         EF = (EFMIN + EFMAX) / 2

#     semi["IINV"][IREG] = IINVSAV
#     return EF

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


# Function to read fort_new.9 file
def read_fort9(filename="fort_new.9"):
    with open(filename, 'r') as file:
        data = file.readlines()

    # 去除注釋部分
    data = [line.split()[0] for line in data if line.strip()
            and not line.lstrip().startswith('C')]

    index = 0

    NPARM = int(data[index].strip())
    index += 1

    parameters = []
    for _ in range(NPARM):
        param = {}
        param["SLOPE"] = float(data[index].strip())
        index += 1
        param["SEPIN"] = float(data[index].strip())
        index += 1
        param["RAD"] = float(data[index].strip())
        index += 1
        param["RAD2"] = float(data[index].strip())
        index += 1
        param["CPot"] = float(data[index].strip())
        index += 1
        param["X0"] = float(data[index].strip())
        index += 1
        param["Y0"] = float(data[index].strip())
        index += 1
        NREG = int(data[index].strip())
        index += 1
        param["NREG"] = NREG

        regions = []
        for _ in range(NREG):
            region = {}
            region["CD"] = float(data[index].strip())
            index += 1
            region["CA"] = float(data[index].strip())
            index += 1
            region["EGAP"] = float(data[index].strip())
            index += 1
            region["DELVB"] = float(data[index].strip())
            index += 1
            region["ED"] = float(data[index].strip())
            index += 1
            region["EA"] = float(data[index].strip())
            index += 1
            region["ACB"] = float(data[index].strip())
            index += 1
            region["AVBH"] = float(data[index].strip())
            index += 1
            region["AVBL"] = float(data[index].strip())
            index += 1
            region["AVBSO"] = float(data[index].strip())
            index += 1
            region["ESO"] = float(data[index].strip())
            index += 1
            region["IDEG"] = int(data[index].strip())
            index += 1
            region["IINV"] = int(data[index].strip())
            index += 1
            regions.append(region)
        param["REGIONS"] = regions

        param["EPSIL"] = float(data[index].strip())
        index += 1
        param["TEM"] = float(data[index].strip())
        index += 1
        param["TK"] = param["TEM"] * 8.617E-5
        param["TK1"] = param["TK"]

        NAR = int(data[index].strip())
        index += 1
        param["NAR"] = NAR

        areas = []
        for _ in range(NAR):
            area = {}
            area["DENS1"] = float(data[index].strip())
            index += 1
            area["EN1"] = float(data[index].strip())
            index += 1
            area["FWHM1"] = float(data[index].strip())
            index += 1
            area["ECENT1"] = float(data[index].strip())
            index += 1
            area["DENS2"] = float(data[index].strip())
            index += 1
            area["EN2"] = float(data[index].strip())
            index += 1
            area["FWHM2"] = float(data[index].strip())
            index += 1
            area["ECENT2"] = float(data[index].strip())
            index += 1
            areas.append(area)
        param["AREAS"] = areas

        param["ISTK"] = int(data[index].strip())
        index += 1

        param["MIRROR"] = int(data[index].strip())
        index += 1

        param["NRIN"] = int(data[index].strip())
        index += 1
        param["NVIN"] = int(data[index].strip())
        index += 1
        param["NSIN"] = int(data[index].strip())
        index += 1
        param["NPIN"] = int(data[index].strip())
        index += 1
        param["SIZE"] = float(data[index].strip())
        index += 1
        if param["SIZE"] <= 0:
            param["DELRIN"] = float(data[index].strip())
            index += 1
            param["DELSIN"] = float(data[index].strip())
            index += 1

        param["IPMAX"] = int(data[index].strip())
        index += 1
        # 將 ITMAX 的值從逗號分隔的字符串轉換為浮點數列表，過濾掉空字符串
        param["ITMAX"] = [int(float(x))
                          for x in data[index].split(',') if x.strip()]
        index += 1

        # 將 EP 的值從逗號分隔的字符串轉換為浮點數列表，過濾掉空字符串
        param["EP"] = [float(x) for x in data[index].split(',') if x.strip()]
        index += 1

        param["NE"] = int(data[index].strip())
        index += 1

        param["IWRIT"] = int(data[index].strip())
        index += 1

        # 檢查下一行是否包含逗號分隔的數值串並處理 NBIAS 和 BBIAS

        param["NBIAS"] = int(data[index].strip())
        index += 1

        # Read BBIAS values from the next line
        # change '2.0, 1.9, 1.8 ...' string to list of floats
        param["BBIAS"] = [float(x)
                          for x in data[index].split(',') if x.strip()]
        index += 1

        param["NUMC"] = int(data[index].strip())
        index += 1
        param["DELPOT"] = float(data[index].strip())
        index += 1
        param["PhiIN"] = float(data[index].strip())
        index += 1
        param["CHI"] = float(data[index].strip())
        index += 1
        param["EFTIP"] = float(data[index].strip())
        index += 1
        param["VACWID"] = float(data[index].strip())
        index += 1
        param["ZVACDEL"] = float(data[index].strip())
        index += 1
        param["EMAX"] = [float(x) for x in data[index].split(',') if x.strip()]
        index += 1
        param["ICOMPIN"] = int(data[index].strip())
        index += 1
        param["BMOD"] = float(data[index].strip())
        index += 1
        param["ANEG"] = float(data[index].strip())
        index += 1
        param["APOS"] = float(data[index].strip())
        index += 1
        param["VSTART"] = float(data[index].strip())
        index += 1
        if param["VSTART"] < 0:
            param["DELS1"] = abs(param["VSTART"]) * param["ANEG"]
        else:
            param["DELS1"] = abs(param["VSTART"]) * param["APOS"]
        param["SEP0"] = param["SEPIN"] - param["DELS1"]

        parameters.append(param)

    return parameters

# Function to process the parameters and initialize common blocks and arrays


def process_parameters(parameters):
    for param in parameters:
        SLOPE = param["SLOPE"]
        SEPIN = param["SEPIN"]
        RAD = param["RAD"]
        RAD2 = param["RAD2"]
        CPot = param["CPot"]
        X0 = param["X0"]
        Y0 = param["Y0"]
        NREG = param["NREG"]

        CC = {
            "CD": np.array([region["CD"] for region in param["REGIONS"]]),
            "CA": np.array([region["CA"] for region in param["REGIONS"]])
        }
        SEMI = {
            "TK": param["TK"],
            "EGAP": np.array([region["EGAP"] for region in param["REGIONS"]]),
            "ED": np.array([region["ED"] for region in param["REGIONS"]]),
            "EA": np.array([region["EA"] for region in param["REGIONS"]]),
            "ACB": np.array([region["ACB"] for region in param["REGIONS"]]),
            "AVB": np.array([np.exp(2.0 * np.log(np.sqrt(region["AVBH"]**3) + np.sqrt(region["AVBL"]**3)) / 3.0) for region in param["REGIONS"]]),
            "IDEG": np.array([region["IDEG"] for region in param["REGIONS"]], dtype=int),
            "IINV": np.array([region["IINV"] for region in param["REGIONS"]], dtype=int),
            "DELVB": np.array([region["DELVB"] for region in param["REGIONS"]]),
            "CD": CC["CD"],
            "CA": CC["CA"]
        }

        print(
            f"RAD, SLOPE, ANGLE = {RAD}, {SLOPE}, {360.0 * np.arctan(1.0 / SLOPE) / PI}")
        print(f"CONTACT POTENTIAL = {CPot}")
        print(f"POSITION OF TIP = {X0}, {Y0}")
        print(f"NUMBER OF DIFFERENT REGIONS OF SEMICONDUCTOR = {NREG}")

        if NREG > NREGDIM:
            print("INPUT NUMBER OF REGIONS > NUMBER OF REGIONS IN PARAMETER STATEMENT")
            print("PROGRAM WILL BE EXITED (TYPE RETURN)")
            input()
            exit()

        for ireg, region in enumerate(param["REGIONS"]):
            if NREG > 1:
                print(f"REGION #{ireg + 1}")

            CC["CD"][ireg] = region["CD"]
            CC["CA"][ireg] = region["CA"]
            SEMI["EGAP"][ireg] = region["EGAP"]
            SEMI["DELVB"][ireg] = region["DELVB"]
            SEMI["ED"][ireg] = region["ED"]
            SEMI["EA"][ireg] = region["EA"]
            SEMI["ACB"][ireg] = region["ACB"]
            AVBH[ireg] = region["AVBH"]
            AVBL[ireg] = region["AVBL"]
            SEMI["AVB"][ireg] = np.exp(
                2.0 * np.log(np.sqrt(AVBH[ireg]**3) + np.sqrt(AVBL[ireg]**3)) / 3.0)
            AVBSO[ireg] = region["AVBSO"]
            ESO[ireg] = region["ESO"]
            SEMI["IDEG"][ireg] = region["IDEG"]
            SEMI["IINV"][ireg] = region["IINV"]

            print(f"DOPING = {CC['CD'][ireg]}, {CC['CA'][ireg]}")
            print(
                f"BAND GAP, VB OFFSET = {SEMI['EGAP'][ireg]}, {SEMI['DELVB'][ireg]}")

            if (CC["CA"][ireg] > CC["CD"][ireg] and (SEMI["IINV"][ireg] == 1 or SEMI["IINV"][ireg] == 3)) or \
               (CC["CD"][ireg] > CC["CA"][ireg] and (SEMI["IINV"][ireg] == 2 or SEMI["IINV"][ireg] == 3)):
                print(
                    "****** WARNING - LIKELY INCOMPATIBLE DOPING AND INVERSION (IINV) PARAMETER")
                print("CONTINUE (y/n) ?")
                ans = input().strip().lower()
                if ans != 'y':
                    exit()

        EPSIL = param["EPSIL"]
        TEM = param["TEM"]
        TK = param["TK"]
        TK1 = param["TK1"]

        NAR = param["NAR"]
        print(f"NUMBER OF DIFFERENT AREAS OF SURFACE STATES = {NAR}")

        if NAR > NARDIM:
            print("INPUT NUMBER OF AREAS > NUMBER OF AREAS IN PARAMETER STATEMENT")
            print("PROGRAM WILL BE EXITED (TYPE RETURN)")
            input()
            exit()

        for iar, area in enumerate(param["AREAS"]):
            if NAR > 1:
                print(f"AREA #{iar + 1}")

            SURF["DENS"][iar, 0] = area["DENS1"]
            SURF["EN"][iar, 0] = area["EN1"]
            SURF["FWHM"][iar, 0] = area["FWHM1"]
            SURF["ECENT"][iar, 0] = area["ECENT1"]
            SURF["DENS"][iar, 1] = area["DENS2"]
            SURF["EN"][iar, 1] = area["EN2"]
            SURF["FWHM"][iar, 1] = area["FWHM2"]
            SURF["ECENT"][iar, 1] = area["ECENT2"]

            print("FIRST DISTRIBUTION OF SURFACE STATES:")
            print(
                f"SURFACE STATE DENSITY, EN = {SURF['DENS'][iar, 0]}, {SURF['EN'][iar, 0]}")
            print(
                f"FWHM, ECENT = {SURF['FWHM'][iar, 0]}, {SURF['ECENT'][iar, 0]}")
            print("SECOND DISTRIBUTION OF SURFACE STATES:")
            print(
                f"SURFACE STATE DENSITY, EN = {SURF['DENS'][iar, 1]}, {SURF['EN'][iar, 1]}")
            print(
                f"FWHM, ECENT = {SURF['FWHM'][iar, 1]}, {SURF['ECENT'][iar, 1]}")

        SURF["ISTK"] = param["ISTK"]

        # 初始化其它參數
        MIRROR = param["MIRROR"]
        if MIRROR == 1:
            print("HORIZONTAL MIRROR PLANE ASSUMED")

        NRIN = param["NRIN"]
        NVIN = param["NVIN"]
        NSIN = param["NSIN"]
        NPIN = param["NPIN"]
        SIZE = param["SIZE"]
        if SIZE <= 0:
            DELRIN = param["DELRIN"]
            DELSIN = param["DELSIN"]

        IPMAX = param["IPMAX"]
        ITMAX = param["ITMAX"]
        EP = param["EP"]
        NE = param["NE"]

        # 找出整體電荷中性水平 EN0
        EN0 = np.zeros(NARDIM)
        for IAR in range(NAR):
            if SURF["DENS"][IAR, 0] == 0 and SURF["DENS"][IAR, 1] == 0:
                EN0[IAR] = 0.0
            elif SURF["DENS"][IAR, 0] == 0:
                EN0[IAR] = SURF["EN"][IAR, 1]
            elif SURF["DENS"][IAR, 1] == 0:
                EN0[IAR] = SURF["EN"][IAR, 0]
            else:
                print("SEARCHING FOR CHARGE NEUTRALITY LEVEL")
                EN0[IAR] = enfind(IAR, SURF["EN"][IAR, 0],
                                  SURF["EN"][IAR, 1], NE, rhos)

            print(f"CHARGE-NEUTRALITY LEVEL = {EN0[IAR]}")

        EN0MAX = np.max(EN0)
        EN0MIN = np.min(EN0)

        IWRIT = param["IWRIT"]
        NBIAS = param["NBIAS"]
        BBIAS = param["BBIAS"]
        NUMC = param["NUMC"]
        DELPOT = param["DELPOT"]
        PhiIN = param["PhiIN"]
        CHI = param["CHI"]
        EFTIP = param["EFTIP"]
        VACWID = param["VACWID"]
        ZVACDEL = param["ZVACDEL"]
        EMAX = param["EMAX"]
        ICOMPIN = param["ICOMPIN"]
        BMOD = param["BMOD"]
        ANEG = param["ANEG"]
        APOS = param["APOS"]
        VSTART = param["VSTART"]
        if VSTART < 0:
            DELS1 = abs(VSTART) * ANEG
        else:
            DELS1 = abs(VSTART) * APOS
        SEP0 = SEPIN - DELS1

        # Find Fermi-level position
        for ireg in range(NREG):
            EF = effind(ireg, SEMI, arho, gsect, rhob,
                        rhocb, rhoa, rhovb, rhod, fjint)
            print(f"REGION TYPE {ireg + 1}, FERMI-LEVEL = {EF}")
            RHOCC = rhocb(ireg, EF, 0, SEMI, fjint)
            RHOVV = rhovb(ireg, EF, 0, SEMI, fjint)
            print(f"CARRIER DENSITY IN CB, VB = {RHOCC}, {RHOVV}")


# Example usage
parameters = read_fort9("fort_new.9")
process_parameters(parameters)
