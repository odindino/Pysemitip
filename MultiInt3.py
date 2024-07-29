import numpy as np
import struct
# 假设这些模块已经存在并且导入
from gsect import gsect
from semirhomult import effind, arho, rhob, rhocb, rhovb, rhoa, rhod, fjint, fd
from surfrhomult import enfind, rhos
from semitip3 import *
from potcut3 import *
from math import atan, sqrt, log, pi, tan, cos, sin, exp
from surfrho import *
from intcurr import*
from contr3 import*

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
CURR = 0.0
CURRV = 0.0
CURRC = 0.0
DELPHI = 0.0
PHI0 = 0.0
IERR = 0
ACB = np.zeros(NREGDIM)

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

def read_fort9(filename="fort_new.9"):
    with open(filename, 'r') as file:
        data = file.readlines()

    data = [line.split()[0] for line in data if line.strip() and not line.lstrip().startswith('C')]

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
        param["EF"] = None
        param["EN0"] = []

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
        param["ITMAX"] = [int(float(x)) for x in data[index].split(',') if x.strip()]
        index += 1
        param["EP"] = [float(x) for x in data[index].split(',') if x.strip()]
        index += 1

        param["NE"] = int(data[index].strip())
        index += 1

        param["IWRIT"] = int(data[index].strip())
        index += 1

        param["NBIAS"] = int(data[index].strip())
        index += 1
        param["BBIAS"] = [float(x) for x in data[index].split(',') if x.strip()]
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

        param["NBARR1"] = 10

        parameters.append(param)

    return parameters

def compute_depletion_width(CD, CA, PotTIP, EPSIL, EPSIL0, E):
    if (CD[0] - CA[0]) == 0.0:
        W = 1.0e10
    else:
        W = 1.0e9 * np.sqrt(2.0 * EPSIL * EPSIL0 * max(1.0, abs(PotTIP)) / (abs(CD[0] - CA[0]) * 1.0e6 * E))
    return W

def semirho(IREG, DELE, ESTART, NE, RHOBTAB):
    RHOBTAB[IREG, :NE] = np.linspace(ESTART, ESTART + DELE * (NE - 1), NE)

def compute_energy_table(EF, PotTIP, EN0MIN, EN0MAX, NE, NREG, RHOBTAB):
    ESTART = min(EF, EF - PotTIP, EN0MIN)
    EEND = max(EF, EF - PotTIP, EN0MAX)
    ETMP = EEND - ESTART
    ESTART -= 2.0 * ETMP
    EEND += 2.0 * ETMP
    DELE = (EEND - ESTART) / float(NE - 1)

    NETMP = int(round((EF - ESTART) / DELE))
    ESTART = EF - (NETMP - 0.5) * DELE
    EEND = ESTART + (NE - 1) * DELE

    for IREG in range(NREG):
        semirho(IREG, DELE, ESTART, NE, RHOBTAB)

    return ESTART, EEND, DELE

def compute_surface_charge_density_table(IAR, DELE, ESTART, NE, NEDIM, RHOSTAB):
    if NE > NEDIM:
        print('*** ERROR - NE > NEDIM; PROGRAM HALTED')
        input('TYPE ENTER TO CONTINUE')
        return

    if ISTK == 1:
        for i in range(NE):
            ef1 = i * DELE + ESTART
            if DENS[1] == 0.0:
                RHOSTAB[i] = rhos1(ef1, DELE)
            elif DENS[0] == 0.0:
                RHOSTAB[i] = rhos2(ef1, DELE)
            else:
                RHOSTAB[i] = rhos(ef1, DELE)
    else:
        if DENS[0] == 0.0 or DENS[1] == 0.0:
            nen = int((EN0 - ESTART) / DELE) + 1
            RHOSTAB[nen] = 0.0
            summ = 0.0
            for i in range(nen + 1, NE):
                ef1 = i * DELE + ESTART
                summ += sigsum(ef1)
                RHOSTAB[i] = summ * DELE
            summ = 0.0
            for i in range(nen - 1, -1, -1):
                ef1 = i * DELE + ESTART
                summ += sigsum(ef1)
                RHOSTAB[i] = summ * DELE
        else:
            nen = int((EN0 - ESTART) / DELE) + 1
            RHOSTAB[nen] = 0.0
            for i in range(nen + 1, NE):
                ef1 = i * DELE + ESTART
                RHOSTAB[i] = rhos(ef1, DELE)
            for i in range(nen - 1, -1, -1):
                ef1 = i * DELE + ESTART
                RHOSTAB[i] = rhos(ef1, DELE)

def compute_bias_and_tip_potential(bias0, BMOD, CPot):
    IMODMAX = 1
    if BMOD == 0.0:
        IMODMAX = -1

    biases = []
    tip_potentials = []

    for imod in range(-1, IMODMAX + 1, 2):
        BIAS = bias0 + imod * BMOD * np.sqrt(2.0)
        PotTIP = BIAS + CPot
        biases.append(BIAS)
        tip_potentials.append(PotTIP)

    return biases, tip_potentials

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
        NBARR1 = param["NBARR1"]
        NLOC = np.zeros(4, dtype=int)
        
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
            "AVB": np.array([exp(2.0 * log(sqrt(region["AVBH"]**3) + sqrt(region["AVBL"]**3)) / 3.0) for region in param["REGIONS"]]),
            "IDEG": np.array([region["IDEG"] for region in param["REGIONS"]], dtype=int),
            "IINV": np.array([region["IINV"] for region in param["REGIONS"]], dtype=int),
            "DELVB": np.array([region["DELVB"] for region in param["REGIONS"]]),
            "CD": CC["CD"],
            "CA": CC["CA"]
        }
        AVBL = np.array([region["AVBL"] for region in param["REGIONS"]])
        AVBH = np.array([region["AVBH"] for region in param["REGIONS"]])
        AVBSO = np.array([region["AVBSO"] for region in param["REGIONS"]])
        ESO = np.array([region["ESO"] for region in param["REGIONS"]])
        # 默认值设置
        FRACZ = param.get("FRACZ", 0)
        E2HH = param.get("E2HH", 0)
        E2SO = param.get("E2SO", 0)
        NEE = param.get("NEE", 0)
        NWK = param.get("NWK", 0)
        EXPANI = param.get("EXPANI", 0)

        # 初始化 TMP 和 IINV 变量
        TMP = np.zeros(NPDIM)
        IINV = np.array([region["IINV"] for region in param["REGIONS"]], dtype=int)

        print(f"RAD, SLOPE, ANGLE = {RAD}, {SLOPE}, {360.0 * atan(1.0 / SLOPE) / pi}")
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
            SEMI["AVB"][ireg] = exp(2.0 * log(sqrt(AVBH[ireg]**3) + sqrt(AVBL[ireg]**3)) / 3.0)
            SEMI["IDEG"][ireg] = region["IDEG"]
            SEMI["IINV"][ireg] = region["IINV"]

            print(f"DOPING = {CC['CD'][ireg]}, {CC['CA'][ireg]}")
            print(f"BAND GAP, VB OFFSET = {SEMI['EGAP'][ireg]}, {SEMI['DELVB'][ireg]}")

            if (CC["CA"][ireg] > CC["CD"][ireg] and (SEMI["IINV"][ireg] == 1 or SEMI["IINV"][ireg] == 3)) or \
               (CC["CD"][ireg] > CC["CA"][ireg] and (SEMI["IINV"][ireg] == 2 or SEMI["IINV"][ireg] == 3)):
                print("****** WARNING - LIKELY INCOMPATIBLE DOPING AND INVERSION (IINV) PARAMETER")
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
            print(f"SURFACE STATE DENSITY, EN = {SURF['DENS'][iar, 0]}, {SURF['EN'][iar, 0]}")
            print(f"FWHM, ECENT = {SURF['FWHM'][iar, 0]}, {SURF['ECENT'][iar, 0]}")
            print("SECOND DISTRIBUTION OF SURFACE STATES:")
            print(f"SURFACE STATE DENSITY, EN = {SURF['DENS'][iar, 1]}, {SURF['EN'][iar, 1]}")
            print(f"FWHM, ECENT = {SURF['FWHM'][iar, 1]}, {SURF['ECENT'][iar, 1]}")

        SURF["ISTK"] = param["ISTK"]

        MIRROR = param["MIRROR"]
        if MIRROR == 1:
            print("HORIZONTAL MIRROR PLANE ASSUMED")

        NRIN = param["NRIN"]
        NVIN = param["NVIN"]
        NSIN = param["NSIN"]
        NPIN = param["NPIN"]
        SIZE = param["SIZE"]
        DELRIN = param.get("DELRIN", 1.0)
        DELSIN = param.get("DELSIN", 1.0)

        IPMAX = param["IPMAX"]
        ITMAX = param["ITMAX"]
        EP = param["EP"]
        NE = param["NE"]

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
                EN0[IAR] = enfind(IAR, SURF["EN"][IAR, 0], SURF["EN"][IAR, 1], NE, rhos)

            param["EN0"].append(EN0[IAR])
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
        NR = param["NRIN"]
        NV = param["NVIN"]
        NS = param["NSIN"]
        NP = param["NPIN"]
        BIAS = param["BBIAS"][0] if param["NBIAS"] > 0 else 0.0
        pot0 = 0.0
        ierr = 0
        IINIT = 1
        DELR = np.zeros(NRDIM)
        DELS = np.zeros(NSDIM)
        DELV = np.zeros(NRDIM)
        DELP = 2.0 * pi / NP
        R = np.zeros(NRDIM)
        S = np.zeros(NSDIM)
        DELXSI = np.zeros(NRDIM)
        VAC = np.zeros((2, NRDIM, NVDIM, NPDIM))
        TIP = np.zeros((NRDIM, NVDIM, NPDIM), dtype=bool)
        SEM = np.zeros((2, NRDIM, NSDIM, NPDIM))
        VSINT = np.zeros((2, NRDIM, NPDIM))
        ETAT, A, Z0, C = 0, 0, 0, 0

        if VSTART < 0:
            DELS1 = abs(VSTART) * ANEG
        else:
            DELS1 = abs(VSTART) * APOS
        SEP0 = param["SEPIN"] - DELS1

        for ireg in range(NREG):
            EF = effind(ireg, SEMI, arho, gsect, rhob, rhocb, rhoa, rhovb, rhod, fjint)
            param["EF"] = EF
            print(f"REGION TYPE {ireg + 1}, FERMI-LEVEL = {EF}")
            RHOCC = rhocb(ireg, EF, 0, SEMI, fjint)
            RHOVV = rhovb(ireg, EF, 0, SEMI, fjint)
            print(f"CARRIER DENSITY IN CB, VB = {RHOCC}, {RHOVV}")

        print(f"NBIAS: {NBIAS}, BBIAS: {BBIAS}")
        for ibias in range(min(NBIAS, len(BBIAS))):
            print(f"ibias: {ibias}")
            bias0 = BBIAS[ibias]
            if bias0 <= 0:
                SEP = SEP0 + ANEG * abs(bias0)
            else:
                SEP = SEP0 + APOS * abs(bias0)
            print(f"SEPARATION = {SEP:.8f}")

            biases, tip_potentials = compute_bias_and_tip_potential(bias0, BMOD, CPot)

            # 初始化 CSAV 相关变量
            CSAV, CSAVE, CSAVL, CSAVV, CSAVVE, CSAVVL, CSAVC, CSAVCE, CSAVCL = 0, 0, 0, 0, 0, 0, 0, 0, 0

            for imod in [-1, 1]:
                for BIAS, PotTIP in zip(biases, tip_potentials):
                    print(f"BIAS, TIP POTENTIAL = {BIAS:.8f}, {PotTIP:.8f}")

                    depletion_width = compute_depletion_width(CC["CD"], CC["CA"], PotTIP, EPSIL, EPSIL0, E)
                    print(f"1-D ESTIMATE OF DEPLETION WIDTH (NM) = {depletion_width:.8f}")

                    RHOBTAB = np.zeros((NREG, NEDIM))
                    ESTART, EEND, DELE = compute_energy_table(EF, PotTIP, EN0MIN, EN0MAX, NE, NREG, RHOBTAB)
                    print(f"ESTART,EEND,NE = {ESTART:.8f},{EEND:.8f},{NE}")
                    print("COMPUTING TABLE OF BULK CHARGE DENSITIES")
                    print("COMPUTING TABLE OF SURFACE CHARGE DENSITIES")

                    NR = NRIN
                    NV = NVIN
                    NS = NSIN

                    if SIZE > 0:
                        DELR0 = RAD
                    if RAD2 != 0:
                        DELR0 = min(RAD2, DELR)
                        DELR0 = min(DELR, depletion_width / NR) * SIZE
            
                    DELS0 = RAD
                    if RAD2 != 0:
                        DELS0 = min(RAD2, DELS)
                        DELS0 = min(DELS, depletion_width / NS) * SIZE
                    else:
                        DELR0 = DELRIN
                        DELS0 = DELSIN
            
                    NP = NPIN

                    if MIRROR == 1:
                        DELP = pi / float(NP)
                    else:
                        DELP = 2 * pi / float(NP)

                    ierr = 0
                    semitip3(SEP, RAD, SLOPE, DELRIN, DELSIN, VAC, TIP, SEM, VSINT, R, S, DELV, DELR, DELXSI, DELP, 
                        NRDIM, NVDIM, NSDIM, NPDIM, NR, NV, NS, NP, PotTIP, IWRIT, ITMAX, EP, IPMAX, pot0, ierr, 
                        IINIT, MIRROR, EPSIL)

                    print(f'RETURN FROM SEMTIP3, NR, NS, NV, IERR = {NR}, {NS}, {NV}, {ierr}')

                    print('COMPUTATION OF CURRENT:')
                    BARR = np.zeros(NVDIM1)
                    PROF = np.zeros(NSDIM)
                    CURRVE, CURRVL, CURRCE, CURRCL, CURRE, CURRL = 0, 0, 0, 0, 0, 0
                    # 确保 icomp 是一个单一的值，而不是数组
                    icomp = int(param["ICOMPIN"])

                    # 打印 icomp 的值以确认
                    print(f"icomp = {icomp}")
                    currv = np.zeros(3)  # 假設需要 3 個元素
                    currv0 = np.zeros(3)  # 假設需要 3 個元素
                    currc = np.zeros(1)  # 假設需要 1 個元素
                    currc0 = np.zeros(1)  # 假設需要 1 個元素
                    curr = np.zeros(1)  # 假設需要 1 個元素
                    curr0 = np.zeros(1)  
                    intcurr(BARR, PROF, NBARR1, NV, NS, NS, NVDIM, NSDIM, S, SEP, BIAS, EF, CHI, EFTIP, CPot, 
                            SEMI["EGAP"], TK, AVBH[0], AVBL[0], AVBSO[0], ACB[0], ESO[0], E2HH, E2SO, NEE, NWK, pot0, 
                            NVDIM1, NVDIM2, NSDIM2, EXPANI, NLOC, CURRVE, CURRVL, CURRCE, CURRCL, CURRE, CURRL, 
                            IWRIT, icomp, TMP, TMP, TMP, TMP, TMP, TMP)

                    if (IINV[0] == 1 or IINV[0] == 3) and CURRVE > 0:
                        CURRVE = 0.0
                        print('VB EXTENDED INVERSION CURRENT SET TO ZERO')
                    if (IINV[0] == 1 or IINV[0] == 3) and CURRVL > 0:
                        CURRVL = 0.0
                        print('VB LOCALIZED INVERSION CURRENT SET TO ZERO')
                    if (IINV[0] == 2 or IINV[0] == 3) and CURRCE < 0:
                        CURRCE = 0.0
                        print('CB EXTENDED INVERSION CURRENT SET TO ZERO')
                    if (IINV[0] == 2 or IINV[0] == 3) and CURRCL < 0:
                        CURRCL = 0.0
                        print('CB LOCALIZED INVERSION CURRENT SET TO ZERO')

                    print(f'valence band current ext, loc = {CURRVE}, {CURRVL}')
                    print(f'conduction band current ext, loc = {CURRCE}, {CURRCL}')

                    CURR = CURRE + CURRL
                    CURRV = CURRVE + CURRVL
                    CURRC = CURRCE + CURRCL

                    if imod == -1:
                        CSAV = CURR
                        CSAVE = CURRE
                        CSAVL = CURRL
                        CSAVV = CURRV
                        CSAVVE = CURRVE
                        CSAVVL = CURRVL
                        CSAVC = CURRC
                        CSAVCE = CURRCE
                        CSAVCL = CURRCL
                    else:
                        COND = (CURR - CSAV) * exp(2.0 * 10.0 * (SEP - SEP0 - DELS1)) / (2.0 * BMOD * sqrt(2.0))
                        CONDE = (CURRE - CSAVE) * exp(2.0 * 10.0 * (SEP - SEP0 - DELS1)) / (2.0 * BMOD * sqrt(2.0))
                        CONDL = (CURRL - CSAVL) * exp(2.0 * 10.0 * (SEP - SEP0 - DELS1)) / (2.0 * BMOD * sqrt(2.0))
                        CONDV = (CURRV - CSAVV) * exp(2.0 * 10.0 * (SEP - SEP0 - DELS1)) / (2.0 * BMOD * sqrt(2.0))
                        CONDVE = (CURRVE - CSAVVE) * exp(2.0 * 10.0 * (SEP - SEP0 - DELS1)) / (2.0 * BMOD * sqrt(2.0))
                        CONDVL = (CURRVL - CSAVVL) * exp(2.0 * 10.0 * (SEP - SEP0 - DELS1)) / (2.0 * BMOD * sqrt(2.0))
                        CONDC = (CURRC - CSAVC) * exp(2.0 * 10.0 * (SEP - SEP0 - DELS1)) / (2.0 * BMOD * sqrt(2.0))
                        CONDCE = (CURRCE - CSAVCE) * exp(2.0 * 10.0 * (SEP - SEP0 - DELS1)) / (2.0 * BMOD * sqrt(2.0))
                        CONDCL = (CURRCL - CSAVCL) * exp(2.0 * 10.0 * (SEP - SEP0 - DELS1)) / (2.0 * BMOD * sqrt(2.0))

                        print(f'BIAS0 = {bias0}, COND = {COND}, CONDE = {CONDE}, CONDL = {CONDL}')
                        print(f'CONDV = {CONDV}, CONDVE = {CONDVE}, CONDVL = {CONDVL}')
                        print(f'CONDC = {CONDC}, CONDCE = {CONDCE}, CONDCL = {CONDCL}')

                    CURR *= exp(2.0 * 10.0 * (SEP - SEP0 - DELS1))
                    CURRE *= exp(2.0 * 10.0 * (SEP - SEP0 - DELS1))
                    CURRL *= exp(2.0 * 10.0 * (SEP - SEP0 - DELS1))
                    CURRV *= exp(2.0 * 10.0 * (SEP - SEP0 - DELS1))
                    CURRVE *= exp(2.0 * 10.0 * (SEP - SEP0 - DELS1))
                    CURRVL *= exp(2.0 * 10.0 * (SEP - SEP0 - DELS1))
                    CURRC *= exp(2.0 * 10.0 * (SEP - SEP0 - DELS1))
                    CURRCE *= exp(2.0 * 10.0 * (SEP - SEP0 - DELS1))
                    CURRCL *= exp(2.0 * 10.0 * (SEP - SEP0 - DELS1))

                    print(f'BIAS = {BIAS}, CURR = {CURR}, CURRE = {CURRE}, CURRL = {CURRL}')
                    print(f'CURRV = {CURRV}, CURRVE = {CURRVE}, CURRVL = {CURRVL}')
                    print(f'CURRC = {CURRC}, CURRCE = {CURRCE}, CURRCL = {CURRCL}')

        # Plot cross-sectional profile
        if IWRIT >= 1:
            with open('output11.txt', 'w') as file11:
                for J in range(NV, 0, -1):
                    file11.write(f"{-J*DELV[0]:.8f} {VAC[0,0,J-1,0]:.8f} {VAC[0,NR-1,J-1,0]:.8f}\n")
                file11.write(f"{0:.8f} {VSINT[0,0,0]:.8f} {VSINT[0,NR-1,0]:.8f}\n")
                for J in range(1, NS+1):
                    file11.write(f"{S[J-1]:.8f} {SEM[0,0,J-1,0]:.8f} {SEM[0,NR-1,J-1,0]:.8f}\n")

        # Plot surface potential and surface charge density
        KPLOT1 = int(round(PhiIN / (DELP * 180.0 / np.pi) + 0.5))
        KPLOT1 = min(max(1, KPLOT1), NP)
        Phi = (KPLOT1 - 0.5) * DELP * 180.0 / np.pi
        if MIRROR == 1:
            KPLOT2 = NP - KPLOT1 + 1
        else:
            KPLOT2 = (KPLOT1 + NP // 2) % NP + 1

        print(f"ACTUAL ANGLE OF CROSS-SECTIONAL PLOT = {Phi}")
        print(f"CORRESPONDING TO ANGULAR GRID LINES {KPLOT2} {KPLOT1}")

        with open('output12.txt', 'w') as file12:
            for I in range(NR, 0, -1):
                file12.write(f"{-R[I-1]:.8f} {VSINT[0,I-1,KPLOT2-1]:.8f} {SEM[0,I-1,NS-1,KPLOT2-1]:.8f}\n")
            for I in range(1, NR+1):
                file12.write(f"{R[I-1]:.8f} {VSINT[0,I-1,KPLOT1-1]:.8f} {SEM[0,I-1,NS-1,KPLOT1-1]:.8f}\n")

        # 关闭文件
        file11.close()
        file12.close()
        ETA1=0
        # Plot contours
        if IWRIT >= 2:
            contr3(ETA1, VAC, TIP, SEM, VSINT, R, S, DELV, NRDIM, NVDIM, NSDIM, NPDIM, NR, NV, NS, NP, NUMC, DELPOT, MIRROR, KPLOT1, KPLOT2)

        # Output entire potential
        if IWRIT >= 3:
            nrecl = 40 + 4 * NR * NP * (NV + NS + 1) + 4 * NR * 2 + 4 * NS
            with open('fort.13', 'wb') as file13:
                data = struct.pack(
                    'i' * 10 + 'f' * 10 + 'f' * (NR * NV * NP + NR * NS * NP + NR * NP + NR + NS + NR),
                    NR, NV, NS, NP, SEP, RAD, RAD2, SLOPE, BIAS, EPSIL,
                    *(VAC[0, i, j, k] for k in range(NP) for j in range(NV) for i in range(NR)),
                    *(SEM[0, i, j, k] for k in range(NP) for j in range(NS) for i in range(NR)),
                    *(VSINT[0, i, k] for k in range(NP) for i in range(NR)),
                    *(R[i] for i in range(NR)),
                    *(S[j] for j in range(NS)),
                    *(DELV[i] for i in range(NR))
                )
                file13.write(data)

        # 按下 Enter 键退出
        input("PRESS THE ENTER KEY TO EXIT")

parameters = read_fort9()
process_parameters(parameters)

