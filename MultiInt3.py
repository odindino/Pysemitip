import math
import numpy as np
from semirhomult import effind, arho, rhob, rhocb, rhovb, rhoa, rhod, fjint, fd
from surfrhomult import enfind, rhos
from gsect import gsect
from semitip import*
from potcut3 import*
from intcurr import*
# Constants
EPSIL0 = 8.854185e-12
E = 1.60210e-19
pi = 4.0 * math.atan(1.0)

# Parameters
NRDIM = int(512)
NVDIM = int(64)
NSDIM = int(512)
NPDIM = int(64)
NVDIM1 = NVDIM + 1
NVDIM2 = int(2048)
NSDIM2 = int(20000)
NEDIM = int(50000)
NREGDIM = int(2)
NARDIM = int(2)
NBARR1 = [0]

# Arrays and Variables
VAC = np.zeros((2, NRDIM, NVDIM, NPDIM), dtype=np.float64)
SEM = np.zeros((2, NRDIM, NSDIM, NPDIM), dtype=np.float64)
VSINT = np.zeros((2, NRDIM, NPDIM), dtype=np.float64)
R = np.zeros(NRDIM, dtype=np.float64)
S = np.zeros(NSDIM, dtype=np.float64)
DELV = np.zeros(NRDIM, dtype=np.float64)
ITMAX = np.full(10, 10, dtype=np.int32)
EP = np.zeros(10, dtype=np.float64)
BBIAS = np.full(1000, 1000, dtype=np.float64)
NLOC = np.zeros(4, dtype=np.int32)
BARR = np.zeros(NVDIM1, dtype=np.float64)
PROF = np.zeros(NSDIM, dtype=np.float64)
AVBL = np.zeros(NREGDIM, dtype=np.float64)
AVBH = np.zeros(NREGDIM, dtype=np.float64)
AVBSO = np.zeros(NREGDIM, dtype=np.float64)
ESO = np.zeros(NREGDIM, dtype=np.float64)
TIP = np.zeros((NRDIM, NVDIM, NPDIM), dtype=bool)

# COMMON blocks equivalent
class SEMI:
    def __init__(self):
        self.TK = 0.0
        self.EGAP = np.zeros(NREGDIM, dtype=np.float64)
        self.ED = np.zeros(NREGDIM, dtype=np.float64)
        self.EA = np.zeros(NREGDIM, dtype=np.float64)
        self.ACB = np.zeros(NREGDIM, dtype=np.float64)
        self.AVB = np.zeros(NREGDIM, dtype=np.float64)
        self.CD = np.zeros(NREGDIM, dtype=np.float64)
        self.CA = np.zeros(NREGDIM, dtype=np.float64)
        self.IDEG = np.zeros(NREGDIM, dtype=np.int32)
        self.IINV = np.zeros(NREGDIM, dtype=np.int32)
        self.DELVB = np.zeros(NREGDIM, dtype=np.float64)
        self.EF = 0.0

SEMI = SEMI()

class PROTRU:
    def __init__(self):
        self.RAD2 = 0.0

PROTRU = PROTRU()

class SURF:
    def __init__(self):
        self.ISTK = 0
        self.TK1 = 0.0
        self.EN0 = np.zeros(NARDIM, dtype=np.float64)
        self.EN = np.zeros((NARDIM, 2), dtype=np.float64)
        self.DENS = np.zeros((NARDIM, 2), dtype=np.float64)
        self.FWHM = np.zeros((NARDIM, 2), dtype=np.float64)
        self.ECENT = np.zeros((NARDIM, 2), dtype=np.float64)
        self.CHARGE_NEUTRALITY_LEVEL = 0.0

SURF = SURF()

class CD:
    def __init__(self):
        self.EF = 0.0
        self.ESTART = 0.0
        self.DELE = 0.0
        self.NE = 0
        self.RHOBTAB = np.zeros((NREGDIM, NEDIM), dtype=np.float64)
        self.RHOSTAB = np.zeros((NARDIM, NEDIM), dtype=np.float64)

CD = CD()

class TIPPOS:
    def __init__(self):
        self.X0 = 0.0
        self.Y0 = 0.0

TIPPOS = TIPPOS()

def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        return False

def read_parameters(filename):
    with open(filename, "r", encoding="utf-8") as file:
        lines = file.readlines()
        # Filter out comment lines and empty lines
        lines = [line.split('*')[0].strip() for line in lines if not line.strip().startswith('*') and line.strip()]
        lines = [','.join(line.split()) for line in lines]  # Replace spaces with single comma
        params = []
        idx = 0
        while idx < len(lines):
            if not lines[idx]:  # Skip empty lines
                idx += 1
                continue
            if idx == 0:
                try:
                    NPARM = int(lines[idx].split(',')[0])
                except ValueError:
                    print(f"Error parsing NPARM at line {idx}: {lines[idx]}")
                    break
                idx += 1
            else:
                try:
                    param_set = {
                        "SLOPE": float(lines[idx].split(',')[0]),
                        "SEPIN": float(lines[idx+1].split(',')[0]),
                        "RAD": float(lines[idx+2].split(',')[0]),
                        "RAD2": float(lines[idx+3].split(',')[0]),
                        "CPOT": float(lines[idx+4].split(',')[0]),
                        "X0": float(lines[idx+5].split(',')[0]),
                        "Y0": float(lines[idx+6].split(',')[0]),
                        "NREG": int(lines[idx+7].split(',')[0]),
                        "DOPING": []
                    }
                    idx += 8
                    for _ in range(param_set["NREG"]):
                        doping_set = {
                            "CD": float(lines[idx].split(',')[0]),
                            "CA": float(lines[idx+1].split(',')[0]),
                            "EGAP": float(lines[idx+2].split(',')[0]),
                            "DELVB": float(lines[idx+3].split(',')[0]),
                            "ED": float(lines[idx+4].split(',')[0]),
                            "EA": float(lines[idx+5].split(',')[0]),
                            "ACB": float(lines[idx+6].split(',')[0]),
                            "AVBH": float(lines[idx+7].split(',')[0]),
                            "AVBL": float(lines[idx+8].split(',')[0]),
                            "AVB": math.exp(2.0 * math.log(math.sqrt(float(lines[idx+7].split(',')[0])**3) + math.sqrt(float(lines[idx+8].split(',')[0])**3)) / 3.0),
                            "AVBSO": float(lines[idx+9].split(',')[0]),
                            "ESO": float(lines[idx+10].split(',')[0]),
                            "IDEG": int(lines[idx+11].split(',')[0]),
                            "IINV": int(lines[idx+12].split(',')[0])
                        }
                        param_set["DOPING"].append(doping_set)
                        idx += 13
                    param_set.update({
                        "EPSIL": float(lines[idx].split(',')[0]),
                        "TEM": float(lines[idx+1].split(',')[0]),
                        "NAR": int(lines[idx+2].split(',')[0]),
                        "SURFACE_STATES": []
                    })
                    idx += 3
                    for _ in range(param_set["NAR"]):
                        surface_state_set = {
                            "DENS1": float(lines[idx].split(',')[0]),
                            "EN1": float(lines[idx+1].split(',')[0]),
                            "FWHM1": float(lines[idx+2].split(',')[0]),
                            "ECENT1": float(lines[idx+3].split(',')[0]),
                            "DENS2": float(lines[idx+4].split(',')[0]),
                            "EN2": float(lines[idx+5].split(',')[0]),
                            "FWHM2": float(lines[idx+6].split(',')[0]),
                            "ECENT2": float(lines[idx+7].split(',')[0])
                        }
                        param_set["SURFACE_STATES"].append(surface_state_set)
                        idx += 8
                    param_set.update({
                        "ISTK": int(float(lines[idx].split(',')[0])),  # Handle ISTK properly
                        "MIRROR": int(float(lines[idx+1].split(',')[0])),  # Handle MIRROR properly
                        "NR": int(float(lines[idx+2].split(',')[0])),  # Handle NR properly
                        "NV": int(float(lines[idx+3].split(',')[0])),  # Handle NV properly
                        "NS": int(float(lines[idx+4].split(',')[0])),  # Handle NS properly
                        "NP": int(float(lines[idx+5].split(',')[0])),  # Handle NP properly
                        "SIZE": float(lines[idx+6].split(',')[0]),
                        "IPMAX": int(float(lines[idx+7].split(',')[0])),  # Handle IPMAX properly
                        "ITMAX": [int(x) for x in lines[idx+8].split(',') if is_number(x.strip())],
                        "EP": [float(x) for x in lines[idx+9].split(',') if is_number(x.strip())],
                        "NE": int(float(lines[idx+10].split(',')[0])),  # Handle NE properly
                        "IWRIT": int(float(lines[idx+11].split(',')[0])),  # Handle IWRIT properly
                        "NBIAS": int(float(lines[idx+12].split(',')[0])),  # Handle NBIAS properly
                        "BBIAS": [float(x) for x in lines[idx+13].split(',') if is_number(x.strip())],
                        "NUMC": int(float(lines[idx+14].split(',')[0])),  # Handle NUMC properly
                        "DELPOT": float(lines[idx+15].split(',')[0]),
                        "PHIIN": float(lines[idx+16].split(',')[0]),
                        "CHI": float(lines[idx+17].split(',')[0]),
                        "EFTIP": float(lines[idx+18].split(',')[0]),
                        "NWK": 1,
                        "NEE": 1,
                        "EXPANI": 1,
                        "FRACZ": float(lines[idx+22].split(',')[0]),
                        "BMOD": float(lines[idx+23].split(',')[0]),
                        "ANEG": float(lines[idx+24].split(',')[0]),
                        "APOS": float(lines[idx+25].split(',')[0]),
                        "VSTART": float(lines[idx+26].split(',')[0])
                    })
                    
                    params.append(param_set)
                    idx += 27
                except Exception as e:
                    print(f"Error parsing parameters at line {idx}: {lines[idx]}")
                    print(e)
                    break
    return params

def main(params):
    semi = SEMI
    surf = SURF

    # Initialize CSAV and related variables
    CSAV = 0.0
    CSAVE = 0.0
    CSAVL = 0.0
    CSAVV = 0.0
    CSAVVE = 0.0
    CSAVVL = 0.0
    CSAVC = 0.0
    CSAVCE = 0.0
    CSAVCL = 0.0
    CURRV = 0.0
    CURRV0 = 0.0
    CURRC = 0.0
    CURRC0 = 0.0
    CURR = 0.0
    CURR0 = 0.0
    CURRVE = 0.0
    CURRVL = 0.0
    CURRCE = 0.0
    CURRCL = 0.0
    CURRE = 0.0
    CURRL = 0.0
    icomp=1
    NWK=1
    NEE=1
    EXPANI=1
    for param in params:
        SLOPE = param["SLOPE"]
        SEPIN = param["SEPIN"]
        RAD = param["RAD"]
        RAD2 = param["RAD2"]
        CPOT = param["CPOT"]
        X0 = param["X0"]
        Y0 = param["Y0"]
        NREG = param["NREG"]

        NWK = param["NWK"]
        NEE = param["NEE"]
        EXPANI = param["EXPANI"]

        print(f"RAD, SLOPE, ANGLE = {RAD:.8f} {SLOPE:.8f} {360.0 * math.atan(1.0 / SLOPE) / pi:.6f}")
        print(f"CONTACT POTENTIAL = {CPOT:.7f}")
        print(f"POSITION OF TIP = {X0:.7f} {Y0:.7f}")
        print(f"NUMBER OF DIFFERENT REGIONS OF SEMICONDUCTOR = {NREG}")

        if NREG > NREGDIM:
            raise ValueError("INPUT NUMBER OF REGIONS > NUMBER OF REGIONS IN PARAMETER STATEMENT")

        for i, doping in enumerate(param["DOPING"]):
            print(f"REGION # {i + 1}")
            print(f"DOPING = {doping['CD']:.8E} {doping['CA']:.7f}")
            print(f"BAND GAP, VB OFFSET = {doping['EGAP']:.7f} {doping['DELVB']:.7f}")

        EPSIL = param["EPSIL"]
        TEM = param["TEM"]
        semi.TK = TEM * 8.617e-5
        surf.CHARGE_NEUTRALITY_LEVEL = param["SURFACE_STATES"][0]["EN1"]
        semi.EF = surf.CHARGE_NEUTRALITY_LEVEL
        NAR = param["NAR"]

        print(f"NUMBER OF DIFFERENT AREAS OF SURFACE STATES = {NAR}")
        if NAR > NARDIM:
            raise ValueError("INPUT NUMBER OF AREAS > NUMBER OF AREAS IN PARAMETER STATEMENT")

        for i, surface_state in enumerate(param["SURFACE_STATES"]):
            print("FIRST DISTRIBUTION OF SURFACE STATES:")
            print(f"SURFACE STATE DENSITY, EN = {surface_state['DENS1']:.8E} {surface_state['EN1']:.8f}")
            print(f"FWHM, ECENT = {surface_state['FWHM1']:.8f} {surface_state['ECENT1']:.7f}")
            print("SECOND DISTRIBUTION OF SURFACE STATES:")
            print(f"SURFACE STATE DENSITY, EN = {surface_state['DENS2']:.7f} {surface_state['EN2']:.7f}")
            print(f"FWHM, ECENT = {surface_state['FWHM2']:.7f} {surface_state['ECENT2']:.7f}")

        ISTK = param["ISTK"]
        MIRROR = param["MIRROR"]
        NR = param["NR"]
        NV = param["NV"]
        NS = param["NS"]
        NP = param["NP"]
        SIZE = param["SIZE"]
        IPMAX = param["IPMAX"]
        ITMAX = param["ITMAX"]
        EP = param["EP"]
        NE = param["NE"]
        IWRIT = param["IWRIT"]
        NBIAS = param["NBIAS"]
        BBIAS = param["BBIAS"]
        NUMC = param["NUMC"]
        DELPOT = param["DELPOT"]
        PHIIN = param["PHIIN"]
        CHI = param["CHI"]
        EFTIP = param["EFTIP"]
        FRACZ = param["FRACZ"]
        BMOD = param["BMOD"]
        ANEG = param["ANEG"]
        APOS = param["APOS"]
        VSTART = param["VSTART"]
        DELRIN = param.get("DELRIN", 1.0)  # 设置默认值
        DELSIN = param.get("DELSIN", 1.0)
        DELXSI = np.zeros(NRDIM)
        IINIT = 1
        DELV = np.zeros(NRDIM)
        DELR = np.zeros(NRDIM)
        DELS = np.zeros(NSDIM)
        DELP = 2.0 * pi / NP  # 计算DELP
        VAC = np.zeros((2, NRDIM, NVDIM, NPDIM), dtype=np.float64)
        SEM = np.zeros((2, NRDIM, NSDIM, NPDIM), dtype=np.float64)
        VSINT = np.zeros((2, NRDIM, NPDIM), dtype=np.float64)
        TIP = np.zeros((NRDIM, NVDIM, NPDIM), dtype=bool)
        BARR = np.zeros(NVDIM1, dtype=np.float64)
        PROF = np.zeros(NSDIM, dtype=np.float64)
        NBARR1 = [0]
        EN0=0.1250000
        if MIRROR == 1:
            print("HORIZONTAL MIRROR PLANE ASSUMED")
            if Y0 != 0.0:
                print("*** WARNING - Y0 <> 0 WITH MIRROR PLANE; WILL SET Y0 TO ZERO")
                Y0 = 0.0
        ireg=1
        EF = effind(ireg, semi, arho, gsect, rhob, rhocb, rhoa, rhovb, rhod, fjint)
        print(f'CHARGE-NEUTRALITY LEVEL = {EN0:.7f}') 
        RHOCC = doping["CD"] * math.exp((semi.EF - semi.TK) / (8.617e-5 * TEM))  # Example calculation
        RHOVV = doping["CA"] * math.exp((semi.TK - semi.EF) / (8.617e-5 * TEM))
        SEP0 = SEPIN - abs(VSTART) * ANEG if VSTART < 0 else SEPIN - abs(VSTART) * APOS
        print(f"REGION TYPE 1, FERMI-LEVEL = {semi.EF:.7f}")
        print(f"CARRIER DENSITY IN CB, VB = {RHOCC:.8E} {RHOVV:.7f}")  # Replace with actual calculations

        for BIAS0 in BBIAS[:1]:  # Limiting to 3 points for demonstration
            SEP = SEP0 + ANEG * abs(BIAS0) if BIAS0 <= 0 else SEP0 + APOS * abs(BIAS0)
            print(f"SEPARATION = {SEP:.8f}")
            BIAS = BIAS0 + (-1 if BMOD == 0 else 1) * BMOD * math.sqrt(2.0)
            POTTIP = BIAS + CPOT
            print(f"BIAS, TIP POTENTIAL = {BIAS:.8f} {POTTIP:.8f}")
            W = 1.0e9 * math.sqrt(2.0 * EPSIL * EPSIL0 * max(1.0, abs(POTTIP)) / (abs(doping["CD"] - doping["CA"]) * 1.0e6 * E)) if (doping["CD"] - doping["CA"]) != 0 else 1.0e10
            print(f"1-D ESTIMATE OF DEPLETION WIDTH (NM) = {W:.8f}")

            ESTART = min(semi.TK, semi.TK - POTTIP, min([surface_state["EN1"] for surface_state in param["SURFACE_STATES"]]))
            EEND = max(semi.TK, semi.TK - POTTIP, max([surface_state["EN1"] for surface_state in param["SURFACE_STATES"]]))
            ETMP = EEND - ESTART
            ESTART = ESTART - 2.0 * ETMP
            EEND = EEND + 2.0 * ETMP
            DELE = (EEND - ESTART) / float(NE - 1)
            NETMP = int((semi.TK - ESTART) / DELE)
            ESTART = semi.TK - (NETMP - 0.5) * DELE
            EEND = ESTART + (NE - 1) * DELE
            print(f"ESTART,EEND,NE = {ESTART:.7f} {EEND:.7f} {NE}")

            print("COMPUTING TABLE OF BULK CHARGE DENSITIES")
            print("COMPUTING TABLE OF SURFACE CHARGE DENSITIES")
            POT0=0
            IERR=0
            ETAT, A, Z0, C, DELR, DELS, DELV, DELP, NR, NS, NV, NP, POT0, IERR,VAC, SEM,VSINT = semitip3(SEP, RAD, SLOPE, DELRIN, DELSIN, VAC, TIP, SEM, VSINT, R, S, DELV, DELR, DELXSI, DELP, 
             NRDIM, NVDIM, NSDIM, NPDIM, NR, NV, NS, NP, BIAS, IWRIT, ITMAX, EP, IPMAX, POT0, IERR, 
             IINIT, MIRROR, EPSIL, DELS)
            POT0 = 0.0  # Initialize POT0
            IERR = 0  # Initialize IERR
            
            BARR, PROF, NBARR1=potcut3(0, VAC, TIP, SEM, VSINT, NRDIM, NVDIM, NSDIM, NPDIM, NV, NS, NP, SEP+RAD2, S, DELV, POT0, BIAS, CHI, CPOT, param["DOPING"][0]["EGAP"], BARR, PROF, NBARR1, NVDIM1, NVDIM2, 0)
          
            NSP = int(round(FRACZ * NS))
            sdepth = (2 * NS * DELS[0] / pi) * math.tan(pi * NSP / (2. * NS))  # Selecting the first element of DELS
            print(f"# GRID POINTS INTO SEMICONDUCTOR USED FOR INTEGRATION = {NSP}")
            print(f"DEPTH INTO SEMICONDUCTOR USED FOR INTEGRATION = {sdepth:.4f}")
            IWRIT1 = IWRIT % 5

            CURRVE, CURRV0, CURRC, CURRC0, CURR, CURR0 = intcurr(
                0, BARR, PROF, NBARR1, NV, NS, NSP, NVDIM, NSDIM, S, SEP, BIAS, semi.EF, CHI, EFTIP,
                CPOT, semi.EGAP[0], semi.TK, AVBH[0], AVBL[0], AVBSO[0], semi.ACB[0], ESO[0], NEE, NWK,EXPANI,
                POT0, NVDIM1, NVDIM2, NSDIM2, NLOC, CURRVE, 0, CURRC, 0, CURR, 0,
                IWRIT, 0, 0, 0, 0, 0, 0, 0
            )
           
            """
             intcurr(impot, barr, prof, nbarr1, nv, ns, nsp, nvdim, nsdim, s, sep, bias, ef, chi, 
             eftip, cpot, egap, tk, avbh, avbl, avbso, acb, eso, nee, nwk, expani ,pot0, nvdim1, nvdim2, nsdim2, nloc, currv, currv0, currc, currc0, 
             curr, curr0, iwrit, icomp, cdesem, cdesurf, cdlsem, cdlsurf, cdevac, cdlvac):
            )
            """
            CURR = CURRE + CURRL
            CURRV = CURRVE + CURRVL
            CURRC = CURRCE + CURRCL

            if semi.IINV[0] in [1, 3] and CURRVE > 0:
                CURRVE = 0.0
                print('VB EXTENDED INVERSION CURRENT SET TO ZERO')

            if semi.IINV[0] in [1, 3] and CURRVL > 0:
                CURRVL = 0.0
                print('VB LOCALIZED INVERSION CURRENT SET TO ZERO')

            print('valence band current ext,loc =', CURRVE, CURRVL)

            if semi.IINV[0] in [2, 3] and CURRCE < 0:
                CURRCE = 0.0
                print('CB EXTENDED INVERSION CURRENT SET TO ZERO')

            if semi.IINV[0] in [2, 3] and CURRCL < 0:
                CURRCL = 0.0
                print('CB LOCALIZED INVERSION CURRENT SET TO ZERO')

            print('conduction band current ext,loc =', CURRCE, CURRCL)

            CURR = CURRE + CURRL
            CURRV = CURRVE + CURRVL
            CURRC = CURRCE + CURRCL

            if BMOD == -1:
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
                COND = (CURR - CSAV)
                CONDE = (CURRE - CSAVE)
                CONDL = (CURRL - CSAVL)
                CONDV = (CURRV - CSAVV)
                CONDVE = (CURRVE - CSAVVE)
                CONDVL = (CURRVL - CSAVVL)
                CONDC = (CURRC - CSAVC)
                CONDCE = (CURRCE - CSAVCE)
                CONDCL = (CURRCL - CSAVCL)
                COND = COND * math.exp(2.0 * 10.0 * (SEP - SEP0 - DELS[0])) / (2.0 * BMOD * math.sqrt(2.0))
                CONDE = CONDE * math.exp(2.0 * 10.0 * (SEP - SEP0 - DELS[0])) / (2.0 * BMOD * math.sqrt(2.0))
                CONDL = CONDL * math.exp(2.0 * 10.0 * (SEP - SEP0 - DELS[0])) / (2.0 * BMOD * math.sqrt(2.0))
                CONDV = CONDV * math.exp(2.0 * 10.0 * (SEP - SEP0 - DELS[0])) / (2.0 * BMOD * math.sqrt(2.0))
                CONDVE = CONDVE * math.exp(2.0 * 10.0 * (SEP - SEP0 - DELS[0])) / (2.0 * BMOD * math.sqrt(2.0))
                CONDVL = CONDVL * math.exp(2.0 * 10.0 * (SEP - SEP0 - DELS[0])) / (2.0 * BMOD * math.sqrt(2.0))
                CONDC = CONDC * math.exp(2.0 * 10.0 * (SEP - SEP0 - DELS[0])) / (2.0 * BMOD * math.sqrt(2.0))
                CONDCE = CONDCE * math.exp(2.0 * 10.0 * (SEP - SEP0 - DELS[0])) / (2.0 * BMOD * math.sqrt(2.0))
                CONDCL = CONDCL * math.exp(2.0 * 10.0 * (SEP - SEP0 - DELS[0])) / (2.0 * BMOD * math.sqrt(2.0))
                print(f"bias0: {BIAS0}, COND: {COND}, CONDE: {CONDE}, CONDL: {CONDL}")
                print(f"bias0: {BIAS0}, CONDV: {CONDV}, CONDVE: {CONDVE}, CONDVL: {CONDVL}")
                print(f"bias0: {BIAS0}, CONDC: {CONDC}, CONDCE: {CONDCE}, CONDCL: {CONDCL}")

            CURR = CURR * math.exp(2.0 * 10.0 * (SEP - SEP0 - DELS[0]))
            CURRE = CURRE * math.exp(2.0 * 10.0 * (SEP - SEP0 - DELS[0]))
            CURRL = CURRL * math.exp(2.0 * 10.0 * (SEP - SEP0 - DELS[0]))
            CURRV = CURRV * math.exp(2.0 * 10.0 * (SEP - SEP0 - DELS[0]))
            CURRVE = CURRVE * math.exp(2.0 * 10.0 * (SEP - SEP0 - DELS[0]))
            CURRVL = CURRVL * math.exp(2.0 * 10.0 * (SEP - SEP0 - DELS[0]))
            CURRC = CURRC * math.exp(2.0 * 10.0 * (SEP - SEP0 - DELS[0]))
            CURRCE = CURRCE * math.exp(2.0 * 10.0 * (SEP - SEP0 - DELS[0]))
            CURRCL = CURRCL * math.exp(2.0 * 10.0 * (SEP - SEP0 - DELS[0]))
            print(f"bias: {BIAS}, CURR: {CURR}, CURRE: {CURRE}, CURRL: {CURRL}")
            print(f"bias: {BIAS}, CURRV: {CURRV}, CURRVE: {CURRVE}, CURRVL: {CURRVL}")
            print(f"bias: {BIAS}, CURRC: {CURRC}, CURRCE: {CURRCE}, CURRCL: {CURRCL}")

if __name__ == "__main__":
    params = read_parameters("fort_new.9Int")
    main(params)
