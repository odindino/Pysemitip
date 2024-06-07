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
    # 初始化 CC 和 SEMI 字典
    CC = {
        "CD": np.zeros(NREGDIM),
        "CA": np.zeros(NREGDIM)
    }
    SEMI = {
        "EGAP": np.zeros(NREGDIM),
        "DELVB": np.zeros(NREGDIM),
        "ED": np.zeros(NREGDIM),
        "EA": np.zeros(NREGDIM),
        "ACB": np.zeros(NREGDIM),
        "AVB": np.zeros(NREGDIM),
        "IDEG": np.zeros(NREGDIM, dtype=int),
        "IINV": np.zeros(NREGDIM, dtype=int)
    }

    for param in parameters:
        SLOPE = param["SLOPE"]
        SEPIN = param["SEPIN"]
        RAD = param["RAD"]
        RAD2 = param["RAD2"]
        CPot = param["CPot"]
        X0 = param["X0"]
        Y0 = param["Y0"]
        NREG = param["NREG"]

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
                # 呼叫 ENFIND 子程序
                EN0[IAR] = (SURF["EN"][IAR, 0] + SURF["EN"]
                            [IAR, 1]) / 2  # 模擬計算

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
        DELS1 = param["DELS1"]
        SEP0 = param["SEP0"]


# Example usage
parameters = read_fort9("fort_new.9")
process_parameters(parameters)
