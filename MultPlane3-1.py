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

CD = {
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

# Data input from fort.9 file


def read_fort9(filename="fort.9"):
    with open(filename, 'r') as file:
        data = file.readlines()

    data = [line.split()[0] for line in data if line.strip()
            and not line.lstrip().startswith('C')]

    # print(data)
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

        parameters.append(param)

    return parameters


# Example usage
parameters = read_fort9("fort.9")
print(parameters)
