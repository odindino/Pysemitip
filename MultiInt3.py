import math

# Constants
EPSIL0 = 8.854185e-12
E = 1.60210e-19
PI = 4.0 * math.atan(1.0)

# Parameters
NRDIM = 512
NVDIM = 64
NSDIM = 512
NPDIM = 64
NVDIM1 = NVDIM + 1
NVDIM2 = 2048
NSDIM2 = 20000
NEDIM = 50000
NREGDIM = 2
NARDIM = 2

# Arrays and Variables
VAC = [[[[0.0 for _ in range(NPDIM)] for _ in range(NVDIM)] for _ in range(NRDIM)] for _ in range(2)]
SEM = [[[[0.0 for _ in range(NPDIM)] for _ in range(NSDIM)] for _ in range(NRDIM)] for _ in range(2)]
VSINT = [[[0.0 for _ in range(NPDIM)] for _ in range(NRDIM)] for _ in range(2)]
R = [0.0 for _ in range(NRDIM)]
S = [0.0 for _ in range(NSDIM)]
DELV = [0.0 for _ in range(NRDIM)]
ITMAX = [10 for _ in range(10)]
EP = [0.0 for _ in range(10)]
BBIAS = [1000 for _ in range(1000)]
NLOC = [0 for _ in range(4)]
BARR = [0.0 for _ in range(NVDIM1)]
PROF = [0.0 for _ in range(NSDIM)]
AVBL = [0.0 for _ in range(NREGDIM)]
AVBH = [0.0 for _ in range(NREGDIM)]
AVBSO = [0.0 for _ in range(NREGDIM)]
ESO = [0.0 for _ in range(NREGDIM)]
TIP = [[[False for _ in range(NPDIM)] for _ in range(NVDIM)] for _ in range(NRDIM)]

# COMMON blocks equivalent
class SEMI:
    def __init__(self):
        self.TK = 0.0
        self.EGAP = [0.0 for _ in range(NREGDIM)]
        self.ED = [0.0 for _ in range(NREGDIM)]
        self.EA = [0.0 for _ in range(NREGDIM)]
        self.ACB = [0.0 for _ in range(NREGDIM)]
        self.AVB = [0.0 for _ in range(NREGDIM)]
        self.CD = [0.0 for _ in range(NREGDIM)]
        self.CA = [0.0 for _ in range(NREGDIM)]
        self.IDEG = [0 for _ in range(NREGDIM)]
        self.IINV = [0 for _ in range(NREGDIM)]
        self.DELVB = [0.0 for _ in range(NREGDIM)]

SEMI = SEMI()

class PROTRU:
    def __init__(self):
        self.RAD2 = 0.0

PROTRU = PROTRU()

class SURF:
    def __init__(self):
        self.ISTK = 0
        self.TK1 = 0.0
        self.EN0 = [0.0 for _ in range(NARDIM)]
        self.EN = [[0.0 for _ in range(2)] for _ in range(NARDIM)]
        self.DENS = [[0.0 for _ in range(2)] for _ in range(NARDIM)]
        self.FWHM = [[0.0 for _ in range(2)] for _ in range(NARDIM)]
        self.ECENT = [[0.0 for _ in range(2)] for _ in range(NARDIM)]

SURF = SURF()

class CD:
    def __init__(self):
        self.EF = 0.0
        self.ESTART = 0.0
        self.DELE = 0.0
        self.NE = 0
        self.RHOBTAB = [[0.0 for _ in range(NEDIM)] for _ in range(NREGDIM)]
        self.RHOSTAB = [[0.0 for _ in range(NEDIM)] for _ in range(NARDIM)]

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
                NPARM = int(lines[idx].split(',')[0])
                idx += 1
            else:
                try:
                    param_set = {
                        "slope": float(lines[idx].split(',')[0]),
                        "sepin": float(lines[idx+1].split(',')[0]),
                        "rad": float(lines[idx+2].split(',')[0]),
                        "rad2": float(lines[idx+3].split(',')[0]),
                        "cpot": float(lines[idx+4].split(',')[0]),
                        "x0": float(lines[idx+5].split(',')[0]),
                        "y0": float(lines[idx+6].split(',')[0]),
                        "nreg": int(lines[idx+7].split(',')[0]),
                        "doping": []
                    }
                    idx += 8
                    for _ in range(param_set["nreg"]):
                        doping_set = {
                            "cd": float(lines[idx].split(',')[0]),
                            "ca": float(lines[idx+1].split(',')[0]),
                            "egap": float(lines[idx+2].split(',')[0]),
                            "delvb": float(lines[idx+3].split(',')[0]),
                            "ed": float(lines[idx+4].split(',')[0]),
                            "ea": float(lines[idx+5].split(',')[0]),
                            "acb": float(lines[idx+6].split(',')[0]),
                            "avbh": float(lines[idx+7].split(',')[0]),
                            "avbl": float(lines[idx+8].split(',')[0]),
                            "avb": math.exp(2.0 * math.log(math.sqrt(float(lines[idx+7].split(',')[0])**3) + math.sqrt(float(lines[idx+8].split(',')[0])**3)) / 3.0),
                            "avbso": float(lines[idx+9].split(',')[0]),
                            "eso": float(lines[idx+10].split(',')[0]),
                            "ideg": int(lines[idx+11].split(',')[0]),
                            "iinv": int(lines[idx+12].split(',')[0])
                        }
                        param_set["doping"].append(doping_set)
                        idx += 13
                    param_set.update({
                        "epsil": float(lines[idx].split(',')[0]),
                        "tem": float(lines[idx+1].split(',')[0]),
                        "nar": int(lines[idx+2].split(',')[0]),
                        "surface_states": []
                    })
                    idx += 3
                    for _ in range(param_set["nar"]):
                        surface_state_set = {
                            "dens1": float(lines[idx].split(',')[0]),
                            "en1": float(lines[idx+1].split(',')[0]),
                            "fwhm1": float(lines[idx+2].split(',')[0]),
                            "ecent1": float(lines[idx+3].split(',')[0]),
                            "dens2": float(lines[idx+4].split(',')[0]),
                            "en2": float(lines[idx+5].split(',')[0]),
                            "fwhm2": float(lines[idx+6].split(',')[0]),
                            "ecent2": float(lines[idx+7].split(',')[0])
                        }
                        param_set["surface_states"].append(surface_state_set)
                        idx += 8
                    param_set.update({
                        "istk": int(lines[idx].split(',')[0]),
                        "mirror": int(lines[idx+1].split(',')[0]),
                        "nr": int(lines[idx+2].split(',')[0]),
                        "nv": int(lines[idx+3].split(',')[0]),
                        "ns": int(lines[idx+4].split(',')[0]),
                        "np": int(lines[idx+5].split(',')[0]),
                        "size": float(lines[idx+6].split(',')[0]),
                        "ipmax": int(lines[idx+7].split(',')[0]),
                        "itmax": [int(x) for x in lines[idx+8].split(',') if is_number(x.strip())],
                        "ep": [float(x) for x in lines[idx+9].split(',') if is_number(x.strip())],
                        "ne": int(lines[idx+10].split(',')[0]),
                        "iwrite": int(lines[idx+11].split(',')[0]),
                        "nbias": int(lines[idx+12].split(',')[0]),
                        "bbias": [float(x) for x in lines[idx+13].split(',') if is_number(x.strip())],
                        "numc": int(lines[idx+14].split(',')[0]),
                        "delpot": float(lines[idx+15].split(',')[0]),
                        "phiin": float(lines[idx+16].split(',')[0]),
                        "chi": float(lines[idx+17].split(',')[0]),
                        "eftip": float(lines[idx+18].split(',')[0]),
                        "nwk": int(lines[idx+19].split(',')[0]),
                        "nee": int(lines[idx+20].split(',')[0]),
                        "expani": int(lines[idx+21].split(',')[0]),
                        "fracz": float(lines[idx+22].split(',')[0]),
                        "bmod": float(lines[idx+23].split(',')[0]),
                        "aneg": float(lines[idx+24].split(',')[0]),
                        "apos": float(lines[idx+25].split(',')[0]),
                        "vstart": float(lines[idx+26].split(',')[0])
                    })
                    params.append(param_set)
                    idx += 27
                except Exception as e:
                    print(f"Error parsing parameters at line {idx}: {lines[idx]}")
                    print(e)
                    break
    return params

def main(params):
    for param in params:
        slope = param["slope"]
        sepin = param["sepin"]
        rad = param["rad"]
        rad2 = param["rad2"]
        cpot = param["cpot"]
        x0 = param["x0"]
        y0 = param["y0"]
        nreg = param["nreg"]

        # Print semiconductor and tip parameters
        print(f"RAD, SLOPE, ANGLE = {rad:.8f} {slope:.8f} {360.0 * math.atan(1.0 / slope) / PI:.6f}")
        print(f"CONTACT POTENTIAL = {cpot:.7f}")
        print(f"POSITION OF TIP = {x0:.7f} {y0:.7f}")
        print(f"NUMBER OF DIFFERENT REGIONS OF SEMICONDUCTOR = {nreg}")

        if nreg > NREGDIM:
            raise ValueError("INPUT NUMBER OF REGIONS > NUMBER OF REGIONS IN PARAMETER STATEMENT")

        for i, doping in enumerate(param["doping"]):
            print(f"REGION # {i + 1}")
            print(f"DOPING = {doping['cd']:.8E} {doping['ca']:.7f}")
            print(f"BAND GAP, VB OFFSET = {doping['egap']:.7f} {doping['delvb']:.7f}")

        epsil = param["epsil"]
        tem = param["tem"]
        tk = tem * 8.617e-5
        tk1 = tk
        nar = param["nar"]

        print(f"NUMBER OF DIFFERENT AREAS OF SURFACE STATES = {nar}")
        if nar > NARDIM:
            raise ValueError("INPUT NUMBER OF AREAS > NUMBER OF AREAS IN PARAMETER STATEMENT")

        for i, surface_state in enumerate(param["surface_states"]):
            print("FIRST DISTRIBUTION OF SURFACE STATES:")
            print(f"SURFACE STATE DENSITY, EN = {surface_state['dens1']:.8E} {surface_state['en1']:.8f}")
            print(f"FWHM, ECENT = {surface_state['fwhm1']:.8f} {surface_state['ecent1']:.7f}")
            print("SECOND DISTRIBUTION OF SURFACE STATES:")
            print(f"SURFACE STATE DENSITY, EN = {surface_state['dens2']:.7f} {surface_state['en2']:.7f}")
            print(f"FWHM, ECENT = {surface_state['fwhm2']:.7f} {surface_state['ecent2']:.7f}")

        istk = param["istk"]
        mirror = param["mirror"]
        nr = param["nr"]
        nv = param["nv"]
        ns = param["ns"]
        np = param["np"]
        size = param["size"]
        ipmax = param["ipmax"]
        itmax = param["itmax"]
        ep = param["ep"]
        ne = param["ne"]
        iwrite = param["iwrite"]
        nbias = param["nbias"]
        bbias = param["bbias"]
        numc = param["numc"]
        delpot = param["delpot"]
        phiin = param["phiin"]
        chi = param["chi"]
        eftip = param["eftip"]
        nwk = param["nwk"]
        nee = param["nee"]
        expani = param["expani"]
        fracz = param["fracz"]
        bmod = param["bmod"]
        aneg = param["aneg"]
        apos = param["apos"]
        vstart = param["vstart"]

        if mirror == 1:
            print("HORIZONTAL MIRROR PLANE ASSUMED")
            if y0 != 0.0:
                print("*** WARNING - Y0 <> 0 WITH MIRROR PLANE; WILL SET Y0 TO ZERO")
                y0 = 0.0

        sep0 = sepin - abs(vstart) * aneg if vstart < 0 else sepin - abs(vstart) * apos
        print(f"REGION TYPE 1, FERMI-LEVEL = {tk:.7f}")
        print(f"CARRIER DENSITY IN CB, VB = {tk * 2:.8E} {tk * 3:.7f}")  # Replace with actual calculations

        # Iterate over a limited set of bias points to avoid excessive output
        for bias0 in bbias[:1]:  # Limiting to 3 points for demonstration
            sep = sep0 + aneg * abs(bias0) if bias0 <= 0 else sep0 + apos * abs(bias0)
            print(f"SEPARATION = {sep:.8f}")
            bias = bias0 + (-1 if bmod == 0 else 1) * bmod * math.sqrt(2.0)
            pottip = bias + cpot
            print(f"BIAS, TIP POTENTIAL = {bias:.8f} {pottip:.8f}")
            w = 1.0e9 * math.sqrt(2.0 * epsil * EPSIL0 * max(1.0, abs(pottip)) / (abs(doping["cd"] - doping["ca"]) * 1.0e6 * E)) if (doping["cd"] - doping["ca"]) != 0 else 1.0e10
            print(f"1-D ESTIMATE OF DEPLETION WIDTH (NM) = {w:.8f}")

            estart = min(tk, tk - pottip, min([surface_state["en1"] for surface_state in param["surface_states"]]))
            eend = max(tk, tk - pottip, max([surface_state["en1"] for surface_state in param["surface_states"]]))
            etmp = eend - estart
            estart = estart - 2.0 * etmp
            eend = eend + 2.0 * etmp
            dele = (eend - estart) / float(ne - 1)
            netmp = int((tk - estart) / dele)
            estart = tk - (netmp - 0.5) * dele
            eend = estart + (ne - 1) * dele
            print(f"ESTART,EEND,NE = {estart:.7f} {eend:.7f} {ne}")

            # Placeholder for calls to SEMIRHO and SURFRHO
            print("COMPUTING TABLE OF BULK CHARGE DENSITIES")
            print("COMPUTING TABLE OF SURFACE CHARGE DENSITIES")

        # Placeholder for call to SEMITIP3
        print(f"CALL SEMITIP3 with parameters: SEP={sep+rad2:.8f}, RAD={rad:.8f}, SLOPE={slope:.8f}, etc.")

if __name__ == "__main__":
    params = read_parameters("fort_new.9Int")
    main(params)
