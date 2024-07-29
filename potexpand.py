import numpy as np

def potexpand(impot, sep, nv, pot0p, s, ns, nsdim, barr, nbarr1, nvdim1, nvdim2, prof, prof2, nsdim2, s2, vacstep, semstep, jsem, nexsem, iwrit):
    kappa, lambda_val = 0.0, 0.0

    # Check for zero or undefined values
    if nv == 0 or vacstep == 0:
        raise ValueError("nv and vacstep must be non-zero and defined")
    '''
    if np.isnan(sep) or np.isnan(nv) or np.isnan(vacstep):
        print(f"Invalid values detected: sep={sep}, nv={nv}, vacstep={vacstep}")
        raise ValueError("sep, nv, and vacstep must not be NaN")
    '''
    """
    print(f"sep: {sep}, nv: {nv}, vacstep: {vacstep}")
    """ #不需要輸出 但可以留下
    # Expand vacuum barrier
    vacstep=1.0
    nexpan = max(1, int((sep / nv) / vacstep))
    if impot == 1:
        nexpan *= 10
    nexvac = nexpan
    if iwrit > 1:
        print('Expansion factor for barrier =', nexpan)
    nbarr2 = nexpan * (nbarr1 - 1) + 1

    # Ensure barr2 has the correct size
    barr2 = np.zeros(nbarr2)
    barr2[nbarr2 - 1] = barr[nbarr1 - 1]
    for j in range(nbarr1 - 2, -1, -1):
        b2 = barr[j + 1]
        b1 = barr[j]
        for k in range(nexpan - 1, -1, -1):
            barr2[j * nexpan + k + 1] = (b2 * float(k) + b1 * float(nexpan - k)) / nexpan
    if iwrit >= 3:
        for i in range(nbarr2):
            print(-(i) * sep / float(nbarr2 - 1), barr2[i])
    if iwrit > 1:
        print('Number of expanded points in vacuum =', nbarr2)
    lambda_val = 3.81**2 * 0.1 * np.log(2.0) / (2.0 * 2.0 * sep)
    if impot == 1:
        for j in range(1, nbarr2 - 1):
            barr2[j] -= 1.15 * lambda_val * (nbarr2 - 1.0)**2 / ((j) * (float(nbarr2) - j))
    if iwrit >= 3:
        for i in range(nbarr2):
            print(-(i) * sep / float(nbarr2 - 1), barr2[i])

    # Expand the potential profile in semiconductor
    for j in range(ns):
        nexsem[j] = 0
    nexpan = max(1, int(2.0 * s[0] / semstep))
    if iwrit > 1:
        print('Initial expansion factor for semiconductor =', nexpan)
    kk = 10  # Start kk from 10
    for j in range(ns):
        if j == 0:
            nexpan = max(1, int(s[0] / semstep))
        else:
            nexpan = max(1, int((s[j] - s[j - 1]) / semstep))
        if nexpan % 2 == 0:
            nexpan += 1

        for k in range(nexpan):
            if kk > nsdim2:
                # Resize the arrays dynamically
                nsdim2 = kk + 100
                jsem = np.resize(jsem, nsdim2)
                prof2 = np.resize(prof2, nsdim2)
                s2 = np.resize(s2, nsdim2)
                nexsem = np.resize(nexsem, nsdim2)
                print(f"Resizing arrays to accommodate {nsdim2} elements.")
            '''
            if j == 0:
                jsem[kk - 1] = j + 1
            else:
                if k <= (nexpan // 2):
                    jsem[kk - 1] = j
                else:
                    jsem[kk - 1] = j + 1
            nexsem[jsem[kk - 1] - 1] += 1
            
            if j == 0:
                prof2[kk - 1] = ((nexpan - k) * pot0p + k * prof[j]) / float(nexpan)
            else:
                prof2[kk - 1] = ((nexpan - k) * prof[j - 1] + k * prof[j]) / float(nexpan)
            if j == 0:
                s2[kk - 1] = ((nexpan - k) * 0.0 + k * s[j]) / float(nexpan)
            else:
                s2[kk - 1] = ((nexpan - k) * s[j - 1] + k * s[j]) / float(nexpan)
            kk += 1
            '''
    ns2 = kk - 1  # Adjust ns2 to the correct size
    if iwrit > 1:
        print('Number of expanded points in semiconductor =', ns2)

    return barr2, prof2, s2, jsem, nexsem, nexvac, ns2

# Example call to potexpand with debugging print statements
sep = 1.0
nv = 1.0
vacstep = 0.1
impot = 1
pot0p = 1.0
s = np.array([1.0, 2.0, 3.0])
ns = len(s)
nsdim = ns
barr = np.array([0.5, 0.7, 0.9])
nbarr1 = len(barr)
nvdim1 = 100
nvdim2 = 200
prof = np.array([0.1, 0.2, 0.3])
prof2 = np.zeros(100)
nsdim2 = 100
s2 = np.zeros(100)
semstep = 0.1
jsem = np.zeros(100, dtype=int)
nexsem = np.zeros(100, dtype=int)
iwrit = 3

# Correct the function call to match the definition
barr2, prof2, s2, jsem, nexsem, nexvac, ns2 = potexpand(impot, sep, nv, pot0p, s, ns, nsdim, barr, nbarr1, nvdim1, nvdim2, prof, prof2, nsdim2, s2, vacstep, semstep, jsem, nexsem, iwrit)

