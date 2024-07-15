import numpy as np

def contr3(eta1, vac, tip, sem, vsint, r, s, delv, nrdim, nvdim, nsdim, npdim, nr, nv, ns, np, numc, delpot, mirror, kplot1, kplot2):
    """
    Translated CONTR3 function from Fortran to Python.
    Draw contours of the potential distributions (3D case).
    """
    nrs = nr // 500
    if nrs == 0:
        nrs = 1

    # Draw TIP
    for i in range(-nr + 1, nr + 1, nrs):
        if i == 0:
            continue
        if i > 0:
            ii = i
            rsav = r[ii - 1]
        else:
            ii = -i
            rsav = -r[ii - 1]
        
        for j in range(1, nv + 1):
            if not tip[ii - 1, j - 1, 0]:
                continue
            print(rsav * np.sqrt(1 - (eta1 * j / nv) ** 2), -j * delv[ii - 1])
            break

    # Search for min, max points in potential
    pmin = 1e10
    pmax = -1e10
    for i in range(1, nr + 1):
        for k in range(1, np + 1, np - 1):
            for j in range(1, nv + 1):
                if pmin > vac[0, i - 1, j - 1, k - 1]:
                    pmin = vac[0, i - 1, j - 1, k - 1]
                if pmax < vac[0, i - 1, j - 1, k - 1]:
                    pmax = vac[0, i - 1, j - 1, k - 1]
            if pmin > vsint[0, i - 1, k - 1]:
                pmin = vsint[0, i - 1, k - 1]
            if pmax < vsint[0, i - 1, k - 1]:
                pmax = vsint[0, i - 1, k - 1]
            for j in range(1, ns + 1):
                if pmin > sem[0, i - 1, j - 1, k - 1]:
                    pmin = sem[0, i - 1, j - 1, k - 1]
                if pmax < sem[0, i - 1, j - 1, k - 1]:
                    pmax = sem[0, i - 1, j - 1, k - 1]

    print('MIN, MAX POTENTIAL VALUES =', pmin, pmax)

    # Draw contours
    if delpot == 0:
        delpot = (pmax - pmin) / (numc + 1)
        print('CONTOUR SPACING =', delpot)

    for i in range(-nr + 1, nr + 1, nrs):
        if i == 0:
            continue
        if i > 0:
            ii = i
            rsav = r[ii - 1]
            kp = kplot1
        else:
            ii = -i
            rsav = -r[ii - 1]
            kp = kplot2

        kdone = np.zeros(numc, dtype=bool)

        for k in range(1, numc + 1):
            for j in range(ns, 0, -1):
                p = k * delpot + pmin
                if j == 1:
                    if (sem[0, ii - 1, j - 1, kp - 1] >= p and vsint[0, ii - 1, kp - 1] <= p) or (sem[0, ii - 1, j - 1, kp - 1] <= p and vsint[0, ii - 1, kp - 1] >= p):
                        print(rsav, s[j - 1])
                        kdone[k - 1] = True
                        break
                else:
                    if (sem[0, ii - 1, j - 1, kp - 1] >= p and sem[0, ii - 1, j - 2, kp - 1] <= p) or (sem[0, ii - 1, j - 1, kp - 1] <= p and sem[0, ii - 1, j - 2, kp - 1] >= p):
                        print(rsav, s[j - 1])
                        kdone[k - 1] = True
                        break

            for k in range(1, numc + 1):
                p = k * delpot + pmin
                if (vsint[0, ii - 1, kp - 1] >= p and vac[0, ii - 1, 0, kp - 1] <= p) or (vsint[0, ii - 1, kp - 1] <= p and vac[0, ii - 1, 0, kp - 1] >= p):
                    if not kdone[k - 1]:
                        print(rsav, 0)
                        kdone[k - 1] = True

            for k in range(1, numc + 1):
                for j in range(1, nv):
                    if tip[ii - 1, j - 1, kp - 1]:
                        continue
                    p = k * delpot + pmin
                    if (vac[0, ii - 1, j - 1, kp - 1] >= p and vac[0, ii - 1, j, kp - 1] <= p) or (vac[0, ii - 1, j - 1, kp - 1] <= p and vac[0, ii - 1, j, kp - 1] >= p):
                        if not kdone[k - 1]:
                            print(rsav * np.sqrt(1 - (eta1 * j / nv) ** 2), -j * delv[ii - 1])
                            kdone[k - 1] = True
                            break

# Sample data preparation
nr = 10
nv = 10
ns = 10
npdim = 5
nrdim = 10
nvdim = 10
nsdim = 10
numc = 5
delpot = 0.5
mirror = 0
kplot1 = 1
kplot2 = 2
eta1 = 1.0
bias = 1.0

vac = np.random.rand(2, nrdim, nvdim, npdim) * 10
tip = np.random.choice([False, True], (nrdim, nvdim, npdim))
sem = np.random.rand(2, nrdim, nsdim, npdim) * 10
vsint = np.random.rand(2, nrdim, npdim) * 10
r = np.linspace(0, 10, nrdim)
s = np.linspace(0, 10, nsdim)
delv = np.linspace(0.1, 1, nrdim)

# Run the translated contr3 function
contr3(eta1, vac, tip, sem, vsint, r, s, delv, nrdim, nvdim, nsdim, npdim, nr, nv, ns, npdim, numc, delpot, mirror, kplot1, kplot2)
