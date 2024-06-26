import numpy as np


def potperiod(sep, nx, ny, nz, xdel, ydel, zdel, vac, tip, sem, vsint, pvac, psem, r, s, delv, delr, dels, delp, nrdim, nvdim, nsdim, npdim, nxdim, nxdim2, nydim, nzdim, nr, nv, ns, np, iwrit, ierr):
    """
    構建平面波計算中使用的等網格間距的周期性勢能，從圓對稱坐標開始，假設X、Y和Z方向的鏡像平面。
    """
    pi = 3.141592654
    delv0 = sep / nv

    if nx > nxdim or ny > nydim or nz > nzdim:
        print('*** ERROR - INSUFFICIENT STORAGE FOR PERIODIC POTENTIAL')
        return

    # 潛勢進入半導體（及其表面）
    for k in range(1, nz + 1):
        z = (k - 1) * zdel
        if z < s[0]:
            fz = z / s[0]
            kk = 0
        else:
            for kk in range(1, ns):
                if s[kk - 1] <= z < s[kk]:
                    break
            else:
                print('*** ERROR - INTERPOLATED S NOT FOUND')
                print(z, ns, kk)
                input('press ENTER to continue')
            fz = (z - s[kk - 1]) / (s[kk] - s[kk - 1])

        for j in range(1, ny + 1):
            y = (j - 1) * ydel
            for i in range(1, 2 * nx - 1):
                x = (i - nx + 1) * xdel
                phi = np.arctan2(y, x)
                ip = int(np.round(phi / delp) + 0.5)
                ip = min(max(ip, 1), np)
                r0 = np.sqrt(x**2 + y**2)
                if i == nx - 1 and j == 1:
                    if kk == 0:
                        psem[nx - 1, 0, k - 1] = ((9.0 * (vsint[0, 0, 0] + vsint[0, 0, np - 1])
                                                   - (vsint[0, 1, 0] + vsint[0, 1, np - 1])) / 16.0) * (1.0 - fz) \
                            + ((9.0 * (sem[0, 0, 0, 0] + sem[0, 0, 0, np - 1])
                                - (sem[0, 1, 0, 0] + sem[0, 1, 0, np - 1])) / 16.0) * fz
                    else:
                        psem[nx - 1, 0, k - 1] = ((9.0 * (sem[0, 0, kk - 1, 0] + sem[0, 0, kk - 1, np - 1])
                                                   - (sem[0, 1, kk - 1, 0] + sem[0, 1, kk - 1, np - 1])) / 16.0) * (1.0 - fz) \
                            + ((9.0 * (sem[0, 0, kk, 0] + sem[0, 0, kk, np - 1])
                                - (sem[0, 1, kk, 0] + sem[0, 1, kk, np - 1])) / 16.0) * fz
                else:
                    if r0 <= r[0]:
                        f = r0 / r[0]
                        if kk == 0:
                            psem[i - 1, j - 1, k - 1] = (((9.0 * (vsint[0, 0, 0] + vsint[0, 0, np - 1])
                                                           - (vsint[0, 1, 0] + vsint[0, 1, np - 1])) / 16.0) * (1.0 - f)
                                                         + vsint[0, 0, ip - 1] * f) * (1.0 - fz) \
                                + (((9.0 * (sem[0, 0, 0, 0] + sem[0, 0, 0, np - 1])
                                     - (sem[0, 1, 0, 0] + sem[0, 1, 0, np - 1])) / 16.0) * (1.0 - f)
                                   + sem[0, 0, 0, ip - 1] * f) * fz
                        else:
                            psem[i - 1, j - 1, k - 1] = (((9.0 * (sem[0, 0, kk - 1, 0] + sem[0, 0, kk - 1, np - 1])
                                                           - (sem[0, 1, kk - 1, 0] + sem[0, 1, kk - 1, np - 1])) / 16.0) * (1.0 - f)
                                                         + sem[0, 0, kk - 1, ip - 1] * f) * (1.0 - fz) \
                                + (((9.0 * (sem[0, 0, kk, 0] + sem[0, 0, kk, np - 1])
                                     - (sem[0, 1, kk, 0] + sem[0, 1, kk, np - 1])) / 16.0) * (1.0 - f)
                                   + sem[0, 0, kk, ip - 1] * f) * fz
                    else:
                        for ii in range(1, nr + 1):
                            if r[ii - 1] <= r0 < r[ii]:
                                break
                        else:
                            print('*** ERROR - INTERPOLATED R NOT FOUND')
                            print('II, NR, R0, R =', ii, nr, r0, r[:10])
                            input('press ENTER to continue')
                        f = (r0 - r[ii - 1]) / (r[ii] - r[ii - 1])
                        if kk == 0:
                            psem[i - 1, j - 1, k - 1] = (vsint[0, ii - 1, ip - 1] * (1.0 - f) + vsint[0, ii, ip - 1] * f) * (1.0 - fz) \
                                + (sem[0, ii - 1, 0, ip - 1] * (1.0 -
                                   f) + sem[0, ii, 0, ip - 1] * f) * fz
                        else:
                            psem[i - 1, j - 1, k - 1] = (sem[0, ii - 1, kk - 1, ip - 1] * (1.0 - f) + sem[0, ii, kk - 1, ip - 1] * f) * (1.0 - fz) \
                                + (sem[0, ii - 1, kk, ip - 1] * (1.0 -
                                   f) + sem[0, ii, kk, ip - 1] * f) * fz

    # 潛勢進入真空
    for k in range(1, nv + 1):
        z = k * delv0
        for j in range(1, ny + 1):
            y = (j - 1) * ydel
            for i in range(1, 2 * nx - 1):
                x = (i - nx + 1) * xdel
                phi = np.arctan2(y, x)
                ip = int(np.round(phi / delp) + 0.5)
                ip = min(max(ip, 1), np)
                r0 = np.sqrt(x**2 + y**2)
                if i == nx - 1 and j == 1:
                    pvac[nx - 1, 0, k - 1] = (9.0 * (vac[0, 0, k - 1, 0] + vac[0, 0, k - 1, np - 1])
                                              - (vac[0, 1, k - 1, 0] + vac[0, 1, k - 1, np - 1])) / 16.0
                else:
                    if r0 <= r[0]:
                        f = r0 / r[0]
                        pvac[i - 1, j - 1, k - 1] = ((9.0 * (vac[0, 0, k - 1, 0] + vac[0, 0, k - 1, np - 1])
                                                     - (vac[0, 1, k - 1, 0] + vac[0, 1, k - 1, np - 1])) / 16.0) * (1.0 - f) \
                            + vac[0, 0, k - 1, ip - 1] * f
                    else:
                        for ii in range(1, nr + 1):
                            if r[ii - 1] <= r0 < r[ii]:
                                break
                        else:
                            print('*** ERROR - INTERPOLATED R NOT FOUND')
                            print('II, NR, R0, R =', ii, nr, r0, r[:10])
                            input('press ENTER to continue')
                        f = (r0 - r[ii - 1]) / (r[ii] - r[ii - 1])
                        delv1 = delv[ii - 1] * (1.0 - f) + delv[ii] * f
                        if z < delv1:
                            fz = z / delv1
                            kk = 0
                        else:
                            for kk in range(1, nv):
                                if kk * delv1 <= z < (kk + 1) * delv1:
                                    break
                            else:
                                print('*** ERROR - INTERPOLATED Z NOT FOUND')
                                input('press ENTER to continue')
                            fz = (z - kk * delv1) / delv1
                        if kk == 0:
                            pvac[i - 1, j - 1, k - 1] = (vsint[0, ii - 1, ip - 1] * (1.0 - f) + vsint[0, ii, ip - 1] * f) * (1.0 - fz) \
                                + (vac[0, ii - 1, 0, ip - 1] * (1.0 -
                                   f) + vac[0, ii, 0, ip - 1] * f) * fz
                        else:
                            pvac[i - 1, j - 1, k - 1] = (vac[0, ii - 1, kk - 1, ip - 1] * (1.0 - f) + vac[0, ii, kk - 1, ip - 1] * f) * (1.0 - fz) \
                                + (vac[0, ii - 1, kk, ip - 1] * (1.0 -
                                   f) + vac[0, ii, kk, ip - 1] * f) * fz

    return
