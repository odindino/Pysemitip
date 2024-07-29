import numpy as np
from math import atan2

def potperiod(SEP, NX, NY, NZ, XDEL, YDEL, ZDEL, VAC, TIP, SEM, VSINT, PVAC, PSEM, R, S, DELV, DELR, DELS, DELP, NRDIM, NVDIM, NSDIM, NPDIM, NXDIM, NXDIM2, NYDIM, NZDIM, NR, NV, NS, NP, IWRIT, IERR):
    """
    構建平面波計算中使用的等網格間距的周期性勢能，從圓對稱坐標開始，假設X、Y和Z方向的鏡像平面。
    """
    PI = 3.141592654
    DELV0 = SEP / NV

    if NX > NXDIM or NY > NYDIM or NZ > NZDIM:
        print('*** ERROR - INSUFFICIENT STORAGE FOR PERIODIC POTENTIAL')
        return

    # 潛勢進入半導體（及其表面）
    for K in range(1, NZ + 1):
        Z = (K - 1) * ZDEL
        if Z < S[0]:
            FZ = Z / S[0]
            KK = 0
        else:
            for KK in range(1, NS):
                if S[KK - 1] <= Z < S[KK]:
                    break
            else:
                print('*** ERROR - INTERPOLATED S NOT FOUND')
                print(Z, NS, KK)
                input('press ENTER to continue')
            FZ = (Z - S[KK - 1]) / (S[KK] - S[KK - 1])

        for J in range(1, NY + 1):
            Y = (J - 1) * YDEL
            for I in range(1, 2 * NX - 1):
                X = (I - NX + 1) * XDEL
                PHI = atan2(Y, X)  # 修正此處
                IP = int(np.round(PHI / DELP) + 0.5)
                IP = min(max(IP, 1), NP)
                R0 = np.sqrt(X**2 + Y**2)
                if I == NX - 1 and J == 1:
                    if KK == 0:
                        PSEM[NX - 1, 0, K - 1] = ((9.0 * (VSINT[0, 0, 0] + VSINT[0, 0, NP - 1])
                                                   - (VSINT[0, 1, 0] + VSINT[0, 1, NP - 1])) / 16.0) * (1.0 - FZ) \
                            + ((9.0 * (SEM[0, 0, 0, 0] + SEM[0, 0, 0, NP - 1])
                                - (SEM[0, 1, 0, 0] + SEM[0, 1, 0, NP - 1])) / 16.0) * FZ
                    else:
                        PSEM[NX - 1, 0, K - 1] = ((9.0 * (SEM[0, 0, KK - 1, 0] + SEM[0, 0, KK - 1, NP - 1])
                                                   - (SEM[0, 1, KK - 1, 0] + SEM[0, 1, KK - 1, NP - 1])) / 16.0) * (1.0 - FZ) \
                            + ((9.0 * (SEM[0, 0, KK, 0] + SEM[0, 0, KK, NP - 1])
                                - (SEM[0, 1, KK, 0] + SEM[0, 1, KK, NP - 1])) / 16.0) * FZ
                else:
                    if R0 <= R[0]:
                        F = R0 / R[0]
                        if KK == 0:
                            PSEM[I - 1, J - 1, K - 1] = (((9.0 * (VSINT[0, 0, 0] + VSINT[0, 0, NP - 1])
                                                           - (VSINT[0, 1, 0] + VSINT[0, 1, NP - 1])) / 16.0) * (1.0 - F)
                                                         + VSINT[0, 0, IP - 1] * F) * (1.0 - FZ) \
                                + (((9.0 * (SEM[0, 0, 0, 0] + SEM[0, 0, 0, NP - 1])
                                     - (SEM[0, 1, 0, 0] + SEM[0, 1, 0, NP - 1])) / 16.0) * (1.0 - F)
                                   + SEM[0, 0, 0, IP - 1] * F) * FZ
                        else:
                            PSEM[I - 1, J - 1, K - 1] = (((9.0 * (SEM[0, 0, KK - 1, 0] + SEM[0, 0, KK - 1, NP - 1])
                                                           - (SEM[0, 1, KK - 1, 0] + SEM[0, 1, KK - 1, NP - 1])) / 16.0) * (1.0 - F)
                                                         + SEM[0, 0, KK - 1, IP - 1] * F) * (1.0 - FZ) \
                                + (((9.0 * (SEM[0, 0, KK, 0] + SEM[0, 0, KK, NP - 1])
                                     - (SEM[0, 1, KK, 0] + SEM[0, 1, KK, NP - 1])) / 16.0) * (1.0 - F)
                                   + SEM[0, 0, KK, IP - 1] * F) * FZ
                    else:
                        for II in range(1, NR + 1):
                            if R[II - 1] <= R0 < R[II]:
                                break
                        else:
                            print('*** ERROR - INTERPOLATED R NOT FOUND')
                            print('II, NR, R0, R =', II, NR, R0, R[:10])
                            input('press ENTER to continue')
                        F = (R0 - R[II - 1]) / (R[II] - R[II - 1])
                        if KK == 0:
                            PSEM[I - 1, J - 1, K - 1] = (VSINT[0, II - 1, IP - 1] * (1.0 - F) + VSINT[0, II, IP - 1] * F) * (1.0 - FZ) \
                                + (SEM[0, II - 1, 0, IP - 1] * (1.0 -
                                   F) + SEM[0, II, 0, IP - 1] * F) * FZ
                        else:
                            PSEM[I - 1, J - 1, K - 1] = (SEM[0, II - 1, KK - 1, IP - 1] * (1.0 - F) + SEM[0, II, KK - 1, IP - 1] * F) * (1.0 - FZ) \
                                + (SEM[0, II - 1, KK, IP - 1] * (1.0 -
                                   F) + SEM[0, II, KK, IP - 1] * F) * FZ

    # 潛勢進入真空
    for K in range(1, NV + 1):
        Z = K * DELV0
        for J in range(1, NY + 1):
            Y = (J - 1) * YDEL
            for I in range(1, 2 * NX - 1):
                X = (I - NX + 1) * XDEL
                PHI = atan2(Y, X)  # 修正此處
                IP = int(np.round(PHI / DELP) + 0.5)
                IP = min(max(IP, 1), NP)
                R0 = np.sqrt(X**2 + Y**2)
                if I == NX - 1 and J == 1:
                    PVAC[NX - 1, 0, K - 1] = (9.0 * (VAC[0, 0, K - 1, 0] + VAC[0, 0, K - 1, NP - 1])
                                              - (VAC[0, 1, K - 1, 0] + VAC[0, 1, K - 1, NP - 1])) / 16.0
                else:
                    if R0 <= R[0]:
                        F = R0 / R[0]
                        PVAC[I - 1, J - 1, K - 1] = ((9.0 * (VAC[0, 0, K - 1, 0] + VAC[0, 0, K - 1, NP - 1])
                                                     - (VAC[0, 1, K - 1, 0] + VAC[0, 1, K - 1, NP - 1])) / 16.0) * (1.0 - F) \
                            + VAC[0, 0, K - 1, IP - 1] * F
                    else:
                        for II in range(1, NR + 1):
                            if R[II - 1] <= R0 < R[II]:
                                break
                        else:
                            print('*** ERROR - INTERPOLATED R NOT FOUND')
                            print('II, NR, R0, R =', II, NR, R0, R[:10])
                            input('press ENTER to continue')
                        F = (R0 - R[II - 1]) / (R[II] - R[II - 1])
                        DELV1 = DELV[II - 1] * (1.0 - F) + DELV[II] * F
                        if Z < DELV1:
                            FZ = Z / DELV1
                            KK = 0
                        else:
                            for KK in range(1, NV):
                                if KK * DELV1 <= Z < (KK + 1) * DELV1:
                                    break
                            else:
                                print('*** ERROR - INTERPOLATED Z NOT FOUND')
                                input('press ENTER to continue')
                            FZ = (Z - KK * DELV1) / DELV1
                        if KK == 0:
                            PVAC[I - 1, J - 1, K - 1] = (VSINT[0, II - 1, IP - 1] * (1.0 - F) + VSINT[0, II, IP - 1] * F) * (1.0 - FZ) \
                                + (VAC[0, II - 1, 0, IP - 1] * (1.0 -
                                   F) + VAC[0, II, 0, IP - 1] * F) * FZ
                        else:
                            PVAC[I - 1, J - 1, K - 1] = (VAC[0, II - 1, KK - 1, IP - 1] * (1.0 - F) + VAC[0, II, KK - 1, IP - 1] * F) * (1.0 - FZ) \
                                + (VAC[0, II - 1, KK, IP - 1] * (1.0 -
                                   F) + VAC[0, II, KK, IP - 1] * F) * FZ

    return
