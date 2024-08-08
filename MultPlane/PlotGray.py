import numpy as np


def plot_gray(dat, ncdim, nrdim, ncol, nrow, ilog, iexpan, shade, ibin, ixbin, iybin, px, py, pz, ipow, ycorrec, nfilt, ixfilt, iyfilt, irfilt, disc, nscal, zsmin, zsmax, igraycut, grang, smin, smax, shad):
    pi = np.pi
    aux1 = np.zeros((512, 512))
    aux2 = np.zeros((512, 512))
    hist = np.zeros(256)
    xcalib = 1.0
    ycorrec = 1.0

    if shade == 7:
        # Split gray scale
        nnrow = nrow
        nncol = ncol
        for j in range(nnrow):
            for i in range(nncol):
                aux1[i, j] = dat[i, j]
                if aux1[i, j] > smax:
                    smax = aux1[i, j]
                if aux1[i, j] < smin:
                    smin = aux1[i, j]
        dscal = 254.0 / (smax - smin)
        for j in range(nnrow):
            for i in range(nncol):
                aux1[i, j] = (aux1[i, j] - smin) * dscal

    elif shade == 6:
        # Local variance
        if ixbin < 1:
            ixbin = 1
        if iybin < 1:
            iybin = 1
        mcol = ixbin
        mrow = iybin
        nncol = ncol - mcol + 1
        nnrow = nrow - mrow + 1
        ioff = (mcol - 1) // 2
        joff = (mrow - 1) // 2
        sumv = 0.0

        for j in range(nnrow):
            for i in range(nncol):
                sum = 0.0
                for jj in range(mrow):
                    for ii in range(mcol):
                        sum += dat[i + ii - 1, j + jj - 2]
                aver = sum / (mrow * mcol)
                sum = 0.0
                for jj in range(mrow):
                    for ii in range(mcol):
                        sum += (dat[i + ii - 1, j + jj - 2] - aver) ** 2
                var = sum / (mrow * mcol)
                sumv += var
                aux1[i, j] = dat[i + ioff, j + joff - 1]
                if disc == -1:
                    aux2[i, j] = var
                else:
                    aux2[i, j] = 254.0 if var > disc else 0.0
                if aux1[i, j] > smax:
                    smax = aux1[i, j]
                if aux1[i, j] < smin:
                    smin = aux1[i, j]
                if aux2[i, j] > smax:
                    smax = aux2[i, j]
                if aux2[i, j] < smin:
                    smin = aux2[i, j]
        sumv /= (nnrow * nncol)

    elif shade == 5:
        # Filter derivative image
        if ixbin < 1:
            ixbin = 1
        if iybin < 1:
            iybin = 1
        mcol = ixbin
        mrow = iybin
        nncol = ncol - mcol + 1
        nnrow = nrow - mrow + 1
        ioff = mcol
        joff = mrow

        for j in range(nnrow):
            for i in range(nncol):
                aux1[i, j] = dat[i + ioff, j + joff - 1]
                sum1, sum2, sum3, sum4 = 0, 0, 0, 0
                for k in range(-((ibin - 1) // 2), (ibin // 2) + 1):
                    sum1 += dat[i + ioff + k - ibin, j + ioff]
                    sum2 += dat[i + ioff + k, j + ioff]
                    sum3 += dat[i + ioff, j + ioff + k - ibin]
                    sum4 += dat[i + ioff, j + ioff + k]
                dx = -(sum2 - sum1) / (ibin * ibin * xcalib)
                dy = -(sum4 - sum3) / (ibin * ibin * ycorrec)
                dz = 1.0
                dmag = np.sqrt(dx ** 2 + dy ** 2 + dz ** 2)
                dx /= dmag
                dy /= dmag
                dz /= dmag
                thetx = np.arccos(dz) * (2.0 / pi) * dx / \
                    np.sqrt(dx ** 2 + dy ** 2)
                thety = np.arccos(dz) * (2.0 / pi) * dy / \
                    np.sqrt(dx ** 2 + dy ** 2)
                ixp = int((1.0 + thetx) * 256.0)
                iyp = int((1.0 + thety) * 256.0)
                aux2[i, j] = 0.0
                for k in range(nfilt):
                    if ixp < (ixfilt[k] - irfilt[k]) or ixp > (ixfilt[k] + irfilt[k]):
                        continue
                    iytemp = np.sqrt(irfilt[k] ** 2 - (ixp - ixfilt[k]) ** 2)
                    if iyp < (iyfilt[k] - iytemp) or iyp > (iyfilt[k] + iytemp):
                        continue
                    aux2[i, j] = 1.0
                if aux1[i, j] > smax:
                    smax = aux1[i, j]
                if aux1[i, j] < smin:
                    smin = aux1[i, j]
                if aux2[i, j] > smax:
                    smax = aux2[i, j]
                if aux2[i, j] < smin:
                    smin = aux2[i, j]
        smin = 0.0
        smax = 1.0

    # Other shading methods can be added similarly following the same pattern.

    # Continue with the rest of the plot_gray logic.
    # Scaling and rendering would need to be implemented depending on the specific needs.
    return aux1, aux2


def dispir(ifile, img):
    xvtop = ['P', '5', ' ', '5', '1', '2', ' ',
             '5', '1', '2', ' ', '2', '5', '5', ' ']
    xvtop[2] = '\n'
    xvtop[6] = '\n'
    xvtop[10] = '\n'
    xvtop[14] = '\n'
    xvtop_header = ''.join(xvtop)

    file_names = ['img1.PGM', 'img2.PGM', 'img3.PGM', 'img4.PGM', 'img5.PGM',
                  'img6.PGM', 'img7.PGM', 'img8.PGM', 'img9.PGM']

    if ifile < 1 or ifile > 9:
        raise ValueError("ifile must be between 1 and 9")

    with open(file_names[ifile - 1], 'wb') as f:
        f.write(xvtop_header.encode('ascii'))
        for j in range(512):
            for i in range(512):
                f.write(img[i, 512 - j - 1])


# Example usage
img = np.chararray((512, 512), itemsize=1)
img[:] = 127  # Example filling
dispir(1, img)
