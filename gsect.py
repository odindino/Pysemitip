def gsect(f, xmin, xmax, ep, *args):
    GS = 0.3819660
    if xmax == xmin or ep == 0:
        return (xmin + xmax) / 2
    if xmax < xmin:
        xmin, xmax = xmax, xmin

    delx = xmax - xmin
    xa = xmin + delx * GS
    fa = f(xa, *args)
    xb = xmax - delx * GS
    fb = f(xb, *args)

    while delx >= ep:
        delxsav = delx
        if fb < fa:
            xmax = xb
            delx = xmax - xmin
            if delx == delxsav:
                result = (xmin + xmax) / 2
                return result
            xb = xa
            fb = fa
            xa = xmin + delx * GS
            fa = f(xa, *args)
        else:
            xmin = xa
            delx = xmax - xmin
            if delx == delxsav:
                result = (xmin + xmax) / 2
                return result
            xa = xb
            fa = fb
            xb = xmax - delx * GS
            fb = f(xb, *args)

    result = (xmin + xmax) / 2
    return result
