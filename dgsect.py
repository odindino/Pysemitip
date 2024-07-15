import numpy as np

def golden_section_search(f, xmin, xmax, ep):
    """
    Golden Section Search to find the minimum of the function f
    in the interval [xmin, xmax] with precision ep.
    
    Parameters:
    f    -- function to minimize
    xmin -- lower bound of the interval
    xmax -- upper bound of the interval
    ep   -- precision of the search
    
    Returns:
    A tuple (xmin, f(xmin)) where xmin is the point at which the minimum is achieved.
    """
    GS = 0.3819660
    
    if xmax == xmin:
        return xmin, f(xmin)
    if ep == 0.0:
        return xmin, f(xmin)
    if xmax < xmin:
        xmin, xmax = xmax, xmin
    
    delx = xmax - xmin
    xa = xmin + delx * GS
    fa = f(xa)
    xb = xmax - delx * GS
    fb = f(xb)
    
    while delx >= ep:
        delxsav = delx
        if fb < fa:
            xmin = xa
            xa = xb
            fa = fb
            delx = xmax - xmin
            if delx == delxsav:
                break
            xb = xmax - delx * GS
            fb = f(xb)
        else:
            xmax = xb
            xb = xa
            fb = fa
            delx = xmax - xmin
            if delx == delxsav:
                break
            xa = xmin + delx * GS
            fa = f(xa)
    
    # Print the final values of xmin and xmax
    print(f"Final xmin: {xmin}, xmax: {xmax}")
    
    # Return the point and value of the minimum found
    if fa < fb:
        return xa, fa
    else:
        return xb, fb


