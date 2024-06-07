import numpy as np
from semirhomult import fd


def sig(iar, index, e):
    # Placeholder for the SIG function implementation
    # This function needs to be implemented based on your requirements
    return 1.0  # example value


def rhos1(iar, ef1, dele, surf):
    sum_val = 0.0
    e = surf["EN"][iar, 0]

    if ef1 == surf["EN"][iar, 0]:
        return sum_val

    if surf["ISTK"] == 0:
        if ef1 < surf["EN"][iar, 0]:
            while e < ef1:
                e += dele
                sum_val += sig(iar, 0, e) * dele
        else:
            while e > ef1:
                e -= dele
                sum_val += sig(iar, 0, e) * dele
        return sum_val

    if ef1 < surf["EN"][iar, 0]:
        while e < ef1 + 10.0 * surf["TK"]:
            e += dele
            sum_val += sig(iar, 0, e) * fd(e, ef1, surf["TK"]) * dele
    else:
        while e > ef1 - 10.0 * surf["TK"]:
            e -= dele
            sum_val += sig(iar, 0, e) * (1.0 - fd(e, ef1, surf["TK"])) * dele

    return sum_val


def rhos2(iar, ef1, dele, surf):
    sum_val = 0.0
    e = surf["EN"][iar, 1]

    if ef1 == surf["EN"][iar, 1]:
        return sum_val

    if surf["ISTK"] == 0:
        if ef1 < surf["EN"][iar, 1]:
            while e < ef1:
                e += dele
                sum_val += sig(iar, 1, e) * dele
        else:
            while e > ef1:
                e -= dele
                sum_val += sig(iar, 1, e) * dele
        return sum_val

    if ef1 < surf["EN"][iar, 1]:
        while e < ef1 + 10.0 * surf["TK"]:
            e += dele
            sum_val += sig(iar, 1, e) * fd(e, ef1, surf["TK"]) * dele
    else:
        while e > ef1 - 10.0 * surf["TK"]:
            e -= dele
            sum_val += sig(iar, 1, e) * (1.0 - fd(e, ef1, surf["TK"])) * dele

    return sum_val


def rhos(iar, ef1, dele, surf):
    return rhos1(iar, ef1, dele, surf) + rhos2(iar, ef1, dele, surf)


def sigsum(iar, ener):
    return sig(iar, 0, ener) + sig(iar, 1, ener)


def enfind(iar, en1, en2, ne, surf):
    en0 = en1
    estart = en1
    dele = abs(en1 - en2) / float(ne)
    if dele == 0:
        return en0
    if en2 < en1:
        dele = -dele

    for ie in range(ne + 1):
        ef1 = estart + ie * dele
        sigtmp = rhos(iar, ef1, abs(dele), surf)
        if dele > 0 and sigtmp <= 0:
            en0 = ef1
            break
        elif dele < 0 and sigtmp >= 0:
            en0 = ef1
            break

    return en0
