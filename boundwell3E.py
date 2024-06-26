import numpy as np

# Define common variables as a dictionary for easy reference in functions
common_vars = {
    'c': None,
    'effm': None,
    'alpha': None,
    'v0': None,
    'func': None
}


def feven(xsi, common_vars):
    c = common_vars['c']
    effm = common_vars['effm']
    alpha = common_vars['alpha']
    eta = np.sqrt(abs(c**2 - xsi**2 / effm))
    return (xsi * np.tan(xsi) / (effm * np.tanh(eta * alpha)) - eta)**2


def fodd(xsi, common_vars):
    c = common_vars['c']
    effm = common_vars['effm']
    alpha = common_vars['alpha']
    eta = np.sqrt(abs(c**2 - xsi**2 / effm))
    return (xsi * np.tanh(eta * alpha) / (np.tan(xsi) * effm) + eta)**2


def feven2(xsi, common_vars):
    c = common_vars['c']
    effm = common_vars['effm']
    alpha = common_vars['alpha']
    v0 = common_vars['v0']
    if v0 >= 0:
        eta = np.sqrt(abs((xsi**2 / effm) - c**2))
    else:
        eta = np.sqrt((xsi**2 / effm) + c**2)
    common_vars['func'] = (xsi * np.tan(xsi) / effm +
                           eta * np.tan(eta * alpha))**2
    return common_vars['func']


def fodd2(xsi, common_vars):
    c = common_vars['c']
    effm = common_vars['effm']
    alpha = common_vars['alpha']
    v0 = common_vars['v0']
    if v0 >= 0:
        eta = np.sqrt(abs((xsi**2 / effm) - c**2))
    else:
        eta = np.sqrt((xsi**2 / effm) + c**2)
    common_vars['func'] = (np.tan(xsi) * effm / xsi +
                           np.tan(eta * alpha) / eta)**2
    return common_vars['func']


def golden_section_search(func, a, b, tol, common_vars):
    invphi = (np.sqrt(5) - 1) / 2
    invphi2 = (3 - np.sqrt(5)) / 2
    h = b - a
    if h <= tol:
        return (a + b) / 2
    n = int(np.ceil(np.log(tol / h) / np.log(invphi)))
    c = a + invphi2 * h
    d = a + invphi * h
    yc = func(c, common_vars)
    yd = func(d, common_vars)
    for _ in range(n - 1):
        if yc < yd:
            b = d
            d = c
            c = a + invphi2 * (b - a)
            yd = yc
            yc = func(c, common_vars)
        else:
            a = c
            c = d
            d = a + invphi * (b - a)
            yc = yd
            yd = func(d, common_vars)
    if yc < yd:
        return (a + d) / 2
    else:
        return (c + b) / 2


def boundener3e(effm, el1, el2, v0, nstate):
    pi = 4.0 * np.arctan(1.0)
    a = el1 / 2.0
    el = el2 / 2.0
    alpha = (el / a) - 1.0
    hbar22m = 197.3**2 / (2.0 * 0.511e6)
    c = np.sqrt(abs(v0) * a**2 / hbar22m)
    n = 0

    common_vars['c'] = c
    common_vars['effm'] = effm
    common_vars['alpha'] = alpha
    common_vars['v0'] = v0

    if v0 < 0:
        return [], [], []

    n = 1 + int(2.0 * c * np.sqrt(effm) / pi)
    specfinal = True

    if n % 2 == 0:
        delxsi0 = np.arctan(effm / (alpha * (n - 1) * pi / 2.0))
        for _ in range(100):
            delxsi1 = np.arctan(
                effm / ((alpha * (n - 1) * pi / 2.0) + delxsi0))
            if abs(delxsi1 - delxsi0) < 1e-6:
                break
            delxsi0 = delxsi1
        else:
            print('*** error - BoundWell, delxsi, too many iterations')
        delxsi = delxsi1
        if (c * np.sqrt(effm) - (n - 1) * pi / 2.0) <= delxsi:
            specfinal = False
            n -= 1

    ener = np.zeros(nstate)
    amp1 = np.zeros(nstate)
    amp2 = np.zeros(nstate)

    for i in range(1, n + 1):
        xsi = i * pi / 2.0
        if i != n or not specfinal:
            for _ in range(1000):
                xsiold = xsi
                eta = np.sqrt(c**2 - xsiold**2 / effm)
                if i % 2 == 1:
                    xsi = np.arctan(
                        effm * eta * np.tanh(eta * alpha) / xsiold) + (i - 1) * pi / 2.0
                else:
                    xsi = np.arctan(
                        effm * eta / (xsiold * np.tanh(eta * alpha))) + (i - 1) * pi / 2.0
                if abs(xsi - xsiold) < 1e-6:
                    break
            else:
                print('*** error - reached max iterations')
        else:
            xmin = (i - 1) * pi / 2.0
            if n % 2 == 0:
                xmin += delxsi
            xmax = c * np.sqrt(effm)
            ep = 1e-6
            if i % 2 == 1:
                xsi = golden_section_search(feven, xmin, xmax, ep, common_vars)
            else:
                xsi = golden_section_search(fodd, xmin, xmax, ep, common_vars)

        ener[i-1] = hbar22m * (xsi / a)**2 / effm
        xsi = a * np.sqrt(effm * ener[i-1] / hbar22m)
        eta = a * np.sqrt(abs(v0 - ener[i-1]) / hbar22m)
        err = False
        if i % 2 == 1:
            if abs((xsi * np.tan(xsi) / effm) / (eta * np.tanh(eta * alpha)) - 1) > 1e-4:
                err = True
        else:
            if abs(xsi * np.tanh(eta * alpha) / (np.tan(xsi) * effm * eta) - 1) > 1e-4:
                err = True
        if err:
            print(f'*** error in solution, i={i}')

        wk = np.sqrt(ener[i-1] * effm / hbar22m)
        wkappa = np.sqrt((v0 - ener[i-1]) / hbar22m)
        if i % 2 == 1:
            atmp = (a + np.sin(2.0 * wk * a) / (2.0 * wk)) + \
                   (np.cos(wk * a) / np.cosh(wkappa * (el - a)))**2 * \
                   (el - a + np.sinh(2.0 * wkappa * (el - a)) / (2.0 * wkappa))
            amp2[i-1] = 1.0 / np.sqrt(atmp)
            amp1[i-1] = amp2[i-1] * np.cos(wk * a) / np.cosh(wkappa * (el - a))
        else:
            atmp = (a - np.sin(2.0 * wk * a) / (2.0 * wk)) + \
                   (np.sin(wk * a) / np.sinh(wkappa * (el - a)))**2 * \
                   (-el + a + np.sinh(2.0 * wkappa * (el - a)) / (2.0 * wkappa))
            amp2[i-1] = 1.0 / np.sqrt(atmp)
            amp1[i-1] = amp2[i-1] * np.sin(wk * a) / np.sinh(wkappa * (el - a))

    xsi = c * np.sqrt(effm) if v0 >= 0 else 0.0
    i = n

    while i < nstate:
        i += 1
        tmp1 = 0.0
        tmp1sav = 0.0
        for _ in range(1000):
            xsiold = xsi
            xsi += 0.01
            if v0 >= 0:
                eta = np.sqrt(abs((xsi**2 / effm) - c**2))
            else:
                eta = np.sqrt((xsi**2 / effm) + c**2)
            tmp1 = xsi * np.tan(xsi) / effm + eta * np.tan(eta * alpha)
            if tmp1 * tmp1sav < 0:
                break
            tmp1sav = tmp1

        ep = 1e-6
        xmin = xsiold
        xmax = xsi
        xsi = golden_section_search(feven2, xmin, xmax, ep, common_vars)

        if common_vars['func'] > 1e-6:
            continue

        ii = n + 2 * (i - n - 1) + 1 + n % 2
        ener[ii-1] = hbar22m * (xsi / a)**2 / effm
        xsi = a * np.sqrt(effm * ener[ii-1] / hbar22m)
        eta = a * np.sqrt(abs(ener[ii-1] - v0) / hbar22m)
        err = False
        if abs((xsi * np.tan(xsi) / effm) / (eta * np.tan(eta * alpha)) - 1) > 1e-4:
            err = True
        if err:
            print(f'*** error in even solution, i={i}')

        wk = np.sqrt(ener[ii-1] * effm / hbar22m)
        wkappa = np.sqrt((ener[ii-1] - v0) / hbar22m)
        atmp = (a + np.sin(2.0 * wk * a) / (2.0 * wk)) + \
               (np.cos(wk * a) / np.cos(wkappa * (el - a)))**2 * \
               (el - a + np.sin(2.0 * wkappa * (el - a)) / (2.0 * wkappa))
        amp2[ii-1] = 1.0 / np.sqrt(atmp)
        amp1[ii-1] = amp2[ii-1] * np.cos(wk * a) / np.cos(wkappa * (el - a))

    xsi = c * np.sqrt(effm) if v0 >= 0 else 0.0
    i = n

    while i < nstate:
        i += 1
        tmp1 = 0.0
        tmp1sav = 0.0
        for _ in range(1000):
            xsiold = xsi
            xsi += 0.01
            if v0 >= 0:
                eta = np.sqrt(abs((xsi**2 / effm) - c**2))
            else:
                eta = np.sqrt((xsi**2 / effm) + c**2)
            tmp1 = np.tan(xsi) * effm / xsi + np.tan(eta * alpha) / eta
            if tmp1 * tmp1sav < 0:
                break
            tmp1sav = tmp1

        ep = 1e-6
        xmin = xsiold
        xmax = xsi
        xsi = golden_section_search(fodd2, xmin, xmax, ep, common_vars)

        if common_vars['func'] > 1e-6:
            continue

        ii = n + 2 * (i - n - 1) + 1 + (n + 1) % 2
        ener[ii-1] = hbar22m * (xsi / a)**2 / effm
        xsi = a * np.sqrt(effm * ener[ii-1] / hbar22m)
        eta = a * np.sqrt(abs(ener[ii-1] - v0) / hbar22m)
        err = False
        if abs((np.tan(xsi) * effm * eta) / (xsi * np.tan(eta * alpha)) - 1) > 1e-4:
            err = True
        if err:
            print(f'*** error in odd solution, i={i}')

        wk = np.sqrt(ener[ii-1] * effm / hbar22m)
        wkappa = np.sqrt((ener[ii-1] - v0) / hbar22m)
        atmp = (a - np.sin(2.0 * wk * a) / (2.0 * wk)) + \
               (np.sin(wk * a) / np.sinh(wkappa * (el - a)))**2 * \
               (-el + a + np.sinh(2.0 * wkappa * (el - a)) / (2.0 * wkappa))
        amp2[ii-1] = 1.0 / np.sqrt(atmp)
        amp1[ii-1] = amp2[ii-1] * np.sin(wk * a) / np.sinh(wkappa * (el - a))

    return ener[:nstate], amp1[:nstate], amp2[:nstate]

# Example usage:
# effm = 1.0
# el1 = 10.0
# el2 = 20.0
# v0 = 1.0
# nstate = 5
# ener, amp1, amp2 = boundener3e(effm, el1, el2, v0, nstate)
# print(ener, amp1, amp2)
