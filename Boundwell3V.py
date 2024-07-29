import numpy as np

def bound_ener3v(effm, el1, el2, v0, n):
    pi = 4.0 * np.arctan(1.0)
    hbar22m = 197.3**2 / (2.0 * 0.511e6)
    
    a = el1 / 2.0
    el = el2 / 2.0
    alpha = (el / a) - 1.0
    c = np.sqrt(v0 * a**2 / hbar22m)
    
    ener = np.zeros(100)
    amp1 = np.zeros(100)
    amp2 = np.zeros(100)
    
    for i in range(1, n + 1):
        xsi = i * pi / 2.0
        for iter in range(1, 1001):
            xsiold = xsi
            eta = np.sqrt(c**2 + xsiold**2 / effm)
            if i % 2 == 1:
                xsi = np.arctan(effm * eta * np.tanh(eta * alpha) / xsiold) + (i - 1) * pi / 2
            else:
                xsi = np.arctan(effm * eta / (xsiold * np.tanh(eta * alpha))) + (i - 1) * pi / 2
            if np.abs(xsi - xsiold) < 1.0e-6:
                break
        else:
            print('*** error - reached max iterations')
        
        ener[i-1] = hbar22m * (xsi / a)**2 / effm
        
        xsi = a * np.sqrt(effm * ener[i-1] / hbar22m)
        eta = a * np.sqrt(abs(v0 + ener[i-1]) / hbar22m)
        err = False
        if i % 2 == 1:
            if np.abs((xsi * np.tan(xsi) * np.tanh(eta * alpha) / effm) / eta - 1) > 1.0e-4:
                err = True
        else:
            if np.abs(xsi * np.tanh(eta * alpha) / (np.tan(xsi) * effm * eta) - 1) > 1.0e-4:
                err = True
        # if err:
        #     print('*** error in solution, i=', i)
        
        wk = np.sqrt(ener[i-1] * effm / hbar22m)
        wkappa = np.sqrt((v0 + ener[i-1]) / hbar22m)
        
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
        
    return ener[:n], amp1[:n], amp2[:n]

# Example usage:
effm = 1.0
el1 = 10.0
el2 = 20.0
v0 = 5.0
n = 10

ener, amp1, amp2 = bound_ener3v(effm, el1, el2, v0, n)
print("ener:", ener)
print("amp1:", amp1)
print("amp2:", amp2)
