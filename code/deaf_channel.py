# The Lapidus-Maier "deaf channel" experiment.
#
# Generalized fractal string with geometric counting function
#     N_L(y) = y^D + eps * Re(y^{D+i*gamma}),   y >= 1     (oscillation at complex dimension D+i*gamma)
# Its frequency (spectral) counting function is the Dirichlet sum
#     N_nu(x) = sum_{k>=1} N_L(x/k)
# and analytically  N_nu = Weyl term + zeta(D)*x^D + eps*Re[ zeta(D+i*gamma) * x^{D+i*gamma} ] + small:
# the spectrum hears the geometric oscillation through the FILTER GAIN |zeta(D+i*gamma)|.
# If gamma is the ordinate of a zeta zero and D = 1/2: gain = 0 -> the spectrum is DEAF to it.
#
# Below: compute N_nu by direct summation (no zeta used anywhere), subtract the exact Weyl
# term, and measure the surviving oscillation amplitude at frequency gamma in log x.

import numpy as np
import math
import mpmath as mp

mp.mp.dps = 15
D = 0.5
eps = 0.5
lam = np.linspace(math.log(1e4), math.log(5e5), 96)   # log x grid
X = np.exp(lam)

def spectrum_residual(gamma):
    c = complex(1/(1-D), 0) + eps*(1/(1-D-1j*gamma))   # Weyl coefficient (both terms)
    res = np.empty(len(X))
    for i, x in enumerate(X):
        k = np.arange(1, int(x) + 1, dtype=np.float64)
        y = x / k
        Nnu = np.sum(y**D) + eps*np.sum((y**D)*np.cos(gamma*np.log(y)))
        res[i] = (Nnu - x*c.real) / x**D           # -> zeta(D) + eps*Re[zeta(D+ig) e^{ig*lam}]
    return res

def osc_amplitude(res, gamma):
    r = res - res.mean()
    w = np.hanning(len(lam))                        # window against leakage
    z = np.trapezoid(r*w*np.exp(-1j*gamma*lam), lam) / np.trapezoid(w, lam) * 2
    return abs(z)

g_zero = 14.134725142                                # first zeta zero ordinate
g_ctrl = 10.0                                        # control frequency (not a zero)

for name, g in [("gamma = 10.0      (NOT a zero)", g_ctrl),
                ("gamma = 14.134725 (zeta zero) ", g_zero)]:
    res = spectrum_residual(g)
    A = osc_amplitude(res, g)
    gain = abs(complex(mp.zeta(mp.mpc(D, g))))
    print(f"{name}:  geometric osc. amplitude = {eps:.2f}")
    print(f"   measured amplitude in SPECTRUM = {A:.4f}   predicted eps*|zeta| = {eps*gain:.4f}")
    print(f"   channel gain |zeta({D}+i*gamma)| = {gain:.4f}\n")

zD = complex(mp.zeta(D))
print(f"sanity: mean of residual = {spectrum_residual(g_ctrl).mean():.4f}  vs  zeta(1/2) = {zD.real:.4f}")
