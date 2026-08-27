# E5 (F1, F2 in the paper): matrix-pencil estimation of the COMPLEX DIMENSIONS of the Cantor string,
# from (a) geometric data (the lengths) and (b) spectral data (the sound).
#
# Cantor string: lengths 3^-n with multiplicity 2^{n-1}.  Exact complex dimensions:
#     omega_k = D + i k p,   D = log2/log3,  p = 2pi/log3     (lattice case: all on one line)
# Geometric measure in lam = log x:   sum_n 2^{n-1} delta(lam - n log3)
#   = (1/(2 log3)) * sum_k e^{omega_k lam}                (exact, by Poisson summation)
# Spectral measure: atoms at lam = log(m * 3^n), weights 2^{n-1}; smooth part
#   = e^lam (Weyl term, total length 1) + sum_k (zeta(omega_k)/(2 log3)) e^{omega_k lam} + ...
# So the SAME pencil pipeline applies; the spectral amplitudes carry the channel gain zeta(omega).

import numpy as np
import math
import mpmath as mp

mp.mp.dps = 15
log3 = math.log(3.0)
D_true = math.log(2.0) / log3
p_true = 2 * math.pi / log3

h, dlt, N, lam0 = 0.10, 0.07, 64, 5.0
lam = lam0 + dlt * np.arange(N)
c0 = 1.0 / (h * math.sqrt(2 * math.pi))
lam_max = lam[-1] + 8 * h

def comb(positions, weights):
    F = np.zeros(N)
    pos = np.asarray(positions); w = np.asarray(weights)
    for i in range(N):
        t = (lam[i] - pos) / h
        m = np.abs(t) < 8.0
        F[i] = c0 * np.sum(w[m] * np.exp(-0.5 * t[m] ** 2))
    return F

# ---------- (a) geometric signal ----------
ns = np.arange(1, 14)
geo = comb(ns * log3, 2.0 ** (ns - 1))

# ---------- (b) spectral signal ----------
pos_s, w_s = [], []
for n in ns:
    m_max = int(math.exp(lam_max) / 3.0 ** n) + 1
    if m_max < 1:
        continue
    m = np.arange(1, m_max + 1, dtype=np.float64)
    pos_s.append(np.log(m) + n * log3)
    w_s.append(np.full(len(m), 2.0 ** (n - 1)))
spec = comb(np.concatenate(pos_s), np.concatenate(w_s)) - np.exp(lam + h * h / 2)  # subtract Weyl term

# ---------- pencil ----------
def pencil(sig, tol):
    L = len(sig) // 2
    H = np.array([[sig[i + j] for j in range(L + 1)] for i in range(len(sig) - L)])
    U, s, Vh = np.linalg.svd(H, full_matrices=False)
    M = int(np.sum(s > tol * s[0]))
    V = Vh[:M, :].conj().T
    z = np.linalg.eigvals(np.linalg.pinv(V[:-1, :]) @ V[1:, :])
    Vd = np.vander(z, N=len(sig), increasing=True).T.astype(complex)
    nrm = np.linalg.norm(Vd, axis=0)
    ce, *_ = np.linalg.lstsq(Vd / nrm, sig.astype(complex), rcond=None)
    return z, np.abs(ce), ce / nrm

def nodes(sig, tol, emin):
    z, energy, c = pencil(sig, tol)
    keep = energy > emin * energy.max()
    om = np.log(z[keep]) / dlt
    cc = c[keep]
    sel = (om.real > 0.2) & (om.real < 1.2) & (om.imag > -0.5)
    om, cc = om[sel], cc[sel]
    o = np.argsort(om.imag)
    return om[o], cc[o]

om_g, c_g = nodes(geo, 1e-9, 1e-5)
om_s, c_s = nodes(spec, 1e-7, 1e-4)

print(f"Cantor string.  Exact: D = {D_true:.6f},  p = {p_true:.6f}\n")
print("(a) complex dimensions from GEOMETRY (the lengths):")
print(f"{'Re (est)':>10} {'Im (est)':>10} {'Im exact':>10} {'k':>3}")
for om in om_g:
    k = round(om.imag / p_true)
    print(f"{om.real:>10.6f} {om.imag:>10.6f} {k * p_true:>10.6f} {k:>3}")

print("\n(b) complex dimensions from the SPECTRUM (the sound):")
print(f"{'Re (est)':>10} {'Im (est)':>10} {'k':>3} {'gain=|c_s/c_g|':>15} {'|zeta(omega_k)|':>16}")
for om, cs in zip(om_s, c_s):
    k = round(om.imag / p_true)
    match = [cg for og, cg in zip(om_g, c_g) if round(og.imag / p_true) == k]
    gain = abs(cs / match[0]) if match else float('nan')
    zval = abs(complex(mp.zeta(mp.mpc(D_true, k * p_true))))
    print(f"{om.real:>10.6f} {om.imag:>10.6f} {k:>3} {gain:>15.4f} {zval:>16.4f}")

print("\nLattice theorem check: every recovered Re should equal D — the Cantor string's own 'RH', provable and measured.")
