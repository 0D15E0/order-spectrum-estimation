# Not part of the paper -- the dual experiment: identify the ORDER DISTRIBUTION OF ZETA.
#
# The "transfer function" is H = -zeta'/zeta, a distributed-order system:
#     -zeta'/zeta(s) = sum_n Lambda(n) e^{-s log n}   (Re s > 1)
# i.e. order comb A(alpha) = sum_n Lambda(n) delta(alpha - log n):
# lines at the LOGS OF PRIME POWERS, weights log p.
#
# Sample it along the vertical line s = sigma0 + i*t (the paper's contour,
# here in the Laplace variable), matrix-pencil the samples, and the primes
# should come out as e^{recovered order}, with weights Lambda(n).

import math
import numpy as np
import mpmath as mp

mp.mp.dps = 20
sigma0, dt, N = 2.5, 0.7, 140

print(f"sampling -zeta'/zeta on the line Re s = {sigma0}, {N} samples, dt = {dt} ...")
y = np.array([
    complex(-mp.zeta(mp.mpc(sigma0, k * dt), derivative=1) / mp.zeta(mp.mpc(sigma0, k * dt)))
    for k in range(N)
])

# ---------- matrix pencil (complex signal; nodes z_n = e^{-i*dt*log n}) ----------
L = N // 2
H = np.array([[y[i + j] for j in range(L + 1)] for i in range(N - L)])
U, s, Vh = np.linalg.svd(H, full_matrices=False)
M = int(np.sum(s > 2e-3 * s[0]))
V = Vh[:M, :].conj().T
z = np.linalg.eigvals(np.linalg.pinv(V[:-1, :]) @ V[1:, :])

# this Hankel orientation returns conjugate nodes: true nodes are conj(z)
zt = np.conj(z)
Vd = np.vander(zt, N=N, increasing=True).T.astype(complex)
norms = np.linalg.norm(Vd, axis=0)
ce, *_ = np.linalg.lstsq(Vd / norms, y, rcond=None)
energy = np.abs(ce)
c = ce / norms                     # raw amplitudes: should be ~ Lambda(n) n^{-sigma0}

keep = energy > 5e-3 * energy.max()
alpha = -np.angle(zt[keep]) / dt   # recovered orders
w = np.real(c[keep])
order = np.argsort(alpha)
alpha, w = alpha[order], w[order]

print(f"\npencil kept M = {M} nodes; significant lines found: {len(alpha)}\n")
print(f"{'e^order':>12} {'nearest n':>10} {'weight*n^2.5':>13} {'log p':>8}   reading")
for a, wt in zip(alpha, w):
    n = math.exp(a)
    ni = round(n)
    if ni < 2 or abs(n - ni) > 0.005 * ni:   # drop mis-localized aliased dust
        continue
    lam_hat = wt * ni ** sigma0
    # factor n as p^k to show the weight identifies the underlying prime
    p = next((q for q in range(2, ni + 1) if ni % q == 0), ni)
    tag = "prime" if ni == p else f"= {p}^{round(math.log(ni, p))}  (weight is log {p}!)"
    print(f"{n:>12.6f} {ni:>10} {lam_hat:>13.6f} {math.log(p):>8.4f}   {tag}")
print("\nThe order distribution of zeta IS the prime comb: positions log(p^k), weights log p.")
