# Not part of the paper -- an exploratory experiment:
# identify the ORDER DISTRIBUTION OF THE PRIMES by matrix pencil on a log grid.
#
# The prime staircase (density form of the explicit formula):
#     sum_n Lambda(n) phi_h(lam - log n)  =  e^{lam + h^2/2}  -  sum_rho e^{rho*lam + rho^2 h^2/2}  -  tiny
# where phi_h = Gaussian of width h (smoothing makes the infinite zero comb
# effectively FINITE: zero rho = 1/2 + i*gamma enters with weight e^{-gamma^2 h^2/2}).
#
# So the residual  y(lam) = e^{lam + h^2/2} - F(lam)  is an exponential sum
#     y(lam_j) = sum_m c_m z_m^j ,   z_m = e^{rho_m * dlam}
# in the sample index j -- exactly eq. (4.10) of the paper, with the zeta zeros
# rho_m as the (complex!) orders. Matrix pencil (Hua-Sarkar) recovers z_m, and
# Prop. 4.1's modulus formula reads off  Re(rho) = ln|z| / dlam :
# the Riemann Hypothesis says every recovered node modulus gives Re(rho) = 1/2.
#
# Input: a prime sieve. No zeta function anywhere.

import math
import numpy as np

# ---------- prime powers up to X ----------
X = 1_700_000
sieve = bytearray([1]) * (X + 1)
sieve[0] = sieve[1] = 0
for i in range(2, int(X ** 0.5) + 1):
    if sieve[i]:
        sieve[i * i :: i] = bytearray(len(sieve[i * i :: i]))
logn, wt = [], []
for p in range(2, X + 1):
    if sieve[p]:
        lp = math.log(p)
        n = p
        while n <= X:
            logn.append(math.log(n))
            wt.append(lp)
            n *= p
logn = np.array(logn)
wt = np.array(wt)
print(f"input: {int(np.sum(np.frombuffer(bytes(sieve), dtype=np.uint8)))} primes, "
      f"{len(logn)} prime powers up to {X}. No zeta anywhere.\n")

# ---------- smoothed prime signal on a uniform log grid ----------
h = 0.16          # Gaussian smoothing width (truncates the zero comb ~ e^{-g^2 h^2/2})
dlam = 0.076      # log-grid step: the "spiral pitch" (|gamma|*dlam < pi for kept zeros)
N = 68            # samples;  N >= 2M with M ~ 12 nodes
lam0 = 8.4
lgrid = lam0 + dlam * np.arange(N)
c0 = 1.0 / (h * math.sqrt(2 * math.pi))

y = np.empty(N)
for j, lam in enumerate(lgrid):
    t = (lam - logn) / h
    m = np.abs(t) < 8.0
    F = c0 * np.sum(wt[m] * np.exp(-0.5 * t[m] ** 2))
    y[j] = math.exp(lam + h * h / 2) - F      # subtract the "steady state" (pole of zeta at s=1)

# ---------- matrix pencil (Hua-Sarkar) ----------
def pencil(sig, rel_tol):
    Np = len(sig)
    L = Np // 2
    H = np.array([[sig[i + j] for j in range(L + 1)] for i in range(Np - L)])
    U, s, Vh = np.linalg.svd(H, full_matrices=False)
    M = int(np.sum(s > rel_tol * s[0]))
    V = Vh[:M, :].conj().T
    z = np.linalg.eigvals(np.linalg.pinv(V[:-1, :]) @ V[1:, :])
    return z, s

def energies(sig, z):
    # least squares on unit-norm node profiles: |coef| = L2 energy carried by each node
    V = np.vander(z, N=len(sig), increasing=True).T.astype(complex)
    norms = np.linalg.norm(V, axis=0)
    c, *_ = np.linalg.lstsq(V / norms, sig.astype(complex), rcond=None)
    return np.abs(c), c / norms

z, s = pencil(y, 3e-11)
print("pencil singular values (model-order indicator, cf. paper Table 2):")
print("  " + "  ".join(f"{v:.2e}" for v in s[:14] / s[0]))
print()

energy, _ = energies(y, z)
keep = energy > 1e-6 * energy.max()
rho = np.log(z[keep]) / dlam
# the critical strip 0 < Re s < 1 is PROVEN (Hadamard / de la Vallee Poussin 1896):
# discarding nodes far outside it is legitimate prior knowledge, not RH
rho = rho[(rho.real > -2) & (rho.real < 3)]
found = sorted((r for r in rho if r.imag > 1.0), key=lambda r: r.imag)

true_g = [14.134725, 21.022040, 25.010858, 30.424876, 32.935062, 37.586178]
print(f"{'Re(rho) recovered':>18} {'Im(rho) recovered':>18} {'true zeta zero':>15} {'Im error':>10}")
for r in found:
    g = min(true_g, key=lambda t: abs(t - r.imag))
    print(f"{r.real:>18.6f} {r.imag:>18.6f} {g:>15.6f} {abs(r.imag - g):>10.2e}")
print("\nRH in the paper's language: every node modulus must give Re(rho) = 1/2.")

# ---------- bonus: DON'T subtract the steady state ----------
y_raw = np.array([math.exp(lam + h * h / 2) for lam in lgrid]) - y   # = F itself
# F contains the main term e^{lam}: the pole of zeta at s = 1 should appear as
# the dominant node with order rho = 1 -- the "steady state" of the primes
# (the exact analogue of y_inf in Proposition 3.2 of the paper).
z2, s2 = pencil(y_raw, 1e-7)
energy2, _ = energies(y_raw, z2)
dominant = z2[np.argmax(energy2)]
rho_dom = np.log(dominant) / dlam
print(f"\nwithout subtraction, dominant recovered order: {rho_dom.real:.8f} "
      f"{rho_dom.imag:+.2e}i   (pole of zeta at s = 1: the primes' steady state)")
