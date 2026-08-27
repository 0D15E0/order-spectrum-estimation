# Appendix — The Order Distribution of the Prime Numbers

> **Integration notes.** Companion draft for *"Identification of the order distribution of linear systems by line-spectrum estimation on logarithmic spirals"*. Written to follow the notation of the main paper: representation (1.1), estimator (1.2), Theorem 3.1, Proposition 3.2, Corollary 3.1, spiral (4.9), exponential-sum model (4.10), Proposition 4.1, matrix pencil [17]. Three new references are appended as [28]–[30]. **Status: not part of the paper.** The paper's E5 calibration benchmark uses the Cantor fractal-string experiments described in `fractal-strings-pencil-note.md` instead; this document explores a distinct, unused benchmark (the primes/zeta function, labeled P1/P2 below) and is kept for the record but is not referenced by the paper. Reproduction scripts: `code/prime_pencil.py`, `code/zeta_comb.py`.

This appendix applies the estimator of Section 4 to the oldest non-commensurate order comb in mathematics. It is a demonstration on a system of independent interest rather than a new identification result: the object recovered is classical, but the estimator turns out to measure precisely the quantity that the Riemann Hypothesis is about, and the two structural obstructions of Section 3 (steady-state bias, readout aliasing) both reappear here in recognizable form. All statements used are unconditional theorems of analytic number theory [28, 29].

## C.1 The zeta function as a distributed-order system

Let $\Lambda$ denote the von Mangoldt function: $\Lambda(n) = \log p$ if $n = p^k$ is a prime power, $\Lambda(n) = 0$ otherwise. For $\operatorname{Re} p > 1$ the logarithmic derivative of the Riemann zeta function is the absolutely convergent Dirichlet series $-\zeta'/\zeta(p) = \sum_{n \ge 2} \Lambda(n)\, n^{-p}$. Under the substitution $s = e^{-p}$ of Section 2 this is the representation (1.1) with

$$
H(s) = \sum_{n \ge 2} \Lambda(n)\, s^{\log n},
\qquad
\mathcal{A}(\alpha) = \sum_{n \ge 2} \Lambda(n)\, \delta(\alpha - \log n).
\tag{C.1}
$$

Thus $-\zeta'/\zeta$ is a distributed-order system whose order distribution is the **prime comb**: spectral lines at the logarithms of the prime powers, the line at $\log p^k$ carrying the weight $\log p$ of its underlying prime. The support is restricted to prime powers — there is no line at $\log 6$ — because the Euler product $\zeta(p) = \prod_p (1 - p^{-p})^{-1}$ collapses the "all integers" comb of $\zeta$ itself to prime powers upon taking the logarithmic derivative; the comb (C.1) is unique factorization written in the order domain. Since $\{\log p\}_{p\ \mathrm{prime}}$ is linearly independent over $\mathbb{Q}$, the comb is genuinely non-commensurate: in the taxonomy of Section 1, the primes are a strictly fractional system.

## C.2 The dual signal and its complex orders

The response side of this system is the Chebyshev step function $\psi(x) = \sum_{p^k \le x} \log p$. The explicit formula of Riemann–von Mangoldt ([28], Ch. 17) gives, for the measure $d\psi$ read in the logarithmic variable $\lambda = \log x$, the density

$$
e^{\lambda} \;-\; \sum_{\rho} e^{\rho\lambda} \;-\; \bigl(e^{2\lambda} - 1\bigr)^{-1},
\tag{C.2}
$$

the sum running over the nontrivial zeros $\rho = \beta + i\gamma$ of $\zeta$ (in the symmetric limit). Equation (C.2) says that on a logarithmic grid the prime staircase is an exponential sum whose **orders are complex**: one real order at $1$ (the pole of $\zeta$) and one order at every zero. Two classical facts frame the experiment: $0 < \beta < 1$ for every zero (Hadamard, de la Vallée Poussin, 1896), and the Riemann Hypothesis (RH) is the statement $\beta = \tfrac12$ for all $\rho$. In the language of this paper: **RH asserts that the order comb of the primes is uniformly damped — all nodes on one vertical line of the order plane.** By Proposition 4.1 the real part of an order is read from a node modulus, so the estimator of Section 4 measures exactly the RH-relevant quantity.

The comb (C.2) is infinite, which Section 7 identifies as the ill-posed regime; a smoothing window restores the finite-comb setting:

**Proposition C.1 (Gaussian finite-comb reduction).** *Let $\varphi_h(t) = (h\sqrt{2\pi})^{-1} e^{-t^2/2h^2}$ and $0 < h \le \tfrac14$. Then for $\lambda \ge 2$*

$$
F_h(\lambda) := \sum_{n \ge 2} \Lambda(n)\, \varphi_h(\lambda - \log n)
= e^{h^2/2}\, e^{\lambda}
- \sum_{\rho} e^{\rho^2 h^2 / 2}\, e^{\rho\lambda}
+ O\!\left(e^{-2\lambda}\right),
\tag{C.3}
$$

*the sum over zeros converging absolutely, since $\bigl| e^{\rho^2 h^2/2} \bigr| = e^{(\beta^2 - \gamma^2) h^2 / 2} \le e^{h^2/2}\, e^{-\gamma^2 h^2 / 2}$ while the zero counting function grows only logarithmically in density.*

*Proof (sketch).* Integrate $\varphi_h(\lambda - \cdot)$ against (C.2); the Gaussian–exponential integrals evaluate to $\int \varphi_h(u)\, e^{-\rho u}\, du = e^{\rho^2 h^2/2}$, and the trivial-zero term $(e^{2t} - 1)^{-1} = \sum_{k \ge 1} e^{-2kt}$ contributes $\sum_k e^{2k^2 h^2 - 2k\lambda} = O(e^{-2\lambda})$ under the stated restrictions on $h$ and $\lambda$. See [28], Ch. 17, for the underlying explicit formula. ∎

The Gaussian weight $e^{-\gamma^2 h^2/2}$ truncates the comb to an effective model order $M(h) \approx \#\{\gamma : \gamma \lesssim h^{-1}\sqrt{2\log(1/\varepsilon)}\}$ at working accuracy $\varepsilon$: the window plays for the zeros the role the finite-comb hypothesis plays in Section 4.

## C.3 The estimator

Fix $\lambda_j = \lambda_0 + j\Delta$, $j = 0, \dots, N-1$, and form from prime data alone the residual

$$
y_j = e^{h^2/2}\, e^{\lambda_j} - F_h(\lambda_j)
= \sum_{\rho} c_\rho\, z_\rho^{\,j} + \varepsilon_j,
\qquad
z_\rho = e^{\rho\Delta}, \quad
c_\rho = e^{\rho^2 h^2/2}\, e^{\rho\lambda_0},
\tag{C.4}
$$

an exponential sum of the form (4.10) with complex orders. Relative to Section 4 the roles of contour and orders are interchanged: there, a complex contour $s = e^{(\sigma + i)u}$ resolved real orders; here a real logarithmic grid (the "real-node Prony" variant noted after Proposition 4.1) resolves complex orders. The matrix pencil [17] applied to (C.4) returns the nodes; by the modulus formula of Proposition 4.1, $\hat\beta = \ln|z_\rho| / \Delta$ with no branch condition, while $\hat\gamma = \arg z_\rho / \Delta$ requires $|\gamma\Delta| < \pi$ for the retained zeros — the analogue of the principal-strip condition there. Since the signal is real, nodes arrive in conjugate pairs.

The subtraction of $e^{h^2/2} e^{\lambda}$ in (C.4) is Corollary 3.1 verbatim: the pole of $\zeta$ at $1$ is the **steady state of the prime system**, and, exactly as in Proposition 3.2, leaving it in place buries the spectrum (see Table C.1, last row, where it is instead *recovered*).

## C.4 Experiment

Primes and prime powers up to $X = 1.7 \times 10^6$ (128,141 primes, 128,430 prime powers) were generated by sieve; parameters $h = 0.16$, $\Delta = 0.076$, $\lambda_0 = 8.4$, $N = 68$, pencil parameter $L = N/2$, quadrature-free evaluation of (C.4) in double precision. **No value of $\zeta$ is computed anywhere.** The singular-value profile of the pencil matrix exhibits the conjugate-pair structure ($1.00,\ 0.97 \mid 4.4,\ 4.4 \times 10^{-2} \mid 4.0,\ 3.7 \times 10^{-3} \mid 8.5,\ 8.4 \times 10^{-5} \mid \cdots$), one pair per zero, mirroring the model-order indicator of Table 2. Nodes were retained by an energy threshold on unit-norm profiles, plus the *a priori* strip condition $-2 < \hat\beta < 3$ justified by the 1896 zero-free bounds (one spurious endpoint node with $\hat\beta \approx 39$, an artifact of the sieve cutoff — see C.6 — is removed by it).

**Table C.1.** P1: order comb of the primes recovered from a sieve to $1.7 \times 10^6$. The first column is the measured quantity: RH predicts $0.5$ in every row.

| $\hat\beta = \ln\lvert z\rvert/\Delta$ | $\hat\gamma = \arg z/\Delta$ | $\gamma$ (true) | $\lvert\hat\gamma - \gamma\rvert$ |
|---|---|---|---|
| 0.500000 | 14.134725 | 14.134725 | 1.4 × 10⁻⁷ |
| 0.500000 | 21.022040 | 21.022040 | 3.6 × 10⁻⁷ |
| 0.500000 | 25.010858 | 25.010858 | 4.3 × 10⁻⁷ |
| 0.500000 | 30.424876 | 30.424876 | 4.2 × 10⁻⁷ |
| 0.500003 | 32.935060 | 32.935062 | 1.8 × 10⁻⁶ |

*Without steady-state subtraction: dominant node $\hat\beta = 0.99999998$, $\hat\gamma = 0$ — the pole of $\zeta$.*

The first column is the point of the experiment: the estimator does not assume where the zeros are; the node moduli **measure** their real parts, and return $\tfrac12$ to six decimals.

## C.5 The reverse direction

Sampling the transform side instead — $-\zeta'/\zeta(\sigma_0 + it)$ on the line $\sigma_0 = \tfrac52$, $N = 140$, $\Delta t = 0.7$, the "$H$ known analytically" regime of Section 3 — and applying the same pencil recovers the comb (C.1) itself:

| $e^{\hat\alpha}$ | nearest $n$ | $\hat\Lambda$ (recovered) | $\Lambda(n)$ (true) | reading |
|---|---|---|---|---|
| 1.999999 | 2 | 0.6933 | log 2 = 0.6931 | prime |
| 3.000009 | 3 | 1.0992 | log 3 = 1.0986 | prime |
| 4.000147 | 4 | 0.6936 | log 2 = 0.6931 | $= 2^2$ — weight is log 2, not log 4 |
| 5.000024 | 5 | 1.6130 | log 5 = 1.6094 | prime |
| 7.000026 | 7 | 1.9464 | log 7 = 1.9459 | prime |
| 8.000241 | 8 | 0.6934 | log 2 = 0.6931 | $= 2^3$ — weight is log 2 |
| 9.000475 | 9 | 1.0982 | log 3 = 1.0986 | $= 3^2$ — weight is log 3 |
| 11.0022 | 11 | 2.4211 | log 11 = 2.3979 | prime |
| 13.0013 | 13 | 2.5486 | log 13 = 2.5649 | prime |

Note the absent line at 6 (and 10, 12, 14, …) and the weights of the lines at 4, 8, 9, which are log 2, log 2, log 3 — not log 4, log 8, log 9: the estimator resolves the prime-power structure of the comb, i.e., **the Euler product is directly visible in the recovered spectrum**. Beyond $n \approx 25$ the weights degrade (the infinite tail of the comb acts as structured noise), in line with the limitation of Section 7. The two runs are dual in the sense of the explicit formula: fed the transform, the estimator returns the primes; fed the primes, it returns the zeros.

## C.6 Scope

P1 is a measurement, and three walls separate it — and any estimator-based approach — from a statement about RH itself.

1. **Finite data constrain finitely many zeros.** The window that makes the comb finite (Proposition C.1) is exactly the operation that erases all zeros beyond $\gamma \sim h^{-1}$, and removing it re-enters the exponentially ill-posed regime of Section 7. In the same vein the sieve bound must satisfy $\log X \gtrsim \lambda_{\max} + 8h$; the run above sits just below this ($\log X = \lambda_{\max} + 5.3h$), the induced endpoint bias appearing as the single out-of-strip node noted in C.4.
2. **A measurement cannot certify exact equality.** The pencil returns $\hat\beta = 0.500000 \pm 10^{-6}$, which is compatible with, but can never prove, $\beta = \tfrac12$; rigorous verification of zeros proceeds by a different, sign-change mechanism and currently covers $|\gamma| \le 3 \times 10^{12}$ [30].
3. **Resolution is logarithmic in the data.** $\hat\beta$ is a damping exponent read over a window of length $\log X$, so each additional digit costs an exponential factor in $X$.

The value of P1 is accordingly not evidentiary but structural: it exhibits the primes as the canonical infinite, exactly structured order comb — the natural boundary object of the method of this paper — and restates the Riemann Hypothesis in the vocabulary of order distributions: **the order comb of the primes is uniformly damped**.

Both experiments are reproduced by `prime_pencil.py` (P1; NumPy only) and `zeta_comb.py` (P2 direction; NumPy and mpmath), included with the code of Appendix B.

## Additional references (renumber on merge)

- [28] H. Davenport, *Multiplicative Number Theory*, 3rd Ed., Graduate Texts in Mathematics **74**, Springer, New York, 2000.
- [29] H. M. Edwards, *Riemann's Zeta Function*, Academic Press, New York, 1974; Dover reprint, 2001.
- [30] D. J. Platt, T. S. Trudgian, The Riemann hypothesis is true up to $3 \cdot 10^{12}$, *Bull. London Math. Soc.* **53** (2021), 792–797.
