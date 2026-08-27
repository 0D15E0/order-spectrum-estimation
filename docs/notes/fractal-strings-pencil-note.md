# Estimating the Complex Dimensions of Fractal Strings by Matrix Pencil

**Ulises Merlan**¹

> **Epistemic status.** Methods note with completed numerical experiments (F1–F3 below, reproduced by `code/cantor_pencil.py` and `code/deaf_channel.py`). The estimation results are validated against exactly known ground truth. Section 6 is a research program (E7), not completed work. Nothing here claims progress on the Riemann Hypothesis; Section 7 states the identifiability boundary explicitly. Companion to the author's paper on order-distribution identification, its experiment E5, and the estimator-criterion note.

## Abstract

A fractal string — a bounded open subset of $\mathbb{R}$ with lengths $\ell_1, \ell_2, \dots$ — is characterized by its geometric zeta function $\zeta_{\mathcal{L}}(s) = \sum_j \ell_j^s$, whose poles, the *complex dimensions*, govern the string's geometric and spectral oscillations (Lapidus–van Frankenhuijsen). We observe that $\zeta_{\mathcal{L}}$ is a distributed-order transfer function in the sense of the companion paper — the order distribution being the comb of log-lengths — and that the matrix-pencil estimator developed there applies verbatim: complex dimensions can be *estimated by super-resolution line-spectrum methods* from geometric or spectral counting data on a logarithmic grid. For the Cantor string, whose complex dimensions are known exactly ($\omega_k = D + ikp$ with $D = \log 2/\log 3$, $p = 2\pi/\log 3$), the estimator recovers eight complex dimensions from geometric data with all real parts equal to $0.630930$ (exact to six decimals) — an instrumental verification of the lattice theorem — and recovers the same dimensions from spectral data, with the ratio of spectral to geometric amplitudes matching $|\zeta(\omega_k)|$ to four decimals in every mode: an *acoustic measurement of the Riemann zeta function along a vertical line*, realizing the spectral factorization $\zeta_\nu = \zeta_{\mathcal{L}} \cdot \zeta$ as a measured transfer function. A third experiment exhibits the Lapidus–Maier "deaf channel": a geometric oscillation at complex dimension $\tfrac12 + i\gamma_1$ (the first zeta zero) enters the spectrum with measured amplitude $0.0007$ against a geometric amplitude of $0.5$, while a control frequency passes with gain equal to $|\zeta|$ to four decimals. We close with a two-part research program: nonlattice strings as measurable instances of non-uniform damping, and Beurling generalized primes with designed zeros as control-group zeta functions for the detection branch of the estimator criterion.

## 1. Introduction

The companion paper identifies the order distribution $\mathcal{A}$ of a system $H(s) = \int \mathcal{A}(\alpha) s^\alpha d\alpha$ by sampling on contours where each order becomes a node of an exponential sum, and estimating the nodes by the matrix pencil of Hua–Sarkar. Its experiment E5 applied the pipeline to the prime-counting signal, whose "order comb" is carried by the nontrivial zeros of $\zeta$.

Fractal geometry supplies a second, and in one respect better, application domain. In the theory of Lapidus and van Frankenhuijsen [3], a fractal string $\mathcal{L} = \{\ell_j\}$ has geometric zeta function $\zeta_{\mathcal{L}}(s) = \sum_j \ell_j^s$ — formally identical to the companion paper's representation (1.1) with order comb at the log-lengths — and the poles of $\zeta_{\mathcal{L}}$, the **complex dimensions**, play exactly the role the zeta zeros play for the primes: they are the frequencies of geometric oscillation, entering tube volumes and counting functions through explicit formulas of Riemann–von Mangoldt type. The paper's *lattice/nonlattice dichotomy* (commensurable versus incommensurable scaling ratios) is precisely the companion paper's commensurate/strictly-fractional taxonomy of order combs.

The advantage over the arithmetic case: fractal geometry offers infinitely many strings whose complex dimensions are **known exactly**, so the estimator can be validated end-to-end — and, unlike $\zeta$, the theory contains systems where the "all dimensions on one line" property *provably holds* (lattice strings) and systems where it *provably fails* (nonlattice strings), giving both phenotypes of the damping question a laboratory realization.

## 2. Background

**Fractal strings and complex dimensions** [1, 2, 3]. For a string with lengths $\ell_j$, the geometric counting function is $N_{\mathcal{L}}(x) = \#\{j : \ell_j^{-1} \le x\}$; the frequency (spectral) counting function of the associated Dirichlet Laplacian is $N_\nu(x) = \sum_j \lfloor \ell_j x \rfloor = \sum_{m \ge 1} N_{\mathcal{L}}(x/m)$. The Dirichlet-convolution structure of the second formula gives the fundamental **spectral factorization**

$$
\zeta_\nu(s) = \zeta_{\mathcal{L}}(s) \cdot \zeta(s):
$$

the sound of a fractal string is its geometry filtered through the Riemann zeta function. Explicit formulas expand both $N_{\mathcal{L}}$ and $N_\nu$ over the complex dimensions $\omega$ (poles of $\zeta_{\mathcal{L}}$); each geometric oscillation $x^\omega$ reaches the spectrum multiplied by the **channel gain** $\zeta(\omega)$.

**The inverse spectral problem** [2]. Lapidus–Maier: one can "hear" Minkowski measurability of strings of dimension $D$ if and only if $\zeta$ has no zeros on the line $\operatorname{Re} s = D$; hence RH is equivalent to solvability of the inverse problem for all $D \in (0,1)$, $D \neq \tfrac12$. In estimation language: **zeros of $\zeta$ are the null frequencies of the geometry-to-spectrum channel**, and the critical line is its only possibly-deaf band.

**The test case.** The Cantor string (complement of the middle-thirds Cantor set in $[0,1]$) has lengths $3^{-n}$ with multiplicities $2^{n-1}$, geometric zeta function $\zeta_{\mathcal{L}}(s) = 3^{-s}/(1 - 2 \cdot 3^{-s})$, and complex dimensions known exactly:

$$
\omega_k = D + ikp, \qquad D = \frac{\log 2}{\log 3} = 0.6309\,2975\ldots, \quad p = \frac{2\pi}{\log 3} = 5.7191\,6993\ldots, \quad k \in \mathbb{Z}.
$$

All on one vertical line — the **lattice case**. For self-similar strings this is a theorem, so the "uniform damping" that is conjectural for $\zeta$ (RH) is certified ground truth here.

## 3. Method

In the logarithmic variable $\lambda = \log x$, the geometric measure of the Cantor string is a weighted comb, and Poisson summation gives its mode expansion **exactly** — with no error term (contrast E5's $O(e^{-2\lambda})$):

$$
\sum_{n \ge 1} 2^{n-1} \delta(\lambda - n \log 3) \;=\; \frac{1}{2\log 3} \sum_{k \in \mathbb{Z}} e^{\omega_k \lambda}.
$$

Gaussian smoothing by $\varphi_h$ multiplies mode $k$ by $e^{\omega_k^2 h^2/2}$ — truncating the tower to an effective model order, exactly as in E5 — and sampling at $\lambda_j = \lambda_0 + j\Delta$ produces the exponential-sum model (4.10) of the companion paper with nodes $z_k = e^{\omega_k \Delta}$. The matrix pencil returns the nodes; $\operatorname{Re}\hat\omega = \ln|z|/\Delta$ (damping = Minkowski-dimension line), $\operatorname{Im}\hat\omega = \arg z/\Delta$ (oscillatory dimensions), valid for $|\operatorname{Im}\omega_k| \Delta < \pi$.

For the spectral run, the smoothed spectral measure has atoms at $\lambda = \log(m \cdot 3^n)$ with weights $2^{n-1}$; its smooth expansion is the Weyl term $e^{\lambda}$ (total length $1$ — the string's steady state, subtracted exactly as in the paper's Corollary 3.1 and E5's pole subtraction) plus the modes $\frac{\zeta(\omega_k)}{2\log 3} e^{\omega_k \lambda}$.

Parameters for both runs: $h = 0.10$, $\Delta = 0.07$, $N = 64$, $\lambda_0 = 5.0$, pencil parameter $N/2$, energy-threshold retention, a-priori filter $\operatorname{Re}\hat\omega \in (0.2, 1.2)$.

## 4. Experiment F1: complex dimensions from geometry

Input: the lengths only. Output:

| $\operatorname{Re}\hat\omega$ | $\operatorname{Im}\hat\omega$ | exact $\operatorname{Im} = kp$ | $k$ |
|---|---|---|---|
| 0.630930 | 0.000000 | 0.000000 | 0 |
| 0.630930 | 5.719202 | 5.719202 | 1 |
| 0.630930 | 11.438403 | 11.438403 | 2 |
| 0.630930 | 17.157605 | 17.157605 | 3 |
| 0.630930 | 22.876807 | 22.876807 | 4 |
| 0.630930 | 28.596009 | 28.596009 | 5 |
| 0.630930 | 34.315211 | 34.315210 | 6 |
| 0.630928 | 40.034412 | 40.034412 | 7 |

Every real part equals $D = \log 2/\log 3$ to the printed six decimals: the **lattice theorem measured by the instrument** — a uniformly damped order comb, this time as certified truth rather than (as in E5) 166-year-old conjecture. A ninth node ($k = 8$) emerges degraded at the smoothing floor $e^{-k^2 p^2 h^2/2}$, exactly where the weight ladder predicts the retained set to end.

## 5. Experiment F2: the same dimensions from the sound — with $\zeta$ as the measured gain

Input: the spectral counting data, Weyl term subtracted. The pencil recovers the same tower (all $\operatorname{Re}\hat\omega = 0.630930$ for $k \le 5$, drifting to $0.630937$ at the floor), and the amplitude ratio of each spectral mode to its geometric counterpart matches the channel gain independently computed from $\zeta$:

| $k$ | measured gain $\lvert c^{\mathrm{spec}}_k / c^{\mathrm{geo}}_k \rvert$ | $\lvert\zeta(D + ikp)\rvert$ |
|---|---|---|
| 0 | 2.1598 | 2.1598 |
| 1 | 0.8594 | 0.8594 |
| 2 | 1.3693 | 1.3693 |
| 3 | 2.0689 | 2.0689 |
| 4 | 1.3275 | 1.3275 |
| 5 | 2.2269 | 2.2269 |
| 6 | 1.9284 | 1.9284 |
| 7 | 1.1398 | 1.1398 |

Four-decimal agreement in every mode: **an acoustic measurement of $|\zeta|$ along the vertical line $\operatorname{Re} s = 0.6309$, obtained from integer counting data**. The spectral factorization $\zeta_\nu = \zeta_{\mathcal{L}} \cdot \zeta$ is here observed as a live transfer function, mode by mode. On this line all gains are nonzero, so the Cantor string is fully hearable, consistent with Lapidus–Maier ($\zeta$ has no zeros on $\operatorname{Re} s = D$ for this $D$ — indeed for any $D \neq \tfrac12$, conditionally on RH).

## 6. Experiment F3: the deaf channel at the critical line

To exhibit the Lapidus–Maier obstruction itself, take a *generalized* fractal string (a measure, in the sense of [3, Ch. 4]) with geometric counting function $N_{\mathcal{L}}(y) = y^{1/2} + \varepsilon \operatorname{Re}(y^{1/2 + i\gamma})$, $\varepsilon = 0.5$: geometry oscillating at complex dimension $\tfrac12 + i\gamma$. Computing $N_\nu(x) = \sum_m N_{\mathcal{L}}(x/m)$ by direct summation (no $\zeta$ anywhere in the forward model), subtracting the exact Weyl term, and measuring the surviving oscillation amplitude at frequency $\gamma$ in $\log x$:

| geometric oscillation | measured in spectrum | predicted $\varepsilon\lvert\zeta(\tfrac12+i\gamma)\rvert$ |
|---|---|---|
| $\gamma = 10.0$ (not a zero), amp. $0.50$ | **0.7747** | 0.7746 |
| $\gamma = 14.134725$ (first zeta zero), amp. $0.50$ | **0.0007** | 0.0000 |

(Sanity: the constant term of the residual returned $-1.4581$ against $\zeta(\tfrac12) = -1.4604$.) A fractal string can hide information from its own sound precisely at the zeta zeros; the estimator's blind spots, plotted against frequency, *are* the zeros. This is the inverse-spectral content of RH made visible in a table: the identifiability of dimension-$\tfrac12$ strings fails exactly on the zero set, and no estimator can be built that hears a frequency whose channel gain is exactly zero.

## 7. Scope: what this does and does not touch

This note is a methods contribution: to our knowledge, super-resolution (Prony-type) estimation of complex dimensions from geometric or spectral data has not been used in the fractal-strings literature, where complex dimensions are computed from analytically known $\zeta_{\mathcal{L}}$. The estimator gives the *inverse* direction — data to dimensions — with quantitative stability inherited from the pencil's perturbation theory, and F2 shows it resolves not only positions but the physically meaningful amplitudes.

It does **not** bear on RH. F3 is the reason, stated as an experiment: the zeros are the kernel of the geometry-to-spectrum map, and information annihilated by a forward map is unrecoverable by any estimator. What the combination legitimately supports is the program of Section 8.

## 8. The E7 program

**E7a — nonlattice strings: non-uniform damping made measurable.** For a self-similar string with incommensurable scaling ratios (e.g., $\zeta_{\mathcal{L}}(s) = (2^{-s}+3^{-s})/(1 - 2^{-s} - 3^{-s})$-type, dimension $D \approx 0.788$ solving $2^{-D} + 3^{-D} = 1$), the nonlattice theorem of [3] says the complex dimensions do **not** lie on one vertical line: their real parts fill an interval quasiperiodically. The pencil applied to such a string's data should measure a *drifting* damping column — the "RH-false" phenotype, on a system where it provably occurs. This exercises, on certified ground truth, exactly the detection capability that the estimator-criterion note's Theorem B invokes; the Diophantine fine structure of the drift (continued-fraction phenomena, per [3, Ch. 3]) gives a second layer of checkable predictions.

**E7b — Beurling control-group zetas.** Beurling generalized prime systems can be *designed* to have zeta functions with prescribed zero configurations, including zeros off the natural critical line (Diamond–Montgomery–Vorhauer [6]). Running the full E5 pipeline plus the validated protocol of the estimator-criterion note on such a system — where a rogue zero genuinely exists — would provide the only possible empirical test of the criterion's falsification branch, which can never be exercised on the true primes. Predicted outcome: the pencil reports a node with $\ln|z|/\Delta$ off the critical value, and the two-window validation gate fails along the Landau-type window sequence.

Both experiments are implementable with the exact tooling of this note (the forward models are sums over designed combs), and both convert "provable but abstract" dichotomies of the theory into instrument readings.

## Reproducibility

`cantor_pencil.py` (F1, F2; NumPy + mpmath for the reference $|\zeta|$ values only — no zeta values enter the forward data) and `deaf_channel.py` (F3), both in `code/`. Parameters as in Section 3; runtime under a minute on commodity hardware.

## References

1. M. L. Lapidus, C. Pomerance, The Riemann zeta-function and the one-dimensional Weyl–Berry conjecture for fractal drums, *Proc. London Math. Soc.* (3) **66** (1993) 41–69.
2. M. L. Lapidus, H. Maier, The Riemann hypothesis and inverse spectral problems for fractal strings, *J. London Math. Soc.* (2) **52** (1995) 15–34.
3. M. L. Lapidus, M. van Frankenhuijsen, *Fractal Geometry, Complex Dimensions and Zeta Functions: Geometry and Spectra of Fractal Strings*, 2nd Ed., Springer Monographs in Mathematics, Springer, New York, 2013.
4. H. Herichi, M. L. Lapidus, *Quantized Number Theory, Fractal Strings and the Riemann Hypothesis*, World Scientific, 2021.
5. Y. Hua, T. K. Sarkar, Matrix pencil method for estimating parameters of exponentially damped/undamped sinusoids in noise, *IEEE Trans. Acoust. Speech Signal Process.* **38** (1990) 814–824.
6. H. G. Diamond, H. L. Montgomery, U. M. A. Vorhauer, Beurling primes with large oscillation, *Math. Ann.* **334** (2006) 1–36.
7. E. C. Titchmarsh, *The Theory of the Riemann Zeta-Function*, 2nd Ed. (revised by D. R. Heath-Brown), Oxford Univ. Press, 1986.

---

¹ Ingeniero Electrónico, Universidad Tecnológica Nacional, Facultad Regional Buenos Aires, Buenos Aires, ARGENTINA. e-mail: ulimerlan@frba.utn.edu.ar. ORCID: 0009-0006-0259-3171
