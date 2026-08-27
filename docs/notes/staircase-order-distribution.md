# The Staircase Form of the Order Distribution: Integration as $1/(iu)$-Filtering

> **Integration notes.** Companion note to the main paper and to
> `order-distribution-of-the-primes.md`. Uses the same numbering convention
> (Def./Eq. references to the main paper's Section 2 unless noted). Not
> referenced by the paper; kept here as a methods note motivated by an
> identity worked out by hand (see `docs/notes/img/` if attached) and by the
> classical contrast between the prime *staircase* $\psi(x)$ and the prime
> *comb* $\Lambda(n)$ (Riemann–von Mangoldt), the same contrast this note
> makes precise and turns into an estimator variant.

## 1. The identity

Let $\mathcal A$ be a tempered distribution on $\mathbb R$ (an order
distribution in the sense of the paper's Def. 2.1) and let

$$
N(\alpha) \;=\; \int_{-\infty}^{\alpha} \mathcal A(\alpha')\, d\alpha'
$$

be its **cumulative distribution**: if $\mathcal A = \sum_m c_m\,
\delta(\alpha - \alpha_m)$ is a finite order comb (Eq. (4.10) of the paper),
$N$ is the right-continuous *staircase* that jumps by $c_m$ at each
$\alpha_m$ — for the prime system of the companion note, $\mathcal A =
\sum \Lambda(n)\,\delta(\alpha-\log n)$ (Eq. (C.1)) and $N$, in the
logarithmic variable, is exactly Chebyshev's $\psi$.

Since $N = \mathcal A * H\!\mathrm{eav}$ (convolution with the Heaviside
step) and the Fourier transform of the Heaviside function is
$\widehat{H\!\mathrm{eav}}(u) = \pi\delta(u) + \dfrac{1}{iu}$, the
convolution theorem gives, with $H(u) := \widehat{\mathcal A}(u)$ the
paper's transfer function evaluated on the unit circle $s=e^{iu}$
(Eq. (2.4)):

$$
\boxed{\ \widehat N(u) \;=\; \pi\, H(0)\, \delta(u) \;+\; \frac{H(u)}{iu}\ } .
\tag{S.1}
$$

*Derivation.* $\widehat{\mathcal A * H\!\mathrm{eav}} = \widehat{\mathcal
A}\cdot \widehat{H\!\mathrm{eav}} = H(u)\bigl(\pi\delta(u) + 1/(iu)\bigr) =
\pi H(0)\delta(u) + H(u)/(iu)$, using $H(u)\delta(u) = H(0)\delta(u)$. $\blacksquare$

On the spiral contour of Section 4 ($s = e^{(\sigma+i)u}$, $\log s =
(\sigma+i)u$), the same computation with $\log s$ playing the role of the
transform variable gives the version actually used for data:

$$
\widetilde N(\log s) \;=\; \pi\, H(0)\, \delta(\log s) \;+\; \frac{H(s)}{\log s} ,
\tag{S.1$'$}
$$

which is well defined and delta-free away from $s=1$ ($\log s = 0$), i.e.
on the entire spiral except its possible passage through the origin of the
$\log s$-plane.

## 2. Why this is not a new fact, but a useful one

**A notational warning first.** In Eq. (S.1), $H(0)$ denotes $H$ evaluated
at $u=0$, i.e. at $s=e^{i\cdot 0}=1$ — the **total order mass** $\int
\mathcal A(\alpha)\,d\alpha$ — and *not* the paper's steady state $H(s{=}0)
= y_\infty$ (Prop. 2.3), which sits at the *opposite* end of the contour
($s=0$, not $s=1$). The two are easy to conflate (both are informally "the
DC term") but are genuinely different numbers; the correct statement is
below.

**(a) It isolates a different, but analogous, singular term.**
Cor. 2.4's steady-state subtraction removes the pole that $y_\infty=H(0)$
(paper's $s=0$) would otherwise place at order $\alpha=-1$ after the
$s$-multiplication undoing the step input; Eq. (S.1)'s delta term instead
isolates the *total mass* of the comb, sitting at frequency $u=0$
($s=1$) of the *integration* operator itself, and is a structurally
separate correction, needed because $1/(iu)$ has its own pole exactly
there. The two corrections address different points of the same contour;
neither is a special case of the other, but both are instances of the same
general phenomenon: an operator with an isolated pole in the transform
variable requires the corresponding sample to be set aside before applying
it elsewhere on the contour.

**(b) It explains, structurally, why $\psi$ is tame and $\Lambda$ is not.**
Dividing by $iu$ is a first-order integrator on the *real line*: it
multiplies the continuous Fourier transform's amplitude at frequency $u$ by
$1/|u|$, attenuating high-frequency order content. This is the abstract
reason the classical explicit formula for $\psi(x)-x$ (companion note,
Eq. (C.2)) is a convergent series with a $\sqrt x$ envelope, while the
corresponding formula for $\Lambda$ as a sum of deltas has no pointwise
meaning at all: same identity, read in the two directions, **on the real
line** — the prime system's natural setting, since $\lambda=\log x$ ranges
over all of $\mathbb R$ and the sum over zeros in (C.2) is a genuine
(conditionally convergent) Fourier-type series, not a finite DFT.

**(c) A caution: this does *not* transfer naively to the paper's
finite-DFT (E1) setting.** For a *periodic* sample vector — the setting of
the paper's circle pipeline, where $u_k = 2\pi k/N$ and the Taylor comb is
read off a length-$N$ DFT — the correct discrete analogue of "integrate to
get a staircase" is a **circular** cumulative sum, whose Fourier multiplier
is $1/(1-e^{-iu})$ (not $1/(iu)$; the two agree only in the small-$u$,
continuum limit) restricted to nonzero bins, and even then a naive
per-bin application does *not* reproduce partial sums of $a_0,\dots,a_n$
correctly without additional boundary bookkeeping, because a finite
circular cumulative sum is only well defined once the sequence's mean
(the true $u=0$ bin) has been handled separately and the result is
inherently linear-in-$n$ if that mean is nonzero. We checked this directly
(§3) and confirm the naive translation fails as stated; the clean
real-line identity (S.1) should be tested on the prime system (real line,
not a period-$N$ circle), not retrofitted onto E1's finite Taylor comb.

**(d) The genuinely useful direction: read $N$ instead of $\mathcal A$ when
the *data itself* is a staircase.** When the measured quantity already is
a counting function or a step response — exactly the prime system's
$\psi(x)$, or any system's step response before differentiation — no
division step is needed at all: transform the staircase directly
(Eq. (4.9)), apply the matrix pencil to recover $\{\rho, c_\rho\}$ from
$N$ itself (as `order-distribution-of-the-primes.md` already does via
Prop. C.1's Gaussian-smoothed $\psi$), and differentiate the recovered
model, if the comb form is wanted, only at the end, in closed form. This
sidesteps resolving deltas from noisy samples entirely — the practical
payoff of the "read the staircase" idea, decoupled from the (incorrect,
see (c)) attempt to derive it as a circular-DFT postprocessing step.

## 3. Numerical check

`code/order_spectrum.py`'s `staircase_transform` implements the continuous
multiplier $H(u)/(iu)$ (with the $u=0$ sample set aside, returned
separately). `experiment_S1` applies it, naively, to E1's circular DFT
samples of $H(s)=1/(s^2+4s+5)$ and compares the result against partial
sums of the known Taylor coefficients: **the comparison fails**, exactly as
flagged in §2(c) above — the printed output includes the mismatch and a
diagnostic using the correct circular kernel $1/(1-e^{-iu})$, which *also*
fails to reproduce partial sums without further boundary corrections, since
$\operatorname{Re}\bigl[1/(1-e^{-iu})\bigr]=\tfrac12$ for every nonzero
$u$, so applying it uniformly just returns (up to a factor) a rescaling of
the original comb, not its cumulative sum. This negative result is kept
deliberately, as a documented pitfall: the identity (S.1) is correct and
useful in its native, real-line setting (§2(b), and the prime system of
the companion note), but does not transplant to a finite periodic DFT by
naive per-bin division.

## 4. Status

This is a real-line Fourier identity (rigorously established in §1) with a
correct qualitative explanation (§2(b)) for the prime staircase versus
comb contrast, and a clearly negative result (§2(c)/§3) warning against a
naive finite-DFT analogue. No new estimator or experiment is validated
here. The natural next step, if pursued, is a direct real-line check using
`order-distribution-of-the-primes.md`'s existing data: compare a
pencil-on-$\psi$ read-out against a pencil-on-$\Lambda$ read-out (P1) of
the *same* sieve data and confirm they return the same zeros $\rho$, which
would be a genuine (not naively-DFT) numerical confirmation of Eq. (S.1)'s
content.
