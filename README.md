# Identification of the order distribution of linear systems by line-spectrum estimation on logarithmic spirals

**Author:** Ing. Ulises Merlan (Ingeniero Electrónico, Universidad Tecnológica
Nacional, Facultad Regional Buenos Aires) — `ulimerlan@frba.utn.edu.ar`

## Abstract

A linear time-invariant system whose transfer function admits the
distributed-order representation $H(s)=\int_{\mathbb{R}} \mathcal{A}(\alpha)\,
s^{\alpha}\, d\alpha$ is characterized by its *order distribution*
$\mathcal{A}$: integer-order systems correspond to Dirac combs supported on
$\mathbb{Z}$ whose weights are Laurent–Taylor coefficients, while
fractional-order systems place mass at non-integer $\alpha$'s. Sampling $H$ on
the unit circle and applying an FFT is a natural closed-form estimator of
$\mathcal{A}$, but it is subject to two structural obstructions: (i) a
truncation bias corrupting every coefficient whenever the response has a
non-zero steady state, and (ii) an aliasing theorem showing that any
unit-circle readout is $2\pi$-periodic, so non-integer orders can never
appear as spectral lines. Both defects are repaired: subtracting the
steady state reduces the integer-order problem to classical Cauchy-integral
coefficient extraction, and replacing the circle by a logarithmic spiral
turns order recovery into a line-spectrum estimation problem solvable by
matrix-pencil methods. As a calibration benchmark, the estimator recovers
the complex dimensions of the Cantor fractal string — known in closed form
— from geometric counting data alone, and, from spectral counting data,
recovers the same dimensions together with channel gains matching
independently computed values of the Riemann zeta function.

## Repository layout

```
paper/
  main.tex          LaTeX source
  main.pdf          Compiled PDF
code/
  order_spectrum.py   Reference implementation (NumPy) of the core
                      estimators; reproduces E1-E4
  cantor_pencil.py    Reproduces E5/F1 and F2 (Cantor string complex
                      dimensions, geometric + spectral)
  deaf_channel.py     Reproduces E5/F3 (the "deaf channel" at a zeta zero)
  prime_pencil.py,
  zeta_comb.py        Earlier, no-longer-referenced experiments (prime
                      counting / Riemann zeta zeros); kept for the record
  results.txt         Captured stdout
docs/
  notes/
    fractal-strings-pencil-note.md
      Methods note for the current E5 (Cantor string) calibration
      experiments referenced by the paper
    order-distribution-of-the-primes.md
      An earlier, unused calibration benchmark exploring the order
      distribution of the primes (prime counting / Riemann zeta zeros);
      kept here for the record but not cited in the paper
    staircase-order-distribution.md
      A methods note on the real-line identity relating an order comb to
      its cumulative (staircase) distribution; includes a documented
      negative result on why it does not transfer naively to the paper's
      finite circular-DFT pipeline (see code/order_spectrum.py's S1)
  figures/            Figures used while developing/validating E1-E4
                      (illustrative, not required to reproduce the
                      paper's numbers)
```

## Reproducing the numerical experiments

```zsh
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
python code/order_spectrum.py      # E1-E4
python code/cantor_pencil.py       # E5, F1 and F2
python code/deaf_channel.py        # E5, F3
```

`order_spectrum.py` and `deaf_channel.py`'s forward data are NumPy-only;
`cantor_pencil.py` and `deaf_channel.py` use `mpmath` only to compute the
independent reference values of $\zeta$ against which the estimator's
output is checked — no zeta values enter the estimator's input.

This reproduces, from closed-form and simulated data:

- **E1** — integer-order Taylor-coefficient recovery for
  $H(s)=1/(s^2+4s+5)$ (corrected vs. naive pipeline, max error $\sim10^{-8}$).
- **E2** — the aliasing signature of the branch cut $s^{-1/2}$ on the unit
  circle (energy spread over all integer bins, no line at $\alpha=-\tfrac12$).
- **E3** — spiral estimator recovery of exact fractional orders
  $\alpha=(-\tfrac12, 1.3)$ from noise-free and noisy samples.
- **E4** — data-driven fractional identification of
  $y(t)=1/\sqrt{\pi t}+2$ directly from time-domain samples.
- **E5** — a calibration benchmark from fractal geometry: complex
  dimensions of the Cantor string recovered from geometric data (F1) and
  from spectral data with channel gains matching $\zeta$ (F2), plus a
  "deaf channel" demonstration at a zeta zero (F3).

## Building the paper

Compiles with any standard TeX distribution (`amsmath`, `amssymb`,
`amsthm`, `booktabs`, `microtype`, `xcolor`, `geometry`, `hyperref`).
Bibliography is inline, no `.bib`/`.bbl` needed:

```zsh
cd paper
tectonic main.tex
```

## Citing

If you use this code or paper, please cite this repository:

```bibtex
@misc{merlan2026orderspectrum,
  author = {Ulises Merlan},
  title  = {Identification of the order distribution of linear systems
            by line-spectrum estimation on logarithmic spirals},
  year   = {2026},
  url    = {https://github.com/0D15E0/order-spectrum-estimation}
}
```

## License

- Code (`code/`): MIT License, see `LICENSE`.
- Paper text (`paper/`): © the author; licensed under CC BY 4.0
  (https://creativecommons.org/licenses/by/4.0/).


