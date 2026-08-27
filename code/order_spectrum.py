"""
Order-spectrum identification of LTI systems from time-domain data.

Implements the corrected pipeline derived from U. Merlan's handwritten notes
(continuous-order Laurent expansion H(s) = int A(alpha) s^alpha d alpha):

1. circle_taylor_coeffs : corrected unit-circle method. From a measured step
   response y(t), subtract the steady state, compute the Laplace transform of
   the residual on the unit circle s = e^(iu), rebuild H(e^(iu)) and FFT to
   obtain the Taylor coefficients of H about s = 0 (the integer-order comb
   weights of A(alpha)). Equivalent to Lyness-Moler (1967) coefficient
   extraction applied to H(s), reconstructed from data as
   y_inf + s*G_A(s) with G_A the truncated transform of y - y_inf.

2. naive_circle_coeffs : the original (uncorrected) pipeline, kept only to
   demonstrate the truncation-bias failure when y(inf) != 0.

3. matrix_pencil : Hua-Sarkar matrix-pencil estimation of exponential sums
   f_j = sum_m c_m z_m^j from uniform samples.

4. spiral_order_estimate : the log-spiral estimator. Sample the transform on
   s = e^{(sigma+i)u} (an arc in the right half-plane when working from
   data); the order content A(alpha) = sum_m c_m delta(alpha - alpha_m)
   becomes an exponential sum in u, and the matrix pencil returns the
   (possibly NON-INTEGER) orders alpha_m and weights c_m. This breaks the
   2*pi-periodicity aliasing that makes the unit circle blind to fractional
   orders in measured data.

5. staircase_transform : implements the real-line identity (docs/notes/
   staircase-order-distribution.md) relating the order comb A(alpha) to its
   cumulative distribution N(alpha) via division by i*u in the Fourier
   variable. experiment_S1 documents, deliberately, why this identity does
   NOT transfer by naive per-bin division to the finite circular DFT used
   in E1 -- see the note for the correct (real-line) setting.

Run "python order_spectrum.py" to reproduce all experiments (E1-E4) reported
in the accompanying paper, plus the staircase-identity note (S1).
"""
import numpy as np

# quadrature helpers

def simpson_uniform(f, h):
    """Composite Simpson rule for samples f on a uniform grid of step h.
    len(f) must be odd."""
    f = np.asarray(f)
    n = f.shape[-1]
    if n % 2 == 0:
        raise ValueError("Simpson rule needs an odd number of samples")
    w = np.ones(n)
    w[1:-1:2] = 4.0
    w[2:-1:2] = 2.0
    return (h / 3.0) * (f @ w)


def laplace_from_samples(t, y, s_points):
    """Truncated Laplace transform int_0^A y(t) e^{-s t} dt by Simpson rule.
    t: uniform grid on [0, A] with an odd number of points.
    y: samples of y on t.
    s_points: array of complex evaluation points.
    Returns complex array, one value per s."""
    if len(t) % 2 == 0:
        raise ValueError("laplace_from_samples needs an odd number of samples")
    h = t[1] - t[0]
    if not np.allclose(np.diff(t), h):
        raise ValueError("laplace_from_samples requires a uniform time grid")
    out = np.empty(len(s_points), dtype=complex)
    for k, s in enumerate(s_points):
        out[k] = simpson_uniform(y * np.exp(-s * t), h)
    return out


# corrected unit-circle pipeline (integer-order comb / Taylor coefficients)
def circle_taylor_coeffs(t, y, n_coeffs, n_circle=64, y_inf=None):
    """Taylor coefficients a_0..a_(n_coeffs-1) of H(s) about s=0, from the
    measured STEP response y(t) of the system H.
    Key correction versus the naive method: the steady state y_inf = H(0) is
    subtracted before transforming, so the integrand decays and the transform
    converges on the whole circle. Then
    G = Laplace of (y - y_inf),
    Y(s) = y_inf / s + G(s),
    H(s) = s Y(s) = y_inf + s G(s),
    and the FFT of H(e^{iu}) over the circle returns the Taylor comb.
    If y_inf is None it is estimated from the mean of the last 2% of samples.
    Returns (real part of the coefficients, max |Im|); the imaginary parts
    are a sanity check and should be ~0 for a real system."""
    if y_inf is None:
        tail = max(2, int(0.02 * len(y)))
        y_inf = float(np.mean(y[-tail:]))
    g = y - y_inf
    u = 2.0 * np.pi * np.arange(n_circle) / n_circle
    s = np.exp(1j * u)
    G = laplace_from_samples(t, g, s)
    H = y_inf + s * G
    a = np.fft.fft(H) / n_circle
    # a_n = (1/N) sum_k H_k e^{-2pi i nk/N}
    coeffs = a[:n_coeffs]
    return np.real(coeffs), float(np.max(np.abs(np.imag(coeffs))))


def naive_circle_coeffs(t, y, n_coeffs, n_circle=64):
    """The UNCORRECTED pipeline of the original notes: transform the raw step
    response on the circle and shift by one order (multiply by s) to undo the
    1/s of the step input. Included only to exhibit the truncation bias
    -y_inf * (-A)^(m+1) / (m+1)! that corrupts every bin when y_inf != 0."""
    u = 2.0 * np.pi * np.arange(n_circle) / n_circle
    s = np.exp(1j * u)
    Y = laplace_from_samples(t, y, s)
    H = s * Y
    a = np.fft.fft(H) / n_circle
    return np.real(a[:n_coeffs])


# 2. matrix-pencil estimation of exponential sums (Hua & Sarkar 1990)
def matrix_pencil(f, n_modes, pencil=None):
    """Estimate {z_m, c_m} such that f_j = sum_m c_m z_m^j, j = 0..N-1.
    f: complex samples on a uniform grid.
    n_modes: number of exponentials M to fit.
    pencil : pencil parameter L (default N//2; N/3 <= L <= N/2 recommended).
    Returns (Z, C, singular_values)."""
    f = np.asarray(f, dtype=complex)
    N = len(f)
    if n_modes < 1:
        raise ValueError("n_modes must be >= 1")
    L = pencil if pencil is not None else N // 2
    if not (n_modes <= L <= N - n_modes):
        raise ValueError(
            f"pencil parameter L={L} outside [{n_modes}, {N - n_modes}] "
            f"for N={N}, M={n_modes}")
    # Hankel data matrix, (N-L) x (L+1)
    Y = np.array([f[i:i + L + 1] for i in range(N - L)])
    _, sv, Vh = np.linalg.svd(Y, full_matrices=False)
    # NB: rows of Y live in span{(z_m^j)_j}; for complex data that span equals
    # span of the CONJUGATES of the right singular vectors, so take Vh.T
    # (transpose without conjugation) to get the shift-invariant subspace.
    V = Vh.T[:, :n_modes]
    V1, V2 = V[:-1, :], V[1:, :]
    z = np.linalg.eigvals(np.linalg.pinv(V1) @ V2)
    # amplitudes by Vandermonde least squares
    Z = z[None, :] ** np.arange(N)[:, None]
    c, *_ = np.linalg.lstsq(Z, f, rcond=None)
    return z, c, sv


def spiral_order_estimate(F, sigma, du, n_modes, u0=0.0, pencil=None):
    """Recover orders alpha_m and weights c_m of F(s) = sum_m c_m s^{alpha_m}
    from samples F_j = F(e^{(sigma+i) u_j}), u_j = u0 + j du (log-spiral).
    on the spiral, F_j = sum_m [c_m e^{(sigma+i) u0 alpha_m}] *
    (e^{(sigma+i) du alpha_m})^j is an exponential sum in j; the matrix pencil
    returns z_m = e^{(sigma+i) du alpha_m}, and alpha_m is read back from
    log z_m. Because sigma != 0, ln|z_m| = sigma du alpha_m is single-valued,
    which resolves the 2*pi branch ambiguity of the imaginary part for ANY
    real alpha (this is precisely the degree of freedom the unit circle,
    sigma = 0, lacks): the modulus selects the branch, then the full
    two-component formula, which is better conditioned, is applied on the
    selected branch.
    Returns (alpha, c, singular_values); alpha should be ~real for a real
    order spectrum -- the imaginary part is a diagnostic."""
    if sigma == 0:
        raise ValueError("spiral estimator requires sigma != 0")
    z, c_raw, sv = matrix_pencil(F, n_modes, pencil)
    # branch selection via the single-valued modulus, then the two-component
    # formula on the selected branch (exact: k recovers the winding number)
    alpha_mod = np.log(np.abs(z)) / (sigma * du)
    k = np.round((du * alpha_mod - np.imag(np.log(z))) / (2.0 * np.pi))
    alpha = (np.log(z) + 2j * np.pi * k) / ((sigma + 1j) * du)
    c = c_raw / np.exp((sigma + 1j) * u0 * alpha)
    order = np.argsort(alpha.real)
    return alpha[order], c[order], sv


def staircase_transform(H_samples, u, sigma=0.0):
    """Staircase-form filter: given samples H(u_j) of the order-distribution
    transfer function (Eq. 2.4) on the circle (sigma=0) or on the spiral's
    log s = (sigma + i) u (Eq. 4.9), return samples of the Fourier/Mellin
    transform of the CUMULATIVE order distribution N(alpha) =
    int_{-inf}^alpha A(alpha') d alpha', via the identity (docs/notes/
    staircase-order-distribution.md, Eq. S.1/S.1'):

        N^(u) = pi H(0) delta(u) + H(u) / (i u)          [circle, sigma=0]
        N~(log s) = pi H(0) delta(log s) + H(s) / log s  [spiral]

    i.e. dividing by (sigma + i) u is an exact integrator in the order
    variable: it turns a comb (deltas) into a staircase (jumps), the two
    being related by the classical Heaviside/(1/(iu)) Fourier pair. The
    u=0 sample, if present, carries the delta weight pi*H(0) and cannot be
    divided; it is returned separately.

    Parameters
    ----------
    H_samples : complex array, H evaluated at s_j = e^{(sigma+i) u_j}.
    u : real array of the same length, the u_j (NOT reduced mod 2*pi for the
        spiral branch; on the circle this is the usual DFT grid).
    sigma : spiral damping parameter (0.0 for the plain unit circle).

    Returns
    -------
    (dc_weight, N_samples) : dc_weight is pi*H(0) read off the u=0 sample if
        present (else 0j); N_samples has the same length as the input, with
        the u=0 entry (if any) set to 0 (its content is entirely in
        dc_weight, matching the delta term of Eq. S.1)."""
    H_samples = np.asarray(H_samples, dtype=complex)
    u = np.asarray(u, dtype=float)
    log_s = (sigma + 1j) * u
    N_samples = np.zeros_like(H_samples)
    dc_weight = 0.0 + 0.0j
    zero_mask = (u == 0.0) if sigma == 0.0 else np.zeros_like(u, dtype=bool)
    nz = ~zero_mask
    N_samples[nz] = H_samples[nz] / log_s[nz]
    if np.any(zero_mask):
        dc_weight = np.pi * H_samples[zero_mask][0]
    return dc_weight, N_samples


# experiments

def _fmt(x, d=12):
    return f"{x:+.{d}f}"


def experiment_E1():
    """Integer-order system from step data: corrected circle pipeline."""
    print("\n" + "=" * 78)
    print("E1  H(s) = 1/(s^2+4s+5), step response, corrected unit-circle pipeline")
    print("=" * 78)
    A, Nt, Nu = 12.0, 65537, 64
    t = np.linspace(0.0, A, Nt)
    y = 0.2 - 0.2 * np.exp(-2 * t) * (2 * np.sin(t) + np.cos(t))
    exact = [1 / 5,
             -4 / 25, 11 / 125, -24 / 625, 41 / 3125, -44 / 15625,
             -29 / 78125]
    est, max_imag = circle_taylor_coeffs(t, y, len(exact), Nu)  # y_inf estimated
    est_k, _ = circle_taylor_coeffs(t, y, len(exact), Nu, y_inf=0.2)  # known y_inf
    naive = naive_circle_coeffs(t, y, len(exact), Nu)
    print(f"window A={A}, time samples Nt={Nt}, circle samples Nu={Nu}")
    print(f"steady state estimated from tail; max |Im a_n| = {max_imag:.2e}\n")
    print(f"{'n':>2} {'exact a_n':>18} {'corrected (est y_inf)':>22} "
          f"{'corrected (known y_inf)':>24} {'naive pipeline':>18}")
    for n, (e, v, vk, nv) in enumerate(zip(exact, est, est_k, naive)):
        print(f"{n:>2} {_fmt(e):>18} {_fmt(v):>22} {_fmt(vk):>24} {nv:>18.6g}")
    err = np.max(np.abs(np.array(est) - np.array(exact)))
    err_k = np.max(np.abs(np.array(est_k) - np.array(exact)))
    print(f"\nmax abs error, corrected (estimated y_inf): {err:.3e}")
    print(f"max abs error, corrected (known y_inf): {err_k:.3e}")
    # the naive method multiplies by s (bin shift), so the m=0 truncation bias
    # -y_inf(-A)^{m+1}/(m+1)! lands at bin 1:
    print("naive-pipeline excess at bin 1: predicted "
          f"{-0.2 * (-A):.3f}, observed {naive[1] - exact[1]:.3f}")
    return est, err


def experiment_E2():
    """Aliasing: unit-circle readout of s^{-1/2} gives an integer comb."""
    print("\n" + "=" * 78)
    print("E2  Aliasing obstruction: F(s)=s^{-1/2} sampled on the unit circle")
    print("-" * 78)
    Nu = 64
    u = 2.0 * np.pi * np.arange(Nu) / Nu  # u in [0, 2pi)
    F = np.exp(-0.5j * u)  # branch s^{-1/2} = e^{-iu/2}
    bins = np.fft.fft(F) / Nu
    print("first bins of the circle FFT vs exact Fourier coefficients of the")
    print("single-valued restriction, c_m = -i/(pi (m+1/2)) aliased over m+1*Nu:")
    print(f"{'m':>3} {'|FFT bin|':>14} {'|c_m| = 1/(pi|m+1/2|)':>24}")
    for m in [0, 1, 2, 3, 5, 10, 20, 31]:
        cm = 1.0 / (np.pi * abs(m + 0.5))
        print(f"{m:>3} {abs(bins[m]):>14.6f} {cm:>24.6f}")
    print("\n-> energy spread over ALL integer bins with ~1/(m+1/2) decay;\n"
          "no line at alpha = -1/2 can ever appear: the circle transform of\n"
          "any data is exactly 2pi-periodic, so its FFT support is Z.")
    return bins


def experiment_E3():
    """Spiral estimator on an exact fractional formula."""
    print("\n" + "-" * 78)
    print("E3  Log-spiral + matrix pencil, exact samples of "
          "H(s) = s^{-1/2} + 0.5 s^{1.3}")
    print("=" * 78)
    sigma, du = 0.2, 0.3
    true_alpha = np.array([-0.5, 1.3])
    true_c = np.array([1.0, 0.5])

    def H_on_spiral(j):
        uj = j * du
        s_log = (sigma + 1j) * uj  # log s along the spiral
        return (true_c[0] * np.exp(true_alpha[0] * s_log)
                + true_c[1] * np.exp(true_alpha[1] * s_log))

    for N, noise in [(8, 0.0), (64, 1e-8), (64, 1e-4)]:
        F = H_on_spiral(np.arange(N))
        if noise > 0:
            rng = np.random.default_rng(1)
            F = F * (1 + noise * (rng.standard_normal(N)
                                   + 1j * rng.standard_normal(N)))
        alpha, c, sv = spiral_order_estimate(F, sigma, du, 2)
        print(f"\nN={N:>3}, relative sample noise {noise:g}:")
        for a_t, c_t, a_e, c_e in zip(true_alpha, true_c, alpha, c):
            print(f"  true alpha={a_t:+.1f}, c={c_t:.1f}  ->  "
                  f"alpha = {a_e.real:+.12f} {a_e.imag:+.1e}i, "
                  f"c = {c_e.real:+.10f}")
        print("  singular values of pencil matrix: "
              + ", ".join(f"{x:.2e}" for x in sv[:4]) + ", ...")
        print(f"  order-position errors: "
              f"{np.max(np.abs(alpha.real - true_alpha)):.2e}")


def experiment_E4():
    """Spiral estimator from TIME-DOMAIN data with a fractional component."""
    print("\n" + "=" * 78)
    print("E4  Data-driven fractional identification: y(t)=1/sqrt(pi*t) + 2,")
    print("    Y(s) = s^{-1/2} + 2 s^{-1}; spiral arc in the right half-plane")
    print("-" * 78)
    # Laplace transform of the data, computed from samples with the
    # substitution t = v^2 (removes the t^{-1/2} singularity):
    # int_0^A y e^{-s t} dt = int_0^sqrt(A) (2/sqrt(pi) + 4v) e^{-s v^2} dv
    A, Nv, N = 120.0, 32769, 64
    v = np.linspace(0.0, np.sqrt(A), Nv)
    base = 2.0 / np.sqrt(np.pi) + 4.0 * v  # smooth integrand factor
    hv = v[1] - v[0]

    u_lo, u_hi = -1.35, 1.35
    u = np.linspace(u_lo, u_hi, N)
    du = u[1] - u[0]
    sigma = 0.15
    s = np.exp((sigma + 1j) * u)  # arc; Re s >= 0.17

    true_alpha = np.array([-1.0, -0.5])
    true_c = np.array([2.0, 1.0])

    def transform(noise=0.0, seed=0):
        w = base.copy()
        if noise > 0:  # relative noise on y samples
            rng = np.random.default_rng(seed)
            # perturb y(t(v)); base = y(t(v)) * 2v with y = 1/sqrt(pi t)+2
            y_tv = np.empty_like(v)
            # placeholder; the v=0 integrand limit is set below (y(0+) is
            # infinite and is not a measurable sample, so it carries no noise)
            y_tv[0] = 0.0
            y_tv[1:] = 1.0 / np.sqrt(np.pi) / v[1:] + 2.0
            y_noisy = y_tv * (1 + noise * rng.standard_normal(Nv))
            w = 2.0 * v * y_noisy
        w[0] = 2.0 / np.sqrt(np.pi)  # limit value at v=0
        F = np.empty(N, dtype=complex)
        for k in range(N):
            F[k] = simpson_uniform(w * np.exp(-s[k] * v * v), hv)
        return F

    for noise in [0.0, 1e-4]:
        F = transform(noise=noise)
        alpha, c, sv = spiral_order_estimate(F, sigma, du, 2, u0=u_lo)
        print(f"\nmeasurement noise (relative, on y samples): {noise:g}")
        print("  singular values of pencil matrix: "
              + ", ".join(f"{x:.2e}" for x in sv[:4]) + ", ...")
        for a_t, c_t, a_e, c_e in zip(true_alpha, true_c, alpha, c):
            print(f"  true alpha={a_t:+.1f}, c={c_t:.1f}  ->  "
                  f"alpha = {a_e.real:+.8f} {a_e.imag:+.1e}i, "
                  f"c = {c_e.real:+.8f}")
        print(f"  order-position errors: "
              f"{np.max(np.abs(alpha.real - true_alpha)):.2e}")


def experiment_S1():
    """Check of the staircase identity (docs/notes/
    staircase-order-distribution.md, Eq. S.1) against E1's finite circular
    DFT setting. THIS IS A DOCUMENTED NEGATIVE RESULT, kept deliberately: the
    identity N^(u) = pi H(0) delta(u) + H(u)/(i u) is a REAL-LINE Fourier
    fact (proved in the note), and does not transfer to a length-N circular
    DFT by naive per-bin division. Two discrete candidates are tried below:
    (i) the naive continuum kernel 1/(i u_k), and (ii) the exact circular
    cumulative-sum kernel 1/(1 - e^{-i u_k}). Neither reproduces the partial
    sums of the known Taylor coefficients, confirming the note's warning
    (Section 2c): (i) has the wrong small-u normalization for a length-N
    sequence, and (ii) has constant real part 1/2 for every nonzero u_k (a
    classical cotangent identity), so applying it uniformly just rescales
    the original comb rather than integrating it -- a circular cumulative
    sum additionally requires removing the sequence's own mean (its true
    u=0 bin) before the kernel is meaningful, and even then involves
    boundary terms beyond a single Fourier multiplier."""
    print("\n" + "=" * 78)
    print("S1  Staircase identity (Eq. S.1): NEGATIVE result on E1's circular DFT")
    print("=" * 78)
    A, Nt, Nu = 12.0, 65537, 64
    t = np.linspace(0.0, A, Nt)
    y = 0.2 - 0.2 * np.exp(-2 * t) * (2 * np.sin(t) + np.cos(t))
    exact_a = np.array([1 / 5,
                         -4 / 25, 11 / 125, -24 / 625, 41 / 3125,
                         -44 / 15625, -29 / 78125])
    y_inf = 0.2
    g = y - y_inf
    u = 2.0 * np.pi * np.arange(Nu) / Nu
    s = np.exp(1j * u)
    G = laplace_from_samples(t, g, s)
    H = y_inf + s * G

    dc_weight, N_naive = staircase_transform(H, u, sigma=0.0)
    bins_naive = np.fft.ifft(N_naive) * Nu

    circ_kernel = np.zeros(Nu, dtype=complex)
    circ_kernel[1:] = 1.0 / (1.0 - np.exp(-1j * u[1:]))
    bins_circ = np.fft.ifft(H * circ_kernel) * Nu

    partial_sums = np.cumsum(exact_a)
    n_show = len(exact_a)
    print(f"pi*H(u=0) = {dc_weight:.6f}  (total comb mass sum a_n over the "
          f"circle's period, NOT the paper's y_inf)")
    print(f"\n{'n':>2} {'partial sum a_0..a_n':>22} "
          f"{'naive 1/(iu) bin':>18} {'circular 1/(1-e^-iu) bin':>26}")
    for n in range(n_show):
        print(f"{n:>2} {partial_sums[n]:>22.10f} "
              f"{bins_naive[n].real:>18.6f} {bins_circ[n].real:>26.6f}")
    print("\n-> neither column tracks the partial sums: as documented in\n"
          "staircase-order-distribution.md Section 2(c)/3, Eq. S.1 is a\n"
          "real-line identity and does not transfer to this finite circular\n"
          "DFT setting by naive per-bin division. Kept as a documented\n"
          "pitfall, not a validated experiment.")


if __name__ == "__main__":
    np.set_printoptions(precision=12, suppress=False)
    experiment_E1()
    experiment_E2()
    experiment_E3()
    experiment_E4()
    experiment_S1()
