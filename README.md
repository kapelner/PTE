Personalized Treatment Evaluator
===

Code for the R package PTE as well as code to duplicate the figures found in 

Kapelner, A, Bleich, J, Cohen, ZD, DeRubeis, RJ and Berk, R (2014) Inference for Treatment Regime Models in Personalized Medicine, arXiv

## Why upgrade to v2.0

v2.0 replaces the R-level bootstrap x cross-validation loop -- the part of
`PTE_bootstrap_inference` that dominates runtime, since it runs thousands of times per call -- with
compiled Rcpp/RcppEigen routines for all three default regression types. The default path now keeps
bootstrap sampling, fold fitting, prediction, and q-score construction in native code, while leaving
the statistical result unchanged (verified against R's own `lm()`/`glm()`/`survreg()` to numerical
precision, plus 1000+ randomized correctness checks against `survival::survfit()` for the
Kaplan-Meier median computation). The timings below were measured sequentially (`num_cores = 1`)
with a native-speed build equivalent to `PTE_NATIVE_SPEED=1 R CMD INSTALL PTE`. They report median
wall-clock time per bootstrap iteration over three runs; the R column uses equivalent R formula
model builders from the v1.7 R-level path, while the Rcpp column uses the default v2.0 path. Both
columns use the same bootstrap resample size: `m_pow_of_n = 1` in v2.0, equivalent to `m_prop = 1`
in v1.7, so each bootstrap sample has `m = n` rows:

| `regression_type` | Default model | R formula path | v2.0 native Rcpp | Speedup |
|---|---|---|---|---|
| `"continuous"`  | OLS (`lm`)               | ~68 ms/iter | ~0.081 ms/iter | **~844x** |
| `"incidence"`   | Logistic (`glm`)         | ~32 ms/iter | ~0.228 ms/iter | **~140x** |
| `"survival"`    | Weibull AFT (`survreg`)  | ~90 ms/iter | ~0.824 ms/iter | **~110x** |

One "iter" is a full bootstrap sample: sampling rows, fitting + predicting every cross-validation
fold, and computing the default q-scores. Continuous and survival use the package's bundled example
datasets; incidence uses the deterministic test fixture because the package does not ship a separate
incidence example dataset. Your own data's speedup will depend on model difficulty, sample size,
predictor count, and whether you install with native CPU flags.

On top of the fit itself, the bootstrap's parallel backend now forks workers (instead of spawning
and serializing data to separate processes) and dispatches one chunk of iterations per worker
instead of one task per bootstrap sample, so wall-clock scaling with `num_cores` is also
substantially better than v1.7's.

## Installing for maximum speed

By default the package compiles portably (no machine-specific CPU flags), which is required for CRAN and safe to redistribute. If you're installing from source on the machine you'll actually run it on and want the fastest possible build, opt into native-CPU optimization:

```
PTE_NATIVE_SPEED=1 R CMD INSTALL PTE
```

This compiles with `-march=native -mtune=native -O3`, tying the resulting binary to the exact CPU that compiled it -- don't use this for a build you plan to share or distribute. Link-time optimization is a separate, independent opt-in (`PTE_NATIVE_LTO=1`) since it has shown non-deterministic link failures on some toolchains; combine both with `PTE_NATIVE_SPEED=1 PTE_NATIVE_LTO=1 R CMD INSTALL PTE` if you've confirmed it's stable on your machine.

**How much does the flag actually buy you?** A naive "install once with the flag, once without"
comparison is misleading on some machines: this development machine's R build has `-O3
-march=native` baked directly into R's own site-wide `Makeconf` (`R CMD config CXXFLAGS`), so a
plain `R CMD INSTALL` *already* gets those flags regardless of `PTE_NATIVE_SPEED`, and neither a
`CXX20FLAGS` environment variable nor `MAKEFLAGS` override it (GNU Make's own variable-precedence
rules let the Makeconf assignment win). To get a real apples-to-apples number, the two fit kernels
were instead compiled standalone with plain `g++` -- fully bypassing R's build system and any site
defaults -- once with `-O2` (a true generic baseline) and once with `-O3 -march=native
-mtune=native`, then loaded and timed directly. Since `m_pow_of_n = 0.75` is now the default (see
below), the per-fold training matrices this exercises are much smaller than the full dataset --
`m = round(n^0.75)` rows per bootstrap replicate, then ~90% of *that* per cross-validation fold
(`n = 154` -> `m = 44` -> ~40 training rows for continuous/incidence; `n = 360` -> `m = 83` -> ~75
training rows for survival) -- so these numbers supersede an earlier version of this table that
benchmarked at the full, un-subsampled `n`:

| Fit kernel | Generic (`-O2`) | `-O3 -march=native -mtune=native` | Speedup |
|---|---|---|---|
| Logistic (incidence fold fit, ~40 rows) | ~0.0045 ms/call | ~0.0038 ms/call | ~1.19x |
| Weibull AFT (survival fold fit, ~75 rows) | ~0.0538 ms/call | ~0.0395 ms/call | ~1.36x |

(Continuous/OLS isn't listed: its default fit path is R's own `.lm.fit()`, not custom compiled
code, so it isn't affected by this package's compile flags either way.) The gain is real but more
modest than it was at full `n` -- for the Weibull AFT fit, `-march=native`'s wider SIMD registers
(confirmed here as `__m256d`/AVX2 vs `__m128d`/SSE2 in the generated code) have proportionally less
raw linear algebra to amortize against per-call fixed overhead once the training matrices shrink
this much under the m-out-of-n default, so its speedup dropped from ~1.84x (at full `n`) to ~1.36x
here; the logistic fit's speedup was essentially unchanged (~1.16x -> ~1.19x, within run-to-run
noise) since it was already a smaller, more overhead-dominated computation even at full `n`.
Whether *your* plain `R CMD INSTALL` already gets some of this for free depends entirely on your
own R installation's default `CXXFLAGS`; check `R CMD config CXXFLAGS` to see which case you're in
before assuming `PTE_NATIVE_SPEED=1` will move the needle by this much on top of your default
build -- and note that a larger `m_pow_of_n` (up to `1`, the classical bootstrap) or a larger `n`
will exercise bigger per-fold matrices and should recover more of the full-`n` speedup shown above.

## The m-out-of-n bootstrap: why, and how to tune it

### The estimand is non-smooth

`PTE_bootstrap_inference`'s reported quantities (`q_adversarial`, `q_average`, `q_best`) are built from
a **treatment recommendation**, not directly from a smooth model parameter. For each held-out subject,
the fitted model produces `est_true` (predicted outcome under the treatment actually given) and
`est_counterfactual` (predicted outcome under the other treatment), and the recommendation is a hard
threshold on their difference:

```r
optimal = est_true > est_counterfactual   # a 0/1 indicator, not a continuous quantity
```

`q_average`, `q_adversarial` and `q_best` are all means computed over subsets defined by this
indicator. An indicator/threshold applied to an estimated quantity is the textbook example of a
**non-smooth (non-Hadamard-differentiable) functional** -- the same kind of object that shows up in
individualized-treatment-rule/policy-value estimation and in other argmax-type statistics. It has a
real consequence: whenever a non-negligible fraction of subjects sit near the decision boundary (an
estimated treatment effect close to zero -- i.e. patients genuinely close to indifferent between the
two arms), small perturbations from resampling flip their recommendation, the functional's derivative
doesn't exist there, and the ordinary bootstrap's usual delta-method justification breaks down.

### Why the ordinary bootstrap can fail here

The classical **n-out-of-n bootstrap** (resample n rows with replacement, size-for-size with the
original data) is consistent for smooth functionals via the delta method, but is not generally
consistent for statistics built from indicators/argmax near a non-regular boundary like this one. In
that regime the bootstrap distribution doesn't reliably approximate the true sampling distribution of
`q_average`/`q_adversarial`/`q_best`, even when the point estimate itself is still fine -- so a naive
bootstrap CI or p-value can be misleading precisely in the clinically interesting case where a
meaningful fraction of patients are near-indifferent between treatments.

### The fix: the m-out-of-n bootstrap

Politis, Romano & Wolf (*Subsampling*, Springer, 1999) show that resampling a smaller number of
observations, `m < n`, restores bootstrap consistency for non-regular functionals like this one, as
long as `m -> infinity` and `m / n -> 0` as `n` grows. This package implements that fix directly:
`PTE_bootstrap_inference`'s `m_pow_of_n` argument sets each bootstrap replicate's resample size
(sampled *with* replacement, i.e. the "m-out-of-n bootstrap" as opposed to without-replacement
subsampling) to

```
m = round(n ^ m_pow_of_n)
```

`m_pow_of_n = 1` recovers the classical, here potentially inconsistent, `n`-out-of-`n` bootstrap.
Whenever `m_pow_of_n < 1`, the reported confidence intervals (both the percentile and basic methods)
are additionally rescaled by `sqrt(m / n)`, as m-out-of-n bootstrap theory requires: the bootstrap
replicates approximate `sqrt(m) * (T*_m - T_n)`, not `sqrt(n) * (T*_m - T_n)`, so using the raw
bootstrap quantiles without this correction would produce intervals of the wrong width. This
rescaling is a no-op when `m_pow_of_n = 1`.

### The default: `m_pow_of_n = 0.75`

`PTE_bootstrap_inference` defaults to `m_pow_of_n = 0.75` (i.e. `m = n^0.75`), the rate Politis,
Romano & Wolf recommend as a practical compromise: large enough that the bootstrap's own Monte Carlo
variance doesn't dominate, small enough to restore consistency near the non-regular boundary. There is
no single population-optimal `m_pow_of_n`, though -- it's dataset-dependent, since it trades off against
how much probability mass actually sits near the decision boundary for your data and model. `0.75` is a
good default, not a guarantee.

### Choosing `m_pow_of_n` for your own data: `select_optimal_m_prop()`

To calibrate `m_pow_of_n` to your specific dataset rather than relying on the default, use
`select_optimal_m_prop()`. It implements the **minimum volatility (MV) method** (Politis, Romano &
Wolf, 1999, Ch. 9; refined by Bickel & Sakov, 2008, *Statistica Sinica* 18(3)): it runs
`PTE_bootstrap_inference` once per candidate `m_pow_of_n` on a grid, measures how much the resulting
confidence interval width varies over a small centered neighborhood of adjacent grid points, and picks
the grid point whose neighborhood is most stable -- the region of the grid where the bootstrap's answer
has stopped reacting to further tuning of `m`.

```r
library(PTE)
data(continuous_example)
X = continuous_example$X
y = continuous_example$y

m_selection = select_optimal_m_prop(
    X, y,
    regression_type = "continuous",
    m_pow_of_n_grid = seq(0.5, 1, by = 0.05),
    B = 250,
    num_cores = 1
)
m_selection$grid_table
```

```
 m_pow_of_n   m estimate ci_width volatility
       0.50  12  0.03786   0.3772         NA
       0.55  16  0.07446   0.4514   0.047147
       0.60  21  0.06503   0.3639   0.044701
       0.65  26  0.09133   0.4236   0.033860
       0.70  34  0.09453   0.3661   0.028761
       0.75  44  0.11820   0.3951   0.014598
       0.80  56  0.14273   0.3835   0.008349
       0.85  72  0.16418   0.3789   0.016426
       0.90  93  0.18828   0.4094   0.029207
       0.95 120  0.19518   0.3510   0.035513
       1.00 154  0.19731   0.4152         NA
```

```r
m_selection$m_pow_of_n_optimal
#> [1] 0.8
```

Here the MV method lands on `m_pow_of_n = 0.80` (`volatility` is minimized at `m = 56`) -- close to,
and consistent with, the package's `0.75` default. `select_optimal_m_prop` is a calibration step: it
uses a small `B` per grid point since its total cost is `length(m_pow_of_n_grid) * B` replicates, so
rerun `PTE_bootstrap_inference` at the selected `m_pow_of_n_optimal` with a production-sized `B` (e.g.
the default `3000`) to get the confidence intervals and p-values you actually report:

```r
pte_results = PTE_bootstrap_inference(
    X, y,
    regression_type = "continuous",
    m_pow_of_n = m_selection$m_pow_of_n_optimal,
    B = 3000
)
pte_results
```

## Developer's note

On a machine whose R installation bakes `-march=native` into its site-wide `Makeconf` (see "How much
does the flag actually buy you?" above), `R CMD check --as-cran` will report a "non-portable
compilation flags used" NOTE even though the package's own build config never requests it. For future
runs, reuse this pattern to keep `-march=native` out of the check:

```
printf 'CFLAGS = -O2\nCXXFLAGS = -O2\n' > /tmp/Makevars.portable
R_MAKEVARS_USER=/tmp/Makevars.portable R CMD check --as-cran PTE_2.0.tar.gz
```
