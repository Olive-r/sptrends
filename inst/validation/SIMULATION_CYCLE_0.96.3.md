# External validation of the simulation and benchmarking cycle

## Scope

This record summarises the independent full run completed for 0.96.3. The
validation script exercised the simulator, factorial design, known-truth
scoring, benchmark orchestration, aggregation and scenario-performance
graphics. All 33 prespecified controls passed.

The run used 1,000 independently generated fields for each spatial covariance
model and 500 paired Monte Carlo replicates in each of eight combinations of
spatial model, temporal AR(1) and trend presence. Paired methods received the
same generated realisation and evaluation domain.

## Spatial and temporal controls

| Process | Target | Observed | Absolute error |
|---|---:|---:|---:|
| Independent unit-lag correlation | 0 | -0.000736 | 0.000736 |
| Exponential unit-lag correlation | 0.65 | 0.648648 | 0.001352 |
| Gaussian unit-lag correlation | 0.65 | 0.652780 | 0.002780 |
| Matérn unit-lag correlation | 0.835004 | 0.829781 | 0.005223 |
| Temporal AR(1) | 0.70 | 0.700570 | 0.000570 |
| Stationary marginal SD | 2.00 | 2.003869 | 0.003869 |

The analytical Matérn correlation also agreed with `fields::Matern()` within
`1e-10`. Gaussian and Student-t marginal distributions, heavier Student-t
tails, exact signal masks, slopes, directions, break truth, reproducibility
and retained generating specifications passed their recorded tolerances.

## Independent Mann-Kendall comparison

The internal MK wrapper and `Kendall::MannKendall()` produced identical
cell-level decisions and directions in all retained paired runs. Consequently,
all confusion-matrix, Type I, Type II, Type III, power, FDR and FWER metrics
were identical. This validates decision extraction and benchmark scoring for
this comparison; it is not a claim that every implementation returns identical
continuous statistics or uses identical conventions.

Without temporal autocorrelation, the empirical cell-wise Type I error was
0.045-0.048. With AR(1) = 0.5 it increased to 0.209-0.218, an expected
demonstration that uncorrected MK is not protected against serial dependence.
Under the simulated signal, within-image power was 0.780-0.784 without AR(1)
and 0.724-0.744 with AR(1) = 0.5.

## Graphical correction

The initial validation figure placed the null scenario (`trend_strength = 0`)
on a power axis, where power is undefined, and pooled temporal AR(1) conditions.
The corrected presentation therefore uses two separate figures:

1. Type I error against AR(1), using only complete-null scenarios.
2. Within-image power against AR(1), using only signal scenarios.

Both are stratified by spatial model and retain both independent MK
implementations. Their exact overlap is evidence of agreement rather than a
missing series.

## Reproducibility and limitations

The original external output retained checks, replicate-level results,
summaries, spatial-correlation estimates and session information. Runtime is
recorded descriptively but is not treated as a universal speed benchmark.
These controls validate the tested data-generating and scoring contracts under
the stated scenarios; they do not establish realism for every environmental
process or validate future methods automatically.
