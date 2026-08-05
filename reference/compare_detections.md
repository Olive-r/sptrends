# Compare detection methods against a known ground truth

Evaluates how well one or more trend-detection methods recover a known
truth – the scientific question this function exists to answer is
simple: *which method actually recovers the true trends better?* A
confusion-matrix comparison is how it answers that question, not the
point of the function itself. Built for
[`sim_trend_stack()`](https://olive-r.github.io/sptrends/reference/sim_trend_stack.md)'s
`true_slope`, but deliberately agnostic to where either side comes from.
`detections` and `ground_truth` are just logical vectors (or objects
that reduce to one): this works equally well comparing
[`workflow_tst()`](https://olive-r.github.io/sptrends/reference/workflow_tst.md)
against a classic Mann-Kendall, a linear model, a method from another
package, or an algorithm you implemented yourself – as long as each is
reduced to a "significant / not significant" call per cell.

## Usage

``` r
compare_detections(
  detections,
  ground_truth,
  metrics = c("sensitivity", "specificity", "precision", "accuracy", "f1", "mcc", "fpr",
    "fdr", "fwer", "type_i", "type_ii", "type_iii", "field_power", "global_power",
    "within_image_power", "directional_power"),
  replicates = FALSE,
  directions = NULL,
  truth_direction = NULL,
  evaluation_mask = NULL,
  verbose = TRUE
)
```

## Arguments

- detections:

  When `replicates = FALSE` (default): a named list, one element per
  method being compared. Each element is the method's significance call:
  a logical vector (`TRUE` = called significant), a single-layer
  [`terra::SpatRaster`](https://rspatial.github.io/terra/reference/SpatRaster-class.html)
  of 0/1 or `TRUE`/`FALSE` (coerced with `!= 0`), or a numeric/integer
  vector coerced the same way. The list names become the `Method`
  column. When `replicates = TRUE`: a list of such named lists, one per
  replicate (e.g. one per random seed) – every inner list must use the
  same method names.

- ground_truth:

  When `replicates = FALSE`: the ground truth, in any of the same forms
  as one element of `detections` – most commonly
  [`sim_trend_stack()`](https://olive-r.github.io/sptrends/reference/sim_trend_stack.md)'s
  `true_slope` (a `SpatRaster`, non-zero where a true trend exists) or
  its raw values. Coerced to logical the same way as `detections`. When
  `replicates = TRUE`: either a single ground truth reused for every
  replicate (the common case, when the same underlying truth is being
  detected repeatedly under different noise draws), or a list the same
  length as `detections` with one ground truth per replicate (for the
  rarer case where the truth itself is re-simulated every time too, e.g.
  [`sim_trend_stack()`](https://olive-r.github.io/sptrends/reference/sim_trend_stack.md)'s
  own `true_slope` with a different `seed` on each call).

- metrics:

  Character vector of derived metrics to include, from
  `c("sensitivity", "specificity", "precision", "accuracy", "f1", "mcc", "fpr", "fdr", "fwer", "type_i", "type_ii", "type_iii", "field_power", "global_power", "within_image_power", "directional_power")`
  (case-insensitive). The omitted default remains the original eight
  confusion-matrix metrics. The confusion matrix counts (`TP`, `FP`,
  `TN`, `FN`) are always included regardless of `metrics` – every
  requested metric is computed from that same confusion matrix, not
  recomputed independently. `"fwer"` is structurally different from the
  other eight – see its own entry in "Value" below – and can only be
  requested when `replicates = TRUE`; it errors otherwise, since a
  single run's own false-positive count is either zero or not, with no
  per-replicate ratio to compute in the first place. `"type_i"` is the
  background false-positive proportion; `"type_ii"` is the missed-signal
  proportion; `"type_iii"` is the wrong-direction proportion among
  detected true signals. `"field_power"` records whether a replicate
  produces at least one rejection anywhere in the evaluation domain. It
  can equal one because of a false positive alone. `"global_power"` is
  stricter: it records whether a replicate detects at least one true
  signal, whereas `"within_image_power"` is the proportion of signal
  cells detected.

- replicates:

  Logical. `FALSE` (default): score a single run, as described above.
  `TRUE`: score every replicate the same way, then aggregate to the mean
  and standard deviation of every column, per method – a single seed can
  favour or disfavour a method by chance, so this is what actually
  supports a claim like "CMK detects more true trends" (see `@examples`
  below). Generating the data and computing each method's detections for
  every replicate still stays as an ordinary loop you write, since that
  step depends entirely on which methods you are comparing and cannot be
  generalised without either guessing at your workflow or adding a
  callback-style interface – `replicates = TRUE` only does the scoring
  and aggregation, for the per-replicate detections/ground-truth you
  already computed.

- directions:

  Optional named list of estimated direction vectors, one per method.
  Alternatively, each detection may be a list containing `significant`
  and `direction`.

- truth_direction:

  Optional known direction vector, normally
  [`sim_trend_stack()`](https://olive-r.github.io/sptrends/reference/sim_trend_stack.md)'s
  `true_direction`. Required for Type III error and directional power.

- evaluation_mask:

  Optional logical vector/raster selecting one common evaluation domain
  for every method. Missing and `FALSE` cells are not scored.
  Method-specific masks are rejected because methods must be compared
  against the same cells and the same ground truth.

- verbose:

  Logical. If `TRUE`, reports progress, elapsed time and the estimated
  time remaining while runs are scored.

## Value

**When `replicates = FALSE`**: a data frame with one row per method:
`Method`, `TP`, `FP`, `TN`, `FN`, and the requested metric columns among
`Sensitivity` (a.k.a. recall or power: `TP / (TP + FN)`), `Specificity`
(`TN / (TN + FP)`), `Precision` (`TP / (TP + FP)`), `Accuracy`
(`(TP + TN) / (TP + FP + TN + FN)` – included because it is widely
expected, though it can be misleading whenever true trends are rare
relative to the whole raster, which is common in these simulations),
`F1` (harmonic mean of precision and sensitivity), `MCC` (Matthews
correlation coefficient – a single balanced summary of all four
confusion-matrix cells at once, generally preferred over `Accuracy` or
`F1` under the class imbalance a low `trend_fraction` produces), and
`FPR` (false positive rate, `FP / (FP + TN)`). `FDR`, in this single-run
table, is the realised false discovery *proportion* (Benjamini &
Hochberg, 1995) in this one result (`FP / (FP + TP)`) – not the false
discovery *rate* itself, which is formally defined as the expectation of
that proportion over many repetitions. `q` in
[`fdr_correction()`](https://olive-r.github.io/sptrends/reference/fdr_correction.md)
bounds that expectation, not the value in any single realisation – see
`replicates = TRUE` below for what actually estimates the expectation
this proportion is a single draw from. By standard convention, this
realised proportion is zero when no discoveries are made; retaining
those zeroes is essential when estimating FDR across replicates. Other
ratios with an undefined denominator (e.g. `Precision` when a method
calls nothing significant, or `MCC` when any confusion-matrix margin is
zero) are returned as `NA`, not `NaN` or `0`. `FieldPower`,
`GlobalPower`, and `WithinImagePower` distinguish, respectively, whether
any cell was rejected, whether at least one true signal cell was
detected, and the fraction of true signal cells detected.

**When `replicates = TRUE`**: a data frame with one row per method:
`Method`, `n_replicates`, and `<column>_mean`/`<column>_sd` for every
column above (`TP_mean`, `TP_sd`, `Sensitivity_mean`, `Sensitivity_sd`,
and so on). `FDR_mean` here is the actual estimate of the false
discovery rate itself (the average of the per-replicate proportion
described above, across replicates) – the quantity `q` in
[`fdr_correction()`](https://olive-r.github.io/sptrends/reference/fdr_correction.md)
is meant to bound, unlike the single-run `FDR` column above. If `"fwer"`
was requested, also a plain `FWER` column (no `_mean`/`_sd` pair): the
proportion of replicates with at least one false positive (`FP > 0`) –
the family-wise error rate (Benjamini & Hochberg, 1995's own point of
contrast for FDR), a single already-aggregated rate rather than a
per-replicate value with its own mean and spread.

## Details

**Function type:** **Validation function** – scores detection methods
against known truth. It is not part of the inferential pipeline.

## Typical use

Single-run validation:

    sim_trend_stack()
        |
    run one or more detection methods on sim$series
        |
    compare_detections()
        |
    one table of cell-wise detection metrics

Replicated validation is an alternative route, not a step after the
single-run call:

    repeat simulation and detection across seeds
        |
    collect detections and ground truths in lists
        |
    compare_detections(replicates = TRUE)
        |
    mean and standard deviation of metrics + optional empirical FWER

See the examples below for both routes worked through in full.

## Methodological details

**Why use simulated data?**

Real environmental datasets have no known "correct" answer, so detection
performance cannot be measured directly – there is no ground truth to
score against for a real NDVI or temperature raster. Simulation studies
([`sim_trend_stack()`](https://olive-r.github.io/sptrends/reference/sim_trend_stack.md))
provide a known ground truth instead, allowing sensitivity, specificity,
precision, and the other metrics below to be computed objectively,
something no amount of visual inspection of a real map can substitute
for.

**A platform for comparing methods, not just running them**

Together,
[`sim_trend_stack()`](https://olive-r.github.io/sptrends/reference/sim_trend_stack.md),
this function, and the two published workflows this package offers
([`workflow_tst()`](https://olive-r.github.io/sptrends/reference/workflow_tst.md)
and
[`workflow_rta()`](https://olive-r.github.io/sptrends/reference/workflow_rta.md))
mean sptrends is not only a tool for running a trend analysis, but also
a platform for validating and comparing *how well* different methods do
it – including methods this package did not itself implement, since both
`detections` and `ground_truth` are accepted in plain, method-agnostic
forms (see the `detections` argument below).

**Methods and metric selection**

All requested metrics derive from the same cell-wise confusion matrix.
Sensitivity and specificity describe conditional detection behaviour;
precision and FDR describe the rejection set; MCC provides a balanced
summary under class imbalance. Empirical FDR and FWER require repeated
simulations and therefore are only interpretable in replicated mode.

**Statistical assumptions and limitations**

Detection and truth must refer to the same cells in the same order and
reduce to binary calls. Results describe performance under the supplied
simulated scenarios; they do not establish performance for every real
environmental process. A single replicate estimates realised
proportions, not repeated-sampling error rates.

**Computational considerations**

Scoring is vectorised over the common evaluation domain and is generally
inexpensive relative to generating data or fitting the compared methods.
Replicated mode retains method-level summaries rather than raster
outputs.

**Quality assurance**

Tests compare confusion counts and every derived metric with
analytically tractable cases, including undefined denominators.
Additional tests cover vector/raster inputs, replicated aggregation,
empirical FDR and FWER, validation failures, method ordering and S3
printing, summaries and plots. In the retained external validation,
`sptrends` and `Kendall` produced identical MK decisions, directions and
derived performance metrics across 4,000 paired Monte Carlo replicates.
See `inst/validation/` and
[`?sptrends`](https://olive-r.github.io/sptrends/reference/sptrends-package.md)
for scope and limitations.

## References

Fawcett, T. (2006) An introduction to ROC analysis. Pattern Recognition
Letters, 27(8), 861-874.
[doi:10.1016/j.patrec.2005.10.010](https://doi.org/10.1016/j.patrec.2005.10.010)

Source of the Matthews correlation coefficient (`MCC`):

- Matthews, B.W. (1975) Comparison of the predicted and observed
  secondary structure of T4 phage lysozyme. Biochimica et Biophysica
  Acta (BBA) - Protein Structure, 405(2), 442-451.
  [doi:10.1016/0005-2795(75)90109-9](https://doi.org/10.1016/0005-2795%2875%2990109-9)

Source of the false discovery rate concept `FDR`/`FDR_mean` estimate,
and of the family-wise error rate `FWER` is contrasted against (see
"Value" above for both):

- Benjamini, Y., & Hochberg, Y. (1995) Controlling the False Discovery
  Rate: A Practical and Powerful Approach to Multiple Testing. Journal
  of the Royal Statistical Society: Series B, 57, 289-300.
  [doi:10.1111/j.2517-6161.1995.tb02031.x](https://doi.org/10.1111/j.2517-6161.1995.tb02031.x)

## See also

Other validation functions:
[`benchmark_methods()`](https://olive-r.github.io/sptrends/reference/benchmark_methods.md),
[`benchmark_summary()`](https://olive-r.github.io/sptrends/reference/benchmark_summary.md),
[`plot_detection_comparison()`](https://olive-r.github.io/sptrends/reference/plot_detection_comparison.md),
[`simulation_design()`](https://olive-r.github.io/sptrends/reference/simulation_design.md)

## Examples

``` r
# Simulated data only, deliberately -- compare_detections() needs a
# KNOWN ground truth to score against (sim$true_slope below), which
# by definition only sim_trend_stack() can provide; real environmental
# data has no equivalent "true" answer to compare a detection against.
# \donttest{
sim <- sim_trend_stack(nrow = 15, ncol = 15, n_time = 15, seed = 1)
#> >> [sim_trend_stack()] elapsed: 0.19 s

# Two variants of the same test: classic Mann-Kendall (no
# neighbourhood averaging) vs. the Contextual version.
trend_mk  <- trend_test(sim$series, method = "MK",
                         report = FALSE, verbose = FALSE)
trend_cmk <- trend_test(sim$series, method = "CMK",
                         report = FALSE, verbose = FALSE)

# Since sim$true_slope is known exactly (this is simulated data), we
# can score each method's raw significance calls against the truth.
compare_detections(
  detections = list(MK = trend_mk$stats$p <= 0.05,
                     CMK = trend_cmk$stats$p <= 0.05),
  ground_truth = sim$true_slope, verbose = FALSE
)
#>   Method TP FP TN  FN Sensitivity Specificity Precision  Accuracy        F1
#> 1     MK 53  1 30 141   0.2731959   0.9677419 0.9814815 0.3688889 0.4274194
#> 2    CMK 67  1 30 127   0.3453608   0.9677419 0.9852941 0.4311111 0.5114504
#>         MCC        FPR        FDR
#> 1 0.1944427 0.03225806 0.01851852
#> 2 0.2349981 0.03225806 0.01470588

# replicates = TRUE: repeat the comparison across several random
# seeds and aggregate -- a single run can favour a method just by
# chance. Both detections and the true_slope ground truth change
# every replicate here, so ground_truth is a list too, not a single
# shared vector.
detections_list <- list()
truths_list <- list()
for (s in 1:10) {
  sim_s <- sim_trend_stack(nrow = 12, ncol = 12, n_time = 12, seed = s)
  mk_s  <- trend_test(sim_s$series, method = "MK",
                       report = FALSE, verbose = FALSE)
  cmk_s <- trend_test(sim_s$series, method = "CMK",
                       report = FALSE, verbose = FALSE)
  detections_list[[s]] <- list(MK = mk_s$stats$p <= 0.05,
                                CMK = cmk_s$stats$p <= 0.05)
  truths_list[[s]] <- sim_s$true_slope
}
#> >> [sim_trend_stack()] elapsed: 0.10 s
#> >> [sim_trend_stack()] elapsed: 0.04 s
#> >> [sim_trend_stack()] elapsed: 0.04 s
#> >> [sim_trend_stack()] elapsed: 0.05 s
#> >> [sim_trend_stack()] elapsed: 0.04 s
#> >> [sim_trend_stack()] elapsed: 0.04 s
#> >> [sim_trend_stack()] elapsed: 0.05 s
#> >> [sim_trend_stack()] elapsed: 0.04 s
#> >> [sim_trend_stack()] elapsed: 0.04 s
#> >> [sim_trend_stack()] elapsed: 0.04 s
compare_detections(detections_list, truths_list, replicates = TRUE,
                   verbose = FALSE)
#>   Method n_replicates TP_mean     TP_sd FP_mean    FP_sd TN_mean    TN_sd
#> 1     MK           10    22.2  7.656515     1.4 1.173788    17.4 2.011080
#> 2    CMK           10    22.4 11.862171     3.0 2.160247    15.8 3.047768
#>   FN_mean     FN_sd Sensitivity_mean Sensitivity_sd Specificity_mean
#> 1   103.0  7.874008        0.1773437     0.06155457          0.92625
#> 2   102.8 12.389960        0.1792087     0.09600344          0.83750
#>   Specificity_sd Precision_mean Precision_sd Accuracy_mean Accuracy_sd
#> 1     0.06022239      0.9396100   0.05615718     0.2750000  0.05260286
#> 2     0.11516896      0.8761194   0.06431335     0.2652778  0.07841723
#>     F1_mean      F1_sd   MCC_mean     MCC_sd FPR_mean     FPR_sd  FDR_mean
#> 1 0.2944415 0.08273766 0.09303124 0.06585726  0.07375 0.06022239 0.0603900
#> 2 0.2884261 0.12094978 0.01317895 0.09155006  0.16250 0.11516896 0.1238806
#>       FDR_sd
#> 1 0.05615718
#> 2 0.06431335

# FWER ("did at least one false positive happen in this run at all?")
# only means anything across many replicates -- request it alongside
# the usual metrics with replicates = TRUE; it errors with
# replicates = FALSE, since a single run cannot estimate a rate.
compare_detections(detections_list, truths_list, replicates = TRUE,
                    metrics = c("sensitivity", "fdr", "fwer"),
                    verbose = FALSE)
#>   Method n_replicates TP_mean     TP_sd FP_mean    FP_sd TN_mean    TN_sd
#> 1     MK           10    22.2  7.656515     1.4 1.173788    17.4 2.011080
#> 2    CMK           10    22.4 11.862171     3.0 2.160247    15.8 3.047768
#>   FN_mean     FN_sd Sensitivity_mean Sensitivity_sd  FDR_mean     FDR_sd FWER
#> 1   103.0  7.874008        0.1773437     0.06155457 0.0603900 0.05615718  0.7
#> 2   102.8 12.389960        0.1792087     0.09600344 0.1238806 0.06431335  1.0
# }
```
