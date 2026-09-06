# Global and local permutation-based spatial autocorrelation inference.

#' @noRd
.prepare_moran_neighbourhood <- function(x, ok, connectivity = "queen") {
  nb <- .prepare_spatial_neighbourhood(x, ok, connectivity)
  list(W = nb$W, nb_count = nb$nb_count, S0 = nb$S0,
       signature = nb$signature)
}

#' @noRd
.spatial_local_statistic <- function(vals, W, ok, method, W_star = NULL) {
  z <- vals
  z[!ok] <- 0
  if (method == "moran") {
    z[ok] <- z[ok] - mean(z[ok])
    m2 <- sum(z^2) / sum(ok)
    return(z * as.numeric(W %*% z) / m2)
  }

  # Local Getis-Ord Gi*: the focal cell is included in W* and the
  # statistic is the local weighted sum divided by the field-wide sum.
  # Its permutation distribution, rather than an analytic normal
  # approximation, supplies the inferential reference distribution.
  if (is.null(W_star)) {
    W_star <- W + Matrix::Diagonal(x = as.numeric(ok))
  }
  as.numeric(W_star %*% z) / sum(z)
}

#' @noRd
.spatial_local_chunk <- function(indices, seeds, x_vals, valid_idx,
                                 valid_vals, W, W_star, ok, method,
                                 has_neighbour, observed, phase,
                                 alternative, null_mean = NULL) {
  N <- length(valid_idx)
  total <- total_sq <- greater <- less <- numeric(N)
  extreme <- numeric(N)

  for (position in seq_along(indices)) {
    i <- indices[position]
    set.seed(seeds[i])
    permuted <- x_vals
    permuted[valid_idx] <- valid_vals[sample.int(N)]
    statistic <- .spatial_local_statistic(
      permuted, W, ok, method, W_star
    )[valid_idx]
    statistic[!has_neighbour] <- NA_real_

    if (phase == 1L) {
      values <- statistic
      values[!is.finite(values)] <- 0
      total <- total + values
      total_sq <- total_sq + values^2
      greater <- greater + ((statistic >= observed) %in% TRUE)
      less <- less + ((statistic <= observed) %in% TRUE)
    } else {
      if (alternative == "two.sided") {
        is_extreme <- abs(statistic - null_mean) >=
          abs(observed - null_mean)
        extreme <- extreme + (is_extreme %in% TRUE)
      }
    }
  }
  list(total = total, total_sq = total_sq, greater = greater,
       less = less, extreme = extreme)
}

#' @noRd
.spatial_raster_from_values <- function(template, ok, values, name) {
  out <- terra::rast(template)
  full <- rep(NA_real_, terra::ncell(template))
  full[ok] <- values
  terra::values(out) <- full
  names(out) <- name
  out
}

#' Permutation-based spatial autocorrelation tests
#'
#' Provides one interface for global and local spatial-autocorrelation
#' analyses of a single-layer raster. Global analysis computes Moran's I
#' or Getis-Ord General G. Local analysis computes a local Moran's I or
#' Getis-Ord Gi* statistic for every valid cell. Local significance is
#' assessed by spatial permutation. Its p-value raster can subsequently
#' be passed to [fdr_correction()] for BH, BKY or BY control.
#'
#' This is a general-purpose spatial diagnostic, not an FDR-specific
#' utility. It can be applied to environmental variables, model
#' residuals, estimated coefficients, test statistics, probabilities
#' or other numeric spatial fields.
#'
#' **Function type:** **Spatial diagnostic function** -- performs global
#' or local spatial-autocorrelation inference on one spatial field. It
#' is independent of the temporal trend workflows.
#'
#' @section Typical use:
#' ```
#' one spatial raster (variable, residual, coefficient, or statistic)
#'     |
#' spatial_autocorrelation(scope = "global" or "local")
#'     |
#' observed statistic + permutation inference + null distribution
#' ```
#' The input is one spatial field, not a raster time series. A trend
#' workflow output can be analysed only after selecting one derived
#' layer, such as a slope or test-statistic raster.
#'
#' @param x A single-layer `terra::SpatRaster` containing finite numeric
#'   values or `NA`.
#' @param method Which statistic to compute. `"moran"` (default): global
#'   Moran's I (Moran, 1950) -- accepts any numeric values (positive,
#'   negative, or mixed), tests whether nearby cells tend to have
#'   *similar* values, whatever those values are. `"getis_ord"`: Getis-Ord
#'   General G (Getis & Ord, 1992) -- **requires non-negative values**
#'   (errors otherwise; see "Methodological background" below for why),
#'   tests specifically whether *high* values cluster near other high
#'   values (or low near low) rather than just "similar near similar" --
#'   a genuinely different question from Moran's I, not a variant of the
#'   same one, and the two do not always agree. Choose Moran's I when the
#'   question is whether neighbouring values are generally similar;
#'   choose General G when the question specifically concerns clustering
#'   of high values.
#' @param scope Spatial scope. `"global"` (default) computes one statistic
#'   for the complete raster. `"local"` computes one local Moran's I or
#'   Getis-Ord Gi* statistic and one permutation p-value per valid cell.
#' @param connectivity `"queen"` (default) or `"rook"`.
#' @param nperm Positive integer. Number of permutations. Start with `99`
#'   to check
#'   the analysis runs correctly, use `999+` for a result to report, and
#'   `9999` for publication-quality analyses when computationally
#'   feasible (p-value resolution is `1/(nperm + 1)`).
#' @param alternative `NULL` (default) selects `"greater"` for global
#'   tests and `"two.sided"` for local tests. May also be supplied
#'   explicitly as `"greater"`, `"two.sided"`, or `"less"`. The global
#'   default reflects the common question of whether nearby observations
#'   are more positively associated than expected under spatial
#'   randomisation.
#' @param alpha One finite number strictly between 0 and 1. Significance
#'   level used only for the exploratory, unadjusted
#'   `significant_raw` map returned by local analysis. It does not
#'   control multiplicity. Use [fdr_correction()] on the returned `p`
#'   raster and specify its target `q` for BH, BKY or BY inference.
#' @param seed One finite numeric value or `NULL`. Random seed for
#'   reproducibility (used for both sequential and parallel execution).
#' @param n_cores Positive integer. Number of cores for the permutation
#'   loop (each permutation is independent). `1` (default): sequential.
#'   `> 1`: uses a
#'   `parallel::makeCluster()` PSOCK cluster. The spatial topology (the
#'   expensive part for large rasters) is computed once and shared, not
#'   recomputed per permutation or per worker.
#' @param precomputed_neighbourhood Optional. Skips recomputing the
#'   spatial adjacency (the same "expensive for large rasters" step
#'   `n_cores` above already shares across permutations) when this
#'   function is called repeatedly on the same raster geometry -- most
#'   usefully, together with [trend_test()]'s own
#'   `precomputed_neighbourhood` argument, when both functions are run on
#'   the same geometry in the same analysis. Accepts the output of
#'   [prepare_cmk_neighbourhood()] directly when it uses the matching
#'   default 3 by 3 adjacency -- despite the name, both functions then
#'   build the identical structure underneath. Broader CMK windows have a
#'   different signature and are rejected rather than silently substituted.
#'   Its `W` component is
#'   validated as a square matrix with one row and column per raster cell.
#'   If `NULL` (typical single-call use), it is computed internally.
#' @param report Logical. If `TRUE` (default), automatically print the
#'   scope-specific summary and draw its plot once the test finishes.
#'   Global output shows the null distribution; local output maps the
#'   statistic, empirical z, raw p-value and exploratory unadjusted
#'   significance decision.
#' @param verbose Logical. Print progress messages and elapsed time.
#'
#' @section Methodological details:
#' **Applications**
#'
#' Typical applications include describing clustering or dispersion in
#' environmental variables, detecting spatial structure in model
#' residuals, and examining the spatial distribution of estimated
#' trends or effect sizes. One optional application is diagnosing
#' dependence in a field of cell-wise inferential results before
#' selecting or interpreting a multiple-testing procedure. Moran's I
#' does not formally test independence or the complete dependence
#' conditions underlying BH, BKY or BY, and must not be used as an
#' automatic selector of an FDR procedure.
#'
#' Spatial autocorrelation is not itself a temporal trend test. The input
#' represents one spatial distribution, not a time series. The optional
#' `moran_check` argument of [fdr_correction()] is only one convenience
#' use of this independent diagnostic and has the same interpretive
#' limitation described above.
#'
#' **Local statistics**
#'
#' With `scope = "local"` and `method = "moran"`, the statistic for cell
#' `i` is `I_i = z_i sum_j(w_ij z_j) / m_2`, where `z_i` is the deviation
#' from the global mean and `m_2 = sum_i(z_i^2) / N`. The binary queen or
#' rook weights are not row-standardised. Under this definition,
#' `sum_i(I_i) = S0 * I_global`, which is checked directly in the test
#' suite. Positive values indicate locally similar deviations; negative
#' values indicate a focal deviation surrounded by deviations of the
#' opposite sign. Significance still depends on the permutation test.
#'
#' With `method = "getis_ord"`, the local Gi* statistic is
#' `G_i* = sum_j(w_ij* x_j) / sum_j(x_j)`, where `w_ii* = 1`: the focal
#' cell is included, which distinguishes Gi* from Gi. High values indicate
#' local concentrations of high raw values. The same non-negative-input
#' requirement as General G applies. This function uses each cell's
#' permutation distribution rather than an analytic normal approximation.
#' The returned `z` raster standardises each observed local statistic by
#' its own permutation mean and standard deviation. It aids comparison
#' across cells with different neighbour counts, but it is not assigned
#' an analytic normal p-value; use the returned permutation inference.
#' Valid cells with no valid spatial neighbour are returned as `NA` for
#' local inference and are excluded from the multiple-testing family.
#'
#' **Methods and method selection**
#'
#' - **Original publications**: Moran (1950) for `method = "moran"`;
#'   Getis & Ord (1992) for `method = "getis_ord"`.
#' - **Main references**: Moran (1950)/Getis & Ord (1992) for the
#'   statistics themselves; Tobler (1970) for the conceptual grounding of
#'   why spatial dependence is the expected default in geographic data;
#'   Hope (1968) for the permutation-based significance test used here in
#'   place of either statistic's own analytic Z approximation. Full
#'   citations under "References" below.
#' - **Why Getis-Ord needs non-negative values, specifically**: General G
#'   is a ratio of a weighted sum of *products* of raw values
#'   (`sum(w_ij * x_i * x_j)`) to the sum of all pairwise products
#'   (`sum(x_i * x_j)`), for `i != j` -- unlike Moran's I, it does not
#'   mean-centre the values first. With mixed-sign data, `x_i * x_j` can
#'   be negative for two high-magnitude values of opposite sign, which
#'   breaks the statistic's actual interpretation (a proportion of
#'   "high-value pairs found near each other") -- not merely unusual
#'   input, a violated precondition. Matches the same restriction in
#'   ArcGIS's High/Low Clustering (Getis-Ord General G) tool and other
#'   established implementations.
#' - **Typical applications**: Moran's I for "are similar values near
#'   each other"; General G for "do high values specifically cluster
#'   together". Inputs may be raw environmental variables, residuals,
#'   coefficients or inferential fields. For mixed-sign statistics,
#'   Moran is usually meaningful directly; General G requires a
#'   scientifically justified non-negative input.
#'
#' **Statistical assumptions and permutation inference**
#'
#' The classic analytic formulas for both statistics' significance need
#' either a normality assumption (the values are normally distributed) or
#' a randomisation assumption (the values are exchangeable), and the
#' resulting Z-score can be inaccurate whenever the actual data violate
#' them (Hope, 1968). Environmental variables, residuals and inferential
#' outputs can all be bounded, skewed or otherwise non-normal. A
#' permutation test builds the null distribution by actually
#' reshuffling the observed values across the raster's cells (preserving
#' their real distribution, skew, and bounds, whatever it is) and
#' recomputing the chosen statistic each time, so significance comes from
#' resampling the data you actually have rather than from an assumption
#' about a distribution it may not follow.
#'
#' **Local inference and multiple testing**
#'
#' A local statistic map without a null distribution is descriptive: it
#' shows where local association is large or small, but it supplies no
#' inferential *p*-value. Consequently, there is no meaningful
#' "alpha without permutations" in this implementation. Local p-values
#' are Monte Carlo p-values obtained by permuting the observed values.
#'
#' The returned `significant_raw` raster applies `p_i <= alpha`, treating
#' every cell as though it were the only hypothesis. It is exploratory
#' and controls neither FDR nor FWER across the map. For multiple-testing
#' inference, pass the returned `p` raster to [fdr_correction()]. That
#' function provides BH, BKY and BY from one independently tested and
#' documented implementation. BH and BKY require independent tests or
#' suitable positive dependence; BY is valid under arbitrary dependence
#' but is usually more conservative. In all three cases, *q* is the target
#' FDR and is not another name for the unadjusted per-cell *alpha*.
#'
#' Every permutation is global and joint: all valid values are reassigned
#' across the raster once, and every local statistic is recomputed from
#' that same reassignment. It is not the conditional LISA permutation
#' used by some software, where the focal value is held fixed and potential
#' neighbours are sampled (Anselin, 1995). The two randomisation
#' hypotheses answer different questions and their cell-wise p-values
#' need not coincide. Permutation validity requires exchangeability under
#' the stated spatial-randomisation null; it is not assumption-free.
#'
#' More permutations improve Monte Carlo resolution and stability; they
#' do not mechanically make a result more significant. The smallest
#' attainable *p*-value is `1 / (nperm + 1)` because both numerator and
#' denominator use the standard plus-one correction (Phipson & Smyth,
#' 2010), so a randomly sampled permutation test never reports zero. A
#' gradual strategy is useful:
#' use about 99 permutations for exploratory checking, 999 or more for a
#' stable analysis, and 9999 when tail resolution is important and the
#' computational cost is acceptable. Results based on only 99
#' permutations should be described as exploratory rather than final.
#'
#' **Computational considerations**
#'
#' Runtime grows with the number of valid cells and permutations; local
#' analysis must update a statistic for every tested cell in every joint
#' randomisation. Reusing `precomputed_neighbourhood` avoids rebuilding
#' topology, while `n_cores > 1` distributes independent permutations.
#' Increasing `nperm` improves p-value resolution but increases runtime
#' approximately proportionally.
#'
#' **Limitations**
#'
#' **Methodological**: the magnitude of either statistic is not
#' universally comparable across different weight matrices (Cliff & Ord,
#' 1981), so this function does not categorise the result into
#' "low/moderate/strong" by default (see [classify_moran()] for an
#' explicitly non-standard convenience label, `method = "moran"` only --
#' General G's natural range and expected value under H0 differ enough
#' from Moran's I that the same low/moderate/strong thresholds would not
#' transfer meaningfully, so no equivalent label is offered for `method =
#' "getis_ord"`). With `scope = "global"`, the result is one number for
#' the complete raster and does not identify where clusters occur.
#' Local Moran's I and Getis-Ord Gi* require one test per cell and,
#' consequently, explicit multiple-testing control for confirmatory maps.
#' Apply BH, BKY or BY with [fdr_correction()] and its target `q`; do not
#' interpret the raw per-cell `alpha` map as family-level inference.
#' This diagnostic is evidence of spatial dependence, not proof
#' that a downstream method's complete dependence assumptions hold.
#'
#' **Implementation**: fixed queen/rook neighbourhood, not a generic
#' distance-band or k-nearest definition. NA cells are excluded, not
#' imputed.
#'
#' **Quality assurance**
#'
#' Moran's I is checked against direct matrix calculations and exact
#' small-raster references. Permutation tests cover reproducible seeds,
#' alternative hypotheses, queen/rook neighbourhoods, missing and
#' constant rasters, reused adjacency objects, and sequential/parallel
#' equivalence. Local Moran's I is checked through its exact algebraic
#' identity with global Moran's I. Gi* is checked against hand-computed
#' small-raster results. Raw permutation p-values, raster geometry,
#' missing values, local S3 methods, and forwarding the local p-value
#' raster to [fdr_correction()] have separate regression tests. BH, BKY
#' and BY are independently tested in that function's own suite. See
#' `?sptrends` for the package-wide release-check protocol.
#'
#' @return For `scope = "global"`, returns an object of class
#'   `c("spatial_autocorrelation_global", "spatial_autocorrelation",
#'   "sptrends")`: a list with `statistic` (observed Moran's I or
#'   Getis-Ord G, depending on `method`), `method`, `sign` (`"positive"`
#'   or `"negative"`, `method = "moran"` only -- General G is always
#'   non-negative for non-negative inputs), `scope`, `p` (permutation
#'   p-value), `alternative`, `nperm`, `null_dist` (the null
#'   distribution), and `N` (valid cell count).
#'   For `scope = "local"`, returns an object of class
#'   `c("spatial_autocorrelation_local", "spatial_autocorrelation",
#'   "sptrends")`. Its `statistic`, `z`, `p`, and `significant_raw`
#'   components are single-layer rasters. `null_mean` and `null_sd`
#'   store cell-wise permutation summaries and the empirical standardised
#'   statistic without retaining the potentially enormous
#'   cell-by-permutation matrix. `permutation_seeds` records the generated
#'   streams for exact auditing of the Monte Carlo calculation. `N` is
#'   the valid-cell count and `N_tested` excludes valid cells without any
#'   valid neighbour; their local outputs are `NA` and
#'   [fdr_correction()] excludes them from the hypothesis family.
#'
#' @references
#' Primary references for the statistics being tested:
#' - Moran, P.A.P. (1950) Notes on Continuous Stochastic Phenomena.
#'   Biometrika, 37(1-2), 17-23. \doi{10.1093/biomet/37.1-2.17}
#' - Getis, A. and Ord, J.K. (1992) The Analysis of Spatial Association by
#'   Use of Distance Statistics. Geographical Analysis, 24(3), 189-206.
#'   \doi{10.1111/j.1538-4632.1992.tb00261.x}
#' - Anselin, L. (1995) Local Indicators of Spatial Association -- LISA.
#'   Geographical Analysis, 27(2), 93-115.
#'   \doi{10.1111/j.1538-4632.1995.tb00338.x}
#' - Ord, J.K. and Getis, A. (1995) Local Spatial Autocorrelation
#'   Statistics: Distributional Issues and an Application. Geographical
#'   Analysis, 27(4), 286-306.
#'   \doi{10.1111/j.1538-4632.1995.tb00912.x}
#'
#' Conceptual grounding for why spatial dependence is the expected default
#' in geographic data (motivates testing for it at all):
#' - Tobler, W. (1970) A Computer Movie Simulating Urban Growth in the
#'   Detroit Region. Economic Geography, 46(sup1), 234-240.
#'   \doi{10.2307/143141}
#'
#' Extended treatment of Moran's I and related spatial autocorrelation
#' statistics:
#' - Cliff, A.D. and Ord, J.K. (1981) Spatial Processes: Models and
#'   Applications. Pion, London.
#'
#' Basis for using permutation rather than the analytic Z approximation to
#' obtain the p-value (see "Why permutation" above):
#' - Hope, A.C. (1968) A simplified Monte Carlo significance test
#'   procedure. Journal of the Royal Statistical Society, Series B,
#'   30(3), 582-598. \doi{10.1111/j.2517-6161.1968.tb00759.x}
#' - Phipson, B. and Smyth, G.K. (2010) Permutation P-values Should Never
#'   Be Zero: Calculating Exact P-values When Permutations Are Randomly
#'   Drawn. Statistical Applications in Genetics and Molecular Biology,
#'   9(1), Article 39. \doi{10.2202/1544-6115.1585}
#'
#' On FDR procedures and the positive-dependence setting this function
#' can diagnose but cannot establish:
#' - Benjamini, Y., & Yekutieli, D. (2001) The control of the false
#'   discovery rate in multiple testing under dependency. Annals of
#'   Statistics, 29(4), 1165-1188. \doi{10.1214/aos/1013699998}
#' - Benjamini, Y. and Hochberg, Y. (1995) Controlling the False
#'   Discovery Rate: A Practical and Powerful Approach to Multiple
#'   Testing. Journal of the Royal Statistical Society: Series B, 57,
#'   289-300. \doi{10.1111/j.2517-6161.1995.tb02031.x}
#' - Benjamini, Y., Krieger, A.M. and Yekutieli, D. (2006) Adaptive
#'   Linear Step-Up Procedures that Control the False Discovery Rate.
#'   Biometrika, 93(3), 491-507. \doi{10.1093/biomet/93.3.491}
#' @examples
#' # Use a small real-data region at the dataset's native resolution so the
#' # example remains suitable for routine package checks.
#' series <- read_ordered_stack(example_data("vhp_ndvi"))
#' complete <- which(stats::complete.cases(
#'   terra::values(series, mat = TRUE)
#' ))
#' centre <- terra::xyFromCell(
#'   series, complete[ceiling(length(complete) / 2)]
#' )
#' resolution <- terra::res(series)
#' region <- terra::ext(c(
#'   centre[1] + c(-6, 6) * resolution[1],
#'   centre[2] + c(-6, 6) * resolution[2]
#' ))
#' series <- terra::crop(series, region, snap = "near")
#' X <- terra::values(series, mat = TRUE)
#' ok <- stats::complete.cases(X)
#' r <- series[[1]]
#' r_values <- terra::values(r, mat = FALSE)
#' r_values[!ok] <- NA_real_
#' terra::values(r) <- r_values
#' \donttest{
#' # I > 0 means neighbouring cells tend to have similar values
#' # (positive spatial autocorrelation); result$p is the permutation
#' # p-value for that observed I.
#' # Nineteen permutations are sufficient for this interface example;
#' # substantive inference should use a larger value.
#' result <- spatial_autocorrelation(
#'   r, nperm = 19, seed = 1, report = FALSE, verbose = FALSE
#' )
#' result
#' plot(result)   # the null distribution, with the observed I marked
#'
#' # Getis-Ord General G needs non-negative values -- abs() of a
#' # (possibly mixed-sign) trend statistic is a common way to get there
#' # when the question is "do high-magnitude values cluster together?"
#' r_abs <- abs(r)
#' result_g <- spatial_autocorrelation(r_abs, method = "getis_ord",
#'                                      nperm = 19, seed = 1,
#'                                      report = FALSE, verbose = FALSE)
#' result_g$statistic
#'
#' # One optional inferential application: describing spatial association
#' # in the p-value raster from trend_test(). This does not prove PRDS
#' # or automatically choose an FDR procedure.
#' trend <- trend_test(series, report = FALSE, verbose = FALSE)
#' spatial_autocorrelation(
#'   trend$stats$p, nperm = 19, seed = 1,
#'   report = FALSE, verbose = FALSE
#' )
#'
#' # --- Local inference, followed by the existing FDR module ---
#' # The statistic is descriptive by itself. The p raster is based on
#' # permutations. significant_raw uses the per-cell alpha without
#' # correcting for the number of cells tested.
#' local <- spatial_autocorrelation(
#'   r, scope = "local", nperm = 19, seed = 1,
#'   report = FALSE, verbose = FALSE
#' )
#' local$statistic
#' local$z
#' local$p
#' local$significant_raw
#'
#' # BH, BKY and BY are all supplied by fdr_correction(); q is the
#' # target FDR, not the raw per-cell alpha above.
#' local_fdr <- fdr_correction(
#'   local$p, method = c("BH", "BKY", "BY"), q = 0.05,
#'   report = FALSE, verbose = FALSE
#' )
#' c(
#'   raw_alpha = terra::global(
#'     local$significant_raw, "sum", na.rm = TRUE
#'   )[1, 1],
#'   fdr_BH = terra::global(
#'     local_fdr$rasters$sig_BH, "sum", na.rm = TRUE
#'   )[1, 1],
#'   fdr_BKY = terra::global(
#'     local_fdr$rasters$sig_BKY, "sum", na.rm = TRUE
#'   )[1, 1],
#'   fdr_BY = terra::global(
#'     local_fdr$rasters$sig_BY, "sum", na.rm = TRUE
#'   )[1, 1]
#' )
#' plot(local)
#' plot(local_fdr)
#'
#' # Both functions build the identical spatial adjacency structure
#' # underneath -- computing it once and reusing it across both calls
#' # (rather than letting each recompute it independently) matters for
#' # large rasters specifically, where that step is the expensive part.
#' nb <- prepare_cmk_neighbourhood(series, ok)
#' spatial_autocorrelation(r, precomputed_neighbourhood = nb,
#'                          nperm = 19, seed = 1,
#'                          report = FALSE, verbose = FALSE)
#' trend_test(series, precomputed_neighbourhood = nb,
#'            report = FALSE, verbose = FALSE)
#' }
#'
#' @family Spatial autocorrelation diagnostic functions
#' @export
spatial_autocorrelation <- function(x, method = c("moran", "getis_ord"),
                                     scope = c("global", "local"),
                                     connectivity = c("queen", "rook"),
                                     nperm = 99,
                                     alternative = NULL,
                                     alpha = 0.05,
                                     seed = NULL, n_cores = 1,
                                     precomputed_neighbourhood = NULL,
                                     report = TRUE,
                                     verbose = TRUE) {
  finish_timer <- .sptrends_elapsed_timer(
    "spatial_autocorrelation()", verbose)
  on.exit(finish_timer(), add = TRUE)
  method <- match.arg(method)
  scope <- match.arg(scope)
  connectivity <- match.arg(connectivity)
  if (is.null(alternative)) {
    alternative <- if (scope == "global") "greater" else "two.sided"
  } else {
    alternative <- match.arg(alternative, c("greater", "two.sided", "less"))
  }
  if (length(alpha) != 1L || !is.numeric(alpha) || !is.finite(alpha) ||
      alpha <= 0 || alpha >= 1) {
    stop("'alpha' must be one finite number strictly between 0 and 1.")
  }
  nperm <- .validate_positive_integer(nperm, "nperm")
  n_cores <- .validate_positive_integer(n_cores, "n_cores")
  seed <- .validate_seed(seed)
  if (!inherits(x, "SpatRaster")) stop("'x' must be a terra SpatRaster.")
  if (terra::nlyr(x) != 1) {
    stop("'x' must have a single layer (use x[[k]] for a stack).")
  }

  x_vals <- terra::values(x, mat = FALSE)
  if (any(!is.na(x_vals) & !is.finite(x_vals))) {
    stop("'x' must contain only finite values or NA.")
  }
  ok <- !is.na(x_vals)
  N <- sum(ok)
  if (N < 3) {
    stop("Too few valid cells to compute spatial autocorrelation ",
         "(practical minimum: 3).")
  }
  if (method == "moran" && stats::var(x_vals[ok]) == 0) {
    stop("'x' has zero variance (every valid cell has the same value) -- ",
         "Moran's I is undefined here, since its denominator is the sum ",
         "of squared deviations from the mean, which is exactly zero for ",
         "a constant raster. There is no spatial pattern to measure when ",
         "there is no variation to begin with.")
  }
  if (method == "getis_ord" && any(x_vals[ok] < 0)) {
    stop("'method = \"getis_ord\"' requires non-negative values -- ",
         "General G is a ratio of products of raw values, and a negative ",
         "input breaks its actual interpretation (see \"Methodological ",
         "background\" in ?spatial_autocorrelation for why). Use ",
         "method = \"moran\" for mixed-sign data, or transform x (e.g. ",
         "abs()) if a non-negative version of the question is what you ",
         "actually want.")
  }
  if (method == "getis_ord" && all(x_vals[ok] == 0)) {
    stop("'x' is exactly zero everywhere -- General G's denominator, the ",
         "sum of all pairwise products, is exactly zero here too, making ",
         "the statistic undefined. There is nothing to measure when ",
         "every value is zero.")
  }
  if (method == "getis_ord" && scope == "global" &&
      sum(x_vals[ok] > 0) < 2L) {
    stop("General G requires at least two positive cells; otherwise ",
         "its pairwise-product denominator is zero.")
  }

  if (verbose) message("Preparing neighbourhood topology...")
  neighbourhood <- .with_timer(
    "Neighbourhood topology",
    if (is.null(precomputed_neighbourhood)) {
      .prepare_moran_neighbourhood(x, ok, connectivity)
    } else {
      precomputed_neighbourhood
    },
    verbose)
  .validate_spatial_neighbourhood(neighbourhood, x, ok, connectivity)
  W <- neighbourhood$W
  # Accept a matching 3 by 3 precomputed object from any of this package's own
  # neighbourhood builders -- .prepare_moran_neighbourhood()'s own
  # (W, S0), prepare_cmk_neighbourhood()'s own (W, nb_count, no S0),
  # or .prepare_spatial_neighbourhood()'s own (all three) -- rather
  # than requiring S0 specifically to already be present. sum(W) for a
  # 0/1 adjacency matrix equals sum(nb_count) exactly (each is the
  # total edge count, just reached by a different route), so deriving
  # it here when missing is not an approximation.
  S0 <- if (!is.null(neighbourhood$S0)) neighbourhood$S0 else sum(W)
  if (S0 == 0) {
    stop("No valid cell has valid neighbours -- spatial autocorrelation ",
         "cannot be computed.")
  }

  if (scope == "local") {
    valid_idx <- which(ok)
    valid_vals <- x_vals[valid_idx]
    W_star <- if (method == "getis_ord") {
      W + Matrix::Diagonal(x = as.numeric(ok))
    } else {
      NULL
    }
    observed_full <- .spatial_local_statistic(
      x_vals, W, ok, method, W_star
    )
    observed <- observed_full[valid_idx]
    has_neighbour <- as.numeric(Matrix::rowSums(W))[valid_idx] > 0
    observed[!has_neighbour] <- NA_real_

    if (!is.null(seed)) set.seed(seed)
    permutation_seeds <- sample.int(
      .Machine$integer.max, nperm, replace = TRUE
    )
    if (verbose) {
      label <- if (method == "moran") "local Moran" else "local Gi*"
      message(sprintf("Computing %s permutation inference (nperm=%d)...",
                      label, nperm))
    }
    n_chunks <- min(n_cores, nperm)
    chunks <- split(
      seq_len(nperm), rep(seq_len(n_chunks), length.out = nperm)
    )
    shared_cl <- .sptrends_shared_cluster(n_cores)
    if (!is.null(shared_cl)) {
      on.exit(parallel::stopCluster(shared_cl), add = TRUE)
    }
    run_pass <- function(phase, null_mean = NULL) {
      run_chunk <- function(indices) {
        .spatial_local_chunk(
          indices, permutation_seeds, x_vals, valid_idx, valid_vals,
          W, W_star, ok, method, has_neighbour, observed, phase,
          alternative, null_mean
        )
      }
      .sptrends_parallel_lapply(
        chunks, run_chunk, n_cores = n_cores,
        export_vars = c(
          "permutation_seeds", "x_vals", "valid_idx", "valid_vals",
          "W", "W_star", "ok", "method", "has_neighbour", "observed",
          "phase", "alternative", "null_mean",
          ".spatial_local_statistic", ".spatial_local_chunk"
        ),
        export_env = environment(), packages = "Matrix",
        shared_cluster = shared_cl
      )
    }
    combine <- function(results, component) {
      Reduce(`+`, lapply(results, `[[`, component))
    }
    first <- run_pass(1L)
    null_mean <- combine(first, "total") / nperm
    null_variance <- if (nperm > 1L) {
      (combine(first, "total_sq") - nperm * null_mean^2) /
        (nperm - 1L)
    } else {
      rep(NA_real_, N)
    }
    null_variance[null_variance < 0] <- 0
    null_sd <- sqrt(null_variance)
    null_mean[!has_neighbour] <- NA_real_
    null_sd[!has_neighbour] <- NA_real_
    z_empirical <- rep(NA_real_, N)
    z_testable <- is.finite(observed) & is.finite(null_sd) & null_sd > 0
    z_empirical[z_testable] <-
      (observed[z_testable] - null_mean[z_testable]) /
      null_sd[z_testable]

    second <- if (alternative == "two.sided") {
      run_pass(2L, null_mean)
    } else {
      NULL
    }
    extreme_count <- switch(
      alternative,
      greater = combine(first, "greater"),
      less = combine(first, "less"),
      two.sided = combine(second, "extreme")
    )
    p_raw <- (extreme_count + 1) / (nperm + 1)
    p_raw[!is.finite(observed)] <- NA_real_

    significant_raw <- p_raw <= alpha

    statistic_name <- if (method == "moran") "local_moran" else "gi_star"
    out <- list(
      statistic = .spatial_raster_from_values(
        x, ok, observed, statistic_name
      ),
      p = .spatial_raster_from_values(x, ok, p_raw, "p"),
      significant_raw = .spatial_raster_from_values(
        x, ok, as.numeric(significant_raw), "significant_raw"
      ),
      method = method, scope = scope, N_tested = sum(has_neighbour),
      alternative = alternative, alpha = alpha,
      nperm = nperm, N = N, null_mean = .spatial_raster_from_values(
        x, ok, null_mean, "null_mean"
      ),
      null_sd = .spatial_raster_from_values(x, ok, null_sd, "null_sd"),
      z = .spatial_raster_from_values(x, ok, z_empirical, "z_empirical"),
      permutation_seeds = permutation_seeds
    )
    class(out) <- c("spatial_autocorrelation_local",
                    "spatial_autocorrelation", "sptrends")
    if (isTRUE(report)) {
      summary(out)
      plot(out)
    }
    return(out)
  }

  if (method == "moran") {
    calc_stat <- function(vals) {
      xbar <- mean(vals[ok])
      z <- vals - xbar
      z[!ok] <- 0
      Wz <- as.numeric(W %*% z)
      (N / S0) * sum(z * Wz) / sum(z^2)
    }
    stat_name <- "Moran's I"
  } else {
    calc_stat <- function(vals) {
      # General G = sum(w_ij * x_i * x_j) / sum(x_i * x_j), i != j (Getis
      # & Ord, 1992) -- reduced to a vectorised closed form: the
      # numerator is x'Wx (W already has a zero diagonal, matching the
      # i != j exclusion); the denominator is the sum of all pairwise
      # products, (sum(x))^2 - sum(x^2), which excludes the i == j terms
      # (x_i^2) from the full cross-product sum the same way.
      z <- vals
      z[!ok] <- 0
      Wz <- as.numeric(W %*% z)
      sum(z * Wz) / (sum(z)^2 - sum(z^2))
    }
    stat_name <- "Getis-Ord General G"
  }

  if (verbose) message(sprintf("Computing observed %s...", stat_name))
  stat_obs <- calc_stat(x_vals)
  sign_ <- if (method == "moran") {
    if (stat_obs >= 0) "positive" else "negative"
  } else {
    NA_character_
  }
  if (verbose) {
    if (method == "moran") {
      message(sprintf("      Observed I = %.4f (sign: %s)", stat_obs, sign_))
    } else {
      message(sprintf("      Observed G = %.4f", stat_obs))
    }
  }

  if (verbose) {
    message(sprintf(
      "Building the null distribution by spatial permutation (nperm=%d)...",
      nperm))
  }
  valid_idx <- which(ok)
  valid_vals <- x_vals[valid_idx]

  if (!is.null(seed)) set.seed(seed)

  if (n_cores <= 1) {
    stat_null <- numeric(nperm)
    pb <- .sptrends_progress(nperm, sprintf("%s permutations", stat_name),
                              verbose)
    for (i in seq_len(nperm)) {
      vals_perm <- x_vals
      vals_perm[valid_idx] <- sample(valid_vals)
      stat_null[i] <- calc_stat(vals_perm)
      .sptrends_progress_step(pb, i)
    }
    .sptrends_progress_close(pb)
  } else {
    if (verbose) {
      message(sprintf("      Parallel permutations over %d cores...", n_cores))
    }
    one_perm <- function(i) {
      vals_perm <- x_vals
      vals_perm[valid_idx] <- sample(valid_vals)
      calc_stat(vals_perm)
    }
    results <- .sptrends_parallel_lapply(
      seq_len(nperm), one_perm, n_cores = n_cores,
      export_vars = c("x_vals", "valid_idx", "valid_vals", "calc_stat", "W",
                       "S0", "N", "ok"),
      export_env = environment(), seed = seed,
      packages = "Matrix")
    stat_null <- unlist(results)
  }

  p <- switch(alternative,
    greater   = (sum(stat_null >= stat_obs) + 1) / (nperm + 1),
    less      = (sum(stat_null <= stat_obs) + 1) / (nperm + 1),
    two.sided = {
      p_g <- (sum(stat_null >= stat_obs) + 1) / (nperm + 1)
      p_l <- (sum(stat_null <= stat_obs) + 1) / (nperm + 1)
      min(1, 2 * min(p_g, p_l))
    }
  )

  if (verbose) message(sprintf("      p-value (%s) = %.4f", alternative, p))

  out <- list(statistic = stat_obs, method = method, scope = scope,
              sign = sign_, p = p,
              alternative = alternative, nperm = nperm,
              null_dist = stat_null, N = N)
  class(out) <- c("spatial_autocorrelation_global",
                  "spatial_autocorrelation", "sptrends")
  if (isTRUE(report)) {
    spatial_autocorrelation_summary(out)
    spatial_autocorrelation_null_plot(out)
  }
  out
}

#' @noRd
.moran_category <- function(I) {
  a <- abs(I)
  if (a < 0.1) "low" else if (a < 0.3) "moderate" else "strong"
}

#' @noRd
.print_spatial_autocorrelation <- function(x, ...) {
  if (identical(x$method, "getis_ord")) {
    cat("<Getis-Ord General G permutation test result>\n")
    cat(sprintf("G = %.4f | p-value (%s, nperm=%d) = %.4f\n",
                x$statistic, x$alternative, x$nperm, x$p))
  } else {
    cat("<Moran's I permutation test result>\n")
    cat(sprintf(
      "I = %.4f (sign: %s, category: %s*) | p-value (%s, nperm=%d) = %.4f\n",
      x$statistic, x$sign, .moran_category(x$statistic), x$alternative,
      x$nperm, x$p))
    cat("* this package's own descriptive convention -- not a recognised\n")
    cat("  disciplinary standard for Moran's I magnitude.\n")
  }
  invisible(x)
}

#' @noRd
.summary_spatial_autocorrelation <- function(object, ...) {
  spatial_autocorrelation_summary(object, ...)
}

#' @noRd
.plot_spatial_autocorrelation <- function(x, ...) {
  spatial_autocorrelation_null_plot(x, ...)
}

# Scope-specific S3 dispatch points keep global and local reporting
# structurally independent behind the common sptrends generics.
#' @noRd
.print_spatial_autocorrelation_global <- function(x, ...) {
  .print_spatial_autocorrelation(x, ...)
}

#' @noRd
.summary_spatial_autocorrelation_global <- function(object, ...) {
  .summary_spatial_autocorrelation(object, ...)
}

#' @noRd
.plot_spatial_autocorrelation_global <- function(x, ...) {
  .plot_spatial_autocorrelation(x, ...)
}

#' @noRd
.print_spatial_autocorrelation_local <- function(x, ...) {
  statistic <- if (x$method == "moran") {
    "Local Moran's I"
  } else {
    "Getis-Ord Gi*"
  }
  significant <- sum(
    terra::values(x$significant_raw, mat = FALSE) == 1,
    na.rm = TRUE
  )
  cat(sprintf("<%s permutation result>\n", statistic))
  cat(sprintf(
    "%d of %d valid cells tested | raw p <= %.3f: %d (unadjusted)\n",
    x$N_tested, x$N, x$alpha, significant
  ))
  cat("Use fdr_correction(x$p, method = c(\"BH\", \"BKY\", \"BY\")).\n")
  invisible(x)
}

#' @noRd
.summary_spatial_autocorrelation_local <- function(object, path = NULL, ...) {
  statistic <- terra::values(object$statistic, mat = FALSE)
  p_raw <- terra::values(object$p, mat = FALSE)
  significant <- terra::values(object$significant_raw, mat = FALSE)
  finite_min <- function(values) {
    values <- values[is.finite(values)]
    if (length(values)) min(values) else NA_real_
  }
  tab <- data.frame(
    method = object$method,
    scope = object$scope,
    alternative = object$alternative,
    N = object$N,
    N_tested = object$N_tested,
    nperm = object$nperm,
    statistic_min = finite_min(statistic),
    statistic_median = stats::median(statistic, na.rm = TRUE),
    statistic_max = -finite_min(-statistic),
    raw_p_min = finite_min(p_raw),
    alpha = object$alpha,
    n_significant_raw = sum(significant == 1, na.rm = TRUE),
    pct_significant_raw = round(
      100 * sum(significant == 1, na.rm = TRUE) / object$N_tested, 2)
  )
  print(tab, row.names = FALSE)
  if (!is.null(path)) {
    utils::write.csv(tab, path, row.names = FALSE)
  }
  invisible(tab)
}

#' @noRd
.plot_spatial_autocorrelation_local <- function(x, path = NULL, ...) {
  layers <- c(x$statistic, x$z, x$p, x$significant_raw)
  names(layers) <- c("Statistic", "Empirical z", "Permutation p",
                     paste("Raw p <=", x$alpha))
  plot_layers <- function() {
    terra::plot(layers, ...)
  }
  plot_layers()
  if (!is.null(path)) .save_current_plot(path, plot_layers)
  invisible(NULL)
}

#' Informal descriptive label for a Moran's I value
#'
#' **Explicit warning**: this is a descriptive convention specific to this
#' package, not a recognised disciplinary standard -- there is no
#' universal, citable threshold for Moran's I magnitude (it depends on the
#' weight matrix used). Use only as an informal aid, never as citable
#' justification in a formal write-up. `method = "moran"` results only --
#' Getis-Ord General G's natural range and expected value under H0 differ
#' enough from Moran's I that these same thresholds would not transfer
#' meaningfully (see `?spatial_autocorrelation`'s "Limitations" section).
#'
#' **Function type:** **Reporting/derived function** -- summarises or
#' plots the output of
#' another function; it does not compute any new statistic. Not
#' exported -- folded into `print()`/`summary()` of a
#' `"spatial_autocorrelation"` object (see `?print.sptrends`), the
#' category now shown there for every `method = "moran"`
#' [spatial_autocorrelation()] result automatically, without a separate
#' call.
#'
#' @param I Numeric. A Moran's I value (e.g. `result$statistic` from a
#'   `method = "moran"` result).
#' @return Invisibly, a character label (`"low"`, `"moderate"`, or
#'   `"strong"`).
#' @family Spatial autocorrelation diagnostic functions
#' @references
#' - Moran, P.A.P. (1950) Notes on Continuous Stochastic Phenomena.
#'   Biometrika, 37(1-2), 17-23. \doi{10.1093/biomet/37.1-2.17}
#' @examples
#' # A moderately positive Moran's I would be labelled by this
#' # function's own convention as follows (see the "Warning" message it
#' # prints about that label not being a disciplinary standard) --
#' # called internally wherever a categorical Moran's I label is
#' # shown, e.g. spatial_autocorrelation()'s own summary output:
#' # classify_moran(0.32)
#'
#' @keywords internal
classify_moran <- function(I) {
  message("Warning: descriptive category of this package's own convention ",
          "-- not backed")
  message("by a recognised disciplinary standard.")
  a <- abs(I)
  category <- if (a < 0.1) "low" else if (a < 0.3) "moderate" else "strong"
  message(sprintf("|I| = %.4f -> category (own convention): %s", a, category))
  invisible(category)
}

#' Summary of a spatial_autocorrelation() result
#'
#' **Function type:** **Reporting/derived function** -- summarises or
#' plots the output of
#' another function; it does not compute any new statistic. Not
#' exported -- called internally by `report = TRUE`, and reachable from
#' outside the package via `summary()`.
#'
#' @param result Output of [spatial_autocorrelation()].
#' @param path Character or `NULL`. If supplied, write the table to this
#'   CSV path.
#' @return Invisibly, a data frame.
#' @family Spatial autocorrelation diagnostic functions
#' @references
#' - Moran, P.A.P. (1950) Notes on Continuous Stochastic Phenomena.
#'   Biometrika, 37(1-2), 17-23. \doi{10.1093/biomet/37.1-2.17}
#' - Getis, A. and Ord, J.K. (1992) The Analysis of Spatial Association by
#'   Use of Distance Statistics. Geographical Analysis, 24(3), 189-206.
#'   \doi{10.1111/j.1538-4632.1992.tb00261.x}
#' @examples
#' r <- read_ordered_stack(example_data("vhp_ndvi"))[[1]]
#' moran_result <- spatial_autocorrelation(r, nperm = 19, seed = 1,
#'                                          verbose = FALSE, report = FALSE)
#'
#' # A plain-language readout of the test above: observed I, its
#' # classification, and whether it is significantly different from
#' # what random spatial arrangement would give. Called internally by
#' # summary() on a spatial_autocorrelation() result -- the public
#' # entry point is:
#' summary(moran_result)
#'
#' @keywords internal
spatial_autocorrelation_summary <- function(result, path = NULL) {
  is_getis <- identical(result$method, "getis_ord")
  stat_label <- if (is_getis) "G" else "I"
  null_min <- min(result$null_dist)
  null_max <- max(result$null_dist)
  null_mean <- mean(result$null_dist)
  null_sd <- stats::sd(result$null_dist)
  p_floor <- min(1, if (result$alternative == "two.sided") {
    2 / (result$nperm + 1)
  } else 1 / (result$nperm + 1))
  z_empirical <- if (is.finite(null_sd) && null_sd > 0) {
    (result$statistic - null_mean) / null_sd
  } else {
    NA_real_
  }

  message(sprintf("N (valid cells): %d", result$N))
  if (is_getis) {
    message(sprintf("Observed G: %.4f", result$statistic))
  } else {
    message(sprintf("Observed I: %.4f (sign: %s)", result$statistic,
                     result$sign))
  }
  category <- if (is_getis) NA_character_ else .moran_category(result$statistic)
  if (!is_getis) {
    message(sprintf(
      paste0("Category (this package's own descriptive convention, not a ",
             "recognised disciplinary standard): %s"), category))
  }
  message(sprintf(
    paste0("Null distribution (%d permutations): min=%.4f, mean=%.4f, ",
           "sd=%.4f, max=%.4f"),
    result$nperm, null_min, null_mean, null_sd, null_max))
  message(sprintf("p-value (%s, nperm=%d): %.4f",
                   result$alternative, result$nperm, result$p))
  message(sprintf("Empirical Z ((%s - null mean) / null sd): %.2f",
                   stat_label, z_empirical))
  if (abs(result$p - p_floor) < 1e-12) {
    message(sprintf(
      "Note: p is at the mathematical floor for nperm=%d (%.4f).",
      result$nperm, p_floor))
  }

  tab <- data.frame(
    metric = c("N", "method", stat_label, "sign", "category", "alternative",
               "nperm", "p_value", "p_floor", "null_min", "null_mean",
               "null_sd", "null_max", "empirical_Z"),
    value = c(result$N, result$method, round(result$statistic, 6),
              ifelse(is_getis, NA_character_, result$sign), category,
              result$alternative, result$nperm, round(result$p, 6),
              round(p_floor, 6), round(null_min, 6), round(null_mean, 6),
              round(null_sd, 6), round(null_max, 6), round(z_empirical, 4)))
  if (!is.null(path)) {
    utils::write.csv(tab, path, row.names = FALSE)
    message(sprintf("Table written to: %s", path))
  }
  invisible(tab)
}

#' Null-distribution plot for a spatial_autocorrelation() result
#'
#' **Function type:** **Reporting/derived function** -- summarises or
#' plots the output of
#' another function; it does not compute any new statistic. Not
#' exported -- called internally by `report = TRUE`, and reachable from
#' outside the package via `plot()`.
#'
#' @inheritParams spatial_autocorrelation_summary
#' @return `NULL`, invisibly.
#' @family Spatial autocorrelation diagnostic functions
#' @references
#' - Moran, P.A.P. (1950) Notes on Continuous Stochastic Phenomena.
#'   Biometrika, 37(1-2), 17-23. \doi{10.1093/biomet/37.1-2.17}
#' - Getis, A. and Ord, J.K. (1992) The Analysis of Spatial Association by
#'   Use of Distance Statistics. Geographical Analysis, 24(3), 189-206.
#'   \doi{10.1111/j.1538-4632.1992.tb00261.x}
#' @examples
#' r <- read_ordered_stack(example_data("vhp_ndvi"))[[1]]
#' moran_result <- spatial_autocorrelation(r, nperm = 19, seed = 1,
#'                                          verbose = FALSE, report = FALSE)
#'
#' # The null distribution (I values from randomly shuffled data), with
#' # the observed I marked -- how unusual is it compared to pure chance?
#' # Called internally by plot() on a spatial_autocorrelation() result
#' # -- the public entry point is:
#' plot(moran_result)
#'
#' @keywords internal
spatial_autocorrelation_null_plot <- function(result, path = NULL) {
  is_getis <- identical(result$method, "getis_ord")
  stat_label <- if (is_getis) "G" else "I"
  stat_name <- if (is_getis) "Getis-Ord General G" else "Moran's I"
  null_range <- range(result$null_dist)
  span <- diff(null_range)
  if (span == 0) span <- 0.01
  padding <- span * 0.15
  xlim_null <- c(null_range[1] - padding, null_range[2] + padding)
  out_of_range <- result$statistic < xlim_null[1] ||
    result$statistic > xlim_null[2]

  plot_null <- function() {
    old_par <- graphics::par(mar = c(5, 4, 4, 3) + 0.1)
    on.exit(graphics::par(old_par), add = TRUE)

    if (out_of_range) {
      graphics::hist(result$null_dist, breaks = 20, col = "steelblue",
           border = "white",
           main = sprintf("Null distribution of %s (spatial permutation)",
                           stat_name),
           sub = sprintf(
             "Observed %s = %.4f (outside this range -- see annotation)",
             stat_label, result$statistic),
           xlab = sprintf("%s (under H0)", stat_label), xlim = xlim_null,
           cex.sub = 0.9)
      usr <- graphics::par("usr")
      right_side <- result$statistic > xlim_null[2]
      arrow <- if (right_side) "->" else "<-"
      x_text <- if (right_side) usr[2] else usr[1]
      txt <- sprintf("%s obs. %s = %.4f", arrow, stat_label,
                     result$statistic)
      graphics::text(x = x_text, y = usr[4] * 0.85, labels = txt,
           col = "red", font = 2, adj = if (right_side) 1 else 0)
    } else {
      graphics::hist(result$null_dist, breaks = 20, col = "steelblue",
           border = "white",
           main = sprintf("Null distribution of %s (spatial permutation)",
                           stat_name),
           xlab = sprintf("%s (under H0)", stat_label),
           xlim = range(c(xlim_null, result$statistic)))
      graphics::abline(v = result$statistic, col = "red", lwd = 2)
      graphics::legend("topright",
             legend = sprintf("Observed %s", stat_label), col = "red",
             lwd = 2, bty = "n", inset = 0.02, cex = 0.9)
    }
  }
  plot_null()

  if (!is.null(path)) .save_current_plot(path, plot_null)
  invisible(NULL)
}
