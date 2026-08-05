# Internal engines for calibrated spatiotemporal simulation -----------------

#' Correlation function used by the formal spatial simulator
#' @noRd
.simulation_correlation <- function(distance, model, rho, range, smoothness) {
  if (model == "independent") return(as.numeric(distance == 0))
  if (model == "exponential") {
    return(exp(log(rho) * distance))
  }
  if (model == "gaussian") {
    return(exp(log(rho) * distance^2))
  }
  scaled <- distance / range
  answer <- distance
  answer[] <- 1
  positive <- scaled > 0
  answer[positive] <- (2^(1 - smoothness) / gamma(smoothness)) *
    scaled[positive]^smoothness *
    besselK(scaled[positive], nu = smoothness)
  answer
}

#' Draw a Gaussian random field by circulant embedding and FFT
#' @noRd
.simulate_gaussian_field_exact <- function(nrow, ncol, model, rho, range,
                                           smoothness, sd) {
  coordinates <- cbind(
    rep(seq_len(nrow), each = ncol),
    rep(seq_len(ncol), times = nrow)
  )
  squared_norm <- rowSums(coordinates^2)
  squared_distance <- outer(squared_norm, squared_norm, "+") -
    2 * tcrossprod(coordinates)
  distance <- sqrt(pmax(squared_distance, 0))
  covariance <- sd^2 * .simulation_correlation(
    distance, model, rho, range, smoothness)
  decomposition <- eigen(covariance, symmetric = TRUE)
  tolerance <- 1e-10 * max(1, max(abs(decomposition$values)))
  if (min(decomposition$values) < -tolerance) {
    stop("The requested spatial covariance matrix is not positive ",
         "semidefinite.")
  }
  eigenvalues <- pmax(decomposition$values, 0)
  draw <- decomposition$vectors %*%
    (sqrt(eigenvalues) * stats::rnorm(nrow * ncol))
  matrix(as.vector(draw), nrow = nrow, ncol = ncol, byrow = TRUE)
}

#' Draw a Gaussian random field by circulant embedding and FFT
#' @noRd
.simulate_gaussian_field <- function(nrow, ncol, model, rho, range,
                                     smoothness, sd, max_expansions = 6L) {
  if (model == "independent") {
    return(matrix(stats::rnorm(nrow * ncol, sd = sd), nrow, ncol))
  }

  base_rows <- max(2L, 2L * nrow)
  base_cols <- max(2L, 2L * ncol)
  for (expansion in 0:max_expansions) {
    embed_rows <- 2^ceiling(log2(base_rows * 2^expansion))
    embed_cols <- 2^ceiling(log2(base_cols * 2^expansion))
    row_distance <- pmin(0:(embed_rows - 1L),
                         embed_rows - 0:(embed_rows - 1L))
    col_distance <- pmin(0:(embed_cols - 1L),
                         embed_cols - 0:(embed_cols - 1L))
    distance <- outer(row_distance^2, col_distance^2, "+")^0.5
    covariance <- sd^2 * .simulation_correlation(
      distance, model, rho, range, smoothness)
    eigenvalues <- Re(stats::fft(covariance))
    tolerance <- 1e-10 * max(1, max(abs(eigenvalues)))
    if (min(eigenvalues) >= -tolerance) {
      eigenvalues[eigenvalues < 0] <- 0
      white <- matrix(stats::rnorm(embed_rows * embed_cols),
                      embed_rows, embed_cols)
      spectrum <- sqrt(eigenvalues) * stats::fft(white)
      field <- Re(stats::fft(spectrum, inverse = TRUE)) /
        (embed_rows * embed_cols)
      return(field[seq_len(nrow), seq_len(ncol), drop = FALSE])
    }
  }
  if (nrow * ncol <= 2500L) {
    return(.simulate_gaussian_field_exact(
      nrow, ncol, model, rho, range, smoothness, sd))
  }
  stop("The requested covariance could not be embedded with a non-negative ",
       "spectrum, and the grid is too large for the exact covariance ",
       "fallback. Increase 'spatial_range', reduce 'spatial_smoothness', ",
       "reduce the grid size, or use another spatial model.")
}

#' Draw one innovation field with the requested marginal distribution
#' @noRd
.simulate_spatial_innovation <- function(nrow, ncol, model, rho, range,
                                         smoothness, sd, dist, df) {
  gaussian <- .simulate_gaussian_field(
    nrow, ncol, model, rho, range, smoothness, sd = 1)
  if (dist == "gaussian") return(sd * gaussian)
  if (df <= 2) {
    stop("'t_df' must be > 2 for the t-distribution to have finite variance.")
  }
  uniforms <- stats::pnorm(gaussian)
  uniforms <- pmin(pmax(uniforms, .Machine$double.eps),
                   1 - .Machine$double.eps)
  sd * sqrt((df - 2) / df) * stats::qt(uniforms, df = df)
}

#' Build an exact geometric signal mask
#' @noRd
.simulation_signal_mask <- function(nrow, ncol, shape, size, location,
                                    angle, axis_ratio, custom_mask = NULL) {
  ncell <- nrow * ncol
  if (shape == "custom") {
    if (is.null(custom_mask)) {
      stop("'custom_mask' is required when trend_shape = 'custom'.")
    }
    if (inherits(custom_mask, "SpatRaster")) {
      custom_mask <- terra::values(custom_mask, mat = FALSE)
    } else if (is.matrix(custom_mask)) {
      if (!all(dim(custom_mask) == c(nrow, ncol))) {
        stop("A matrix 'custom_mask' must have dimensions nrow by ncol.")
      }
      custom_mask <- as.vector(t(custom_mask))
    }
    if (length(custom_mask) != ncell) {
      stop("'custom_mask' must contain exactly nrow * ncol values.")
    }
    return(!is.na(custom_mask) & custom_mask != 0)
  }

  rows <- rep(seq_len(nrow), each = ncol)
  cols <- rep(seq_len(ncol), times = nrow)
  if (location == "centre") {
    centre_row <- floor((nrow + 1) / 2)
    centre_col <- floor((ncol + 1) / 2)
  } else {
    centre_row <- sample.int(nrow, 1L)
    centre_col <- sample.int(ncol, 1L)
  }
  if (length(size) == 1L) size <- rep(size, 2L)
  height <- size[[1]]
  width <- size[[2]]
  theta <- angle * pi / 180
  row_offset <- rows - centre_row
  col_offset <- cols - centre_col
  rotated_col <- col_offset * cos(theta) + row_offset * sin(theta)
  rotated_row <- -col_offset * sin(theta) + row_offset * cos(theta)

  if (shape %in% c("block", "square", "rectangle") &&
      angle %% 180 == 0) {
    height <- max(1L, as.integer(round(height)))
    width <- max(1L, as.integer(round(width)))
    row_start <- floor(centre_row - (height - 1L) / 2)
    col_start <- floor(centre_col - (width - 1L) / 2)
    return(rows >= row_start & rows < row_start + height &
             cols >= col_start & cols < col_start + width)
  }
  if (shape %in% c("block", "square", "rectangle")) {
    return(abs(rotated_row) <= height / 2 &
             abs(rotated_col) <= width / 2)
  }
  if (shape == "ellipse") {
    semi_row <- max(height / 2, .Machine$double.eps)
    semi_col <- max(width * axis_ratio / 2, .Machine$double.eps)
    return((rotated_row / semi_row)^2 + (rotated_col / semi_col)^2 <= 1)
  }
  rep(TRUE, ncell)
}
