#' Wrapper for `nprobust::lpbwselect`
#'
#' This function wraps `nprobust::lpbwselect` and returns the selected
#' optimal bandwidth as a function of data and kernel.
#' Only supports local linear regression. Supports using a
#' subset of the data to estimate the optimal bandwidth by `n_real`.
#' That is, the bandwidth is equal to `h * (n_real / length(x))^-0.2`.
#' If `length(x) < 40`, then the bandwidth
#' `h <- 1.06 * sd(x) * (length(x) + 3)^-0.2`
#' is returned, and a warning raised.
#'
#' @param x A vector of the x-values.
#' @param y A vector of the y-values.
#' @param n_real An integer equal to the sample size. Optional and defaults to `length(x)`.
#' @param kernel A string denoting a kernel type. See `nprobust::lpbwselect`.
#' Optional and defaults to `epa` for the Epanechnikov kernel.
#' @param bwselect A string denoting a bandwidth selection method.
#' See `nprobust::lpbwselect`. Defaults to `imse-dpi`.
#' @param ... Additional arguments to be passed to `nprobust::lpbwselect`.
#'
#' @return A list with the following elements:
#' \itemize{
#'  \item{`h`} The selected bandwidth.
#'  \item{`kernel`} The kernel used.
#' }
#'
#' @examples
#' x <- runif(100)
#' y <- x + rnorm(100)
#' bw_result <- llr_bwselect_wrapper(x, y)

#' @export
llr_bwselect_wrapper <- function(x, y, n_real = NULL, kernel = "epa", bwselect = "imse-dpi", ...) {
    exponent <- -0.2
    if (length(x) > 40) {
        obj <- nprobust::lpbwselect(
            x = x,
            y = y,
            p = 1,
            deriv = 0,
            kernel = kernel,
            bwselect = bwselect,
            interior = TRUE,
            ...
        )
        h <- getElement(obj$bws[1, ], "h")
    } else {
        warning("Too few observations for bandwidth selection, using Silverman's rule for LLR.")
        h <- 1.06 * sd(x) * (length(x) + 3)^exponent
    }

    if (is.null(n_real)) {
        n_real <- length(x)
    }


    h_const <- h / (length(x)^exponent)
    h_new <- h_const * (n_real^exponent)

    return(list(h = h_new, kernel = kernel))
}

#' Use a subset of the observations to select bandwidth
#'
#' This function is a wrapper for `llr_bwselect_wrapper` that
#' first picks `frac` of the observations at random and
#' pass to `llr_bwselect_wrapper`.
#' @param x A vector of the x-values.
#' @param y A vector of the y-values.
#' @param frac A number between 0 and 1 denoting the fraction of observations to use.
#' @param ... Additional arguments to be passed to `llr_bwselect_wrapper`.
#' @return A list with the following elements:
#' \itemize{
#'  \item{`h`} The selected bandwidth.
#'  \item{`kernel`} The kernel used.
#' }
#' @examples
#' x <- runif(10000)
#' y <- x + rnorm(10000)
#' bw_result <- fast_bandwidth_select(x, y, frac = 0.2)
#' @export
fast_bandwidth_select <- function(x, y, frac = 0.2, ...) {
    n <- length(x)
    idx <- sample(n, floor(frac * n), replace = FALSE)
    results <- llr_bwselect_wrapper(x[idx], y[idx], n_real = n, ...)
    return(results)
}

#' Local linear regression
#'
#' This function implements local linear regression with a given
#' bandwidth and evaluation points. The result is linearly interpolated
#' from the evaluation points to all x-values.
#'
#' @param x A vector of the x-values.
#' @param y A vector of the y-values.
#' @param h A number denoting the bandwidth.
#' @param evals A vector of evaluation points. Optional and defaults to x.
#' @param kernel A string denoting a kernel type. See `nprobust::lpbwselect`.
#' Defaults to "epa" for the Epanechnikov kernel.
#' @return A list with the following elements:
#' \itemize{
#'     \item{`f_hat`} A vector that contains the estimated value
#'      of the conditional expectation at each x-value.
#'    \item{`eff_sample_size`} A number equal to the average effective
#'  sample size averaged over evaluation points. The effective sample
#'  size is defined as the inverse sum of squares of the weights used
#'  at the evaluation point. That is, if `f_hat(x) = sum_i w_i(x) y_i`,
#' then `eff_sample_size_i = 1 / sum_i w_i(x)^2`.
#' }
#'
#' @examples
#' x <- runif(100)
#' y <- x + rnorm(100)
#' bw_result <- llr_bwselect_wrapper(x, y)
#' llr_result <- local_linear_smoothing(x, y, bw_result$h)
#' @export
local_linear_smoothing <- function(x, y, h, evals = NULL, kernel = "epa") {
    kernels <- list(
        epa = function(x) pmax(0, (1 - x^2) * 3 / 4),
        gau = function(x) exp(-(x^2) / 2) / sqrt(2 * pi),
        tri = function(x) pmax(0, 1 - abs(x)),
        uni = function(x) ifelse(abs(x) <= 1, 0.5, 0)
    )

    if (is.null(evals)) {
        evals <- x
        interp <- FALSE
    } else {
        interp <- TRUE
    }
    if (!(kernel %in% names(kernels))) {
        warning("Unknown kernel, using Epanechnikov kernel instead.")
        kernel <- "epa"
    }

    # matrix of xx = [1, x, x^2, ..., x^p]
    xx <- matrix(rep(x, 2), ncol = 2)^matrix(rep(0:1, each = length(x)), nrow = length(x))

    # analogous to xx for evaluation points
    evals_xx <- matrix(rep(evals, 2), ncol = 2)^matrix(
        rep(0:1, each = length(evals)),
        nrow = length(evals)
    )

    xxt_matrices <- outer_broadcast(xx)
    dist_mat <- outer(evals, x, FUN = "-")
    kernel_mat <- matrix(kernels[[kernel]](dist_mat / h), nrow = dim(dist_mat)[1])
    llr_point <- llr_point_maker(xxt_matrices, xx, kernel_mat, evals_xx, y)

    result <- t(sapply(seq_along(evals), llr_point))
    f_hat_eval <- unlist(result[, "f_hat"])
    eff_sample_size <- mean(unlist(result[, "eff_sample_size"]))

    if (interp) {
        f_hat <- approx(x = evals, y = f_hat_eval, xout = x)$y
    } else {
        f_hat <- f_hat_eval
    }

    return(list(f_hat = f_hat, eff_sample_size = eff_sample_size))
}


outer_broadcast <- function(x) {
    # for x an n x p array, return an n x p x p array where each slice is x[, i] %*% t(x[, i])
    n <- dim(x)[1]
    p <- dim(x)[2]
    x3d1 <- array(x, dim = c(n, p, 1))
    x3d1 <- array(rep(x3d1, times = p, along = 3), dim = c(n, p, p))
    x3d2 <- aperm(x3d1, c(1, 3, 2))

    return(x3d1 * x3d2)
}


llr_point_maker <- function(xxt_matrices, xx, kernel_mat, evals_xx, y) {
    llr_point <- function(idx) {
        eval_point <- evals_xx[idx, ]
        weight <- kernel_mat[idx, ]

        # p x p
        singular <- FALSE
        tryCatch(
            expr = {
                sum_w_xxt_inverse <- solve(colSums(xxt_matrices * weight, dims = 1))
            },
            error = function(e) {
                singular <<- TRUE
            }
        )
        if (singular) {
            return(list(f_hat = NA, eff_sample_size = NA))
        }

        w_xt <- t(weight * xx) # p x n

        hat_mat_row <- t(eval_point) %*% (sum_w_xxt_inverse %*% w_xt) # 1 x n

        f_hat <- hat_mat_row %*% y
        eff_sample_size <- 1 / sum(hat_mat_row^2)

        return(list(f_hat = f_hat, eff_sample_size = eff_sample_size))
    }
    return(llr_point)
}

#' Local linear regression with bandwidth estimated on
#' a subset of the data and a small number of evaluation points
#'
#' This function implements local linear regression with (a) bandwidth
#'  selected by `fast_bandwidth_select` and (b) a small number of evaluation points.
#' The estimated function is a linear interpolation of the estimated
#'  values at the evaluation points.
#'
#' @param x A numeric vector of x-values
#' @param y A numeric vector of y-values
#' @param frac The fraction of the data to use for bandwidth selection.
#' See `fast_bandwidth_select` for details.
#' @param kernel The kernel to use. See `nprobust::lpbwselect` for details. Defaults to `epa`.
#' @param bwselect The bandwidth selection method to use.
#' See `nprobust::lpbwselect` for details. Defaults to `imse-dpi`.
#' @param ngrid The number of evaluation points to use. Defaults to 500.
#' @param ... Additional arguments to pass to `fast_bandwidth_select`.
#' @return A list with the following elements:
#' \itemize{
#'     \item{`f_hat`} A vector that contains the estimated value
#'      of the conditional expectation at each x-value.
#'    \item{`eff_sample_size`} A number equal to the average effective
#'  sample size averaged over evaluation points. The effective sample
#'  size is defined as the inverse sum of squares of the weights used
#'  at the evaluation point. That is, if
#' `f_hat(x) = sum_i w_i(x) y_i`,
#' then `eff_sample_size_i = 1 / sum_i w_i(x)^2`.
#' }
#' @examples
#' x <- runif(10000)
#' y <- x + rnorm(10000)
#' llr_result <- local_linear_regression_fast(x, y)

#' @export
local_linear_regression_fast <- function(
    x, y, frac = 0.2,
    kernel = "epa",
    bwselect = "imse-dpi",
    ngrid = 500,
    ...) {
    h <- fast_bandwidth_select(x, y, frac = frac, kernel = kernel, bwselect = bwselect, ...)$h
    evals <- seq(min(x), max(x), length.out = ngrid)
    return(local_linear_smoothing(x = x, y = y, h = h, evals = evals, kernel = kernel))
}

#' Local linear regression
#'
#' This function implements local linear regression .
#'
#' @param x A numeric vector of x-values
#' @param y A numeric vector of y-values
#' @param kernel The kernel to use. See `nprobust::lpbwselect` for details. Defaults to `epa`.
#' @param bwselect The bandwidth selection method to use.
#' See `nprobust::lpbwselect` for details. Defaults to `imse-dpi`.
#' @param ... Additional arguments to pass to `fast_bandwidth_select`.
#' @return A list with the following elements:
#' \itemize{
#'     \item{`f_hat`} A vector that contains the estimated value
#'      of the conditional expectation at each x-value.
#'    \item{`eff_sample_size`} A number equal to the average effective
#'  sample size averaged over evaluation points. The effective sample
#'  size is defined as the inverse sum of squares of the weights used
#'  at the evaluation point. That is, if
#' `f_hat(x) = sum_i w_i(x) y_i`,
#' then `eff_sample_size_i = 1 / sum_i w_i(x)^2`.}
#' @examples
#' x <- runif(10000)
#' y <- x + rnorm(10000)
#' llr_result <- local_linear_regression(x, y)
#' @export
local_linear_regression <- function(
    x, y,
    kernel = "epa",
    bwselect = "imse-dpi",
    ...) {
    h <- fast_bandwidth_select(x, y, kernel = kernel, frac = 1, bwselect = bwselect, ...)$h
    return(local_linear_smoothing(x = x, y = y, h = h, kernel = kernel))
}

#' Estimation of conditional moments of `theta | sigma` using local linear regression
#'
#' This function estimates the conditional mean and variance of `theta | sigma`
#'  using local linear regression.
#'
#' @param estimates A numeric vector of noisy estimates (`Y`).
#' @param standard_errors A numeric vector of standard errors (`sigma`).
#' @param morphed_standard_errors A numeric vector of transformed
#' standard errors (e.g. `log(sigma)`) to estimate the conditional moments with.
#' Defaults to `standard_errors`.
#' @param kernel The kernel to use. See `nprobust::lpbwselect` for details. Defaults to `epa`.
#' @param bwselect The bandwidth selection method to use. Defaults to `imse-dpi`.
#' @param fast Whether to use `local_linear_regression_fast` or `local_linear_regression`.
#' @param variance_fit_type The type of variance estimator to use. The choices are
#' `squared_residual` (the default), `fit_conditional`, and `difference_of_squares`.
#' \itemize{
#'     \item{`squared_residual`} The variance estimator is an analogue estimator of
#'          ```Var(theta | sigma) = E[(Y - E[Y | sigma])^2 | sigma] - sigma^2```.
#'     \item{`fit_conditional`} The variance estimator is an analogue estimator of
#'          ```Var(theta | sigma) = E[(Y - E[Y | sigma])^2 - sigma^2 | sigma]```.
#'     \item{`difference_of_squares`} The variance estimator is an analogue estimator of
#'          ```Var(theta | sigma) = E[Y^2 | sigma^2] - E[E[Y | sigma]^2 | sigma] - sigma^2```.
#' }
#'
#' @param truncate Whether to truncate based on
#' the Kubokawa, Robert, and Saleh (Canadian Journal of Statistics, 1993) heuristic.
#' Defaults to `TRUE`.
#' If TRUE, then the conditional variance is truncated below at
#' `2/(eff_sample_size + 2) * E[(Y - E[Y | sigma])^2 | sigma]`,
#' where `E[(Y - E[Y | sigma])^2 | sigma]` is an estimated analogue
#' Otherwise, the conditional variance is truncated below at zero.
#' @param ... Additional arguments to pass to
#' `local_linear_regression_fast` or `local_linear_regression`.
#' @return A list with the following elements:
#' \itemize{
#'    \item{`conditional_mean`} A numeric vector of the estimated conditional mean
#'   \item{`conditional_std`} A numeric vector of the estimated conditional standard deviation
#'  \item{`eff_sample_size`} A number equal to the average effective
#' sample size averaged over evaluation points. See `local_linear_regression`.
#' }
#' @examples
#' n <- 10000
#' standard_errors <- 0.1 + runif(n)
#' thetas <- standard_errors + 0.3 * rnorm(n)
#' estimates <- thetas + standard_errors * rnorm(n)
#' conditional_moments <- local_linear_regression_conditional_moments(
#'     estimates = estimates, standard_errors = standard_errors
#' )
#' @export
local_linear_regression_conditional_moments <- function(
    estimates,
    standard_errors,
    morphed_standard_errors = NULL,
    kernel = "epa",
    bwselect = "imse-dpi",
    fast = TRUE,
    variance_fit_type = "squared_residual",
    truncate = TRUE,
    ...) {
    if (is.null(morphed_standard_errors)) {
        morphed_standard_errors <- standard_errors
    }
    if (fast) {
        regression <- local_linear_regression_fast
    } else {
        regression <- local_linear_regression
    }
    conditional_mean <- regression(
        x = morphed_standard_errors, y = estimates, kernel = kernel, bwselect = bwselect, ...
    )$f_hat

    if (variance_fit_type == "fit_conditional") {
        difference <- (estimates - conditional_mean)^2 - standard_errors^2
        smoothed_cond_var_result <- regression(
            morphed_standard_errors, difference,
            kernel = kernel, bwselect = bwselect, ...
        )
        smoothed_cond_var <- smoothed_cond_var_result$f_hat
        effective_sample_size <- smoothed_cond_var_result$eff_sample_size
    } else if (variance_fit_type == "difference_of_squares") {
        squared <- estimates^2
        smoothed_squared_result <- regression(
            morphed_standard_errors,
            squared,
            kernel = kernel,
            bwselect = bwselect,
            ...
        )
        smoothed_squared_residual <- smoothed_squared_result$f_hat - conditional_mean^2
        effective_sample_size <- smoothed_squared_result$eff_sample_size
    } else {
        if (variance_fit_type != "squared_residual") {
            warning("variance_fit_type must be one of `squared_residual`, `fit_conditional`, or `difference_of_squares`. Defaulting to `squared_residual`.")
        }
        squared_residual <- (estimates - conditional_mean)^2
        smoothed_squared_residual_result <- regression(
            morphed_standard_errors,
            squared_residual,
            kernel = kernel,
            bwselect = bwselect,
            ...
        )
        smoothed_squared_residual <- smoothed_squared_residual_result$f_hat
        effective_sample_size <- smoothed_squared_residual_result$eff_sample_size
    }

    if (variance_fit_type == "fit_conditional") {
        conditional_var <- pmax(smoothed_cond_var, 0)
        smoothed_squared_residual <- conditional_var + standard_errors^2
    } else {
        conditional_var <- pmax(smoothed_squared_residual - standard_errors^2, 0)
    }

    if (truncate) {
        truncation_factor <- 2 / (effective_sample_size + 2)
        truncation_point <- truncation_factor * smoothed_squared_residual
    }

    conditional_std <- (pmax(conditional_var, truncation_point))^0.5
    return(list(conditional_mean = conditional_mean, conditional_std = conditional_std))
}
