#' Conditional Location Scale Empirical Bayes (CLOSE)
#'
#' This function takes in estimates and standard errors and computes
#' CLOSE-GAUSS and CLOSE-NPMLE empirical Bayes posterior means.
#'
#' @param estimates A vector of noisy estimates for some underlying parameters. It is assumed
#' that these estimates are Gaussian with standard deviation equal to `standard_errors`
#' and means equal to the underlying parameters.
#' @param standard_errors A vector of standard errors for the estimates.
#' @param conditional_mean A vector of estimated conditional mean.
#' `E[estimates | standard_error]` for the parameters. If `NULL`, this will be estimated by
#' `local_linear_regression_conditional_moments`.
#' @param conditional_std A vector of estimated conditional standard deviation
#' for the parameters: `Sqrt(Var(estimates | standard_errors) - standard_errors^2)`.
#' If `NULL`, this will be estimated by `local_linear_regression_conditional_moments`.
#' @param truncate A small number. If the estimated conditional standard deviation is
#' less than truncate, the observations are ignored.
#' @param grid A grid of values to use for the NPMLE estimate of the conditional distribution
#' of the parameters. If `NULL`, a default grid is used. The default grid is
#' ```
#' grid <- seq(-6, 6, length.out=400)
#' minx <- min(-7, min(transformed_estimates[good_obs]))
#' max <- max(7, max(transformed_estimates[good_obs]))
#' grid <- c(seq(minx, -6, length.out=100)[1:99], grid, seq(6, max, length.out=100)[2:100])
#' ```
#' @param ... Additional arguments to be passed to `local_linear_regression_conditional_moments`.
#' In particular, if one wishes to estimate the conditional moment based on some transformation of
#' `standard_errors`, one should pass in, e.g., `morphed_standard_errors = log(standard_errors)`
#'  here.
#' @return A list with the following elements:
#' \itemize{
#'   \item{`close_gauss`} The posterior means according to CLOSE-GAUSS
#'   \item{`close_npmle`} The posterior means according to CLOSE-NPMLE
#'   \item{`conditional_mean`}  The estimated conditional mean of the parameters
#'   \item{`conditional_std`}  The estimated conditional standard deviation of the parameters
#'   \item{`npmle`}  The metadata on NPMLE returned by `REBayes::GLmix`
#' }
#' @examples
#' n <- 10000
#' standard_errors <- 0.1 + runif(n)
#' thetas <- standard_errors + 0.3 * rnorm(n)
#' estimates <- thetas + standard_errors * rnorm(n)
#' close_results <- compute_close(estimates, standard_errors)
#'
#' @export
compute_close <- function(
    estimates, standard_errors, conditional_mean = NULL,
    conditional_std = NULL, truncate = 1e-8, grid = NULL, ...) {
    if (is.null(conditional_mean) || is.null(conditional_std)) {
        result <- local_linear_regression_conditional_moments(
            estimates = estimates, standard_errors = standard_errors, ...
        )
        conditional_mean <- result$conditional_mean
        conditional_std <- result$conditional_std
    }

    close_gauss_weight <- (conditional_std^2) / (conditional_std^2 + standard_errors^2)
    close_gauss <- estimates * close_gauss_weight + conditional_mean * (1 - close_gauss_weight)

    transformed_estimates <- (estimates - conditional_mean) / conditional_std
    transformed_ses <- standard_errors / conditional_std

    good_obs <- (conditional_std > truncate) & (!is.na(conditional_std))

    # empty array of size length(estimates)
    tau_posterior <- estimates * 0

    if (is.null(grid)) {
        # Default grid
        grid <- seq(-6, 6, length.out = 400)
        minx <- min(-7, min(transformed_estimates[good_obs]))
        max <- max(7, max(transformed_estimates[good_obs]))
        grid <- c(
            seq(minx, -6, length.out = 100)[1:99], grid, seq(6, max, length.out = 100)[2:100]
        )
    }

    npmle <- REBayes::GLmix(
        x = transformed_estimates[good_obs],
        sigma = transformed_ses[good_obs], v = grid
    )
    tau_posterior[good_obs] <- npmle$dy
    close_npmle <- tau_posterior * conditional_std + conditional_mean

    return(
        list(
            close_gauss = close_gauss,
            close_npmle = close_npmle,
            conditional_mean = conditional_mean,
            conditional_std = conditional_std,
            npmle = npmle
        )
    )
}
