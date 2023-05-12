# CLOSE:  Conditional Location-Scale Empirical Bayes

The `close` R package implements the empirical Bayes posterior mean estimators developed
in [Chen (2023)](https://arxiv.org/abs/2212.14444).

## Installation

```R
library(devtools)
install_github("close", "jiafengkevinchen")
```

## Summary of method

Let $Y_i \sim N(\theta_i, \sigma_i^2)$ denote some noisy estimates for parameters $\theta_i$
with known standard errors $\sigma_i$. We observe $(Y_i, \sigma_i)$ and would like to make
decisions based on $\theta_i$.

CLOSE assumes that $\theta \mid \sigma$ is a location-scale family with location parameter
$m_0$, scale parameter $s_0$, and shape parameter $G_0$.
$$ P(\theta_i < t \mid \sigma_i) = G_0\left( \frac{t - m_0(\sigma_i)}{s_0(\sigma_i)} \right).$$

CLOSE proposes to estimate $m_0, s_0$ via nonparametric regression.
The method CLOSE-NPMLE additionally estimates $G_0$ via
nonparametric maximum likelihood (see [`REBayes::GLmix`](https://www.jstatsoft.org/article/view/v082i08)). The method CLOSE-Gauss assumes $G_0$ is $N(0,1)$.

### Usage

The main function is `compute_close`. The user may specify estimated $m_0$ and $s_0$.
If they are not specified, then they are estimated by local linear regression specified in
`local_linear_regression_conditional_moments`.

As an example
```R
library(close)  # Depends on REBayes and nprobust. In particular, REBayes require a MOSEK installation

set.seed(123)
n <- 10000
standard_errors <- 0.1 + runif(n)

# thetas are correlated with standard errors
thetas <- standard_errors + 0.3 * rnorm(n)

# Y ~ N(theta, sigma^2)
estimates <- thetas + standard_errors * rnorm(n)

# close_results$close_gauss contains estimated posterior means via CLOSE-Gauss
# close_results$close_npmle contains estimated posterior means via CLOSE-NPMLE
close_results <- compute_close(estimates, standard_errors)
```
