# High-dimensional binomial/logistic simulation; base R only.
# References:
# Breheny & Huang (2011), Annals of Applied Statistics 5:232-253.
# https://doi.org/10.1214/10-AOAS388 (Section 5.2, Setting 1).
# Fan & Song (2010), Annals of Statistics 38:3567-3604.
# https://doi.org/10.1214/10-AOS798 (Sections 7 and 7.1).
# Huang, Ma & Zhang (2008), Iowa Technical Report 392, Section 4.2.
# Only bh2011() below follows the first paper's independent-design,
# four-signal structure; its default p=1000 and amplitude=1 are user choices.
# The other presets are proposed extensions, not exact paper reproductions.

simulate_binomial <- function(n = 200L, p = 1000L,
                              beta_active = rep(c(1, -1), 5),
                              active = seq_along(beta_active),
                              design = c("ar1", "independent", "equicor", "block"),
                              rho = 0.5, block_size = 20L,
                              intercept = 0, prevalence = NULL,
                              signal_sd = NULL, trials = 1L,
                              n_test = 2000L, seed = NULL) {
  design <- match.arg(design)
  is_count <- function(z, lower) {
    length(z) == 1L && is.finite(z) && z >= lower && z == floor(z)
  }
  stopifnot(is_count(n, 1), is_count(p, 1), is_count(n_test, 0),
            is_count(block_size, 1), is_count(trials, 1),
            length(rho) == 1L, is.finite(rho), rho >= 0, rho < 1,
            length(intercept) == 1L, is.finite(intercept),
            length(active) == length(beta_active),
            all(is.finite(beta_active)), all(is.finite(active)),
            all(active == floor(active)), all(active >= 1 & active <= p),
            !anyDuplicated(active))
  if (!is.null(seed)) set.seed(seed)
  beta <- numeric(p)
  beta[active] <- beta_active

  # Work only with the active covariance matrix: no p-by-p allocation.
  if (length(active)) {
    sigma <- switch(design,
      independent = diag(length(active)),
      ar1 = rho ^ abs(outer(active, active, "-")),
      equicor = matrix(rho, length(active), length(active)),
      block = rho * outer((active - 1L) %/% block_size,
                          (active - 1L) %/% block_size, "=="))
    diag(sigma) <- 1
    variance <- drop(crossprod(beta_active, sigma %*% beta_active))
  } else {
    variance <- 0
  }
  if (!is.null(signal_sd)) {
    stopifnot(length(signal_sd) == 1L, is.finite(signal_sd),
              signal_sd > 0, variance > 0)
    beta <- beta * signal_sd / sqrt(variance)
    variance <- signal_sd ^ 2
  }

  # Exact population calibration up to numerical integration error,
  # valid for these Gaussian designs. Never calibrate on the test sample.
  if (!is.null(prevalence)) {
    stopifnot(length(prevalence) == 1L, is.finite(prevalence),
              prevalence > 0, prevalence < 1)
    population_mean <- function(a) {
      if (variance == 0) return(plogis(a))
      integrate(function(z) plogis(a + sqrt(variance) * z) * dnorm(z),
                -Inf, Inf, rel.tol = 1e-9)$value
    }
    bound <- abs(qlogis(prevalence)) + 12 * sqrt(variance) + 40
    intercept <- uniroot(function(a) population_mean(a) - prevalence,
                         c(-bound, bound), tol = 1e-9)$root
  }

  draw_x <- function(m) {
    x <- matrix(rnorm(m * p), nrow = m, ncol = p)
    if (design == "ar1" && p > 1L) {
      for (j in 2:p) x[, j] <- rho * x[, j - 1L] + sqrt(1 - rho^2) * x[, j]
    } else if (design == "equicor") {
      x <- sqrt(1 - rho) * x + sqrt(rho) * rnorm(m)
    } else if (design == "block") {
      for (start in seq.int(1L, p, by = block_size)) {
        idx <- seq.int(start, min(start + block_size - 1L, p))
        x[, idx] <- sqrt(1 - rho) * x[, idx, drop = FALSE] + sqrt(rho) * rnorm(m)
      }
    }
    x
  }
  draw_data <- function(m) {
    x <- draw_x(m)
    eta <- as.vector(intercept + x %*% beta)
    prob <- plogis(eta)
    y <- rbinom(m, size = trials, prob = prob)
    list(x = x, y = y, prob = prob, eta = eta,
         y_pm1 = if (trials == 1L) 2L * y - 1L else NULL,
         response = cbind(success = y, failure = trials - y))
  }
  train <- draw_data(n)
  test <- if (n_test > 0L) draw_data(n_test) else NULL
  list(train = train, test = test, beta = beta, intercept = intercept,
       support = which(beta != 0),
       settings = list(n = n, p = p, design = design, rho = rho,
                       block_size = block_size, trials = trials,
                       signal_sd = sqrt(variance), prevalence = prevalence))
}

bh2011 <- function(n = 400L, p = 1000L, amplitude = 1,
                   n_test = 2000L, seed = NULL) {
  simulate_binomial(n = n, p = p,
                    beta_active = amplitude * c(1, 1, -1, -1),
                    design = "independent", n_test = n_test, seed = seed)
}

simulation_setting <- function(setting = c("independent", "ar1", "block", "imbalanced"),
                               n = 200L, p = 1000L, amplitude = 1,
                               n_test = 2000L, seed = NULL) {
  setting <- match.arg(setting)
  simulate_binomial(n = n, p = p,
    beta_active = amplitude * rep(c(1, -1), 5),
    design = switch(setting, independent = "independent", ar1 = "ar1",
                    block = "block", imbalanced = "ar1"),
    rho = 0.5, block_size = 20L,
    prevalence = if (setting == "imbalanced") 0.2 else 0.5,
    n_test = n_test, seed = seed)
}

# Examples (uncomment to execute):
# dat <- simulation_setting("ar1", seed = 2026)
# x <- dat$train$x
# y <- dat$train$y
# beta_true <- dat$beta
# dat_paper <- bh2011(seed = 2026)
# dat_imb <- simulation_setting("imbalanced", seed = 2026)
# dat_grouped <- simulate_binomial(trials = 10L, seed = 2026)
# For glmnet grouped binomial: use dat_grouped$train$response as y.
# Repeated runs: set seed once outside the loop and omit per-call seed.
# set.seed(2026)
# for (b in seq_len(200L)) {
#   dat <- simulation_setting("ar1", n_test = 2000L)
#   # Fit ALL competing methods to this same dat; retain only metrics.
# }
# Fit standardization using training data only and transform test accordingly.
# Returned beta and intercept always refer to the original X scale.
