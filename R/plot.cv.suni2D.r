#' Plot a two-dimensional cross-validation surface for
#'   a \code{cv.suni2D} object.
#'
#' Show the mean validation loss, plus and minus one standard error, and the
#'   selected penalty pairs. Use a curve when there is only one
#'   second-coordinate value.
#'
#' @details Transformed coordinates are matched to the stored traversal path.
#'   Missing or inconsistent grid metadata produces an error. Zero tau values
#'   are mapped to negative infinity on the log axis, with a warning in the
#'   transformed-coordinate branch.
#'
#' @param x A \code{cv.suni2D} object.
#' @param use_axis \code{"original"} uses the fitted parameterization;
#'   \code{"zeta_tau"} uses transformed coordinates from \code{penalty.path}
#'   if the \code{method} of the fitted object is not "zeta_tau". Default
#'   is \code{"original"}.
#'   Zeta, tau, lambda, sum, and rate axes are logarithmic; alpha is linear on
#'   its original axis.
#' @param ... Additional arguments passed to \code{plotly::plot_ly} in the
#'   surface branch; unused by the single-curve branch.
#'
#' @return For a surface, prints and invisibly returns a plotly widget. For a
#'   curve, draws a base graphics plot and invisibly returns \code{NULL}.
#' @examples
#' set.seed(1)
#' nobs <- as.integer(50)
#' nvars <- as.integer(6)
#' x <- matrix(rnorm(nobs * nvars), nrow = nobs, ncol = nvars)
#' y <- drop(x[, 1] - x[, 2] + rnorm(nobs))
#' cvfit <- cv.suni2D(
#'   x, y, nfolds = 3, method = "lambda_alpha",
#'   nlambda = 8, alpha = c(0.6, 0.75, 0.9)
#' )
#'
#' # Plot using the fitted model's lambda-alpha coordinates (the default).
#' plot(cvfit)
#'
#' # Plot the same CV surface using the equivalent zeta-tau coordinates.
#' plot(cvfit, use_axis = "zeta_tau")
#' @import plotly
#' @importFrom graphics legend lines par
#' @export
plot.cv.suni2D <- function(
  x,
  use_axis = c("original", "zeta_tau"), ...
) {
  cvobj <- x
  method <- cvobj$suni2D.fit$method
  use_axis <- match.arg(use_axis)
  z <- cvobj$cvm
  uppersurf <- cvobj$cvupper
  lowersurf <- cvobj$cvlower
  penalty1 <- cvobj[[1L]]
  penalty2 <- cvobj[[2L]]
  is_surface <- length(penalty2) > 1L

  y.min <- cvobj$cv.min[[2]]
  y.1se <- cvobj$cv.1se[[2]]
  x.min <- cvobj$cv.min[[1]]
  x.1se <- cvobj$cv.1se[[1]]

  if (method == "zeta_tau") {
    xlab <- "log(zeta)"
    ylab <- "log(tau)"
    x <- log(cvobj$zeta)
    y <- log(cvobj$tau)
    y.min <- log(y.min)
    y.1se <- log(y.1se)
    x.min <- log(x.min)
    x.1se <- log(x.1se)
  } else if (use_axis == "original") {
    xlab <- paste0("log(", strsplit(method, "_")[[1]][1], ")")
    ylab <- switch(method,
      lambda_alpha = strsplit(method, "_")[[1]][2],
      sum_rate = paste0("log(", strsplit(method, "_")[[1]][2], ")")
    )
    x <- log(cvobj[[1]])
    y <- switch(method,
      lambda_alpha = cvobj[[2]],
      sum_rate = log(cvobj[[2]])
    )
    x.min <- log(x.min)
    x.1se <- log(x.1se)
    y.min <- switch(method,
      lambda_alpha = y.min,
      sum_rate = log(y.min)
    )
    y.1se <- switch(method,
      lambda_alpha = y.1se,
      sum_rate = log(y.1se)
    )
  } else {
    xlab <- "log(zeta)"
    ylab <- "log(tau)"

    penalty.path <- cvobj$suni2D.fit$penalty.path
    if (is.null(penalty.path) || length(dim(penalty.path)) != 2L ||
      ncol(penalty.path) < 4L) {
      stop(
        "suni2D.fit$penalty.path must be a matrix with at least four columns",
        call. = FALSE
      )
    }
    if (nrow(z) != length(penalty2) || ncol(z) != length(penalty1)) {
      stop("the CV surface dimensions do not match the penalty grid", call. = FALSE)
    }

    ## cvm rows correspond to penalty2 (alpha/rate), while its columns
    ## correspond to penalty1 (lambda/sum).  Zeta and tau both depend on
    ## the two original penalties, so plotly needs a coordinate matrix for
    ## each rather than two independent coordinate vectors.
    x <- matrix(NA_real_, nrow(z), ncol(z), dimnames = dimnames(z))
    y <- matrix(NA_real_, nrow(z), ncol(z), dimnames = dimnames(z))
    match_penalty <- function(values, target) {
      exact <- which(values == target)
      if (length(exact)) {
        return(exact)
      }

      finite <- is.finite(values) & is.finite(target)
      distance <- abs(values - target)
      candidates <- which(finite &
        distance <= sqrt(.Machine$double.eps) *
          pmax(1, abs(values), abs(target)))
      if (!length(candidates)) {
        return(integer(0))
      }

      ## Return every exact duplicate of the nearest grid value, but do not
      ## merge distinct grid values merely because both are within tolerance.
      nearest <- values[candidates[which.min(distance[candidates])]]
      which(values == nearest)
    }

    for (i in seq_len(nrow(penalty.path))) {
      row.id <- match_penalty(penalty2, penalty.path[i, 2L])
      col.id <- match_penalty(penalty1, penalty.path[i, 1L])
      if (!length(row.id) || !length(col.id)) {
        stop("suni2D.fit$penalty.path does not match the CV penalty grid", call. = FALSE)
      }
      x[row.id, col.id] <- penalty.path[i, 3L]
      y[row.id, col.id] <- penalty.path[i, 4L]
    }

    x <- log(x)
    if (any(y == 0, na.rm = TRUE)) {
      warning("tau contains 0, which is regular lasso; plot ignores 0 case.")
    }
    y <- ifelse(y == 0, -Inf, log(y))

    selected_coordinates <- function(selected1, selected2) {
      row.id <- match_penalty(penalty2, selected2)[1L]
      col.id <- match_penalty(penalty1, selected1)[1L]
      coordinates <- c(x[row.id, col.id], y[row.id, col.id])
      if (length(coordinates) != 2L || anyNA(coordinates)) {
        stop(
          "the selected CV penalty pair is missing from suni2D.fit$penalty.path",
          call. = FALSE
        )
      }
      coordinates
    }

    min.coordinates <- selected_coordinates(x.min, y.min)
    one.se.coordinates <- selected_coordinates(x.1se, y.1se)
    x.min <- min.coordinates[1L]
    y.min <- min.coordinates[2L]
    x.1se <- one.se.coordinates[1L]
    y.1se <- one.se.coordinates[2L]
  }

  if (is_surface) {
    fig <- plot_ly(
      x = x,
      y = y,
      z = z,
      type = "surface",
      colorscale = "Viridis",
      name = "CVM",
      ...
    ) %>%
      add_surface(
        x = x,
        y = y,
        z = uppersurf,
        opacity = 0.3,
        showscale = FALSE,
        name = "Upper CI",
        colorscale = list(c(0, 1), c("red", "red"))
      ) %>%
      add_surface(
        x = x,
        y = y,
        z = lowersurf,
        opacity = 0.3,
        showscale = FALSE,
        name = "Lower CI",
        colorscale = list(c(0, 1), c("blue", "blue"))
      ) %>%
      # inherit = FALSE keeps surface-only attributes (e.g. colorscale) set in
      # plot_ly() from propagating to the scatter3d traces
      add_trace(
        type = "scatter3d",
        mode = "lines",
        x = c(x.min, x.min),
        y = c(y.min, y.min),
        z = c(min(lowersurf, na.rm = TRUE), max(uppersurf, na.rm = TRUE)),
        line = list(color = "green", width = 6),
        name = "min line",
        inherit = FALSE
      ) %>%
      add_trace(
        type = "scatter3d",
        mode = "lines",
        x = c(x.1se, x.1se),
        y = c(y.1se, y.1se),
        z = c(min(lowersurf, na.rm = TRUE), max(uppersurf, na.rm = TRUE)),
        line = list(color = "orange", width = 6),
        name = "1se line",
        inherit = FALSE
      ) %>%
      layout(
        scene = list(
          xaxis = list(title = xlab),
          yaxis = list(title = ylab),
          zaxis = list(title = "cvm")
        )
      )

    print(fig)
    invisible(fig)
  } else {
    y_range <- range(c(z, as.vector(uppersurf), as.vector(lowersurf)), na.rm = TRUE)
    plot(as.vector(x), as.vector(z),
      type = "l", lwd = 2, col = "blue",
      xlab = xlab, ylab = "Mean CV Error",
      main = "Cross-Validation Curve", ylim = y_range
    )
    lines(as.vector(x), as.vector(uppersurf), col = "gray", lty = 2)
    lines(as.vector(x), as.vector(lowersurf), col = "gray", lty = 2)
    abline(v = x.min, col = "red", lty = 3, lwd = 2)
    abline(v = x.1se, col = "green", lty = 3, lwd = 2)

    legend("topright",
      legend = c("CV Error", "Upper/Lower CI", "min", "1se"),
      col = c("blue", "gray", "red", "green"), lty = c(1, 2, 3, 3), lwd = 2
    )
    invisible()
  }
}
