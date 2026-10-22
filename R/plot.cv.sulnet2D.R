#######################################################################
# This function is adapted/modified based on the plot.cv function
# from the `gcdnet` package:

#' Plot a two-dimensional cross-validation surface for
#'   a \code{cv.sulnet2D} object.
#'
#' Show the mean validation loss, plus and minus one standard error, and the
#'   selected penalty pairs. Use a curve when there is only one
#'   second-coordinate value.
#'
#' @details The horizontal axis is log lambda and the second axis is alpha or
#'   lamPos.
#'
#' @param x A \code{cv.sulnet2D} object.
#' @param ... Additional arguments passed to \code{plotly::plot_ly} in the
#'   surface branch; unused by the single-curve branch.
#'
#' @return For a surface, prints and invisibly returns a plotly widget. For a
#'   curve, draws a base graphics plot and invisibly returns \code{NULL}.
#' @import plotly
#' @importFrom graphics legend lines par
#' @export
plot.cv.sulnet2D <- function(x, ...) {
  cvobj <- x
  xlab <- "log(Lambda)"
  x <- log(cvobj$lambda)
  if (cvobj$sulnet2D.fit$use_alpha) {
    y <- cvobj$alpha
    ylab <- "alpha"
    y.min <- cvobj$cv.min$alpha
    y.1se <- cvobj$cv.1se$alpha
  } else {
    y <- cvobj$lamPos
    ylab <- "lamPos"
    y.min <- cvobj$cv.min$lamPos
    y.1se <- cvobj$cv.1se$lamPos
  }
  z <- cvobj$cvm

  if (length(y) > 1) {
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
        z = cvobj$cvupper,
        opacity = 0.3,
        showscale = FALSE,
        name = "Upper CI",
        colorscale = list(c(0, 1), c("red", "red"))
      ) %>%
      add_surface(
        z = cvobj$cvlower,
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
        x = log(c(cvobj$cv.min$lambda, cvobj$cv.min$lambda)),
        y = c(y.min, y.min),
        z = c(min(cvobj$cvlower, na.rm = TRUE), max(cvobj$cvupper, na.rm = TRUE)),
        line = list(color = "green", width = 6),
        name = "min line",
        inherit = FALSE
      ) %>%
      add_trace(
        type = "scatter3d",
        mode = "lines",
        x = log(c(cvobj$cv.1se$lambda, cvobj$cv.1se$lambda)),
        y = c(y.1se, y.1se),
        z = c(min(cvobj$cvlower, na.rm = TRUE), max(cvobj$cvupper, na.rm = TRUE)),
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
    y_range <- range(c(z, as.vector(cvobj$cvupper), as.vector(cvobj$cvlower)), na.rm = TRUE)
    plot(x, as.vector(z),
      type = "l", lwd = 2, col = "blue",
      xlab = "log(Lambda)", ylab = "Mean CV Error",
      main = "Cross-Validation Curve", ylim = y_range
    )
    lines(x, as.vector(cvobj$cvupper), col = "gray", lty = 2)
    lines(x, as.vector(cvobj$cvlower), col = "gray", lty = 2)
    abline(v = log(cvobj$cv.min$lambda), col = "red", lty = 3, lwd = 2)
    abline(v = log(cvobj$cv.1se$lambda), col = "green", lty = 3, lwd = 2)

    legend("topright",
      legend = c("CV Error", "Upper/Lower CI", "Lambda.min", "Lambda.1se"),
      col = c("blue", "gray", "red", "green"), lty = c(1, 2, 3, 3), lwd = 2
    )
    invisible()
  }
}
