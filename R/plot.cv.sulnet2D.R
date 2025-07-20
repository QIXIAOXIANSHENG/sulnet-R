########################################################################
## This function is adapted/modified based on the plot.cv function
## from the `glmnet` package:
##
## Jerome Friedman, Trevor Hastie, Robert Tibshirani (2010).
## Regularization Paths for Generalized Linear Models via Coordinate Descent.
##   Journal of Statistical Software, 33(1), 1-22.
##   URL: https://www.jstatsoft.org/v33/i01/.

##' Plot the cross-validation curve produced by cv.sulnet
##'
##' Plots the cross-validation surface, and upper and lower standard deviation
##' surface.
##'
##' A plotly object is returned.
##'
##' @param x fitted \code{\link{cv.sulnet2D}} object
##' @param \dots other graphical parameters to plotly
##' @author Yi Yang, Yuwen Gu and Hui Zou\cr
##'
##' Maintainer: Yi Yang <yi.yang6@mcgill.ca>
##' @seealso \code{\link{cv.sulnet2D}}.
##' @references Yang, Y. and Zou, H. (2012).
##'   "An Efficient Algorithm for Computing The HHSVM and Its Generalizations."
##'   \emph{Journal of Computational and Graphical Statistics}, 22, 396-415.\cr
##'   BugReport: \url{https://github.com/emeryyi/sulnet}\cr
##'
##'   Gu, Y., and Zou, H. (2016).
##'   "High-dimensional generalizations of asymmetric least squares regression and their applications."
##'   \emph{The Annals of Statistics}, 44(6), 2661–2694.\cr
##'
##'   Friedman, J., Hastie, T., and Tibshirani, R. (2010).
##'   "Regularization paths for generalized linear models via coordinate descent."
##'   \emph{Journal of Statistical Software, 33, 1.}\cr
##'   \url{https://www.jstatsoft.org/v33/i01/}
##' @keywords models regression
##' @examples
##'
##' # fit an elastic net penalized logistic regression with lambda2 = 1 for the
##' # L2 penalty. Use the logistic loss as the cross validation prediction loss.
##' # Use five-fold CV to choose the optimal lambda for the L1 penalty.
##' data(FHT)
##' set.seed(2011)
##' cv=cv.sulnet2D(FHT$x, FHT$y)
##' plot(cv)
##'
##' @import plotly
##' @export
plot.cv.sulnet2D <- function(x,  ...) {
  cvobj <- x
  xlab <- "log(Lambda)"
  x = log(cvobj$lambda)
  if(cvobj$sulnet2D.fit$use_alpha){
    y = cvobj$alpha
    ylab = "alpha"
    y.min = cvobj$cv.min$alpha.min
    y.1se = cvobj$cv.1se$alpha.1se
  }else{
    y = cvobj$lamPos
    ylab = "lamPos"
    y.min = cvobj$cv.min$lamPos.min
    y.1se = cvobj$cv.1se$lamPos.1se
  }
  z = cvobj$cvm

  if(length(y)>1){
    fig <- suppressWarnings({

      fig <- plot_ly(
        x = x,
        y = y,
        z = z,
        type = "surface",
        colorscale = "Viridis",
        name = "CVM",
        ...
      )%>%
        add_surface(
          z = cvobj$cvupper,
          opacity = 0.3,
          showscale = FALSE,
          name = "Upper CI",
          colorscale = list(c(0,1), c("red","red"))
        ) %>%
        add_surface(
          z = cvobj$cvlower,
          opacity = 0.3,
          showscale = FALSE,
          name = "Lower CI",
          colorscale = list(c(0,1), c("blue","blue"))
        )%>%
        add_trace(
          type = 'scatter3d',
          mode = 'lines',
          x = log(c(cvobj$cv.min$lambda.min, cvobj$cv.min$lambda.min)),
          y = c(y.min, y.min),
          z = c(min(cvobj$cvlower), max(cvobj$cvupper)),
          line = list(color = 'green', width = 6),
          name = "min line"
        )%>%
        add_trace(
          type = 'scatter3d',
          mode = 'lines',
          x = log(c(cvobj$cv.1se$lambda.1se, cvobj$cv.1se$lambda.1se)),
          y = c(y.1se, y.1se),
          z = c(min(cvobj$cvlower), max(cvobj$cvupper)),
          line = list(color = 'orange', width = 6),
          name = "1se line"
        ) %>%
        layout(
          scene = list(
            xaxis = list(title = xlab),
            yaxis = list(title = ylab),
            zaxis = list(title = "cvm")
          )
        )

    })

    suppressWarnings(print(fig))
    invisible(fig)
  }else{
    y_range <- range(c(z, as.vector(cvobj$cvupper), as.vector(cvobj$cvlower)), na.rm = TRUE)
    plot(x, as.vector(z), type = "l", lwd = 2, col = "blue",
         xlab = "log(Lambda)", ylab = "Mean CV Error",
         main = "Cross-Validation Curve", ylim = y_range)
    lines(x, as.vector(cvobj$cvupper), col = "gray", lty = 2)
    lines(x, as.vector(cvobj$cvlower), col = "gray", lty = 2)
    abline(v = log(cvobj$cv.min$lambda.min), col = "red", lty = 3, lwd = 2)
    abline(v = log(cvobj$cv.1se$lambda.1se), col = "green", lty = 3, lwd = 2)

    legend("topright", legend = c("CV Error", "Upper/Lower CI", "Lambda.min", "Lambda.1se"),
           col = c("blue", "gray", "red", "green"), lty = c(1, 2, 3, 3), lwd = 2)
    invisible()
  }



}
