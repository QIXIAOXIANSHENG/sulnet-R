##' @importFrom methods rbind2
##' @import Matrix
##' @export
coef.sunipath <- function(object, penalty1 = NULL, penalty2 = NULL, type = c("coefficients", "nonzero", "inner"), ...) {
  type <- match.arg(type)

  using <- object$use
  pen <- switch(using,
    zeta_tau = pen_zt(object, penalty1, penalty2),
    lambda_alpha = pen_la(object, penalty1, penalty2),
    sum_rate = pen_sr(object, penalty1, penalty2)
  )

  if (type == "inner") {
    theta0 <- Diagonal(x = pen$pen1$frac) %*%
      object$inner$theta0[pen$pen1$left, pen$pen2$left, drop = FALSE] %*%
      Diagonal(x = pen$pen2$frac) +
      Diagonal(x = pen$pen1$frac) %*%
      object$inner$theta0[pen$pen1$left, pen$pen2$right, drop = FALSE] %*%
      Diagonal(x = 1 - pen$pen2$frac) +
      Diagonal(x = 1 - pen$pen1$frac) %*%
      object$inner$theta0[pen$pen1$right, pen$pen2$left, drop = FALSE] %*%
      Diagonal(x = pen$pen2$frac) +
      Diagonal(x = 1 - pen$pen1$frac) %*%
      object$inner$theta0[pen$pen1$right, pen$pen2$right, drop = FALSE] %*%
      Diagonal(x = 1 - pen$pen2$frac)
    nz <- length(pen$pen1$frac)
    W <- Matrix::sparseMatrix(
      i = c(pen$pen1$left, pen$pen1$right),
      j = rep(seq_len(nz), 2L),
      x = c(pen$pen1$frac, 1 - pen$pen1$frac),
      dims = c(ncol(object$inner$theta[[1L]]), nz)
    )
    ## Interpolate coefficients and combine with intercept
    theta <- lapply(seq_along(pen$pen2$frac), function(j) {
      t1 <- pen$pen2$left[j]
      t2 <- pen$pen2$right[j]
      b <- pen$pen2$frac[j]
      g <- (b * object$inner$theta[[t1]] +
        (1 - b) * object$inner$theta[[t2]]) %*% W
      rbind(theta0[, j], g)
    })
    names(theta) <- paste0(
      strsplit(using, "_")[[1]][2],
      seq_along(theta)
    )
    return(theta)
  } else {
    gamma0 <- Diagonal(x = pen$pen1$frac) %*%
      object$gamma0[pen$pen1$left, pen$pen2$left, drop = FALSE] %*%
      Diagonal(x = pen$pen2$frac) +
      Diagonal(x = pen$pen1$frac) %*%
      object$gamma0[pen$pen1$left, pen$pen2$right, drop = FALSE] %*%
      Diagonal(x = 1 - pen$pen2$frac) +
      Diagonal(x = 1 - pen$pen1$frac) %*%
      object$gamma0[pen$pen1$right, pen$pen2$left, drop = FALSE] %*%
      Diagonal(x = pen$pen2$frac) +
      Diagonal(x = 1 - pen$pen1$frac) %*%
      object$gamma0[pen$pen1$right, pen$pen2$right, drop = FALSE] %*%
      Diagonal(x = 1 - pen$pen2$frac)

    nz <- length(pen$pen1$frac)
    W <- Matrix::sparseMatrix(
      i = c(pen$pen1$left, pen$pen1$right),
      j = rep(seq_len(nz), 2L),
      x = c(pen$pen1$frac, 1 - pen$pen1$frac),
      dims = c(ncol(object$gamma[[1L]]), nz)
    )
    ## Interpolate coefficients and combine with intercept
    gamma <- lapply(seq_along(pen$pen2$frac), function(j) {
      t1 <- pen$pen2$left[j]
      t2 <- pen$pen2$right[j]
      b <- pen$pen2$frac[j]
      g <- (b * object$gamma[[t1]] +
        (1 - b) * object$gamma[[t2]]) %*% W
      g <- rbind(gamma0[, j], g)
      colnames(g) <- paste0(strsplit(using, "_")[[1]][1], seq_len(ncol(g)))
      g
    })
    names(gamma) <- paste0(
      strsplit(using, "_")[[1]][2],
      seq_along(gamma)
    )
    if (type == "coefficients") {
      return(gamma)
    }
    if (type == "nonzero") {
      return(lapply(gamma, function(x) {
        nonzero(x[-1, , drop = FALSE], bystep = TRUE)
      }))
    }
  }
}

pen_zt <- function(object, penalty1, penalty2) {
  if (!is.null(penalty1)) {
    pen1 <- lambda.interp(object$zeta, penalty1)
  } else {
    pen1 <- list(
      left = seq_along(object$zeta),
      frac = rep(1, length(object$zeta)),
      right = seq_along(object$zeta)
    )
  }
  if (!is.null(penalty2)) {
    pen2 <- lambda.interp(object$tau, penalty2)
  } else {
    pen2 <- list(
      left = seq_along(object$tau),
      frac = rep(1, length(object$tau)),
      right = seq_along(object$tau)
    )
  }
  return(list(pen1 = pen1, pen2 = pen2))
}
pen_la <- function(object, penalty1, penalty2) {
  if (!is.null(penalty1)) {
    pen1 <- lambda.interp(object$lambda, penalty1)
  } else {
    pen1 <- list(
      left = seq_along(object$lambda),
      frac = rep(1, length(object$lambda)),
      right = seq_along(object$lambda)
    )
  }
  if (!is.null(penalty2)) {
    pen2 <- lambda.interp(object$alpha, penalty2)
  } else {
    pen2 <- list(
      left = seq_along(object$alpha),
      frac = rep(1, length(object$alpha)),
      right = seq_along(object$alpha)
    )
  }
  return(list(pen1 = pen1, pen2 = pen2))
}
pen_sr <- function(object, penalty1, penalty2) {
  if (!is.null(penalty1)) {
    pen1 <- lambda.interp(object$sum, penalty1)
  } else {
    pen1 <- list(
      left = seq_along(object$sum),
      frac = rep(1, length(object$sum)),
      right = seq_along(object$sum)
    )
  }
  if (!is.null(penalty2)) {
    pen2 <- lambda.interp(object$rate, penalty2)
  } else {
    pen2 <- list(
      left = seq_along(object$rate),
      frac = rep(1, length(object$rate)),
      right = seq_along(object$rate)
    )
  }
  return(list(pen1 = pen1, pen2 = pen2))
}
