##' @export
predict.sunipath_2 <- function(object, newx, s = NULL, type = c("link"),
                               alpha, ...) {
  type <- match.arg(type)
  b0 <- lapply(object$b0, function(x){
    b = t(as.matrix(x))
    rownames(b) = "(Intercept)"
    b
  })
  nbeta <- lapply(seq_along(b0),function(x){
    rbind2(b0[[x]],object$beta[[x]])
  })
  names(nbeta) = names(b0)
  if (!is.null(s)) {
    lambda <- object$lambda
    lamlist <- lambda.interp(lambda, s)

    nbeta = lapply(nbeta, function(x){
      vnames <- dimnames(x)[[1]]
      dimnames(x) <- list(NULL, NULL)

      x <- x[, lamlist$left, drop = FALSE] %*%
        Diagonal(x = lamlist$frac) +
        x[, lamlist$right, drop = FALSE] %*%
        Diagonal(x = 1-lamlist$frac)
      dimnames(x) <- list(vnames, paste(seq(along = s)))
      x
    })
  }

  if(!missing(alpha)){
    nbeta_2 = list()
    if(object$use_alpha){
      for(a in alpha){
        if(a %in% object$alpha){
          nbeta_2[[paste0("alpha_",a)]] = nbeta[[paste0("alpha_",a)]]
        }
      }
    }else{
      for(a in alpha){
        if(a %in% object$lamPos){
          nbeta_2[[paste0("lamPos_",a)]] = nbeta[[paste0("lamPos_",a)]]
        }
      }
    }
    if(length(nbeta_2) != 0){
      nbeta = nbeta_2
    }
  }

  nnewx = as.matrix(cbind2(1, newx))
  nfit <- lapply(seq_along(nbeta), function(x){as.matrix(nnewx %*% nbeta[[x]])})
  names(nfit) = names(nbeta)
  nfit
}
