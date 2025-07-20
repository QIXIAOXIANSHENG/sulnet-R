##' Fits the regularization paths for large margin classifiers
##'
##' Fits a regularization path for large margin classifiers at a sequence of
##' regularization parameters lambda.
##'
##' @param x matrix of predictors, of dimension \eqn{N \times p}{N*p}; each row
##'   is an observation vector.
##' @param y response variable. This argument should be a two-level factor for
##'   classification.
##' @param nlambda the number of \code{lambda} values - default is 100.
##' @param method a character string specifying the loss function to use, valid
##'   options are: \itemize{ \item \code{"ls"} least square loss. \item \code{"suni"}
##'   soft unilasso loss. } Default is \code{"ls"}.
##' @param lambda.factor The factor for getting the minimal lambda in
##'   \code{lambda} sequence, where \code{min(lambda)} = \code{lambda.factor} *
##'   \code{max(lambda)}, where \code{max(lambda)} is the smallest value of
##'   \code{lambda} for which all coefficients are zero. The default depends on
##'   the relationship between \eqn{N} (the number of rows in the matrix of
##'   predictors) and \eqn{p}{p} (the number of predictors). If \eqn{N > p}, the
##'   default is \code{0.0001}, close to zero. If \eqn{N<p}, the default is
##'   \code{0.01}. A very small value of \code{lambda.factor} will lead to a
##'   saturated fit. It takes no effect if there is user-defined \code{lambda}
##'   sequence.
##' @param lambda a user supplied \code{lambda} sequence. Typically, by leaving
##'   this option unspecified users can have the program compute its own
##'   \code{lambda} sequence based on \code{nlambda} and \code{lambda.factor}.
##'   Supplying a value of \code{lambda} overrides this. It is better to supply
##'   a decreasing sequence of \code{lambda} values than a single (small) value,
##'   if not, the program will sort user-defined \code{lambda} sequence in
##'   decreasing order automatically.
##' @param lambda2 regularization parameter \eqn{\lambda_2}{lambda2} for the
##'   quadratic penalty of the coefficients.
##' @param pf L1 penalty factor of length \eqn{p}{p} used for adaptive LASSO or
##'   adaptive elastic net. Separate L1 penalty weights can be applied to each
##'   coefficient of \eqn{\beta}{beta} to allow differential L1 shrinkage. Can
##'   be 0 for some variables, which implies no L1 shrinkage, and results in
##'   that variable always being included in the model. Default is 1 for all
##'   variables (and implicitly infinity for variables listed in
##'   \code{exclude}).
##' @param pf2 L2 penalty factor of length \eqn{p}{p} used for adaptive LASSO or
##'   adaptive elastic net. Separate L2 penalty weights can be applied to each
##'   coefficient of \eqn{\beta}{beta} to allow differential L2 shrinkage. Can
##'   be 0 for some variables, which implies no L2 shrinkage. Default is 1 for
##'   all variables.
##' @param exclude indices of variables to be excluded from the model. Default
##'   is none. Equivalent to an infinite penalty factor.
##' @param dfmax limit the maximum number of variables in the model. Useful for
##'   very large \eqn{p}, if a partial path is desired. Default is \eqn{p+1}.
##' @param pmax limit the maximum number of variables ever to be nonzero. For
##'   example once \eqn{\beta} enters the model, no matter how many times it
##'   exits or re-enters model through the path, it will be counted only once.
##'   Default is \code{min(dfmax*1.2,p)}.
##' @param standardize logical flag for variable standardization, prior to
##'   fitting the model sequence. If \code{TRUE}, \code{x} matrix is normalized
##'   such that \code{x} is centered (i.e.
##'   \eqn{\sum^N_{i=1}x_{ij}=0}{sum(Xj)=0}), and sum squares of each column
##'   \eqn{\sum^N_{i=1}x_{ij}^2/N=1}{<Xj,Xj>/N=1}. If \code{x} matrix is
##'   standardized, the ending coefficients will be transformed back to the
##'   original scale. Default is \code{FALSE}.
##' @param intercept logical flag to indicate whether to include or exclude the
##'   intercept in the model.
##' @param eps convergence threshold for coordinate majorization descent. Each
##'   inner coordinate majorization descent loop continues until the relative
##'   change in any coefficient (i.e., \eqn{\max_j(\beta_j^{new}
##'   -\beta_j^{old})^2}{max(j)(beta_new[j]-beta_old[j])^2}) is less than
##'   \code{eps}. For HHSVM, i.e., \code{method="hhsvm"}, it is
##'   \eqn{\frac{2}{\delta}\max_j(\beta_j^{new}-\beta_j^{old})^2}{2*max(j)
##'   (beta_new[j]-beta_old[j])^2/delta}. For expectile regression,
##'   i.e., \code{method="er"}, it is \eqn{2\max(1-\omega,\omega)
##'   \max_j(\beta_j^{new}-\beta_j^{old})^2}{2*max(1-omega,omega)*max(j)
##'   (beta_new[j]-beta_old[j])^2}. Defaults value is \code{1e-8}.
##' @param maxit maximum number of outer-loop iterations allowed at fixed lambda
##'   value. Default is 1e6. If models do not converge, consider increasing
##'   \code{maxit}.
##' @param lamPos the penalty term for the negativity constraint in soft unilasso
##'   model. Default is 0.1.
##' @param loo the logical variable determining if the leave-one-out fits are used
##'   in the soft unilasso model. Default is \code{TRUE}.
##' @param negOnly the logical variable determining if only the second step in the
##'   soft unilasso model is fitted. Default is \code{FALSE}.
##' @param alpha the proportion of penalty term split between lasso penalty \code{lambda}
##'   and soft-unilasso penalty \code{lamPos}. If \code{alpha} is not \code{NULL}, then
##'   \code{alpha} will be used to generate different \code{lamPos}. Otherwise,
##'   \code{lamPos} will be used. Default is \code{seq(0, 0.5, length.out = 100)}
##' @return An object with S3 class \code{\link{sulnet}}. \item{call}{the call
##'   that produced this object} \item{b0}{intercept sequence of length
##'   \code{length(lambda)}} \item{beta}{a \code{p*length(lambda)} matrix of
##'   coefficients, stored as a sparse matrix (\code{dgCMatrix} class, the
##'   standard class for sparse numeric matrices in the \code{Matrix} package.).
##'   To convert it into normal type matrix use \code{as.matrix()}.}
##'   \item{lambda}{the actual sequence of \code{lambda} values used}
##'   \item{df}{the number of nonzero coefficients for each value of
##'   \code{lambda}.} \item{dim}{dimension of coefficient matrix (ices)}
##'   \item{npasses}{total number of iterations (the most inner loop) summed
##'   over all lambda values} \item{jerr}{error flag, for warnings and errors, 0
##'   if no error.}
##'
##' @details Note that the objective function in \code{sulnet} is \deqn{Loss(y,
##'   X, \beta)/N + \lambda_1\Vert\beta\Vert_1 +
##'   0.5\lambda_2\Vert\beta\Vert_2^2}{Loss(y, X, beta))/N + lambda1 * |beta| +
##'   0.5 * lambda2 * |beta|^2} where the penalty is a combination of L1 and L2
##'   term. Users can specify the loss function to use, options include
##'   Huberized squared hinge loss, Squared hinge loss, least square loss,
##'   logistic regression and expectile regression loss. Users can also tweak
##'   the penalty by choosing different \eqn{lambda2} and penalty factor.
##'
##'   For computing speed reason, if models are not converging or running slow,
##'   consider increasing \code{eps}, decreasing \code{nlambda}, or increasing
##'   \code{lambda.factor} before increasing \code{maxit}.
##'
##'   \strong{FAQ:}
##'
##'   \bold{Question: }\dQuote{\emph{I couldn't get an idea how to specify an
##'   option to get adaptive LASSO, how to specify an option to get elastic net
##'   and adaptive elastic net? Could you please give me a quick hint?}}
##'
##'   \bold{Answer: } \code{lambda2} is the regularize parameter for L2 penalty
##'   part. To use LASSO, set \code{lambda2=0}. To use elastic net, set
##'   \code{lambda2} as nonzero.
##'
##'   \code{pf} is the L1 penalty factor of length \eqn{p}{p} (\eqn{p}{p} is the
##'   number of predictors). Separate L1 penalty weights can be applied to each
##'   coefficient to allow differential L1 shrinkage. Similiarly \code{pf2} is
##'   the L2 penalty factor of length \eqn{p}{p}.
##'
##'   To use adaptive LASSO, you should set \code{lambda2=0} and also specify
##'   \code{pf} and \code{pf2}. To use adaptive elastic net, you should set
##'   \code{lambda2} as nonzero and specify \code{pf} and \code{pf2},
##'
##'   For example:
##'
##'   \preformatted{
##'     library('sulnet')
##'     # Dataset N = 100, p = 10
##'     x_log <- matrix(rnorm(100*10),100,10)
##'     y_log <- sample(c(-1,1),100,replace=TRUE)
##'
##'     # LASSO
##'     m <- sulnet(x=x_log,y=y_log,lambda2=0,method="log")
##'     plot(m)
##'
##'     # elastic net with lambda2 = 1
##'     m <- sulnet(x=x_log,y=y_log,lambda2=1,method="log")
##'     plot(m)
##'
##'     # adaptive lasso with penalty factor
##'     # pf = 0.5 0.5 0.5 0.5 0.5 1.0 1.0 1.0 1.0 1.0
##'     m <- sulnet(x=x_log,y=y_log,lambda2=0,method="log",
##'                 pf=c(rep(0.5,5),rep(1,5)))
##'     plot(m)
##'
##'     # adaptive elastic net with lambda2 = 1 and penalty factor pf =
##'     # c(rep(0.5,5),rep(1,5)) pf2 = 3 3 3 3 3 1 1 1 1 1
##'     m <- sulnet(x=x_log,y=y_log,lambda2=1,method="log",
##'                 pf=c(rep(0.5,5),rep(1,5)),
##'                 pf2 = c(rep(3,5),rep(1,5)))
##'     plot(m)
##'   }
##'   \bold{Question: }\dQuote{\emph{what is the meaning of the parameter
##'   \code{pf}? On the package documentation, it said \code{pf} is the penalty
##'   weight applied to each coefficient of beta?}}
##'
##'   \bold{Answer: } Yes, \code{pf} and \code{pf2} are L1 and L2 penalty factor
##'   of length \eqn{p}{p} used for adaptive LASSO or adaptive elastic net. 0
##'   means that the feature (variable) is always excluded, 1 means that the
##'   feature (variable) is included with weight 1.
##'
##'   \bold{Question: }\dQuote{\emph{Does sulnet deal with both continuous and
##'   categorical response variables?}}
##'
##'   \bold{Answer: } Yes, both are supported, you can use a continuous type
##'   response variable with the least squares regression loss, or a categorical
##'   type response with losses for classification problem.
##'
##'   \bold{Question: }\dQuote{\emph{Why does predict function not work? predict
##'   should return the predicted probability of the positive class. Instead I
##'   get:}}
##'   \preformatted{
##'     Error in as.matrix(as.matrix(cbind2(1, newx)) \%*\% nbeta):
##'     error in evaluating the argument 'x' in selecting a method for function 'as.matrix':
##'     Error in t(.Call(Csparse_dense_crossprod, y, t(x))):
##'     error in evaluating the argument 'x' in selecting a method for function 't':
##'     Error: Cholmod error 'X and/or Y have wrong dimensions' at
##'       file ../MatrixOps/cholmod_sdmult.c, line 90?
##'   }
##'
##'   \dQuote{\emph{Using the Arcene dataset and executing the following code
##'   will give the above error:}}
##'   \preformatted{
##'     library(sulnet)
##'     arc <- read.csv("arcene.csv", header=FALSE)
##'     fit <- sulnet(arc[,-10001], arc[,10001], standardize=FALSE,
##'                   method="logit")
##'     pred <- rnorm(10000)
##'     predict(fit, pred, type="link")
##'   }
##'
##'   \bold{Answer: } It is actually NOT a bug of sulnet. When make prediction
##'   using a new matrix x, each observation of x should be arranged as a row of
##'   a matrix. In your code, because "pred" is a vector, you need to convert
##'   "pred" into a matrix, try the following code:
##'   \preformatted{
##'      pred <- rnorm(10000)
##'      pred <- matrix(pred,1,10000)
##'      predict(fit, pred, type="link")
##'   }
##'
##' @author Yi Yang, Yuwen Gu and Hui Zou\cr
##'   Maintainer: Yi Yang <yi.yang6@mcgill.ca>
##'
##' @seealso \code{plot.sulnet}
##'
##' @references Yang, Y. and Zou, H. (2012).
##'   "An Efficient Algorithm for Computing The HHSVM and Its Generalizations."
##'   \emph{Journal of Computational and Graphical Statistics}, 22, 396-415.\cr
##'   BugReport: \url{https://github.com/emeryyi/sulnet}\cr
##'
##'   Gu, Y., and Zou, H. (2016).
##'   "High-dimensional generalizations of asymmetric least squares regression and their applications."
##'   \emph{The Annals of Statistics}, 44(6), 2661–2694.\cr
##'
##' @keywords models regression
##'
##' @examples
##'
##' data(FHT)
##' # 1. solution paths for the LASSO penalized least squares.
##' # To use LASSO set lambda2 = 0.
##'
##' m1 <- sulnet(x = FHT$x, y = FHT$y_reg, lambda2 = 0, method = "ls")
##' plot(m1)
##'
##'
##' @export
##'
##' @useDynLib sulnet, .registration = TRUE
##'
sulnet <- function(x, y, nlambda = 100,
                   method = c("ls", "suni"),
                   lambda.factor = ifelse(nobs < nvars, 0.01, 1e-04),
                   lambda = NULL, lambda2 = 0, pf = rep(1, nvars),
                   pf2 = rep(1, nvars), exclude, dfmax = nvars + 1,
                   pmax = min(dfmax * 1.2, nvars), standardize = FALSE,
                   intercept = TRUE, eps = 1e-08, maxit = 1e+06, lamPos = 0.1,
                   loo = TRUE, alpha = 0.2, negOnly = FALSE) {
################################################################################
  ## data setup
  method <- match.arg(method)
  this.call <- match.call()
  y <- drop(y)
  x <- as.matrix(x)
  np <- dim(x)
  nobs <- as.integer(np[1])
  nvars <- as.integer(np[2])
  vnames <- colnames(x)
  if (is.null(vnames))
    vnames <- paste("V", seq(nvars), sep = "")
  if (NROW(y) != nobs)
    stop("x and y have different number of observations")
  if (NCOL(y) > 1L) stop("Multivariate response is not supported now")
################################################################################
  ## parameter setup
  if (length(pf) != nvars)
    stop("The size of L1 penalty factor must be same as the number of input variables")
  if (length(pf2) != nvars)
    stop("The size of L2 penalty factor must be same as the number of input variables")
  if (lambda2 < 0)
    stop("lambda2 must be non-negative")
  maxit <- as.integer(maxit)
  lam2 <- as.double(lambda2)
  pf <- as.double(pf)
  pf2 <- as.double(pf2)
  isd <- as.integer(standardize)
  intr <- as.integer(intercept)
  eps <- as.double(eps)
  dfmax <- as.integer(dfmax)
  pmax <- as.integer(pmax)
  lamPos <- as.double(lamPos)
  ignore_lamPos = FALSE
  if(!is.null(alpha)){
    alpha <- as.double(alpha)
    ignore_lamPos = TRUE
  }
  if (!missing(exclude)) {
    jd <- match(exclude, seq(nvars), 0)
    if (!all(jd > 0))
      stop("Some excluded variables out of range")
    jd <- as.integer(c(length(jd), jd))
  } else jd <- as.integer(0)
################################################################################
  ## lambda setup
  nlam <- as.integer(nlambda)
  if (is.null(lambda)) {
    if (lambda.factor >= 1)
      stop("lambda.factor should be less than 1")
    flmin <- as.double(lambda.factor)
    ulam <- double(1)
  } else {
    ## flmin=1 if user define lambda
    flmin <- as.double(1)
    if (any(lambda < 0))
      stop("lambdas should be non-negative")
    ulam <- as.double(rev(sort(lambda)))
    nlam <- as.integer(length(lambda))
  }
################################################################################
  fit <- switch(method,
                suni = sunipath(x, y, nlam, flmin, ulam, isd, intr, eps, dfmax,
                                  pmax, jd, pf, pf2, maxit, lam2,lamPos, loo, negOnly,
                                  nobs, nvars, vnames, alpha, ignore_lamPos),
                ls = lspath(x, y, nlam, flmin, ulam, isd, intr, eps, dfmax,
                            pmax, jd, pf, pf2, maxit, lam2, nobs, nvars, vnames))
  if (is.null(lambda))
    fit$lambda <- lamfix(fit$lambda)
  fit$call <- this.call
################################################################################
  class(fit) <- c("sulnet", class(fit))
  fit
}
