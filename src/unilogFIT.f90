! Fit a univariate logistic regression with labels y in {-1, 1} by minimizing
! (1/nobs)*SUM(LOG(1 + EXP(-y*(beta0 + beta*x)))) with Newton's method and
! step halving. x is standardized internally; beta0 and beta are returned on
! the original scale. prob(i) is the fitted P(y(i) = 1). npass counts Newton
! iterations; convergence is declared when the full Newton step on the
! standardized scale is below eps.
!    jerr = 0    => converged
!    jerr = 7777 => x is constant; beta = 0, beta0 is the intercept-only fit,
!                   and prob is the fraction of y equal to 1
!    jerr = 8001 => some y(i) is not -1 or 1; prob = 0
!    jerr = 8002 => y has a single class; the MLE does not exist,
!                   beta0 = beta = 0, and prob is its limit (all 0 or all 1)
!    jerr = -1   => not converged in maxit iterations, singular Hessian, or
!                   failed step halving (typically complete separation);
!                   the last iterate and its prob are returned
SUBROUTINE unilogFIT(nobs, x, y, xsd, beta0, beta, maxit, eps, npass, jerr)
    ! -----------------------------------------------------------------------
    IMPLICIT NONE
    ! -------- INPUT VARIABLES -------- !
    INTEGER, INTENT(IN) :: nobs
    INTEGER, INTENT(IN) :: maxit
    DOUBLE PRECISION, INTENT(IN) :: x(nobs)
    DOUBLE PRECISION, INTENT(IN) :: y(nobs)
    DOUBLE PRECISION, INTENT(IN) :: eps
    ! -------- OUTPUT VARIABLES -------- !
    INTEGER, INTENT(OUT) :: npass
    INTEGER, INTENT(OUT) :: jerr
    DOUBLE PRECISION, INTENT(OUT) :: beta0
    DOUBLE PRECISION, INTENT(OUT) :: beta
    DOUBLE PRECISION, INTENT(INOUT) :: xsd
    ! -------- LOCAL DECLARATIONS -------- !
    INTEGER :: k
    INTEGER :: npos
    LOGICAL :: conv
    DOUBLE PRECISION :: xmean
    DOUBLE PRECISION :: b0
    DOUBLE PRECISION :: b
    DOUBLE PRECISION :: g0
    DOUBLE PRECISION :: g1
    DOUBLE PRECISION :: h00
    DOUBLE PRECISION :: h01
    DOUBLE PRECISION :: h11
    DOUBLE PRECISION :: det
    DOUBLE PRECISION :: d0
    DOUBLE PRECISION :: d1
    DOUBLE PRECISION :: step
    DOUBLE PRECISION :: loss
    DOUBLE PRECISION :: lossnew
    DOUBLE PRECISION :: xs(nobs)
    DOUBLE PRECISION :: eta(nobs)
    DOUBLE PRECISION :: etanew(nobs)
    DOUBLE PRECISION :: e(nobs)
    DOUBLE PRECISION :: p(nobs)
    DOUBLE PRECISION :: w(nobs)
    DOUBLE PRECISION :: prob(nobs)

    npass = 0
    jerr = 0
    beta0 = 0.0D0
    beta = 0.0D0
    prob = 0.0D0

    ! -------- CHECK RESPONSE -------- !
    IF (ANY(y /= 1.0D0 .AND. y /= -1.0D0)) THEN
        jerr = 8001
        RETURN
    END IF
    npos = COUNT(y > 0.0D0)
    IF (npos == 0 .OR. npos == nobs) THEN
        ! prob = DBLE(npos) / DBLE(nobs)
        jerr = 8002
        RETURN
    END IF

    ! -------- INTERCEPT-ONLY START -------- !
    b0 = LOG(DBLE(npos) / DBLE(nobs - npos))
    b = 0.0D0
    IF (ALL(x == x(1))) THEN
        beta0 = b0
        ! prob = DBLE(npos) / DBLE(nobs)
        jerr = 7777
        RETURN
    END IF

    ! -------- STANDARDIZE x -------- !
    xmean = SUM(x) / DBLE(nobs)
    xs = x - xmean
    xsd = SQRT(SUM(xs * xs) / DBLE(nobs))
    xs = xs / xsd

    eta = b0 + b * xs
    loss = SUM(MAX(-y * eta, 0.0D0) + LOG(1.0D0 + EXP(-ABS(y * eta)))) / DBLE(nobs)

    ! -------- NEWTON ITERATIONS -------- !
    conv = .FALSE.
    DO WHILE (npass < maxit)
        npass = npass + 1
        ! p = 1/(1 + EXP(y*eta)) = probability of the opposite label
        e = EXP(-ABS(y * eta))
        p = MERGE(e / (1.0D0 + e), 1.0D0 / (1.0D0 + e), y * eta >= 0.0D0)
        w = p * (1.0D0 - p)
        ! gradient and Hessian of the mean loss in (b0, b)
        g0 = -SUM(y * p) / DBLE(nobs)
        g1 = -SUM(y * p * xs) / DBLE(nobs)
        h00 = SUM(w) / DBLE(nobs)
        h01 = SUM(w * xs) / DBLE(nobs)
        h11 = SUM(w * xs * xs) / DBLE(nobs)
        det = h00 * h11 - h01 * h01
        IF (det <= EPSILON(1.0D0) * h00 * h11) EXIT
        d0 = -(h11 * g0 - h01 * g1) / det
        d1 = -(h00 * g1 - h01 * g0) / det
        ! -------- STEP HALVING -------- !
        ! the slack keeps rounding noise in loss from rejecting full Newton
        ! steps near the minimum, where the true decrease is below precision
        step = 1.0D0
        DO k = 1, 30
            etanew = (b0 + step * d0) + (b + step * d1) * xs
            lossnew = SUM(MAX(-y * etanew, 0.0D0) + &
                & LOG(1.0D0 + EXP(-ABS(y * etanew)))) / DBLE(nobs)
            IF (lossnew <= loss * (1.0D0 + 1.0D-10)) EXIT
            step = 0.5D0 * step
        END DO
        IF (.NOT. (lossnew <= loss * (1.0D0 + 1.0D-10))) EXIT
        b0 = b0 + step * d0
        b = b + step * d1
        eta = etanew
        loss = lossnew
        IF (MAX(ABS(d0), ABS(d1)) < eps) THEN
            conv = .TRUE.
            EXIT
        END IF
    END DO
    IF (.NOT. conv) jerr = -1

    ! ! -------- FITTED PROBABILITIES -------- !
    ! ! eta always matches the returned (b0, b); prob = 1/(1 + EXP(-eta))
    ! e = EXP(-ABS(eta))
    ! prob = MERGE(1.0D0 / (1.0D0 + e), e / (1.0D0 + e), eta >= 0.0D0)

    ! -------- BACK TO ORIGINAL SCALE -------- !
    beta = b / xsd
    beta0 = b0 - beta * xmean

END SUBROUTINE unilogFIT


SUBROUTINE logitIntercept(nobs, x, y, w, maxit, eps, w0)
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: nobs, maxit
  DOUBLE PRECISION, INTENT(IN) :: x(nobs), y(nobs), w, eps
  DOUBLE PRECISION, INTENT(OUT) :: w0

  INTEGER :: i, iter, npos
  DOUBLE PRECISION :: q, logitq, tmin, tmax
  DOUBLE PRECISION :: lower, upper, eta, p
  DOUBLE PRECISION :: grad, hess, wnew, tol

  w0 = 0.0D0
  IF (nobs <= 0) RETURN

  ! Proportion of positive responses (y = 1).
  npos = COUNT(y > 0.0D0)

  ! A finite MLE requires both response classes.
  ! Use a smoothed proportion if only one class is present.
  IF (npos == 0 .OR. npos == nobs) THEN
    q = (DBLE(npos) + 0.5D0) / (DBLE(nobs) + 1.0D0)
  ELSE
    q = DBLE(npos) / DBLE(nobs)
  END IF

  logitq = LOG(q) - LOG(1.0D0 - q)
  tol = MAX(eps, EPSILON(1.0D0))

  ! Initialize the intercept and its search interval.
  tmin = MIN(w * MINVAL(x), w * MAXVAL(x))
  tmax = MAX(w * MINVAL(x), w * MAXVAL(x))

  lower = logitq - tmax
  upper = logitq - tmin

  w0 = logitq - w * SUM(x) / DBLE(nobs)

  ! Safeguarded Newton-Raphson iteration.
  DO iter = 1, maxit

    grad = 0.0D0
    hess = 0.0D0

    DO i = 1, nobs

      eta = w0 + w * x(i)

      ! Numerically stable sigmoid.
      IF (eta >= 0.0D0) THEN
        p = 1.0D0 / (1.0D0 + EXP(-eta))
      ELSE
        p = EXP(eta) / (1.0D0 + EXP(eta))
      END IF

      grad = grad + (p - q)
      hess = hess + p * (1.0D0 - p)

    END DO

    grad = grad / DBLE(nobs)
    hess = hess / DBLE(nobs)

    ! Check convergence of the score equation.
    IF (ABS(grad) <= tol * q * (1.0D0 - q)) EXIT

    ! Update the bracket.
    IF (grad > 0.0D0) THEN
      upper = w0
    ELSE
      lower = w0
    END IF

    ! Newton step with bisection fallback.
    IF (hess > 0.0D0) THEN
      wnew = w0 - grad / hess
    ELSE
      wnew = 0.5D0 * lower + 0.5D0 * upper
    END IF

    IF (wnew <= lower .OR. wnew >= upper) THEN
      wnew = 0.5D0 * lower + 0.5D0 * upper
    END IF

    ! Check convergence of the intercept.
    IF (ABS(wnew - w0) <= tol * MAX(1.0D0, ABS(w0))) THEN
      w0 = wnew
      EXIT
    END IF

    w0 = wnew

  END DO

END SUBROUTINE logitIntercept

