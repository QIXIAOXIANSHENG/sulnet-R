! DESCRIPTION:
!
!    Functions standard and chkvars are minor modifications from the `glmnet` package:
!
!    Jerome Friedman, Trevor Hastie, Robert Tibshirani (2010).
!    Regularization Paths for Generalized Linear Models via Coordinate Descent.
!    Journal of Statistical Software, 33(1), 1-22.
!    URL: https://www.jstatsoft.org/v33/i01/.
!
!    Function loofit is a minor modification from the `uniLasso` package:
!
!    Chatterjee, S., Hastie, T., & Tibshirani, R. (2025).
!    Univariate‑Guided Sparse Regression (arXiv:2501.18360v9).
!    arXiv.
!    URL: https://doi.org/10.48550/arXiv.2501.18360
!
! --------------------------------------------------------------------------
! standard: An auxiliary function for standardize x matrix.
! --------------------------------------------------------------------------
!
! USAGE:
!
! CALL standard(nobs, nvars, x, ju, isd, intr, xmean, xnorm, maj)
!
! INPUT ARGUMENTS:
!
!    nobs = number of observations
!    nvars = number of predictor variables
!    x(nobs, nvars) = matrix of predictors, of dimension N * p; each row is an
!                     observation vector.
!    ju(nvars) = flag of predictor variables
!                ju(j) = 0 => this predictor has zero variance
!                ju(j) = 1 => this predictor does not have zero variance
!    isd = standarization flag:
!          isd = 0 => do not standardize predictor variables
!          isd = 1 => standardize predictor variables
!          NOTE: no matter isd is 1 or 0, matrix x is always centered by column.
!                That is, col.mean(x) = 0.
!
! OUTPUT:
!
!    x(nobs, nvars) = standarized matrix x
!    xmean(nvars) = column mean of x matrix
!    xnorm(nvars) = column standard deviation of x matrix
!    maj(nvars) = column variance of x matrix
!
! --------------------------------------------------------------------------
! chkvars: An auxiliary function for variable check.
! --------------------------------------------------------------------------
!
! USAGE:
!
! CALL chkvars(nobs, nvars, x, ju)
!
! INPUT ARGUMENTS:
!
!    nobs = number of observations
!    nvars = number of predictor variables
!    x(nobs, nvars) = matrix of predictors, of dimension N * p; each row is an
!                     observation vector.
!    y(nobs) = response variable. This argument should be a two-level factor
!              {-1, 1} for classification.
!
! OUTPUT:
!
!    ju(nvars) = flag of predictor variables
!                ju(j) = 0 => this predictor has zero variance
!                ju(j) = 1 => this predictor does not have zero variance
!
! --------------------------------------------------------------------------
! gridIndex: An auxiliary function for ordering the (zeta, tau) grid and
!            finding the warm start of each grid point.
! --------------------------------------------------------------------------
!
! USAGE:
!
! CALL gridIndex(nzeta, ntau, uzeta, utau, ind)
!
! INPUT ARGUMENTS:
!
!    nzeta = number of zeta values
!    ntau = number of tau values
!    uzeta(nzeta) = zeta values, positive and ordered from large to small
!    utau(ntau) = tau values, positive and ordered from large to small
!    ind(nzeta * ntau, 4) = zero matrix
!
! OUTPUT:
!
!    ind(nzeta * ntau, 4) = grid index matrix; rows follow the order of
!                             expand.grid(uzeta, utau), i.e.
!                             row k = (j - 1) * nzeta + i
!       ind(k, 1:2) = (i, j), index of the current zeta and tau
!       ind(k, 3:4) = index of the zeta and tau, among the pairs in rows
!                       1, ..., k - 1 with zeta >= uzeta(i) and tau >= utau(j),
!                       that has the smallest squared Euclidean distance to
!                       (log(uzeta(i)), log(utau(j))) on the log scale;
!                       (0, 0) for the first row
!
! --------------------------------------------------------------------------
! priorStandard: An auxiliary function for the summary statistics of a
!                pair of vectors x and y.
! --------------------------------------------------------------------------
!
! USAGE:
!
! CALL priorStandard(n, x, y, xmean, ymean, sdratio, xycor)
!
! INPUT ARGUMENTS:
!
!    n = length of x and y
!    x(n) = first vector
!    y(n) = second vector
!
! OUTPUT:
!
!    xmean = mean of x
!    ymean = mean of y
!    sdratio = standard deviation ratio S_x / S_y
!    xycor = covariance of x and y
!    NOTE: variances and covariance use the divisor n (not n - 1), the same
!          as in standard. sdratio does not depend on this choice.
!
! --------------------------------------------------------------------------
! median: An auxiliary function for the median of a vector x.
! --------------------------------------------------------------------------
!
! USAGE:
!
! CALL median(x, n, ju, medx)
!
! INPUT ARGUMENTS:
!
!    x(n) = input vector; it is not modified
!    n = length of x
!    ju(n) = inclusion flags (0/1); x(i) is skipped if ju(i) = 0
!
! OUTPUT:
!
!    medx = median of x(ju /= 0); with m = COUNT(ju /= 0), the average of
!           the two middle values when m is even, the same as median() in R.
!           medx = 0 if m < 1.
!    NOTE: uses Wirth's selection algorithm on a copy of x(ju /= 0), O(m) on
!          average.
!


! --------------------------------------------------------------------------
SUBROUTINE standard(nobs, nvars, x, ju, isd, intr, xmean, xnorm, maj)
    ! ------------------------------------------------------------------------
    IMPLICIT NONE
    ! -------- INPUT VARIABLES -------- !
    INTEGER :: nobs
    INTEGER :: nvars
    INTEGER :: isd
    INTEGER :: intr
    INTEGER :: ju(nvars)
    DOUBLE PRECISION :: xmsq
    DOUBLE PRECISION :: xvar
    DOUBLE PRECISION :: x(nobs, nvars)
    DOUBLE PRECISION :: xmean(nvars)
    DOUBLE PRECISION :: xnorm(nvars)
    DOUBLE PRECISION :: maj(nvars)
    ! -------- LOCAL DECLARATIONS -------- !
    INTEGER :: j
    ! -------- STANDARDIZATION -------- !
    IF (intr == 0) THEN
        DO j = 1, nvars
            IF (ju(j) == 1) THEN
                xmean(j) = 0.0D0
                maj(j) = DOT_PRODUCT(x(:, j), x(:, j))/nobs
                IF (isd == 1) THEN
                    xmsq = (SUM(x(:, j))/nobs)**2
                    xvar = maj(j) - xmsq
                    xnorm(j) = SQRT(xvar)
                    x(:, j) = x(:, j)/xnorm(j)
                    maj(j) = 1.0D0 + xmsq/xvar
                END IF
            END IF
        END DO
    ELSE
        DO j = 1, nvars
            IF (ju(j) == 1) THEN
                xmean(j) = SUM(x(:, j))/nobs ! MEAN
                x(:, j) = x(:, j) - xmean(j)
                maj(j) = DOT_PRODUCT(x(:, j), x(:, j))/nobs
                IF (isd == 1) THEN
                    xnorm(j) = SQRT(maj(j)) ! STANDARD DEVIATION
                    x(:, j) = x(:, j)/xnorm(j)
                    maj(j) = 1.0D0
                END IF
            END IF
        END DO
    END IF
END SUBROUTINE standard

! -------------------------------------------------------------------------
SUBROUTINE chkvars(nobs, nvars, x, ju)
    ! -----------------------------------------------------------------------
    IMPLICIT NONE
    ! -------- INPUT VARIABLES -------- !
    INTEGER :: nobs
    INTEGER :: nvars
    INTEGER :: ju(nvars)
    DOUBLE PRECISION :: x(nobs, nvars)
    ! -------- LOCAL DECLARATIONS -------- !
    INTEGER :: i
    INTEGER :: j
    DOUBLE PRECISION :: t
    ! -------- VARIABLE CHECKING -------- !
    DO j = 1, nvars
        ju(j) = 0
        t = x(1, j)
        DO i = 2, nobs
            IF (x(i, j) /= t) THEN
                ju(j) = 1
                EXIT
            END IF
        END DO
    END DO
END SUBROUTINE chkvars


! -------------------------------------------------------------------------
SUBROUTINE gridIndexZT(nzeta, ntau, uzeta, utau, ind, diffzeta, difftau)

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: nzeta, ntau
    DOUBLE PRECISION, INTENT(IN) :: uzeta(nzeta), utau(ntau)
    DOUBLE PRECISION, INTENT(IN) :: diffzeta, difftau
    INTEGER, INTENT(INOUT) :: ind(nzeta*ntau, 4)

    INTEGER :: izeta, itau, k
    DOUBLE PRECISION :: dzeta, dtau
    DOUBLE PRECISION :: big

    big = HUGE(1.0D0)

    ind = 0
    k = 0

    DO itau = 1, ntau

        DO izeta = 1, nzeta

            k = k + 1

            ! Current grid point
            ind(k, 1) = izeta
            ind(k, 2) = itau

            ! The largest zeta and largest tau point has no predecessor
            IF (izeta == 1 .AND. itau == 1) THEN
                ind(k, 3) = 0
                ind(k, 4) = 0
                CYCLE
            END IF

            ! Only the tau-direction predecessor exists
            IF (izeta == 1) THEN
                ind(k, 3) = izeta
                ind(k, 4) = itau - 1
                CYCLE
            END IF

            ! Only the zeta-direction predecessor exists
            IF (itau == 1) THEN
                ind(k, 3) = izeta - 1
                ind(k, 4) = itau
                CYCLE
            END IF

            ! If both grid spacings are known, no distance calculation is needed
            IF (diffzeta >= 0.0D0 .AND. difftau >= 0.0D0) THEN

                dzeta = diffzeta
                dtau  = difftau

            ELSE

                ! Calculate the actual log-distance in the zeta direction
                IF (uzeta(izeta) == uzeta(izeta-1)) THEN

                    dzeta = 0.0D0

                ELSE IF (uzeta(izeta) > 0.0D0 .AND. &
                         uzeta(izeta-1) > 0.0D0) THEN

                    dzeta = ABS(LOG(uzeta(izeta)) - &
                                LOG(uzeta(izeta-1)))

                ELSE

                    dzeta = big

                END IF

                ! Calculate the actual log-distance in the tau direction
                IF (utau(itau) == utau(itau-1)) THEN

                    dtau = 0.0D0

                ELSE IF (utau(itau) > 0.0D0 .AND. &
                         utau(itau-1) > 0.0D0) THEN

                    dtau = ABS(LOG(utau(itau)) - &
                               LOG(utau(itau-1)))

                ELSE

                    dtau = big

                END IF

            END IF

            ! Choose the closer predecessor
            IF (dzeta <= dtau) THEN

                ind(k, 3) = izeta - 1
                ind(k, 4) = itau

            ELSE

                ind(k, 3) = izeta
                ind(k, 4) = itau - 1

            END IF

        END DO

    END DO

END SUBROUTINE gridIndexZT


SUBROUTINE gridIndexSR(ns, nr, us, ur, ind)

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: ns, nr
    DOUBLE PRECISION, INTENT(IN) :: us(ns), ur(nr)
    INTEGER, INTENT(OUT) :: ind(ns*nr, 4)

    INTEGER :: i, j, k, l
    INTEGER :: row, bestk, bestl
    DOUBLE PRECISION :: lambda1c, lambda2c
    DOUBLE PRECISION :: lambda1p, lambda2p
    DOUBLE PRECISION :: dist, bestdist
    DOUBLE PRECISION :: tol

    ! us must be ordered from large to small.
    ! ur can be ordered arbitrarily, although increasing order is recommended.
    !
    ! For each current point (us(i), ur(j)), search all points from
    ! stronger S levels and select the closest componentwise-stronger
    ! point in the (log(lambda1), log(lambda2)) space.
    !
    ! lambda1 = S / (1 + r)
    ! lambda2 = S * r / (1 + r)

    ind = 0
    tol = 100.0D0 * EPSILON(1.0D0)

    DO i = 1, ns

        DO j = 1, nr

            row = (i - 1) * nr + j

            ind(row, 1) = i
            ind(row, 2) = j

            ! No stronger S level is available for the first layer.
            IF (i == 1) THEN
                ind(row, 3) = 0
                ind(row, 4) = 0
                CYCLE
            END IF

            lambda1c = us(i) / (1.0D0 + ur(j))
            lambda2c = us(i) * ur(j) / (1.0D0 + ur(j))

            bestdist = HUGE(1.0D0)
            bestk = 0
            bestl = 0

            ! Search from the nearest stronger S layer first.
            DO k = i - 1, 1, -1

                DO l = 1, nr

                    lambda1p = us(k) / (1.0D0 + ur(l))
                    lambda2p = us(k) * ur(l) / (1.0D0 + ur(l))

                    ! The warm-start point must be at least as strongly
                    ! penalized in both directions.
                    IF (lambda1p >= lambda1c * (1.0D0 - tol) .AND. &
                        lambda2p >= lambda2c * (1.0D0 - tol)) THEN

                        dist = (LOG(lambda1p) - LOG(lambda1c))**2 + &
                               (LOG(lambda2p) - LOG(lambda2c))**2

                        IF (dist < bestdist) THEN
                            bestdist = dist
                            bestk = k
                            bestl = l
                        END IF

                    END IF

                END DO

            END DO

            ind(row, 3) = bestk
            ind(row, 4) = bestl

        END DO

    END DO

END SUBROUTINE gridIndexSR

SUBROUTINE gridIndexLA(nlambda, nalpha, ulambda, ualpha, ind)

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: nlambda, nalpha
    DOUBLE PRECISION, INTENT(IN) :: ulambda(nlambda)
    DOUBLE PRECISION, INTENT(IN) :: ualpha(nalpha)
    INTEGER, INTENT(OUT) :: ind(nlambda*nalpha, 4)

    INTEGER :: i, j, k, l
    INTEGER :: row, bestk, bestl
    DOUBLE PRECISION :: lambda1c, lambda2c
    DOUBLE PRECISION :: lambda1p, lambda2p
    DOUBLE PRECISION :: dist, bestdist
    DOUBLE PRECISION :: tol

    ! ulambda must be ordered from large to small.
    ! ualpha is assumed to satisfy 0.5 <= alpha <= 1.0.
    !
    ! For each current point (lambda(i), alpha(j)), search all points
    ! from stronger total-penalty levels and select the closest
    ! componentwise-stronger point.
    !
    ! lambda1 = lambda * (1 - alpha)
    ! lambda2 = lambda * alpha
    !
    ! When alpha = 1, lambda1 = 0 and a logarithmic distance cannot
    ! be used for lambda1. In that case, the lambda1 difference is
    ! normalized by the current total penalty.

    ind = 0
    tol = 100.0D0 * EPSILON(1.0D0)

    DO i = 1, nlambda

        DO j = 1, nalpha

            row = (i - 1) * nalpha + j

            ind(row, 1) = i
            ind(row, 2) = j

            ! No stronger total-penalty level is available for the first layer.
            IF (i == 1) THEN
                ind(row, 3) = 0
                ind(row, 4) = 0
                CYCLE
            END IF

            lambda1c = ulambda(i) * (1.0D0 - ualpha(j))
            lambda2c = ulambda(i) * ualpha(j)

            bestdist = HUGE(1.0D0)
            bestk = 0
            bestl = 0

            ! Search from the nearest stronger lambda layer first.
            DO k = i - 1, 1, -1

                DO l = 1, nalpha

                    lambda1p = ulambda(k) * (1.0D0 - ualpha(l))
                    lambda2p = ulambda(k) * ualpha(l)

                    ! The warm-start point must be at least as strongly
                    ! penalized in both directions.
                    IF (lambda1p >= lambda1c - tol * ulambda(i) .AND. &
                        lambda2p >= lambda2c * (1.0D0 - tol)) THEN

                        IF (lambda1c > 0.0D0) THEN

                            dist = (LOG(lambda1p) - LOG(lambda1c))**2 + &
                                   (LOG(lambda2p) - LOG(lambda2c))**2

                        ELSE

                            ! Boundary case alpha = 1.
                            dist = (lambda1p / ulambda(i))**2 + &
                                   (LOG(lambda2p) - LOG(lambda2c))**2

                        END IF

                        IF (dist < bestdist) THEN
                            bestdist = dist
                            bestk = k
                            bestl = l
                        END IF

                    END IF

                END DO

            END DO

            ind(row, 3) = bestk
            ind(row, 4) = bestl

        END DO

    END DO

END SUBROUTINE gridIndexLA

! -------------------------------------------------------------------------
SUBROUTINE priorStandard(n, x, y, xmean, ymean, sdratio, xycor)
    ! -----------------------------------------------------------------------
    IMPLICIT NONE
    ! -------- INPUT VARIABLES -------- !
    INTEGER, INTENT(IN) :: n
    DOUBLE PRECISION, INTENT(IN) :: x(n)
    DOUBLE PRECISION, INTENT(IN) :: y(n)
    DOUBLE PRECISION, INTENT(IN) :: ymean
    ! -------- OUTPUT VARIABLES -------- !
    DOUBLE PRECISION, INTENT(OUT) :: xmean
    DOUBLE PRECISION, INTENT(OUT) :: sdratio
    DOUBLE PRECISION, INTENT(OUT) :: xycor
    ! -------- LOCAL DECLARATIONS -------- !
    INTEGER :: i
    DOUBLE PRECISION :: dx
    DOUBLE PRECISION :: dy
    DOUBLE PRECISION :: sxx
    DOUBLE PRECISION :: syy
    DOUBLE PRECISION :: sxy
    ! -------- MEANS -------- !
    xmean = SUM(x)/n
    ! ymean = SUM(y)/n
    ! -------- CENTERED CROSS PRODUCTS -------- !
    sxx = 0.0D0
    syy = 0.0D0
    sxy = 0.0D0
    DO i = 1, n
        dx = x(i) - xmean
        dy = y(i) - ymean
        sxx = sxx + dx * dx
        syy = syy + dy * dy
        sxy = sxy + dx * dy
    END DO
    ! -------- SD RATIO AND COVARIANCE -------- !
    sdratio = SQRT(sxx/syy)
    xycor = sxy/SQRT(sxx*syy)
END SUBROUTINE priorStandard


! -------------------------------------------------------------------------
SUBROUTINE median(x, n, ju, medx)
    ! -----------------------------------------------------------------------
    IMPLICIT NONE
    ! -------- INPUT VARIABLES -------- !
    INTEGER, INTENT(IN) :: n
    INTEGER, INTENT(IN) :: ju(n)
    DOUBLE PRECISION, INTENT(IN) :: x(n)
    ! -------- OUTPUT VARIABLES -------- !
    DOUBLE PRECISION, INTENT(OUT) :: medx
    ! -------- LOCAL DECLARATIONS -------- !
    INTEGER :: i
    INTEGER :: j
    INTEGER :: k
    INTEGER :: l
    INTEGER :: r
    INTEGER :: m
    DOUBLE PRECISION :: pivot
    DOUBLE PRECISION :: tmp
    DOUBLE PRECISION, ALLOCATABLE :: a(:)
    ! -------- NUMBER OF INCLUDED VALUES -------- !
    m = COUNT(ju /= 0)
    ! -------- EMPTY INPUT -------- !
    IF (m < 1) THEN
        medx = 0.0D0
        RETURN
    END IF
    ! -------- WORKING COPY OF x(ju /= 0) -------- !
    ALLOCATE (a(1:m))
    a = PACK(x, ju /= 0)
    ! -------- SELECT THE k-TH SMALLEST VALUE (WIRTH) -------- !
    ! on exit: a(1:k-1) <= a(k) <= a(k+1:m)
    k = (m + 1)/2
    l = 1
    r = m
    DO WHILE (l < r)
        pivot = a(k)
        i = l
        j = r
        DO
            DO WHILE (a(i) < pivot)
                i = i + 1
            END DO
            DO WHILE (pivot < a(j))
                j = j - 1
            END DO
            IF (i <= j) THEN
                tmp = a(i)
                a(i) = a(j)
                a(j) = tmp
                i = i + 1
                j = j - 1
            END IF
            IF (i > j) EXIT
        END DO
        IF (j < k) l = i
        IF (k < i) r = j
    END DO
    ! -------- MEDIAN -------- !
    IF (MOD(m, 2) == 1) THEN
        medx = a(k)
    ELSE
        ! the (k+1)-th smallest value is the minimum of a(k+1:m)
        medx = 0.5D0 * (a(k) + MINVAL(a(k+1:m)))
    END IF
    DEALLOCATE (a)
END SUBROUTINE median
