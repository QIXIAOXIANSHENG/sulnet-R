! Shared logistic grid solver for suniZTlogNET, suniLAlogNET and
! suniSRlogNET. Each cell is fitted by IRLS as in glmnet's lognet: an outer
! loop forms the quadratic approximation of the mean logistic loss at the
! current fit, and coordinate descent solves the penalized weighted least
! squares problem. Along the first coordinate, each second-coordinate level
! stops early once its fit saturates, as glmnet stops its lambda path.

! Traverse the grid in the order of ind (columns 1:2 current indices, 3:4
! seed indices). Row l of pen holds the multipliers of traversal cell l: the
! positive threshold of predictor k is pen(l, 1)*pf(k) and the negative
! threshold is pen(l, 2)*pf3(k) + pen(l, 3)*pf(k).
!
! A second-coordinate level stops after a cell whose fit is saturated
! (fitted probabilities numerically 0 or 1), or, once the level has mnl
! computed cells, after a cell whose deviance ratio exceeds devmax. The
! level's remaining cells get ntheta = -1, as cells skipped for pmax do, and
! are reported as uncomputed. glmnet's third rule, a deviance ratio that
! stops growing, is not used: with a nearly unpenalized side (alpha = 1, a
! large rate) the deviance can stay flat for several cells until the other
! side's threshold drops, and then rise again.
!
! maxit is a pass budget for each second-coordinate level. A cell that does
! not converge within its level's remaining budget is left uncomputed and
! stops its level, as glmnet returns the solutions for the larger lambdas;
! the traversal continues with the other levels. jerr is then -l for the
! first such cell l (replacing an earlier pmax code), otherwise -10000 - l
! for the first cell skipped for pmax.
SUBROUTINE sunilogpath(lam2, nobs, nvars, x, y, ju, pf, pf2, pf3, pmax, &
    & n1, n2, pen, ind, eps, maxit, intr, napenalty, theta0, theta, &
    & itheta, ntheta, npass, jerr)

    IMPLICIT NONE

    DOUBLE PRECISION, PARAMETER :: devmax = 0.999D0
    INTEGER, PARAMETER :: mnl = 5

    INTEGER :: nobs
    INTEGER :: nvars
    INTEGER :: pmax
    INTEGER :: n1
    INTEGER :: n2
    INTEGER :: maxit
    INTEGER :: intr
    INTEGER :: napenalty
    INTEGER :: npass
    INTEGER :: jerr
    INTEGER :: ju(nvars)
    INTEGER :: ind(n1*n2, 4)
    INTEGER :: itheta(pmax, n1, n2)
    INTEGER :: ntheta(n1, n2)

    DOUBLE PRECISION :: lam2
    DOUBLE PRECISION :: eps
    DOUBLE PRECISION :: x(nobs, nvars)
    DOUBLE PRECISION :: y(nobs)
    DOUBLE PRECISION :: pf(nvars)
    DOUBLE PRECISION :: pf2(nvars)
    DOUBLE PRECISION :: pf3(nvars)
    DOUBLE PRECISION :: pen(n1*n2, 3)
    DOUBLE PRECISION :: theta0(1, n1, n2)
    DOUBLE PRECISION :: theta(pmax, n1, n2)

    DOUBLE PRECISION :: thr
    DOUBLE PRECISION :: nulldev
    DOUBLE PRECISION :: ybar
    DOUBLE PRECISION :: dev
    DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: t

    INTEGER :: l
    INTEGER :: j
    INTEGER :: k
    INTEGER :: i1
    INTEGER :: i2
    INTEGER :: p1
    INTEGER :: p2
    INTEGER :: ni
    INTEGER :: status
    INTEGER :: npass0
    INTEGER :: ierr
    INTEGER, DIMENSION(:), ALLOCATABLE :: mm
    INTEGER, DIMENSION(:), ALLOCATABLE :: nlevel
    INTEGER, DIMENSION(:), ALLOCATABLE :: lpass
    INTEGER, DIMENSION(:), ALLOCATABLE :: stopped

    ALLOCATE (t(0:nvars), STAT=jerr)
    ALLOCATE (mm(1:nvars), STAT=ierr)
    jerr = jerr + ierr
    ALLOCATE (nlevel(1:n2), STAT=ierr)
    jerr = jerr + ierr
    ALLOCATE (lpass(1:n2), STAT=ierr)
    jerr = jerr + ierr
    ALLOCATE (stopped(1:n2), STAT=ierr)
    jerr = jerr + ierr
    IF (jerr /= 0) RETURN

    ! convergence threshold: eps times the null deviance per observation
    ! (y is coded -1/+1; 2*log(2) without intercept or with a single class)
    nulldev = 2.0D0*LOG(2.0D0)
    IF (intr == 1) THEN
        ybar = DBLE(COUNT(y > 0.0D0))/DBLE(nobs)
        IF (ybar > 0.0D0 .AND. ybar < 1.0D0) THEN
            nulldev = -2.0D0*(ybar*LOG(ybar) + (1.0D0 - ybar)*LOG(1.0D0 - ybar))
        END IF
    END IF
    thr = eps*nulldev

    nlevel = 0
    lpass = 0
    stopped = 0

    DO l = 1, n1*n2
        i1 = ind(l, 1)
        i2 = ind(l, 2)
        IF (stopped(i2) == 1) THEN
            ntheta(i1, i2) = -1
            napenalty = l
            CYCLE
        END IF
        ! Seed from the recorded predecessor; if it was not computed (pmax
        ! skip or a stopped level), from the previous first-coordinate value
        ! of this level; otherwise start from zero.
        p1 = ind(l, 3)
        p2 = ind(l, 4)
        IF (p1 > 0 .AND. p2 > 0) THEN
            IF (ntheta(p1, p2) < 0) p1 = 0
        END IF
        IF (p1 <= 0 .OR. p2 <= 0) THEN
            p1 = i1 - 1
            p2 = i2
            IF (p1 > 0) THEN
                IF (ntheta(p1, p2) < 0) p1 = 0
            END IF
        END IF
        t = 0.0D0
        mm = 0
        ni = 0
        IF (p1 > 0 .AND. p2 > 0) THEN
            ni = ntheta(p1, p2)
            t(0) = theta0(1, p1, p2)
            DO j = 1, ni
                k = itheta(j, p1, p2)
                itheta(j, i1, i2) = k
                mm(k) = j
                t(k) = theta(j, p1, p2)
            END DO
        END IF

        npass0 = npass
        CALL sunilogcell(lam2, nobs, nvars, x, y, ju, pf, pf2, pf3, pmax, &
            & pen(l, 1), pen(l, 2), pen(l, 3), thr, &
            & npass + maxit - lpass(i2), intr, t, ni, itheta(:, i1, i2), &
            & mm, npass, dev, status)
        lpass(i2) = lpass(i2) + npass - npass0

        IF (status == 3) THEN
            IF (jerr == 0 .OR. jerr < -10000) jerr = -l
            stopped(i2) = 1
            ntheta(i1, i2) = -1
            napenalty = l
            CYCLE
        END IF
        IF (status == 2) THEN
            ! skip this cell and continue the grid; R reports the first skip
            IF (jerr == 0) jerr = -10000 - l
            ntheta(i1, i2) = -1
            napenalty = l
            CYCLE
        END IF
        IF (ni > 0) theta(1:ni, i1, i2) = t(itheta(1:ni, i1, i2))
        ntheta(i1, i2) = ni
        theta0(1, i1, i2) = t(0)
        napenalty = l

        ! -------- PATH TERMINATION -------- !
        nlevel(i2) = nlevel(i2) + 1
        IF (status == 1) stopped(i2) = 1
        IF (nlevel(i2) >= mnl .AND. 1.0D0 - dev/nulldev > devmax) THEN
            stopped(i2) = 1
        END IF
    END DO
    DEALLOCATE (t, mm, nlevel, lpass, stopped)
    RETURN
END SUBROUTINE sunilogpath


! Fit one cell by IRLS from the warm start in t, whose active set is
! m(1:ni) with mm(k) the slot of predictor k. Thresholds are
! ppos*pf(k) on the positive side and pneg3*pf3(k) + pneg1*pf(k) on the
! negative side; the quadratic coordinate denominator adds lam2*pf2(k).
! Fitted probabilities beyond pmin of 0 or 1 are set to 0 or 1 and get zero
! weight. Convergence is measured, as in glmnet, by the largest change of a
! coefficient squared and weighted by its current curvature. On return dev is
! the deviance per observation and status is 0 (converged), 1 (saturated:
! the total weight fell to vmin), 2 (more than pmax active predictors) or
! 3 (npass exceeded maxpass).
SUBROUTINE sunilogcell(lam2, nobs, nvars, x, y, ju, pf, pf2, pf3, pmax, &
    & ppos, pneg3, pneg1, thr, maxpass, intr, t, ni, m, mm, npass, dev, status)

    IMPLICIT NONE

    DOUBLE PRECISION, PARAMETER :: pmin = 1.0D-9

    INTEGER :: nobs
    INTEGER :: nvars
    INTEGER :: pmax
    INTEGER :: maxpass
    INTEGER :: intr
    INTEGER :: ni
    INTEGER :: npass
    INTEGER :: status
    INTEGER :: ju(nvars)
    INTEGER :: m(pmax)
    INTEGER :: mm(nvars)

    DOUBLE PRECISION :: lam2
    DOUBLE PRECISION :: ppos
    DOUBLE PRECISION :: pneg3
    DOUBLE PRECISION :: pneg1
    DOUBLE PRECISION :: thr
    DOUBLE PRECISION :: dev
    DOUBLE PRECISION :: x(nobs, nvars)
    DOUBLE PRECISION :: y(nobs)
    DOUBLE PRECISION :: pf(nvars)
    DOUBLE PRECISION :: pf2(nvars)
    DOUBLE PRECISION :: pf3(nvars)
    DOUBLE PRECISION :: t(0:nvars)

    DOUBLE PRECISION :: fmax
    DOUBLE PRECISION :: vmin
    DOUBLE PRECISION :: q
    DOUBLE PRECISION :: xmz
    DOUBLE PRECISION :: d
    DOUBLE PRECISION :: dif
    DOUBLE PRECISION :: oldt
    DOUBLE PRECISION :: u
    DOUBLE PRECISION :: pospenalty
    DOUBLE PRECISION :: negpenalty
    DOUBLE PRECISION :: mrg
    DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: eta
    DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: v
    DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: r
    DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: xv
    DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: oldtheta

    INTEGER :: i
    INTEGER :: j
    INTEGER :: k
    INTEGER :: vrg

    ALLOCATE (eta(1:nobs), v(1:nobs), r(1:nobs), xv(1:nvars), &
        & oldtheta(0:nvars))
    fmax = LOG(1.0D0/pmin - 1.0D0)
    vmin = (1.0D0 + pmin)*pmin*(1.0D0 - pmin)
    status = 0
    xv = 0.0D0
    oldtheta = 0.0D0

    eta = t(0)
    DO j = 1, ni
        k = m(j)
        IF (t(k) /= 0.0D0) eta = eta + x(:, k)*t(k)
    END DO

    ! -------- IRLS LOOP -------- !
    DO
        ! quadratic approximation at eta; v holds the weights p(1-p)/n and
        ! r the weighted working residuals (y01 - p)/n
        DO i = 1, nobs
            IF (eta(i) < -fmax) THEN
                q = 0.0D0
            ELSE IF (eta(i) > fmax) THEN
                q = 1.0D0
            ELSE
                q = 1.0D0/(1.0D0 + EXP(-eta(i)))
            END IF
            v(i) = q*(1.0D0 - q)/DBLE(nobs)
            r(i) = (0.5D0*(1.0D0 + y(i)) - q)/DBLE(nobs)
        END DO
        xmz = SUM(v)
        IF (xmz <= vmin) THEN
            status = 1
            EXIT
        END IF
        DO k = 1, nvars
            IF (ju(k) /= 0) xv(k) = DOT_PRODUCT(v, x(:, k)**2)
        END DO
        IF (intr == 1) oldtheta(0) = t(0)
        IF (ni > 0) oldtheta(m(1:ni)) = t(m(1:ni))
        ! -------- MIDDLE LOOP -------- !
        DO
            npass = npass + 1
            dif = 0.0D0
            DO k = 1, nvars
                IF (ju(k) /= 0) THEN
                    oldt = t(k)
                    u = DOT_PRODUCT(r, x(:, k)) + xv(k)*t(k)
                    pospenalty = ppos*pf(k)
                    negpenalty = pneg3*pf3(k) + pneg1*pf(k)
                    IF (u > pospenalty) THEN
                        t(k) = (u - pospenalty)/(xv(k) + pf2(k)*lam2)
                    ELSE IF (u < -negpenalty) THEN
                        t(k) = (u + negpenalty)/(xv(k) + pf2(k)*lam2)
                    ELSE
                        t(k) = 0.0D0
                    END IF
                    d = t(k) - oldt
                    IF (ABS(d) > 0.0D0) THEN
                        dif = MAX(dif, xv(k)*d**2)
                        r = r - d*v*x(:, k)
                        IF (mm(k) == 0) THEN
                            ni = ni + 1
                            IF (ni > pmax) EXIT
                            mm(k) = ni
                            m(ni) = k
                        END IF
                    END IF
                END IF
            END DO
            IF (ni > pmax) THEN
                status = 2
                EXIT
            END IF
            IF (intr == 1) THEN
                d = SUM(r)/xmz
                IF (d /= 0.0D0) THEN
                    t(0) = t(0) + d
                    r = r - d*v
                    dif = MAX(dif, xmz*d**2)
                END IF
            END IF
            IF (dif < thr) EXIT
            IF (npass > maxpass) THEN
                status = 3
                EXIT
            END IF
            ! -------- INNER LOOP -------- !
            DO
                npass = npass + 1
                dif = 0.0D0
                DO j = 1, ni
                    k = m(j)
                    oldt = t(k)
                    u = DOT_PRODUCT(r, x(:, k)) + xv(k)*t(k)
                    pospenalty = ppos*pf(k)
                    negpenalty = pneg3*pf3(k) + pneg1*pf(k)
                    IF (u > pospenalty) THEN
                        t(k) = (u - pospenalty)/(xv(k) + pf2(k)*lam2)
                    ELSE IF (u < -negpenalty) THEN
                        t(k) = (u + negpenalty)/(xv(k) + pf2(k)*lam2)
                    ELSE
                        t(k) = 0.0D0
                    END IF
                    d = t(k) - oldt
                    IF (ABS(d) > 0.0D0) THEN
                        dif = MAX(dif, xv(k)*d**2)
                        r = r - d*v*x(:, k)
                    END IF
                END DO
                IF (intr == 1) THEN
                    d = SUM(r)/xmz
                    IF (d /= 0.0D0) THEN
                        t(0) = t(0) + d
                        r = r - d*v
                        dif = MAX(dif, xmz*d**2)
                    END IF
                END IF
                IF (dif < thr) EXIT
                IF (npass > maxpass) THEN
                    status = 3
                    EXIT
                END IF
            END DO
            IF (status == 3) EXIT
        END DO
        IF (status /= 0) EXIT
        eta = t(0)
        DO j = 1, ni
            k = m(j)
            IF (t(k) /= 0.0D0) eta = eta + x(:, k)*t(k)
        END DO
        !-------- FINAL CHECK -------- !
        vrg = 1
        IF (intr == 1) THEN
            IF (xmz*(t(0) - oldtheta(0))**2 >= thr) vrg = 0
        END IF
        DO j = 1, ni
            k = m(j)
            IF (xv(k)*(t(k) - oldtheta(k))**2 >= thr) THEN
                vrg = 0
                EXIT
            END IF
        END DO
        IF (vrg == 1) EXIT
    END DO

    ! deviance per observation, log(1 + exp(-mrg)) evaluated stably
    dev = 0.0D0
    IF (status <= 1) THEN
        DO i = 1, nobs
            mrg = y(i)*eta(i)
            dev = dev + MAX(-mrg, 0.0D0) + LOG(1.0D0 + EXP(-ABS(mrg)))
        END DO
        dev = 2.0D0*dev/DBLE(nobs)
    END IF
    DEALLOCATE (eta, v, r, xv, oldtheta)
    RETURN
END SUBROUTINE sunilogcell
