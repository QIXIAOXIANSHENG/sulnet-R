! Construct the univariate transformed design used by the grid solvers.
! Nonzero beta_ju entries mark supplied priors; other entries use univariate
! slopes. Scale priors by medr/medpj, then return beta, predictor offsets
! beta0, and fit(:,j) = x(:,j)*beta(j) + beta0(j). No LOO step is performed.
SUBROUTINE betaFIT(nobs, nvars, x, y, beta0, beta, beta_ju, medpj, fit)
    ! -----------------------------------------------------------------------
    IMPLICIT NONE
    ! -------- INPUT VARIABLES -------- !
    INTEGER :: nobs
    INTEGER :: nvars
    INTEGER, INTENT(IN) :: beta_ju(nvars)
    DOUBLE PRECISION :: x(nobs, nvars)
    DOUBLE PRECISION :: y(nobs)
    DOUBLE PRECISION, INTENT(IN) :: medpj
    ! -------- OUTPUT VARIABLES -------- !
    DOUBLE PRECISION :: beta0(nvars)
    DOUBLE PRECISION :: beta(nvars)
    DOUBLE PRECISION :: fit(nobs, nvars)
    ! -------- LOCAL DECLARATIONS -------- !
    INTEGER :: j
    ! INTEGER :: lenju
    INTEGER :: ju(nvars)
    ! INTEGER :: ierr
    ! DOUBLE PRECISION :: beta0_temp
    ! DOUBLE PRECISION :: beta_temp
    DOUBLE PRECISION :: r
    DOUBLE PRECISION :: sdratio(nvars)
    DOUBLE PRECISION :: sjdsy
    DOUBLE PRECISION :: medr
    DOUBLE PRECISION :: c
    ! DOUBLE PRECISION :: Ri(nobs, nvars)
    ! DOUBLE PRECISION :: x_new(nobs, nvars)
    ! DOUBLE PRECISION :: maj(nvars)
    ! DOUBLE PRECISION :: xnorm(nvars)
    DOUBLE PRECISION :: xmj
    DOUBLE PRECISION :: ymean
    DOUBLE PRECISION :: xmean(nvars)
    ! DOUBLE PRECISION :: Z(nobs, nvars)
    ! DOUBLE PRECISION :: XY(nobs, nvars)
    ! DOUBLE PRECISION :: beta0_mat(nobs, nvars)
    ! DOUBLE PRECISION :: beta_mat(nobs, nvars)
    DOUBLE PRECISION :: xycor(nvars)
    ! DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: pj

    ! lenju = SUM(beta_ju)
    ! IF(lenju > 0) THEN
    !     ALLOCATE(pj(lenju), STAT = ierr)
    !     IF(ierr /=0 ) RETURN
    ! ELSE
    !     ALLOCATE(pj(1:nvars), STAT = ierr)
    !     IF(ierr /=0 ) RETURN
    !     lenju = nvars
    ! END IF
    CALL chkvars(nobs, nvars, x, ju)

    ymean = SUM(y)/DBLE(nobs)

    DO j = 1, nvars
        xmj = 0.0D0
        r = 0.0D0
        sjdsy = 0.0D0
        IF (ju(j) /= 0) CALL priorStandard(nobs, x(:,j), y, xmj, ymean, sjdsy, r)
        sdratio(j) = sjdsy
        xycor(j) = r
        xmean(j) = xmj
    END DO

    CALL median(ABS(xycor), nvars, ju, medr)
    c = medr / medpj

    DO j = 1, nvars
        IF(beta_ju(j) /= 0) THEN
            beta(j) = beta(j) * c
        ELSE
            beta(j) = sdratio(j) * xycor(j)
        END IF

        fit(:, j) = x(:, j) * beta(j)
        beta0(j) = ymean - SUM(fit(:, j)) / DBLE(nobs)
        fit(:, j) = fit(:, j) + beta0(j)
    END DO

END SUBROUTINE betaFIT

! Logistic version of betaFIT with y in {-1, 1}. jerr(j) is the unilogFIT code
! for predictor j: -1 (no convergence, typically separation) keeps the last
! iterate; 8001 or 8002 (invalid or single-class y) returns immediately.
SUBROUTINE betalogFIT(nobs, nvars, x, y, beta0, beta, beta_ju, medpj, fit, jerr)
    ! -----------------------------------------------------------------------
    IMPLICIT NONE
    ! -------- INPUT VARIABLES -------- !
    INTEGER :: nobs
    INTEGER :: nvars
    INTEGER, INTENT(IN) :: beta_ju(nvars)
    DOUBLE PRECISION :: x(nobs, nvars)
    DOUBLE PRECISION :: y(nobs)
    DOUBLE PRECISION, INTENT(IN) :: medpj
    ! -------- OUTPUT VARIABLES -------- !
    DOUBLE PRECISION :: beta0(nvars)
    DOUBLE PRECISION :: beta(nvars)
    DOUBLE PRECISION :: fit(nobs, nvars)
    INTEGER :: jerr(nvars)
    ! -------- LOCAL DECLARATIONS -------- !
    INTEGER :: j
    INTEGER :: ju(nvars)
    INTEGER, PARAMETER :: maxit = 1D5
    INTEGER :: npass
    DOUBLE PRECISION, PARAMETER :: eps = 1.0D-8
    DOUBLE PRECISION :: r(nvars)
    DOUBLE PRECISION :: medr
    DOUBLE PRECISION :: c
    DOUBLE PRECISION :: b
    DOUBLE PRECISION :: b0
    DOUBLE PRECISION :: xnorm(nvars)
    DOUBLE PRECISION :: xsd
    DOUBLE PRECISION :: w0

    r = 0.0D0
    jerr = 0

    CALL chkvars(nobs, nvars, x, ju)
    DO j = 1, nvars
        b = 0.0D0
        b0 = 0.0D0
        npass = 0
        xsd = 0.0D0
        IF (ju(j) /= 0) CALL unilogFIT(nobs, x(:, j), y, xsd, b0, b, maxit, eps, npass, jerr(j))
        IF (jerr(j) == 8001 .OR. jerr(j) == 8002) RETURN
        beta0(j) = b0
        xnorm(j) = xsd
        r(j) = b*xsd
    END DO
     
    CALL median(ABS(r), nvars, ju, medr)
    c = medr / medpj

    DO j = 1, nvars
        w0 = 0.0D0
        IF(ju(j) /= 0) THEN
            IF(beta_ju(j) /= 0) THEN
                beta(j) = beta(j) * c / xnorm(j)
                CALL logitIntercept(nobs, x(:, j), y, beta(j), maxit, eps, w0)
                beta0(j) = w0
            ELSE
                beta(j) = r(j) / xnorm(j)
            END IF
            fit(:, j) = x(:, j) * beta(j) + beta0(j)
        END IF
    END DO
END SUBROUTINE betalogFIT



