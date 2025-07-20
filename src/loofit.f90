! DESCRIPTION:
!
!    Function loofit is a minor modification from the `uniLasso` package:
!
!    Chatterjee, S., Hastie, T., & Tibshirani, R. (2025).
!    Univariate‑Guided Sparse Regression (arXiv:2501.18360v9).
!    arXiv.
!    URL: https://doi.org/10.48550/arXiv.2501.18360
!
! --------------------------------------------------------------------------
! loofit: An auxiliary function for leave-one-out cross-validation.
! --------------------------------------------------------------------------
!
! USAGE:
!
! CALL loofit(nobs, nvars, x, y, loo, beta0, beta, fit)
!
! INPUT ARGUMENTS:
!
!    nobs = number of observations
!    nvars = number of predictor variables
!    x(nobs, nvars) = matrix of predictors, of dimension N * p; each row is an
!                     observation vector.
!    y(nobs) = response variable. This argument should be a two-level factor
!              {-1, 1} for classification.
!    loo = logical flag for leave-one-out univariate fit; if true, then compute
!          the leave-one-out univariate fit, otherwise compute the standard fit.
!
!
! OUTPUT:
!
!    beta0(nvars) = intercepts for each univariate predictor variable
!    beta(nvars) = coefficients for each univariate predictor variable
!    fit(nobs, nvars) = fitted values for each observation and predictor; if loo
!                      is true, then this is the leave-one-out fit, otherwise
!                      this is the standard fit.
! --------------------------------------------------------------------------


! --------------------------------------------------------------------------
SUBROUTINE loofit(nobs, nvars, x, y, loo, beta0, beta, fit)
    ! -----------------------------------------------------------------------
    IMPLICIT NONE
    ! -------- INPUT VARIABLES -------- !
    INTEGER :: nobs
    INTEGER :: nvars
    DOUBLE PRECISION :: x(nobs, nvars)
    DOUBLE PRECISION :: y(nobs)
    LOGICAL, INTENT(IN) :: loo
    ! -------- OUTPUT VARIABLES -------- !
    DOUBLE PRECISION :: beta0(nvars)
    DOUBLE PRECISION :: beta(nvars)
    DOUBLE PRECISION :: fit(nobs, nvars)
    ! -------- LOCAL DECLARATIONS -------- !
    INTEGER :: j
    INTEGER :: ju(nvars)
    DOUBLE PRECISION :: beta0_temp
    DOUBLE PRECISION :: beta_temp
    DOUBLE PRECISION :: Ri(nobs, nvars)
    DOUBLE PRECISION :: x_new(nobs, nvars)
    DOUBLE PRECISION :: maj(nvars)
    DOUBLE PRECISION :: xnorm(nvars)
    DOUBLE PRECISION :: xmean(nvars)
    DOUBLE PRECISION :: Z(nobs, nvars)
    DOUBLE PRECISION :: XY(nobs, nvars)
    DOUBLE PRECISION :: beta0_mat(nobs, nvars)
    DOUBLE PRECISION :: beta_mat(nobs, nvars)

    call chkvars(nobs, nvars, x, ju)

    ! -------- DEFINE VARIABLES -------- !
    x_new = x
    DO j = 1, nvars
        Z(:, j) = y
    END DO

    ! -------- CENTERING AND STANDARDIZATION -------- !
    CALL standard(nobs, nvars, x_new, ju, 1, 1, xmean, xnorm, maj)
    XY = x_new * Z

    ! -------- COMPUTE BETA -------- !
    DO j = 1, nvars
        IF (ju(j) == 1) THEN

            beta_temp = SUM(XY(:, j))/nobs
            beta0_temp = SUM(Z(:, j))/nobs

            beta_mat(:, j) = beta_temp
            beta0_mat(:, j) = beta0_temp
            beta0(j) = beta0_temp
            beta(j) = beta_temp
        ELSE
            beta_mat(:, j) =0.0D0
            beta0_mat(:, j) = 0.0D0
            beta(j) = 0.0D0
            beta0(j) = 0.0D0
        END IF
    END DO

    beta = beta / xnorm
    beta0 = beta0 - beta*xmean/xnorm


    ! -------- LEAVE-ONE-OUT FITTING -------- !
    IF (loo) THEN
        Ri = DBLE(nobs) * (Z - beta0_mat - x_new * beta_mat) / &
             (DBLE(nobs) - 1.0D0 - x_new ** 2)

        fit = Z - Ri
    ELSE
        fit = beta0_mat + x_new * beta_mat
    END IF

END SUBROUTINE loofit



