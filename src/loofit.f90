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

    beta0 = beta0 - beta*xmean / xnorm
    beta = beta / xnorm
    


    ! -------- LEAVE-ONE-OUT FITTING -------- !
    IF (loo) THEN
        Ri = DBLE(nobs) * (Z - beta0_mat - x_new * beta_mat) / &
             (DBLE(nobs) - 1.0D0 - x_new ** 2)

        fit = Z - Ri
    ELSE
        fit = beta0_mat + x_new * beta_mat
    END IF

END SUBROUTINE loofit


! --------------------------------------------------------------------------
SUBROUTINE loofit_ts(nobs, nobs0, nvars, x0, y0, xA, yA, loo, beta0, beta, fit)
    ! -----------------------------------------------------------------------
    IMPLICIT NONE
    ! -------- INPUT VARIABLES -------- !
    INTEGER :: nobs
    INTEGER :: nobs0
    INTEGER :: nvars
    DOUBLE PRECISION :: x0(nobs0, nvars)
    DOUBLE PRECISION :: y0(nobs0)
    DOUBLE PRECISION :: xA(nobs, nvars)
    DOUBLE PRECISION :: yA(nobs)
    LOGICAL, INTENT(IN) :: loo
    ! -------- OUTPUT VARIABLES -------- !
    DOUBLE PRECISION :: beta0(nvars)
    DOUBLE PRECISION :: beta(nvars)
    DOUBLE PRECISION :: fit(nobs, nvars)
    ! -------- LOCAL DECLARATIONS -------- !
    INTEGER :: j
    INTEGER :: ju(nvars)
    INTEGER :: juA(nvars)
    DOUBLE PRECISION :: beta0_temp
    DOUBLE PRECISION :: beta_temp
    DOUBLE PRECISION :: Ri(nobs, nvars)
    DOUBLE PRECISION :: x_new0(nobs0, nvars)
    DOUBLE PRECISION :: x_newA(nobs, nvars)
    DOUBLE PRECISION :: maj(nvars)
    DOUBLE PRECISION :: xnorm0(nvars)
    DOUBLE PRECISION :: xmean0(nvars)
    DOUBLE PRECISION :: xnormA(nvars)
    DOUBLE PRECISION :: xmeanA(nvars)
    DOUBLE PRECISION :: ynew(nobs)
    DOUBLE PRECISION :: Z0(nobs0, nvars)
    DOUBLE PRECISION :: ZA(nobs, nvars)
    DOUBLE PRECISION :: XY0(nobs0, nvars)
    DOUBLE PRECISION :: XYA(nobs, nvars)
    DOUBLE PRECISION :: beta0_mat(nobs, nvars)
    DOUBLE PRECISION :: beta_mat(nobs, nvars)

    call chkvars(nobs0, nvars, x0, ju)

    ! -------- DEFINE VARIABLES -------- !
    x_new0 = x0
    x_newA = xA
    ynew = yA
    DO j = 1, nvars
        Z0(:, j) = y0
        ZA(:, j) = yA
    END DO
    call chkvars(nobs, nvars, x_newA, juA)

    ! -------- CENTERING AND STANDARDIZATION -------- !
    CALL standard(nobs, nvars, x_newA, juA, 1, 1, xmeanA, xnormA, maj)
    CALL standard(nobs0, nvars, x_new0, ju, 1, 1, xmean0, xnorm0, maj)
    XY0 = x_new0 * Z0
    XYA = x_newA * ZA

    ! -------- COMPUTE BETA -------- !
    DO j = 1, nvars
        IF (ju(j) == 1) THEN

            beta_temp = SUM(XY0(:, j))/nobs0
            beta0_temp = SUM(Z0(:, j))/nobs0
            beta0(j) = beta0_temp
            beta(j) = beta_temp
        ELSE
            beta(j) = 0.0D0
            beta0(j) = 0.0D0
        END IF
        IF (juA(j) == 1) THEN
            beta_temp = SUM(XYA(:, j))/nobs
            beta0_temp = SUM(ZA(:, j))/nobs
            beta_mat(:, j) = beta_temp
            beta0_mat(:, j) = beta0_temp
        ELSE
            beta_mat(:, j) =0.0D0
            beta0_mat(:, j) = 0.0D0
        END IF
    END DO

    beta = beta / xnorm0
    beta0 = beta0 - beta*xmean0 / xnorm0


    ! -------- LEAVE-ONE-OUT FITTING -------- !
    IF (loo) THEN
        Ri = DBLE(nobs) * (ZA - beta0_mat - x_newA * beta_mat) / &
             (DBLE(nobs) - 1.0D0 - x_newA ** 2)

        fit = ZA - Ri
    ELSE
        fit = beta0_mat + x_newA * beta_mat
    END IF

END SUBROUTINE loofit_ts


! --------------------------------------------------------------------------
SUBROUTINE loofit_s(nobs, nobs0, nvars, x0, y0, xA, yA, loo, beta0, beta, fit)
    ! -----------------------------------------------------------------------
    IMPLICIT NONE
    ! -------- INPUT VARIABLES -------- !
    INTEGER :: nobs
    INTEGER :: nobs0
    INTEGER :: nvars
    DOUBLE PRECISION :: x0(nobs0, nvars)
    DOUBLE PRECISION :: y0(nobs0)
    DOUBLE PRECISION :: xA(nobs, nvars)
    DOUBLE PRECISION :: yA(nobs)
    LOGICAL, INTENT(IN) :: loo
    ! -------- OUTPUT VARIABLES -------- !
    DOUBLE PRECISION :: beta0(nvars)
    DOUBLE PRECISION :: beta(nvars)
    DOUBLE PRECISION :: fit(nobs0, nvars)
    ! -------- LOCAL DECLARATIONS -------- !
    INTEGER :: j
    INTEGER :: ju(nvars)
    INTEGER :: juA(nvars)
    DOUBLE PRECISION :: beta0_temp
    DOUBLE PRECISION :: beta_temp
    DOUBLE PRECISION :: Ri(nobs, nvars)
    DOUBLE PRECISION :: x_new0(nobs0, nvars)
    DOUBLE PRECISION :: x_newA(nobs, nvars)
    DOUBLE PRECISION :: maj(nvars)
    DOUBLE PRECISION :: xnorm0(nvars)
    DOUBLE PRECISION :: xmean0(nvars)
    DOUBLE PRECISION :: xnormA(nvars)
    DOUBLE PRECISION :: xmeanA(nvars)
    DOUBLE PRECISION :: ynew(nobs)
    DOUBLE PRECISION :: Z0(nobs0, nvars)
    DOUBLE PRECISION :: ZA(nobs, nvars)
    DOUBLE PRECISION :: XY0(nobs0, nvars)
    DOUBLE PRECISION :: XYA(nobs, nvars)
    DOUBLE PRECISION :: beta0_mat(nobs, nvars)
    DOUBLE PRECISION :: beta_mat(nobs, nvars)
    DOUBLE PRECISION :: fit_temp(nobs, nvars)

    call chkvars(nobs0, nvars, x0, ju)

    ! -------- DEFINE VARIABLES -------- !
    x_new0 = x0
    x_newA = xA
    ynew = yA
    DO j = 1, nvars
        Z0(:, j) = y0
        ZA(:, j) = yA
    END DO
    call chkvars(nobs, nvars, x_newA, juA)

    ! -------- CENTERING AND STANDARDIZATION -------- !
    CALL standard(nobs, nvars, x_newA, juA, 1, 1, xmeanA, xnormA, maj)
    CALL standard(nobs0, nvars, x_new0, ju, 1, 1, xmean0, xnorm0, maj)
    XY0 = x_new0 * Z0
    XYA = x_newA * ZA

    ! -------- COMPUTE BETA -------- !
    DO j = 1, nvars
        IF (ju(j) == 1) THEN

            beta_temp = SUM(XY0(:, j))/nobs0
            beta0_temp = SUM(Z0(:, j))/nobs0
            beta0(j) = beta0_temp
            beta(j) = beta_temp
        ELSE
            beta(j) = 0.0D0
            beta0(j) = 0.0D0
        END IF
        IF (juA(j) == 1) THEN
            beta_temp = SUM(XYA(:, j))/nobs
            beta0_temp = SUM(ZA(:, j))/nobs
            beta_mat(:, j) = beta_temp
            beta0_mat(:, j) = beta0_temp
        ELSE
            beta_mat(:, j) =0.0D0
            beta0_mat(:, j) = 0.0D0
        END IF
    END DO

    beta = beta / xnorm0
    beta0 = beta0 - beta*xmean0 / xnorm0


    ! -------- LEAVE-ONE-OUT FITTING -------- !
    IF (loo) THEN
        Ri = DBLE(nobs) * (ZA - beta0_mat - x_newA * beta_mat) / &
             (DBLE(nobs) - 1.0D0 - x_newA ** 2)

        fit_temp = ZA - Ri
        
    ELSE
        fit_temp = beta0_mat + x_newA * beta_mat
    END IF
    fit = fit_temp(1:nobs0, :)

END SUBROUTINE loofit_s

! --------------------------------------------------------------------------
SUBROUTINE loofit_st(nobs0, nvars, nobsA, x0, y0, xA, yA, beta0, beta, fit)
    ! -----------------------------------------------------------------------
    IMPLICIT NONE
    ! -------- INPUT VARIABLES -------- !
    INTEGER :: nobs0
    INTEGER :: nvars
    INTEGER :: nobsA
    DOUBLE PRECISION :: x0(nobs0, nvars)
    DOUBLE PRECISION :: y0(nobs0)
    DOUBLE PRECISION :: xA(nobsA, nvars)
    DOUBLE PRECISION :: yA(nobsA)
    ! -------- OUTPUT VARIABLES -------- !
    DOUBLE PRECISION :: beta0(nvars)
    DOUBLE PRECISION :: beta(nvars)
    DOUBLE PRECISION :: fit(nobs0, nvars)
    ! -------- LOCAL DECLARATIONS -------- !
    INTEGER :: j
    INTEGER :: ju(nvars)
    INTEGER :: juA(nvars)
    DOUBLE PRECISION :: beta0_temp
    DOUBLE PRECISION :: beta_temp
    DOUBLE PRECISION :: x_new0(nobs0, nvars)
    DOUBLE PRECISION :: x_newA(nobsA, nvars)
    DOUBLE PRECISION :: maj(nvars)
    DOUBLE PRECISION :: xnorm0(nvars)
    DOUBLE PRECISION :: xmean0(nvars)
    DOUBLE PRECISION :: xnormA(nvars)
    DOUBLE PRECISION :: xmeanA(nvars)
    DOUBLE PRECISION :: ynew(nobsA)
    DOUBLE PRECISION :: Z0(nobs0, nvars)
    DOUBLE PRECISION :: ZA(nobsA, nvars)
    DOUBLE PRECISION :: XY0(nobs0, nvars)
    DOUBLE PRECISION :: XYA(nobsA, nvars)
    DOUBLE PRECISION :: beta0_mat(nobs0, nvars)
    DOUBLE PRECISION :: beta_mat(nobs0, nvars)

    call chkvars(nobs0, nvars, x0, ju)

    ! -------- DEFINE VARIABLES -------- !
    x_new0 = x0
    x_newA = xA
    ynew = yA
    DO j = 1, nvars
        Z0(:, j) = y0
        ZA(:, j) = yA
    END DO
    call chkvars(nobsA, nvars, x_newA, juA)

    ! -------- CENTERING AND STANDARDIZATION -------- !
    CALL standard(nobsA, nvars, x_newA, juA, 1, 1, xmeanA, xnormA, maj)
    CALL standard(nobs0, nvars, x_new0, ju, 1, 1, xmean0, xnorm0, maj)
    XY0 = x_new0 * Z0
    XYA = x_newA * ZA

    ! -------- COMPUTE BETA -------- !
    DO j = 1, nvars
        IF (ju(j) == 1) THEN

            beta_temp = SUM(XY0(:, j))/nobs0
            beta0_temp = SUM(Z0(:, j))/nobs0
            beta0(j) = beta0_temp
            beta(j) = beta_temp
        ELSE
            beta(j) = 0.0D0
            beta0(j) = 0.0D0
        END IF
        IF (juA(j) == 1) THEN
            beta_temp = SUM(XYA(:, j))/nobsA
            beta0_temp = SUM(ZA(:, j))/nobsA
            beta_mat(:, j) = beta_temp
            beta0_mat(:, j) = beta0_temp
            x0(:, j) = (x0(:, j) - xmeanA(j)) / xnormA(j)
        ELSE
            beta_mat(:, j) =0.0D0
            beta0_mat(:, j) = 0.0D0
        END IF
    END DO

    beta = beta / xnorm0
    beta0 = beta0 - beta*xmean0 / xnorm0


    fit = beta0_mat + x0 * beta_mat

END SUBROUTINE loofit_st


SUBROUTINE loofit_binom(nobs, nvars, x, y, loo, beta0, beta, fit)
    ! -----------------------------------------------------------------------
    USE, INTRINSIC :: ieee_arithmetic
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
    INTEGER :: i
    INTEGER :: iter
    INTEGER :: niter
    INTEGER :: ju(nvars)
    DOUBLE PRECISION :: mus(nobs)
    DOUBLE PRECISION :: musmat(nobs, nvars)
    DOUBLE PRECISION :: etas(nobs)
    DOUBLE PRECISION :: w(nobs)
    DOUBLE PRECISION :: wmat(nobs, nvars)
    DOUBLE PRECISION :: wmatkeep(nobs, nvars)
    DOUBLE PRECISION :: wxx(nobs, nvars)
    DOUBLE PRECISION :: sws(nobs, nvars)
    DOUBLE PRECISION :: swxx(nobs, nvars)
    DOUBLE PRECISION :: ztemp(nobs)
    DOUBLE PRECISION :: beta0_temp
    DOUBLE PRECISION :: beta_temp
    DOUBLE PRECISION :: Ri(nobs, nvars)
    DOUBLE PRECISION :: x_new(nobs, nvars)
    DOUBLE PRECISION :: xm(nobs, nvars)
    DOUBLE PRECISION :: xxp(nvars)
    DOUBLE PRECISION :: maj(nvars)
    DOUBLE PRECISION :: xnorm(nvars)
    DOUBLE PRECISION :: xmean(nvars)
    DOUBLE PRECISION :: Z(nobs, nvars)
    DOUBLE PRECISION :: XY(nobs, nvars)
    DOUBLE PRECISION :: beta0_mat(nobs, nvars)
    DOUBLE PRECISION :: beta_mat(nobs, nvars)

    call chkvars(nobs, nvars, x, ju)
    

    ! -------- DEFINE VARIABLES -------- !
    niter = 2
    mus = (y + 0.5D0) / 2.0D0
    w = mus * (1.0D0 - mus)
    etas = LOG(mus / (1.0D0 - mus))
    ztemp = etas + (y - mus) / w
    DO j = 1, nvars
        Z(:, j) = ztemp
        wmatkeep(:, j) = w
    END DO
    wmat = wmatkeep / DBLE(nobs)
    x_new = x * wmat

    ! --------- START OF WLSU ---------- !
    ! CALL standard(nobs, nvars, x_new, ju, 0, 0, xmean, xnorm, maj)
    Do j = 1, nvars
        xmean(j) = SUM(x_new(:, j)) / SUM(wmat(:, 1))
    END DO
    DO j = 1, nobs
        xm(j, :) = x(j, :) - xmean
    END DO
    XY = xm * wmat * Z
    DO j = 1, nvars
        xxp(j) = SUM(wmat(:, j) * xm(:, j) * xm(:, j))
    END DO
    ! -------- COMPUTE BETA 1st -------- !
    DO j = 1, nvars
        IF (ju(j) == 1) THEN

            beta_temp = SUM(XY(:, j))/xxp(j)
            beta0_temp = SUM(wmat(:, j) * Z(:, j)) / SUM(w) / DBLE(nobs)

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

    beta0 = beta0 - beta*xmean
    beta = beta
    
    fit = beta0_mat + xm * beta_mat

    ! -------- END OF WLSU --------!

    ! -------- COMPUTE BETA 2nd -------- !
    DO iter = 1, niter
        musmat = 1.0D0 / (1.0D0 + EXP(-fit))
        wmat = musmat * (1.0D0 - musmat)
        DO j = 1, nvars
            Z(:, j) = fit(:, j) + (y - musmat(:, j)) / wmat(:, j)
        END DO
        wmat = wmat / DBLE(nobs)
        x_new = x * wmat
        ! CALL standard(nobs, nvars, x_new, ju, 0, 1, xmean, xnorm, maj)
        Do j = 1, nvars
            xmean(j) = SUM(x_new(:, j)) / SUM(wmat(:, 1))
        END DO
        DO j = 1, nobs
            xm(j, :) = x(j, :) - xmean
        END DO
        XY = xm * wmat * Z
        DO j = 1, nvars
            xxp(j) = SUM(wmat(:, j) * xm(:, j) * xm(:, j))
        END DO
        ! -------- COMPUTE BETA 1st -------- !
        DO j = 1, nvars
            IF (ju(j) == 1) THEN

                beta_temp = SUM(XY(:, j))/xxp(j)
                beta0_temp = SUM(wmat(:, j) * Z(:, j)) / SUM(w) / DBLE(nobs)

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

        beta0 = beta0 - beta*xmean
        beta = beta
        
        fit = beta0_mat + xm * beta_mat
    END DO


    ! -------- LEAVE-ONE-OUT FITTING -------- !
    IF (loo) THEN
        wmat = wmatkeep / DBLE(nobs)
        DO j = 1, nvars
            wxx(:, j) = wmat(:, j) * xm(:, j) * xm(:, j)
            ! xxp(j) = SUM(wxx(:, j)) = swxx
            sws(:, j) = SUM(wmat(:, j))
            swxx(:, j) = SUM(wxx(:, j))
        END DO
        
        Ri = wmat / sws + wxx / swxx

        Ri = (fit - Ri * Z) / (1.0D0 - Ri)
        DO j = 1, nvars
            DO i = 1, nobs
                IF (ieee_is_nan(Ri(i,j))) THEN
                    Ri(i,j) = fit(i,j)
                END IF
            END DO
        END DO
        fit = Ri
    END IF

END SUBROUTINE loofit_binom

! SUBROUTINE wlsu(nobs, nvars, x, wmat, y, beta0, beta, fit, xmean)
!     IMPLICIT NONE
!     ! -------- INPUT VARIABLES -------- !
!     INTEGER :: nobs
!     INTEGER :: nvars
!     DOUBLE PRECISION :: x(nobs, nvars)
!     DOUBLE PRECISION :: y(nobs)
!     DOUBLE PRECISION :: xmean(nvars)
!     DOUBLE PRECISION :: beta0(nvars)
!     DOUBLE PRECISION :: beta(nvars)
!     DOUBLE PRECISION :: fit(nobs, nvars)
!     DOUBLE PRECISION :: wmat(nobs, nvars)

