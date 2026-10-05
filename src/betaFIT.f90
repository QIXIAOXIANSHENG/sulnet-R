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











    ! call chkvars(nobs, nvars, x, ju)

    ! ! -------- DEFINE VARIABLES -------- !
    ! x_new = x
    ! DO j = 1, nvars
    !     Z(:, j) = y
    ! END DO

    ! ! -------- CENTERING AND STANDARDIZATION -------- !
    ! CALL standard(nobs, nvars, x_new, ju, 1, 1, xmean, xnorm, maj)
    ! XY = x_new * Z

    ! ! -------- COMPUTE BETA -------- !
    ! DO j = 1, nvars
    !     IF (ju(j) == 1) THEN

    !         beta_temp = SUM(XY(:, j))/nobs
    !         beta0_temp = SUM(Z(:, j))/nobs

    !         beta_mat(:, j) = beta_temp
    !         beta0_mat(:, j) = beta0_temp
    !         beta0(j) = beta0_temp
    !         beta(j) = beta_temp
    !     ELSE
    !         beta_mat(:, j) =0.0D0
    !         beta0_mat(:, j) = 0.0D0
    !         beta(j) = 0.0D0
    !         beta0(j) = 0.0D0
    !     END IF
    ! END DO

    ! beta0 = beta0 - beta*xmean / xnorm
    ! beta = beta / xnorm
    


    ! ! -------- LEAVE-ONE-OUT FITTING -------- !
    ! IF (loo) THEN
    !     Ri = DBLE(nobs) * (Z - beta0_mat - x_new * beta_mat) / &
    !          (DBLE(nobs) - 1.0D0 - x_new ** 2)

    !     fit = Z - Ri
    ! ELSE
    !     fit = beta0_mat + x_new * beta_mat
    ! END IF

END SUBROUTINE betaFIT