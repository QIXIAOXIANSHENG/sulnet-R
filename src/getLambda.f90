SUBROUTINE getlambda(nobs, nvars, nlam, ulam, x, y, pf, flmin)
    IMPLICIT NONE
    DOUBLE PRECISION, PARAMETER :: big = 9.9E30
    DOUBLE PRECISION, PARAMETER :: mfl = 1.0E-6

    INTEGER :: nobs
    INTEGER :: nvars
    INTEGER :: nlam
    

    DOUBLE PRECISION :: x(nobs, nvars)
    DOUBLE PRECISION :: y(nobs)
    DOUBLE PRECISION :: ulam(nlam)
    DOUBLE PRECISION :: flmin
    DOUBLE PRECISION :: pf(nvars)

    INTEGER :: l,j
    INTEGER :: ju(nvars)
    DOUBLE PRECISION :: u
    DOUBLE PRECISION:: alf
    DOUBLE PRECISION :: al
    DOUBLE PRECISION :: altemp(nlam)
    DOUBLE PRECISION :: xmean(nvars)
    DOUBLE PRECISION :: xnorm(nvars)
    DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: maj

    maj = 2.0D0*maj
    IF (flmin < 1.0D0) THEN
        flmin = MAX(mfl, flmin)
        alf = flmin**(1.0D0/(DBLE(nlam) - 1.0D0))
    END IF

    CALL standard(nobs, nvars, x, ju, 0, 1, xmean, xnorm, maj)

    DO l = 1, nlam
        ! -------- COMPUTING LAMBDA -------- !
        IF (flmin >= 1.0D0) THEN
            al = ulam(l)
            altemp(l) = al
        ELSE
            IF (l > 2) THEN
                al = al*alf
                altemp(l) = al
            ELSE IF (l == 1) THEN
                al = big
                altemp(l) = al
            ELSE IF (l == 2) THEN
                al = 0.0D0
                DO j = 1, nvars
                    IF (ju(j) /= 0) THEN
                        IF (pf(j) > 0.0D0) THEN
                            u = DOT_PRODUCT(y, x(:, j))
                            al = MAX(al, ABS(u)/pf(j))
                        END IF
                    END IF
                END DO
                al = al*alf/nobs
                altemp(l) = al
            END IF
        END IF
    END DO

    IF (flmin < 1.0D0) THEN
        ulam = altemp
        altemp = LOG(altemp)
        ulam(1) = EXP(2 * altemp(2) - altemp(3))
        ! ulam(1) = 100.0D0
    END IF
END SUBROUTINE getlambda