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
    INTEGER :: ierr
    INTEGER :: isd
    INTEGER :: intr
    INTEGER :: ju(nvars)
    DOUBLE PRECISION :: u
    DOUBLE PRECISION:: alf
    DOUBLE PRECISION :: al = 0.0D0
    DOUBLE PRECISION :: altemp(nlam)
    DOUBLE PRECISION :: xmean(nvars)
    DOUBLE PRECISION :: xnorm(nvars)
    DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: maj

    ALLOCATE (maj(1:nvars), STAT=ierr)
    maj = 2.0D0*maj
    IF (flmin < 1.0D0) THEN
        flmin = MAX(mfl, flmin)
        alf = flmin**(1.0D0/(DBLE(nlam) - 1.0D0))
    END IF

    CALL chkvars(nobs, nvars, x, ju)

    isd = 0
    intr = 1
    ! CALL DBLEPR("x0", -1, x(1,:), 14)
    CALL standard(nobs, nvars, x, ju, 0, 1, xmean, xnorm, maj)
    ! CALL DBLEPR("x", -1, x(1,:), 14)
    ! CALL DBLEPR("y", -1, y(1:14), 14)
    ! CALL DBLEPR("xmean", -1, xmean, 14)

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
                ! CALL DBLEPR("altemp", -1, altemp(1), 1)
            ELSE IF (l == 2) THEN
                al = 0.0D0
                DO j = 1, nvars
                    IF (ju(j) /= 0) THEN
                        IF (pf(j) > 0.0D0) THEN
                            u = DOT_PRODUCT(y, x(:, j))
                            ! CALL DBLEPR("x", -1, x(1:14,j), 14)
                            ! CALL INTPR("j", -1, j, 1)
                            ! RETURN
                            ! CALL DBLEPR("u", -1, u, 1)
                            al = MAX(al, ABS(u)/pf(j))
                            ! CALL DBLEPR("al", -1, al, 1)
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
    DEALLOCATE (maj)
END SUBROUTINE getlambda