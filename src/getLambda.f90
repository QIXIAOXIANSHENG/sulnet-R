SUBROUTINE getlambdagauss(nobs, nvars, nlam, ulam, x, y, vp, flmin)
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
    DOUBLE PRECISION :: vp(nvars)

    INTEGER :: l, j
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
                        IF (vp(j) > 0.0D0) THEN
                            u = DOT_PRODUCT(y, x(:, j))
                            ! CALL DBLEPR("x", -1, x(1:14,j), 14)
                            ! CALL INTPR("j", -1, j, 1)
                            ! RETURN
                            ! CALL DBLEPR("u", -1, u, 1)
                            al = MAX(al, ABS(u)/vp(j))
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
        ulam(1) = EXP(2*altemp(2) - altemp(3))
        ! ulam(1) = 100.0D0
    END IF
    DEALLOCATE (maj)
END SUBROUTINE getlambdagauss

SUBROUTINE getlambdabinom(nobs, nvars, nlam, ulam, x, y, vp, flmin)
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
    DOUBLE PRECISION :: vp(nvars)

    INTEGER :: l, j
    INTEGER :: ierr
    INTEGER :: isd
    INTEGER :: intr
    INTEGER :: ju(nvars)
    DOUBLE PRECISION:: alf
    DOUBLE PRECISION :: al = 0.0D0
    DOUBLE PRECISION :: altemp(nlam)
    DOUBLE PRECISION :: xmean(nvars)
    DOUBLE PRECISION :: xnorm(nvars)
    DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: maj

    ALLOCATE (maj(1:nvars), STAT=ierr)
    maj = 0.25D0*maj
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
                        IF (vp(j) > 0.0D0) THEN
                            al = MAX(al, ABS(DOT_PRODUCT(y*0.5D0, &
                                 & x(:, j)))/vp(j))
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
        ulam(1) = EXP(2*altemp(2) - altemp(3))
        ! ulam(1) = 100.0D0
    END IF
    DEALLOCATE (maj)
END SUBROUTINE getlambdabinom


SUBROUTINE getzt(nobs, nvars, nzeta, ntau, uzeta, utau, x, y, vp, flminz, flmint)
    IMPLICIT NONE
    DOUBLE PRECISION, PARAMETER :: big = 9.9E30
    DOUBLE PRECISION, PARAMETER :: mfl = 1.0E-6

    INTEGER :: nobs
    INTEGER :: nvars
    INTEGER :: nzeta
    INTEGER :: ntau

    DOUBLE PRECISION :: x(nobs, nvars)
    DOUBLE PRECISION :: y(nobs)
    DOUBLE PRECISION :: uzeta(nzeta)
    DOUBLE PRECISION :: utau(ntau)
    DOUBLE PRECISION :: flminz
    DOUBLE PRECISION :: flmint
    DOUBLE PRECISION :: vp(nvars)
    
    INTEGER :: l, j
    INTEGER :: ierr
    INTEGER :: isd
    INTEGER :: intr
    INTEGER :: ju(nvars)
    DOUBLE PRECISION :: zf
    DOUBLE PRECISION :: tf
    DOUBLE PRECISION :: z = 0.0D0
    DOUBLE PRECISION :: t = 0.0D0
    DOUBLE PRECISION :: ztemp(nzeta)
    DOUBLE PRECISION :: ttemp(ntau)
    DOUBLE PRECISION :: u
    DOUBLE PRECISION :: pos
    DOUBLE PRECISION :: neg
    DOUBLE PRECISION :: xmean(nvars)
    DOUBLE PRECISION :: xnorm(nvars)
    DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: maj

    ALLOCATE (maj(1:nvars), STAT=ierr)
    maj = 2.0D0*maj
    IF (flminz < 1.0D0) THEN
        flminz = MAX(mfl, flminz)
        zf = flminz**(1.0D0/(DBLE(nzeta) - 1.0D0))
    END IF
    IF (flmint < 1.0D0) THEN
        flmint = MAX(mfl, flmint)
        tf = flmint**(1.0D0/(DBLE(ntau) - 1.0D0))
    END IF

    CALL chkvars(nobs, nvars, x, ju)

    isd = 0
    intr = 1
    ! CALL DBLEPR("x0", -1, x(1,:), 14)
    CALL standard(nobs, nvars, x, ju, 0, 1, xmean, xnorm, maj)
    ! CALL DBLEPR("x", -1, x(1,:), 14)
    ! CALL DBLEPR("y", -1, y(1:14), 14)
    ! CALL DBLEPR("xmean", -1, xmean, 14)
    IF (flminz < 1.0D0 .OR. flmint < 1.0D0) THEN
        ztemp = 0.0D0
        ttemp = 0.0D0
        DO j = 1, nvars
            IF (ju(j) /= 0) THEN
                IF (vp(j) > 0.0D0) THEN
                    u = DOT_PRODUCT(y, x(:,j))
                    IF (u > 0.0D0) THEN
                        pos = u/vp(j)
                        z = MAX(z, pos)
                    ELSE IF (u < 0.0D0) THEN
                        neg = ABS(u)/vp(j)
                        t = MAX(t, neg)
                    END IF
                END IF
            END IF
        END DO
        ! z = z*zf/nobs
        ! t = t*tf/nobs
    END IF

    IF (flminz >= 1.0D0) THEN
    ELSE
        DO l = 1, nzeta
            IF(l>2) THEN
                pos = pos*zf
                ztemp(l) = pos
            ELSE IF (l == 1) THEN
                ztemp(l) = big
            ELSE IF (l == 2) THEN
                pos = MAX(z, t)
                pos = pos*zf/nobs
                ztemp(l) = pos
            END IF
        END DO
        uzeta = ztemp
        ztemp = LOG(ztemp)
        uzeta(1) = EXP(2*ztemp(2) - ztemp(3))
    END IF

    IF (flmint >= 1.0D0) THEN
    ELSE
        DO l = 1, ntau
            IF(l>2) THEN
                neg = neg*tf
                ttemp(l) = neg
            ELSE IF (l == 1) THEN
                ttemp(l) = big
            ELSE IF (l == 2) THEN
                neg = ABS(z-t)
                neg = neg*tf/nobs
                ttemp(l) = neg
            END IF
        END DO
        utau = ttemp
        ttemp = LOG(ttemp)
        utau(1) = EXP(2*ttemp(2) - ttemp(3))
    END IF
    DEALLOCATE (maj)
END SUBROUTINE getzt


