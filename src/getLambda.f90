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


SUBROUTINE getzt(nobs, nvars, nzeta, ntau, uzeta, utau, x, y, vp, flminz, flmint, isd, intr, diffz, difft)
    IMPLICIT NONE
    DOUBLE PRECISION, PARAMETER :: mfl = 1.0E-6

    INTEGER :: nobs
    INTEGER :: nvars
    INTEGER :: nzeta
    INTEGER :: ntau
    INTEGER :: isd
    INTEGER :: intr

    DOUBLE PRECISION :: x(nobs, nvars)
    DOUBLE PRECISION :: y(nobs)
    DOUBLE PRECISION :: uzeta(nzeta)
    DOUBLE PRECISION :: utau(ntau)
    DOUBLE PRECISION :: flminz
    DOUBLE PRECISION :: flmint
    DOUBLE PRECISION :: diffz
    DOUBLE PRECISION :: difft
    DOUBLE PRECISION :: vp(nvars)
    
    INTEGER :: l, j
    INTEGER :: ierr

    INTEGER :: ju(nvars)
    DOUBLE PRECISION :: zf
    DOUBLE PRECISION :: tf
    DOUBLE PRECISION :: z
    DOUBLE PRECISION :: t
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
        zf = 1.0D0
        IF (nzeta > 1) zf = flminz**(1.0D0/(DBLE(nzeta) - 1.0D0))
    ELSE
        diffz = -1.0D0
    END IF
    IF (flmint < 1.0D0) THEN
        flmint = MAX(mfl, flmint)
        tf = 1.0D0
        IF (ntau > 1) tf = flmint**(1.0D0/(DBLE(ntau) - 1.0D0))
    ELSE
        difft = -1.0D0
    END IF

    CALL chkvars(nobs, nvars, x, ju)

    ! isd = 0
    ! intr = 1
    ! CALL DBLEPR("x0", -1, x(1,:), 14)
    CALL standard(nobs, nvars, x, ju, isd, intr, xmean, xnorm, maj)
    ! CALL DBLEPR("x", -1, x(1,:), 14)
    ! CALL DBLEPR("y", -1, y(1:14), 14)
    ! CALL DBLEPR("xmean", -1, xmean, 14)
    IF (flminz < 1.0D0 .OR. flmint < 1.0D0) THEN
        z = 0.0D0
        t = 0.0D0
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

    ! Geometric grids starting at the smallest all-zero zeta. The first zeta
    ! is raised by a relative 1e-9 so that rounding in the solver cannot leave
    ! a tiny nonzero coefficient in the first row. Grids of any length >= 1.
    IF (flminz < 1.0D0) THEN
        pos = MAX(z, t)/DBLE(nobs)
        DO l = 1, nzeta
            uzeta(l) = pos*zf**(l - 1)
        END DO
        uzeta(1) = uzeta(1)*(1.0D0 + 1.0D-9)
        diffz = -LOG(zf)
    END IF

    ! The tau grid starts at |z - t|, or at the zeta scale when the largest
    ! positive and negative scores (nearly) coincide.
    IF (flmint < 1.0D0) THEN
        neg = ABS(z - t)
        IF (neg <= 1.0D-8*MAX(z, t)) neg = MAX(z, t)
        neg = neg/DBLE(nobs)
        DO l = 1, ntau
            utau(l) = neg*tf**(l - 1)
        END DO
        difft = -LOG(tf)
    END IF
    DEALLOCATE (maj)
END SUBROUTINE getzt

SUBROUTINE getla(nobs, nvars, nlambda, nalpha, ulambda, ualpha,&
    & x, y, vp, flmin, flmina, isd, intr)
    IMPLICIT NONE
    DOUBLE PRECISION, PARAMETER :: mfl = 1.0E-6

    INTEGER :: i
    INTEGER :: nobs
    INTEGER :: nvars
    INTEGER :: nlambda
    INTEGER :: nalpha
    INTEGER :: isd
    INTEGER :: intr

    DOUBLE PRECISION :: x(nobs, nvars)
    DOUBLE PRECISION :: y(nobs)
    DOUBLE PRECISION :: ulambda(nlambda)
    DOUBLE PRECISION :: ualpha(nalpha)
    DOUBLE PRECISION :: flmin
    DOUBLE PRECISION :: flmina
    DOUBLE PRECISION :: vp(nvars)

    DOUBLE PRECISION :: diffz
    DOUBLE PRECISION :: difft
    DOUBLE PRECISION :: maxa
    DOUBLE PRECISION :: mina
    DOUBLE PRECISION :: max_l
    DOUBLE PRECISION :: max_a
    DOUBLE PRECISION :: min_a
    DOUBLE PRECISION :: min_l


    diffz = 0.0D0
    difft = 0.0D0

    CALL getzt(nobs, nvars, nlambda, nalpha, ulambda, ualpha &
    &, x, y, vp, flmin, flmina, isd, intr, diffz, difft)

    
    IF(flmina < 1.0D0) THEN
        max_l = MAXVAL(ulambda)
        max_a = MAXVAL(ualpha)
        min_a = MINVAL(ualpha)
        min_l = MINVAL(ulambda)
        maxa = MAX(0.5D0, min_a / (max_l + max_a))
        mina = MIN(1.0D0-mfl, &
            & MAX(max_l, max_a)/(MAX(max_l, max_a)+MIN(min_l, min_a)))
        IF (nalpha == 1) THEN
            ualpha(1) = mina
        ELSE
            DO i = 1, nalpha
                ualpha(i) = mina + (maxa - mina) * DBLE(i - 1) / DBLE(nalpha - 1)
            END DO
            ! exact endpoint, so the grid can be passed back in (alpha >= 0.5)
            ualpha(nalpha) = maxa
        END IF
    END IF
    IF(flmin < 1.0D0) THEN
        ulambda = ulambda * 2.0D0
    END IF
    
    
END SUBROUTINE getla


SUBROUTINE getsr(nobs, nvars, nsum, nrate, usum, urate,&
    & x, y, vp, flmins, flminr, isd, intr)
    IMPLICIT NONE
    DOUBLE PRECISION, PARAMETER :: mfl = 1.0E-6

    INTEGER :: i
    INTEGER :: nobs
    INTEGER :: nvars
    INTEGER :: nsum
    INTEGER :: nrate
    INTEGER :: isd
    INTEGER :: intr

    DOUBLE PRECISION :: x(nobs, nvars)
    DOUBLE PRECISION :: y(nobs)
    DOUBLE PRECISION :: usum(nsum)
    DOUBLE PRECISION :: urate(nrate)
    DOUBLE PRECISION :: flmins
    DOUBLE PRECISION :: flminr
    DOUBLE PRECISION :: vp(nvars)

    DOUBLE PRECISION :: diffz
    DOUBLE PRECISION :: difft
    DOUBLE PRECISION :: maxr
    DOUBLE PRECISION :: minr
    DOUBLE PRECISION :: max_s
    DOUBLE PRECISION :: max_r
    DOUBLE PRECISION :: min_s
    DOUBLE PRECISION :: min_r
    DOUBLE PRECISION :: step


    diffz = 0.0D0
    difft = 0.0D0

    CALL getzt(nobs, nvars, nsum, nrate, usum, urate &
    &, x, y, vp, flmins, flminr, isd, intr, diffz, difft)

    
    IF(flminr < 1.0D0) THEN
        max_s = MAXVAL(usum)
        max_r = MAXVAL(urate)
        min_r = MINVAL(urate)
        min_s = MINVAL(usum)
        maxr = MAX(0.5D0, min_r / (max_s + max_r))
        minr = MIN(1.0D0-mfl, &
            & MAX(max_s, max_r)/(MAX(max_s, max_r)+MIN(min_s, min_r)))
        maxr = maxr / (1 - maxr)
        minr = minr / (1 - minr)
        IF (nrate == 1) THEN
            urate(1) = minr
        ELSE
            step = (LOG(maxr) - LOG(minr)) / DBLE(nrate-1)
            DO i = 1, nrate
                urate(i) = minr * EXP(DBLE(i-1) * step)
            END DO
            ! exact endpoint, so the grid can be passed back in (rate >= 1)
            urate(nrate) = maxr
        END IF
    END IF
    IF(flmins < 1.0D0) THEN
        usum = usum * 2.0D0
    END IF
    
    
END SUBROUTINE getsr


SUBROUTINE getztlog(nobs, nvars, nzeta, ntau, uzeta, utau, x, y, vp, flminz, flmint, isd, intr, diffz, difft)
    IMPLICIT NONE
    DOUBLE PRECISION, PARAMETER :: mfl = 1.0E-6

    INTEGER :: nobs
    INTEGER :: nvars
    INTEGER :: nzeta
    INTEGER :: ntau
    INTEGER :: isd
    INTEGER :: intr

    DOUBLE PRECISION :: x(nobs, nvars)
    DOUBLE PRECISION :: y(nobs)
    DOUBLE PRECISION :: uzeta(nzeta)
    DOUBLE PRECISION :: utau(ntau)
    DOUBLE PRECISION :: flminz
    DOUBLE PRECISION :: flmint
    DOUBLE PRECISION :: diffz
    DOUBLE PRECISION :: difft
    DOUBLE PRECISION :: vp(nvars)
    
    INTEGER :: l, j
    INTEGER :: ierr

    INTEGER :: ju(nvars)
    DOUBLE PRECISION :: zf
    DOUBLE PRECISION :: tf
    DOUBLE PRECISION :: z
    DOUBLE PRECISION :: t
    DOUBLE PRECISION :: u
    DOUBLE PRECISION :: pos
    DOUBLE PRECISION :: neg
    DOUBLE PRECISION :: xmean(nvars)
    DOUBLE PRECISION :: xnorm(nvars)
    DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: maj

    ALLOCATE (maj(1:nvars), STAT=ierr)
    maj = 0.25D0*maj
    IF (flminz < 1.0D0) THEN
        flminz = MAX(mfl, flminz)
        zf = 1.0D0
        IF (nzeta > 1) zf = flminz**(1.0D0/(DBLE(nzeta) - 1.0D0))
    ELSE
        diffz = -1.0D0
    END IF
    IF (flmint < 1.0D0) THEN
        flmint = MAX(mfl, flmint)
        tf = 1.0D0
        IF (ntau > 1) tf = flmint**(1.0D0/(DBLE(ntau) - 1.0D0))
    ELSE
        difft = -1.0D0
    END IF

    CALL chkvars(nobs, nvars, x, ju)

    ! isd = 0
    ! intr = 1
    ! CALL DBLEPR("x0", -1, x(1,:), 14)
    CALL standard(nobs, nvars, x, ju, isd, intr, xmean, xnorm, maj)
    ! CALL DBLEPR("x", -1, x(1,:), 14)
    ! CALL DBLEPR("y", -1, y(1:14), 14)
    ! CALL DBLEPR("xmean", -1, xmean, 14)
    IF (flminz < 1.0D0 .OR. flmint < 1.0D0) THEN
        z = 0.0D0
        t = 0.0D0
        DO j = 1, nvars
            IF (ju(j) /= 0) THEN
                IF (vp(j) > 0.0D0) THEN
                    u = DOT_PRODUCT(y*0.5D0, x(:, j))
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

    ! Same grids as getzt: see the comments there.
    IF (flminz < 1.0D0) THEN
        pos = MAX(z, t)/DBLE(nobs)
        DO l = 1, nzeta
            uzeta(l) = pos*zf**(l - 1)
        END DO
        uzeta(1) = uzeta(1)*(1.0D0 + 1.0D-9)
        diffz = -LOG(zf)
    END IF

    IF (flmint < 1.0D0) THEN
        neg = ABS(z - t)
        IF (neg <= 1.0D-8*MAX(z, t)) neg = MAX(z, t)
        neg = neg/DBLE(nobs)
        DO l = 1, ntau
            utau(l) = neg*tf**(l - 1)
        END DO
        difft = -LOG(tf)
    END IF
    DEALLOCATE (maj)
END SUBROUTINE getztlog

SUBROUTINE getlalog(nobs, nvars, nlambda, nalpha, ulambda, ualpha,&
    & x, y, vp, flmin, flmina, isd, intr)
    IMPLICIT NONE
    DOUBLE PRECISION, PARAMETER :: mfl = 1.0E-6

    INTEGER :: i
    INTEGER :: nobs
    INTEGER :: nvars
    INTEGER :: nlambda
    INTEGER :: nalpha
    INTEGER :: isd
    INTEGER :: intr

    DOUBLE PRECISION :: x(nobs, nvars)
    DOUBLE PRECISION :: y(nobs)
    DOUBLE PRECISION :: ulambda(nlambda)
    DOUBLE PRECISION :: ualpha(nalpha)
    DOUBLE PRECISION :: flmin
    DOUBLE PRECISION :: flmina
    DOUBLE PRECISION :: vp(nvars)

    DOUBLE PRECISION :: diffz
    DOUBLE PRECISION :: difft
    DOUBLE PRECISION :: maxa
    DOUBLE PRECISION :: mina
    DOUBLE PRECISION :: max_l
    DOUBLE PRECISION :: max_a
    DOUBLE PRECISION :: min_a
    DOUBLE PRECISION :: min_l


    diffz = 0.0D0
    difft = 0.0D0

    CALL getztlog(nobs, nvars, nlambda, nalpha, ulambda, ualpha &
    &, x, y, vp, flmin, flmina, isd, intr, diffz, difft)

    
    IF(flmina < 1.0D0) THEN
        max_l = MAXVAL(ulambda)
        max_a = MAXVAL(ualpha)
        min_a = MINVAL(ualpha)
        min_l = MINVAL(ulambda)
        maxa = MAX(0.5D0, min_a / (max_l + max_a))
        mina = MIN(1.0D0-mfl, &
            & MAX(max_l, max_a)/(MAX(max_l, max_a)+MIN(min_l, min_a)))
        IF (nalpha == 1) THEN
            ualpha(1) = mina
        ELSE
            DO i = 1, nalpha
                ualpha(i) = mina + (maxa - mina) * DBLE(i - 1) / DBLE(nalpha - 1)
            END DO
            ! exact endpoint, so the grid can be passed back in (alpha >= 0.5)
            ualpha(nalpha) = maxa
        END IF
    END IF
    IF(flmin < 1.0D0) THEN
        ulambda = ulambda * 2.0D0
    END IF
    
    
END SUBROUTINE getlalog

SUBROUTINE getsrlog(nobs, nvars, nsum, nrate, usum, urate,&
    & x, y, vp, flmins, flminr, isd, intr)
    IMPLICIT NONE
    DOUBLE PRECISION, PARAMETER :: mfl = 1.0E-6

    INTEGER :: i
    INTEGER :: nobs
    INTEGER :: nvars
    INTEGER :: nsum
    INTEGER :: nrate
    INTEGER :: isd
    INTEGER :: intr

    DOUBLE PRECISION :: x(nobs, nvars)
    DOUBLE PRECISION :: y(nobs)
    DOUBLE PRECISION :: usum(nsum)
    DOUBLE PRECISION :: urate(nrate)
    DOUBLE PRECISION :: flmins
    DOUBLE PRECISION :: flminr
    DOUBLE PRECISION :: vp(nvars)

    DOUBLE PRECISION :: diffz
    DOUBLE PRECISION :: difft
    DOUBLE PRECISION :: maxr
    DOUBLE PRECISION :: minr
    DOUBLE PRECISION :: max_s
    DOUBLE PRECISION :: max_r
    DOUBLE PRECISION :: min_s
    DOUBLE PRECISION :: min_r
    DOUBLE PRECISION :: step


    diffz = 0.0D0
    difft = 0.0D0

    CALL getztlog(nobs, nvars, nsum, nrate, usum, urate &
    &, x, y, vp, flmins, flminr, isd, intr, diffz, difft)

    
    IF(flminr < 1.0D0) THEN
        max_s = MAXVAL(usum)
        max_r = MAXVAL(urate)
        min_r = MINVAL(urate)
        min_s = MINVAL(usum)
        maxr = MAX(0.5D0, min_r / (max_s + max_r))
        minr = MIN(1.0D0-mfl, &
            & MAX(max_s, max_r)/(MAX(max_s, max_r)+MIN(min_s, min_r)))
        maxr = maxr / (1 - maxr)
        minr = minr / (1 - minr)
        IF (nrate == 1) THEN
            urate(1) = minr
        ELSE
            step = (LOG(maxr) - LOG(minr)) / DBLE(nrate-1)
            DO i = 1, nrate
                urate(i) = minr * EXP(DBLE(i-1) * step)
            END DO
            ! exact endpoint, so the grid can be passed back in (rate >= 1)
            urate(nrate) = maxr
        END IF
    END IF
    IF(flmins < 1.0D0) THEN
        usum = usum * 2.0D0
    END IF
    
    
END SUBROUTINE getsrlog