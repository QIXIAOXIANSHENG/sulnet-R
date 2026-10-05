SUBROUTINE suniZTNET(lam2, nobs, nvars, x, y, jd, pf, pf2, pf3, pmax, nzeta,&
& ntau, uzeta, utau, eps, isd, intr, maxit, napenalty, theta0, &
& theta, itheta, ntheta, apenalty, npass, diffzeta, difftau, jerr)
    IMPLICIT NONE

    INTEGER :: nobs
    INTEGER :: nvars
    INTEGER :: pmax
    INTEGER :: nzeta
    INTEGER :: ntau
    INTEGER :: isd
    INTEGER :: intr
    INTEGER :: napenalty
    INTEGER :: npass
    INTEGER :: jerr
    INTEGER :: maxit
    INTEGER :: jd(*)
    INTEGER :: itheta(pmax, nzeta, ntau)
    INTEGER :: ntheta(nzeta, ntau)
    DOUBLE PRECISION :: lam2
    DOUBLE PRECISION :: eps
    DOUBLE PRECISION :: diffzeta
    DOUBLE PRECISION :: difftau
    DOUBLE PRECISION :: x(nobs, nvars)
    DOUBLE PRECISION :: y(nobs)
    DOUBLE PRECISION :: pf(nvars)
    DOUBLE PRECISION :: pf2(nvars)
    DOUBLE PRECISION :: pf3(nvars)
    DOUBLE PRECISION :: uzeta(nzeta)
    DOUBLE PRECISION :: utau(ntau)
    DOUBLE PRECISION :: theta0(1, nzeta, ntau)
    DOUBLE PRECISION :: theta(pmax, nzeta, ntau)
    DOUBLE PRECISION :: apenalty(nzeta * ntau, 2)

    INTEGER :: j
    INTEGER :: l
    INTEGER :: nk
    INTEGER :: iz
    INTEGER :: it
    INTEGER :: ierr
    INTEGER, DIMENSION(:), ALLOCATABLE :: ju
    INTEGER, DIMENSION(:,:), ALLOCATABLE :: ind
    DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: xmean
    DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: xnorm
    DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: maj

    ALLOCATE (ju(1:nvars), STAT=ierr)
    jerr = jerr + ierr
    ALLOCATE (xmean(1:nvars), STAT=ierr)
    jerr = jerr + ierr
    ALLOCATE (maj(1:nvars), STAT=ierr)
    jerr = jerr + ierr
    ALLOCATE (xnorm(1:nvars), STAT=ierr)
    jerr = jerr + ierr
    ALLOCATE (ind(1:(nzeta * ntau), 1:4), STAT=ierr)
    jerr = jerr + ierr
    IF (jerr /= 0) RETURN
    CALL chkvars(nobs, nvars, x, ju)
    IF (jd(1) > 0) ju(jd(2:(jd(1) + 1))) = 0
    IF (MAXVAL(ju) <= 0) THEN
        jerr = 7777
        RETURN
    END IF
    IF (MAXVAL(pf) <= 0.0D0) THEN
        jerr = 10000
        RETURN
    END IF
    IF (MAXVAL(pf2) <= 0.0D0) THEN
        jerr = 10000
        RETURN
    END IF
    pf = MAX(0.0D0, pf)
    pf2 = MAX(0.0D0, pf2)
    pf3 = MAX(0.0D0, pf3)
    CALL standard(nobs, nvars, x, ju, isd, intr, xmean, xnorm, maj)
    CALL gridIndexZT(nzeta, ntau, uzeta, utau, ind, diffzeta, difftau)

    CALL suniZTNETpath(lam2, maj, nobs, nvars, x, y, ju, pf, pf2, pf3, &
        & pmax, nzeta, ntau, uzeta, utau, eps, maxit, napenalty, theta0, &
        & theta, itheta, ntheta, apenalty, npass, ind, intr, jerr)
    IF (jerr > 0) RETURN

    ! -------- ORGANIZE THETA AFTERWARDS -------- !
    DO l = 1, napenalty
        iz = ind(l, 1)
        it = ind(l, 2)
        nk = ntheta(iz, it)
        IF (isd == 1 .AND. nk > 0) THEN
            DO j = 1, nk
                theta(j, iz, it) = theta(j, iz, it)/xnorm(itheta(j, iz, it))
            END DO
        END IF
        IF (nk>0) THEN
            theta0(1, iz, it) = theta0(1, iz, it) - &
            & DOT_PRODUCT(theta(1:nk, iz, it), xmean(itheta(1:nk, iz, it)))
        END IF
    END DO
    DEALLOCATE (ju, xmean, maj, xnorm, ind)
    RETURN
END SUBROUTINE suniZTNET


SUBROUTINE suniZTNETpath(lam2, maj, nobs, nvars, x, y, ju, pf, pf2, pf3, &
    & pmax, nzeta, ntau, uzeta, utau, eps, maxit, napenalty, theta0, theta, &
    & itheta, ntheta, apenalty, npass, ind, intr, jerr)

    IMPLICIT NONE

    INTEGER :: nobs
    INTEGER :: nvars
    INTEGER :: pmax
    INTEGER :: nzeta
    INTEGER :: ntau
    INTEGER :: maxit
    INTEGER :: napenalty
    INTEGER :: npass
    INTEGER :: intr
    INTEGER :: jerr
    INTEGER :: ju(nvars)
    INTEGER :: itheta(pmax, nzeta, ntau)
    INTEGER :: ntheta(nzeta, ntau)
    INTEGER :: ind(nzeta*ntau, 4)

    DOUBLE PRECISION :: lam2
    DOUBLE PRECISION :: maj(nvars)
    DOUBLE PRECISION :: x(nobs, nvars)
    DOUBLE PRECISION :: y(nobs)
    DOUBLE PRECISION :: pf(nvars)
    DOUBLE PRECISION :: pf2(nvars)
    DOUBLE PRECISION :: pf3(nvars)
    DOUBLE PRECISION :: uzeta(nzeta)
    DOUBLE PRECISION :: utau(ntau)
    DOUBLE PRECISION :: eps
    DOUBLE PRECISION :: theta0(1, nzeta, ntau)
    DOUBLE PRECISION :: theta(pmax, nzeta, ntau)
    DOUBLE PRECISION :: apenalty(nzeta*ntau, 2)

    DOUBLE PRECISION :: d
    DOUBLE PRECISION :: dif
    DOUBLE PRECISION :: oldt
    DOUBLE PRECISION :: u
    DOUBLE PRECISION :: zeta
    DOUBLE PRECISION :: tau
    DOUBLE PRECISION :: negpenalty
    DOUBLE PRECISION :: pospenalty
    DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: t
    DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: oldtheta
    DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: r

    INTEGER :: k
    INTEGER :: j
    INTEGER :: l
    INTEGER :: vrg
    INTEGER :: ierr
    INTEGER :: ni
    INTEGER :: iz
    INTEGER :: it
    INTEGER :: prevz
    INTEGER :: prevt
    INTEGER :: npenalty
    INTEGER, DIMENSION(:), ALLOCATABLE :: mm

    npenalty = nzeta * ntau

    ALLOCATE (t(0:nvars), STAT=jerr)
    ALLOCATE (oldtheta(0:nvars), STAT = ierr)
    jerr = jerr + ierr
    ALLOCATE (mm(1:nvars), STAT=ierr)
    jerr = jerr + ierr
    ALLOCATE (r(1:nobs), STAT=ierr)
    jerr = jerr + ierr
    IF (jerr /= 0) RETURN

    maj = 2.0D0*maj

    DO l = 1, npenalty
        iz = ind(l, 1)
        it = ind(l, 2)
        prevz = ind(l, 3)
        prevt = ind(l, 4)
        zeta = uzeta(iz)
        tau = utau(it)
        t = 0.0D0
        oldtheta = 0.0D0
        mm = 0
        r = y
        ni = 0
        IF(l > 1) THEN
            ni = ntheta(prevz, prevt)
            t(0) = theta0(1, prevz, prevt)
            r = r - t(0)
            DO j = 1, ni
                k = itheta(j, prevz, prevt)
                itheta(j, iz, it) = k
                mm(k) = j
                t(k) = theta(j, prevz, prevt)
                r = r - x(:, k)*t(k)
            END DO
        END IF
        ! -------- OUTER LOOP -------- !
        DO
            IF (intr == 1) oldtheta(0) = t(0)
            IF (ni > 0) oldtheta(itheta(1:ni, iz, it)) = t(itheta(1:ni, iz, it))
            ! -------- MIDDLE LOOP -------- !
            DO
                npass = npass + 1
                dif = 0.0D0
                DO k = 1, nvars
                    IF (ju(k) /= 0) THEN
                         oldt = t(k)
                        u = DOT_PRODUCT(r, x(:, k))
                        u = maj(k) * t(k) + u/DBLE(nobs)
                        pospenalty = zeta * pf(k)
                        negpenalty = tau * pf3(k) + pospenalty
                        IF (u > pospenalty) THEN
                            t(k) = (u - pospenalty) / &
                                & (maj(k) + pf2(k)*lam2)
                        ELSE IF (u < -negpenalty) THEN
                            t(k) = (u + negpenalty) / &
                                & (maj(k) + pf2(k)*lam2)
                        ELSE
                            t(k) = 0.0D0
                        END IF
                        d = t(k) - oldt
                        IF (ABS(d) > 0.0D0) THEN
                            dif = MAX(dif, d**2)
                            r = r - x(:, k)*d
                            IF (mm(k) == 0) THEN
                                ni = ni + 1
                                IF (ni > pmax) EXIT
                                mm(k) = ni
                                itheta(ni, iz, it) = k
                            END IF
                        END IF
                    END IF
                END DO
                IF (ni > pmax) EXIT
                IF (intr == 1) THEN
                    d = SUM(r)/DBLE(nobs)
                    IF (d /= 0.0D0) THEN
                        t(0) = t(0) + d
                        r = r - d
                        dif = MAX(dif, d**2)
                    END IF
                END IF
                IF (dif < eps) EXIT
                IF (npass > maxit) THEN
                    jerr = -l
                    RETURN
                END IF
                ! -------- INNER LOOP -------- !
                DO
                    npass = npass + 1
                    dif = 0.0D0
                    DO j = 1, ni
                        k = itheta(j, iz, it)
                        oldt = t(k)
                        u = DOT_PRODUCT(r, x(:, k))
                        u = maj(k) * t(k) + u/DBLE(nobs)
                        pospenalty = zeta * pf(k)
                        negpenalty = tau * pf3(k) + pospenalty
                        IF (u > pospenalty) THEN
                            t(k) = (u - pospenalty)/(maj(k) + pf2(k)*lam2)
                        ELSE IF (u < -negpenalty) THEN
                            t(k) = (u + negpenalty)/(maj(k) + pf2(k)*lam2)
                        ELSE
                            t(k) = 0.0D0
                        END IF
                        d = t(k) - oldt
                        IF (ABS(d) > 0.0D0) THEN
                            dif = MAX(dif, d**2)
                            r = r - x(:, k)*d
                        END IF
                    END DO
                    IF (intr == 1) THEN
                        d = SUM(r)/DBLE(nobs)
                        IF (d /= 0.0D0) THEN
                            t(0) = t(0) + d
                            r = r - d
                            dif = MAX(dif, d**2)
                        END IF
                    END IF
                    IF (dif < eps) EXIT
                    IF (npass > maxit) THEN
                        jerr = -l
                        RETURN
                    END IF
                END DO
            END DO
            IF (ni > pmax) EXIT
            !-------- FINAL CHECK -------- !
            vrg = 1
            IF ((t(0) - oldtheta(0))**2 >= eps) vrg = 0
            DO j = 1, ni
                IF ((t(itheta(j, iz, it)) - oldtheta(itheta(j, iz, it)))**2 >= eps) THEN
                    vrg = 0
                    EXIT
                END IF
            END DO
            IF (vrg == 1) EXIT
        END DO
        ! -------- FINAL UPDATE & SAVE RESULTS -------- !
        IF (ni > pmax) THEN
            jerr = -10000 - l
            EXIT
        END IF
        IF (ni > 0) theta(1:ni, iz, it) = t(itheta(1:ni, iz, it))
        ntheta(iz, it) = ni
        theta0(1, iz, it) = t(0)
        apenalty(l, 1) = zeta
        apenalty(l, 2) = tau
        napenalty = l
    END DO
    DEALLOCATE (t, oldtheta, mm, r)
    RETURN
END SUBROUTINE suniZTNETpath

            