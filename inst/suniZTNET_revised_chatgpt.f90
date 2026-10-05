SUBROUTINE suniZTNET(lam2, nobs, nvars, x, y, jd, pf, pf2, dfmax, pmax, nzeta, &
    & ntau, uzeta, utau, eps, isd, intr, maxit, napenalty, theta0, theta, &
    & itheta, ntheta, apenalty, npass, diffzeta, difftau, jerr)

    IMPLICIT NONE

    INTEGER :: nobs
    INTEGER :: nvars
    INTEGER :: dfmax
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
    DOUBLE PRECISION :: uzeta(nzeta)
    DOUBLE PRECISION :: utau(ntau)
    DOUBLE PRECISION :: theta0(1, nzeta, ntau)
    DOUBLE PRECISION :: theta(pmax, nzeta, ntau)
    DOUBLE PRECISION :: apenalty(nzeta*ntau, 2)

    INTEGER :: j
    INTEGER :: l
    INTEGER :: nk
    INTEGER :: ierr
    INTEGER, DIMENSION(:), ALLOCATABLE :: ju
    INTEGER, DIMENSION(:,:), ALLOCATABLE :: ind

    DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: xmean
    DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: xnorm
    DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: maj

    jerr = 0
    napenalty = 0
    npass = 0
    theta0 = 0.0D0
    theta = 0.0D0
    itheta = 0
    ntheta = 0
    apenalty = 0.0D0

    ALLOCATE (ju(1:nvars), STAT=ierr)
    jerr = jerr + ierr
    ALLOCATE (xmean(1:nvars), STAT=ierr)
    jerr = jerr + ierr
    ALLOCATE (maj(1:nvars), STAT=ierr)
    jerr = jerr + ierr
    ALLOCATE (xnorm(1:nvars), STAT=ierr)
    jerr = jerr + ierr
    ALLOCATE (ind(1:(nzeta*ntau), 1:4), STAT=ierr)
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

    CALL standard(nobs, nvars, x, ju, isd, intr, xmean, xnorm, maj)
    CALL gridIndexZT(nzeta, ntau, uzeta, utau, ind, diffzeta, difftau)

    CALL suniZTNETpath(lam2, maj, nobs, nvars, x, y, ju, pf, pf2, dfmax, &
        & pmax, nzeta, ntau, uzeta, utau, eps, maxit, napenalty, theta0, &
        & theta, itheta, ntheta, apenalty, npass, ind, intr, jerr)

    IF (jerr > 0) RETURN

    DO l = 1, napenalty
        nk = ntheta(ind(l, 1), ind(l, 2))

        IF (isd == 1 .AND. nk > 0) THEN
            DO j = 1, nk
                theta(j, ind(l, 1), ind(l, 2)) = &
                    & theta(j, ind(l, 1), ind(l, 2)) / &
                    & xnorm(itheta(j, ind(l, 1), ind(l, 2)))
            END DO
        END IF

        IF (nk > 0) THEN
            theta0(1, ind(l, 1), ind(l, 2)) = &
                & theta0(1, ind(l, 1), ind(l, 2)) - &
                & DOT_PRODUCT(theta(1:nk, ind(l, 1), ind(l, 2)), &
                & xmean(itheta(1:nk, ind(l, 1), ind(l, 2))))
        END IF
    END DO

    DEALLOCATE (ju, xmean, maj, xnorm, ind)
    RETURN

END SUBROUTINE suniZTNET


SUBROUTINE suniZTNETpath(lam2, maj, nobs, nvars, x, y, ju, pf, pf2, dfmax, &
    & pmax, nzeta, ntau, uzeta, utau, eps, maxit, napenalty, theta0, theta, &
    & itheta, ntheta, apenalty, npass, ind, intr, jerr)

    IMPLICIT NONE

    INTEGER :: nobs
    INTEGER :: nvars
    INTEGER :: dfmax
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
    DOUBLE PRECISION :: negpenalty
    DOUBLE PRECISION :: pospenalty
    DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: t
    DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: oldtheta
    DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: r
    DOUBLE PRECISION, DIMENSION(:,:,:), ALLOCATABLE :: rmat

    INTEGER :: k
    INTEGER :: j
    INTEGER :: l
    INTEGER :: izeta
    INTEGER :: itau
    INTEGER :: prevz
    INTEGER :: prevt
    INTEGER :: vrg
    INTEGER :: ierr
    INTEGER :: ni
    INTEGER :: npenalty
    INTEGER, DIMENSION(:), ALLOCATABLE :: m
    INTEGER, DIMENSION(:), ALLOCATABLE :: mm

    npenalty = nzeta*ntau

    ALLOCATE (t(0:nvars), STAT=ierr)
    jerr = jerr + ierr
    ALLOCATE (oldtheta(0:nvars), STAT=ierr)
    jerr = jerr + ierr
    ALLOCATE (m(1:nvars), STAT=ierr)
    jerr = jerr + ierr
    ALLOCATE (mm(1:nvars), STAT=ierr)
    jerr = jerr + ierr
    ALLOCATE (r(1:nobs), STAT=ierr)
    jerr = jerr + ierr
    ALLOCATE (rmat(1:nobs, 1:nzeta, 1:ntau), STAT=ierr)
    jerr = jerr + ierr
    IF (jerr /= 0) RETURN

    t = 0.0D0
    oldtheta = 0.0D0
    m = 0
    mm = 0
    r = y
    rmat = 0.0D0

    npass = 0
    napenalty = 0
    maj = 2.0D0*maj

    DO l = 1, npenalty

        izeta = ind(l, 1)
        itau = ind(l, 2)

        ! Restore the complete state of the selected warm-start point.
        t = 0.0D0
        m = 0
        mm = 0

        IF (l == 1) THEN
            ni = 0
            r = y
        ELSE
            prevz = ind(l, 3)
            prevt = ind(l, 4)

            ni = ntheta(prevz, prevt)
            r = rmat(:, prevz, prevt)
            t(0) = theta0(1, prevz, prevt)

            IF (ni > 0) THEN
                m(1:ni) = itheta(1:ni, prevz, prevt)
                t(m(1:ni)) = theta(1:ni, prevz, prevt)

                DO j = 1, ni
                    mm(m(j)) = j
                END DO
            END IF
        END IF

        ! The outer/middle/inner loops below follow lslassoNETpath.
        DO
            IF (intr == 1) oldtheta(0) = t(0)
            IF (ni > 0) oldtheta(m(1:ni)) = t(m(1:ni))

            DO
                npass = npass + 1
                dif = 0.0D0

                DO k = 1, nvars
                    IF (ju(k) /= 0) THEN
                        oldt = t(k)
                        u = DOT_PRODUCT(r, x(:, k))
                        u = maj(k)*t(k) + u/DBLE(nobs)

                        pospenalty = uzeta(izeta)*pf(k)
                        negpenalty = utau(itau) + pospenalty

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
                                m(ni) = k
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

                DO
                    npass = npass + 1
                    dif = 0.0D0

                    DO j = 1, ni
                        k = m(j)
                        oldt = t(k)
                        u = DOT_PRODUCT(r, x(:, k))
                        u = maj(k)*t(k) + u/DBLE(nobs)

                        pospenalty = uzeta(izeta)*pf(k)
                        negpenalty = utau(itau) + pospenalty

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

            vrg = 1
            IF ((t(0) - oldtheta(0))**2 >= eps) vrg = 0

            DO j = 1, ni
                IF ((t(m(j)) - oldtheta(m(j)))**2 >= eps) THEN
                    vrg = 0
                    EXIT
                END IF
            END DO

            IF (vrg == 1) EXIT
        END DO

        IF (ni > pmax) THEN
            jerr = -10000 - l
            EXIT
        END IF

        ! Save this grid point before moving to another branch.
        IF (ni > 0) THEN
            itheta(1:ni, izeta, itau) = m(1:ni)
            theta(1:ni, izeta, itau) = t(m(1:ni))
        END IF

        ntheta(izeta, itau) = ni
        theta0(1, izeta, itau) = t(0)
        rmat(:, izeta, itau) = r

        apenalty(l, 1) = uzeta(izeta)
        apenalty(l, 2) = utau(itau)
        napenalty = l

        ! dfmax is intentionally not used as a global stopping rule here.
        ! Unlike a one-dimensional decreasing lambda path, this two-dimensional
        ! traversal is not totally ordered by penalty strength. Stopping the
        ! entire grid after one point exceeds dfmax can skip valid stronger
        ! points on another branch.

    END DO

    DEALLOCATE (t, oldtheta, m, mm, r, rmat)
    RETURN

END SUBROUTINE suniZTNETpath
