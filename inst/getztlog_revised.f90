! Revised getztlog (drop-in replacement for the one in src/getlambda.f90).
! Changes relative to the original:
!   1. z and t are no longer initialized in their declarations. An
!      initializer implies SAVE, so the original carried z and t over from
!      the previous call and returned MAX(old, new) maxima on later calls.
!   2. The gradient is computed on centered but unscaled x. suniZTlogNET
!      divides pf and pf3 by xnorm when isd = 1, so its entry thresholds are
!      on the unscaled x; scaling x here put uzeta(1) off by xnorm(j).
!   3. The grids are built as top*factor**(l-1) instead of starting from
!      big and extrapolating the first value through LOG. Same values and
!      diffz/difft, but no LOG(0) = NaN when the top value is 0 (e.g.
!      z == t gives utau(1) = 0), and nzeta, ntau < 3 no longer read past
!      the end of ztemp/ttemp.
!   4. maj is no longer scaled while still uninitialized (it is unused).
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

    INTEGER :: ju(nvars)
    DOUBLE PRECISION :: zf
    DOUBLE PRECISION :: tf
    DOUBLE PRECISION :: z
    DOUBLE PRECISION :: t
    DOUBLE PRECISION :: u
    DOUBLE PRECISION :: xmean(nvars)
    DOUBLE PRECISION :: xnorm(nvars)
    DOUBLE PRECISION :: maj(nvars)

    zf = 1.0D0
    tf = 1.0D0
    IF (flminz < 1.0D0) THEN
        flminz = MAX(mfl, flminz)
        IF (nzeta > 1) zf = flminz**(1.0D0/(DBLE(nzeta) - 1.0D0))
    ELSE
        diffz = -1.0D0
    END IF
    IF (flmint < 1.0D0) THEN
        flmint = MAX(mfl, flmint)
        IF (ntau > 1) tf = flmint**(1.0D0/(DBLE(ntau) - 1.0D0))
    ELSE
        difft = -1.0D0
    END IF

    CALL chkvars(nobs, nvars, x, ju)

    ! Never scale x here (isd only matters through pf/xnorm in suniZTlogNET).
    CALL standard(nobs, nvars, x, ju, 0, intr, xmean, xnorm, maj)

    IF (flminz < 1.0D0 .OR. flmint < 1.0D0) THEN
        ! Largest positive and negative null-model gradients, y in {-1, 1}.
        z = 0.0D0
        t = 0.0D0
        DO j = 1, nvars
            IF (ju(j) /= 0) THEN
                IF (vp(j) > 0.0D0) THEN
                    u = DOT_PRODUCT(y*0.5D0, x(:, j))
                    IF (u > 0.0D0) THEN
                        z = MAX(z, u/vp(j))
                    ELSE IF (u < 0.0D0) THEN
                        t = MAX(t, ABS(u)/vp(j))
                    END IF
                END IF
            END IF
        END DO
        z = z/DBLE(nobs)
        t = t/DBLE(nobs)
    END IF

    IF (flminz < 1.0D0) THEN
        uzeta(1) = MAX(z, t)
        DO l = 2, nzeta
            uzeta(l) = uzeta(l - 1)*zf
        END DO
        diffz = -LOG(zf)
    END IF

    IF (flmint < 1.0D0) THEN
        utau(1) = ABS(z - t)
        DO l = 2, ntau
            utau(l) = utau(l - 1)*tf
        END DO
        difft = -LOG(tf)
    END IF
END SUBROUTINE getztlog
