MODULE coxfunc
    IMPLICIT NONE
CONTAINS

    FUNCTION nonzero(n, v) RESULT(res)
        IMPLICIT NONE
        INTEGER :: n
        DOUBLE PRECISION :: v(n)

        INTEGER :: res, i
        res = 0

        DO i = 1, n
            IF (v(i) /= 0.0D0) THEN
                res = 1
                EXIT
            END IF
        END DO 
    END FUNCTION nonzero   

    FUNCTION risk(nobs, nvars, nk, d, dk, f, e, kp, jp, u) RESULT(res)
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: nobs, nvars, nk
        INTEGER, INTENT(IN) :: kp(nk), jp(nobs)
        DOUBLE PRECISION, INTENT(IN) :: d(nobs), dk(nk), f(nobs)
        DOUBLE PRECISION, INTENT(INOUT) :: u(nk)
        DOUBLE PRECISION, INTENT(IN) :: e(nobs)

        DOUBLE PRECISION :: res

        INTEGER :: nvars_local

        nvars_local = nvars
        CALL usk(nobs, nk, kp, jp, e, u)
        u = LOG(u)
        res = DOT_PRODUCT(d, f) - DOT_PRODUCT(dk, u)
    END FUNCTION risk

    SUBROUTINE usk(nobs, nk, kp, jp, e, u)
        IMPLICIT NONE

        INTEGER, INTENT(IN) :: nobs, nk
        INTEGER, INTENT(IN) :: kp(nk), jp(nobs)
        DOUBLE PRECISION, INTENT(IN) :: e(nobs)
        DOUBLE PRECISION, INTENT(OUT) :: u(nk)

        INTEGER :: k, j, j1, j2
        DOUBLE PRECISION :: h

        h = 0.0D0

        DO k = nk, 1, -1
            j2 = kp(k)
            j1 = 1
            IF (k > 1) j1 = kp(k-1) + 1

            DO j = j2, j1, -1
            h = h + e(jp(j))
            END DO

            u(k) = h
        END DO

        RETURN
    END SUBROUTINE usk

    SUBROUTINE vars(nobs, nvars, x, w, ixx, v)
        INTEGER :: nobs, nvars
        INTEGER :: ixx(nvars)
        DOUBLE PRECISION :: x(nobs, nvars), w(nobs)
        DOUBLE PRECISION :: v(nvars)

        INTEGER :: j

        DO j = 1, nvars
            IF (ixx(j) > 0) THEN
                v(j) = DOT_PRODUCT(w, x(:, j)**2)
            END IF
        END DO

    END SUBROUTINE vars

END MODULE coxfunc