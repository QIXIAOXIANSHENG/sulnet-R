SUBROUTINE suniZTlogNET(lam2, nobs, nvars, x, y, jd, pf, pf2, pf3, pmax, nzeta,&
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
    IF (MAXVAL(pf3) <= 0.0D0) THEN
        jerr = 10000
        RETURN
    END IF
    pf = MAX(0.0D0, pf)
    pf2 = MAX(0.0D0, pf2)
    pf3 = MAX(0.0D0, pf3)
    CALL standard(nobs, nvars, x, ju, isd, intr, xmean, xnorm, maj)
    CALL gridIndexZT(nzeta, ntau, uzeta, utau, ind, diffzeta, difftau)
    CALL suniZTlogNETpath(lam2, nobs, nvars, x, y, ju, pf, pf2, pf3, &
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
END SUBROUTINE suniZTlogNET


! Positive threshold zeta*pf and negative threshold zeta*pf + tau*pf3,
! solved cell by cell by the shared IRLS grid solver in sunilogpath.f90.
SUBROUTINE suniZTlogNETpath(lam2, nobs, nvars, x, y, ju, pf, pf2, pf3, &
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

    DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: pen

    INTEGER :: l
    INTEGER :: iz
    INTEGER :: it

    ALLOCATE (pen(1:(nzeta*ntau), 1:3), STAT=jerr)
    IF (jerr /= 0) RETURN

    DO l = 1, nzeta*ntau
        iz = ind(l, 1)
        it = ind(l, 2)
        pen(l, 1) = uzeta(iz)
        pen(l, 2) = utau(it)
        pen(l, 3) = uzeta(iz)
        apenalty(l, 1) = uzeta(iz)
        apenalty(l, 2) = utau(it)
    END DO
    CALL sunilogpath(lam2, nobs, nvars, x, y, ju, pf, pf2, pf3, pmax, &
        & nzeta, ntau, pen, ind, eps, maxit, intr, napenalty, theta0, theta, &
        & itheta, ntheta, npass, jerr)
    DEALLOCATE (pen)
    RETURN
END SUBROUTINE suniZTlogNETpath
