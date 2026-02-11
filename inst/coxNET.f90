SUBROUTINE coxNET(parm, nobs, nvars, x, y, d, g, w, jd, vp, cl, &
    & ne, nx, nlam, flmin, ulam, thr, maxit, isd, lmu, ca, ia, &
    & nin, dev0, dev, alm, nlp, jerr)
    IMPLICIT NONE
    ! ---------Input Variables-------------
    INTEGER          :: nobs, nvars, ne, nx, nlam
    INTEGER          :: maxit, isd, lmu, nlp, jerr

    DOUBLE PRECISION :: parm, flmin, thr, dev0
    DOUBLE PRECISION :: x(nobs,nvars)
    DOUBLE PRECISION :: y(nobs), d(nobs), g(nobs), w(nobs)
    DOUBLE PRECISION :: vp(nvars)
    DOUBLE PRECISION :: ulam(nlam)
    DOUBLE PRECISION :: ca(nx,nlam)
    DOUBLE PRECISION :: dev(nlam), alm(nlam)
    DOUBLE PRECISION :: cl(2,nvars)

    INTEGER          :: jd(*)
    INTEGER          :: ia(nx), nin(nlam)

    DOUBLE PRECISION, ALLOCATABLE :: xs(:), ww(:), vq(:)
    INTEGER,          ALLOCATABLE :: ju(:)
    
    INTEGER :: j, k, nk
    DOUBLE PRECISION :: sw


    IF (MAXVAL(vp) <= 0.0D0) THEN
        jerr = 10000
        RETURN
    END IF
    ALLOCATE (ww(1:nobs), stat=jerr)
    IF (jerr /= 0) THEN
        RETURN 
    END IF
    ALLOCATE (ju(1:nvars), stat=jerr)
    IF (jerr /= 0) THEN
        RETURN 
    END IF
    ALLOCATE (vq(1:nvars), stat=jerr)
    IF (jerr /= 0) THEN
        RETURN 
    END IF
    IF (isd > 0) THEN
        ALLOCATE (xs(1:nvars), stat=jerr)
        IF (jerr /= 0) THEN
            RETURN 
        END IF 
    END IF
    CALL chkvars(nobs, nvars, x, ju)
    IF (jd(1) > 0) THEN
        ju(jd(2:(jd(1) + 1))) = 0
    END IF
    IF (maxval(ju) <= 0) THEN
        jerr = 7777
        RETURN
    END IF
    vq = MAX(0.0D0, vp)
    vq = vq * nvars / SUM(vq)
    ww = MAX(0.0D0, w)
    sw = SUM(ww)
    IF (sw <= 0.0) THEN
        jerr = 9999
        RETURN
    END IF
    ww = ww / sw
    CALL cstandard(nobs, nvars, x, ww, ju, isd, xs)
    IF (isd > 0) THEN
        DO j = 1, nvars
            cl(:, j) = cl(:, j) * xs(j)
        END DO
    END IF
    CALL coxnet1(parm, nobs, nvars, x, y, d, g, ww, ju, vq, cl, ne, nx, nlam, flmin, ulam, &
        & thr, isd, maxit, lmu, ca, ia, nin, dev0, dev, alm, nlp, jerr)
    IF (jerr > 0) THEN
        RETURN
    END IF
    dev0 = 2.0 * sw * dev0
    IF (isd > 0) THEN
        DO k = 1, lmu
            nk = nin(k)
            ca(1:nk, k) = ca(1:nk, k)/xs(ia(1:nk))
        END DO
    END IF
    DEALLOCATE (ww, ju, vq)
    IF (isd > 0) THEN
        DEALLOCATE (xs)
    END IF
    RETURN
END SUBROUTINE coxNET

SUBROUTINE cstandard(nobs, nvars, x, w, ju, isd, xs)
    IMPLICIT NONE
    INTEGER :: nobs, nvars, isd
    INTEGER, DIMENSION(nvars) :: ju
    DOUBLE PRECISION, DIMENSION(nobs, nvars) :: x
    DOUBLE PRECISION, DIMENSION(nobs) :: w
    DOUBLE PRECISION, DIMENSION(nvars) :: xs
    
    INTEGER :: j
    DOUBLE PRECISION :: xm

    DO j = 1, nvars
        IF (ju(j) /= 0) THEN
            xm = DOT_PRODUCT(w, x(:, j))
            x(:, j) = x(:, j) - xm
            IF (isd > 0) THEN
                xs(j) = SQRT(DOT_PRODUCT(w, x(:, j)**2))
                x(:, j) = x(:, j) / xs(j)
            END IF
        END IF
    END DO
    RETURN
END SUBROUTINE cstandard

SUBROUTINE coxnet1(parm, nobs, nvars, x, y, d, g, q, ju, vp, cl, ne, nx, nlam, flmin, &
                  & ulam, cthri, isd, maxit, lmu, ao, m, kin, dev0, dev, alm, nlp, jerr)
    USE coxfunc
    IMPLICIT NONE

    INTEGER :: nobs, nvars, ne, nx, nlam, isd, maxit, lmu, nlp, jerr
    DOUBLE PRECISION :: parm, flmin, cthri, dev0

    DOUBLE PRECISION :: x(nobs, nvars), y(nobs), q(nobs), d(nobs), g(nobs)
    DOUBLE PRECISION :: vp(nvars), ulam(nlam)
    DOUBLE PRECISION :: ao(nx, nlam), dev(nlam), alm(nlam), cl(2, nvars)

    INTEGER :: ju(nvars), m(nx), kin(nlam)

    DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: w, dk, v, xs, wr
    DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: a, as, f, dq
    DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: e, uu, ga
    INTEGER, DIMENSION(:), ALLOCATABLE :: jp, kp, mm, ixx

    DOUBLE PRECISION :: sml, eps, big, devmax, pmin, exmx, alpha, oma, &
        & al, t0, fmax, r0, rr, alf, cthr, al0, sa, omal, tlam, &
        & eqs, u, at, del, dli
    INTEGER :: mnlam, itrace, nlm, nk, i, mnl, l, nin, j, ilm, k, &
        & ix, me

    isd = isd * 1
    sml = 1.0D-2
    eps = 1.0D-6
    big = 9.9D35
    mnlam = 5
    devmax = 0.99
    pmin = 1.0D-9
    exmx = 250.0D0
    itrace = 0

    ALLOCATE (e(1:nobs), STAT=jerr)
    IF (jerr /= 0) RETURN

    ALLOCATE (uu(1:nobs), STAT=jerr)
    IF (jerr /= 0) RETURN

    ALLOCATE (f(1:nobs), STAT=jerr)
    IF (jerr /= 0) RETURN

    ALLOCATE (w(1:nobs), STAT=jerr)
    IF (jerr /= 0) RETURN

    ALLOCATE (v(1:nvars), STAT=jerr)
    IF (jerr /= 0) RETURN

    ALLOCATE (a(1:nvars), STAT=jerr)
    IF (jerr /= 0) RETURN

    ALLOCATE (as(1:nvars), STAT=jerr)
    IF (jerr /= 0) RETURN

    ALLOCATE (xs(1:nvars), STAT=jerr)
    IF (jerr /= 0) RETURN

    ALLOCATE (ga(1:nvars), STAT=jerr)
    IF (jerr /= 0) RETURN

    ALLOCATE (ixx(1:nvars), STAT=jerr)
    IF (jerr /= 0) RETURN

    ALLOCATE (jp(1:nobs), STAT=jerr)
    IF (jerr /= 0) RETURN

    ALLOCATE (kp(1:nobs), STAT=jerr)
    IF (jerr /= 0) RETURN

    ALLOCATE (dk(1:nobs), STAT=jerr)
    IF (jerr /= 0) RETURN

    ALLOCATE (wr(1:nobs), STAT=jerr)
    IF (jerr /= 0) RETURN

    ALLOCATE (dq(1:nobs), STAT=jerr)
    IF (jerr /= 0) RETURN

    ALLOCATE (mm(1:nvars), STAT=jerr)
    IF (jerr /= 0) RETURN

    CALL groups(nobs, y, d, q, nk, kp, jp, t0, jerr)

    IF (jerr /= 0) THEN
        DEALLOCATE (e, uu, w, dk, v, xs, f, wr, a, as, jp, kp, dq, mm, ga, ixx)
        RETURN
    END IF
    alpha = parm
    oma = 1.0 - alpha
    nlm = 0
    ixx = 0
    al = 0.0
    dq = d*q
    CALL died(nobs, nk, dq, kp, jp, dk)
    a = 0.0D0
    f(1) = 0.0D0
    fmax = LOG(HUGE(f(1)) * 0.1D0)
    IF (nonzero(nobs, g) /= 0) THEN
        f = g - DOT_PRODUCT(q, g)
        e = q * EXP(SIGN(MIN(ABS(f), fmax), f))
    ELSE
        f = 0.0
        e = q
    END IF
    r0 = risk(nobs, nvars, nk, dq, dk, f, e, kp, jp, uu)
    rr = -(DOT_PRODUCT(dk(1:nk), LOG(dk(1:nk))) + r0)
    dev0 = rr

    DO i = 1, nobs
        IF ((y(i) >= t0) .and. (q(i) >= 0.0)) EXIT
        w(i) = 0.0
        wr(i) = w(i)
    END DO

    CALL outer(nobs, nk, dq, dk, kp, jp, e, wr, w, jerr, uu)
    IF (jerr /= 0) THEN
        DEALLOCATE (e, uu, w, dk, v, xs, f, wr, a, as, jp, kp, dq, mm, ga, ixx)
        RETURN
    END IF
    alf = 1.0D0
    IF (flmin < 1.0) THEN
        eqs = MAX(eps, flmin)
        alf = eqs ** (1.0 / (nlam - 1))
    END IF
    m = 0
    mm = 0
    nlp = 0
    nin = nlp
    mnl = min(mnlam, nlam)
    as = 0.0
    cthr = cthri*dev0
    DO j = 1, nvars
        IF (ju(j) /= 0) THEN
            ga(j) = ABS(DOT_PRODUCT(wr, x(:, j)))
        END IF
    END DO

    DO ilm = 1, nlam
        al0 = al
        IF (flmin < 1.0) THEN
            IF (ilm == 1) THEN
                al = big
            ELSE IF (ilm == 2) THEN
                al0 = 0.0D0
                DO j = 1, nvars
                    IF (ju(j) /= 0) THEN
                        IF (vp(j) > 0.0) THEN
                            al0 = MAX(al0, ga(j) / vp(j))
                        END IF
                    END IF
                END DO
                al0 = al0 / MAX(parm, 1.0D-3)
                al = alf * al0
            ELSE
                al = al * alf
            END IF
        ELSE
            al = ulam(ilm)
        END IF

        sa = alpha * al
        omal = oma * al
        tlam = alpha * (2.0 * al - al0)

        DO k = 1, nvars
            IF (ixx(k) == 1) EXIT
            IF (ju(k) == 0) EXIT
            IF (ga(k) > tlam * vp(k)) THEN
                ixx(k) = 1
            END IF
        END DO
        DO !10360
            DO !10371
                IF (nlp > maxit) THEN
                    jerr = -ilm
                    RETURN
                END IF
                IF (nin > 0) THEN
                    as(m(1:nin)) = a(m(1:nin))
                END IF
                CALL vars(nobs, nvars, x, w, ixx, v)
                DO !10401
                    nlp = nlp + 1
                    dli = 0.0
                    DO j = 1, nvars
                        IF (ixx(j) == 0) CYCLE
                        u = a(j)*v(j) + DOT_PRODUCT(wr, x(:, j))
                        IF (ABS(u) > vp(j) * sa) THEN
                            at = MAX(cl(1, j), min(cl(2, j), SIGN(ABS(u) - vp(j)*sa, u)/(v(j) + vp(j)*omal)))
                        ELSE
                            at = 0.0D0
                        END IF
                        IF (at == a(j)) CYCLE
                        del = at - a(j)
                        a(j) = at
                        dli = max(dli, v(j)*del**2)
                        wr = wr - del*w*x(:, j)
                        f = f + del*x(:, j)
                        IF (mm(j) /= 0) CYCLE
                        nin = nin + 1
                        IF (nin > nx) EXIT
                        mm(j) = nin
                        m(nin) = j
                    END DO
                    IF (nin > nx) EXIT
                    IF (dli .lt. cthr) EXIT
                    IF (nlp > maxit) THEN
                        jerr = -ilm
                        RETURN
                    END IF
                    DO !10511
                        nlp = nlp + 1
                        dli = 0.0
                        DO l = 1, nin !10521
                            j = m(l)
                            u = a(j)*v(j) + DOT_PRODUCT(wr, x(:, j))
                            IF (ABS(u) > vp(j) * sa) THEN
                                at = MAX(cl(1, j), MIN(cl(2, j), SIGN(ABS(u) - vp(j)*sa, u)/(v(j) + vp(j)*omal)))
                            ELSE
                                at = 0.0D0
                            END IF
                            IF (at == a(j)) CYCLE
                            del = at - a(j)
                            a(j) = at
                            dli = MAX(dli, v(j)*del**2)
                            wr = wr - del*w*x(:, j)
                            f = f + del*x(:, j)
                        END DO
                        IF (dli < cthr) EXIT
                        IF (nlp <= maxit) CYCLE
                        jerr = -ilm
                        RETURN
                    END DO
                END DO
                IF (nin <= nx) THEN !10372
                    e = q*EXP(SIGN(MIN(ABS(f), fmax), f))
                    CALL outer(nobs, nk, dq, dk, kp, jp, e, wr, w, jerr, uu)
                    IF (jerr /= 0) THEN
                        jerr = jerr - ilm
                        DEALLOCATE (e, uu, w, dk, v, xs, f, wr, a, as, jp, kp, dq, mm, ga, ixx)
                        RETURN
                    END IF
                    ix = 0
                    DO j = 1, nin
                        k = m(j)
                        IF (v(k)*(a(k) - as(k))**2 < cthr) CYCLE
                        ix = 1
                        EXIT
                    END DO
                    IF (ix /= 0) EXIT
                    DO k = 1, nvars
                        IF (ixx(k) == 1) CYCLE
                        IF (ju(k) == 0) CYCLE
                        ga(k) = ABS(DOT_PRODUCT(wr, x(:, k)))
                        IF (ga(k) <= sa*vp(k)) CYCLE
                        ixx(k) = 1
                        ix = 1
                    END DO
                END IF
            END DO
            IF (ix /= 1) EXIT
        END DO
        IF (nin > nx) THEN
            jerr = -10000 - ilm
            EXIT
        END IF
        IF (nin > 0) THEN
            ao(1:nin, ilm) = a(m(1:nin))
        END IF
        kin(ilm) = nin
        alm(ilm) = al
        lmu = ilm
        dev(ilm) = (risk(nobs, nvars, nk, dq, dk, f, e, kp, jp, uu) - r0)/rr
        IF (ilm < mnl) EXIT
        IF (flmin >= 1.0) EXIT
        me = 0
        DO j = 1, nin
            IF (ao(j, ilm) /= 0.0) THEN
                me = me + 1
            END IF
        END DO
        IF (me > ne) EXIT
        IF ((dev(ilm) - dev(ilm - mnl + 1))/dev(ilm) < sml) EXIT
        IF (dev(ilm) > devmax) EXIT
    END DO
    g = f
    DEALLOCATE (e, uu, w, dk, v, xs, f, wr, a, as, jp, kp, dq, mm, ga, ixx)
    RETURN
END SUBROUTINE coxnet1




SUBROUTINE groups(nobs, y, d, q, nk, kp, jp, t0, jerr)
    IMPLICIT NONE
    INTEGER :: nobs, nk
    INTEGER :: jerr

    DOUBLE PRECISION :: y(nobs), d(nobs), q(nobs)
    INTEGER :: jp(nobs)
    INTEGER :: kp(*)

    DOUBLE PRECISION :: t0


    INTEGER :: j, nj, j0
    DOUBLE PRECISION :: yk

    DO j = 1, nobs
        jp(j) = j
    END DO

    CALL psort7(y, jp, 1, nobs)
    nj = 0
    DO j = 1, nobs
        IF (q(jp(j)) > 0) THEN 
            nj = nj + 1
            jp(nj) = jp(j)
        END IF
    END DO
    IF (nj == 0) THEN
        jerr = 20000
        RETURN
    END IF

    j = 1
    DO 
        IF (j > nj) EXIT
        IF (d(jp(j)) > 0.0) EXIT
        j = j + 1
    END DO

    IF (j >= nj - 1) THEN
        jerr = 30000
        RETURN
    END IF

    t0 = y(jp(j))
    j0 = j - 1
    IF (j0 > 0) THEN
        DO
            IF (y(jp(j0)) < t0) EXIT
            j0 = j0 - 1
            IF (j0 == 0) EXIT
        END DO
        IF (j0 > 0) THEN
            nj = nj - j0
            DO j = 1, nj
                jp(j) = jp(j + j0)
            END DO
        END IF
        
    END IF
    jerr = 0
    nk = 0
    yk = t0
    j = 2

    DO
        DO
            IF (d(jp(j)) > 0.0) EXIT 
            IF (y(jp(j)) > yk) EXIT
            j = j + 1
            IF (j > nj) EXIT
        END DO
        nk = nk + 1
        kp(nk) = j - 1
        IF (j > nj) EXIT
        IF (j == nj) THEN
            nk = nk + 1
            kp(nk) = nj
            EXIT
        END IF
        yk = y(jp(j))
        j = j + 1
    END DO
RETURN
END SUBROUTINE groups

SUBROUTINE psort7(v, a, ii, jj)
    IMPLICIT NONE
    INTEGER :: ii, jj
    INTEGER :: a(jj)
    DOUBLE PRECISION :: v(jj)
    INTEGER :: iu(20), il(20)
    INTEGER :: t, tt
    INTEGER :: m, i, j, k
    INTEGER :: cycle90goto, ij, l
    DOUBLE PRECISION :: vt, vtt
    m = 1
    i = ii
    j = jj
    cycle90goto = 0
    DO

        IF (i < j .or. cycle90goto == 20) THEN
            DO ! 20 - 90 LOOP
                k = i
                ij = (j + i)/2
                t = a(ij)
                vt = v(t)
                IF (v(a(i)) > vt) THEN
                    a(ij) = a(i)
                    a(i) = t
                    t = a(ij)
                    vt = v(t)
                END IF
                l = j
                IF (v(a(j)) < vt) THEN
                    a(ij) = a(j)
                    a(j) = t
                    t = a(ij)
                    vt = v(t)
                END IF
                IF (v(a(i)) > vt) THEN
                    a(ij) = a(i)
                    a(i) = t
                    t = a(ij)
                    vt = v(t)
                END IF
                
                DO 
                    l = l - 1
                    if (v(a(l)) <= vt) EXIT
                END DO
                tt = a(l)
                vtt = v(tt)

                DO
                    k = k + 1
                    if (v(a(k)) >= vt) EXIT
                END DO

                DO
                    IF (k > l) EXIT
                    a(l) = a(k)
                    a(k) = tt
                    DO 
                        l = l - 1
                        if (v(a(l)) <= vt) EXIT
                    END DO
                    tt = a(l)
                    vtt = v(tt)

                    DO
                        k = k + 1
                        if (v(a(k)) >= vt) EXIT
                    END DO
                END DO
                
                IF (l - i > j - k) THEN
                    il(m) = i
                    iu(m) = l
                    i = k
                    m = m + 1
                ELSE
                    il(m) = k
                    iu(m) = j
                    j = l
                    m = m + 1
                END IF
                
                IF (j - i <= 10) EXIT
            END DO
            
            IF (i == ii) CYCLE
            
            i = i - 1
            DO
                i = i + 1
                IF (i /= j) THEN
                    t = a(i + 1)
                    vt = v(t)
                ELSE
                    EXIT ! GOTO 80
                END IF
                IF (v(a(i)) <= vt) CYCLE
                k = i
                DO
                    a(k + 1) = a(k)
                    k = k - 1
                    IF (vt .lt. v(a(k))) CYCLE
                    a(k + 1) = t
                END DO

            END DO
        END IF
        
        IF (i >= j .or. cycle90goto == 80) THEN ! 80
            m = m - 1
            IF (m == 0) RETURN
            i = il(m)
            j = iu(m)
            IF (j - i > 10) THEN
                cycle90goto = 20
            ELSE IF (i == ii) THEN
                cycle90goto = 10
            ELSE
                cycle90goto = 0
            END IF
        END IF

        IF (cycle90goto /= 0) CYCLE
        
        i = i - 1
        DO
            i = i + 1
            IF (i /= j) THEN
                t = a(i + 1)
                vt = v(t)
            ELSE
                cycle90goto = 80
                EXIT ! GOTO 80
            END IF
            IF (v(a(i)) <= vt) CYCLE
            k = i
            DO
                a(k + 1) = a(k)
                k = k - 1
                IF (vt .lt. v(a(k))) CYCLE
                a(k + 1) = t
            END DO

        END DO
    END DO
END SUBROUTINE psort7

SUBROUTINE died(nobs, nk, d, kp, jp, dk)
    IMPLICIT NONE
    INTEGER :: nobs, nk
    INTEGER :: kp(nk), jp(nobs)
    DOUBLE PRECISION :: d(nobs), dk(nk)

    INTEGER :: k

    dk(1) = SUM(d(jp(1:kp(1))))
    DO k = 2, nk
        dk(k) = SUM(d(jp((kp(k - 1) + 1):kp(k))))
    END DO
    RETURN
END SUBROUTINE died





SUBROUTINE outer(nobs, nk, d, dk, kp, jp, e, wr, w, jerr, u)
    USE coxfunc
    IMPLICIT NONE

    INTEGER :: nobs, nk, jerr
    INTEGER, DIMENSION(nk)   :: kp
    INTEGER, DIMENSION(nobs) :: jp

    DOUBLE PRECISION, DIMENSION(nobs) :: d, e, wr, w, u
    DOUBLE PRECISION, DIMENSION(nk)   :: dk
    DOUBLE PRECISION :: b, c

    INTEGER :: j, k, i, j1, j2

    call usk(nobs, nk, kp, jp, e, u)
    b = dk(1)/u(1)
    c = dk(1)/u(1)**2
    jerr = 0

    DO j = 1, kp(1)
        i = jp(j)
        w(i) = e(i)*(b - e(i)*c)
        IF (w(i) > 0.0) THEN
            wr(i) = d(i) - e(i)*b
        ELSE
            jerr = -30000
            RETURN
        END IF
    END DO

    DO k = 2, nk
        j1 = kp(k - 1) + 1
        j2 = kp(k)
        b = b + dk(k) / u(k)
        c = c + dk(k) / u(k) ** 2
        DO j = j1, j2
            i = jp(j)
            w(i) = e(i)*(b - e(i)*c)
            IF (w(i) > 0.0) THEN
                wr(i) = d(i) - e(i)*b
            ELSE
                jerr = -30000
                RETURN
            END IF
        END DO
    END DO
    RETURN
END SUBROUTINE outer

