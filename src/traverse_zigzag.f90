SUBROUTINE traverse_zigzag(nrow, ncol, lamrow, lamcol, nlam, lamoutput, lamindex)
    IMPLICIT NONE
    ! -------- INPUT VARIABLES -------- !
    INTEGER, INTENT(IN) :: nrow, ncol, nlam
    DOUBLE PRECISION, INTENT(IN) :: lamrow(nrow), lamcol(ncol)
    ! -------- OUTPUT VARIABLES -------- !
    DOUBLE PRECISION, INTENT(OUT) :: lamoutput(2, nlam)
    DOUBLE PRECISION, INTENT(OUT) :: lamindex(2, nlam)
    ! -------- LOCAL DECLARATIONS -------- !
    INTEGER :: i, j

    !-------- CHECK INPUTS --------!
    IF (nlam /= nrow * ncol) THEN
        PRINT *, "Error: nlam does not match the product of nrow and ncol."
        RETURN
    END IF

    DO i = 1, nrow
        IF(MOD(i, 2) == 1) THEN
            ! Odd row: traverse left to right
            DO j = 1, ncol
                lamoutput(1, (i - 1) * ncol + j) = lamrow(i)
                lamoutput(2, (i - 1) * ncol + j) = lamcol(j)
                lamindex(1, (i - 1) * ncol + j) = i
                lamindex(2, (i - 1) * ncol + j) = j
            END DO
        ELSE
            ! Even row: traverse right to left
            DO j = ncol, 1, -1
                lamoutput(1, (i - 1) * ncol + (ncol - j + 1)) = lamrow(i)
                lamoutput(2, (i - 1) * ncol + (ncol - j + 1)) = lamcol(j)
                lamindex(1, (i - 1) * ncol + (ncol - j + 1)) = i
                lamindex(2, (i - 1) * ncol + (ncol - j + 1)) = j
            END DO
        END IF
    END DO


END SUBROUTINE traverse_zigzag


SUBROUTINE traverse_monoexpand(nrow, ncol, lamrow, lamcol, nlam, lamoutput, lamindex)
    IMPLICIT NONE
    ! -------- INPUT VARIABLES -------- !
    INTEGER, INTENT(IN) :: nrow, ncol, nlam
    DOUBLE PRECISION, INTENT(IN) :: lamrow(nrow), lamcol(ncol)
    ! -------- OUTPUT VARIABLES -------- !
    DOUBLE PRECISION, INTENT(OUT) :: lamoutput(2, nlam)
    DOUBLE PRECISION, INTENT(OUT) :: lamindex(2, nlam)
    ! -------- LOCAL DECLARATIONS -------- !
    INTEGER :: queue(2, nlam)
    INTEGER :: front, rear, count
    INTEGER :: r, c, new_r, new_c
    LOGICAL :: visited(nrow, ncol)

    !-------- CHECK INPUTS --------!
    IF (nlam /= nrow * ncol) THEN
        PRINT *, "Error: nlam does not match the product of nrow and ncol."
        RETURN
    END IF

    !-------- INITIALIZE --------!
    visited = .FALSE.
    queue = 0
    front = 1
    rear = 1
    count = 0

    !-------- START FROM (1, 1) --------!
    queue(1, rear) = 1
    queue(2, rear) = 1
    visited(1,1) = .TRUE.

    DO WHILE (front <= rear .AND. count < nlam)
        r = queue(1, front)
        c = queue(2, front)
        front = front + 1

        count = count + 1
        lamoutput(1, count) = lamrow(r)
        lamoutput(2, count) = lamcol(c)
        lamindex(1, count) = DBLE(r)
        lamindex(2, count) = DBLE(c)

        ! Try to add (r+1, c)
        new_r = r + 1
        new_c = c
        IF (new_r <= nrow .AND. .NOT. visited(new_r, new_c)) THEN
            rear = rear + 1
            queue(1, rear) = new_r
            queue(2, rear) = new_c
            visited(new_r, new_c) = .TRUE.
        END IF

        ! Try to add (r, c+1)
        new_r = r
        new_c = c + 1
        IF (new_c <= ncol .AND. .NOT. visited(new_r, new_c)) THEN
            rear = rear + 1
            queue(1, rear) = new_r
            queue(2, rear) = new_c
            visited(new_r, new_c) = .TRUE.
        END IF
    END DO

END SUBROUTINE traverse_monoexpand