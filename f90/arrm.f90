SUBROUTINE arrm(p,d,nd)
     
    !     SCALED ARITHMETIC ROUTINES--ARRANGING ROUTINE
 
 
    DOUBLE PRECISION, INTENT(IN OUT)         :: p(3)
    DOUBLE PRECISION, INTENT(OUT)            :: d(3)
    INTEGER, INTENT(OUT)                     :: nd(3)
    DOUBLE PRECISION :: dx, px
 
 
    DO  i=1,3
        IF(d(i) == 0.0D0) CYCLE
10      IF(DABS(d(i)) >= 1.0) GO TO 20
        d(i) = d(i)*10.0
        nd(i) = nd(i)-1
        GO TO 10
20      IF(DABS(d(i)) < 10.0) CYCLE
        d(i) = d(i)*0.1
        nd(i) = nd(i)+1
        GO TO 20
    END DO
    IF(nd(1) > nd(2) .AND. nd(2) > nd(3)) RETURN
    IF(nd(1) > nd(2) .AND. nd(1) > nd(3)) GO TO 110
    IF(nd(2) - nd(3) < 0) THEN
        GO TO    50
    ELSE IF (nd(2) - nd(3) == 0) THEN
        GO TO    40
    ELSE
        GO TO    80
    END IF
40  IF(DABS(d(2)) >= DABS(d(3))) GO TO 80
50  IF(nd(1) - nd(3) < 0) THEN
        GO TO    70
    ELSE IF (nd(1) - nd(3) == 0) THEN
        GO TO    60
    ELSE
        GO TO   110
    END IF
60  IF(DABS(d(1)) >= DABS(d(3))) GO TO 110
70  nx=nd(1)
    dx=d(1)
    px=p(1)
    nd(1)=nd(3)
    d(1)=d(3)
    p(1)=p(3)
    nd(3)=nx
    d(3)=dx
    p(3)=px
    GO TO 110
80  IF(nd(1) - nd(2) < 0) THEN
        GO TO   100
    ELSE IF (nd(1) - nd(2) == 0) THEN
        GO TO    90
    ELSE
        GO TO   110
    END IF
90  IF(DABS(d(1)) >= DABS(d(2))) GO TO 110
100 nx=nd(1)
    dx=d(1)
    px=p(1)
    nd(1)=nd(2)
    d(1)=d(2)
    p(1)=p(2)
    nd(2)=nx
    d(2)=dx
    p(2)=px
110 IF(nd(2) - nd(3) < 0) THEN
        GO TO   130
    ELSE IF (nd(2) - nd(3) == 0) THEN
        GO TO   120
    ELSE
        GO TO   140
    END IF
120 IF(DABS(d(2)) >= DABS(d(3))) RETURN
130 nx=nd(2)
    dx=d(2)
    px=p(2)
    nd(2)=nd(3)
    d(2)=d(3)
    p(2)=p(3)
    nd(3)=nx
    d(3)=dx
    p(3)=px

140 RETURN
END SUBROUTINE arrm
