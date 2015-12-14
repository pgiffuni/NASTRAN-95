SUBROUTINE cdetm2 (p,d,ip,pr,pi,dr,di,ips1)
     
    !     ARRANGES  P,D,IP  IN ORDER BY MAGNITUDE OF DETERMINANT
 
 
    DOUBLE PRECISION, INTENT(IN)             :: p(6)
    DOUBLE PRECISION, INTENT(IN)             :: d(6)
    INTEGER, INTENT(IN)                      :: ip(6)
    DOUBLE PRECISION, INTENT(OUT)            :: pr(3)
    DOUBLE PRECISION, INTENT(OUT)            :: pi(3)
    DOUBLE PRECISION, INTENT(OUT)            :: dr(3)
    DOUBLE PRECISION, INTENT(OUT)            :: di(3)
    INTEGER, INTENT(OUT)                     :: ips1(3)
    INTEGER :: ips(3)
    DOUBLE PRECISION :: d1,d2,d3,dd(3), d4
    EQUIVALENCE      (d1,dd(1)),(d2,dd(2)),(d3,dd(3))
 
    d1 = d(1)*d(1) + d(2)*d(2)
    d2 = d(3)*d(3) + d(4)*d(4)
    d3 = d(5)*d(5) + d(6)*d(6)
    DO  i = 1,3
        dd(i) = DSQRT(dd(i))
    END DO
    DO  i = 2,6,2
        k = i/2
        ips(k)  = ip(i)
        ips1(k) = ip(i)
    END DO
 
    !     SAVE STUFF IN OUTPUT AREAS
 
    DO  i = 1,3
        pr(i) = p(2*i-1)
        pi(i) = p(2*i  )
        dr(i) = d(2*i-1)
        di(i) = d(2*i  )
    END DO
 
    !     SCALE  MAGNITUDES
 
    DO  i = 1,3
40      IF (dd(i) > 10.0D0) GO TO 60
50      IF (dd(i) <  1.0D0) GO TO 70
        CYCLE
60      dd(i)  = dd(i)*0.1D0
        ips(i) = ips(i) + 1
        GO TO 40
70      dd(i)  = dd(i)*10.0D0
        ips(i) = ips(i) - 1
        GO TO 50
    END DO
 
    !     START COMPARISON TESTS
 
    IF (ips(1) > ips(2) .AND. ips(2) > ips(3)) GO TO 190
    IF (ips(1) > ips(2) .AND. ips(1) > ips(3)) GO TO 160
    IF (ips(2)-ips(3) < 0) THEN
        GO TO   100
    ELSE IF (ips(2)-ips(3) == 0) THEN
        GO TO    90
    ELSE
        GO TO   130
    END IF
90  IF (d2 >= d3) GO TO 130
100 IF (ips(1)-ips(3) < 0) THEN
        GO TO   120
    ELSE IF (ips(1)-ips(3) == 0) THEN
        GO TO   110
    ELSE
        GO TO   160
    END IF
110 IF (d1 >= d3) GO TO 160
120 is1 = 1
    is2 = 3
    ASSIGN 160 TO isret
    GO TO 200
130 IF (ips(1)-ips(2) < 0) THEN
        GO TO   150
    ELSE IF (ips(1)-ips(2) == 0) THEN
        GO TO   140
    ELSE
        GO TO   160
    END IF
140 IF (d1 >= d2) GO TO 160
150 is1 = 1
    is2 = 2
    ASSIGN 160 TO isret
    GO TO 200
160 IF (ips(2)-ips(3) < 0) THEN
        GO TO   180
    ELSE IF (ips(2)-ips(3) == 0) THEN
        GO TO   170
    ELSE
        GO TO   190
    END IF
170 IF (d2 >= d3) GO TO 190
180 is1 = 2
    is2 = 3
    ASSIGN 190 TO isret
    GO TO 200
190 RETURN
 
    !      SWITCHES VALUES ON IS1, IS2
 
200 nx = ips(is1)
    ips(is1) = ips(is2)
    ips(is2) = nx
    nx       = ips1(is1)
    ips1(is1)= ips1(is2)
    ips1(is2)= nx
    d4       = pr(is1)
    pr(is1)  = pr(is2)
    pr(is2)  = d4
    d4       = pi(is1)
    pi(is1)  = pi(is2)
    pi(is2)  = d4
    d4       = dr(is1)
    dr(is1)  = dr(is2)
    dr(is2)  = d4
    d4       = di(is1)
    di(is1)  = di(is2)
    di(is2)  = d4
    d4       = dd(is1)
    dd(is1)  = dd(is2)
    dd(is2)  = d4
    GO TO isret, (160,190)

END SUBROUTINE cdetm2
