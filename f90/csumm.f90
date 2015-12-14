SUBROUTINE  csumm(d1,d2,id1,d3,d4,id2,d5,d6,id5)
     
    !     ADDS  D1+D2 TO  D3+D4 SCALING OUTPUT
 
 
    DOUBLE PRECISION, INTENT(IN OUT)         :: d1
    DOUBLE PRECISION, INTENT(IN OUT)         :: d2
    INTEGER, INTENT(IN OUT)                  :: id1
    DOUBLE PRECISION, INTENT(IN OUT)         :: d3
    DOUBLE PRECISION, INTENT(IN OUT)         :: d4
    INTEGER, INTENT(IN OUT)                  :: id2
    DOUBLE PRECISION, INTENT(OUT)            :: d5
    DOUBLE PRECISION, INTENT(OUT)            :: d6
    INTEGER, INTENT(OUT)                     :: id5
    DOUBLE PRECISION :: t1,t2,t3,t4
    mult = IABS(id1-id2)
    IF(mult <= 38) factor = 10.0**mult
    t1 =d1
    t2 =d2
    t3 =d3
    t4 =d4
    id5 =id1
    IF(id1-id2 < 0) THEN
        GO TO    30
    ELSE IF (id1-id2 == 0) THEN
        GO TO    50
    ELSE
        GO TO    20
    END IF
30  IF(mult > 38) GO TO 40
    t3 =t3*factor
    t4 =t4*factor
    GO TO 50
20  IF(mult > 38) GO TO 35
    t1 = t1*factor
    t2 = t2*factor
    id5= id2
    GO TO 50
35  d5 = d3
    d6 =d4
    id5 = id2
    GO TO 70
40  d5 = d1
    d6 = d2
    GO TO 70
50  d5 = t1 +t3
    d6 = t2 + t4
70  RETURN
    ENTRY csqrtn(d1,d2,id1,d3,d4,id2)
 
    !     COMPUTES COMPLEX SQRT = SCALED
 
    id2 = id1
    d3=d1
    d4= d2
    IF( MOD(id1,2) == 0) GO TO 100
    id2 = id2 -1
    IF(id2 < 0) GO TO 105
101 d3 = d3*10.0
    d4 =d4*10.0
100 id2 = id2/2
    t1 =DSQRT(d3*d3 +d4*d4)
    t2 = DSQRT( DABS(d3+t1)/2.0)
    t3 = DSQRT(DABS(-d3+t1)/2.0)
    d3 =t2
    d4 = t3
    IF(d2 == 0.0D0) GO TO 70
    d4 =DSIGN(t3,d2)
    GO TO 70
 
    !     NEGATIVE EXPONENT
 
105 id2 = id2+1
    GO TO 101
 
    !     SCALES DETERMINANT
 
    ENTRY cdetm3(d1,d2,id1)
    t1 = DMAX1(DABS(d1),DABS(d2))
    IF(t1 == 0.0D0) GO TO 70
4125 IF(t1 > 10.0D0) GO TO 4153
4126 IF(t1 < 1.0D0) GO TO 4140
    GO TO 70
4153 d1 = d1*0.1D0
    d2 = d2*0.1D0
    t1 = t1*0.1D0
    id1 = id1+1
    GO TO 4125
4140 d1 = d1*10.0D0
    d2 = d2*10.0D0
    t1 = t1*10.0D0
    id1 = id1-1
    GO TO 4126

END SUBROUTINE  csumm
