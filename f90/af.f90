SUBROUTINE af (f,n,a,b,c,c1,c2,c3,t1,t3,t5,jump)
     
    !     THIS AREA INTEGRATION ROUTINE IS USED IN TRIM6, TRPLT1 AND TRSHL
    !     IT COMPUTES THE F FUNCTION, AND CONSTANTS C1, C2, C3
 
    !     FAC ARE THE FACTORIALS 1 THRU 36
    !     B   IS  DISTANCE OF GRID POINT 1
    !     A   IS  DISTANCE OF GRID POINT 3
    !     C   IS  DISTANCE OF GRID POINT 5
    !     T1  IS  ASSOCIATIVE VARIABLE AT GRID POINT 1
    !     T3  IS  ASSOCIATIVE VARIABLE AT GRID POINT 3
    !     T5  IS  ASSOCIATIVE VARIABLE AT GRID POINT 5
    !     N   IS  DIMENSION OF AREA FUNCTION F
 
    REAL, INTENT(OUT)                        :: f(n,n)
    INTEGER, INTENT(IN)                      :: n
    REAL, INTENT(IN)                         :: a
    REAL, INTENT(IN)                         :: b
    REAL, INTENT(IN)                         :: c
    REAL, INTENT(OUT)                        :: c1
    REAL, INTENT(OUT)                        :: c2
    REAL, INTENT(OUT)                        :: c3
    REAL, INTENT(IN OUT)                     :: t1
    REAL, INTENT(IN)                         :: t3
    REAL, INTENT(IN OUT)                     :: t5
    INTEGER, INTENT(IN OUT)                  :: jump
 
    DOUBLE PRECISION :: fac(20), temp
    DATA fac / 1.d0,1.d0, 2.d0,6.d0, 2.4D1,  1.2D2, 7.2D2, 5.04D3,  &
        4.032D4,       3.6288D5,     3.6288D6,     3.99168D7,  &
        4.790016D8,    6.227021D9,   8.7178291D10, 1.307674D12,  &
        2.092279D13,   3.556874D14,  6.402374D15,  1.216451D17/
 
    IF (jump > 0) GO TO 30
    IF (n > 18) STOP 'in AF'
    DO  i=1,n
        DO  j=1,n
            f(i,j)=0.0
        END DO
    END DO
    DO  i=1,n
        i1=i
        DO  j=1,i
            temp = DBLE(c**j) * fac(i1) / fac(i+2)
            temp = DBLE(a**i1-(-b)**i1) * temp * fac(j)
            f(i1,j) = SNGL(temp)
            i1=i1-1
        END DO
    END DO
    IF (jump < 0) RETURN
 
30  ab=a-b
    IF (a == b .AND. a /= 0.0) ab=a+b
    IF (ab == 0.0) CALL mesage (-37,0,0)
    c1=(t1*a-t3*b)/ab
    c2=(t3-t1)/ab
    c3=(t5-c1)/c

    RETURN
END SUBROUTINE af
