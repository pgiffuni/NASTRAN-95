SUBROUTINE cxloop (x,y,n)
     
    DOUBLE PRECISION, INTENT(OUT)            :: x(1)
    DOUBLE PRECISION, INTENT(IN)             :: y(1)
    INTEGER, INTENT(IN)                      :: n
 
    DOUBLE PRECISION :: xx(2) , yy(2), mpy(2)
 
    nn = n + n
    DO  i = 1,nn
        x(i) = y(i)
    END DO
    RETURN
    ENTRY cloop( xx, yy, mpy, m)
    mm = m+m
    DO  i = 1,mm,2
        xx(i) = xx(i) - mpy(1)*yy(i) + mpy(2) * yy(i+1)
        xx(i+1) = xx(i+1) -mpy(2)*yy(i) -mpy(1)*yy(i+1)
    END DO

    RETURN
END SUBROUTINE cxloop
