SUBROUTINE csqrtx(xx,y)
    !*******
    !     ROUTINE TO FIND THE COMPLEX SQUARE ROOT OF X AND STORE IT IN Y
    !*******
 
 
    DOUBLE PRECISION, INTENT(IN)             :: xx(2)
    DOUBLE PRECISION, INTENT(OUT)            :: y(2)
    DOUBLE PRECISION :: x(2), r
    x(1) = xx(1)
    x(2) = xx(2)
    r = DSQRT(x(1)**2+x(2)**2)
    y(1) = DSQRT(DABS(x(1)+r)/2.)
    y(2) = DSQRT(DABS(-x(1)+r)/2.)
    IF(x(2) == 0.0D0) RETURN
    y(2) = DSIGN(y(2),x(2))

    RETURN
END SUBROUTINE csqrtx
