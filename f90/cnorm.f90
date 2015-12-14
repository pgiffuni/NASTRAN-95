SUBROUTINE cnorm(x,div,y)
     
    !     CNORM WILL NORMALIZE X TO THE MAXIMUM ELEMENT EQUAL TO A MODULUS
    !     OF ONE AND STORE THE DIVISOR IN MAX (X MAY BE COMPLEX)

    DOUBLE PRECISION, INTENT(IN OUT)         :: x(1)
    DOUBLE PRECISION, INTENT(OUT)            :: div(2)
    DOUBLE PRECISION, INTENT(IN OUT)         :: y(1)
    DOUBLE PRECISION :: MAX,temp, SIGN,cosang,xo,d,r,ri
    COMMON /system/  ibuf,nout
    COMMON /cinvpx/  filek(7)
    COMMON /cinvxx/  dum(30),ind1,iter
    EQUIVALENCE      (ncol,filek(2))
 
    ncol2 = ncol + ncol
    MAX   = 0.d0
    SIGN  = 1.0D0
    ind   = 0

    DO  i = 1,ncol2,2
        IF (x(i)**2+x(i+1)**2 <= MAX) CYCLE
        MAX = x(i)**2 + x(i+1)**2
        ind = i
    END DO

    IF (ind  ==   0) GO TO 80
    IF (iter ==   1) GO TO 60
    IF (ind == ind1) GO TO 40
    CALL sswtch (12,idiag)
    IF (idiag ==  0) GO TO 40
    WRITE  (6,30) ind,ind1
30  FORMAT (10H change     ,2I5)
40 CONTINUE

   d  = x(ind)**2 + x(ind+1)**2
   r  = (x(ind1)*x(ind) + x(ind1+1)*x(ind+1))/d
   ri = (x(ind1+1)*x(ind) - x(ind1)*x(ind+1))/d
   cosang = xo*r/DSQRT(r**2 + ri**2)
   IF (DABS(cosang+1.d0) <= 0.1D0) SIGN = -1.0D0
60 i  = ind
   div(1) = x(i  )*SIGN
   div(2) = x(i+1)*SIGN
   ind1 = ind
   MAX  = 1.0D0/MAX

   DO  i = 1,ncol2,2
       temp   = (x(i)*div(1)+x(i+1)*div(2))*MAX
       x(i+1) = (x(i+1)*div(1)-x(i)*div(2))*MAX
       x(i) = temp
   END DO

   xo = x(ind)
   RETURN
 
80 WRITE  (nout,90)
90 FORMAT (//5X,37HCONOR  received a vector of all zeros)
   CALL mesage (-37,0,0)

   RETURN
END SUBROUTINE cnorm
