SUBROUTINE rand3 (f,s,q,n)
     
!     COMPUTES  MEAN RESPONSE  Q
 
 
 REAL, INTENT(IN)                         :: f(1)
 REAL, INTENT(IN)                         :: s(1)
 REAL, INTENT(OUT)                        :: q(2)
 INTEGER, INTENT(IN)                      :: n
 INTEGER :: check   ,NAME(2)
 
 COMMON /condas/ phi     ,twopi    ,radeg    ,degra    , s4pisq
 DATA    NAME  / 4HRAND  ,4H3      /
 
!     F IS ARRAY OF FREQUENCIES
!     S IS ARRAY OF POWER SPECTRAL DENSITY FUNCTIONS
!     Q IS MEAN RESPONSE
!     N IS NUMBER OF FREQUENCIES
 
 sum1 = 0.0
 nn   = n - 1
 sum  = 0.0
 DO  i = 1,nn
   df   = f(i+1) - f(i)
   sum  = sum + (s(i)+s(i+1))*df
   fi   = f(i   )*f(i  )
   fi1  = f(i+1 )*f(i+1)
   fii1 = 2.*f(i)*f(i+1)
   alp  = (3.*fi+fii1+fi1)/6.
   bta  = (fi+fii1+3.*fi1)/6.
   sum1 = sum1 + (alp*s(i)+bta*s(i+1))*df
 END DO
 sum  = SQRT(sum*0.5)
 sum1 = SQRT(sum1*.5)
 q(1) = sum
 q(2) = 0.0
 q1   = q(1)
 IF (q1 /= 0.0) q(2) = sum1/q1
 check= 123456789
 RETURN
 
!     AUTOCORRALATION FUNCTION
 
 ENTRY rand4 (f,s,tau,r,n)
!     =========================
 
!     COMPUTES  AUTOCORRALATION FUNCTION R  AT TIME TAU
!     WHERE F,S AS ABOVE. IF TAU = 0.0  R = Q*Q
 
 IF (check /= 123456789) CALL mesage (-37,0,NAME)
 IF (tau == 0.0) GO TO 30
 nn  = n - 1
 a   = 2.0*phi*tau
 b   = 1.0/a
 sum = 0.0
 DO  i = 1,nn
   sum = sum + b*(s(i+1)-s(i))/(f(i+1)-f(i))*(COS(a*f(i+1))  &
       - COS(a*f(i))) + s(i+1)*SIN(a*f(i+1))-s(i)*SIN(a*f(i))
 END DO
 sum = sum*b
 r   = sum
 GO TO 40
 30 r   = q1*q1
 40 RETURN
END SUBROUTINE rand3
