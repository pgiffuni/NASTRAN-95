SUBROUTINE alg18 (beta1,beta2,i1,i2,fact,x0,y0,s0,xr,y1,x1,y2,rdi  &
        us,s,c1)
     
 
 
 REAL, INTENT(IN OUT)                     :: beta1
 REAL, INTENT(IN OUT)                     :: beta2
 INTEGER, INTENT(IN)                      :: i1
 INTEGER, INTENT(IN)                      :: i2
 REAL, INTENT(IN)                         :: fact
 REAL, INTENT(IN)                         :: x0
 REAL, INTENT(IN)                         :: y0
 REAL, INTENT(IN)                         :: s0
 REAL, INTENT(IN)                         :: xr
 REAL, INTENT(OUT)                        :: y1
 REAL, INTENT(OUT)                        :: x1
 REAL, INTENT(OUT)                        :: y2
 REAL, INTENT(IN OUT)                     :: us
 REAL, INTENT(OUT)                        :: s(80)
 REAL, INTENT(IN)                         :: c1
 
 
 delx=xr/FLOAT(i2-i1)
 xx=x0
 i3=i1+1
 IF (beta1 == beta2) GO TO 20
 y1=-COS(beta1/c1)/(SIN(beta1/c1)-SIN(beta2/c1))
 x1=SIN(beta1/c1)/(SIN(beta1/c1)-SIN(beta2/c1))
 y2=TAN((beta1+beta2)/(2.0*c1))
 rdius=ABS(1.0/(SIN(beta1/c1)-SIN(beta2/c1)))
 y2=y2*fact+y0
 y1=y1*fact+y0
 x1=x1*fact+x0
 rdius=rdius*fact
 DO  j=i3,i2
   xx=xx+delx
   phi1=ATAN(-1./SQRT(rdius**2-(xx-x1)**2)*(xx-x1))
   IF ((beta1-beta2) < 0.0) phi1=-phi1
   phi2=ABS(beta1/c1-phi1)
   s(j)=rdius*phi2+s0
 END DO
 RETURN
 20    am=TAN(beta1/c1)
 DO  j=i3,i2
   xx=xx+delx
   s(j)=(xx-x0)*SQRT(am*am+1.0)+s0
 END DO
 y2=am*(xx-x0)+y0
 RETURN
END SUBROUTINE alg18
