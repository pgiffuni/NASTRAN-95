SUBROUTINE fvrs1b (base,w1,nf)
     
!     SUBROUTINE TO COMPUTE BASE(FI)(3X1) FOR MODFRL=FALSE
 
 
 COMPLEX, INTENT(OUT)                     :: base(3,nf)
 REAL, INTENT(IN)                         :: w1(nf)
 INTEGER, INTENT(IN)                      :: nf
 COMPLEX :: z1,p
 
 
 
 COMMON /condas/ pi,twopi,radeg,degra,s4piq
 COMMON /BLANK / dum(5),it(6),dum1(3)
 
 DO  k=1,nf
   f=w1(k)/twopi
   LT=1
   lp=2
   DO  i=1,3
     IF(it(LT) == -1)GO TO 90
     CALL tab(it(LT),f,xo)
     IF(it(lp) == -1)GO TO 40
     CALL tab(it(lp),f,phi)
     rad=phi*degra
     z1=CMPLX(0.0,rad)
     p=CEXP(z1)
     GO TO 50
     40  p=(1.0,0.0)
     50  base(i,k)=xo*p
     GO TO 95
     90  base(i,k)=(0.0,0.0)
     95  LT=LT+2
     lp=lp+2
   END DO
 END DO
 RETURN
END SUBROUTINE fvrs1b
