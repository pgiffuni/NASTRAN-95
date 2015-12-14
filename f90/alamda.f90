SUBROUTINE alamda (arg,y,blamda)
     
!     SUBROUTINE FOR COMPUTING LAMDA
 
 
 REAL, INTENT(IN)                         :: arg
 REAL, INTENT(IN OUT)                     :: y
 COMPLEX, INTENT(OUT)                     :: blamda
 COMPLEX :: ai,c1
 
 COMMON/blk1/scrk,sps,sns,dstr,ai,pi,del,sigma,beta,res
 
 scrk1= ABS(scrk)
 arg1=  ABS(arg)
 s1=(arg-del)*sps+sigma
 IF( scrk1 > arg1) GO TO 10
 gam=SQRT(arg**2-scrk**2)
 c1=COS(gam*(sns-y))-CEXP(ai*s1)*COS(gam*y)
 c2=COS(sns*gam)-COS(s1)
 blamda=c1/c2
 RETURN
 10 CONTINUE
 gam=SQRT(scrk**2-arg**2)
 c1=COSH(gam*(sns-y))-CEXP(ai*s1)*COSH(gam*y)
 c2=COSH(sns*gam)-COS(s1)
 blamda=c1/c2
 RETURN
END SUBROUTINE alamda
