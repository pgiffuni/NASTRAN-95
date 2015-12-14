SUBROUTINE akappa (arg,bkappa)
     
!     SUBROUTINE FOR COMPUTING KAPPA
 
 
 REAL, INTENT(IN)                         :: arg
 REAL, INTENT(OUT)                        :: bkappa
 COMPLEX :: ai
 
 COMMON/blk1/scrk,sps,sns,dstr,ai,pi,del,sigma,beta,res
 
 scrk1 =  ABS (scrk)
 arg1  =  ABS (arg)
 IF ( scrk1 > arg1)  GO TO 10
 gam=SQRT(arg**2-scrk**2)
 s1=sns*gam
 c1=beta*gam*SIN(s1)
 c2=COS(s1)-COS((arg-del)*sps+sigma)
 bkappa=c1/c2
 RETURN
 10 CONTINUE
 gam=SQRT(scrk**2-arg**2)
 s1=sns*gam
 c1=-beta*gam*SINH(s1)
 c2=COSH(s1)-COS((arg-del)*sps+sigma)
 bkappa=c1/c2
 RETURN
END SUBROUTINE akappa
