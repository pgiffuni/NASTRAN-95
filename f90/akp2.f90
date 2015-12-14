SUBROUTINE akp2
     
 COMPLEX :: ai
 
 COMMON/blk1/scrk,sps,sns,dstr,ai,pi,del,sigma,beta,res
 
 gam=SQRT(del**2-scrk**2)
 s1=sns*gam
 c1=(sigma-s1)/2.0
 c2=(sigma+s1)/2.0
 dgda=del/gam
 d1=sps/2.0
 d2=sns/2.0*dgda
 dc1da=d1-d2
 dc2da=d1+d2
 res=1.0/gam*dgda+sns*COS(s1)/SIN(s1)*dgda  &
     -COS(c1)/SIN(c1)*dc1da-COS(c2)/SIN(c2)*dc2da
 RETURN
END SUBROUTINE akp2
