SUBROUTINE akapm (arg,bkpm)
     
!     SUBROUTINE FOR COMPUTING KAPPA MINUS
 
 
 COMPLEX, INTENT(IN)                      :: arg
 COMPLEX, INTENT(OUT)                     :: bkpm
 COMPLEX :: c1,ai,c1test,bsycon, at2,at3,alp0,alp,aln
 CHARACTER (LEN=23) :: ufm
 COMMON /xmssg / ufm
 COMMON /system/ sysbuf,ibbout
 COMMON /blk1  / scrk,sps,sns,dstr,ai,pi,del,sigma,beta,res
 COMMON /blk2  / bsycon
 
 c1   = CEXP(-ai*arg/2.0*(sps-sns))
 gam0 = sps*del - sigma
 pi2  = 2.0*pi
 s1   = sps/(dstr**2)
 s2   = sns/dstr
 c2q  = gam0/dstr - scrk
 c3q  = gam0/dstr + scrk
 nn   = 0
 csec = c2q*c3q
 IF (csec < 0.0) nn = 1
 t1   = gam0*s1
 t2   = s2*SQRT(ABS(csec))
 IF (c2q < 0.0 .AND. c3q < 0.0) t2 =-t2
 IF (nn == 0) alp0 = t1 + t2
 IF (nn == 1) alp0 = CMPLX(t1,t2)
 c1   = c1*(1.0-arg/alp0)
 a1   = pi2/(sps-sns)
 a2   =-a1
 b1   = gam0/(sps-sns)
 c1test = 0.0
 DO  i = 1,200
   r    = i
   gamp = pi2*r + gam0
   gamn =-pi2*r + gam0
   c2p  = gamp/dstr - scrk
   c2q  = gamp/dstr + scrk
   c2n  = gamn/dstr - scrk
   c3q  = gamn/dstr + scrk
   nn   = 0
   csec = c2p*c2q
   IF (csec < 0.0) nn = 1
   t1   = gamp*s1
   t2   = s2*SQRT(ABS(csec))
   IF (c2p < 0.0 .AND. c2q < 0.0) t2 =-t2
   IF (nn == 0) alp = t1 + t2
   IF (nn == 1) alp = CMPLX(t1,t2)
   nn   = 0
   csec = c2n*c3q
   IF (csec < 0.0) nn = 1
   t1   = gamn*s1
   t2   = s2*SQRT(ABS(csec))
   IF (c2n < 0.0 .AND. c3q < 0.0) t2 =-t2
   IF (nn == 0) aln = t1 + t2
   IF (nn == 1) aln = CMPLX(t1,t2)
   at2  = (alp-a1*r-b1)/(a1*r+b1-arg)
   at3  = (aln-a2*r-b1)/(a2*r+b1-arg)
   c1   = c1*(1.0+at2)*(1.0+at3)
   IF (cabs((c1-c1test)/c1) < 0.0009) GO TO 50
   c1test = c1
 END DO
 GO TO 70
 50 CONTINUE
 c1   = c1*b1/(arg-b1)*CSIN(pi/a1*(arg-b1))/(SIN(pi*b1/a1))
 c1   = c1*bsycon
 bkpm = c1
 RETURN
 
 70 WRITE  (ibbout,80) ufm
 80 FORMAT (a23,' - AMG MODULE -SUBROUTINE AKAPM')
 CALL mesage (-61,0,0)
 RETURN
END SUBROUTINE akapm
