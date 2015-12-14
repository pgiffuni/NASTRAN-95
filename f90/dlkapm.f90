SUBROUTINE dlkapm (arg,blkapm)
     
!     SUBROUTINE FOR COMPUTING LOGARITHMIC DERIVATIVE OF KAPPA MINUS
 
 
 COMPLEX, INTENT(IN)                      :: arg
 COMPLEX, INTENT(OUT)                     :: blkapm
 COMPLEX :: ai,c1,d1,d2,c1test, e1,alp0,alp,aln
 CHARACTER (LEN=23) :: ufm
 COMMON /xmssg / ufm
 COMMON /system/ sysbuf,ibbout
 COMMON /blk1  / scrk,sps,sns,dstr,ai,pi,del,sigma,beta,res
 
 c1   =-ai/2.0*(sps-sns)
 pi2  = 2.0*pi
 s1   = sps/(dstr**2)
 s2   = sns/dstr
 gam0 = sps*del - sigma
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
 c1   = c1 + 1.0/(arg-alp0)
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
   e1   = a1*r + b1 - arg
   d1   = (alp-a1*r-b1)/e1
   d2   = d1/e1
   c1   = c1 + 1.0/(1.0+d1)*d2
   e1   = a2*r + b1 - arg
   d1   = (aln-a2*r-b1)/e1
   d2   = d1/e1
   c1   = c1 + 1.0/(1.0+d1)*d2
   IF (cabs((c1-c1test)/c1) < 0.0006) GO TO 50
   c1test = c1
 END DO
 GO TO 70
 50 CONTINUE
 e1   = arg - b1
 b    = pi/a1
 c1   = c1 - 1.0/e1 + b*CCOS(b*e1)/(CSIN(b*e1))
 blkapm = c1
 RETURN
 
 70 WRITE  (ibbout,80) ufm
 80 FORMAT (a23,' - AMG MODULE -SUBROUTINE DLKAPM')
 CALL mesage (-61,0,0)
 RETURN
END SUBROUTINE dlkapm
