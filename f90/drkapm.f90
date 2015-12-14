SUBROUTINE drkapm (arg,indx,reslt)
     
!     THIS SUBROUTINE COMPUTES THE DERVIATIVE OF KAPPA MINUS
 
 
 COMPLEX, INTENT(IN)                      :: arg
 INTEGER, INTENT(IN)                      :: indx
 COMPLEX, INTENT(OUT)                     :: reslt
 COMPLEX :: ai, bsycon,c1,c2,c2test,at2,at3,alp0, alp,aln
 CHARACTER (LEN=23) :: ufm
 COMMON /xmssg / ufm
 COMMON /system/ sysbuf,ibbout
 COMMON /blk1  / scrk,sps,sns,dstr,ai,pi,del,sigma,beta,res
 COMMON /blk2  / bsycon
 
 pi2  = 2.0*pi
 a1   = pi2/(sps-sns)
 a2   =-a1
 gam0 = sps*del - sigma
 b1   = gam0/(sps-sns)
 c1   = CEXP(-ai*arg/2.0*(sps-sns))
 c2q  = gam0/dstr - scrk
 c3q  = gam0/dstr + scrk
 s1   = sps/(dstr**2)
 s2   = sns/dstr
 nn   = 0
 csec = c2q*c3q
 IF (csec < 0.0) nn = 1
 t1   = gam0*s1
 t2   = s2*SQRT(ABS(csec))
 IF (c2q < 0.0 .AND. c3q < 0.0) t2 =-t2
 IF (nn == 0) alp0 = t1 + t2
 IF (nn == 1) alp0 = CMPLX(t1,t2)
 rindx = indx
 IF (indx == 0) GO TO 10
 c2   = c1*b1/alp0*CSIN(pi/a1*(arg-b1))/(a1*rindx+b1-arg)*  &
     (1.0+(alp0-b1)/(b1-arg))/(SIN(pi*b1/a1))*bsycon
 GO TO 20
 10 CONTINUE
 c2   = c1*b1/alp0*CSIN(pi/a1*(arg-b1))/((b1-alp0)*SIN(pi*b1/a1))* bsycon
 20 CONTINUE
 c2test = 0.0
 DO  i = 1,200
   r    = i
   IF (indx < 0 .AND. ABS(rindx) == r) CYCLE
   IF (indx > 0 .AND. rindx == r) CYCLE
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
   c2   = c2*(1.0+at2)*(1.0+at3)
   IF (cabs((c2-c2test)/c2) < 0.0009) GO TO 40
   c2test = c2
 END DO
 GO TO 50
 40 CONTINUE
 reslt = c2
 RETURN
 
 50 CONTINUE
 WRITE  (ibbout,60) ufm
 60 FORMAT (a23,' - AMG MODULE -SUBROUTINE DRKAPM')
 CALL mesage (-61,0,0)
 RETURN
END SUBROUTINE drkapm
