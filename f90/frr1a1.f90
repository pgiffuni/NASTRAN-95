SUBROUTINE frr1a1 (rz,cz,ib,reb,ceb)
     
 
 
 REAL, INTENT(IN OUT)                     :: rz
 REAL, INTENT(IN OUT)                     :: cz
 INTEGER, INTENT(IN)                      :: ib
 REAL, INTENT(OUT)                        :: reb
 REAL, INTENT(OUT)                        :: ceb
 COMPLEX :: z,sum,zk,term
 
 z  = CMPLX(rz,cz)
 IF (cabs(z) < .1) GO TO 100
 zk  = CMPLX(1.,0.)
 n   = ib
 bf  = 1.
 bf1 = 0.
 sum = CMPLX(0.,0.)
 DO  i = 1,n
   sum = sum + zk/CMPLX(bf,0.)
   zk  = zk*z
   bf1 = bf1 + 1.
   bf  = bf*bf1
 END DO
 zk  = CMPLX(bf,0.)/zk*(CEXP(z)-sum)
 reb = REAL(zk)
 ceb = AIMAG(zk)
 RETURN
 
 100 CONTINUE
 zk  = z
 den = FLOAT(ib) + 1.
 sum = CMPLX(1.,0.)
 DO  i = 1,30
   term= zk/den
   sum = sum + term
   IF (cabs(term) < 1.e-9) EXIT
   zk  = zk*z
   den = den*(FLOAT(ib)+ FLOAT(i+1))
 END DO
 200 reb = REAL(sum)
 ceb = AIMAG(sum)
 RETURN
END SUBROUTINE frr1a1
