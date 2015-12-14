SUBROUTINE mbgate(ntote,dphite,n,ywte,q,q1,q2,kte,kte1,kte2)
     
!     SUM ON TRAILING EDGE
 
 
 INTEGER, INTENT(IN)                      :: ntote
 COMPLEX, INTENT(IN)                      :: dphite(3,n)
 INTEGER, INTENT(IN OUT)                  :: n
 REAL, INTENT(IN OUT)                     :: ywte(1)
 COMPLEX, INTENT(OUT)                     :: q(1)
 COMPLEX, INTENT(OUT)                     :: q1(1)
 COMPLEX, INTENT(OUT)                     :: q2(1)
 INTEGER, INTENT(IN)                      :: kte(1)
 INTEGER, INTENT(IN)                      :: kte1(1)
 INTEGER, INTENT(IN)                      :: kte2(1)
 LOGICAL :: cntrl1,cntrl2
 
 COMPLEX :: dphi
 COMMON /mboxc/ njj ,crank1,crank2,cntrl1,cntrl2
 COMMON /mboxa/ x(12),y(12)
 
 DO      j = 1 , ntote
   dphi  =  dphite(1,j) * 0.5 * AMIN0(j,2)
   IF(cntrl1.AND.ywte(j) >= y(7).AND.ywte(j) <= y(11)) GO TO 1100
   IF(cntrl2.AND.ywte(j) > y(11).AND.ywte(j) <= y(12)) GO TO 1150
   isp=kte(j)
   q(isp) = dphi
   GO TO 1300
   1100 isp=kte1(j)
   q1(isp) = dphi
   GO TO 1300
   1150 isp=kte2(j)
   q2(isp) = dphi
   1300 CONTINUE
 END DO
 RETURN
END SUBROUTINE mbgate
