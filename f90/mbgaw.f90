SUBROUTINE mbgaw(boxl,dphi,ws,paw,paf1,paf2,q,q1,q2,m,kc,kc1,kc2)
     
!     MAIN PLANE BOXES
!     (NEW MSC METHOD USED)
 
 
 
 REAL, INTENT(IN)                         :: boxl
 COMPLEX, INTENT(IN)                      :: dphi
 COMPLEX, INTENT(OUT)                     :: ws
 REAL, INTENT(IN)                         :: paw
 REAL, INTENT(IN)                         :: paf1
 REAL, INTENT(IN)                         :: paf2
 COMPLEX, INTENT(OUT)                     :: q(1)
 COMPLEX, INTENT(OUT)                     :: q1(1)
 COMPLEX, INTENT(OUT)                     :: q2(1)
 INTEGER, INTENT(IN OUT)                  :: m
 INTEGER, INTENT(IN OUT)                  :: kc
 INTEGER, INTENT(IN OUT)                  :: kc1
 INTEGER, INTENT(IN OUT)                  :: kc2
 
 ws  =  (-0.5 * AMIN0(m,2) * boxl ) * dphi
 IF( paw < 0.005 ) GO TO 120
 q(kc) = paw * ws
 120 IF( paf1 < 0.005 ) GO TO 140
 q1(kc1) = paf1 * ws
 140 IF( paf2 < 0.005 ) GO TO 300
 q2(kc2) = paf2 * ws
 300 CONTINUE
 RETURN
END SUBROUTINE mbgaw
