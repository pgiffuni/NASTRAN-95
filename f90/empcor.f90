SUBROUTINE empcor(mt1x,mt2x,pt,pc,frsrow,midrow,lasrow,nx,a,z)
     
!     EMPTY CORE OF A TRIANGULAR MATRIX
 
 
 INTEGER, INTENT(IN)                      :: mt1x
 INTEGER, INTENT(IN)                      :: mt2x
 INTEGER, INTENT(IN)                      :: pt
 INTEGER, INTENT(IN)                      :: pc
 INTEGER, INTENT(IN)                      :: frsrow
 INTEGER, INTENT(IN)                      :: midrow
 INTEGER, INTENT(IN)                      :: lasrow
 INTEGER, INTENT(IN)                      :: nx
 REAL, INTENT(IN OUT)                     :: a(1)
 REAL, INTENT(IN OUT)                     :: z(1)
 INTEGER :: row,mcb(7)
 
 
 
!     MT1      FIRST PART OF THE MATRIX (UP TO ROW -MIDROW-).
!     MT2      REST OF THE MATRIX.
!     PT       PRECISION OF THE MATRIX ON TAPE.
!     PC       ......... .. ... ...... IN CORE.
!     FRSROW   FIRST ROW IN CORE.
!     LAST     LAST  ... .. CORE.
!     N        SIZE OF THE COMPLETE MATRIX.
!     A        LOCATION OF THE COMPLETE MATRIX.
 
 COMMON /packx/it1,it2,ii,jj,incr
 DATA  mcb /7*0/
 
 mt1 = mt1x
 mt2 = mt2x
 n   = nx
 mt  = mt1
 IF(frsrow > midrow .AND. mt2 /= 0) mt = mt2
 na  =1
 incr = 1
 it1 = pc
 it2 = pt
 jj  = n
 DO  row =  frsrow,lasrow
   ii = row
   CALL pack(a(na),mt,mcb)
   IF( row == n) GO TO 110
   na = na + pc* (n-row+1)
   IF( row /= midrow .OR. mt2 == 0) CYCLE
   CALL CLOSE(mt,1)
   mt = mt2
   CALL gopen(mt,z,1)
 END DO
 GO TO 115
 
!     END OF CORE DUMP
 
 110 CALL CLOSE(mt,1)
 115 RETURN
END SUBROUTINE empcor
