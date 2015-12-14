INTEGER FUNCTION filcor(mt1x,mt2x,   pc,frsrow,midrow,nx,a,nza,z)
     
!     FILL CORE WITH A TRIANGULAR MATRIX
 
 
 INTEGER, INTENT(IN)                      :: mt1x
 INTEGER, INTENT(IN)                      :: mt2x
 INTEGER, INTENT(IN)                      :: pc
 INTEGER, INTENT(IN)                      :: frsrow
 INTEGER, INTENT(IN)                      :: midrow
 INTEGER, INTENT(IN)                      :: nx
 REAL, INTENT(OUT)                        :: a(1)
 INTEGER, INTENT(IN)                      :: nza
 REAL, INTENT(IN OUT)                     :: z(1)
 
 
 
 COMMON /unpakx/it1,ii,jj,incr1
 
!     MT1      FIRST PART OF THE MATRIX (UP TO ROW -MIDROW-).
!     MT2      REST OF THE MATRIX
!     PC       PRECISION OF THE MATRIX IN CORE
!     NX       COLUMN SIZE OF THE MATRIX
!     A        STORAGE FOR THE MATRIX
!     Z        BUFFER FOR GINO
!     FRSROW   FIRST ROW OF THE MATRIX TO BE READ
!     ANSWER   LAST ROW READ
 
 mt1 = mt1x
 mt2 = mt2x
 n   = nx
 mt = mt1
 lasrow = frsrow-1
 it1 = pc
 incr1 =1
 jj = n
 IF( lasrow >= midrow .AND. mt2 /= 0) mt = mt2
 nn = nza/pc
 na = 0
 
!     READ IN EACH ROW
 
 105 IF (na + n -lasrow > nn) GO TO 115
 lasrow  = lasrow +1
 i = pc*na +1
 ii = lasrow
 CALL unpack(*106,mt,a(i))
 GO TO 107
 106 k  = i +lasrow*pc-1
 DO  j = i,k
   a(j) =0.0
 END DO
 107 IF (lasrow == n) GO TO 110
 na = na + (n-lasrow +1)
 IF( lasrow /= midrow .OR. mt2 == 0) GO TO 105
 CALL CLOSE(mt,1)
 mt = mt2
 CALL gopen(mt,z,0)
 GO TO 105
 
!     END OF ROUTINE
 
 110 CALL CLOSE(mt,1)
 115 filcor = lasrow
 RETURN
END FUNCTION filcor
