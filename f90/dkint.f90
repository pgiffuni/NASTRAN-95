DOUBLE PRECISION FUNCTION dkint(i,j,a,b,iv,iw,r,z)
     
 INTEGER, INTENT(IN OUT)                  :: i
 INTEGER, INTENT(IN OUT)                  :: j
 DOUBLE PRECISION, INTENT(IN)             :: a
 DOUBLE PRECISION, INTENT(IN)             :: b
 INTEGER, INTENT(IN)                      :: iv
 INTEGER, INTENT(IN)                      :: iw
 DOUBLE PRECISION, INTENT(IN)             :: r(1)
 DOUBLE PRECISION, INTENT(IN OUT)         :: z(1)
 DOUBLE PRECISION :: bint, c1p, c2p, c1, c2, aw, dkj, dkef
 DOUBLE PRECISION :: sp1
 
 
 bint = 0.0D0
 iw1 = iw + 1
 c1p = b
 c2p = a
 c1 = c1p
 c2 = c2p
 aw = 0.0D0
 IF( r(i) /= 0.0D0 .AND. r(j) /= 0.0D0 ) aw = DLOG(r(j)/r(i))
 DO  it = 1,iw1
   ic = iw - it + 1
   IF (ic == 0) c1 = 1.0D0
   IF (it == 1) c2 = 1.0D0
   
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   
! THE FOLLOWING CODE REPLACES DOUBLE PRECISION FUNCTION DKEF
   
   IF(it == 1) GO TO 20
   in=1
   id=1
   DO  k=2,it
     in=in*(iw-k+2)
     id=id*(k-1)
   END DO
   dkef=in/id
   GO TO 30
   20 dkef=1.0D0
   30 CONTINUE
   
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   
   
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   
! THE FOLLOWING CODE REPLACES DOUBLE PRECISION FUNCTION DKJ
   
   is1 = ic+iv+1
   IF(is1 == 0) GO TO 60
   sp1=is1
   dkj=( r(j)**is1 - r(i)**is1 )  /  sp1
   GO TO 70
   60 dkj = aw
   70 CONTINUE
   
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   
   bint = bint + c1 ** ic *dkj * c2 ** (it - 1) * dkef
   c1 = c1p
   c2 = c2p
 END DO
 aw = iw
 bint = bint / aw
 dkint = bint
 RETURN
END FUNCTION dkint
