SUBROUTINE sdumx2
     
!     DELETE ANY OF THE FOLLOW ENTRY POINT IF A SUBROUTINE OF THE SAME
!     NAME ALREADY EXISTS
 
 INTEGER :: ii(9),kk(9)
 CHARACTER (LEN=25) :: uwm
 CHARACTER (LEN=23) :: ufm
 COMMON /xmssg / ufm,uwm
 COMMON /system/ ibuf,nout
 DATA    ii    / 9*0/,   jj /4HSDUM/,     kk /  &
     2H12,2H22,2H32,2H42,2H52,2H62,2H72,2H82,2H92 /
 
 GO TO 30
 
 
 ENTRY sdum92
!     ============
 
 j = 9
 GO TO 10
 
 
 ENTRY sdum82
!     ============
 
 j = 8
 GO TO 10
 
 
 ENTRY sdum72
!     ============
 
 j = 7
 GO TO 10
 
 
 ENTRY sdum62
!     ============
 
 j = 6
 GO TO 10
 
 
 ENTRY sdum52
!     ============
 
 j = 5
 GO TO 10
 
 
 ENTRY sdum42
!     ============
 
 j = 4
 GO TO 10
 
 
 ENTRY sdum32
!     ============
 
 j = 3
 GO TO 10
 
 
 ENTRY sdum22
!     ============
 
 j = 2
 GO TO 10
 
 
 ENTRY sdum12
!     ============
 
 j = 1
!     GO TO 10
 
 10 IF (ii(j) /= 0) GO TO 30
 ii(j)  = 1
 WRITE  (nout,20) uwm,jj,kk(j)
 20 FORMAT (a25,' 2182, SUBROUTINE ',2A4,' IS DUMMY.  ONLY ONE OF ',  &
     'THESE MESSAGES WILL APPEAR PER OVERLAY OF THIS DECK.')
 30 RETURN
END SUBROUTINE sdumx2
