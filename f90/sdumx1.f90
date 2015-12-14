SUBROUTINE sdumx1
     
!     DELETE ANY OF THE FOLLOW ENTRY POINT IF A SUBROUTINE OF THE SAME
!     NAME ALREADY EXISTS
 
 INTEGER :: ii(9),kk(9)
 CHARACTER (LEN=25) :: uwm
 CHARACTER (LEN=23) :: ufm
 COMMON /xmssg / ufm,uwm
 COMMON /system/ ibuf,nout
 DATA    ii    / 9*0/,   jj /4HSDUM/,     kk /  &
     2H11,2H21,2H31,2H41,2H51,2H61,2H71,2H81,2H91 /
 
 GO TO 30
 
 
 ENTRY sdum91
!     ============
 
 j = 9
 GO TO 10
 
 
 ENTRY sdum81
!     ============
 
 j = 8
 GO TO 10
 
 
 ENTRY sdum71
!     ============
 
 j = 7
 GO TO 10
 
 
 ENTRY sdum61
!     ============
 
 j = 6
 GO TO 10
 
 
 ENTRY sdum51
!     ============
 
 j = 5
 GO TO 10
 
 
 ENTRY sdum41
!     ============
 
 j = 4
 GO TO 10
 
 
 ENTRY sdum31
!     ============
 
 j = 3
 GO TO 10
 
 
 ENTRY sdum21
!     ============
 
 j = 2
 GO TO 10
 
 
 ENTRY sdum11
!     ============
 
 j = 1
!     GO TO 10
 
 10 IF (ii(j) /= 0) GO TO 30
 ii(j)  = 1
 WRITE  (nout,20) uwm,jj,kk(j)
 20 FORMAT (a25,' 2182, SUBROUTINE ',2A4,' IS DUMMY.  ONLY ONE OF ',  &
     'THESE MESSAGES WILL APPEAR PER OVERLAY OF THIS DECK.')
 30 RETURN
END SUBROUTINE sdumx1
