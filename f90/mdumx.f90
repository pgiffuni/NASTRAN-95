SUBROUTINE mdumx
     
!     DELETE ANY OF THE FOLLOW ENTRY POINT IF A SUBROUTINE OF THE SAME
!     NAME ALREADY EXISTS
 
 INTEGER :: ii(9),kk(9)
 CHARACTER (LEN=25) :: uwm
 CHARACTER (LEN=23) :: ufm
 COMMON /xmssg / ufm,uwm
 COMMON /system/ ibuf,nout
 DATA    ii    / 9*0/,   jj /4HMDUM/,     kk /  &
     1H1,1H2,1H3,1H4,1H5,1H6,1H7,1H8,1H9 /
 
 GO TO 30
 
 
 ENTRY mdum9
!     ===========
 
 j = 9
 GO TO 10
 
 
 ENTRY mdum8
!     ==========
 
 j = 8
 GO TO 10
 
 
 ENTRY mdum7
!     ==========
 
 j = 7
 GO TO 10
 
 
 ENTRY mdum6
!     ==========
 
 j = 6
 GO TO 10
 
 
 ENTRY mdum5
!     ==========
 
 j = 5
 GO TO 10
 
 
 ENTRY mdum4
!     ==========
 
 j = 4
 GO TO 10
 
 
 ENTRY mdum3
!     ==========
 
 j = 3
 GO TO 10
 
 
 ENTRY mdum2
!     ==========
 
 j = 2
 GO TO 10
 
 
 ENTRY mdum1
!     ==========
 
 j = 1
!     GO TO 10
 
 10 IF (ii(j) /= 0) GO TO 30
 ii(j)  = 1
 WRITE  (nout,20) uwm,jj,kk(j)
 20 FORMAT (a25,' 2182, SUBROUTINE ',2A4,' IS DUMMY.  ONLY ONE OF ',  &
     'THESE MESSAGES WILL APPEAR PER OVERLAY OF THIS DECK.')
 30 RETURN
END SUBROUTINE mdumx
