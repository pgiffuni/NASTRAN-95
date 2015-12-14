SUBROUTINE kdumx
     
!     DELETE ANY OF THE FOLLOW ENTRY POINT IF A SUBROUTINE OF THE SAME
!     NAME ALREADY EXISTS
 
 INTEGER :: ii(9),kk(9)
 CHARACTER (LEN=25) :: uwm
 CHARACTER (LEN=23) :: ufm
 COMMON /xmssg / ufm,uwm
 COMMON /system/ ibuf,nout
 DATA    ii    / 9*0/,   jj /4HKDUM/,     kk /  &
     1H1,1H2,1H3,1H4,1H5,1H6,1H7,1H8,1H9 /
 
 GO TO 30
 
 
 ENTRY kdum9
!     ===========
 
 j = 9
 GO TO 10
 
 
 ENTRY kdum8
!     ==========
 
 j = 8
 GO TO 10
 
 
 ENTRY kdum7
!     ==========
 
 j = 7
 GO TO 10
 
 
 ENTRY kdum6
!     ==========
 
 j = 6
 GO TO 10
 
 
 ENTRY kdum5
!     ==========
 
 j = 5
 GO TO 10
 
 
 ENTRY kdum4
!     ==========
 
 j = 4
 GO TO 10
 
 
 ENTRY kdum3
!     ==========
 
 j = 3
 GO TO 10
 
 
 ENTRY kdum2
!     ==========
 
 j = 2
 GO TO 10
 
 
 ENTRY kdum1
!     ==========
 
 j = 1
!     GO TO 10
 
 10 IF (ii(j) /= 0) GO TO 30
 ii(j)  = 1
 WRITE  (nout,20) uwm,jj,kk(j)
 20 FORMAT (a25,' 2182, SUBROUTINE ',2A4,' IS DUMMY.  ONLY ONE OF ',  &
     'THESE MESSAGES WILL APPEAR PER OVERLAY OF THIS DECK.')
 30 RETURN
END SUBROUTINE kdumx
