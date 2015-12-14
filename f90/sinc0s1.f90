SUBROUTINE sinc0s1 (row,sick, d,o,COS)
!                    =
!     SUBROUTINE SICOX (D,O,COS)
 
!     THIS ROUTINE WAS CALLED SICOX BEFORE, WITH ENTRY POINT SINCAS
!                                                                =
!     THIS ROUTINE IS CALLED ONLY BY TRIDI SUBROUTINE, WHICH IS CALLED
!     ONLY BY VALVEC
 
!     IT CALCULATES SINES AND COSINES FOR GIVENS TRIDIAGONALIZATION
 
 
 INTEGER, INTENT(IN)                      :: row
 INTEGER, INTENT(OUT)                     :: sick
 REAL, INTENT(IN OUT)                     :: d(1)
 REAL, INTENT(OUT)                        :: o(1)
 REAL, INTENT(OUT)                        :: COS(1)
 INTEGER :: rowp2
 REAL :: z
 COMMON /givn  /  title(100),n
 
!     D   = DIAGONAL AND SINES.
!     O   = OFF-DIAGONAL.
!     COS = COSINES.
 
!     RETURN
 
 
!     ENTRY SINCAS (ROW,SICK)
!     =======================
 
!     CALCULATE THE SINES AND COSINES OF ROW -ROW-.
 
 sick  = 0
 rowp2 = row + 2
 DO  i = rowp2,n
   IF (d(i) == 0.0) GO TO 101
   
!     CALCULATE THE ROTATION.
   
   sick = 1
   z    = SQRT(d(i)**2 + d(row+1)**2)
   d(i) = d(i)/z
   COS(i) = d(row+1)/z
   d(row+1) = z
   IF (COS(i) >= 0.0) CYCLE
   COS(i) = ABS(COS(i))
   d(i)   = -d(i)
   d(row+1) = -d(row+1)
   CYCLE
   
!     NO ROTATION.
   
   101 COS(i) = 1.0
 END DO
 o(row) = d(row+1)
 RETURN
END SUBROUTINE sinc0s1
