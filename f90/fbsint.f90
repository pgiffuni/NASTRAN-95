SUBROUTINE fbsint (x,y)
     
!     GIVEN THE DECOMPOSITION OF A REAL SYMMETRIC MATRIX, FBSINT WILL
!     SOLVE A SYSTEM OF SIMULTANEOUS LINEAR EQUATIONS BY FORWARD-
!     BACKWARD SUBSTITUTION
 
!     THIS ROUTINE IS SUITABLE FOR BOTH SINGLE AND DOUBLE PRECISION
!     OPERATION
 
 
 REAL, INTENT(IN)                         :: x(1)
 REAL, INTENT(OUT)                        :: y(1)
 INTEGER :: filel     ,iblk(15)
 
 COMMON /infbsx/ filel(7)
 COMMON /fbsx  / lfile(7)
 EQUIVALENCE     (filel(3) ,nrow)
 
 nrow2 = nrow
 IF (filel(5) == 2) nrow2 = 2*nrow
 DO  i = 1,nrow2
   y(i) = x(i)
 END DO
 DO  i = 1,7
   lfile(i) = filel(i)
 END DO
 CALL REWIND (filel)
 CALL skprec (filel,1)
 iblk(1) = filel(1)
 IF (filel(5) == 1) CALL fbs1 (iblk,y,y,nrow2)
 IF (filel(5) == 2) CALL fbs2 (iblk,y,y,nrow2)
 RETURN
END SUBROUTINE fbsint
