SUBROUTINE angtrs (theta,k,trans)
!     &    ENTRY ANGTRD
 
!     ROUTINE TO CALCULATE AND OUTPUT THE INPLANE ROTATION
!     TRANSFORMATION IN 3-D USING THE ANGLE OF ROTATION.
 
!     IF K=1, TRANS OR TRAND WILL BE TRANSPOSED AND THEN RETURNED.
 
!     SINGLE PRECISION -
 
 
 REAL, INTENT(IN OUT)                     :: theta
 INTEGER, INTENT(IN OUT)                  :: k
 REAL, INTENT(OUT)                        :: trans(9)
 
 DOUBLE PRECISION :: trand(9), thetad
 
 DO  i = 1,9
   trans(i) = 0.0
 END DO
 
 trans(1) =  COS(theta)
 trans(2) =  SIN(theta)
 trans(4) = -trans(2)
 trans(5) =  trans(1)
 trans(9) =  1.0
 
 IF (k /= 1) GO TO 30
 trans(2) = -trans(2)
 trans(4) = -trans(4)
 RETURN
 
 ENTRY angtrd (thetad,k,trand)
!     =============================
 
!     DOUBLE PRECISION -
 
 DO  i = 1,9
   trand(i) = 0.0D0
 END DO
 
 trand(1) =  DCOS(thetad)
 trand(2) =  DSIN(thetad)
 trand(4) = -trand(2)
 trand(5) =  trand(1)
 trand(9) =  1.0D0
 
 IF (k /= 1) GO TO 30
 trand(2) = -trand(2)
 trand(4) = -trand(4)
 30 RETURN
END SUBROUTINE angtrs
