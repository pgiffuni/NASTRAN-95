SUBROUTINE cfe1ao (tpose,v1,v2,v3,zb)
!*******
!     CFE1AO IS A SINGLE PRECISION ROUTINE WHICH PERFORMS THE OPERATION
!     (A) OR (A)-TRANSPOSE FOR THE COMPLEX FEER METHOD. THIS OPERATION
!     IS CALLED THE EIGENMATRIX MULTIPLICATION.
!*******
!     DEFINITION OF INPUT AND OUTPUT PARAMETERS
!*******
!     TPOSE    = .FALSE. --- PERFORM OPERATION (A)
!              = .TRUE.  --- PERFORM OPERATION (A)-TRANSPOSE
!     V1       = INPUT  VECTOR
!     V2       = OUTPUT VECTOR
!     V3       = INPUT WORKING SPACE (FOR INTERNAL USE)
!     ZB       = INPUT GINO BUFFER
!*******
 
 LOGICAL, INTENT(IN OUT)                  :: tpose(1)
 REAL, INTENT(IN)                         :: v1(1)
 REAL, INTENT(OUT)                        :: v2(1)
 REAL, INTENT(IN OUT)                     :: v3(1)
 REAL, INTENT(IN OUT)                     :: zb(1)
 
 LOGICAL :: no b     ,qpr
 REAL :: lambda
 COMMON  /feeraa/  ik(7)    ,im(7)    ,ib(7)    ,dumaa(117) ,mcblmb(7)
 COMMON  /feerxc/  lambda(4),dum01(2) ,nord     ,idiag  &
     ,epsdum(2),northo   ,nord2    ,nord4  &
     ,nordp1   ,dum02(2) ,no b     ,dum03(4) ,qpr
 COMMON  /system/  ksys     ,nout
!*******
 IF (qpr) WRITE (nout,8881) tpose
 8881 FORMAT(1H0,12HENTER cfe1ao,8X,11HTRANSPOSE =,l2)
 IF (tpose(1)) GO TO 50
!*******
!     PERFORM OPERATION (A)  = EIGENMATRIX MULTIPLICATION
!*******
 IF ( no b ) GO TO 30
!*******
!     MULTIPLY LOWER HALF OF INPUT VECTOR BY MASS MATRIX
!*******
 CALL cfe1my (tpose(1),v1(nordp1),v3(1),im(1),zb(1))
 IF (qpr) WRITE (nout,8882) (v3(i),i=1,nord)
 8882 FORMAT(3H --,32(4H----),/(1H ,6E21.13))
!*******
!     MULTIPLY UPPER HALF OF INPUT VECTOR BY -(LAMBDA*M+B)
!*******
 CALL cfe1my (tpose(1),v1(1),v3(nordp1),mcblmb(1),zb(1))
 IF (qpr) WRITE (nout,8882) (v3(i),i=nordp1,nord2)
!*******
!     CALCULATE RIGHT-HAND SIDE OF SWEEP EQUATION
!*******
 DO   i = 1,nord
   j = nord+i
   v2(i) = -v3(i) + v3(j)
 END DO
 IF (qpr) WRITE (nout,8882) (v2(i),i=1,nord)
!*******
!     PERFORM FORWARD AND BACKWARD SWEEPS
!     (GENERATES UPPER HALF OF OUTPUT VECTOR)
!*******
 CALL cf1fbs (tpose(1),v2(1),zb(1))
 IF (qpr) WRITE (nout,8882) (v2(i),i=1,nord)
!*******
!     COMPUTE LOWER HALF OF OUTPUT VECTOR
!*******
 DO   i = 1,nord,2
   j = i+1
   ni = nord+i
   nj = ni+1
   v2(ni) = v1(i) + lambda(1)*v2(i) - lambda(3)*v2(j)
   v2(nj) = v1(j) + lambda(1)*v2(j) + lambda(3)*v2(i)
 END DO
 IF (qpr) WRITE (nout,8882) (v2(i),i=nordp1,nord2)
 GO TO 200
!*******
!     DAMPING MATRIX ABSENT
!*******
!     MULTIPLY INPUT VECTOR BY MASS MATRIX
!*******
 30 CALL cfe1my (tpose(1),v1(1),v2(1),im(1),zb(1))
 DO   i = 1,nord2
   v2(i) = -v2(i)
 END DO
 IF (qpr) WRITE (nout,8882) (v2(i),i=1,nord2)
!*******
!     PERFORM FORWARD AND BACKWARD SWEEPS
!*******
 CALL cf1fbs (tpose(1),v2(1),zb(1))
 IF (qpr) WRITE (nout,8882) (v2(i),i=1,nord2)
 GO TO 200
!*******
!     PERFORM OPERATION (A)-TRANSPOSE  = TRANSPOSED EIGENMATRIX
!                                                   MULTIPLICATION
!*******
 50 IF ( no b ) GO TO 90
!*******
!     CALCULATE RIGHT-HAND SIDE OF SWEEP EQUATION
!*******
 DO   i = nordp1,nord2,2
   j = i+1
   ni = i-nord
   nj = ni+1
   v3(i) = v1(ni) + lambda(1)*v1(i) - lambda(3)*v1(j)
   v3(j) = v1(nj) + lambda(1)*v1(j) + lambda(3)*v1(i)
 END DO
 IF (qpr) WRITE (nout,8882) (v3(i),i=nordp1,nord2)
!*******
!     PERFORM BACKWARD AND FORWARD SWEEPS
!*******
 CALL cf1fbs (tpose(1),v3(nordp1),zb(1))
 IF (qpr) WRITE (nout,8882) (v3(i),i=nordp1,nord2)
!*******
!     MULTIPLY SWEEP OUTPUT VECTOR BY -(LAMBDA*M+B)-TRANSPOSE
!*******
 CALL cfe1my (tpose(1),v3(nordp1),v3(1),mcblmb(1),zb(1))
 IF (qpr) WRITE (nout,8882) (v3(i),i=1,nord)
!*******
!     COMPUTE UPPER HALF OF OUTPUT VECTOR
!*******
 DO   i = 1,nord
   j = nord+i
   v2(i) = v1(j) + v3(i)
 END DO
 IF (qpr) WRITE (nout,8882) (v2(i),i=1,nord)
!*******
!     MULTIPLY SWEEP OUTPUT VECTOR BY TRANSPOSED MASS MATRIX
!     (GENERATES NEGATIVE OF LOWER HALF OF OUTPUT VECTOR)
!*******
 CALL cfe1my (tpose(1),v3(nordp1),v2(nordp1),im(1),zb(1))
 DO   i = nordp1,nord2
   v2(i) = -v2(i)
 END DO
 IF (qpr) WRITE (nout,8882) (v2(i),i=nordp1,nord2)
 GO TO 200
!*******
!     DAMPING MATRIX ABSENT
!*******
!     PERFORM BACKWARD AND FORWARD SWEEPS
!*******
 90 DO   i = 1,nord2
   v3(i) = v1(i)
 END DO
 CALL cf1fbs (tpose(1),v3(1),zb(1))
 IF (qpr) WRITE (nout,8882) (v3(i),i=1,nord2)
!*******
!     MULTIPLY SWEEP OUTPUT VECTOR BY TRANSPOSED MASS MATRIX
!*******
 CALL cfe1my (tpose(1),v3(1),v2(1),im(1),zb(1))
 DO   i = 1,nord2
   v2(i) = -v2(i)
 END DO
 IF (qpr) WRITE (nout,8882) (v2(i),i=1,nord2)
 200 RETURN
END SUBROUTINE cfe1ao
