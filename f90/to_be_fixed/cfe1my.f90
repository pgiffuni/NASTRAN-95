SUBROUTINE cfe1my (tpose,y,x,FILE,buf)
!*******
!     CFE1MY FORMS THE COMPLEX SINGLE PRECISION MATRIX
!     PRODUCT X = M*Y FOR THE COMPLEX FEER METHOD
!*******
!     DEFINITION OF INPUT AND OUTPUT PARAMETERS
!*******
!     TPOSE    = .FALSE. --- USE MATRIX M
!              = .TRUE.  --- USE MATRIX M-TRANSPOSE
!     Y        = INPUT  VECTOR
!     X        = OUTPUT VECTOR
!     FILE     = INPUT MATRIX CONTROL BLOCK FOR THE
!                REQUIRED MATRIX
!     BUF      = INPUT REQUIRED GINO BUFFER
!*******
 
 LOGICAL, INTENT(IN OUT)                  :: tpose(1)
 REAL, INTENT(IN)                         :: y(1)
 REAL, INTENT(OUT)                        :: x(1)
 INTEGER, INTENT(IN)                      :: FILE(7)
 INTEGER, INTENT(IN OUT)                  :: buf(1)
 
 INTEGER :: eol      ,diag
 
 COMMON  /names /  rd       ,rdrew    ,wrt      ,wrtrew  &
     ,rew      ,norew    ,eofnrw   ,rsp ,rdp      ,csp      ,cdp      ,sqr  &
     ,rect     ,diag     ,lowtri   ,uprtri ,sym      ,row      ,identy
 COMMON  /zntpkx/  da(4)    ,ii       ,eol
 
 ncol2 = FILE(2)+FILE(2)
 IF (FILE(4) == identy) GO TO 50
 CALL gopen (FILE(1),buf(1),rdrew)
 DO   i = 1,ncol2
   x(i) = 0.
 END DO
 IF (FILE(4) == diag) GO TO 40
 IF (tpose(1)) GO TO 31
!*******
!     GENERAL MATRIX*VECTOR PRODUCT
!*******
 DO   i = 1,ncol2,2
   j = i+1
   IF (y(i) == 0..AND.y(j) == 0.) GO TO 25
   CALL intpk(*30,FILE(1),0,csp,0)
   22 CALL zntpki
   jj = ii+ii
   ii = jj-1
   x(ii) = x(ii) + da(1)*y(i) - da(2)*y(j)
   x(jj) = x(jj) + da(1)*y(j) + da(2)*y(i)
   IF (eol == 0.0) THEN
     GO TO    22
   ELSE
     GO TO    30
   END IF
   25 CALL skprec (FILE(1),1)
 END DO
 GO TO 80
!*******
!     GENERAL MATRIX-TRANSPOSE*VECTOR PRODUCT
!*******
 31 DO   i = 1,ncol2,2
   j = i+1
   CALL intpk(*36,FILE(1),0,csp,0)
   32 CALL zntpki
   jj = ii+ii
   ii = jj-1
   x(i) = x(i) + da(1)*y(ii) - da(2)*y(jj)
   x(j) = x(j) + da(1)*y(jj) + da(2)*y(ii)
   IF (eol == 0.0) THEN
     GO TO    32
   ELSE
     GO TO    36
   END IF
 END DO
 GO TO 80
!*******
!     MATRIX IS DIAGONAL
!*******
 40  CALL intpk(*80,FILE(1),0,csp,0)
 45 CALL zntpki
 jj = ii+ii
 ii = jj-1
 x(ii) = da(1)*y(ii) - da(2)*y(jj)
 x(jj) = da(1)*y(jj) + da(2)*y(ii)
 IF (eol == 0.0) THEN
   GO TO    45
 ELSE
   GO TO    80
 END IF
!*******
!     MATRIX IS IDENTITY
!*******
 50 DO   i = 1,ncol2
   x(i) = y(i)
 END DO
 GO TO 90
 80 CALL CLOSE (FILE(1),rew)
 90 RETURN
END SUBROUTINE cfe1my
