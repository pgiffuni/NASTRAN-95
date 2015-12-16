SUBROUTINE matvc2 (y,x,filea,buf)
     
!     MATVC2 WILL FORM THE PRODUCT X = X + A*Y WHERE A IS A MATRIX
!     AND Y IS A VECTOR
 
!     THIS ROUTINE IS SUITABLE FOR DOUBLE PRECISION OPERATION
 
 
 DOUBLE PRECISION, INTENT(IN)             :: y(1)
 DOUBLE PRECISION, INTENT(OUT)            :: x(1)
 INTEGER, INTENT(IN)                      :: filea(7)
 REAL, INTENT(IN OUT)                     :: buf(1)
 INTEGER :: sub(2)   ,diag     ,rdp     ,eol
 DOUBLE PRECISION :: a        ,da
 
 COMMON   /zntpkx/  a(2)      ,ii       ,eol
!     COMMON   /DESCRP/  LENGTH    ,MAJOR(1)
 COMMON   /names /  rd        ,rdrew    ,wrt      ,wrtrew   ,  &
     rew       ,norew    ,eofnrw   ,rsp      ,  &
     rdp       ,csp      ,cdp      ,sqr      ,  &
     rect      ,diag     ,iwtri    ,uprtri   , sym       ,row      ,identy
 COMMON   /trdxx /  idum(27)  ,iopen
 EQUIVALENCE        (a(1),da)
 DATA      sub   /  4HMATV, 4HC2   /
 
 IF (filea(1) == 0) RETURN
 ncol = filea(2)
 IF (filea(4) == identy) GO TO 60
 IF (iopen == 1) GO TO 5
 CALL OPEN (*90,filea(1),buf,rdrew)
 5 CALL fwdrec (*100,filea(1))
 IF (filea(4) == diag) GO TO 40
 
!     MATRIX IS FULL
 
 DO  i = 1,ncol
   IF (y(i) == 0.0D0) GO TO 20
   CALL intpk (*30,filea(1),0,rdp,0)
   10 CALL zntpki
   x(ii) = da*y(i) + x(ii)
   IF (eol == 0.0) THEN
     GO TO    10
   ELSE
     GO TO    30
   END IF
   20 CALL fwdrec (*100,filea(1))
   30 CONTINUE
 END DO
 GO TO 80
 
!     MATRIX IS DIAGONAL
 
 40 CALL intpk (*80,filea(1),0,rdp,0)
 50 CALL zntpki
 x(ii) = y(ii)*da +x(ii)
 IF (eol == 0.0) THEN
   GO TO    50
 ELSE
   GO TO    80
 END IF
 
!     MATRIX IS THE IDENTITY
 
 60 DO  i = 1,ncol
   x(i) = y(i) + x(i)
 END DO
 RETURN
 
 80 CALL REWIND (filea(1))
 IF (iopen == 0) CALL CLOSE (filea(1),rew)
 RETURN
 
 90 no = -1
 GO TO 110
 100 no = -2
 110 CALL mesage (no,filea(1),sub(1))
 RETURN
END SUBROUTINE matvc2
