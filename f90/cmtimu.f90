SUBROUTINE cmtimu (y,x,FILE,buf)
     
!     CM TIM U FORMS THE MATRIX PRODUCT X = M*Y WHERE ALL MAY BE COMPLEX
 
 
 DOUBLE PRECISION, INTENT(IN)             :: y(1)
 DOUBLE PRECISION, INTENT(OUT)            :: x(1)
 INTEGER, INTENT(IN)                      :: FILE(1)
 INTEGER, INTENT(IN OUT)                  :: buf(1)
 INTEGER :: diag      ,eol      ,eor      ,filem(7) ,  &
     filek     , filemm   , NAME(2)
 DOUBLE PRECISION :: da
 COMMON   /cinvpx/  filek(7)  ,filemm(7)
 COMMON   /names /  rd        ,rdrew    ,wrt      ,wrtrew   ,  &
     rew       ,norew    ,eofnrw   ,rsp      ,  &
     rdp       ,csp      ,cdp      ,sqr      ,  &
     rect      ,diag     ,lowtri   ,uprtri   , sym       ,row      ,identy
 COMMON   /cinvxx/  dum(21)   ,nzero
 COMMON   /zntpkx/  da(2)     ,ii       ,eol      ,eor
!     COMMON   /DESCRP/  LENGTH    ,MAJOR(1)
 EQUIVALENCE        (ncol,filek(2))
 DATA      NAME  /  4HCMTI    ,4HMU     /
 
 IF (FILE(1) == 0) GO TO 5
 
!     USE MATRIX OTHER THAN THE MASS MATRIX
 
 DO  i = 1,7
   filem(i) = FILE(i)
 END DO
 GO TO 8
 
!     USE MASS MATRIX
 
 5 DO  i = 1,7
   filem(i) = filemm(i)
 END DO
 8 CONTINUE
 ncol2 = ncol + ncol
 IF (filem(4) == identy) GO TO 50
 nzero = 0
 CALL gopen (filem(1),buf,rdrew)
 DO  i = 1,ncol2
   x(i) = 0.d0
 END DO
 IF (filem(4) == diag) GO TO 40
 
!     MASS MATRIX IS NOT DIAGONAL
 
 DO  i = 1,ncol2,2
   IF (y(i) == 0.d0 .AND. y(i+1) == 0.d0) GO TO 25
   CALL intpk (*30,filem(1),0,cdp,0)
   22 CALL zntpki
   IF (ii == i) nzero = nzero + 1
   ii = ii+ii-1
   x(ii  ) = x(ii  ) + da(1)*y(i  )-da(2)*y(i+1)
   x(ii+1) = x(ii+1) + da(1)*y(i+1)+da(2)*y(i  )
   IF (eol == 0) IF (eor) 30,22,30
   CYCLE
   25 CALL fwdrec (*80,filem(1))
   30 CONTINUE
 END DO
 GO TO 80
 
!     FILE ERROR
 
!  35 J = -1
!     GO TO 37
!  36 J = -2
!  37 CALL MESAGE (J,FILEM(1),NAME)
 
!     MASS MATRIX IS DIAGONAL
 
 40 CALL intpk (*80,filem(1),0,cdp,0)
 45 CALL zntpki
 ii = ii + ii - 1
 x(ii  ) = y(ii)*da(1) - y(ii+1)*da(2)
 x(ii+1) = y(ii)*da(2) + y(ii+1)*da(1)
 nzero = nzero + 1
 IF (eol == 0) IF (eor) 80,45,80
 GO TO 80
 
!     MASS MATRIX IS THE IDENTY
 
 50 DO  i = 1,ncol2
   x(i) = y(i)
 END DO
 nzero = 0
 RETURN
 
 80 CALL CLOSE (filem(1),rew)
 nzero = 0
 RETURN
END SUBROUTINE cmtimu
