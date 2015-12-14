SUBROUTINE axis10 (x1,y1,x2,y2,penden,opt)
     
!     (X1,Y1) = STARTING POINT OF THE AXIS.
!     (X2,Y2) = TERMINAL POINT OF THE AXIS.
!     PENDEN  = PEN NUMBER OR LINE DENSITY.
!     OPT     = -1 TO INITIATE  THE AXIS MODE.
!     ...     = +1 TO TERMINATE THE AXIS MODE.
!     ...     =  0 TO DRAW AN AXIS.
 
 
 REAL, INTENT(IN)                         :: x1
 REAL, INTENT(IN)                         :: y1
 REAL, INTENT(IN)                         :: x2
 REAL, INTENT(IN)                         :: y2
 INTEGER, INTENT(IN)                      :: penden
 INTEGER, INTENT(IN)                      :: opt
 INTEGER :: optx,a(6),axis
 REAL :: xy(2,2)
 COMMON /pltdat/ skpplt(2),xymin(2),xymax(2)
 DATA    optx  / -1 /
 DATA    axis  /  6 /
 
 IF (optx >= 0) optx = opt
 IF (opt < 0.0) THEN
   GO TO   200
 ELSE IF (opt == 0.0) THEN
   GO TO   100
 ELSE
   GO TO   150
 END IF
 100 xy(1,1) = x1
 xy(2,1) = y1
 xy(1,2) = x2
 xy(2,2) = y2
 DO  j = 1,2
   DO  i = 1,2
     IF (xy(i,j) < xymin(i)) xy(i,j) = xymin(i)
     IF (xy(i,j) > xymax(i)) xy(i,j) = xymax(i)
   END DO
 END DO
 
!     DRAW THE AXIS.
 
 a(1) = axis
 a(2) = penden
 DO  j = 1,2
   a(2*j+1) = xy(1,j) + .1
   a(2*j+2) = xy(2,j) + .1
 END DO
 IF (optx == 0) GO TO 120
 
!     INITIATE THE AXIS MODE.
 
 a(1) = a(1) + 10
 optx = 0
 
!     DRAW THE LINE.
 
 120 CALL wplt10 (a,0)
 GO TO 200
 
 
!     TERMINATE THE LINE MODE.
 
 150 CALL wplt10 (a,1)
 optx = -1
 
 200 RETURN
END SUBROUTINE axis10
