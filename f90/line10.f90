SUBROUTINE line10 (x1,y1,x2,y2,penden,opt)
     
!     X1,Y1  = STARTING POINT OF THE LINE
!     X2,Y2  = TERMINAL POINT OF THE LINE
!     PENDEN = PEN NUMBER OR LINE DENSITY
!     OPT    = -1 TO INITIATE  THE LINE MODE
!            = +1 TO TERMINATE THE LINE MODE
!            = 0 TO DRAW A LINE
 
 
 REAL, INTENT(IN OUT)                     :: x1
 REAL, INTENT(IN OUT)                     :: y1
 REAL, INTENT(IN OUT)                     :: x2
 REAL, INTENT(IN OUT)                     :: y2
 INTEGER, INTENT(IN)                      :: penden
 INTEGER, INTENT(IN)                      :: opt
 INTEGER :: optx,a(6)
 DATA    optx,line / -1, 5  /
 
 IF (optx >= 0) optx = opt
 IF (opt < 0.0) THEN
   GO TO   200
 ELSE IF (opt == 0.0) THEN
   GO TO   100
 ELSE
   GO TO   150
 END IF
 100 a(1) = line
 a(2) = penden
 a(3) = IFIX(x1+.1)
 a(4) = IFIX(y1+.1)
 a(5) = IFIX(x2+.1)
 a(6) = IFIX(y2+.1)
 IF (optx == 0) GO TO 120
 
!     INITIATE THE LINE MODE.
 
 a(1) = a(1) + 10
 optx = 0
 
!     DRAW THE LINE.
 
 120 CALL wplt10 (a,0)
 GO TO 200
 
!     TERMINATE THE LINE MODE.
 
 150 CALL wplt10 (a,1)
 optx = -1
 
 200 RETURN
END SUBROUTINE line10
