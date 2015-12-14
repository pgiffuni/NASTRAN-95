SUBROUTINE line (x1,y1,x2,y2,penx,opt)
     
!     (X1,Y1) = STARTING POINT OF THE LINE
!     (X2,Y2) = TERMINAL POINT OF THE LINE
!     PENX    = PEN NUMBER OR DENSITY (DEPENDING ON PLOTTER)
!     OPT     = -1  TO INITIATE  THE LINE MODE
!             = +1  TO TERMINATE THE LINE MODE
!             =  0  TO DRAW A LINE.
 
 
 REAL, INTENT(IN)                         :: x1
 REAL, INTENT(IN)                         :: y1
 REAL, INTENT(IN)                         :: x2
 REAL, INTENT(IN)                         :: y2
 INTEGER, INTENT(IN)                      :: penx
 INTEGER, INTENT(IN OUT)                  :: opt
 INTEGER :: pen, ploter,tra1,tra2
 REAL :: xy(2,2),infnty
 COMMON /pltdat/ model,ploter,reg(2,2),skpplt(14),skpa(6),npens
 DATA    infnty/ 1.e+10 /
 
 IF (opt /= 0) GO TO 220
 slp = infnty
 b   = 0.
 IF (x1 == x2) GO TO 10
 slp = (y2-y1)/(x2-x1)
 b   =  y1 - slp*x1
 10 xy(1,1) = x1
 xy(2,1) = y1
 xy(1,2) = x2
 xy(2,2) = y2
 
!     CHECK TO SEE IF AN END OF THE LINE IS OUTSIDE THE PLOT REGION.
 
 20 DO  j = 1, 2
   DO  i = 1,2
     IF (xy(i,j) < reg(i,1) .OR. xy(i,j) > reg(i,2)) GO TO 40
   END DO
 END DO
 GO TO 210
 40 DO  i = 1,2
   IF (xy(i,1) < reg(i,1) .AND. xy(i,2) < reg(i,1)) GO TO 230
   IF (xy(i,1) > reg(i,2) .AND. xy(i,2) > reg(i,2)) GO TO 230
 END DO
 
!     AN END IS OUTSIDE THE REGION, BUT NOT THE ENTIRE LINE. FIND THE
!     END POINTS OF THE PORTION OF THE LINE WITHIN THE REGION.
 
 j = 1
 60 i = 1
 70 IF (xy(i,j) >= reg(i,1)) GO TO 130
 ASSIGN 120 TO tra2
 SELECT CASE ( i )
   CASE (    1)
     GO TO 100
   CASE (    2)
     GO TO 110
 END SELECT
 100 ASSIGN 350 TO tra1
 x = reg(1,1)
 GO TO 300
 110 ASSIGN 310 TO tra1
 y = reg(2,1)
 GO TO 300
 120 xy(1,j) = x
 xy(2,j) = y
 
 130 IF (xy(i,j) <= reg(i,2)) GO TO 170
 ASSIGN 160 TO tra2
 SELECT CASE ( i )
   CASE (    1)
     GO TO 140
   CASE (    2)
     GO TO 150
 END SELECT
 140 ASSIGN 350 TO tra1
 x = reg(1,2)
 GO TO 300
 150 ASSIGN 310 TO tra1
 y = reg(2,2)
 GO TO 300
 160 xy(1,j) = x
 xy(2,j) = y
 170 i = i + 1
 IF (i == 2) GO TO 70
 j = j + 1
 IF (j == 2) GO TO 60
 
!     MAKE SURE THE LINE SEGMENT IS WITHIN THE PLOT REGION.
 
 DO  j = 1,2
   DO  i = 1,2
     IF (xy(i,j)+.1 < reg(i,1) .OR. xy(i,j)-.1 > reg(i,2)) GO TO 400
   END DO
 END DO
 
!     FIND THE CORRECT PEN NUMBER FOR THIS PLOTTER.
 
 210 pen = penx
 pen = pen - npens*((pen-1)/npens)
 
!     DRAW THE LINE.
 
 220 CALL line10 (xy(1,1),xy(2,1),xy(1,2),xy(2,2),pen,opt)
 GO TO 400
 
 230 ifl = 0
 DO  j = 1, 2
   DO  m = 1, 2
     IF (ABS(xy(i,j)-reg(i,m)) > 1.0E-8) CYCLE
     ifl = 1
     xy(i,j) = reg(i,m)
   END DO
 END DO
 IF (ifl > 0) THEN
   GO TO    20
 ELSE
   GO TO   400
 END IF
 
 
!     CALCULATE THE EQUATION OF THE LINE TO BE DRAWN.
 
 300 GO TO tra1, (310,350)
 
!     GIVEN Y, CALCULATE X.
 
 310 IF (slp == infnty) GO TO 330
 IF (slp ==     0.) GO TO 320
 x = (y-b)/slp
 GO TO 340
 320 x = infnty
 GO TO 340
 330 x = x1
 340 GO TO tra2, (120,160)
 
!     GIVEN X, CALCULATE Y.
 
 350 IF (slp == infnty) GO TO 370
 IF (slp ==     0.) GO TO 360
 y = slp*x + b
 GO TO 380
 360 y = y1
 GO TO 380
 370 y = infnty
 380 GO TO tra2, (120,160)
 
 400 RETURN
END SUBROUTINE line
