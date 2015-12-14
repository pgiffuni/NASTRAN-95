SUBROUTINE idplot (idx)
     
 
 INTEGER, INTENT(OUT)                     :: idx
 COMMON /output/ skpout(32,6),id(32)
 COMMON /pltdat/ skpplt(2),xymin(2),xymax(2),axymax(2),edge(12)  &
     ,      skpa(3),cntx,cnty,skpb(4),pltype
 INTEGER :: pltype
 
 INTEGER :: BLANK
 REAL :: SAVE(2,4)
 DATA BLANK,linsiz / 1H ,3 /
 
!     DOES A PLOT ID EXIST AT ALL
 
 idx = 1
 DO  i = 1,20
   IF (id(i) /= BLANK)  GO TO 102
 END DO
 idx = 0
 GO TO 200
 
 102 DO  i = 1,2
   SAVE(i,1) = xymin(i)
   xymin(i) = 0.
   SAVE(i,2) = xymax(i)
   xymax(i) = axymax(i)+edge(i)
   SAVE(i,3) = axymax(i)
   axymax(i) = xymax(i)
   SAVE(i,4) = edge(i)
   edge(i) = 0.
 END DO
 nlines = (axymax(2)-7.*cnty) / FLOAT(2*linsiz) + .1
 IF (IABS(pltype) /= 1)  GO TO 122
 
!     FILL TOP HALF OF PLOT WITH X-AXIS LINES ALL THE WAY ACROSS.
 
 CALL axis (0,0,0,0,0,-1)
 DO  i = 1,nlines
   y = xymax(2) - FLOAT((i-1)*linsiz)
   CALL axis (xymin(1),y,xymax(1),y,1,0)
 END DO
 
!     PRINT THE PLOT ID 2 TIMES IN THE MIDDLE OF THE PLOT.
 
 CALL PRINT (0,0,0,0,0,-1)
 x = xymin(1) + AMAX1(0.,(axymax(1)-80.*cntx)/2.)
 yy = y-cnty
 DO  i = 1,2
   y = yy - cnty*FLOAT(i-1)
   CALL PRINT (x,y,1,id,20,0)
 END DO
 
!     FILL BOTTOM HALF OF PLOT WITH X-AXIS LINES ALL THE WAY ACROSS.
 
 CALL axis (0,0,0,0,0,-1)
 DO  i = 1,nlines
   y = xymin(2) + FLOAT((i-1)*linsiz)
   CALL axis (xymin(1),y,xymax(1),y,1,0)
 END DO
 CALL axis (0,0,0,0,0,1)
 GO TO 125
 
!     NOT A CRT PLOTTER. TYPE THE ID ONCE AT THE BOTTOM OF THE PAPER.
 
 122 CALL PRINT (0,0,0,0,0,-1)
 x = xymin(1) + AMAX1(0.,(axymax(1)-80.*cntx)/2.)
 y = 0.
 IF (pltype < 0)  y=cnty/2.
 CALL PRINT (x,y,1,id,20,0)
 
!     END OF ID PLOT. PUT BLANKS IN THE PLOT ID.
 
 125 CALL PRINT (0,0,0,0,0,1)
 DO  i = 1,20
   id(i) = BLANK
 END DO
 DO  i = 1,2
   xymin(i) = SAVE(i,1)
   xymax(i) = SAVE(i,2)
   axymax(i) = SAVE(i,3)
   edge(i) = SAVE(i,4)
 END DO
 
 200 RETURN
END SUBROUTINE idplot
