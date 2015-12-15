SUBROUTINE optpx1 (*,stor,nogo,nen,loc1)
     
!     PROCESS PID DATA ON PLIMIT CARD
 
 INTEGER, INTENT(IN OUT)                  :: stor(15)
 INTEGER, INTENT(OUT)                     :: nogo
 INTEGER, INTENT(IN)                      :: nen
 INTEGER, INTENT(IN OUT)                  :: loc1
 INTEGER :: sysbuf,outtap,ycor,thru,nam(2),x(7),iy(1)
 CHARACTER (LEN=23) :: ufm
 COMMON /xmssg / ufm
 COMMON /BLANK / skp1(5),ycor
 COMMON /zzzzzz/ core(1)
 COMMON /system/ sysbuf,outtap
 EQUIVALENCE     (core(1),x(1)),(x(7),iy(1))
 DATA    thru  / 4HTHRU  /
 
 nam(1) = stor(1)
 nam(2) = stor(2)
 IF (stor(6) == thru) GO TO 100
 
!     USER SPECIFIED BY EXPLICIT ID-S
 
 CALL sort (0,0,1,1,stor(5),5)
 
!     CREATE PSEUDO THRU RANGE
!     LOCATE FIRST NONZERO
 
 DO  l = 5,9
   IF (stor(l) /= 0) GO TO 30
 END DO
 CALL page2 (-2)
 WRITE  (outtap,20) ufm,nam
 20 FORMAT (a23,' 2293, NO PID ENTRIES ON PLIMIT CARD (',2A4,2H).)
 nogo = nogo + 1
 GO TO 110
 
!     LOOP ON ENTRIES
 
 30 CONTINUE
 i1 = stor(l)
 i3 = 1
 35 i2 = stor(l+1)
 IF (l-9 < 0) THEN
   GO TO    40
 ELSE IF (l-9 == 0) THEN
   GO TO    60
 ELSE
   GO TO   130
 END IF
 40 IF (i2-i1-i3 < 0) THEN
   GO TO    80
 ELSE IF (i2-i1-i3 == 0) THEN
   GO TO    50
 ELSE
   GO TO    60
 END IF
 
!     THRU CAN BE EXPANDED
 
 50 l  = l  + 1
 i3 = i3 + 1
 GO TO 35
 
!     PUT OUT I1,I2
 
 60 stor(1) = i1
 stor(2) = stor(l)
 IF (loc1+3+nen > ycor) GO TO 120
 CALL bishel (*80,stor,nen,4,iy(loc1))
 70 l = l + 1
 IF (l-9 > 0) THEN
   GO TO   110
 ELSE
   GO TO    30
 END IF
 
!     DUPLICATE ENTRIES FOUND
 
 80 CALL page2 (-2)
 WRITE  (outtap,90) ufm,i1,i2,nam
 90 FORMAT (a23,' 2294, DUPLICATE',i8,' THRU',i8,' RANGE FOR ELEMENT',  &
     1X,2A4,' REJECTED PLIMIT. SCAN CONTINUED.')
 nogo = nogo + 1
 GO TO 70
 
!     USER SPECIFIED BY USING THRU
 
 100 l = 8
 stor(9) = stor(8)
 stor(8) = stor(5)
 GO TO 30
 
!     THIS PLIMIT FINISHED
 
 110 CONTINUE
 RETURN
 
!     INSUFFICIENT CORE
 
 120 CONTINUE
 stor(1) = nam(1)
 stor(2) = nam(2)
 RETURN 1
 
 130 CALL mesage (-7,0,nam)
 GO TO 120
 
END SUBROUTINE optpx1
