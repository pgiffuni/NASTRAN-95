SUBROUTINE intfbs (dx,dy,iobuf)
     
!     GIVEN THE TRIANGULAR FACTORS FOR A GENERAL MATRIX, INTFBS WILL
!     PERFORM THE FORWARD-BACKWARD SUBSTITUTION NECESSARY TO SOLVE
!     A SYSTEM OF EQUATIONS
 
!     DEFINITION OF INPUT PARAMETERS
 
!     FILEL    =  MATRIX CONTROL BLOCK FOR THE LOWER TRIANGLE L
!     FILEU    =  MATRIX CONTROL BLOCK FOR THE UPPER TRIANGLE U
!     DX       =  THE LOAD VECTOR B
!     DY       =  THE SOLUTION VECTOR X
!     IOBUF    =  THE INPUT BUFFER
 
!     NAMED COMMONS
 
 
 REAL, INTENT(IN)                         :: dx(1)
 REAL, INTENT(OUT)                        :: dy(1)
 INTEGER, INTENT(IN OUT)                  :: iobuf(1)
 INTEGER :: filel     ,fileu    ,typear   ,rdp      , parm(4)   ,rsp      ,eol
 
 COMMON   /names /  rd        ,rdrew    ,wrt      ,wrtrew   ,  &
     rew       ,norew    ,eofnrw   ,rsp      ,  &
     rdp       ,csp      ,cdp      ,sqr      ,  &
     rect      ,diag     ,lowtri   ,uprtri   , sym       ,row      ,identy
 COMMON   /zntpkx/  a(4)      ,ii       ,eol
 COMMON   /infbsx/  filel(7)  ,fileu(7)
 COMMON   /trdxx /  idummy(27),iopen
 EQUIVALENCE        (a(1),da) , (filel(3) ,nrow)
 DATA      parm(3), parm(4)   /4HINTF   ,4HBS   /
 
 
!     TRANSFER THE LOAD VECTOR TO THE SOLUTION VECTOR
 
 DO  i = 1,nrow
   dy(i)  = dx(i)
 END DO
 typear = rsp
 
!     OPEN FILE FOR THE LOWER TRIANGLE
 
 parm(2) = filel(1)
 IF (iopen == -10) GO TO 15
 IF (iopen ==   0) CALL OPEN (*200,filel(1),iobuf,rdrew)
 CALL fwdrec (*210,filel(1))
 
!     BEGIN FORWARD PASS
 
 15 j = 1
 20 CALL intpk (*90,filel(1),0,typear,0)
 30 IF (eol == 0.0) THEN
   GO TO    40
 ELSE
   GO TO   220
 END IF
 40 CALL zntpki
 IF (j-ii < 0) THEN
   GO TO    80
 ELSE IF (j-ii == 0) THEN
   GO TO    50
 ELSE
   GO TO    30
 END IF
 
!     PERFORM THE REQUIRED ROW INTERCHANGE
 
 50 in1   = j + IFIX(a(1))
 dtemp = dy(j)
 dy(j) = dy(in1)
 dy(in1) = dtemp
 60 IF (eol == 0.0) THEN
   GO TO    70
 ELSE
   GO TO    90
 END IF
 70 CALL zntpki
 80 dy(ii) = dy(ii) - dy(j)*da
 GO TO 60
 90 j = j + 1
 IF (j < nrow) GO TO 20
 CALL REWIND (filel(1))
 IF (iopen ==   0) CALL CLOSE (filel(1),rew)
 IF (iopen == -10) CALL skprec (filel,1)
 
!     BEGIN BACKWARD PASS
 
 ioff = fileu(7) - 1
 parm(2) = fileu(1)
 IF (iopen == -10) GO TO 95
 IF (iopen ==   0) CALL OPEN (*200,fileu(1),iobuf,rdrew)
 CALL fwdrec (*210,fileu(1))
 95 j = nrow
 100 CALL intpk (*220,fileu(1),0,typear,0)
 IF (eol == 0.0) THEN
   GO TO   120
 ELSE
   GO TO   220
 END IF
 120 CALL zntpki
 i = nrow - ii + 1
 IF (i /= j) GO TO 170
 
!     DIVIDE BY THE DIAGONAL
 
 dy(i) = dy(i)/da
 
!     SUBTRACT OFF REMAINING TERMS
 
 140 IF (i > j) GO TO 120
 IF (eol == 0.0) THEN
   GO TO   160
 ELSE
   GO TO   190
 END IF
 160 CALL zntpki
 i = nrow - ii + 1
 170 in1 = i
 in2 = j
 IF (i < j) GO TO 180
 k   = in1
 in1 = in2 - ioff
 in2 = k
 180 dy(in1) = dy(in1) - dy(in2)*da
 GO TO 140
 190 j = j - 1
 IF (j > 0) GO TO 100
 CALL  REWIND (fileu(1))
 IF (iopen ==   0) CALL CLOSE (fileu(1),rew)
 IF (iopen == -10) CALL skprec (fileu,1)
 RETURN
 
 200 parm(1) = -1
 GO TO 230
 210 parm(1) = -2
 GO TO 230
 220 parm(1) = -5
 230 CALL mesage (parm(1),parm(2),parm(3))
 RETURN
END SUBROUTINE intfbs
