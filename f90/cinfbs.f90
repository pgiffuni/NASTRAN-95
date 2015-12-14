SUBROUTINE cinfbs (dx,dy,iobuf)
     
!     CINVFB DOES THE FORWARD AND BACKWARD PASS FOR COMPLEX INVERSE POWE
 
 
 DOUBLE PRECISION, INTENT(IN)             :: dx(1)
 DOUBLE PRECISION, INTENT(OUT)            :: dy(1)
 INTEGER, INTENT(IN OUT)                  :: iobuf(1)
 INTEGER :: NAME(2)   ,typear   ,cdp      , eol
 DOUBLE PRECISION :: da       ,dtemp
!     COMMON   /DESCRP/  LENGTH    ,MAJOR
 COMMON   /names /  rd        ,rdrew    ,wrt      ,wrtrew   ,  &
     rew       ,norew    ,eofnrw   ,rsp      , rdp       ,csp      ,cdp
 COMMON   /zntpkx/  da(2)     ,ii       ,eol
 COMMON   /cinfbx/  ifill(7)  ,ifilu(7)
 EQUIVALENCE        (ifill(3),nrow)
 DATA      NAME  /  4HCINF, 4HBS   /
 
!     TRANSFER THE LOAD VECTOR TO THE SOLUTION VECTOR
 
 typear = cdp
 nrow2  = nrow + nrow
 DO  i = 1,nrow2
   dy(i) = dx(i)
 END DO
 
!     BEGIN FORWARD PASS
 
!     CALL GOPEN (IFILL(1),IOBUF,RDREW)
 j = 1
 100 CALL intpk (*200,ifill(1),0,typear,0)
 110 IF (eol == 0.0) THEN
   GO TO   120
 ELSE
   GO TO  3010
 END IF
 120 CALL zntpki
 IF (j-ii < 0) THEN
   GO TO   184
 ELSE IF (j-ii == 0) THEN
   GO TO   130
 ELSE
   GO TO   110
 END IF
 
!     PERFORM THE REQUIRED ROW INTERCHANGE
 
 130 in1 = (j+IFIX(SNGL(da(1))))*2 - 1
 dtemp     = dy(2*j-1)
 dy(2*j-1) = dy(in1)
 dy(in1)   = dtemp
 dtemp     = dy(2*j)
 dy(2*j)   = dy(in1+1)
 dy(in1+1) = dtemp
 160 IF (eol == 0.0) THEN
   GO TO   170
 ELSE
   GO TO   200
 END IF
 170 CALL zntpki
 184 dy(2*ii-1) = dy(2*ii-1) - dy(2*j-1)*da(1) + dy(2*j)*da(2)
 dy(2*ii  ) = dy(2*ii  ) - dy(2*j-1)*da(2) - dy(2*j)*da(1)
 GO TO 160
 200 j = j + 1
 IF (j < nrow) GO TO 100
 CALL REWIND (ifill(1))
 
!     BEGIN BACKWARD PASS
 
 ioff = ifilu(7) - 1
 j = nrow
 210 CALL intpk (*3020,ifilu(1),0,typear,0)
 IF (eol == 0.0) THEN
   GO TO   230
 ELSE
   GO TO  3020
 END IF
 230 CALL zntpki
 i = nrow - ii + 1
 IF (i /= j) GO TO 275
 
!     DIVIDE BY THE DIAGONAL
 
 dtemp     = (dy(2*i-1)*da(1)+dy(2*i)*da(2))/(da(1)**2+da(2)**2)
 dy(2*i  ) = (dy(2*i)*da(1)-dy(2*i-1)*da(2))/(da(1)**2+da(2)**2)
 dy(2*i-1) = dtemp
 
!     SUBTRACT OFF REMAINING TERMS
 
 255 IF (i > j) GO TO 230
 IF (eol == 0.0) THEN
   GO TO   270
 ELSE
   GO TO   300
 END IF
 270 CALL zntpki
 i = nrow - ii + 1
 275 in1 = i
 in2 = j
 IF (i < j) GO TO 279
 k   = in1
 in1 = in2 - ioff
 in2 = k
 279 in1 = in1 + in1 - 1
 in2 = in2 + in2 - 1
 dy(in1  ) = dy(in1  ) - dy(in2)*da(1) + dy(in2+1)*da(2)
 dy(in1+1) = dy(in1+1) - dy(in2)*da(2) - dy(in2+1)*da(1)
 GO TO 255
 300 j = j - 1
 IF (j > 0) GO TO 210
 CALL REWIND (ifilu(1))
 RETURN
 
 3010 ifile = ifill(1)
 GO TO 3040
 3020 ifile = ifilu(1)
 3040 CALL mesage (-5,ifile,NAME)
 RETURN
END SUBROUTINE cinfbs
