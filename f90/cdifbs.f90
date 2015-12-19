SUBROUTINE cdifbs(dz,buf)
     
!     SUBROUTINE TO DO THE FBS PASS TO FIND THE LEFT EIGENVECTOR FOR
!     THE TRANSPOSED MATRIX
 
 
 DOUBLE PRECISION, INTENT(IN OUT)         :: dz(1)
 REAL, INTENT(IN OUT)                     :: buf(1)
 INTEGER :: uprtri    ,eol      ,NAME(2)
 DOUBLE PRECISION :: dtemp     , da
 
!     COMMON   /DESCRP/  LENGTH    ,MAJOR
 COMMON   /zntpkx/  da(2)     ,ii       ,eol
 COMMON   /names/   rd        ,rdrew    ,wrt      ,wrtrew   ,  &
     rew       ,norew    ,eofnrw   ,rsp      , rdp       ,csp      ,cdp
 COMMON   /cdcmpx/  idumm(20) ,iof
 COMMON   /cinvpx/  idum(36)  ,scrfil(11)
 EQUIVALENCE        (scrfil(6),uprtri)  ,(scrfil(8),lowtri) , (idum(2),nrow)
 DATA      NAME  /  4HCDIF,4HBS  /
 
 CALL sswtch (12,idiag)
 
!     BEGIN THE FORWARD PASS USING THE UPPER TRIANGLE
 
 ioff = iof - 1
 CALL gopen (uprtri,buf,rdrew)
 nrow2 = nrow + nrow
 DO  i = 1,nrow
   j = i + i
   CALL intpk (*100,uprtri,0,cdp,0)
   10 CALL zntpki
   IF (ii-i < 0) THEN
     GO TO    30
   ELSE IF (ii-i == 0) THEN
     GO TO    20
   ELSE
     GO TO    40
   END IF
   
!     DIVIDE BY THE DIAGONAL
   
   20 dtemp = (dz(j-1)*da(1)+dz(j  )*da(2))/(da(1)**2 + da(2)**2)
   dz(j) = (dz(j  )*da(1)-dz(j-1)*da(2))/(da(1)**2 + da(2)**2)
   dz(j-1) = dtemp
   GO TO 90
   
!     SUBTRACT OFF NORMAL TERM
   
   30 dz(j-1) = dz(j-1) - dz(2*ii-1)*da(1) + dz(2*ii)*da(2)
   dz(j  ) = dz(j  ) - dz(2*ii-1)*da(2) - dz(2*ii)*da(1)
   GO TO 90
   
!     SUBTRACT OFF ACTIVE COLUMN TERMS
   
   40 k = (i-ioff)*2
   dz(2*ii-1) = dz(2*ii-1) - dz(k-1)*da(1) + dz(k  )*da(2)
   dz(2*ii  ) = dz(2*ii  ) - dz(k  )*da(1) - dz(k-1)*da(2)
   90 IF (eol == 0.0) THEN
     GO TO    10
   ELSE
     GO TO   100
   END IF
   100 CONTINUE
 END DO
 CALL CLOSE (uprtri,rew)
 
!     BEGIN BACKWARD PASS USING THE LOWER TRIANGLE
 
 CALL gopen (lowtri,buf,rdrew)
 CALL skprec (lowtri,nrow)
 DO  i = 1,nrow
   CALL bckrec (lowtri)
   intchn = 0
   CALL intpk (*200,lowtri,0,cdp,0)
   j = (nrow-i+1)*2
   120 CALL zntpki
   IF (ii /= nrow-i+1) GO TO 150
   IF (ii < j/2) GO TO 1010
   
!     PERFORM THE INTERCHANGE
   
   intchn = IFIX(SNGL(da(1)))*2
   IF (idiag /= 0) WRITE (6,131) i,intchn
   131 FORMAT (5H i = ,i5,10HINTCHNG =  ,i5)
   GO TO 190
   130 in1   = j + intchn
   dtemp = dz(j)
   dz(j) = dz(in1)
   dz(in1) = dtemp
   dtemp   = dz(j-1)
   dz(j-1  ) = dz(in1-1)
   dz(in1-1) = dtemp
   GO TO 200
   150 CONTINUE
   dz(j-1) = dz(j-1) - dz(2*ii-1)*da(1) + dz(2*ii)*da(2)
   dz(j  ) = dz(j  ) - dz(2*ii-1)*da(2) - dz(2*ii)*da(1)
   190 IF (eol == 0.0) THEN
     GO TO   120
   END IF
   195 IF (intchn > 0) THEN
     GO TO   130
   END IF
   200 CALL bckrec (lowtri)
 END DO
 CALL CLOSE (lowtri,rew)
 RETURN
 
 1010 CALL mesage (-7,lowtri,NAME)
 RETURN
END SUBROUTINE cdifbs
