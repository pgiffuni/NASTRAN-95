SUBROUTINE cf2fbs (tpose,xout,iobuf)
!*******
!     CF2FBS PERFORMS THE DOUBLE-PRECISION FORWARD AND BACKWARD SWEEPS
!     FOR THE COMPLEX FEER METHOD. THESE SWEEPS CONSTITUTE THE
!     OPERATIONAL INVERSE (MATRIX INVERSION).
!*******
!     DEFINITION OF INPUT AND OUTPUT PARAMETERS
!*******
!     TPOSE    = .FALSE. --- PERFORM OPERATION L * U
!              = .TRUE.  --- PERFORM OPERATION U-TRANSPOSE * L-TRANSPOSE
!     XOUT     = INPUT VECTOR GETS TRANSFORMED TO OUTPUT VECTOR
!     IOBUF    = INPUT  GINO BUFFER
!*******
 
 LOGICAL, INTENT(IN OUT)                  :: tpose(1)
 DOUBLE PRECISION, INTENT(IN OUT)         :: xout(1)
 INTEGER, INTENT(IN OUT)                  :: iobuf(1)
 DOUBLE PRECISION :: dtemp    , da       ,unidum
 INTEGER :: NAME(2)  , eol      ,cdp
 LOGICAL :: symmet   ,qpr
 COMMON  /feeraa/  aadum(117),mcblt(7),mcbut(7)
 COMMON  /feerxc/  xcd01(4) ,symmet   ,xcd02(9) ,nswp ,xcd03(6) ,qpr
 COMMON  /zntpkx/  da(2)    ,ii       ,eol
 COMMON  /names /  rd       ,rdrew    ,wrt      ,wrtrew  &
     ,rew      ,norew    ,eofnrw   ,rsp ,rdp      ,csp      ,cdp
 COMMON  /system/  ksystm   ,nout
 EQUIVALENCE       (aadum(42),iscr6)
 DATA   NAME       /4HCF2F,4HBS  /
 
 IF (qpr) WRITE (nout,8887) tpose(1),symmet,nswp,iscr6
 8887 FORMAT(1H0,12HENTER cf2fbs,8X,11HTRANSPOSE =,l2,l9,2I10)
 junk = 0
 IF (tpose(1) .AND. .NOT.symmet) GO TO 399
!*******
!     BELOW FOR OPERATION L * U
!     (LOGIC COPIED FROM SUBROUTINE CINFBS)
!*******
!     BEGIN FORWARD PASS USING THE LOWER TRIANGLE
!*******
 CALL gopen (mcblt(1),iobuf(1),rdrew)
 j = 1
 100 CALL intpk(*200,mcblt(1),0,cdp,0)
 110 IF (eol == 0.0) THEN
   GO TO   120
 ELSE
   GO TO  3010
 END IF
 120 CALL zntpki
 IF (qpr) WRITE (nout,8882) da,ii,eol,j
 8882          FORMAT(1H ,4HDA =,2D16.8,4X,4HII =,i6,  &
     4X,5HEOL =,i2,4X,3HJ =,i4)
 IF (j-ii < 0) THEN
   GO TO   184
 ELSE IF (j-ii == 0) THEN
   GO TO   130
 ELSE
   GO TO   110
 END IF
!*******
!     PERFORM THE REQUIRED ROW INTERCHANGE
!*******
 130 in1 = ( j + IFIX(SNGL(da(1))) )*2 - 1
 IF (qpr) WRITE (nout,8883) in1,eol
 8883          FORMAT(1H ,3X,5HIN1 =,i6,4X,5HEOL =,i2)
 in2 = in1+1
 j2 = 2*j
 unidum = xout(j2)
 xout(j2) = xout(in2)
 xout(in2) = unidum
 j2 = j2-1
 unidum = xout(j2)
 xout(j2) = xout(in1)
 xout(in1) = unidum
 160 IF (eol == 0.0) THEN
   GO TO   170
 ELSE
   GO TO   200
 END IF
 170 CALL zntpki
 IF (qpr) WRITE (nout,8882) da,ii,eol,j
 184 ii2 = 2*ii
 ii1 = ii2-1
 j2 = 2*j
 j1 = j2-1
 xout(ii1) = xout(ii1) - da(1)*xout(j1) + da(2)*xout(j2)
 xout(ii2) = xout(ii2) - da(2)*xout(j1) - da(1)*xout(j2)
 GO TO 160
 200 j = j+1
 IF (j < nswp) GO TO 100
 CALL CLOSE (mcblt(1),rew)
!*******
!     BEGIN BACKWARD PASS USING THE UPPER TRIANGLE
!*******
 ioff = mcbut(7)-1
 IF (qpr) WRITE (nout,8866) ioff,mcblt,mcbut
 8866          FORMAT(1H ,15(1X,i7))
 CALL gopen (mcbut(1),iobuf(1),rdrew)
 j = nswp
 210  CALL intpk(*3020,mcbut(1),0,cdp,0)
 IF (eol == 0.0) THEN
   GO TO   230
 ELSE
   GO TO  3020
 END IF
 230 CALL zntpki
 IF (qpr) WRITE (nout,8882) da,ii,eol,j
 i = nswp - ii + 1
 IF (i /= j) GO TO 275
!*******
!     DIVIDE BY THE DIAGONAL
!*******
 i2 = 2*i
 i1 = i2-1
 unidum = 1.d0/(da(1)**2+da(2)**2)
 dtemp = (da(1)*xout(i1)+da(2)*xout(i2))*unidum
 xout(i2) = (da(1)*xout(i2)-da(2)*xout(i1))*unidum
 xout(i1) = dtemp
 IF (qpr) WRITE (nout,8884)
 8884          FORMAT(1H ,6X,8HDIAGONAL)
!*******
!     SUBTRACT OFF REMAINING TERMS
!*******
 255 IF (i > j) GO TO 230
 IF (eol == 0.0) THEN
   GO TO   270
 ELSE
   GO TO   300
 END IF
 270 CALL zntpki
 IF (qpr) WRITE (nout,8882) da,ii,eol,j
 i = nswp - ii + 1
 275 in1 = i
 in2 = j
 IF (i < j) GO TO 279
 k = in1
 in1 = in2-ioff
 junk = 1
 IF (in1 <= 0) GO TO 3020
 in2 = k
 279 in1 = 2*in1
 in2 = 2*in2
 ii1 = in1-1
 ii2 = in2-1
 IF (qpr) WRITE (nout,8820) i,j,ii1,ii2
 8820       FORMAT(1H ,3HI =,i6,6X,3HJ =,i6,6X,5HII1 =,i6,6X,5HII2 =,i6)
 xout(ii1) = xout(ii1) - da(1)*xout(ii2) + da(2)*xout(in2)
 xout(in1) = xout(in1) - da(2)*xout(ii2) - da(1)*xout(in2)
 GO TO 255
 300 j = j-1
 IF (j > 0) GO TO 210
 CALL CLOSE (mcbut(1),rew)
 GO TO 4000
!*******
!     BELOW FOR OPERATION U-TRANSPOSE * L-TRANSPOSE
!     (LOGIC COPIED FROM SUBROUTINE CDIFBS)
!*******
!     BEGIN THE FORWARD PASS USING THE UPPER TRIANGLE
!*******
 399 ioff = mcbut(7)-1
 IF (qpr) WRITE (nout,2216) ioff
 2216          FORMAT(1H ,30X,6HIOFF =,i10)
 mcsave = mcbut(1)
 mcbut(1) = iscr6
 CALL gopen (mcbut(1),iobuf(1),rdrew)
 DO   i = 1,nswp
   IF (qpr) WRITE (nout,2218) i
   2218          FORMAT(1H ,12HLOOP INDEX =,i6)
   j = i+i
   CALL intpk(*500,mcbut(1),0,cdp,0)
   410 CALL zntpki
   IF (qpr) WRITE (nout,2224) ii,eol,da
   2224          FORMAT(1H ,4HII =,i14,6X,5HEOL =,i2, 8X,4HDA =,2D16.8)
   IF (ii-i < 0) THEN
     GO TO   430
   ELSE IF (ii-i == 0) THEN
     GO TO   420
   ELSE
     GO TO   440
   END IF
!*******
!     DIVIDE BY THE DIAGONAL
!*******
   420 i1 = j-1
   unidum = 1.d0/(da(1)**2+da(2)**2)
   dtemp = (xout(i1)*da(1) + xout(j)*da(2))*unidum
   xout(j) = (xout(j)*da(1) - xout(i1)*da(2))*unidum
   xout(i1) = dtemp
   IF (qpr) WRITE (nout,8884)
   GO TO 490
!*******
!     SUBTRACT OFF NORMAL TERM
!*******
   430 i2 = ii+ii
   i1 = i2-1
   j1 = j-1
   xout(j1) = xout(j1) - xout(i1)*da(1) + xout(i2)*da(2)
   xout(j) = xout(j) - xout(i1)*da(2) - xout(i2)*da(1)
   GO TO 490
!*******
!     SUBTRACT OFF ACTIVE COLUMN TERMS
!*******
   440 k = (i-ioff)*2
   junk = 1
   in1 = k
   IF (in1 <= 0) GO TO 3020
   i2 = ii+ii
   i1 = i2-1
   j1 = k-1
   xout(i1) = xout(i1) - xout(j1)*da(1) + xout(k)*da(2)
   xout(i2) = xout(i2) - xout(k)*da(1) - xout(j1)*da(2)
   490 IF (eol == 0.0) THEN
     GO TO   410
   ELSE
     GO TO   500
   END IF
 END DO
 CALL CLOSE (mcbut(1),rew)
 mcbut(1) = mcsave
!*******
!     BEGIN BACKWARD PASS USING THE LOWER TRIANGLE
!*******
 CALL gopen (mcblt(1),iobuf(1),rdrew)
 CALL skprec (mcblt(1),nswp)
 DO   i = 1,nswp
   IF (qpr) WRITE (nout,2218) i
   CALL bckrec (mcblt(1))
   intchn = 0
   CALL intpk(*600,mcblt(1),0,cdp,0)
   j = (nswp-i+1)*2
   520 CALL zntpki
   IF (qpr) WRITE (nout,2224) ii,eol,da
   IF (ii /= nswp-i+1) GO TO 550
   IF (ii < j/2) GO TO 3010
!*******
!     PERFORM THE INTERCHANGE
!*******
   intchn = IFIX(SNGL(da(1)))*2
   IF (qpr) WRITE (nout,2226) intchn
   2226          FORMAT(1H ,4X,11HINTERCHANGE,i6)
   GO TO 590
   530 in1 = j+intchn
   IF (qpr) WRITE (nout,2232) j,intchn,in1
   2232          FORMAT(1H ,15X,3I6)
   dtemp = xout(j)
   xout(j) = xout(in1)
   xout(in1) = dtemp
   j1 = j-1
   i1 = in1-1
   dtemp = xout(j1)
   xout(j1) = xout(i1)
   xout(i1) = dtemp
   GO TO 600
   550 j1 = j-1
   i2 = ii+ii
   i1 = i2-1
   xout(j1) = xout(j1) - xout(i1)*da(1) + xout(i2)*da(2)
   xout(j) = xout(j) - xout(i1)*da(2) - xout(i2)*da(1)
   590 IF (eol == 0.0) THEN
     GO TO   520
   END IF
   595 IF (intchn > 0) THEN
     GO TO   530
   END IF
   600 CALL bckrec (mcblt(1))
 END DO
 CALL CLOSE (mcblt(1),rew)
 GO TO 4000
 3010 j = mcblt(1)
 GO TO 3040
 3020 j = mcbut(1)
 3040 CALL mesage (-5,j,NAME)
 4000 CONTINUE
 IF (qpr.AND.junk == 0) WRITE (nout,5516)
 5516          FORMAT(1H0,30X,13HIOFF NOT used,/1H )
 IF (qpr.AND.junk /= 0) WRITE (nout,5518)
 5518          FORMAT(1H0,30X,13HIOFF was used,/1H )
 RETURN
END SUBROUTINE cf2fbs
