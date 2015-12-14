SUBROUTINE invfbs (dx,dy,iobuf)
     
!     DOUBLE PRECISION VERSION
 
!     INVFBS IS A SPECIAL FORWARD-BACKWARD SUBSTITUTION ROUTINE FOR
!     INVPWR. IT OPERATES ON CONJUNCTION WITH SDCOMP.
!     THE ARITHMETIC PRECISION IS THAT OF THE INPUT FILE
 
!     FILEL    =  MATRIX CONTROL BLOCK FOR THE LOWER TRIANGLE L
!     FILEU    =  MATRIX CONTROL BLOCK FOR THE UPPER TRIANGLE U
!     DX       =  THE LOAD VECTOR B
!     DY       =  THE SOLUTION VECTOR X
!     IOBUF    =  THE INPUT BUFFER
 
!     COMMENT FROM G.CHAN/UNISYS, 6/89
!     IF LOAD IS SUDDENLY INCREADED TO A LARGE VALUE, THE VAX MACHINE
!     MAY BLOW ITS TOP (ARITHMETIC FAULT, FLOATING OVERFLOW) BECAUSE
!     VAX DOUBLE PRECISION REAL NUMBERS ARE LIMITED TO 10**38, SAME
!     LIMIT AS THE SINGLE PRECISION REAL NUMBERS. OTHER MACHINES ALLOW
!     MUCH LARGER LIMITS FOR DOUBLE PRECISION NUMBERS.
 
 
 DOUBLE PRECISION, INTENT(IN)             :: dx(1)
 DOUBLE PRECISION, INTENT(OUT)            :: dy(1)
 INTEGER, INTENT(IN OUT)                  :: iobuf(1)
 INTEGER :: filel     ,fileu    ,typear   ,rdp      ,  &
     parm(4)   ,eol      ,ijj(2)
 DOUBLE PRECISION :: da       ,dtemp    , djj       ,dyj      ,epsi
 
 CHARACTER (LEN=25) :: sfm
 CHARACTER (LEN=29) :: uim
 CHARACTER (LEN=25) :: uwm
 CHARACTER (LEN=23) :: ufm
 COMMON   /xmssg /  ufm       ,uwm      ,uim      ,sfm
 COMMON   /machin/  mach
 COMMON   /system/  ibuf      ,nout
!     COMMON   /DESCRP/  LENGTH    ,MAJOR
 COMMON   /names /  rd        ,rdrew    ,wrt      ,wrtrew   ,  &
     rew       ,norew    ,eofnrw   ,rsp      ,  &
     rdp       ,csp      ,cdp      ,sqr      ,  &
     rect      ,diag     ,lowtri   ,uprtri   , sym       ,row      ,identy
 COMMON   /TYPE  /  rc(2)     ,nwds(4)
 COMMON   /zntpkx/  a(4)      ,ii       ,eol
 COMMON   /infbsx/  filel(7)  ,fileu(7)
 COMMON   /trdxx /  idummy(27),iopen
 EQUIVALENCE        (a(1),da) ,(filel(3),nrow)    ,(djj,ijj(1))
 DATA      epsi  /  1.0D-24   /
 DATA      parm(3), parm(4)   /4HINVF,4HBS  /
 
 
!     TRANSFER THE LOAD VECTOR TO THE SOLUTION VECTOR
 
 DO  i = 1,nrow
   dy(i)  = dx(i)
 END DO
 typear = rdp
 
!     OPEN FILE FOR THE LOWER TRIANGLE
!     IOPEN WAS SET TO -20 BY STEP2
 
 parm(2) = filel(1)
 IF (iopen ==  -20) CALL fwdrec (*360,filel(1))
 IF (filel(7) < 0) GO TO 300
 
!     NASTRAN ORIGINAL CODE
 
!     BEGIN FORWARD PASS
 
 j = 1
 20 CALL intpk (*100,filel(1),0,typear,0)
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
 
 50 in1   = j+IFIX(SNGL(da))
 dtemp = dy(j)
 dy(j) = dy(in1)
 dy(in1) = dtemp
 60 IF (eol == 0.0) THEN
   GO TO    70
 ELSE
   GO TO   100
 END IF
 70 CALL zntpki
 80 IF (mach /= 5 .OR.  &
     (DABS(da) < 1.d+19 .AND. DABS(dy(j)) < 1.d+19)) GO TO 90
 x1 = ALOG10(ABS(SNGL(da)))
 x2 = ALOG10(ABS(SNGL(dy(j))))
 IF (x1+x2 > 38.) GO TO 200
 90 dy(ii) = dy(ii) - dy(j)*da
 GO TO 60
 100 j = j + 1
 IF (j < nrow) GO TO 20
 CALL REWIND (filel(1))
 IF (iopen == -20) GO TO 110
 CALL skprec (filel,1)
 
!     BEGIN BACKWARD PASS
 
 110 ioff    = fileu(7) - 1
 parm(2) = fileu(1)
 IF (iopen == -20) CALL fwdrec (*360,fileu(1))
 j = nrow
 120 CALL intpk (*220,fileu(1),0,typear,0)
 IF (eol /= 0) GO TO 220
 130 CALL zntpki
 i = nrow - ii + 1
 IF (i /= j) GO TO 150
 
!     DIVIDE BY THE DIAGONAL
 
 dy(i) = dy(i)/da
 
!     SUBTRACT OFF REMAINING TERMS
 
 140 IF (i   > j) GO TO 130
 IF (eol /= 0) GO TO 180
 CALL zntpki
 i   = nrow - ii + 1
 150 in1 = i
 in2 = j
 IF (i < j) GO TO 160
 k   = in1
 in1 = in2 - ioff
 in2 = k
 160 IF (mach /= 5 .OR.  &
     (DABS(da) < 1.d+19 .AND. DABS(dy(in2)) < 1.d+19)) GO TO 170
 x1 = ALOG10(ABS(SNGL(da)))
 x2 = ALOG10(ABS(SNGL(dy(in2))))
 IF (x1+x2 > 38.) GO TO 200
 170 dy(in1) = dy(in1) - dy(in2)*da
 GO TO 140
 180 j = j - 1
 IF (j > 0) GO TO 120
 CALL REWIND (fileu)
 IF (iopen == -20) RETURN
 CALL skprec (fileu,1)
 GO TO 450
 
 200 WRITE  (nout,210) sfm,parm(1),parm(2)
 210 FORMAT (a25,' FROM ',2A4,'- SOLUTION VECTOR VALUE OVERFLOWS,',/5X,  &
     'POSSIBLY DUE TO SUDDEN INCREASE OF LARGE LOAD VECTOR OR ',  &
     'OTHER INPUT CONDITION')
 GO TO 420
 220 parm(1) = -5
 GO TO 440
 
 
!     NEW METHOD
!     FILEL HAS BEEN RE-WRITTEN FORWARD FIRST THAN BACKWARD BY UNPSCR
!     IN INVP3)
 
!     THE LOAD VECTOR DX WILL BE DESTROYED IN THIS NEW METHOD
 
!     FORWARD SWEEP DIRECTLY ON SOLUTION VECTOR DY
 
 300 ifile  =-filel(7)
 parm(2)= ifile
 nwd    = nwds(filel(5))
 IF (filel(4) /= 2) GO TO 400
 ifw = +1
 CALL REWIND (ifile)
 CALL skprec (ifile,1)
 CALL READ (*360,*370,ifile,dx,2,0,i)
 ntms = 0
 DO  j = 1,nrow
   djj = dx(ntms+1)
   ii  = ijj(1)
   jj  = ijj(2)
   IF (ii /= j) GO TO 380
   ntms = jj - ii + 1
   ji = ntms*nwd + 2
   CALL READ (*360,*370,ifile,dx,ji,0,i)
   IF (ntms <= 1) GO TO 320
   dyj = dy(j)
   IF (DABS(dyj) < epsi) GO TO 320
   DO  i = 2,ntms
     ii = ii + 1
     dy(ii) = dy(ii) + dx(i)*dyj
   END DO
   320 dy(j) = dy(j)/dx(1)
 END DO
 
!     BACKWARD SUBSTITUTION OMIT DIAGONAL
 
 ifw = -1
 IF (nrow == 1) GO TO 450
 j  = nrow
 DO  jx = 1,nrow
   djj = dx(ntms+1)
   ii  = ijj(1)
   jj  = ijj(2)
   IF (ii /= j) GO TO 380
   ntms = jj - ii + 1
   ji = ntms*nwd + 2
   CALL READ (*360,*370,ifile,dx,ji,0,i)
   IF (ntms <= 1) GO TO 340
   DO  i = 2,ntms
     ii = ii + 1
     dy(j) = dy(j) + dx(i)*dy(ii)
   END DO
   340 j = j - 1
 END DO
 GO TO 450
 
!     ERROR
 
 360 parm(1) = -2
 GO TO 440
 370 parm(1) = -3
 GO TO 440
 380 WRITE  (nout,390) ifw,ii,j
 390 FORMAT ('   ERROR IN INVFBS.   IFW),II,J =',i3,1H),2I6)
 GO TO 420
 400 WRITE  (nout,410) filel(4)
 410 FORMAT ('0*** FILEL MATRIX IN WRONG FORMAT. UNPSCR FLAG =',i3)
 420 parm(1) = -37
 440 CALL mesage (parm(1),parm(2),parm(3))
 
 450 RETURN
END SUBROUTINE invfbs
