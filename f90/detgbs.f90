SUBROUTINE detgbs (iy,iobuf,kcount)
     
!     DETGFBS IS A SPECIAL VERSION OF THE DETFBS ROUTINE AND IS USED BY
!     THE REAL DETERMINANT METHOD FOR UNSYMMETRIC DECOMPOSITION.
!     IT IS SUITABLE FOR BOTH SINGLE AND DOUBLE PRECISION OPERATION.
 
 
!     DEFINITION OF PARAMETERS
!     ------------------------
 
!     FILEU  = MATRIX CONTROL BLOCK FOR THE UPPER TRIANGLE
!     FILEV  = SAME AS FILEU
!     FILEVT = MATRIX CONTROL BLOCK FOR THE TRANSPOSE OF THE UPPER
!              TRIANGLE
!     X, DX  = THE SOLUTION VECTOR
!     Y, DY  = REGION USED FOR UNPACKING
!     IY     = POINTER TO Y (DY) RELATIVE TO X (DX)
!     IOBUF  = THE INPUT BUFFER
!     NROW   = MATRIX SIZE
!     KCOUNT = EIGENVALUE COUNTER
 
 
 INTEGER, INTENT(IN)                      :: iy
 INTEGER, INTENT(IN OUT)                  :: iobuf(7)
 INTEGER, INTENT(IN)                      :: kcount
 INTEGER :: fileu(7),parm(4) , option , fc  &
     ,                scr3    ,scr4    ,scr6    ,eol
 REAL :: x(1)    ,y(1)
 DOUBLE PRECISION :: dx(1)   ,dy(1)   ,dxmin   ,dsdiag ,da(2)
 COMMON /zzzzzz/  core(1)
 COMMON /detmx /  dum3(36),ipdeta
 COMMON /names /  rd      ,rdrew   ,wrt     ,wrtrew  ,rew
 COMMON /regean/  dum1(23),scr3    ,scr4    ,dum2(11),scr6
 COMMON /reigkr/  option
 COMMON /unpakx/  itypex  , iunpak ,junpak  ,incr
 COMMON /zntpkx/  a(4)    , ii     ,eol
 COMMON /dcompx/  fa(7)   , fl(7)  ,fc(7)
 EQUIVALENCE      (core(1),x(1),dx(1),y(1),dy(1))  &
     ,                (xmin,dxmin) ,   (sdiag,dsdiag)  &
     ,                (a(1),da(1) )
 DATA    parm(3), parm(4) / 4HDETG, 4HFBS  /
 
 fileu(1) = fc(1)
 CALL rdtrl ( fileu )
 itypex = fileu(5)
 nrow   = fileu(2)
 ioff   = fileu(7)-1
 iprec  = 1
 IF ( itypex == 2 ) iprec = 2
 INDEX  = -1
 incr   = 1
 nfile  = fileu(1)
 INDEX  = 1
 lcore  = ipdeta - iy*itypex - 1
 IF (lcore < 0) CALL mesage (-8,0,parm(3))
 nfile  = fileu(1)
 parm(2) = nfile
 CALL gopen (nfile,iobuf,rdrew)
 xmin = 1.0E20
 IF (itypex /= 1) dxmin = 1.0D20
 DO  i = 1,nrow
   iunpak = i
   junpak = i
   ind    = nrow - i + 1
   IF (itypex /= 1) GO TO 70
   CALL unpack (*400,nfile,x(ind))
   IF (xmin > ABS(x(ind))) xmin = ABS(x(ind))
   CYCLE
   70 CALL unpack (*400,nfile,dx(ind))
   IF (dxmin > DABS(dx(ind))) dxmin = DABS(dx(ind))
 END DO
 IF (itypex == 1 .AND. xmin /= 0.0  ) GO TO 120
 IF (itypex /= 1 .AND. dxmin /= 0.0D0) GO TO 120
 xmin = 1.0E20
 IF (itypex /= 1) dxmin = 1.0D20
 DO  i = 1,nrow
   IF (itypex /= 1) GO TO 90
   IF (x(i) == 0.0) CYCLE
   IF (xmin > ABS(x(i))) xmin = ABS(x(i))
   CYCLE
   90 IF (dx(i) == 0.0D0) CYCLE
   IF (dxmin > DABS(dx(i))) dxmin = DABS(dx(i))
 END DO
 IF (itypex /= 1) GO TO 110
 IF (xmin > 1.0E-8) xmin = 1.0E-8
 GO TO 120
 110 IF (dxmin > 1.0D-8) dxmin = 1.0D-8
 
!     BUILD LOAD VECTOR FOR BACKWARD PASS
 
 120 sdiag = 1.0
 IF (itypex /= 1) dsdiag = 1.0D0
 DO  i = 1,nrow
   anum = (-1)**(i*kcount)
   ai   = i
   aden = 1.0 + (1.0 - ai/nrow)*kcount
   avalue = anum/aden
   IF (itypex /=    1) GO TO 140
   130 x(i) = xmin*avalue/sdiag
   CYCLE
   140 dx(i) = dxmin*avalue/dsdiag
 END DO
 
!     BEGIN BACKWARD PASS
 
 CALL REWIND ( fileu )
 CALL skprec ( fileu, 1 )
 j = nrow
 390 CALL intpk (*650,fileu(1),0,iprec,0)
 IF (eol == 0.0) THEN
   GO TO   410
 ELSE
   GO TO   650
 END IF
 410 CALL zntpki
 i = nrow - ii + 1
 IF (i /= j) GO TO 510
 
!     DIVIDE BY THE DIAGONAL
 
 in1 = i
 k   = 0
 420 SELECT CASE ( iprec )
   CASE (    1)
     GO TO 430
   CASE (    2)
     GO TO 440
 END SELECT
 430 CONTINUE
 IF ( a(1) >= 0.0 .AND. ABS( a(1)) < xmin ) a(1) =  xmin
 IF ( a(1) < 0.0 .AND. ABS( a(1)) < xmin ) a(1) = -xmin
 x(in1) = x(in1)/a(1)
 GO TO 470
 440 CONTINUE
 IF ( da(1) >= 0.0D0 .AND. DABS(da(1)) < dxmin) da(1) =  dxmin
 IF ( da(1) < 0.0D0 .AND. DABS(da(1)) < dxmin) da(1) = -dxmin
 dx(in1) = dx(in1)/da(1)
 470 GO TO 490
 
!     SUBTRACT OFF REMAINING TERMS
 
 480 IF (i > j) GO TO 410
 490 IF (eol == 0.0) THEN
   GO TO   500
 ELSE
   GO TO   590
 END IF
 500 CALL zntpki
 i   = nrow - ii + 1
 510 in1 = i
 in2 = j
 IF (i < j) GO TO 530
 k   = in1
 in1 = in2 - ioff
 in2 = k
 530 SELECT CASE ( iprec )
   CASE (    1)
     GO TO 540
   CASE (    2)
     GO TO 550
 END SELECT
 540 x(in1) = x(in1) - a(1)*x(in2)
 GO TO 580
 550 dx(in1) = dx(in1) - dx(in2)*da(1)
 580 in1 = in1 + nrow
 in2 = in2 + nrow
 GO TO 480
 590 j = j - 1
 IF (j > 0) GO TO 390
! END OF BACKWARD SUBSTITUTION, NEGATE TERMS AND RETURN
 SELECT CASE ( iprec )
   CASE (    1)
     GO TO 600
   CASE (    2)
     GO TO 620
 END SELECT
 600 DO  k = 1, nrow
   x(k)  = -x(k)
 END DO
 GO TO 700
 620 DO  k = 1, nrow
   dx(k) = -dx(k)
 END DO
 700 CONTINUE
 CALL CLOSE ( fileu, rew )
 RETURN
 
!     ATTEMPT TO OPERATE ON SINGULAR MATRIX
 
 400 parm(1) = -5
 CALL mesage (parm(1),parm(2),parm(3))
 650 CALL mesage ( -5    ,parm(2),parm(3))
 RETURN
END SUBROUTINE detgbs
