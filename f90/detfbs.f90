SUBROUTINE detfbs (iy,iobuf,fileu,nrow,kcount)
     
!     DETFBS IS A SPECIAL VERSION OF THE GFBS ROUTINE AND IS USED BY
!     THE REAL DETERMINANT METHOD.  IT IS SUITABLE FOR BOTH SINGLE
!     AND DOUBLE PRECISION OPERATION.
 
 
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
 INTEGER, INTENT(IN)                      :: fileu(7)
 INTEGER, INTENT(IN)                      :: nrow
 INTEGER, INTENT(IN)                      :: kcount
 INTEGER :: parm(4) , option ,sdet ,filev ,  &
     filevt  ,scr3    ,scr4    ,scr6   ,scr
 REAL :: x(1)    ,y(1)
 DOUBLE PRECISION :: dx(1)   ,dy(1)   ,dxmin   ,dsdiag
 COMMON /zzzzzz/  core(1)
 COMMON /detmx /  dum3(36),ipdeta
 COMMON /names /  rd      ,rdrew   ,wrt     ,wrtrew  ,rew
 COMMON /regean/  dum1(23),scr3    ,scr4    ,dum2(11),scr6
 COMMON /reigkr/  option
 COMMON /trnspx/  filev(7),filevt(7),lcore  ,ncr     ,scr(2)
 COMMON /unpakx/  itypex  , iunpak ,junpak  ,incr
 EQUIVALENCE      (core(1),x(1),dx(1),y(1),dy(1))    ,  &
     (xmin,dxmin) ,   (sdiag,dsdiag)
 DATA    sdet  /  4HSDET  /
 DATA    parm(3), parm(4) / 4HDETF, 4HBS   /
 
!WKBI SPR 94011 10/94
 IF ( option /= sdet ) GO TO 1000
 itypex = fileu(5)
 INDEX  = -1
 incr   = 1
 nfile  = fileu(1)
 IF (option == sdet) GO TO 30
 INDEX  = 1
 lcore  = ipdeta - iy*itypex - 1
 IF (lcore < 0) CALL mesage (-8,0,parm(3))
 ncr = 2
 scr(1) = scr3
 scr(2) = scr4
 DO  i = 1,7
   filev(i)  = fileu(i)
   filevt(i) = fileu(i)
 END DO
 30 filevt(1) = scr6
 nfile  = filevt(1)
 IF (itypex == 1) CALL trnsp ( y(iy))
 IF (itypex /= 1) CALL trnsp (dy(iy))
 IF (itypex == 1) GO TO 50
 ASSIGN 230 TO isd
 ASSIGN 260 TO ius
 40 parm(2) = nfile
 CALL gopen (nfile,iobuf,rdrew)
 GO TO 60
 50 ASSIGN 240 TO isd
 ASSIGN 270 TO ius
 GO TO 40
 60 xmin = 1.0E20
 IF (itypex /= 1) dxmin = 1.0D20
 DO  i = 1,nrow
   iunpak = 0
   IF (itypex /= 1) GO TO 70
   CALL unpack (*400,nfile,x(i))
   IF (xmin > ABS(x(i))) xmin = ABS(x(i))
   CYCLE
   70 CALL unpack (*400,nfile,dx(i))
   IF (dxmin > DABS(dx(i))) dxmin = DABS(dx(i))
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
   IF (option /= sdet) GO TO 130
   sdiag = x(i)
   IF (x(i) >= 0.0 .AND. ABS(x(i)) < xmin) sdiag = xmin
   IF (x(i) < 0.0 .AND. ABS(x(i)) < xmin) sdiag =-xmin
   130 x(i) = xmin*avalue/sdiag
   CYCLE
   140 IF (option /= sdet) GO TO 150
   dsdiag = dx(i)
   IF (dx(i) >= 0.0 .AND. DABS(dx(i)) < dxmin) dsdiag = dxmin
   IF (dx(i) < 0.0 .AND. DABS(dx(i)) < dxmin) dsdiag =-dxmin
   150 dx(i) = dxmin*avalue/dsdiag
 END DO
 
 
!     BEGIN BACKWARD PASS
 
 DO  i = 1,nrow
   iunpak = 0
   j = nrow - i + 1
   CALL bckrec (nfile)
   IF (itypex == 1) CALL unpack (*400,nfile,y(iy))
   IF (itypex /= 1) CALL unpack (*400,nfile,dy(iy))
   CALL bckrec (nfile)
   ising = 0
   k = junpak - iunpak + iy
   GO TO isd, (230,240)
   
!     DIVIDE BY THE DIAGONAL TERM
   
   200 IF (option == sdet) CYCLE
   IF (dy(k) >= 0.0D0 .AND. DABS(dy(k)) < dxmin) dy(k) = dxmin
   IF (dy(k) < 0.0D0 .AND. DABS(dy(k)) < dxmin) dy(k) =-dxmin
   dx(j) = dx(j)/dy(k)
   CYCLE
   210 IF (option == sdet) CYCLE
   IF (y(k) >= 0.0 .AND. ABS(y(k)) < xmin) y(k) = xmin
   IF (y(k) < 0.0 .AND. ABS(y(k)) < xmin) y(k) =-xmin
   x(j) = x(j)/y(k)
   CYCLE
   220 k = k - 1
   junpak = junpak - 1
   IF (k < iy) GO TO 280
   GO TO isd, (230,240)
   230 IF (dy(k) == 0.0D0) GO TO 220
   IF (junpak - j < 0) THEN
     GO TO   280
   ELSE IF (junpak - j == 0) THEN
     GO TO   200
   ELSE
     GO TO   250
   END IF
   240 IF (y(k) == 0.0) GO TO 220
   IF (junpak - j < 0) THEN
     GO TO   280
   ELSE IF (junpak - j == 0) THEN
     GO TO   210
   END IF
   250 GO TO ius, (260,270)
   260 dx(j) = dx(j) - INDEX*dx(junpak)*dy(k)
   GO TO 220
   270 x(j) = x(j) - INDEX*x(junpak)*y(k)
   GO TO 220
   280 IF (ising == 0) GO TO 400
 END DO
 
 IF (option == sdet) GO TO 340
 IF (itypex ==    1) GO TO 320
 DO  i = 1,nrow
   dx(i) = -dx(i)
 END DO
 GO TO 340
 320 DO  i = 1,nrow
   x(i) = -x(i)
 END DO
 340 CALL CLOSE (nfile,rew)
!WKBI 10/94 SPR94011
 700 CONTINUE
 RETURN
!WKBNB 10/94 SPR94011
 1000  CALL detgbs( iy, iobuf, kcount )
 GO TO 700
!WKBNE 10/94 SPR94011
 
!     ATTEMPT TO OPERATE ON SINGULAR MATRIX
 
 400 parm(1) = -5
 CALL mesage (parm(1),parm(2),parm(3))
 RETURN
END SUBROUTINE detfbs
