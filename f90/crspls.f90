SUBROUTINE crspls (*,jump,mu,bp,rs,again,n23)
     
!     THIS ROUTINE HANDLES CRBE3 AND CRSPLINE RIGID ELEMENTS
!     CALLED ONLY BY CRIGGP SUBROUTINE
 
!     SINGLE PRECISION VERSION
 
 IMPLICIT INTEGER (a-z)
 INTEGER, INTENT(IN OUT)                  :: jump
 INTEGER, INTENT(OUT)                     :: mu
 INTEGER, INTENT(IN)                      :: bp
 INTEGER, INTENT(OUT)                     :: rs(3)
 LOGICAL, INTENT(IN OUT)                  :: again
 INTEGER, INTENT(OUT)                     :: n23
 LOGICAL :: debug
 INTEGER :: mcode(2),NAME(2),sild(6),retn1,retn2,retn3,retn4,retn
 REAL :: z(1),zk,wt,dl,coeff
 REAL :: x1,x2,x3,y1,y2,y3,z1,z2,z3,a(3),b(3),c(3),d(9),  &
     LEN,leng,one,zero,half,eps,espx,ans,di,fac,ln3,  &
     t(36),tx(36),knn(36),gnn(36),unn(36),znn(36), snn(36),x(36),y(36),w(6)
 CHARACTER (LEN=23) :: ufm
 COMMON /xmssg /  ufm
 COMMON /system/  sysbuf,nout
 COMMON /gp4fil/  geomp,bgpdt,cstm,rgt
 COMMON /gp4prm/  buf(20),buf1,buf2,buf3,buf4,knkl1,two16,nogo, gpoint,kn
 COMMON /crsply/  x1,y1,z1,x2,y2,z2,x3,y3,z3
 COMMON /zzzzzz/  iz(1)
 EQUIVALENCE      (z(1),iz(1)), (wt,iwt ), (dl,idl ),  &
     (x1  ,a(1) ), (x2,b(1)), (x3,c(1))
 DATA    one,     zero,    half,    eps,     times, debug  /  &
     1.0,     0.0,     0.5,     1.0E-10, 0,     .false./
 DATA    cm,cn,   nogox,   mask15,  NAME            /  &
     6 ,12,   0,       32767,   4HCRSP,4HLS     /
 
 IF (again) RETURN 1
 IF (debug) WRITE (nout,10) kn,knkl1,bp,geomp,bgpdt,cstm,rgt,jump
 10 FORMAT ('0  CRSPLS DEBUG- KN,KNKL1,BP,GEOMP,BGPDT,CSTM,RGT,JUMP=',  &
     /3X,8I7)
 kn2 = kn/2
 CALL sswtch (38,l38)
 
!     UNIT MATRIX UNN
 
 DO  i = 2,35
   unn( i) = zero
 END DO
 unn( 1) = one
 unn( 8) = one
 unn(15) = one
 unn(22) = one
 unn(29) = one
 unn(36) = one
 
!     JUMP = 1 FOR CRBE3 DATA,  JUMP = 2 FOR CRSPLINE DATA
 
 IF (jump == 2) GO TO 400
 
!     READ CRBE3 DATA ON INPUT FILE
!     =============================
 
!     CLEAR WORKING SPACE
!     READ INPUT CARD, SAVE BEGINING POINTER, BEGN, AND
!     COUNT NUMBER OF WORDS READ, NWDS, IN FIRST PASS
!     (EACH INPUT CARD WILL BE READ TWICE)
 
 begn = 3
 30 pass = 1
 DO  i = 1,36
   t(i)   = zero
   knn(i) = zero
 END DO
 CALL READ (*1300,*1300,geomp,buf,3,0,flag)
 nwds = 3
 IF (debug .OR. l38 == 1) WRITE (nout,50) buf(1)
 50 FORMAT (5X,'ELEMENT',i8,' IS BEING PROCESSED')
 eid    = buf(1)
 gpoint = buf(2)
 ASSIGN 60   TO retn
 ASSIGN 1000 TO retn1
 kx = buf(3)
 nm = cn
 GO TO 950
 60 refg = k
 sil  = gpoint
 DO  i = 1,6
   sild(i) = sil + i - 1
 END DO
 x2 = z(k+1)
 y2 = z(k+2)
 z2 = z(k+3)
 
!     READ WEIGHT FACTORS AND COMPONENTS.
!     GENERATE WEIGHT VECTOR W
 
 80 CALL READ (*1110,*1110,geomp,iwt,1,0,flag)
 IF (pass  ==  1) nwds = nwds + 1
 IF (iwt   == -2) GO TO 170
 IF (iwt   == -3) GO TO 240
 CALL READ (*1110,*1110,geomp,comp,1,0,flag)
 IF (pass == 1) nwds = nwds + 1
 ASSIGN 90 TO retn1
 kx = comp
 nm = cm
 GO TO 950
 90 DO  i = 1,6
   w(i) = zero
   IF (buf(cm+i) /= 0) w(i) = wt
 END DO
 
!     READ GRID POINT, GET TRANSFORMATION MATRIX, AND SUMMING UP
!     WT MATRIX, AND FINALLY KNN MATRIX
 
 110 CALL READ (*1110,*1110,geomp,grid,1,0,flag)
 IF (pass ==  1) nwds = nwds + 1
 IF (grid == -1) GO TO 80
 ASSIGN 120 TO retn
 gpoint = grid
 GO TO 1000
 120 x1 = z(k+1)
 y1 = z(k+2)
 z1 = z(k+3)
 ASSIGN 850 TO retn2
 ASSIGN 130 TO retn3
 zk = z(k)
 GO TO 800
 130 CALL gmmats (t,6,6,0, unn,6,6,0, x)
 IF (pass == 2) GO TO 270
 DO  i = 1,36
   tx(i) = x(i)
 END DO
 l = 0
 DO  i = 1,6
   DO  j = 1,6
     x(l+j) = x(l+j)*w(i)
   END DO
   l = l + 6
 END DO
 CALL gmmats (tx,6,6,-1, x,6,6,0, knn)
 
!     REPEAT FOR MORE GRID POINT
 
 GO TO 110
 
!     UM SET WAS SPECIFIED BY USER. REBUILD SILD WITH THE UM SET, AND
!     CHECK TOTAL NUMBER OF COMPONENTS FOR POSSIBLE ERROR
 
 170 IF (pass == 2) GO TO 310
 jj = 1
 180 CALL READ (*1110,*1110,geomp,grid,1,0,flag)
 nwds = nwds + 1
 IF (grid == -3) GO TO 240
 CALL READ (*1110,*1110,geomp,comp,1,0,flag)
 nwds = nwds + 1
 ASSIGN 190 TO retn1
 kx = comp
 nm = cm
 GO TO 950
 190 gpoint = grid
 ASSIGN 200 TO retn
 GO TO 1000
 200 DO  i = 1,6
   IF (buf(cm+ i) == 0) CYCLE
   IF (jj > 6) GO TO 1160
   210 IF (buf(cn+jj) /= 0) GO TO 220
   jj = jj + 1
   IF (jj > 6) GO TO 1160
   GO TO 210
   220 sild(jj) = gpoint + i - 1
   jj = jj + 1
 END DO
 GO TO 180
 240 IF (pass == 2) GO TO 320
 
!     STORE DIAG TERMS WITH -1.
!     ADD DEPENDENT SIL TO THE END OF OPEN CORE VIA MU POINTER
 
 DO  i = 1,6
   IF (buf(cn+i) == 0) CYCLE
   mcode(1) = sil + i - 1
   mcode(2) = sild(i)
   coeff    = -1.
   CALL WRITE (rgt,mcode,2,0)
   CALL WRITE (rgt,coeff,1,0)
   iz(mu) = mcode(2)
   mu = mu - 1
 END DO
 
!     GET MATRIX READY FOR SECOND PASS, IN TX
 
 sing = -1
 CALL invers (6,knn,6,0,0,LEN,sing,x)
 IF (sing == 2) GO TO 1120
 ASSIGN 260 TO retn2
 zk = z(refg)
 GO TO 800
 260 CALL gmmats (knn,6,6,0, t,6,6,0, tx)
 
!     BACK RECORD FOR 2ND PASS
!     SKIP TO WHERE WEIGHT FACTORS BEGIN
 
 CALL bckrec (geomp)
 pass = 2
 i = begn + 3
 CALL READ (*1110,*1110,geomp,buf,-i,0,flag)
 GO TO 80
 
!     INSERT THIS GRID MPC EQUATIONS
 
 270 CALL gmmats (tx,6,6,0, x,6,6,1, knn)
 DO  i = 1,6
   DO  j = 1,31,6
     l = i + j - 1
     knn(l) = knn(l)*w(i)
   END DO
 END DO
 DO  i = 1,6
   IF (buf(cn+i) == 0) CYCLE
   sil = sild(i)
   l = (i-1)*6
   DO  j = 1,6
     IF (buf(cm+j) == 0) CYCLE
     ans = knn(l+j)
     IF (ans == zero) CYCLE
     mcode(1) = gpoint + j - 1
     mcode(2) = sil
     coeff    = ans
     CALL WRITE (rgt,mcode,2,0)
     CALL WRITE (rgt,coeff,1,0)
   END DO
 END DO
 GO TO 110
 
!     SKIP TO END OF CARD
 
 310 CALL READ (*1110,*1110,geomp,j,1,0,flag)
 IF (j /= -3) GO TO 310
 
!     UPDATE BEGIN POINTER, AND RETURN FOR ANOTHER RBE3 CARD
 
 320 begn = begn + nwds
 GO TO 30
 
 
!     READ CRSPLINE DATA ON INPUT FILE
!     ================================
 
!     INPUT DATA WILL BE SAVED IN RS ARRAY
!     3 WORDS SAVED FOR EACH GRID - BGPDT POINTER, COMPONENT, AND SIL
 
 400 CALL READ (*1300,*1300,geomp,buf,3,0,flag)
 eid = buf(1)
 idl = buf(2)
 rs(1) = buf(3)
 rs(2) =-1
 rs(3) = 0
 IF (debug .OR. l38 == 1) WRITE (nout,50) buf(1)
 k = 4
 410 CALL READ (*1110,*1110,geomp,rs(k),2,0,flag)
 IF (rs(k) == -1) GO TO 420
 rs(k+2) = 0
 k = k + 3
 IF (k > mu) CALL mesage (-8,0,NAME)
 GO TO 410
 
!     END OF INPUT FOR THIS RIGID ELEMENT, NOW COMPUTE LENGTH, INTERNAL
!     NUMBER (BGPDT POINTER), AND CHANGE GRID TO SIL
 
 420 IF (k < 8) GO TO 1100
 IF (debug) CALL bug1 ('RS-     ',310,rs,k)
 iend = k - 1
 LEN  = zero
 ASSIGN 430 TO retn
 
!     DO 460 I = 1,IEND,3
 i = 1
 425 gpoint = rs(i)
 GO TO 1000
 
!     UPON RETURN FROM 1000, K IS BGPDT AND GPOINT IS SIL
 
 430 rs(i  ) = k
 rs(i+2) = gpoint
 IF (debug) WRITE (nout,440) i,gpoint,k,z(k+1)
 440 FORMAT (/10X,'@430  I, NEW GPOINT & K=',i4,2I6,e11.3)
 
 IF (i /= 1) GO TO 450
 x1 = z(k+1)
 y1 = z(k+2)
 z1 = z(k+3)
 GO TO 460
 450 x2 = z(k+1)
 y2 = z(k+2)
 z2 = z(k+3)
 LEN= LEN + SQRT((x2-x1)**2 + (y2-y1)**2 + (z2-z1)**2)
 x1 = x2
 y1 = y2
 z1 = z2
 460 i  = i + 3
 IF (i < iend) GO TO 425
 
 di = LEN*dl
 is = 1
 IF (.NOT.debug) GO TO 480
 CALL bug1 ('RS-     ',345,rs,iend)
 WRITE  (nout,470) LEN,di
 470 FORMAT ('0  LEN,DI/@470 =',2E14.5)
 
!     COMPUTATION FOR EACH SEPARATED SPLINE
!     SET NUMBER OF SEGMENTS, NS
 
 480 ns = 0
 DO  i = is,iend,3
   IF (rs(i+1) == 0) GO TO 500
   ns = ns + 1
 END DO
 
!     IB = BEGIN,  IE = END,  IS = PRESENT SEGMENT
!     ND = NUMBER OF DEPENDENT POINTS
 
!     ZERO MTRAIX WORKING SPACE KNN ,GNN, AND T
 
 500 ie = i
 ib = is
 nd = ns - 1
 DO  i = 1,36
   knn(i) = zero
   gnn(i) = zero
   t(i)   = zero
 END DO
 
!     COMPUTE FLEXIBILITY MATRIX ZNN AND ITS INVERSE KNN, FOR EACH
!     SPLINE SEGMENT.
 
!     DO 540 I = 1,NS
 i = 0
 515 i = i + 1
 i1 = rs(is  )
 i2 = rs(is+3)
 ASSIGN 520 TO retn4
 GO TO 900
 520 IF (nogox == 1) GO TO 540
 x2 = z(i1+1)
 y2 = z(i1+2)
 z2 = z(i1+3)
 x1 = z(i2+1)
 y1 = z(i2+2)
 z1 = z(i2+3)
 x2 = (x2+x1)*half
 y2 = (y2+y1)*half
 z2 = (z2+z1)*half
 i1 = rs(ie)
 x1 = z(i1+1)
 y1 = z(i1+2)
 z1 = z(i1+3)
 
!     FORM UNN USING BASIC UNN MATRIX
!     DO NOT DESTROY RIGID TRANSFER MATRIX
 
 ASSIGN 530 TO retn3
 GO TO 850
 530 CALL gmmats (unn,6,6,0, znn,6,6,0, snn)
 
!     SUM INTO KNN
 
 CALL gmmats (snn,6,6,-2, unn,6,6,1, knn)
 540 is = is + 3
 IF (i < ns) GO TO 515
 
 IF (nogox == 1) GO TO 730
 
!     INVERT KNN
 
 sing = -1
 CALL invers (6,knn,6,0,0,LEN,sing,snn)
 IF (sing == 2) GO TO 1120
 
!     LOOP FOR FINAL CONSTRAINT EQUATIONS
 
 is = ib
 js = rs(is)
 ii = 0
 545 ii = ii + 1
 i1 = rs(is)
 id = is + 3
 i2 = rs(id)
 ASSIGN 550 TO retn4
 GO TO 900
 550 x1 = z(i2+1)
 y1 = z(i2+2)
 z1 = z(i2+3)
 x2 = z(i1+1)
 y2 = z(i1+2)
 z2 = z(i1+3)
 x3 = z(i2+1)
 y3 = z(i2+2)
 z3 = z(i2+3)
 x2 = (x2+x3)*half
 y2 = (y2+y3)*half
 z2 = (z2+z3)*half
 
!     Y I+1 I X   S I+1 S
 
 ASSIGN 560 TO retn3
 GO TO 850
 560 CALL gmmats (unn,6,6,0, znn,6,6,0, snn)
 CALL gmmats (snn,6,6,0, unn,6,6,1,   y)
 x2 = z(i1+1)
 y2 = z(i1+2)
 z2 = z(i1+3)
 
!     S I+1 I X GIN
 
 ASSIGN 570 TO retn3
 GO TO 850
 570 CALL gmmats (unn,6,6,0, gnn,6,6,0, snn)
 i3 = rs(ie)
 x3 = z(i3+1)
 y3 = z(i3+2)
 z3 = z(i3+3)
 x2 = z(js+1)
 y2 = z(js+2)
 z2 = z(js+3)
 
!     GNN = G I+1 N
 
 ASSIGN 580 TO retn3
 GO TO 860
 580 CALL gmmats (y  ,6,6,0, unn,6,6,1, znn)
 CALL gmmats (znn,6,6,0, knn,6,6,0, gnn)
 DO  j = 1,36
   gnn(j) = gnn(j) + snn(j)
 END DO
 
!     Y = G I+1 1
 
 ASSIGN 600 TO retn3
 GO TO 870
 600 CALL gmmats (gnn,6,6,0, unn,6,6,0, snn)
 ASSIGN 610 TO retn3
 GO TO 850
 610 DO  j = 1,36
   y(j) = unn(j) - snn(j)
 END DO
 
!     TRANSFORM TO GLOBAL AND STORE ANSWERS IN Y AND SNN
 
 ASSIGN 630 TO retn2
 zk = z(i2)
 GO TO 800
 630 CALL gmmats (t,6,6,1, y  ,6,6,0, snn)
 CALL gmmats (t,6,6,1, gnn,6,6,0, znn)
 ASSIGN 640 TO retn2
 zk = z(js)
 GO TO 800
 640 CALL gmmats (snn,6,6,0, t,6,6,0, y)
 ASSIGN 650 TO retn2
 zk = z(i3)
 GO TO 800
 650 CALL gmmats (znn,6,6,0, t,6,6,0, snn)
 
!     Y = G I 1  SNN = G I N
 
 ASSIGN 660 TO retn1
 kx = rs(id+1)
 nm = cm
 GO TO 950
 
!     ADD DEPENDENT TO LIST AND MPC EQUATIONS TO RGT
 
 660 IF (.NOT.debug) GO TO 680
 WRITE (nout,670) y
 WRITE (nout,670) snn
 670 FORMAT ('0  CRSPLS/@670',/,(2X,10E12.4))
 680 DO  j = 1,6
   IF (buf(cm+j) == 0) CYCLE
   
!     SELF TERM FOR DEPENDENT SIL
   
   sil = rs(id+2) + j - 1
   mcode(1) = sil
   mcode(2) = sil
   coeff    = -1.
   CALL WRITE (rgt,mcode,2,0)
   CALL WRITE (rgt,coeff,1,0)
   iz(mu) = mcode(2)
   mu = mu - 1
   IF (ii >= mu) CALL mesage (-8,0,NAME)
   ll = (j-1)*6
   
!     END ONE DEPENDENT
   
   DO  l = 1,6
     ans = y(ll+l)
     
!     TEST FOR COMPUTED ZERO
     
     espx = eps
     IF (j > 3 .AND. l < 4) espx = espx/leng
     IF (j < 4 .AND. l > 3) espx = espx*leng
     IF (ABS(ans) < espx) CYCLE
     mcode(1) = rs(ib+2) + l - 1
     mcode(2) = sil
     coeff    = ans
     CALL WRITE (rgt,mcode,2,0)
     CALL WRITE (rgt,coeff,1,0)
   END DO
   
!     END N INDEPENDENT
   
   DO  l = 1,6
     ans = snn(ll+l)
     
!     TEST FOR COMPUTED ZERO
     
     espx = eps
     IF (j > 3 .AND. l < 4) espx = espx/leng
     IF (j < 4 .AND. l > 3) espx = espx*leng
     IF (ABS(ans) < espx) CYCLE
     mcode(1) = rs(ie+2) + l - 1
     mcode(2) = sil
     coeff    = ans
     CALL WRITE (rgt,mcode,2,0)
     CALL WRITE (rgt,coeff,1,0)
   END DO
 END DO
 
 is = is + 3
 IF (ii < nd) GO TO 545
 
!     END BIG DO (720) LOOP
 
 730 IF (ie+2 >= iend) GO TO 400
 is = ie
 rs(is+1) = -1
 GO TO 480
 
!     ----------------------------------------------------
 
!     INTERNAL ROUTINE TO BUILD 6X6 BASIC TO GLOBAL MATRIX
!     (T = 0 ON ENTRY)
 
 800 CALL transs (zk,d)
 j = 1
 DO  i = 1,15,6
   t(i   ) = d(j  )
   t(i+ 1) = d(j+1)
   t(i+ 2) = d(j+2)
   t(i+21) = d(j  )
   t(i+22) = d(j+1)
   t(i+23) = d(j+2)
   j = j + 3
 END DO
 GO TO retn2, (260,850,630,640,650)
 
!     INTERNAL ROUTINE TO MAKE RIGID BODY TRANSFER MATRIX FOR CRSPLINES
!     (UNN = IDENTITY MATRIX ON ENTRY)
 
 850 unn( 5) = a(3) - b(3)
 unn( 6) = b(2) - a(2)
 unn(10) = b(3) - a(3)
 unn(12) = a(1) - b(1)
 unn(16) = a(2) - b(2)
 unn(17) = b(1) - a(1)
 GO TO 880
 860 unn( 5) = c(3) - a(3)
 unn( 6) = a(2) - c(2)
 unn(10) = a(3) - c(3)
 unn(12) = c(1) - a(1)
 unn(16) = c(2) - a(2)
 unn(17) = a(1) - c(1)
 GO TO 880
 870 unn( 5) = c(3) - b(3)
 unn( 6) = b(2) - c(2)
 unn(10) = b(3) - c(3)
 unn(12) = c(1) - b(1)
 unn(16) = c(2) - b(2)
 unn(17) = b(1) - c(1)
 880 GO TO retn3, (130,530,560,570,580,600,610)
 
!     INTERNAL ROUTINE TO FORM FLEXIBILITY MATRIX FOR CRSPLINE
 
 900 DO  i = 1,36
   znn(i) = zero
 END DO
 x1 = z(i1+1)
 y1 = z(i1+2)
 z1 = z(i1+3)
 x2 = z(i2+1)
 y2 = z(i2+2)
 z2 = z(i2+3)
 leng = SQRT((x2-x1)**2 + (y2-y1)**2 + (z2-z1)**2)
 IF (leng == zero) GO TO 930
 znn(22) = leng
 znn(29) = leng
 znn(36) = leng
 fac = leng/12.0*((3.0*di**2)/(2.0*leng**2)-one)
 ln3 = leng**3/12.0
 znn( 1) = ln3 + fac*(x2-x1)**2
 znn( 2) = fac *(x2-x1)*(y2-y1)
 znn( 3) = fac *(x2-x1)*(z2-z1)
 znn( 7) = znn(2)
 znn( 8) = ln3 + fac*(y2-y1)**2
 znn( 9) = fac*(y2-y1)*(z2-z1)
 znn(13) = znn(3)
 znn(14) = znn(9)
 znn(15) = ln3 + fac*(z2-z1)**2
 920 GO TO retn4, (520,550)
 930 CALL mesage (30,31,eid)
 nogox = 1
 GO TO 920
 
!     INTERNAL ROUTINE TO ABSTRACT CODED DOF
 
 950 DO  i = 1,6
   buf(nm+i) = 0
 END DO
 IF (kx <= 0) GO TO 980
 DO  i = 1,6
   k1 = kx/10
   k2 = kx - k1*10
   IF (k2 > 6) GO TO 980
   buf(nm+k2) = k2
   IF (k1  == 0) GO TO 980
   kx = k1
 END DO
 980 GO TO retn1, (90,190,1000,660)
 
!     INTERNAL ROUTINE TO PERFORM BINARY SEARCH IN EQEXIN AND
!     CONVERT THE EXTERNAL NUMBER TO A SIL VALUE
 
 1000 klo = 0
 khi = kn2
 lastk = 0
 1010 k = (klo+khi+1)/2
 IF (lastk == k) GO TO 1140
 lastk = k
 IF (gpoint-iz(2*k-1) < 0.0) THEN
   GO TO  1020
 ELSE IF (gpoint-iz(2*k-1) == 0.0) THEN
   GO TO  1040
 ELSE
   GO TO  1030
 END IF
 1020 khi = k
 GO TO 1010
 1030 klo = k
 GO TO 1010
 1040 k = iz(2*k)
 gpoint = iz(k+2*kn)
 k = (k-1)*4 + bp
 IF (gpoint+5 > mask15) n23 = 3
 GO TO retn, (60,120,200,430)
 
!     ERROR MESSAGES
 
 1100 msg = 131
 GO TO 1130
 1110 CALL mesage (-3,geomp,NAME)
 1120 msg = 38
 1130 CALL mesage (30,msg,eid)
 GO TO 1180
 1140 WRITE  (nout,1150) ufm,gpoint,eid
 1150 FORMAT (a23,', UNDEFINED GRID POINT',i9,' SPECIFIED BY RIGID ',  &
     'ELEMENT ID',i9)
 times = times + 1
 IF (times > 50) CALL mesage (-37,0,NAME)
 GO TO 1180
 1160 WRITE  (nout,1170) ufm,eid
 1170 FORMAT (a23,', RIGID ELEMENT CRBE3',i9,' HAS ILLEGAL UM SET ',  &
     'SPECIFICATION')
 GO TO 1190
 
 1180 nogo  = 1
 nogox = 0
 SELECT CASE ( jump )
   CASE (    1)
     GO TO 30
   CASE (    2)
     GO TO 400
 END SELECT
 
!     REPOSITION GEOMP FILE FOR NEXT CRBE3 INPUT CARD
 
 1190 nogo  = 1
 nogox = 0
 CALL bckrec (geomp)
 i = begn + 1
 CALL READ (*1110,*1110,geomp,j,-i,0,flag)
 1200 CALL READ (*1110,*1110,geomp,j, 1,0,flag)
 i = i + 1
 IF (j /= -3) GO TO 1200
 begn = i
 GO TO 30
 
 1300 IF (nogox == 1) nogo = 1
 
 RETURN
END SUBROUTINE crspls
