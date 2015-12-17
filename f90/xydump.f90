SUBROUTINE xydump (outfil,TYPE)
     
 
 INTEGER, INTENT(IN OUT)                  :: outfil
 INTEGER, INTENT(IN)                      :: TYPE
 LOGICAL :: punch      ,plot       ,PRINT      ,outopn    ,  &
     null       ,random     ,ones       ,paplot    ,  &
     on         ,ok         ,intore     ,dec
 INTEGER :: oneone(2)  ,eor        , tcurve    ,  &
     xaxis      ,yaxis      ,ytaxis     ,ybaxis    ,  &
     curve      ,ycycle(2)  ,xcycle     ,center    ,  &
     z          ,limit(2,3) ,steps      ,vector    ,  &
     FILE       ,subc       ,vecid      ,buf       , begin      , two1
 REAL :: ymin(2)    ,ymax(2)    ,value(60)  ,rz(1)     ,  &
     rbuf(100)  ,idoutr(300),ylimit(2,3)
 COMMON /machin/ machx
 COMMON /BLANK / blkcom     ,vari(3)    ,nframe     ,ncard
 COMMON /two   / two1(32)
 COMMON /zzzzzz/ z(1)
 COMMON /xywork/ FILE       ,tcurve(32) ,ntops      ,PRINT     ,  &
     ifile      ,xaxis(32)  ,nbots      ,plot      ,  &
     vector     ,yaxis(32)  ,vecid(5)   ,punch     ,  &
     major      ,ytaxis(32) ,subc(5)    ,center    ,  &
     random     ,ybaxis(32) ,idin(153)  ,buf(100)  ,  &
     ivalue(60) ,iat        ,idout(300) ,outopn    ,  &
     steps      ,nat        ,paplot
 EQUIVALENCE     (limit(1,1),ylimit(1,1)) ,(z(1),rz(1)) ,  &
     (buf(1),rbuf(1)) ,(idout(1),idoutr(1)) , (ivalue(1),value(1))
 DATA   oneone / 1,1 /, eor/ 1 /, noeor/ 0 /
 DATA   npaplt / 0   /
 
!     SET MINIMUM X-DIFFERENCE
 
!     BUT FIRST CONVERT X FROM INTERGER TO REAL IF NECESSARY.
 
 intore = .false.
 dec = machx == 5 .OR. machx == 6 .OR. machx == 21
 j   = 1
 is1 = steps - 1
 
!     NOW SEARCH LIST FOR FIRST NON-ZERO ENTRY
 
 10 IF (z(iat+j) /= 0) GO TO 20
 j = j + 1
 IF (j > is1) GO TO 50
 GO TO 10
 
!              UNIVAC             CDC             CRAY
 20 IF (machx == 3 .OR. machx == 4 .OR. machx == 12)  THEN
      IF (IABS(z(iat+j))-two1(2) > 0.0) THEN
        GO TO    50
      ELSE
        GO TO    40
      END IF
    END IF 
   
!     IBM, VAX, UNIX
   
   IF (.NOT.dec .AND. IABS(z(iat+j)) > two1(9)) GO TO 50
   IF (     dec .AND. (z(iat+j) < 1 .OR. z(iat+j) > 127)) GO TO 50
   40 intore = .true.
   IF (j == 1) rz(iat+j) = z(iat+j)
   
   50 ok = .false.
   DO  i = 1,is1
     j = iat + i
     IF (intore) rz(j+1) = z(j+1)
     diff = rz(j+1) - rz(j)
     IF (.NOT.ok) GO TO 60
     IF (diff == 0.0) CYCLE
     xinc = AMIN1(xinc,diff)
     CYCLE
     60 IF (diff == 0.0) CYCLE
     xinc = diff
     ok   = .true.
   END DO
   IF (.NOT.ok) xinc = 1.0
   
!     SET XMIN AND XMAX FOR ALL DATA
   
   xcycle    = 0
   ycycle(1) = 0
   ycycle(2) = 0
   xmin = rz(iat+1)
   j    = iat + steps
   xmax = rz(j)
   
!     REDUCE THESE LIMITS TO USER SPECIFIED LIMITS
   
   IF (ivalue(1) /= 1) xmin = (value(1))
   IF (ivalue(2) /= 1) xmax = (value(2))
   
!     FURTHER EXPAND XLIMITS TO INCLUDE Y-AXIS INTERCEPT
   
   IF (ivalue(9) == 1) GO TO 80
   IF (ivalue(36) == 1 .AND. value(9) <= 0.0) GO TO 90
   xmin = AMIN1(xmin,value(9))
   xmax = AMAX1(xmax,value(9))
   
!     IF X-DIRECTION IS LOG AND XMIN IS NEGATIVE OR ZERO, SET YMIN
!     EQUAL TO THE SMALLEST NON-ZERO POSITIVE VALUE
   
   80 IF (ivalue(36) /= 1) GO TO 130
   90 IF (xmin    >  0.0) GO TO 120
   DO  i = 1,steps
     j = iat + i
     IF (rz(j) > 0.0) GO TO 110
   END DO
   xmin = 1.0
   xmax = 10.
   GO TO 120
   110 xmin = rz(j)
   120 CALL xylog (xmin,xmax,xcycle)
   
!     SWITCH XMIN AND XMAX (SAFETY CHECK) IF NECESSARY
   
   130 IF (xmin <= xmax) GO TO 140
   temp = xmin
   xmin = xmax
   xmax = temp
   
!     USING XMIN AND XMAX AS LIMITS DETERMINE Y-LIMITS FOR TOP AND
!     BOTTOM.
   
!     I1 = FIRST STEP WITHIN XMIN TO XMAX
!     I2 = LAST  STEP WITHIN XMIN TO XMAX
   
!     FIRST FIND I1 AND I2
   
   140 DO  i = 1,steps
     j = iat + i
     IF (xmin <= rz(j) .AND. rz(j) <= xmax) GO TO 160
   END DO
   i1 = 0
   GO TO 180
   160 i1 = i
   j  = iat + steps + 1
   DO  i = 1,steps
     j = j - 1
     IF (xmin <= rz(j) .AND. rz(j) <= xmax) GO TO 190
   END DO
   180 i2 = 0
   GO TO 200
   190 i2 = j - iat
   200 IF (i1 /= 0) GO TO 210
   
!     FIND FOLLOWING VALUES FOR CURVES AS A GROUP
   
!    YLIMIT(1,1)=YMIN TOP, YLIMIT(1,2)=YMAX TOP, YLIMIT(1,3)=MIN POS TOP
!    YLIMIT(2,1)=YMIN BOT, YLIMIT(2,2)=YMAX BOT, YLIMIT(2,3)=MIN POS BOT
   
   ymin(1) = 0.0
   ymin(2) = 0.0
   ymax(1) = 10.
   ymax(2) = 10.
   GO TO 330
   210 m = 1
   IF (nbots /= 0) m = 2
   begin = iat
   DO  i = 1,m
     limit(i,1) = 1
     limit(i,2) = 1
     limit(i,3) = 1
     DO  j = 1,ntops
       k  = j*steps + begin
       j1 = k + i1
       j2 = k + i2
       IF (limit(i,1) /= 1) GO TO 240
       
!     FIND FIRST NON-INTEGER 1 VALUE
       
       DO  k = j1,j2
         IF (z(k) /= 1) GO TO 230
       END DO
       CYCLE
       230 ylimit(i,1) = rz(k)
       ylimit(i,2) = rz(k)
       240 DO  k = j1,j2
         IF (z(k) == 1) CYCLE
         ylimit(i,1) = AMIN1(rz(k),ylimit(i,1))
         ylimit(i,2) = AMAX1(rz(k),ylimit(i,2))
         IF (rz(k) <= 0.0) CYCLE
         IF (limit(i,3) == 1) ylimit(i,3) = rz(k)
         ylimit(i,3) = AMIN1(ylimit(i,3),rz(k))
       END DO
     END DO
     begin = center
     
!     DEFAULT YLIMITS IF ALL CURVES NULL
     
     IF (limit(i,1) /= 1) GO TO 270
     ylimit(i,1) = 0.0
     ylimit(i,2) = 100.
     270 IF (limit(i,3) == 1) ylimit(i,3) = 10.0
     
   END DO
   
!     SET FINAL Y-LIMITS FOR UPPER AND LOWER CURVES
   
   
!     K=1 IMPLIES WHOLE CURVES
!     K=2 IMPLIES UPPER AND LOWER CURVES
   
   k = 1
   IF (nbots > 0) k = 2
   DO  i = 1,k
     ymin(i) = ylimit(i,1)
     ymax(i) = ylimit(i,2)
     
!     REDUCE THESE CURVE LIMITS TO LIMITS SET BY USER
     
     itemp = 2*(i+k)
     IF (ivalue(itemp-1) /= 1) ymin(i) = (value(itemp-1))
     IF (ivalue(itemp  ) /= 1) ymax(i) = (value(itemp  ))
     
!     FURTHER EXPAND LIMITS TO INCLUDE X-AXIS
     
     itemp = i + k
     IF (ivalue(itemp+8) == 1) GO TO 290
     IF (ivalue(itemp+35) == 1 .AND. value(itemp+8) <= 0.e0) GO TO 300
     ymin(i) = AMIN1(ymin(i),value(itemp+8))
     ymax(i) = AMAX1(ymax(i),value(itemp+8))
     
!     IF Y-DIRECTION IS LOG AND YMIN IS NEGATIVE OR ZERO SET YMIN
!     EQUAL TO SMALLEST POSITIVE CURVE VALUE WITHIN XLIMITS
     
     290 IF (ivalue(itemp+35) /= 1) GO TO 310
     300 IF (ymin(i) <= 0.0) ymin(i) = ylimit(i,3)
     CALL xylog (ymin(i),ymax(i),ycycle(i))
     
!     SWITCH YMIN AND YMAX (SAFETY CHECK) IF NECESSARY
     
     310 IF (ymin(i) <= ymax(i)) CYCLE
     temp    = ymin(i)
     ymin(i) = ymax(i)
     ymax(i) = temp
   END DO
   
!     ALL CURVE LIMITS HAVE NOW BEEN SET FOR THIS FRAME
   
   
!     OUTPUT EACH CURVE AND AN IDOUT RECORD IF PLOTS = .TRUE.
   
!     FILL IDOUT
   
   330 DO  i = 1,300
     idout(i) = 0
   END DO
   IF (plot .AND. outopn) nframe = nframe + 1
   idout(1) = subc(FILE)
   idout(2) = nframe
   idout(6) = vector
   idout(9) = ivalue(45)
   IF (ivalue(43) == 0) value(43) = 1.0
   idout(43)  = ivalue(43)
   idoutr(10) = xinc
   idout(245) = TYPE
   idout(246) = steps
   idoutr(282)= value(57)
   IF (idoutr(282) < 1.0) idoutr(282) = 1.0
   idout(283) = ivalue(50)
   IF (ivalue(47) == 3) idout(283) = ivalue(41)
   idout(284) = ivalue(47)
   idout(285) = ivalue(48)
   idout(286) = ivalue(49)
   idout(287) = ivalue(46)
   idout( 44) = ivalue(58)
   idout( 45) = ivalue(59)
   IF (PRINT) idout(288) = 1
   IF (plot ) idout(289) = 1
   IF (.NOT.paplot) GO TO 350
   IF (.NOT.plot) idout(289) = -1
   IF (     plot) idout(289) =  2
   npaplt = npaplt+1
   idout(281) = npaplt
   350 on = .false.
   IF (plot .OR. paplot) on = .true.
   IF (punch) idout(290) = 1
   DO  i = 51,146
     idout(i) = idin(i)
   END DO
   
!     BRANCH ON TOP, BOTTOM, OR WHOLE CURVE (FIRST WILL BE TOP OR WHOLE)
   
   i = 3
   IF (z(i) == 0 .OR. random) GO TO 400
   
!     TOP CURVE ID
   
   curve = 0
   idout(7) = 1
   idout(8) = 1
   idoutr(11) = xmin
   idoutr(12) = xmax
   idoutr(13) = ymin(1)
   idoutr(14) = ymax(1)
   iflag = 0
   IF (intore) iflag = 1
   CALL xytics (idout(15),idoutr(15),ivalue(17),idout(11) ,  &
       idout(12),ivalue(21),xcycle,iflag)
   CALL xytics (idout(23),idoutr(23),ivalue(19),idout(13),  &
       idout(14),ivalue(23),ycycle(1),0)
   idout(31) = ivalue(34) + ivalue(25)
   idout(32) = ivalue(34) + ivalue(26)
   idout(33) = ivalue(34) + ivalue(29)
   idout(34) = ivalue(34) + ivalue(30)
   idout(35) = xcycle
   idout(36) = ycycle(1)
   idout(37) = ivalue(15)
   idout(38) = ivalue(11)
   IF (idout(38)  == 1) idoutr(38) = 0.0
   IF (idoutr(38) < ymin(1)) idout(37) = 0
   idout(39) = ivalue(14)
   idout(40) = ivalue( 9)
   IF (idout(40)  == 1) idoutr(40) = 0.0
   IF (idoutr(40) < xmin) idout(39) = 0
   idout(41) = ivalue(40)
   idout(243) = ivalue(53)
   idout(244) = ivalue(54)
   DO  i=1,32
     idout(i+146) = tcurve(i)
     idout(i+178) = xaxis(i)
     idout(i+210) = ytaxis(i)
   END DO
   GO TO 420
   
!     BOTTOM CURVE ID (SET ONLY VALUES THAT CHANGE FROM THE TOP CURVES)
   
   380 curve = 0
   idout(7)   = -1
   idoutr(13) = ymin(2)
   idout(8)   = 1
   idoutr(14) = ymax(2)
   CALL xytics (idout(23),idoutr(23),ivalue(20),idout(13),  &
       idout(14),ivalue(24),ycycle(2),0)
   idout(31) = ivalue(35) + ivalue(25)
   idout(32) = ivalue(35) + ivalue(26)
   idout(33) = ivalue(35) + ivalue(31)
   idout(34) = ivalue(35) + ivalue(32)
   idout(36) = ycycle(2)
   idout(37) = ivalue(16)
   idout(38) = ivalue(12)
   IF (idout(38)  == 1) idoutr(38) = 0.0
   IF (idoutr(38) < ymin(2)) idout(37) = 0
   idout(243) = ivalue(55)
   idout(244) = ivalue(56)
   DO  i = 1,32
     idout(i+146) = tcurve(i)
     idout(i+178) = xaxis(i)
     idout(i+210) = ybaxis(i)
   END DO
   ipair = center + steps
   GO TO 430
   
!     WHOLE CURVE ID
   
   400 curve = 0
   idout(7) = 0
   idout(8) = 1
   idoutr(11) = xmin
   idoutr(12) = xmax
   idoutr(13) = ymin(1)
   idoutr(14) = ymax(1)
   iflag = 0
   IF (intore) iflag = 1
   CALL xytics (idout(15),idoutr(15),ivalue(17),idout(11),  &
       idout(12),ivalue(21),xcycle,iflag)
   CALL xytics (idout(23),idoutr(23),ivalue(18),idout(13),  &
       idout(14),ivalue(22),ycycle(1),0)
   idout(31) = ivalue(33) + ivalue(25)
   idout(32) = ivalue(33) + ivalue(26)
   idout(33) = ivalue(33) + ivalue(27)
   idout(34) = ivalue(33) + ivalue(28)
   idout(35) = xcycle
   idout(36) = ycycle(1)
   idout(37) = ivalue(13)
   idout(38) = ivalue(10)
   IF (idout(38)  == 1) idout(38) = 0.0
   IF (idoutr(38) < ymin(1)) idout(37) = 0
   idout(39) = ivalue(14)
   idout(40) = ivalue( 9)
   IF (idout(40)  == 1) idoutr(40) = 0.0
   IF (idoutr(40) < xmin) idout(39) = 0
   idout(41 ) = ivalue(40)
   idout(243) = ivalue(51)
   idout(244) = ivalue(52)
   DO  i=1,32
     idout(i+146) = tcurve(i)
     idout(i+178) = xaxis(i)
     idout(i+210) = yaxis(i)
   END DO
   GO TO 420
   
!     IDOUT IS COMPLETE   OUTPUT CURVES
   
   420 ASSIGN 590 TO icont
   ipair = iat + steps
   n = 1
   
   430 mcount = 0
   DO  m = 1,nat,3
     mcount = mcount + 1
     
!     CURVE NUMBER, ID, COMPONENT
     
     idout(4) = z(m)
     itemp = m + n
     idout(5) = z(itemp)
     IF (idout(5) /= 1000) curve = curve + 1
     idout(3) = curve
     
!     MEAN RESPONSE IN PLACE OF SUBCASE IF RANDOM
     
     IF (random) idout(1) = z(itemp+1)
     
!     SET NUMBER OF ZERO CROSSINGS IF RANDOM
     
     IF (random) idout(42) = buf(mcount+20)
     
!     COMPUTE Y1 = YMIN  AND Y2 = YMAX  FOR ALL DATA FOR THIS CURVE
     
     begin = ipair + mcount*steps - steps
     null  = .true.
     DO  k = 1,steps
       i = begin + k
       IF (z(i) == 1) CYCLE
       IF (.NOT.null  ) GO TO 440
       nx1 = k
       nx2 = k
       y1  = rz(i)
       y2  = rz(i)
       null= .false.
       CYCLE
       440 IF (rz(i) >= y1) GO TO 450
       y1  = rz(i)
       nx1 = k
       CYCLE
       450 IF (rz(i) <= y2) CYCLE
       y2  = rz(i)
       nx2 = k
     END DO
     
     IF (.NOT.null) GO TO 470
     idoutr(297) = 0.0
     idoutr(298) = 0.0
     idoutr(299) = 0.0
     idoutr(300) = 0.0
     GO TO 480
     470 nx1 = nx1 + iat
     nx2 = nx2 + iat
     idoutr(297) = y1
     idoutr(298) = rz(nx1)
     idoutr(299) = y2
     idoutr(300) = rz(nx2)
     
!     COMPUTE Y1 AND Y2 FOR DATA BETWEEN XMIN AND XMAX
     
     480 null = .true.
     IF (i1 == 0) GO TO 520
     DO  k = i1,i2
       i = begin + k
       IF (z(i) == 1) CYCLE
       IF (.NOT.null  ) GO TO 490
       nx1 = k
       nx2 = k
       y1  = rz(i)
       y2  = rz(i)
       null= .false.
       CYCLE
       490 IF (rz(i) >= y1) GO TO 500
       y1  = rz(i)
       nx1 = k
       CYCLE
       500 IF (rz(i) <= y2) CYCLE
       y2  = rz(i)
       nx2 = k
     END DO
     IF (.NOT.null) GO TO 530
     520 idoutr(293) = 0.0
     idoutr(294) = 0.0
     idoutr(295) = 0.0
     idoutr(296) = 0.0
     GO TO 540
     530 nx1 = nx1 + iat
     nx2 = nx2 + iat
     idoutr(293) = y1
     idoutr(294) = rz(nx1)
     idoutr(295) = y2
     idoutr(296) = rz(nx2)
     
     540 idoutr(291) = rz(iat+1)
     itemp = iat + steps
     idoutr(292) = rz(itemp)
     
!     IDOUT IS COMPLETE FOR THIS CURVE
     
     IF (idout(5) /= 0 .AND. idout(5) /= 1000)  &
         CALL xyout (-1,idout(1),idoutr(1))
     IF (on) CALL WRITE (outfil,idout(1),300,eor)
     idout(8) = 0
     
!     DUMP ALL PAIRS TO PRINTER AND PUNCH,  THOSE IN RANGE TO PLOTTER
     
     y1 = idoutr(13)
     y2 = idoutr(14)
     IF (on) CALL WRITE (outfil,oneone(1),2,noeor)
     ones = .true.
     IF (idout(5) == 1000) GO TO 570
     DO  k = 1,steps
       i = begin + k
       j = iat + k
       buf(1) = z(j)
       buf(2) = z(i)
       IF (z(i) == 1) CYCLE
       IF (k < i1 .OR. k > i2) CYCLE
       IF (PRINT .OR. punch) CALL xyout (1,buf(1),rbuf(1))
       IF (rz(i) < y1 .OR. rz(i) > y2) GO TO 550
       IF (on) CALL WRITE (outfil,buf(1),2,noeor)
       ones = .false.
       CYCLE
       550 IF (ones) CYCLE
       IF (on) CALL WRITE (outfil,oneone(1),2,noeor)
       ones = .true.
     END DO
     570 IF (on) CALL WRITE (outfil,buf(1),0,eor)
   END DO
   
   GO TO icont, (590,600)
   
!     DO BOTTOM CURVES IF ANY
   
   590 ASSIGN 600 TO icont
   n = 2
   IF (idout(7) > 0) GO TO 380
   600 RETURN
 END SUBROUTINE xydump
