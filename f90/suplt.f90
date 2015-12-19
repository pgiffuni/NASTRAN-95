SUBROUTINE suplt (iz,iy,x,u,gplst,pen,deform)
     
!     TO CREATE A SET OF UNIQUE LINES TO BE PLOTTED.  IN ADDITION, TO
!     AVOID SKIPPING ALL OVER THE PLOT FOR EACH LINE.
 
!     INPUT (SIMULAR TO -GPCT-)
!       NGRID - NUMBER OF INTERNAL GRID POINTS.
!       IY    - 1 THRU NGRID - POINTERS TO FIRST CONNECTION OF THE
!                              INTERNAL GRID MATCHING THIS INDEX.
!                              IF THE GRID HAS NO ENTRIES IZ(I+1) WILL
!                              HAVE THE SAME POINTER VALUE.
!             - NGRID+1      - POINTER TO END-OF-RECORD.
!       IZ    - CONNECTING INTERNAL GRIDS.  POINTER FOR NEXT GRID
!                              DETERMINES LAST ENTRY.  ENTRIES PUSHED
!                              DOWN AND  -1  ADDED AT END AS EACH ENTRY
!                              IS USED.
 
!       NTAB  = TOTAL ENTRY COUNTER IN GPCT
!       ID1   = START OF CURRENT -LINE-
!       ID2   = END   OF CURRENT -LINE-
!       ID3   = START OF LAST -LINE-
!       ID4   = END   OF LAST -LINE-
 
 
 INTEGER, INTENT(IN OUT)                  :: iz(1)
 INTEGER, INTENT(IN)                      :: iy(1)
 REAL, INTENT(IN)                         :: x(3,1)
 REAL, INTENT(IN)                         :: u(2,1)
 INTEGER, INTENT(IN OUT)                  :: gplst(1)
 INTEGER, INTENT(IN OUT)                  :: pen
 INTEGER, INTENT(IN OUT)                  :: deform
 INTEGER :: nm(5),m(71), m1(17),m2(17),m3(11),m4(17),m5(9),ERR(4),lm(5),  &
     limt(2)
 
 COMMON /system/ sysbuf,iout
 COMMON /BLANK / ngrid,skp1(19),merr
 EQUIVALENCE     (m1(1),m( 1)), (m2(1),m(18)), (m3(1),m(35)),  &
     (m4(1),m(46)), (m5(1),m(63))
 DATA    nm    / 17,17,11,17, 9 /, lm    /  1,18,35,46,63 /,  &
     m1    / 4H(35X, 4H,26H, 4HSUPL, 4HT re, 4HJECT, 4HED p,  &
     4HLOT., 4H piv, 4HOT,i, 4H8,26, 4HH is, 4H zer,  &
     4HO OR, 4H sam, 4HE as, 4H ent, 4HRY.)/,  &
     m2    / 4H(35X, 4H,54H, 4HSUPL, 4HT re, 4HJECT, 4HED p,  &
     4HLOT., 4H neg, 4HATIV, 4HE nu, 4HMBER, 4H ent,  &
     4HRIES, 4H - n, 4H3,n4, 4H =,2, 4HI10)/,  &
     m3    / 4H(35X, 4H,31H, 4HUNEX, 4HPECT, 4HED e, 4HOF i,  &
     4HN su, 4HPLT , 4H- pi, 4HVOT,, 4HI10)/,  &
     m4    / 4H(35X, 4H,11H, 4HSUPL, 4HT-en, 4HTRY,, 4HI10,,  &
     4H22H , 4HFOR , 4HPIVO, 4HT no, 4HT fo, 4HUND ,  &
     4H(,2I, 4H6,8H, 4H) ra, 4HNGE., 4H)   /,  &
     m5    / 4H(35X, 4H,24H, 4HNO e, 4HLEME, 4HNTS , 4HIN t,  &
     4HHIS , 4HSET., 4H)   /
 
 linesp = 0
 
!     LOCATE FIRST PIVOT (ID1) WITH ODD NUMBER OF ENTRIES
 
 id2 = 0
 id1 = 0
 DO  i = 1,ngrid
   ll  = iy(i+1) - iy(i)
   IF (ll == 0) CYCLE
   
!     IN CASE AN ODD NO. ENTRIES ISN'T FOUND THE DEFAULT IS FIRST PIVOT
   
   IF (id2 == 0) id2 = i
   IF (MOD(ll,2) == 0) CYCLE
   id1 = i
   i1  = iy(i)
   id2 = iz(i1)
   GO TO 20
 END DO
 
 20 ntab = iy(ngrid+1) - iy(1)
 
!     NO ELEMENTS IN THE SET
 
 IF (ntab == 0) GO TO 440
 
!     SEE IF ANY ODD ENTRIES FOUND
 
 limt(1) = 0
 limt(2) = ngrid + 1
 IF (id1 /= 0) GO TO 140
 id1 = id2
 i1  = iy(id1)
 id2 = iz(i1)
 GO TO 140
 
!     START OF LOOP AFTER FIRST -LINE-
 
 30 IF (n4 < 0) THEN
   GO TO   410
 ELSE IF (n4 == 0) THEN
   GO TO   260
 END IF
 
!     LAST END HAS ENTRY TO CONTINUE FROM
 
 40 id1 = id4
 
!     INTERNAL SEARCH FOR FIRST GPCT ENTRY ABOVE AND BELOW THE PIVOT
!     VALUE FOR THE PIVOT           *** M I N / M A X ***
 
 50 i1 = iy(id1)
 il =-100000
 ih = 100000
 j1 = iy(id1+1) - 1
 
 DO  i = i1,j1
   IF (iz(i)     < 0) THEN
     GO TO   100
   ELSE IF (iz(i)     == 0) THEN
     GO TO   400
   END IF
   60 IF (iz(i)-id1 < 0) THEN
     GO TO    70
   ELSE IF (iz(i)-id1 == 0) THEN
     GO TO   400
   ELSE
     GO TO    90
   END IF
   70 il = iz(i)
 END DO
 GO TO 100
 
 90 ih = iz(i)
 
!     DETERMINE  WHICH IS CLOSER TO PIVOT
 
 100 i = id1 - il
 j = ih  - id1
 IF (j-i < 0) THEN
   GO TO   130
 ELSE IF (j-i == 0) THEN
   GO TO   110
 ELSE
   GO TO   120
 END IF
 
!     EQUAL DISTANT, GO TO SAME DIRECTION AS BEFORE
 
 110 IF (id4-id3 < 0) THEN
   GO TO   120
 ELSE IF (id4-id3 == 0) THEN
   GO TO   400
 ELSE
   GO TO   130
 END IF
 
!     ID2 IS LESSOR ID
 
 120 id2 = il
 GO TO 140
 
!     ID2 IS GREATER ID
 
 130 id2 = ih
 
!     OUTPUT THE LINE -
!     NOTE THAT ID4 MAY BE RESET AT 320 SO DONT TAKE SHORTCUTS
 
 140 CONTINUE
 i  = IABS(gplst(id1))
 j  = IABS(gplst(id2))
 IF (deform /= 0) GO TO 160
 x1 = x(2,i)
 y1 = x(3,i)
 x2 = x(2,j)
 y2 = x(3,j)
 GO TO 170
 160 x1 = u(1,i)
 y1 = u(2,i)
 x2 = u(1,j)
 y2 = u(2,j)
 170 CONTINUE
 CALL line (x1,y1,x2,y2,pen,0)
 linesp = linesp + 1
 
!     REMOVE ENTRIES FROM CORE, LEFT SHIFT AS NEEDED AND PUT -1 AT THE
!     END OF THE TABLE.  PLACE THE NUMBER OF ENTRIES LEFT IN N3 AND N4.
!     DECREMENT THE TOTAL NUMBER OF ENTRIES BY 2. SET ID3 AND ID4.
 
 IF (ntab <= 2) GO TO 460
 kk = 0
 j1 = i1
 j2 = iy(id1+1) - 1
 ll = id2
 180 ipar = j1
 il = 0
 
 DO  i = j1,j2
   IF (iz(i)-ll < 0) THEN
     GO TO   210
   ELSE IF (iz(i)-ll == 0) THEN
     GO TO   200
   END IF
   190 iz(ipar) = iz(i)
   GO TO 220
   
!     COMPONENT TO BE ELIMINATED HAS BEEN FOUND
   
   200 il = 1
   CYCLE
   210 IF (iz(i) < 0) THEN
     GO TO   240
   ELSE IF (iz(i) == 0) THEN
     GO TO   400
   END IF
   220 ipar = ipar + 1
 END DO
 
 240 iz(ipar) = -ll
 IF (il == 0) GO TO 430
 IF (kk /= 0) GO TO 250
 kk = 2
 id3= id1
 n3 = ipar - j1
 ll = id2
 j1 = iy(ll  )
 j2 = iy(ll+1) - 1
 ll = id1
 GO TO 180
 250 ntab= ntab - 2
 id4 = id2
 n4  = ipar - j1
 GO TO 30
 
!     CASE ID4 HAS NO MORE ENTRIES. CHECK IF ID3 CAN BE PIVOT
 
 260 IF (n3 < 0) THEN
   GO TO   410
 ELSE IF (n3 == 0) THEN
   GO TO   280
 END IF
 
!     NONZERO - ID3 IS TO BE ID1
 
 270 id1 = id3
 GO TO 50
 
!     ID3 AND ID4 ARE NULL.  GO TO CLOSEST END OF TABLE FROM ID4
 
 280 i = ngrid - id4
 j = 1
 IF (i > id4) j = -1
 l = (j+2)/2 + 1
 lim = limt(l)
 LEN = id4
 
 ASSIGN 310 TO iret1
 kk = id4
 290 kk = kk + j
 IF (kk == lim) GO TO iret1, (310,420)
 ipar = 2
 ASSIGN 300 TO iret
 GO TO 370
 
!     CHECK IF ANY ENTRIES FOUND
 
 300 IF (ipar == 0) GO TO 290
 
!     ENTRY FOUND
 
 LEN = kk + j
 id4 = kk
 n4  = ipar
 GO TO 320
 
!     THAT END OF TABLE FAILED - TRY OTHER END
 
 310 j   = -j
 limt(l) = LEN
 l   = (j+2)/2 + 1
 lim = limt(l)
 kk  = id4
 ASSIGN 420 TO iret1
 GO TO 290
 
!     AN ENTRY WAS FOUND - CHECK FOR ODD NUMBER OF ENTRIES FOR PIVOT
 
 320 IF (MOD(ipar,2) == 1) GO TO 360
 
!     NOT AN ODD NUMBER OF ENTRIES FOR ID4.  CHECK GPCT ENTRIES
!     FOR ONLY ONE ENTRY.
 
 ASSIGN 340 TO iret
 
 ih = j2
 il = j1 - 1
 330 il = il + 1
 kk = iz(il)
 IF (kk <= 0) GO TO 40
 ipar = 1
 GO TO 370
 340 CONTINUE
 IF (ipar ==  1) GO TO 360
 IF (il   < ih) GO TO 330
 
!     PIVOT NOW DETERMINED
 
 GO TO 40
 360 id1 = kk
 GO TO 50
 
 
!     INTERNAL ROUTINE TO DETERMINE NUMBER OF ENTRIES FOR PIVOT
 
!     INPUT
!        IPAR = 1 -- 0,1 OR MORE THAN 1 ENTRY RETURN
!             = 2 -- ACTUAL NUMBER OF ENTRIES RETURN
!        KK   = ID OF PIVOT
 
!     OUTPUT
!        IPAR = DESIRED NUMBER OF ENTRIES
!        KK   = SAME AS INPUT
!        J1   = POINTER TO 1ST LOCATION
!        J2   = NOT NECESSARILY LAST LOCATION (I.E. IPAR INPUT AS 1)
 
 370 j1 = iy(kk  )
 j2 = iy(kk+1) - 1
 IF (ipar == 1) j2 = MIN0(j1+2,j2)
 ipar = 0
 IF (j2-j1 < 0) GO TO 390
 
 DO  i = j1,j2
   IF (iz(i) < 0) THEN
     GO TO   390
   ELSE IF (iz(i) == 0) THEN
     GO TO   400
   ELSE
     GO TO   380
   END IF
 380 ipar = ipar + 1
 END DO
 390 GO TO iret, (300,340)
 
!     ERROR MESSAGES
 
 400 ERR(1) = 1
 ERR(2) = id1
 k = 1
 GO TO 450
 410 ERR(1) = 2
 ERR(2) = n3
 ERR(3) = n4
 k = 2
 GO TO 450
 420 ERR(1) = 1
 ERR(2) = id4
 k = 3
 GO TO 450
 430 ERR(1) = 3
 ERR(2) = ll
 ERR(3) = j1
 ERR(4) = j2
 k = 4
 GO TO 450
 440 ERR(1) = 0
 k = 5
 
 450 i = lm(k)
 CALL wrtprt (merr,ERR,m(i),nm(k))
 IF (k == 5) GO TO 530
 
!     CONVERT TABLE TO ORIGINAL VALUES UNLESS THIS IS THE LAST CALL
 
 460 il = ngrid  + 1
 i  = iy(il) - 1
 IF (deform /= 0) GO TO 530
 DO  j = 1,i
   iz(j) = IABS(iz(j))
 END DO
 
 DO  j1 = 1,ngrid
   IF (iy(j1) == iy(j1+1)) CYCLE
   i = iy(j1)
   l = i
   n = iy(j1+1) - 1
   IF (i+1 > n) CYCLE
   
!     SHUTTLE EXCHANGE
!     (NOTE FROM G.CHAN/UNISYS  10/1990
!     THERE ARE MORE THAN JUST A SHUTTLE SORTING HERE. REPLACING THE
!     SHUTTLE EXCHANGE METHOD BY SORT ROUTINE, WHICH USES A MUCH FASTER
!     TECHNEQUE, DOES NOT WORK HERE)
   
   490 IF (iz(i) <= iz(i+1)) GO TO 510
   k = iz(i+1)
   iz(i+1) = iz(i)
   iz(i  ) = k
   j = i
   500 IF (j == l) GO TO 510
   IF (iz(j) >= iz(j-1)) GO TO 510
   k = iz(j)
   iz(j  ) = iz(j-1)
   iz(j-1) = k
   j = j - 1
   GO TO 500
   510 IF (i >= n-1) CYCLE
   i = i + 1
   GO TO 490
   
 END DO
 
!     A NONSTANDARD RETURN COULD BE ADDED HERE.  BAD PLOT RESULTS IF
!     THIS ROUTINE FAILS.  THE FRAME WILL BE PRESENT HOWEVER
 
 530 CONTINUE
 RETURN
END SUBROUTINE suplt
