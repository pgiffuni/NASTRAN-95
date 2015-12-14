SUBROUTINE gp3d
     
!     GP3D CREATES THE ETT (ELEMENT TEMPERATURE TABLE)
 
!     THE GPTT AS PREPARED BY GP3B COMES TO THIS ROUTINE VIA SCRATCH
!     DATA SET 1.
 
!     DATA IN THE GPTT IS USED TOGETHER WITH DATA OBTAINED FROM TEMPP1,
!     TEMPP2, TEMPP3, AND TEMPRB CARDS WHICH RESIDE ON GEOM3.
 
 LOGICAL :: anygpt   ,anyet    ,lflag    ,any      ,heat
 INTEGER :: geom3    ,eqexin   ,geom2    ,slt      ,ett      ,  &
     scr1     ,scr2     ,buf1     ,buf2     ,buf      ,  &
     FILE     ,cardid   ,carddt   ,STATUS   ,pload2   ,  &
     tempd    ,tempp1   ,tempp2   ,tempp3   ,temprb   ,  &
     rd       ,rdrew    ,wrt      ,wrtrew   ,rew      ,  &
     norew    ,z        ,flag     ,twoi     ,defalt   ,  &
     nam(2)   ,record   ,gptrec   ,setid    ,outpt    ,  &
     sysbuf   ,outwds   ,ectwds   ,elem     ,buf3
 REAL :: rz(1)    ,rbuf(50) ,tgrid(32)
 CHARACTER (LEN=23) :: ufm
 COMMON /xmssg / ufm
 COMMON /BLANK / nograv   ,noload   ,notemp
 COMMON /system/ ksystm(63)
 COMMON /names / rd       ,rdrew    ,wrt      ,wrtrew   ,rew      , norew
 COMMON /gp3com/ geom3    ,eqexin   ,geom2    ,slt      ,ett      ,  &
     scr1     ,scr2     ,buf1     ,buf2     ,buf(50)  ,  &
     cardid(60),idno(30),carddt(60),mask(60),STATUS(60)  &
     ,               ntypes   ,ipload   ,igrav    ,pload2(2),load(2)  ,  &
     nopld2   ,temp(2)  ,tempd(2) ,tempp1(2),tempp2(2), tempp3(2),temprb(2),buf3
 COMMON /gpta1 / nelem    ,last     ,incr     ,elem(1)
 COMMON /zzzzzz/ z(1)
 EQUIVALENCE     (rz(1),z(1)), (rbuf(1),buf(1)), (defalt,deftmp),  &
     (ksystm(1),sysbuf), (ksystm(2),outpt), (ksystm(56),iheat)
 DATA   nam    / 4HGP3D,4H     /
 
!                 +---------------------+
!     OPEN CORE   I                     I  Z(ILIST) = Z(1)
!                 I  ET SET-LIST        I
!     DESIGN FOR  I  2 WDS/ENTRY        I
!                 I                     I  Z(NLIST)
!     GP3D        +---------------------+
!                 I                     I  Z(IGPTT)
!                 I  GPT SET-LIST       I
!                 I  3 WDS/ENTRY        I
!                 I                     I  Z(NGPTT)
!                 +---------------------+
!                 I                     I  Z(IGPT) *
!                 I  GPTT DATA          I           *
!                 I  FOR CURRENT SETID  I            *
!                 I  2 WDS/ENTRY        I             *
!                 I                     I  Z(NGPT)     *
!                 +---------------------+               *
!                 I                     I  Z(IET1)       *
!                 I  2-DIMEN EL-TEMP    I                * THIS SPACE IS
!                 I  FOR CURRENT SETID  I                * DYNAMIC FOR
!                 I  7 WDS/ENTRY        I                * EACH SET OF
!                 I                     I  Z(NET1)       * TEMPERATURE
!                 +---------------------+                * DATA.
!                 I                     I  Z(IET2)       *
!                 I  1-DIMEN EL-TEMP    I                *
!                 I  FOR CURRENT SETID  I               *
!                 I  15 WDS/ENTRY       I              *
!                 I                     I  Z(NET2)    *
!                 +---------------------+            *
!                 I/////////////////////I           *
!                 I/////////////////////I          *
!                 +---------------------+
!                 I                     I  Z(BUF1)
!                 I  BUFFER 2           I
!                 I                     I
!                 +---------------------+
!                 I                     I  Z(BUF2)
!                 I  BUFFER 1           I
!                 I                     I  Z(KORSZ)
!                 +---------------------+
 
 
 
!     OPEN GEOM3, AND SCR1. READ IN TEMPP1, TEMPP2, TEMPP3, TEMPRB CARDS
!     CONVERT AND WRITE THEM OUT ON SCR2.
 
 heat  = .false.
 IF (iheat == 1) heat = .true.
 lflag = .false.
 j     = -1
 nwords= 8
 ilist = 1
 nlist = 0
 FILE  = geom3
 any   = .false.
 CALL preloc (*820,z(buf1),geom3)
 FILE  = scr2
 CALL OPEN (*820,scr2,z(buf2),wrtrew)
 
!     PICK UP TEMPP1 CARDS
 
 FILE = geom3
 CALL locate (*20,z(buf1),tempp1,flag)
 any  = .true.
 ASSIGN 10 TO iretrn
 buf(7) = 0
 buf(8) = 1
 10 CALL READ (*840,*20,geom3,buf,6,0,flag)
 GO TO 170
 
!     PICK UP TEMPP2 CARDS
 
 20 CALL locate (*40,z(buf1),tempp2,flag)
 any = .true.
 ASSIGN 30 TO iretrn
 30 CALL READ (*840,*40,geom3,buf,8,0,flag)
 GO TO 170
 
!     PICK UP TEMPP3 CARDS (CONVERT THESE TO LOOK LIKE TEMPP1 CARDS)
 
 40 CALL locate (*140,z(buf1),tempp3,flag)
 any = .true.
 ASSIGN 50 TO iretrn
 50 CALL READ (*840,*140,geom3,buf,24,0,flag)
 n  = 25
 DO  i = 1,11
   n  = n - 2
   IF (rbuf(n) /= 0.0 .OR. rbuf(n+1) /= 0.0) EXIT
 END DO
 70 n  = n/2
 t1 = rbuf(4)
 t2 = rbuf(2*n+2)
 IF (n == 1) GO TO 100
 h  = rbuf(2*n+1) - rbuf(3)
 sum= 0.0
 n  = n-1
 DO  i = 1,n
   twoi = 2*i
   factor = rbuf(twoi+3) - rbuf(twoi+1)
   IF (factor <= 0.0) GO TO 120
   sum  = sum + (rbuf(twoi+2) + rbuf(twoi+4))*factor
 END DO
 tbar = sum/(2.0*h)
 hover2 = h/2.0
 sum  = 0.0
 DO  i = 1,n
   twoi = 2*i
   sum  = sum + (rbuf(twoi+3) - rbuf(twoi+1)    )*(3.0*  &
       (rbuf(twoi+1) - rbuf(3) - hover2)* (rbuf(twoi+4) + rbuf(twoi+2)    ) +  &
       (rbuf(twoi+2) + 2.0*rbuf(twoi+4))* (rbuf(twoi+3) - rbuf(twoi+1)   ))
 END DO
 tprime = 2.0*sum/h**3
 GO TO 110
 
 100 tbar = rbuf(4)
 tprime = 0.0
 
 110 rbuf(3) = tbar
 rbuf(4) = tprime
 rbuf(5) = t1
 rbuf(6) = t2
 buf(7)  = 0
 buf(8)  = 1
 GO TO 170
 
!     BAD DATA ON A TEMPP3 CARD
 
 120 WRITE  (outpt,130) ufm,buf(1),buf(2)
 130 FORMAT (a23,' 4010, TEMPP3 BULK DATA CARD WITH SET ID =',i8,  &
     ' AND ELEMENT ID =',i8, /27X,  &
     'DOES NOT HAVE ASCENDING VALUES SPECIFIED FOR Z.')
 lflag = .true.
 GO TO 50
 
!     END OF 8 WORD CARDS.  WRITE EOR ON SCR2 AND DO TEMPRB CARDS NOW.
 
 140 CALL WRITE (scr2,0,0,1)
 nwords = 16
 CALL locate (*160,z(buf1),temprb,flag)
 any = .true.
 ASSIGN 150 TO iretrn
 150 CALL READ (*840,*160,geom3,buf,16,0,flag)
 GO TO 170
 
!     WRITE EOR ON SCR2. SCR2 THEN WILL HAVE 2 RECORDS (1 OR BOTH EMPTY)
 
 160 CALL WRITE (scr2,0,0,1)
 CALL CLOSE (geom3,rew )
 CALL CLOSE (scr2 ,rew )
 GO TO 230
 
!     INTERNAL SUBROUTINE TO BUILD SET LIST FROM TEMPERATURE CARD DATA
!     FIND SET-ID OR ADD IT TO LIST IN SORT, BUMP COUNT AND WRITE CARD.
 
 170 IF (j == -1) GO TO 210
 IF (buf(1) == z(j)) GO TO 180
 IF (buf(1) > z(j) .AND. j == nlist-1) GO TO 210
 
!     LOOK FOR MATCHING SETID OR FIND WHERE NEW SETID BELONGS
 
 CALL bisloc (*190,buf(1),z(ilist),2,nlist/2,j)
 
!     MATCH WAS FOUND  (ILIST ASSUMED TO BE EQUAL TO 1)
 
 180 z(j+1) = z(j+1) + nwords
 GO TO 220
 
!     ADD THIS NEW SETID INTO LIST
 
 190 IF (buf(1) > z(j)) j = j + 2
 
!     PUSH Z(J) THRU Z(NLIST) DOWN TWO WORDS TO MAKE ROOM FOR NEW SETID
 
 i = nlist + 2
 DO  k = j,nlist
   z(i) = z(i-2)
   i = i - 1
 END DO
 GO TO 211
 
!     ADD NEW SETID TO LIST
 
 210 j = j + 2
 211 z(j)  = buf(1)
 nlist = nlist + 2
 z(j+1)= nwords
 
!     WRITE OUT THE DATA CARD ON THE SCRATCH FILE FOR LATER USE
 
 220 CALL WRITE (scr2,buf,nwords,0)
 GO TO iretrn, (10,30,50,150)
 
!     READ IN GPTT HEADER RECORD FROM SCR1
 
 230 igptt = nlist + 1
 ngptt = igptt
 FILE  = scr1
 IF (notemp /= 1) GO TO 250
 CALL OPEN (*820,scr1,z(buf1),rdrew)
 CALL READ (*840,*240,scr1,z(igptt),buf2-igptt,1,flag)
 CALL mesage (-8,0,nam)
 240 ngptt = nlist + flag
 igptt = igptt + 2
 nsets = (ngptt - igptt + 1)/3
 
!     DETERMINE NUMBER OF RECORDS OF EXTERNAL INDEX TEMP DATA
!     FOLLOWING HEADER RECORD.
 
 irecs = 0
 IF (nsets > 0) THEN
   GO TO   241
 ELSE
   GO TO   247
 END IF
 241 DO  i = igptt,ngptt,3
   irecs = MAX0(z(i+2),irecs)
 END DO
 247 CONTINUE
 CALL CLOSE (scr1,norew)
 
!     OPEN ETT, PUT OUT HEADER RECORD WITH THE 3 WORD SET ENTRIES.
 
 250 IF (notemp /= 1 .AND. .NOT.any) GO TO 810
 notemp = 1
 FILE   = ett
 CALL OPEN (*820,ett,z(buf2),wrtrew)
 CALL fname (ett,buf)
 CALL WRITE (ett,buf,2,0)
 list1 = ilist
 list2 = igptt
 record= 0
 260 IF (list1 > nlist-1 .AND. list2 <= ngptt-2) GO TO 290
 IF (list1 <= nlist-1 .AND. list2 > ngptt-2) GO TO 270
 IF (list1 > nlist-1 .AND. list2 > ngptt-2) GO TO 330
 
 IF (z(list1) - z(list2) < 0.0) THEN
   GO TO   270
 ELSE IF (z(list1) - z(list2) == 0.0) THEN
   GO TO   280
 ELSE
   GO TO   290
 END IF
 
!     SET-ID OF LIST1 IS .LT. SET-ID OF LIST2 OR LIST2 IS ALL USED.
 
 270 buf(1) = z(list1)
 buf(2) = -1
 list1  = list1 + 2
 GO TO 300
 
!     SET-ID OF LIST1 IS .EQ. SET-ID OF LIST2.
 
 280 buf(1) = z(list2  )
 buf(2) = z(list2+1)
 list1  = list1 + 2
 list2  = list2 + 3
 GO TO 300
 
!     SET-ID OF LIST2 IS .LT. SET-ID OF LIST1 OR LIST1 IS ALL USED.
 
 290 buf(1) = z(list2  )
 buf(2) = z(list2+1)
 list2  = list2 + 3
 IF (z(list2-1) == 0) GO TO 310
 
!     WRITE 3-WORD SET-ID ENTRY IN HEADER
 
 300 record = record + 1
 buf(3) = record
 GO TO 320
 310 buf(3) = 0
 320 CALL WRITE (ett,buf,3,0)
 GO TO 260
 
!     HEADER RECORD IS COMPLETE.  WRITE EOR AND CLOSE WITH NOREWIND.
 
 330 CALL WRITE (ett,0,0,1)
 CALL CLOSE (ett,norew)
 
!     FOR EACH SET DEFINED IN THE EL-TEMP SET LIST AND OR THE GRID-TEMP
!     SET LIST PASS GEOM2 USING LOCATE FOR ALL THE ELEMENTS FOR
!     WHICH ETT TEMP DATA OUTPUT IS POSSIBLE.
!     IF ANY ELEMENTS CONCERNED ARE PRESENT THEN SELECT FROM THE TEMP
!     DATA AVAILABLE THAT WHICH IS APPLICABLE AND OUTPUT THE DATA ON THE
!     ETT IN THE FOLLOWING FORMAT.
 
!     CONTENTS OF 1 RECORD OF THE OUTPUT FILE ETT. EACH RECORD CONTAINS
!     DATA FOR 1 SET.
 
!         SET-ID
!         ELEMENT TYPE          * * * * * * * * * *
!         NUMBER OF TEMPERATURE DATA VALUES/EL-ID  *
!         EL-ID          *                          *
!         TEMP-VALUE      *                          *
!             .           * EL-ID                    *
!             .           * ENTRY                    *
!             .           *                          *  ELEMENT-TYPE
!         LAST-TEMP-VALUE*                           *     ENTRY
!               *             (1 OR MORE EL-ID       *
!               *              ENTRIES PER EL-TYPE   *   (1 OR MORE
!               *              ENTRY)                *    PER RECORD)
!         EL-ID          *                           *
!         TEMP-VALUE      *                          *
!             .           * EL-ID                    *
!             .           * ENTRY                    *
!             .           *                         *
!         LAST-TEMP-VALUE*                         *
!         0                     * * * * * * * * * *
 
!     IN THE ABOVE IF THE ELEMENT HAS NO SPECIAL DATA, A NEGATIVE
!     ELEMENT ID IS INSERTED FOLLOWED BY NO TEMPERATURE DATA.
 
!     NOW GATHER THE DATA AVAILABLE FOR A SET FROM SCR1 AND OR SCR2.
 
 gptrec = 1
 list1  = ilist
 list2  = igptt
 340 anygpt = .false.
 anyet  = .false.
 igpt   = 0
 ngpt   = 0
 iet1   = 0
 net1   = 0
 iet2   = 0
 net2   = 0
 IF (list1 > nlist-1) GO TO 350
 IF (list2 <= ngptt-2) GO TO 360
 GO TO 370
 350 IF (list2 <= ngptt-2) GO TO 390
 GO TO 770
 
 360 IF (z(list1) - z(list2) < 0.0) THEN
   GO TO   370
 ELSE IF (z(list1) - z(list2) == 0.0) THEN
   GO TO   380
 ELSE
   GO TO   390
 END IF
 
!     NEXT SET-ID HAS ONLY EL-TEMP DATA
 
 370 setid  = z(list1)
 defalt = -1
 anyet  = .true.
 nwords = z(list1+1)
 list1  = list1 + 2
 GO TO 400
 
!     NEXT SET-ID HAS BOTH GRID-TEMP AND EL-TEMP DATA
 
 380 setid  = z(list2  )
 defalt = z(list2+1)
 anyet  = .true.
 inrec  = z(list2+2)
 IF (inrec > 0) anygpt = .true.
 nwords = z(list1+1)
 list1  = list1 + 2
 list2  = list2 + 3
 GO TO 400
 
!     NEXT SET-ID HAS ONLY GRID-TEMP DATA
 
 390 setid  = z(list2  )
 defalt = z(list2+1)
 inrec  = z(list2+2)
 IF (inrec > 0) anygpt = .true.
 list2  = list2 + 3
 GO TO 400
 
!     AT THIS POINT READ IN ANY GRID-TEMP DATA AND/OR ANY EL-TEMP DATA.
!     SORT THE EL-TEMP DATA ON EL-ID. THE GRID-TEMP DATA IS SORTED ON
!     GRIDS
 
 400 igpt = ngptt + 1
 ngpt = igpt
 IF (.NOT.anygpt) GO TO 460
 FILE = scr1
 CALL OPEN (*820,scr1,z(buf1),rd)
 
!     POSITION GPTT TO DESIRED GRID-POINT-TEMP SET AND READ IT IN.
 
 move = inrec - gptrec
 IF (move < 0) THEN
   GO TO   410
 ELSE IF (move == 0) THEN
   GO TO   440
 ELSE
   GO TO   420
 END IF
 410 CALL REWIND (scr1)
 move = inrec
 420 DO  i = 1,move
   CALL fwdrec (*840,scr1)
 END DO
 440 gptrec = inrec + 1
 CALL READ (*840,*450,scr1,z(igpt),buf2-igpt,1,flag)
 CALL mesage (-8,0,nam)
 450 ngpt = igpt + flag - 1
 CALL CLOSE (scr1,norew)
 
!     READ IN EL-TEMP DATA PERTAINING TO THIS SET-ID
 
 460 IF (.NOT.anyet) GO TO 520
 IF (ngpt+nwords >= buf2) CALL mesage (-8,0,nam)
 FILE = scr2
 CALL OPEN (*820,scr2,z(buf1),rdrew)
 iet1 = ngpt + 1
 net1 = ngpt
 470 CALL READ (*840,*490,scr2,buf,8,0,flag)
 IF (buf(1) /= setid) GO TO 470
 DO  i = 2,8
   net1 = net1 + 1
   z(net1) = buf(i)
 END DO
 nwords  = nwords - 8
 IF (nwords /= 0) GO TO 470
 CALL fwdrec (*820,scr2)
 490 iet2 = net1 + 1
 net2 = net1
 500 CALL READ (*840,*520,scr2,buf,16,0,flag)
 IF (buf(1) /= setid) GO TO 500
 DO  i = 2,16
   net2 = net2 + 1
   z(net2) = buf(i)
 END DO
 nwords  = nwords - 16
 IF (nwords /= 0) GO TO 500
 
!     ALL DATA IS NOW IN CORE FOR THIS SET-ID
 
 520 CALL CLOSE (scr2,rew)
 IF (.NOT.anyet .AND. .NOT.anygpt) GO TO 340
 
!     SORT THE 7-WORD TEMP CARDS ON ID AND CHECK FOR DUPLICATE ID S
!     AMONG ALL THE ELEMENT TEMPERATURE DATA
 
 IF (iet1 < net1) CALL sort (0,0, 7,1,z(iet1),net1-iet1+1)
 IF (iet2 < net2) CALL sort (0,0,15,1,z(iet2),net2-iet2+1)
 
 let1 = (net1 - iet1 + 1)/7
 let2 = (net2 - iet2 + 1)/15
 lgpt = (ngpt - igpt + 1)/2
 lflag = .false.
 IF (let1 <= 1) GO TO 560
 id = z(iet1)
 j  = iet1 + 7
 DO  i = j,net1,7
   IF (id /= z(i)) GO TO 540
   
!     ERROR - TWO OR MORE ID-S EQUAL IN TEMPERATURE DATA WITHIN A SET.
   
   WRITE  (outpt,530) ufm,setid,id
   530 FORMAT (a23,' 4011, ELEMENT TEMPERATURE SET',i9,' CONTAINS ',  &
       'MULTIPLE TEMPERATURE DATA SPECIFIED FOR ELEMENT ID',i9)
   lflag = .true.
   540 id = z(i)
 END DO
 560 IF (let2 <= 1) GO TO 590
 id = z(iet2)
 j  = iet2 + 15
 DO  i = j,net2,15
   IF (id /= z(i)) GO TO 570
   WRITE (outpt,530) ufm,setid,id
   lflag = .true.
   570 id = z(i)
 END DO
 
!     OPEN GEOM2, PREPARE TO PASS GEOM2, AND OUTPUT A RECORD OF THE ETT.
 
 590 FILE = geom2
 CALL preloc (*820,z(buf1),geom2)
 
!     OPEN ETT TO PUT OUT DATA-RECORD FOR THIS SET AND WRITE SETID,
 
 FILE = ett
 CALL OPEN (*820,ett,z(buf2),wrt)
 CALL WRITE (ett,setid,1,0)
 
!     RUN THROUGH POSSIBLE TEMPERATURE DEPENDENT ELEMENTS ON GEOM2.
 
 FILE = geom2
 595 CALL ectloc (*760,FILE,buf,i)
 
!     OK DATA FOR A CARD TYPE HAS BEEN FOUND.  WRITE EL-TYPE AND
!     DATA FOR A CARD TYPE FOUND.
 
 buf(1) = elem(i+2)
 buf(2) = elem(i+14) - 1
 ieltyp = buf(1)
 
!     WRITE ELEMENT TYPE HEADER
 
 CALL WRITE (ett,buf,2,0)
 IF (elem(i+13) == 0) GO TO 740
 jtemp  = elem(i+13)
 outwds = elem(i+14)
 ectwds = elem(i+ 5)
 igrid  = elem(i+12)
 ngrid  = igrid + elem(i+9) - 1
 fgrids = 0.0
 600 CALL READ (*840,*740,geom2,buf,ectwds,0,flag)
 
!     ON FIRST PASS COUNT NUMBER OF NON-ZERO GRIDS
 
 IF (fgrids == 0.0) THEN
   GO TO   601
 ELSE
   GO TO   605
 END IF
 601 DO  j = igrid,ngrid
   IF (buf(j) /= 0) fgrids = fgrids + 1.0
 END DO
 605 CONTINUE
 
!     SELECT DATA TO BE OUTPUT
 
 IF (.NOT.anyet) GO TO 650
 SELECT CASE ( jtemp )
   CASE (    1)
     GO TO 610
   CASE (    2)
     GO TO 620
   CASE (    3)
     GO TO 650
   CASE (    4)
     GO TO 650
 END SELECT
 
!     1 - DIMENSIONAL ELEMENT-TEMP DATA MAY BE AVAIL.
 
 610 IF (let2 < 1) GO TO 650
 CALL bisloc (*650,buf(1),z(iet2),15,let2,j)
 j = iet2 + j
 
!     AVERAGE T-BAR-A AND T-BAR-B IF THIS IS A ROD, CONROD, OR TUBE
 
 IF (ieltyp /= 1 .AND. ieltyp /= 3 .AND. ieltyp /= 10) GO TO 630
 rbuf(2) = (rz(j) + rz(j+1))/2.0
 GO TO 730
 
!     2 - DIMENSIONAL ELEMENT-TEMP DATA MAY BE AVAIL.
 
 620 IF (let1 < 1) GO TO 650
 CALL bisloc (*650,buf(1),z(iet1),7,let1,j)
 j = iet1 + j
 630 DO  k = 2,outwds
   buf(k) = z(j)
   j = j + 1
 END DO
 GO TO 730
 
!     CHECK FOR GRID-POINT-TEMP-DATA
 
 650 IF (.NOT.anygpt) GO TO 700
 
!     GRID-POINT-TEMP-DATA IS AVAILABLE FOR SOME OR ALL GRID POINTS.
 
 any   = .false.
 rtemp = 0.0
 ii    = 0
 DO  k = igrid,ngrid
   ii = ii + 1
   IF (buf(k) == 0.0) THEN
     GO TO   665
   END IF
   655 CALL bisloc (*660,buf(k),z(igpt),2,lgpt,j)
   j  = igpt + j
   rtemp = rtemp + rz(j)
   IF (ii > 32) CALL mesage (-61,0,0)
   tgrid(ii) = rz(j)
   any = .true.
   CYCLE
   660 IF (defalt == -1) GO TO 710
   rtemp = rtemp + deftmp
   tgrid(ii) = deftmp
   CYCLE
   
!     UNDEFINED GRID-POINT
   
   665 tgrid(ii) = 0
 END DO
 
!     IF NOTHING BUT DEFAULT DATA THEN WRITE NOTHING SINCE THE
!     DEFAULT IS IN THE HEADER RECORD.
 
 IF (.NOT.any) GO TO 735
 
!     IF BAR ELEMENT PUT GRID TEMPS INTO BUFFER FOR T-BAR-A AND T-BAR-B
 
 IF (ieltyp /= 34) GO TO 675
 rbuf(2) = tgrid(1)
 rbuf(3) = tgrid(2)
 j = 4
 GO TO 676
 
 675 rbuf(2) = rtemp/fgrids
 j = 3
 IF (jtemp == 4) j = 2
 
 676 IF (jtemp < 3) GO TO 690
 DO  k = 1,ii
   rbuf(j) = tgrid(k)
   j = j + 1
 END DO
 690 IF (j > outwds) GO TO 730
 buf(j) = 0
 j = j + 1
 GO TO 690
 
!     NO GRID-POINT-TEMP-DATA.  VERIFY THAT THERE IS A DEFAULT TEMP.
 
 700 IF (defalt /= -1) GO TO 735
 
!     ERROR NO TEMP DATA OR DEFALT OF ANY KIND FOR THIS ID.
 
 710 lflag = .true.
 WRITE  (outpt,720) ufm,setid,buf(1)
 720 FORMAT (a23,' 4012, THERE IS NO ELEMENT, GRID POINT, OR DEFAULT',  &
     ' TEMPERATURE DATA FOR', /30X,'TEMPERATURE SET',i12,  &
     ', WITH RESPECT TO ELEMENT ID =',i8)
 GO TO 735
 
!     OUTPUT ELEMENT-TEMPERATURE DATA FOR 1 ELEMENT OF THIS TYPE IN SET
 
 730 CALL WRITE (ett,buf,outwds,0)
 GO TO 600
 
!     OUTPUT A NEGATIVE ELEMENT ID SINCE THERE IS NO DATA AVAILABLE.
 
 735 id = -buf(1)
 CALL WRITE (ett,id,1,0)
 GO TO 600
 
!     END OF ELEMENTS FOR THIS EL-TYPE.  WRITE ZERO ON ETT
 
 740 CALL WRITE (ett,0,1,0)
 GO TO 595
 760 CONTINUE
 
!     ETT-RECORD IS COMPLETE FOR THIS SET. WRITE EOR AND PROCESS NEXT
!     SET.
 
 CALL WRITE (ett,0,0,1)
 CALL CLOSE (ett,norew)
 GO TO 340
 
!     ETT IS COMPLETE
 
 770 IF (lflag) CALL mesage (-61,0,0)
 
!     WRITE TRAILER FOR ETT
 
 buf(1) = ett
 buf(7) = 7
 DO  i = 2,6
   buf(i) = 0
 END DO
 
!     OPEN ETT AND APPEND GPTT SECTION OF TEMP DATA IN INTERNAL NOTATION
 
 FILE = ett
 CALL OPEN (*820,ett,z(buf2),wrt)
 IF (.NOT.anygpt .AND. .NOT.heat) GO TO 800
 
!     OPEN SCR1 AND SKIP THE TEMPERATURE DATA HAVING EXTERNAL INDICES
 
 FILE = scr1
 CALL gopen (scr1,z(buf1),rdrew)
 IF (irecs > 0) THEN
   GO TO   780
 ELSE
   GO TO   790
 END IF
 780 DO  i = 1,irecs
   CALL fwdrec (*840,scr1)
 END DO
 
!     COPY BALANCE OF SCR1 TO ETT
 
 790 CALL READ  (*800,*795,scr1,z,buf2-1,0,flag)
 CALL WRITE (ett,z,buf2-1,0)
 GO TO 790
 795 CALL WRITE (ett,z,flag,1)
 GO TO 790
 800 CALL CLOSE (scr1,rew)
 CALL CLOSE (ett, rew)
 CALL wrttrl (buf)
 
!     THERE WAS NO GPTT DATA AND ALSO NO ETT DATA. THUS RETURN HAVING
!     CREATED NO ETT DATA SET.
 
 810 RETURN
 
!     ERROR CONDITIONS ON FILES
 
 
!     FILE NOT IN FIST OR PURGED
 
 820 j = -1
 GO TO 850
 
!     EOF HIT WHILE READING FILE
 
 840 j = -2
 850 CALL mesage (j,FILE,nam)
 RETURN
END SUBROUTINE gp3d
