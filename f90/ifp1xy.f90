SUBROUTINE ifp1xy (card,xycard)
     
!     THIS ROUTINE PROCESSES THE XRCARD IMAGES OF THE XY-PLOT CONTROL
!     CARDS AND CREATES THE -XYCDB- FILE WHICH IS OPENED AND CLOSED BY
!     THE CALLING ROUTINE.
 
!     THE ARGUMENT -CARD- IS = 1 ON THE FIRST CALL TO THIS ROUTINE
!                            = 0 ON OTHER CALLS WHEN AN IMAGE IS SENT
!                            =-1 ON LAST CALL AND NO IMAGE IS SENT
 
!     TWO RECORDS WILL BE FORMED BY THIS ROUTINE.
!     THE FIRST RECORD HAS XY-PLOT, XY-PRINT, AND XY-PUNCH DATA AND IS
!     USED BY THE XYTRAN MODULE.
!     THE SECOND RECORD IS A SORTED NX6 MATRIX STORED BY ROWS. EACH ROW
!     CONTAINS THE FOLLOWING.
 
!          1 - SUBCASE ID OR 0 INDICATING ALL.
!          2 - VECTOR CODE NUMBER E.G. DISP,STRESS,SPCF, ETC.
!          3 - POINT OR ELEMENT ID NUMBER.
!          4 - COMPONENT NUMBER.
!          5 - TYPE OF PLOT (1=RESP,2=AUTO,3=PSDF)
!          6 - DESTINATION CODE 1-7 (BIT1=PRINT,BIT2=PLOT,BIT3=PUNCH).
!              CODE 8 ADDED - BIT 4 PAPERPLOT
 
 
 INTEGER, INTENT(OUT)                     :: card
 INTEGER, INTENT(IN)                      :: xycard(1)
 EXTERNAL        lshift,rshift,andf,orf,complf
 LOGICAL :: contin,slash,pairs,ofbcd,tapbit,xycm,bit64
 CHARACTER (LEN=23) :: ufm
 INTEGER :: BLANK, bcd, buff

 COMMON /xmssg / ufm
 DIMENSION       buf(10),modid(2), kword(3),rword(14),  &
     iword(26),bword(16),icse(400),incard(20),  &
     buff(150),subcas(200),z(1),corey(771)
     
 COMMON /machin/ mach,ihalf
 COMMON /system/ ksystm(21),ilink,skp63(63),intr
 COMMON /ifp1a / dum3(3),nwpc,dummy(11),a377
 COMMON /xifp1 / BLANK,bit64
!     COMMON /ZZIFP1/ ICSE(400),INCARD(20),BUFF(150),SUBCAS(200),Z(1)
 COMMON /zzzzzz/ corex(1)
 COMMON /ifpx0 / lbd,lcc,bits(1)
 EQUIVALENCE     (ksystm( 1),sysbuf    ), (ksystm( 2),l      ),  &
     (ksystm( 3),nogo      ), (ksystm( 9),nlpp   ),  &
     (ksystm(12),line      ), (ksystm(21),iresrt ),  &
     (icse(1)   ,corex(1   ), corey(1)           ),  &
     (incard(1) ,corey(401)), (buff(1),corey(421)),  &
     (subcas(1) ,corey(571)), (z(1)   ,corey(771))
 DATA    nrword/ 14 /,  niword / 26    /, nbword  / 16 /
 DATA    ilnk  / 4HNS01        /
 DATA    kword / 4HFILM, 4HPAPE, 4HBOTH/
 DATA    rword / 4HXMIN, 4HXMAX, 4HYMIN, 4HYMAX, 4HYTMI, 4HYTMA,  &
     4HYBMI, 4HYBMA, 4HYINT, 4HXINT, 4HYTIN, 4HYBIN, 4HXPAP, 4HYPAP/
 DATA    iword / 4HXDIV, 4HYDIV, 4HYTDI, 4HYBDI, 4HXVAL, 4HYVAL,  &
     4HYTVA, 4HYBVA, 4HUPPE, 4HLOWE, 4HLEFT, 4HRIGH,  &
     4HTLEF, 4HTRIG, 4HBLEF, 4HBRIG, 4HALLE, 4HTALL,  &
     4HBALL, 4HCURV, 4HDENS, 4HCAME, 4HPENS, 4HSKIP, 4HCSCA, 4HCOLO/
 DATA    bword / 4HXAXI, 4HYAXI, 4HXTAX, 4HXBAX, 4HXLOG, 4HYLOG,  &
     4HYTLO, 4HYBLO, 4HXGRI, 4HYGRI, 4HXTGR, 4HXBGR,  &
     4HPLOT, 4HYTGR, 4HYBGR, 4HLONG/
 DATA    clea  / 4HCLEA /,   yes   / 4HYES  /,   no    / 4HNO   /,  &
     t1    / 4HT1   /,   r1    / 4HR1   /,   t1rm  / 4HT1RM /,  &
     t2    / 4HT2   /,   r2    / 4HR2   /,   t2rm  / 4HT2RM /,  &
     t3    / 4HT3   /,   r3    / 4HR3   /,   t3rm  / 4HT3RM /,  &
     t1ip  / 4HT1IP /,   r1rm  / 4HR1RM /,   r1ip  / 4HR1IP /,  &
     t2ip  / 4HT2IP /,   r2rm  / 4HR2RM /,   r2ip  / 4HR2IP /,  &
     t3ip  / 4HT3IP /,   r3rm  / 4HR3RM /,   r3ip  / 4HR3IP /,  &
     xypl  / 4HXYPL /,   xypu  / 4HXYPU /,   xypr  / 4HXYPR /,  &
     slas  / 4H/    /,   thru  / 4HTHRU /,   fram  / 4HFRAM /,  &
     xy    / 4HXY   /,   auto  / 4HAUTO /,   resp  / 4HRESP /
 DATA    psdf  / 4HPSDF /,   vdum  / 4HVDUM /,   disp  / 4HDISP /,  &
     velo  / 4HVELO /,   svel  / 4HSVEL /,   elst  / 4HELST /,  &
     acce  / 4HACCE /,   spcf  / 4HSPCF /,   sacc  / 4HSACC /,  &
     oloa  / 4HOLOA /,   load  / 4HLOAD /,   stre  / 4HSTRE /,  &
     nonl  / 4HNONL /,   subc  / 4HSUBC /,   forc  / 4HFORC /,  &
     sdis  / 4HSDIS /,   elfo  / 4HELFO /,   xtit  / 4HXTIT /,  &
     ytit  / 4HYTIT /,   ytti  / 4HYTTI /,   tcur  / 4HTCUR /,  &
     ybti  / 4HYBTI /,   xype  / 4HXYPE /,   vect  / 4HVECT /,  &
     plt1  / 4HPLT1 /,   plt2  / 4HPLT2 /,   eor   /  1     /,  &
     xypa  / 4HXYPA /,   xycm  / .false./,   noeor /  0     /
 
 DATA    iden  / 4HDENS /,   iequal/ 4H=    /,   oparen/ 4H(    /,  &
     FILE  / 4HXYCD /,   vg    / 2HVG   /,   imodel/ 4HMODE /,  &
     REAL  / -2     /,   inte  / -1     /,   cont  / 0      /,  &
     g     / 1HG    /,   f     / 1HF    /
 
 bitwrd = lbd + 1
 n = 1
 IF (intr <= 1 .AND. ilink == ilnk) GO TO 5
 incard(1) = xycard(1)
 CALL xrcard (buff,149,xycard)
 buff(150) = rshift(complf(0),1)
 a377 = buff(150)
 FILE = 301
 5 CONTINUE
 IF (card < 0.0) THEN
   GO TO   710
 ELSE IF (card == 0.0) THEN
   GO TO    20
 END IF
 
!     FIRST CALL AND FIRST CARD IMAGE.
 
 10 iat    = 0
 card   = 0
 plots  = 0
 ploter = 0
 sdrbit = 0
 vdrbit = 0
 binplt = 0
 contin = .false.
 a777   = complf(0)
 icore  = korsz(z) - 2*sysbuf - nwpc - 1
 
!     RETURNING WITH ANOTHER CARD IMAGE
 
 20 IF (buff(n) == a377) RETURN
 
 IF (.NOT.contin) GO TO 30
 contin = .false.
 GO TO icont, (370,410,430,460,520,540,570,680,640)
 
!     BEGIN PROCESSING NON-CONTINUATION CARD (MUST BEGIN WITH BCD FIELD)
 
 30 iwrd = incard(1)
 IF (bit64) CALL mvbits (BLANK,0,32,iwrd,0)
 IF (iwrd == xtit .OR. iwrd == ytit .OR. iwrd == ytti .OR.  &
     iwrd == ybti .OR. iwrd == tcur) GO TO 70
 IF (buff(n) == 0) RETURN
 
 IF (buff(n) < 0) GO TO 730
 bcd = buff(n+1)
 IF (bit64) CALL mvbits (BLANK,0,32,bcd,0)
 DO  i = 1,nrword
   IF (bcd == rword(i)) GO TO 120
 END DO
 
 DO  i = 1,niword
   IF (bcd == iword(i)) GO TO 80
 END DO
 
 DO  i = 1,nbword
   IF (bcd == bword(i)) GO TO 130
 END DO
 
 IF (bcd == clea .OR. bcd == vdum) GO TO 110
 IF (bcd == xype .OR. bcd == xypl .OR. bcd == xypr .OR.  &
     bcd == xypu .OR. bcd == xypa) GO TO 140
 GO TO 750
 
!     TITLE CARD
 
 70 CALL WRITE (FILE,incard(1),1,noeor)
 CALL WRITE (FILE,buff(1),32,noeor )
 RETURN
 
!     VERB FOLLOWED BY AN INTEGER VALUE
!     ON CAMERA CARD BCD ALSO ACCEPTED
 
 80 n = n + 2*buff(n) + 1
 IF (i == 22 .AND. buff(n) /= inte) GO TO 81
 IF (buff(n) /= inte) GO TO 770
 IF (i == 26) GO TO 95
 IF (buff(n+1) >= 0 .AND. i <= 8) GO TO 90
 IF (i <= 8) GO TO 770
 IF (buff(n+1) == 0 .OR. i > 19) GO TO 90
 buff(n+1) = buff(n+1)/IABS(buff(n+1))
 90 buf(1) = bcd
 buf(2) = buff(n+1)
 100 CALL WRITE (FILE,buf(1),2,noeor)
 RETURN
 
 95 buf(1) = bcd
 buf(2) = buff(n+1)
 buf(3) = buff(n+3)
 CALL WRITE (FILE,buf(1),3,noeor)
 RETURN
 
 110 CALL WRITE (FILE,bcd,1,noeor)
 RETURN
 
 81 iwrd = buff(n-2)
 IF (bit64) CALL mvbits (BLANK,0,32,iwrd,0)
 DO  i = 1,3
   IF (iwrd /= kword(i)) CYCLE
   buff(n+1) = i
   GO TO 90
 END DO
 GO TO 770
 
!     VERB FOLLOWED BY A REAL VALUE
 
 120 n = n + 2*buff(n) + 1
 IF (buff(n) /= REAL) GO TO 770
 GO TO 90
 
!     VERB FOLLOWED BY BCD YES OR NO, UNLESS BCD = PLOT...
 
 130 IF (i == 13) GO TO 138
 n = n + 2*buff(n) - 2
 j = n
 
!     SEARCH FOR EQUAL SIGN
 
 132 iwrd = buff(n)
 IF (bit64) CALL mvbits (BLANK,0,32,iwrd,0)
 IF (iwrd == iequal) GO TO 133
 n = n - 2
 IF (n > 0) GO TO 132
 n = j
 133 CONTINUE
 i = -1
 134 iwrd = buff(n+1)
 IF (bit64) CALL mvbits (BLANK,0,32,iwrd,0)
 IF (iwrd == yes) i = 1
 IF (iwrd == no ) i = 0
 IF (i < 0) GO TO 135
 buf(1) = bcd
 buf(2) = i
 GO TO 100
 135 IF (i < -3) GO TO 136
 i = i - 1
 n = n + 1
 GO TO 134
 136 n = j
 GO TO 750
 
!     PLOTTER SPECIFICATION CARD LOGIC
 
 138 IF (buff(n+3) == a777) n = n + 2
 n    = n + 2
 nmod = n + 3
 iwrd = buff(nmod)
 IF (bit64) CALL mvbits (BLANK,0,32,iwrd,0)
 IF (iwrd == imodel) nmod = nmod + 2
 IF (iwrd ==   iden) GO TO 147
 modid(1) = 0
 modid(2) = 0
 IF (buff(nmod) == a377) GO TO 147
 IF (buff(nmod) ==   -1) modid(1) = buff(nmod+1)
 IF (buff(nmod) /=   -1) modid(1) = buff(nmod  )
 nmod = nmod + 2
 IF (buff(nmod) ==    1) nmod = nmod + 1
 iwrd = buff(nmod)
 IF (bit64) CALL mvbits (BLANK,0,32,iwrd,0)
 IF (iwrd       == iden) GO TO 147
 IF (buff(nmod) == a377) GO TO 147
 IF (buff(nmod) ==   -1) modid(2) = buff(nmod+1)
 IF (buff(nmod) /=   -1) modid(2) = buff(nmod  )
 147 CALL fndplt (ploter,model,modid(1))
 buf(1) = bcd
 buf(2) = orf(lshift(ploter,ihalf),model+100)
 binplt = binplt + 1
 GO TO 100
 
!     PRINT, PLOT, OR PUNCH COMMAND CARD
 
 140 xtype = 0
 TYPE  = 0
 xvect = 0
 vector= 0
 PRINT = 0
 plot  = 0
 punch = 0
 paplot= 0
 slash = .false.
 n1    = 2
 n2    = 2*buff(n) + n
 
!     PROCESS ALL WORDS
 
 DO  i = n1,n2,2
   bcd = buff(i)
   IF (bcd == a777) GO TO 350
   IF (bit64) CALL mvbits (BLANK,0,32,bcd,0)
   IF (bcd == xypl) GO TO 150
   IF (bcd == xypr) GO TO 160
   IF (bcd == xypu) GO TO 170
   IF (bcd == xype) GO TO 359
   IF (bcd == xypa) GO TO 175
   IF (bcd == resp) GO TO 180
   IF (bcd == auto) GO TO 190
   IF (bcd == psdf) GO TO 200
   IF (bcd == subc) GO TO 220
   IF (bcd == disp) GO TO 230
   IF (bcd == vect) GO TO 230
   IF (bcd == velo) GO TO 240
   IF (bcd == acce) GO TO 250
   IF (bcd == spcf) GO TO 260
   IF (bcd == load) GO TO 270
   IF (bcd == stre) GO TO 280
   IF (bcd == forc) GO TO 290
   IF (bcd == sdis) GO TO 300
   IF (bcd == svel) GO TO 310
   IF (bcd == sacc) GO TO 320
   IF (bcd == nonl) GO TO 330
   IF (bcd == elfo) GO TO 290
   IF (bcd == elst) GO TO 280
   IF (bcd == oloa) GO TO 270
   IF (bcd == vg  ) GO TO 270
   n = i - 1
   GO TO 750
   150 plot  = 2
   plots = 1
   IF (ploter /= 0) GO TO 359
   ploter = 1
   model  =-1
   buf(1) = bword(13)
   buf(2) = orf(lshift(ploter,ihalf),model+100)
   binplt = binplt + 1
   CALL WRITE (FILE,buf(1),2,noeor)
   GO TO 359
   160 PRINT = 1
   GO TO 359
   170 punch = 4
   GO TO 359
   175 paplot = 1
   GO TO 359
   180 TYPE = 1
   GO TO 210
   190 TYPE = 3
   GO TO 210
   200 TYPE = 2
   GO TO 210
   210 IF (xtype /= 0) GO TO 790
   xtype = 1
   220 CYCLE
   230 vector = 1
   GO TO 291
   240 vector = 2
   GO TO 291
   250 vector = 3
   GO TO 291
   260 vector = 4
   GO TO 291
   270 vector = 5
   GO TO 291
   280 vector = 6
   GO TO 291
   290 vector = 7
   291 sdrbit = 16
   GO TO 340
   300 vector = 8
   GO TO 331
   310 vector = 9
   GO TO 331
   320 vector = 10
   GO TO 331
   330 vector = 11
   331 vdrbit = 2
   GO TO 340
   340 IF (xvect /= 0) GO TO 790
   xvect = 1
   CYCLE
   
!     DELIMETER HIT OF SOME KIND. IGNORE IF NOT LAST WORD OF BCD GROUP.
   
   350 IF (i /= n2-1) CYCLE
   iwrd = buff(i+1)
   IF (bit64) CALL mvbits (BLANK,0,32,iwrd,0)
   IF (iwrd == slas) slash = .true.
   IF (.NOT.slash) GO TO 810
   CYCLE
   359 xycm = .true.
 END DO
 
!     WRITE PLOT CONTROL INFORMATION
 
 buf(1) = xy
 buf(2) = PRINT
 buf(3) = plot
 buf(4) = punch
 IF (paplot == 1) plot = 2
 destin = PRINT + plot + punch
 IF (TYPE == 0) TYPE = 1
 buf(5) = TYPE
 IF (vector == 0) GO TO 1030
 buf(6) = vector
 buf(7) = paplot
 CALL WRITE (FILE,buf(1),7,noeor)
 
!     ALL WORDS PROCESSED. IF SLASH HAS NOT BEEN HIT, START READING
!     SUBCASE NUMBERS.
 
 nsubs = 0
 n = n2 + 1
 IF (slash) GO TO 490
 
!     FORM LIST OF SUBCASES, MAXIMUM OF 200 FOR THIS COMMAND CARD.
 
 370 IF (buff(n) /= cont) GO TO 380
 ASSIGN 370 TO icont
 GO TO 700
 
 380 subcas(1) = 0
 IF (buff(n) /= inte) GO TO 830
 
!     SUBCASES ARE NOT APPLICABLE IN AUTO AND PSDF
 
 IF (TYPE /= 1) GO TO 850
 
 390 nsubs = nsubs + 1
 IF (nsubs   > 200) GO TO 890
 IF (buff(n+1) <= 0) GO TO 910
 subcas(nsubs) = buff(n+1)
 400 n = n + 2
 410 IF (buff(n) == a377) GO TO 1080
 IF (buff(n) /= cont) GO TO 420
 ASSIGN 410 TO icont
 GO TO 700
 
 420 IF (buff(n) /= inte) GO TO 430
 IF (subcas(nsubs)-buff(n+1) < 0.0) THEN
   GO TO   390
 ELSE IF (subcas(nsubs)-buff(n+1) == 0.0) THEN
   GO TO   400
 ELSE
   GO TO   910
 END IF
 430 IF (buff(n) /= cont) GO TO 440
 ASSIGN 430 TO icont
 GO TO 700
 
 440 IF (buff(n) < 0) GO TO 830
 iwrd = buff(n+2)
 IF (bit64) CALL mvbits (BLANK,0,32,iwrd,0)
 IF (iwrd /= slas) GO TO 450
 slash = .true.
 n = n + 3
 GO TO 490
 
 450 IF (iwrd /= thru) GO TO 830
 n = n + 3
 460 IF (buff(n) /= cont) GO TO 470
 ASSIGN 460 TO icont
 GO TO 700
 
 470 IF (buff(n  ) /=         inte ) GO TO 830
 IF (buff(n+1) < subcas(nsubs)) GO TO 910
 IF (buff(n+1) == subcas(nsubs)) GO TO 400
 480 nsubs = nsubs + 1
 IF (nsubs > 200) GO TO 890
 subcas(nsubs) = subcas(nsubs-1) + 1
 IF (subcas(nsubs) < buff(n+1)) GO TO 480
 GO TO 400
 
!     SLASH HIT. BEGIN PROCESSING FRAME DATA. FIRST WRITE SUBCASE
!     NUMBERS.
 
 490 CALL WRITE (FILE,nsubs,1,noeor)
 IF (nsubs /= 0) CALL WRITE (FILE,subcas(1),nsubs,noeor)
 IF (nsubs == 0) subcas(1) = 0
 IF (nsubs == 0) nsubs = 1
 500 slash  = .false.
 CALL WRITE (FILE,fram,1,noeor)
 pairs  = .false.
 ncurve = 0
 520 IF (buff(n) /= cont) GO TO 530
 ASSIGN 520 TO icont
 GO TO 700
 
 530 IF (buff(n) /= inte) GO TO 830
 buf(1) = buff(n+1)
 buf(2) = 0
 buf(3) = 0
 idcom  = 0
 ncurve = ncurve + 1
 IF (buf(1) <= 0) GO TO 930
 
!     GET COMPONENT. POSITIVE INTEGER.
!     MAY BE T1,T2,T3,R1,R2,R3 ETC. IF THE VECTOR IS NOT STRESS OR FORCE
 
 n = n + 2
 540 IF (buff(n) /= cont) GO TO 550
 ASSIGN 540 TO icont
 GO TO 700
 
 550 iwrd = buff(n+2)
 IF (bit64) CALL mvbits (BLANK,0,32,iwrd,0)
 IF (buff(n) > 0 .AND. iwrd == oparen) GO TO 560
 GO TO 830
 
 560 ofbcd = .false.
 IF (buff(n) > 1) GO TO 590
 
!     FALL HERE AND A POSITIVE INTEGER COMPONENT IS EXPECTED.
 
 n = n + 3
 570 IF (buff(n) /= cont) GO TO 580
 ASSIGN 570 TO icont
 GO TO 700
 
 580 IF (buff(n  ) /= inte) GO TO 830
 IF (buff(n+1) <=    0) GO TO 950
 ofbcd  =.false.
 compon = buff(n+1)
 n = n + 2
 GO TO 620
 
!     FALL HERE AND A BCD COMPONENT IS EXPECTED. T1,T2,T3,R1,R2,R3
 
 590 n1  = n + 3
 600 n   = n + 2*buff(n) + 1
 610 bcd = buff(n1)
 IF (bit64) CALL mvbits (BLANK,0,32,bcd,0)
 IF (bcd == BLANK) GO TO 615
 IF (vector == 6 .OR. vector == 7) GO TO 970
 615  ofbcd  = .true.
 compon = 3
 IF (bcd == t1 .OR. bcd == t1rm) GO TO 620
 IF (bcd == g) GO TO 620
 compon = 4
 IF (bcd == t2 .OR. bcd == t2rm) GO TO 620
 IF (bcd == f) GO TO 620
 compon = 5
 IF (bcd == t3 .OR. bcd == t3rm) GO TO 620
 compon = 6
 IF (bcd == r1 .OR. bcd == r1rm) GO TO 620
 compon = 7
 IF (bcd == r2 .OR. bcd == r2rm) GO TO 620
 compon = 8
 IF (bcd == r3 .OR. bcd == r3rm) GO TO 620
 compon = 9
 IF (bcd == t1ip) GO TO 620
 compon = 10
 IF (bcd == t2ip) GO TO 620
 compon = 11
 IF (bcd == t3ip) GO TO 620
 compon = 12
 IF (bcd == r1ip) GO TO 620
 compon = 13
 IF (bcd == r2ip) GO TO 620
 compon = 14
 IF (bcd == r3ip) GO TO 620
 compon = 1000
 IF (bcd == BLANK) GO TO 620
 GO TO 990
 
 620 idcom = idcom + 1
 buf(idcom+1) = compon
 
!     CHECK RANGE OF COMPONENT
 
 IF (compon == 1000) GO TO 631
 IF ((TYPE == 2 .OR. TYPE == 3) .AND. (compon < 3 .OR.compon > 8)  &
     .AND. (vector /= 6 .AND. vector /= 7)) GO TO 1130
 IF ((compon < 3 .OR. compon > 14) .AND.  &
     (vector /= 6 .AND. vector /= 7)) GO TO 1150
 IF (nogo /= 0) GO TO 631
 
!     ADD THIS COMPONENT-ID TO XY-MASTER SET IN OPEN CORE.
 
 DO  i = 1,nsubs
   IF (iat+6 > icore) GO TO 1090
   z(iat+1) = subcas(i)
   z(iat+2) = vector
   z(iat+3) = buf(1)
   z(iat+4) = compon
   z(iat+5) = TYPE
   z(iat+6) = destin
   iat = iat + 6
 END DO
 
!     PROCEED TO NEXT COMPONENT OR ID OF THIS FRAME
 
 631 IF (ncurve == 1 .AND. idcom == 2) pairs = .true.
 IF (pairs .AND. (TYPE == 2 .OR. TYPE == 3)) GO TO 1110
 IF (.NOT.pairs .AND. idcom == 2) GO TO 1050
 IF (idcom > 2) GO TO 1050
 IF (.NOT.ofbcd  ) GO TO 640
 IF (n1 >=  n-2) GO TO 640
 n1   = n1 + 2
 iwrd = buff(n1+1)
 IF (bit64) CALL mvbits (BLANK,0,32,iwrd,0)
 IF (iwrd /= slas) GO TO 610
 slash = .true.
 GO TO 670
 
!     IS NEXT FIELD AN INTEGER FOLLOWED BY AN OPAREN
 
 640 IF (buff(n) /= cont) GO TO 650
 ASSIGN 640 TO icont
 GO TO 700
 
 650 IF (buff(n  ) /= inte) GO TO 660
 IF (buff(n+2) == a377) GO TO 580
 iwrd = buff(n+4)
 IF (bit64) CALL mvbits (BLANK,0,32,iwrd,0)
 IF (iwrd == oparen) GO TO 670
 GO TO 580
 660 IF (buff(n) == a377) GO TO 670
 iwrd = buff(n+2)
 IF (bit64) CALL mvbits (BLANK,0,32,iwrd,0)
 IF (buff(n) <= 0 .OR. iwrd == slas) GO TO 670
 n1 = n + 1
 GO TO 600
 670 IF (pairs .AND. idcom == 1) GO TO 1050
 CALL WRITE (FILE,buf(1),3,noeor)
 IF (.NOT.slash .AND. buff(n) == inte) GO TO 520
 buf(1) = -1
 buf(2) = -1
 buf(3) = -1
 CALL WRITE (FILE,buf(1),3,noeor)
 IF (buff(n) == a377) RETURN
 
 680 IF (buff(n) /= cont) GO TO 690
 ASSIGN 680 TO icont
 GO TO 700
 690 IF (slash) GO TO 500
 iwrd = buff(n+2)
 IF (bit64) CALL mvbits (BLANK,0,32,iwrd,0)
 IF (iwrd /= slas) GO TO 830
 n = n + 2*buff(n) + 1
 GO TO 500
 
!     RETURN FOR A CONTINUATION CARD
 
 700 contin = .true.
 RETURN
 
!     NO MORE CARDS AVAILABLE. RAP IT UP IF NO ERROR. WRITE XY-SET
!     RECORD
 
 710 IF (contin) GO TO 1010
 CALL WRITE (FILE,z(1),0,eor)
 IF (iat == 0) GO TO 720
 j = 7
 DO  i = 1,6
   j = j - 1
   CALL sort (0,0,6,-j,z(1),iat)
 END DO
 720 CALL WRITE (FILE,z(1),iat,eor)
 
!     SET CARD = 0 IF NO PLOTS
!     SET CARD = 1 IF PLOTS
 
 card = plots
 
!     SET RESTART BITS FOR VDR AND SDR
 
 IF (iresrt < 0) bits(bitwrd) = orf(bits(bitwrd),vdrbit+sdrbit)
 
!     CHECK FOR COMMAND OP CARD
 
 IF (.NOT.xycm) GO TO 1170
 
!     CHECK PLOT TAPE BITS
 
 IF (plots == 0) RETURN
 
!     CHECK FOR TAPE SETUPS
 
 IF (binplt /= 0 .AND. .NOT.tapbit(plt1) .AND. .NOT.tapbit(plt2))  &
     CALL ifp1d (-618)
 RETURN
 
!     FATAL ERROR CONDITIONS
 
 730 j = 675
 WRITE  (l,740) ufm,j
 740 FORMAT (a23,i4,', ABOVE CARD DOES NOT BEGIN WITH A NON-NUMERIC ',  &
     'WORD.')
 GO TO 2000
 750 j = 676
 WRITE  (l,760) ufm,j,buff(n+1),buff(n+2)
 760 FORMAT (a23,i4,1H,,2A4,' IS NOT RECOGNIZED AS AN XYPLOT COMMAND ',  &
     'CARD OR PARAMETER.')
 GO TO 2000
 770 j = 677
 WRITE  (l,780) ufm,j
 780 FORMAT (a23,i4,', ILLEGAL VALUE SPECIFIED.')
 GO TO 2000
 790 j = 678
 WRITE  (l,800) ufm,j,buff(i),buff(i+1)
 800 FORMAT (a23,i4,1H,,2A4,' CONTRADICTS PREVIOUS DEFINITION.')
 GO TO 2000
 810 j = 679
 WRITE  (l,820) ufm,j,buff(i+1)
 820 FORMAT (a23,i4,1H,,a4,' DELIMITER ILLEGALLY USED.')
 GO TO 2000
 830 IF (buff(n) == REAL) GO TO 850
 IF (buff(n) == inte) GO TO 870
 j = 680
 WRITE  (l,840) ufm,j,buff(n+1),buff(n+2)
 840 FORMAT (a23,i4,1H,,2A4,' IS ILLEGAL IN STATEMENT.')
 GO TO 2000
 850 j = 681
 WRITE  (l,860) ufm,j,buff(n+1)
 860 FORMAT (a23,i4,1H,,e16.8,' IS ILLEGAL IN STATEMENT.')
 GO TO 2000
 870 j = 682
 WRITE  (l,880) ufm,j,buff(n+1)
 880 FORMAT (a23,i4,1H,,i10,' IS ILLEGAL IN STATEMENT.')
 GO TO 2000
 890 j = 683
 WRITE  (l,900) ufm,j
 900 FORMAT (a23,i4,', TOO MANY SUBCASES. MAXIMUM = 200 ON ANY ONE XY',  &
     '-OUTPUT COMMAND CARD.')
 GO TO 2000
 910 j = 684
 WRITE  (l,920) ufm,j
 920 FORMAT (a23,i4,', SUBCASE-ID IS LESS THAN 1 OR IS NOT IN ',  &
     'ASCENDING ORDER.')
 GO TO 2000
 930 j = 685
 WRITE  (l,940) ufm,j,buf(1)
 940 FORMAT (a23,i4,1H,,i12,' = POINT OR ELEMENT ID IS ILLEGAL (LESS ',  &
     'THAN 1).')
 GO TO 2000
 950 j = 686
 WRITE  (l,960) ufm,j
 960 FORMAT (a23,i4,', NEGATIVE OR ZERO COMPONENTS ARE ILLEGAL.')
 GO TO 2000
 970 j = 687
 WRITE  (l,980) ufm,j
 980 FORMAT (a23,i4,', ALPHA-COMPONENTS ARE NOT PERMITTED FOR STRESS ',  &
     'OR FORCE XY-OUTPUT REQUESTS.')
 GO TO 2000
 990 j = 688
 WRITE  (l,1000) ufm,j,bcd
 1000 FORMAT (a23,i4,1H,,a4,' COMPONENT NAME NOT RECOGNIZED.')
 GO TO 2000
 1010 j = 689
 WRITE  (l,1020) ufm,j
 1020 FORMAT (a23,i4,', LAST CARD ENDED WITH A DELIMITER BUT NO ',  &
     'CONTINUATION CARD WAS PRESENT.')
 GO TO 2000
 1030 j = 690
 WRITE  (l,1040) ufm,j
 1040 FORMAT (a23,i4,', TYPE OF CURVE WAS NOT SPECIFIED. (E.G. ',  &
     'DISPLACEMENT, STRESS, ETC.).')
 GO TO 2000
 1050 j = 691
 WRITE  (l,1060) ufm,j
 1060 FORMAT (a23,i4,', MORE THAN 2 OR UNEQUAL NUMBER OF COMPONENTS ',  &
     'FOR ID-S WITHIN A SINGLE FRAME.')
 GO TO 2000
 1070 FORMAT (a23,i4,', XY-OUTPUT COMMAND IS INCOMPLETE.')
 1080 j = 692
 WRITE  (l,1070) ufm,j
 GO TO 2000
 1090 j = 693
 WRITE  (l,1100) ufm,j
 1100 FORMAT (a23,i4,', INSUFFICIENT CORE FOR SET TABLE.')
 icrq = (nsubs-i+1) * 6
 WRITE  (l,1101) icrq
 1101 FORMAT (5X,8HAT least,i8,19H more words needed.)
 GO TO 2000
 1110 j = 694
 WRITE  (l,1120) ufm,j
 1120 FORMAT (a23,i4,', AUTO OR PSDF REQUESTS MAY NOT USE SPLIT FRAME',  &
     ', THUS ONLY ONE COMPONENT PER ID IS PERMITTED.')
 GO TO 2000
 1130 j = 695
 WRITE  (l,1140) ufm,j,compon
 1140 FORMAT (a23,i4,', COMPONENT VALUE =',i8,', IS ILLEGAL FOR AUTO ',  &
     'OR PSDF VECTOR REQUESTS.')
 GO TO 2000
 1150 j = 696
 WRITE  (l,1160) ufm,j,compon
 1160 FORMAT (a23,i4,', COMPONENT VALUE =',i8,', IS ILLEGAL FOR VECTOR',  &
     ' TYPE SPECIFIED.')
 GO TO 2000
 1170 j = 697
 WRITE  (l,1180) ufm,j
 1180 FORMAT (a23,i4,', XYPLOT, XYPRINT, XYPUNCH, XYPEAK, OR XYPAPLOT',  &
     /5X,' COMMAND CARD NOT FOUND IN XY PLOTTER OUT PUT PACKAGE.')
 GO TO 2000
 2000 nogo = 1
 line = line + 2
 IF (line >= nlpp) CALL page
 RETURN
END SUBROUTINE ifp1xy
