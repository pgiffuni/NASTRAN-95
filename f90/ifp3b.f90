SUBROUTINE ifp3b
     
!        CARDS           TYPE         REC.ID-BIT  CARDS-FILE, CARDS-FILE
!    === =======         ===========  ==========  ==========  ==========
!     1  AXIC     -----  AX.SY.SHELL     515- 5
!     2  CCONEAX  -----  AX.SY.SHELL    8515-85  CCONE-GEOM2,
!     3  FORCEAX  -----  AX.SY.SHELL    2115-21  FORCE-GEOM3,
!     4  FORCE    -----  STANDARD       4201-42  FORCE-GEOM3,
!     5  GRAV     -----  STANDARD       4401-44   GRAV-GEOM3,
!     6  LOAD     -----  STANDARD       4551-61   LOAD-GEOM3,
!     7  MOMAX    -----  AX.SY.SHELL    3815-38  MOMNT-GEOM3,
!     8  MOMENT   -----  STANDARD       4801-48  MOMNT-GEOM3,
!     9  MPCADD   -----  STANDARD       4891-60 MPCADD-GEOM4,
!    10  MPCAX    -----  AX.SY.SHELL    4015-40    MPC-GEOM4,
!    11  OMITAX   -----  AX.SY.SHELL    4315-43   OMIT-GEOM4,
!    12  POINTAX  -----  AX.SY.SHELL    4915-49    MPC-GEOM4, GRID-GEOM1
!    13  PRESAX   -----  AX.SY.SHELL    5215-52  PLOAD-GEOM3,
!    13+ RFORCE   -----  STANDARD       5509-55 RFORCE-GEOM3,
!    14  RINGAX   -----  AX.SY.SHELL    5615-56    SPC-GEOM4, GRID-GEOM1
!    15  SECTAX   -----  AX.SY.SHELL    6315-63    MPC-GEOM4, GRID-GEOM1
!    16  SEQGP    -----  STANDARD       5301-53  SEQGP-GEOM1,
!    17  SPCADD   -----  STANDARD       5491-59 SPCADD-GEOM4,
!    18  SPCAX    -----  AX.SY.SHELL    6215-62    SPC-GEOM4,
!    19  SUPAX    -----  AX.SY.SHELL    6415-64 SUPORT-GEOM4,
!    20  TEMPAX   -----  AX.SY.SHELL    6815-68   TEMP-GEOM3,
!    21  TEMPD    -----  STANDARD       5641-65  TEMPD-GEOM3,
 
 IMPLICIT INTEGER (a-z)
 EXTERNAL        lshift    ,andf      ,orf
 LOGICAL :: nogo      ,recoff    ,ifpdco
 REAL :: nphi      ,nphi1     ,nisq      ,ni        ,  &
     raddeg    ,rz        ,t1        ,t2        ,  &
     t3        ,coef      ,a1        ,a2        ,  &
     a3        ,a4        ,angle     ,gc        , sum       ,consts
 DIMENSION       geom(4)   ,z(13)
 DIMENSION       isystm(175)
 CHARACTER (LEN=25) :: sfm
 CHARACTER (LEN=29) :: uim
 CHARACTER (LEN=25) :: uwm
 CHARACTER (LEN=23) :: ufm
 COMMON /xmssg / ufm       ,uwm       ,uim       ,sfm
 COMMON /BLANK / bottom
 COMMON /system/ ibufsz    ,nout      ,noflag    ,dumdum(8) ,  &
     nlines    ,ddd(14)   ,mn        ,dum50(50) , ipiez
 COMMON /machin/ mach      ,ihalf
 COMMON /two   / two(32)
 COMMON /condas/ consts(5)
 COMMON /ifp3lv/            recid(3)  ,recid1(3) ,recidx(3) ,  &
     iend      ,REC(3)    ,rec1(3)   ,trail(7)  ,  &
     it        ,axtrl(7)  ,openfl(6) ,n         ,  &
     a1        ,csid      ,ni        ,nisq      ,  &
     a2        ,ibuff1    ,ibuff2    ,ibuff3    ,  &
     a3        ,buff      ,nogo      ,op        ,  &
     a4        ,iheadr    ,ibitr     ,ifile     ,  &
     noreg     ,last      ,ierrtn    ,icont     ,  &
     noaxic    ,ringid    ,outbuf    ,veor      ,  &
     istart    ,iretrn    ,flag      ,iamt      ,  &
     sum       ,ibit      ,setid     ,sorc      ,  &
     ibegin    ,mpcon     ,nwords    ,nnn       ,  &
     angle     ,k3or6     ,nphi1     ,zpt       ,  &
     nmove     ,csset     ,nopont    ,non       ,  &
     iphi      ,recoff    ,nphi      ,n3or5     ,  &
     ion       ,nplus1    ,nosect    ,coef      ,  &
     ipt       ,compon    ,icore     ,iscrat    ,  &
     icore1    ,ncards    ,i1        ,iat       ,  &
     i2        ,t1        ,t2        ,nfile     , nadd      ,ncard
 COMMON /ifp3cm/ FILE(6)   ,iname(12) ,cdtype(50),axic1(3)  ,  &
     cconex(3) ,forcex(3) ,force(3)  ,grav(3)   ,  &
     load(3)   ,momax(3)  ,moment(3) ,mpcadd(3) ,  &
     mpcax(3)  ,omitax(3) ,pointx(3) ,presax(3) ,  &
     ringax(3) ,sectax(3) ,seqgp(3)  ,spcax(3)  ,  &
     supax(3)  ,tempax(3) ,tempd(3)  ,pload(3)  ,  &
     mpc(3)    ,spc(3)    ,grid(3)   ,suport(3) ,  &
     neg111(3) ,t65535(3) ,temp(3)   ,omit(3)   ,  &
     spcadd(3) ,one       ,zero      ,iheadb(96), ctriaa(3) ,ctrapa(3) ,iconso
 COMMON /output/ dummy(96) ,ihead(96)
 COMMON /ifpdta/ dum(521)  ,gc(7)     ,ll(6)
 COMMON /zzzzzz/ rz(1)
 EQUIVALENCE     (consts(4),raddeg )  , (z(1)    ,rz(1)  )  ,  &
     (geom(1)  ,FILE(1))  , (scrtch  ,FILE(5))  ,  &
     (axic     ,FILE(6))                        ,  &
     (noeor    ,inprwd    ,  zero            )  ,  &
     (eor      ,clorwd    ,  outrwd  ,one    )
 EQUIVALENCE     (ibufsz   ,isystm(1))
 DATA    ifist / 4HFIST/   ,i3,i4,i5  /  3,4,5   /
 
 
!     GEOM4 PROCESSING
!     ================
 
!     OPEN GEOM4
 
 ifile= geom(4)
 i    = 4
 op   = outrwd
 buff = ibuff2
 ASSIGN 20 TO iretrn
 GO TO 1340
 
!     SPCADD OR MPCADD CARDS
!     ======================
 
 20 ASSIGN 30 TO icont
 REC(1)  = mpcadd(1)
 REC(2)  = mpcadd(2)
 REC(3)  = mpcadd(3)
 rec1(1) = mpcax(1)
 rec1(2) = mpcax(2)
 rec1(3) = mpcax(3)
 21 ASSIGN 28 TO iheadr
 GO TO 1470
 
!     MANDATORY SPCADD AND MPCADD CARDS.
 
 28 z( 1) = 100000101
 z( 2) = 101
 z( 3) = -1
 z( 4) = 200000102
 z( 5) = 102
 z( 6) = -1
 z( 7) = 100000000
 z( 8) = 101
 z( 9) = -1
 z(10) = 200000000
 z(11) = 102
 z(12) = -1
 IF (nogo) GO TO 22
 CALL WRITE (geom(4),z(1),12,noeor)
 22 CALL locate (*23,z(ibuff1),REC(1),flag)
 
!     READ AN OPEN ENDED SPCADD OR MPCADD CARD INTO CORE.
 
 i = 1
 27 CALL READ (*1540,*23,axic,z(i),1,noeor,iamt)
 IF (z(i) < 0.0) THEN
   GO TO    25
 END IF
 24 i = i + 1
 IF ((i+1) > icore) GO TO 1580
 GO TO 27
 
!     COMPLETE CARD IS AT HAND
 
 25 z(i) = 101
 i    = i + 1
 z(i) = -1
 z(1) = z(1) + 100000000
 IF (nogo) GO TO icont, (30,610)
 26 CALL WRITE (geom(4),z(1),i,noeor)
 IF (z(i-1) == 102) GO TO 23
 z(i-1) = 102
 z(1  ) = z(1) + 100000000
 GO TO 26
 
!     ALL SPCADD OR MPCADD CARDS COMPLETE.
!     NOW CREATE SPCADD OR MPCADD FROM SPCAX OR MPCAX
!     CARDS RESPECTIVELY.
 
 23 irec   = REC(1)
 REC(1) = rec1(1)
 REC(2) = rec1(2)
 REC(3) = rec1(3)
 CALL locate (*35,z(ibuff1),REC(1),flag)
 
!     OK SPCAX OR MPCAX RECORD EXISTS.
 
 ilast = -1
 38 z(4)  = -1
 CALL READ (*1540,*35,axic,z(2),1,noeor,iamt)
 
!     MPCAX CARDS ARE OPEN ENDED
!     SPCAX CARDS ARE 5 WORDS LONG.
 
 IF (z(2) == ilast) GO TO 47
 ilast = z(2)
 
!     CREATE TWO SPCADD OR MPCADD CARDS.
 
 z(3) = 101
 z(1) = z(2) + 100000000
 IF (nogo) GO TO icont, (30,610)
 33 CALL WRITE (geom(4),z(1),4,noeor)
 IF (z(3) == 102) GO TO 47
 z(3) = 102
 z(1) = z(1) + 100000000
 GO TO 33
 
!     READ UP TO NEXT CARD
 
 47 CALL READ (*1540,*35,axic,z(1),4,noeor,iamt)
 IF (REC(1) == spcax(1) .OR. z(1) == (-1)) GO TO 38
 GO TO 47
 
!     ALL CARDS COMPLETE.
!     WRITE EOR AND PUT BITS IN TRAILER.
 
 35 iamt = 0
 ASSIGN 37 TO iretrn
 IF (irec == spcadd(1)) GO TO 39
 REC(1) = mpcadd(1)
 REC(2) = mpcadd(2)
 REC(3) = mpcadd(3)
 GO TO 1300
 39 REC(1) = spcadd(1)
 REC(2) = spcadd(2)
 REC(3) = spcadd(3)
 GO TO 1300
 37 GO TO icont, (30,610)
 
!     MPCAX CARD
!     ==========
 
 30 mpcon  = 0
 REC(1) = mpc(1)
 REC(2) = mpc(2)
 REC(3) = mpc(3)
 recoff = .false.
 last   = -1
 ncard  = 10
 nwords = 0
 CALL locate (*130,z(ibuff1),mpcax(1),flag)
 
!     WRITE RECORD HEADER
 
 recoff = .true.
 ASSIGN 40 TO iheadr
 GO TO 1470
 
 40 mpcon = 1
 last  = 0
 
!     READ SET ID
 
 50 CALL READ (*1540,*120,axic,setid,1,noeor,iamt)
 IF (setid > 100) GO TO 130
 nwords = nwords + 1
 IF (nogo) GO TO 60
 CALL WRITE (geom(4),setid,1,noeor)
 
!     READ 4-WORDS SETS UNTIL -1,-1,-1,-1 ENCOUNTERED...
 
 60 CALL READ (*1540,*100,axic,z(1),4,noeor,iamt)
 nwords = nwords + 4
 IF (z(4) == -1) GO TO 90
 
!     CHECK HARMONIC NUMBER
 
 nnn = z(2)
 ASSIGN 70 TO ierrtn
 GO TO 1420
 
!     CHECK RING ID
 
 70 nnn = z(1)
 ASSIGN 80 TO ierrtn
 GO TO 1440
 
 80 z(2) = z(1) + (z(2)+1)*1000000
 IF (nogo) GO TO 60
 CALL WRITE (geom(4),z(2),3,noeor)
 GO TO 60
 
!     END OF EQUATION
 
 90 IF (nogo) GO TO 50
 CALL WRITE (geom(4),neg111(1),3,noeor)
 GO TO 50
 100 CALL page2 (3)
 imsg = 1063
 WRITE  (nout,105) sfm,imsg
 105 FORMAT (a25,i5)
 WRITE  (nout,110) sfm,imsg
 110 FORMAT (5X,50HEOR on axic FILE WHILE reading mpcax card records.)
 nogo = .true.
 GO TO 1530
 120 last = 1
 
!     FIRST NWORDS HAVE BEEN PROCESSED OF MPCAX CARDS UNLESS
!     LAST = 1, IN WHICH CASE ALL MPCAX CARDS ARE COMPLETE.
!     GO NOW TO THE S-SET MPC CARD-GENERATION FOR POINTAX CARDS
!     IF LAST = -1, THERE ARE NO MPCAX CARDS.
 
 
!     S-SET MPC-S FROM POINTAX CARDS
!     ==============================
 
 130 REC(1) = pointx(1)
 REC(2) = pointx(2)
 REC(3) = pointx(3)
 ncard  = 12
 n3or5  = 3
 k3or6  = 6
 sorc   = 101
 ASSIGN 380 TO icont
!     TURN NOPONT OR NOSECT ON IF POINTAX OR SECTAX CARDS EXIST RESPECT.
 
 ibit = pointx(2)
 ASSIGN 140 TO ibitr
 GO TO 1460
 140 nopont = non
 ibit = sectax(2)
 ASSIGN 150 TO ibitr
 GO TO 1460
 150 nosect = non
 
 IF (nopont == 0) THEN
   GO TO   370
 END IF
 
 160 CALL locate (*370,z(ibuff1),REC(1),flag)
 mpcon = 1
 IF (recoff) GO TO 170
 
!     WRITE RECORD HEADER
 
 recoff = .true.
 REC(1) = mpc(1)
 REC(2) = mpc(2)
 REC(3) = mpc(3)
 ASSIGN 170 TO iheadr
 GO TO 1470
 
 170 CALL READ (*1540,*370,axic,z(1),n3or5,noeor,iamt)
 
!     CHECK RING ID FOR S-SET PASS ONLY FOR POINTAX AND SECTAX CARDS.
!     NO CHECK WILL BE MADE IN THE GRID CARD GENERATION AREA.
 
!     IF (SORC .EQ. 102) GO TO 785
 nnn = z(2)
 ASSIGN 180 TO ierrtn
 GO TO 1440
 
 180 iat = n3or5 + 1
 DO  i = 1,k3or6
   z(iat) = sorc
   z(iat+1) = z(1)
   z(iat+2) = i
   rz(iat+3)= -1.0
   IF (nogo) GO TO 190
   CALL WRITE (geom(4),z(iat),4,noeor)
   190 DO  j = 1,nplus1
     
!     COMPUTE COEFFICIENT.
     
     ni = j - 1
     IF (n3or5 == 5) GO TO 240
     
!     POINTAX CARD COEFFICIENTS
     
     t1 = ni*rz(i3)*raddeg
     IF (andf(i,1) > 0.0) THEN
       GO TO   200
     ELSE
       GO TO   210
     END IF
     
!     ODD I
     
     200 IF (sorc - 101 > 0.0) THEN
       GO TO   230
     ELSE
       GO TO   220
     END IF
     
!     EVEN I
     
     210 IF (sorc - 101 > 0.0) THEN
       GO TO   220
     ELSE
       GO TO   230
     END IF
     
     220 coef = SIN(t1)
     GO TO 340
     
     230 coef = COS(t1)
     IF (sorc == 101) coef = -coef
     IF (ni == 0.0 .AND. sorc == 101) coef = 1.0
     GO TO 340
     
!     SECTAX CARD COEFFICIENTS
     
     240 t1 = ni*rz(i4)*raddeg
     t2 = ni*rz(i5)*raddeg
     IF (i >=  4) GO TO 245
     IF (andf(i,1) > 0.0) THEN
       GO TO   280
     ELSE
       GO TO   250
     END IF
     245 IF (andf(i,1) > 0.0) THEN
       GO TO   250
     ELSE
       GO TO   280
     END IF
     
!     EVEN I
     
     250 IF (sorc == 101) GO TO 290
     260 IF (ni == 0) THEN
       GO TO   320
     END IF
     270 t3 = t2
     t2 = COS(t1)
     t1 = COS(t3)
     GO TO 310
     
!     ODD I
     
     280 IF (sorc == 101) GO TO 260
     290 IF (ni == 0) THEN
       GO TO   330
     END IF
     300 t1 = SIN(t1)
     t2 = SIN(t2)
     310 coef = rz(i3)*(t2-t1)/ni
     IF (sorc == 101 .AND. (i == 2 .OR. i == 5)) coef = -coef
     GO TO 340
     320 coef = 0.0
     GO TO 340
     330 coef = rz(i3)*(rz(i5)-rz(i4))*raddeg
     
     340 z(iat  ) = z(2) + j*1000000
     z(iat+1) = i
     rz(iat+2)= coef
     IF (nogo) CYCLE
     CALL WRITE (geom(4),z(iat),3,noeor)
   END DO
   IF (nogo) CYCLE
   CALL WRITE (geom(4),neg111(1),3,noeor)
 END DO
 GO TO 170
 
 370 GO TO icont, (380,390,400,410)
 
!     S-SET MPC-S FROM SECTAX CARDS
!     =============================
 
!     DO SECTAX CARDS FOR S-SET.
 
 380 REC(1) = sectax(1)
 REC(2) = sectax(2)
 REC(3) = sectax(3)
 n3or5  = 5
 k3or6  = 6
 sorc   = 101
 ncard  = 15
 ASSIGN 390 TO icont
 IF (nosect == 0) THEN
   GO TO   390
 ELSE
   GO TO   160
 END IF
 
!     C-SET MPC-S FROM POINTAX CARDS
!     ==============================
 
 
 390 REC(1) = pointx(1)
 REC(2) = pointx(2)
 REC(3) = pointx(3)
 n3or5  = 3
 k3or6  = 6
 sorc   = 102
 ASSIGN 400 TO icont
 IF (nopont == 0) THEN
   GO TO   400
 ELSE
   GO TO   160
 END IF
 
!     C-SET MPC-S FROM SECTAX CARDS
!     =============================
 
 400 REC(1) = sectax(1)
 REC(2) = sectax(2)
 REC(3) = sectax(3)
 n3or5  = 5
 k3or6  = 6
 sorc   = 102
 ASSIGN 410 TO icont
 IF (nosect == 0) THEN
   GO TO   410
 ELSE
   GO TO   160
 END IF
 
!     BALANCE OF MPCAX CARDS
 
 410 IF (last == 0) THEN
   GO TO   420
 ELSE
   GO TO   510
 END IF
 420 CALL locate (*510,z(ibuff1),mpcax(1),flag)
 ncard = 10
 IF (nwords == 0) GO TO 440
 DO  i = 1,nwords
   CALL READ (*1540,*470,axic,z(1),1,noeor,iamt)
 END DO
 
!     NOW POSITIONED AT POINT LEFT OFF AT ABOVE.
 
 440 CALL READ (*1540,*510,axic,setid,1,noeor,iamt)
 IF (setid < 101) GO TO 470
 IF (setid > 102) GO TO 448
 nogo = .true.
 CALL page2(3)
 imsg = 366
 WRITE  (nout,445) ufm,imsg
 445 FORMAT (a23,i5)
 WRITE  (nout,442)
 442 FORMAT (5X,'SPCAX OR MPCAX CARD HAS A SETID = 101 OR 102.  101 ',  &
     'AND 102 ARE SYSTEM ID-S RESERVED FOR SINE AND COSINE SETS')
 448 IF (nogo) GO TO 450
 CALL WRITE (geom(4),setid,1,noeor)
 450 CALL READ (*1540,*100,axic,z(1),4,noeor,iamt)
 IF (z(4) == (-1)) GO TO 500
 
!     CHECK HARMONIC NUMBER
 
 nnn = z(2)
 ASSIGN 460 TO ierrtn
 GO TO 1420
 
!     CHECK RING ID
 
 460 nnn = z(1)
 ASSIGN 490 TO ierrtn
 GO TO 1440
 470 CALL page2 (3)
 imsg = 1063
 WRITE  (nout,105) sfm,imsg
 WRITE  (nout,480) cdtype(19),cdtype(20)
 480 FORMAT (5X,'EOR ON AXIC FILE WHILE READING ',2A4,'CARD RECORDS.')
 nogo = .true.
 GO TO 1530
 
 490 z(2) = z(1) + (z(2)+1)*1000000
 IF (nogo) GO TO 450
 CALL WRITE (geom(4),z(2),3,noeor)
 GO TO 450
 
!     END OF EQUATION
 
 500 IF (nogo) GO TO 440
 CALL WRITE (geom(4),neg111(1),3,noeor)
 GO TO 440
 
!     AT 713(?) WRITE EOR AND PUT BITS IN TRAILER.
 
 510 IF (mpcon == 0) THEN
   GO TO   530
 END IF
 520 iamt   = 0
 REC(1) = mpc(1)
 REC(2) = mpc(2)
 REC(3) = mpc(3)
 ASSIGN 530 TO iretrn
 GO TO 1300
 
!     OMITAX CARDS
 
 530 REC(1) = omitax(1)
 REC(2) = omitax(2)
 REC(3) = omitax(3)
 ncard  = 11
 rec1(1)= omit(1)
 rec1(2)= omit(2)
 rec1(3)= omit(3)
 ASSIGN 600 TO icont
 540 CALL locate (*590,z(ibuff1),REC(1),flag)
 IF (nogo) GO TO 550
 CALL WRITE (geom(4),rec1(1),3,noeor)
 550 CALL READ (*1540,*580,axic,z(1),3,noeor,iamt)
 
!     CHECK HARMONIC NUMBER
 
 nnn = z(2)
 ASSIGN 560 TO ierrtn
 GO TO 1420
 
!     CHECK RING ID
 
 560 nnn = z(1)
 ASSIGN 570 TO ierrtn
 GO TO 1440
 
 570 z(2) = z(1) + (z(2)+1)*1000000
 IF (ifpdco(z(3))) GO TO 571
 DO  l2 = 1,6
   IF (ll(l2) == 0) CYCLE
   z(3) = ll(l2)
   IF (nogo) GO TO 550
   CALL WRITE (geom(4),z(2),2,noeor)
 END DO
 GO TO 550
 571 nogo = .true.
 CALL page2 (3)
 imsg = 367
 WRITE  (nout,445) ufm,imsg
 WRITE  (nout,573) z(3),cdtype(2*ncard-1),cdtype(2*ncard)
 573 FORMAT (5X,'COMPONENT SPECIFICATION',i8,4H on ,2A4, ' CARD IS INCORRECT')
 GO TO 550
 
!     WRITE EOR AND PUT BITS IN TRAILER
 
 580 iamt   = 0
 REC(1) = rec1(1)
 REC(2) = rec1(2)
 REC(3) = rec1(3)
 ASSIGN 590 TO iretrn
 GO TO 1300
 590 GO TO icont, (600,870)
 
!     SPCADD CARD
!     ===========
 
 600 REC(1)  = spcadd(1)
 REC(2)  = spcadd(2)
 REC(3)  = spcadd(3)
 rec1(1) = spcax(1)
 rec1(2) = spcax(2)
 rec1(3) = spcax(3)
 ASSIGN 610 TO icont
 GO TO 21
 
!     SPCAX CARD
!     ==========
 
 610 REC(1) = spc(1)
 REC(2) = spc(2)
 REC(3) = spc(3)
 
!     RECORD HEADER FOR SPC-S
 
 ASSIGN 620 TO iheadr
 GO TO 1470
 
 620 last  = -1
 ncard = 18
 CALL locate (*670,z(ibuff1),spcax(1),flag)
 last   = 0
 nwords = 0
 630 CALL READ (*1540,*660,axic,z(1),5,noeor,iamt)
 IF (z(1) > 100) GO TO 670
 nwords = nwords + 5
 
!     ALTER CARD JUST READ AND OUTPUT
 
!     CHECK HARMONIC NUMBER
 
 nnn = z(3)
 ASSIGN 640 TO ierrtn
 GO TO 1420
 
!     CHECK RING ID
 
 640 nnn = z(2)
 ASSIGN 650 TO ierrtn
 GO TO 1440
 
 650 z(2) = z(2) + (z(3)+1)*1000000
 z(3) = z(4)
 z(4) = z(5)
 IF (nogo) GO TO 630
 
 CALL WRITE (geom(4),z(1),4,noeor)
 GO TO 630
 660 last = 1
 
!     FIRST NWORDS HAVE BEEN PROCESSED OF SPCAX CARDS
!     UNLESS LAST = 1, IN WHICH CASE ALL SPCAX CARDS ARE COMPLETE.
!     IF LAST = -1, THERE ARE NO SPCAX CARDS
 
!     S-SET AND C-SET SPC-S FROM RINGAX CARDS
!     =======================================
 
 670 sorc   = 101
 ncard  = 14
 compon = 135
 IF (iconso == 1) compon = 13
 ASSIGN 750 TO icont
 680 CALL locate (*760,z(ibuff1),ringax(1),flag)
 690 CALL READ (*1540,*740,axic,z(1),4,noeor,iamt)
 
 IF (sorc == 102) GO TO 730
 
!     GIVE RING CARD A CHECK FOR MINIMUM DATA.
 
!     CHECK RING ID
 
 nnn = z(1)
 ASSIGN 700 TO ierrtn
 GO TO 1440
 
!     CHECK FOR NON-ZERO RADIUS
 
 700 IF (rz(i3-1) == 0.0) THEN
   GO TO   710
 ELSE
   GO TO   730
 END IF
 710 CALL page2 (3)
 imsg = 368
 WRITE  (nout,445) ufm,imsg
 WRITE  (nout,720) z(1)
 720 FORMAT (5X,'RINGAX CARD WITH RING ID =',i10,' HAS A ZERO RADIUS',  &
     ' SPECIFIED.')
 nogo = .true.
 730 z(4) = 0
 z(3) = compon
 z(2) = z(1) + 1000000
 z(1) = sorc
 IF (nogo) GO TO 690
 CALL WRITE (geom(4),z(1),4,noeor)
 GO TO 690
 
 740 GO TO icont, (750,770)
 750 sorc   = 102
 compon = 246
 IF (iconso == 1) compon = 2
 
!     KEEP DOF 4 FOR PIEZOELECTRIC PROBLEM
 
 IF (ipiez == 1) compon = 26
 ASSIGN 770 TO icont
 GO TO 680
 
!     MISSING REQUIRED CARD
 
 760 ASSIGN 770 TO ierrtn
 GO TO 1510
 
!     BALANCE OF SPCAX CARDS
 
 770 IF (last == 0) THEN
   GO TO   780
 ELSE
   GO TO   830
 END IF
 780 CALL locate (*830,z(ibuff1),spcax(1),flag)
 ncard = 18
 IF (nwords == 0) GO TO 800
 DO  i = 1,nwords,5
   CALL READ (*1540,*840,axic,z(1),5,noeor,iamt)
 END DO
 
!     NOW POSITIONED AT POINT LEFT OFF AT ABOVE...
 
 800 CALL READ (*1540,*830,axic,z(1),5,noeor,iamt)
 IF (z(1) < 101) GO TO 840
 IF (z(1) > 102) GO TO 808
 nogo = .true.
 CALL page2 (3)
 imsg = 366
 WRITE (nout,445) ufm,imsg
 WRITE (nout,442)
 
!     CHECK HARMONIC NUMBER
 
 808 nnn = z(3)
 ASSIGN 810 TO ierrtn
 GO TO 1420
 
!     RING ID CHECK
 
 810 nnn = z(2)
 ASSIGN 820 TO ierrtn
 GO TO 1440
 
 820 z(2) = z(2) + (z(3)+1)*1000000
 z(3) = z(4)
 z(4) = z(5)
 IF (nogo) GO TO 800
 CALL WRITE (geom(4),z(1),4,noeor)
 GO TO 800
 
!     WRITE EOR AND PUT BITS IN THE TRAILER
 
 830 iamt = 0
 ASSIGN 860 TO iretrn
 GO TO 1300
 840 CALL page2 (3)
 imsg = 1063
 WRITE (nout,105) sfm,imsg
 WRITE (nout,480) cdtype(35),cdtype(36)
 nogo = .true.
 GO TO 1530
 
!     SUPAX CARDS
!     ===========
 
 860 REC(1)  = supax(1)
 REC(2)  = supax(2)
 REC(3)  = supax(3)
 ncard   = 19
 rec1(1) = suport(1)
 rec1(2) = suport(2)
 rec1(3) = suport(3)
 ASSIGN 870 TO icont
 GO TO 540
 
!     CLOSE GEOM4
 
 870 i = 4
 ASSIGN 880 TO iretrn
 GO TO 1380
 
 
!     GEOM1 PROCESSING
!     ================
 
!     OPEN GEOM1
 
 880 ifile = geom(1)
 i    = 1
 op   = outrwd
 buff = ibuff2
 ASSIGN 890 TO iretrn
 GO TO 1340
 
!     GRID CARDS FROM POINTAX AND SECTAX CARDS
 
!     NOPONT = 0 OR 1, DEPENDING ON THE PRESSENCE OF POINTAX CARDS
!     NOSECT = 0 OR 1, DEPENDING ON THE PRESSENCE OF SECTAX  CARDS
 
!     RECORD HEADER FOR GRID CARDS
 
 890 REC(1) = grid(1)
 REC(2) = grid(2)
 REC(3) = grid(3)
 ASSIGN 900 TO iheadr
 GO TO 1470
 
 900 IF (nosect == 0) THEN
   GO TO   910
 ELSE
   GO TO   920
 END IF
 910 IF (nopont == 0) THEN
   GO TO  1110
 ELSE
   GO TO   980
 END IF
 920 IF (nopont == 0) THEN
   GO TO   940
 END IF
 
!     LOCATE SECTAX CARDS, READ SECTAX, CONVERT TO GRID, PUT ON NFILE
 
 930 nfile = scrtch
 
!     OPEN SCRTCH FILE
 
 i   = 5
 op  = outrwd
 buff= ibuff3
 ASSIGN 950 TO iretrn
 GO TO 1340
 
 940 nfile = geom(1)
 
 950 icard = 15
 CALL locate (*1090,z(ibuff1),sectax(1),flag)
 960 CALL READ (*1540,*970,axic,z(1),5,noeor,iamt)
 z(2) = 0
 z(6) = csid
 z(7) = 0
 z(8) = 0
 IF (nogo) GO TO 960
 CALL WRITE (nfile,z(1),8,noeor)
 GO TO 960
 970 IF (nopont == 0) THEN
   GO TO  1110
 END IF
 980 icard = 12
 CALL locate (*1090,z(ibuff1),pointx(1),flag)
 
!     READ POINT CARD CONVERT TO GRID CARD AND PUT OUT ON GEOM(1)
!     MERGING GRID CARDS FROM SCRTCH IF NOSECT IS NON-ZERO
 
 IF (nosect == 0) THEN
   GO TO  1000
 END IF
 990 IF (nogo  ) GO TO 1110
 CALL CLOSE (scrtch,clorwd)
 CALL OPEN (*1570,scrtch,z(ibuff3),inprwd)
 CALL READ (*1050,*1050,scrtch,z(9),8,noeor,iamt)
 1000 CALL READ (*1540,*1070,axic,z(1),3,noeor,iamt)
 
!     CONVERT POINTAX CARD
 
 z(2)   = 0
 rz(i4) = 0.0
 rz(i5) = 0.0
 z(6)   = csid
 z(7)   = 0
 z(8)   = 0
 IF (nosect == 0) THEN
   GO TO  1020
 END IF
 1010 IF (z(1) >= z(9)) GO TO 1030
 1020 zpt = 1
 GO TO 1040
 1030 zpt = 9
 1040 IF (nogo) GO TO 1110
 CALL WRITE (geom(1),z(zpt),8,noeor)
 IF (zpt == 1) GO TO 1000
 CALL READ (*1050,*1050,scrtch,z(9),8,noeor,iamt)
 IF (nopont == 0) THEN
   GO TO  1040
 ELSE
   GO TO  1010
 END IF
 1050 nosect = 0
 
!     CLOSE SCRTCH
 
 i = 5
 ASSIGN 1060 TO iretrn
 GO TO 1380
 1060 IF (nopont == 0) THEN
   GO TO  1110
 ELSE
   GO TO  1020
 END IF
 
 1070 IF (nosect == 0) THEN
   GO TO  1110
 END IF
 1080 zpt = 9
 nopont = 0
 GO TO 1040
 
 1090 CALL page2 (3)
 imsg = 1064
 WRITE  (nout,105) sfm,imsg
 WRITE  (nout,1100) cdtype(2*icard-1),cdtype(2*icard)
 1100 FORMAT (5X,2A4,' CARD COULD NOT BE LOCATED ON AXIC FILE AS ',  &
     'EXPECTED.')
 nogo = .true.
 GO TO 1110
 
!     GRID CARDS FROM RING CARDS
 
!     COPY RINGAX CARDS INTO CORE AND TO SCRTCH IF CORE IS EXCEEDED.
 
 1110 CALL locate (*1240,z(ibuff1),ringax(1),flag)
 nwords = (icore/4)*4 - 12
 ibegin = 13
 iscrat = 0
 CALL READ (*1540,*1140,axic,z(13),nwords,noeor,iamt)
 
!     FALL HERE IMPLIES CORE IS FULL.. SPILL BALANCE TO SCRTCH FILE.
 
 ion    = 0
 iscrat = 0
 IF (nogo) GO TO 1240
 CALL OPEN (*1570,scrtch,z(ibuff3),outrwd)
 1120 CALL READ (*1540,*1130,axic,z(1),8,noeor,iamt)
 ion = 1
 CALL WRITE (scrtch,z(1),8,noeor)
 GO TO 1120
 1130  IF ((iamt/4)*4 /= iamt) GO TO 1230
 IF (ion == 0 .AND. iamt == 0) GO TO 1160
 iscrat = 1
 IF (nogo) GO TO 1240
 CALL WRITE (scrtch,z(1),iamt,eor)
 CALL CLOSE (scrtch,clorwd)
 GO TO 1160
 
 1140 IF ((iamt/4)*4 /= iamt) GO TO 1230
 nwords = iamt
 
!     NWORDS-WORDS ARE IN CORE AND IF ISCRAT = 1 THERE IS
!     A RECORD OF RINGAX CARDS ON SCRTCH FILE ALSO
 
!     NOW MAKE N+1 PASSES THROUGH THE RING CARDS
 
 1160 IF (iscrat == 0) THEN
   GO TO  1180
 END IF
 1170 IF (nogo  ) GO TO 1240
 CALL OPEN (*1570,scrtch,z(ibuff3),inprwd)
 1180 z(2) = 0
 z(5) = 0
 z(6) = csid
 z(8) = 0
 ncards = nwords/4
 
!     27TH WORD OF SYSTEM IS PACKED AND HOLDS NUMBER OF RINGS AND HARMS
 
 mn   = nplus1
 isystm(161) = ncards
 nadd = 0
 DO  i = 1,nplus1
   nadd = nadd + 1000000
   ipt  = ibegin - 4
   
!     PASS THROUGH THE INCORE CARDS
   
   DO  j = 1,ncards
     ipt  = ipt + 4
     z(1) = z(ipt) + nadd
     z(3) = z(ipt+1)
     z(4) = z(ipt+2)
     z(7) = z(ipt+3)
     IF (nogo) CYCLE
     CALL WRITE (geom(1),z(1),8,noeor)
   END DO
   
!     PASS THROUGH SCRTCH CARDS IF ANY
   
   IF (nogo  ) CYCLE
   IF (iscrat == 0) THEN
     GO TO  1220
   END IF
   1200 CALL READ (*1540,*1210,scrtch,z(9),4,noeor,iamt)
   z(1) = z(9) + nadd
   z(3) = z(10)
   z(4) = z(11)
   z(7) = z(12)
   CALL WRITE (geom(1),z(1),8,noeor)
   GO TO 1200
   
   1210 CALL REWIND (scrtch)
   1220 CONTINUE
 END DO
 
!     PUT BITS IN TRAILER AND WRITE EOR FOR GRID CARDS
 
 iamt   = 0
 REC(1) = grid(1)
 REC(2) = grid(2)
 REC(3) = grid(3)
 ASSIGN 1240 TO iretrn
 GO TO 1300
 1230 ncard  = 14
 ASSIGN 1240 TO ierrtn
 GO TO 1490
 
!     SEQGP CARD
!     ==========
 
 1240 REC(1) = seqgp(1)
 REC(2) = seqgp(2)
 REC(3) = seqgp(3)
 ASSIGN 1250 TO iretrn
 GO TO 1260
 
!     CLOSE GEOM1
 
 1250 i = 1
 ASSIGN 1530 TO iretrn
 GO TO 1380
 
 
!     UTILITY SECTION FOR IFP3
!     AXIS-SYMETRIC-CONICAL-SHELL DATA GENERATOR.
!     ==========================================
 
!     COMMON CODE FOR TRANSFER OF RECORD FROM AXIC FILE TO SOME
!     OTHER FILE
 
 1260 CALL locate (*1330,z(ibuff1),REC(1),flag)
 IF (nogo) GO TO 1330
 CALL WRITE (ifile,REC(1),3,noeor)
 1290 CALL READ (*1540,*1300,axic,z(1),icore,noeor,iamt)
 iamt = icore
 CALL WRITE (ifile,z(1),iamt,noeor)
 GO TO 1290
 1300 IF (nogo) GO TO 1330
 CALL WRITE (ifile,z(1),iamt,eor)
 
!     PUT BITS IN TRAILER
 
 i1 = (REC(2)-1)/16 + 2
 i2 =  REC(2)-(i1-2)*16 + 16
 trail(i1) = orf(trail(i1),two(i2))
 
 1330 GO TO iretrn, (590,530,610,30,1240,1250,860,37)
 
!     OPEN A FILE AND GET THE TRAILER
 
 1340 IF (nogo) GO TO 1350
 CALL OPEN (*1360,FILE(i),z(buff),op)
 openfl(i) = 1
 IF (i > 4) GO TO 1350
 
!     WRITE THE HEADER RECORD
 
 CALL WRITE (FILE(i),iname(2*i-1),2,eor)
 trail(1) = FILE(i)
 CALL rdtrl (trail(1))
 
 1350 GO TO iretrn, (890,950,20)
 
 1360 CALL page2 (3)
 imsg = 1061
 WRITE  (nout,105 ) sfm,imsg
 WRITE  (nout,1370) FILE(i),iname(2*i-1),iname(2*i),ifist
 1370 FORMAT (5X,11HFILE NUMBER ,i4,3H ( ,2A4,12H) is NOT in ,a4)
 nogo = .true.
 GO TO 1530
 
!     CLOSE A FILE
 
 1380 IF (openfl(i) == 0.0) THEN
   GO TO  1410
 END IF
 1390 IF (i > 4) GO TO 1400
 CALL WRITE (FILE(i),t65535(1),3,eor)
 1400 CALL CLOSE (FILE(i),clorwd)
 openfl(i) = 0
 IF (i > 4) GO TO 1410
 CALL wrttrl (trail(1))
 1410 GO TO iretrn, (880,1060,1530)
 
!     HARMONIC NUMBER,  ON CARD TYPE ..... IS OUT OF RANGE 0 TO 998
 
 1420 IF (nnn < 999 .AND. nnn >= 0) GO TO ierrtn, (70,460,560,640,810)
 CALL page2 (3)
 imsg = 364
 WRITE  (nout,445 ) ufm,imsg
 WRITE  (nout,1430) nnn,cdtype(2*ncard-1),cdtype(2*ncard)
 1430 FORMAT (5X,'HARMONIC NUMBER',i6,4H on ,2A4,' CARD OUT OF 0 TO ',  &
     '998 ALLOWABLE RANGE')
 nogo = .true.
 GO TO ierrtn, (70,460,560,640,810)
 
!     RING ID OUT PERMISSABLE RANGE OF 1 TO 999999
 
 1440 IF (nnn > 0 .AND. nnn <= 999999)  &
     GO TO ierrtn, (80,180,490,570,650,700,820)
 CALL page2 (3)
 imsg = 365
 WRITE  (nout,445 ) ufm,imsg
 WRITE  (nout,1450) nnn,cdtype(2*ncard-1),cdtype(2*ncard)
 1450 FORMAT (5X,'RING ID',i10,4H on ,2A4,' CARD OUT OF 0 TO 999999',  &
     ' ALLOWABLE RANGE')
 nogo = .true.
 GO TO ierrtn, (80,180,490,570,650,700,820)
 
!     CHECK BIT-IBIT IN TRAILER AND RETURN NON = ZERO OR NON-ZERO
 
 1460 i1 = (ibit-1)/16  +  2
 i2 = ibit - (i1-2)*16 + 16
 non = andf(axtrl(i1),two(i2))
 GO TO ibitr, (140,150)
 
!     WRITE 3 WORD RECORD HEADER
 
 1470 IF (nogo) GO TO 1480
 CALL WRITE (ifile,REC(1),3,noeor)
 1480 GO TO iheadr, (40,170,620,900,28)
 
!     END-OF-RECORD ON AXIC FILE.
 
 1490 CALL page2 (3)
 imsg = 1063
 WRITE (nout,105) sfm,imsg
 WRITE (nout,480) cdtype(2*ncard-1),cdtype(2*ncard)
 nogo = .true.
 GO TO ierrtn, (1240)
 
!     MISSING REQUIRED CARD
 
 1510 CALL page2 (3)
 imsg = 362
 WRITE  (nout,445 ) ufm,imsg
 WRITE  (nout,1520) cdtype(2*ncard-1),cdtype(2*ncard)
 1520 FORMAT (5X,'MINIMUM PROBLEM REQUIRES ',2A4,' CARD.  NONE FOUND.')
 nogo = .true.
 GO TO ierrtn, (770)
 
!     RETURN TO IFP3
 
 1530 RETURN
 
!     EOF ENCOUNTERED READING AXIC FILE
 
 1540 nfile = axic
 CALL page2 (3)
 imsg = 3002
 WRITE  (nout,105 ) sfm,imsg
 WRITE  (nout,1560) iname(11),iname(12),nfile
 1560 FORMAT (5X,'EOF ENCOUNTERED WHILE READING DATA SET ',2A4,' (FILE',  &
     i4,') IN SUBROUTINE IFP3B')
 nogo = .true.
 GO TO 1530
 
 1580 CALL page2 (3)
 imsg = 363
 WRITE  (nout,445 ) ufm,imsg
 WRITE  (nout,1590)
 1590 FORMAT (5X,'INSUFFICIENT CORE TO PROCESS AXIC DATA IN SUBROUTINE',  &
     'IFP3B')
 nogo = .true.
 GO TO 1530
 
 1570 i = 5
 GO TO 1360
END SUBROUTINE ifp3b
