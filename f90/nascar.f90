SUBROUTINE nascar
     
!     NASCAR READS THE NASTRAN CARD (IF PRESENT) AND CALLS TTLPGE.
 
 IMPLICIT INTEGER (a-z)
 EXTERNAL        rshift,orf,complf
 INTEGER :: hdg(14),nstrn(2),bdt(7),files(2), modcom(9),keywds(2,17),buf(75)
 REAL :: s1,rtolel
 CHARACTER (LEN=16) :: s2
 CHARACTER (LEN=25) :: uwm
 CHARACTER (LEN=23) :: ufm
 COMMON /xmssg / ufm,uwm
 COMMON /machin/ mach
 COMMON /system/ system(100)
 COMMON /BLANK / flag,card(20)
 COMMON /output/ pghdg(96)
 COMMON /xfist / nfist,lfist,fist(2)
 COMMON /xpfist/ npfist
 COMMON /lhpwx / lhpw(4),mxfl
 EQUIVALENCE     (system( 1),sysbuf), (system( 2),outtap),  &
     (system( 3),nogo  ), (system( 4),intap ),  &
     (system( 7),logfl ), (system(20),pltflg),  &
     (system(29),maxfil), (system(30),maxopn),  &
     (system(34),idrum ), (system(57),modcom(1)),  &
     (system(70),itolel,rtolel), (system(77),bandit)
 
 DATA nstrn  /4HNAST,  4HRAN /
 DATA files  /4HFILE,  1HS   /, BLANK / 1H  /
 DATA lkeywd /   17 /
 DATA keywds /4HBUFF,  4HSIZE, 4HCONF,  4HIG  ,  &
     4HMAXF,  4HILES, 4HMAXO,  4HPEN ,  &
     4HSYST,  4HEM  , 4HKON3,  4H60  ,  &
     4HNLIN,  4HES  , 4HTITL,  4HEOPT,  &
     4HMODC,  4HOM  , 4HHICO,  4HRE  ,  &
     4HDRUM,  4H    , 4HTRAC,  4HKS  ,  &
     4HSTST,  4H    , 4HBAND,  4HIT  ,  &
     4HBULK,  4HDATA, 4HPLOT,  4HOPT ,  &
     4HLOGF,  4HL   /
 DATA hdg    /4HN a ,4HS t ,4HR a ,4HN  s,4H y s,4H t e,4H m  ,  &
     4HP a ,4HR a ,4HM e ,4HT e ,4HR  e,4H c h,4H o  /
 DATA bdt    /4HTCRI,4HTMTH,4HTMPC,4HTDEP,4HTPCH,4HTRUN,4HTDIM/
 DATA topt   /  -9  /
!     DATA ADD    /4H@ASG,4H,T  ,4HLOG-,4HFILE,4H.,F ,4H .  /
 DATA s1,s2  /4HWORD, ' OF /SYSTEM/ IS '/
 
 
!     CONMSG (BCD7,1,1)
 mask1 = complf(0)
 mask2 = rshift(mask1,1)
 
!     CALL NSINFO TO OPEN NAINFO FILE AND PICK UP ANY PRESET SYSTEM
!     PARAMETERS FROM THE SECOND SECTION OF THAT FILE
 
 j = 2
 CALL nsinfo (j)
 IF (j /= 2) topt = j
 
!     READ FIRST CARD IN DATA STREAM AND CALL XRCARD TO CONVERT IT.
!     IF INPUT CARD IS BLANK, READ NEXT CARD
 
 10 CALL xread (*3500,card)
 IF (card(1) == BLANK .AND. card(2) == BLANK .AND. card(3) == BLANK  &
     .AND. card(5) == BLANK .AND. card(7) == BLANK) GO TO 10
 CALL xrcard (buf,75,card)
 flag = 1
 IF (buf(1) < 0.0) THEN
   GO TO  4000
 ELSE IF (buf(1) == 0.0) THEN
   GO TO    15
 ELSE
   GO TO    20
 END IF
 15 IF (nogo == 0) GO TO 10
 nogo = 0
 GO TO 4000
 
!     IF CARD IS NASTRAN PARAMETER CARD, ECHO IT.
 
 20 IF (buf(2) /= nstrn(1) .OR. buf(3) /= nstrn(2)) GO TO 4000
 DO  i = 1,14
   pghdg(i+2) = hdg(i)
 END DO
 IF (system(11) <= 0) CALL page1
 WRITE  (outtap,30) (card(i),i=1,20)
 30 FORMAT (5X,20A4)
 
!     RETURN IF NO KEYWORD ON NASTRAN CARD
 
 IF (buf(4) == mask2) GO TO 4000
 flag = 0
 
!     IDENTIFY KEYWORDS AND BRANCH TO APPROPRIATE CODE.
 
 j  = 4
 35 jn = 2*buf(1) + 1
 j1 = 1
 GO TO 50
 40 IF (buf(j1) < 0.0) THEN
   GO TO    85
 ELSE IF (buf(j1) == 0.0) THEN
   GO TO    45
 ELSE
   GO TO    50
 END IF
 45 CALL xread  (*3500,card)
 CALL xrcard (buf,75,card)
 WRITE (outtap,30) card
 IF (buf(1) == 0) GO TO 45
 j = 2
 GO TO 35
 50 IF (buf(j1) == mask2) GO TO 4000
 DO  i = 1,lkeywd
   IF (buf(j) == keywds(1,i) .AND. buf(j+1) == keywds(2,i)) GO TO 110
 END DO
 IF (buf(j) == keywds(1,14)) GO TO 100
 IF (buf(j) == files(1)) GO TO 3000
 IF (buf(j) /= BLANK) GO TO 60
 j = j + 2
 IF (buf(j+2) == mask2) GO TO 4000
 IF (buf(j) == 0) GO TO 45
 IF (j < jn) GO TO 50
 j1 = jn + 1
 IF (buf(j1) == mask2) GO TO 4000
 jn = 2*buf(j1) + 1
 j  = j1 + 1
 GO TO 50
 
!     PRINT MESSAGE FOR UNIDENTIFIED KEYWORD.
 
 60 CONTINUE
 WRITE  (outtap,65) ufm,buf(j),buf(j+1)
 65 FORMAT (a23,' 17, UNIDENTIFIED NASTRAN CARD KEYWORD ',2A4,  &
     '.  ACCEPTABLE KEYWORDS FOLLOW ---', /1H0 )
 DO  i = 1,lkeywd
   WRITE (outtap,75) keywds(1,i),keywds(2,i)
 END DO
 75 FORMAT (5X,2A4)
 WRITE  (outtap,80) (bdt(i),i=1,7)
 80 FORMAT (7(5X,4HBAND,a4), /5X,'FILES (MUST BE LAST IN INPUT LIST)')
 nogo = 1
 GO TO 4000
 
!     PRINT MESSAGE FOR BAD FORMAT.
 
 85 WRITE  (outtap,90) ufm
 90 FORMAT (a23,' 43, INCORRECT FORMAT FOR NASTRAN CARD.')
 nogo = 1
 GO TO 4000
 
!     . CHECK FOR LEGAL REAL NUMBER...
 
 95 CONTINUE
 IF (buf(j1-2) /= -2) GO TO 85
 IF (buf(j1-2) == -2) GO TO 120
 IF (i == 11) GO TO 1100
 
!     . BANDIT KEYWORDS.
 
 100 i = 1400
 k = buf(j+1)
 
!     KEYWORD FOUND.
 
 110 CONTINUE
 j1 = jn + 1
 param = buf(j1+1)
 j1 = j1 + 2
 IF (buf(j1) /= mask2) jn = 2*buf(j1) + j1
 j  = j1 + 1
 IF (buf(j1-2) /= -1) GO TO 95
 120 CONTINUE
 IF (i == 1400) GO TO 1400
 SELECT CASE ( i )
   CASE (    1)
     GO TO  150
   CASE (    2)
     GO TO  200
   CASE (    3)
     GO TO  300
   CASE (    4)
     GO TO  400
   CASE (    5)
     GO TO  500
   CASE (    6)
     GO TO  600
   CASE (    7)
     GO TO  700
   CASE (    8)
     GO TO  800
   CASE (    9)
     GO TO  900
   CASE (   10)
     GO TO 1000
   CASE (   11)
     GO TO 1100
   CASE (   12)
     GO TO 1200
   CASE (   13)
     GO TO 1300
   CASE (   14)
     GO TO 1450
   CASE (   15)
     GO TO 1500
   CASE (   16)
     GO TO 1600
   CASE (   17)
     GO TO 1700
 END SELECT
 
!     BUFFSIZE
 
 150 CONTINUE
 sysbuf = param
 GO TO 40
 
!     IGNORE THE CONFIG PARAMETER
 
 200 CONTINUE
 GO TO 40
 
!     MAXFILES UPPER LIMIT
 
 300 CONTINUE
 m = mxfl
 IF (param <= m) GO TO 320
 WRITE  (outtap,310) m
 310 FORMAT (' *** MAXFILES IS RESET TO THE LIMIT OF 74')
 param = m
 320 maxfil = param
 GO TO 40
 
!     MAXOPEN
 
 400 CONTINUE
 IF (param <= maxfil) GO TO 420
 IF (param >   mxfl) GO TO 430
 WRITE  (outtap,410) param
 410 FORMAT (' *** MAXOPEN EXCEEDS MAXFILES. MAXFILES IS AUTOMATICALLY'  &
     ,       ' EXPANDED TO',i4)
 maxfil = param
 420 maxopn = param
 GO TO 40
 
 430 m = mxfl
 WRITE  (outtap,440) m,m
 440 FORMAT (' *** MAXOPEN EXCEEDS MAXFILES LIMIT OF ',i3,'.  BOTH ',  &
     'MAXOPEN AND MAXFILES ARE RESET ONLY TO ',i3,' EACH')
 maxfil = m
 maxopn = m
 GO TO 40
 
!     SYSTEM
 
 500 CONTINUE
 IF (param <=  0) GO TO 85
 IF (param /= 24) GO TO 510
 WRITE (outtap,505)
 505 FORMAT ('0*** FATAL, USER SHOULD NOT CHANGE THE 24TH WORD OF ',  &
     '/SYSTEM/')
 nogo = 1
 GO TO 40
 
 510 IF (buf(j1) < 0.0) THEN
   GO TO   530
 ELSE IF (buf(j1) == 0.0) THEN
   GO TO    40
 END IF
 520 j1 = jn + 1
 j  = j1 + 1
 IF (buf(j1) == mask2) GO TO 85
 jn = j1 + 2*buf(j1)
 GO TO 510
 530 IF (buf(j1) == -2 .AND. param /= 70) GO TO 85
 
!     IGNORE THE CONFIG PARAMETER
 
 IF (param /= 28) system(param) = buf(j)
 
!     SYSTEM WORD ECHO
 
 IF (param >= 10) GO TO 531
 IF (param ==  1) WRITE (outtap,541) s1,param,s2
 IF (param ==  2) WRITE (outtap,542) s1,param,s2
 IF (param ==  3) WRITE (outtap,543) s1,param,s2
 IF (param ==  4) WRITE (outtap,544) s1,param,s2
 IF (param ==  7) WRITE (outtap,547) s1,param,s2
 IF (param == 7 .AND. mach == 3) WRITE (outtap,548)
 IF (param ==  9) WRITE (outtap,549) s1,param,s2
 GO TO 590
 531 k = param/10
 IF (k <= 9) THEN
    SELECT CASE ( k )
     CASE (    1)
       GO TO 532
     CASE (    2)
       GO TO 532
     CASE (    3)
       GO TO 533
     CASE (    4)
       GO TO 534
     CASE (    5)
       GO TO 534
     CASE (    6)
       GO TO 536
     CASE (    7)
       GO TO 536
     CASE (    8)
       GO TO 536
     CASE (    9)
       GO TO 536
   END SELECT
 END IF
 WRITE (outtap,540) s1,param,s2
 GO TO 590
 532 IF (param == 20) WRITE (outtap,550) s1,param,s2
 IF (param == 28) WRITE (outtap,558) s1,param,s2
 IF (param == 29) WRITE (outtap,559) s1,param,s2
 GO TO 590
 533 IF (param == 30) WRITE (outtap,560) s1,param,s2
 IF (param == 31) WRITE (outtap,561) s1,param,s2
 IF (param == 34 .AND. mach /= 4) WRITE (outtap,564) s1,param,s2
 IF (param == 34 .AND. mach == 4) WRITE (outtap,565) s1,param,s2
 GO TO 590
 534 IF (param == 42) WRITE (outtap,572) s1,param,s2
 IF (param == 45) WRITE (outtap,575) s1,param,s2
 IF (param == 57) WRITE (outtap,577) s1,param,s2
 IF (param == 58) WRITE (outtap,578) s1,param,s2
 IF (param == 59) WRITE (outtap,579) s1,param,s2
 GO TO 590
 536 IF (param >= 60 .AND. param <= 65) WRITE (outtap,577) s1,param,s2
 IF (param == 70) WRITE (outtap,580) s1,param,s2
 IF (param == 77) WRITE (outtap,587) s1,param,s2
 GO TO 590
 540 FORMAT (5X,a4,i3,a16,'NOT AVAILABLE. INPUT IGNORED')
 541 FORMAT (5X,a4,i3,a16,'GINO BUFFER SIZE')
 542 FORMAT (5X,a4,i3,a16,'OUTPUT UNIT')
 543 FORMAT (5X,a4,i3,a16,'NOGO FLAG')
 544 FORMAT (5X,a4,i3,a16,'INPUT UNIT')
 547 FORMAT (5X,a4,i3,a16,'NO. OF CONSOLE LOG MESSAGES')
 548 FORMAT (1H+,31X,'. (95 MAX.)')
 549 FORMAT (5X,a4,i3,a16,'NO. OF LINES PER PAGE. MINIMUM 10')
 550 FORMAT (5X,a4,i3,a16,'PLOT OPTION')
 558 FORMAT (5X,a4,i3,a16,'MACHINE CONFIGURATION (IGNORED)')
 559 FORMAT (5X,a4,i3,a16,'MAX FILES')
 560 FORMAT (5X,a4,i3,a16,'MAX FILES OPEN')
 561 FORMAT (5X,a4,i3,a16,'HI-CORE')
 564 FORMAT (5X,a4,i3,a16,'DRUM FLAG')
 565 FORMAT (5X,a4,i3,a16,'NOS/NOS-BE FLAG')
 572 FORMAT (5X,a4,i3,a16,'SYSTEM RELEASE DATE')
 575 FORMAT (5X,a4,i3,a16,'TAPE BIT')
 577 FORMAT (5X,a4,i3,a16,'DATA EXTRACTED FROM ADUM CARDS')
 578 FORMAT (5X,a4,i3,a16,'MPYAD METHOD SELECTION')
 579 FORMAT (5X,a4,i3,a16,'PLOT TAPE TRACK SPEC')
 580 FORMAT (5X,a4,i3,a16,'SMA1 SINGULAR TOLERANCE')
 587 FORMAT (5X,a4,i3,a16,'BANDIT/BULKDATA FLAG')
 
!     SET BOTTOM LIMIT OF 10 TO NUMBER OF LINES PER PAGE
!     AND FOR UNIVAC ONLY, LIMIT THE CONSOLE LOG MESSAGES TO 95 MAXIMUM
 
 590 IF (param == 9 .AND. system(9) < 10) system(9) = 10
 IF (mach == 3 .AND. param == 7 .AND. system(7) > 95) system(7) = 95
 j1 = j1 + 2
 j  = j1 + 1
 IF (buf(j1) == mask2) GO TO 4000
 jn = j1 + 2*buf(j1)
 GO TO 40
 
!     KON360/HICORE
 
 600 CONTINUE
 system(31) = param
 GO TO 40
 
!     NLINES - BOTTOM-LIMITED TO 10
 
 700 CONTINUE
 system(9) = param
 IF (system(9) < 10) system(9) = 10
 GO TO 40
 
!     TITLEOPT
 
 800 CONTINUE
 topt = param
 IF (mach == 3 .AND. topt <= -2) logfl = 3
 GO TO 40
 
!     MODCOM COMMUNICATION AREA
 
 900 CONTINUE
 IF (param <= 0) GO TO 85
 910 IF (buf(j1) < 0.0) THEN
   GO TO   930
 ELSE IF (buf(j1) == 0.0) THEN
   GO TO    40
 END IF
 920 j1 = jn + 1
 j  = j1 + 1
 IF (buf(j1) == mask2) GO TO 85
 jn = j1 + 2*buf(j1)
 GO TO 910
 930 modcom(param) = buf(j)
 j1 = j1 + 2
 j  = j1 + 1
 IF (buf(j1) == mask2) GO TO 4000
 jn = j1 + 2*buf(j1)
 GO TO 40
 
!     HICORE = LENGTH OF CORE ON UNIVAC, VAX, AND UNIX
 
 1000 CONTINUE
 system(31) = param
 GO TO 40
 
!     UNIVAC - DRUM ALLOCATION, 1 BY POSITIONS, 2 BY TRACKS
!              DEFAULT IS 1,  150 POSITIONS  (GOOD FOR LARGE JOB)
!              IF DRUM IS 2, 1280 TRKS. IS ASSIGNED (SUITABLE FOR
!                 SMALLER JOB)
 
!     CDC - IDRUM (34TH WORD OF /SYSTEM/) IS LENGTH OF FET + DUMMY INDEX
 
 1100 idrum = param
 GO TO 40
 
!     PLOT TAPE TRACK SIZE    TRACK=7 IMPLIES 7 TRACK
!                             TRACK=9 IMPLIES 9 TRACK
 
 1200 IF (param /= 7 .AND. param /= 9) GO TO 1250
 IF (param == 7) system(59) = 1
 IF (param == 9) system(59) = 2
 GO TO 40
 1250 WRITE (outtap,1480) uwm,param,keywds(1,12),keywds(2,12)
 nogo = 1
 GO TO 40
 
!     . ELEMENT SINGULARITY TOLERANCE (A REAL S.P. NUMBER)...
 
 1300 itolel = param
 IF (buf(j1-2) == -1) rtolel = itolel
 GO TO 40
 
!     BANDIT (77TH WORD OF SYSTEM)
!     BANDIT KEYWORDS (DEFAULT VALUES IN BRACKETS, SET BY BGRID ROUTINE)
!        BANDTCRI = (1),2,3,4     CRITERION
!        BANDTMTH = 1,2,(3)       METHOD
!        BANDTMPC = (0),1,2       MPC EQUS. AND RIGID ELEMENTS
!        BANDTDEP = (0),1         DEPENDANT GRID
!        BANDTPCH = (0),1         PUNCH SEQGP CARDS
!        BANDTRUN = (0),1         RUN/SEQGP
!        BANDTDIM = (0),1,2,...,9 SCRATCH ARRAY DIMENSION
!        BANDIT   = -1,(0)        BANDIT SKIP FLAG
!     WHERE,
!        CRITERION = 1, USE RMS WAVEFRONT TO DETERMINE BEST RESULT,
!                  = 2, BANDWIDTH,  =3, PROFILE, OR  =4, MAX WAVEFRONT
!        METHOD    = 1, CM METHOD IS USED,   3, GPS, OR     2, BOTH
!        MPC       = 0, MPC'S AND RIGID ELEM ARE NOT CONSIDERED
!                  = 1, MPC'S AND RIGID ELEM ARE USED IN RESEQUENCING
!                  = 2, ONLY RIGID ELEMENTS ARE USED IN RESEQUENCING
!        DEPEND    = 0, DEPENDANT GRID IS OMITTED IN RESEQUENCING
!                       IF MPC IS NON-ZERO
!                  = 1, DEPENDANT GRIDS ARE INCLUDED
!        PUNCH     = 0, NO SEQGP CARDS PUNCHED
!                  = 1, PUNCH OUT BANDIT GENERATED SEQGP CARDS AND
!                       TERMINATE NASTRAN JOB
!        RUN/SEQGP = 0, BANDIT WOULD QUIT IF THERE IS ONE OR MORE SEQGP
!                       CARD IN THE INPUT DECK
!                  = 1, TO FORCE BANDIT TO BE EXECUTED EVEN IF SEQGP
!                       CARDS ARE PRESENT
!        DIM       = 1,2,...,N, TO SET THE SCRATCH AREA, USED ONLY IN
!                       GPS METHOD, TO N*100. (N IS 9 OR LESS)
!                  = 0, DIMENSION IS SET TO 150
!        BANDIT    =-1, BANDIT COMPUTATION IS SKIPPED UNCONDITIONALLY
!                  = 0, BANDIT WOULD BE EXECUTED IF BULK DATA CONTAINS
!                       NO INPUT ERROR
 
 1400 CONTINUE
 IF (bandit < 0) GO TO 40
 IF (k == bdt(7) .AND. param >= 100) param = param/100
 IF (param < 0 .OR. param > 9) GO TO 1470
 DO  i = 1,7
   IF (k /= bdt(i)) CYCLE
   k = param*10**(i-1)
   GO TO 1430
 END DO
 GO TO 60
 1430 bandit = bandit + k
 GO TO 40
 1450 IF (param < 0) bandit = -1
 IF (param <= 0) GO TO 40
 k = keywds(2,14)
 1470 WRITE  (outtap,1480) uwm,param,keywds(1,14),k
 1480 FORMAT (a25,' 65, ILLEGAL VALUE OF ',i7,' IN NASTRAN ',2A4, ' CARD')
 GO TO 40
 
!     BULK DATA CHECK ONLY
!     TO TERMINATE JOB AFTER BULK DATA CHECK, AND SKIP OVER BANDIT
!     (OPTION TO PRINTOUT TIME CONSTANTS IN /NTIME/, IF BULKDATA=-3)
 
 1500 CONTINUE
 IF (param  /=  0) bandit = -2
 IF (bandit == -2) maxfil = 23
 IF (param  == -3) bandit = -3
 GO TO 40
 
!     PLOT OPTIONS -
 
!     PLTFLG   BULKDATA    PLOT COMMANDS       ACTION TAKEN
!     -----  ------- --   -------------  -------------------------------
!      0      NO ERROR    NO ERROR       EXECUTES ALL LINKS, NO PLOTS
!             NO ERROR       ERROR       STOPS AFTER LNK1 DATA CHECK
!                ERROR  ERR OR NO ERR    STOPS AFTER LINK1 CHECK
!      1      NO ERROR    NO ERROR       GO, ALL LINKS AND PLOTS
!             NO ERROR       ERROR       STOP AFTER LINK1 DATA CHECK
!                ERROR    NO ERROR       STOP AFTER LINK1 DATA CHECK
!                ERROR       ERROR       STOP AFTER LINK1 DATA CHECK
!      2      NO ERROR    NO ERROR       STOP AFTER UNDEFORM PLOT/LINK2
!             NO ERROR       ERROR       STOP AFTER LINK1 DATA CHECK
!                ERROR    NO ERROR       STOP AFTER UNDEFORM PLOT/LINK2
!                ERROR       ERROR       STOP AFTER LINK1 DATA CHECK
!      3      (ERROR OR   (ERROR OR      (ATTEMPT TO PLOT UNDEFORM MODEL
!             NO ERROR)   NO ERROR)      THEN STOP/LINK2)
!      4      NO ERROR    NO ERROR       GO, ALL LINKS AND PLOTS
!             NO ERROR       ERROR       STOP AFTER UNDEFORM PLOT/LINK2
!                ERROR    NO ERROR       STOP AFTER UNDEFORM PLOT/LINK2
!                ERROR       ERROR       STOP AFTER LINK1 DATA CHECK
!      5      NO ERROR    NO ERROR       GO, ALL LINKS AND PLOTS
!             NO ERROR       ERROR       GO, ALL LINKS BUT NO PLOTS
!                ERROR    NO ERROR       STOP AFTER UNDEFORM PLOT/LINK2
!                ERROR       ERROR       STOP AFTER LINK1 DATA CHECK
!     PLTFLG 0 OR 1 IS SET BY THE PRESENCE OF THE PLOT TAPE.
!     PLTFLG WILL BE RESET TO POSITIVE IN IFP1
!     CUT MAXFIL TO HALF AND SKIP BANDIT IF PLOT OPTION IS 2 OR 3
 
 1600 CONTINUE
 IF (param < 2 .OR. param > 5) GO TO 1650
 IF (param >= 2) pltflg = -param
 IF (pltflg /= -2 .AND. pltflg /= -3) GO TO 40
 maxfil = 24
 bandit = -1
 GO TO 40
 1650 WRITE (outtap,1480) uwm,param,keywds(1,16),keywds(2,16)
 GO TO 40
 
!     LOGFL = LOGFILE MESSAGE CONTROL ON UNIVAC 1100
 
 1700 CONTINUE
 logfl = param
 GO TO 40
 
!     FILES
 
 3000 IF (buf(j+2) /= mask1) GO TO 85
 IF (j+4 >= jn) GO TO 85
 IF (buf(j+4) == BLANK .AND. buf(j+6) == mask1) GO TO 3010
 j = j + 4
 khr = 0
 GO TO 3020
 3010 j = j + 8
 khr = 7
 3020 IF (buf(j) == mask1 .OR. khr == 1) GO TO 3090
 DO  ii = 1,npfist
   IF (buf(j) == fist(2*ii-1)) GO TO 3060
 END DO
 DO  i = 1,lkeywd
   IF (buf(j) == keywds(1,i) .AND. buf(j+1) == keywds(2,i)) GO TO 110
 END DO
 IF (buf(j) /= BLANK) WRITE (outtap,3050) uwm,buf(j)
 3050 FORMAT (a25,' 64, ',a4,' IS NOT DEFINED AS A NASTRAN FILE AND ',  &
     'WILL BE IGNORED.')
 GO TO 3070
 3060 ixx = 2**(ii-1)
 system(45) = orf(system(45),ixx)
 3070 j = j + 2
 khr = khr + 1
 IF (j < jn) GO TO 3020
 j1 = jn + 1
 IF (buf(j1) == mask2) GO TO 4000
 IF (buf(j1) /= 0) GO TO 85
 3080 CALL xread  (*3500,card)
 CALL xrcard (buf,75,card)
 WRITE (outtap,30) (card(i),i=1,20)
 IF (buf(1) == 0) GO TO 3080
 j  = 2
 j1 = 1
 jn = 2*buf(1) + 1
 GO TO 3020
 3090 j = j + 2
 GO TO 40
 
 
!     END-OF-FILE ENCOUNTERED ON INPUT FILE
 
 3500 WRITE  (outtap,3600) ufm,intap
 3600 FORMAT (a23,' 74, EOF ENCOUNTERED ON UNIT ',i4,  &
     ' WHILE READING THE INPUT DATA IN SUBROUTINE NASCAR')
 CALL mesage (-61,0,0)
 
 
!     GENERATE TITLE PAGE
 
 4000 DO  i = 1,14
   pghdg(i+2) = BLANK
 END DO
 CALL ttlpge (topt)
 
 RETURN
END SUBROUTINE nascar
