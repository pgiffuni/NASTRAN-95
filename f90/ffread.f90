SUBROUTINE ffread (*,card)
     
    !     THIS ROUTINE READS INPUT CARDS IN FREE FIELD OR FIXED FIELD
    !     FORMATS.
 
    !     IF READFILE COMMAND IS ENCOUNTERED, IT SWITCHE THE INPUT FILE TO
    !     THE ONE SPECIFIED BY READFILE UNTIL EOF IS REACHED. THEN IT
    !     SWITCHES BACK TO THE NORMAL CARD READER. NESTED READFILE IS
    !     ALLOWED.
 
    !     IT ALSO PRINTS THE INPUT CARDS IF UNSORTED ECHO FLAG IS ONE
 
    !     ALL INTEGERS, BCD, AND REAL NUMBERS ARE LEFT ADJUSTED BEFORE
    !     RETURNING TO THE CALLER, XSORT2
 
    !     IN BULK DATA SECTION -
    !     ALL INTEGERS ARE LIMITED TO 8 DIGITS. REAL NUMBERS CAN BE UP TO 12
    !     DIGITS IF INPUT VIA FREE-FIELD, OR UP TO 8 DIGITS IF FIXED-FIELD.
    !     ALL REAL NUMBER MUST HAVE A DECIMAL POINT.  10E-6 OR 1+7 ARE
    !     NOT ACCEPTABLE
 
    !     THREE WORDS ARE ATTACHED TO THE END OF AN INPUT CARD TO BE USED
    !     FOR ALPHA-NUMERIC SORTING
 
    CHARACTER (LEN=8), INTENT(OUT)           :: card(10)
    LOGICAL :: fp,       star,     pct,      notyet,   twodot
    INTEGER :: ffflag,   inflag,   none,     screen,   prom,  &
               univc(11),xsort,    wasff
    INTEGER :: ib,       ic,       ie,       is,       il,  &
               ir,       id,       ip,       im,       ig,  &
               ia,       ih,       pt,       sp,       aii,  &
               a1,       dot,      at,       a(94)
    CHARACTER (LEN=1) :: cb,       cc,       ce,       cs,       cl,  &
                         cr,       cd,       cp,       cm,       cg,  &
                         ca,       ch,       ct,       c1,       c(80),  &
                         cx(94),   tmp,      qmark,    dotc
    CHARACTER (LEN=4) :: prompt,   on,       off,      yes,      temp4, echo
    CHARACTER (LEN=5) :: a5,       seqgp,    seqep
    CHARACTER (LEN=8) :: blank,    a8(10),   a81,      cancel,  &
                         save,     rdfl,     skfl,     dend,     dbgn,  &
                         temp,     from,     spill,    list,     help,  &
                         stop,     scale8,   scale1,   noprt,    slash, a8x(12)
    CHARACTER (LEN=48) :: a48

    COMMON /xreadx/ screen,   loop,     kount,    prom,     notyet,  &
                    star,     pct,      jc(9),    l(9),     rc(9), f(9)
    COMMON /qmarkq/ qmark,    tmp(8),   spill,    save(10)
    COMMON /xechox/ ffflag,   iechou,   iechos,   iechop,   xsort,  &
                    wasff,    ncard,    dum(2),   noecho
    COMMON /xxread/ inflag,   insave,   loop4,    ibmcdc,   ierr
    COMMON /machin/ mchn
    COMMON /system/ ibuf,     nout,     nogo,     in

    EQUIVALENCE     (c(1),cx(1),a8(1),a8x(1),a5,a48,a81),  (kkf,fkk),  &
                    (temp4,temp,tmp(1)),(a1,a(1))
 
    DATA            none,     prompt,   on,       off,      yes     /  &
                  4HNONE,   'PROM',   'ON, ',   'OFF,',   'YES,'    /
    DATA            blank,    dend,     dbgn,     from,     slash   /  &
                    '      ','$   END ','$   ...',' FROM-', '/    ' /
    DATA            ct,       xxxx,     cancel,   list,     lout    /  &
                    '.',      4HXXXX,   'CANCEL', 'LIST',   3       /
    DATA            rdfl,               skfl,     dotc,     echo    /  &
                    'READFILE',        'SKIPFILE','.',      'ECHO'  /
    DATA            help,     iwo,      scale8,   scale1,   stop    /  &
                    'HELP',   60,      'SCALE/8','SCALE/10','STOP'  /
    DATA            a,        ib,       noprt,    seqgp,    seqep   /  &
                    94*1H ,   0,       'NOPRINT,','SEQGP',  'SEQEP' /
    DATA            cb , cc , ce , cs , cl , cr , cd , cp , cm , cg /  &
                    ' ', ',', '=', '*', '(', ')', '$', '+', '-', '%'/
    DATA            l12, l94/ 10, 80 /, ca, ch, at / '/', '!', 2H@  /
    DATA            univc   / 4H*ADD,   4H,E  ,   8*4H    ,  4H .   /

    !     THIS ROUTINE IS A PREPROCESSOR FOR THE XREAD ROUTINE IN NASTRAN
    !     WRITTEN BY G. CHAN/SPERRY,  APRIL 1985

    !     FFFLAG IN /XECHOX/ MUST BE SET TO 1234 FOR FREE-FIELD INPUT.
    !     IECHOS IS SET TO -2 IN STAND-ALONE VERSION.
    !     MUST RESERVE 43 WORDS IN SEMINT FOR /XREADX/ IN ALL MACHINES.

    !     FREE FIELD INPUT IS TRIGGERED BY THE PRESENCE OF COMMA (,) OR
    !     EQUAL SIGN (=) IN COLS. 1 THRU 10, AND AFTER BEGIN BULK CARD WAS
    !     READ.

    !     FFREAD IS DESIGNED TO BE USER FRIENDLY -
    !     UNDER NO CIRCUMSTANCES SHOULD THE USER BE KICKED OUT OF THE
    !     COMPUTER DUE TO HIS OR HER STUPID INPUT ERROR(S).

    !     DURING FREE-FIELD INPUT SESSION, FOUR CONTROL CARDS ARE ALLOWED -

    !        ECHO  = SORT, UNSORT, BOTH, NONE, PUNCH, LINK1
    !        PROMPT= ON, OFF, YES    (YES = ON + GENERATED CARD ECHO)
    !        CANCEL= N, TO CANCEL N PREVIOUSLY GENERATED LINES
    !        LIST  =-N, TO   LIST N PREVIOUSLY GENERATED LINES
    !        (CANCEL AND LIST ARE AVAILABLE ONLY IN STAND-ALONE VERSION AND
    !         A SAVE FILE HAS BEEN REQUESTED)

    !     WRITTEN BY G.CHAN/UNISYS ON A COLD DECEMBER MORNING, 1983
    !     REFERENCE - CHAN, G.C.: 'COSMIC/NASTRAN FREE-FIELD INPUT',
    !                 12TH NASTRAN USERS' COLLOQUIUM, MAY 1984

    !     THIS ROUTINE WILL HANDLE COMPUTER WORD OF ANY SIZE, 32,36,60,64
    !     BITS, UPPER CASE AND LOWER CASE ASCII AND EBCDIT CHARACTER SETS.

    !     VAX AND UNIX ONLY -
    !     (UNIVAC TOO, ONLY IF OPEN STATEMENT IS USED FOR LOGICAL UNIT 5)
    !     DURING FREE-FIELD SESSION, 94 COLUMNS, INSTEAD OF REGULARLY 80,
    !     ARE ALLOWED FOR AN INPUT CARD COMING FROM CARD READER OR READFILE
    !     (A MAXINUM OF 94 COLUMNS IS ALLOWED IN PRINT FORMAT 310)

    !     THIS ROUTINE CALLS THE FOLLOWING SUPPORTING SUBROUTINES FOR BCD
    !     (LEFT ADJUSTED), INTEGER, AND F.P. NUMBER CONVERSION -

    !        INT 2 K8  - DECODES INTEGER TO A8 CHAR.
    !        FP  2 K8  - DECODES F.P. NUMBER TO A8 CHAR.
    !        NK1 2 IF  - ENCODES N(A1) CHAR. TO INTEGER OR F.P. NUMBER
    !        NK1 2 K8  - ENCODES N(A1) CHARS. TO A A8 CHAR. WORD
    !        K8  2 INT - DECODES A8 CHAR. TO INTEGER
    !        K8  2 FP  - DECODES A8 CHAR. TO F.P. NUMBER
    !        UPCASE    - REPLACES ANY LOWER-CASE LETTER BY ITS UPPER CASE

    !     THIS ROUTINE WILL ALSO HANDLE 'READFILE' AND 'SKIPFILE' CARDS.
    !     FILE NAME IS LIMITED UP TO 48 CHARACTERS,  8/91

    !     THIS ROUTINE TRIES NOT TO USE SYSTEM ENCODE/DECODE FUNCTIONS,
    !     SHIFT, AND ANY NON-STANDARD CHARACTER FUNCTIONS.


    !     INPUT FILE LOGIC:

    !     IN UNIVAC, INPUT CARDS ARE READ FROM CARD READER INFLAG, UNIT 5.
    !     ALL OTHER INPUT FILES, NESTED OR NOT, ARE DYNAMICALLY INSERTED IN-
    !     TO INPUT STREAM (WITH THE E-O-F MARK STRIPPED OFF), AND READ INTO
    !     COMPUTER SYSTEM FROM UNIT 5 ALSO. IF AN E-O-F MARK ENCOUNTERED
    !     BEFORE ENDDATA CARD, IT IS FATAL. INFLAG=TWO=IN=5

    !     IN ALL OTHER MACHINES, INPUT CARDS ARE READ FROM CARD READER
    !     INFLAG, UNIT 5. WHEN A READFILE CARD IS ENCOUNTERED, DATA ARE READ
    !     INTO COMPUTER SYSTEM FROM UNIT INFLAG, WHICH BEGINS AT 60;
    !          INFLAG = IWO = 60 FOR THE FIRST FILE
    !          INFLAG = 61 FOR THE SECOND FILE
    !          INFLAG = 62 FOR THE THIRD  FILE, ETC.
    !     (NOTE, SINCE NASTRAN USES READFILE INTERNALLY TO READ RIGID FORMAT
    !     FILE, NESTED READFILE IS NOT UNCOMMON)
    !     WHEN E-O-F IS ENCOUNTERED, CURRENT FILE IS CLOSED AND INFLAG IS
    !     DECREASE BY 1. INFLAG IS SET TO ZERO WHEN INFLAG .LE. IWO (END
    !     OF CURRENT NESTED FILE OPERATION). NEXT READFILE, NESTED OR NOT,
    !     IS ALLOWED.

    !     ADD READFILE,NOPRINT OPTION.  2/2/1989
    !     LAST REVISED, 8/1989, IMPROVED EFFICIENCY BY REDUCING CHARACTER
    !     OPERATIONS (VERY IMPORTANT FOR CDC MACHINE)
    !     8/93, LIBERAL READFILE NOPRINT FORMATS:
    !           READFILE,NOPRINT  FILENAME
    !           READFILE,NOPRINT, FILENAME
    !           READFILE NOPRINT  FILENAME
    !           READFILE(NOPRINT) FILENAME
    !           (EMBEDDED BLANK, COMMA, BRACKETS, AND EQUAL-SIGN ALLOWED)
    !           READFILE = FILENAME

    !     INITIALIZE THE FOLLOWING ITEMS SO THAT COMPILER WILL NOT COMPLAIN

    DATA   c1,i,ii,jj,kk / ' ', 4*0 /

    mach = mchn
    IF (mach == 12) mach = 4
    IF (mach < 5) GO TO 40
    l12 = 12
    l94 = 94
40  IF (ib /= 0) GO TO 50
    CALL k2b (cb,ib,1)
    CALL k2b (cc,ic,1)
    CALL k2b (ce,ie,1)
    CALL k2b (cs,is,1)
    CALL k2b (cl,il,1)
    CALL k2b (cr,ir,1)
    CALL k2b (cd,id,1)
    CALL k2b (cp,ip,1)
    CALL k2b (cm,im,1)
    CALL k2b (cg,ig,1)
    CALL k2b (ca,ia,1)
    CALL k2b (ch,ih,1)
    CALL k2b (ct,pt,1)
    CALL k2b (dotc,dot,1)
    CALL khrfn1 (univc(1),1,at,1)

50  IF (kount  /= 0) GO TO 300
60  IF (inflag == 0) IF (ffflag-1234) 80,200,80
    READ (inflag,65,END=150) (a8x(j),j=1,l12)
65  FORMAT (11A8,a6)
    !     NCARD = NCARD + 1
    IF (iechos == -2) WRITE (lout,65) a8x
    IF (a81 == rdfl) GO TO 4500
    IF (a81 == skfl .AND. a8(2) == BLANK) GO TO 130
    IF (ffflag == 1234) GO TO 240
    DO  i = 1,10
        card(i) = a8(i)
        SAVE(i) = a8(i)
    END DO
    GO TO 2800

    !     10A8 INPUT

80  READ (in,90,END=150) card
90  FORMAT (10A8)
    ncard = ncard + 1
    IF (iechos == -2) WRITE (lout,90) card

    IF (card(1) == skfl .AND. card(2) == BLANK) GO TO 130
    IF (card(1) /= rdfl) GO TO 2000
    DO  i = 1,10
        a8(i) = card(i)
    END DO
    CALL k2b (a8,a,80)
    GO TO 350

    !     IT IS A SKIPFILE CARD - TO SKIP TO THE END OF INPUT FILE

130 IF (inflag == 0) GO TO 5200
140 READ (inflag,90,END=5100) card
    GO TO 140

    !     CLOSE FILE, AND SET INFLAG BACK TO ZERO, OR PREVIOUS FILE OPENED

150 IF (mach >= 5) GO TO 154
    SELECT CASE ( mach )
        CASE (    1)
            GO TO 154
        CASE (    2)
            GO TO 154
        CASE (    3)
            GO TO 158
        CASE (    4)
            GO TO 152
    END SELECT
152 IF (inflag ==   0) RETURN 1
    IF (inflag >= iwo) REWIND inflag
    ierr = ierr + 1
    IF (ierr-15 > 0) THEN
        GO TO  3070
    ELSE
        GO TO   156
    ENDIF
154 IF (inflag ==   0) RETURN 1
156 CLOSE (UNIT=inflag)
158 IF (inflag == 0) RETURN 1
    inflag = inflag - 1
    IF (inflag <= iwo) inflag = 0
    card(1) = dend
    card(2) = rdfl
    DO  j = 3,10
        card(j) = BLANK
    END DO
    IF (iechos == -2) GO TO 60
    CALL page2 (-2)
    noecho = noecho - 1
    IF (noecho >= 0) WRITE (nout,165) noecho
165 FORMAT (12X,1H(,i4,' CARDS READ)')
    WRITE  (nout,460) card
    noecho = 0
    GO TO 60

170 loop  = 0
    loop4 = loop - 4
    kount = 0
    star  = .false.
    pct   = .false.
    notyet= .false.
    DO  j = 1,9
        l(j)  = 0
        f(j)  = 0.0
    END DO
    IF (inflag-iwo < 0) THEN
        GO TO   200
    ELSE
        GO TO    60
    ENDIF

    !     FREE FIELD INPUT

190 WRITE (nout,3020)
    ierr = ierr + 1
    IF (ierr > 3) GO TO 3070
    WRITE (screen,3060) a8
    IF (mach == 4 .AND. in == 5) REWIND in
200 IF (prom /= 0) WRITE (screen,210)
210 FORMAT (7H enter )
    READ   (in,220,END=190) (cx(j),j=1,l94)
220 FORMAT (94A1)
    ncard = ncard + 1
    lash  = 0
240 CONTINUE
    IF (iechos == -2) WRITE (lout,220) cx
    CALL k2b (a8,a,l94)
    IF (a1 == id) GO TO 280

    IF (a81 == rdfl) GO TO 350
    IF (a81 == skfl .AND. a8(2) == BLANK) GO TO 130
    IF (ffflag  == 1234) GO TO 260
    DO  i = 1,10
        card(i) = a8(i)
    END DO
    GO TO 2800
260 wasff = +1
    DO  i = 1,10
        IF (a(i) == ic .OR. a(i) == ie) GO TO 300
    END DO
280 wasff = -1
    IF (iechou == 0 .OR. xsort == 0) GO TO 288
    CALL page2 (-1)
    WRITE  (nout,285) a
285 FORMAT (30X,94A1)
288 IF (a1 == id) GO TO 60
    j = 0
    DO  i = 1,10
        IF (a8(i) /= BLANK) j = 1
        card(i) = a8(i)
    END DO
    loop  = -1
    loop4 = loop - 4
    IF (j == 0 .AND. iechos == -2) GO TO 4700
    GO TO 2000

300 IF (iechos == -2) GO TO 340
    IF (iechou == 0 .OR. kount >= 1) GO TO 320
    CALL page2 (-1)
    WRITE  (nout,310) a
310 FORMAT (30X,4H-ff-,4X,94A1)
320 IF (loop == -1) GO TO 340
    DO  j = 1,10
        card(j) = SAVE(j)
    END DO
340 IF (kount /= 0) GO TO 900
350 ke = 0
    k  = 0
    DO  j = 1,l94
        aii = a(j)
        IF (aii /= ib) GO TO 360
        IF (ke  ==  0) CYCLE
        IF (a(ke ) == ic .OR. a(ke ) == il) CYCLE
        IF (a(j+1) == ic .OR. a(j+1) == ib) CYCLE
        IF (a(j+1) == ir .AND.     k == 1) GO TO 370
        aii = ic
360     IF (aii == id) GO TO 390
        ke = ke + 1
        a(ke) = aii
        c(ke) = c(j)
        IF (aii == ic) c(ke) = cc
        IF (aii == il) k = k + 1
        IF (aii == ir) k = k - 1
        IF (k-1 > 0) THEN
            GO TO  5000
        ELSE
            GO TO   380
        ENDIF
370     k = 0
380 CONTINUE
    END DO
    IF (k  > 0) GO TO 5000
    IF (ke == 0) GO TO 4700
390 IF (a(ke) == ic) GO TO 400
    ke = ke + 1
    a(ke) = ic
    c(ke) = cc
400 IF (a81 /= rdfl) GO TO 520

    !     IT IS A READFILE CARD -
    !     CHECK NOPRINT OPTION, SET NOECHO = 1, IF FOUND.
    !     LOOK FOR FILE NAME. SET INFLAG TO UNIT IWO (OR IWO+ IF NESTED
    !     READFFILE), AND OPEN USERS FILE (NOT MEMBER OF A FILE AS IN IBM)

    !     READFILE FORMAT - '(', ')', ',', AND '=' ARE IGNORED.

    noecho = 0
    noec = 0
    i    = 9
405 a(1) = ib
    c(1) = cb
    c(8) = cc
    j = 0
410 i = i + 1
    IF (i > l94) GO TO 480
    aii = a(i)
    IF (aii == ib) GO TO 415
    IF (aii == il .OR. aii == ir .OR. aii == ic .OR. aii == ie) THEN
        IF (noec > 0) THEN
            GO TO   415
        ELSE
            GO TO   410
        ENDIF
    ENDIF
    j = j + 1
    IF (j > 48) GO TO 480
    a(j) = aii
    c(j) = c(i)
    IF (j /= 7 .OR. a81 /= noprt) GO TO 410
    noecho = 1
    noec = 1
    GO TO 405
415 IF (j ==  0) GO TO 410
    IF (j >= 60) GO TO 422
    j1 = j + 1
    DO  i = j1,60
        c(i) = cb
        a(i) = ib
    END DO
422 IF (mach == 3) GO TO 425
    IF (inflag < iwo) inflag = iwo - 1
    inflag = inflag + 1
    !WKBI 8/94 ALPHA-VMS
    IF ( mach == 21 ) GO TO 423
    IF (ibmcdc == 0) OPEN(UNIT=inflag,FILE=a8(1),STATUS='OLD',ERR=470)
    IF (ibmcdc /= 0) OPEN(UNIT=inflag,FILE=a48  ,STATUS='OLD',ERR=470)
    !WKBNB 8/94 ALPHA-VMS
    GO TO 424
423 indx = INDEX( a48, ' ' )
    a48(indx:indx) = '.'
    OPEN(UNIT=inflag,FILE=a48,STATUS='OLD',ERR=470)
424 CONTINUE
    !WKBNE 8/94 ALPHA-VMS

    IF (mach == 4) REWIND inflag
    GO TO 450

    !     UNIVAC - USE SYSTEM FACSF ROUTINE, SO THAT IT CAN READ A FILE OR
    !              AN ELEMENT OF A FILE. INPUT UNIT IWO IS NOT USED
    !              MAKE SURE FILE NAME CONTAINS A DOT
  
425 k = 0
    DO  i = 1,48
        IF (a(i) == dot) GO TO 440
        IF (a(i) /=  ib) k = 1
        IF (k == 1 .AND. a(i) == ib) GO TO 435
    END DO
    i = 49
435 a(i) = dot
440 inflag = in
    iwo    = in
    READ (a48,445) (univc(i),i=3,14)
445 FORMAT (12A4)
    i = facsf(univc)
    IF (i /= 0) GO TO 470
  
450 card(1) = dbgn
    card(2) = rdfl
    card(3) = from
    DO  j = 4,10
        card(j) = a8(j-3)
    END DO
    IF (iechos == -2) GO TO 465
    CALL page2 (-1)
    WRITE  (nout,460) card
460 FORMAT (5H0*** ,10A8)
    GO TO 60
465 prom = +1
    GO TO 60
  
470 WRITE  (nout,475) inflag,(a(i),i=1,j)
<<<<<<< HEAD
475 FORMAT (//,' *** CAN NOT OPEN FILE (UNIT=',i3,4H) - ,94A1)
    GO TO 500
480 j = j - 1
    WRITE  (nout,485) (a(i),i=1,j)
485 FORMAT (//,' *** FILE NAME ERROR - ',48A1)
    IF (j >= 48) WRITE (nout,490)
490 FORMAT (5X,'FILE NAME EXCEEDS 48 CHARACTERS')
500 nogo = 1
    IF (mach == 3 .OR. mach >= 5) WRITE (nout,505)
505 FORMAT (5X,'SUGGESTION- CHECK USER ID OR QUALIFIER')
=======
475 FORMAT (//,29H *** can NOT OPEN FILE (UNIT=,i3,4H) - ,94A1)
    GO TO 500
480 j = j - 1
    WRITE  (nout,485) (a(i),i=1,j)
485 FORMAT (//,23H *** FILE NAME error - ,48A1)
    IF (j >= 48) WRITE (nout,490)
490 FORMAT (5X,31HFILE NAME exceeds 48 characters)
500 nogo = 1
    IF (mach == 3 .OR. mach >= 5) WRITE (nout,505)
505 FORMAT (5X,38HSUGGESTION- check user id OR qualifier)
>>>>>>> d6697cb3550dda9a69b01e86a7002d41d2d4f575
    inflag = inflag - 1
    IF (inflag <= iwo) inflag = 0
    card(1) = BLANK
    card(2) = BLANK
    RETURN
  
    !     HERE WE GO
  
520 kk = 0
    ii = 0
    jj = 0
    twodot = .false.
530 iisave = jj - 2
540 jj = ii + 1
550 ii = ii + 1
    IF (ii > ke) GO TO 1500
    aii = a(ii)
    IF (aii == ih) GO TO 540
    IF (aii == ie) GO TO 700
    IF (jj  >  1) GO TO 580
    IF ((star .OR. pct) .AND. loop /= -1) WRITE (nout,560)
560 FORMAT (' *** PREVIOUS CARD SETTING UP FOR DUPLICATION IS NOW ',  &
        'ABANDONNED')
    kount = 0
    loop  = 0
    star  =.false.
    pct   =.false.
    notyet=.false.
    DO  j = 1,9
        l(j) = 0
        f(j) = 0.0
    END DO
580 IF (aii == ic) GO TO 600
    IF (aii == ia) GO TO 650
    IF (aii == ir) GO TO 1300
    IF (aii == is .OR. aii == ig) GO TO 1000
    IF (aii == il) GO TO 5400
    GO TO 550
  
    ! ... COMMA (,):
  
600 kk = kk + 1
    IF (kk == 1 .OR. kk == 10) GO TO 620
    je = ii - 1
    IF (je <= jj) GO TO 620
    i = 0
    DO  j = jj,je
        IF (a(j) == pt) i = i + 1
    END DO
    IF (i <= 1) GO TO 620
    IF (a5 /= seqgp .AND. a5 /= seqep) GO TO 4400
    twodot =.true.
    loop =-1
620 CALL nk12k8 (*3200,c(jj),ii-jj,card(kk),1)
    GO TO 530
  
    ! ... ECHO OR PROMPT:
  
630 CALL nk12k8 (*3200,c(jj),ii-jj,temp,1)
    IF (temp == cancel .OR. temp == list) GO TO 1600
    IF (temp4 ==   echo) GO TO 4600
    IF (temp4 /= prompt) GO TO 3000
    CALL nk12k8 (*3200,c(ii+1),4,temp,-1)
    IF (temp4 /= on .AND. temp4 /= off .AND. temp4 /= yes) GO TO 3000
    IF (temp4 == on ) prom =-1
    IF (temp4 == off) prom = 0
    IF (temp4 == yes) prom =+1
    GO TO 60
  
    ! ... SLASH (/):
  
650 IF (iisave <= 0) GO TO 660
    a(ii) = ih
    c(ii) = ch
    ii = ii + 1
    IF (a(ii) /= ic) GO TO 655
    a(ii) = ih
    c(ii) = ch
655 ii = iisave - 1
    GO TO 540
660 IF (lash == 0 .AND. kk == 0) GO TO 680
    j = kk + 1
    WRITE  (nout,670) j
<<<<<<< HEAD
670 FORMAT (' *** ILLEGAL USE OF SLASH IN FIELD',i3)
=======
670 FORMAT (34H *** illegal use of slash in field,i3)
>>>>>>> d6697cb3550dda9a69b01e86a7002d41d2d4f575
    GO TO 540
  
    !     A DELETE CARD (/) READ
  
680 lash = +1
    GO TO 530
  
    ! ... EQUAL (=):
  
700 IF (jj /= ii) GO TO 630
    kk = kk + 1
    ii = ii + 1
    IF (ii > ke) GO TO 3600
    aii = a(ii)
    IF (aii == il) GO TO 750
    IF (aii == ie) GO TO 730
    IF (aii == ic) GO TO 530
    GO TO 3600
  
730 kk = 10
    IF (twodot) GO TO 2400
    IF (loop > 0) THEN
        GO TO   850
    ELSE
        GO TO  2000
    ENDIF
  
    ! ... DUPLICATE WITH INCREMENT, =(N):
  
750 IF (kk /= 1) GO TO 3600
    jj = ii + 1
800 ii = ii + 1
    IF (ii  > ke) GO TO 3600
    aii = a(ii)
    IF (aii == ir) GO TO 820
    IF (aii == ic .OR. aii == is .OR. aii == ie) GO TO 3000
    GO TO 800
820 INT = 1
    CALL nk12if (*3800,c(jj),ii-jj,loop,INT)
    IF (loop <= 0) GO TO 4100
    loop4 = loop - 4
    ii = ii + 1
    IF (ii+1 < ke) GO TO 530
    IF (.NOT.star .AND. .NOT.pct) GO TO 3300
850 kount = 0
    IF (.NOT.notyet) GO TO 900
    notyet = .false.
    DO  kk = 2,9
        IF (l(kk) == NONE) GO TO 860
        IF (f(kk) /= xxxx) GO TO 870
        f(kk) = 0.0
        i = (l(kk)-jc(kk))/loop
        IF (i*loop+jc(kk) /= l(kk)) GO TO 4200
        l(kk) = i
        CYCLE
860     l(kk) = 0
        f(kk) = (f(kk)-rc(kk))/FLOAT(loop)
        CYCLE
870     IF (l(kk) /=   0) jc(kk) = jc(kk) - l(kk)
        IF (f(kk) /= 0.0) rc(kk) = rc(kk) - f(kk)
    END DO
900 kount = kount + 1
    IF (kount > loop) GO TO 170
    DO  kk = 2,9
        IF (l(kk) == 0) GO TO 920
        jc(kk) = jc(kk) + l(kk)
        CALL int2k8 (*3200,jc(kk),card(kk))
        CYCLE
920     IF (f(kk) == 0.0) CYCLE
        rc(kk) = rc(kk) + f(kk)
        CALL fp2k8 (*3000,rc(kk),card(kk))
    END DO
    IF (prom < 0 .AND. kount == loop) WRITE (screen,970) loop,card
970 FORMAT (/,i5,' ADDITIONAL CARDS WERE GENERATED.  LAST CARD WAS-',  &
        /1X,10A8)
    GO TO 2000
  
    ! ... STAR (*), OR PERCENTAGE (%):
  
1000 sp = aii
    ii = ii + 1
    IF (a(ii) /= il) GO TO 4000
    jj = ii + 1
    fp =.false.
    IF (star .OR. pct) GO TO 1030
    DO  k = 1,9
        l(k) = 0
        f(k) = 0.0
    END DO
1030 IF (sp == is) star =.true.
    IF (sp == ig) pct  =.true.
1050 ii = ii + 1
    aii= a(ii)
    IF (ii > ke .OR. aii == ic) GO TO 4000
    IF (aii == pt) fp =.true.
    IF (ii > jj .AND. (aii == ip .OR. aii == im)) fp =.true.
    IF (aii /= ir) GO TO 1050
    IF (ii  <= jj) GO TO 4000
    kk = kk + 1
    IF (fp) GO TO 1070
    INT = 1
    CALL nk12if (*3800,c(jj),ii-jj,l(kk),INT)
    CALL k82int (*3100,SAVE(kk),8,jc(kk),INT)
1060 IF (sp  == ig) GO TO 1120
    IF (loop > 0) GO TO 1100
    jc(kk) = jc(kk) + l(kk)
    CALL int2k8 (*3200,jc(kk),card(kk))
    GO TO 1100
1070 INT =-1
    CALL nk12if (*3900,c(jj),ii-jj,kkf,INT)
    f(kk) = fkk
    CALL k82fp  (*3100,SAVE(kk),8,rc(kk),INT)
1080 IF (sp  == ig) GO TO 1150
    IF (loop > 0) GO TO 1100
    rc(kk) = rc(kk) + f(kk)
    CALL fp2k8 (*3000,rc(kk),card(kk))
1100 ii = ii + 1
    GO TO 530
  
1120 IF (loop > 0) GO TO 1130
    f(kk) = xxxx
    GO TO 1160
1130 i = (l(kk)-jc(kk))/loop
    IF (i*loop+jc(kk) /= l(kk)) GO TO 4200
    l(kk) = i
    GO TO 1100
1150 IF (loop > 0) GO TO 1180
    l(kk)  = NONE
1160 notyet =.true.
    GO TO 1100
1180 f(kk) = (f(kk)-rc(kk))/FLOAT(loop)
    GO TO 1100
  
    ! ... RIGHT BRACKET ):
  
1300 IF (kk   ==  0) GO TO 1450
    IF (ii+1 >= ke) GO TO 3400
    aii = a(ii+1)
    IF (aii == is .OR. aii == ie) GO TO 3400
    j  = 10
    INT= 1
    IF (aii /= ip) CALL nk12if (*3900,c(jj),ii-jj,j,INT)
    IF (j <= 0 .OR. j > 10) GO TO 3700
    IF (j <= kk) GO TO 1400
    kk = kk + 1
    DO  k = kk,j
        card(k) = BLANK
    END DO
    kk = j
1400 IF (a(ii+1) == ic) ii = ii + 1
    jj = ii + 1
1420 ii = ii + 1
    IF (ii    > ke) GO TO 1430
    IF (a(ii) /= ic) GO TO 1420
1430 CALL nk12k8 (*3000,c(jj),ii-jj,card(j),1)
    IF (kk < 10) IF (ii-ke) 530,1500,1500
    GO TO 730
1450 kk = 1
    card(kk) = SAVE(10)
    ii = ii + 1
    IF (ii > ke .OR. a(ii) /= ic) GO TO 3000
    GO TO 530
  
    ! ... END OF CARD READ
  
1500 IF (kk-10 < 0) THEN
        GO TO  1550
    ELSE IF (kk-10 == 0) THEN
        GO TO   730
    ELSE
        GO TO  3500
    ENDIF
1550 kk = kk + 1
    card(kk) = BLANK
    IF (kk < 10) GO TO 1550
    GO TO 730
  
    ! ... CANCEL = N, LIST = +N
  
1600 IF (iechos /= -2) GO TO 5300
    card(1) = temp
    jj = ii + 1
1650 ii = ii + 1
    IF (a(ii) /= ic) GO TO 1650
    INT = 1
    CALL nk12if (*3800,c(jj),ii-jj,jc(1),INT)
    IF (temp == cancel .AND. jc(1) <= 0) GO TO 3800
    IF (temp ==  list .AND. jc(1) <= 0) GO TO 3800
    card(3) = temp
    GO TO 2800
  
    !     PREPARE TO RETURN
  
2000 IF (notyet) GO TO 60
  
    ! ... UPDATE CONTINUATION FIELDS IF WE ARE IN A DUPLICATION LOOP
  
    IF (loop == -1) GO TO 2400
    IF (kount == 0 .AND. .NOT.star) GO TO 2400
    kk = 10
    IF (SAVE(kk) == BLANK) GO TO 2300
2100 temp = SAVE(kk)
    IF (tmp(1) /= cp) GO TO 2300
    jj = 0
    DO  i = 3,8
        IF (tmp(i) == cm) jj = i
        IF (tmp(i) == cb) GO TO 2200
    END DO
    i = 9
2200 IF (jj == 0) GO TO 2300
    INT = 1
    CALL nk12if (*4800,tmp(jj+1),i-jj-1,j,INT)
    IF (mach == 3) GO TO 2230
    j = j + 1
    CALL int2k8 (*3800,j,tmp(jj+1))
    GO TO 2270
  
    ! ... UNIVAC USES NEXT 5 CARDS INSTEAD OF THE 3 ABOVE
  
2230 CALL int2k8 (*3800,j,spill)
    j = 9 - jj
    DO  i = 1,j
        tmp(jj+i) = tmp(8+i)
    END DO
2270 j = 9
    IF (tmp(j) /= cb) GO TO 4900
    card(kk) = temp
2300 IF (kk == 1) GO TO 2400
    kk = 1
    GO TO 2100
  
2400 IF (ffflag /= 1234) GO TO 2700
    IF (lash == +1) card(1) = slash
    IF (prom /= +1) GO TO 2500
    IF (kount < 7  .OR. kount > loop4) WRITE (screen,2450) card
    IF (kount == 7 .AND. kount <= loop4) WRITE (screen,2460)
2450 FORMAT (1X,10A8)
2460 FORMAT (9X,1H.,2(/,9X,1H.))
2500 IF (loop == -1) GO TO 2700
    DO  kk = 1,10
        SAVE(kk) = card(kk)
    END DO
2700 IF (card(1) == help .AND. card(2) == BLANK .AND. iechos == -2)  &
        CALL ffhelp (*60,*2900,2)
    IF (card(1) == STOP .AND. card(2) == BLANK .AND. iechos /= -2) GO TO 2900
    IF (card(1) /= scale8 .AND. card(1) /= scale1) GO TO 2800
    IF (card(1) == scale8) WRITE (nout,2710) (i,i=1,10)
    IF (card(1) == scale1) WRITE (nout,2720) (i,i=1,8 )
2710 FORMAT (/1X,10(i5,3X),/1X,5('--------++++++++'))
2720 FORMAT (/1X,     8I10,/1X,8('1234567890'))
    GO TO 60
  
2800 RETURN
2900 STOP
  
    !     ERRORS
  
3000 WRITE  (screen,3020)
3020 FORMAT (31H *** card error - INPUT ignored)
3050 IF (iechos == -2) GO TO 170
    IF (ierr   <= 15) WRITE (screen,3060) a8
3060 FORMAT (5X,1H',10A8,1H',/)
    nogo = 1
    ierr = ierr + 1
    IF (ierr < 30) GO TO 170
3070 WRITE  (screen,3080)
<<<<<<< HEAD
3080 FORMAT ('0*** JOB TERMINATED DUE TO TOO MANY INPUT ERRORS')
    STOP
3100 je = ii - 1
    WRITE  (screen,3150) kk,card(kk),(a(j),j=jj,je)
3150 FORMAT (5X,'FIELD',i3,' (',a8,') OF PREVIOUS CARD SHOULD NOT BE ',  &
=======
3080 FORMAT (48H0*** job terminated due TO too many INPUT errors)
    STOP
3100 je = ii - 1
    WRITE  (screen,3150) kk,card(kk),(a(j),j=jj,je)
3150 FORMAT (5X,5HFIELD,i3,2H (,a8,') OF PREVIOUS CARD SHOULD NOT BE ',  &
>>>>>>> d6697cb3550dda9a69b01e86a7002d41d2d4f575
        'USED FOR', /5X,'INCREMENTATION (BY ',8A1, ').  ZERO IS ASSUMED')
    IF (INT > 0) jc(kk) = 0
    IF (INT < 0) rc(kk) = 0.0
    IF (INT) 1080,3000,1060
3200 je = ii - 1
    WRITE  (screen,3250) kk,(a(j),j=jj,je)
3250 FORMAT (5X,'FIELD',i3,' IS TOO LONG. ONLY 8 DIGITS ALLOWED - ', 16A1)
    GO TO  3000
3300 WRITE  (screen,3350)
<<<<<<< HEAD
3350 FORMAT (5X,'PREVIOUS CARD WAS NOT SET UP FOR DUPLICATION')
    GO TO  3000
3400 WRITE  (screen,3450) a8
3450 FORMAT (' *** INDEX ERROR.  NO VALUE AFTER )')
    GO TO  3050
3500 WRITE  (screen,3550)
3550 FORMAT (' *** INPUT ERROR - TOO MANY FIELDS.  REPEAT INPUT')
    GO TO  3050
3600 WRITE  (screen,3650)
3650 FORMAT (' *** INPUT ERROR AFTER EQUAL SIGN (=)')
=======
3350 FORMAT (5X,44HPREVIOUS card was NOT set up for duplication)
    GO TO  3000
3400 WRITE  (screen,3450) a8
3450 FORMAT (35H *** INDEX error.  no value after ))
    GO TO  3050
3500 WRITE  (screen,3550)
3550 FORMAT (49H *** INPUT error - too many fields.  repeat INPUT)
    GO TO  3050
3600 WRITE  (screen,3650)
3650 FORMAT (37H *** INPUT error after equal SIGN (=))
>>>>>>> d6697cb3550dda9a69b01e86a7002d41d2d4f575
    IF (iechos == -2) GO TO 60
    WRITE  (screen,3060) a8
    nogo = 1
    GO TO  60
3700 WRITE  (screen,3750)
3750 FORMAT (5X,'INDEX ERROR BEFORE RIGHT BRACKET )')
    GO TO  3050
3800 je = ii - 1
    WRITE  (screen,3850) (a(j),j=jj,je)
<<<<<<< HEAD
3850 FORMAT (5X,'INVALID INTEGER - ',16A1)
    GO TO  3000
3900 je = ii - 1
    WRITE  (screen,3950) (a(j),j=jj,je)
3950 FORMAT (5X,'INVALID F.P. NUMBER - ',16A1)
    GO TO  3000
4000 WRITE  (screen,4050)
4050 FORMAT (' *** INPUT ERROR AFTER STAR (*), OR PERCENT (%)')
    GO TO  3050
4100 WRITE  (screen,4150)
4150 FORMAT (' *** ZERO LOOP COUNT.  NO CARDS GENERATED')
    GO TO  3050
4200 WRITE  (screen,4250) kk,l(kk),jc(kk),loop
4250 FORMAT (5X,'FIELD',i3,' (',i8,'-',i8,') IS NOT DIVIDABLE BY',i4,  &
        /5X,'RESUME INPUT',/)
=======
3850 FORMAT (5X,18HINVALID INTEGER - ,16A1)
    GO TO  3000
3900 je = ii - 1
    WRITE  (screen,3950) (a(j),j=jj,je)
3950 FORMAT (5X,22HINVALID f.p. NUMBER - ,16A1)
    GO TO  3000
4000 WRITE  (screen,4050)
4050 FORMAT (47H *** INPUT error after star (*), OR percent (%))
    GO TO  3050
4100 WRITE  (screen,4150)
4150 FORMAT (41H *** zero loop count.  no cards generated)
    GO TO  3050
4200 WRITE  (screen,4250) kk,l(kk),jc(kk),loop
4250 FORMAT (5X,5HFIELD,i3,2H (,i8,1H-,i8,21H) is NOT dividable by,i4,  &
        /5X,12HRESUME INPUT,/)
>>>>>>> d6697cb3550dda9a69b01e86a7002d41d2d4f575
4300 IF (iechos /= -2) nogo = 1
    DO  j = 1,10
        card(j) = SAVE(j)
    END DO
    GO TO  60
4400 WRITE  (screen,4450) (a(j),j=jj,je)
<<<<<<< HEAD
4450 FORMAT (5X,'MORE THAN ONE DEC. PT.,  - ',16A1)
    GO TO  3000
4500 WRITE  (screen,4550)
4550 FORMAT (' *** WARNING- NESTED READFILE OPERATION')
    GO TO  350
4600 WRITE  (screen,4650)
4650 FORMAT (' *** SO BE IT.  TO RUN NASTRAN LINK1 ONLY ***',/)
    GO TO  60
4700 WRITE  (screen,4750)
4750 FORMAT (' *** BLANK LINE IGNORED')
    GO TO  60
4800 WRITE  (screen,4850) temp
4850 FORMAT (' *** INTEGER ERROR IN CONTINUATION ID - ',a8)
    IF (iechos /= -2) WRITE (screen,3060) a8
    GO TO  4300
4900 WRITE  (screen,4950) (tmp(j),j=1,9)
4950 FORMAT (' *** CONTINUATION FIELD TOO LONG - ',9A1, /5X,  &
        'LAST GENERATED CARD WAS -',/)
    WRITE  (screen,2450) SAVE
    GO TO  4300
5000 WRITE  (screen,5050)
5050 FORMAT (' *** TOO MANY LEFT BRACKETS')
    GO TO  3050
5100 WRITE  (nout,5150)
5150 FORMAT (/,' *** EOF ENCOUNTERED')
    IF (mach == 4 .AND. inflag == 5) REWIND inflag
    GO TO  60
5200 WRITE  (nout,5250)
5250 FORMAT (/,' *** SKIPFILE IGNORED.  FILE HAS NOT BEEN OPENED')
    GO TO  60
5300 WRITE  (nout,5350)
5350 FORMAT (/,' *** FEATURE NOT AVAILABLE')
    IF (iechos /= -2) WRITE (screen,3060) a8
    GO TO  60
5400 WRITE  (nout,5450)
5450 FORMAT (/,' *** LEFT BRACKET ENCOUNTERED WITHOUT FIRST PRECEEDED BY ', &
               "'=', '*', OR '%'")
=======
4450 FORMAT (5X,27HMORE than one dec. pt.,  - ,16A1)
    GO TO  3000
4500 WRITE  (screen,4550)
4550 FORMAT (39H *** warning- nested readfile operation)
    GO TO  350
4600 WRITE  (screen,4650)
4650 FORMAT (45H *** so be it.  TO run nastran link1 only ***,/)
    GO TO  60
4700 WRITE  (screen,4750)
4750 FORMAT (23H *** BLANK line ignored)
    GO TO  60
4800 WRITE  (screen,4850) temp
4850 FORMAT (40H *** INTEGER error in continuation id - ,a8)
    IF (iechos /= -2) WRITE (screen,3060) a8
    GO TO  4300
4900 WRITE  (screen,4950) (tmp(j),j=1,9)
4950 FORMAT (35H *** continuation field too long - ,9A1, /5X,  &
        25HLAST generated card was -,/)
    WRITE  (screen,2450) SAVE
    GO TO  4300
5000 WRITE  (screen,5050)
5050 FORMAT (27H *** too many left brackets)
    GO TO  3050
5100 WRITE  (nout,5150)
5150 FORMAT (/,20H *** eof encountered )
    IF (mach == 4 .AND. inflag == 5) REWIND inflag
    GO TO  60
5200 WRITE  (nout,5250)
5250 FORMAT (/,48H *** skipfile ignored.  FILE has NOT been OPENED)
    GO TO  60
5300 WRITE  (nout,5350)
5350 FORMAT (/,26H *** feature NOT available)
    IF (iechos /= -2) WRITE (screen,3060) a8
    GO TO  60
5400 WRITE  (nout,5450)
5450 FORMAT (/,73H *** left bracket encountered without first preceeded &
        &by '=', '*', OR '%')
>>>>>>> d6697cb3550dda9a69b01e86a7002d41d2d4f575
    GO TO  3000
  
END SUBROUTINE ffread
