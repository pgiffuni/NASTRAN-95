SUBROUTINE xcsa
     
    !     XCSA READS AND PROCESSES THE NASTRAN EXECUTIVE CONTROL DECK.
 
    IMPLICIT INTEGER (a-z)
    EXTERNAL        lshift,rshift,andf,orf,complf
    LOGICAL :: tapbit

    DIMENSION alter(2),apptyp(4),bgnal(2),cend(2),diagx(11), &
              dmapbf(1),ectt(51),endal(2),hdg(19),iptdic(1), &
              iufile(2),iz(2),nxptdc(2),nxcsa(2),osolu(2),   &
              outcrd(200),solrec(6),solu(12),solnm3(7,11),   &
              solnms(7,31),solnm1(7,10),solnm2(7,10),        &
              solnmx(6), xalt(2),xsys(100)

    INTEGER :: insert(4), DELETE(9), altrbs, alnogo
    INTEGER :: altfil, erralt, altopn
    CHARACTER (LEN=25) :: sfm
    CHARACTER (LEN=29) :: uim
    CHARACTER (LEN=25) :: uwm
    CHARACTER (LEN=23) :: ufm

    COMMON /xmssg / ufm,uwm,uim,sfm
    COMMON /machin/ mach,ijhalf(3),mchnam
    COMMON /sem   / mskdum(3),links(15)
    COMMON /system/ ibufsz,outtap,xnogo,intape,sy5,sy6,logfl,sy8,   &
                    nlpp,sy10,sy11,nlines,sy13,sy14,idate(3),sy18,  &
                    iecho,sy20,apprch,sy22,sy23,icfiat,rfflag,      &
                    sy26(11),lu,sy38,nbpc,nbpw,ncpw,sy42(13),prec,  &
                    sy56(13),isubs,sy70(9),switch(3),icpflg,sy83(2),&
                    sy85,intra,sy87(5),ldict
    COMMON /xechox/ dum9(9),noecho
    COMMON /xrgdxx/ irestr,nsubst
    COMMON /altrxx/ altfil, newalt,alnogo
    COMMON /resdic/ irdict, iropen
    COMMON /xoldpt/ itop,ibot,ldic,nrlfl,iseqno
    COMMON /xxfiat/ ixxfat(1)
    COMMON /xpfist/ ixpfst
    COMMON /xfist / ifist(1)
    COMMON /xfiat / ifiat(1)
    COMMON /zzzzzz/ gbuff(1)
    COMMON /BLANK / zcom,card(20)
    COMMON /stapid/ tapid(6),otapid(6)
    COMMON /stime / time
    COMMON /l15_l8/ l15,l8,l13
    COMMON /xlink / lxlink,maxlnk
    COMMON /output/ pghdg1(32),pghdg2(32), pghdg3(32),  &
                    pghdg4(32),pghdg5(32), pghdg6(32)

    EQUIVALENCE     (ibufsz   ,xsys(1)  ), (mask    ,maskhi  ),  &
                    (ectt(16) ,bgnal(1) ), (ectt(25),endal(1)),  &
                    (ectt(13) ,cend(1)  ), (ectt(34),id      ),  &
                    (solrec(1),apprec   ), (solrec(2),rstrt  ),  &
                    (solrec(3),alter(1) ), (solrec(5),solu(1)),  &
                    (gbuff(1) ,dmapbf( 1), iptdic(1)), (solnms(1, 1)&
                    ,solnm1(1,1)),  &
                    (solnms(1,11),solnm2(1,1)), (solnms(1,21),solnm3(1,1))

    DATA  apptyp                                                   /  &
          4HDMAP,   4HDISP,   4HHEAT,   4HAERO                     /
    DATA  BLANK,    ixdmap,   nsubs,    renter,   dolsin           /  &
          1H ,      4HXDMA,   4HSUBS,   4HREEN,   4H$              /
    DATA  iyes,     no,       idisk,    ptape,    optape,   dmend  /  &
          4HYES ,   4HNO  ,   4HDISK,   4HNPTP,   4HOPTP,   4HEND  /
    DATA  iufile,   xalt,               nxptdc,             intgr  /  &
          2*0,      4HXALT,   4HER  ,   4HXPTD,   4HC   ,   -1     /
    DATA  nxcsa,              diagx                                /  &
          4HXCSA,   4H    ,   4,9,14,17,23,24,25,28,29,30,31       /
    DATA  appdmp,   apphea,   appaer,   numapp,   solrec           /  &
          1,        3,        4,        4,        0,1,0,0,0,0      /
    DATA  soluf,    osolu,    icold,    ignore,   outcrd           /  &
          0,        2*0,      1,        0,        3,199*4H         /
    DATA  plot,     prnt,     both,     inp9  ,   notalt           /  &
          4HPLOT,   4HPRIN,   4HBOTH,   4HINP9,   0                /
    DATA  mask /    32767 /
    !                     32767 = O77777 = 2**15-1 = MASK HI
    DATA  lectt,    ectt /    51,  &
          4HTIME,4H    ,0   ,   4HAPP ,4H    ,0   ,   4HCHKP,4HNT  ,0,  &
          4HREST,4HART ,0   ,   4HCEND,4H    ,0   ,   4HALTE,4HR   ,0,  &
          4HSOL ,4H    ,0   ,   4HBEGI,4HN   ,0   ,   4HENDA,4HLTER,0,  &
          4HDIAG,4H    ,0   ,   4HUMF ,4H    ,0   ,   4HID  ,4H    ,1,  &
          4HUMFE,4HDIT ,0   ,   4HPREC,4H    ,0   ,   4HINTE,4HRACT,0,  &
          4HINSE,4HRT  ,0   ,   4HDELE,4HTE  ,0/
    DATA ileft  /4H(    /
    DATA altopn / 0     /
    DATA hdg/4HN A ,4HS T ,4HR A ,4HN   ,4H E X,4H E C,4H U T,4H I V,  &
        4H E  ,4H  C ,4HO N ,4HT R ,4HO L ,4H   D,4H E C,4H K  ,4H  E , 4HC H ,4HO   /
    DATA nsolnm /26/

    DATA solnm1 / 4HSTAT,4HICS  , 4H    ,4H     , 4H    ,4H     ,  1 ,  &
                  4HINER,4HTIA  , 4HRELI,4HEF   , 4H    ,4H     ,  2 ,  &
                  4HNORM,4HAL   , 4HMODE,4HS    , 4H    ,4H     ,  3 ,  &
                  4HDIFF,4HEREN , 4HSTIF,4HFNES , 4H    ,4H     ,  4 ,  &
                  4HBUCK,4HLING , 4H    ,4H     , 4H    ,4H     ,  5 ,  &
                  4HPIEC,4HEWIS , 4HLINE,4HAR   , 4H    ,4H     ,  6 ,  &
                  4HDIRE,4HCT   , 4HCOMP,4HLEX  , 4HEIGE,4HNVAL ,  7 ,  &
                  4HDIRE,4HCT   , 4HFREQ,4HUENC , 4HRESP,4HONSE ,  8 ,  &
                  4HDIRE,4HCT   , 4HTRAN,4HSIEN , 4HRESP,4HONSE ,  9 ,  &
                  4HMODA,4HL    , 4HCOMP,4HLEX  , 4HEIGE,4HNVAL , 10 /
    DATA solnm2 / 4HMODA,4HL    , 4HFREQ,4HUENC , 4HRESP,4HONSE , 11 ,  &
                  4HMODA,4HL    , 4HTRAN,4HSIEN , 4HRESP,4HONSE , 12 ,  &
                  4HSTEA,4HDY   , 4HSTAT,4HE    , 4H    ,4H     ,  3 ,  &
                  4HTRAN,4HSIEN , 4H    ,4H     , 4H    ,4H     ,  9 ,  &
                  4HMODE,4HS    , 4H    ,4H     , 4H    ,4H     ,  3 ,  &
                  4HREAL,4H     , 4HEIGE,4HNVAL , 4H    ,4H     ,  3 ,  &
                  4HMODA,4HL    , 4HFLUT,4HTER  , 4HANAL,4HYSIS , 10 ,  &
                  4HMODA,4HL    , 4HAERO,4HELAS , 4HRESP,4HONSE , 11 ,  &
                  4HNORM,4HAL   , 4HMODE,4HS    , 4HANAL,4HYSIS , 13 ,  &
                  4HSTAT,4HICS  , 4HCYCL,4HIC   , 4HSYMM,4HETRY , 14 /
    DATA solnm3 / 4HMODE,4HS    , 4HCYCL,4HIC   , 4HSYMM,4HETRY , 15 ,  &
                  4HSTAT,4HIC   , 4HAERO,4HTHER , 4HMOEL,4HASTI , 16 ,  &
                  4HBLAD,4HE    , 4HCYCL,4HIC   , 4HMODA,4HL    ,  9 ,  &
                  4HDYNA,4HMIC  , 4HDESI,4HGN A , 4HNALY,4HSIS  , 17 ,  &
                  4HDIRE,4HCT   , 4HFORC,4HED V , 4HIBRA,4HTION , 18 ,  &
                  4HMODA,4HAL   , 4HFORC,4HED V , 4HIBRA,4HTION , 19 ,  &
                  4H****,4H**** , 4H****,4H**** , 4H****,4H**** ,  0 ,  &
                  4H****,4H**** , 4H****,4H**** , 4H****,4H**** ,  0 ,  &
                  4H****,4H**** , 4H****,4H**** , 4H****,4H**** ,  0 ,  &
                  4H****,4H**** , 4H****,4H**** , 4H****,4H**** ,  0 ,  &
                  4H****,4H**** , 4H****,4H**** , 4H****,4H**** ,  0 /

    !     SET UP DATA IN COMMON

    itop   = 0
    ibot   = 0
    ldic   = 0
    nrlfl  = 0
    iseqno = 0
    altfil = 301
    newalt = 0
    alnogo = 0
    erralt = 0
    nscr   = 315
    irestr = 0
    nsubst = 0
    nwpc   = 18
    drecsz = 0


    !     INITIALIZE MACHINE DEPENDENT CONSTANTS

    !     ALLON  = O777777777777  ALL BITS ON
    !     ISIGN  = O400000000000  SIGN ON ONLY
    !     MASK5  = O500000000000  SIGN AND NEXT BIT ON
    !     ENDCD  = O377777777777  ALL BITS ON EXCEPT SIGN
    !     MHIBYT = O770000000000  MASK IN HIGH ORDER BYTE

    ISIGN  = lshift(1,nbpw-1)
    mask5  = orf(ISIGN,rshift(ISIGN,1))
    allon  = complf(0)
    mhibyt = lshift(allon,(ncpw-1)*nbpc)
    endcd  = rshift(allon,1)
    j      = diagx(2)*5 - 1
    card(j  ) = xsys(j)
    card(j+1) = khrfn1(bnk,1,xsys(j),2)
    CALL na12if (*1420,card(j),2,s7,1)
    IF (s7 /= 0) i7 = mach*100

    !     DETERMINE OPEN CORE SIZE AND ALLOCATE BUFFER AREA

    dmapbs = korsz(gbuff) - 2*ibufsz
    altrbs = dmapbs + ibufsz
    CALL waltim (timew)
    timew = MOD(timew,10000000)

    !     LOAD PAGE HEADING IN /OUTPUT/

    j = 32
    DO  i = 1,j
        pghdg1(i) = BLANK
        pghdg2(i) = BLANK
        pghdg3(i) = BLANK
        pghdg4(i) = BLANK
        pghdg5(i) = BLANK
        pghdg6(i) = BLANK
    END DO
    DO  i = 1,19
        pghdg3(i+1) = hdg(i)
    END DO
    CALL page

    !     CARD PREPARATION

    n7 = i7 + s7
    i7 = i7/100
    n7 = n7 - 2*i7
    m7 = card(lectt+9)
    j  = IABS(m7)
    i  = 3
    IF (m7 < 0 .AND. MOD(j,10) == 7) i = 4
    IF (j/10 == n7 .AND. xsys(17)-i <= s7) card(lectt+2) = icold
    card(lectt+11) = khrfn1(card(lectt+11),2,xalt(1),3)
    card(lectt+13) = khrfn1(card(lectt+13),1,nxcsa(1),1)
    card(lectt+14) = khrfn1(card(lectt+14),2,idisk,1)

    !     WRITE DUMMY ID FILE ON PROBLEM TAPE IN CASE OF ID CONTROL CARD
    !     ERROR.

    nogo   = xnogo
    xnogo  = 0
    oldalt = 0

    !     READ CONTROL CARD AND PROCESS

20  IF (altopn <= 0) ASSIGN 70 TO irtn1
30  nlines = nlines + 1
    IF (nlines >= nlpp) CALL page
    IF (zcom /= 0) GO TO 40
    CALL xread (*1232,card)

    !     ECHO CARD
    !     (NOECHO IS SET BY SEMDBD AND READFILE OF FFREAD)

40  zcom = 0
    IF (noecho /= 0) GO TO 52
    WRITE  (outtap,50) card
50  FORMAT (5X,20A4)
    GO TO 55
52  noecho = noecho + 1
    nlines = nlines - 1

    !     CHECK FOR COMMENT CARD

55  IF (khrfn1(BLANK,1,card(1),1) == dolsin) GO TO 30

    !     CALL RMVEQ TO REPLACE ONE EQUAL SIGN BY ONE BLANK
    !     IF CARD IS NOT WITHIN ALTER RANGE

    !CCCC   NEXT LINE CAUSE ERROR IN READING RESTART DICTIONARY. POSITION
    !CCCC   PROBLEM
    !CCCC
    !CCCC      IF (NOTALT .EQ. 0) CALL RMVEQ (CARD)
    CALL xrcard (outcrd,200,card)

    !     CHECK FOR ERROR DETECTED BY XRCARD

    IF (xnogo == 0) GO TO 60
    IF (nogo  == 0) nogo = 1
    xnogo = 0
    GO TO 30

    !     CHECK FOR BLANK CARD

60  IF (outcrd(1) == 0) GO TO 30
    GO TO irtn1, (70,270,370,510)
70  j = 0
    DO  i = 1,lectt,3
        j = j + 1
        IF (outcrd(2) == ectt(i) .AND. outcrd(3) == ectt(i+1)) GO TO 90
    END DO
    IF (outcrd(2) == ixdmap) GO TO 400
    IF (ignore == 0) GO TO 690
    GO TO 20

    !     HAS THIS TYPE CARD ALREADY BEEN PROCESSED

90  ignore = 0
    IF (ectt(i+2) < 0 .AND. outcrd(2) == ectt(28)) ectt(i+2) = 0
    !                                               DIAG
    IF (ectt(i+2) < 0.0) THEN
        GO TO   720
    END IF
100 ectt(i+2) = orf(ectt(i+2),mask5)
    SELECT CASE ( j )
        CASE (    1)
            GO TO 110
        CASE (    2)
            GO TO  120
        CASE (    3)
            GO TO  140
        CASE (    4)
            GO TO  210
        CASE (    5)
            GO TO  570
        CASE (    6)
            GO TO  330
        CASE (    7)
            GO TO  390
        CASE (    8)
            GO TO  400
        CASE (    9)
            GO TO 1180
        CASE (   10)
            GO TO  480
        CASE (   11)
            GO TO 460
        CASE (   12)
            GO TO  530
        CASE (   13)
            GO TO  560
        CASE (   14)
            GO TO  565
        CASE (   15)
            GO TO  555
        CASE (   16)
            GO TO  330
        CASE (   17)
            GO TO  330
    END SELECT


    !     NOW PROCESS TIME CARD

110 imhere = 110
    IF (outcrd(4) /= -1 .OR. outcrd(5) <= 0) GO TO 760
    time = outcrd(5)*60
    GO TO 20


    !     NOW PROCESS APPROACH CARD

    120 DO  jj = 1,numapp
        apprch = jj
        apprec = jj
        IF (outcrd(4) == apptyp(jj)) GO TO 132
    END DO
    imhere = 130
    GO TO 760

    !     CHECK FOR SUBSTRUCTURE ANALYSIS

132 IF (outcrd(6) /= nsubs) GO TO 20
    isubs = apprch
    IF (outcrd(8) /= -1) GO TO 20
    isubs = isubs + 10*outcrd(9)
    GO TO 20


    !     NOW PROCESS CHKPNT CARD

140 IF (outcrd(4) == no .OR. outcrd(6) == no) GO TO 20

    !     CHECK FOR ILLEGAL FORMAT

    imhere = 140
    IF (outcrd(4) /= iyes .AND. outcrd(6) /= iyes) GO TO 750
    icpflg = 1
    IF (outcrd(6) == idisk) GO TO 20
    ASSIGN 150 TO l
    idfist = ptape

    !     CHECKPOINT FLAG IS ON,MAKE SURE NEW PROBLEM TAPE IS ON
    !     PHYSICAL TAPE DRIVE

    GO TO 160
150 IF (nostup /= 0) GO TO 790
    GO TO 20

    !     CHECK TAPE SETUP

160 IF (tapbit(idfist)) GO TO 190

    !     TAPE NOT SETUP

    nostup = 1
    GO TO 200
190 CONTINUE

    !     TAPE SETUP

    nostup = 0
    !     GO TO L, (150,470)
200 GO TO l, (150)


    !     NOW PROCESS RESTART CARD

210 ngino  = optape
    irestr = 1

    !     SET UNSORTED AND SORTED BULK DATA OUTPUT (ECHO = BOTH)
    !     AS THE DEFAULT FOR RESTART RUNS

    iecho = 3
    CALL OPEN (*850,optape,gbuff(dmapbs+1),0)
    CALL READ (*1350,*1350,optape,otapid,6,0,flgwrd)
    CALL READ (*1350,*222,optape,timex,1,1,flgwrd)
    GO TO 225
222 outcrd(21) = 0
    timex = 0

    !     COMPARE ID OF OLD PTAPE WITH THAT ON RSTART CARD

225 rstrt = 2

    !     UNPACK DATE

    i     = lshift(otapid(5),7)
    iyear = rshift(andf(i,maskhi),7)
    i     = rshift(i,6)
    iday  = rshift(andf(i,maskhi),9)
    i     = rshift(i,5)
    imnth = rshift(andf(i,maskhi),10)
    jj    = outcrd(1)*2 - 2
    DO  jk = 1,jj
        IF (otapid(jk) /= outcrd(jk+3)) GO TO 820
    END DO
    IF (outcrd( 9) == 0 .AND. outcrd(14) == 0 .AND. outcrd(19) == 0) GO TO 235
    IF (imnth /= outcrd(9) .OR. iday /= outcrd(14) .OR.  &
        iyear /= outcrd(19)) GO TO 820
235 CONTINUE
    IF (outcrd(21) == 0) timex = 0
    IF (timex /= outcrd(21)) GO TO 820

    !     MAKE SURE CORRCET REEL IS MOUNTED

    IF (otapid(6) == 1) GO TO 240
    GO TO 820

    !     GET OLD SOLUTION NUMBER

240 CALL skpfil (optape,1)
    CALL READ  (*1350,*1350,optape,osolu,1,0,flgwrd)
    IF (osolu(1) == xalt(1)  ) oldalt = oldalt + 1
    IF (osolu(1) == nxptdc(1)) oldalt = oldalt + 1
    IF (osolu(1) /= nxcsa(1) ) GO TO 240
    CALL fwdrec (*1350,optape)
    CALL READ (*1350,*1350,optape,0,-4,0,flgwrd)
    CALL READ (*1350,*1350,optape,osolu,2,1,flgwrd)
    CALL skpfil (optape,1)
    CALL CLOSE  (optape,2)

    !     LOAD PROBLEM TAPE DICTIONARY

    icrdct = 0
    iseqno = 0
    itop = drecsz + 1
    ldic = korsz(iptdic(itop)) - ibufsz
    ibot = itop - 3

    !     ZERO FIRST PTDIC ENTRY IN CASE THERE ARE NO ENTRIES

    iptdic(itop  ) = 0
    iptdic(itop+1) = 0
    iptdic(itop+2) = 0

    !     SET ITOPX SO THAT FIRST XVPS ENTRY IN PTDIC WILL BE PRESERVED

    itopx  = itop + 3
260 icrdct = 1 + icrdct

    !     READ IN NEXT CONTROL CARD

    ASSIGN 270 TO irtn1
    GO TO 30
270 IF (outcrd(1) /=     -1) GO TO 320
    IF (outcrd(2) /= icrdct) GO TO 1210
    IF (outcrd(3) ==      5) GO TO 310
    IF (outcrd(3) ==  endcd) GO TO 320
    IF (outcrd(3) >      3) GO TO 310

    !     CHECK FORMAT

    imhere = 275
    IF (outcrd(3) /= 3 .OR. outcrd(10) /= -1 .OR. outcrd(12) /= 2 .OR.  &
        outcrd(17) /= -1 .OR. outcrd(19) /= 2 .OR. outcrd(24) /= -1) GO TO 760

    !     PACK FLAGS/REEL/FILE

    flags = 0
    IF (outcrd(11) >= 4) flags = ISIGN
    reel = orf(lshift(outcrd(18),16),outcrd(25))

    !     SEE IF FILE IS ALREADY IN PTDIC - IF IT IS, PUT LATEST REEL/FILE
    !     NO. IN EXISTING ENTRY

    IF (ibot < itopx) GO TO 290
    DO  k = itopx,ibot,3
        IF (iptdic(k) == outcrd(4) .AND. iptdic(k+1) == outcrd(5)) GO TO 300
    END DO

    !     FILE NOT IN PTDIC - MAKE NEW ENTRY

290 ibot = ibot + 3

    !     CHECK FOR OVERFLOW

    IF (ibot+3-itop > ldic) GO TO 1260
    k = ibot
    iptdic(k  ) = outcrd(4)
    iptdic(k+1) = outcrd(5)
300 iptdic(k+2) = orf(flags,reel)
    GO TO 260

    !     THIS IS A REENTRY CARD - LOAD DMAP INSTRUCTION NO. IN ISEQNO

310 imhere = 310
    IF (outcrd(4) /= renter .OR. outcrd(14) /= -1) GO TO 760
    iseqno = lshift(outcrd(15),16)
    GO TO 260

    !     DICTIONARY PROCESSED - COPY ONTO NEW PROBLEM TAPE.
    !     THERE MUST ALWAYS BE AT LEAST ONE ENTRY IN PTDIC

320 IF (ibot < itop) ibot = itop
    ngino = ptape
    imhere= 320
    CALL OPEN (*1320,ptape,gbuff(dmapbs+1),3)

    !     RECORD 1 = ID

    CALL WRITE (ptape,nxptdc,2,1)

    !     RECORD 2 = CONTENTS OF IPTDIC

    CALL WRITE (ptape,iptdic(itop),ibot+3-itop,1)
    CALL eof   (ptape)
    CALL CLOSE (ptape,2)
    IF (outcrd(3) == endcd) GO TO 20
    GO TO 70


    !     PROCESS ALTER CONTROL CARDS

330 ASSIGN 370 TO irtn1
    IF (ectt(27) < 0) GO TO 30
    notalt = 1
    imhere = 330
    ngino = altfil
    CALL OPEN (*1320,altfil,gbuff(altrbs+1),1)
    altopn = 1
    IF (j == 16) GO TO 3605
    IF (j == 17) GO TO 3655
340 IF (outcrd(6) /= endcd) GO TO 350
    outcrd(6) = intgr
    outcrd(7) = 0
350 imhere = 350
    IF (outcrd(4) /= intgr .OR. outcrd(6) /= intgr .OR. outcrd(5) <= 0  &
        .OR. outcrd(7) < 0) GO TO 750
    IF (outcrd(7) > 0 .AND. outcrd(8) /= endcd) GO TO 750


    alter(1) = outcrd(5)
    alter(2) = outcrd(7)


    !     WRITE ALTER PARAMETERS ONTO THE ALTER SCRATCH FILE
    !     AND FOLLOW IT BY THE CARD IMAGE

    CALL WRITE (altfil, alter,  2, 1)
    CALL WRITE (altfil, card , 18, 1)

    !     READ NEXT CARD INTO CORE

    GO TO 30

    !     PROCESS INSERT CONTROL CARDS HERE

3605 insert(1) = outcrd(4)
    insert(2) = outcrd(5)
    insert(3) = 1
    insert(4) = 0
    IF (outcrd(6) == allon .AND. outcrd(7) == ileft .AND.  &
        outcrd(8) == intgr) GO TO 3610
    jn = 7
    IF (outcrd(6) == intgr) GO TO 3615
    IF (outcrd(6) == endcd) GO TO 3620
    GO TO 750
3610 IF (outcrd(9) <= 0) GO TO 750
    insert(3) = outcrd(9)
    jn = 11
    IF (outcrd(10) == intgr) GO TO 3615
    IF (outcrd(10) == endcd) GO TO 3620
    GO TO 750
3615 insert(4) = outcrd(jn)
    IF (outcrd(jn+1) /= endcd) GO TO 750

    !     WRITE INSERT PARAMETERS ONTO THE ALTER SCRATCH FILE
    !     AND FOLLOW IT BY THE CARD IMAGE

3620 CALL WRITE (altfil, insert,  4, 1)
    CALL WRITE (altfil, card  , 18, 1)
    newalt = 1
    GO TO 30

    !     PROCESS DELETE CONTROL CARDS HERE

3655 DELETE(1) = outcrd(4)
    DELETE(2) = outcrd(5)
    DELETE(3) = 1
    DELETE(4) = 0
    DELETE(5) = 0
    IF (outcrd(6) == allon .AND. outcrd(7) == ileft .AND.  &
        outcrd(8) == intgr) GO TO 3660
    jn = 7
    jnx = 7
    IF (outcrd(6) == intgr) GO TO 3665
    IF (outcrd(6) == endcd) GO TO 3670
    jnx = 6
    GO TO 3675
3660 IF (outcrd(9) <= 0) GO TO 750
    DELETE(3) = outcrd(9)
    jn = 11
    jnx = 11
    IF (outcrd(10) == intgr) GO TO 3665
    IF (outcrd(10) == endcd) GO TO 3670
    IF (outcrd(10) >     0) GO TO 3675
    GO TO 750
3665 DELETE(4) = outcrd(jn)
    jn = jn + 1
    jnx = jn + 1
    IF (outcrd(jn) == endcd) GO TO 3670
    IF (outcrd(jn) >     0) GO TO 3675
    GO TO 750

    !     WRITE DELETE PARAMETERS ONTO THE ALTER SCRATCH FILE
    !     AND FOLLOW IT BY THE CARD IMAGE

3670 CALL WRITE (altfil, DELETE,  5, 1)
    CALL WRITE (altfil, card  , 18, 1)
    newalt = 1
    GO TO 30

3675 jn = jnx
    DELETE(5) = 1
    DELETE(6) = outcrd(jn  )
    DELETE(7) = outcrd(jn+1)
    DELETE(8) = 1
    DELETE(9) = 0
    jn = jn + 2
    jnx = jn + 3
    IF (outcrd(jn  ) == allon .AND. outcrd(jn+1) == ileft .AND.  &
        outcrd(jn+2) == intgr) GO TO 3680
    jnx = jn + 1
    IF (outcrd(jn) == intgr) GO TO 3685
    IF (outcrd(jn) == endcd) GO TO 3690
    GO TO 750
3680 jn = jnx
    IF (outcrd(jn) <= 0) GO TO 750
    DELETE(8) = outcrd(jn)
    jn = jn + 1
    jnx = jn + 1
    IF (outcrd(jn) == intgr) GO TO 3685
    IF (outcrd(jn) == endcd) GO TO 3690
    GO TO 750
3685 DELETE(9) = outcrd(jnx)
    IF (outcrd(jnx+1) /= endcd) GO TO 750

    !     WRITE DELETE PARAMETERS ONTO THE ALTER SCRATCH FILE
    !     AND FOLLOW IT BY THE CARD IMAGE

3690 CALL WRITE (altfil, DELETE,  9, 1)
    CALL WRITE (altfil, card  , 18, 1)
    newalt = 1
    GO TO 30
370 CONTINUE

    !     CHECK FOR CEND CARD TO PREVENT STREAMING THRU BULK DATA

    IF (outcrd(2) == cend(1) .AND. outcrd(3) == cend(2)) GO TO 910

    !     CHECK FOR ANOTHER ALTER CARD

    IF (outcrd(2) == bgnal(1) .AND. outcrd(3) == bgnal(2)) GO TO 340

    !     CHECK FOR ANOTHER INSERT CARD

    IF (outcrd(2) == ectt(46) .AND. outcrd(3) == ectt(47)) GO TO 3605

    !     CHECK FOR ANOTHER DELETE CARD

    IF (outcrd(2) == ectt(49) .AND. outcrd(3) == ectt(50)) GO TO 3655

    !     CHECK FOR ENDALTER CARD

    IF (outcrd(2) /= endal(1) .OR.  outcrd(3) /= endal(2)) GO TO 380

    !     ENDALTER ENCOUNTERED

    IF (ectt(27) < 0) GO TO 720
    ectt(27) = orf (ectt(27), mask5)
    CALL eof (altfil)
    CALL CLOSE (altfil,2)
    altopn = -1
    notalt = 0
    GO TO 20



    !     WRITE DMAP INSTRUCTION ON THE ALTER SCRATCH FILE

380 IF (ectt(27) < 0) GO TO 30
    CALL WRITE (altfil, card, 18, 1)
    GO TO 30


    !     NOW PROCESS SOL CONTROL CARD

390 soluf = 1

    !     =====================================
    !     ECTT(I+2) = 0
    !     DO 2000 JJ = 1,12
    !2000 SOLU(JJ) = 0
    !     WRITE  (6,2001)
    !2001 FORMAT (16H0+++ OUTCARD +++)
    !     JJ = 1
    !2002 WRITE  (6,2003) JJ,OUTCRD(JJ)
    !2003 FORMAT (20X,I5,5X,O20)
    !     IF (OUTCRD(JJ) .EQ. ENDCD) GO TO 2004
    !     JJ = JJ + 1
    !     GO TO 2002
    !2004 CONTINUE
    !     =====================================

    IF (outcrd(1) == 1) GO TO 395

    DO  jj = 1,6
        solnmx(jj) = BLANK
    END DO
    jk = 2*outcrd(1) + 3
    solnmx(1) = outcrd(4)
    solnmx(2) = outcrd(5)
    IF (outcrd(1) == 2 .OR. outcrd(7) == BLANK) GO TO 392
    solnmx(3) = outcrd(6)
    solnmx(4) = outcrd(7)
    IF (outcrd(1) == 3 .OR. outcrd(9) == BLANK) GO TO 392
    solnmx(5) = outcrd(8)
    solnmx(6) = outcrd(9)
    392 loop394: DO  jj = 1,nsolnm
        DO  k  = 1,6
            IF (solnmx(k) /= solnms(k,jj)) CYCLE loop394
        END DO
        solu(1) = solnms(7,jj)
        GO TO 396
    END DO loop394
    iufile(1) = outcrd(4)
    iufile(2) = outcrd(5)
    solu(1)   = 0
    GO TO 396

395 imhere = 395
    IF (outcrd(4) /= -1) GO TO 750
    jk = 7
    solu(1) = outcrd(5)
    IF (outcrd(6) == 1) jk = jk + 3
    IF (outcrd(6) == 2) jk = jk + 5

396 CONTINUE
    rfflag = solu(1)
    IF (outcrd(jk-1) == endcd) GO TO 399
    imhere = 397
    jj = 1
397 jj = jj + 1
    IF (jj > 12) GO TO 750
    IF (outcrd(jk-1) /= -1) GO TO 750
    nsubst = jj
    solu(jj) = outcrd(jk)
    IF (outcrd(jk+1) == endcd) GO TO 399
    jk = jk + 2
    GO TO 397
399 CONTINUE

    !     ===========================================
    !2005 FORMAT (1H0,100(1H+)/1H0/1H0)
    !     WRITE  (6,2006)
    !2006 FORMAT (13H0+++ SOLU +++)
    !     JJ = 1
    !2007 IF (SOLU(JJ).EQ.0 .AND. JJ.GT.2) GO TO 2009
    !     WRITE  (6,2008) JJ,SOLU(JJ)
    !2008 FORMAT (20X,I5,5X,I10)
    !     JJ = JJ + 1
    !     GO TO 2007
    !2009 CONTINUE
    !     WRITE (6,2005)
    !     ===========================================

    GO TO 20


    !     B E G I N  CONTROL CARD
    !     PROCESS DMAP SEQUENCE

400 jj = 0
    WRITE  (outtap,410)
410 FORMAT (5X,'(SEE NASTRAN SOURCE PROGRAM COMPILATION FOR LISTING ',  &
        'OF DMAP SEQUENCE)')
    DO  jk = 1,nwpc
        jj = jj + 1
        dmapbf(jj) = card(jk)
    END DO
430 CALL xread (*1232,card)
    DO  jk = 1,nwpc
        jj = jj + 1
        dmapbf(jj) = card(jk)
    END DO
    IF (jj > dmapbs) GO TO 1290

    !     CHECK FOR END OR CEND CARD

    CALL xrcard (outcrd,200,card)

    !     CHECK FOR ERROR DETECTED BY XRCARD

    IF (xnogo == 0) GO TO 450
    WRITE (outtap,50) card
    IF (nogo  == 0) nogo = 1
    xnogo = 0
    GO TO 430
450 IF (outcrd(2) == cend(1) .AND. outcrd(3) == cend(2)) GO TO 940
    IF (outcrd(2) /= dmend) GO TO 430
    WRITE (outtap,50) card
    drecsz = jj
    GO TO 20


    !     NOW PROCESS UMF CARD
    !     CHECK FORMAT

460 WRITE  (outtap,465) uwm,ectt(i),ectt(i+1)
465 FORMAT (a25,', ',2A4,' CARD IS NO LONGER AVAILABLE')
    GO TO 20

! 460 IMHERE = 460
!     IF (OUTCRD(4).NE.INTGR .OR. OUTCRD(6).NE.INTGR .OR.
!    1    OUTCRD(5).LE.    0 .OR. OUTCRD(7).LT.   0) GO TO 750

!     SET UNSORTED AND SORTED BULK DATA OUTPUT (ECHO = BOTH)
!     AS THE DEFAULT FOR RUNS USING THE UMF

!     IECHO = 3

!     MAKE SURE UMF TAPE IS SETUP

!     ASSIGN 470 TO L
!     IDFIST = NUMF
!     GO TO 160
! 470 IF (NOSTUP .NE. 0) GO TO 970

!     MAKE SURE CORRECT UMF TAPE IS MOUNTED

!     NGINO = NUMF
!     IMHERE= 470
!     CALL OPEN  (*1320,NUMF,GBUFF(DMAPBS+1),0)
!     CALL READ  (*1350,*1350,NUMF,UMFID,1,0,FLGWRD)
!     CALL SKPFIL (NUMF,1)
!     CALL CLOSE (NUMF,2)
!     IF (UMFID .NE. OUTCRD(5)) GO TO 1000
!     UMFID = OUTCRD(7)
!     GO TO 20


!     PROCESS DIAG CARD
!     ALLOW MULTIPLE DIAG CARDS TO BE PROCESSED.

480 CONTINUE
    i = 2
490 i = i + 2
    IF (outcrd(i) ==     0) GO TO 505
    IF (outcrd(i) /= intgr) GO TO 520

    !     SET SENSE SWITCH BITS. (DIAG 1 THRU 48, BIT COUNTS 0 THRU 47)
    !     BITS 49 THRU 63 ARE RESERVED FOR LINK NO.  (-1 THRU -15)

    jj = outcrd(i+1)
    !WKBD IF (JJ .GT. 63-MAXLNK) GO TO 503
    !WKBD IF (JJ.GE.-MAXLNK .AND. JJ.LE.-1) JJ = 63 - MAXLNK - JJ
    IF (jj > 31) GO TO 500
    switch(1) = orf(lshift(1,jj-1),switch(1))

    !     TURN ON DIAG 14 IF DIAG 25 HAS BEEN REQUESTED

    IF (jj == 25) switch(1) = orf(lshift(1,13),switch(1))
    GO TO 503
500 IF (jj == 42 .AND. mach > 5) WRITE (outtap,501) uwm,mchnam
501 FORMAT (a25,', DIAG 42 IS UNSUPPORTED IN ALL UNIX MACHINES, ',  &
        'INCLUDING ',a6,' ***')
    jj = jj - 31
    switch(2) = orf(lshift(1,jj-1),switch(2))
503 CONTINUE
    GO TO 490

    !     DIAG CONTINUED ON NEXT CARD - READ IN NEXT CARD

505 ASSIGN 510 TO irtn1
    GO TO 30
510 IF (outcrd(2) == cend(1) .AND. outcrd(3) == cend(2)) GO TO 570
    i = -1
    GO TO 490

    !     SHOULD BE END OF LOGICAL DIAG CARD

520 imhere = 520
    IF (outcrd(i) /= endcd) GO TO 750
    !IBMDB 5/95
    !      SWITCH(3) = ORF(SWITCH(3),SWITCH(1))
    !      SWITCH(1) = 0
    !      CALL PRESSW (LINKS(1),I)

    !     RE-ACTIVATE THOSE LINK1 SPECIAL DIAGS IN DIAGX LIST IF NECESSARY

    !      IF (SWITCH(1) .EQ. SWITCH(3)) GO TO 527
    !      DO 525 I = 1,11
    !      JJ = DIAGX(I) - 1
    !      SWITCH(1) = ORF(ANDF(LSHIFT(1,JJ),SWITCH(3)),SWITCH(1))
    !  525 CONTINUE
    !      IF (SWITCH(1) .NE. SWITCH(3)) CALL PRESSW (RENTER,I)
    !IBMDE 5/95
527 CALL sswtch (15,l15)
    CALL sswtch (8 ,l_8)
    CALL sswtch (13,l13)
    GO TO 20


    !     NOW PROCESS ID CARD
    !     CHECK FORMAT - MUST BE AT LEAST 3 BCD FIELDS

530 imhere = 530
    IF (outcrd(1) < 3) GO TO 750

    !     MAKE SURE ID CARD IS FIRST CONTROL CARD
    !     IF ID CARD WAS IN ERROR CONTROL WILL STILL RETURN TO HERE

    531 DO  i = 1,lectt,3
        IF (ectt(i+2) < 0 .AND. ectt(i) /= id) GO TO 1060
    END DO
    IF (logfl <= 0) CALL logfil (card)
    DO  jj = 1,4
        tapid(jj) = outcrd(jj+3)
    END DO

    !      PACK DATE -

    imnth = lshift(idate(1),14)
    iday  = lshift(idate(2),8)
    iyear = idate(3)
    tapid(5) = orf(imnth,orf(iday,iyear))

    !     REEL NO. TO TAPID

    tapid(6) = 1

    !     OUTPUT IF ON NEW PROBLEM TAPE

    ngino = ptape
    CALL OPEN  (*1320,ptape,gbuff(dmapbs+1),1)
    CALL WRITE (ptape,tapid,6,0)
    CALL WRITE (ptape,timew,1,1)
    CALL eof   (ptape)
    CALL CLOSE (ptape,2)
    GO TO 20


    !     PROCESS INTERACTIVE CARD
    !     SET INTRA TO NEGATIVE IN BATCH RUN (I.E. PRE-INTERACTIVE RUN)
    !     INTRA WILL BE RESET TO POSITIVE IN AN ON-LINE INTERACTIVE RUN

    !     CHECK FORMAT AND FILE ASSIGNMENT

555 intra = 0
    DO  jj = 4,9
        IF (outcrd(jj) == plot) intra = orf(intra,1)
        IF (outcrd(jj) == prnt) intra = orf(intra,2)
        IF (outcrd(jj) == both) intra = orf(intra,3)
    END DO
    IF (intra == 0) GO TO 700
    intra = -intra
    jj = 1
    IF (mach == 3) CALL facil (inp9,jj)
    IF (jj   == 2) GO TO 1250
    GO TO 20


    !     UMFEDIT CARD FOUND - SET EDTUMF FLAG

560 WRITE (outtap,465) uwm,ectt(i),ectt(i+1)
    !     EDTUMF = 1
    GO TO 20


    !     PROCESS PREC CARD

565 imhere = 565
    IF (outcrd(5) /= 1 .AND. outcrd(5) /= 2) GO TO 750
    prec = outcrd(5)
    GO TO 20

    !     CEND CARD FOUND - NO MORE CONTROL CARDS TO PROCESS


    !     SET APP DEFAULT TO 'DISPLACEMENT' AND TIME TO 10 MINUTES

570 IF (apprch /= 0) GO TO 572
    apprch  = 2
    apprec  = 2
    WRITE  (outtap,571)
571 FORMAT ('0*** APP  DECLARATION CARD MISSING.  DISPLACEMENT IS ',  &
        'SELECTED BY DEFAULT')
572 IF (time > 0) GO TO 575
    time = 300
    WRITE  (outtap,573)
573 FORMAT ('0*** TIME  CARD MISSING. MAXIMUM EXECUTION TIME IS SET ',  &
        'TO 5 MINUTES BY DEFAULT')

    !     CALL NSINFO TO PRINT DIAG48, OR
    !     PRINT THE FOLLOWING MESSAGE OUT ONLY IF THE JOB IS RUN ON THE SAME
    !     YEAR OF THE RELEASE DATE, AND USER DOES NOT MAKE A DIAG48 REQUEST

    !     DIAG48 TEXT IS STORED IN 4TH SECTION OF THE NASINFO FILE


575 CALL sswtch (48,jj)
    IF (jj /= 1) GO TO 576
    CALL nsinfo (4)
    GO TO 580
576 jj = idate(3)
    jj = MOD(jj,100)
    CALL int2a8 (*577,jj,iz(1))
577 IF (iz(1) == sy42(3)) WRITE (outtap,578) uim
578 FORMAT (//,a29,', TURN DIAG 48 ON FOR NASTRAN RELEASE NEWS, ',  &
        'DIAG DEFINITION, NEW DMAP', /9X,  &
        'MODULES AND NEW BULKDATA CARDS INFORMATION')

    !     CLOSE NASINFO FILE IF IT EXISTS
    !     AND RESET THE 37TH WORD OF /SYSTEM/ BACK TO ZERO

580 IF (lu /= 0) CLOSE (UNIT=lu)
    lu = 0

    !     NOW MAKE SURE ALL NECESSARY CARDS HAVE BEEN FOUND

    DO  i = 1,lectt,3
        test = andf(ectt(i+2),mask)
        IF (test > 0) IF (ectt(i+2)) 590,1090,1090
590 CONTINUE
    END DO

    !     SET APPRCH NEGATIVE FOR RESTART

    IF (rstrt /= icold) apprch = -apprch
    IF (soluf == 1 .AND. drecsz /= 0) GO TO 1120
    IF (soluf == 0 .AND. drecsz == 0) GO TO 1150
    !     IF (RSTRT.NE.ICOLD .AND. UMFID.NE.0) GO TO 1030


600 IF (nogo > 1) GO TO 1380

    !     WRITE XCSA CONTROL FILE ONTO PROBLEM TAPE
    !     FIRST RECORD IS HEADER RECORD CONTAINING A SINGLE WORD (XCSA)

    IF (apprec == appdmp) GO TO 610

    !     IF APPROACH IS HEAT ADD TWENTY THREE TO SOLUTION

    IF (apprec == apphea) solu(1) = solu(1) + 23

    !     IF APPROACH IS AEROELASTIC ADD THIRTY TO SOLUTION

    IF (apprec == appaer) solu(1) = solu(1) + 30
    GO TO 612
610 ngino = ptape
    imhere= 610
    CALL OPEN  (*1320,ptape,gbuff(dmapbs+1),3)
    CALL WRITE (ptape,nxcsa,2,1)

    !     DIS OLD PT HAVE AN ALTER FILE AND/OR CKPT DIST

    solrec(4) = oldalt

    !     WRITE SIX-WORD CONTROL FILE RECORD

    CALL WRITE (ptape,solrec,6,1)
    CALL eof   (ptape)
    CALL CLOSE (ptape, 3)
    IF (apprec /= appdmp) GO TO 640
612 ngino = nscr
    imhere= 612
    CALL OPEN (*1320,nscr,gbuff(dmapbs+1),1)
    IF (apprec == appdmp) GO TO 620

    !     APPROACH IS RIGID FORMAT
    !     WRITE RIGID FORMAT AND MED TABLES ONTO SCRATCH FILE

    isize = korsz (dmapbf(1)) - ibufsz
    IF (altopn == 0) GO TO 614
    IF (erralt == 0) GO TO 613
    newalt = 0
613 IF (newalt == 0) GO TO 614
    isize = isize - ibufsz
    ngino = altfil
    CALL OPEN (*1320, altfil, gbuff(altrbs+1), 3)
614 CALL xrgdfm (solu,osolu,apprec,iufile,dmapbf,isize,nscr,nogo)
    IF (xnogo == 0) GO TO 615
    IF (nogo  == 0) nogo = 1
    xnogo = 0
615 CONTINUE
    IF (nogo > 1) GO TO 1380
    CALL CLOSE (nscr, 1)
    solrec(3) = 0
    IF (altopn == 0) GO TO 610
    IF (erralt == 1) GO TO 610
    solrec(3) = 1
    ngino = ptape
    CALL OPEN (*1320, ptape,  gbuff(dmapbs+1), 3)
    ngino = altfil
    CALL OPEN (*1320, altfil, gbuff(altrbs+1), 0)
    CALL dmpalt (isize, dmapbf, ptape)
    CALL eof (ptape)
    CALL CLOSE (ptape,  2)
    CALL CLOSE (altfil, 1)
    IF (alnogo == 0) GO TO 610
    IF (nogo < 2) nogo = 2
    GO TO 610

    !     APPROACH IS DMAP
    !     WRITE DMAP SEQUENCE ONTO SCRATCH FILE FROM OPEN CORE

620 CALL WRITE (nscr,dmapbf,drecsz,1)
630 CALL CLOSE (nscr,1)
640 CONTINUE

    !     PUNCH RESTART CARD IF CHECKPOINT FLAG IS SET.

    IF (icpflg == 0) GO TO 660
    !      IF (IROPEN .EQ. 1) GO TO 6405
    !      OPEN (UNIT=4, FILE=DSNAMES(4), STATUS='UNKNOWN')
    !      IROPEN = 1
    WRITE (irdict,641) (tapid(i),i=1,4),(idate(j),j=1,3),timew
641 FORMAT (9HRESTART  ,2A4,1H,,2A4,1H,,i2,1H/,i2,1H/,i2,1H,,i8,1H,)
    CALL sswtch (9,diag09)
    IF (diag09 == 1) GO TO 660
    CALL page
    WRITE  (outtap,651) (tapid(i),i=1,4),(idate(j),j=1,3),timew
651 FORMAT ('0ECHO OF FIRST CARD IN CHECKPOINT DICTIONARY TO BE ',  &
        'PUNCHED OUT FOR THIS PROBLEM', /  &
        14H0   restart   ,2A4,1H,,2A4,1H,,i2,1H/,i2,1H/,i2,1H,,i8,1H,)
660 xnogo = nogo
    RETURN

    !     ERROR MESSAGES

    !     USER  FATAL MESSAGES

670 nlines = nlines + 2
    IF (nlines >= nlpp) CALL page
    IF (nogo   <    1) nogo = 1
    ignore = 1
    GO TO irtn2, ( 700, 730, 770, 800, 830, 860,      920, 950,  &
        1070,1100,1130,1160,1190,1220,1234)

690 ASSIGN 700 TO irtn2
    msgnum = 505
    GO TO 670
700 WRITE  (outtap,710) ufm,msgnum,outcrd(2),outcrd(3)
710 FORMAT (a23,i5,', CONTROL CARD ',2A4,11H is illegal)
    GO TO 20

720 ASSIGN 730 TO irtn2
    msgnum = 506
    GO TO 670
730 WRITE  (outtap,740) ufm,msgnum,outcrd(2),outcrd(3)
740 FORMAT (a23,i5,', CONTROL CARD ',2A4,11H duplicated)
    GO TO 20

750 CONTINUE
    erralt = 1
760 ASSIGN 770 TO irtn2
    msgnum = 507
    GO TO 670
770 WRITE  (outtap,780) ufm,msgnum,imhere
780 FORMAT (a23,i5,', ILLEGAL SPECIFICATION OR FORMAT ON PRECEDING ',  &
        'CARD.', /5X,'IMHERE =',i5)
    IF (outcrd(2) == ectt(34) .AND. outcrd(3) == ectt(35)) GO TO 531
    GO TO 20

790 ASSIGN 800 TO irtn2
    msgnum = 508
    GO TO 670
800 WRITE  (outtap,810) ufm,msgnum
810 FORMAT (a23,i5,', PROBLEM TAPE MUST BE ON PHYSICAL TAPE FOR ',  &
        'CHECK POINTING')
    ignore = 0
    icpflg = 0
    GO TO 20

820 ASSIGN 830 TO irtn2
    msgnum = 509
    GO TO 670
830 WRITE  (outtap,840) ufm,msgnum,(otapid(i),i=1,4),imnth,iday,  &
        iyear,timex,otapid(6)
840 FORMAT (a23,i5,', WRONG OLD TAPE MOUNTED.', /30X,  &
        23H OLD problem tape id = ,2A4,1H,,2A4,1H,,i2,1H/,i2,1H/,  &
        i2,1H,,2X,i8,1H,,5X,10HREEL no. =,i4)
    GO TO 1410

850 ASSIGN 860 TO irtn2
    msgnum = 512
    GO TO 670
860 WRITE  (outtap,870) ufm,msgnum
870 FORMAT (a23,i5,', OLD PROBLEM TAPE IS MISSING AND IS NEEDED FOR ',  &
        'RESTART')
    nogo = 3
    GO TO 20


910 ASSIGN 920 TO irtn2
    msgnum = 514
    GO TO 670
920 WRITE  (outtap,930) ufm,msgnum
930 FORMAT (a23,i5,', ENDALTER CARD IS MISSING')
    IF (nogo < 2) nogo = 2
    GO TO 570

940 ASSIGN 950 TO irtn2
    msgnum = 515
    GO TO 670
950 WRITE  (outtap,960) ufm,msgnum
960 FORMAT (a23,i5,', END INSTRUCTION MISSING IN DMAP SEQUENCE')
    IF (nogo < 2) nogo = 2
    GO TO 570

1060 ASSIGN 1070 TO irtn2
    msgnum = 519
    GO TO 670
1070 WRITE  (outtap,1080) ufm,msgnum
1080 FORMAT (a23,i5,', ID CARD MUST PRECEDE ALL OTHER CONTROL CARDS')
    nogo = 3
    GO TO 20

1090 ASSIGN 1100 TO irtn2
    msgnum = 520
    GO TO 670
1100 WRITE  (outtap,1110) ufm,msgnum,ectt(i),ectt(i+1)
1110 FORMAT (a23,i5,', CONTROL CARD ',2A4,' IS MISSING')
    ectt(i+2) = orf(ectt(i+2),mask5)
    IF (ectt(i) /= ectt(4)) GO TO 570

    !     MISSING CARD IS APP

    IF (nogo < 2) nogo = 2
    GO TO 570

1120 ASSIGN 1130 TO irtn2
    msgnum = 521
    GO TO 670
1130 WRITE  (outtap,1140) ufm,msgnum
1140 FORMAT (a23,i5,', SPECIFY A SOLUTION OR A DMAP SEQUENCE BUT NOT ',  &
        'BOTH')
    IF (nogo < 2) nogo = 2
    GO TO 1380

1150 ASSIGN 1160 TO irtn2
    msgnum = 522
    GO TO 670
1160 WRITE  (outtap,1170) ufm,msgnum
1170 FORMAT (a23,i5,', NEITHER A SOL CARD NOR A DMAP SEQUENCE WAS ',  &
        'INCLUDED')
    IF (nogo < 2) nogo = 2
    GO TO 1380

1180 ASSIGN 1190 TO irtn2
    notalt = 0
    msgnum = 523
    GO TO 670
1190 WRITE  (outtap,1200) ufm,msgnum
1200 FORMAT (a23,i5,', ENDALTER CARD OUT OF ORDER')
    GO TO 20

1210 ASSIGN 1220 TO irtn2
    msgnum = 526
    GO TO 670
1220 WRITE  (outtap,1230) ufm,msgnum
1230 FORMAT (a23,i5,', CHECKPOINT DICTIONARY OUT OF SEQUENCE - ',  &
        'REMAINING RESTART CARDS IGNORED')
    GO TO 20
1232 ASSIGN 1234 TO irtn2
    msgnum = 529
    GO TO 670
1234 WRITE  (outtap,1236) ufm,msgnum
1236 FORMAT (a23,i5,', MISSING CEND CARD.')
    nogo = 3
    GO TO 1380

    !     SYSTEM FATAL MESSAGES

1240 nlines = nlines +2
    IF (nlines >= nlpp) CALL page
    IF (nogo   <    2) nogo = 2
    ignore = 1
    GO TO irtn3, (1255,1270,1300,1330,1360)

1250 ASSIGN 1255 TO irtn3
    msgnum = 530
    GO TO 1240
1255 WRITE  (outtap,1256) sfm,msgnum
1256 FORMAT (a25,i5,2H, , /5X,'INP9 FILE WAS NOT ASSIGNED FOR ',  &
        'NASTRAN INTERACTIVE POST-PROCESSOR',/)
    GO TO 20
1260 ASSIGN 1270 TO irtn3
    msgnum = 510
    GO TO 1240
1270 WRITE  (outtap,1280) sfm,msgnum
1280 FORMAT (a25,i5,', CHECKPOINT DICTIONARY EXCEEDS CORE SIZE - ',  &
        'REMAINING RESTART CARDS IGNORED')
    GO TO 20

1290 ASSIGN 1300 TO irtn3
    msgnum = 511
    GO TO 1240
1300 WRITE  (outtap,1310) sfm,msgnum
1310 FORMAT (a25,i5,', DMAP SEQUENCE EXCEEDS CORE SIZE - ',  &
        'REMAINING DMAP INSTRUCTIONS IGNORED')
    IF (nogo < 2) nogo = 2
    GO TO 20

1320 ASSIGN 1330 TO irtn3
    msgnum = 524
    GO TO 1240
1330 WRITE  (outtap,1340) sfm,msgnum,ngino,imhere
1340 FORMAT (a25,i5,', ALTERNATE RETURN TAKEN WHEN OPENING FILE ',a4,  &
        3X,1H-,i3)
    nogo = 3
    GO TO 1410

1350 ASSIGN 1360 TO irtn3
    msgnum = 525
    GO TO 1240
1360 WRITE  (outtap,1370) sfm,msgnum,ngino
1370 FORMAT (a25,i5,', ILLEGAL FORMAT ENCOUNTERED WHILE READING FILE ', a4)
    nogo = 3
    GO TO 1410

1380 SELECT CASE ( nogo )
        CASE (    1)
            GO TO 600
        CASE (    2)
            GO TO 1400
        CASE (    3)
            GO TO 1390
    END SELECT

    !     NOGO = 3 - TERMINATE JOB HERE
1390 icpflg = 0
    CALL mesage (-61,0,0)

    !     NOGO = 2 - PUT IN DUMMY CONTROL FILE ON PROBLEM TAPE
1400 ngino = ptape
    CALL CLOSE (ptape,1)
    CALL OPEN  (*1320,ptape,gbuff(dmapbs+1),0)
    CALL skpfil(ptape,1)
    CALL CLOSE (ptape,2)
    CALL OPEN  (*1320,ptape,gbuff(dmapbs+1),3)
    CALL WRITE (ptape,nxcsa,2,1)
    solu(1) = 0
    solu(2) = 0
    apprch  = appdmp
    IF (rstrt /= icold) apprch = -apprch
    CALL WRITE (ptape,solrec,6,1)
    CALL eof   (ptape)
    CALL CLOSE (ptape,3)
    GO TO 640

    !     XCSA HAS BEEN DISASTERED - GET DUMP AND QUIT.
1410 icpflg = 0
1420 CALL mesage (-37,0,nxcsa)

    RETURN
END SUBROUTINE xcsa
