SUBROUTINE exio1
     
!     EXIO1 SERVICES INTERNAL FORMAT FUNCTIONS FOR EXIO.
 
 EXTERNAL lshift   ,rshift   ,andf     ,orf
 LOGICAL :: first    ,opnsof   ,ditup    ,mdiup    ,nxtup     ,  &
     nxtrst   ,tapbit
 INTEGER :: dry      ,cor(1)   ,device   ,uname    ,pos       ,  &
     datype   ,pdate    ,ptime    ,time     ,sec       ,  &
     hours    ,ssname   ,savrec   ,hdrec    ,rewi2     ,  &
     buf      ,sysbuf   ,date     ,rd       ,rdrew     ,  &
     wrt      ,wrtrew   ,rew      ,eofnrw   ,filnam    ,  &
     filsiz   ,STATUS   ,passwd   ,blksiz   ,dirsiz    ,  &
     supsiz   ,avblks   ,dit      ,ditpbn   ,ditlbn    ,  &
     ditsiz   ,ditnsb   ,ditbl    ,z        ,tape      ,  &
disk     ,sofin    ,sofout   ,check    ,APPEND    ,  &
    comprs   ,rewi     ,eqf      ,all      ,tables    ,  &
    phase3   ,dump     ,restor   ,whole(2) ,xitems(50),  &
    subr(2)  ,BLANK    ,sof      ,eoi      ,hdr       ,  &
    q        ,qqqq     ,xxxx     ,scr1     ,srd       ,  &
    swrt     ,rshift   ,andf     ,sofsiz   ,orf       ,  &
    buf1     ,buf2     ,buf3     ,buf4     ,UNIT      ,  &
    rc       ,flag     ,oldtsz   ,buf5     ,scr2      ,  &
    head1    ,head2    ,inblk(15),outblk(15)
CHARACTER (LEN=27) :: swm
CHARACTER (LEN=25) :: sfm
CHARACTER (LEN=29) :: uim
CHARACTER (LEN=31) :: sim
CHARACTER (LEN=25) :: uwm
CHARACTER (LEN=23) :: ufm
COMMON  /xmssg /   ufm      ,uwm      ,uim      ,sfm       , swm      ,sim
COMMON  /machin/   mach     ,ihalf    ,jhalf
COMMON  /BLANK /   dry      ,xmach    ,device(2),uname(2)  ,  &
    formt(2) ,mode(2)  ,pos(2)   ,datype(2) ,  &
    names(10),pdate    ,ptime    ,time(3)   ,  &
    ssname(2),savrec(9),hdrec(10),buf(10)
COMMON  /system/   sysbuf   ,nout     ,x1(6)    ,nlpp      ,  &
    x2(2)    ,line     ,x3(2)    ,date(3)   , x4(21)   ,nbpc     ,nbpw     ,ncpw
COMMON  /names /   rd       ,rdrew    ,wrt      ,wrtrew    ,  &
    rew      ,norew    ,eofnrw
COMMON  /sofcom/   nfiles   ,filnam(10)         ,filsiz(10),  &
    STATUS   ,passwd(2),first    ,opnsof
COMMON  /sys   /   blksiz   ,dirsiz   ,supsiz   ,avblks    , noblks   ,ifrst
COMMON  /itemdt/   nitem    ,items(7,1)
COMMON  /output/   head1(96),head2(96)
COMMON  /sof   /   dit      ,ditpbn   ,ditlbn   ,ditsiz    ,  &
    ditnsb   ,ditbl    ,io       ,iopbn     ,  &
    iolbn    ,iomode   ,ioptr    ,iosind    ,  &
    ioitcd   ,ioblk    ,mdi      ,mdipbn    ,  &
    mdilbn   ,mdibl    ,nxt      ,nxtpbn    ,  &
    nxtlbn   ,nxttsz   ,nxtfsz(10)          ,  &
    nxtcur   ,ditup    ,mdiup    ,nxtup     , nxtrst
COMMON  /zzzzzz/   z(1)
EQUIVALENCE        (cor(1)  ,z(1))
EQUIVALENCE        (time(1) ,hours)   ,(time(2) ,MIN)      , (time(3) ,sec)
DATA     tape     ,disk     ,sofin    ,sofout   ,check     /  &
    4HTAPE   ,4HDISK   ,4HSOFI   ,4HSOFO   ,4HCHEC    /,  &
APPEND   ,comprs   ,norewi   ,rewi     ,eqf       /  &
    4HAPPE   ,4HCOMP   ,4HNORE   ,4HREWI   ,4HEOF     /,  &
    all      ,matric   ,tables   ,phase3   ,dump      /  &
    4HALL    ,4HMATR   ,4HTABL   ,4HPHAS   ,4HDUMP    /,  &
    restor   ,whole              ,rewi2               /  &
    4HREST   ,4HWHOL   ,4HESOF   ,4HND                /,  &
    subr               ,BLANK    ,sof      ,eoi       /  &
    4HEXIO   ,4H1      ,4H       ,4HSOF    ,4HEOI     /,  &
    id       ,hdr      ,q        ,qqqq     ,xxxx      /  &
    4H$id$   ,4H$hd$   ,4HQ      ,4HQQQQ   ,4HXXXX    /
DATA     scr1     ,scr2     ,srd      ,swrt     ,iz2       /  &
    301      ,302      ,1        ,2        ,2         /

!     INITIALIZE

IF (nitem > 50) CALL errmkn (25,10)
lcore = korsz(z)
buf1  = lcore- sysbuf + 1
buf2  = buf1 - sysbuf - 1
buf3  = buf2 - sysbuf
buf4  = buf3 - sysbuf
buf5  = buf4 - sysbuf
lcore = buf5 - 1
ncore = lcore
nos   = 0
idm   = 1
IF (lcore <= 0) CALL mesage (-8,0,subr)
IF (mode(1) /= restor) CALL sofopn (z(buf1),z(buf2),z(buf3))
UNIT  = uname(1)

!     CHECK TAPE BIT IF DEVICE=TAPE

IF (device(1) == disk .OR. mode(1) == comprs .OR. mode(1) == APPEND) GO TO 10
IF (device(1) /= tape) GO TO 1810
IF (.NOT.  tapbit(UNIT)) GO TO 1800

!     SET REWIND VARIABLE

!     IF SOFOUT COMMAND POSITION TO END-OF-FILE IF REQUESTED

!     IF POSITION = REWIND AND WE ARE WRITING THEN BCKREC OVER LAST EOF

!     IF POSITION = EOF AND WE ARE WRITING THEN BCKREC FIRST TO INSURE
!     WE ARE INFRONT OF AND EOF AND THEN SEARCH FOR EOF

10 ipos = -1
IF (pos(1) == norewi .OR. pos(1) == eqf) ipos = 2
IF (pos(1) == rewi) ipos = 0
IF (ipos < 0) GO TO 1830
IF (mode(1) == dump .OR. mode(1) == restor) ipos = 0
IF (ipos /= 0) GO TO 20
head2(13) = rewi
head2(14) = rewi2
20 IF (mode(1) /= sofout) GO TO 60
IF (ipos == 0) GO TO 60
CALL OPEN (*1860,UNIT,z(buf4),rd)
CALL bckrec (UNIT)
IF (pos(1) == norewi) GO TO 50
30 CALL fwdrec (*40,UNIT)
GO TO 30
40 CALL bckrec (UNIT)
50 CALL CLOSE (UNIT,norew)

!     BRANCH ON MODE OF OPERATION

60 IF (mode(1) == sofout .OR. mode(1) ==  dump) GO TO 70
IF (mode(1) == sofin  .OR. mode(1) == restor) GO TO 370
IF (mode(1) == check ) GO TO 1160
IF (mode(1) == APPEND) GO TO 1220
IF (mode(1) == comprs) GO TO 1500
GO TO 1820


!     **********************   W R I T E   **********************

!     OPEN FILE AND WRITE 9 WORD ID RECORD

70 CALL OPEN (*1860,UNIT,z(buf4),wrtrew+ipos)
CALL waltim (sec)
hours= sec/3600
sec  = MOD(sec,3600)
MIN  = sec/60
sec  = MOD(sec,60)
hdrec(1) = id
hdrec(2) = passwd(1)
hdrec(3) = passwd(2)
DO  i = 1,3
  hdrec(i+3) = date(i)
  hdrec(i+6) = time(i)
END DO
CALL WRITE (UNIT,hdrec,9,1)
CALL page
WRITE (nout,2120) uim,passwd,date,time
line = line + 1

!     WRITE DIT AND MDI CONTROL WORDS

n = ditsiz/2
CALL WRITE (UNIT,n,1,0)
DO  i = 1,n
  CALL fdit(i,j)
  CALL WRITE (UNIT,cor(j),2,0)
  CALL fmdi  (i,j)
  CALL WRITE (UNIT,cor(j+1),2,0)
END DO
CALL WRITE (UNIT,0,0,1)
CALL WRITE (UNIT,eoi,1,1)
IF (mode(1) /= dump) GO TO 110


!     DUMP FORM --

!     COPY OUT ALL SOF SUPERBLOCKS WHICH HAVE BEEN USED WITHOUT REGARD
!     TO THE DATA SEQUENCE OR CONTENT.


DO  i = 1,noblks
  CALL sofio (srd,i,z(buf1))
  CALL WRITE (UNIT,z(buf1+3),blksiz,0)
END DO
CALL WRITE (UNIT,0,0,1)
CALL CLOSE (UNIT,rew)
WRITE (nout,2130) uim,noblks,nxttsz,uname
GO TO 1740

!     STANDARD FORM --

!     COPY OUT EACH SUBSTRUCTURE/ITEM WITH ITS DATA IN THE CORRECT
!     SEQUENCE.

!     SETUP THE ARRAY XITEMS OF NAMES OF ITEMS TO BE COPIED.

110 IF (datype(1) /= all) GO TO 130
nitems = nitem
DO  i = 1,nitem
  xitems(i) = items(1,i)
END DO
GO TO 200
130 IF (datype(1) /= tables) GO TO 150
nitems = 0
DO  i = 1,nitem
  IF (items(2,i) > 0) CYCLE
  nitems = nitems + 1
  xitems(nitems) = items(1,i)
END DO
GO TO 200
150 IF (datype(1) /= matric) GO TO 170
nitems = 0
DO  i = 1,nitem
  IF (items(2,i) <= 0) CYCLE
  nitems = nitems + 1
  xitems(nitems) = items(1,i)
END DO
GO TO 200
170 IF (datype(1) /= phase3) GO TO 190
nitems = 0
DO  i = 1,nitem
  IF (andf(items(7,i),8) == 0) CYCLE
  nitems = nitems + 1
  xitems(nitems) = items(1,i)
END DO
GO TO 200
190 nitems = 2
xitems(1) = datype(1)
xitems(2) = datype(2)
IF (xitems(2) == BLANK) nitems = 1

!     LOOP OVER SUBSTRUCTURE NAMES.  FOR EACH SUBSTRUCTURE, WRITE OUT
!     THE NITEMS IN XITEMS.

200 iss = 0
210 iss = iss + 1
IF (names(1) /= whole(1) .OR. names(2) /= whole(2)) GO TO 220

!     WRITE ALL SUBSTRUCTURES IN THE RESIDENT SOF.

IF (iss > ditsiz/2) GO TO 360
CALL fdit (iss,i)
IF (cor(i) == BLANK) GO TO 210
ssname(1) = cor(i  )
ssname(2) = cor(i+1)
GO TO 230

!     WRITE ONLY THOSE SUBSTRUCTURES IN THE PARAMETER LIST

220 IF (iss > 5) GO TO 360
IF (names(2*iss-1) == xxxx) GO TO 210
ssname(1) = names(2*iss-1)
ssname(2) = names(2*iss  )

!     LOOP OVER ALL ITEMS OF THIS SUBSTRUCTURE.

230 DO  item = 1,nitems
  kdh = ittype(xitems(item))
  IF (kdh == 1) GO TO 260
  CALL sfetch (ssname,xitems(item),srd,rc)
  SELECT CASE ( rc )
    CASE (    1)
      GO TO 260
    CASE (    2)
      GO TO 240
    CASE (    3)
      GO TO 350
    CASE (    4)
      GO TO 250
    CASE (    5)
      GO TO 250
  END SELECT
  240 line = line + 2
  IF (line > nlpp) CALL page
  WRITE (nout,2160) uwm,ssname,xitems(item)
  CYCLE
  250 line = line+2
  IF (line > nlpp) CALL page
  CALL smsg (rc-2,xitems(item),ssname)
  CYCLE
  
!     WRITE SUBSTRUCTURE/ITEM HEADER RECORD
  
  260 CALL waltim (sec)
  hours = sec/3600
  sec   = MOD(sec,3600)
  MIN   = sec/60
  sec   = MOD(sec,60)
  hdrec(1) = hdr
  hdrec(2) = ssname(1)
  hdrec(3) = ssname(2)
  hdrec(4) = xitems(item)
  DO  i = 1,3
    hdrec(i+4) = date(i)
    hdrec(i+7) = time(i)
  END DO
  IF (kdh == 1) GO TO 310
  CALL WRITE (UNIT,hdrec,10,1)
  
!     COPY DATA
  
  280 CALL suread (z(1),lcore,nwds,rc)
  SELECT CASE ( rc )
    CASE (    1)
      GO TO 290
    CASE (    2)
      GO TO 300
    CASE (    3)
      GO TO 340
  END SELECT
  290 CALL WRITE (UNIT,z,lcore,0)
  GO TO 280
  300 CALL WRITE (UNIT,z,nwds,1)
  GO TO 280
  
!     COPY MATRIX DATA ITEMS
  
  310 ifile = scr1
  CALL mtrxi (scr1,ssname,xitems(item),0,rc)
  SELECT CASE ( rc )
    CASE (    1)
      GO TO 320
    CASE (    2)
      GO TO 240
    CASE (    3)
      GO TO 350
    CASE (    4)
      GO TO 250
    CASE (    5)
      GO TO 250
    CASE (    6)
      GO TO 2010
  END SELECT
  320 CALL WRITE (UNIT,hdrec,10,1)
  z(1) = scr1
  CALL rdtrl (z(1))
  CALL WRITE (UNIT,z(iz2),6,1)
  CALL OPEN  (*2010,scr1,z(buf5),rdrew)
  CALL cpyfil(scr1,UNIT,z,lcore,icount)
  CALL CLOSE (scr1,1)
  
!     WRITE END-OF-ITEM RECORD AND USER MESSAGE
  
  340 CALL WRITE (UNIT,eoi,1,1)
  line = line + 1
  IF (line > nlpp) CALL page
  WRITE (nout,2170) uim,ssname,xitems(item),sof,UNIT,date,time
  350 CONTINUE
END DO

!     BOTTOM OF LOOP OVER SUBSTRUCTURES

GO TO 210

!     ALL SUBSTRUCTURE/ITEMS HAVE NOW BEEN COPIED.  CLOSE WITH EOF AND
!     NO REWIND (IN CASE MORE DATA TO FOLLOW).

!     WRITE EOF FOR NOU BECAUSE LEVEL 16 GINO OPT=3 DOESN T AS
!     ADVERTISED

360 CALL eof (UNIT)
CALL CLOSE (UNIT,eofnrw)
GO TO 1740

!     ***********************   R E A D  ************************

!     BRANCH FOR RESTORE OR STANDARD READ

370 IF (mode(1) /= restor) GO TO 400

!     RESTORE FORM --

!     COPY EACH LOGICAL RECORD ON THE EXTERNAL FILE INTO CONSEQUTIVE,
!     CONTIGUOUS BLOCKS ON THE RESIDENT SOF.

!     MAKE SURE THE RESIDENT SOF IS EMPTY.

IF (STATUS /= 0) GO TO 1840
CALL sofopn (z(buf1),z(buf2),z(buf3))
CALL sofcls

!     OPEN FILE AND READ THE ID RECORD

CALL OPEN (*1860,UNIT,z(buf4),rdrew)
CALL READ (*1850,*1850,UNIT,hdrec,9,1,flag)
IF (hdrec(1) /= id) GO TO 1850
CALL page
line = line+1
WRITE (nout,2120) uim,(hdrec(i),i=2,9)
CALL fwdrec (*1870,UNIT)
CALL fwdrec (*1870,UNIT)

!     BEGIN DATA TRANSFER

i = 1
380 CALL READ (*1870,*390,UNIT,z(buf1+3),blksiz,0,flag)
CALL sofio (swrt,i,z(buf1))
i = i+1
GO TO 380

!     RESTORE COMPLETE.  CLOSE FILE AND GIVE USER THE NEWS.

390 CALL CLOSE (UNIT,rew)
i = i - 1
WRITE (nout,2200) uim,i
GO TO 1750

!     STANDARD FORM -

!     COPY IN EACH INDIVIDUAL SUBSTRUCTURE/ITEM.

400 iss = 0

!     SETUP ARRAY OF NAMES OF ITEMS TO BE COPIED.

IF (datype(1) /= all) GO TO 420
nitems = nitem
DO  i = 1,nitem
  xitems(i) = items(1,i)
END DO
GO TO 490
420 IF (datype(1) /= tables) GO TO 440
nitems = 0
DO  i = 1,nitem
  IF (items(2,i) > 0) CYCLE
  nitems = nitems + 1
  xitems(nitems) = items(1,i)
END DO
GO TO 490
440 IF (datype(1) /= matric) GO TO 460
nitems = 0
DO  i = 1,nitem
  IF (items(2,i) <= 0) CYCLE
  nitems = nitems + 1
  xitems(nitems) = items(1,i)
END DO
GO TO 490
460 IF (datype(1) /= phase3) GO TO 480
nitems = 0
DO  i = 1,nitem
  IF (andf(items(7,i),8) == 0) CYCLE
  nitems = nitems + 1
  xitems(nitems) = items(1,i)
END DO
GO TO 490
480 nitems = 2
xitems(1) = datype(1)
xitems(2) = datype(2)
IF (xitems(2) == BLANK) nitems = 1

!     DETERMINE NUMBER OF SUBSTRUCTURE/ITEMS TO BE COPIED AND INITIALIZE
!     COUNTER.

490 jcopy = 0
ncopy = 0
IF (names(1) == whole(1) .AND. names(2) == whole(2)) GO TO 510
DO  i = 1,5
  IF (names(2*i-1) /= xxxx) ncopy = ncopy + 1
END DO
ncopy = ncopy*nitems
IF (pdate /= 0) ncopy = 1

!     OPEN THE EXTERNAL FILE AND READ THE IDENTIFICATION OR HEADER
!     RECORD.
!     REMEMBER IT IN CASE THE USER HAS REQUESTED A SUBSTRUCTURE/ITEM
!     WHICH IS NOT PRESENT ON THE FILE.

510 CALL page
CALL OPEN (*1860,UNIT,z(buf4),rdrew+ipos)
520 CALL READ (*530,*540,UNIT,hdrec,10,1,lrec1)
lrec1 = 10
GO TO 540
530 CALL REWIND (UNIT)
GO TO 520
540 DO  i = 1,lrec1
  buf(i) = hdrec(i)
END DO
IF (hdrec(1) /=  id) GO TO 560
GO TO 610
560 IF (hdrec(1) /= hdr) GO TO 1850
GO TO 610

!     SCAN THROUGH THE EXTERNAL TAPE.  FOR EACH SUBSTRUCTURE/ITEM
!     ENCOUNTERED, CHECK TO SEE IF IT SHOULD BE READ.  THEN, EITHER
!     READ OR SKIP IT.

!     FOR EACH SUBSTRUCTURE/ITEM WHICH IS READ, SAVE THE HEADER RECORD
!     IN OPEN CORE.  WHEN DUPLICATES ARE FOUND, AND THE DATE AND TIME
!     PARAMETERS HAVE NOT BEEN SET, ISSUE A WARNING AND USE THE MOST
!     RECENT.

!     IF THE DATE AND TIME PARAMETERS ARE NON-ZERO, READ ONLY THE
!     SUBSTRUCTURE/ITEM WHICH HAS MATCHING VALUES AND IGNORE THE
!     SUBSTRUCTURE AND ITEM NAME PARAMETERS.

!     READ AN IDENTIFICATION OR HEADER RECORD

570 CALL READ (*580,*590,UNIT,buf,10,1,flag)
GO TO 590
580 IF (names(1) == whole(1) .AND. names(2) == whole(2)) GO TO 1150
CALL REWIND (UNIT)
GO TO 570

!     CHECK IT AGAINST THE FIRST RECORD READ.  IF IT MATCHES, THE ENTIRE
!     TAPE HAS BEEN SCANNED, BUT NOT ALL ITEMS WERE FOUND.

590 DO  i = 1,lrec1
  IF (buf(i) /= hdrec(i)) GO TO 610
END DO
GO TO 1080

!     IF THAT WAS AN ID RECORD, ISSUE MESSAGE AND GO BACK TO READ THE
!     IMMEDIATELY FOLLOWING HEADER RECORD.

610 IF (buf(1) /= id) GO TO 620

!     READ OLD DIT AND MDI DATA

CALL READ (*1870,*1880,UNIT,nos,1,0,flag)
lcore= ncore - 4*nos
idm  = lcore + 1
IF (lcore <= 0) GO TO 1890
nos4 = nos*4
CALL READ (*1870,*1880,UNIT,z(idm),nos4,1,flag)
CALL fwdrec (*1870,UNIT)
line = line + 1
IF (line > nlpp) CALL page
WRITE (nout,2120) uim,(buf(i),i=2,9)
GO TO 570

!     READ OR SKIP THE SUBSTRUCTURE/ITEM DATA.

620 IF (pdate /= 0) GO TO 820
IF (names(1) == whole(1) .AND. names(2) == whole(2)) GO TO 680
DO  i = 1,5
  IF (names(2*i-1) == xxxx) CYCLE
  IF (buf(2) == names(2*i-1) .AND. buf(3) == names(2*i)) GO TO 640
END DO
GO TO 660
640 DO  i = 1,nitems
  IF (buf(4) == xitems(i)) GO TO 680
END DO

!     SKIP -

660 CALL rectyp (UNIT,irec)
IF (irec == 0) GO TO 670

!     STRING RECORD - SKIP IT

CALL fwdrec (*1870,UNIT)
GO TO 660

!     NORMAL GINO RECORD - CHECK IF EOI

670 CALL READ (*1870,*660,UNIT,i,1,1,flag)
IF (i-eoi == 0) THEN
  GO TO   570
ELSE
  GO TO   660
END IF

!     READ -

!     CHECK HEADER RECORDS SAVED IN CORE FOR DUPLICATE

680 IF (iss == 0) GO TO 850
loop700:  DO  i = 1,iss
  jss = 10*(i-1)
  DO  j = 1,3
    IF (buf(j+1) /= z(jss+j)) CYCLE loop700
  END DO
  GO TO 710
END DO loop700
GO TO 850

!     DUPLICATE SUBSTRUCTURE/ITEM ENCOUNTER.  USE MOST RECENT.

710 IF (z(jss+10) /= 0) GO TO 780
line = line+3
IF (line > nlpp) CALL page

!     CHECK YEAR, MONTH, DAY, HOUR, MINUTE, SECOND

IF (z(jss+6)-buf( 7) < 0.0) THEN
  GO TO   800
ELSE IF (z(jss+6)-buf( 7) == 0.0) THEN
  GO TO   720
ELSE
  GO TO   770
END IF
720 IF (z(jss+4)-buf( 5) < 0.0) THEN
  GO TO   800
ELSE IF (z(jss+4)-buf( 5) == 0.0) THEN
  GO TO   730
ELSE
  GO TO   770
END IF
730 IF (z(jss+5)-buf( 6) < 0.0) THEN
  GO TO   800
ELSE IF (z(jss+5)-buf( 6) == 0.0) THEN
  GO TO   740
ELSE
  GO TO   770
END IF
740 IF (z(jss+7)-buf( 8) < 0.0) THEN
  GO TO   800
ELSE IF (z(jss+7)-buf( 8) == 0.0) THEN
  GO TO   750
ELSE
  GO TO   770
END IF
750 IF (z(jss+8)-buf( 9) < 0.0) THEN
  GO TO   800
ELSE IF (z(jss+8)-buf( 9) == 0.0) THEN
  GO TO   760
ELSE
  GO TO   770
END IF
760 IF (z(jss+9)-buf(10) < 0.0) THEN
  GO TO   800
ELSE IF (z(jss+9)-buf(10) == 0.0) THEN
  GO TO   780
END IF

!     MOST RECENT VERSION IS THE ONE ALREADY READ.  THEREFORE, SKIP THE
!     ONE ON TAPE.

770 WRITE (nout,2210) uwm,buf(2),buf(3),buf(4),uname,(buf(i),i=5,10)
780 CALL rectyp (UNIT,irec)
IF (irec == 0) GO TO 790

!     STRING RECORD - SKIP IT

CALL fwdrec (*1870,UNIT)
GO TO 780

!     NORMAL GINO RECORD - CHECK IF EOI

790 CALL READ (*1870,*780,UNIT,i,1,1,flag)
IF (i-eoi == 0) THEN
  GO TO   570
ELSE
  GO TO   780
END IF

!     MOST RECENT VERSION IS ON TAPE.  REPLACE OLDER VERSION ALREADY
!     READ.

800 WRITE (nout,2210) uwm,buf(2),buf(3),buf(4),uname,(z(jss+i),i=4,9)
DO  i = 1,9
  z(jss+i) = buf(i+1)
END DO
jcopy = jcopy - 1
CALL DELETE (buf(2),buf(4),rc)
GO TO 870

!     IF DATE AND TIME PARAMETERS WERE INVOKED, CHECK THEM.

820 IF (MOD(pdate,100)       == buf( 7) .AND.  &
    pdate/10000          == buf( 5) .AND. MOD(pdate,10000)/100 == buf( 6) .AND.  &
    ptime/10000          == buf( 8) .AND. MOD(ptime,10000)/100 == buf( 9) .AND.  &
    MOD(ptime,100)       == buf(10)) GO TO 870

!     DATE AND TIME DONT MATCH. SKIP THIS SUBSTRUCTURE/ITEM.

830 CALL rectyp (UNIT,irec)
IF (irec == 0) GO TO 840

!     STRING RECORD - SKIP IT

CALL fwdrec (*1870,UNIT)
GO TO 830

!     NORMAL GINO RECORD - CHECK IF EOI

840 CALL READ (*1870,*830,UNIT,i,1,1,flag)
IF (i-eoi == 0) THEN
  GO TO   570
ELSE
  GO TO   830
END IF

!     NO DUPLICATE.  ADD THIS HEADER TO THOSE IN CORE.

850 IF (10*(iss+1) > lcore) GO TO 1890
DO  i = 1,9
  z(10*iss+ i) = buf(i+1)
END DO
z(10*iss+10) = 0
iss = iss+1

!     FETCH THE ITEM ON THE SOF.

870 rc  = 3
kdh = ittype(buf(4))
IF (kdh == 1) GO TO 970
CALL sfetch (buf(2),buf(4),swrt,rc)
IF (rc == 3) GO TO 930
line = line + 2
IF (line > nlpp) CALL page
SELECT CASE ( rc )
  CASE (    1)
    GO TO 880
  CASE (    2)
    GO TO 930
  CASE (    3)
    GO TO 930
  CASE (    4)
    GO TO 890
  CASE (    5)
    GO TO 900
END SELECT

!     ITEM ALREADY EXISTS.

880 WRITE (nout,2220) uwm,buf(2),buf(3),buf(4)
z(10*iss) = 1
GO TO 910

!     SUBSTRUCTURE DOES NOT EXIST.  ADD IT TO THE SOF HIERARCHY.

890 CALL exlvl (nos,z(idm),buf(2),z(10*iss+1),lcore-10*iss)
GO TO 870

!     INVALID ITEM NAME

900 CALL smsg (3,buf(4),buf(2))

!     BECAUSE OF ERRORS, NO COPY.  SKIP DATA.

910 CALL rectyp (UNIT,irec)
IF (irec == 0) GO TO 920

!     STRING RECORD - SKIP IT

CALL fwdrec (*1870,UNIT)
GO TO 910

!     NORMAL GINO RECORD - CHECK IF EOI

920 CALL READ (*1870,*910,UNIT,i,1,1,flag)
IF (i-eoi == 0) THEN
  GO TO   570
ELSE
  GO TO   910
END IF

!     COPY THE DATA FROM THE GINO FILE TO THE SOF.

930 i = 10*iss + 1
j = lcore - i + 1
IF (j < 2) GO TO 1890
940 CALL READ (*1870,*950,UNIT,z(i),j,0,flag)
rc = 1
CALL suwrt (z(i),j,rc)
GO TO 940
950 IF (z(i) == eoi) GO TO 960
rc = 2
CALL suwrt (z(i),flag,rc)
GO TO 940
960 rc = 3
CALL suwrt (0,0,rc)
GO TO 1070

!     COPY MATRIX DATA FROM THE GINO FILE TO THE SOF.

970 ifile = scr2
i = 10*iss + 1
j = lcore - i + 1
IF (j < 7) GO TO 1890
CALL READ (*2020,*2030,UNIT,z(i+1),6,1,nw)
z(i) = scr2
CALL wrttrl (z(i))
inblk(1)  = UNIT
outblk(1) = scr2
CALL OPEN (*2010,scr2,z(buf5),wrtrew)
980 CALL rectyp (UNIT,itype)
IF (itype /= 0) GO TO 1010
990 CALL READ (*2010,*1000,UNIT,z(i),j,0,nw)
CALL WRITE (scr2,z(i),j,0)
GO TO 990
1000 IF (z(i) == eoi) GO TO 1020
CALL WRITE (scr2,z(i),nw,1)
GO TO 980
1010 CALL cpystr (inblk,outblk,0,0)
GO TO 980
1020 CALL CLOSE (scr2,1)
1030 CALL mtrxo (scr2,buf(2),buf(4),0,rc)
SELECT CASE ( rc )
  CASE (    1)
    GO TO 1040
  CASE (    2)
    GO TO 1070
  CASE (    3)
    GO TO 1070
  CASE (    4)
    GO TO 1050
  CASE (    5)
    GO TO 1060
  CASE (    6)
    GO TO 2010
END SELECT

!     ITEM ALREADY EXISTS

1040 line = line + 2
IF (line > nlpp) CALL page
WRITE (nout,2220) uwm,buf(2),buf(3),buf(4)
z(10*iss) = 1
GO TO 570

!     SUBSTRUCTURE DOES NOT EXIST - ADD IT TO THE SOF HIERARCHY

1050 CALL exlvl (nos,z(idm),buf(2),z(10*iss+1),lcore-10*iss)
GO TO 1030

!     ILLEGAL ITEM NAME

1060 line = line + 2
IF (line > nlpp) CALL page
CALL smsg (3,buf(4),buf(2))
GO TO 570

!     ITEM COPIED - PRINT MESSAGE

1070 line = line + 1
IF (line > nlpp) CALL page
WRITE (nout,2170) uim,buf(2),buf(3),buf(4),UNIT,sof, (buf(i),i=5,10)

!     INCREMENT NUMBER OF ITEMS COPIED.  IF NOT ALL ARE COPIED, LOOP
!     BACK TO FIND NEXT SUBSTRUCTURE/ITEM ON THE EXTERNAL FILE TO BE
!     COPIED.

jcopy = jcopy + 1
IF (ncopy-jcopy == 0) THEN
  GO TO  1150
ELSE
  GO TO   570
END IF

!     THE ENTIRE EXTERNAL FILE HAS NOW BEEN SCANNED, BUT NOT ALL ITEMS
!     WERE FOUND.  WARN USER OF EACH ITEM NOT FOUND.

!     SKIP REMAINDER OF CURRENT ITEM SO FILE IS PROPERLY POSITIONED
!     FOR NEXT EXECUTION OF MODULE.

1080 DO  i = 1,9,2
  IF (names(i) == xxxx) CYCLE
  loop1110:  DO  item = 1,nitems
    IF (iss == 0) GO TO 1100
    DO  j = 1,iss
      jss = 10*(j-1)
      IF (names(i) == z(jss+1) .AND. names(i+1) == z(jss+2) .AND.  &
          xitems(item) == z(jss+3)) CYCLE loop1110
    END DO
    1100 line = line + 2
    IF (line > nlpp) CALL page
    WRITE (nout,2230) uwm,names(i),names(i+1),xitems(item),uname
  END DO loop1110
END DO
1130 CALL rectyp (UNIT,irec)
IF (irec == 0) GO TO 1140

!     STRING RECORD - SKIP IT

CALL fwdrec (*1870,UNIT)
GO TO 1130

!     NORMAL GINO RECORD - CHECK IF EOI

1140 CALL READ (*1870,*1130,UNIT,i,1,1,flag)
IF (i-eoi == 0) THEN
  GO TO  1150
ELSE
  GO TO  1130
END IF

!     READ OPERATION COMPLETE

1150 CALL CLOSE (UNIT,norew)
GO TO 1740

!     *********************   C H E C K   ***************************

!     REWIND THE EXTERNAL FILE AND PRINT A LIST OF ALL SUBSTRUCTURE/
!     ITEMS ON IT WITH THE DATE AND TIME WHEN THEY WERE WRITTEN THERE.

1160 CALL OPEN (*1860,UNIT,z(buf4),rdrew)
CALL page
WRITE (nout,2240) uim,uname
line = line + 1
CALL READ (*1870,*1880,UNIT,buf,9,1,flag)
GO TO 1180
1170 CALL READ (*1210,*1180,UNIT,buf,10,1,flag)
line = line + 1
IF (line > nlpp) CALL page
WRITE (nout,2250) (buf(i),i=2,10)
GO TO 1190
1180 line = line + 1
IF (line > nlpp) CALL page
WRITE (nout,2120) uim,(buf(i),i=2,9)
1190 CALL rectyp (UNIT,irec)
IF (irec == 0) GO TO 1200

!     STRING RECORD - SKIP IT

CALL fwdrec (*1870,UNIT)
GO TO 1190

!     NORMAL GINO RECORD - CHECK IF EOI

1200 CALL READ (*1870,*1190,UNIT,i,1,1,flag)
IF (i-eoi == 0) THEN
  GO TO  1170
ELSE
  GO TO  1190
END IF
1210 CALL bckrec (UNIT)
CALL CLOSE (UNIT,norew)
GO TO 1740

!     ********************   A P P E N D   ***************************

!     ADD AN EXISTING SOF IN ITS RANDOM ACCESS FORM TO THE RESIDENT SOF.
!     THE MDI AND DIT OF THE EXTERNAL SOF ARE MERGED INTO THOSE OF THE
!     RESIDENT SOF.  THE NXT OF THE EXTERNAL SOF IS INCREMENTED BY THE
!     NUMBER OF BLOCKS IN THE RESIDENT SOF.  THE COMMON BLOCKS /SYS/,
!     /SOF/, AND /SOFCOM/ ARE UPDATED AND WRITTEN TO THE FIRST PHYSICAL
!     BLOCK ON EACH FILE OF THE RESIDENT SOF BY SOFCLS.  NOTE THAT NO
!     USER DATA IS ACTUALLY MOVED.

!     FIRST, ADD THE EXTERNAL SOF TO /SOFCOM/ SO THAT SOFIO CAN BE USED
!     TO READ IT.

1220 IF (nfiles < 10) GO TO 1230
WRITE (nout,2260) uwm,uname
WRITE (nout,2270)
GO TO 1910
1230 nfiles = nfiles + 1
filnam(nfiles) = UNIT
filsiz(nfiles) = 4
nsave = noblks + 1

!     READ THE FIRST PHYSICAL BLOCK OF THE EXTERNAL SOF AND SEE THAT IT
!     IS COMPATIBLE WITH THE RESIDENT SOF.
incblk =-4
DO  i = 1,nfiles
  incblk = incblk+filsiz(i)
END DO
CALL sofio (srd,incblk+1,z(buf4))

!     PASSWORD CHECK

IF (z(buf4+3) == datype(1) .AND. z(buf4+4) == datype(2)) GO TO 1250
WRITE (nout,2260) uwm,uname
WRITE (nout,2310)
incblk =-1

!     FILE SEQUENCE NUMBER CHECK

1250 IF (z(buf4+5) == 1) GO TO 1260
WRITE (nout,2260) uwm,uname
WRITE (nout,2280)
incblk =-1

!     NUMBER OF EXTERNAL FILES CHECK

1260 IF (z(buf4+6) == 1) GO TO 1270
WRITE (nout,2260) uwm,uname
WRITE (nout,2290)
incblk =-1

!     BLOCKSIZE CHECK

1270 IF (z(buf4+27) == blksiz) GO TO 1280
WRITE (nout,2260) uwm,uname
WRITE (nout,2300) blksiz,z(buf4+27)
incblk =-1
1280 IF (incblk < 0) GO TO 1490

!     COMPLETE THE UPDATING OF THE COMMON BLOCKS

filsiz(nfiles) = z(buf4+17)
avblks = avblks + z(buf4+30)
nxtcur = 1
nxtrst =.true.
nxtfsz(nfiles) = z(buf4+36)
j = nfiles-1
nxttsz = 0
DO  i = 1,j
  nxttsz = nxttsz + nxtfsz(i)
END DO
oldtsz = nxttsz + 1
nxttsz = nxttsz + z(buf4+35)

!     READ THE DIT OF THE EXTERNAL SOF AND ADD EACH SUBSTRUCTURE THERE
!     TO THE DIT OF THE RESIDENT SOF.  KEEP A TABLE IN OPEN CORE OF TWO
!     WORDS PER SUBSTRUCTURE -

!     (1)  SUBSTRUCTURE NUMBER FROM THE EXTERNAL SOF.
!     (2)  NEW SUBSTRUCTURE NUMBER ON THE RESIDENT SOF.

nold = z(buf4+32)
IF (2*nold > lcore) GO TO 1890
iss = 1
k   = 1
kdit = z(buf4+33) + incblk
kmdi = z(buf4+34) + incblk
1300 CALL sofio (srd,kdit,z(buf4))
DO  i = 1,blksiz,2
  ssname(1) = z(buf4+i+2)
  ssname(2) = z(buf4+i+3)
  IF (ssname(1) == BLANK) CYCLE
  1320 CALL fdsub (ssname,j)
  IF (j == -1) GO TO 1330
  
!     DUPLICATE NAME ON RESIDENT SOF.  PREFIX IT WITH -Q- AND TRY AGAIN.
  
  WRITE (nout,2320) uwm,ssname
  CALL prefix (q,ssname)
  IF (ssname(2) /= qqqq) GO TO 1320
  WRITE (nout,2330)
  z(iss  ) = (i+1)/2
  z(iss+1) = 0
  iss = iss + 2
  GO TO 1340
  1330 CALL crsub (ssname,j)
  z(iss  ) = k
  z(iss+1) = j
  iss = iss + 2
  k   = k + 1
  1340 IF (iss/2 >= nold) GO TO 1380
END DO

!     GET THE NEXT BLOCK OF THE DIT FROM THE EXTERNAL SOF

CALL fnxt (kdit,j)
IF (MOD(kdit,2) == 1) GO TO 1360
i = andf(rshift(cor(j),ihalf),jhalf)
GO TO 1370
1360 i = andf(cor(j),jhalf)
1370 kdit = i + incblk
GO TO 1300

!     THE DIT OF THE EXTERNAL SOF HAS NOW BEEN MERGED WITH THE DIT OF
!     THE RESIDENT SOF.  NOW MERGE THE MDI

1380 iss = 0
1390 CALL sofio (srd,kmdi,z(buf4))
DO  i = 1,blksiz,dirsiz
  IF (blksiz-i+1 < dirsiz) CYCLE
  iss = iss  + 1
  jmdi= buf4 + i + 1
  CALL bisloc (*1900,iss,z,2,nold,k)
  CALL fmdi (z(k+1),jrmdi)
  
!     PUT THE CONVERTED SUBSTRUCTURE INDICES IN THE FIRST TWO WORDS OF
!     THE MDI OF THE RESIDENT SOF.
  
  DO  j = 1,6
    mask = lshift(1023,10*((j-1)/2))
!                   1023 = 2*10-1, LEFT SHIFT 0, 10, AND 20 BITS
    
    k = MOD(j-1,2) + 1
    jss = andf(z(jmdi+k),mask)
    IF (jss == 0) CYCLE
    CALL bisloc (*1900,jss,z,2,nold,k)
    jss = z(k+1)
    cor(jrmdi+k) = andf(cor(jrmdi+k),lshift(jss,10*((j-1)/2)))
  END DO
  
!     INCREMENT THE BLOCK INDICES OF THE ITEMS IN THIS MDI DIRECTORY BY
!     THE NUMBER OF BLOCKS ON THE RESIDENT SOF.
  
  DO  j = ifrst,dirsiz
    IF (andf(z(jmdi+j),jhalf) == 0) CYCLE
    cor(jrmdi+j) = z(jmdi+j) + incblk
  END DO
  IF (iss == nold) GO TO 1450
END DO

!     GET THE NEXT BLOCK OF THE MDI FROM THE EXTERNAL SOF.

CALL fnxt (kmdi,j)
IF (MOD(kmdi,2) == 1) GO TO 1430
i = andf(rshift(cor(j),ihalf),jhalf)
GO TO 1440
1430 i = andf(cor(j),jhalf)
1440 kmdi = i + incblk
GO TO 1390

!     THE MDI OF THE EXTERNAL SOF HAS NOW BEEN MERGED WITH THE MDI OF
!     THE RESIDENT SOF.  NOW UPDATE THE NXT OF THE EXTERNAL SOF.

1450 n = blksiz
knxt = incblk + 2
incblk = orf(incblk,lshift(incblk,ihalf))
DO  i = oldtsz,nxttsz
  CALL sofio (srd,knxt,z(buf4))
  IF (i-oldtsz+1 == nxtfsz(nfiles))  &
      n = (MOD(filsiz(nfiles)-2,supsiz)+1)/2 + 1
  DO  j = 1,n
    z(buf4+j+2) = z(buf4+j+2) + incblk
  END DO
  CALL sofio (swrt,knxt,z(buf4))
  knxt = knxt + supsiz
END DO

!     RELEASE THE BLOCKS USED BY THE MDI AND DIT OF THE EXTERNAL SOF.
!     (THIS WILL CAUSE THE EXTERNAL SOF TO BE UNUSEABLE IN ITS ORIGINAL
!     FORM.)

incblk = andf(incblk,jhalf)
CALL sofio (srd,incblk+1,z(buf4))
kdit = z(buf4+33) + incblk
kmdi = z(buf4+34) + incblk
CALL retblk (kdit)
CALL retblk (kmdi)

!     WRITE ON ALL BLOCKS BETWEEN THE HIGHEST BLOCK WRITTEN ON THE
!     ORIGINAL RESIDENT SOF AND THE FIRST BLOCK OF THE APPENDED SOF.
!     THIS IS REQUIRED TO AVOID DATA TRANSMISSION ERRORS.

n = filsiz(nfiles-1)
DO  i = nsave,n
  CALL sofio (swrt,nsave,z(buf4))
END DO

!     SOFCLS WILL UPDATE THE FIRST PHYSICAL BLOCK ON EACH SOF UNIT.

CALL sofcls

!     APPEND OPERATION COMPLETED SUCCESSFULLY.  TELL USER THE NEWS.

WRITE (nout,2340) uim,uname
n = sofsiz(n)
WRITE (nout,2360) uim,avblks,n
GO TO 1750

!     APPEND OPERATION ABORTED.  RESTORE THE COMMON BLOCKS FOR THE
!     RESIDENT SOF.

1490 first =.true.
opnsof=.false.
CALL sofopn (z(buf1),z(buf2),z(buf3))
GO TO 1900

!     ********************   C O M P R E S S   **********************

!     FOR EACH SUBSTRUCTURE IN THE DIT, COPY EACH ITEM WHICH EXISTS OR
!     PSEUDO-EXISTS TO SCR1 AND DELETE THE ITEM ON THE SOF.  THEN COPY
!     ALL ITEMS BACK.  ALL INTERMEDIATE FREE BLOCKS WILL THUS BE
!     ELIMINATED AND THE DATA FOR ANY ONE ITEM WILL BE STORED ON
!     CONTIGUOUS BLOCKS.

!     THE FORMAT OF THE SCRATCH FILE IS --

!                                      +------------+
!     SUBSTRUCTURE NAME (2 WORDS)      I            I+
!     ITEM NAME (1 WORD)               I HEADER     I +
!     PSEUDO FLAG -- 2 FOR PSEUDO-ITEM I RECORD     I  +
!                    3 FOR REAL DATA   I            I   +   REPEATED
!                                      +------------+    +  FOR EACH
!     DATA -- 1 SOF GROUP PER          I DATA       I   +   SUBS./ITEM
!             GINO LOGICAL RECORD      I RECORDS    I  +
!                                      +------------+ +
!     END OF ITEM FLAG (1 WORD)        I EOI RECORD I+
!                                      +------------+

1500 UNIT = scr1
CALL OPEN (*1860,scr1,z(buf4),wrtrew)

!     COPY OUT DIT AND MDI INFORMATION

iss = 0
DO  k = 1,ditsiz,2
  iss = iss + 1
  CALL fdit  (iss,j)
  CALL WRITE (scr1,cor(j),2,0)
  CALL fmdi  (iss,j)
  CALL WRITE (scr1,cor(j+1),2,0)
END DO
CALL WRITE (scr1,0,0,1)

!     COPY OUT SUBSTRUCTURE ITEMS

iss = 0
DO  k = 1,ditsiz,2
  iss = iss + 1
  CALL fdit (iss,j)
  ssname(1) = cor(j  )
  ssname(2) = cor(j+1)
  IF (ssname(1) == BLANK) CYCLE
  DO  item = 1,nitem
    kdh = items(2,item)
    IF (kdh == 1) GO TO 1570
    CALL sfetch (ssname,items(1,item),srd,rc)
    SELECT CASE ( rc )
      CASE (    1)
        GO TO 1540
      CASE (    2)
        GO TO 1530
      CASE (    3)
        GO TO 1590
      CASE (    4)
        GO TO 1520
      CASE (    5)
        GO TO 1520
    END SELECT
    1520 CALL smsg (rc-2,items(1,item),ssname)
    CYCLE
    
!     ITEM PSEUDO-EXISTS.  WRITE PSEUDO-HEADER RECORD AND EOI RECORD.
    
    1530 CALL WRITE (scr1,ssname,2,0)
    CALL WRITE (scr1,items(1,item),1,0)
    CALL WRITE (scr1,2,1,1)
    CALL WRITE (scr1,eoi,1,1)
    CYCLE
    
!     ITEM EXISTS.  COPY IT OUT.
    
    1540 CALL WRITE (scr1,ssname,2,0)
    CALL WRITE (scr1,items(1,item),1,0)
    CALL WRITE (scr1,3,1,1)
    1550 CALL suread(z,lcore,n,rc)
    IF (rc > 1) GO TO 1560
    CALL WRITE (scr1,z,lcore,0)
    GO TO 1550
    1560 CALL WRITE (scr1,z,n,1)
    IF (rc == 2) GO TO 1550
    
!     END OF ITEM HIT.  WRITE EOI RECORD
    
    CALL WRITE (scr1,eoi,1,1)
    CYCLE
    
!     PROCESS MATRIX ITEMS
    
    1570 CALL mtrxi (scr2,ssname,items(1,item),0,rc)
    ifile = scr2
    SELECT CASE ( rc )
      CASE (    1)
        GO TO 1580
      CASE (    2)
        GO TO 1530
      CASE (    3)
        GO TO 1590
      CASE (    4)
        GO TO 1520
      CASE (    5)
        GO TO 1520
      CASE (    6)
        GO TO 2010
    END SELECT
    1580 CALL WRITE (scr1,ssname,2,0)
    CALL WRITE (scr1,items(1,item),1,0)
    CALL WRITE (scr1,3,1,1)
    CALL OPEN  (*2010,scr2,z(buf5),rdrew)
    z(1) = scr2
    CALL rdtrl (z(1))
    CALL WRITE (scr1,z(iz2),6,1)
    CALL cpyfil(scr2,scr1,z,lcore,icount)
    CALL WRITE (scr1,eoi,1,1)
    CALL CLOSE (scr2,1)
    1590 CONTINUE
  END DO
END DO

!     COPY ALL ITEMS BACK TO THE SOF

CALL CLOSE (scr1,rew)
CALL OPEN  (*1860,scr1,z(buf4),rdrew)

!     RE-INITIALIZE THE SOF, THEN RESTORE THE OLD DIT AND MDI

CALL sofcls
STATUS= 0
first =.true.
CALL sofopn (z(buf1),z(buf2),z(buf3))
CALL page
iss = 0
1610 CALL READ (*1870,*1620,scr1,buf,4,0,flag)
iss = iss + 1
IF (buf(1) == BLANK) GO TO 1610
CALL crsub (buf,i)
CALL fmdi  (i,j)
cor(j+1) = buf(3)
cor(j+2) = buf(4)
mdiup = .true.
GO TO 1610

!     READ HEADER RECORD AND FETCH THE SOF ITEM

1620 CALL READ (*1730,*1880,scr1,buf,4,1,flag)
kdh = ittype(buf(3))
IF (kdh == 1) GO TO 1660
CALL sfetch (buf,buf(3),2,buf(4))

!     COPY THE DATA

1630 CALL READ (*1870,*1640,scr1,z,lcore,0,flag)
IF (z(1) == eoi) GO TO 1650
CALL suwrt (z,lcore,1)
GO TO 1630
1640 IF (z(1) == eoi) GO TO 1650
CALL suwrt (z,flag,2)
GO TO 1630

!     EOI FOUND

1650 CALL suwrt (0,0,3)
GO TO 1720

!     COPY IN MATRIX ITEMS

1660 CALL OPEN (*2010,scr2,z(buf5),wrtrew)
CALL READ (*1870,*1880,scr1,z(iz2),6,1,nw)
z(1) = scr2
CALL wrttrl (z(1))
inblk(1)  = scr1
outblk(1) = scr2
1670 CALL rectyp (scr1,itype)
IF (itype /= 0) GO TO 1700
1680 CALL READ (*1870,*1690,scr1,z,lcore,0,nw)
CALL WRITE (scr2,z,lcore,0)
GO TO 1680
1690 IF (z(1) == eoi) GO TO 1710
CALL WRITE (scr2,z,nw,1)
GO TO 1670
1700 CALL cpystr (inblk,outblk,0,0)
GO TO 1670

!     EOI FOUND

1710 CALL CLOSE (scr2,1)
CALL mtrxo (scr2,buf,buf(3),0,rc)
1720 CONTINUE
line = line + 1
IF (line > nlpp) CALL page
WRITE (nout,2350) uim,buf(1),buf(2),buf(3)
GO TO 1620

!     COMPRESS COMPLETE

1730 CALL CLOSE (scr1,rew)

!     **********************   C O D A   ************************

!     NORMAL TERMINATION

1740 CALL sofcls
1750 RETURN

!     ERRORS CAUSING MODULE AND/OR JOB TERMINATION

1800 WRITE (nout,2100) uwm,uname
GO TO 1910
1810 WRITE (nout,2110) uwm,device
GO TO 1910
1820 WRITE (nout,2140) uwm,mode
GO TO 1910
1830 WRITE (nout,2150) uwm,pos
GO TO 1910
1840 WRITE (nout,2180) uwm
GO TO 1910
1850 WRITE (nout,2190) swm,uname
CALL CLOSE (UNIT,norew)
GO TO 1910

1860 n = -1
GO TO 2000
1870 n = -2
GO TO 2000
1880 n = -3
GO TO 2000
1890 n = 8
GO TO 2000
1900 n = -61
GO TO 2000
1910 CALL sofcls
dry = -2
WRITE (nout,2370) sim
RETURN

2000 CALL sofcls
CALL mesage (n,UNIT,subr)
dry = -2
WRITE (nout,2370) sim
RETURN

2010 n = -1
GO TO 2040
2020 n = -2
GO TO 2040
2030 n = -3
2040 CALL sofcls
CALL mesage (n,ifile,subr)
RETURN

!     TEXT OF ERROR MESSAGES

2100 FORMAT (a25,' 6334, EXIO DEVICE PARAMETER SPECIFIES TAPE, BUT ',  &
    'UNIT ',2A4,' IS NOT A PHYSICAL TAPE')
2110 FORMAT (a25,' 6335, ',2A4,' IS AN INVALID DEVICE FOR MODULE EXIO')
2120 FORMAT (a29,' 6336, EXIO FILE IDENTIFICATION.  PASSWORD= ',2A4,  &
    '  DATE=',i3,1H/,i2,1H/,i2,7H  time=,i3,1H.,i2,1H.,i2)
2130 FORMAT (a29,' 6337,',i6,' BLOCKS (',i4,' SUPERBLOCKS) OF THE SOF',  &
    ' SUCCESSFULLY DUMPED TO EXTERNAL FILE ',2A4)
2140 FORMAT (a25,' 6338, ',2A4,' IS AN INVALID MODE PARAMETER FOR ',  &
    'MODULE EXIO')
2150 FORMAT (a25,' 6339, ',2A4,' IS AN INVALID FILE POSITIONING ',  &
    'PARAMETER FOR MODULE EXIO')
2160 FORMAT (a25,' 6340, SUBSTRUCTURE ',2A4,' ITEM ',a4,  &
    ' PSEUDOEXISTS ONLY AND CANNOT BE COPIED OUT BY EXIO')
2170 FORMAT (a29,' 6341, SUBSTRUCTURE ',2A4,' ITEM ',a4,  &
    ' SUCCESSFULLY COPIED FROM ',a4,' TO ',a4,2H (,  &
    i2,1H/,i2,1H/,i2,2H, ,i2,1H.,i2,1H.,i2,1H))
2180 FORMAT (a25,' 6342, SOF RESTORE OPERATION FAILED.  THE RESIDENT ',  &
    'SOF IS NOT EMPTY')
2190 FORMAT (a27,' 6343, ',2A4,' IS NOT AN EXTERNAL SOF FILE')
2200 FORMAT (a29,' 6344, SOF RESTORE OF ',i6,' BLOCKS SUCCESSFULLY ',  &
    'COMPLETED')
2210 FORMAT (a25,' 6345, SUBSTRUCTURE ',2A4,' ITEM ',a4,  &
    ' IS DUPLICATED ON EXTERNAL FILE ',2A4, /32X,  &
    'OLDER VERSION (',i2,1H/,i2,1H/,i2,2H, ,i2,1H.,i2,1H.,i2, ') IS IGNORED')
2220 FORMAT (a25,' 6346, SUBSTRUCTURE ',2A4,' ITEM ',a4,  &
    ' NOT COPIED.  IT ALREADY EXISTS ON THE SOF')
2230 FORMAT (a25,' 6348, SUBSTRUCTURE ',2A4,' ITEM ',a4,  &
    ' NOT FOUND ON EXTERNAL FILE ',2A4)
2240 FORMAT (a29,' 6349, CONTENTS OF EXTERNAL SOF FILE ',2A4,' FOLLOW')
2250 FORMAT (5X,'SUBSTRUCTURE ',2A4,5X,'ITEM ',a4,10X,5HDATE ,i2,1H/,  &
    i2,1H/,i2,10X,5HTIME ,i2,1H.,i2,1H.,i2)
2260 FORMAT (a25,' 6350, SOF APPEND OF FILE ',2A4,' FAILED')
2270 FORMAT (32X,'TOO MANY PHYSICAL SOF UNITS. MAXIMUM ALLOWED IS 10')
2280 FORMAT (32X,'THE SEQUENCE NUMBER OF THE EXTERNAL SOF FILE IS NOT', ' 1')
2290 FORMAT (32X,'THE EXTERNAL SOF FILE MUST CONSIST OF ONLY ONLY ONE',  &
    ' PHYSICAL UNIT')
2300 FORMAT(32X,45HTHE EXTERNAL sof has incompatible BLOCK size., /  &
    32X, 32HBLOCK size of the resident sof = ,i5, /  &
    32X, 32HBLOCK size of the EXTERNAL sof = ,i5 )
2310 FORMAT (32X,17HINVALID password.)
2320 FORMAT (a25,' 6351, DUPLICATE SUBSTRUCTURE NAME ',2A4,  &
' FOUND DURING SOF APPEND OF FILE ',2A4, /32X,  &
    'THE SUBSTRUCTURE WITH THIS NAME ON THE FILE BEING ',  &
    'APPENDED WILL BE PREFIXED WITH Q')
2330 FORMAT (1H0,31X, 37HPREFIX failed.  substructure ignored.)
2340 FORMAT (a29,' 6352, EXTERNAL SOF FILE ',2A4,  &
    ' SUCCESSFULLY APPENDED TO THE RESIDENT SOF')
2350 FORMAT (a29,' 6353, SUBSTRUCTURE ',2A4,' ITEM ',a4,  &
    ' HAS BEEN SUCCESSFULLY COMPRESSED')
2360 FORMAT (a29,' 6354, THERE ARE',i7,' FREE BLOCKS (',i9,  &
    ' WORDS) ON THE RESIDENT SOF')
2370 FORMAT (a31,' 6355, EXIO TERMINATED WITH ERRORS.  DRY RUN MODE ',  &
    'ENTERED')
END SUBROUTINE exio1
