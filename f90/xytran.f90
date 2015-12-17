SUBROUTINE xytran
     
 IMPLICIT INTEGER (a-z)
 LOGICAL :: vgp,random,outopn,PRINT,plot,paplot,oompp,oomcp, punch
 INTEGER :: word(58),namev(11),files(11),subcas(200),  &
     NAME(2),majid(11),routin(2),headsv(96), xycard(20),openf(5),indb(5)
 REAL :: temp,temp1,rbuf(100),rz(1),value(60)
 CHARACTER (LEN=25) :: uwm
 CHARACTER (LEN=23) :: ufm
 COMMON /xmssg / ufm,uwm
 COMMON /BLANK / blkcom,idum1,ipset,ipset2,nframe,ncard
 COMMON /system/ sysbuf,l,nogo,nin,ksystm(81),intr
 COMMON /output/ ihead(96)
 COMMON /zzzzzz/ z(1)
 COMMON /xywork/ FILE,tcurve(32),ntops,PRINT,ifile,xaxis(32),  &
     nbots,plot,vector,yaxis(32),vecid(5),punch,  &
     major,ytaxis(32),subc(5),center,random,ybaxis(32),  &
     idin(153),buf(100),ivalue(60),iat,idout(300), outopn,steps,nat,paplot,knt
 EQUIVALENCE     (z(1),rz(1)),(buf(1),rbuf(1)),(ivalue(1),value(1))
 
 DATA    STOP  / 4HSTOP /, GO    / 4HGO   /, vdum  / 4HVDUM /,  &
     xy    / 4HXY   /, fram  / 4HFRAM /, clea  / 4HCLEA /,  &
     tcur  / 4HTCUR /, xaxi  / 4HXTIT /, yaxi  / 4HYTIT /,  &
     ytax  / 4HYTTI /, ybax  / 4HYBTI /, BLANK / 4H     /, pset  / 4HPSET /
 
 DATA    eor   / 1  /, noeor /0  /, outrwd/1/, inprwd/0/, rewd/1/
 DATA    xycdb / 101/, outfil/201/, indb  / 102,103,104,105,106 /
 DATA    nwords/ 58 /, routin/4HXYTR, 4HAN  /, rand  / 4HRAND   /
 DATA    vg    / 4HVG  /, i3 / 3 /
 
 DATA    word  / 4HXMIN, 4HXMAX, 4HYMIN, 4HYMAX, 4HYTMI, 4HYTMA, 4HYBMI,  &
     4HYBMA, 4HXINT, 4HYINT, 4HYTIN, 4HYBIN, 4HXAXI, 4HYAXI,  &
     4HXTAX, 4HXBAX, 4HXDIV, 4HYDIV, 4HYTDI, 4HYBDI, 4HXVAL,  &
     4HYVAL, 4HYTVA, 4HYBVA, 4HUPPE, 4HLOWE, 4HLEFT, 4HRIGH,  &
     4HTLEF, 4HTRIG, 4HBLEF, 4HBRIG, 4HALLE, 4HTALL, 4HBALL,  &
     4HXLOG, 4HYLOG, 4HYTLO, 4HYBLO, 4HCURV, 4HDENS, 4H....,  &
     4H...., 4H...., 4HSKIP, 4HCAME, 4HPLOT, 4HXPAP, 4HYPAP,  &
     4HPENS, 4HXGRI, 4HYGRI, 4HXTGR, 4HYTGR, 4HXBGR, 4HYBGR, 4HCSCA, 4HCOLO/
 
!     DATA FOR THE 11 VECTOR TYPES POSSIBLE
 
!                                                         BASIC
!              VECTOR-NAME         RESIDENT-FILE       MAJOR - ID
!          ******************     ***************   ****************
 DATA namev( 1) / 4HDISP /,  files( 1) / 3 /,  majid( 1) /  1 /
 DATA namev( 2) / 4HVELO /,  files( 2) / 3 /,  majid( 2) / 10 /
 DATA namev( 3) / 4HACCE /,  files( 3) / 3 /,  majid( 3) / 11 /
 DATA namev( 4) / 4HSPCF /,  files( 4) / 2 /,  majid( 4) /  3 /
 DATA namev( 5) / 4HLOAD /,  files( 5) / 1 /,  majid( 5) /  2 /
 DATA namev( 6) / 4HSTRE /,  files( 6) / 4 /,  majid( 6) /  5 /
 DATA namev( 7) / 4HFORC /,  files( 7) / 5 /,  majid( 7) /  4 /
 DATA namev( 8) / 4HSDIS /,  files( 8) / 1 /,  majid( 8) / 15 /
 DATA namev( 9) / 4HSVEL /,  files( 9) / 1 /,  majid( 9) / 16 /
 DATA namev(10) / 4HSACC /,  files(10) / 1 /,  majid(10) / 17 /
 DATA namev(11) / 4HNONL /,  files(11) / 2 /,  majid(11) / 12 /
 DATA namevg    / 4H vg  /
 
!     - IDOUT DATA RECORD DISCRIPTION -
 
!     WORD    TYPE   DISCRIPTION
!     ==================================================================
!       1     I/R    SUBCASE ID OR IF RANDOM THE MEAN RESPONSE
!       2      I     FRAME NUMBER
!       3      I     CURVE NUMBER
!       4      I     POINT-ID OR ELEMENT-ID
!       5      I     COMPONENT NUMBER
!       6      I     VECTOR NUMBER  1 THRU 11
 
!       7      I     1 -- CURVE USES TOP HALF OF FRAME
!                    0 -- CURVE USES FULL FRAME
!                   -1 -- CURVE USES LOWER HALF OF FRAME
 
!       8      I     0 -- AXIS,TICS,LABELS,VALUES, ETC. HAVE BEEN DRAWN
!                         AND THIS CURVE IS TO BE SCALED AND PLOTTED
!                         IDENTICALLY AS LAST EXCEPT FOR CURVE SYMBOLS.
!                    1 -- AXIS, TICS, LABELS, SCALEING, ETC. ARE TO BE
!                         PERFORMED OR COMPUTED AND IF IDOUT(7)=0 OR 1
!                         A SKIP TO NEW FRAME IS TO BE MADE.
 
!       9      I     NUMBER OF BLANK FRAMES BETWEEN FRAMES (FRAME-SKIP)
!      10      R     MINIMUM X-INCREMENT
!      11      R     XMIN  *
!      12      R     XMAX   *   DEFINES ACTUAL LIMITS OF DATA OF THIS
!      13      R     YMIN   *   UPPER, LOWER, OR FULL FRAME CURVE.
!      14      R     YMAX  *
!      15      R     ACTUAL VALUE OF FIRST TIC                 *
!      16      R     ACTUAL INCREMENT TO SUCCESSIVE TICS        *
!      17      I     ACTUAL MAXIMUM VALUE OF FRAME               *  X-
!      18      I     MAXIMUM NUMBER OF DIGITS IN ANY PRINT-VALUE  * DIRE
!      19      I     + OR - POWER FOR PRINT VALUES                * TICS
!      20      I     TOTAL NUMBER OF TICS TO PRINT THIS EDGE     *
!      21      I     VALUE PRINT SKIP  0,1,2,3---               *
!      22      I     SPARE                                     *
!      23      R     *
!      24      R      *
!      25      I       *
!      26      I        *  SAME AS  15 THRU 22
!      27      I        *  BUT FOR  Y-DIRECTION TICS
!      28      I       *
!      29      I      *
!      30      I     *
!      31      I     TOP EDGE TICS   **   EACH OF 31 THRU 34 MAY BE
!      32      I     BOTTOM EDGE TICS **  LESS THAN 0 -- TICS W/O VALUES
!      33      I     LEFT EDGE TICS   **  EQUAL TO  0 -- NO TICS HERE
!      34      I     RIGHT EDGE TICS **   GREATER   0 -- TICS W VALUES
 
!      35      I     0 -- X-DIRECTION IS LINEAR
!                    GREATER THAN 0 - NUMBR OF CYCLES AND X-DIREC IS LOG
!      36      I     0 -- Y-DIRECTION IS LINEAR
!                    GREATER THAN 0 - NUMBR OF CYCLES AND Y-DIREC IS LOG
!      37      I     0 -- NO X-AXIS
!                    1 -- DRAW X-AXIS
 
!      38      R     X-AXIS  Y-INTERCEPT
 
!      39      I     0 -- NO Y-AXIS
!                    1 -- DRAW Y-AXIS
 
!      40      R     Y-AXIS  X-INTERCEPT
 
!      41      I     LESS THAN 0 ----- PLOT SYMBOL FOR EACH CURVE POINT.
!                                      SELECT SYMBOL CORRESPONDING TO
!                                      CURVE NUMBER IN IDOUT(3)
!                    EQUAL TO  0 ----- CONNECT POINTS BY LINES WHERE
!                                      POINTS ARE CONTINUOUS I.E.(NO
!                                      INTEGER 1 PAIRS)
!                    GREATER THAN 0 -- DO BOTH OF ABOVE
 
!      42
!       .
!       .
!      50
!      51     BCD    TITLE(32)
!       .     BCD    SUBTITLE(32)
!       .     BCD    LABEL(32)
!       .     BCD    CURVE TITLE(32)
!       .     BCD    X-AXIS TITLE(32)
!     242     BCD    Y-AXIS TITLE(32)
!     243      I     XGRID LINES   0=NO   1=YES
!     244      I     YGRID LINES   0=NO   1=YES
!     245      I     TYPE OF PLOT  1=RESPONSE, 2=PSDF, 3=AUTO
!     246      I     STEPS
!       .
!       .
!     281      I     PAPLOT FRAME NUMBER
!     282      R     CSCALE (REAL NUMBER)
!     283      I     PENSIZE OR DENSITY
!     284      I     PLOTTER (LSHIFT 16) AND MODEL NUMBER.
!     285      R     INCHES PAPER X-DIRECTION
!     286      R     INCHES PAPER Y-DIRECTION
!     287      I     CAMERA FOR SC4020 LESS THAN 0=35MM, 0=F80,
!                                        GREATER 0=BOTH
!     288      I     PRINT FLAG  **
!     289      I     PLOT  FLAG  ** 0=NO, +=YES (PLOT- 2=BOTH, -1=PAPLT)
!     290      I     PUNCH FLAG  **
!     291      R     X-MIN OF ALL DATA
!     292      R     X-MAX OF ALL DATA
!     293      R     Y-MIN WITHIN X-LIMITS OF FRAME
!     294      R     X-VALUE AT THIS Y-MIN
!     295      R     Y-MAX WITHIN X-LIMITS OF FRAME
!     296      R     X-VALUE AT THIS Y-MAX
!     297      R     Y-MIN FOR ALL DATA
!     298      R     X-VALUE AT THIS Y-MIN
!     299      R     Y-MAX FOR ALL DATA
!     300      R     X-VALUE AT THIS Y-MAX
!     ==================================================================
 
!     SAVE OUTPUT HEADING
 
 DO  i = 1,96
   headsv(i) = ihead(i)
 END DO
 
!     ALLOCATE CORE AND OPEN DATA BLOCKS
 
 oompp = .false.
 vgp   = .false.
 oomcp = .false.
 random= .false.
 ifle  = xycdb
 core  = korsz(z) - 1
 DO  i = 1,32
   tcurve(i) = BLANK
   xaxis(i)  = BLANK
   yaxis(i)  = BLANK
   ytaxis(i) = BLANK
   ybaxis(i) = BLANK
 END DO
 DO  i = 1,5
   subc(i) = 1
 END DO
 nsubs  = 0
 core   = core - sysbuf
 IF (core < 0) GO TO 825
 intrwd = inprwd
 IF (intr <= 0) GO TO 35
 intrwd = outrwd
 xycdb  = 301
 35 CALL OPEN (*835,xycdb,z(core+1),intrwd)
 IF (intr <= 0) GO TO 65
 card = 1
 WRITE (l,900)
 40 DO  ij = 1,20
   xycard(ij) = BLANK
 END DO
 CALL xread (*43,xycard)
 IF (xycard(1) == STOP) GO TO 1500
 IF (xycard(1) ==   GO) card = -1
 CALL ifp1xy (card,xycard)
 IF (xycard(1) == GO) GO TO 50
 card = 0
 IF (nogo == 0) GO TO 45
 nogo = 0
 43 WRITE (l,902)
 45 WRITE (l,910) xycard
 GO TO 40
 50 CALL CLOSE (xycdb,rewd)
 IF (intr > 10) l = 1
 CALL OPEN (*835,xycdb,z(core+1),inprwd)
 65 IF (intr <= 0) CALL fwdrec (*80,xycdb)
 outopn = .false.
 IF (blkcom == rand) random = .true.
 IF (blkcom ==   vg) vgp = .true.
 IF (blkcom ==   vg) namev(5) = namevg
 
 core = core - sysbuf
 DO  i = 1,5
   openf(i) = -1
   IF (core < 0) GO TO 235
   
   CALL OPEN (*70,indb(i),z(core),inprwd)
   openf(i) = 0
   vecid(i) = 0
   core = core - sysbuf
   70 CONTINUE
 END DO
 
 core = core + sysbuf - 1
 
!     NOTE - OUTPUT DATA BLOCKS WILL BE OPENED WHEN AND IF REQUIRED
 
 
 
!     READ FIRST BCD WORD FROM -XYCDB- THEN GO INITIALIZE DATA
 
 bcd = clea
 GO TO 800
 80 ier = 2
 GO TO 237
 90 ier = 3
 GO TO 237
 
!     BRANCH ON BCD WORD
 
 100 IF (bcd == xy  ) GO TO 230
 IF (bcd == tcur) GO TO 180
 IF (bcd == xaxi) GO TO 190
 IF (bcd == yaxi) GO TO 200
 IF (bcd == ytax) GO TO 210
 IF (bcd == ybax) GO TO 220
 
!     SET SINGLE VALUE FLAGS. READ IN VALUE
 
 IF (bcd == clea) GO TO 150
 IF (bcd == vdum) GO TO 820
 CALL READ (*80,*90,xycdb,ival,1,noeor,flag)
 DO  i = 1,nwords
   IF (bcd == word(i)) GO TO 130
 END DO
 
!     WORD NOT RECOGNIZED
 
 CALL page2 (2)
 WRITE (l,120) uwm,bcd
 GO TO 140
 
!     KEY WORD FOUND
 
 130 IF (bcd /= word(58)) GO TO 135
 ivalue(i) = ival
 CALL READ (*80,*90,xycdb,ival,1,noeor,flag)
 ivalue(i+1) = ival
 GO TO 140
 135 ivalue(i) = ival
 
!     READ NEXT BCD WORD
 
 140 CALL READ (*80,*240,xycdb,bcd,1,noeor,flag)
 GO TO 100
 
!     CLEAR ALL VALUES SET AND RESTORE DEFAULTS
 
 150 DO  i = 1,12
   ivalue(i) = 1
 END DO
 DO  i = 13,nwords
   IF (i /= 47) ivalue(i) = 0
 END DO
 DO  i = 25,32
   ivalue(i) = 1
 END DO
 
!     DEFAULT CAMERA TO BOTH
 
 ivalue(46) = 3
 GO TO 140
 
!     SET TITLES
 
 180 CALL READ (*80,*90,xycdb,tcurve(1),32,noeor,flag)
 GO TO 140
 190 CALL READ (*80,*90,xycdb,xaxis(1),32,noeor,flag)
 GO TO 140
 200 CALL READ (*80,*90,xycdb,yaxis(1),32,noeor,flag)
 GO TO 140
 210 CALL READ (*80,*90,xycdb,ytaxis(1),32,noeor,flag)
 GO TO 140
 220 CALL READ (*80,*90,xycdb,ybaxis(1),32,noeor,flag)
 GO TO 140
 
!     XY-COMMAND OPERATIONS HIT
 
 230 CALL READ (*80,*90,xycdb,buf(1),7,noeor,flag)
 IF (buf(6) /= 0) paplot = .true.
 IF (buf(6) /= 0) oompp  = .true.
 IF (buf(2) /= 0) oomcp  = .true.
 IF (buf(1) /= 0) PRINT  = .true.
 IF (buf(2) /= 0) plot   = .true.
 kasknt = 0
 IF (outopn) GO TO 280
 IF (.NOT.plot .AND. .NOT.paplot) GO TO 280
 
!     OPEN OUTPUT PLOT DATA BLOCK
 
 core = core - sysbuf
 IF (core > 0) GO TO 260
 235 ier = 8
 ifle = -core
 237 CALL mesage (ier,ifle,routin)
 
!     CLOSE ANY OPEN FILES AND RETURN
 
 240 CALL CLOSE (xycdb,rewd)
 DO  i = 1,5
   CALL CLOSE (indb(i),rewd)
 END DO
 IF (.NOT.outopn) RETURN
 
!     NO CAMERA PLOTS SO DONT WRITE TRAILER
 
 IF (.NOT. oomcp) GO TO 255
 buf(1) = outfil
 buf(2) = 9999999
 CALL wrttrl (buf(1))
 255 CALL CLOSE  (outfil,rewd)
 GO TO 830
 
 260 CALL OPEN  (*270,outfil,z(core+1),outrwd)
 CALL fname (outfil,NAME(1))
 CALL WRITE (outfil,NAME(1),2,eor)
 outopn = .true.
 GO TO 280
 
!     ERROR,  PLOTS REQUESTED AND OUTFIL PURGED.  DO ALL ELSE.
 
 270 CALL page2 (2)
 WRITE (l,290) uwm,outfil
 plot = .false.
 
 280 IF (buf(3) /= 0) punch = .true.
 TYPE   = buf(4)
 vector = buf(5)
 nsubs  = buf(7)
 knt    = 0
 IF (nsubs > 0) CALL READ (*80,*90,xycdb,subcas(1),nsubs,noeor, flag)
 IF (nsubs > 0) CALL sort (0,0,1,1,subcas(1),nsubs)
 IF (random .AND. TYPE /= 2 .AND. TYPE /= 3) GO TO 380
 IF ((.NOT.random) .AND. (TYPE == 2    .OR.  TYPE == 3) ) GO TO 380
 IF ((.NOT.random) .AND. ipset == pset .AND. vector > 7) GO TO 380
 IF ((.NOT.random) .AND. ipset /= pset .AND. vector <= 7) GO TO 380
 
!     INITIALIZE DATA BLOCK POINTERS FOR THIS VECTOR
 
 FILE = files(vector)
 
!     CHECK FOR RANDOM
 
 IF (random .AND. TYPE == 3) FILE = 2
 IF (random .AND. TYPE == 2) FILE = 1
 ifile = indb(FILE)
 IF (openf(FILE) < 0.0) THEN
   GO TO   360
 ELSE
   GO TO   400
 END IF
 
!     EOR HIT ON IFILE.  SHOULD NOT HAVE HAPPENED
 
 330 ier = 3
 GO TO 355
 
!     EOF HIT ON IFILE.  SHOULD NOT HAVE HAPPENED
 
 350 ier = 2
 355 CALL mesage (ier,ifile,routin)
 openf(FILE) = -1
 
!     FILE IFILE IS NOT SATISFACTORY
 
 360 CALL fname (ifile,buf(1))
 CALL page2 (3)
 WRITE (l,370) uwm,buf(1),buf(2),namev(vector)
 
!     SKIP OVER ANY AND ALL FRAME DATA FOR THIS CARD.
 
 380 CALL READ (*80,*240,xycdb,bcd,1,noeor,flag)
 IF (bcd /= fram) GO TO 800
 390 CALL READ (*80,*90,xycdb,buf(1),3,noeor,flag)
 IF (buf(1) /= -1) GO TO 390
 GO TO 380
 
!     CHECK TO SEE IF THIS FILES SUBCASE IS TO BE OUTPUT
 
 400 CONTINUE
 IF (openf(FILE) < 0.0) THEN
   GO TO   360
 ELSE IF (openf(FILE) == 0.0) THEN
   GO TO   401
 ELSE
   GO TO   402
 END IF
 401 CONTINUE
 CALL fwdrec (*350,ifile)
 CALL READ (*350,*330,ifile,idin(1),20,eor,flag)
 CALL READ (*350,*403,ifile,idin(1),-core,eor,flag)
 GO TO 235
 403 CALL bckrec (ifile)
 CALL bckrec (ifile)
 size  = flag/idin(10)
 ktype = (idin(2)/1000)*1000
 openf(FILE) = 1
 402 CONTINUE
 kasknt = kasknt + 1
 IF (nsubs == 0) GO TO 415
 subc(FILE) = subcas(kasknt)
 GO TO 420
 415 subc(FILE) = 0
 
!     NOW READY TO PROCEED WITH DATA SELECTION
 
 420 CALL READ (*80,*240,xycdb,bcd,1,noeor,flag)
 IF (bcd /= fram) GO TO 800
 
!     READ IN THE ID-COMP-COMP SETS AND SORT ON ID-S.
 
 knt  = 0
 itry = 0
 iat  = 0
 430 CALL READ (*80,*90,xycdb,z(iat+1),3,noeor,flag)
 IF (z(iat+1) == -1) GO TO 440
 iat  = iat + 3
 GO TO 430
 
!     SORT ON ID-S
 
 440 CALL sort (0,0,3,1,z(1),iat)
 450 icore = core - iat
 
!     COMPUTE FINAL REGIONS
 
 nslots= iat/3
 nat   = iat
 IF (z(i3) > 0 .AND. .NOT.random) nslots = nslots + nslots
 554 steps = size
 IF (.NOT.vgp) GO TO 559
 itemp = 0
 nuq   = 0
 DO  i = 1,nat,3
   IF (z(i) == itemp) CYCLE
   nuq   = nuq + 1
   itemp = z(i)
 END DO
 steps = steps*nuq
 
!     SET CORE TO 1
 
 j = iat + 1
 n = j + MIN0(icore,(nslots+1)*steps)
 DO  i = j,n
   z(i) = 1
 END DO
 559 CONTINUE
 IF (steps*(nslots+1) <= icore) GO TO 580
 CALL page2 (4)
 WRITE (l,570) uwm,z(iat-2),z(iat-1),z(iat)
 icrq = steps*(nslots+1) - icore
 WRITE (l,571) icrq
 nslots = nslots - 1
 IF (z(i3) > 0 .AND. .NOT.random) nslots = nslots - 1
 nat = nat - 3
 IF (nslots > 0) GO TO 554
 GO TO 420
 580 ntops = nslots/2
 nbots = ntops
 IF (z(i3) > 0 .AND. .NOT.random) GO TO 590
 ntops = nslots
 nbots = 0
 590 CONTINUE
 center = iat + ntops*steps
 
!     GET CURVE DATA
 
 major = ktype + majid(vector)
 i2    = 0
 ifcrv =-1
 istsv = 0
 idtot = nat/3
 
!     I1 = 1-ST ROW OF NEXT ID
!     I2 = LAST ROW OF NEXT ID
 
 630 i1   = i2 + 1
 nbeg = 3*i1 - 3
 IF (nbeg >= nat) GO TO 780
 idz = nbeg + 1
 id  = z(idz)
 i2  = i1
 640 IF (i2 >= idtot .OR. random) GO TO 650
 IF (z(3*i2+1) /= id) GO TO 650
 i2  = i2 + 1
 GO TO 640
 
!     FIND THIS ID ON IFILE
 
 650 CALL xyfind (*350,*330,*660,majid(1),idz)
 knt = -1
 IF (itry == 0 .AND. subc(FILE) == -1) GO TO 661
 
!     THIS IS THE WAY OUT FOR ALL SUBCASE REQUEST
 
 IF (itry /= 0 .AND. subc(FILE) == -1) GO TO 415
 ktype = (idin(2)/1000)*1000
 IF (ktype == 2000 .OR. ktype == 3000) GO TO 690
 CALL page2 (2)
 WRITE (l,310) uwm
 GO TO 360
 
!     ID NOT FOUND. PRINT MESSAGE AND SHRINK LIST.
 
 
!     SUBCASE REQUEST EITHER SUBCASE NOT FOUND OR POINT NOT FOUND
 
 660 IF (knt == -1) idz = -1
 IF (idz /= -1) GO TO 784
 CALL page2 (3)
 WRITE (l,530) uwm,id,namev(vector),ifile
 WRITE (l,635) subc(FILE)
 knt = 0
 IF (nat/3 <= i2 .AND. i1 == 1) GO TO 784
 GO TO 666
 
!     NSUBS = 0 AND POINT NOT FOUND START FRAME OVER
 
 661 CALL page2 (3)
 WRITE (l,530) uwm,id,namev(vector),ifile
 CALL REWIND (ifile)
 subc(FILE) = 0
 knt = 0
 IF (nat/3 > i2) GO TO 666
 IF (i1    ==  1) GO TO 415
 666 CONTINUE
 i13 = 3*i1 - 3
 i23 = 3*i2 + 1
 IF (i23 >= nat) GO TO 680
 DO  i = i23,nat
   i13 = i13 + 1
   z(i13) = z(i)
 END DO
 680 idtot = idtot - (i2-i1) - 1
 i2  = i1 - 1
 nat = i13
 IF (idz == -1 .AND. i1 /= 1 .AND. .NOT.vgp) GO TO 630
 iat = nat
 GO TO 450
 
!     ID FOUND. READ DATA AND DISTRIBUTE INTO SLOTS.
 
 690 nwds  = idin(10)
 istep = 0
 ifcrv = ifcrv + 1
 IF (vgp) istep = istsv
 700 CALL READ (*350,*630,ifile,buf(1),nwds,noeor,flag)
 istep = istep + 1
 IF (istep > steps) GO TO 700
 itemp = iat + istep
 istsv = istep
 IF (.NOT.vgp) GO TO 709
 IF (ifcrv == 0) GO TO 709
 
!     SORT X AND MOVE Y TO PROPER SLOTS
 
 IF (rbuf(1) >= rz(itemp-1)) GO TO 709
 n = istep - 1
 DO  i = 1,n
   IF (rbuf(1) < rz(iat+i)) GO TO 707
 END DO
 GO TO 709
 707 istep = i
 j  = istep
 itemp = iat + istep
 n  = nslots + 1
 jj = istsv  - 1
 DO  i = 1,n
   item  = iat + (i-1)*steps + j
   temp1 = rz(item)
   z(item) = 1
   DO  ij = j,jj
     item = iat + (i-1)*steps + ij + 1
     temp = rz(item)
     rz(item) = temp1
     temp1 = temp
   END DO
 END DO
 709 rz(itemp) = rbuf(1)
 
!     DISTRIBUTE DATA
 
 DO  i = i1,i2
   place = i*steps + istep
   
!     TOP CURVE
   
   comp = z(3*i-1)
   
!     SET MEAN RESPONSE IF RANDOM
   
   IF (random) z(3*i) = idin(8)
   
!     SET NUMBER OF ZERO CROSSINGS IF RANDOM
   
   IF (random) buf(i+20) = idin(9)
   IF (comp == 1000) GO TO 745
   IF (comp ==    0) CYCLE
   IF (random) comp = 2
   IF (comp <= nwds) GO TO 740
   z(3*i-1) = 0
   CALL page2 (2)
   WRITE (l,730) uwm,comp,id
   GO TO 750
   
   740 itemp = iat + place
   z(itemp) = buf(comp)
   GO TO 750
   745 itemp = iat + place
   z(itemp) = 1
   
!     BOTTOM CURVE IF DOUBLE FRAME
   
   750 IF (random) CYCLE
   comp = z(3*i)
   IF (comp == 1000) GO TO 765
   IF (comp ==    0) CYCLE
   IF (comp <= nwds) GO TO 760
   z(3*i) = 0
   CALL page2 (2)
   WRITE (l,730) uwm,comp,id
   CYCLE
   
   760 itemp = center + place
   z(itemp) = buf(comp)
   CYCLE
   765 itemp = center + place
   z(itemp) = 1
 END DO
 istep = istsv
 GO TO 700
 
!     ALL DATA IS NOW IN SLOTS. INTEGER 1-S REMAIN IN VACANT SLOTS.
 
 780 IF (nsubs /= 0) GO TO 783
 subc(FILE) = idin(4)
 783 CONTINUE
 CALL xydump (outfil,TYPE)
 knt = 1
 IF (nsubs /= 0) GO TO 784
 subc(FILE) = 0
 itry = itry + 1
 GO TO 450
 784 CONTINUE
 IF (kasknt < nsubs) GO TO 785
 kasknt = 0
 GO TO 402
 785 kasknt = kasknt + 1
 subc(FILE) = subcas(kasknt)
 DO  i = 1,5
   vecid(i) = 0
 END DO
 GO TO 450
 
!     INITIALIZE PARAMETERS
 
 800 plot  = .false.
 punch = .false.
 PRINT = .false.
 paplot= .false.
 DO  i = 1,5
   vecid(i) = 0
 END DO
 GO TO 100
 
!     VALUE DUMP
 
 820 CONTINUE
 GO TO 140
 
!     INTERACTIVE STOP INITIATED HERE.
 
 1500 nogo = 1
 RETURN
 
!     INSUFFICIENT CORE
 
 825 CALL mesage (8,-core,routin)
 
!     CALL THE PRINTER-PLOTTER IF ANY REQUESTS FOR PRINTER-PLOTTER
 
 830 IF (oompp) CALL xyprpl
 
!     RESTORE OUTPUT HEADING AND RETURN
 
 835 DO  i = 1,96
   ihead(i) = headsv(i)
 END DO
 RETURN
 
 
 120 FORMAT (a25,' 975, XYTRAN DOES NOT RECOGNIZE ',a4, ' AND IS IGNORING')
 290 FORMAT (a25,' 976, OUTPUT DATA BLOCK',i4,' IS PURGED.',  &
     '  XYTRAN WILL PROCESS ALL REQUESTS OTHER THAN PLOT')
 310 FORMAT (a25,' 977, FOLLOWING NAMED DATA-BLOCK IS NOT IN SORT-II',  &
     ' FORMAT')
 370 FORMAT (a25,' 978', /5X,'XYTRAN MODULE FINDS DATA-BLOCK(',2A4,  &
     ') PURGED, NULL, OR INADEQUATE, AND IS IGNORING XY-OUTPUT',  &
     ' REQUEST FOR -',a4,'- CURVES')
 530 FORMAT (a25,' 979, AN XY-OUTPUT REQUEST FOR POINT OR ELEMENT ID',  &
     i10, /5X,1H-,a4,'- CURVE IS BEING PASSED OVER.  THE ID ',  &
     'COULD NOT BE FOUND IN DATA BLOCK',i10)
 570 FORMAT (a25,' 980, INSUFFICIENT CORE TO HANDLE ALL DATA FOR ALL ',  &
     'CURVES OF THIS FRAME', /5X,' ID =',i10,2(' COMPONENT =',  &
     i4,5X),' DELETED FROM OUTPUT')
 571 FORMAT (5X,'ADDITIONAL CORE NEEDED =',i9,' WORDS.')
 635 FORMAT (5X,'SUBCASE',i10 )
 730 FORMAT (a25,' 981, COMPONENT =',i10,' FOR ID =',i10,  &
     ' IS TOO LARGE. THIS COMPONENTS CURVE NOT OUTPUT')
 
 900 FORMAT ('  ENTER XYPLOT DEFINITION OR GO TO PLOT OR STOP TO EXIT')
 902 FORMAT ('  BAD CARD TRY AGAIN')
 910 FORMAT (20A4)
END SUBROUTINE xytran
