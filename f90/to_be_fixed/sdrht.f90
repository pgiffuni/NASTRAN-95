SUBROUTINE sdrht
     
!     SPECIAL FLUX-DATA-RECOVERY MODULE FOR HBDY ELEMENTS IN HEAT
!     TRANSFER.
 
!     DMAP CALLING SEQUENCE.
 
!     SDRHT SIL,USET,UGV,OEF1,SLT,EST,DIT,QGE,DLT,/OEF1X/V,N,TABS $
 
 LOGICAL :: havids,cardin,lhbdy,transt,found,mch521
 INTEGER :: tablst(13),z,buf(50),sysbuf,rd,rdrew,wrt,wrtrew,  &
     clsrew,cls,subr(2),NAME(2),eor,idpos(3),outpt,  &
     sltyps,ldword(16),ug,oef1,slt,est,dit,qge,dlt,  &
     oef1x,buf1,buf2,buf3,core,pass,hbdytp,mcbugv(7),  &
     FILE,eltype,estwds,ecpt(100),sltat,sltrec,mcb(7), eol,gsize
 REAL :: rz(1),rbuf(50),grids(6)
 CHARACTER (LEN=27) :: swm
 CHARACTER (LEN=25) :: sfm
 CHARACTER (LEN=29) :: uim
 CHARACTER (LEN=25) :: uwm
 CHARACTER (LEN=23) :: ufm
 COMMON /xmssg / ufm,uwm,uim,sfm,swm
 COMMON /condas/ consts(5)
 COMMON /system/ ksystm(65)
 COMMON /names / rd,rdrew,wrt,wrtrew,clsrew,cls
 COMMON /zntpkx/ ai(4),irow,eol
 COMMON /zzzzzz/ z(1)
 COMMON /machin/ mach
 COMMON /BLANK / tabs
 EQUIVALENCE     (ksystm(1),sysbuf),(ksystm(2),outpt),  &
     (consts(2),twopi ),(consts(3),raddeg), (z(1),rz(1)), (buf(1),rbuf(1))
 DATA    tablst / 4, 1105,11,1, 1205,12,2, 1305,13,3, 1405,14,4  /
 DATA    lentry / 14   /, eor,noeor  / 1, 0 /, subr/ 4HSDRH,2HT  /
 DATA    idpos  / 2,1,5/, hbdytp/ 52 /,     sltyps / 16 /
 DATA    ldword / 6,6,4,4,6,6,2,5,5,6,6,7,2,2,5,5  /
 DATA    ug,oef1, slt,est,dit,qge,dlt,oef1x /  &
     103,104, 105,106,107,108,109, 201  /
 DATA    grids  / 1.0, 2.0, 2.0, 3.0, 4.0, 2.0     /
 
!     SET UP CORE AND BUFFERS
 
 buf1   = korsz(z) - sysbuf - 2
 buf2   = buf1 - sysbuf - 2
 buf3   = buf2 - sysbuf - 2
 core   = buf3 - 1
 idrec  = 1
 ieltab = 1
 pass   = 0
 havids = .false.
 cardin = .false.
 mch521 = mach == 5 .OR. mach == 21
 
!     OPEN INPUT FORCES -OEF1- AND OUTPUT FORCES -OEF1X-.
 
 CALL OPEN  (*1190,oef1,z(buf1),rdrew)
 CALL OPEN  (*1200,oef1x,z(buf2),wrtrew)
 CALL fname (oef1x,NAME)
 CALL WRITE (oef1x,NAME,2,eor)
 CALL fwdrec (*1220,oef1)
 
!     COPY RECORD PAIRS OF DATA FROM OEF1 TO OEF1X UNTIL HBDY DATA IS
!     DETECTED.
 
 10 lcore = core - idrec
 IF (lcore < 300) CALL mesage (-8,0,subr)
 FILE = oef1
 CALL READ (*1180,*30,oef1,z(idrec),lcore,noeor,iflag)
 WRITE  (outpt,20) swm
 20 FORMAT (a27,' 3063, INPUT FORCES DATA BLOCK HAS INCORRECT DATA.')
 GO TO 1220
 
!     MODIFY ID-RECORD IF THIS IS FOR HBDY ELEMENTS.
 
 30 IF (z(idrec+2) /= hbdytp) GO TO 40
 lhbdy = .true.
 GO TO 50
 40 lhbdy = .false.
 GO TO 90
 
!     SET CONSTANTS FROM OEF1 ID RECORD.
 
 50 IF (z(idrec)/10 == 6) GO TO 70
 loadid = z(idrec+7)
 transt = .false.
 GO TO 80
 70 loadid = z(idrec+7)
 time   = rz(idrec+4)
 transt = .true.
 80 lwords = z(idrec+9)
 z(idrec+9) = 5
!/////
!     CALL BUG (4HTRAN ,90,TRANST,1)
!/////
 90 CALL WRITE (oef1x,z(idrec),iflag,eor)
 IF (lhbdy) GO TO 120
 
!     NOT AN HBDY ELEMENT TYPE THUS COPY DATA ACROSS.
 
 100 CALL READ  (*1260,*110,oef1,z(idrec),lcore,noeor,iflag)
 CALL WRITE (oef1x,z(idrec),lcore,noeor)
 GO TO 100
 110 CALL WRITE (oef1x,z(idrec),iflag,eor)
 GO TO 10
 
!     HBDY ELEMENT DATA ENCOUNTERED.
 
 120 pass = pass + 1
 
!     ON FIRST PASS ELEMENT-DATA-TABLE IS FORMED.
 
!     EACH ENTRY WILL CONTAIN,
 
!      1) ELEMENT-ID (HBDY).
!      2) FLUX-RADIATION TERM FOR THIS ELEMENT.
!      3) FLUX-X FROM OEF1 DATA.
!      4) APPLIED LOAD (USING SLT DATA).
!      5) HBDY ELEMENT TYPE (1 TO 6).
!      6) HBDY AREA FACTOR.
!      7) ALPHA VALUE.
!      8) V1(1) *
!      9) V1(2)  *  VECTOR-V1
!     10) V1(3) *
!     11) V2(1) *
!     12) V2(2)  *  VECTOR-V2
!     13) V2(3) *
!     14) OUTPUT ID*10 + DEVICE CODE (FROM OEF1).
 
!     ON PASSES OTHER THAN THE FIRST ONLY THE FLUX-X VALUE IS EXTRACTED
!     FROM OEF1-HBDY-ENTRIES.
 
 jeltab = ieltab - 1
 
!     INPUT = ID*10+CODE, NAME1,NAME2,GRD-X,GRD-Y,GRD-Z,FLUX-X,FLUX-Y,
!             FLUX-Z   (TOTAL OF 9 WORDS)
 
 130 CALL READ (*1260,*150,oef1,buf,lwords,noeor,iflag)
 IF (pass > 1) GO TO 140
 z(jeltab+1) = buf(1)/10
 
!     STORE (FLUX-X) AND (OUTPUT ID*10 + DEVICE CODE).
 
 140 z(jeltab+ 3) = buf(7)
 rz(jeltab+2) = 0.0
 z(jeltab+14) = buf(1)
 jeltab = jeltab + lentry
 IF (jeltab > core) CALL mesage (-8,0,subr)
 GO TO 130
 
!     END OF DATA.
 
 150 IF (pass > 1) GO TO 160
 neltab = jeltab
 NUMBER = (neltab - ieltab + 1)/lentry
 
!     OPEN UG FILE FOR INPUT OF UG VECTORS.
 
 CALL OPEN (*1220,ug,z(buf3),rdrew)
 CALL fwdrec (*1220,ug)
 mcbugv(1) = ug
 CALL rdtrl (mcbugv)
 gsize = mcbugv(3)
 GO TO 180
 160 IF (jeltab == neltab) GO TO 180
 WRITE  (outpt,170) swm,jeltab,neltab
 170 FORMAT (a27,' 3064, INCONSISTANT HBDY DATA RECORDS.  ',2I20)
 GO TO 1220
 
!     ALL DATA FROM OEF1 IS AT HAND NOW.
!     FILES ARE CLOSED WITHOUT REWIND.
 
 180 CALL CLOSE (oef1,cls)
 CALL CLOSE (oef1x,cls)
 
!     ALLOCATE UG VECTOR SPACE ON FIRST PASS.
 
 IF (pass /= 1) GO TO 190
 iug  = neltab + 1
 nug  = neltab + gsize
 iugz = iug - 1
 IF (nug + 5 > core) CALL mesage (-8,0,subr)
 
!     BRING NEXT DISPLACEMENT VECTOR INTO CORE.
 
 190 DO  i = iug,nug
   rz(i) = tabs
 END DO
 
 CALL intpk (*220,ug,0,1,0)
 210 CALL zntpki
 kk = iugz + irow
 rz(kk) = rz(kk) + ai(1)
 IF (eol > 0.0) THEN
   GO TO   220
 ELSE
   GO TO   210
 END IF
 
!     RAISE VECTOR RESULT TO 4TH POWER
 
 220 DO  i = iug,nug
   rz(i) = rz(i)**4
 END DO
 
!     IF TRANSIENT PROBLEM SKIP ACCELERATION AND VELOCITY VECTORS
 
 IF (.NOT. transt) GO TO 250
 DO  i = 1,2
   CALL fwdrec (*250,ug)
 END DO
 250 CONTINUE
!/////
!     CALL BUG (4HUG4  ,200,Z(IUG),NUG-IUG+1)
!/////
 
!     IF NONLINEAR PROBLEM, COMPUTE FLUX RADIATION TERMS.
 
!                                            T           4
!     (FLUX-RADIATION         ) = (Q           )(U +TABS)
!                    EL-SUBSET      G,EL-SUBSET   G
 
 FILE = qge
 CALL OPEN (*261,qge,z(buf1),rdrew)
 IF (pass == 1) GO TO 260
 CALL fwdrec (*1260,qge)
 GO TO 280
 261 iqgid = nug + 1
 nqgid = nug
 idrec = nqgid
 GO TO 311
 260 CALL READ (*1260,*1270,qge,buf,-2,noeor,iflag)
 
!     ON FIRST PASS PICK UP ELEMENT ID LIST.
 
 iqgid = nug + 1
 CALL READ (*1260,*270,qge,z(iqgid),core-iqgid,noeor,iflag)
 CALL mesage (-7,0,subr)
 270 nqgid = nug + iflag
 idrec = nqgid + 1
!/////
!     CALL BUG (4HQGID,410,Z(IQGID),NQGID-IQGID+1)
!/////
 
!     EACH FLUX-RADIATION TERM IN THE ELEMENT TABLE IS CREATED BY
!     FORMING THE DOT-PRODUCT OF THE COLUMN OF -QGE- HAVING THE
!     SAME ELEMENT-ID WITH THE -UG- VECTOR IN CORE.
 
 280 DO  i = iqgid,nqgid,1
   CALL intpk (*310,qge,0,1,0)
   
!     FIND OUT IF ID OF THIS VECTOR IS IN ELEMENT TABLE.
   
   kid = z(i)
   CALL bisloc (*300,kid,z(ieltab),lentry,NUMBER,jpoint)
   jword = ieltab + jpoint
   z(jword) = 0
   
!     FORM DOT PRODUCT
   
   290 CALL zntpki
   k = iugz + irow
   rz(jword) = rz(jword) - ai(1)*rz(k)
   IF (eol > 0.0) THEN
     GO TO   310
   ELSE
     GO TO   290
   END IF
   
!     ID OF THIS COLUMN NOT IN ELEMENT TABLE
   
   300 CALL zntpki
   IF (eol > 0.0) THEN
     GO TO   310
   ELSE
     GO TO   300
   END IF
 END DO
!/////
!     CALL BUG (4HELTB ,440,Z(IELTAB),NELTAB-IELTAB+1)
!/////
 CALL CLOSE (qge,clsrew)
 
!     ON FIRST PASS, EST IS PASSED AND HBDY ELEMENTS CALLED.
 
 311 CONTINUE
 IF (pass > 1) GO TO 380
 FILE = est
 CALL gopen (est,z(buf1),rdrew)
 next = ieltab
 
!     READ THE ELEMENT TYPE
 
 320 CALL READ (*350,*1270,est,eltype,1,noeor,iflag)
 IF (eltype == hbdytp) GO TO 330
 CALL fwdrec (*1260,est)
 GO TO 320
 
!     HBDY ELEMENT-SUMMARY-TABLE DATA FOUND.
 
 330 estwds = 53
 340 CALL READ (*1260,*1270,est,ecpt,estwds,noeor,iflag)
 
!     CHECK TO SEE IF THIS ELEMENT IS IN OUTPUT SET.
 
 IF (z(next) - ecpt(1) < 0.0) THEN
   GO TO   350
 ELSE IF (z(next) - ecpt(1) == 0.0) THEN
   GO TO   370
 ELSE
   GO TO   340
 END IF
 350 WRITE  (outpt,360) swm,z(next)
 360 FORMAT (a27,' 3065, THERE IS NO EST DATA FOR HBDY ELEMENT ID =', i10)
 GO TO 1220
 
!     THIS ELEMENT IS IN TABLE.
 
 370 CALL hbdy (ecpt,ecpt,2,rbuf,buf)
 
!     PLANT HBDY OUTPUTS INTO TABLE.
 
 z(next+ 4) = ecpt(2)
 z(next+ 5) = buf(2)
 z(next+ 6) = ecpt(17)
 z(next+ 7) = buf(11)
 z(next+ 8) = buf(12)
 z(next+ 9) = buf(13)
 z(next+10) = buf(14)
 z(next+11) = buf(15)
 z(next+12) = buf(16)
 next = next + lentry
 IF (next < neltab) GO TO 340
 CALL CLOSE (est,clsrew)
 
!     LOAD SET PROCESSING IF LOAD-SET-ID IS NON-ZERO.
 
 380  IF (loadid > 0) THEN
   GO TO   390
 ELSE
   GO TO   944
 END IF
 
!     OPEN SLT FOR LOAD DATA.
 
 
 390 FILE = slt
 CALL OPEN (*1250,slt,z(buf1),rdrew)
 IF (havids) CALL fwdrec (*1260,slt)
 IF (havids) GO TO 810
 havids = .true.
 ildid  = nqgid + 1
 nldid  = ildid - 1
 
!     IDS OF LOAD SETS NOT IN CORE THUS BRING IN IDS FROM HEADER RECORD.
 
 imast  = ildid
 nsets  = 0
 nldset = 3
 CALL READ (*1260,*1270,slt,buf,-2,noeor,iflag)
 400 IF (nldid+5 > core) CALL mesage (-8,0,subr)
 CALL READ (*1260,*410,slt,z(nldid+1),1,noeor,iflag)
 nsets = nsets + 1
 z(nldid +2) = 1
 z(nldid +3) = z(nldid+1)
 rz(nldid+4) = 1.0
 z(nldid +5) = nsets
 nldid = nldid + 5
 GO TO 400
 
!     IF TRANSIENT PROBLEM THEN DLT OPERATIONS BEGIN
 
 410 IF (.NOT.transt) GO TO 800
 FILE = dlt
 CALL OPEN (*1250,dlt,z(buf2),rdrew)
 CALL READ (*1260,*1270,dlt,buf,3,noeor,iflag)
 m = buf(3)
 found = .false.
 IF (m <= 0) GO TO 430
 DO  i = 1,m
   CALL READ (*1260,*1270,dlt,buf,1,noeor,iflag)
   IF (buf(1) == loadid) found = .true.
 END DO
 
!     NOW READ RLOAD1, RLOAD2, TLOAD1, AND TLOAD2 IDS.
 
 430 irtids = nldid + 1
 CALL READ (*1260,*440,dlt,z(irtids),core-irtids,noeor,iflag)
 CALL mesage (-7,0,subr)
 GO TO 1220
 440 nrtids = nldid + iflag
 
!     IF LOADID WAS FOUND AMONG THE DLOAD IDS, SEARCH IS NOW MADE IN
!     RECORD 1 OF THE DLT FOR THAT ID, AND ITS SUB-IDS.
 
 jj1   = ildid
 jj2   = nldid
 ildid = nrtids + 1
 nldid = nrtids + 2
 z(ildid  ) = loadid
 z(ildid+1) = 0
 IF (.NOT. found) GO TO 520
 
!     READ A MASTER DLOAD SET-ID.
 
 450 CALL READ (*1260,*1270,dlt,buf,2,noeor,iflag)
 IF (buf(1) == loadid) GO TO 470
 
!     SKIP SUB-ID DATA OF THIS MASTER
 
 460 CALL READ (*1260,*1270,dlt,buf,2,noeor,iflag)
 IF (buf(2) < 0.0) THEN
   GO TO   450
 ELSE
   GO TO   460
 END IF
 
!     MASTER-ID FOUND.  BUILD LOAD-SET-ID TABLE.
 
 470 factor = rbuf(2)
 472 IF (nldid+11 <= core) GO TO 480
 CALL mesage (8,0,subr)
 GO TO 1220
 480 CALL READ (*1260,*1270,dlt,buf,2,noeor,iflag)
 IF (buf(2) > 0.0) THEN
   GO TO   490
 ELSE
   GO TO   540
 END IF
 490 z(ildid +1) = z(ildid+1) + 1
 z(nldid +1) = buf(2)
 rz(nldid+2) = 0.0
 z(nldid +3) = 0
 rz(nldid+4) = rbuf(1)*factor
 z(nldid +5) = 0
 nldid = nldid + 11
 GO TO 472
 
!     LOADID NOT AMONG DLOADS FOR THIS TRANSIENT PROBLEM
 
 520 z(ildid+1) = 1
 IF (nldid+13 <= core) GO TO 530
 CALL mesage (8,0,subr)
 GO TO 1220
 530 z(nldid +1) = z(ildid)
 rz(nldid+2) = 0.0
 z(nldid +3) = 0
 rz(nldid+4) = 1.0
 z(nldid +5) = 0
 nldid = nldid + 11
 
!     IF THERE ARE ANY DLOAD CARDS AT ALL THEN BALANCE OF (OR ALL OF)
!     RECORD 1 IS NOW SKIPPED.
 
 540 IF (m > 0) CALL fwdrec (*1260,dlt)
 
!     NOW PICKING UP DATA NEEDED OF DYNAMIC LOAD SET RECORDS.
 
 k1 = ildid + 2
 k2 = nldid
 DO  i = irtids,nrtids
   
!     READ THE LOAD TYPE
   
   CALL READ (*1260,*1270,dlt,buf,2,noeor,iflag)
   IF (buf(1) /= 3 .AND. buf(1) /= 4) GO TO 570
   
!     CHECK AND SEE IF THIS TLOAD ID IS AMONG THE SUB-IDS
   
   DO  j = k1,k2,11
     IF (z(j) == z(i)) GO TO 560
   END DO
   GO TO 570
   
!     YES THIS RECORD IS NEEDED.  THUS PUT ITS DATA IN TABLE.
   
   560 z(j+4) = buf(1)
   
!     SLT ID INTO TABLE
   
   z(j) = -buf(2)
   
!     SET SLT RECORD NUMBER
   
   k = 0
   DO  l = jj1,jj2,5
     k = k + 1
     IF (z(l) == buf(2)) GO TO 566
   END DO
   k = 0
   566 z(j+2) = k
   CALL READ (*1260,*1270,dlt,z(j+5),6,eor,iflag)
   IF (buf(1) == 3) z(j+6) = 0
   CYCLE
   570 CALL fwdrec (*1260,dlt)
 END DO
 
!     CHECK IS NOW MADE TO INSURE ALL SUB-IDS RECEIVED DLT DATA.
 
 
!     SET SLT IDS POSITIVE
 
 DO  i = k1,k2,11
   z(i) = IABS(z(i))
 END DO
 DO  i = k1,k2,11
   IF (z(i+4) > 0.0) THEN
     GO TO   610
   END IF
   
!     ERROR
   
   590 WRITE  (outpt,600) uwm,z(i)
   600 FORMAT (a25,' 3066, THERE IS NO TLOAD1 OR TLOAD2 DATA FOR LOAD-',  &
       'ID =',i9)
 END DO
 
 CALL CLOSE (dlt,clsrew)
 nldset = 11
 
!     SORT SUB-ID TABLE ON SLT RECORD NUMBERS.
 
 CALL sort (0,0,11,3,z(k1),k2-k1+1)
!/////
!     CALL BUG (4HTABL,640,Z(K1),K2-K1+1)
 
!     CONSTRUCTION OF TABLE-ID LIST.
 
 itabid = nldid + 1
 ntabid = nldid + 1
 
!     FIRST GET TABLE ID-S PRESENT IN THE SUB-ID TABLE.
 
 loop680:  DO  i = k1,k2,11
   
!     CHECK FOR OTHER THAN TLOAD1 TYPE CARD
   
   IF (z(i+4) /= 3) CYCLE loop680
   
!     CHECK FOR ID IN TABLE.
   
   IF (ntabid <= itabid) GO TO 660
   DO  j = itabid,ntabid
     IF (z(i+5) == z(j)) CYCLE loop680
   END DO
   660 ntabid = ntabid + 1
   IF (ntabid > core) CALL mesage (-8,0,subr)
   z(ntabid) = z(i+5)
 END DO loop680
 
!     NOW PASS SLT AND GET ANY TABLE IDS PRESENT IN QVECT PORTION OF
!     RECORDS WE WILL BE USING.  (SLT IS CURRENTLY POSITIONED AT FIRST
!     RECORD.)
 
 sltat = 1
 FILE  = slt
 DO  i = k1,k2,11
   ngo = z(i+2) - sltat
   IF (ngo < 0) THEN
     GO TO   780
   ELSE IF (ngo == 0) THEN
     GO TO   710
   END IF
   690 DO  j = 1,ngo
     CALL fwdrec (*1260,slt)
   END DO
   sltat = sltat + ngo
   
!     LOOK FOR QVECT CARDS.
   
   710 CALL READ (*1260,*770,slt,buf,2,noeor,iflag)
   itype  = buf(1)
   ncards = buf(2)
   nwords = ldword(itype)
   IF (itype /= 16) GO TO 760
   
!     QVECT CARDS FOUND
   
   IF (ncards <= 0) GO TO 710
   DO  j = 1,ncards
     CALL READ (*1260,*1270,slt,buf,nwords,noeor,iflag)
     loop740:  DO  k = 2,4
       l = numtyp(buf(k))
       IF (mch521 .AND. buf(k) > 16000 .AND. buf(k) <= 99999999) l= 1
       IF (buf(k) <= 0 .OR. l /= 1) CYCLE loop740
       
!     TABLE ID FOUND.  ADD TO LIST IF NOT YET IN.
       
       IF (ntabid <= itabid) GO TO 730
       DO  l = itabid,ntabid
         IF (buf(k) == z(l)) CYCLE loop740
       END DO
       730 ntabid = ntabid + 1
       IF (ntabid > core) CALL mesage (-8,0,subr)
       z(ntabid) = buf(k)
     END DO loop740
   END DO
   GO TO 710
   760 IF (itype > 16) GO TO 1170
   nwdcrd = -nwords*ncards
   CALL READ (*1260,*1270,slt,buf,nwdcrd,noeor,iflag)
   GO TO 710
   770 sltat = sltat + 1
 END DO
 numtab   = ntabid - itabid
 z(itabid)= numtab
 numtab   = z(itabid)
 
!     TABLE-ID LIST COMPLETE. NOW SORT IT AND PRIME TAB ROUTINE.
 
 CALL REWIND (slt)
 CALL fwdrec (*1260,slt)
!/////
!     CALL BUG (4HTBID,555,Z(ITABID),NTABID-ITABID+1)
!/////
 idit = ntabid + 1
 ndit = idit
 lz   = core - idit
 IF (lz > 10) GO TO 790
 CALL mesage (8,0,subr)
 GO TO 1220
 790 CONTINUE
 IF (numtab == 0) GO TO 792
 CALL sort (0,0,1,1,z(itabid+1),numtab)
 CALL pretab(dit,z(idit),z(idit),z(buf2),lz,lused,z(itabid),tablst)
 ndit = idit + lused
!/////
!     CALL BUG (4HDITS,557,Z(IDIT),NDIT-IDIT+1)
!/////
 792 CONTINUE
 idrec = ndit + 1
 GO TO 810
 
!     DETERMINE IF -LOADID- IS IN LIST OF LOAD SET IDS.
 
 800 idrec = nldid + 1
!/////
!     CALL BUG (4HLD1   ,360,Z(ILDID),NLDID-ILDID+1)
!/////
 810 CONTINUE
 nmast = nldid
 j = ildid
 820 IF (j <= nldid) GO TO 910
 
!     LOAD SET ID LIST EXHAUSTED.
!     BRING IN ANY LOAD CARDS IF NOT YET IN.
 
 IF (cardin) GO TO 920
 
!     THE LOAD-SET-ID TABLE HAS THE FOLLOWING FORMAT.
 
!     MASTER ID                                ******       Z(ILDID)
!     NUMBER OF SUB-IDS FOR THIS MASTER              *
!     SUB-ID             **  3-WORDS  ***             *
!     SCALE FACTOR = F(T)  * ONLY IF     *             *
!     SLT RECORD NUMBER  **  STATICS      * 11 WORDS   * REPEATS FOR
!     CONSTANT SCALE FACTOR               * REPEATS    * EACH MASTER
!     TYPE OF TLOAD =(3 OR 4)             * FOR EACH   * ID PRESENT
!     TYPE3 TABLE ID (OR) TYPE4 T1        * SUB-ID     *
!           0                   T2        * OF THIS    *
!           0                   OMEGA     * MASTER     *
!           0                   PHI       * ID         *
!           0                   N        *             *
!           0                   ALPHA ***              *
!              .                                       *
!              .                                      *
!              .                                     *
!                                              ******
!             ...                    ...
!             ...                    ...
!             ...                    ...                    Z(NLDID)
 
 
 cardin = .true.
 
!     FORWARD SLT TO LOAD CARD RECORD.
 
 IF (nsets > 0) THEN
   GO TO   830
 ELSE
   GO TO   850
 END IF
 830 DO  i = 1,nsets
   CALL fwdrec (*1250,slt)
 END DO
 
!     READ AND ENTER MASTER ID INTO TABLE
 
 850 IF (nldid+2 > core) CALL mesage (-8,0,subr)
 CALL READ (*1260,*900,slt,z(nldid+1),2,noeor,iflag)
 scale  = rz(nldid+2)
 nldid  = nldid + 2
 jcount = nldid
 z(jcount) = 0
 
!     READ THE (SID, SCALE-FACTOR)  PAIRS FOR THIS ID.
 
 860 IF (nldid+3 > core) CALL mesage (-8,0,subr)
 CALL READ (*1260,*1270,slt,z(nldid+1),2,noeor,iflag)
 IF (z(nldid+1) == -1) GO TO 890
 
!     MULTIPLY SUBID SCALE FACTOR BY MASTER SCALE FACTOR.
 
 rz(nldid+2) = rz(nldid+2)*scale
 
!     DETERMIND SLT RECORD NUMBER OF THIS SUB ID.
 
 krec = 1
 DO  i = imast,nmast,nldset
   IF (z(nldid+1) == z(i)) GO TO 880
   krec = krec + 1
 END DO
 krec = 0
 880 z(nldid+3) = krec
 nldid = nldid + 3
 z(jcount) = z(jcount) + 1
 GO TO 860
 
!     SORT ALL SUB-ID 3 WORD GROUPS ON SLT RECORD NUMBER.
 
 890 CALL sort (0,0,3,3,z(jcount+1),nldid-jcount)
 GO TO 850
 
!     REPOSITION SLT TO BEGINNING OF FIRST SLT RECORD
 
 900 idrec = nldid + 1
!/////
!     CALL BUG (4HLDID,460,Z(ILDID),NLDID - ILDID+1)
!/////
 CALL REWIND (slt)
 CALL fwdrec (*1260,slt)
 GO TO 820
 
!     CONTINUE SEARCH FOR LOADID
 
 910 IF (loadid == z(j)) GO TO 940
 
!     POSITION -J- TO NEXT LOAD-SET-ID IN TABLE
 
 j = j + nldset*z(j+1) + 2
 GO TO 820
 
!     -LOADID- NOT FOUND ANYWHERE.
 
 920 WRITE  (outpt,930) uwm,loadid
 930 FORMAT (a25,' 3067, LOAD SET ID =',i9,' IS NOT PRESENT.')
 GO TO 1220
 
!     MATCH ON MASTER ID HAS BEEN FOUND.
 
 940 nloads = z(j+1)
 iload  = j + 2
 nload  = iload + nldset*nloads - 1
 
!     PROCESS ALL THE LOAD RECORDS FOR THIS MASTER-ID
 
 sltat = 1
!/////
!     CALL BUG (4HLOAD ,500,Z(ILOAD),NLOAD-ILOAD+1)
!/////
 
!     INITIALIZE APPLIED LOAD TO 0.0 FOR ALL ELEMENTS IN TABLE
 
 944  CONTINUE
 DO  i = ieltab,neltab,lentry
   rz(i+3) = 0.0
 END DO
 IF (loadid <= 0) GO TO 1140
 DO  i = iload,nload,nldset
   factor = rz(i+1)
   IF (.NOT.transt) GO TO 960
   
!     FACTOR HAS TO BE FOUND AS F(TIME)
   
   IF (z(i+4) == 4) GO TO 950
   CALL tab (z(i+5),time,yvalue)
   factor = rz(i+3)*yvalue
   GO TO 960
   950 tt = time - rz(i+5)
   IF (tt == 0.0) GO TO 951
   IF (tt <= 0.0 .OR. time >= rz(i+6)) GO TO 955
   factor = rz(i+3)*EXP(rz(i+10)*tt)*(tt**rz(i+9))*COS(twopi*  &
       rz(i+7)*tt + rz(i+8)/raddeg)
   GO TO 960
   951 IF (rz(i+9) /= 0.0) GO TO 955
   factor = COS(twopi*rz(i+7))
   GO TO 960
   955 factor = 0.0
   960 sltrec = z(i+2)
   IF (sltrec <= 0 .OR. factor == 0.0) CYCLE
   
!     POSITION SLT TO RECORD DESIRED.
   
   980 ngo = sltrec - sltat
   IF (ngo < 0) THEN
     GO TO   990
   ELSE IF (ngo == 0) THEN
     GO TO  1020
   ELSE
     GO TO  1000
   END IF
   
!     NEED TO BACK UP ON SLT.
   
   990 CALL bckrec (slt)
   sltat = sltat - 1
   GO TO 980
   
!     NEED TO GO FORWARD ON SLT
   
   1000 DO  j = 1,ngo
     CALL fwdrec (*1260,slt)
   END DO
   sltat = sltat + ngo
   
!     SLT IS NOW POSITIONED TO LOAD RECORD DESIRED.
   
   
!     GENERATE LOADS FOR THOSE ELEMENTS IN THE TABLE USING ONLY QBDY1,
!     QBDY2, AND QVECT CARDS.
   
   1020 CALL READ (*1260,*1130,slt,buf,2,noeor,iflag)
   itype = buf(1)
   IF (itype <= sltyps) GO TO 1040
   WRITE  (outpt,1030) swm,itype
   1030 FORMAT (a27,' 3068, UNRECOGNIZED CARD TYPE =',i9,  &
       ' FOUND IN -SLT- DATA BLOCK.')
   GO TO 1220
   1040 ncards = buf(2)
   IF (ncards > 0) THEN
     GO TO  1050
   ELSE
     GO TO  1020
   END IF
   1050 nwords = ldword(itype)
   IF (itype >= 14 .AND. itype <= 16) GO TO 1060
   nwdcrd = -nwords*ncards
   CALL READ (*1260,*1270,slt,buf,nwdcrd,noeor,iflag)
   GO TO 1020
   1060 itype = itype - 13
   jid = idpos(itype)
   DO  k = 1,ncards
     
!     READ A QBDY1, QBDY2, OR QVECT ENTRY.
     
     CALL READ (*1260,*1270,slt,buf,nwords,noeor,iflag)
     
!     CHECK FOR ID IN THE TABLE (OTHERWISE SKIP).
     
     CALL bisloc (*1120,buf(jid),z(ieltab),lentry,NUMBER,jpoint)
     kk = ieltab + jpoint
     
!     THIS ELEMENT IS IN TABLE, THUS COMPUTE AND SUM IN THE LOAD.
     
     SELECT CASE ( itype )
       CASE (    1)
         GO TO 1070
       CASE (    2)
         GO TO 1080
       CASE (    3)
         GO TO 1090
     END SELECT
     1070 rz(kk+2) = rz(kk+2) + rz(kk+4)*rbuf(1)*factor
     CYCLE
     
     1080 ktype    = z(kk+3)
     rz(kk+2) = rz(kk+2) + factor*rz(kk+4)*(rbuf(2)+rbuf(3) + rbuf(4) +  &
         rbuf(5))/grids(ktype)
     CYCLE
     
     
!     CALL TAB IF E1,E2,E3 OF QVECT DATA ARE TABLE ID-S IMPLYING
!     TIME DEPENDENCE
     
     1090 IF (.NOT.transt) GO TO 1099
     DO  kkk = 2,4
       l = numtyp(buf(kkk))
       IF (mch521 .AND. buf(kkk) > 16000 .AND. buf(kkk) <= 99999999) l = 1
       IF (buf(kkk) <= 0 .OR. l /= 1) CYCLE
       CALL tab (buf(kkk),time,yvalue)
       rbuf(kkk) = yvalue
     END DO
     1099 ktype = z(kk+3)
     c = rbuf(2)*rz(kk+6) + rbuf(3)*rz(kk+7) + rbuf(4)*rz(kk+8)
     IF (ktype == 6) GO TO 1110
     IF (c > 0.0) THEN
       GO TO  1120
     END IF
     1100 rz(kk+2) = rz(kk+2) - c*rz(kk+4)*rz(kk+5)*rbuf(1)*factor
     CYCLE
     1110 rz(kk+2) = rz(kk+2) + factor*rz(kk+4)*rbuf(1)*rz(kk+5)*SQRT(c*c +  &
         (rbuf(2)*rz(kk+9) + rbuf(3)*rz(kk+10) + rbuf(4)*rz(kk+11))**2)
     
   END DO
   
 END DO
 CALL CLOSE (slt,clsrew)
!/////
!     CALL BUG (4HTELT,670,Z(IELTAB),NELTAB-IELTAB+1)
!/////
 
!     ELEMENT TABLE IS NOW COMPLETE FOR OUTPUT.
 
 1140 FILE = oef1x
 CALL OPEN (*1250,oef1x,z(buf1),wrt)
 DO  i = ieltab,neltab,lentry
   buf( 1) = z(i+13)
   rbuf(2) = rz(i+3)
   rbuf(3) = rz(i+2)
   rbuf(4) = rz(i+1)
   rbuf(5) = rbuf(2) + rbuf(3) + rbuf(4)
   CALL WRITE (oef1x,buf(1),5,noeor)
 END DO
 CALL WRITE (oef1x,0,0,eor)
 FILE = oef1
 CALL OPEN (*1250,oef1,z(buf2),rd)
 GO TO 10
 1170 WRITE (outpt,1030)itype
 GO TO 1220
 
!     ALL PROCESSING COMPLETE.
 
 1180 mcb(1) = oef1
 CALL rdtrl (mcb)
 mcb(1) = oef1x
 CALL wrttrl (mcb)
 GO TO 1220
 
!     ERROR CONDITIONS.
 
 1190 RETURN
 
 1200 WRITE  (outpt,1210) uwm
 1210 FORMAT (a25,' 3069, OUTPUT DATA BLOCK FOR FORCES IS PURGED.')
 1220 CALL CLOSE (oef1,clsrew)
 CALL CLOSE (oef1x,clsrew)
 CALL CLOSE (ug,clsrew)
 CALL CLOSE (est,clsrew)
 CALL CLOSE (slt,clsrew)
 CALL CLOSE (dlt,clsrew)
 GO TO 1190
 1250 n = 1
 GO TO 1280
 1260 n = 2
 GO TO 1280
 1270 n = 3
 1280 CALL mesage (n,FILE,subr)
 GO TO 1220
END SUBROUTINE sdrht
