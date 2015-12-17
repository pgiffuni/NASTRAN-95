SUBROUTINE rmg
     
!     RADIATION MATRIX GENERATOR MODULE.
 
!     DMAP CALLING SEQUENCE
 
!     RMG    EST,MATPOOL,GPTT,KGGX/RGG,QGE,KGG/C,Y,TABS/C,Y,SIGMA/
!            V,N,NLR/V,N,LUSET $
 
!     THIS MODULE COMPUTES AND OUTPUTS DATA IN SINGLE OR DOUBLE
!     PRECISION BASED ON -PRECIS-.
 
 LOGICAL :: nogo     ,DOUBLE   ,lrad
 INTEGER :: buf(10)  ,subr(2)  ,radlst(2),radmtx(2),hbdytp  ,  &
     est      ,gptt     ,rgg      ,qge      ,scrt1   ,  &
     scrt2    ,scrt3    ,scrt4    ,scrt5    ,scrt6   ,  &
     sysbuf   ,outpt    ,tset     ,precis   ,unout   ,  &
     unirow   ,unnrow   ,unincr   ,pkin     ,pkout   ,  &
     pkirow   ,pknrow   ,pkincr   ,rd       ,rdrew   ,  &
     wrt      ,wrtrew   ,clsrew   ,cls      ,elem    ,  &
     z        ,core     ,buf1     ,buf2     ,buf3    ,  &
     flag     ,words    ,eor      ,mcb1(7)  ,mcb2(7) ,  &
     mcb3(7)  ,eltype   ,estwds   ,ecpt(100),rcol    ,  &
     dcol     ,rx       ,dx       ,sqr      ,FILE    ,  &
     dirgg    ,dnrgg    ,idata(16),BLOCK(20),radchk  ,  &
     mcb(7)   ,NAME(2)  ,radtyp(2),block2(20)
 REAL :: rz(2)    ,rbuf(10) ,rdata(16),ai(4)
 DOUBLE PRECISION :: dett     ,mindia   ,dsumfa   ,DO(2)    ,di(2)   ,  &
     dz(1)    ,dtemp2   ,dvalue
 CHARACTER (LEN=27) :: swm
 CHARACTER (LEN=25) :: sfm
 CHARACTER (LEN=29) :: uim
 CHARACTER (LEN=25) :: uwm
 CHARACTER (LEN=23) :: ufm
 COMMON /xmssg /  ufm      ,uwm      ,uim      ,sfm      ,swm
 COMMON /system/  ksystm(65)
 COMMON /names /  rd       ,rdrew    ,wrt      ,wrtrew   ,clsrew  ,  &
     cls      ,skip5(5) ,sqr
 COMMON /zblpkx/  ao(4)    ,irow
 COMMON /packx /  pkin     ,pkout    ,pkirow   ,pknrow   ,pkincr
 COMMON /unpakx/  unout    ,unirow   ,unnrow   ,unincr
 COMMON /gfbsx /  jl(7)    ,ju(7)    ,jb(7)    ,jx(7)    ,nzzz    ,  &
     ipr      ,isgn
 COMMON /dcompx/  ia(7)    ,il(7)    ,iu(7)    ,isr1     ,isr2    ,  &
     isr3     ,dett     ,ipow     ,nzz      ,mindia  , ib       ,ibbar
 COMMON /gpta1 /  nelems   ,last     ,incr     ,elem(1)
 COMMON /zzzzzz/  z(1)
 COMMON /BLANK /  tabs     ,sigma    ,nlr      ,luset
 EQUIVALENCE      (z(1),rz(1),dz(1) ),(buf(1),rbuf(1)   ),  &
     (DO(1)   ,ao(1)   ),(di(1) ,ai(1)     ),  &
     (idata(1),rdata(1)),(defalt,ideflt    ),  &
     (ksystm( 1),sysbuf),(ksystm( 2),outpt ),  &
     (ksystm(10),tset  ),(ksystm(55),iprec ),  &
     (ksystm(57),myradm),(ksystm(58),radchk)
 
!     MYRADM  = 1  IMPLIES SYMMETRIC SCRIPT-AF INPUT
!     RADCHK NE 0  REQUESTS DIAGNOSTIC PRINTOUT OF AREAS AND VIEW FACTOR
!     MYRADM  = 2  IMPLIES UNSYMMETRIC SCRIPT-AF INPUT
 
 DATA     subr  / 4HRMG ,4H     /
 DATA     radtyp/ 4H    ,4H  un /
 DATA     radlst/ 2014,  20     /
 DATA     radmtx/ 3014,  30     /
 DATA     hbdytp/ 52            /
 DATA     noeor / 0  /, eor / 1 /
 DATA     est   , matpol, gptt, kggx, rgg, qge, kgg   /  &
     101   , 102   , 103 , 104 , 201, 202, 203   /
 DATA     scrt1 , scrt2, scrt3, scrt4, scrt5, scrt6   /  &
     301   , 302  , 303  , 304  , 305  , 306     /
 
!     DEFINITION OF CORE AND BUFFER POINTERS
 
 CALL delset
 scrt1  = 301
 precis = 2
 IF (iprec /= 2) precis = 1
 core = korsz(z)
 buf1 = core - sysbuf - 2
 buf2 = buf1 - sysbuf - 2
 buf3 = buf2 - sysbuf - 2
 core = buf3 - 1
 IF (core < 100) CALL mesage (-8,0,subr)
 nogo   = .false.
 DOUBLE = .false.
 IF (precis == 2) DOUBLE = .true.
 IF (myradm == 1 .OR. myradm == 2) WRITE (outpt,5) uwm, radtyp(myradm)
 5 FORMAT (a25,' 2358, ',a4,'SYMMETRIC SCRIPT-AF MATRIX (HREE) ',  &
     'ASSUMED IN RADMTX')
 
!     OPEN MATPOOL DATA BLOCK.
 
 FILE = matpol
 CALL preloc (*1100,z(buf1),matpol)
 
!     LOCATE RADLST DATA
 
 CALL locate (*1090,z(buf1),radlst,flag)
 
!     BUILD ELEMENT DATA TABLE.  -LENTRY- WORDS PER ELEMENT ID PRESENT
!     IN RADLST.
 
!     EACH ENTRY CONTAINS THE FOLLOWING OR MORE
 
!     WORD  1 = ELEMENT ID OF HBDY ELEMENT
!     WORD  2 = DIAGONAL MATRIX ELEMENT A-SUB-I
!     WORD  3 = DIAGONAL MATRIX ELEMENT E-SUB-I
!     WORD  4 = ELEMENT FA SUM (USED FOR RADMTX CHECK)
!     WORD  5 = SIL-1
!     WORD  6 = SIL-2
!     WORD  7 = SIL-3
!     WORD  8 = SIL-4
!     WORD  9 = GIJ-1  (GIJ TERMS MAY BE 2 WORDS EACH IF DOUBLE PREC)
!     WORD 10 = GIJ-2
!     WORD 11 = GIJ-3
!     WORD 12 = GIJ-4
 
 
 lentry = 8 + 4*precis
 ieltab = 1
 idxm8  = ieltab - lentry - 1
 neltab = ieltab - 1
 10 IF (neltab+lentry > core) CALL mesage (-8,0,subr)
 CALL READ (*1110,*1120,matpol,z(neltab+1),1,noeor,words)
 IF (z(neltab+1) > 0.0) THEN
   GO TO    20
 ELSE
   GO TO    30
 END IF
 20 z(neltab+2) = 0
 z(neltab+3) = 0
 neltab = neltab + lentry
 GO TO 10
 
!     ALL RADLST DATA NOW IN CORE.
!     (POSITION TO END OF RECORD ON  MATPOOL)
 
 30 CALL READ (*1110,*50,matpol,buf,1,eor,words)
 WRITE  (outpt,40) swm
 40 FORMAT (a27,' 3071, EXTRA DATA IN RADLST RECORD OF MATPOOL DATA ',  &
     'BLOCK IGNORED.')
 
!     LOCATE RADMTX DATA
 
 50 NE = (neltab-ieltab+1)/lentry
 CALL locate (*135,z(buf1),radmtx,flag)
 lrad = .true.
 
!     READ IN RADMTX DATA.  FOR LOWER TRIANGLE COLUMNS PRESENT
!     ENTRY WORDS 2 AND 3 IN -ELTAB- WILL BE USED TO STORE FIRST
!     AND LAST LOCATIONS OF LOWER TRIANGLE COLUMN.  ZEROS IMPLY COLUMN
!     IS NULL.
 
 irad = neltab + 1
 
!     READ COLUMN INDEX
 
 60 CALL READ (*1110,*140,matpol,INDEX,1,noeor,words)
 
!     MAXIMUM NUMBER OF INPUT TERMS FOR THIS COLUMN. (LOWER TRIANGLE)
 
 MAX = NE - INDEX + 1
 IF (myradm == 2) MAX = NE
 
!     SET -IDX- TO ELTAB ENTRY
 
 idx = idxm8 + INDEX*lentry
 
!     READ IN COLUMN ELEMENTS IF ANY
 
 n = 0
 70 CALL READ (*1110,*1120,matpol,z(irad),1,noeor,words)
 IF (z(irad) == -1) GO TO 100
 n = n + 1
 irad = irad + 1
 IF (irad > core) CALL mesage (-8,0,subr)
 IF (n <= MAX) GO TO 70
 
!     TOO MANY COLUMN ELEMENTS INPUT
 
 irad = irad - 1
 
!     SKIP TO END OF COLUMN
 
 80 CALL READ (*1110,*1120,matpol,idum,1,noeor,words)
 IF (idum /= -1) GO TO 80
 WRITE  (outpt,90) uwm,INDEX,NE
 90 FORMAT (a25,' 3072, TOO MANY MATRIX VALUES INPUT VIA RADMTX BULK',  &
     ' DATA FOR COLUMN',i9,1H., /5X,'EXTRA VALUES IGNORED AS ',  &
     'MATRIX SIZE IS DETERMINED TO BE OF SIZE',i9,  &
     ' FROM RADLST COUNT OF ELEMENT ID-S.')
 
!     ALL DATA FOR LOWER TRIANGLE PORTION OF COLUMN IS IN CORE.
!     (BACK UP OVER ANY ZEROS)
 
 100 IF (n > 0) THEN
   GO TO   110
 ELSE
   GO TO    60
 END IF
 110 IF (z(irad-1) == 0.0) THEN
   GO TO   120
 ELSE
   GO TO   130
 END IF
 120 n = n - 1
 irad = irad - 1
 GO TO 100
 
!     SET FIRST AND LAST POINTERS
 
 130 z(idx+2) = irad - n
 z(idx+3) = irad - 1
 
!     GO READ NEXT COLUMN
 
 GO TO 60
 
!     NULL RADMTX ASSUMED
 
 135 lrad = .false.
 
!     RADMTX IS COMPLETELY IN CORE IN TEMPORARY SPECIAL PACKED FORM.
 
!     NOW PACK OUT EACH COLUMN OF MATRIX F TO SCRATCH 1
 
 140 CALL CLOSE (matpol,clsrew)
 IF (myradm == 1 .OR. myradm == 2) scrt1 = 303
 CALL gopen (scrt1,z(buf1),wrtrew)
 CALL makmcb (mcb1,scrt1,NE,sqr,precis)
 DO  jcol = 1,NE
   
!     INITIALIZE PACKING OF COLUMN -JCOL-
   
   CALL bldpk (1,precis,scrt1,0,0)
   
!     PACK OUT ELEMENTS OF COLUMN -JCOL-
   
   inxcol = idxm8 + jcol*lentry
   
!     SET FA SUM TO ZERO FOR CURRENT COLUMN.
   
   sumfa = 0.0
   IF (.NOT.lrad) GO TO 205
   DO  irow = 1,NE
     
!     LOCATE ELEMENT ROW-IROWK, COL-JCOL.
     
     IF (irow >= jcol .OR. myradm == 2) GO TO 180
     
!     HERE IF ABOVE THE DIAGONAL
!     ELEMENT DESIRED IS IN COLUMN -IROW- IN CORE AND POSITION
!     (JCOL-IROW+1) OF THE LOWER TRIANGLE PORTION.
     
     idx = idxm8 + irow*lentry
     i1  = z(idx+2)
     IF (i1 > 0) THEN
       GO TO   150
     ELSE
       GO TO   200
     END IF
     150 i2  = z(idx+3)
     ipos= jcol - irow + i1
     160 IF (ipos > i2) CYCLE
     IF (rz(ipos) < 0.0) THEN
       GO TO   162
     ELSE IF (rz(ipos) == 0.0) THEN
       GO TO   200
     ELSE
       GO TO   170
     END IF
     162 WRITE  (outpt,164) uwm,jcol,irow,rz(ipos)
     164 FORMAT (a25,' 2359, COL',i6,', ROW',i6,  &
         ' OF RADMTX IS NEGATIVE (',e14.6,').')
     170 ao(1) = rz(ipos)
     IF (myradm == 1 .OR. myradm == 2) ao(1) = -sigma*rz(ipos)
     IF (jcol == irow .AND. (myradm == 1 .OR. myradm == 2)) GO TO 175
     sumfa = sumfa + rz(ipos)
     175 CALL zblpki
     CYCLE
     
!     HERE IF BELOW OR ON DIAGONAL.
!     ELEMENT DESIRED IS IN COLUMN -JCOL- IN POSITION (IROW-JCOL+I1)
     
     180 idx= inxcol
     i1 = z(idx+2)
     IF (i1 > 0) THEN
       GO TO   190
     ELSE
       GO TO   200
     END IF
     190 i2 = z(idx+3)
     ipos = irow - jcol + i1
     IF (myradm == 2) ipos = irow + i1 - 1
     GO TO 160
     200 CONTINUE
   END DO
   
!     COMPLETE COLUMN
   
   205 CALL bldpkn (scrt1,0,mcb1)
   
!     SAVE COLUMN FA SUM IN ELTAB FOR AWHILE.
   
   rz(inxcol+4) = sumfa
   210 CONTINUE
 END DO
 
!     PACKED MATRIX IS COMPLETE
 
 CALL wrttrl (mcb1)
 CALL CLOSE (scrt1,clsrew)
!/////
!     CALL DMPFIL (-SCRT1,Z(NELTAB+1),CORE-NELTAB-2)
!     CALL BUG (10HF-MATRIX    ,210,0,1)
!/////
 
!     OUTPUT OF ELEMENT-ID LIST TO QGE HEADER RECORD IS PERFORMED AT
!     THIS TIME.
 
 FILE = qge
 CALL OPEN  (*1130,qge,z(buf1),wrtrew)
 CALL fname (qge,NAME)
 CALL WRITE (qge,NAME,2,noeor)
 DO  i = ieltab,neltab,lentry
   CALL WRITE (qge,z(i),1,noeor)
 END DO
 CALL WRITE (qge,0,0,eor)
 CALL CLOSE (qge,cls)
 
!     OPEN EST AND PROCESS EST ELEMENT DATA OF ONLY THE HBDY ELEMENTS
!     WHOSE ELEMENT ID-S ARE IN THE RADLST.  I.E. NOW IN THE RDLST TABLE
 
 FILE = est
 CALL gopen (est,z(buf1),rdrew)
 GO TO 230
 
!     LOCATE HBDY ELEMENT TYPE RECORD
 
 220 CALL fwdrec (*1110,est)
 
!     READ ELEMENT TYPE
 
 230 CALL READ (*300,*1120,est,eltype,1,noeor,words)
 IF (eltype /= hbdytp) GO TO 220
 
!     NOW POSITIONED TO READ EST DATA FOR HBDY ELEMENT.
 
 j = (eltype-1)*incr
 estwds = elem(j+12)
 lost = 0
 
!     READ EST FOR ONE ELEMENT
 
 240 CALL READ (*1110,*300,est,ecpt,estwds,noeor,words)
 
!     FIND ID IN LIST
 
 DO  i = ieltab,neltab,lentry
   IF (ecpt(1) == z(i)) GO TO 260
 END DO
 GO TO 240
 
!     ELEMENT ID IS IN LIST
 
 260 CALL hbdy (ecpt,ecpt,1,rdata,idata)
 
!     ON RETURN TAKE ELEMENT OUTPUTS AND PLANT THEM IN ALL ENTRIES
!     HAVING THIS SAME ID.
 
 iadd = 4*precis + 7
 DO  j = ieltab,neltab,lentry
   IF (ecpt(1) /= z(j)) CYCLE
   
!     CHECK TO SEE IF SUM FA/A EQUALS 1.0 FOR THIS ELEMENT.
   
   IF (rdata(2) > 1.0E-10) GO TO 261
   check = 9999999.
   GO TO 263
   261 check = rz(j+3)/rdata(2)
   IF (myradm == 1 .OR. myradm == 2) check = check/rdata(3)
   IF (check > 0.99) GO TO 262
   lost = lost + 1
   262 IF (check < 1.01) GO TO 266
   263 WRITE  (outpt,264) ufm,z(j),check,rdata(2)
   264 FORMAT (a23,' 2360, TOTAL VIEW FACTOR (FA/A), FOR ELEMENT',i9,  &
       ' IS',1P,e14.6,', (ELEMENT AREA IS ',1P,e14.5,').')
   nogo = .true.
   266 IF (check < 1.01 .AND. radchk /= 0) WRITE (outpt,267) uim,z(j),  &
       check,rdata(2)
   267 FORMAT (a29,' 2360, TOTAL VIEW FACTOR (FA/A), FOR ELEMENT',i9,  &
       ' IS ',1P,e14.6,', (ELEMENT AREA IS ',1P,e14.5,')')
   z(j  ) = idata(1)
   z(j+1) = idata(2)
   z(j+2) = idata(3)
   z(j+3) = idata(4)
   z(j+4) = idata(5)
   z(j+5) = idata(6)
   z(j+6) = idata(7)
   z(j+7) = idata(8)
   IF (DOUBLE) GO TO 270
   rz(j+ 8) = rdata( 9)
   rz(j+ 9) = rdata(10)
   rz(j+10) = rdata(11)
   rz(j+11) = rdata(12)
   GO TO 280
   270 dx = j/2 + 1
   dz(dx+4) = rdata( 9)
   dz(dx+5) = rdata(10)
   dz(dx+6) = rdata(11)
   dz(dx+7) = rdata(12)
   280 z(j) = -z(j)
 END DO
 GO TO 240
 
!     ALL ELEMENTS PROCESSED.
 
 300 CALL CLOSE (est,clsrew)
 IF (lost > 0) WRITE (outpt,302) uim,lost
 302 FORMAT (a29,' 2361, ',i4,' ELEMENTS HAVE A TOTAL VIEW FACTOR (FA',  &
     '/A) LESS THAN 0.99 , ENERGY MAY BE LOST TO SPACE.')
 
!     CHECK TO SEE IF ALL ELEMENTS WERE PROCESSED.
 
!/////
!     CALL BUG (4HELTB ,270,Z(IELTAB),NELTAB-IELTAB+1)
!/////
 DO  i = ieltab,neltab,lentry
   IF (z(i) > 0.0) THEN
     GO TO   320
   END IF
   310 z(i) = -z(i)
   CYCLE
   320 nogo = .true.
   WRITE  (outpt,330) ufm,z(i)
   330 FORMAT (a23,' 3073, NO -HBDY- ELEMENT SUMMARY DATA IS PRESENT ',  &
       'FOR ELEMENT ID =',i9, /5X, 'WHICH APPEARS ON A -RADLST- BULK DATA CARD.')
 END DO
 IF (nogo) CALL mesage (-61,0,0)
 IF (myradm == 1 .OR. myradm == 2) GO TO 345
 
!     FORMATION OF THE Y MATRIX.  MATRIX F IS STORED ON SCRATCH 1
 
!         Y    = -F  (1.0 - E )  +  A
!          IJ      IJ        J       I
 
!         A  IS ADDED IN ONLY TO THE DIAGONAL TERMS I.E. I = J
!          I
 
!     MATRIX Y WILL BE STORED ON SCRATCH 2.
 
 
!     OPEN SCRATCH 1 FOR MATRIX F COLUMN UNPACKING.
 
 CALL gopen (scrt1,z(buf1),rdrew)
 
!     OPEN SCRATCH 2 FOR MATRIX Y COLUMN PACKING
 
 CALL gopen  (scrt2,z(buf2),wrtrew)
 CALL makmcb (mcb2,scrt2,NE,sqr,precis)
 
!     SET UP VECTOR CORE (INSURE EVEN BOUNDARY)
 
 345 icol = MOD(neltab,2) + neltab + 1
 rcol = icol
 dcol = icol/2 + 1
 ncol = icol + precis*NE - 1
 IF (ncol > core) CALL mesage (-8,0,subr)
 IF (myradm == 1 .OR. myradm == 2) GO TO 465
 meltab = ieltab - lentry - 1
 
!     SETUP /PACKX/ FOR PACKING COLUMNS OF Y (SCRATCH 2)
 
 pkin   = precis
 pkout  = precis
 pkirow = 1
 pknrow = NE
 pkincr = 1
 
!     SETUP /UNPAKX/ FOR UNPACKING COLUMNS OF F (SCRATCH 1)
 
 unout  = precis
 unirow = 1
 unnrow = NE
 unincr = 1
 DO  i = 1,NE
   meltab = meltab + lentry
   rx = rcol
   dx = dcol
   
!     UNPACK A COLUMN OF F INTO CORE.
   
   CALL unpack (*350,scrt1,z(icol))
   GO TO 370
   350 DO  j = icol,ncol
     z(j) = 0
   END DO
   
!     COMPUTE THE Y-COLUMN
   
   370 DO  irow = 1,NE
     IF (DOUBLE) GO TO 380
     
!     REAL COMPUTATION
     
     rz(rx) = -rz(rx)*(1.0E0 - rz(meltab+3))
     IF (irow == i) rz(rx) = rz(rx) + rz(meltab+2)
     rx = rx + 1
     CYCLE
     
!     DOUBLE PRECISION COMPUTATION
     
     380 dz(dx) = -dz(dx)*(1.0D0 - DBLE(rz(meltab+3)))
     IF (irow == i) dz(dx) = dz(dx) + DBLE(rz(meltab+2))
     dx = dx + 1
   END DO
   
!     PACK COLUMN OUT
   
   mcbsav  = mcb2(6)
   mcb2(6) = 0
   CALL pack (z(icol),scrt2,mcb2)
   IF (mcb2(6) > 0) THEN
     GO TO   420
   END IF
   400 nogo = .true.
   WRITE  (outpt,410) ufm,i
   410 FORMAT (a23,' 3074, COLUMN',i9,' OF THE Y MATRIX IS NULL.')
   420 mcb2(6) = MAX0(mcb2(6),mcbsav)
   
 END DO
 IF (nogo) CALL mesage (-61,0,subr)
 CALL CLOSE (scrt1,clsrew)
 CALL wrttrl (mcb2)
 CALL CLOSE (scrt2,clsrew)
!/////
!     CALL DMPFIL (-SCRT2,Z(ICOL),CORE-ICOL-1)
!     CALL BUG (10HY-MATRIX    ,400,0,1)
!/////
 
!     NOW SOLVING FOR MATRIX X ON SCRATCH-3
 
!     (Y) (X) = (F)
 
!     F IS ON SCRATCH 1
!     Y IS ON SCRATCH 2
 
 
!     SETUP /DCOMPX/
 
 ia(1) = scrt2
 il(1) = 201
 iu(1) = 203
 il(5) = precis
 isr1  = scrt4
 isr2  = scrt5
 isr3  = scrt6
 CALL rdtrl (ia)
 nzz   = korsz(z(icol))
 ib    = 0
 ibbar = 0
 CALL decomp (*440,z(icol),z(icol),z(icol))
 GO TO 460
 440 WRITE  (outpt,450) ufm
 450 FORMAT (a23,' 3075, INTERMEDIATE MATRIX Y IS SINGULAR.')
 CALL mesage (-61,0,subr)
 
!     SETUP /GFBSX/
 
 460 jl(5) = il(5)
 ju(7) = iu(7)
 jl(1) = 201
 ju(1) = 203
 jb(1) = scrt1
 jx(1) = scrt3
 ipr   = precis
!//// WHAT ABOUT IDET
 isgn  = 1
 nzzz  = nzz
 jl(3) = NE
 jx(5) = precis
 CALL rdtrl (jb(1))
 CALL gfbs (z(icol),z(icol))
 jx(3) = NE
 jx(4) = sqr
 CALL wrttrl (jx)
!/////
!     CALL DMPFIL (-SCRT3,Z(ICOL),CORE-ICOL-1)
!     CALL BUG (10HX-MATRIX     ,438,0,1)
!/////
 
!     FORMATION OF THE R MATRIX (TO BE STORED ON SCRATCH 1)
 
!          R    =(-SIGMA*E *A *E *X  ) + (SIGMA*E *A )
!           IJ            J  I  I  IJ            J  I
 
!     (TERM2 IS ADDED IN ONLY WHEN I = J)
 
!     IF MYRADM = 1 OR 2    , RADMTX MULTIPLIED BY -SIGMA IS ON SCRT3
!     MATRIX X IS ON SCRATCH 3
 
 
!     OPEN SCRATCH 3 FOR MATRIX X COLUMN UNPACKING.
 
 465 CALL gopen (scrt3,z(buf3),rdrew)
 
!     THE FOLLOWING CARD IS NEEDED IF DIRECT SCRIPT-F INPUT IS USED
 
 scrt1 = 301
 
!     OPEN SCRATCH 1 FOR MATRIX R COLUMN PACKING.
 
 FILE = scrt1
 CALL gopen  (scrt1,z(buf1),wrtrew)
 CALL makmcb (mcb1,scrt1,NE,sqr,precis)
 meltab = ieltab - lentry - 1
 
!     SETUP /PACKX/ FOR PACKING COLUMNS OF R (SCRATCH 1)
 
 pkin   = precis
 pkout  = precis
 pkirow = 1
 pknrow = NE
 pkincr = 1
 
!     SETUP /UNPAKX/ FOR UNPACKING COLUMNS OF X (SCRATCH 3)
 
 unout  = precis
 unirow = 1
 unnrow = NE
 unincr = 1
 
 idx1 = ieltab - lentry
 DO  icolum = 1,NE
   dsumfa = 0.
   sumfa  = 0.
   meltab = meltab + lentry
   
!     COMPUTE CONSTANT FOR COLUMN
   
!     TEMP1 = SIGMA*E
!                    J
   
   temp1 = sigma*rz(meltab+3)
   rx = rcol
   dx = dcol
   
!     UNPACK A COLUMN OF X INTO CORE.
   
   CALL unpack (*470,scrt3,z(icol))
   GO TO 490
   470 DO  j = icol,ncol
     z(j) = 0
   END DO
   
!     COMPUTE THE R-COLUMN
   
   490 idx2 = idx1
   DO  irow = 1,NE
     idx2 = idx2 + lentry
     IF (DOUBLE) GO TO 500
     
!     REAL COMPUTATION.
     
     IF (myradm == 1 .OR. myradm == 2) GO TO 495
     temp2  = temp1*rz(idx2+1)
     rz(rx) =-temp2*rz(idx2+2)*rz(rx)
     IF (irow == icolum) rz(rx) = rz(rx) + temp2
     495 IF (irow /= icolum) sumfa  = sumfa  + rz(rx)
     rx = rx + 1
     CYCLE
     
!     DOUBLE PRECISON COMPUTATION
     
     500 IF (myradm == 1 .OR. myradm == 2) GO TO 505
     dtemp2 = DBLE(temp1)*DBLE(rz(idx2+1))
     dz(dx) =-dtemp2*DBLE(rz(idx2+2))*dz(dx)
     IF (irow == icolum) dz(dx) = dz(dx) + dtemp2
     505 IF (irow /= icolum) dsumfa = dsumfa + dz(dx)
     dx = dx + 1
     
   END DO
   
!     PACK COLUMN OF R OUT
   
   IF (myradm /= 1 .AND. myradm /= 2) GO TO 515
   IF (DOUBLE) dz(dx-1-NE+icolum) = -dsumfa
   IF (.NOT.DOUBLE) rz(icolum+rx-1-NE) = -sumfa
   515 CALL pack (z(icol),scrt1,mcb1)
 END DO
 CALL wrttrl (mcb1)
 CALL CLOSE (scrt1,clsrew)
 CALL CLOSE (scrt3,clsrew)
!/////
!     CALL DMPFIL (-SCRT1,Z(ICOL),CORE-ICOL-1)
!     CALL BUG (10HR-MATRIX    ,490,0,1)
!/////
 
!     ALL OF THE HBDY ELEMENTS OF THE RADLST HAVE
!     HAD THEIR G TERMS COMPUTED, THESE G TERMS MAY BE INSERTED INTO
!     THE FULL MATRIX G.
 
!     GOING THROUGH THE RADLST TABLE WE HAVE EACH ELEMENT ENTRY FORMING
!     A COLUMN OF G WITH THE G TERMS OF THE RESPECTIVE ENTRY BEING
!     ENTERED INTO THE COLUMN AT THE SIL LOCATIONS.  (THE SILS WERE
!     PLACED IN THE RADLST ENTRY EARLIER)
 
 
!     AS THE X MATRIX STORED ON SCRATCH 3 IS NO LONGER NEEDED
!     WE WILL USE SCRATCH 3 FOR THE G MATRIX NOW.
 
 CALL gopen  (scrt3,z(buf3),wrtrew)
 CALL makmcb (mcb3,scrt3,luset,2,precis)
 
!     LOOP ON THE RADLST TABLE
 
 DO  i = ieltab,neltab,lentry
   
!     BEGIN PACKING A COLUMN OUT
   
   CALL bldpk (precis,precis,scrt3,0,0)
   
!     PACK 1 TO 4 TERMS OUT.
   
   i1 = i + 4
   i2 = i + 7
   DO  j = 1,4
     
!     PICKING THE SMALLEST SIL NOT ZERO FOR THE NEXT TERM OUT
     
     isil = 0
     DO  l = i1,i2
       IF (z(l) > 0.0) THEN
         GO TO   530
       ELSE
         GO TO   560
       END IF
       530 IF (isil > 0) THEN
         GO TO   540
       ELSE
         GO TO   550
       END IF
       540 IF (z(l)-isil > 0.0) THEN
         GO TO   560
       END IF
       550 isil = z(l)
       k = l
       560 CONTINUE
     END DO
     
!     ZERO SIL IMPLYS OUT OF VALUES
     
     IF (isil > 0) THEN
       GO TO   570
     ELSE
       GO TO   590
     END IF
     
!     PACK OUT TERM  (MAY BE SINGLE OR DOUBLE PRECISON)
     
     570 irow = z(k)
     z(k) = 0
     
!     RESET K TO GIJ TERM PTR.
     
     kk = k + 4
     IF (DOUBLE) kk = kk + k - i1
     ao(1) = rz(kk  )
     ao(2) = rz(kk+1)
     CALL zblpki
   END DO
   
!     COMPLETE THE COLUMN
   
   590 CALL bldpkn (scrt3,0,mcb3)
 END DO
 
!     G MATRIX IS COMPLETE ON SCRATCH 3.
 
 CALL wrttrl (mcb3)
 CALL CLOSE  (scrt3,clsrew)
!/////
!     CALL DMPFIL (-SCRT3,Z(ICOL),CORE-ICOL-1)
!     CALL BUG (10HG-MATRIX     ,570,0,1)
!/////
 
!     FORM OUTPUT MATRIX  (Q  ) = (G)(R )
!                           GE         E
 
 
!     ALL CORE AT THIS POINT IS AVAILABLE THUS OPEN CORE FOR SSG2B
!     WHICH IS IN /SSGB2/ MAY BE AT THE SAME LEVEL AS
!     /RMGZZZ/.  SSG2B IS THE DRIVER FOR MPYAD.
 
 CALL ssg2b (scrt3,scrt1,0,scrt5,0,precis,1,scrt2)
 
!                                        T
!     FORM OUTPUT MATRIX  (R  ) = (Q  )(G )
!                           GG      GE
 
 
!     THE MATRIX G IS FIRST TRANSPOSED.
 
!     MATRIX G IS ON SCRATCH-3.  MATRIX G TRANSPOSE WILL BE ON SCRATCH-2
 
!     OPEN CORE /DTRANX/ FOR TRANP1 MAY BE AT SAME LEVEL AS /RMGZZZ/.
 
 CALL tranp1 (scrt3,scrt2,4,scrt4,scrt6,scrt1,rgg,0,0,0,0)
!/////
!     CALL DMPFIL (-SCRT2,Z(ICOL),CORE-ICOL-1)
!     CALL BUG (10HG-TRANSP    ,570,0,1)
!/////
 
!     SSG2B MAY BE CALLED NOW TO COMPUTE (R  )
!                                          GG
 
 CALL ssg2b (scrt5,scrt2,0,rgg,0,precis,1,scrt1)
 
!     QGE WAS PLACED ON SCRT5.  NOW COPY IT TO QGE (WHERE THE HEADER
!     RECORD HAS BEEN SPECIALLY PREPARED EARLIER) .
 
 FILE = qge
 CALL OPEN (*1130,qge,z(buf1),wrt)
 FILE = scrt5
 CALL gopen  (scrt5,z(buf2),rdrew)
 CALL cpyfil (scrt5,qge,z,core,icount)
 mcb(1) = scrt5
 CALL rdtrl (mcb)
 mcb(1) = qge
 CALL wrttrl (mcb)
 CALL CLOSE (scrt5,clsrew)
 CALL CLOSE (qge,clsrew)
 
!                    1      3
!     FORM  S   = 4(U  + T )  THIS IS ACTUALLY A DIAGONAL MATRIX.
!            GG      G    A
 
!     NOW ALLOCATE S   DIAGONAL MATRIX SPACE AND STORE -TABS- EVERYWHERE
!                   GG
 
 
 isgg = 1
 nsgg = precis*luset
 IF (nsgg > core) CALL mesage (-8,0,subr)
 IF (DOUBLE) GO TO 620
 
!     REAL VECTOR
 
 DO  i = isgg,nsgg
   rz(i) = tabs
 END DO
 GO TO 640
 
!     DOUBLE PRECISION VECTOR
 
 620 dx = isgg/2 + 1
 ndx = dx + luset - 1
 DO  i = dx,ndx
   dz(i) = tabs
 END DO
 
!     IF -TSET- IS SPECIFIED THEN THAT SET OF TEMPERATURES IS ADDED TO
!     THE UG VECTOR IN CORE.
 
 640 IF (tset > 0.0) THEN
   GO TO   650
 ELSE
   GO TO   900
 END IF
 
!     TSET IS REQUESTED
 
 650 FILE = gptt
 CALL OPEN (*1130,gptt,z(buf1),rdrew)
 
!     DETERMINE NUMBER OF RECORDS IN ELEMENT TEMPERATURE SECTION TO
!     SKIP OVER. (FIRST SKIP THE NAME IN HEADER)
 
 CALL READ (*1110,*1120,gptt,buf,2,noeor,flag)
 
!     LOOK FOR REQUESTED TSET POINTERS AND REPOSITION GPTT.
 
 NUMBER =  0
 numtst = -1
 660 CALL READ (*1110,*670,gptt,buf,3,noeor,flag)
 IF (buf(3) > NUMBER) NUMBER = buf(3)
 IF (tset /= buf(1)) GO TO 660
 
!     BUF(1)=SET-ID, BUF(2)=-1 OR DEFAULT TEMP, BUF(3)=GPTT DATA RECORD.
 
 defalt = rbuf(2)
 numtst = buf(3)
 GO TO 660
 
!     CHECK FOR TSET NOT FOUND.
 
 670 IF (numtst == -1) GO TO 1170
 
!     ADD SKIP COUNTS (EL. RECORDS + DUPE HEADER + TEMP SET -1)
 
 NUMBER = NUMBER + numtst
 
!     NO NEED TO DO FURTHER I/O IF TSET IS ALL DEFAULT TEMPS.
 
 IF (numtst == 0) NUMBER = 0
 IF (NUMBER) 740,740,720
 720 DO  i = 1,NUMBER
   CALL fwdrec (*1110,gptt)
 END DO
 
!     TEMPERATURE DATA IS IN PAIRS OF INTERNAL ID AND TEMPERATURE.
 
 
!     AT THIS POINT THE GRID POINT TEMPERATUE DATA IS ADDED INTO THE SGG
!     DIAGONAL HELD IN CORE.
 
 740 nsil = 1
 rx = isgg - 1
 dx = isgg/2
 ASSIGN 750 TO iretrn
 IF (NUMBER) 790,790,750
 750 CALL READ (*1110,*870,gptt,buf,2,noeor,flag)
 760 IF (buf(1)-nsil < 0.0) THEN
   GO TO   770
 ELSE IF (buf(1)-nsil == 0.0) THEN
   GO TO   820
 ELSE
   GO TO   800
 END IF
 770 WRITE  (outpt,780) sfm
 780 FORMAT (a25,' 3076, GPTT DATA IS NOT IN SORT BY INTERNAL ID.')
 CALL mesage (-61,0,subr)
 
!     ADD DEFAULT TEMPERATURE (IF ONE EXISTS) TO THOSE POINTS NOT HAVING
!     AN EXPLICIT TEMPERATURE DEFINED.
 
 790 buf(1) = luset + 1
 800 IF (ideflt /= -1) GO TO 830
 WRITE  (outpt,810) ufm,nsil
 810 FORMAT (a23,' 3077, THERE IS NO GRID POINT TEMPERATURE DATA OR ',  &
     'DEFAULT TEMPERATURE DATA FOR SIL POINT',i9, /5X,  &
     'AND POSSIBLY OTHER POINTS.')
 CALL mesage (-61,0,subr)
 820 value = rbuf(2)
 GO TO 840
 830 value = defalt
 ASSIGN 880 TO iretrn
 isil = nsil
 ksil = buf(1) - 1
 nsil = buf(1)
 GO TO 845
 840 isil = nsil
 ksil = buf(1)
 nsil = buf(1) + 1
 845 DO  i = isil,ksil
   IF (i > luset) GO TO 890
   IF (DOUBLE) GO TO 850
   rx = rx + 1
   rz(rx) = rz(rx) + value
   CYCLE
   850 dx = dx + 1
   dz(dx) = dz(dx) + DBLE(value)
 END DO
 GO TO iretrn, (750,890,880)
 870 ASSIGN 890 TO iretrn
 buf(1) = luset
 value  = defalt
 GO TO 840
 880 ASSIGN 750 TO iretrn
 GO TO 760
 
!     ALL TEMPERATURE DATA HAS BEEN ADDED IN.
 
 890 CALL CLOSE (gptt,clsrew)
!/////
!     CALL BUG (4HTMPS,890,Z(ISGG),NSGG-ISGG+1)
!/////
 
!     NOW CUBE EACH TERM AND THEN MULTIPLY EACH TERM BY 4.0
 
 900 IF (DOUBLE) GO TO 920
 
!     REAL COMPUTATION
 
 DO  i = isgg,nsgg
   rz(i) = 4.0*(rz(i)**3)
 END DO
 GO TO 940
 
!     DOUBLE PRECISION COMPUTATION
 
 920 dx  = isgg/2 + 1
 ndx = dx + luset - 1
 DO  i = dx,ndx
   dz(i) = 4.0D0*(dz(i)**3)
 END DO
 
!     ALLOCATION OF CORE FOR A COLUMN OF MATRIX RGG.
 
 940 irgg  = nsgg + 1
 nrgg  = irgg + precis*luset - 1
 dirgg = irgg/2 + 1
 dnrgg = dirgg + luset - 1
 IF (nrgg > core) CALL mesage (-8,0,subr)
 
!                                   X
!     FORM OUTPUT MATRIX  (K  ) = (K  ) + (R  )(S  )
!                           GG      GG      GG   GG
 
 
!     THE DIAGONAL MATRIX (S  ) RESIDES IN CORE FROM Z(ISGG) TO Z(NSGG)
!                           GG
 
!     Z(IRGG) TO Z(NRGG) WILL BE USED TO HOLD A COLUMN OF R.
 
!       X
!     (K  ) WILL BE UNPACKED INCREMENTALLY AND ADDED INTO THE COLUMN
!       GG
 
!     OF R, AFTER THAT COLUMN OF R HAS BEEN MULTIPLIED BY THE RESPECTIVE
 
!     DIAGONAL ELEMENT OF (S  ).
!                           GG
 
!/////
!     CALL BUG (4HSGG  ,829,Z(ISGG),NSGG-ISGG+1)
!/////
 CALL gopen (rgg,z(buf1),rdrew)
 CALL gopen (kggx,z(buf2),rdrew)
 CALL gopen (kgg,z(buf3),wrtrew)
 CALL makmcb (mcb1,kgg,luset,sqr,precis)
 
!     SET UP /PACKX/ FOR PACKING COLUMN OF KGG OUT.
 
 pkin   = precis
 pkout  = precis
 pkirow = 1
 pknrow = luset
 pkincr = 1
 rx = isgg - 1
 dx = isgg/2
 
!     LOOP THROUGH -LUSET- COLUMNS TO BE OUTPUT.
 
 DO  i = 1,luset
   IF (DOUBLE) GO TO 950
   rx = rx + 1
   value = rz(rx)
   GO TO 960
   950 dx = dx + 1
   dvalue = dz(dx)
   
!     UNPACK A COLUMN OF R
   
   960 DO  j = irgg,nrgg
     z(j) = 0
   END DO
   
!     -UNPACK- CAN NOT BE USED HERE DUE TO UNPACKING OF KGGX BELOW.
   
   CALL intpk  (*990,rgg,block2,precis,1)
   984 CALL intpki (ai,iirow,rgg,block2,ieol)
   IF (DOUBLE) GO TO 985
   k = irgg - 1 + iirow
   rz(k) = rz(k) + ai(1)
   IF (ieol > 0) THEN
     GO TO   990
   ELSE
     GO TO   984
   END IF
   985 k = dirgg - 1 + iirow
   dz(k) = dz(k) + di(1)
   IF (ieol > 0) THEN
     GO TO   990
   ELSE
     GO TO   984
   END IF
   
!     MULTIPLY RGG COLUMN BY DIAGONAL ELEMENT OF SGG.
   
   990 IF (DOUBLE) GO TO 1010
   
!     REAL COMPUTATION
   
   DO  j = irgg,nrgg
     rz(j) = rz(j)*value
   END DO
   GO TO 1030
   
!     DOUBLE PRECISION COMPUTATION
   
   1010 DO  j = dirgg,dnrgg
     dz(j) = dz(j)*dvalue
   END DO
   
!     INCREMENTAL UNPACK OF A COLUMN OF KGGX.
!     ADD TO MODIFIED COLUMN OF RGG IN CORE, AND THEN
!     BLAST PACK OUT FURTHER MODIFIED COLUMN AS A COLUMN OF KGG.
   
!     START UNPACKING COLUMN OF KGGX
   
   1030 CALL intpk  (*1070,kggx,BLOCK,precis,1)
   1040 CALL intpki (ai,iirow,kggx,BLOCK,ieol)
   
!     ADD VALUE IN
   
   IF (iirow > luset) GO TO 1050
   IF (DOUBLE) GO TO 1060
   
!     REAL ADD IN
   
   k = irgg - 1 + iirow
   rz(k) = rz(k) + ai(1)
   1050 IF (ieol > 0) THEN
     GO TO  1070
   ELSE
     GO TO  1040
   END IF
   
!     DOUBLE PRECISION ADD IN
   
   1060 k = dirgg - 1 + iirow
   dz(k) = dz(k) + di(1)
   IF (ieol > 0) THEN
     GO TO  1070
   ELSE
     GO TO  1040
   END IF
   
!     PACK OUT COMPLETED COLUMN.
   
   1070 CALL pack (z(irgg),kgg,mcb1)
 END DO
 CALL wrttrl (mcb1)
 CALL CLOSE (kgg,clsrew )
 CALL CLOSE (kggx,clsrew)
 CALL CLOSE (rgg,clsrew )
 
!     ALL PROCESSING COMPLETED.
 
 nlr = +1
 RETURN
 1090 CALL CLOSE (matpol,clsrew)
 1100 nlr = -1
 RETURN
 
!     ERROR CONDITIONS
 
 
!     END OF FILE
 
 1110 j = -2
 GO TO 1140
 
!     END OF RECORD
 
 1120 j = -3
 GO TO 1140
 
!     UNDEFINED FILE
 
 1130 j = -1
 1140 CALL mesage (j,FILE,subr)
 
!     GPTT DATA MISSING FOR SET -TSET-.
 
 1170 WRITE  (outpt,1180) ufm,tset
 1180 FORMAT (a23,' 3078, NO GPTT DATA IS PRESENT FOR TEMPERATURE SET ',  &
     i8,1H.)
 CALL mesage (-61,0,subr)
 
!     NO HBDY ELEMENTS
 
 
 RETURN
END SUBROUTINE rmg
