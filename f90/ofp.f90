SUBROUTINE ofp
     
!     THE OUTPUT FILE PROCESSOR
 
!     THIS SUBROUTINE IS THE MAIN AND ONLY DRIVER.
!     OFP1 OUTPUTS HEADINGS ONLY.
 
 IMPLICIT INTEGER (a-z)
 LOGICAL :: axic,fluid,temper,onefil,headng,solset,elemen,  &
     pnched,heat,gpfb,ese,dummy,gpst,eor,pack,strain
 INTEGER :: REAL(10),imag(5),isave(20),gse(4),i15blk(2),  &
     filex(6),b(23,4),FMT(300),iout(100),tsave(96), scan(2),id(50),buff(1),of(56)
 REAL :: freal(10),fimag(2),out(100)
 DOUBLE PRECISION :: dout(50)
!WKBI
 CHARACTER (LEN=1) :: cfmt(300)
 CHARACTER (LEN=27) :: swm
 CHARACTER (LEN=25) :: sfm
 CHARACTER (LEN=29) :: uim
 CHARACTER (LEN=25) :: uwm
 CHARACTER (LEN=23) :: ufm
 COMMON /xmssg /  ufm,uwm,uim,sfm,swm
 COMMON /BLANK /  icard,option(2)
 COMMON /system/  ksystm(65)
 COMMON /ofpcom/  temper,mpunch
 COMMON /ofpbd1/  d(1)
 COMMON /ofpbd5/  esingl(64),e(1)
!ZZ   COMMON /ZZOFPX/  CORE(1)
 COMMON /zzzzzz/  core(20000)
 COMMON /ofp1id/  idm(2)
 COMMON /output/  head(96)
 EQUIVALENCE      (ksystm( 1),sysbuf), (ksystm( 2),l     ),  &
     (ksystm( 9),maxlns), (ksystm(12),line  ),  &
     (ksystm(38),axif  ), (ksystm(56),itherm),  &
     (freal (1),REAL(1)), (fimag (1),imag(1)),  &
     (id    (3),ieltyp ), (iout(1),out(1),dout(1)),  &
     (l1, of(1),core(1)), (l2,of(2)), (l3,of(3)),  &
     (l4, of(4)), (l5,of(5)), (id(1),of(6)), (buff(1),of(56))
!WKBI
 EQUIVALENCE      ( cfmt, FMT )
 DATA    pe    /  4H1P,e /,   pf   / 4H0P,f /
 DATA    e236  /  4H23.6 /,   f236 / 4H14.1 /
 DATA    e156  /  4H15.6 /,   f156 / 4H6.1  /
 DATA    i8    /  4H,i8, /,   i12  / 4H,i11 /
 DATA    i2x   /  4H2X   /,   i2xx / 4H,2X  /
 DATA    i1x   /  4H(1X  /,   i1xx / 4H,1X  /
 DATA    istar /  4H,1H* /,   i15x / 4H/14X /
 DATA    ih0   /  4H/1H0 /,   i1h0 / 4H(1H0 /
 DATA    i9x   /  4H,9X  /,   i6x  / 4H,6X  /
 DATA    f174  /  4H17.4 /
 DATA    static,  reigen  ,   freq ,  trans  , bk1  , ceigen  /  &
     1     ,  2       ,   5    ,  6      , 8    , 9       /
 DATA    a4    ,  comma   ,   cparen, oparen  /  &
     4HA4  ,  4H,     ,   4H)   , 4H(     /
 DATA    eend  /  195    /,   i15blk/ 4HA4,  , 4H11X /,  &
     gse   /  4HG     ,   4HS   , 4HE    , 4HM   /
 DATA    filex /  101, 102,  103,104, 105    , 106   /
 DATA    iblank,  e9pt1  /    4H    , 195    /
 DATA    iheat /  4HHEAT /,   center/ 4HTER  /
 DATA    phase /  4H1P9E /,   scan  / 4HSCAN , 4HNED /
 DATA    hex1  ,  hex2, hex3 /4HHEX1, 4HHEX2 , 4HHEX3/
 
 
!     THE FOLLOWING ARE ZERO POINTERS TO THE DATA-BLOCK AND LINE-SET
!     SELECTION LISTS.  TO THIS THE SUBSET SELECTION POINTER IS ADDED.
!     THE SUBSET SELECTION POINTER IS BASED ON INFORMATION IN THE ID
!     RECORD.
 
!     A MINUS ONE IN THE FOLLOWING ARRAY INDICATES AN UNDEFINED OUTPUT.
 
!            S O R T - I    S O R T - I I
!          *************** ****************        ZERO-BASE
!          REAL    COMPLEX REAL    COMPLEX    POINTERS INTO C-ARY
!     *************************************   *******************
!     DISPLACEMENT VECTOR
 DATA b( 1,1),b( 1,2),b( 1,3),b( 1,4) /    0,   4,   2,   6 /
 
!     LOAD VECTOR
 DATA b( 2,1),b( 2,2),b( 2,3),b( 2,4) /    8,  12,  10,  14 /
 
!     SPCF VECTOR
 DATA b( 3,1),b( 3,2),b( 3,3),b( 3,4) /   16,  20,  18,  22 /
 
!     ELEMENT FORCE   (ZERO POINTERS INTO OVERLAY BLOCK DATA)
 DATA b( 4,1),b( 4,2),b( 4,3),b( 4,4) /    0,   0,   0,   0 /
 
!     ELEMENT STRESS  (ZERO POINTERS FOR OVERLAY BLOCK DATA)
 DATA b( 5,1),b( 5,2),b( 5,3),b( 5,4) /    0,   0,   0,   0 /
 
!     EIGENVALUE SUMMARY
 DATA b( 6,1),b( 6,2),b( 6,3),b( 6,4) /   38,  39,  -1,  -1 /
 
!     EIGENVECTOR
 DATA b( 7,1),b( 7,2),b( 7,3),b( 7,4) /   72,  40,  -1,  -1 /
 
!     GPST
 DATA b( 8,1),b( 8,2),b( 8,3),b( 8,4) /   73,  -1,  -1,  -1 /
 
!     EIGENVALUE ANALYSIS SUMMARY
 DATA b( 9,1),b( 9,2),b( 9,3),b( 9,4) /   64,  68,  -1,  -1 /
 
!     VELOCITY VECTOR
 DATA b(10,1),b(10,2),b(10,3),b(10,4) /   24,  30,  26,  32 /
 
!     ACCELERATION VECTOR
 DATA b(11,1),b(11,2),b(11,3),b(11,4) /   25,  34,  27,  36 /
 
!     NON-LINEAR-FORCE VECTOR
 DATA b(12,1),b(12,2),b(12,3),b(12,4) /   28,  -1,  29,  -1 /
 
!     GRID-POINT-WEIGHT-OUTPUT
 DATA b(13,1),b(13,2),b(13,3),b(13,4) /   -1,  -1,  -1,  -1 /
 
!     EIGENVECTOR (SOLUTION SET FROM VDR)
 DATA b(14,1),b(14,2),b(14,3),b(14,4) /   -1,  60,  -1,  62 /
 
!     DISP-VECTOR (SOLUTION SET FROM VDR)
 DATA b(15,1),b(15,2),b(15,3),b(15,4) /   42,  44,  43,  46 /
 
!     VELO-VECTOR (SOLUTION SET FROM VDR)
 DATA b(16,1),b(16,2),b(16,3),b(16,4) /   48,  50,  49,  52 /
 
!     ACCE-VECTOR (SOLUTION SET FROM VDR)
 DATA b(17,1),b(17,2),b(17,3),b(17,4) /   54,  56,  55,  58 /
 
!     ELEMENT STRAIN ENERGY (FROM GPFDR)
 DATA b(18,1),b(18,2),b(18,3),b(18,4) /   74,  -1,  -1,  -1 /
 
!     GRID POINT FORCE BALANCE (FROM GPFDR)
 DATA b(19,1),b(19,2),b(19,3),b(19,4) /   76,  -1,  -1,  -1 /
 
!     MPCFORCE VECTOR
 DATA b(20,1),b(20,2),b(20,3),b(20,4) /   78,  -1,  -1,  -1 /
 
!     ELEMENT STRAIN/CURVATURE (ZERO POINTER FOR OVERLAY BLOCK DATA)
 DATA b(21,1),b(21,2),b(21,3),b(21,4) /    0,  -1,  -1,  -1 /
 
!     STRESSES IN LAYERED COMPOSITE ELEMENTS (ZERO POINTER)
 DATA b(22,1),b(22,2),b(22,3),b(22,4) /    0,   0,   0,   0 /
 
!     FORCES IN LAYERED COMPOSITE ELEMENTS   (ZERO POINTER)
 DATA b(23,1),b(23,2),b(23,3),b(23,4) /    0,   0,   0,   0 /
!     ************************************   *******************
 
!     SAVE OLD TITLES WHATEVER THEY BE AND RESTORE BEFORE RETURNING
 
 CALL totape (3,buff(1))
 heat = .false.
 IF (itherm /= 0) heat = .true.
 option(1) = 0
 IF (heat) option(1) = iheat
 onefil = .false.
 GO TO 10
 
 
 ENTRY ofpdmp (ifile1)
!     =====================
 
 onefil = .true.
 10 DO  i = 1,96
   tsave(i) = head(i)
 END DO
 
 icore = korsz(buff)
 IF (icore >= sysbuf) GO TO 40
 WRITE  (6,30) uwm,icore,sysbuf
 30 FORMAT (a25,' 2043, OFP HAS INSUFFICIENT CORE FOR ONE GINO ',  &
     'BUFFER ****    OFP NOT EXECUTED.')
 RETURN
 
 40 line  = 0
 ifile = 0
 
!     LOOP FOR 6 FILES
 
 50 ifile = ifile + 1
 IF (onefil .AND. ifile > 1) GO TO 2060
 FILE = filex(ifile)
 IF (onefil) FILE = ifile1
 CALL OPEN (*2050,FILE,buff(1),0)
 from = 55
 CALL fwdrec (*2020,FILE)
 60 CALL READ (*2040,*2040,FILE,  id(1),50,0,flag)
 CALL READ (*2040,*2040,FILE,head(1),96,1,flag)
 axic   = .false.
 temper = .false.
 dummy  = .false.
 gpst   = .false.
 sort   =  1
 pnched = .false.
 headng = .false.
 gpfb   = .false.
 ese    = .false.
 strain = .false.
 
!     COMPUTE I AND J, THE B ARRAY SUBSCRIPTS
 
 j = id(2)/1000
 i = id(2) - j*1000
 j = j + 1
 IF (i /= 4 .AND. i /= 5 .AND. i /= 21) GO TO 70
 icurv = id(3)/1000
 id(3) = id(3) - 1000*icurv
 
 70 pack   = .false.
 solset = .false.
 fluid  = .false.
 IF (axif /= 0) fluid = .true.
 elemen = .false.
 iapp = id(1)/10
 nadd = 1
 from = 75
 IF (j > 4) GO TO 2020
 from = 77
 IF (j < 0) THEN
   GO TO  2020
 END IF
 80 from = 80
 IF (i < 1 .OR. i > 23) GO TO 2020
 IF (j > 2) sort = 2
 IF (j /= 3 .OR. iapp /= static) GO TO 120
 IF (head(74) == scan(1) .AND. head(75) == scan(2)) GO TO 100
 DO  ihd = 65,96
   head(ihd) = iblank
 END DO
 GO TO 120
 100 DO  ihd = 65,72
   IF (ihd >= 68) head(ihd+22) = iblank
   head(ihd) = iblank
 END DO
 120 GO TO (150,150,150,230,240,270,150,290,300,150,  &
     150,150,340,380,380,380,380,390,400,220, 410,420,420), i
 150 pack = .true.
 IF (id(3) == 1000) axic = .true.
 SELECT CASE ( i )
   CASE (    1)
     GO TO 200
   CASE (    2)
     GO TO 210
   CASE (    3)
     GO TO 220
   CASE (    4)
     GO TO 160
   CASE (    5)
     GO TO 160
   CASE (    6)
     GO TO 160
   CASE (    7)
     GO TO 280
   CASE (    8)
     GO TO 160
   CASE (    9)
     GO TO 160
   CASE (   10)
     GO TO 310
   CASE (   11)
     GO TO 320
   CASE (   12)
     GO TO 330
 END SELECT
 160 CALL mesage (-61,0,0)
 
!     DISPLACEMENT VECTOR
 
 200 IF (j == 3 .AND. iapp == trans) nadd = 7
 IF (option(1) /= iheat) GO TO 500
 IF (i == 1 .AND. (j == 1 .OR. j == 3)) temper =.true.
 GO TO 500
 
!     LOAD VECTOR
 
 210 IF (j == 3 .AND. iapp == trans) nadd = 7
 GO TO 500
 
!     SPCF VECTOR, MPCF VECTOR
 
 220 IF (j == 3 .AND. iapp == trans) nadd = 7
!WKBI 11/93 SPR93007
 pack = .true.
 GO TO 500
 
!     ELEMENT FORCE, ELEMENT STRESS
 
 230 CONTINUE
 240 from = 240
 IF (id(3) < 1 .OR.  id(3) > 100) GO TO 2020
 IF (id(3) > 52 .AND. id(3) < 62) dummy = .true.
 elemen = .true.
 iopt   = 2
 IF (icurv > 0 .AND. j == 1) GO TO 260
 IF (icurv > 0 .AND. (j == 2 .OR. j == 4)) GO TO 250
 nadd = 6*(id(3)-1) + 1
 IF (j == 2 .OR. j == 4) nadd = nadd*2 - 1
 GO TO 500
 250 from = 250
 IF (icurv > 1) GO TO 2020
 nadd = 0
 IF (id(3) ==  6) nadd =  1
 IF (id(3) == 17) nadd = 13
 IF (id(3) == 18) nadd = 25
 IF (id(3) == 19) nadd = 37
 GO TO 500
 
!     ELEMENT STRESS IN MATERIAL COORDINATE SYSTEM
 
 260 from = 260
 IF (icurv > 2) GO TO 2020
 nadd = 0
 IF (id(3) ==  6) nadd =  1
 IF (id(3) == 17) nadd =  7
 IF (id(3) == 18) nadd = 13
 IF (id(3) == 19) nadd = 19
 IF (icurv ==  2) nadd = 25
 GO TO 500
 
!     EIGENVALUE SUMMARY
 
 270 CONTINUE
 GO TO 500
 
!     EIGENVECTOR
 
 280 CONTINUE
 GO TO 500
 
!     GPST
 
 290 gpst = .true.
 GO TO 500
 
!     EIGENVALUE ANALYSIS SUMMARY
!       ID(3) = 1  DETERMINANT METHOD TABLE
!       ID(3) = 2  INVERSE POWER TABLE
!       ID(3) = 3  DETERMINANT METHOD SWEPT FUNCTION DATA VECTORS
!       ID(3) = 4  UPPER HESSENBERG METHOD TABLE
 
 300 nadd = 6*(id(3)-1) + 1
 from = 300
 IF (id(3) > 4) GO TO 2020
 GO TO 500
 
!     VELOCITY VECTOR
 
 310 CONTINUE
 GO TO 500
 
!     ACCELERATION VECTOR
 
 320 CONTINUE
 GO TO 500
 
!     NON-LINERAR FORCE VECTOR
 
 330 solset = .true.
 GO TO 500
 
!     GRID-POINT-WEIGHT-OUTPUT
!     (FROM = 345 AND 355 ARE SETUP IN OFPGPW)
 
 340 from = 340
 IF (j > 1) GO TO 2020
 CALL ofpgpw (*2020,FILE,dout,from)
 GO TO 60
 
!     EIGENVECTOR, DISPLACEMENT, VELOCITY, ACCELERATION
!     (VDR OUTPUT ONLY)
 
 380 pack   = .true.
 solset = .true.
 GO TO 500
 
!     ELEMENT STRAIN ENERGY.
 
 390 ese  = .true.
 iopt =  3
 GO TO 500
 
!     GRID POINT FORCE BALANCE.
 
 400 gpfb = .true.
 iopt =  4
 lastid = 0
 GO TO 500
 
!     ELEMENT STRAIN/CURVATURE
 
 410 from = 410
 IF (id(3) /= 6 .AND. id(3) /= 17 .AND. id(3) /= 18 .AND. &
!WKBR NCL93012 3/94     1    ID(3).NE.19) GO TO 2020  &
 id(3) /= 19 .AND. id(3) /= 64 .AND. id(3) /= 83) GO TO 2020
 from = 415
 IF (icurv > 2) GO TO 2020
 strain = .true.
 elemen = .true.
 iopt = 2
 nadd = 0
 IF (id(3) ==  6) nadd =  1
 IF (id(3) == 17) nadd =  7
 IF (id(3) == 18) nadd = 13
 IF (id(3) == 19) nadd = 19
!WKBNB NCL93012 3/94
 IF (id(3) == 64) nadd = 55
 IF (id(3) == 83) nadd = 61
!WKBNE NCL93012 3/94
 IF (icurv ==  1) nadd = nadd + 24
 IF (icurv ==  2) nadd = 49
 GO TO 500
 
!     STRESSES AND FORCES IN LAYERED COMPOSITE ELEMENTS
 
 420 CALL ofcomp (*60,FILE,j,ieltyp,iapp,headng,pnched,i)
 GO TO 60
 
 500 from = 500
 IF (b(i,j) == -1) GO TO 2020
 IF (pack ) iopt = 1
 point = nadd + b(i,j)*6
 
!     IS THIS MAGNITUDE / PHASE OUTPUT
 
 IF (id(9) == 3 .AND. (iapp == freq .OR. iapp == ceigen)) point = point + 6
 
 IF (strain) GO TO 660
 IF (elemen) GO TO 510
 
!     CALL NON-STRESS AND NON-FORCE OVERLAY.
 
 CALL ofpmis (ix,l1,l2,l3,l4,l5,point)
 from = 505
 GO TO 690
 
!     CALL PARTICULAR STRESS OR FORCE OVERLAY CONSIDERING
!     REAL, COMPLEX, SORT1, SORT2.
 
 510 IF (dummy) GO TO 580
 IF (icurv <= 0) GO TO 515
 IF (j == 1) GO TO 650
 IF (j == 2) GO TO 670
 IF (j == 4) GO TO 680
 515 itype = j + 4*(5-i)
 SELECT CASE ( itype )
   CASE (    1)
     GO TO 520
   CASE (    2)
     GO TO 530
   CASE (    3)
     GO TO 540
   CASE (    4)
     GO TO 560
   CASE (    5)
     GO TO 570
   CASE (    6)
     GO TO 610
   CASE (    7)
     GO TO 620
   CASE (    8)
     GO TO 640
 END SELECT
 
 520 CALL ofprs1 (ix,l1,l2,l3,l4,l5,point)
 from = 525
 GO TO 690
 
 530 CALL ofpcs1 (ix,l1,l2,l3,l4,l5,point)
 from = 535
 GO TO 690
 
 540 IF (iapp /= static) GO TO 550
 CALL ofrs2s (ix,l1,l2,l3,l4,l5,point)
 from = 545
 GO TO 690
 550 CALL ofprs2 (ix,l1,l2,l3,l4,l5,point)
 from = 555
 GO TO 690
 
 560 CALL ofpcs2 (ix,l1,l2,l3,l4,l5,point)
 from = 565
 GO TO 690
 
 570 CALL ofprf1 (ix,l1,l2,l3,l4,l5,point)
 from = 575
 IF (.NOT.heat .OR. id(3) == 82) GO TO 690
 IF (id(10) /= -9) GO TO 600
 l2 = 405
 l4 = 0
 l5 = 406
 id(10) = 9
 GO TO 700
 
 580 CALL odum (1,ix,itype,nmult,nlines,id)
 dummy = .false.
 from = 580
 GO TO 690
 
!     REAL FORCE SORT 1 (HEAT)
 
 600 l2 = 297
 IF (id(10) == 5) l2 = 302
 l4 = 0
 l5 = 298
 IF (id(10) == 5) l5 = 300
 GO TO 700
 
 610 CALL ofpcf1 (ix,l1,l2,l3,l4,l5,point)
 from = 615
 IF (heat) GO TO 2020
 GO TO 690
 
 620 IF (iapp /= static) GO TO 630
 CALL ofrf2s (ix,l1,l2,l3,l4,l5,point)
 from = 625
 GO TO 690
 630 CALL ofprf2 (ix,l1,l2,l3,l4,l5,point)
 from = 635
 IF (.NOT.heat .OR. id(3) == 82) GO TO 690
 
!     REAL FORCE SORT 2 (HEAT)
 
 l1 = 108
 l2 = 297
 IF (id(10) == 5) l2 = 302
 l5 = 299
 IF (id(10) == 5) l5 = 301
 GO TO 700
 
 640 CALL ofpcf2 (ix,l1,l2,l3,l4,l5,point)
 from = 645
 IF (heat) GO TO 2020
 GO TO 690
 
 650 CALL ofpss1 (ix,l1,l2,l3,l4,l5,point)
 from = 655
 GO TO 690
 660 CALL ofpsn1 (ix,l1,l2,l3,l4,l5,point)
 from = 665
 GO TO 690
 
 670 CALL ofpcc1 (ix,l1,l2,l3,l4,l5,point)
 from = 675
 GO TO 690
 
 680 CALL ofpcc2 (ix,l1,l2,l3,l4,l5,point)
 from = 685
 
 690 IF (ix == 0) GO TO 2000
 
!     IF THERMAL DISPLACEMENTS IN -HEAT- PROBLEMS, CHANGE HEADING
!     FROM  DISPLACEMENT TO TEMPERATURE
 
 700 IF (temper .AND. l2 == 1) l2 = 253
 IF (temper .AND. l3 == 1) l3 = 253
 
!     HEAT PROBLEMS REAL-SORT1-VECTORS ONLY
 
 IF (heat .AND. pack .AND. j == 1) l5 = 296
 IF (heat .AND. sort == 2 .AND. .NOT.elemen) l5 = 303
 IF (axic) l4 = -1
 IF (axic) l5 = 203
 IF (axic .AND. iapp == trans  .AND. j == 3) l5 = 402
 IF (axic .AND. iapp == static .AND. j == 3) l5 = 403
 IF (axic .AND. iapp == freq   .AND. j == 4) l5 = 404
 IF (j /= 1) GO TO 710
 IF (iapp == trans .AND. i /= 8) l1 = 106
 IF ((iapp == reigen .OR. iapp == bk1) .AND. i /= 6 .AND. i /= 8  &
     .AND. i /= 9) l1 = 102
 GO TO 720
 710 IF (j /= 2) GO TO 720
 IF (iapp == ceigen .AND. i /= 6 .AND. i /= 9) l1 = 110
 720 idd = 0
 IF (sort == 1) GO TO 730
 idd   = id(5)
 itemp = idd/10
 device= idd - 10*itemp
 device= MOD(device,8)
 idd   = itemp
 id(5) = idd
 IF (heat .AND. .NOT.elemen .AND. sort == 2) l5 = 303
 IF (iapp == static) idd = -1
 
!     SORT2 HARMONIC VECTOR OUTPUT
 
 730 IF (.NOT.pack .OR. .NOT.fluid .OR. sort == 1) GO TO 750
 IF (id(5) < 500000) GO TO 740
 IF (l1 == 107) l1 = 229
 740 fluid = .false.
 IF (ese  .OR.   gpfb) GO TO 800
 750 IF (pack .OR. elemen) GO TO 800
 CALL ofp1
 headng = .true.
 
!     OUTPUT THE DATA BLOCK
 
 800 eor = .false.
 IF (elemen  .AND. id(3) == 35) axic = .true.
 IF (elemen  .AND. id(3) == 70) axic = .true.
 IF (elemen  .AND. id(3) == 71) axic = .true.
 IF (axic) solset = .true.
 
!     D(IX) CONTINS TWO VALUES IN A PACKED 4 DIGIT NUMBER.
!     THE RIGHT 2 DIGITS GIVES THE NUMBER OF LINES FORMAT PRODUCES.
!     THE LEFT 2 DIGITS GIVES THE NUMBER OF DATA VECTORS PER LINE.
!     IF  THE LEFT 2 DIGITS ARE 0 OR NULL, 1 VECTOR IS ASSUMED.
 
 IF (heat ) GO TO 810
 IF (dummy) GO TO 820
 nmult  = d(ix)/100
 nlines = d(ix) - nmult*100
 GO TO 820
 810 nmult  = 1
 nlines = 1
 820 IF (nmult == 0) nmult = 1
 nwds  = id(10)*nmult
 nwdsav= nwds
 maxn  = maxlns - nlines
 
 IF (nwds == 0) GO TO 60
 900 IF (eor) GO TO 60
 from  = 900
 CALL READ (*2020,*910,FILE,iout(1),nwds,0,flag)
 GO TO 920
 910 IF (flag == 0) GO TO 60
 IF (flag == 0) GO TO 60
 nwds  = flag
 nwdsav= nwds
 eor   = .true.
 920 i1    = 0
 IF (axic) GO TO 930
 igse  = iout(2)
 IF (pack .OR. gpst) iout(2) = gse(igse)
 IF (ese  .OR. gpfb) GO TO 930
 IF (.NOT.pack .AND. .NOT.elemen) GO TO 1030
 930 IF (sort == 2) GO TO 990
 incr = id(10)
 i  = 1
 k1 = 1
 940 itemp  = iout(i)/10
 device = iout(i) - 10*itemp
 device = MOD(device,8)
 iout(i)= itemp
 IF (device < 4) GO TO 950
 CALL ofppun (iout(i),iout(i),incr,iopt,idd,pnched)
 device = device - 4
 950 IF (device > 0) GO TO 980
 
!     ELIMINATE VECTOR FROM MULTIPLE VECTOR PER LINE
 
 nwds = nwds - incr
 IF (nwds > i) GO TO 960
 device = 1
 IF (nwds > 0) GO TO 990
 IF (eor) GO TO 60
 nwds = nwdsav
 GO TO 900
 960 k1 = k1 + incr
 k2 = k1 + incr - 1
 jj = i - 1
 DO  j = k1,k2
   jj = jj + 1
   iout(jj) = iout(j)
 END DO
 GO TO 940
 980 i = i + incr
 IF (i <= nwds) GO TO 940
 990 IF (device < 4) GO TO 1020
 IF (elemen) GO TO 1000
 CALL ofppun (iout(1),iout(1),nwds,iopt,idd,pnched)
 GO TO 1020
 
!     SORT 2 ELEMENT PUNCH
 
 1000 incr = id(10)
 DO  jj = 1,nwds,incr
   CALL ofppun (iout(jj),iout(jj),incr,iopt,idd,pnched)
 END DO
 1020 IF (device /= 1 .AND. device /= 5) GO TO 900
 1030 IF (.NOT.pack .OR. axic) GO TO 1100
 IF (fluid .AND. sort == 1 .AND. iout(1) >= 500000) GO TO 1800
 1040 IF (iout(2) /= gse(1)) GO TO 1720
 
!     BUILD FORMAT CHECKING DATA FOR SPECIAL CASES.
 
 1100 i = 1
 IF (heat .AND. elemen .AND. id(3) /= 82) GO TO 1400
 IF (dummy) GO TO 1640
 FMT(1) = oparen
 ifmt = 1
 j = ix + 1
 GO TO 1120
 1110 j = j + 1
 ifmt = ifmt + 1
 FMT(ifmt) = comma
 
!     IF K IS NEGATIVE THEN BUILDING BLOCK IS NOT FOR A VARIABLE.
!     IN THIS CASE THEN K IS ACTUAL POINTER TO BE USED IN THE ESINGL ARR
 
 1120 k = d(j)
 IF (k < 0) THEN
   GO TO  1130
 ELSE IF (k == 0) THEN
   GO TO  1300
 ELSE
   GO TO  1140
 END IF
 1130 k = -k
 ifmt = ifmt + 1
 FMT(ifmt) = esingl(k)
 GO TO 1110
 
!     CHECK FOR  SPECIAL PACKING FORMATS
 
 1140 k = 5*k - 5
 IF (.NOT. axic) GO TO 1150
 IF (k /= 200 .AND. k /= 275) GO TO 1160
 IF (iout(2) == iblank) GO TO 1230
 GO TO 1240
 1150 IF (.NOT.pack .OR. iout(2) == gse(1)) GO TO 1160
 IF ((i >= i1 .AND. i <= 8) .OR. (i >= i2 .AND. i <= 14)) GO TO 1250
 
!     IF SOLSET AND K=0 OR K=80 OR K=365 OR K=75 USE I15BLK IF INTEGER 1
 
 1160 IF (.NOT. solset) GO TO 1170
 IF (k /= 0 .AND. k /= 80 .AND. k /= 365 .AND. k /= 75) GO TO 1170
 IF (iout(i) /= 1) GO TO 1170
 iout(i) = iblank
 ifmt = ifmt + 2
 FMT(ifmt-1) = i15blk(1)
 FMT(ifmt  ) = i15blk(2)
 IF (axic) FMT(ifmt) = i2x
 GO TO 1260
 
!     CHECK FOR  0.0 ON AN E-FORMAT
 
 1170 IF (k < eend) IF (out(i)) 1230,1240,1230
 IF (k ==  440) IF (out(i)) 1230,1240,1230
 
!     CHECK FOR MID-EDGE OR CENTER STRESS POINTS ON ISOPARAMETRIC
!     SOLID ELEMENTS
 
 IF ((k /= 390 .AND. k /= 395) .OR. iout(i) /= 0) GO TO 1190
 IF (k == 395) GO TO 1180
 iout(i) = center
 GO TO 1240
 1180 iout(i) = iblank
 GO TO 1240
 
!     CHECK FOR SPECIAL INTEGER ON E9.1 FORMAT ONLY
 
 1190 IF (k /= e9pt1 .OR. iout(i) /= 1) GO TO 1200
 iout(i) = iblank
 GO TO 1240
 
!     CHECK FOR SPECIAL GPST FORMATS
 
 1200 IF (iout(i) /= 0 .OR. k < 301 .OR. k > 325) GO TO 1210
 iout(i) = iblank
 GO TO 1240
 
!     CHECK FOR HARMONIC NUMBER OR POINT ANGLE
 
 1210 IF (k /= 355 .AND. k /= 360 .AND. k /= 445) GO TO 1220
 IF (iout(i) <= 0 .OR. iout(i) >= 1000) GO TO 1220
 iout(i) = iout(i) - 1
 GO TO 1240
 
!     CHECK FOR PHASE ANGLE ON STRESSES WITH TRAPAX AND TRIAAX ELEMENTS
 
!     COMMENTS FROM G.CHAN/UNISYS   1/93
!     FMT AND PHASE ARE LOCAL. I SEE NOBODY SETTING UP FMT() TO PHASE.
!     PHASE IS '1PE9'. IN ANSI FORTRAN STANDARD, A COMMA IS NEEDED
!     BETWEEN P AND E IF PHASE IS REALLY USED IN SETTING UP FMT().
 
 1220 IF (k /= 430 .OR. id(9) /= 3 .OR. iapp /= freq .OR. FMT(ifmt-4)  &
      /= phase) GO TO 1230
 GO TO 1240
 
!     NO OTHER SPECIAL CHECKS AT THIS TIME
 
!     *** ADD FORMAT BLOCKS ***
 
!     STANDARD BLOCKS
 
 1230 ifmt = ifmt + 2
 FMT(ifmt-1) = e(k+1)
 FMT(ifmt  ) = e(k+2)
 GO TO 1260
 
!     ALTERNATE BLOCKS
 
 1240 ifmt = ifmt + 3
 FMT(ifmt-2) = e(k+3)
 FMT(ifmt-1) = e(k+4)
 FMT(ifmt  ) = e(k+5)
 GO TO 1260
 
!     SPECIAL BLOCKS FOR PACKED OUTPUT
 
 1250 ifmt = ifmt + 1
 FMT(ifmt) = a4
 GO TO 1260
 
 1260 i = i + 1
 GO TO 1110
 
!     OUTPUT THE LINE OR LINES OF DATA WITH THE NEW FORMAT
 
 1300 FMT(ifmt) = cparen
 IF (line > maxn .OR. .NOT.headng) CALL ofp1
 headng = .true.
 nwds = nwdsav
 
!     IF GRID-POINT-FORCE-BALANCE ENTRY, BLANK OUT NONEXISTENT (ZERO)
!     ELEMENT ID-S.
 
 IF (.NOT.gpfb) GO TO 1330
 IF (iout(2) == 0) THEN
   GO TO  1310
 ELSE
   GO TO  1320
 END IF
 1310 iout(2) = iblank
 FMT( 9) = a4
 FMT(10) = i9x
 
!     ALSO, FOR GPFB, SET FORMAT TO SPACE TWO LINES ON NEW POINT-ID.
 
 1320 IF (iout(1) == lastid) GO TO 1330
 lastid = iout(1)
 FMT(2) = ih0
 line   = line + 2
 1330 CONTINUE
!hgs revixe call to use new ofppnt: gfortran doesn't like variable format
!     ofppnt now calls forwrt to perform write
 CALL ofppnt (iout,nwds,FMT)
 GO TO 1700
 
!     ELEMENT FORCES IN HEAT PROBLEMS
 
 1400 IF (line > maxn .OR. .NOT.headng) CALL ofp1
 headng = .true.
 IF (sort == 2) GO TO 1520
 
!     BRANCH ON SPECIAL HBDY OUTPUT
 
 IF (nwds    == 5) GO TO 1460
 IF (iout(5) == 1) GO TO 1440
 IF (iout(6) == 1) GO TO 1420
 IF (iout(2) == hex1 .OR. iout(2) == hex2) GO TO 1480
 IF (iout(2) == hex3) GO TO 1500
!RPKR THE FOLLOWING LINE HAS BEEN CHANGED
!RPKR WRITE  (L,1410) IOUT(1),(OUT(I),I=2,9)
 WRITE  (l,1410) (iout(i),i=1,3), (out(i),i=4,9)
 1410 FORMAT (i14,4X,2A4,6(1P,e17.6))
 GO TO 1700
!RPKR THE FOLLOWING LINE HAS BEEN CHANGED
!1420 WRITE  (L,1430) IOUT(1),OUT(2),OUT(3),OUT(4),OUT(5),OUT(7),OUT(8)
 1420 WRITE  (l,1430) iout(1),iout(2),iout(3),out(4),out(5), out(7),out(8)
 1430 FORMAT (i14,4X,2A4,2(1P,e17.6),17X,2(1P,e17.6))
 GO TO 1700
!RPKR THE FOLLOWING LINE HAS BEEN CHANGED
!1440 WRITE  (L,1450) IOUT(1),OUT(2),OUT(3),OUT(4),OUT(7)
 1440 WRITE  (l,1450) iout(1),iout(2),iout(3),out(4),out(7)
 1450 FORMAT (i14,4X,2A4,1P,e17.6,34X,1P,e17.6)
 GO TO 1700
 
 1460 WRITE  (l,1470) iout(1),out(2),out(3),out(4),out(5)
 1470 FORMAT (18X,i18,4(1P,e18.6))
 GO TO 1700
 
 1480 WRITE  (l,1490) (iout(i),i=1,3),(out(i),i=4,9)
 1490 FORMAT (i14,1X,a4,i7,6(1P,e17.6))
 GO TO 1700
 1500 WRITE  (l,1510) (iout(i),i=1,3),(out(i),i=4,9)
 1510 FORMAT (i14,2X,a4,2X,a4,6(1P,e17.6))
 GO TO 1700
 
!     BRANCH ON SPECIAL HBDY FORCES
 
 1520 IF (nwds    == 5) GO TO 1580
 IF (iout(5) == 1) GO TO 1560
 IF (iout(6) == 1) GO TO 1540
 IF (iout(2) == hex1 .OR. iout(2) == hex2) GO TO 1600
 IF (iout(2) == hex3) GO TO 1620
!RPKR THE FOLLOWING LINE HAS BEEN CHANGED
!RPKR WRITE  (L,1530) (OUT(I),I=1,9)
 WRITE  (l,1530) out(1),iout(2),iout(3),(out(i),i=4,9)
 1530 FORMAT (1P,e14.6,4X,2A4,6(1P,e17.6))
 GO TO  1700
!RPKR THE FOLLOWING LINE HAS BEEN CHANGED
!1540 WRITE  (L,1550) (OUT(KK),KK=1,5),OUT(7),OUT(8)
 1540 WRITE  (l,1550) out(1),iout(2),iout(3),out(4),out(5), out(7),out(8)
 1550 FORMAT (1P,e14.6,4X,2A4,2(1P,e17.6),17X,2(1P,e17.6))
 GO TO  1700
!RPKR THE FOLLOWING LINE HAS BEEN CHANGED
!1560 WRITE  (L,1570) (OUT(KK),KK= 1,4),OUT(7)
 1560 WRITE  (l,1570) out(1),iout(2),iout(3),out(4),out(7)
 1570 FORMAT (1P,e14.6,4X,2A4,1P,e17.6,34X,1P,e17.6)
 GO TO  1700
 1580 WRITE  (l,1590) (out(kk),kk=1,5)
 1590 FORMAT (18X,5(1P,e18.6))
 GO TO  1700
!RPKR THE FOLLOWING LINE HAS BEEN CHANGED
!1600 WRITE  (L,1610) OUT(1),OUT(2),IOUT(3),(OUT(I),I=4,9)
 1600 WRITE  (l,1610) out(1),iout(2),iout(3),(out(i),i=4,9)
 1610 FORMAT (1P,e14.6,1X,a4,i7,6(1P,e17.6))
 GO TO  1700
!RPKR THE FOLLOWING LINE HAS BEEN CHANGED
!1620 WRITE  (L,1630) (OUT(I),I=1,9)
 1620 WRITE  (l,1630) out(1),iout(2),iout(3),(out(i),i=4,9)
 1630 FORMAT (1P,e14.6,2X,a4,2X,a4,6(1P,e17.6))
 GO TO  1700
 
!     DUMMY ELEMENT
 
 1640 IF (line > maxn .OR. .NOT.headng) CALL odum (2,l,itype,iapp,0,id)
 headng = .true.
 nwds   = nwdsav
 CALL odum (3,l,itype,iapp,nwds,iout)
 GO TO 1700
 
 1700 line = line + nlines
 IF (eor ) GO TO 60
 IF (axic) GO TO 900
 IF (.NOT.pack .OR. iout(2) == gse(1) .OR. i1 == 9) GO TO 900
 
!     TRANSFER THE SAVED BLOCK
 
 DO  i = 1,nwds
   iout(i) = isave(i)
 END DO
 IF (.NOT.fluid) GO TO 1040
 IF (iout(1) >= 500000) GO TO 1800
 GO TO 1040
 
!     SPECIAL ROUTINE TO PACK SCALAR OR EXTRA POINTS OUTPUT..
!     PACKING IS PERFORMED ONLY WHEN IDS ARE SEQUENTIAL,
!     AND THE TYPE REMAINS THE SAME.
 
 1720 i = 1
 grdpt = iout(1)
 TYPE  = iout(2)
 1730 IF (i == 6) GO TO 1760
 1740 from = 1740
 CALL READ (*2020,*1790,FILE,isave(1),nwds,0,flag)
 igse = isave(2)
 IF (pack) isave(2) = gse(igse)
 IF (sort == 2) GO TO 1750
 itemp  = isave(1)/10
 device = isave(1) - 10*itemp
 device = MOD(device,8)
 isave(1) = itemp
 1750 IF (device >= 4) CALL ofppun (isave(1),isave(1),nwds,iopt,idd, pnched)
 IF (device /= 1 .AND. device /= 3 .AND. device /= 5 .AND.  &
     device /= 7) GO TO 1740
 j = grdpt + i
 IF (fluid .AND. isave(1) >= 500000) GO TO 1760
 IF (isave(2) /= TYPE .OR. isave(1) /= (grdpt+i)) GO TO 1760
 
!     PACK THIS VECTOR INTO LINE OF DATA
!     IF COMPLEX TWO LINES OF DATA
!     IMAGINARY PART WILL BE PACKED EVEN IF IT DOES NOT EXIST.
 
 i = i + 1
 iout(i+2) = isave(3)
 iout(i+8) = isave(9)
 GO TO 1730
 
!     PUT BLANKS IN ANY OPEN SLOTS
 
 1760 j = i + 3
 IF (j > 8) GO TO 1780
 DO  k  = j,8
   iout(k  ) = iblank
   iout(k+6) = iblank
 END DO
 
 1780 i1 = j
 i2 = j + 6
 GO TO 1100
 
 1790 eor = .true.
 GO TO 1760
 
!     SPECIAL LOGIC FOR SORT-1 VECTOR OUTPUT IN A FLUID PROBLEM FOR
!     HARMONIC POINTS ONLY
 
 1800 oldhrm = -1
 l5  = 230
 line= maxn + 1
 k   = 0
 eor = .false.
 GO TO 1820
 1810 from = 1810
 CALL READ (*2020,*1840,FILE,iout(1),nwds,0,flag)
 
!     PUNCH PROCESSING
 
 itemp  = iout(1)/10
 device = iout(1) - 10*itemp
 iout(1)= itemp
 IF (device < 4) GO TO 1820
 device = MOD(device,8)
 CALL ofppun (iout(1), iout(1),incr,iopt,idd,pnched)
 device = device - 4
 IF (device <= 0) GO TO 1810
 
!     DECODE THE HARMONIC
 
 1820 itemp = MOD(iout(1),500000)
 nharm = (iout(1)-itemp)/500000
 iout(1) = itemp
 IF (oldhrm == -1) oldhrm = nharm
 IF (nharm /= oldhrm .OR. k >= 5) GO TO 1850
 k = k + 1
 1830 REAL(2*k-1) = iout(1)
 REAL(:: 2*k  ) = iout(3)
 imag(  k  ) = iout(9)
 GO TO 1810
 
!     OUTPUT THE LINE OF DATA
 
 1840 eor = .true.
 IF (k <= 0) GO TO 60
 
!     BUILD THE FORMAT
 
 1850 FMT(1) = i1x
 IF (nlines > 1) FMT(1) = i1h0
 FMT(2) = i12
 FMT(3) = i2xx
 ifmt   = 3
 
!     ADD STAR IF THIS IS AN UN-SYMETRIC HARMONIC
 
 IF (MOD(oldhrm,2) == 0) GO TO 1860
 FMT(3) = istar
 FMT(4) = i1xx
 ifmt   = 4
 
!     VARIABLES IN MAIN LINE
 
 1860 DO  i = 1,k
   FMT(ifmt+1) = i8
   IF (freal(2*i) == 0.0) THEN
     GO TO  1880
   END IF
   1870 FMT(ifmt+2) = pe
   FMT(ifmt+3) = e156
   ifmt = ifmt + 3
   CYCLE
   1880 FMT(ifmt+2) = pf
   FMT(ifmt+3) = f156
   FMT(ifmt+4) = i9x
   ifmt = ifmt + 4
 END DO
 
!     VARIABLES IN SECOND LINE IF COMPLEX
 
 IF (nlines <= 1) GO TO 1940
 ifmt = ifmt + 1
 FMT(ifmt) = i15x
 DO  i = 1,k
   ifmt = ifmt + 1
   FMT(ifmt) = comma
   IF (fimag(i) == 0.0) THEN
     GO TO  1900
   ELSE
     GO TO  1910
   END IF
   1900 FMT(ifmt+1) = pf
   FMT(ifmt+2) = f236
   FMT(ifmt+3) = i9x
   ifmt = ifmt + 3
   CYCLE
   1910 IF (l3 == 126) GO TO 1920
   FMT(ifmt+1) = pe
   FMT(ifmt+2) = e236
   ifmt = ifmt + 2
   CYCLE
   1920 FMT(ifmt+1) = pf
   FMT(ifmt+2) = f174
   FMT(ifmt+3) = i6x
   ifmt = ifmt + 3
 END DO
 
!     COMPLETE FORMAT
 
 1940 FMT(ifmt+1) = cparen
 IF (line > maxn) CALL ofp1
 line = line + nlines
 k2   = 2*k
 iharm= (oldhrm-1)/2
 IF (nlines <= 1) GO TO  1950
!RPKR THE FOLLOWING LINE HAS BEEN REPLACED BY TWO LINES
!RPKR WRITE (L,FMT) IHARM,(FREAL(I),I=1,K2),(FIMAG(I),I=1,K)
 WRITE (l,FMT) iharm,(REAL(i),freal(i+1),i=1,k2,2), (fimag(i),i=1,k)
 GO TO 1960
!RPKR THE FOLLOWING LINE HAS BEEN REPLACED
!1950 WRITE (L,FMT) IHARM,(FREAL(I),I=1,K2)
 1950 WRITE (l,FMT) iharm,(REAL(i),freal(i+1),i=1,k2,2)
 1960 k = 1
 oldhrm = nharm
 IF (.NOT. eor) GO TO 1830
 GO TO 60
 
!     ERROR CONDITION THIS FILE
 
 2000 WRITE  (l,2010) swm,ix,point,from
 2010 FORMAT (a27,', OFP BLOCK DATA ROUTINES UNAVAILABLE FOR THIS ',  &
     'ELEMENT.',11X,'IX,POINT,FROM =,',3I5)
 2020 WRITE  (l,2030) uwm,from
 2030 FORMAT (a25,' 3030, OFP UNABLE TO PROCESS DATA BLOCK.  A TABLE ',  &
     'PRINT OF THE DATA BLOCK FOLLOWS.   FROM =',i5,'/OFP')
 CALL CLOSE  (FILE,1)
 CALL tabprt (FILE)
 GO TO 2050
 
!     CLOSE FILE UP
 
 2040 CALL CLOSE (FILE,1)
 2050 IF (ifile == 6) GO TO 2060
 GO TO 50
 
!     RESTORE TITLES TO WHATEVER THEY WERE AT ENTRY TO OFP
 
 2060 DO  i = 1,96
   head(i) = tsave(i)
 END DO
 RETURN
END SUBROUTINE ofp
