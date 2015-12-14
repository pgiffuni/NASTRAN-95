SUBROUTINE ofppun (ibuf,buf,nwds,iopt,idd,pnched)
     
!     MAIN OFP PUNCH ROUTINE FOR PUNCHING OF DATA LINES ONLY
 
!  $MIXED_FORMATS
 
 
 INTEGER, INTENT(IN)                      :: ibuf(nwds)
 REAL, INTENT(IN)                         :: buf(nwds)
 INTEGER, INTENT(IN OUT)                  :: nwds
 INTEGER, INTENT(IN OUT)                  :: iopt
 INTEGER, INTENT(IN)                      :: idd
 LOGICAL, INTENT(OUT)                     :: pnched
 LOGICAL :: temper
 INTEGER :: vector,id(50),of(56)
 REAL :: rid(50)
 COMMON /system/ sysbuf,l,dum53(53),itherm,dum34(34),lpch
 COMMON /output/ hd(96)
!     COMMON /ZZOFPX/ L1,L2,L3,L4,L5,ID(50)
 COMMON /zzzzzz/ core(1)
 COMMON /BLANK / icard
 COMMON /ofpcom/ temper, m
 COMMON /gpta1 / nelm,last,incr,ie(25,1)
 EQUIVALENCE     (rid(1),id(1),of(6)), (l1,of(1),core(1)),  &
     (l2,of(2)), (l3,of(3)), (l4,of(4)), (l5,of(5))
 DATA    vector, idtemp / 1, 0 /
 
 
 IF (.NOT. pnched) GO TO 700
 20 IF (nwds <  0) GO TO 1710
 
!     FIRST CARD OUT
 
 icard = icard + 1
 IF (iopt == vector) GO TO 200
 
!     GENERAL 1-ST CARD (FIRST WORD OF BUF ASSUMED INTEGER)
 
 n = MIN0(4,nwds)
 IF (idd < 0) THEN
   GO TO    30
 ELSE IF (idd == 0) THEN
   GO TO    90
 ELSE
   GO TO    40
 END IF
 30 IF (idd == -1) GO TO 90
 40 SELECT CASE ( n )
   CASE (    1)
     GO TO 50
   CASE (    2)
     GO TO 60
   CASE (    3)
     GO TO 70
   CASE (    4)
     GO TO 80
 END SELECT
 50 WRITE (lpch,440,ERR=180) buf(1),icard
 GO TO 180
 60 WRITE (lpch,450,ERR=180) buf(1),buf(2),icard
 GO TO 180
 70 WRITE (lpch,460,ERR=180) buf(1),buf(2),buf(3),icard
 GO TO 180
 80 WRITE (lpch,470,ERR=180) buf(1),buf(2),buf(3),buf(4),icard
 GO TO 180
 90 SELECT CASE ( n )
   CASE (    1)
     GO TO 100
   CASE (    2)
     GO TO 110
   CASE (    3)
     GO TO 120
   CASE (    4)
     GO TO 130
 END SELECT
 100 WRITE (lpch,400) ibuf(1),icard
 GO TO 180
 110 WRITE (lpch,410,ERR=180) ibuf(1),buf(2),icard
 GO TO 180
 120 WRITE (lpch,420,ERR=180) ibuf(1),buf(2),buf(3),icard
 GO TO 180
 
!     CHECK FOR THERMAL FORCES FOR ISOPARAMETRICS
 
 130 IF (itherm == 0 .OR. m /= 4) GO TO 150
 IF (id(3) < 65 .OR. id(3) > 67) GO TO 150
 WRITE  (lpch,140) ibuf(1),buf(2),ibuf(3),buf(4),icard
 140 FORMAT (i10,8X,a4,14X,i10,8X,1P,e18.6,i8)
 GO TO 180
 
!     CHECK FOR INTEGER IN SECOND ARGUMENT ALSO.
 
 150 IF (m == 19) GO TO 170
 IF (numtyp(buf(2)) <= 1) GO TO 160
 WRITE (lpch,430,ERR=180) ibuf(1),buf(2),buf(3),buf(4),icard
 GO TO 180
 160 WRITE (lpch,500,ERR=180) ibuf(1),ibuf(2),buf(3),buf(4),icard
 GO TO 180
 170 WRITE (lpch,510) ibuf(1),ibuf(2),buf(3),buf(4),icard
 GO TO 180
 180 nword = 4
 GO TO 230
 
!     VECTOR 1-ST CARD (FIRST WORD INTEGER, SECOND WORD BCD)
 
 200 IF (temper) GO TO 280
 IF (idd /= 0 .AND. idd /= -1) GO TO 210
 WRITE (lpch,520,ERR=220) ibuf(1),buf(2),buf(3),buf(4),buf(5),icard
 GO TO 220
 210 WRITE (lpch,530,ERR=220)  buf(1),buf(2),buf(3),buf(4),buf(5),icard
 220 nword = 5
 
!     CONTINUATION CARDS IF ANY.
 
 230 IF (nword >= nwds) GO TO 1710
 icard = icard + 1
 nword = nword + 3
 IF (nword <= nwds) GO TO 250
 nword = nword - 1
 IF (nword == nwds) GO TO 240
 nword = nword - 1
 
!     1 WORD OUT
 
 WRITE (lpch,610,ERR=1710) buf(nword),icard
 GO TO 1710
 
!     2 WORDS OUT
 
 240 WRITE (lpch,600,ERR=1710) buf(nword-1),buf(nword),icard
 GO TO 1710
 
!     3 WORDS OUT
 
 250 IF (ibuf(nword-1) == vector) GO TO 260
 IF (ibuf(nword  ) == vector) GO TO 270
 WRITE (lpch,590,ERR=230) buf(nword-2),buf(nword-1),buf(nword), icard
 GO TO 230
 260 WRITE (lpch,620) buf(nword-2),buf(nword),icard
 GO TO 230
 270 WRITE (lpch,600) buf(nword-2),buf(nword-1),icard
 GO TO 230
 
!     SPECIAL PUNCH ONLY WHEN TEMPER FLAG IS ON IN A -HEAT- FORMULATION.
 
 280 ic1 = ibuf(1)
 IF (idd == 0 .OR. idd == -1) GO TO 290
 idtemp = idtemp + 1
 ic1 = idd
 290 CONTINUE
 WRITE  (lpch,300) idtemp,ic1,buf(3),icard
 300 FORMAT (8HTEMP*    ,i16,i16,1P,e16.6,16X,i8)
 GO TO 1710
 
 400 FORMAT (i10,62X,i8)
 410 FORMAT (i10,8X,1P,e18.6,36X,i8)
 420 FORMAT (i10,8X,2(1P,e18.6),18X,i8)
 430 FORMAT (i10,8X,3(1P,e18.6),i8)
 440 FORMAT (1P,e18.6,54X,i8)
 450 FORMAT (2(1P,e18.6),36X,i8)
 460 FORMAT (3(1P,e18.6),18X,i8)
 470 FORMAT (4(1P,e18.6),i8)
 500 FORMAT (i10,8X,i10,8X,2(1P,e18.6),i8)
 510 FORMAT (i10,8X,i10,8X,2A4,28X,i8)
 520 FORMAT (i10,7X,a1,3(1P,e18.6),i8)
 530 FORMAT (1P,e16.6,1X,a1,3(1P,e18.6),i8)
 590 FORMAT (6H-cont-,12X,3(1P,e18.6),i8)
 600 FORMAT (6H-cont-,12X,2(1P,e18.6),18X,i8)
 610 FORMAT (6H-cont-,12X,1P,e18.6,36X,i8)
 620 FORMAT (6H-cont-,12X,1P,e18.6,18X,1P,e18.6,i8)
 
 
!     PUNCH HEADING CARDS
 
 
!     TITLE,SUBTITLE,AND LABEL
 
 700 DO  i = 1,3
   icard = icard + 1
   SELECT CASE ( i )
     CASE (    1)
       GO TO 710
     CASE (    2)
       GO TO 720
     CASE (    3)
       GO TO 730
   END SELECT
   710 WRITE (lpch,750) (hd(j),j= 1,15),icard
   CYCLE
   720 WRITE (lpch,760) (hd(j),j=33,47),icard
   CYCLE
   730 WRITE (lpch,770) (hd(j),j=65,79),icard
 END DO
 
 750 FORMAT (10H$title   =,15A4,2X,i8)
 760 FORMAT (10H$subtitle=,15A4,2X,i8)
 770 FORMAT (10H$label   =,15A4,2X,i8)
 
 ktype = id(2)/1000
 m = id(2) - (ktype)*1000
 IF (m < 1 .OR. m > 19) GO TO 1200
 icard = icard + 1
 SELECT CASE ( m )
   CASE (    1)
     GO TO 780
   CASE (    2)
     GO TO 790
   CASE (    3)
     GO TO 800
   CASE (    4)
     GO TO 810
   CASE (    5)
     GO TO 900
   CASE (    6)
     GO TO 1170
   CASE (    7)
     GO TO 910
   CASE (    8)
     GO TO 1170
   CASE (    9)
     GO TO 1170
   CASE (   10)
     GO TO 920
   CASE (   11)
     GO TO 930
   CASE (   12)
     GO TO 940
   CASE (   13)
     GO TO 1170
   CASE (   14)
     GO TO 950
   CASE (   15)
     GO TO 960
   CASE (   16)
     GO TO 970
   CASE (   17)
     GO TO 980
   CASE (   18)
     GO TO 990
   CASE (   19)
     GO TO 1000
 END SELECT
 780 WRITE (lpch,1010) icard
 GO TO 1200
 790 WRITE (lpch,1020) icard
 GO TO 1200
 800 WRITE (lpch,1030) icard
 GO TO 1200
 810 WRITE (lpch,1040) icard
 GO TO 1200
 
!     PUNCH ELEMENT STRESS OR GRID POINT STRESS HEADING LINE
 
 900 IF (l2 /= 378) WRITE(lpch,1050) icard
 IF (l2 == 378) WRITE(lpch,1060) icard
 GO TO 1200
 910 WRITE (lpch,1070) icard
 GO TO 1200
 920 WRITE (lpch,1080) icard
 GO TO 1200
 930 WRITE (lpch,1090) icard
 GO TO 1200
 940 WRITE (lpch,1100) icard
 GO TO 1200
 950 WRITE (lpch,1110) icard
 GO TO 1200
 960 WRITE (lpch,1120) icard
 GO TO 1200
 970 WRITE (lpch,1130) icard
 GO TO 1200
 980 WRITE (lpch,1140) icard
 GO TO 1200
 990 WRITE (lpch,1150) icard
 GO TO 1200
 1000 WRITE (lpch,1160) icard
 GO TO 1200
 
 1010 FORMAT (14H$displacements,58X,i8)
 1020 FORMAT (7H$oloads,65X,i8)
 1030 FORMAT (5H$spcf,67X,i8)
 1040 FORMAT (15H$element forces,57X,i8)
 1050 FORMAT (17H$element stresses,55X,i8)
 1060 FORMAT (24H$stresses at grid points,48X,i8)
 1070 FORMAT (12H$eigenvector,60X,i8)
 1080 FORMAT (9H$velocity,63X,i8)
 1090 FORMAT (13H$acceleration,59X,i8)
 1100 FORMAT (18H$non-linear-forces,54X,i8)
 1110 FORMAT (27H$eigenvector (solution set),45X,i8)
 1120 FORMAT (29H$displacements (solution set),43X,i8)
 1130 FORMAT (24H$velocity (solution set),48X,i8)
 1140 FORMAT (28H$acceleration (solution set),43X,i8)
 1150 FORMAT (23HELEMENT strain energies ,49X,i8)
 1160 FORMAT (24HGRID point force balance ,48X,i8)
 1170 icard = icard - 1
 
!     REAL, REAL/IMAGINARY, MAGNITUDE/PHASE
 
 1200 icard = icard + 1
 IF (ktype < 1 .OR. ktype == 2) GO TO 1210
 IF (id(9) == 3) GO TO 1230
 GO TO 1220
 1210 WRITE (lpch,1240) icard
 GO TO 1300
 1220 WRITE (lpch,1250) icard
 GO TO 1300
 
 1230 WRITE  (lpch,1260) icard
 1240 FORMAT (12H$REAL output,60X,i8)
 1250 FORMAT (22H$REAL-imaginary output, 50X,i8)
 1260 FORMAT (23H$magnitude-phase output,49X,i8)
 
!     SUBCASE NUMBER FOR SORT1 OUTPUT, OR
!     SUBCASE NUMBER FOR SORT2, FREQUENCY AND TRANSIENT RESPONSE ONLY
 
 1300 IF (ktype <= 1) GO TO 1310
 iapp = id(1)/10
 IF (iapp /= 5 .AND. iapp /= 6) GO TO 1400
 1310 icard = icard + 1
 WRITE  (lpch,1320) id(4),icard
 1320 FORMAT (13H$subcase id =,i12,47X,i8)
 
!     IF ELEMENT STRESS OR FORCE PUNCH ELEMENT TYPE NUMBER
 
 1400 IF (m /= 4 .AND. m /= 5) GO TO 1500
 icard = icard + 1
 id3   = id(3)
 IF (l2 /= 378) WRITE (lpch,1410) id3,ie(1,id3),ie(2,id3),icard
 IF (l2 == 378) WRITE (lpch,1420) icard
 1410 FORMAT (15H$element TYPE =,i12,3X,1H(,2A4,1H),32X,i8)
 1420 FORMAT (38H$punched in material coordinate system,34X,i8)
 
!     PUNCH EIGENVALUE, FREQUENCY, POINT OR ELEMENT ID, OR TIME
 
 1500 iapp = id(1)/10
 IF (iapp < 1 .OR. iapp > 10) GO TO 1700
 SELECT CASE ( iapp )
   CASE (    1)
     GO TO 1590
   CASE (    2)
     GO TO 1510
   CASE (    3)
     GO TO 1590
   CASE (    4)
     GO TO 1590
   CASE (    5)
     GO TO 1550
   CASE (    6)
     GO TO 1570
   CASE (    7)
     GO TO 1590
   CASE (    8)
     GO TO 1510
   CASE (    9)
     GO TO 1510
   CASE (   10)
     GO TO 1590
 END SELECT
 
!     PUNCH EIGENVALUE
 
 1510 icard = icard + 1
 IF (ktype == 1) GO TO 1530
 WRITE  (lpch,1520,ERR=1700) rid(6),id(5),icard
 1520 FORMAT (13H$eigenvalue =,e15.7,2X,6HMODE =,i6,30X,i8)
 GO TO 1700
 1530 WRITE  (lpch,1540,ERR=1700) rid(6),rid(7),id(5),icard
 1540 FORMAT (15H$eigenvalue = (,e15.7,1H,,e15.7,8H) mode =,i6,12X,i8)
 GO TO 1700
 
!     FREQUENCY OR TIME, POINT OR ELEMENT ID
 
 1550 IF (ktype > 1) GO TO 1590
 icard = icard + 1
 WRITE  (lpch,1560,ERR=1700) rid(5),icard
 1560 FORMAT (12H$frequency =,e16.7,44X,i8)
 GO TO 1700
 1570 IF (ktype > 1) GO TO 1590
 icard = icard + 1
 WRITE  (lpch,1580,ERR=1700) rid(5),icard
 1580 FORMAT (7H$time =,e16.7,49X,i8)
 GO TO 1700
 1590 IF (ktype <= 1) GO TO 1700
 icard = icard + 1
 IF (m == 4 .OR. m == 5) GO TO 1610
 WRITE  (lpch,1600) id(5),icard
 1600 FORMAT (11H$point id =,i12,49X,i8)
 GO TO 1700
 1610 WRITE  (lpch,1620) id(5),icard
 1620 FORMAT (13H$element id =,i10,49X,i8)
 
!     CARD HEADING COMPLETE
 
 1700 pnched = .true.
 IF (.NOT.temper) GO TO 20
 idtemp = idtemp + 1
 IF (idd > 0) idtemp = 0
 GO TO 20
 
 1710 RETURN
END SUBROUTINE ofppun
