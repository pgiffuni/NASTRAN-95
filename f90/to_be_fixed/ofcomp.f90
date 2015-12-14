SUBROUTINE ofcomp (*,FILE,TYPE,eltyp,iapp,headng,pnched,FORM)
     
!     OFP ROUTINE TO HANDLE PRINT AND PUNCH OF LAYERED COMPOSITE
!     ELEMENT STRESSES AND FORCES.  CURRENTLY, THIS INVOLVES ONLY
!     THE CQUAD4 AND CTRIA3 ELEMENTS.
 
!     FILE     = OUTPUT FILE UNDER PROCESSING
!     TYPE     = TYPE OF DATA-  REAL   , SORT 1       = 1
!                               COMPLEX, SORT 1       = 2
!                               REAL   , SORT 2       = 3
!                               COMPLEX, SORT 2       = 4
!     ELTYP    = ELEMENT TYPE-  QUAD4                 = 64
!                               TRIA3                 = 83
!     IAPP     = SOLUTION TYPE
!     HEADNG   = INDICATES PRINT HEADINGS ARE DONE FOR A PAGE
!     PNCHED   = INDICATES PUNCH HEADINGS ARE DONE
!     FORM     = DATA TYPE-     STRESSES              = 22
!                               FORCES                = 23
!                               STRAIN                = 21
 
 
 , INTENT(IN OUT)                         :: *
 INTEGER, INTENT(IN OUT)                  :: FILE
 INTEGER, INTENT(IN OUT)                  :: TYPE
 INTEGER, INTENT(IN OUT)                  :: eltyp
 INTEGER, INTENT(IN OUT)                  :: iapp
 LOGICAL, INTENT(OUT)                     :: headng
 LOGICAL, INTENT(OUT)                     :: pnched
 INTEGER, INTENT(IN OUT)                  :: FORM
 EXTERNAL        andf
 LOGICAL :: heat, cmpxdt,sort1,sort2, magpha, quad4,tria3,stress,force,strn
 INTEGER :: ist(86), flag,nout,punch,buf(86),ibuf(3),  &
     device,andf,head, static,freq,  &
     ceig,ititle(32),isubtl(32),label(32),elemid,  &
     failth,hill(2),hoffmn(2),tsaiwu(2),stresf(2),  &
     strain(2),ifail(2),blnk,astr,subst(3), id(50),of(58)
!     INTEGER         REIG,TRANS,BK1,ELEC
 REAL :: rst(86),rid(50),bufr(86),rbuf(3)
!     REAL            HARMON,PANGLE,BUFF(1)
 CHARACTER (LEN=5) :: t3q4,t3,q4
 COMMON /BLANK / icard
!     COMMON /ZZOFPX/ L1,L2,L3,L4,L5,ID(50),HARMON,PANGLE,BUFF(1)
 COMMON /zzzzzz/ core(1)
 COMMON /output/ head(96)
 COMMON /system/ ksystm(100)
 EQUIVALENCE     (ist(1)    ,rst(1) ), (id(1)     ,rid(1)  ),  &
     (buf(1)    ,bufr(1)), (ibuf(1)   ,rbuf(1) ),  &
     (ifail(1)  ,failmx ), (ifail(2)  ,maxflg  ),  &
     (ksystm(2) ,nout   ), (ksystm(9) ,maxlns  ),  &
     (ksystm(12),line   ), (ksystm(33),iflg    ),  &
     (ksystm(56),itherm ), (ksystm(69),isubs   ),  &
     (ksystm(91),punch  ), (head( 1) ,ititle(1)),  &
     (head(65) ,label(1)), (head(33) ,isubtl(1)),  &
     (l1, of(1) ,core(1)), (l2,of(2)),(l3,of(3)),  &
     (id(1)     ,of  (6)), (l4,of(4)),(l5,of(5))
!     EQUIVALENCE     (HARMON    ,OF (56)), (PANGLE   ,OF   (57)),
!    1                (BUFF(1)   ,OF (58))
 
 DATA  static,freq,ceig    / 1 , 5 , 9      /
!     DATA  REIG,TRANS,BK1,ELEC / 2 , 6 , 8 , 11 /
 DATA  hill  ,       hoffmn,       tsaiwu,       stresf       /  &
     4H   h,4HILL ,4HHOFF,4HMAN ,4HTSAI,4H-wu ,4H str,4HESS /
 DATA  strain       /4H str,4HAIN /
 DATA  blnk  ,astr  /4H    ,4H  * /
 DATA  subst        /4HSUBS,4HTRUC,4HTURE /
 DATA  t3q4, t3, q4 /' ', 'TRIA3', 'QUAD4'/
 
!     INITIALIZE
 
 cmpxdt = TYPE == 2 .OR. TYPE == 4
 sort1  = TYPE  <= 2
 sort2  = TYPE  > 2
 heat   = itherm == 1
 magpha = id(9) == 3 .AND. (iapp == freq .OR. iapp == ceig)
 quad4  = eltyp == 64
 tria3  = eltyp == 83
 stress = FORM  == 22
 force  = FORM  == 23
 strn   = FORM  == 21
 IF (heat .OR. sort2 .OR. cmpxdt) GO TO 1800
 
!     GET THE DEVICE CODE IF SORT=2,  1=PRINT  2=POST  4=PUNCH
 
 IF (sort1) GO TO 10
 idd = id(5)/10
 device = id(5) - 10*idd
 idevce = device
 id(5)  = idd
 elemid = idd
 10 CONTINUE
 
!     GET THE NUMBER OF OUTPUT WORDS PER ELEMENT.
 
 nwds = id(10)
 IF (nwds == 0) GO TO 1800
 IF (force) GO TO 40
 
!     ********************
!     ******* READ *******
!     ********************
 
 20 CALL READ (*1910,*1800,FILE,ist(1),3,0,flag)
 IF (sort1) elemid = ist(1)
 IF (sort2) time = rst(1)
 nlayer = ist(2)
 failth = ist(3)
 iply = 0
 30 iply = iply + 1
 IF (iply > nlayer) GO TO 20
 
 40 CALL READ (*1910,*1900,FILE,ist(1),nwds,0,flag)
 IF (stress .AND. iply == nlayer) CALL READ (*1910,*1910,FILE,ifail,2,0,flag)
 IF (force) elemid = ist(1)
 
!     GET THE DEVICE CODE IF SORT=1,   1=PRINT  2=POST  4=PUNCH
 
 IF (sort2) GO TO 100
 IF (stress .AND. iply > 1) GO TO 100
 itemp  = elemid / 10
 device = elemid - 10*itemp
 idevce = device
 elemid = itemp
 
!     *********************
!     ******* PUNCH *******
!     *********************
 
 100 IF (device < 4) GO TO 820
 
!     TAKE OUT INDEX FAILURE FLAGS FOR STRESSES
 
 numwds = nwds
 IF (stress) numwds = numwds - 2
 DO  ii=1,nwds
   buf(ii) = ist(ii)
 END DO
 IF (force) GO TO 120
 buf(6) = buf(7)
 buf(7) = buf(8)
 buf(8) = buf(9)
 120 CONTINUE
 
 IF (pnched) GO TO 500
 
!     PUNCH HEADINGS - TITLE, SUBTITLE, AND LABEL
 
 icard = icard + 1
 WRITE  (punch,130) (ititle(j),j=1,15),icard
 icard = icard + 1
 WRITE  (punch,140) (isubtl(j),j=1,15),icard
 icard = icard + 1
 WRITE  (punch,150) ( label(j),j=1,15),icard
 130 FORMAT (10H$title   =,15A4,2X,i8)
 140 FORMAT (10H$subtitle=,15A4,2X,i8)
 150 FORMAT (10H$label   =,15A4,2X,i8)
 
!     IF SUBSTRUCTURE (PHASE2) EXTRACTED ALSO SUBS-NAME AND COMPONENT
 
 IF (isubs == 0) GO TO 170
 IF (isubtl(20) /= subst(1) .OR. isubtl(21) /= subst(2) .OR.  &
     isubtl(22) /= subst(3)) GO TO 170
 icard = icard + 1
 WRITE (punch,160) (isubtl(j),j=20,26),icard
 icard = icard + 1
 WRITE  (punch,160) ( label(j),j=20,26),icard
 160 FORMAT (1H$,7A4,43X,i8)
 
 170 icard = icard + 1
 IF (stress) WRITE (punch,190) icard
 IF (force ) WRITE (punch,180) icard
 180 FORMAT (15H$element forces,57X,i8)
 190 FORMAT (17H$element stresses,55X,i8)
 
!     REAL, REAL/IMAGINARY, MAGNITUDE/PHASE
 
 icard = icard + 1
 IF (cmpxdt) GO TO 200
 WRITE  (punch,220) icard
 GO TO 250
 200 IF (magpha) GO TO 210
 WRITE  (punch,230) icard
 GO TO 250
 210 WRITE  (punch,240) icard
 220 FORMAT (12H$REAL output,60X,i8)
 230 FORMAT (22H$REAL-imaginary output,50X,i8)
 240 FORMAT (23H$magnitude-phase output,49X,i8)
 
!     SUBCASE OR ELEMENT ID
 
 250 icard = icard + 1
 IF (sort2) GO TO 260
 WRITE  (punch,280) id(4),icard
 GO TO 270
 260 WRITE  (punch,290) elemid,icard
 270 CONTINUE
 280 FORMAT (13H$subcase id =,i12,47X,i8)
 290 FORMAT (13H$element id =,i10,49X,i8)
 
!     PUNCH ELEMENT TYPE NUMBER,
!     IT IS SWITCHED TO MATCH THOSE OF POST PROCESSOR.
 
 icard  = icard + 1
 ieltyp = id(3)
 t3q4   = t3
 IF (ieltyp == 64) t3q4 = q4
 WRITE  (punch,300) ieltyp,t3q4,icard
 300 FORMAT (15H$element TYPE =,i12,4H   (,a5,1H),37X,i8)
 
!     EIGENVALUE, FREQUENCY, OR TIME
 
 SELECT CASE ( iapp )
   CASE (    1)
     GO TO 480
   CASE (    2)
     GO TO 400
   CASE (    3)
     GO TO 480
   CASE (    4)
     GO TO 480
   CASE (    5)
     GO TO 440
   CASE (    6)
     GO TO 450
   CASE (    7)
     GO TO 480
   CASE (    8)
     GO TO 400
   CASE (    9)
     GO TO 400
   CASE (   10)
     GO TO 480
   CASE (   11)
     GO TO 480
 END SELECT
 
!     PUNCH EIGENVALUE
 
 400 icard = icard + 1
 IF (sort1 .AND. cmpxdt) GO TO 410
 WRITE  (punch,420) rid(6),id(5),icard
 GO TO 480
 410 WRITE  (punch,430) rid(6),rid(7),id(5),icard
 GO TO 480
 420 FORMAT (13H$eigenvalue =,e15.7,2X,6HMODE =,i6,30X,i8)
 430 FORMAT (15H$eigenvalue = (,e15.7,1H,,e15.7,8H) mode =,i6,12X,i8)
 
!     FREQUENCY OR TIME
 
 440 IF (sort2) GO TO 480
 icard = icard + 1
 WRITE  (punch,460) rid(5),icard
 GO TO 480
 450 IF (sort2) GO TO 480
 icard = icard + 1
 WRITE  (punch,470) rid(5),icard
 460 FORMAT (12H$frequency =,e16.7,44X,i8)
 470 FORMAT (7H$time =,e16.7,49X,i8)
 
 480 pnched = .true.
 
!     PUNCH HEADINGS COMPLETE
 
 500 icard = icard + 1
 
!     ELEMENT STRESSES,  FIRST SUB-RECORD
 
 IF (force) GO TO 570
 IF (iply <= 1) GO TO 520
 WRITE  (punch,510) buf(1),bufr(2),bufr(3),icard
 510 FORMAT (6H-cont-,12X,i10,8X,2(1P,e18.6),i8)
 GO TO 560
 
 520 IF (sort2 .AND. iapp /= static) GO TO 540
 
!     FIRST CARD BEGINS WITH AN INTEGER
 
 WRITE  (punch,530) elemid,buf(1),bufr(2),bufr(3),icard
 530 FORMAT (i10,8X,i10,8X,2(1P,e18.6),i8)
 GO TO 560
 
!     FIRST CARD BEGINS WITH A REAL
 
 540 WRITE  (punch,550) time,buf(1),bufr(2),bufr(3),icard
 550 FORMAT (1P,e18.6,i10,8X,2(1P,e18.6),i8)
 560 nword = 3
 GO TO 620
 
!     ELEMENT FORCES,  FIRST SUB-RECORD
 
 570 IF (sort2 .AND. iapp /= static) GO TO 590
 
!     FIRST CARD BEGINS WITH AN INTEGER
 
 WRITE  (punch,580) buf(1),bufr(2),bufr(3),bufr(4),icard
 580 FORMAT (i10,8X,3(1P,e18.6),i8)
 GO TO 610
 
!     FIRST CARD BEGINS WITH A REAL
 
 590 WRITE  (punch,600) bufr(1),bufr(2),bufr(3),bufr(4),icard
 600 FORMAT (4(1P,e18.6),i8)
 610 nword = 4
 
 620 length = 8
 
!     SUBSEQUENT SUB-RECORDS
 
 700 left = numwds - nword
 IF (left > 0) GO TO 710
 IF (sort1) GO TO 810
 GO TO 820
 
!     PUNCH THE SUB-RECORDS
 
 710 IF (nword >= length) GO TO 700
 icard = icard + 1
 nword = nword + 3
 jout  = 3
 IF (nword <= length) GO TO  720
 nword = nword - 1
 jout  = 2
 IF (nword == length) GO TO  720
 nword = nword - 1
 jout  = 1
 
 720 jj = nword - jout + 1
 DO  ii = 1,jout
   ibuf(ii) = buf(jj)
   jj = jj + 1
 END DO
 SELECT CASE ( jout )
   CASE (    1)
     GO TO 740
   CASE (    2)
     GO TO 760
   CASE (    3)
     GO TO 780
 END SELECT
 
!     1 WORD OUT
 
 740 WRITE  (punch,750) rbuf(1),icard
 750 FORMAT (6H-cont-,12X,1P,e18.6,36X,i8)
 GO TO 800
 
!     2 WORDS OUT
 
 760 IF (iply < nlayer) WRITE (punch,770) rbuf(1),rbuf(2),icard
 IF (iply == nlayer) WRITE (punch,775) rbuf(1),rbuf(2),rbuf(3), icard
 770 FORMAT (6H-cont-,12X,1P,e18.6,0P,f18.4,18X,i8)
 775 FORMAT (6H-cont-,12X,1P,e18.6, 2(0P,f18.4),i8)
 GO TO 800
 
!     3 WORDS OUT
 
 780 WRITE  (punch,790) rbuf(1),rbuf(2),rbuf(3),icard
 790 FORMAT (6H-cont-,12X,1P,e18.6,0P,f18.4,1P,e18.6,i8)
 800 IF (jout < 3) GO TO 700
 GO TO 710
 
!     END OF PUNCH, SEE IF PRINT IS REQUESTED
 
 810 idevce = device - 4
 820 IF (andf(idevce,1) /= 0) GO TO 900
 IF (stress) GO TO 30
 GO TO 40
 
!     *********************
!     ******* PRINT *******
!     *********************
 
!     WRITE TITLES IF HAVE NOT DONE SO YET
 
 900 icheck = 0
 IF (line <= maxlns-2 .AND. headng) GO TO 910
 iflg = 1
 CALL page1
 headng = .true.
 icheck = 1
 
!     *** PRINT OF ELEMENT STRESSES ***
 
 910 IF (force) GO TO 1500
 
!     BRANCH ON TYPE OF OUTPUT
 
 SELECT CASE ( TYPE )
   CASE (    1)
     GO TO 920
   CASE (    2)
     GO TO 1400
   CASE (    3)
     GO TO 1410
   CASE (    4)
     GO TO 1420
 END SELECT
 
!     *** REAL, SORT 1 ***
 
 920 IF (icheck == 0) GO TO 1200
 SELECT CASE ( iapp )
   CASE (    1)
     GO TO 960
   CASE (    2)
     GO TO 930
   CASE (    3)
     GO TO 960
   CASE (    4)
     GO TO 960
   CASE (    5)
     GO TO 960
   CASE (    6)
     GO TO 940
   CASE (    7)
     GO TO 960
   CASE (    8)
     GO TO 950
   CASE (    9)
     GO TO 960
   CASE (   10)
     GO TO 960
   CASE (   11)
     GO TO 960
 END SELECT
 
 930 WRITE  (nout,970) id(5),rid(8),rid(6)
 GO TO 1010
 940 WRITE  (nout,980) rid(5)
 GO TO 1010
 950 WRITE  (nout,990) rid(6)
 GO TO 1010
 960 WRITE  (nout,1000)
 970 FORMAT (6X,'MODE NUMBER = ',i4,26X,'FREQUENCY = ',1P,e13.6,26X,  &
     'EIGENVALUE = ',1P,e13.6)
 980 FORMAT (6X,6HTIME =,1P,e14.6)
 990 FORMAT (6X,12HEIGENVALUE =,1P,e14.6)
 1000 FORMAT (1H )
 
 1010 CONTINUE
 IF (quad4) GO TO 1020
 IF (tria3) GO TO 1030
 GO TO 1050
 1020 WRITE  (nout,1070)
 GO TO 1050
 1030 WRITE  (nout,1080)
 GO TO 1050
 1050 WRITE  (nout,1100)
 WRITE  (nout,1110)
 1070 FORMAT (20X,'S T R E S S E S   I N   L A Y E R E D   ',  &
     'C O M P O S I T E   E L E M E N T S   ( Q U A D 4 )')
 1080 FORMAT (20X,'S T R E S S E S   I N   L A Y E R E D   ',  &
     'C O M P O S I T E   E L E M E N T S   ( T R I A 3 )')
 1100 FORMAT ('0 ELEMENT',3X,'PLY *STRESSES IN FIBER AND MATRIX',  &
     ' DIRECTIONS*  *DIRECT FIBER *  *INTER-LAMINAR STRESS',  &
     'ES*  * SHEAR BOND  *   *MAXIMUM*')
 1110 FORMAT (4X, 'ID', 6X, 'ID  *  NORMAL-1', 6X, 'NORMAL-2', 6X,  &
     'SHEAR-12 *  *FAILURE INDEX*  *SHEAR-1Z',6X,'SHEAR-2Z*',  &
     '  *FAILURE INDEX*   * INDEX *',/)
 
!     WRITE THE DATA
!     BUT FIRST, MODIFY THE FAILURE INDEX FLAGS FROM INTEGER TO BCD
 
 1200 IF (ist( 6) == 0) ist( 6) = blnk
 IF (ist( 6) == 1) ist( 6) = astr
 IF (ist(10) == 0) ist(10) = blnk
 IF (ist(10) == 1) ist(10) = astr
 
 IF (iply > 1) GO TO 1220
 WRITE  (nout,1210) elemid,ist(1),(rst(k),k=2,5),ist(6),  &
     (rst(k),k=7,9),ist(10)
 1210 FORMAT (1H0,i8,2X,i4,3(1P,e14.5),2X,0P,f10.3,a4,2(1P,e14.5),  &
     0P,f10.3,a4)
 nlines = 3
 GO TO 1730
 
 1220 WRITE  (nout,1230) ist(1),(rst(k),k=2,5),ist(6), (rst(k),k=7,9),ist(10)
 1230 FORMAT (11X,i4,3(1P,e14.5),2X,0P,f10.3,a4,2(1P,e14.5),0P,f10.3,a4)
 nlines = 1
 IF (iply < nlayer) GO TO 1730
 
!     IF THE LAST LAYER, CHECK THE MAXIMUM FAILURE INDEX
 
 nlines = 2
 IF (maxflg == 0) maxflg = blnk
 IF (maxflg == 1) maxflg = astr
 IF (failth /= 0) THEN
    SELECT CASE ( failth )
     CASE (    1)
       GO TO 1250
     CASE (    2)
       GO TO 1260
     CASE (    3)
       GO TO 1270
     CASE (    4)
       GO TO 1280
     CASE (    5)
       GO TO 1290
   END SELECT
 END IF
 failmx = 0.0
 WRITE  (nout,1240) failmx
 1240 FORMAT (1H ,116X,0P,f10.3)
 GO TO 1730
 1250 WRITE  (nout,1300) hill(1),hill(2),failmx,maxflg
 GO TO 1730
 1260 WRITE  (nout,1300) hoffmn(1),hoffmn(2),failmx,maxflg
 GO TO 1730
 1270 WRITE  (nout,1300) tsaiwu(1),tsaiwu(2),failmx,maxflg
 GO TO 1730
 1280 WRITE  (nout,1300) stresf(1),stresf(2),failmx,maxflg
 GO TO 1730
 1290 WRITE  (nout,1300) strain(1),strain(2),failmx,maxflg
 1300 FORMAT (1H ,41X,2A4,'FAILURE THEORY WAS USED FOR THIS ELEMENT.',  &
     26X,0P,f10.3,a4)
 GO TO 1730
 
!     *** COMPLEX, SORT 1 ***
 
 1400 GO TO 1800
 
!     *** REAL, SORT 2 ***
 
 1410 GO TO 1800
 
!     *** COMPLEX, SORT 2 ***
 
 1420 GO TO 1800
 
!     *** PRINT OF ELEMENT FORCES ***
 
 1500 CONTINUE
 
!     BRANCH ON TYPE OF OUTPUT
 
 SELECT CASE ( TYPE )
   CASE (    1)
     GO TO 1510
   CASE (    2)
     GO TO 1700
   CASE (    3)
     GO TO 1710
   CASE (    4)
     GO TO 1720
 END SELECT
 
!     *** REAL, SORT 1 ***
 
 1510 IF (icheck == 0) GO TO 1670
 SELECT CASE ( iapp )
   CASE (    1)
     GO TO 1550
   CASE (    2)
     GO TO 1520
   CASE (    3)
     GO TO 1550
   CASE (    4)
     GO TO 1550
   CASE (    5)
     GO TO 1550
   CASE (    6)
     GO TO 1530
   CASE (    7)
     GO TO 1550
   CASE (    8)
     GO TO 1540
   CASE (    9)
     GO TO 1550
   CASE (   10)
     GO TO 1550
   CASE (   11)
     GO TO 1550
 END SELECT
 
 1520 WRITE (nout,970) id(5),rid(8),rid(6)
 GO TO 1560
 1530 WRITE (nout,980) rid(5)
 GO TO 1560
 1540 WRITE (nout,990) rid(6)
 GO TO 1560
 1550 WRITE (nout,1000)
 
 1560 IF (quad4) GO TO 1570
 IF (tria3) GO TO 1580
 GO TO 1600
 1570 WRITE  (nout,1620)
 GO TO 1600
 1580 WRITE  (nout,1630)
 GO TO 1600
 1600 WRITE  (nout,1650)
 WRITE  (nout,1660)
 1620 FORMAT (22X,'F O R C E S   I N   L A Y E R E D   C O M P O S ',  &
     'I T E   E L E M E N T S   ( Q U A D 4 )'/)
 1630 FORMAT (22X,'F O R C E S   I N   L A Y E R E D   C O M P O S ',  &
     'I T E   E L E M E N T S   ( T R I A 3 )'/)
 1650 FORMAT (6X,'ELEMENT',18X,'- MEMBRANE  FORCES -',22X,'- BENDING',  &
     '   MOMENTS -',11X,'- TRANSVERSE SHEAR FORCES -')
 1660 FORMAT (8X,'ID',16X,2HFX,12X,2HFY,12X,3HFXY,11X,  &
     2HMX,12X,2HMY,12X,3HMXY,11X,2HVX,12X,2HVY)
 
!     WRITE THE DATA
 
 1670 WRITE  (nout,1680) elemid,(rst(k),k=2,9)
 1680 FORMAT (1H0,4X,i8,6X,8(1X,1P,e13.5))
 nlines = 2
 GO TO 1730
 
!     *** COMPLEX, SORT 1 ***
 
 1700 GO TO 1800
 
!     *** REAL, SORT 2 ***
 
 1710 GO TO 1800
 
!     *** COMPLEX, SORT 2 ***
 
 1720 GO TO 1800
 
!     DONE WITH ONE ENTRY, GO BACK AND READ ANOTHER ONE.
 
 1730 line = line + nlines
 IF (stress) GO TO 30
 GO TO 40
 
 1800 CONTINUE
 RETURN
 
 1900 IF (force) RETURN
 1910 CONTINUE
 RETURN 1
 
 
END SUBROUTINE ofcomp
