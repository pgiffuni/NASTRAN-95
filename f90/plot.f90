SUBROUTINE plot (mode,buf1,b1,setid,deflst,nofind)
     
!     THIS PLOT ROUTINE IS CALLED ONLY BY PARAM
 
 
 INTEGER, INTENT(IN OUT)                  :: mode
 INTEGER, INTENT(IN OUT)                  :: buf1
 INTEGER, INTENT(IN)                      :: b1
 INTEGER, INTENT(IN)                      :: setid(1)
 INTEGER, INTENT(OUT)                     :: deflst(2)
 INTEGER, INTENT(OUT)                     :: nofind
 EXTERNAL          andf
 LOGICAL :: tapbit   ,stress   ,disp
 INTEGER :: andf     ,anydef   ,awrd(2)  ,bfrms    ,  &
     bufsiz   ,casecc   ,d1       ,d2       ,defbuf   ,defile ,  &
     defid    ,DIRECT   ,dtype    , eor      ,elset  ,  &
     ERR(17)  ,gpset    ,oes1     ,origin   ,for      ,parm   ,  &
     pbufsz   ,pcon     ,pedge    ,plabel   ,pltbuf   ,pltnum ,  &
     plttap   ,plttyp   ,porig    ,ppen     ,prnt     ,prject ,  &
     pset     ,psymm    ,pshape   ,psymbl   ,pvectr   ,rew    ,  &
      setd     ,skpttl   ,stereo   ,subc(2)  ,subcas ,  &
     tra      ,where    ,word     ,NAME(2)  ,skplod   ,thlid  ,  &
     fscale   ,fvp      ,offscl   ,org
 INTEGER :: nf(2)    ,f1(10)   ,f2(20)   ,msg1(19) ,msg2(17) ,mf4(6) ,  &
     msg7(13) ,mf3(3,3) ,used(10)
 INTEGER :: all      ,both     ,contur   ,defo     ,elem     ,epid   ,  &
     keywd    ,grid     ,gspc     ,lag(2)   ,magc(2)  ,TO     ,  &
     mf1(2,5) ,mf2(2,5) ,poin     ,rang     ,rqst(17) ,thru   , time
 REAL :: frr(17)  ,maxdef
 DOUBLE  PRECISION dwrd
 CHARACTER (LEN=25) :: uwm
 CHARACTER (LEN=23) :: ufm
 COMMON /xmssg /   ufm      ,uwm
 COMMON / pltscr / ncor, pltsc(50)
 COMMON /BLANK /   ngp      ,skp11    ,nsets    ,prnt     ,pltnum ,  &
     ngpset   ,anydef   ,skp12(3) ,parm     ,gpset    ,elset  ,  &
     casecc   ,skp21(3) ,defile(3),merr     ,setd     ,skp31  , oes1
 COMMON /system/   bufsiz   ,nout     ,dummy(66),isubs
 COMMON /output/   title(32,3)
 COMMON /pltdat/   skpplt(20),skpa(10),plttap
 COMMON /xxparm/   pbufsz   ,camera   ,bfrms    ,skpcam(3),  &
     penpap(30),scale(4),defmax   ,view(15) ,vanpnt(8),prject ,  &
     projct   ,for      ,org      ,norg     ,origin(11),  &
     skppar(77),ncntr   ,cntr(50) ,icntvl   ,where    ,DIRECT ,  &
     subcas   ,flag     ,DATA     ,skp19(19),icolor   ,skp235 , offscl
 COMMON /drwdat/   pset     ,plabel   ,porig    ,ppen     ,pshape ,  &
     psymbl(2),psymm(6) ,pvectr   ,pcon     ,pedge    ,offlag
 COMMON /plothd/   iused
 EQUIVALENCE       (ERR(1),frr(1))  , (word,awrd(1),iwrd,fwrd,dwrd)
 EQUIVALENCE       (fscale,scale(3)), (fvp,vanpnt(1))
 EQUIVALENCE       (skp19(1), lasset )
 DATA eor , inprew,norew,rew,skpttl,skplod / 1000000,0,2,1,37,5 /
 DATA subc/ 4HSUBC, 4HASES/
 DATA NAME/ 4H  pl, 4HOT  /
 
 DATA nf  / 10, 20 /      ,  &
     f1  / 4H(49X, 4H,4HP, 4HLOT,, 4HI9,2, 4HX,16, 4HHUND, 4HEFOR,  &
     4HMED , 4HSHAP, 4HE)   /       ,  &
     f2  / 4H(10X, 4H,4HP, 4HLOT,, 4HI5,3, 4HX,2(, 4HA4,a, 4H3),i,  &
     4H6,10, 4HH - , 4HSUBC, 4HASE,, 4HI8,3, 4HH - , 4H,a4,,  &
     4H1P,e, 4H15.6, 4H,1X,, 4H6A4,, 4HE11., 4H3)  /
 
!     DATA FOR FORMAT F2 - ORDER CORRESPONDING TO DTYPE, +10=VEL,+20=ACC
 
 DATA mf1 / 4HSTAT, 2HIC  , 4HFREQ, 1H.   , 4HTRAN,2HS.  ,  &
     4HMODA, 1HL   , 4HCMOD, 2HAL  /       ,  &
     mf2 / 4HDEFO, 3HRM. , 4HVELO, 1H.   , 4HACCE,2HL.  ,  &
     4HSTRE, 2HSS  , 4HSTRA, 2HIN  /
 DATA imod, load  /         4HMODE, 4HLOAD/
 DATA mf3 / 4H- fr, 4HEQUE, 4HNCY , 4H- ei, 4HGENV,4HALUE,  &
     4H- ti, 2HME  , 1H    /       ,  &
     mf4 / 4H pha, 4HSE l, 4HAG  , 4H mag, 4HNITU,2HDE  /
 
 DATA nmsg1,nmsg2,  nmsg7 / 19, 17, 13/   ,  &
     msg1/ 4H(33X, 4H,26H, 4HAN u, 4HNREC, 4HOGNI, 4HZABL, 4HE op,  &
     4HTION, 4H (,2, 4HA4,3, 4H1H) , 4HWAS , 4HDETE, 4HCTED,  &
     4H on , 4HA -p, 4HLOT-, 4H car, 4HD)  /       ,  &
     msg2/ 4H(34X, 4H,21H, 4HA no, 4HN-ex, 4HISTE, 4HNT o, 4HRIGI,  &
     4HN,i7, 4H,31H, 4H  is, 4H spe, 4HCIFI, 4HED o, 4HN a ,  &
     4H-plo, 4HT- c, 4HARD)/       ,  &
     msg7/ 4H(33X, 4H,41H, 4H*** , 4HINCO, 4HMPLE, 4HTE p, 4HLOT ,  &
     4HDUE , 4HTO i, 4HNPUT, 4H OR , 4HFILE, 4H.)  /
 
!     SET OPTIONS - FOLLOWING THE SET REQUEST(S)
 
 DATA rqst/ 4HSET , 4HORIG, 4HSHAP, 4HSYMB, 4HLABE, 4HVECT, 4HDENS,  &
     4HPEN , 4HSYMM, 4HANTI, 4HMAXI, 4HOUTL, 4HHIDD, 4HSHRI,  &
     4HNOFI, 4HFILL, 4HOFFS/
 
 DATA used/ 4H(49X, 4H,6HO, 4HRIGI, 4HN,i7, 4H,19H, 4H  us, 4HED i,  &
     4HN th, 4HIS p, 4HLOT)/
 
!     THE FOLLOWING ARE POSSIBLE OPTIONS ON THE PLOT CARD
 
 DATA defo/ 4HDEFO/, lorig/ 0     /,  &
     all / 3HALL /, TO   / 2HTO  /, thru/ 4HTHRU/, rang / 4HRANG/,  &
     time/ 4HTIME/, both / 4HBOTH/, grid/ 4HGRID/, poin / 4HPOIN/,  &
     elem/ 4HELEM/, gspc / 4HGSPC/, lag / 4HPHAS , 4HLAG        /,  &
     magc/ 4HMAGN , 4HIT.        /, epid/ 4HEPID/,contur/ 4HCONT/
 
 ncntr  = 10
 icntvl = 1
 where  = 1
 lasset = 0
 DIRECT = 2
 ncor   = 50
 DO  i = 1, 50
   pltsc(i) = 0
   cntr(i)  = 0
 END DO
 pltbuf = b1 - pbufsz
 defbuf = pltbuf - bufsiz
 IF (defbuf <= 0) GO TO 1400
 v1     =-1.e+30
 v2     =+1.e+30
 ph     = 0.0
 mag    = 0
 pcon   = 0
 loadid = 0
 lpcon  = 0
 flag   = 0.0
 subcas = 0
 defid  = 0
 disp   =.false.
 stress =.false.
 twopi  = 8.0*ATAN(1.0)
 ndef   = 0
 nogo   = 0
 CALL rdmodx (parm,mode,word)
 
 10 IF (mode <= 0) CALL rdmode (*10,*20,*40,mode,word)
 20 CALL rdword (mode,word)
 
!     CHECK FOR A DEFORMATION TYPE
 
 DO  dtype = 1,5
   IF (word == mf1(1,dtype)) GO TO 50
 END DO
 40 dtype = 0
 IF (word /= contur .OR. mode >= eor) GO TO 180
 pcon  = 1
 plttyp= 1
 GO TO 90
 
!     DEFORMATION TYPE SPECIFIED. CHECK IF ALL ARE TO BE PLOTTED
 
 50 plttyp = 1
 IF (mode <= 0) CALL rdmode (*120,*60,*110,mode,word)
 60 CALL rdword (mode,word)
 DO  plttyp = 1,3
   IF (word == mf2(1,plttyp)) GO TO 80
 END DO
 plttyp = 1
 IF (word /= contur) GO TO 110
 pcon = 1
 GO TO 80
 
!     ACCEL, VELOCITY ONLY ALLOWED FOR TRANS OR FREQUENCY RESPONSE.
!     NOTE THAT A COMPLEX  IGENVALUE WOULD BE NEEDED FOR -CMODAL-
 
 80 IF ((dtype == 2 .OR. dtype == 3) .OR. plttyp == 1) GO TO 90
 ERR(1) = 2
 ERR(2) = awrd(1)
 ERR(3) = awrd(2)
 CALL wrtprt (merr,ERR,msg1,nmsg1)
 plttyp = 1
 90 IF (mode <= 0) CALL rdmode (*120,*100,*110,mode,word)
 100 CALL rdword (mode,word)
 110 ndef = 1
 deflst(1) = all
 GO TO 180
 
!     THE DEFORMATIONS MAY BE EXPLICITLY LISTED AND/OR A RANGE MAY BE
!     LISTED (I.E., N1,N2 AND/OR N1 -TO/THRU- N2)
 
 120 ASSIGN 130 TO tra
 GO TO 1450
 130 ndef = ndef + 1
 deflst(ndef) = iwrd
 CALL rdmode (*1450,*140,*170,mode,word)
 140 CALL rdword (mode,word)
 IF (mode /= 0 .OR. (word /= TO .AND. word /= thru)) GO TO 170
 ASSIGN 150 TO tra
 CALL rdmode (*1450,*160,*170,mode,word)
 150 deflst(ndef+1) = TO
 deflst(ndef+2) = iwrd
 ndef = ndef + 2
 CALL rdmode (*120,*160,*170,mode,word)
 160 CALL rdword (mode,word)
 170 IF (ndef /= 1 .OR. deflst(1) /= 0) GO TO 180
 ndef = 2
 deflst(2) = all
 
!     ALL THE LISTED DEFORMATION ID-S HAVE BEEN READ
 
 180 deflst(ndef+1) = 0
 IF (mode >= eor) GO TO 340
 
!     TEST FOR CONTOUR REQUEST
 
 190 IF (word /= contur) GO TO 240
 IF (pcon ==      0) GO TO 220
 200 ERR(2) = awrd(1)
 ERR(3) = awrd(2)
 210 ERR(1) = 2
 CALL wrtprt (merr,ERR,msg1,nmsg1)
 GO TO 320
 
 220 pcon = 1
 IF (dtype == 0) plttyp = 1
 IF (ndef  /= 1) GO TO 230
 ndef = 0
 GO TO 90
 230 IF (mode > 0) GO TO 320
 ERR(2) = subc(1)
 ERR(3) = subc(2)
 GO TO 210
 
!     TEST FOR RANGE / TIME  (UNITS=LAMDA,F, OR TIME)
 
 240 IF (word /= rang .AND. word /= time) GO TO 270
 IF (pcon == 0 .AND. dtype == 1) GO TO 200
 ASSIGN 250 TO tra
 IF (mode > 0) GO TO 200
 CALL rdmode (*1490,*330,*340,mode,word)
 250 v1 = fwrd
 ASSIGN 260 TO tra
 CALL rdmode (*1490,*330,*340,mode,word)
 260 v2 = fwrd
 GO TO 320
 
!     TEST FOR PHASE LAG (COMPLEX DATA)
 
 270 IF (word /= lag(1)) GO TO 310
 IF (dtype /= 2 .AND. dtype /= 5) GO TO 200
 ASSIGN 300 TO tra
 280 IF (mode <= 0) CALL rdmode (*1490,*290,*340,mode,word)
 290 CALL rdword (mode,word)
 IF (word == lag(2)) GO TO 280
 GO TO 340
 300 IF (mag == 0) ph = fwrd
 GO TO 320
 
!     TEST FOR MAGNITUDE (COMPLEX DATA)
 
 310 IF (word /= magc(1)) GO TO 340
 IF (dtype /= 2 .AND. dtype /= 5) GO TO 200
 IF (ph == 0.0) mag = 1
 GO TO 320
 
 320 IF (mode <= 0) CALL rdmode (*320,*330,*340,mode,word)
 330 CALL rdword (mode,word)
 GO TO 190
 
!     READ THE REST OF THE PLOT CARD INTO STORAGE - DEFLST(N1-N2)
 
 340 n1 = ndef + 1
 n2 = n1 + 1
 IF (mode < eor) GO TO 350
 deflst(n2) = mode
 n2 = n2 + 1
 GO TO 400
 350 n  = 0
 360 deflst(n2+1) = awrd(1)
 deflst(n2+2) = awrd(2)
 n2 = n2 + 2
 n  = n  + 1
 IF (mode == 0) GO TO 370
 CALL rdword (mode,word)
 GO TO 360
 370 n2 = n2 + 1
 deflst(n1+1) = n
 380 CALL READ (*1520,*390,parm,deflst(n2),defbuf-n2+1,0,n)
 GO TO 1400
 390 n2 = n2 + n
 
!     SAVE LENGTH OF OPEN CORE USED IN IUSED FOR HDPLOT
 
 iused = n2 + nsets
 IF (deflst(n2-1) == 0) GO TO 380
 400 n2 = n2 - 1
 
!     INITIATE THE PLOTS OF THE REQUESTED DEFORMATIONS.
 
 nplots = 0
 IF (prnt < 0) GO TO 420
 IF (dtype == 0 .AND. pcon == 0) GO TO 410
 anydef = 1
 GO TO 1430
 
!     DO THE UNDEFORMED PLOT
 
 410 defid  = 0
 defbuf = defbuf + bufsiz
 IF (isubs == 0 .AND. .NOT.tapbit(plttap)) GO TO 1520
 GO TO 700
 420 IF (dtype == 0 .AND. pcon == 0) GO TO 1430
 
!     DO THE DEFORMED PLOT
 
!     STRESS IS TRUE IF CONTOUR REQUEST IS FOR STRESS
 
 lpcon = pcon
 IF (.NOT.tapbit(plttap)) GO TO 1520
 IF (pcon /= 0 .AND. icntvl <= 9) stress = .true.
 IF (pcon /= 0 .AND. icntvl > 13) stress = .true.
 IF ((pcon /= 0 .AND. (icntvl > 9.AND.icntvl < 14)) .OR.  &
     dtype /= 0) disp = .true.
 IF (.NOT.disp) GO TO 470
 mdef = defile(1)
 IF (dtype > 1) mdef = defile(2)
 IF (dtype > 0) GO TO 460
 
!     USER SPECIFIED CONTOUR DISP AND NOT THE TYPE
!     USE FIRST NON-NULL FILE
 
 430 CALL OPEN (*440,mdef,deflst(defbuf),inprew)
 CALL skprec (mdef,1)
 GO TO 450
 440 IF (mdef == defile(2)) CALL mesage (-1,mdef,NAME)
 mdef = defile(2)
 GO TO 430
 
!     SET DTYPE BY MFILE
 
 450 CALL READ (*1390,*1390,mdef,ERR(1),2,0,i)
 mfile = MOD(ERR(2),10)
 dtype = mfile
 CALL CLOSE (mdef,rew)
 460 CONTINUE
 
!     CALCULATE HEADER WORD 2 NEEDED FOR PLOT FILE CHECK
 
 mfile = dtype
 IF (dtype == 3) mfile = 3 + (plttyp-1)*10
 
!     OPEN OES1 AND MDEF
 
 IF (.NOT.disp) GO TO 470
 CALL OPEN (*1430,mdef,deflst(defbuf),inprew)
 CALL skprec (mdef,1)
 470 IF (.NOT.stress) GO TO 500
 CALL OPEN (*1390,oes1,deflst(b1),inprew)
 CALL skprec (oes1,1)
 IF (.NOT.disp) plttyp = 4
 CALL fread  (oes1,i,1,0)
 CALL bckrec (oes1)
 i = i/10
 japp = i
 IF (dtype /= 0) GO TO 475
 IF (i == 1 .OR. i == 3 .OR. i == 4 .OR. i == 7 .OR. i == 10) dtype = 1
 IF (i == 2 .OR. i == 8) dtype = 4
 IF (i == 6) dtype = 3
 
!     FOR STRESS PLOTS SET -FLAG- SO FNDSET KNOWS WHICH WORD TO COMPARE
 
 475 IF (dtype == 1) GO TO 480
 IF (dtype > 1) flag = 1.0
 IF (dtype > 3) flag = 2.0
 480 IF (dtype == 0) GO TO 1410
 IF (.NOT.disp) defbuf = defbuf + bufsiz
 
!     READ THE PLOT TITLES FOR EACH DEFORMED SHAPE TO BE DRAWN
 
 500 pcon = lpcon
 IF (.NOT.disp) GO TO 540
 510 CALL READ  (*1385,*1385,mdef,defid,1,0,i)
 CALL fread (mdef,n,1,0)
 IF (n == mfile) GO TO 515
 CALL skprec (mdef,1)
 GO TO 530
 515 CONTINUE
 CALL fread (mdef,loadid,1,0)
 CALL fread (mdef,value, 1,1)
 IF (value < v1 .OR. value > v2) GO TO 530
 DATA   = value
 subcas = defid
 n = 1
 520 IF (deflst(n) == all) GO TO 540
 CALL intlst (deflst,n,i,d1,d2)
 IF (defid >= d1 .AND. defid <= d2) GO TO 540
 IF (n < n1) GO TO 520
 530 CALL skprec (mdef,1)
 GO TO 510
 
!     POSITION OES1 IF NEEDED
 
 540 IF (.NOT.stress) GO TO 660
 IF (nplots /= 0) CALL OPEN (*1390,oes1,deflst(b1),norew)
 550 CALL READ (*1385,*1385,oes1,iapp,1,0,i)
 
!     VERIFY OES1 IS FOR CURRENT DTYPE
 
 iapp = iapp/10
 IF (iapp /= japp) GO TO 1385
 CALL fread (oes1,0,-2,0)
 CALL fread (oes1,i,1,0)
 IF (.NOT.disp ) GO TO 570
 IF (i /= defid) GO TO 620
 570 subcas = i
 v = value
 CALL fread (oes1,ERR(1),4,0)
 IF (dtype == 1) GO TO 575
 IF (dtype >= 4) GO TO 580
 
!     TRANSIENT
 
 v = frr(1)
 
!     STATICS
 
 575 j = ERR(4)
 GO TO 590
 
!     MODAL
 
 580 j = ERR(1)
 v = frr(2)
 IF (dtype == 4 .AND. iapp == 2) v = SQRT(ABS(v))/twopi
 590 IF (.NOT.disp) GO TO 600
 
!     ACCOUNT FOR ROUNDOFF
 
 IF (ABS(v-value) > 1.0E-6) GO TO 620
 DATA = value
 GO TO 650
 600 IF (v < v1 .OR. v > v2) GO TO 620
 DATA = v
 n = 1
 610 IF (deflst(n) == all) GO TO 650
 CALL intlst (deflst,n,i,d1,d2)
 IF (subcas >= d1 .AND. subcas <= d2) GO TO 650
 IF (n < n1) GO TO 610
 
!     WRONG CASE
 
 620 CALL fwdrec (*1410,oes1)
 CALL fwdrec (*1410,oes1)
 GO TO 550
 
!     LOCATED CASE TO PLOT
 
 650 CALL bckrec (oes1)
 loadid = j
 defid  = subcas
 value  = DATA
 
 660 CALL gopen (casecc,deflst(buf1),inprew)
 670 CALL READ  (*690,*690,casecc,n,1,0,i)
 IF (n == defid) GO TO 675
 CALL fread (casecc,0,0,1)
 GO TO 670
 675 CALL fread (casecc,0,-skplod,0)
 CALL fread (casecc,thlid,1,0)
 IF (loadid == 0) loadid = thlid
 skpttl = 31
 CALL fread (casecc,0,-skpttl,0)
 CALL fread (casecc,title,3*32,0)
 CALL CLOSE (casecc,rew)
 GO TO 700
 690 CALL CLOSE (casecc,rew)
 IF (.NOT.disp) GO TO 550
 CALL fread (mdef,0,0,1)
 GO TO 510
 
!     IDENTIFY THE PLOT
 
 700 pltnum = pltnum + 1
 IF (stress) CALL CLOSE (oes1,norew)
 CALL sopen (*1430,plttap,deflst(pltbuf),pbufsz)
 ncntr = -IABS(ncntr)
 IF (nplots == 0) CALL pltopr
 nplots = nplots + 1
 stereo = 0
 mtyp   = 0
 ERR(2) = pltnum
 IF (.NOT.(disp .OR. stress)) GO TO 720
 ERR(3) = mf1(1,dtype)
 ERR(4) = mf1(2,dtype)
 IF (icntvl == 20) plttyp = 4
 ERR(5) = mf2(1,plttyp)
 ERR(6) = mf2(2,plttyp)
 ERR(7) = defid
 ERR(8) = loadid
 ERR(9) = load
 IF (dtype /= 1) GO TO 710
 ERR(1) = 8
 GO TO 730
 710 ERR(1) = 12
 IF (dtype > 3) ERR(9) = imod
 frr(10) = value
 mtyp = 1
 IF (dtype == 3) mtyp = 3
 IF (dtype == 4 .AND. loadid < 0) mtyp = 2
 IF (mtyp == 2) ERR(8) = -loadid
 ERR(11) = mf3(1,mtyp)
 ERR(12) = mf3(2,mtyp)
 ERR(13) = mf3(3,mtyp)
 IF (dtype == 3 .OR. dtype == 4) GO TO 730
 ERR(1) = 15
 m = 0
 IF (mag /= 0) m = 3
 ERR(14) = mf4(m+1)
 ERR(15) = mf4(m+2)
 ERR(16) = mf4(m+3)
 IF (mag /= 0) GO TO 730
 ERR(1)  = 16
 frr(17) = ph
 GO TO 730
 720 ERR(1) = 1
 CALL wrtprt (merr,ERR,f1,nf(1))
 GO TO 740
 730 CALL wrtprt (merr,ERR,f2,nf(2))
 740 CALL stplot (pltnum)
 CALL head (dtype,plttyp,mtyp,ERR)
 
!     PLOT EACH SET REQUESTED. INTERPRET THE ASSOCIATED REQUESTS.
 
 750 CALL rdmody (deflst(n1+1),mode,word)
 mode   = 0
 maxdef = 0.
 porig  = 1
 ppen   = 1
 pset   = 0
 760 plabel = -1
 pcon   = lpcon
 pshape = 1
 pvectr = 0
 offlag = 0
 pedge  = 0
 psymbl(1) = 0
 psymbl(2) = 0
 psymm(1) = 1
 psymm(2) = 1
 psymm(3) = 1
 psymm(4) = 1
 psymm(5) = 1
 psymm(6) = 1
 780 IF (mode <= 0) CALL rdmode (*780,*790,*1180,mode,word)
 790 CALL rdword (mode,word)
 
!     CHECK FOR THE KEYWORD. THIS MAY BE FOLLOWED BY QUALIFIERS
 
 800 CONTINUE
 DO  keywd = 1,17
   IF (word == rqst(keywd)) GO TO 804
 END DO
 GO TO 1170
 804 SELECT CASE ( keywd )
   CASE (    1)
     GO TO 1080
   CASE (    2)
     GO TO  910
   CASE (    3)
     GO TO  960
   CASE (    4)
     GO TO  990
   CASE (    5)
     GO TO  830
   CASE (    6)
     GO TO 1060
   CASE (    7)
     GO TO  810
   CASE (    8)
     GO TO  810
   CASE (    9)
     GO TO 1020
   CASE (   10)
     GO TO 1020
   CASE (   11)
     GO TO 880
   CASE (   12)
     GO TO 1140
   CASE (   13)
     GO TO 1148
   CASE (   14)
     GO TO 1142
   CASE (   15)
     GO TO 1175
   CASE (   16)
     GO TO  805
   CASE (   17)
     GO TO 1160
 END SELECT
 
!             SET ORIG SHAP SYMB LABE VECT DENS  PEN SYMM ANTI
!    1       MAXI OUTL HIDD SHRI NOFI FILL OFFS
 
!     FILL ELEMENTS BY SET HERE
!     FILL PRESENTLY DOES NOT WORK TOGETHER WITH SHRINK AND HIDDEN
 
 805 ppen  = ppen + 31
 pedge = 100
 GO TO 780
 
!     DENSITY I, PEN I
 
 810 IF (mode /= 0) GO TO 1170
 ASSIGN 820 TO tra
 GO TO 1440
 820 ppen = iwrd
 GO TO 780
 
!     LABEL GRID / ELEMENTS
 
 830 plabel = 0
 IF (mode <= 0) CALL rdmode (*780,*840,*1180,mode,word)
 840 CALL rdword (mode,word)
 IF (word == both) GO TO 870
 IF (word == elem) GO TO 860
 IF (word /= grid) GO TO 872
 IF (mode <= 0) CALL rdmode (*780,*850,*1180,mode,word)
 850 CALL rdword (mode,word)
 IF (word-poin == 0.0) THEN
   GO TO   780
 ELSE
   GO TO   800
 END IF
 860 plabel = 3
 GO TO 780
 870 plabel = 6
 GO TO 780
 872 IF (word == gspc) plabel = 1
 IF (word == epid) plabel = 4
 IF (plabel /=  0) GO TO 780
 GO TO 800
 
!     MAXIMUM DEFORMATION X.X
 
 880 CONTINUE
 ASSIGN 900 TO tra
 IF (mode <= 0) CALL rdmode (*1490,*890,*1180,mode,word)
 890 CALL rdword (mode,word)
 IF (word /= defo .OR. mode /= 0) GO TO 800
 GO TO 1480
 900 maxdef = ABS(fwrd)
 GO TO 780
 
!     ORIGIN I
 
 910 IF (mode /= 0) GO TO 1170
 ASSIGN 920 TO tra
 GO TO 1440
 920 DO  i = 1,org
   IF (origin(i) == iwrd) GO TO 940
 END DO
 IF (stereo /= 0) GO TO 780
 ERR(1) = 1
 ERR(2) = iwrd
 CALL wrtprt (merr,ERR,msg2,nmsg2)
 GO TO 780
 940 porig  = i
 GO TO 780
 
!     SHAPE
 
 960 IF (pedge /= 0) GO TO 1170
 IF ((.NOT.(disp .OR. stress) .AND. dtype /= 0)) GO TO 1170
 IF (.NOT.disp) GO TO 780
 pshape = 2
 DO  i = 1,ndef
   IF (deflst(i) == 0) GO TO 980
 END DO
 GO TO 780
 980 pshape = 3
 GO TO 780
 
!     SYMBOL I,I
 
 990 psymbl(1) = 1
 IF (mode /= 0) GO TO 1170
 ASSIGN 1010 TO tra
 i = 0
 1000 i = i + 1
 GO TO 1440
 1010 psymbl(i) = iwrd
 IF (i-2 < 0) THEN
   GO TO  1000
 ELSE
   GO TO   780
 END IF
 
!     SYMMETRY B / ANTISYMMETRY B
 
 1020 n = 1
 IF (keywd == 10) n = -1
 IF (mode  <=  0) GO TO 1170
 CALL rdword (mode,word)
 CALL intvec (word)
 IF (word < 1 .OR. word > 7) GO TO 1170
 DO  i = 1,3
   psymm(i) = 1
   IF (andf(word,2**(i-1)) /= 0) psymm(i) = -1
   psymm(i+3) = n*psymm(i)
 END DO
 GO TO 780
 
!     VECTOR B
 
 1060 IF (.NOT.disp .OR. mode == 0) GO TO 1170
 CALL rdword (mode,word)
 pvectr = word
 GO TO 780
 
!     SET -  SAVE FIRST ENCOUNTERED, DO PLOT WHEN EOR OR ANOTHER SET
 
 1080 IF (mode /= 0) GO TO 1170
 ASSIGN 1090 TO tra
 GO TO 1440
 1090 iwrd = IABS(iwrd)
 DO  i = setd,nsets
   IF (iwrd == setid(i)) GO TO 1120
 END DO
 IF (stereo /= 0) GO TO 1110
 WRITE  (nout,1105) ufm,iwrd
 1105 FORMAT (a23,' 700, SET',i9,' REQUESTED ON PLOT CARD HAS NOT BEEN',  &
     ' DEFINED.')
 nogo = 1
 1110 iwrd = setd
 GO TO 1130
 1120 iwrd = i
 1130 IF (pset /= 0) GO TO 1180
 pset = iwrd
 GO TO 780
 
!     OUTLINE
 
 1140 IF (pshape /= 1) GO TO 1170
 IF (pcon   == 0) GO TO 780
 pedge = 1
 GO TO 1149
 
!     SHRINK
 
 1142 IF (pedge /= 2) pedge = 75
 IF (pedge == 2) pedge = 75 + 200
!                           SHRINK + HIDDEN
 
 IF (mode > 0) GO TO 780
 CALL rdmode (*1144,*1143,*1180,mode,word)
 1143 CALL rdword (mode,word)
 GO TO 1149
 1144 IF (mode == -2 .AND. fwrd > 0.0 .AND. fwrd <= 1.0) GO TO 1147
 WRITE  (nout,1145) uwm
 1145 FORMAT (a25,', INPUT VALUE ERROR FOR SHRINK.  0.85 IS SUBSTITUED')
 IF (mode == -1) WRITE (nout,1146) iwrd
 1146 FORMAT (5X,'FOR INTEGER VALUE',i5)
 fwrd = 0.85
 1147 j = fwrd*100
 IF (j <  10) j =  10
 IF (j > 100) j = 100
 IF (pedge /= 2) pedge = j
 IF (pedge == 2) pedge = j + 200
!                          SHRINK + HIDDEN
 
 GO TO 1149
 
!     HIDDEN
 
 1148 IF (pedge < 10) pedge = 2
 IF (pedge >= 10 .AND. pedge <= 100) pedge = 200 + pedge
!                                              HIDDEN + SHRINK
 1149 IF (.NOT.disp) GO TO 780
 DO  i = 1,ndef
   IF (deflst(i) == 0) GO TO 1155
 END DO
 pshape = 2
 GO TO 780
 1155 pshape = 3
 GO TO 780
 
!     OFFSET n
!     TURN OFFSET PLOT ON  IF n IS .GE. 0. +n IS MAGNIFYING FACTOR
!     TURN OFFSET PLOT OFF IF n IS .LT. 0
 
 
 1160 IF (mode /= 0) GO TO 1170
 ASSIGN 1165 TO tra
 GO TO 1440
 1165 offscl = iwrd
 IF (offscl >= 0) pedge = 3
 GO TO 780
 
!     UNRECOGNIZABLE OPTION ON THE -PLOT- CARD.
 
 1170 IF (stereo /= 0) GO TO 780
 ERR(1) = 2
 ERR(2) = awrd(1)
 ERR(3) = awrd(2)
 CALL wrtprt (merr,ERR,msg1,nmsg1)
 GO TO 780
 
!     NOFIND
 
!     COMMENTS FROM G.CHAN/UNISYS  11/1990
!     THE 'NOFIND' FEATURE IN NASTRAN PLOTTING COMMANDS IS REALLY NOT
!     NEEDED. IT ONLY LIMITS TO PREVIOUS PLOT CASE. THE FOLLOWING TWO
!     EXAMPLES GIVE EXACTLY THE SAME RESULT IN $ PLOT 2
 
!     $ PLOT 1                           $ PLOT 1
!     FIND SCALE, ORIGIN 100, SET 2      FIND SCALE, ORIGIN 100, SET 2
!     PLOT ORIGIN 100                    PLOT ORIGIN 100
!     $ PLOT 2                           $ PLOT 2
!     PLOT ORIGIN 100                    PLOT NOFIND
!       :
!     (NOTE - ORIGIN 100 IS STILL AVAILABLE
!      IN ANY FOLLOWING PLOT)
!     $ PLOT N
!     PLOT ORIGIN 100
 
 1175 nofind = +1
 IF (lorig == 0) GO TO 1530
 porig  = lorig
 GO TO 780
 
 
 1180 IF (nofind >= 0) GO TO 1185
 IF (fscale /= 0 .OR. for /= 0) GO TO 1182
 IF (prject == 1 .OR. fvp == 0) GO TO 1185
 1182 forg  = 1
 fscale= 1
 isetd = setd
 setd  = MAX0(setd,pset)
 modex = mode
 mode  = -1
 org   = MAX0(1,org)
 CALL find (mode,buf1,b1,setid,deflst)
 nofind= +1
 setd  = isetd
 mode  = modex
 
!     PLOT THIS SET
 
 1185 IF (.NOT.disp) GO TO 1210
 IF (pvectr /= 0 .OR. pshape /= 1 .OR. pedge /= 0) GO TO 1210
 IF (pcon /= 0 .AND. icntvl > 9) GO TO 1210
 IF (pcon /= 0 .AND. icntvl > 13) GO TO 1210
 
!     CREATE A DEFAULT OF SHAPE OR SHAPE + UNDERLAY
 
 DO  i = 1,ndef
   IF (deflst(i) == 0) GO TO 1200
 END DO
 pshape = 2
 GO TO 1210
 1200 pshape = 3
 1210 pset = MAX0(pset,setd)
 
!     DEFAULT OF FIRST DEFINED SET WILL BE USED
 
 CALL gopen  (gpset,deflst(b1),inprew)
 CALL skprec (gpset,pset)
 CALL fread  (gpset,ngpset,1,0)
 
!     TEST FOR CORE NEEDED FOR BOTH UNDEF, DEFOR PLOTS, GRID INDEX
 
 i1 = n2 + ngp + 1
 
!     UNDEFORMED COORDINATES
 
 i2 = i1 + 3*ngpset
 
!     DEFORMATION VALUES
 
 i3 = i2 + 3*ngpset
 
!     REDUCE CORE FOR UNDEFORMED PLOTS
 
 IF (disp) GO TO 1230
 i3 = i2
 n  = 0
 GO TO 1240
 
!     DEFORMED PLOTS NEED X-Y LOCATIONS OF RESULTANT DEFLECTIONS ON
!     FRAME
 
 1230 n = 2*ngpset
 
 1240 IF (i3+n-1 >= defbuf) GO TO 1400
 iused = MAX0(i3+n-1,iused+ngp)
 
 CALL fread (gpset,deflst(n2+1),ngp,0)
 CALL CLOSE (gpset,rew)
 CALL fndset (deflst(n2+1),deflst(i1),buf1-n2,0)
 
 CALL gopen (elset,deflst(b1),inprew)
 IF (pset == 1) GO TO 1280
 CALL skprec (elset,pset-1)
 
 1280 IF (.NOT.stress) GO TO 1290
 IF (icntvl < 4 .OR. DIRECT /= 2) GO TO 1290
 i = b1 + bufsiz
 CALL CLOSE (parm,norew)
 CALL gopen (oes1,deflst(i),norew)
 
 CALL rotat (elset,buf1-n2,deflst(n2+1),deflst(i1))
 
 CALL CLOSE (oes1,norew)
 CALL gopen (parm,deflst(i),norew)
 
 1290 IF (.NOT.disp) GO TO 1320
 
!     CONVERSION FOR ACCEL OR VELOCITY
 
 conv = 1.0
 IF (plttyp == 1) GO TO 1310
 IF (plttyp == 3 .OR. plttyp == 4) GO TO 1300
 
!     VELOCITY
 
 conv = value*twopi
 GO TO 1310
 
!     ACCEL
 
 1300 conv = (value*twopi)**2
 1310 i = 3*bufsiz + b1
 ph1 = ph * twopi/360.0
 CALL getdef (mdef,ph1,mag,conv,plttyp,deflst(i),deflst(n2+1), deflst(i2))
!                  FILE PH  MAG   W   RESP   BUF(1)     GPLST
!                  DEFLECTION
 
!     PRINT THE MAXIMUM FOUND ON THE PLOT FILE
 
 IF (mode >= eor .AND. icolor == 0) CALL head (0,0,-1,defmax)
 ASSIGN 1320 TO incom
 IF (maxdef /= 0.0) defmax = maxdef
 IF (defmax == 0.0 .OR. scale(4) == 0.0) GO TO 1420
 
!                GPLST       ,X         ,U         ,S         ,
 1320 CALL draw (deflst(n2+1),deflst(i1),deflst(i2),deflst(i3),  &
     disp,stereo,defbuf-(i3+n),buf1-n2)
 
!     NOTE - THE NEXT TO LAST ARGUMENT, DEFBUF-(I3+N), IS THE SIZE OF
!            AVAILABLE OPEN CORE. IT IS NOT A POINTER, AND IT IS NOT AN
!            OPEN CORE ARRAY
 
!     OPEN CORE /ZZPLOT/
!     SETID NSETS NDOF      NGP 3*NGPSET 3*NGPSET SCRATCH  N
!     -----+-----+----+----+---+--------+--------+-------+--+--+-+-+-+-+
!          !          N1   N2  I1 (X)   I2 (U)   I3 (S)   DEFBUF ..BUF..
!          !(DEFLST)         /
!                       (GPLST)                      N=2*NGPSET
 
 CALL CLOSE (elset,rew)
 IF (mode >= eor) GO TO 1360
 IF (.NOT.disp) GO TO 1350
 CALL bckrec (mdef)
 1350 pset = iwrd
 IF (.NOT.stress) GO TO 760
 
!     POSITION OES1
 
 i = 1
 ASSIGN 1360 TO incom
 CALL fndset (deflst(n2+1),deflst(i1),buf1-n2,i)
 IF (i == 1) GO TO 760
 GO TO 1420
 
!     END OF A DEFORMATION
 
 1360 CALL stplot (-1)
 IF (prject /= 3 .OR. stereo /= 0) GO TO 1380
 stereo = 1
 CALL sopen (*1430,plttap,deflst(pltbuf),pbufsz)
 j = bfrms
 bfrms = 2
 CALL stplot (pltnum)
 bfrms  = j
 pltnum = pltnum + 1
 IF (.NOT.disp) GO TO 1370
 CALL bckrec (mdef)
 1370 IF (.NOT.stress) GO TO 750
 
!     POSITION OES1
 
 i = 1
 ASSIGN 1360 TO incom
 CALL fndset (deflst(n2+1),deflst(i1),buf1-n2,i)
 IF (i /= 1) GO TO 1420
 GO TO 750
 1380 IF (disp .OR. stress) GO TO 500
 
!     END OF THIS PLOT CARD.
 
 1385 IF (stress) CALL CLOSE (oes1,rew)
 1390 IF (disp  ) CALL CLOSE (mdef,rew)
 GO TO 1430
 
!     INSUFFICIENT CORE TO START PROCESSING
 
 1400 CALL mesage (-8,defbuf,NAME)
 
 1410 CONTINUE
 GO TO 1385
 
!     INCOMPLETE PLOT RESULTED
 
 1420 ERR(1) = 0
 CALL wrtprt (merr,ERR,msg7,nmsg7)
 GO TO incom, (1360,1320)
 
!     FINISHING ONE PLOT
!     ECHO OUT WHICH ORIGIN WAS USED
 
 1430 IF (nogo  /= 0) CALL mesage (-61,0,0)
 IF (porig == 0) GO TO 1550
 ERR(1) = 1
 ERR(2) = origin(porig)
 CALL wrtprt (merr,ERR,used,10)
 CALL WRITE (merr,0,0,1)
 lorig = porig
 porig = 0
 GO TO 1550
 
!     READ AN INTEGER VALUE FROM THE -PLOT- CARD
 
 1440 CALL rdmode (*1450,*790,*1180,mode,word)
 1450 IF (mode == -1) GO TO 1470
 IF (mode == -4) GO TO 1460
 iwrd = fwrd
 GO TO 1470
 1460 iwrd = dwrd
 1470 GO TO tra, (130,150,820,920,1010,1090,1165)
 
!     READ A REAL VALUE FROM THE -PLOT- CARD
 
 1480 CALL rdmode (*1490,*790,*1180,mode,word)
 1490 IF (mode == -4) GO TO 1500
 IF (mode == -1) fwrd = iwrd
 GO TO 1510
 1500 fwrd = dwrd
 1510 GO TO tra, (250,260,900,300)
 
 1520 WRITE  (nout,1525) ufm,plttap
 1525 FORMAT (a23,' 702, PLOT FILE ',a4,' DOES NOT EXIST.')
 nogo = 1
 GO TO 1390
 1530 WRITE  (nout,1535) uwm,lorig
 1535 FORMAT (a25,' 704, NO PREVIOUS PLOT TO INITIATE NOFIND OPERATION')
 
 1550 RETURN
END SUBROUTINE plot
