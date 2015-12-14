SUBROUTINE ifp1
     
!     READS AND INTERPRETS CASE CONTROL DECK FOR NASTRAN
 
 IMPLICIT INTEGER (a-z)
 EXTERNAL        rshift,andf,orf,complf
 LOGICAL :: tapbit,setcd,plotcd,bit64
 REAL :: symseq(360),xcore(1),xintcd
 DIMENSION       case(200,2),xcase(200,2),nifp(2),casen(11),  &
     NAME(2),ttlcd(9),cccd(9),cccds(54),xyprm(5),  &
     outop(15),isubc(5),outpch(13),core(7),corey(401)
 CHARACTER (LEN=23) :: ufm
 COMMON /xmssg / ufm
 COMMON /output/ title(32),subtit(32),label(32),head1(32),  &
     head2(32),head3(32),pltid(32)
 COMMON /system/ sysbuf,otpe,nogo,intp,mpcn,spcn,logfl,loadnn,nlpp,  &
     stftem,ipage,line,tline,maxlin,date(3),tim,iecho,  &
     splots,app,idum,lsystm,dumms(16),nbpw,dummy(28),  &
     isubs,dumz(16),intra,dmz(4),lpch
 COMMON /zzzzzz/ corex(1)
 COMMON /xifp1 / BLANK,bit64
 COMMON /ifp1a / scr1,casecc,is,nwpc,ncpw4,nmodes,icc,nset,  &
     nsym,zzzzbb,istr,isub,lencc,iben,equal,IEOR
 COMMON /ifp1hx/ msst,misset(20)
 COMMON /xsortx/ ibuf41
 EQUIVALENCE     (corex(1) ,xcase(1,1) , case(1,1),corey(1)),  &
     (core(1)  ,xcore( 1)  , corey(401)        ),  &
     (nono     ,outop(15)) , (iaxic  ,dumms(4) ),  &
     (iaxif    ,dumms(15)) , (set    ,cccd( 7) ),  &
     (plot     ,ttlcd(4) ) , (xypl   ,xyprm(1) ),  &
     (outp     ,cccd( 1) ) , (begi   ,cccd( 2) ),  &
     (both     ,outop(1) ) , (NONE   ,outop(2) )
 DATA    nifp  / 4H  IF, 4HP1  /
 DATA    casen / 4HC a , 4HS e , 4H   c, 4H o n, 4H t r, 4H o l,  &
     4H   d, 4H e c, 4H k  , 4H e c, 4H h o  /
 DATA    ibob,   isymcm, loadn , iout2 , inomor          /  &
     0,      0,      1,      0,      0               /
 DATA    outpch/ 11,18 , 21,24 , 27,30 , 33,36 , 152,155,158,168, 171   /
 DATA    blank4, card  , coun  , t     , equal1, nptp   ,dol1   /  &
     4H    , 4HCARD, 4HCOUN, 4HT   , 4H=   , 4HNPTP ,4H$    /
 DATA    NAME  / 4HCASE, 4HCC  /
 DATA    ttlcd / 4HTITL, 4HSUBT, 4HLABE, 4HPLOT, 4HXTIT, 4HYTIT,  &
     4HTCUR, 4HYTTI, 4HYBTI /
 DATA    cccd  / 4HOUTP, 4HBEGI, 4HSYM , 4HSUBC, 4HSYMC, 4HREPC,  &
     4HSET , 4HNCHE, 4HSCAN /
 DATA    cccds / 4HMPC , 4HSPC , 4HLOAD, 4HNLLO, 4HDEFO, 4HTEMP,  &
     4HDLOA, 4HMETH, 4HFREQ, 4HIC  , 4HDISP, 4HVECT,  &
     4HPRES, 4HTHER, 4HSTRE, 4HELST, 4HELFO, 4HFORC,  &
     4HACCE, 4HVELO, 4HSPCF, 4HMAXL, 4HTSTE, 4HSYMS,  &
     4HSUBS, 4HECHO, 4HMODE, 4HLINE, 4HDSCO, 4HK2PP,  &
     4HM2PP, 4HB2PP, 4HTFL , 4HFMET, 4HOFRE, 4HOTIM,  &
     4HCMET, 4HSDAM, 4HSDIS, 4HSVEC, 4HSVEL, 4HSACC,  &
     4HNONL, 4HPLCO, 4HAXIS, 4HHARM, 4HRAND, 4HOLOA,  &
     4HGPFO, 4HESE , 4HMPCF, 4HAERO, 4HGUST, 4HSTRA/
 DATA    all   / 4HALL /, cosi / 4HCOSI/
 DATA    defa  / 4HDEFA/, mat  / 4HMATE/
 DATA    om    / 4HOM  /, oneb / 4H1   /
 DATA    pcdb  / 4HPCDB/, plt1 / 4HPLT1/
 DATA    plt2  / 4HPLT2/, sine / 4HSINE/
 DATA    xycb  / 4HXYCD/, xyou / 4HXYOU/
 DATA    ptit  / 4HPTIT/, flui / 4HFLUI/
 DATA    symm  / 4HSYMM/, anti / 4HANTI/
 DATA    anom  / 4HANOM/
 DATA    xyprm / 4HXYPL, 4HXYPR, 4HXYPU, 4HXYPE, 4HXYPA /
 DATA    outop / 4HBOTH, 4HNONE, 4HUNSO, 4HSORT, 4HPUNC, 4HPRIN,  &
     4HREAL, 4HIMAG, 4HPHAS, 4HNOPR, 4HMAXS, 4HVONM, 4HEXTR, 4HLAYE, 4HNONO/
 
!     INITIALIZE
 
 icc    = 1
 icnt   = 0
 nset   = 0
 nsym   = 0
 isub   = 1
 msst   = 0
 org    = 0
 porg   =-1
 istr   = 1
 ncpw4  = 4
 nwpc   = 20
 jumph  = 0
 npch   = 0
 nogopc = 0
 scr1   = 301
 setcd  =.false.
 plotcd =.false.
 BLANK  = blank4
 bit64  = nbpw == 64
 casecc = NAME(1)
 zzzzbb = 0
 zzzzbb = khrfn1(zzzzbb,1,zzzzbb,4)
 equal  = khrfn1(zzzzbb,1,equal1,1)
 dol    = khrfn1(zzzzbb,1,dol1  ,1)
 iben   = khrfn1(zzzzbb,1,BLANK ,1)
 is     = 9999999
 IEOR   = rshift(complf(0),1)
 nmodes = 1
 lencc  = 200
 DO  j = 1,2
   DO  i = 1,lencc
     case(i,j) = 0
   END DO
 END DO
 case(166,1) = lencc
 DO  j = 1,2
   DO  i = 1,96
     case(i+38,j) = BLANK
   END DO
 END DO
 DO  i = 1,5
   isubc(i) = BLANK
 END DO
 nz = korsz(core) - nwpc - 1
 
!     BLANK TITLE
 
 DO  i = 1,96
   title(i) = BLANK
 END DO
 DO  i = 1,11
   head1(i+9) = casen(i)
 END DO
 head2(  4) = card
 head3(  4) = coun
 head3(  5) = t
 
 i81 = nwpc + 1
 
!     READ IN DATA-- STORE TITLE CARDS
 
 nz   = nz  - sysbuf
 icrq = i81 - nz
 IF (i81 > nz) GO TO 330
 CALL OPEN (*300,scr1,corex(nz+1),1)
 80 CALL xread (*2000,core(1))
 CALL WRITE (scr1,core(1),nwpc,0)
 IF (ibuf41 == -1) GO TO 80
 
!     IS THIS A TITLE SUBTITLE,LABEL,ETC CARD
 
 CALL ifp1f (*80,iword,i2)
 IF (bit64) CALL mvbits (BLANK,0,32,iword,0)
 ASSIGN 80 TO iret1
 istr = 0
 isub = 1
 DO  i = 1,6
   IF (iword == cccd(i)) THEN
      SELECT CASE ( i )
       CASE (    1)
         GO TO 145
       CASE (    2)
         GO TO  340
       CASE (    3)
         GO TO  140
       CASE (    4)
         GO TO  140
       CASE (    5)
         GO TO  140
       CASE (    6)
         GO TO  140
     END SELECT
   END IF
!                                    OUTP BEGI SYM  SUBC SYMC REPC
   
 END DO
 IF (inomor == 1) GO TO 80
 DO  i = 1,3
   IF (iword == ttlcd(i)) THEN
      SELECT CASE ( i )
       CASE (    1)
         GO TO 110
       CASE (    2)
         GO TO  120
       CASE (    3)
         GO TO  130
     END SELECT
   END IF
!                                     TITL SUBT LABE
   
 END DO
 GO TO 80
 110 IF (logfl <= 0) CALL logfil (core(1))
 115 itype = 1
 GO TO 150
 120 itype = 2
 GO TO 150
 130 itype = 3
 GO TO 150
 131 itype = 7
 GO TO 150
 
!     STOP TITLE SEARCH
 
 140 inomor = 1
 GO TO 80
 
!     IDENTIFY PLOT PACKETS
 
 145 CALL xrcard (core(i81),nz,core(1))
 temp = core(i81+5)
 IF (bit64) CALL mvbits (BLANK,0,32,temp,0)
 IF (temp == plot) GO TO 146
 IF (temp == xypl .OR. temp == xyou) GO TO 140
 GO TO 80
 
!     SET PLOT FLAG
 
 146 case(135,1) = 1
 GO TO 140
 
!     FIND EQUAL SIGN COPY REMAINING DATA ON CARD
 
 150 CALL ifp1g (itype,case(1,1),isub)
 GO TO iret1, (80,350)
 
!     FILE ERRORS
 
 300 ip1 = -1
 301 CALL mesage (ip1,FILE,nifp)
 RETURN
 
 310 ip1 = -2
 GO TO 301
 320 ip1 = -3
 GO TO 301
 330 ip1 = -8
 FILE = icrq
 GO TO 301
 340 CALL CLOSE (scr1,1)
 
!     START BUILDING RECORDS
 
 CALL page
 nwdsc  = nwpc + 1
 ASSIGN 350 TO iret1
 ihowdy = 1
 nsym   = 0
 nsyms  = 0
 iun    = 0
 ixypl  = 0
 icasec = 0
 istr   = 1
 nsub   = 0
 msst   = 0
 ibuf1  = nz + 1
 FILE   = scr1
 CALL OPEN (*300,scr1,corex(ibuf1),0)
 nz     = nz - sysbuf
 ibuf2  = nz + 1
 FILE   = casecc
 IF (isubs == 0) GO TO 603
 
!     IN SUBSTRUCTURES, THE CASECC FILE CONTAINS DATA ON THE FRONT.
!     SKIP FILE BEFORE WRITING.
 
 CALL OPEN (*603,FILE,corex(ibuf2),3)
 CALL WRITE (FILE,NAME,2,1)
 350 FILE  = scr1
 icont = 0
 icrq  = i81 - nz
 IF (i81 > nz) GO TO 330
 351 CONTINUE
 CALL READ (*310,*320,scr1,core(1),nwpc,0,flag)
 WRITE  (otpe,360) icc,(core(i),i=1,nwpc)
 360 FORMAT (11X,i8,6X,20A4)
 icc  = icc  + 1
 line = line + 1
 IF (line >= nlpp) CALL page
 IF (dol == khrfn1(zzzzbb, 1,core(1),1))  GO TO 350
 
!     IS THIS TITLE SUBTITLE OR LABEL CARD
 
 CALL ifp1f (*350,iword,i2)
 IF (bit64) CALL mvbits (BLANK,0,32,iword,0)
 DO  i = 1,4
   IF (iword == ttlcd(i) .AND. ibob+ixypl == 0) SELECT CASE ( i )
     CASE (    1)
       GO TO 115
     CASE (    2)
       GO TO  120
     CASE (    3)
       GO TO  130
     CASE (    4)
       GO TO  131
   END SELECT
!                TITL SUBT LABE PLOT
   
 END DO
 IF (iword == ptit .AND. ibob == 1) GO TO 1838
 IF (ixypl /= 1) GO TO 374
 DO  i = 5,9
   IF (iword == ttlcd(i)) GO TO 1838
 END DO
 374 CALL xrcard (core(i81),nz,core(1))
 IF (icont == 1) GO TO 650
 
 IF (bit64) CALL mvbits (BLANK,0,32,core(i81+1),0)
 IF (core(i81+1) == outp) GO TO 590
 IF (core(i81+1) == begi) GO TO 1320
 IF (ibob  == 1) GO TO 1500
 IF (ixypl == 1) GO TO 1836
 IF (core(i81) < 0) GO TO 380
 iword = core(i81+1)
 DO  i = 3,9
   io = i - 2
   IF (iword == cccd(i)) SELECT CASE ( io )
     CASE (    1)
       GO TO 580
     CASE (    2)
       GO TO  1060
     CASE (    3)
       GO TO  1560
     CASE (    4)
       GO TO  1720
     CASE (    5)
       GO TO  1050
     CASE (    6)
       GO TO  791
     CASE (    7)
       GO TO  1055
   END SELECT
!                SYM  SUBC  SYMC  REPC  SET  NCHE  SCAN
   
 END DO
 
 
!     FIND VALUE AFTER EQUAL SIGN
 
 l = 2*IABS(core(i81)) + i81
 DO  i = i81,l
   temp = core(i)
   IF (bit64) CALL mvbits (BLANK,0,32,temp,0)
   IF (temp == equal1) GO TO 377
 END DO
 il = -617
 GO TO 1291
 377 i1 = i + 1
 IF (i == l) i1 = i1 + 1
 
 iword = core(i81+1)
 IF (bit64) CALL mvbits (BLANK,0,32,iword,0)
 DO  i = 1,54
   IF (iword == cccds(i))
   
!               MPC   SPC   LOAD  NLLO  DEFO  TEMP  DLOA  METH  FREQ  &
   GO TO ( 400,  430,  440,  460,  540,  690,  550,  760,  560,
   
!                IC   DISP  VECT  PRES  THER  STRE  ELST  ELFO  FORC  &
   570,  770,  770,  770,  770,  780,  780,  790,  790,
   
!               ACCE  VELO  SPCF  MAXL  TSTE  SYMS  SUBS  ECHO  MODE  &
   800,  810,  820,  610,  620,  630,  630, 1420, 1490,
   
!               LINE  DSCO  K2PP  M2PP  B2PP  TFL   FMET  OFRE  OTIM  &
   1630, 1660, 1680, 1700, 1710, 1730, 1880, 1740, 1740,
   
!               CMET  SDAM  SDIS  SVEC  SVEL  SACC  NONL  PLCO  AXIS  &
   1750, 1760, 1780, 1780, 1790, 1800, 1810, 1665, 1850,
   
!               HARM  RAND  OLOA  GPFO  ESE   MPCF  AERO  GUST  STRA  &
   1860, 1870,  480, 1890, 1900,  405, 1910, 1930, 1950), i
   
 END DO
 
!     UNABLE TO FIND CARD TYPE
 
 380 CALL ifp1d (-601)
 iun  = iun + 1
 IF (iun < 10) GO TO 350
 
!     ASSUME BEGIN BULK MISSING
 
 CALL ifp1d (-611)
 GO TO 1320
 
!     MPC CARD FOUND
 
 400 ik = 2
 GO TO 490
 
!     MPCFORCE CARD
 
 405 ik = 173
 GO TO 830
 
!     TOO MANY SPECIFICATIONS
 
 410 CALL ifp1d  (602)
 GO TO iret, (720,500,860)
 
!     SPC CARD DETECTED
 
 430 ik = 3
 GO TO 490
 
!     LOAD SET SELECTION
 
 440 ik = 4
 GO TO 490
 
!     PNL FOR VDR
 
 460 ik = 10
 GO TO 830
 
!     OUTPUT LOAD SET
 
 480 ik = 17
 GO TO 830
 490 IF (core(i1) <= 0) CALL ifp1d (-617)
 491 ASSIGN 500 TO iret
 
!     SKIP CHECK FOR HARMONIC AS DEFAULT IS NON-ZERO
 
 IF (case(ik,isub) /= 0) GO TO 410
 500 case(ik,isub) = core(i1)
 501 IF (core(i1-1) /= -1) GO TO 520
 
!     CHECK FOR END OF DATA
 
 IF (core(i1+1) == IEOR) GO TO 350
 
!     DATA CARD DID NOT END PROPERLY
 
 503 CONTINUE
 il = -603
 GO TO 1291
 
!     NO INTEGER IN INTEGER FIELD
 
 520 il = -604
 GO TO 1291
 
!     DEFORMATION SET
 
 540 ik = 6
 GO TO 490
 
!     DLOAD CARD
 
 550 ik = 13
 GO TO 490
 
!     FREQUENCY CARD
 
 560 ik = 14
 GO TO 490
 
!     IC CARD
 
 570 ik = 9
 GO TO 490
 
!     SYM CARD
 
 580 nsym = nsym + 1
 IF (nsym-361 < 0) THEN
   GO TO   585
 ELSE IF (nsym-361 == 0) THEN
   GO TO   586
 ELSE
   GO TO  1070
 END IF
 585 symseq(nsym) = 1.0
 GO TO 1070
 586 CALL ifp1d (-633)
 GO TO 1070
 
!     OUTPUT
 
 590 iout2 = 1
 
!     BLANK CHECK
 
 temp = core(i81+5)
 IF (bit64) CALL mvbits (BLANK,0,32,temp,0)
 IF (core(i81+3) == IEOR .AND. core(i81) == 1) GO TO 350
 IF (temp == plot) GO TO 600
 IF (ibob == 1 .AND. .NOT.(setcd.AND.plotcd)) CALL ifp1d (-631)
 IF (temp == xypl .OR. temp == xyou) GO TO 1830
 il = -617
 GO TO 1291
 600 ibob = 1
 
!     TURN ON TRAIL BITS FOR PLOT
 
 core(1) = pcdb
 core(2) = 0
 core(3) = 0
 core(4) = 0
 core(5) = 7777
 core(6) = 0
 core(7) = 0
 CALL wrttrl (core(1))
 
!     CHECK FOR PRESENCE OF PLOT TAPE
!     (SPLOTS COULD BE SET ALREADY BY NASTRAN PLTFLG CARD)
 
 IF (isubs == 0 .AND. .NOT.tapbit(plt1) .AND. .NOT.tapbit(plt2))  &
     CALL ifp1d (-618)
 IF (splots == 0) splots = 1
 IF (splots < 0) splots =-splots
 ASSIGN 605 TO iret3
 GO TO 1321
 
!     CLOSE OPEN STUFF
 
 605 IF (ixypl /= 1) GO TO 602
 
!     TERMINANT XY PACKAGE
 
 ihowdy = -1
 CALL ifp1xy (ihowdy,xintcd)
 CALL CLOSE (xycb,1)
 ixypl = 0
 602 CALL CLOSE (casecc,1)
 
!     OPEN  PCDB
 
 FILE = pcdb
 
!     OPEN WRITE FILE
 
 603 CALL gopen (FILE,corex(ibuf2),1)
 GO TO 350
 
!     MAXLINES CARD
 
 610 maxlin = core(i1)
 GO TO 501
 
!     TIME STEP CARD
 
 620 ik = 38
 GO TO 490
 
!     SYMSEQ AND SUBSEQ
 
 630 IF (isymcm /= 0) GO TO 631
 
!     SYMSEQ  CARD WITHOUT SYMCOM
 
 il = -605
 GO TO 1291
 631 nsymsq = 1
 nsym = 1
 650 IF (nsym-361 < 0) THEN
   GO TO   655
 ELSE IF (nsym-361 == 0) THEN
   GO TO   665
 ELSE
   GO TO   660
 END IF
 655 symseq(nsym) = xcore(i1)
 660 IF (core(i1+1) < 0.0) THEN
   GO TO   670
 ELSE IF (core(i1+1) == 0.0) THEN
   GO TO   680
 ELSE
   GO TO   350
 END IF
 665 CALL ifp1d (-633)
 GO TO 660
 
!     CHECK FOR END OF DATA
 
 670 IF (core(i1+1) == IEOR) GO TO 350
 nsym = nsym + 1
 i1   = i1 + 2
 GO TO 650
 
!     CONTINUATION CARD
 
 680 icont = 1
 nsym  = nsym + 1
 i1    = i81  + 1
 GO TO 351
 
!     TEMPERATURE CARD
 
 690 IF (core(i81) == 2) GO TO 710
 temp = core(i81+5)
 IF (bit64) CALL mvbits (BLANK,0,32,temp,0)
 IF (temp == both) GO TO 710
 IF (temp ==  mat) GO TO 730
 
!     THERMAL LOAD
 
 700 ik = 7
 GO TO 490
 
!     THERMAL + STIFFNESS
 
 710 ASSIGN 720 TO iret
 720 case(8,isub) = core(i1)
 stftem = core(i1)
 IF (isub /= 1) GO TO 740
 GO TO 700
 
!     STIFNESS LOAD
 
 730 ik = 8
 stftem = core(i1)
 IF (isub /= 1) GO TO 740
 GO TO 490
 
!     THERMAL REQUEST AT SUBCASE LEVEL
 
 740 il = 606
 GO TO 1291
 
!     METHOD
 
 760 ik = 5
 GO TO 490
 
!     DISP(PLOT,1) CARD
 
 770 ik = 20
 GO TO 830
 
!     STRESS CARD
 
 780 ik = 23
 GO TO 830
 
!     ELFORCE CARD
 
 790 ik = 26
 GO TO 830
 
!     NCHECK CARD
 
 791 ik = 146
 IF (core(i81  ) ==  1) GO TO 793
 IF (core(i81+5) == -1) GO TO 792
 il = -617
 GO TO 1291
 792 case(ik,isub) = core(i81+6)
 IF (core(i81+7) /= IEOR) GO TO 503
 GO TO 350
 793 case(ik,isub) = 5
 GO TO 350
 
!     ACC
 
 800 ik = 29
 GO TO 830
 
!     VEL CARD
 
 810 ik = 32
 GO TO 830
 
!     SPC FORC
 
 820 ik = 35
 GO TO 830
 
!     OUTPUT SPECIFICATION
!     STRESS AND FORCE FLAGS MAY BE PRE-SET TO 2 (NOPRINT) BY IFP1H
 
 830 ASSIGN 860 TO iret
 IF ((ik == 23 .OR. ik == 26) .AND. case(ik+1,isub) == 2) GO TO 860
 IF (case(ik,isub) /= 0) GO TO 410
 
!     FIND EQUAL SIGN
 
 860 ido = core(i81)
 case(ik+1,isub) = 0
 case(ik+2,isub) = 1
 DO  i = 1,ido
   ii   = i81 + 2*i
   temp = core(ii)
   IF (bit64) CALL mvbits (BLANK,0,32,temp,0)
   IF (temp == equal1) EXIT
   iwrd = core(ii-1)
   DO  io = 4,14
     IF (bit64) CALL mvbits (BLANK,0,32,iwrd,0)
     iop = io - 3
     IF (iwrd == outop(io)) GO TO  &
         (940, 890, 900, 910, 920, 930, 905, 950, 943, 950, 946), iop
!         SORT PUNC PRIN REAL IMAG PHAS NOPR MAXS VONM EXTR LAYE
     
   END DO
   CYCLE
   
!     PUNCH
   
   890 case(ik+1,isub) = case(ik+1,isub) + 4
   CYCLE
   
!     PRINT
   
   900 case(ik+1,isub) = case(ik+1,isub) + 1
   CYCLE
   
!     COMPUTE BUT NO PRINT
!     DEVICE CODE IS 2 (AND SUBPRESS PRINT CODE 1)
   
   905 case(ik+1,isub) = case(ik+1,isub) - MOD(case(ik+1,isub),2) + 2
   CYCLE
   
!     REAL PRINT OUT FORMAT
   
   910 ii = 1
   GO TO 931
   
!     REAL AND IMAGINARY
   
   920 ii = 2
   GO TO 931
   
!     MAGNITUE AND PHASE ANGLE
   
   930 ii = 3
   931 case(ik+2,isub) = ISIGN(ii,case(ik+2,isub))
   CYCLE
   
!     SORT TWO REQUEST
!     (COMMENTS FORM G.C.  7/1989
!     SINCE OES2L FILE HAS NOT BEEN IMPLEMENTED IN ALL DMAPS, SORT2
!     STRESS REQUEST ON LAYERED ELEMENTS IS NOT AVAILABLE)
   
   940 temp = core(ii)
   IF (bit64) CALL mvbits (BLANK,0,32,temp,0)
   IF (temp == oneb) CYCLE
   IF (ik == 23 .AND. case(183,isub) >= 2) CALL ifp1d (-645)
   case(ik+2,isub) = -IABS(case(ik+2,isub))
   CYCLE
   
!     VON MISES STRESS
!     (183 WORD ON CASECC, FIRST RIGHT-MOST BIT)
   
   943 case(183,isub) = orf(case(183,isub),1)
   CYCLE
   
!     LAYER STRESSES FOR COMPOSITE ELEMENTS
!     (183 WORD ON CASECC, SECOND RIGHT-MOST BIT)
!     (SORT2 STRESS REQUEST ON LAYERED ELEMENTS NOT AVAILABLE)
   
   946 IF (ik /= 23) CALL ifp1d (-646)
   IF (ik == 23 .AND. case(25,isub) < 0) CALL ifp1d (-645)
   case(183,isub) = orf(case(183,isub),2)
   
 END DO
 960 IF (case(ik+1,isub) == 0) case(ik+1,isub) = 1
 IF (core(ii+1) /=  0) GO TO 962
 CALL ifp1d (610)
 GO TO 970
 962 temp = core(ii+1)
 IF (bit64) CALL mvbits (BLANK,0,32,temp,0)
 IF (temp == all) GO TO 970
 IF (temp == NONE .OR. temp == nono) GO TO 980
 IF (core(ii+1) == -1) GO TO 964
 il = -617
 GO TO 1291
 964 i1 = ii + 2
 GO TO 990
 
!     ALL SPECIFIED -- SET SET NO. MINUS
 
 970 case(ik,isub) = -1
 GO TO 1042
 
!     NONE SPECIFIED
 
 980 case(ik,isub) = NONE
 GO TO 1042
 
!     FIND SET NUMBER
 
 990 IF (nset /= 0) GO TO 1020
 
!     UNDEFINED SET ID ON CARD
 
 1000 CALL ifp1d (-608)
 GO TO 350
 1020 jj = nwdsc
 DO  il = 1,nset
   IF (core(jj) == core(i1)) GO TO 1040
   jj = jj + core(jj+1) + 3
 END DO
 GO TO 1000
 1040 case(ik,isub) =core(i1)
 1042 IF (core(ii+3) /= IEOR) GO TO 503
 GO TO 350
 
!     SET CARD
 
 1050 nset = nset + 1
 CALL ifp1c (i81,nz)
 GO TO 350
 
!     SCAN CARD
 
 1055 CALL ifp1h (i81,nz,jumph)
 GO TO 350
 
!     SUBCASE
 
 1060 temp = core(i81+2)
 IF (bit64) CALL mvbits (BLANK,0,32,temp,0)
 IF (temp   == om) GO TO 1560
 IF (isymcm ==  1) GO TO 1330
 nsym  = 0
 nsyms = nsyms + 1
 IF (nsyms-361 < 0) THEN
   GO TO  1062
 ELSE IF (nsyms-361 == 0) THEN
   GO TO  1064
 ELSE
   GO TO  1070
 END IF
 1062 symseq(nsyms) = 1.0
 GO TO 1070
 1064 CALL ifp1d (-633)
 1070 ASSIGN 350 TO iret3
 IF (isub == 2) GO TO 1080
 isub  = 2
 loadn = core(i81+4)
 CALL ifp1f (*350,iword,i2)
 IF (bit64) CALL mvbits (BLANK,0,32,iword,0)
 DO  i = 1,5
   isubc(i) = core(i2)
   i2 = i2 + 1
 END DO
 IF (core(i81+3) + 1 == 0.0) THEN
   GO TO   350
 ELSE
   GO TO  1290
 END IF
 
!     TURN STRESS AND FORCE NO-PRINT FLAGS ON IF INTERACTIVE FLAG IS ON
 
 1080 IF (intra < 2) GO TO 1085
 case(24,isub) = orf(case(24,isub),8)
 case(27,isub) = orf(case(27,isub),8)
 
 1085 case(1,isub) = loadn
 IF (core(i81+4) <= loadn+nmodes-1) GO TO 1310
 loadn = core(i81+4)
 IF (core(i81+3) /= -1) GO TO 1310
 1090 IF (case(137,1) ==  0) case(137,1) = 1
 CALL ifp1e (isubc(1),symseq,nwdsc,i81,icaste)
 stftem = icaste
 nsub = nsub + nmodes
 
!     CHECK SET NOS. THAT WERE SPECIFIED AFTER SCAN CARDS
 
!     FORM G.C./UNISYS   4/1990
!     IFP1H IS BY-PASSING THIS NEW CODE HERE (MSST=0) BECAUSE SET DATA
!     IS NOT AVAILABLE HERE. SAVE THIS CODE FOR FURTHER INVESTIGATION.
 
 IF (msst == 0) GO TO 1281
 mm = 0
 ll = lencc + case(lencc,isub) + 1
 DO  m = 1,msst
   i  = ll
   mset = misset(m)
   
!     WRITE (6,2345) MSET,MSST,LL
   
   1091 iset = case(i,isub)
   
!     LX1 = I - 3
!     LX2 = I + 3
!     WRITE (6,6789) ISET,(CASE(LX,ISUB),LX=LX1,LX2)
   
   IF (iset == 0) CYCLE
   IF (mset-iset == 0) THEN
     GO TO  1093
   END IF
   1092 i = i + case(i+1,isub)
   IF (i >= 400) CYCLE
   GO TO 1091
   1093 misset(m) = 0
   mm = mm + 1
 END DO
 IF (mm == msst) GO TO 1281
 DO  m = 1,msst
   IF (misset(m) == 0) CYCLE
   WRITE  (otpe,1095) ufm,misset(m)
   1095 FORMAT (a23,' 608A, UNIDENTIFIED SET',i8,' WAS REQUESTED FOR ',  &
       'SCAN')
   nogo = 1
 END DO
 
 1281 GO TO iret3, (350,1370,605,1835)
 
!     SUBCASE ID MISSING
 
 1290 il = -609
 loadn = case(1,2)
 1291 CALL ifp1d (il)
 GO TO 350
 1310 CALL ifp1d (-609)
 loadn = case(1,2)
 GO TO 1090
 
!     BEGIN BULK
 
 1320 ASSIGN 1370 TO iret3
 1321 core(i81+3) = -1
 core(i81+4) = 9999999
 IF (icasec == 1) GO TO 1281
 icasec = 1
 IF (isymcm == 1) GO TO 1330
 nsym = 0
 GO TO 1080
 
!     PUT OUT SUBCOM OR SYMCOM RECORD
 
 1330 isymcm = 0
 1340 IF (nsymsq /= 0 .OR. nsym /= 0) GO TO 1360
 
!     NO SUBSEQ OR SYMSEQ CARD
 
 nsym = nsyms
 
 1360 nsymsq = 0
 case(lencc,2) = MAX0(nsym,0)
 case(16,2) = nsym
 GO TO 1080
 1370 CALL CLOSE (scr1,1)
 IF (ibob /= 1 .AND. ixypl /= 1) CALL CLOSE (casecc,1)
 IF (ibob == 1) CALL CLOSE (pcdb,1)
 IF (ibob == 1 .AND. .NOT.(setcd.AND.plotcd)) CALL ifp1d (-631)
 IF (ixypl /= 1) GO TO 1371
 
!     TERMINATE XYPLOT PACKAGE
 
 ihowdy = -1
 CALL ifp1xy (ihowdy,xintcd)
 CALL CLOSE  (xycb,1)
 
!     PUT CASECC ON NPTP
 
 1371 CONTINUE
 FILE  = casecc
 CALL OPEN (*300,casecc,corex(ibuf1),0)
 FILE  = nptp
 maxcc = 0
 CALL OPEN (*300,nptp,corex(ibuf2),3)
 1380 CALL READ (*1400,*1390,casecc,core(1),nz,0,flag)
 icrq  = nz
 GO TO 330
 1390 CALL WRITE (nptp,core(1),flag,1)
 maxcc = MAX0(maxcc,flag)
 
!     CHECK ANY PUNCH REQUEST  ON OUTPUT DATA BLOCKS
 
 IF (npch == 1 .OR. flag < 166) GO TO 1380
 DO  i = 1,13
   j = outpch(i)
   IF (andf(core(j),4) /= 0) GO TO 1395
 END DO
 GO TO 1380
 1395 npch = 1
 GO TO 1380
 1400 CALL CLOSE (casecc,1)
 CALL eof (nptp)
 CALL CLOSE (nptp,2)
 IF (splots < 0) splots = 0
 
!     IF THIS IS A RESTART  SET CHANGE FLAGS IN IFP1B
 
 IF (app < 0) CALL ifp1b
 IF (iun /= 0) CALL ifp1d (-612)
 CALL makmcb (core,casecc,nsub,0,0)
 core(2) = nsub
 core(4) = maxcc
 CALL wrttrl (core)
 
!     SET NOGO FLAG TO -9 IF ERROR IN BULKDATA AND PLOT COMMANDS
!     SET NOGO FLAG TO POSITIVE IF ERROR IN BULKDATA, AND NOT IN PLOT
!     SET NOGO FLAG TO NEGATIVE IF NO ERROR IN BULKDATA, BUT IN PLOT
!     PUNCH AN IDENTIFICATION CARD IF PUNCH IS REQUESTED ON OUTPUT DATA,
!     AND PRINT SCAN KEYWORDS IF ERROR FLAG (JUMPH) WAS TURNED ON
 
 IF (nogo /= 0 .AND. nogopc == -1) nogo = -9
 IF (nogo == 0) nogo = nogopc
 IF (npch == 1) WRITE (lpch,1415) (title(j),j=1,17)
 1415 FORMAT (2H$ ,17A4)
 IF (jumph == 1) CALL ifp1h (0,0,2)
 RETURN
 
!     ECHO REQUEST
 
 1420 iecho = 0
 ido = core(i81) - 2
 DO  i = 1,ido
   iwrd = core(i1)
   IF (bit64) CALL mvbits (BLANK,0,32,iwrd,0)
   DO  io = 1,5
     IF (iwrd == outop(io)) THEN
        SELECT CASE ( io )
         CASE (    1)
           GO TO 1435
         CASE (    2)
           GO TO  1480
         CASE (    3)
           GO TO  1440
         CASE (    4)
           GO TO  1430
         CASE (    5)
           GO TO  1431
       END SELECT
     END IF
!                                     BOTH  NONE  UNSO  SORT  PUNC
     
   END DO
   IF (iwrd == outop(15)) GO TO 1470
   CALL ifp1d (629)
   GO TO 1432
   
!     SORTED ECHO
   
   1430 CONTINUE
   IF (andf(iecho,2) /= 0) CALL ifp1d (629)
   1432 iecho = orf(iecho,2)
   GO TO 1450
   
!     PUNCH ECHO
   
   1431 CONTINUE
   IF (andf(iecho,4) /= 0) CALL ifp1d (629)
   iecho = orf(iecho,4)
   npch  = 1
   GO TO 1450
   
!     BOTH ECHO
   
   1435 CONTINUE
   IF (andf(iecho,3) /= 0) CALL ifp1d (629)
   iecho = orf(iecho,3)
   GO TO 1450
   
!     UNSORTED ECHO
   
   1440 CONTINUE
   IF (andf(iecho,1) /= 0) CALL ifp1d (629)
   iecho = orf(iecho,1)
   1450 i1 = i1 + 2
 END DO
 
 GO TO 350
 
!     NONO ECHO - ABSOLUTELY NO ECHO, NO EVEN IN RESTART
 
 1470 io = 16
 
!     NONE ECHO
 
 1480 CONTINUE
 IF (iecho /= 0 .OR. i < ido) CALL ifp1d (630)
 iecho = -1
 IF (io == 16) iecho = -2
 GO TO 350
 
!     LOOP CONTROL FOR EIGENVALUE
 
 1490 nmodes = core(i1)
 GO TO 350
 
!     PLOT DATA FOR BO BATA
 
 1500 i1 = i81
 
!     TEST FOR REQUIRED PLOT AND SET CARDS IN STRUCTURE PLOT OUTPUT PKG
 
 temp = core(i81+2)
 IF (bit64) CALL mvbits (BLANK,0,32,temp,0)
 IF (core(i81+1) == plot .AND. temp == BLANK) plotcd =.true.
 IF (core(i81+1) == set) setcd  = .true.
 
!     TEST FOR XYPLOT COMMAND CARDS IN STRUCTURE PLOT OUTPUT PACKAGE
 
 iwrd = core(i81+1)
 DO  i = 1,5
   IF (iwrd == xyprm(i)) CALL ifp1d (-632)
 END DO
 
!     TEST FORMAT OF PLOT COMMAND CARDS
 
 i = nogo
 nogo = 0
 CALL ifp1pc (i81,icnt,xintcd,org,porg)
 IF (nogo /= 0) nogopc = -1
 nogo = i
 
!     COMPUTE LENGTH OF RECORD
 
 ik = 0
 1510 IF (core(i1) < 0.0) THEN
   GO TO  1520
 ELSE IF (core(i1) == 0.0) THEN
   GO TO  1550
 ELSE
   GO TO  1530
 END IF
 1520 CONTINUE
 ip = 2
 GO TO 1540
 1530 IF (core( i1) == IEOR) GO TO 1550
 ip = 2*core(i1) + 1
 1540 ik = ik + ip
 i1 = i1 + ip
 GO TO 1510
 1550 CALL WRITE (pcdb,core(i81),ik+1,1)
 GO TO 350
 
!     PLOT TITLE CARD
 
 1555 core(i81  ) = 10
 core(i81+1) = iword
 core(i81+2) = BLANK
 CALL ifp1g (itype,core(i81+3),1)
 core(i81+21) = 9999999
 ik = 21
 GO TO 1550
 
!     SYMCOM OR SUBCOM CARD
 
 1560 IF (isymcm == 0) GO TO 1570
 ASSIGN 350 TO iret3
 GO TO 1340
 1570 isymcm = 1
 nsymsq = 0
 GO TO 1070
 
!     LINE CARD - NLPP BOTTOM-LIMITED TO 10
 
 1630 CONTINUE
 IF (core(i1-1) /= -1) GO TO 520
 IF (IABS(core(i1)) > 0) nlpp = IABS(core(i1))
 IF (nlpp < 10) nlpp = 10
 GO TO 350
 
!     DIFFERENTIAL STIFFNESS OR PIECEWISE LINEAR COEFFICIENT SET
 
 1660 ik = 138
 GO TO 1670
 1665 ik = 164
 1670 temp = core(i1)
 IF (bit64) CALL mvbits (BLANK,0,32,temp,0)
 IF (temp /= defa) GO TO 490
 core(i1  ) = -1
 core(i1+1) = IEOR
 core(i1-1) = -1
 GO TO 491
 
!     K2PP
 
 1680 ik = 139
 1690 case(ik  ,isub) = core(i1  )
 case(ik+1,isub) = core(i1+1)
 GO TO 350
 
!     M2PP
 
 1700 ik = 141
 GO TO 1690
 
!     B2PP
 
 1710 ik = 143
 GO TO 1690
 
!     REPRINT OF ABOVE CASE
 
 1720 nsym = -1
 IF(isub /= 2) CALL ifp1d (-607)
 GO TO 1560
 
!     TRANSFER FUNCTION SELECTION
 
 1730 ik = 15
 GO TO 490
 
!     OUTPUT FREQUENCY LIST SET
 
 1740 ik = 145
 temp = core(i1)
 IF (bit64) CALL mvbits (BLANK,0,32,temp,0)
 IF (temp /= all) GO TO 830
 core(i1  ) = -1
 core(i1-1) = -1
 core(i1+1) = IEOR
 GO TO 491
 
!     COMPLEX EIGENVALUE METHOD
 
 1750 ik = 148
 GO TO 490
 
!     STRUCTURAL DAMPING TABLE
 
 1760 ik = 149
 GO TO 490
 
!     INERTIA RELIEF SET SELECTION
 
!1770 IK = 150
!     GO TO 490
 
!     ANALYSIS SET FOR VDR
 
 1780 ik = 151
 GO TO 830
 
!     ANALYSIS VELOCITY
 
 1790 ik = 154
 GO TO 830
 
!     ANALYSIS ACCELERATION
 
 1800 ik = 157
 GO TO 830
 
!     NON LINEAR FORCE VECTOR FOR TRANSIENT ANALYSIS
 
 1810 ik = 160
 GO TO 490
 
!     X-Y PLOTTER PACKAGE
 
 1830 ASSIGN 1835 TO iret3
 GO TO 1321
 1835 CALL CLOSE (casecc,1)
 DO  i = 2,6
   core(i) = 0
 END DO
 core(1) = xycb
 core(7) = 1
 CALL wrttrl (core(1))
 
!     OPEN XYCB
 
 IF (ibob /= 1) GO TO 1837
 CALL CLOSE (pcdb,1)
 ibob  = 0
 1837 FILE  = xycb
 ixypl = 1
 i81   = nwpc + 1
 GO TO 603
 
!     AXIS TITLE CARDS
 
 1838 itype = 8
 core(1) = iword
 DO  i = 1,32
   k = i81 + i - 1
   core(k) = BLANK
 END DO
 IF (ibob == 1) GO TO 1555
 CALL ifp1g (itype,core(i81),1)
 
!     PROCESS XYPLOTTER CARD
 
 1836 CALL ifp1xy (ihowdy,xintcd)
 GO TO 350
 
!     DELETE SETS FOR FORCE
 
!1840 IK = 161
!     GO TO 490
 
!     AXISYM CARD
 
 1850 temp = core(i1)
 IF (bit64) CALL mvbits (BLANK,0,32,temp,0)
 IF (temp == sine) GO TO 1851
 IF (temp == cosi) GO TO 1852
 IF (temp == flui) GO TO 1852
 IF (temp == symm) GO TO 1853
 IF (temp == anti) GO TO 1854
 IF (temp == anom) GO TO 1855
 
!     ILLEGAL  SPECIFICATION
 
 il = -617
 GO TO 1291
 1851 case(136,isub) = 1
 iaxic = 1
 GO TO 350
 1852 case(136,isub) = 2
 IF (temp == cosi) iaxic = 1
 IF (temp == flui) iaxif = 1
 GO TO 350
 1853 case(136,isub) = -2
 GO TO 1856
 1854 case(136,isub) = -1
 GO TO 1856
 1855 case(136,isub) = -30
 GO TO 350
 1856 temp = core(i1+1)
 IF (bit64) CALL mvbits (BLANK,0,32,temp,0)
 IF (temp == anom) case(136,isub) = case(136,isub)*10
 GO TO 350
 
!     HARMONIC SELECTOR
 
 1860 ik = 137
 temp = core(i1)
 IF (bit64) CALL mvbits (BLANK,0,32,temp,0)
 IF (temp ==  all) GO TO 1861
 IF (temp == NONE) GO TO 1862
 core(i1) = core(i1) + 1
 GO TO 490
 1861 core(i1) = -1
 GO TO 1863
 1862 case(137,1)= 0
 core(i1  ) = 0
 1863 core(i1-1) = -1
 core(i1+1) = IEOR
 GO TO 491
 
!     RANDOM SET SELECTION
 
 1870 ik = 163
 GO TO 490
 
!     FMETHOD
 
 1880 ik = 165
 GO TO 490
 
!     GRID POINT FORCE REQUEST
 
 1890 ik = 167
 GO TO 830
 
!     ELEMENT STRAIN ENERGY
 
 1900 ik = 170
 GO TO 830
 
!     AEROFORCE OUTPUT REQUEST
 
 1910 ik = 176
 GO TO 830
 
!     AEROELASTIC GUST LOAD REQUEST
 
 1930 ik = 179
 GO TO 490
 
!     STRAIN CARD
!     (180 THRU 182 WORDS OF CASECC)
 
 1950 ik = 180
 GO TO 830
 
!     EOF ON INPUT UNIT
 
 2000 CALL ifp1d  (-624)
 CALL mesage (-37,0,nifp)
 RETURN
END SUBROUTINE ifp1
