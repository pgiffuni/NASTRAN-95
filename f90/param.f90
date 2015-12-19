SUBROUTINE param (setid,xx,buf4)
     
!     THIS PARAM ROUTINE IS CALLED ONLY BY DPLOT, WHICH IS THE DRIVER
!     OF THE PLOT MODULE           ==== == =====
 
!     THE DRIVER FOR THE PARAM MODULE IS QPARAM
 
 
 INTEGER, INTENT(IN OUT)                  :: setid(1)
 INTEGER, INTENT(IN OUT)                  :: xx(1)
 INTEGER, INTENT(IN)                      :: buf4
 LOGICAL :: test
 INTEGER :: buf1     ,bufsiz   ,title   ,  &
     prnt    ,parm     ,pltbuf   ,camera   ,bframs   ,pltmod  ,  &
     tapden  ,pensiz   ,penclr   ,paptyp   ,axis     ,daxis   ,  &
     fvp     ,prject   ,for      ,org      ,origin   ,ploter  ,  &
     word    ,awrd(2)  ,ERR(3)   ,BLANK    ,pltnam(2),eor     ,  &
     tra     ,where    ,DIRECT   ,fscale   ,pltype   ,fpltit  ,  &
     pltitl  ,savtit(96),msg1(20),msg2(20) ,msg4(22) ,msg5(16),  &
     anti    ,axisd(7) ,both     ,bpi      ,by       ,colo    ,  &
     comm    ,defo     ,dens     ,disp     ,even     ,film    ,  &
     fram    ,hmode    ,hplot(2) ,hkey(19) ,hx       ,oesx    ,  &
     nkwd(3) ,icnda(20),plan     ,poin     ,pape     ,sepa    ,  &
     size    ,stre     ,symm     ,TYPE     ,z1       ,z2      ,  &
     fill    ,color    ,layer    ,oes1     ,oes1l    ,onrgy1
 REAL :: maxdef
 DOUBLE PRECISION :: dwrd
 COMMON /system/  ksystm(65)
 COMMON /output/  title(96)
 COMMON /BLANK /  skp11(3) ,prnt     ,skp12(6) ,parm     ,skp2(9) ,  &
     merr     ,skpit(2) ,oesx
 COMMON /xxparm/  pltbuf   ,camera   ,bframs   ,pltmod(2),tapden  ,  &
     npens    ,papsiz(2),paptyp(2),pensiz(8)         ,  &
     penclr(8,2),penpap ,scale(2) ,fscale   ,maxdef  ,  &
     defmax   ,axis(3)  ,daxis(3) ,vangle(3)         ,  &
     skpvue(6),fvp      ,vanpnt(5),d02      ,d03     ,  &
     prject   ,s0s      ,for      ,org      ,norg    ,  &
     origin(11),edge(11,4),xy(11,3),ncntr   ,cntr(50),  &
     icntvl   ,where    ,DIRECT   ,subcas   ,flag    ,  &
     DATA     ,lasset   ,fpltit   ,pltitl(17)        , color    ,layer
 COMMON /pltdat/  model    ,ploter   ,skpplt(17),chrscl  ,skpa(2) ,  &
     cntsin   ,skpd1(6) ,pltype    ,skpd2(3),cntin3
 EQUIVALENCE      (ksystm(1),bufsiz) ,(pape,hkey(10))    ,  &
     (word,awrd(1),dwrd ,fwrd,iwrd)
 
!     THE FOLLOWING ARE THE ALLOWABLE FIRST WORDS ON THE LOGICAL CARD.
!     THE PROJECTION DETERMINES HOW MANY WORDS ARE CHECKED.
 
!     OES1   IS THE NORMAL STRESS FILE, 111
!     OES1L  IS THE LAYER COMPOSITE STRESS FILE, 112
!     ONRGY1 IS THE ELEMENT  STRAIN ENERGY FILE, 113
 
 DATA oes1 , oes1l , onrgy1 /111,112,113   /
 DATA nkwd / 17,19,19/     , blank4 /4H    /
 DATA hkey / 4HFIND, 4HVIEW, 4HAXES, 4HMAXI, 4HORTH, 4HPERS,  &
     4HSTER, 4HCONT, 4HCAME, 4HPAPE, 4HPEN , 4HBLAN,  &
     4HORIG, 4HSCAL, 4HCSCA, 4HPROJ, 4HPTIT, 4HOCUL, 4HVANT/
 
!     THE FOLLOWING ARE RECOGNIZABLE PARAMETERS
 
 DATA axisd/ 2HMZ   , 2HMY, 2HMX, 0, 1HX,  1HY,     1HZ        /,  &
     anti / 4HANTI/, both/ 4HBOTH/, bpi / 4HBPI /, by  /4HBY  /,  &
     colo / 4HCOLO/, defo/ 4HDEFO/, dens/ 4HDENS/, film/4HFILM/,  &
     fram / 4HFRAM/,hmode/ 4HMODE/,hplot/ 4HPLOT,  4HTER      /,  &
     hx   / 4HX   /, plan/ 4HPLAN/, poin/ 4HPOIN/, sepa/4HSEPA/,  &
     size / 4HSIZE/, symm/ 4HSYMM/, TYPE/ 4HTYPE/, & 
!     CONTOUR PLOTTING
  disp / 4HDISP/, stre/ 4HSTRE/, even/ 4HEVEN/, laye/ 4HLAYE/,  &
     list / 4HLIST/, z1  / 2HZ1  /, z2  / 2HZ2  /, MAX / 3HMAX /,  &
     mid  / 3HMID /, comm/ 4HCOMM/, loca/ 4HLOCA/, fill/ 4HFILL/
 
 DATA  icnda  /4HMAJP, 4HMINP, 4HMAXS, 4HXNOR, 4HYNOR, 4HZNOR,  &
     4HXYSH, 4HXZSH, 4HYZSH, 4HXDIS, 4HYDIS, 4HZDIS, 4HMAGN,  &
     4HNRM1, 4HNRM2, 4HSH12, 4HSH1Z, 4HSH2Z, 4HBDSH, 4HSTRA/
 
 DATA  eor   , BLANK/ 1000000, 1H  /,  &
     nmsg5 , msg5 / 16,4H(25X, 4H,31H, 4HMORE, 4H tha, 4HN 50,  &
     4H con, 4HTOUR,   4HS sp, 4HECIF, 4HIED,, 4H1P,e, 4H14.6,  &
     4H,9H , 4HREJE,   4HCTED, 4H)   /
 DATA  nmsg1 / 20   /
 DATA  msg1  / 4H(34X  ,4H,45H   ,4HAN a   ,4HTTEM   ,4HPT h   ,  &
     4HAS b  ,4HEEN    ,4HMADE   ,4H TO    ,4HDEFI   ,  &
     4HNE m  ,4HORE    ,4HTHAN   ,4H ,i2   ,4H,17H   ,  &
     4H dis  ,4HTINC   ,4HT OR   ,4HIGIN   ,4HS)     /
 DATA  nmsg2 / 20   /
 DATA  msg2  / 4H(30X  ,4H,34H   ,4HAN u   ,4HNREC   ,4HOGNI   ,  &
     4HZABL  ,4HE pl   ,4HOT p   ,4HARAM   ,4HETER   ,  &
     4H (,2  ,4HA4,2   ,4H9H)    ,4HHAS    ,4HBEEN   ,  &
     4H det  ,4HECTE   ,4HD -    ,4HIGNO   ,4HRED)   /
 DATA  nmsg4 / 22   /
 DATA  msg4  / 4H(25X  ,4H,4HP   ,4HEN ,   ,4HI4,6   ,4H9H i   ,  &
     4HS no  ,4HT a    ,4HLEGA   ,4HL pe   ,4HN nu   ,  &
     4HMBER  ,4H for   ,4H thi   ,4HS pl   ,4HOTTE   ,  &
     4HR. p  ,4HEN 1   ,4H wil   ,4HL be   ,4H red   , 4HEFIN  ,4HED.)    /
 DATA  test  / .false. /
 
 
!     COMMENTS FROM G.CHAN/UNISYS ABOUT THE NOFIND FLAG       11/1990
 
!     THE NOFIND FLAG WAS TOO CONFUSING BEFORE. I'M SETTING THE NEW RULE
!     HERE
 
!     NOFIND FLAG IS USED IN PARAM AND PLOT ROUTINES ONLY. ITS USE IS
!     TO INDICATE WHETHER SUBROUTINE FIND SHOULD BE CALLED.
!     (SUBROUTINE FIND COMPUTES THE NEW ORIGIN, FRAME SIZE, NEW VIEW,
!     VANTAGE POINT ETC. DUE TO CERTAIN PLOT PARAMETERS).
!     NOFIND FLAG CAN BE SET BY USER VIA THE FIND AND NOFIND COMMANDS,
!     OR IT IS SET AUTOMATICALLY BY THIS PARAM SUBROUTINE.
 
!      NOFIND                  ACTION
!     --------    ----------------------------------------------------
!        -1       FIND ROUTINE SHOULD BE CALLED IN NEXT OPPORTUNITY
!                 BEFORE THE ACTUAL PLOTTING
!        +1       (1) A NOFIND CARD WAS ENCOUNTERED. USER WANTS TO KEEP
!                 ALL PARAMETERS AS IN THE PREVIOUS PLOT CASE, OR
!                 (2) FIND ROUTINE WAS JUST CALLED. PROGRAM SHOULD NOT
!                 CALL FIND AGAIN
!         0       THE CURRENT STATUS OF ALL PARAMETERS THAT WERE FOUND
!                 BY PREVIOUS FIND REMAIN UNCHANGED. HOWEVER, ANY
!                 CHANGE IN THE PLOT PARAMETERS BY THE USER (SCALE,
!                 CSCALE, VIEW, VENTAGE POINT, REGION, ORIGIN, PLOTTER,
!                 MAX.DEFORMATION, PROJECTION AND PAPER SIZE) WILL
!                 CHANGE NOFIND FLAG TO -1
 
!     IF A FIND COMMAND IS ENCOUNTERED, SUBROUTINE FIND IS CALLED
!     IMMEDIATELY AND UNCONDISIONALLY, THEN NOFIND FLAG IS SET TO +1
 
!     IF USER HAS ALREADY ONE OR MORE ORIGINS, AND IF HE USES A FIND
!     CARD TO FIND ANOTHER ORIGIN, BUT THE NEXT PLOT CARD DOES NOT USE
!     THIS NEWLY DEFINED ORIGIN, A WARNING MESSAGE SHOULD BE ISSUED TO
!     INFORM THE USER THAT THE DEFAULT ORIGIN, WHICH IS THE FIRST
!     DEFINDED ORIGIN, IS GOING TO BE USED, NOT THE ONE HE JUST DEFINED
 
 nofind = -1
 lasset = 0
 CALL pltset
 buf1 = buf4 + 3*bufsiz
 
!     SAVE THE TITLE, SUBTITLE AND LABEL IF DEFORMED PLOTS ...
 
 IF (prnt >= 0) GO TO 30
 DO  i = 1,96
   savtit(i) = title(i)
 END DO
 20 nofind = 0
 30 CALL rdmodx (parm,mode,word)
 40 CALL READ (*1800,*1800,parm,mode,1,0,i)
 IF (mode < 0) THEN
   GO TO    50
 ELSE IF (mode == 0) THEN
   GO TO    40
 ELSE
   GO TO    60
 END IF
 50 i = 1
 IF (mode == -4) i = 2
 CALL fread (parm,0,-i,0)
 GO TO 40
 60 IF (mode < eor) GO TO 70
 CALL fread (parm,0,0,1)
 GO TO 40
 70 mode = mode + 1
 CALL rdword (mode,word)
 CALL rdword (mode,word)
 IF (awrd(1) /= hplot(1)) GO TO 160
 IF (awrd(2) ==    BLANK) GO TO 110
 IF (awrd(2) == hplot(2)) GO TO 900
 GO TO 1750
 
!     FIND
 
 100 CALL find (mode,buf1,buf4,setid,xx)
 nofind = +1
 IF (mode >= 0) GO TO 30
 mode = modex
 GO TO 130
 
!     PLOT
 
 110 IF (test) GO TO 130
 
!         WHEN PLOTTER OR PROJECTION WERE HIT
!              FSCALE=FOR=FVP=1
!              PROJECTION=KWRD-4, SOME NUMBER
!         WHEN SCALE IS HIT,       FSCALE SET TO 0
!         WHEN VANTAGE POINT IS HEIT, FVP SET TO 0
!         WHEN ORIGIN IS HIT,         ORG SET TO 0
 
 IF (fscale /= 0 .OR. for /= 0) GO TO 120
 IF (prject == 1 .OR. fvp == 0) GO TO 130
 120 modex = mode
 mode  = -1
 org   = MAX0(1,org)
 GO TO 100
 130 CALL plot (mode,buf1,buf4,setid,xx,nofind)
 oesx = oes1
 IF (nofind == -1) org = MAX0(1,org)
 GO TO 20
 
!     PLOT PARAMETER CARD.
 
 140 IF (mode <= 0) CALL rdmode (*140,*150,*40,mode,word)
 150 CALL rdword (mode,word)
 160 i = nkwd(prject)
 DO  kwrd = 1,i
   IF (hkey(kwrd) == word) GO TO 200
 END DO
 GO TO 1750
 
 200 SELECT CASE ( kwrd )
   CASE (    1)
     GO TO 100
   CASE (    2)
     GO TO  1230
   CASE (    3)
     GO TO   250
   CASE (    4)
     GO TO   500
   CASE (    5)
     GO TO   230
   CASE (    6)
     GO TO   230
   CASE (    7)
     GO TO   230
   CASE (    8)
     GO TO  1300
   CASE (    9)
     GO TO   400
   CASE (   10)
     GO TO   700
   CASE (   11)
     GO TO 800
   CASE (   12)
     GO TO   440
   CASE (   13)
     GO TO   600
   CASE (   14)
     GO TO  1120
   CASE (   15)
     GO TO  1700
   CASE (   16)
     GO TO  1100
   CASE (   17)
     GO TO  1720
   CASE (   18)
     GO TO   520
   CASE (   19)
     GO TO  1200
 END SELECT
 
!           FIND  VIEW  AXES  MAXI  ORTH  PERS  STER  CONT  CAME  PAPE
!    1       PEN  BLAN  ORIG  SCAL  CSCA  PROJ  PTIT  OCUL  VANT
 
 
!     RECHECK IF PROJECTION CARD
 
 210 DO  kwrd = 5,7
   IF (word == hkey(kwrd)) GO TO 230
 END DO
 GO TO 1750
 
!     PROJECTION
 
 230 prject    = kwrd-4
 vangle(1) = 0.
 vangle(2) =-1.e10
 vangle(3) = 34.27
 fscale    = 1
 fvp = 1
 for = 1
 IF (nofind == 0) nofind = -1
 CALL rdword (mode,word)
 IF (word /= hkey(16)) GO TO 140
 
!     READ SECOND WORD OF ORTHO.,PERS.,OR STERO. SHOULD BE PROJECTION
 
 IF (org == 0) GO TO 140
 DO  i  = 1,org
   edge(i,1) = 0.
   edge(i,2) = 0.
   edge(i,3) = 1.
   edge(i,4) = 1.
 END DO
 org = 0
 GO TO 140
 
!     AXES
 
 250 DO  j = 1,3
   IF (mode == 0) CALL rdmode (*140,*260,*40,mode,word)
   260 CALL rdword (mode,word)
   DO  i = 1,7
     IF (word == axisd(i)) GO TO 280
   END DO
   GO TO 310
   280 axis(j) = i - 4
 END DO
 IF (mode == 0) CALL rdmode (*320,*300,*320,mode,word)
 300 CALL rdword (mode,word)
 310 IF (word == anti) GO TO 330
 320 k = 1
 GO TO 340
 330 k = -1
 340 DO  j = 1,3
   daxis(j) = k*axis(j)
 END DO
 IF (mode >= eor) GO TO 40
 IF (mode < 0 .OR. word == symm .OR. word == anti) GO TO 140
 GO TO 160
 
!     CAMERA
 
 400 ASSIGN 420 TO tra
 IF (mode <= 0) CALL rdmode (*1910,*410,*40,mode,word)
 410 CALL rdword (mode,word)
 n = 2
 IF (word == film) n = 1
 IF (word == pape) n = 2
 IF (word == both) n = 3
 IF (n == 0) THEN
   GO TO  1750
 ELSE
   GO TO   430
 END IF
 420 n = iwrd
 430 camera = n
 GO TO 140
 
!     BLANK FRAMES
 
 440 IF (mode == 0) GO TO 1750
 CALL rdword (mode,word)
 IF (word /= fram .OR. mode /= 0) GO TO 1750
 ASSIGN 450 TO tra
 GO TO 1900
 450 bframs = iwrd
 GO TO 140
 
!     MAXIMUM DEFORMATION
 
 500 IF (mode <= 0) GO TO 1750
 CALL rdword (mode,word)
 IF (word /= defo .OR. mode /= 0) GO TO 1750
 ASSIGN 510 TO tra
 GO TO 1940
 510 maxdef = fwrd
 GO TO 140
 
!     OCULAR SEPARATION
 
 520 IF (mode <= 0) GO TO 1750
 CALL rdword (mode,word)
 IF (word /= sepa .OR. mode /= 0) GO TO 1750
 ASSIGN 530 TO tra
 GO TO 1940
 530 s0s = fwrd
 GO TO 140
 
!     ORIGIN
 
 600 IF (mode /= 0) GO TO 1750
 ASSIGN 610 TO tra
 GO TO 1900
 
!     ORIGIN ID
 
 610 id = iwrd
 ASSIGN 620 TO tra
 GO TO 1940
 
!     HORIZONTAL LOCATION (LEFT EYE - STEREO)
 
 620 x = fwrd*cntsin
 ASSIGN 630 TO tra
 GO TO 1940
 
!     VERTICAL LOCATION
 
 630 y = fwrd*cntsin
 IF (org == 0) GO TO 670
 DO  j = 1,org
   IF (origin(j) == id) GO TO 680
 END DO
 IF (org < norg) GO TO 670
 IF (prnt <   0) GO TO 650
 ERR(1) = 1
 ERR(2) = norg
 CALL wrtprt (merr,ERR,msg1,nmsg1)
 650 org = norg
 DO  i = 1,2
   edge(org+1,i+0) = 0.
   edge(org+1,i+2) = 1.
 END DO
 670 org = org + 1
 j   = org
 origin(j) = id
 IF (nofind == 0) nofind = -1
 680 xy(j,1) = x
 xy(j,3) = y
 for = 0
 ASSIGN 690 TO tra
 GO TO 1940
 
!     HORIZONTAL LOCATION (RIGHT EYE - STEREO)
 
 690 xy(j,2) = fwrd*cntsin
 GO TO 140
 
!     PAPER SIZE, TYPE
 
 700 IF (mode <= 0) GO TO 1750
 CALL rdword (mode,word)
 IF (word == TYPE) GO TO 760
 IF (word /= size .OR. mode /= 0) GO TO 1750
 ASSIGN 710 TO tra
 GO TO 1940
 710 x = fwrd
 CALL rdmode (*730,*720,*40,mode,word)
 720 CALL rdword (mode,word)
 IF (word /= by .AND. word /= hx) GO TO 1750
 IF (mode /= 0) GO TO 1750
 730 ASSIGN 740 TO tra
 GO TO 1940
 740 papsiz(1) = x
 papsiz(2) = fwrd
 CALL pltset
 CALL rdmode (*140,*750,*40,mode,word)
 750 CALL rdword (mode,word)
 IF (word /= TYPE) GO TO 160
 
!     PAPER TYPE
 
 760 IF (mode == 0) GO TO 1750
 CALL rdword (mode,word)
 paptyp(1) = awrd(1)
 paptyp(2) = awrd(2)
 IF (mode > 0) THEN
   GO TO   700
 ELSE
   GO TO   140
 END IF
 
!     PEN SIZE / COLOR
 
 800 IF (mode /= 0) GO TO 1750
 ASSIGN 810 TO tra
 GO TO 1900
 810 IF (iwrd /= 1 .AND. iwrd <= npens) GO TO 820
 ERR(1) = 1
 ERR(2) = iwrd
 CALL wrtprt (merr,ERR,msg4,nmsg4)
 iwrd = 1
 820 id = iwrd
 830 CALL rdmode (*140,*840,*40,mode,word)
 840 CALL rdword (mode,word)
 IF (word == size) GO TO 850
 IF (word /= colo) GO TO 160
 IF (mode ==    0) GO TO 1750
 CALL rdword (mode,word)
 penclr(id,1) = awrd(1)
 penclr(id,2) = awrd(2)
 IF (mode < 0) THEN
   GO TO   140
 ELSE IF (mode == 0) THEN
   GO TO   830
 ELSE
   GO TO   840
 END IF
 
!     PEN SIZE
 
 850 IF (mode /= 0) GO TO 1750
 ASSIGN 860 TO tra
 GO TO 1900
 860 pensiz(id) = iwrd
 GO TO 830
 
!     PLOTTER
 
 900 IF (mode == 0) GO TO 1750
 CALL rdword (mode,word)
 pltnam(1) = awrd(1)
 pltnam(2) = awrd(2)
 pltmod(1) = 0
 pltmod(2) = 0
 camera = 2
 fscale = 1
 fvp = 1
 for = 1
 IF (org == 0) GO TO 920
 DO  i  = 1,org
   edge(i,1) = 0.
   edge(i,2) = 0.
   edge(i,3) = 1.
   edge(i,4) = 1.
 END DO
 org = 0
 
!     CHECK FOR A MODEL NUMBER
 
 920 ASSIGN 960 TO tra
 j = 1
 IF (mode <= 0) CALL rdmode (*1910,*930,*970,mode,word)
 930 CALL rdword (mode,word)
 IF (word == dens ) GO TO 970
 IF (word /= hmode) GO TO 960
 940 IF (mode <= 0) CALL rdmode (*1910,*950,*970,mode,word)
 950 CALL rdword (mode,word)
 IF (word == dens) GO TO 970
 960 pltmod(j) = word
 j = j + 1
 IF (j == 2) GO TO 940
 970 CALL fndplt (id,n,pltmod)
 ploter = id
 model  = n
 CALL pltset
 IF (word == dens) GO TO 1000
 IF (mode >=  eor) GO TO 40
 
!     TAPE DENSITY ON PLOTTER CARD
 
 980 IF (mode <= 0) CALL rdmode (*980,*990,*40,mode,word)
 990 CALL rdword (mode,word)
 1000 IF (word /= dens) GO TO 160
 IF (mode /=    0) GO TO 140
 ASSIGN 1010 TO tra
 GO TO 1900
 1010 tapden = iwrd
 CALL rdmode (*140,*1020,*40,mode,word)
 1020 CALL rdword (mode,word)
 IF (word == bpi) GO TO 140
 GO TO 160
 
!     PROJECTION PLANE SEPARATION
 
 1100 IF (mode == 0) GO TO 1750
 CALL rdword (mode,word)
 
!     USER MAY HAVE REVERSE ENGLISH
 
 IF (word /= plan) GO TO 210
 IF (mode ==    0) GO TO 1750
 CALL rdword (mode,word)
 IF (mode /= 0 .OR. word /= sepa) GO TO 1750
 ASSIGN 1110 TO tra
 GO TO 1940
 1110 IF (prject == 2) d02 = fwrd
 IF (prject == 3) d03 = fwrd
 GO TO 140
 
!     SCALE
 
 1120 IF (mode /= 0) GO TO 1750
 ASSIGN 1130 TO tra
 GO TO 1940
 1130 IF (fwrd  == 0.) GO TO 1140
 IF (prject /= 3) scale(1) = cntsin*fwrd
 IF (prject == 3) scale(1) = cntin3*fwrd
 1140 fscale = 0
 ASSIGN 1150 TO tra
 GO TO 1940
 1150 IF (fwrd  /= 0.) scale(2) = fwrd
 IF (nofind == 0) nofind = -1
 GO TO 140
 
!     VANTAGE POINT
 
 1200 IF (mode == 0) GO TO 1750
 CALL rdword (mode,word)
 IF (word /= poin .OR. mode /= 0) GO TO 1750
 ASSIGN 1220 TO tra
 j = 0
 1210 j = j + 1
 IF (j == 3) j = 4
 IF (prject == 3 .AND. j == 6) j = 3
 GO TO 1940
 1220 vanpnt(j) = fwrd
 IF ((prject /= 3 .AND. j /= 5) .OR. (prject == 3 .AND. j /= 3)) GO TO 1210
 fvp =  0
 IF (nofind == 0) nofind = -1
 GO TO 140
 
!     VIEW
 
 1230 IF (mode /= 0) GO TO 1750
 ASSIGN 1250 TO tra
 j = 4
 1240 j = j - 1
 GO TO 1940
 1250 vangle(j) = fwrd
 IF (nofind == 0) nofind = -1
 IF (j-1 == 0) THEN
   GO TO   140
 ELSE
   GO TO  1240
 END IF
 
!     CONTOUR
 
!     RESTORE DEFAULTS
 
 1300 icntvl = 1
 ncntr  = 10
 color  = 0
 layer  = 0
 where  = 1
 DIRECT = 2
 cntr(1)= 0.0
 cntr(2)= 0.0
 
!     FLAG AND LASSET SET IN PLOT AND CONPLT
 
 1310 IF (mode <= 0) CALL rdmode (*1310,*1320,*40,mode,word)
 1320 CALL rdword (mode,word)
 IF (word == colo .OR. word == fill .OR. word == laye) GO TO 1340
 IF (word /= even) GO TO 1370
 ASSIGN 1330 TO tra
 GO TO 1900
 1330 ncntr = MIN0 (50,iwrd)
 GO TO 1310
 1340 IF (word == colo) ASSIGN 1350 TO tra
 IF (word == fill) ASSIGN 1360 TO tra
 IF (word == laye) GO TO 1600
 GO TO 1900
 1350 color = iwrd
 GO TO 1310
 1360 color = -iwrd
 GO TO 1310
 
 1370 IF (word /= list) GO TO 1500
 IF (mode >    0) GO TO 1580
 ncntr = 0
 ASSIGN 1390 TO tra
 1380 CALL rdmode (*1950,*1320,*40,mode,word)
 1390 IF (ncntr < 50) GO TO 1400
 IF (prnt  <  0) GO TO 1380
 ERR(1) = 1
 ERR(2) = iwrd
 CALL wrtprt (merr,ERR,msg5,nmsg5)
 GO TO 1380
 1400 ncntr = ncntr + 1
 cntr(ncntr) = fwrd
 GO TO 1380
 
 1500 IF (word == z1  ) GO TO 1510
 IF (word == z2  ) GO TO 1520
 IF (word == MAX ) GO TO 1530
 IF (word == mid ) GO TO 1540
 IF (word == comm) GO TO 1550
 IF (word == disp) GO TO 1310
 IF (word == stre) GO TO 1310
 IF (word /= loca) GO TO 1560
 DIRECT = 1
 GO TO 1310
 1510 where  = 1
 GO TO 1310
 1520 where  =-1
 GO TO 1310
 1530 where  = 2
 GO TO 1310
 1540 where  = 3
 GO TO 1310
 1550 DIRECT = 2
 GO TO 1310
 
 1560 DO  j = 1,20
   IF (word == icnda(j)) GO TO 1590
 END DO
 1580 IF (prnt < 0) GO TO 1310
 ERR(1) = 2
 ERR(2) = awrd(1)
 ERR(3) = awrd(2)
 CALL wrtprt (merr,ERR,msg2,nmsg2)
 GO TO 1310
 
 1590 icntvl = j
 
!     SET STRESS FILE TO STRAIN FILE
 
 IF (icntvl == 20) oesx = onrgy1
 GO TO 1310
 
!     ASSIGN LAYER NUMBER HERE FOR COMPOSITS
 
 1600 ASSIGN 1610 TO tra
 
!     SET STRESS FILE TO LAYER STRESS
 
 oesx = oes1l
 GO TO 1900
 1610 layer = iwrd
 GO TO 1310
 
!     CSCALE
 
 1700 IF (mode /= 0) GO TO 1750
 ASSIGN 1710 TO tra
 GO TO 1940
 1710 chrscl = fwrd
 IF (nofind ==   0) nofind = -1
 IF (chrscl < 1.0) chrscl = 1.0
 CALL pltset
 GO TO 140
 
!     PTITLE
 
 1720 fpltit = 1
 DO  i = 1,17
   pltitl(i) = blank4
 END DO
 j = color
 DO  i = 1,17,2
   CALL rdword (mode,word)
   pltitl(i  ) = awrd(1)
   pltitl(i+1) = awrd(2)
   IF (mode == 0) GO TO 140
 END DO
 color = j
 IF (mode /= 0) CALL rdword (mode,word)
 GO TO 140
 
!     UNRECOGNIZABLE PLOT PARAMETER.
 
 1750 IF (prnt < 0) GO TO 140
 ERR(1) = 2
 ERR(2) = awrd(1)
 ERR(3) = awrd(2)
 CALL wrtprt (merr,ERR,msg2,nmsg2)
 GO TO 140
 
!     END OF PLOT INPUT
 
 1800 IF (prnt >= 0) GO TO 1820
 DO  i = 1,96
   title(i) = savtit(i)
 END DO
 1820 CONTINUE
 RETURN
 
 
!     READ AN INTEGER ON A PARAMETER CARD
 
 1900 CALL rdmode (*1910,*140,*40,mode,word)
 1910 IF (mode == -1) GO TO 1930
 IF (mode == -4) GO TO 1920
 iwrd = fwrd
 GO TO 1930
 1920 iwrd = dwrd
 1930 GO TO tra, (420,450,610,810,860,1330,1350,1360,1610,960,1010)
 
!     READ A DECIMAL NUMBER ON A PARAMETER CARD
 
 1940 CALL rdmode (*1950,*140,*40,mode,word)
 1950 IF (mode == -4) GO TO 1960
 IF (mode /= -1) GO TO 1970
 fwrd = iwrd
 GO TO 1970
 1960 fwrd = dwrd
 1970 GO TO tra, ( 510, 530, 620, 630, 690, 710, 740,1110,1130,1150,  &
     1220,1250,1390,1710)
 
END SUBROUTINE param
