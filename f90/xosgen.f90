SUBROUTINE xosgen
     
!     THE PURPOSE OF THIS ROUTINE IS TO GENERATE THE OSCAR ARRAY.
 
!          ... DESCRIPTION OF PROGRAM VARIABLES ...
!     IENDF  = FLAG SIGNALING END OF DMAP SEQUENCE.
!     LDEF   = SCRATCH USED IN SCANNING LBLTBL TABLE.
!     LBLTOP = TOP OF LBLTBL ARRAY.
!     LBLBOT = BOTTOM OF LBLTBL ARRAY.
!     LSTLBL = POINTER TO LAST LABEL ENTRY MADE IN LBLTBL.
!     LSTPAR = POINTER TO LAST PARAMETER NAME ENTRY MADE IN LBLTBL.
!     NAMTBL = NAME CONVERSION TABLE FOR TYPE E NAMES.
!     IEXFLG = FLAG INDICATING LAST OSCAR ENTRY WAS EXIT.
!     IOSPNT = POINTER TO NEXT AVAILABLE WORD IN OSCAR ENTRY.
!     NOSPNT = POINTER TO DATA BLOCK NAME COUNT IN OSCAR ENTRY.
!     NTYPEE = TABLE CONTAINING TYPE E DMAP NAMES
!     IPRCFO = POINTER TO LAST TYPE F OR O OSCAR ENTRY.
!     NDIAG1 = NAME OF THE DIAGNOSTIC O/P PROCESSOR
!     ITYPE  = TABLE FOR TRANSLATING TYPE CODES TO WORD LENGTH
!     VARFLG = FLAG INDICATING VARIABLE FOUND IN EQUIV OR PURGE
!              INSTRUCTION.
 
 IMPLICIT INTEGER (a-z)
 EXTERNAL        lshift,rshift,andf,orf
 LOGICAL :: skip
 DIMENSION       prechk(2),xdmap(2),declar(3),fparam(3),  &
     dmpcrd(1),nskip(5,2),cdcomp(3),namtbl(12),  &
     itype(6),med(1),lbltbl(1),oscar(1),os(5)
 COMMON /xfiat / fiat(3)
 COMMON /system/ bufsz,optape,nogo,dum(20),icfiat,junk(54), iswtch(3),icpflg
 COMMON /moddmp/ iflg(6),namopt(26)
 COMMON /xgpic / icold,islsh,iequl,nblank,nxequi,  &
     nmed,nsol,ndmap,nestm1,nestm2,nexit,  &
     nbegin,nend,njump,ncond,nrept,ntypee(9),  &
     maskhi,masklo,isgnon,nosgn,iallon,masks(1)
 COMMON /zzzzzz/ core(1)
 COMMON /xgpi2 / lmpl,mplpnt,mpl(1)
!WKBR COMMON /XGPI3 / PVT(2)
 COMMON /xgpi3 / pvt(200)
 COMMON /xgpi4 / irturn,insert,iseqn,dmpcnt,  &
     idmpnt,dmppnt,bcdcnt,length,icrdtp,ICHAR,newcrd, modidx,ldmap,isavdw,dmap(1)
 COMMON /xgpi5 / iapp,start,alter(2),sol,subset,iflag,iestim,  &
     icftop,icfpnt,lctlfl,ictlfl(1)
 COMMON /xgpi6 / medtp,fnmtp,cnmtp,medpnt,lmed,iplus,diag14,diag17,  &
     diag4,diag25,ifirst,ibuff(20)
 COMMON /xgpi7 / fpnt,lfile,FILE(1)
 COMMON /xgpid / icst,iunst,imst,ihapp,idsapp,idmapp,  &
     isave,itape,iappnd,intgr,losgn
 COMMON /xgpie / nscr
 COMMON /xvps  / vps(2)
!WKBR COMMON /XCEITB/ CEITBL(2)
 COMMON /xceitb/ ceitbl(42)
 COMMON /xoldpt/ xx(4),seqno
 COMMON /autocm/ preflg,nnames,prenam(100)
 COMMON /autosm/ nwords,savnam(100)
 COMMON /passer/ istopf,modnam
 
!     EQUIVALENCE     (NTYPEE(1),NTIME ), (NTYPEE(2),NSAVE )
!    1                (NTYPEE(3),NOUTPT), (NTYPEE(4),NCHKPT)
!    2                (NTYPEE(5),NPURGE), (NTYPEE(6),NEQUIV)
!    3                (NTYPEE(7),NCPW  ), (NTYPEE(8),NBPC  )
!    4                (NTYPEE(9),NWPC  )
 EQUIVALENCE     (namtbl(9),nxpurg)
 EQUIVALENCE     (oscar (1),dmpcrd(1),lbltbl(1),med(1),os(5)),  &
     (core(1),os(1),loscar), (os(2),osprc), (os(3),osbot), (os(4),ospnt)
 
 DATA    xchk  / 4HXCHK/
 DATA    itype / 1,1,2,2,2,4/
 DATA    iprcfo/ 0     /, iendf / 0/
 DATA    nfile / 4HFILE/
 DATA    nvps  / 4HVPS /
 DATA    prechk/ 4HPREC,  4HHK  /,xdmap / 4HXDMA, 4HP     /
 DATA    nceit1/ 4HCEIT/, nceit2/ 4HBL  /
 DATA    nlblt1/ 4HLBLT/, nlblt2/ 4HBL  /
 DATA    declar/ 4HBEGI,  4HLABE, 4HFILE/
 DATA    fparam/ 4HTAPE,  4HAPPE, 4HSAVE/
 DATA    namtbl/ 4HXTIM,  4HE   , 4HXSAV, 4HE   , 4HXUOP, 4H    ,  &
     4HXCHK,  4H    , 4HXPUR, 4HGE  , 4HXEQU, 4HIV  /
 DATA    nskip / 10*0 /, cdcomp / 4HCOMP, 4HON  , 4HOFF   /
 
!     INITIALIZE
 
 ifirst = 0
 osbot  = 1
 nwords = 0
 lookup = 0
 preflg = 0
 ivrept = 0
 ilevel = 0
 skip   =.false.
 ospnt  = osbot
 oscar(osbot  ) = 0
 oscar(osbot+1) = 1
 
!     FOR RESTART ALLOW CHECKPOINT AND JUMP ENTRIES TO BE INSERTED IN
!     OSCAR BY XGPI.
 
 IF (start == icst) GO TO 10
 oscar(osbot+1) = 3
 
!     ALLOCATE 50 WORDS IN OPEN CORE FOR LBLTBL AND SET LBLTBL
!     PARAMETERS.
 
 10 lblbot = loscar
 lbltop = loscar - 50
 loscar = lbltop - 1
 lstlbl = lbltop - 4
 lstpar = lblbot + 1
 
!     INITIALIZE DMPCRD ARRAY FOR RIGID FORMAT
 
 icrdtp = loscar
 
!     ****************************************
!     PREPARE TO PROCESS NEXT DMAP INSTRUCTION
!     ****************************************
 
 100 dmpcnt = dmpcnt + 1
 IF (iapp == idmapp) GO TO 110
 medpnt = med(medtp+1)*(dmpcnt - 1) + medtp + 2
 IF (med(medtp) < dmpcnt .AND. iapp /= idmapp) GO TO 2390
 110 newcrd =-1
 insert = 0
 
!     SEE IF DMAP INSTRUCTION IS TO BE DELETED OR INSERTED
 
 IF (alter(1) == 0 .OR. alter(1) > dmpcnt) GO TO 130
 IF (alter(1) <= dmpcnt .AND. alter(2) >= dmpcnt) GO TO 150
 IF (alter(2) == 0) GO TO 120
 
!     JUST FINISHED DELETING, SET INSERT AND ALTER FOR INSERTING
 
 alter(1) = alter(2)
 alter(2) = 0
 120 IF (alter(1) /= dmpcnt-1) GO TO 130
 insert = 1
 dmpcnt = dmpcnt - 1
 GO TO 160
 
!     GET NEXT DMAP INSTRUCTION
!     FOR RIGID FORMAT SEE IF OSCAR ENTRY IS PART OF SUBSET
 
 130 IF (iapp == idmapp) GO TO 160
 i = med(medtp+1)
 DO  j = 1,i
   k = medpnt + j -1
   IF (med(k) /= 0) GO TO 160
 END DO
 
!     SET INSERT FLAG TO NO PRINT
 
 150 insert = -2
 GO TO 310
 
!     CHECK FOR CONDITIONAL COMPILATION END
 
 160 IF (ilevel <= 0) GO TO 190
 DO  i = 1,ilevel
   IF (IABS(nskip(i,1)) < 99999) nskip(i,1) = nskip(i,1) - 1
 END DO
 IF (nskip(ilevel,1) == -1) GO TO 180
 IF (skip) insert = insert - 2
 GO TO 190
 180 skip   =.false.
 ilevel = ilevel - 1
 
 190 IF (lookup /= 1 .OR. preflg == 0) GO TO 200
 preflg = -preflg
 CALL autock (ospnt)
 200 modnam = 1
 lookup = 0
 CALL xscndm
 modnam = 0
 SELECT CASE ( irturn )
   CASE (    1)
     GO TO 2120
   CASE (    2)
     GO TO 210
   CASE (    3)
     GO TO 2120
   CASE (    4)
     GO TO 100
   CASE (    5)
     GO TO 2060
 END SELECT
 210 IF (.NOT.skip) GO TO 220
 
!     CHECK LABELS EVEN IF CONDITIONAL COMPILATION
 
 IF (dmap(dmppnt) == declar(2)) GO TO 1270
 GO TO 310
 
!     FIND MPL ENTRY AND BRANCH ON TYPE
 
 220 mplpnt = 1
 modidx = 1
 IF (dmap(dmppnt) == prechk(1) .AND. dmap(dmppnt+1) == prechk(2)) GO TO 1490
 IF (dmap(dmppnt) == xdmap(1)  .AND. dmap(dmppnt+1) == xdmap(2)) GO TO 1570
 IF (dmap(dmppnt) == cdcomp(1) .AND. (dmap(dmppnt+1) == cdcomp(2)  &
     .OR. dmap(dmppnt+1) == cdcomp(3))) GO TO 1740
 230 IF (mpl(mplpnt+1) == dmap(dmppnt) .AND. mpl(mplpnt+2) ==  &
     dmap(dmppnt+1)) GO TO 240
 
!     CHECK FOR ERROR IN MPL TABLE
 
 IF (mpl(mplpnt) < 1 .OR. mpl(mplpnt) > lmpl) GO TO 2140
 mplpnt = mpl(mplpnt) + mplpnt
 modidx = 1 + modidx
 IF (mplpnt-lmpl < 0) THEN
   GO TO   230
 ELSE
   GO TO  2130
 END IF
 
!     GET FORMAT TYPE FROM MPL AND BRANCH
 
 240 i = mpl(mplpnt + 3)
 IF (i < 1 .OR. i > 5) GO TO 2140
 SELECT CASE ( i )
   CASE (    1)
     GO TO 400
   CASE (    2)
     GO TO 400
   CASE (    3)
     GO TO 500
   CASE (    4)
     GO TO 800
   CASE (    5)
     GO TO 1200
 END SELECT
 
!     *****************************************************
!     RETURN HERE AFTER DMAP INSTRUCTION HAS BEEN PROCESSED
!     *****************************************************
 
!     CHECK FOR FATAL ERROR
 
 300 IF (nogo == 2) GO TO 2060
 
!     CHECK FOR END OF DMAP SEQUENCE.
 
 IF (iendf /= 0) GO TO 1900
 
!     CHECK FOR $ ENTRY IN DMAP AND GET NEXT DMAP INSTRUCTION
 
 310 CALL xscndm
 SELECT CASE ( irturn )
   CASE (    1)
     GO TO 320
   CASE (    2)
     GO TO 320
   CASE (    3)
     GO TO 320
   CASE (    4)
     GO TO 100
   CASE (    5)
     GO TO 2060
 END SELECT
 320 IF (nogo == 0 .AND. insert >= 0) GO TO 2160
 GO TO 310
 
!     ********************************************
!     GENERATE OSCAR ENTRY WITH TYPE F OR O FORMAT
!     ********************************************
 
!     GENERATE LINK HEADER SECTION
 
 400 CALL xlnkhd
 iprcfo = ospnt
 
!     GENERATE I/P FILE SECTION
 
 CALL xipfl
 SELECT CASE ( irturn )
   CASE (    1)
     GO TO 410
   CASE (    2)
     GO TO 2100
 END SELECT
 
!     SAVE POINTER TO O/P FILE SECTION
 
 410 j = ospnt + oscar(ospnt)
 
!     GENERATE O/P FILE SECTION
 
 CALL xopfl
 SELECT CASE ( irturn )
   CASE (    1)
     GO TO 420
   CASE (    2)
     GO TO 2110
 END SELECT
 
!     NUMBER OF SCRATCH FILES TO OSCAR
 
 420 i = ospnt + oscar(ospnt)
 oscar(i) = mpl(mplpnt)
 
!     INCREMENT OSCAR WORD COUNT AND MPLPNT
 
 oscar(ospnt) = 1 + oscar(ospnt)
 mplpnt = 1 + mplpnt
 
!     GENERATE PARAMETER SECTION
 
 CALL xparam
 SELECT CASE ( irturn )
   CASE (    1)
     GO TO 430
   CASE (    2)
     GO TO 2060
 END SELECT
 
!     CONTINUE COMPILATION
!     ZERO INTERNAL CHECKPOINT FLAG IN OSCAR ENTRY FOR TYPE F ENTRY
 
 430 IF (andf(oscar(ospnt+2),maskhi) == 2) GO TO 440
 i = ospnt + oscar(ospnt)
 oscar(i) = 0
 oscar(ospnt) = 1 + oscar(ospnt)
 440 CONTINUE
 IF (nwords == 0) GO TO 450
 CALL autosv
 nwords = 0
 450 IF (preflg == 0 .OR. istopf == 0) GO TO 460
 CALL autock (istopf)
 460 CONTINUE
 GO TO 300
 
!     ***************************************
!     GENERATE OSCAR ENTRY WITH TYPE C FORMAT
!     ***************************************
 
!     GENERATE LINK HEADER SECTION
 
 500 CALL xlnkhd
 
!     UPDATE OSCAR ENTRY WORD COUNT TO INCLUDE VALUE SECTION.
 
 oscar(ospnt) = 7
 
!     CHECK FOR END CARD
 
 IF (oscar(ospnt+3) /= nend) GO TO 510
 oscar(ospnt+3) = nexit
 iendf = 1
 
!     SET EXECUTE FLAG IN OSCAR FOR END
 
 oscar(ospnt+5) = orf(isgnon,oscar(ospnt+5))
 
!     GET NEXT ENTRY IN DMAP
 
 510 CALL xscndm
 SELECT CASE ( irturn )
   CASE (    1)
     GO TO 2160
   CASE (    2)
     GO TO 520
   CASE (    3)
     GO TO 630
   CASE (    4)
     GO TO 630
   CASE (    5)
     GO TO 2060
 END SELECT
 
!     IF NEXT DMAP ENTRY IS BCD IT SHOULD BE LABEL NAME FOR BRANCH
!     DMAP INSTRUCTION.
 
 520 IF (oscar(ospnt+3) == nexit) GO TO 2160
 
!     SEARCH LABEL TABLE FOR LABEL NAME
 
 IF (lstlbl < lbltop) GO TO 540
 DO  j = lbltop,lstlbl,4
   IF (dmap(dmppnt) == lbltbl(j) .AND. dmap(dmppnt+1) == lbltbl(j+1))  &
       GO TO 550
 END DO
 
!     NAME NOT FOUND IN TABLE
 
 540 ldef = 0
 GO TO 560
 
!     NOW SEE IF LABEL HAS BEEN REFERENCED
 
 550 IF (lbltbl(j+3) == 0) GO TO 580
 ldef = lbltbl(j+2)
 
!     MAKE NEW ENTRY IN LABEL TABLE, CHECK FOR TABLE OVERFLOW
 
 560 ASSIGN 570 TO irturn
 IF (lstlbl+8 >= lstpar) GO TO 2220
 570 lstlbl = lstlbl + 4
 j = lstlbl
 lbltbl(j  ) = dmap(dmppnt  )
 lbltbl(j+1) = dmap(dmppnt+1)
 lbltbl(j+2) = ldef
 580 lbltbl(j+3) = ospnt
 
!     GET NEXT ENTRY FROM DMAP, ENTRY IS $ FOR JUMP,NAME FOR COND,
!     VALUE FOR REPT.
 
 CALL xscndm
 SELECT CASE ( irturn )
   CASE (    1)
     GO TO 2160
   CASE (    2)
     GO TO 600
   CASE (    3)
     GO TO 720
   CASE (    4)
     GO TO 590
   CASE (    5)
     GO TO 2060
 END SELECT
 
!     DMAP INSTRUCTION IS JUMP
 
 590 oscar(ospnt+6) = 0
 IF (oscar(ospnt+3) == njump) GO TO 300
 GO TO 2160
 
!     COND DMAP INSTRUCTION, ENTER PARAMETER NAME IN LABEL TABLE.
 
 600 IF (oscar(ospnt+3) /= nrept) GO TO 610
 ivrept =  1
 GO TO 640
 610 IF (oscar(ospnt+3) /= ncond) GO TO 2160
 ASSIGN 620 TO irturn
 IF (lstpar-8 <= lstlbl) GO TO 2220
 620 lstpar = lstpar - 4
 lbltbl(lstpar  ) = dmap(dmppnt  )
 lbltbl(lstpar+1) = dmap(dmppnt+1)
 lbltbl(lstpar+2) = ospnt + 6
 lbltbl(lstpar+3) = ospnt
 GO TO 300
 
!     EXIT DMAP INSTRUCTION, SET EXECUTE FLAG AND OSCAR VALUE SECTION.
 
 630 IF (oscar(ospnt+3) /= nexit) GO TO 2160
 IF (dmap(dmppnt) /= intgr) dmap(dmppnt+1) = 0
 dmap(dmppnt  ) = intgr
 dmap(dmppnt+2) = rshift(iallon,1)
 
!     ENTER LOOP COUNT IN CEITBL FOR REPT AND EXIT INSTRUCTIONS
 
 640 ceitbl(2) = ceitbl(2) + 4
 IF (ceitbl(2) > ceitbl(1)) GO TO 2280
 
!     I = POINTER TO LOOP COUNT IN CEITBL ENTRY
 
 i = ceitbl(2) - 2
 IF (ivrept == 0) GO TO 700
 
!     PROCESS VARIABLE REPT INSTRUCTION - FIND PARAM IN VPS
 
 kdh = 3
 650 IF (dmap(dmppnt) == vps(kdh) .AND. dmap(dmppnt+1) == vps(kdh+1))  &
     GO TO 660
 kdh = kdh + andf(vps(kdh+2),maskhi) + 3
 IF (kdh - vps(2) < 0) THEN
   GO TO   650
 ELSE
   GO TO   670
 END IF
 
!     PARAMETER FOUND
 
 660 IF (andf(rshift(vps(kdh+2),16),15) /= 1) GO TO 2210
 ceitbl(i) = lshift(kdh,16)
 ceitbl(i) = orf(ceitbl(i),isgnon)
 GO TO 710
 
!     CHECK PVT FOR PARAMETER
 
 670 kdh = 3
 680 length = andf(pvt(kdh+2),nosgn)
 length = itype(length)
 IF (dmap(dmppnt) == pvt(kdh) .AND. dmap(dmppnt+1) == pvt(kdh+1)) GO TO 690
 kdh = kdh + length + 3
 IF (kdh - pvt(2) < 0) THEN
   GO TO   680
 ELSE
   GO TO  2200
 END IF
 690 IF (length /= itype(1)) GO TO 2210
 ceitbl(i) = lshift(pvt(kdh+3),16)
 GO TO 710
 700 ceitbl(i) = lshift(dmap(dmppnt+1),16)
 
!     FIRST WORD OF CEITBL ENTRY CONTAINS OSCAR RECORD NUMBERS OF
!     BEGINNING AND END OF LOOP
 
 710 ceitbl(i-1) = iseqn
 ivrept = 0
 
!     OSCAR VALUE SECTION CONTAINS POINTER TO LOOP COUNT IN CEITBL ENTRY
 
 oscar(ospnt+6) = i
 GO TO 300
 
!     REPT DMAP INSTRUCTION, COUNT TO VALUE SECTION.
 
 720 IF (oscar(ospnt+3) == nrept) GO TO 640
 GO TO 2160
 
!     ***************************************
!     GENERATE OSCAR ENTRY WITH TYPE E FORMAT
!     ***************************************
 
!     PREFIX MODULE NAME WITH AN X
 
 800 DO  i = 1,6
   IF (ntypee(i) == dmap(dmppnt)) EXIT
 END DO
 820 i = 2*i - 1
 dmap(dmppnt  ) = namtbl(i  )
 dmap(dmppnt+1) = namtbl(i+1)
 
!     GENERATE LINK HEADER FOR OSCAR
 
 IF (i == 9 .OR. i == 11) lookup = 1
 os2b4 = osprc
 CALL xlnkhd
 
!     BRANCH ON DMAP NAME AND GENERATE VALUE/OUTPUT SECTION OF OSCAR
 
 i = (i+1)/2
 SELECT CASE ( i )
   CASE (    1)
     GO TO 830
   CASE (    2)
     GO TO 860
   CASE (    3)
     GO TO 990
   CASE (    4)
     GO TO 990
   CASE (    5)
     GO TO 990
   CASE (    6)
     GO TO 990
 END SELECT
 
!     EXTIME ENTRY, CHECK ESTIM IN CONTROL FILE
 
 830 oscar(ospnt+5) = andf(oscar(ospnt+5),nosgn)
 IF (iestim == 0) GO TO 300
 
!     GET TIME SEGMENT NAME
 
 CALL xscndm
 SELECT CASE ( irturn )
   CASE (    1)
     GO TO 2370
   CASE (    2)
     GO TO 840
   CASE (    3)
     GO TO 2370
   CASE (    4)
     GO TO 2370
   CASE (    5)
     GO TO 2060
 END SELECT
 840 i = iestim + ictlfl(iestim) - 1
 j = iestim + 1
 DO  k = j,i,2
   IF (dmap(dmppnt) == ictlfl(k) .AND. dmap(dmppnt+1) == ictlfl(k+1))  &
       oscar(ospnt+5) = orf(oscar(ospnt+5),isgnon)
 END DO
 GO TO 300
 
!     XSAVE ENTRY, ENTER POINTERS IN VALUE SECTION OF OSCAR.
 
 860 i = ospnt + oscar(ospnt)
 oscar(i) = 0
 k = i - 1
 
!     GET PARAMETER NAME FROM DMAP.
 
 870 CALL xscndm
 SELECT CASE ( irturn )
   CASE (    1)
     GO TO 2260
   CASE (    2)
     GO TO 880
   CASE (    3)
     GO TO 2260
   CASE (    4)
     GO TO 930
   CASE (    5)
     GO TO 2060
 END SELECT
 
!     FIND PARAMETER IN VPS AND ENTER POINTER TO VALUE IN OSCAR.
 
 880 k = k + 2
 oscar(i  ) = oscar(i) + 1
 oscar(k  ) = 0
 oscar(k+1) = 0
 j = 3
 890 IF (vps(j) == dmap(dmppnt) .AND. vps(j+1) == dmap(dmppnt+1)) GO TO 900
 l = andf(vps(j+2),maskhi)
 j = j + l + 3
 IF (j < vps(2)) GO TO 890
 
!     PARAMETER NOT IN VPS - ERROR
 
 GO TO 2270
 
!     PARAMETER FOUND IN VPS
 
 900 oscar(k) = j + 3
 
!     SEE IF PARAMETER WAS ALREADY SAVED
 
 j  = i + 1
 j1 = k - 2
 IF (j1 < j) GO TO 870
 DO  l = j,j1,2
   IF (oscar(l) == oscar(k)) GO TO 920
 END DO
 GO TO 870
 
!     PARAMETER DUPLICATED
 
 920 k = k - 2
 oscar(i) = oscar(i) - 1
 GO TO 2150
 
 
!     END OF SAVE PARAMETER NAME LIST, INCREMENT OSCAR WORD COUNT.
 
 930 oscar(ospnt) = oscar(ospnt) + 2*oscar(i) + 1
 
!     GET PARAMETER VALUE DISPLACEMENT IN COMMON FROM PRECEDING
!     OSCAR ENTRY.
 
 iosdav = osprc
 IF (oscar(osprc+3) == xchk) osprc = os2b4
 IF (andf(oscar(osprc+2),maskhi) > 2) GO TO 2420
 
!     J = OSCAR POINTER TO BEGINNING OF PARAMETER SECTION.
 
 j = osprc + 6 + 3*oscar(osprc+6) + 1
 IF (andf(oscar(osprc+2),maskhi) == 1) j = j + 1 + 3*oscar(j)
 j = j + 1
 
!     N1 = PARAMETER COUNT,N2=PARAMETER DISPLACEMENT IN COMMON,
!     N3 = OSCAR POINTER TO PARAMETER ENTRIES IN PRECEDING OSCAR ENTRY.
 
 n3 = j + 1
 n1 = oscar(j)
 n2 = 1
 
!     SCAN PARAMETER LIST OF PRECEDING OSCAR ENTRY
 
 DO  m = 1,n1
   l = andf(oscar(n3),nosgn)
   IF (oscar(n3) > 0) GO TO 970
   n3 = n3 + 1
   
!     VARIABLE PARAMETER, COMPARE VPS POINTER WITH XSAVE VPS POINTERS.
   
   i1 = i + 1
   DO  k1 = i1,k,2
     IF (oscar(k1) == l) GO TO 950
   END DO
   GO TO 960
   950 oscar(k1+1) = n2
   960 l = andf(vps(l-1),maskhi)
   GO TO 980
   
!     CONSTANT PARAMETER, INCREMENT N2, N3
   
   970 n3 = n3 + l + 1
   980 n2 = n2 + l
 END DO
 
!     PARAMETER SECTION SCANNED, CHECK EXSAVE PARAMETER LIST FOR
!     PARAMETERS NOT FOUND IN PRECEDING OSCAR.
 
 GO TO 2290
 
!     XUOP,XCHK,XPURGE,OR XEQUIV OSCAR ENTRY - GENERATE FILE NAME LIST.
 
 990 nospnt = ospnt + oscar(ospnt)
 iprime = 1
 iospnt = nospnt + 1
 oscar(nospnt) = 0
 
!     GET NEXT ENTRY FROM DMAP CARD
 
 1000 CALL xscndm
 SELECT CASE ( irturn )
   CASE (    1)
     GO TO 1040
   CASE (    2)
     GO TO 1010
   CASE (    3)
     GO TO 2160
   CASE (    4)
     GO TO 1080
   CASE (    5)
     GO TO 2060
 END SELECT
 
!     DMAP ENTRY IS DATA BLOCK NAME, STORE IN OSCAR
 
 1010 oscar(iospnt  ) = dmap(dmppnt  )
 oscar(iospnt+1) = dmap(dmppnt+1)
 
!     MAKE SURE FILE IS NOT BLANK
 
 IF (oscar(iospnt) == nblank) GO TO 1000
 
!     FOR CHKPNT - MAKE SURE FILE IS NOT OUTPUT BY USER I/P PROCESSOR
 
 IF (oscar(ospnt+3) /= namtbl(7)) GO TO 1030
 m = fiat(3)*icfiat - 2
 DO  j = 4,m,icfiat
   IF (oscar(iospnt) == fiat(j+1) .AND. oscar(iospnt+1) == fiat(j+2))  &
       GO TO 2400
 END DO
 1030 iospnt = iospnt + 2
 oscar(nospnt) = 1 + oscar(nospnt)
 
!     INSERT EXTRA WORD INTO OSCAR FOR EACH PRIMARY DATA BLOCK IN
!     EQUIV STATEMENT
 
 IF (oscar(ospnt+3) /= namtbl(11) .OR. oscar(ospnt+4) /= namtbl(12)  &
     ) GO TO 1000
 IF (iprime == 0) GO TO 1000
 oscar(iospnt) = 0
 iospnt = iospnt + 1
 iprime = 0
 GO TO 1000
 
!     DMAP ENTRY IS OPERATOR, CHECK FOR / OPERATOR
 
 1040 IF ((dmap(dmppnt+1) /= islsh) .OR. (oscar(ospnt+3) /= nxequi .AND.  &
     oscar(ospnt+3) /= nxpurg)) GO TO 2160
 
!     OSCAR ENTRY IS XEQUIV OR XPURGE
 
 varflg = 0
 IF (oscar(ospnt+3) == nxpurg) GO TO 1050
 IF (oscar(nospnt)  < 2     ) GO TO 2160
 
!     GET PARAMETER NAME AND ENTER INTO LBLTBL
 
 1050 CALL xscndm
 SELECT CASE ( irturn )
   CASE (    1)
     GO TO 1110
   CASE (    2)
     GO TO 1060
   CASE (    3)
     GO TO 2160
   CASE (    4)
     GO TO 2160
   CASE (    5)
     GO TO 2060
 END SELECT
 1060 varflg = 1
 IF (dmap(dmppnt) == nblank) GO TO 1100
 ASSIGN 1070 TO irturn
 IF (lstpar-8 <= lstlbl) GO TO 2220
 1070 lstpar = lstpar - 4
 lbltbl(lstpar  ) = dmap(dmppnt  )
 lbltbl(lstpar+1) = dmap(dmppnt+1)
 lbltbl(lstpar+2) = iospnt
 lbltbl(lstpar+3) = ospnt
 idlhss = 2*oscar(nospnt)+oscar(ospnt) + 2
 IF (oscar(ospnt+3) == namtbl(11)) idlhss = idlhss + 1
 oscar(ospnt) = idlhss
 
!     CHECK FOR POSSIBILITY OF ANOTHER DATA BLOCK NAME LIST.
 
 CALL xscndm
 SELECT CASE ( irturn )
   CASE (    1)
     GO TO 990
   CASE (    2)
     GO TO 2160
   CASE (    3)
     GO TO 2160
   CASE (    4)
     GO TO 300
   CASE (    5)
     GO TO 2060
 END SELECT
 
!     END OF DMAP INSTRUCTION, INCREMENT OSCAR WORD COUNT IF NOT XEQUIV
!     OR XPURGE.
 
 1080 IF (oscar(ospnt+3) /= nxequi .AND. oscar(ospnt+3) /= nxpurg) GO TO 1090
 oscar(iospnt) = -1
 idlhss = 2*oscar(nospnt) + oscar(ospnt) + 2
 IF (oscar(ospnt+3) == namtbl(11)) idlhss = idlhss + 1
 oscar(ospnt) = idlhss
 GO TO 300
 1090 oscar(ospnt) = 2*oscar(nospnt) + oscar(ospnt) + 1
 
!     ELIMINATE ENTRY IF NOTHING CHECKPOINTED.
 
 IF (oscar(nospnt) == 0) osbot = osprc
 GO TO 300
 1100 CALL xscndm
 SELECT CASE ( irturn )
   CASE (    1)
     GO TO 1110
   CASE (    2)
     GO TO 2160
   CASE (    3)
     GO TO 2160
   CASE (    4)
     GO TO 2160
   CASE (    5)
     GO TO 2060
 END SELECT
 1110 IF ((dmap(dmppnt+1) /= islsh) .OR. (oscar(ospnt+3) /= nxequi .AND.  &
     oscar(ospnt+3) /= nxpurg)) GO TO 2160
 oscar(iospnt) = -1
 idlhss = 2*oscar(nospnt) + oscar(ospnt) + 2
 IF (oscar(ospnt+3) == namtbl(11)) idlhss = idlhss + 1
 oscar(ospnt) = idlhss
 GO TO 990
 
!     *******************************
!     DMAP INSTRUCTION IS DECLARATIVE
!     *******************************
 
!     PUT DUMMY ENTRY IN OSCAR FOR DIAGNOSTIC USE.
 
 1200 j = osbot  + oscar(osbot)
 oscar(j+3) = dmap(dmppnt)
 oscar(j+4) = dmap(dmppnt+1)
 oscar(j+5) = dmpcnt
 CALL xlnkhd
 
!     NOW PROCESS INSTRUCTION
 
 DO  j = 1,3
   IF (dmap(dmppnt) == declar(j)) THEN
      SELECT CASE ( j )
       CASE (    1)
         GO TO 1220
       CASE (    2)
         GO TO 1270
       CASE (    3)
         GO TO 1350
     END SELECT
   END IF
 END DO
 
!     BEGIN DECLARATIVE - PREPARE TO PROCESS NEXT DMAP INSTRUCTION
 
 1220 INDEX = 1
 1230 IF (ifirst > 0) GO TO 1250
 IF (diag14 == 0 .AND. diag17 == 0) GO TO 1250
 ifirst = 1
 CALL xgpimw (5,18,dmpcnt,ibuff)
 1240 IF (start /= icst) CALL xgpimw (10,0,0,0)
 1250 IF (INDEX >    1) GO TO 300
 1260 CALL xscndm
 SELECT CASE ( irturn )
   CASE (    1)
     GO TO 1260
   CASE (    2)
     GO TO 1260
   CASE (    3)
     GO TO 1260
   CASE (    4)
     GO TO 300
   CASE (    5)
     GO TO 2060
 END SELECT
 
!     LABEL DECLARATIVE - GET LABEL NAME
 
 1270 CALL xscndm
 SELECT CASE ( irturn )
   CASE (    1)
     GO TO 2170
   CASE (    2)
     GO TO 1280
   CASE (    3)
     GO TO 2170
   CASE (    4)
     GO TO 2170
   CASE (    5)
     GO TO 2060
 END SELECT
 
!     CHECK IF LABEL IS FOR CONDITIONAL COMPILATION
 
 1280 CONTINUE
 IF (dmap(dmppnt) /= nskip(ilevel,1) .OR. dmap(dmppnt+1) /=  &
     nskip(ilevel,2)) GO TO 1290
 ilevel = ilevel - 1
 skip   = .false.
 GO TO 300
 1290 IF (skip) GO TO 300
 
!     SCAN LABEL TABLE FOR LABEL NAME
 
 IF (lstlbl < lbltop) GO TO 1310
 DO  j = lbltop,lstlbl,4
   IF (dmap(dmppnt) == lbltbl(j) .AND. dmap(dmppnt+1) == lbltbl(j+1))  &
       GO TO 1340
 END DO
 
!     NAME NOT IN LABEL TABLE, MAKE NEW ENTRY
 
 1310 ASSIGN 1320 TO irturn
 IF (lstlbl+8 >= lstpar) GO TO 2220
 1320 lstlbl = lstlbl + 4
 j = lstlbl
 lbltbl(j  ) = dmap(dmppnt  )
 lbltbl(j+1) = dmap(dmppnt+1)
 lbltbl(j+3) = 0
 1330 lbltbl(j+2) = iseqn + 1
 GO TO 300
 
!     LABEL NAME FOUND IN LABEL TABLE, DEF ENTRY SHOULD BE ZERO
 
 1340 IF (lbltbl(j+2) == 0) THEN
   GO TO  1330
 ELSE
   GO TO  2250
 END IF
 
!     FILE DECLARATIVE
!     SET FILE NAME FLAG
!     DO NOT PROCESS FILE DECLARATION WHEN EXECUTE FLAG IS OFF ON
!     MODIFIED RESTART.
 
 1350 IF (start == imst .AND. oscar(ospnt+5) >= 0) GO TO 1260
 1360 i = 1
 1370 CALL xscndm
 SELECT CASE ( irturn )
   CASE (    1)
     GO TO 1380
   CASE (    2)
     GO TO 1410
   CASE (    3)
     GO TO 2170
   CASE (    4)
     GO TO 300
   CASE (    5)
     GO TO 2060
 END SELECT
 
!     DELIMITER ENCOUNTERED
 
 1380 IF (dmap(dmppnt+1) == islsh) GO TO 1390
 IF (dmap(dmppnt+1) == iequl) GO TO 1400
 GO TO 2170
 
!     DELIMITER IS /, TEST FILE NAME FLAG
 
 1390 IF (i /= 0) GO TO 2170
 GO TO 1360
 
!     DELIMITER IS =, TURN OFF FILE NAME FLAG
 
 1400 i = 0
 GO TO 1370
 
!     NAME ENCOUNTERED - TEST FILE NAME FLAG
 
 1410 IF (i == 0) GO TO 1430
 
!     FILE NAME - ENTER IN FILE TABLE
 
 fpnt = fpnt + 3
 IF (fpnt > lfile-2) GO TO 2410
 FILE(fpnt  ) = dmap(dmppnt  )
 FILE(fpnt+1) = dmap(dmppnt+1)
 
!     PUT FILE NAME INTO LABEL TABLE FOR DMAP XREF
 
 ASSIGN 1420 TO irturn
 IF (lstlbl+8 >= lstpar) GO TO 2220
 1420 lstlbl = lstlbl + 4
 lbltbl(lstlbl  ) = FILE(fpnt  )
 lbltbl(lstlbl+1) = FILE(fpnt+1)
 lbltbl(lstlbl+2) = iseqn
 lbltbl(lstlbl+3) = -1
 GO TO 1370
 
!     FILE PARAMETER FOUND - ENTER APPROPRIATE CODE IN FILE TABLE
 
 1430 DO  j = 1,3
   IF (dmap(dmppnt) == fparam(j)) THEN
      SELECT CASE ( j )
       CASE (    1)
         GO TO 1450
       CASE (    2)
         GO TO 1460
       CASE (    3)
         GO TO 1470
     END SELECT
   END IF
 END DO
 GO TO 2160
 
!     TAPE PARAM
 
 1450 fcode = itape
 GO TO 1480
 
!     APPEND PARAM
 
 1460 fcode = iappnd
 GO TO 1480
 
!     SAVE PARAM
 
 1470 fcode = isave
 
!     PUT CODE IN FILE TABLE
 
 1480 FILE(fpnt+2) = orf(FILE(fpnt+2),fcode)
 GO TO 1370
 
!     PROCESS PRECHK CARD
 
 1490 INDEX = 3
 CALL xscndm
 SELECT CASE ( irturn )
   CASE (    1)
     GO TO 2160
   CASE (    2)
     GO TO 1500
   CASE (    3)
     GO TO 2160
   CASE (    4)
     GO TO 2160
   CASE (    5)
     GO TO 2160
 END SELECT
 
!     TEST FOR  ALL  OPTION OR BLANK
 
 1500 IF (dmap(dmppnt) == nblank) GO TO 1490
 preflg = 1
 nnames = 0
 IF (dmap(dmppnt) == namopt(23)) GO TO 1520
 IF (dmap(dmppnt) ==      nend ) GO TO 1550
 
!     LIST HAS BEEN FOUND, STORE IN /AUTOCM/
 
 1510 nnames = nnames + 1
 IF (nnames > 50) GO TO 2180
 prenam(2*nnames-1) = dmap(dmppnt  )
 prenam(2*nnames  ) = dmap(dmppnt+1)
 CALL xscndm
 SELECT CASE ( irturn )
   CASE (    1)
     GO TO 2160
   CASE (    2)
     GO TO 1510
   CASE (    3)
     GO TO 2160
   CASE (    4)
     GO TO 1560
   CASE (    5)
     GO TO 2060
 END SELECT
 
!     ALL  OPTION FOUND, LOOK FOR  EXCEPT
 
 1520 CALL xscndm
 SELECT CASE ( irturn )
   CASE (    1)
     GO TO 2160
   CASE (    2)
     GO TO 1530
   CASE (    3)
     GO TO 2160
   CASE (    4)
     GO TO 1530
   CASE (    5)
     GO TO 2060
 END SELECT
 1530 IF (dmap(dmppnt) == namopt(25) .AND. dmap(dmppnt+1) == namopt(26))  &
     GO TO 1540
 preflg = 2
 GO TO 1560
 1540 preflg = 3
 CALL xscndm
 SELECT CASE ( irturn )
   CASE (    1)
     GO TO 2160
   CASE (    2)
     GO TO 1510
   CASE (    3)
     GO TO 2160
   CASE (    4)
     GO TO 1560
   CASE (    5)
     GO TO 2060
 END SELECT
 1550 preflg = 0
 1560 IF (icpflg /= 0) GO TO 1240
 preflg = 0
 GO TO 300
 
!     PROCESS XDMAP INSTRUCTION
 
 1570 iold = diag14
 1580 CALL xscndm
 SELECT CASE ( irturn )
   CASE (    1)
     GO TO 2160
   CASE (    2)
     GO TO 1610
   CASE (    3)
     GO TO 2160
   CASE (    4)
     GO TO 1590
   CASE (    5)
     GO TO 2060
 END SELECT
 1590 INDEX = 2
 IF (iold == 0 .OR. ifirst == 0) GO TO 1230
 IF (start /= icst) WRITE (optape,1600) iplus,iplus
 1600 FORMAT (a1,2X,a1)
 GO TO 300
 1610 IF (dmap(dmppnt) == nblank) GO TO 1580
 
!     HAVE LOCATED AN XDMAP OPTION
 
 DO  k = 1,22,2
   IF (dmap(dmppnt) == namopt(k) .AND. dmap(dmppnt+1) == namopt(k+1))  &
       GO TO 1630
 END DO
 GO TO 2190
 1630 kk = k/2 + 1
 SELECT CASE ( kk )
   CASE (    1)
     GO TO 1580
   CASE (    2)
     GO TO 1640
   CASE (    3)
     GO TO 1710
   CASE (    4)
     GO TO 1660
   CASE (    5)
     GO TO 1650
   CASE (    6)
     GO TO 1680
   CASE (    7)
     GO TO 1690
   CASE (    8)
     GO TO 1700
   CASE (    9)
     GO TO 1580
   CASE (   10)
     GO TO 1670
   CASE (   11)
     GO TO 1580
 END SELECT
 1640 iflg(1) = 0
 GO TO 1580
 1650 IF (diag14 == 1) GO TO 1580
 iflg(3) = 0
 diag14  = 0
 GO TO 1580
 1660 IF (diag14 == 1) GO TO 1580
 iflg(3) = 1
 diag14  = 2
 GO TO 1580
 1670 IF (diag4 == 1) GO TO 1580
 iflg(6) = 1
 diag4   = 1
 GO TO 1580
 1680 IF (diag17 == 1) GO TO 1580
 iflg(4) = 1
 diag17  = 2
 GO TO 1580
 1690 IF (diag17 == 1) GO TO 1580
 iflg(4) = 0
 diag17  = 0
 GO TO 1580
 1700 IF (diag25 == 1) GO TO 1580
 iflg(5) = 1
 diag25  = 1
 GO TO 1580
 
!     CODE TO PROCESS  ERR  OPTION
 
 1710 CALL xscndm
 SELECT CASE ( irturn )
   CASE (    1)
     GO TO 1720
   CASE (    2)
     GO TO 2160
   CASE (    3)
     GO TO 2160
   CASE (    4)
     GO TO 2160
   CASE (    5)
     GO TO 2060
 END SELECT
 1720 IF (dmap(dmppnt+1) /= iequl) GO TO 2160
 CALL xscndm
 SELECT CASE ( irturn )
   CASE (    1)
     GO TO 2160
   CASE (    2)
     GO TO 2160
   CASE (    3)
     GO TO 1730
   CASE (    4)
     GO TO 2160
   CASE (    5)
     GO TO 2060
 END SELECT
 1730 iflg(2) = dmap(dmppnt+1)
 IF (iflg(2) < 0 .OR. iflg(2) > 2) GO TO 2190
 GO TO 1580
 
!     PROCESS CONDCOMP INSTRUCTION
 
 1740 IF (ilevel >= 5) GO TO 2160
 ion = 0
 IF (dmap(dmppnt+1) == cdcomp(2)) ion = 1
 CALL xscndm
 SELECT CASE ( irturn )
   CASE (    1)
     GO TO 2160
   CASE (    2)
     GO TO 1750
   CASE (    3)
     GO TO 1760
   CASE (    4)
     GO TO 2160
   CASE (    5)
     GO TO 2060
 END SELECT
 
!     LABEL SPECIFIED FOR END
 
 1750 nskip(ilevel+1,1) = dmap(dmppnt  )
 nskip(ilevel+1,2) = dmap(dmppnt+1)
 GO TO 1770
 
!     INSTRUCTION COUNT GIVEN FOR END
 
 1760 CONTINUE
 IF (dmap(dmppnt+1) < 0) GO TO 2160
 nskip(ilevel+1,1) = dmap(dmppnt+1)
 
!     GET LABEL AND LOOK FOR IT IN PVT
 
 1770 CALL xscndm
 SELECT CASE ( irturn )
   CASE (    1)
     GO TO 2160
   CASE (    2)
     GO TO 1780
   CASE (    3)
     GO TO 2160
   CASE (    4)
     GO TO 2160
   CASE (    5)
     GO TO 2060
 END SELECT
 1780 ilevel = ilevel + 1
 kdh = 3
 1790 length = andf(pvt(kdh+2),nosgn)
 length = itype(length)
 IF (dmap(dmppnt) == pvt(kdh) .AND. dmap(dmppnt+1) == pvt(kdh+1)) GO TO 1810
 kdh = kdh + length + 3
 IF (kdh - pvt(2) < 0) THEN
   GO TO  1790
 END IF
 
!     PARAMETER NOT FOUND - ASSUME FALSE VALUE
 
 1800 IF (ion == 0) GO TO 300
 GO TO 1820
 
!     CHECK IF VALUE IS FALSE
 
 1810 pvt(kdh+2) = orf(pvt(kdh+2),isgnon)
 IF (andf(pvt(kdh+2),nosgn) /=  1) GO TO 2160
 IF (pvt(kdh+3) < 0 .AND. ion == 1) GO TO 300
 IF (pvt(kdh+3) >= 0 .AND. ion == 0) GO TO 300
 1820 skip = .true.
 GO TO 300
 
!     ***********************************************************
!     DMAP INSTRUCTIONS ALL PROCESSED - PREPARE OSCAR FOR PHASE 2
!     ***********************************************************
 
!     CHECK FOR DISCREPENCY BETWEEN RIGID FORMAT AND MED TABLE.
 
 1900 IF (med(medtp) /= dmpcnt .AND. iapp /= idmapp) GO TO 2390
 
!     USE LBLTBL PARAMETER NAMES TO UPDATE VALUE SECTIONS OF TYPE C AND
!     E OSCAR ENTRIES.
 
 1910 IF (lstpar >= lblbot) GO TO 1990
 
!     FIND PARAMETER NAME IN VPS
 
 k = 3
 1920 IF (lbltbl(lstpar) == vps(k) .AND. lbltbl(lstpar+1) == vps(k+1))  &
     GO TO 1930
 k = k + andf(vps(k+2),maskhi) + 3
 IF (k - vps(2) < 0) THEN
   GO TO  1920
 ELSE
   GO TO  1950
 END IF
 
!     NAME FOUND IN VPS, VPS POINTER TO OSCAR VALUE SECTION.
 
 1930 i = lbltbl(lstpar+2)
 oscar(i) = k + 3
 
!     GET NEXT ENTRY FROM LBLTBL
 
 1940 lstpar = lstpar + 4
 GO TO 1910
 
!     SEARCH PVT TABLE FOR PARAMETER. IF FOUND ENTER PARAMETER IN VPS.
 
 1950 k1 = 3
 1960 length = andf(pvt(k1+2),nosgn)
 length = itype(length)
 IF (lbltbl(lstpar) == pvt(k1) .AND. lbltbl(lstpar+1) == pvt(k1+1))  &
     GO TO 1970
 k1 = k1 + length + 3
 IF (k1-pvt(2) < 0) THEN
   GO TO  1960
 ELSE
   GO TO  2310
 END IF
 1970 k = vps(2) + 1
 pvt(k1+2) = orf(pvt(k1+2),isgnon)
 vps(2) = k + 2 + length
 IF (vps(2) >= vps(1)) GO TO 2380
 k2 = length + 3
 DO  m = 1,k2
   j  = k  + m - 1
   j1 = k1 + m - 1
   vps(j) = pvt(j1)
 END DO
 GO TO 1930
 
!     USE LBLTBL ENTRIES TO LOAD SEQUENCE NOS. INTO VALUE SECTION OF
!     TYPE C OSCAR ENTRIES.
 
 1990 lblerr = 0
 lstlsv = lstlbl
 2000 IF (lstlbl < lbltop) GO TO 2050
 IF (lbltbl(lstlbl+2) == 0) GO TO 2330
 
!     IGNORE FILE NAMES IN LBLTBL USED FOR XREF
 
 2010 IF (lbltbl(lstlbl+3) < 0) THEN
   GO TO  2040
 ELSE IF (lbltbl(lstlbl+3) == 0) THEN
   GO TO  2360
 END IF
 2020 i = lbltbl(lstlbl+3) + 6
 IF (oscar(i-3) == ncond .OR. oscar(i-3) == njump) GO TO 2030
 j = oscar(i)
 
!     LABEL NAME TO WORDS 3 AND 4 OF CEITBL ENTRY
 
 ceitbl(j+1) = lbltbl(lstlbl  )
 ceitbl(j+2) = lbltbl(lstlbl+1)
 
!     OSCAR RECORD NO. OF BEGIN LOOP TO FIRST WORD OF CEITBL ENTRY
 
 ceitbl(j-1) = orf(lshift(lbltbl(lstlbl+2),16),ceitbl(j-1))
 2030 oscar(i)    = orf(lshift(lbltbl(lstlbl+2),16),oscar(i))
 
!     GET NEXT LBLTBL ENTRY.
 
 2040 lstlbl = lstlbl - 4
 GO TO 2000
 
!     NORMAL RETURN -     DUMP LBLTBL ONTO SCRATCH FOR DMAP XREF
!                         THEN DELETE LBLTBL AND DMPCRD ARRARYS
!                         FROM OPEN CORE
 
 2050 lstlbl = lstlsv
 2060 loscar = lblbot
 idpbuf = korsz(oscar) - 2*bufsz
 CALL CLOSE (nscr,1)
 lstlbl = lstlbl - lbltop + 4
 IF (lstlbl < 0) lstlbl = 0
 RETURN
 
!     DIAGNOSTIC MESSAGES -
 
!     DMAP INPUT FILE ERROR
 
 2100 CALL xgpidg (-10,ospnt,0,0)
 GO TO 410
 
!     DMAP OUTPUT FILE ERROR
 
 2110 CALL xgpidg (-11,ospnt,0,0)
 GO TO 420
 
!     NO MACRO INSTRUCTION NAME ON DMAP CARD.
 
 2120 CALL xgpidg (12,0,dmpcnt,0)
 GO TO 300
 
!     NO MPL ENTRY FOR THIS DMAP MACRO INSTRUCTION
 
 2130 CALL xgpidg (13,0,dmppnt,dmpcnt)
 GO TO 300
 
!     MPL TABLE INCORRECT
 
 2140 CALL xgpidg (49,0,0,0)
 GO TO 2500
 
!     DUPLICATE PARAMETER NAMES (WARNING)
 
 2150 CALL xgpidg (-2,ospnt,dmap(dmppnt),dmap(dmppnt+1))
 GO TO 870
 
!     DMAP FORMAT ERROR
 
 2160 CALL xgpidg (16,ospnt,0,0)
 GO TO 300
 2170 j = osbot + oscar(osbot) + 6
 CALL xgpidg (16,j,0,0)
 GO TO 300
 
!     PRECHK NAME LIST OVERFLOW
 
 2180 CALL xgpidg (55,0,0,0)
 GO TO 2500
 
!     ILLEGAL OPTION ON XDMAP CARD
 
 2190 CALL xgpidg (56,0,0,0)
 GO TO 300
 
!     VARIABLE REPT INSTRUCTION ERRORS
 
 2200 CALL xgpidg (58,0,0,0)
 GO TO 300
 2210 CALL xgpidg (57,0,0,0)
 GO TO 300
 
!     LBLTBL OVERFLOWED - ALLOCATE 50 MORE WORDS FOR IT.
 
 2220 icrdtp = icrdtp - 50
 IF (icrdtp < oscar(osbot)+osbot) GO TO 2240
 loscar = loscar - 50
 
!     MOVE LABEL NAME PORTION OF LBLTBL
 
 jx = lstlbl + 3
 DO  ix = lbltop,jx
   iy = ix - 50
   lbltbl(iy) = lbltbl(ix)
 END DO
 lbltop = lbltop - 50
 lstlbl = lstlbl - 50
 GO TO irturn, (570,620,1070,1320,1420)
 
!     LABEL TABLE OVERFLOW, DISCONTINUE COMPILATION
 
 2240 CALL xgpidg (14,nlblt1,nlblt2,dmpcnt)
 GO TO 2500
 
!     LABEL IS MULTIPLY DEFINED
 
 2250 CALL xgpidg (19,dmpcnt,dmppnt,0)
 GO TO 300
 
!     ILLEGAL CHARACTERS IN DMAP SAVE PARAMETER NAME LIST
 
 2260 CALL xgpidg (20,ospnt,oscar(i)+1,0)
 GO TO 870
 
!     XSAVE PARAMETER NAME NOT ON PRECEDING DMAP CARD
 
 2270 CALL xgpidg (21,ospnt,dmap(dmppnt),dmap(dmppnt+1))
 GO TO 870
 
!     CEITBL OVERFLOW, DISCONTINUE COMPILATION
 
 2280 CALL xgpidg (14,nceit1,nceit2,dmpcnt)
 GO TO 2500
 
!     CHECK FOR XSAVE PARAMETERS NOT ON PRECEDING DMAP CARD
 
 2290 i1 = i + 2
 k  = k + 1
 DO  k1 = i1,k,2
   IF (oscar(k1) > 0 .OR. oscar(k1-1) == 0) CYCLE
   j = oscar(k1-1)
   CALL xgpidg (21,ospnt,vps(j-3),vps(j-2))
 END DO
 GO TO 300
 
!     PARAMETER NOT DEFINED FOR USE IN COND, PURGE OR EQUIV INSTRUCTIONS
 
 2310 CALL xgpidg (25,lbltbl(lstpar+3),lbltbl(lstpar),lbltbl(lstpar+1))
 GO TO 1940
 
!     LABEL NOT DEFINED
 
 2320 CALL xgpidg (26,lbltbl(lstlbl+3),lbltbl(lstlbl),lbltbl(lstlbl+1))
 nogo = 1
 GO TO 2040
 
!     CHECK FOR LABEL DEFINED
 
 2330 DO  j = lbltop,lstlbl,4
   IF (lbltbl(j) == lbltbl(lstlbl) .AND. lbltbl(j+1) ==  &
       lbltbl(lstlbl+1) .AND. lbltbl(j+2) > 0) GO TO 2350
 END DO
 GO TO 2320
 2350 lbltbl(lstlbl+2) = lbltbl(j+2)
 GO TO 2010
 
!     LABEL NOT REFERENCED - WARNING ONLY
 
 2360 CALL xgpidg (-27,lbltbl(lstlbl+2),lbltbl(lstlbl),lbltbl(lstlbl+1))
 GO TO 2040
 
!     TIME SEGMENT NAME INCORRECT - WARNING ONLY
 
 2370 CALL xgpidg (-17,ospnt,0,0)
 GO TO 300
 
!     VPS TABLE OVERFLOWED
 
 2380 CALL xgpidg (14,nvps,nblank,0)
 GO TO 2500
 
!     DMAP SEQUENCE DOES NOT CORRESPOND TO MED TABLE
 
 2390 CALL xgpidg (39,0,0,0)
 GO TO 2500
 
!     WARNING - CANNOT CHECKPOINT USER INPUT
 
 2400 CALL xgpidg (-48,ospnt,oscar(iospnt),oscar(iospnt+1))
 GO TO 1030
 
!     OVERFLOWED FILE TABLE
 
 2410 CALL xgpidg (14,nfile,nblank,0)
 GO TO 2500
 
!     SAVE OUT OF POSITION
 
 2420 CALL xgpidg (61,ospnt,0,0)
 ospnt = iosdav
 osprc = os2b4
 GO TO 300
 
!     RETURN WHEN XGPI HAS BEEN DISASTERED.
 
 2500 nogo = 2
 RETURN
END SUBROUTINE xosgen
