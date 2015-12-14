SUBROUTINE xgpi
     
!     THE PURPOSE OF XGPI IS TO INITIALIZE AND CALL THE FOLLOWING
!     SUBROUTINES - XOSGEN AND XFLORD.
 
 IMPLICIT INTEGER (a-z)
 EXTERNAL        lshift,rshift,andf,orf,complf
 LOGICAL :: lnogo
 INTEGER :: cnm(1),fnm(1),ptdic(1)
 DIMENSION       icpdpl(1),isol(1),icf(1),iccnam(1),ibufr(1),  &
     ibf(8),ioshdr(2),itype(6),itrl(7),dmpcrd(1),  &
     med(1),nxgpi(2),nxptdc(2),oscar(1),os(5)
 CHARACTER (LEN=25) :: sfm
 CHARACTER (LEN=29) :: uim
 CHARACTER (LEN=25) :: uwm
 CHARACTER (LEN=23) :: ufm
 COMMON /xmssg / ufm,uwm,uim,sfm
 COMMON /xoldpt/ ptdtop,ptdbot,lptdic,nrlfl,seqno
 COMMON /stapid/ tapid(6),otapid(6)
 COMMON /xgpi5 / isol,start,alter(2),sol,subset,iflag,iestim,  &
     icftop,icfpnt,lctlfl,ictlfl(1)
 COMMON /moddmp/ iflg(6)
 COMMON /xgpi6 / medtp,fnmtp,cnmtp,medpnt,lmed,iplus,diag14,  &
     diag17,diag4,diag25,ifirst,ibuff(20)
 COMMON /xgpi8 / icptop,icpbot,lcpdpl
 COMMON /ifpx0 / lbd,lcc,mjmsk(1)
 COMMON /ifpx1 / ncds,mjcd(1)
 COMMON /two   / two(32)
 COMMON /xmdmsk/ nmskcd,nmskfl,nmskrf,medmsk(7)
 COMMON /system/ ibufsz,optape,nogo,sys4,mpc,spc,sys7,load,sys9(2),  &
     pagect,sys12(7),iecho,sys20,apprch,sys22(2),  &
     icfiat,sys25,cppgct,sys27(42),sscell,sys70(7), bandit,sys78(4),icpflg
 COMMON /l15 l8/ l15,l8
 COMMON /xgpi4 / irturn,insert,iseqn,dmpcnt,  &
     idmpnt,dmppnt,bcdcnt,length,icrdtp,ICHAR,newcrd, modidx,licf,isavdw,dmap(1)
 
!                  ** CONTROL CARD NAMES **
!                  ** DMAP    CARD NAMES **
 COMMON /xgpic / icold,islsh,iequl,nblank,nxequi,  &
     nmed,nsol,ndmap,nestm1,nestm2,nexit,  &
     nbegin,nend,njump,ncond,nrept,ntime,nsave,noutpt,  &
     nchkpt,npurge,nequiv,ncpw,nbpc,nwpc,  &
     maskhi,masklo,isgnon,nosgn,iallon,masks(1)
 COMMON /zzzzzz/ core(1)
 COMMON /xgpi2 / libf,mplpnt,mpl(1)
!WKBR COMMON /XGPI3 / PVT(6)
 COMMON /xgpi3 / pvt(200)
 COMMON /xdpl  / dpl(3)
 COMMON /xvps  / vps(4)
 COMMON /xfist / ifist(1)
 COMMON /xfiat / ifiat(3)
!WKBR COMMON /XCEITB/ CEITBL(2)
 COMMON /xceitb/ ceitbl(42)
 COMMON /xgpid / icst,iunst,imst,ihapp,idsapp,idmapp,  &
     isave,itape,modflg,intgr,losgn, noflgs,seteor,eotflg,ieqflg,  &
     cpntry(7),jmp(7)
 COMMON /xgpie / nscr
 EQUIVALENCE     (loscar,os(1),core(1)), (osprc,os(2)),  &
     (osbot ,os(3)), (ospnt,os(4)), (oscar(1),os(5)),  &
     (oscar(1),med(1),fnm(1),cnm(1),icpdpl(1)),  &
     (oscar(1),ibufr(1),dmpcrd(1),ptdic(1))
 EQUIVALENCE     (dmap(1),icf(1)), (nmed,iccnam(1)),  &
     (mpl(1),ibf(1)), (dpl(1),ndpfil),  &
     (dpl(2),maxdpl), (dpl(3),lstdpl), (nogo,lnogo)
 
!           ** DEFINITION OF PROGRAM VARIABLES **
!     LICF   = NUMBER OF WORDS IN ICF ARRAY
!     WRNGRL = COUNTER FOR NUMBER OF TIMES WRONG REEL WAS MOUNTED.
!     FILCON = FLAG INDICATING FILE IS CONTINUED ON NEXT REEL.
!     IDPFCT = DATA POOL FILE NUMBER OF OSCAR FILE
!     IOSHDR = ARRAY CONTAINING HEADER RECORD FOR XOSCAR FILE IN IDP.
!     PTFCT  = PROBLEM TAPE FILE POSITION
!     EORFLG = END OF RECORD FLAG
!     ICST   = COLD START FLAG
!     IUNST  = UNMODIFIED RESTART
!     IMST   = MODIFIED   RESTART
 
!           ** VARIABLES USED IN GINO CALLS **
!     NPTBUF = NEW PROBLEM TAPE BUFFER AREA
!     IOPBUF = OLD PROBLEM TAPE BUFFER AREA
!     IDPBUF = DATA POOL FILE BUFFER AREA
!     NPTWRD = NUMBER OF WORDS READ FROM NEW PROBLEM TAPE
!     IOPWRD = NUMBER OF WORDS READ FROM OLD PROBLEM TAPE
!     IDPWRD = NUMBER OF WORDS READ FROM DATA POOL FILE
 
!           ** SYMBOLS EQUATED TO CONSTANTS **
!     NPT    = NEW PROBLEM TAPE GINO I.D. NAME (NPTP)
!     IOP    = OLD PROBLEM TAPE GINO I.D. NAME (OPTP)
!     IDP    = DATA POOL FILE GINO I.D. NAME   (POOL)
!     NSCR   = SCRATCH FILE USED FOR RIGID FORMAT. DATA IN NSCR WAS
!              PASSED OVER BY XCSA. IT MUST BE THE LAST SCRATCH FILE
!              IN LINK1, AND NOT TO BE OVER WRITTEN BY XSORT2
!              (CURRENTLY, NSCR = 315)
 
 DATA    jcard ,jfile  /  4HCARD,4HFILE /, idp   /4HPOOL /
 DATA    npvt  /4HPVT  /, ixtim /4HXTIM /, npt   /4HNPTP /,  &
     iop   /4HOPTP /, idpwrd/0      /, iopwrd/0      /,  &
     ndpl  /4HDPL  /, ioshdr/4HXOSC  , 4HAR          /,  &
     nptwrd/0      /, filcon/0      /, nxvps /4HXVPS /
 DATA    itype /1,1,2,2 , 2,4           /
 DATA    nxcsa /4HXCSA /, nxaltr/4HXALT /, nparam/216    /
 DATA    nxptdc/4HXPTD  , 4HIC  /,nxgpi /  4HXGPI,4H     /
 
 
!     LOAD COMMON AREAS AND PERFORM INITIAL CALCULATIONS
 
 CALL xgpidd
 CALL xmpldd
 CALL xlnkdd
 
!     INITIALIZE
 
 nscr  = 315
 CALL sswtch ( 4,diag4 )
 CALL sswtch (14,diag14)
 CALL sswtch (17,diag17)
 CALL sswtch (25,diag25)
 IF (diag14 == 1) iflg(3) = 1
 IF (diag17 == 1) iflg(4) = 1
 IF (diag4  == 1) iflg(6) = 1
 IF (diag4  /= 0) iflg(5) = orf(iflg(5),lshift(1,16))
 IF (diag25 == 1) iflg(5) = 1
 
!     SET DMAP COMPILER DEFAULT OPTION TO LIST FOR
!     APPROACH DMAP RUNS, RESTART RUNS AND SUBSTRUCTURE RUNS
!     RESET TO NO LIST IF ECHO=NONO (IECHO=-2)
 
 IF (iecho /= -2 .AND. (apprch < 2 .OR. sscell /= 0)) iflg(3) = 1
 IF (diag14 == 0 .AND. iflg(3) == 1) diag14 = 2
 IF (iecho == -2) diag14 = 0
 CALL xgpimw (1,0,0,0)
 
 CALL xgpibs
 IF (nogo > 1) GO TO 2210
 
!     SET UP GINO BUFFER AREAS FOR OLD PROBLEM TAPE,NEW PROBLEM TAPE
!     AND DATA POOL TAPE.
 
 loscar = korsz(ibufr)
 nptbuf = loscar - ibufsz
 
!     OLD PROBLEM TAPE AND NEW PROBLEM TAPE SHARE BUFFER
 
 iopbuf = nptbuf
 idpbuf = nptbuf - ibufsz
 loscar = idpbuf - 1
 
!     ALLOW MINIMAL SIZE FOR MED ARRAY RESIDING IN OPEN CORE.
!     WE WILL EXPAND MED IF NECESSARY.
 
 medtp = loscar
 lmed  = 1
 IF (loscar < 1) GO TO 2080
 
!     OPEN NEW PROBLEM TAPE AS INPUT FILE
 
 CALL OPEN (*1900,npt,ibufr(nptbuf),0)
 
!     NUMBER OF FILE ON NPT + 1
 
 nrlfl = lshift(tapid(6),16) + 5
 
!     FILE POSITION OF IOP AT ENTRY TO XGPI
 
 ptfct = lshift(otapid(6),16) + 4
 
!     FIND XCSA FILE ON NEW PROBLEM TAPE
 
 nam1 = nxcsa
 nam2 = nblank
 
!     SKIP HEADER FILE
 
 10 CALL skpfil (npt,1)
 CALL READ (*1950,*1950,npt,icf,1,1,nptwrd)
 
!     CHECK FOR ALTER FILE
!     SET DIAG14 TO 10 IF ALTER CARDS ARE PRESENT. DIAG14 WOULD BE
!     CHANGED TO 11 IF DMAP CONTAINS POTENTIAL FATAL ERROR. IN SUCH
!     CASE, DMAP LISTING WILL BE PRINTED.
 
 IF (icf(1) /= nxaltr) GO TO 15
 nrlfl = nrlfl + 1
 IF (diag14 == 0)                 diag14 = 10
 
!     CHECK FOR CHECKPOINT DICTIONARY FILE
 
 15 IF (icf(1) == nxptdc(1)) nrlfl = nrlfl + 1
 
!     CHECK FOR CONTROL FILE
 
 IF (icf(1) /= nxcsa) GO TO 10
 
!     PROBLEM TAPE IS POSITIONED AT EXECUTIVE CONTROL FILE.
 
 icfpnt = icftop
 
!     READ THE SIX-WORD DATA RECORD
 
 CALL READ (*1950,*20,npt,isol,7,1,nptwrd)
 GO TO 1950
 20 CALL CLOSE (npt,1)
 IF (IABS(apprch) == 1) GO TO 620
 
!     FILL MED ARRAY
 
 medtp = 1
 
!     SET VALUE FOR NUMBER OF WORDS PER MED ENTRY
 
 med(medtp+1) = 1
 IF (start /= icst) med(medtp+1) = nmskcd + nmskfl + nmskrf
 
 CALL gopen (nscr,ibufr(nptbuf),0)
 lloscr = loscar - 2
 
!     READ THE MED TABLE
 
 CALL READ (*1960,*30,nscr,med(medtp+2),lloscr,1,lmed)
 GO TO 2090
 
!     SET VALUE FOR NUMBER OF DMAP INSTRUCTIONS
 
 30 med(medtp) = lmed/med(medtp+1)
 
!     CHECK FOR ILLEGAL NUMBER OF WORDS IN MED TABLE RECORD
 
 IF (start /= icst .AND. lmed /= med(medtp)*med(medtp+1)) GO TO 1980
 
!     SET THE POINTERS TO THE FILE NAME AND CARD NAME TABLES
 
 fnmtp = medtp + lmed + 2
 cnmtp = fnmtp
 IF (start == icst) GO TO 600
 lloscr = lloscr - lmed
 
!     READ THE FILE NAME TABLE
 
 CALL skprec (nscr,1)
 jtype = jfile
 CALL READ (*1970,*40,nscr,med(fnmtp+1),lloscr,1,lmed)
 GO TO 2090
 
!     SET THE VALUE FOR THE NUMBER OF ENTRIES IN THE FILE NAME TABLE
 
 40 med(fnmtp) = lmed/3
 
!     CHECK FOR ILLEGAL NUMBER OF WORDS IN FILE NAME TABLE RECORD
 
 IF (lmed /= 3*med(fnmtp)) GO TO 1990
 
!     CHECK FOR ILLEGAL BIT NUMBERS IN FILE NAME TABLE
 
 istrbt = 31*nmskcd + 1
 iendbt = 31*(nmskcd+nmskfl)
 DO  j = 3,lmed,3
   IF (med(fnmtp+j) < istrbt .OR. med(fnmtp+j) > iendbt) GO TO 2000
 END DO
 
!     RESET THE POINTER FOR THE CARD NAME TABLE
 
 cnmtp  = fnmtp + 3*fnm(fnmtp) + 1
 lloscr = lloscr - lmed
 
!     READ THE CARD NAME TABLE
 
 CALL skprec (nscr,-2)
 jtype = jcard
 CALL READ (*1970,*60,nscr,med(cnmtp+1),lloscr,1,lmed)
 GO TO 2090
 
!     SET THE VALUE FOR THE NUMBER OF ENTRIES IN THE CARD NAME TABLE
 
 60 med(cnmtp) = lmed/3
 
!     CHECK FOR ILLEGAL NUMBER OF WORDS IN CARD NAME TABLE RECORD
 
 IF (lmed /= 3*med(cnmtp)) GO TO 1990
 
!     CHECK FOR ILLEGAL BIT NUMBERS IN CARD NAME TABLE
 
 istrbt = 1
 iendbt = 31*nmskcd
 DO  j = 3,lmed,3
   IF (med(cnmtp+j) < istrbt .OR. med(cnmtp+j) > iendbt) GO TO 2000
 END DO
 
!     RESTART - CHECK MEDMSK TABLE
!     IF MEDMSK WORD(S), CORRESPONDING TO RIGID FORMAT SWITCH, IS(ARE)
!     NON-ZERO, SOLUTION HAS BEEN CHANGED.
!     RESET ENTRY SEQUENCE NO. TO INFINITE IF SOLUTION IS CHANGED.
 
 nmask = med(medtp+1)
 ibegn = nmskcd + nmskfl + 1
 DO  i = ibegn,nmask
   IF (medmsk(i) == 0) CYCLE
   seqno = masklo
   start = imst
 END DO
 
!     SEE IF ANY BULK DATA OR CASE CONTROL CARDS HAVE BEEN MODIFIED.
 
 bgnmsk = 1
 endmsk = lbd + lcc
 
!     TURN OFF BIT IN MJMSK ARRAY IF THE CORRESPONDING CARD NAME
!     IS NOT IN THE CARD NAME RESTART TABLE
 
 i1 = cnmtp + 1
 i2 = i1 + 3*cnm(cnmtp) - 3
 DO  lx = bgnmsk,endmsk
   IF (mjmsk(lx) == 0) CYCLE
   l = lx - bgnmsk + 1
   loop100:  DO  l1 = 2,32
     IF (andf(mjmsk(lx),two(l1)) == 0) CYCLE loop100
     
!     IGNORE BIT IF IT CORRESPONDS TO QOUT$ OR BOUT$
     
     IF (lx == lbd+2 .AND. (l1 == 3 .OR. l1 == 4)) CYCLE loop100
     i = 62*(l-1) + 2*(l1-2) + 1
     DO  ii = i1,i2,3
       IF (mjcd(i) == cnm(ii) .AND. mjcd(i+1) == cnm(ii+1)) CYCLE loop100
     END DO
     ii = complf(two(l1))
     mjmsk(lx) = andf(mjmsk(lx),ii)
   END DO loop100
 END DO
 IF (start == imst) GO TO 130
 
!     DETERMINE TYPE OF RESTART
 
 INDEX = 0
 iend  = lbd
 DO  l = bgnmsk,iend
   IF (mjmsk(l) == 0) CYCLE
   INDEX = 1
   EXIT
 END DO
 130 l = lbd + 1
 IF (start == imst) GO TO 160
 IF (INDEX ==    1) GO TO 150
 IF (mjmsk(l) == 0) GO TO 140
 
!     CHECK FOR NOLOOP$ AND LOOP$
!                                          2**21
 IF (mjmsk(l) /= 1 .AND. mjmsk(l) /= two(11)) GO TO 150
 
!     CHECK FOR GUST$
!                         2**30
 140 IF (mjmsk(l+1) < two(2)) GO TO 170
 150 start = imst
 
!     TURN ON POUT$ IF QOUT$ IS ON
!                         2**29                                  2**14
 160 IF (andf(mjmsk(l+1),two(3)) /= 0) mjmsk(l)=orf(mjmsk(l),two(18))
 
!     TURN ON AOUT$ IF BOUT$ IS ON
!                         2**28                                  2**22
 IF (andf(mjmsk(l+1),two(4)) /= 0) mjmsk(l)=orf(mjmsk(l),two(10))
 
!     TURN OFF BOUT$ AND QOUT$
!                 2**28    2**29
 ii = complf(two(4) + two(3))
 mjmsk(l+1) = andf(mjmsk(l+1),ii)
 
!     TURN OFF NOLOOP$ FOR UNMODIFIED RESTARTS
 
 170 IF (start == iunst .AND. mjmsk(lbd+1) == 1) mjmsk(lbd+1) = 0
 180 CALL page1
 IF (start /= iunst) GO TO 200
 WRITE  (optape,190) uim
 190 FORMAT (a29,' 4143, THIS IS AN UNMODIFIED RESTART.')
 bandit = -1
 IF (apprch == -1) GO TO 700
 GO TO 600
 200 CALL page2 (-2)
 IF (seqno /= masklo) WRITE (optape,210) uim
 IF (seqno == masklo) WRITE (optape,220) uim
 210 FORMAT (a29,' 4144, THIS IS A MODIFIED RESTART.')
 220 FORMAT (a29,' 4145, THIS IS A MODIFIED RESTART INVOLVING RIGID ',  &
     'FORMAT SWITCH.')
 ibulk = 0
 icase = 0
 DO  l = 1,lbd
   IF (mjmsk(l) == 0) CYCLE
   ibulk = 1
   EXIT
 END DO
 240 lbd1 = lbd + 1
 lbdlcc = lbd + lcc
 DO  l = lbd1,lbdlcc
   IF (mjmsk(l) == 0) CYCLE
   icase = 1
   EXIT
 END DO
 260 IF (ibulk /= 0 .OR. icase /= 0) GO TO 290
 IF (seqno == masklo) GO TO 270
 WRITE (optape,460)
 CALL mesage (-61,0,0)
 270 WRITE  (optape,280) uim
 280 FORMAT (a29,'. THERE ARE NO CASE CONTROL OR BULK DATA DECK ',  &
     'CHANGES AFFECTING THIS RESTART.')
 GO TO 600
 290 CALL page2 (-4)
 WRITE  (optape,300) uim
 300 FORMAT (a29,'. CASE CONTROL AND BULK DATA DECK CHANGES AFFECTING',  &
     ' THIS RESTART ARE INDICATED BELOW.',/)
 DO  llx = 1,2
   IF (llx == 1) GO TO 360
   CALL page2 (-3)
   WRITE  (optape,310) uim
   310 FORMAT (a29,'. EFFECTIVE BULK DATA DECK CHANGES', /1X,32(1H-))
   IF (ibulk /= 0) GO TO 330
   CALL page2 (-3)
   WRITE  (optape,320)
   320 FORMAT (//,' NONE',/)
   CYCLE
   330 CALL page2 (-3)
   IF (apprch /= -1) WRITE (optape,340)
   IF (apprch == -1) WRITE (optape,350)
   340 FORMAT (//,' MASK WORD - BIT POSITION - CARD/PARAM NAME - PACKED',  &
       ' BIT POSITION',/)
   350 FORMAT (//,' MASK WORD - BIT POSITION - CARD/PARAM NAME',/)
   lim1 = 1
   lim2 = lbd
   GO TO 410
   360 CALL page2 (-3)
   WRITE  (optape,370) uim
   370 FORMAT (a29,'. EFFECTIVE CASE CONTROL DECK CHANGES', /1X,35(1H-))
   IF (icase /= 0) GO TO 380
   CALL page2 (-3)
   WRITE (optape,320)
   CYCLE
   380 CALL page2 (-3)
   IF (apprch /= -1) WRITE (optape,390)
   IF (apprch == -1) WRITE (optape,400)
   390 FORMAT (//,' MASK WORD - BIT POSITION ---- FLAG NAME ---- PACKED',  &
       ' BIT POSITION',/)
   400 FORMAT (//,' MASK WORD - BIT POSITION ---- FLAG NAME',/)
   lim1 = lbd1
   lim2 = lbdlcc
   410 DO  l = lim1,lim2
     IF (mjmsk(l) == 0) CYCLE
     CALL page2 (-1)
     WRITE  (optape,420) l
     420 FORMAT (1X,i5)
     loop480:  DO  k = 2,32
       IF (andf(mjmsk(l),two(k)) == 0) CYCLE loop480
       
!     GET CORRESPONDING CARD NAME FROM MAIN CARD TABLE
       
       i  = 62*(l-1) + 2*(k-2) + 1
       kz = k - 1
       CALL page2 (-1)
       IF (apprch /= -1) GO TO 430
       WRITE (optape,440) kz,mjcd(i),mjcd(i+1)
       CYCLE loop480
       
!     SEARCH RIGID FORMAT CARD NAME RESTART TABLE FOR A MATCH
       
       430 DO  ii = i1,i2,3
         IF (mjcd(i) /= cnm(ii) .OR. mjcd(i+1) /= cnm(ii+1)) CYCLE
         
!     CARD NAME FOUND - SET BIT IN MEDMSK
         
         WRITE  (optape,440) kz,mjcd(i),mjcd(i+1),cnm(ii+2)
         440 FORMAT (17X,i3,11X,2A4,14X,i3)
         l1 = (cnm(ii+2)-1)/31
         ll = l1 + 1
         kk = cnm(ii+2) - 31*l1 + 1
         medmsk(ll) = orf(medmsk(ll),two(kk))
         CYCLE loop480
       END DO
       WRITE  (optape,460) sfm
       460 FORMAT (a25,' 4146, LOGIC ERROR IN SUBROUTINE XGPI WHILE ',  &
           'PROCESSING DATA CHANGES FOR MODIFIED RESTART.')
       WRITE  (optape,470) mjcd(i),mjcd(i+1),(cnm(ll),cnm(ll+1), ll=i1,i2,3)
       470 FORMAT (/10X,2A4, //,10(4X,2A4))
       CALL mesage (-61,0,0)
     END DO loop480
   END DO
 END DO
 IF (apprch == -1) GO TO 700
 
!     MOVE MED AND FILE NAME TABLES TO BOTTOM OF OPEN CORE.
 
 600 CALL CLOSE (nscr,1)
 lmed = cnmtp - medtp
 DO  i = 1,lmed
   ll = medtp + lmed - i
   m  = loscar - i + 1
   med(m) = med(ll)
 END DO
 medtp  = loscar - lmed + 1
 fnmtp  = medtp + med(medtp)*med(medtp+1) + 2
 loscar = medtp - 1
 
!     DETERMINE TYPE OF RESTART IF IT IS A RESTART OF A DMAP RUN
 
 620 IF (apprch /=  -1) GO TO 700
 IF (mjmsk(lbd+1) == 0) GO TO 630
 
!     CHECK FOR NOLOOP$ AND LOOP$
!                                                  2**21
 IF (mjmsk(lbd+1) /= 1 .AND. mjmsk(lbd+1) /= two(11)) GO TO 650
 mjmsk(lbd+1) = 0
 
!     CHECK FOR GUST$
!                           2**30
 630 IF (mjmsk(lbd+2) >= two(2)) GO TO 650
 DO  l = 1,lbd
   IF (mjmsk(l) /= 0) GO TO 650
 END DO
 GO TO 180
 650 start = imst
 seqno = lshift(1,16)
 GO TO 180
 
!     CONTROL FILE LOADED, LOAD PVT TABLE
!     BUMP NUMBER OF FILES IF OLD PROBLEM TAPE HAD ALTERS
 
 700 ptfct = ptfct + alter(2)
 itrl(1) = nparam
 CALL rdtrl (itrl(1))
 IF (itrl(2) <= 0) GO TO 760
 CALL OPEN (*1900,nparam,ibufr(nptbuf),0)
 CALL READ (*760,*710,nparam,pvt(6),2,1,nptwrd)
 710 IF (pvt(6) /= npvt) GO TO 1950
 i = 3
 
!      LOAD PVT VALUES INTO PVT TABLES
 
 720 CALL READ (*740,*730,nparam,pvt(i),pvt(1)-i+1,0,nptwrd)
 GO TO 2020
 730 i = i + nptwrd
 GO TO 720
 740 pvt(2) = i - 1
 CALL CLOSE (nparam,1)
 
!     ELIMINATE TRAILER SO FILE WILL BE DELETED
 
 DO  i = 2,7
   itrl(i) = 0
 END DO
 CALL wrttrl (itrl(1))
 760 CONTINUE
 IF (start == icst) GO TO 1000
 IF (apprch == -1 .AND. start == imst) GO TO 1000
 
!     INITIALIZE VPS TABLE FOR RESTART
!     GET FIRST ENTRY IN CHECKPOINT DICTIONARY
 
 ptdtop = 1
 ASSIGN 770 TO irturn
 GO TO 1090
 770 i = ptdtop
 IF (ptdic(ptdtop) /= nxvps) GO TO 1000
 
!     FIRST ENTRY IN CHECKPOINT DICTIONARY IS XVPS - GET FILE OFF OF OLD
!     PROBLEM TAPE, OPTP
 
 CALL OPEN (*1910,iop,ibufr(iopbuf),2)
 
!     CHECK TO SEE IF OLD RESTART TAPE HAS PVT  J = 0 WITHOUT PVT
 
 j = andf(maskhi,ptdic(ptdtop+2)) - (andf(maskhi,ptfct)+1)
 ptfct = ptfct + j
 CALL skpfil (iop,j)
 CALL READ (*2060,*780,iop,vps(3),2,1,iopwrd)
 780 IF (vps(3) /= nxvps .OR. vps(4) /= nblank) GO TO 2060
 j = vps(1)
 CALL READ (*2060,*790,iop,vps,j,1,iopwrd)
 790 CALL skpfil (iop,1)
 CALL CLOSE  (iop,2)
 ptfct  = ptfct + 1
 vps(1) = j
 
!     FOR RESTART COMPARE PVT VALUES WITH VPS VALUES. IF NOT EQUAL SET
!     MODFLG INVPS ENTRY.
 
 IF (pvt(2) <= 2) GO TO 850
 i = 3
 800 j = 3
 810 IF (pvt(2) < j) GO TO 840
 IF (pvt(j) == vps(i) .AND. pvt(j+1) == vps(i+1)) GO TO 820
 jj = andf(pvt(j+2),nosgn)
 j  = j + itype(jj) + 3
 GO TO 810
 
!     FOUND VARIABLE IN PVT TABLE
 
 820 l = andf(vps(i+2),maskhi)
 pvt(j+2) = orf(pvt(j+2),isgnon)
 DO  ll = 1,l
   ii = i + ll + 2
   jj = j + ll + 2
   vps(i+2) = orf(vps(i+2),modflg)
   vps(ii ) = pvt(jj)
 END DO
 840 i = i + andf(vps(i+2),maskhi) + 3
 IF (i < vps(2)) GO TO 800
 850 i = lbd + lcc + 1
 iparpt = mjmsk(i)
 iparw1 = (iparpt-1)/31 + 1
 iparw2 = lbd
 iparbt = MOD(iparpt-1,31) + 2
 idelet = 0
 DO  j1 = iparw1,iparw2
   IF (mjmsk(j1) /= 0) GO TO 870
 END DO
 idelet = 1
 GO TO 1000
 870 DO  j1 = iparw1,iparw2
   IF (mjmsk(j1) == 0) GO TO 910
   DO  i1 = iparbt,32
     IF (andf(mjmsk(j1),two(i1)) == 0) CYCLE
     nampt = 2*(31*(j1-1)+i1-1) - 1
     i2 = 3
     880 IF (mjcd(nampt) /= vps(i2) .OR. mjcd(nampt+1) /= vps(i2+1))  &
         GO TO 890
     IF (andf(vps(i2+2),two(2)) /= 0) CYCLE
     vps(i2  ) = nblank
     vps(i2+1) = nblank
     CYCLE
     890 i2 = i2 + andf(vps(i2+2),maskhi) + 3
     IF (i2 < vps(2)) GO TO 880
   END DO
   910 iparbt = 2
 END DO
 
!     DMAP SEQUENCE COMPILATION - PHASE 1
!     ***********************************
 
!     GENERATE OSCAR
!     POSITION NEW PROBLEM TAPE AT ALTER FILE IF IT EXISTS
 
 1000 IF (alter(1) == 0) GO TO 1030
 nam1 = nxaltr
 nam2 = nblank
 CALL OPEN (*1900,npt,ibufr(nptbuf),0)
 1010 CALL skpfil (npt,1)
 CALL READ (*1950,*1020,npt,icf,2,1,nptwrd)
 1020 IF (icf(1) /= nxaltr) GO TO 1010
 
!     ALTER FILE FOUND - INITIALIZE ALTER CELLS
 
 CALL READ (*1950,*1950,npt,alter,2,1,nptwrd)
 1030 CALL OPEN (*2010,nscr,ibufr(idpbuf),0)
 CALL xgpimw (1,1,0,0)
 CALL xosgen
 IF (start == icst) GO TO 1050
 DO  i = 1,nmask
   medmsk(i) = 0
 END DO
 1050 IF (alter(1) == 0) GO TO 1060
 CALL CLOSE (npt,2)
 1060 CONTINUE
 IF (pvt(2) <= 2) GO TO 1080
 j = 5
 1070 IF (pvt(2) < j) GO TO 1080
 IF (pvt(j) >= 0) CALL xgpidg (-54,0,pvt(j-2),pvt(j-1))
 jj = andf(pvt(j),nosgn)
 j  = j + itype(jj) + 3
 GO TO 1070
 1080 IF (nogo == 2) GO TO 2210
 CALL CLOSE (nscr,1)
 IF (start /= icst) CALL xgpimw (2,0,0,0)
 CALL xgpimw (1,0,0,0)
 
!     ALLOW MINIMAL SIZE FOR PTDIC ARRAY IN OPEN CORE.
!     WE WILL EXPAND IF THIS IS RESTART.
 
 ptdtop = oscar(osbot) + osbot
 ptdbot = ptdtop
 lptdic = 3
 ASSIGN 1130 TO irturn
 
 1090 IF (start == icst) GO TO 1130
 
!     RESTART - LOAD OLD PROBLEM TAPE DICTIONARY INTO OPEN CORE.
 
 CALL OPEN (*1900,npt,ibufr(nptbuf),0)
 
!     FIND XPTDIC ON NEW PROBLEM TAPE
 
 nam1 = nxptdc(1)
 nam2 = nxptdc(2)
 1100 CALL skpfil (npt,1)
 CALL READ (*1950,*1110,npt,ptdic(ptdtop),2,1,nptwrd)
 1110 IF (ptdic(ptdtop) ==     nxcsa) GO TO 1950
 IF (ptdic(ptdtop) /= nxptdc(1)) GO TO 1100
 
!     FOUND XPTDIC
 
 lptdic = loscar - ptdtop
 CALL READ (*1950,*1120,npt,ptdic(ptdtop),lptdic,1,nptwrd)
 GO TO 2030
 1120 ptdbot = ptdtop + nptwrd - 3
 CALL CLOSE (npt,1)
 GO TO irturn, (770,1130)
 
!     IF BOTH DIAGS 14 AND 20 ARE ON, TERMINATE JOB
 
 1130 IF (diag14 /= 1) GO TO 1200
 CALL sswtch (20,i)
 IF (i == 0) GO TO 1200
 WRITE  (optape,1140)
 1140 FORMAT (//' *** JOB TERMINATED BY DIAG 20',//)
 CALL pexit
 
!     DMAP SEQUENCE COMPILATION - PHASE 2
!     ***********************************
 
!     COMPUTE NTU AND LTU FOR DATA SETS IN OSCAR
 
 1200 IF (nogo /= 0 .AND. start /= icst .AND. ptdtop == ptdbot) GO TO 2210
 CALL xflord
 IF (diag14 == 11) GO TO 2120
 IF (nogo /= 0 .OR. lnogo) GO TO 2210
 IF (diag4 /= 0) CALL dumper
 
!     PURGE ALL FILES IN FIAT TABLE THAT HAVE NOT BEEN GENERATED BY
!     IFP SUBROUTINE
 
 i = ifiat(1)*icfiat - 2
 DO  k = 4,i,icfiat
   IF (ifiat(k+1) == 0) GO TO 1210
   IF (ifiat(k+3) /= 0 .OR. ifiat(k+4) /= 0 .OR. ifiat(k+5) /= 0) CYCLE
   IF (icfiat == 11  .AND. (ifiat(k+8) /= 0 .OR. ifiat(k+9) /= 0 .OR.  &
       ifiat(k+10) /= 0)) CYCLE
   
!     FILE NOT GENERATED - PURGE IT.
   
   k1 = ifiat(3)*icfiat + 4
   ifiat(3)  = ifiat(3) + 1
   ifiat(k1) = orf( andf(ifiat(k),masklo),maskhi)
   ifiat(k ) = andf(ifiat(k),orf(maskhi,losgn))
   ifiat(k1+1) = ifiat(k+1)
   ifiat(k1+2) = ifiat(k+2)
   
!     MAKE SURE NO RESIDUE LEFT IN FIAT TABLE
   
   1210 j1 = k + 1
   j2 = k + icfiat - 1
   DO  j = j1,j2
     ifiat(j) = 0
   END DO
   
 END DO
 
!     WRITE OSCAR ON DATA POOL FILE.
 
!     PUT OSCAR NAME IN DPL AND ASSIGN FILE NO.
 
 lstdpl = lstdpl + 1
 i = lstdpl*3 + 1
 dpl(i  ) = ioshdr(1)
 dpl(i+1) = ioshdr(2)
 dpl(i+2) = ndpfil
 ndpfil   = 1 + ndpfil
 
!     WRITE OSCAR HEADER RECORD
!     POSITION FILE
 
 IF (ndpfil == 2) GO TO 1240
 CALL OPEN (*1940,idp,ibufr(idpbuf),0)
 CALL skpfil (idp,ndpfil-2)
 CALL CLOSE  (idp,2)
 1240 idpfct = ndpfil - 1
 CALL OPEN (*1940,idp,ibufr(idpbuf),3)
 CALL WRITE (idp,ioshdr,2,1)
 
!     IF CHECKPOINT AND RESTART FLAGS ARE ON INSERT CHECKPOINT ENTRY IN
!     OSCAR TO SAVE FILES LISTED IN ICPDPL TABLE
 
 IF (start == icst) GO TO 1290
 IF (icpbot >= icptop .AND. icpflg /= 0) GO TO 1250
 cpntry(6) = 1
 CALL WRITE (idp,cpntry,6,1)
 GO TO 1270
 
!     CHECKPOINT ALL FILES LISTED IN ICPDPL
 
 1250 cpntry(7) = (icpbot - icptop + 3)/3
 cpntry(1) = 7 + cpntry(7)*2
 
!     FOR UNMODIFIED RESTART - DMAP SEQUENCE NO. OF THIS INITIAL
!     CHECKPOINT MUST = REENTRY POINT - 1
 
 IF (start == iunst) cpntry(6) = orf(isgnon,rshift(andf(seqno,masklo),16)-1)
 CALL WRITE (idp,cpntry,7,0)
 DO  i = icptop,icpbot,3
   CALL WRITE (idp,icpdpl(i),2,0)
 END DO
 CALL WRITE (idp,0,0,1)
 
!     FOR RESTART - INSERT JUMP IN OSCAR TO POSITION OSCAR AT CORRECT
!     REENTRY POINT
!     FOR MODIFIED RESTART - START AT FIRST EXECUTABLE MODULE
 
 1270 IF (start == imst) jmp(6) = 1
 
!     SEE IF RE-ENTRY POINT IS WITHIN BOUNDS UNLESS SOLUTION CHANGED.
 
 IF (andf(seqno,masklo) == masklo) GO TO 1280
 i = andf(seqno,maskhi)
 IF (i > oscar(osbot+1) .OR. i == 0) GO TO 2110
 jmp(7) = lshift(i,16)
 1280 CALL WRITE (idp,jmp,7,1)
 1290 ospnt = 1
 
!     WRITE NEXT OSCAR ENTRY ON DATA POOL TAPE
 
 1300 CALL WRITE (idp,oscar(ospnt),oscar(ospnt),1)
 IF (oscar(ospnt+3) == ixtim) GO TO 1330
 i = andf(oscar(ospnt+2),maskhi)
 IF (i > 2 .OR. oscar(ospnt+5) >= 0) GO TO 1340
 
!     MAKE SURE SYSTEM HAS ENOUGH FILES AVAILABLE TO HANDLE MODULE
!     REQUIREMENTS.
!     COUNT NUMBER OF I/P AND O/P FILES NEEDED
 
 j1 = 2
 IF (i == 2) j1 = 1
 k = 0
 l = ospnt + 6
 DO  j = 1,j1
   l2 = oscar(l)*3 - 2 + l
   l1 = l + 1
   IF (oscar(l1-1) == 0) GO TO 1320
   DO  l = l1,l2,3
     IF (oscar(l) /= 0) k = k + 1
   END DO
   1320 l = l2 + 3
 END DO
 
!     ADD ON NUMBER OF SCRATCH FILES NEEDED
 
 k = k + oscar(l)
 IF (ifiat(1) < k) GO TO 2070
 GO TO 1340
 
!     OSCAR ENTRY IS XTIME, COMPUTE ROUGH TIME ESTIMATES FOR MODULES IN
!     TIME SEGMENT, AND
!     WRITE XTIME HEADER AND TIME ESTIMATES ONTO DATA POOL
!     (THIS SECTION TEMPORARILY OMITTED)
 
 1330 GO TO 1340
 
!     INCREMENT OSPNT AND CHECK FOR END OF OSCAR
 
 1340 ospnt = ospnt + oscar(ospnt)
 IF (ospnt-osbot > 0.0) THEN
   GO TO  1350
 ELSE
   GO TO  1300
 END IF
 1350 CALL eof (idp)
 IF (start == icst) GO TO 1800
 
 
!     *** RESTART ***
 
 IF (icpbot < icptop) GO TO 1800
 
!     LIST ICPDPL CONTENTS
 
 CALL xgpimw (8,icptop,icpbot,icpdpl)
 
!     ELIMINATE PURGED FILES FROM ICPDPL
 
 i1 = icptop
 DO  i = i1,icpbot,3
   IF (andf(icpdpl(i+2),maskhi) /= 0) GO TO 1410
   icptop = icptop + 3
 END DO
 1410 IF (icpbot < icptop) GO TO 1800
 CALL CLOSE (idp,2)
 ib1s   = idpbuf
 idpbuf = icpbot + 3
 iopbuf = idpbuf + ibufsz
 CALL gopen (idp,ibufr(idpbuf),3)
 
!     TRANSFER CHECKPOINT INFO FROM OLD PROBLEM TAPE TO DATA POOL TAPE
 
 k = lstdpl*3 + 4
 CALL OPEN (*1910,iop,ibufr(iopbuf),2)
 DO  i = icptop,icpbot,3
   dpl(k+2) = 0
   IF (andf(icpdpl(i+2),noflgs) > ptfct) GO TO 1420
   
!     FILE IS EQUIVALENCED TO PREVIOUS ENTRY IN DPL
   
   ndpfil = ndpfil - 1
   dpl(k+2) = dpl(k-1)
   GO TO 1570
   
!     MAKE SURE CORRECT REEL IS MOUNTED FOR OLD PROBLEM TAPE
   
   1420 IF (andf(andf(noflgs,masklo),icpdpl(i+2)) == andf(masklo,ptfct))  &
       GO TO 1480
   
!     ** NEW REEL NEEDED **
!     MOUNT REEL SPECIFIED BY ICPDPL ENTRY
   
   otapid(6) = rshift(andf(noflgs,icpdpl(i+2)),16)
   wrngrl = 0
   
!     SEND OPERATOR MESSAGE
   
   1430 CALL xeot (iop,rshift(ptfct,16),otapid(6),ibufr(iopbuf))
   CALL OPEN (*1910,iop,ibufr(iopbuf),0)
   CALL READ (*2050,*1440,iop,ibf,libf,0,iopwrd)
   
!     SEE THAT CORRECT REEL HAS BEEN MOUNTED.
   
   1440 DO  ii = 1,6
     IF (otapid(ii) /= ibf(ii)) GO TO 1460
   END DO
   GO TO 1470
   1460 wrngrl = wrngrl + 1
   IF (wrngrl < 2) GO TO 1430
   GO TO 2100
   
!     CORRECT REEL MOUNTED - CARRY ON
   
   1470 CALL skpfil (iop,1)
   ptfct = lshift(otapid(6),16) + 1
   IF (filcon == 0.0) THEN
     GO TO  1480
   ELSE
     GO TO  1560
   END IF
   
!     WRITE FILE ON DATA POOL
   
   1480 CALL skpfil (iop,andf(maskhi,icpdpl(i+2))-(andf(maskhi,ptfct)+1))
   
!     CHECK FOR CORRECT FILE
   
!     5 OR 8 WORDS (DEPEND ON ICFIAT VALUE OF 8 OR 11) WRITTEN TO IOP
!     BY XCHK OF PREVIOUS CHECKPOINT RUN.
!     IF ICFIAT=11, READ 5 WORDS HERE FIRST, AND CHECK IF THERE ARE 3
!     MORE WORDS BEHIND.  I.E. OPTP MAY BE WRITTEN WITH A 5-WORD RECORD
!     IF ICFIAT= 8, READ 5 WORDS
   
   IF (icfiat == 11) GO TO 1490
   CALL READ (*2050,*2050,iop,ibf,5,1,iopwrd)
   ibf(8) = 0
   GO TO 1510
   1490 ibf(8) = -999
   CALL READ (*2050,*2050,iop,ibf(1),5,0,iopwrd)
   CALL READ (*2050,*1510,iop,ibf(6),3,1,iopwrd)
   
   1510 DO  ii = i,icpbot,3
     IF (ibf(1) == icpdpl(ii) .AND. ibf(2) == icpdpl(ii+1)) GO TO 1530
   END DO
   GO TO 2050
   
!     A 5-WORD RECORD READ, EXPANDED (THE TRAILERS) TO 8 WORDS
   
   1530 IF (ibf(8) /= -999) GO TO 1540
   ibf(8) = andf(ibf(5),65535)
   ibf(7) = rshift(ibf(5),16)
   ibf(6) = andf(ibf(4),65535)
   ibf(5) = rshift(ibf(4),16)
   ibf(4) = andf(ibf(3),65535)
   ibf(3) = rshift(ibf(3),16)
   
!     COPY FILE TO POOL
   
   1540 CALL WRITE (idp,ibf,icfiat-3,1)
   1560 CALL cpyfil (iop,idp,ibf,libf,iopwrd)
   dpl(k+2) = dpl(k+2) + iopwrd/1000 + 1
   
!     FILE ALL ON DATA POOL TAPE
   
   CALL eof (idp)
   filcon = 0
   
!     MAKE DPL ENTRY FOR ICPDPL ENTRY
   
   dpl(k+2) = orf(orf(lshift(dpl(k+2),16),ndpfil), andf(icpdpl(i+2),ieqflg))
   1570 dpl(k  ) = icpdpl(i  )
   dpl(k+1) = icpdpl(i+1)
   IF (l8 /= 0) CALL conmsg (dpl(k),2,0)
   k = k + 3
   ndpfil = ndpfil + 1
   lstdpl = 1 + lstdpl
   IF (lstdpl > maxdpl) GO TO 2040
   ptfct = andf(noflgs,icpdpl(i+2))
 END DO
 
!     FILES ALL COPIED OVER FROM OLD PROBLEM TAPE TO DATA POOL TAPE.
 
 CALL CLOSE (iop,1)
 
!     SEE IF XVPS IS ON POOL TAPE
 
 k = lstdpl*3 + 1
 l = ndpfil
 IF (dpl(k) /= nxvps) GO TO 1590
 
!     VPS FILE IS LAST ENTRY IN DPL - DELETE ENTRY
 
 lstdpl = lstdpl - 1
 ndpfil = ndpfil - 1
 j = k
 GO TO 1620
 
!     VPS FILE IS NOT LAST ENTRY IN DPL - SEARCH DPL FOR IT
 
 1590 DO  j = 4,k,3
   IF (dpl(j) == nxvps) GO TO 1610
 END DO
 
!     NO RESTART VPS TABLE
 
 GO TO 1800
 
!     XVPS FOUND - ZERO NAME WHEN NOT LAST ENTRY IN DPL
 
 1610 dpl(j  ) = 0
 dpl(j+1) = 0
 
!     XVPS FILE FOUND IN DPL - POSITION POOL TAPE AND INITIALIZE
!     VPS TABLE WITH CHECKPOINT VALUES
 
 1620 CALL CLOSE (idp,3)
 CALL OPEN (*1940,idp,ibufr(idpbuf),2)
 CALL skpfil (idp,andf(dpl(j+2),maskhi)-l-1)
 nam1 = nxvps
 nam2 = nblank
 CALL skpfil (idp,1)
 CALL READ (*1930,*1630,idp,ibf,libf,1,idpwrd)
 1630 IF (ibf(1) /= nxvps) GO TO 1930
 CALL READ (*1930,*1640,idp,ibf,libf,1,idpwrd)
 
!     COMPARE RESTART PARAMETER NAMES WITH VPS NAMES
 
 1640 k = 3
 1650 j = 3
 IF (andf(vps(k+2),modflg) == modflg) GO TO 1730
 1660 IF (ibf(2) < j) GO TO 1730
 IF (ibf(j) == vps(k) .AND. ibf(j+1) == vps(k+1)) GO TO 1670
 j = j + ibf(j+2) + 3
 GO TO 1660
 
!     PARAMETER NAMES MATCH AND MODFLG NOT ON - INITIALIZE VPS WITH
!     RESTART VALUE.
 
 1670 l = ibf(j+2)
 IF (idelet == 1) GO TO 1710
 iparbt = MOD(iparpt-1,31) + 2
 DO  jjj = iparw1,iparw2
   IF (mjmsk(jjj) == 0) GO TO 1690
   DO  iii = iparbt,32
     IF (andf(mjmsk(jjj),two(iii)) == 0) CYCLE
     nampt = 2*(31*(jjj-1) + iii - 1) - 1
     IF (mjcd(nampt) == vps(k) .AND. mjcd(nampt+1) == vps(k+1)) GO TO 1730
   END DO
   1690 iparbt = 2
 END DO
 1710 DO  m = 1,l
   j1 = m + 2 + j
   k1 = m + 2 + k
   vps(k1) = ibf(j1)
 END DO
 
!     CLEAR FLAGS AND TYPE CODE IN VPS ENTRY AND GET NEXT ENTRY.
 
 1730 vps(k+2) = andf(vps(k+2),maskhi)
 k = k + vps(k+2) + 3
 IF (k < vps(2)) GO TO 1650
 
!     FOR UNMODIFIED RESTART LOAD CEITBL FROM LAST CHECKPOINT
 
 CALL READ (*1930,*1740,idp,ibf,libf,1,idpwrd)
 1740 IF (start == imst) GO TO 1770
 k1 = ceitbl(2)
 j1 = ibf(2)
 
!     FOR RESTART INITIALIZE REPT LOOP COUNTS WITH CHECKPOINT INFO
 
 DO  j = 3,j1,4
   DO  k = 3,k1,4
     IF (ceitbl(k+2) == ibf(j+2) .AND. ceitbl(k+3) == ibf(j+3) .AND.  &
         ibf(j+2) /= 0) ceitbl(k+1) = ibf(j+1)
   END DO
 END DO
 
!     FOR BOTH MOD AND UNMOD RESTART - LOAD VARIOUS CELLS IN /SYSTEM/
!     WITH LAST CHECKPOINT INFO
 
 1770 CALL READ (*1790,*1780,idp,ibf,libf,1,idpwrd)
 1780 mpc  = ibf(5)
 spc  = ibf(6)
 load = ibf(8)
 1790 CONTINUE
 CALL CLOSE (idp,1)
 idpbuf = ib1s
 
 
!     POSITION DATA POOL TAPE AT FIRST OSCAR ENTRY
 
 1800 CALL CLOSE (idp,1)
 
!     *** FIRST, PRODUCE DMAP XREF IF REQUESTED
 
 CALL OPEN (*1940,idp,ibufr(idpbuf),2)
 CALL skpfil (idp,idpfct-1)
 CALL fwdrec (*1920,idp)
 IF (andf(iflg(5),1) /= 0) CALL oscxrf (idpfct-1,idpbuf-1)
 CALL CLOSE (idp,2)
 
!     WRITE VPS TABLE ON NEW PROBLEM TAPE IF CHECKPOINT FLAG ES SET
!     CLEAR FLAGS IN VPS
 
 k = 3
 1810 vps(k+2) = andf(vps(k+2),maskhi)
 k = k + vps(k+2) + 3
 IF (k < vps(2)) GO TO 1810
 IF (icpflg == 0) GO TO 1820
 
!     POSITION TAPE FOR WRITING XVPS
 
 CALL OPEN (*1900,npt,ibufr(nptbuf),0)
 CALL skpfil (npt,andf(nrlfl,maskhi)-1)
 CALL CLOSE (npt,2)
 CALL OPEN (*1900,npt,ibufr(nptbuf),3)
 ibf(1) = nxvps
 ibf(2) = nblank
 CALL WRITE (npt,ibf,2,1)
 CALL WRITE (npt,vps,vps(2),1)
 
!     WRITE CEITBL TABLE ON NEW PROBLEM TAPE
 
 CALL WRITE (npt,ceitbl,ceitbl(2),1)
 CALL eof (npt)
 CALL CLOSE (npt,2)
 
!     INITIALIZE CHECKPOINT PARAMETERS FOR XCHK AND XCEI ROUTINES
 
 ptdic(ptdtop  ) = nxvps
 ptdic(ptdtop+1) = nblank
 ptdic(ptdtop+2) = nrlfl
 nrlfl = nrlfl + 1
 seqno = 1
 
!     WRITE NEW DICTIONARY ON XPTD
 
 CALL OPEN  (*1900,nxptdc,ibufr(nptbuf),1)
 CALL WRITE (nxptdc,nxptdc,2,1)
 CALL WRITE (nxptdc,nrlfl, 2,1)
 CALL WRITE (nxptdc,ptdic(ptdtop),3,1)
 CALL CLOSE (nxptdc,1)
 
!     PUNCH DICTIONARY ENTRY FOR XVPS TABLE
 
 nfile = andf(maskhi,ptdic(ptdtop+2))
 1820 CONTINUE
 IF (nogo /= 0 .OR. lnogo) GO TO 2210
 CALL xgpimw (9,nfile,icpflg,ifiat)
 cppgct = pagect
 IF (iflg(1) == 0) CALL pexit
 
!     TERMINATE RUN IF ANY OF THE DIAG (17, 25, 28, OR 30) AND DIAG 20
!     ARE REQUESTED SIMULTANEOUSLY
 
 CALL sswtch (20,j)
 IF (j == 0) RETURN
 CALL sswtch (28,i)
 CALL sswtch (30,j)
 IF (diag17+diag25+i+j == 0) RETURN
 WRITE  (optape,1830)
 1830 FORMAT (10X,'JOB TERMINATED BY DIAG 20')
 CALL pexit
 
!     E R R O R    M E S S A G E S
 
!     UNEXPECTED END OF TAPE ON NEW PROBLEM TAPE
 
 1900 CALL xgpidg (28,0,0,0)
 GO TO 2200
 
!     UNEXPECTED END OF TAPE ON OLD PROBLEM TAPE
 
 1910 CALL xgpidg (29,0,0,0)
 GO TO 2200
 
!     CANNOT FIND FILE ON DATA POOL TAPE
 
 1920 nam1 = ioshdr(1)
 nam2 = ioshdr(2)
 1930 CALL xgpidg (24,nam1,nam2,0)
 GO TO 2200
 
!     UNEXPECTED END OF TAPE ON DATA POOL TAPE
 
 1940 CALL xgpidg (30,0,0,0)
 GO TO 2200
 
!     CONTROL FILE INCOMPLETE OR MISSING ON NEW PROBLEM TAPE.
 
 1950 CALL xgpidg (31,nam1,nam2,0)
 GO TO 2200
 
!     MED TABLE RECORD MISSING ON SCRATCH FILE
 
 1960 CALL xgpidg (69,nxgpi(1),nxgpi(2),0)
 GO TO 2200
 
!     CARD OR FILE NAME TABLE RECORD MISSING ON SCRATCH FILE
 
 1970 CALL xgpidg (70,nxgpi(1),nxgpi(2),jtype)
 GO TO 2200
 
!     ILLEGAL NUMBER OF WORDS IN MED TABLE RECORD
 
 1980 CALL xgpidg (71,lmed,0,0)
 GO TO 2200
 
!     ILLEGAL NUMBER OF WORDS IN CARD OR FILE NAME TABLE RECORD
 
 1990 CALL xgpidg (72,lmed,jtype,0)
 GO TO 2200
 
!     ILLEGAL BIT NUMBERS IN CARD OR FILE NAME TABLE
 
 2000 CALL xgpidg (73,jtype,0,0)
 GO TO 2200
 
!     SCRATCH FILE CONTAINING DMAP DATA COULD NOT BE OPENED
 
 2010 CALL xgpidg (33,nxgpi(1),nxgpi(2),0)
 GO TO 2200
 
!     PVT TABLE OVERFLOW
 
 2020 CALL xgpidg (14,npvt,nblank,0)
 GO TO 2210
 
!     XPTDIC OVERFLOWED
 
 2030 CALL xgpidg (14,nxptdc(1),nxptdc(2),0)
 GO TO 2200
 
!     DPL TABLE OVERFLOW
 
 2040 CALL xgpidg (14,ndpl,nblank,0)
 GO TO 2210
 
!     CANNOT FIND FILE ON OLD PROBLEM TAPE
 
 2050 CALL xgpidg (36,icpdpl(i),icpdpl(i+1),0)
 GO TO 2210
 2060 CALL xgpidg (36,ptdic(i),ptdic(i+1),0)
 GO TO 2210
 
!     NOT ENOUGH FILES AVAILABLE FOR MODULE REQUIREMENTS.
 
 2070 CALL xgpidg (-37,ospnt,k,ifiat(1))
 GO TO 1340
 
!     NOT ENOUGH CORE FOR GPI TABLES
 
 2080 CALL xgpidg (38,-loscar,0,0)
 GO TO 2200
 
!     MED TABLE OVERFLOW
 
 2090 CALL xgpidg (14,nmed,nblank,0)
 GO TO 2200
 
!     INCORRECT OLD PROBLEM TAPE MOUNTED
 
 2100 CALL xgpidg (35,0,0,0)
 GO TO 2210
 
!     REENTRY POINT NOT WITHIN BOUNDS
 
 2110 CALL xgpidg (46,0,0,0)
 GO TO 2210
 
!     USER DMAP ALTER CONTAINS ERROR, DIAG 14 FLAG IS NOT REQUESTED, AND
!     ECHO IS NOT 'NONO', PRINT RIGID FORMAT BEFORE QUITTING
 
 2120 IF (iecho /= -2) CALL xgpimw (13,0,0,core)
 GO TO 2210
 
!     TERMINATE JOB IF NOGO = 1
 
 2200 nogo = 2
 2210 WRITE  (optape,2220)
 2220 FORMAT (//5X,'*** JOB TERMINATED DUE TO ABOVE ERRORS')
 CALL mesage (-37,0,nxgpi)
 RETURN
END SUBROUTINE xgpi
