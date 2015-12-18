SUBROUTINE xflord
     
!     THE PURPOSE OF THIS ROUTINE IS TO COMPUTE THE LTU (LAST TIME USED)
!     VALUE AND THE NTU (NEXT TIME USED) VALUE FOR THE INPUT AND OUTPUT
!     FILE SECTIONS OF THE OSCAR ENTRIES.
 
!          ... DESCRIPTION OF PROGRAM VARIABLES ...
!     LPTOP  = POINTER/SEQUENCE NUMBER OF FIRST ENTRY IN A DMAP LOOP.
!     LPBOT  = LAST ENTRY IN A LOOP.
!     IOPNT  = POINTER TO FILE NAME IN I/O SECTION OF OSCAR ENTRY.
!     LPORD  = POINTER TO IORDNL TABLE ENTRY CORRESPONDING TO LPTOP.
!     IORDNO = FILE ORDINAL NUMBER
 
 IMPLICIT INTEGER (a-z)
 EXTERNAL        lshift,rshift,andf,orf,complf
 DIMENSION       ptdic(1),xnam(12),itmp(1),iordnl(1800),icpdpl(1),  &
     oscar(2),os(5)
 COMMON /system/ bufsz,optape,nogo,dum1(20),icfiat,dum2(57),icpflg
 COMMON /xgpi4 / irturn,insert,iseqn,dmpcnt,  &
     idmpnt,dmppnt,bcdcnt,length,icrdtp,ICHAR,newcrd, modidx,ldmap,isavdw,dmap(1)
 COMMON /xgpic / icold,islsh,iequl,nblank,nxequi,  &
     ndiag,nsol,ndmap,nestm1,nestm2,nexit,  &
     nbegin,nend,njump,ncond,nrept,ntime,nsave,noutpt, nchkpt,npurge,nequiv,  &
     ncpw,nbpc,nwpc, maskhi,masklo,isgnon,nosgn,iallon,masks(1)
 COMMON /zzzzzz/ core(1)
 COMMON /xgpi2 / lmpl,mplpnt,mpl(1)
 COMMON /xdpl  / dpl(3)
 COMMON /xgpi5 / iapp,start,alter(2),sol,subset,iflag,  &
     iestim,icftop,icfpnt,lctlfl,ictlfl(1)
 COMMON /xgpi6 / medtp,dum5(5),diag14
 COMMON /xgpi7 / fpnt,lfile,FILE(1)
 COMMON /xgpi8 / icptop,icpbot,lcpdpl
 COMMON /xfiat / ifiat(3)
 COMMON /xfist / ifist(1)
 COMMON /xgpid / icst,iunst,imst,ihapp,idsapp,idmapp,  &
     isave,itape,iappnd,intgr,losgn,noflgs
 COMMON /xoldpt/ ptdtop,ptdbot,lptdic,nrlfl,seqno
 COMMON /two   / two(4)
 EQUIVALENCE     (core(1),loscar  , os(1)),  &
     (osprc  ,os(2)  ), (osbot   ,os(3)),  &
     (ospnt  ,os(4)  ), (oscar(1),os(5),ptdic(1)),  &
     (dmap(1),itmp(1)), (oscar(1),icpdpl(1)),  &
     (lmpl   ,lordnl ), (mplpnt,iorbot), (dpl(1) ,ndpfil ), (dpl(2),maxdpl),  &
     (dpl(3) ,lstdpl ), (two(4),reuse )
 DATA    nordn1/ 4HIORD/,  nordn2/ 4HNL  /,  &
     nxvps / 4HXVPS/,  ncpdp1/ 4HICPD/, ncpdp2/4HPL  /,  &
     xnam  / 4HXTIM ,  4HE   , 4HXSAV , 4HE   ,  &
     4HXUOP ,  4H    , 4HXCHK , 4H    , 4HXPUR ,  4HGE  , 4HXEQU , 4HIV  /
 DATA    nthpas/ 0     /,  dlyerr/ 0     /
 
 OR(i,j)  = orf(i,j)
 AND(i,j) = andf(i,j)
 compl(l) = complf(l)
 
!     USE AREA IN OPEN CORE BETWEEN PTDIC AND MED ARRAYS FOR STORING
!     MISSING FILE DATA
 
 iflag  = 0
 icptop = ptdbot + 3
 icpbot = icptop - 3
 lcpdpl = medtp  - icptop
 IF (start == imst) dlyerr = 1
 irentr = AND(maskhi,seqno)
 idmpct = rshift(seqno,16)
 
!     PREPARE FOR NTH PASS THRU OSCAR
!     *******************************
 
 10 IF (nogo > 1) GO TO 960
 ospnt = 1
 osprc = ospnt
 iorbot= 0
 ifeq  = 0
 
!     INCREMENT NUMBER OF PASSES MADE THRU OSCAR
 
 nthpas = 1 + nthpas
 
!     ENTER DPL FILE NAMES IN IORDNL TABLE
 
 i = lstdpl*3 + 1
 idpl = i
 IF (lstdpl == 0) GO TO 30
 DO  k = 4,i,3
   iorbot = iorbot + 4
   iordnl(iorbot  ) = dpl(k  )
   iordnl(iorbot+1) = dpl(k+1)
   iordnl(iorbot+2) = 0
   iordnl(iorbot+3) = 0
 END DO
 
!     ENTER FIAT NAMES IN IORDNL TABLE
 
 30 i = ifiat(3)*icfiat - 2
 DO  k = 4,i,icfiat
   IF (ifiat(k+1) == 0) CYCLE
   iorbot = iorbot + 4
   ifiat(k) = OR(lshift(iorbot,16),AND(ifiat(k),orf(maskhi,losgn)))
   iordnl(iorbot  ) = ifiat(k+1)
   iordnl(iorbot+1) = ifiat(k+2)
   iordnl(iorbot+2) = 0
   iordnl(iorbot+3) = 0
 END DO
 
!     FOR UNMODIFIED RESTART BEGIN OSCAR PROCESSING AT RE-ENTRY POINT IF
!     THIS IS FIRST PASS THRU OSCAR
 
 IF (nthpas > 1) GO TO 60
 IF (start /= iunst .OR. irentr == 0) GO TO 60
 DO  j = 1,irentr
   IF (oscar(ospnt+1) >= irentr) GO TO 60
   osprc = ospnt
   ospnt = ospnt + oscar(ospnt)
 END DO
 
!     GET NEXT OSCAR ENTRY
!     ********************
 
!     BRANCH ON OSCAR ENTRY TYPE IF EXECUTE FLAG IS UP
 
 60 IF (oscar(ospnt+5) >= 0) GO TO 70
 i = AND(oscar(ospnt+2),maskhi)
 lstbot = iorbot
 SELECT CASE ( i )
   CASE (    1)
     GO TO 310
   CASE (    2)
     GO TO 390
   CASE (    3)
     GO TO 520
   CASE (    4)
     GO TO 80
 END SELECT
 
!     GET NEXT OSCAR ENTRY
 
 70 IF (ospnt >= osbot) GO TO 650
 IF (oscar(ospnt+5) < 0 .AND. AND(oscar(ospnt+2),maskhi) <= 2) osprc = ospnt
 ospnt = ospnt + oscar(ospnt)
 GO TO 60
 
!     PROCESS TYPE E OSCAR ENTRY
!     **************************
 
!     BRANCH ON NAME
 
 80 DO  i = 1,11,2
   IF (oscar(ospnt+3) /= xnam(i)) CYCLE
   j = (i+1)/2
   SELECT CASE ( j )
     CASE (    1)
       GO TO 70
     CASE (    2)
       GO TO 70
     CASE (    3)
       GO TO 100
     CASE (    4)
       GO TO 100
     CASE (    5)
       GO TO 70
     CASE (    6)
       GO TO 190
   END SELECT
 END DO
 
!     ENTRY IS XUOP OR XCHK - MAKE SURE FILES HAVE BEEN DEFINED OR
!     PREPURGED.
 
 100 i1 = ospnt + 7
 i2 = oscar(ospnt+6)*2 + i1 - 2
 iopnt = i1
 110 IF (iorbot <= 0) GO TO 130
 DO  j = 4,iorbot,4
   IF (oscar(iopnt) /= iordnl(j) .OR. oscar(iopnt+1) /= iordnl(j+1)) CYCLE
   IF (start /= iunst .OR. j > idpl) GO TO 170
   nnfind = -1
   CALL xfldef (oscar(iopnt),oscar(iopnt+1),nnfind)
   GO TO 170
 END DO
 
!     FILE NOT IN ORDNAL TABLE - SEE IF IT IS IN PREVIOUS PURGE OR
!     EQUIV ENTRY
 
 130 k1 = 2
 k1 = oscar(k1)
 k2 = oscar(ospnt+1) - 1
 kk = 1
 DO  k = k1,k2
   IF (oscar(kk+3) /= xnam(9) .AND. oscar(kk+3) /= xnam(11)) GO TO 160
   
!     PURGE OR EQUIV ENTRY FOUND - SEARCH FOR FILE NAME MATCH
   
   l1 = kk + 7
   l3 = kk + oscar(kk)
   
!     GET FIRST/NEXT FILE LIST
   
   140 l2 = oscar(l1-1)*2 + l1 - 2
   incrlp = 2
   IF (oscar(kk+3) /= xnam(11)) GO TO 145
   l2 = l2 + 1
   incrlp = 3
   145 CONTINUE
   DO  l = l1,l2,incrlp
     IF (oscar(l) == oscar(iopnt) .AND. oscar(l+1) == oscar(iopnt+1))  &
         GO TO 180
     IF (l == l1+incrlp) GO TO 153
   END DO
   GO TO 159
   153 l4 = l1 + incrlp
   incrlp = 2
   l4 = l4 + incrlp
   DO  l = l4,l2,incrlp
     IF (oscar(l) == oscar(iopnt) .AND. oscar(l+1) == oscar(iopnt+1))  &
         GO TO 180
   END DO
   159 l1 = l2 + 4
   IF (l1 < l3) GO TO 140
   160 kk = oscar(kk) + kk
 END DO
 
!     FILE IS NOT PURGED OR DEFINED - SEE IF IT IS ON PROBLEM TAPE
 
 nofind = -1
 GO TO 450
 
!     FILE IS IN ORDNAL TABLE - ENTER RANGE
 
 170 IF (iordnl(j+3) < 0) THEN
   GO TO   180
 END IF
 175 iordnl(j+3) = lshift(oscar(ospnt+1),16)
 180 iopnt = iopnt + 2
 IF (iopnt <= i2) GO TO 110
 GO TO 70
 
!     PROCESS EQUIV INSTRUCTION
 
 190 l1 = ospnt + 7
 nwdh = oscar(ospnt) - 6
 230 ndatab = oscar(l1-1)
 iprime = 0
 DO  khr = 1,ndatab
   
!     CHECK FOR DATA BLOCK IN IORDNL
   
   IF (iorbot <= 0) GO TO 200
   DO  i = 4,iorbot,4
     IF (iordnl(i) /= oscar(l1) .OR. iordnl(i+1) /= oscar(l1+1)) CYCLE
     IF (start /= iunst .OR. i > idpl) GO TO 210
     nnfind = -1
     CALL xfldef (oscar(l1),oscar(l1+1),nnfind)
     GO TO 210
   END DO
   
!     FILE NOT IN IORDNL, SEE IF ON PTDIC OR REGEN
   
   200 IF (start == icst .OR. iprime /= 0) GO TO 220
   nofind = 1
   CALL xfldef (oscar(l1),oscar(l1+1),nofind)
   IF (nofind < 0) THEN
     GO TO    10
   ELSE IF (nofind == 0) THEN
     GO TO   215
   END IF
   220 IF (dlyerr /= 0 .OR. iprime /= 0) GO TO 215
   
!     PRIMARY EQUIV FILE NOT DEFINED
   
   CALL xgpidg (32,ospnt,oscar(l1),oscar(l1+1))
   GO TO 210
   
!     PUT FILE IN IORDNL, FLAG FOURTH WORD FOR EQUIV
   
   215 iorbot = iorbot + 4
   IF (iorbot - lordnl < 0) THEN
     GO TO   225
   ELSE
     GO TO   780
   END IF
   225 iordnl(iorbot  ) = oscar(l1)
   iordnl(iorbot+1) = oscar(l1+1)
   iordnl(iorbot+2) = 0
   iordnl(iorbot+3) = isgnon
   210 IF (iprime /= 0) GO TO 211
   lstuse = AND(maskhi,iordnl(i+2))
   IF (lstuse == 0) GO TO 212
   ntu = OR(oscar(ospnt+1),AND(iordnl(i+2),itape))
   oscar(lstuse) = OR(AND(oscar(lstuse),masklo),ntu)
   212 iordnl(i+2) = OR(oscar(l1+2),AND(iordnl(i+2),itape))
   iordnl(i+3) = lshift(oscar(ospnt+1),16)
   oscar(l1+2) = OR(AND(oscar(l1+2),maskhi),lshift(i,16))
   211 l1 = l1 + 2
   IF (iprime == 0) l1 = l1 + 1
   iprime = 1
 END DO
 nwdh = nwdh - 2*ndatab - 3
 IF (nwdh <= 0) GO TO 70
 l1 = l1 + 2
 GO TO 230
 
!     PROCESS TYPE F OSCAR ENTRY
!     **************************
 
!     SCAN OSCAR OUTPUT FILE SECTION,ENTER NAMES IN IORDNL TABLE.
 
 310 k = ospnt + 6
 k = oscar(k)*3   + 2 + k
 i = oscar(k-1)*3 - 3 + k
 iopnt = k
 ASSIGN 380 TO irturn
 
!     GET FIRST/NEXT FILE NAME FROM OSCAR
 
 320 IF (oscar(iopnt) == 0) GO TO 380
 
!     SEE IF FILE NAME IS ALREADY IN ORDNAL TABLE
 
 IF (iorbot <= 0) GO TO 340
 DO  k = 4,iorbot,4
   IF (iordnl(k) /= oscar(iopnt) .OR. iordnl(k+1) /= oscar(iopnt+1)) CYCLE
   IF (start /= iunst .OR. k > idpl) GO TO 345
   nnfind = -1
   CALL xfldef (oscar(iopnt),oscar(iopnt+1),nnfind)
   GO TO 345
 END DO
 GO TO 340
 345 IF (iordnl(k+3) < 0) THEN
   GO TO   346
 ELSE
   GO TO   820
 END IF
 346 kxt = k
 GO TO 347
 
!     INCREMENT TO NEXT ORDNAL ENTRY AND ENTER FILE NAME AND LU POINTER.
 
 340 iorbot = iorbot + 4
 IF (iorbot - lordnl < 0) THEN
   GO TO   350
 ELSE
   GO TO   780
 END IF
 350 iordnl(iorbot  ) = oscar(iopnt  )
 iordnl(iorbot+1) = oscar(iopnt+1)
 
!     SEE IF TAPE FLAG IS SET FOR THIS FILE
 
 kxt = iorbot
 347 lstuse = iopnt + 2
 IF (AND(oscar(ospnt+2),maskhi) > 2) lstuse=0
 IF (fpnt < 1) GO TO 370
 DO  k = 1,fpnt,3
   IF (oscar(iopnt) == FILE(k) .AND. oscar(iopnt+1) == FILE(k+1))  &
       lstuse = OR(lstuse,AND(FILE(k+2),itape))
 END DO
 370 iordnl(kxt+2) = lstuse
 iordnl(kxt+3) = lshift(oscar(ospnt+1),16)
 
!     IORDNL POINTER  TO OSCAR IF TYPE F OR O FORMAT
 
 IF (AND(oscar(ospnt+2),maskhi) <= 2)  &
     oscar(iopnt+2) = OR(lshift(kxt,16),AND(oscar(iopnt+2),maskhi))
 GO TO irturn, (380,510)
 
!     O/P FILE PROCESSED  -  INCREMENT TO NEXT O/P FILE
 
 380 iopnt = iopnt + 3
 IF (iopnt <= i) GO TO 320
 
!     OUTPUT SECTION SCANNED, NOW SCAN INPUT FILE SECTION OF OSCAR.
 
!     PROCESS TYPE F OR O OSCAR ENTRY
!     *******************************
 
!     SCAN OSCAR INPUT FILE SECTION,ENTER RANGES IN IORDNL TABLE.
 
 390 k = ospnt + 7
 i = oscar(k-1)*3 -3 + k
 iopnt = k
 
!     GET FIRST/NEXT FILE NAME FROM OSCAR
 
 400 IF (oscar(iopnt) == 0) GO TO 510
 nofind = 1
 ASSIGN 510 TO irturn
 
!     NOW SCAN IORDNAL TABLE FOR FILE NAME
 
 j1 = lstbot
 IF (j1 <= 0) GO TO 440
 DO  j = 4,j1,4
   IF (oscar(iopnt) /= iordnl(j) .OR. oscar(iopnt+1) /= iordnl(j+1)) CYCLE
   IF (start /= iunst .OR. j > idpl) GO TO 420
   nnfind = -1
   CALL xfldef (oscar(iopnt),oscar(iopnt+1),nnfind)
   GO TO 420
 END DO
 GO TO 440
 
!     FOUND FILE IN IORDNL TABLE - ENTER NTU AND TAPE FLAG INTO
!     OSCAR ENTRY POINTED TO BY IORDNL ENTRY
 
 420 lstuse = AND(maskhi,iordnl(j+2))
 IF (lstuse == 0) GO TO 430
 ntu = OR(oscar(ospnt+1),AND(iordnl(j+2),itape))
 oscar(lstuse) = OR(AND(oscar(lstuse),masklo),ntu)
 
!     SET RANGE AND LASTUSE POINTER IN IORDNAL ENTRY
 
 430 nofind = -1
 iordnl(j+2) = OR(iopnt+2,AND(iordnl(j+2),itape))
 iordnl(j+3) = lshift(oscar(ospnt+1),16)
 
!     LINK OSCAR I/P FILE TO IORDNL ENTRY
 
 oscar(iopnt+2) = OR(AND(oscar(iopnt+2),maskhi),lshift(j,16))
 
!     I/P FILE PROCESSED - MAKE SURE IT WAS DEFINED
 
 440 IF (nofind < 0) THEN
   GO TO   510
 END IF
 
!     I/P FILE NOT DEFINED
 
 450 IF (start == icold) GO TO 470
 
!     RESTART - SEE IF FILE IS ON PROBLEM TAPE OR CAN BE REGENERATED
!     BY RE-EXECUTING SOME MODULES.
 
 CALL xfldef (oscar(iopnt),oscar(iopnt+1),nofind)
 IF (nofind < 0) THEN
   GO TO    10
 ELSE IF (nofind == 0) THEN
   GO TO   500
 END IF
 
!     ERROR - FILE NOT DEFINED(PUT OUT MESSAGE AT END OF XFLORD)
!     SEE IF FILE IS ALREADY IN ICPDPL TABLE
 
 470 IF (dlyerr /=      0) GO TO 500
 IF (icpbot < icptop) GO TO 490
 DO  l = icptop,icpbot,3
   IF (oscar(iopnt) == icpdpl(l) .AND. oscar(iopnt+1) == icpdpl(l+1))  &
       GO TO 500
 END DO
 
!     ENTER FILE IN ICPDPL TABLE
 
 490 icpbot = icpbot + 3
 IF (icpbot+3-icptop > lcpdpl) GO TO 830
 icpdpl(icpbot  ) = oscar(iopnt  )
 icpdpl(icpbot+1) = oscar(iopnt+1)
 icpdpl(icpbot+2) = -ospnt
 
!     ENTER FILE IN ORDNAL TABLE IF NOT CHKPNT MODULE
 
 500 IF (oscar(ospnt+3) /= xnam(7)) GO TO 340
 GO TO 180
 
!     CHECK FOR ANOTHER I/P FILE
 
 510 iopnt = iopnt + 3
 IF (iopnt <= i) GO TO 400
 
!     INPUT FILE SECTION SCANNED,GET NEXT OSCAR ENTRY.
 
 GO TO 70
 
!     PROCESS TYPE C OSCAR ENTRY
!     **************************
 
!     CHECK FOR LOOPING
 
 520 lptop = rshift(oscar(ospnt+6),16)
 IF ((nexit == oscar(ospnt+3)) .OR. (oscar(ospnt+1) < lptop)) GO TO 70
 
!     FIND BEGINNING OF LOOP AND ADJUST IORDNL RANGES INSIDE LOOP.
 
 lpbot = ospnt
 i     = oscar(ospnt+1)
 ospnt = 1
 j1    = oscar(ospnt+1)
 DO  j = j1,i
   IF (oscar(ospnt+1) == lptop) GO TO 540
   ospnt = oscar(ospnt) + ospnt
 END DO
 540 lptop = ospnt
 
!     LOOP TOP FOUND - IF UNMODIFIED RESTART,EXECUTE ALL MODULES INSIDE
!     LOOP.
 
 IF (oscar(lptop+5) < 0 .OR. start /= iunst) GO TO 570
 
!     MAKE SURE FIRST INSTRUCTION IN LOOP IS NOT CHKPNT
 
 IF (oscar(lptop+3) == xnam(7) .AND. oscar(lptop+4) == xnam(8)) GO TO 790
 
!     EXECUTE FLAGS NOT ALL SET - SET FLAGS AND BEGIN OSCAR SCAN AGAIN
 
 550 j1 = oscar(lptop+1)
 DO  j = j1,i
   IF (oscar(ospnt+3) == xnam(7) .AND. oscar(ospnt+4) == xnam(8)  &
       .AND. icpflg == 0) GO TO 560
   IF (oscar(ospnt+5) < 0) GO TO 560
   IF (iflag == 1) GO TO 5510
   iflag = 1
   CALL page1
   CALL xgpimw (11,idmpct,0,0)
   5510 CALL xgpimw (4,0,0,oscar(ospnt))
   oscar(ospnt+5) = OR(isgnon,oscar(ospnt+5))
   560 ospnt = oscar(ospnt) + ospnt
 END DO
 GO TO 10
 
!     EXTEND RANGE OF FILES DEFINED OUTSIDE OF LOOP IF USED INSIDE LOOP
!     GET FIRST/NEXT OSCAR ENTRY INSIDE LOOP
 
 570 ospnt = lptop
 j1 = oscar(lptop+1)
 j2 = oscar(lpbot+1)
 DO  j = j1,j2
   IF (AND(oscar(ospnt+2),maskhi) > 2) GO TO 640
   
!     GET FIRST/NEXT I/P FILE OF OSCAR ENTRY
   
   k1 = ospnt + 7
   k2 = oscar(k1-1)*3 - 3 + k1
   loop630:  DO  k = k1,k2,3
     IF (oscar(k) == 0) CYCLE loop630
     
!     SEE IF FILE SAVE IS ON
     
     IF (fpnt < 1) GO TO 590
     DO  l = 1,fpnt,3
       IF (oscar(k) /= FILE(l) .OR. oscar(k+1) /= FILE(l+1)) CYCLE
       IF (AND(isave,FILE(l+2)) == isave) GO TO 620
       GO TO 590
     END DO
     
!     FILE SAVE FLAG NOT ON - SEE IF I/P FILE IS GENERATED INSIDE LOOP
     
     590 l1 = oscar(ospnt+1)
     
!     GET FIRST/NEXT OSCAR ENTRY INSIDE LOOP
     
     n = lptop
     DO  l = j1,l1
       IF (AND(oscar(n+2),maskhi) /= 1 .OR. oscar(n+5) >= 0) GO TO 610
       
!     GET FIRST/NEXT O/P FILE
       
       m1 = oscar(n +6)*3 + n + 8
       m2 = oscar(m1-1)*3 - 3 + m1
       DO  m = m1,m2,3
         IF (oscar(m) == 0) CYCLE
         IF (oscar(m) == oscar(k) .AND. oscar(m+1) == oscar(k+1))  &
             CYCLE loop630
       END DO
       610 n = oscar(n) + n
     END DO
     
!     EXTEND I/P FILE RANGE TO END OF LOOP
     
     620 n = rshift(oscar(k+2),16)
     iordnl(n+3) = lshift(i,16)
   END DO loop630
   IF (start /= iunst) GO TO 640
   
!     FOR UNMODIFIED RESTART, MARK ALL OUTPUT FILES WITHIN THE
!     LOOP AND BEFORE THE RE-ENTRY POINT FOR REUSE
   
   kk1 = k1
   IF (oscar(kk1-6) >= irentr) GO TO 640
   IF (AND(oscar(kk1-5),maskhi) /= 1) GO TO 640
   k1 = k2 + 4
   k2 = 3*oscar(k1-1) - 3 + k1
   DO  k = k1,k2,3
     IF (oscar(k) == 0) CYCLE
     nofind = -1
     CALL xfldef (oscar(k),oscar(k+1),nofind)
   END DO
   640 ospnt = oscar(ospnt) + ospnt
 END DO
 
!     LOOP SCANNED, GET NEXT OSCAR ENTRY AFTER LOOP ENTRIES
 
 ospnt = lpbot
 GO TO 70
 
!     OSCAR HAS BEEN PROCESSED
!     ************************
 
 650 IF (dlyerr == 0) GO TO 653
 dlyerr = 0
 GO TO 10
 
!     SET  NTU = LTU FOR LAST REFERENCE TO EACH FILE IN OSCAR.
 
 653 DO  i = 4,iorbot,4
   lstuse = AND(iordnl(i+2),maskhi)
   IF (lstuse == 0) CYCLE
   ntu = OR(AND(itape,iordnl(i+2)),rshift(iordnl(i+3),16))
   oscar(lstuse) = OR(ntu,AND(oscar(lstuse),masklo))
 END DO
 
!     SEARCH FILE TABLE FOR FILES WITH APPEND OR SAVE FLAG UP
 
 IF (fpnt < 1) GO TO 690
 DO  j = 1,fpnt,3
   IF (AND(FILE(j+2),iappnd) == 0 .AND. AND(FILE(j+2),isave) == 0) CYCLE
   
!     FOR RESTART, MARK APPEND AND SAVE FILES FOR REUSE
   
   nofind = -1
   CALL xfldef (FILE(j),FILE(j+1),nofind)
   IF (AND(FILE(j+2),isave) /= 0) CYCLE
   
!     APPEND FLAG SET - FIND CORRESPONDING IORDNL ENTRY AND SET FLAG
   
   DO  i = 4,iorbot,4
     IF (iordnl(i) == FILE(j) .AND. iordnl(i+1) == FILE(j+1))  &
         iordnl(i+3) = OR(iappnd,iordnl(i+3))
   END DO
 END DO
 
!     STORE LTU IN OSCAR FILE ENTRIES
 
 690 ospnt = 1
 700 IF (oscar(ospnt+5) >= 0 .OR. AND(oscar(ospnt+2),maskhi) > 2) GO TO 730
 k = ospnt + 7
 j = 1
 IF (AND(oscar(ospnt+2),maskhi) == 1) j = 2
 DO  l = 1,j
   
   i = oscar(k-1)*3 - 3 + k
   DO  iopnt = k,i,3
     IF (oscar(iopnt) == 0) CYCLE
     j1  = rshift(oscar(iopnt+2),16)
     ltu = AND(oscar(iopnt+2),OR(losgn,maskhi))
     oscar(iopnt+2) = OR(ltu,iordnl(j1+3))
   END DO
   k  = i + 4
 END DO
 730 IF (oscar(ospnt+3) /= xnam(11)) GO TO 735
 i  = oscar(ospnt) - 6
 k  = ospnt + 7
 733 j1 = rshift(oscar(k+2),16)
 ltu= AND(oscar(k+2),OR(losgn,maskhi))
 oscar(k+2) = OR(ltu,iordnl(j1+3))
 i  = i - 2*oscar(k-1) - 3
 IF (i <= 0) GO TO 735
 k  = k + 2*oscar(k-1) + 3
 GO TO 733
 735 ospnt = ospnt + oscar(ospnt)
 IF (ospnt - osbot > 0.0) THEN
   GO TO   740
 ELSE
   GO TO   700
 END IF
 
!     STORE LTU IN FIAT ENTRIES
 
 740 i = ifiat(3)*icfiat - 2
 DO  k = 4,i,icfiat
   IF (ifiat(k+1) == 0) CYCLE
   j = rshift(AND(ifiat(k),masklo),16)
   
!     SEE IF FILE HAS BEEN REFERENCED
   
   IF (AND(iordnl(j+3),compl(iappnd)) /= 0) GO TO 760
   
!     FILE NOT USED - DROP IT FROM FIAT
   
   ifiat(k) = AND(ifiat(k),OR(maskhi,losgn))
   k1 = k + 1
   k2 = k + icfiat - 3
   DO  kk = k1,k2
     ifiat(kk) = 0
   END DO
   CYCLE
   760 ltu = AND(ifiat(k),OR(OR(isgnon,losgn),maskhi))
   ifiat(k) = OR(ltu,iordnl(j+3))
 END DO
 GO TO 840
 
!     ERROR MESSAGES
!     **************
 
!     IORDNL TABLE OVERFLOW
 
 780 CALL xgpidg (14,nordn1,nordn2,AND(oscar(ospnt+5),nosgn))
 GO TO 960
 
!     CHKPNT IS FIRST INSTRUCTION IN LOOP
 
 790 CALL xgpidg (47,lptop,0,0)
 oscar(lptop+5) = OR(oscar(lptop+5),isgnon)
 GO TO 550
 
!     FILE APPEARS MORE THAN ONCE AS OUTPUT
 
!     SUPPRESS MESSAGE ONCE IF FILE IS INITIALLY UNDEFINED
 
 820 IF (icpbot < icptop) GO TO 8220
 DO  ii = icptop,icpbot,3
   IF (oscar(iopnt) /= icpdpl(ii) .OR. oscar(iopnt+1) /= icpdpl(ii+1) ) CYCLE
   IF (icpdpl(ii+2) >= 0) EXIT
   icpdpl(ii+2) = -icpdpl(ii+2)
   GO TO 346
 END DO
 8220 CALL xgpidg (-45,ospnt,oscar(iopnt),oscar(iopnt+1))
 GO TO 346
 
!     ICPDPL TABLE OVERFLOW
 
 830 CALL xgpidg (14,ncpdp1,ncpdp2,0)
 GO TO 960
 
!     CHECK ICPDPL TABLE FOR UNDEFINED FILES
 
 840 IF (icpbot < icptop) GO TO 860
 DO  i = icptop,icpbot,3
   CALL xgpidg (-22,IABS(icpdpl(i+2)),icpdpl(i),icpdpl(i+1))
 END DO
 
!     IF DIAG 14 IS NOT ON, AND THERE ARE UNDEFINED FILES FROM USER'S
!     ALTER (DIAG14 IS SET TO 10 BY XGPI AT THIS TIME), SET DIAG14 TO 11
!     TO FLAG XGPI TO PRINT THE DMAP COMPILE LISTING.
 
!     IF DIAG 14 IS ON, THE DMAP LISTING IS ALREADY PRINTTED BY XSCNDM,
!     SHICH IS CALLED BY XOSGEN. XOSGEN IS CALLED BY XGPI BEFORE THIS
!     XFLORD IS CALLED (ALSO BY XGPI)
 
 IF (diag14 ==  10) diag14 = 11
 IF (start /= icst) GO TO 865
 GO TO 960
 
!     NO UNDEFINED FILES - CHECK FOR RESTART
 
 860 IF (start == icst) GO TO 960
 
!     RESTART - USE LAST XVPS ENTRY IN PTDIC FOR RESTART.
!     EXCLUDE FIRST NXVPS ENTRY
 
 865 ptdtop = ptdtop + 3
 nofind = -1
 CALL xfldef (nxvps,nblank,nofind)
 ptdtop = ptdtop - 3
 
!     OVERLAY PTDIC TABLE WITH ICPDPL TABLE
 
 icptop = ptdtop
 icpbot = icptop - 3
 lcpdpl = lptdic
 
!     SCAN PTDIC FOR REUSE FLAGS
 
 DO  j = ptdtop,ptdbot,3
   IF (AND(ptdic(j+2),reuse) == 0) CYCLE
   
!     REUSE FLAG UP - ENTER FILE IN ICPDPL
   
   icpbot = icpbot + 3
   icpdpl(icpbot  ) = ptdic(j  )
   icpdpl(icpbot+1) = ptdic(j+1)
   icpdpl(icpbot+2) = ptdic(j+2)
 END DO
 
!     ORDER FILES IN ICPDPL BY REEL/FILE NUMBER
 
 IF (icpbot < icptop) GO TO 960
 
!     DO NOT DISTURB EXISTING ORDER
 
 IF (icpbot == icptop) GO TO 900
 k = icptop
 881 l = k
 882 IF (AND(icpdpl(k+2),noflgs) <= AND(icpdpl(k+5),noflgs)) GO TO 890
 
!     SWITCH
 
 DO  m = 1,3
   j = k + m + 2
   itmp(1) = icpdpl(j)
   icpdpl(j) = icpdpl(j-3)
   icpdpl(j-3) = itmp(1)
 END DO
 k = k - 3
 IF (k >= icptop) GO TO 882
 890 k = l + 3
 IF (k < icpbot) GO TO 881
 900 CONTINUE
 
!     ENTER PURGED FILE IN FIAT IF THERE IS NO POSSIBLE WAY TO GENERATE
!     FILE
 
 j1 = 2
 j1 = oscar(j1)
 j2 = oscar(osbot+1)
 loop950:  DO  i = icptop,icpbot,3
   IF (AND(icpdpl(i+2),maskhi) /= 0) EXIT loop950
   ospnt = 1
   DO  j = j1,j2
     IF (AND(maskhi,oscar(ospnt+2)) > 2 .OR. oscar(ospnt+5) >= 0) GO TO 940
     
!     SEE IF PURGED FILE IS IN I/P SECTION
     
     k1 = ospnt + 7
     k2 = oscar(k1-1)*3 - 3 + k1
     DO  k = k1,k2,3
       IF (oscar(k) == icpdpl(i) .AND. oscar(k+1) == icpdpl(i+1)) GO TO 930
     END DO
     
!     PURGED FILE IS NOT IN I/P SECTION - SEARCH O/P SECTION FOR IT.
     
     IF (AND(maskhi,oscar(ospnt+2)) /= 1) GO TO 940
     k1 = oscar(ospnt+6)*3  + ospnt + 8
     k2 = oscar(k1-1)*3 - 3 + k1
     DO  k = k1,k2,3
       IF (oscar(k) == icpdpl(i) .AND. oscar(k+1) == icpdpl(i+1))  &
           CYCLE loop950
     END DO
     GO TO 940
     
!     PURGED FILE FIRST USED AS INPUT - THEREFORE IT CANNOT BE GENERATED
!     ENTER PURGED FILE IN FIAT
     
     930 l = ifiat(3)*icfiat + 4
     ifiat(3  ) = ifiat(3) + 1
     ifiat(l  ) = OR(maskhi,oscar(k+2))
     ifiat(l+1) = oscar(k  )
     ifiat(l+2) = oscar(k+1)
     CYCLE loop950
     940 ospnt = oscar(ospnt) + ospnt
   END DO
 END DO loop950
 960 RETURN
END SUBROUTINE xflord
