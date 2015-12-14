SUBROUTINE xchk
     
!     THE PURPOSE OF THIS ROUTINE IS TO SAVE ON THE NEW PROBLEM NPTP
!     TAPE ALL FILES REQUESTED BY XCHK OSCAR ENTRY TOGETHER WITH ANY
!     OTHER DATA NECESSARY FOR RESTART.
 
!          ... DEFINITION OF PROGRAM VARIABLES ...
!     NPTPNT = POINTER TO GINO BUFFER FOR NEW PROBLEM NPTP TAPE
!     DPPNT  = POINTER TO GINO BUFFER FOR DATA POOL TAPE
!     FPNT   = POINTER TO GINO BUFFER FOR FILES LISTED IN FIAT TABLE
!     IOBUF  = INPUT/OUTPUT BUFFER AREA
!     IOPNT  = POINTER TO IOBUF
!     LIOBUF = LENGTH OF IOBUF
!     DICT   = PRELIMINARY FILE DICTIONARY
!     FDICT  = FINAL FILE DICTIONARY TO BE WRITTEN ON NEW PROBLEM TAPE
!     LDC    = POINTER TO LAST DICT ENTRY MADE.
!     DCPNT  = POINTER TO DICT ENTRY BEING SCANNED.
!     NPTFN  = NEW PROBLEM TAPE (NPTP) FILE NUMBER TO BE ASSIGNED
!     UCBPNT = UCB POINTER FOUND IN FIAT ENTRIES
!     MINFN  = SMALLEST DATA POOL FILE NUMBER
!     DPFCT  = DATA POOL FILE POSITION
!     OSCFN  = DATA POOL FILE NUMBER OF OSCAR FILE
!     EORFLG = END OF RECORD FLAG
!     PURGE  = TABLE OF PURGED CHECKPOINT FILES
!     LPURGE = LENGTH OF PURGE TABLE
!     PRGPNT = POINTER TO LAST PURGE ENTRY
!     REELCT = KEEPS TRACK OF HOW MANY PROBLEM TAPE REELS A FILE IS
!              USING
!     EQFLG  = EQUIVALENCE FLAG
!     DPLFLG = DATA POOL FLAG
!     EOTFLG = END OF TAPE FLAG
!     SETEOR = END OF RECORD FLAG SET
!     FNASS  = NPTP FILE NUMBER ASSIGNED FLAG
!     MASKHI = MASK FOR ALL BITS EXCEPT LOWEST ORDER 16 BITS OF A WORD.
!     NOFLGS = MASK FOR ALL FLAG BITS
!     ALLON  = ALL BITS ON
!     PTDIC  = ARRAY CONTAINING CHECKPOINT DICTIONARY
!     SEQNO  = SEQUENCE NO. OF LAST PTDIC ENTRY THAT WAS PUNCHED OUT.
!     NRLFL  = NEXT REEL/FILE NO. TO BE USED IN PTDIC
!     PTDTOP = POINTER TO FIRST WORD OF FIRST ENTRY IN PTDIC
!     PTDBOT = POINTER TO FIRST WORD OF LAST  ENTRY IN PTDIC
!     LCPTP  = POINTER TO FIRST WORD OF FIRST ENTRY OF NEW GROUP OF
!              ENTRIES TO BE PUT IN PTDIC.
!     LPTDIC = LENGTH (IN WORDS) OF PTDIC
 
 IMPLICIT INTEGER (a-z)
 EXTERNAL        lshift,rshift,andf,orf,complf
 DIMENSION       blkcnt(90),dcparm(2),head(2),purge(100),svfst(2),  &
     pghdg(1),hdg(32),dict(400),fdict(400),ptdic(1),  &
     nxptdc(2),iobuf(1),nxchk(2),nvps(2)
 CHARACTER (LEN=25) :: sfm
 CHARACTER (LEN=29) :: uim
 CHARACTER (LEN=25) :: uwm
 CHARACTER (LEN=23) :: ufm
 COMMON /xmssg / ufm,uwm,uim,sfm
 COMMON /xfist / fist(2)
 COMMON /xpfist/ ipfst
 COMMON /oscent/ oscar(7)
!WKBR COMMON /XCEITB/ CEITBL(2)
 COMMON /xceitb/ ceitbl(42)
 COMMON /xfiat / fiat(3)
 COMMON /xdpl  / dpl(3)
 COMMON /xvps  / vps(2)
 COMMON /zzzzzz/ gbuf(1)
 COMMON /system/ zsys(91)
 COMMON /machin/ mach
 COMMON /output/ pghdg
 COMMON /stapid/ tapid(6)
 COMMON /resdic/ irdict,iropen
!WKBI
!      INCLUDE 'NASNAMES.COM'
 EQUIVALENCE     (zsys( 1),bufsz ),(zsys( 2),otpe  ),  &
     (zsys( 9),nlpp  ),(zsys(11),npages), (zsys(12),nlines),(zsys(24),icfiat),  &
     (zsys(26),cppgct),(zsys(40),nbpw  )
 EQUIVALENCE     (dcparm(1),nrlfl),(dcparm(2),seqno),  &
     (gbuf(1),iobuf(1),ptdic(1))
 DATA    nptp  / 4HNPTP/
 DATA    dpt   / 4HPOOL/
 DATA    nblank/ 4H    /
 DATA    noscar/ 4HXOSC/
 DATA    nxchk / 4HXCHK,4H    /
 DATA    nvps  / 4HXVPS,4H    /
 DATA    nxptdc/ 4HXPTD,4HIC  /, dcparm/4H(non,4HE)  /
 DATA    hdg   / 4H    ,4HADDI,4HTION,4HS TO,4H che,4HCKPO,4HINT ,  &
     4HDICT,4HIONA,4HRY  ,22*4H    /
 DATA    blkcnt/ 90*0  /,      limit /90       /
 DATA    lpurge/ 100   /
 
!     INITIALIZE
 
 filcnt = 0
 reelct = 0
 ldc    =-2
 recsz  = 0
 prgpnt =-1
 CALL sswtch (9,diag09)
 IF (mach < 5) CALL xflszd (0,blksiz,0)
 
!     MASKHI - O000000077777
 maskhi = 32767
 
!     DPLFLG - O004000000000
 dplflg = lshift(1,29)
 
!     SETEOR - O004000000000
 seteor = dplflg
 
!     FNASS  - O010000000000
 fnass  = lshift(1,30)
 
!     EOTFLG - O010000000000
 eotflg = fnass
 
!     ALLON  - O777777777777
 allon  = complf(0)
 
!     NOSGN  - O377777777777
 nosgn  = rshift(allon,1)
 
!     EQFLG  - O400000000000
 eqflg  = complf(nosgn)
 
!     NOFLGS - O003777777777
 noflgs = rshift(allon,nbpw-29)
 
 
!     FIND OSCAR FILE NUMBER IN DPL
 
 j1 = dpl(3)*3 + 1
 DO  j = 4,j1,3
   IF (dpl(j) == noscar) EXIT
 END DO
 20 oscfn = andf(dpl(j+2),maskhi)
 dpfct = oscfn
 
!     ALLOCATE CORE FOR GINO BUFFERS
 
 nptpnt = korsz(gbuf) - bufsz - 1
 dppnt  = nptpnt - bufsz
 fpnt   = dppnt  - bufsz
 IF (fpnt < 1) CALL mesage (-8,0,nxchk)
 
!     INITIALIZE PTDIC PARAMETERS AND LOAD CHECKPOINT DICTIONARY
 
 ngino = nxptdc(1)
 CALL OPEN (*905,nxptdc,gbuf(nptpnt),0)
 CALL READ (*970,*30,nxptdc,dcparm,2,1,recsz)
 30 IF (dcparm(1) /= nxptdc(1)) GO TO 970
 CALL READ (*970,*35,nxptdc,dcparm,2,1,recsz)
 35 ptdtop = 1
 lptdic = nptpnt - ptdtop
 CALL READ (*970,*40,nxptdc,ptdic(ptdtop),lptdic,1,recsz)
 GO TO 940
 40 ptdbot = recsz  + ptdtop - 3
 iopnt  = ptdbot + 6
 liobuf = fpnt - iopnt
 IF (liobuf < 1) CALL mesage (-8,0,nxchk)
 CALL CLOSE (nxptdc,1)
 lcptp = ptdbot + 3
 
!     SAVE CHECKPOINT DMAP SEQ. NO. AND RECORD NO.
 
 ptdic(lcptp  ) = nblank
 ptdic(lcptp+1) = nblank
 ptdic(lcptp+2) = orf(oscar(2),lshift(andf(maskhi,oscar(6))+1,16))
 nptfn = nrlfl
 
!     GET FIRST/NEXT FILE NAME FROM OSCAR ENTRY
 
 i1 = oscar(7)*2 + 6
 loop200:  DO  i = 8,i1,2
   
!     SEE IF FILE IS ALREADY IN DICT
   
   IF (oscar(i) == nvps(1) .AND. oscar(i+1) == nvps(2)) CYCLE loop200
   IF (ldc < 0) GO TO 110
   DO  j = 1,ldc,3
     IF (dict(j) == oscar(i) .AND. dict(j+1) == oscar(i+1)) CYCLE loop200
   END DO
   
!     CHECK FIAT TABLE FOR FILE NAME
   
   110 j1 = fiat(3)*icfiat - 2
   DO  j = 4,j1,icfiat
     IF (oscar(i) == fiat(j+1) .AND. oscar(i+1) == fiat(j+2)) GO TO 120
   END DO
   GO TO 160
   
!     FILE IS IN FIAT - ENTER FILE AND ALL EQUIVALENCED FILES IN DICT
   
   120 IF (andf(fiat(j),maskhi) == maskhi) GO TO 155
   
!     FILE NOT PURGED - CHECK FIAT TRAILER WORDS TO INSURE THAT FILE HAS
!     BEEN GENERATED
   
   IF (fiat(j+3) /= 0 .OR. fiat(j+4) /= 0 .OR. fiat(j+5) /= 0) GO TO 125
   IF (icfiat == 11 .AND. (fiat(j+8) /= 0 .OR. fiat(j+9) /= 0 .OR.  &
       fiat(j+10) /= 0)) GO TO 125
   GO TO 155
   125 IF (fiat(j) < 0) GO TO 145
   ldc = ldc + 3
   dict(ldc  ) = fiat(j+1)
   dict(ldc+1) = fiat(j+2)
   dict(ldc+2) = orf(lshift(j,16),andf(fiat(j),maskhi))
   
!     DESTROY ANY EQUIVS TO THIS FILE
   
!     FIND LAST DICTIONARY REFERENCE TO THIS DATA BLOCK NAME
   
   DO  j = ptdtop,ptdbot,3
     k = ptdbot - (j-ptdtop)
     IF (dict(ldc) == ptdic(k) .AND. ptdic(k+1) == dict(ldc+1)) GO TO 132
   END DO
   GO TO 140
   
!     FILE EXISTS IN DICTIONARY SEE IF IT IS EQUIVED
   
   132 CONTINUE
   IF (andf(ptdic(k+2),eqflg) == 0) GO TO 140
   
!     FILE IS EQUIVED.  PURGE ALL SUBSEQUENT ENTRIES FOR THIS FILE
   
   IF (k == ptdbot) GO TO 140
   DO  j = k,ptdbot,3
     IF (ptdic(j+2) /= ptdic(k+2)) CYCLE
     
!     PURGE FILE
     
     prgpnt = prgpnt + 2
     IF (lpurge < prgpnt+1) GO TO 960
     purge(prgpnt  ) = ptdic(j  )
     purge(prgpnt+1) = ptdic(j+1)
   END DO
   140 CONTINUE
   CYCLE loop200
   145 k = andf(fiat(j),orf(maskhi,eqflg))
   DO  j = 4,j1,icfiat
     IF (andf(fiat(j),orf(maskhi,eqflg)) /= k) CYCLE
     ldc = ldc + 3
     
!     EQUIVALENCED FILE FOUND
     
     dict(ldc  ) = fiat(j+1)
     dict(ldc+1) = fiat(j+2)
     
!     ENTER EQUIVALENCE FLAG, FIAT POINTER AND UCB POINTER IN DICT
     
     dict(ldc+2) = orf(lshift(j,16),k)
   END DO
   CYCLE loop200
   
!     ENTER PURGED FILE IN PURGE TABLE
   
   155 prgpnt = prgpnt + 2
   IF (lpurge < prgpnt+1) GO TO 960
   purge(prgpnt  ) = oscar(i  )
   purge(prgpnt+1) = oscar(i+1)
   CYCLE loop200
   
!     SEE IF FILE IS IN DPL
   
   160 j1 = dpl(3)*3 + 1
   DO  j = 4,j1,3
     IF (oscar(i) == dpl(j) .AND. oscar(i+1) == dpl(j+1)) GO TO 180
   END DO
   GO TO 155
   
!     FILE IS IN DPL - ENTER FILE AND ALL EQUIVALENCED FILES IN DICT
   
   180 k = andf(dpl(j+2),maskhi)
   dpfct = MIN0(oscfn,k)
   DO  j = 4,j1,3
     IF (andf(dpl(j+2),maskhi) /= k) CYCLE
     ldc = ldc + 3
     
!     EQUIVALENCED FILE FOUND
     
     dict(ldc  ) = dpl(j  )
     dict(ldc+1) = dpl(j+1)
     
!     ENTER EQUIVALENCE FLAG, DPLFLG AND FILE NO. IN DICT
     
     dict(ldc+2) = orf(dplflg,andf(dpl(j+2),orf(maskhi,eqflg)))
   END DO
 END DO loop200
 
!     MOVE DICT ENTRIES TO FDICT TABLE
!     GET FIRST NEXT/ENTRY IN DICT
 
 IF (ldc < 1) GO TO 400
 DO  i = 1,ldc,3
   
!     IF DICT ENTRY IS EQUIVALENCED - SEE IF IT IS IN PTDIC
   
   IF (andf(dict(i+2),fnass) == fnass) CYCLE
   IF (dict(i+2) > 0) GO TO 225
   
!     SEARCH BACKWARD FOR PREVIOUS ENTRY
   
   DO  j = ptdtop,ptdbot,3
     k = ptdbot - (j-ptdtop)
     IF (ptdic(k) == dict(i) .AND. ptdic(k+1) == dict(i+1) .AND.  &
         ptdic(k+2) /= 0) GO TO 215
   END DO
   GO TO 225
   
!     DICT ENTRY IS IN PTDIC
   
   215 fdict(i  ) = dict(i  )
   fdict(i+1) = dict(i+1)
   fdict(i+2) = orf(ptdic(k+2),eqflg)
   ucbpnt     = dict(i+2)
   dict(i+2)  = fnass
   
!     ENTER PTDIC FILE NUMBER IN FDICT ENTRIES THAT ARE EQUIVALENCED TO
!     PTDIC ENTRY
   
   ucbpnt = andf(ucbpnt,orf(maskhi,dplflg))
   DO  j = 1,ldc,3
     IF (andf(dict(j+2),orf(maskhi,dplflg)) /= ucbpnt) CYCLE
     fdict(j  ) = dict(j  )
     fdict(j+1) = dict(j+1)
     fdict(j+2) = orf(eqflg,ptdic(k+2))
     dict(j+2)  = fnass
   END DO
   
!     MOVE DICT ENTRY TO FDICT IF NOT ALREADY MOVED
   
   225 IF (andf(dict(i+2),fnass) == fnass) CYCLE
   fdict(i  ) = dict(i  )
   fdict(i+1) = dict(i+1)
   fdict(i+2) = dict(i+2)
   IF (andf(dict(i+2),dplflg) == dplflg) CYCLE
   
!     DICT ENTRY IS FIAT FILE - ENTER NPTP FILE NO. IN FDICT
   
   fdict(i+2) = orf(andf(fdict(i+2),eqflg),nptfn)
   dict (i+2) = orf(dict(i+2),fnass)
   IF (dict(i+2) > 0) GO TO 295
   
!     FILE IS EQUIVALENCED - ENTER NPTP FILE NO. IN FDICT FOR FILES THAT
!     THIS ENTRY IS EQUIVALENCED TO.
   
   ucbpnt = andf(dict(i+2),maskhi)
   j1 = i + 3
   IF (j1 > ldc) GO TO 295
   DO  j = j1,ldc,3
     IF (andf(dict(j+2),maskhi) /= ucbpnt) CYCLE
     fdict(j  ) = dict(j  )
     fdict(j+1) = dict(j+1)
     fdict(j+2) = fdict(i+2)
     dict(j+2)  = orf(dict(j+2),fnass)
   END DO
   295 nptfn = 1 + nptfn
 END DO
 
!     NOW ASSIGN NPTP FILE NUMBERS TO DATA POOL FILES IN SAME ORDER THAT
!     FILES APPEAR ON DATA POOL TAPE.
 
 310 minfn = rshift(allon,1)
 
!     GET FIRST/NEXT DICT ENTRY
 
 DO  i = 1,ldc,3
   IF (andf(dict(i+2),fnass) == fnass) CYCLE
   minfn = MIN0(minfn,andf(dict(i+2),maskhi))
 END DO
 IF (minfn == rshift(allon,1)) GO TO 400
 DO  i = 1,ldc,3
   IF (andf(dict(i+2),fnass)  == fnass) CYCLE
   IF (andf(dict(i+2),maskhi) /= minfn) CYCLE
   fdict(i+2) = orf(nptfn,andf(fdict(i+2),eqflg))
   dict (i+2) = orf(dict(i+2),fnass)
 END DO
 nptfn = nptfn + 1
 GO TO 310
 
!     OPEN DATA POOL TAPE SO IT IS POSITIONED BEFORE FIRST FILE TO
!     CHECKPOINT.
 
 400 IF (dpfct < oscfn) GO TO 401
 j = 2
 dpfct = oscfn
 GO TO 402
 401 j = 0
 dpfct = 1
 402 NAME  = dpt
 CALL OPEN (*905,dpt,gbuf(dppnt),j)
 NAME  = nptp
 
!     OPEN NEW PROBELM NPTP TAPE FOR WRITE
 
 CALL OPEN (*905,nptp,gbuf(nptpnt),3)
 
!     MAKE TEMPORARY ENTRY IN FIST FOR FIAT FILES
 
 ifstmp   = 2*ipfst + 3
 svfst(1) = fist(ifstmp  )
 svfst(2) = fist(ifstmp+1)
 fist(2)  = ipfst + 1
 fist(ifstmp) = 301
 
!     WRITE FILES ON NEW PROBLEM NPTP TAPE AS SPECIFIED IN FDICT.
 
 n1 = nptfn - 1
 
!     GET FIRST/NEXT FDICT ENTRY
 
 n = nrlfl
 IF (ldc < 1 .OR. n1 < n) GO TO 615
 405 DO  i = 1,ldc,3
   IF (andf(fdict(i+2),noflgs) == n) GO TO 415
 END DO
 
!     FDICT ENTRIES SHOULD ALL BE COPIED - MAKE SURE ALL IS O.K.
 
 DO  i = 1,ldc,3
   IF (andf(fdict(i+2),noflgs) > n) GO TO 920
 END DO
 nptfn = n
 GO TO 615
 
!     THIS FDICT ENTRY IS NEXT TO GO ON NEW PROBLEM NPTP TAPE.
 
 415 IF (andf(dict(i+2),dplflg) == dplflg) GO TO 450
 
!     FILE IS IN FIAT TABLE
 
 k = rshift(andf(noflgs,dict(i+2)),16)
 IF (dict(i+2) > 0) GO TO 418
 
!     GET SMALLEST FIAT POINTER FOR EQUIVALENCED FIAT FILES
 
 DO  ii = 1,ldc,3
   IF (andf(dict(ii+2),dplflg) == dplflg) CYCLE
   IF (andf(dict(i+2),maskhi)  == andf(dict(ii+2),maskhi))  &
       k = MIN0(rshift(andf(noflgs,dict(ii+2)),16),k)
 END DO
 
!     INSERT FIAT POINTER IN TEMPORARY FIST ENTRY
 
 418 fist(ifstmp+1) = k - 1
 
!     READ FIRST 2 WORDS OF DATA BLOCK, CHECK NAME AND WRITE TO NEW
!     PROBLEM NPTP TAPE SPECIAL HEADER AND 3 OR 6 TRAILER WORDS
!     (TOTAL OF 5 OR 8 WORDS IN THIS NPTP RECORD)
 
 ngino = fist(ifstmp)
 CALL OPEN (*900,ngino,gbuf(fpnt),0)
 filcnt = filcnt + 1
 IF (filcnt > limit) GO TO 990
 CALL xflszd (-1,blkcnt(filcnt),ngino)
 CALL READ (*930,*930,ngino,head,2,0,recsz)
 DO  j = i,ldc,3
   IF (head(1) == fdict(j) .AND. head(2) == fdict(j+1) .AND.  &
       fdict(j+2) == fdict(i+2)) GO TO 445
 END DO
 GO TO 930
 445 CALL WRITE (nptp,head,2,0)
 IF (icfiat == 11) GO TO 447
 CALL WRITE (nptp,fiat(k+3),3,1)
 GO TO 448
 447 CALL WRITE (nptp,fiat(k+3),3,0)
 CALL WRITE (nptp,fiat(k+8),3,1)
 
!     COPY ENTIRE FILE ONTO NEW PROBLEM NPTP TAPE USING CPYFIL
 
 448 CALL WRITE  (nptp,head,2,0)
 CALL cpyfil (ngino,nptp,iobuf(iopnt),liobuf,recsz)
 CALL CLOSE  (ngino,1)
 GO TO 600
 
!     FILE IS ON POOL -- POSITION POOL AND COPY FILE USING CPYFIL
 
 450 ngino = dpt
 k = andf(dict(i+2),maskhi)
 CALL skpfil (dpt,k-dpfct)
 dpfct  = k + 1
 filcnt = filcnt + 1
 IF (filcnt > limit) GO TO 990
 CALL xflszd (k,blkcnt(filcnt),0)
 CALL cpyfil (dpt,nptp,iobuf(iopnt),liobuf,recsz)
 
!     GET NEXT FDICT ENTRY
 
 600 CALL eof (nptp)
 n = n + 1
 IF (n <= n1) GO TO 405
 
!     RESTORE FIST ENTRY
 
 fist(ifstmp  ) = svfst(1)
 fist(ifstmp+1) = svfst(2)
 
!     WRITE VPS TABLE ONTO NEW PROBLEM NPTP TAPE
!     MAKE ENTRY IN FDICT FOR VPS TABLE
 
 615 ldc = ldc + 3
 fdict(ldc  ) = nvps(1)
 fdict(ldc+1) = nblank
 fdict(ldc+2) = nptfn
 eorflg = seteor
 i = ldc
 CALL WRITE (nptp,nvps,5,1)
 CALL WRITE (nptp,vps,vps(2),1)
 
!     WRITE CEITBL TABLE ONTO PROBLEM TAPE
 
 CALL WRITE (nptp,ceitbl,ceitbl(2),1)
 
!     WRITE /SYSTEM/ ONTO PROBLEM TAPE
 
 CALL WRITE (nptp,bufsz,20,1)
 CALL eof   (nptp)
 CALL CLOSE (nptp,2)
 
!     POSITION DATA POOL TAPE AT CORRECT OSCAR ENTRY FOR RETURN TO XSEM
 
 IF (dpfct == oscfn) GO TO 675
 CALL REWIND (dpt)
 IF (oscfn > 1) CALL skpfil (dpt,oscfn-1)
 j1 = oscar(2)
 DO  j = 1,j1
   CALL fwdrec (*910,dpt)
 END DO
 675 CALL CLOSE (dpt,2)
 
!     UPDATE PTDIC AND ASSOCIATED VARIABLES
 
 nrlfl  = nptfn + 1
 ptdbot = lcptp
 loop690:  DO  i = 1,ldc,3
   DO  j = ptdtop,ptdbot,3
     
!     SCAN PTDIC TO SEE IF FILE IS ALREADY THERE
     
     IF (fdict(i) == ptdic(j) .AND. fdict(i+1) == ptdic(j+1) .AND.  &
         fdict(i+2) == ptdic(j+2)) CYCLE loop690
   END DO
   
!     ENTER FILE IN PTDIC
   
   ptdbot = ptdbot + 3
   ptdic(ptdbot  ) = fdict(i  )
   ptdic(ptdbot+1) = fdict(i+1)
   ptdic(ptdbot+2) = fdict(i+2)
 END DO loop690
 
!     PUT PURGED FILES IN PTDIC
 
 IF (prgpnt < 1) GO TO 800
 loop710:  DO  i = 1,prgpnt,2
   DO  j = ptdtop,ptdbot,3
     IF (purge(i) == ptdic(j) .AND. purge(i+1) == ptdic(j+1) .AND.  &
         ptdic(j+2) == 0) CYCLE loop710
   END DO
   ptdbot = ptdbot + 3
   ptdic(ptdbot  ) = purge(i)
   ptdic(ptdbot+1) = purge(i+1)
   ptdic(ptdbot+2) = 0
 END DO loop710
 
!     CHECK FOR PTDIC OVERFLOW
 
 800 IF (ptdbot+3-ptdtop > lptdic) GO TO 940
 
 
!     PUNCH AND PRINT LATEST ENTRIES IN PTDIC
!     INITIALIZE PAGE HEADING AND CHECK PAGE COUNT
 
 IF (diag09 == 1) GO TO 802
 DO  i = 1,32
   pghdg(i+ 96) = hdg(i)
   pghdg(i+128) = nblank
   pghdg(i+160) = nblank
 END DO
 IF (cppgct /= npages) CALL page
 802 CONTINUE
 i1 = ((lcptp  - ptdtop)/3) + 1
 i2 = ((ptdbot - ptdtop)/3) + 1
 DO  i = i1,i2
   j1 = (i-1)*3 + ptdtop
   j2 = j1 + 2
   
!     SEPARATE FLAGS, REEL NO., FILE NO.
   
   nflags = 0
   IF (ptdic(j2) < 0) nflags = 4
   nflags = orf(nflags,rshift(andf(ptdic(j2),nosgn),29))
   nreel  = rshift(andf(ptdic(j2),noflgs),16)
   nfile  = andf(ptdic(j2),maskhi)
   seqno  = 1 + seqno
   IF (ptdic(j1) == nblank) GO TO 805
!      IF (IROPEN .EQ. 1) GO TO 815
!      OPEN (UNIT=4, FILE=DIC, STATUS='UNKNOWN')
!      IROPEN = 1
   815 WRITE  (irdict,820) seqno,ptdic(j1),ptdic(j1+1),nflags,nreel,nfile
   820 FORMAT (i10,4H,   ,2A4,12H,   flags = ,i1,11H,   reel = ,i2,  &
       11H,   FILE = ,i6)
   IF (diag09 == 1) CYCLE
   nlines = nlines + 1
   IF (mach < 5 .AND. nfile /= 0 .AND. ptdic(j1) /= nvps(1))  &
       nlines = nlines + 1
   IF (nlines >= nlpp) CALL page
   WRITE  (otpe,821) seqno,ptdic(j1),ptdic(j1+1),nflags,nreel,nfile
   821 FORMAT (1H ,i9 ,4H,   ,2A4,12H,   flags = ,i1,11H,   reel = ,i2,  &
       11H,   FILE = ,i6)
   IF (mach < 5 .AND. nfile /= 0 .AND. ptdic(j1) /= nvps(1))  &
       WRITE (otpe,822) ptdic(j1),ptdic(j1+1),blkcnt(i-i1),blksiz
   822 FORMAT (13X,6H FILE ,2A4, 9H contains, i10,  &
       28H blocks, each BLOCK contains,i5,7H words.)
   CYCLE
   805 CONTINUE
!      IF (IROPEN .EQ. 1) GO TO 8055
!      OPEN (UNIT=4, FILE=DIC, STATUS='UNKNOWN')
!      IROPEN = 1
!8055  CONTINUE
   WRITE  (irdict,806) seqno,nreel
   806 FORMAT (i10,36H,   reenter at dmap sequence NUMBER ,i5)
   IF (diag09 == 1) CYCLE
   nlines = nlines + 2
   IF (nlines >= nlpp) CALL page
   WRITE  (otpe,807) seqno,nreel
   807 FORMAT (1H ,/1H ,i9,36H,   reenter at dmap sequence NUMBER ,i5)
 END DO
 
!     WRITE PTDIC ONTO XPTD
 
 ngino = nxptdc(1)
 CALL OPEN  (*905,nxptdc,gbuf(nptpnt),1)
 CALL WRITE (nxptdc,nxptdc,2,1)
 CALL WRITE (nxptdc,dcparm,2,1)
 CALL WRITE (nxptdc,ptdic(ptdtop),ptdbot+3-ptdtop,1)
 CALL CLOSE (nxptdc,1)
 cppgct = npages
 
 fist(2) = ipfst
 RETURN
 
 
!     ERRORS -
 
 900 n = 1101
 ASSIGN 901 TO RETURN
 GO TO 980
 901 WRITE  (otpe,902) fdict(i),fdict(i+1)
 902 FORMAT (4X,26HCOULD NOT OPEN FILE NAMED ,2A4)
 GO TO 995
 
 905 n = 1102
 ASSIGN 906 TO RETURN
 GO TO 985
 906 WRITE (otpe,902) ngino,nblank
 GO TO 995
 
 910 n = 1103
 ASSIGN 911 TO RETURN
 GO TO 985
 911 WRITE  (otpe,912)
 912 FORMAT (4X,43HUNABLE TO POSITION DATA pool tape correctly )
 GO TO 995
 
 920 n = 1104
 ASSIGN 921 TO RETURN
 GO TO 985
 921 WRITE  (otpe,922)
 922 FORMAT (4X,24HFDICT table is incorrect )
 GO TO 995
 
 930 n = 1105
 ASSIGN 931 TO RETURN
 GO TO 980
 931 WRITE  (otpe,932) fdict(i),fdict(i+1),head(1),head(2)
 932 FORMAT (4X,29HCANNOT find DATA BLOCK NAMED ,2A4,17H header record = ,2A4)
 GO TO 995
 
 940 n = 1106
 ASSIGN 941 TO RETURN
 GO TO 985
 941 WRITE  (otpe,942)
 942 FORMAT (4X,32HCHECKPOINT dictionary overflowed)
 GO TO 995
 
 960 n = 1108
 ASSIGN 961 TO RETURN
 GO TO 985
 962 FORMAT (4X,22HPURGE table overflowed)
 961 WRITE  (otpe,962)
 GO TO 995
 
 970 n = 1109
 ASSIGN 971 TO RETURN
 GO TO 985
 971 WRITE (otpe,932) nxptdc,dcparm
 GO TO 995
 
!     USER FATAL ERROR
 
 980 WRITE  (otpe,981) ufm,n
 981 FORMAT (a23,i5)
 GO TO 987
 
!     SYSTEM FATAL ERROR
 
 985 CALL page2 (3)
 WRITE  (otpe,986) sfm,n
 986 FORMAT (a25,i5)
 987 GO TO RETURN, (901,906,911,921,931,941,961,971)
 
 990 WRITE  (otpe,991) sfm
 991 FORMAT (a25,', BLKCNT ARRAY EXCEEDED IN XCHK')
 
 995 CALL mesage (-37,0,nxchk)
 RETURN
END SUBROUTINE xchk
