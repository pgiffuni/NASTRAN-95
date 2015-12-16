SUBROUTINE xgpibs
     
!     PURPOSE OF THIS ROUTINE IS TO INITIALIZE MACHINE DEPENDENT
!     CONSTANTS FOR XGPI AND ASSOCIATED ROUTINES AND TO INITIALIZE
!     THE MODULE LINK TABLE.
 
 IMPLICIT INTEGER (a-z)
 EXTERNAL        lshift,rshift,andf,orf,complf
 DIMENSION       opbuff(1),pghdg(113),hdg1(32),hdg2(32),lnkedt(15),  &
     enddta(2),opncor(1),utilty(1),nwptyp(6),ll(15),  &
     NONE(2),lnkspc(1),inbuff(20),modnam(2),os(2)
 CHARACTER (LEN=27) :: swm
 CHARACTER (LEN=25) :: sfm
 CHARACTER (LEN=29) :: uim
 CHARACTER (LEN=25) :: uwm
 CHARACTER (LEN=23) :: ufm
 COMMON /xmssg / ufm,uwm,uim,sfm,swm
 COMMON /machin/ ijhalf(4),mchnam
 COMMON /system/ xsys(90),lpch
 COMMON /zzzzzz/ core(1)
 COMMON /xgpi2 / lmpl,mplpnt,mpl(1)
 COMMON /xgpi2x/ ixx(1)
 COMMON /xlink / lxlink,maxlnk,mlink(1)
 COMMON /xlkspc/ llink,link(1)
 COMMON /output/ pghdg
 COMMON /lhpwx / lhpw(2),nwpic
 COMMON /xgpic / a(22),ncpw,nbpc,nwpc,maskhi,masklo,isgnon,nosgn,  &
     iallon,masks(1)
 COMMON /xgpid / b(7),itape,iappnd,intgr,losgn,noflgs,seteor,  &
     eotflg,ieqflg,cpntry(7),jmp(7)
 COMMON /xgpi6 / d(5),iplus
 EQUIVALENCE     (xsys( 2) ,optap )  ,(xsys(9) ,nlpp ) ,  &
     (xsys(12) ,nlines)  ,(xsys(4) ,intap) ,  &
     (opncor(1),lnkspc(1),opbuff(1),os(2)  ,utilty(1))  &
     ,               (core(1)  ,os(1)    ,lopncr)
 DATA    modnam/ 4HCHKP,4HNT  /
 DATA    DELETE/ 4HDELE       /,  xnone /4HNONE/
 DATA    enddta/ 4HENDD,4HATA /,  dolsgn/4H$   /
 
!             NWPTYP = NUMBER OF WORDS PER PARAMETER TYPE CODE
!                        INT,   REAL,   BCD,   D.P.,   CMPLX,  D.P.CMPLX
 DATA    nwptyp/      1,      1,     2,      2,       2,     4   /
 DATA    nblank/ 4H    /, NONE/4H(non,4HE)  /
 DATA    lnkedt/ 4H  1 ,4H  2 ,4H  3 ,4H  4 ,4H  5 ,4H  6 ,4H  7 ,  &
     4H  8 , 4H  9 ,4H 10 ,4H 11 ,4H 12 ,4H 13 ,4H 14 ,4H 15 /
 DATA    hdg1  / 4HMODU,4HLE -,4H dma,4HP na,4HME -,4H MOD,4HULE ,  &
     4HENTR, 4HY - ,4HLINK,4HS mo,4HDULE,4H res,4HIDES,4H in ,  &
     4HON  , 16*4H        /
 DATA    hdg2  / 4HINDE,4HX   ,4H of ,4HMODU,4HLE  ,4H poi,4HNT n,  &
     4HAME , 24*4H        /
 
!     INITIALIZE MACHINE DEPENDENT CONSTANTS FOR XGPI
!     SEE SUBROUTINE XGPIDD FOR DESCRIPTION OF CONSTANTS.
 
!     INITIALIZE  /XGPIC/
 
!     NCPW   = NUMBER OF CHARACTERS PER WORD
!     NBPC   = NUMBER OF BITS PER CHARACTER
!     NWPC   = NUMBER OF WORDS PER CARD = NWPIC
!                 7094         360        1108             6600
!    MASKLO = 017777600000, 7FFF0000, 017777600000, 00000000017777600000
!    ISGNON = 400000000000, 80000000, 400000000000, 40000000000000000000
!    NOSGN  = 377777777777, 7FFFFFFF, 377777777777, 37777777777777777777
!    IALLON = 777777777777, FFFFFFFF, 777777777777, 77777777777777777777
 
!    MASKHI = MASK FOR LOW ORDER 16 BITS AND SIGN BIT = 32767,
!             INITIALIZED IN XGPIDD
 
 ncpw   = xsys(41)
 nbpc   = xsys(39)
 nwpc   = nwpic
 masklo = lshift(maskhi,16)
 isgnon = lshift(1,xsys(40)-1)
 nosgn  = complf(isgnon)
 iallon = complf(0)
 
!     GENERATE MASKS ARRAY
!     MASK IS IN 4 PARTS - MASK DESCRIPTION WILL BE GIVEN IN TERMS OF
!                          IBM 360
!     PART 1 - FFOOOOOO,OOFFOOOO,OOOOFFOO,OOOOOOFF
!     PART 2 - COMPLEMENT OF PART 1
!     PART 3 - FFFFFFFF,OOFFFFFF,OOOOFFFF,OOOOOOFF
!     PART 4 - COMPLEMENT OF PART 3
 
 mhibyt = lshift(iallon,nbpc*(ncpw-1))
 DO  j = 1,ncpw
   masks(j) = rshift(mhibyt,nbpc*(j-1))
   j2 = j + ncpw
   masks(j2) = complf(masks(j))
   j3 = 2*ncpw + j
   masks(j3) = rshift(iallon,nbpc*(j-1))
   j4 = 3*ncpw + j
   masks(j4) = complf(masks(j3))
 END DO
 
!     INITIALIZE  /XGPID/
 
!                 7094         360        1108             6600
!    ITAPE  = 000000100000, 00008000, 000000100000, 00000000000000100000
!    IAPPND = 010000000000, 40000000, 010000000000, 00000000010000000000
!    INTGR  = 400000000001, 80000001, 400000000001, 40000000000000000001
!    LOSGN  = 000000100000, 00008000, 000000100000, 00000000000000100000
!    NOFLGS = 000377777777, 03FFFFFF, 000377777777, 00000000000377777777
!    SETEOR = 004000000000, 20000000, 004000000000, 00000000004000000000
!    EOTFLG = 010000000000, 40000000, 010000000000, 00000000010000000000
!    IEQFLG = 400000000000, 80000000, 400000000000, 40000000000000000000
!    CPNTRY(3) = CHKPNT MODULE INDEX/TYPE CODE
!    NTRY(6)= 400000000001, 80000001, 400000000001, 40000000000000000001
!    JMP(3) = JUMP MODULE INDEX/TYPE CODE
 
 itape  = lshift(1,15)
 iappnd = lshift(1,30)
 intgr  = orf(isgnon,1)
 losgn  = lshift(1,15)
 noflgs = rshift(iallon,xsys(40)-26)
 seteor = lshift(1,29)
 eotflg = lshift(1,30)
 ieqflg = isgnon
 
!     PRINT MPL CONTENTS IF DIAG 31 IS ON
 
 ASSIGN 40 TO irtn
 CALL sswtch (31,l)
 IF (l /= 0) CALL mplprt
 
!     GET CHKPNT MODULE INDEX
 
 20 modidx = 1
 mplpnt = 1
 30 IF (mpl(mplpnt+1) == modnam(1) .AND. mpl(mplpnt+2) == modnam(2))  &
     GO TO irtn, (40,50)
 modidx = modidx + 1
 mplpnt = mplpnt + mpl(mplpnt)
 IF (mplpnt > lmpl .OR. mpl(mplpnt) < 1) GO TO 1240
 GO TO 30
 40 cpntry(3) = lshift(modidx,16) + 4
 
!     GET JUMP MODULE INDEX
 
 ASSIGN 50 TO irtn
 modnam(1) = jmp(4)
 modnam(2) = jmp(5)
 GO TO 20
 50 jmp(3)    = lshift(modidx,16) + 3
 cpntry(6) = orf(isgnon,1)
 jmp(6)    = cpntry(6)
 
!     COMPUTE LENGTH OF OPENCORE (SUBTRACT OFF SOME FOR UTILITY BUFFERS)
 
 lopncr = korsz(opncor) - xsys(1) - 1
 utltop = lopncr + 1
 utlbot = utltop + xsys(1) - 1
 
!     INITIALIZE  /XGPI2/ (I.E. MPL TABLE)
 
!     LOAD FLOATING POINT NUMBERS INTO MPL FROM ARRAY IN /XGPI2X/
 
 mplpnt = 1
 60 IF (mpl(mplpnt) < 4) GO TO 150
 IF (mpl(mplpnt+3) < 1 .OR. mpl(mplpnt+3) > 2) GO TO 150
 
!     MPL ENTRY HAS MODULE TYPE CODE 1 OR 2 - PROCESS PARAMETER SECTION.
 
 i = mplpnt + 7
 
!     CHECK FOR END OF MPL ENTRY
 
 70 IF (i >= mplpnt+mpl(mplpnt)) GO TO 150
 
!     CHECK VALIDITY OF PARAMETER TYPE CODE
 
 j = IABS(mpl(i))
 IF (j < 1 .OR. j > 6) GO TO 1230
 l = 1
 
!     SEE IF PARAMETER VALUE FOLLOWS TYPE CODE.
 
 IF (mpl(i) < 0) GO TO 100
 
!     GET LENGTH OF PARAMETER VALUE TO BE LOADED.
 
 l = nwptyp(j)
 
!     A VALUE FOLLOWS IF TYPE CODE IS INTEGER OR BCD - OTHERWISE AN
!     INDEX INTO A TABLE CONTAINING THE VALUE FOLLOWS THE TYPE CODE.
 
 IF (j == 1 .OR. j == 3) GO TO 90
 
!     GET INDEX INTO VALUE TABLE - NOTE INDEX MUST BE CONVERTED FROM
!     DOUBLE PRECISION INDEX TO ONE DIMENSIONAL INDEX.
 
 m  = mpl(i+1)*2 - 1
 DO  k = 1,l
   n  = k + m - 1
   k1 = i + k
   mpl(k1) = ixx(n)
 END DO
 90 i  = i + 1
 
!     INCREMENT TO NEXT PARAMETER TYPE CODE.
 
 100 i  = i + l
 GO TO 70
 
!     GET NEXT MPL ENTRY
 
 150 IF (mpl(mplpnt)+mplpnt > lmpl) GO TO 160
 mplpnt = mplpnt + mpl(mplpnt)
 IF (mpl(mplpnt) < 1) GO TO 1240
 GO TO 60
 160 CONTINUE
 
!     INITIALIZE /XLINK/
 
!     MAXLNK = MAXIMUM NUMBER OF LINKS THAT CAN BE HANDLED. IF MAXLNK IS
!              INCREASED THEN LNKEDT TABLE MUST BE INCREASED.
!              (MAXLNK WAS SET IN SEMDBD ROUTINE)
 
!     MOVE LINK TABLE INTO OPEN CORE
 
 lnktop = 1
 lnkbot = llink + lnktop - 5
 DO  j = 1,llink
   lnkspc(j) = link(j)
 END DO
 
!     UPDATE LNKSPC TABLE IF SENSE SWITCH 29 IS ON
 
 CALL sswtch (29,l)
 IF (l == 0) GO TO 600
 ASSIGN 280 TO irtn
 
!     PROCESS INPUT CARD (NOTE-DO NOT USE VARIABLES I,J OR M)
 
 210 CALL page1
 nlines = nlines + 2
 WRITE  (optap,220)
 220 FORMAT (42H0LINK specification table update deck echo )
 230 nlines = nlines + 1
 IF (nlines >= nlpp) GO TO 210
 CALL xread (*240,inbuff)
 GO TO 260
 240 CALL page2 (2)
 WRITE  (optap,250) ufm
 250 FORMAT (a23,' 220, MISSING ENDDATA CARD.')
 GO TO 1250
 260 CONTINUE
 WRITE  (optap,270) inbuff
 270 FORMAT (5X,20A4)
 
!     CHECK FOR COMMENT CARD
 
 IF (khrfn1(0,1,inbuff(1),1) == khrfn1(0,1,dolsgn,1)) GO TO 230
 
!     CONVERT CARD IMAGE
 
 CALL xrcard (utilty(utltop),utlbot-utltop+1,inbuff)
 IF (utilty(utltop) == 0) GO TO 230
 
!     CHECK FOR ENDDATA CARD
 
 IF (utilty(utltop+1) == enddta(1) .AND.  &
     utilty(utltop+2) == enddta(2)) GO TO 380
 GO TO irtn, (280,330)
 
!     CHECK FORMAT OF CARD
 
 280 IF (utilty(utltop) < 2) GO TO 1220
 
!     SEE IF MODULE NAME IS IN LNKSPC TABLE
 
 DO  i = lnktop,lnkbot,5
   IF (lnkspc(i  ) == utilty(utltop+1) .AND.  &
       lnkspc(i+1) == utilty(utltop+2)) GO TO 300
 END DO
 
!     MODULE IS NOT IN LNKSPC - MAKE NEW ENTRY
 
 lnkbot = lnkbot + 5
 IF (lnkbot > lopncr) GO TO 1200
 i = lnkbot
 
!     TRANSFER MODULE NAME AND ENTRY POINT TO LNKSPC
 
 300 lnkspc(i  ) = utilty(utltop+1)
 lnkspc(i+1) = utilty(utltop+2)
 lnkspc(i+2) = utilty(utltop+3)
 lnkspc(i+3) = utilty(utltop+4)
 
!     CHECK FOR DELETE OR NONE
 
 IF (utilty(utltop  ) ==      2) GO TO 320
 IF (utilty(utltop+5) == DELETE) GO TO 310
 IF (utilty(utltop+5) /=  xnone) GO TO 1220
 
!     MODULE HAS NO ENTRY POINT
 
 lnkspc(i+2) = NONE(1)
 lnkspc(i+3) = NONE(2)
 m = 0
 j = 7
 IF (utilty(utltop+7) /= -1) j = 9
 GO TO 330
 
!     MODULE IS TO BE DELETED
 
 310 lnkspc(i) = 0
 GO TO 370
 
!     GENERATE A LINK FLAG WORD
 
 320 m = 0
 j = 5
 
!     CHECK MODE WORD
 
 330 k = utltop + j
 IF (utilty(k) < 0.0) THEN
   GO TO   340
 ELSE IF (utilty(k) == 0.0) THEN
   GO TO   350
 ELSE
   GO TO   360
 END IF
 
!     INTEGER FOUND
 
 340 IF (utilty(k) /= -1) GO TO 1220
 m = orf(m,lshift(1,utilty(k+1)-1))
 j = j + 2
 GO TO 330
 
!     CONTINUE MODE FOUND
 
 350 j = 1
 ASSIGN 330 TO irtn
 GO TO 230
 
!     END OF INSTRUCTION FOUND
 
!     TRANSFER GENERATED LINK WORD TO LNKSPC ENTRY
 
 360 IF (utilty(k) /= nosgn) GO TO 1220
 j = i + 4
 lnkspc(j) = m
 
!     PROCESS NEXT INPUT CARD
 
 370 ASSIGN 280 TO irtn
 GO TO 230
 
!     PUNCH OUT LNKSPC TABLE IF SENSE SWITCH 28 IS ON.
 
 380 CALL sswtch (28,l)
 IF (l == 0) GO TO 600
 
!     ELIMINATE DELETED LNKSPC ENTRIES
 
 390 DO  i = lnktop,lnkbot,5
   IF (lnkspc(i) == 0) GO TO 410
 END DO
 GO TO 430
 410 k = i + 4
 n = lnkbot - 1
 DO  m = i,k
   n = n + 1
   lnkspc(m) = lnkspc(n)
 END DO
 lnkbot = lnkbot - 5
 GO TO 390
 430 CALL page2 (2)
 WRITE  (optap,440)
 440 FORMAT (98H0***user requests link specification table be punched &
     &out for use in recompiling SUBROUTINE xlnkdd )
 WRITE  (lpch,450)
 450 FORMAT (70(1H*),/38HLINK spec. table for SUBROUTINE xlnkdd )
 j  = lnkbot - lnktop + 5
 n  = j/90
 WRITE  (lpch,460) j
 460 FORMAT (6X,16HDIMENSION link (,i4,1H))
 k  = 90
 IF (n == 0) GO TO 490
 DO  i = 1,n
   i10= i/10
   i1 = i - 10*i10
   WRITE  (lpch,470) i10,i1,k
   470 FORMAT (5X,2H1,,9X,4HLINK,2I1,1H(,i4,1H))
 END DO
 490 k  = MOD(j,90)
 i  = n + 1
 i10= i/10
 i1 = i - 10*i10
 IF (k > 0) WRITE (lpch,470) i10,i1,k
 WRITE  (lpch,500) j
 500 FORMAT (6X,28HCOMMON/xlkspc/ llink, klink(,i4,1H),/,  &
     6X,34HEQUIVALENCE (link(   1),link01(1)) )
 IF (k > 0) n = n + 1
 IF (n < 2) GO TO 530
 DO  i = 2,n
   i10= i/10
   i1 = i - 10*i10
   k  = 90*(i-1) + 1
   WRITE  (lpch,510) k,i10,i1
   510 FORMAT (5X,2H1,,11X,6H(link( ,i4,6H),link ,2I1,4H(1)))
 END DO
 530 CONTINUE
 j  = lnktop - 1
 m  = 0
 540 j  = j + 1
 m  = m + 1
 m10= m/10
 m1 = m - 10*m10
 k  = MIN0(j+89,lnkbot+4)
 WRITE  (lpch,550) m10,m1,(lnkspc(i),i=j,k)
 550 FORMAT (6X,9HDATA link,2I1,1H/, /,  &
     5X,4H1 4H,a4,3H,4H,a4,4H, 4H,a4,3H,4H,a4,1H,,i6,/,  &
     (5X,4H1,4H,a4,3H,4H,a4,4H, 4H,a4,3H,4H,a4,1H,,i6))
 WRITE  (lpch,560)
 560 FORMAT (5X,2H1/)
 j = k
 IF (j < lnkbot+4) GO TO 540
 j = lnkbot - lnktop + 5
 WRITE  (lpch,570) j
 570 FORMAT (6X,8HLLINK = ,i4)
 
!     INITIALIZE PAGE HEADING
 
 600 DO  i = 1,32
   pghdg(i+ 96) = hdg1(i)
   pghdg(i+128) = hdg2(i)
   pghdg(i+160) = nblank
 END DO
 pghdg(  113) = mchnam
 nlines = nlpp
 
!     INITIALIZE O/P BUFFER PARAMETERS - O/P BUFFERS ARE IN OPEN CORE
 
 opbtop = lnkbot + 5
 nxtlin = opbtop - 20
 
!     GET FIRST/NEXT MPL ENTRY
 
 mplpnt = 1
 modidx = 1
 
!     CHECK FOR DECLARATIVE OR NULL ENTRY
 
 620 IF (mpl(mplpnt+3) > 4 .OR. mpl(mplpnt+3) < 1) GO TO 630
 GO TO 700
 630 IF (mpl(mplpnt) < 1) GO TO 800
 mplpnt = mplpnt + mpl(mplpnt)
 modidx = modidx + 1
 IF (mplpnt < lmpl) GO TO 620
 GO TO 800
 
!     PREPARE TO GENERATE NEXT LINE OF OUTPUT
 
 700 nxtlin = nxtlin + 20
 i  = nxtlin + 19
 IF (i > lopncr) GO TO 1240
 DO  j = nxtlin,i
   opbuff(j) = nblank
 END DO
 
!     MODULE INDEX INTO WORD 1 OF O/P ENTRY
 
 opbuff(nxtlin) = modidx
 
!     DMAP NAME TO WORDS 2,3 OF O/P ENTRY
 
 opbuff(nxtlin+1) = mpl(mplpnt+1)
 opbuff(nxtlin+2) = mpl(mplpnt+2)
 
!     GET ENTRY POINT NAME AND ENTER IN WORDS 4,5 OF O/P ENTRY
 
 opbuff(nxtlin+3) = NONE(1)
 opbuff(nxtlin+4) = NONE(2)
 DO  i = lnktop,lnkbot,5
   IF (lnkspc(i) == mpl(mplpnt+1) .AND. lnkspc(i+1) == mpl(mplpnt+2))  &
       GO TO 730
 END DO
 GO TO 630
 730 opbuff(nxtlin+3) = lnkspc(i+2)
 opbuff(nxtlin+4) = lnkspc(i+3)
 
!     EXAMINE LINK FLAG
 
 l = lnkspc(i+4)
 DO  j = 1,maxlnk
   IF (andf(l,lshift(1,j-1)) == 0) CYCLE
   
!     MODULE IS IN LINK J - SET BIT J IN MAIN LINK TABLE AND O/P BUFFER
!     MAKE SURE LINK TABLE IS LONG ENOUGH.
   
   IF (lxlink < modidx) GO TO 1210
   mlink(modidx) = orf(mlink(modidx),lshift(1,j-1))
   k = nxtlin + j + 4
   opbuff(k) = lnkedt(j)
 END DO
 GO TO 630
 
!     SEE IF O/P BUFFER IS TO BE PRINTED (I.E. SENSE SWITCH 31 IS ON)
 
 800 CALL sswtch (31,l)
 IF (l /= 0) GO TO 810
 
!     PRINT O/P BUFFER IF LINK DRIVER PUNCHED O/P  REQUESTED(I.E. SENSE
!     SWITCH 30 IS ON)
 
 CALL sswtch (30,l)
 IF (l == 0) RETURN
 810 DO  i = opbtop,nxtlin,20
   nlines = nlines + 1
   IF (nlines >= nlpp) CALL page
   j = i + 19
   WRITE  (optap,820) (opbuff(k),k=i,j)
   820 FORMAT (5X,i6,3X,2A4,4X,2A4,7X,15A4)
 END DO
 
!     SEE IF ANY DRIVERS SHOULD BE PUNCHED  (I.E. SENSE SWITCH 30 ON)
 
 CALL sswtch (30,l)
 IF (l == 0) RETURN
 CALL page1
 nlines = nlines + 2
 WRITE  (optap,910)
 910 FORMAT ('0USER REQUESTS PUNCHED OUTPUT FOR THE FOLLOWING LINK ',  &
     'DRIVER SUBROUTINES')
 WRITE  (lpch,920)
 920 FORMAT (70(1H*), /,' INSERT FOLLOWING FORTRAN CODE IN RESPECTIVE',  &
     ' LINK DRIVER ROUTINES')
 DO  j = 1,maxlnk
   CALL sswtch (j,l)
   IF (l == 0) CYCLE
   j10= j/10
   j1 = j - 10*j10
   WRITE  (lpch,930) j10,j1,j
   930 FORMAT (70(1H*), /6X,15HSUBROUTINE xsem,2I1, /6X,12HDATA thislk,  &
       /,i2,1H/)
   WRITE  (lpch,940) j10,j1
   940 FORMAT (6X,21HDATA subnam/4HXSEM,4H,2I1,3H  /)
   nlines = nlines + 2
   IF (nlines >= nlpp) CALL page
   WRITE  (optap,950) j10,j1
   950 FORMAT (9H0    xsem,2I1)
   
!     USER REQUESTS PUNCHED O/P FOR LINK J
!     SEARCH LINK TABLE FOR MODULES RESIDING IN LINK J
   
   frstin = 0
   l      = 0
   lastin = 0
   nxtgrp = 1000
   DO   i = 1,lxlink
     ll(l+1) = 940
     IF (andf(mlink(i),lshift(1,j-1)) /= 0) ll(l+1) = 2000 + i
     IF (i - lastin < 0) THEN
       GO TO   980
     ELSE IF (i - lastin == 0) THEN
       GO TO   990
     END IF
     960 IF (frstin <= 0 .AND. ll(l+1) == 940) CYCLE
     frstin = i
     lastin = MIN0(i+180,lxlink)
     lstgrp = nxtgrp
     nxtgrp = nxtgrp - 5
     970 FORMAT (i5,8H GO TO (  )
     980 l = l + 1
     IF (ll(l) /= 940) GO TO 1000
     
!     ONLY TWO CONSECUTIVE BRANCHES TO 940 IN COMPUTED  -GO TO -
     
     IF (last+2 >= i) GO TO 1010
     lastin = last
     l      = MAX0(0,l-1+lastin-i)
     990 ll(15) = ll(l+1)
     IF (l > 0) THEN
       GO TO  1020
     ELSE
       GO TO  1050
     END IF
     1000 last = i
     1010 IF (l < 13) CYCLE
     1020 IF (frstin == ll(1)-2000) WRITE (lpch,970) nxtgrp
     lk = MIN0(l,10)
     WRITE  (lpch,1030) (ll(k),k=1,lk)
     1030 FORMAT (5X,1H1,10(i5,1H,))
     l  = l - lk
     DO  k = 1,l
       ll(k) = ll(k+10)
     END DO
     IF (i < lastin) CYCLE
     1050 last = nxtgrp + 15
     IF (i == lxlink) last = 970
     IF (frstin == lastin) GO TO 1070
     frstin = frstin - 1
     WRITE  (lpch,1060) ll(15),lstgrp,lastin,last,frstin,nxtgrp
     1060 FORMAT (5X, 1H1, i5, 4H ),i ,/,  &
     i5, 14H IF (modx  >  ,i3, 8H) GO TO ,i5, /, &
     6X, 10HI = modx - , i3, /,  &
     6X, 17HIF (i ) 940, 940, ,i5)
     GO TO 1090
     1070 WRITE  (lpch,1080) lstgrp,frstin,ll(15),last
     1080 FORMAT (i5,12H IF (modx - ,i3,7H ) 940,,i5,1H,,i5)
     1090 nxtgrp = last
     frstin = -1
   END DO
   
!     PUNCH OUT GO TO AND IF STATEMENTS FOR LAST GROUP OF MODULES IN
!     LINK J.
   
   IF (frstin /= 0) GO TO 1120
   
!     CANNOT FIND ANY MODULES IN THIS LINK
   
   nlines = nlines + 2
   WRITE  (optap,1110) j
   1110 FORMAT (1H0,10X,29HTHERE are no modules in link  ,i3)
   CYCLE
   1120 IF (last /= 970) WRITE (lpch,1130) nxtgrp
   1130 FORMAT (i5,34H IF (modx - lxlink ) 940, 940, 970 )
   
!     SEARCH O/P BUFFER FOR MODULES RESIDING IN LINK J
   
   DO  i = opbtop,nxtlin,20
     k = i + 4 + j
     IF (opbuff(k) == nblank) CYCLE
     
!     THIS MODULE IS IN LINK J - PUNCH OUT CALL AND GO TO STATEMENT
     
     n = 2000 + opbuff(i)
     WRITE  (lpch,1150) n,opbuff(i+3),opbuff(i+4)
     1150 FORMAT (i5,1X,5HCALL ,2A4,/6X,8HGO TO 10)
   END DO
 END DO
 j = llink/8
 IF (j > lxlink) CALL page2 (-3)
 IF (j > lxlink) WRITE (optap,1180) swm,j,lxlink
 1180 FORMAT (a27,' 54, THE NUMBER OF MODULES SPECIFIED IN THE LINK ',  &
     'SPECIFICATION TABLE,',i5, /20X,'EXCEEDS THE ALLOWABLE ',  &
     'NUMBER SPECIFIED BY SEMDBD,',i5,1H.)
 CALL pexit
 
!     ERROR MESSAGES -
 
!     NOT ENOUGH OPEN CORE
 
 1200 CALL xgpidg (51,lnkbot-lopncr,0,0)
 GO TO 1250
 
!     NAMED COMMON /XLINK/ IS TOO SMALL
 
 1210 CALL xgpidg (52,0,0,0)
 RETURN
 
!     INCORRECT FORMAT IN ABOVE CARD.
 
 1220 CALL xgpidg (53,0,0,0)
 ASSIGN 280 TO irtn
 GO TO 230
 
!     ERROR IN PARAMETER SECTION OF MPL TABLE
 
 1230 CALL xgpidg (49,mplpnt,mpl(mplpnt+1),mpl(mplpnt+2))
 GO TO 150
 
!     FATAL ERROR IN MPL TABLE
 
 1240 CALL xgpidg (49,mplpnt,mpl(mplpnt+1),mpl(mplpnt+2))
 GO TO 1250
 
!     FATAL ERROR EXIT
 
 1250 xsys(3) = 3
 RETURN
END SUBROUTINE xgpibs
