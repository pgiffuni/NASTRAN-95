SUBROUTINE sdr2b
     
!     SDR2B PROCESSES THE EST. FOR EACH ELEMENT IN THE MASTER SET,
!     PRELIMINARY COMPUTATIONS ARE MADE. IF THE PROBLEM CONTAINS EXTRA
!     POINTS, SIL NOS. ARE CONVERTED TO SILD NOS. THE DATA IS WRITTEN
!     ON ESTA FOR INPUT TO SDR2D WHERE FINAL STRESS AND FORCE RECOVERY
!     COMPUTATIONS ARE MADE.
 
 
 IMPLICIT INTEGER (a-z)
 LOGICAL :: anyout,axic  ,heat  ,reject,strain
 INTEGER :: NAME(2)
 INTEGER :: mmre(2)
 REAL :: scrtch,zz(1) ,bufr(1)
 DIMENSION       kdefrm(2)
 CHARACTER (LEN=27) :: swm
 CHARACTER (LEN=25) :: sfm
 CHARACTER (LEN=29) :: uim
 CHARACTER (LEN=25) :: uwm
 CHARACTER (LEN=23) :: ufm
 COMMON /xmssg / ufm   ,uwm   ,uim   ,sfm   ,swm
 COMMON /zzzzzz/ z(1)
 COMMON /system/ ksystm(63)
 COMMON /BLANK / app(2),sort2 ,idummy(7)    ,strain
 COMMON /sdr2x1/ ieigen,ieldef,itload,isymfl,iloads,idispl,istr  ,  &
     ielf  ,iacc  ,ivel  ,ispcf ,ittl  ,ilsym
 COMMON /sdr2x2/ casecc,cstm  ,mpt   ,dit   ,eqexin,sil   ,gptt  ,  &
                 edt   ,bgpdt ,pg    ,qg    ,ugv   ,est   ,phig  ,  &
                 eigr  ,opg1  ,oqg1  ,ougv1 ,oes1  ,oef1  ,pugv1 ,  &
                 oeigr ,ophig ,pphig ,esta  ,gptta ,harms
 COMMON /hmatdd/ ihmat ,nhmat ,mptmpt,idit
 COMMON /gpta1 / nelem ,last  ,incr  ,elem(1)
 COMMON /sdr2x4/ nam(2), end  ,mset  ,icb(7),ocb(7),mcb(7),dtype(8), &
                 icstm ,ncstm ,ivec  ,ivecn ,temp  ,deform, file,    &
                 buf1  ,buf2  ,buf3  ,buf4  ,buf5  ,any   ,all   ,   &
                 tloads,eldef ,symflg,branch,ktype ,loads ,spcf  ,  &
                 displ ,vel   ,acc   ,stress,force ,kwdest,kwdedt,  &
                 kwdgpt,kwdcc ,nrigds,sta(2),rei(2),ds0(2),ds1(2),  &
                 frq(2),trn(2),bk0(2),bk1(2),cei(2),pla(22)      ,  &
                 nrings,nharms,axic  ,knset ,isopl ,strspt,ddrmm , isopl8
COMMON /sdr2x5/ buf(100),bufa(100)  ,bufb(4176)
COMMON /sdr2x6/ scrtch(300)
COMMON /names / rd    ,rdrew ,wrt   ,wrtrew,clsrew
EQUIVALENCE     (ksystm( 1),sysbuf) ,(ksystm( 2),ioutpt),  &
                (ksystm(55),iprec ) ,(ksystm(56),itherm),  &
                (z(1)      ,zz( 1)) ,(bufr(1)   ,buf(1))
DATA   NAME   / 4HSDR2, 4HB   /      ,star  / 4H* *      /
DATA   kdefrm / 104,1/
DATA   iz1st  / 1    /
!WKBI 7/94 SPR 94007
DATA   mmre   / 4HMMRE, 4HIGEN /
!            IZ1ST  IS THE START OF OPEN CORE AVAILABLE


!     IF APPROACH IS COMPLEX EIGENVALUES, FREQUENCY OR TRANSIENT
!     RESPONSE, TEST FOR EXTRA POINTS. IF PRESENT, READ EQUIVALENCE
!     TABLE (SIL,SILD) INTO CORE.

CALL delset
heat = .false.
IF (itherm /= 0) heat = .true.
isopl = 0
icstm = iz1st
m8    =-8
noep  = 0
!WKBR 7/94 SPR 94007
!     IF (APP(1).EQ.CEI(1) .OR. APP(1).EQ.FRQ(1) .OR. APP(1).EQ.TRN(1))
!    1    GO TO 20
IF (app(1) == cei(1) .OR. app(1) == frq(1) .OR. app(1) == trn(1)  &
    .OR. app(1) == mmre(1) )GO TO 20
GO TO 40
20 icb(1) = sil
CALL rdtrl (icb)
noep = icb(3)
IF (noep == 0) GO TO 40
FILE = sil
CALL OPEN (*560,sil,z(buf1),rdrew)
CALL fwdrec (*570,sil)
CALL fwdrec (*570,sil)
CALL READ (*570,*30,sil,z,buf2,1,nsil)
CALL mesage (m8,0,nam)
30 CALL CLOSE (sil,clsrew)
knsil = nsil/2
icstm = nsil + 1
IF (nsil < mset) GO TO 40
mset = buf2 - 1
all  = 1

!     READ THE CSTM INTO CORE (IF PRESENT).

40 ncstm = 0
FILE  = cstm
CALL OPEN (*60,cstm,z(buf1),rdrew)
CALL fwdrec (*570,cstm)
CALL READ (*570,*50,cstm,z(icstm),buf2-icstm,1,ncstm)
CALL mesage (m8,0,nam)
50 CALL CLOSE (cstm,clsrew)
CALL pretrs (z(icstm),ncstm)
60 imat = icstm + ncstm
IF (imat < mset) GO TO 70
mset = buf2 - 1
all  = 1

!     READ MATERIAL PROPERTY DATA INTO CORE.

70 n1mat = buf2 - imat
IF (.NOT.heat) GO TO 77

!     FOR HEAT PROBLEMS ONLY, -HMAT- ROUTINE IS USED.

ihmat = imat
nhmat = buf1 + sysbuf
mptmpt= mpt
idit  = dit
CALL prehma (z)
n2mat = nhmat - ihmat+1 - 2*(sysbuf+1)
GO TO 78

77 CALL premat (z(imat),z(imat),z(buf1),n1mat,n2mat,mpt,dit)
78 IF (imat+n2mat < mset) GO TO 80
mset = buf2 - 1
all  = 1

!     OPEN EST AND ESTA.

80 FILE = est
CALL OPEN (*620,est,z(buf1),rdrew)
CALL fwdrec (*570,est)
FILE = esta
CALL OPEN (*560,esta,z(buf2),wrtrew)
FILE   = est
kwdest = 0
kwdedt = 0
kwdgpt = 0

!     READ ELEMENT TYPE. SET PARAMETERS AS A FUNCTION OF ELEM TYPE.

90 CALL READ (*430,*580,est,eltype,1,0,flag)
IF (eltype < 1 .OR. eltype > nelem) GO TO 3800
anyout = .false.
ipr = iprec
IF (ipr /= 1) ipr = 0
jltype = 2*eltype - ipr
ielem  = (eltype-1)*incr
nwds   = elem(ielem+12)
nwdsa  = elem(ielem+17)
IF (heat) nwdsa = 142
ngps   = elem(ielem+10)

!     READ DATA FOR AN ELEMENT.
!     DETERMINE IF ELEMENT BELONGS TO MASTER SET.

100 CALL READ (*570,*420,est,buf,nwds,0,flag)
DO  i = 1,nwds
  scrtch(100+i) = bufr(i)
END DO
strspt = 0
isopl  =-1
idsave = buf(1)
IF (all /= 0) GO TO 110
itabl = mset
kn  = knset
l   = 1
n12 = 1
ASSIGN 100 TO ret1
IF (.NOT. axic) GO TO 630

!     DECODE ELEMENT ID SINCE THIS IS A CONICAL SHELL PROBLEM

buf(1) = buf(1)/1000
GO TO 630

!     CALL APPROPRIATE ELEMENT SUBROUTINE.

110 CONTINUE
buf(1) = idsave

IF (.NOT.strain) GO TO 112

!     IF THE STRAIN FLAG IS TURNED ON, IGNORE ALL ELEMENTS
!WKBR NCL93012 3/94 EXCEPT CTRIA1, CTRIA2, CQUAD1 AND CQUAD2 ELEMENTS
!     EXCEPT CTRIA1, CTRIA2, CTRIA3, CQUAD1, CQUAD2 AND CQUAD4 ELEMENTS

IF (eltype == 6 .OR. eltype == 17 .OR. eltype == 18 .OR. &
!WKBR NCL93012 3/94     1    ELTYPE.EQ.19) GO TO 112  &
    eltype == 19 .OR. eltype == 64 .OR. eltype == 83) GO TO 112
WRITE  (ioutpt,111) swm,elem(ielem+1),elem(ielem+2)
111 FORMAT (a27,', STRAIN REQUEST FOR ',2A4,' ELEMENTS WILL', /5X,  &
    'NOT BE HONORED AS THIS OUTPUT IS NOT DEFINED FOR THIS ', 'ELEMENT TYPE.')
CALL fwdrec (*570,est)
GO TO 420

112 IF (heat) GO TO 389
local = jltype - 100
IF (local > 0) THEN
  GO TO   115
END IF

!     PAIRED -GO TO- ENTRIES PER ELEMENT SINGLE/DOUBLE PRECISION

!             1 CROD      2 C.....    3 CTUBE     4 CSHEAR    5 CTWIST
114 GO TO (120,  120,  380,  380,  140,  140,  150,  150,  160,  160  &
!             6 CTRIA1    7 CTRBSC    8 CTRPLT    9 CTRMEM   10 CONROD  &
,      180,  180,  190,  190,  200,  200,  210,  210,  120,  120      &
!            11 ELAS1    12 ELAS2    13 ELAS3    14 ELAS4    15 CQDPLT  &
,      220,  220,  230,  230,  240,  240,  250,  250,  270,  270      &
!            16 CQDMEM   17 CTRIA2   18 CQUAD2   19 CQUAD1   20 CDAMP1  &
,      280,  280,  290,  290,  300,  300,  310,  310,  380,  380      &
!            21 CDAMP2   22 CDAMP3   23 CDAMP4   24 CVISC    25 CMASS1  &
,      380,  380,  380,  380,  380,  380,  380,  380,  380,  380      &
!            26 CMASS2   27 CMASS3   28 CMASS4   29 CONM1    30 CONM2  &
,      380,  380,  380,  380,  380,  380,  380,  380,  380,  380      &
!            31 PLOTEL   32 C.....   33 C.....   34 CBAR     35 CCONE  &
,      380,  380,  380,  380,  380,  380,  330,  330,  340,  340      &
!            36 CTRIARG  37 CTRAPRG  38 CTORDRG  39 CTETRA   40 CWEDGE  &
,      350,  350,  360,  360,  370,  370,  371,  371,  372,  372      &
!            41 CHEXA1   42 CHEXA2   43 CFLUID2  44 CFLUID3  45 CFLUID4  &
,      373,  373,  374,  374,  380,  380,  380,  380,  380,  380      &
!            46 CFLMASS  47 CAXIF2   48 CAXIF3   49 CAXIF4   50 CSLOT3  &
,      380,  380,  375,  375,  376,  376,  377,  377,  378,  378      &
 ), jltype


!            51 CSLOT4   52 CHBDY    53 CDUM1    54 CDUM2    55 CDUM3
115 GO TO (379,  379,  380,  380,  451,  451,  452,  452,  453,  453      &
!            56 CDUM4    57 CDUM5    58 CDUM6    59 CDUM7    60 CDUM8  &
,      454,  454,  455,  455,  456,  456,  457,  457,  458,  458      &
!            61 CDUM9    62 CQDMEM1  63 CQDMEM2  64 CQUAD4   65 CIHEX1  &
,      459,  459,  460,  460,  461,  461,  462,  462,  383,  383      &
!            66 CIHEX2   67 CIHEX3   68 CQUADTS  69 CTRIATS  70 CTRIAAX  &
,      383,  383,  383,  383,  465,  465,  466,  466,  467,  467      &
!            71 CTRAPAX  72 CAERO1   73 CTRIM6   74 CTRPLT1  75 CTRSHL  &
,      468,  468,  380,  380,  469,  469,  470,  470,  471,  471      &
!            76 CFHEX1   77 CFHEX2   78 CFTETRA  79 CFWEDGE  80 CIS2D8  &
,      380,  380,  380,  380,  380,  380,  380,  380,  472,  472      &
!            81 CELBOW   82 CFTUBE   83 CTRIA3  &
,      473,  473,  380,  380,  463,  463  &
 ), local

120 CALL srod1
GO TO 390
140 CALL stube1
GO TO 390
150 k = 4
GO TO 170
160 k = 5
170 CALL spanl1 (k)
GO TO 390
180 k = 1
GO TO 320
190 CALL strbs1 (0)
GO TO 390
200 CALL strpl1
GO TO 390
210 CALL strme1 (0)
GO TO 390
220 k = 1
GO TO 260
230 k = 2
GO TO 260
240 k = 3
GO TO 260
250 k = 4
260 CALL selas1 (k)
GO TO 390
270 CALL sqdpl1
GO TO 390
280 CALL sqdme1
GO TO 390
290 k = 2
GO TO 320
300 k = 4
GO TO 320
310 k = 3
320 CALL strqd1 (k)
GO TO 390
330 CALL sbar1
GO TO 390
340 CALL scone1
GO TO 390
350 CALL strir1
GO TO 390
360 CALL strap1
GO TO 390
370 CALL stord1
GO TO 390
371 CALL ssold1 (1)
GO TO 390
372 CALL ssold1 (2)
GO TO 390
373 CALL ssold1 (3)
GO TO 390
374 CALL ssold1 (4)
GO TO 390
375 k = 0
GO TO 381
376 k = 1
GO TO 381
377 k = 2
381 CALL saxif1 (k)
GO TO 390
378 k = 0
GO TO 382
379 k = 1
382 CALL sslot1 (k)
GO TO 390
383 CONTINUE
CALL sihex1 (eltype-64,strspt,nip)
IF (strspt >= nip**3+1) strspt = 0
GO TO 390
451 CALL sdum11
GO TO 391
452 CALL sdum21
GO TO 391
453 CALL sdum31
GO TO 391
454 CALL sdum41
GO TO 391
455 CALL sdum51
GO TO 391
456 CALL sdum61
GO TO 391
457 CALL sdum71
GO TO 391
458 CALL sdum81
GO TO 391
459 CALL sdum91
GO TO 391
460 CALL sqdm11
GO TO 390
461 CALL sqdm21
GO TO 390
462 CALL squd41
GO TO 390
463 CALL stri31
GO TO 390
465 CONTINUE
GO TO 390
466 CONTINUE
GO TO 390
467 CALL strax1
GO TO 390
468 CALL stpax1
GO TO 390
469 CALL strm61
GO TO 390
470 CALL strp11
GO TO 390
471 CALL  strsl1
GO TO 390
472 CALL ss2d81
isopl8 = 8
GO TO 390
473 CALL selbo1
GO TO 390

!     ELEMENT UNDEFINE TO SDR2BD

3800 WRITE  (ioutpt,385) star,star,eltype
GO TO 388
380 WRITE  (ioutpt,385) swm,elem(ielem+1),elem(ielem+2),eltype
385 FORMAT (a27,' 2184,  STRESS OR FORCE REQUEST FOR ELEMENT ',2A4,  &
    ' (NASTRAN ELEM. TYPE =',i4,1H), /5X,'WILL NOT BE HONORED'  &
    ,       ' AS THIS ELEMENT IS NOT A STRUCTURAL ELEMENT.')
388 CALL fwdrec (*570,est)
GO TO 420

!     HEAT PROBLEMS (ALL ELEMENTS).

389 CALL sdhtf1 (eltype,reject)
IF (eltype < 65 .OR.  eltype > 67) GO TO 3890
IF (eltype == 65 .AND. strspt >= 9) strspt = 0
IF (strspt >= 21) strspt = 0
3890 CONTINUE
IF (.NOT.reject) GO TO 390
CALL fwdrec (*570,est)
GO TO 420

!     IF EXTRA POINTS PRESENT, CONVERT SIL NOS. TO SILD NOS.

391 nwdsa = elem(ielem+17)
390 IF (noep == 0) GO TO 410
n = ngps + 101
itabl = 1
kn  = knsil
n12 = 2
ASSIGN 740 TO ret1
l = 102
!WKBNB 7/94 SPR 94006
! REMOVE COMPONENT FROM SIL AND THEN ADD AFTER SILD NUMBER FOUND FOR
! CELAS1 AND CELAS2 ELEMENTS-SEE SUBROUTINE SELAS1
IF ( eltype == 11 ) GO TO 392
IF ( eltype == 12 ) GO TO 393
GO TO 394
! SET SIL NUMBER TO SIL OF GRID POINT WITHOUT COMPONENT CODE INCLUDED FOR
! CELAS1 SO SIL NUMBER CAN BE FOUND IN SILD
392 buf(l)   = buf(2)
buf(l+1) = buf(3)
GO TO 394
! SET SIL NUMBER TO SIL OF GRID POINT WITHOUT COMPONENT CODE INCLUDED FOR
! CELAS2 SO SIL NUMBER CAN BE FOUND IN SILD
393 buf(l)   = buf(3)
buf(l+1) = buf(4)
GO TO 394
394 CONTINUE
!WKBNE 7/94 SPR 94006
IF (buf(l) == 0) GO TO 400
GO TO 630
400 l = l + 1
!WKBR 7/94 SPR 94006 IF (L      .GT. N) GO TO 410
IF (l      > n) GO TO 401
IF (buf(l) == 0) GO TO 400
GO TO 630
!WKBNB 7/94 SPR94006
401 CONTINUE
IF ( eltype == 11 ) GO TO 402
IF ( eltype == 12 ) GO TO 403
GO TO 404
! ADD COMPONENT CODES FOR SILD NUMBERS FOR CELAS1
402 IF ( buf(4) /= 0 ) buf(102) = buf(102) + buf(4) - 1
IF ( buf(5) /= 0 ) buf(103) = buf(103) + buf(5) - 1
GO TO 404
! ADD COMPONENT CODES FOR SILD NUMBERS FOR CELAS2
403 IF ( buf(5) /= 0 ) buf(102) = buf(102) + buf(5) - 1
IF ( buf(6) /= 0 ) buf(103) = buf(103) + buf(6) - 1
404 CONTINUE
!WKBNE 7/94 SPR 94006

!     WRITE ELEMENT COMPUTATIONS ON ESTA. GO TO READ ANOTHER ELEMENT.

410 IF (anyout) GO TO 411
CALL WRITE (esta,eltype,1,0)
kwdest = kwdest + 2
anyout = .true.
411 CALL WRITE (esta,bufa,nwdsa,0)

!     DIAG 20 OUTPUT ONLY

!     CALL BUG (4HESTA,0,BUFA,NWDSA)

kwdest = kwdest + nwdsa
IF (strspt == 0) GO TO 100
strspt = strspt + 1
GO TO 112

!     CLOSE RECORD FOR CURRENT ELEMENT TYPE.
!     GO TO READ ANOTHER ELEM TYPE.

420 IF (anyout) CALL WRITE (esta,0,0,1)
GO TO 90

!     CLOSE FILES.

430 CALL CLOSE (est ,clsrew)
CALL CLOSE (esta,clsrew)

!     IF ELEMENT DEFORMATIONS, DETERMINE MAXIMUM NO. OF
!     WORDS IN ANY ONE DEFORMATION SET.

IF (eldef == 0) RETURN
CALL preloc (*620,z(buf1),edt)
CALL locate (*600,z(buf1),kdefrm,flag)
id = 0
k  = 0
530 CALL READ (*600,*550,edt,buf,3,0,flag)
IF (buf(1) == id) GO TO 540
kwdedt = MAX0(kwdedt,k)
k  = 3
id = buf(1)
GO TO 530
540 k  = k + 3
GO TO 530
550 kwdedt = MAX0(kwdedt,k)
CALL CLOSE (edt,clsrew)
RETURN


!     FATAL FILE ERRORS.

560 n = -1
GO TO 590
570 n = -2
GO TO 590
580 n = -3
GO TO 590
590 CALL mesage (n,FILE,nam)

!     ABNORMAL RETURN FROM SDR2B.

600 CALL CLOSE (edt,clsrew)
eldef = 0
GO TO 620
620 CALL mesage (30,79,0)
stress = 0
force  = 0
any    = 0
RETURN

!     BINARY SEARCH ROUTINE

630 klo = 1
khi = kn
640 k   = (klo+khi+1)/2
650 kx  = itabl + n12*(k-1)
IF (buf(l)-z(kx) < 0.0) THEN
  GO TO   660
ELSE IF (buf(l)-z(kx) == 0.0) THEN
  GO TO   720
ELSE
  GO TO   670
END IF
660 khi = k
GO TO 680
670 klo = k
680 IF (khi-klo-1  < 0) THEN
  GO TO   730
ELSE IF (khi-klo-1  == 0) THEN
  GO TO   690
ELSE
  GO TO   640
END IF
690 IF (k == klo) GO TO 700
k = klo
GO TO 710
700 k = khi
710 klo = khi
GO TO 650
720 IF (n12 == 1) GO TO 110
buf(l) = z(kx+1)
GO TO 400
730 GO TO ret1, (100,740)
740 CALL mesage (-61,0,NAME)
GO TO 740

END SUBROUTINE sdr2b
