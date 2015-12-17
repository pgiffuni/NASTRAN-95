SUBROUTINE scan
     
!     THIS IS THE MAIN DRIVER FOR THE OUTPUT SCAN MODULE - SCAN
 
!     THIS SCAN MODULE CAN BE CALLED DIRECTLY FROM ALL RIGID FORMATS, OR
!     BY USER DMAP ALTER. THE CALLING INSTRUCTIONS ARE
 
!     (THREE INPUT FILES IF CALLED BY RIGID FORMAT VIA SCAN INPUT CARDS)
!     (1) FORCE AND STRESS SCAN -
!     SCAN  CASECC,OESI,OEFI/OESFI/*RF*  $    (WHERE I=1, OR 2)
!           OR
!     SCAN  CASECC,OESI,OEFI/OESFI/*OLI* $    FOR ON-LINE SCAN
 
!         . IF INPUT FILES ARE OES1, OEF1, SORT1 TYPE DATA ARE SCANNED
!         . IF INPUT FILES ARE OES2, OEF2, SORT2 TYPE DATA ARE SCANNED
 
!     (ONE INPUT FILE ONLY IF CALLED BY USER VIA DMAP ALTER)
!     (2) STRESS SCAN -
!     SCAN, ,OESI, /OESFI/C,N,ELEMENT/C,N,ICOMP/C,N,NTOP/C,N,AMAX/
!           C,N,AMIN/C,N,IBEG/C,N,IEND/C,N,ICOMPX $
!     OR (3) FORCE SCAN -
!     SCAN, ,,OEFI /OESFI/C,N,ELEMENT/C,N,ICOMP/C,N,NTOP/C,N,AMAX/
!           C,N,AMIN/C,N,IBEG/C,N,IEND/C,N,ICOMPX $
 
!         . FOR SORT1 TYPE DATA, OESI AND OEFI ARE OES1 AND OEF1, AND
!           IBEG AND IEND ARE RANGE OF ELEMENT IDS TO BE SCANNED
!         . FOR SORT2 TYPE DATA, OESI AND OEFI ARE OES2 AND OEF2, AND
!           IBEG AND IEND ARE RANGE OF SUBCASE IDS TO BE SCANNED
!         . IF IBEG AND IEND ARE NOT GIVEN, ALL IDS IMPLY
 
!         . OESB1, OESC1, OEFB1, AND OEFC1 CAN BE USED IN LIEU OF OES1
!           AND OEF1. SIMILARLY, OESC2 AND OEFC2  FOR OES2 AND OEF2
 
!     INPUT  FILES  - CASECC, OES1, OEF1, (OR OES2, OEF2)
!                     (OESB1, OESC1, OEFB1, OEFC1, OESB2, OEFB2 CAN BE
!                     USED INSTEAD)
!     OUTPUT FILE   - OESF1 (OR OESF2)
!     SCRATCH FILE  - SCR1
 
!     THIS SCAN MODULE SHOULD BE FOLLOWED BY OFP TO PRINT SCAN RESULTS
!     OFP  OESFI,,,,, //S,N,CARDNO $
 
!     PARAMETERS -
 
!           ELEMENT - ELEMENT NAME IN BCD.  E.G. BAR, CBAR, QUAD2, ETC.
!           ICOMP   - THE OUTPUT FIELD NO. (BY COLUMN, 1 THRU 31) OF
!                     OUTPUT LISTING.
!           ICOMPX  - OUTPUT FIELD NO. CONTINUATION (FROM 32 THRU 62)
!           NTOP    - TOP N VALUES TO BE OUTPUT.  DEFAULT=20
!      AMAX-AMIN    - SCAN VALUES OUTSIDE THIS MAX-MIN RANGE, DEFAULT=0.
!      IBEG,IEND    - SEE EXPLANATION ABOVE
 
!     DEFINITION OF SOME LOCAL VARIABLES
 
!           DEBUG   - USED FOR LOCAL DEBUG
!           S OR F  - STRESS OR FORCE SCAN FLAG
!           NSCAN   - NO. OF SCAN INPUT CARDS IN CASECC
!           SUBC    - CURRENT SUBCASE ID
!           NZ      - TOP OF OPEN CORE, JUST BELOW GINO BUFFERS
!           LCORE   - AVAILABLE CORE FOR STRSCN ROUTINE
!           IOPEN   - INPUT  FILE STATUS FLAG, .T. FOR OPEN, .F. NOT
!           JOPEN   - OUTPUT FILE STATUS FLAG, .T. FOR OPEN, .F. NOT
!           KOPEN   - SCR1   FILE STATUS FLAG, .T. FOR OPEN, .F. NOT
!           LOPEN   - CASECC FILE STATUS FLAG, .T. FOR OPEN, .F. NOT
!           ISET    - SCAN ONLY BY THE SPECIFIED SET OF ELEM. IDS
!                   - ALL IS IMPLIED IF ISET IS NOT GIVEN
!                   - USED ONLY IF SCAN IS CALLED FROM RIGID FORMAT
!      IDUPL,INC    - SET UP COMPONENT FIELDS TO BE REPEATEDLY SCANNED
!                     IDUPL TIMES, WITH FIELD INCREMENT BY INC (RF ONLY)
!      LBEG,LEND    - A LIST OF TO-BE-SCANNED ELEMENT IDS, STORED IN
!                     Z(LBEG) THRU Z(LEND).
!                   - NO SUCH LIST EXISTS IF LBEG.GT.LEND OR LBEG=LEND=0
!           IOPT    - DATA SCAN BY AMAX AND AMIN IF IOPT=1, BY NTOP IF 2
!           ISORT   - SET TO 1 (BY STRSCN) IF DATA TYPE IS IN SORT1
!                     FORMAT, AND SET TO 2 IF SORT2
 
!     WRITTEN BY G.CHAN/SPERRY      OCTOBER 1984
 
!     THIS ROUTINE OPENS AND CLOSES ALL INPUT AND OUTPUT FILES.
!     IT SETS UP THE SCANNING PARAMETERS AND CALL STRSCN TO SCAN THE
!     OUTPUT STRESS OR FORCE DATA
 
!     THE SCAN INPUT CARDS OPERATE VERY SIMILARY TO THE ELEMENT STRESS
!     OR FORCE CARDS. THEY CAN BE PLACED ABOVE ALL SUBCASES, OR INSIDE
!     ANY SUBCASE LEVEL, OR BOTH
!     HOWEVER, UNLIKE THE STRESS OR FORCE CARDS, MULTI-SCAN CARDS ARE
!     ALLOWED, AND THEY DO NOT EXCLUDE ONE ANOTHER.
 
!     MODIFIED IN 10/1989, TO ALLOW SETS TO BE DEFINED BEFORE OR AFTER
!     SCAN CARDS IN CASE CONTROL SECTION
!     (CURRENTLY, THIS MODIFICATION IS OK, BUT IFP1/IFP1H DO NOT ALLOW
!     SET TO BE DEFINED AFTER SCAN. IN FACT, IFP1 DOES NOT ALLOW SET TO
!     BE DEFINED AFTER ANY GUY WHO USES THE SET)
 
 LOGICAL :: debug,    iopen,    jopen,    kopen,    lopen
!WKBI  1/4/94 SPR93010 & 93011
 LOGICAL :: stress,   force,    layerd
!WKBI  1/4/94 SPR93010 & 93011
 INTEGER :: quad4,    tria3
!RLBR 12/29/93 SPR 93010 & 93011
!     INTEGER         CASECC,   OESI,     OEFI,     OESFI,    SCR1,
 INTEGER :: casecc,   oesi(2),  oefi(2),  oesfi(2), scr1,  &
     oufile,   FILE,     sorf,     z(166),   nam(2),  &
     e,        eor,      subc,     osubc,    oel
!RLBNB 12/29/93 SPR 93010 & 93011
 INTEGER :: jelt(2)
!RLBNE 12/29/93 SPR 93010 & 93011
 
 CHARACTER (LEN=27) :: swm
 CHARACTER (LEN=25) :: sfm
 CHARACTER (LEN=29) :: uim
 CHARACTER (LEN=25) :: uwm
 CHARACTER (LEN=23) :: ufm
 COMMON /xmssg / ufm,      uwm,      uim,      sfm,      swm
 COMMON /BLANK / ielt(2),  icomp,    ntop,     amax,     amin,  &
     ibeg,     iend,     icompx
 COMMON /system/ ibuf,     nout,     skp(83),  intra
 COMMON /names / rd,       rdrew,    wrt,      wrtrew,   rew, norew,    eofnrw
 COMMON /gpta1 / nelem,    last,     incr,     e(1)
 COMMON /xscanx/ infile,   oufile,   lcore,    lbeg,     lend,  &
     iopen,    jopen,    iel,      iopt,     iset,  &
     isort,    itrl3,    subc,     osubc,    oel, &
!WKBR 1/4/94 SPR93010 & 93011     3                DEBUG  &
 debug,    lloop,    quad4,    tria3,    stress, force,    layerd
 COMMON /zzzzzz/ core(1)
 EQUIVALENCE     (imax,amax),        (imin,amin),  &
     (idupl,ibeg),       (inc,iend), (core(1),z(1))
!RLBDB 12/29/93 SPR 93010 & 93011
!     DATA            CASECC,   OESI,     OEFI,     OESFI,    SCR1    /
!    1                101,      102,      103,      201,      301     /
!RLBDE 12/29/93 SPR 93010 & 93011
!RLBNB 12/29/93 SPR 93010 & 93011
 DATA            casecc,   oesi(1),  oefi(1),  oesi(2),  oefi(2),  &
     oesfi(1), oesfi(2), scr1    / 101,      102,      103,      104,      105,  &
     201,      202,      301     /
!RLBNE 12/29/93 SPR 93010 & 93011
 DATA            nam,                llc,      eor,      irf     /  &
     4HSCAN,   4H    ,   4HC   ,   1,        4HRF    /
 DATA            iol1,     iol2      / 4HOL1 ,   4HOL2     /
 
 debug = .false.
!WKBNB 1/4/94 SPR93011 & 93010
 quad4 = 0
 tria3 = 0
 
!     ALLOCATE OPEN CORE
 
!RLBNB 12/29/93 SPR 93010 & 93011
 lloop = 1
 jelt(1) = ielt(1)
 jelt(2) = ielt(2)
 10  CONTINUE
!RLBNB 12/29/93 SPR 93010 & 93011
 nz    = korsz(z)
 ibuf1 = nz - ibuf + 1
 ibuf2 = ibuf1 - ibuf
 ibuf3 = ibuf2 - ibuf
 nz    = ibuf3 - 1
 lcore = ibuf2 - 1
 iopen =.false.
 jopen =.false.
 kopen =.false.
 lopen =.false.
 
!     OPEN CASECC AND CHECK SCAN DATA
 
 iset = 0
 IF (ielt(1) /= irf) iset = -2
 IF (ielt(1) == iol1 .OR. ielt(1) == iol2) iset = -3
 IF (iset == -2) GO TO 40
 FILE = casecc
 CALL OPEN (*310,casecc,z(ibuf1),rdrew)
 lopen = .true.
 CALL fwdrec (*320,casecc)
 IF (iset == -3) GO TO 40
 30   CALL READ (*80,*80,casecc,z(1),200,1,l)
 lencc = z(166)
 nscan = z(lencc-1)
 IF (nscan == 0) GO TO 30
 
!     CHECK THE PRESENCE OF STRESS AND/OR FORCE FILE.
!     QUIT IF BOTH ARE PURGED
 
 40   ioes  = 1
 ioef  = 1
!RLBDB 12/29/93 SPR 93010 & 93011
!     Z( 1) = OESI
!     Z(11) = OEFI
!RLBDE 12/29/93 SPR 93010 & 93011
!RLBNB 12/29/93 SPR 93010 & 93011
 z( 1) = oesi(lloop)
 z(11) = oefi(lloop)
!RLBNE 12/29/93 SPR 93010 & 93011
 CALL rdtrl (z( 1))
 CALL rdtrl (z(11))
 IF (z( 1) < 0) ioes = 0
 IF (z(11) < 0) ioef = 0
 IF (ioes+ioef == 0 .AND. iset /= -3) GO TO 300
 
!     OPEN OUTPUT FILE OESFI
 
!RLBDB 12/29/93 SPR 93010 & 93011
!     FILE = OESFI
!     OUFILE = OESFI
!     CALL FNAME (OESFI,Z)
!     CALL OPEN  (*310,OESFI,Z(IBUF2),WRTREW)
!     CALL WRITE (OESFI,Z,2,EOR)
!RLBDE 12/29/93 SPR 93010 & 93011
!RLBNB 12/29/93 SPR 93010 & 93011
 FILE = oesfi(lloop)
 oufile = oesfi(lloop)
 CALL fname (oufile,z)
 CALL OPEN  (*310,oufile,z(ibuf2),wrtrew)
 CALL WRITE (oufile,z,2,eor)
!RLBNE 12/29/93 SPR 93010 & 93011
 jopen  =.true.
 itrl3  = 0
 lx     =-1
 IF (ielt(1) == iol2) lx = -2
 IF (iset == -3) CALL onlins (*280,lx)
 IF (iset /= -2) GO TO 90
 
!     SCAN CALLED BY USER VIA DMAP ALTER (ISET=-2)
!     ============================================
 
 ls   = lcore
 lbeg = 0
 lend = 0
 
!     CHECK USER DMAP ERROR, SET IOPT FLAG, AND INITIALIZE ISCAN ARRAY
!     FOR COMPONENT SPECIFIED.
 
 IF (ioes+ioef > 1) GO TO 400
 IF (amin  >  amax) GO TO 410
 IF (icomp <=     1) GO TO 420
 IF ((amax == 0. .AND. amin == 0.) .AND. ntop == 0) GO TO 430
 IF ((amax /= 0. .OR.  amin /= 0.) .AND. ntop /= 0) GO TO 440
 IF ((ibeg == 0 .AND. iend /= 0) .OR. ibeg > iend .OR.  &
     (ibeg /= 0 .AND. iend == 0)) GO TO 460
 IF ( ibeg == 0 .AND. iend == 0 ) ibeg = -1
 iopt = 1
 IF (ntop > 0) iopt = 2
 
!     DETERMINE ELEMENT TYPE, DROP THE FIRST LETTER C IF NECESSARY
 
 z(1) = irf
 z(2) = irf
 IF (khrfn2(ielt(1),1,1) /= llc) GO TO 50
 z(1) = khrfn3(nam(2),ielt(1),1,1)
 z(1) = khrfn1(z(1),4,ielt(2),1  )
 z(2) = khrfn3(nam(2),ielt(2),1,1)
 50   DO  i = 1,last,incr
   IF (ielt(1) == e(i) .AND. ielt(2) == e(i+1)) GO TO 70
   IF (   z(1) == e(i) .AND.    z(2) == e(i+1)) GO TO 70
 END DO
 GO TO 450
 70   iel = e(i+2)
 
!     SPECIAL HANDLING OF THE QUAD4 AND TRIA3 ELEMENT, STRESS ONLY
!     (THE 2ND, 3RD, 9TH, AND 13TH WORDS IN OES1/OES1L FILES ARE
!     NOT PRINTED. THE 9TH AND 13TH WORDS MAY BE BLANKS OR ASTERISKS)
 
 IF ((iel /= 64 .AND. iel /= 83) .OR. ioes == 0) GO TO 75
!WKBD 1/3/94 SPR93011 & 93011      ICOMP = ICOMP + 2
!WKBD 1/3/94 SPR93010 & 93011      IF (ICOMP .GT. 8) ICOMP = ICOMP + 1
 
!     OPEN INPUT FILE
 
!RLBDB 12/29/93 SPR 93010 & 93011
!75   INFILE = OESI
!     IF (IOES .EQ. 0) INFILE = OEFI
!RLBDE 12/29/93 SPR 93010 & 93011
!RLBNB 12/29/93 SPR 93010 & 93011
 75   infile = oesi(lloop)
 stress = .true.
 force  = .false.
 IF (ioes /= 0) GO TO 76
 stress = .false.
 force  = .true.
 infile = oefi(lloop)
!RLBNE 12/29/93 SPR 93010 & 93011
 76   FILE = infile
 CALL OPEN (*340,infile,z(ibuf1),rdrew)
 iopen = .true.
 
! ... NEXT I/O OPERATION ON INFILE WILL BE IN SUBROUTINE STRSCN
 
!     ALL SET TO GO
 
 j = 1
 IF (ioes == 0) j = 2
 CALL strscn (j)
 GO TO 280
 
 80   CALL CLOSE (casecc,rew)
 lopen = .false.
 RETURN
 
 
!     SCAN IS CALLED BY RIGID FORMAT (ISET .GE. -1)
!     OR CALLED BY INTERACTIVE MODE  (ISET .EQ. -3)
!     =============================================
 
 90   ls = nz
 
!     OPEN SCR1 FILE, SEPERATE SCAN DATA FROM SET DATA IN CASECC, AND
!     SAVE THE COMPLETE SCAN DATA IN SCR1 FILE.
 
 FILE = scr1
 CALL OPEN (*310,scr1,z(ibuf3),wrtrew)
 kopen =.true.
 nscan = 0
 ncase = 0
 nxx   = nz
 IF (intra <= 0) GO TO 95
 nxx = 198
 l   = lx
 IF (lx > 0) GO TO 110
 95   FILE = casecc
 CALL REWIND (casecc)
 CALL fwdrec (*320,casecc)
 
!     READ CASECC AND PROCESS ALL SUBCASES
 
 100  CALL READ (*210,*110,casecc,z(1),nxx,1,l)
 IF (nxx >= 200) GO TO 380
 110  ncase = ncase + 1
 lencc = z(166)
 nscan = z(lencc-1)
 lsem  = z(lencc)
 subc  = z(1)
 
!     PICK UP ALL THE SET ID'S AND THEIR LOCATIONS IN Z ARRAY, Z(L1)
!     THRU Z(LL). SORT, AND CHECK DUPLICATE
 
 jmp = 0
 ii  = lencc + lsem
 l1  = l + 1
 ll  = l
 115  ii  = ii + jmp
 IF (ii >= l) GO TO 120
 jmp = z(ii+2) + 2
 IF (z(ii+1) >= 10000000 .AND. jmp == 8) GO TO 115
 z(ll+1) = z(ii+1)
 z(ll+2) = ii
 ll = ll + 2
 GO TO 115
 120  lll1 = ll - l1 + 1
 ll2  = lll1/2
 IF (debug) WRITE (nout,125) (z(i),i=l1,ll)
 125  FORMAT (' ...SET/@125',/,(10X,i8,' @',i6))
 
 jmp = 0
 ii  = lencc + lsem
 kk  = nz
 IF (ll2 <= 1) GO TO 140
 CALL sort (0,0,2,1,z(l1),lll1)
 j = l1 + 2
 DO  i = j,ll,2
   IF (z(i) == z(i-2)) WRITE (nout,600) uwm,z(i)
 END DO
 
!     PROCESS THE SCAN CARDS
 
!     PICK UP SCAN 8 WORD ARRAY, AND PICK UP SET DATA
!     WRITE TO SCR1 A RECORD (OF EACH SUBCASE) OF THE SCAN INPUT DATA
!     IN REVERSE ORDER (FIRST SCAN CARD LAST, AS SET UP BY CASECC)
 
 140  ii = ii + jmp
 IF (ii >= l) GO TO 190
 jmp = z(ii+2) + 2
 IF (z(ii+1) < 10000000 .OR. jmp /= 8) GO TO 140
 ie  = 0
 iset= z(ii+4)
 IF (iset == -1) GO TO 160
 IF (lll1 <=  0) GO TO 470
 CALL bisloc (*470,iset,z(l1),2,ll2,i)
 ib = z(i+l1) + 2
 ie = z(ib)
 IF (debug) WRITE (nout,145) iset,i,ib,ie
 145  FORMAT (' @145, SET',i8,' FOUND.  I,IB,IE =',3I6)
 kk = kk - ie
 DO  i = 1,ie
   z(kk+i) = z(ib+i)
 END DO
 160  kk = kk - 9
 DO  i = 1,8
   z(kk+i) = z(ii+i)
 END DO
 z(kk+9) = 0
 idupl   = z(kk+8)
 IF (idupl == 0) GO TO 180
!WKBD 1/3/94 SPR93010 & 93011      INC = IDUPL/100
!WKBD 1/3/94 SPR93010 & 93011      Z(KK+8) = MOD(IDUPL,100)
!WKBNB 1/3/94 SPR93010 & 93011
 inc = MOD ( idupl, 100 )
 z(kk+8) = idupl / 100
!WKBNE 1/3/94 SPR93010 & 93011
 z(kk+9) = inc
 180  z(kk+2) = z(kk+2) + 1 + ie
 
!     HERE AT THE TAIL END OF OPEN CORE, WE ACCUMULATE ANOTHER RECORD
!     OF A SCAN DATA SET
!        WORD 1,  10000000 FOR STRESS, OR 20000000 FOR FORCE
!             2,  NO. OF WORDS OF THIS DATA SET (SCAN + SET)
!                 (FIRST 2 WORDS NOT INCLUDED)
!             3,  ELEMENT TYPE NUMERIC CODE
!             4,  SET-ID, OR -1
!             5,  COMPONENT CODE, ICOMP
!             6,  NTOP, OR AMAX
!             7,  -1,   OR AMIN
!             8,  COMPONENT - DUPLICATION, OR ZERO
!             9,  COMPONENT - INCREMENT,   OR ZERO
!        10-END,  SET DATA
!     REPEAT FOR ANOTHER SCAN CARD
 
 
!     SPECIAL HANDLING OF THE QUAD4 AND TRIA3 ELEMENT, STRESS ONLY
!     (THE 2ND, 3RD, 9TH,  AND 13TH WORDS IN OES1/OES1L FILES ARE
!     NOT PRINTED. THE 9TH AND 13TH WORDS MAY BE BLANKS OR ASTERISKS)
!WKBI 12/93 SPR93010 & 93011
!     ABOVE IS TRUE ONLY FOR LAMINATED QUAD4 AND TRIA3)
 
!WKBD 12/31/93 SPR93010 & 93011
!     IF ((Z(KK+3).NE.64 .AND. Z(KK+3).NE.83) .OR. Z(KK+1).NE.10000000)
 IF ((z(kk+3) /= 64 .AND. z(kk+3) /= 83) .OR. z(kk+8) == 0) GO TO 140
!WKBDB 1/3/94 SPR93010 & 93011
!      Z(KK+5) = Z(KK+5) + 2
!      IF (Z(KK+5) .GT. 8) Z(KK+5) = Z(KK+5) + 1
!      IF (Z(KK+9) .NE. 0) Z(KK+9) = Z(KK+9) + 2
!WKBDE 1/3/94 SPR93010 & 93011
 GO TO 140
 
!     AT THE END OF EACH SUBCASE, WE COMPUTE THE TOTAL LENGTH OF THIS
!     SCAN DATA ARRAY, AND WRITE THE ARRAY OUT TO SCR1.  ONE RECORD PER
!     SUBCASE
 
 190  kk = kk - 2
 IF (kk < ll) GO TO 610
 ie = nz - kk
 z(kk+1) = subc
 z(kk+2) = ie - 2
 CALL WRITE (scr1,z(kk+1),ie,1)
 l  = kk + 1
 nn = 200
 IF (debug) WRITE (nout,200) nn,(z(j),j=l,nz)
 200  FORMAT (/,11H scan/debug,i3,  (/2X,13I9))
 IF (intra <= 0 .OR. lx < 200) GO TO 100
 
!     THUS, END OF THE PREPARATION PHASE.  CLOSE CASECC AND SCR1
 
 210  CALL CLOSE (casecc,rew)
 CALL CLOSE (scr1  ,rew)
 kopen =.false.
 lopen =.false.
 
!     NOW, SET UP 2 LOOPS FOR STRESS (10000000) AND FORCE (20000000)
!     OUTPUT SCAN
 
 sorf = 30000000
 220  sorf = sorf - 10000000
 IF (debug) WRITE (nout,225) sorf
 225  FORMAT (///,18H processing series,i15 /1X,8(4H====),/)
 IF (iopen) CALL CLOSE (infile,rew)
 iopen = .false.
 IF (sorf == 10000000 .AND. ioes == 0) GO TO 220
 IF (sorf == 20000000 .AND. ioef == 0) GO TO 220
 IF (sorf <= 0) GO TO 280
 
!     OPEN INPUT FILES
 
!RLBDB 12/29/93 SPR 93010 & 93011
!     INFILE = OESI
!     IF (SORF .GE. 20000000) INFILE=OEFI
!RLBDE 12/29/93 SPR 93010 & 93011
!RLBNB 12/29/93 SPR 93010 & 93011
 infile = oesi(lloop)
 stress = .true.
 force  = .false.
 IF (sorf < 20000000) GO TO 226
 stress = .false.
 force  = .true.
 infile=oefi(lloop)
!RLBNE 12/29/93 SPR 93010 & 93011
 226   FILE = infile
 CALL OPEN (*310,infile,z(ibuf1),rdrew)
 iopen = .true.
! ... NEXT I/O OPERATION ON INFILE WILL BE IN SUBROUTINE STRSCN
 
!     NOW, LOAD THE SCAN DATA PREVIOUSLY SAVED IN SCR1, TO THE TAIL END
!     OF THE OPEN CORE.
!     ONE OR MORE SCAN CARDS MAY BE PRESENT IN  ONE SUBCASE
!     SET UP POINTERS IN FRONT OF THE SCAN DATA, SO THAT FIRST SCAN
!     INPUT CARD WILL BE PROCESS FIRST, SECOND CARD SECOND, ETC.
!     NOTE - USE SUBCASE 1 SCAN DATA IF OUTPUT IS SORT 2 TYPE
!            (IF SUBCASE 1 DOES NOT HAVE SCAN DATA, USE NEXT SUBCASE)
 
 FILE = scr1
 IF (.NOT.kopen) CALL OPEN (*310,scr1,z(ibuf3),rdrew)
 IF (     kopen) CALL REWIND (scr1)
 kopen =.true.
 isort = 0
 osubc = 0
 oel   = 0
 
 DO  ii = 1,ncase
   IF (isort == 2) GO TO 220
   CALL READ (*320,*330,scr1,z(1),2,0,l)
   j = z(2)
   IF (j == 0) GO TO 260
   subc = z(1)
   ls   = nz - j
   CALL READ (*320,*330,scr1,z(ls+1),j,1,l)
   LE = ls
   i  = ls
   230  z(ls) = i
   ls = ls - 1
   i  = i + z(i+2) + 2
   IF (i < nz) GO TO 230
   lcore = ls
   j  = ls + 1
   kk = 230
   IF (debug) WRITE (nout,200) kk,subc,(z(i),i=j,nz)
   
!     NOW IS THE TIME TO SET THE SCAN PARAMETERS FOR EACH SCAN CARD
!     WITHIN A SUBCASE, AND CALL STRSCN TO SCAN THE OUTPUT DATA
   
   i = ls
   240  i = i + 1
   IF (i > LE) CYCLE
   ib = z(i)
   IF (z(ib+1) /= sorf) GO TO 240
   jmp   = z(ib+2)
   iel   = z(ib+3)
! ONLY QUAD4 (=64) AND TRIA3 (=83) ARE VALID FOR LLOOP=2
   IF ( lloop == 2 .AND. iel /= 64 .AND. iel /= 83 ) GO TO 240
   iset  = z(ib+4)
   icomp = z(ib+5)
   ntop  = z(ib+6)
   imax  = z(ib+6)
   imin  = z(ib+7)
   idupl = z(ib+8)
   inc   = z(ib+9)
   iopt  = 1
   IF (imin == -1) iopt = 2
   IF (iopt /=  2) ntop = 0
   lbeg = lcore
   lend = lcore - 1
   IF (iset == -1) GO TO 250
   lbeg = ib + 10
   lend = ib + jmp + 2
   250  j    = (iel-1)*incr
   ielt(1) = e(j+1)
   ielt(2) = e(j+2)
   IF (debug) WRITE (nout,255) ielt,(z(ib+j),j=3,9),iopt,lbeg,lend, ii,subc
   255  FORMAT (/5X,16HDEBUG/scan255 - ,2A4,/5X,12I9)
   CALL strscn (sorf/10000000)
   IF (iopt < 0) GO TO 480
   GO TO 240
   260  CALL fwdrec (*320,scr1)
 END DO
 
!     GO BACK TO PROCESS NEXT INPUT FILE
 
 GO TO 220
 
!     ALL SCAN DONE.  WRITE OUTPUT FILE TRAILERS AND CLOSE ALL FILES
 
 280  IF (itrl3 <= 0) GO TO 300
!RLBR 12/29/93 SPR 93010 & 93011
!     Z(1) = OESFI
 z(1) = oesfi(lloop)
 z(2) = 1
 z(3) = itrl3
 DO  i = 4,7
   z(i) = 0
 END DO
 CALL wrttrl (z(1))
 
 300  IF (iopen) CALL CLOSE (infile,rew)
 IF (jopen) CALL CLOSE (oufile,rew)
 IF (kopen) CALL CLOSE (scr1  ,rew)
 IF (lopen) CALL CLOSE (casecc,rew)
!RLBNE 12/29/93 SPR 93010 & 93011
 IF (lloop == 2) GO TO 305
 lloop = 2
 ielt(1) = jelt(1)
 ielt(2) = jelt(2)
 GO TO 10
 305  CONTINUE
 IF ( quad4 == -1 ) WRITE ( nout, 605 ) 'quad4'
 IF ( tria3 == -1 ) WRITE ( nout, 605 ) 'tria3'
 605  FORMAT(//' SCAN MODULE DID NOT FIND ELEMENT ',a5,  &
     ' IN USER OUTPUT REQUESTS.',/  &
     ,' POSSIBLY WRONG COMPONENT SPECIFIED FOR LAYERED OR '  &
     ,'NON-LAYERED CASE',//)
!RLBNE 12/29/93 SPR 93010 & 93011
 RETURN
 
!     FILE ERRORS
 
 310  j = -1
 GO TO 350
 320  j = -2
 GO TO 350
 330  j = -3
 GO TO 350
 340  CONTINUE
 GO TO 70
 350  CALL mesage (j,FILE,nam)
 380  j = -8
 GO TO 350
 
!     ERROR MESSAGES
 
 400  WRITE (nout,500)
 GO TO 490
 410  WRITE (nout,510)
 GO TO 490
 420  WRITE (nout,520)
 GO TO 490
 430  WRITE (nout,530)
 GO TO 490
 440  WRITE (nout,540)
 GO TO 490
 450  WRITE (nout,550) ielt
 GO TO 490
 460  WRITE (nout,560) sfm,ielt,ibeg,iend
 GO TO 490
 470  WRITE (nout,570) uwm,iset
 GO TO 140
 480  WRITE (nout,580) iopt
 490  WRITE (nout,590) swm
 GO TO 280
 
 500  FORMAT (//5X,48HONLY one INPUT FILE allowed from scan dmap alter)
 510  FORMAT (//5X,21HAMAX-amin range error)
 520  FORMAT (//5X,35HFIELD component specification error)
 530  FORMAT (//5X,30HNO amax-amin OR ntop specified)
 540  FORMAT (//5X,46HSPECIFY either amax-amin OR ntop, but NOT both,  &
     /5X,21H(ntop=20  by default))
 550  FORMAT (//5X,22HELEMENT mis-spelled - ,2A4)
560  FORMAT (a25,' - SCANNING ',2A4,' ELEMENT. IBEG-IEND OUT OF RANGE',  &
    '.  SCAN ABORTED')
570  FORMAT (a25,' FROM SCAN, SET',i9,' NOT FOUND')
580  FORMAT (//5X,44HUSER error.  illegal INPUT FILE sent TO scan,i6)
590  FORMAT (a27,' FROM SCAN.  CASE ABORTED ***')
600  FORMAT (a25,' FROM SCAN, DUPLICATE SET',i9)

610  CALL mesage (8,0,nam)
RETURN
END SUBROUTINE scan
