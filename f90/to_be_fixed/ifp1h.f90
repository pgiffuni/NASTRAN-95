SUBROUTINE ifp1h (i81,nz,j400)
     
!     THIS ROUTINE PROCESSES THE SCAN CARD IN CASE CONTROL SECTION
 
!     WRITTEN BY G.CHAN/SPERRY,  OCTOBER 1984
 
!     PROGRAM METHOD
 
!     A 'SCAN(HELP)' INPUT CARD WILL SET J400 TO 2, AND ANY ERROR IN
!     A SCAN INPUT CARD WILL SET J400 TO 1. NON-ZERO J400 WILL CAUSE
!     SCAN COMPONENT KEY-WORDS TO BE PRINTED.
 
!     THE SCAN INPUT CARDS, AND THEIR DATA, ARE DECODED AND SAVED IN
!     CASECC FILE AS SETS OF PSEUDO SET COMMANDS (SET ID OF 10000000 FOR
!     STRESS, AND 20000000 FOR FORCE). IN THIS WAY, THE SCAN CARDS CAN
!     BE USED IN ALL SUBCASE LEVELS, OR ABOVE-SUBCASE LEVEL, SIMILAR TO
!     THE ELEM. STRESS AND ELEM. FORCE CARDS IN THE CASE CONTROL SECTION
!     HOWEVER, MULTIPLE SCAN CARDS CAN BE USED IN ALL SUBCASE LEVELS,
!     AND WITHIN EACH SUBCASE
 
!     ELEM. NAME CAN BE SPECIFIED WITH OR WITHOUT THE LEADING LETTER C
!     E.G.  BAR, CBAR, QUAD2, CQUAD2
 
!     SCAN COMPONENTS CAN BE REQUESTED BY ACTUAL OUTPUT COLUMN NUMBER(S)
!     OR BY COMPONENT KEYWORD(S)
!     IF THE ACTUAL OUTPUT COLUMN IS NOT IN THE SAME WORD ORDER AS IN
!     THE OUTPUT PRINT FILE (E.G. OES1L FOR THE QUAD4 LAYER), THE ACTUAL
!     COLUMN COUNT AS IT APPEARS IN THE PRINTOUT, IS USED HERE. ANY
!     DISCREPANCY SHOULD BE HANDLED BY SCAN OR STRSCN ROUTINES.
 
!     A LIST OF KEYWORDS WILL BE PRINTED AUTOMATICALLY IF ELEM. NAME OR
!     COMPONENT KEYWORD ARE MISSPELLED OR MISSSING
 
!     THIS LIST IS ALSO PRINTED IF A  SCAN (HELP) CARD IS IN INPUT DECK
 
!     THIS ROUTINE MAY ISSUE THE FOLLOWING ERROR MESSAGES -
 
!        604 - NON-INTEGER IN INTEGER FIELD
!        608 - SET NOT DEFINED
!        617 - IMPROPER FORMAT
!        634 - KEYWORD INSIDE BRACKET IS ILLEGAL OR MISSPELLED
!        635 - ONLY ONE SET-ID ALLOWED IN A SCAN CARD
!        636 - EXTRA VALUE ENCOUNTERED
!        637 - ILLEGAL COMPONENT SPECIFIED
!        638 - COMPONENT LIMIT OF 31 IS EXCEEDED
!        639 - SET ID ERROR (REQUESTED BEFORE EQUAL SIGN OR SPLITTED ID)
!        640 - TOO MANY COMPONENTS BY NAME
!        641 - -MAX EXCEEDS +MAX
!        642 - COMPONENT NAME NOT AVAILABLE FOR ELEMENT SELECTED
!        643 - SCAN BY STRESS OR FORCE ONLY
!        644 - WARNING MESSAGE FOR POSSIBLE INSUFFICIENT CORE
!        909 - CORE ARRAY NOT INITIALIZED CORRECTLY, OR MZERO IS NOT SET
!              IN AGREEMENT WITH XRCARD
 
!     EXAMPLE - TO ADD A NEW ELEMENT TO THE SCAN MODULE   BY G.C. 7/89
!           1.  INCREASE COMP DIMENSION TO ALLOW NEW COMPONENT WORDS
!               IF THEY DO NOT ALREADY EXIST.
!           2.  EXPAND THE SP-ARRAY IF NECESSARY. INCREASE NCOMP BY
!               THE NUMBER OF NEW WORDS ADDED
!           3.  REACTIVATE THE CORRESPONDING WORD IN ETAB THAT POINTS
!               TO THE NEW ARRAY IN TAB
!           4.  IF SP-ARRAY IS USED, MAKE SURE THAT THE COMPONENT WORDS
!               ARE PROPERLY PROCESSED, IN STATEMENT NOS. 110-120
!           5.  SET THE CODED WORDS IN TAB. SEE COMMENTS FUTHER DOWN
!           6.  PREPARE FORMAT FOR COMPONENT WORDS PRINT OUT (FMT 690)
!               UPDATE ISP-ARRAY IN CASE SP-ARRAY WAS USED PREVIOUSLY
 
 
 INTEGER, INTENT(IN OUT)                  :: i81
 INTEGER, INTENT(OUT)                     :: nz
 INTEGER, INTENT(OUT)                     :: j400
 LOGICAL :: debug,   bit64
 INTEGER :: core(1), scr1,    nam(2),  stress,  force,  &
     seti,    BLANK,   e,       ERR,     equal
 INTEGER :: SAVE(5), etab(90),tab(10,17),       comp(2,60),  &
     comp1(2,19),      comp2(2,19),      isp(10),  &
     tab1(10,9),       tab2(10,8),       sp(30),  &
     comp3(2,19),      comp4(2,3),       icse(400), corey(401)
 DIMENSION       rcore(1),ll(4),   cc(4),   keywds(3)
 REAL :: bcd(2,3)
 CHARACTER (LEN=23) :: ufm
 COMMON /xmssg / ufm
 COMMON /system/ ibuf,    nout,    nogo,    skip(5), nlpp, more(2), line
 COMMON /machin/ mach
 COMMON /gpta1 / nelem,   last,    incr,    e(1)
 COMMON /xifp1 / BLANK,   bit64
 COMMON /ifp1a / scr1,    casecc,  is,      nwpc,    ncpw,  &
     nmodes,  icc,     nset,    dummy(3),isub, lencc,   iblnk,   iequal,  IEOR
 COMMON /ifp1hx/ msst,    misset(1)
 COMMON /zzzzzz/ corex(1)
 EQUIVALENCE     (core(1)  ,rcore(1)   ,  corey(401)           ),  &
     (corex(1) ,corey(1)   ,  icse(1)              ),  &
     (comp(1,1),comp1(1,1)), (comp4(1,1),comp(1,58)),  &
     (tab1(1,1) ,tab(1, 1)), (comp2(1,1),comp(1,20)),  &
     (tab2(1,1),tab(1,10)) , (comp3(1,1),comp(1,39)),  &
     (BLANK    ,xblank   ) , (stress    ,bcd(1,1)  ),  &
     (force    ,bcd(1,2) ) , (shea      ,comp(1, 9)),  &
     (norm     ,comp(1,9)) , (mome      ,comp(1,22))
 DATA    ncomp,  equal,    lll,    seti,   debug    /  &
     58,     4H=   ,   4HL   , 4HSET , .false.  /
 DATA    llc,    comma,    mzero,  bcd                          /  &
     4HC   , 4H,   ,   -0    , 4HSTRE,2HSS,4HFORC,1HE,2*1H  /
 DATA    nam   / 4HIFP1,   4HH   /
 DATA    comp1 / 4HAXIA,   4HL   , 4HTORS,   4HIONA,  &
     4HRADI,   4HAL  , 4HNORM,   4HAL  ,  &
     4HPRIN,   4HCIPA, 4HMAJO,   4HR   ,  &
     4HMINO,   4HR   , 4HBEND,   4HING ,  &
     4HNORM,   4H-x  ,
!    9                or          -U  , AL-1, AL-X  &
 4HNORM,   4H-y  ,
!    O                or          -V  , AL-2, AL-Y  &
 4HNORM,   4H-z  , 4HSHEA,   4HR   ,
!    2                or          R-1Z, R-41  &
 4HSHEA,   4HR-xy,
!    3                or          R-ZR, R-X , R-U , R-12  &
 4HSHEA,   4HR-yz,
!    4                or          R-RT, R-Y , R-V , R-23  &
 4HSHEA,   4HR-zx,
!    5                or          R-ZT, R-UV, R-2Z, R-34  &
 4HMAX-,   4HSHR , 4HSHR-,   4HFORC,  &
     4HOCT-,   4HSHR , 4HSA-m,   4HAX  /
 DATA    comp2 / 4HSB-m,   4HAX  , 4HMOME,   4HNT  ,  &
     4HMOME,   4HNT-a,
!    2                or          NT-X, NT-U , NT-1  &
 4HMOME,   4HNT-b,
!    3                or          NT-Y, NT-V , NT-2  &
 4HCURV,   4H    , 4HTORQ,   4HUE  ,  &
     4HCIRC,   4HUM  , 4HTWIS,   4HT   ,  &
     4HMARG,   4HIN  , 4HMAX ,   4H    ,  &
     4HMEAN,   4H    , 4HAVG ,   4H    ,  &
     4HMEM-,   4HT   , 4HMEM-,   4HC   ,  &
     4HFLEX,   4H-t  , 4HFLEX,   4H-c  ,  &
     4HPRIN,   4HC-a , 4HPRIN,   4HC-b ,  &
     4HPRIN,   4HC-c /
 DATA    comp3 / 4HEFOR,   4HCE  , 4HFORC,   4HE-1 ,
!                     or          E-12,  &
 4HFORC,   4HE-2 ,
!                     or          E-23,  &
 4HFORC,   4HE-3 ,
!                     or          E-34,  &
 4HFORC,   4HE-4 ,
!                     or          E-41,  &
 4HKICK,   4H-for, 4HSIG-,   4HX   ,  &
     4HSIG-,   4HY   , 4HTAU-,   4HXY  ,  &
     4HHELP,   4H    , 4HON-l,   4HINE ,  &
     4HFX+f,   4HY   , 4HFXY ,   4H    ,  &
     4HMX+m,   4HY   , 4HMXY ,   4H    ,  &
     4HVX+v,   4HY   , 4HKICK,   4H on1,  &
     4HKICK,   4H on2, 4HKICK,   4H on3/
 DATA    comp4 / 4HKICK,   4H on4, 4H .. ,   4H .. ,  &
     4H .. ,   4H .. /
 DATA    sp    / 4HR-zr,   4HR-u ,   4HR-rt,   4HR-v ,   4HR-zt,  &
     4HR-uv,   4HNT-x,   4HNT-u,   4HNT-y,   4HNT-v,  &
     4H-u  ,   4H-v  ,   4HR-x ,   4HR-y ,   4HR-41,  &
     4HR-12,   4HR-23,   4HR-34,   4HNT-1,   4HNT-2,  &
     4HR-1Z,   4HR-2Z,   4HAL-1,   4HAL-2,   4HAL-x,  &
     4HAL-y,   4HE-12,   4HE-23,   4HE-34,   4HE-41/
 DATA    etab  / 1, -02,   1,   2,   2,   3,   3,   3,   4,   1,  &
     6,   6,   6, -14,   3,   4,   3,   3,   3, -20,  &
     -21, -22, -23, -24, -25, -26, -27, -28, -29, -30,  &
     -31, -32, -33,   7,   8,   9,  10,  11, -39, -40,  &
     -41, -42, -43, -44, -45, -46, -47, -48, -49, -50,  &
     -51, -52, -53, -54, -55, -56, -57, -58, -59, -60,  &
     -61,   4,   4,  15,  12,  12,  13,  14,  14, -70,  &
     -71, -72, -73, -74, -75, -76, -77, -78, -79,   5,  &
     7, -82,  15, -84, -85, -86, -87, -88, -89, -90/
 DATA    tab1  /
!    1. ROD, TUB, CONROD  &
 01000002, 02000004, 28000503,        0,        0,  &
     -01000002,-25000003,        0,        0,        0,
!    2. SHEAR, TWIST  &
 16000002, 28000004, 31000003, 29000002,        0,  &
     -40000002,-41000003,-22000002,-23000003,        0,
!    3. TRIA1, TRIA2, QUAD1, QUAD2, TRBSC, TRPLT, QDPLT  &
 09001103, 10001204, 13001305, 06001507, 07001608,  &
     16001709,-22000002,-23000003,-13000005,-14000006,
!    4. TRMEM, QDMEM, QDMEM1, DQMEM2  &
 09000002, 10000003, 13000004, 06000006, 07000007,  &
     16000008,-40000403,-41000605,-42000807,-43000902,
!   **  CONTINUE...  &
 -55000010,-56000012,-57000014,-58000016,-13000011,  &
     -14000013,-15000015,-12000017,        0,        0,
!    6. ELAS1, ELAS2, ELAS3, IS2D8  &
 18000002,-26000002,        0,        0,        0,  &
     -40000904,-41000603,-42000805,-43000702,        0,
!    7. BAR, ELBOW  &
 19000807, 20001514, 28001609, 01000006,        0,  &
     -01000008,-25000009,-12000605,-22000302,-23000504,
!    8. CONEAX  &
 09041852, 10051852, 15061852, 06081852, 07091852,  &
     16101852,-22000003,-23000004,-13000006,-14000007,
!    9. TRIARG  &
 03000002, 26000003, 01000004, 12000005,        0,  &
     -03020353,-26030353,-01040353,        0,        0/
 DATA    tab2  /
!    O. TRAPRG  &
 03020455, 26030455, 01040455, 12050455, 17060455,  &
     -03020354,-26030354,-01040354,        0,        0,
!    1. TORDRG  &
 32020553, 33030553, 34040553, 35050553, 17060553,  &
     -03000802,-26000903,-01001004,-21001105,-24001307,
!   12. IHEX1, IHEX2  &
 09032258, 13042258, 36052258, 30092258, 10112258,  &
     14122258, 37132258, 11172258, 15182258, 38192258,
!   13. IHEX3  &
 09032382, 13042382, 36052382, 30092382, 10122382,  &
     14132382, 37142382, 11182382, 15192382, 38202382,
!   14. TRIAAX, TRAPAX  &
 03030853, 01040853, 26050853, 33060853, 34070853,  &
     35080853,-03030453,-26040453,-01050453,        0,
!   15. QUAD4, TRIA3 (GENERAL)  &
 09000003, 10000004, 13000005, 06000007, 07000008,  &
     16000009,-50000302,-51000004,-52000605,-53000007, -54000908,
!   **. QUAD4, TRIA3 (LAYER), 9 DIGIT CODE  &
 81030899, 82040899, 84050899, 83070899,  &
     85080899,        0,        0,        0,        0,
!   17.  &
 10*0/
 
!     FIRST 2 DIGITS IN A TAB ITEM ARE COMPONENT POINTER, POINTING TO
!     THE BCD WORDS IN COMP ARRAY. POSITIVE FOR STRESS, AND NEGATIVE FOR
!     FORCE DATA (WITH SOME EXCEPTIONS). THIS POINTER IS USED ONLY
!     LOCALLY WITHIN THIS SUBROUTINE.
!     NEXT 3 NUMBERS (2 DIGITS EACH) ARE POINTERS TO THE FIELD NOS.
 
!     SPECIAL CASE -
!     IF LAST FIELD IS GREATER THAN 50 THEN, THIS LAST FIELD MINUS 50 IS
!     THE REPEAT FLAG. IF LAST FIELD IS 99, WE HAVE AN OPEN-END REPEAT.
!     IF LAST FIELD IS GREATER THAN 50, NEXT TO LAST FIELD IS FIELD
!     INCREMENT, AND THE FIELD IN FRONT IS THE FIRST STARTING COLUMN TO
!     BE SCANNED.
!     THE QUAD4/TRIA3 LAYER HAS COMPONENT INDICES 81 THRU 85
!     (THUS - IN FUTURE EXPANSION, ARRAY COMP SHOULD NOT EXCEED 80)
 
!     E.G.  TAB(3,1) =  09 00 11 03
!                       09          = NORMAL-X (STRESS)
!                          00       = SKIP
!                             11 03 = SCAN BY 3RD AND 11TH FIELDS
 
!     E.G.  TAB(9,8) = -01 04 03 54
!                      -01          = AXIAL (FORCE)
!                                54 = REPEAT 4 TIMES
!                             03    = INCREASE BY 3 ON EACH REPEAT
!                          04       = SCAN BY 4, 7, 10, AND 13TH FIELDS
 
 IF (j400 == 2) GO TO 400
 IF (mach == 2 .OR. mach >= 5) mzero = -1
 CALL sswtch (20,j)
 IF (j == 1) debug = .true.
 ERR   = -909
 nscan = lencc - 1
 IF (core(i81+3) /= mzero) GO TO 300
 iscan = i81
 iwds  = i81 + 1
 iisub = i81 + 2
 ielem = i81 + 3
 iset  = i81 + 4
 icomp = i81 + 5
!     +MAX  = I81 + 6 = TOP N
!     -MAX  = I81 + 7
 irept = i81 + 8
 
!     NOTE - THE IISUB WORD WILL BE DROPPED WHEN THESE WORDS ARE
!            TRANSFERRED TO CASECC
 
 jcomp = iwds
 iend  = i81 + core(i81)*2 - 1
 core(iscan) = 0
 core(iisub) = isub
 core(ielem) = 0
 core(iset ) = 0
 core(jcomp) = 0
 nwdss = 0
 nwdsf = 0
 ieq   = 0
 MAX   = 0
 MIN   = 0
 nsv   = 0
 nrp   = 0
 ii    = i81 + 3
 
 10   ii = ii + 2
 jj = core(ii  )
 kk = core(ii+1)
 jx = jj
 IF (.NOT.bit64) GO TO 15
!WKBD 3/94      CALL MVBITS (BLANK,0,32,JX,0)
!WKBD 3/94      CALL MVBITS (BLANK,0,32,KK,0)
 15   IF (jj == IEOR) GO TO 200
 IF (ii > iend) IF (jj) 150,190,20
 GO TO 30
 
!     DECODE BCD WORD
 
 20   iend = ii + jj*2 - 1
 ii   = ii - 1
 GO TO 10
 
!     LOOK FOR EQUAL SIGN OR SET
 
 30   ERR = -617
 IF (jj /= mzero) GO TO 40
 IF (kk /= equal) GO TO 300
 ieq = ieq + 1
 IF (ieq > 1) GO TO 300
 core(icomp+1) = 0
 core(icomp+2) = 0
 GO TO 10
 40   IF (jx == seti) GO TO 130
 
!     LOOK FOR STRESS OR FORCE
 
 IF (ieq == 1) GO TO 300
 IF (jx /= stress) GO TO 50
 core(iscan) = core(iscan) + 10000000
 GO TO 10
 50   IF (jx /= force) GO TO 60
 core(iscan) = core(iscan) + 20000000
 GO TO 10
 
!     LOOK FOR ELEMENT, DROP THE FIRST LETTER C IF NECESSARY
 
 60   IF (core(ielem) /= 0) GO TO 100
 jc = nam(1)
 kc = jc
 IF (khrfn2(jx,1,1) /= llc) GO TO 70
 jc = khrfn3(BLANK,jx,1,1)
 kc = khrfn3(BLANK,kk,1,1)
 jc = khrfn1(jc,4,kk,1)
 70   j  = 1
 DO  i = 1,nelem
   IF (jx == e(j) .AND. kk == e(j+1)) GO TO 90
   IF (jc == e(j) .AND. kc == e(j+1)) GO TO 90
   j = j + incr
 END DO
 GO TO 100
 90   core(ielem) = i
 nwdss = e(j+17)
 nwdsf = e(j+18)
 GO TO 10
 
!     LOOK FOR COMPONENT
 
 100  DO  i = 1,ncomp
   IF (jx == comp(1,i) .AND. kk == comp(2,i)) GO TO 120
 END DO
 ERR = -634
 i = 0
 
!     SP ARRAYS
!        1     2     3     4     5     6     7     8     9     10
!        R-ZR  R-U   R-RT  R-V   R-ZT  R-UV  NT-X  NT-U  NT-Y  NT-V
!        11    12    13    14    15    16    17    18    19    20
!        -U    -V    R-X   R-Y   R-41  R-12  R-23  R-34  NT-1  NT-2
!        21    22    23    24    25    26    27    28    29    30
!        R-1Z  R-2Z  AL-1  AL-2  AL-X  AL-Y  E-12  E-23  E-34  E-41
 
 IF (jx /=  force) GO TO 115
 IF (kk == sp(27)) i = 40
!                 FORCE-12 (USED IN QDMEM2)
 IF (kk == sp(28)) i = 41
!                 FORCE-23
 IF (kk == sp(29)) i = 42
!                 FORCE-34
 IF (kk == sp(30)) i = 43
!                 FORCE-41
 IF (i == 0) GO TO 300
 115  IF (jx /= norm .AND. jx /= shea .AND. jx /= mome) GO TO 300
 IF (kk == sp( 1) .OR. kk == sp( 2) .OR. kk == sp(13)) i = 13
!         SHEAR-ZR          SHEAR-U           SHEAR-X
 IF (kk == sp( 3) .OR. kk == sp( 4) .OR. kk == sp(14)  .OR.  &
     kk == sp(17)) i = 14
!         SHEAR-RT          SHEAR-V           SHEAR-6
!         SHEAR-23
 IF (kk == sp( 5) .OR. kk == sp( 6) .OR. kk == sp(18)) i = 15
!         SHEAR-ZT          SHEAR-UV          SHEAR-34
 IF (i /= 0) GO TO 120
 IF (kk == sp( 7) .OR. kk == sp( 8)) i = 22
!         MOMENT-X          MOMENT-U
 IF (kk == sp( 9) .OR. kk == sp(10)) i = 23
!         MOMENT-Y          MOMENT-V
 IF (kk == sp(15) .OR. kk == sp(19)) i = 12
!         SHEAR-41          MOMENT-1
 IF (kk == sp(11) .OR. kk == sp(25)) i = 9
!         NORM-U            NORNAL-X
 IF (kk == sp(12) .OR. kk == sp(26)) i = 10
!         NORM-V            NORNAL-Y
 IF (i /= 0) GO TO 120
 IF (kk == sp(16) .OR. kk == sp(20)) i = 13
!         SHEAR-12          MOMENT-2
 
!     SECOND SET KEYWORDS FOR QUAD4/TRIA3 LAYER, 81 AND HIGHER
!     (THE GENERAL QUAD4/TRIA3 KEYWORDS ARE BELOW 80)
 
 IF (i /= 0) GO TO 120
 IF (kk == sp(23)) i = 81
!         NORAML-1
 IF (kk == sp(24)) i = 82
!         NORMAL-2
 IF (kk == sp(16)) i = 84
!         SHEAR-12
 IF (kk == sp(21)) i = 83
!         SHEAR-1Z
 IF (kk == sp(22)) i = 85
!         SHEAR-2Z
 IF (i == 0) GO TO 300
 
 120  IF (i == 48) GO TO 320
 IF (i == 49) GO TO 900
 ERR = -640
 IF (nsv > 5) GO TO 300
 IF (nsv <= 0) GO TO 125
 DO  j = 1,nsv
   IF (SAVE(j) == i) GO TO 127
 END DO
 125  nsv = nsv + 1
 SAVE(nsv) = i
 
!     TWO WORDS, PRINCIPAL AND TORSIONAL, HAVE A LETTER L TOO LONG
 
!     LLL IS 4HL   , BLANK FILL
!     IBLNK IS 1H  , ZERO  FILL
 
 127  iword = core(ii+2)
!WKBD 3/94      IF (BIT64) CALL MVBITS (BLANK,0,32,IWORD,0)
 IF (ii < iend .AND. iword == lll .AND. core(ii+3) == iblnk) ii = ii + 2
 GO TO 10
 
!     PROCESS SET
 
 130  ERR = -635
 IF (core(iset) /= 0) GO TO 300
 ERR = -639
 IF (ieq /= 1  .OR. core(ii+2) /= -1) GO TO 300
 core(iset) = core(ii+3)
 ii = ii + 2
 j  = nwpc + 1 + icse(lencc)
 DO  i = 1,nset
   IF (core(iset) == core(j)) GO TO 10
   j = j + core(j+1) + 3
 END DO
 ERR  = -608
 msst = msst + 1
 IF (msst >  0) GO TO 300
 misset(msst) = core(iset)
 ERR = 0
 GO TO 10
 
!     NUMERIC DATA
 
 150  IF (jj == -2) GO TO 170
 IF (ieq == 1) GO TO 160
 
!     INTEGER BEFORE EQUAL SIGN = COMPONENT(S)
 
 ERR = -637
 IF (jj /= -1 .OR. kk <= 1) GO TO 300
 ERR = -638
 IF (kk > 31) GO TO 300
 core(jcomp) = core(jcomp) + 2**(kk-1)
 GO TO 10
 
!     INTEGER AFTER EQUAL SIGN = TOP N
 
 160  ERR = -608
 IF (ieq /= 1) GO TO 300
 ERR = -636
 IF (MAX-1) 180,300,300
 
!     F.P. DATA = +MAX OR -MAX
 
 170  ERR = -608
 IF (ieq /= 1) GO TO 300
 MIN = 1
 ERR = -636
 IF (MAX >= 2) GO TO 300
 180  MAX = MAX + 1
 core(icomp+MAX) = kk
 ERR = -641
 IF (MAX == 2 .AND. rcore(icomp+2) > rcore(icomp+1)) GO TO 300
 GO TO 10
 
!     READ CONTINUATION CARD
 
 190  CALL READ (*290,*290,scr1,core(1),nwpc,0,flag)
 WRITE (nout,310) icc,(core(i),i=1,nwpc)
 icc  = icc  + 1
 line = line + 1
 IF (line > nlpp) CALL page
 ii = i81 + 8
 nz = nz - ii
 CALL xrcard (core(ii),nz,core(1))
 ii = ii - 2
 iend = ii
 GO TO 10
 
!     SCAN CARD COMPLETED
 
 200  IF (nogo /= 0) GO TO 330
 ERR = -643
 IF (core(iscan) /= 10000000 .AND. core(iscan) /= 20000000) GO TO 300
 core(icomp) = core(jcomp)
 core(iwds ) = 6
 IF (core(iset ) == 0) core(iset) = -1
 IF (core(ielem) == 0 .AND. core(icomp+2) == 0) core(ielem) = -1
 IF (MAX == 0) core(icomp+1) = 20
 IF (MAX <= 1 .AND. MIN == 0) core(icomp+2) = -1
 IF (core(icomp+2) /= -1) GO TO 205
 
!     COMPUTE HOW HIGH TOPN CAN GO ASSUMING LINK14 HAS AN OPEN CORE SIZE
!     SAME AS THAT OF LINK1
 
 IF (core(iscan) >= 20000000) nwdss = nwdsf
 IF (2*nwdss*core(icomp+1) > korsz(icse(1))) CALL ifp1d (644)
 
!     CONVERT NAMED COMPONENT TO FIELD NO.
 
 205  IF (nsv == 0 .OR. core(ielem) == -1) GO TO 250
 i = core(ielem)
 i = etab(i)
 ERR = -642
 IF (i < 0) GO TO 300
 ie = 10
 IF (i == 14) ie = 16
 loop240:  DO  k = 1,nsv
   IF (core(iscan) == 20000000) SAVE(k) = -SAVE(k)
   DO  j = 1,ie
     ii = tab(j,i)/1000000
     IF (ii == 0) CYCLE
     IF (SAVE(k) == ii) GO TO 220
   END DO
   
!     5 SPECIAL CASES WHERE TAB ARRAY OF 10 IS NOT LONG ENOUGH
!     SET THE 11TH ITEM OF ECAH OF THESE 3 CASES
   
   ii = 0
   IF (i == 3 .AND. SAVE(k) == 27) ii = -27000004
!         TRIA1         TWIST (MOMENT)
   IF (i == 12 .AND. SAVE(k) == 18) ii = +18102258
!         IHEX1         OCT-SHR (STRESS)    +18102270 (IHEX2)
   IF (i == 13 .AND. SAVE(k) == 18) ii = +18102382
!         IHEX3         OCT-SHR (STRESS)
   IF (i == 12 .AND. SAVE(k) == 16) ii = +16102270
!         IHEX2         MAX-SHR (STRESS)
   IF (i == 13 .AND. SAVE(k) == 16) ii = +16102382
!         IHEX2         MAX-SHR (STRESS)
   ERR = -637
   IF (ii > 0) THEN
     GO TO   225
   ELSE
     GO TO   300
   END IF
   220  ii = IABS(tab(j,i))
   225  ii = MOD(ii,1000000)
   DO  jj = 1,3
     kk = MOD(ii,100)
     IF (core(icomp) /= 1) GO TO 223
     core(icomp) = 0
     IF (ielem == 66) kk = 70
!                    IHEX2
     IF (ielem == 69) kk = 54
!                    TRIATS
     nrp = nrp + kk
     kk  = 0
     223  IF (kk <= 50) GO TO 227
     nrp = (kk-50)*100
!     NRP = 4900 FOR OPEN-END REPEAT FLAG
     kk  = 1
     227  IF (kk > 0) core(icomp) = core(icomp) + 2**(kk-1)
     ii  = ii/100
     IF (ii == 0) CYCLE loop240
   END DO
 END DO loop240
 
!     NRP/100 IS REPEAT FLAG, AND MOD(NRP,100) IS INCREMENT
 
 250  IF (nogo /= 0) GO TO 330
 core(irept) = nrp
 
!     FINAL ERROR CHECK
 
 IF (core(ielem) == 0 .OR. core(icomp) == 0) j400 = 1
 ERR = -617
 j = 0
 DO  i = 1,8
   IF (core(i81+j) == 0) GO TO 300
   j = j + 1
 END DO
 
!     ALL GO WELL, RE-SET PARAMETERS
!     NOTE - THE (LENCC-1) WORD OF CASECC RECORDS THE NO. OF SCAN CARDS
 
 nset = nset + 1
 ii = (isub-1)*lencc
 icse(nscan+ii) = icse(nscan+ii) + 1
!WKBD IF (DEBUG) CALL BUG1 ('IFP1H',270,CORE(I81),9)
 i81 = i81 + 9
 
!     TURN ON STRESS OR FORCE OUTPUT FLAGS IF THEY ARE NOT ALREADY DONE
!     BY THE USER.   SET OUTPUT OPTIONS TO - ALL, NOPRINT, AND REAL
!     (WORD 23 ON CASECC IS STRESS OUTPUT FLAG, AND
!      WORD 26 ON CASECC IS FORCE  OUTPUT FLAG)
 
 IF (core(iscan) == 20000000 .OR. icse(23+ii+1) /= 0) GO TO 280
 icse(23+ii  ) =-1
 icse(23+ii+1) = 2
 icse(23+ii+2) = 1
 280  IF (core(iscan) /= 20000000 .OR. icse(26+ii+1) /= 0) GO TO 330
 icse(26+ii  ) =-1
 icse(26+ii+1) = 2
 icse(26+ii+2) = 1
 GO TO 330
 
 290  CALL mesage (-1,scr1,nam)
 300  CALL ifp1d (ERR)
 310  FORMAT (11X,i8,6X,20A4)
 IF (icse(nscan) < 0) GO TO 330
 IF (ERR /= -634 .AND. ERR /= -637 .AND. ERR /= -642) GO TO 330
 icse(nscan) = -10000
 320  j400 = 1
 330  RETURN
 
!     PRINT OUT SCAN COMPONENT KEYWORDS
 
 400  CALL page1
 ii = 20
 GO TO 810
 410  ii = 0
 WRITE  (nout,420)
 420  FORMAT (46H0*** component keywords for the scan operation, //5X,  &
     59HFORCE/stress    keyword        component (output field no.),  &
     /5X,15(4H----),/)
 lline = 15
 425  FORMAT (/5X,17HROD, tube, conrod,/)
 GO TO 700
 430  isp(8) = 19
 isp(9) = 20
 lline  = 11
 435  FORMAT (/5X,12HSHEAR, twist,/)
 GO TO 700
 440  isp( 7) = 7
 isp( 8) = 9
 isp( 9) = 13
 isp(10) = 14
 lline   = 14
 445  FORMAT (/5X,47HTRIA1, tria2, quad1, quad2, trbsc, trplt, qdplt,/)
 GO TO 700
 450  WRITE  (nout,455)
 455  FORMAT (10X,'FORCE      TWIST',15X,'4')
 line = line + 1
 isp( 7) = 27
 isp( 8) = 28
 isp( 9) = 29
 isp(10) = 30
 lline = 13
 460  FORMAT (/5X,28HTRMEM, qdmem, qdmem1, qdmem2,/)
 GO TO 700
 470  lline = 8
 GO TO 700
 480  lline = 9
 490  FORMAT (/5X,26HELAS1, elas2, elas3, is2d8,/)
 GO TO 700
 500  lline = 12
 510  FORMAT (/5X,10HBAR, elbow,/)
 GO TO 700
 520  isp(1) = 11
 isp(2) = 12
 isp(3) = 6
 isp(7) = 8
 isp(8) = 10
 lline  = 13
 530  FORMAT (/5X, 6HCONEAX,/)
 GO TO 700
 540  lline = 10
 550  FORMAT (/5X, 6HTRIARG,/)
 GO TO 700
 560  lline = 11
 570  FORMAT (/5X, 6HTRAPRG,/)
 GO TO 700
 580  lline = 13
 590  FORMAT (/5X, 6HTORDRG,/)
 GO TO 700
 600  lline = 14
 610  FORMAT (/5X,12HIHEX1, ihex2,/)
 GO TO 700
 620  WRITE  (nout,625)
 625  FORMAT (10X,'STRESS     MAX-SHR',12X,'10, 32, 54, 76 ... ETC',  &
     /10X,'STRESS     OCT-SHR',12X,'10, 32, 54, 76 ... ETC')
 line  = line + 1
 lline = 14
 630  FORMAT (/5X, 6HIHEX3 ,/)
 GO TO 700
 640  WRITE  (nout,645)
 645  FORMAT (10X,'STRESS     MAX-SHR',12X,'10, 33, 56, 79 ... 746',  &
     /10X,'STRESS     OCT-SHR',12X,'10, 33, 56, 79 ... 746')
 line  = line + 1
 lline = 12
 650  FORMAT (/5X,14HTRIAAX, trapax,/)
 GO TO 700
 660  isp(1) = 25
 isp(2) = 26
 lline  = 19
 665  FORMAT (/5X,12HQUAD4, tria3,/)
 GO TO 700
 
!   . QUAD4/TRIA3 LAYER
 
 670  isp(2) = 23
 isp(3) = 24
 isp(4) = 16
 isp(5) = 21
 isp(6) = 22
 lline  = 5
 GO TO 700
 
 680  WRITE  (nout,685)
 685  FORMAT (1X)
 GO TO 840
 
 700  ii = ii + 1
 CALL page2 (lline)
 SELECT CASE ( ii )
   CASE (    1)
     GO TO 701
   CASE (    2)
     GO TO 702
   CASE (    3)
     GO TO 703
   CASE (    4)
     GO TO 704
   CASE (    5)
     GO TO 720
   CASE (    6)
     GO TO 706
   CASE (    7)
     GO TO 707
   CASE (    8)
     GO TO 708
   CASE (    9)
     GO TO 709
   CASE (   10)
     GO TO 710
   CASE (   11)
     GO TO 711
   CASE (   12)
     GO TO 712
   CASE (   13)
     GO TO 713
   CASE (   14)
     GO TO 714
   CASE (   15)
     GO TO 715
   CASE (   16)
     GO TO 720
   CASE (   17)
     GO TO 717
   CASE (   18)
     GO TO 700
   CASE (   19)
     GO TO 700
   CASE (   20)
     GO TO 720
 END SELECT
 701  WRITE (nout,425)
 GO TO 720
 702  WRITE (nout,435)
 GO TO 720
 703  WRITE (nout,445)
 GO TO 720
 704  WRITE (nout,460)
 GO TO 720
 706  WRITE (nout,490)
 GO TO 720
 707  WRITE (nout,510)
 GO TO 720
 708  WRITE (nout,530)
 GO TO 720
 709  WRITE (nout,550)
 GO TO 720
 710  WRITE (nout,570)
 GO TO 720
 711  WRITE (nout,590)
 GO TO 720
 712  WRITE (nout,610)
 GO TO 720
 713  WRITE (nout,630)
 GO TO 720
 714  WRITE (nout,650)
 GO TO 720
 715  WRITE (nout,665)
 GO TO 720
 717  lline = 0
 GO TO 700
 
 720  DO  i = 1,10
   jj = tab(i,ii)
   IF (jj == 0) CYCLE
   bcd(1,3) = bcd(1,1)
   bcd(2,3) = bcd(2,1)
   IF (jj >= 0) GO TO 725
   bcd(1,3) = bcd(1,2)
   bcd(2,3) = bcd(2,2)
   725  jj = IABS(jj)
   
!                     +-------------------+    LL(1) = XX
!     JJ = TAB(I,J)=  ! CC ! ZZ ! YY ! XX !    LL(2) = YY
!                     +-------------------+    LL(3) = ZZ
!                                              LL(4) = CC
   DO  j = 1,4
     ll(j) = MOD(jj,100)
     jj = jj/100
   END DO
   jj = ll(4)
   
!     QUAD4/TRIA3 LAYER IF JJ IS 81 THRU 85
   
   IF (jj == 81 .OR.  jj == 82) jj = 9
   IF (jj >= 83 .AND. jj <= 85) jj = 12
   
   keywds(1) = comp(1,jj)
   keywds(2) = comp(2,jj)
   keywds(3) = BLANK
   IF (ii == 4 .OR. ii == 16) GO TO  735
   IF (jj == 2 .OR. jj == 5) keywds(3) = lll
   IF (jj < 9 .OR. jj > 30) GO TO 740
   IF (jj == 11 .OR. (jj >= 16 .AND. jj <= 21)) GO TO 740
   735  j = isp(i)
   IF (j > 0) keywds(2) = sp(j)
   740  IF (ll(1) > 50) GO TO 745
   ll(4) = 0
   idupl = 0
   jj    = 3
   GO TO 760
   745  idupl = ll(1) - 50
   inc   = ll(2)
   jj = MIN0(idupl,4)
   kk = ll(3)
   DO  j = 1,jj
     ll(j) = kk
     kk = kk+inc
   END DO
   kk = inc*idupl + ll(1)
   760  DO  j = 1,jj
     IF (ll(j) == 0) GO TO 770
     cc(j) = comma
   END DO
   j  = jj + 1
   770  jj = j  - 1
   cc(jj) = xblank
   WRITE  (nout,775) bcd(1,3),bcd(2,3),keywds,(ll(j),cc(j),j=1,jj)
   775  FORMAT (10X,a4,a2,5X,2A4,a1,9X,4(i3,a1))
   IF (idupl <= 4) CYCLE
   IF (ii /= 12 .AND. ii /= 14 .AND. ii /= 16) WRITE (nout,780) kk
   IF (ii == 12 .OR.  ii == 14 .OR.  ii == 16) WRITE (nout,785)
   780  FORMAT (1H+,54X,3H...,i4)
   785  FORMAT (1H+,54X,3H...,5H etc.)
 END DO
 810  DO  j = 1,10
   isp(j) = 0
 END DO
 SELECT CASE ( ii )
   CASE (    1)
     GO TO 430
   CASE (    2)
     GO TO 440
   CASE (    3)
     GO TO 450
   CASE (    4)
     GO TO 470
   CASE (    5)
     GO TO 480
   CASE (    6)
     GO TO 500
   CASE (    7)
     GO TO 520
   CASE (    8)
     GO TO 540
   CASE (    9)
     GO TO 560
   CASE (   10)
     GO TO 580
   CASE (   11)
     GO TO 600
   CASE (   12)
     GO TO 620
   CASE (   13)
     GO TO 640
   CASE (   14)
     GO TO 660
   CASE (   15)
     GO TO 670
   CASE (   16)
     GO TO 680
   CASE (   17)
     GO TO 840
   CASE (   18)
     GO TO 840
   CASE (   19)
     GO TO 840
   CASE (   20)
     GO TO 410
 END SELECT
 840  WRITE  (nout,850)
 850  FORMAT (//5X,'USE OUTPUT FIELD NUMBER(S) TO SPECIFY COMPONENT(S)',  &
     'FOR ELEMENTS OR KEYWORDS', /5X,'NOT LISTED ABOVE',/)
 RETURN
 
!     ON-LINE
 
 900  WRITE  (nout,910) ufm
 910  FORMAT (a23,', SCAN ON-LINE OPTION IS NOT AVAILABLE IN THIS ',  &
     'NASTRAN RELEASE')
 nogo = 1
 RETURN
END SUBROUTINE ifp1h
