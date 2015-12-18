SUBROUTINE onlins (*,lx)
     
!     ON-LINE SCAN ROUTINE, CALLED ONLY BY SCAN
 
!     WRITTEN FEBY G.CHAN/SPERRY,  FEB. 1986
 
 INTEGER, INTENT(IN OUT)                  :: lx
 EXTERNAL        lshift,   rshift,   andf,     orf,      complf
 LOGICAL :: debug
 INTEGER :: NAME(2),  card(20), iz(1)
 REAL :: r(2),     z(30)
 COMMON /machin/ mach
 COMMON /BLANK / ielt(2),  icomp,    ntop,     amax,     amin,  &
     ibeg,     iend,     icompx
 COMMON /system/ ibuf,     outtap,   nogo,     in,       dum(74),  &
     swtch1,   jdum(6),  intra
 COMMON /xscanx/ skip(2),  lcore,    lbeg,     lend,     dumm(2),  &
     iel,      iopt,     iset,     isort,    idum(4), debug
 COMMON /ifp1a / scr1,     casecc,   is,       nwpc,     ncpw,  &
     nmodes,   icc,      nset,     dummy(3), isub,  &
     lencc,    iblnk,    iequal,   IEOR
 
!            /ZZIFP1/ IS THE OPEN CORE FOR SCAN
 COMMON /zzzzzz/ lcse(400),core(1)
 EQUIVALENCE     (iz(1),lcse(1))
 EQUIVALENCE     (imax,amax),        (imin,amin),  &
     (idupl,ibeg),       (inc,iend), (card(1),core(1)),  (iz(1),z(1))
 DATA    BLANK , equal ,   STOP  ,   all   ,   NAME             /  &
     4H    , 4H=   ,   4HSTOP,   4HALL ,   4HONLI,   4HNS   /
 DATA    lu    , debug1,   debug2,   debug3,   i0               /  &
     1     , 4HDEBU,   4HG on,   4HG of,   0                /
 
!     INITIALIZE /IFP1A/
 
 scr1   = 301
 casecc = 101
 is     = 0
 nwpc   = 20
 ncpw   = 4
 nmodes = 0
 icc    = 0
 isub   = 1
 iblnk  = BLANK
 iequal = equal
 IEOR   = complf(0)
 IEOR   = rshift(IEOR,1)
 
!     SET INTERACTIVE FLAG TO POSITIVE, A SIGNAL TO SCAN, TOTAPE, IFP1C
 
 intra  = IABS(intra)
 IF (intra == 0) intra = 1
 
 icomp = lx
 nwpc1 = nwpc + 1
 nout  = outtap
 WRITE  (nout,10)
 10   FORMAT (///1X,'*** SCAN INTERACTIVE INPUT ***')
 
!     READ CASECC FILE AND SAVE DATA IN LCSE, ONE SUBCASE AT A TIME
!     SAVE SET DATA IN CORE BEGIN AT CORE(BGN)
 
 15   lcse(166) = 200
 lcse(199) = 0
 lcse(200) = 0
 nz    = korsz(core(1)) - 3*ibuf - 1
 nz    = MIN0(nz,lcore)
 iscan = 0
 nset  = 0
 i81   = nwpc1
 subid =-1
 lx    = 0
 IF (icomp == -2) GO TO 30
 
!     NO QUESTION ASKED IF SORT2 DATA TYPE IS USED.
 
 lx = 1
 20   WRITE  (nout,25)
 25   FORMAT (//,' ENTER SUBCASE ID (DEFAULT=FIRST SUBCASE)')
 READ   (in,26) r
 26   FORMAT (2A4)
 CALL a82int (*20,r,8,subid,i)
 IF (subid ==  0) subid = -1
 IF (intra > 10) WRITE (lu,27) subid
 27   FORMAT (///3X,'SUBCASE ID',i8)
 30   jj = 1
 CALL REWIND (casecc)
 CALL fwdrec (*110,casecc)
 32   jj = jj + 1
 CALL READ (*110,*110,casecc,lcse(jj),1,0,i)
 IF (subid == -1) subid = lcse(jj)
 IF (lcse(jj) == subid) GO TO 35
 CALL fwdrec (*110,casecc)
 GO TO 32
 35   lcse(1) = lcse(jj)
 CALL READ (*110,*125,casecc,lcse(2),199,0,i)
 lencc = lcse(166)
 lsem  = lcse(lencc)
 nset  = lcse(lencc-1)
 IF (lsem > 0) CALL READ (*110,*125,casecc,core(i81),lsem,0,i)
 i81   = i81 + lsem
 bgn   = i81
END   = i81
37   CALL READ (*40,*40,casecc,core(i81),2,0,i)
jmp   = core(i81+1)
core(i81+2) = jj
i81  = i81 + 3
CALL READ (*110,*125,casecc,core(i81),jmp,0,i)
nset = nset + 1
i81  = i81 + jmp
GO TO 37

!     SET CARD

40   WRITE  (nout,43)
43   FORMAT (//,' ENTER A BLANK, OR A SET CARD (SEE USER MANUAL P. ',  &
    '2.3-44)', /,' E.G.  SET 101 = 1, 5 THRU 20')
45   core(i81) = IEOR
nogo  = 0
CALL xread (*40,card)
IF (card(1) == BLANK .AND. card(2) == BLANK) GO TO 60
WRITE (lu,77) card
IF (card(1) /= debug1) GO TO 46
j  = lshift(1,20)
IF (card(2) == debug2) swtch1 = orf(j,swtch1)
j  = complf(j)
IF (card(2) == debug3) swtch1 = andf(j,swtch1)
debug = .false.
IF (card(2) == debug2) debug = .true.
GO TO 40
46   ib  = i81
nzz = nz - i81
CALL xrcard (core(i81),nzz,card(1))
IF (core(i81+8) /= all) GO TO 47
core(i81  ) = core(i81+4)
core(i81+1) = 1
core(i81+2) = jj
core(i81+3) =-1
i81 = i81 + 4
GO TO 50
47   icc = 1
CALL ifp1c (i81,nzz)

!     CONTINUATION CARDS FOR SET ARE READ IN BY IFP1C

IF (nogo == 0) GO TO 50
i81  = ib
GO TO 40
50   nset = nset + 1
WRITE  (nout,52) core(ib)
52   FORMAT (/,' THIS NEW SET',i6,' IS DEFINED FOR LOCAL USE ONLY',  &
    //,' ENTER A BLANK, OR ANOTHER SET CARD')
kk = 55
IF (debug) WRITE (6,55) kk,i81
55   FORMAT ('   ONLINS/',i2,4X,'I81 =',i7)
GO TO 45

!     SET DATA - FROM CORE(BGN) THRU CORE(END)

60   END = i81 - 1
nzz = nz - i81

!     SCAN CARD

70   WRITE  (nout,72)
72   FORMAT (//,' ENTER A BLANK, OR A SCAN CARD (SEE USER MANUAL P.2.3-  &
    41A',   /,'  E.G. SCAN (STRESS,CBAR,AXIAL,SA/MAX) = 15, SET 102',  &
    /,'       SCAN (FORCE,3,ROD,2) = +2000.,-1500.', /,'       SCAN (HELP)' )

75   jumph = 0
CALL xread (*70,card)
IF (card(1) == STOP  .AND. card(2) == BLANK) GO TO 135
IF (card(1) == BLANK .AND. card(2) == BLANK) GO TO 90
WRITE  (lu,77) card
77   FORMAT (20A4)
ib = i81
CALL xrcard (core(i81),nzz,card(1))
CALL ifp1h (i81,nzz,jumph)
IF (nogo  /= 0) GO TO 80
IF (jumph == 0) GO TO 82
CALL ifp1h (0,0,2)
80   i81 = ib
IF (nogo == 0) THEN
  GO TO    75
ELSE
  GO TO    70
END IF

82   j = core(ib)
IF (iscan == 0) iscan = j
IF (iscan == j) iscan = 30000000
WRITE  (nout,85)
85   FORMAT (/,' ENTER A BLANK, OR ANOTHER SCAN CARD')
kk = 87
IF (debug) WRITE (6,55) kk,i81
GO TO 75

!     MOVE SET AND SCAN DATA TO THE END OF CASECC ARRAY IN /ZZIFP1/
!     THEN, MOVE THE ENTIRE CASECC DATA (SET AND SCAN INCLUDED) TO
!     THE END OF THE OPEN CORE. FINALLY, MOVE THE SAME DATA BLOCK
!     TO THE BEGINNING OF THE OPEN CORE SPACE IN /ZZSCAN/ FOR SCAN
!     OPERATION

90   l   = lencc
IF (i81 <= nwpc1) GO TO 100
j   = bgn + 2
i81 = i81 - 1
DO  i = nwpc1,i81
  IF (i /= j) GO TO 92
  j   = j + core(j-1) + 3
  CYCLE
  92   l   = l + 1
  lcse(l) = core(i)
END DO
j   = lcore
DO  i = 1,l
  lcse(j) = lcse(i)
  j   = j - 1
END DO
IF (i > j) CALL mesage (+8,0,NAME)
j   = lcore
DO  i = 1,l
  z(i) = lcse(j)
  j   = j - 1
END DO
IF (debug) WRITE (6,99) (z(i),i=1,l)
99   FORMAT (//,' Z(1...200+) =', (/4X,10I7))
100  IF (lx > 0) lx = l

IF (iscan == 20000000) GO TO 103
IF (z(25) == 0) GO TO 140

!     STRESS SCAN

z(24) =-1
z(25) = 1
z(26) = 1
103  IF (iscan /= 20000000) GO TO 105
IF (z(28) == 0) GO TO 150

!     FORCE SCAN

z(27) =-1
z(28) = 1
z(29) = 1
105  IF (intra > 10) outtap = lu
RETURN

110  jj = jj - 1
WRITE  (nout,115) subid,(z(i),i=1,jj)
115  FORMAT (//,' SUBCASE',i5,' NOT FOUND',  &
    //,' EXISTING SUBCASES ARE -', (/5X,10I7))
GO TO 15

125  CALL mesage (+2,casecc,NAME)
GO TO 105
135  RETURN 1

140  WRITE  (nout,145)
145  FORMAT (//,' STRESS OUTPUT FILE NOT AVAILABLE FOR SCAN',//)
GO TO 75
150  WRITE  (nout,155)
155  FORMAT (//,' FORCE  OUTPUT FILE NOT AVAILABLE FOR SCAN',//)
GO TO 75
END SUBROUTINE onlins
