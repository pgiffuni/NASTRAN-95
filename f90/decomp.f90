SUBROUTINE decomp (*,ix,x,dx)
     
!     DECOMP WILL DECOMPOSE A REAL UNSYMETRIC MATRIX INTO A UNIT LOWER
!     TRIANGULAR MATRIX AND AN UPPER TRIANGULAR MATRIX,USING PARTIAL
!     PIVOTING WITHIN THE LOWER BAND
 
!     DEFINITION OF INPUT PARAMETERS
 
!     FILEA    =  MATRIX CONTROL BLOCK FOR THE INPUT  MATRIX A
!     FILEL    =  MATRIX CONTROL BLOCK FOR THE OUTPUT MATRIX L
!     FILEU    =  MATRIX CONTROL BLOCK FOR THE OUTPUT MATRIX U
!     SR1FIL   =  SCRATCH FILE
!     SR2FIL   =  SCRATCH FILE
!     SR3FIL   =  SCRATCH FILE
!     NX       =  NUMBER OF CELLS OF CORE AVAILABLE AT IX
!     DET      =  CELL WHERE THE DETERMINATE OF A WILL BE STORED
!     POWER    =  SCALE FACTOR TO BE APPLIED TO THE DETERMINATE
!                 (DETERMINATE = DET*10**POWER)
!     MINDIA   =  CELL WHERE THE VALUE OF THE MINIMUM DIAGONAL WILL BE
!                 SAVED
!     IX       =  BLOCK OF CORE AVAILABLE AS WORKING STORAGE TO DECOMP
!     X        =  SAME BLOCK AS IX, BUT TYPED REAL
!     DX       =  SAME BLOCK AS IX, BUT TYPED DOUBLE PRECISION
 
 
 INTEGER, INTENT(OUT)                     :: ix(1)
 REAL, INTENT(OUT)                        :: x(1)
 DOUBLE PRECISION, INTENT(OUT)            :: dx(1)
 INTEGER :: filea     ,filel    ,fileu    ,power    ,  &
     sysbuf    ,forma    ,typea    ,rdp      ,  &
     typel     ,eol      ,parm(5)  ,bufa     ,  &
     outbuf    ,sr1buf   ,sr2buf   ,sr3buf   ,  &
     b         ,bbar     ,c        ,cbar     ,  &
     bbar1     ,r        ,ccount   ,cbcnt    ,  &
scrflg    ,END      ,bbbar    ,bbbar1   ,  &
    count     ,sr2fl    ,sr3fl    ,sr1fil   ,  &
    sr2fil    ,sr3fil   ,sqr      ,sym      , flag      ,itran(4)
DOUBLE PRECISION :: dz        ,da       ,det      ,MAX      , mindia    , dtrn

CHARACTER (LEN=29) :: uim
CHARACTER (LEN=25) :: uwm
CHARACTER (LEN=23) :: ufm
COMMON   /xmssg /  ufm       ,uwm      ,uim
COMMON   /dcompx/  filea(7)  ,filel(7) ,fileu(7) ,sr1fil   ,  &
    sr2fil    ,sr3fil   ,det      ,power    ,  &
    nx        ,mindia   ,b        ,bbar     , c         ,cbar     ,r
COMMON   /system/  sysbuf    ,nout
COMMON   /names /  rd        ,rdrew    ,wrt      ,wrtrew   ,  &
    rew       ,norew    ,eofnrw   ,rsp      ,  &
    rdp       ,csp      ,cdp      ,sqr      ,  &
    rect      ,diag     ,lowtri   ,uprtri   , sym       ,row      ,ident
COMMON   /zntpkx/  a(4)      ,ii       ,eol
COMMON   /zblpkx/  z(4)      ,jj
COMMON   /unpakx/  itypex    ,ixy      ,jxy      ,incrx
COMMON   /packx /  itype1    ,itype2   ,iy       ,jy       , incry
EQUIVALENCE        (da,a(1))           ,(dz,z(1))          ,  &
    (forma,filea(4))    ,(typea,filea(5))   ,  &
    (ncol,filea(3))     ,(typel,filel(5))
EQUIVALENCE        (itran(1),itrn)     ,(itran(2),jtrn)    , (itran(3),dtrn)
DATA      parm(3), parm(4)/ 4HDECO,4HMP   /
DATA      ibegn  / 4HBEGN /, iend /4HEND  /

!     AT LAST, THE START OF THE PROGRAM

IF ((forma /= sqr .AND. forma /= sym) .OR. typea > rdp) GO TO 1660

!     BUFFER ALLOCATION

bufa   = nx     - sysbuf
ibufl  = bufa   - sysbuf
outbuf = ibufl  - sysbuf
sr1buf = outbuf - sysbuf
sr2buf = sr1buf - sysbuf
sr3buf = sr2buf - sysbuf
icrq   =-sr3buf
IF (icrq > 0) GO TO 1668
det    = 1.d0
power  = 0
mindia = 1.d+25
iterm  = 0
IF (filea(1) < 0) iterm = 1
filea(1) = IABS(filea(1))

!     WRITE THE HEADER RECORD ON THE OUTPUT TAPES AND INITIALIZE THE
!     TRAILER RECORDS.

CALL gopen (filel,ix(ibufl),wrtrew)
parm(2) = sr2fil
CALL OPEN  (*1670,sr2fil,ix(outbuf),wrtrew)
CALL fname (fileu(1),x(1))
CALL WRITE (sr2fil,x(1),2,1)
filel(3) = ncol
filel(4) = 4
filel(2) = 0
filel(6) = 0
filel(7) = 0
fileu(2) = 0
fileu(3) = ncol
fileu(4) = 5
fileu(6) = 0
fileu(7) = 0
filea(5) = 2
IF (ncol > 2 ) GO TO 10
imhere = 9
CALL onetwo (*1710,ix(1),x(1),dx(1),iterm)

!     CALL GENVEC TO PICK B,BBAR,C,CBAR, AND R

RETURN
10 IF (b > 0 .AND. bbar > 0) GO TO 15
imhere = 10
CALL genvec (*1710,ix(bufa),filea(1),nx,ix(1),ncol,b,bbar,c,cbar, r,1)
15 CONTINUE
bbar1  = bbar + 1
bbbar  = MIN0(b+bbar,ncol)
bbbar1 = bbbar - 1
scrflg = 0
IF (r < bbbar1) scrflg = 1
IF (scrflg == 0) GO TO 20
icrq = (bbbar1-r)*2*bbar
CALL page2 (3)
WRITE  (nout,2000) uim,icrq
2000 FORMAT (a29,' 2177, SPILL WILL OCCUR IN UNSYMMETRIC DECOMPOSITION'  &
    ,      /,i10,' ADDITIONAL MEMORY WORDS NEEDED TO STAY IN CORE.')

!     INITIALIZE POINTERS TO SPECIFIC AREAS OF CORE

20 i1   = 1
i1sp = (i1+bbar*r)*2 - 1
ipak = i1 + bbar*r + bbbar/2 + 1
i2   = ipak
i3sp = (i2 + MIN0(ncol,bbbar+bbar))*2 - 1
i3   = i2  + MIN0(ncol,bbbar+bbar) + c
i4sp = i3sp + (bbar+2)*c*2
i4   = i3 + bbar1*c + cbar
i5   = i4 + bbbar*cbar
i6sp = (i5+c*cbar)*2 - 1
i7sp = i6sp + cbar
END  = i7sp + c
parm(5) = ibegn
CALL conmsg (parm(3),3,0)

!     DEFINITION OF KEY PROGRAM PARAMETERS

!     I1     =  POINTER TO AREA WHERE COMPLETED COLUMNS OF L ARE STORED
!     I1SP   =  POINTER TO AREA WHERE THE PERMUTATION INDEXES ARE STORED
!     IPAK   =  POINTER TO AREA WHERE COLUMNS WILL BE PACKED FROM
!     I2     =  POINTER TO AREA WHERE THE NEXT COLUMN OF A IS STORED
!     I3     =  POINTER TO AREA WHERE ACTIVE COLUMNS ARE STORED
!     I4     =  POINTER TO AREA WHERE ACTIVE ROWS ARE STORED
!     I5     =  POINTER TO AREA WHERE INTERACTION ELEMENTS ARE STORED
!     I6SP   =  POINTER TO AREA WHERE SEQUENCED ACTIVE ROW INDICES
!               ARE STORED
!     I7SP   =  POINTER TO AREA WHERE SEQUENCED ACTIVE COLUMN INDICES
!               ARE STORED
!     B      =  UPPER HALF-BAND
!     BBAR   =  LOWER HALF-BAND
!     C      =  NUMBER OF ACTIVE COLUMNS
!     CBAR   =  NUMBER OF ACTIVE ROWS
!     R      =  NUMBER OF COLUMNS OF L THAT CAN BE STORED IN CORE
!     JPOS   =  CURRENT PIVOTAL COLUMN INDEX
!     JPOSL  =  NEXT COLUMN OF L TO BE WRITTEN OUT
!     LCOL   =  NUMBER OF COLUMNS OF L CURRENTLY STORED IN CORE OR ON
!               SCRATCH FILES
!     CCOUNT =  CURRENT NUMBER OF ACTIVE COLUMNS
!     CBCNT  =  CURRENT NUMBER OF ACTIVE ROWS
!     ITRN   =  ROW INDEX OF NEXT ACTIVE COLUMN ELEMENT
!     JTRN   =  COLUMN INDEX  OF NEXT ACTIVE COLUMN ELEMENT
!     IOFF   =  ROW POSITION OF THE FIRST ELEMENT IN AREA II
!     ITERM  =  IF NONZERO, TERMINATE BEFORE THE RE-WRITE
!     NCOL   =  SIZE OF THE INPUT MATRIX
!     BBBAR  =  B + BBAR
!     BBAR1  =  BBAR + 1
!     BBBAR1 =  B+BBAR - 1
!     SCRFLG =  NONZERO MEANS SPILL

!     ****************************************************************
!     RE-WRITE THE UPPER TRIANGLE OF ACTIVE ELEMENTS IN THE TRANSPOSED
!     ORDER
!     ****************************************************************

parm(2) = filea(1)
CALL OPEN (*1670,filea(1),ix(bufa),rdrew)
ccount = 0
IF (c == 0) GO TO 40
CALL transp (ix(1),x(1),nx,filea(1),b,sr1fil)

!     ZERO CORE

40 DO  i = 1,END
x(i) = 0.
END DO
IF (c == 0) GO TO 260

!     ****************************************************************
!     OPEN THE FILE CONTAINING THE TRANSPOSED ACTIVE ELEMENTS AND READ I
!     THE FIRST BBAR + 1 ROWS
!     ****************************************************************

parm(2) = sr1fil
CALL OPEN (*1670,sr1fil,ix(sr1buf),rd)
k = 0
60 CALL READ (*1680,*1690,sr1fil,itran(1),4,0,flag)
IF (itrn > 0) GO TO 70
CALL CLOSE (sr1fil,rew)
GO TO 140
70 IF (itrn > k+1) GO TO 130

!     DETERMINE IF COLUMN IS ALREADY ACTIVE

IF (jtrn <= bbbar) GO TO 60
kk  = 0
80 in1 = i3sp + kk
IF (ix(in1) == jtrn) GO TO 90
kk  = kk + 1
IF (kk-c < 0) THEN
  GO TO    80
ELSE IF (kk-c == 0) THEN
  GO TO   100
ELSE
  GO TO  1700
END IF

!     ADD IN ACTIVE ELEMENT TO EXISTING COLUMN

90 in1 = i3 + kk*bbar1 + k
dx(in1) = dtrn
GO TO 60

!     CREATE NEW ACTIVE COLUMN

100 ccount = ccount + 1
kk  = 0
110 in1 = i3sp + kk
IF (ix(in1) == 0) GO TO 120
kk  = kk + 1
IF (kk-c < 0) THEN
  GO TO   110
ELSE
  GO TO  1700
END IF
120 ix(in1) = jtrn
in1     = in1 + c
ix(in1) = k+1
in1     = i3 + kk*bbar1 + k
dx(in1) = dtrn
GO TO 60
130 k = k + 1
IF (k-bbar1 < 0) THEN
  GO TO    70
ELSE IF (k-bbar1 == 0) THEN
  GO TO   140
ELSE
  GO TO  1700
END IF

!     SET INDEXES IN AREA VII TO POINT TO THE ACTIVE COLUMNS IN SEQUENCE

140 ASSIGN 260 TO kk
150 in1 = i7sp
k   = 0
160 in2 = i3sp + k
IF (ix(in2) < 0) THEN
  GO TO  1700
ELSE IF (ix(in2) == 0) THEN
  GO TO   180
ELSE
  GO TO   190
END IF
170 in1 = in1 + 1
180 k   = k + 1
IF (k-c < 0) THEN
  GO TO   160
ELSE IF (k-c == 0) THEN
  GO TO   250
ELSE
  GO TO  1700
END IF
190 IF (in1 /= i7sp) GO TO 200
ix(in1) = k
GO TO 170
200 kkk = 0
210 in3 = in1 -kkk
IF (in3 > i7sp) GO TO 220
ix(in3) = k
GO TO 170
220 in4 = i3sp + ix(in3-1)
IF (ix(in2)-ix(in4) < 0) THEN
  GO TO   240
ELSE IF (ix(in2)-ix(in4) == 0) THEN
  GO TO  1700
END IF
230 ix(in3) = k
GO TO 170
240 ix(in3) = ix(in3-1)
kkk = kkk + 1
GO TO 210
250 GO TO kk, (260,1560)
260 CONTINUE

!     INITIALIZE

sr2fl = fileu(1)
sr3fl = sr3fil
jpos  = 1
parm(2) = filea(1)
CALL fwdrec (*1680,filea(1))
lcol  = 0
cbcnt = 0
jposl = 0
270 IF (jpos > ncol) GO TO 1650
!****************************************************************
!     READ NEXT COLUMN OF A INTO AREA II
!****************************************************************
ioff   = MAX0(1,jpos-bbbar1)
count  = cbcnt
imhere = 275
CALL intpk (*1710,filea(1),0,rdp,0)
k = 1
IF (jpos > bbbar) k = jpos - b + 1
280 IF (eol == 0.0) THEN
  GO TO   290
ELSE
  GO TO   400
END IF
290 CALL zntpki
IF (ii < k) GO TO 280
k = jpos + bbar
300 IF (ii > k) GO TO 330

!     READ ELEMENTS WITHIN THE BAND INTO AREA II

in1 = i2 - ioff + ii
dx(in1) = da
310 IF (eol == 0.0) THEN
  GO TO   320
ELSE
  GO TO   400
END IF
320 CALL zntpki
GO TO 300

!     TAKE CARE OF ACTIVE ELEMENTS BELOW THE BAND

330 kk  = 0
340 in1 = i4sp + kk
IF (ix(in1)-ii == 0) THEN
  GO TO   360
END IF
350 kk  = kk + 1
IF (kk-cbar < 0) THEN
  GO TO   340
ELSE IF (kk-cbar == 0) THEN
  GO TO   370
ELSE
  GO TO  1700
END IF

!     ADD IN ACTIVE ELEMENT TO EXISTING ROW

360 in1 = i4 + (kk+1)*bbbar - 1
dx(in1) = da
GO TO 310

!     CREATE NEW ACTIVE ROW

370 kk  = 0
380 in1 = i4sp + kk
IF (ix(in1) == 0) GO TO 390
kk  = kk + 1
IF (kk-cbar < 0) THEN
  GO TO   380
ELSE
  GO TO  1700
END IF
390 ix(in1) = ii
in1     = in1 + cbar
ix(in1) = jpos
in1     = i4 + (kk+1)*bbbar - 1
dx(in1) = da
cbcnt   = cbcnt + 1
GO TO 310

!     ARRANGE ACTIVE ROW INDEXES IN SEQUENCE AND STORE THEM IN AREA VI

400 IF (count == cbcnt) GO TO 500
in1 = i6sp
k   = 0
410 in2 = i4sp + k
IF (ix(in2) < 0) THEN
  GO TO  1700
ELSE IF (ix(in2) == 0) THEN
  GO TO   430
ELSE
  GO TO   440
END IF
420 in1 = in1 + 1
430 k   = k + 1
IF (k-cbar < 0) THEN
  GO TO   410
ELSE IF (k-cbar == 0) THEN
  GO TO   500
ELSE
  GO TO  1700
END IF
440 IF (in1 /= i6sp) GO TO 450
ix(in1) = k
GO TO 420
450 kk  = 0
460 in3 = in1 - kk
IF (in3 > i6sp) GO TO 470
ix(in3) = k
GO TO 420
470 in4 = i4sp + ix(in3-1)
IF (ix(in2)-ix(in4) < 0) THEN
  GO TO   490
ELSE IF (ix(in2)-ix(in4) == 0) THEN
  GO TO  1700
END IF
480 ix(in3) = k
GO TO 420
490 ix(in3) = ix(in3-1)
kk = kk + 1
GO TO 460
500 CONTINUE

!     TEST FOR POSSIBLE MERGING BETWEEN AN INACTIVE-ACTIVE COLUMN AND
!     THE CURRENT PIVOTAL COLUMN

IF (ccount == 0) GO TO 600
in1 = ix(i7sp) + i3sp
IF (ix(in1)-jpos < 0) THEN
  GO TO  1700
ELSE IF (ix(in1)-jpos == 0) THEN
  GO TO   510
ELSE
  GO TO   600
END IF

!     MERGE ACTIVE COLUMN AND CURRENT PIVOTAL COLUMN AND ZERO THAT
!     ACTIVE COLUMN IN AREA III

510 ix(in1) = 0
in1     = in1 + c
ix(in1) = 0
in1     = i3 + ix(i7sp)*bbar1
ccount  = ccount - 1
kk  = 0
520 in2 = in1 + kk
in3 = i2 + kk
dx(in3) = dx(in3) + dx(in2)
dx(in2) = 0.d0
kk = kk + 1
IF (kk-bbar1 < 0) THEN
  GO TO   520
ELSE IF (kk-bbar1 == 0) THEN
  GO TO   530
ELSE
  GO TO  1700
END IF

!     MERGE INTERACTION ELEMENTS

530 CONTINUE
IF (cbcnt == 0) GO TO 580
in1 = i5 + ix(i7sp)*cbar
k   = 0
540 in2 = i4sp + k
IF (ix(in2) == 0) GO TO 560
in3 = in1 + k
IF (dx(in3) == 0.d0) GO TO 560
IF (ix(in2) > jpos+bbar) GO TO 570

!     STORE ELEMENT WITHIN THE LOWER BAND

in2 = i2 + ix(in2) - ioff
dx(in2) = dx(in2) - dx(in3)
550 dx(in3) = 0.d0
560 k = k + 1
IF (k-cbar < 0) THEN
  GO TO   540
ELSE IF (k-cbar == 0) THEN
  GO TO   580
ELSE
  GO TO  1700
END IF

!     STORE ELEMENT IN THE ACTIVE ROW

570 in2 = i4 + (k+1)*bbbar - 1
dx(in2) = dx(in2) - dx(in3)
dx(in3) = 0.d0
GO TO 550

!     MOVE THE POINTERS IN AREA VII UP ONE

580 in1 = i7sp + ccount - 1
DO  i = i7sp,in1
  ix(i    ) = ix(i+1)
END DO
ix(in1+1) = 0
600 IF(lcol == 0)GO TO 820

!     ****************************************************************
!     OPERATE ON THE CURRENT COLUMN OF A BY ALL PREVIOUS COLUMNS OF L,
!     MAKING NOTED INTERCHANGES AS YOU GO
!     ****************************************************************

IF (scrflg == 0) GO TO 630
IF (lcol-(r-1) < 0) THEN
  GO TO   630
ELSE IF (lcol-(r-1) == 0) THEN
  GO TO   620
END IF
610 parm(2) = sr2fl
CALL OPEN (*1670,sr2fl,ix(sr2buf),rd)
620 parm(2) = sr3fl
CALL OPEN (*1670,sr3fl,ix(sr3buf),wrtrew)
630 ll   = 0
lll  = 0
llll = 0

!     PICK UP INTERCHANGE INDEX FOR COLUMN JPOSL + LL + 1

640 in1    = i1sp + ll
intchn = ix(in1)
in2    = i2 + ll
IF (intchn == 0) GO TO 650

!     PERFORM ROW INTERCHANGE

in1 = in2 + intchn
da  = dx(in1)
dx(in1) = dx(in2)
dx(in2) = da
650 CONTINUE

!     COMPUTE THE CONTRIBUTION FROM THAT COLUMN

END = MIN0(bbar1,ncol-(jposl+ll))
END = END - 1
IF (dx(in2) == 0.0) THEN
  GO TO   710
END IF
660 in1 = i1 + lll*bbar
CALL dloop (dx(in2+1),dx(in1),-dx(in2),END)
IF (cbcnt == 0) GO TO 710

!     TEST TO SEE IF AN INACTIVE-ACTIVE ROW CONTRIBUTION SHOULD BE
!     ADDED IN

kkk = 0
680 in3 = i6sp + kkk
in1 = ix(in3) + i4sp
IF (ix(in1) > jpos+bbar) GO TO 710
kk = in1 + cbar
IF (ix(kk) > jposl+ll+1) GO TO 700
IF (ix(in1)-jposl-bbar1 <= ll) GO TO 700

!     ADD IN EFFECT OF THE INACTIVE-ACTIVE ROW

in4 = i2 + ix(in1) - ioff
k   = jposl + bbbar - jpos  + ll + i4 + ix(in3)*bbbar
dx(in4) = dx(in4) - dx(k)*dx(in2)
700 kkk = kkk + 1
IF (kkk < cbcnt) GO TO 680
710 ll  = ll  + 1
lll = lll + 1
IF (ll == lcol) GO TO 770
IF (ll-r+1 < 0) THEN
  GO TO   640
ELSE IF (ll-r+1 == 0) THEN
  GO TO   720
ELSE
  GO TO   750
END IF
720 IF (r == bbbar1) GO TO 640
in1  = i1  + ll*bbar
740 icrq = in1 + bbar*2 - 1 - sr3buf
IF (icrq > 0) GO TO 1668
ibbar2 = bbar*2
CALL READ (*1680,*1690,sr2fl,dx(in1),ibbar2,0,flag)
GO TO 640
750 in1 = i1 + (lll-1)*bbar
IF (ll == r .AND. lcol == bbbar1) GO TO 760
CALL WRITE (sr3fl,dx(in1),2*bbar,0)
760 lll = lll - 1
GO TO 740
770 CONTINUE

!     COMPUTE ELEMENTS FOR THE ACTIVE ROWS

IF (cbcnt == 0) GO TO 820
k   = 0
780 in1 = i4sp + k
IF (ix(in1) > jpos+bbar) GO TO 800
790 k   = k + 1
IF (k-cbar < 0) THEN
  GO TO   780
ELSE IF (k-cbar == 0) THEN
  GO TO   820
ELSE
  GO TO  1700
END IF
800 in1 = in1 + cbar
IF (ix(in1) == jpos) GO TO 790
kkk = MAX0(0,bbbar-jpos+ix(in1)-1)
in2 = i4  + k*bbbar - 1
in3 = i2  + kkk - 1 - MAX0(0,bbbar-jpos)
in1 = in2 + bbbar
in2 = in2 + kkk
810 in2 = in2 + 1
kkk = kkk + 1
in3 = in3 + 1
dx(in1) = dx(in1)-dx(in2)*dx(in3)
IF (kkk-bbbar1 < 0) THEN
  GO TO   810
ELSE IF (kkk-bbbar1 == 0) THEN
  GO TO   790
ELSE
  GO TO  1700
END IF

!     SEARCH THE LOWER BAND FOR THE MAXIMUM ELEMENT AND INTERCHANGE
!     ROWS TO BRING IT TO THE DIAGONAL

820 k   = 1
in1 = i2 + jpos - ioff
MAX = DABS(dx(in1))
END = MIN0(bbar1,ncol-jpos+1)
intchn = 0
IF (END == 1) GO TO 860
830 in2 = in1 + k
IF (DABS(dx(in2)) > MAX) GO TO 850
840 k = k + 1
IF (k-END < 0) THEN
GO TO   830
ELSE IF (k-END == 0) THEN
GO TO   860
ELSE
GO TO  1700
END IF
850 MAX = DABS(dx(in2))
intchn = k
GO TO 840

860 IF (intchn == 0) GO TO 870

!     INTERCHANGE ROWS IN AREA II

det = -det

MAX = dx(in1)
in2 = in1 + intchn
dx(in1) = dx(in2)
dx(in2) = MAX

!     STORE THE PERMUTATION INDEX

in2 = i1sp + lcol
ix(in2) = intchn

!     DIVIDE THE LOWER BAND BY THE DIAGONAL ELEMENT

870 imhere = 870
IF (dx(in1) == 0.d0) GO TO 1710
MAX = 1.d0/dx(in1)
mindia = DMIN1(DABS(dx(in1)),mindia)
880 IF (DABS(det) <= 10.d0) GO TO 890
det   = det/10.d0
power = power + 1
GO TO 880
890 IF (DABS(det) >= .1D0) GO TO 900
det   = det*10.d0
power = power - 1
GO TO 890
900 det = det*dx(in1)
k   = 1
END = MIN0(bbar1,ncol-jpos+1)
IF (END == 1) GO TO 920
910 in2 = in1 + k
dx(in2) = dx(in2)*MAX
k = k + 1
IF (k-END < 0) THEN
GO TO   910
ELSE IF (k-END == 0) THEN
GO TO   920
ELSE
GO TO  1700
END IF
920 IF (cbcnt == 0) GO TO 940

!     DIVIDE THE ACTIVE ROWS BY THE DIAGONAL

k   = 0
in1 = i4 + bbbar1
930 dx(in1) = dx(in1)*MAX
in1 = in1 + bbbar
k   = k + 1
IF (k-cbar < 0) THEN
  GO TO   930
ELSE IF (k-cbar == 0) THEN
  GO TO   940
ELSE
  GO TO  1700
END IF
940 CONTINUE

!     INTERCHANGE ACTIVE COLUMNS AND ADD IN EFFECT OF THE COLUMN OF L
!     ABOUT TO BE WRITTEN OUT

IF (ccount == 0) GO TO 990
IF (jpos < bbbar) GO TO 990
intch = ix(i1sp)
k   = 0
950 in1 = i3sp + k
IF (intch == 0) GO TO 960
in1 = i3  + k*bbar1
in2 = in1 + intch
da  = dx(in1)
dx(in1) = dx(in2)
dx(in2) = da
960 kk  = 1
in2 = i1 - 1
in1 = i3 + k*bbar1
IF (dx(in1) == 0.d0) GO TO 980
970 in3 = in1 + kk
in4 = in2 + kk
dx(in3) = dx(in3) - dx(in1)*dx(in4)
kk = kk + 1
IF (kk-bbar1 < 0) THEN
  GO TO   970
ELSE IF (kk-bbar1 == 0) THEN
  GO TO   980
ELSE
  GO TO  1700
END IF
980 k = k + 1
IF (k-c < 0) THEN
  GO TO   950
ELSE IF (k-c == 0) THEN
  GO TO   990
ELSE
  GO TO  1700
END IF

!     WRITE OUT THE NEXT COLUMN OF U AND THE ROW OF ACTIVE ELEMENTS

990 parm(2) = sr2fil
CALL bldpk (rdp,typel,sr2fil,0,0)
in1 = i2
jj  = ioff
imhere = 1030
1000 dz = dx(in1)
IF (dz == 0.0) THEN
  GO TO  1020
END IF
1010 CALL zblpki
1020 in1 = in1 + 1
jj  = jj  + 1
IF (jj-jpos > 0) THEN
  GO TO  1030
ELSE
  GO TO  1000
END IF
1030 IF (dx(in1-1) == 0.0) THEN
  GO TO  1710
END IF
1040 CONTINUE

!     PACK ACTIVE COLUMN ELEMENTS ALSO

IF (ccount ==   0) GO TO 1080
IF (jpos < bbbar) GO TO 1080
k   = 0
1050 in1 = i7sp + k
in2 = ix(in1) + i3sp
GO TO 1070
1060 k = k + 1
IF (k-ccount < 0) THEN
  GO TO  1050
ELSE IF (k-ccount == 0) THEN
  GO TO  1080
ELSE
  GO TO  1700
END IF
1070 in3 = i3 + ix(in1)*bbar1
dz  = dx(in3)
IF (dz == 0.d0) GO TO 1060
jj = ix(in2)
CALL zblpki
GO TO 1060
1080 CALL bldpkn (sr2fil,0,fileu)

!     COMPUTE ACTIVE ROW-COLUMN INTERACTION

IF (ccount == 0 .OR. cbcnt == 0) GO TO 1130
IF (jpos < bbbar) GO TO 1130
k = 0
1090 CONTINUE
in1 = i3 + k*bbar1
IF (dx(in1) == 0.d0) GO TO 1120
kk  = 0
1100 in2 = i4sp + kk
in2 = i4 + kk*bbbar
IF (dx(in2) == 0.d0) GO TO 1110
in3 = i5 + k*cbar + kk
dx(in3) = dx(in3)+dx(in2)*dx(in1)
1110 kk = kk + 1
IF (kk-cbar < 0) THEN
  GO TO  1100
ELSE IF (kk-cbar == 0) THEN
  GO TO  1120
ELSE
  GO TO  1700
END IF
1120 k  = k  + 1
IF (k-c < 0) THEN
  GO TO  1090
ELSE IF (k-c == 0) THEN
  GO TO  1130
ELSE
  GO TO  1700
END IF

!     MOVE ELEMENTS IN AREA III UP ONE CELL

1130 IF (ccount == 0) GO TO 1180
IF (jpos < bbbar) GO TO 1180
k   = 0
1140 in1 = i3sp + k
IF (ix(in1) == 0) GO TO 1170
kk  = 0
in1 = i3  + k*(bbar1)
1150 in2 = in1 + kk
dx(in2) = dx(in2+1)
kk = kk + 1
IF (kk-bbar < 0) THEN
  GO TO  1150
ELSE IF (kk-bbar == 0) THEN
  GO TO  1160
ELSE
  GO TO  1700
END IF
1160 dx(in2+1) = 0.d0
1170 k  = k + 1
IF (k-c < 0) THEN
  GO TO  1140
ELSE IF (k-c == 0) THEN
  GO TO  1180
ELSE
  GO TO  1700
END IF

!     DETERMINE IF A COLUMN OF L CAN BE WRITTEN OUT

1180 IF (lcol-bbbar1 < 0) THEN
  GO TO  1360
END IF

!     OUTPUT A COLUMN OF L

1190 parm(2) = filel(1)
jposl   = jposl + 1
CALL bldpk (rdp,typel,filel(1),0,0)

!     STORE THE PERMUTATION INDEX AS THE DIAGONAL ELEMENT

jj = jposl
dz = ix(i1sp)
CALL zblpki
k  = 0
1200 jj = jposl + k + 1
in2= i1  + k
dz = dx(in2)
IF (dz == 0.0) THEN
  GO TO  1220
END IF
1210 CALL zblpki
1220 k = k + 1
IF (k-bbar < 0) THEN
  GO TO  1200
ELSE IF (k-bbar == 0) THEN
  GO TO  1230
ELSE
  GO TO  1700
END IF

!     PACK ACTIVE ROW ELEMENTS ALSO

1230 IF (cbcnt == 0) GO TO 1270
k   = 0
1240 in1 = i6sp + k
in2 = i4 + ix(in1)*bbbar
in1 = ix(in1) + i4sp
jj  = ix(in1)
dz  = dx(in2)
IF (dz == 0.d0) GO TO 1260
CALL zblpki
1260 k = k + 1
IF (k-cbcnt < 0) THEN
  GO TO  1240
ELSE IF (k-cbcnt == 0) THEN
  GO TO  1270
ELSE
  GO TO  1700
END IF
1270 CALL bldpkn (filel,0,filel)

!     MOVE PERMUTATION INDICES OVER ONE ELEMENT

END = i1sp + lcol
DO  i = i1sp,END
ix(i) = ix(i+1)
END DO

!     MOVE ELEMENTS IN AREA I OVER ONE COLUMN

k = 0
IF (scrflg == 0) GO TO 1300
CALL CLOSE (sr2fl,rew)
IF (r > 2) GO TO 1300
icrq = i1 + bbar*2 - 1 - sr3buf
IF (icrq > 0) GO TO 1668
CALL OPEN (*1670,sr2fl,ix(sr2buf),rd)
ibbar2 = 2*bbar
CALL READ (*1680,*1690,sr2fl,dx(i1),ibbar2,0,flag)
GO TO 1350
1300 in1 = i1  + k*bbar
in2 = in1 + bbar
CALL xloop (dx(in1),dx(in2),bbar)
k = k + 1
IF (k-r+2 < 0) THEN
  GO TO  1300
ELSE IF (k-r+2 == 0) THEN
  GO TO  1330
ELSE
  GO TO  1350
END IF
1330 IF (r-bbbar1 < 0.0) THEN
  GO TO  1340
ELSE IF (r-bbbar1 == 0.0) THEN
  GO TO  1300
ELSE
  GO TO  1700
END IF
1340 icrq = in2 + bbar*2 - 1 - sr3buf
IF (icrq > 0) GO TO 1668
CALL OPEN (*1670,sr2fl,ix(sr2buf),rd)
ibbar2 = bbar*2
CALL READ (*1680,*1690,sr2fl,dx(in2),ibbar2,0,flag)
1350 lcol = lcol - 1

!     STORE CURRENT COLUMN OF L

1360 IF (cbcnt == 0) GO TO 1410

!     MOVE ELEMENTS IN AREA IV UP ONE CELL

k   = 0
1370 in1 = i4sp + k
IF (ix(in1) == 0) GO TO 1400
kk  = 0
in1 = i4 + k*bbbar
1380 in2 = in1 + kk
dx(in2) = dx(in2+1)
kk = kk + 1
IF (kk-bbbar1 < 0) THEN
  GO TO  1380
ELSE IF (kk-bbbar1 == 0) THEN
  GO TO  1390
ELSE
  GO TO  1700
END IF
1390 dx(in2+1) = 0.d0
1400 k = k + 1
IF (k-cbar < 0) THEN
  GO TO  1370
ELSE IF (k-cbar == 0) THEN
  GO TO  1410
ELSE
  GO TO  1700
END IF
1410 IF (scrflg /= 0) GO TO 1440

!     STORE COLUMN IN CORE

1420 in1 = i1 + lcol*bbar
END = MIN0(bbar,ncol-jpos)
IF (END == 0) GO TO 1470
k   = 0
in3 = i2 + jpos - ioff + 1
1430 in2 = in1 + k
in4 = in3 + k
dx(in2) = dx(in4)
k = k + 1
IF (k-END < 0) THEN
GO TO  1430
ELSE IF (k-END == 0) THEN
GO TO  1470
ELSE
GO TO  1700
END IF

!     STORE COLUMN ON THE SCRATCH FILE

1440 IF (lcol-r+1 < 0) THEN
  GO TO  1420
ELSE IF (lcol-r+1 == 0) THEN
  GO TO  1460
END IF
1450 in1 = i1 + (lll-1)*bbar
CALL WRITE (sr3fl,dx(in1),bbar*2,0)
1460 in1 = i2 + jpos - ioff + 1
CALL WRITE (sr3fl,dx(in1),bbar*2,0)

!     CLOSE SCRATCH FILES AND SWITCH THE POINTERS TO THEM

CALL CLOSE (sr3fl,rew)
CALL CLOSE (sr2fl,rew)
in1   = sr2fl
sr2fl = sr3fl
sr3fl = in1
1470 lcol  = lcol + 1
IF (c == 0) GO TO 1560
IF (jpos < bbbar) GO TO 1560

!     READ IN THE NEXT ROW OF ACTIVE COLUMN ELEMENTS

count = ccount
IF (itrn < 0) GO TO 1560
1480 IF (itrn > jpos-b+2) GO TO 1550

!     TEST TO SEE IF COLUMN IS ALREADY ACTIVE

k   = 0
1490 in1 = i3sp + k
IF (ix(in1) == jtrn) GO TO 1530
k   = k + 1
IF (k-c < 0) THEN
  GO TO  1490
ELSE IF (k-c == 0) THEN
  GO TO  1500
ELSE
  GO TO  1700
END IF

!     CREATE A NEW ACTIVE COLUMN

1500 k   = 0
1510 in1 = i3sp + k
IF (ix(in1) == 0) GO TO 1520
k   = k + 1
IF (k-c < 0) THEN
  GO TO  1510
ELSE
  GO TO  1700
END IF
1520 ix(in1) = jtrn
in1     = in1 + c
ix(in1) = itrn
in1     = i3 + (k+1)*bbar1 - 1
dx(in1) = dtrn
ccount  = ccount + 1
GO TO 1540

!     STORE ELEMENT IN EXISTING COLUMN

1530 in1 = i3 + (k+1)*bbar1 - 1
dx(in1) = dx(in1) + dtrn
1540 CALL READ (*1680,*1690,sr1fil,itran,4,0,flag)
IF (itrn > 0) GO TO 1480
CALL CLOSE (sr1fil,rew)
1550 IF (ccount == count) GO TO 1560

!     RE-ARRANGE INDEXES IN SEQUENTIAL ORDER

ASSIGN 1560 TO kk
GO TO 150
1560 CONTINUE
jpos = jpos + 1

!     ZERO AREA II

END = i2 + MIN0(jpos-ioff+bbar-1,ncol-1)
DO  i = i2,END
dx(i) = 0.d0
END DO

!      TEST TO SEE IF ROW INTERACTION ELEMENTS WILL MERGE INTO AREA III

IF (cbcnt  == 0) GO TO 270
IF (ccount == 0) GO TO 1620
IF (jpos-1 < bbbar) GO TO 270
in1 = i4sp
k   = 0
1590 in2 = in1 + k
IF (ix(in2) == jpos-b+1) GO TO 1600
k   = k + 1
IF (k < cbar) GO TO 1590
GO TO 270
1600 in1 = i5 + k
in2 = i3 + bbar
k   = 0
1610 dx(in2) = dx(in2)-dx(in1)
dx(in1) = 0.d0
in2 = in2 + bbar1
in1 = in1 + cbar
k   = k + 1
IF (k < c) GO TO 1610

!      TEST TO SEE IF ACTIVE ROW HAS BEEN ELIMINATED

1620 in1 = ix(i6sp) + i4sp
IF (ix(in1)-jposl-bbar1 == 0) THEN
  GO TO  1630
ELSE
  GO TO   270
END IF

!     ELIMINATE THE ACTIVE ROW

1630 ix(in1) = 0
in1     = in1 + cbar
ix(in1) = 0
cbcnt   = cbcnt - 1

!     MOVE INDEXES IN AREA VI UP ONE

in1 = i6sp + cbcnt - 1
DO  i = i6sp,in1
  ix(i    ) = ix(i+1)
END DO
ix(in1+1) = 0
GO TO 270

!     FINISH WRITING OUT THE COMPLETED COLUMNS OF L

1650 CONTINUE
CALL CLOSE (sr1fil,rew)
CALL CLOSE (filel,norew)
CALL CLOSE (sr2fil,norew)
parm(5) = iend
CALL conmsg (parm(3),3,0)
CALL finwrt (iterm,scrflg,sr2fl,jposl,i1sp,bbar,i1,cbcnt,ipak,r,  &
    bbbar1,bbbar,i6sp,i4,i4sp,ix,dx,x,lcol)
fileu(7) = bbbar
RETURN

!     ERROR EXITS

1660 parm(1) = -7
GO TO 1720
1668 parm(1) = -8
parm(2) = icrq
GO TO 1720
1670 parm(1) = -1
GO TO 1720
1680 parm(1) = -2
GO TO 1720
1690 parm(1) = -3
GO TO 1720
1700 parm(1) = -25
GO TO 1720

!     SINGULAR MATRIX - CLOSE ALL FILES AND RETURN TO USER

1710 CALL CLOSE (filea(1),rew)
CALL CLOSE (filel(1),rew)
CALL CLOSE (fileu(1),rew)
CALL CLOSE (sr1fil,rew)
CALL CLOSE (sr2fil,rew)
CALL CLOSE (sr3fil,rew)
WRITE  (nout,1715) imhere
1715 FORMAT (/60X,'DECOMP/IMHERE@',i5)
!WKBA 4/95 SPR94018
fileu(7) = bbbar
RETURN 1
1720 CALL mesage (parm(1),parm(2),parm(3))
RETURN
END SUBROUTINE decomp
