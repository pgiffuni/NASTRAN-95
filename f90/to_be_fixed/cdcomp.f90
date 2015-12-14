SUBROUTINE cdcomp (*,ix,x,dx)
     
!     CDCOMP WILL DECOMPOSE A COMPLEX UNSYMETRIC MATRIX INTO A UNIT
!     LOWER TRIANGULAR MATRIX AND AN UPPER TRIANGULAR MATRIX,USING
!     PARTIAL PIVOTING WITHIN THE LOWER BAND
 
!     IMPORTANT - CALLER MUST FIRST INITIALIZE B AND/OR BBAR IN /CDCMPX/
 
!     DEFINITION OF INPUT PARAMETERS
 
!     FILEA  =  MATRIX CONTROL BLOCK FOR THE  INPUT MATRIX A
!     FILEL  =  MATRIX CONTROL BLOCK FOR THE OUTPUT MATRIX L
!     FILEU  =  MATRIX CONTROL BLOCK FOR THE OUTPUT MATRIX U
!     SR1FIL =  SCRATCH FILE
!     SR2FIL =  SCRATCH FILE
!     SR3FIL =  SCRATCH FILE
!     NX     =  NUMBER OF CELLS OF CORE AVAILABLE AT IX
!     DET    =  CELL WHERE THE DETERMINATE OF A WILL BE STORED
!     POWER  =  SCALE FACTOR TO BE APPLIED TO THE DETERMINATE
!               ( DETERMINATE = DET*10**POWER )
!     MINDIA =  CELL WHERE THE VALUE OF THE MINIMUM DIAGONAL WILL BE
!               STORED
!     IX     =  BLOCK OF CORE AVAILABLE AS WORKING STORAGE TO DECOMP
!     X      =  SAME BLOCK AS IX, BUT TYPED REAL
!     DX     =  SAME BLOCK AS IX, BUT TYPED DOUBLE PRECISION
 
 
 
 , INTENT(IN OUT)                         :: *
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
    sr2fil    ,sr3fil   ,sqr      ,sym      , flag      ,itran(6)
DOUBLE PRECISION :: dz(2)     ,da(2)    ,det      ,MAX(2)   ,  &
    mindia    , dtrn(2)  ,dx1      , dx2       ,limit

CHARACTER (LEN=29) :: uim
CHARACTER (LEN=25) :: uwm
CHARACTER (LEN=23) :: ufm
COMMON   /xmssg /  ufm       ,uwm      ,uim
COMMON   /cdcmpx/  filea(7)  ,filel(7) ,fileu(7) ,sr1fil   ,  &
    sr2fil    ,sr3fil   ,det(2)   ,power    ,  &
    nx        ,mindia   ,b        ,bbar     , c         ,cbar     ,r
COMMON   /system/  sysbuf    ,nout
COMMON   /names /  rd        ,rdrew    ,wrt      ,wrtrew   ,  &
    rew       ,norew    ,eofnrw   ,rsp      ,  &
    rdp       ,csp      ,cdp      ,sqr      ,  &
    rect      ,diag     ,lowtri   ,uprtri   , sym       ,row      ,ident
COMMON   /zntpkx/  a(4)      ,ii       ,eol
!     COMMON   /DESCRP/  LENGTH    ,MAJOR
COMMON   /zblpkx/  z(4)      ,jj
COMMON   /unpakx/  itypex    ,ixy      ,jxy      ,incrx
COMMON   /packx /  itype1    ,itype2   ,iy       ,jy       , incry
EQUIVALENCE        (da(1),a(1))        ,(dz(1),z(1))       ,  &
    (forma,filea(4))    ,(typea,filea(5))   ,  &
    (ncol,filea(3))     ,(typel,filel(5))
EQUIVALENCE        (itran(1),itrn)     ,(itran(2),jtrn)    ,  &
    (itran(3),dtrn(1))
DATA      parm(3), parm(4)/  4HCDCO,4HMP  /
DATA      ibegn  / 4HBEGN /, iend  /4HEND /
DATA      limit  / 1.0D-38/

!     BUFFER ALLOCATION

bufa   = nx     - sysbuf
ibufl  = bufa   - sysbuf
outbuf = ibufl  - sysbuf
sr1buf = outbuf - sysbuf
sr2buf = sr1buf - sysbuf
sr3buf = sr2buf - sysbuf
icrq   =-sr3buf
IF (icrq > 0) GO TO 1715
det(1) = 1.d0
det(2) = 0.d0
power  = 0
mindia = 1.d+25
iterm  = 0
IF (filea(1) < 0) iterm = 1
filea(1) = IABS(filea(1))

!     WRITE THE HEADER RECORD ON THE OUTPUT TAPES AND INITIALIZE THE
!     TRAILER RECORDS.

CALL gopen (filel,ix(ibufl),wrtrew)
parm(2) = sr2fil
CALL OPEN (*1680,sr2fil,ix(outbuf),wrtrew)
CALL fname (fileu(1),x(1))
CALL WRITE (sr2fil,x(1),2,1)
filel(2) = 0
filel(3) = ncol
filel(4) = 4
filel(6) = 0
filel(7) = 0
fileu(2) = 0
fileu(3) = ncol
fileu(4) = 5
fileu(6) = 0
fileu(7) = 0
IF (ncol > 2) GO TO 10
CALL com12 (*1720,ix(1),x(1),dx(1),iterm)
parm(5) = iend
CALL conmsg (parm(3),3,0)
RETURN

!     CALL GENVEC TO PICK B, BBAR, C, CBAR, AND R

10 IF (b <= 0 .OR. bbar <= 0) CALL genvec (*1720,ix(bufa),filea(1),  &
    nx,ix(1),ncol,b,bbar,c,cbar,r,2)
bbar1  = bbar + 1
bbbar  = MIN0(b+bbar,ncol)
bbbar1 = bbbar - 1
scrflg = 0
IF (r < bbbar1) scrflg = 1
IF (scrflg == 0) GO TO 20
icrq = (bbbar1-r)*4*bbar
CALL page2 (2)
WRITE  (nout,15) uim,icrq
15 FORMAT (a29,' 2177, SPILL WILL OCCUR IN COMPLEX UNSYMMETRIC ',  &
    'DECOMPOSITION.', /i10, ' ADDITIONAL WORDS NEEDED TO STAY IN CORE.')

!     INITIALIZE POINTERS TO SPECIFIC AREAS OF CORE

20 i1   = 1
ipak = i1 + 2*bbar*r + bbbar/2 + 1
i1sp = bbar*r*4 + 1
i2   = ipak
i3sp = (i2+ 2*MIN0(ncol,bbbar+bbar))*2 - 1
i3   = i2 + 2*MIN0(ncol,bbbar+bbar) + c
i4sp = i3sp + (bbar+2)*c*4 - 2*c
i4   = i3 + 2*bbar1*c + cbar
i5   = i4 + 2*bbbar*cbar
i6sp = (i5 + 2*c*cbar)*2 - 1
i7sp = i6sp + cbar
parm(5) = ibegn
CALL conmsg (parm(3),3,0)
END  = i7sp + c

!     DEFINITION OF KEY PROGRAM PARAMETERS

!     I1     =  POINTER TO AREA WHERE COMPLETED COLUMNS OF L ARE STORE
!     I1SP   =  POINTER TO AREA WHERE THE PERMUTATION INDEXES ARE STOR
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
!     BBBAR1 =  B + BBAR - 1
!     SCRFLG =  NONZERO MEANS SPILL

!     RE-WRITE THE UPPER TRIANGLE OF ACTIVE ELEMENTS IN THE TRANSPOSED
!     ORDER

parm(2) = filea(1)
CALL OPEN (*1680,filea(1),ix(bufa),rdrew)
ccount = 0
IF (c == 0) GO TO 40
CALL ctrnsp (ix(1),x(1),nx,filea(1),b,sr1fil)

!     ZERO CORE

40 DO  i = 1,END
x(i) = 0.
END DO
IF (c == 0) GO TO 260

!     OPEN THE FILE CONTAINING THE TRANSPOSED ACTIVE ELEMENTS AND READ I
!     THE FIRST BBAR + 1 ROWS

parm(2) = sr1fil
CALL OPEN (*1680,sr1fil,ix(sr1buf),rd)
k = 0
60 CALL READ (*1690,*1700,sr1fil,itran(1),6,0,flag)
IF (itrn > 0) GO TO 70
CALL CLOSE (sr1fil,rew)
GO TO 140
70 IF (itrn > k+1) GO TO 130

!     DETERMINE IF COLUMN IS ALREADY ACTIVE

IF (jtrn <= bbbar) GO TO 60
kk  = 0
80 in1 = i3sp + kk
IF (ix(in1) == jtrn) GO TO 90
kk = kk + 1
IF (kk-c < 0) THEN
  GO TO    80
ELSE IF (kk-c == 0) THEN
  GO TO   100
ELSE
  GO TO  1710
END IF

!     ADD IN ACTIVE ELEMENT TO EXISTING COLUMN

90 in1 = i3 + 2*kk*bbar1 + k + k
dx(in1  ) = dtrn(1)
dx(in1+1) = dtrn(2)
GO TO 60

!     CREATE NEW ACTIVE COLUMN

100 ccount = ccount + 1
kk  = 0
110 in1 = i3sp + kk
IF (ix(in1) == 0) GO TO 120
kk = kk + 1
IF (kk-c < 0) THEN
  GO TO   110
ELSE
  GO TO  1710
END IF
120 ix(in1) = jtrn
in1     = in1 + c
ix(in1) = k  + 1
in1     = i3 + 2*kk*bbar1 + k + k
dx(in1  ) = dtrn(1)
dx(in1+1) = dtrn(2)
GO TO 60
130 k = k + 1
IF (k-bbar1 < 0) THEN
  GO TO    70
ELSE IF (k-bbar1 == 0) THEN
  GO TO   140
ELSE
  GO TO  1710
END IF

!     SET INDEXES IN AREA VII TO POINT TO THE ACTIVE COLUMNS IN SEQUENCE

140 ASSIGN 260 TO kk
150 in1 = i7sp
k   = 0
160 in2 = i3sp + k
IF (ix(in2) < 0) THEN
  GO TO  1710
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
  GO TO  1710
END IF
190 IF (in1 /= i7sp) GO TO 200
ix(in1) = k
GO TO 170
200 kkk = 0
210 in3 = in1 - kkk
IF (in3 > i7sp) GO TO 220
ix(in3) = k
GO TO 170
220 in4 = i3sp + ix(in3-1)
IF (ix(in2)-ix(in4) < 0) THEN
  GO TO   240
ELSE IF (ix(in2)-ix(in4) == 0) THEN
  GO TO  1710
END IF
230 ix(in3) = k
GO TO 170
240 ix(in3) = ix(in3-1)
kkk = kkk + 1
GO TO 210
250 GO TO kk, (260,1570)
260 CONTINUE

!     INITIALIZE

sr2fl = fileu(1)
sr3fl = sr3fil
jpos  = 1
parm(2) = filea(1)
CALL fwdrec (*1690,filea(1))
lcol  = 0
cbcnt = 0
jposl = 0
270 IF (jpos > ncol) GO TO 1670

!     READ NEXT COLUMN OF A INTO AREA II

ioff  = MAX0(1,jpos-bbbar1)
count = cbcnt
CALL intpk (*1720,filea(1),0,cdp,0)
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

in1 = i2 + 2*(ii-ioff)
dx(in1  ) = da(1)
dx(in1+1) = da(2)
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
  GO TO  1710
END IF

!     ADD IN ACTIVE ELEMENT TO EXISTING ROW

360 in1 = i4 + 2*(kk+1)*bbbar - 2
dx(in1  ) = da(1)
dx(in1+1) = da(2)
GO TO 310

!     CREATE NEW ACTIVE ROW

370 kk  = 0
380 in1 = i4sp + kk
IF (ix(in1) == 0) GO TO 390
kk  = kk + 1
IF (kk-cbar < 0) THEN
  GO TO   380
ELSE
  GO TO  1710
END IF
390 ix(in1) = ii
in1 = in1 + cbar
ix(in1) = jpos
in1 = i4 + (kk+1)*bbbar*2 - 2
dx(in1  ) = da(1)
dx(in1+1) = da(2)
cbcnt = cbcnt + 1
GO TO 310

!     ARRANGE ACTIVE ROW INDEXES IN SEQUENCE AND STORE THEM IN AREA VI

400 IF (count == cbcnt) GO TO 500
in1 = i6sp
k   = 0
410 in2 = i4sp + k
IF (ix(in2) < 0) THEN
  GO TO  1710
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
  GO TO  1710
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
  GO TO  1710
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
  GO TO  1710
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
in1     = i3 + ix(i7sp)*bbar1*2
ccount  = ccount - 1
kk  = 0
520 in2 = in1 + kk + kk
in3 = i2  + kk + kk
dx(in3  ) = dx(in3  ) + dx(in2)
dx(in3+1) = dx(in3+1) + dx(in2+1)
dx(in2  ) = 0.d0
dx(in2+1) = 0.d0
kk = kk + 1
IF (kk-bbar1 < 0) THEN
  GO TO   520
ELSE IF (kk-bbar1 == 0) THEN
  GO TO   530
ELSE
  GO TO  1710
END IF

!     MERGE INTERACTION ELEMENTS

530 CONTINUE
IF (cbcnt == 0) GO TO 580
in1 = i5 + 2*ix(i7sp)*cbar
k   = 0
540 in2 = i4sp + k
IF (ix(in2) == 0) GO TO 560
in3 = in1 + 2*k
IF (dx(in3) == 0.d0 .AND. dx(in3+1) == 0.d0) GO TO 560
IF (ix(in2) > jpos+bbar) GO TO 570

!     STORE ELEMENT WITHIN THE LOWER BAND

in2 = i2 + 2*(ix(in2)-ioff)
dx(in2  ) = dx(in2  ) - dx(in3)
dx(in2+1) = dx(in2+1) - dx(in3+1)
550 dx(in3  ) = 0.d0
dx(in3+1) = 0.d0
560 k = k + 1
IF (k-cbar < 0) THEN
  GO TO   540
ELSE IF (k-cbar == 0) THEN
  GO TO   580
ELSE
  GO TO  1710
END IF

!     STORE ELEMENT IN THE ACTIVE ROW

570 in2 = i4 + 2*(k+1)*bbbar - 2
dx(in2+1) = dx(in2+1) - dx(in3+1)
dx(in3+1) = 0.d0
dx(in2) = dx(in2) - dx(in3)
dx(in3) = 0.d0
GO TO 550

!     MOVE THE POINTERS IN AREA VII UP ONE

580 in1 = i7sp + ccount - 1
DO  i = i7sp,in1
  ix(i    ) = ix(i+1)
END DO
ix(in1+1) = 0
600 IF (lcol == 0) GO TO 830

!     OPERATE ON THE CURRENT COLUMN OF A BY ALL PREVIOUS COLUMNS OF L,
!     MAKING NOTED INTERCHANGES AS YOU GO

IF (scrflg == 0) GO TO 630
IF (lcol-(r-1) < 0) THEN
  GO TO   630
ELSE IF (lcol-(r-1) == 0) THEN
  GO TO   620
END IF
610 parm(2) = sr2fl
CALL OPEN (*1680,sr2fl,ix(sr2buf),rd)
620 parm(2) = sr3fl
CALL OPEN (*1680,sr3fl,ix(sr3buf),wrtrew)
630 ll   = 0
lll  = 0
llll = 0

!     PICK UP INTERCHANGE INDEX FOR COLUMN JPOSL + LL + 1

640 in1    = i1sp + ll
intchn = ix(in1)
in2    = i2 + ll + ll
IF (intchn == 0) GO TO 650

!     PERFORM ROW INTERCHANGE

in1     = in2 + 2*intchn
da(1)   = dx(in1)
dx(in1) = dx(in2)
dx(in2) = da(1)
da(1)   = dx(in1+1)
dx(in1+1) = dx(in2+1)
dx(in2+1) = da(1)
650 CONTINUE

!     COMPUTE THE CONTRIBUTION FROM THAT COLUMN

END = MIN0(bbar1,ncol-(jposl+ll))
IF (dx(in2) == 0.d0 .AND. dx(in2+1) == 0.d0) GO TO 720
in1 = i1 + 2*lll*bbar
CALL cloop (dx(in2+2),dx(in1),dx(in2),END-1)
IF (cbcnt == 0) GO TO 720

!     TEST TO SEE IF AN INACTIVE-ACTIVE ROW CONTRIBUTION SHOULD BE
!     ADDED IN

kkk = 0
690 in3 = i6sp + kkk
in1 = ix(in3) + i4sp
IF (ix(in1) > jpos+bbar) GO TO 720
kk  = in1 + cbar
IF (ix(kk) > jposl+ll+1) GO TO 710
IF (ix(in1)-jposl-bbar1 <= ll) GO TO 710

!     ADD IN EFFECT OF THE INACTIVE-ACTIVE ROW

in4 = i2 + 2*(ix(in1)-ioff)
k   = i4 + 2*(jposl+bbbar - jpos+ll + ix(in3)*bbbar)
dx(in4  ) = dx(in4  ) - dx(k)*dx(in2)   + dx(k+1)*dx(in2+1)
dx(in4+1) = dx(in4+1) - dx(in2+1)*dx(k) - dx(in2)*dx(k+1)
710 kkk = kkk + 1
IF (kkk < cbcnt) GO TO 690
720 ll  = ll  + 1
lll = lll + 1
IF (ll == lcol) GO TO 780
IF (ll-r+1 < 0) THEN
  GO TO   640
ELSE IF (ll-r+1 == 0) THEN
  GO TO   730
ELSE
  GO TO   760
END IF
730 IF (r == bbbar1) GO TO 640
in1 = i1 + 2*ll*bbar
750 icrq = in1 + bbar*4 - 1 - sr3buf
IF (icrq > 0) GO TO 1715
ibbar4 = bbar*4
CALL READ (*1690,*1700,sr2fl,dx(in1),ibbar4,0,flag)
GO TO 640
760 in1 = i1 + (lll-1)*bbar*2
IF (ll == r .AND. lcol == bbbar1) GO TO 770
CALL WRITE (sr3fl,dx(in1),4*bbar,0)
770 lll = lll - 1
GO TO 750
780 CONTINUE

!     COMPUTE ELEMENTS FOR THE ACTIVE ROWS

IF (cbcnt == 0) GO TO 830
k   = 0
790 in1 = i4sp + k
IF (ix(in1) > jpos+bbar) GO TO 810
800 k = k + 1
IF (k-cbar < 0) THEN
  GO TO   790
ELSE IF (k-cbar == 0) THEN
  GO TO   830
ELSE
  GO TO  1710
END IF
810 in1 = in1 + cbar
IF (ix(in1) == jpos) GO TO 800
kkk = MAX0(0,bbbar-jpos+ix(in1)-1)
in2 = i4  + 2*k*bbbar - 2
in3 = i2  + 2*(kkk-1-MAX0(0,bbbar-jpos))
in1 = in2 + 2*bbbar
in2 = in2 + 2*kkk
820 in2 = in2 + 2
kkk = kkk + 1
in3 = in3 + 2
dx(in1  ) = dx(in1  ) - dx(in2)*dx(in3)   + dx(in2+1) *dx(in3+1)
dx(in1+1) = dx(in1+1) - dx(in2+1)*dx(in3) - dx(in2)*dx(in3+1)
IF (kkk-bbbar1 < 0) THEN
  GO TO   820
ELSE IF (kkk-bbbar1 == 0) THEN
  GO TO   800
ELSE
  GO TO  1710
END IF

!     SEARCH THE LOWER BAND FOR THE MAXIMUM ELEMENT AND INTERCHANGE
!     ROWS TO BRING IT TO THE DIAGONAL

830 k   = 1
in1 = i2 + (jpos-ioff)*2
dx1 = 0.d0
dx2 = 0.d0
IF (DABS(dx(in1  )) > limit) dx1 = dx(in1  )**2
IF (DABS(dx(in1+1)) > limit) dx2 = dx(in1+1)**2
MAX(1) = dx1 + dx2
intchn = 0
END = MIN0(bbar1,ncol-jpos+1)
IF (END == 1) GO TO 870
840 in2 = in1 + k + k
dx1 = 0.d0
dx2 = 0.d0
IF (DABS(dx(in2  )) > limit) dx1 = dx(in2  )**2
IF (DABS(dx(in2+1)) > limit) dx2 = dx(in2+1)**2
dx2 = dx2 + dx1
IF (dx2 > MAX(1)) GO TO 860
850 k = k + 1
IF (k-END < 0) THEN
GO TO   840
ELSE IF (k-END == 0) THEN
GO TO   870
ELSE
GO TO  1710
END IF
860 MAX(1) = dx2
intchn = k
GO TO 850

870 IF (intchn == 0) GO TO 880

!     INTERCHANGE ROWS IN AREA II

det(1)  =-det(1)
det(2)  =-det(2)

MAX(1)  = dx(in1)
in2     = in1+2*intchn
dx(in1) = dx(in2)
dx(in2) = MAX(1)
MAX(1)  = dx(in1+1)
dx(in1+1) = dx(in2+1)
dx(in2+1) = MAX(1)

!     STORE THE PERMUTATION INDEX

in2 = i1sp + lcol
ix(in2) = intchn

!     DIVIDE THE LOWER BAND BY THE DIAGONAL ELEMENT

880 dx1 = 0.d0
dx2 = 0.d0
IF (DABS(dx(in1  )) > limit) dx1 = dx(in1  )**2
IF (DABS(dx(in1+1)) > limit) dx2 = dx(in1+1)**2
da(1) = dx1 + dx2
IF (da(1) == 0.d0) GO TO 1720
MAX(1) = dx(in1  )/da(1)
MAX(2) =-dx(in1+1)/da(1)
mindia = DMIN1(DSQRT(da(1)),mindia)
da(1)  = DMAX1(DABS(det(1)),DABS(det(2)))
890 IF (da(1) <= 10.d0) GO TO 900
det(1) = det(1)*.1D0
det(2) = det(2)*.1D0
da(1)  = da(1) *.1D0
power  = power + 1
GO TO 890
900 IF (da(1) >= .1D0) GO TO 910
det(1) = det(1)*10.d0
det(2) = det(2)*10.d0
da(1)  = da(1) *10.d0
power  = power - 1
GO TO 900
910 da(1)  = det(1)*dx(in1) - det(2)*dx(in1+1)
det(2) = det(2)*dx(in1) + det(1)*dx(in1+1)
det(1) = da(1)
k   = 1
END = MIN0(bbar1,ncol-jpos+1)
IF (END == 1) GO TO 930
920 in2 = in1 + k + k
da(1)     = dx(in2)*MAX(1) - dx(in2+1)*MAX(2)
dx(in2+1) = dx(in2)*MAX(2) + dx(in2+1)*MAX(1)
dx(in2) = da(1)
k = k + 1
IF (k-END < 0) THEN
GO TO   920
ELSE IF (k-END == 0) THEN
GO TO   930
ELSE
GO TO  1710
END IF
930 IF (cbcnt == 0) GO TO 950

!     DIVIDE THE ACTIVE ROWS BY THE DIAGONAL

k   = 0
in1 = i4 + 2*bbbar1
940 da(1)     = dx(in1)*MAX(1) - dx(in1+1)*MAX(2)
dx(in1+1) = dx(in1)*MAX(2) + dx(in1+1)*MAX(1)
dx(in1) = da(1)
in1 = in1 + 2*bbbar
k   = k + 1
IF (k-cbar < 0) THEN
  GO TO   940
ELSE IF (k-cbar == 0) THEN
  GO TO   950
ELSE
  GO TO  1710
END IF
950 CONTINUE

!     INTERCHANGE ACTIVE COLUMNS AND ADD IN EFFECT OF THE CURRENT COLUMN

IF (ccount ==   0) GO TO 1000
IF (jpos < bbbar) GO TO 1000
intch = ix(i1sp)
k   = 0
960 in1 = i3sp + k
IF (intch == 0) GO TO 970
in1 = i3  + 2*k*bbar1
in2 = in1 + intch + intch
da(1)   = dx(in1)
dx(in1) = dx(in2)
dx(in2) = da(1)
da(1)   = dx(in1+1)
dx(in1+1) = dx(in2+1)
dx(in2+1) = da(1)
970 kk  = 1
in2 = i1 - 2
in1 = i3 + 2*k*bbar1
IF (dx(in1) == 0.d0 .AND. dx(in1+1) == 0.d0) GO TO 990
980 in3 = in1 + kk + kk
in4 = in2 + kk + kk
dx(in3  ) = dx(in3  ) - dx(in1)*dx(in4  ) + dx(in1+1)*dx(in4+1)
dx(in3+1) = dx(in3+1) - dx(in1)*dx(in4+1) - dx(in1+1)*dx(in4)
kk = kk + 1
IF (kk-bbar1 < 0) THEN
  GO TO   980
ELSE IF (kk-bbar1 == 0) THEN
  GO TO   990
ELSE
  GO TO  1710
END IF
990 k = k + 1
IF (k-c < 0) THEN
  GO TO   960
ELSE IF (k-c == 0) THEN
  GO TO  1000
ELSE
  GO TO  1710
END IF

!     WRITE OUT THE NEXT COLUMN OF U AND THE ROW OF ACTIVE ELEMENTS

1000 parm(2) = sr2fil
CALL bldpk (cdp,typel,sr2fil,0,0)
in1 = i2
jj  = ioff
1010 dz(1) = dx(in1)
dz(2) = dx(in1+1)
IF (dz(1) == 0.d0 .AND. dz(2) == 0.d0) GO TO 1030
CALL zblpki
1030 in1 = in1 + 2
jj  = jj  + 1
IF (jj-jpos > 0) THEN
  GO TO  1040
ELSE
  GO TO  1010
END IF
1040 IF (dx(in1-2) == 0.d0 .AND. dx(in1-1) == 0.d0) GO TO 1720

!     PACK ACTIVE COLUMN ELEMENTS ALSO

IF (ccount ==   0) GO TO 1090
IF (jpos < bbbar) GO TO 1090
k   = 0
1060 in1 = i7sp + k
in2 = ix(in1) + i3sp
GO TO 1080
1070 k   = k + 1
IF (k-ccount < 0) THEN
  GO TO  1060
ELSE IF (k-ccount == 0) THEN
  GO TO  1090
ELSE
  GO TO  1710
END IF
1080 in3 = i3 + 2*(ix(in1)*bbar1)
dz(1) = dx(in3  )
dz(2) = dx(in3+1)
IF (dz(1) == 0.d0 .AND.  dz(2) == 0.d0) GO TO 1070
jj = ix(in2)
CALL zblpki
GO TO 1070
1090 CALL bldpkn (sr2fil,0,fileu)

!     COMPUTE ACTIVE ROW-COLUMN INTERACTION

IF (ccount == 0 .OR. cbcnt == 0) GO TO 1140
IF (jpos < bbbar) GO TO 1140
k = 0
1100 CONTINUE
in1 = i3 + 2*k*bbar1
IF (dx(in1) == 0.d0 .AND. dx(in1+1) == 0.d0) GO TO 1130
kk  = 0
1110 in2 = i4 + 2*kk*bbbar
IF (dx(in2) == 0.d0 .AND. dx(in2+1) == 0.d0) GO TO 1120
in3 = i5 + 2*(k*cbar+kk)
dx(in3  ) = dx(in3  ) + dx(in2)*dx(in1)   - dx(in2+1)*dx(in1+1)
dx(in3+1) = dx(in3+1) + dx(in2)*dx(in1+1) + dx(in2+1)*dx(in1)
1120 kk = kk + 1
IF (kk-cbar < 0) THEN
  GO TO  1110
ELSE IF (kk-cbar == 0) THEN
  GO TO  1130
ELSE
  GO TO  1710
END IF
1130 k = k + 1
IF (k-c < 0) THEN
  GO TO  1100
ELSE IF (k-c == 0) THEN
  GO TO  1140
ELSE
  GO TO  1710
END IF

!     MOVE ELEMENTS IN AREA III UP ONE CELL

1140 IF (ccount ==   0) GO TO 1190
IF (jpos < bbbar) GO TO 1190
k   = 0
1150 in1 = i3sp + k
IF (ix(in1) == 0) GO TO 1180
kk  = 0
in1 = i3  + 2*k*bbar1
1160 in2 = in1 + kk + kk
dx(in2  ) = dx(in2+2)
dx(in2+1) = dx(in2+3)
kk = kk + 1
IF (kk-bbar < 0) THEN
  GO TO  1160
ELSE IF (kk-bbar == 0) THEN
  GO TO  1170
ELSE
  GO TO  1710
END IF
1170 dx(in2+2) = 0.d0
dx(in2+3) = 0.d0
1180 k = k + 1
IF (k-c < 0) THEN
  GO TO  1150
ELSE IF (k-c == 0) THEN
  GO TO  1190
ELSE
  GO TO  1710
END IF


!     DETERMINE IF A COLUMN OF L CAN BE WRITTEN OUT

1190 IF (lcol-bbbar1 < 0) THEN
  GO TO  1370
END IF

!     OUTPUT A COLUMN OF L

1200 parm(2) = filel(1)
jposl   = jposl + 1
CALL bldpk (cdp,typel,filel(1),0,0)

!     STORE THE PERMUTATION INDEX AS THE DIAGONAL ELEMENT

jj    = jposl
dz(1) = ix(i1sp)
dz(2) = 0.d0
CALL zblpki
k  = 0
1210 jj = jposl + k + 1
in2 = i1 + k + k
dz(1) = dx(in2  )
dz(2) = dx(in2+1)
IF (dz(1) == 0.d0 .AND. dz(2) == 0.d0) GO TO 1230
CALL zblpki
1230 k = k + 1
IF (k-bbar < 0) THEN
  GO TO  1210
ELSE IF (k-bbar == 0) THEN
  GO TO  1240
ELSE
  GO TO  1710
END IF

!     PACK ACTIVE ROW ELEMENTS ALSO

1240 IF (cbcnt == 0) GO TO 1280
k   = 0
1250 in1 = i6sp + k
in2 = i4 + (ix(in1)*bbbar)*2
in1 = ix(in1) + i4sp
jj  = ix(in1)
dz(1) = dx(in2  )
dz(2) = dx(in2+1)
IF (dz(1) == 0.d0 .AND. dz(2) == 0.d0) GO TO 1270
CALL zblpki
1270 k = k + 1
IF (k-cbcnt < 0) THEN
  GO TO  1250
ELSE IF (k-cbcnt == 0) THEN
  GO TO  1280
ELSE
  GO TO  1710
END IF
1280 CALL bldpkn (filel,0,filel)

!     MOVE PERMUTATION INDICES OVER ONE ELEMENT

END = i1sp + lcol
DO  i = i1sp,END
ix(i) = ix(i+1)
END DO

!     MOVE ELEMENTS IN AREA I OVER ONE COLUMN

k = 0
IF (scrflg == 0) GO TO 1310
CALL CLOSE (sr2fl,rew)
CALL OPEN (*1680,sr2fl,ix(sr2buf),rd)
IF (r > 2) GO TO 1310
icrq = i1 + bbar*4 - 1 - sr3buf
IF (icrq > 0) GO TO 1715
ibbar4 = bbar*4
CALL READ (*1690,*1700,sr2fl,dx(i1),ibbar4,0,flag)
GO TO 1360
1310 in1 = i1  + k*bbar*2
in2 = in1 + bbar + bbar
CALL cxloop (dx(in1),dx(in2),bbar)
k = k + 1
IF (k-r+2 < 0) THEN
  GO TO  1310
ELSE IF (k-r+2 == 0) THEN
  GO TO  1340
ELSE
  GO TO  1360
END IF
1340 IF (r-bbbar1 < 0.0) THEN
  GO TO  1350
ELSE IF (r-bbbar1 == 0.0) THEN
  GO TO  1310
ELSE
  GO TO  1710
END IF
1350 icrq = in2 + bbar*4 - 1 - sr3buf
IF (icrq > 0) GO TO 1715
ibbar4 = bbar*4
CALL READ (*1690,*1700,sr2fl,dx(in2),ibbar4,0,flag)
1360 lcol = lcol - 1

!     STORE CURRENT COLUMN OF L

1370 IF (cbcnt == 0) GO TO 1420

!     MOVE ELEMENTS IN AREA IV UP ONE CELL

k   = 0
1380 in1 = i4sp + k
IF (ix(in1) == 0) GO TO 1410
kk  = 0
in1 = i4  + 2*k*bbbar
1390 in2 = in1 + kk + kk
dx(in2  ) = dx(in2+2)
dx(in2+1) = dx(in2+3)
kk = kk + 1
IF (kk-bbbar1 < 0) THEN
  GO TO  1390
ELSE IF (kk-bbbar1 == 0) THEN
  GO TO  1400
ELSE
  GO TO  1710
END IF
1400 dx(in2+2) = 0.d0
dx(in2+3) = 0.d0
1410 k = k + 1
IF (k-cbar < 0) THEN
  GO TO  1380
ELSE IF (k-cbar == 0) THEN
  GO TO  1420
ELSE
  GO TO  1710
END IF
1420 IF (scrflg /= 0) GO TO 1450

!     STORE COLUMN IN CORE

1430 in1 = i1 + 2*lcol*bbar
END = MIN0(bbar,ncol-jpos)
IF (END == 0) GO TO 1480
k   = 0
in3 = i2  + 2*(jpos-ioff+1)
1440 in2 = in1 + k + k
in4 = in3 + k + k
dx(in2  ) = dx(in4  )
dx(in2+1) = dx(in4+1)
k = k + 1
IF (k-END < 0) THEN
GO TO  1440
ELSE IF (k-END == 0) THEN
GO TO  1480
ELSE
GO TO  1710
END IF

!     STORE COLUMN ON THE SCRATCH FILE

1450 IF (lcol-r+1 < 0) THEN
  GO TO  1430
ELSE IF (lcol-r+1 == 0) THEN
  GO TO  1470
END IF
1460 in1 = i1 + (lll-1)*bbar*2
CALL WRITE (sr3fl,dx(in1),bbar*4,0)
1470 in1 = i2 + 2*(jpos-ioff+1)
CALL WRITE (sr3fl,dx(in1),bbar*4,0)

!     CLOSE SCRATCH FILES AND SWITCH THE POINTERS TO THEM

CALL CLOSE (sr3fl,rew)
CALL CLOSE (sr2fl,rew)
in1   = sr2fl
sr2fl = sr3fl
sr3fl = in1
1480 lcol  = lcol + 1
IF (c == 0) GO TO 1570
IF (jpos < bbbar) GO TO 1570

!     READ IN THE NEXT ROW OF ACTIVE COLUMN ELEMENTS

count = ccount
IF (itrn < 0) GO TO 1570
1490 IF (itrn > jpos-b+2) GO TO 1560

!     TEST TO SEE IF COLUMN IS ALREADY ACTIVE

k   = 0
1500 in1 = i3sp + k
IF (ix(in1) == jtrn) GO TO 1540
k   = k + 1
IF (k-c < 0) THEN
  GO TO  1500
ELSE IF (k-c == 0) THEN
  GO TO  1510
ELSE
  GO TO  1710
END IF

!     CREATE A NEW ACTIVE COLUMN

1510 k   = 0
1520 in1 = i3sp + k
IF (ix(in1) == 0) GO TO 1530
k   = k + 1
IF (k-c < 0) THEN
  GO TO  1520
ELSE
  GO TO  1710
END IF
1530 ix(in1) = jtrn
in1 = in1 + c
ix(in1) = itrn
in1 = i3 + 2*(k+1)*bbar1 - 2
dx(in1  ) = dtrn(1)
dx(in1+1) = dtrn(2)
ccount = ccount + 1
GO TO 1550

!     STORE ELEMENT IN EXISTING COLUMN

1540 in1 = i3 + 2*(k+1)*bbar1 - 2
dx(in1  ) = dx(in1  ) + dtrn(1)
dx(in1+1) = dx(in1+1) + dtrn(2)
1550 CALL READ (*1690,*1700,sr1fil,itran(1),6,0,flag)
IF (itrn > 0) GO TO 1490
CALL CLOSE (sr1fil,rew)
1560 IF (ccount == count) GO TO 1570

!     RE-ARRANGE INDEXES IN SEQUENTIAL ORDER

ASSIGN 1570 TO kk
GO TO 150
1570 CONTINUE
jpos = jpos + 1

!     ZERO AREA II

END = i2 + 2*MIN0(jpos-ioff+bbar-1,ncol-1) + 1
DO  i = i2,END
dx(i) = 0.d0
END DO

!      TEST TO SEE IF ROW INTERACTION ELEMENTS WILL MERGE INTO AREA III

IF (cbcnt  == 0) GO TO 270
IF (ccount == 0) GO TO 1640
IF (jpos-1 < bbbar) GO TO 270
in1 = i4sp
k   = 0
1600 in2 = in1 + k
IF (ix(in2) == jpos-b+1) GO TO 1610
k   = k + 1
IF (k < cbar) GO TO 1600
GO TO 270
1610 in1 = i5 + k + k
in2 = i3 + bbar + bbar
k   = 0
1620 dx(in2  ) = dx(in2  ) - dx(in1)
dx(in2+1) = dx(in2+1) - dx(in1+1)
dx(in1  ) = 0.d0
dx(in1+1) = 0.d0
in2 = in2 + bbar1 + bbar1
in1 = in1 + cbar  + cbar
k   = k + 1
IF (k < c) GO TO 1620

!     TEST TO SEE IF A ACTIVE ROW HAS BEEN ELIMINATED

1640 in1 = ix(i6sp) + i4sp
IF (ix(in1)-jposl-bbar1 == 0) THEN
  GO TO  1650
ELSE
  GO TO   270
END IF

!     ELIMINATE THE ACTIVE ROW

1650 ix(in1) = 0
in1     = in1 + cbar
ix(in1) = 0
cbcnt   = cbcnt - 1

!     MOVE INDEXES IN AREA VI UP ONE

in1 = i6sp + cbcnt - 1
DO  i = i6sp,in1
  ix(i) = ix(i+1)
END DO
ix(in1+1) = 0
GO TO 270

!     FINISH WRITING OUT THE COMPLETED COLUMNS OF L

1670 CALL CLOSE (sr1fil,rew)
CALL CLOSE (filel,norew)
CALL CLOSE (sr2fil,norew)
CALL comfin (iterm,scrflg,sr2fl,jposl,i1sp,bbar,i1,cbcnt,ipak,r,  &
    bbbar1,bbbar,i6sp,i4,i4sp,ix,dx,x,lcol)
parm(5) = iend
CALL conmsg (parm(3),3,0)
fileu(7) = bbbar
RETURN

!     ERROR EXITS

1680 parm(1) = -1
GO TO 1730
1690 parm(1) = -2
GO TO 1730
1700 parm(1) = -3
GO TO 1730
1710 parm(1) = -25
GO TO 1730
1715 parm(1) = -8
parm(2) = icrq
GO TO 1730

!     SINGULAR MATRIX - CLOSE ALL FILES AND RETURN TO USER

1720 CALL CLOSE (filea(1),rew)
CALL CLOSE (filel(1),rew)
CALL CLOSE (fileu(1),rew)
CALL CLOSE (sr1fil,rew)
CALL CLOSE (sr2fil,rew)
CALL CLOSE (sr3fil,rew)
!WKBR SPR94018 4/95      FILEU(2) = BBBAR
fileu(7) = bbbar
RETURN 1

1730 CALL mesage (parm(1),parm(2),parm(3))
RETURN
END SUBROUTINE cdcomp
