SUBROUTINE opt2a (ip,el,iel,pr,ipr,rr)
     
 
 INTEGER, INTENT(IN)                      :: ip(2,1)
 REAL, INTENT(IN)                         :: el(1)
 INTEGER, INTENT(IN)                      :: iel(1)
 REAL, INTENT(OUT)                        :: pr(1)
 INTEGER, INTENT(IN)                      :: ipr(1)
 REAL, INTENT(OUT)                        :: rr(1)
 LOGICAL :: first,unsafe
 INTEGER :: count,etyp, iz(10),NAME(2),  &
     oes1,outtap,pest,pstres,ptelt,zcor,oldtyp,eid(20), plus(5),iy(1)
 REAL :: y(1),parm(8)
 CHARACTER (LEN=29) :: uim
 CHARACTER (LEN=25) :: uwm
 CHARACTER (LEN=23) :: ufm
 COMMON /xmssg / ufm,uwm,uim
 COMMON /BLANK / skp(2),count,skq(2),kore,skr(2),nwdse,nwdsp,sks,  &
     oes1,skt(3),nelw,nprw,sku,ntotl,conv
 COMMON /optpw2/ zcor,z(16)
 COMMON /zzzzzz/ core(1)
 COMMON /names / nrd,noeor,nwrt,next
 COMMON /system/ sysbuf,outtap
!     EQUIVALENT ARE  (EL,IEL), (PR,IPR)
 EQUIVALENCE     (z(1),iz(1)), (core(1),parm(1),MAX), (iy(1),y(1),parm(8))
 DATA    NAME  / 4H opt,4H2A   /
 DATA    plus  / 4H    , 4H+   , 4H++  , 4H+++ , 4H++++ /
 
 nelr  = 0
 NE    = 0
 ptelt = 0
 idel  = 0
 kel   = kore
 kconv = 0
 conv  = 1.0
 icp   = ntotl - 4
 first =.true.
 
!     READ HEADER, ODD RECORDS
 
 GO TO 10
 5 CALL fread (oes1,0,0,next)
 10 CALL READ  (*630,*100,oes1,z(1),10,next,i)
 etyp  = iz(3)
 nesw  = iz(10)
 oldtyp= ptelt
 ptelt = iy(etyp)
 IF (ptelt > 0) GO TO 15
 
!     ELEMENT TYPE NOT TO OPTIMIZE
 
 GO TO 5
 15 IF (ptelt >= oldtyp .OR. oldtyp == 0) GO TO 20
 IF (kel /= -1) kel = kore
 IF (NE  ==  0) GO TO 16
 CALL page2 (1)
 WRITE (outtap,580) (eid(j),j=1,NE)
 NE = 0
 16 WRITE  (outtap,17)
 17 FORMAT (/5X,15HNEXT subcase...)
 
!     SET POINTERS TO ELEMENT TYPE AND PROPERTIES IN CORE.
!     L = LOCATION OF FIRST, M = MAX LOCATION
 
 20 lel = ip(1,ptelt)
 mel = ip(1,ptelt+1) - 1
 IF (mel <= lel) GO TO 5
 loce  = lel
 locp1 = ip(2,ptelt) - 1
 IF (nesw > zcor) GO TO 70
 
!     SEQUENTIALLY READ ONE ELEMENT FROM EVEN NUMBERED RECORDS.
!     LOCE IS CURRENT ELEMENT TO COMPARE TO.
 
 30 CALL READ (*90,*10,oes1,z(1),nesw,noeor,i)
 ides = iz(1)/10
 50 IF (ides == iel(loce)) GO TO 110
 
!     SCAN THE CORE FILE UNTIL ELEMENT ID .GT. IDES
 
 IF (ides < iel(loce)) GO TO 30
 
!     CORE ELEMENT NOT TO BE OPTIMIZED
 
 loce = loce + nwdse
 IF (loce < mel) GO TO 50
 
!     END OF ELEMENT SEARCH FOR THIS TYPE (EOR NOT READ)
 
 GO TO 5
 
!     ELEMENT TYPE EXCEEDS CORE
 
 70 ier = -8
 ifle = nesw - zcor
 GO TO 105
 
!     ILLEGAL EOF, EOR
 
 90 ier = -2
 GO TO 101
 100 ier = -3
 101 ifle = oes1
 
 105 CALL mesage (ier,ifle,NAME)
 
!     PROCES THIS ELEMENT
 
 110 CONTINUE
 nelr = nelr + 1
 locp = iel(loce+4) + locp1
 pest = ipr(locp+1)/100
 mest = ipr(locp+1) - pest*100
 rc   = 1.0
 x1a  = 0.0
 x2a  = 0.0
 e1   = 999.
 unsafe = .false.
 
 SELECT CASE ( ptelt )
   CASE (    1)
     GO TO 160
   CASE (    2)
     GO TO 160
   CASE (    3)
     GO TO 180
   CASE (    4)
     GO TO 150
   CASE (    5)
     GO TO 150
   CASE (    6)
     GO TO 150
   CASE (    7)
     GO TO 140
   CASE (    8)
     GO TO 140
   CASE (    9)
     GO TO 140
   CASE (   10)
     GO TO 120
   CASE (   11)
     GO TO 130
   CASE (   12)
     GO TO 140
   CASE (   13)
     GO TO 140
   CASE (   14)
     GO TO 140
   CASE (   15)
     GO TO 170
   CASE (   16)
     GO TO 150
   CASE (   17)
     GO TO 140
   CASE (   18)
     GO TO 120
   CASE (   19)
     GO TO 140
   CASE (   20)
     GO TO 140
 END SELECT
 
!     ROD, TUBE
 
 120 limit  = 1
 pstres = 4
 ASSIGN 121 TO iret
 GO TO 500
 121 limit  = 2
 pstres = 2
 ASSIGN  540 TO iret
 GO TO 500
 
!     SHEAR
 
 130 limit  = 1
 pstres = 2
 ASSIGN  540 TO iret
 GO TO 500
 
!     TRBSC, TRPLT, QDPLT, TRIA1, TRIA2, TRIA3, QUAD1, QUAD2, QUAD4
 
 140 IF (mest == 1) GO TO 144
 limit  = 2
 pstres = 7
 ASSIGN 141 TO iret
 GO TO 500
 141 pstres = 8
 ASSIGN 142 TO iret
 GO TO 500
 142 pstres = 15
 ASSIGN 143 TO iret
 GO TO 500
 143 pstres = 16
 ASSIGN 144 TO iret
 x1a = AMAX1(ABS(z( 7)),ABS(z( 8)))
 x2a = AMAX1(ABS(z(15)),ABS(z(16)))
 x1a = AMAX1(x1a,x2a)
 k   = 0
 IF (x1a == ABS(z(8)) .OR. x1a == ABS(z(15))) k = 1
 x1a = z( 7+k)
 x2a = z(16-k)
 GO TO 500
 144 IF (mest == 2) GO TO 540
 limit  = 1
 pstres = 9
 ASSIGN 145 TO iret
 GO TO 500
 145 pstres = 17
 ASSIGN  540 TO iret
 GO TO 500
 
!     TRMEM, QDMEM, QDMEM1, QDMEM2
 
 150 IF (mest == 1) GO TO 152
 limit  = 2
 pstres = 6
 ASSIGN 151 TO iret
 GO TO 500
 151 pstres = 7
 ASSIGN 152 TO iret
 GO TO 500
 152 IF (mest == 2) GO TO 30
 limit  = 1
 pstres = 8
 ASSIGN  540 TO iret
 GO TO 500
 
!     BAR, ELBOW
 
 160 limit  = 2
 pstres = 7
 x2a = ABS(z(7))
 ASSIGN 161 TO iret
 GO TO 500
 161 pstres = 8
 x1a = ABS(z(8))
 ASSIGN 162 TO iret
 GO TO 500
 162 pstres = 14
 ASSIGN 163 TO iret
 GO TO 500
 163 pstres = 15
 ASSIGN  540 TO iret
 GO TO 500
 
!     TRIM6
 
 170 IF (iel(loce) == idel) GO TO 172
 idel = iel(loce)
 icp  = icp + 4
 IF (kel /= -1 .AND. icp >= kel) CALL mesage (-8,0,NAME)
 iy(icp) = locp
 iy(icp+4) =-1
 172 k  = 0
 m1 =-1
 DO  i = 1,3
   m1 = m1 + 7
   ii = 3 + loce
   s1s = 0.0
   s3s = 0.0
   IF (mest /= 2) s3s = ABS(z(m1+2)/el(ii))
   ii = ii - 2
   IF (z(m1) < 0.0) ii = ii + 1
   IF (mest /= 1) s1s = ABS(z(m1)/el(ii))
   ii = 1 + loce
   IF (z(m1+1) < 0.0) ii = ii + 1
   s2s = ABS(z(m1+1)/el(ii))
   s13 = AMAX1(s1s,s2s)
   s13 = AMAX1(s13,s3s)
   y(icp+i) = AMAX1(y(icp+i),s13)
   pr(locp+4) = AMAX1(pr(locp+4),s13)
   e1 = ABS(s13) - 1.0
   IF (ABS(e1) <= parm(2)) k = k + 1
 END DO
 ASSIGN 540 TO iret
 IF (k-3 < 0) THEN
   GO TO   550
 ELSE
   GO TO   520
 END IF
 
!     IS2D8
 
 180 m1  = 1
 s1s = 0.0
 s2s = 0.0
 s3s = 0.0
 DO  m = 1,8
   m1 = m1 + 5
   ii = 3 + loce
   IF (mest /= 2) s3s = AMAX1(s3s,ABS(z(m1+2)/el(ii)))
   ii = ii - 2
   IF (z(m1) < 0.0) ii = ii + 1
   IF (mest /= 1) s1s = AMAX1(s1s,ABS(z(m1)/el(ii)))
   ii = 1 + loce
   IF (z(m1+1) < 0.0) ii = ii + 1
   s2s = AMAX1(s2s,ABS(z(m1+1)/el(ii)))
   s13 = AMAX1(s1s,s2s)
   s13 = AMAX1(s13,s3s)
 END DO
 e1 = ABS(s13) - 1.0
 pr(locp+4) = AMAX1(pr(locp+4),s13)
 ASSIGN 540 TO iret
 GO TO 520
 
!     FUNCTION E1  -  RATIO STRESS MINUS LIMIT DIVIDED BY LIMIT,
!     WITH RESET OF -ALPHA-
!     LOCP   = POINTER TO PID OF PROPERTY.
!     LOCE   = POINTER TO EID OF ELEMENT.
!     LIMIT  = 1=SHEAR, 2= COMPRESSION/TENSION.
!     PSTRES = CORRESPONDING STRESS, POINTER TO Z ARRAY.
 
 500 ii = 3 + loce
 IF (limit == 1) GO TO 510
 ii = ii - 2
 IF (z(pstres) < 0.0) ii = ii + 1
 510 IF (el(ii) <= 0.0) GO TO 530
 
!     POSITIVE LIMIT
 
 pr(locp+4) = AMAX1(pr(locp+4),ABS(z(pstres)/el(ii)))
 
!                                        I
!                  NEGATIVE E1, SAFE     I    POSITIVE E1, UNSAFE
!                                        I
!   --+------+------+------+------+------+------+------------------- E1
!    UL     4P     3P     2P      P      0      P  (WHERE P=PARM(2),
!      ++++    +++    ++     +    I             I        UL=UNLOADED)
!            OVER DESIGNED        I REGION WHEREI  UNDER DESIGNED
!            REGION               I  AE1 .LE. P I          REGION
!                      (UNSAFE=.FALSE.)         I  (UNSAFE=.TRUE.)
 
 e1 = ABS(z(pstres)/el(ii)) - 1.0
 520 IF (e1 > parm(2)) unsafe = .true.
 IF (unsafe) kel = -1
 ae1 = AMIN1(ae1,ABS(e1))
 530 GO TO iret, (121,141,142,143,144,145,151,152,161,162,163,540)
 
 540 x1 = ABS(x1a)
 x2 = ABS(x2a)
 IF (x1 == 0.0 .OR. x2 == 0.0) GO TO 550
 x1a= AMIN1(x1a,x2a)
 x1 = AMIN1(x1,x2)/AMAX1(x1,x2)
 x1 = SIGN(x1,x1a)
 IF (ABS(x1) > 1.0E-8) rc = x1
 
!     SAVE IN RR AN EMPIRICAL ALPHA MODIFIER FOR SPEEDY CONVERGENCE
 
 550 irr = (locp+nwdsp)/nwdsp
 rr(irr) = rc
 
 IF (unsafe) GO TO 30
 
!     PRINT ELEMENT IDS THAT HAVE CONVERGED, OR OVER DESIGNED
 
 IF (.NOT.first) GO TO 570
 first = .false.
 CALL page2 (-3)
 WRITE  (outtap,560) uim
 560 FORMAT (a29,' 2304A, THE FOLLOWING ELEMENTS EITHER CONVERGED (NO',  &
     ' PLUS) OR OVER-DESIGNED (PLUS(ES))',/5X,'IN ONE OR MORE ',  &
     'SUBCASES,  (EACH PLUS INDICATES AN INCREMENTAL PERCENTAGE'  &
     ,      ' OF OVER-DESIGN BASED ON CONVERGENCE CRITERION, EPS)',/)
 570 xstar = (pr(locp+4)-1.0) - parm(2)
 j  = IFIX(ABS(xstar)/parm(2))
 IF (j > 3) j = 3
 ii = 1
 IF (pr(locp+4) < 1.0E-8) ii = 0
 IF (ii == 0) j = 4
 eid(NE+1) = iel(loce)
 eid(NE+2) = plus(j+1)
 NE = NE + 2
 IF (NE < 20) GO TO 590
 NE = 0
 CALL page2 (1)
 WRITE  (outtap,580) eid
 580 FORMAT (5X,10(i8,a4))
 590 IF (kel == -1) GO TO 30
 kel = kel - 1
!WKBR 9/93 IZK = IZ(KEL)
 izk = iy(kel)
 IF (pr(locp+3) < 1.0E-6) ii = 0
 IF (j > 0 .AND. izk == -1 .AND. ii /= 0) kconv = kconv - 1
 IF (ii == 0) GO TO 600
 IF (ae1 > parm(2)) GO TO 30
 600 IF (iel(loce) == izk) GO TO 30
 IF (ae1 <= parm(2) .AND. izk == -1) GO TO 610
 IF (ii == 0 .AND. izk == -1) GO TO 30
!WKBR 9/93 IZ(KEL) = IEL(LOCE)
 iy(kel) = iel(loce)
!WKBR 9/93 IF (II .EQ. 0) IZ(KEL) = -1
 IF (ii == 0) iy(kel) = -1
 kconv = kconv + 1
 GO TO 30
!WKBR 9/93  610 IZ(KEL) = IEL(LOCE)
 610 iy(kel) = iel(loce)
 GO TO 30
 
!     EOF
 
 630 CONTINUE
 IF (NE > 0) WRITE (outtap,580) (eid(j),j=1,NE)
 
!     IF KEL=-1 HERE, OR
!     IF NUMBER OF ELEMENTS CONVERGED, KORE-KEL, IS LESS THAN NUMBER OF
!     ELEMENTS IN THE PROBLEM, NELW/NWDSE, CONVERGENCE IS INCOMPLETE
 
 IF (kel == -1) GO TO 650
 IF (kconv < nelw/nwdse) GO TO 650
!WKBR CALL PAGE (-4)
 CALL page2 (-4)
 WRITE  (outtap,640) uim
 640 FORMAT (a29,' 2304B, CONVERGENCE ACHIEVED FOR ALL ELEMENTS ',  &
     'REQUESTED, AND IN ALL SUBCASE(S)', /5X, 'FULLY-STRESSED DESIGN COMPLETED',/)
 conv = 2.0
 GO TO 670
 
!     IF NELR IS ZERO, NO ELEMENT MATCH MADE
 
 650 IF (nelr > 0) GO TO 670
 CALL page2 (-2)
 WRITE  (outtap,660) ufm
 660 FORMAT (a23,' 2295, NO ELEMENTS EXIST FOR OPTIMIZATION.')
 count = MAX + 1
 
 670 RETURN
END SUBROUTINE opt2a
