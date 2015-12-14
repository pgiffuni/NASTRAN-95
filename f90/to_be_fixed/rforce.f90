SUBROUTINE rforce (lcore)
     
!     COMPUTES STATIC LOADS DUE TO ROTATING COORDINATE SYSTEMS
 
 
 INTEGER, INTENT(IN)                      :: lcore
 EXTERNAL        rshift,andf
 LOGICAL :: nonshl,cupmas
 INTEGER :: FILE,slt,bgpdt,OLD,icard(6),sysbuf,NAME(2), strtmn,andf,rshift
 REAL :: mt(3,3),mtr(3,3),mr(3,3)
 DIMENSION       card(6),ra(4),wb(3),wg(3),ri(4),xm(6,6),iy(7)
 DIMENSION       isystm(175)
 COMMON /machin/ mach,ihalf,jhalf
 COMMON /condas/ pi,twophi,radeg,degra,s4pisq
 COMMON /unpakx/ it1,ii,jj,incr
 COMMON /xcstm / ti(3,3)
 COMMON /tranx / ix(5),TO(3,3)
 COMMON /zzzzzz/ z(1)
 COMMON /system/ sysbuf,dumy(25),mn
 COMMON /zntpkx/ a(4),irow,ieol,IEOR
 COMMON /loadx / lc,slt,bgpdt,OLD,nn(11),mgg
 EQUIVALENCE     (icard(1),card(1)), (ir,ri(1)), (ira,ra(1))
 EQUIVALENCE     (sysbuf  ,isystm(1))
 DATA    NAME  / 4HRFOR,4HCE  /
 
!     DEFINITION OF VARIABLES
 
!     SLT      STATIC LOAD TABLE
!     BGPDT    BASIC GRID POINT DEFINITION TABLE
!     MGG      MASS  MATRIX
!     FILE     FILE NAME FOR ERROR MESAGES
!     CARD     CARD IMAGE OF RFORCE CARD
!     RA       BGPDT ENTRY FOR AXIAL GRID POINT
!     WB       OMEGA-S IN BASIC COORDINATES
!     II       SIL OF CURRENT  POINT
!     IT1      UNPACK TYPE(REAL)
!     INCR     INCREMENT( TO ROW STORE COLUMNS)
!     RI       BGPDT ENTRY FOR CURRENT GRID POINT
!     WG       OMEGA-S IN GLOBAL COORDINANTS AT CURRENT GRID POINT
!     XM       6X6 DIAGONAL PARTION OF MGG
!     MT       3X3 PARTITION OF  MGG
!     MR       3X3 PARTITION OF  MGG
!     MTR      3X3 PARTITION OF  MGG
!     OLD      CURRENT POSITION OF BGPDT  0 IMPLIES BEGINNING
 
 
!     BRING IN CARD IMAGE
 
 CALL fread (slt,card,6,0)
 
!     FIND LOCATION OF AXIAL GRID POINT
 
 DO  i = 1,3
   ra(i+1) = 0.0
 END DO
 IF (icard(1) == 0) GO TO 30
 CALL fndpnt (ra(1),icard(1))
 
!     CHECK FOR GRID POINT
 
 IF (ira /= -1) GO TO 30
 DO  i = 1,3
   ra(i+1) = 0.0
 END DO
 30 CALL REWIND (bgpdt)
 CALL skprec (bgpdt,1)
 
!     CONVERT WI'S TO BASIC COORDINANTS
 
 DO  i = 4,6
   wb(i-3) = card(i)*twophi*card(3)
 END DO
 IF (icard(2) == 0) GO TO 60
 CALL fdcstm (icard(2))
 CALL mpyl (TO,wb,3,3,1,wg)
 DO  i = 1,3
   wb(i) = wg(i)
 END DO
 
!     OPEN MASS MATRIX
 
 60 CONTINUE
 j   = lcore - sysbuf
 IF (j > 0) GO TO 65
 icrrqd = IABS(j) + 1
 CALL mesage (-8, icrrqd, NAME)
 65 CALL gopen (mgg,z(j),0)
 it1 = 1
 
!     TEST FOR COUPLED MASS
 
 iy(1) = mgg
 CALL rdtrl (iy)
 cupmas = .false.
 IF (iy(6) == 1) GO TO 90
 IF (iy(6) > 6) cupmas = .true.
 IF (cupmas) GO TO 90
 incr = 0
 ncol = iy(2)
 DO  i = 1,ncol
   ii = 0
   CALL unpack (*70,mgg,a)
   IF (jj-ii > 6) cupmas = .true.
   IF (cupmas) EXIT
 END DO
 80 CALL REWIND (mgg)
 CALL skprec (mgg,1)
 90 ii = 1
 incr = 6
 
!     TEST FOR CONICAL SHELL PROBLEM
 
 nonshl = .true.
 IF (mn == 0) GO TO 100
 nonshl = .false.
 nharms = mn
 nrings = isystm(161)
 iy(1)  = bgpdt
 CALL rdtrl (iy)
 strtmn = iy(2) - nharms*nrings
 iptax  = 0
 kountm = 0
 
!     BRING IN BGPDT
 
 100 FILE = bgpdt
 CALL READ (*410,*330,bgpdt,ri(1),4,0,iflag)
 
!     TEST FOR CONICAL SHELL PROCESSING
 
 IF (nonshl) GO TO 120
 iptax = iptax + 1
 IF (iptax < strtmn) GO TO 110
 kountm = kountm + 1
 IF (kountm <= nrings) GO TO 240
 GO TO 330
 
 110 IF (ir /= -1) CALL skprec (mgg,6)
 
!     CHECK FOR SCALAR POINT
 
 120 CONTINUE
 IF (ir /= -1) GO TO 130
 CALL skprec (mgg,1)
 ii = ii + 1
 GO TO 100
 
!     TEST FOR COUPLED MASS PROCESSING
 
 130 IF (cupmas) GO TO 250
 
!     CONVERT WB'S TO GLOBAL COORDINATES AT RI
 
 DO  i = 1, 3
   wg(i) = wb(i)
 END DO
 IF (ir == 0) GO TO 150
 CALL basglb (wb(1),wg(1),ri(2),ir)
 
!     BRING IN  6X6  ON DIAGONAL OF MASS MATRIX
 
 150 jj = ii + 5
 DO  j = 1,6
   DO  i = 1,6
     xm(i,j) = 0.0
   END DO
 END DO
 DO  i = 1,6
   CALL unpack (*170,mgg,xm(i,1))
 END DO
 
!     MOVE  6X6 TO PARTITIONS
 
 DO  i = 1,3
   DO  j = 1,3
     mt(j,i) = xm(j,i)
     mr(j,i) = xm(j+3,i+3)
     mtr(j,i)= xm(j+3,i)
   END DO
 END DO
 
!     COMPUTE WBX(RI-RA)
 
 DO  i = 1,3
   xm(i,1) = ri(i+1) - ra(i+1)
 END DO
 CALL cross (wb(1),xm(1,1),xm(1,3))
 DO  i = 1,3
   xm(i,1) = xm(i,3)
 END DO
 IF (ir == 0) GO TO 210
 CALL mpyl (ti(1,1),xm(1,1),3,3,1,xm(1,3))
 210 CONTINUE
 
!     COMPUTE MOMENTS
 
 CALL mpyl  (mr(1,1),wg(1),3,3,1,xm(1,1))
 CALL cross (xm(1,1),wg(1),xm(1,2))
 CALL mpylt (mtr(1,1),xm(1,3),3,3,1,xm(1,1))
 CALL cross (xm(1,1),wg,xm(1,4))
 j = ii + 2
 DO  i = 1,3
   j = j + 1
   z(j) = z(j) + xm(i,2) + xm(i,4)
 END DO
 
!     COMPUTE FORCES
 
 CALL mpyl  (mtr(1,1),wg(1),3,3,1,xm(1,1))
 CALL cross (xm(1,1),wg(1),xm(1,2))
 CALL mpyl  (mt(1,1),xm(1,3),3,3,1,xm(1,1))
 CALL cross (xm(1,1),wg,xm(1,4))
 j = ii - 1
 DO  i = 1,3
   j = j + 1
   z(j) = z(j) + xm(i,4) + xm(i,2)
 END DO
 
!     BUMP  II
 
 ii  = ii + 6
 GO TO 100
 
!     CONICAL SHELL PROCESSING
!     COMPUTE A = R*WB**2
 
 240 xm(2,3) = 0.0
 xm(3,3) = 0.0
 xm(1,3) = ri(2)*wb(2)*wb(2)
 GO TO 290
 
!     COUPLED MASS PROCESSING
!     COMPUTE -WB*(WB*(RI - RA))
 
 250 DO  i = 1, 3
   xm(i,1) = ri(i+1) - ra(i+1)
 END DO
 CALL cross (wb(1),xm(1,1),xm(1,3))
 CALL cross (xm(1,3),wb(1),xm(1,1))
 IF (ir == 0) GO TO 270
 CALL basglb (xm(1,1),xm(1,3),ri(2),ir)
 GO TO 290
 270 DO  i = 1, 3
   xm(i,3) = xm(i,1)
 END DO
 
!     COMPUTE F = M*A
 
 290 i1 = 1
 DO  i = 1, 3
   CALL intpk (*320,mgg,0,i1,0)
   IF (xm(i,3) == 0.0) GO TO 310
   300 CALL zntpki
   z(irow) = z(irow) + a(1)*xm(i,3)
   IF (ieol /= 1) GO TO 300
   CYCLE
   310 CALL skprec (mgg,1)
 END DO
 CALL skprec (mgg,3)
 GO TO 100
 
!     EOR IN BGPDT
 
 330 CALL CLOSE  (mgg,1)
 CALL REWIND (bgpdt)
 OLD = 0
 CALL skprec (bgpdt,1)
 RETURN
 
!     FILE ERRORS
 
 400 CALL mesage (ip1,FILE,NAME(1))
 410 ip1 = -2
 GO TO 400
END SUBROUTINE rforce
