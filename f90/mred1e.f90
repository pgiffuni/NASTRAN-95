SUBROUTINE mred1e
     
!     THIS SUBROUTINE GENERATES THE RIGID BODY MATRIX DMX IF FREEBODY
!     MODES ARE REQUESTED FOR THE MRED1 MODULE.
 
!     INPUT DATA
!     SOF -    BGSS   - BASIC GRID POINT IDENTIFICATION TABLE
!              EQSS   - SUBSTRUCTURE EQUIVALENCE TABLE
!              CSTM   - COORDINATE SYSTEM TRANSFORMATION MATRIX
 
!     OUTPUT DATA
!     GINO -   SCR1   - SCRATCH FILE HOLDING UNTRANSPOSED DMX MATRIX
!              DMR    - RIGID BODY MATRIX
 
!     PARAMETERS
!     INPUT  - DRY    - MODULE OPERATION FLAG
!              RGRID  - FREEBODY MODES FLAGS
!                       RGRID(1) .EQ. INTERNAL GRID POINT IDENTIFICATION
!                                     NUMBER (SET IN MRED1C)
!                       RGRID(2) .EQ. NUMBER OF THE CONTRIBUTING
!                                     SUBSTRUCTURE (SET IN MRED1)
!              KORBGN - BEGINNING ADDRESS OF OPEN CORE
!              KORLEN - LENGTH OF OPEN CORE
!              RGRID0 - FREE BODY MODE BASIC COORDINATES
!     OTHERS - NBGSS  - NUMBER OF INTERNAL GRID IDENTIFICATION POINTS
!              LOCBGS - BEGINNING ADDRESS OF BGSS DATA
!              LOCSTM - BEGINNING ADDRESS OF CSTM DATA
!              LOCSIL - BEGINNING ADDRESS OF SIL DATA
!              SMALD  - MATRIX OF COORDINATE LOCATION DIFFERENCES (3X3)
!              TI     - MATRIX OF COORDINATE TRANSFORMATIONS (3X3)
!              BIGD   - PARTITIONED MATRIX OF TRANSFORMATIONS (6X6)
 
!                                                  T
!              TTD    - TEMPORARY MATRIX HOLDING (T SMALD) (3X3)
!              KOMPNT - ARRAY HOLDING DECODED SIL COMPONENTS
 
 INTEGER :: oldnam,dry,gbuf1,gbuf2,rgrid,rname,z,typin,typpck,  &
     typunp,scr1,dmr,tittd,tiijd1,ttdijd,zeroij,tiijd2, dmrnam
 DIMENSION       modnam(2),itrlr(7),ti(9),smald(9),bigd(36),ttd(9),  &
     kompnt(32),rz(1),dmrnam(2)
 COMMON /BLANK / oldnam(2),dry,idum1(6),gbuf1,gbuf2,idum2(3),  &
     korlen,idum3(6),rgrid(2),rname(2),idum4,korbgn,  &
     ncsubs,idum5(3),nsil,idum6(4),rgrid0(3)
 COMMON /zzzzzz/ z(1)
 COMMON /packx / typin,typpck,irowp,lrowp,incrp
 COMMON /unpakx/ typunp,irowup,lrowup,incrup
 EQUIVALENCE     (rz(1),z(1))
 DATA    modnam/ 4HMRED,4H1E   /
 DATA    nhbgss, nhcstm,nheqss /4HBGSS,4HCSTM,4HEQSS/
 DATA    scr1  , dmr /301,204  /
 
!     TEST FOR MODULE ERRORS
 
 IF (dry == -2) GO TO 190
 
!     TEST FOR FREEBODY MODES REQUEST
 
 IF (rgrid(1) == -1) GO TO 190
 
!     READ BGSS DATA
 
 it = 1
 CALL sfetch (oldnam,nhbgss,1,itest)
 IF (itest == 3) GO TO 240
 IF (itest == 4) GO TO 250
 IF (itest == 5) GO TO 260
 CALL suread (z(korbgn),-1,nwdsrd,itest)
 
!     EXTRACT SUBSTRUCTURE IP DATA
 
 nbgss  = z(korbgn+2)
 locbgs = korbgn
 CALL suread (z(korbgn),-2,nwdsrd,itest)
 korbgn = korbgn + nwdsrd
 
!     READ CSTM DATA
 
 locstm = korbgn
 IF (korlen <= locstm) GO TO 200
 it = 2
 CALL sfetch (oldnam,nhcstm,1,itest)
 IF (itest == 3) GO TO 30
 IF (itest == 4) GO TO 250
 IF (itest == 5) GO TO 260
 CALL suread (z(locstm),-2,nwdsrd,itest)
 CALL pretrs (z(locstm+3),nwdsrd-4)
 
!     CHECK FOR BASIC COORDINATES
 
 30 DO  i = 1, 3
   rgrid0(i) = 0.0
 END DO
 IF (rgrid(1) == 0) GO TO 60
 
!     EXTRACT FREEBODY BASIC COORDINATES
 
 locrgr = locbgs + (4*(rgrid(1)-1))
 DO  i = 1,3
   rgrid0(i) = rz(locrgr+i)
 END DO
 
!     OPEN SCRATCH FILE
 
 60 ifile = scr1
 itrlr(1) = ifile
 CALL OPEN (*210,scr1,z(gbuf2),1)
 typin = 1
 typpck= 1
 irowp = 1
 lrowp = 6
 incrp = 1
 
!     OPEN EQSS FILE AND CHECK OPEN CORE LENGTH
 
 it = 3
 CALL sfetch (oldnam,nheqss,1,itest)
 IF (itest == 3) GO TO 240
 IF (itest == 4) GO TO 250
 IF (itest == 5) GO TO 260
 locsil = locstm + nwdsrd
 CALL suread (z(locsil),-1,nwdsrd,itest)
 IF (korlen <= locsil) GO TO 240
 
!     READ UP TO SIL DATA
 
 IF (korlen <= 2*nsil) GO TO 240
 DO  i = 1,ncsubs
   CALL suread (z(locsil),-1,nwdsrd,itest)
   IF (korlen <= locsil+nwdsrd) GO TO 240
 END DO
 
!     GENERATE SMALD MATRIX (3X3)
 
!                **                               **
!                *                                 *
!                *    0.0      DELTA(Z)  -DELTA(Y) *
!                *                                 *
!        SMALD = * -DELTA(Z)     0.0      DELTA(X) *
!                *                                 *
!                *  DELTA(Y)  -DELTA(X)     0.0    *
!                *                                 *
!                **                               **
 
 DO  i = 1,nbgss
   ii       = 4*(i-1)
   smald(1) = 0.0
   smald(2) = rz(locbgs+ii+3) - rgrid0(3)
   smald(3) =-rz(locbgs+ii+2) + rgrid0(2)
   smald(4) =-smald(2)
   smald(5) = 0.0
   smald(6) = rz(locbgs+ii+1) - rgrid0(1)
   smald(7) =-smald(3)
   smald(8) =-smald(6)
   smald(9) = 0.0
   
!     SELECT TI, TTD MATRIX GENERATION
   
   IF (z(locbgs+ii) < 0.0) THEN
     GO TO   120
   ELSE IF (z(locbgs+ii) == 0.0) THEN
     GO TO    85
   END IF
   
!     GENERATE TI, TTD MATRICES (3X3)
!     (CID .GT. 0)
   
   80   CALL transs (z(locbgs+ii),ti)
   CALL gmmats (ti,3,3,0,smald,3,3,1,ttd)
   GO TO 95
   
!     GENERATE TI, TTD MATRICES (3X3)
!     (CID .EQ. 0)
   
   85 DO  j = 1,3
     DO  k = 1,3
       l = k + 3*(j-1)
       ti(l) = 0.0
       IF (j == k) ti(l) = 1.0
       ttd(l) = smald(l)
     END DO
   END DO
   
!     GENERATE BIGD MATRIX (6X6)
   
!               **            **
!               *    .         *
!               *  T .  T      *
!               * T  . T SMALD *
!               *    .         *
!        BIGD = *..............*
!               *    .         *
!               *    .    T    *
!               * 0  .  T      *
!               *    .         *
!               **            **
   
   95 DO  j = 1,3
     DO  k = 1,3
       tittd  = k + 3*(j-1)
       tiijd1 = k + 6*(j-1)
       ttdijd = tiijd1 + 3
       zeroij = tiijd1 + 18
       tiijd2 = tiijd1 + 21
       bigd(tiijd1) = ti(tittd)
       bigd(ttdijd) = ttd(tittd)
       bigd(zeroij) = 0.0
       bigd(tiijd2) = ti(tittd)
     END DO
   END DO
   
!     EXTRACT ROWS OF BIGD CORRESPONDING TO ACTIVE SIL COMPONENTS
   
   CALL suread (z(locsil),2,nwdsrd,itest)
   icode = z(locsil+1)
   CALL decode (icode,kompnt,nwdsd)
   DO  j = 1,nwdsd
     irowd = 1 + 6*kompnt(j)
     CALL pack (bigd(irowd),scr1,itrlr)
   END DO
   CYCLE
   
!     SCALAR POINT ADDS NULL COLUMN TO BIGD
!     (CID .LT. 0)
   
   120 DO  j = 1,6
     bigd(j) = 0.0
   END DO
   irowp = 1
   CALL pack (bigd(1),scr1,itrlr)
 END DO
 CALL CLOSE (scr1,1)
 itrlr(3) = lrowp
 
!     READ SCR1 INTO TRANSPOSED FORM
 
 CALL OPEN (*210,scr1,z(gbuf1),0)
 typunp = 1
 irowup = 1
 lrowup = 6
 incrup = itrlr(2)
 kolmns = itrlr(2)
 korbgn = locbgs
 IF (korlen <= korbgn+lrowp*kolmns) GO TO 240
 DO  i = 1,kolmns
   CALL unpack (*150,scr1,z(korbgn))
   GO TO 170
   150 j = korbgn
   DO  k = 1,6
     rz(j) = 0.0
     j = j + incrup
   END DO
   170 korbgn = korbgn + 1
 END DO
 CALL CLOSE (scr1,1)
 
!     PLACE TRANSPOSED BIGD ONTO DMR OUTPUT FILE
 
 ifile  = dmr
 CALL OPEN (*210,dmr,z(gbuf2),1)
 CALL fname (dmr,dmrnam)
 CALL WRITE (dmr,dmrnam,2,1)
 locdmr = locbgs
 lrowp  = kolmns
 iform  = 2
 CALL makmcb (itrlr,dmr,lrowp,iform,typin)
 DO  i = 1,6
   CALL pack (z(locdmr),dmr,itrlr)
   locdmr = locdmr + kolmns
 END DO
 CALL CLOSE (dmr,1)
 CALL wrttrl (itrlr)
 190 RETURN
 
!     PROCESS SYSTEM FATAL ERRORS
 
 200 imsg  =-8
 ifile = 0
 GO TO 230
 210 imsg = -1
 230 CALL sofcls
 CALL mesage (imsg,ifile,modnam)
 GO TO 190
 
!     PROCESS MODULE FATAL ERRORS
 
 240 imsg = -1
 GO TO 270
 250 imsg = -2
 GO TO 270
 260 imsg = -3
 270 IF (it-2 < 0) THEN
   GO TO   280
 ELSE IF (it-2 == 0) THEN
   GO TO   290
 ELSE
   GO TO   300
 END IF
 280 CALL smsg (imsg,nhbgss,oldnam)
 GO TO 190
 
 290 CALL smsg (imsg,nhcstm,oldnam)
 GO TO 190
 
 300 CALL smsg (imsg,nheqss,oldnam)
 GO TO 190
 
END SUBROUTINE mred1e
