SUBROUTINE mred2d
     
!     THIS SUBROUTINE CALCULATES THE MODAL MASS AND STIFFNESS MATRICES
!     IF USERMODE = TYPE2 FOR THE MRED2 MODULE.
 
!     INPUT DATA
!     GINO - USETMR  - USET TABLE FOR REDUCED SUBSTRUCTURE
!            LAMAMR  - EIGENVALUE TABLE FOR SUBSTRUCTURE BEING REDUCED
!            PHISS   - EIGENVECTORS FOR SUBSTRUCTURE BEING REDUCED
!            QSM     - MODEL REACTION MATRIX
!            PAA     - SUBSTRUCTURE LOAD MATRIX
 
!     OUTPUT DATA
!     GINO - KHH     - REDUCED STIFFNESS MATRIX
!            MHH     - REDUCED MASS MATRIX
!            PHH     - REDUCED LOAD MATRIX
!     SOF  - HORG    - H TRANSFORMATION MATRIX
!            KMTX    - STIFFNESS MATRIX FOR REDUCED SUBSTRUCTURE
!            MMTX    - MASS MATRIX FOR REDUCED SUBSTRUCTURE
!            PVEC    - LOAD MATRIX FOR REDUCED SUBSTRUCTURE
!            PAPP    - APPENDED LOAD MATRIX FOR REDUCED SUBSTRUCTURE
!            POVE    - INTERNAL POINT LOADS FOR ORIGINAL SUBSTRUCTURE
!            POAP    - INTERNAL POINTS APPENDED LOADS FOR ORIGINAL
!                       SUBSTRUCTURE
 
!     PARAMETERS
!     INPUT - DRY    - MODULE OPERATION FALG
!             GBUF   - GINO BUFFERS
!             INFILE - INPUT FILE NUMBERS
!             OTFILE - OUTPUT FILE NUMBERS
!             ISCR   - SCRATCH FILE NUMBERS
!             KORLEN - LENGTH OF OPEN CORE
!             KORBGN - BEGINNING ADDRESS OF OPEN CORE
!             OLDNAM - NAME OF SUBSTRUCTURE BEING REDUCED
!             NEWNAM - NAME OF REDUCED SUBSTRUCTURE
!             USRMOD - USERMODE FLAG
 
 INTEGER :: dry,popt,gbuf1,sbuf1,sbuf2,sbuf3,otfile,oldnam,  &
     usrmod,z,ul,ua,uf,us,un,ub,typin,typout,typea,  &
     typeb,paa,phh,rprtn,cprtn,bbzero,zero,usetmr,snb, papp
 DIMENSION       itrlr(7),itrlr1(7),itrlr2(7),modnam(2),itmlst(6),  &
     BLOCK(11),isub(4),itmnam(2),rz(1)
 CHARACTER (LEN=23) :: ufm
 COMMON /xmssg / ufm
 COMMON /BLANK / idum1,dry,popt,gbuf1,idum3(2),sbuf1,sbuf2,sbuf3,  &
     infile(12),otfile(6),iscr(10),korlen,korbgn,  &
     oldnam(2),newnam(2),idum4(4),usrmod,idum2(5), nmodes
 COMMON /zzzzzz/ z(1)
 COMMON /two   / itwo(32)
 COMMON /bitpos/ idum5(5),ul,ua,uf,us,un,idum6(10),ub
 COMMON /packx / typin,typout,irow,nrow,incr
 COMMON /system/ idum7,iprntr
 EQUIVALENCE     (usetmr,infile(5)),(kbb,infile(6)),  &
     (mbb,infile(7)),(paa,infile(10)),  &
     (khh,otfile(1)),(mhh,otfile(2)),(phh,otfile(5)),  &
     (rprtn,iscr(8)),(cprtn,iscr(8)),(k,iscr(3)),  &
     (bbzero,iscr(9)),(m,iscr(10)),(zero,iscr(3)),  &
     (rz(1),z(1)),(typea,BLOCK(1)),(typeb,BLOCK(7))
 DATA    modnam/ 4HMRED,4H2D  /
 DATA    papp  / 4HPAPP/
 DATA    mred2 / 27    /
 DATA    itmlst/ 4HKMTX,4HMMTX,4HPVEC,4HPAPP,4HPOVE,4HPOAP/
 
!     CHECK USERMODE OPTION FLAG
 
 IF (dry == -2) GO TO 400
 
!     COUNT NUMBER OF FREE, FIXED POINTS WITHIN BOUNDARY SET
 
 itrlr(1) = usetmr
 CALL rdtrl(itrlr)
 ifile = usetmr
 IF (itrlr(1) < 0) GO TO 270
 luset = itrlr(3)
 IF ((korbgn + luset) >= korlen) GO TO 280
 CALL gopen (usetmr,z(gbuf1),0)
 CALL READ  (*260,*270,usetmr,z(korbgn),luset,0,nwdsrd)
 CALL CLOSE (usetmr,1)
 nuf   = 0
 nus   = 0
 snb   = itwo(us) + itwo(un) + itwo(ub)
 lafnb = itwo(ul) + itwo(ua) + itwo(uf) + itwo(un) + itwo(ub)
 DO  i = 1,luset
   IF (z(korbgn+i-1) == lafnb) nuf = nuf + 1
   IF (z(korbgn+i-1) ==   snb) nus = nus + 1
 END DO
 
!     IF FIXED SET, COMPUTE GS MATRIX
 
 IF (nus == 0) GO TO 20
 CALL mred2i (1,0,0)
 
!     IF FREE SET, PARTITION PHISS
 
 20 IF (nuf == 0) GO TO 50
 CALL mred2j (nuf,n2)
 
!     FORM HK MATRIX
 
 CALL mred2l (nuf,n2,nus,ufbits)
 50 CALL mred2m (nuf,n2,nus)
 
!     COMPUTE K MATRIX
 
 CALL mred2n
 
!     COMPUTE HM MATRIX
 
 CALL mred2o (nus)
 
!     OUTPUT HORG
 
 CALL mred2p (nus,nuf,n2)
 
!     PROCESS STIFFNESS, MASS MATRICES
!     II = 1, PROCESS STIFFNESS MATRIX
!     II = 2, PROCESS MASS MATRIX
 
 IF (dry == -2) GO TO 240
 CALL setlvl (newnam,1,oldnam,itest,mred2)
 IF (itest == 8) GO TO 380
 DO  ii = 1,2
   itrlr1(1) = kbb
   km   = k
   kmhh = khh
   IF (ii == 1) GO TO 60
   itrlr1(1) = mbb
   km   = m
   kmhh = mhh
   60 kmbb = itrlr1(1)
   CALL rdtrl (itrlr1)
   IF (itrlr1(1) < 0) GO TO 160
   CALL sofcls
   
!     FORM MERGE VECTOR
   
   jrow   = itrlr1(3)
   kolumn = itrlr1(2)
   itrlr2(1) = km
   CALL rdtrl (itrlr2)
   nrow   = itrlr2(3)
   kolmns = itrlr2(2)
   DO  i = 1,nrow
     rz(korbgn+i-1) = 0.0
     IF (i > jrow) rz(korbgn+i-1) = 1.0
   END DO
   iform  = 7
   typin  = 1
   typout = 1
   irow   = 1
   incr   = 1
   CALL makmcb (itrlr1,rprtn,nrow,iform,typin)
   CALL gopen (rprtn,z(gbuf1),1)
   CALL pack (z(korbgn),rprtn,itrlr1)
   CALL CLOSE (rprtn,1)
   CALL wrttrl (itrlr1)
   
!     MERGE (K,M)BB MATRIX WITH ZERO MATRICES
   
   isub(1) = kolumn
   isub(2) = kolmns - kolumn
   isub(3) = jrow
   isub(4) = nrow - jrow
   itype   = 1
   CALL gmmerg (bbzero,kmbb,0,0,0,rprtn,rprtn,isub,itype,z(korbgn), korlen)
   
!     FORM STIFFNESS, MASS MATRICES
   
!                                  **           **
!                                  *         .   *
!        **       **   **     **   * (K,M)BB . 0 *
!        *         *   *       *   *         .   *
!        * (K,M)HH * = * (K,M) * + *.............*
!        *         *   *       *   *         .   *
!        **       **   **     **   *    0    . 0 *
!                                  *         .   *
!                                  **           **
   DO  i = 1,11
     BLOCK(i) = 0.0
   END DO
   BLOCK(2) = 1.0
   BLOCK(8) = 1.0
   typea = itrlr2(5)
   typeb = itrlr1(5)
   iop   = 1
   CALL ssg2c (km,bbzero,kmhh,iop,BLOCK)
   CALL sofopn (z(sbuf1),z(sbuf2),z(sbuf3))
   GO TO 170
   
!     NO BB MATRIX PARTITION
   
   160 kmhh = km
   
!     STORE MATRIX ON SOF
!     II = 1, STORE KHH AS KMTX
!     II = 2, STORE MHH AS MMTX
   
   170 item      = itmlst(ii)
   itmnam(1) = newnam(1)
   itmnam(2) = newnam(2)
   CALL mtrxo (kmhh,newnam,item,0,itest)
   IF (itest /= 3) GO TO 300
 END DO
 
!     PROCESS LOAD DATA
 
 itrlr1(1) = paa
 CALL rdtrl (itrlr1)
 IF (itrlr1(1) < 0) GO TO 240
 
!     EXPAND PAA FOR MODAL DOF
 
!                  **   **
!                  *     *
!        **   **   * PAA *
!        *     *   *     *
!        * PHH * = *.....*
!        *     *   *     *
!        **   **   *  0  *
!                  *     *
!                  **   **
 
 nrow = itrlr1(3) + n2
 IF (n2 == 0) nrow = nrow + (nmodes - nuf)
 iform  = 7
 typin  = 1
 typout = 1
 irow   = 1
 incr   = 1
 CALL makmcb (itrlr2,cprtn,nrow,iform,typin)
 DO  i = 1,nrow
   rz(korbgn+i-1) = 0.0
   IF (i > itrlr1(3)) rz(korbgn+i-1) = 1.0
 END DO
 CALL gopen (cprtn,z(gbuf1),1)
 CALL pack (z(korbgn),cprtn,itrlr2)
 CALL CLOSE (cprtn,1)
 CALL wrttrl (itrlr2)
 
!     MERGE PAA WITH ZERO MATRIX
 
 isub(3) = itrlr1(3)
 isub(4) = n2
 IF (n2 == 0) isub(4) = nmodes - nuf
 itype   = 1
 CALL gmmerg (phh,paa,0,0,0,0,cprtn,isub,itype,z(korbgn),korlen)
 
!     SAVE PHH AS PVEC OR PAPP ON SOF
 
 item = itmlst(3)
 IF (popt == papp) item = itmlst(4)
 itmnam(1) = newnam(1)
 itmnam(2) = newnam(2)
 CALL mtrxo (phh,newnam,item,0,itest)
 IF (itest /= 3) GO TO 300
 
!     STORE NULL MATRIX AS POVE OR POAP ON SOF
 
 iform  = 2
 kolmns = itrlr1(2)
 nrow   = n2
 IF (n2 == 0) nrow = nmodes - nuf
 CALL makmcb (itrlr2,zero,nrow,iform,typin)
 CALL gopen (zero,z(gbuf1),1)
 DO  i = 1,kolmns
   DO  j = 1,nrow
     rz(korbgn+j-1) = 0.0
   END DO
   CALL pack (z(korbgn),zero,itrlr2)
 END DO
 CALL CLOSE (zero,1)
 CALL wrttrl (itrlr2)
 item = itmlst(5)
 IF (popt == papp) item = itmlst(6)
 itmnam(1) = oldnam(1)
 itmnam(2) = oldnam(2)
 CALL mtrxo (zero1,oldnam,item,0,itest)
 IF (itest /= 3) GO TO 300
 240 CONTINUE
 GO TO 400
 
!     PROCESS SYSTEM FATAL ERRORS
 
 260 imsg = -2
 GO TO 290
 270 imsg = -3
 GO TO 290
 280 imsg = -8
 ifile = 0
 290 CALL sofcls
 CALL mesage (imsg,ifile,modnam)
 GO TO 400
 
!     PROCESS MODULE FATAL ERRORS
 
 300 SELECT CASE ( itest )
   CASE (    1)
     GO TO 310
   CASE (    2)
     GO TO 320
   CASE (    3)
     GO TO 330
   CASE (    4)
     GO TO 340
   CASE (    5)
     GO TO 350
   CASE (    6)
     GO TO 370
 END SELECT
 310 imsg = -9
 GO TO 390
 320 imsg = -11
 GO TO 390
 330 imsg = -1
 GO TO 360
 340 imsg = -2
 GO TO 360
 350 imsg = -3
 360 CALL smsg (imsg,item,itmnam)
 GO TO 400
 370 imsg = -10
 GO TO 390
 380 WRITE  (iprntr,385) ufm
 385 FORMAT (a23,' 6518, ONE OF THE COMPONENT SUBSTRUCTURES HAS BEEN ',  &
     'USED IN A PREVIOUS COMBINE OR REDUCE.')
 dry = -2
 GO TO 400
 390 dry = -2
 CALL smsg1 (imsg,item,itmnam,modnam)
 400 RETURN
END SUBROUTINE mred2d
