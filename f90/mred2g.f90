SUBROUTINE mred2g (kode)
     
!     THIS SUBROUTINE CALCULATES THE FINAL STRUCTURAL MATRICES FOR THE
!     MRED2 MODULE.
 
!     INPUT DATA
!     GINO -   KBB    - STIFFNESS PARTITION MATRIX
!              KIB    - KIB  STIFFNESS PATTITION MATRIX
!              HIE    - HIE  PARTITION MATRIX
!              KII    - KII  PARTITION MATRIX
!              HGH    - HORG PARTITION MATRIX
!              MAA    - MASS INPUT MATRIX
!              BAA    - DAMPING INPUT MATRIX
!              K4AA   - STIFFNESS INPUT MATRIX
!              PAA    - LOADS INPUT MATRIX
!     SOF  -   GIMS   - G TRANSFORMATION MATRIX
 
!     OUTPUT DATA
!     GINO -   KHH    - STIFFNESS MATRIX
!              MHH    - MASS MATRIX
!              BHH    - DAMPING MATRIX
!              K4HH   - K4HH  MATRIX
!              PHH    - LOADS MATRIX
!     SOF  -   KMTX   - STIFFNESS MATRIX
!              MMTX   - MASS  MATRIX
!              PVEC   - LOADS MATRIX
!              PAPP   - APPENDED LOADS MATRIX
!              BMTX   - DAMPING MATRIX
!              K4MX   - K4MX STIFFNESS MATRIX
 
!     PARAMETERS
!     INPUT  - POPT   - LOADS OPTION FLAG
!              GBUF   - GINO BUFFERS
!              INFILE - INPUT   FILE NUMBERS
!              OTFILE - OUTPUT  FILE NUMBERS
!              ISCR   - SCRATCH FILE NUMBERS
!              KORLEN - LENGTH OF OPEN CORE
!              KORBGN - BEGINNING ADDRESS OF OPEN CORE
!              OLDNAM - NAME OF SUBSTRUCTURE BEING REDUCED
!     OTHERS - PAA    - LOADS INPUT FILE NUMBER
!              KHH    - STIFFNESS I
!              KHH    - STIFFNESS OUTPUT FILE NUMBER
!              POVE   - LOADS OUTPUT FILE NUMBER
!              UPRT   - PARTITION VECTOR FILE NUMBER
!              ZEROEB - ZERO PARTITION FILE NUMBER
!              KBB    - KBB INPUT FILE NUMBER
!              ZEROBE - ZERO PARTITION MATRIX
!              KIB    - KIB INPUT FILE NUMBER
!              KII    - KII INPUT FILE NUMBER
!              KBARBB - KBARBB FILE NU BER
!              GIB    - GIB INPUT FILE NUMBER
!              KEE    - KEE FILE NUMBER
!              HGH    - HORG INPUT FILE NUMBER
 
 
 INTEGER, INTENT(IN OUT)                  :: kode
 LOGICAL :: frebdy,bounds,modes,ponly
 INTEGER :: dry,popt,gbuf1,sbuf1,sbuf2,sbuf3,otfile,oldnam,z,  &
     t,signab,signc,prec,scr,typin,typout,un,ub,ui,  &
     fuset,prec3,zerobe,zeroeb,blanks,papp,paa,pove,  &
     uprt,gib,hgh,hie,rprtn,cprtn,usetmr,dblkor,eqst
 DOUBLE PRECISION :: dz
 DIMENSION        modnam(2),itrlr1(7),itrlr2(7),itrlr3(7),isub(4),  &
     itmlst(11),itmnam(2),rz(1),dz(1)
 CHARACTER (LEN=23) :: ufm
 COMMON /xmssg /  ufm
 COMMON /BLANK /  idum1,dry,popt,gbuf1,idum2(2),sbuf1,sbuf2,sbuf3,  &
     infile(12),otfile(6),iscr(10),korlen,korbgn,  &
     oldnam(2),newnam(2),frebdy,idum6(5),bounds,modes, idum7(4),ponly,lstzwd
 COMMON /zzzzzz/  z(1)
 COMMON /system/  idum3,iprntr
 COMMON /mpyadx/  itrlra(7),itrlrb(7),itrlrc(7),itrlrd(7),nz,t,  &
     signab,signc,prec,scr,dumm
 COMMON /packx /  typin,typout,irow,nrow,incr
 COMMON /bitpos/  idum4(9),un,idum5(10),ub,ui
 COMMON /patx  /  lcore,nsub(3),fuset
 COMMON /mpy3tl/  jtrlra(7),jtrlrb(7),jtrlre(7),jtrlrc(7),jscr(3),  &
     lkore,icode,prec3,dummy(13)
 EQUIVALENCE      (eqst,infile(5)),(usetmr,infile(5)),  &
     (paa,infile(10)),(khh,otfile(1)),(pove,otfile(6))
 EQUIVALENCE      (zerobe,iscr(1)),(uprt,iscr(1)),(kib,iscr(2)),  &
     (zeroeb,iscr(3)),(kii,iscr(3)),(kbb,iscr(1)),  &
     (gib,iscr(4)),(kbarbb,iscr(5)),(kee,iscr(6)),  &
     (hie,iscr(7)),(hgh,iscr(8)),(rprtn,iscr(2)),  &
     (cprtn,iscr(4)),(rz(1),z(1)),(dz(1),z(1))
 DATA    modnam/  4HMRED,4H2G  /,      papp  ,blanks/4HPAPP,4H    /
 DATA    itmlst/  4HKMTX,4HMMTX,4HBMTX,4HK4MX,4HPVEC,4HPAPP,4HPOVE,  &
     4HGIMS,4HHORG,4HPOAP,4HUPRT/
 DATA    mred2 /  27 /
 
!     SELECT OPERATION
!     KODE = 1, NO SETLVL, NO STIFFNESS CALCULATIONS
!     KODE = 2, SETLVL, STIFFNESS CALCULATIONS
!     KODE = 3, NO SETLVL, NO STIFFNESS CALCULATIONS
!     KODE = 4, SETLVL, NO STIFFNESS CALCULATIONS
 
 IF (dry == -2) GO TO 300
 SELECT CASE ( kode )
   CASE (    1)
     GO TO 90
   CASE (    2)
     GO TO 1
   CASE (    3)
     GO TO 90
   CASE (    4)
     GO TO 1
 END SELECT
 
!     SET UP NEW SUBSTRUCTURE
 
 1 IF (bounds .OR. modes) GO TO 5
 numb = 1
 CALL setlvl (newnam,numb,oldnam,itest,mred2)
 IF (itest == 8) GO TO 290
 5 IF (kode  == 4) GO TO 90
 
!     FORM PRELIMINARY STIFFNESS CALCULATION
 
!                                      T
!        **      **   **   **   **   ** **   **
!        *        *   *     *   *     * *     *
!        * KBARBB * = * KBB * + * GIB * * KIB *
!        *        *   *     *   *     * *     *
!        **      **   **   **   **   ** **   **
 
 itrlr1(1) = kbb
 CALL rdtrl (itrlr1)
 item = itmlst(8)
 itmnam(1) = oldnam(1)
 itmnam(2) = oldnam(2)
 CALL softrl (oldnam,item,itrlr2)
 itest = itrlr2(1)
 IF (itest /= 1) GO TO 200
 CALL mtrxi (gib,oldnam,item,0,itest)
 IF (itest /= 1) GO TO 200
 CALL sofcls
 itrlr2(1) = gib
 CALL rdtrl (itrlr2)
 itrlr3(1) = kib
 CALL rdtrl (itrlr3)
 DO  i = 1,7
   itrlra(i) = itrlr2(i)
   itrlrb(i) = itrlr3(i)
   itrlrc(i) = itrlr1(i)
 END DO
 iform = 6
 iprc  = 1
 ityp  = 0
 IF ((itrlra(5) == 2) .OR. (itrlra(5) == 4)) iprc = 2
 IF ((itrlrb(5) == 2) .OR. (itrlrb(5) == 4)) iprc = 2
 IF ((itrlrc(5) == 2) .OR. (itrlrc(5) == 4)) iprc = 2
 IF (itrlra(5) >= 3) ityp = 2
 IF (itrlrb(5) >= 3) ityp = 2
 IF (itrlrc(5) >= 3) ityp = 2
 itype = iprc + ityp
 CALL makmcb (itrlrd,kbarbb,itrlr1(3),iform,itype)
 t      = 1
 signab = 1
 signc  = 1
 prec   = 0
 scr    = iscr(9)
 dblkor = korbgn/2 + 1
 nz = lstzwd - (2*dblkor - 1)
 CALL mpyad  (dz(dblkor),dz(dblkor),dz(dblkor))
 CALL wrttrl (itrlrd)
 kbarow = itrlrd(3)
 kcol   = itrlrd(2)
 
!     FORM PRELIMINARY STIFFNESS CALCULATION
 
!                         T
!        **   **   **   ** **   ** **   **
!        *     *   *     * *     * *     *
!        * KEE * = * HIE * * KII * * HIE *
!        *     *   *     * *     * *     *
!        **   **   **   ** **   ** **   **
 
 itrlr1(1) = hie
 itrlr2(1) = kii
 CALL rdtrl (itrlr1)
 CALL rdtrl (itrlr2)
 DO  i = 1,7
   jtrlra(i) = itrlr1(i)
   jtrlrb(i) = itrlr2(i)
   jtrlre(i) = 0
 END DO
 iprc = 1
 ityp = 0
 IF (jtrlra(5) == 2 .OR. jtrlra(5) == 4) iprc = 2
 IF (jtrlrb(5) == 2 .OR. jtrlrb(5) == 4) iprc = 2
 IF (jtrlra(5) >= 3 .OR. jtrlrb(5) >= 3) ityp = 2
 itype = iprc + ityp
 CALL makmcb (jtrlrc,kee,itrlr1(2),iform,itype)
 jscr(1) = iscr(9)
 jscr(2) = iscr(2)
 jscr(3) = iscr(1)
 lkore   = nz
 icode   = 0
 prec3   = 0
 CALL mpy3dr (dz(dblkor))
 CALL wrttrl (jtrlrc)
 keerow = jtrlrc(3)
 keecol = jtrlrc(2)
 
!     GENERATE MERGE PARTITION VECTOR
 
 nrow = kcol + keecol
 DO  i = 1,nrow
   rz(korbgn+i-1) = 0.0
   IF (i > kcol) rz(korbgn+i-1) = 1.0
 END DO
 typin  = 1
 typout = 1
 irow   = 1
 incr   = 1
 iform  = 7
 CALL makmcb (itrlr1,rprtn,nrow,iform,typin)
 CALL gopen  (rprtn,z(gbuf1),1)
 CALL pack   (rz(korbgn),rprtn,itrlr1)
 CALL CLOSE  (rprtn,1)
 CALL wrttrl (itrlr1)
 
!     FORM STIFFNESS MATRIX
 
!                  **            **
!                  *        .     *
!        **   **   * KBARBB .  0  *
!        *     *   *        .     *
!        * KHH * = *..............*
!        *     *   *        .     *
!        **   **   *   0    . KEE *
!                  *        .     *
!                  **            **
 
 isub(1) = kcol
 isub(2) = keecol
 isub(3) = kbarow
 isub(4) = keerow
 iform   = 6
 CALL gmmerg (khh,kbarbb,0,0,kee,rprtn,rprtn,isub,iform,z(korbgn), korlen)
 
!     STORE KHH AS KMTX ON SOF
 
 CALL sofopn (z(sbuf1),z(sbuf2),z(sbuf3))
 itmnam(1) = newnam(1)
 itmnam(2) = newnam(2)
 CALL mtrxo (khh,newnam,itmlst(1),0,itest)
 item = itmlst(1)
 IF (itest /= 3) GO TO 200
 GO TO 100
 
!     LOCATE HGH MATRIX
 
 90 CALL mtrxi (hgh,oldnam,itmlst(9),0,itest)
 item = itmlst(9)
 itmnam(1) = oldnam(1)
 itmnam(2) = oldnam(2)
 IF (itest /= 1) GO TO 200
 100 signab = 1
 signc  = 1
 scr    = iscr(1)
 dblkor = korbgn/2 + 1
 lcore  = lstzwd - (2*dblkor - 1)
 
!     GENERATE MATRICES REQUESTED
!     I = 2, GENERATE MHH MATRIX
!     I = 3, GENERATE BHH MATRIX
!     I = 4, GENERATE K4HH MATRIX
!     I = 5, GENERATE PHH MATRIX
 
 DO  i = 2,5
   itrlr1(1) = infile(i+5)
   CALL rdtrl (itrlr1)
   IF (itrlr1(1) < 0) CYCLE
   CALL sofcls
   
!     CALCULATE MATRIX REQUIRED
   
!                                T
!        **          **   **   ** **          ** **   **
!        *            *   *     * *            * *     *
!        * (M,B,K4)HH * = * HGH * * (M,B,K4)AA * * HGH *
!        *            *   *     * *            * *     *
!        **          **   **   ** **          ** **   **
   
!                         T
!        **   **   **   ** **   **
!        *     *   *     * *     *
!        * PHH * = * HGH * * PAA *
!        *     *   *     * *     *
!        **   **   **   ** **   **
   
   itrlr2(1) = hgh
   CALL rdtrl (itrlr2)
   item = itmlst(i)
   IF (i == 5 .AND. popt == papp) item = itmlst(6)
   DO  j = 1,7
     jtrlra(j) = itrlr2(j)
     jtrlrb(j) = itrlr1(j)
     jtrlre(j) = 0
   END DO
   iform = 6
   iprc  = 1
   ityp  = 0
   IF (jtrlra(5) == 2 .OR. jtrlra(5) == 4) iprc = 2
   IF (jtrlrb(5) == 2 .OR. jtrlrb(5) == 4) iprc = 2
   IF (jtrlra(5) >= 3 .OR. jtrlrb(5) >= 3) ityp = 2
   itype = iprc + ityp
   CALL makmcb (jtrlrc,otfile(i),itrlr2(2),iform,itype)
   jscr(1) = iscr(9)
   jscr(2) = iscr(2)
   jscr(3) = iscr(1)
   icode   = 0
   IF (i == 5) icode = 1
   prec3 = 0
   CALL mpy3dr (dz(dblkor))
   CALL wrttrl (jtrlrc)
   
!     STORE MATRIX ON SOF
!     I = 2, STORE MHH AS MMTX
!     I = 3, STORE BHH AS BMTX
!     I = 4, STORE K4HH AS K4MX
!     I = 5, STORE PHH AS PVEC OR PAPP
   
   CALL sofopn (z(sbuf1),z(sbuf2),z(sbuf3))
   itmnam(1) = newnam(1)
   itmnam(2) = newnam(2)
   CALL mtrxo (otfile(i),newnam,item,0,itest)
   IF (itest /= 3) GO TO 200
 END DO
 
!     TEST FOR LOAD PROCESSING
 
 IF (popt == blanks) GO TO 190
 IF (.NOT. ponly) GO TO 184
 itrlr1(1) = eqst
 CALL rdtrl (itrlr1)
 nsub(1 )  = itrlr1(6)
 nsub(2)   = itrlr1(7)
 item      = itmlst(11)
 itmnam(1) = oldnam(1)
 itmnam(2) = oldnam(2)
 CALL mtrxi (uprt,oldnam,item,0,itest)
 IF (itest /= 1) GO TO 200
 GO TO 188
 
!     PARTITION PAA VECTOR
 
 184 lcore = korlen
 fuset = usetmr
 CALL calcv (uprt,un,ui,ub,z(korbgn))
 188 CONTINUE
 CALL gmprtn (paa,pove,0,0,0,0,uprt,nsub(1),nsub(2),z(korbgn), korlen)
 
!     SAVE POVE AS POVE OR POAP ON SOF
 
 IF (modes) GO TO 190
 item = itmlst(7)
 IF (popt == papp) item = itmlst(10)
 CALL mtrxo (pove,oldnam,item,0,itest)
 IF (itest /= 3) GO TO 200
 190 CONTINUE
 GO TO 300
 
!     PROCESS MODULE ERRORS
 
 200 SELECT CASE ( itest )
   CASE (    1)
     GO TO 210
   CASE (    2)
     GO TO 220
   CASE (    3)
     GO TO 230
   CASE (    4)
     GO TO 240
   CASE (    5)
     GO TO 250
   CASE (    6)
     GO TO 270
 END SELECT
 210 imsg = -9
 GO TO 280
 220 imsg = -11
 GO TO 280
 230 imsg = -1
 GO TO 260
 240 imsg = -2
 GO TO 260
 250 imsg = -3
 260 CALL smsg (imsg,item,itmnam)
 GO TO 300
 270 imsg = -10
 280 dry  = -2
 CALL smsg1 (imsg,item,itmnam,modnam)
 GO TO 300
 290 WRITE  (iprntr,295) ufm
 295 FORMAT (a23,' 6518, ONE OF THE COMPONENT SUBSTRUCTURES HAS BEEN ',  &
     'USED IN A PREVIOUS COMBINE OR REDUCE.')
 dry = -2
 300 RETURN
END SUBROUTINE mred2g
