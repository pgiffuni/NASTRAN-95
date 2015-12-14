SUBROUTINE mred2f
     
!     THIS SUBROUTINE COMPUTES THE FREEBODY EFFECTS FOR THE MRED2
!     MODULE.
 
!     INPUT DATA
!     GINO   - MAA    - SUBSTRUCTURE MASS MATRIX
!              DMR    - FREEBODY MATRIX
!     SOF    - GIMS   - G TRANSFORMATION MATRIX FOR BOUNDARY POINTS OF
!                       ORIGINAL SUBSTRUCTURE
 
!     OUTPUT DATA
!     GINO   - HGH    - HORG PARTITION MATRIX
!     SOF    - HORG   - H TRANSFORMATION MATRIX FOR ORIG. SUBSTRUCTURE
 
!     PARAMETERS
!     INPUT  - GBUF   - GINO BUFFERS
!              INFILE - INPUT FILE NUMBERS
!              ISCR   - SCRATCH FILE NUMBERS
!              KORLEN - LENGTH OF OPEN CORE
!              KORBGN - BEGINNING ADDRESS OF OPEN CORE
!              FREBDY - FREEBODY MODES OPTION FLAG
!     OTHERS - RPRTN  - ROW PARTITIONING VECTOR FILE NUMBER
!              LII    - LII PARTITION MATRIX FILE NUMBER (ISCR11)
!              IDENT  - IDENTITY MATRIX FILE NUMBER
!              ZERO   - ZERO MATRIX FILE NUMBER
!              HIE    - HIE PARTITION MATRIX FILE NUMBER
!              HIR    - HIR PARTITION MATRIX FILE NUMBER
!              HIRSCR - HIR SCRATCH PARTITION MATRIX FILE NUMBER
!              FBR    - FBR PARTITION MATRIX FILE NUMBER
!              FIR    - FIR PARTITION MATRIX FILE NUMBER
!              GIB    - GIMS INPUT FILE NUMBER
!              CPRTN  - COLUMN PARTITIONING VECTOR FILE NUMBER
!              HIM    - HIM PARTITION MATRIX FILE NUMBER
!              HGH    - HORG MATRIX FILE NUMBER
 
 LOGICAL :: frebdy,bounds
 INTEGER :: dry,gbuf1,gbuf2,sbuf1,sbuf2,sbuf3,oldnam,z,t,  &
     signab,signc,precmp,scr,fuset,precfb,SIGN,typinp,  &
     typeop,typinu,un,ub,ui,dmr,far,fir,farind,zero,  &
     rprtn,hie,hir,cprtn,him,hgh,gib,hirscr,usetmr, dblkor,sglkor
 DOUBLE           PRECISION dz,dhirmg
 DIMENSION modnam (2),itrlr1(7),itrlr2(7),rz(1),isub(4),itmlst(4), dz(1)
 COMMON /BLANK /  idum1,dry,idum7,gbuf1,gbuf2,idum2,sbuf1,sbuf2,  &
     sbuf3,infile(12),otfile(6),iscr(10),korlen,  &
     korbgn,oldnam(2),idum4(2),frebdy,idum8(5),bounds, idum9(6),lstzwd,iscr11
 COMMON /zzzzzz/  z(1)
 COMMON /mpyadx/  itrlra(7),itrlrb(7),itrlrc(7),itrlrd(7),nzmpy,t,  &
     signab,signc,precmp,scr
 COMMON /bitpos/  idum5(9),un,idum6(10),ub,ui
 COMMON /patx  /  lcore,nsub(3),fuset
 COMMON /fbsx  /  jtrlrl(7),jtrlru(7),jtrlrb(7),jtrlrx(7),nzfbs, precfb,SIGN
 COMMON /packx /  typinp,typeop,irowp,nrowp,incrp
 COMMON /unpakx/  typinu,irowu,nrowu,incru
 COMMON /system/  idum3,iprntr
 EQUIVALENCE      (usetmr,infile(5)),(maa,infile(7)),  &
     (dmr,infile(11)),(rprtn,iscr(9)),(ident,iscr(5)),  &
     (cprtn,iscr(10)),(pprtn,iscr(4)),(rz(1),z(1)),  &
     (dz(1),z(1)),(gib,iscr(4)),(lii ,iscr11),  &
     (hirscr,iscr(5)),(hgh,iscr(8)),(zero,iscr(6)),  &
     (him,iscr(8)),(hie,iscr(7)),(hir ,iscr(9)), (far,iscr(9)),(fir,iscr(10))
 DATA    modnam/  4HMRED,4H2F  /
 DATA    farind,  iscr7 ,iscr8 /6, 307, 308  /
 DATA    itmlst/  4HGIMS,4HHORG,4HUPRT,4HLMTX/
 
!     TEST FREEBODY MODES CALCULATION FLAG
 
 IF (dry == -2) GO TO 300
 itrlr2(1) = dmr
 CALL rdtrl (itrlr2)
 IF (itrlr2(1) < 0) GO TO 110
 
!     COMPUTE FREEBODY MATRIX
 
!        **   **   **   ** **   **
!        *     *   *     * *     *
!        * FAR * = * MAA * * DMR *
!        *     *   *     * *     *
!        **   **   **   ** **   **
 
 CALL sofcls
 frebdy = .true.
 itrlr1(1) = maa
 CALL rdtrl (itrlr1)
 DO  i = 1,7
   itrlra(i) = itrlr1(i)
   itrlrb(i) = itrlr2(i)
   itrlrc(i) = 0
 END DO
 iform = 2
 iprc  = 1
 ityp  = 0
 IF (itrlra(5) == 2 .OR. itrlra(5) == 4) iprc = 2
 IF (itrlrb(5) == 2 .OR. itrlrb(5) == 4) iprc = 2
 IF (itrlra(5) >= 3) ityp = 2
 IF (itrlrb(5) >= 3) ityp = 2
 itype = iprc + ityp
 CALL makmcb (itrlrd,far,itrlr1(3),iform,itype)
 t     = 0
 signab= 1
 signc = 1
 prec  = 0
 scr   = iscr(4)
 dblkor= 1 + korbgn/2
 nzmpy = lstzwd - 2*dblkor - 1
 CALL mpyad  (dz(dblkor),dz(dblkor),dz(dblkor))
 CALL wrttrl (itrlrd)
 
!     PARTITION FAR INTO BOUNDARY, INTERIOR POINTS
 
!                  **   **
!                  *     *
!        **   **   * FBR *
!        *     *   *     *
!        * FAR * = *.....*
!        *     *   *     *
!        **   **   * FIR *
!                  *     *
!                  **   **
 
 lcore = nzmpy
 fuset = usetmr
 CALL calcv (pprtn,un,ui,ub,z(korbgn))
 CALL gmprtn (far,fir,0,0,0,0,pprtn,nsub(1),nsub(2),z(korbgn), korlen)
 
!     CALCULATE FREEBODY TRANSFORMATION MATRIX
 
!                       T
!        **   ** **   ** **   **    **   **
!        *     * *     * *     *    *     *
!        * LII * * LII * * HIR * = -* FIR *
!        *     * *     * *     *    *     *
!        **   ** **   ** **   **    **   **
 
 IF (.NOT.bounds) GO TO 20
 item = itmlst(4)
 CALL softrl (oldnam,item,jtrlrl)
 itest = jtrlrl(1)
 IF (itest /= 1) GO TO 20
 jtrlrl(1) = lii
 CALL sofopn (z(sbuf1),z(sbuf2),z(sbuf3))
 CALL mtrxi (lii,oldnam,item,0,itest)
 IF (itest /= 1) GO TO 210
 CALL sofcls
 GO TO 30
 20 jtrlrl(1) = lii
 CALL rdtrl (jtrlrl)
 30 jtrlrb(1) = fir
 CALL rdtrl (jtrlrb)
 iform = 2
 iprc  = 1
 ityp  = 0
 IF (jtrlrl(5) == 2 .OR. jtrlrl(5) == 4) iprc = 2
 IF (jtrlrb(5) == 2 .OR. jtrlrb(5) == 4) iprc = 2
 IF (jtrlrl(5) >= 3) ityp = 2
 IF (jtrlrb(5) >= 3) ityp = 2
 itype = iprc + ityp
 CALL makmcb (jtrlrx,hir,jtrlrb(3),iform,itype)
 nzfbs  = nzmpy
 precfb = itype
 SIGN   = -1
 CALL fbs (z(korbgn),z(korbgn))
 CALL wrttrl (jtrlrx)
 
!     UNPACK HIR COLUMNS FOR SCALING
 
 typinu = jtrlrx(5)
 irowu  = 1
 nrowu  = jtrlrx(3)
 incru  = jtrlrx(5)
 typinp = jtrlrx(5)
 typeop = jtrlrx(5)
 irowp  = 1
 nrowp  = jtrlrx(3)
 incrp  = jtrlrx(5)
 CALL gopen (hir,z(gbuf1),0)
 iform  = jtrlrx(4)
 CALL makmcb (itrlr1,hirscr,jtrlrx(3),iform,jtrlrx(5))
 CALL gopen (hirscr,z(gbuf2),1)
 sglkor = 2*dblkor - 1
 DO  i = 1,farind
   CALL unpack (*60,hir,dz(dblkor))
   
!     CALCULATE MAGNITUDE OF HIR
   
   IF (jtrlrx(5) == 2) GO TO 42
   hirmag = rz(sglkor)
   IF (nrowu == 1) GO TO 50
   DO  j = 2,nrowu
     IF (ABS(rz(sglkor+j-1)) > ABS(hirmag)) hirmag = rz(sglkor+j-1)
   END DO
   GO TO 50
   42 dhirmg = dz(dblkor)
   IF (nrowu == 1) GO TO 50
   DO  j = 2,nrowu
     IF (DABS(dz(dblkor+j-1)) > DABS(dhirmg)) dhirmg =dz(dblkor+j-1)
   END DO
   
!     SCALE HIR COLUMN
   
   50 IF (jtrlrx(5) == 2) GO TO 54
   DO  j = 1,nrowu
     rz(sglkor+j-1) = rz(sglkor+j-1)/hirmag
   END DO
   GO TO 80
   54 DO  j = 1,nrowu
     dz(dblkor+j-1) = dz(dblkor+j-1)/dhirmg
   END DO
   GO TO 80
   
!     NULL COLUMN
   
   60 IF (jtrlrx(5) == 2) GO TO 74
   DO  j = 1,nrowu
     rz(sglkor+j-1) = 0.0
   END DO
   GO TO 80
   74 DO  j = 1,nrowu
     dz(dblkor+j-1) = 0.0D0
   END DO
   
!     PACK HIR COLUMN
   
   80 CALL pack (dz(dblkor),hirscr,itrlr1)
 END DO
 CALL CLOSE (hirscr,1)
 CALL CLOSE (hir,1)
 CALL wrttrl (itrlr1)
 isub(1) = itrlr1(2)
 
!     SET UP MERGE COLUMN PARTITION VECTOR
 
 itrlr2(1) = him
 CALL rdtrl (itrlr2)
 i = itrlr1(2) + itrlr2(2)
 isub(2) = itrlr2(2)
 DO  j = 1,i
   rz(korbgn+j-1) = 0.0
   IF (j > isub(1)) rz(korbgn+j-1) = 1.0
 END DO
 typinp = 1
 typeop = 1
 irowp  = 1
 nrowp  = i
 incrp  = 1
 iform  = 7
 CALL makmcb (itrlr2,rprtn,nrowp,iform,typinp)
 CALL gopen (rprtn,z(gbuf1),1)
 CALL pack(rz (korbgn),rprtn,itrlr2)
 CALL CLOSE (rprtn,1)
 CALL wrttrl (itrlr2)
 
!     MERGE FREEBODY, MODAL TRANSFORMATION MATRICES
 
!        **   **   **         **
!        *     *   *     .     *
!        * HIE * = * HIR . HIM *
!        *     *   *     .     *
!        **   **   **         **
 
 IF (hie /= him) GO TO 105
 hie = iscr8
 hgh = iscr7
 105 itype = 1
 IF (i /= itrlr2(3)) itype = 2
 CALL gmmerg (hie,hirscr,0,him,0,rprtn,0,isub,itype,z(korbgn), korlen)
 CALL sofopn (z(sbuf1),z(sbuf2),z(sbuf3))
 GO TO 120
 
!     FREEBODY MODES NOT REQUESTED
 
 110 hie = him
 IF (hie == iscr7) hgh = iscr8
 IF (hie == iscr8) hgh = iscr7
 
!     FORM HGH MATRIX
 
!                  **         **
!                  *     .     *
!        **   **   *  I  .  0  *
!        *     *   *     .     *
!        * HGH * = *...........*
!        *     *   *     .     *
!        **   **   * GIB . HIE *
!                  *     .     *
!                  **         **
 
 120 CALL softrl (oldnam,itmlst(2),itrlr1)
 IF (itrlr1(1) == 1) GO TO 190
 
!     GENERATE IDENTITY MATRIX
 
 CALL softrl (oldnam,itmlst(1),itrlr1)
 itest = itrlr1(1)
 item  = itmlst(1)
 IF (itest /= 1) GO TO 210
 typinp = 1
 typeop = itrlr1(5)
 irowp  = 1
 nrowp  = itrlr1(2)
 incrp  = 1
 iform  = 8
 ii = itrlr1(2)
 CALL makmcb (itrlr1,ident,nrowp,iform,typeop)
 CALL gopen  (ident,z(gbuf1),1)
 DO  i = 1,ii
   DO  j = 1,ii
     rz(korbgn+j-1) = 0.0
     IF (i == j) rz(korbgn+j-1) = 1.0
   END DO
   CALL pack (rz(korbgn),ident,itrlr1)
 END DO
 CALL CLOSE (ident,1)
 CALL wrttrl (itrlr1)
 
!     SET UP MERGE ROW PARTITION VECTOR
 
 itrlr1(1) = hie
 CALL rdtrl (itrlr1)
 iter  = itrlr1(2)
 nrowp = ii + iter
 DO  i = 1,nrowp
   rz(korbgn+i-1) = 0.0
   IF (i > ii) rz(korbgn+i-1) = 1.0
 END DO
 typinp = 1
 typeop = 1
 incrp  = 1
 iform  = 7
 CALL makmcb (itrlr2,rprtn,nrowp,iform,typinp)
 CALL gopen (rprtn,z(gbuf1),1)
 CALL pack (rz(korbgn),rprtn,itrlr2)
 CALL CLOSE (rprtn,1)
 CALL wrttrl (itrlr2)
 nrows = nrowp
 
!     SET UP MERGE COLUMN PARTITION VECTOR
 
 item = itmlst(3)
 CALL mtrxi (cprtn,oldnam,item,0,itest)
 IF (itest /= 1) GO TO 210
 
!     SET UP GIB MATRIX
 
 CALL mtrxi (gib,oldnam,itmlst(1),0,itest)
 item = itmlst(1)
 IF (itest /= 1) GO TO 210
 
!     MERGE ALL STRUCTURAL REDUCTION TRANSFORMATION MATRICES
 
 isub(1) = ii
 isub(2) = iter
 isub(3) = itrlr1(3)
 isub(4) = ii
 itype   = 1
 IF (nrows /= nrowp) itype = 2
 CALL gmmerg (hgh,gib,ident,hie,0,rprtn,cprtn,isub,itype,z(korbgn), korlen)
 
!     SAVE HGH ON SOF AS HORG MATRIX
 
 CALL mtrxo (hgh,oldnam,itmlst(2),0,itest)
 item = itmlst(2)
 IF (itest /= 3) GO TO 210
 190 CONTINUE
 GO TO 300
 
!     PROCESS MODULE FATAL ERRORS
 
 210 SELECT CASE ( itest )
   CASE (    1)
     GO TO 220
   CASE (    2)
     GO TO 230
   CASE (    3)
     GO TO 240
   CASE (    4)
     GO TO 250
   CASE (    5)
     GO TO 260
   CASE (    6)
     GO TO 280
 END SELECT
 220 imsg = -9
 GO TO 290
 230 imsg = -11
 GO TO 290
 240 imsg = -1
 GO TO 270
 250 imsg = -2
 GO TO 270
 260 imsg = -3
 270 CALL smsg (imsg,item,oldnam)
 GO TO 300
 280 imsg = -10
 290 dry  = -2
 CALL smsg1 (imsg,item,oldnam,modnam)
 300 RETURN
END SUBROUTINE mred2f
