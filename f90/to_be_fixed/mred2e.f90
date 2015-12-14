SUBROUTINE mred2e
     
!     THIS SUBROUTINE CALCULATES THE MODAL TRANSFORMATION MATRIX FOR THE
!     MRED2 MODULE.
 
!     INPUT DATA
!     GINO   - LAMAMR - EIGENVALUE TABLE FOR SUBSTRUCTURE BEING REDUCED
!              PHISS  - EIGENVECTOR MATRIX FOR SUBSTRUCTURE BEING REDUCE
!     SOF    - GIMS   - G TRANSFORMATION MATRIX FOR ORIGINAL SUBSTRUCTUR
 
!     OUTPUT DATA
!     GINO   - HIM    - HIM MATRIX PARTITION
 
!     PARAMETERS
!     INPUT  - GBUF   - GINO BUFFERS
!              INFILE - INPUT FILE NUMBERS
!              OTFILE - OUTPUT FILE NUMBERS
!              ISCR   - SCRATCH FILE NUMBERS
!              KORLEN - LENGTH OF OPEN CORE
!              KORBGN - BEGINNING ADDRESS OF OPEN CORE
!              OLDNAM - NAME OF SUBSTRUCTURE BEING REDUCED
!              NMAX   - MAXIMUM NUMBER OF FREQUENCIES TO BE USED
!     OUTPUT - MODUSE - BEGINNING ADDRESS OF MODE USE DESCRIPTION ARRAY
!              MODLEN - LENGTH OF MODE USE ARRAY
!              NFOUND - NUMBER OF MODAL POINTS FOUND
!     OTHERS - HIMPRT - HIM PARTITION VECTOR
!              PPRTN  - PHISS MATRIX PARTITION VECTOR
!              PHIAM  - PHIAM MATRIX PARTITION
!              PHIBM  - PHIBM MATRIX PARTITION
!              PHIIM  - PHIIM MATRIX PARTITION
!              IPARTN - BEGINNING ADDRESS OF PHISS PARTITION VECTOR
!              LAMAMR - LAMAMR INPUT FILE NUMBER
!              PHISS  - PHISS INPUT FILE NUMBER
!              PPRTN  - PARTITION VECTOR FILE NUMBER
!              HIMPRT - HIM PARTITION VECTOR FILE NUMBER
!              GIB    - GIB INPUT FILE NUMBER
!              PHIAM  - PHIAM PARTITION MATRIX FILE NUMBER
!              PHIBM  - PHIBM PARTITION MATRIX FILE NUMBER
!              PHIIM  - PHIIM PARTITION MATRIX FILE NUMBER
!              HIM    - HIM INPUT FILE NUMBER
!              HIMSCR - HIM SCRATCH INPUT FILE NUMBER
 
 LOGICAL :: frebdy
 INTEGER :: dry,gbuf1,gbuf2,gbuf3,sbuf1,sbuf2,sbuf3,otfile,  &
     oldnam,z,typin,typep,fuset,un,ub,ui
 INTEGER :: t,signab,signc,prec,scr,rule,typeu,phiss
 INTEGER :: pprtn,gib,phiam,phibm,phiim,him,himscr,himprt,  &
     usetmr,himtyp,fbmods,dblkor,sglkor,dicore
 DOUBLE PRECISION :: dz,dhimsm,dhimag
 DIMENSION        modnam(2),rz(1),itrlr(7),dz(1)
 CHARACTER (LEN=23) :: ufm
 COMMON /xmssg /  ufm
 COMMON /BLANK /  idum1,dry,idum6,gbuf1,gbuf2,gbuf3,sbuf1,sbuf2,  &
     sbuf3,infile(12),otfile(6),iscr(10),korlen,  &
     korbgn,oldnam(2),idum4(2),frebdy,range(2),nmax,  &
     idum5(5),moduse,nfound,modlen,idum2,lstzwd
 COMMON /zzzzzz/  z(1)
 COMMON /packx /  typin,typep,irowp,nrowp,incrp
 COMMON /patx  /  lcore,nsub(3),fuset
 COMMON /mpyadx/  itrlra(7),itrlrb(7),itrlrc(7),itrlrd(7),nz,t,  &
     signab,signc,prec,scr
 COMMON /bitpos/  idum3(9),un,idum7(10),ub,ui
 COMMON /parmeg/  ia(7),ia11(7),ia21(7),ia12(7),ia22(7),lcr,rule
 COMMON /unpakx/  typeu,irowu,nrowu,incru
 COMMON /system/  idum8,iprntr
 EQUIVALENCE      (lamamr,infile(2)),(phiss,infile(3)), (usetmr,infile(5))
 EQUIVALENCE      (gib,iscr(8)),(pprtn,iscr(5)), (him,iscr(8)),  &
     (himprt,iscr(9)),(phibm,iscr(9))
 EQUIVALENCE      (rz(1),z(1)),(dz(1),z(1))
 DATA    modnam/  4HMRED,4H2E  /
 DATA    epslon,  iscr4,fbmods /1.0E-03,304,6/
 DATA    item  /  4HGIMS       /
 
!     READ LAMAMR FILE
 
 IF (dry == -2) GO TO 300
 kore  = korbgn
 ifile = lamamr
 CALL gopen (lamamr,z(gbuf1),0)
 CALL fwdrec (*170,lamamr)
 iter = 0
 2 CALL READ (*160,*4,lamamr,z(korbgn),7,0,nwds)
 
!     REJECT MODES WITH NO ASSOCIATED VECTORS
 
 IF (rz(korbgn+5) <= 0.0) GO TO 2
 korbgn = korbgn + 7
 IF (korbgn >= korlen) GO TO 180
 iter = iter + 1
 GO TO 2
 4 CALL CLOSE (lamamr,1)
 
!     ZERO OUT PARTITIONING VECTOR AND SET UP MODE USE DESCRIPTION
!     RECORD
 
 modext  = korbgn
 itrlr(1)= phiss
 CALL rdtrl (itrlr)
 itphis  = itrlr(2)
 nrows   = itrlr(3)
 IF ((3*itphis)+modext >= korlen) GO TO 180
 lamlen  = 7*itphis
 nnmax   = MIN0(nmax,itphis)
 moduse  = modext + itphis
 ipartn  = modext + 2*itphis
 modlen  = itphis
 DO  i = 1,itphis
   z(modext+i-1) = 0
   z(moduse+i-1) = 3
   rz(ipartn+i-1) = 0.0
 END DO
 
!     SELECT DESIRED MODES
 
 korbgn = modext + 3*itphis
 IF (korbgn >= korlen) GO TO 180
 nfound = 0
 DO  i = 1,itphis
   j = 4 + 7*(i-1)
   IF (rz(kore+j) <= range(1) .OR. rz(kore+j) >= range(2)) CYCLE
   
!     REMOVE MODES WITH NEGATIVE EIGENVALUES
   
   IF (rz(kore+j-2) < 0.0) CYCLE
   z(modext+nfound) = i
   nfound = nfound + 1
   z(moduse +i-1) = 1
   rz(ipartn+i-1) = 1.0
 END DO
 
!     PACK OUT PARTITIONING VECTOR
 
 typin = 1
 typep = 1
 irowp = 1
 nrowp = itrlr(2)
 incrp = 1
 iform = 2
 CALL makmcb (itrlr,pprtn,nrowp,iform,typin)
 CALL gopen (pprtn,z(gbuf1),1)
 CALL pack (rz(ipartn),pprtn,itrlr)
 CALL CLOSE (pprtn,1)
 CALL wrttrl (itrlr)
 
!     PARTITION PHISS MATRIX
 
!        **     **   **         **
!        *       *   *   .       *
!        * PHISS * = * 0 . PHIAM *
!        *       *   *   .       *
!        **     **   **         **
 
 nsub(1) = itphis - nfound
 nsub(2) = nfound
 nsub(3) = 0
 lcore   = korlen - korbgn
 icore   = lcore
 
!     TEST FOR ALL MODES
 
 IF (nsub(1) == 0) GO TO 32
 phiam = iscr(8)
 CALL gmprtn (phiss,0,0,phiam,0,pprtn,0,nsub(1),nsub(2),z(korbgn), icore)
 
!     PARTITION PHIAM MATRIX
 
!                    **     **
!                    *       *
!        **     **   * PHIBM *
!        *       *   *       *
!        * PHIAM * = *.......*
!        *       *   *       *
!        **     **   * PHIIM *
!                    *       *
!                    **     **
 
 GO TO 34
 32 phiam = phiss
 34 CONTINUE
 
!     CALCULATE THE VECTOR MAGNITUDE
 
 IF (korbgn+nrows >= korlen) GO TO 180
 CALL gopen (phiam,z(gbuf1),0)
 typeu = 1
 irowu = 1
 nrowu = nrows
 incru = 1
 DO  i = 1,nfound
   l     = ipartn + i - 1
   rz(l) = 0.0
   CALL unpack (*40,phiam,rz(korbgn))
   DO  j = 1,nrows
     k     = korbgn + j - 1
     rz(l) = rz(l) + rz(k)**2
   END DO
 END DO
 CALL CLOSE (phiam,1)
 fuset = usetmr
 CALL calcv (pprtn,un,ui,ub,z(korbgn))
 
!     TEST FOR NULL B SET
 
 itrlr(1) = pprtn
 CALL rdtrl (itrlr)
 IF (itrlr(6) > 0) GO TO 44
 phiim   = phiam
 ia11(1) = phiam
 CALL rdtrl (ia11)
 DO  i = 1,7
   ia21(i) = 0
 END DO
 GO TO 55
 44 CONTINUE
 phiim = iscr(7)
 CALL gmprtn (phiam,phiim,phibm,0,0,0,pprtn,nsub(1),nsub(2), z(korbgn),icore)
 jhim = 0
 
!     COMPUTE MODAL TRANSFORMATION MATRIX
 
!        **   **   **     **   **   ** **     **
!        *     *   *       *   *     * *       *
!        * HIM * = * PHIIM * - * GIB * * PHIBM *
!        *     *   *       *   *     * *       *
!        **   **   **     **   **   ** **     **
 
 CALL mtrxi (gib,oldnam,item,0,itest)
 IF (itest /= 1) GO TO 200
 CALL softrl (oldnam,item,itrlr)
 itest = itrlr(1)
 IF (itest /= 1) GO TO 200
 DO  i = 1, 7
   itrlra(i) = itrlr(i)
   itrlrb(i) = ia21(i)
   itrlrc(i) = ia11(i)
 END DO
 itrlra(1) = gib
 himscr = iscr(4)
 iform  = 2
 iprc   = 1
 ityp   = 0
 IF (itrlra(5) == 2 .OR. itrlra(5) == 4) iprc = 2
 IF (itrlrb(5) == 2 .OR. itrlrb(5) == 4) iprc = 2
 IF (itrlrc(5) == 2 .OR. itrlrc(5) == 4) iprc = 2
 IF (itrlra(5) >= 3) ityp = 2
 IF (itrlrb(5) >= 3) ityp = 2
 IF (itrlrc(5) >= 3) ityp = 2
 itype  = iprc + ityp
 CALL makmcb (itrlrd,himscr,itrlr(3),iform,itype)
 CALL sofcls
 t      = 0
 signab = -1
 signc  = 1
 prec   = 0
 scr    = iscr(6)
 dblkor = korbgn/2 + 1
 nz     = lstzwd - ((2*dblkor)-1)
 CALL mpyad (dz(dblkor),dz(dblkor),dz(dblkor))
 CALL wrttrl (itrlrd)
 CALL sofopn (z(sbuf1),z(sbuf2),z(sbuf3))
 i      = itrlrd(2)
 ii     = itrlrd(3)
 iform  = itrlrd(4)
 himtyp = itrlrd(5)
 GO TO 60
 
!     PHIBM IS NULL, HIM = PHIIM
 
 55 himscr = phiim
 i      = ia11(2)
 ii     = ia11(3)
 iform  = ia11(4)
 himtyp = ia11(5)
 jhim   = 1
 
!     TEST SELECTED MODES
 
 60 ncore  = i
 IF (korbgn+ncore >= korlen) GO TO 180
 typin  = himtyp
 typep  = himtyp
 irowp  = 1
 nrowp  = ii
 incrp  = 1
 irowu  = 1
 CALL gopen (himscr,z(gbuf1),0)
 CALL makmcb (itrlr,him,ii,iform,himtyp)
 CALL gopen (him,z(gbuf3),1)
 nfound = 0
 iter   = i
 dblkor = korbgn/2 + 1
 sglkor = 2*dblkor - 1
 IF (himtyp == 1) dicore = (sglkor+ii)/2 + 1
 IF (himtyp == 2) dicore = dblkor + ii
 icore  = 2*dicore - 1
 
!     UNPACK HIM COLUMN
 
 DO  i = 1,iter
   
!     LIMIT VECTORS TO NMAX
   
   IF (nfound < nnmax) GO TO 65
   j      = z(modext+i-1) + moduse - 1
   z(j)   = 3
   CYCLE
   65 typeu  = himtyp
   incru  = 1
   nrowu  = ii
   ihim   = nrowu
   CALL unpack (*90,himscr,dz(dblkor))
   
!     SAVE LARGEST HIM COLUMN VALUE AND CALCULATE MAGNITUDE OF HIM,
!     COLUMN
   
   IF (himtyp == 2) GO TO 74
   itype  = 0
   himsum = 0.0
   himmag = 0.0
   DO  j = 1,ihim
     IF (ABS(rz(sglkor+j-1)) >= ABS(himmag)) himmag = rz(sglkor+j-1)
     himsum = himsum + (rz(sglkor+j-1)**2)
   END DO
   GO TO 78
   74 itype  = 2
   dhimsm = 0.0D0
   dhimag = 0.0D0
   DO  j = 1,ihim
     IF (DABS(dz(dblkor+j-1)) >= DABS(dhimag)) dhimag = dz(dblkor+j-1)
     dhimsm = dhimsm + dz(dblkor+j-1)**2
   END DO
   himsum = dhimsm
   78 IF (jhim == 1) GO TO 95
   phimsm = rz(ipartn+i-1)
   IF (phimsm <= 0.0) GO TO 90
   pmsm   = phimsm*epslon*epslon
   IF (himsum >= pmsm) GO TO 95
   
!     REJECT MODE
   
   90 j = z(modext+i-1)
   z(moduse+j-1) = 2
   CYCLE
   
!     USE MODE
   
   95 nfound = nfound + 1
   
!     SCALE HIM COLUMN
   
   IF (himtyp == 2) GO TO 104
   DO  j = 1,ihim
     rz(sglkor+j-1) = rz(sglkor+j-1)/himmag
   END DO
   GO TO 130
   104 DO  j = 1,ihim
     dz(dblkor+j-1) = dz(dblkor+j-1)/dhimag
   END DO
   
!     PACK HIM COLUMN
   
   130 nrowp = nrowu
   CALL pack(dz (dblkor),him,itrlr)
 END DO
 CALL CLOSE (him,1)
 IF (jhim == 0) CALL CLOSE (phiim,1)
 CALL CLOSE (himscr,1)
 CALL wrttrl (itrlr)
 korbgn = kore
 IF (jhim == 1) himscr = iscr4
 
!     TEST NUMBER OF MODAL POINTS
 
 modal = itrlr(2)
 IF (frebdy) modal = modal + fbmods
 IF (modal <= itrlr(3)) GO TO 300
 WRITE  (iprntr,145) ufm,oldnam,modal,itrlr(3)
 145 FORMAT (a23,' 6633, FOR SUBSTRUCTURE ',2A4,' THE TOTAL NUMBER OF',  &
     ' MODAL COORDINATES (',i8,1H), /30X,  &
     'IS LARGER THAN THE NUMBER OF INTERNAL DOF (',i8,2H).)
 dry = -2
 GO TO 300
 
!     PROCESS SYSTEM FATAL ERRORS
 
 160 imsg = -2
 GO TO 190
 170 imsg = -3
 GO TO 190
 180 imsg = -8
 ifile = 0
 190 CALL sofcls
 CALL mesage (imsg,ifile,modnam)
 GO TO 300
 
!     PROCESS MODULE FATAL ERRORS
 
 200 SELECT CASE ( itest )
   CASE (    1)
     GO TO 210
   CASE (    2)
     GO TO 210
   CASE (    3)
     GO TO 220
   CASE (    4)
     GO TO 230
   CASE (    5)
     GO TO 240
   CASE (    6)
     GO TO 260
 END SELECT
 210 imsg = -11
 GO TO 270
 220 imsg = -1
 GO TO 250
 230 imsg = -2
 GO TO 250
 240 imsg = -3
 250 CALL smsg (imsg,item,oldnam)
 GO TO 300
 260 imsg = -10
 270 CALL smsg1 (imsg,item,oldnam,modnam)
 dry = -2
 300 RETURN
END SUBROUTINE mred2e
