SUBROUTINE cmrd2d (iter)
     
    !     THIS SUBROUTINE CALCULATES THE MODAL TRANSFORMATION MATRIX FOR THE
    !     CMRED2 MODULE.
 
    !     INPUT  DATA
    !     GINO - LAMAMR - EIGENVALUE TABLE FOR SUBSTRUCTURE BEING REDUCED
    !            PHISSR - RIGHT EIGENVECTOR MATRIX FOR SUBSTRUCTURE BEING
    !                     REDUCED
    !            PHISSL - LEFT EIGENVECTOR MATRIX FOR SUBSTRUCTURE BEING
    !                     REDUCED
    !     SOF  - GIMS   - G TRANSFORMATION MATRIX FOR ORIGINAL SUBSTRUCTURE
 
    !     OUTPUT DATA
    !     GINO - HIM    - MODAL TRANSFORMATION MATRIX
 
    !     PARAMETERS
    !     INPUT- GBUF   - GINO BUFFERS
    !            INFILE - INPUT FILE NUMBERS
    !            OTFILE - OUTPUT FILE NUMBERS
    !            ISCR   - SCRATCH FILE NUMBERS
    !            KORLEN - LENGTH OF OPEN CORE
    !            KORBGN - BEGINNING ADDRESS OF OPEN CORE
    !            OLDNAM - NAME OF SUBSTRUCTURE BEING REDUCED
    !            NMAX   - MAXIMUM NUMBER OF FREQUENCIES TO BE USED
    !     OUTPUT-MODUSE - BEGINNING ADDRESS OF MODE USE DESCRIPTION ARRAY
    !            NFOUND - NUMBER OF MODAL POINTS FOUND
    !            MODLEN - LENGTH OF MODE USE ARRAY
    !     OTHERS-HIMPRT - HIM PARTITION VECTOR
    !            PPRTN  - PHISS MATRIX PARTITION VECTOR
    !            PHIAM  - PHIAM MATRIX PARTITION
    !            PHIBM  - PHIBM MATRIX PARTITION
    !            PHIIM  - PHIIM MATRIX PARTITION
    !            IPARTN - BEGINNING ADDRESS OF PHISS PARTITION VECTOR
    !            LAMAMR - LAMAMR INPUT FILE NUMBER
    !            PHISS  - PHISS INPUT FILE NUMBER
    !            PPRTN  - PARTITION VECTOR FILE NUMBER
    !            HIMPRT - HIM PARTITION VECTOR FILE NUMBER
    !            GIB    - GIB INPUT FILE NUMBER
    !            PHIAM  - PHIAM PARTITION MATRIX FILE NUMBER
    !            PHIBM  - PHIBM PARTITION MATRIX FILE NUMBER
    !            PHIIM  - PHIIM PARTITION MATRIX FILE NUMBER
    !            HIM    - HIM INPUT FILE NUMBER
    !            HIMSCR - HIM SCRATCH INPUT FILE NUMBER
 
 
    INTEGER, INTENT(IN OUT)                  :: iter
    LOGICAL :: modes
    INTEGER :: dry,gbuf1,gbuf2,gbuf3,sbuf1,sbuf2,sbuf3,otfile,  &
        oldnam,z,typin,typep,fuset,t,signab,signc,prec, scr,un,ub,ui,rule,typeu,  &
        phiss,pprtn,gib,phiam,phibm,phiim,him,himprt,  &
        phissr,phissl,gibbar,himbar,himscr,usetmr,himtyp, dblkor,sglkor,dicore
    DOUBLE PRECISION :: dz,dhimsm,dhimag,dphim,dhimg
    DIMENSION       modnam(2),rz(1),itrlr(7),dz(1)
    CHARACTER (LEN=23) :: ufm
    COMMON /xmssg / ufm
    COMMON /BLANK / idum1,dry,idum6,gbuf1,gbuf2,gbuf3,sbuf1,sbuf2,  &
        sbuf3,infile(11),otfile(6),iscr(11),korlen,korbgn,  &
        oldnam(2),idum4(3),range(2),nmax,idum5,modes,  &
        idum8,moduse,nfound,modlen,idum9,lstzwd
    COMMON /zzzzzz/ z(1)
    COMMON /packx / typin,typep,irowp,nrowp,incrp
    COMMON /patx  / lcore,nsub(3),fuset
    COMMON /mpyadx/ itrlra(7),itrlrb(7),itrlrc(7),itrlrd(7),nz,t,  &
        signab,signc,prec,scr
    COMMON /bitpos/ idum3(9),un,idum7(10),ub,ui
    COMMON /parmeg/ ia(7),ia11(7),ia21(7),ia12(7),ia22(7),lcr,rule
    COMMON /unpakx/ typeu,irowu,nrowu,incru
    COMMON /system/ idum2,iprntr
    EQUIVALENCE     (lamamr,infile(2)),(phissr,infile(3)),  &
        (phissl,infile(4)),(usetmr,infile(6)),  &
        (phiam,iscr(8)),(himscr,iscr(7)),(phibm,iscr(9)),  &
        (gib,iscr(8)),(gibbar,iscr(11)),(phiim,iscr(6)),  &
        (himprt,iscr(7)),(himbar,iscr(8)),(pprtn,iscr(7)),  &
        (him,iscr(10)),(rz(1),z(1)),(dz(1),z(1))
    DATA    modnam/ 4HCMRD,4H2D  /
    DATA    epslon/ 1.0E-03/
    DATA    item  / 4HGIMS /
    DATA    iscr7 / 307    /
 
    !     READ LAMA FILE
 
    IF (dry == -2) RETURN
    kore  = korbgn
    ifile = lamamr
    CALL gopen (lamamr,z(gbuf1),0)
    CALL fwdrec (*170,lamamr)
    lamwds = 6
    IF (modes) lamwds = 7
    it = 0
2   CALL READ (*160,*4,lamamr,z(korbgn),lamwds,0,nwds)
    korbgn = korbgn + 6
    IF (korbgn >= korlen) GO TO 180
    it = it + 1
    GO TO 2
4   CALL CLOSE (lamamr,1)
 
    !     ZERO OUT PARTITIONING VECTOR AND SET UP MODE USE DESCRIPTION
    !     RECORD
 
    modext   = korbgn
    itrlr(1) = phissr
    IF (iter == 2) itrlr(1) = phissl
    CALL rdtrl (itrlr)
    itphis = itrlr(2)
    IF (3*itphis+modext >= korlen) GO TO 180
    lamlen = lamwds*itphis
    nnmax  = MIN0(nmax,itphis)
    moduse = modext + itphis
    ipartn = modext + 2*itphis
    modlen = itphis
    DO  i = 1,itphis
        z(moduse+i-1) = 3
        z(modext+i-1) = 0
        rz(ipartn+i-1) = 0.0
    END DO
 
    !     SELECT DESIRED MODES
 
    korbgn = modext + 3*itphis
    nfound = 0
    DO  i = 1,itphis
        IF (nfound == nnmax) GO TO 30
        j = 3 + lamwds*(i-1)
        IF (rz(kore+j) <= range(1) .OR. rz(kore+j) >= range(2)) CYCLE
        z(modext+nfound) = i
        nfound = nfound + 1
        z(moduse+i-1)  = 1
        rz(ipartn+i-1) = 1.0
    END DO
 
    !     PACK OUT PARTITIONING VECTOR
 
30  typin = 1
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
    korbgn = korbgn - itphis
 
    !     PARTITION PHISS(R,L) MATRICES
 
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
    phiss   = phissr
    IF (iter == 2) phiss = phissl
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
 
    fuset = usetmr
    CALL calcv (pprtn,un,ui,ub,z(korbgn))
    CALL gmprtn (phiam,phiim,phibm,0,0,0,pprtn,nsub(1),nsub(2), z(korbgn),icore)
    khim = 0
    IF (ia21(6) == 0) GO TO 55
 
    !     COMPUTE MODAL TRANSFORMATION MATRIX
 
    !        **   **   **     **   **   ** **     **
    !        *     *   *       *   *     * *       *
    !        * HIM * = * PHIIM * - * GIB * * PHIBM *
    !        *     *   *       *   *     * *       *
    !        **   **   **     **   **   ** **     **
 
    IF (iter == 2) GO TO 40
    CALL softrl (oldnam,item,itrlr)
    itest = itrlr(1)
    IF (itest /= 1) GO TO 200
    CALL mtrxi (gib,oldnam,item,0,itest)
    IF (itest /= 1) GO TO 200
    itrlr(1) = gib
    GO TO 45
40  itrlr(1) = gibbar
    CALL rdtrl (itrlr)
    45 DO  i = 1, 7
        itrlra(i) = itrlr(i)
        itrlrb(i) = ia21(i)
        itrlrc(i) = ia11(i)
    END DO
    iform = 2
    iprc  = 1
    ityp  = 0
    IF (itrlra(5) == 2 .OR. itrlra(5) == 4) iprc = 2
    IF (itrlrb(5) == 2 .OR. itrlrb(5) == 4) iprc = 2
    IF (itrlrc(5) == 2 .OR. itrlrc(5) == 4) iprc = 2
    IF (itrlra(5) >= 3) ityp = 2
    IF (itrlrb(5) >= 3) ityp = 2
    IF (itrlrc(5) >= 3) ityp = 2
    itype = iprc + ityp
    CALL makmcb (itrlrd,himscr,itrlr(3),iform,itype)
    CALL sofcls
    t      = 0
    signab =-1
    signc  = 1
    prec   = 0
    scr    = iscr(7)
    dblkor = korbgn/2 + 1
    nz     = lstzwd - 2*dblkor - 1
    CALL mpyad  (dz(dblkor),dz(dblkor),dz(dblkor))
    CALL wrttrl (itrlrd)
    CALL sofopn (z(sbuf1),z(sbuf2),z(sbuf3))
    i      = itrlrd(2)
    ii     = itrlrd(3)
    iform  = itrlrd(4)
    himtyp = itrlrd(5)
    GO TO 60
 
    !     PHIBM IS NULL, HIM = PHIIM
 
55  himscr = phiim
    i      = ia11(2)
    ii     = ia11(3)
    iform  = ia11(4)
    himtyp = ia11(5)
    khim   = 1
    dblkor = korbgn/2 + 1
 
    !     TEST SELECTED MODES
 
60  ncore = 4*ii
    IF (khim == 0) ncore = ncore + 4*ia11(3)
    IF (korbgn+ncore >= korlen) GO TO 180
    typin = himtyp
    typep = himtyp
    irowp = 1
    nrowp = ii
    incrp = 1
    irowu = 1
    jhim  = him
    IF (iter == 2) jhim = himbar
    CALL gopen (himscr,z(gbuf1),0)
    IF (khim == 0) CALL gopen (phiim,z(gbuf2),0)
    CALL makmcb (itrlr,jhim,ii,iform,himtyp)
    CALL gopen  (jhim,z(gbuf3),1)
    nfound = 0
    it     = i
    dblkor = korbgn/2 + 1
    sglkor = 2*dblkor - 1
    IF (himtyp == 3) dicore = ((sglkor + 2*ii)/2) + 1
    IF (himtyp == 4) dicore = dblkor + 2*ii
    icore = 2*dicore - 1
 
    !     UNPACK HIM AND PHIIM COLUMNS
 
    DO  i = 1,it
        typeu = himtyp
        incru = 1
        nrowu = ii
        ihim  = nrowu
        CALL unpack (*110,himscr,dz(dblkor))
        IF (khim == 1) GO TO 70
        typeu = ia11(5)
        incru = 1
        nrowu = ia11(3)
        iphim = nrowu
        CALL unpack (*90,phiim,dz(dicore))
   
        !     SAVE LARGEST HIM COLUMN VALUE AND CALCULATE MAGNITUDE OF HIM,
        !     PHIIM COLUMNS
   
70      IF (himtyp == 4) GO TO 74
        itype  = 0
        himsum = 0.0
        himmag = 0.0
        DO  j = 1,ihim
            k = 1 + 2*(j-1)
            himag = SQRT((rz(sglkor+k-1)**2) + (rz(sglkor+k)**2))
            IF (himag >= himmag) himmag = himag
            himsum = himsum + (rz(sglkor+k-1)**2) + (rz(sglkor+k)**2)
        END DO
        GO TO 78
74      itype  = 1
        dhimsm = 0.0D0
        dhimag = 0.0D0
        DO  j = 1,ihim
            k = 1 + 2*(j-1)
            dhimg = DSQRT((dz(dblkor+k-1)**2) + (dz(dblkor+k)**2))
            IF (dhimg >= dhimag) dhimag = dhimg
            dhimsm = dhimsm + (dz(dblkor+k-1)**2) + (dz(dblkor+k)**2)
        END DO
78      IF (khim    == 1) GO TO 95
        IF (ia11(5) == 4) GO TO 82
        itype  = itype + 1
        phimsm = 0.0
        DO  j = 1,iphim
            k = 1 + 2*(j-1)
            phimsm = phimsm + (rz(icore+k-1)**2) + (rz(icore+k)**2)
        END DO
        GO TO 85
82      itype = itype + 2
        dphim = 0.0D0
        DO  j = 1,iphim
            k = 1 + 2*(j-1)
            dphim = dphim + (dz(dicore+k-1)**2) + (dz(dicore+k)**2)
        END DO
   
        !     TEST FOR INCLUSION
   
85      SELECT CASE ( itype )
            CASE (    1)
                GO TO 86
            CASE (    2)
                GO TO 87
            CASE (    3)
                GO TO 88
            CASE (    4)
                GO TO 89
        END SELECT
86      IF (phimsm == 0.0) GO TO 90
        IF (SQRT(himsum)/SQRT(phimsm) >= epslon) GO TO 95
        GO TO 90
87      IF (dphim == 0.0) GO TO 90
        IF (SQRT(himsum)/DSQRT(dphim) >= epslon) GO TO 95
        GO TO 90
88      IF (phimsm == 0.0) GO TO 90
        IF (DSQRT(dhimsm)/SQRT(phimsm) >= epslon) GO TO 95
        GO TO 90
89      IF (dphim == 0.0D0) GO TO 90
        IF (DSQRT(dhimsm)/DSQRT(dphim) >= epslon) GO TO 95
   
        !     REJECT MODE
   
90      j = z(modext+i-1)
        z(moduse+j-1) = 2
        CYCLE
   
        !     USE MODE
   
95      nfound = nfound + 1
   
        !     SCALE HIM COLUMN
   
        ihim = 2*ihim
        IF (himtyp == 4) GO TO 104
        DO  j = 1,ihim
            rz(sglkor+j-1) = rz(sglkor+j-1)/himmag
        END DO
        GO TO 130
        104 DO  j = 1,ihim
            dz(dblkor+j-1) = dz(dblkor+j-1)/dhimag
        END DO
        GO TO 130
   
        !     NULL COLUMN
   
110     ihim = 2*ihim
        IF (himtyp == 4) GO TO 114
        DO  j = 1,ihim
            rz(sglkor+j-1) = 0.0
        END DO
        GO TO 130
        114 DO  j = 1,ihim
            dz(dblkor+j-1) = 0.0D0
        END DO
   
        !     PACK HIM COLUMN
   
130     nrowp = nrowu
        CALL pack (dz(dblkor),jhim,itrlr)
    END DO
    CALL CLOSE (jhim,1)
    IF (khim == 0) CALL CLOSE (phiim,1)
    CALL CLOSE (himscr,1)
    CALL wrttrl (itrlr)
    korbgn = kore
    IF (khim == 1) himscr = iscr7
    RETURN
 
    !     PROCESS SYSTEM FATAL ERRORS
 
160 imsg = -2
    GO TO 190
170 imsg = -3
    GO TO 190
180 imsg = -8
    ifile = 0
190 CALL sofcls
    CALL mesage (imsg,ifile,modnam)
    RETURN
 
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
210 WRITE (iprntr,900) ufm,modnam,item,oldnam
    dry = -2
    RETURN
 
220 imsg = -1
    GO TO 250
230 imsg = -2
    GO TO 250
240 imsg = -3
250 CALL smsg (imsg,item,oldnam)
    RETURN
 
260 WRITE (iprntr,901) ufm,modnam,item,oldnam
    dry = -2
    RETURN
 
900 FORMAT (a23,' 6215, MODULE ',2A4,' - ITEM ',a4,  &
        ' OF SUBSTRUCTURE ',2A4,' PSEUDO-EXISTS ONLY.')
901 FORMAT (a23,' 6632, MODULE ',2A4,' - NASTRAN MATRIX FILE FOR I/O',  &
        ' OF SOF ITEM ',a4,', SUBSTRUCRURE ',2A4,', IS PURGED.')
 
END SUBROUTINE cmrd2d
