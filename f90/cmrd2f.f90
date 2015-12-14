SUBROUTINE cmrd2f (kode)
     
    !     THIS SUBROUTINE CALCULATES THE FINAL STRUCTURAL MATRICES FOR THE
    !     CMRED2 MODULE.
 
    !     INPUT  DATA
    !     GINO - KBB    - STIFFNESS PARTITION MATRIX
    !            KIB    - KIB STIFFNESS PATTITION MATRIX
    !            HIE    - HIE PARTITION MATRIX
    !            KII    - KII PARTITION MATRIX
    !            HGH    - HORG PARTITION MATRIX
    !            MAA    - MASS INPUT MATRIX
    !            BAA    - DAMPING INPUT MATRIX
    !            K4AA   - STIFFNESS INPUT MATRIX
    !            PAA    - LOADS INPUT MATRIX
    !     SOF  - GIMS   - G TRANSFORMATION MATRIX
 
    !     OUTPUT DATA
    !     GINO - KHH    - STIFFNESS MATRIX
    !            MHH    - MASS MATRIX
    !            BHH    - DAMPING MATRIX
    !            K4HH   - K4HH MATRIX
    !            PHH    - LOADS MATRIX
    !     SOF  - KMTX   - STIFFNESS MATRIX
    !            MMTX   - MASS MATRIX
    !            PVEC   - LOADS MATRIX
    !            PAPP   - APPENDED LOADS MATRIX
    !            BMTX   - DAMPING MATRIX
    !            K4MX   - K4MX STIFFNESS MATRIX
 
    !     PARAMETERS
    !     INPUT- POPT   - LOADS OPTION FLAG
    !            GBUF   - GINO BUFFERS
    !            INFILE - INPUT FILE NUMBERS
    !            OTFILE - OUTPUT FILE NUMBERS
    !            ISCR   - SCRATCH FILE NUMBERS
    !            KORLEN - LENGTH OF OPEN CORE
    !            KORBGN - BEGINNING ADDRESS OF OPEN CORE
    !            OLDNAM - NAME OF SUBSTRUCTURE BEING REDUCED
    !     OTHERS-PAA    - LOADS INPUT FILE NUMBER
    !            KHH    - STIFFNESS OUTPUT FILE NUMBER
    !            POVE   - LOADS OUTPUT FILE NUMBER
    !            UPRT   - PARTITION VECTOR FILE NUMBER
    !            ZEROMB - ZERO PARTITION FILE NUMBER
    !            KBB    - KBB INPUT FILE NUMBER
    !            ZEROBM - ZERO PARTITION MATRIX
    !            KIB    - KIB INPUT FILE NUMBER
    !            KII    - KII INPUT FILE NUMBER
    !            KBARBB - KBARBB FILE NU BER
    !            GIB    - GIB INPUT FILE NUMBER
    !            KMM    - KMM FILE NUMBER
    !            HGH    - HORG INPUT FILE NUMBER
 
 
    INTEGER, INTENT(IN OUT)                  :: kode
    LOGICAL :: symtry,modes,ponly
    INTEGER :: dry,popt,gbuf1,sbuf1,sbuf2,sbuf3,otfile,oldnam,z,  &
        t,signab,signc,prec,scr,typin,typout,un,ub,ui,  &
        fuset,prec3,paa,him,pove,uprt,gib,gibbar,hghbar,  &
        hgh,usetmr,cmred2,papp,blanks,dblkor,himbar,eqst
    DOUBLE PRECISION :: dz
    DIMENSION       modnam(2),itrlr1(7),itrlr2(7),itrlr3(7),isub(4),  &
        itmlst(12),itmnam(2),rz(1),dz(1)
    CHARACTER (LEN=23) :: ufm
    COMMON /xmssg / ufm
    COMMON /BLANK / idum1,dry,popt,gbuf1,idum2(2),sbuf1,sbuf2,sbuf3,  &
        infile(11),otfile(6),iscr(11),korlen,korbgn,  &
        oldnam(2),newnam(2),symtry,idum6(4),modes, idum7(4),ponly,lstzwd
    COMMON /zzzzzz/ z(1)
    COMMON /mpyadx/ itrlra(7),itrlrb(7),itrlrc(7),itrlrd(7),nz,t,  &
        signab,signc,prec,scr
    COMMON /packx / typin,typout,irow,nrow,incr
    COMMON /bitpos/ idum4(9),un,idum5(10),ub,ui
    COMMON /patx  / lcore,nsub(3),fuset
    COMMON /system/ idum3,iprntr
    COMMON /mpy3tl/ jtrlra(7),jtrlrb(7),jtrlre(7),jtrlrc(7),jscr(3),  &
        lkore,icode,prec3
    EQUIVALENCE     (eqst,infile(5)),(usetmr,infile(6)),  &
        (paa,infile(11)),(khh,otfile(1)),(pove,otfile(6)),  &
        (kbb,iscr(1)),(kib,iscr(2)),(kii,iscr(4)),  &
        (him,iscr(10)),(uprt,iscr(1)),(himbar,iscr(8)),  &
        (kbarbb,iscr(5)),(kmm,iscr(6)),(gib,iscr(3)),  &
        (gibbar,iscr(11)),(hghbar,iscr(9)),(hgh,iscr(8)),  &
        (rprtn,iscr(1)),(rz(1),z(1)),(dz(1),z(1))
    DATA    modnam/ 4HCMRD,4H2F  /, papp  / 4HPAPP/, blanks/ 4H    /
    DATA    cmred2/ 26    /
    DATA    itmlst/ 4HKMTX,4HHORG,4HHLFT,4HMMTX,4HBMTX,4HK4MX,4HPVEC,  &
        4HPAPP,4HPOVE,4HGIMS,4HPOAP,4HUPRT/
 
    !     SELECT OPERATION MODE
 
    IF (dry == -2) RETURN
    IF (ponly .OR. dry == 0) GO TO 90
 
    !     SET UP NEW SUBSTRUCTURE
 
    IF (modes) GO TO 1
    numb = 1
    CALL setlvl (newnam,numb,oldnam,itest,cmred2)
    IF (itest == 8) GO TO 290
 
    !     CHECK FOR STIFFNESS MATRIX GENERATION
 
1   itrlr1(1) = khh
    CALL rdtrl (itrlr1)
    IF (itrlr1(1) < 0) GO TO 90
 
    !     FORM PRELIMINARY STIFFNESS CALCULATION
 
    !                                           T
    !        **      **   **   **   **        ** **   **
    !        *        *   *     *   *          * *     *
    !        * KBARBB * = * KBB * + * GIB(BAR) * * KIB *
    !        *        *   *     *   *          * *     *
    !        **      **   **   **   **        ** **   **
 
    itrlr1(1) = kbb
    CALL rdtrl (itrlr1)
    IF (symtry) GO TO 2
    itrlr2(1) = gibbar
    CALL rdtrl (itrlr2)
    GO TO 4
2   item  = itmlst(10)
    CALL softrl (oldnam,item,itrlr2)
    itest = itrlr2(1)
    itmnam(1) = oldnam(1)
    itmnam(2) = oldnam(2)
    IF (itest /= 1) GO TO 200
    CALL mtrxi (gib,oldnam,item,0,itest)
    IF (itest /= 1) GO TO 200
    itrlr2(1) = gib
4   itrlr3(1) = kib
    CALL rdtrl (itrlr3)
    DO  i = 1,7
        itrlra(i) = itrlr2(i)
        itrlrb(i) = itrlr3(i)
        itrlrc(i) = itrlr1(i)
    END DO
    iform = 1
    iprc  = 1
    ityp  = 0
    IF (itrlra(5) == 2 .OR. itrlra(5) == 4) iprc = 2
    IF (itrlrb(5) == 2 .OR. itrlrb(5) == 4) iprc = 2
    IF (itrlrc(5) == 2 .OR. itrlrc(5) == 4) iprc = 2
    IF (itrlra(5) >= 3) ityp = 2
    IF (itrlrb(5) >= 3) ityp = 2
    IF (itrlrc(5) >= 3) ityp = 2
    itype  = iprc + ityp
    CALL makmcb (itrlrd,kbarbb,itrlr1(3),iform,itype)
    t      = 1
    signab = 1
    signc  = 1
    prec   = 0
    scr    = iscr(7)
    scr    = iscr(1)
    CALL sofcls
    dblkor = korbgn/2 + 1
    nz     = lstzwd - (2*dblkor - 1)
    CALL mpyad (dz(dblkor),dz(dblkor),dz(dblkor))
    CALL wrttrl (itrlrd)
    kbarow = itrlrd(3)
    kcol   = itrlrd(2)
 
    !     FORM PRELIMINARY STIFFNESS CALCULATION
 
    !                              T
    !        **   **   **        ** **   ** **   **
    !        *     *   *          * *     * *     *
    !        * KMM * = * HIM(BAR) * * KII * * HIM *
    !        *     *   *          * *     * *     *
    !        **   **   **        ** **   ** **   **
 
    itrlr1(1) = kii
    itrlr2(1) = him
    CALL rdtrl (itrlr1)
    CALL rdtrl (itrlr2)
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
    CALL makmcb (itrlrd,iscr(2),itrlr2(3),iform,itype)
    prec  = 0
    t     = 0
    signab= 1
    signc = 1
    scr   = iscr(1)
    CALL mpyad (dz(dblkor),dz(dblkor),dz(dblkor))
    CALL wrttrl (itrlrd)
    itrlr1(1) = him
    IF (.NOT.symtry) itrlr1(1) = himbar
    CALL rdtrl (itrlr1)
    DO  i = 1,7
        itrlra(i) = itrlr1(i)
        itrlrb(i) = itrlrd(i)
    END DO
    iform = 1
    iprc  = 1
    ityp  = 0
    IF (itrlra(5) == 2 .OR. itrlra(5) == 4) iprc = 2
    IF (itrlrb(5) == 2 .OR. itrlrb(5) == 4) iprc = 2
    IF (itrlra(5) >= 3) ityp = 2
    IF (itrlrb(5) >= 3) ityp = 2
    itype = iprc + ityp
    CALL makmcb (itrlrd,kmm,itrlr1(2),iform,itype)
    t    = 1
    prec = 0
    CALL mpyad (dz(dblkor),dz(dblkor),dz(dblkor))
    CALL wrttrl (itrlrd)
    kmmrow = itrlrd(3)
    kmmcol = itrlrd(2)
 
    !     GENERATE MERGE PARTITION VECTOR
 
    nrow = kcol + kmmcol
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
    CALL gopen (rprtn,z(gbuf1),1)
    CALL pack (rz(korbgn),rprtn,itrlr1)
    CALL CLOSE (rprtn,1)
    CALL wrttrl (itrlr1)
 
    !     FORM STIFFNESS MATRIX
 
    !                  **            **
    !                  *        .     *
    !        **   **   * KBARBB .  0  *
    !        *     *   *        .     *
    !        * KHH * = *..............*
    !        *     *   *        .     *
    !        **   **   *   0    . KMM *
    !                  *        .     *
    !                  **            **
 
    isub(1) = kcol
    isub(2) = kmmcol
    isub(3) = kbarow
    isub(4) = kmmrow
    itype   = 1
    CALL gmmerg (khh,kbarbb,0,0,kmm,rprtn,rprtn,isub,itype,z(korbgn), korlen)
 
    !     STORE KHH AS KMTX ON SOF
 
    CALL sofopn (z(sbuf1),z(sbuf2),z(sbuf3))
    CALL mtrxo (khh,newnam,itmlst(1),0,itest)
    item      = itmlst(1)
    itmnam(1) = newnam(1)
    itmnam(2) = newnam(2)
    IF (itest /= 3) GO TO 200
 
    !     LOCATE HGH MATRIX
    !       KODE .EQ. 0, BOTH HORG, HLFT ON SOF
    !       KODE .EQ. 1, HORG CALCULATED, HLFT ON SOF
    !       KODE .EQ. 2, HORG ON SOF, HLFT CALCULATED
    !       KODE .EQ. 3, BOTH HORG, HLFT CALCULATED
 
90  item      = itmlst(2)
    itmnam(1) = oldnam(1)
    itmnam(2) = oldnam(2)
    CALL mtrxi (hgh,oldnam,item,0,itest)
    IF (itest /= 1) GO TO 200
    IF (kode > 1 .OR. symtry) GO TO 100
    item = itmlst(3)
    CALL mtrxi (hghbar,oldnam,item,0,itest)
    IF (itest /= 1) GO TO 200
100 signab = 1
    signc  = 1
    scr    = iscr(1)
    dblkor = korbgn/2 + 1
    nz     = lstzwd - (2*dblkor - 1)
    itmnam(1) = newnam(1)
    itmnam(2) = newnam(2)
 
    !     GENERATE MATRICES REQUESTED
    !        I .EQ. 2, GENERATE MHH MATRIX
    !        I .EQ. 3, GENERATE BHH MATRIX
    !        I .EQ. 4, GENERATE K4HH MATRIX
    !        I .EQ. 5, GENERATE PHH MATRIX
 
    DO  i = 2,5
        itrlr1(1) = infile(i+6)
        CALL rdtrl (itrlr1)
        IF (itrlr1(1) < 0) CYCLE
        CALL sofcls
   
        !     CALCULATE MATRIX REQUIRED
   
        !                                     T
        !        **          **   **        ** **          ** **   **
        !        *            *   *          * *            * *     *
        !        * (M,B,K4)HH * = * HGH(BAR) * * (M,B,K4)AA * * HGH *
        !        *            *   *          * *            * *     *
        !        **          **   **        ** **          ** **   **
   
        !                              T
        !        **   **   **        ** **   **
        !        *     *   *          * *     *
        !        * PHH * = * HGH(BAR) * * PAA *
        !        *     *   *          * *     *
        !        **   **   **        ** **   **
   
        itrlr2(1) = hgh
        CALL rdtrl (itrlr2)
        IF (i == 5) GO TO 112
        DO  j = 1,7
            itrlra(j) = itrlr1(j)
            itrlrb(j) = itrlr2(j)
            itrlrc(j) = 0
        END DO
        iform = 2
        IF (itrlr1(3) == itrlr2(2)) iform = 1
        iprc  = 1
        ityp  = 0
        IF (itrlr1(5) == 2 .OR. itrlr1(5) == 4) iprc = 2
        IF (itrlr2(5) == 2 .OR. itrlr2(5) == 4) iprc = 2
        IF (itrlr1(5) >= 3) ityp = 2
        IF (itrlr2(5) >= 3) ityp = 2
        itype = iprc + ityp
        CALL makmcb (itrlrd,iscr(2),itrlr1(3),iform,itype)
        prec  = 0
        t     = 0
        signab= 1
        signc = 1
        scr   = iscr(1)
        CALL mpyad (dz(dblkor),dz(dblkor),dz(dblkor))
        CALL wrttrl (itrlrd)
        item = itmlst(i+2)
        GO TO 116
        112 DO  j = 1,7
            itrlrd(j) = itrlr1(j)
        END DO
        item = itmlst(7)
        IF (popt == papp) item = itmlst(8)
116     itrlr2(1) = hgh
        IF (.NOT. symtry) itrlr2(1) = hghbar
        CALL rdtrl (itrlr2)
        DO  j = 1,7
            itrlra(j) = itrlr2(j)
            itrlrb(j) = itrlrd(j)
        END DO
        iform = 1
        iprc  = 1
        ityp  = 0
        IF (itrlrd(5) == 2 .OR. itrlrd(5) == 4) iprc = 2
        IF (itrlr2(5) == 2 .OR. itrlr2(5) == 4) iprc = 2
        IF (itrlrd(5) >= 3) ityp = 2
        IF (itrlr2(5) >= 3) ityp = 2
        itype = iprc + ityp
        CALL makmcb (itrlrd,otfile(i),itrlr2(2),iform,itype)
        t    = 1
        prec = 0
        CALL mpyad (dz(dblkor),dz(dblkor),dz(dblkor))
        CALL wrttrl (itrlrd)
   
        !     STORE MATRIX ON SOF
        !        I .EQ. 2, STORE MHH AS MMTX
        !        I .EQ. 3, STORE BHH AS BMTX
        !        I .EQ. 4, STORE K4HH AS K4MX
        !        I .EQ. 5, STORE PHH AS PVEC OR PAPP
   
        CALL sofopn (z(sbuf1),z(sbuf2),z(sbuf3))
        CALL mtrxo (otfile(i),newnam,item,0,itest)
        IF (itest /= 3) GO TO 200
    END DO
 
    !     TEST FOR LOAD PROCESSING
 
    IF (popt == blanks) GO TO 190
    itmnam(1) = oldnam(1)
    itmnam(2) = oldnam(2)
    IF (.NOT.ponly) GO TO 184
    itrlr1(1) = eqst
    CALL rdtrl (itrlr1)
    nsub(1) = itrlr1(6)
    nsub(2) = itrlr1(7)
    item = itmlst(12)
    CALL mtrxi (uprt,oldnam,item,0,itest)
    IF (itest /= 1) GO TO 200
    GO TO 188
 
    !     PARTITION PAA VECTOR
 
184 lcore = korlen
    fuset = usetmr
    CALL calcv (uprt,un,ui,ub,z(korbgn))
188 CALL gmprtn (paa,pove,0,0,0,0,uprt,nsub(1),nsub(2),z(korbgn), korlen)
 
    !     SAVE POVE AS POVE OR POAP ON SOF
 
    IF (modes) GO TO 190
    item = itmlst(9)
    IF (popt == papp) item = itmlst(11)
    CALL mtrxo (pove,oldnam,item,0,itest)
    IF (itest /= 3) GO TO 200
190 CONTINUE
    RETURN
 
    !     PROCESS MODULE ERRORS
 
200 SELECT CASE ( itest )
        CASE (    1)
            GO TO 210
        CASE (    2)
            GO TO 210
        CASE (    3)
            GO TO 210
        CASE (    4)
            GO TO 220
        CASE (    5)
            GO TO 230
        CASE (    6)
            GO TO 250
    END SELECT
210 WRITE (iprntr,900) ufm,modnam,item,itmnam
    dry = -2
    RETURN
 
220 imsg = -2
    GO TO 240
230 imsg = -3
240 CALL smsg (imsg,item,itmnam)
    RETURN
 
250 WRITE (iprntr,901) ufm,modnam,item,itmnam
    dry = -2
    RETURN
 
290 WRITE (iprntr,902) ufm
    dry = -2
    RETURN
 
900 FORMAT (a23,' 6211, MODULE ',2A4,' - ITEM ',a4,  &
        ' OF SUBSTRUCTURE ',2A4,' HAS ALREADY BEEN WRITTEN.')
901 FORMAT (a23,' 6632, MODULE ',2A4,' - NASTRAN MATRIX FILE FOR I/O',  &
        ' OF SOF ITEM ',a4,', SUBSTRUCTURE ',2A4,', IS PURGED.')
902 FORMAT (a23,' 6518, ONE OF THE COMPONENT SUBSTRUCTURES HAS BEEN ',  &
        'USED IN A PREVIOUS COMBINE OR REDUCE.')
 
END SUBROUTINE cmrd2f
