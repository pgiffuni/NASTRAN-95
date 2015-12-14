SUBROUTINE cmrd2e (iter)
     
    !     THIS SUBROUTINE CALCULATES THE H TRANSFORMATION MATRIX FOR THE
    !     CMRED2 MODULE.
 
    !     INPUT  DATA
    !     GINO - HIM     - MODAL TRANSFORMATION MATRIX
    !     SOF  - GIMS    - G TRANSFORMATION MATRIX FOR BOUNDARY POINTS OF
    !                      ORIGINAL SUBSTRUCTURE
 
    !     OUTPUT  DATA
    !     GINO  - HGH    - HORG PARTITION MATRIX
    !     SOF   - HORG   - H TRANSFORMATION MATRIX FOR ORIGINAL SUBSTRUCTURE
 
    !     PARAMETERS
    !     INPUT - GBUF   - GINO BUFFERS
    !             INFILE - INPUT FILE NUMBERS
    !             OTFILE - OUTPUT FILE NUMBERS
    !             ISCR   - SCRATCH FILE NUMBERS
    !             KORLEN - LENGTH OF OPEN CORE
    !             KORBGN - BEGINNING ADDRESS OF OPEN CORE
    !             OLDNAM - NAME OF SUBSTRUCTURE BEING REDUCED
    !     OTHERS- HIM    - HIM PARTITION MATRIX FILE NUMBER (RIGHT SIDE)
    !             HGH    - HORG MATRIX FILE NUMBER (RIGHT SIDE)
    !             GIB    - GIMS INPUT FILE NUMBER (RIGHT SIDE)
    !             HIMBAR - HIM PARTITION MATRIX FILE NUMBER (LEFT SIDE)
    !             HGHBAR - HGH PARTITION MATRIX FILE NUMBER (LEFT SIDE)
    !             GIBBAR - GIB PARTITION MATRIX FILE NUMBER (LEFT SIDE)
    !             UPRT   - USET PARTITIONING VECTOR FILE NUMBER
 
 
    INTEGER, INTENT(IN OUT)                  :: iter
    INTEGER :: dry,gbuf1,gbuf2,z,typinp,typeop,typeu,him,hgh,  &
        gib,himbar,hghbar,gibbar,uprt,himrl,hghrl,gibrl,  &
        dblkor,gibtyp,himtyp,sglkor,dicore
    DOUBLE PRECISION :: dz
    DIMENSION        modnam(2),itrlr1(7),itrlr2(7),rz(1),itmlst(4),  &
        dz(1),itrlr3(7)
    CHARACTER (LEN=23) :: ufm
    COMMON /xmssg /  ufm
    COMMON /BLANK /  idum1,dry,idum7,gbuf1,gbuf2,idum2(4),infile(11),  &
        otfile(6),iscr(11),korlen,korbgn,oldnam(2)
    COMMON /zzzzzz/  z(1)
    COMMON /packx /  typinp,typeop,irowp,nrowp,incrp
    COMMON /unpakx/  typeu,irowu,nrowu,incru
    COMMON /system/  idum3,iprntr
    EQUIVALENCE      (him,iscr(10)),(gib,iscr(6)),(uprt,iscr(7)),  &
        (gibbar,iscr(11)),(hghbar,iscr(9)),(hgh,iscr(9)),  &
        (himbar,iscr(8)),(rz(1),z(1)),(dz(1),z(1))
    DATA    modnam/  4HCMRD,4H2E  /
    DATA    itmlst/  4HHORG,4HHLFT,4HGIMS,4HUPRT/
 
    !     SET UP ROW PARTITION
 
    IF (dry == -2) RETURN
    item = itmlst(4)
    CALL mtrxi (uprt,oldnam,item,0,itest)
    IF (itest /= 1) GO TO 210
    CALL softrl (oldnam,item,itrlr1)
    IF (korbgn+itrlr1(3) >= korlen) GO TO 270
    typeu = itrlr1(5)
    irowu = 1
    nrowu = itrlr1(3)
    incru = 1
    CALL gopen (uprt,z(gbuf1),0)
    CALL unpack (*5,uprt,rz(korbgn))
    GO TO 15
    5 DO  i = 1, nrowu
        rz(korbgn+i-1) = 0.0
    END DO
15  CALL CLOSE (uprt,1)
    luprt  = nrowu
    kore   = korbgn
    korbgn = korbgn + luprt
 
    !     GET GIB MATRIX
 
    IF (iter == 2) GO TO 20
    item = itmlst(3)
    CALL softrl (oldnam,item,itrlr1)
    IF (itest /= 1) GO TO 210
    CALL mtrxi (gib,oldnam,item,0,itest)
    IF (itest /= 1) GO TO 210
    itrlr1(1) = gib
    gibrl = gib
    GO TO 30
20  itrlr1(1) = gibbar
    CALL rdtrl (itrlr1)
    gibrl = gibbar
 
    !     SET UP HGH TRAILER
 
30  hghrl = hgh
    IF (iter == 2) hghrl = hghbar
    nrows1 = itrlr1(3)
    kols1  = itrlr1(2)
    gibtyp = itrlr1(5)
    himrl  = him
    IF (iter == 2) himrl = himbar
    itrlr2(1) = himrl
    CALL rdtrl (itrlr2)
    nrows2 = itrlr2(3)
    kols2  = itrlr2(2)
    himtyp = itrlr2(5)
    iform  = 2
    IF (itrlr1(2)+itrlr1(3) == itrlr2(2)+itrlr2(3)) iform = 1
    iprc   = 1
    ityp   = 0
    IF (itrlr1(5) == 2 .OR. itrlr1(5) == 4) iprc = 2
    IF (itrlr2(5) == 2 .OR. itrlr2(5) == 4) iprc = 2
    IF (itrlr1(5) >= 3) ityp = 2
    IF (itrlr2(5) >= 3) ityp = 2
    itype = iprc + ityp
    CALL makmcb (itrlr3,hghrl,luprt,iform,itype)
 
    !     SET UP PACK/UNPACK PARAMETERS
 
    typeop = itrlr3(5)
    irowp  = 1
    nrowp  = itrlr1(2) + itrlr1(3)
    incrp  = 1
    incru  = 1
    dblkor = korbgn/2 + 1
    sglkor = 2*dblkor - 1
 
    !     FORM HGH MATRIX
 
    !                  **         **
    !                  *     .     *
    !        **   **   *  I  .  0  *
    !        *     *   *     .     *
    !        * HGH * = *...........*
    !        *     *   *     .     *
    !        **   **   * GIB . HIM *
    !                  *     .     *
    !                  **         **
 
    CALL gopen (hghrl,z(gbuf1),1)
 
    !     PROCESS GIB MATRIX
 
    typeu  = itrlr1(5)
    nrowu  = itrlr1(3)
    typinp = itrlr1(5)
    nrows  = itrlr1(3)
    IF (itrlr1(5) > 2) nrows = 2*itrlr1(3)
    IF (itrlr1(5) == 1 .OR. itrlr1(5) == 3) dicore = (sglkor+nrows)/2 + 1
    IF (itrlr1(5) == 2 .OR. itrlr1(5) == 4) dicore = dblkor + nrows
    icore  = 2*dicore - 1
    IF (dicore+nrows >= korlen) GO TO 270
    CALL gopen (gibrl,z(gbuf2),0)
    DO  i = 1,kols1
        k  = 0
        kk = 0
        CALL unpack (*40,gibrl,dz(dblkor))
        GO TO 50
   
        !     NULL GIB COLUMN
   
40      SELECT CASE ( gibtyp )
            CASE (    1)
                GO TO 42
            CASE (    2)
                GO TO 46
            CASE (    3)
                GO TO 42
            CASE (    4)
                GO TO 46
        END SELECT
        42 DO  j = 1,nrows
            rz(sglkor+j-1) = 0.0
        END DO
        GO TO 50
        46 DO  j = 1,nrows
            dz(dblkor+j-1) = 0.0D0
        END DO
   
        !     MOVE GIB DATA
   
        50 DO  j = 1,luprt
            IF (rz(kore+j-1) == 1.0) GO TO 70
            kk = kk + 1
            l  = 1 + 2*(kk-1)
            ll = 1 + 2*( j-1)
            SELECT CASE ( gibtyp )
                CASE (    1)
                    GO TO 62
                CASE (    2)
                    GO TO 64
                CASE (    3)
                    GO TO 66
                CASE (    4)
                    GO TO 68
            END SELECT
62          rz(icore+j-1) = rz(sglkor+kk-1)
            CYCLE
64          dz(dicore+j-1) = dz(dblkor+kk-1)
            CYCLE
66          rz(icore+ll-1) = rz(sglkor+l-1)
            rz(icore+ll  ) = rz(sglkor+l)
            CYCLE
68          dz(dicore+ll-1) = dz(dblkor+l-1)
            dz(dicore+ll  ) = dz(dblkor+l)
            CYCLE
     
            !     MOVE IDENTITY MATRIX DATA
     
70          k = k + 1
            l = 1 + 2*(j-1)
            SELECT CASE ( gibtyp )
                CASE (    1)
                    GO TO 72
                CASE (    2)
                    GO TO 74
                CASE (    3)
                    GO TO 76
                CASE (    4)
                    GO TO 78
            END SELECT
72          rz(icore+j-1) = 0.0
            IF (k == i) rz(icore+j-1) = 1.0
            CYCLE
74          dz(dicore+j-1) = 0.0D0
            IF (k == i) dz(dicore+j-1) = 1.0D0
            CYCLE
76          rz(icore+l-1) = 0.0
            IF (k == i) rz(icore+l-1) = 1.0
            rz(icore+l) = 0.0
            CYCLE
78          dz(dicore+l-1) = 0.0D0
            IF (k == i) dz(dicore+l-1) = 1.0D0
            dz(dicore+l) = 0.0D0
        END DO
        CALL pack (dz(dicore),hghrl,itrlr3)
    END DO
    CALL CLOSE (gibrl,1)
 
    !     PROCESS HIM MATRIX
 
    typeu  = itrlr2(5)
    nrowu  = itrlr2(3)
    typinp = itrlr2(5)
    nrows  = itrlr2(3)
    IF (itrlr2(5) > 2) nrows = 2*itrlr2(3)
    IF (itrlr2(5) == 2 .OR. itrlr2(5) == 4) dicore = (sglkor+nrows)/2 + 1
    IF (itrlr2(5) == 1 .OR. itrlr2(5) == 3) dicore = dblkor + nrows
    icore = 2*dicore - 1
    IF (dicore+nrows >= korlen) GO TO 270
    CALL gopen (himrl,z(gbuf2),0)
    DO  i = 1,kols2
        kk = 0
        CALL unpack (*100,himrl,dz(dblkor))
        GO TO 110
   
        !     NULL HIM COLUMN
   
100     SELECT CASE ( himtyp )
            CASE (    1)
                GO TO 102
            CASE (    2)
                GO TO 106
            CASE (    3)
                GO TO 102
            CASE (    4)
                GO TO 106
        END SELECT
        102 DO  j = 1,nrows
            rz(sglkor+j-1) = 0.0
        END DO
        GO TO 110
        106 DO  j = 1,nrows
            dz(dblkor+j-1) = 0.0D0
        END DO
   
        !     MOVE HIM MATRIX DATA
   
        110 DO  j = 1,luprt
            IF (rz(kore+j-1) == 1.0) GO TO 130
            kk = kk + 1
            l  = 1 + 2*(kk-1)
            ll = 1 + 2*( j-1)
            SELECT CASE ( himtyp )
                CASE (    1)
                    GO TO 122
                CASE (    2)
                    GO TO 124
                CASE (    3)
                    GO TO 126
                CASE (    4)
                    GO TO 128
            END SELECT
122         rz(icore+j-1) = rz(sglkor+kk-1)
            CYCLE
124         dz(dicore+j-1) = dz(dblkor+kk-1)
            CYCLE
126         rz(icore+ll-1) = rz(sglkor+l-1)
            rz(icore+ll  ) = rz(sglkor+l)
            CYCLE
128         dz(dicore+ll-1) = dz(dblkor+l-1)
            dz(dicore+ll  ) = dz(dblkor+l)
            CYCLE
     
            !     MOVE ZERO MATRIX DATA
     
130         l = 1 + 2*(j-1)
            SELECT CASE ( himtyp )
                CASE (    1)
                    GO TO 132
                CASE (    2)
                    GO TO 134
                CASE (    3)
                    GO TO 136
                CASE (    4)
                    GO TO 138
            END SELECT
132         rz(icore+j-1) = 0.0
            CYCLE
134         dz(dicore+j-1) = 0.0D0
            CYCLE
136         rz(icore+l-1) = 0.0
            rz(icore+l  ) = 0.0
            CYCLE
138         dz(dicore+l-1) = 0.0D0
            dz(dicore+l  ) = 0.0D0
        END DO
        CALL pack (dz(dicore),hghrl,itrlr3)
    END DO
    CALL CLOSE (himrl,1)
    CALL CLOSE (hghrl,1)
    CALL wrttrl (itrlr3)
    korbgn = kore
 
    !     SAVE HGH ON SOF AS H(ORG,LFT) MATRIX
 
    item = itmlst(iter)
    CALL mtrxo (hghrl,oldnam,item,0,itest)
    IF (itest /= 3) GO TO 210
    RETURN
 
    !     PROCESS MODULE FATAL ERRORS
 
210 SELECT CASE ( itest )
        CASE (    1)
            GO TO 220
        CASE (    2)
            GO TO 220
        CASE (    3)
            GO TO 220
        CASE (    4)
            GO TO 230
        CASE (    5)
            GO TO 240
        CASE (    6)
            GO TO 260
    END SELECT
220 WRITE (iprntr,900) ufm,modnam,item,oldnam
    dry = -2
    RETURN
 
230 imsg = -2
    GO TO 250
240 imsg = -3
250 CALL smsg(imsg,item,oldnam)
    RETURN
 
260 WRITE (iprntr,901) ufm,modnam,item,oldnam
    dry = -2
    RETURN
 
    !     PROCESS SYSTEM FATAL ERRORS
 
270 imsg = -8
    ifile = 0
    CALL sofcls
    CALL mesage (imsg,ifile,modnam)
    RETURN
 
900 FORMAT (a23,' 6211, MODULE ',2A4,' - ITEM ',a4,  &
        ' OF SUBSTRUCTURE ',2A4,' HAS ALREADY BEEN WRITTEN.')
901 FORMAT (a23,' 6632, MODULE ',2A4,' - NASTRAN MATRIX FILE FOR I/O',  &
        ' OF SOF ITEM ',a4,', SUBSTRUCTURE ',2A4,', IS PURGED.')
 
END SUBROUTINE cmrd2e
