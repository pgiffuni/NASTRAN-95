SUBROUTINE adr
     
    !     AERODYNAMIC DATA RECOVERY   -  FORCE OUTPUT BY SET SELECTION
 
    !     DMAP
    !     FLUTTER
    !     ADR  CPHIH1,CASEZZ,QKHL,CLAMAL1,SPLINE,SILA,USETA/PKF/C,N,BOV/C,
    !          N,MACH=0.0/C,N,APP $
    !     DYNAMICS
    !     ADR  UHVT1,CASECC,QKHL,TOL1,SPLINE,SILA,USETA/PKF/V,N,BOV/C,Y,
    !          MACH=0.0/C,N,APP $

    IMPLICIT NONE
    INTEGER :: sysbuf,out   ,casecc,disp  ,qkhl  ,load  ,spline,sila  ,  &
               useta ,iz(1) ,pkf   ,app   ,flut  ,freq  ,scr1  ,scr2  ,  &
               scr3  ,scr4  ,mcb(7),i     ,iaero ,ibuf1 ,ibuf2 ,ibuf3 ,  &
               ii    ,incr  ,incr1 ,inn   ,iout  ,ipa   ,ipd   ,ipq   ,  &
               iti   ,ito   ,j     ,k     ,l     ,m     ,nam   ,ncol  ,  &
               ncore ,next  ,nfreq ,nload ,nn    ,nnn   ,nns1  ,nogo  ,  &
               nrow  ,nterma,ntermd,nterms,nw

    REAL    :: mach  ,bov   ,pi    ,twopi ,z

    CHARACTER (LEN=29) :: uim
    CHARACTER (LEN=25) :: uwm
    CHARACTER (LEN=23) :: ufm

    COMMON /xmssg / ufm   ,uwm   ,uim
    COMMON /blank / bov   ,mach  ,app
    COMMON /system/ sysbuf,out
    COMMON /condas/ pi    ,twopi
    COMMON /unpakx/ iout  ,inn   ,nnn   ,incr1
    COMMON /packx / iti   ,ito   ,ii    ,nn    ,incr
    COMMON /zzzzzz/ z(1)

    EQUIVALENCE     (z(1),iz(1))

    DATA    iaero / 176/
    DATA    flut  / 4HFLUT/, freq /4HFREQ/
    DATA    disp  / 101/ , casecc /102/ , qkhl /103/ , load /104/
    DATA    spline/ 105/ , sila   /106/ , useta/107/ , pkf  /201/
    DATA    scr1  / 301/ , scr2   /302/ , scr3 /303/ , scr4 /304/
 
 
    !     BUILD    P    =  Q    *  U
    !               KF      KH      H
    !     WHERE  QKH INTERPOLATED FOR A EIGENVALUE OR FREQUENCY - MACH DEP.
    !            UH  - EIGENVALUE OR FREQUENCY
 
 
    !     INITIALIZE  - LOOK FOR A REQUEST
 
    IF (app == flut .OR. app == freq) GO TO 5
    GO TO 1000
5   ncore  = korsz(z)
    ibuf1  = ncore - sysbuf
    CALL OPEN (*1000,casecc,iz(ibuf1),0)
    CALL fwdrec (*1000,casecc)
    CALL READ (*1000,*10,casecc,z,ibuf1,0,nw)
10  IF (iz(iaero) == 0) GO TO 1000
    CALL CLOSE (casecc,1)
 
    !     BUILD INTERPOLATED MATRIX FROM QHKL ON SCR1
    !     DEPENDENT LIST
    !     IF CLAMAL1 PICK UP FREQUENCY FROM OFP TABLE
    !     IF TOL1    PICK UP FREQUENCY FROM HEADER
    !     INDEPENDENT LIST ON QKHL
 
    CALL OPEN (*1000,load,iz(ibuf1),0)
    IF (app == flut) GO TO 30
 
    !     TOL1 = LOAD
 
    mcb(1) = casecc
    CALL rdtrl (mcb)
    CALL READ (*1000,*1000,load,iz,-2,0,nfreq)
    CALL READ (*1000,*20,load,iz,ibuf1,0,nfreq)
    GO TO 999
20  nload = mcb(2)
    GO TO 60
 
    !     CLAMAL1 = LOAD
 
30  CALL fwdrec (*1000,load)
    CALL fwdrec (*1000,load)
    CALL READ (*1000,*40,load,iz,ibuf1,0,nfreq)
    GO TO 999
40  nfreq = nfreq/6
    IF(bov == 0.0) GO TO 997
    DO  i = 1,nfreq
        k     =  i*6 - 1
        z(i)  = z(k)/(twopi*bov)
    END DO
    nload = 1
 
    !     CALL ADRI TO BUILD  (AFTER ADRI FREQUENCY*2PI*BOV IS IN Z AT EVERY
    !     OTHER SLOT 0.0 ,W FOR NFREQ*2
 
60  CALL CLOSE (load,1)
    CALL adri (z,nfreq,ncore,qkhl,scr1,scr2,scr3,scr4,nrow,ncol,nogo)
    IF (nogo /= 0) GO TO 1000
 
    !     SCR1 NOW HAS QKH INTERPOLATED    NROW*NCOL(ROW5)  NFREQ(COLUMNS)
 
    ipq   = nfreq*2 + 1
 
    !     BUILD PKF
 
    iout  = 3
    iti   = 3
    ito   = 3
    incr  = 1
    incr1 = 1
    mcb(1)= disp
    CALL rdtrl (mcb)
    IF (mcb(1) <    0) GO TO 1000
    IF (mcb(3) /= ncol) GO TO 998
    nns1  = nrow*ncol
    ii    = 1
    nn    = nrow
    inn   = 1
    ibuf2 = ibuf1 - sysbuf
    CALL gopen (pkf,z(ibuf2),1)
    ibuf3 = ibuf2 - sysbuf
    CALL gopen (disp,z(ibuf3),0)
    CALL gopen (scr1,z(ibuf1),0)
    mcb(1) = pkf
    mcb(2) = 0
    mcb(3) = nn
    mcb(6) = 0
    mcb(7) = 0
    nterms = nns1*2
    ntermd = ncol*2
    nterma = nrow*2
    ipd    = ipq + nterms
    ipa    = ipd + ntermd
    next   = ipa + nterma
    IF (next > ibuf3) GO TO 999
    DO  i = 1,nload
        DO  j = 1,nfreq
     
            !     UNPACK INTERPOLATED MATRIX COLUMN THEN DISP VECTOR  MULTIPLY AND
            !     PACK OUT
     
            nnn = nns1
            CALL unpack (*70,scr1,z(ipq))
     
            !     MULTIPLY BACK BY FREQUENCY (K)
     
            DO  l = 1,nterms,2
                m = j*2
                z(ipq+l) = z(ipq+l)*z(m)
            END DO
            GO TO 75
70          CALL zeroc (z(ipq),nterms)
75          nnn = ncol
            CALL unpack (*80,disp,z(ipd))
            GO TO 90
80          CALL zeroc (z(ipd),ntermd)
90          CALL gmmatc (z(ipd),1,ncol,0,z(ipq),ncol,nrow,0,z(ipa))
            CALL pack (z(ipa),pkf,mcb)
        END DO
        IF (i == nload) CYCLE
        CALL REWIND (scr1)
        CALL skprec (scr1,1)
    END DO
    CALL CLOSE  (scr1,1)
    CALL CLOSE  (disp,1)
    CALL CLOSE  (pkf, 1)
    CALL wrttrl (mcb)
    CALL dmpfil (-pkf,z(ipq),ibuf3-ipq)
 
    !     PUT FREQUENCY BACK TO ORIGINAL VALUE
 
    DO  i = 1,nfreq
        z(i) = z(i*2)/(twopi*bov)
    END DO
 
    !     PRINT RESULTS
 
    CALL adrprt (casecc,pkf,spline,sila,useta,z,nfreq,ncore,nload)
 
    !     STOP  CLOSE ALL POSSIBLE OPENS
 
1000 CALL CLOSE (casecc,1)
    CALL CLOSE (load  ,1)
    CALL CLOSE (pkf   ,1)
    CALL CLOSE (disp  ,1)
    RETURN
 
    !     ERROR MESSAGES
 
999 CALL mesage (8,0,nam)
    GO TO 1000
998 CALL mesage (7,0,nam)
    GO TO 1000
997 WRITE  (out,9970) uim
9970 FORMAT (a29,' 2272, NO FLUTTER CALCULATIONS CAN BE MADE IN ',  &
        'MODULE ADR SINCE BOV = 0.0.')
    GO TO 1000
END SUBROUTINE adr
