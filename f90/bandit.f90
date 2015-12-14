SUBROUTINE bandit
     
    !     BANDIT - A COMPUTER PROGRAM TO RE-SEQUENCE MATRIX BY BANDWIDTH,
    !              PROFILE, AND WAVEFRONT METHODS FOR NASTRAN.
 
    !     THIS PROGRAM GENERATES THE RE-SEQUENCE CARDS, SEQGP (AFTER GEOM1,
    !     GEOM2, AND GEOM4 DATA BLOCKS ARE ASSEMBLED), AND ADD THESE CARDS
    !     TO THE END OF GEOM1 FILE.
 
    !     HOWEVER, IF THE ORIGINAL NASTRAN INPUT DECK CONTAINS ONE OR MORE
    !     SEQGP CARD, BANDIT WILL BE AUTOMATICALLY SKIPPED.
 
    !     ******************************************************************
 
    !     ACKNOWLEDGEMENT:
 
    !     THE ORIGINAL BANDIT PROGRAM (VERSION 9, DEC. 1978, DISTRIBUTED BY
    !     COSMIC  NO. DOD-0034) WAS WRITTEN BY G. C. EVERTINE OF NAVAL SHIP
    !     RESEARCH AND DEVELOPMENT CENTER (NSRDC), BETHESDA, MD.
 
    !     THE FOLLOWING SUBROUTINES WERE WRITTEN BY E. CUTHILL AND J. MCKEE
    !     OF NSRDC
    !     - CTHMCK,DEGREE,DIAM,IDIST,KOMPNT,MAXDGR,MINDEG,RELABL
 
    !     THE FOLLOWING SUBROUTINES WERE WRITTEN BY N. GIBBS, W. POOLE,
    !     P. STOCKMEYER, AND H. CRANE OF THE COLLEGE OF WILLIAM AND MARY
    !     - DGREE,FNDIAM,GIBSTK,NUMBER,PIKLVL,RSETUP,SORTDG,SORT2,TREE.
    !     (THESE ROUTINES AND CTHMCK WERE MODIFIED BY G. C. EVERSTINE.)
 
    !     ******************************************************************
 
    !     ONLY HALF OF THE ORIGINAL BANDIT PROGRAM WAS ADOPTED IN THIS
    !     NASTRAN VERSION BY G. C. CHAN OF SPERRY, HUNTSVILLE, AL., 1982
 
    !     THE ORIGINAL BANDIT ROUTINES WERE UPDATED TO MEET NASTRAN
    !     PROGRAMMING STYLE AND STANDARD.
    !     NASTRAN GINO FILES AND GINO I/O ARE USED INSTEAD OF FORTRAN FILES
    !     AND FORTRAN READ/WRITE
    !     THE INTEGER PACK AND UNPACK ROUTINES, BPACK AND BUNPK, WERE RE-
    !     WRITTEN TO ALLOW COMMON USAGE FOR IBM, CDC, UNIVAC AND VAX MACH.
 
    !     ROUTINES BANDIT, SCHEME, BREAD, BGRID, BSEQGP, AND TIGER WERE
    !     COMPLETELY RE-WRITTEN.
    !     (SCHEME WAS FORMALLY CALLED NASNUM, AND CTHMCK WAS SCHEME)
 
    !     ******************************************************************
 
    !     THIS NASTRAN VERSION DOES NOT USE $-OPTION CARDS AS IN THE CASE OF
    !     ORIGINAL BANDIT PROGRAM.
 
    !     THE FOLLOWING 'OPTIONS' ARE PRE-SELECTED -
 
    !        $ADD        (NOT USE)            $INSERT     (NOT USE)
    !        $APPEND     (NOT USE)            $METHOD     (GPS    )
    !        $CONFIG     (NOT USE)            $MPC        (NO     )
    !        $CRITERION  (RMS    )            $NASTRAN    (NOT USE)
    !        $DEGREE     (NOT USE)            $PLUS       (NOT USE)
    !        $DIMENSION  (NOT USE)            $PRINT      (MIN    )
    !        $ELEMENTS   (NOT USE)            $PUNCH      (NONE   )
    !        $FRONTAL    (NOT USE)            $SEQUENCE   (YES    )
    !        $GRID       (NOT USE)            $SPRING     (NO     )
    !        $HICORE     (NOT USE)            $TABLE      (NO     )
    !        $IGNORE     (NOT USE)            $START      (NOT USE)
 
    !     ******************************************************************
 
    INTEGER :: hicore,   geom1,    geom2,    geom4,    scr1,  &
              z,        sub(3),   END, rd,       rdrew,    wrt,      wrtrew,   rew
 
    COMMON /machin/  machin
    COMMON /banda /  ibuf1,    nompc,    nodep,    nopch,    norun,  &
        method,   icrit,    ngpts(2)
    COMMON /bandb /  nbitin,   kore,     ifl,      ngrid,    ipass,  &
        nw,       kdim,     nbpw,     irept
    COMMON /bandd /  korig,    knew,     iop,      inp,      ncm,  &
        nzero,    nel,      neq,      neqr
    COMMON /bandg /  dum3g(3)
    COMMON /bands /  nn,       mm,       ih,       ib,       maxgrd,  &
        maxdeg,   kmod,     mach,     mindeg,   nedge, mask
    COMMON /bandw /  dum4w(4), i77,      dum2w(2)
    COMMON /system/  ibuf,     nout,     nogo,     is(97)
    COMMON /geomx /  geom1,    geom2,    geom4,    scr1
    COMMON /names /  rd,       rdrew,    wrt,      wrtrew,   rew
    COMMON /zzzzzz/  z(1)
 
    EQUIVALENCE      (hicore,is(28))
    DATA             sub  /    4HBAND,   4HIT  ,   4HBEGN  /
    DATA             END,      iquit /   4HEND ,   4HQUIT  /

    !     INITIALIZE PROGRAM PARAMETERS

    !     NOMPC =  0, MPC'S AND RIGID ELEM. ARE NOT USED IN BANDIT COMPUTATI
    !           = +1, ONLY RIGID ELEMENTS ARE USED IN BANDIT RESEQUENCING
    !           = +2, BOTH MPC'S  AND RIGID ELEM. ARE USED IN BANDIT
    !           = +3, ONLY MPC'S, NOT RIGID ELEM. ARE USED IN BANDIT
    !     NODEP = +1, MPC DEPENDENT PTS. ARE TO BE REMOVED FROM COMPUTATION
    !           = -1, MPC DEPENDENT PTS. ARE NOT TO BE REMOVED.
    !                 (NOTE - NODEP DICTATES ALSO THE DEPENDENT GRIDS OF
    !                         THE RIGID ELEMENTS)
    !     NOPCH = +1, PUNCH OUT SEQGP CARDS
    !           = -1, NO SEQGP CARDS PUNCHED
    !     NORUN = +1, BANDIT WILL RUN EVEN SEQGP CARDS ARE PRESENT
    !           = -1, BANDIT IS SKIPPED IF ONE OR MORE SEQGP CARD IS
    !                 PRESENT IN  THE INPUT DECK
    !     METHOD= -1, CM METHOD ONLY
    !           =  0, BOTHE CM AND GPS METHODS ARE USED
    !           = +1, USE GPS METHOD ONLY
    !     ICRIT =     RE-SEQUENCING CRITERION
    !           =  1, RMS WAVEFRONT
    !           =  2, BANDWIDTH
    !           =  3, PROFILE
    !           =  4, MAX WAVEFRONT

    nzero = 0
    i77   = 77
    nompc = 0
    nodep =-1
    nopch =-1
    norun =-1
    method=+1
    kdim  = 1
    icrit = 1
    irept = 0

    !     THE ABOVE DEFAULT VALUES CAN BE RESET BY THE NASTRAN CARD.
    !     (SEE SUBROUTINE NASCAR BANDIT FLAG FOR MORE DETAILS)
    !     ******************************************************************

    CALL conmsg (sub,3,0)
    nbpw = is(37)
    mach = machin
    kore = korsz(z(1))
    ibuf1= kore - ibuf - 2
    kore = ibuf1 - 1

    !     CALL BGRID TO GET THE NO. OF GRID POINTS IN THE PROBLEM, SET
    !     THE INTEGER PACKING CONSTANT, NW, AND COMPUTE MAXGRD AND MAXDEG.
    !     BANDIT QUITS IF PROBLEM IS TOO SMALL TO BE WORTHWHILE.

5   irept = irept + 1
    CALL bgrid
    IF (ngrid < 15) GO TO 30
    kdim4 = kdim*4
    ii3   = 2*maxgrd

    !     PARTITION OPEN CORE FOR SCHEME COMPUTATION.

    k2 =  1 + kdim4
    k3 = k2 + 2*ii3  + 2
    IF (method <= 0 .AND. maxdeg > maxgrd) k3 = k3 + maxdeg - maxgrd
    k4 = k3 + maxgrd + 1
    k5 = k4 + maxgrd
    k6 = k5 + maxgrd + 1
    k7 = k6 + maxgrd
    k8 = k7 + maxdeg
    k1 = k8 + maxdeg + nw
    k9 = k1 + maxgrd*maxdeg/nw
    IF (k9 > kore) CALL mesage (-8,k9-kore,sub)

    !     READ BULK DATA, SET UP CONNECTION TABLE, AND RESEQUENCE NODES.

    CALL scheme(z(k1),z(k2),ii3,z(k3),z(k4),z(k5),z(k6),z(k7),z(k8),z)
    IF (ngrid == -1) CALL sptchk
    IF (irept ==  2) GO TO 5
    IF (ngrid < 0) THEN
        GO TO    20
    ELSE IF (ngrid == 0) THEN
        GO TO    30
    END IF

    !     JOB DONE.

10  sub(3) = END
    GO TO 40

    !     NO BANDIT RUN.

20  nogo = 1
30  sub(3) = iquit
40  CALL conmsg (sub,3,0)

    RETURN
END SUBROUTINE bandit
