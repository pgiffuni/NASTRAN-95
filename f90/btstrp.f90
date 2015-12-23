SUBROUTINE btstrp
     
    !     BASED ON MACHINE NUMBER, BTSTRP WILL DEFINE ALL THE
    !     MACHINE-DEPENDENT CONSTANTS NEEDED IN NASTRAN. THESE CONSTANTS
    !     ARE SAVED IN LABEL COMMONS /SYSTEM/, /LHPWX/, /MACHIN/ & /CHMACH/
 
    !     SEE ALSO  PRTPRM, SDCMPS, SDR2E, AND UPCASE WHEN NASTRAN SOURCE
    !     CODE IS PORTED TO OTHER (NEW) MACHINE
 
    !                ===
    EXTERNAL lshift, rshift, andf, complf

    CHARACTER (LEN=1) :: mchnam*11, machos*7, comput(22)*11, compos(22)*7

    INTEGER           :: sysbuf, outtap ,  two   , complf ,rshift ,  &
                         fcb   , order  ,idate(3), abcd   ,ak     ,  &
                         recl  , andf   ,sperlk  , qp     ,m1(110),  &
                         m2(110) ,mconst(220), highpw
    REAL :: xx     ,yy

    COMMON /machin/ machx  ,ihalf  ,jhalf  ,lqro
    COMMON /chmach/ mchnam ,machos
    COMMON /sem   / a      ,mask2  ,mask3  ,lnknos(15)
    COMMON /lhpwx / lowpw  ,highpw ,nwpic  ,nudflw    ,mxfl, kshift ,mtisa
    COMMON /system/ b(100)
    COMMON /two   / two(32),mzero
    COMMON /xxread/ dum(3) ,ibmcdc

    EQUIVALENCE (b( 1),sysbuf) ,(b(22),linkno) ,(b(40),nbpw ) ,  &
                (b( 2),outtap) ,(b(41),ncpw ) , (b(42),idate(1)),  &
                (b( 4),intp  ) ,(b(34),idrum ) ,(b(55),iprec) ,  &
                (b( 9),nlpp  ) ,(b(39),nbpc  ) ,(b(91),lpch ) ,  &
                (b(92),ldict ) ,(b(95),sperlk) ,(machx,mach ) ,  &
                (m1(1),mconst(1))     ,(m2(1)  ,mconst(111) )
 
    !  DEFINE SYSTEM RELEASE DATE
 
    DATA imnth, iyr1, iyr2 /4HAPR., 3H 19, 2H95 /
    DATA xx   , yy    / 1.2E-38, 0.3E-38       /
    DATA mvax , abcd  , ka   /  1H1,   4HABCD  ,4HA           /
 
    !     MACH   = MACHX = HOST MACHINE
    !              ANY SUBROUTINE, THAT USES 'MACHX' INSTEAD OF 'MACH' IN
    !              LABEL COMMON /MACHIN/, CONTAINES MACHINE CONSTANTS THAT
    !              ARE USED LOCALLY.
    !     NMACH  = NUMBER OF MACHINES
    !     MCONST = ARRAY CONTAINING MACHINE DEPENDENT CONSTANTS
 
 
    !     COSMIC/NASTRAN SUPPORTS ONLY MACHINES 2, 5, 6, 7, 8, 9, 10, 16,
    !     21, & 22.   CONSTANTS FOR OTHER MACHINES MAY NOT BE EXACT
 
    !     -MACHINE-    IBM/  UNIVAC   CDC   DEC/    DEC/   SUN   IBM/    HP
    !           DUMMY   MVS    FTN   FTN5    VMS  ULTRIX SOLARIS  AIX    UX
    !     MACH = -1-  ---2-  ---3-  ---4-  ---5-  ---6-  ---7-  ---8-  ---9-
 
    !           SGI    MAC   CRAY  CONVEX   NEC  FUJITSU  SUN  AMDAHL  PRIME
    !           IRIS        UNICOS                       SUNOS
    !          --10-  --11-  --12-  --13-  --14-  --15-  --16-  --17-  --18-
 
    !          INTEL          DEC/   DEC/
    !          LINUX  DUMMY OPENVMS  OSF
    !          --19-  --20-  --21-  --22-
 
    !  MACHINE NAMES
    DATA comput/ 'DUMMY      ', 'IBM        ', 'UNIVAC     ',  &
                 'CDC        ', 'DEC-VAX    ', 'DEC-MIPS   ',  &
                 'SUN        ', 'IBM RS6000 ', 'HP         ',  &
                 'SGI        ', 'MACINTOCH  ', 'CRAY       ',  &
                 'CONVEX     ', 'NEC        ', 'FUJITSU    ',  &
                 'SUN        ', 'AMDAHL     ', 'PRIME      ',  &
                 'INTEL x86  ', 'DUMMY      ', 'DEC-ALPHA  ',  &
                 'DEC-ALPHA  '/
 
    !  MACHINE OPERATING SYSTEM
    DATA compos/ '       ', 'MVS    ', 'FTN    ',  &
                 'FTN5   ', 'VMS    ', 'ULTRIX ',  &
                 'SOLARIS', 'AIX    ', 'HP-UX  ',  &
                 'IRIX   ', '       ', 'UNICOS ',  &
                 '       ', '       ', '       ',  &
                 'SUNOS  ', '       ', '       ',  &
                 'LINUX  ', '       ', 'OPENVMS',  &
                 'OSF    '/
 
    DATA    nmach / 22 /,    m1/ &
        !     SYSBUF  =   LENGTH OF NASTRAN I/O BUFFER
        200,  4100,   871,  1042,  1028,  1028,  1028,  1028,  1028,  &
        1028,  1028,  2052,  1028,  2052,  2052,  1028,  1028,  1028,  &
        1028,  1028,  1028,  1028, &
        !     INTP(X100)  =  FORTRAN UNIT NO. FOR INPUT DATA
        !     OUTTAP      =  FORTRAN UNIT NO. FOR PRINTED OUTPUT
        506,  506,  506,  506,  506,  506,  506,  506,  506,  &
        506,  506,  506,  506,  506,  506,  506,  506,  506,  &
        506,  506,  506,  506, &
        !     NLPP(X100)  =  NUMBER OF LINES PRINTED PER PAGE
        !     NWPIC       =  NUMBER OF WORDS PER INPUT CARD, USED ONLY IN XGPIBS
        5000, 5518, 5518, 4208, 5518, 5518, 5518, 5518, 5518,  &
        5518, 5518, 5509, 5518, 5500, 5500, 5500, 5500, 5500,  &
        5518, 5500, 5518, 5500, &
        !     NBPC(X100)  =  NUMBER OF BITS PER CHARACTER
        !     NBPW        =  NUMBER OF BITS PER WORD
        636,  832,  936,  660,  832,  832,  832,  832,  832,  &
        832,  832,  864,  832,  864,  864,  832,  832,  832,  &
        832,  832,  832,  832, &
        !     IPREC(X100) =  PRECISION (1 = S.P., 2 = D.P.)
        !     RECL(X10)   =  DIRECT FILE RECORD LENGTH (USED IN FORTRAN OPEN
        !                    STATEMENT) BY WORDS (= 1), OR BYTE (= NCPW)
        !     QP          =  REAL*16 PRECISION FLAG (1 = YES, 0 = NO)
        !WKBR3    2 0 0, 2 4 0, 2 1 1, 1 1 0, 2 1 1, 2 1 0, 2 4 0, 2 4 1, 2 4 1,  &
        200, 240, 210, 110, 210, 210, 240, 240, 240,  &
        210, 200, 180, 240, 100, 100, 200, 200, 200,  &
        240, 200, 210, 200 /
 
    DATA          m2/ &
        !     LPCH(X100)  =  FORTRAN UNIT NO. FOR PUNCHED OUTPUT
        !     LDICT       =  FORTRAN UNIT NO. FOR RESTART DICTIONARY PUNCH
        703,  707,  103,  707,  104,  104,  104,  104,  104,  &
        104,  104,  104,  104,  104,  104,  104,  104,  104,  &
        104,  104,  104,  104, &
        !     LOWPW, HIGHPW = MACHINE NUMERIC RANGE FOR S. P. REAL NUMBER,
        !     USED ONLY BY RCARD, RCARD2, XRCARD AND YRCARD
        38,    75,    38,   321,    38,    38,    38,    38,    38,  &
        38,    38,  2465,    38,     0,     0,     0,     0,     0,  &
        38,     0,    38,     0, &
        !     NUDFLW(X100) =  FLOATING NUMBER UNDERFLOW CONTROL
        !                     (USED ONLY BY FQRW AND FQRWV)
        !     MXFL         =  MAXINUM FILES MAXFIL CAN REQUEST VIA THE NASTRAN
        !                     CARD, USED ONLY IN NASCAR
        1650, 1650, 1849, 1475,  875, 1675, 1675, 1675, 1675,  &
        1675, 1675, 1675, 1675, 1675, 1675, 1675, 1675, 1675,  &
        1675, 1675,  975, 1675, &
        !     KSHIFT  =  SHIFT COUNTS USED IN A DIVIDE TO CONVERT A GINO LOC
        !                RETURNED FROM SAVPOS TO GINO BLOCK NUMBER, USED IN EMA
        1,  4096,   4096,262144,  4096,  4096, 4096,  4096,  4096,  &
        4096,  4096,   4096,  4096,     0,     0,    0,     0,     0,  &
        4096,     0,   4096,     0, &
        !     MANTISSA BITS, USED ONLY IN SDCMPS
        000, 2426, 2760, 4896, 2355, 2355, 2352, 2355, 2355,  &
        2355, 2355, 4896, 2352, 4896, 4896, 2355, 2355, 2355,  &
        2355,  000, 2355,  000/
 
    !     DEFINE SYSTEM (42), SYSTEM(43), SYSTEM(44)
    idate(1) = imnth
    idate(2) = iyr1
    idate(3) = iyr2
 
    !     MACHINE TYPE IS SET HERE
    !     +++++++++++++++++++++++++++++++
100 mach = 19
    mchnam = comput(mach)
    machos = compos(mach)
    sysbuf = mconst(mach)
    ibmcdc = 1
    IF (mach == 2 .OR. mach == 4) ibmcdc = 0
 
    i  = mach + nmach
    intp   = mconst(i)/100
    outtap = MOD(mconst(i),100)
 
    i  = i + nmach
    nlpp   = mconst(i)/100
    nwpic  = MOD(mconst(i),100)
 
    i  = i + nmach
    nbpc   = mconst(i)/100
    nbpw   = MOD(mconst(i),100)
 
    i  = i + nmach
    iprec  = mconst(i)/100
    RECL   = MOD(mconst(i),100)/10
    qp     = MOD(mconst(i),10)
 
 
    i  = i + nmach
    lpch   = mconst(i)/100
    ldict  = MOD(mconst(i),100)
 
    !     MACHINE S.P. RANGE
    i  = i + nmach
    highpw = mconst(i)
    lowpw  = 1 - highpw
    IF (mach == 2) lowpw = -78
    IF (mach == 4) lowpw = -292
 
    !     FLOATING NUMBER UNDERFLOW CONTROL
    !     MAXINUM FILES FOR MAXFIL CHECK
 
    i  = i + nmach
    nudflw = mconst(i)/100
    mxfl   = MOD(mconst(i),100)
 
    !     SHIFT COUNTER FOR EMA SUBROUTINE
    i  = i + nmach
    kshift = mconst(i)
 
    !     MANTISSA BITS
    i  = i + nmach
    mtisa  = mconst(i)/100
    IF (iprec == 2) mtisa = MOD(mconst(i),100)
 
    !     NUMBER OF BITS PER HALF WORD, USED MAINLY FOR INTEGER PACKING
    !     IHALF = NBPW/2
    !     JHALF = 2**IHALF - 1
    ihalf = 16
    jhalf = 65535
 
    !     NUMBER OF CHARACTERS PER WORD
    ncpw = nbpw/nbpc
 
    !     ZERO FIELD KA, AK AND GENERATE A MASK FOR FIRST BYTE
    ak   = khrfn1(0,1,ka,4)
    ka   = khrfn1(0,1,ka,1)
    i    = 2**nbpc - 1
    mask = lshift(i,nbpw-nbpc)
 
    !     CHECK BCD WORD (NOT CHARACTER WORD) STORING ORDER.
    !     IF 'ABCD' IS STORED INTERNALLY IN A-B-C-D ORDER, SET ORDER TO 0,
    !     OR  IF IT IS STORED IN REVERSED ORDER, D-C-B-A,  SET ORDER TO 1
     i = andf(abcd,mask)
    order = 0
    IF (nbpw < 60 .AND. i /= ka .AND. i /= ak) order = 1
 
    !     CHECK SYSTEM LOC OR %LOC FUNCTION.
    !     IF SYSTEM LOC FUNCTION IS WORD COUNT, SET LOCF TO 1
    !     IF SYSTEM LOC FUNCTION IS BYTE COUNT, SET LOCF TO NCPW
    lqro = 1000
    i    = locfx(b(11)) - locfx(b(1))
    locf = i/10
 
    !     MERGE LOCF, QP, RECL, AND ORDER INTO LQRO
    lqro = locf*1000 + qp*100 + RECL*10 + order
 
    !     GENERATE MASKS
    mask2  = complf(lshift(2**nbpc-1, nbpw-4*nbpc))
    mask3  = rshift(complf(0),1)
    mzero  = lshift(1,nbpw-1)
    two(1) = lshift(1,31)
 
    CALL cnstdd
    linkno = lnknos(1)
 
    RETURN
END SUBROUTINE btstrp
