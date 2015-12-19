SUBROUTINE outpt5
     
!     DRIVER OF OUTPUT5 MODULE
!     COPIES UP TO 5 GINO DATA BLOCKS TO TAPE, BY FORTRAN WRITE,
!     FORMATTED (ASCII), OR UNFORMATTED (BINARY)
 
!     THIS MODULE HAS BEEN EXPANDED TO INCLUDE TABLE DATA BLOCKS.
!     ORIGINALLY IT HANDLES ONLY MATRIX DATA BLOCKS.  G.CHAN/MAY 88
 
!     ==== TABLE  ====
!     OUTPT5 CALLS TABLE5 TO PROCESS TABLE DATA BLOCKS
!      . UNFORMATTED (BINARY) OR FORMATTED (UNDER P4 CONTROL)
!      . IF BINARY, EACH RECORD IS WRITTEN OUT BY -
!           WRITE (OUT) L,(Z(J),J=1,L)
!      . IF FORMATTED,  5 BYTES ARE USED FOR BCD WORD,
!                      10 BYTES FOR INTEGER,
!                      15 BYTES FOR REAL, S.P. OR D.P.
!      . A HEADER RECORD, WHICH CONFORMS TO OUTPT5 HEADER STANDARD,
!           IS WRITTEN OUT FIRST, PRECEEDING THE TABLE DATA RECORDS.
 
!     ==== MATRIX ====
!     COPY GINO MATRIX DATA BLOCK(S) ONTO FORTRAN UNIT IN
!      . UNPACKED BANDED RECORD
!      . BANDED COLUMN  RECORD (FIRST TO LAST NON-ZERO ELEMENTS),
!      . UNFORMATTED (BINARY) OR FORMATTED
!      . SINGLE PRECISION OR DOUBLE, REAL OR COMPLEX DATA
!      . OUTPUT FORTRAN TAPE INPI (I=T,1,2,..,9) FOR UNIVAC, IBM, VAX
!                    OR TAPE UTI  (I=  1,2,..,5) FOR CDC
!        (DEFAULT=INP1, UNIT 15, OR UT1, UNIT 11)
 
!     THIS MODULE HANDLES ONLY MATRIX DATA BLOCKS, NOT TRUE ANY MORE
 
!     UNFORMATTED RECORDS CAN ONLY BE USED BY THE SAME COMPUTER SYSTEM,
!     WHILE FORMATTED RECORDS CAN BE USED ACROSS COMPUTER BOUNDARY
!     (E.G. WRITTEN BY CDC MACHINE AND READ BY IBM) AND ALSO, CAN BE
!     EDITED BY SYSTEM EDITOR, OR PRINTED OUT BY SYSTEM PRINT COMMAND.
 
!     CALL TO THIS MODULE IS
 
!     OUTPUT5  IN1,IN2,IN3,IN4,IN5//C,N,P1/C,N,P2/C,N,P3/C,N,P4
!                                  /C,N,T1/C,N,T2/C,N,T3/C,N,T4... $
 
!              P1=+N, SKIP FORWARD N MATRIX DATA BLOCKS OR TABLES BEFORE
!                     WRITE. (EXCEPT THE FIRST HEADER RECORD. EACH
!                     MATRIX DATA BLOCK OR TABLE, PRECEEDED BY A HEADER
!                     RECORD, IS A COMPLETE MATRIX OR TABLE, MADE UP OF
!                     MANY PHYSICAL RECORDS.
!                     SKIP TO THE END OF TAPE IF P1 EXCEEDS THE
!                     NO. OF DATA BLOCKS AVAILABLE ON THE OUTPUT FILE)
!              P1= 0, NO ACTION TAKEN BEFORE WRITE. (DEFAULT)
!              P1=-1, FORTRAN TAPE IS REWOUND, A TAPE HEADER RECORD IS
!                     WRITTEN TO TAPE. DATA IN FIRST GINO DATA BLOCK IS
!                     COPIED TO TAPE, FOLLOWED BY 4 MORE GINO DATA
!                     BLOCKS IF THEY ARE PRESENT.
!                     AT END, NO EOF WRITTEN, AND TAPE NOT REWOUND
!              P1=-3, THE NAMES OF ALL DATA BLOCKS ON FORTRAN TAPE
!                     ARE PRINTED AND WRITE OCCURS AT THE END OF TAPE
!              P1=-9, WRITE AN INTERNAL END-OF-FILE RECORD, FOLLOWED BY
!                     A SYSTEM ENDFILE MARK, AND REWIND FORTRAN TAPE
!              P2  IS THE FORTRAN UNIT NO. ON WHICH THE DATA BLOCKS WILL
!                     BE WRITTEN.  DEFAULT IS 15 (INP1 FOR UNIVAC, IBM,
!                     VAX), OR UNIT 11 (UT1 FOR CDC)
 
!              P3  IS TAPE ID IF GIVEN BY USER. DEFAULT IS XXXXXXXX
 
!              P4= 0, OUTPUT FILE IS FORTRAN WRITTEN, UNFORMATTED
!              P4= 1, OUTPUT FILE IS FORTRAN WRITTEN, FORMATTED
!                     (BCD IN 2A4, INTEGER IN I8, REAL IN 10E13.6 AND
!                      D.P. IN 5D26.17)
!              P4= 2, SAME AS P4=1, EXECPT 5E26.17 IS USED FOR S.P. REAL
!                     DATA. P4=2 IS USED ONLY IN MACHINES WITH LONG WORD
!                     FOR ACCURACY (60 OR MORE BITS PER WORD)
 
!              TI     10 WORD ARRAY USED ONLY BY TABLE BLOCK DATA.
!                     TO OVERRIDE AUTOMATIC FORMAT TYPE SETTING.
 
!     OUTPT5 LOGIC -
!                                                       (P4=0)   (P4=1)
!     RECORD  WORD        CONTENTS                      BINARY   FORMAT
!     ------  ----  --------------------------------   -------  -------
!        0            TAPE HEADER RECORD -
!              1,2    TAPEID                             2*BCD      2A4
!              3,4    MACHINE (2ND WORD BLANK)           2*BCD      2A4
!              5-7    DATE                               3*INT      3I8
!               8     SYSTEM BUFFSIZE                      INT       I8
!               9     P4 (0,1, OR 2)                       INT       I8
!      1A,1B%         FIRST MATRIX HEADER RECORD -
!               1     ZERO                                 INT       I8
!              2,3    ONE,ONE                            2*INT      2I8
!               4     D.P. ZERO                           F.P.   D26.17
!              5-10   MATRIX TRAILER                     6*INT      6I8
!                     (COL,ROW,FORM,TYPE,MAX,DENSITY)
!             11,12   DMAP NAME OF FIRST INPUT MATRIX    2*BCD      2A4
!      2A,2B    1     1 (FIRST MATRIX COLUMN ID)           INT       I8
!               2     COLUMN LOC. OF FIRST NON-ZERO ELEM.  INT       I8
!               3     COLUMN LOC. OF LAST  NON-ZERO ELEM.  INT       I8
!              1-W    FIRST BANDED COLUMN DATA            F.P.     (**)
!                     (W=WORD3-WORD2)
!      3A,3B    1     2 (SECOND MATRIX COLUMN ID)          INT       I8
!              2-3    FIRST AND LAST NON-ZERO ELEM LOC.  2*INT      2I8
!              1-W    SECOND BANDED COLUMN DATA           F.P.     (**)
!      4A,4B   1-3    THIRD  MATRIX COLUMN, SAME FORMAT  3*INT      3I8
!              1-W    AS RECORD 1                         F.P.     (**)
!        :      :       :
!      ZA,ZB    1     (A NULL COLUMN ID)                   INT       I8
!              2,3    1,1                                2*INT      2I8
!               1     0.0                                 F.P.     (**)
!        :      :       :
!      MA,MB   1-3    LAST MATRIX COLUMN, SAME AS REC #2 3*INT      3I8
!              1-W    LAST BANDED COLUMN DATA             F.P.     (**)
 
!      SA,SB    :     SECOND MATRIX HEADER RECORD   3*INT+F.P. 3I8+D26.
!                                                       +2*BCD   +2*BCD
!                                                       +6*INT     +6I8
!    S+1A,S+1B 1-W    FIRST THRU LAST COLS OF 2ND MATRIX
!        :      :     REPEAT FOR MORE MATRICES
!        :      :     (UP TO 5 MATRIX DATA BLOCKS PER ONE OUTPUT FILE)
 
!    EOFA,EOFB  1     -1                                   INT       I8
!              2,3    1,1                                2*INT      2I8
!               1     D.P. ZERO                           F.P.   D26.17
 
!                                                               - NOTE -
!                                                  BCD AND INTEGERS IN 8
!                                         SINGLE PRECISION REAL IN  13.6
!                                         DOUBLE PRECISION DATA IN 26.17
!                                         S.P. LOGN WORD MACHINE   26.17
 
!     WHERE   %  RECORDS A AND B ARE 2 (OR MORE) RECORDS ON FORMATTED
!                OUTPUT FILE, WHILE
!                A & B ARE 1 CONTINUOUS RECORD IN UNFORMATTED TAPE
!           (**) IS (10E13.6) FOR S.P.REAL, OR (5D26.17) FOR D.P. DATA.
!                OR (5E26.17) FOR S.P. AND D.P. DATA (P4=2 ONLY)
!     NOTE -
!     NO SYSTEM END-OF-FILE MARK WRITTEN BETWEEN MATRICES.
 
!     TO READ BINARY TAPE              TO READ FORMATTED TAPE
!     ----------------------------     --------------------------------
!                LOGICAL SP,DP
!                INTEGER COL,ROW,FORM,TYPE,DENS,FILE(2),IZ(M,N)
!    *                   TAPEID(2),MAC(2),DATE(3),BUFSZ,P4
!                DOUBLE PRECISION DZ(M/2,N/2),DTEMP
!                COMMON  /ZZZZZZ/ Z(M,N)
!                EQUIVALENCE      (Z,IZ,DZ)
!                DATA     SP,DP / .TRUE.,.FALSE./
!     READ (TAPE,ERR=7)                READ (TAPE,10,ERR=7)
!    *           TAPEID,MAC,DATE,BUFSZ,P4
!   1            K = 0
!   2            K = K + 1
!     READ (TAPE,ERR=7,END=3) I,JB,    IF (SP) READ (TAPE,8,ERR=7,END=3)
!    *                        JE               I,JB,JE,( Z(J,K),J=JB,JE)
!                                      IF (DP) READ (TAPE,9,ERR=7,END=3)
!    *                                         I,JB,JE,(DZ(J,K),J=JB,JE)
!                IF (I)   3,            4,     6
!C                      EOF,MATRIX-HEADER,COLUMN
!   3            CONTINUE
!C               (EOF ENCOUNTERED, COMPLETE TAPE READ)
!                CALL EXIT
!   4            BACKSPACE TAPE
!                                      BACKSPACE TAPE
!C               (MATRIX-HEADER READ)
!     READ (TAPE) J,J,J,               READ (TAPE,11) J,J,J
!    *           DTEMP,COL,ROW,FORM,TYPE,MAX,DENS,FILE
!                DP = .FALSE.
!                IF (TYPE.EQ.2 .OR. TYPE.EQ.4) DP=.TRUE.
!                SP = .NOT.DP
!                JTYP = TYPE
!                IF (TYPE .EQ. 3) JTYP = 2
!                IF (COL*JTYP.GT.M .OR. ROW*JTYP.GT.N) STOP 'Z DIM ERR'
!                J = COL*ROW*JTYP
!                DO 5 I = 1,J
!   5            Z(I,1) = 0.0
!                GO TO 1
!   6            CONTINUE
!C               (A COLUMN OF MATRIX READ)
!                IF (I .NE. K) STOP 'COLUMN COUNTER MISSMATCH'
!                GO TO 2
!   7            STOP 'READ ERROR. CHECK TAPE FORMAT TYPE'
!                                     8 FORMAT (3I8,/,(10E13.6))
!                                     9 FORMAT (3I8,/,(5D26.17))
!                                    10 FORMAT (4A4,5I8)
!                                    11 FORMAT (3I8,/,D26.17,6I8,2A4)
!C             FOR LONG WORD MACHINE  8 FORMAT (3I8,/,(5E26.17))
 
!     SEE SUBROUTINE INPTT5 FOR MORE COMPREHENSIVE DETAILS IN RECOVERING
!     MATRIX DATA FROM THE TAPE GENERATED IN THIS OUTPT5 ROUTINE.
!     OR SUBROUTINE TABLE-V FOR TABLE DATA BLOCK RECOVERY.
 
!     WRITTEN BY G.CHAN/UNISYS   1987
 
 IMPLICIT INTEGER (a-z)
 LOGICAL :: p40,p40s,p40d,p41,p41s,p41d,p41c,complx
 INTEGER :: trl(9),NAME(2),tapeid(2),subnam(2),ti,dt(3), fn(3,10),NONE(2)
 REAL :: rz(1),x,zero
 DOUBLE PRECISION :: dz(1),dx,dzero
 CHARACTER (LEN=8) :: binary,formtd,bf
 CHARACTER (LEN=25) :: sfm
 CHARACTER (LEN=29) :: uim
 CHARACTER (LEN=25) :: uwm
 CHARACTER (LEN=23) :: ufm
!WKBNB
 CHARACTER (LEN=80) :: dsnames
 COMMON /dsname/  dsnames(80)
!WKBNE
 COMMON /xmssg /  ufm,uwm,uim,sfm
 COMMON /BLANK /  p1,p2,p3(2),p4,ti(10)
 COMMON /machin/  mach,ijhalf(3),mchnam
 COMMON /system/  ibuf,nout,dum6(7),line,dumm4(4),date(3),  &
     dum22(22),nbpw,dum50(50),lpch
 COMMON /zzzzzz/  iz(1)
 COMMON /unpakx/  ityp,ii,jj,incr
 EQUIVALENCE      (rz(1),dz(1),iz(1))
 DATA    binary,           formtd,         subnam             /  &
     'BINARY  ',        'FORMATTD',     4HOUTP, 2HT5       /
 DATA    zero,    dzero,   izero,  one,    mone,   fn         /  &
     0.0,     0.0D0,   0,      1,      -1,     30*4H      /
 DATA    mtrx,    tble,    BLANK / 4HMTRX, 4HTBLE, 4H         /
 DATA    NONE  /  4H (no,  4HNE) /
 
!     IF MACHINE IS CDC OR UNIVAC, CALL CDCOPN OR UNVOPN TO OPEN OUTPUT
!     FILE, A FORMATTED SEQUENTIAL TAPE.  NO CONTROL WORDS ARE TO BE
!     ADDED TO EACH FORMATTED RECORD. RECORD LENGTH IS 132 CHARACTERS,
!     AN ANSI STANDARD.
 
!WKBD IF (MACH .EQ. 3) CALL UNVOPN (P2)
!WKBD IF (MACH .EQ. 4) CALL CDCOPN (P2)
 bf = binary
 IF (p4 >= 1) bf = formtd
 CALL page
 CALL page2 (1)
 WRITE  (nout,3) uim,bf,p1
 3    FORMAT (a29,', MODULE OUTPUT5 CALLED BY USER DMAP ALTER, ON ',a8,  &
     ' TAPE,', /5X,'WITH FOLLOWING REQUEST  (P1=',i2,1H))
 IF (p1 == -9) WRITE (nout,4)
 IF (p1 == -3) WRITE (nout,5)
 IF (p1 == -1) WRITE (nout,6)
 IF (p1 ==  0) WRITE (nout,7)
 IF (p1 >  0) WRITE (nout,8) p1
 4    FORMAT (5X,'WRITE AN INTERNAL E-O-F RECORD, FOLLOWED BY A SYSTEM',  &
     ' E-O-F MARK, AND REWIND OUTPUT TAPE')
 5    FORMAT (5X,'REWIND TAPE, PRINT DATA BLOCK NAMES AND THEN WRITE ',  &
     'AFTER THE LAST DATA BLOCK ON TAPE')
 6    FORMAT (5X,'REWIND, WRITE A TAPE HEADER RECORD, THEN FOLLOWED BY '  &
     ,       'DATA BLOCKS WRITING.',/5X,'AT END, NO EOF AND NO REWIND')
 7    FORMAT (5X,'DATA BLOCKS ARE WRITTEN STARTING AT CURRENT TAPE ',  &
     'POSITION. AT END, NO EOF AND NO REWIND')
 8    FORMAT (5X,'SKIP FORWARD',i4,' DATA BLOCKS BEFORE WRITING (TAPE ',  &
     'HEADER RECORD NOT COUNTED AS A DATA BLOCK).', /5X,  &
     'NO REWIND BEFORE SKIPPING. AT END, NO EOF AND NO REWIND')
 
 buf1 = korsz(rz(1)) - ibuf - 1
 IF (buf1 <= 0) CALL mesage (-8,0,subnam)
 out = p2
 wrt = 0
 lfn = -1
 IF (p1 == -3) lfn = 0
 
!     SET P4 FLAGS
 
!     SET P40  TO .TRUE. IF USER SPECIFIES P4 TO ZERO (BINARY)
!     SET P41  TO .TRUE. IF USER SPECIFIES P4 TO ONE  (FORMATTED)
!     SET P40D TO .TRUE. IF P40 IS TRUE AND DATA IS IN D.P.
!     SET P40S TO .TRUE. IF P40 IS TRUE AND DATA IS IN S.P.
!     SET P41D TO .TRUE. IF P41 IS TRUE AND DATA IS IN D.P.
!     SET P41S TO .TRUE. IF P41 IS TRUE AND DATA IS IN S.P.
!     SET P41C TO .TRUE. IF P4=2, AND RESET P41S AND P41D TO .FALSE.
 
 p40d = .false.
 p41s = .false.
 p41d = .false.
 p41c =  p4 == 2 .AND. nbpw >= 60
 p41  = .false.
 IF (p4 >= 1) p41 = .true.
!WKBNB
 CLOSE ( UNIT=out )
 IF ( p4 /= 0  ) GO TO 1
 OPEN ( UNIT=out, FILE=dsnames(out), FORM='UNFORMATTED' ,STATUS='UNKNOWN' )
 GO TO 2
 1     CONTINUE
 OPEN ( UNIT=out, FILE=dsnames(out), STATUS='UNKNOWN' )
 2     CONTINUE
!WKBNE
 IF (p41c) GO TO 10
 p41s =  p41
 IF (p41 ) p41d = .NOT.p41s
 10   p40  = .NOT.p41
 p40s =  p40
 IF (p40) p40d = .NOT.p40s
 IF (p1 /= -9) GO TO 20
 
!     FINAL CALL TO OUTPUT5
 
 IF (p40) WRITE (out    ) mone,one,one,dzero
 IF (p41) WRITE (out,290) mone,one,one,dzero
 ENDFILE out
 REWIND out
 RETURN
 
 20   IF (p1 == -3) GO TO 60
 IF (p1 == -1) GO TO 180
 IF (p1 < 0.0) THEN
   GO TO    30
 ELSE IF (p1 == 0.0) THEN
   GO TO   190
 ELSE
   GO TO    65
 END IF
 
 30   WRITE  (nout,35) ufm,p1
 35   FORMAT (a23,' 4120, MODULE OUTPUT5 - ILLEGAL VALUE FOR FIRST ',  &
     'PARAMETER = ',i8)
 40   ERR = -37
 50   CALL mesage (ERR,INPUT,subnam)
 RETURN
 
!     OLD TAPE. CHECK TAPE ID
 
 60   REWIND out
 65   IF (p40) READ (out,    END=150) tapeid,NAME,dt,i,k
 IF (p41) READ (out,185,END=150) tapeid,NAME,dt,i,k
 IF (tapeid(1) == p3(1) .AND. tapeid(2) == p3(2)) GO TO 70
 WRITE  (nout,67) tapeid,p3
 67   FORMAT ('0*** WRONG TAPE MOUNTED - TAPEID =',2A4,', NOT ',2A4)
 GO TO 40
 70   CALL page2 (6)
 WRITE  (nout,75) tapeid,NAME,dt,i
 75   FORMAT (/5X,'MODULE OUTPUT5 IS PROCESSING TAPE ',2A4, /5X,  &
     'WRITTEN BY ',2A4, /5X,'ON ',i2,1H/,i2,1H/,i2,  /5X, 'BUFFSIZE USED =',i7,/)
 IF (k == 0) WRITE (nout,80) binary
 IF (k >= 1) WRITE (nout,80) formtd
 80   FORMAT (5X,'ORIGINAL TAPE IS ',a8)
 IF (k == p4) GO TO 90
 WRITE  (nout,85) ufm,p4
 85   FORMAT (a23,', THE 4TH PARAMETER TO OUTPUT5 DOES NOT AGREE WITH ',  &
     'ORIG. TAPE FORMAT    P4=',i5,/)
 CALL mesage (-37,0,subnam)
 
!     TO SKIP P1 MATRIX DATA BLOCKS OR TABLES ON THE OLD OUTPUT FILE
!     OR TO TABULATE TAPE CONTENTS IF P1 = -3
 
 90   lfn = 0
 100  IF (p40 ) READ (out,    ERR=160,END=150) nc,jb,je
 IF (p41s) READ (out,280,ERR=100,END=150) nc,jb,je,( x,j=jb,je)
 IF (p41c) READ (out,285,ERR=100,END=150) nc,jb,je,( x,j=jb,je)
 IF (p41d) READ (out,290,ERR=100,END=150) nc,jb,je,(dx,j=jb,je)
 IF (nc < 0) THEN
   GO TO   140
 ELSE IF (nc == 0) THEN
   GO TO   120
 ELSE
   GO TO   100
 END IF
 110  IF (p40 ) READ (out,    ERR=160,END=150) l
 IF (p41 ) READ (out,115,ERR=100,END=150) l,(table,j=1,l)
 115  FORMAT (i10,24A5,/,(26A5))
 IF (l < 0) THEN
   GO TO   140
 ELSE IF (l == 0) THEN
   GO TO   120
 ELSE
   GO TO   110
 END IF
 120  IF (p1 /= -3 .AND. lfn >= p1) GO TO 140
 lfn = lfn + 1
 BACKSPACE out
 IF (p41) BACKSPACE out
 IF (p40) READ (out    ) i,i,i,dx,j,j,j,j,k,k,fn(1,lfn),fn(2,lfn)
 IF (p41) READ (out,250) i,i,i,dx,j,j,j,j,k,k,fn(1,lfn),fn(2,lfn)
 IF (k > 0 .AND. j >= 1 .AND. j <= 4) GO TO 130
 fn(3,lfn) = tble
 GO TO 110
 130  fn(3,lfn) = mtrx
 IF (p40) GO TO 100
 p41s = .false.
 p41d = .false.
 p41c =  p4 == 2 .AND. nbpw >= 60
 IF (p41c) GO TO 100
 IF (j == 1 .OR. j == 3) p41s = .true.
 p41d = .NOT.p41s
 GO TO 100
 140  IF (p41) BACKSPACE out
 150  BACKSPACE out
 IF (p1 == -3 .AND. lfn > 0) GO TO 430
 GO TO 200
 
 160  WRITE  (nout,170) uwm,tapeid
 170  FORMAT (a25,' FROM OUTPUT5 MODULE. ERROR WHILE READING ',2A4)
 GO TO 40
 
!     NEW TAPE (P1=-1)
 
!     WRITE A TAPE IDENTIFICATION RECORD (NOTE -THIS IS THE ONLY TIME
!     A TAPE HEADER RECORD IS WRITTEN)
 
 180  IF (p1 /= -1) GO TO 200
 REWIND out
 trl(1) = p3(1)
 trl(2) = p3(2)
 trl(3) = mchnam
 trl(4) = BLANK
 trl(5) = date(1)
 trl(6) = date(2)
 trl(7) = date(3)
 IF (p40) WRITE (out    ) (trl(j),j=1,7),ibuf,p4
 IF (p41) WRITE (out,185) (trl(j),j=1,7),ibuf,p4
 185  FORMAT (4A4,5I8)
 190  lfn = 0
 
!     COPY MATRICES OR TABLES OUT TO TAPE
 
 200  DO  mx = 1,5
   INPUT = mx + 100
   CALL fname (INPUT,NAME)
   IF (NAME(1) == NONE(1) .AND. NAME(2) == NONE(2)) GO TO 390
   trl(1) = INPUT
   CALL rdtrl (trl)
   IF (trl(1) <= 0) GO TO 390
   IF (trl(1) > 0) GO TO 220
   CALL page2 (3)
   WRITE  (nout,210) INPUT,NAME
   210  FORMAT (/5X,'INPUT FILE ',2A4,'(',i3,') IS PURGED. NO DATA ',  &
       'TRANSFERRED TO OUTPUT FILE')
   CYCLE
   220  IF (trl(4) > 8 .OR. trl(5) > 4 .OR. trl(6) <= 0 .OR. trl(7) <= 0  &
       ) CALL table5 (*400,INPUT,out,trl,buf1,wrt,lfn,fn)
   col  = trl(2)
   row  = trl(3)
   TYPE = trl(5)
   complx = .false.
   IF (TYPE >= 3) complx = .true.
   
!     CHECK FOR NULL MATRIX
   
   IF (row == 0 .OR. col == 0 .OR. TYPE == 0) GO TO 380
   
!     SET FLAGS FOR FORMATTED OR UNFORMATTED WRITE, SINGLE OR DOUBLE
!     PRECISION DATA, THEN WRITE THE MATRIX HEADER WITH PROPER FORMAT.
!     MATRIX HEADER CONSISTS OF ONE SCRATCH WORD, ORIGINAL MATRIX
!     TRAILER, AND MATRIX DMAP NAME
   
   p40s = .false.
   p40d = .false.
   p41s = .false.
   p41d = .false.
   p41c =  p4 == 2 .AND. nbpw >= 60
   IF (p41) GO TO 230
   IF (TYPE == 1 .OR. TYPE == 3) p40s = .true.
   p40d = .NOT.p40s
   GO TO 240
   230  IF (p41c) GO TO 240
   IF (TYPE == 1 .OR. TYPE == 3) p41s = .true.
   p41d = .NOT.p41s
   240  IF (p40) WRITE (out    ) izero,one,one,dzero,(trl(k),k=2,7),NAME
   IF (p41) WRITE (out,250) izero,one,one,dzero,(trl(k),k=2,7),NAME
   250  FORMAT (3I8,/,d26.17,6I8,2A4)
   wrt = 1
   
!     OPEN INPUT DATA BLOCK AND SAVE DMAP NAME IN FN ARRAY
   
   ERR = -1
   CALL OPEN (*50,INPUT,rz(buf1),0)
   CALL fwdrec (*50,INPUT)
   IF (lfn == -1 .OR. lfn >= 10) GO TO 260
   lfn = lfn + 1
   fn(1,lfn) = NAME(1)
   fn(2,lfn) = NAME(2)
   fn(3,lfn) = mtrx
   
!     UNPACK A MATRIX COLUMN, AND WRITE TO OUTPUT FILE THE BANDED DATA
!     (FROM FIRST TO LAST NON-ZERO ELEMENTS)
   
   260  ityp = TYPE
   incr = 1
   DO  nc = 1,col
     ii = 0
     jj = 0
     CALL unpack (*300,INPUT,rz)
     jb = ii
     je = jj
     nwds = jj - ii + 1
     IF (.NOT.complx) GO TO 270
     nwds = nwds + nwds
     je   = nwds + jb - 1
     270  IF (nwds > buf1) CALL mesage (-8,0,subnam)
     IF (p40s) WRITE (out) nc,jb,je,(rz(j),j=1,nwds)
     IF (p40d) WRITE (out) nc,jb,je,(dz(j),j=1,nwds)
     IF (p41s) WRITE (out,280,ERR=480) nc,jb,je,(rz(j),j=1,nwds)
     IF (p41c) WRITE (out,285,ERR=480) nc,jb,je,(rz(j),j=1,nwds)
     IF (p41d) WRITE (out,290,ERR=480) nc,jb,je,(dz(j),j=1,nwds)
     280  FORMAT (3I8,/,(10E13.6))
     285  FORMAT (3I8,/,(5E26.17))
     290  FORMAT (3I8,/,(5D26.17))
     CYCLE
     
!     A NULL COLUMN
     
     300  je = 1
     IF (complx) je = 2
     IF (p40s) WRITE (out    ) nc,one,je,( zero,i=1,je)
     IF (p40d) WRITE (out    ) nc,one,je,(dzero,i=1,je)
     IF (p41s) WRITE (out,280) nc,one,je,( zero,i=1,je)
     IF (p41c) WRITE (out,285) nc,one,je,( zero,i=1,je)
     IF (p41d) WRITE (out,290) nc,one,je,(dzero,i=1,je)
   END DO
   
!     CLOSE INPUT DATA BLOCK WITH REWIND.
   
   CALL CLOSE (INPUT,1)
   CALL page2 (10)
   WRITE  (nout,350) NAME,out,(trl(j),j=2,5),ibuf
   350  FORMAT (/5X,'MODULE OUTPUT5 UNPACKED MATRIX DATA BLOCK ',2A4,  &
       ' AND WROTE IT OUT TO', /5X,'FORTRAN UNIT',i4,  &
       ', IN BANDED DATA FORM (FIRST TO LAST NON-ZERO ELEMENTS)',  &
       /9X,'NO. OF COLS =',i8, /9X,'NO. OF ROWS =',i8, /16X,  &
       'FORM =',i8, /16X,'TYPE =',i8, /5X,'SYSTEM BUFFSIZE =',i8)
   IF (p40 ) WRITE (nout,360)
   IF (p41s) WRITE (nout,365)
   IF (p41c) WRITE (nout,370)
   IF (p41d) WRITE (nout,375)
   360  FORMAT (5X,'IN FORTRAN BINARY RECORDS')
   365  FORMAT (5X,'IN FORTRAN FORMATTED RECORDS - (3I8,/,(10E13.6))')
   370  FORMAT (5X,'IN FORTRAN FORMATTED RECORDS - (3I8,/,(5E26.17))')
   375  FORMAT (5X,'IN FORTRAN FORMATTED RECORDS - (3I8,/,(5D26.17))')
   CYCLE
   
!     NULL MATRIX, OR GINO DATA BLOCK IS NOT A MTRIX FILE
   
   380  CALL page2 (5)
   WRITE  (nout,385) uwm,NAME
   385  FORMAT (a25,' FROM OUTPUT5 MODULE. ',2A4,' IS EITHER A NULL ',  &
       'MATRIX OR NOT A MATRIX DATA BLOCK', /5X,  &
       'NO DATA WERE COPIED TO OUTPUT FILE',/)
   CYCLE
   
   390  trl(1) = INPUT + 1
   CALL rdtrl (trl)
   IF (trl(1) > 0) WRITE (nout,395) uwm,INPUT,NAME
   395  FORMAT (a25,' FROM OUTPUT5 MODULE. INPUT DATA BLOCK',i5,2H, ,2A4,  &
       ' IS EITHER PURGED OR DOES NOT EXIST')
   400  CONTINUE
 END DO
 
 IF (wrt == 0) WRITE (nout,410) uwm
 410  FORMAT (a25,' FROM OUTPUT5 MODULE. NO DATA BLOCK WRITTEN TO ',  &
     'OUTPUT FILE')
 ENDFILE out
 BACKSPACE out
 IF (p1 == -3) GO TO 460
 
!     PRINT LIST OF DATA BLOCKS ON FORTRAN TAPE (P1=-3).
 
 IF (lfn <= 0) RETURN
 430  CALL page2 (lfn+10)
 WRITE  (nout,440) out,mchnam,bf,(j,fn(1,j),fn(2,j),fn(3,j), j=1,lfn)
 440  FORMAT (/5X,'SUMMARY FROM OUTPUT5 MODULE', //16X,'DATA BLOCKS ',  &
     'WRITTEN TO FORTRAN UNIT',i4, /17X,'(BY ',a4,' MACHINE, ',  &
     a8,' RECORDS)', ///22X,'FILE',8X,'NAME',8X,'TYPE' /17X,  &
     9(4H----), /,(22X,i3,9X,2A4,4X,a4))
 IF (p1 == -3) GO TO 200
 IF (p40) GO TO 460
 CALL page2 (2)
 WRITE  (nout,450)
 450  FORMAT (/5X,'THIS FORMATTED OUTPUT FILE CAN BE VIEWED OR EDITED',  &
     ' VIA SYSTEM EDITOR',/)
 
 460  IF (mach == 3) CALL unvcls (p2)
 IF (mach == 4) CALL cdccls (p2)
 RETURN
 
!     WRITE ERROR
 
 480  WRITE  (nout,490) sfm
 490  FORMAT (a25,' IN WRITING OUTPUT FILE', /5X,'IBM USER - CHECK FILE'  &
     ,      ' ASSIGNMENT FOR DCB PARAMETER OF 132 BYTES')
 CALL mesage (-37,0,subnam)
END SUBROUTINE outpt5
