SUBROUTINE inptt5
     
!     DRIVER OF INPUTT5 MODULE
!     THIS MODULE HANDLES BOTH TABLE AND MATRIX DATA BLOCKS     5/88
 
!     THIS IS A COMPANION MODULE TO OUTPUT5
 
!     ==== TABLE  ====
!     CALLS TABLE-V ROUTINE TO COPY FROM A FORTRAN UNIT (FORMATTED OR
!     BINARY TAPE) TABLE DATA TO NASTRAN GINO TABLE DATA BLOCKS
 
!     ==== MATRIX ====
!     COPIES FROM A FORTRAN UNIT (BINARY OR FORMATTED TAPE) OF BANDED
!     MATRICES ONTO NASTRAN GINO MATRIX DATA BLOCKS, IN GINO PACKED
!     FORMAT
 
!     UNFORMATTED RECORDS CAN ONLY BE USED BY THE SAME COMPUTER SYSTEM,
!     WHILE FORMATTED RECORDS CAN BE USED ACROSS COMPUTER BOUNDARY
!     (E.G. WRITTEN BY CDC MACHINE AND READ BY IBM), AND ASLO, IT CAN
!     BE EDITED BY SYSTEM EDITOR, OR PRINTED OUT BY SYSTEM PRINT COMMAND
 
!     ******************************************************************
!     *                                                                *
!     *                       - IMPORTANT -                            *
!     *                                                                *
!     *  IF USER ASSEMBLES HIS OWN MATRIX IN INPUTT5 FORMAT, AND USES  *
!     *  INPUTT5 MODULE TO READ IT INTO NASTRAN, BE SURE THAT THE      *
!     *  DENSITY TERM (DENS) OF THE MATRIX TRAILER IS SET TO NON-ZERO  *
!     *  (NEED NOT BE EXACT) AND THE PRECISION TERM (TYPE) IS 1,2,3,   *
!     *  OR 4. OTHERWISE, HIS MATRIX WILL BE TREATED AS TABLE AND      *
!     *  EVERYTHING GOES HAYWIRE.                                      *
!     *                                                                *
!     ******************************************************************
 
!     CALL TO THIS MODULE IS
 
!     INPUTT5  /O1,O2,O3,O4,O5/C,N,P1/C,N,P2/C,N,P3/C,N,P4 $
 
!              P1=+N, SKIP FORWARD N MATRIX DATA BLOCKS OR TABLES BEFORE
!                     COPYING (EXCEPT THE FIRST HEADER RECORD. EACH
!                     MATRIX DATA BLOCK OR TABLE, PRECEEDED BY A HEADER
!                     RECORD, IS A COMPLETE MATRIX OR TABLE, MADE UP OF
!                     MANY PHYSICAL RECORDS.
!                     SKIP TO THE END OF TAPE IF P1 EXCEEDS THE NO. OF
!                     DATA BLOCKS AVAILABLE ON THE OUTPUT TAPE
!                     NO REWIND BEFORE SKIPPING)
!              P1= 0, NO ACTION TAKEN BEFORE COPYING. (DEFAULT)
!                     HOWEVER, IF TAPE IS POSITIONED AT THE BEGINNING,
!                     THE TAPE ID RECORD IS SKIPPED FIRST.
!              P1=-1, INPUT TAPE IS REWOUND, AND TAPEID CHECKED. IF
!                     OUTPUT GINO FILES ARE PRESENT, DATA FROM TAPE ARE
!                     THEN COPIED TO GINO FILES - IN PACKED MATRIX FORM
!                     IF MATRIX DATA, OR TABLE FORM IF TABLE DATA.
!              P1=-3, TAPE IS REWOUND AND READ. THE NAMES OF ALL DATA
!                     BLOCKS ON FORTRAN TAPE ARE PRINTED. AT END, TAPE
!                     IS REWOUND AND POSITIONED AFTER TAPE HEADER RECORD
!                     (NOTE - SERVICE UP TO 15 FILE NAMES ON ONE INPUT
!                     TAPE. AND THE 'AT END' TREATMENT IS NOT THE SAME
!                     AS IN OUTPUT5)
!              P1=-4  THRU -8 ARE NOT USED
!              P1=-9, REWIND TAPE
 
!              P2  IS THE FORTRAN UNIT NO. ON WHICH THE DATA BLOCKS WILL
!                     WRITTEN.  DEFAULT IS 16 (INP2 FOR UNIVAC,IBM,VAX),
!                     OR 12 (UT2 FOR CDC)
 
!              P3  IS TAPE ID IF GIVEN BY USER. DEFAULT IS XXXXXXXX
 
!              P4=0, OUTPUT TAPE IS FORTRAN WRITTEN, UNFORMATTED RECORDS
!              P4=1, OUTPUT TAPE IS FORTRAN WRITTED, FORMATTED
!                    BCD IN 2A4, INTEGER IN I8, S.P. REAL IN 10E13.6,
!                    AND D.P. IN 5D26.17.
!              P4=2, SAME AS P4=1, EXECPT FORMAT 5E26.17 IS USED FOR
!                    S.P. REAL DATA. (THIS OPTION IS USED ONLY IN
!                    MACHINES WITH 60 OR MORE BITS PER WORD)
 
 
!     CONTENTS OF INPUT TAPE, AS WRITTEN BY OUTPUT5
!                                                       (P4=0)   (P4=1)
!     RECORD  WORD        CONTENTS                      BINARY   FORMAT
!     ------  ----  --------------------------------   -------  -------
!        0            TAPE HEADER RECORD -
!              1,2    TAPEID                             2*BCD      2A4
!              3,4    MACHINE                            2*BCD      2A4
!              5-7    DATE                               3*INT      3I8
!               8     BUFFSIZE                             INT       I8
!               9     0 (BINARY), OR 1 OR 2 (FORMATTED)    INT       I8
!       1/2@          FIRST MATRIX HEADER RECORD -
!               1     ZERO                                 INT       I8
!              2,3    1,1                                2*INT      2I8
!               1     DUMMY (D.P.)                        F.P.   D26.17
!              2-7    MATRIX TRAILER                     6*INT      6I8
!                     (COL,ROW,FORM,TYPE,MAX,DENS)
!              8-9    MATRIX DMAP NAME                   2*BCD      2A4
!       3/4     1     1 (FIRST COLUMN ID)                  INT       I8
!               2     LOC. OF FIRST NON-ZERO ELEMENT, L1   INT       I8
!               3     LOC. OF LAST  NON-ZERO ELEMENT, L2   INT       I8
!              1-W    FIRST MATRIX COLUMN DATA            F.P.     (**)
!                     (W=L2-L1+1)
!       5/6     1     2 (SECOND COLUMN ID)                 INT       I8
!              2-3    LOC. OF FIRST AND LAST NON-ZERO    2*INT      2I8
!                     ELEMENTS
!              1-W    SECOND MATRIX COLUMN DATA           F.P.     (**)
!       7/8    1-3    THIRD MATRIX COLUMN, SAME FORMAT   3*INT      3I8
!              1-W    AS RECORD 1                         F.P.     (**)
!        :      :       :
!       M/M+1  1-3    LAST MATRIX COLUMN, SAME FORMAT    3*INT      3I8
!                     AS RECORD 1                         F.P.     (**)
!     M+2/M+3  1-3    SECOND MATRIX HEADER RECORD        3*INT      3I8
!               1     DUMMY                               F.P.     (**)
!              2-7    MATRIX TRAILER                     6*INT      6I8
!              8,9    MATRIX DMAP NAME                   2*BCD      2A4
!     M+4-N     :     FIRST THRU LAST COLUMNS OF MATRIX  3*INT      3I8
!                                                        +F.P.    +(**)
!        :      :     REPEAT FOR 3RD,4TH,5TH MATRICES
!        :      :     (UP TO 5 MATRIX DATA BLOCKS PER ONE OUTPUT TAPE)
 
!       EOF    1-3    -1,1,1                              3*INT     3I8
!               1     ZEROS (D.P.)                         F.P.  D26.17
 
!     @  RECORDS 1/2 (3/4, 5/6, ETC) ARE TWO RECORDS IN THE FORMATTED
!        TAPE, AND ARE PHYSICALLY ONE RECORD IN THE BINARY TAPE (AND
!        THE WORD COUNT SHOULD BE ADDED)
!     ** IS (10E13.6) FOR S.P.REAL OR (5D26.17) FOR D.P.DATA
!        (5E26.17) FOR LONG WORD MACHINE
 
!                                                               - NOTE -
!                                                  BCD AND INTEGERS IN 8
!                                                     S.P. REAL IN  13.7
!                                                     D.P. DATA IN 26.17
!                                             LONG WORD MACHINE IN 26.17
 
!     NO SYSTEM END-OF-FILE MARK WRITTEN BETWEEN MATRICES
!     EXCEPT FOR THE TAPE HEADER RECORD, AND THE MATRIX HEADERS, THE
!     ENTIRE FORMATTED INPUT TAPE CAN BE READ BY A STANDARD FORMAT
!     (3I8,/,(10E13.6)), (3I8,/,(5D26.17)), OR (3I8,/,(5E26.17))
 
!     ALSO, USER MAY OR MAY NOT CALL OUTPUT5 WITH P1=-9 TO WRITE AN
!     'OUPUT5 E-O-F' MARK ON TAPE. THIS CAUSED PROBLEM BEFORE.
 
!     THE PROCEDURE TO READ AND/OR WRITE THE TAPE IS COMMONLY USED
!     AMONG INPUTT5, OUTPUT5, AND DUMOD5. ANY PROCEDURE CHANGE SHOULD
!     BE MADE TO ALL THREE SUBROUTINES.
 
!     WRITTEN BY G.CHAN/UNISYS   1987
!     MAJOR REVISED 12/1992 BY G.C.
 
 IMPLICIT INTEGER (a-z)
 LOGICAL :: opn,p40,p40s,p40d,p41,p41s,p41d,p41c,debug
 INTEGER :: NAME(2),tapeid(2),mac(2),subnam(2),dt(3), iz(7),fn(3,15),bk
 REAL :: rz,x
 DOUBLE PRECISION :: dz(7),dx
 CHARACTER (LEN=8) :: binary,formtd,bf
!WKBI
 CHARACTER (LEN=5) :: z5(1)
 CHARACTER (LEN=25) :: sfm
 CHARACTER (LEN=29) :: uim
 CHARACTER (LEN=25) :: uwm
 CHARACTER (LEN=23) :: ufm
!WKBNB
 CHARACTER (LEN=80) :: dsnames
 COMMON /dsname/  dsnames(80)
!WKBNE
 COMMON /xmssg /  ufm,uwm,uim,sfm
 COMMON /BLANK /  p1,p2,p3(2),p4
 COMMON /input5/  mcb(1),col,row,FORM,TYPE,MAX,dens
 COMMON /machin/  mach
 COMMON /system/  ibuf,nout,nogo,dum36(36),nbpw
 COMMON /zzzzzz/  rz(1)
 COMMON /packx /  typin,typout,ii,jj,incr
!WKBR EQUIVALENCE      (RZ(1),IZ(1),DZ(1))
 EQUIVALENCE      (rz(1),iz(1),dz(1),z5)
 DATA    binary,  formtd,       subnam,         fn,bk    /  &
     'BINARY','FORMATTD',    4HINPT, 2HT5,   46*2H    /
 DATA    mtrx,    tble,skip   / 4HMTRX, 4HTBLE, 4HSKIP   /
 DATA    debug /  .false.     /
 
!     IF MACHINE IS CDC OR UNIVAC, CALL CDCOPN OR UNVOPN TO OPEN OUTPUT
!     FILE, A SEQUENTIAL FORMATTED TAPE. NO CONTROL WORDS ARE TO BE
!     ADDED TO EACH SEQUENTIAL RECORD. RECORD LENGTH IS 132 CHARACTERS,
!     AN ANSI STANDARD.
 
 bf = binary
 IF (p4 >= 1) bf = formtd
 CALL page1
 WRITE  (nout,10) uim,bf,p1
 10 FORMAT (a29,', MODULE INPUTT5 CALLED BY USER DMAP ALTER, ON ',a8,  &
     ' INPUT FILE,',/5X,'WITH THE FOLLOWING REQUEST.  (P1=', i2,1H))
 IF (p1 == -9) WRITE (nout,20)
 IF (p1 == -3) WRITE (nout,30)
 IF (p1 == -1) WRITE (nout,40)
 IF (p1 ==  0) WRITE (nout,50)
 IF (p1 >  0) WRITE (nout,60) p1
 20 FORMAT (5X,'REWIND TAPE ONLY')
 30 FORMAT (5X,'REWIND AND READ TAPE. PRINT ALL DATA BLOCK NAMES ON ',  &
     'TAPE. AT END, TAPE IS REWOUND', /5X,'AND POSITIONED ',  &
     'PASS TAPE HEADER RECORD')
 40 FORMAT (5X,'REWIND, POSITION PAST TAPE HEADER RECORD, THEN READ ',  &
     'TAPE. AT END, NO REWIND')
 50 FORMAT (5X,'READ TAPE STARTING AT CURRENT POSITION, OR POSITION ',  &
     'PAST THE TAPE HEADER RECORD (FIRST USE OF TAPE).', /5X,  &
     ' NO REWIND AT BEGINNING AND AT END')
 60 FORMAT (5X,'SKIP FORWARD',i4,' DATA BLOCKS (NOT COUNTING TAPE ',  &
'HEADER RECORD) BEFORE READING, AT END NO REWIND')

buf1 = korsz(rz(1)) - ibuf - 1
IF (buf1 <= 0) CALL mesage (-8,0,subnam)
INPUT= p2
opn  = .false.
ll   = 0
p41  =.false.
IF (p4 >=  1) p41 =.true.
p40  =.NOT.p41
p40s =.false.
p41s =.false.
p40d = p40
p41d = p41
p41c = p4 == 2 .AND. nbpw >= 60
IF (p41c) p40d = .false.
IF (p41c) p41d = .false.
col12= 0
p1n  = p1
IF (p1 <  0) p1n = 0
!WKBNB
CLOSE( UNIT=INPUT )
IF ( p4 /= 0  ) GO TO 62
OPEN ( UNIT=INPUT, FILE=dsnames(INPUT), FORM='UNFORMATTED' ,STATUS='UNKNOWN' )
GO TO 65
62    CONTINUE
OPEN ( UNIT=INPUT, FILE=dsnames(INPUT), STATUS='UNKNOWN')
65    CONTINUE
!WKBNE
IF (p1 /= -9) GO TO 70
REWIND INPUT
GO TO 1000

70 DO  i = 1,15
  fn(3,i) = bk
END DO
IF (p1 >= -1 .OR. p1 == -3 .OR. p1 == -9) GO TO 200

WRITE  (nout,90) ufm,p1
90 FORMAT (a23,', MODULE INPUTT5 - ILLEGAL VALUE FOR FIRST PARAMETER'  &
    ,      ' = ',i8, /5X,'ONLY -9, -3 AND GREATER THAN -1 ALLOWED')
100 ERR = -37
120 CALL mesage (ERR,output,subnam)
RETURN

200 IF (p1 == 0) GO TO 500

!     CHECK TAPE ID

REWIND INPUT
ERR = -1
IF (p40) READ (INPUT,    END=420) tapeid,mac,dt,i,k
IF (p41) READ (INPUT,210,END=420) tapeid,mac,dt,i,k
210 FORMAT (4A4,5I8)
IF (tapeid(1) == p3(1) .AND. tapeid(2) == p3(2)) GO TO 230
WRITE  (nout,220) tapeid,p3,mac,dt
220 FORMAT ('0*** WRONG TAPE MOUNTED - TAPEID =',2A4,', NOT ',2A4,  &
    /5X,'MACHINE=',2A4,' DATE WRITTEN-',i4,1H/,i2,1H/,i2)
IF (p1 == -1) GO TO 100
230 IF (k  == p4) GO TO 250
WRITE  (nout,240) uwm,p4
240 FORMAT (a25,', MODULE INPUTT5 4TH PARAMETER SPECIFIED WRONG TAPE',  &
    ' FORMAT.   P4=',i5, /5X,  &
    'INPUTT5 WILL RESET P4 AND TRY TO READ THE TAPE AGAIN.',/)
p4  = k
p40 =.NOT.p40
p41 =.NOT.p41
250 CALL page2 (4)
WRITE  (nout,260) tapeid,mac,dt,i
260 FORMAT (/5X,'MODULE INPUTT5 IS NOW PROCESSING TAPE ',2A4,  &
    ' WHICH WAS WRITTEN BY ',2A4,'MACHINE', /5X,  &
    'ON',i4,1H/,i2,1H/,i2,4X,'SYSTEM BUFFSIZE=',i8)
IF (p40) WRITE (nout,270)
IF (p41) WRITE (nout,280)
270 FORMAT (5X,'TAPE IN BINARY RECORDS',/)
280 FORMAT (5X,'TAPE IN FORMATTED RECORDS',/)
ll = 0
IF (p1 > 0 .OR. p1 == -3) GO TO 300
IF (p1 == -1) GO TO 510
imhere = 290
WRITE  (nout,290) sfm,imhere,p1
290 FORMAT (a25,' @',i5,i10)
GO TO 100

!     TO SKIP P1 MATRIX DATA BLOCKS OR TABLES ON INPUT TAPE (P1 = +N)
!     OR PRINT CONTENTS OF INPUT TAPE (P1 = -3)

300 IF (p40 ) READ (INPUT,    ERR=390,END=420) nc,jb,je
IF (p41s) READ (INPUT,520,ERR=390,END=420) nc,jb,je,( x,j=jb,je)
IF (p41c) READ (INPUT,525,ERR=390,END=420) nc,jb,je,( x,j=jb,je)
IF (p41d) READ (INPUT,530,ERR=390,END=420) nc,jb,je,(dx,j=jb,je)
IF (debug .AND. (nc <= 15 .OR. nc >= col12)) WRITE (nout,540) nc,jb,je,ll
IF (nc < 0) THEN
  GO TO   360
ELSE IF (nc == 0) THEN
  GO TO   340
ELSE
  GO TO   300
END IF

310 IF (p40) READ (INPUT,    ERR=390,END=420) l
IF (p41) READ (INPUT,320,ERR=390,END=420) l,(tabel,j=1,l)
320 FORMAT (i10,24A, /,(26A5))
IF (debug) WRITE (nout,330) l,ll
330 FORMAT (30X,'L AND LL=',2I6)
imhere = 330
IF (l < 0) THEN
  GO TO   360
ELSE IF (l == 0) THEN
  GO TO   340
ELSE
  GO TO   310
END IF
340 IF (p1 /= -3 .AND. ll >= p1) GO TO 360
imhere = 340
ll = ll + 1
BACKSPACE INPUT
IF (p41) BACKSPACE INPUT
IF (ll > 15) GO TO 370
IF (p40) READ (INPUT    ) i,i,i,dx,j,j,j,j,k,k,fn(1,ll),fn(2,ll)
IF (p41) READ (INPUT,560) i,i,i,dx,j,j,j,j,k,k,fn(1,ll),fn(2,ll)
IF (p1 /= -3 .OR. ll <= p1) fn(3,ll) = skip
IF (k > 0 .AND. j >= 1 .AND. j <= 4) GO TO 350

!     FILE IS A TABLE

IF (ll > p1) fn(3,ll) = tble
imhere = 345
GO TO 310

!     FILE IS A MATRIX

350 IF (ll > p1) fn(3,ll) = mtrx
IF (p40) GO TO 300
p41s = .false.
p41d = .false.
p41c =  p4 == 2 .AND. nbpw >= 60
IF (p41c) GO TO 300
IF (j == 1 .OR. j == 3) p41s = .true.
p41d = .NOT.p41s
GO TO 300

360 IF (p1 == -3) GO TO 900
IF (p41) BACKSPACE INPUT
BACKSPACE INPUT
GO TO 510

370 WRITE  (nout,380) uim
380 FORMAT (a29,', INPUTT5, WITH P1= -3, CAN ONLY PRINT UP TO 15 ',  &
    ' FILE NAMES ON ONE INPUT TAPE.', /5X,'TAPE IS POSITIONED',  &
    ' AFTER THE 15TH FILE')
ll = ll - 1
GO TO 920

390 WRITE  (nout,400) ufm,p3,ll,nc,imhere
400 FORMAT (a23,', TAPE ERROR DURING READ/INPUTT5  ',2A4, /5X,  &
    'LL,NC =',2I5,'   IMHERE =',i5)
imhere = 405
IF (p41 .AND. mach == 2) WRITE (nout,410) imhere
410 FORMAT (/5X,'IBM USER - CHECK FILE ASSIGNMENT FOR DCB PARAMETER ',  &
    'OF 132 BYTES',i15)
GO TO  100
420 IF (p1 == -3) GO TO 440
WRITE  (nout,430) ufm,p3,imhere,ll,nc
430 FORMAT (a23,', EOF ENCOUNTERED ON INPUT TAPE ',2A4,5X,  &
    'IMHERE,LL,NC =',3I5)
IF (p1 /= -3) nogo = 1
GO TO  900
440 WRITE  (nout,450) uwm,p3
450 FORMAT (a25,', EOF ENCOUNTERED ON INPUT TAPE ',2A4,'. TAPE DOES ',  &
    'NOT CONTAIN AN ''OUTPUT5 E-O-F'' MARK')
IF (debug) WRITE (nout,460) imhere,ll,nc
460 FORMAT (5X,'IMHERE,LL,NC =',3I5)
GO TO  900

!     P1 = 0,
!     MUST SKIP TAPE HEADER RECORD IF CURRENT TAPE POSITION IS AT THE
!     VERY BEGINNING

500 ll = 0
imhere = 500
IF (p40) READ (INPUT,    ERR=770,END=420) tapeid
IF (p41) READ (INPUT,210,ERR=770,END=420) tapeid
IF (tapeid(1) /= p3(1) .OR. tapeid(2) /= p3(2)) BACKSPACE INPUT

!     COPY MATRIX TO TAPE

imhere = 510
510 IF (p40s) READ(INPUT,    ERR=770,END=910) nc,jb,je,(rz(j),j=jb,je)
IF (p40d) READ(INPUT,    ERR=770,END=910) nc,jb,je,(dz(j),j=jb,je)
IF (p41s) READ(INPUT,520,ERR=770,END=910) nc,jb,je,(rz(j),j=jb,je)
IF (p41c) READ(INPUT,525,ERR=770,END=910) nc,jb,je,(rz(j),j=jb,je)
IF (p41d) READ(INPUT,530,ERR=770,END=910) nc,jb,je,(dz(j),j=jb,je)
520 FORMAT (3I8,/,(10E13.6))
525 FORMAT (3I8,/,(5E26.17))
530 FORMAT (3I8,/,(5D26.17))
IF (debug .AND. (nc <= 15 .OR. nc >= col12))  &
    WRITE (nout,540) nc,jb,je,ll,imhere
540 FORMAT (30X,'NC,JB,JE,LL=',5I6,'=IMHERE')
IF (nc < 0) THEN
  GO TO   800
ELSE IF (nc == 0) THEN
  GO TO   550
ELSE
  GO TO   700
END IF
!             EOF, MATRIX-HEADER, COLUMN-DATA

!     MATRIX OR TABLE HEADER

550 IF (opn) GO TO 810
ll = ll + 1
IF (ll > 15) GO TO 370
BACKSPACE INPUT
IF (p41) BACKSPACE INPUT
j = -1
IF (p40) READ (INPUT,    ERR=570) k,j,j,dx,col,row,FORM,TYPE,  &
    MAX,dens,fn(1,ll),fn(2,ll)
IF (p41) READ (INPUT,560,ERR=570) k,j,j,dx,col,row,FORM,TYPE,  &
    MAX,dens,fn(1,ll),fn(2,ll)
560 FORMAT (3I8,/,d26.17,6I8,2A4)
col12 = col - 12
IF (col12 < 0) col12 = 0
IF (.NOT.debug) GO TO 590
570 WRITE  (nout,580) col,row,FORM,TYPE,MAX,dens,dx,fn(1,ll),fn(2,ll)
580 FORMAT (' COL,ROW,FORM,TYPE,MAX,DENS,DX,FILE=',6I6,d12.3,3X,2A4)
IF (j == -1) CALL mesage (-37,0,subnam)

590 IF (k == 0 .AND. (dens == 0 .OR. TYPE <= 0 .OR. TYPE > 4)) &
!WKBR1    CALL TABLE V (*510,INPUT,LL,MCB,FN(1,LL),P4,BUF1,RZ)  &
CALL tablev (*510,INPUT,ll,mcb,fn(1,ll),p4,buf1,z5)

fn(3,ll) = mtrx
p40s = .false.
p40d = .false.
p41s = .false.
p41d = .false.
p41c =  p4 == 2 .AND. nbpw >= 60
IF (p41c) GO TO 610
IF (p41 ) GO TO 600
IF (TYPE == 1 .OR. TYPE == 3) p40s = .true.
p40d = .NOT.p40s
GO TO 610
600 IF (TYPE == 1 .OR. TYPE == 3) p41s = .true.
p41d = .NOT.p41s
610 IF (debug) WRITE (nout,620) p40,p40s,p40d,p41,p41s,p41d,p41c
620 FORMAT ('0  P40,P40S,P40D,P41,P41S,P41D,P41C = ',7L4)
typin  = TYPE
typout = TYPE
jtyp   = TYPE
IF (TYPE == 3) jtyp = 2
ii   = 1
jj   = row
incr = 1
nwds = row*jtyp
IF (nwds > buf1) CALL mesage (-8,0,subnam)

!     OPEN GINO FILE FOR OUTPUT

imhere = 640
IF (p1 == -3) GO TO 640
rowx   = row
formx  = FORM
output = 200 + ll - p1n
mcb(1) = output
CALL rdtrl (mcb(1))
IF (mcb(1) <= 0) GO TO 750
ERR  = -1
CALL OPEN  (*120,output,rz(buf1),1)
CALL fname (output,NAME)
CALL WRITE (output,NAME,2,1)
opn  = .true.
col  = 0
row  = rowx
FORM = formx
TYPE = typout
MAX  = 0
dens = 0
nck  = 0
GO TO 510

640 WRITE (nout,290) sfm,imhere,p1
CALL mesage (-37,0,subnam)

!     RECOVER INPUT MATRIX, AND WRITE IT OUT BY COLUMN

700 imhere = 700
IF (p1 == -3) GO TO 510
nck = nck + 1
IF (nc /= nck) GO TO 390
IF (jb <= 1) GO TO 720
jb  = (jb-1)*jtyp
DO  j = 1,jb
  rz(j) = 0.0
END DO
720 IF (je >= nwds) GO TO 740
je  = (je*jtyp) + 1
DO  j = je,nwds
  rz(j) = 0.0
END DO
740 CALL pack (rz,output,mcb)
GO TO 510

!     OUTPUT FILE PURGED, SKIP FORWARD FOR NEXT MATRIX ON TAPE

750 IF (p40 ) READ (INPUT    ,ERR=390,END=420) nc,jb,je
IF (p41s) READ (INPUT,520,ERR=390,END=420) nc,jb,je,( x,j=jb,jb)
IF (p41c) READ (INPUT,525,ERR=390,END=420) nc,jb,je,( x,j=jb,jb)
IF (p41d) READ (INPUT,530,ERR=390,END=420) nc,jb,je,(dx,j=jb,jb)
IF (nc > 0) GO TO 750
CALL page2 (2)
WRITE  (nout,760) uwm,fn(1,ll),fn(2,ll)
760 FORMAT (a25,', OUTPUT FILE PURGED.  ',2A4,' FROM INPUT TAPE NOT ',  &
    'COPIED')
!     LL = LL + 1
GO TO 550

770 imhere = -imhere
WRITE  (nout,400) ufm,p3,ll,nc,imhere
WRITE  (nout,780) p40,p41,p40s,p40d,p41s,p41d,p41c
780 FORMAT ('  P40,P41,P40S,P40D,P41S,P41D,P41C =',7L2)
imhere = 770
IF (p41 .AND. mach == 2) WRITE (nout,410) imhere
GO TO 750

!     END OF MATRIX ENCOUNTERED. CLOSE GINO DATA BLOCK WITH REWIND.

800 IF (.NOT.opn) GO TO 840
810 CALL CLOSE (output,1)
opn = .false.
IF (FORM >= 1 .AND. FORM <= 6) GO TO 820
FORM = 1
IF (col /= row) FORM = 2
820 CALL wrttrl (mcb)
CALL fname (output,NAME)
CALL page2 (10)
WRITE  (nout,830) fn(1,ll),fn(2,ll),INPUT,NAME,(mcb(j),j=1,7)
830 FORMAT (/5X,'MATRIX DATA BLOCK ',2A4,' WAS SUCESSFULLY RECOVERED',  &
    ' FROM FORTRAN UNIT',i4,' TO ',2A4, /8X,'GINO UNIT =',i8,  &
    /6X,'NO. OF COLS =',i8, /6X,'NO. OF ROWS =',i8,  /13X,  &
    'FORM =',i8, /13X,'TYPE =',i8, /3X,'NON-ZERO WORDS =',i8, /10X,'DENSITY =',i8)
840 imhere = 840
IF (ll >=  5+p1n) GO TO 860
IF (nc < 0) THEN
  GO TO   850
ELSE IF (nc == 0) THEN
  GO TO   550
ELSE
  GO TO   390
END IF
850 IF (p1 == -3) GO TO 1000
GO TO 900
860 BACKSPACE INPUT
IF (p41) BACKSPACE INPUT
GO TO 920

!     IF NC = -2, THIS IS AN ELEM/GRID ID RECORD WRITTEN BY DUMOD5

900 IF (nc == -2) GO TO 510
910 nc = -3
IF (opn) GO TO 800
IF (fn(3,ll) == bk) ll = ll - 1
IF (ll <= 0) GO TO 970

!     PRINT LIST OF DATA BLOCKS ON FORTRAN TAPE (P1=-3).

920 CALL page2 (ll+9)
WRITE  (nout,930)
IF (p1 /= -3) WRITE (nout,940) INPUT
IF (p1 == -3) WRITE (nout,950) INPUT
WRITE  (nout,960) mac,bf,(j,fn(1,j),fn(2,j),fn(3,j),j=1,ll)
930 FORMAT (/5X,'SUMMARY FROM INPUTT5 MODLUE -')
940 FORMAT (/34X,'FILES RECOVERED FROM FORTRAN UNIT',i5)
950 FORMAT (/34X,'FILE CONTENTS ON FORTRAN UNIT',i5)
960 FORMAT (28X,'(WRITTEN BY ',2A4,' MACHINE ',a8,' RECORDS)', //37X,  &
    'FILE',8X,'NAME',8X,'TYPE', /33X,9(4H----), /, (37X,i3,7X,2A4,6X,a4))
IF (nogo == 1) GO TO 100

IF (p1  /= -3) GO TO 1000
REWIND INPUT
IF (p40) READ (INPUT)
IF (p41) READ (INPUT,210)
GO TO 1000

970 IF (p1 == -3) WRITE (nout,980) uim,INPUT
980 FORMAT (a29,' FROM INPUTT5 MODULE, INPUT TAPE (FORTRAN UNIT',i5,  &
    ') CONTAINS NO DATA BLOCK')

1000 IF (mach == 3) CALL unvcls (p2)
IF (mach == 4) CALL cdccls (p2)
RETURN
END SUBROUTINE inptt5
