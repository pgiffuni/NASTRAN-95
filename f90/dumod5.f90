SUBROUTINE dumod5
     
!     MSFC ROUTINE, TO CONVERT NASTRAN TABULAR DATA BLOCKS INTO 2-
!     DIMENSIONAL DATA BLOCKS (S.P. REAL ONLY) FOR CONVENIENCE IN
!     MANIPULATION AND OUTPUT, SPECIALLY TO BE USED WITH OUTPUT5 AND
!     INPUT5.
 
!     THIS VERSION WAS MODIFIED BY R. MOORE/MSFC IN JAN. 1989
!     TO ALLOW SELECTION OF EITHER 8 OR 16 VALUES PER ELEMENT BY
!     USING A 7TH PARAMETER ON DMAP
 
!     DUMMOD5  T1,T2,T3,T4,T5/O1,O2,O3,O4,O5/C,N,P1/C,N,P2/C,N,P3
!              C,N,P4/C,N,P5/C,N,Q/C,N,R $
 
!     TI  = INPUT GINO FILE, OEF1, OQG1 OR SIMILAR TYPE OF TABULAR
!           DATA BLOCKS
!     OI  = OUTPUT GINO DATA BLOCK, PACKED, BUT NOT QUITE A REGULAR
!           NASTRAN MATRIX BLOCK, SEE PICTURE BELOW
!           IF OI IS PURGED (NOT PRESENT), MATRIX BLOCK IS WRITTEN OUT
!           TO FORTRAN UNIT 15 (INP1) DIRECTLY, IN BINARY RECORDS,
!           BANDED MATRIX FORM (FROM FIRST NON-ZERO TO LAST NON-ZERO
!           ELEMENTS), D.R. OR S.P.
!     PI  = TI TABLE IS MAPPED INTO A PI X 8 2-DIMENSIONAL BLOCKS.
!           EACH BLOCK IS PACKED AS A COLUMN OF A MATRIX
!     Q   = ELEMENT/GRID POINT ID PRINT-PUNCH CONTROL
!         = -1, NO PRINT AND NO PUNCH
!         =  0, PRINT ONLY, NO PUNCH
!         = +1, BOTH PRINT AND PUNCH
!         = /2/ CONTENTS OF OUTPUT TAPE, INP1, WILL BE PRINTED OUT
!     R   = SWITCH TO CHANGE FROM 8 TO 16 VALUES IN TABLE MAPPING
!           DEFAULT = 0 WHICH SETS TO 8.    R = 1 SETS IT TO 16
 
!     CDC USER ONLY - FORTRAN UNIT 11 (UT1) IS USED INSTEAD OF 15 (INP1)
 
 
!           |<------ 8 OR 16 ------->|
!           ==========================
!         / I                        I \
!        /  I------- TABULAR --------I  \
!       P1  I         DATA           I  BLOCK 1 (MATRIX COLUMN 1)
!        \  I-------- BLOCKS --------I  /
!         \ I                        I /
!           ==========================
!         / I                        I \
!        /  I------------------------I  \
!       P1  I                        I  BLOCK 2 (MATRIX COLUMN 2)
 
!     WRITTEN BY SOMEBODY FOR MARSHALL SPACE FLIGHT CENTER (MSFC).
!     MODIFIED BY G.CHAN/UNISYS TO EMPLOY OPEN-CORE SPACE INSTEAD OF
!     THE FIXED DIMENSION ARRAYS, AND TO EXPAND FROM ONE INPUT DATA
!     BLOCK TO FIVE. IF A CORRESPONDING OUTPUT FILE IS MISSING OR
!     PURGED, THE DATA BLOCKS ARE WRITTEN DIRECTLY TO FORTRAN TAPE
!     (UNIT 15, INP1) USING OUTPUT5 BINARY FORMAT.
 
!     CONTENTS OF INP1 TAPE IF IT IS WRITTEN -
 
!         RECORD   WORD     CONTENT                           TYPE
!         ------  ------   ----------------------------------------
!            0              TAPE HEADER RECORD
!                   1-2     'XXXXXXXX', TAPE ID              2*BCD
!                   3-4     MACHINE TYPE                     2*BCD
!                   5-7     DATE                             3*INT
!                    8      SYSTEM BUFFSIZE                    INT
!                    9      0 (BINARY TAPE)                    INT
!            1              FIRST MATRIX HEADER
!                    1      0                                  INT
!                   2,3     1,1                              2*INT
!                    4      A DOUBLE PRECISION ZERO           D.P.
!                   5-10    6 WORDS FROM MATRIX TRAILER      6*INT
!                           (COL,ROW,FORM,TYPE,MAX,DENSITY-
!                            TYPE=1 OR 3, DENSITY=1)
!                  11-12    MATRIX DMAP NAME                 2*BCD
!            2       1      1 (FIRST COLUMN ID)                INT
!                    2      LOCATION OF FIST NON-ZERO ELEMENT  INT
!                    3      LOCATION OF LAST NON-ZERO ELEMENT  INT
!                   4-N     S.P. DATA                         REAL
!            3       1      2 (SECOND COLUMN ID)               INT
!                   2-N     SAME AS RECORD 1
!            :      1-N     REPEAT FOR MORE COLUMNS
 
!            X       1      X (X-TH COLUMN ID, A NUL COLUMN)   INT
!                   2-3     1,1                                INT
!                   4-5     0.0,0.0                           REAL)
 
!            M      1-N     LAST COLUMN, SAME AS RECORD 1
!           M+1      1      -1 (ELEM) OR -2 (GRID)             INT
!                    2      1                                  INT
!                    3      LENGTH OF ELEM./GRID ID LIST, L    INT
!                  4-L+4    LIST OF ELEMENT OR GRID IDS        INT
 
!           M+2             SECOND MATRIX HEADER
!            :       :      REPEAT 1 THRU (M+1) FOR THE SECOND MATRIX
 
!            :       :      REPEAT, UP TO 5 OUTPUT DATA BLOCKS PER TAPE
 
!     COMMENTS FROM G.C. -
!     (1) THIS MODULE IS VERY LIMITED IN SCOPE. IT HANDLES ONLY SOME
!         SPECIAL TYPES OF TABULAR INPUT DATA BLOCKS. THE (PI X 8) MATRI
!         SPACE IS FOR PRINT/PUNCH PURPOSE. THE ORIGINAL PROGRAM SEEMS
!         TO BE WRITTEN TO MEET A PARTICULAR JOB REQUIREMENT.
 
!     (2) CURRENT MODULE HANDLES ONLY SINGLE PRECISION DATA
 
!     (3) THE PROCEDURE TO READ AND/OR WRITE THE TAPE IS COMMONLY USED
!         AMONG INPUTT5, OUTPUT5, AND DUMMOD5. ANY PROCEDURE CHANGE
!         SHOULD BE MADE TO ALL THREE MODULES.
 
 IMPLICIT INTEGER (a-z)
 LOGICAL :: NONE,     debug
 INTEGER :: NAME(2),  mcb(7),   trl(7),  iz(8),  temp(10),  &
     eg(2),    ir(5001), id(5001),unvc(2),mt(2),  &
     infile(2),outfil(2),date(3), SAVE(2,5)
 REAL :: z,        epsi
 DOUBLE PRECISION :: dzero,    dtemp
 CHARACTER (LEN=29) :: uim
 CHARACTER (LEN=25) :: uwm
 CHARACTER (LEN=23) :: ufm
!WKBNB
 CHARACTER (LEN=80) :: dsnames
 COMMON /dsname/ dsnames(80)
!WKBNE
 COMMON /xmssg / ufm,      uwm,      uim
 COMMON /zzzzzz/ z(1)
 COMMON /machin/ mach,     ijhalf(3),mchnam
 COMMON /system/ ibuf,     nout,     dumm(88),lpch
 COMMON /packx / typin,    typout,   ii,      jj,     incr
 COMMON /BLANK / p(5),     q,        r
 EQUIVALENCE     (z(1),iz(1)),       (date(1),dumm(13))
!WKBI
 DATA    ifirst/0/
 DATA    tape,   irdlmt,   id,       im, ie,  xx,     epsi    /  &
     15,     5000,     5001*0,   1H-,1H=, 4HXXXX, 1.0E-30 /
 DATA    zero,   one,      eg,                NAME            /  &
     0,      1,        4HELEM,   4HGRID,  4HDUMO, 4HD5    /
 DATA    unvc,   mt  /     4HUNIV,   4HAC  ,  2*4H            /
 DATA    debug,  dzero,    SAVE  /   .false., 0.d0,   10*1H   /
 
 IF (mach == 12) tape = 11
 CALL page
 WRITE  (nout,5) p,q,r
 5 FORMAT ('0*** MODULE DUMMOD5 CALLED BY USER DMAP ALTER.', /5X,  &
     'PARAMETERS ARE    P=',5(i5,1H,),5X,'Q=',i5,5X,'R=',i4,/)
 i6_or_8 = 8
 IF (r == 1) i6_or_8 = 16
 incr  = 1
 typin = 1
 typout= 1
 ii    = 1
 tapx  =-1
 tapp  =-1
 core  = korsz(z)
 buf1  = core - ibuf + 1
 buf2  = buf1 - ibuf
 core  = buf2 - 1
 half  = core/2
 half1 = half + 1
!WKBNB
 IF ( ifirst /= 0 ) GO TO 1
 CLOSE ( UNIT=tape)
 OPEN ( UNIT=tape, FILE=dsnames(tape), FORM='UNFORMATTED', STATUS='UNKNOWN' )
 ifirst = 1
 1     CONTINUE
!WKBNE
 
 DO  loop = 1,5
   INPUT = 100 + loop
   outpt = 200 + loop
   trl(1)= INPUT
   CALL rdtrl (trl(1))
   IF (trl(1) <= 0) CYCLE
   CALL fname (INPUT,infile)
   
!     INPUT DATA PRECISION TYPE IS S.P. ONLY
   
   TYPE = 1
   
   IF (p(loop) <= 0) p(loop) = pv
   pv = p(loop)
   jj = p(loop)*i6_or_8
   DO  j = 1,jj
     z(j+half) = 0.0
   END DO
   CALL gopen (INPUT,z(buf1),0)
   mcb(1) = outpt
   CALL rdtrl (mcb)
   NONE = .false.
   IF (mcb(1) <= 0) NONE = .true.
   IF (NONE) GO TO 15
   CALL gopen  (outpt,z(buf2),1)
   CALL fname  (outpt,outfil)
   CALL makmcb (mcb,outpt,0,2,1)
   GO TO 20
   15 tapx = tapx + 1
   IF (tapx <= 0) GO TO 20
   SAVE(1,tapx) = infile(1)
   SAVE(2,tapx) = infile(2)
   20 i    = 1
   nxzh = 0
   nxir = 0
   CALL READ (*290,*30,INPUT,temp,10,1,m)
   nwds = temp(10)
   neltp= temp( 3)
!     IF (NELTP.GE.11 .AND. NELTP.LE.14) GO TO 320
!               CELAS1            CELAS4
   GO TO 60
   30 CALL mesage (-37,0,NAME)
   40 CALL READ (*290,*50,INPUT,temp,10,1,m)
   nwds = temp(10)
   IF (temp(3) /= neltp) GO TO 60
   GO TO 130
   50 CALL mesage (-61,INPUT,NAME)
!  60 IF (TEMP(3).GE.11 .AND. TEMP(3).LE.14) GO TO 320
!                 CELAS1              CELAS4
   60 CONTINUE
   newlt = temp(3)
   nwds1 = nwds - 1
   nwds2 = nwds - 2
   DO  l = 1,jj
     z(l) = 0.0
   END DO
   DO  l = 1,irdlmt
     ir(l) = 0
   END DO
   CALL READ (*330,*350,INPUT,ir(1),1,0,m)
   kount = 0
   DO  jsq = 1,irdlmt
     kount = kount + 1
     loc = nwds1*jsq - nwds2
     CALL READ (*330,*350,INPUT,z(loc),nwds1,0,m)
!     LAST = LOC + NWDS1 - 1
     last = kount*i6_or_8
     CALL READ (*330,*100,INPUT,ir(jsq+1),1,0,m)
   END DO
   100 m   = nwds*kount
   ijk = 0
   DO  j = 1,m,nwds
     ijk = ijk + 1
     nrop  = (ir(ijk)-1)/10
     locid = nxir + ijk
     id(locid) = nrop*100 + newlt
     loca = (ijk*i6_or_8) - (i6_or_8 -1) + nxzh
     lj  = nwds1*ijk - nwds1
     kk  = loca + nwds + half
     IF (kk > core) CALL mesage (-8,0,NAME)
     DO  jm = 1,nwds1
       z(loca+jm-1+half) = z(lj+jm)
     END DO
   END DO
   nxir = nxir + jsq
   nxzh = nxzh + last
   GO TO 40
   130 IF (q < 1) GO TO 150
   is = im
   kk = half + nxzh
   WRITE  (nout,140) is,i,(z(j),j=half1,kk)
   140 FORMAT ('  COLUMN',a1,i5, /,(2X,8E16.6))
   150 i = i + 1
   IF (NONE) GO TO 180
   CALL pack (z(half1),outpt,mcb)
   GO TO 270
   160 IF (tapx > 0) GO TO 170
   
!     WRITE TAPE HEADER AND MATRIX HEADER
!     (CURRENTLY, OUTPUT TAPE IS WRITTEN OUT IN SINGLE PRECISION ONLY)
!     CHANGE IN 89 VERSION -
!     MUST SET MATRIX DENSITY IN MATRIX TRAILER TO NON-ZERO IF INPUT5
!     IS TO BE USED
   
   tapx = 1
   SAVE(1,tapx) = infile(1)
   SAVE(2,tapx) = infile(2)
   mt(1) = mchnam
   IF (mach /= 3) GO TO 162
   mt(1) = unvc(1)
   mt(2) = unvc(2)
   162 WRITE (tape) xx,xx,mt,date,ibuf,zero
   IF (debug) WRITE (nout,165) xx,xx,mt,date,ibuf,zero
   165 FORMAT ('0+++TAPE HEADER/DUMMOD5-',/3X,2A4,1X,2A4,3I4,2I6)
   170 IF (tapx == tapp) GO TO 190
   tapp   = tapx
   trl(5) = typout
   trl(7) = 1
   WRITE (tape) zero,one,one,dzero,(trl(k),k=2,7),infile
   IF (debug) WRITE (nout,175) zero,one,one,dzero,(trl(k),k=2,7), infile
   175 FORMAT (' +++MATRIX HEADER/DUMMOD5- ',3I5,d8.0,6I5,1X,2A4)
   GO TO 190
   
   180 ASSIGN 270 TO retn
   190 DO  jb = ii,jj
!WKBNB 8/94 ALPHA-VMS
     itype = numtyp( z(jb+half) )
     IF ( itype <= 1 ) CYCLE
!WKBNE 8/94 ALPHA-VMS
     IF (ABS(z(jb+half)) > epsi) GO TO 210
   END DO
   WRITE (tape) i,one,one,(zero,j=1,TYPE)
   IF (debug) WRITE (nout,205) i,one,one,(zero,j=1,TYPE)
   205 FORMAT (' +++ZEROS/DUMMOD5- ',7I5)
   GO TO 265
   210 je = jj
   DO  j = ii,jj
!WKBNB 8/94 ALPHA-VMS
     itype = numtyp( z(je+half) )
     IF ( itype <= 1 ) GO TO 220
!WKBNE 8/94 ALPHA-VMS
     IF (ABS(z(je+half)) > epsi) GO TO 230
     220 je = je - 1
   END DO
   230 SELECT CASE ( TYPE )
     CASE (    1)
       GO TO 260
     CASE (    2)
       GO TO 240
     CASE (    3)
       GO TO 240
     CASE (    4)
       GO TO 250
   END SELECT
   240 IF (MOD(jb,2) == 0) jb = jb - 1
   IF (MOD(je,2) == 1) je = je + 1
   GO TO 260
   250 j = MOD(jb,4)
   IF (j == 0) j = 4
   jb = jb - j + 1
   j  = MOD(je,4)
   IF (j == 0) j = 4
   je = je - j + 4
   260 WRITE (tape) i,jb,je,(z(j+half),j=jb,je)
   IF (debug) WRITE (nout,262) i,jb,je
   262 FORMAT (' +++DATA RECORD/DUMMOD5- ',3I5)
   265 GO TO retn, (270,370)
   
   270 DO  l = 1,jj
     z(l+half) = 0.0
   END DO
   nxzh = 0
   nxir = 0
   GO TO 60
   290 IF (q < 0) GO TO 300
   is = ie
   kk = half + nxzh
   WRITE (nout,140) is,i,(z(j),j=half1,kk)
   300 ASSIGN 370 TO retn
   IF (NONE) GO TO 160
   CALL pack (z(half1),outpt,mcb)
   mcb(3) = jj
   CALL wrttrl (mcb)
   IF (q == 2) WRITE (nout,310) (mcb(j),j=1,5)
   310 FORMAT (/2X,'MCB=',6I8)
   GO TO 370
! 320 CALL READ (*330,*40 ,INPUT,IR(1),1,0,M)
!     CALL READ (*330,*350,INPUT, Z(1),1,0,M)
!     Z(1) = 0.0
!     GO TO 320
   330 WRITE  (nout,340) infile
   340 FORMAT (/5X,'*** EOF ENCOUNTERED ON INPUT ',2A4,' DATA BLOCK')
   GO TO 440
   350 WRITE  (nout,360) infile
   360 FORMAT (/5X,'*** INPUT ',2A4,'DATA BLOCK IS EMPTY')
   GO TO 440
   370 IF (.NOT.NONE) WRITE (nout,380) uim,infile,outfil
   380 FORMAT (a29,', MODULE DUMMOD5 SUCCESSFULLY PROCESSED TABULAR ',  &
       'DATA FROM ',2A4,' TO DATA BLOCK ',2A4, /5X, 'IN GINO PACKED FORM')
   IF (NONE) WRITE (nout,390) uim,infile,tape
   390 FORMAT (a29,', MODULE DUMMOD5 SUCCESSFULLY COPIED TABULAR DATA ',  &
       'FROM ',2A4,' TO OUTPUT TAPE', /5X,  &
       '(FORTRAN UNIT',i4,') IN BANDED MATRIX FORM')
   IF (q > 0) WRITE (lpch,400) (id(j),j=1,nxir)
   400 FORMAT (8I10)
   l = eg(1)
   IF (newlt > 0) GO TO 420
   l = eg(2)
   DO  j = 1,nxir
     id(j) = id(j)/100
   END DO
   420 WRITE  (nout,430) l,infile,(id(j),j=1,nxir)
   430 FORMAT (//5X,a4,'-ID ARRAY FOLLOWS/FROM ',2A4, (/5X,15I8))
   IF (.NOT.NONE) GO TO 440
   i = -1
   IF (newlt == 0) i = -2
   WRITE (tape) i,one,nxir,(id(j),j=1,nxir)
   IF (debug) WRITE (nout,435) i,one,nxir
   435 FORMAT (' +++ELEM/GRID ID RECORD/DUMMOD5- ',3I5)
   440 CONTINUE
   CALL CLOSE (INPUT,1)
   IF (.NOT.NONE) CALL CLOSE (outpt,1)
 END DO
 
 IF (tapx <= 0) GO TO 590
 WRITE  (nout,455) uim,tape,(SAVE(1,j),SAVE(2,j),j=1,tapx)
 455 FORMAT (a29,', FOLLOWING DATA BLOCKS WERE COPIED TO FORTRAN UNIT',  &
     i3,' BY MODULE DUMMOD5', /5X,  &
     'USING UNFORMATTED (BINARY) WRITE', /6X,5(2A4,3X))
 ENDFILE tape
 REWIND tape
 
!     TO READ THE OUTPUT TAPE, Q=/2/
 
 IF (IABS(q) < 2) GO TO 590
 CALL page1
 k = 1
 READ (tape,END=575) mcb,j,i
 WRITE  (nout,460) mcb,j
 460 FORMAT (//,'  TAPEID=',2A4,'   FROM ',a4,a2,' MACHINE,  DATE',i5,  &
     1H/,i2,1H/,i2,'  BINARY TAPE.   BUFFSIZE=',i7//)
 470 READ (tape,END=580) i,jb,je,(z(j),j=jb,je)
 IF (i < 0) THEN
   GO TO   560
 ELSE IF (i == 0) THEN
   GO TO   480
 ELSE
   GO TO   500
 END IF
 480 BACKSPACE tape
 READ (tape,END=580) i,jb,je,dtemp,(iz(j),j=1,8)
 WRITE  (nout,490) k,iz(7),iz(8),(iz(j),j=1,6)
 490 FORMAT (//,'  DATA BLOCK',i3,3X,2A4,'  TRAILER=',6I5)
 k = k + 1
 GO TO 470
 500 WRITE  (nout,510) i,jb,je,(z(j),j=jb,je)
 510 FORMAT (//,'  COLUMN RECORD =',i3,'   JB,JE =',2I5,/,(1X,10E13.6))
 GO TO 470
 
 560 l = eg(-i)
 WRITE  (nout,570) l,(iz(j),j=jb,je)
 570 FORMAT (//2X,a4,'-ID LIST -',/,(1X,10I10))
 GO TO 470
 575 WRITE  (nout,577)
 577 FORMAT (//,'  EMPTY TAPE')
 580 REWIND tape
 590 CONTINUE
 RETURN
END SUBROUTINE dumod5
