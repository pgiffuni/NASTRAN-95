SUBROUTINE partn1
     
!     THIS IS THE DMAP MODULE PARTN WHICH PARTITIONS A MATRIX -A- INTO
!     FOUR PARTITIONS, SOME OR ALL OF WHICH MAY BE PURGED.
 
 
!                             **                  **
!                             *       I            *
!                             *  A11  I    A12     *
!          **   **            *       I            *
!          *     *            * ------+----------- *
!          *  A  *  BECOMES   *       I            *
!          *     *            *       I            *
!          **   **            *  A21  I    A22     *
!                             *       I            *
!                             **                  **
 
 
!     BASED ON ROW PARTITION MATRIX -RP- AND COLUMN PARTITION MATRIX
!     -CP-
 
!     DMAP SEQUENCE.
 
!     PARTN A,CP,RP/A11,A21,A12,A22/V,Y,SYM/V,Y,TYPE/V,Y,F11/V,Y,F21/
!                                   V,Y,F12/V,Y,F22 $
 
 IMPLICIT INTEGER (a-z)
 EXTERNAL        rshift,andf
 LOGICAL :: cphere,rphere,cpnull,rpnull
 DIMENSION       mcb(7,4),BLOCK(80),subr(2),aij(4),buffs(5),  &
     head(2),mcba(7),refus(3)
 CHARACTER (LEN=27) :: swm
 CHARACTER (LEN=25) :: sfm
 CHARACTER (LEN=29) :: uim
 CHARACTER (LEN=25) :: uwm
 CHARACTER (LEN=23) :: ufm
 COMMON /xmssg / ufm,uwm,uim,sfm,swm
 COMMON /mahcin/ machx
 COMMON /zntpkx/ elem(4),row,eol
 COMMON /names / rd,rdrew,wrt,wrtrew,clsrew,cls
 COMMON /system/ sysbuf,outpt,xxx(37),nbpw
 COMMON /prtmrg/ cpsize,rpsize,cpones,rpones,cpnull,rpnull,cphere,  &
     rphere,icp,ncp,irp,nrp
 COMMON /zzzzzz/ z(1)
 COMMON /BLANK / sym,TYPE,FORM(4),cpcol,rpcol,ireqcl
 DATA    subr  / 4HPART,4HN1    /, a,cp,rp /101,102,103/,  &
     aij   / 201,202,203,204/
 DATA    naform/ 4HFORM /, natype/ 4HTYPE /, refus/2*3H   ,3HREF /
 DATA    eor   / 1      /
 
 core = korsz(z)
 buffs(1) = core - sysbuf - 2
 DO  i = 2,5
   buffs(i) = buffs(i-1) - sysbuf - 2
 END DO
 core = buffs(5) - 1
 IF (core < 10) CALL mesage (-8,0,subr)
 
!     OPEN MATRIX TO BE PARTITIONED.  IF PURGED RETURN IS MADE
 
 buff = buffs(5)
 CALL OPEN (*425,a,z(buff),rdrew)
 CALL skprec (a,1)
 mcba(1) = a
 CALL rdtrl (mcba)
 inform = mcba(4)
 
!     CALL TO PARTN2 WILL PROCESS -CP- AND -RP- INTO BIT STRINGS AND
!     DETERMINE SIZES OF THE PARTITIONS.
 
 buff = buffs(4)
 CALL partn2 (cp,rp,core,z(buff))
 
!     IF RPSIZE OR CPSIZE ARE 0 THEY ARE SET EQUAL TO THE RESPECTIVE
!     SIZE OF A
 
 IF (cpsize == 0) cpsize = mcba(2)
 IF (rpsize == 0) rpsize = mcba(3)
 
!     MATRIX COMPATIBILITY CHECKS
 
 IF (rpsize == mcba(3) .AND. cpsize == mcba(2)) GO TO 40
 WRITE  (outpt,30) swm,mcba(3),mcba(2),rpsize,cpsize
 30 FORMAT (a27,' 2166, MATRIX TO BE PARTITIONED IS OF SIZE',i10,  &
     ' ROWS BY',i10,' COLUMNS.', /5X,'ROW PARTITION SIZE IS',  &
     i10,' COLUMN PARTITION SIZE IS',i10,' (INCOMPATIBLE).')
 
!     PREPARE OUTPUT DATA BLOCKS AS REQUIRED.
 
 40 cpzero = mcba(2) - cpones
 rpzero = mcba(3) - rpones
 
!     CHECK OF TYPE PARAMETER
 
 ntype = mcba(5)
 IF (ntype == TYPE) GO TO 60
 IF (TYPE  ==    0) GO TO 54
 IF (TYPE < 0 .OR. TYPE > 4) GO TO 52
 WRITE  (outpt,50) swm,natype,TYPE,refus(1),subr,ntype
 50 FORMAT (a27,' 2163, REQUESTED VALUE OF ',a4,i10,2X,a3,  &
     'USED BY ',2A4,'. LOGICAL CHOICE IS',i10)
 ntype = TYPE
 GO TO 60
 52 WRITE (outpt,50) swm,natype,TYPE,refus(3),subr,ntype
 54 TYPE = ntype
 
 60 DO  i = 1,4
   FILE = aij(i)
   mcb(1,i) = 0
   cols = cpzero
   rows = rpzero
   IF (i == 3 .OR. i == 4) cols = cpones
   IF (i == 2 .OR. i == 4) rows = rpones
   
!     IF ROWS OR COLS EQUAL ZERO NOTHING IS WRITTEN ON THIS PARTITION
   
   IF (rows == 0 .OR. cols == 0) CYCLE
   buff = buffs(i)
   CALL OPEN (*140,FILE,z(buff),wrtrew)
   CALL fname (FILE,head)
   CALL WRITE (FILE,head,2,eor)
   
!     CHECK OF THE FORM PARAMETER
   
   nform = FORM(i)
   IF (nform < 1 .OR. nform > 8) GO TO 110
   SELECT CASE ( nform )
     CASE (    1)
       GO TO 70
     CASE (    2)
       GO TO 130
     CASE (    3)
       GO TO 100
     CASE (    4)
       GO TO 70
     CASE (    5)
       GO TO 70
     CASE (    6)
       GO TO 70
     CASE (    7)
       GO TO 100
     CASE (    8)
       GO TO 70
   END SELECT
   
!     FORM IMPLIES SQUARE
   
   70 IF (rows == cols) GO TO 130
   80 WRITE  (outpt,90) swm,head,nform,rows,cols
   90 FORMAT (a27,' 2168, THE FORM PARAMETER AS GIVEN TO THE PARTITION',  &
       'ING MODULE FOR SUB-PARTITION ',2A4, /5X,'IS INCONSISTANT'  &
       ,       ' WITH ITS SIZE.  FORM =',i9,' SIZE =',i9,' ROWS BY',i9,  &
       ' COLUMNS.')
   GO TO 130
   
!     DIAGONAL OR ROW MATRIX
   
   100 IF (cols == 1) GO TO 130
   GO TO 80
   
!     NO FORM SPECIFIED THUS IT IS SQUARE IF ROWS = COLS OR RECTANGULAR
!     OTHERWISE.
   
   110 nform = 2
   IF (rows == cols) nform = 1
   IF (sym < 0 .AND. inform == 6 .AND. nform == 1 .AND.  &
       (i == 1 .OR. i == 4)) nform = 6
   IF (FORM(i) == 0) GO TO 128
   jj = 1
   IF (FORM(i) < 1 .OR. FORM(i) > 8) jj = 3
   WRITE (outpt,50) swm,naform,FORM(i),refus(jj),subr,nform
   IF (jj /= 3) nform = FORM(i)
   128 FORM(i) = nform
   
!     TRAILER INITIALIZATION.  BLDPKN WILL SET MCB(2) AND MCB(6) LATER.
   
   130 CALL makmcb (mcb(1,i),FILE,rows,nform,ntype)
   140 CONTINUE
 END DO
 
!     ROW PARTITIONING BIT STRING IS AT THIS POINT CONVERTED TO A CORE
!     VECTOR ONE WORD PER BIT.  EACH WORD CONTAINS THE ROW NUMBER OF THE
!     PARTITION TO WHICH THE ELEMENT OF -A- IS TO BE MOVED TO.  IF THE
!     NUMBER IS NEGATIVE THE ELEMENT IS MOVED TO THE LOWER PARTITIONS
!     AND IF THE NUMBER IS POSITIVE THE ELEMENT IS MOVED TO THE UPPER
!     PARTITION
 
 iz = nrp + 1
 nz = iz  + rpsize - 1
 IF (nz+nbpw > core) CALL mesage (-8,0,subr)
 IF (.NOT.rpnull .AND. rpones /= 0) GO TO 160
 k  = 0
 DO  i = iz,nz
   k  = k + 1
   z(i) = k
 END DO
 GO TO 210
 160 jz   = iz - 1
 zero = 0
 ones = 0
 
!     NOTE THIS LOGIC WORKS ON CRAY WITH 48 OF 64 BIT INTEGER WORD
 
 DO  i = irp,nrp
   DO  j = 1,nbpw
     shift = nbpw - j
     bit   = rshift(z(i),shift)
     jz    = jz + 1
     IF (andf(bit,1) == 0.0) THEN
       GO TO   170
     ELSE
       GO TO   180
     END IF
     170 zero  = zero + 1
     z(jz) = zero
     CYCLE
     180 ones  = ones - 1
     z(jz) = ones
   END DO
 END DO
 
!     LOOP ON ALL THE COLUMNS OF -A-.
 
 210 izm1  = iz - 1
 DO  i = 1,cpsize
   IF (cpnull) GO TO 220
   il1   = i - 1
   bitwd = il1/nbpw + icp
   shift = nbpw - MOD(il1,nbpw) - 1
   bit   = rshift(z(bitwd),shift)
   IF (andf(bit,1) == 0.0) THEN
     GO TO   220
   ELSE
     GO TO   230
   END IF
   
!     ZERO-S COLUMN (LEFT PARTITIONS A11 AND A21)
   
   220 ifile  = 1
   iblock = 1
   GO TO 240
   
!     ONE-S COLUMN (RIGHT PARTITIONS A12 AND A22)
   
   230 ifile  = 3
   iblock = 41
   GO TO 240
   
!     START COLUMNS OF THE 2 AIJ PARTITIONS.
   
   240 kfile  = ifile
   kblock = iblock
   m = 0
   DO  j = 1,2
     IF (mcb(1,kfile) > 0) THEN
       GO TO   250
     ELSE
       GO TO   260
     END IF
     250 CALL bldpk (ntype,mcb(5,kfile),mcb(1,kfile),BLOCK(kblock),1)
     m = 1
     260 kfile  = kfile + 1
     kblock = kblock + 20
   END DO
   IF (m == 0) THEN
     GO TO   390
   END IF
   
!     START THE I-TH COLUMN OF THE MATRIX BEING PARTITIONED -A-.
   
   280 CALL intpk (*350,a,0,ntype,0)
   
!     LOOP ON NON-ZEROS OF THE COLUMN
   
   290 IF (eol > 0.0) THEN
     GO TO   350
   END IF
   
!     PICK UP A NON-ZERO ELEMENT
   
   300 CALL zntpki
   
!     DETERMINE ROW POSITION AND FILE DESTINATION.
   
   l = izm1 + row
   IF (z(l) < 0.0) THEN
     GO TO   320
   END IF
   
!     ZERO-S ROW PARTITION.
   
   310 jrow   = z(l)
   kfile  = ifile
   kblock = iblock
   GO TO 330
   
!     ONE-S ROW PARTITION.
   
   320 jrow   = -z(l)
   kfile  = ifile  + 1
   kblock = iblock + 20
   
!     OUTPUT THE ELEMENT.
   
   330 IF (mcb(1,kfile) > 0) THEN
     GO TO   340
   ELSE
     GO TO   290
   END IF
   340 CALL bldpki (elem,jrow,mcb(1,kfile),BLOCK(kblock))
   GO TO 290
   
!     COMPLETE COLUMNS OF THE 2 AIJ PARTITIONS BEING WORKED ON.
   
   350 kfile  = ifile
   kblock = iblock
   DO  j = 1,2
     IF (mcb(1,kfile) > 0) THEN
       GO TO   360
     ELSE
       GO TO   370
     END IF
     360 CALL bldpkn (mcb(1,kfile),BLOCK(kblock),mcb(1,kfile))
     370 kfile  = kfile  + 1
     kblock = kblock + 20
   END DO
   CYCLE
   
!     COLUMN NOT BEING OUTPUT TO ANY PARTITIONS AT ALL THUS SKIP IT.
   
   390 CALL skprec (a,1)
   
 END DO
 
!     WRAP UP.
 
 CALL CLOSE (a,clsrew)
 DO  i = 1,4
   IF (mcb(1,i) > 0) THEN
     GO TO   410
   ELSE
     GO TO   420
   END IF
   410 CALL wrttrl (mcb(1,i))
   CALL CLOSE (mcb(1,i),clsrew)
   420 CONTINUE
 END DO
 425 RETURN
END SUBROUTINE partn1
