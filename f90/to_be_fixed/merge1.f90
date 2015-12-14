SUBROUTINE merge1
     
!     THIS IS THE DMAP MODULE MERGE WHICH MERGES 1 TO 4 PARTITIONS
!     A11, A21, A12, A22, INTO A SINGLE MATRIX -A-.
 
!          **                  **           **                  **
!          *       I            *           *                    *
!          *  A11  I    A12     *           *                    *
!          *       I            *           *                    *
!          * ------+----------- *  BECOMES  *          A         *
!          *       I            *           *                    *
!          *       I            *           *                    *
!          *  A21  I    A22     *           *                    *
!          *       I            *           *                    *
!          **                  **           **                  **
 
!     BASED ON THE ZEROS AND NON-ZEROS IN THE ROW PARTITIONING VECTOR
!     -RP- AND THE COLUMN PARTITIONING VECTOR -CP-.
 
!     DMAP CALLING SEQUENCE.
 
!     MERGE  A11,A21,A12,A22,CP,RP/ A /V,Y,SYM  /V,Y,TYPE/V,Y,FORM/
!                                      V,Y,CPCOL/V,Y,RPCOL  $
 
 IMPLICIT INTEGER (a-z)
 EXTERNAL        rshift  ,andf
 LOGICAL :: cpnull  ,rpnull   ,cphere   ,rphere   ,only    , pass
 DIMENSION       subr(2) ,head(2)  ,aij(4)   ,mcb(7,4) ,mcba(7) ,  &
     elem1(4),elem2(4 ),refus(3) ,BLOCK(80)
 CHARACTER (LEN=27) :: swm
 CHARACTER (LEN=25) :: sfm
 CHARACTER (LEN=29) :: uim
 CHARACTER (LEN=25) :: uwm
 CHARACTER (LEN=23) :: ufm
 COMMON /xmssg / ufm     ,uwm      ,uim      ,sfm      ,swm
 COMMON /system/ sysbuf  ,outpt    ,xxx(37)  ,nbpw
 COMMON /names / rd      ,rdrew    ,wrt      ,wrtrew   ,clsrew  , cls
 COMMON /zblpkx/ elem(4) ,row
 COMMON /prtmrg/ cpsize  ,rpsize   ,cpones   ,rpones   ,cpnull  ,  &
     rpnull  ,cphere   ,rphere   ,icp      ,ncp     , irp     ,nrp
 COMMON /zzzzzz/ z(1)
 COMMON /BLANK / sym     ,TYPE     ,FORM     ,cpcol    ,rpcol   ,  &
     dumfor(3)         ,ireqcl
 DATA    subr  / 4HMERG  ,4HE1   /, eor    / 1      /
 DATA    aij   / 101,102 ,103,104/, cp,rp  / 105,106/     ,a /201/
 DATA    naform/ 4HFORM /,natype /  4HTYPE /,refus/2*3H   ,3HREF /
 
!     OPEN MATRICES TO BE MERGED.  IF ALL ARE PURGED, RETURN IS MADE.
 
 core = korsz(z)
 m = 0
 DO  i = 1,4
   kfile = aij(i)
   mcb(1,i) = kfile
   CALL rdtrl (mcb(1,i))
   IF (mcb(1,i) > 0) THEN
     GO TO    10
   ELSE
     GO TO    20
   END IF
   10 buff = core - sysbuf - 2
   core = buff - 1
   IF (core < 10) CALL mesage (-8,0,subr)
   CALL OPEN (*20,kfile,z(buff),rdrew)
   CALL skprec (kfile,1)
   m = 1
 END DO
 IF (m == 0) RETURN
 buff = core - sysbuf - 2
 core = buff - 1
 IF (core < 10) CALL mesage (-8,0,subr)
 CALL OPEN (*440,a,z(buff),wrtrew)
 CALL CLOSE (a,clsrew)
 
!     CALL TO PARTN2 WILL PROCESS -CP- AND -RP- INTO BIT STRINGS AND
!     DETERMINE SIZES OF PARTITIONS REQUIRED.
 
!     STANDARDISE BLANK COMMON FOR PARTN2 CALLS FROM MERGE1-PARTN1
 
 dumfor(2) = cpcol
 dumfor(3) = rpcol
 CALL partn2 (cp,rp,core,z(buff))
 cpcol = dumfor(2)
 rpcol = dumfor(3)
 
!     IF CPSIZE OR RPSIZE IS 0 AS A RESULT OF A NULL VECTOR (PURGED
!     VECTOR) THERE SIZE IS ESTIMATED HERE FROM THEIR RESPECTIVE
!     PARTITIONS.
 
 IF (cpsize /= 0) GO TO 24
 IF (mcb(1,1) > 0) THEN
   GO TO    23
 ELSE
   GO TO    25
 END IF
 23 cpsize = mcb(2,1)
 GO TO 24
 
 25 IF (mcb(1,2) > 0) THEN
   GO TO    26
 ELSE
   GO TO    24
 END IF
 26 cpsize = mcb(2,2)
 GO TO 24
 
 24 IF (rpsize /= 0) GO TO 29
 IF (mcb(1,1) > 0) THEN
   GO TO    22
 ELSE
   GO TO    21
 END IF
 22 rpsize = mcb(3,1)
 GO TO 29
 
 21 IF (mcb(1,3) > 0) THEN
   GO TO    27
 ELSE
   GO TO    29
 END IF
 27 rpsize = mcb(3,3)
 
!     MATRIX COMPATIBILITY CHECKS.
 
 29 cpzero = cpsize - cpones
 rpzero = rpsize - rpones
 ipr    = 1
 irlcx  = 0
 DO  i = 1,4
   IF (mcb(1,i) > 0) THEN
     GO TO    30
   ELSE
     GO TO    70
   END IF
   30 cols = mcb(2,i)
   rows = mcb(3,i)
   icol = cpzero
   irow = rpzero
   IF (mcb(5,i) == 2 .OR. mcb(5,i) == 4) ipr   = 2
   IF (mcb(5,i) == 3 .OR. mcb(5,i) == 4) irlcx = 2
   IF (i == 3 .OR. i == 4) icol = cpones
   IF (i == 2 .OR. i == 4) irow = rpones
   IF (icol > 0) THEN
     GO TO    40
   ELSE
     GO TO    70
   END IF
   40 IF (irow > 0) THEN
     GO TO    50
   ELSE
     GO TO    70
   END IF
   
!     CHECK PARTITION SIZE WITH PARTITIONING VECTOR DEMANDS.
   
   50 IF (rows == irow .AND. cols == icol) CYCLE
   WRITE  (outpt,60) swm,aij(i),rows,cols,irow,icol
   60 FORMAT (a27,' 2161, PARTITION FILE',i4,' IS OF SIZE',i10,  &
       ' ROWS BY',i10,' COLUMNS.', /5X,'PARTITIONING VECTORS ',  &
       'INDICATE THAT THIS PARTITION SHOULD BE OF SIZE',i10,  &
       ' ROWS BY',i10,' COLUMNS FOR A SUCCESSFUL MERGE.')
 END DO
 
!     CHECK OF FORM VALUE.
 
 nform = FORM
 IF (nform < 1 .OR. nform > 8) GO TO 120
 SELECT CASE ( nform )
   CASE (    1)
     GO TO 80
   CASE (    2)
     GO TO 140
   CASE (    3)
     GO TO 110
   CASE (    4)
     GO TO 80
   CASE (    5)
     GO TO 80
   CASE (    6)
     GO TO 80
   CASE (    7)
     GO TO 110
   CASE (    8)
     GO TO 80
 END SELECT
 
!     FORM = SQUARE
 
 80 IF (cpsize == rpsize) GO TO 140
 90 WRITE  (outpt,100) swm,nform,rpsize,cpsize
 100 FORMAT (a27,' 2162, THE FORM PARAMETER AS GIVEN TO THE MERGE ',  &
     'MODULE IS INCONSISTANT WITH THE SIZE OF THE', /5X,  &
     'MERGED MATRIX, HOWEVER IT HAS BEEN USED.  FORM =',i9,  &
     ' SIZE =',i10,' ROWS BY',i10,' COLUMNS.')
 GO TO 140
 110 IF (cpsize == 1) GO TO 140
 GO TO 90
 120 nform = 2
 IF (rows /= cols .AND. cpsize /= rpsize) GO TO 122
 nform = 1
 IF (sym < 0) nform= 6
 122 IF (FORM == 0 .OR. FORM == nform) GO TO 132
 WRITE  (outpt,130) swm,naform,FORM,refus(3),subr,nform
 130 FORMAT (a27,' 2163, REQUESTED VALUE OF ',a4,i10,2X,a3,'USED BY ',  &
     2A4,'. LOGICAL CHOICE IS',i10)
 132 FORM = nform
 
!     CHECK PARAMETER -TYPE-
 
 140 ntype =   irlcx + ipr
 IF (ntype == TYPE) GO TO 160
 IF (TYPE  ==    0) GO TO 154
 IF (TYPE < 0 .OR. TYPE > 4) GO TO 152
 WRITE (outpt,130) swm,natype,TYPE,refus(1),subr,ntype
 ntype = TYPE
 GO TO 160
 152 WRITE (outpt,130) swm,natype,TYPE,refus(3),subr,ntype
 154 TYPE = ntype
 
!     THE ROW PARTITIONING BIT STRING IS AT THIS POINT CONVERTED TO A
!     CORE VECTOR ONE WORD PER BIT.  EACH WORD CONATINS THE ACTUAL ROW
!     POSITION THE SUB-PARTITON ELEMENT WILL OCCUPY IN THE MERGED
!     MATRIX.
 
 160 iz = nrp + 1
 nz = iz + rpsize - 1
 IF (nz > core) CALL mesage (-8,0,subr)
 IF (.NOT.rpnull .AND. rpones /= 0) GO TO 180
 k = 0
 DO  i = iz,nz
   k = k + 1
   z(i) = k
 END DO
 GO TO 240
 180 k = 0
 zero = iz - 1
 ones = zero + rpzero
 DO  i = irp,nrp
   DO  j = 1,nbpw
     shift = nbpw - j
     bit   = rshift(z(i),shift)
     k = k + 1
     IF (k -  rpsize > 0) THEN
       GO TO   240
     END IF
     190 IF (andf(bit,1) == 0.0) THEN
       GO TO   200
     ELSE
       GO TO   210
     END IF
     200 zero = zero + 1
     z(zero) = k
     CYCLE
     210 ones = ones + 1
     z(ones) = k
   END DO
 END DO
 
!     OPEN OUTPUT FILE AND FILL MCB.
 
 240 CALL OPEN (*440,a,z(buff),wrtrew)
 CALL fname (a,head)
 CALL WRITE (a,head,2,eor)
 CALL makmcb (mcba,a,rpsize,nform,ntype)
 
!     MERGE OPERATIONS.  LOOPING ON OUTPUT COLUMNS OF -A-.
 
 i1 = iz - 1
 i2 = i1 + rpzero
 DO  i = 1,cpsize
   
!     START A COLUMN OUT ON -A-
   
   CALL bldpk (ntype,ntype,a,0,0)
   IF (cpnull) GO TO 250
   il1   = i - 1
   bitwd = il1/nbpw + icp
   shift = nbpw - MOD(il1,nbpw) - 1
   bit   = rshift(z(bitwd),shift)
   IF (andf(bit,1) == 0.0) THEN
     GO TO   250
   ELSE
     GO TO   260
   END IF
   
!     ZERO-S COLUMN (LEFT PARTITONS A11 AND A21 USED THIS PASS)
   
   250 ifile  = 1
   iblock = 1
   GO TO 270
   
!     ONE-S COLUMN (RIGHT PARTITIONS A12 AND A22 USED THIS PASS)
   
   260 ifile  = 3
   iblock = 41
   GO TO 270
   
!     START UNPACKING COLUMN OF EACH PARTITION BEING USED THIS PASS.
   
   270 kfile  = ifile
   kblock = iblock
   mpart  = 0
   DO  j = 1,2
     IF (mcb(1,kfile) > 0) THEN
       GO TO   280
     ELSE
       GO TO   290
     END IF
     280 CALL intpk (*290,mcb(1,kfile),BLOCK(kblock),ntype,1)
     mpart  = mpart  + j
     290 kfile  = kfile  + 1
     kblock = kblock + 20
   END DO
   IF (mpart > 0) THEN
     GO TO   310
   ELSE
     GO TO   420
   END IF
   
!     UNPACK NON-ZEROS FROM EACH OF THE TWO PARTITIONS AS NEEDED UNTIL
!     BOTH PARTITIONS HAVE THIS COLUMN EXHAUSED.
   
   310 eol1  = 1
   eol2  = 1
   nam1  = mcb(1,ifile)
   nam2  = mcb(1,ifile+1)
   ibloc1= iblock
   ibloc2= iblock + 20
   IF (mpart == 1 .OR. mpart == 3) eol1 = 0
   IF (mpart > 1) eol2 = 0
   pass = .false.
   only = .false.
   IF (eol1 > 0.0) THEN
     GO TO   340
   END IF
   320 IF (eol2 > 0.0) THEN
     GO TO   330
   ELSE
     GO TO   350
   END IF
   330 only = .true.
   GO TO 350
   340 only = .true.
   GO TO 360
   
!     UNPACK A NON-ZERO FROM THE ZEROS PARTITION
   
   350 CALL intpki (elem1,irow1,nam1,BLOCK(ibloc1),eol1)
   
!     SET OUTPUT ROW POSITION
   
   jrow  = i1 + irow1
   ipos1 = z(jrow)
   IF (only) GO TO 380
   IF (pass) GO TO 370
   
!     UNPACK A NON-ZERO FROM THE ONE-S PARTITION
   
   360 CALL intpki (elem2,irow2,nam2,BLOCK(ibloc2),eol2)
   
!     SET OUTPUT ROW POSITION
   
   jrow  = i2 + irow2
   ipos2 = z(jrow)
   IF (only) GO TO 400
   pass  = .true.
   
!     OK COMING HERE MEANS THERE IS ONE ELEMENT FORM EACH PARTITION
!     AVAILABLE FOR OUTPUT.  THUS OUTPUT THE ONE WITH THE LOWEST
!     OUTPUT ROW NUMBER.
   
   370 IF (ipos2 < ipos1) GO TO 400
   
!     OUTPUT ELEMENT FROM ZERO-S PARTITION.
   
   380 row     = ipos1
   elem(1) = elem1(1)
   elem(2) = elem1(2)
   elem(3) = elem1(3)
   elem(4) = elem1(4)
   CALL zblpki
   IF (eol1 > 0.0) THEN
     GO TO   390
   ELSE
     GO TO   350
   END IF
   390 IF (only) GO TO 420
   only = .true.
   GO TO 400
   
!     OUTPUT ELEMENT FROM ONES-PARTITION.
   
   400 row     = ipos2
   elem(1) = elem2(1)
   elem(2) = elem2(2)
   elem(3) = elem2(3)
   elem(4) = elem2(4)
   CALL zblpki
   IF (eol2 > 0.0) THEN
     GO TO   410
   ELSE
     GO TO   360
   END IF
   410 IF (only) GO TO 420
   only = .true.
   GO TO 380
   
!     COMPLETE THE COLUMN BEING OUTPUT
   
   420 CALL bldpkn (a,0,mcba)
 END DO
 
!     MERGE IS COMPLETE.  WRAP UP.
 
 CALL CLOSE (a,clsrew)
 CALL wrttrl (mcba)
 440 DO  i = 1,4
   IF (mcb(1,i) > 0) THEN
     GO TO   450
   ELSE
     GO TO   460
   END IF
   450 CALL CLOSE (mcb(1,i),clsrew)
 END DO
 RETURN
 
END SUBROUTINE merge1
