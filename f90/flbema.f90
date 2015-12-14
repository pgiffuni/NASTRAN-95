SUBROUTINE flbema (TYPE)
     
!     ASSEMBLES THE AF OR DKGG MATRIX UTITLIZING THE ELEMENT
!     MATRICES GENERATED IN FLBEMG
 
!     TYPE = 1  AFF MATRIX
!     TYPE = 2  DKGG MATRIX
 
 
 INTEGER, INTENT(IN)                      :: TYPE
 LOGICAL :: error    ,skip
 INTEGER :: geom2    ,ect      ,bgpdt    ,sil      ,mpt      ,  &
     geom3    ,cstm     ,uset     ,eqexin   ,usetf    ,  &
     usets    ,af       ,dkgg     ,fbelm    ,frelm    ,  &
     conect   ,afmat    ,afdict   ,kgmat    ,kgdict   ,  &
      outmat   ,xmat     ,xdict    ,z        ,  &
     FILE     ,NAME(2)  ,mcb(7)   ,alloc(3) ,dict(2)  ,  &
     typin    ,typout   ,rowsil(4),colsil(12)         ,  &
     optc     ,optw     ,rd       ,rdrew    ,wrt      ,  &
     wrtrew   ,rew      ,norew    ,terms(288)
 CHARACTER (LEN=25) :: sfm
 CHARACTER (LEN=29) :: uim
 CHARACTER (LEN=25) :: uwm
 CHARACTER (LEN=23) :: ufm
 COMMON /xmssg / ufm      ,uwm      ,uim      ,sfm
 COMMON /flbfil/ geom2    ,ect      ,bgpdt    ,sil      ,mpt      ,  &
     geom3    ,cstm     ,uset     ,eqexin   ,usetf    ,  &
     usets    ,af       ,dkgg     ,fbelm    ,frelm    ,  &
     conect   ,afmat    ,afdict   ,kgmat    ,kgdict
 COMMON /zzzzzz/ z(1)
 COMMON /flbptr/ error    ,icore    ,lcore    ,ibgpdt   ,nbgpdt   ,  &
     isil     ,nsil     ,igrav    ,ngrav    ,igrid    ,  &
     ngrid    ,ibuf1    ,ibuf2    ,ibuf3    ,ibuf4    , ibuf5
 COMMON /system/ sysbuf   ,nout
 COMMON /names / rd       ,rdrew    ,wrt      ,wrtrew   ,rew      , norew
 DATA    NAME  / 4HFLBE,4HMA   /
 
 
!     ASSIGN FILES DEPENDING ON TYPE
 
 SELECT CASE ( TYPE )
   CASE (    1)
     GO TO 2
   CASE (    2)
     GO TO 4
 END SELECT
 
!     AF MATRIX
 
 2 outmat = af
 xmat   = afmat
 xdict  = afdict
 GO TO 6
 
!     DKGG MATRIX
 
 4 outmat = dkgg
 xmat   = kgmat
 xdict  = kgdict
 
!     ALLOCATE COLUMN POINTER VECTOR IN TOP OF CORE
 
 6 mcb(1) = uset
 CALL rdtrl(mcb)
 luset = mcb(3)
 icol  = 1
 ncol  = luset
 DO  i = 1,ncol
   z(i) = 0
 END DO
 
!     INITILIZE OPEN AND CLOSE OPTIONS
 
 optw = wrtrew
 optc = norew
 
!     POSITION CONNECT FILE TO PROPER RECORD
 
 FILE = conect
 CALL OPEN (*1001,conect,z(ibuf1),rdrew)
 IF (TYPE == 2) CALL skpfil (conect,1)
 CALL fwdrec (*1002,conect)
 CALL CLOSE (conect,norew)
 
!     INITIALIZE PACK - UNPACK DATA
 
 typin  = 2
 typout = 2
 mcb(1) = outmat
 mcb(2) = 0
 mcb(3) = luset
 mcb(4) = 3 - TYPE
 mcb(5) = typout
 mcb(6) = 0
 mcb(7) = 0
 
!     SET UP CORE POINTERS
 
 icore = ncol  + 1
 lcore = ibuf2 - 1
 ncore = lcore - icore
 IF (ncore < 200) GO TO 1008
 
 skip  = .false.
 ilcol = 0
 
 
!     ALLOCATE ALL AVALABLE CORE FOR THIS PASS BY USE OF CONECT FILE
 
 30 ifcol = ilcol + 1
 jcore = icore
 FILE  = conect
 
 CALL gopen (conect,z(ibuf1),rd)
 
 IF (skip) GO TO 60
 50 CALL READ (*70,*1008,conect,alloc,3,1,n)
 
 60 isil     = alloc(1)
 z(isil)  = jcore
 z(jcore) = jcore + 1
 jcore    = jcore + 1 + alloc(2) + 2*alloc(3)
 IF(jcore > lcore ) GO TO 80
 ilcol = isil
 GO TO 50
 
!     END OF RECORD ON CONECT - ALL COLUMNS ALLOCATED
 
 70 ilcol = luset
 optc  = rew
 GO TO 90
 
!     INSUFFICIENT CORE FOR NEXT COLUMN - SET FLAG TO SAVE CURRENT
!     CONECT ALLOCATION RECORD
 
 80 skip = .true.
 
 90 CALL CLOSE (conect,optc)
 
!     OPEN DICTIONARY AND MATRIX FILES AND PREPARE TO MAKE PASS
 
 CALL gopen (xdict,z(ibuf1),rdrew)
 CALL gopen (xmat,z(ibuf2),rdrew)
 icpos = 0
 
!     READ XDICT ENTRY AND DETERMINE IF COLUMN IS IN CORE FOR THIS
!     PASS
 
 100 FILE = xdict
 CALL READ (*1002,*200,xdict,dict,2,0,n)
 isil = dict(1)
 IF (isil < ifcol .OR. isil > ilcol) GO TO 100
 
!     THE COLUMN IS IN CORE - OBTAIN MATRIX DATA FROM XMAT FILE IF
!     WE DO NOT ALREADY HAVE IT
 
 IF (dict(2) == icpos) GO TO 150
 icpos = dict(2)
 FILE  = xmat
 CALL filpos (xmat,icpos)
 CALL READ (*1002,*1003,xmat,rowsil,4,0,n)
 CALL READ (*1002,*1003,xmat,colsil,4,0,n)
 nrow  = 4
 IF (rowsil(4) < 0) nrow = 3
 ncol  = 4
 IF (colsil(4) < 0) ncol = 3
 CALL READ (*1002,*110,xmat,terms,289,0,nwds)
 icode = 1
 GO TO 8010
 
!     EXPAND COLSIL TO INCLUDE ALL SILS
 
 110 IF(nwds < 162) GO TO 130
 DO  i = 1,4
   j = 4 - i
   colsil(3*j+1) = colsil(j+1)
   colsil(3*j+2) = colsil(j+1) + 1
   colsil(3*j+3) = colsil(j+1) + 2
 END DO
 ncol   = ncol * 3
 130 ntpers = 2
 IF(nwds < 54) GO TO 150
 ntpers = 6
 
!     LOCATE POSITION OF MATRIX TERMS FOR DESIRED SIL
 
 150 DO  kcol = 1,ncol
   IF (colsil(kcol) == isil) GO TO 170
 END DO
 icode = 2
 GO TO 8010
 
 170 iloc = (kcol-1)*nrow*ntpers + 1
 
!     EXTRACT MATRIX TERMS AND STORE THEM IN CORE
 
 icode = 3
 jcore = z(isil)
 IF (jcore == 0) GO TO 8010
 kcore = z(jcore)
 DO  i = 1,nrow
   z(kcore) = rowsil(i)
   IF (ntpers == 2) z(kcore) = -rowsil(i)
   kcore = kcore + 1
   DO  j = 1,ntpers
     z(kcore) = terms(iloc)
     iloc  = iloc  + 1
     kcore = kcore + 1
   END DO
 END DO
 z(jcore) = kcore
 
 GO TO 100
 
!     END OF FILE ON XDICT - PREPARE TO PACK OUT COLUMNS IN CORE
 
 200 CALL CLOSE (xdict,optc)
 CALL CLOSE (xmat,optc)
 CALL gopen (outmat,z(ibuf1),optw)
 
!     PACK OUT COLUMNS
 
 DO  i = ifcol,ilcol
   CALL bldpk (typin,typout,outmat,0,0)
   IF (z(i) == 0) GO TO 210
   
   iloc = z(i) + 1
   nloc = z(iloc-1) - iloc
   CALL pakcol (z(iloc),nloc)
   
   210 CALL bldpkn (outmat,0,mcb)
 END DO
 
 CALL CLOSE (outmat,optc)
 
!     RETURN FOR ADDITIONAL PASS IF MORE NONZERO COLUMNS REMAIN
 
 optw = wrt
 IF (ilcol < luset) GO TO 30
 
!     ALL COLUMNS PROCESSED - WRITE TRAILER AND RETURN
 
 CALL wrttrl (mcb)
 RETURN
 
!     ERROR CONDITIONS
 
 1001 n = -1
 GO TO 1100
 1002 n = -2
 GO TO 1100
 1003 n = -3
 GO TO 1100
 1008 n = -8
 GO TO 1100
 
 1100 CALL mesage (n,FILE,NAME)
 
 8010 WRITE  (nout,9010) sfm,icode
 9010 FORMAT (a25,' 8010, LOGIC ERROR IN SUBROUTINE FLBEMA - CODE',i3/)
 n = -61
 GO TO 1100
END SUBROUTINE flbema
