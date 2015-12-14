SUBROUTINE partn2 (cp,rp,core,buf)
     
!     THIS IS AN INITIALIZATION ROUTINE FOR PARTN1 AND MERGE1.
!     IT CALLS PARTN3 TO BUILD THE BIT STRINGS FROM THE PARTITIONING
!     VECTORS -CP- AND -RP- AND SETS DEFAULT OPTIONS FOR -CP- AND  -RP-
!     BASED ON -SYM-.
 
 
 
 INTEGER, INTENT(IN OUT)                  :: cp
 INTEGER, INTENT(IN OUT)                  :: rp
 INTEGER, INTENT(IN OUT)                  :: core
 INTEGER, INTENT(IN OUT)                  :: buf(4)
 LOGICAL :: cpnull  ,rpnull  ,cphere  ,rphere
 INTEGER :: subr(2) , cpsize  ,rpsize  ,cpones  ,rpones  ,z       ,  &
     sym     ,TYPE    ,FORM    ,sysbuf  ,outpt   , cpcol   ,rpcol
 CHARACTER (LEN=27) :: swm
 CHARACTER (LEN=25) :: sfm
 CHARACTER (LEN=29) :: uim
 CHARACTER (LEN=25) :: uwm
 CHARACTER (LEN=23) :: ufm
 COMMON /xmssg / ufm     ,uwm     ,uim     ,sfm     ,swm
 COMMON /system/ sysbuf  ,outpt
 COMMON /prtmrg/ cpsize  ,rpsize  ,cpones  ,rpones  ,cpnull  ,  &
     rpnull  ,cphere  ,rphere  ,icp     ,ncp     , irp     ,nrp
 COMMON /zzzzzz/ z(1)
 COMMON /BLANK / sym     ,TYPE    ,FORM(4) ,cpcol   ,rpcol   , ireqcl
 DATA    subr  / 4HPART  ,4HN2    /
 
 
!              I             I             I                I
!       SYM    I  RP-PURGED  I  CP-PURGED  I NEITHER-PURGED I
!     ---------+-------------+-------------+----------------+--------
!              I             I             I                I
!       .LT.0  I  RP IS SET  I  CP IS SET  I  RP MUST HAVE  I
!              I  = TO CP    I  = TO RP    I  SAME ONES-S   I
!              I             I             I  COUNT AS CP   I
!     ---------+-------------+-------------+----------------+--------
!              I             I             I                I
!       .GE.0  I  RP IS SET  I  CP IS SET  I  USE CP AND RP I
!              I  TO ALL 0   I  TO ALL 0   I                I
!              I             I             I                I
 
!     IN ALL CASES, RESULTANT -CP- AND -RP- DIMENSIONS MUST BE
!     COMPATIBLE TO THOSE OF  -A-
 
 
!     CONVERT COLUMN PARTITIONING VECTOR TO BIT STRING.
 
 icp = 1
 ireqcl = cpcol
 CALL partn3 (cp,cpsize,cpones,icp,ncp,cphere,buf,core)
 cpcol = ireqcl
 IF (cphere) GO TO 10
 irp = 1
 GO TO 20
 10 irp = ncp + 1
 
!     CONVERT ROW PARTITIONING VECTOR TO BIT STRING.
 
 20 ireqcl = rpcol
 CALL partn3 (rp,rpsize,rpones,irp,nrp,rphere,buf,core)
 rpcol = ireqcl
 
!     BRANCH ON SYMMETRIC OR  NON-SYMMETRIC DMAP VARIABLE SYM
 
 cpnull = .false.
 rpnull = .false.
 IF (sym < 0.0) THEN
   GO TO    30
 ELSE
   GO TO   140
 END IF
 
!     DMAP USER CLAIMS SYMMETRIC INPUT AND OUTPUT
 
 30 IF (cphere) GO TO 70
 
!     -CP- IS PURGED.  CHECK FOR -RP- PURGED (ERROR), AND IF NOT SET
!     -CP- BITS EQUAL TO -RP- BITS
 
 IF (rphere) GO TO 60
 
!     BOTH -RP- AND -CP- PURGED AND -A- IS NOT PURGED (ERROR).
 
 40 WRITE  (outpt,50) sfm
 50 FORMAT (a25,' 2170, BOTH THE ROW AND COLUMN PARTITIONING VECTORS',  &
     ' ARE PURGED AND ONLY ONE MAY BE.')
 CALL mesage (-61,0,subr)
 
!     SET CP-ONES = RP-ONES BY SIMPLE EQUIVALENCE OF CORE SPACE
 
 60 icp = irp
 ncp = nrp
 cpones = rpones
 cpsize = rpsize
 GO TO 170
 
!     -CP- IS NOT PURGED.  IF -RP- IS PURGED IT IS SET EQUAL TO -CP-.
 
 70 IF (rphere) GO TO 80
 irp = icp
 nrp = ncp
 rpones = cpones
 rpsize = cpsize
 GO TO 170
 
!     BOTH -RP- AND -CP- ARE PRESENT AND SINCE USER HAS SPECIFIED A
!     SYMMETRIC OUTPUT PARTITION IS DESIRED THE NUMBER OF
!     NON-ZEROS IN-CP- MUST EQUAL THE NUMBER OF NON-ZEROS IN -RP- FOR NO
!     ERROR HERE.
 
 80 IF (cpones == rpones .AND. cpsize == rpsize) GO TO 100
 WRITE  (outpt,90) swm,cp,rp
 90 FORMAT (a27,' 2171, SYM FLAG INDICATES TO THE PARTITION OR MERGE',  &
     ' MODULE THAT A SYMMETRIC MATRIX IS TO BE', /5X,  &
 ' OUTPUT.  THE PARTITIONING VECTORS',2I4,' HOWEVER DO NOT',  &
       ' CONTAIN AN IDENTICAL NUMBER OF ZEROS AND NON-ZEROS.')
   
!     CHECK FOR ORDER OF ONES IN ROW AND COLUMN PARTITIONING VECTOR.
   
   100 IF (cpsize /= rpsize) GO TO 170
   j = irp
   DO  i = icp,ncp
     IF (z(i) /= z(j)) GO TO 120
     j = j + 1
   END DO
   GO TO 170
   
!     ROW AND COLUMN PARTITIONING VECTORS DO NOT HAVE SAME ORDER.
   
   120 WRITE  (outpt,130) swm
   130 FORMAT (a27,' 2172, ROW AND COLUMN PARTITIONING VECTORS DO NOT ',  &
         'HAVE IDENTICAL ORDERING OF ZERO', /5X,' AND NON-ZERO ',  &
         'ELEMENTS, AND SYM FLAG INDICATES THAT A SYMMETRIC ',  &
         'PARTITION OR MERGE IS TO BE PERFORMED.')
     GO TO 170
     
!     DMAP USER DOES NOT REQUIRE SYMMETRY
     
     140 IF (cphere) GO TO 160
     
!     -CP- IS PURGED.  THUS -RP- MUST BE PRESENT FOR NO ERROR.
     
     IF (rphere) GO TO 150
     GO TO 40
     
!     SET CP-ONES EQUAL TO 0 AND CPSIZE = 0
     
     150 cpnull = .true.
     cpsize = 0
     cpones = 0
     GO TO 170
     
!     -CP- NOT PURGED.  IF -RP- IS PURGED SET IT NULL.
     
     160 IF (rphere) GO TO 170
     rpnull = .true.
     nrp    = irp - 1
     rpsize = 0
     rpones = 0
     170 RETURN
   END SUBROUTINE partn2
