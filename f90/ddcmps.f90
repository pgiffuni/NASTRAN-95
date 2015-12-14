SUBROUTINE ddcmps
     
!     DDCMPS IS THE DMAP DRIVER FOR SDCMPS
 
!     SDCMPS   USET,GPL,SIL,KAA/LLL,ULL/SYM=0/DIAGCK=0/DIAGET=20/
!              PDEFCK=0/SING=0/SET=L/CHLSKY=0/DET=0.0D0/MINDIA=0.0D0/
!              POWER=0/SUBNAM=NONE
 
!     SYM      =  1 - USE SYMMETRIC DECOMPOSITION
!                 0 - CHOOSE WHICH DECOMPOSITION BASED ON INPUT MATRIX
!                -1 - USE UNSYMETRIC DECOMPOSITION
!     DIAGCK   =  DIAGONAL SINGULARITY CHECK              (SDCMPS)
!                 - = NO CHECK                            (SDCMPS)
!                 0 = NONFATAL                            (SDCMPS)
!                 + = MAX ALLOWED FATAL                   (SDCMPS)
!     DIAGET   =  DIAGONAL SINGULARITY ERROR TOLERANCE.   (SDCMPS)
!     PDEFCK   =  POSITIVE DEFINATE CHECK                 (SDCMPS)
!                 - = NO CHECK                            (SDCMPS)
!                 0 = NONFATAL                            (SDCMPS)
!                 + = MAX ALLOWED FATAL                   (SDCMPS)
!     SING     =  SINGULARITY OUTPUT FLAG
!                 1 = OK
!                 0 = NONCONSERVATIVE OR ES FAILURE
!                -1 = SINGULAR OR LIMITS EXCEEDED
!     SET      =  SET MATRIX BELONGS TO                   (SDCMPS)
!     CHLSKY   =  1 USE CHOLESKY DECOMPOSITION LLL = C
!     DET      =  DETERMINANT OF KAA
!     MINDIA   =  MINIMUM DIAGONAL OF ULL
!     POWER    =  SCALE FACTOR FOR DET
!     SUBNAM   =  SUBSTRUCTURE NAME                       (SDCMPS)
 
 LOGICAL :: opnscr   ,first
 INTEGER :: buf6     ,chlsky   ,diagck   ,diaget   ,NAME(2) ,  &
     nam(2)   ,outpt    ,parm     ,pdefck   ,power   ,  &
     rect     ,set      ,sing     ,sqr      ,sym     , ull
 REAL :: zz(1)    ,zzz(1)   ,zzzz(1)  ,zm(1)
 DOUBLE PRECISION :: cdet     ,cmndia   ,mindia   ,sdetc    ,minds   ,  &
     ddet     ,dmndia   ,sdet
 CHARACTER (LEN=27) :: swm
 CHARACTER (LEN=25) :: sfm
 CHARACTER (LEN=29) :: uim
 CHARACTER (LEN=25) :: uwm
 CHARACTER (LEN=23) :: ufm
 COMMON /xmssg /  ufm      ,uwm      ,uim      ,sfm      ,swm
 COMMON /BLANK /  isym     ,diagck   ,diaget   ,pdefck   ,sing    ,  &
     set(2)   ,chlsky   ,det(2)   ,mindia   ,power   , subnam(2)
 COMMON /sdcq  /  nerr(2)   ,noglev   ,buf6     ,iscmsg   ,iscdia ,  &
     istscr    ,kpdfck   ,kdgck    ,kdget    ,kprec  , parm(4)   ,opnscr   ,first
 COMMON /sfact /  ifila(7)  ,ifill(7) ,ifilu(7) ,kscr1    ,  &
     kscr2     ,nz       ,sdet     ,sdetc    ,kpow   , kscr3     ,minds    ,ichlk
 COMMON /dcompx/  ia(7)     ,il(7)    ,iu(7)    ,iscr1    ,  &
     iscr2     ,iscr3    ,ddet     ,ipow     , nzz       ,dmndia   ,ib
 COMMON /cdcmpx/  ja(7)     ,jl(7)    ,ju(7)    ,jscr1    ,  &
     jscr2     ,jscr3    ,cdet(2)  ,jpow     , nzzz      ,cmndia   ,jb
 COMMON /names /  knames(19)
 COMMON /system/  ksystm(69)
 COMMON /zzzzzz/  z(1)
 EQUIVALENCE      (zz(1),z(1))
 EQUIVALENCE      (zzz(1),z(1))
 EQUIVALENCE      (zzzz(1),z(1))
 EQUIVALENCE      (zm(1),z(1))
 EQUIVALENCE      (ksystm( 1),nbufsz) ,(ksystm( 2),outpt) ,  &
     (knames(12),sqr   ) ,(knames(13),rect ) , (knames(17),sym   )
 DATA    luset ,  lgpl ,lsil ,kaa  ,lll  ,ull  ,lscr1,lscr2,lscr3 /  &
     101   ,  102  ,103  ,104  ,201  ,202  ,301  ,302  ,303   /
 DATA    lscr4 ,  lscr5,lscr6/ 304   ,  305  ,306  /
 DATA    NAME  /  4HDDCM, 4HPS   /
 DATA    nam   /  4HSDCM, 4HPS   /
 
!     NOTE SYM DECOMP DOES NOT OUTPUT  ULL
 
 
 opnscr = .false.
 first  = .true.
 sing   = 1
 ja(1)  = kaa
 CALL rdtrl (ja)
 IF (ja(1) < 0) GO TO 490
 iform  = ja(4)
 IF (isym < 0) THEN
   GO TO    10
 ELSE IF (isym == 0) THEN
   GO TO    50
 ELSE
   GO TO    30
 END IF
 10 IF (iform /= sym) GO TO 20
 CALL page2 (2)
 WRITE  (outpt,15) swm,nam
 15 FORMAT (a27,' 2340, MODULE ',2A4,' HAS BEEN REQUESTED TO DO ',  &
       'UNSYMETRIC DECOMPOSITION OF A SYMETRIC MATRIX')
   20 iform = rect
   IF (ja(2) == ja(3)) iform = sqr
   GO TO 50
   
   30 IF (iform == sym) GO TO 50
   CALL page2 (2)
   WRITE  (outpt,40) swm,nam
   40 FORMAT (a27,' 2341, MODULE ',2A4,' HAS BEEN FURNISHED A SQUARE ',  &
       'MATRIX MARKED UNSYMETRIC FOR SYMETRIC DECOMPOSITION.')
   iform = sym
   50 isym  = -1
   IF (iform == sym) isym = 1
   ja(4) = iform
   i = 0
   IF (ja(2) == ja(3)) GO TO 60
   CALL page2 (2)
   i = 1
   WRITE  (outpt,55) swm,nam
   55 FORMAT (a27,' 2375, MODULE ',2A4,' HAS BEEN REQUESTED TO ',  &
       'DECOMPOSE A RECTANGULAR MATRIX')
   60 CONTINUE
   IF (isym < 0) GO TO 200
   
!     SET UP CALL TO SDCOMP
   
   IF (i /= 0) GO TO 500
   ifila(1) = kaa
   CALL rdtrl (ifila)
   ifill(1) = lll
   ifilu(1) = lscr4
   kscr1 = lscr1
   kscr2 = lscr2
   kscr3 = lscr3
   ifill(5) = ifila(5)
   ichlk = chlsky
   IF (ifila(5) <= 2) GO TO 100
   nz = korsz (z)
   CALL sdcomp (*400,z,z,z)
   GO TO 130
   100 nz = korsz(zzzz)
   iscmsg = lscr5
   iscdia = lscr6
   kpdfck = pdefck
   kdgck  = diagck
   kdget  = diaget
   CALL sdcmps (zzzz,zzzz,zzzz)
   IF (nerr(1)+nerr(2) == 0) GO TO 110
   buf6 = korsz(zm) - 2*nbufsz + 1
   IF (buf6+nbufsz <= 0) GO TO 510
   CALL sdcmm (zm,set(1),ifila(2),ifila(1),luset,lgpl,lsil,subnam)
   sing = 0
   
!     ONLY ES CHECK AND NONCONSERVATIVE COLUMN CAN EXIT WITH SING = 1
!     OR IF USER DESIRES TO CONTINUE
   
   IF (noglev > 0) sing = -1
   110 CONTINUE
   IF (parm(1) /= 0) CALL mesage (parm(1),parm(2),parm(3))
   130 det(1) = sdet
   det(2) = sdetc
   mindia = minds
   power  = kpow
   ifill(2) = ifila(2)
   ifill(3) = ifila(3)
   ifill(4) = 4
   IF (sing >= 0) CALL wrttrl (ifill)
   GO TO 410
   
!     SET UP CALL TO DECOMP
   
   200 CONTINUE
   IF (ja(5) > 2) GO TO 300
   ia(1) = kaa
   CALL rdtrl (ia)
   il(1) = lll
   iu(1) = ull
   nzz   = korsz(zz)
   iscr1 = lscr1
   iscr2 = lscr2
   iscr3 = lscr3
   ib    = 0
   il(5) = 2
   CALL decomp (*400,zz,zz,zz)
   iu(5) = 2
   il(4) = 4
   iu(4) = 5
   il(3) = il(2)
   iu(3) = iu(2)
   det(1)= ddet
   det(2)= 0.0
   power = ipow
   mindia= dmndia
   CALL wrttrl (iu)
   CALL wrttrl (il)
   GO TO 410
   
!     SET UP CALL TO CDCOMP
   
   300 CONTINUE
   jl(1) = lll
   ju(1) = ull
   jscr1 = lscr1
   jscr2 = lscr2
   jscr3 = lscr3
   nzzz  = korsz(zzz)
   jl(5) = 4
   jb    = 0
   CALL cdcomp (*400,zzz,zzz,zzz)
   ju(5) = 4
   jl(4) = 4
   ju(4) = 5
   jl(3) = jl(2)
   ju(3) = ju(2)
   det(1)= cdet(1)
   det(2)= cdet(2)
   mindia= cmndia
   power = jpow
   CALL wrttrl (jl)
   CALL wrttrl (ju)
   GO TO 410
   
   400 sing   = -1
   det(1) = 0.0
   det(2) = 0.0
   power  = 0
   mindia = 0.0
   410 RETURN
   
!     ERROR  MESSAGES
   
!     PURGED INPUT
   
   490 parm(1) = -1
   parm(2) = kaa
   GO TO 520
   
!     NUMBER ROWS.NE.COLUMNS
   
   500 parm(1) = -16
   parm(2) = kaa
   GO TO 520
   
!     INSUFFICIENT CORE
   
   510 parm(1) = -8
   parm(2) = -buf6 - nbufsz
   520 parm(3) = NAME(1)
   parm(4) = NAME(2)
   GO TO 110
 END SUBROUTINE ddcmps
