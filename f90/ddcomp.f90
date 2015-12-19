SUBROUTINE ddcomp
     
!     DDCOMP IS THE DMAP DRIVER FOR DECOMP
 
!     DECOMP    KAA/LLL,ULL/SYM/CHLSKY/MINDIA/DET/POWER/SING $
 
!        SYM    =  1 - USE SYMMETRIC DECOMPOSITION
!                  0 - CHOOSE WHICH DECOMPOSITION BASED ON INPUT MATRIX
!                 -1 - USE UNSYMETRIC DECOMPOSITION
!        CHLSKY =  1 USE CHOLESKY DECOMPOSITION LLL = C
!        DET    =  DETERMINANT OF KAA
!        POWER  =  SCALE FACTOR FOR DET
!        MINDIA =  MINIMUM DIAGONAL OF ULL
!        SING   = -1 SINGULAR MATRIX
 
 INTEGER :: ull       ,sym      ,power    ,sing     ,  &
     chlsky    ,NAME(2)  ,sqr      ,rect     , outpt     ,upper
 REAL :: zz(1)     ,zzz(1)
 DOUBLE PRECISION :: cdet      ,cmndia   ,mindia   ,sdetc    ,  &
     minds     ,ddet     ,dmndia   ,sdet
 CHARACTER (LEN=27) :: swm
 CHARACTER (LEN=25) :: sfm
 CHARACTER (LEN=29) :: uim
 CHARACTER (LEN=25) :: uwm
 CHARACTER (LEN=23) :: ufm 
 COMMON /xmssg /  ufm       ,uwm      ,uim      ,sfm      , swm
 COMMON /BLANK /  isym      ,chlsky                       ,  &
     mindia    ,det(2)   ,power    ,sing
 COMMON /sfact /  ifila(7)  ,ifill(7) ,ifilu(7) ,kscr1    ,  &
     kscr2     ,nz       ,sdet     ,sdetc    ,  &
     kpow      ,kscr3    ,minds    ,ichlk
 COMMON /dcompx/  ia(7)     ,il(7)    ,iu(7)    ,iscr1    ,  &
     iscr2     ,iscr3    ,ddet     ,ipow     , nzz       ,dmndia   ,ib
 COMMON /cdcmpx/  ja(7)     ,jl(7)    ,ju(7)    ,jscr1    ,  &
     jscr2     ,jscr3    ,cdet(2)  ,jpow     , nzzz      ,cmndia   ,jb
 COMMON /names /  knames(19)
 COMMON /system/  ksystm(65)
 COMMON /zzzzzz/  z(1)
 EQUIVALENCE      (zz(1),z(1))
 EQUIVALENCE      (zzz(1),z(1))
 EQUIVALENCE      (ksystm( 2),outpt)  ,(knames(12),sqr)   ,  &
     (knames(13),rect )  ,(knames(17),sym)   ,  &
     (knames(16),upper)  ,(knames(15),lower)
 DATA    kaa,     lll,   ull,    lscr1,  lscr2,  lscr3,  lscr4 /  &
     101,     201,   202,    301  ,  302  ,  303  ,  304   /
 DATA    NAME  /  4HDDCO,4HMP    /
 
 sing  = 0
 ja(1) = kaa
 CALL rdtrl (ja)
 iform = ja(4)
 IF (isym < 0) THEN
   GO TO    10
 ELSE IF (isym == 0) THEN
   GO TO    50
 ELSE
   GO TO    30
 END IF
 10 IF (iform == sym) WRITE (outpt,20) swm,NAME
 20 FORMAT (a27,' 2340, MODULE ',2A4,' HAS BEEN REQUESTED TO DO ',  &
       'UNSYMMETRIC DECOMPOSITION OF A SYMMETRIC MATRIX')
   iform = rect
   IF (ja(2) == ja(3)) iform = sqr
   GO TO 50
   30 IF (ja(2) == ja(3) .AND. iform /= sym) WRITE (outpt,40) swm,NAME
   40 FORMAT (a27,' 2341, MODULE ',2A4,'HAS BEEN FURNISHED A SQUARE ',  &
       'MATRIX MARKED UNSYMMETRIC FOR SYMMETRIC DECOMPOSITION.')
   iform = sym
   50 isym  = -1
   IF (iform == sym) isym = 1
   ja(4) = iform
   IF (isym < 0) GO TO 200
   
!     SET UP CALL TO SDCOMP
   
   ifila(1) = kaa
   CALL rdtrl (ifila)
   ifill(1) = lll
   ifilu(1) = lscr4
   kscr1    = lscr1
   kscr2    = lscr2
   kscr3    = lscr3
   nz       = korsz(z)
   ifill(5) = ifila(5)
   ichlk    = chlsky
   CALL sdcomp (*400,z,z,z)
   det(1)   = sdet
   det(2)   = sdetc
   mindia   = minds
   power    = kpow
   ifill(2) = ifila(2)
   ifill(3) = ifila(3)
   ifill(4) = lower
   CALL wrttrl (ifill)
   RETURN
   
!     SET UP CALL TO DECOMP
   
   200 CONTINUE
   IF (ja(5) > 2) GO TO 300
   ia(1)  = kaa
   CALL rdtrl (ia)
   il(1)  = lll
   iu(1)  = ull
   nzz    = korsz(zz)
   iscr1  = lscr1
   iscr2  = lscr2
   iscr3  = lscr3
   ib     = 0
   il(5)  = 2
   CALL decomp (*400,zz,zz,zz)
   iu(5)  = 2
   il(4)  = lower
   iu(4)  = upper
   il(3)  = il(2)
   iu(3)  = iu(2)
   det(1) = ddet
   det(2) = 0.0
   power  = ipow
   mindia = dmndia
   CALL wrttrl (iu)
   CALL wrttrl (il)
   RETURN
   
!     SET UP CALL TO CDCOMP
   
   300 CONTINUE
   jl(1)  = lll
   ju(1)  = ull
   jscr1  = lscr1
   jscr2  = lscr2
   jscr3  = lscr3
   nzzz   = korsz(zzz)
   jl(5)  = 4
   jb     = 0
   CALL cdcomp (*400,zzz,zzz,zzz)
   ju(5)  = 4
   jl(4)  = lower
   ju(4)  = upper
   jl(3)  = jl(2)
   ju(3)  = ju(2)
   det(1) = cdet(1)
   det(2) = cdet(2)
   mindia = cmndia
   power  = jpow
   CALL wrttrl (jl)
   CALL wrttrl (ju)
   RETURN
   
   400 sing   = -1
   det(1) = 0.0
   det(2) = 0.0
   power  = 0
   mindia = 0.0
   CALL fname (kaa,ja(1))
   WRITE  (outpt,410) uim,ja(1),ja(2)
   410 FORMAT (a29,' FORM DECOMP MODULE. MATRIX ',2A4,' IS SINGULAR')
   RETURN
 END SUBROUTINE ddcomp
