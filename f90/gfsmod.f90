SUBROUTINE gfsmod
     
!     THIS ROUTINE PERFORMS THE MODAL FORMULATION OF THE
!     FLUID / STRUCTURE MATRICES
 
 EXTERNAL        andf
 INTEGER :: scr1     ,scr2     ,scr3     ,scr4     ,scr5     ,  &
     scr6     ,scr7     ,scr8     ,scr9     ,scr10    ,  &
     axy      ,afry     ,kyy      ,dkaa     ,dkfrfr   ,  &
     usetf    ,phia     ,phix     ,ac       ,pout     ,  &
     lama     ,kmat     ,mmat     ,gia      ,pvec     ,  &
     ident    ,uset     ,usetd    ,kc       ,h        ,  &
     comptp   ,azy      ,ahy      ,ahj      ,ajh      ,  &
     kjj      ,ayh      ,kjjl     ,gjh      ,mzz      ,  &
     kzz      ,mhhbar   ,phiar    ,kzzbar   ,khhbar   ,  &
     gyh      ,mcb(7)   ,sfbit    ,FILE     ,um       ,  &
     uz       ,unz      ,ufr      ,uh       ,uy       ,  &
     uf       ,us       ,ui       ,z        ,sysbuf   ,  &
     two      ,typin    ,typout   ,badd(11) ,NAME(2)  ,  &
     phixr    ,andf     ,namex(2)
 REAL :: rz(1)
 DOUBLE PRECISION :: dbadd(5)
 CHARACTER (LEN=25) :: uwm
 CHARACTER (LEN=23) :: ufm
 COMMON /xmssg / ufm      ,uwm
 COMMON /BLANK / nograv   ,nofree   ,kcomp    ,comptp   ,FORM     , llmode
 COMMON /system/ sysbuf   ,nout
 COMMON /patx  / lcore    ,nsub0    ,nsub1    ,nsub2    ,uset
 COMMON /zzzzzz/ z(1)
 COMMON /packx / typin    ,typout   ,ii       ,nn       ,incr
 COMMON /two   / two(32)
 COMMON /bitpos/ unz      ,uz       ,um       ,uh       ,bit1(3)  ,  &
     uf       ,us       ,bit2(15) ,uy       ,ufr      , bit3(2)  ,ui
 COMMON /gfsmox/ axy      ,afry     ,kyy      ,dkaa     ,dkfrfr   ,  &
     usetf    ,phia     ,phix     ,lama     ,  &
     kmat     ,mmat     ,gia      ,pout     ,  &
     scr1     ,scr2     ,scr3     ,scr4     ,scr5     ,  &
     scr6     ,scr7     ,scr8     ,  &
     lmodes   ,nmodes   ,ibuf     ,sfbit    ,badd     , NAME
 EQUIVALENCE     (badd(2),dbadd(1)) ,(rz(1),z(1)) ,  &
     (scr1,usetd) ,(scr2,pvec,ident,kjjl) ,  &
     (scr3,azy,ahj,kjj,gjh) ,(scr4,ajh,khhbar,gyh) ,  &
     (scr5,ac,ayh,mzz,kzzbar) ,(scr6,kzz) ,  &
     (scr7,kc,ahy) ,(scr8,h) ,(scr9,mmat) , (scr10,gia,mhhbar)
!     DATA            AXY      ,AFRY     ,KYY      ,DKAA     ,DKFRFR   ,
!    1                USETF    ,PHIA     ,PHIX     ,LAMA     ,
!    2                KMAT     ,MMAT     ,GIA      ,POUT     ,
!    3                SCR1     ,SCR2     ,SCR3     ,SCR4     ,SCR5     ,
!    4                SCR6     ,SCR7     ,SCR8     /
!    5                101      ,102      ,103      ,104      ,105      ,
!    6                111      ,112      ,113      ,114      ,
!    7                201      ,202      ,203      ,204      ,
!    8                301      ,302      ,303      ,304      ,305      ,
!    9                306      ,307      ,308      /
 
!     DATA    BADD  / 11*0   /
 DATA    namex / 4HGFSM , 4HOD   /
 
 axy    = 101
 afry   = 102
 kyy    = 103
 dkaa   = 104
 dkfrfr = 105
 usetf  = 111
 phia   = 112
 phix   = 113
 lama   = 114
 kmat   = 201
 mmat   = 202
 gia    = 203
 pout   = 204
 scr1   = 301
 scr2   = 302
 scr3   = 303
 scr4   = 304
 scr5   = 305
 scr6   = 306
 scr7   = 307
 scr8   = 308
 NAME(1)= namex(1)
 NAME(2)= namex(2)
 DO  i = 1,11
   badd(i)= 0
 END DO
 
 
 phiar = scr4
 phixr = scr2
 
 lcore = korsz(z(1))
 ibuf  = lcore - sysbuf - 1
 IF (ibuf < 0) GO TO 1008
 
!     CREATE A DUMMY USET VECTOR TOR USE WITH THE MODAL DISPLACEMENTS
 
!     BIT POSITIONS WILL BE
 
!     UM  - MODAL POINT  UZ + UNZ
!     UZ  - DESIRED MODAL POINT
!     UNZ - MODAL POINT TO BE SKIPPED
!     UFR - FREE SURFACE POINT
!     UH  - UFR + UZ
 
!     SET MODAL DISPLACEMENTS
 
 FILE   = phix
 mcb(1) = phix
 CALL rdtrl (mcb)
 IF (mcb(1) < 0) GO TO 1001
 nmodes = mcb(2)
 IF (llmode > nmodes .OR. llmode == 0) llmode = -1
 lmodes = llmode
 IF (lmodes <= 0) lmodes = nmodes
 IF (lmodes <= 0 .OR. lmodes > nmodes) lmodes = nmodes
 izm  = two(uz)  + two(um) + two(uh)
 inzm = two(unz) + two(um)
 IF (ibuf <= nmodes) GO TO 1008
 DO  i = 1,nmodes
   z(i) = izm
   IF (i > lmodes) z(i) = inzm
 END DO
 
!     SET FREE SURFACE DISPLACEMENTS
 
 ifr  = two(ufr) + two(uh)
 lvec = nmodes
 IF (nofree < 0) THEN
   GO TO    45
 END IF
 20 CALL gopen (usetf,z(ibuf),0)
 30 CALL READ (*40,*40,usetf,ibit,1,0,n)
 IF (andf(ibit,two(ufr)) == 0) GO TO 30
 lvec = lvec + 1
 IF (lvec >= ibuf) GO TO 1008
 z(lvec) = ifr
 GO TO 30
 
 40 CALL CLOSE (usetf,1)
 
!     WRITE DUMMY USETD FILE
 
 45 CALL gopen (usetd,z(ibuf),1)
 CALL WRITE (usetd,z(1),lvec,1)
 CALL CLOSE (usetd,1)
 mcb(1) = usetd
 mcb(2) = lvec
 DO  i = 3,7
   mcb(i) = 0
 END DO
 CALL wrttrl (mcb)
 
!     EXTRACT THE DESIRED MODES FORM THE PHIX MATRIX
 
 uset = usetd
 IF (lmodes >= nmodes) GO TO 70
 CALL calcv (pvec,um,uz,unz,z(1))
 CALL gfsptn (phix,phixr,0,0,0,pvec,0)
 GO TO 80
 
 70 phixr = phix
 
!     TRANSFORM THE FLUID STRUCTURE AREA MATRIX
 
 80 CALL ssg2b (phixr,axy,0,azy,1,2,1,scr5)
 
!     IF FREE SURFACE POINTS EXIST - MERGE THEM WITH THE TRANSFORMED
!     AREA MATRIX
 
 IF (nofree < 0) THEN
   GO TO   100
 END IF
 90 CALL calcv (pvec,uh,uz,ufr,z(1))
 CALL gfsmrg (ahy,azy,afry,0,0,0,pvec)
 GO TO 110
 
 100 CALL gfswch (ahy,azy)
 
!     DETERMINE IF ANY SINGLE POINT CONSTRAINTS EXIST ON THE FLUID
 
 110 uset = usetf
 CALL calcv (pvec,uy,uf,us,z(1))
 nuy = nsub0 + nsub1
 sfbit = 1
 IF (nsub1 == 0) sfbit = 0
 
!     IF SPC POINTS EXIST ON THE FLUID - PARTITION THEM OUT OF THE
!     FLUID AREA AND STIFFNESS MATRICES
 
 IF (sfbit == 0.0) THEN
   GO TO   130
 END IF
 120 CALL gfsptn (ahy,ahj,0,0,0,pvec,0)
 CALL gfstrn (ahj,ajh,scr5,scr6)
 CALL gfsptn (kyy,kjj,0,0,0,pvec,pvec)
 GO TO 160
 
!     IF NO SPC POINTS EXIST ON THE FLUID, CONSTRAIN THE FIRST FLUID
!     POINT TO REMOVE POTENTIAL SINGULARITIES
 
 130 IF (comptp > 0) WRITE (nout,140) uwm
 140 FORMAT (a25,' 8015. THE PURELY INCOMPRESSIBLE METHOD IS AVAIL',  &
     'ABLE ONLY WITH THE DIRECT FORMULATION.')
 CALL gfsspc (nuy,pvec)
 nsub0 = nuy - 1
 nsub1 = 1
 CALL gfsptn (kyy,kjj,0,0,0,pvec,pvec)
 
!     GENERATE THE H TRANSFORMATION MATRIX
 
 CALL gfsh (nuy,h)
 CALL gfstrn (ahy,ayh,scr2,scr6)
 CALL ssg2b (h,ayh,0,ajh,0,2,1,scr6)
 
!     GENERATE THE COMPRESSIBLITY MATRIX
 
 CALL gfscom (ahy,nuy,kc,ident,ac,scr6)
 
!     SOLVE FOR THE INITIAL PRESSURE TRANSFORMATION MATRIX
 
 160 CALL factor (kjj,kjjl,scr5,scr6,scr9,scr10)
 CALL ssg3a (0,kjjl,ajh,gjh,scr5,scr6,-1,0)
 
!     FOR COMPUTER CORE CONSERVATION REASON, THE REST OF GFSMOD IS
!     MOVED TO GFSMO2, WHICH CAN BE SEGMENTED IN PARALLEL WITH GFSMOD.
 
 RETURN
 
!     ERROR EXITS
 
 1001 n = -1
 GO TO 9999
 1008 n = -8
 
 9999 CALL mesage (n,FILE,NAME)
 RETURN
END SUBROUTINE gfsmod
