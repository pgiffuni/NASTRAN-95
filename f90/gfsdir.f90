SUBROUTINE gfsdir
     
!     THIS ROUTINE PERFORMS THE DIRECT FORMULATION OF THE
!     FLUID/STRUCTURE MATRICES
 
 EXTERNAL      andf
 INTEGER :: scr1     ,scr2     ,scr3     ,scr4     ,scr5     ,  &
     scr6     ,scr7     ,scr8     ,axy      ,afry     ,  &
     kyy      ,dkaa     ,dkfrfr   ,kaa      ,maa      ,  &
     gm       ,GO       ,usets    ,usetf    ,kmat     ,  &
     mmat     ,gia      ,pvec     ,ident    ,kjjl     ,  &
     anybar   ,afy      ,awy      ,scr9     ,kaabar   ,  &
     amy      ,aaybar   ,awj      ,kjj      ,gjw      ,  &
     any      ,aoy      ,ajw      ,ac       ,gyw      ,  &
     aay      ,kwwbar   ,mwwbar   ,kc       ,h        ,  &
     uset     ,mcb(7)   ,mbit     ,sbit     ,sfbit    ,  &
     obit     ,um       ,us       ,uo       ,ug       ,  &
     un       ,ua       ,uf       ,uy       ,uab      ,  &
     ufr      ,ui       ,z        ,sysbuf   ,two      ,  &
     FILE     ,typin    ,typout   ,badd(11) ,NAME(2)  ,  &
     ur       ,usg      ,usb      ,ul       ,ux       ,  &
     uz       ,bit      ,ayw      ,mt       ,hc       , comptp   ,andf
 REAL :: rz(1)    ,kcomp    ,rbadd(12)
 DOUBLE PRECISION :: dbadd(5)
 
!     MODULE PARAMETERS
 
 COMMON / BLANK  /       nograv   ,nofree   ,kcomp   ,comptp    ,  &
     FORM     ,lmodes
 
!     SYSTEM COMMON
 
 COMMON / system /       sysbuf
 
!     CALCV COMMON BLOCK
 
 COMMON / patx   /       lcore    ,nsub0    ,nsub1    ,nsub2    , uset
 
!     OPEN CORE
 
 COMMON / zzzzzz /       z(1)
 
!     PACK COMMON BLOCKS
 
 COMMON / zblpkx /       a(4)     ,irow
 COMMON / packx  /       typin    ,typout   ,ii       ,nn       , incr
 
!     POWERS OF TWO
 
 COMMON / two    /       two(32)
 
!     USET BIT POSITIONS
 
 COMMON / bitpos /       um       ,uo       ,ur       ,usg      ,  &
     usb      ,ul       ,ua       ,uf       ,  &
     us       ,un       ,ug       ,bit(12)  ,  &
     ux       ,uy       ,ufr      ,uz       , uab      ,ui
 
!     SCRATCH FILE ASSIGNMENTS
 
 EQUIVALENCE   (badd(1),rbadd(2)) , (dbadd(1),rbadd(3)) , (rz(1), z(1)),  &
     (scr1 , pvec , ident , kjjl) , (scr2 , anybar , afy , awy ) ,  &
     (scr3 , amy , aaybar , awj , gjw , kjj) , (scr4 , aoy , ajw , gyw)  ,  &
     (scr5 , aay , kwwbar , ayw ) , (scr6 , ac) ,  &
     (scr7 , kc  , mt) , (scr8 , h)  ,  &
     (scr9 , gia)
 
!     GINO FILE ASSIGNMENTS
 
 DATA          axy      ,afry     ,kyy      ,dkaa     ,dkfrfr   ,  &
     kaa      ,maa      ,gm       ,GO       ,usets    ,  &
     usetf    ,kmat     ,mmat     ,hc       ,  &
     gia      ,scr1     ,scr2     ,scr3     ,scr4     ,  &
     scr5     ,scr6     ,scr7     ,scr8               /  &
     101      ,102      ,103      ,104      ,105      ,  &
     106      ,107      ,108      ,109      ,110      ,  &
     111      ,201      ,202      ,205      ,  &
     203      ,301      ,302      ,303      ,304      ,  &
     305      ,306      ,307      ,308                /
 
 DATA   NAME / 4HGFSD   ,4HIR     /
 DATA   badd / 11*0     /
 
 
 any    = scr4
 kaabar = scr2
 mwwbar = scr6
 
 lcore = korsz(z(1))
 ibuf  = lcore - sysbuf - 1
 IF (ibuf < 0) GO TO 1008
 
!     REDUCE FLUID / STRUCTURE AREA MATRIX.  MATRIX IS TREATED AS
!     A LOAD VECTOR
 
 mcb(1) = usets
 CALL rdtrl (mcb)
 mbit = andf(mcb(5),two(um))
 sbit = andf(mcb(5),two(us))
 obit = andf(mcb(5),two(uo))
 
 uset = usets
 
!     PARTITION OUT MULTIPOINT CONSTRAINTS
 
 IF (mbit == 0) THEN
   GO TO    20
 END IF
 10 CALL calcv  (pvec,ug,un,um,z(1))
 CALL gfsptn (axy,anybar,amy,0,0,0,pvec)
 CALL ssg2b  (gm,amy,anybar,any,1,2,1,scr1)
 GO TO 30
 
 20 any = axy
 
!     PARTITION OUT SINGLE POINT CONSTRAINTS
 
 30 IF (sbit == 0.0) THEN
   GO TO    50
 END IF
 40 CALL calcv  (pvec,un,uf,us,z(1))
 CALL gfsptn (any,afy,0,0,0,0,pvec)
 GO TO 60
 
 50 CALL gfswch (afy,any)
 
!     PARTITION OUT OMITS
 
 60 IF (obit == 0.0) THEN
   GO TO    80
 END IF
 70 CALL calcv  (pvec,uf,ua,uo,z(1))
 CALL gfsptn (afy,aaybar,aoy,0,0,0,pvec)
 CALL ssg2b  (GO,aoy,aaybar,aay,1,2,1,scr1)
 GO TO 85
 
 80 CALL gfswch (aay,afy)
 
!     IF FREE SURFACE POINTS EXIST - MERGE THEM WITH THE REDUCED
!     AREA MATRIX
 
 85 uset = usetf
 IF (nofree < 0) THEN
   GO TO   100
 END IF
 90 CALL calcv  (pvec,ua,uab,ufr,z(1))
 CALL gfsmrg (awy,aay,afry,0,0,0,pvec)
 GO TO 110
 
 100 CALL gfswch (awy,aay)
 
!     DETERMINE IF ANY SINGLE POINT CONSTRAINTS EXIST ON THE FLUID
 
 110 CALL calcv (pvec,uy,uf,us,z(1))
 nuy   = nsub0 + nsub1
 sfbit = 1
 IF (nsub1 == 0) sfbit = 0
 
!     IF SPC POINTS EXIST ON THE FLUID - PARTITION THEM OUT OF
!     THE FLUID AREA AND STIFFNESS MATRIX
 
 IF (sfbit == 0.0) THEN
   GO TO   130
 END IF
 120 CALL gfsptn (awy,awj,0,0,0,pvec,0)
 CALL gfstrn (awj,ajw,scr2,scr5)
 CALL gfsptn (kyy,kjj,0,0,0,pvec,pvec)
 GO TO 170
 
!     NO SPC POINTS EXIST ON THE FLUID
 
!     CONSTRAIN THE FIRST FLUID POINT TO REMOVE ANY POTENTIAL
!     SINGULARITIES
 
 130 CALL gfsspc (nuy,pvec)
 nsub0 = nuy - 1
 nsub1 = 1
 CALL gfsptn (kyy,kjj,0,0,0,pvec,pvec)
 
!     GENERATE THE H TRANSFORMATION MATRIX
 
 CALL gfsh   (nuy,h)
 CALL gfstrn (awy,ayw,scr1,scr6)
 CALL ssg2b  (h,ayw,0,ajw,0,2,1,scr6)
 
!     CHECK COMPRESSIBLITY TYPE
 
 IF (comptp > 0) GO TO 140
 
!     A SPRING WILL BE GENERATED TO COUPLE THE STRUCTURE AND THE
!     FREE SURFACE TO RESTRICT VOLUME CHANGES
 
!     COMPUTE THE COMPRESSIBLITY MATRIX WHICH CONTAINS THIS SPRING
 
 CALL gfscom (awy,nuy,kc,ident,ac,scr5)
 GO TO 170
 
!     PURELY INCOMPRESSIBLE APPROACH - A CONSTRAINT EQUATION IS
!     GENERATED TO RESTRICT VOLUME CHANGE
 
!     GENERATE HC MATRIX WHICH CONTAINS THE CONSTRAINT
 
 140 CALL gfshc (awy,nuy,hc,ident,ac,mrow)
 
!     SOLVE FOR THE INITIAL PRESSURE TRANSFORMATION MATRIX
 
 170 CALL factor (kjj,kjjl,scr2,scr5,scr6,scr9)
 CALL ssg3a  (0,kjjl,ajw,gjw,scr5,scr6,-1,0)
 
!     IF GRAVITY EXISTS - ADD THE ADDITIONAL STIFFNESS
 
 IF (nograv < 0) THEN
   GO TO   190
 END IF
 180 badd (1) = 2
 dbadd(1) = 1.0D0
 badd (7) = 2
 dbadd(4) = 1.0D0
 CALL ssg2c (kaa,dkaa,kaabar,0,badd)
 GO TO 200
 
 190 kaabar = kaa
 
!     IF FREE SURFACE EXISTS - MERGE THE STIFFNESS TO SOLUTION SIZE
!     AND EXPAND THE MASS MATRIX
 
 200 IF (nofree < 0) THEN
   GO TO   220
 END IF
 210 CALL calcv  (pvec,ua,uab,ufr,z(1))
 CALL gfsmrg (kwwbar,kaabar,0,0,dkfrfr,pvec,pvec)
 CALL gfsmrg (mwwbar,maa,0,0,0,pvec,pvec)
 GO TO 230
 
 220 CALL gfswch (kwwbar,kaabar)
 mwwbar = maa
 
!     COMPUTE THE FINAL MASS MATRIX
!     FOR COMPTP = 1 THIS MATRIX IS NOT THE FINAL ONE
 
 230 CALL ssg2b (ajw,gjw,mwwbar,mmat,1,2,1,scr2)
 
!     COMPUTE THE FINAL STIFFNESS MATRIX
 
 IF (sfbit == 0.0) THEN
   GO TO   240
 ELSE
   GO TO   260
 END IF
 240 IF (comptp > 0) GO TO 250
 
!     ADD IN THE SPRING FACTOR KC
 
 badd (1) = 2
 dbadd(1) = 1.0D0
 badd (7) = 2
 dbadd(4) = 1.0D0
 CALL ssg2c (kwwbar,kc,kmat,0,badd)
 GO TO 270
 
!     APPLY THE CONSTRAINT EQUATION TO STIFFNESS AND MASS FOR
!     THE INCOMPRESSIBLE APPROACH
 
 250 CALL ssg2b (hc,kwwbar,0,scr2,1,2,1,scr1)
 CALL ssg2b (scr2,hc,0,kmat,0,2,1,scr1)
 CALL ssg2b (hc,mmat,0,scr2,1,2,1,scr1)
 CALL ssg2b (scr2,hc,0,mt,0,2,1,scr1)
 
!     ADD 1.0 TO THE NULL COLUMN IN THE MASS MATRIX TO PREVENT
!     SINGULATITIES
 
 CALL gfsmt (mt,mmat,mrow)
 GO TO 270
 
 260 CALL gfswch (kmat,kwwbar)
 
!     TRANSFORM THE FINAL PRESSURE TRANSFORMATION MATRIX OR IF
!     SPC POINTS EXIST ON THE FLUID MERGE IN ZEROS
 
 270 IF (sfbit == 0.0) THEN
   GO TO   280
 ELSE
   GO TO   300
 END IF
 280 CALL ssg2b (h,gjw,0,gyw,1,2,1,scr5)
 GO TO 310
 
 300 CALL calcv  (pvec,uy,uf,us,z(1))
 CALL gfsmrg (gyw,gjw,0,0,0,0,pvec)
 
!     PARTITON OUT THE FREE SURFACE POINTS
 
 310 IF (nofree < 0) THEN
   GO TO   330
 END IF
 320 CALL calcv  (pvec,uy,ufr,ui,z(1))
 CALL gfsptn (gyw,0,gia,0,0,0,pvec)
 RETURN
 
 330 CALL gfswch (gia,gyw)
 RETURN
 
!     ERROR CONDITIONS
 
 1008 n = -8
 CALL mesage (n,FILE,NAME)
 RETURN
END SUBROUTINE gfsdir
