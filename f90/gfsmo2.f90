SUBROUTINE gfsmo2
     
!     THIS ROUTINE IS THE CONTINUATION OF GFSMOD
 
 REAL :: rz(1)     ,eigval(7)
 
 DOUBLE PRECISION :: dbadd(5)
 
 INTEGER :: scr1     ,scr2     ,scr3     ,scr4     ,scr5  &
     ,scr6     ,scr7     ,scr8     ,scr9     ,scr10  &
     ,axy      ,afry     ,kyy      ,dkaa     ,dkfrfr  &
     ,usetf    ,phia     ,phix     ,ac       ,pout  &
     ,lama     ,kmat     ,mmat     ,gia      ,pvec  &
     ,ident    ,uset     ,usetd    ,kc       ,h  &
     ,comptp   ,azy      ,ahy      ,ahj      ,ajh  &
     ,kjj      ,ayh      ,kjjl     ,gjh      ,mzz  &
     ,kzz      ,mhhbar   ,phiar    ,kzzbar   ,khhbar  &
     ,gyh      ,mcb(7)   ,sfbit    ,FILE     ,um  &
     ,uz       ,unz      ,ufr      ,uh       ,uy  &
     ,uf       ,us       ,ui       ,z        ,sysbuf  &
     ,two      ,typin    ,typout   ,badd(11) ,NAME(2)
 
!     MODULE PARAMETERS
 
 COMMON /BLANK/          nograv   ,nofree   ,kcomp   ,comptp ,FORM     ,llmode
 
!     SYSTEM COMMON
 
 COMMON / system /       sysbuf   ,nout
 
!     CALCV COMMON BLOCK
 
 COMMON / patx /         lcore    ,nsub0    ,nsub1    ,nsub2 ,uset
 
!     OPEN CORE
 
 COMMON / zzzzzz /        z(1)
 
!     PACK COMMON BLOCK
 
 COMMON / packx /        typin    ,typout   ,ii       ,nn ,incr
 
!     POWERS OF TWO
 
 COMMON / two /          two(32)
 
!     USET BIT POSITIONS - SOME OF THESE ARE USED JUST HERE
 
 COMMON / bitpos /       unz      ,uz       ,um       ,uh  &
     ,bit1(3)  ,uf       ,us       ,bit2(15) ,uy       ,ufr      ,bit3(2)  ,ui
 
!     LOCAL VARIABLES FOR GFSMO1 AND GFSMO2
 
 COMMON /gfsmox/axy     ,afry     ,kyy      ,dkaa     ,dkfrfr  &
     ,usetf    ,phia     ,phix     ,lama ,kmat     ,mmat     ,gia      ,pout  &
     ,scr1     ,scr2     ,scr3     ,scr4     ,scr5 ,scr6     ,scr7     ,scr8  &
     ,lmodes   ,nmodes   ,ibuf     ,sfbit    ,badd ,NAME
 
!     SCRATCH FILE ASSIGNMENTS
 
 EQUIVALENCE   ( badd(2) , dbadd(1) ) ,( rz(1) , z(1) )  &
     ,( scr1 , usetd ) ,( scr2 , pvec , ident , kjjl )  &
     ,( scr3 , azy , ahj , kjj , gjh ) ,( scr4 , ajh , khhbar , gyh )  &
     ,( scr5 , ac , ayh , mzz , kzzbar ) ,( scr6 , kzz )  &
     ,( scr7 , kc , ahy ) ,( scr8 , h )  &
     ,( scr9 , mmat ) ,( scr10 , gia , mhhbar )
 
 
!***********************************************************************
 
 
!     GET THE GENERALIZED STIFFNESS AND MASS FOR THE DESIRED MODES
!     FROM THE LAMA DATA BLOCK
 
 IF(2*lmodes >= ibuf) GO TO 1008
 CALL gopen(lama,z(ibuf),0)
 FILE = lama
 CALL fwdrec(*1002,lama)
 igk = 1
 igm = lmodes + 1
 DO  i=1,lmodes
   165 CALL READ(*1001,*1002,lama,eigval,7,0,n)
   IF(eigval(6) == 0.0) GO TO 165
   rz(igk) = eigval(7)
   igk = igk + 1
   rz(igm) = eigval(6)
   igm = igm + 1
 END DO
 CALL CLOSE(lama,1)
 
!     GENERATE THE DIAGONAL MODAL STIFFNESS MATRIX
 
 i1 = 1
 i2 = lmodes
 CALL makmcb(mcb,kzz,lmodes,6,2)
 CALL gopen(kzz,z(ibuf),1)
 typin = 1
 typout = 2
 incr = 1
 DO  i=i1,i2
   ii = i
   nn = i
   CALL pack(rz(i),kzz,mcb)
 END DO
 CALL CLOSE(kzz,1)
 CALL wrttrl(mcb)
 
!     GENERATE THE DIAGANOL MODAL MASS MATRIX
 
 i1 = lmodes + 1
 i2 = 2 * lmodes
 CALL makmcb(mcb,mzz,lmodes,6,2)
 CALL gopen(mzz,z(ibuf),1)
 DO  i=i1,i2
   ii = i - lmodes
   nn = ii
   CALL pack(rz(i),mzz,mcb)
 END DO
 CALL CLOSE(mzz,1)
 CALL wrttrl(mcb)
 
!     IF A FREE SURFACE EXISTS - EXPAND THE MASS MATRIX
!     THE PARTITIONING VECTOR WILL BE SAVED FOR DMAP USE
 
 IF(nofree < 0) THEN
   GO TO   210
 END IF
 200 uset = usetd
 CALL calcv(pout,uh,uz,ufr,z(1))
 nsub0s = nsub0
 nsub1s = nsub1
 CALL gfsmrg(mhhbar,mzz,0,0,0,pout,pout)
 GO TO 220
 
 210 CALL gfswch(mhhbar,mzz)
 
!     COMPUTE THE FINAL MASS MATRIX
 
 220 CALL ssg2b(ajh,gjh,mhhbar,mmat,1,2,1,scr2)
 
!     IF GRAVITY EXISTS - TRANSFORM THE ADDITIONAL STIFFNESS AND
!     ADD IT IN.  BE SURE TO USE ONLY THOSE MODES REQUESTED IN
!     THE TRANSFORMATION FROM PHIA
 
 IF(nograv < 0) THEN
   GO TO   260
 END IF
 230 uset = usetd
 IF(lmodes >= nmodes) GO TO 240
 CALL calcv(pvec,um,uz,unz,z(1))
 CALL gfsptn(phia,phiar,0,0,0,pvec,0)
 GO TO 250
 
 240 phiar = phia
 
 250 CALL ssg2b(phiar,dkaa,0,scr2,1,2,1,scr5)
 CALL ssg2b(scr2,phiar,kzz,kzzbar,0,2,1,scr10)
 GO TO 270
 
 260 CALL gfswch(kzz,kzzbar)
 
!     IF A FREE SURFACE EXISTS - MERGE THE FREE SURFACE STIFFNESS IN
 
 270 IF(nofree < 0) THEN
   GO TO   290
 END IF
 280 nsub0 = nsub0s
 nsub1 = nsub1s
 CALL gfsmrg(khhbar,kzzbar,0,0,dkfrfr,pout,pout)
 GO TO 300
 
 290 CALL gfswch(khhbar,kzzbar)
 
!     COMPUTE THE FINAL STIFFNESS MATRIX BY ADDING IN COMPRESSIBILITY
!     IF IT EXISTS
 
 300 IF(sfbit == 0.0) THEN
   GO TO   310
 ELSE
   GO TO   320
 END IF
 310 badd(1) = 2
 dbadd(1) = 1.0D0
 badd(7) = 2
 dbadd(4) = 1.0D0
 CALL ssg2c(khhbar,kc,kmat,0,badd)
 
 GO TO 330
 320 CALL gfswch(khhbar,kmat)
 
!     TRANSFORM THE FINAL PRESSURE TRANSFORMATION MATRIX OR IF SPC
!     POINTS EXIST ON THE FLUID MERGE IN ZEROS
 
 330 uset = usetf
 IF(sfbit == 0.0) THEN
   GO TO   340
 ELSE
   GO TO   350
 END IF
 340 CALL ssg2b(h,gjh,0,gyh,1,2,1,scr5)
 GO TO 360
 
 350 CALL calcv(pvec,uy,uf,us,z(1))
 CALL gfsmrg(gyh,gjh,0,0,0,0,pvec)
 
!     PARTITION OUT THE FREE SURFACE POINTS
 
 360 IF(nofree < 0) THEN
   GO TO   380
 END IF
 370 CALL calcv(pvec,uy,ufr,ui,z(1))
 CALL gfsptn(gyh,0,gia,0,0,0,pvec)
 RETURN
 
 380 CALL gfswch(gyh,gia)
 RETURN
 
!     ERROR EXITS
 
 1001 n = -1
 GO TO 9999
 1002 n = -2
 GO TO 9999
 1008 n = -8
 
 9999 CALL mesage(n,FILE,NAME)
 RETURN
END SUBROUTINE gfsmo2
