SUBROUTINE trlg
     
!     THIS IS THE MODULE DRIVER FOR  TRLG(TRANSIENT LOAD GENERATOR)
 
!     INPUTS(14)
!       CASEXX      CASECONTROL
!       USETD
!       DLT         DYNAMIC LOAD TABLE
!       SLT         STATIC  LOAD TABLE
!       BGPDT       BASIC GRID POINT DEFINITION TABLE
!       SIL         SCALAR INDEX LIST
!       CSTM        COORDINATE SYSTEMS
!       TRL         TRANSIENT RESPONSE LIST
!       DIT         DIRECT INPUT TABELS
!       GMD
!       GOD
!       PHIDH
!       EST
!       MGG         MASS MATRIX FOR GRAVITY LOADS
!       MPT
!     OUTPUTS(6)
!       PPO
!       PSO
!       PDO
!       PD
!       PH
!       TOL
!     PARAMETERS
!      IP1 = -1     IF (AP = AD)
!      NCOL.LE.0    NO CONTINUE MODE (TO = 0.0)
!      NCOL.GT.0    CONTINUE MODE (TO = LAST TIME)
 
!     SCRATCHES (9)
 
 INTEGER :: casexx,usetd,dlt,slt,bgpdt,sil,cstm,trl,dit,gmd,god,  &
     phidh,ppo,pso,pdo,pd,ph,tol,scr1,scr2,scr3,scr4,scr5,  &
     scr6,scr7,ap,as,ad,ah,tmldtb,fct,fco,scr8,est,scr9,mcb(7)
 
 COMMON /BLANK/ ip1,ncol
 
 DATA  casexx,usetd,dlt,slt,bgpdt,sil,cstm,trl,dit,gmd,god,phidh/  &
     101   , 102 ,103,104,105  ,106,107 ,108,109,110,111,112  /
 DATA  est,mgg,mpt / 113,114,115 /
 DATA  ppo,pso,pdo,pd ,ph ,tol / 201,202,203,204,205,206 /
 DATA  scr1,scr2,scr3,scr4,scr5,scr6,scr7,scr8,scr9 /  &
     301 ,302 ,303 ,304 ,305, 306 ,307 ,308 ,309  /
 
!     FORM AP MATRIX AND EXTRACT LOAD TABLES
 
 ap = scr1
 tmldtb = scr2
 CALL trlga (casexx,usetd,dlt,slt,bgpdt,sil,cstm,ap,tmldtb,itrl,  &
     scr3,scr4,scr5,est,scr6,mgg,scr7,mpt)
 
!     REDUCE TRANSFORMATION MATRIX
 
 as = scr3
 ad = scr4
 ah = scr5
 mcb(1) = ap
 CALL rdtrl (mcb)
 IF (mcb(2) <= 0) GO TO 10
 CALL trlgb (usetd,ap,gmd,god,phidh,as,ad,ah,ip1,scr6,scr7,scr8, scr9)
 
!     PRODUCE TIME FUNCTION MATRIX
 
 10 CONTINUE
 fct = scr6
 fco = scr7
 CALL  trlgc (tmldtb,trl,dit,itrl,fct,fco,tol,ip2)
 IF (mcb(2) <= 0) GO TO 20
 IF (ip2  ==  -1) fco = fct
 
!     COMPUTE LOAD FACTORS
 
 CALL trlgd (fct,fco,ap,as,ad,ah,ppo,pso,pdo,pd,ph,ip1,scr2,ip2)
 20 CONTINUE
 RETURN
END SUBROUTINE trlg
