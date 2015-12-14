SUBROUTINE ddr2
     
!     DYNAMIC DATA RECOVERY--PART 2 --MODE ACCELERATION
 
!     DMAP SEQUENCE
 
!     INPUTS = 9
 
!     USETD,VUD,PD,K2DD,B2DD,MDD,FRL,LLL,DM
 
!     OUTPUTS = 3
 
!     UAV,UEV,PAF
 
!     SCRATCHES = 6
 
!     PARAMETERS 1 BCD, 3INTEGERS
 
 INTEGER :: usetd,pd,b2dd,frl,dm,uav,paf, scr2,scr3,scr4,scr5,scr6,scr7,  &
     TYPE,react,tran,uset,vud,pad,uev,pl
 COMMON /BLANK / TYPE(2),noue,react,frqset
 COMMON /bitpos/ um,uo,ur,usg,usb,ul,ua,uf,us,un,ug,ue,up,une,ufe, ud
 COMMON /patx  / lc,n,no,n4,uset
 COMMON /zzzzzz/ core(1)
 DATA    usetd , vud, pd, k2dd,b2dd,mdd,frl,lll,dm  /  &
     101   , 102,103,  104, 105,106,107,108,109 /
 DATA    uav   , uev, paf, tran   / 201   , 202, 203, 4HTRAN /
 DATA    scr2  , scr3,scr4,scr5,scr6,scr7,pad /  &
     302   , 303 , 304, 305, 306, 301,302 /
 
 
 lc   = korsz(core)
 vud  = 102
 scr7 = 301
 uset = usetd
 pl   = scr6
 isol = scr7
 IF (noue >= 0) GO TO 10
 pad  = paf
 10 CONTINUE
 IF (TYPE(1) /= tran) scr7 = uav
 IF (TYPE(1) /= tran .AND. react < 0 .AND. noue >= 0) scr7 = vud
 
!     MODE ACCELERATION
 
!     FORM PAD
 
 
 CALL ddr1a (pd,k2dd,b2dd,mdd,vud,pad,frl,frqset,scr3,scr4,scr5,  &
     scr6,TYPE(1),scr7)
 
!     DISP ON SCR7 IN TRANSIENT
 
 IF (noue < 0) GO TO 50
 CALL calcv (scr3,ud,ua,ue,core(1))
 CALL ssg2a (vud,scr4,uev,scr3)
 
!     UA IS ON SCR4
 
 vud = scr4
 
!     BREAK UP PAD
 
 CALL ssg2a (pad,paf,scr5,scr3)
 50 IF (react >= 0) GO TO 90
 
!     UR NULL
 
 IF (TYPE(1) /= tran) scr7 = isol
 IF (TYPE(1) /= tran .AND. noue < 0) scr7 = uav
 CALL ssg3a (0,lll,paf,scr7,scr3,scr6,-1,0)
 60 IF (TYPE(1) /= tran) GO TO 80
 
!     MERGE RECALCULATED SOLUTIONS AND ACCEL AND VELOCITY
 
 isol = uav
 IF (noue < 0) GO TO 70
 isol = scr5
 70 CALL ddr1b (vud,scr7,isol)
 
!     BUILD UP TO DSIZE  ADDING IN UEV
 
 80 IF (noue < 0) GO TO 30
 CALL sdr1b (scr4,isol,uev,uav,ud,ua,ue,usetd,0,0)
 30 RETURN
 
!     FREE BODY PROBLEM
 
 90 CALL calcv (scr3,ua,ul,ur,core(1))
 
!     PARTITION PAF AND UA
 
 CALL ssg2a (paf,pl,scr5,scr3)
 ivec = vud
 IF (TYPE(1) == tran) ivec = scr7
 CALL ssg2a (ivec,scr2,scr5,scr3)
 
!     UR IS  ON SCR5
 
 CALL ssg3a (0,lll,pl,scr3,scr2,scr6,-1,0)
 CALL ssg2b (dm,scr5,scr3,scr4,0,2,1,scr6)
 CALL sdr1b (scr3,scr4,scr5,scr7,ua,ul,ur,usetd,0,0)
 GO TO 60
END SUBROUTINE ddr2
