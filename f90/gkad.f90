SUBROUTINE gkad
     
!     GENERAL K ASSEMBLER DIRECT
 
!     INPUT = 10,  USETD,GM,GO,KAA,BAA,MAA,K4AA,K2PP,M2PP,B2PP
!     OUTPUT = 8,  KDD,BDD,MDD,GMD,GOD,K2DD,M2DD,B2DD
!     SCRATCHES = 6
!     PARAMETERS 3 BCD, 3 REAL, 11 INTERGER
!     - TYPE,APP,FORM, G,W3,W4, NOK2PP,MOM2PP,NOB2PP,MULTI,SINGLE,OMIT,
!       NOUE,NOK4GG,NOBGG,NOKMGG,MODACC
 
 
 INTEGER :: TYPE(2),app(2),FORM(2),iblock(11),blck(12),mcb(7),  &
     tran,forc,omit,baa,b2pp,b2dd,b1dd,scr1,scr2,scr3,  &
     scr4,scr5,scr6,gm,GO,god,gmd,bdd,usetd,single
 DOUBLE PRECISION :: BLOCK(5)
 COMMON /BLANK / TYPE,app,FORM, g,w3,w4, ik2pp,im2pp,ib2pp,multi,  &
     single,omit,noue,nok4gg,nobgg,nokmgg,modacc
 COMMON /bitpos/ um,uo,ur,usg,usb,ul,ua,uf,us,un,ug,ue,up,une,ufe, ud
 EQUIVALENCE     (iblock(1),blck(2)),(BLOCK(1),blck(3))
 DATA    usetd , gm, GO, kaa,baa,maa,k4aa,k2pp,m2pp,b2pp /  &
     101   , 102,103,104,105,106,107, 108, 109, 110  /
 DATA    kdd   , bdd,mdd,gmd,god,k2dd,m2dd,b2dd /  &
     201   , 202,203,204,205,206, 207, 208  /
 DATA    scr1  , scr2,scr3,scr4,scr5,scr6 / 301   , 302, 303, 304, 305, 306  /
 DATA    forc  , tran,modal / 4HFORC,4HTRAN,4HMODA     /
 DATA    BLOCK(1),BLOCK(2),BLOCK(4),BLOCK(5),iblock(1),iblock(7) /  &
     1.0D0 , 0.0D0,    1.0D0,   0.0D0,   2,        2         /
 DATA    xnum  , mcb / 1.0,7*0  /  ,iblock(6) /    -1  /
 
 
 kdd   = 201
 bdd   = 202
 mdd   = 203
 k2dd  = 206
 m2dd  = 207
 b2dd  = 208
 k1dd  = 302
 m1dd  = 303
 b1dd  = 304
 k41dd = 305
 scr3  = 303
 scr4  = 304
 IF (noue > 0) GO TO 10
 
!     NO E-S A = 1DD
 
 k1dd  = kaa
 b1dd  = baa
 m1dd  = maa
 k41dd = k4aa
 10 IF (TYPE(1) == tran) GO TO 20
 
!     COMPLEX EIGENVALUE OR FREQUENCY RESPONSE - SET UP FOR FINAL ADD
 
 IF (ib2pp < 0) b1dd = bdd
 IF (im2pp < 0) m1dd = mdd
 GO TO 50
 
!     TRANSIENT ANALYSIS - SETUP FOR FINAL ADD
 
 20 IF (ik2pp < 0) k1dd = kdd
 IF (im2pp < 0) m1dd = mdd
 IF (w3 /=  0.0) GO TO 30
 g  = 0.0
 w3 = 1.0
 30 IF (w4 /= 0.0) GO TO 50
 w4   = 1.0
 xnum = 0.0
 50 IF (app(1) /= forc) GO TO 60
 
!     FORCE APPROACH P = D
 
 k2dd = k2pp
 b2dd = b2pp
 m2dd = m2pp
 GO TO 140
 
!     DISPLACEMENT APPROACH - REDUCE P TO D
 
!     IF MODAL DO NOT MAKE KDD AND BDD
 
 60 IF (FORM(1) /= modal) GO TO 70
 kdd  = 0
 k1dd = 0
 bdd  = 0
 b1dd = 0
 70 IF (noue < 0) GO TO 100
 
!     BUILD GMD AND GOD
 
!     M-S PRESENT
 
 IF (multi >= 0) CALL gkad1a (usetd,gm,gmd,scr1,ue,un,une)
 
!     0-S PRESENT
 
 IF (omit >= 0) CALL gkad1a (usetd,GO,god,scr1,ue,ua,ud)
 
 100 IF (multi < 0 .AND. single < 0 .AND. omit < 0) GO TO 130
 IF (im2pp < 0 .AND. ib2pp < 0 .AND. ik2pp < 0) GO TO 130
 
!     REDUCE 2PP-S TO 2DD-S
 
 CALL gkad1c (gmd,god,scr1,scr2,scr3,scr4,scr5,scr6,usetd)
 IF (ik2pp >= 0) CALL gkad1d (k2pp,k2dd)
 IF (im2pp >= 0) CALL gkad1d (m2pp,m2dd)
 IF (ib2pp >= 0) CALL gkad1d (b2pp,b2dd)
 130 IF (FORM(1) == modal .AND. modacc < 0) GO TO 180
 IF (noue < 0) GO TO 140
 
!     EXPAND AA-S TO DD SET
 
 CALL gkad1b (usetd,kaa,maa,baa,k4aa,k1dd,m1dd,b1dd,k41dd,ua,ue, ud,scr1)
 140 IF (TYPE (1) == tran) GO TO 190
 
!     FREQUENCY RESPONSE OR COMPLEX EIGENVALUE
 
 IF (b1dd == bdd .OR. nobgg < 0 .OR. FORM(1) == modal) GO TO 150
 CALL ssg2c (b1dd,b2dd,bdd,1,iblock(1))
 150 IF (m1dd == mdd .OR. nokmgg < 0) GO TO 160
 CALL ssg2c (m1dd,m2dd,mdd,1,iblock(1))
 160 IF (k1dd == kdd .OR. FORM(1) == modal .OR. nokmgg < 0) GO TO 180
 iblock(1) = 4
 BLOCK(2)  = g
 IF (nok4gg < 0) scr4 = kdd
 
!     DETERMINE IF KDD IS REAL OR IMAGINARY  (COMPLEX EIGEN)
 
 mcb(1) = k2dd
 CALL rdtrl (mcb(1))
 IF (g /= 0.0 .OR. nok4gg > 0 .OR. mcb(5) > 2) GO TO 170
 iblock(1) = 2
 iblock(7) = 2
 170 CALL ssg2c (k1dd,k2dd,scr4,1,iblock)
 IF (nok4gg < 0) GO TO 180
 BLOCK(1) = 0.0D0
 BLOCK(2) = 1.0D0
 CALL ssg2c (k41dd,scr4,kdd,1,iblock(1))
 180 RETURN
 
!     TRANSIENT ANALYSIS
 
 190 iblock(1) = 2
 iblock(7) = 2
 IF (k1dd == kdd .OR. nokmgg < 0) GO TO 200
 CALL ssg2c (k1dd,k2dd,kdd,1,iblock(1))
 200 IF (m1dd == mdd .OR. nokmgg < 0) GO TO 210
 CALL ssg2c (m1dd,m2dd,mdd,1,iblock(1))
 210 IF (b1dd == bdd) GO TO 180
 BLOCK(1) = g/w3
 BLOCK(4) = xnum/w4
 IF (g == 0.0 .AND. xnum == 0.0 .AND. nobgg < 0 .AND. ib2pp < 0) GO TO 180
 IF (nobgg < 0 .AND. ib2pp < 0) scr3 = bdd
 CALL ssg2c (k1dd,k41dd,scr3,1,iblock(1))
 IF (scr3 == bdd) GO TO 180
 BLOCK(1) = 1.0D0
 BLOCK(4) = 1.0D0
 CALL ssg2c (b1dd,b2dd,scr5,1,iblock(1))
 CALL ssg2c (scr5,scr3,bdd, 1,iblock(1))
 GO TO 180
END SUBROUTINE gkad
