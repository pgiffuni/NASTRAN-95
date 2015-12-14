SUBROUTINE sdr1
     
 EXTERNAL        andf
 INTEGER :: andf,ue,reig,uset,pg,ulv,uoov,ys,GO,gm,ps,qr,  &
     ugvx,pgx,qsx,um,uo,ur,us,ia(7),  &
     ua1,uf1,un1,ug1,up,une,ufe,ud,ua,uf,un,ug,dyna
 COMMON /two   / two1(32)
 COMMON /BLANK / APPEND,itype(2)
 COMMON /bitpos/ um,uo,ur,usg,usb,ul,ua1,uf1,us,un1,ug1,ue,up,une, ufe,ud
 COMMON /system/ isys(54),iprec,iheat
 EQUIVALENCE     (isys(25),irfno)
 DATA    dyna  , reig / 4HDYNA, 4HREIG /
 
 ium   = 304
 iscr6 = 306
 iur   = 0
 iuo   = 304
 ipvect= 301
 ius   = 304
 uset  = 101
 pg    = 102
 ulv   = 103
 uoov  = 104
 ys    = 105
 GO    = 106
 gm    = 107
 ps    = 108
 kss   = 110
 qr    = 111
 ugvx  = 201
 pgx   = 202
 qsx   = 203
 kfs = 109
 iua = 302
 iuf = 303
 iun = 302
 iug = 306
 iscr5 = 0
 
!     COPY PG ONTO PGX
 
 CALL sdr1a (pg,pgx)
 
!     SET FLAGS TO CONTROL LOGIC
 
 ia(1) = uset
 CALL rdtrl (ia(1))
 IF (ia(1) <= 0) RETURN
 iomt   = andf(ia(5),two1(uo))
 noue   = andf(ia(5),two1(ue))
 isng   = andf(ia(5),two1(us))
 ireact = andf(ia(5),two1(ur))
 imulti = andf(ia(5),two1(um))
 itran  = 1
 
!     TEST FOR DYNAMICS OR STATICS
 
 IF (noue /= 0 .OR. itype(1) == dyna) GO TO 10
 
!     STATICS
 
 ua = ua1
 uf = uf1
 un = un1
 ug = ug1
 GO TO 20
 
!     DYNAMICS
 
 10 ug = up
 un = une
 uf = ufe
 ua = ud
 IF (iheat /= 0) itran = 0
 IF (irfno == 9) itran = 0
 20 CONTINUE
 
!     IF REAL EIGENVALUE,BUCKLING,OR DYNAMICS PROBLEM UR = 0
 
 IF (itype(1) == dyna .OR. itype(1) == reig) GO TO 70
 IF (ireact == 0) THEN
   GO TO    70
 END IF
 
!     REACTIONS
 
 40 CALL sdr1b (ipvect,ulv,iur,iua,ua,ul,ur,uset,0,0)
 iscr5 = 305
 IF (isng == 0) THEN
   GO TO    60
 END IF
 50 CONTINUE
 CALL sdr1b (ipvect,qr,0,iscr5,uf,ur,ul,uset,0,0)
 iug = iua
 GO TO 80
 
!     REACTS BUT NO SINGLES - MAKE QG
 
 60 CALL sdr1b (ipvect,qr,0,iscr5,ug,ur,ul,uset,0,0)
 CALL sdr1a (iscr5,qsx)
 GO TO 80
 
!     NO REACT
 
 
!     NON STATICS APPROACH
 
 70 iua = ulv
 80 IF (iomt == 0) THEN
   GO TO   100
 END IF
 
!     OMITTED POINTS
 
 90 CALL ssg2b (GO,iua,uoov,iuo,0,iprec,1,iscr6)
 CALL sdr1b (ipvect,iua,iuo,iuf,uf,ua,uo,uset,0,0)
 iug = iuf
 GO TO 110
 
!     NO OMITTED POINTS
 
 100 isav = iuf
 iuf  = iua
 iun  = isav
 110 IF (isng == 0) THEN
   GO TO   180
 END IF
 
!     SINGLE POINT CONSTRAINTS
 
 
!     TEST FOR PRESENCE OF YS VECTOR
 
 120 ia(1) = ys
 CALL rdtrl (ia(1))
 IF (ia(1) < 0 .OR. ia(6) == 0) GO TO 130
 CALL sdr1b (ipvect,iuf,ys,iun,un,uf,us,uset,1,ius)
 
!     IUS CONTAINS EXPANDED YS FROM SPC
 
 
!     IS QS REWUESTED
 
 ia(1) = qsx
 CALL rdtrl (ia(1))
 IF (ia(1) <= 0) GO TO 190
 
!     COMPUTE QS
 
 CALL ssg2b (kss,ius,ps,ipvect,0,iprec,2,iscr6)
 CALL ssg2b (kfs,iuf,ipvect,ius,1,iprec,1,iscr6)
 IF (imulti /= 0 .AND. ireact /= 0) GO TO 160
 CALL sdr1b (ipvect,ius,iscr5,iscr6,ug,us,uf,uset,0,0)
 CALL sdr1a (iscr6,qsx)
 GO TO 190
 
!     NO YS VECTOR
 
 130 CALL sdr1b (ipvect,iuf,0,iun,un,uf,us,uset,0,0)
 ia(1) = qsx
 CALL rdtrl (ia(1))
 IF (ia(1) <= 0) GO TO 190
 
!     COMPUTE QS = KFS T*UF
 
 iuf1 = iuf
 IF (itype(1) /= dyna) GO TO 140
 
!     EXPAND  KFS TO  D SET
 
 IF (noue == 0) GO TO 140
 CALL sdr1b (ipvect,kfs,0,ius,uf,uf1,ue,uset,0,0)
 kfs = ius
 
!     IF TRANSIENT STRIP VELOCITY AND ACCERERATION FROM IUF
 
 140 CALL sdr1d (ps,iuf,qsx,itran)
 IF (itran == 1) GO TO 150
 iuf1 = qsx
 150 CALL ssg2b (kfs,iuf1,ps,ipvect,1,iprec,2,iscr6)
 IF (imulti /= 0 .AND. ireact /= 0 .AND. itype(1) /= dyna .AND.  &
     itype(1) /= reig) GO TO 170
 CALL sdr1b (ius,ipvect,iscr5,iscr6,ug,us,uf,uset,0,0)
 CALL sdr1a (iscr6,qsx)
 GO TO 190
 160 CALL sdr1b (ipvect,ius,iscr5,iscr6,un,us,uf,uset,0,0)
 CALL sdr1b (ipvect,iscr6,0,ius,ug,un,um,uset,0,0)
 CALL sdr1a (ius,qsx)
 GO TO 190
 170 CALL sdr1b (ius,ipvect,iscr5,iuf,un,us,uf,uset,0,0)
 CALL sdr1b (ius,iuf,0,ipvect,ug,un,um,uset,0,0)
 CALL sdr1a (ipvect,qsx)
 GO TO 190
 
!     NO SINGLE POINT CONSTRAINTS
 
 180 iug = iun
 iun = iuf
 
 190 IF (imulti == 0) THEN
   GO TO   200
 ELSE
   GO TO   210
 END IF
 
!     NO MULTI POINT CONSTRAINTS
 
 200 iug = iun
 GO TO 220
 
!     MULTI POINT CONSTRAINTS
 
 210 iug = iscr6
 CALL ssg2b (gm,iun,0,ium,0,iprec,1,iscr6)
 CALL sdr1b (ipvect,iun,ium,iug,ug,un,um,uset,0,0)
 220 CALL sdr1a (iug,ugvx)
 RETURN
END SUBROUTINE sdr1
