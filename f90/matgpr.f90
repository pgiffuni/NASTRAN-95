SUBROUTINE matgpr
     
!     DMAP FOR MATGPR MODULE
 
!     MATGPR    GPL,USET,SIL,KFS//C,N,F/C,N,S/C,N,PRTOPT/
!                                 C,N,FILTER/C,N,FLTRFLAG $
 
!     THIS MODULE ENHANCED BY P.R.PAMIDI/RPK CORPORATION, 3/1988
 
 EXTERNAL        andf
 INTEGER :: gpl,uset,sil,im(7),two1,andf,sysbuf,core,BLANK,  &
     otpe,tycomp,scalar,comps(6),exid,prbuf(15),  &
     head2(32),iprbf(4),ICHAR(17),prbufc(5)
 INTEGER :: NAME(2),extra,hset
 REAL :: a(4),prbufx(5),xxbuf(15)
 CHARACTER (LEN=25) :: uwm
 CHARACTER (LEN=23) :: ufm
 COMMON /xmssg / ufm,uwm
 COMMON /condas/ idum(2),raddeg
 COMMON /system/ sysbuf,otpe,inx(6),nlpp,inx1(2),line
 COMMON /zzzzzz/ core(1)
 COMMON /bitpos/ ibits(32),ICHAR
 COMMON /output/ head(96),label(96)
 COMMON /zntpkx/ ia(4),ii,ieol,IEOR
 COMMON /BLANK / iiset(2),kkset(2),ipopt(2),filter,iflflg
 COMMON /two   / two1(32)
 EQUIVALENCE     (xxbuf(1),prbuf(1))
 EQUIVALENCE     (prbufc(1),prbufx(1))
 EQUIVALENCE     (ia(1), a(1))
 DATA    gpl,uset,sil,matrx    / 101   ,102   ,103   ,104     /
 DATA    scalar,comps,nline    / 4H s  ,4HT1  ,4HT2  ,4HT3    ,  &
     4HR1  ,4HR2  ,4HR3  ,15      /
 DATA    NAME  / 4HMATG, 4HPR  /
 DATA    null  / 4HNULL        /
 DATA    BLANK,  extra,  hset  / 4H    ,4H e  ,4H h           /
 DATA    ihset / 4HH           /
 DATA    iallp / 4HALLP        /
 DATA    head2 / 4H    , 4HPOIN, 4HT   , 4H    , 4H   v, 4HALUE, 4H    ,  &
     4H poi, 4HNT  , 4H    , 4H    , 4HVALU, 4HE   , 4H  p0,  &
     4HINT , 4H    , 4H    , 4H val, 4HUE  , 4H   p, 4HOINT,  &
     4H    , 4H    , 4H  va, 4HLUE , 4H    , 4H poi, 4HNT  ,  &
     4H    , 4H    , 4HVALU, 4HE   /
 
 
 iset = iiset(1)
 kset = kkset(1)
 inlopt = 0
 IF (ipopt(1) == null) inlopt = 1
 IF (filter   ==  0.0) GO TO 5
 iflag = 1
 IF (filter < 0.0) iflag = 2
 IF (iflflg /=   0) iflag = iflag + 2
 5 im(1) = matrx
 CALL rdtrl (im(1))
 IF (im(1) < 0) GO TO 380
 
!     CONVERT BCD TO BIT POSITION IN USET
 
 DO  i = 1,32
   IF (ICHAR(i) == iset) GO TO 20
 END DO
 IF (iset /= ihset) GO TO 15
 iset = -1
 GO TO 21
 15 WRITE  (otpe,16) uwm,iiset
 16 FORMAT (a25,', UNKNOWN SET ',2A4,' SPECIFIED FOR THE FIRST PARA',  &
     'METER OF THE MATGPR MODULE.  MODULE NOT EXECUTED.')
 RETURN
 
 20 iset = ibits(i)
 21 CONTINUE
 DO  i = 1,32
   IF (ICHAR(i) == kset) GO TO 40
 END DO
 kset = iset
 GO TO 50
 40 kset = ibits(i)
 50 CONTINUE
 lcore = korsz(core) - sysbuf
 ibuf  = lcore + 1
 IF (iset+kset == -2)  GO TO 51
 CALL gopen (gpl,core(ibuf),0)
 CALL READ  (*460,*60,gpl,core(1),lcore,0,lgpl)
 CALL CLOSE (gpl,1)
 GO TO 500
 
!     NSET ONLY  NO GPL,USET,  ETC.
 
 51 lgpl  = 0
 luset = 0
 lsil  = 0
 iuset = 1
 isil  = 1
 GO TO 81
 60 CALL CLOSE (gpl,1)
 lcore = lcore - lgpl
 CALL gopen (uset,core(ibuf),0)
 iuset = lgpl + 1
 CALL READ  (*480,*70,uset,core(lgpl+1),lcore,0,luset)
 CALL CLOSE (uset,1)
 GO TO 500
 70 CALL CLOSE (uset,1)
 lcore = lcore - luset
 CALL gopen (sil,core(ibuf),0)
 isil  = lgpl + luset + 1
 CALL READ  (*490,*80,sil,core(isil),lcore,0,lsil)
 CALL CLOSE (sil,1)
 GO TO 500
 80 CALL CLOSE (sil,1)
 k = isil + lsil
 lcore   = lcore - lsil - 1
 core(k) = luset + 1
 
!     LOAD HEADER FOR PAGES
 
 lsil = lsil + 1
 81 CONTINUE
 DO  i = 1,96
   label(i) = BLANK
 END DO
 DO  i = 1,32
   k = 32 + i
   label(k) = head2(i)
 END DO
 ncol  = im(2)
 CALL fname (matrx,label(4))
 CALL gopen (matrx,core(ibuf),0)
 ie    = ibits(12)
 inull = 0
 loop  = 0
 icmpx = 1
 IF (im(5) > 2) icmpx = 3
 IF (iset /= -1) mask  = two1(iset)
 IF (kset /= -1) mask1 = two1(kset)
 muset = 0
 jc    = 0
 iksil = 1
 l     = 1
 ASSIGN 210 TO iout
 CALL page
 
!     START LOOP ON EACH COLUMN
 
 110 loop = loop + 1
 CALL intpk (*390,matrx,0,icmpx,0)
 IF (inull /= 0) GO TO 400
 120 CONTINUE
 IF (inlopt == 1) GO TO 359
 
!     CHECK FOR HSET
 
 121 IF (iset  ==   -1) GO TO 150
 IF (muset == loop) GO TO 160
 130 jc = jc + 1
 IF (jc > luset) GO TO 150
 kk = lgpl + jc
 IF (andf(core(kk),mask) == 0.0) THEN
   GO TO   130
 END IF
 
!     FOUND COLUMN IN USET
 
 140 muset = muset + 1
 GO TO 121
 
!     COLUMN NOT IN USET
 
 150 iprbf(l  ) = loop
 iprbf(l+1) = hset
 GO TO 200
 
!     JC IS INDEX OF NON-ZERO IN G SET-- SOOK UP SIL
 
 160 IF (iksil == lsil+1) GO TO 150
 kk = isil + iksil
 IF (jc < core(kk)) GO TO 170
 iksil = iksil + 1
 GO TO 160
 170 icomp = jc - core(kk-1) + 1
 IF (icomp /= 1) GO TO 180
 
!     CHECK FOR SCALAR POINT
 
 IF (core(kk)-core(kk-1) > 1) GO TO 180
 tycomp = scalar
 
!     CHECK FOR EXTRA
 
 kk = lgpl + jc
 IF (andf(core(kk),two1(ie)) == 0.0) THEN
   GO TO   190
 END IF
 171 tycomp = extra
 GO TO 190
 180 tycomp = comps(icomp)
 190 exid   = core(iksil)
 iprbf(l+1) = tycomp
 iprbf(l  ) = exid
 200 GO TO iout, (210,420,430)
 210 WRITE  (otpe,220)loop,iprbf(1),iprbf(2)
 220 FORMAT ('0COLUMN',i8,2H (,i8,1H-,a2,2H).)
 line  = line + 2
 IF (line >= nlpp) CALL page
 jj    = 0
 kuset = 0
 ksil  = 1
 ipb   = 1
 ipbc  = 1
 iend  = 0
 230 IF (ieol == 0) THEN
   GO TO   240
 ELSE
   GO TO   350
 END IF
 240 CALL zntpki
 
!     CHECK FILTER
 
 IF (filter == 0.0) GO TO 246
 
!     FILTER IS NON-ZERO
 
 value = a(1)
 IF (icmpx == 3) value = SQRT(a(1)*a(1) + a(2)*a(2))
 SELECT CASE ( iflag )
   CASE (    1)
     GO TO 241
   CASE (    2)
     GO TO 242
   CASE (    3)
     GO TO 243
   CASE (    4)
     GO TO 244
 END SELECT
 
 241 IF (ABS(value) < filter) GO TO 230
 GO TO 246
 242 IF (ABS(value) > ABS(filter)) GO TO 230
 GO TO 246
 243 IF (value < filter .AND. value > 0.0) GO TO 230
 GO TO 246
 244 IF (value > filter .AND. value < 0.0) GO TO 230
 
!     CHECK FOR HSET
 
 246 IF (kset == -1) GO TO 306
 
!     LOOK UP ROW IN USET
 
 250 IF (kuset > luset+1) GO TO 500
 IF (kuset ==      ii) GO TO 280
 260 jj = jj + 1
 
!     PROTECT AGINST NO BITPOS OR NO USET
 
 IF (jj > luset) GO TO 306
 kk = lgpl + jj
 IF (andf(core(kk),mask1) == 0.0) THEN
   GO TO   260
 END IF
 
!     FOUND ELEMENT IN USET
 
 270 kuset = kuset + 1
 GO TO 250
 
!     JJ IS INDEX OF NON-ZERO IN G SET - NOW SEARCH SIL FOR JJ
 
 280 IF (ksil == lsil+1) GO TO 510
 kk = isil + ksil
 IF (jj < core(kk)) GO TO 290
 ksil = ksil + 1
 GO TO 280
 290 icomp = jj - core(kk-1) + 1
 IF (icomp /= 1) GO TO 300
 
!     CHECK FOR SCALAR POINT
 
 IF (core(kk)-core(kk-1) > 1) GO TO 300
 tycomp = scalar
 
!     CHECK FOR EXTRA POINT
 
 kk = lgpl + jj
 IF (andf(core(kk),two1(ie)) == 0.0) THEN
   GO TO   310
 END IF
 
!     EXTRA POINT
 
 305 tycomp = extra
 GO TO 310
 
!     H POINT
 
 306 tycomp = hset
 exid   = ii
 GO TO 311
 300 tycomp = comps(icomp)
 310 exid   = core(ksil)
 311 IF (ipb >= nline) GO TO 330
 320 prbuf(ipb  ) = exid
 prbuf(ipb+1) = tycomp
 IF (icmpx == 1) GO TO 325
 IF (ipopt(1) /= iallp) GO TO 325
 amag = SQRT(a(1)*a(1) + a(2)*a(2))
 IF (amag == 0.0) GO TO 325
 a(2) = ATAN2(a(2),a(1))*raddeg
 IF (a(2) < -0.00005) a(2) = a(2) + 360.0
 a(1) = amag
 325 prbuf(ipb+2) = ia(1)
 prbufc(ipbc) = ia(2)
 ipbc = ipbc + 1
 ipb  = ipb  + 3
 GO TO 230
 330 ipb1 = ipb  - 1
 ipbc = ipbc - 1
 WRITE  (otpe,340) (prbuf(i),prbuf(i+1),xxbuf(i+2),i=1,ipb1,3)
 340 FORMAT (5X,5(1X,i8,1X,1A2,1X,1P,e12.5))
 line = line + 1
 IF (icmpx == 1) GO TO 343
 WRITE  (otpe,341) (prbufx(i),i=1,ipbc)
 341 FORMAT (5X,5(13X,1P,e12.5))
 WRITE  (otpe,342)
 342 FORMAT (1H )
 line = line + 2
 343 CONTINUE
 ipbc = 1
 ipb  = 1
 IF (line >= nlpp) CALL page
 IF (iend == 1) GO TO 360
 GO TO 320
 
!     END OF COLUMN
 
 350 iend = 1
 IF (ipb == 1) GO TO 360
 GO TO 330
 359 CALL fwdrec (*510,matrx)
 360 IF (loop  /= ncol) GO TO 110
 IF (inull /=    0) GO TO 450
 370 CALL CLOSE (matrx,1)
 380 RETURN
 
 390 IF (inull /= 0) GO TO 360
 inull = 1
 ibegn = loop
 GO TO 360
 400 ifin  = loop - 1
 inull = 0
 410 loops = loop
 loop  = ibegn
 ASSIGN 420 TO iout
 GO TO 121
 420 l = 3
 loop = ifin
 ASSIGN 430 TO iout
 GO TO 121
 430 ASSIGN 210 TO iout
 l = 1
 loop = loops
 WRITE  (otpe,440) ibegn,iprbf(1),iprbf(2),ifin,iprbf(3),iprbf(4)
 440 FORMAT ('0COLUMNS',i8,2H (,i8,1H-,a2,6H) thru,i8,2H (,i8,1H-,a2,  &
     11H) are null.)
 line = line + 2
 IF (line >= nlpp) CALL page
 IF (ifin /= ncol) GO TO 120
 GO TO 370
 450 ifin = loop
 GO TO 410
 460 in = gpl
 470 CALL mesage (-2,in,NAME)
 480 in = uset
 GO TO 470
 490 in = sil
 GO TO 470
 500 CALL mesage (8,0,NAME)
 GO TO 370
 510 CALL mesage (7,0,NAME)
 GO TO 370
END SUBROUTINE matgpr
