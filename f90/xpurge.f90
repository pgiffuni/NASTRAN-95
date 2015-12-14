SUBROUTINE xpurge
     
!     THIS SUBROUTINE PURGES AND EQUATES FILES WITHIN FIAT AND DPD
 
 IMPLICIT INTEGER (a-z)
 EXTERNAL        lshift    ,andf    ,orf
 DIMENSION       purge1( 2),ddbn( 1),dfnu( 1),fcum( 1),fcus( 1),  &
     fdbn  ( 1),fequ( 1),FILE( 1),fknd( 1),fmat( 1),  &
     fntu  ( 1),fpun( 1),fon ( 1),ford( 1),minp( 1),  &
     mlsn  ( 1),mout( 1),mscr( 1),sal ( 1),sdbn( 1), sntu  ( 1),sord( 1)
 CHARACTER (LEN=25) :: sfm
 CHARACTER (LEN=29) :: uim
 CHARACTER (LEN=25) :: uwm
 CHARACTER (LEN=23) :: ufm
 COMMON /xmssg / ufm,uwm,uim,sfm
 COMMON /system/ ibufsz,outtap,dum(21),icfiat,dumm(14),nbpc,nbpw, ncpw
 COMMON /oscent/ x(1)
 COMMON /xfiat / fiat(7)
 COMMON /xfist / fist
 COMMON /xdpl  / dpd(6)
 COMMON /xvps  / vps(1)
 COMMON /xsfa1 / md(401),sos(1501),comm(20),xf1at(5)
 COMMON /ipurge/ j,k,nsav,prisav,hold
 EQUIVALENCE              (dpd  (1),dnaf    ),(dpd  (2),dmxlg   ),  &
     (dpd  (3),dculg   ),(dpd  (4),ddbn (1)),(dpd  (6),dfnu (1)),  &
     (fiat (1),funlg   ),(fiat (2),fmxlg   ),(fiat (3),fculg   ),  &
     (fiat (4),fequ (1)),(fiat (4),FILE (1)),(fiat (4),ford (1)),  &
     (fiat (5),fdbn (1)),(fiat (7),fmat (1)),(md   (1),mlgn    ),  &
     (md   (2),mlsn (1)),(md   (3),minp (1)),(md   (4),mout (1)),  &
     (md   (5),mscr (1)),(sos  (1),slgn    ),(sos  (2),sdbn (1)),  &
     (sos  (4),sal  (1)),(sos  (4),sntu (1)),(sos  (4),sord (1)),  &
     (xf1at(1),fntu (1)),(xf1at(1),fon  (1)),(xf1at(2),fpun (1)),  &
     (xf1at(3),fcum (1)),(xf1at(4),fcus (1)),(xf1at(5),fknd (1))
 EQUIVALENCE              (comm (1),almsk   ),(comm (2),apndmk  ),  &
     (comm (3),cursno  ),(comm (4),entn1   ),(comm (5),entn2   ),  &
     (comm (6),entn3   ),(comm (7),entn4   ),(comm (8),flag    ),  &
     (comm (9),fnx     ),(comm(10),lmsk    ),(comm(11),lxmsk   ),  &
     (comm(12),macsft  ),(comm(13),rmsk    ),(comm(14),rxmsk   ),  &
     (comm(15),s       ),(comm(16),scornt  ),(comm(17),tapmsk  ),  &
     (comm(18),thcrmk  ),(comm(19),zap     )
 DATA  purge1  /4HXPUR   ,4HGE     /
 
 CALL xsfadd
 k = -1
 GO TO 10
 
 
 ENTRY xequiv
!     ============
 
 CALL xsfadd
 k = +1
 prisav= 0
 secchn= 0
 10 entn1 = icfiat
 lmt1  = funlg* entn1
 lmt2  = lmt1 + 1
 lmt3  = fculg* entn1
 IF (fculg >= fmxlg) GO TO 610
 nfculg= lmt3 + 1
 incr  = 1
 
!     S = O 400000000000     Z 80000000
 
 20 s = lshift(1,nbpw-1)
 
!     INITIALIZE FOR FIRST SET OF DATA BLOCKS
 
 nwds = x(1)
 i = 7
 100 ndbs = x(i)
 
!     FIND POSITION OF VPS POINTER WORD
 
 jpt = i + 2*ndbs + 1
 IF (k == 1) jpt = jpt + 1
 ivps  = x(jpt)
 iexec = -1
 IF (ivps > 0) iexec = vps(ivps)
 
!     TEST CONDITIONAL INDICATOR (NOT HERE, BELOW TO PERMIT UNPURGE-
!     UNEQUIV
 
 GO TO 200
 
!     INCREMENT AND LOOK AT NEXT SET OF DATA BLOCKS
 
 150 i = jpt + 1
 IF (i < nwds) GO TO 100
 RETURN
 
!     TEST FOR PURGE OR EQUIV
 
 200 j = i + 1
 IF (k > 0) GO TO 400
 
!     PURGE LOGIC FOLLOWS
 
 220 xj1 = x(j  )
 xj2 = x(j+1)
 DO  n = 1,lmt3,entn1
   IF (xj1 /= fdbn(n) .OR. xj2 /= fdbn(n+1)) CYCLE
   IF (n  <= lmt1) GO TO 240
   IF (iexec >= 0) GO TO 230
   FILE(n) = zap
   GO TO 280
   
   230 IF (andf(rmsk,FILE(n)) /= zap) GO TO 300
   
!     UNPURGE (CLEAR THE ENTRY)
   
   lmt4 = n + entn1 - 1
   DO  m = n,lmt4
     FILE(m) = 0
   END DO
   GO TO 300
   240 IF (iexec >= 0) GO TO 300
   hold = andf(rxmsk,FILE(n))
   lmt4 = n + entn1 - 1
   DO  m = n,lmt4
     FILE(m) = 0
   END DO
   FILE(n) = hold
   
   EXIT
 END DO
 270 IF (iexec >= 0) GO TO 300
 FILE(nfculg  ) = zap
 fdbn(nfculg  ) = xj1
 fdbn(nfculg+1) = xj2
 fculg = fculg + incr
 IF (fculg >= fmxlg) GO TO 620
 nfculg = nfculg + entn1
 
 280 CALL xpolck (xj1,xj2,fn,l)
 IF (fn == 0) GO TO 300
 ddbn(l  ) = 0
 ddbn(l+1) = 0
 
 300 j = j + 2
 IF (j-2 == i+1 .AND. k > 0) j = j + 1
 IF (j < jpt) GO TO 220
 GO TO 150
 
!     EQUIV LOGIC FOLLOWS
 
 400 xj1 = x(j  )
 xj2 = x(j+1)
 DO  n = 1,lmt3,entn1
   IF (xj1 /= fdbn(n) .OR. xj2 /= fdbn(n+1)) CYCLE
   IF (j /= i+1) GO TO 420
   
!     PRIMARY
   
   prisav = andf(rxmsk,FILE(n))
   IF (iexec >= 0) GO TO 550
   
!     IF PRIMARY FILE IS PURGED OR HAS ZERO TRAILERS, PURGE SECONDARYS
   
   IF (andf(rmsk,prisav) == zap) GO TO 300
   IF (fmat(n) /= 0 .OR.  fmat(n+1) /= 0 .OR. fmat(n+2) /= 0) GO TO 405
   IF (entn1 == 11 .AND. (fmat(n+5) /= 0 .OR. fmat(n+6) /= 0 .OR.  &
       fmat(n+7) /= 0)) GO TO 405
   GO TO 300
   405 fequ(n) = orf(s,fequ(n))
   nsav = n
   CALL xpolck (xj1,xj2,fnsav,lsav)
   IF (fnsav /= 0) dfnu(lsav) = orf(s,dfnu(lsav))
   
!     IF PRIMARY FILE CONTAINS OTHER UNEQUIV D.B.- CLEAR THEM
   
   DO  jx = 1,lmt3,entn1
     IF (fequ(jx) < 0 .OR. jx == n) CYCLE
     IF (prisav /= andf(rxmsk,FILE(jx))) CYCLE
     lmt4 = jx + entn1 - 1
     DO  m = jx,lmt4
       FILE(m) = 0
     END DO
     IF (jx <= lmt1) FILE(jx) = prisav
   END DO
   GO TO 550
   
!     SECONDARY
   
   420 IF (iexec >= 0) GO TO 425
   IF (prisav == 0) GO TO 435
   IF (FILE(n) < 0 .AND. andf(rxmsk,FILE(n)) /= prisav)  &
       secchn = andf(rmsk,FILE(n))
   IF (n <= lmt1) GO TO 430
   FILE(n  ) = orf(andf(lxmsk,FILE(n)),prisav)
   fequ(n  ) = orf(s,fequ(n))
   fmat(n  ) = fmat(nsav  )
   fmat(n+1) = fmat(nsav+1)
   fmat(n+2) = fmat(nsav+2)
   IF (entn1 /= 11) GO TO 480
   fmat(n+5) = fmat(nsav+5)
   fmat(n+6) = fmat(nsav+6)
   fmat(n+7) = fmat(nsav+7)
   GO TO 480
   425 IF (fequ(n) >= 0 .OR. prisav /= andf(rxmsk,FILE(n))) GO TO 550
   
!     UNEQUIV (CLEAR SEC EQUIV ENTRY)
   
   lmt4 = n + entn1 - 1
   DO  m = n,lmt4
     FILE(m) = 0
   END DO
   IF (n <= lmt1) FILE(n) = prisav
   CALL xpolck (xj1,xj2,fn,l)
   IF (fn == 0) GO TO 550
   ddbn(l  ) = 0
   ddbn(l+1) = 0
   GO TO 550
   430 FILE(nfculg) = orf(andf(lxmsk,FILE(n)),prisav)
   435 hold = andf(rxmsk,FILE(n))
   lmt4 = n + entn1 -1
   DO  m = n,lmt4
     FILE(m) = 0
   END DO
   IF (n <= lmt1) FILE(n) = hold
   IF (prisav == 0) GO TO 480
   GO TO 470
   
!     FILE IS NOT IN FIAT -- CHECK PARM FOR TYPE OF EQUIV
   
 END DO
 IF (iexec < 0) GO TO 458
 
!     ELIMINATE EQUIV FILES -- CHECK FOR PRIMARY FILE
 
 IF (j /= i+1) GO TO 455
 
!     PRIMARY FILE
 
 CALL xpolck (xj1,xj2,fnsav,lsav)
 
!     LEAVE EQUIV FLAG FOR XDPH
 
 GO TO 550
 
!     SECONDARY FILE --  BREAK EQUIV
 
 455 CALL xpolck (xj1,xj2,snsav,ssav)
 
!     CHECK IF FILE EXISTS AND IS EQUIVED TO PRIMARY FILE
 
 IF (snsav == 0 .OR. fnsav /= snsav) GO TO 550
 ddbn( ssav  ) = 0
 ddbn( ssav+1) = 0
 GO TO 550
 
!     CHECK FOR PRIMARY FILE
 
 458 IF (j /= i+1) GO TO 460
 
!     -IF PRIMARY, IT MUST BE ON POOL
 
 CALL xpolck (xj1,xj2,fnsav,lsav)
 IF (fnsav == 0) GO TO 300
 dfnu(lsav) = orf(s,dfnu(lsav))
 GO TO 550
 
!     -IF SECONDARY, WAS PRIMARY IN FIAT
 
 460 IF (prisav == 0) GO TO 480
 
!     -PRIMARY WAS IN FIAT, SET UP SECONDARY IN FIAT
 
 FILE(nfculg  ) = prisav
 470 fequ(nfculg  ) = orf(s,fequ(nfculg))
 fdbn(nfculg  ) = xj1
 fdbn(nfculg+1) = xj2
 fmat(nfculg  ) = fmat(nsav  )
 fmat(nfculg+1) = fmat(nsav+1)
 fmat(nfculg+2) = fmat(nsav+2)
 IF (entn1 /= 11) GO TO 475
 fmat(nfculg+5) = fmat(nsav+5)
 fmat(nfculg+6) = fmat(nsav+6)
 fmat(nfculg+7) = fmat(nsav+7)
 475 fculg = fculg + incr
 IF (fculg >= fmxlg) GO TO 630
 nfculg = nfculg + entn1
 
!     WAS SECONDARY FILE IN FIAT ALREADY EQUIV TO OTHERS
 
 480 IF (secchn == 0) GO TO 490
 
!     SEC. FILE WAS EQUIV - DRAG ALONG ALL EQUIVS
 
 DO  ij = 1,lmt3,entn1
   IF (FILE(ij) >= 0) CYCLE
   IF (ij       == n) CYCLE
   IF (andf(rmsk,FILE(ij)) /= secchn) CYCLE
   
!     CREATE AN ENTRY IN OSCENT TO EXPLICITLY EQUIV THIS DB
   
   m1   = nwds + 1
   nwds = nwds + 6
   x(m1  ) = 2
   x(m1+1) = xj1
   x(m1+2) = xj2
   IF (k /= 1) GO TO 482
   x(m1+3) = 0
   nwds = nwds + 1
   m1   = m1   + 1
   482 CONTINUE
   x(m1+3) = fdbn(ij  )
   x(m1+4) = fdbn(ij+1)
   x(m1+5) = ivps
 END DO
 
!     IS SECONDARY FILE ON POOL
 
 490 CALL xpolck (xj1,xj2,fn,l)
 IF (fn == 0) GO TO 500
 
!     WAS SEC. FILE ON POOL ALREADY EQUIV TO OTHERS
 
 IF (dfnu(l) >= 0 .OR. fnsav == fn) GO TO 495
 
!     SEC. FILE ON POOL WAS EQUIV - DRAG ALONG ALL EQUIVS
 
 lmt4 = dculg* entn4
 m    = lmt4 + 1
 DO  ij = 1,lmt4,entn4
   IF (dfnu(ij) >= 0  .OR. ij == l) CYCLE
   IF (andf(rmsk,dfnu(ij)) /= fn) CYCLE
   IF (fnsav == 0) GO TO 491
   ddbn(m  ) = ddbn(ij  )
   ddbn(m+1) = ddbn(ij+1)
   dfnu(m  ) = dfnu(lsav)
   dculg = dculg + 1
   IF (dculg > dmxlg) GO TO 910
   m = m + entn4
   GO TO 493
   
!     CREATE AN ENTRY IN OSCENT TO EXPLICITLY EQUIV THIS DB
   
   491 m1   = nwds + 1
   nwds = nwds + 6
   x(m1  ) = 2
   x(m1+1) = xj1
   x(m1+2) = xj2
   IF (k /= 1) GO TO 492
   x(m1+3) = 0
   nwds = nwds + 1
   m1   = m1 + 1
   492 CONTINUE
   x(m1+3) = ddbn(ij  )
   x(m1+4) = ddbn(ij+1)
   x(m1+5) = ivps
   493 ddbn(ij  ) = 0
   ddbn(ij+1) = 0
 END DO
 
!     IF SECONDARY FILE IS ON POOL AND PRIMARY IS NOT - DELETE SEC REF
 
 495 IF (fnsav /= 0) GO TO 530
 ddbn(l  ) = 0
 ddbn(l+1) = 0
 GO TO 550
 
!     IF SECONDARY FILE IS NOT ON POOL AND PRIMARY IS - ADD SEC REF
 
 500 IF (fnsav == 0) GO TO 550
 520 m = dculg*entn4 + 1
 ddbn(m  ) = xj1
 ddbn(m+1) = xj2
 dfnu(m  ) = dfnu(lsav)
 dculg = dculg + 1
 IF (dculg > dmxlg) GO TO 910
 GO TO 550
 
!     BOTH PRIMARY AND SECONDARY ON POOL - IF NOT SAME FILE,
!     DELETE OLD SEC REF AND ADD NEW SEC REF
 
 530 IF (fnsav == fn) GO TO 550
 ddbn(l  ) = 0
 ddbn(l+1) = 0
 GO TO 520
 
 550 j = j + 2
 IF (j-2 == i+1) j = j + 1
 secchn = 0
 IF (j < jpt) GO TO 400
 prisav = 0
 GO TO 150
 
!     POTENTIAL FIAT OVERFLOW- LOOK FOR OTHER SLOTS IN FIAT TAIL
 
 610 ASSIGN 20 TO iback
 GO TO 640
 620 ASSIGN 280 TO iback
 GO TO 640
 630 ASSIGN 480 TO iback
 640 IF (fculg > fmxlg) GO TO 900
 DO  nn = lmt2,lmt3,entn1
   IF (FILE(nn) < 0 .OR. andf(zap,FILE(nn)) == zap) CYCLE
   IF (fmat(nn) /= 0 .OR. fmat(nn+1) /= 0 .OR. fmat(nn+2) /= 0) CYCLE
   IF (entn1 == 11 .AND. (fmat(nn+5) /= 0 .OR. fmat(nn+6) /= 0 .OR.  &
       fmat(nn+7) /= 0)) CYCLE
   nfculg = nn
   incr   = 0
   GO TO 660
 END DO
 nfculg = nfculg + entn1
 incr   = 1
 660 GO TO iback, (20,280,480)
 
!     ERROR MESSAGES
 
 900 WRITE  (outtap,901) sfm
 901 FORMAT (a25,' 1201, FIAT OVERFLOW.')
 GO TO 1000
 910 WRITE  (outtap,911) sfm
 911 FORMAT (a25,' 1202, DPL OVERFLOW.')
 1000 CALL mesage (-37,0,purge1)
 RETURN
END SUBROUTINE xpurge
