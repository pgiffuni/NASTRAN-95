SUBROUTINE xpunp
     
!     THIS SUBROUTINE POOLS AND UNPOOLS FILES AS PRESCRIBED BY XFIAT
 
 IMPLICIT INTEGER (a-z)
 EXTERNAL        lshift    ,andf    ,orf
 DIMENSION       header( 8),head( 2),npunp( 2),BLOCK(1000),  &
     ddbn  ( 1),dfnu( 1),fcum ( 1),fcus (   1),  &
     fdbn  ( 1),fequ( 1),FILE ( 1),fknd (   1),  &
     fmat  ( 1),fntu( 1),fpun ( 1),fon  (   1),  &
     ford  ( 1),minp( 1),mlsn ( 1),mout (   1),  &
     mscr  ( 1),sal ( 1),sdbn ( 1),sntu (   1), sord  ( 1)
 CHARACTER (LEN=25) :: sfm
 CHARACTER (LEN=29) :: uim
 CHARACTER (LEN=25) :: uwm
 CHARACTER (LEN=23) :: ufm
 COMMON /xmssg / ufm,uwm,uim,sfm
 COMMON /xfiat / fiat(7)
 COMMON /xpfist/ pfist
 COMMON /xfist / fist(2)
 COMMON /xdpl  / dpd(6)
 COMMON /system/ ibufsz,outtap
 COMMON /zzzzzz/ buf1(1)
 COMMON /xsfa1 / md(401),sos(1501),comm(20),xf1at(5)
 EQUIVALENCE               (dpd  (1),dnaf    ),(dpd  (2),dmxlg   ),  &
     (dpd  (3),dculg   ),(dpd  (4),ddbn (1)),(dpd  (6),dfnu (1)),  &
     (fiat (1),funlg   ),(fiat (2),fmxlg   ),(fiat (3),fculg   ),  &
     (fiat (4),fequ (1)),(fiat (4),FILE (1)),(fiat (4),ford (1)),  &
     (fiat (5),fdbn (1)),(fiat (7),fmat (1)),(md   (1),mlgn    ),  &
     (md   (2),mlsn (1)),(md   (3),minp (1)),(md   (4),mout (1)),  &
     (md   (5),mscr (1)),(sos  (1),slgn    ),(sos  (2),sdbn (1)),  &
     (sos  (4),sal  (1)),(sos  (4),sntu (1)),(sos  (4),sord (1)),  &
     (xf1at(1),fntu (1)),(xf1at(1),fon  (1)),(xf1at(2),fpun (1)),  &
     (xf1at(3),fcum (1)),(xf1at(4),fcus (1)),(xf1at(5),fknd (1))
 EQUIVALENCE               (comm (1),almsk   ),(comm (2),apndmk  ),  &
     (comm (3),cursno  ),(comm (4),entn1   ),(comm (5),entn2   ),  &
     (comm (6),entn3   ),(comm (7),entn4   ),(comm (8),flag    ),  &
     (comm (9),fnx     ),(comm(10),lmsk    ),(comm(11),lxmsk   ),  &
     (comm(12),macsft  ),(comm(13),rmsk    ),(comm(14),rxmsk   ),  &
     (comm(15),s       ),(comm(16),scornt  ),(comm(17),tapmsk  ),  &
     (comm(18),thcrmk  ),(comm(19),zap     )
 DATA  n/1000/  ,pool/4HPOOL/  ,entn5/2/,  npunp/4HXPUN,4HP   /
 
!     ENTRY SIZE NUMBERS,  1=FIAT, 4=DPD
 
 isw1 = 0
 isw2 = 0
 entn1x  = entn1- 1
 fist(2) = 1 + pfist
 
!     COMPUTE INDEX FOR DUMMY ENTRY IN FIST
 
 fstidx = fist(2)*2 + 1
 lmt3   = fculg* entn1
 
!     CHECK FOR ANY FILES TO POOL
 
 fist(fstidx) = 101
 loop360:  DO  i = 1,lmt3,entn1
   IF (fpun(i) >= 0) CYCLE loop360
   nn = andf(almsk,fpun(i))
   fpun(i)= 0
   IF (fmat(i) /= 0 .OR.  fmat(i+1) /= 0 .OR. fmat(i+2) /= 0) GO TO 105
   IF (entn1 == 11 .AND. (fmat(i+5) /= 0 .OR. fmat(i+6) /= 0 .OR.  &
       fmat(i+7) /= 0)) GO TO 105
   nn = 1
   GO TO 268
   105 CALL xpolck (fdbn(i),fdbn(i+1),fn,nx)
   IF (fn == 0) GO TO 110
   j = nx
   GO TO 268
   110 IF (isw1 /= 0) GO TO 220
   isw1 = 1
   CALL OPEN (*900,pool,buf1,2)
   CALL xfilps (dnaf)
   CALL CLOSE (pool,2)
   CALL OPEN (*900,pool,buf1,3)
   fnx = dnaf
   220 fist(fstidx+1) = i + entn5
   CALL OPEN (*900,101,buf1(ibufsz+1),0)
   ncnt = 0
   
!     WRITE SPECIAL FILE HEADER RECORD -- XPOOL DICT NAME    ( 2 WORDS )
!                                       + DATA BLOCK TRAILER ( 3 WORDS
!                                                         OR   6 WORDS )
   
   CALL WRITE (pool,fdbn(i),2,0)
   IF (entn1 == 11) GO TO 230
   CALL WRITE (pool,fmat(i),3,1)
   GO TO 240
   230 CALL WRITE (pool,fmat(i  ),3,0)
   CALL WRITE (pool,fmat(i+5),3,1)
   
!     READ AND WRITE 1ST 2 WORDS OF DATA BLOCK HEADER.
!     THEN CALL CPYFIL TO COPY REMAINDER OF FILE.
   
   240 CALL READ (*910,*920,101,head,2,0,flag)
   CALL WRITE (pool,head,2,0)
   CALL cpyfil (101,pool,BLOCK,n,flag)
   ncnt = andf(lxmsk,lshift(flag/1000+1,16))
   CALL eof (pool)
   CALL CLOSE (101,1)
   
!     ADD FILE NAME OF FILE JUST POOLED TO DPD
   
   j = dculg* entn4+ 1
   dculg = dculg+ 1
   IF (dculg > dmxlg) GO TO 700
   dfnu(j  ) = orf(dnaf,ncnt)
   ddbn(j  ) = fdbn(i  )
   ddbn(j+1) = fdbn(i+1)
   CALL sswtch (3,l)
   IF (l /= 1)  GO TO 267
   CALL page2 (-2)
   WRITE  (outtap,266) ddbn(j),ddbn(j+1),head(1),head(2)
   266 FORMAT (16H0POOL FILE NAME , 2A4, 17H DATA BLOCK NAME ,2A4)
   267 dnaf = dnaf + 1
   fnx  = fnx  + 1
   268 hold = andf(rxmsk,FILE(i))
   lmt4 = i + entn1x
   DO  kk = i,lmt4
     FILE(kk) = 0
   END DO
   FILE(i)  = hold
   fdbn(i)  = almsk
   
!     CHECK FOR EQUIV FILES
   
   IF (nn == 1) CYCLE loop360
   
!     THERE ARE EQUIV FILES
   
   dfnu(j) = orf(s,dfnu(j))
   dfnusv  = dfnu(j)
   DO  k = 1,lmt3,entn1
     IF (fequ(k) >= 0 .OR. i == k) CYCLE
     IF (andf(rmsk,FILE(i)) /= andf(rmsk,FILE(k))) CYCLE
     
!     THIS IS AN EQUIV FILE
     
     CALL xpolck (fdbn(k),fdbn(k+1),fn,nx)
     IF (fn == 0) GO TO 272
     IF (dfnu(nx) == dfnusv) GO TO 277
     ddbn(nx  ) = 0
     ddbn(nx+1) = 0
     272 j = j + entn4
     dculg = dculg + 1
     IF (dculg > dmxlg) GO TO 700
     dfnu(j  ) = dfnusv
     ddbn(j  ) = fdbn(k)
     ddbn(j+1) = fdbn(k+1)
     lmt4= k+ entn1x
     DO  kk = k,lmt4
       FILE(kk) =  0
     END DO
     277 nn = nn - 1
     IF (nn == 1) CYCLE loop360
   END DO
   GO TO 930
 END DO loop360
 IF (isw1 == 0) GO TO 400
 CALL CLOSE (pool,1)
 fnx = 1
 
!     CHECK FOR ANY FILES TO UNPOOL
 
 400 fist(fstidx) = 201
 405 fn = dnaf
 DO  i = 1,lmt3,entn1
   IF (fpun(i) <= 0 .OR. fpun(i) >= fn) CYCLE
   fn = fpun(i)
   ii = i
 END DO
 IF (fn == dnaf) GO TO 570
 fpun(ii) = 0
 IF (isw2 /= 0) GO TO 470
 isw2 = 1
 CALL OPEN (*900,pool,buf1,0)
 fnx = 1
 470 CALL xfilps (fn)
 fnx= fn
 fist(fstidx+1) = ii + entn5
 CALL OPEN (*900,201,buf1(ibufsz+1),1)
 
!     READ SPECIAL FILE HEADER RECORD AND, IF DIAG 3 IS ON, PRINT MSG
 
 CALL READ (*910,*920,pool,header,entn1-3,1,flag)
 CALL sswtch (3,l)
 IF (l /= 1) GO TO 500
 CALL page2 (-2)
 WRITE  (outtap,501) fdbn(ii),fdbn(ii+1),header(1),header(2)
 501 FORMAT (17H0XUNPL-dict NAME ,2A4, 16H pool FILE NAME , 2A4 )
 
!     COPY FILE USING CPYFIL
 
 500 CALL cpyfil (pool,201,BLOCK,n,flag)
 CALL CLOSE (201,1)
 fnx = fnx + 1
 fmat(ii  ) = header(3)
 fmat(ii+1) = header(4)
 fmat(ii+2) = header(5)
 IF (entn1 /= 11) GO TO 510
 fmat(ii+5) = header(6)
 fmat(ii+6) = header(7)
 fmat(ii+7) = header(8)
 
!     IS FILE EQUIVALENCED
 
 510 IF (fequ(ii) >= 0) GO TO 405
 
!     YES, COPY SAME TRAILER INTO ALL EQUIV FILES
 
 hold = andf(rmsk,FILE(ii))
 DO  j = 1,lmt3,entn1
   IF (fequ(j) >= 0 .OR. ii == j) CYCLE
   IF (hold /= andf(rmsk,FILE(j))) CYCLE
   fmat(j  ) = header(3)
   fmat(j+1) = header(4)
   fmat(j+2) = header(5)
   IF (entn1 /= 11) CYCLE
   fmat(j+5) = header(6)
   fmat(j+6) = header(7)
   fmat(j+7) = header(8)
 END DO
 GO TO 405
 570 IF (isw2 == 0) GO TO 600
 CALL CLOSE (pool,1)
 fnx = 1
 600 CONTINUE
 RETURN
 
 
 700 WRITE  (outtap,701)
 701 FORMAT (1H0,23X,19H 1031, dpl overflow)
 GO TO 1000
 900 WRITE  (outtap,901)
 901 FORMAT (1H0,23X,62H 1032, pool OR FILE being pooled/un-pooled could &
      &not be OPENED)
 GO TO 1000
 910 WRITE  (outtap,911)
 911 FORMAT (1H0,23X,39H 1033, illegal eof on FILE being pooled)
 GO TO 1000
 920 WRITE  (outtap,921)
 921 FORMAT (1H0,23X,39H 1034, illegal eor on FILE being pooled)
 GO TO 1000
 930 WRITE  (outtap,931)
 931 FORMAT (1H0,23X,33H 1035, equiv indicated,NONE found)
 1000 CALL page2 (-4)
 WRITE  (outtap,1001) sfm
 1001 FORMAT (a25,1H.)
 CALL mesage (-37,0,npunp)
 RETURN
END SUBROUTINE xpunp
