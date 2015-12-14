SUBROUTINE dpd4
     
!     DPD4 ASSEMBLES THE NON-LINEAR FORCING TABLE (NLFT)
!     AND THE TRANSIENT RESPONSE LIST (TRL).
 
 EXTERNAL        andf
 INTEGER :: gpl   ,sil   ,uset  ,usetd ,gpld  ,sild  ,dpool ,  &
     dlt   ,frl   ,tfl   ,trl   ,psdl  ,eed   ,scr1  ,  &
     scr2  ,scr3  ,scr4  ,buf   ,buf1  ,buf2  ,buf3  ,  &
     buf4  ,flag  ,FILE  ,epoint,seqep ,z     ,loads ,  &
     andf  ,dload ,freq1 ,freq  ,tic   ,tstep ,tf    ,  &
     psd   ,eigr  ,eigb  ,eigc  ,ngrid ,eqdyn ,sdt   , ud    ,ue    ,two
 DIMENSION       buf(24)   ,epoint(2)    ,seqep(2)     ,mcb(7)   ,  &
     nam(2)    ,loads(32)    ,dload(2)     ,freq1(2) ,  &
     freq(2)   ,zz(1)        ,bufr(20)     ,nolin(21),  &
     tic(2)    ,tstep(2)     ,tf(2)        ,psd(2)   ,  &
     msg(3)    ,eigr(2)      ,eigb(2)      ,eigc(2)
 COMMON /BLANK / luset ,lusetd,notfl ,nodlt ,nopsdl,nofrl ,nonlft,  &
     notrl ,noeed
 COMMON /names / rd    ,rdrew ,wrt   ,wrtrew,clsrew
 COMMON /dpdcom/ dpool ,gpl   ,sil   ,uset  ,gpld  ,sild  ,usetd ,  &
     dlt   ,frl   ,nlft  ,tfl   ,trl   ,psdl  ,eed   ,  &
     scr1  ,scr2  ,scr3  ,scr4  ,buf   ,buf1  ,buf2  ,  &
     buf3  ,buf4  ,epoint,seqep ,l     ,kn    ,neqdyn,  &
     loads ,dload ,freq1 ,freq  ,nolin ,nogo  ,  &
     msg   ,tic   ,tstep ,tf    ,psd   ,eigr  ,eigb  ,  &
     eigc  ,mcb   ,nam   ,eqdyn ,sdt   ,ineq
 COMMON /two   / two(32)
 COMMON /bitpos/ um    ,uo    ,ur    ,usg   ,usb   ,ul    ,ua    ,  &
     uf    ,us    ,un    ,ug    ,ue    ,up    ,une   , ufe   ,ud
 COMMON /zzzzzz/ z(1)
 EQUIVALENCE     (z(1) ,zz(1)),(buf(1),bufr(1)),(msg(2),ngrid)
 DATA            nolinr/ 7 /
 
!     INITIALIZE POINTERS. OPEN SCR1. OPEN DYNAMICS POOL.
 
 inolin = neqdyn + 2
 j = inolin
 mskud = two(ud)
 mskue = two(ue)
 incore = 0
 ii = 1
 i  = 1
 msg(1) = 67
 CALL preloc (*2001,z(buf1),dpool)
 CALL OPEN (*2001,scr1,z(buf2),wrtrew)
 ineq = 0
 
!     LOCATE NOLINI CARD. IF PRESENT, TURN NONLFT FLAG OFF,
 
 1320 CALL locate (*1358,z(buf1),nolin(i),flag)
 nonlft = 1
 nwdin  = nolin(i+2)
 
!     READ A NOLINI CARD. CONVERT POINTS ON CARD TO SIL NOS.
!     STORE DATA IN CORE. IF SPILL, WRITE ON SCRATCH FILE.
 
 1340 CALL READ (*2002,*1358,dpool,buf,nwdin,0,flag)
 msg(3) = 100000000*ii + buf(1)
 IF (ii >= 5) IF (ii-6)  1350,  1354,  1342
!                             NOLIN5,NFTUBE,NOLIN6
 iii = ii
 IF (buf(6) < 10) GO TO 1341
 iii = ii + 4
 buf(6) = buf(6) - 10
 1341 IF (ii /= 2) GO TO 1343
 IF (buf(8) < 10) GO TO 1343
 buf(8) = buf(8) - 10
 IF(iii == 2) iii = 10
 IF(iii == 6) iii = 9
 GO TO 1343
 1342 iii = 13
 IF (buf(6) < 10) GO TO 1343
 iii = 14
 buf(6) = buf(6) - 10
 1343 l = 2
 CALL dpdaa
 buf(3) = buf(2)
 l = 5
 CALL dpdaa
 l = 7
 IF (ii == 2) CALL dpdaa
 buf(6) = buf(7)
 buf(2) = iii
 1344 nn = 6
 1345 IF (incore  /= 0) GO TO 1348
 IF (j+nn >= buf2) GO TO 1347
 DO  k = 1,nn
   z(j) = buf(k)
   j = j + 1
 END DO
 GO TO 1340
 1347 CALL WRITE (scr1,z(inolin),j-inolin,0)
 incore = 1
 1348 CALL WRITE (scr1,buf,nn,0)
 GO TO 1340
 
!     SPECIAL HANDLING OF NOLIN5 CARD
!     CARD FORMAT AS RECEIVED FROM IFP
!        SID  AA   AB   FAB  EA/TEA  EB/TEB  ALPA/TALPA  ALPB/TALPB
!        GA1  GA2  GA3  GA4  GB1     GB2     GB3         GB4
 
!     WE CONVERT THIS CARD INTO THE FOLLOWING 6-WORD ENTRY FORMAT
 
!        SID  12  SILA1  AA          SILA2  AB
!        SID  12  SILA3  FAB         SIL4   0
!        SID  12  SILB1  EA/TEA      SILB2  EB/TEB
!        SID  12  SILB3  ALPA/TALPA  SILB4  ALPB/TALPB
 
 1350 l = 23
 kk= 16
 DO  k = 1,8
   buf(l+1) = 0
   buf(l  ) = buf(kk)
   IF (buf(l) /= 0) CALL dpdaa
   kk = kk - 1
   l  = l -2
 END DO
 buf(24) = buf( 8)
 buf(22) = buf( 7)
 buf(18) = buf( 6)
 buf(16) = buf( 5)
 buf(12) = 0
 buf(10) = buf( 4)
 buf( 6) = buf( 3)
 buf( 4) = buf( 2)
 buf( 3) = buf( 9)
 buf( 5) = buf(11)
 buf( 9) = buf(13)
 buf(11) = buf(15)
 buf(17) = buf(19)
 DO  k = 1,24,6
   buf(k  ) = buf(1)
   buf(k+1) = 12
 END DO
 nn = 24
 GO TO 1345
 
 1354 l = 7
 buf(7) = buf(2)
 buf(8) = 1
 CALL dpdaa
 buf(3) = buf(7)
 buf(7) = buf(3)
 buf(8) = 1
 CALL dpdaa
 buf(5) = buf(7)
 buf(6) = buf(5)
 buf(2) = 11
 msg(3) = buf(1)
 GO TO 1344
 
!     HERE WHEN ALL CARDS OF CURRENT TYPE HAVE BEEN READ.
!     TEST FOR ALL CARDS READ.
 
 1358 i  = i + 3
 ii = ii+ 1
 IF (ii <= nolinr) GO TO 1320
 CALL WRITE (scr1,0,0,1)
 CALL CLOSE (scr1,clsrew)
 IF (nonlft == -1) GO TO 1400
 
!     SORT THE DATA ON SET ID.
 
 IF (incore /= 0) GO TO 1362
 nnolin = j - 6
 n = j - inolin
 GO TO 1364
 1362 CALL OPEN (*2001,scr1,z(buf2),rdrew)
 CALL READ (*2002,*1363,scr1,z,buf1,1,n)
 CALL mesage (-8,0,nam)
 1363 CALL CLOSE (scr1,clsrew)
 inolin = 1
 nnolin = n - 5
 1364 CALL sort (0,0,6,1,z(inolin),n)
 
!     READ USETD INTO CORE.
 
 FILE = usetd
 CALL OPEN (*2001,usetd,z(buf2),rdrew)
 CALL fwdrec (*2002,usetd)
 iusetd = nnolin + 7
 CALL READ (*2002,*1365,usetd,z(iusetd),buf2-iusetd,1,n)
 CALL mesage (-8,0,nam)
 1365 CALL CLOSE (usetd,clsrew)
 
!     OPEN THE NLFT. WRITE SET IDS IN HEADER RECORD.
 
 FILE = nlft
 CALL OPEN (*1392,nlft,z(buf2),wrtrew)
 CALL fname (nlft,buf)
 CALL WRITE (nlft,buf,2,0)
 z(nnolin+6) = 0
 DO  i = inolin,nnolin,6
   IF (z(i+6) /= z(i)) CALL WRITE (nlft,z(i),1,0)
 END DO
 CALL WRITE (nlft,0,0,1)
 
!     WRITE ONE RECORD PER SET. WITHIN EACH SET, SORT DATA ON SIL NO.
!     CONVERT SIL NOS. TO SIL NOS. IN UD AND UE SETS
 
 i = inolin
 1381 j = i
 1382 IF (z(i+6) /= z(i)) GO TO 1383
 i = i + 6
 GO TO 1382
 1383 n = i + 6 - j
 
! ... THE FOLLOWING SORT WAS REMOVED DUE TO THE INSTALLATION OF NOLIN5
!     CALL SORT (0,0,6,3,Z(J),N)
 
!WKBR SPR94005 6/94   DO 1387 KC = J,I,6
 DO  k = j,i,6
   buf(1) = z(k+1)
   buf(2) = z(k+2)
   buf(4) = z(k+3)
   buf(5) = z(k+4)
   buf(8) = z(k+5)
   buf(9) = 0
   DO  kk = 2,8,3
     IF (kk >= 8 .AND. buf(1) /= 2 .AND. buf(1) /= 6 .AND. buf(1) /= 9  &
         .AND.  buf(1) /= 10.AND.kk == 8) CYCLE
     k1 = 0
     k2 = 0
     nusetd = iusetd + buf(kk) - 1
     IF (nusetd < iusetd) GO TO 1385
     DO  kkk = iusetd,nusetd
       buf(10) = z(kkk)
       IF (andf(buf(10),mskud) /= 0) k1 = k1 + 1
       IF (andf(buf(10),mskue) /= 0) k2 = k2 + 1
     END DO
     1385 buf(kk  ) = k1
     buf(kk+1) = k2
     IF (nusetd < iusetd) CYCLE
     IF (andf(buf(10),mskue) == 0) buf(kk+1) = 0
     IF (andf(buf(10),mskud) /= 0) CYCLE
     nogo = 1
     buf(1) = z(k)
     buf(2) = k1
     CALL mesage (30,93,buf)
   END DO
   buf(7) = buf(8)
   buf(8) = buf(9)
   CALL WRITE (nlft,buf,8,0)
 END DO
 CALL WRITE (nlft,0,0,1)
 i = i + 6
 IF (z(i) /= 0) GO TO 1381
 
!     CLOSE FILE AND WRITE TRAILER.
 
 CALL CLOSE (nlft,clsrew)
 mcb(1) = nlft
 mcb(2) = (nnolin-inolin)/6 + 1
 CALL wrttrl (mcb)
 IF (incore /= 0) ineq = 0
 GO TO 1400
 1392 nonlft =-1
 
!     LOCATE TIC CARDS IN DYNAMICS POOL.
 
 1400 notrl =-1
 notic = 0
 notstp= 0
 CALL locate (*1500,z(buf1),tic,flag)
 notrl = 1
 
!     OPEN SCR1. INITIALIZE TO READ TIC CARDS.
 
 FILE = scr1
 CALL OPEN (*2001,scr1,z(buf2),wrtrew)
 itic = neqdyn + 2
 nset = buf3 - 1
 j    = nset
 l    = 2
 msg(1) = 69
 id   = 0
 
!     READ A TIC CARD. IF SET ID IS DIFFERENT, STORE IT IN LIST.
!     IF NOT FIRST CARD, SORT DATA ON SIL NO. AND WRITE IT IN SCR1.
 
 1420 CALL READ (*2002,*1440,dpool,buf,5,0,flag)
 IF (buf(1) == id) GO TO 1430
 IF (id == 0) GO TO 1421
 n = i - itic
 CALL sort (0,0,3,1,z(itic),n)
 CALL WRITE (scr1,z(itic),n,1)
 1421 id = buf(1)
 z(j) = id
 j  = j - 1
 i  = itic
 msg(3) = id
 
!     CONVERT POINT AND COMPONENT TO SIL NO.
!     STORE SIL NO., UO, VO IN CORE.
 
 1430 CALL dpdaa
 z(i  ) = buf(2)
 z(i+1) = buf(4)
 z(i+2) = buf(5)
 i = i + 3
 IF (i < j) GO TO 1420
 CALL mesage (-8,0,nam)
 
!     HERE WHEN LAST CARD READ - SORT AND WRITE LAST RECORD.
 
 1440 n = i - itic
 CALL sort (0,0,3,1,z(itic),n)
 CALL WRITE (scr1,z(itic),n,1)
 CALL CLOSE (scr1,clsrew)
 iset = j + 1
 
!     OPEN TRL. WRITE SET IDS IN HEADER.
 
 FILE = trl
 CALL OPEN (*1493,trl,z(buf2),wrtrew)
 CALL fname (trl,buf)
 n = nset - iset + 1
 buf(3) = n
 notic  = n
 CALL WRITE (trl,buf,3,0)
 i  = iset
 j  = nset
 1451 id = z(j)
 z(j) = z(i)
 z(i) = id
 i  = i + 1
 j  = j - 1
 IF (i < j) GO TO 1451
 CALL WRITE (trl,z(iset),n,0)
 
!     READ USETD INTO CORE.
!     COMPUTE NO. OF POINTS UN UD SET. WRITE NO. AS LAST WORD OF HEADER.
 
 1460 FILE = usetd
 CALL OPEN (*2001,usetd,z(buf3),rdrew)
 CALL fwdrec (*2002,usetd)
 iusetd = 1
 ineq   = 0
 CALL READ (*2002,*1462,usetd,z(iusetd),buf3-iusetd,1,n)
 CALL mesage (-8,0,nam)
 1462 CALL CLOSE (usetd,clsrew)
 nusetd = iusetd + n - 1
 k = 0
 DO  i = iusetd,nusetd
   IF (andf(z(i),mskud) /= 0) k = k + 1
 END DO
 CALL WRITE (trl,k,1,1)
 IF (notic == 0) GO TO 1481
 
!     READ SCR1. CONVERT SIL NO. TO AN SIL NO. IN THE D-SET.
!     WRITE TRL ONE RECORD PER SET.
 
 FILE = scr1
 kset = iset
 CALL OPEN (*2001,scr1,z(buf3),rdrew)
 1475 k = 0
 ipoint = iusetd
 1471 CALL READ (*1474,*1473,scr1,buf,3,0,flag)
 nusetd = iusetd + buf(1) - 1
 DO  i = ipoint,nusetd
   IF (andf(z(i),mskud) /= 0) k = k + 1
 END DO
 buf(1) = k
 IF (andf(z(nusetd),mskud) /= 0) GO TO 1476
 nogo = 1
 CALL mesage (30,133,z(kset))
 1476 CALL WRITE (trl,buf,3,0)
 ipoint = nusetd + 1
 GO TO 1471
 1473 CALL WRITE (trl,0,0,1)
 kset = kset + 1
 GO TO 1475
 1474 CALL CLOSE (scr1,clsrew)
 
!     IF TSTEP CARDS PRESENT, COPY THEM ONTO TRL.
 
 CALL locate (*1490,z(buf1),tstep,flag)
 1481 CALL READ (*2002,*1483,dpool,buf,1,0,flag)
 notstp = notstp + 1
 CALL WRITE (trl,buf,1,0)
 1482 CALL READ (*2002,*2003,dpool,buf,3,0,flag)
 IF (buf(1) == -1) GO TO 1485
 CALL WRITE (trl,buf,3,0)
 GO TO 1482
 1485 CALL WRITE (trl,0,0,1)
 GO TO 1481
 1483 CONTINUE
 
!     CLOSE FILES AND WRITE TRAILER.
 
 1490 CALL CLOSE (trl,clsrew)
 mcb(1) = trl
 mcb(2) = notic
 mcb(3) = notstp
 CALL wrttrl (mcb)
 1492 CALL CLOSE (dpool,clsrew)
 RETURN
 
 1493 notrl = -1
 GO TO 1492
 
!     HERE IF NO TIC CARDS - LOCATE TSTEP CARDS IN DYNAMICS POOL.
!     IF ABSENT, RETURN. OTHERWISE OPEN TRL AND WRTIE HEADER.
 
 1500 CALL locate (*1492,z(buf1),tstep,flag)
 notrl = 1
 FILE  = trl
 CALL OPEN (*1493,trl,z(buf2),wrtrew)
 CALL fname (trl,buf)
 buf(3) = 0
 CALL WRITE (trl,buf,3,0)
 GO TO 1460
 
!     FATAL FILE ERRORS
 
 2001 n = -1
 GO TO 2005
 2002 n = -2
 GO TO 2005
 2003 n = -3
 2005 CALL mesage (n,FILE,nam)
 RETURN
END SUBROUTINE dpd4
