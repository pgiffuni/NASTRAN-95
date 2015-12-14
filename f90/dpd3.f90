SUBROUTINE dpd3
     
!     DPD3 ASSEMBLES THE FREQUENCY RESPONSE LIST (FRL)
!     AND THE POWER SPECTRAL DENSITY LIST (PSDL).
 
 INTEGER :: gpl   ,sil   ,uset  ,usetd ,gpld  ,sild  ,dpool ,  &
     dlt   ,frl   ,tfl   ,trl   ,psdl  ,eed   ,scr1  ,  &
     scr2  ,scr3  ,scr4  ,buf   ,buf1  ,buf2  ,buf3  ,  &
     buf4  ,flag  ,FILE  ,epoint,seqep ,z     ,loads ,  &
     randt2,dload ,freq1 ,freq  ,tic   ,tstep ,tf    ,  &
     psd   ,eigr  ,eigb  ,eigc  ,ngrid ,eqdyn ,sdt   , freq2 ,randps,randt1
 DIMENSION       buf(24)   ,epoint(2)    ,seqep(2)     ,mcb(7)   ,  &
     nam(2)    ,loads(32)    ,dload(2)     ,freq1(2) ,  &
     freq(2)   ,zz(1)        ,bufr(20)     ,nolin(21),  &
     tic(2)    ,tstep(2)     ,tf(2)        ,psd(2)   ,  &
     msg(3)    ,eigr(2)      ,eigb(2)      ,eigc(2)  ,  &
     freq2(2)  ,randps(2)    ,randt1(2)    ,randt2(2)
 COMMON /condas/ consts(5)
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
 COMMON /zzzzzz/ z(1)
 EQUIVALENCE     (consts(2),twopi), (z(1),zz(1)), (buf(1),bufr(1)),  &
     (msg(2),ngrid)
 DATA    freq2 , randps , randt1 , randt2  /  &
     1107,11, 2107,21, 2207,22, 2307,23 /
 
!     OPEN DYNAMICS POOL. SET POINTERS.
 
 FILE = dpool
 CALL preloc (*2001,z(buf1),dpool)
 nofrq1 = 0
 nofrq2 = 0
 nofrq  = 0
 ifrq1  = 1
 ifrq2  = ifrq1
 ifrq   = ifrq1
 i = ifrq1
 j = i
 
!     READ FREQ1 CARDS. CONVERT F1 AND DELTA F TO RADIANS.
 
 CALL locate (*1265,z(buf1),freq1,flag)
 nofrq1 = 1
 1261 CALL READ (*2002,*1262,dpool,z(i),4,0,flag)
 zz(i+1) = twopi*zz(i+1)
 zz(i+2) = twopi*zz(i+2)
 i = i + 4
 GO TO 1261
 1262 nfrq1 = i - 4
 ifrq2 = i
 ifrq  = i
 j = i
 
!     READ FREQ2 CARDS. CONVERT FREQUENCIES TO RADIANS.
 
 1265 CALL locate (*1270,z(buf1),freq2,flag)
 nofrq2 = 1
 1266 CALL READ (*2002,*1267,dpool,z(i),4,0,flag)
 zz(i+1) = twopi*zz(i+1)
 zz(i+2) = twopi*zz(i+2)
 i = i + 4
 GO TO 1266
 1267 nfrq2 = i - 4
 ifrq  = i
 j = i
 
!     READ FREQ CARDS. CONVERT FREQUENCIES TO RADIANS.
 
 1270 CALL locate (*1274,z(buf1),freq,flag)
 nofrq = 1
 1271 CALL READ (*2002,*1274,dpool,z(j+1),1,0,flag)
 j = j + 2
 1272 CALL READ (*2002,*2003,dpool,z(j),1,0,flag)
 IF (z(j) == -1) GO TO 1273
 zz(j) = twopi*zz(j)
 j = j + 1
 GO TO 1272
 1273 z(i) = j - (i+1)
 i = j
 GO TO 1271
 
!     TEST FOR ANY FREQ TYPE CARDS.
 
 1274 nofrl = nofrq1 + nofrq2 + nofrq
 IF (nofrl /= 0) GO TO 1280
 GO TO 1276
 1275 ineq  = 0
 1276 nofrl =-1
 GO TO 1310
 
!     COLLECT LIST OF FREQUENCY SET IDS AND POINTERS TO CARDS.
!     SORT THIS LIST ON SET ID.
 
 1280 ilist = j + 1
 i = ilist
 IF (nofrq1 == 0) GO TO 1282
 
!     FOR FREQ1 SET STORE SET ID, POINTER TO SET, 0.
 
 DO  k = ifrq1,nfrq1,4
   z(i  ) = z(k)
   z(i+1) = k
   z(i+2) = 0
   i = i + 3
 END DO
 nlist = i - 3
 1282 IF (nofrq2 == 0) GO TO 1287
 
!     FOR FREQ2 SET STORE SET ID, POINTER TO SET, -1.
 
 DO  k = ifrq2,nfrq2,4
   z(i  ) = z(k)
   z(i+1) = k
   z(i+2) =-1
   i = i + 3
 END DO
 nlist = i - 3
 1287 IF (nofrq == 0) GO TO 1285
 
!     FOR FREQ SET STORE SET ID, POINTER TO SET, NO. OF WORDS IN SET.
 
 j = ifrq
 1283 n = z(j)
 IF (n == -1) GO TO 1284
 j = j + 1
 z(i  ) = z(j)
 z(i+1) = j
 z(i+2) = n
 i = i + 3
 j = j + n
 GO TO 1283
 1284 nlist = i - 3
 1285 n = i - ilist
 CALL sort (0,0,3,1,z(ilist),n)
 
!     OPEN THE FRL. WRITE NAME + SET IDS IN HEADER.
 
 FILE = frl
 CALL OPEN (*1275,frl,z(buf2),wrtrew)
 CALL fname (frl,buf)
 CALL WRITE (frl,buf,2,0)
 DO  i = ilist,nlist,3
   buf(1) = z(i)
   CALL WRITE (frl,buf,1,0)
 END DO
 CALL WRITE (frl,0,0,1)
 
!     WRITE THE FRL ONE RECORD PER FREQUENCY SET.
!     CONVERT FREQ1 SETS TO LOOK LIKE FREQ SETS.
!     CONVERT FREQ2 SETS TO LOOK LIKE FREQ SETS.
 
 DO  i = ilist,nlist,3
   j = z(i+1)
   n = z(i+2)
   IF (n < 0) THEN
     GO TO  1304
   ELSE IF (n == 0) THEN
     GO TO  1301
   END IF
   
!     FREQ SET ---  SORT FREQUENCY LIST AND DISCARD ANY DUPLICATES.
!     THEN WRITE FREQUENCIES ON THE FRL
   
   1303 n = n - 1
   IF (n == 1) GO TO 1307
   CALL sortf (0,0,1,1,z(j+1),n)
   j1 = j + 2
   jn = j + n
   ix = j + 1
   DO  jx = j1,jn
     IF (z(jx) == z(ix)) CYCLE
     ix = ix + 1
     z(ix) = z(jx)
   END DO
   n = ix - j
   1307 CALL WRITE (frl,z(j+1),n,1)
   CYCLE
   
!     FREQ1 SET-- FORM F = F0 + (I-1)*DELTA F, WHERE I = 1 THRU N+1.
   
   1301 f0   = zz(j+1)
   delf = zz(j+2)
   n  = z(j+3) + 1
   fi = 0.
   DO  k = 1,n
     f = f0 + fi*delf
     CALL WRITE (frl,f,1,0)
     fi = fi + 1.0
   END DO
   CALL WRITE (frl,0,0,1)
   CYCLE
   
!     FREQ2 SET-- FORM F = F0*10.0**((I-1)*DELTA)
!     WHERE DELTA = (LOG10(FE/F0))/N AND I = 1 THRU N+1.
   
   1304 f0 = zz(j+1)
   fe = zz(j+2)
   n  =  z(j+3)
   fn = n
   delta = (ALOG10(fe/f0))/fn
   fi = 0.
   n  = n + 1
   DO  k = 1,n
     f  = f0*10.0**(fi*delta)
     CALL WRITE (frl,f,1,0)
     fi = fi + 1.0
   END DO
   CALL WRITE (frl,0,0,1)
 END DO
 
!     CLOSE FRL AND WRITE TRAILER.
 
 mcb(1) = frl
 mcb(2) = (nlist-ilist)/3 + 1
 CALL wrttrl (mcb)
 CALL CLOSE (frl,clsrew)
 ineq = 0
 
!     OPEN PSDL. IF PURGED, BYPASS PSDL PROCESSING.
!     OTHERWISE, LOCATE RANDPS CARDS. IF ABSENT, BYPASS PSDL PROCESSING.
 
 1310 FILE = psdl
 CALL OPEN (*1381,psdl,z(buf2),wrtrew)
 CALL locate (*1381,z(buf1),randps,flag)
 
!     READ RANDPS CARDS INTO CORE.
 
 irps = 1
 FILE = dpool
 CALL READ (*2002,*1322,dpool,z(irps),buf2-irps,1,nrps)
 GO TO 2004
 1322 irt1 = irps + nrps
 irt2 = irt1
 i    = irt1
 j    = i
 nort1= 0
 nort2= 0
 
!     READ RANDT1 CARDS.
 
 CALL locate (*1340,z(buf1),randt1,flag)
 CALL READ (*2002,*1332,dpool,z(irt1),buf2-irt1,1,nort1)
 GO TO 2004
 1332 irt2 = irt1 + nort1
 nrt1 = irt2 - 4
 i = irt2
 j = i
 
!     READ RANDT2 CARDS.
 
 1340 CALL locate (*1350,z(buf1),randt2,flag)
 nort2 = 1
 1341 CALL READ (*2002,*1350,dpool,z(j+1),1,0,flag)
 j = j + 2
 1342 CALL READ (*2002,*2003,dpool,z(j),1,0,flag)
 IF (z(j) == -1) GO TO 1343
 j = j + 1
 IF (j < buf2) GO TO 1342
 GO TO 2004
 1343 z(i) = j - (i+1)
 i = j
 GO TO 1341
 
!     COLLECT LIST OF RANDT1 AND RANDT2 SET IDS AND POINTERS TO DATA.
 
 1350 nort = nort1 + nort2
 IF (nort == 0) GO TO 1360
 ilist = j + 1
 i = ilist
 IF (nort1 == 0) GO TO 1352
 
!     FOR RANDT1 SETS STORE SET ID, POINTER TO SET, 0.
 
 DO  k = irt1,nrt1,4
   z(i  ) = z(k)
   z(i+1) = k
   z(i+2) = 0
   i = i + 3
 END DO
 nlist = i - 3
 IF (i >  buf2) GO TO 2004
 1352 IF (nort2 == 0) GO TO 1355
 
!     FOR RANDT2 SETS STORE SET ID, POINTER TO SET, NO. OF WORDS IN SET.
 
 j = irt2
 1353 n = z(j)
 IF (n == -1) GO TO 1354
 z(i  ) = z(j)
 z(i+1) = j
 z(i+2) = n
 i = i + 3
 j = j + n
 IF (i < buf2) GO TO 1353
 GO TO 2004
 1354 nlist = i - 3
 
!     SORT LIST ON SET ID.
 
 1355 n = i - ilist
 CALL sort (0,0,3,1,z(ilist),n)
 
!     WRITE SET IDS FOR RANDT1 AND RANDT2 CARDS IN HEADER RECORD OF
!     PSDL. THEN WRITE RANDPS DATA AS FIRST RECORD OF PSDL.
 
 1360 CALL fname (psdl,buf)
 CALL WRITE (psdl,buf,2,0)
 IF (nort == 0) GO TO 1362
 DO  i = ilist,nlist,3
   CALL WRITE (psdl,z(i),1,0)
 END DO
 1362 CALL WRITE (psdl,0,0,1)
 CALL WRITE (psdl,z(irps),nrps,1)
 IF (nort == 0) GO TO 1380
 
!     WRITE ONE RECORD ON PSDL FOR EACH RANDT1 OR RANDT2 SET.
 
 DO  i = ilist,nlist,3
   j = z(i+1)
   n = z(i+2)
   IF (n == 0) GO TO 1372
   
!     RANDT2 SET--  SORT DATA AND DISCARD ANY DUPLICATES. THEN WRITE SET
   
   n = n - 1
   IF (n == 1) GO TO 1376
   CALL sortf (0,0,1,1,z(j+1),n)
   j1 = j + 2
   jn = j + n
   ix = j + 1
   DO  jx = j1,jn
     IF (z(jx) == z(ix)) CYCLE
     ix = ix + 1
     z(ix) = z(jx)
   END DO
   n = ix - j
   1376 CALL WRITE (psdl,z(j+1),n,1)
   CYCLE
   
!     RANDT1 SET-- WRITE TI = T0 + (I-1)*DELTA T, WHERE I = 1 THRU N+1.
   
   1372 n  = z(j+1)
   fn = n
   delt = (zz(j+3)-zz(j+2))/fn
   t0 = zz(j+2)
   fi = 0.
   n  = n + 1
   DO  k = 1,n
     ti = t0 + fi*delt
     CALL WRITE (psdl,ti,1,0)
     fi = fi + 1.0
   END DO
   CALL WRITE (psdl,0,0,1)
 END DO
 
!     CLOSE FILES, WRITE TRAILER AND EXIT.
 
 1380 mcb(1) = psdl
 mcb(2) = (nlist-ilist)/3 + 1
!      2147483647  = 2**31 - 1
 IF (nort == 0) mcb(2) = 2147483647
 CALL wrttrl (mcb)
 ineq  = 0
 nopsdl= 1
 1381 CALL CLOSE (dpool,clsrew)
 CALL CLOSE (psdl ,clsrew)
 RETURN
 
!     FATAL FILE ERRORS
 
 2001 n = -1
 GO TO 2005
 2002 n = -2
 GO TO 2005
 2003 n = -3
 GO TO 2005
 2004 n = -8
 2005 CALL mesage (n,FILE,nam)
 RETURN
END SUBROUTINE dpd3
