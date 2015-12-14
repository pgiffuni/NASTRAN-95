SUBROUTINE dpd2
     
!     DPD2 ASSEMBLES THE DYNAMIC LOADS TABLE (DLT).
 
 INTEGER :: gpl   ,sil   ,uset  ,usetd ,gpld  ,sild  ,dpool ,  &
     dlt   ,frl   ,tfl   ,trl   ,psdl  ,eed   ,scr1  ,  &
     scr2  ,scr3  ,scr4  ,buf   ,buf1  ,buf2  ,buf3  ,  &
     buf4  ,flag  ,FILE  ,epoint,seqep ,z     ,loads ,  &
     sdt   ,dload ,freq1 ,freq  ,tic   ,tstep ,tf    ,  &
     psd   ,eigr  ,eigb  ,eigc  ,ngrid ,eqdyn ,scr   , bufx
 DIMENSION       buf(24)      ,epoint(2)    ,seqep(2)     ,mcb(7),  &
     nam(2),loads(32)    ,dload(2)     ,freq1(2)     ,  &
     freq(2)      ,zz(1) ,bufr(20)     ,nolin(21)    ,  &
     tic(2),tstep(2)     ,tf(2) ,psd(2),msg(3),eigr(2)  &
     ,               eigb(2)      ,eigc(2)      ,scr(4),bufx(3)
 COMMON /BLANK / luset ,lusetd,notfl ,nodlt ,nopsdl,nofrl ,nonlft,  &
     notrl ,noeed ,nosdt ,noue
 COMMON /system/ idummy(55)   ,ithrml
 COMMON /names / rd    ,rdrew ,wrt   ,wrtrew,clsrew
 COMMON /dpdcom/ dpool ,gpl   ,sil   ,uset  ,gpld  ,sild  ,usetd ,  &
     dlt   ,frl   ,nlft  ,tfl   ,trl   ,psdl  ,eed   ,  &
     scr1  ,scr2  ,scr3  ,scr4  ,buf   ,buf1  ,buf2  ,  &
     buf3  ,buf4  ,epoint,seqep ,l     ,kn    ,neqdyn,  &
     loads ,dload ,freq1 ,freq  ,nolin ,nogo  ,msg   ,  &
     tic   ,tstep ,tf    ,psd   ,eigr  ,eigb  ,eigc  ,  &
     mcb   ,nam   ,eqdyn ,sdt   ,ineq
 COMMON /zzzzzz/ z(1)
 EQUIVALENCE     (z(1),zz(1)), (buf(1),bufr(1)), (msg(2),ngrid),  &
     (scr1,scr(1)),(buf2 ,bufx(1))
 
!     OPEN DYNAMICS POOL. SET POINTERS TO LOOP THRU DAREA, DELAY
!     AND DPHASE TABLES.
 
 FILE = dpool
 CALL preloc (*2001,z(buf1),dpool)
 ii  = 1
 iii = 1
 itabl = neqdyn + 2
 l   = 2
 j   = buf4 - 1
 msg(1) = 66
 
!     LOCATE CARD TYPE. IF PRESENT--
!     STORE POINTER TO 1ST TABLE NO. IN LOADS TABLE, OPEN SCRATCH FILE
!     FOR TABLES, SET ID = 0.
 
 1110 CALL locate (*1141,z(buf1),loads(ii),flag)
 loads(ii+2) = j
 FILE = scr(iii)
 CALL OPEN (*2001,FILE,z(buf2),wrtrew)
 id   = 0
 
!     READ A CARD. IF TABLE NO. IS DIFFERENT, STORE TABLE NO. IN TABLE
!     LIST.  IF NOT FIRST CARD, SORT TABLE ON SIL NO. AND WRITE ON
!     SCRATCH FILE.
 
 1120 CALL READ (*2002,*1140,dpool,buf,4,0,flag)
 IF (buf(1) == id) GO TO 1130
 IF (id == 0) GO TO 1122
 n = i - itabl
 CALL sort (0,0,2,1,z(itabl),n)
 CALL WRITE (FILE,z(itabl),n,1)
 1122 id = buf(1)
 z(j) = id
 j = j - 1
 i = itabl
 msg(3) = id
 
!     CONVERT POINT AND COMPONENT TO SIL NO.
!     STORE SIL NO. AND VALUE IN CORE.
 
 1130 CALL dpdaa
 z(i  ) = buf(2)
 z(i+1) = buf(4)
 i = i + 2
 IF (i < j) GO TO 1120
 CALL mesage (-8,0,nam)
 
!     HERE WHEN LAST CARD OF CURRENT TYPE HAS BEEN READ--
!     SORT AND WRITE LAST RECORD. CLOSE SCRATCH FILE.  STORE
!     NUMBER OF TABLES IN TABLE LIST. TEST FOR ALL CARD TYPES PROCESSED.
 
 1140 n = i - itabl
 CALL sort (0,0,2,1,z(itabl),n)
 CALL WRITE (FILE,z(itabl),n,1)
 CALL CLOSE (FILE,clsrew)
 loads(ii+3) = loads(ii+2) - j
 1141 ii  = ii  + 4
 iii = iii + 1
 IF (iii <= 3) GO TO 1110
 
!     SET POINTERS TO LOOP THRU RLOAD1,2 AND TLOAD1,2 CARDS
 
 ncore = j
 j     = 1
 iii   = 1
 ineq  = 0
 
!     LOCATE A CARD TYPE. IF PRESENT--
!     READ ALL CARDS OF TYPE INTO CORE.
 
 1160 CALL locate (*1165,z(buf1),loads(ii),flag)
 m = loads(ii+2)
 1161 z(j) = iii
 CALL READ (*2002,*1165,dpool,z(j+1),m,0,flag)
 j = j + 11
 IF (j < ncore) GO TO 1161
 CALL mesage (-8,0,nam)
 
!     TEST FOR ALL CARD TYPES PROCESSED.
!     IF SO, SORT CARDS ON LOAD SET ID.
 
 1165 ii  = ii  + 4
 iii = iii + 1
 IF (iii <= 4) GO TO 1160
 n   = j - 1
 IF (n /= 0) GO TO 1166
 CALL CLOSE (dpool,clsrew)
 RETURN
 
 1166 CALL sort (0,0,11,2,z,n)
 nlist = j - 11
 
!     LOCATE DLOAD CARDS ON DYNAMICS POOL.
!     IF PRESENT READ INTO CORE. SORT EACH DLOAD CARD ON REFERENCED SET
!     ID.
 
 nodld  = 0
 CALL locate (*1174,z(buf1),dload,flag)
 idload = j
 i = idload
 j = i
 1171 CALL READ (*2002,*1174,dpool,z(j+1),2,0,flag)
 j = j + 3
 nodld = nodld + 1
 1172 CALL READ (*2002,*2003,dpool,z(j),2,0,flag)
 IF (z(j) == -1) GO TO 1173
 j = j + 2
 IF (j >= ncore) CALL mesage (-8,0,nam)
 GO TO 1172
 1173 n = j - (i+3)
 CALL sort (0,0,2,2,z(i+3),n)
 
!     CHECK FOR DLOAD SET ID UNIQUENESS
 
 DO  kk = 2,n,2
   jj = i + 2 + kk
   IF (kk >= n) CYCLE
   IF (z(jj) /= z(jj+2)) CYCLE
   nogo = 1
   msg(2) = z(i+1)
   msg(3) = z(jj)
   CALL mesage (30,135,msg(2))
 END DO
 z(i) = n + 2
 i = j
 GO TO 1171
 1174 CALL CLOSE (dpool,clsrew)
 
!     OPEN THE DLT. WRITE NAME IN HEADER RECORD.
!     THEN WRITE NO. OF DLOAD CARDS FOLLOWED BY DLOAD SET IDS.
!     THEN WRITE SET IDS FOR EACH RECORD OF THE DLT (FOLLOWING DLOAD
!     RECORD)
 
 FILE = dlt
 CALL OPEN (*1249,dlt,z(buf1),wrtrew)
 CALL fname (dlt,buf)
 buf(3) = nodld
 CALL WRITE (dlt,buf,3,0)
 IF (nodld == 0) GO TO 1182
 i = idload
 j = 1
 1181 CALL WRITE (dlt,z(i+1),1,0)
 i = i + z(i) + 1
 j = j + 1
 IF (j <= nodld) GO TO 1181
 
!     CHECK DLOAD SID  VS  RLOAD1,2 AND TLOAD1,2 FOR UNIQUENESS
 
 i = idload
 DO  jj = 1,nodld
   itemp = z(i+1)
   DO  kk = 1,nlist,11
     IF (itemp /= z(kk+1)) CYCLE
     nogo = 1
     msg(2) = itemp
     CALL mesage (30,136,msg(2))
   END DO
   i = i + z(i) + 1
 END DO
 1182 DO  i = 1,nlist,11
   buf(1) = z(i+1)
   
!     CHECK FOR UNIQUE SET IDS ON TLOAD1,2 AND RLOAD1,2 CARDS  THEN WRIT
   
   IF (i >= nlist) GO TO 1184
   IF (z(i+1) /= z(i+12)) GO TO 1184
   nogo = 1
   msg(2) = itemp
   CALL mesage (30,136,msg(2))
   1184 CALL WRITE (dlt,buf,1,0)
 END DO
 CALL WRITE (dlt,0,0,1)
 
!     IF DLOAD CARDS PRESENT, WRITE THE DLOAD RECORD.
 
 IF (nodld == 0) GO TO 1200
 buf(1) = -1
 buf(2) = -1
 i = idload
 j = 1
 1191 n = z(i)
 CALL WRITE (dlt,z(i+1),n,0)
 CALL WRITE (dlt,buf,2,0)
 i = i + n + 1
 j = j + 1
 IF (j <= nodld) GO TO 1191
 CALL WRITE (dlt,0,0,1)
 
!     INITIALIZE TO LOOP THRU ALL LOAD SETS. THE REMAINDER OF THE DLT
!     WILL CONSIST OF ONE LOGICAL RECORD PER LOAD SET.
 
 1200 i = 1
 
!     WRITE FIXED SECTION OF DLT RECORD.
 
 1205 buf(1) = z(i  )
 buf(2) = z(i+2)
 
!     SAVE INFORCED MOTION FLAG ON TLOAD CARDS
 
 IF (z(i) < 3 .OR. z(i) > 4) GO TO 1206
 iemf = z(i+4)
 z(i+4) = 0
 1206 CONTINUE
 CALL WRITE (dlt,buf,2,0)
 CALL WRITE (dlt,z(i+5),6,0)
 
!     POSITION SCRATCH FILES TO SELECTED TABLES.
 
 idarea = 0
 DO  j = 1,3
   buf(2*j-1) = 16777215
!                  16777215 =2**24 - 1
   k = i + j
   buf(j+16) = z(k+1)
   IF (buf(j+16) == 0) CYCLE
   jj = loads(4*j-1)
   nn = loads(4*j  )
   IF (nn == 0) GO TO 1212
   DO  nx = 1,nn
     IF (z(jj) == buf(j+16)) GO TO 1213
     jj = jj - 1
   END DO
   1212 IF (ithrml /= 1 .OR. j /= 1) GO TO 1300
   idarea  = -1
   buf(17) = 0
   CYCLE
   1300 buf(10) = z(i+1)
   buf(11) = buf(j+16)
   buf(11) = buf(11) + 100000000*j
   nogo    = 1
   CALL mesage (30,71,buf(10))
   buf(j+16) = 0
   CYCLE
   1213 nn   = nx - 1
   FILE = scr(j)
   ibuf = bufx(j)
   CALL OPEN (*2001,FILE,z(ibuf),rdrew)
   IF (nn == 0) CYCLE
   DO  nx = 1,nn
     CALL fwdrec (*2002,FILE)
   END DO
 END DO
 
!     INITIALIZE TABLE READ.
 
 buf(14) = buf(17)
 buf(15) = buf(18)
 buf(16) = buf(19)
 
!     READ AN ENTRY FROM APPROPRIATE TABLE/S).
!     IF ALL ENTRIES HAVE BEEN READ, GO TO CLOSE DLT RECORD.
 
 1220 DO  j = 1,3
   IF (ithrml /= 1 .OR. j /= 1) GO TO 1320
   IF (idarea == 0) GO TO 1320
   IF (idarea == -2) GO TO 1221
   idarea  = -2
   buf(1)  = 1
   buf(2)  = 0
   buf(14) = 0
   1320 IF (buf(j+13) == 0) CYCLE
   FILE = scr(j)
   j2   = 2*j
   CALL READ (*2002,*1221,FILE,buf(j2-1),2,0,flag)
   CYCLE
   1221 buf(2*j-1) = 16777215
   buf(j+13)  = 0
 END DO
 IF (buf(1)+buf(3)+buf(5) == 3*16777215) GO TO 1240
 
!     SELECT MINIMUM SIL NO(S) AND FORMAT OUTPUT.
 
 DO  j = 1,6
   buf(j+10) = 0
 END DO
 buf(7) = 1
 buf(8) = 2
 buf(9) = 3
 IF (buf(1) > buf(3)) GO TO 1232
 
!     1 .LE. 2--COMPARE 2 TO 3. IF 2 .GT. 3, SWITCH 2 AND 3.
 
 IF (buf(3) <= buf(5)) GO TO 1234
 k = buf(8)
 buf(8) = buf(9)
 buf(9) = k
 GO TO 1233
 
!     1 .GT. 2--SWITCH 1 AND 2 THEN COMPARE 2 AND 3. IF 2 .GT. 3, SWITCH
 
 1232 k = buf(7)
 buf(7) = buf(8)
 buf(8) = k
 IF (buf(1) <= buf(5)) GO TO 1234
 k = buf(8)
 buf(8) = buf(9)
 buf(9) = k
 
!     COMPARE 1 TO 2--IF 1 .GT. 2, SWITCH 1 AND 2.
 
 1233 k = buf(7)
 l = buf(8)
 IF (buf(2*k-1) <= buf(2*l-1)) GO TO 1234
 buf(7) = l
 buf(8) = k
 
!     PICK UP 1. SET TO READ 1.
 
 1234 k = buf(7)
 buf(  10) = buf(2*k-1)
 buf(k+10) = buf(2*k  )
 buf(k+13) = k
 
!     IF 1 .EQ. 2, PICK UP 2 AND SET TO READ 2.
 
 l = buf(8)
 IF (buf(2*k-1) /= buf(2*l-1)) GO TO 1235
 buf(l+10) = buf(2*l)
 buf(l+13) = l
 
!     IF 1 .EQ. 2 .EQ. 3, PICK UP 3 AND SET TO READ 3.
 
 m = buf(9)
 IF (buf(2*l-1) /= buf(2*m-1)) GO TO 1235
 buf(m+10) = buf(2*m)
 buf(m+13) = m
 
!     WRITE SIL NO., A, TAU, THETA. THEN GO TO READ ANOTHER TABLE
!     ENTRY(S).
 
 1235 IF (z(i) < 3 .OR. z(i) > 4) GO TO 1236
 buf(13) = iemf
 1236 CALL WRITE (dlt,buf(10),4,0)
 GO TO 1220
 
!     CLOSE DLT RECORD,  CLOSE TABLES AND TEST FOR COMPLETION OF DLT.
 
 1240 CALL WRITE (dlt,0,0,1)
 DO  j = 1,3
   IF (buf(j+16) /= 0) CALL CLOSE (scr(j),clsrew)
 END DO
 i = i + 11
 IF (i <= nlist) GO TO 1205
 
!     CLOSE DLT, WRITE TRAILER AND RETURN.
 
 CALL CLOSE (dlt,clsrew)
 mcb(1) = dlt
 mcb(2) = dlt
 CALL wrttrl (mcb)
 nodlt = 1
 1249 RETURN
 
!     FATAL FILE ERRORS
 
 2001 n= -1
 GO TO 2005
 2002 n= -2
 GO TO 2005
 2003 n= -3
 2005 CALL mesage (n,FILE,nam)
 RETURN
END SUBROUTINE dpd2
