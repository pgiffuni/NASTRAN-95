SUBROUTINE dpd5
     
!     DPD5 ASSEMBLEMS
!     (1) THE EIGENVALUE EXTRACTION DATA BLOCK (EED), AND
!     (2) THE TRANSFER FUNCTION LIST (TFL).
 
!     REVISED  9/1989, BY G.C./UNISYS
!     NO COLUMN AND ROW WORD PACKING IN TFL FILE FOR MACHINES WITH 32
!     BIT WORD SIZE, OR LESS
 
 EXTERNAL        andf  ,orf   ,lshift
 LOGICAL :: pack
 INTEGER :: gpl   ,sil   ,uset  ,usetd ,gpld  ,sild  ,dpool  ,  &
     dlt   ,frl   ,tfl   ,trl   ,psdl  ,eed   ,scr1   ,  &
     scr2  ,scr3  ,scr4  ,buf   ,buf1  ,buf2  ,buf3   ,  &
     buf4  ,flag  ,FILE  ,epoint,seqep ,z     ,loads  ,  &
     orf   ,dload ,freq1 ,freq  ,tic   ,tstep ,tf     ,  &
     psd   ,eigr  ,eigb  ,eigc  ,ngrid ,eqdyn ,sdt    ,  &
     ua    ,ud    ,eigp  ,andf  ,two
 INTEGER :: imsg(2)
 DIMENSION       buf(24)   ,epoint(2)    ,seqep(2)     ,mcb(7)    ,  &
     nam(2)    ,loads(32)    ,dload(2)     ,freq1(2)  ,  &
     freq(2)   ,zz(1)        ,bufr(20)     ,nolin(21) ,  &
     tic(2)    ,tstep(2)     ,tf(2)        ,psd(2)    ,  &
     msg(3)    ,eigr(2)      ,eigb(2)      ,eigc(2)   , eigp(2)
 COMMON /machin/ mach  ,ihalf ,jhalf
 COMMON /BLANK / luset ,lusetd,notfl ,nodlt ,nopsdl,nofrl ,nonlft ,  &
     notrl ,noeed ,nosdt ,noue
 COMMON /names / rd    ,rdrew ,wrt   ,wrtrew,clsrew
 COMMON /dpdcom/ dpool ,gpl   ,sil   ,uset  ,gpld  ,sild  ,usetd  ,  &
     dlt   ,frl   ,nlft  ,tfl   ,trl   ,psdl  ,eed    ,  &
     scr1  ,scr2  ,scr3  ,scr4  ,buf   ,buf1  ,buf2   ,  &
     buf3  ,buf4  ,epoint,seqep ,l     ,kn    ,neqdyn ,  &
     loads ,dload ,freq1 ,freq  ,nolin ,nogo  ,  &
     msg   ,tic   ,tstep ,tf    ,psd   ,eigr  ,eigb   ,  &
     eigc  ,mcb   ,nam   ,eqdyn ,sdt   ,ineq
 COMMON /bitpos/ um    ,uo    ,ur    ,usg   ,usb   ,ul    ,ua     ,  &
     uf    ,us    ,un    ,ug    ,ue    ,up    ,une    , ufe   ,ud
 COMMON /setup / ifile(6)
 COMMON /zzzzzz/ z(1)
 COMMON /two   / two(32)
 COMMON /system/ ibuf  ,nout
 EQUIVALENCE     (z(1) ,zz(1)), (buf(1) ,bufr(1)), (msg(2),ngrid)
 DATA    eigp  / 257   ,4/
 
!     (1) PROCESS EDD
!     ===============
 
!     OPEN EED AND WRITE HEADER.
!     INITIALIZE TO LOOP THROUGH EIG CARDS.
 
!     OPEN DYNAMICS POOL.
 
 FILE = dpool
 CALL preloc (*310,z(buf1),dpool)
 
 FILE = eed
 CALL OPEN (*170,eed,z(buf2),wrtrew)
 FILE = dpool
 CALL fname (eed,buf)
 CALL WRITE (eed,buf,2,1)
 in = 0
 DO  j = 2,7
   mcb(j) = 0
 END DO
 l = 12
 msg(1) = 75
 
!     LOCATE EIGB CARDS IN DYNAMICS POOL. IF PRESENT, TURN NOEED FLAG
!     OFF, WRITE ID ON EED AND TURN ON TRAILER BIT.
 
 CALL locate (*30,z(buf1),eigb,flag)
 noeed = 1
 CALL WRITE (eed,eigb,2,0)
 CALL WRITE (eed,0,1,0)
 j = (eigb(2)-1)/16
 k =  eigb(2) - 16*j
 mcb(j+2) = orf(mcb(j+2),two(k+16))
 ASSIGN 23 TO nback
 l = 12
 mask = two(ua)
 
!     READ EIGB CARDS. IF GRID NO. IS PRESENT, CONVERT TO SIL VALUE.
!     WRITE DATA ON EED.
 
 22 CALL READ (*320,*24,dpool,buf,18,0,flag)
 GO TO 120
 23 CALL WRITE (eed,buf,12,0)
 CALL WRITE (eed,buf(14),6,0)
 GO TO 22
 24 CALL WRITE (eed,0,0,1)
 
!     LOCATE EIGC CARDS IN DYNAMICS POOL. IF PRESENT, TURN OFF NOEED
!     FLAG, WRITE ID ON EED AND TURN ON TRL BIT.
 
 30 CALL locate (*80,z(buf1),eigc,flag)
 noeed = 1
 CALL WRITE (eed,eigc,2,0)
 CALL WRITE (eed,0,1,0)
 j = (eigc(2)-1)/16
 k =  eigc(2) - 16*j
 mcb(j+2) = orf(mcb(j+2),two(k+16))
 ASSIGN 50 TO nback
 l = 6
 mask = two(ud)
 
!     READ EIGC CARDS. IF GRID NO. IS PRESENT, CONVERT TO SIL VALUE.
!     WRITE DATA ON EED.
 
 40 CALL READ (*320,*70,dpool,buf,10,0,flag)
 GO TO 120
 50 CALL WRITE (eed,buf,7,0)
 CALL WRITE (eed,buf(8),3,0)
 60 CALL READ  (*320,*320,dpool,buf,7,0,flag)
 CALL WRITE (eed,buf,7,0)
 IF (buf(1) /= -1) GO TO 60
 GO TO 40
 70 CALL WRITE (eed,0,0,1)
 
!     LOCATE EIGP CARDS. IF PRESENT, TURN NOEED FLAG OFF,
!     WRITE ID ON EED AND TURN ON TRAILER BIT. COPY DATA ON EED.
 
 80 CALL locate (*89,z(buf1),eigp,flag)
 noeed = 1
 CALL WRITE (eed,eigp,2,0)
 CALL WRITE (eed,0,1,0)
 j = (eigp(2)-1)/16
 k =  eigp(2) - 16*j
 mcb(j+2) = orf(mcb(j+2),two(k+16))
 81 CALL READ  (*320,*82,dpool,buf,4,0,flag)
 CALL WRITE (eed,buf,4,0)
 GO TO 81
 82 CALL WRITE (eed,0,0,1)
 
!     LOCATE EIGR CARDS IN DYNAMICS POOL. IF PRESENT, TURN OFF NOEED
!     FLAG, WRITE ID ON EED AND TURN ON TRL BIT.
 
 89 CALL locate (*160,z(buf1),eigr,flag)
 noeed = 1
 CALL WRITE (eed,eigr,2,0)
 CALL WRITE (eed,0,1,0)
 j = (eigr(2)-1)/16
 k =  eigr(2) - 16*j
 mcb(j+2) = orf(mcb(j+2),two(k+16))
 ASSIGN 100 TO nback
 l = 12
 mask = two(ua)
 
!     READ EIGR CARDS. IF GRID NO. IS PRESENT, CONVERT TO SIL VALUE.
!     WRITE DATA ON EED.
 
 90 CALL READ (*320,*110,dpool,buf,18,0,flag)
 GO TO 120
 100 CALL WRITE (eed,buf,12,0)
 CALL WRITE (eed,buf(14),6,0)
 GO TO 90
 110 CALL WRITE (eed,0,0,1)
 GO TO 160
 
!     CODE TO CONVERT GRID NO. AND COMPOIENT CODE TO SIL NO.
!     SIL NO. IS IN A SET FOR EIGR, EIGB - IN D SET FOR EIGC.
!     WRITE DATA ON EED.
 
 120 IF (buf(l) == 0) GO TO nback, (23,50,100)
 IF (in /= 0) GO TO 140
 FILE = usetd
 CALL OPEN (*310,usetd,z(buf3),rdrew)
 CALL fwdrec (*320,usetd)
 iusetd = neqdyn+2
 CALL READ (*320,*130,usetd,z(iusetd),buf3-iusetd,1,n)
 CALL mesage (-8,0,nam)
 130 CALL CLOSE (usetd,clsrew)
 in = 1
 140 imsg(1) = buf(1)
 imsg(2) = buf(l)
 CALL dpdaa
 nusetd = iusetd + buf(l) - 1
 buf(l) = 0
 DO  j = iusetd,nusetd
   IF (andf(z(j),mask) /= 0) buf(l)= buf(l) + 1
 END DO
 IF (andf(z(nusetd),mask) /= 0) GO TO nback, (23,50,100)
 nogo = 1
 CALL mesage (30,107,imsg)
 GO TO nback, (23,50,100)
 
!     COMPLETE EIG CARD PROCESSING.
 
 160 CONTINUE
 CALL CLOSE (eed,clsrew)
 mcb(1) = eed
 CALL wrttrl (mcb)
 
 
!     (2) PRECESS TFL FILE
!     ====================
 
!     SELECT PACK OR NO-PACK LOGIC
 
 170 pack = .true.
 i45 = 4
 i23 = 3
 IF (ihalf > 16) GO TO 175
 pack = .false.
 i45 = 5
 i23 = 2
 175 CONTINUE
 
!     OPEN TFL. WRITE HEADER. INITIALIZE TO READ TF CARDS.
 
 DO  j = 2,7
   mcb(j) = 0
 END DO
 CALL locate (*300,z(buf1),tf,flag)
 notfl = 0
 FILE  = tfl
 CALL OPEN (*300,tfl,z(buf2),wrtrew)
 CALL fname (tfl,buf)
 CALL WRITE (tfl,buf,2,1)
 msg(1) = 68
 l   = 2
 id  = 0
 itfl= neqdyn + 2
 i   = itfl
 isw = 0
 last= 0
 
!     READ FIXED SECTION OF TF CARD. CONVERT GRID POINT AND COMP. TO
!     SIL NO. TEST FOR NEW TRANSFER FUNCTION SET.
 
 190 CALL READ (*320,*200,dpool,buf,6,0,flag)
 msg(3) = buf(1)
 CALL dpdaa
 irow = buf(2)
 IF (buf(1) == id) GO TO 250
 notfl = notfl + 1
 IF (id /= 0) GO TO 210
 id = buf(1)
 GO TO 250
 
!     SORT TRANSFER EQUATIONS AND WRITE ON TFL ONE RECORD PER TRANSFER
!     FUNCTION SET. FIRST WORD OF RECORD IS SETID.
 
 200 last = 1
 210 CALL WRITE (tfl,id,1,0)
 IF (isw == 0) GO TO 220
 CALL WRITE (scr1,0,0,1)
 CALL CLOSE (scr1,clsrew)
 CALL OPEN  (*310,scr1,z(buf2),rdrew)
 ifile(1) = scr2
 ifile(2) = scr3
 ifile(3) = scr4
 n = buf3 - itfl
 IF (     pack) CALL sorti  (scr1,tfl,4,1,z(itfl),n)
 IF (.NOT.pack) CALL sorti2 (scr1,tfl,5,1,z(itfl),n)
 CALL CLOSE  (scr1,clsrew)
 GO TO 230
 220 n = i - itfl
 IF (     pack) CALL sorti  (0,0,4,1,z(itfl),n)
 IF (.NOT.pack) CALL sorti2 (0,0,5,1,z(itfl),n)
 CALL WRITE  (tfl,z(itfl),n,1)
 230 i  = itfl
 id = buf(1)
 isw= 0
 IF (last /= 0) GO TO 290
 GO TO 250
 
!     READ GRID POINT, COMP AND VALUES. CONVERT POINT AND COMP. TO SIL
!     NO. STORE IN CORE. IF SPILL, WRITE ON SCR1.
 
 240 CALL READ (*320,*310,dpool,buf(2),5,0,flag)
 IF (buf(2) == -1) GO TO 190
 CALL dpdaa
 
!     INTEGER PACKING LOGIC (MACHINES WITH 36 BITS WORDS, OR MORE) -
!     PACK COLN AND ROW INTO ONE WORD IF BOTH CAN BE STORED IN HALF WORD
!     THEN FOLLOWED BY 3 COEFFICIENTS, TOTALLY 4 WORDS
 
!     NON-INTEGER PACKING LOGIC (MACHINES WITH 32 BITS WORDS) -
!     THE COLUMN AND ROW ARE NOT PACKED, AND THEREFORE NOT BOUNED TO
!     65536 SIZE LIMIT. 1ST WORD IS COLUMN, 2ND WORD IS ROW, THEN
!     FOLLOWED BY 3 COEFFICIENTS, TOTALLY 5 WORDS
!     THE SUBROUTINE SORTI2 IS USED WHEN SORTING TFL BY 2 KEY WORDS
 
 250 IF (.NOT.pack) GO TO 252
 IF (buf(2) >= jhalf .OR. irow >= jhalf) GO TO 340
 buf(3) = lshift(buf(2),ihalf)
 buf(3) = orf(buf(3),irow)
 GO TO 255
 252 buf(3) = irow
 255 IF (isw /= 0) GO TO 280
 IF (i+i45 > buf3) GO TO 270
 DO  j = i23,6
   z(i) = buf(j)
   i = i + 1
 END DO
 GO TO 240
 270 isw = 1
 FILE= scr1
 CALL OPEN (*310,scr1,z(buf3),wrtrew)
 n = i - itfl
 CALL WRITE (scr1,z(itfl),n,0)
 280 CALL WRITE (scr1,buf(i23),i45,0)
 GO TO 240
 
!     HERE WHEN ALL TRANSFER FUNCTION SETS COMPLETE.
!     CLOSE FILE AND WRITE TRAILER.
 
 290 CALL CLOSE (tfl,clsrew)
 mcb(2) = notfl
 mcb(1) = tfl
 CALL wrttrl (mcb)
 300 CALL CLOSE  (dpool,clsrew)
 RETURN
 
!     FATAL ERRORS
 
 310 n = -1
 GO TO 330
 320 n = -2
 330 CALL mesage (n,FILE,nam)
 340 WRITE  (nout,350) ihalf,buf(2),irow
 350 FORMAT ('0*** COLUMN OR ROW SIL NO. EXCEEDS',i3,' BITS WORD ',  &
     'PACKING LIMIT',2I9)
 CALL mesage (-37,nam,nam)
 RETURN
END SUBROUTINE dpd5
