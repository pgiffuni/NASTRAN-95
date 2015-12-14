SUBROUTINE mtrxin
     
!     TWO CAPABILITIES -
 
!     (1) TO PROVIDE DIRECT INPUT MATRICES CAPABILITY, IN DYNAMIC RIGID
!         FORMATS, AND
!     (2) TO CONVERT DMIG TYPE MATRICES TO NASTRAN MATRIX FORMAT.
 
!     REVISED  1/90 BY G.CHAN/UNISYS
!     NO INTEGER ROW AND COLUMN PACKING FOR 32-BIT WORD MACHINE, AND
!     REAL AND COMPLEX DMIG MATRIX GENERATION FROM DMIG INPUT CARDS
!     WITH DOUBLE PRECISION DATA
 
 IMPLICIT INTEGER (a-z)
 EXTERNAL         orf,lshift,rshift,andf
 LOGICAL :: pack
 REAL :: x,alpha(4),beta(4),bufr(13)
 DOUBLE PRECISION :: xd(2),bufd(2)
 DIMENSION        buf(20),mcb(50),dmig(2),nam(2),BLOCK(81),db(13),  &
     filei(3),bufi(3),filea(7),fileb(7),filec(7)
 CHARACTER (LEN=23) :: ufm
 COMMON /xmssg /  ufm
 COMMON /machin/  mach,ihalf,jhalf
 COMMON /BLANK /  luset,nomat(3) /zzzzzz/  z(1)  &
     /system/  sysbuf,nout,xx(5),loadnn /TYPE  /  prc(2),nwds(4)  &
     /zblpkx/  x(4),row /names /  rd,rdrew,wrt,wrtrew,clsrew  &
     /setup /  ifile(6) /saddx /  nomats,nz,mcbs(67)
 EQUIVALENCE      (buf(1)  ,bufr(1) ), (x(1)    ,xd(1)   ),  &
     (filei(1),file1   ), (filei(2),file2   ),  &
     (bufi(1) ,buf3    ), (bufi(2) ,buf4    ),  &
     (nomat(1),nomat1  ), (nomat(2),nomat2  ),  &
     (filei(3),file3   ), (bufi(3) ,buf5    ),  &
     (nomat(3),nomat3  ), (bufd(1) ,db(13)  ), (buf(1)  ,db(2)   )
 EQUIVALENCE      (mcbs( 1),filea(1)), (mcbs( 8),typalp  ),  &
     (mcbs( 9),alpha(1)), (mcbs(13),fileb(1)),  &
     (mcbs(20),typbet  ), (mcbs(21),beta(1) ), (mcbs(61),filec(1))
 DATA    mcb   /  201   ,9*0   ,202   ,9*0   ,203   ,29*0  /,  &
     casecc,  mpool ,eqex  ,tfpool                     /  &
     101   ,  102   ,103   ,105                        /,  &
     scr1  ,  scr2  ,scr3  ,scr4  ,scr5  ,scr6  ,scr7  /  &
     301   ,  302   ,303   ,304   ,305   ,306   ,307   /,  &
     BLOCK /  81*0         /, nam   /  4HMTRX,4HIN  /,  &
     dmig  /  114   ,1     /, nfiles/  21           /,  &
     imat1 ,  imat2 ,imat3 ,itf /  139,   141,  143, 15/
 
!     PERFORM GENERAL INITIALIZATION
 
 buf1  = korsz(z) - sysbuf - 2
 buf2  = buf1 - sysbuf
 buf3  = buf2 - sysbuf
 buf4  = buf3 - sysbuf
 buf5  = buf4 - sysbuf
 nomat1= -1
 nomat2= -1
 nomat3= -1
 mask16= jhalf
 
!     IF MACHINE IS MORE THEN 32 BITS PER WORD, WE USE PACKING LOGIC
!     OTHERWISE, WE DO NOT PACK ROW AND COLUMN INDICES INTO ONE WORD
 
 pack  = .false.
 IF (ihalf > 16) pack = .true.
 DO  i = 1,nfiles,10
   j1 = i + 1
   jn = i + 9
   DO  j = j1,jn
     mcb(j)  = 0
   END DO
 END DO
 tfset   = 0
 typalp  = 1
 alpha(1)= 1.0
 typbet  = 1
 beta(1) = 1.0
 nogo    = 0
 
!     OPEN MPOOL. IF PURGED, SET FLAG.
 
 FILE   = mpool
 nompoo = 0
 nodmig = 0
 CALL preloc (*30,z(buf1),mpool)
 nompoo = 1
 CALL locate (*30,z(buf1),dmig,flag)
 nodmig = 1
 
!     READ CASE CONTROL RECORD.
!     SET NAMES OF REQUESTED MATRICES.
!     IF CASE CONTROL IS PURGED, SET NAMES OF REQUESTED MATRICES EQUAL
!     NAMES OF OUTPUT DATA BLOCKS.
 
 30 FILE = casecc
 CALL OPEN (*90,casecc,z(buf2),rdrew)
 CALL fwdrec (*680,casecc)
 CALL READ (*650,*50,casecc,z,buf2,1,flag)
 CALL mesage (-8,0,nam)
 50 CONTINUE
 CALL CLOSE (casecc,clsrew)
 tfset = z(itf)
 IF (z(imat1) == 0) GO TO 70
 nomat1 = 1
 mcb(8) = z(imat1  )
 mcb(9) = z(imat1+1)
 70 IF (z(imat2) == 0) GO TO 80
 nomat2  = 1
 mcb(18) = z(imat2  )
 mcb(19) = z(imat2+1)
 80 IF (z(imat3) == 0) GO TO 110
 nomat3  = 1
 mcb(28) = z(imat3  )
 mcb(29) = z(imat3+1)
 GO TO 110
 90 DO  i = 1,21,10
   mcb(31)  = mcb(i)
   CALL rdtrl (mcb(31))
   IF (mcb(31) == mcb(i)) CALL fname (mcb(31),mcb(i+7))
 END DO
 GO TO 290
 
!     IF TRANSFER FUNCTION MATRICES EXIST, BUILD THEM IN MATRIX FORMAT.
!     WRITE THEM ON 201,202,203 IF NO DMIG MATRICES TO ADD, OTHERWISE,
!     WRITE THEM ON SCR5,SCR6,SCR7.
!     IF NO TRANSFER FUNCTION MATRICES AND NO DMIG MATRICES, EXIT.
 
 110 IF (nomat1+nomat2+nomat3 == -3) GO TO 114
 IF (nodmig == 0) THEN
   GO TO   630
 ELSE
   GO TO   116
 END IF
 114 nodmig = 0
 116 IF (nodmig == 0 .AND. tfset == 0) GO TO 650
 IF (tfset == 0) GO TO 290
 file1 = scr5
 
!     TEST FOR PURGED OUTPUT DATA BLOCKS.
 
 DO  i = 1,21,10
   CALL rdtrl (mcb(i))
 END DO
 IF (nomat1 == -1) file1 = mcb( 1)
 file2 = scr7
 IF (nomat3 == -1) file2 = mcb(21)
 file3 = scr6
 IF (nomat2 == -1) file3 = mcb(11)
 nomat1 = 1
 nomat2 = 1
 nomat3 = 1
 
!     OPEN TFPOOL AND POSITION TO REQUESTED SET.
!     IF SET NOT IN TFPOOL, QUEUE MESSAGE AND TURN ON NOGO FLAG.
 
 FILE = tfpool
 CALL OPEN (*140,tfpool,z(buf2),rdrew)
 130 CALL fwdrec (*140,tfpool)
 CALL READ (*140,*140,tfpool,buf,1,0,flag)
 IF (buf(1) == tfset) GO TO 150
 GO TO 130
 140 nogo = 1
 buf(1) = tfset
 buf(2) = 0
 CALL mesage (30,74,buf)
 IF (dmig(1) == tfset) GO TO 150
 CALL CLOSE (tfpool,clsrew)
 GO TO 290
 
!     OPEN OUTPUT FILES. WRITE HEADER RECORDS.
 
 150 DO  i = 1,3
   
!     CHECK FOR PURGED OUTPUT DATA BLOCKS.
   
   IF (filei(i) > 0) GO TO 152
   nomat(i) = -1
   CYCLE
   152 FILE = filei(i)
   bufx = bufi(i)
   CALL gopen (FILE,z(bufx),wrtrew)
 END DO
 
!     PACK MATRICES ONTO OUTPUT FILES.
 
 ncol = luset
 icol = 1
 jsw  = 0
 isw  = 0
 i45  = 5
 IF (pack) i45 = 4
 i12  = i45 - 3
 180 DO  i = 1,3
   IF (filei(i) <= 0) CYCLE
   CALL bldpk (1,1,filei(i),BLOCK(20*i-19),1)
 END DO
 IF (isw /= 0) GO TO 210
 200 IF (jsw /= 0) GO TO 240
 CALL READ (*680,*260,tfpool,buf,i45,0,flag)
 isw = 1
 col = buf(1)
 row = buf(2)
 IF (.NOT.pack) GO TO 210
 col = rshift(buf(1),ihalf)
 row = andf(buf(1),mask16)
 210 IF (col > icol) GO TO 240
 DO  i = 1,3
   IF (filei(i) <= 0) CYCLE
   CALL bldpki (buf(i+i12),row,filei(i),BLOCK(20*i-19))
 END DO
 isw = 0
 GO TO 200
 240 DO  i = 1,3
   IF (filei(i) <= 0) CYCLE
   CALL bldpkn (filei(i),BLOCK(20*i-19),BLOCK(7*i+54))
 END DO
 icol = icol + 1
 IF (icol <= ncol) GO TO 180
 GO TO 270
 260 jsw = 1
 GO TO 240
 
!     CLOSE FILES AND WRITE TRAILERS. IF NO DMIG MATRICES, RETURN
 
 270 CALL CLOSE (tfpool,clsrew)
 DO  i = 1,3
   IF (filei(i) <= 0) CYCLE
   i7 = 7*i
   BLOCK(i7+54) = filei(i)
   BLOCK(i7+56) = ncol
   BLOCK(i7+57) = 1
   BLOCK(i7+58) = 1
   CALL CLOSE  (filei(i),1)
   CALL wrttrl (BLOCK(i7+54))
 END DO
 IF (nodmig == 0) GO TO 650
 
!     READ EQUIVALENCE TABLE INTO CORE
 
 290 FILE = eqex
 CALL gopen (eqex,z(buf2),0)
 CALL skprec (eqex,1)
 CALL READ (*680,*300,eqex,z,buf2,1,neqex)
 CALL mesage (-8,0,nam)
 300 CALL CLOSE (eqex,clsrew)
 kn = neqex/2
 nn = neqex - 1
 DO  i = 1,nn,2
   z(i+1) = z(i+1)/10
 END DO
 
!     READ MATRIX HEADER INFORMATION.
!     LOOK UP MATRIX NAME IN NAME LIST. IF ABSENT, SKIP MATRIX.
 
 320 CALL READ (*680,*630,mpool,buf,9,0,flag)
 
!     BUF(5)= INPUT MATRIX TYPE, BUF(6)= OUTOUT MATRIX TYPE
!     BUF(1) AND BUF(2) ARE MATRIX NAME FROM DMIG CARDS.
 
 k    = buf(6)
 prec = prc(k)
 k    = buf(5)
 iprc = MOD(k,2)
 nwd  = nwds(k)
 nwd1 = nwd + 1
 i11  = 11
 IF (pack) GO TO 325
 i11  = 10
 nwd1 = nwd + 2
 325 i10  = i11 - 1
 DO  i = 1,nfiles,10
   IF (mcb(i+7) == buf(1) .AND. mcb(i+8) == buf(2)) GO TO 360
 END DO
 340 CALL fread (mpool,buf,2,0)
 IF (buf(1) == -1) GO TO 320
 350 CALL fread (mpool,buf,2,0)
 IF (buf(1) == -1) GO TO 340
 CALL fread (mpool,buf,-nwd,0)
 GO TO 350
 
!     OPEN SCRATCH FILE. SET POINTERS.
 
 360 iptr = i
 FILE = mcb(iptr)
 mcb(iptr+2) = luset
 mcb(iptr+3) = buf(4)
 mcb(iptr+4) = buf(6)
 mcb(iptr+9) = 1
 iqq = (iptr-1)/10
 nomat(iqq+1) = +1
 isw = 0
 imtrx = neqex + 1
 i = imtrx
 
 
 CALL OPEN (*670,scr1,z(buf2),wrtrew)
 
!     READ COLUMN GRID AND COMPONENT, AND CHECK DUPLICATE.
!     CONVERT GRID AND COMPONENT TO SIL NO.
 
!     REMOVE DUPLICATE CHECK ADDED HERE IN 91 VERSION. CHECKING SHOULD
!     BE DONE EARLY IN IFP MODULE, AND NOT HERE. REMOVED ALSO NOGOX AND
!     ITS ASSOCIATED LINES.   3/93
 
 370 CALL fread (mpool,buf(10),2,0)
 IF (buf(10) == -1) GO TO 450
 
 
 ASSIGN 380 TO ret
 GO TO 710
 380 col = z(2*k)
 IF (buf(11) /= 0) col = col + buf(11) - 1
 IF (pack) col = lshift(col,ihalf)
 
!     READ A COLUMN OF THE MATRIX.
!     STORE IN CORE OR ON SCRATCH FILE IF TOO BIG FOR CORE.
 
 390 CALL fread (mpool,buf(10),2,0)
 IF (buf(10) == -1) GO TO 370
 ASSIGN 400 TO ret
 GO TO 710
 400 row = z(2*k)
 IF (buf(11) /= 0) row = row + buf(11) - 1
 buf(11) = row
 buf(10) = col
 IF (pack) buf(11) = row + col
 CALL fread (mpool,buf(12),nwd,0)
 IF (isw == 0) GO TO 420
 410 CALL WRITE (scr1,buf(i11),nwd1,0)
 GO TO 390
 420 IF (i+nwd1 < buf2) GO TO 430
 isw = 1
 CALL WRITE (scr1,z(imtrx),i-imtrx,0)
 GO TO 410
 430 DO  j = 1,nwd1
   z(i) = buf(j+i10)
   i = i + 1
 END DO
 GO TO 390
 
!     SORT MATRIX.
 
 450 IF (isw   == 0) GO TO 460
 CALL WRITE (scr1,0,0,1)
 CALL CLOSE (scr1,clsrew)
 CALL OPEN  (*670,scr1,z(buf2),rdrew)
 ifile(1) = scr2
 ifile(2) = scr3
 ifile(3) = scr4
 IF (     pack) CALL sorti  (scr1,0,nwd1,1,z(imtrx),buf2-imtrx)
 IF (.NOT.pack) CALL sorti2 (scr1,0,nwd1,1,z(imtrx),buf2-imtrx)
 CALL CLOSE (scr1,clsrew)
 GO TO 470
 460 n = i - imtrx
 nmtrx = i - nwd1
 CALL CLOSE (scr1,clsrew)
 IF (     pack) CALL sorti  (0,0,nwd1,1,z(imtrx),n)
 IF (.NOT.pack) CALL sorti2 (0,0,nwd1,1,z(imtrx),n)
 
!     OPEN OUTPUT FILE. WRITE HEADER RECORD
!     IF SORTED MATRIX NOT IN CORE, OPEN FILE WITH MATRIX.
 
 470 IF (tfset /= 0) FILE = scr1
 CALL OPEN  (*670,FILE,z(buf2),wrtrew)
 CALL fname (FILE,buf(19))
 CALL WRITE (FILE,buf(19),2,1)
 IF (isw /= 0) CALL OPEN (*670,ifile(6),z(buf3),rdrew)
 
!     PACK MATRIX ONTO OUTPUT FILE.
 
 ncol = luset
 j    = imtrx
 icol = 1
 jsw  = 0
 490 CALL bldpk (buf(6),buf(6),FILE,0,0)
 IF (jsw /=   0) GO TO 540
 500 IF (j > nmtrx) GO TO 570
 IF (isw ==   0) GO TO 510
 CALL READ (*680,*580,ifile(6),buf(i11),nwd1,0,flag)
 GO TO 530
 510 DO  k = 1,nwd1
   buf(k+i10) = z(j)
   j   = j + 1
 END DO
 530 col = buf(10)
 row = buf(11)
 IF (.NOT.pack) GO TO 540
 col = rshift(buf(11),ihalf)
 row = andf(buf(11),mask16)
 540 IF (col > icol) GO TO 590
 jsw = 0
 IF (prec == 2) GO TO 550
 IF (iprc == 0) GO TO 545
 x(1) = bufr(12)
 x(2) = bufr(13)
 GO TO 560
 545 x(1) = bufd(1)
 x(2) = bufd(2)
 GO TO 560
 550 IF (iprc == 0) GO TO 555
 xd(1) = bufr(12)
 xd(2) = bufr(13)
 GO TO 560
 555 xd(1) = bufd(1)
 xd(2) = bufd(2)
 560 CALL zblpki
 GO TO 500
 570 CALL bldpkn (FILE,0,mcb(iptr))
 icol = icol + 1
 IF (icol <= ncol) GO TO 490
 GO TO 600
 580 j = nmtrx + 1
 GO TO 570
 590 jsw = 1
 GO TO 570
 600 CALL CLOSE (FILE,clsrew)
 IF (isw /= 0) CALL CLOSE (ifile(6),clsrew)
 CALL wrttrl (mcb(iptr))
 
!     IF TRANSFER FUNCTION MATRICES ARE TO BE ADDED, CALL MATRIX ADD
!     ROUTINE THEN RETURN TO READ NEXT MATRIX IN THE MATRIX POOL.
 
 IF (tfset == 0) GO TO 320
 j = 2
 IF (iptr ==  1) j = 1
 IF (iptr == 11) j = 3
 DO  i = 1,7
   k = iptr + i - 1
   filea(i) = mcb(k)
   k = 7*j + i
   fileb(i) = BLOCK(k+53)
   filec(i) = 0
 END DO
 filea(1) = scr1
 filec(1) = mcb(iptr)
 filec(2) = ncol
 filec(3) = ncol
 filec(4) = filea(4)
 filec(5) = filea(5)
 nz = buf1 - imtrx
 nomats = 2
 k = orf(imtrx,1)
 CALL sadd (z(k),z(k))
 CALL wrttrl (filec)
 GO TO 320
 
!     TEST FOR ALL REQUESTED MATRICES FOUND.
 
 630 DO  i = 1,nfiles,10
   IF (mcb(i+7) == 0 .OR. mcb(i+9) /= 0) CYCLE
   CALL mesage (30,70,mcb(i+7))
   nogo = 1
 END DO
 650 IF (nompoo /= 0) CALL CLOSE (mpool,clsrew)
 IF (nogo   /= 0) CALL mesage (-61,0,nam)
 RETURN
 
!     FATAL ERRORS
 
 670 n = -1
 GO TO 700
 680 n = -2
 700 CALL mesage (n,FILE,nam)
 
!     BINARY SEARCH ROUTINE
 
 710 klo = 1
 khi = kn
 720 k = (klo+khi+1)/2
 730 IF (buf(10)-z(2*k-1) < 0.0) THEN
   GO TO   740
 ELSE IF (buf(10)-z(2*k-1) == 0.0) THEN
   GO TO   810
 ELSE
   GO TO   750
 END IF
 740 khi = k
 GO TO 760
 750 klo = k
 760 IF (khi-klo-1 < 0) THEN
   GO TO   800
 ELSE IF (khi-klo-1 == 0) THEN
   GO TO   770
 ELSE
   GO TO   720
 END IF
 770 IF (k == klo) GO TO 780
 k = klo
 GO TO 790
 780 k = khi
 790 klo = khi
 GO TO 730
 800 nogo = 1
 buf(11) = 0
 810 GO TO ret, (380,400)
END SUBROUTINE mtrxin
