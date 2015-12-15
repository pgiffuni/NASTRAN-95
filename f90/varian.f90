SUBROUTINE varian
     
!     VARIANCE ANALYSIS POST PROCESSOR MODULE
 
!     INPUTS--O1,O2,O3,O4,O5 OR EDT
 
!     OUTPUTS--O1O,O2O,O3O,O4O,O5O
 
!     PARAMETERS--OP--BCD'DER' OR 'VAR'
!                 DELTA--REAL--DEFAULT=1.0
 
 LOGICAL :: tapbit
 INTEGER :: it1,it2,it3,it4,it5,ot1,ot2,ot3,ot4,ot5,NAME(2),mcb(7),  &
     sysbuf,scr1,scr2,der,var,FILE,itf(5),ito(5),varl(3), derl(3),varan(2),iz(260)
 REAL :: z(200)
 COMMON /system/ sysbuf,skip(91),jrun
 COMMON /BLANK / iop(2),delta
 COMMON /zzzzzz/ rz(1)
 EQUIVALENCE (it1,itf(1)),(itf(2),it2) ,(itf(3),it3),(itf(4),it4),  &
     (itf(5),it5),(ito(1) ,ot1),(ito(2),ot2),(ito(3),ot3),  &
     (ito(4),ot4),(ito(5),ot5) ,(z(1),iz(1),rz(1))
 DATA it1,it2,it3,it4,it5,ot1,ot2,ot3,ot4,scr1,ot5,inpt, scr2 /  &
     101,102,103,104,105,201,202,203,204,301 ,205,4HINPT,302 /
 DATA  der   , var   , derl  ,                 NAME             /  &
     4HDER , 4HVAR , 4HDERI, 4HVATI, 4HVE  , 4HVARI, 4HAN     /
 DATA  varl  ,                 iblnk , varan   ,       mcb      /  &
     4HVARI, 4HANCE, 4H    , 4H    , 4202, 42,       7*0      /
 
 ibuf1 = korsz(z(1))-sysbuf
 ibuf2 = ibuf1-sysbuf
 ibuf3 = ibuf2-sysbuf
 nz    = ibuf3-1
 IF (nz <= 0) CALL mesage(-8,0,NAME)
 IF (.NOT.tapbit(inpt)) CALL mesage(-7,0,NAME)
 CALL int2a8 (*5,jrun,iz)
 5 njrun=iz(1)
 nw=1
 IF (iop(1) == var) GO TO 300
 IF (iop(1) /= der) RETURN
 
!     DERIVATIVES SECTION
 
 IF (jrun /= 0) GO TO 30
 
!     COPY INPUT FILES TO INPT TAPE
 
 CALL OPEN (*900,inpt,iz(ibuf1),1)
 i=1
 ASSIGN 20 TO iret
 10 FILE = itf(i)
 GO TO 700
 20 i = i + 1
 IF (i <= 5) GO TO 10
 CALL CLOSE (inpt,2)
 RETURN
 
!     COMPUTE DERIVATIVES  DJ = (OJ - O0)/DELTA
 
 30 CONTINUE
 CALL OPEN (*900,inpt,iz(ibuf1),0)
 
 DO  i = 1,5
   ifound = 0
   CALL fwdrec (*220,inpt)
   CALL OPEN (*230,itf(i),iz(ibuf2),0)
   CALL fwdrec (*910,itf(i))
   CALL gopen (ito(i),iz(ibuf3),1)
   FILE = itf(i)
   40 ASSIGN 60 TO irtn
   50 CALL READ (*220,*920,inpt,iz(1),146,1,iflag)
   GO TO irtn, (60,70)
   60 CALL READ (*230,*920,FILE,iz(147),146,1,iflag)
   
!     CHECK FOR MATCH ON SUBCASE
   
   70 IF (iz(4)-iz(150) < 0) THEN
     GO TO    80
   ELSE IF (iz(4)-iz(150) == 0) THEN
     GO TO   100
   ELSE
     GO TO    90
   END IF
   
!     NEED NEW INPT RECORD
   
   80 CALL fwdrec (*910,inpt)
   ASSIGN 70 TO irtn
   GO TO 50
   
!     NEED NEW FILE RECORD
   
   90 CALL fwdrec (*910,FILE)
   GO TO 60
   
!     CHECK FOR MATCH ON TIME, FREQ ETC
   
   100 IF (z(5)-z(151) < 0.0) THEN
     GO TO    80
   ELSE IF (z(5)-z(151) == 0.0) THEN
     GO TO   110
   ELSE
     GO TO    90
   END IF
   
!     CHECK FOR MATCH ON ELTYPE
   
   110 IF (iz(3)-iz(149) < 0) THEN
     GO TO    80
   ELSE IF (iz(3)-iz(149) == 0) THEN
     GO TO   120
   ELSE
     GO TO    90
   END IF
   
!     WE GOT ONE
   
   120 CONTINUE
   iz(257) = derl(1)
   iz(258) = derl(2)
   iz(259) = derl(3)
   iz(260) = njrun
   ifound  = 1 + ifound
   CALL WRITE (ito(i),iz(147),146,1)
   nrec=iz(10)
   130 ASSIGN 150 TO irtn1
   140 CALL READ (*910,*190,inpt,iz(1),nrec,0,iflag)
   id1 = iz(1)/10
   GO TO irtn1, (150,160)
   150 CALL READ (*910,*200,FILE,iz(nrec+1),nrec,0,iflag)
   id2 = iz(nrec+1)/10
   ASSIGN 160 TO irtn1
   160 IF (id1-id2 < 0) THEN
     GO TO   140
   ELSE IF (id1-id2 == 0) THEN
     GO TO   170
   ELSE
     GO TO   150
   END IF
   
!     POINT CHECKS
   
   170 CONTINUE
   DO  j=2,nrec
     itype = numtyp(iz(j))
     IF (itype /= 2 .AND. itype /= 0) CYCLE
     z(nrec+j) = (z(nrec+j) - z(j))/delta
   END DO
   CALL WRITE (ito(i),iz(nrec+1),nrec,0)
   GO TO 130
   
!     END OF DATA RECORD
   
   190 CALL fwdrec (*910,FILE)
   GO TO 210
   200 CALL fwdrec (*910,inpt)
   210 CALL WRITE (ito(i),0,0,1)
   GO TO 40
   
!     EOF ON INPT
   
   220 GO TO 240
   
!     EOF ON FILE
   
   230 CALL skpfil (inpt,1)
   240 CALL CLOSE (FILE,1)
   CALL CLOSE (ito(i),1)
   mcb(1)=ito(i)
   mcb(2)=ifound
   IF (ifound /= 0) CALL wrttrl (mcb)
 END DO
 
!     SKIP OVER OLD DERIVATIVES
 
 i = 5*jrun - 5
 CALL skpfil (inpt,i)
 CALL CLOSE (inpt,2)
 CALL gopen (inpt,iz(ibuf1),3)
 i = 1
 ASSIGN 270 TO iret
 260 FILE = ito(i)
 mcb(1) = FILE
 CALL rdtrl (mcb)
 IF (mcb(2) /= 0) GO TO 700
 CALL eof (inpt)
 270 i = i + 1
 IF (i <= 5) GO TO 260
 CALL CLOSE (inpt,2)
 280 RETURN
 
!     VARIANCE SECTION
 
 300 IF (jrun == 0) RETURN
 
!     SEE IF VARIANCE IS TO BE COMPUTED
 
 CALL preloc (*280,iz(ibuf1),it1)
 CALL locate (*320,iz(ibuf1),varan,iflag)
 
!     READ IN VARIANCES
 
 CALL READ (*910,*310,it1,iz(1),nz,0,iflag)
 CALL mesage (-8,0,NAME)
 310 IF (iflag-1 == jrun) GO TO 330
 320 CALL CLOSE (it1,1)
 RETURN
 
!     SET UP FOR VARIANCES
 
 330 CALL CLOSE (it1,1)
 CALL OPEN (*900,inpt,iz(ibuf1),0)
 CALL skpfil (inpt,5)
 in1 = inpt
 io1 = scr1
 DO  i = 1,jrun
   IF (i == jrun) GO TO 340
   CALL OPEN (*900,io1,iz(ibuf2),1)
   340 CONTINUE
   IF (i == 1) GO TO 350
   CALL OPEN (*900,in1,iz(ibuf3),0)
   350 CONTINUE
   DO  j = 1,5
     IF (i /= jrun) GO TO 360
     
!     FIX UP FOR WRITING ON OUTPUT FILES
     
     CALL OPEN (*610,ito(j),iz(ibuf2),1)
     CALL fname (ito(j),mcb)
     CALL WRITE (ito(j),mcb,2,1)
     ifound = 0
     io1 = ito(j)
     360 CONTINUE
     CALL fwdrec (*590,inpt)
     370 ASSIGN 400 TO irtn2
     380 CALL READ (*600,*920,in1,iz(jrun+1),146,1,iflag)
     IF (i /= 1) GO TO 390
     iz(jrun+111) = varl(1)
     iz(jrun+112) = varl(2)
     iz(jrun+113) = varl(3)
     iz(jrun+114) = iblnk
     GO TO 460
     
!     CHECK FOR MATCH
     
     390 GO TO irtn2, (400,410)
     400 CALL READ (*580,*920,inpt,iz(jrun+147),146,1,iflag)
     410 IF (iz(jrun+4)-iz(jrun+150) < 0) THEN
       GO TO   420
     ELSE IF (iz(jrun+4)-iz(jrun+150) == 0) THEN
       GO TO   440
     ELSE
       GO TO   430
     END IF
     420 CALL fwdrec (*910,in1)
     ASSIGN 410 TO irtn2
     GO TO 380
     430 CALL fwdrec (*910,inpt)
     GO TO 400
     440 IF ( z(jrun+5)- z(jrun+151) < 0.0) THEN
       GO TO   420
     ELSE IF ( z(jrun+5)- z(jrun+151) == 0.0) THEN
       GO TO   450
     ELSE
       GO TO   430
     END IF
     450 IF (iz(jrun+3)-iz(jrun+149) < 0) THEN
       GO TO   420
     ELSE IF (iz(jrun+3)-iz(jrun+149) == 0) THEN
       GO TO   460
     ELSE
       GO TO   430
     END IF
     
!     MATCH
     
     460 CALL WRITE (io1,iz(jrun+1),146,1)
     nrec = iz(jrun+10)
     m = jrun + nrec
     470 ASSIGN 490 TO irtn3
     480 CALL READ (*910,*550,in1,iz(jrun+1),nrec,0,iflag)
     IF (i == 1) GO TO 510
     id1 = iz(jrun+1)/10
     GO TO irtn3, (490,500)
     490 CALL READ (*910,*560,inpt,iz(m+1),nrec,0,iflag)
     id2 = iz(m+1) /10
     ASSIGN 500 TO irtn3
     500 IF (id1-id2 < 0) THEN
       GO TO   480
     ELSE IF (id1-id2 == 0) THEN
       GO TO   510
     ELSE
       GO TO   490
     END IF
     
!     POINT MATCH
     
     510 CONTINUE
     IF (i == jrun) ifound = ifound +1
     DO  k =2,nrec
       itype = numtyp(iz(jrun+k))
       IF (itype /= 2 .AND. itype /= 0) CYCLE
       IF (i /= 1) GO TO 520
       z(jrun+k) = (z(jrun+k)*z(1))**2
       GO TO 530
       520 z(jrun+k) = z(jrun+k) + (z(m+k)*z(i))**2
       530 IF (i /= jrun) CYCLE
       z(jrun+k) = SQRT(z(jrun+k))
     END DO
     CALL WRITE (io1,iz(jrun+1),nrec,0)
     GO TO 470
     
!     END OF DATA ON IN1
     
     550 IF (i == 1) GO TO 570
     CALL fwdrec (*910,inpt)
     GO TO 570
     560 CALL fwdrec (*910,in1)
     570 CALL WRITE (io1,0,0,1)
     GO TO 370
     
!     EOF ON INPT
     
     580 CALL skpfil (in1,1)
     590 CALL eof (io1)
     IF (i /= jrun) CYCLE
     CALL CLOSE (io1,1)
     mcb(1) = io1
     mcb(2) = ifound
     CALL wrttrl (mcb)
     CYCLE
     600 IF (i == 1) GO TO 590
     CALL skpfil (inpt,1)
     GO TO 590
     610 CONTINUE
   END DO
   
!     SWITCH FILES
   
   IF (i /= jrun) CALL CLOSE (io1,1)
   IF (i /=    1) CALL CLOSE (in1,1)
   j=in1
   in1=io1
   io1=j
   IF (i == 1) io1 = scr2
 END DO
 CALL CLOSE (inpt,1)
 jrun = 9999999
 RETURN
 
!     INTERNAL ROUTINE TO COPY FILES
 
 700 CONTINUE
 CALL OPEN (*730,FILE,iz(ibuf2),0)
 710 IEOR = 1
 CALL READ (*730,*720,FILE,iz(1),nz,0,iread)
 IEOR = 0
 720 CALL WRITE (inpt,iz(1),iread,IEOR)
 GO TO 710
 730 CALL eof (inpt)
 CALL CLOSE (FILE,1)
 GO TO iret, (20,270)
 
!     ERROR MESSAGES
 
 900 ip1 = -1
 GO TO 930
 910 ip1 = -2
 GO TO 930
 920 ip1 = -3
 930 CALL mesage (ip1,FILE,NAME)
 STOP
END SUBROUTINE varian
