SUBROUTINE xcei
     
 IMPLICIT INTEGER (a-z)
 EXTERNAL        lshift,rshift,andf,orf
 DIMENSION       dcparm(2),nxcei(2),idic(1),nxptdc(2),contrl(4)
 COMMON /xvps  / vps(1)
!WKBR COMMON /XCEITB/ CEITBL(2)
 COMMON /xceitb/ ceitbl(42)
 COMMON /oscent/ buf(7)
 COMMON /zzzzzz/ databf(1)
 COMMON /system/ bfsz,isysou,dum21(21),icfiat,dum57(57),icpflg
 COMMON /xfiat / ifiat(3)
 COMMON /xdpl  / idpl(3)
 EQUIVALENCE     (databf(1),idic(1))
 DATA    nxptdc/ 4HXPTD,4HIC  /
 DATA    nxcei / 4HXCEI,4H    /
 DATA    noscar/ 4HXOSC/
 DATA    pool  / 4HPOOL/
 DATA    contrl/ 4HJUMP,4HREPT,4HCOND,4HEXIT/
 DATA    nblank/ 4H    /
 DATA    mask1 / 65535 /, noflgs / 536870911/
 
!     MASK1  = 000000177777 =     65536 = 2**16-1
!     NOFLGS = 003777777777 = 536870911 = 2**29-1
!     MASK   = 017777600000
!     LPFLG  = 010000000000
 
 mask  = lshift(mask1,16)
 lpflg = lshift(1,30)
 CALL OPEN (*310,pool,databf,2)
 
!     DETERMINE WHICH TYPE OF CONTROL REQUEST
 
 DO  j = 1,4
   IF (buf(4) == contrl(j)) THEN
      SELECT CASE ( j )
       CASE (    1)
         GO TO 150
       CASE (    2)
         GO TO 110
       CASE (    3)
         GO TO 250
       CASE (    4)
         GO TO 270
     END SELECT
   END IF
 END DO
 CALL mesage (-61,0,0)
 
!     PROCESS  JUMP CONTROL REQUEST
 
 30 IF (newsq > buf(2)) GO TO 60
 
!     MUST BACKSPACE WITHIN OSCAR FILE
!     DUE TO GINO TECHNIQUES IT IS USUALLY FASTER TO REWIND AND FORWARD
!     REC RATHER THAN BACKREC
 
 CALL REWIND (pool)
 
!     POSITION POOL TAPE AT BEGINNING OF OSCAR FILE
 
 jj = idpl(3)*3 + 1
 DO  j = 4,jj,3
   IF (idpl(j) == noscar) GO TO 50
 END DO
 CALL mesage (-61,0,0)
 50 CALL skpfil (pool,andf(idpl(j+2),mask1)-1)
 newsq = newsq - 1
 GO TO 70
 
!     MUST FORWARD REC WITHIN OSCAR FILE
 
 60 newsq = newsq - buf(2) - 1
 IF (newsq == 0) GO TO 260
 70 DO  i = 1,newsq
   CALL fwdrec (*290,pool)
 END DO
 
!     CHECK FOR REPEAT INSTRUCTION
 
 IF (buf(4) == contrl(2)) GO TO 260
 
!     JUMP REQUEST - CHECK FOR JUMP OUT OF LOOPS
 
 newsq = rshift(andf(buf(7),mask),16)
 kk = 3
 ceitbx = 0
 100 ceitbx = 4 + ceitbx
 IF (ceitbx > ceitbl(2)) GO TO 260
 IF (andf(ceitbl(ceitbx-1),lpflg) == 0 .OR. ceitbl(ceitbx+1) == 0) GO TO 100
 nbegn = rshift(andf(ceitbl(ceitbx-1),noflgs),16)
 nend  = andf(mask1,ceitbl(ceitbx-1))
 IF (newsq < nbegn .OR. newsq > nend) GO TO 130
 GO TO 100
 
!     PROCESS  REPEAT CONTROL REQUEST
 
 110 kk = 1
 120 ceitbx = andf(buf(7),mask1)
 IF (ceitbl(ceitbx) < 0.0) THEN
   GO TO   121
 ELSE
   GO TO   122
 END IF
 
!      NEGATIVE ENTRY IMPLIES VARIABLE REPT INSTRUCTION
!      FIND VALUE IN VPS AND UPDATE CEITBL
 
 121 ivpspt = rshift(andf(ceitbl(ceitbx),noflgs),16)
 loop   = andf(ceitbl(ceitbx),mask1)
 ivpspt = vps(ivpspt+3)
 ceitbl(ceitbx) = orf(lshift(ivpspt,16),loop)
 122 CONTINUE
 
!     CHECK FOR END OF LOOP
 
 mxloop = rshift(andf(ceitbl(ceitbx),noflgs),16)
 loop   = andf(ceitbl(ceitbx),mask1)
 IF (mxloop > loop) GO TO 140
 
!     REPEATS FINISHED - ZERO LOOP COUNT AND TURN OFF LOOP FLAG
 
 130 ceitbl(ceitbx  ) = andf(ceitbl(ceitbx  ),mask  )
 ceitbl(ceitbx-1) = andf(ceitbl(ceitbx-1),noflgs)
 SELECT CASE ( kk )
   CASE (    1)
     GO TO 260
   CASE (    2)
     GO TO 280
   CASE (    3)
     GO TO 100
 END SELECT
 
!     ANOTHER  TIME THRU - INCREMENT COUNTER BY 1
 
 140 ceitbl(ceitbx) = ceitbl(ceitbx) + 1
 
!     SET LOOP FLAG IN WORD 1 OF CEITBL ENTRY
 
 ceitbl(ceitbx-1) = orf(ceitbl(ceitbx-1),lpflg)
 SELECT CASE ( kk )
   CASE (    1)
     GO TO 150
   CASE (    2)
     GO TO 260
 END SELECT
 150 newsq = rshift(andf(buf(7),mask),16)
 
!     MAKE SURE WE ARE LOOPING
 
 IF (newsq >= buf(2)) GO TO 30
 
!     IF CHECKPOINTING - BACKUP PROBLEM TAPE DICTIONARY TO BEGINNING OF
!     LOOP
 
 IF (icpflg == 0) GO TO 210
 
!     READ IN CHECKPOINT DICTIONARY
 
 itop = 2*bfsz + 1
 ldic = korsz(idic(itop))
 CALL OPEN (*310,nxptdc,databf(bfsz+1),0)
 CALL READ (*300,*160,nxptdc,dcparm,2,1,nrecsz)
 160 IF (nxptdc(1) /= dcparm(1)) CALL mesage (-61,0,0)
 CALL READ (*300,*170,nxptdc,dcparm,2,1,nrecsz)
 170 CALL READ (*300,*180,nxptdc,idic(itop),ldic,1,nrecsz)
 GO TO 310
 180 ibot = nrecsz + itop - 3
 CALL CLOSE (nxptdc,1)
 j = ibot
 DO  i = itop,ibot,3
   IF (idic(i) /= nblank) CYCLE
   IF (andf(idic(i+2),mask1) < newsq) CYCLE
   j = i - 3
   EXIT
 END DO
 200 ibot = j
 
!     WRITE IDIC ON NEW PROBLEM TAPE
 
 CALL OPEN  (*310,nxptdc,databf(bfsz+1),1)
 CALL WRITE (nxptdc,nxptdc,2,1)
 CALL WRITE (nxptdc,dcparm,2,1)
 CALL WRITE (nxptdc,idic(itop),ibot+3-itop,1)
 CALL CLOSE (nxptdc,1)
 
!     SCAN FIAT FOR FILES REGENERATED NEXT TIME THRU LOOP.
 
 210 j  = ifiat(3)*icfiat - 2
 jj = idpl(3) *3 + 1
 DO  i = 4,j,icfiat
   IF (rshift(andf(ifiat(i),noflgs),16) >= buf(2)) CYCLE
   IF (rshift(andf(ifiat(i),noflgs),16) ==      0) CYCLE
   IF (andf(rshift(ifiat(i),30),1)      /=      0) CYCLE
   
!     LTU IS LESS THAN LOOP END - CLEAR FIAT TRAILER
   
   ifiat(i+ 3) = 0
   ifiat(i+ 4) = 0
   ifiat(i+ 5) = 0
   IF (icfiat == 8) GO TO 212
   ifiat(i+ 8) = 0
   ifiat(i+ 9) = 0
   ifiat(i+10) = 0
   
!     IF EQUIV, REMOVE ENTIRE ENTRY FROM FIAT
!     REMOVE ENTIRE ENTRY FROM FIAT TO FORCE REALLOCATION
   
   212 ihold = andf(mask1,ifiat(i))
   ifiat(i  ) = 0
   ifiat(i+1) = 0
   ifiat(i+2) = 0
   IF (i < ifiat(1)*icfiat) ifiat(i) = ihold
   
!     ZERO FILE NAME IF IN DPL
   
   DO  ii = 4,jj,3
     IF (idpl(ii) == ifiat(i+1) .AND. idpl(ii+1) == ifiat(i+2)) GO TO 230
   END DO
   CYCLE
   230 idpl(ii  ) = 0
   idpl(ii+1) = 0
 END DO
 GO TO 30
 
!     PROCESS  CONDITIONAL CONTROL REQUEST
 
 250 ceitbx = andf(buf(7),mask1)
 IF (vps(ceitbx) < 0) GO TO 150
 260 CALL CLOSE (pool,2)
 RETURN
 
!     PROCESS EXIT  CONTROL REQUESTS
 
 270 kk = 2
 IF (buf(7) /= contrl(4)) GO TO 120
 280 CALL pexit
 290 CALL mesage (-2,pool  ,nxcei)
 300 CALL mesage (-2,nxptdc,nxcei)
 310 CALL mesage (-61,0,0)
 RETURN
END SUBROUTINE xcei
