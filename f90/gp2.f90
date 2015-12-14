SUBROUTINE gp2
     
!     GP2 BUILDS THE ELEMENT CONNECTION TABLE (ECT).
!     STRUCTURAL ELEMENT CONNECTION CARDS ARE ON GEOM2.
!     EACH EXTERNAL GRID PT. NO. IS CONVERTED TO AN INTERNAL INDEX.
!     IN ADDITION, GENERAL ELEMENT CARDS ARE READ AND
!     EXTERNAL GRID NUMBERS ARE CONVERTED TO INTERNAL NUMBERS.
 
 
 INTEGER :: elem  ,sysbuf,buf1  ,buf2  ,eqexin,rd    ,rdrew ,  &
     wrt   ,wrtrew,clsrew,cls   ,ect   ,geomp ,b     ,  &
     FILE  ,z     ,genel ,geom2 ,ret   ,ret1  ,gp2h  , cbar  ,cbeam ,buf3  ,two
 DIMENSION       b(34) ,gp2h(2)      ,mcb(7)       ,genel(2)
 COMMON /BLANK / noect
 COMMON /zzzzzz/ z(1)
 COMMON /gpta1 / nelem ,last  ,incr  ,elem(1)
 COMMON /system/ sysbuf,junk(36)     ,iaxif ,nbpc  ,nbpw
 COMMON /names / rd    ,rdrew ,wrt   ,wrtrew,clsrew,cls
 COMMON /setup / nfile(6)
 COMMON /two   / two(32)
 EQUIVALENCE     (geomp,geom2)
 
!     INPUT  DATA FILES
 DATA   geom2,eqexin / 101,102 /
 
!     OUTPUT DATA FILES
 DATA   ect / 201 /
 
!     MISC   DATA
 DATA   gp2h/ 4HGP2 ,4H    /, cbar / 4HBAR /, cbeam / 4HBEAM /
 
!     GENEL DATA CARDS PROCESSED BY GP2 IN ADDITION TO ELEMENTS.
 DATA  genel / 4301, 43 /
 
 
!     PERFORM GENERAL INITIALIZATION
 
 CALL delset
 buf1  = korsz(z) - sysbuf - 2
 buf2  = buf1 - sysbuf
 noect = -1
 buf3  = buf2 - sysbuf
 mcb(1)= geom2
 CALL rdtrl (mcb)
 
!     READ EQEXIN INTO CORE
 
 FILE = eqexin
 CALL OPEN (*580,eqexin,z(buf1),rdrew)
 CALL fwdrec (*590,eqexin)
 CALL READ (*590,*30,eqexin,z,buf2,1,n)
 CALL mesage (-8,0,gp2h)
 30 CALL CLOSE (eqexin,clsrew)
 kn = n/2
 n1 = n + 1
 
!     OPEN GEOM2. IF PURGED, RETURN.
!     OTHERWISE, OPEN ECT AND WRITE HEADER RECORD.
 
 nogeo2 = 0
 CALL preloc (*50,z(buf1),geom2)
 nogeo2 = 1
 GO TO 60
 50 RETURN
 
 60 noect = 1
 nogo  = 0
 FILE  = ect
 CALL OPEN (*580,ect,z(buf2),wrtrew)
 CALL fname (ect,b)
 CALL WRITE (ect,b,2,1)
 
!     READ 3-WORD ID FROM GEOM2. SEARCH ELEMENT TABLE FOR MATCH.
!     IF FOUND, BRANCH TO ELEMENT CODE. IF NOT FOUND, SEARCH GENEL
!     TABLE  FOR MATCH. IF FOUND BRANCH TO APPROPRIATE CODE. IF NOT
!     FOUND, SKIP RECORD AND CONTINUE.
 
 70 CALL READ (*460,*600,geom2,b,3,0,flag)
 DO  i = 1,last,incr
   IF (elem(i+3) == b(1)) GO TO 120
 END DO
 IF (genel(1) == b(1)) GO TO 110
 CALL fwdrec (*460,geom2)
 GO TO 70
 110 k = (i+1)/2
 GO TO 280
 
!     WRITE 3-WORD ID ON ECT. READ ALL CARDS FOR ELEMENT AND
!     CONVERT EXTERNAL GRID NOS. TO INTERNAL NOS.  WRITE ENTRIES ON ECT
!     DIRECTLY AFTER CONVERSION.
 
 120 ASSIGN 170 TO ret
 ASSIGN 630 TO ret1
 CALL WRITE (ect,b,3,0)
 m  = elem(i+5)
 lx = elem(i+12)
 mm = lx + elem(i+9)
 NAME = elem(i)
 ii   = n1
 FILE = geom2
 150 CALL READ (*590,*270,FILE,b,m,0,flag)
 
!     CHECK LATER TO SEE IF RESTRICTION APPLIES TO AXIF PROBLEMS
 
 IF (iaxif /= 0) GO TO 155
 IF (nbpw <= 32 .AND. b(1) > 16777215) GO TO 670
!                                  16777215 = 2**24 - 1
 155 l = lx
 160 IF (b(l) /= 0) GO TO 470
 170 l= l + 1
 IF (l    <    mm) GO TO 160
 IF (NAME == cbeam) GO TO 180
 IF (NAME /=  cbar) GO TO 200
 
!     SPECIAL PROCESSING FOR BAR AND BEAM ELEMENTS
 
 IF (b(8) == 1) GO TO 200
 ASSIGN 190 TO ret
 l = 5
 GO TO 470
 180 IF (b(8) == 0) GO TO 200
 ASSIGN 190 TO ret
 l = 8
 GO TO 470
 190 ASSIGN 170 TO ret
 
 200 CALL WRITE (ect,b,m,0)
 GO TO 150
 
!     CURRENT ELEMENT IS COMPLETE
 
 270 CALL WRITE (ect,0,0,1)
 GO TO 70
 
!     GENERAL ELEMENTS-- WRITE 3-WORD ID ON ECT. READ ALL GENELS,
!     CONVERT EXTERNAL GRID NOS. TO INTERNAL NOS. AND WRITE THEM ON ECT.
 
 280 CALL WRITE (ect,b,3,0)
 FILE = geom2
 l = 2
 ASSIGN 310 TO ret
 ASSIGN 640 TO ret1
 290 ijk = 0
 CALL READ (*590,*360,geom2,b,1,0,flag)
 CALL WRITE (ect,b,1,0)
 300 CALL READ (*590,*600,geom2,b(2),2,0,flag)
 IF (b(2) == -1) GO TO 320
 GO TO 470
 310 CALL WRITE (ect,b(2),2,0)
 GO TO 300
 320 nud = b(3)
 IF (ijk /= 0) GO TO 330
 nui = b(3)
 ijk = 1
 GO TO 310
 330 CALL WRITE (ect,b(2),2,0)
 CALL READ (*590,*600,geom2,ijk1,1,0,flag)
 CALL WRITE (ect,ijk1,1,0)
 ncore = buf2 - n1
 nz = (nui*(nui+1))/2
 nread = 0
 340 n= MIN0(ncore,nz-nread)
 CALL READ (*590,*600,geom2,z(n1),n,0,flag)
 CALL WRITE (ect,z(n1),n,0)
 nread = nread + n
 IF (nread < nz) GO TO 340
 CALL READ (*590,*600,geom2,ijk,1,0,flag)
 CALL WRITE (ect,ijk,1,0)
 IF (ijk == 0) GO TO 290
 ns = nui*nud
 nread = 0
 350 n= MIN0(ncore,ns-nread)
 CALL READ (*590,*600,geom2,z(n1),n,0,flag)
 CALL WRITE (ect,z(n1),n,0)
 nread = nread + n
 IF (nread < ns) GO TO 350
 GO TO 290
 360 CALL WRITE (ect,0,0,1)
 GO TO 70
 
!     CLOSE FILES, WRITE TRAILER AND RETURN.
 
 460 CALL CLOSE (geom2,clsrew)
 CALL CLOSE (ect  ,clsrew)
 mcb(1) = geom2
 CALL rdtrl (mcb)
 mcb(1) = ect
 CALL wrttrl (mcb)
 IF (nogo /= 0) CALL mesage (-61,0,0)
 RETURN
 
 
!     INTERNAL BINARY SEARCH ROUTINE
!     ==============================
 
 470 klo = 1
 khi = kn
 igrid = b(l)
 480 k = (klo+khi+1)/2
 490 IF (igrid-z(2*k-1) < 0) THEN
   GO TO   500
 ELSE IF (igrid-z(2*k-1) == 0) THEN
   GO TO   560
 ELSE
   GO TO   510
 END IF
 500 khi = k
 GO TO 520
 510 klo = k
 520 IF (khi-klo-1 < 0) THEN
   GO TO   570
 ELSE IF (khi-klo-1 == 0) THEN
   GO TO   530
 ELSE
   GO TO   480
 END IF
 530 IF (k == klo) GO TO 540
 k = klo
 GO TO 550
 540 k = khi
 550 klo = khi
 GO TO 490
 560 b(l) = z(2*k)
 GO TO ret,  (170,310,190)
 570 GO TO ret1, (630,640)
 
 
!     FATAL ERROR MESSAGES
 
 580 j = -1
 GO TO 610
 590 j = -2
 GO TO 610
 600 j = -3
 610 CALL mesage (j,FILE,gp2h)
 630 k = 7
 GO TO 660
 640 k = 61
 660 b(2) = igrid
 CALL mesage (30,k,b)
 nogo = 1
 GO TO ret, (170,310)
 670 nogo = 1
 CALL mesage (30,138,b)
 GO TO 155
END SUBROUTINE gp2
