SUBROUTINE gp3b
     
!     GP3B BUILDS THE GRID POINT TEMPERATURE TABLE (GPTT).
!     TEMPD AND TEMP CARDS ARE READ.
!     THE GPTT HEADER CONTAINS THE FILE NAME PLUS 3 WORDS FOR EACH
!     TEMPERATURE SET.
!       WORD 1 = TEMP SET ID.
!       WORD 2 = DEFAULT TEMP OR -1 IF NO DEFAULT TEMP.
!       WORD 3 = RECORD NO. (AFTER HEADER RECORD) OF TEMPERATURE DATA
!                FOR THE SET, OR
!                ZERO IF ONLY A DEFAULT TEMP IS DEFINED FOR THE SET.
!     DATA RECORDS OF THE GPTT CONSIST OF PAIRS OF EXTERNAL INDEX AND
!     TEMPERATURE. EACH DATA RECORD IS SORTED ON EXTERNAL INDEX.
 
!     AN IDENTICAL SET OF RECORDS WITH INTERNAL INDICES IS APPENDED AT
!     THE END OF THE GPTT.
 
 
 LOGICAL :: intern
 INTEGER :: geomp ,eqexin,slt   ,gptt  ,scr1  ,buf1  ,buf2  ,  &
     buf   ,temp  ,tempd ,FILE  ,flag  ,z     ,rd    ,  &
     rdrew ,wrtrew,wrt   ,clsrew,nam(2),geom3 ,ett   ,  &
     tempp1,tempp2,tempp3,temprb,buf3  ,tempg ,tempp4
 CHARACTER (LEN=23) :: ufm
 COMMON /xmssg / ufm
 COMMON /BLANK / nograv,noload,notemp
 COMMON /gp3com/ geom3 ,eqexin,geom2 ,slt   ,ett   ,scr1  ,scr2  ,  &
     buf1  ,buf2  ,buf(50)      ,cardid(60)   ,idno(30)  &
     ,               carddt(60)   ,mask(60)     ,STATUS(60)   ,ntypes,  &
     ipload,igrav ,pload2(2)    ,load(2)      ,nopld2,  &
     temp(2)      ,tempd(2)     ,tempp1(2)           ,  &
     tempp2(2)    ,tempp3(2)    ,temprb(2)    ,buf3  ,  &
     pload3(2)    ,ipld3        ,tempg(2)            , tempp4(2)
 COMMON /names / rd    ,rdrew ,wrt   ,wrtrew,clsrew
 COMMON /zzzzzz/ z(1)
 COMMON /system/ sysbuf,nout
 EQUIVALENCE     (geom3,geomp),(gptt,scr1)
 DATA    nam   / 4HGP3B,4H    /
 
!     TURN NODEF FLAG ON
 
 id    = 0
 nodef = 0
 
!     READ EQEXIN INTO CORE
 
 FILE = eqexin
 CALL OPEN (*400,eqexin,z(buf2),rdrew)
 CALL fwdrec (*410,eqexin)
 CALL READ (*410,*10,eqexin,z,buf3,1,neqx)
 CALL mesage (-8,0,nam)
 10 CALL CLOSE (eqexin,clsrew)
 kn = neqx/2
 itempd = neqx + 1
 itabl  = itempd
 
!     READ TEMPERATURE DEFAULT CARDS (IF PRESENT)
 
 FILE = geomp
 CALL preloc (*460,z(buf1),geomp)
 CALL locate (*40,z(buf1),tempd,flag)
 i = itempd
 nodef  = 1
 notemp = 1
 20 CALL READ (*410,*30,geomp,z(i),2,0,flag)
 i = i + 2
 GO TO 20
 30 itabl  = i
 ntempd = i - 2
 n = itabl - itempd
 CALL sort (0,0,2,1,z(itempd),n)
 
!     READ TEMP CARDS.  DETERMINE NO. OF TEMP SETS
!     FOR EACH SET ID, LOOK UP THE DEFAULT TEMPERATURE
!     WRITE SET ID, DEFAULT TEMP (OR -1) AND RECORD NUMBER
!     OF THE TEMPERATURE DATA (OR 0) IN THE GPTT HEADER
 
 40 j = 0
 k = itempd
 i = itabl
 l = 1
 FILE = geomp
 CALL locate (*270,z(buf1),temp,flag)
 notemp = 1
 FILE   = gptt
 CALL OPEN  (*400,gptt,z(buf2),wrtrew)
 CALL fname (gptt,buf)
 CALL WRITE (gptt,buf,2,0)
 
!     OPEN ETT AS TEMPORARY SCRATCH TO FORM IDENTICAL FILE WITH
!     INTERNAL NOTATION
 
 FILE = ett
 CALL OPEN  (*400,ett,z(buf3),wrtrew)
 CALL fname (ett,buf)
 CALL WRITE (ett,buf,2,0)
 FILE = geomp
 50 CALL READ (*410,*110,geomp,buf,3,0,flag)
 j = j + 1
 IF (id == buf(1)) GO TO 50
 id = buf(1)
 z(i) = j
 i = i + 1
 IF (nodef ==  0) GO TO 80
 60 IF (k > ntempd) GO TO 80
 IF (id-z(k) < 0) THEN
   GO TO    80
 ELSE IF (id-z(k) == 0) THEN
   GO TO    90
 END IF
 70 buf(1) = z(k  )
 buf(2) = z(k+1)
 buf(3) = 0
 CALL WRITE (gptt,buf,3,0)
 CALL WRITE (ett ,buf,3,0)
 k = k + 2
 GO TO 60
 80 buf(2) = -1
 GO TO 100
 90 buf(2) = z(k+1)
 k = k + 2
 100 buf(3) = l
 buf(1) = id
 l = l + 1
 CALL WRITE (gptt,buf,3,0)
 CALL WRITE (ett ,buf,3,0)
 j = 0
 GO TO 50
 110 IF (nodef ==  0) GO TO 130
 IF (k > ntempd) GO TO 130
 buf(3) = 0
 DO  l = k,ntempd,2
   buf(1) = z(l  )
   buf(2) = z(l+1)
   CALL WRITE (ett ,buf,3,0)
   CALL WRITE (gptt,buf,3,0)
 END DO
 130 CALL WRITE (gptt,0,0,1)
 CALL WRITE (ett ,0,0,1)
 CALL bckrec (geomp)
 n = i
 z(n) = j + 1
 i = itabl + 1
 
!     READ EACH TEMP SET
!     SORT ON EXTERNAL INDEX AND WRITE ON GPTT
 
 ifile  = gptt
 intern = .false.
 isave  = i
 nogo   = 0
 140 CALL READ (*410,*420,geomp,0,-3,0,flag)
 n1 = n + 1
 150 j  = n1
 nx = z(i)
 ni = 1
 160 CALL READ (*410,*420,geomp,buf,3,0,flag)
 IF (intern) GO TO 300
 170 z(j  ) = buf(2)
 z(j+1) = buf(3)
 j  = j + 2
 IF (j >= buf3) GO TO 430
 ni = ni + 1
 IF (ni <= nx) GO TO 160
 nx = j - n1
 CALL sort (0,0,2,1,z(n1),nx)
 
!     TEST FOR UNIQUENESS OF POINT AND TEMPERATURE
 
 khi = j  - 1
 klo = n1 + 2
 k   = j
 IF (klo >= khi) GO TO 210
 k   = klo
 DO  j = klo,khi,2
   IF (z(j) /= z(j-2)) GO TO 190
   
!     NOT FATAL IF SAME TEMPERATURE
   
   IF (z(j+1) /= z(j-1)) nogo = nogo + 1
   IF (intern) CYCLE
   CALL page2 (2)
   WRITE  (nout,180) ufm,z(j-1),z(j+1),z(j)
   180 FORMAT (a23,' 2100, TEMPERATURE SPECIFIED HAS ',1P,e10.3,4H AND,  &
       1P,e10.3,' FOR GRID',i9)
   CYCLE
   
!     VALID TEMPERATURE
   
   190 z(k  ) = z(j  )
   z(k+1) = z(j+1)
   k  = k + 2
 END DO
 
 210 nx = k - n1
 CALL WRITE (ifile,z(n1),nx,1)
 i  = i + 1
 IF (i <= n) GO TO 150
 
!     NOW DO SAME AS ABOVE WITH OUTPUT IN INTERNAL INDEX NOTATION.
 
 IF (nogo /= 0) CALL mesage (-61,nogo,0)
 IF (intern) GO TO 220
 CALL bckrec (geomp)
 intern = .true.
 ifile  = ett
 i = isave
 GO TO 140
 
!     NOW APPEND ENTIRE ETT FILE TO GPTT FILE
 
 220 FILE = ett
 CALL CLOSE (ett,clsrew)
 CALL OPEN  (*400,ett,z(buf3),rdrew)
 230 CALL READ  (*250,*240,ett,z,buf3-1,0,flag)
 CALL WRITE (gptt,z,buf3-1,0)
 GO TO 230
 240 CALL WRITE (gptt,z,flag,1)
 GO TO 230
 250 CALL CLOSE (gptt,clsrew)
 CALL CLOSE (ett,clsrew)
 260 CALL CLOSE (geomp,clsrew)
 GO TO 460
 
!     NO TEMP CARDS PRESENT. IF NO DEFAULT CARDS, NO GPTT.
!     OTHERWISE, GPTT IS COMPRISED ONLY OF DEFAULT TEMPERATURES.
!     WRITE THE SET IDS AND DEFAULT TEMPS IN THE HEADER RECORD.
 
 270 IF (nodef == 0) GO TO 260
 FILE = gptt
 CALL OPEN  (*400,gptt,z(buf2),wrtrew)
 CALL fname (gptt,buf)
 CALL WRITE (gptt,buf,2,0)
 FILE = ett
 CALL OPEN  (*400,ett,z(buf3),wrtrew)
 CALL fname (ett,buf)
 CALL WRITE (ett,buf,2,0)
 buf(3) =  0
 DO  k = itempd,ntempd,2
   buf(1) = z(k  )
   buf(2) = z(k+1)
   CALL WRITE (gptt,buf,3,0)
 END DO
 CALL WRITE (ett ,buf,3,0)
 CALL WRITE (gptt,0,0,1)
 CALL WRITE (ett ,0,0,1)
 GO TO 220
 
!     INTERNAL BINARY SEARCH ROUTINE.
 
 300 klo = 1
 khi = kn
 310 k = (klo+khi+1)/2
 320 IF (buf(2)-z(2*k-1) < 0.0) THEN
   GO TO   330
 ELSE IF (buf(2)-z(2*k-1) == 0.0) THEN
   GO TO   390
 ELSE
   GO TO   340
 END IF
 330 khi = k
 GO TO 350
 340 klo = k
 350 IF (khi -klo-1 < 0) THEN
   GO TO   440
 ELSE IF (khi -klo-1 == 0) THEN
   GO TO   360
 ELSE
   GO TO   310
 END IF
 360 IF (k == klo) GO TO 370
 k = klo
 GO TO 380
 370 k = khi
 380 klo = khi
 GO TO 320
 390 buf(2) = z(2*k)
 GO TO 170
 
!     FATAL ERROR MESAGES
 
 400 j = -1
 GO TO 450
 410 j = -2
 GO TO 450
 420 j = -3
 GO TO 450
 430 j = -8
 GO TO 450
 440 CALL mesage (-30,9,buf)
 450 CALL mesage (j,FILE,nam)
 
 460 RETURN
END SUBROUTINE gp3b
