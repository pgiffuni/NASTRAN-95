SUBROUTINE gp3c
     
!     GP3C EXECUTES ONLY IF PLOAD2 AND/OR PLOAD3 CARDS ARE PRESENT. ITS
!     FUNCTION IS TO --
!     (1) PROCESS PLOAD2 CARDS SO THAT THEIR FORMAT IS IDENTICAL TO
!         PLOAD CARDS.  IF A PLOAD RECORD EXISTS ON GEOM3, PLOAD2 DATA
!         IS APPENDED TO THE DATA, SORTED, AND ALL RESULTING PLOAD DATA
!         IS WRITTEN ON SCR2.
!     (2) PROCESS PLOAD3 CARDS SO THAT ALL PRESSURES APPLIED TO AN ISO-
!         PARAMETRIC SOLID ARE GATHERED IN ONE ENTRY AND SORTED BY THE
!         FACE NUMBER TO WHICH THE PRESSURE IS APPLIED.  THE SORTED
!         PRESSURES AND GRID POINT NUMBERS FOR EACH ELEMENT ARE WRITTEN
!         ON SCR2.
 
 
 EXTERNAL        andf
 INTEGER :: geom3 ,geom2 ,scr2  ,buf1  ,buf2  ,buf   ,clsrew,  &
     wrtrew,rdrew ,elem  ,FILE  ,nam(2),pload2,z     ,  &
     cardid,pload3,pl2   ,pl3   ,pld(3),d1    ,d2    ,  &
     andf  ,faces(6,12)  ,pl3err(14)
 REAL :: rz(1),p(12)
 COMMON /gp3com/ geom3 ,eqexin,geom2 ,slt   ,gptt  ,scr1  ,scr2   ,  &
     buf1  ,buf2  ,buf(50)      ,cardid(60)   ,idno(30)  &
     , carddt(60)   ,mask(60)     ,STATUS(60)   ,ntypes ,  &
     ipload,igrav ,pload2(2)    ,load(2)      ,nopld2 ,  &
     temp(2)      ,tempd(2)     ,tempp1(2)    ,  &
     tempp2(2)    ,tempp3(2)    ,temprb(2)    ,buf3   , pload3(2)    ,ipld3
 COMMON /names / rd    ,rdrew ,wrt   ,wrtrew,clsrew
 COMMON /zzzzzz/ z(1)
 COMMON /gpta1 / nelem ,last  ,incr  ,elem(1)
 COMMON /two   / two(32)
 COMMON /system/ sysbuf,nout
 EQUIVALENCE     (rz(1),z(1))
 
!  FACE                 IHEX1             IHEX2             IHEX3
!   NO                D1      D2        D1      D2        D1      D2
 DATA    faces/   1,      3,        1,      5,        1,      7,  &
     2,      4,        3,      7,        4,     10,  &
     1,      6,        1,     15,        1,     24,  &
     2,      5,        3,     13,        4,     21,  &
     2,      7,        3,     17,        4,     27,  &
     3,      6,        5,     15,        7,     24,  &
     3,      8,        5,     19,        7,     30,  &
     4,      7,        7,     17,       10,     27,  &
     1,      8,        1,     19,        1,     30,  &
     4,      5,        7,     13,       10,     21,  &
     5,      7,       13,     17,       21,     27,  &
     6,      8,       15,     19,       24,     30/
 
 DATA   n3304,n3305,pl3err/4H3304, 4H3305, 4H0***, 4H use, 4HR fa,  &
     4HTAL , 4HMESS, 4HAGE , 4H330*, 4H, pl,  &
     4HOAD3, 4H car, 4HD fr, 4HOM l, 4HOAD , 4HSET /
 DATA   nam  /    4HGP3C,4H      /
 
!     CHECK TRAILER BITS FOR PRESENCE OF PLOAD2 AND PLOAD3 CARDS.
!     IF NONE EXIST, RETURN.  OTHERWISE, BRANCH AND INITIALIZE TO
!     PROCESS ONE OF THESE CARD TYPES.
 
 nogo = 0
 pl2  = 0
 pl3  = 0
 j    = (pload2(2)-1)/16
 k    = pload2(2) - 16*j
 IF (andf(buf(j+2),two(k+16)) /= 0) pl2 = 1
 j    = (pload3(2)-1)/16
 k    = pload3(2) - 16*j
 IF (andf(buf(j+2),two(k+16)) /= 0) pl3 = 1 - 2*pl2
 FILE = scr2
 IF (pl2-pl3 /= 0) CALL OPEN (*210,scr2,z(buf2),wrtrew)
 IF (pl2 == 0) GO TO 15
 nopld2 = 1
 pld(1) = pload2(1)
 pld(2) = pload2(2)
 pld(3) = 24
 incrd  = 3
 incl   = 6
 idl    = 2
 GO TO 10
 15 IF (pl3 == 0) GO TO 196
 nopld2 = nopld2 + 2
 pld(1) = pload3(1)
 pld(2) = pload3(2)
 pld(3) = 255
 incrd  = 5
 incl   = 39
 idl    = 1
 
!     READ PLOAD2 OR PLOAD3 CARDS INTO CORE IN AN EXPANDED FORMAT.
!     SET THE SET ID NEGATIVE TO INDICATE THE CARD IS NOT YET CONVERTED.
 
 10 i = 1
 FILE = geom3
 CALL preloc (*210,z(buf1),geom3)
 CALL locate (*230,z(buf1),pld,flag)
 IF (pl2 /= 1) GO TO 20
 pld(1) = cardid(ipload  )
 pld(2) = cardid(ipload+1)
 20 CALL READ (*220,*30,geom3,z(i),incrd,0,flag)
 z(i) = -z(i)
 IF (pl2 ==  1) GO TO 29
 IF (i < incl) GO TO 25
 DO  j = 2,i,incl
   k = j
   IF (z(j) /= z(i+2)) CYCLE
   IF (z(j-1) == z(i)) GO TO 22
 END DO
 25 p(1)   = rz(i+1)
 z(i+1) = z(i+2)
 rz(i+2)= p(1)
 z(i+14)= z(i+3)
 z(i+15)= z(i+4)
 z(i+3) =-1
 GO TO 29
 22 j = k + 2
 23 IF (z(j) == -1) GO TO 24
 j = j + 1
 IF (j <= k+12) GO TO 23
 GO TO 25
 24 rz(j) = rz(i+1)
 IF (j < k+12) z(j+1) = -1
 j = k + 15 + 2*(j-k-2)
 z(j  ) = z(i+3)
 z(j+1) = z(i+4)
 GO TO 20
 29 z(i+incl-1) = 0
 i = i+incl
 IF (i < buf2) GO TO 20
 CALL mesage (-8,0,nam)
 30 CALL CLOSE (geom3,clsrew)
 npld2 = i - incl
 nwds  = i - 1
 
!     POSITION TO FIRST DATA RECORD ON GEOM2.
 
 FILE = geom2
 CALL OPEN (*130,geom2,z(buf1),rdrew)
 CALL fwdrec (*220,geom2)
 
!     READ 3-WORD RECORD ID. LOOK FOR ID IN ELEM TABLE.
!     IF NOT THERE, SKIP RECORD.
!     IF PROCESSING PLOAD2, AND NOT A TWO-DIMENSIONAL ELEMENT, SKIP REC.
!     IF PROCESSING PLOAD3, AND NOT AN ISOPARAMETRIC ELEMENT, SKIP REC.
!     OTHERWISE,  INITIALIZE PARAMETERS.
 
 50 CALL READ (*130,*50,geom2,buf,3,0,flag)
 DO  i = 1,last,incr
   IF (buf(1) == elem(i+3)) GO TO 80
 END DO
 70 CALL fwdrec (*220,geom2)
 GO TO 50
 80 ngps  = elem(i+9)
 itype = elem(i+2)
 
!   . IF ELEMENT TYPE IS 68 (QUADTS) THEN USE FIRST FOUR  GRID POINTS
!   . IF ELEMENT TYPE IS 69 (TRIATS) THEN USE FIRST THREE GRID POINTS
 
 IF (itype == 68 .OR. itype == 69) ngps = ngps/2
 IF (pl2 == 1 .AND. (ngps < 3 .OR. ngps > 4)) GO TO 70
 IF (pl3 == 1 .AND. (itype < 65 .OR. itype > 67)) GO TO 70
 itype  = 2*(itype-64) - 1
 nwdect = elem(i+5)
 j1 = elem(i+12)
 j2 = j1 + ngps - 1
 
!     READ EACH ELEMENT IN RECORD. LOOK FOR ELEMENT ID MATCH IN PLOAD2
!     OR PLOAD3 LIST.  IF FOUND, SET THE SET ID POSITIVE TO INDICATE
!     ENTRY IS CONVERTED.
 
 90 CALL READ (*220,*50,geom2,buf,nwdect,0,flag)
 DO  i = 1,npld2,incl
   IF (z(i) > 0) CYCLE
   IF (z(i+idl) /= buf(1)) CYCLE
   z(i) = -z(i)
   ix   = i
   IF (pl3 == 1) GO TO 300
   
!     PLACE GRID POINT NUMBERS FROM ELEMENT CARD IN PLOAD2 ENTRY TO
!     MAKE IT LOOK LIKE PLOAD CARD.
   
   DO  j = j1,j2
     z(ix+2) = buf(j)
     ix = ix + 1
   END DO
   CYCLE
   
!     FIND THE DIAGONALS ON THE PLOAD3 CARD ON THE ELEMENT CARD TO
!     DETERMINE THE FACES TO WHICH THE PRESSURES ARE APPLIED.  SORT
!     THE PRESSURES BY FACE NUMBER AND APPEN+ THE GRID POINT NUMBERS
!     FROM THE ELEMENT CARD TO THE PLOAD3 ENTRY.
   
   300 np = 0
   DO  j = 1,12
     IF (z(i+j+1) == -1) EXIT
     np = np + 1
     p(j) = rz(i+j+1)
   END DO
   315 DO  j = 1,6
     rz(i+j)  = 0.0
   END DO
   DO  j = 1,np
     k   = i + 14 + 2*(j-1)
     id1 = z(k  )
     id2 = z(k+1)
     DO  k = j1,j2
       IF (id1 == buf(k)) GO TO 324
     END DO
     GO TO 335
     324 id1 = k - j1 + 1
     DO  k = j1,j2
       IF (id2 == buf(k)) GO TO 328
     END DO
     GO TO 335
     328 id2 = k - j1 + 1
     d1  = MIN0(id1,id2)
     d2  = MAX0(id1,id2)
     DO  k = 1,12
       nface = (k+1)/2
       IF (d1 /= faces(itype  ,k)) CYCLE
       IF (d2 == faces(itype+1,k)) GO TO 340
     END DO
     335 nogo = 1
     pl3err(7) = n3305
     WRITE (nout,420) pl3err,z(i),buf(1)
     CYCLE
     340 rz(i+nface) = rz(i+nface)+p(j)
   END DO
   ix = ix + 7
   DO  j = j1,j2
     z(ix) = buf(j)
     ix = ix + 1
   END DO
   IF (ix+1-i > 39) CYCLE
   k = i + 38
   DO  j = ix,k
     z(j) = 0
   END DO
 END DO
 GO TO 90
 
!     HERE WHEN END-OF-FILE ON GEOM2 IS ENCOUNTERED.
!     MAKE SURE ALL PLOAD2 OR PLOAD3 ENTRIES HAVE BEEN CONVERTED.
 
 130 CALL CLOSE (geom2,clsrew)
 DO  i = 1,npld2,incl
   IF (z(i) > 0) CYCLE
   nogo = 1
   buf(1) = -z(i)
   buf(2) =  z(i+idl)
   IF (pl2 == 1) CALL mesage (30,105,buf)
   pl3err(7) = n3304
   IF (pl3 == 1) WRITE (nout,410) pl3err,buf(1),buf(2)
 END DO
 IF (nogo /= 0) CALL mesage (-61,0,0)
 IF (pl3  == 1) GO TO 190
 
!     LOCATE PLOAD RECORD ON GEOM3. IF PRESENT, READ PLOAD DATA INTO
!     CORE (AFTER PLOAD2 DATA) AND SORT COMBINED DATA ON SET ID.
 
 CALL preloc (*210,z(buf1),geom3)
 CALL locate (*180,z(buf1),cardid(ipload),flag)
 i = npld2 + 6
 160 CALL READ (*220,*170,geom3,z(i),6,0,flag)
 i = i + 6
 IF (i < buf2) GO TO 160
 CALL mesage (-8,0,nam)
 170 npld2 = i - 6
 nwds  = i - 1
 CALL sort (0,0,6,1,z,nwds)
 180 CALL CLOSE (geom3,clsrew)
 
!     WRITE DATA ON SCR2, SET FLAG TO INDICATE AND RETURN.
 
 190 CALL WRITE (scr2,pld,3,0)
 CALL WRITE (scr2,z,nwds,1)
 IF (pl2 /= 1) GO TO 196
 195 pl2 = -pl2
 pl3 = -pl3
 GO TO 15
 196 CALL CLOSE (scr2,clsrew)
 RETURN
 
!     ERROR MESSAGES.
 
 200 CALL mesage (n,FILE,nam)
 210 n = -1
 GO TO 200
 220 n = -2
 GO TO 200
 
!     ABNORMAL RETURN.
 
 230 IF (pl3 < 0) GO TO 195
 CALL CLOSE (geom3,clsrew)
 RETURN
 
!     PLOAD3 CARD ERRORS
 
 410 FORMAT (14A4,i9,' REFERENCES MISSING OR NON-ISOPARAMETRIC ELEMENT' ,i9)
 420 FORMAT (14A4,i9,' HAS INVALID GRID POINT NUMBERS FOR ELEMENT',i9)
END SUBROUTINE gp3c
