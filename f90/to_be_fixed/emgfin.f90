SUBROUTINE emgfin
     
!     THIS ROUTINE OF THE -EMG- MODULE WRAPS UP THE WORK OF THE MODULE.
 
 LOGICAL :: error, heat, linear
 INTEGER :: cls, clsrew, rdrew, wrtrew, date, mcb(7), precis,  &
     est, dictn, flags, scr3, scr4, vafile, sub(2), ix(6), sil(32)
 REAL :: inpi(10), z(2), rx(200), corey(201)
 COMMON /BLANK / nokmb(3), nok4gg, cmass, dummy(11), volume, surfac
 COMMON /hmatdd/ skp(4), linear
 COMMON /emgprm/ icore, jcore, ncore, dum12(12), flags(3), precis,  &
     error, heat, icmbar, lcstm, lmat, lhmat
 COMMON /names / rd, rdrew, wrt, wrtrew, clsrew, cls
 COMMON /emgfil/ est, cstm, mpt, dit, geom2, matrix(3), dictn(3)
 COMMON /output/ head(96)
 COMMON /system/ ibuf, nout, skip6(6), nlpp, skip2(2), line, sk1p2(2), date(3)
 COMMON /machin/ mach
 COMMON /zzzzzz/ corex(1)
 EQUIVALENCE     (corex(1),corey(1),rx(1),ix(1)), (z(1),corey(201))
 DATA    vafile, scr4  / 207,    304           /
 DATA    d2,     d3    / 4H2-d , 4H3-d         /
 DATA    sub,    scr3  / 4HEMGF, 4HIN   , 303  /
 DATA    inpi  / 4HINPT, 4HINP1, 4HINP2, 4HINP3, 4HINP4,  &
     4HINP5, 4HINP6, 4HINP7, 4HINP8, 4HINP9/
 
!     CLOSE ALL FILES, EXCEPT SCR4
 
 DO  i = 1,3
   nokmb(i) = -flags(i) - 1
   IF (flags(i)+1 == 0) flags(i) = 0
   IF (flags(i)   == 0) nokmb(i) =-1
   CALL CLOSE (matrix(i),clsrew)
   CALL CLOSE (dictn(i),clsrew)
 END DO
 CALL CLOSE (est,clsrew)
 
!     HEAT ONLY - SET NONILINEAR FLAG BASED ON VALUE PREVIOUSLY SET BY
!     HMAT ROUTINE
 
 IF (heat .AND. .NOT.linear) nokdgg = +1
 IF (heat .AND.      linear) nokdgg = -1
 
!  WRITE TRAILERS FOR FILES PREPARED.
 
 IF (error) GO TO 340
 DO  i = 1,3
   
!     PRECISION IS STORED IN FIRST DATA WORD OF TRAILER.
   
   IF (flags(i) == 0) CYCLE
   mcb(1) = matrix(i)
   CALL rdtrl (mcb)
   IF (mcb(1) > 0) THEN
     GO TO    20
   ELSE
     GO TO    40
   END IF
   20 mcb(2) = precis
   mcb(3) = 0
   mcb(4) = 0
   mcb(5) = 0
   mcb(6) = 0
   mcb(7) = 0
   CALL wrttrl (mcb)
   
   mcb(1) = dictn(i)
   CALL rdtrl (mcb)
   IF (mcb(1) > 0) THEN
     GO TO    30
   ELSE
     GO TO    40
   END IF
   30 mcb(2) = precis
   mcb(3) = 0
   mcb(4) = 0
   mcb(5) = 0
   mcb(6) = 0
   mcb(7) = 0
   CALL wrttrl (mcb)
 END DO
 IF (volume <= 0.0 .AND. surfac <= 0.0) GO TO 330
 
!     COMPUTE AND PRINT VOLUMES AND SURFACE AREAS FOR THE 2-D AND 3-D
!     ELEM. IF USER REQUESTED VIA PARAM CARD.
 
 CALL CLOSE (scr4,clsrew)
 ibuf1 = icore + 200
 ibuf2 = ibuf1 + ibuf
 CALL OPEN (*350,scr4,z(ibuf1),rdrew)
 tvol2 = 0.0
 tvol3 = 0.0
 tmas2 = 0.0
 tmas3 = 0.0
 nrec  = 0
 line  = nlpp
 inp   = 0
 
!     CHECK ANY REQUEST TO SAVE VOLUME AND AREA COMPUTATIONS ON OUTPUT
!     FILE SET INP TO APPROPRIATE VALUE IF IT IS AN INPI FILE
 
 mcb(1) = vafile
 CALL rdtrl (mcb(1))
 IF (mcb(1) <= 0) GO TO 70
 CALL fname (vafile,z(1))
 DO  i = 1,10
   IF (z(1) == inpi(i)) GO TO 60
 END DO
 GO TO 65
 60 inp = i + 13
 IF (inp == 14 .AND. mach == 2) inp = 24
 vafile = scr3
 mcb(1) = scr3
 65 CALL OPEN (*360,vafile,z(ibuf2),wrtrew)
 CALL WRITE (vafile,z(1),    2,0)
 CALL WRITE (vafile,head(1),96,0)
 CALL WRITE (vafile,date(1), 3,1)
 nrec = 1
 70 CALL READ (*210,*90,scr4,rx,201,1,i)
 WRITE  (nout,80)
 80 FORMAT (' *** WARNING,   RX TOO SMALL IN EMGFIN ***')
 90 ngpt = ix(6)
 IF (ngpt <    3) GO TO 70
 IF (line < nlpp) GO TO 130
 line = 5
 CALL page1
 WRITE  (nout,100) (i,i=1,6)
 100 FORMAT (17X,'V O L U M E S,  M A S S E S,  A N D  S U R F A C E ',  &
     ' A R E A S  O F  2-  A N D  3-  D  E L E M E N T S',  &
     ///10X,7HELEMENT,8X,3HEID,8X,6HVOLUME,7X,4HMASS,1X,  &
     6(3X,7HSURFACE,i2), /10X, 29(4H----),/)
 IF (volume <= 0.0) WRITE (nout,110)
 IF (surfac <= 0.0) WRITE (nout,120)
 IF (volume <= 0.0 .OR. surfac <= 0.0) line=line+2
 110 FORMAT (10X,42H(no mass AND volume computation requested),/)
 120 FORMAT (10X,39H(no surface area computation requested),/)
 130 l = 5
 
!     ENTRIES IN RX ARRAY, AS SAVED IN SCR4 BY KTRIQD,KTETRA,IHEXI,
!     EMGPRO
!        RX( 1),RX(2) = ELEMENT BCD NAME
!        IX( 3) = ELEMENT ID
!        RX( 4) = VOLUME (SOLID), OR THICKNESS (PLATE)
!        RX( 5) = TOTAL MASS (SOLID), OR DENSITY (PLATE)
!        IX( 6) = NO. OF GRID POINTS, = NGPT
!        IX(7)...IX(6+NGPT) = SIL OF THE GRID POINTS
!        RX( 7+NPGT) = CID OF 1ST GRID POINT
!        RX( 8+NPGT) = X COORD. OF 1ST GRID POINT
!        RX( 9+NPGT) = Y COORD. OF 1ST GRID POINT
!        RX(10+NPGT) = Z COORD. OF 1ST GRID POINT
!        IX(11+NPGT...) = REPEAT FOR OTHER GRID POINTS
 
!     CALL SFAREA TO COMPUTE AREAS, 6 VALUES ARE RETURNED IN RX(6...11)
!     AND NO. OF SURFACES IN NGPT
!     VOLUME AND MASS ARE ALSO COMPUTED FOR THE PLATE ELEMENTS.
 
 DO  i = 1,ngpt
   sil(i) = ix(6+i)
 END DO
 ln = ngpt
 CALL sfarea (ln,rx,ix(ngpt+7))
 l = 5
 IF (surfac > 0.0) l = 5 + ln
 IF (volume > 0.0) WRITE (nout,160) (ix(i),i=1,3),(rx(i),i=4,l)
 IF (volume <= 0.0) WRITE (nout,170) (ix(i),i=1,3),(rx(i),i=6,l)
 160 FORMAT (10X,2A4,i10, 2X,8E12.4)
 170 FORMAT (10X,2A4,i10,26X,6E12.4)
 line = line + 1
 
 IF (nrec == 0) GO TO 190
 nrec = nrec + 1
 ix(5) = (ln*100) + ngpt
 CALL WRITE (vafile,rx(1),l,0)
 n4 = ngpt*4
 n7 = n4 + 7
 j  = 1
 DO  i = 7,n7,4
   ix(i) = sil(j)
   j = j + 1
 END DO
 CALL WRITE (vafile,rx(ngpt+7),n4,1)
 
!     A RECORD IS SAVED IN VAFILE FOR EACH ELEM., HAVING THE FOLLOWING
!     DATA
 
!        WORDS  1,2 ELEMENT BCD NAME
!                 3 ELEMENT ID, INTEGER
!                 4 VOLUME (3-D ELEMS), ZERO (2-D ELEMS), REAL
!                 5 (NO. OF SURFACES, N)*100 + NO. OF GRID PTS, INTEGER
!                 6 AREA OF FIRST SURFACE, REAL
!        7 THRU 5+N REPEAT FOR N SURFACES, REAL
!             5+N+1 SIL OF FIRST GRID POINT, INTEGER
!         5+N+2,3,4 X,Y,Z COORDINATES OF THE FIRST GRID POINT, REAL
!               ... REPEAT LAST 4 WORDS FOR OTHER GRID POINTS, REAL
 
 190 IF (volume <= 0.0) GO TO 70
 IF (ngpt   >   1) GO TO 200
 tvol2 = tvol2 + rx(4)
 tmas2 = tmas2 + rx(5)
 GO TO 70
 200 tvol3 = tvol3 + rx(4)
 tmas3 = tmas3 + rx(5)
 GO TO 70
 210 CALL CLOSE (scr4,clsrew)
 IF (nrec == 0) GO TO 230
 CALL CLOSE (vafile,clsrew)
 mcb(2) = nrec
 DO  i = 3,7
   mcb(i) = 0
 END DO
 CALL wrttrl (mcb(1))
 230 IF (volume <= 0.0) GO TO 330
 IF (tvol2  > 0.0) WRITE (nout,240) tvol2,tmas2,d2
 IF (tvol3  > 0.0) WRITE (nout,240) tvol3,tmas3,d3
 240 FORMAT (/6X,24H* total volume AND mass=,2E12.4,3H  (,a4, 9HELEMENTS))
 IF (nrec <= 0) GO TO 330
 
!     IF OUTPUT FILE REQUESTED BY USER IS AN INPI FILE, COPY FROM VAFILE
!     TO INPI, A FORTRAN WRITTEN BINARY FILE
 
 IF (inp == 0) GO TO 280
 CALL OPEN (*360,vafile,z(ibuf2),rdrew)
 260 CALL READ (*280,*270,vafile,z(1),ibuf2,1,j)
 GO TO 370
 270 WRITE (inp) (z(i),i=1,j)
 GO TO 260
 280 CALL CLOSE (vafile,clsrew)
 CALL fname (207,z(1))
 WRITE  (nout,300) z(1),z(2),date
 300 FORMAT ('0*** VOLUMES AND EXTERNAL SURFACE AREAS WERE SAVED IN ',  &
     'OUTPUT FILE ',2A4,4H on ,i2,1H/,i2,3H/19,i2)
 IF (inp == 0) WRITE (nout,310)
 IF (inp /= 0) WRITE (nout,320) inp
 310 FORMAT (1H+,91X,21H(a gino written FILE))
 320 FORMAT (1H+,91X,28H(a fortran binary FILE, UNIT,i3,1H))
 330 volume = 0.0
 surfac = 0.0
 RETURN
 
 340 IF (volume /= 0. .OR. surfac /= 0.) CALL CLOSE (scr4,clsrew)
 GO TO 330
 350 CALL mesage (-1,scr4,sub)
 360 j = -1
 GO TO 380
 370 j = -8
 380 CALL mesage (j,vafile,sub)
 RETURN
END SUBROUTINE emgfin
