SUBROUTINE ssgslt (slt,newslt,est)
     
!     THIS SUBROUTINE OF THE SSG1 MODULE COPIES THE SLT TO ANOTHER
!     FILE.  IN THE COPYING PROCESS ANY -QVOL-, -QBDY1-, -QBDY2-, OR
!     -QVECT- EXTERNAL LOAD TYPE DATA FOUND WILL BE ALTERED SO AS TO
!     REPLACE THEIR ELEMENT ID REFERENCES WITH THE APPROPRIATE SILS, AND
!     MISC. CONSTANTS.  THE EXTERNAL LOADS WILL BE PREPARED AS USUAL FOR
!     THESE AND OTHER LOAD CARD TYPES VIA SUBROUTINE EXTERN.
 
 
 INTEGER, INTENT(IN)                      :: slt
 INTEGER, INTENT(IN)                      :: newslt
 INTEGER, INTENT(IN)                      :: est
 IMPLICIT INTEGER (a-z)
 LOGICAL :: any,nogo,bgcore,bgopen
 REAL :: area,hc1,hc2,hc3,q0,piovr4,xx,yy,zz, rbuf(50),rz(1),recpt(100)
 INTEGER :: mcb(7),ecpt(100),buf(50),subr(2),TYPE(25,4)
 CHARACTER (LEN=25) :: sfm
 CHARACTER (LEN=29) :: uim
 CHARACTER (LEN=25) :: uwm
 CHARACTER (LEN=23) :: ufm
 COMMON /xmssg / ufm,uwm,uim,sfm
 COMMON /BLANK / nrowsp
 COMMON /gpta1 / nelem,last,incr,NE(1)
 COMMON /system/ ksystm(65)
 COMMON /names / rd,rdrew,wrt,wrtrew,clsrew,cls
 COMMON /zzzzzz/ z(1)
 EQUIVALENCE     (ksystm(1),sysbuf  ), (ksystm(2),outpt  ),  &
     (ecpt(1)  ,recpt(1)), (buf(1)   ,rbuf(1)), (z(1)     ,rz(1)   )
 DATA    subr  / 4HSSGS,4HLT    /    , noeor,eor/ 0, 1 /
 DATA    bgpdt / 102/
 DATA    piovr4/ 0.7853981634E0 /
 DATA    cbar  / 34 /
 DATA    crod  /  1 /
 DATA    conrod/ 10 /
 DATA    ctube /  3 /
 DATA    ctrmem/  9 /
 DATA    ctria1/  6 /
 DATA    ctria2/ 17 /
 DATA    cqdmem/ 16 /
 DATA    cqdmm1/ 62 /
 DATA    cqdmm2/ 63 /
 DATA    cquad4/ 64 /
 DATA    ctria3/ 83 /
 DATA    cquad1/ 19 /
 DATA    cquad2/ 18 /
 DATA    ctrirg/ 36 /
 DATA    ctrprg/ 37 /
 DATA    ctetra/ 39 /
 DATA    cwedge/ 40 /
 DATA    chexa1/ 41 /
 DATA    chexa2/ 42 /
 DATA    chbdy / 52 /
 DATA    cihex1/ 65 /
 DATA    cihex2/ 66 /
 DATA    cihex3/ 67 /
 
 DATA    ntypes/ 25 /
 
!          SLT            NEWSLT         FLAG FOR       DATA
!          WORDS-IN       WORDS-OUT      SPEC-PROC      CORE-LOCAT
!          ==========     ==========     ==========     ==========
!     FORCE
 DATA TYPE( 1,1)/ 6/,TYPE( 1,2)/ 6/,TYPE( 1,3)/ 0/,TYPE( 1,4)/ 0/
!     MOMENT
 DATA TYPE( 2,1)/ 6/,TYPE( 2,2)/ 6/,TYPE( 2,3)/ 0/,TYPE( 2,4)/ 0/
!     FORCE1
 DATA TYPE( 3,1)/ 4/,TYPE( 3,2)/ 4/,TYPE( 3,3)/ 0/,TYPE( 3,4)/ 0/
!     MOMNT1
 DATA TYPE( 4,1)/ 4/,TYPE( 4,2)/ 4/,TYPE( 4,3)/ 0/,TYPE( 4,4)/ 0/
!     FORCE2
 DATA TYPE( 5,1)/ 6/,TYPE( 5,2)/ 6/,TYPE( 5,3)/ 0/,TYPE( 5,4)/ 0/
!     MOMNT2
 DATA TYPE( 6,1)/ 6/,TYPE( 6,2)/ 6/,TYPE( 6,3)/ 0/,TYPE( 6,4)/ 0/
!     SLOAD
 DATA TYPE( 7,1)/ 2/,TYPE( 7,2)/ 2/,TYPE( 7,3)/ 0/,TYPE( 7,4)/ 0/
!     GRAV
 DATA TYPE( 8,1)/ 5/,TYPE( 8,2)/ 5/,TYPE( 8,3)/ 0/,TYPE( 8,4)/ 0/
!     PLOAD
 DATA TYPE( 9,1)/ 5/,TYPE( 9,2)/ 5/,TYPE( 9,3)/ 0/,TYPE( 9,4)/ 0/
!     RFORCE
 DATA TYPE(10,1)/ 6/,TYPE(10,2)/ 6/,TYPE(10,3)/ 0/,TYPE(10,4)/ 0/
!     PRESAX
 DATA TYPE(11,1)/ 6/,TYPE(11,2)/ 6/,TYPE(11,3)/ 0/,TYPE(11,4)/ 0/
!     QHBDY
 DATA TYPE(12,1)/ 7/,TYPE(12,2)/ 7/,TYPE(12,3)/ 0/,TYPE(12,4)/ 0/
!     QVOL
 DATA TYPE(13,1)/ 2/,TYPE(13,2)/12/,TYPE(13,3)/ 1/,TYPE(13,4)/ 0/
!     QBDY1
 DATA TYPE(14,1)/ 2/,TYPE(14,2)/10/,TYPE(14,3)/ 1/,TYPE(14,4)/ 0/
!     QBDY2
 DATA TYPE(15,1)/ 5/,TYPE(15,2)/10/,TYPE(15,3)/ 1/,TYPE(15,4)/ 0/
!     QVECT
 DATA TYPE(16,1)/ 5/,TYPE(16,2)/19/,TYPE(16,3)/ 1/,TYPE(16,4)/ 0/
!     PLOAD3
 DATA TYPE(17,1)/38/,TYPE(17,2)/38/,TYPE(17,3)/ 0/,TYPE(17,4)/ 0/
!     PLOAD1
 DATA TYPE(18,1)/ 7/,TYPE(18,2)/ 7/,TYPE(18,3)/ 0/,TYPE(18,4)/ 0/
!     PLOADX
 DATA TYPE(19,1)/ 5/,TYPE(19,2)/ 5/,TYPE(19,3)/ 0/,TYPE(19,4)/ 0/
!     SPCFLD  (WORDS OUT IS A DUMMY VALUE-IT WILL REALLY BE 3*NROWSP)
 DATA TYPE(20,1)/ 5/,TYPE(20,2)/ 4/,TYPE(20,3)/ 0/,TYPE(20,4)/ 0/
!     CEMLOOP
 DATA TYPE(21,1)/12/,TYPE(21,2)/12/,TYPE(21,3)/ 0/,TYPE(21,4)/ 0/
!     GEMLOOP  ,BOTH INPUT AND OUTPUT ARE DUMMY.
 DATA TYPE(22,1)/ 5/,TYPE(22,2)/ 4/,TYPE(22,3)/ 0/,TYPE(22,4)/ 0/
!     MDIPOLE (OUTPUT VALUE IS A DUMMY)
 DATA TYPE(23,1)/ 9/,TYPE(23,2)/ 5/,TYPE(23,3)/ 0/,TYPE(23,4)/ 0/
!     REMFLUX    (OUTPUT VALUE IS A DUMMY)
 DATA TYPE(24,1)/ 5/,TYPE(24,2)/ 5/,TYPE(24,3)/ 0/,TYPE(24,4)/ 0/
!     PLOAD4
 DATA TYPE(25,1)/11/,TYPE(25,2)/11/,TYPE(25,3)/ 0/,TYPE(25,4)/ 0/
 
!                     SLT                         NEWSLT
!     CARD=TYPE       WORDS IN                    WORDS OUT
!     =========       ========                    =========
 
!     QVOL=13         1 = QV                      1 = NUM-POINTS(1 TO 8)
!                     2 = ELEMENT ID              2 = ELEMENT ID
!                                                 3 THRU 10 = 8 SILS
!                                                 11 = COEFICIENT
!                                                 12 = TYPE 1 = 1 DIMEN
!                                                           2 = 2 DIMEN
!                                                           3 = BELL-EL
!                                                           4 = SOLID
 
!     QBDY1=14        1 = Q0                      1 = TYPE (1 TO 5)
!                     2 = ELEMENT ID              2 = ELEMENT ID
!                                                 3 THRU  6 = 4 SILS
!                                                 7 THRU 10 = 4 COEFS.
 
!     QBDY2=15        1 = ELEMENT ID              1 = ELEMENT ID
!                     2 = Q01                     2 = TYPE (1 TO 5)
!                     3 = Q02                     3 THRU  6 = 4 SILS
!                     4 = Q03                     7 THRU 10 = 4 COEFS.
!                     5 = Q04
 
!     QVECT=16        1 = Q0                      1 THRU 4 = 4 SILS
!                     2 = E1                      5 = ELEMENT ID
!                     3 = E2                      6 = TYPE (1 TO 5)
!                     4 = E3                      7 THRU 10 = 4 COEFS.
!                     5 = ELEMENT ID              11 = E1
!                                                 12 = E2
!                                                 13 = E3
!                                                 14 THRU 16 = V1 VECTOR
!                                                 17 THRU 19 = V2 VECTOR
 
 
 
!     SPCFLD=20       1 = CID                     1 THRU 3*NROWSP=
!                     2 = HCX                     TOTAL HC VALUES AT
!                     3 = HCY                     THE GRID POINTS
!                     4 = HCZ
!                     5 = GRID ID OR -1
 
 
!     CEMLOOP=21                                  SAME AS FOR
!     GEMLOOP=22                                  SPCFLD
 
 
 
!     MDIPOLE=23      1  =CID                     SAME AS
!                     2-4=LOCATION OF DIPOLE      SPCFLD
!                     5-7=DIPOLE MOMENT
!                     8  =MIN. DISTANCE
!                     9  =MAX. DISTANCE
 
!     REMFLUX=24      SAME INPUT AS               1 THRU 3*(NO. OF
!                     SPCFLD EXCEPT               ELEMENTS)= TOTAL
!                     WORD 5 IS ELEMENT ID        REMANENT FLUX DENSITY
!                                                 FOR EACH ELEMNT IN
!                                                 ORDER ON EST
 
!     THE ELEMENT ID MUST REMAIN IN THE SAME LOCATION ON OUTPUT.
 
 
!     SET UP CORE AND BUFFERS. (PG BUFFER IS OPEN IN SSG1)
 
 bgcore=.false.
 bgopen=.false.
 nogo  = .false.
 buf1  = korsz(z) - 2*sysbuf - 2
 buf2  = buf1 - sysbuf - 2
 buf3  = buf2 - sysbuf - 2
 core  = buf3 - 1
 IF (core < 100) CALL mesage (-8,0,subr)
 
!     OPEN SLT, AND NEWSLT.  COPY HEADER RECORD ACROSS.
 
 CALL OPEN (*650,slt,z(buf1),rdrew)
 CALL OPEN (*660,newslt,z(buf2),wrtrew)
 CALL READ (*670,*10,slt,z,core,eor,iwords)
 CALL mesage (-8,0,subr)
 10 CALL fname (newslt,z)
 CALL WRITE (newslt,z,iwords,eor)
 
!     READ TRAILER OF SLT AND GET COUNT OF LOAD SET RECORDS.
 
 mcb(1) = slt
 CALL rdtrl (mcb)
 nrecs  = mcb(2)
 mcb(1) = newslt
 CALL wrttrl (mcb)
 
!     PROCESSING OF LOAD SET RECORDS IF ANY.
 
 IF (nrecs > 0) THEN
   GO TO    20
 ELSE
   GO TO   620
 END IF
 20 i = 1
 21 CONTINUE
 any  = .false.
 ncore= nrecs+2
 ifirst = 0
 
!     ZERO OUT EXISTANCE FLAGS
 
 DO  k = 1,ntypes
   TYPE(k,4) = 0
 END DO
 
!     READ CARD TYPE AND COUNT OF CARDS
 
 30 CALL READ (*670,*110,slt,buf,2,noeor,iwords)
 35 CONTINUE
 itype  = buf(1)
 entrys = buf(2)
 
!     CHECK FOR KNOWN TYPE
 
 IF (itype <= ntypes .AND. itype > 0) GO TO 50
 WRITE  (outpt,40) sfm,itype
 40 FORMAT (a25,' 3094, SLT LOAD TYPE',i9,' IS NOT RECOGNIZED.')
 CALL mesage (-61,0,subr)
 
!     CHECK FOR SPECIAL PROCESSING.
 
 50 incnt = TYPE(itype,1)
 
!     IF TYPE IS CEMLOOP,SPCFLD,MDIPOLE, OR GEMLOOP,GO TO 800 FOR
!     SPECIAL PROCESSING. GO TO 1000 FOR REMFLUX PROCESSING
 
 IF (itype >= 20 .AND. itype <= 23) GO TO 800
 IF (itype == 24) GO TO 1000
 outcnt = TYPE(itype,2)
 iflag  = TYPE(itype,3)
 
 IF (iflag > 0) THEN
   GO TO    60
 ELSE
   GO TO    90
 END IF
 
!     OK BRING DATA INTO CORE.
 
 60 jcore = ncore + 1
 TYPE(itype,4) = jcore
 z(jcore  ) = itype
 z(jcore+1) = entrys
 jcore = jcore + 2
 ncore = jcore + entrys*outcnt - 1
 IF (ncore > core) CALL mesage (-8,0,subr)
 DO  j = jcore,ncore
   z(j) = 1
 END DO
 kcore = jcore
 
!     READ IN THE LOAD ENTRIES.
 
 DO  j = jcore,ncore,outcnt
   CALL fread (slt,z(j),incnt,0)
 END DO
 id = 2
 IF (itype == 15) id = 1
 IF (itype == 16) id = 5
 CALL sort (0,0,outcnt,id,z(kcore),ncore-kcore+1)
 any = .true.
 GO TO 30
 
!     NO SPECIAL PROCESSING OF THIS LOAD TYPE THUS JUST COPY IT ACROSS.
 
 90 CALL WRITE (newslt,buf,2,noeor)
 DO  j = 1,entrys
   CALL fread (slt,buf,incnt,0)
   CALL WRITE (newslt,buf,outcnt,noeor)
 END DO
 GO TO 30
 
!     ALL DATA NOW IN CORE FOR THIS LOAD SET.
 
 110 IF (.NOT. any) GO TO 600
 
!     THE EST IS NOW PROCESSED FOR ELEMENT TYPES CHECKED BELOW.
 
 CALL gopen (est,z(buf3),rdrew)
 
!     READ ELEMENT TYPE
 
 120 CALL READ (*500,*690,est,eltype,1,noeor,iwords)
 IF (eltype == cbar  ) GO TO 140
 IF (eltype == crod  ) GO TO 150
 IF (eltype == conrod) GO TO 150
 IF (eltype == ctube ) GO TO 160
 IF (eltype == ctrmem) GO TO 170
 IF (eltype == ctria1) GO TO 180
 IF (eltype == ctria2) GO TO 170
 IF (eltype == ctria3) GO TO 175
 IF (eltype == cqdmem) GO TO 190
 IF (eltype == cqdmm1) GO TO 190
 IF (eltype == cqdmm2) GO TO 190
 IF (eltype == cquad1) GO TO 200
 IF (eltype == cquad2) GO TO 190
 IF (eltype == cquad4) GO TO 195
 IF (eltype == ctrirg) GO TO 210
 IF (eltype == ctrprg) GO TO 220
 IF (eltype == ctetra) GO TO 230
 IF (eltype == cwedge) GO TO 240
 IF (eltype == chexa1) GO TO 250
 IF (eltype == chexa2) GO TO 250
 IF (eltype == cihex1) GO TO 252
 IF (eltype == cihex2) GO TO 254
 IF (eltype == cihex3) GO TO 256
 IF (eltype == chbdy ) GO TO 360
 130 CALL fwdrec (*700,est)
 GO TO 120
 
!     BAR
 
 140 estwds = 42
 grid1  = 2
 points = 2
 iarea  = 17
 itype  = 1
 GO TO 260
 
!     ROD AND CONROD
 
 150 estwds = 17
 grid1  = 2
 points = 2
 iarea  = 5
 itype  = 1
 GO TO 260
 
!     TUBE
 
 160 estwds = 16
 grid1  = 2
 points = 2
 iarea  = 5
 itype  = 1
 GO TO 260
 
!     TRMEM AND TRIA2
 
 170 estwds = 21
 grid1  = 2
 points = 3
 iarea  = 7
 itype  = 2
 GO TO 260
 
!     TRIA3
 
 175 estwds = 39
 grid1  = 2
 points = 3
 iarea  = 7
 itype  = 2
 GO TO 260
 
!     TRIA1
 
 180 estwds = 27
 grid1  = 2
 points = 3
 iarea  = 7
 itype  = 2
 GO TO 260
 
!     QDMEM AND QUAD2
 
 190 estwds = 26
 grid1  = 2
 points = 4
 iarea  = 8
 itype  = 2
 GO TO 260
 
!     QUAD4
 
 195 estwds = 45
 grid1  = 2
 points = 4
 iarea  = 8
 itype  = 2
 GO TO 260
 
!     QUAD1
 
 200 estwds = 32
 grid1  = 2
 points = 4
 iarea  = 8
 itype  = 2
 GO TO 260
 
!     TRIRG
 
 210 estwds = 19
 grid1  = 2
 points = 3
 iarea  = 0
 itype  = 3
 GO TO 260
 
!     TRAPRG
 
 220 estwds = 24
 grid1  = 2
 points = 4
 iarea  = 0
 itype  = 3
 GO TO 260
 
!     TETRA
 
 230 estwds = 23
 grid1  = 3
 points = 4
 iarea  = 0
 itype  = 4
 GO TO 260
 
!     WEDGE
 
 240 estwds = 33
 grid1  = 3
 points = 6
 iarea  = 0
 itype  = 4
 GO TO 260
 
!     HEXA1 AND HEXA2
 
 250 estwds = 43
 grid1  = 3
 points = 8
 iarea  = 0
 itype  = 4
 GO TO 260
 
!     IHEX1
 
 252 estwds = 55
 grid1  = 3
 points = 8
 iarea  = 0
 itype  = 4
 GO TO 260
 
!     IHEX2 AND IHEX3 ARE NOT IMPLEMENTED DUE TO
!        1. ECPT ARRAY TOO SMALL IN THIS ROUTINE
!        2. QVOL ROUTINE CAN NOT HANDLE SOLID ELEMENTS HAVING MORE THAN
!           8 GRID POINTS
 
!     IHEX2
 
 254 estwds = 127
 grid1  = 3
 points = 20
 iarea  = 0
 itype  = 4
 GO TO 130
 
!     IHEX3
 
 256 estwds = 199
 grid1  = 3
 points = 32
 iarea  = 0
 itype  = 4
 GO TO 130
 
!     MISC. ELEMENTS OF EST FILE.  DO QVOL REFERENCES.
 
 260 IF (TYPE(13,4) == 0) GO TO 130
 idqvol = TYPE(13,4) + 3
 entrys = z(idqvol-2)
 iwords = 12
 j1 = idqvol
 j2 = j1 + entrys*iwords
 idptr = 2
 ASSIGN 280 TO iretrn
 270 CALL READ (*700,*120,est,ecpt,estwds,noeor,iflag)
 
!     LOOK FOR THIS ELEMENT ID AMONG QVOL DATA.
 
 CALL bisloc (*270,ecpt(1),z(idqvol),12,entrys,jpoint)
 
!     MATCH FOUND ON ID. COMPUTE ZERO POINTER TO ZERO WORD OF QVOL ENTRY
 
 INDEX = idqvol + jpoint - 3
 GO TO 710
 
!     IF COEFICIENT IS NOT AN INTEGER 1 ENTRY HAS BEEN ALTERED BEFORE.
 
 280 DO  INDEX = k1,k2,iwords
   IF (z(INDEX+11) == 1) GO TO 300
   WRITE  (outpt,290) uwm,eltype,ecpt(1),z(i+2)
   290 FORMAT (a25,' 3095, ELEMENT TYPE',i9,' WITH ID =',i9,  &
       ', REFERENCED BY A QVOL CARD IN LOAD SET',i9,1H,, /5X,  &
       'IS NOT BEING USED FOR INTERNAL HEAT GENERATION IN THIS ',  &
       'LOAD SET BECAUSE ANOTHER ELEMENT TYPE WITH THE SAME ID',  &
       /5X,'HAS ALREADY BEEN USED.')
   CYCLE
   
!     ALTER ENTRY IN PLACE (NOTE THE CONVERSION TABLE ABOVE)
   
!     GET AREA FACTOR FROM ECPT AND REVISE ENTRY.
   
   300 IF (iarea == 0) GO TO 310
   area = recpt(iarea)
   IF (eltype == ctube) area = piovr4*(area**2-(area-2.*recpt(6))**2)
   GO TO 320
   310 area = 1.0
   320 i1 = INDEX + 3
   i2 = INDEX + 10
   DO  j = i1,i2
     z(j) = 0
   END DO
   i2 = i1 + points - 1
   igrid = grid1
   DO  j = i1,i2
     z(j)  = ecpt(igrid)
     igrid = igrid+1
   END DO
   rz(INDEX+11) = rz(INDEX+1)*area
   z(INDEX + 1) = points
   z(INDEX +12) = itype
 END DO
 GO TO 270
 
!     HBDY ELEMENTS OF EST FILE.  DO QBDY1, QBDY2, AND QVECT REFERENCES.
 
 
!     BUF(3) IS SET TO 0 AS A FLAG TO TELL IF HBDY HAS BEEN CALLED FOR
!     THIS ELEMENT.
 
 360 IF (TYPE(14,4)+TYPE(15,4)+TYPE(16,4) == 0) GO TO 130
 idbdy1 = TYPE(14,4) + 3
 idbdy2 = TYPE(15,4) + 2
 idqvec = TYPE(16,4) + 6
 qbdy1s = z(idbdy1-2)
 qbdy2s = z(idbdy2-1)
 qvects = z(idqvec-5)
 estwds = 53
 
!     READ AN HBDY ELEMENT ECPT FROM THE EST.
 
 370 CALL READ (*700,*120,est,ecpt,estwds,noeor,iflag)
 
!     LOOK FOR ID AMONG QBDY1 DATA
 
 IF (TYPE(14,4) == 0) GO TO 410
 iwords = 10
 j1 = idbdy1
 j2 = j1 + qbdy1s*iwords
 idptr = 2
 ASSIGN 380 TO iretrn
 CALL bisloc (*410,ecpt(1),z(idbdy1),10,qbdy1s,jpoint)
 
!     MATCH FOUND.  CHECK FOR PREVIOUS REFERENCE.
 
 INDEX = idbdy1 + jpoint - 3
 GO TO 710
 380 DO  INDEX = k1,k2,iwords
   IF (z(INDEX+10) == 1) GO TO 390
   WRITE  (outpt,382) ufm,ecpt(1)
   382 FORMAT (a23,' 2362, CHBDY CARDS WITH DUPLICATE IDS FOUND IN EST,',  &
       ' CHBDY ID NUMBER =',i9)
   nogo = .true.
   GO TO 650
   
!     ALTER ENTRY FOR OUTPUT.  GET AREA FACTORS FOR HBDY ELEMENT.
   
   390 CALL hbdy (ecpt,ecpt,2,rbuf,buf)
   z(INDEX +3) = buf(3)
   z(INDEX +4) = buf(4)
   z(INDEX +5) = buf(5)
   z(INDEX +6) = buf(6)
   rz(INDEX+7) = rbuf(7)*rz(INDEX+1)
   rz(INDEX+8) = rbuf(8)*rz(INDEX+1)
   rz(INDEX+9) = rbuf(9)*rz(INDEX+1)
   rz(INDEX+10)= rbuf(10)*rz(INDEX+1)
   z(INDEX +1) = ecpt(2)
 END DO
 
!     LOOK FOR ID AMONG QBDY2 DATA.
 
 410 IF (TYPE(15,4) == 0) GO TO 450
 iwords = 10
 j1 = idbdy2
 j2 = j1 + qbdy2s*iwords
 idptr = 1
 ASSIGN 420 TO iretrn
 CALL bisloc (*450,ecpt(1),z(idbdy2),10,qbdy2s,jpoint)
 
!     MATCH FOUND.  CHECK FOR PREVIOUS REFERENCE.
 
 INDEX = idbdy2 + jpoint - 2
 GO TO 710
 420 DO  INDEX = k1,k2,iwords
   IF (z(INDEX+10) == 1) GO TO 430
   WRITE (outpt,382) ufm,ecpt(1)
   nogo = .true.
   GO TO 650
   
!     ALTER ENTRY FOR OUTPUT.  GET AREA FACTORS FOR HBDY ELEMENT.
   
   430 CALL hbdy (ecpt,ecpt,2,rbuf,buf)
   rz(INDEX+7) = rbuf(7)*rz(INDEX+2)
   rz(INDEX+8) = rbuf(8)*rz(INDEX+3)
   rz(INDEX+9) = rbuf(9)*rz(INDEX+4)
   rz(INDEX+10)= rbuf(10)*rz(INDEX+5)
   z(INDEX +3) = buf(3)
   z(INDEX +4) = buf(4)
   z(INDEX +5) = buf(5)
   z(INDEX +6) = buf(6)
   z(INDEX +2) = ecpt(2)
 END DO
 
!     LOOK FOR ID AMONG QVECT DATA
 
 450 IF (TYPE(16,4) == 0) GO TO 490
 iwords = 19
 j1 = idqvec
 j2 = j1 + qvects*iwords
 idptr = 5
 ASSIGN 460 TO iretrn
 CALL bisloc (*490,ecpt(1),z(idqvec),19,qvects,jpoint)
 
!     MATCH FOUND.  CHECK FOR PREVIOUS REFERENCE.
 
 INDEX = idqvec + jpoint - 6
 GO TO 710
 460 DO  INDEX = k1,k2,iwords
   IF (z(INDEX+19) == 1) GO TO 470
   WRITE (outpt,382) ufm,ecpt(1)
   nogo = .true.
   GO TO 650
   
!     ALTER ENTRY FOR OUTPUT.  GET AREA FACTORS FOR HBDY ELEMENT.
   
   470 CALL hbdy (ecpt,ecpt,3,rbuf,buf)
   rz(INDEX+11) = rz(INDEX+2)
   rz(INDEX+12) = rz(INDEX+3)
   rz(INDEX+13) = rz(INDEX+4)
   rz(INDEX+14) = rbuf(11)
   rz(INDEX+15) = rbuf(12)
   rz(INDEX+16) = rbuf(13)
   rz(INDEX+17) = rbuf(14)
   rz(INDEX+18) = rbuf(15)
   rz(INDEX+19) = rbuf(16)
   q0           = rz(INDEX+1)
   z(INDEX + 1) = buf(3)
   z(INDEX + 2) = buf(4)
   z(INDEX + 3) = buf(5)
   z(INDEX + 4) = buf(6)
   z(INDEX + 6) = ecpt(2)
   rz(INDEX+ 7) = rbuf(7)*q0
   rz(INDEX+ 8) = rbuf(8)*q0
   rz(INDEX+ 9) = rbuf(9)*q0
   rz(INDEX+10) = rbuf(10)*q0
 END DO
 490 GO TO 370
 
!     EST HAS BEEN PASSED FOR ALL ELEMENTS.  NOW OUTPUT DATA TO NEWSLT.
 
 500 CALL CLOSE (est,clsrew)
 DO  j = 13,16
   jcore = TYPE(j,4)
   IF (jcore > 0) THEN
     GO TO   510
   ELSE
     GO TO   590
   END IF
   510 nwords = z(jcore+1)*TYPE(j,2) + 2
   
!     INSURE THAT ALL ENTRYS WERE MODIFIED.
!     CHECK WORD 7 FOR NO INTEGER 1 IN TYPES 14,15, AND 16.
!     CHECK WORD 11 FOR NO INTEGER 1 IN TYPE 13.
   
   k = 8
   IF (j == 13) k = 12
   i1 = jcore + k
   i2 = i1 + nwords - 3
   outcnt = TYPE(j,2)
   DO  l = i1,i2,outcnt
     IF (z(l) /= 1) CYCLE
     k = j - 12
     SELECT CASE ( k )
       CASE (    1)
         GO TO 520
       CASE (    2)
         GO TO 530
       CASE (    3)
         GO TO 540
       CASE (    4)
         GO TO 550
     END SELECT
     520 id = l - 9
     GO TO 560
     530 id = l - 5
     GO TO 560
     540 id = l - 6
     GO TO 560
     550 id = l - 2
     GO TO 560
     560 WRITE  (outpt,570) ufm,z(id)
     570 FORMAT (a23,' 3096, ELEMENT ID =',i9,' AS REFERENCED ON A QVOL, ',  &
         'QBDY1, QBDY2, OR QVECT LOAD CARD,', /5X,'COULD NOT BE ',  &
         'FOUND AMONG ACCEPTABLE ELEMENTS FOR THAT LOAD TYPE.')
     nogo = .true.
   END DO
   CALL WRITE (newslt,z(jcore),nwords,noeor)
 END DO
 
!     COMPLETE THIS LOAD SET RECORD ON -NEWSLT-.
 
 600 CALL WRITE (newslt,0,0,eor)
 i = i+1
 IF (i <= nrecs) GO TO 21
 
!     COPY BALANCE OF DATA ON -SLT- TO -NEWSLT- WHATEVER IT BE.
 
 620 CALL READ (*640,*630,slt,z,core,noeor,iwords)
 CALL WRITE (newslt,z,core,noeor)
 GO TO 620
 630 CALL WRITE (newslt,z,iwords,eor)
 GO TO 620
 
!     NEWSLT IS COMPLETE.
 
 640 CALL CLOSE (slt,clsrew)
 CALL CLOSE (newslt,clsrew)
 650 IF (nogo) CALL mesage (-61,0,subr)
 RETURN
 
!     FATAL FILE ERRORS
 
 660 CALL mesage (-1,newslt,subr)
 670 CALL mesage (-2,slt   ,subr)
 690 CALL mesage (-3,est   ,subr)
 700 CALL mesage (-2,est   ,subr)
 701 CALL mesage (-2,bgpdt ,subr)
 RETURN
 
!     INTERNAL ROUTINE TO FIND THE START AND END OF ENTRYS HAVING THE
!     SAME ID IN A GIVEN CARD-TYPE SET.
 
 
!     BACK UP TO FIRST ENTRY OF THIS ID.
 
 710 jindex = INDEX + idptr - iwords
 720 IF (jindex < j1) GO TO 730
 IF (z(jindex) /= ecpt(1)) GO TO 730
 jindex = jindex - iwords
 GO TO 720
 730 k1 = jindex + iwords - idptr
 
!     FIND LAST ENTRY OF THIS ID.
 
 jindex = k1 + iwords + idptr
 740 IF (jindex >= j2) GO TO 750
 IF (z(jindex) /= ecpt(1)) GO TO 750
 jindex = jindex + iwords
 GO TO 740
 750 k2 = jindex - iwords - idptr
 GO TO iretrn, (280,380,420,460)
 
!     SPECIAL PROCESSING FOR SPCFLD,CEMLOOP,MDIPOLE, AND GEMLOOP. SET UP
!     A VECTOR FOR ALL BUT SPCFLD CARDS, COMPUTE FIELD AT EACH POINT
!     IN BGPDT. WHEN FINISHED, ALL THE E AND M CARD TYPES WILL BE
!     ACCUMULATED INTO ONE SPCFLD-LIKE CARD WITH FIELD VALUSS AT EACH
!     POINT
 
 800 IF (ifirst == 1) GO TO 811
 ifirst = 1
 jcore1 = ncore+1
 jcoren = ncore+3*nrowsp
 IF (jcoren > core) CALL mesage (-8,0,subr)
 
 DO  j1 = jcore1,jcoren
   rz(j1) = 0.
 END DO
 
!     ALL E AND M CARDS WILL BE COMBINED INTO ONE LOGICAL CARD OF
!     TYPE=20, 3*NROWSP VALUES HCX,HCY,HCZ AT EACH POINT IN THE MODEL.
!     FOR CEMLOOP AND GEMLOOP, WE MUST PICK UP BGPDT FOREACH POINT AND
!     COMPUTE FIELD FOR EACH LOOP
! *** 10/1/80 WE MUST ALSO FIND HC AT INTEGRATION POINTS AND CENTROIDS.
!     SO ALSO COPY SLT INFO TO NEWSLT FOR USE IN EANDM
 
 
!     1ST OCCURRENCE OF A CARD TYPE. CHECK ON TYPE
 
 811 jtype = itype-19
 SELECT CASE ( jtype )
   CASE (    1)
     GO TO 812
   CASE (    2)
     GO TO 840
   CASE (    3)
     GO TO 840
   CASE (    4)
     GO TO 840
 END SELECT
 
!     SPCFLD
 
 812 buf(1) = 20
 buf(2) = 1
 CALL WRITE (newslt,buf,2,0)
 DO  j1 = 1,entrys
   
!     READ ONE SPCFLD CARD
   
   CALL fread (slt,buf,5,0)
   IF (buf(5) /= -1) GO TO 825
   
!     BUF(1)=CID WHICH IS ASSUMED TO BE 0 FOR NOW
   
!     ALL GRIDS GET HC
   
   DO  j2 = jcore1,jcoren,3
     rz(j2  ) = rz(j2  )+rbuf(2)
     rz(j2+1) = rz(j2+1)+rbuf(3)
     rz(j2+2) = rz(j2+2)+rbuf(4)
   END DO
   CYCLE
   
   825 isub = ncore+3*buf(5)-2
   rz(isub  ) = rz(isub  )+rbuf(2)
   rz(isub+1) = rz(isub+1)+rbuf(3)
   rz(isub+2) = rz(isub+2)+rbuf(4)
 END DO
 
!     DONE WITH ALL SPCFLD CARDS IN THIS LOAD SET. CHECK FOR OTHER CARD
!     TYPES IN THIS LOAD SET
 
 CALL WRITE (newslt,rz(jcore1),3*nrowsp,0)
 GO TO 910
 
!     CEMLOOP,GEMLOOP, OR MDIPOLE
!     CHECK FOR ENOUGH CORE TO READ IN BGPDT. IF NOT, READ ONE POINT AT
!     A TIME
 
 840 IF (bgopen) GO TO 850
 
!     IF MODCOM(9) IS NOT SET TO NONZERO, THEN WE WILL NOT COMPUTE HCFLD
!     AT GRID POINTS FOR COILS, ETC.(ONLY SPCFLD) SINCE IT TAKES TIME
!     AND IS NOT NEEDED IN ANY SUBSEQUENT COMPUTATION. (ONLY SPCFLD INFO
!     IS NEEDED LATER. ALL OTHER HC INFO IS COMPUTED LATER) IF MODCOM(9)
!     IS SET TO NONZERO, HCFLD IS COMPUTED AT THE POINTS FOR ALL LOAD
!     TYPES AND CAN BE PRINTED FOR INFORMATIONAL PURPOSES IF DESIRED.
 
 IF (ksystm(65) == 0) GO TO 850
 CALL gopen (bgpdt,z(buf3),0)
 mcb(1) = bgpdt
 CALL rdtrl (mcb)
 npts   = mcb(2)
 bgcore = .true.
 bgopen = .true.
 IF (jcoren+4*npts > core) bgcore=.false.
 next = jcoren+4*npts
 IF (.NOT.bgcore) next=jcoren
 IF (bgcore) CALL fread (bgpdt,z(jcoren+1),4*npts,0)
 850 CONTINUE
 CALL WRITE (newslt,buf,2,0)
 
 DO  j1 = 1,entrys
   
!     READ CEMLOOP, GEMLOOP, OR MDIPOLE ENTRY
   
   iwords = 12
   IF (itype == 22) iwords = 48
   IF (itype == 23) iwords = 9
   CALL fread (slt,buf,iwords,0)
   CALL WRITE (newslt,buf,iwords,0)
   IF (ksystm(65) == 0) CYCLE
   
   
!     DO THIS LOOP FOR ALL POINTS
   
   DO  kk = 1,npts
     IF (bgcore) GO TO 880
     CALL fread (bgpdt,buf,4,0)
     IF (buf(1) == -1) CYCLE
     xx = rbuf(2)
     yy = rbuf(3)
     zz = rbuf(4)
     GO TO 885
     880 jcor = jcoren+4*kk
     IF (z(jcor-3) == -1) CYCLE
     xx = rz(jcor-2)
     yy = rz(jcor-1)
     zz = rz(jcor )
     885 IF (itype == 21) GO TO 886
     IF (itype == 23) GO TO 888
     CALL geloop (rbuf,buf,xx,yy,zz,hc1,hc2,hc3)
     GO TO 887
     886 CALL axloop (rbuf,buf,xx,yy,zz,hc1,hc2,hc3)
     GO TO 887
     888 CALL dipole (rbuf,buf,xx,yy,zz,hc1,hc2,hc3)
     887 isub = ncore+3*kk-2
     rz(isub  ) = rz(isub  )+hc1
     rz(isub+1) = rz(isub+1)+hc2
     rz(isub+2) = rz(isub+2)+hc3
     
!     GO BACK FOR ANOTHER POINT
     
   END DO
   IF (bgcore) CYCLE
   CALL REWIND (bgpdt)
   CALL fwdrec (*701,bgpdt)
   
!     GET ANOTHER LOOP OR DIPOLE
   
 END DO
 
!     CHECK IF NEXT CARD TYPE IS 21 ,22, OR 23. CARD TYPES ON SLT ARE
!     IN INCREASING CARD TYPE). IF SO, STAY HERE. OTHERWISE, WRITE OUT
!     ALL CARD TYPES GENERATING AN SPCFLD-TYPE CARD AND GOING ONTO
!     HCFLDS MUST HAVE CONSECUTIVE TYPE NUMBERS FOR THIS SPECIAL
!     PROCESSING THE GENERATED SPCFLD AND GO BACK TO NORMAL PROCESSING
 
 910 CALL READ (*670,*920,slt,buf,2,noeor,iwords)
 itype  = buf(1)
 entrys = buf(2)
 IF (buf(1) >= 20 .AND. buf(1) <= 23) GO TO 811
 IEOR = 0
 GO TO 930
 920 IEOR = 1
 930 buf(1) =-20
 buf(2) = 1
 CALL WRITE (newslt,buf,2,0)
 CALL WRITE (newslt,rz(jcore1),3*nrowsp,0)
 IF (bgopen) CALL CLOSE (bgpdt,1)
 IF (IEOR == 1) GO TO 110
 GO TO 35
 
!     REMFLUX PROCESSING. CREATE A VECTOR OF ORDER 3N,N=NUMBER OF
!     ELEMENTS IN MODEL,N IS 1ST TRAILER WORD OF EST. THE VECTOR
!     CONTAINS TOTAL BX,BY,BZ FROM ALL REMFLUX CARDS FOR EACH ELEMENT
!     IN THE ORDER OF ELEMENTS ON EST
 
 1000 CALL gopen (est,z(buf3),0)
 mcb(1) = est
 CALL rdtrl (mcb)
 nel    = mcb(2)
 jcore1 = ncore+1
 jcoren = ncore+3*nel
 jcorex = jcoren+5*entrys
 IF (jcorex > core) CALL mesage (-8,0,subr)
 
 nels = 0
 DO  j1 = jcore1,jcoren
   rz(j1) = 0.
 END DO
 
!     READ ALL REMFLUX CARDS
 
 CALL fread (slt,rz(jcoren+1),5*entrys,0)
 
 1020 CALL READ (*1050,*690,est,eltype,1,0,iwords)
 idx = (eltype-1)*incr
 estwds = NE(idx+12)
 1025 CALL READ (*700,*1020,est,elid,1,0,iwords)
 nels = nels+1
 isub = ncore+3*nels-2
 CALL fread (est,dum,-estwds+1,0)
 
!     CHECK FOR THIS ELID AMONG THE REMFLUX CARDS
 
 DO  j1 = 1,entrys
   isub1 = jcoren+5*j1
   IF (z(isub1) ==   -1) GO TO 1030
   IF (elid /= z(isub1)) CYCLE
   
!     MATCH-STORE THIS PERM MAG
   
   1030 rz(isub  ) = rz(isub  )+rz(isub1-3)
   rz(isub+1) = rz(isub+1)+rz(isub1-2)
   rz(isub+2) = rz(isub+2)+rz(isub1-1)
 END DO
 
!     READ ANOTHER ELEMENT ID
 
 GO TO 1025
 
!     EST EXHAUSTED
 
 1050 CALL CLOSE (est,1)
 buf(1) = 24
 buf(2) = 1
 CALL WRITE (newslt,buf,2,0)
 CALL WRITE (newslt,rz(jcore1),3*nel,0)
 GO TO 30
END SUBROUTINE ssgslt
