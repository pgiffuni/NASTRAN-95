SUBROUTINE ifp5
     
    !     ACOUSTIC CAVITY PREFACE ROUTINE
 
    !     THIS PREFACE MODULE OPERATES ON ACOUSTIC-CAVITY-ANALYSIS DATA
    !     CARDS WHICH AT THIS POINT ARE IN THE FORM OF IFP-OUTPUT IMAGES ON
    !     THE AXIC DATA BLOCK.
 
    !     THE FOLLOWING LIST GIVES THE CARD IMAGES IFP5 WILL LOOK FOR ON THE
    !     AXIC OR GEOM2 DATA BLOCKS,  THE CARD IMAGES IFP5 WILL GENERATE OR
    !     MODIFY, AND THE DATA BLOCKS ONTO WHICH THE GENERATED OR MODIFIED
    !     CARD IMAGES WILL BE PLACED.
 
    !      IFP5 INPUT         IFP5 OUTPUT        DATA BLOCK
    !      CARD IMAGE         CARD IMAGE         EFFECTED
    !     ------------       -----------        ----------
    !      AXSLOT/AXIC        -NONE-             -NONE-
    !      CAXIF2/GEOM2       PLOTEL             GEOM2
    !      CAXIF3/GEOM2       PLOTEL             GEOM2
    !      CSLOT3/GEOM2       PLOTEL             GEOM2
    !      CSLOT4/GEOM2       PLOTEL             GEOM2
    !      CAXIF4/GEOM2       PLOTEL             GEOM2
    !      GRIDF/AXIC         GRID               GEOM1
    !      GRIDS/AXIC         GRID               GEOM1
    !      SLBDY/AXIC         CELAS2             GEOM2
 
    !     SOME OF THE ABOVE OUTPUT DATA CARDS ARE A FUNCTION OF MORE THAN
    !     ONE INPUT DATA CARDS
 
    LOGICAL :: g1eof    ,g2eof    ,plotel   ,any
    INTEGER :: sysbuf   ,output   ,rd       ,rdrew    ,cls      ,  &
        buf(24)  ,z        ,wrt      ,wrtrew   ,clsrew   ,  &
        core     ,subr(2)  ,flag     ,  &
        caxif(6) ,cslot(4) ,grid(2)  ,celas2(2),plotls(2),  &
        card(10) ,eor      ,grids(2) ,gridf(2) ,slbdy(2) ,  &
        words    ,buf1     ,buf2     ,buf3     ,buf4     ,  &
        FILE     ,geom1    ,geom2    ,scrt1    ,scrt2    ,  &
        axslot(2),msg1(2)  ,msg2(2)  ,axic     ,entrys
    REAL :: rz(4)    ,rbuf(24) ,rr(3)    ,zz(3)    ,ww(3)    ,  &
        l1       ,l2       ,l3       ,l1l2     ,kf       ,  &
        LE       ,lc       ,rcard(10)

    COMMON /condas/ consts(5)
    COMMON /system/ ksystm(65)
    COMMON /names / rd,rdrew ,wrt      ,wrtrew   ,clsrew   ,cls
    COMMON /zzzzzz/ z(1)

    EQUIVALENCE (consts(2),twopi)  ,(ksystm(1),sysbuf),  &
                (ksystm(2),output) ,(z(1),rz(1)), (buf(1),rbuf(1))   ,(card(1),rcard(1))
    DATA    axslot/ 1115,   11 /
    DATA    slbdy / 1415,   14 /
    DATA    caxif / 2108,   21 , 2208,   22,  &
                    2308,   23 /
    DATA    celas2/  701,    7 /
    DATA    cslot / 4408,   44 , 4508,   45 /
    DATA    grid  / 4501,   45 /
    DATA    grids / 1315,   13 /
    DATA    gridf / 1215,   12 /
    DATA    plotls/ 5201,   52 /
    DATA    subr  / 4HIFP5, 4H    /
    DATA    eor   , noeor / 1, 0  /
 
    !     NOTE...  SCRATCH2 IN IFP5 AS IN IFP4 IS EQUIVALENCED TO THE
    !     -FORCE- DATA BLOCK.
 
    DATA axic, geom1, geom2, scrt1, scrt2 / 215, 201, 208, 301, 213 /
    DATA msg1/ 4HIFP5, 4HBEGN/, msg2 /4HIFP5, 4HEND /

    !     DEFINE CORE AND BUFFER POINTERS

    CALL conmsg (msg1,2,0)
    core = korsz(z)
    buf1 = core - sysbuf - 2
    buf2 = buf1 - sysbuf - 2
    buf3 = buf2 - sysbuf - 2
    buf4 = buf3 - sysbuf - 2
    core = buf4 - 1
    icrq = 100  - core
    IF (core < 100) GO TO 980
    plotel = .false.

    !     OPEN AXIC DATA BLOCK. (IF NAME IS NOT IN FIST RETURN - NO MESSAGE)

    CALL preloc (*910,z(buf1),axic)

    !     PICK UP THE AXSLOT CARD AND SAVE THE VALUES ON IT.
    !     RHOD, BD, N, WD, MD (FATAL ERROR IF NOT PRESSENT)

    CALL locate (*10,z(buf1),axslot,flag)
    CALL READ (*920,*30,axic,z(1),6,eor,words)
10  CALL ifp5a (1)
    WRITE  (output,20)
20  FORMAT (' AXSLOT DATA CARD IS NOT PRESENT OR IS INCORRECT.')

    !     SET VALUES FOR CONTINUING DATA CHECK

    rhod = 0.0
    bd   = 0.0
    n    = 0
    wd   = 1.0
    md   = 0
    GO TO 40
30  IF (words /= 5) GO TO 10
    rhod = rz(1)
    bd   = rz(2)
    j    = 3
    n    = z(j)
    wd   = rz(4)
    md   = z(j+2)

    !     READ GRIDS DATA CARDS INTO CORE FROM AXIC DATA BLOCK.

40  igrids = 1
    ngrids = igrids - 1
    CALL locate (*70,z(buf1),grids,flag)
    CALL READ (*920,*60,axic,z(igrids),core,eor,words)
    CALL ifp5a (2)
    WRITE  (output,50)
50  FORMAT (49H insufficient core TO hold all grids card images.)
    WRITE  (output,431) core
    GO TO 70
60  ngrids = ngrids + words

    !     READ GRIDF DATA CARDS INTO CORE FROM AXIC DATA BLOCK.

70  igridf = ngrids + 1
    ngridf = igridf - 1
    CALL locate (*100,z(buf1),gridf,flag)
    CALL READ (*920,*90,axic,z(igridf),core-ngrids,eor,words)
    CALL ifp5a (3)
    WRITE  (output,80)
80  FORMAT (49H insufficient core TO hold all gridf card images.)
    icrq = core - ngrids
    WRITE (output,431) icrq
    GO TO 100
90  ngridf = ngridf + words

    !     INSERT DEFAULT SLOT WIDTH INTO ANY GRIDS IMAGE HAVING NONE
    !     SPECIFIED BY THE USER.

100 IF (ngrids < igrids) GO TO 170
    DO  i = igrids,ngrids,5
        IF (z(i+3) == 1) rz(i+3) = wd
    END DO

    !     CREATE A GRIDF CARD FOR EACH GRIDS DATA CARD THAT HAS A NON-ZERO
    !     IDF

    DO  i = igrids,ngrids,5
        IF (z(i+4) > 0.0) THEN
            GO TO   130
        ELSE
            GO TO   140
        END IF
130     ngridf = ngridf + 3
        IF (ngridf > core) GO TO 150
        z(ngridf-2) = z(i+4)
        z(ngridf-1) = z(i+1)
        z(ngridf  ) = z(i+2)
140 CONTINUE
    END DO
    GO TO 170
150 CALL ifp5a (4)
    WRITE  (output,160)
160 FORMAT (' INSUFFICIENT CORE TO HOLD ALL GRIDF CARD IMAGES BEING ',  &
        'CREATED INTERNALLY DUE TO GRIDS CARDS SPECIFYING AN IDF.')
    icrq = ngridf - core
    WRITE (output,431) icrq
    ngridf = ngridf - 3

    !     SORT THE GRIDF CARDS ON THEIR ID.

170 IF (ngridf > igridf) CALL sort (0,0,3,1,z(igridf),ngridf-igridf+1)

    !     OPEN GEOM1 AND SCRATCH1, COPY HEADER REC FROM GEOM1 TO SCRATCH1.

    CALL ifp4c (geom1,scrt1,z(buf2),z(buf3),g1eof)

    !     COPY ALL DATA FROM GEOM1 TO SCRATCH1 UP TO FIRST GRID CARD.

    CALL ifp4b (geom1,scrt1,any,z(ngridf+1),core-ngridf,grid,g1eof)
    FILE = geom1

    !     CREATE GRID CARDS FROM GRIDS AND GRIDF CARDS.
    !     MERGE THESE INTO EXISTING GRID CARDS CHECKING FOR DUPLICATE ID-S.

    igf  = igridf
    idgf = 0
    igs  = igrids
    idgs = 0
    IF (igf < ngridf) idgf = z(igf)
    IF (igs < ngrids) idgs = z(igs)
    card(2) = 0
    card(6) =-1
    card(7) = 0
    card(8) = 0

    !     READ A GRID CARD INTO BUF.

    IF (.NOT.any) GO TO 190
180 CALL READ (*940,*190,geom1,buf,8,noeor,words)
    idg = buf(1)
    GO TO 200
190 idg = 0

    !     DETERMINE WHETHER GRID, GRIDF, OR GRIDS CARD IS TO OUTPUT NEXT.

200 IF (  idg     > 0) THEN
        GO TO   250
    END IF
210 IF (  idgf    > 0) THEN
        GO TO   230
    END IF
220 IF (  idgs    > 0) THEN
        GO TO   370
    ELSE
        GO TO   390
    END IF
230 IF (  idgs    > 0) THEN
        GO TO   240
    ELSE
        GO TO   360
    END IF
240 IF (idgf-idgs < 0) THEN
        GO TO   360
    ELSE IF (idgf-idgs == 0) THEN
        GO TO   330
    ELSE
        GO TO   370
    END IF
250 IF (  idgf    > 0) THEN
        GO TO   280
    END IF
260 IF (  idgs    > 0) THEN
        GO TO   270
    ELSE
        GO TO   350
    END IF
270 IF (idg -idgs < 0) THEN
        GO TO   350
    ELSE IF (idg -idgs == 0) THEN
        GO TO   330
    ELSE
        GO TO   370
    END IF
280 IF (idg -idgf < 0) THEN
        GO TO   310
    ELSE IF (idg -idgf == 0) THEN
        GO TO   330
    END IF
290 IF (  idgs    > 0) THEN
        GO TO   300
    ELSE
        GO TO   360
    END IF
300 IF (idgf-idgs < 0) THEN
        GO TO   360
    ELSE IF (idgf-idgs == 0) THEN
        GO TO   330
    ELSE
        GO TO   370
    END IF
310 IF (  idgs    > 0) THEN
        GO TO   320
    ELSE
        GO TO   360
    END IF
320 IF (idg -idgs < 0) THEN
        GO TO   350
    ELSE IF (idg -idgs == 0) THEN
        GO TO   330
    ELSE
        GO TO   370
    END IF

    !     ERROR - DUPLICATE ID-S ENCOUNTERED

330 CALL ifp5a (10)
    WRITE  (output,340) idg,idgs,idgf
340 FORMAT (' ONE OF THE FOLLOWING NON-ZERO IDENTIFICATION NUMBERS ',  &
        'APPEARS ON SOME COMBINATION', /,' OF GRID, GRIDS, OR ',  &
        'GRIDF BULK DATA CARDS.',3(6H   id=,i12))
    IF (idg == idgf) GO TO 350
    IF (idg == idgs) GO TO 350
    GO TO 370

    !     OUTPUT GRID CARD AND READ ANOTHER

350 CALL WRITE (scrt1,buf,8,noeor)
    GO TO 180

    !     OUTPUT A GRID FROM GRIDF CARD.

360 card(1)  = idgf
    rcard(3) = rz(igf+1)
    rcard(4) = rz(igf+2)
    rcard(5) = 0.0
    igf = igf + 3
    IF (igf > ngridf) idgf = 0
    IF (idgf /= 0) idgf = z(igf)
    GO TO 380

    !     OUTPUT A GRID FROM GRIDS CARD.

370 card(1)  = idgs
    rcard(3) = rz(igs+1)
    rcard(4) = rz(igs+2)
    rcard(5) = rz(igs+3)
    igs = igs + 5
    IF (igs > ngrids) idgs = 0
    IF (idgs /= 0) idgs = z(igs)
380 CALL WRITE (scrt1,card,8,noeor)
    GO TO 200

    !     ALL GRID CARDS HAVE BEEN OUPTUT, WRITE EOR.

390 CALL WRITE (scrt1,0,0,eor)

    !     COPY BALANCE OF GEOM1 TO SCRT1, WRAP UP AND COPY BACK.

    CALL ifp4b (geom1,scrt1,any,z(igridf),core-igridf,-1,g1eof)

    !     SLBDY CARD IMAGES ARE NOW PROCESSED AND A BOUNDARY TABLE IS FORMED
    !     IN CORE.  EACH ENTRY IN THE TABLE CONTAINS,

    !              IDS , IDS   , IDS   , RHO, M
    !                 I     I-1     I+1

    !     IDS    = -1 IF IDS  IS THE FIRST ID ON SLBDY CARD.
    !        I-1            I

    !     IDS    = -1 IF IDS  IS THE LAST ID ON SLBDY CARD.
    !        I+1            I

    islbdy = ngrids + 1
    nslbdy = islbdy - 1
    CALL locate (*440,z(buf1),slbdy,flag)
400 CALL READ (*920,*440,axic,buf,2,noeor,words)
    rho = rbuf(1)
    m   = buf(2)
    idsl1 = -1
    CALL READ (*920,*930,axic,ids,1,noeor,words)
410 CALL READ (*920,*930,axic,idsp1,1,noeor,words)

    !     PLACE 5 WORD ENTRY INTO CORE

    nslbdy = nslbdy + 5
    IF (nslbdy > core) GO TO 420
    z(nslbdy-4) = ids
    z(nslbdy-3) = idsl1
    z(nslbdy-2) = idsp1
    rz(nslbdy-1) = rho
    z(nslbdy   ) = m
    idsl1 = ids
    ids   = idsp1
    IF (idsp1+1 == 0) THEN
        GO TO   400
    ELSE
        GO TO   410
    END IF

    !     OUT OF CORE

420 CALL ifp5a (5)
    WRITE  (output,430)
430 FORMAT (' INSUFFICIENT CORE TO CONSTRUCT ENTIRE BOUNDARY TABLE ',  &
        'FOR SLBDY CARDS PRESENT.')
    icrq = nslbdy - core
    WRITE  (output,431) icrq
431 FORMAT (5X,24HADDITIONAL core needed =,i8,7H words.)
    nslbdy = nslbdy - 5

    !     SKIP BALANCE OF SLBDY DATA.

    CALL READ (*920,*440,axic,buf,1,eor,words)

    !     SORT BOUNDARY TABLE ON IDS . (FIRST WORD OF EACH ENTRY)
    !                               I

440 IF (nslbdy > islbdy) CALL sort (0,0,5,1,z(islbdy),nslbdy-islbdy+1)
    !/////
    !     CALL BUG (10H BOUNDRY      ,440,Z(ISLBDY),NSLBDY-ISLBDY+1)

    !     OPEN GEOM2, OPEN SCRATCH2, COPY HEADER REC FROM GEOM2 TO SCRATCH2.

    FILE = geom2
    CALL ifp4c (geom2,scrt2,z(buf2),z(buf3),g2eof)

    !     OPEN SCRATCH1, FOR TEMPORARY OUTPUT OF PLOTEL IMAGES CREATED FROM
    !     CAXIF2, CAXIF3, CAXIF4, CSLOT3, AND CSLOT4 CARDS.

    FILE = scrt1
    CALL OPEN (*960,scrt1,z(buf4),wrtrew)

    !     CREATE PLOTEL IMAGES FROM CAXIF2, CAXIF3, AND CAXIF4 AT THIS TIME

    FILE = geom2
    DO  i = 1,3
        ibase = (i-1)*1000000
        IF (i == 3) ibase = 4000000
        k  = 2*i - 1
        k4 = i + 5
  
        !     CHECK TRAILER TO SEE IF CAXIF(I+1) EXISTS
  
        CALL ifp4f (caxif(k+1),geom2,any)
        IF (.NOT.any) CYCLE
  
        !     COPY ALL DATA FROM GEOM2 TO SCRATCH2 UP TO FIRST CAXIF(I+1) IMAGE.
  
        CALL ifp4b(geom2,scrt2,any,z(nslbdy+1),core-nslbdy,caxif(k),g2eof)
        IF (.NOT.any) GO TO 1610
  
        !     COPY EACH IMAGE TO SCRATCH2 AND CREATE PLOTELS AT SAME TIME.
  
460     CALL READ (*940,*480,geom2,buf,k4,noeor,words)
        CALL WRITE (scrt2,buf,k4,noeor)
        nlines = i + 1
        IF (i == 1) nlines = 1
        DO  j = 1,nlines
            card(1) = buf(1) + ibase + j*1000000
            card(2) = buf(j+1)
            jj = j + 1
            IF (jj > nlines .AND. nlines /= 1) jj = 1
            card(3) = buf(jj+1)
            CALL WRITE (scrt1,card,3,noeor)
        END DO
        plotel = .true.
        GO TO 460
  
        !     END OF RECORD HIT ON GEOM2.  COMPLETE RECORD ON SCRATCH2
  
480     CALL WRITE (scrt2,0,0,eor)
    END DO

    !     COPY ALL DATA FROM GEOM2 TO SCRATCH2 UP TO FIRST CELAS2 CARD
    !     IMAGE.

    CALL ifp4b (geom2,scrt2,any,z(nslbdy+1),core-nslbdy,celas2,g2eof)

    !     COPY ANY CELAS2 DATA CARDS, MAKE SURE ALL ID ARE LESS THAN
    !     10000001.

    IF (.NOT.any) GO TO 540
510 CALL READ (*940,*540,geom2,buf,8,noeor,words)
    IF (buf(1) < 10000001) GO TO 530
    CALL ifp5a (6)
    WRITE  (output,520) buf(1)
520 FORMAT (' CELAS2 DATA CARD HAS ID =',i14,  &
        ', WHICH IS GREATER THAN 10000000,', /,' AND 10000000 IS THE ',  &
        'LIMIT FOR CELAS2 ID WITH ACOUSTIC ANALYSIS DATA CARDS PRESENT')
530 CALL WRITE (scrt2,buf,8,noeor)
    GO TO 510

    !     OUTPUT THREE CELAS2 IMAGES FOR EACH ENTRY IN THE BOUNDARY TABLE.

540 IF (nslbdy < islbdy) GO TO 800
    entrys = (ngrids-igrids+1)/5
    !/////
    !     CALL BUG(10H BOUNDRY      ,540,Z(ISLBDY),NSLBDY-ISLBDY+1)
    !     CALL BUG(10H GRIDS        ,540,Z(IGRIDS),NGRIDS-IGRIDS+1)
    ide = 10000000
    loop790:  DO  i = islbdy,nslbdy,5
  
        !     FIND  R, Z, W FOR IDS , IDS   , IDS    RESPECTIVELY.
        !                          I     I-1     I+1
  
        k   = 0
        is1 = i
        is3 = i + 2
        DO  j = is1,is3
            k = k + 1
            IF (z(j)   < 0.0) THEN
                GO TO   570
            ELSE IF (z(j)   == 0.0) THEN
                GO TO   580
            END IF
550         IF (entrys > 0.0) THEN
                GO TO   560
            ELSE
                GO TO   580
            END IF
560         kid = z(j)
            CALL bisloc (*580,kid,z(igrids),5,entrys,jpoint)
            ntemp = igrids + jpoint
    
            !     NTEMP NOW POINTS TO THE SECOND WORD OF THE GRIDS ENTRY HAVING
            !     THE ID SPECIFIED BY Z(J).  (1ST,2ND,OR 3RD ID IN SLBDY ENTRY)
    
    
            !     NO CELAS2 CARDS ARE GENERATED IF GRIDS FOR IDS  HAS NO IDF.
            !                                                   I
    
            IF (k == 1) idf = z(ntemp+3)
            IF (k == 1 .AND. idf <= 0) CYCLE loop790
            rr(k) = rz(ntemp  )
            zz(k) = rz(ntemp+1)
            ww(k) = rz(ntemp+2)
            CYCLE
    
            !     IDS = -1
    
570         rr(k) = rr(1)
            zz(k) = zz(1)
            ww(k) = ww(1)
            CYCLE
    
            !     IDS COULD NOT BE FOUND IN GRIDS ENTRYS.
    
580         CALL ifp5a (7)
            WRITE  (output,590) z(j)
590         FORMAT (11H slbdy id =,i12,  &
                ' DOES NOT APPEAR ON ANY GRIDS DATA CARD.')
            rr(k) = 0.0
            zz(k) = 0.0
            ww(k) = 0.0
        END DO
  
        !     COMPUTE GEOMETRY AND OTHER DATA.
  
        l1 = SQRT((zz(3)-zz(1))**2 + (rr(3)-rr(1))**2)
        l2 = SQRT((zz(2)-zz(1))**2 + (rr(2)-rr(1))**2)
        l3 = SQRT((zz(3)-zz(2))**2 + (rr(3)-rr(2))**2)/2.0
  
        l1l2 = (l1 + l2)*4.0
        IF (l1l2 > 0) THEN
            GO TO   630
        END IF
  
        !     ERROR, ZERO OR NEGATIVE LENGTH
  
610     CALL ifp5a (8)
        WRITE  (output,620) z(i),z(i+1),z(i+2)
620     FORMAT (' ONE OR MORE OF THE FOLLOWING ID-S NOT EQUAL TO -1 HAVE',  &
            ' INCORRECT OR NO GEOMETRY DATA',/3(10X,4HID = ,i10))
        CYCLE loop790
  
        !     COMPUTE W-BAR AND R-BAR
  
630     wbar = (l1*ww(3) + l2*ww(2))/l1l2 + 0.75*ww(1)
        rbar = (l1*rr(3) + l2*rr(2))/l1l2 + 0.75*rr(1)
  
        IF (wbar   == 0.0) THEN
            GO TO   610
        END IF
640     IF (rbar   == 0.0) THEN
            GO TO   610
        END IF
650     IF (z(i+4) == 0.0) THEN
            GO TO   740
        END IF
  
        !     COMPUTE BETA,LC
  
660     beta = (twopi*rbar)/(FLOAT(z(i+4))*wbar)
        IF (beta - 1.0 > 0.0) THEN
            GO TO   670
        ELSE
            GO TO   610
        END IF
670     bl1 = beta - 1.0
        bp1 = beta + 1.0
        lc  = wbar/twopi
        lc  = lc*((beta+1.0/beta)*ALOG(bp1/bl1) + 2.0*ALOG(bp1*bl1/(4.0*beta)))
        term = 0.01*wbar
        LE  = AMAX1(lc,term)
        IF (LE) 680,610,680
680     IF (rz(i+3) == 0.0) THEN
            GO TO   690
        ELSE
            GO TO   710
        END IF
690     CALL ifp5a (9)
        WRITE  (output,700) z(i)
700     FORMAT (' RHO AS SPECIFIED ON SLBDY OR AXSLOT CARD IS 0.0 FOR ID',  &
            ' =',i12)
        CYCLE loop790
  
        !     FIND F  = M, IF N=0 OR N=M/2  OTHERWISE F  = M/2
        !           I                                  I
  
710     IF (n == 0 .OR. 2*n == z(i+4)) GO TO 720
        fi = FLOAT(z(i+4))/2.0
        GO TO 730
720     fi = z(i+4)
730     kf = (wbar*l3*fi)/(rz(i+3)*LE)
        GO TO 750
  
        !     M = 0, THUS K  = 0.0
        !                  F
  
740     kf = 0.0
  
        !                           N WBAR
        !                      SIN( ------ )
        !                           2 RBAR
        !     COMPUTE  ALPHA = --------------
        !                           N WBAR
        !                         ( ------ )
        !                           2 RBAR
  
750     term = (FLOAT(n)*wbar)/(2.0*rbar)
        IF (term == 0.0) THEN
            GO TO   770
        END IF
760     alpha = SIN(term)/term
        GO TO 780
770     alpha = 1.0
  
        !  OUTPUT THE 3 CELAS2 CARDS
  
780     buf( 1) = ide + 1
        rbuf(2) = kf*(1.0 - alpha)
        buf( 3) = z(i)
        buf( 4) = 0
        buf( 5) = 1
        buf( 6) = 0
        buf( 7) = 0
        buf( 8) = 0
        buf( 9) = ide + 2
        rbuf(10)= kf*alpha
        buf(11) = z(i)
        buf(12) = idf
        buf(13) = 1
        buf(14) = 1
        buf(15) = 0
        buf(16) = 0
        buf(17) = ide + 3
        rbuf(18)= kf*alpha*(alpha - 1.0)
        buf(19) = idf
        buf(20) = 0
        buf(21) = 1
        buf(22) = 0
        buf(23) = 0
        buf(24) = 0
        CALL WRITE (scrt2,buf,24,noeor)
        ide = ide + 3
    END DO loop790

    !     COMPLETE THE CELAS2 RECORD.

800 CALL WRITE (scrt2,0,0,eor)

    !     CREATE PLOTEL IMAGES FROM CSLOT3, AND CSLOT4 AT THIS TIME IF ANY

    DO  i = 1,2
        ibase = 3000000*i + 5000000
        k  = 2*i - 1
        k6 = i + 7
  
        !     CHECK TRAILER BIT TO SEE IF CSLOT(I+2) EXISTS.
  
        CALL ifp4f (cslot(k+1),geom2,any)
        IF (.NOT.any) CYCLE
  
        !     COPY ALL DATA FROM GEOM2 TO SCRATCH2 UP TO FIRST CSLOT(I+2) IMAGE.
  
        CALL ifp4b (geom2,scrt2,any,z(igrids),core-igrids,cslot(k),g2eof)
        IF (.NOT.any) GO TO 1610
  
        !     COPY EACH IMAGE TO SCRATCH2 AND CREATE PLOTELS AT SAME TIME.
  
810     CALL READ (*940,*830,geom2,buf,k6,noeor,words)
        CALL WRITE (scrt2,buf,k6,noeor)
        nlines = i + 2
        DO  j = 1,nlines
            card(1) = buf(1) + ibase + j*1000000
            card(2) = buf(j+1)
            jj = j + 1
            IF (jj > nlines) jj = 1
            card(3) = buf(jj+1)
            CALL WRITE (scrt1,card,3,noeor)
        END DO
        plotel = .true.
        GO TO 810
  
        !     END OF RECORD ON GEOM2.  COMPLETE RECORD ON SCRATCH2.
  
830     CALL WRITE (scrt2,0,0,eor)
    END DO

    !     APPEND PLOTELS ON SCRATCH1 TO ANY PLOTELS ON GEOM2.
    !     MAKE SURE ALL PLOTEL ID-S ARE .LE. 1000000
    !     /// ID CHECK NOT IN YET.
    !     POSITION TO PLOTELS ON GEOM2 IF ANY ARE ON SCRATCH1

    IF (.NOT. plotel) GO TO 900
    CALL ifp4b (geom2,scrt2,any,z(igrids),core-igrids,plotls,g2eof)
    IF (.NOT.any) GO TO 870

    !     BLAST COPY PLOTELS FROM GEOM2 TO SCRATCH2

850 CALL READ (*940,*860,geom2,z(igrids),core-igrids,noeor,words)
    CALL WRITE (scrt2,z(igrids),core-igrids,noeor)
    GO TO 850
860 CALL WRITE (scrt2,z(igrids),words,noeor)

    !     CLOSE AND OPEN SCRATCH1 CONTAINING GENERATED PLOTEL IMAGES.

870 FILE = scrt1
    CALL CLOSE (scrt1,clsrew)
    CALL OPEN (*960,scrt1,z(buf4),rdrew)

    !     BLAST COPY PLOTELS FROM SCRATCH1 TO SCRATCH2.

880 CALL READ (*890,*890,scrt1,z(igrids),core-igrids,noeor,words)
    CALL WRITE (scrt2,z(igrids),core-igrids,noeor)
    GO TO 880
890 CALL WRITE (scrt2,z(igrids),words,eor)
900 CALL CLOSE (scrt1,clsrew)

    !     ALL PROCESSING OF GEOM2 IS COMPLETE SO COPY BALANCE OF GEOM2 TO
    !     SCRATCH2, WRAP UP, AND COPY BACK.

    CALL ifp4b (geom2,scrt2,any,z(igrids),core-igrids,-1,g2eof)

    !     ALL PROCESSING COMPLETE.

910 CALL CLOSE (axic,clsrew)
    CALL conmsg (msg2,2,0)
    RETURN

    !     END OF FILE ON AXIC.

920 FILE = axic
    GO TO 940

    !     END OF RECORD ON AXIC

930 FILE = axic
    ier  = -3
    GO TO 2000

    !     END OF FILE OR END OF RECORD ON -FILE-.

940 ier = -2
    GO TO 2000

    !     FILE NOT IN FIST

960 ier = -1
    GO TO 2000

    !     INSUFFICIENT CORE

980 ier = -8
    FILE = icrq
    GO TO 2000

    !     BISLOC EXIT

1610 ier = -37

2000 CALL mesage (ier,FILE,subr)
    RETURN
END SUBROUTINE ifp5
