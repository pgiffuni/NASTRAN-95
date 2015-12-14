SUBROUTINE premat (iz,rz,bfr,nimat,n2mat,mptf,ditf)
     
!     REVISED 7/92, BY G.CHAN, NEW REFERENCE TO OPEN CORE ARRAY, SUCH
!     THAT THE SOURCE CODE IS UP TO ANSI FORTRAN 77 STANADARD
 
 
 INTEGER, INTENT(IN OUT)                  :: iz(1)
 REAL, INTENT(IN OUT)                     :: rz(1)
 INTEGER, INTENT(IN OUT)                  :: bfr(1)
 INTEGER, INTENT(IN)                      :: nimat
 INTEGER, INTENT(OUT)                     :: n2mat
 INTEGER, INTENT(IN)                      :: mptf
 INTEGER, INTENT(IN OUT)                  :: ditf
 LOGICAL :: part1 ,pla   ,tdep
 INTEGER :: qmat1 ,qmat2 ,qmatx ,flag  ,  &
     mat1(2)      ,matt1(2)     ,dit   ,buf(3),tempid,  &
     nam(2),back  ,ret   ,mat2(2)      ,matt2(2)     ,  &
     ret1  ,pass  ,tablid,tablei(16)   ,mats1(2)     ,  &
     mat3(2)      ,matt3(2)     ,qmat3 ,qmat8 ,mat8(2)
 INTEGER :: qmatf ,elemid,qmtpz1,qmtpz2,qmat6 ,sysbuf,matf(2)
 INTEGER :: z     ,offset
 REAL :: nu    , x(27) ,y(25) ,nux   ,nuxx  ,j11   ,  &
     j12   ,j22   ,nuxy3 ,nuyz3 ,nuzx3 ,matset,zz(1)
 DIMENSION       matpz1(2)    ,matpz2(2)    ,mttpz1(2)    ,  &
     mttpz2(2)    ,bufpz(51)    ,mat6(2)      ,  &
     matt6(2)     ,buftm6(39)   ,xy(108)      , ib(46)
 COMMON /zzzzzz/ z(1)
 COMMON /system/ sysbuf,nout  ,skp(7),tempid
 COMMON /names / rd    ,rdrew ,wrt   ,wrtrew,clsrew
 COMMON /matin / matid ,inflag,temp  ,plaarg,sinth ,costh
 COMMON /matout/ e     ,g     ,nu    ,rho   ,alpha ,TO    ,GE    ,  &
     sigmat,sigmac,sigmas,SPACE(15)    ,tdep  , dum26(26)
 
!     COMMON FOR PIEZOELECTRIC MATERIAL
 COMMON /matpz / pzout(51)
 
!     ISOPARAMETRIC MATERIALS
!     COMMON /MATISO/ G11   ,G12   ,G13   ,G14   ,G15   ,G16   ,G12   ,
!                     G22,..,G56   ,G66   ,RHO   ,
!                     AXX   ,AYY   ,AZZ   ,AXY   ,AYZ   ,AZX   ,TREF  ,
!                     GE    ,IER
 
 COMMON /matiso/ bufm6(46)
 EQUIVALENCE     (e    ,buf(1),y(1) ), (y( 1),g11   ,ex3  ,plaans ,indstr),  &
     (y( 2),g12   ,ey3  ,icell2), (y( 3),g13   ,ez3  ),  &
     (y( 4),g22   ,nuxy3), (y( 5),g23   ,nuyz3),  &
     (y( 6),g33   ,nuzx3), (y( 7),rhoy  ,rho3 ),  &
     (y( 8),alph1 ,gxy3 ), (y( 9),alph2 ,gyz3 ),  &
     (y(10),alph12,gzx3 ), (y(11),toy   ,ax3  ),  &
     (y(12),gey   ,ay3  ), (y(13),sigty ,az3  ),  &
     (y(14),sigcy ,tref3), (y(15),sigsy ,ge3  )
 EQUIVALENCE     (y(16),j11   ,g113 ), (y(17),j12   ,g123 ),  &
     (y(18),j22   ,g133 ), (y(19),       g223 ),  &
     (y(20),       g233 ), (y(21),       g333 ),  &
     (y(22),     sigty3 ), (y(23),     sigcy3 ),  &
     (y(24),     sigsy3 ), (y(25),     matset )
 EQUIVALENCE     (x( 1)    , ex     ), (x( 2)    , gx     ),  &
     (x( 3)    , nux    ), (x( 4)    , rhox   ),  &
     (x( 5)    , alphx  ), (x( 6)    , tox    ),  &
     (x( 7)    , gex    ), (x( 8)    , sigtx  ),  &
     (x( 9)    , sigcx  ), (x(10)    , sigsx  )
 EQUIVALENCE     (bufm6(1) , ib(1)  )
 EQUIVALENCE     (zz(1)    , z(1)   )
 
!     DATA DEFINING CARDS TO BE READ
 
 DATA    mat1   ,  kmat1  ,  lmat1  /  103, 1,   11,     31 /  ,  &
     mat2   ,  kmat2  ,  lmat2  /  203, 2,   16,     46 /  ,  &
     mat3   ,  kmat3  ,  lmat3  / 1403,14,   16,     46 /  ,  &
     mat6   ,  kmat6  ,  lmat6  / 2503,25,   40,    118 /  ,  &
     mat8   ,  kmat8  ,  lmat8  /  603, 6,   18,     52 /  ,  &
     matt1                      /  703, 7               /  ,  &
     matt2                      /  803, 8               /  ,  &
     matt3                      / 1503,15               /  ,  &
     matt6                      / 2603,26               /  ,  &
     mats1                      /  503, 5               /  ,  &
     matpz1 ,  kmtpz1 ,  lmtpz1 / 1603,16,   15,     43 /  ,  &
     matpz2 ,  kmtpz2 ,  lmtpz2 / 1703,17,   52,    154 /  ,  &
     mttpz1                     / 1803,18               /  ,  &
     mttpz2                     / 1903,19               /  ,  &
     matf   ,  kmatf  ,  lmatf  / 5110,51,    3,      3 /
 DATA    ntypes,   tablei /  8,  &
     105 , 1,  205 , 2,  305 , 3,  405 , 4,  3105,31,  3205,32, 3305,33,  3405,34/
 DATA    nam    /  4HMAT ,   4H     /
 
!     MAT1 AND MAT2 HAVE ONE EXTRA WORD AT END
 
 DATA nwmat1 /  12 /  ,   nwmat2 / 17  /
 
!      - DATA IN /MATOUT/ IN VARIOUS MAT FORMATS, AND INFLAGS -
 
!     FORMAT  MAT1    MAT2                     MAT3      MAT6  MAT8 MATP
!     INFLAG=   1   2,3,12    4      5    6,8     7   11   10    12
!     -WORD- -----  ------  --- ------ ------ ----- ---- ---- ----- ----
!       1       E     G11   RHO INDSTR PLAANS    EX         :    E1    E
!       2       G     G12             ICDLL/8    EY         :  NU12    E
!       3      NU     G13                        EZ              E2
!       4     PHO     G22                      NUXY  RHO        G12
!       5    ALPH     G23                      NUYZ             G2Z
!       6      T0     G33                      NUZX             G1Z
!       7      GE     RHO                       RHO             RHO
!       8    SIGT   ALPH1    (X/N INDICATES     GXY           ALPH1
!       9    SIGC   ALPH2     ITEM X IS FOR     GYZ           ALPH2
!      10    SIGS  ALPH12     INFLAG=N ONLY)    GZX              TO
!      11              TO                        AX              TL
!      12              GE                        AY              CL
!      13            SIGT                        AZ              TT
!      14            SIGC                        TO              CT
!      15            SIGS                        GE              IS
!      16         E/12 J11/3                     G11             GE
!      17         E/12 J12/3                     G12            F12
!      18              J22/3                     G13
!      19                                        G22
!      20                                        G23
!      21                                        G33
!      22                                       SIGT
!      23                                       SIGC
!      24                                       SIGS
!      25  MATSET  MATSET                     MATSET          MATSET
!      26    TDEP
!       :
 
!     PERFORM GENERAL INITIALIZATION
 
 qmat1  = 0
 qmat2  = 0
 qmat3  = 0
 qmat6  = 0
 qmat8  = 0
 qmatf  = 0
 qmatx  = 0
 qmtpz1 = 0
 qmtpz2 = 0
 pla    = .false.
 IF (ditf < 0) pla = .true.
 part1  = .true.
 i      =-1
 mpt    = mptf
 dit    = IABS(ditf)
 offset = locfx(iz(1)) - locfx(z(1))
 IF (offset < 0) CALL errtrc ('premat  ',10)
 n1mat  = nimat + offset
 
!     READ MAT1,MAT2 AND MAT3 CARDS.  SPREAD FORMAT SO THAT MATTI AND
!     MATSI TEMPERATURE AND STRESS-STRAIN TABLE NUMBERS CAN BE MERGED
 
 CALL preloc (*350,bfr,mpt)
 imat1 = 1 + offset
 i     = 1 + offset
 CALL locate (*60,bfr,mat1,flag)
 qmat1 = 1
 imhere= 30
 30 CALL READ (*1350,*50,mpt,buf,nwmat1,0,flag)
 z(i) = buf(1)
 i    = i + 1
 DO  j = 2,kmat1
   z(i  ) = buf(j)
   z(i+1) = 0
   z(i+2) = 0
   IF (i > n1mat) GO TO 1398
   i = i + 3
 END DO
 GO TO 30
 50 nmat1 = i - lmat1
 60 imat2 = i
 CALL locate (*100,bfr,mat2,flag)
 qmat2 = 1
 imhere= 70
 70 CALL READ (*1350,*90,mpt,buf,nwmat2,0,flag)
 z(i) = buf(1)
 i    = i + 1
 DO  j = 2,kmat2
   z(i  ) = buf(j)
   z(i+1) = 0
   z(i+2) = 0
   IF (i > n1mat) GO TO 1398
   i = i + 3
 END DO
 GO TO 70
 90 nmat2 = i - lmat2
 100 imat3 = i
 CALL locate (*131,bfr,mat3,flag)
 qmat3 = 1
 imhere= 110
 110 CALL READ (*1350,*130,mpt,buf,kmat3,0,flag)
 z(i) = buf(1)
 i    = i + 1
 DO  j = 2,kmat3
   z(i  ) = buf(j)
   z(i+1) = 0
   z(i+2) = 0
   IF (i > n1mat) GO TO 1398
   i = i + 3
 END DO
 GO TO 110
 130 nmat3  = i - lmat3
 131 imtpz1 = i
 CALL locate (*135,bfr,matpz1,flag)
 qmtpz1 = 1
 imhere = 132
 132 CALL READ (*1350,*134,mpt,buf,kmtpz1,0,flag)
 z(i) = buf(1)
 i    = i + 1
 DO  j = 2,kmtpz1
   z(i  ) = buf(j)
   z(i+1) = 0
   z(i+2) = 0
   IF (i > n1mat) GO TO 1398
   i = i + 3
 END DO
 GO TO 132
 134 nmtpz1 = i - lmtpz1
 135 imtpz2 = i
 CALL locate (*140,bfr,matpz2,flag)
 qmtpz2 = 1
 imhere = 136
 136 CALL READ (*1350,*138,mpt,buf,kmtpz2,0,flag)
 z(i) = buf(1)
 i    = i + 1
 DO  j = 2,kmtpz2
   z(i  ) = buf(j)
   z(i+1) = 0
   z(i+2) = 0
   IF (i > n1mat) GO TO 1398
   i = i + 3
 END DO
 GO TO 136
 138 nmtpz2 = i - lmtpz2
 140 imat6  = i
 CALL locate (*144,bfr,mat6,flag)
 qmat6  = 1
 imhere = 141
 141 flag   = 0
 CALL READ (*1350,*143,mpt,buf,kmat6,0,flag)
 z(i) = buf(1)
 i    = i + 1
 DO  j = 2,kmat6
   z(i  ) = buf(j)
   z(i+1) = 0
   z(i+2) = 0
   IF (i > n1mat) GO TO 1398
   i = i + 3
 END DO
 GO TO 141
 143 nmat6 = i - lmat6
 IF (flag /= 0) GO TO 1570
 144 imat8 = i
 CALL locate (*1444,bfr,mat8,flag)
 qmat8 = 1
 imhere= 1441
 1441 flag  = 0
 CALL READ (*1350,*1443,mpt,buf,kmat8,0,flag)
 z(i) = buf(1)
 i    = i + 1
 DO  j = 2,kmat8
   z(i  ) = buf(j)
   z(i+1) = 0
   z(i+2) = 0
   IF (i > n1mat) GO TO 1398
   i = i + 3
 END DO
 GO TO 1441
 1443 nmat8 = i - lmat8
 IF (flag /= 0) GO TO 1570
 1444 imatf = i
 CALL locate (*149,bfr,matf,flag)
 qmatf  = 1
 imhere = 145
 145 CALL READ (*1350,*148,mpt,buf,kmatf,0,flag)
 z(i  ) = buf(1)
 z(i+1) = buf(2)
 z(i+2) = buf(3)
 i = i + 3
 GO TO 145
 148 nmatf = i - lmatf
 149 ilist = i
 IF (i > n1mat) GO TO 1398
 CALL CLOSE (mpt,clsrew)
 
!     IF TEMPERATURE OR PLA PROBLEM, READ THE MATTI OR MATSI CARDS.
!     MERGE MATSI AND MATTI DATA WITH MATI DATA.
!     SAVE A LIST OF TABLES REFERENCED.
 
 IF (     pla .AND. tempid /= 0) GO TO 1540
 IF (.NOT.pla .AND. tempid == 0) GO TO 350
 CALL preloc (*350,bfr,mpt)
 IF (tempid /= 0) GO TO 160
 nx = 3
 CALL locate (*150,bfr,mats1,flag)
 qmatx = 1
 ASSIGN  150 TO back
 ASSIGN 1420 TO ret1
 ASSIGN  820 TO pass
 n = kmat1
 GO TO 910
 150 CONTINUE
 160 nx = 2
 CALL locate (*170,bfr,matt1,flag)
 qmatx = 1
 ASSIGN  170 TO back
 ASSIGN 1450 TO ret1
 ASSIGN  820 TO pass
 n = kmat1
 GO TO 910
 170 CALL locate (*180,bfr,matt2,flag)
 qmatx = 1
 ASSIGN  180 TO back
 ASSIGN 1460 TO ret1
 ASSIGN  850 TO pass
 n = kmat2
 GO TO 910
 180 CALL locate (*181,bfr,matt3,flag)
 qmatx = 1
 ASSIGN  181 TO back
 ASSIGN 1520 TO ret1
 ASSIGN  880 TO pass
 n = kmat3
 GO TO 910
 181 CALL locate (*182,bfr,mttpz1,flag)
 qmatx = 1
 ASSIGN  182 TO back
 ASSIGN 1551 TO ret1
 ASSIGN  901 TO pass
 n = kmtpz1
 GO TO 910
 182 CALL locate (*183,bfr,mttpz2,flag)
 qmatx = 1
 ASSIGN  183 TO back
 ASSIGN 1552 TO ret1
 ASSIGN  904 TO pass
 n = kmtpz2
 GO TO 910
 183 CALL locate (*190,bfr,matt6,flag)
 qmatx = 1
 ASSIGN  190 TO back
 ASSIGN 1560 TO ret1
 ASSIGN  907 TO pass
 n = 31
 GO TO 910
 190 itabl = i
 imhere= 190
 IF (i > n1mat) GO TO 1398
 nlist = itabl - 11
 CALL CLOSE (mpt,clsrew)
 
!     IF ANY MATTI OR MATSI CARDS WERE READ, FORM A SORTED LIST OF TABLE
!     NUMBERS REFERENCED ON THESE CARDS. THEN, DISCARD ANY DUPLICATES IN
!     THE LIST SO THAT THE LIST CONTAINS UNIQUE TABLE NOS. TO BE READ.
 
 IF (qmatx == 0) GO TO 350
 DO  ii = ilist,nlist,11
   MIN = 999999999
   DO  jj = ii,nlist,11
     IF (z(jj) >= MIN) CYCLE
     MIN = z(jj)
     jx  = jj
   END DO
   z(jx) = z(ii)
   z(ii) = MIN
 END DO
 z(itabl) = 0
 jj = ilist
 DO  ii = ilist,nlist,11
   IF (z(ii+11) == z(ii)) CYCLE
   z(jj) = z(ii)
   jj = jj + 11
 END DO
 itabl = jj
 nlist = jj - 11
 
!     READ THE DIT BY TABLE TYPE. FOR EACH TABLE IN THE DIT, LOOK UP IN
!     TABLE NO. LIST TO DETERMINE IF THE TABLE IS REQUIRED FOR PROBLEM
!     SOLUTION. IF NOT, SKIP THE TABLE. IF SO, READ THE TABLE INTO CORE
!     AND STORE POINTERS TO THE FIRST AND LAST ENTRIES IN THE TABLE AND
!     THE TYPE OF TABLE. THIS INFORMATION IS STORED IN THE TABLE NO.
!     LIST
 
 CALL preloc (*1370,bfr,dit)
 i = itabl
 j = 1
 ASSIGN 260 TO ret
 ASSIGN 280 TO ret1
 240 jj = j + j - 1
 CALL locate (*290,bfr,tablei(jj),flag)
 250 CALL READ (*1380,*290,dit,buf,8,0,flag)
 nwds = 2
 IF (j == 4 .OR. j == 8) nwds = 1
 tablid = buf(1)
 GO TO 960
 260 z(l+1) = j
 IF (j > 4) z(l+1) = j - 4
 z(l+2) = i
 imhere = 270
 270 CALL READ (*1380,*1390,dit,z(i),nwds,0,flag)
 IF (z(i) == -1) GO TO 300
 i = i + nwds
 IF (i > n1mat) GO TO 1398
 GO TO 270
 280 CALL READ (*1380,*1390,dit,buf,nwds,0,flag)
 IF (buf(1) == -1) GO TO 250
 GO TO 280
 290 j = j + 1
 IF (j <= ntypes) GO TO 240
 CALL CLOSE (dit,clsrew)
 GO TO 330
 300 z(l+3) = i - nwds
 
!     STORE THE PARAMETERS ON THE TABLEI CARD IN LOCATIONS
!     Z(L+4),Z(L+5),...,Z(L+10)
 
 DO  k = 2,8
   lx = l + k
   z(lx+2) = buf(k)
 END DO
 
!     IF THIS TABLE IS A POLYNOMIAL (TABLE4), EVALUATE THE END POINTS
!     AND STORE AT ZZ(L+8) AND ZZ(L+9)
 
 IF (j /= 4) GO TO 250
 xx = (zz(l+6) - zz(l+4))/zz(l+5)
 ASSIGN 1330 TO igoto
 GO TO 1280
 320 ASSIGN 1340 TO igoto
 xx = (zz(l+7) - zz(l+4))/zz(l+5)
 GO TO 1280
 
!     TEST TO FOR ALL REFERENCED TABLES IN CORE
 
 330 flag = 0
 DO  l = ilist,nlist,11
   IF (z(l+1) /= 0) CYCLE
   flag = 1
   buf(1) = z(l)
   buf(2) = 0
   CALL mesage (30,41,buf)
 END DO
 IF (flag /= 0) CALL mesage (-37,0,nam)
 
!     WRAP UP PREMAT
 
 350 n2mat  = i + 1 - offset
 matido = 0
 sintho = 2.0
 costho = 2.0
 inflgo = 0
 part1  =.false.
 mapck  =+999
 RETURN
 
!     THE FOLLOWING POINTERS AND FLAGS ARE SET IN PREMAT FOR USE BY MAT-
!     QMAT1 = 0, NO MAT1 TABLE, = 1, MAT1 TABLE PRESENT
!     IMAT1 = POINTER TO FIRST ENTRY IN MAT1 TABLE
!     LMAT1 = LENGTH  OF EACH  ENTRY IN MAT1 TABLE
!     NAMT1 = POINTER TO LAST  ENTRY IN MAT1 TABLE
!     QMAT2,  IMAT2, LMAT2 AND NMAT2 ARE DEFINED AS ABOVE FOR MAT2 TABLE
!     QMATX = 0, NO TEMP OR STRESS TABLES PRESENT, = 1, OTHERWISW
!     ILIST = POINTER TO FIRST ENTRY IN TABLE LIST
!     NLIST = POINTER TO  LAST ENTRY IN TABLE LIST
 
!     THE TABLE LIST HAS 11 WORDS PER ENTRY AS FOLLOWS--
!      1. TABLE NUMBER (I.E. ID NO.)
!      2. TABLE TYPE (I.E. 1 = TABLE1, 2 = TABLE2 ETC.)
!      3. POINTER TO FIRST ENTRY IN TABLE
!      4. POINTER TO  LAST ENTRY IN TABLE
!      5. THRU 11. PARAMETERS ON TABLEI CARD
!     MATIDO = OLD MATERIAL ID (INITIALIZED TO  0 BY PREMAT)
!     SINTHO = OLD SIN THETA   (INITIALIZED TO 2. BY PREMAT)
!     INFLGO = OLD INFLAG      (INITIALIZED TO  0 BY PREMAT)
 
 
 
 ENTRY mat (elemid)
!     ==================
 
!     IF MAPCK .NE. +999 PREMAT HAS BEEN CORRUPTED.  (OVERLAY ERROR)
 
 IF (mapck == +999) GO TO 355
 WRITE  (nout,353) mapck
 353 FORMAT (//,' *** PREMAT OVERLEY ERROR',i12)
 CALL errtrc ('PREMAT  ',353)
 
 
!     INFLAG DESCRIBES PROCESSING REQUESTED BY CALLER
 
 355 SELECT CASE ( inflag )
   CASE (    1)
     GO TO 360
   CASE (    2)
     GO TO 400
   CASE (    3)
     GO TO 480
   CASE (    4)
     GO TO 560
   CASE (    5)
     GO TO 620
   CASE (    6)
     GO TO 640
   CASE (    7)
     GO TO 680
   CASE (    8)
     GO TO 2000
   CASE (    9)
     GO TO 2200
   CASE (   10)
     GO TO 2400
   CASE (   11)
     GO TO 2600
   CASE (   12)
     GO TO 2700
 END SELECT
 
!     INFLAG = 1 MEANS CALLER WANTS ONLY MAT1 PROPERTIES IN MAT1 FORMAT
!     IF NO TEMPERATURE DEPENDENT PROPERTIES AND MATID = OLD MATID,
!     RETURN SINCE DATA IS ALREADY IN MATOUT
 
 360 IF (tempid == 0 .AND. matid == matido .AND. inflag == inflgo .AND.  &
     .NOT. pla) RETURN
 matido = matid
 inflgo = inflag
 tdep   =.false.
 
!     LOOK UP MATID IN MAT1 TABLE
 
 ASSIGN  380 TO ret
 ASSIGN 1480 TO ret1
 GO TO 820
 
!     PICK UP MATERIAL PROPERTIES FROM MAT1 ENTRY.
 
 380 i = k + 1
 j = 1
 ASSIGN 390 TO back
 GO TO 980
 390 y(j) = prop
 i = i + 3
 j = j + 1
 IF (j < kmat1) GO TO 980
 RETURN
 
!     INFLAG = 2 MEANS CALLER WANTS MAT2 FORMAT WHETHER PROPERTIES ARE
!     DEFINED IN MAT1 OR MAT2 TABLE.
!     IF NO TEMPERATURE DEPENDENT PROPERTIES AND MATID = OLD MATID AND
!     SIN THETA = OLD SIN THETA, RETURN
 
 400 IF (tempid == 0 .AND.  matid == matido .AND. sinth == sintho .AND.  &
     .NOT.pla    .AND. inflag == inflgo .AND. costh == costho) RETURN
 inflgo = inflag
 matido = matid
 sintho = sinth
 costho = costh
 
!     LOOK UP MATID IN MAT1 TABLE
 
 410 ASSIGN 420 TO ret1
 ASSIGN 430 TO ret
 GO TO 820
 
!     MATID NOT IN MAT1 TABLE, LOOK UP IN MAT2 TABLE
!     - IF NOT PRESENT, FATAL ERROR IF INFLAG = 2
!     - IF NOT PRESENT, SEARCH MAT8 TABLE IF INFLAG = 12
 
 420 ASSIGN 450 TO ret
 ASSIGN 425 TO ret1
 GO TO 850
 425 IF (inflgo == 12) GO TO 2710
 GO TO 1480
 
!     MATID FOUND IN MAT1 TABLE.
!     COMPUTE G MATRIX FROM MAT1 PROPERTIES.
!     COMPLETE REMAINDER OF OUTPUT BUFFER IN MAT2 FORMAT.
 
 430 i = k + 1
 j = 1
 ASSIGN 440 TO back
 mmat = 1
 GO TO 980
 440 x(j) = prop
 i = i + 3
 j = j + 1
 IF (j < kmat1) GO TO 980
 nuxx  = 1.0 - nux**2
 g11   = ex/nuxx
 g12   = nux*g11
 g13   = 0.
 g22   = g11
 g23   = 0.
 g33   = gx
 rhoy  = rhox
 alph1 = alphx
 alph2 = alphx
 alph12= 0.
 toy   = tox
 gey   = gex
 sigty = sigtx
 sigcy = sigcx
 sigsy = sigsx
 IF (inflgo == 12) GO TO 2701
 RETURN
 
!     MATID FOUND IN MAT2 TABLE.
!     PLACE PROPERTIES IN OUTPUT AREA IN MAT2 FORMAT
!     THEN TEST FOR TRANSFORMATION. IF IDENTITY, RETURN. 3THERWISE,
!     PERFORM  U(T)*G*U .
 
 450 i = k + 1
 j = 1
 mmat = 2
 ASSIGN 460 TO back
 GO TO 980
 460 y(j) = prop
 i = i + 3
 j = j + 1
 IF (j    < kmat2) GO TO 980
 IF (inflgo ==  12) GO TO 2705
 IF (sinth  == 0.0) GO TO 470
 IF (ABS(sinth**2 + costh**2 - 1.0) > .0001) GO TO 1485
 
!     TRANSFORM G , THE MATERIAL STIFFNESS PROPERTY MATRIX.
!                M                   T
!                            G   =  U * G * U
!                             E          M
 
 x( 1) = costh**2
 x( 2) = sinth**2
 x( 3) = costh*sinth
 x( 4) = x(2)
 x( 5) = x(1)
 x( 6) =-x(3)
 x( 7) = 2.0*x(6)
 x( 8) =-x(7)
 x( 9) = x(1) - x(2)
 x(10) = g11
 x(11) = g12
 x(12) = g13
 x(13) = g12
 x(14) = g22
 x(15) = g23
 x(16) = g13
 x(17) = g23
 x(18) = g33
 CALL gmmats (x(10),3,3,0,x( 1),3,3,0,x(19))
 CALL gmmats (x( 1),3,3,1,x(19),3,3,0,x(10))
 g11   = x(10)
 g12   = x(11)
 g13   = x(12)
 g22   = x(14)
 g23   = x(15)
 g33   = x(18)
 
!     COMPUTE THE TRANSFORMED TEMPERATURE EXPANSION VECTOR
!               (ALPHA) = (U)*(ALPHA)
!                                    M
 
 x(3)  = -x(3)
 x(6)  = -x(6)
 x(7)  = -x(7)
 x(8)  = -x(8)
 CALL gmmats (x(1),3,3,0, y(8),3,1,0, x(10))
 alph1 = x(10)
 alph2 = x(11)
 alph12= x(12)
 470 IF (inflag == 7) GO TO 813
 RETURN
 
!     INFLAG = 3 IMPLIES THE CALLER WANTS
!             (1) ONLY J11, J12 AND J22, AND
!             (2) THE FIRST 15 LOCATIONS OF /MATOUT/ TO BE UNDISTURBED.
 
 480 IF (matid == matido .AND. inflag == inflgo .AND. .NOT.pla) RETURN
 IF (matid  /= matido) GO TO 490
 IF (inflgo /=      2) GO TO 490
 IF (mmat-2 < 0) THEN
   GO TO   530
 ELSE
   GO TO   540
 END IF
 
!     SEARCH MAT1 TABLE FOR MATID
 
 490 ASSIGN 500 TO ret1
 ASSIGN 510 TO ret
 GO TO 820
 
!     MATID NOT IN MAT1 TABLE. LOOK IN MAT2 TABLE - ERROR IF NOT PRESENT
 
 500 ASSIGN 1480 TO ret1
 ASSIGN  540 TO ret
 GO TO 850
 510 i = k + 4
 ASSIGN 520 TO back
 GO TO 980
 520 j11 = prop
 j12 = 0.0
 j22 = prop
 GO TO 550
 530 j11 = gx
 j12 = 0.0
 j22 = gx
 GO TO 550
 540 j11 = 0.0
 j12 = 0.0
 j22 = 0.0
 550 inflgo = inflag
 matido = matid
 RETURN
 
!     INFLAG = 4 MEANS CALLER DESIRES ONLY THE DENSITY PROPERTY (RHO)
!     LOOK UP MATID IN MAT1 TABLE.
 
 560 IF (tempid == 0 .AND. matid == matido .AND. inflag == inflgo .AND.  &
     .NOT.pla) RETURN
 ASSIGN 580 TO ret
 ASSIGN 570 TO ret1
 GO TO 820
 
!     MATID NOT IN MAT1 TABLE, LOOK UP IN MAT2 TABLE - ERROR IF NOT
!     PRESENT
 
 570 ASSIGN  610 TO ret
 ASSIGN 1480 TO ret1
 GO TO 850
 
!     MATID FOUND IN MAT1 TABLE. PICK UP RHO
 
 580 i = k + 10
 590 ASSIGN 600 TO back
 GO TO 980
 600 y(1)   = prop
 matido = matid
 inflgo = inflag
 RETURN
 
!     MATID FOUND IN MAT2 TABLE. PICK UP RHO.
 
 610 i = k + 19
 GO TO 590
 
!     INFLAG = 5, USED ONLY IN MODULE PLA1, DETERMINES IF THE MAT CARD
!     REFERENCED IS A MAT1 WITH E, YOUNGS MODULUS, DEFINED AS STRESS
!     DEPENDENT.  IF IT IS STRESS DEPENDENT, INDSTR, THE FIRST WORD OF
!     THE /MATOUT/ BLOCK IS SET = 1.  IF NOT STRESS DEPENDENT, INDSTR
!     IS SET = 0 ONLY MAT1 CARDS ARE ADMISSIBLE FOR THIS TEST.
 
 620 IF (pla .AND. matid == matido .AND. inflag == inflgo) RETURN
 matido = matid
 inflgo = inflag
 ASSIGN 630 TO ret
 ASSIGN 635 TO ret1
 indstr = 0
 GO TO 820
 
!     TEST TO SEE IF THE MATERIAL PROPERTY E IS DEPENDENT ON A TABLE OF
!     STRAIN VS. STRESS (EPSILON VS. SIGMA)
 
 630 tablid = z(k+3)
 IF (tablid /= 0) indstr = 1
 635 RETURN
 
!     INFLAG = 6, USED ONLY IN SUBROUTINES PLA3 AND PLA4, ACCEPTS
!     EPSILON - STRAIN - IN THE /MATIN/ BLOCK (PLAARG) AND LOOKS-UP
!     SIGMA   - STRESS - AND STORES THIS VALUE IN PLAANS IN /MATOUT/.
!     ONLY MAT1 AND MATS1 CARDS ARE ADMISSIBLE FOR THIS INFLAG.
 
 640 ASSIGN  650 TO ret
 ASSIGN 1500 TO ret1
 matido = matid
 inflgo = inflag
 GO TO 820
 650 tablid = z(k+3)
 IF (tablid <= 0) GO TO 1510
 xx = plaarg
 ASSIGN  660 TO ret
 ASSIGN 1490 TO ret1
 GO TO 960
 660 itype = z(l+1)
 IF (itype /= 1) GO TO 1530
 ASSIGN 670 TO iret
 GO TO 1080
 670 plaans = prop
 RETURN
 
!     INFLAG = 7, USED CURRENTLY ONLY BY BELL AEROSYSTEMS ELEMENTS,
!     IMPLIES THE USER WANTS HIS DATA IN MAT3 FORMAT.  IF THE MATID IS
!     FOUND IN THE MAT1 SET, THE DATA IS STORED IN MAT3 FORMAT.  IF NOT
!     FOUND IN THE MAT1 SET, THE MAT3 SET IS SEARCHED. IF NOT FOUND IN
!     THE MAT3 SET THE MAT2 SET IS SEARCHED. IF NOT FOUND HERE, A FATAL
!     ERROR EXISTS.
 
 680 IF (tempid == 0 .AND. inflag == inflgo .AND. matid == matido) RETURN
 inflgo = inflag
 matido = matid
 ASSIGN 690 TO ret
 685 CONTINUE
 ASSIGN 790 TO ret1
 GO TO 820
 690 i = k + 1
 j = 1
 ASSIGN 700 TO back
 GO TO 980
 700 SELECT CASE ( j )
   CASE (    1)
     GO TO 710
   CASE (    2)
     GO TO 720
   CASE (    3)
     GO TO 730
   CASE (    4)
     GO TO 740
   CASE (    5)
     GO TO 750
   CASE (    6)
     GO TO 760
   CASE (    7)
     GO TO 770
   CASE (    8)
     GO TO 771
   CASE (    9)
     GO TO 772
   CASE (   10)
     GO TO 773
 END SELECT
 710 ex3   = prop
 ey3   = prop
 ez3   = prop
 GO TO 780
 720 gxy3  = prop
 gyz3  = prop
 gzx3  = prop
 GO TO 780
 730 nuxy3 = prop
 nuyz3 = prop
 nuzx3 = prop
 GO TO 780
 740 rho3  = prop
 GO TO 780
 750 ax3   = prop
 ay3   = prop
 az3   = prop
 GO TO 780
 760 tref3 = prop
 GO TO 780
 770 ge3   = prop
 GO TO 780
 771 sigty3 = prop
 GO TO 780
 772 sigcy3 = prop
 GO TO 780
 773 sigsy3 = prop
 matset = 1.0
 RETURN
 
 780 j = j + 1
 i = i + 3
 GO TO 980
 
!     SEARCH FOR MATID IN THE MAT3 SET
 
 790 ASSIGN 800 TO ret
 ASSIGN 811 TO ret1
 GO TO 880
 
!     PICK UP MATERIAL PROPERTIES FROM MAT3 ENTRY
 
 800 i = k + 1
 j = 1
 ASSIGN 810 TO back
 GO TO 980
 810 y(j) = prop
 i = i + 3
 j = j + 1
 IF (j < kmat3) GO TO 980
 matset = 3.0
 RETURN
 
!     SEARCH FOR MATID IN THE MAT2 SET
 
 811 ASSIGN  812 TO ret
 ASSIGN 1480 TO ret1
 GO TO 850
 
!     GO TO INFLAG = 2 CODE TO PICK UP MAT2 PROPERTIES
 
 812 GO TO 450
 813 sigty3 = sigty
 sigcy3 = sigcy
 sigsy3 = sigsy
 tref3  = toy
 ge3    = gey
 ax3    = alph1
 ay3    = alph2
 az3    = alph12
 g113   = g11
 g123   = g12
 g133   = g13
 g223   = g22
 g233   = g23
 g333   = g33
 matset = 2.0
 RETURN
 
!     INFLAG = 8 IS USED ONLY BY TWO-DIMENSIONAL ELEMENTS IN PIECEWISE
!     LINEAR ANALYSIS.  HERE WE PERFORM AN INVERSE INTERPOLATION TO
!     OBTAIN STRAIN (EPS) GIVEN STRESS (TAU)
 
 2000 ASSIGN 2010 TO ret
 ASSIGN 1500 TO ret1
 matido = matid
 inflgo = inflag
 yy     = plaarg
 GO TO 820
 2010 tablid = z(k+3)
 IF (tablid <= 0) GO TO 1510
 ASSIGN 2020 TO ret
 ASSIGN 1490 TO ret1
 GO TO 960
 2020 itype = z(l+1)
 IF (itype /= 1) GO TO 1530
 
!     ROUTINE TO PERFORM INVERSE LINEAR INTERPOLATION OR EXTRAPOLATION.
 
 itabl = z(l+2)
 ntabl = z(l+3)
 up    = 1.0
 IF (zz(itabl) > zz(itabl+2)) up = -1.0
 kxx1 = itabl
 IF ((yy - zz(itabl+1))*up < 0.0) GO TO 2180
 kxx1 = ntabl - 2
 IF ((yy - zz(ntabl+1))*up <= 0.0) GO TO 2030
 IF (zz(ntabl+1) == zz(ntabl-1)) GO TO 2180
 2030 klo = 1
 khi = (ntabl - itabl)/2 + 1
 2090 kx  = (klo + khi + 1)/2
 kxx = (kx - 1)*2 + itabl
 IF ((yy - zz(kxx+1))*up < 0.0) THEN
   GO TO  2100
 ELSE IF ((yy - zz(kxx+1))*up == 0.0) THEN
   GO TO  2150
 ELSE
   GO TO  2110
 END IF
 2100 khi = kx
 GO TO 2120
 2110 klo = kx
 2120 IF (khi-klo /= 1) GO TO 2090
 kxx1 = 2*(klo-1) + itabl
 IF (kxx ==       kxx1) GO TO 2130
 IF (yy  == zz(kxx1+3)) GO TO 2140
 2130 plaans = (yy - zz(kxx1+1))*(zz(kxx1+2) - zz(kxx1))/(zz(kxx1+3)  &
     - zz(kxx1+1)) + zz(kxx1)
 2135 icell2 = 0
 RETURN
 
 2140 kxx = kxx1 + 2
 2150 IF (yy == zz(kxx-1)) GO TO 2160
 IF (yy == zz(kxx+3)) GO TO 2170
 plaans = zz(kxx)
 GO TO 2135
 2160 plaans = (zz(kxx) + zz(kxx-2))/2.0
 GO TO 2135
 2170 plaans = (zz(kxx) + zz(kxx+2))/2.0
 GO TO 2135
 
!     YY IS OUT OF THE RANGE OF THE FUNCTION, SET THE SECOND CELL OF
!     /MATOUT/ EQUAL TO ONE.
 
 2180 plaans = 0.0
 icell2 = 1
 RETURN
 
!     INFLAG = 9 IS USED ONLY BY TRAPAX AND TRIAAX WHEN PIEZOELECTRIC
!     MATERIALS ARE SELECTED.  WANT MATERIALS RETURNED INTO MATPZ2
!     FORMAT.
 
!     MATPZ1 CODE TRANSFORMS 1,2,3 MATERIAL DIRECTIONS INTO Z, THETA,
!     R = 0 DIRECTIONS, RESPECTIVLELY, AND INTERCHANGES 4TH AND 6TH ROWS
!     AND COLUMNS TO ACCOUNT FOR DIFFERENT SHEAR ORDERING.
!     ELEMENT ROUTINE WILL TRANSFORM FOR R-POLARIZATION
!     MATPZ2 CODE ASSUMES USER HAS PERFO-MED ALL TRANSFORMATIONS AS
!     EXPLAINED FOR MATPZ1
 
 2200 IF (tempid == 0 .AND. inflag == inflgo .AND. matid == matido) RETURN
 inflgo = inflag
 matido = matid
 
!     LOOK UP MATID IN MATPZ1 TABLE
 
 ASSIGN 2210 TO ret1
 ASSIGN 2220 TO ret
 GO TO 901
 
!     NOT IN MATPZ1, LOOK AT MATPZ2
 
 2210 ASSIGN 2300 TO ret
 ASSIGN 685 TO ret1
 GO TO 904
 
!     FOUND IN MATPZ1 - PUT OUT LIKE MATPZ2
 
 2220 i = k + 1
 j = 1
 ASSIGN 2230 TO back
 GO TO 980
 2230 bufpz(j) = prop
 i = i + 3
 j = j + 1
 IF (j < kmtpz1)  GO TO 980
 epso = 8.854E-12
 DO  ijk = 1,8
   bufpz(ijk) = bufpz(ijk)*1.e-12
 END DO
 se1   = (bufpz(4)-bufpz(1))*2.*bufpz(5)**2 - bufpz(2)*  &
     (bufpz(4)**2-bufpz(1)**2)
 se2   = 2.*bufpz(5)**2 - bufpz(2)*(bufpz(4) + bufpz(1))
 IF (se1 == 0. .OR. se2 == 0.)GO TO 1556
 ce11  =-(bufpz(5)**2-bufpz(1)*bufpz(2))/se1
 ce12  = (bufpz(5)**2-bufpz(2)*bufpz(4))/se1
 ce13  =  bufpz(5)/se2
 ce33  =-(bufpz(4)+bufpz(1))/se2
 ce44  = 1./bufpz(3)
 ce66  = 0.5/(bufpz(1) - bufpz(4))
 e15   = bufpz(8)*ce44
 e31   = bufpz(6)*(ce11+ce12) + bufpz(7)*ce13
 e33   = bufpz(6)*ce13*2. + bufpz(7)*ce33
 eps11 = bufpz(9)*epso
 eps33 = bufpz(10)*epso
 DO  ijk = 4,44
   pzout(ijk)= 0.
 END DO
 pzout( 1) = ce33
 pzout( 2) = ce13
 pzout( 3) = ce13
 pzout( 7) = ce11
 pzout( 8) = ce12
 pzout(12) = ce11
 pzout(16) = ce44
 pzout(19) = ce44
 pzout(21) = ce66
 pzout(22) = e33
 pzout(23) = e31
 pzout(24) = e31
 pzout(31) = e15
 pzout(38) = e15
 pzout(40) = eps33
 pzout(43) = eps11
 pzout(45) = eps11
 pzout(46) = bufpz(11)
 pzout(47) = bufpz(12)
 pzout(48) = bufpz(12)
 pzout(49) = bufpz(12)
 pzout(50) = bufpz(13)
 pzout(51) = bufpz(14)
 matset = 4.0
 RETURN
 
!     FOUND IN MATPZ2 FORMAT
 
 2300 i = k + 1
 j = 1
 ASSIGN 2310 TO back
 GO TO 980
 2310 pzout(j) = prop
 i = i + 3
 j = j + 1
 IF (j < kmtpz2) GO TO 980
 matset = 5.0
 RETURN
 
!     INFLAG = 10, USED CURRENTLY ONLY BY ISOPARAMETRIC SOLIDS IHEX1,2,3
!     IMPLIES CALLER WANTS HIS DATA IN MAT6 FORMAT STORED IN MATISO.
!     MATERIALS COULD BE ON MAT1 OR ON MAT6. IN EITHER CASE,MATERIALS
!     WILL BE COMPUTED FOR MAT6 OUTPUT. IF NOT FOUND ON MAT1 OR MAT6,
!     FATAL.
 
 2400 IF (tempid == 0 .AND. matid == matido .AND. inflag == inflgo) RETURN
 inflgo = inflag
 matido = matid
 tdep   =.false.
 
!     LOOK UP MATID IN MAT1 TABLE
 
 ASSIGN 2420 TO ret1
 ASSIGN 2430 TO ret
 GO TO 820
 
!     MATID NOT IN MAT1. CHECK MAT6
 
 2420 ASSIGN 2470 TO ret
 ASSIGN 1480 TO ret1
 GO TO 907
 
!     MATID FOUND IN MAT1 TABLE. COMPUTE G MATRIX,ETC.
 
 2430 ib(46) = 1
 i = k + 1
 j = 1
 ASSIGN 2440 TO back
 GO TO 980
 2440 x(j) = prop
 i = i + 3
 j = j + 1
 IF (j < kmat1) GO TO 980
 dd = (1.+nux)*(1.-2.*nux)
 IF (dd /= 0.) GO TO 2450
 ib(46) = 0
 RETURN
 
 2450 dd   =  ex*(1.-nux)/dd
 ddn1 = nux/(1.-nux)
 ddn2 = 0.5*(1.-2.*nux)/(1.-nux)
 DO  ijkl = 1,45
   bufm6(ijkl) = 0.
 END DO
 bufm6( 1) = dd
 bufm6( 2) = dd*ddn1
 bufm6( 3) = bufm6(2)
 bufm6( 7) = bufm6(2)
 bufm6( 8) = bufm6(1)
 bufm6( 9) = bufm6(2)
 bufm6(13) = bufm6(2)
 bufm6(14) = bufm6(2)
 bufm6(15) = bufm6(1)
 bufm6(22) = dd*ddn2
 bufm6(29) = bufm6(22)
 bufm6(36) = bufm6(22)
 bufm6(37) = rhox
 bufm6(38) = alphx
 bufm6(39) = alphx
 bufm6(40) = alphx
 bufm6(44) = tox
 bufm6(45) = gex
 RETURN
 
!     MATID FOUND IN MAT6 TABLE. PUT PROPERTIES IN MAT6 FORMAT AND
!     TRANSFORM USING DIRECTION COSINES
 
 2470 ib(46) = 6
 i = k + 1
 j = 1
 ASSIGN 2480 TO back
 GO TO 980
 2480 buftm6(j) = prop
 i = i + 3
 j = j + 1
 IF (j < kmat6) GO TO 980
 
!     PUT SYMMETRIC PORTION OF G INTO A FULL 6 X 6 AND CREATE A 6 X 6
!     DIRECTION COSINE MATRIX BY COOK PP. 212-213. THEN TRANSFORM
!     (U-TRANSPOSE)*G*U
 
 kkk = 0
 lll = 0
 DO  iii = 1,6
   DO  jjj = iii,6
     kkk = kkk + 1
     lll = lll + 1
     xy(lll) = buftm6(kkk)
     IF (jjj == iii) CYCLE
     l5   = 5*(jjj-iii)
     isub = lll + l5
     xy(isub) = xy(lll)
   END DO
   lll = lll + iii
 END DO
 xl1 = buftm6(31)
 xm1 = buftm6(32)
 xn1 = buftm6(33)
 xl2 = buftm6(34)
 xm2 = buftm6(35)
 xn2 = buftm6(36)
 xl3 = buftm6(37)
 xm3 = buftm6(38)
 xn3 = buftm6(39)
 xy(37) = xl1**2
 xy(38) = xm1**2
 xy(39) = xn1**2
 xy(40) = xl1*xm1
 xy(41) = xm1*xn1
 xy(42) = xn1*xl1
 xy(43) = xl2**2
 xy(44) = xm2**2
 xy(45) = xn2**2
 xy(46) = xl2*xm2
 xy(47) = xm2*xn2
 xy(48) = xn2*xl2
 xy(49) = xl3**2
 xy(50) = xm3**2
 xy(51) = xn3**2
 xy(52) = xl3*xm3
 xy(53) = xm3*xn3
 xy(54) = xn3*xl3
 xy(55) = xl1*xl2*2.
 xy(56) = xm1*xm2*2.
 xy(57) = xn1*xn2*2.
 xy(58) = xl1*xm2 + xl2*xm1
 xy(59) = xm1*xn2 + xm2*xn1
 xy(60) = xn1*xl2 + xn2*xl1
 xy(61) = xl2*xl3*2.
 xy(62) = xm2*xm3*2.
 xy(63) = xn2*xn3*2.
 xy(64) = xl2*xm3 + xl3*xm2
 xy(65) = xm2*xn3 + xm3*xn2
 xy(66) = xn2*xl3 + xn3*xl2
 xy(67) = xl3*xl1*2.
 xy(68) = xm3*xm1*2.
 xy(69) = xn3*xn1*2.
 xy(70) = xl3*xm1 + xl1*xm3
 xy(71) = xm3*xn1 + xm1*xn3
 xy(72) = xn3*xl1 + xn1*xl3
 
 CALL gmmats (xy(1),6,6,0,xy(37),6,6,0,xy(73))
 CALL gmmats (xy(37),6,6,1,xy(73),6,6,0,bufm6(1))
 
!     MUST ALSO TRANSFORM THERMAL EXPANSION VECOT= (U-INVERSE)*ALPHA
!     BY COOK P.212, THE INVERSE OF U IS THE TRANSPOSE OF THE
!     MATRIX WHICH TRANSFORMS STRESSES
 
 kkk = 72
 DO  iii = 1,6
   DO  jjj = 1,36,6
     kkk = kkk + 1
     lll = jjj + iii + 35
     xy(kkk) = xy(lll)
   END DO
 END DO
 DO  iii = 75,87,6
   DO  jjj = 1,3
     kkk = iii + jjj
     xy(kkk) = xy(kkk)*0.5
   END DO
 END DO
 DO  iii = 90,102,6
   DO  jjj = 1,3
     kkk = iii + jjj
     xy(kkk) = xy(kkk)*2.0
   END DO
 END DO
 
 CALL gmmats (xy(73),6,6,0,buftm6(23),6,1,0,bufm6(38))
 
 bufm6(37) = buftm6(22)
 bufm6(44) = buftm6(29)
 bufm6(45) = buftm6(30)
 RETURN
 
!     INFLAG = 11 IS USED ONLY BY A HYDROELASTIC ANALYSIS TO FIND THE
!     DENSITY FOR THREE DIMENSIONAL FLUID ELEMENTS FROM MATF CARDS.
 
 2600 IF (qmatf == 0) GO TO 1480
 DO  k = imatf,nmatf,lmatf
   IF (z(k) == matid) GO TO 2620
 END DO
 GO TO 1480
 2620 rho = zz(k+1)
 RETURN
 
!     INFLAG = 12 IS USED ONLY BY SHELL ELEMENTS QUAD4 AND TRIA3.
!     MAT1 IS FIRST SEARCHED, IF NOT FOUND, MAT2 IS SEARCHED. IF FOUND
!     IN EITHER CASE, /MATOUT/ WILL BE FILLED WITH MAT2 FORMAT DATA.
!     IF NOT FOUND IN MAT1 OR MAT2, MAT8 IS SEARCHED AND MAT8 FORMAT IS
!     USED IN /MATOUT/. FATAL ERROR IF MAT8 IS NOT FOUND.
 
 2700 IF (tempid == 0 .AND. matid == matido .AND. inflag == inflgo .AND.  &
     sinth == sintho .AND. costh == costho .AND. .NOT.pla)  RETURN
 inflgo = inflag
 matido = matid
 sintho = sinth
 costho = costh
 
!     GO TO INFLAG = 2 CODE TO PICK UP MAT1 OR MAT2 PROPERTIES
!     SET MATSET TO 1.0 IF PROPERTY DATA COMES FROM MAT1, OR
!     TO 2.0 IF FROM MAT2
 
 GO TO 410
 
 2701 matset = 1.0
 y(16)  = ex
 y(17)  = ex
 RETURN
 2705 matset = 2.0
 RETURN
 
!     NOT FOUND IN MAT1 AND MAT2.  LOOK FOR MAT8, ERROR IF NOT FOUND
 
 2710 IF (qmat8 == 0) GO TO 1480
 DO  k = imat8,nmat8,lmat8
   IF (z(k) == matid) GO TO 2730
 END DO
 GO TO 1480
 2730 i = k + 1
 j = 1
 ASSIGN 2740 TO back
 GO TO 980
 
!     OUTPUT IN MAT8 FORMAT AND SET MATSET TO 8.0
 
 2740 x(j) = prop
 i = i + 3
 j = j + 1
 IF (j < kmat8) GO TO 980
 DO  k = 1,17
   y(k) = x(k)
 END DO
 y(2) = x(3)
 y(3) = x(2)
 y(5) = x(6)
 y(6) = x(5)
 matset = 8.0
 RETURN
 
 
!     INTERNAL ROUTINE TO SEARCH FOR MATERIAL IN MAT1 TABLE
 
 820 IF (qmat1 == 0) GO TO 840
 DO  k = imat1,nmat1,lmat1
   IF (z(k) == matid) GO TO ret, ( 380,430,510,580,630,650,690,930,  &
       2010,2430)
 END DO
 840 GO TO ret1, (420,500,570,635,790,1420,1450,1480,1500,2420)
 
!     INTERNAL ROUTINE TO SEARCH FOR MATERIAL IN MAT2 TABLE
 
 850 IF (qmat2 == 0) GO TO 870
 DO  k = imat2,nmat2,lmat2
   IF (z(k) == matid) GO TO ret, (930,450,610,540,812)
 END DO
 870 GO TO ret1, (425,1460,1480)
 
!     INTERNAL ROUTINE TO SEARCH FOR MATERIAL IN MAT3 TABLE.
 
 880 IF (qmat3 == 0) GO TO 900
 DO  k = imat3,nmat3,lmat3
   IF (z(k) == matid) GO TO ret, (800,930)
 END DO
 900 GO TO ret1, (811)
 
!     PIEZOELECTRIC MATERIALS
 
 901 IF (qmtpz1 == 0)  GO TO 903
 DO  k = imtpz1,nmtpz1,lmtpz1
   IF (z(k) == matid) GO TO ret, (2220,930)
 END DO
 903 GO TO ret1, (1551,2210)
 904 IF (qmtpz2 == 0)  GO TO 906
 DO  k = imtpz2,nmtpz2,lmtpz2
   IF (z(k) == matid) GO TO ret, (2300,930)
 END DO
 906 GO TO ret1, (1552,1480,685)
 
!     SEARCH FOR MATERIAL IN MAT6 TABLE(ISOPARAMETRIC SOLIDS)
 
 907 IF (qmat6 == 0) GO TO 909
 DO  k = imat6,nmat6,lmat6
   IF (z(k) == matid) GO TO ret, (2470,930)
 END DO
 909 GO TO ret1, (1560,1480)
 
!     INTERNAL ROUTINE TO READ MATXI CARDS, MERGE DATA IN MATI TABLE
!     AND STORE TABLE IDS IN CORE.
 
 910 ASSIGN 930 TO ret
 920 CALL READ (*1350,*950,mpt,buf,n,0,flag)
 matid = buf(1)
 GO TO pass, (820,850,880,901,904,907)
 930 DO  j = 2,n
   IF (buf(j) == 0) CYCLE
   jx    = k + 3*(j-2) + nx
   z(jx) = buf(j)
   z(i ) = buf(j)
   z(i+1)= 0
   i     = i + 11
 END DO
 GO TO 920
 950 GO TO back, (150,170,180,181,182,183,190)
 
!     INTERNAL ROUTINE TO SEARCH FOR A TABLE IN THE TABLE LIST
 
 960 DO  l = ilist,nlist,11
   IF (z(l) == tablid) GO TO ret, (260,660,990,1030,2020)
 END DO
 GO TO ret1, (280,1490,1520)
 
!     ROUTINE TO TEST FOR DEPENDENCE OF A MATERIAL PROPERTY ON
!     TEMPERATURE OR STRESS. IF DEPENDENT, APPROPRIATE TABLE LOOK UP
!     PROCEDURE IS EMPLOYED. IN EITHER CASE, THE PROPERTY IS RETURNED
!     IN PROP.
 
 980 IF (qmatx == 0) GO TO 1060
 flag   = 0
 tablid = z(i+1)
 IF (elemid < 0) GO TO 1020
 IF (tablid == 0) GO TO 1060
 xx   = temp
 tdep =.true.
 flag = 1
 ASSIGN  990 TO ret
 ASSIGN 1490 TO ret1
 GO TO 960
 990 ASSIGN 1010 TO ret
 1000 itype = z(l+1)
 SELECT CASE ( itype )
   CASE (    1)
     GO TO 1180
   CASE (    2)
     GO TO 1200
   CASE (    3)
     GO TO 1230
   CASE (    4)
     GO TO 1240
 END SELECT
 1010 propt = prop
 GO TO 1070
 
!     SINCE THIS IS NOT A PIECEWISE LINEAR ANALYSIS PROBLEM, NO STRESS
!     DEPENDENT MATERIAL PROPERTIES ARE ALLOWED.  IF AND WHEN THIS
!     RESTRICTION IS LIFTED THE FOLLOWING CODE CAN BE IMPLEMENTED.
!     CURRENTLY A TRANSFER IS ALWAYS MADE TO STATEMENT 1060, SINCE THE
!     ELEMENT ID. IS ALWAYS POSITIVE.
 
 1020 IF (pla) GO TO 1550
 IF (elemid > 0) GO TO 1060
 tablid = z(i+2)
 IF (tablid == 0) GO TO 1050
 ASSIGN 1030 TO ret
 ASSIGN 1490 TO ret1
 GO TO 960
 1030 ASSIGN 1040 TO ret
 xx = plaarg
 GO TO 1000
 1040 IF (flag /= 0) prop = prop*propt
 GO TO 1070
 1050 IF (flag /= 0) GO TO 1070
 1060 prop = zz(i)
 1070 GO TO back, ( 390,440,460,520,600,700,810,2230,2310,2440,2480, 2740)
 
!     ROUTINE TO PERFORM LINEAR INTERPOLATION FOR FUNCTION IN TABLE.
!     L POINTS TO THE ENTRY IN THE TABLE LIST WHICH DEFINES THE TABLE.
!     ARGUMENT IS XX. FUNCTION VALUE IS RETURNED IN PROP. EXTRAPOLATION
!     IS MADE IF XX IS OUTSIDE THE LIMITS OF THE TABLE.
 
 1080 itabl = z(l+2)
 ntabl = z(l+3)
 up    = 1.0
 IF (zz(itabl) > zz(itabl+2)) up = -1.0
 kxx1 = itabl
 IF ((xx - zz(itabl))*up <= 0.) GO TO 1130
 kxx1 = ntabl - 2
 IF ((xx - zz(ntabl))*up >= 0.) GO TO 1130
 klo = 1
 khi = (ntabl-itabl)/2 + 1
 1090 kx  = (klo+khi+1)/2
 kxx = (kx-1)*2 + itabl
 IF ((xx - zz(kxx))*up < 0.0) THEN
   GO TO  1100
 ELSE IF ((xx - zz(kxx))*up == 0.0) THEN
   GO TO  1150
 ELSE
   GO TO  1110
 END IF
 1100 khi = kx
 GO TO 1120
 1110 klo = kx
 1120 IF (khi-klo /= 1) GO TO 1090
 kxx1 = (klo-1)*2 + itabl
 IF (kxx == kxx1) GO TO 1130
 IF (xx  == zz(kxx1+2)) GO TO 1140
 1130 prop = (xx - zz(kxx1))*(zz(kxx1+3) - zz(kxx1+1))/(zz(kxx1+2)  &
     - zz(kxx1)) + zz(kxx1+1)
 GO TO iret, (670,1190,1220)
 1140 kxx = kxx1 + 2
 1150 IF (xx == zz(kxx-2)) GO TO 1160
 IF (xx == zz(kxx+2)) GO TO 1170
 prop = zz(kxx+1)
 GO TO iret, (670,1190,1220)
 1160 prop = (zz(kxx-1) + zz(kxx+1))/2.0
 GO TO iret, (670,1190,1220)
 1170 prop = (zz(kxx+1) + zz(kxx+3))/2.0
 
!     TABLE TYPE = 1
!     ARGUMENT = XX
 
 1180 ASSIGN 1190 TO iret
 GO TO 1080
 1190 GO TO ret, (1010,1040)
 
!     TABLE TYPE = 2
!     ARGUMENT = (XX-X1)
 
 1200 xx = xx - zz(l+4)
 1210 ASSIGN 1220 TO iret
 GO TO 1080
 1220 prop = zz(i)*prop
 GO TO ret, (1010,1040)
 
!     TABLE TYPE = 3
!     ARGUMENT = (XX-X1)/X2
 
 1230 xx = (xx - zz(l+4))/zz(l+5)
 GO TO 1210
 
!     TABLE TYPE = 4
!     PERFORM POLYNOMIAL INTERPOLATION
 
 
!     NOTE...
!         ZZ(L+4) = X1
!         ZZ(L+5) = X2
!         ZZ(L+6) = X3
!         ZZ(L+7) = X4
!         ZZ(L+8) = F((X3-X1)/X2)
!         ZZ(L+9) = F((X4-X1)/X2)
!         WHERE X1 AND X2 ARE TRANSLATION AND SCALE FACTORS RESPECTIVELY
!         AND X3 AND X4 (X3 .LT. X4) ARE THE END POINTS OF THE
!         INTERVAL OVER WHICH THE POLYNOMIAL IS DEFINED.
 
 1240 factor = zz(i)
 
!     DETERMINE THE ARGUMENT XX
 
 xx = (xx - zz(l+4))/zz(l+5)
 IF   (xx - (zz(l+6) - zz(l+4))/zz(l+5) > 0.0) THEN
   GO TO  1260
 END IF
 1250 prop = zz(l+8)
 GO TO 1310
 1260 IF   (xx - (zz(l+7) - zz(l+4))/zz(l+5) < 0.0) THEN
   GO TO  1280
 END IF
 1270 prop = zz(l+9)
 GO TO 1310
 1280 nn   = z(l+3)
 prop = zz(nn)
 1290 IF (nn <= z(l+2)) GO TO 1300
 prop = prop*xx + zz(nn-1)
 nn   = nn - 1
 GO TO 1290
 1300 IF (part1) GO TO igoto, (1330,1340)
 1310 prop = prop*factor
 GO TO ret, (1010,1040)
 1330 zz(l+8) = prop
 GO TO 320
 1340 zz(l+9) = prop
 GO TO 250
 
!     FATAL ERROR MESSAGES
 
 1350 n = -2
 dit = mpt
 GO TO 1400
 1370 n = -1
 GO TO 1400
 1380 n = -2
 GO TO 1400
 1390 n = -3
 GO TO 1400
 1385 WRITE  (nout,1386) imhere,i,n1mat,offset,nimat
 1386 FORMAT ('0*** NIMAT SPACE TOO SMALL.  ERROR AT',i5,'/PREMAT', /5X,  &
     'I,N1MAT,OFFSET,NIMAT =',3I12,i7,/)
 GO TO 1472
 1398 IF (nimat <= 2*sysbuf+4) GO TO 1385
 n = -8
 dit = i - n1mat
 1400 CALL mesage (n,dit,nam)
 n = 16
 buf(1) = 0
 buf(2) = 0
 GO TO 1470
 1420 n =  17
 1430 buf(1) = matid
 buf(2) = 0
 GO TO 1470
 1450 n =  19
 GO TO 1430
 1460 n =  20
 GO TO 1430
 
 1470 CALL sswtch (20,j)
 IF (j == 0) GO TO 1475
 WRITE  (nout,1471) buf(1),buf(2)
 1471 FORMAT (' PREMAT/1471 - BUF(1),BUF(2) =',2I10)
 1472 CALL errtrc ('MAT     ',1472)
 
 1475 CALL mesage (-30,n,buf)
 
 1480 n = 42
 buf(1) = elemid
 buf(2) = matid
 GO TO 1470
 1485 n = 103
 buf(1) = 0
 buf(2) = 0
 GO TO 1470
 1490 n = 112
 buf(1) = tablid
 GO TO 1470
 1500 n = 113
 GO TO 1430
 1510 n = 116
 buf(1) = matid
 buf(2) = tablid
 GO TO 1470
 1520 n = 114
 GO TO 1430
 1530 n = 115
 buf(1) = tablid
 buf(2) = itype
 GO TO 1470
 1540 buf(1) = tempid
 1541 buf(2) = 0
 n = 117
 GO TO 1470
 1550 buf(1) = elemid
 GO TO 1541
 1551 buf(2) = 1
 GO TO 1553
 1552 buf(2) = 2
 1553 n = 216
 1554 buf(1) = matid
 GO TO 1470
 1556 n = 214
 GO TO 1554
 1560 n = 217
 GO TO 1554
 1570 n = 219
 GO TO 1554
END SUBROUTINE premat
