SUBROUTINE ihexd (TYPE)
     
!     DOUBLE PRECISION VERSION
 
!     THIS ROUTINE PROCESSES IHEX1, IHEX2, AND IHEX3 ELEMENT DATA TO
!     PRODUCE STIFFNESS AND MASS MATRICES.  IF THE HEAT TRANSFER OPTION
!     IS ON, CONDUCTIVITY AND CAPACITY MATRICES ARE PRODUCED.  IF THE
!     DISPLACEMENT VECTOR POINTER IS NON-ZERO, THE DIFFERENTIAL
!     STIFFNESS MATRIX ONLY IS PRODUCED.
 
!           TYPE = 1    IHEX1
!           TYPE = 2    IHEX2
!           TYPE = 3    IHEX3
 
!           THE EST ENTRIES ARE
 
!     NAME  ----------INDEX----------   DESCRIPTION
!            IHEX1    IHEX2    IHEX3
 
!     EID        1        1        1    ELEMENT ID NO.
!     SIL      2-9     2-21     2-33    SCALAR INDEX LIST
!     MID       10       22       34    MATERIAL ID NO.
!     CID       11       23       35    MATERIAL COORD. SYSTEM ID NO.
!     NIP       12       24       36    NO. INTEGRATION POINTS PER EDGE
!     MAXAR     13       25       37    MAX ASPECT RATIO
!     ALFA      14       26       38    MAX ANGLE FOR NORMALS
!     BETA      15       27       39    MAX ANGLE FOR MIDSIDE POINTS
!     BGPDT  16-47   28-107   40-167    BASIC GRID POINT DATA
!     GPT    48-55  108-127  168-199    GRID POINT TEMPERATURES
 
!     - INSTALLATION NOTE --
!     GPTLD IS SUPPOSED TO CONTAIN GRID POINT TEMPERATURE LOADS FOR
!     COMPUTING DIFFERENTIAL STIFFNESS.  FOR INSTALLATION, GPTLD MUST
!     BE LOADED WITH DATA BY EMG.  IF GPTLD(1)=-1, NO TEMP LOAD IS
!     ASSUMED.
 
 
 INTEGER, INTENT(IN)                      :: TYPE
 LOGICAL :: anis       ,rect       ,tdep       ,diag       ,  &
     mtdep      ,heat1      ,nogo       ,nocstm
 INTEGER :: heat       ,eid        ,sil(1)     ,scr4       ,  &
      jz(1)      ,cid        ,iest(1)    ,  &
     bcord      ,bgpdt      ,gpt        ,nc(8)      ,  &
     edge       ,face       ,ib(46)     ,elno(3)    ,  &
     excd(3)    ,twins(9)   ,rvrs(5)    ,iwork(1)   ,  &
     back       ,otpt       ,ugv        ,cdamp      , dict(40)
 REAL :: nu         ,kheat      ,maxar      ,dmaxar(3)  ,  &
     dalfa(3)   ,dbeta(2)   ,evec(3,12) ,work(66)   ,  &
     vn(3,2)    ,gptld(32)  ,bcd2(3)
 DOUBLE PRECISION :: z(1)       ,jacob(3,3) ,detj       ,s(4)       ,  &
     h(4)       ,gauss(8)   ,sfact      ,part(3,3)  ,  &
     e1         ,e2         ,e3         ,tf(3,3)    ,  &
     tk(3,3)    ,prt1       ,sig(6)     ,sx         ,  &
     sy         ,sz         ,sxy        ,syz        ,  &
     szx        ,str(18)    ,c(3,3)     ,temp
 DOUBLE PRECISION :: gmat(36)   ,dalpha(6)  ,store(45)  ,tvol
 CHARACTER (LEN=23) :: ufm
 COMMON /xmssg /  ufm
 COMMON /matin/   mid        ,inflag     ,temp
 COMMON /matout/  e          ,g          ,nu         ,rho        ,  &
     talpha     ,tref       ,cdamp      ,SPACE(18)  , mtdep
 COMMON /matiso/  bufm6(46)
!     OMMON  /MATISO/  G11,G12,G13,...,G46,G56,G66,RHO,AXX,AYY,AZZ,AXY,
!                      AYZ,AZX,TREF,GE,IER
 COMMON /BLANK /  skip16(16) ,volume     ,surfac
 COMMON /hmtout/  kheat(6)   ,cp
!     COMMON /EMG***/  ...,UGV,...
 
!     - INSTALLATION NOTE --
!     UGV POINTS TO BEGINNING OF SINGLE PRECISION GLOBAL DISPLACEMENT
!     VECTOR IN OPEN CORE ARRAY RZ.
 
 COMMON /emgprm/  iext       ,izs        ,nzs         ,dum(12)   ,  &
     kgg1       ,mgg1       ,bgg1        ,iprec     , nogo       ,heat1
 
!     SZ IS OPEN CORE.  USE ONLY RZ(IZS) TO RZ(NZS).
 
 COMMON /zzzzzz/  rz(1)
 COMMON /emgest/  est(200)
 COMMON /system/  sysbuf     ,otpt       ,sys1(7)    ,mtemp
 COMMON /emgdic/  spac(2)    ,ngrids     ,spac1      ,iestid
 EQUIVALENCE      (z(1),jz(1),rz(1))     ,(eid,est(1),iest(1))   ,  &
     (sil(1),est(2))        ,(work(1),iwork(1))     ,  &
     (sig(1),sx)            ,(sig(2),sy)            ,  &
     (sig(3),sz)            ,(sig(4),sxy)           ,  &
     (sig(5),syz)           ,(sig(6),szx)           , (dstld,idstld)
 EQUIVALENCE      (work(1),evec(1,1))    ,(work(37),vn(1,1))     ,  &
     (work(43),nc(1))
 EQUIVALENCE      (work(1),jacob(1,1))   ,(work(19),h(1))        ,  &
     (work(27),s(1))        ,(work(35),part(1,1))   ,  &
     (work(53),sig(1))      ,(work(1),c(1,1))
 EQUIVALENCE      (work(1),tf(1,1))      ,(work(35),tk(1,1))
 EQUIVALENCE      (ib(1),bufm6(1))
 DATA    scr4  /  304 /
 DATA    bcd1  ,  bcd2/ 4HCIHE, 4HX1  , 4HX2  , 4HX3   /
 DATA    dmaxar,  dalfa,dbeta / 5.0      ,10.0       ,15.0       ,  &
     45.0      ,45.0       ,45.0       , 45.0       ,45.0       /
 DATA    dtor  ,  gauss /0.017453292519943E0, 0.577350269189626D0,  &
     0.555555555555556D0, 0.774596669241483D0,  &
     0.888888888888889D0, 0.347854845137454D0,  &
     0.861136311594053D0, 0.652145154862546D0,  &
     0.339981043584856D0/
 DATA    ihex,elno       /4HIHEX,4H ele,4HMENT,4H no./
 DATA    bar,balfa,bbeta /4H  ar,4HALFA,4HBETA/
 DATA    excd            /4H exc,4HEEDE,4HD.  /
 DATA    rvrs            /4HREVE,4HRSED,4H num,4HBERI,4HNG. /
 DATA    twins           /4HCOOR,4HDINA,4HTES ,4HOF t,4HWO p,  &
     4HOINT,4HS ar,4HE sa,4HME. /
 DATA    nerr1,nerr2     /3301, 3302 /
 
!     FOR DOUBLE PRECISION, OPEN CORE POINTERS MUST BE MODIFIED
 
 iz = izs/2 + 1
 nz = nzs/2 + 1
 
!     THIS ROUTINE OPERATES IN DOUBLE PRECISION.
!     EMGOUT WILL PRODUCE THE REQUIRED MATRIX IN THE REQUESTED PRECISION
 
!     ALLOCATE LARGE ARRAYS IN OPEN CORE
 
 ngp  = 12*TYPE - 4
 heat = 0
 kgg  = 0
 mgg  = 0
 IF (heat1) heat = 1
 IF (kgg1 /= 0) kgg = 1
 IF (mgg1 /= 0) mgg = 1
 ngrids = ngp
 ugv  = 0
 ngg  = 3*ngp
 dict(1) = iestid
 dict(2) = 1
 IF (.NOT.heat1) GO TO 5
 dict(3) = ngp
 dict(4) = 1
 GO TO 30
 5 dict(3) = ngg
 dict(4) = 7
 IF (kgg <= 0) GO TO 10
 ik = iz + 3*ngg
 nk = ik - 1 + (ngg+1)*ngg/2
 GO TO 20
 10 ik = iz
 nk = ik + 3*ngg - 1
 im = nk + 1
 nm = (ngp+1)*ngp/2 + nk
 GO TO 40
 20 nm = nk
 IF (mgg <= 0) GO TO 40
 im = nk + 1
 nm = nk + (ngp+1)*ngp/2
 GO TO 40
 30 ik = iz + 17
 nk = ik - 1 + ngp**2
 im = nk + 1
 nm = im - 1 + ngp**2
 ngg= ngp
 40 in = nm + 1
 ig = in + ngp
 ix = ig + 3*ngp
 nd = nm + 9*ngp
 IF (ugv == 0) GO TO 50
 id = nd + 1
 nd = id + ngg - 1
 50 IF (nd <= nz) GO TO 100
 WRITE (otpt,7100) ufm,nerr1,ihex,TYPE,elno,eid
 nogo = .true.
 
!     ***** OPEN CORE MAP *****
 
!     DOUBLE PRECISION Z(1)
!     COMMON /EMGZZZ/  Z
 
!     NGG = ORDER OF ELEMENT MATRIX
 
!     INDEX      STIFFNESS             MASS                HEAT
!                AND MASS              ONLY              TRANSFER
 
!     IZ    NGG BY 3 PARTITION  NGG BY 3 PARTITION  FOUR WORD COORDINATE
!           OF MATRIX           OF MATRIX           VECTOR.  INPUT TO
!                                                   TRANSD
 
!     IZ+2                                          TRANSFORMED THERMAL
!                                                   CONDUCTANCE MATRIX
 
!     IT                                            MATERIAL TRANSFOR-
!                                                   MATION MATRIX
 
!     IK    SYMMETRIC HALF OF   SAME AS IZ          FULL CONDUCTANCE
!           STIFFNESS
 
!     IM    SYMMETRIC HALF OF   SYMMETRIC HALF OF   FULL CAPACITANCE
!           MASS                MASS
 
!     IN    --------------------SHAPE FUNCTIONS-------------------------
 
!     IG    --------------------D(SHAPE)/D(GREEK)-----------------------
 
!     IX    --------------------D(SHAPE)/D(BASIC XYZ)-------------------
 
!     ID    DISPLACEMENT
!           VECTOR IN BASIC
!           COORDINATES
 
!     CHECK GEOMETRY.  THE FOLLOWING CHECKS ARE MADE
!           1.  ASPECT RATIO
!           2.  ANGLES BETWEEN NORMALS OF SUB-TRIANGLES ON EACH FACE
!           3.  ANGLES BETWEEN VECTORS BETWEEN POINTS ALONG EACH EDGE
!           4.  REVERSE SEQUENCING
!           5.  DUPLICATE COORDINATE VALUES
 
!     FETCH EPT DATA, COMPUTE EST POINTERS
 
 100 mid  = 10 + 12*(TYPE-1)
 cid  = iest(mid+1)
 nip  = iest(mid+2)
 maxar= est(mid+3)
 alfa = est(mid+4)
 beta = est(mid+5)
 bgpdt= mid + 6
 gpt  = bgpdt + ngp*4
 mid  = iest(mid)
 IF (nip < 2 .OR. nip > 4) nip = TYPE/2 + 2
 IF (maxar <= 0.0) maxar = dmaxar(TYPE)
 IF (alfa  < 0.0) alfa  = dalfa(TYPE)
 IF (beta < 0.0 .AND. TYPE /= 1) beta = dbeta(TYPE-1)
 alfa = COS(dtor*alfa)
 beta = COS(dtor*beta)
 IF (ugv == 0) GO TO 105
 
!     TRANSFORM DISPLACEMENT VECTOR TO BASIC COORDINATES
!     MULTIPLY BY 1/4 TO AVOID MULTIPLYING STRAIN-DISPLACEMENT
!     RELATIONS BY 1/2 UNDER THE INTEGRAL.  DITTO FOR LOADING TEMP-S.
 
 dstld = gptld(1)
 DO  i = 1,ngp
   m = bgpdt + 4*i - 4
   j = ugv + sil(i) - 1
   k = id + 3*i - 3
   IF (iest(m) == 0) GO TO 102
   CALL transd (est(m),tk)
   DO  l = 1,3
     z(iz+l-1) = DBLE(rz(j+l-1)*0.25)
   END DO
   CALL gmmatd (tk,3,3,0,z(iz),3,1,0,z(n))
   gptld(i) = 0.25*gptld(i)
   CYCLE
   102 DO  l = 1,3
     z(n+l-1) = DBLE(rz(j+l-1)*0.25)
   END DO
   gptld(i) = 0.25*gptld(i)
 END DO
 
!     REARRANGE BGPDT
 
 105 DO  i = 1,ngp
   jz(izs+i) = iest(bgpdt+i*4-4)
 END DO
 bcord = gpt - 3
 DO  i = 2,ngp
   DO  j = 1,3
     k = bgpdt + 4*(ngp-i) + 4 - j
     bcord = bcord - 1
     est(bcord) = est(k)
   END DO
 END DO
 DO  i = 2,ngp
   iest(bgpdt+i-1) = jz(izs+i)
 END DO
 
!     IF COMPUTING DIFFERENTIAL STIFFNESS, SKIP CHECKS
 
 IF (ugv > 0) GO TO 500
 
!     FIND 8 POINTERS TO CORNER COORDINATES IN EST
 
!     EDGE        CORNERS
!       1         1     2
!       2         2     3
!       3         3     4
!       4         4     1
!       5         1     5
!       6         2     6
!       7         3     7
!       8         4     8
!       9         5     6
!      10         6     7
!      11         7     8
!      12         8     5
 
 nc(1) = bcord
 j = 3*TYPE
 SELECT CASE ( TYPE )
   CASE (    1)
     GO TO 140
   CASE (    2)
     GO TO 150
   CASE (    3)
     GO TO 160
 END SELECT
 140 nc(5) = bcord + 12
 GO TO 170
 150 nc(5) = bcord + 36
 GO TO 170
 160 nc(5) = bcord + 60
 170 DO  i = 2,4
   nc(i  ) = nc(i-1) + j
   nc(i+4) = nc(i+3) + j
 END DO
 
!     COMPUTE 12 EDGE VECTORS, FIND SMALLEST AND LARGEST MAGNITUDES
 
 i = 0
 j = 1
 smag = 1.0E20
 bmag = 0.0
 DO  edge = 1,12
   SELECT CASE ( edge )
     CASE (    1)
       GO TO 190
     CASE (    2)
       GO TO 190
     CASE (    3)
       GO TO 190
     CASE (    4)
       GO TO 200
     CASE (    5)
       GO TO 210
     CASE (    6)
       GO TO 190
     CASE (    7)
       GO TO 190
     CASE (    8)
       GO TO 190
     CASE (    9)
       GO TO 220
     CASE (   10)
       GO TO 190
     CASE (   11)
       GO TO 190
     CASE (   12)
       GO TO 200
   END SELECT
   190 i = i + 1
   j = j + 1
   l = nc(i) - 1
   m = nc(j) - 1
   GO TO 230
   200 l = m
   m = nc(j-3) - 1
   GO TO 230
   210 i = 0
   j = 4
   GO TO 190
   220 i = 4
   j = 5
   GO TO 190
   230 tmag = 0.0
   DO  k = 1,3
     evec(k,edge) = est(m+k) - est(l+k)
     tmag = tmag + evec(k,edge)**2
   END DO
   IF (tmag < smag) smag = tmag
   IF (tmag > bmag) bmag = tmag
 END DO
 
!     CHECK ASPECT RATIO
 
 IF (smag > 0.0) GO TO 260
 smag = 1.0E-10
 260 IF (bmag/smag <= maxar**2) GO TO 265
 WRITE (otpt,7200) ufm,nerr2,ihex,TYPE,elno,eid,bar,excd
 nogo = .true.
 
!     CHECK ANGLES BETWEEN FACE NORMALS
 
!     FACE              CORNERS
!       1         1     4     3     2
!       2         1     2     6     5
!       3         2     3     7     6
!       4         3     4     8     7
!       5         4     1     5     8
!       6         5     6     7     8
 
 265 DO  face = 1,6
   SELECT CASE ( face )
     CASE (    1)
       GO TO 270
     CASE (    2)
       GO TO 280
     CASE (    3)
       GO TO 290
     CASE (    4)
       GO TO 290
     CASE (    5)
       GO TO 300
     CASE (    6)
       GO TO 310
   END SELECT
   270 i = 1
   j = 4
   k = 3
   l = 2
   GO TO 320
   280 i = 1
   j = 6
   k = 9
   l = 5
   GO TO 320
   290 i = i + 1
   j = j + 1
   k = k + 1
   l = l + 1
   GO TO 320
   300 i = 4
   j = 5
   k = 12
   l = 8
   GO TO 320
   310 i = 12
   j = 9
   k = 10
   l = 11
   320 DO  n = 1,2
     vn(1,1) = evec(2,i)*evec(3,j) - evec(3,i)*evec(2,j)
     vn(2,1) = evec(3,i)*evec(1,j) - evec(1,i)*evec(3,j)
     vn(3,1) = evec(1,i)*evec(2,j) - evec(2,i)*evec(1,j)
     vn(1,2) = evec(2,k)*evec(3,l) - evec(3,k)*evec(2,l)
     vn(2,2) = evec(3,k)*evec(1,l) - evec(1,k)*evec(3,l)
     vn(3,2) = evec(1,k)*evec(2,l) - evec(2,k)*evec(1,l)
     smag = 0.0
     bmag = 0.0
     tmag = 0.0
     DO  m = 1,3
       smag = smag + vn(m,1)**2
       bmag = bmag + vn(m,2)**2
       tmag = vn(m,1)*vn(m,2) + tmag
     END DO
     smag = SQRT(smag*bmag)
     IF (smag == 0.0) GO TO 335
     
!     EPSILON INTRODUCED TO OVERCOME ROUNDOUT ERROR
     
     IF (tmag/smag >= 0.99*alfa) GO TO 335
     WRITE (otpt,7200) ufm,nerr2,ihex,TYPE,elno,eid,balfa,excd
     nogo = .true.
     335 m = i
     i = l
     l = k
     k = j
     j = m
   END DO
 END DO
 
!     CHECK MID-EDGE POINTS
 
 IF (TYPE == 1) GO TO 455
 m = 1
 DO  edge = 1,12
   SELECT CASE ( edge )
     CASE (    1)
       GO TO 370
     CASE (    2)
       GO TO 370
     CASE (    3)
       GO TO 370
     CASE (    4)
       GO TO 370
     CASE (    5)
       GO TO 380
     CASE (    6)
       GO TO 390
     CASE (    7)
       GO TO 390
     CASE (    8)
       GO TO 390
     CASE (    9)
       GO TO 400
     CASE (   10)
       GO TO 370
     CASE (   11)
       GO TO 370
     CASE (   12)
       GO TO 370
   END SELECT
   370 i = nc(m)
   j = i + 3
   k = j + 3
   l = k + 3
   m = m + 1
   IF (edge /= 4 .AND. edge /= 12) GO TO 410
   IF (TYPE == 2) k = nc(m-4)
   IF (TYPE == 3) l = nc(m-4)
   GO TO 410
   380 m = 0
   390 m = m + 1
   i = nc(m)
   j = i + 12*TYPE - 3*(m-1)*(TYPE-1)
   k = j + 12
   k = k + 3*(m-1)*(3-TYPE)
   l = nc(m+4)
   GO TO 410
   400 m = 5
   GO TO 370
   410 smag = 0.0
   bmag = 0.0
   tmag = 0.0
   DO  n = 1,3
     vn(n,1) = est(j+n-1) - est(i+n-1)
     vn(n,2) = est(k+n-1) - est(j+n-1)
     tmag = tmag + vn(n,1)*vn(n,2)
     smag = smag + vn(n,1)**2
     bmag = bmag + vn(n,2)**2
   END DO
   smag = SQRT(smag*bmag)
   IF (smag == 0.0) GO TO 430
   IF (tmag/smag >= beta) GO TO 430
   GO TO 445
   430 IF (TYPE == 2) CYCLE
   tmag = 0.0
   smag = 0.0
   DO  n = 1,3
     vn(n,1) = est(l+n-1) - est(k+n-1)
     tmag = tmag + vn(n,1)*vn(n,2)
     smag = smag + vn(n,1)**2
   END DO
   smag = SQRT(smag*bmag)
   IF (smag == 0.0) CYCLE
   IF (tmag/smag >= beta) CYCLE
   445 WRITE (otpt,7200) ufm,nerr2,ihex,TYPE,elno,eid,bbeta,excd
   nogo = .true.
 END DO
 
!     CHECK FOR LEFT-HANDED ELEMENT COORDINATE SYSTEM
 
!     VOL = EVEC(5)*(EVEC(1) X -EVEC(4))
 
 455 vn(1,1) = evec(2,4)*evec(3,1) - evec(3,4)*evec(2,1)
 vn(2,1) = evec(3,4)*evec(1,1) - evec(1,4)*evec(3,1)
 vn(3,1) = evec(1,4)*evec(2,1) - evec(2,4)*evec(1,1)
 tmag    = 0.0
 DO  i = 1,3
   tmag = tmag + evec(i,5)*vn(i,1)
 END DO
 IF (tmag > 0.0) GO TO 470
 WRITE (otpt,7200) ufm,nerr2,ihex,TYPE,elno,eid,rvrs
 nogo = .true.
 
!     CHECK FOR DUPLICATE COORDINATE VALUES
 
 470 l = ngp - 1
 DO  i = 1,l
   m = bcord + 3*(i-1)
   k = i + 1
   DO  j = k,ngp
     n = bcord + 3*(j-1)
     IF (est(m  ) /= est(n  )) CYCLE
     IF (est(m+1) /= est(n+1)) CYCLE
     IF (est(m+2) /= est(n+2)) CYCLE
     WRITE (otpt,7200) ufm,nerr2,ihex,TYPE,elno,eid,twins
     nogo = .true.
   END DO
 END DO
 
!     IF NOGO FLAG ON, DON T COMPUTE ELEMENT MATRICES
 
 IF (nogo) RETURN
 
!     INITIALIZE FOR NUMERICAL INTEGRATION
 
!     ABSCISSAE AND WEIGHT COEFFICIENTS FOR GAUSSIAN QUADRATURE
 
 500 i = nip - 1
 SELECT CASE ( i )
   CASE (    1)
     GO TO 510
   CASE (    2)
     GO TO 520
   CASE (    3)
     GO TO 530
 END SELECT
 510 h(1) = 1.0
 s(1) = gauss(1)
 h(2) = 1.0
 s(2) =-gauss(1)
 GO TO 540
 520 h(1) = gauss(2)
 s(1) = gauss(3)
 h(2) = gauss(4)
 s(2) = 0.0
 h(3) = gauss(2)
 s(3) =-gauss(3)
 GO TO 540
 530 h(1) = gauss(5)
 s(1) = gauss(6)
 h(2) = gauss(7)
 s(2) = gauss(8)
 h(3) = gauss(7)
 s(3) =-gauss(8)
 h(4) = gauss(5)
 s(4) =-gauss(6)
 
!     GENERATE TABLE OF EQUIVALENTS IN SIL ARRAY SO MATRIX WILL BE
!     ORDERED ACCORDING TO INCREASING SIL NUMBERS
 
 540 i = -ngp
 545 j =  0
 DO  k = 1,ngp
   IF (sil(k) < j) CYCLE
   j = sil(k)
   l = k
 END DO
 sil(l) = i
 i = i + 1
 IF (i < 0) GO TO 545
 DO  i = 1,ngp
   sil(i) = -sil(i)
 END DO
 
!     NOW SIL(I) = PARTITION NUMBER OF ELEMENT GRID POINT I
 
!     ZERO OUT OPEN CORE FOR MATRIX SUMMATION
 
 DO  i = ik,nm
   z(i) = 0.0
 END DO
 
!     BRANCH ON HEAT TRANSFER FLAG
 
 IF (heat == 1) GO TO 3000
 
!     FETCH MATERIAL PROPERTIES
 
!     =============================================================
!     THIS SECTION OF CODE MUST BE UPDATED WHEN GENERAL ANISOTROPIC
!     MATERIAL IS ADDED.
 
!     TEST FOR ANISOTROPIC MATERIAL
 
 inflag = 10
 anis   =.false.
 
!     TEST FOR RECTANGULAR COORDINATE SYSTEM IN WHICH THE ANISOTROPIC
!     MATERIAL IS DEFINED
 
 rect = .true.
!     ===============================================================
 
!     CHECK FOR TEMPERATURE DEPENDENCE
 
 tdep = .true.
 DO  i = 2,ngp
   IF (est(gpt) /= est(gpt+i-1)) GO TO 630
 END DO
 tdep = .false.
 630 temp = est(gpt)
 CALL mat (eid)
 IF (.NOT.mtdep) tdep = .false.
 IF (ib(46) == 6) anis = .true.
 IF (kgg <= 0) GO TO 1000
 
!     IF ISOTROPIC, TEMPERATURE INDEPENDENT MATERIAL, COMPUTE CONSTANTS
 
 IF (anis .OR. tdep) GO TO 1000
 IF (ib(46) /=  0) GO TO 640
 WRITE (otpt,7300) ufm,mid,eid
 nogo = .true.
 RETURN
 
!     SET UP FOR EASY MULTIPLICATION IF MATERIALS ARE ON MAT1
 
 640 e1 = bufm6(1)
 e2 = bufm6(2)
 e3 = bufm6(22)
 
!     ============================================================
!     CODE TO TRANSFORM GENERAL ANISOTROPIC MATERIAL PROPERTIES TO
!     BASIC COORDINATE SYSTEM MUST BE ADDED HERE.
!     ============================================================
 
!     ALL SET TO BEGIN INTEGRATION LOOPS.  DO IT.
 
 1000 tvol = 0.0D+0
 DO  i = 1,nip
   DO  j = 1,nip
     DO  k = 1,nip
       
!     GENERATE SHAPE FUNCTIONS AND JACOBIAN MATRIX INVERSE
       
       CALL ihexsd (TYPE,z(in),z(ig),jacob,detj,eid,s(i),s(j),s(k),  &
           est(bcord))
       IF (detj /= 0.0) GO TO 1010
       
!     BAD ELEMENT IF FALL HERE.  JACOBIAN MATRIX WAS SINGULAR.
       
       nogo = .true.
       RETURN
       
       1010 sfact = h(i)*h(j)*h(k)*detj
       tvol  = tvol + sfact
       IF (kgg <= 0) GO TO 1015
       
!     STIFFNESS
       
!     COMPUTE STRAIN-DISPLACEMENT RELATIONS
       
!     MUST REVERSE CALLING ORDER SINCE MATRICES ARE STORED BY COLUMNS
       
       CALL gmmatd (z(ig),ngp,3,0,jacob,3,3,0,z(ix))
       
!     IF MATERIAL IS TEMPERATURE DEPENDENT, MUST COMPUTE TEMPERATURE
!     AT THIS INTEGRATION POINT AND FETCH MATERIAL PROPERTIES AGAIN
       
       1015 IF (.NOT. tdep) GO TO 1030
       temp = 0.0
       DO  l = 1,ngp
         temp = temp + z(in+l-1)*est(gpt+l-1)
       END DO
       CALL mat (eid)
       IF (kgg <= 0) GO TO 1100
       IF (anis) GO TO 1040
       IF (ib(46) /= 0) GO TO 1025
       WRITE (otpt,7300) ufm,mid,eid
       nogo = .true.
       RETURN
       
       1025 e1 = bufm6(1)
       e2 = bufm6(2)
       e3 = bufm6(22)
       GO TO 1100
       1030 IF (kgg <= 0) GO TO 1100
       
!     IF MATERIAL IS ANISOTROPIC AND NOT DEFINED IN RECTANGULAR COOR-
!     DINATE SYSTEM, MUST TRANSFORM TO BASIC COORDINATE SYSTEM AT THIS
!     INTEGRATION POINT
!     IN THIS VERSION, ANISOTROPIC MATERIAL SYSTEMS MUST BE RECTANGULAR.
!     THEREFORE, NO FURTHER TRANSFORMATIONS ARE NECESSARY
       
       
!     ================================================================
!     THIS CODE MUST BE COMPLETED WHEN GENERAL ANISOTROPIC MATERIAL IS
!     ADDED
       
       IF (.NOT.anis) GO TO 1100
       1040 CONTINUE
       
!     INSERT GLOBAL TO BASIC TRANSFORMATION OPERATIONS HERE FOR
!     ANISOTROPIC MATERIAL MATRIX
!     =============+==================================================
       
       DO  ijk = 1,36
         gmat(ijk) = bufm6(ijk)
       END DO
       IF (rect) GO TO 1100
       
!     MATERIAL HAS BEEN EVALUATED FOR THIS INTEGRATION POINT WHEN
!     FALL HERE.
       
       1100 IF (ugv == 0) GO TO 1170
       
!     COMPUTE STRESSES FOR DIFFERENTIAL STIFFNESS MATRIX
       
!     THERMAL EFFECTS
       
       IF (idstld == -1) GO TO 1120
       temp = 0.0
       DO  l = 1,ngp
         temp = temp + z(in+l-1)*DBLE(gptld(l))
       END DO
       temp = temp - DBLE(tref)
       IF (anis) GO TO 1115
       sig(1) =-DBLE(talpha)*(e1+2.0*e2)*temp
       sig(2) = sig(1)
       sig(3) = sig(1)
       sig(4) = 0.0
       sig(5) = 0.0
       sig(6) = 0.0
       GO TO 1140
!     ===========================================================
       1115 CONTINUE
       
!     ADD THERMAL STRESS COMPUTATIONS FOR ANISOTROPIC MATERIAL
       
!     STORE ALPHA IN DOUBLE PRECISION
       
       DO  ijk = 1,6
         dalpha(ijk) = bufm6(ijk+37)
       END DO
       
       CALL gmmatd (gmat,6,6,0, dalpha,6,1,0,sig)
       DO  ijk = 1,6
         sig(ijk) = -sig(ijk)*temp
       END DO
       GO TO 1140
!     ===========================================================
       1120 DO  l = 1,6
         sig(l) = 0.0
       END DO
       
!     DISPLACEMENT EFFECTS, COMPUTE STRESS MATRIX AND MULTIPLY BY DISPL.
       
       1140 str(12) = 0.0
       str(13) = 0.0
       str(17) = 0.0
       DO  l = 1,ngp
         ii = ix + 3*l - 4
         IF (anis) GO TO 1145
         str( 1) = e1*z(ii+1)
         str( 2) = e2*z(ii+2)
         str( 3) = e2*z(ii+3)
         str( 4) = e2*z(ii+1)
         str( 5) = e1*z(ii+2)
         str( 6) = e2*z(ii+3)
         str( 7) = e2*z(ii+1)
         str( 8) = e2*z(ii+2)
         str( 9) = e1*z(ii+3)
         str(10) = e3*z(ii+2)
         str(11) = e3*z(ii+1)
         str(14) = e3*z(ii+3)
         str(15) = e3*z(ii+2)
         str(16) = e3*z(ii+3)
         str(18) = e3*z(ii+1)
         GO TO 1150
!     =========================================================
         
         1145 CONTINUE
         
!     ADD STRESS MATRIX COMPUTATION FOR ANISOTROPIC MATERIAL
         
         DO  ijk = 1,18
           store(ijk) = 0.d0
         END DO
         store( 1) = z(ii+1)
         store( 5) = z(ii+2)
         store( 9) = z(ii+3)
         store(10) = z(ii+2)
         store(11) = z(ii+1)
         store(14) = z(ii+3)
         store(15) = z(ii+2)
         store(16) = z(ii+3)
         store(18) = z(ii+1)
         
         CALL gmmatd (gmat,6,6,0,store(1),6,3,0,str)
         
!     ============================================================
         
         1150 CALL gmmatd (str,6,3,-2,z(id+3*l-3),3,1,0,sig)
       END DO
       str(1) = sx
       sx = sx + sy
       sy = sy + sz
       sz = sz + str(1)
       
!     NOW BEGIN LOOPS OVER GRID POINTS ALONG ROWS AND COLUMNS
       
       1170 DO  n = 1,ngp
         DO  m = n,ngp
           
!     COMPUTE PARTITION FOR POINTWISE ROW M AND COLUMN N
           
           IF (kgg <= 0) GO TO 1300
           IF (.NOT.anis ) GO TO 1200
           
!     =================================================================
!     MUST ADD CODE TO COMPUTE THE CONTRIBUTION TO THE STIFFNESS MATRIX
!     FOR ANISOTROPIC MATERIAL HERE
!     =================================================================
           
           1200 IF (sil(m) >= sil(n)) GO TO 1210
           
!     MUST COMPUTE TRANSPOSE OF THIS PARTITION FOR SUMMATION IN ELEMENT
!     MATRIX
           
           mz = ix + (n-1)*3
           nz = ix + (m-1)*3
           GO TO 1220
           1210 mz = ix + (m-1)*3
           nz = ix + (n-1)*3
           1220 IF (ugv == 0) GO TO 1222
           
!     DIFFERENTIAL STIFFNESS
           
           DO  l = 1,3
             DO  inc = 1,3
               c(l,inc) = z(mz+inc-1)*z(nz+l-1)
             END DO
           END DO
           part(1,1) = sx*c(2,2) + syz*(c(2,3)+c(3,2)) + sz*c(3,3)
           part(2,2) = sy*c(3,3) + szx*(c(3,1)+c(1,3)) + sx*c(1,1)
           part(3,3) = sz*c(1,1) + sxy*(c(1,2)+c(2,1)) + sy*c(2,2)
           part(2,1) =-sx*c(2,1) + sxy*c(3,3) -syz*c(1,3) - szx*c(2,3)
           part(3,1) =-sz*c(3,1) - sxy*c(3,2) -syz*c(2,1) + szx*c(2,2)
           part(1,2) =-sx*c(1,2) + sxy*c(3,3) -syz*c(3,1) - szx*c(3,2)
           part(3,2) =-sy*c(3,2) - sxy*c(3,1) +syz*c(1,1) - szx*c(1,2)
           part(1,3) =-sz*c(1,3) - sxy*c(2,3) -syz*c(1,2) + szx*c(2,2)
           part(2,3) =-sy*c(2,3) - sxy*c(1,3) +syz*c(1,1) - szx*c(2,1)
           GO TO 1228
           
!     ELASTIC STIFFNESS
           
           1222 IF (.NOT.anis) GO TO 1226
           
!     STORE CI MATRIX
           
           DO  ijk = 1,18
             store(ijk) = 0.d0
           END DO
           store( 1) = z(mz  )
           store( 4) = z(mz+1)
           store( 6) = z(mz+2)
           store( 8) = z(mz+1)
           store(10) = z(mz  )
           store(11) = z(mz+2)
           store(15) = z(mz+2)
           store(17) = z(mz+1)
           store(18) = z(mz  )
           
           CALL gmmatd (store(1),3,6,0,gmat(1),6,6,0,store(19))
           
!     STORE CJ
           
           DO  ijk = 1,18
             store(ijk) = 0.d0
           END DO
           store( 1) = z(nz  )
           store( 5) = z(nz+1)
           store( 9) = z(nz+2)
           store(10) = z(nz+1)
           store(11) = z(nz  )
           store(14) = z(nz+2)
           store(15) = z(nz+1)
           store(16) = z(nz+2)
           store(18) = z(nz  )
           
           CALL gmmatd (store(19),3,6,0,store(1),6,3,0,store(37))
           ijkl = 0
           DO  ijk = 1,3
             DO  ijl = 1,3
               ijkl = ijkl + 1
               part(ijk,ijl) = store(ijkl+36)
             END DO
           END DO
           GO TO 1228
           1226 part(1,1) = e1*z(nz)*z(mz) + e3*(z(nz+1)*z(mz+1) +z(nz+2)*z(mz+2))
           part(2,2) = e1*z(nz+1)*z(mz+1) + e3*(z(nz)*z(mz) +z(nz+2)*z(mz+2))
           part(3,3) = e1*z(nz+2)*z(mz+2) + e3*(z(nz)*z(mz) +z(nz+1)*z(mz+1))
           part(2,1) = e2*z(nz  )*z(mz+1) + e3*z(nz+1)*z(mz  )
           part(3,1) = e2*z(nz  )*z(mz+2) + e3*z(nz+2)*z(mz  )
           part(1,2) = e2*z(nz+1)*z(mz  ) + e3*z(nz  )*z(mz+1)
           part(3,2) = e2*z(nz+1)*z(mz+2) + e3*z(nz+2)*z(mz+1)
           part(1,3) = e2*z(nz+2)*z(mz  ) + e3*z(nz  )*z(mz+2)
           part(2,3) = e2*z(nz+2)*z(mz+1) + e3*z(nz+1)*z(mz+2)
           
!     ADD STIFFNESS PARTITION TO ELEMENT MATRIX
           
!     COMPUTE INDEX INTO OPEN CORE WHERE PART(1,1) IS TO BE ADDED.
           
           1228 IF (sil(m)-sil(n) < 0.0) THEN
             GO TO  1230
           ELSE IF (sil(m)-sil(n) == 0.0) THEN
             GO TO  1240
           ELSE
             GO TO  1250
           END IF
           1230 mz = sil(n)
           nz = sil(m)
           diag = .false.
           GO TO 1260
           1240 mz = sil(m)
           nz = sil(n)
           diag = .true.
           GO TO 1260
           1250 mz = sil(m)
           nz = sil(n)
           diag = .false.
           
!     COLUMN NUMBER
           
           1260 l = (nz-1)*3 + 1
           
!     INCREMENT BETWEEN COLUMNS
           
           inc = ngg - l
           
!     FIRST WORD OF COLUMN
           
           l = ik + ((l-1)*l)/2 + (inc+1)*(l-1)
           
!     WORD IN COLUMN FOR THIS ROW
           
           l = l + 3*(mz-nz)
           
!     ADD PARTITION
           
           DO  nz = 1,3
             DO  mz = 1,3
               IF (diag .AND. mz < nz) CYCLE
               z(l+mz-1) = z(l+mz-1) + part(mz,nz)*sfact
             END DO
             l = l + inc
             inc = inc - 1
           END DO
           1300 IF (mgg <= 0) CYCLE
           
!     MASS
           
!     COMPUTE TERM FOR MASS MATRIX
           
           rho = bufm6(37)
           mz  = sil(m)
           nz  = sil(n)
           IF (mz >= nz) GO TO 1310
           mz  = sil(n)
           nz  = sil(m)
           
!     COMPUTE INDEX INTO OPEN CORE FOR THIS MASS TERM
           
           1310 l = (nz*(nz+1))/2 + (nz-1)*(ngp-nz) + mz - nz + im - 1
           
!     COMPUTE AND ADD MASS TERM TO ELEMENT MATRIX
           
           z(l) = z(l) + DBLE(rho)*sfact*z(in+m-1)*z(in+n-1)
         END DO
       END DO
     END DO
   END DO
 END DO
 
!     END OF INTEGRATION LOOPS
 
 icode = 7
 
!     LOOK FOR NON-BASIC COORDINATE SYSTEM
 
 nocstm = .false.
 DO  i = 1,ngp
   IF (iest(bgpdt+i-1) /= 0) GO TO 2005
 END DO
 nocstm = .true.
 GO TO 2065
 
!     RESTORE GRID POINT DATA TO ORIGINAL FORM FOR DOING TRANSFORM
!     TO GLOBAL COORDINATES
 
!     FIRST, TRANSFER IT TO OPEN CORE AT IN
 
 2005 k = (in-1)*2 + 1
 j = ngp*4
 DO  i = 1,j
   rz(k+i-1) = est(bgpdt+i-1)
 END DO
 
!     NOW MOVE IT BACK AND REARRANGE IT
 
 DO  i = 1,ngp
   iest(bgpdt+4*i-4) = jz(k+i-1)
   DO  j = 1,3
     est(bgpdt+4*i-4+j) = rz(k+ngp+3*i+j-4)
   END DO
 END DO
 
!     FETCH GLOBAL TO BASIC TRANSFORMATION MATRICES
 
 DO  i = 1,ngp
   j = in + (i-1)*9
   CALL transd (est(bgpdt+4*i-4),z(j))
 END DO
 IF (kgg <= 0) GO TO 2110
 
!     TRANSFORM STIFFNESS TO GLOBAL COORDINATES
 
 i = 0
 2026 i = i + 1
 icp = sil(i)
 
!     COLUMN INDICES
 
 k   = (icp-1)*3 + 1
 inc = ngg - k + 1
 l   = ik + ((k-1)*k)/2 + inc*(k-1)
 m   = l + inc
 n   = m + inc - 1
 
!     TRANSFORMATION MATRIX INDEX
 
 igcs = iest(bgpdt+4*i-4)
 nz = in + (i-1)*9
 IF (igcs == 0) GO TO 2028
 
!     TERMS ON DIAGONAL PARTITION
 
 ASSIGN 2028 TO back
 GO TO 6000
 
!     OFF-DIAGONAL PARTITIONS
 
 2028 l = l + 3
 m = m + 2
 n = n + 1
 irp = icp + 1
 IF (irp > ngp) GO TO 2060
 mz = nz
 DO  j = irp,ngp
   DO  k = 1,ngp
     IF (j == sil(k)) EXIT
   END DO
   2031 IF (igcs /= 0) GO TO 2032
   IF (iest(bgpdt+4*k-4) == 0) GO TO 2045
   2032 nz = in + (k-1)*9
   DO  k = 1,3
     tk(k,1) = 0.0
     tk(k,2) = 0.0
     tk(k,3) = 0.0
     DO  ii = 1,3
       tk(k,1) = tk(k,1) + z(l+ii-1)*z(nz+3*ii+k-4)
       tk(k,2) = tk(k,2) + z(m+ii-1)*z(nz+3*ii+k-4)
       tk(k,3) = tk(k,3) + z(n+ii-1)*z(nz+3*ii+k-4)
     END DO
   END DO
   DO  k = 1,3
     z(l+k-1) = 0.0
     z(m+k-1) = 0.0
     z(n+k-1) = 0.0
     DO  ii = 1,3
       z(l+k-1) = z(l+k-1) + tk(k,ii)*z(mz+3*ii-3)
       z(m+k-1) = z(m+k-1) + tk(k,ii)*z(mz+3*ii-2)
       z(n+k-1) = z(n+k-1) + tk(k,ii)*z(mz+3*ii-1)
     END DO
   END DO
   2045 l = l + 3
   m = m + 3
   n = n + 3
 END DO
 2060 IF (i < ngp) GO TO 2026
 
!     BUILD STIFFNESS PARTITIONS AND PASS TO EMGOUT
 
 2065 idon = 0
 DO  i = 1,ngp
   IF (i == ngp) idon = 1
   DO  j = 1,3
     
!     COLUMN NUMBER
     
     k = (i-1)*3 + j
     
!     NUMBER OF TERMS TO FETCH TO COMPLETE THIS COLUMN IN PARTITION
     
     l = k - 1
     IF (l == 0) GO TO 2075
     
!     FETCH TERMS AND LOAD INTO J-TH COLUMN OF PARTITION
     
     n   = ik + l
     inc = ngg - 1
     DO  m = 1,l
       z(iz+ngg*j-ngg+m-1) = z(n)
       n   = n + inc
       inc = inc - 1
     END DO
     
!     FILL OUT PARTITION WITH COLUMNS OF STIFFNESS MATRIX
     
!     COMPUTE INDEX IN OPEN CORE OF FIRST TERM OF COLUMN K
     
     2075 n = ik + ((k-1)*k)/2 + (ngg-k+1)*(k-1)
     
!     INSERT THIS COLUMN IN PARTITION
     
     DO  m = k,ngg
       z(iz+ngg*j-ngg+m-1) = z(n)
       n = n + 1
     END DO
   END DO
   dict(5) = ib(45)
   CALL emgout (z(iz),z(iz),3*ngg,idon,dict,1,2)
 END DO
 
!     EXPAND AND TRANSFORM MASS MATRIX AND PASS TO EMGOUT
 
 IF (mgg <= 0) GO TO 2400
 2110 idon = 0
 DO  i = 1,ngp
   IF (i == ngp) idon = 1
   DO  j = 1,ngp
     
!     COMPUTE INDEX INTO OPEN CORE FOR MASS TERM
     
     k = i
     l = j
     IF (i <= j) GO TO 2115
     k = j
     l = i
     2115 n = ((k-1)*k)/2 + (k-1)*(ngp-k+1) + l - k + im
     
!     MULTIPLY GLOBAL TO BASIC TRANSFORMATIONS
     
     m = iz - ngg+3*j - 4
     IF (i == j .OR. nocstm) GO TO 2116
     IF (iest(bgpdt+4*i-4) /= 0) GO TO 2118
     IF (iest(bgpdt+4*j-4) /= 0) GO TO 2118
     2116 z(m+ngg  +1) = z(n)
     z(m+ngg  +2) = 0.0
     z(m+ngg  +3) = 0.0
     z(m+ngg*2+1) = 0.0
     z(m+ngg*2+2) = z(n)
     z(m+ngg*2+3) = 0.0
     z(m+ngg*3+1) = 0.0
     z(m+ngg*3+2) = 0.0
     z(m+ngg*3+3) = z(n)
     CYCLE
     2118 DO  k = 1,ngp
       IF (i == sil(k)) mz = in + 9*(k-1)
       IF (j == sil(k)) nz = in + 9*(k-1)
     END DO
     CALL gmmatd (z(mz),3,3,1,z(nz),3,3,0,tf)
     
!     MULTIPLY BY MASS SCALAR FOR THIS 3 BY 3 PARTITION AND STORE
!     IN NGG BY 3 PARTITION
     
     DO  k = 1,3
       DO  l = 1,3
         z(m+ngg*l+k) = tf(k,l)*z(n)
       END DO
     END DO
   END DO
   dict(5) = 0
   CALL emgout (z(iz),z(iz),3*ngg,idon,dict,2,2)
 END DO
 
!     SAVE ELEMENT BCD NAME, ID, VOLUME, MASS, NO. OF GRID POINTS, AND
!     GRID POINT DATA IN SCR4 IF USER REQUESTED VOLUME/AREA PRINTOUT
!     (NOTE - MAKE SURE THE GRID POINT DATA, BGPDT, IS IN ISTS ORIGIANL
!      FORM)
 
 2400 IF (volume <= 0.0 .AND. surfac <= 0.0) GO TO 5000
 il = iz*2
 rz(il+1) = bcd1
 rz(il+2) = bcd2(TYPE)
 jz(il+3) = eid
 rz(il+4) = tvol*volume
 rz(il+5) = tvol
 IF (rho > 0.0) rz(il+5) = tvol*rho
 jz(il+6) = ngp
 k = il + 6
 DO  i = 1,ngp
   k = k + 1
   rz(k) = est(1+i)
 END DO
 IF (surfac <= 0.0) GO TO 2460
 IF (.NOT.nocstm) GO TO 2440
 l = bgpdt + ngp
 DO  i = 1,ngp
   k = k + 1
   jz(k) = iest(bgpdt+i-1)
   DO  j = 1,3
     k = k + 1
     rz(k) = est(l)
     l = l + 1
   END DO
 END DO
 GO TO 2460
 2440 j = ngp*4
 DO  i = 1,j
   k = k + 1
   rz(k) = est(bgpdt+i-1)
 END DO
 2460 l = k - il
 CALL WRITE (scr4,rz(il+1),l,1)
 GO TO 5000
 
!     HEAT TRANSFER SECTION
 
 3000 inflag = 3
 CALL hmat (eid)
 anis = .false.
 IF (kgg <= 0) GO TO 3100
 
!     CHECK FOR ANISOTROPY
 
 IF (kheat(1) /= kheat(4) .OR. kheat(1) /= kheat(6)) GO TO 3010
 IF (kheat(2) /= 0.0 .OR. kheat(3) /= 0.0 .OR. kheat(5) /= 0.0) GO TO 3010
 GO TO 3100
 3010 anis = .true.
 it = iz + 8
 
!     CHECK FOR RECTANGULAR COORDINATE SYSTEM FOR MATERIAL
 
 rect = .true.
 IF (cid == 0) GO TO 3100
 jz(izs) = cid
 DO  i = 1,3
   rz(izs+i) = est(bcord+i-1)
 END DO
 CALL transd (rz(izs),z(it))
 DO  i = 1,3
   rz(izs+i) = -rz(izs+i)
 END DO
 CALL transd (rz(izs),z(in))
 DO  i = 1,9
   IF (z(it+i-1) /= z(in+i-1)) rect = .false.
 END DO
 
!     IF NOT DEFINED IN A RECTANGULAR SYSTEM, MUST TRANSFORM INSIDE
!     INTEGRATION LOOPS
 
 IF (.NOT.rect) GO TO 3100
 
!     TRANSFORM MATERIAL MATRIX TO BASIC SYSTEM
 
 DO  i = 1,6
   z(iz+i+1) = DBLE(kheat(i))
 END DO
 l  = iz + 2
 m  = l + 3
 n  = m + 2
 nz = it
 ASSIGN 3100 TO back
 GO TO 6000
 
!     ANISOTROPIC CONDUCTIVITY MATERIAL MATRIX NOW STORED AT RZ(IZ+2)
!     TO RZ(IZ+7)
 
!     ALL SET FOR DOING INTEGRATION.  DO IT.
 
 3100 i = 0
 3101 i = i + 1
 j = 0
 3102 j = j + 1
 k = 0
 3103 k = k + 1
 
!     GENERATE SHAPE FUNCTIONS AND JACOBIAN MATRIX INVERSE
 
 CALL ihexsd (TYPE,z(in),z(ig),jacob,detj,eid,s(i),s(j),s(k), est(bcord))
 IF (detj /= 0.0) GO TO 3110
 
!     FALL HERE IF JACOBIAN MATRIX WAS SINGULAR
 
 nogo = .true.
 RETURN
 
 3110 sfact = h(i)*h(j)*h(k)*detj
 IF (kgg <= 0) GO TO 3120
 
!     COMPUTE DERIVATIVES OF SHAPE FUNCTION W.R.T. BASIC SYSTEM.
 
!     MUST REVERSE CALLING ORDER SINCE MATRICES ARE STORED BY COLUMNS
 
 CALL gmmatd (z(ig),ngp,3,0,jacob,3,3,0,z(ix))
 
!     IF MATERIAL IS ANISOTROPIC AND NOT DEFINED IN A RECTANGULAR
!     CORDINATE SYSTEM, MUST TRANSFORM TO BASIC SYSTEM AT THIS
!     INTEGRATION POINT
 
 3120 IF (.NOT.anis) GO TO 3160
 IF (rect) GO TO 3160
 
!     COMPUTE BASIC COORDINATES VECTOR AT THIS POINT
 
 DO  l = 1,3
   rz(izs+l) = 0.0
 END DO
 DO  l = 1,ngp
   DO  m = 1,3
     rz(izs+m) = rz(izs+m) + z(in+l-1)*est(bcord + 3*l+m-4)
   END DO
 END DO
 
!     FETCH TRANSFORMATION AND CONDUCTIVITY MATRICES AND PERFORM
!     TRANSFORMATION OPERATIONS
 
 CALL transd (rz(izs),z(it))
 DO  l = 1,6
   z(iz+l+1) = DBLE(kheat(l))
 END DO
 nz = it
 l  = iz + 2
 m  = l + 3
 n  = m + 2
 ASSIGN 3160 TO back
 GO TO 6000
 
!     MATERIAL HAS BEEN EVALUATED FOR THIS INTEGRATION POINT WHEN
!     FALL HERE
 
!     NOW BEGIN LOOPS OVER GRID POINTS ALONG ROWS AND COLUMNS
 
 3160 DO  n = 1,ngp
   DO  m = n,ngp
     
!     COMPUTE 1 BY 1 PARTITION FOR ROW M AND COLUMN N
     
     IF (kgg <= 0) GO TO 3210
     
!     CONDUCTIVITY
     
     IF (anis) GO TO 3180
     
!     ISOTROPIC CASE
     
     prt1 = 0.0
     DO  l = 1,3
       prt1 = prt1 + z(ix+3*m+l-4)*z(ix+3*n+l-4)
     END DO
     prt1 = sfact*DBLE(kheat(1))*prt1
     GO TO 3190
     
!     ANISOTROPIC CASE
     
     3180 l  = ix + 3*(m-1)
     e1 = z(l)*z(iz+2) + z(l+1)*z(iz+3) + z(l+2)*z(iz+4)
     e2 = z(l)*z(iz+3) + z(l+1)*z(iz+5) + z(l+2)*z(iz+6)
     e3 = z(l)*z(iz+4) + z(l+1)*z(iz+6) + z(l+2)*z(iz+7)
     l  = ix + 3*(n-1)
     prt1 = sfact*(z(l)*e1 + z(l+1)*e2 + z(l+2)*e3)
     
!     COMPUTE INDEX INTO OPEN CORE FOR THIS TERM
     
     3190 l  = sil(m)
     mz = sil(n)
     IF (l <= mz) GO TO 3200
     l  = mz
     mz = sil(m)
     3200 l  = (l-1)*ngg + mz + ik - 1
     
!     ADD TERM TO MATRIX
     
     z(l) = z(l) + prt1
     
!     CAPACITANCE
     
     3210 IF (mgg <= 0) CYCLE
     
!     COMPUTE INDEX INTO OPEN CORE FOR THIS TERM
     
     l  = sil(m)
     mz = sil(n)
     IF (l <= mz) GO TO 3215
     l  = mz
     mz = sil(m)
     3215 l  = (l-1)*ngg + mz + im - 1
     
!     COMPUTE AND ADD TERM
     
     z(l) = z(l) + sfact*DBLE(cp)*z(in+m-1)*z(in+n-1)
   END DO
 END DO
 IF (k < nip) GO TO 3103
 IF (j < nip) GO TO 3102
 IF (i < nip) GO TO 3101
 
!     END OF HEAT TRANSFER INTEGRATION LOOPS
 
 icode = 1
 
!     FILL IN THE UPPER TRIANGLES OF THE MATRICES
 
 IF (kgg <= 0) GO TO 4010
 mz = ik
 GO TO 4020
 4010 IF (mgg <= 0) GO TO 4040
 mz = im
 4020 l  = ngg - 1
 DO  i = 1,l
   j  = i + 1
   DO  k = j,ngg
     m  = (i-1)*ngg + k + mz - 1
     n  = (k-1)*ngg + i + mz - 1
     z(n) = z(m)
   END DO
 END DO
 IF (mz == ik) GO TO 4010
 
!     PASS MATRICES TO EMGOUT
 
 4040 k = ngg**2
 dict(5) = 0
 IF (kgg > 0) CALL emgout (z(ik),z(ik),k,1,dict,1,2)
 IF (mgg > 0) CALL emgout (z(im),z(im),k,1,dict,3,2)
 
!     ALL DONE, NO ERRORS
 
 5000 RETURN
 
 
!     INTERNAL SUBROUTINE
 
!     TRANSFORM COORDINATE SYSTEM OF SYMMETRIC HALF OF A 3 BY 3 MATRIX
 
 6000 tk(1,1) = z(nz  )*z(l  )  + z(nz+3)*z(l+1)  + z(nz+6)*z(l+2)
 tk(2,1) = z(nz+1)*z(l  )  + z(nz+4)*z(l+1)  + z(nz+7)*z(l+2)
 tk(3,1) = z(nz+2)*z(l  )  + z(nz+5)*z(l+1)  + z(nz+8)*z(l+2)
 tk(1,2) = z(nz  )*z(l+1)  + z(nz+3)*z(m  )  + z(nz+6)*z(m+1)
 tk(2,2) = z(nz+1)*z(l+1)  + z(nz+4)*z(m  )  + z(nz+7)*z(m+1)
 tk(3,2) = z(nz+2)*z(l+1)  + z(nz+5)*z(m  )  + z(nz+8)*z(m+1)
 tk(1,3) = z(nz  )*z(l+2)  + z(nz+3)*z(m+1)  + z(nz+6)*z(n  )
 tk(2,3) = z(nz+1)*z(l+2)  + z(nz+4)*z(m+1)  + z(nz+7)*z(n  )
 tk(3,3) = z(nz+2)*z(l+2)  + z(nz+5)*z(m+1)  + z(nz+8)*z(n  )
 z(l  )  = z(nz  )*tk(1,1) + z(nz+3)*tk(1,2) + z(nz+6)*tk(1,3)
 z(l+1)  = z(nz  )*tk(2,1) + z(nz+3)*tk(2,2) + z(nz+6)*tk(2,3)
 z(l+2)  = z(nz  )*tk(3,1) + z(nz+3)*tk(3,2) + z(nz+6)*tk(3,3)
 z(m  )  = z(nz+1)*tk(2,1) + z(nz+4)*tk(2,2) + z(nz+7)*tk(2,3)
 z(m+1)  = z(nz+1)*tk(3,1) + z(nz+4)*tk(3,2) + z(nz+7)*tk(3,3)
 z(n  )  = z(nz+2)*tk(3,1) + z(nz+5)*tk(3,2) + z(nz+8)*tk(3,3)
 GO TO back, (2028,3100,3160)
 
 7100 FORMAT (a23,i5,2H, ,a4,i1,3A4,i9,' INSUFFICIENT CORE TO COMPUTE',  &
     ' ELEMENT MATRIX')
 7200 FORMAT (a23,i5,2H, ,a4,i1,3A4,i9,3X,18HILLEGAL geometry, ,9A4)
 7300 FORMAT (a23,' 4005. AN ILLEGAL VALUE OF -NU- HAS BEEN SPECIFIED ',  &
     'UNDER MATERIAL ID =',i10,17H for element id =,i10)
 
END SUBROUTINE ihexd
