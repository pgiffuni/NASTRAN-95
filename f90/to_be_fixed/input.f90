SUBROUTINE INPUT
     
!     INPUT I1,I2,I3,I4,I5/O1,O2,O3,O4,O5/C,N,-V1-/C,N,-V2-/C,N,-V3- $
 
 
 EXTERNAL        orf
 LOGICAL :: inopen(5)
 INTEGER :: filin(5),filout(5),FILE,fil,hfil(3,5),sperlk,  &
     modcom(9),rdflg,r,r1,r2,parama,paramb,paramc,  &
     param1,paramn,t(7),mnam(2),eee(3),ix(1),two,orf,  &
     k(100),kt(20),i1t(20),j1t(20),i2t(20),j2t(20),  &
     cord2c(2),kdn(2,20),kl(20),kno(20),kdsort(270), k1(100),k2(100),k3(70)
 REAL :: lambda,qk(100)
 CHARACTER (LEN=8) :: e1,e2,e3,chr
 CHARACTER (LEN=25) :: sfm
 CHARACTER (LEN=29) :: uim
 CHARACTER (LEN=25) :: uwm
 CHARACTER (LEN=23) :: ufm
 COMMON /xmssg / ufm,uwm,uim,sfm
 COMMON /machin/ mach
 COMMON /system/ ksystm(100)
 COMMON /BLANK / parama,paramb,paramc
 COMMON /two   / two(32)
 COMMON /zzzzzz/ x(1)
 COMMON /condas/ pi,twopi,radeg,degra,s4pisq
 EQUIVALENCE     (ksystm( 1),nbuf), (ksystm( 2),nout  ),  &
     (ksystm( 4),nin ), (ksystm(12),nlines),  &
     (ksystm(57),modcom(1)), (ksystm(95),sperlk), (qk(1),k(1)), (x(1),ix(1)),  &
     (kdsort(  1),k1(1)), (kdsort(101),k2(1)), (kdsort(201),k3(1))
 
 DATA    cord2c / 2001,20        /
 DATA    mnam   / 4HINPU, 4HT    /
 DATA    filin  / 101,102,103,104,105/
 DATA    filout / 201,202,203,204,205/
 DATA    eee    / 3*2147483647 /
 DATA    param1 , paramn  / 1, 8 /
DATA    e1,e2  / 'ENDDATA ', 'END DATA' /, e3 / 'ENDATA  ' /
DATA    inopen / 5*.false./

!     SORTSEQUENCE (INTERNALSEQUENCEID)

!                1    2    3    4    5    6    7    8    9    0

DATA  k1/ 116, 115,   2, 211,  58,  57,  59,  61,  60,  62  &
    , 169, 215, 216, 221, 144, 214, 137, 104, 134, 105  &
    , 135, 106, 136, 165, 213, 114, 233, 113, 181, 189  &
    , 191,   3, 185, 186, 184, 188, 187, 177, 178, 176  &
    , 172, 182, 190, 170, 151, 161,  56,  70,  83,  85  &
    ,   4,  78,  79,  77,  82,  81,  68,  69,  67,  63  &
    ,  71,  84,  54,  55,  49,  50,  51,  52,  23,  24  &
    ,  25,  26,  36,  37,  38,  39, 122, 123,  80,  76  &
    ,  89, 148, 138, 121, 101,  98,  99,   5,   1, 127  &
    , 128, 145, 227, 228, 229, 230, 231, 235,   6,   7 /

!                1    2    3    4    5    6    7    8    9    0

DATA  k2/   8, 129,   9,  75, 219,  10,  27,  28,  29,  30  &
    ,  31,  32,  33,  34,  35, 152, 153, 154,  92,  94  &
    , 183, 124,  91, 102, 110, 109, 140, 141, 142, 143  &
    , 205, 206, 223, 224, 210, 240, 241, 242, 243, 225  &
    , 226, 244,  96,  13, 203,  22, 150, 217, 139, 146  &
    , 220, 171, 209, 179, 234, 107, 133, 100, 155, 156  &
    , 157, 222, 158, 159, 160, 111, 245, 247, 249, 251  &
    , 253, 255, 257, 259, 261, 246, 248, 250,  16,  21  &
    , 149,  88,  90,  95, 164, 252, 254, 256, 130, 202  &
    , 232, 258, 260, 262, 199, 200, 201, 166, 167,  97 /

!                1    2    3    4    5    6    7    8    9    0

DATA  k3/ 236, 237, 238, 239, 117, 112, 204, 180,  40,  41  &
    ,  42,  14,  17, 108,  11,  12,  74,  86,  44,  45  &
    ,  93, 103,  15,  18,  19,  20,  72,  73, 118, 119  &
    , 212,  43, 194, 125, 126, 162, 131, 132, 192, 193  &
    , 195, 196, 197, 198, 207, 208, 120, 147,  64, 173  &
    ,  46,  47,  48, 163, 168, 218,  87,  53,  65, 174  &
    ,  66, 175,   0,   0,   0,   0,   0,   0,   0,   0 /

!                1    2    3    4    5    6    7    8    9    0

DATA  kl/ 901,1301,5501,5481,4501,2408, 501,5301,2801,3301,  &
    5401,5551,3001,5001,5008,5561,2001,   0,   0,   0/
DATA  kt/   2,   2,   1,   1,   1,   1,   2,   1,   1,   1,  &
    1,   1,   2,   1,   1,   1,   1,   0,   0,   0/

DATA i1t/   2,   2,   5,   5,   4,   3,   2,   5,   3,   4,  &
    5,   5,   3,   5,   5,   6,   3,   0,   0,   0/
DATA j1t/   9,  13,   7,  10,  13,   8,   5,   5,  12,   1,  &
    6,   1,  14,   2,   2,  12,   4,   0,   0,   0/

DATA i2t/   7,   7,   0,   0,   0,   0,   7,   0,   0,   0,  &
    0,   0,   7,   0,   0,   0,   0,   0,   0,   0/
DATA j2t/   8,  10,   0,   0,   0,   0,   9,   0,   0,   0,  &
    0,   0,  16,   0,   0,   0,   0,   0,   0,   0/

DATA kno/  76,  68,  16,  12,   1, 180,  72,   4,  57,  52,  &
    25, 105,  48,  15, 258, 215,   9,   0,   0,   0/

DATA       kdn(1, 1),kdn(2, 1) / 4HCELA , 4HS4   /,  &
    kdn(1, 2),kdn(2, 2) / 4HCMAS , 4HS4   /,  &
    kdn(1, 3),kdn(2, 3) / 4HSPC  , 4H     /,  &
    kdn(1, 4),kdn(2, 4) / 4HSPC1 , 4H     /,  &
    kdn(1, 5),kdn(2, 5) / 4HGRID , 4H     /,  &
    kdn(1, 6),kdn(2, 6) / 4HCBAR , 4H     /,  &
    kdn(1, 7),kdn(2, 7) / 4HCDAM , 4HP4   /,  &
    kdn(1, 8),kdn(2, 8) / 4HSEQG , 4HP    /,  &
    kdn(1, 9),kdn(2, 9) / 4HCQUA , 4HD1   /,  &
    kdn(1,10),kdn(2,10) / 4HCTRI , 4HA1   /,  &
    kdn(1,11),kdn(2,11) / 4HSLOA , 4HD    /,  &
    kdn(1,12),kdn(2,12) / 4HSPOI , 4HNT   /,  &
    kdn(1,13),kdn(2,13) / 4HCROD , 4H     /,  &
    kdn(1,14),kdn(2,14) / 4HOMIT , 4H     /,  &
    kdn(1,15),kdn(2,15) / 4HCNGR , 4HNT   /,  &
    kdn(1,16),kdn(2,16) / 4HASET , 4H     /,  &
    kdn(1,17),kdn(2,17) / 4HXXXX , 4H     /,  &
    kdn(1,18),kdn(2,18) / 4HXXXX , 4H     /,  &
    kdn(1,19),kdn(2,19) / 4HXXXX , 4H     /,  &
    kdn(1,20),kdn(2,20) / 4HXXXX , 4H     /


lf(i,j,n) = i + n*(j-1)

IF (param1 <= parama .AND. parama <= paramn) GO TO 20
WRITE  (nout,10) ufm,parama
10 FORMAT (a23,' 1738, UTILITY MODULE INPUT FIRST PARAMETER VALUE - '  &
    ,       i20,' OUT OF RANGE')
GO TO 9999

20 kor  = 10*nbuf + 1
nkor = korsz(x) - 10*nbuf
IF (nkor <= 0) CALL mesage (-8,nkor,mnam)
CALL page1
nlines = nlines + 8
WRITE  (nout,1)
1 FORMAT (//20X,'* U T I L I T Y   M O D U L E   I N P U T *',///,  &
    20X,'INPUT DATA ECHO (DATA READ VIA FORTRAN, REMEMBER TO ',  &
    'RIGHT ADJUST)', ///20X,'*   1  **   2  **   3  **   4  ',  &
    '**   5  **   6  **   7  **   8  **   9  **  10  *' ,///)
iox = 0
ioy = 0
IF (mach < 5 .OR. sperlk /= 0) GO TO 100

!     ON VAX-11/780 OR UNIX MACHINES, SEARCH FOR END OF BULK DATA DECK.

60 READ (nin,70,END=80) chr
70 FORMAT (a8)
IF (chr == e1 .OR. chr == e2 .OR. chr == e3) GO TO 100
GO TO 60

!     ENDDATA CARD NOT FOUND

80 WRITE  (nout,90) ufm
90 FORMAT (a23,' - "ENDDATA" CARD NOT FOUND BY INPUT MODULE')
CALL mesage (-37,0,mnam)

100 SELECT CASE ( parama )
  CASE (    1)
    GO TO 1000
  CASE (    2)
    GO TO 2000
  CASE (    3)
    GO TO 3000
  CASE (    4)
    GO TO 4000
  CASE (    5)
    GO TO 5000
  CASE (    6)
    GO TO 6000
  CASE (    7)
    GO TO 7000
  CASE (    8)
    GO TO 8000
END SELECT


!     PARAMA = 1 LAPLACE NETWORK

!     INPUT, ,,,,/,G2,,G4,/C,N,1/C,N,1 $         STATICS
!     INPUT, ,GEOM2,,GEOM4,/,G2,,G4,/C,N,1/C,N,1 $         STATICS
!     INPUT, ,,,,/,G2,,,/C,N,1/C,N,2 $          REAL-EIG W/O MASS COUPL
!     INPUT, ,GEOM2,,,/,G2,,,/C,N,1/C,N,2 $     REAL-EIG W/O MASS COUPL
!     INPUT, ,,,,/,G2,,,/C,N,1/C,N,3 $          REAL-EIG WITH MASS COUPL
!     INPUT, ,GEOM2,,,/,G2,,,/C,N,1/C,N,3 $     REAL-EIG WITH MASS COUPL


1000 SELECT CASE ( paramb )
  CASE (    1)
    GO TO 1100
  CASE (    2)
    GO TO 1200
  CASE (    3)
    GO TO 1300
END SELECT

1100 READ   (nin,1110) n,zk,u
1110 FORMAT (i8,2E8.0)
CALL page2 (-1)
WRITE  (nout,1111) n,zk,u
1111 FORMAT (21X,i8,1P,2E8.1,0P,f8.5)

ASSIGN 1140 TO r2
GO TO 1205
1140 CONTINUE

!     G4

ifil = 4
ASSIGN 1181 TO r
GO TO 9100

!     SPC

1181 ic = 3
ASSIGN 1182 TO r
GO TO 9200
1182 k(1) = 1000 + n
k(3) = 0
k(4) = 0
DO  i = 2,n
  k(2) = i
  CALL WRITE (FILE,k,4,0)
END DO
DO  i = 2,n
  k(2)  = lf(1,i,n1)
  qk(4) = u
  CALL WRITE (FILE,k,4,0)
  k(2)  = k(2) + n
  k(4)  = 0
  CALL WRITE (FILE,k,4,0)
END DO
DO  i = 2,n
  k(2)  = n*n1 + i
  CALL WRITE (FILE,k,4,0)
END DO
ASSIGN 1190 TO r1
GO TO 9600
1190 RETURN


1200 READ   (nin,1201) n,zk,zm
1201 FORMAT (i8,2E8.0)
CALL page2 (-1)
WRITE (nout,1111) n,zk,zm

1204 ASSIGN 1299 TO r2

1205 n1  = n + 1
nm1 = n - 1

!     G2

ifil = 2
ASSIGN 1211 TO r
GO TO 9100

!     CELAS4

1211 ic = 1
ASSIGN 1213 TO r
GO TO 9200
1213 qk(2) = zk
DO  j = 2,n
  DO  i = 1,n
    k(1) = lf(i,j,n1)
    k(3) = k(1)
    k(4) = k(3) + 1
    IF (paramb /= 1 .AND. i == 1) k(3) = 0
    IF (paramb /= 1 .AND. i == n) k(4) = 0
    CALL WRITE (FILE,k,4,0)
  END DO
END DO
DO  j = 1,n
  DO  i = 2,n
    k(3) = lf(i,j,n1)
    k(4) = k(3) + n1
    k(1) = k(3) + 1000000
    IF (paramb /= 1 .AND. j == 1) k(3) = 0
    IF (paramb /= 1 .AND. j == n) k(4) = 0
    CALL WRITE (FILE,k,4,0)
  END DO
END DO
ASSIGN 1216 TO r1
GO TO 9650
1216 IF (paramb == 1) GO TO 1240

!     CMASS4

ic = 2
ASSIGN 1218 TO r
GO TO 9200
1218 qk(2) = zm
k(4)  = 0
DO  j = 2,n
  DO  i = 2,n
    k(3) = lf(i,j,n1)
    k(1) = k(3) + 2000000
    CALL WRITE (FILE,k,4,0)
  END DO
END DO
IF (paramb == 3) GO TO 1230
1220 ASSIGN 1240 TO r1
GO TO 9650

1230 qk(2) = -f*zm
DO  j = 2,n
  DO  i = 1,n
    k(3) = lf(i,j,n1)
    k(1) = k(3) + 3000000
    k(4) = k(3) + 1
    IF (i == 1) k(3) = 0
    IF (i == n) k(4) = 0
    CALL WRITE (FILE,k,4,0)
  END DO
END DO
DO  j = 1,n
  DO  i = 2,n
    k(3) = lf(i,j,n1)
    k(4) = k(3) + n1
    k(1) = k(3) + 4000000
    IF (j == 1) k(3) = 0
    IF (j == n) k(4) = 0
    CALL WRITE (FILE,k,4,0)
  END DO
END DO
qk(2) = -f*zm/2.0
DO  j = 1,n
  DO  i = 1,n
    k(3) = lf(i,j,n1)
    k(1) = k(3) + 5000000
    k(4) = k(3) + n1 + 1
    IF (i == 1 .OR. j == 1) k(3) = 0
    IF (i == n .OR. j == n) k(4) = 0
    IF (k(3) /= 0 .OR. k(4) /= 0) CALL WRITE (FILE,k,4,0)
  END DO
END DO
DO  j = 1,n
  DO  i = 1,n
    k(3) = lf(i,j,n1)
    k(1) = k(3) + 6000000
    k(4) = k(3) + n1
    k(3) = k(3) + 1
    IF (i == n .OR. j == 1) k(3) = 0
    IF (i == 1 .OR. j == n) k(4) = 0
    IF (k(3) /= 0 .OR. k(4) /= 0) CALL WRITE (FILE,k,4,0)
  END DO
END DO
GO TO 1220
1240 IF (modcom(1) /= 0) GO TO 1295

!     DO NOT GENERATE CNGRNT DATA FOR N LESS THAN 3.

IF (n < 3) GO TO 1295

!     CNGRNT

ic = 15
ASSIGN 1245 TO r
GO TO 9200
1245 DO  j = 2,n
  DO  i = 1,n
    IF (paramb /= 1 .AND. (i == 1 .OR. i == n)) CYCLE
    k(1) = lf(i,j,n1)
    CALL WRITE (FILE,k,1,0)
  END DO
END DO
DO  j = 1,n
  IF (paramb /= 1 .AND. (j == 1 .OR. j == n)) CYCLE
  DO  i = 2,n
    k(1) = lf(i,j,n1) + 1000000
    CALL WRITE (FILE,k,1,0)
  END DO
END DO
k(1) = -1
CALL WRITE (FILE,k,1,0)
IF (paramb == 1) GO TO 1259
DO  j = 2,n
  DO  i = 1,n,nm1
    k(1) = lf(i,j,n1)
    CALL WRITE (FILE,k,1,0)
  END DO
END DO
DO  j = 1,n,nm1
  DO  i = 2,n
    k(1) = lf(i,j,n1) + 1000000
    CALL WRITE (FILE,k,1,0)
  END DO
END DO
k(1) = -1
CALL WRITE (FILE,k,1,0)
1259 CONTINUE
IF (paramb == 1) GO TO 1285
DO  j = 2,n
  DO  i = 2,n
    k(1) = lf(i,j,n1) + 2000000
    CALL WRITE (FILE,k,1,0)
  END DO
END DO
k(1) = -1
CALL WRITE (FILE,k,1,0)
IF (paramb == 2) GO TO 1285
DO  j = 2,n
  DO  i = 2,nm1
    k(1) = lf(i,j,n1) + 3000000
    CALL WRITE (FILE,k,1,0)
  END DO
END DO
DO  j = 2,nm1
  DO  i = 2,n
    k(1) = lf(i,j,n1) + 4000000
    CALL WRITE (FILE,k,1,0)
  END DO
END DO
k(1) = -1
CALL WRITE (FILE,k,1,0)
DO  j = 2,n
  DO  i = 1,n,nm1
    k(1) = lf(i,j,n1) + 3000000
    CALL WRITE (FILE,k,1,0)
  END DO
END DO
DO  j = 1,n,nm1
  DO  i = 2,n
    k(1) = lf(i,j,n1) + 4000000
    CALL WRITE (FILE,k,1,0)
  END DO
END DO
k(1) = -1
CALL WRITE (FILE,k,1,0)
DO  l = 1,2
  DO  j = 2,nm1
    DO  i = 2,nm1
      k(1) = lf(i,j,n1) + 1000000*l + 4000000
      CALL WRITE (FILE,k,1,0)
    END DO
  END DO
END DO
k(1) = -1
CALL WRITE (FILE,k,1,0)
DO  j = 1,n
  DO  i = 1,n
    IF (i /= 1 .AND. i /= n .AND. j /= 1 .AND. j /= n) CYCLE
    IF (j == 1 .AND. i == n .OR.  j == n .AND. i == 1) CYCLE
    k(1) = lf(i,j,n1) + 5000000
    CALL WRITE (FILE,k,1,0)
  END DO
END DO
DO  j = 1,n
  DO  i = 1,n
    IF (i /= 1 .AND. i /= n .AND. j /= 1 .AND. j /= n) CYCLE
    IF (j == 1 .AND. i == 1 .OR.  j == n .AND. i == n) CYCLE
    k(1) = lf(i,j,n1) + 6000000
    CALL WRITE (FILE,k,1,0)
  END DO
END DO
k(1) = -1
CALL WRITE (FILE,k,1,0)
1285 ASSIGN 1290 TO r1
GO TO 9600
1290 GO TO r2, (1140,1299)
1295 ASSIGN 1290 TO r
GO TO 9500
1299 RETURN

1300 READ   (nin,1301) n,zk,zm,f
1301 FORMAT (i8,3E8.0)
CALL page2 (-1)
WRITE (nout,1111) n,zk,zm,f
GO TO 1204


!     PARAMA = 2 RECTANGULAR FRAME MADE FROM BAR-S OR ROD-S

!     INPUT, ,,,,/G1,G2,,,/C,N,2/C,N,I/C,N,J $
!     INPUT GEOM1,GEOM2,,,/G1,G2,,,/C,N,2/C,N,I/C,N,J $
!            I=1  REGULAR BANDING
!            I=2  DOUBLE BANDING
!            I=3  ACTIVE COLUMN BANDING
!            I=4  REVERSE DOUBLE BANDING
!            J=0  BAR CONFIGURATION
!            J=1  ROD CONFIGURATION 1  (DIAGONALS IN BOTH DIRECTIONS)
!            J=2  ROD CONFIGURATION 2  (DIAGONALS IN LR TO UL DIRECTN)
!            J=3  ROD CONFIGURATION 3  (STATICALLY DETERMINATE)


2000 READ   (nin,2001) nx,ny,dx,dy,ip,lambda
2001 FORMAT (2I8,2E8.0,i8,e8.0)
CALL page2 (-1)
WRITE  (nout,2002) nx,ny,dx,dy,ip,lambda
2002 FORMAT (21X,2I8,1P,2E8.1,i8,1P,2E8.1)
nx1 = nx + 1
ny1 = ny + 1
ASSIGN 2295 TO r2

!     G1

2005 ifil = 1
ASSIGN 2010 TO r
GO TO 9100

!     GRID

2010 ic = 5
ASSIGN 2015 TO r
GO TO 9200
2015 qk(5) = 0.0
k(2)  = 0
k(6)  = 0
k(8)  = 0
sl  = SIN(degra*lambda)
cl  = COS(degra*lambda)
ddy = dy*cl
jj  = -1
DO  j = 1,ny1
  jj  = jj + 1
  IF (jj > ioy) jj = 0
  qk(4) = ddy*FLOAT(j-1)
  xo  = FLOAT(j-1)*sl
  ii  = -1
  DO  i = 1,nx1
    ii  = ii + 1
    IF (ii > iox) ii = 0
    k(1)  = lf(i,j,nx1)
    qk(3) = dx*FLOAT(i-1) + xo
    IF (ii == 0 .OR. jj == 0) GO TO 2018
    
    k(7) = 6
    
    GO TO 2020
    2018 k(7) = ip
    2020 CALL WRITE (FILE,k,8,0)
  END DO
END DO
2020  CONTINUE
ASSIGN 2050 TO r1
GO TO 9650

2050 IF (paramb == 1) GO TO 2290
ic = 8
ASSIGN 2060 TO r
GO TO 9200
2060 SELECT CASE ( paramb )
  CASE (    1)
    GO TO 2290
  CASE (    2)
    GO TO 2200
  CASE (    3)
    GO TO 2300
  CASE (    4)
    GO TO 2350
END SELECT

!     DOUBLE BANDING

2200 IF (MOD(ny,2) == 0) GO TO 2212
kk = ny1/2
nn = 1
GO TO 2214
2212 kk = ny1/2 + 1
nn = 2
2214 ij = ny1*nx1
DO  j = 1,ij
  k(1) = j
  iw = MOD(j,nx1)
  IF (iw == 0) iw = nx1
  il = (j-1)/nx1 + 1
  ilmk = il - kk
  SELECT CASE ( nn )
    CASE (    1)
      GO TO 2220
    CASE (    2)
      GO TO 2230
  END SELECT
  2220 IF (ilmk > 0) THEN
    GO TO  2222
  END IF
  2221 ill = -2*ilmk + 1
  GO TO 2223
  2222 ill  = 2*ilmk
  2223 k(2) = lf(iw,ill,nx1)
  GO TO 2240
  2230 IF (ilmk < 0) THEN
    GO TO  2231
  ELSE
    GO TO  2232
  END IF
  2231 ill  = -2*ilmk
  GO TO 2233
  2232 ill  = 2*ilmk + 1
  2233 k(2) = lf(iw,ill,nx1)
  2240 CALL WRITE (FILE,k,2,0)
END DO

2270 ASSIGN 2290 TO r1
GO TO 9650
2290 GO TO r2, (2295,3010)
2295 ASSIGN 2400 TO r
GO TO 9500

!     ACTIVE COLUMNS BANDING

2300 ij  = nx1*ny1
IF (MOD(ny,2) == 0) GO TO 2311
kk  = ij/2
kkk = 0
nn  = 1
GO TO 2315
2311 kk  = (ny/2+1)*nx1
kkk = ij - kk
nn  = 2
2315 DO  j = 1,ij
  k(1) = j
  SELECT CASE ( nn )
    CASE (    1)
      GO TO 2320
    CASE (    2)
      GO TO 2330
  END SELECT
  2320 IF (j-kk > 0) THEN
    GO TO  2322
  END IF
  2321 k(2) = j + kk
  GO TO 2340
  2322 k(2) = j - kk
  GO TO 2340
  2330 IF (j-kkk > 0) THEN
    GO TO  2332
  END IF
  2331 k(2) = j + kk
  GO TO 2340
  2332 k(2) = j - kkk
  2340 CALL WRITE (FILE,k,2,0)
END DO
GO TO 2270

!     REVERSE DOUBLE BANDING

2350 ij = nx1*ny1
IF (MOD(nx,2) == 0) GO TO 2360
kk = nx1/2
nn = 1
GO TO 2370
2360 kk = nx1/2 + 1
nn = 2
2370 DO  j = 1,ij
  k(1) = j
  iw   = MOD(j,nx1)
  IF (iw == 0) iw = nx1
  il   = (j-1)/nx1 + 1
  iwmk = iw - kk
  SELECT CASE ( nn )
    CASE (    1)
      GO TO 2380
    CASE (    2)
      GO TO 2385
  END SELECT
  2380 IF (iwmk > 0) THEN
    GO TO  2382
  END IF
  2381 iww  = -2*iwmk + 1
  GO TO 2383
  2382 iww  = 2*iwmk
  2383 k(2) = lf(il,iww,ny1)
  GO TO 2390
  2385 IF (iwmk < 0) THEN
    GO TO  2386
  ELSE
    GO TO  2387
  END IF
  2386 iww = -2*iwmk
  GO TO 2388
  2387 iww  = 2*iwmk + 1
  2388 k(2) = lf(il,iww,ny1)
  2390 CALL WRITE (FILE,k,2,0)
END DO
GO TO 2270

!     G2

2400 ifil = 2
ASSIGN 2410 TO r
GO TO 9100
2410 IF (paramc /= 0) GO TO 2700

!     CBAR

ic = 6
ASSIGN 2420 TO r
GO TO 9200
2420 k(2)  = 101
qk(5) = 0.0
qk(6) = 0.0
qk(7) = 1.0
k(8)  = 1
DO  i = 9,16
  k(i) = 0
END DO
DO  j = 1,ny1
  DO  i = 1,nx
    k(1) = lf(i,j,nx1)
    k(3) = k(1)
    k(4) = k(1) + 1
    CALL WRITE (FILE,k,16,0)
  END DO
END DO
DO  j = 1,ny
  DO  i = 1,nx1
    k(3) = lf(i,j,nx1)
    k(4) = k(3) + nx1
    k(1) = k(3) + 1000000
    CALL WRITE (FILE,k,16,0)
  END DO
END DO
2470 ASSIGN 2600 TO r1
GO TO 9650
2600 IF (modcom(1) /= 0) GO TO 2695

!     CNGRNT    (OUT OF SEQUENCE FOR CROD CASES)

ic = 15
ASSIGN 2610 TO r
GO TO 9200
2610 DO  j = 1,ny1
  DO  i = 1,nx
    k(1) = lf(i,j,nx1)
    CALL WRITE (FILE,k,1,0)
  END DO
END DO
k(1) = -1
CALL WRITE (FILE,k,1,0)
DO  j = 1,ny
  DO  i = 1,nx1
    k(1) = lf(i,j,nx1) + 1000000
    CALL WRITE (FILE,k,1,0)
  END DO
END DO
k(1) = -1
CALL WRITE (FILE,k,1,0)
IF (paramc == 0) GO TO 2680
loop2655:  DO  j = 1,ny
  DO  i = 1,nx
    k(1) = lf(i,j,nx1)*2 + 1999999
    CALL WRITE (FILE,k,1,0)
    IF (paramc == 3 .AND. j > 1) CYCLE loop2655
  END DO
END DO loop2655
k(1) = -1
CALL WRITE (FILE,k,1,0)
IF (paramc /= 1) GO TO 2680
DO  j = 1,ny
  DO  i = 1,nx
    k(1) = lf(i,j,nx1)*2 + 2000000
    CALL WRITE (FILE,k,1,0)
  END DO
END DO
k(1) = -1
CALL WRITE (FILE,k,1,0)
2680 ASSIGN 2690 TO r1
GO TO 9600
2690 RETURN
2695 ASSIGN 2690 TO r
GO TO 9500

!     CROD

2700 ic = 13
ASSIGN 2710 TO r
GO TO 9200
2710 k(2) = 101
DO  j = 1,ny1
  DO  i = 1,nx
    k(1) = lf(i,j,nx1)
    k(3) = k(1)
    k(4) = k(3) + 1
    CALL WRITE (FILE,k,4,0)
  END DO
END DO
DO  j = 1,ny
  DO  i = 1,nx1
    k(3) = lf(i,j,nx1)
    k(4) = k(3) + nx1
    k(1) = k(3) + 1000000
    CALL WRITE (FILE,k,4,0)
  END DO
END DO
loop2750:  DO  j = 1,ny
  DO  i = 1,nx
    k(3) = lf(i,j,nx1) + 1
    k(4) = k(3) + nx
    k(1) = 2*k(3) + 1999997
    CALL WRITE (FILE,k,4,0)
    IF (paramc == 3 .AND. j > 1) CYCLE loop2750
    IF (paramc /= 1) CYCLE
    k(1) = k(1) + 1
    k(3) = k(3) - 1
    k(4) = k(4) + 1
    CALL WRITE (FILE,k,4,0)
  END DO
END DO loop2750
GO TO 2470


!     PARAMA = 3 RECTANGULAR PLATE MADE FROM QUAD1-S

!     INPUT, ,,,,/G1,G2,,G4,/C,N,3/C,N,I $
!     INPUT GEOM1,GEOM2,,GEOM4,/G1,G2,,G4,/C,N,3/C,N,I $
!              I=1  REGULAR BANDING
!              I=2  DOUBLE BANDING
!              I=3  ACTIVE COLUMN BANDING
!              I=4  REVERSE DOUBLE BANDING


3000 READ  (nin,3001) nx,ny,dx,dy,ip,lambda,th
3001 FORMAT (2I8,2E8.0,i8,2E8.0)
CALL page2 (-2)
WRITE  (nout,2002) nx,ny,dx,dy,ip,lambda,th
READ   (nin,3002) iy0,ix0,iyl,ixw,iox,ioy
3002 FORMAT (6I8)
WRITE  (nout,3003) iy0,ix0,iyl,ixw,iox,ioy
3003 FORMAT (21X,6I8)
nx1 = nx + 1
ny1 = ny + 1

!     GRID

ASSIGN 3010 TO r2
GO TO 2005
3010 ASSIGN 3020 TO r
GO TO 9500

!     G2

3020 ifil = 2
ASSIGN 3030 TO r
GO TO 9100

!     CQUAD1

3030 IF (parama == 4) GO TO 4100
ic = 9
ASSIGN 3040 TO r
GO TO 9200
3040 k(2)  = 101
qk(7) = th
DO  j = 1,ny
  DO  i = 1,nx
    k(1) = lf(i,j,nx1)
    k(3) = k(1)
    k(4) = k(3) + 1
    k(6) = k(1) + nx1
    k(5) = k(6) + 1
    CALL WRITE (FILE,k,7,0)
  END DO
END DO
ASSIGN 3061 TO r1
GO TO 9650
3061 IF (modcom(1) /= 0) GO TO 3066

!     CNGRNT (OUT OF SEQUENCE)

ic = 15
ASSIGN 3062 TO r
GO TO 9200
3062 DO  j = 1,ny
  DO  i = 1,nx
    k(1) = lf(i,j,nx1)
    CALL WRITE (FILE,k,1,0)
  END DO
END DO
k(1) = -1
CALL WRITE (FILE,k,1,0)
3065 ASSIGN 3070 TO r1
GO TO 9600
3066 ASSIGN 3070 TO r
GO TO 9500

!     SPC-S AND OMIT-S

3070 IF (ip+iy0+ix0+iyl+ixw+iox+ioy == 0) GO TO 3090

!     G4

ifil = 4
ASSIGN 3080 TO r
GO TO 9100

!     SPC

3071 ic = 3
ASSIGN 3072 TO r
GO TO 9200
3072 k(1) = 1000*nx + ny
k(4) = 0
DO  i = 1,nx1
  k(2) = i
  k(3) = iy0
  IF (i ==    1) k(3) = iunion(iy0,ix0)
  IF (i ==  nx1) k(3) = iunion(iy0,ixw)
  IF (k(3) /= 0) CALL WRITE (FILE,k,4,0)
END DO
DO  i = 2,ny
  k(2) = lf(1,i,nx1)
  k(3) = ix0
  IF (k(3) /= 0) CALL WRITE (FILE,k,4,0)
  k(2) = k(2) + nx
  k(3) = ixw
  IF (k(3) /= 0) CALL WRITE (FILE,k,4,0)
END DO
k(2) = nx1*ny
DO  i = 1,nx1
  k(2) = k(2) + 1
  k(3) = iyl
  IF (i ==    1) k(3) = iunion(iyl,ix0)
  IF (i ==  nx1) k(3) = iunion(iyl,ixw)
  IF (k(3) /= 0) CALL WRITE (FILE,k,4,0)
END DO
ASSIGN 3089 TO r1
GO TO 9650
3080 IF (iox+ioy == 0) GO TO 3071

!     OMIT

ic = 14
ASSIGN 3081 TO r
GO TO 9200
3081 DO  i = 2,12,2
  k(i) = i/2
END DO
jj = 0
DO  j = 2,ny
  jj = jj + 1
  IF (jj > ioy) GO TO 3086
  ii = 0
  DO  i = 2,nx
    ii = ii + 1
    IF (ii > iox) GO TO 3084
    k(1) = lf(i,j,nx1)
    DO   l = 3,11,2
      k(l) = k(1)
    END DO
    
    CALL WRITE (FILE,k,10,0)
    
    CYCLE
    3084 ii = 0
  END DO
  CYCLE
  3086 jj = 0
END DO
ASSIGN 3088 TO r1
GO TO 9650
3088 IF (ip+iy0+ix0+iyl+ixw > 0) GO TO 3071
3089 ASSIGN 3090 TO r
GO TO 9500
3090 RETURN


!     PARAMA = 4 RECTANGULAR PLATE MADE FROM TRIA1-S

!     INPUT, ,,,,/G1,G2,,G4,/C,N,4/C,N,I/C,N,J $
!     INPUT GEOM1,GEOM2,,GEOM4,/G1,G2,,G4,/C,N,4/C,N,I/C,N,J $
!            I=1  REGULAR BANDING
!            I=2  DOUBLE BANDING
!            I=3  ACTIVE COLUMN BANDING
!            I=4  REVERSE DOUBLE BANDING
!            J=1  TRIANGLE CONFIGURATION OPTION NO. 1    (LL TO UR)
!            J=2  TRIANGLE CONFIGURATION OPTION NO. 2    (LR TO UL)


4000 GO TO 3000

!     CTRIA1

4100 ic = 10
ASSIGN 4200 TO r
GO TO 9200
4200 k(2)  = 101
qk(6) = th
DO  j = 1,ny
  DO  i = 1,nx
    k(3) = lf(i,j,nx1)
    k(4) = k(3) + 1
    k(1) = 2*k(3) - 1
    SELECT CASE ( paramc )
      CASE (    1)
        GO TO 4300
      CASE (    2)
        GO TO 4400
    END SELECT
    4300 k(5) = k(4) + nx1
    CALL WRITE (FILE,k,6,0)
    k(1) = k(1) + 1
    k(4) = k(3) + nx1
    k(3) = k(5)
    k(5) = k(4) - nx1
    CALL WRITE (FILE,k,6,0)
    CYCLE
    4400 k(5) = k(3) + nx1
    CALL WRITE (FILE,k,6,0)
    k(1) = k(1) + 1
    k(3) = k(5) + 1
    k(4) = k(5)
    k(5) = k(3) - nx1
    CALL WRITE (FILE,k,6,0)
  END DO
END DO
ASSIGN 4550 TO r1
GO TO 9650
4550 IF (modcom(1) /= 0) GO TO 3066

!     CNGRNT (OUT OF SEQUENCE)

ic = 15
ASSIGN 4600 TO r
GO TO 9200
4600 DO  j = 1,ny
  DO  i = 1,nx
    k(1) = lf(i,j,nx1)*2 - 1
    CALL WRITE (FILE,k,1,0)
  END DO
END DO
k(1) = -1
CALL WRITE (FILE,k,1,0)
DO  j = 1,ny
  DO  i = 1,nx
    k(1) = lf(i,j,nx1)*2
    CALL WRITE (FILE,k,1,0)
  END DO
END DO
k(1) = -1
CALL WRITE (FILE,k,1,0)
GO TO 3065


!     PARAMA = 5 N-SEGMENT STRING

!     INPUT, ,,,,/,G2,,,/C,N,5 $
!     INPUT, ,GEOM2,,,/,G2,,,/C,N,5 $


5000 READ   (nin,5010) n,xk1,xk2,xm,xb
5010 FORMAT (i8,4E8.0)
CALL page2 (-1)
WRITE  (nout,5011) n,xk1,xk2,xm,xb
5011 FORMAT (21X,i8,1P,4E8.1)
n1  = n + 1
nm1 = n - 1

!     G2

ifil = 2
ASSIGN 5100 TO r
GO TO 9100
5100 IF (xb == 0.0) GO TO 5140

!     CDAMP4

ic = 7
ASSIGN 5110 TO r
GO TO 9200
5110 qk(2) = xb
k(4)  = 0
DO  i = 2,n
  k(1)  = i + 2000000
  k(3)  = i
  CALL WRITE (FILE,k,4,0)
END DO
ASSIGN 5140 TO r1
GO TO 9650

!     CELAS4

5140 ic = 1
ASSIGN 5160 TO r
GO TO 9200
5160 qk(2) = xk1
DO  i = 1,n
  k(1) = i
  k(3) = i
  k(4) = i + 1
  IF (i == 1) k(3) = 0
  IF (i == n) k(4) = 0
  CALL WRITE (FILE,k,4,0)
END DO
IF (xk2 /= 0.0) GO TO 5190
5175 ASSIGN 5210 TO r1
GO TO 9650

5190 qk(2) = xk2
k(4)  = 0
DO  i = 2,n
  k(1)  = i + 3000000
  k(3)  = i
  CALL WRITE (FILE,k,4,0)
END DO
GO TO 5175

5210 IF (xm == 0.0) GO TO 5260

!     CMASS4

ic = 2
ASSIGN 5220 TO r
GO TO 9200
5220 qk(2) = xm
k(4)  = 0
DO  i = 2,n
  k(1)  = i + 1000000
  k(3)  = i
  CALL WRITE (FILE,k,4,0)
END DO
ASSIGN 5260 TO r1
GO TO 9650
5260 IF (modcom(1) /= 0) GO TO 5750
IF (n <= 2) GO TO 5750
IF (n == 3 .AND. xm == 0.0 .AND. xb == 0.0 .AND. xk2 == 0.0) GO TO 5750

!     CNGRNT

ic = 15
ASSIGN 5300 TO r
GO TO 9200
5300 IF (n == 3) GO TO 5400
DO  i = 2,nm1
  k(1) = i
  CALL WRITE (FILE,k,1,0)
END DO
k(1) = -1
CALL WRITE (FILE,k,1,0)
5400 IF (xm == 0.0) GO TO 5500
DO  i = 2,n
  k(1) = i + 1000000
  CALL WRITE (FILE,k,1,0)
END DO
k(1) = -1
CALL WRITE (FILE,k,1,0)
5500 IF (xb == 0.0) GO TO 5600
DO  i = 2,n
  k(1) = i + 2000000
  CALL WRITE (FILE,k,1,0)
END DO
k(1) = -1
CALL WRITE (FILE,k,1,0)
5600 IF (xk2 == 0.0) GO TO 5700
DO  i = 2,n
  k(1) = i + 3000000
  CALL WRITE (FILE,k,1,0)
END DO
k(1) = -1
CALL WRITE (FILE,k,1,0)
5700 ASSIGN 5900 TO r1
GO TO 9600
5750 ASSIGN 5900 TO r
GO TO 9500
5900 RETURN


!     PARAMA = 6 N-CELL BAR

!     INPUT, ,,,,/G1,G2,,G4,/C,N,6 $
!     INPUT GEOM1,GEOM2,,GEOM4,/G1,G2,,G4,/C,N,6 $


6000 READ   (nin,6010) n,xl,ip,iflg,ig0,m,iox
6010 FORMAT (i8,e8.0,5I8)
CALL page2 (-1)
WRITE  (nout,6011) n,xl,ip,iflg,ig0,m,iox
6011 FORMAT (21X,i8,1P,e8.1,5I8)
n1 = n + 1

!     G1

ifil = 1
ASSIGN 6100 TO r
GO TO 9100

!     GRID

6100 ic = 5
ASSIGN 6200 TO r
GO TO 9200
6200 k(2)  = 0
qk(4) = 0.0
qk(5) = 0.0
k(6)  = 0
k(8)  = 0
ii    = 0
DO  i = 1,n1
  ii    = ii + 1
  k(1)  = i
  qk(3) = xl*FLOAT(i-1)/FLOAT(n)
  IF (i == 1 .OR. ii > iox) GO TO 6280
  k(7) = 0
  GO TO 6300
  6280 k(7) = ip
  ii   = 0
  6300 CALL WRITE (FILE,k,8,0)
END DO
ASSIGN 6600 TO r1
GO TO 9600

!     G2

6600 ifil = 2
ASSIGN 6610 TO r
GO TO 9100

!     CBAR

6610 ic = 6
ASSIGN 6620 TO r
GO TO 9200
6620 k(2) = 101
k(8) = iflg
k(5) = ig0
k(6) = 0
k(7) = 0
IF (iflg == 2) GO TO 6635
READ   (nin,6630) x1,x2,x3
6630 FORMAT (3E8.0)
CALL page2 (-1)
WRITE  (nout,6631) x1,x2,x3
6631 FORMAT (21X,1P,3E8.1)
qk(5) = x1
qk(6) = x2
qk(7) = x3
GO TO 6640
6635 GO TO 9907
6640 DO  i = 9,16
  k(i) = 0
END DO
DO  i = 1,n
  k(1) = i
  k(3) = i
  k(4) = i + 1
  CALL WRITE (FILE,k,16,0)
END DO
IF (m <= 0 .OR. m > n-1) GO TO 6670
k(2) = 102
k(3) = 2
DO  i = 1,m
  k(1) = n + i
  k(4) = n - i + 2
  CALL WRITE (FILE,k,16,0)
END DO
6670 ASSIGN 6680 TO r1
GO TO 9650
6680 IF (modcom(1) /= 0) GO TO 6694

!     CNGRNT

ic = 15
ASSIGN 6685 TO r
GO TO 9200
6685 DO  i = 1,n
  k(1) = i
  CALL WRITE (FILE,k,1,0)
END DO
k(1) = -1
CALL WRITE (FILE,k,1,0)
ASSIGN 6695 TO r1
GO TO 9600
6694 ASSIGN 6695 TO r
GO TO 9500
6695 IF (iox == 0) RETURN

!     G4

ifil = 4
ASSIGN 6700 TO r
GO TO 9100

!     OMIT

6700 ic = 14
ASSIGN 6710 TO r
GO TO 9200
6710 DO  i = 2,12,2
  k(i) = i/2
END DO
ii = 0
DO  i = 2,n
  ii = ii + 1
  IF (ii > iox) GO TO 6716
  k(1) = i
  DO  l = 3,11,2
    k(l) = k(1)
  END DO
  CALL WRITE (FILE,k,12,0)
  CYCLE
  6716 ii = 0
END DO
ASSIGN 6730 TO r1
GO TO 9600
6730 RETURN


!     PARAMA = 7 FULL MATRIX AND OPTIONAL UNIT LOAD

!     INPUT, ,,,,/,G2,G3,,G5/C,N,7 $
!     INPUT, ,GEOM2,GEOM3,,/,G2,G3,,G5/C,N,7 $


7000 READ   (nin,7001) n,nsload
7001 FORMAT (2I8)
CALL page2 (-1)
WRITE  (nout,7002) n,nsload
7002 FORMAT (21X,2I8)
n1 = n + 1

!     G2

ifil = 2
ASSIGN 7010 TO r
GO TO 9100

!     CELAS4

7010 ic = 1
ASSIGN 7011 TO r
GO TO 9200
7011 qk(2) = 1.0
ii = 0
DO  i = 1,n
  IF (i > 1) ii = ii + n1 - i
  DO  j = i,n
    k(1) = ii + j
    k(3) = i
    k(4) = j
    IF (i == j) k(4) = 0
    CALL WRITE (FILE,k,4,0)
  END DO
END DO
ASSIGN 7030 TO r1
GO TO 9650
7030 IF (modcom(1) /= 0) GO TO 7036

!     DO NOT GENERATE CNGRNT DATA FOR N LESS THAN 3.

IF (n < 3) GO TO 7036

!     CNGRNT

ic = 15
ASSIGN 7032 TO r
GO TO 9200
7032 ii = 0
DO  i = 1,n
  IF (i > 1) ii = ii + n1 - i
  k(1) = ii + i
  CALL WRITE (FILE,k,1,0)
END DO
k(1) = -1
CALL WRITE (FILE,k,1,0)
ii = 0
DO  i = 1,n
  IF (i > 1) ii = ii + n1 - i
  DO  j = i,n
    IF (j == i) CYCLE
    k(1) = ii + j
    CALL WRITE (FILE,k,1,0)
  END DO
END DO
k(1) = -1
CALL WRITE (FILE,k,1,0)
ASSIGN 7037 TO r1
GO TO 9600
7036 ASSIGN 7037 TO r
GO TO 9500
7037 IF (nsload == 0) GO TO 7070

!     G3

ifil = 3
ASSIGN 7040 TO r
GO TO 9100

!     SLOAD

7040 ic = 11
ASSIGN 7050 TO r
GO TO 9200
7050 k(1)  = n
qk(3) = 1.0
DO  i = 1,n
  k(2)  = i
  CALL WRITE (FILE,k,3,0)
END DO
ASSIGN 7070 TO r1
GO TO 9600

!     G5

7070 ifil = 5
ASSIGN 7080 TO r
GO TO 9100

!     SPOINT

7080 ic = 12
ASSIGN 7090 TO r
GO TO 9200
7090 DO  i = 1,n
  CALL WRITE (FILE,i,1,0)
END DO
ASSIGN 7110 TO r1
GO TO 9600
7110 NE = n*n1/2
CALL page2 (-2)
WRITE  (nout,7201) n,NE
7201 FORMAT ('0*INPUT* FULL MATRIX OF ORDER',i9,'   GENERATED WITH',  &
    i9,'  ELEMENTS')
IF (nsload == 0) GO TO 7203
CALL page2 (-2)
WRITE  (nout,7202) n
7202 FORMAT ('0*INPUT* LOAD SET',i9,'  GENERATED')
7203 RETURN


!     PARAMA = 8 N-SPOKE WHEEL

!     INPUT, ,,,,/G1,G2,,,/C,N,8 $
!     INPUT GEOM1,GEOM2,,,/G1,G2,,,/C,N,8 $


8000 READ   (nin,8010) n,xl,ip,iflg,ig0,icen
8010 FORMAT (i8,e8.0,4I8)
CALL page2 (-1)
WRITE  (nout,8011) n,xl,ip,iflg,ig0,icen
8011 FORMAT (21X,i8,1P,e8.1,4I8)
n1 = n + 1

!     G1

ifil = 1
ASSIGN 8100 TO r
GO TO 9100
8100  CONTINUE

!     LOCATE AND COPY CORD2C CARD FROM THE FIRST INPUT FILE

ibuf = (ifil+4)*nbuf + 1
CALL preloc (*9908,x(ibuf),filin(ifil))
CALL locate (*9909,x(ibuf),cord2c,qk(3))
CALL READ (*9908,*8120,filin(ifil),qk(4),13,0,iflag)
CALL CLOSE (filin(ifil),1)
inopen(ifil) = .false.
8120  CONTINUE
ic = 17
ASSIGN 8200 TO r
GO TO 9200
8200  CONTINUE
CALL WRITE (FILE,qk(4),13,0)
ASSIGN 8250 TO r1
GO TO 9650

!     GRID

8250  CONTINUE
ic = 5
ASSIGN 8260 TO r
GO TO 9200
8260  CONTINUE
k(2)  = 2
!     QK(2) = QK(5)    THIS WILL ASSIGN REFERENCE NUMBER ON CORD2C CARD
!                      TO THE GRID POINTS

IF (n > 0 .AND. n < 256) GO TO 8050
CALL page2 (-2)
WRITE  (nout,8030) uwm
8030 FORMAT (a25,' 2369, WHEEL MUST HAVE FEWER THAN 256 SPOKES. ',  &
    'INPUT MODULE RESETTING TO 255')
n = 255
8050 n1 = n + 1
qk(3) = xl
qk(5) = 0.0
k(6)  = 2
k(8)  = 0
k(7)  = ip
DO  i = 1,n
  k(1)  = i
  qk(4) = 360.0/FLOAT(n)*FLOAT(i-1)
  CALL WRITE (FILE,k,8,0)
END DO
k(1)  = n1
k(2)  = 0
qk(3) = 0.0
qk(4) = 0.0
k(6)  = 0
IF (icen /= 0) k(7) = icen
CALL WRITE (FILE,k,8,0)
ASSIGN 8600 TO r1
GO TO 9600

!     G2

8600 ifil = 2
ASSIGN 8610 TO r
GO TO 9100

!     CBAR

8610 ic = 6
ASSIGN 8620 TO r
GO TO 9200
8620 k(2) = 101
k(8) = iflg
k(5) = ig0
k(6) = 0
k(7) = 0
IF (iflg == 2) GO TO 8635
READ   (nin,8630) x1,x2,x3
8630 FORMAT (3E8.0)
CALL page2 (-1)
WRITE  (nout,8631) x1,x2,x3
8631 FORMAT (21X,1P,3E8.1)
qk(5) = x1
qk(6) = x2
qk(7) = x3
GO TO 8640
8635 GO TO 9907
8640 DO  i = 9,16
  k(i) = 0
END DO
DO  i = 1,n
  k(1) = i
  k(3) = i
  k(4) = i + 1
  IF (k(4) == n1) k(4) = 1
  CALL WRITE (FILE,k,16,0)
END DO
k(4) = n1
DO  i = 1,n
  k(1) = n + i
  k(3) = i
  CALL WRITE (FILE,k,16,0)
END DO
ASSIGN 8950 TO r1
GO TO 9650
8950 ASSIGN 8900 TO r
GO TO 9500
8900 RETURN


!     UTILITY I/O ROUTINES


9100 FILE = filout(ifil)
ibuf = (ifil-1)*nbuf + 1
CALL gopen (FILE,x(ibuf),1)
t(1) = FILE
DO  j = 2,7
  t(j) = 0
END DO
fil  = filin(ifil)
ibuf = (ifil+4)*nbuf + 1
IF (parama - 8 == 0.0) THEN
  GO TO  9190
END IF
9115 CONTINUE
CALL OPEN (*9130,fil,x(ibuf),0)
inopen(ifil) = .true.
t(1) = fil
CALL rdtrl (t)
t(1) = FILE
CALL skprec (fil,1)
CALL fread (fil,hfil(1,ifil),3,0)
DO  j = 1,3
  IF (hfil(j,ifil) /= eee(j)) GO TO 9190
END DO
GO TO 9904
9130 inopen(ifil) = .false.
9190 GO TO r, (1181,1211,2010,2410,3030,3080,5100,6100,6610,6700,  &
    7010,7040,7080,8100,8610)


9200 IF (inopen(ifil)) GO TO 9230
9210 CALL WRITE (FILE,kl(ic),1,0)
CALL WRITE (FILE,16*(i1t(ic)-2)+j1t(ic),1,0)
CALL WRITE (FILE,kno(ic),1,0)
GO TO r, (1182,1213,1218,1245,2015,2060,2420,2610,2710,3040,  &
    3062,3072,3081,4200,4600,5110,5160,5220,5300,6200,  &
    6620,6685,6710,7011,7032,7050,7090,8200,8620,8260)
9230 knoic = kno(ic)
ksrt  = kdsort(knoic)
9235 knox  = hfil(3,ifil)
ksrtx = kdsort(knox)
IF (ksrt < ksrtx) GO TO 9210
IF (ksrt == ksrtx) GO TO 9906
CALL WRITE (FILE,hfil(1,ifil),3,0)
9240 CALL READ (*9902,*9250,fil,x(kor),nkor,0,rdflg)
CALL WRITE (FILE,x(kor),nkor,0)
GO TO 9240
9250 CALL WRITE (FILE,x(kor),rdflg,1)
CALL fread (fil,hfil(1,ifil),3,0)
DO  j = 1,3
  IF (hfil(j,ifil) /= eee(j)) GO TO 9235
END DO
inopen(ifil) = .false.
CALL CLOSE (fil,1)
GO TO 9210

9400 ktt  = kt(ic)
i1tt = i1t(ic)
j1tt = j1t(ic) + 16
t(i1tt) = orf(t(i1tt),two(j1tt))
IF (ktt == 1) GO TO 9450
i2tt = i2t(ic)
j2tt = j2t(ic) + 16
t(i2tt) = orf(t(i2tt),two(j2tt))
9450 GO TO r, (9620,9670)


9500 IF (inopen(ifil)) GO TO 9520
CALL WRITE (FILE,eee,3,1)
GO TO 9510
9505 CALL CLOSE (fil,1)
inopen(ifil) = .false.
9510 CALL wrttrl (t)
CALL CLOSE (FILE,1)
GO TO r, (1290,2400,2690,3020,3070,3090,5900,6695,7037,8900, 9630)
9520 CALL WRITE (FILE,hfil(1,ifil),3,0)
9525 CALL READ  (*9505,*9530,fil,x(kor),nkor,0,rdflg)
CALL WRITE (FILE,x(kor),nkor,0)
GO TO 9525
9530 CALL WRITE (FILE,x(kor),rdflg,1)
GO TO 9525

9600 CALL WRITE (FILE,0,0,1)
ASSIGN 9620 TO r
GO TO 9400
9620 ASSIGN 9630 TO r
GO TO 9500
9630 GO TO r1, (1190,1290,2690,3070,5900,6600,6695,6730,7037,7070,  &
    7110,8600,8900)

9650 CALL WRITE (FILE,0,0,1)
ASSIGN 9670 TO r
GO TO 9400
9670 GO TO r1, (1216,1240,2050,2290,2600,3061,3088,3089,4550,5140,  &
    5210,5260,6680,7030,8250,8950)

!     DIAGNOSTIC PROCESSING

9900 CALL mesage (m,FILE,mnam)
9902 m = -2
GO TO 9900
9904 WRITE  (nout,9954) sfm
9954 FORMAT (a25,' 1742, NO DATA PRESENT')
GO TO 9999
9906 WRITE  (nout,9956) ufm,kdn(1,ic),kdn(2,ic)
9956 FORMAT (a23,' 1744, DATA CARD(S) -',2A4,'- GENERATED BY UTILITY',  &
    ' MODULE INPUT NOT ALLOWED TO APPEAR IN BULK DATA')
GO TO 9999
9907 WRITE  (nout,9957) ufm
9957 FORMAT (a23,' 1745, UTILITY MODULE CANNOT HANDLE THE IFLG=2 CASE',  &
    ' SINCE THERE IS NO WAY TO GENERATE GRID POINT G0')
GO TO 9999
9908 m = -1
GO TO 9900
9909 WRITE  (nout,9959) ufm
9959 FORMAT (a23,' 1746, COORDINATE SYSTEM NOT DEFINED ON A CORD2C', ' CARD')

9999 m = -61
CALL page2 (-2)
GO TO 9900

END SUBROUTINE INPUT
