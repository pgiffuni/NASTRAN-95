SUBROUTINE dtria (iopt)
     
!     THIS ROUTINE GENERATES THE FOLLOWING
 
!     THREE 6X6 DIFFERENTIAL STIFFNESS MATRIX PARTITION FOR ONE PIVOT
!     POINT FOR A TRIA1, TRIA2 OR TORA3 ELEMENT.
 
 
!     CALLS FROM THIS ROUTINE ARE MADE TO
!           DTRBSC - BASIC BENDING TRI. ROUTINE.
!           DTRMEM - TRIANGLULAR MEMBRANE ROUTINE
!           TRANSD - SUPPLIES 3X3 TRANSFORMATIONS
!           INVERD - MATRIX INVERSION ROUTINE
!           GMMATD - GENERAL MATRIX MULITPLY AND TRANSPOSE ROUTINE
!           DS1B   - INSERTION ROUTINE
 
 
!          IOPT  = 1               2           3
!     ECPT INDEX   TRIA1           TRIA2       TRIA3       TRMEM
!     **********   *********       ********    ********    ********
!           1      EL ID           EL ID       EL ID        EL ID
!           2      SIL1            SIL1        SIL1         SIL1
!           3      SIL2            SIL2        SIL2         SIL2
!           4      SIL3            SIL3        SIL3         SIL3
!           5      THETA           THETA       MEM T1       THETA
!           6      MAT ID 1        MAT ID      MEM T2       MAT ID
!           7      T1              T           MEM T3       T
!           8      MAT ID 2        NSM         THETA        NSM
!           9      INERTIA I       CID1        FLAG FOR 8   CID1
!          10      MAT ID 3        X1          GRD OFFSET   X1
!          11      T2              Y1          MAT ID1      Y1
!          12      NSM             Z1          THICKNESS    Z1
!          13      Z1              CID2        MAT ID2      CID2
!          14      Z2              X2          INERTIA I    X2
!          15      CID1            Y2          MAT ID 3     Y2
!          16      X1              Z2          TS/T         Z2
!          17      Y1              CID3        NSM          CID3
!          18      Z1              X3          Z1           X3
!          19      CID2            Y3          Z2           Y3
!          20      X2              Z3          MAT ID 4     Z3
!          21      Y2              EL TEMP     THETA        EL TEMP
!          22      Z2                          FLAG FOR 21  EL DEFORM
!          23      CID3                        INTEGRATION  LOAD TEMP
!          24      X3              U1          STRESS ANGLE U1
!          25      Y3              V1          FLAG FOR 24  V2
!          26      Z3              W1          ZOFF1        W3
!          27      EL TEMP         U2          CID1         U2
!          28      EL DEFORM       V2          X1           V2
!          29      EL LOAD TEMP    W2          Y1           W2
!          30      U1 -DISP FOR U1 U3          Z1           U3
!          31      V1 -DISP FOR V1 V3          CID2         V3
!          32      W1 -DISP FOR Z1 W3          X2           W3
!          33      U2 -DISP FOR X2             Y2
!          34      V2 -DISP FOR Y2             Z2
!          35      W2 -DISP FOR Z2             CID3
!          36      U3 -DISP FOR X3             X3
!          37      V3 -DISP FOR Y3             Y3
!          38      W3 -DISP FOR Z3             Z3
!          39                                  EL TEMP
!          40
!          41
!          42                                  U1
!          43                                  V1
!          44                                  W1
!          45                                  U2
!          46                                  V2
!          47                                  W2
!          48                                  U3
!          49                                  V3
!          50                                  W3
 
 
 INTEGER, INTENT(IN OUT)                  :: iopt
 INTEGER :: subsca        ,subscb        ,subscc        , cid1
 DOUBLE PRECISION :: r             ,d1            ,habc          ,  &
     temp          ,d2            ,hinv          ,  &
     ksum          ,ivect         ,g             ,  &
     v             ,jvect         ,e             ,  &
     vv            ,kvect         ,tite          ,  &
     xsubb         ,temp9         ,tjte          ,  &
     xsubc         ,prod9         ,arr9          ,  &
     ysubc         ,u1            ,array9        ,  &
     t             ,u2            ,temp18        ,  &
     a             ,temp1         ,prod12        ,  &
     c1            ,temp2         ,hq            ,  &
     c2            ,l1            ,y1            ,  &
     x1            ,l2            ,y2            ,  &
     x2            ,s1            ,determ        ,  &
     s2            ,kout          ,s             , requiv
 DOUBLE PRECISION :: sigx          ,sigy          ,sigxy         ,  &
     stres         ,dumtwo
 DIMENSION necpt(100)    ,m(9)          ,requiv(8)     ,  &
     hq(12)        ,prod12(12)    ,habc(18)      ,  &
     g(36)         ,tite(18)      ,tjte(18)      ,  &
     kout(36)      ,temp18(18)    ,v1(3)         ,  &
     v2(3)         ,v3(3)         ,d1(3)         , d2(3)
 CHARACTER (LEN=29) :: uim
 CHARACTER (LEN=25) :: uwm, sfm
 CHARACTER (LEN=23) :: ufm
 COMMON /xmssg /  ufm           ,uwm           ,uim           , sfm
 COMMON /matin /  matid,inflag  ,eltemp,stress ,sinth,costh
 COMMON /matout/  g11,g12,g13   ,g22,g23,g33   ,rho,alpha1    ,  &
     alpha2,alp12  ,t_sub_0       ,g_sub_e       ,  &
     sigten,sigcom ,sigshe,g2x211 ,g2x212 ,g2x222
 COMMON /ds1aaa/  npvt          ,icstm         ,ncstm
 COMMON /ds1aet/  ecpt(100)
 COMMON /ds1adp/  a(54)         ,s(18)         ,hinv(36)      ,  &
     t(9)          ,temp9(9)      ,prod9(9)      ,  &
     arr9(9)       ,array9(9)     ,  &
     e(18)         ,temp          ,temp1         ,  &
     temp2         ,l1            ,l2            ,  &
     s1            ,s2            ,c1            ,  &
     c2            ,x1            ,x2            ,  &
     y1            ,y2            ,dumtwo(2)     ,  &
     determ        ,sigx          ,sigy          ,  &
     sigxy         ,xsubb         ,xsubc         ,  &
     ysubc         ,stres(3)      ,ksum(63)      ,  &
     ivect(3)      ,jvect(3)      ,kvect(3)      , r(2,4)        ,  &
     v(2)          ,vv(2)         ,u1            ,  &
     u2            ,npoint        ,km            ,  &
     subsca        ,subscb        ,subscc        ,  &
     npivot        ,ipvt          ,theta         ,  &
     nsubb         ,nsubc         ,ising         ,  &
     npt1          ,sinang        ,cosang
 COMMON /condas/  pi            ,twopi         ,radeg         ,  &
     degra         ,s4pisq
 COMMON /system/  ibuff         ,nout          ,nogo
 EQUIVALENCE (necpt(1),ecpt(1)) , (prod12(1),a(13))      ,  &
     (habc(1),a(25))    , (tite(1),a(37))        ,  &
     (tjte(1),s( 1))    , (kout(1),a(1))         ,  &
     (temp18(1),hinv(1)), (v1(1),ecpt(66))       ,  &
     (v2(1),ecpt(70))   , (v3(1),ecpt(74))       ,  &
     (requiv(1),r(1,1)) , (d1(1),a(1))           ,  &
     (d2(1),a(4))       , (hq(1),a(1))
 
 
 DATA     m    /  1,2,4,  2,3,4,  3,1,4 /,  cid1  / 65        /
 
 
!     THE ECPT DATA IS COPIED TO ECPT(PLUS 50)
!     THE DATA IN ECPT(BELOW 50) IS THEN PUT INTO TRMEM FORMAT TO BE
!     USED BY DTRMEM
!     THE DATA IN ECPT(ABOVE 50, SPECIALLY 51 THRU 62, 65 THRU 88) IS
!     PUT INTO TRIA1 FORMAT, WHICH WILL BE USED BY DTRBSC AND LOCALLY
 
 icid = cid1 - 4
 DO  i = 1,50
   ecpt(i+50) = ecpt(i)
 END DO
 SELECT CASE ( iopt )
   CASE (    1)
     GO TO 15
   CASE (    2)
     GO TO 25
   CASE (    3)
     GO TO 35
 END SELECT
 
!     TRIA1
 
 15 j = 15
 DO  i = 9,32
   ecpt(i) = ecpt(j)
   j = j + 1
 END DO
 GO TO 60
 
!     TRIA2
 
 25 ecpt(58) = ecpt(6)
 ecpt(59) =(ecpt(7)**3)/12.0
 ecpt(60) = ecpt(6)
 ecpt(61) = ecpt(7)
 
 j = 9
 DO  i = 65,88
   ecpt(i) = ecpt(j)
   j = j + 1
 END DO
 GO TO 60
 
!     TRIA3
 
!     IF NECPT(9)=0, ECPT(8) IS MATERIAL PROPERTY ORIENTAION ANGLE THETA
!     IF NECPT(9).NE.0, NECPT(8) IS MATERIAL COORDINATE SYSTEM ID. IN
!     THIS CASE, WE CAN NOT CONTINUE (NEED MORE STUFFS TO COMPUTE THETA,
!     SEE SHCSGD)
 
 35 IF (necpt(9) /= 0) GO TO 410
 ecpt(5) = ecpt( 8)
 ecpt(6) = ecpt(11)
 ecpt(7) = ecpt(12)
 j = 27
 DO  i = 9,32
   ecpt(i) = ecpt(j)
   j = j + 1
 END DO
 
 ecpt(55) = ecpt(58)
 j = 61
 DO  i = 56,60
   ecpt(i) = ecpt(j)
   j = j + 1
 END DO
 ecpt(61) = ecpt(62)
 j = 77
 DO  i = 65,88
   ecpt(i) = ecpt(j)
   j = j + 1
 END DO
 
 60 theta  = ecpt(5)*degra
 sinang = SIN(theta)
 cosang = COS(theta)
 sinth  = sinang
 costh  = cosang
 
 CALL dtrmem (2)
 
!     SIGX, SIGY , SIGXY ARE NOW AVAILABLE. SAVE THEM.
 
 stres(1) = sigx
 stres(2) = sigy
 stres(3) = sigxy
 
 eltemp = ecpt(21)
 
!     DETERMINE PIVOT POINT NUMBER
 
 DO  i = 1,3
   IF (npvt /= necpt(i+1)) CYCLE
   npivot = i
   GO TO 80
 END DO
 RETURN
 
!     FALL THRU ABOVE LOOP IMPLIES ERROR CONDITION
 
 80 CONTINUE
 
!     FORMATION OF THE R-MATRIX CONTAINING COORDINATES OF THE
!     SUB TRIANGLES. (2X4) FOR TRIANGULAR PLATE. (COLUMN 4 BLANK)
!     FORMATION ALSO OF THE I,J, AND K VECTORS USED IN THE E-MATRIX.
 
!     ZERO OUT R-MATRIX
 
 DO  i = 1,8
   requiv(i) = 0.0D0
 END DO
 
 DO  i = 1,3
   d2(i) = DBLE(v2(i)) - DBLE(v1(i))
   d1(i) = DBLE(v3(i)) - DBLE(v1(i))
 END DO
 
!     X2  GOES IN R(1,2)
 
 r(1,2) = DSQRT(d2(1)**2 + d2(2)**2 + d2(3)**2)
 DO  i = 1,3
   ivect(i) = d2(i)/r(1,2)
 END DO
 
!     NON-NORMALIZED K-VECTOR
 
 kvect(1) = ivect(2)*d1(3) - d1(2)*ivect(3)
 kvect(2) = ivect(3)*d1(1) - d1(3)*ivect(1)
 kvect(3) = ivect(1)*d1(2) - d1(1)*ivect(2)
 
!     Y3 GOES INTO R(2,3)
 
 r(2,3) = DSQRT(kvect(1)**2 + kvect(2)**2 + kvect(3)**2)
 DO  i = 1,3
   kvect(i) = kvect(i)/r(2,3)
 END DO
 
!     J-VECTOR = K X I  VECTORS
 
 jvect(1) = kvect(2)*ivect(3) - ivect(2)*kvect(3)
 jvect(2) = kvect(3)*ivect(1) - ivect(3)*kvect(1)
 jvect(3) = kvect(1)*ivect(2) - ivect(1)*kvect(2)
 
!     NORMALIZE J VECTOR TO MAKE SURE
 
 temp = DSQRT(jvect(1)**2 + jvect(2)**2 + jvect(3)**2)
 DO  i = 1,3
   jvect(i) = jvect(i)/temp
 END DO
 
!     X3 GOES INTO R(1,3) = D1 DOT IVECT
 
 r(1,3) = d1(1)*ivect(1) + d1(2)*ivect(2) + d1(3)*ivect(3)
 
!     CENTROID POINT GOES INTO R(1,4) AND R(2,4)
 
 r(1,4) = (r(1,2) + r(1,3))/3.0D0
 r(2,4) = r(2,3)/3.0D0
 
 
!     THE COORDINATES AND CENTROID OF THE PLATE IN THE ELEMENT
!     SYSTEM ARE STORED IN THE R-MATRIX WHERE THE COLUMN DENOTES THE
!     POINT AND THE ROW DENOTES THE X OR Y COORDINATE FOR ROW 1 OR
!     ROW 2 RESPECTIVELY.
 
 
!     SET UP THE M-MATRIX FOR MAPPING TRIANGLES, IN DATA STATEMENT.
 
!     ZERO OUT THE KSUM MATRIX FOR 63 AND THE GSUM MATRIX FOR 36
 
 DO  i = 1,63
   ksum(i) = 0.0D0
 END DO
 DO  i = 1,36
   g(i) = 0.0D0
 END DO
 
 DO  j = 1,3
   km = 3*j - 3
   subsca = m(km+1)
   subscb = m(km+2)
   subscc = m(km+3)
   
   DO  i = 1,2
     v(i)  = r(i,subscb) - r(i,subsca)
     vv(i) = r(i,subscc) - r(i,subsca)
   END DO
   xsubb = DSQRT(v(1)**2 + v(2)**2)
   u1    = v(1)/xsubb
   u2    = v(2)/xsubb
   xsubc = u1*vv(1) + u2*vv(2)
   ysubc = u1*vv(2) - u2*vv(1)
   
   sinth = sinang*u1 - cosang*u2
   costh = cosang*u1 + sinang*u2
   IF (ABS(sinth) < 1.0E-06) sinth = 0.0
   
!     AT THIS POINT, XSUBB, XSUBC, YSUBC ARE AT HAND FOR TRIANGLE -J-
   
   c2   = u1**2
   s2   = u2**2
   l1   = u1*u2
   sigx = c2*stres(1) + s2*stres(2) + 2.0D0*l1*stres(3)
   sigy = s2*stres(1) + c2*stres(2) - 2.0D0*l1*stres(3)
   sigxy=-l1*stres(1)+ l1*stres(2) + (c2-s2)*stres(3)
   ipvt = 0
   DO  i = 1,2
     npoint = km + i
     IF (m(npoint) == npivot) ipvt = i
   END DO
   CALL dtrbsc (1,ipvt)
   
!     NOW WE HAVE 6 MATRICES STORED AT A(1) TO A(54)- HIA,HIB,HIC
!                                                     HAC,HBC,HCC
   
!     NOW ADD CERTAIN OF THESE INTO THE SUMMED MATRICES
   
   
!     SET UP OF T-MATRIX
   
   t(1) = 1.0D0
   t(2) = 0.0D0
   t(3) = 0.0D0
   t(4) = 0.0D0
   t(5) = u1
   t(6) = u2
   t(7) = 0.0D0
   t(8) =-u2
   t(9) = u1
   
   DO  i = 1,3
     CALL gmmatd (t(1),3,3,1, a(9*i+19),3,3,0, temp9(1))
     CALL gmmatd (temp9(1),3,3,0, t(1),3,3,0,  prod9(1))
     
!     ADD THIS PRODUCT IN NOW.
!     COMPUTE POINTER TO KSUM MATRIX DESIRED.  (ZERO POINTER)
     
     npoint = km + i
     npoint = 9*m(npoint) + 18
     
     DO  k = 1,9
       nsubc  = npoint + k
       ksum(nsubc) = ksum(nsubc) + prod9(k)
     END DO
   END DO
   IF (ipvt == 0) GO TO 220
   DO  i = 1,2
     npoint = km +i
     npoint = 9*m(npoint) -9
     
!     TRANSFORM
     
     CALL gmmatd (t(1),3,3,1, a(9*i-8),3,3,0, temp9(1))
     CALL gmmatd (temp9(1),3,3,0, t(1),3,3,0, prod9(1))
     
!     INSERT
     
     DO  k = 1,9
       nsubc = k + npoint
       ksum(nsubc) = ksum(nsubc) + prod9(k)
     END DO
   END DO
   220 CONTINUE
   
!     FORM HQ (2X6)
   
   temp1 = xsubb - xsubc
   temp2 = ysubc**2
   l1 = DSQRT(xsubc**2 + temp2)
   l2 = DSQRT(temp1**2 + temp2)
   s1 = xsubc/l1
   s2 = temp1/l2
   c1 = ysubc/l1
   c2 = ysubc/l2
   x1 = xsubc/2.0D0
   y1 = ysubc/2.0D0
   x2 = (xsubb+xsubc)/2.0D0
   y2 = y1
   hq( 1) =-xsubc*c1
   hq( 2) = x1*s1 - y1*c1
   hq( 3) = 2.0D0*y1*s1
   hq( 4) =-3.0D0*x1*x1*c1
   hq( 5) = y1*(2.0D0*x1*s1 - y1*c1)
   hq( 6) = 3.0D0*y1*y1*s1
   hq( 7) = 2.0D0*x2*c2
   hq( 8) = x2*s2 + y2*c2
   hq( 9) = 2.0D0*y2*s2
   hq(10) = 3.0D0*x2*x2*c2
   hq(11) = y2*(2.0D0*x2*s2 + y2*c2)
   hq(12) = 3.0D0*y2*y2*s2
   
!                      I                    -1
!     COMPUTE (H       I  H     )  = (HQ)(H)    STORE IN PROD12
!               PSI,B  I   PSI,C
!                      I
   
   
   CALL gmmatd (hq(1),2,6,0, hinv(1),6,6,0, prod12(1))
   
   
!     COMPUTE (H     ) = -(PROD12)(S)
!               PSI,A
   
   CALL gmmatd (prod12(1),2,6,0, s(1),6,3,0, habc(1))
   
   habc(1) = -habc(1)
   habc(2) = -habc(2) + s1
   habc(3) = -habc(3) + c1
   habc(4) = -habc(4)
   habc(5) = -habc(5) + s2
   habc(6) = -habc(6) - c2
   
!     SPLIT (H     ) AND (H     )    PARTITION
!             PSI,B        PSI,C
   
   habc( 7) = prod12( 1)
   habc( 8) = prod12( 2)
   habc( 9) = prod12( 3)
   habc(10) = prod12( 7)
   habc(11) = prod12( 8)
   habc(12) = prod12( 9)
   habc(13) = prod12( 4)
   habc(14) = prod12( 5)
   habc(15) = prod12( 6)
   habc(16) = prod12(10)
   habc(17) = prod12(11)
   habc(18) = prod12(12)
   
!     MAP  H , H , AND H  INTO THE G-MATRICES.
!           A   B       C
   
!     TRIANGLE NUMBER = J, THE THREE POINTS ARE SUBSCA, SUBSCB, SUBSCC.
   
   DO  i = 1,3
     
!     POINTER TO H  = 6*I-6
!                 I
     
     
!     TRANSFORM H SUB I
     
     CALL gmmatd (habc(6*i-5),2,3,0, t(1),3,3,0, temp9(1))
     
     
     npoint = km + i
     npoint = 9*m(npoint) - 9
     
!     J = 1    ROW 1 OF H INTO ROW 1 OF G.
!              ROW 2 OF H INTO ROW 2 OF G.
!     J = 2    ROW 1 OF H INTO ROW 2 OF G.
!              ROW 2 OF H INTO ROW 3 OF G.
!     J = 3    ROW 1 OF H INTO ROW 3 OF G.
!              ROW 2 OF H INTO ROW 1 OF G.
     
     IF (j-2 < 0) THEN
       GO TO   240
     ELSE IF (j-2 == 0) THEN
       GO TO   230
     ELSE
       GO TO   260
     END IF
     
     230 npoint = npoint + 3
     240 DO  k = 1,6
       npoint = npoint + 1
       g(npoint) = g(npoint) + temp9(k)
     END DO
     CYCLE
     260 g(npoint+7) = g(npoint+7) + temp9(1)
     g(npoint+8) = g(npoint+8) + temp9(2)
     g(npoint+9) = g(npoint+9) + temp9(3)
     g(npoint+1) = g(npoint+1) + temp9(4)
     g(npoint+2) = g(npoint+2) + temp9(5)
     g(npoint+3) = g(npoint+3) + temp9(6)
     
   END DO
   
   
!     END OF LOOP FOR BASIC TRIANGLES
   
 END DO
 
 
!     FILL E-MATRIX
 
 DO  i = 1,18
   e(i) = 0.0D0
 END DO
 e( 1) = kvect(1)
 e( 4) = kvect(2)
 e( 7) = kvect(3)
 e(11) = ivect(1)
 e(14) = ivect(2)
 e(17) = ivect(3)
 e(12) = jvect(1)
 e(15) = jvect(2)
 e(18) = jvect(3)
 
!              T
!     FORM   T   E      STORE IN TITE-MATRIX (6X3)
!             I
 
 IF (necpt(4*npivot+icid) == 0) GO TO 300
 CALL transd (necpt(4*npivot+icid),t(1))
 CALL gmmatd (t(1),3,3,1, e( 1),3,3,0, tite( 1))
 CALL gmmatd (t(1),3,3,1, e(10),3,3,0, tite(10))
 GO TO 320
 
 300 DO  k = 1,18
   tite(k) = e(k)
 END DO
 
!     SOLVE NOW FOR
 
!       E                   T     T                       T
!    (K  ) = (K  ) - (TERM ) (K  ) - (K  )(TERM ) + (TERM )(K  )(TERM )
!      IJ      IJ         I    J4      I4      J         I   44      J
 
!                           -1                               I=NPIVOT
!      WHERE  (TERM ) = (G )  (G ) ,I=NPIVOT                 J=1,2,3
!                  I      4     I
 
!                           -1
!             (TERM ) = (G )  (G ) ,J=1,2,3 AS ABOVE
!                  J      4     J
 
!     AND WITH TRANSFORMATIONS
 
!       G        T      E   T
!    (K  ) = (C ) (E)(K  )(E )(C )
!      IJ      I       IJ       J
 
 
!     COMPUTE  (TERM        )  STORE IN PROD9
!                   I=NPIVOT
 
!                  -1
!     FIRST GET (G )
!                 4
 
 320 CONTINUE
 
!     NO NEED TO COMPUTE DETERMINANT SINCE IT IS NOT USED SUBSEQUENTLY.
 
 ising = -1
 CALL inverd (3,g(28),3,prod9,0,determ,ising,temp9)
 
 CALL gmmatd (g(28),3,3,0, g(9*npivot-8),3,3,0, prod9(1))
 
!                       T
!     GET  (TERM        )(K  ) -(K  )  STORE IN TEMP9
!               I=NPIVOT   44     I4
 
 CALL gmmatd (prod9(1),3,3,1, ksum(55),3,3,0, temp9(1))
 DO  k = 1,9
   npoint = 9*npivot + 18 + k
   temp9(k) = temp9(k) - ksum(npoint)
 END DO
 
 
!     THE TWO COMMON PRODUCTS ARE NOW AT HAND IN PROD9 AND TEMP9.
 
 DO  j = 1,3
   
!                   T     T
!     (TERM        ) (K  )    STORE IN ARR9
!          I=NPIVOT    J4
   
   CALL gmmatd (prod9(1),3,3,1, ksum(9*j+19),3,3,1, arr9(1))
   
!     SUBTRACT FROM (K  )
!                     IJ
   
   nbegin = 9*j - 9
   DO  i = 1,9
     npoint = nbegin + i
     ksum(npoint) = ksum(npoint) - arr9(i)
   END DO
   
   
!      COMPUTE  (TERM )  STORE IN ARR9
!                   J
   
   CALL gmmatd (g(28),3,3,0, g(9*j-8),3,3,0, arr9(1))
   
!                            T
!     COMPUTE ((TERM        )(K  ) -(K  )) (TERM ) = (TEMP9)(ARR9)
!                   I=NPOINT   44     I4        J
   
   CALL gmmatd (temp9(1),3,3,0, arr9(1),3,3,0, array9(1))
   
!     ADD TO K
!             IJ
   
   DO  i = 1,9
     npoint = nbegin + i
     ksum(npoint) = ksum(npoint) + array9(i)
   END DO
   
!       E
!     K    COMPLETE
!      IJ
   
!     TRANSFORM NOW, AND INSERT.
   
   
!     TRANSFORMATIONS AND INSERTION
   
   IF (necpt(4*j+icid) == 0) GO TO 370
   CALL transd (necpt(4*j+icid),t(1))
   CALL gmmatd (t(1),3,3,1, e( 1),3,3,0, tjte( 1))
   CALL gmmatd (t(1),3,3,1, e(10),3,3,0, tjte(10))
   GO TO 390
   
   370 DO  k = 1,18
     tjte(k) = e(k)
   END DO
   390 CALL gmmatd (ksum(nbegin+1),3,3,0, tjte(1),6,3,1, temp18(1))
   CALL gmmatd (tite(1),6,3,0, temp18(1),3,6,0, kout(1))
   CALL ds1b (kout(1),necpt(j+1))
 END DO
 RETURN
 
!     COULD NOT DO IT
 
 410 WRITE  (nout,420) sfm
 420 FORMAT (a25,', DEFFICIENT SOURCE CODE IN DTRIA TO HANDLE CTRIA3 ',  &
     'ELEMENT WITH MATERIAL', /5X,  &
     'PROPERTY COORD. SYSTEM. ANGLE MUST BE SPECIFIED')
 nogo = 1
 RETURN
END SUBROUTINE dtria
