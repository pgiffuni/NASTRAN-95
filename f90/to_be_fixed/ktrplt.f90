SUBROUTINE ktrplt
     
!     THIS ROUTINE GENERATES THE FOLLOWING
 
!                        3-6X6 STIFFNESS MATRICES WITH RESPECT
!                        TO ONE PIVOT POINT OF A TRIANGULAR PLATE
!                        ELEMENT.
 
!         REF. FMMS-55  NOVEMBER 1ST, 1967
 
!         CALLS FROM THIS ROUTINE ARE MADE TO
!                        KTRBSC - BASIC BENDING TRI. ROUTINE.
!                        TRANSD - SUPPLIES 3X3 TRANSFORMATIONS
!                        INVERD - MATRIX INVERSION ROUTINE
!                        SMA1B  - INSERTION ROUTINE
!                        GMMATD - GENERAL MATRIX MULITPLY AND
!                                 TRANSPOSE ROUTINE
!                        MESAGE - ERROR MESSAGE WRITER
 
 
 INTEGER :: subsca        ,subscb        ,subscc
 DOUBLE PRECISION :: r             ,d1            ,habc  &
     ,temp          ,d2            ,hinv ,ksum          ,ivect         ,g  &
     ,v             ,jvect         ,e ,vv            ,kvect         ,tite  &
     ,xsubb         ,temp9         ,tjte ,xsubc         ,prod9         ,arr9  &
     ,ysubc         ,u1            ,array9  &
     ,t             ,u2            ,temp18  &
     ,a             ,temp1         ,prod12 ,c1            ,temp2         ,hq  &
     ,c2            ,l1            ,y1 ,x1            ,l2            ,y2  &
     ,x2            ,s1            ,determ  &
     ,s2            ,kout          ,s ,requiv
!     ******************************************************************
 
!     ECPT LISTS AS OF AUGUST 4, 1967
 
!                 DEFINITION
!       ECPT      TRI.PLATE AND BASIC BENDING TRI.
!     ******************************************************************
!     ECPT( 1) = ELEMENT ID         INTEGER
!     ECPT( 2) = GRID PT. A         INTEGER
!     ECPT( 3) = GRID PT. B         INTEGER
!     ECPT( 4) = GRID PT. C         INTEGER
!     ECPT( 5) = THETA              REAL
!     ECPT( 6) = MAT ID 1           INTEGER
!     ECPT( 7) = I  MOM. OF INERT.  REAL
!     ECPT( 8) = MAT ID 2           INTEGER
!     ECPT( 9) = T2                 REAL
!     ECPT(10) = NON-STRUCT. MASS   REAL
!     ECPT(11) = Z1                 REAL
!     ECPT(12) = Z2                 REAL
!     ECPT(13) = COORD. SYS. ID 1   INTEGER
!     ECPT(14) = X1                 REAL
!     ECPT(15) = Y1                 REAL
!     ECPT(16) = Z1                 REAL
!     ECPT(17) = COORD. SYS. ID 2   INTEGER
!     ECPT(18) = X2                 REAL
!     ECPT(19) = Y2                 REAL
!     ECPT(20) = Z2                 REAL
!     ECPT(21) = COORD. SYS. ID 3   INTEGER
!     ECPT(22) = X3                 REAL
!     ECPT(23) = Y3                 REAL
!     ECPT(24) = Z3                 REAL
!     ECPT(25) = ELEMENT TEMP       REAL
!     ******************************************************************
 DIMENSION necpt(100)    ,m(9)          ,requiv(8)  &
     ,hq(12)        ,prod12(12)    ,habc(18)  &
     ,g(36)         ,tite(18)      ,tjte(18)  &
     ,kout(36)      ,temp18(18)    ,v1(3) ,v2(3)         ,v3(3)         ,r(2,4)  &
     ,d1(3)         ,d2(3)
 
 COMMON /condas/    consts(5)
 COMMON /matin /    matid         ,inflag        ,eltemp  &
     ,stress        ,sinth         ,costh
 COMMON /matout/    g11,g12,g13   ,g22,g23,g33   ,rho  &
     ,alpha1        ,alpha2        ,alp12 ,t sub 0       ,g sub e       ,sigten  &
     ,sigcom        ,sigshe        ,g2x211 ,g2x212        ,g2x222
 COMMON /sma1io/    dum1(10)      ,ifkgg         ,dum2(1)  &
     ,if4gg         ,dum3(23)
 COMMON /sma1cl/    iopt4         ,k4ggsw        ,npvt  &
     ,dumcl(7)      ,link(10)      ,idetck ,dodet         ,nogo
 COMMON /sma1et/    ecpt(100)
 COMMON /sma1dp/ a(81)         ,s(18)         ,t(9)  &
     ,temp9(9)      ,prod9(9)      ,arr9(9)  &
     ,array9(9)     ,hinv(36)      ,ksum(63)  &
     ,xsubb         ,xsubc         ,ysubc ,e(18)         ,temp          ,l1  &
     ,l2            ,s1            ,s2 ,c1            ,c2            ,x1  &
     ,x2            ,y1            ,y2  &
     ,temp1         ,temp2         ,dumtwo(2) ,determ  &
     ,npoint        ,km            ,subsca  &
     ,subscb        ,subscc        ,npivot ,theta         ,nsubc         ,ising  &
     ,npt1          ,v(2)          ,vv(2)  &
     ,ivect(3)      ,jvect(3)      ,kvect(3)  &
     ,u1            ,u2            ,sinang ,cosang
 
 EQUIVALENCE       (consts(4),degra)
 EQUIVALENCE (necpt(1),ecpt(1))  &
     ,(prod12(1),a(13)) ,(habc(1),a(25))  &
     ,(tite(1),a(37)) ,(tjte(1),a(55))  &
     ,(kout(1),a(1)) ,(temp18(1),hinv(1))  &
     ,(v1(1),ecpt(14)) ,(v2(1),ecpt(18))  &
     ,(v3(1),ecpt(22)) ,(requiv(1),r(1,1))  &
     ,(d1(1),a(1)) ,(d2(1),a(4))  &
     ,(hq(1),a(1))
 
 DATA  m/ 1,2,4,   2,3,4,   3,1,4 /
 
 eltemp = ecpt(25)
!     DETERMINE PIVOT POINT NUMBER
 
 DO  i=1,3
   IF( npvt /= necpt(i+1) ) CYCLE
   npivot = i
   GO TO 20
 END DO
 
 
!     FALL THRU ABOVE LOOP IMPLIES ERROR CONDITION
 CALL mesage(-30,34,ecpt(1))
 
 20 theta = ecpt(5) * degra
 sinang = SIN( theta )
 cosang = COS( theta )
!     FORMATION OF THE R-MATRIX CONTAINING COORDINATES OF THE
!     SUB TRIANGLES. (2X4) FOR TRIANGULAR PLATE. (COLUMN 4 BLANK)
!     FORMATION ALSO OF THE I,J, AND K VECTORS USED IN THE E-MATRIX.
 
!     ZERO OUT R-MATRIX
 DO  i=1,8
   requiv(i)=0.0D0
 END DO
 
 DO  i=1,3
   d2(i) = DBLE( v2(i) ) - DBLE( v1(i) )
   d1(i) = DBLE( v3(i) ) - DBLE( v1(i) )
 END DO
 
!     X2  GOES IN R(1,2)
 r(1,2) = DSQRT ( d2(1)**2  +  d2(2)**2  +  d2(3)**2 )
 IF (r(1,2) < 1.0D-7) GO TO 370
 DO  i=1,3
   ivect(i) = d2(i) / r(1,2)
 END DO
 
!     NON-NORMALIZED K-VECTOR
 kvect(1) = ivect(2) * d1(3)    -    d1(2) * ivect(3)
 kvect(2) = ivect(3) * d1(1)    -    d1(3) * ivect(1)
 kvect(3) = ivect(1) * d1(2)    -    d1(1) * ivect(2)
 
!     Y3 GOES INTO R(2,3)
 r(2,3) = DSQRT ( kvect(1)**2 + kvect(2)**2 + kvect(3)**2 )
 IF (r(2,3) < 1.0D-7) GO TO 370
 DO  i=1,3
   kvect(i) = kvect(i) / r(2,3)
 END DO
 
!     J-VECTOR = K X I  VECTORS
 jvect(1) = kvect(2) * ivect(3) - ivect(2) * kvect(3)
 jvect(2) = kvect(3) * ivect(1) - ivect(3) * kvect(1)
 jvect(3) = kvect(1) * ivect(2) - ivect(1) * kvect(2)
!     NORMALIZE J VECTOR TO MAKE SURE
 temp = DSQRT ( jvect(1)**2 + jvect(2)**2 + jvect(3)**2 )
 IF (temp < 1.0D-7) GO TO 370
 DO  i=1,3
   jvect(i) = jvect(i) / temp
 END DO
!     X3 GOES INTO R(1,3) = D1 DOT IVECT
 r(1,3) = d1(1) * ivect(1)  +  d1(2) * ivect(2) + d1(3) * ivect(3)
 
!     CENTROID POINT GOES INTO R(1,4) AND R(2,4)
 r(1,4) = (  r(1,2) + r(1,3)  ) / 3.0D0
 r(2,4) = r(2,3) / 3.0D0
!     ******************************************************************
!            THE COORDINATES AND CENTROID OF THE PLATE IN THE ELEMENT
!     SYSTEM ARE STORED IN THE R-MATRIX WHERE THE COLUMN DENOTES THE
!     POINT AND THE ROW DENOTES THE X OR Y COORDINATE FOR ROW 1 OR
!     ROW 2 RESPECTIVELY.
!     ******************************************************************
 
!     SET UP THE M-MATRIX FOR MAPPING TRIANGLES, IN DATA STATEMENT.
 
!     ******************************************************************
!     ZERO OUT THE KSUM MATRIX FOR 63 AND THE GSUM MATRIX FOR 36...
 
 DO  i=1,63
   ksum(i) = 0.0D0
 END DO
 DO  i=1,36
   g(i) = 0.0D0
 END DO
 
 
 DO  j=1,3
   km = 3*j - 3
   
   subsca = m(km+1)
   subscb = m(km+2)
   subscc = m(km+3)
   
   DO  i=1,2
     v(i) = r(i,subscb) - r(i,subsca)
     vv(i)= r(i,subscc) - r(i,subsca)
   END DO
   xsubb = DSQRT ( v(1)**2 + v(2)**2 )
   u1 = v(1) / xsubb
   u2 = v(2) / xsubb
   xsubc = u1 * vv(1) + u2 * vv(2)
   ysubc = u1 * vv(2) - u2 * vv(1)
   
   sinth = sinang * u1  -  cosang * u2
   costh = cosang * u1  +  sinang * u2
   IF(ABS(sinth) < 1.0E-06) sinth = 0.0E0
   
!     AT THIS POINT, XSUBB, XSUBC, YSUBC ARE AT HAND FOR
!     TRIANGLE -J-
   
   CALL ktrbsc( 2 )
!                         U
!     NOW HAVE AT HAND  K    I,J, =1,2,3.   9-3X3 MATRICES STORED AT
!                        IJ                 A(1) THROUGH A(81).
   
!           -1
!     ALSO H   (6X6) AT A(145) TO A(181) AND S (6X3) AT A(82) TO A(99)
   
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
   
   DO  i=1,3
     CALL gmmatd( t(1),3,3,1,  a(27*i-8),3,3,0,  temp9(1) )
     CALL gmmatd( temp9(1),3,3,0,  t(1),3,3,0,  prod9(1) )
     
!     ADD THIS PRODUCT IN NOW.
!     COMPUTE POINTER TO KSUM MATRIX DESIRED.  (ZERO POINTER)
     npoint = km + i
     npoint = 9*m(npoint) + 18
     
     DO  k=1,9
       nsubc  = npoint + k
       ksum(nsubc) = ksum(nsubc) + prod9(k)
     END DO
   END DO
   DO  k=1,2
     npoint = km + k
     IF( m(npoint) /= npivot ) CYCLE
     CALL gmmatd( t(1),3,3,1,  a(36*k-35),3,3,0,  temp9(1) )
     CALL gmmatd( temp9(1),3,3,0,  t(1),3,3,0,  prod9(1) )
     
!     COMPUTE POINTER TO KSUM MATRIX (ZERO POINTER)
     
     npoint = 9 * npivot - 9
     DO  i=1,9
       nsubc = npoint + i
       ksum(nsubc) = ksum(nsubc) + prod9(i)
     END DO
     
     CALL gmmatd(t(1),3,3,1,  a(18*k-8),3,3,0,  temp9(1) )
     CALL gmmatd( temp9(1),3,3,0,  t(1),3,3,0,  prod9(1) )
     
!     COMPUTE ZERO POINTER TO KSUM MATRIX DESIRED
     
     npoint = km + 3 - k
     npoint = 9 * m(npoint) - 9
     DO  i=1,9
       nsubc = npoint + i
       ksum(nsubc) = ksum(nsubc) +  prod9(i)
     END DO
   END DO
   
!     FORM HQ (2X6)
   
   temp1 = xsubb - xsubc
   temp2 = ysubc ** 2
   l1 = DSQRT( xsubc**2 + temp2 )
   l2 = DSQRT( temp1**2 + temp2 )
   s1 = xsubc / l1
   s2 = temp1 / l2
   c1 = ysubc / l1
   c2 = ysubc / l2
   x1 = xsubc / 2.0D0
   y1 = ysubc / 2.0D0
   x2 = (xsubb + xsubc) / 2.0D0
   y2 = y1
   hq( 1) = -xsubc * c1
   hq( 2) = x1 * s1 - y1 * c1
   hq( 3) = 2.0D0 * y1 * s1
   hq( 4) = -3.0D0 * x1 * x1 * c1
   hq( 5) = y1 * (2.0D0 * x1 * s1 - y1 * c1 )
   hq( 6) = 3.0D0 * y1 * y1 * s1
   hq( 7) = 2.0D0 * x2 * c2
   hq( 8) = x2 * s2  + y2 * c2
   hq( 9) = 2.0D0 * y2 * s2
   hq(10) = 3.0D0 * x2 * x2 * c2
   hq(11) = y2 * ( 2.0D0 * x2 * s2 + y2 * c2 )
   hq(12) = 3.0D0 * y2 * y2 * s2
   
!                      I                    -1
!     COMPUTE (H       I  H     )  = (HQ)(H)    STORE IN PROD12
!               PSI,B  I   PSI,C
!                      I
   
   
   CALL gmmatd( hq(1),2,6,0, hinv(1),6,6,0,  prod12(1) )
   
   
!     COMPUTE (H     ) = -(PROD12)(S)
!               PSI,A
   
   CALL gmmatd( prod12(1),2,6,0, s(1),6,3,0,  habc(1) )
   
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
   
   DO  i=1,3
     
!     POINTER TO H  = 6*I-6
!                 I
     
     
!     TRANSFORM H SUB I
     
     CALL gmmatd( habc(6*i-5),2,3,0,  t(1),3,3,0,  temp9(1) )
     
     
     npoint = km + i
     npoint = 9*m(npoint) - 9
     
!     J = 1    ROW 1 OF H INTO ROW 1 OF G.
!              ROW 2 OF H INTO ROW 2 OF G.
!     J = 2    ROW 1 OF H INTO ROW 2 OF G.
!              ROW 2 OF H INTO ROW 3 OF G.
!     J = 3    ROW 1 OF H INTO ROW 3 OF G.
!              ROW 2 OF H INTO ROW 1 OF G.
     
     IF( j-2  < 0) THEN
       GO TO   170
     ELSE IF ( j-2  == 0) THEN
       GO TO   160
     ELSE
       GO TO   190
     END IF
     
     160 npoint = npoint + 3
     170 DO  k=1,6
       npoint = npoint + 1
       g(npoint) = g(npoint) + temp9(k)
     END DO
     CYCLE
     190 g(npoint + 7) = g(npoint + 7) + temp9(1)
     g(npoint + 8) = g(npoint + 8) + temp9(2)
     g(npoint + 9) = g(npoint + 9) + temp9(3)
     g(npoint + 1) = g(npoint + 1) + temp9(4)
     g(npoint + 2) = g(npoint + 2) + temp9(5)
     g(npoint + 3) = g(npoint + 3) + temp9(6)
     
   END DO
   
   
!     END OF LOOP FOR BASIC TRIANGLES
   
 END DO
!     ******************************************************************
 
!     FILL E-MATRIX
 
 DO  i=1,18
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
 
 IF( necpt(4*npivot+9) == 0 ) GO TO 230
 CALL transd( necpt(4*npivot+9), t(1) )
 CALL gmmatd( t(1),3,3,1,  e( 1),3,3,0,  tite( 1) )
 CALL gmmatd( t(1),3,3,1,  e(10),3,3,0,  tite(10) )
 GO TO 250
 230 DO  k=1,18
   tite(k) = e(k)
 END DO
 
!     SOLVE NOW FOR ....
 
!    E                   T     T                       T
! (K  ) = (K  ) - (TERM ) (K  ) - (K  )(TERM ) + (TERM )(K  )(TERM )
!   IJ      IJ         I    J4      I4      J         I   44      J
 
!                        -1                               I=NPIVOT
! WHERE... (TERM ) = (G )  (G ) ,I=NPIVOT                 J=1,2,3
!               I      4     I
 
!                        -1
!          (TERM ) = (G )  (G ) ,J=1,2,3 AS ABOVE
!               J      4     J
 
!     AND WITH TRANSFORMATIONS....
 
!    G        T      E   T
! (K  ) = (C ) (E)(K  )(E )(C )
!   IJ      I       IJ       J
 
 
!     COMPUTE  (TERM        )  STORE IN PROD9
!                   I=NPIVOT
 
!                   -1
!     FIRST GET (G )
!                 4
 
 250 CONTINUE
!     NO NEED TO COMPUTE DETERMINANT SINCE IT IS NOT USED SUBSEQUENTLY.
 ising = -1
 CALL inverd( 3,g(28),3,prod9,0,determ,ising,temp9 )
 
!     CHECK FOR SINGULARITY. ISING=2 IMPLIES SINGULARITY.
 SELECT CASE ( ising )
   CASE (    1)
     GO TO 270
   CASE (    2)
     GO TO 260
 END SELECT
 260 CALL mesage(30,36,ecpt(1))
 
!  SET FLAG FOR FATAL ERROR WHILE ALLOWING ERROR MESSAGES TO ACCUMULATE
 
 nogo=1
 RETURN
 
 270 CALL gmmatd ( g(28),3,3,0,  g(9*npivot-8),3,3,0,  prod9(1) )
 
!                       T
!     GET  (TERM        )(K  )  STORE IN TEMP9
!               I=NPIVOT   44
 
 CALL gmmatd( prod9(1),3,3,1,  ksum(55),3,3,0,  temp9(1) )
 
 
 
!     THE TWO COMMON PRODUCTS ARE NOW AT HAND IN PROD9 AND TEMP9.
 
 DO  j=1,3
   
!                   T     T
!     (TERM        ) (K  )    STORE IN ARR9
!          I=NPIVOT    J4
   
   CALL gmmatd( prod9(1),3,3,1,  ksum(9*j+19),3,3,1,  arr9(1) )
   
!     SUBTRACT FROM (K  )
!                     IJ
   
   nbegin = 9*j-9
   DO  i=1,9
     npoint = nbegin + i
     ksum(npoint) = ksum(npoint) - arr9(i)
   END DO
   
   
!      COMPUTE  (TERM )  STORE IN ARR9
!                   J
   
   CALL gmmatd( g(28),3,3,0,  g(9*j-8),3,3,0,  arr9(1) )
   
   
!     GET  (K  )(TERM )  STORE IN ARRAY9
!            I4      J
   
   CALL gmmatd( ksum(9*npivot+19),3,3,0,  arr9(1),3,3,0,  array9(1))
   
!     SUBTRACT FROM KIJ
   
   DO  i=1,9
     npoint = nbegin + i
     ksum(npoint) = ksum(npoint) - array9(i)
   END DO
   
!                           T
!     COMPUTE  (TERM        )(K  )(TERM ) = (TEMP9)(ARR9)
!                   I=NPOINT   44      J
   
   CALL gmmatd( temp9(1),3,3,0,  arr9(1),3,3,0,  array9(1) )
   
!     ADD TO K
!             IJ
   
   DO  i=1,9
     npoint = nbegin + i
     ksum(npoint) = ksum(npoint) + array9(i)
   END DO
   
!       E
!     K    COMPLETE
!      IJ
   
!     TRANSFORM NOW, AND INSERT.
   
   
!     TRANSFORMATIONS AND INSERTION
   
   IF( necpt(4*j+9) == 0) GO TO 310
   CALL transd( necpt(4*j+9), t(1) )
   CALL gmmatd( t(1),3,3,1,  e( 1),3,3,0,  tjte( 1) )
   CALL gmmatd( t(1),3,3,1,  e(10),3,3,0,  tjte(10) )
   GO TO 330
   310 DO  k=1,18
     tjte(k) = e(k)
   END DO
   330 CALL gmmatd( ksum(nbegin+1),3,3,0,  tjte(1),6,3,1,  temp18(1) )
   CALL gmmatd ( tite(1),6,3,0,  temp18(1),3,6,0,  kout(1))
   CALL sma1b( kout(1), necpt(j+1), -1, ifkgg, 0.0D0 )
   temp = gsube
   IF( iopt4  == 0) THEN
     GO TO   360
   END IF
   340 IF( gsube  == 0.0) THEN
     GO TO   360
   END IF
   350 CALL sma1b( kout(1), necpt(j+1), -1, if4gg, temp )
   k4ggsw = 1
   
 END DO
 RETURN
 370 CALL mesage(30,26,ecpt(1))
!  SET FLAG FOR FATAL ERROR WHILE ALLOWING ERROR MESSAGES TO ACCUMULATE
 
!  NOTE - NON-FATAL MESSAGE HERE MAY INDUCE PARTITION ERROR 3111 LATER
!         IN EMGOUT
 
 nogo=1
 RETURN
END SUBROUTINE ktrplt
