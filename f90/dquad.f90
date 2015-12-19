SUBROUTINE dquad (itype)
     
!     THIS ROUTINE GENERATES THE FOLLOWING
 
!     FOUR 6X6 DIFFERENTIAL STIFFNESS MATRICES FOR ONE PIVOT POINT OF
!     A QUADRILATERAL
 
 
!     CALLS FROM THIS ROUTINE ARE MADE TO
!           DTRBSC - BASIC BENDING TRI. ROUTINE.
!           DTRMEM - TRIANGULAR MEMBRANE ROUTINE
!           TRANSD - SUPPLIES 3X3 TRANSFORMATIONS
!           GMMATD - GENERAL MATRIX MULITPLY AND TRANSPOSE ROUTINE
!           DS1B   - INSERTION ROUTINE
 
 
!        ITYPE    = 1             2                       4
!     ECPT INDEX    QUAD1         QUAD2        TRMEM      QUAD4
!     **********    *******       *******      *******    ********
!          1        EL. ID.       EL. ID.      EL. ID.    EL. ID
!          2        SIL1          SIL1         SIL1       SIL1
!          3        SIL2          SIL2         SIL2       SIL2
!          4        SIL3          SIL3         SIL3       SIL3
!          5        SIL4          SIL4         THETA      SIL4
!          6        THETA         THETA        MAT. ID.   MEM.T1
!          7        MAT. ID. 1    MAT. ID.     T          MEM.T2
!          8        T1            T            NSM        MEM.T3
!          9        MAT. ID. 2    NSM          CID1       MEM.T4
!         10        INERTIA I     CID1         X1         THETA
!         11        MAT ID  3     X1           Y1         FLAG FOR 10
!         12        T2            Y1           Z1         GRD OFFSET
!         13        NSM           Z1           CID2       MAT. ID 1
!         14        Z1            CID2         X2         THICKNESS
!         15        Z2            X2           Y2         MAT. ID 2
!         16        CID1          Y2           Z2         INERTIA I
!         17        X1            Z2           CID3       MAT. ID 3
!         18        Y1            CID3         X3         TS/T
!         19        Z1            X3           Y3         NSM
!         20        CID2          Y3           Z3         Z1
!         21        X2            Z3           EL TEMP    Z2
!         22        Y2            CID4         EL DEFORM  MAT. ID 4
!         23        Z2            X4           LOAD TEMP  THETA
!         24        CID3          Y4           U1         FLAG FOR 23
!         25        X3            Z4           V1         INTEGRATION
!         26        Y3            EL TEMP      W1         STRESS ANGLE
!         27        Z3            EL DEFORM    U2         FLAG FOR 26
!         28        CID4          LOAD TEMP    V2         ZOFF1
!         29        X4            U1           W2         CID1
!         30        Y4            V1           U3         X1
!         31        Z4            W1           V3         Y1
!         32        EL TEMP       U2           W3         Z1
!         33        EL DEFORM     V2                      CID2
!         34        LOAD TEMP     W2                      X2
!         35        U1            U3                      Y2
!         36        V1            V3                      Z2
!         37        W1            W3                      CID3
!         38        U2            U4                      X3
!         39        V2            V4                      Y3
!         40        W2            W4                      Z3
!         41        U3                                    CID4
!         42        V3                                    X4
!         43        W3                                    Y4
!         44        U4                                    Z4
!         45        V4                                    EL TEMP
!         46        W4
!         47
!         48                                              U1
!         49                                              V1
!         50                                              W1
!         51                                              U2
!         52                                              V2
!         53                                              W2
!         54                                              U3
!         55                                              V3
!         56                                              W3
!         57                                              U4
!         58                                              V4
!         59                                              W4
 
 
 INTEGER, INTENT(IN OUT)                  :: itype
 INTEGER :: subsca        ,subscb        ,subscc
 DOUBLE PRECISION :: kout          ,tite          ,dpdum        ,  &
     tjte          ,dpdum2        ,ivect        ,  &
     d1            ,jvect         ,d2           ,  &
     kvect         ,a1            ,ksum         ,  &
     t             ,xsubb         ,v            ,  &
     xsubc         ,vv            ,ysubc        ,  &
     prod9         ,temp          ,temp9        ,  &
     u1            ,h             ,u2           ,  &
     e             ,a             ,temp18       ,  &
     requiv        ,r             ,sigxy        , sigx          ,sigy
 DIMENSION       m(12)         ,necpt(100)    ,requiv(8)    ,  &
     vq1(3),vq2(3) ,vq3(3),vq4(3) ,a(1)
 CHARACTER (LEN=25) :: sfm
 CHARACTER (LEN=29) :: uim
 CHARACTER (LEN=25) :: uwm
 CHARACTER (LEN=23) :: ufm
 COMMON /xmssg / ufm           ,uwm           ,uim          , sfm
 COMMON /condas/ consts(5)
 COMMON /system/ ibuff         ,nout          ,nogo
 COMMON /matin / matid         ,inflag        ,eltemp       ,  &
     stress        ,sinth         ,costh
 COMMON /matout/ g11           ,g12           ,g13          ,  &
     g22           ,g23           ,g33          ,  &
     rho           ,alpha1        ,alpha2       ,  &
     alp12         ,t_sub_0       ,g_sub_e      ,  &
     sigten        ,sigcom        ,sigshe       ,  &
     g2x211        ,g2x212        ,g2x222
 COMMON /ds1aaa/ npvt          ,icstm         ,ncstm
 COMMON /ds1aet/ ecpt(100)
 COMMON /ds1adp/ kout(36)      ,tite(18)      ,tjte(18)     ,  &
     temp18(18)    ,d1(3)         ,d2(3)        ,  &
     a1(3)         ,v(2)          ,vv(2)        ,  &
     prod9(9)      ,temp9(9)      ,h            ,  &
     u1            ,u2            ,dpdum(1)     ,  &
     temp          ,dpdum2(43)    ,e(18)        ,  &
     sigx          ,sigy          ,sigxy        ,  &
     xsubb         ,xsubc         ,ysubc        ,  &
     ksum(36)      ,t(9)          ,ivect(3)     ,  &
     jvect(3)      ,kvect(3)      ,r(2,4)       ,  &
     sp1(2)        ,theta         ,sinang       ,  &
     cosang        ,km            ,nbegin       ,  &
     jnot          ,npivot        ,nsubc        ,  &
     ising         ,subsca        ,subscb       ,  &
     subscc        ,npoint        ,ipvt
 EQUIVALENCE     (consts(4),degra) , (necpt(1),ecpt(1))     ,  &
     (requiv(1),r(1,1)), (vq1(1),ecpt(17))      ,  &
     (vq2(1),ecpt(21)) , (vq3(1),ecpt(25))      ,  &
     (vq4(1),ecpt(29)) , (a(1),kout(1))
 DATA     m   /  2, 4, 1,   3, 1, 2,   4, 2, 3,   1, 3, 4   /
 
 
!     IF ITYPE = 2, QUAD2 EST DATA IS MOVED AND STORED IN QUAD1 FORMAT
!     IF ITYPE = 4, QUAD4 EST DATA IS MOVED AND STORED IN QUAD1 FORMAT
 
 IF (itype == 4) GO TO 15
 IF (itype /= 2) GO TO 20
 
 DO  i = 10,40
   npoint = 50 - i
   ecpt(npoint+6) = ecpt(npoint)
 END DO
 
 ecpt( 9) = ecpt(7)
 ecpt(10) =(ecpt(8)**3.0)/12.0
 ecpt(11) = ecpt(7)
 ecpt(12) = ecpt(8)
 GO TO 20
 
!     QUAD4
 
!     IF NECPT(11)=0, ECPT(10) IS THE MATERIAL PROPERTY ORIENTAION
!     ANGLE THETA. IF IT IS NOT, NECPT(10) IS MATERIAL COORDINATE
!     SYSTEM ID. IN THIS CASE, WE CAN NOT CONTINUE
 
 15 IF (necpt(11) /= 0) GO TO 350
 ecpt(6) = ecpt(10)
 ecpt(7) = ecpt(13)
 ecpt(8) = ecpt(14)
 ecpt(9) = ecpt(15)
 ecpt(10)= ecpt(16)
 ecpt(11)= ecpt(17)
 ecpt(12)= ecpt(14)
 DO  i = 16,46
   ecpt(i) = ecpt(i+13)
 END DO
 20 IF (ecpt(8) == 0.0) RETURN
 
!     CALL BUG (4HQDET,5,ECPT,52-6*ITYPE)
 
!     DETERMINE PIVOT POINT NUMBER
 
 DO  i = 1,4
   IF (npvt /= necpt(i+1)) CYCLE
   npivot = i
   GO TO 40
 END DO
 RETURN
 
 40 theta  = ecpt(6)*degra
 sinang = SIN(theta)
 cosang = COS(theta)
 
 IF (npivot-2 > 0) THEN
   GO TO    60
 END IF
 50 jnot = npivot + 2
 GO TO 70
 60 jnot = npivot - 2
 
!     FORMATION OF THE R-MATRIX CONTAINING COORDINATES OF THE
!     SUB TRIANGLES.  (2X4) FOR QUADRILATERAL PLATE...
!     FORMATION ALSO OF THE I,J, AND K VECTORS USED IN THE E-MATRIX.
 
!     ZERO OUT R-MATRIX
 
 70 DO  i = 1,8
   requiv(i) = 0.0D0
 END DO
 
 DO  i = 1,3
   d1(i) = DBLE(vq3(i)) - DBLE(vq1(i))
   d2(i) = DBLE(vq4(i)) - DBLE(vq2(i))
   a1(i) = DBLE(vq2(i)) - DBLE(vq1(i))
 END DO
 
!     NON-NORMALIZED K-VECTOR = D1 CROSS D2
 
 kvect(1) = d1(2)*d2(3) - d2(2)*d1(3)
 kvect(2) = d1(3)*d2(1) - d2(3)*d1(1)
 kvect(3) = d1(1)*d2(2) - d2(1)*d1(2)
 
!     NORMALIZE K-VECTOR
 
 temp = DSQRT(kvect(1)**2 + kvect(2)**2 + kvect(3)**2)
 IF (temp == 0.0D0) CALL mesage (-30,26,ecpt(1))
 DO  i = 1,3
   kvect(i) = kvect(i)/temp
 END DO
 
!     COMPUTE H = (A1 DOT KVECT) / 2
 
 temp = (a1(1)*kvect(1) + a1(2)*kvect(2) + a1(3)*kvect(3))/2.0D0
 
!     I-VECTOR =(A1) - H*(KVECT)    NON-NORMALIZED
 
 DO  i = 1,3
   ivect(i) = a1(i) - temp*kvect(i)
 END DO
 
!     NORMALIZE I-VECTOR
 
 temp = DSQRT(ivect(1)**2 + ivect(2)**2 + ivect(3)**2)
 IF (temp == 0.0D0) CALL mesage (-30,26,ecpt(1))
 DO  i = 1,3
   ivect(i) = ivect(i)/temp
 END DO
 
!     J-VECTOR = K CROSS I, AND X3 CALCULATION
 
 jvect(1) = kvect(2)*ivect(3) - ivect(2)*kvect(3)
 jvect(2) = kvect(3)*ivect(1) - ivect(3)*kvect(1)
 jvect(3) = kvect(1)*ivect(2) - ivect(1)*kvect(2)
 
!     NORMALIZE J VECTOR TO MAKE SURE
 
 temp =  DSQRT(jvect(1)**2 + jvect(2)**2 + jvect(3)**2)
 IF (temp == 0.0D0) CALL mesage (-30,26,ecpt(1))
 DO  i = 1,3
   jvect(i) = jvect(i)/temp
 END DO
 
!     X3 GOES INTO R(1,3) = D1 DOT IVECT
 
 r(1,3) = d1(1)*ivect(1) + d1(2)*ivect(2) + d1(3)*ivect(3)
 
!     X2 GOES INTO R(1,2) AND Y3 GOES INTO R(2,3)
 
 r(1,2) = a1(1)*ivect(1) + a1(2)*ivect(2) + a1(3)*ivect(3)
 r(2,3) = d1(1)*jvect(1) + d1(2)*jvect(2) + d1(3)*jvect(3)
 
!     X4 GOES INTO R(1,4) AND Y4 GOES INTO R(2,4)
 
 r(1,4) = d2(1)*ivect(1) + d2(2)*ivect(2) + d2(3)*ivect(3) + r(1,2)
 r(2,4) = d2(1)*jvect(1) + d2(2)*jvect(2) + d2(3)*jvect(3)
 
!     AT THIS POINT, THE COORDINATES OF THE PLATE IN THE ELEMENT
!     SYSTEM ARE STORED IN THE R-MATRIX WHERE THE COLUMN DENOTES THE
!     POINT AND THE ROW DENOTES THE X OR Y COORDINATE FOR ROW 1 OR
!     ROW 2 RESPECTIVELY.
 
!     SET UP THE M-MATRIX FOR MAPPING TRIANGLES, IN DATA STATEMENT.
 
!     COMPUTE SUB-TRIANGLE COORDINATES
 
!     ZERO OUT KSUM MATRICES
 
 DO  i = 1,36
   ksum(i) = 0.0D0
 END DO
 
 eltemp = ecpt(32)
 
!     MOVE ECPT INTO POSITIONS 51-93
 
 DO  i = 1,46
   ecpt(i+50) = ecpt(i)
 END DO
 
!     MOVE MISCELLANEOUS VARIABLES INTO TRMEM FORMAT
 
 ecpt( 6) = ecpt( 7)
 ecpt( 7) = ecpt( 8)
 ecpt(21) = ecpt(32)
 ecpt(22) = ecpt(33)
 ecpt(23) = ecpt(34)
 
 DO  j = 1,4
   IF (j == jnot) CYCLE
   km   = 3*j - 3
   ipvt = 0
   DO  i = 1,3
     npoint = km+i
     nsubc  = m(npoint)
     IF (nsubc == npivot) ipvt = i
     necpt(i+1) = necpt(nsubc+51)
     DO  k = 1,4
       npoint = 4*(nsubc-1) + k + 65
       subsca = 4*(i-1) + k + 8
       ecpt(subsca) = ecpt(npoint)
     END DO
     DO  k = 1,3
       npoint = 3*(nsubc-1) + k + 84
       subsca = 3*(i-1) + k + 23
       ecpt(subsca) = ecpt(npoint)
     END DO
   END DO
   IF (ipvt == 0) CYCLE
   
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
   
   sinth = sinang*u1 - cosang*u2
   costh = cosang*u1 + sinang*u2
   IF (ABS(sinth) < 1.0E-06) sinth = 0.0
   
!     AT THIS POINT, XSUBB, XSUBC, YSUBC ARE AT HAND FOR TRIANGLE -J-
   
   CALL dtrmem (3)
   CALL dtrbsc (2,ipvt)
   
!     NOW WE HAVE AT HAND  K   I=NPIVOT,J=1,2,3   THREE 6X6 MATRICES
!                           IJ
!                                STORED AT  A(1) THROUGH A(27)
   
!     MAP THE THE 3X3 S FOR THE PIVOT ROW INTO THE SUMMATION ARRAYS
   
   DO  i = 1,3
     npoint = 9*i - 8
     
     CALL gmmatd (t,3,3,1, a(npoint),3,3,0, temp9)
     CALL gmmatd (temp9,3,3,0, t,3,3,0, prod9)
     
!     ADD THIS PRODUCT IN NOW.
     
     npoint = km + i
     npoint = 9*m(npoint) - 9
     DO  k = 1,9
       npoint = npoint + 1
       ksum(npoint) = ksum(npoint) + prod9(k)/2.0D0
     END DO
   END DO
   
 END DO
 
!     CALL BUG (4HQDKD,220,KSUM,72)
 
!     FILL E-MATRIX
 
 DO  i = 1,18
   e(i)  = 0.0D0
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
 
 IF (necpt(4*npivot + 62) == 0) GO TO 260
 CALL transd (necpt(4*npivot+62),t)
 CALL gmmatd (t,3,3,1, e( 1),3,3,0, tite( 1))
 CALL gmmatd (t,3,3,1, e(10),3,3,0, tite(10))
 GO TO 290
 
 260 DO  k = 1,18
   tite(k) = e(k)
 END DO
 
!     RESTORE ECPT FOR CKECKOUT
 
 DO  k = 1,46
   ecpt(k) = ecpt(k+50)
 END DO
 
 290 DO  j = 1,4
   
!     TRANSFORMATIONS AND INSERTION
   
   IF (necpt(4*j+62) == 0) GO TO 300
   CALL transd (necpt(4*j+62),t)
   CALL gmmatd (t,3,3,1, e(1),3,3,0,  tjte( 1))
   CALL gmmatd (t,3,3,1, e(10),3,3,0, tjte(10))
   GO TO 320
   
   300 DO  k = 1,18
     tjte(k) = e(k)
   END DO
   320 CALL gmmatd (ksum(9*j-8),3,3,0, tjte,6,3,1, temp18(1))
   CALL gmmatd (tite(1),6,3,0, temp18(1),3,6,0, kout(1))
   CALL ds1b (kout,necpt(j+51))
 END DO
 RETURN
 
!     COULD NOT CONTINUE
 
 350 WRITE  (nout,360) sfm
 360 FORMAT (a25,', DEFFICIENT SOURCE CODE IN DQUAD TO HANDLE CQUAD4 ',  &
     'ELEMENT WITH MATERIAL', /5X,  &
     'PROPERTY COORD. SYSTEM. ANGLE MUST BE SPECIFIED')
 nogo = 1
 RETURN
END SUBROUTINE dquad
