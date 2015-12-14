SUBROUTINE ektrms (ntype)
     
!     TRIANGULAR MEMBRANE ELEMENT
 
!     ECPT LIST
!                                                      THIS
!       ECPT       DESCRIPTION                         ROUTINE   TYPE
!     ******************************************************************
!       ECPT( 1) = ELEMENT ID                          NECPT(1)  INTEGER
!       ECPT( 2) = GRID POINT A                        NGRID(1)  INTEGER
!       ECPT( 3) = GRID POINT B                        NGRID(2)  INTEGER
!       ECPT( 4) = GRID POINT C                        NGRID(3)  INTEGER
!       ECPT( 5) = THETA = ANGLE OF MATERIAL           ANGLE     REAL
!       ECPT( 6) = MATERIAL ID                         MATID     INTEGER
!       ECPT( 7) = T                                   T         REAL
!       ECPT( 8) = NON-STRUCTURAL MASS                 FMU       REAL
!       ECPT( 9) = COORD. SYSTEM ID 1                  NECPT(9)  INTEGER
!       ECPT(10) = X1                                  X1        REAL
!       ECPT(11) = Y1                                  Y1        REAL
!       ECPT(12) = Z1                                  Z1        REAL
!       ECPT(13) = COORD. SYSTEM ID 2                  NECPT(13) INTEGER
!       ECPT(14) = X2                                  X2        REAL
!       ECPT(15) = Y2                                  Y2        REAL
!       ECPT(16) = Z2                                  Z2        REAL
!       ECPT(17) = COORD. SYSTEM ID 3                  NECPT(17) INTEGER
!       ECPT(18) = X3                                  X3        REAL
!       ECPT(19) = Y3                                  Y3        REAL
!       ECPT(20) = Z3                                  Z3        REAL
!       ECPT(21) = ELEMENT TEMPERATURE                 ELTEMP    REAL
 
!     IF NTYPE = 0  COMPLETE MEMBRANE COMPUTATION IS PERFORMED
!     IF NTYPE = 1  9  3X3 MATRICES FOR THE GRID POINTS  IN ECPT
 
 
 INTEGER, INTENT(IN OUT)                  :: ntype
 LOGICAL :: nogo,heat
 INTEGER :: necpt(6)
 REAL :: k,kout,ecpt(21),matbuf,kij,g(9),c(18),tt(2),ti(9), tempar(27)
 COMMON /emgtrx/ a(225),prod9(9),temp9(9),xsubb,xsubc,ysubc,dict5,  &
     e(18),k(324),kout(324),kij(81)
 COMMON /system/ ksystm(60)
 COMMON /condas/ consts(5)
 COMMON /emgprm/ dum(19),nogo,heat
 COMMON /emgest/ mecpt(1),ngrid(3),angle,matid1,t,fmu,dummy1,x1,y1,  &
     z1,dummy2,x2,y2,z2,dummy3,x3,y3,z3,dumb(80)
 COMMON /matin / matid,inflag,eltemp,stress,sinth,costh
 COMMON /matout/ g11,g12,g13,g22,g23,g33,rho,alpha1,alpha2,alp12,  &
     t_sub_0,g sub e,sigten,sigcom,sigshe, g2x211,g2x212,g2x222
 COMMON /hmtout/ matbuf(4)
 EQUIVALENCE     (consts(4),degra),(ecpt(1),necpt(1),mecpt(1)),  &
     (ksystm(2),ioutpt),(ksystm(56),iheat)
 
!     SET UP THE E MATRIX WHICH IS (3X2) FOR THE TRI-MEMBRANE
 
!     E(1), E(3), E(5) WILL BE THE I-VECTOR
!     E(2), E(4), E(6) WILL BE THE J-VECTOR
!     E(7), E(8), E(9) WILL BE THE K-VECTOR NOT USED IN E FOR MEMBRANE
 
!     FIRST FIND I-VECTOR = RSUBB - RSUBA  (NON-NORMALIZED)
 
 e(1) = x2 - x1
 e(3) = y2 - y1
 e(5) = z2 - z1
 
!     NOW FIND LENGTH = X-SUB-B   COORD. IN ELEMENT SYSTEM
 
 xsubb = SQRT(e(1)**2 + e(3)**2 + e(5)**2)
 IF (xsubb <= 1.e-06)  GO TO 7770
 
!  20 NOW NORMALIZE I-VECTOR WITH X-SUB-B
 
 e(1) = e(1)/xsubb
 e(3) = e(3)/xsubb
 e(5) = e(5)/xsubb
 
!     HERE WE NOW TAKE RSUBC - RSUBA AND STORE TEMPORARILY IN
!     E(2), E(4), E(6) WHICH IS WHERE THE J-VECTOR WILL FIT LATER
 
 e(2) = x3 - x1
 e(4) = y3 - y1
 e(6) = z3 - z1
 
!     X-SUB-C  =  I . (RSUBC - RSUBA), THUS
 
 xsubc = e(1)*e(2) + e(3)*e(4) + e(5)*e(6)
 
!     AND CROSSING THE I-VECTOR TO (RSUBC-RSUBA) GIVES THE K-VECTOR
!     (NON-NORMALIZED)
 
 e(7) = e(3)*e(6) - e(5)*e(4)
 e(8) = e(5)*e(2) - e(1)*e(6)
 e(9) = e(1)*e(4) - e(3)*e(2)
 
!     THE LENGTH OF THE K-VECTOR IS NOW FOUND AND EQUALS Y-SUB-C
!     COORD. IN ELEMENT SYSTEM
 
 ysubc = SQRT(e(7)**2 + e(8)**2 + e(9)**2)
 IF (ysubc <= 1.e-06) GO TO 7780
 
!  25 NOW NORMALIZE K-VECTOR WITH YSUBC JUST FOUND
 
 e(7) = e(7)/ysubc
 e(8) = e(8)/ysubc
 e(9) = e(9)/ysubc
 
!     J VECTOR = K CROSS I
!     STORE IN THE SPOT FOR J
 
 e(2) = e(5)*e(8) - e(3)*e(9)
 e(4) = e(1)*e(9) - e(5)*e(7)
 e(6) = e(3)*e(7) - e(1)*e(8)
 
!     AND JUST FOR COMPUTER EXACTNESS NORMALIZE J-VECTOR TO MAKE SURE.
 
 temp = SQRT(e(2)**2 + e(4)**2 + e(6)**2)
 IF (temp == 0.0) GO TO 7790
 e(2) = e(2)/temp
 e(4) = e(4)/temp
 e(6) = e(6)/temp
 
!     VOLUME OF ELEMENT, THETA, MU, LAMDA, AND DELTA
 
 vol    = xsubb*ysubc*t/2.
 reelmu = 1./xsubb
 flamda = 1./ysubc
 delta  = xsubc/xsubb - 1.
 
!     NOW FORM THE  C MATRIX   (3X6) PARTITIONED AS FOLLOWS HERE.
!         CSUBA = (3X2) STORED IN C( 1) THRU C( 6) BY ROWS
!         CSUBB = (3X2) STORED IN C( 7) THRU C(12) BY ROWS
!         CSUBC = (3X2) STORED IN C(13) THRU C(18) BY ROWS
 
 c(1)  =-reelmu
 c(2)  = 0.
 c(3)  = 0.
 c(4)  = flamda*delta
 c(5)  = c(4)
 c(6)  =-reelmu
 c(7)  = reelmu
 c(8)  = 0.
 c(9)  = 0.
 c(10) =-flamda*reelmu*xsubc
 c(11) = c(10)
 c(12) = reelmu
 c(13) = 0.
 c(14) = 0.
 c(15) = 0.
 c(16) = flamda
 c(17) = flamda
 c(18) = 0.
 
 IF (ntype == 1) GO TO 30
 
 theta = angle*degra
 sinth = SIN(theta)
 costh = COS(theta)
 30 IF (ABS(sinth) < 1.0E-06) sinth = 0.0
 
!     BRANCH ON -HEAT- PROBLEM AT THIS POINT.
 
 IF (heat) GO TO 300
 eltemp = ecpt(21)
 matid  = matid1
 inflag = 2
 CALL mat (ecpt(1))
 
!     FILL G-MATRIX WITH OUTPUT FROM MAT ROUTINE
 
 g(1) = g11
 g(2) = g12
 g(3) = g13
 g(4) = g12
 g(5) = g22
 g(6) = g23
 g(7) = g13
 g(8) = g23
 g(9) = g33
 
!     AT 50 G, E, AND C MATRICES ARE COMPLETE
 
!     AT THIS POINT THE FOLLOWING EQUATION CAN BE SOLVED FOR K-SUB-IJ
 
!                     T        T             T
!       K   = VOL . T  * E * C  * G * C  * E  * T
!        IJ          I        I        J         J
 
!     T-SUB-I WILL BE USED IN THE ABOVE ONLY IF THE PIVOT COORDINATE
!     SYSTEM ID IS NOT ZERO, OTHERWISE IT IS ASSUMED TO BE THE
!     IDENTITY MATRIX.
 
!     THE I SUBSCRIPT IMPLIES THE PIVOT POINT  1,2, OR 3 (ELEMENT SYST)
!     THE J SUBSCRIPT IMPLIES  1 THRU 3  FOR EACH CALL TO THIS ROUTINE.
 
!     DO COMPUTATIONS FOR EACH POINT IN ECPT LIST
 
 DO  i = 1,81
   kij(i) = 0.
 END DO
 DO  npvt = 1,3
   ka = 4*npvt + 5
   npoint = 6*npvt - 5
   
!                     T
!     COMPUTE   E * C   * G       AND STORE IN TEMPAR(1 THRU 9)
!                    I
   
   CALL gmmats (e,3,2,0, c(npoint),3,2,1, tempar(10))
   CALL gmmats (tempar(10),3,3,0, g,3,3,0, tempar(1))
   
!     NCOM WILL ALWAYS POINT TO THE COMMON 3 X 3 PRODUCT ABOVE
!     NPT1 WILL POINT TO FREE WORKING SPACE LENGTH 9
   
   ncom = 1
   npt1 = 10
   
!     MULTIPLY COMMON PRODUCT BY SCALER VOL
   
   DO  i = 1,9
     tempar(i) = tempar(i)*vol
   END DO
   
!     CHECK FOR PIVOT  CSID = 0,  IF ZERO SKIP TRANSFORMATION TSUBI.
   
   IF (necpt(ka) == 0) GO TO 110
   
!     NOT-ZERO THUS GET TI
   
   CALL transs (necpt(ka),ti)
   
!     INTRODUCE TI INTO THE COMMON PRODUCT AND STORE AT
!     TEMPAR(10 THRU 18)
   
   CALL gmmats (ti,3,3,1, tempar(1),3,3,0, tempar(10))
   
!     COMMON PRODUCT NOW STARTS AT TEMPAR(10) THUS CHANGE NCOM AND NPT1
   
   ncom = 10
   npt1 =  1
   
!  80 NOW HAVE COMMON PRODUCT STORED BEGINNING TEMPAR(NCOM),  (3X3).
!     NPT1 POINTS TO FREE WORKING SPACE LENGTH 9.
   
!     PROCEED NOW AND RUN OUT THE 3 6X6 MATRICES KIJ-SUB-1,2,3.
   
   110 nsave  = npt1
   npoint = (npvt-1)*27
   
!     INSERT G INTO TEMPAR
   
   DO  i = 1,9
     tempar(i+18) = g(i)
   END DO
   DO  i = 1,3
     CALL gmmats (c(6*i-5),3,2,0, e,3,2,1, tempar(nsave))
     
!     NPT2 IS SET TO POINT TO THE BEGINNING OF THE PRODUCT  C * E * T
!                                                            J       J
     npt2 = nsave
     npt1 = 19
     
!     CHECK FOR ZERO CSID IN WHICH CASE TJ IS NOT NEEDED
     
     IF (necpt(4*i +5) == 0) GO TO 120
     
!     COMMING HERE IMPLIES NEED FOR TJ
!     WILL STORE TJ IN TI
     
     CALL transs (necpt(4*i+5),ti)
     CALL gmmats (tempar(npt2),3,3,0, ti,3,3,0, tempar(19))
     npt1 = npt2
     npt2 = 19
     
!  60 AT THIS POINT COMPLETE COMPUTATION FOR  K-SUB-I,J
     
     120 CALL gmmats (tempar(ncom),3,3,0, tempar(npt2),3,3,0, tempar(npt1))
     npt36 = npt1 + 35
     
     DO  j = 1,9
       npoint = npoint + 1
       npt2   = npt1 + j - 1
       kij(npoint) = tempar(npt2)
     END DO
   END DO
 END DO
 
 dict5 = gsube
 RETURN
 
!     HEAT PROBLEM LOGIC PICKS UP HERE.  CALL HMAT FOR MATERIAL DATA.
 
 300 inflag = 2
 matid  = necpt(6)
 eltemp = ecpt(21)
 CALL hmat (necpt)
 g(1) = matbuf(1)
 g(2) = matbuf(2)
 g(3) = matbuf(2)
 g(4) = matbuf(3)
 
!     CONDENSE C MATRIX FOR HEAT PROBLEM (FORMED ABOVE)  C IS (2X3)
 
 c(2) = c(4)
 c(3) = c(7)
 c(4) = c(10)
 c(5) = c(13)
 c(6) = c(16)
 
!     DETERMINE THE PIVOT POINT.
 
 kq   = 3
 kmax = kq*3
 DO  i = 1,kmax
   kij(i) = 0.
 END DO
 DO  npvt = 1,3
   
!     PIVOT C MATRIX TIMES VOLUME (STORED INTO TT(1) AND TT(2).)
   
   tt(1) = vol*c(2*npvt-1)
   tt(2) = vol*c(2*npvt  )
   
!     OUTPUT THE CONDUCTIVITY MATRICES
   
   npoint = (npvt-1)*kq
   
   DO  i = 1,3
     n2 = 2*i
     n1 = n2 - 1
     tempar(1) = (g(1)*c(n1) + g(2)*c(n2))*tt(1)  +  &
         (g(3)*c(n1) + g(4)*c(n2))*tt(2)
     
!     SUB-TRIANGLE (RETURN 3X3-S AS ABOVE IN STIFFNESS PORTION)
     
     kij(npoint+1) = tempar(1)
     npoint = npoint + 1
   END DO
 END DO
 RETURN
 
!     ERROR  EXITS
 
 7770 CALL mesage (30,31,necpt(1))
 7777 nogo = .true.
 RETURN
 7780 CALL mesage (30,32,necpt(1))
 GO TO 7777
 7790 CALL mesage (30,26,necpt(1))
 GO TO 7777
 
END SUBROUTINE ektrms
