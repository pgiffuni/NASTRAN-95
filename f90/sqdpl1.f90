SUBROUTINE sqdpl1
     
!     PHASE I OF STRESS DATA RECOVERY FOR TRI OR QUAD PLATE.
 
!     OUTPUTS FROM THIS PHASE FOR USE IN PHASE II ARE THE FOLLOWING.
 
!     1) ELEMENT ID
!     2) 4 SILS
!     3) I
!     4) Z1 AND Z2
!     5) 4  5X6 S-SUB-I ARRAYS
!     6) 3X1 S SUB T MATRIX
!     THUS, 131 WORDS FOR QUAD-PLATE
 
!     ECPT LISTS AS OF AUGUST 4, 1967
 
!                 DEFINITION                   DEFINITION
!       ECPT      BSC.BEND.TRI.-----TYPE       QUAD.PLT.---------TYPE
!     --------   --------------------------    -------------------------
!     ECPT( 1) = ELEMENT ID         INTEGER ** ELEMENT           INTEGER
!     ECPT( 2) = GRID PT. A         INTEGER ** GRID PT.A         INTEGER
!     ECPT( 3) = GRID PT. B         INTEGER ** GRID PT.B         INTEGER
!     ECPT( 4) = GRID PT. C         INTEGER ** GRID PT.C         INTEGER
!     ECPT( 5) = THETA              REAL    ** GRID PT.D         INTEGER
!     ECPT( 6) = MAT ID 1           INTEGER ** THETA             REAL
!     ECPT( 7) = I  MOM. OF INERT.  REAL    ** MAT ID 1          INTEGER
!     ECPT( 8) = MAT ID 2           INTEGER ** I  MOM. OF INERT. REAL
!     ECPT( 9) = T2                 REAL    ** MAT ID 2          INTEGER
!     ECPT(10) = NON-STRUCT. MASS   REAL    ** T2                REAL
!     ECPT(11) = Z1                 REAL    ** NON-STRUCT. MASS  REAL
!     ECPT(12) = Z2                 REAL    ** Z1                REAL
!     ECPT(13) = COORD. SYS. ID 1   INTEGER ** Z2                REAL
!     ECPT(14) = X1                 REAL    ** COORD. SYS. ID 1  INTEGER
!     ECPT(15) = Y1                 REAL    ** X1                REAL
!     ECPT(16) = Z1                 REAL    ** Y1                REAL
!     ECPT(17) = COORD. SYS. ID 2   INTEGER ** Z1                REAL
!     ECPT(18) = X2                 REAL    ** COORD. SYS. ID 2  INTEGER
!     ECPT(19) = Y2                 REAL    ** X2                REAL
!     ECPT(20) = Z2                 REAL    ** Y2                REAL
!     ECPT(21) = COORD. SYS. ID 3   INTEGER ** Z2                REAL
!     ECPT(22) = X3                 REAL    ** COORD. SYS. ID 3  INTEGER
!     ECPT(23) = Y3                 REAL    ** X3                REAL
!     ECPT(24) = Z3                 REAL    ** Y3                REAL
!     ECPT(25) = ELEMENT TEMP       REAL    ** Z3                REAL
!     ECPT(26) =                            ** COORD. SYS. ID 4  INTEGER
!     ECPT(27) =                            ** X4                REAL
!     ECPT(28) =                            ** Y4                REAL
!     ECPT(29) =                            ** Z4                REAL
!     ECPT(30) =                            ** ELEMENT TEMP      REAL
 
 INTEGER :: subsca,subscb,subscc
 REAL :: ivect,jvect,kvect,d(9)
 DIMENSION       necpt(100),m(12),vq1(3),vq2(3),vq3(3),vq4(3), requiv(10)
 COMMON /condas/ consts(5)
 COMMON /sdr2x5/ ecpt(100),ph1out(128),st(3)
 COMMON /sdr2x6/ a(45),temp15(15),prod15(15),t(9),tite(18),v(25),  &
     d1(3),d2(3),spdum1(18),u1,u2,sinang,cosang,  &
     ssum(60),r(2,5),xsubb,xsubc,ysubc,e(18),temp,  &
     vv1(2),vv2(2),h,a1(3),npoint,spdum2(5),ivect(3),  &
     jvect(3),kvect(3),spdum3(15),theta,nsubc,  &
     spdum4(1),subsca,subscb,subscc,spdum5(2),xc,yc, spdum6(5)
 COMMON /matin / matid,inflag,eltemp,stress,sinth,costh
 COMMON /matout/ g11,g12,g13,g22,g23,g33,rho,alpha(3)
 EQUIVALENCE     (consts(4),degra),(ecpt(1),necpt(1)),  &
     (vq1(1),ecpt(15)),(vq2(1),ecpt(19)), (vq3(1),ecpt(23)),(vq4(1),ecpt(27)),  &
     (requiv(1),r(1,1))
 DATA    m     / 2,4,1,  3,1,2,  4,2,3,  1,3,4 /
 
 idsave = necpt(7)
 eye    = ecpt(8)
 theta  = ecpt(6)*degra
 sinang = SIN(theta)
 cosang = COS(theta)
 
!     FORMATION OF THE R-MATRIX CONTAINING COORDINATES OF THE
!     SUB TRIANGLES. (2X5) FOR QUADRILATERAL PLATE.
!     FORMATION ALSO OF THE I,J, AND K VECTORS USED IN THE E-MATRIX.
 
!     ZERO OUT R-MATRIX
 
 DO  i = 1,10
   requiv(i) = 0.0
 END DO
 
!     SHIFT ECPT UP TO MATCH STRBS1 FOR CERTAIN VARIABLES.
 
 DO  i = 6,12
   ecpt(i) = ecpt(i+1)
 END DO
 
 DO  i = 1,3
   d1(i) = vq3(i) - vq1(i)
   d2(i) = vq4(i) - vq2(i)
   a1(i) = vq2(i) - vq1(i)
 END DO
 
!     NON-NORMALIZED K-VECTOR = D1 CROSS D2
 
 kvect(1) = d1(2)*d2(3) - d2(2)*d1(3)
 kvect(2) = d1(3)*d2(1) - d2(3)*d1(1)
 kvect(3) = d1(1)*d2(2) - d2(1)*d1(2)
 
!     NORMALIZE K-VECTOR
 
 temp = SQRT(kvect(1)**2 + kvect(2)**2 + kvect(3)**2)
 DO  i = 1,3
   kvect(i) = kvect(i)/temp
 END DO
 
!     COMPUTE H = (A1 DOT KVECT)/2
 
 temp = (a1(1)*kvect(1) + a1(2)*kvect(2) + a1(3)*kvect(3))/2.0
 
!     I-VECTOR =(A1) - H*(KVECT)    NON-NORMALIZED
 
 DO  i = 1,3
   ivect(i) = a1(i) - temp*kvect(i)
 END DO
 
!     NORMALIZE I-VECTOR
 
 temp =  SQRT(ivect(1)**2 + ivect(2)**2 + ivect(3)**2)
 DO  i = 1,3
   ivect(i) = ivect(i)/temp
 END DO
 
!     J-VECTOR = K X I  VECTORS
 
 jvect(1) = kvect(2)*ivect(3) - ivect(2)*kvect(3)
 jvect(2) = kvect(3)*ivect(1) - ivect(3)*kvect(1)
 jvect(3) = kvect(1)*ivect(2) - ivect(1)*kvect(2)
 
!     NORMALIZE J VECTOR TO MAKE SURE
 
 temp =  SQRT(jvect(1)**2 + jvect(2)**2 + jvect(3)**2)
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
 
!     STRESS CALCULATION POINT WHICH IS THE DIAGONALS INTERSECTION.
 
 ftemp  = r(1,3)*r(2,4) + r(2,3)*(r(1,2)-r(1,4))
 IF (ftemp == 0.0) CALL mesage (-30,26,ecpt(1))
 r(1,5) = r(1,2)*r(1,3)*r(2,4)/ftemp
 r(2,5) = r(1,2)*r(2,3)*r(2,4)/ftemp
 
!     CHECK OF 4 POINTS FOR ANGLE GREATER THAN OR EQUAL TO 180 DEGREES.
 
 IF (r(2,3) <= 0.0 .OR. r(2,4) <= 0.0) GO TO 90
 temp = r(1,2) - (r(1,2)-r(1,3))*r(2,4)/r(2,3)
 IF (r(1,4) >= temp) GO TO 90
 temp = r(2,3)*r(1,4)/r(2,4)
 IF (r(1,3) > temp) GO TO 100
 90 CALL mesage (-30,35,ecpt(1))
 
!     SET UP THE M-MATRIX FOR MAPPING TRIANGLES, IN DATA STATEMENT
 
!     COMPUTE SUB-TRIANGLE COORDINATES
!     CALL BASIC BENDING ROUTINE FOR ALL SUB-TRIANGLES.
 
 100 eltemp = ecpt(30)
 DO  i = 1,60
   ssum(i) = 0.0
 END DO
 
 DO  j = 1,4
   km     = 3*j - 3
   subsca = m(km+1)
   subscb = m(km+2)
   subscc = m(km+3)
   
   DO  i = 1,2
     vv1(i) = r(i,subscb) - r(i,subsca)
     vv2(i) = r(i,subscc) - r(i,subsca)
   END DO
   xsubb  = SQRT(vv1(1)**2 + vv1(2)**2)
   u1     = vv1(1)/xsubb
   u2     = vv1(2)/xsubb
   xsubc  = u1*vv2(1) + vv2(2)*u2
   ysubc  = u1*vv2(2) - vv2(1)*u2
   
   xc    = SQRT((r(1,subsca)-r(1,5))**2 + (r(2,subsca)-r(2,5))**2)
   yc    = 0.0
   
   sinth = sinang*u1 - cosang*u2
   costh = cosang*u1 + sinang*u2
   IF (ABS(sinth) < 1.0E-06) sinth = 0.0
   
!     AT THIS POINT, XSUBB, XSUBC, YSUBC ARE AT HAND FOR
!     TRIANGLE -J-
   
   CALL strbs1 (1)
   
!     RETURNING FROM STRBS1 THE FOLLOWING QUANTITIES ARE AT HAND.
   
!       S   , S   , S   , EACH 5X3.   45 WORDS STORED IN A( 1)...A(45)
!        A     B     C
   
   
!     SET UP OF T-MATRIX
   
   t(1) = 1.0
   t(2) = 0.0
   t(3) = 0.0
   t(4) = 0.0
   t(5) = u1
   t(6) = u2
   t(7) = 0.0
   t(8) =-u2
   t(9) = u1
   
!     SET UP V-MATRIX PER FMMS 51-A
   
   v( 1) = u1*u1*0.25
   v( 2) = u2*u2*0.25
   v(11) = u1*u2*0.25
   v( 3) =-v(11)*2.0
   v( 4) = 0.0
   v( 5) = 0.0
   v( 6) = v(2)
   v( 7) = v(1)
   v( 8) =-v(3)
   v( 9) = 0.0
   v(10) = 0.0
   v(12) =-v(11)
   v(13) = v(1) - v(2)
   v(14) = 0.0
   v(15) = 0.0
   v(16) = 0.0
   v(17) = 0.0
   v(18) = 0.0
   v(19) = u1*0.25
   v(20) =-u2*0.25
   v(21) = 0.0
   v(22) = 0.0
   v(23) = 0.0
   v(24) =-v(20)
   v(25) = v(19)
   
!     ADD IN S , S , S   TO THE 4 5X3 SSUM MATRICES
!             A   B   C
   
   DO  i = 1,3
     CALL gmmats (v,5,5,0, a(15*i-14),5,3,0, temp15)
     CALL gmmats (temp15,5,3,0, t,3,3,0, prod15)
     
!     POINTER TO SSUM MATRIX
     
     npoint = km + i
     npoint = 15*m(npoint) - 15
     DO  k = 1,15
       nsubc = npoint + k
       ssum(nsubc) = ssum(nsubc) + prod15(k)
     END DO
   END DO
   
 END DO
 
!     FILL E-MATRIX
 
 DO  i = 1,18
   e( i) = 0.0
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
 
 DO  i = 1,4
   
!     DO WE NEED TRANSFORMATION T
!                                I
   nsubc = 4*i + 10
   IF (necpt(nsubc) == 0) GO TO 180
   CALL transs (necpt(nsubc),t)
   CALL gmmats (t,3,3,1, e( 1),3,3,0, tite( 1))
   CALL gmmats (t,3,3,1, e(10),3,3,0, tite(10))
   GO TO 200
   
   180 DO  k = 1,18
     tite(k) = e(k)
   END DO
   
   200 CALL gmmats (ssum(15*i-14),5,3,0, tite,6,3,1, ph1out(30*i-21))
   
 END DO
 
!     I,Z1,Z2,ELEM ID, 4 SILS FOR PHASE 2
 
 ph1out(1) = ecpt( 1)
 ph1out(2) = ecpt( 2)
 ph1out(3) = ecpt( 3)
 ph1out(4) = ecpt( 4)
 ph1out(5) = ecpt( 5)
 ph1out(6) = ecpt( 7)
 ph1out(7) = ecpt(11)
 ph1out(8) = ecpt(12)
 
!     GET S SUB T MATRIX
 
 matid   = idsave
 ecpt(8) = eye
 stress  = 0
 sinth   = sinang
 costh   = cosang
 inflag  = 2
 CALL mat (ecpt(1))
 d(1) = g11*ecpt(8)
 d(2) = g12*ecpt(8)
 d(3) = g13*ecpt(8)
 d(4) = d(2)
 d(5) = g22*ecpt(8)
 d(6) = g23*ecpt(8)
 d(7) = d(3)
 d(8) = d(6)
 d(9) = g33*ecpt(8)
 CALL gmmats (d(1),3,3,0, alpha(1),3,1,0, st(1))
 
!     ALL PHASE ONE COMPLETE
 
 RETURN
END SUBROUTINE sqdpl1
