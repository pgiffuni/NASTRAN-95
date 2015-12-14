SUBROUTINE pkqdpl
     
!     THIS ROUTINE CALCULATES THE STIFFNESS MATRIX FOR QUAD-PLATES IN
!     PLA4
 
!     THIS ROUTINE GENERATES THE FOLLOWING
 
!     FOUR 6X6 STIFFNESS MATRICES WITH RESPECT TO ONE PIVOT POINT OF A
!     QUADRILATERAL PLATE ELEMENT.
 
!     REF.  FMMS-44   JULY  18, 1967   TRI.BENDING ELEMENT STIFF.
!           FMMS-48   AUGUST 1, 1967   QUAD. BENDING ELEMENT STIFF.
 
!     CALLS FROM THIS ROUTINE ARE MADE TO
!           PKTRBS - BASIC BENDING TRIANGLE
!           TRANSD - SUPPLIES 3X3 TRANSFORMATIONS
!           PLA4B  - INSERTION ROUTINE
!           GMMATD - GENERAL MATRIX MULITPLY AND TRANSPOSE ROUTINE
!           MESAGE - ERROR MESSAGE WRITER
 
 
!     ECPT LISTS AS OF AUGUST 4, 1967
 
!                 DEFINITION                   DEFINITION
!       ECPT      BSC.BEND.TRI.-----TYPE       QUAD.PLT.---------TYPE
!     ==================================================================
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
 DOUBLE PRECISION :: kout,tite,tjte,dpdum1,dpdum2,d1,d2,ivect,jvect,  &
                     kvect,a1,ksum,t,v,vv,xsubb,xsubc,ysubc,prod9,  &
                     temp,temp9,h,u1,u2,e,a,temp18,requiv,r
 DIMENSION        m(12),necpt(100),requiv(8),vq1(3),vq2(3),vq3(3), vq4(3),a(1)
 COMMON /condas/  consts(5)
 COMMON /matin /  matid,inflag,eltemp,stress,sinth,costh
 COMMON /matout/  g11,g12,g13,g22,g23,g33,rho,alpha1,alpha2,alp12,  &
                  t_sub_0,g_sub_e,sigten,sigcom,sigshe, g2x211,g2x212,g2x222
 COMMON /pla42c/  npvt,dum1(148),nogo
 COMMON /pla4es/  ecpt(100)
 COMMON /pla42d/  kout(36),tite(18),tjte(18),temp18(18),dpdum1(54),  &
                  ivect(3),jvect(3),kvect(3),d1(3),d2(3),a1(3),  &
                  t(9),v(2),vv(2),h,u1,u2,r(2,4),ksum(36),  &
                  dpdum2(3),prod9(9),temp9(9),xsubb,xsubc,ysubc,  &
                  e(18),temp,sp1(28),sp2(2),km,nbegin,jnot,npivot,  &
                  theta,nsubc,ising,subsca,subscb,subscc,sinang, cosang,npoint
 EQUIVALENCE      (consts(4),degra),(necpt(1),ecpt(1)),  &
                  (r(1,1),requiv(1)),(vq1(1),ecpt(15)), (vq2(1),ecpt(19)),(vq3(1),ecpt(23)),  &
                  (vq4(1),ecpt(27)),(a(1),kout(1))
 DATA    m     /  2,4,1, 3,1,2, 4,2,3, 1,3,4 /
 
!     DETERMINE PIVOT POINT NUMBER
 
 DO  i = 1,4
   IF (npvt /= necpt(i+1)) CYCLE
   npivot = i
   GO TO 20
 END DO
 
!     FALL THRU ABOVE LOOP IMPLIES ERROR CONDITION
 
 CALL mesage (-30,34,ecpt(1))
 
 20 theta  = ecpt(6)*degra
 sinang = SIN(theta)
 cosang = COS(theta)
 
 IF (npivot-2 > 0) THEN
   GO TO    40
 END IF
 30 jnot = npivot + 2
 GO TO 50
 40 jnot = npivot - 2
 
!     FORMATION OF THE R-MATRIX CONTAINING COORDINATES OF THE
!     SUB TRIANGLES.  (2X4) FOR QUADRILATERAL PLATE...
!     FORMATION ALSO OF THE I,J, AND K VECTORS USED IN THE E-MATRIX.
 
!     ZERO OUT R-MATRIX
 
 50 DO  i = 1,8
   requiv(i) = 0.0D0
 END DO
 
!     SHIFT ECPT UP TO MATCH PKTRBS FOR CERTAIN VARIABLES.
 
 DO  i = 6,12
   ecpt(i) = ecpt(i+1)
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
 IF (temp == 0.0D0) GO TO 330
 DO  i = 1,3
   kvect(i) = kvect(i)/temp
 END DO
 
!     COMPUTE H = (A1 DOT KVECT)/2
 
 temp = (a1(1)*kvect(1) + a1(2)*kvect(2) + a1(3)*kvect(3))/2.0D0
 
!     I-VECTOR =(A1) - H*(KVECT)    NON-NORMALIZED
 
 DO  i = 1,3
   ivect(i) = a1(i) - temp*kvect(i)
 END DO
 
!     NORMALIZE I-VECTOR
 
 temp = DSQRT(ivect(1)**2 + ivect(2)**2 + ivect(3)**2)
 IF (temp == 0.0D0) GO TO 330
 DO  i = 1,3
   ivect(i) = ivect(i)/temp
 END DO
 
!     J-VECTOR = K CROSS I, AND X3 CALCULATION
 
 jvect(1) = kvect(2)*ivect(3) - ivect(2)*kvect(3)
 jvect(2) = kvect(3)*ivect(1) - ivect(3)*kvect(1)
 jvect(3) = kvect(1)*ivect(2) - ivect(1)*kvect(2)
 
!     NORMALIZE J VECTOR TO MAKE SURE
 
 temp =  DSQRT(jvect(1)**2 + jvect(2)**2 + jvect(3)**2)
 IF (temp == 0.0D0) GO TO 330
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
 
!     CHECK OF 4 POINTS FOR ANGLE GREATER THAN OR EQUAL TO 180 DEGREES.
 
 IF (r(2,3) <= 0.0D0 .OR. r(2,4) <= 0.0D0) GO TO 140
 temp = r(1,2) - (r(1,2)-r(1,3))*r(2,4)/r(2,3)
 IF (r(1,4) >= temp) GO TO 140
 temp = r(2,3)*r(1,4)/ r(2,4)
 IF (r(1,3) > temp) GO TO 150
 140 CALL mesage (30,35,ecpt(1))
 
!     SET FLAG FOR FATAL ERROR WHILE ALLOWING ERROR MESSAGES TO
!     ACCUMULATE
 
 nogo = 1
 RETURN
 
! 140 AT 140 THE COORDINATES OF THE PLATE IN THE ELEMENT
!     SYSTEM ARE STORED IN THE R-MATRIX WHERE THE COLUMN DENOTES THE
!     POINT AND THE ROW DENOTES THE X OR Y COORDINATE FOR ROW 1 OR
!     ROW 2 RESPECTIVELY.
 
 
!     SET UP THE M-MATRIX FOR MAPPING TRIANGLES, IN DATA STATEMENT.
 
 
!     COMPUTE SUB-TRIANGLE COORDINATES
 
!     ZERO OUT KSUM MATRICES
 
 150 DO  i = 1,36
   ksum(i) = 0.0D0
 END DO
 
 DO  j = 1,4
   IF (j == jnot) CYCLE
   km = 3*j - 3
   subsca = m(km+1)
   subscb = m(km+2)
   subscc = m(km+3)
   
   DO  i = 1,2
     v(i) = r(i,subscb) - r(i,subsca)
     vv(i)= r(i,subscc) - r(i,subsca)
   END DO
   xsubb = DSQRT(v(1)**2 + v(2)**2)
   u1 = v(1)/xsubb
   u2 = v(2)/xsubb
   xsubc = u1*vv(1) + u2*vv(2)
   ysubc = u1*vv(2) - u2*vv(1)
   
   sinth = sinang*u1  -  cosang*u2
   costh = cosang*u1  +  sinang*u2
   IF (ABS(sinth) < 1.0E-06) sinth = 0.0E0
   
!     AT THIS POINT, XSUBB, XSUBC, YSUBC ARE AT HAND FOR
!     TRIANGLE -J-
   
   CALL pktrbs (1)
!                         U
!     NOW HAVE AT HAND  K    I,J, =1,2,3.   9-3X3 MATRICES STORED AT
!                        IJ                 A(1) THROUGH A(81).
   
!     MAP THE 3 3X3-S FOR THE PIVOT ROW INTO THE SUMMATION ARRAYS...
   
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
   
!     FIND WHICH POINT OF THE SUBTRIANGLE IS ALSO THE PIVOT OF THE
!     QUADRILATERAL
   
   DO  i = 1,3
     npoint = km + i
     IF (m(npoint) /= npivot) CYCLE
     nbegin = 27*i - 27
     GO TO 190
   END DO
   
   190 DO  i = 1,3
     npoint = nbegin + 9*i - 8
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
 
 IF (necpt(4*npivot+10) == 0) GO TO 240
 CALL transd (necpt(4*npivot+10),t)
 CALL gmmatd (t,3,3,1, e( 1),3,3,0, tite( 1))
 CALL gmmatd (t,3,3,1, e(10),3,3,0, tite(10))
 GO TO 260
 
 240 DO  k = 1,18
   tite(k) = e(k)
 END DO
 
 260 DO  j = 1,4
   
!     TRANSFORMATIONS AND INSERTION
   
   IF (necpt(4*j+10) == 0) GO TO 270
   CALL transd (necpt(4*j+10),t)
   CALL gmmatd (t,3,3,1,  e(1),3,3,0, tjte(1 ))
   CALL gmmatd (t,3,3,1, e(10),3,3,0, tjte(10))
   GO TO 290
   270 DO  k = 1,18
     tjte(k) = e(k)
   END DO
   290 CALL gmmatd (ksum(9*j-8),3,3,0, tjte,6,3,1, temp18(1))
   CALL gmmatd (tite(1),6,3,0, temp18(1),3,6,0, kout(1))
   CALL pla4b (kout(1),necpt(j+1))
   
 END DO
 RETURN
 
 330 CALL mesage (30,26,ecpt(1))
 
!     SET FLAG FOR FATAL ERROR WHILE ALLOWING ERROR MESSAGES TO
!     ACCUMULATE
 
 nogo = 1
 
 RETURN
END SUBROUTINE pkqdpl
