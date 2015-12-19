SUBROUTINE sqdme1
     
!          ECPT                        ECPT
!       RECEIVED BY                 REQUIRED BY
!         SQDME1                      STRME1
!     -----------------------     --------------------------
!     ECPT( 1) = EL. ID           ECPT( 1) = EL. ID
!     ECPT( 2) = GRD. PT. A       ECPT( 2) = GRD. PT. A
!     ECPT( 3) = GRD. PT. B       ECPT( 3) = GRD. PT. B
!     ECPT( 4) = GRD. PT. C       ECPT( 4) = GRD. PT. C
!     ECPT( 5) = GRD. PT. D       ECPT( 5) = THETA
!     ECPT( 6) = THETA            ECPT( 6) = MATERIAL ID
!     ECPT( 7) = MATERIAL ID      ECPT( 7) = T
!     ECPT( 8) = T                ECPT( 8) = NON-STRUCT. MASS
!     ECPT( 9) = NON-STRUCT. MASS ECPT( 9) = COORD. SYS. ID 1
!     ECPT(10) = COORD. SYS. ID 1 ECPT(10) = X1
!     ECPT(11) = X1               ECPT(11) = Y1
!     ECPT(12) = Y1               ECPT(12) = Z1
!     ECPT(13) = Z1               ECPT(13) = COORD. SYS. ID 2
!     ECPT(14) = COORD. SYS. ID 2 ECPT(14) = X2
!     ECPT(15) = X2               ECPT(15) = Y2
!     ECPT(16) = Y2               ECPT(16) = Z2
!     ECPT(17) = Z2               ECPT(17) = COORD. SYS. ID 3
!     ECPT(18) = COORD. SYS. ID 3 ECPT(18) = X3
!     ECPT(19) = X3               ECPT(19) = Y3
!     ECPT(20) = Y3               ECPT(20) = Z3
!     ECPT(21) = Z3               ECPT(21) = ELEMENT TEMPERATURE
!     ECPT(22) = COORD. SYS. ID 4
!     ECPT(23) = X4
!     ECPT(24) = Y4
!     ECPT(25) = Z4
!     ECPT(26) = ELEMENT TEMPERATURE
 
!     NOTE. THE FOLLOWING ARE INTEGERS - GRID POINTS, MAT ID, EL.ID,
!                                        COORD. SYS. IDS.
!           ALL OTHERS ARE REAL IN THE ECPT.
 
 INTEGER :: necpt(100)
 REAL :: ivec,jvec,kvec
 DIMENSION       m(12),r(6),ngrid(4),coord(16),ssubt(3),s(27)
 COMMON /condas/ consts(5)
 COMMON /sdr2x6/ dummy(100),sum(36),stemp(9),d1(3),d2(3),a1(3),  &
     a2(3),a3(3),a4(3),ivec(3),jvec(3),kvec(3),vecl,h,  &
     v(8),ecptsa(36),st(3),ncoord,npoint,nsub1,nsub2,  &
     nsub3,t(9),cosang,sinang,u1,u2,dumy(61)
 COMMON /sdr2x5/ ecpt(100),ph1out(100),forvec(25)
 COMMON /matin / matid,inflag,eltemp,stress,sinth,costh
 EQUIVALENCE     (consts(4),degra),(necpt(1),ecpt(1)),  &
     (r(1),ivec(1)),(ngrid(1),ecptsa(2)),  &
     (coord(1),ecptsa(10)),(s(1),ph1out(10)), (ssubt(1),ph1out(7))
 DATA    m     / 1, 2, 4, 2, 3, 1, 3, 4, 2, 4, 1, 3 /
 
 
 angl   = ecpt(6)*degra
 cosang = COS(angl)
 sinang = SIN(angl)
 
!     VECTORS D1 AND D2  FMMS-46 PAGE 6
!     A1 A2 A3 A4
 
 DO  i = 1,3
   d1(i) = ecpt(i+18) - ecpt(i+10)
   d2(i) = ecpt(i+22) - ecpt(i+14)
   a1(i) = ecpt(i+14) - ecpt(i+10)
   a2(i) = ecpt(i+18) - ecpt(i+14)
   a3(i) = ecpt(i+22) - ecpt(i+18)
   a4(i) = ecpt(i+10) - ecpt(i+22)
 END DO
 
!     K-VECTOR = NORMALIZED D1 CROSS D2
 
 kvec(1) = d1(2)*d2(3) - d1(3)*d2(2)
 kvec(2) = d1(3)*d2(1) - d1(1)*d2(3)
 kvec(3) = d1(1)*d2(2) - d1(2)*d2(1)
 vecl    = SQRT(kvec(1)**2 + kvec(2)**2 + kvec(3)**2)
 IF (vecl < 1.0E-06) CALL mesage (-30,26,ecpt(1))
 kvec(1) = kvec(1)/vecl
 kvec(2) = kvec(2)/vecl
 kvec(3) = kvec(3)/vecl
 
!     I-VECTOR = NORMALIZED A SUB 12 - H * KVECTOR
!     GET H FIRST = (A SUB 12 DOT KVECTOR)/2
 
 h = (a1(1)*kvec(1) + a1(2)*kvec(2) + a1(3)*kvec(3))/2.0
 
 ivec(1) = a1(1) - h*kvec(1)
 ivec(2) = a1(2) - h*kvec(2)
 ivec(3) = a1(3) - h*kvec(3)
 vecl    = SQRT(ivec(1)**2 + ivec(2)**2 + ivec(3)**2)
 IF (vecl < 1.0E-06) CALL mesage (-30,26,ecpt(1))
 ivec(1) = ivec(1)/vecl
 ivec(2) = ivec(2)/vecl
 ivec(3) = ivec(3)/vecl
 
!     J-VECTOR = K CROSS I
 
 jvec(1) = kvec(2)*ivec(3) - kvec(3)*ivec(2)
 jvec(2) = kvec(3)*ivec(1) - kvec(1)*ivec(3)
 jvec(3) = kvec(1)*ivec(2) - kvec(2)*ivec(1)
 
 vecl    = SQRT(jvec(1)**2 + jvec(2)**2 + jvec(3)**2)
 jvec(1) = jvec(1)/vecl
 jvec(2) = jvec(2)/vecl
 jvec(3) = jvec(3)/vecl
 
 
 v(1) = 1.0
 v(2) = 0.0
 
!     R ARRAY IS EQUIVALENCED TO IVECTOR AND JVECTOR
 
 CALL gmmats (r,2,3,0, a2,3,1,0, v(3))
 CALL gmmats (r,2,3,0, a3,3,1,0, v(5))
 CALL gmmats (r,2,3,0, a4,3,1,0, v(7))
 
!     NORMALIZE THE 4 2X1 V ARRAYS
 
 DO  i = 1,4
   vecl = SQRT(v(2*i-1)**2 + v(2*i)**2)
   IF (vecl < 1.0E-10) CALL mesage (-30,26,ecpt(1))
   v(2*i-1) = v(2*i-1)/vecl
   v(2*i  ) = v(2*i  )/vecl
 END DO
 
!     MAPPING MATRIX M IS IN DATA STATEMENT.
 
!     NOW MAKE 4 CALLS TO STRME1 WHICH WILL RETURN
!     S , S , S , S , T SUB 0
!      A   B   C   T
 
!     SAVE GRID SILS AND COORDINATE SYSTEMS.
 
 DO  i = 1,36
   ecptsa(i) = ecpt(i)
 END DO
 
 ecpt(6) = ecpt(7)
 ecpt(7) = ecpt(8)
 ecpt(8) = ecpt(9)
 
!     ZERO OUT SUM MATRICES
 
 DO  i = 1,36
   sum(i) = 0.0
 END DO
 st(1) = 0.0
 st(2) = 0.0
 st(3) = 0.0
 
 ecpt(21) = ecpt(26)
 
 DO  i = 1,4
   
!     POINTER TO THE SILS IN THE MAPPING MATRIX
   
   ncoord = 8
   npoint = 3*i - 3
   DO  j = 2,4
     npoint = npoint + 1
     nsub1  = m(npoint)
     DO  k = 1,4
       nsub3  = 4*nsub1 - 4 + k
       ncoord = ncoord + 1
       ecpt(ncoord) = coord(nsub3)
     END DO
     necpt(j) = ngrid(nsub1)
   END DO
   
!     SET UP T MATRIX FOR THIS TRIANGLE.  T IS 3X3
   
   u1   = v(2*i-1)
   u2   = v(2*i  )
   
   t(1) = u1**2
   t(2) = u2**2
   t(7) = u1*u2
   t(3) =-2.0*t(7)
   t(4) = t(2)
   t(5) = t(1)
   t(6) =-t(3)
   t(8) =-t(7)
   t(9) = t(1) - t(2)
   
!     COMPUTE NET SINTH AND COSTH FOR ANISOTROPIC POSSIBILITY
   
   sinth = sinang*u1 - cosang*u2
   costh = cosang*u1 + sinang*u2
   
   CALL strme1 (1)
   
   
!     NOW TRANSFORM AND ADD THE S MATRICES INTO THE RESPECTIVE SUM
!     MATRICES.
   
   DO  j = 1,3
     
!     POINTER TO TRIANGLE I ROW IN THE MAPPING MATRIX
     
     npoint = 3*i - 3
     
!     TRANSFORM S
     
     CALL gmmats (t,3,3,0, s(9*j-8),3,3,0, stemp)
     
!     ADD STEMP INTO RESPECTIVE KSUM POSITIONS
     
!     ZERO POINTER INTO KSUM MATRICES
     
     nsub1 = npoint + j
     nsub1 = m(nsub1)*9 - 9
     DO  k = 1,9
       nsub1 = nsub1 + 1
       sum(nsub1) = sum(nsub1) + stemp(k)
     END DO
   END DO
   
!     TRANSFORM AND ADD IN S SUB T
   
   CALL gmmats (t,3,3,0, ssubt, 3,1,0, stemp)
   st(1) = st(1) + stemp(1)
   st(2) = st(2) + stemp(2)
   st(3) = st(3) + stemp(3)
 END DO
 
!     ALL MATRICES COMPLETE
 
!     FILL OUTPUT BLOCK
 
 DO  i = 1,5
   ph1out(i) = ecptsa(i)
 END DO
 ph1out(7) = st(1)*0.25
 ph1out(8) = st(2)*0.25
 ph1out(9) = st(3)*0.25
 DO  i = 1,36
   ph1out(i+9) = 0.25*sum(i)
 END DO
 
!     PHASE 1 COMPLETE OUTPUT BLOCK CONTAINS 45 WORDS
 
 RETURN
END SUBROUTINE sqdme1
