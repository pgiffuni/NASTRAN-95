SUBROUTINE qdmem(t,core)
     
 REAL, INTENT(IN)                         :: t(1)
 REAL, INTENT(IN OUT)                     :: core(1)
 REAL :: ivec,jvec,kvec,ngrid
 DIMENSION m(12),r(6),ngrid(4),coord(16)
 DIMENSION d1(3),d2(3),a1(3),a2(3),a3(3),a4(3),ivec(3),jvec(3),  &
     kvec(3),v(8),ecptsa(36)
 COMMON /condas/    pi       ,twopi    ,radeg    ,degra    , s4pisq
 COMMON /trimex/ ecpt(26)
 COMMON /matin / matid,inflag,eltemp,stress,sinth,costh
 EQUIVALENCE (r(1),ivec(1)),(ngrid(1),ecptsa(2)),  &
     (coord(1),ecptsa(10)),(r(4),jvec(1))
 
 DATA  m / 1, 2, 4, 2, 3, 1, 3, 4, 2, 4, 1, 3 /
!     ******************************************************************
!          ECPT                       ECPT
!       RECEIVED BY                REQUIRED BY
!         SQDME1                     STRME1
!     ******************************************************************
!     ECPT( 1) = EL. ID          ECPT( 1) = EL. ID
!     ECPT( 2) = GRD. PT. A      ECPT( 2) = GRD. PT. A
!     ECPT( 3) = GRD. PT. B      ECPT( 3) = GRD. PT. B
!     ECPT( 4) = GRD. PT. C      ECPT( 4) = GRD. PT. C
!     ECPT( 5) = GRD. PT. D      ECPT( 5) = THETA
!     ECPT( 6) = THETA           ECPT( 6) = MATERIAL ID
!     ECPT( 7) = MATERIAL ID     ECPT( 7) = T
!     ECPT( 8) = T               ECPT( 8) = NON-STRUCT. MASS
!     ECPT( 9) = NON-STRUCT. MASSECPT( 9) = COORD. SYS. ID 1
!     ECPT(10) = COORD. SYS. ID 1ECPT(10) = X1
!     ECPT(11) = X1              ECPT(11) = Y1
!     ECPT(12) = Y1              ECPT(12) = Z1
!     ECPT(13) = Z1              ECPT(13) = COORD. SYS. ID 2
!     ECPT(14) = COORD. SYS. ID 2ECPT(14) = X2
!     ECPT(15) = X2              ECPT(15) = Y2
!     ECPT(16) = Y2              ECPT(16) = Z2
!     ECPT(17) = Z2              ECPT(17) = COORD. SYS. ID 3
!     ECPT(18) = COORD. SYS. ID 3ECPT(18) = X3
!     ECPT(19) = X3              ECPT(19) = Y3
!     ECPT(20) = Y3              ECPT(20) = Z3
!     ECPT(21) = Z3              ECPT(21) = ELEMENT TEMPERATURE
!     ECPT(22) = COORD. SYS. ID 4    NOTE. THE FOLLOWING ARE INTEGERS...
!     ECPT(23) = X4                  GRID POINTS, MAT ID, EL.ID,
!     ECPT(24) = Y4                  COORD. SYS. IDS.
!     ECPT(25) = Z4                  ALL OTHERS ARE REAL IN THE ECPT.
!     ECPT(26) = ELEMENT TEMPERATURE
!     ******************************************************************
 
 
!     VECTORS D1 AND D2  FMMS-46 PAGE 6
!     A1 A2 A3 A4
 
 DO  i=1,3
   d1(i) = ecpt(i + 18) - ecpt(i + 10)
   d2(i) = ecpt(i + 22) - ecpt(i + 14)
   a1(i) = ecpt(i + 14) - ecpt(i + 10)
   a2(i) = ecpt(i + 18) - ecpt(i + 14)
   a3(i) = ecpt(i + 22) - ecpt(i + 18)
   a4(i) = ecpt(i + 10) - ecpt(i + 22)
 END DO
 
!     K-VECTOR = NORMALIZED D1 CROSS D2
 
 kvec(1) = d1(2) * d2(3) - d1(3) * d2(2)
 kvec(2) = d1(3) * d2(1) - d1(1) * d2(3)
 kvec(3) = d1(1) * d2(2) - d1(2) * d2(1)
 vecl = SQRT ( kvec(1)**2 + kvec(2)**2 + kvec(3)**2 )
 IF(vecl == 0.0) CALL mesage(-30,26,ecpt(1))
 kvec(1) = kvec(1)/vecl
 kvec(2) = kvec(2)/vecl
 kvec(3) = kvec(3)/vecl
 
!     I-VECTOR = NORMALIZED A SUB 12 - H * KVECTOR
!     GET H FIRST = ( A SUB 12 DOT KVECTOR)/2
 
 h = (a1(1)*kvec(1) + a1(2)*kvec(2) + a1(3)*kvec(3))/2.0E0
 
 ivec(1) = a1(1) - h * kvec(1)
 ivec(2) = a1(2) - h * kvec(2)
 ivec(3) = a1(3) - h * kvec(3)
 vecl = SQRT ( ivec(1)**2 + ivec(2)**2 + ivec(3)**2 )
 IF(vecl == 0.0) CALL mesage(-30,26,ecpt(1))
 ivec(1) = ivec(1)/vecl
 ivec(2) = ivec(2)/vecl
 ivec(3) = ivec(3)/vecl
 
!     J-VECTOR = K CROSS I
 
 jvec(1) = kvec(2) * ivec(3) - kvec(3) * ivec(2)
 jvec(2) = kvec(3) * ivec(1) - kvec(1) * ivec(3)
 jvec(3) = kvec(1) * ivec(2) - kvec(2) * ivec(1)
 
 vecl = SQRT(jvec(1)**2 + jvec(2)**2 + jvec(3)**2)
 jvec(1) = jvec(1)/vecl
 jvec(2) = jvec(2)/vecl
 jvec(3) = jvec(3)/vecl
 
 theta=ecpt(6)*degra
 sinang=SIN(theta)
 cosang=COS(theta)
 
 v(1) = 1.0E0
 v(2) = 0.0E0
 
!     R ARRAY IS EQUIVALENCED TO IVECTOR AND JVECTOR
 
 CALL gmmats(r,2,3,0,  a2,3,1,0,  v(3))
 CALL gmmats(r,2,3,0,  a3,3,1,0,  v(5))
 CALL gmmats(r,2,3,0,  a4,3,1,0,  v(7))
 
!     NORMALIZE THE 4 2X1 V ARRAYS
 
 DO  i=1,4
   vecl = SQRT ( v(2*i-1)**2 + v(2*i)**2 )
   IF(vecl == 0.0) CALL mesage(-30,26,ecpt(1))
   v(2*i-1) = v(2*i-1)/vecl
   v(2*i  ) = v(2*i  )/vecl
 END DO
 
!     MAPPING MATRIX M IS IN DATA STATEMENT.
 
!     NOW MAKE 4 CALLS TO STRME1 WHICH WILL RETURN
 
!     SAVE GRID SILS AND COORDINATE SYSTEMS.
 
 
!     REDUCE THICKNESS BY 0.5
 
 ecpt(8) = ecpt(8)/2.0
 DO  i=1,36
   ecptsa(i) = ecpt(i)
 END DO
 
 ecpt(6) = ecpt(7)
 ecpt(7) = ecpt(8)
 ecpt(8) = ecpt(9)
 
 ecpt(21) = ecpt(26)
 
 DO  i=1,4
   
!     POINTER TO THE SILS IN THE MAPPING MATRIX
   ncoord = 8
   npoint = 3*i-3
   tbar = t(1)
   DO  j=2,4
     npoint = npoint + 1
     nsub1 = m(npoint)
     DO  k=1,4
       nsub3 = 4*nsub1 - 4 + k
       ncoord = ncoord + 1
       ecpt(ncoord) = coord(nsub3)
     END DO
     ecpt(j) = ngrid(nsub1)
   END DO
   
!     SET UP T MATRIX FOR THIS TRIANGLE.  T IS 3X3
   
   u1 = v(2*i-1)
   u2 = v(2*i  )
   
   
!     COMPUTE NET SINTH AND COSTH FOR ANISOTROPIC POSSIBILITY
   
   sinth = sinang * u1 - cosang * u2
   costh = cosang * u1 + sinang * u2
   
   CALL trimem(1,tbar,core)
 END DO
 RETURN
END SUBROUTINE qdmem
