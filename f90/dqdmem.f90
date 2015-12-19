SUBROUTINE dqdmem
     
 
!     QUADRILATERAL MEMBRANE ROUTINE FOR DIFFERENTIAL STIFFNESS..
 
 REAL :: ivec,jvec,kvec
 
 DIMENSION          m(12)         ,necpt(5)
 
 COMMON /condas/ consts(5)
 COMMON /ds1aet/    ecpt(100)
 COMMON /matin / matid,inflag,eltemp,stress,sinth,costh
 COMMON /ds1aaa/ npvt,icstm,ncstm ,                  dumcl(32)          ,nogo
 COMMON /ds1adp/    dummy(400)    ,ivec(3) ,ngrid(4)      ,jvec(3)  &
     ,coord(16)     ,kvec(3) ,sdisp(12)     ,pvec(3)  &
     ,vsubk(3)      ,npt1 ,jnot          ,npt2  &
     ,npivot        ,npt3 ,mpoint        ,nsubsc  &
     ,v(3)          ,u1 ,si(3)         ,u2  &
     ,vecl          ,angl ,sinang        ,cosang
 
 EQUIVALENCE ( consts(4) , degra  )
 EQUIVALENCE        (necpt(1),ecpt(1))
 
 DATA  m / 1, 2, 4, 2, 3, 1, 3, 4, 2, 4, 1, 3 /
!     ******************************************************************
!          ECPT                       ECPT
!       RECEIVED BY                REQUIRED BY
!         DQDMEM                     DTRMEM
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
!     ECPT(21) = Z3              ECPT(21) = ELEMENT TEMP.
!     ECPT(22) = COORD. SYS. ID 4ECPT(22) = EL. DEF.
!     ECPT(23) = X4              ECPT(23) = LDTEMP
!     ECPT(24) = Y4              ECPT(24) = XT 1
!     ECPT(25) = Z4              ECPT(25) = YT 1
!     ECPT(26) = ELEMENT TEMP.   ECPT(26) = ZT 1
!     ECPT(27) = EL. DEF.        ECPT(27) = XT 2
!     ECPT(28) = LDTEMP          ECPT(28) = YT 2
!     ECPT(29) = XT 1            ECPT(29) = ZT 2
!     ECPT(30) = YT 1            ECPT(30) = XT 3
!     ECPT(31) = ZT 1            ECPT(31) = YT 3
!     ECPT(32) = XT 2            ECPT(32) = ZT 3
!     ECPT(33) = YT 2
!     ECPT(34) = ZT 2
!     ECPT(35) = XT 3
!     ECPT(36) = YT 3
!     ECPT(37) = ZT 3
!     ECPT(38) = XT 4
!     ECPT(39) = YT 4
!     ECPT(40) = ZT 4
!     ******************************************************************
 
!     THE FOLLOWING COMPUTATION IS PERFORMED FOR USE WITH THE
!     COMPUTATION OF SINTH AND COSTH BELOW (ANISOTROPIC MATERIAL
!     POSSIBILITY)  NOTE  FMMS-46 PAGE -9-
 
 angl = ecpt(6) * degra
 cosang = COS( angl )
 sinang = SIN( angl )
 ivec(1) = ecpt(15) - ecpt(11)
 ivec(2) = ecpt(16) - ecpt(12)
 ivec(3) = ecpt(17) - ecpt(13)
 vecl = SQRT( ivec(1)**2 + ivec(2)**2 + ivec(3)**2 )
 IF (vecl == 0.0E0) GO TO 150
 ivec(1) = ivec(1)/vecl
 ivec(2) = ivec(2)/vecl
 ivec(3) = ivec(3)/vecl
 vsubk(1) =ivec(2) *(ecpt(25)-ecpt(13))-ivec(3)*(ecpt(24)-ecpt(12))
 vsubk(2) =ivec(3) *(ecpt(23)-ecpt(11))-ivec(1)*(ecpt(25)-ecpt(13))
 vsubk(3) =ivec(1) *(ecpt(24)-ecpt(12))-ivec(2)*(ecpt(23)-ecpt(11))
 vecl = SQRT(vsubk(1)**2 + vsubk(2)**2 + vsubk(3)**2 )
 IF (vecl == 0.0E0) GO TO 150
 kvec(1) = vsubk(1)/vecl
 kvec(2) = vsubk(2)/vecl
 kvec(3) = vsubk(3)/vecl
 jvec(1) = kvec(2) * ivec(3) - kvec(3) * ivec(2)
 jvec(2) = kvec(3) * ivec(1) - kvec(1) * ivec(3)
 jvec(3) = kvec(1) * ivec(2) - kvec(2) * ivec(1)
 DO  i=1,3
   pvec(i) = cosang * ivec(i) + sinang * jvec(i)
 END DO
 
 
!     SAVE COORDINATE SYSTEMS, GRID POINT SIL NUMBERS, AND DISP VECTOR.
 
 ngrid(1) = necpt(2)
 ngrid(2) = necpt(3)
 ngrid(3) = necpt(4)
 ngrid(4) = necpt(5)
 
 DO  i=1,16
   coord(i) = ecpt(i + 9)
 END DO
 
 DO  i=1,12
   sdisp(i) = ecpt(i+28)
 END DO
 
!     NOTE. COORD 1, 5, 9, AND 13  ARE INTEGER CSID NUMBERS.
 
!     CORRECT ECPT FOR MEMBRANE USE
 ecpt(5) = ecpt(6)
 ecpt(6) = ecpt(7)
 ecpt(7) = ecpt(8)/2.0E0
 ecpt(8) = ecpt(9)
 ecpt(21) = ecpt(26)
 ecpt(22) = ecpt(27)
 ecpt(23) = ecpt(28)
 
!     FOR EACH TRIANGLE THEN THE THREE GRID POINTS AND COORDINATES
!     ARE INSERTED INTO THE ECPT BEFORE THE CALL TO KTRMEM.
 
!     FILL MAP MATRIX  (PERFORMED IN DATA STATEMENT - DO NOT ALTER)
!              A              B              C
!           M1 = 1         M2 = 2         M3 = 4      (TRIANGLE    I)
 
!           M4 = 2         M5 = 3         M6 = 1      (TRIANGLE   II)
 
!           M7 = 3         M8 = 4         M9 = 2      (TRIANGLE  III)
 
!           M10= 4         M11= 1         M12= 3      (TRIANGLE   IV)
 
!     ******************************************************************
!     FIND WHICH POINT IS THE PIVOT POINT.
 DO  i=1,4
   IF(npvt /= ngrid(i)) CYCLE
   npivot = i
   GO TO 50
 END DO
 
!     FALL THRU ABOVE LOOP IMPLIES AN ERROR CONDITION.
 
 CALL mesage(-30,34,ecpt(1))
 
!     COMPUTE JNOT WHICH EQUALS THE ONE TRIANGLE OF THE FOUR NOT USED
!     AND THUS NOT COMPUTED FOR THE PIVOT POINT IN QUESTION.  (NOTE THE
!     ROWS OF THE MAPPING MATRIX ABOVE AND THE TRIANGLE NUMBERS)
 
 50 IF(npivot - 2 > 0) THEN
   GO TO    70
 END IF
 60 jnot = npivot + 2
 GO TO 80
 70 jnot = npivot - 2
 
 
 80 DO  j=1,4
   IF (j == jnot) GO TO 140
   
!     FILL IN ECPT FOR TRIANGLE J
   mpoint = 3*j - 3
   DO  i=1,3
     npt1 = mpoint + i
     nsubsc = m(npt1)
     necpt(i+1) = ngrid(nsubsc)
     
     npt1 = 3*nsubsc - 3
     npt3 = 3 * i + 20
     DO  k=1,3
       npt2 = npt1 + k
       npt3 = npt3 + 1
       ecpt(npt3) = sdisp(npt2)
     END DO
     
     npt1 = 4*nsubsc - 4
     DO  k=1,4
       npt2 = npt1 + k
       npt3 = 4*i + 4 + k
       ecpt(npt3) = coord(npt2)
     END DO
   END DO
   
!     ECPT IS COMPLETE FOR TRIANGLE J
   
!     SET UP SINTH AND COSTH FOR THIS SUB TRIANGLE
   
   IF( j /= 1 ) GO TO 120
   sinth = sinang
   costh = cosang
   GO TO 130
   
!     NOTE FMMS-46 PAGE-9 FOR FOLLOWING
   
   120 v(1) = ecpt(14) - ecpt(10)
   v(2) = ecpt(15) - ecpt(11)
   v(3) = ecpt(16) - ecpt(12)
   vecl = SQRT( v(1)**2 + v(2)**2 + v(3)**2 )
   IF (vecl == 0.0E0) GO TO 150
   u1 = ( v(1)*pvec(1) + v(2)*pvec(2) + v(3)*pvec(3) )/vecl
   si(1) = v(2) * pvec(3) - v(3) * pvec(2)
   si(2) = v(3) * pvec(1) - v(1) * pvec(3)
   si(3) = v(1) * pvec(2) - v(2) * pvec(1)
   u2 = ( si(1)*kvec(1) + si(2)*kvec(2) + si(3)*kvec(3) )/vecl
   vecl = SQRT( u1**2 + u2**2 )
   u1 = u1 / vecl
   u2 = u2 / vecl
   sinth = sinang * u1 - cosang * u2
   costh = cosang * u1 + sinang * u2
   130 IF( ABS(sinth) < 1.0E-06 ) sinth = 0.0E0
   
   CALL dtrmem(   1   )
   
!     INSERTIONS ARE PERFORMED BY DTRMEM
   
   140 CONTINUE
 END DO
 
!     ******************************************************************
 
 RETURN
 150 CALL mesage(30,26,ecpt(1))
 
!  SET FLAG FOR FATAL ERROR WHILE ALLOWING ERROR MESSAGES TO ACCUMULATE
 
 nogo=1
 RETURN
END SUBROUTINE dqdmem
