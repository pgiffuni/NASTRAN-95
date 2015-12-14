SUBROUTINE q2trmd(ra,rb,rc,alpha,isinth,icosth,gsube,it,  &
        ierror,iopt,kmat,pmat,smat,zmat)
!*****
!  SUB-TRIANGLE COMPUTATION ROUTINE FOR THE QDMEM2 ELEMENT
 
!  ON INPUT
!  ========
!            RA,RB,RC = 3 (3X1) COORDINATE VECTORS FOR TRIANGLE
!            IOPT     = 1  CALL FROM STIFFNESS GENERATION MODULE
!                     = 2  CALL FROM STATIC LOAD MODULE
!                     = 3  CALL FROM STRESS RECOVERY MODULE
!            ALPHA    = 3X1 VECTOR APPROPRIATE FOR CALL
!            ISINTH   = SIN OF MATERIAL ANGLE(WHOLE - ELEMENT)
!            ICOSTH   = COS OF MATERIAL ANGLE(WHOLE - ELEMENT)
!            GSUBE    = MATERIAL MATRIX (3X3)
!            IT       = THICKNESS OF ELEMENT
 
!  ON OUTPUT
!  =========
!            IERROR = 0  IF NO ERROR
!                   = 1  IF BAD ELEMENT GEOMETRY
 
!            KMAT,PMAT,SMAT,ZMAT = FOLLOWING PER IOPT VALUE SENT
 
 
!            IOPT=1
!            ------
!            KMAT = 7 (3X3)-S = KCA,KCB,KCC,KAA,KAB,KBA,KBB
!            PMAT = UNCHANGED
!            SMAT = UNCHANGED
!            ZMAT = UNCHANGED
 
!            IOPT=2
!            ------
!            KMAT = 3 (3X3)-S = KCA,KCB,KCC
!            PMAT = 3 (3X1)-S = PA,PB,PC
!            SMAT = UNCHANGED
!            ZMAT = UNCHANGED
 
!            IOPT=3
!            ------
!            KMAT = 7 (3X3)-S = KCA,KCB,KCC,KAA,KAB,KBA,KBB
!            PMAT = 3 (3X1)-S = PTA,PTB,PTC
!            SMAT = 3 (3X3)-S = SA,SB,SC
!            ZMAT = 3 (3X1)-S = ZA,ZB,ZC
 
!*****
 
 DOUBLE PRECISION, INTENT(IN)             :: ra(3)
 DOUBLE PRECISION, INTENT(IN)             :: rb(3)
 DOUBLE PRECISION, INTENT(IN)             :: rc(3)
 DOUBLE PRECISION, INTENT(IN)             :: alpha(3)
 DOUBLE PRECISION, INTENT(IN)             :: isinth
 DOUBLE PRECISION, INTENT(IN)             :: icosth
 DOUBLE PRECISION, INTENT(IN OUT)         :: gsube(9)
 DOUBLE PRECISION, INTENT(IN)             :: it
 INTEGER, INTENT(OUT)                     :: ierror
 INTEGER, INTENT(IN)                      :: iopt
 DOUBLE PRECISION, INTENT(IN OUT)         :: kmat(1)
 DOUBLE PRECISION, INTENT(OUT)            :: pmat(1)
 DOUBLE PRECISION, INTENT(IN OUT)         :: smat(1)
 DOUBLE PRECISION, INTENT(IN OUT)         :: zmat(1)
 DOUBLE PRECISION :: e(9)     ,iareat   ,i33
 DOUBLE PRECISION :: mag      ,ivec(3)
 DOUBLE PRECISION :: v12(3)   ,iarea    ,jvec(3)
 DOUBLE PRECISION :: v13(3)   ,ixsubb   ,kvec(3)
 DOUBLE PRECISION :: ic       , ca(6)    ,ixsubc   ,tm(9)
 DOUBLE PRECISION :: is       , cb(6)    ,iysubc   ,tm3(3)
 DOUBLE PRECISION :: cc(6)    ,c(3,6)   ,alp(3)
 DOUBLE PRECISION :: temp9(9) ,hi(27)   ,hitge(9)
 
 DOUBLE PRECISION :: dadotb
 
 INTEGER :: ipart(3,3)
 
 EQUIVALENCE (c(1,1),ca(1)), (c(1,3),cb(1)), (c(1,5),cc(1))
 EQUIVALENCE (e(1),ivec(1)), (e(4),jvec(1)), (e(7),kvec(1))
 
 DATA ipart/ 28,46,1,   37,55,10,   0,0,19 /
 
!     V    = R   - R   ,      V    = R   - R
!      12     B     A          13     C     B
 
 DO  i = 1,3
   v12(i) = rb(i) - ra(i)
   v13(i) = rc(i) - ra(i)
 END DO
 
!     KVEC(UN-NORMALIZED)  =  V     X  V
!                              12       13
 
 CALL daxb( v12, v13, kvec )
 mag = DSQRT( dadotb(kvec,kvec) )
 IF( mag  > 0) THEN
   GO TO    20
 ELSE
   GO TO   190
 END IF
 
!     NORMALIZE  K-VECTOR, AND AREA
 
 20 kvec(1) = kvec(1) / mag
 kvec(2) = kvec(2) / mag
 kvec(3) = kvec(3) / mag
 iarea = 0.50D0 * mag
 
!     I-VECTOR = V   (NORMALIZED) THUS
!                 12
 
 mag = DSQRT( dadotb( v12, v12 ) )
 IF( mag  > 0) THEN
   GO TO    30
 ELSE
   GO TO   190
 END IF
 30 ivec(1) = v12(1) / mag
 ivec(2) = v12(2) / mag
 ivec(3) = v12(3) / mag
 ixsubb = mag
 
!     J-VECTOR = K-VECTOR CROSS I-VECTOR THUS
 
 CALL daxb( kvec, ivec, jvec )
 
!     MATERIAL COEFFICIENTS C AND S    U,V,W = I-VECTOR
 
 mag = DSQRT( ivec(1)**2 + ivec(2)**2 )
 IF( mag <= 0.d0 ) GO TO 190
 ic =(ivec(1)*icosth + ivec(2)*isinth)/mag
 is =(ivec(1)*isinth - ivec(2)*icosth)/mag
 
!     X = MAGNITUDE OF V  , X = I-VEC DOT V  , Y = J-VEC DOT V
!      B                12   C             13   C             13
 
 ixsubc = dadotb( ivec, v13 )
 iysubc = dadotb( jvec, v13 )
 IF( ixsubb  == 0) THEN
   GO TO   190
 END IF
 40 IF( iysubc  == 0) THEN
   GO TO   190
 END IF
 
 50 ca(1) = -1.0D0 / ixsubb
 ca(2) = 0.0D0
 i33 = 1.0D0 / iysubc
 ca(3) = i33 * (ixsubc/ixsubb - 1.0D0)
 ca(4) = 0.0D0
 ca(5) = ca(3)
 ca(6) = ca(1)
 
 cb(1) = -ca(1)
 cb(2) = 0.0D0
 cb(3) = - i33   * (ixsubc / ixsubb )
 cb(4) = 0.0D0
 cb(5) = cb(3)
 cb(6) = cb(1)
 
 cc(1) = 0.0D0
 cc(2) = 0.0D0
 cc(3) = i33
 cc(4) = 0.0D0
 cc(5) = i33
 cc(6) = 0.0D0
 
!     FORM MATERIAL-ORIENTATION-TRANSFORMATION-MATRIX  (BY-ROWS)
 
 tm(1) = ic * ic
 tm(2) = is * is
 tm(3) = ic * is
 tm(4) = tm(2)
 tm(5) = tm(1)
 tm(6) = -tm(3)
 tm(7) = 2.0D0 * tm(6)
 tm(8) = -tm(7)
 tm(9) = tm(1) - tm(2)
 iareat= iarea * it
 
!     IF SSG CALL MULTIPLY ALPHA(T-TO) VECTOR BY IAREAT
 
 IF( iopt /= 2 ) GO TO 60
 alp(1) = alpha(1) * iareat
 alp(2) = alpha(2) * iareat
 alp(3) = alpha(3) * iareat
 
!     IF SDR CALL COMPUTE AREA   = X  * T
!                                   B
 60 IF( iopt /= 3 ) GO TO 70
 tm3(1) = tm(3) * it
 tm3(2) = tm(6) * it
 tm3(3) = tm(9) * it
 
!     SET FIRST PARTITION ROW TO COMPUTE FOR STIFFNESS MATRICES.
 
 70 irow1 = 1
 IF( iopt == 2 ) irow1 = 3
!*****
!           M
!     H  = T  C  E
!      I       I
 
!*****
 DO  i = 1,3
   CALL gmmatd( tm,3,3,0, c(1,2*i-1),2,3,1, temp9 )
   CALL gmmatd( temp9,3,2,0, e,2,3,0, hi(9*i-8) )
 END DO
!*****
!     FORM OUTPUTS FOR POINTS I = A,B,C
!*****
 DO  i = 1,3
   
!              T
!     HITGE= H  G
!             I  E
   
   CALL gmmatd( hi(9*i-8),3,3,1, gsube,3,3,0, hitge )
   
!     STIFFNESS MATRIX CALCULATIONS
   
!     ONLY KAA,KAB     ARE FORMED.  OUTPUT ORDER WITH EACH 3X3 STORED
!          KBA,KBB                  BY ROWS =
!          KCA,KCB,KCC              KCA,KCB,KCC,KAA,KAB,KBA,KBB
   
   IF( i < irow1 ) GO TO 150
   kk = 0
   DO  j = 1,3
     ipartn = ipart(i,j)
     IF( ipartn  > 0) THEN
       GO TO    90
     ELSE
       GO TO   140
     END IF
     90 DO  k = 1,9
       kk = kk + 1
       temp9(k) = hi(kk)*iareat
     END DO
     CALL gmmatd( hitge,3,3,0, temp9,3,3,0, kmat(ipartn) )
   END DO
   150 SELECT CASE ( iopt )
     CASE (    1)
       GO TO 180
     CASE (    2)
       GO TO 160
     CASE (    3)
       GO TO 170
   END SELECT
!****
!  SSG LOAD GENERATION CALL ADDITIONAL DATA TO OUTPUT.
   
!  ONLY PA,PB,PC ARE FORMED.
!*****
   160 CALL gmmatd( hitge,3,3,0,   alp,3,1,0, pmat(3*i-2) )
   CYCLE
!*****
!  SDR ADDITIONAL PHASE-1 STRESS OUTPUTS
!*****
   170 jpartn = 9*i-8
   CALL gmmatd( gsube,3,3,0, hi(jpartn),3,3,0, smat(jpartn) )
   ipartn = 3*i - 2
   CALL gmmatd( hitge,3,3,0, alpha,3,1,0, pmat(ipartn)  )
   CALL gmmatd( tm3,3,1,1, smat(jpartn),3,3,0, zmat(ipartn) )
   DO  j=1,3
     k = ipartn + j - 1
     pmat(k) = pmat(k)*iareat
   END DO
   
 END DO
 ierror = 0
 RETURN
!*****
!  ERROR CONDITION, BAD GEOMETRY.
!*****
 190 ierror = 1
 RETURN
END SUBROUTINE q2trmd
