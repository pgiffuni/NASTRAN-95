SUBROUTINE mcone
     
!     MASS MATRIX GENERATION FOR AXIS-SYMETRIC CONICAL SHELL ELEMENT
 
!     ECPT( 1) = ELEMENT ID             INTEGER        ECT
!     ECPT( 2) = SIL PT A               INTEGER        ECT
!     ECPT( 3) = SIL PT B B             INTEGER        ECT
!     ECPT( 4) = MATID 1                INTEGER        EPT
!     ECPT( 5) = T   (MEMBRANE THICK)   REAL           EPT
!     ECPT( 6) = MATID 2                INTEGER        EPT
!     ECPT( 7) = I   (MOM.OF INERTIA)   REAL           EPT
!     ECPT( 8) = MATID 3                INTEGER        EPT
!     ECPT( 9) = TS  (SHEAR THICKNESS)  REAL           EPT
!     ECPT(10) = NON-STRUCTURAL-MASS    REAL           EPT
!     ECPT(11) = Z1                     REAL           EPT
!     ECPT(12) = Z2                     REAL           EPT
!     ECPT(13) = PHI  1                 REAL           EPT
!     ECPT(14) = PHI  2                 REAL           EPT
!     ECPT(15) = PHI  3                 REAL           EPT
!     ECPT(16) = PHI  4                 REAL           EPT
!     ECPT(17) = PHI  5                 REAL           EPT
!     ECPT(18) = PHI  6                 REAL           EPT
!     ECPT(19) = PHI  7                 REAL           EPT
!     ECPT(20) = PHI  8                 REAL           EPT
!     ECPT(21) = PHI  9                 REAL           EPT
!     ECPT(22) = PHI 10                 REAL           EPT
!     ECPT(23) = PHI 11                 REAL           EPT
!     ECPT(24) = PHI 12                 REAL           EPT
!     ECPT(25) = PHI 13                 REAL           EPT
!     ECPT(26) = PHI 14                 REAL           EPT
!     ECPT(27) = COORD. SYS. ID PT.1    INTEGER        BGPDT
!     ECPT(28) = RADIUS PT. 1           REAL           BGPDT
!     ECPT(29) = DISTANCE TO PT.1       REAL           BGPDT
!     ECPT(30) = NULL                   REAL           BGPDT
!     ECPT(31) = COORD. SYS. ID PT.2    INTEGER        BGPDT
!     ECPT(32) = RADIUS PT 2            REAL           BGPDT
!     ECPT(33) = DISTANCE TO PT. 2      REAL           BGPDT
!     ECPT(34) = NULL                   REAL           BGPDT
!     ECPT(35) = ELEMENT TEMPERATURE    REAL           GEOM3
 
 INTEGER :: necpt(100)
 REAL :: l,mu
 DOUBLE PRECISION :: mass
 COMMON /condas/  pi, twopi, radeg, degra, s4pisq
 COMMON /matin /  matid, inflag, eltemp
 COMMON /matout/  rho
 COMMON /sma2et/  ecpt(100)
 COMMON /sma2io/  dum4(10), ifmgg, dum5(25)
 COMMON /sma2cl/  dum3(2),npvt
 COMMON /sma2dp/  mass(36), temp, l, term, m1
 EQUIVALENCE      (ra,ecpt(28)), (rb,ecpt(32)), (za,ecpt(29)),  &
     (zb,ecpt(33)), (t,ecpt(5))  , (mu,ecpt(10)), (necpt(1),ecpt(1))
 
 l = SQRT((rb-ra)**2 + (zb-za)**2)
 
!     NEXT LINE WAS REMOVED BY M.H./NAVY. ERROR FOR CONICAL SHELL MASS
 
 
 temp = rb/6.0 + ra/3.0
 
 IF (t == 0.0) THEN
   GO TO    40
 END IF
 30 inflag = 4
 matid  = necpt(4)
 eltemp = ecpt(35)
 CALL mat (necpt(1))
 40 DO  i = 1,36
   mass(i) = 0.0D0
 END DO
 term = pi*l*temp*(rho*t + mu)
 IF (necpt(1)-(necpt(1)/1000)*1000-1 == 0) term = term*2.0
 mass( 1) = term
 mass( 8) = term
 mass(15) = term
 m1 = -1
 CALL sma2b (mass(1),npvt,m1,ifmgg,0.0D0)
 RETURN
END SUBROUTINE mcone
