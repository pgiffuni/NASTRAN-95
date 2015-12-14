SUBROUTINE lprops (g)
!     &    ENTRY LPROPD (D)
 
!     THIS ROUTINE RETURNS INTRINSIC LAYER PROPERTIES FOR
!     ALL LAYERS REFERENCING MAT1, MAT2 OR MAT8 PROPERTY
!     ENTRIES IN A STANDARD FORMAT AS REQUIRED FOR FILE PCOMPS
 
 
 REAL, INTENT(OUT)                        :: g(25)
 REAL :: mtype,nu12,nu21
 DOUBLE PRECISION :: d(25),donst
 COMMON  /matout/ rmtout(25)
 EQUIVALENCE      (rmtout(1),e1),(rmtout(2),nu12),(rmtout(3),e2)
 
!     SINGLE PRECISION -
 
 DO  i=1,25
   g(i) = 0.0
 END DO
 mtype = rmtout(25)
 mtyp = IFIX(mtype+.05) - 2
 IF (mtyp < 0) THEN
   GO TO    20
 ELSE IF (mtyp == 0) THEN
   GO TO    30
 ELSE
   GO TO    60
 END IF
 
!****
!     ISOTROPIC MATERIALS, MAT1 IN MAT2 FORMAT
 
!****  LAYER PROPERTY MATRIX
 
 20 CONTINUE
 
!****
!     ANISOTROPIC MATERIALS, MAT2
 
!****  LAYER PROPERTY MATRIX
 
 30 DO  i=1,3
   g(i) = rmtout(i)
 END DO
 g(5) = rmtout(4)
 g(6) = rmtout(5)
 g(9) = rmtout(6)
 g(4) = g(2)
 g(7) = g(3)
 g(8) = g(6)
 
!**** TRANSVERSE SHEAR FLEXIBILITY MATRIX
 DO  i=10,13
   ii = i - 9
   g(i) = rmtout(ii)
 END DO
 g(12) = g(11)
 
!**** THERMAL COEFFICIENTS OF EXPANSION
 g(14) = rmtout( 8)
 g(15) = rmtout( 9)
 g(16) = rmtout(10)
 
!**** ULTIMATE STRENGTH VALUES
 g(17) = rmtout(13)
 g(18) = rmtout(13)
 g(19) = rmtout(14)
 g(20) = rmtout(14)
 g(21) = rmtout(15)
 g(22) = 0.0
 
!*** RHO, TREF, GE
 g(23) = rmtout( 7)
 g(24) = rmtout(11)
 g(25) = rmtout(12)
 GO TO 160
 
!****
!     ORTHOTROPIC MATERIALS, MAT8
 
!****  LAYER PROPERTY MATRIX
 
 60 nu21 = nu12 * e2 / e1
 const= 1.0 - (nu21*nu12)
 g(1) = e1 / const
 g(2) = nu12 * e2 / const
 g(5) = e2 / const
 g(4) = g(2)
 g(9) = rmtout(4)
 
!**** TRANSVERSE SHEAR FLEXIBILITY MATRIX
 g(10) = rmtout(6)
 g(13) = rmtout(5)
 
!**** THERMAL COEFFICIENTS OF EXPANSION
 g(14) = rmtout(8)
 g(15) = rmtout(9)
 
!**** ULTIMATE STRENGTH VALUES
 g(17) = rmtout(11)
 g(18) = rmtout(12)
 g(19) = rmtout(13)
 g(20) = rmtout(14)
 g(21) = rmtout(15)
 g(22) = rmtout(17)
 
!*** RHO, TREF, GE
 g(23) = rmtout( 7)
 g(24) = rmtout(10)
 g(25) = rmtout(16)
 GO TO 160
 
 ENTRY lpropd (d)
!     ================
 
!     DOUBLE PRECISION -
 
 DO  i=1,25
   d(i) = 0.0D0
 END DO
 mtype = rmtout(25)
 mtyp  = IFIX(mtype+.05) - 2
 IF (mtyp < 0) THEN
   GO TO   110
 ELSE IF (mtyp == 0) THEN
   GO TO   120
 ELSE
   GO TO   150
 END IF
 
!****
!     ISOTROPIC MATERIALS, MAT1 IN MAT2 FORMAT
 
!****  LAYER PROPERTY MATRIX
 
 110 CONTINUE
 
!****
!     ANISOTROPIC MATERIALS, MAT2
 
!****  LAYER PROPERTY MATRIX
 
 120 DO  i=1,3
   d(i) = rmtout(i)
 END DO
 d(5) = rmtout(4)
 d(6) = rmtout(5)
 d(9) = rmtout(6)
 d(4) = d(2)
 d(7) = d(3)
 d(8) = d(6)
 
!**** TRANSVERSE SHEAR FLEXIBILITY MATRIX
 DO  i=10,13
   ii = i - 9
   d(i) = rmtout(ii)
 END DO
 d(12) = d(11)
 
!**** THERMAL COEFFICIENTS OF EXPANSION
 d(14) = rmtout( 8)
 d(15) = rmtout( 9)
 d(16) = rmtout(10)
 
!**** ULTIMATE STRENGTH VALUES
 d(17) = rmtout(13)
 d(18) = rmtout(13)
 d(19) = rmtout(14)
 d(20) = rmtout(14)
 d(21) = rmtout(15)
 d(22) = 0.0D0
 
!*** RHO, TREF, GE
 d(23) = rmtout( 7)
 d(24) = rmtout(11)
 d(25) = rmtout(12)
 GO TO 160
 
!****
!     ORTHOTROPIC MATERIALS, MAT8
 
!****  LAYER PROPERTY MATRIX
 
 150 nu21 = nu12 * e2 / e1
 donst= 1.0D0 - DBLE(nu21*nu12)
 d(1) = e1 / donst
 d(2) = nu12 * e2 / donst
 d(5) = e2 / donst
 d(4) = d(2)
 d(9) = rmtout(4)
 
!**** TRANSVERSE SHEAR FLEXIBILITY MATRIX
 d(10) = rmtout(6)
 d(13) = rmtout(5)
 
!**** THERMAL COEFFICIENTS OF EXPANSION
 d(14) = rmtout(8)
 d(15) = rmtout(9)
 
!**** ULTIMATE STRENGTH VALUES
 d(17) = rmtout(11)
 d(18) = rmtout(12)
 d(19) = rmtout(13)
 d(20) = rmtout(14)
 d(21) = rmtout(15)
 d(22) = rmtout(17)
 
!*** RHO, TREF, GE
 d(23) = rmtout( 7)
 d(24) = rmtout(10)
 d(25) = rmtout(16)
 
 160 RETURN
END SUBROUTINE lprops
