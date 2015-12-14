SUBROUTINE q4gmgs (mid,factor,g)
!     &    ENTRY Q4GMGD (MID,FACTOD,D)
 
 
!     MATERIAL PROPERTY MATRIX GENERATOR ROUTINE FOR QUAD4 ELEMENT
 
!     THIS ROUTINE BUILDS THE MATERIAL PROPERTY MATRIX, G, USING THE
!     OUTPUT OF SUBROUTINE 'MAT' (/MATOUT/).
 
!     ALL THE MATERIAL OPTIONS, ISOTROPIC, ORTHOTROPIC, AND ANISOTROPIC,
!     ARE AVAILABLE.
 
!     OUTPUT WILL BE G(9) OR D(9) FOR MID1, MID2 AND MID4.
!     FOR MID3,  G(4) OR D(4) IS SENT BACK.
 
 
 INTEGER, INTENT(OUT)                     :: mid
 REAL, INTENT(IN)                         :: factor
 REAL, INTENT(OUT)                        :: g(9)
 REAL :: const,mtype,nu12,nu21
 DOUBLE PRECISION :: factod,d(9),donst
 COMMON  /matout/ rmtout(25)
 EQUIVALENCE      (rmtout(1),e1),(rmtout(2),nu12),(rmtout(3),e2)
 
!     SINGLE PRECISION SECTION -
 
 DO  i=1,9
   g(i) = 0.0
 END DO
 mtype= rmtout(25)
 mtyp = IFIX(mtype+.05) - 2
 IF (mtyp < 0) THEN
   GO TO    20
 ELSE IF (mtyp == 0) THEN
   GO TO    30
 ELSE
   GO TO    80
 END IF
 
!     ISOTROPIC MATERIALS (MAT1)
 
 20 IF (mid /= 3) GO TO 40
 g(1) = rmtout(6)
 g(4) = g(1)
 GO TO 100
 
!     ANISOTROPIC MATERIALS (MAT2)
 
 30 IF (mid == 3) GO TO 60
 40 DO  i=1,3
   g(i) = rmtout(i)
 END DO
 g(4) = g(2)
 g(5) = rmtout(4)
 g(6) = rmtout(5)
 g(7) = g(3)
 g(8) = g(6)
 g(9) = rmtout(6)
 GO TO 100
 
 60 DO  i=1,4
   g(i) = rmtout(i)
 END DO
 g(3) = g(2)
 GO TO 100
 
!     ORTHOTROPIC MATERIALS (MAT8)
 
 80 IF (mid == 3) GO TO 90
 nu21 = nu12 * e2 / e1
 const= 1.0 - (nu21*nu12)
 g(1) = e1 / const
 g(2) = nu12 * e2 / const
 g(4) = g(2)
 g(5) = e2 / const
 g(9) = rmtout(4)
 GO TO 100
 
 90 g(1) = rmtout(6)
 g(4) = rmtout(5)
 IF (g(1) == 0.0 .AND. g(4) == 0.0) GO TO 120
 
!     STANDARD RETURN
 
 100 DO  i=1,9
   g(i) = g(i)*factor
 END DO
 GO TO 310
 
!     FATAL RETURN
 
 120 mid = -mid
 GO TO 310
 
 ENTRY q4gmgd (mid,factod,d)
!     ===========================
 
 DO  i=1,9
   d(i) = 0.0D0
 END DO
 mtype= rmtout(25)
 mtyp = IFIX(mtype+.05) - 2
 IF (mtyp < 0) THEN
   GO TO   210
 ELSE IF (mtyp == 0) THEN
   GO TO   220
 ELSE
   GO TO   270
 END IF
 
!     ISOTROPIC MATERIALS (MAT1)
 
 210 IF (mid /= 3) GO TO 230
 d(1) = rmtout(6)
 d(4) = d(1)
 GO TO 290
 
!     ANISOTROPIC MATERIALS (MAT2)
 
 220 IF (mid == 3) GO TO 250
 230 DO  i=1,3
   d(i) = rmtout(i)
 END DO
 d(4) = d(2)
 d(5) = rmtout(4)
 d(6) = rmtout(5)
 d(7) = d(3)
 d(8) = d(6)
 d(9) = rmtout(6)
 GO TO 290
 
 250 DO  i=1,4
   d(i) = rmtout(i)
 END DO
 d(3) = d(2)
 GO TO 290
 
!     ORTHOTROPIC MATERIALS (MAT8)
 
 270 IF (mid == 3) GO TO 280
 nu21 = nu12 * e2 / e1
 donst= 1.0D0 - DBLE(nu21*nu12)
 d(1) = e1 / donst
 d(2) = nu12 * e2 / donst
 d(4) = d(2)
 d(5) = e2 / donst
 d(9) = rmtout(4)
 GO TO 290
 
 280 d(1) = rmtout(6)
 d(4) = rmtout(5)
 IF (d(1) == 0.0D0 .AND. d(4) == 0.0D0) GO TO 120
 
!     STANDARD RETURN
 
 290 DO  i=1,9
   d(i) = d(i)*factod
 END DO
 
 310 RETURN
END SUBROUTINE q4gmgs
