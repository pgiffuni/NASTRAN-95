SUBROUTINE q2bcd (est,planar,rmat,et,ierror)
     
!     BASIC CALCULATIONS ARE PERFORMED FOR THE QDMEM2 ELEMENT IN THIS
!     ROUTINE (DOUBLE-PRECISION VERSION)
 
 
 REAL, INTENT(IN)                         :: est(1)
 LOGICAL, INTENT(OUT)                     :: planar
 DOUBLE PRECISION, INTENT(OUT)            :: rmat(3,5)
 DOUBLE PRECISION, INTENT(OUT)            :: et(3,3)
 INTEGER, INTENT(OUT)                     :: ierror
 
 
 DOUBLE PRECISION :: mag    ,d12(3) ,g1(3) ,iarea  ,d13(3),grid(3,5),  &
     g2(3)  ,itwoh  ,d24(3),vec(3) ,g3(3) , g5(3)  ,g4(3)  ,dadotb
 EQUIVALENCE      (grid(1,1),g1(1)),(grid(1,2),g2(1)),  &
     (grid(1,3),g3(1)),(grid(1,4),g4(1)), (grid(1,5),g5(1))
 
!     MOVE GRID COORDINATES AND MAKE DOUBLE-PRECISION IF THIS IS THE
!     DOUBLE-PRECISION VERSION.
 
 DO  i = 1,3
   g1(i) = est(i+10)
   g2(i) = est(i+14)
   g3(i) = est(i+18)
   g4(i) = est(i+22)
 END DO
 
!     FORM  D   , D   AND  D   VECTORS
!            13    24       12
 
 DO  i = 1,3
   d12(i) = g2(i) - g1(i)
   d13(i) = g3(i) - g1(i)
   d24(i) = g4(i) - g2(i)
 END DO
 
!     NVEC = D13 CROSS D24 = K-VECTOR (UN-NORMALIZED)
 
 CALL daxb (d13,d24,vec)
 mag   = DSQRT(dadotb(vec,vec))
 iarea = 0.5D0*mag
 
!     NORMALIZE K-VECTOR
 
 IF (mag > 0) THEN
   GO TO    30
 ELSE
   GO TO   100
 END IF
 30 et(1,3) = vec(1)/mag
 et(2,3) = vec(2)/mag
 et(3,3) = vec(3)/mag
 
!     H = .5 * (D   DOT K-VEC)
!                12
 
 itwoh = dadotb(d12,et(1,3))
 
!     I-VECTOR (UN-NORMALIZED) = (D  ) - 2 H (K-VECTOR)
!                                  12
 
 DO  i = 1,3
   vec(i) = d12(i) - itwoh*et(i,3)
 END DO
 mag = DSQRT(dadotb(vec,vec))
 
!     NORMALIZE I-VECTOR
 
 IF (mag > 0) THEN
   GO TO    50
 ELSE
   GO TO   100
 END IF
 50 et(1,1) = vec(1)/mag
 et(2,1) = vec(2)/mag
 et(3,1) = vec(3)/mag
 
!     JVEC = KVEC CROSS IVEC
 
 CALL daxb (et(1,3),et(1,1),et(1,2))
 
!     FILL THE SUB-TRIANGLE ELEMENT COORDINATE MATRIX
 
 DO  i = 1,3
   g5(i) = 0.25D0*(g1(i) + g2(i) + g3(i) + g4(i))
 END DO
 rmat(1,1) = 0.0D0
 rmat(2,1) = 0.0D0
 rmat(3,1) =-itwoh/2.0D0
 DO  i = 2,5
   vec(1) = grid(1,i) - g1(1)
   vec(2) = grid(2,i) - g1(2)
   vec(3) = grid(3,i) - g1(3)
   CALL gmmatd (et,3,3,0, vec,3,1,0, rmat(1,i))
   rmat(1,i) = rmat(1,i) + rmat(1,1)
   rmat(2,i) = rmat(2,i) + rmat(2,1)
   rmat(3,i) = rmat(3,i) + rmat(3,1)
 END DO
 
!     SET PLANAR FLAG .TRUE. OR .FALSE.
 
 IF ((itwoh/2.0D0)**2/iarea <= 0.01D0) GO TO 80
 planar = .false.
 GO TO 90
 80 planar = .true.
 
!     ALL BASIC CALCULATIONS NOW COMPLETE
 
 90 ierror = 0
 RETURN
 
!     ERROR CONDITION, BAD ELEMENT GEOMETRY.
 
 100 ierror = 1
 RETURN
END SUBROUTINE q2bcd
