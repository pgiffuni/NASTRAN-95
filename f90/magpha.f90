SUBROUTINE magpha( a, b )
!*****
! THIS SUBROUTINE FORMS THE MAGNITUDE OF (A,B) AND STORES IT IN A...
! THE PHASE OF (X=A, Y=B) IS THEN FORMED AND THE RESULT STORED IN B...
!*****
 
 REAL, INTENT(IN OUT)                     :: a
 REAL, INTENT(IN OUT)                     :: b
 COMMON /condas/ consts(5)
 
 EQUIVALENCE ( consts(3) , radeg  )
 
 value = SQRT( a**2 + b**2 )
 IF( value  == 0.0) THEN
   GO TO    20
 END IF
 10 phase = ATAN2( b, a ) * radeg
 GO TO 30
 20 phase = 0.0E0
 GO TO 40
 30 IF( phase < (-0.00005E0) ) phase = phase + 360.0E0
 40 a = value
 b = phase
 RETURN
END SUBROUTINE magpha
