SUBROUTINE curvps (sigs, prin)
!*****
!  COMPUTES PRINCIPAL STRESSES OR STRAINS AND ANGLE OF MAXIMUM.
!*****
 
 
 REAL, INTENT(IN)                         :: sigs(3)
 REAL, INTENT(OUT)                        :: prin(4)
 
 
 temp = sigs(1) - sigs(2)
 prin(4) = SQRT( (temp/2.0)**2 + sigs(3)**2 )
 delta = ( sigs(1) + sigs(2) ) / 2.0
 prin(2) = delta + prin(4)
 prin(3) = delta - prin(4)
 delta = 2.0 * sigs(3)
 IF( ABS(delta) < 1.0E-15 .AND. ABS(temp) < 1.0E-15 ) GO TO 100
 prin(1) = ATAN2( delta, temp )*28.6478898E0
 RETURN
 
 100 prin(1) = 0.0
 RETURN
END SUBROUTINE curvps
