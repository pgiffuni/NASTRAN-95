SUBROUTINE smccrs ( temp, zil, ilim, zol )
     
 
 REAL, INTENT(OUT)                        :: temp( ilim )
 REAL, INTENT(IN)                         :: zil( ilim )
 INTEGER, INTENT(IN)                      :: ilim
 REAL, INTENT(IN)                         :: zol
 
 DO  i = 1, ilim
   temp( i ) = temp( i ) + zil( i ) * zol
 END DO
 RETURN
END SUBROUTINE smccrs


