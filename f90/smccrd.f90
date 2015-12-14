SUBROUTINE smccrd ( dtemp, zil, ilim, zol )
     
 
 DOUBLE PRECISION, INTENT(OUT)            :: dtemp( ilim )
 DOUBLE PRECISION, INTENT(IN)             :: zil( ilim )
 INTEGER, INTENT(IN)                      :: ilim
 DOUBLE PRECISION, INTENT(IN)             :: zol
 
 DO  i = 1, ilim
   dtemp( i ) = dtemp( i ) + zil( i ) * zol
 END DO
 RETURN
END SUBROUTINE smccrd


