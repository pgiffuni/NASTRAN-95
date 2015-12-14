SUBROUTINE smcccd ( dtemp, zil, ilim, zol )
     
 
 REAL, INTENT(OUT)                        :: dtemp
 DOUBLE COMPLEX      dtem, INTENT(IN)     :: zil( ilim )
 INTEGER, INTENT(IN)                      :: ilim
 DOUBLE COMPLEX      dtem, INTENT(IN)     :: zol
 DOUBLE COMPLEX      dtemp( ilim )
 DO  i = 1, ilim
   dtemp( i ) = dtemp( i ) + zil( i ) * zol
 END DO
 RETURN
END SUBROUTINE smcccd


