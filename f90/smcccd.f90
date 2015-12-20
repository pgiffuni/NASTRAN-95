SUBROUTINE smcccd ( dtemp, zil, ilim, zol )

 INTEGER, INTENT(IN)                      :: ilim
 DOUBLE COMPLEX, INTENT(IN)               :: zil( ilim ),zol
 DOUBLE COMPLEX   ::  dtemp( ilim )

 DO  i = 1, ilim
   dtemp( i ) = dtemp( i ) + zil( i ) * zol
 END DO

 RETURN
END SUBROUTINE smcccd


