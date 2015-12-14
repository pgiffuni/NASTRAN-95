SUBROUTINE smcccs ( ctemp, zil, ilim, zol )
     
 
 COMPLEX, INTENT(OUT)                     :: ctemp( ilim )
 COMPLEX, INTENT(IN)                      :: zil( ilim )
 INTEGER, INTENT(IN)                      :: ilim
 COMPLEX, INTENT(IN)                      :: zol
 
 DO  i = 1, ilim
   ctemp( i ) = ctemp( i ) + zil( i ) * zol
 END DO
 RETURN
END SUBROUTINE smcccs


