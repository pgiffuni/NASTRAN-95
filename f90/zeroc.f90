SUBROUTINE zeroc(iz,n)
     
!     SET AND ARRAY TO ZERO
 
 
 
 INTEGER, INTENT(OUT)                     :: iz(n)
 INTEGER, INTENT(IN)                      :: n
 
 
 DO  i=1,n
   iz(i) = 0
 END DO
 RETURN
END SUBROUTINE zeroc
