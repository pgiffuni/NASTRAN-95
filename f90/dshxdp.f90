SUBROUTINE dshxdp ( iarr, LEN )
     
 
 INTEGER, INTENT(IN OUT)                  :: iarr( 10000)
 INTEGER, INTENT(IN)                      :: LEN
 
 WRITE ( 6, 901 ) (iarr(i),i=1,LEN )
 901   FORMAT(200(8(1X,z8),/))
 RETURN
END SUBROUTINE dshxdp
