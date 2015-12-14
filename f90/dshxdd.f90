SUBROUTINE dshxdd ( ii,iarr, LEN)
     
 INTEGER, INTENT(OUT)                     :: ii
 INTEGER, INTENT(IN)                      :: iarr( 10000)
 INTEGER, INTENT(IN)                      :: LEN
 
 COMMON / dsbuff / iibuff(2048)
 
 DO  k = 1, LEN
   iibuff(k) = iarr(k)
 END DO
 WRITE ( 6, 901 ) ii,(iibuff(i),i=1,LEN )
 901   FORMAT(i5,200(8(1X,z8),/))
 RETURN
END SUBROUTINE dshxdd
