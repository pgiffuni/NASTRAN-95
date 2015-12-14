SUBROUTINE dbmmov ( index1, index2, no )
     INCLUDE 'ZZZZZZ.COM'
 DO  i = 1, no
   mem( index2+i-1 ) = mem( index1+i-1 )
 END DO
 RETURN
END SUBROUTINE dbmmov
