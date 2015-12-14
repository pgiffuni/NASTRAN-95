SUBROUTINE logfil (line)
     
 
 INTEGER, INTENT(IN OUT)                  :: line(18)
 
 
 COMMON /logout/ lout
 
 WRITE (lout, 2000) line
 RETURN
 
 2000 FORMAT (1X, 18A4)
 
END SUBROUTINE logfil
