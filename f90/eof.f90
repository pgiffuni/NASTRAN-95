SUBROUTINE eof ( FILE )
     INCLUDE 'DSIOF.COM'
 
 
 INTEGER, INTENT(IN)                      :: FILE
 
 NAME   = FILE
 iretrn = 0
 CALL dsgefl
 CALL dsefwr
 CALL dssdcb
 RETURN
END SUBROUTINE eof
